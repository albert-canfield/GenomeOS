# SPDX-License-Identifier: AGPL-3.0-or-later
"""Run scripts/enhancer_targets_all.py over every chromosome, one after another.

Albert's instruction of 2026-09-13: let the running chr21 job finish, then continue chromosome by
chromosome until every one is scored. This driver changes nothing in the scoring script (owned by the
area I decoding lane); it only sequences it. The order is smallest first, as that lane asked, so the
shared quota turns into finished chromosomes sooner.

    uv run python scripts/enhancer_targets_all_chain.py            # wait for chr21, then chr1..chrY
    uv run python scripts/enhancer_targets_all_chain.py --dry-run  # say what it would do

Rules it keeps:
- one chromosome at a time: every run spends the same AlphaGenome quota, so two in parallel only race;
- a chromosome whose result says `"complete": true` is skipped, so the chain resumes after any restart
  (the scoring script itself resumes element by element from its cache);
- a chromosome another process is already scoring is waited on, never started twice;
- each chromosome is started through the job registry (`jobs.start`), so the Progress tab shows it and
  the registry refuses to start a twin of a live run;
- a run that ends without completing is retried after a pause;
- a completed chromosome's summary result is committed alone through a private index and pushed; a
  result still carrying the per-element table (the old shape, megabytes) is refused and logged, never
  committed. The tables stay local under data/knowledge/alphagenome/all_elements.
Progress goes to data/jobs/enhancer_targets_all_chain.log and each chromosome's own log beside it.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

from genomeos import jobs
from genomeos.jobs import heartbeat

JOB = "enhancer_targets_all_chain"
RESULTS = Path("data/results")
JOBS = Path("data/jobs")
FIRST = "chr21"  # already running when the chain was asked for
# smallest first (the decoding lane's order of 2026-09-13); chr21 is the run already in progress
ORDER = [
    "chr22", "chrY", "chr19", "chr20", "chr18", "chr17", "chr16", "chr15", "chr14", "chr13", "chr12",
    "chr11", "chr10", "chr9", "chr8", "chrX", "chr7", "chr6", "chr5", "chr4", "chr3", "chr2", "chr1",
]  # fmt: skip
POLL = 60
RETRY_PAUSE = 600
MAX_RETRIES = 20
PIDFILE = JOBS / f"{JOB}.pid"


def log(msg: str) -> None:
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    print(line, flush=True)
    JOBS.mkdir(parents=True, exist_ok=True)
    with open(JOBS / f"{JOB}.log", "a") as fh:
        fh.write(line + "\n")


def complete(chrom: str) -> bool:
    p = RESULTS / f"enhancer_targets_all_{chrom}.json"
    if not p.exists():
        return False
    try:
        return bool(json.loads(p.read_text()).get("complete"))
    except (OSError, ValueError):
        return False


def scored(chrom: str) -> tuple[int, int] | None:
    p = RESULTS / f"enhancer_targets_all_{chrom}.json"
    try:
        d = json.loads(p.read_text())
        return int(d.get("scored", 0)), int(d.get("elements_total", 0))
    except (OSError, ValueError):
        return None


def running(chrom: str) -> bool:
    """Is some process already scoring this chromosome (started by the registry, a person or this chain)?"""
    # match the interpreter running the script, not any shell whose command line mentions it (a peer's
    # wait loop can contain the same words)
    pattern = f"[Pp]ython[^ ]* scripts/enhancer_targets_all.py --chrom {chrom}( |$)"
    out = subprocess.run(["pgrep", "-f", pattern], capture_output=True, text=True)
    return out.returncode == 0 and bool(out.stdout.strip())


def wait_for(chrom: str, why: str) -> None:
    last = 0.0
    while running(chrom) or (chrom == FIRST and not complete(chrom)):
        heartbeat(JOB)
        if time.time() - last > 1800:
            last = time.time()
            sc = scored(chrom)
            state = f"{sc[0]:,}/{sc[1]:,} scored" if sc else "no result yet"
            log(f"waiting on {chrom} ({why}): {state}")
        time.sleep(POLL)


def committable(chrom: str) -> bool:
    """A complete summary in the committed shape: `complete` true and no per-element table inside."""
    p = RESULTS / f"enhancer_targets_all_{chrom}.json"
    try:
        d = json.loads(p.read_text())
    except (OSError, ValueError):
        return False
    return bool(d.get("complete")) and "elements" not in d


def git(*args: str, env: dict | None = None, timeout: int = 900) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *args], capture_output=True, text=True, env=env, timeout=timeout)


def commit_result(chrom: str) -> bool:
    """Commit one chromosome's summary alone, through a private index, and push. False on any refusal."""
    path = f"data/results/enhancer_targets_all_{chrom}.json"
    if not committable(chrom):
        log(f"{chrom}: result not in the committed shape (per-element table inside?); not committed")
        return False
    d = json.loads((RESULTS / f"enhancer_targets_all_{chrom}.json").read_text())
    s_ = d.get("summary") or {}
    named = s_.get("fraction_with_target")
    msg = (
        f"All-enhancer deletion scoring, {chrom} complete: {d.get('scored', 0):,} elements inside a node "
        f"deleted in AlphaGenome with the effect per cell line"
        + (f", {named:.1%} name a gene" if isinstance(named, float) else "")
        + "; the summary is committed and the per-element table stays local "
        "(data/knowledge/alphagenome/all_elements)"
    )
    index = JOBS / f"{JOB}.index"
    for attempt in range(3):
        env = {**os.environ, "GIT_INDEX_FILE": str(index.resolve())}
        old = git("rev-parse", "refs/heads/dev").stdout.strip()
        steps = [
            git("read-tree", old, env=env),
            git("add", path, env=env),
        ]
        if any(r.returncode for r in steps):
            log(f"{chrom}: staging failed: {' '.join(r.stderr.strip() for r in steps)[:200]}")
            return False
        if not git("diff", "--cached", "--quiet", old, env=env).returncode:
            log(f"{chrom}: summary already committed")
            return True
        tree = git("write-tree", env=env).stdout.strip()
        new = git("commit-tree", tree, "-p", old, "-m", msg).stdout.strip()
        if not new:
            log(f"{chrom}: commit-tree failed")
            return False
        moved = git("update-ref", "refs/heads/dev", new, old)
        index.unlink(missing_ok=True)
        if moved.returncode == 0:
            git("read-tree", "HEAD")  # keep the shared index in step with the branch (LESSONS)
            log(f"{chrom}: committed {new[:7]}")
            push = git("push", "origin", "dev")
            if push.returncode == 0:
                log(f"{chrom}: pushed")
            else:
                log(f"{chrom}: push refused, left for the next push: {push.stderr.strip()[-200:]}")
            return True
        log(f"{chrom}: dev moved while committing (attempt {attempt + 1}); rebuilding on the new tip")
    return False


def run_one(chrom: str) -> bool:
    """Score one chromosome to completion through the registry, retrying runs that end early."""
    name = f"enhancer_targets_all_{chrom}"
    for attempt in range(1, MAX_RETRIES + 1):
        if complete(chrom):
            break
        if running(chrom):
            wait_for(chrom, "already running elsewhere")
            continue
        log(f"{chrom}: starting through the registry (attempt {attempt})")
        jobs.start(name)
        time.sleep(POLL)
        while running(chrom):
            heartbeat(JOB)
            time.sleep(POLL)
        log_tail = (JOBS / f"{name}.log").read_text(errors="replace").splitlines()[-3:]
        if any("AlphaGenome is disabled" in line for line in log_tail):
            log(f"{chrom}: AlphaGenome is disabled; the chain stops here")
            return False
        if complete(chrom):
            break
        pause = RETRY_PAUSE // 60
        log(f"{chrom}: run ended before completing ({' | '.join(log_tail)[-160:]}); retry in {pause} min")
        time.sleep(RETRY_PAUSE)
    if not complete(chrom):
        log(f"{chrom}: not complete after {MAX_RETRIES} attempts; the chain stops here")
        return False
    sc = scored(chrom)
    log(f"{chrom}: complete, {sc[0]:,} elements" if sc else f"{chrom}: complete")
    commit_result(chrom)
    return True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    todo = [c for c in ORDER if c != FIRST and not complete(c)]
    if args.dry_run:
        print(f"first: {FIRST} complete={complete(FIRST)} running={running(FIRST)} scored={scored(FIRST)}")
        print(f"then, in order: {' '.join(todo)}")
        print(f"already complete: {' '.join(c for c in ORDER if c != FIRST and complete(c)) or 'none'}")
        return 0
    if PIDFILE.exists():
        try:
            other = int(PIDFILE.read_text().strip())
            os.kill(other, 0)
            print(f"another chain is running (pid {other}); not starting a second one")
            return 1
        except (ValueError, ProcessLookupError, PermissionError):
            pass
    PIDFILE.write_text(str(os.getpid()))
    try:
        log(f"chain started: wait for {FIRST}, then {' '.join(todo)}")
        wait_for(FIRST, "the run in progress")
        log(f"{FIRST}: complete, starting the chain")
        for chrom in todo:
            if not run_one(chrom):
                return 1
        log("chain done: every chromosome scored")
        return 0
    finally:
        PIDFILE.unlink(missing_ok=True)


if __name__ == "__main__":
    sys.exit(main())
