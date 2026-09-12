# SPDX-License-Identifier: AGPL-3.0-or-later
"""Run scripts/enhancer_targets_all.py over every chromosome, one after another.

Albert's instruction of 2026-09-13: let the running chr21 job finish, then continue with chr1, chr2,
chr3 and on until every chromosome is scored. This driver changes nothing in the scoring script (owned
by the area I decoding lane); it only sequences it.

    uv run python scripts/enhancer_targets_all_chain.py            # wait for chr21, then chr1..chrY
    uv run python scripts/enhancer_targets_all_chain.py --dry-run  # say what it would do

Rules it keeps:
- one chromosome at a time: every run spends the same AlphaGenome quota, so two in parallel only race;
- a chromosome whose result says `"complete": true` is skipped, so the chain resumes after any restart
  (the scoring script itself resumes element by element from its cache);
- a chromosome another process is already scoring is waited on, never started twice;
- a run that ends without completing is retried after a pause; AlphaGenome being disabled stops the chain.
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

from genomeos.jobs import heartbeat

JOB = "enhancer_targets_all_chain"
RESULTS = Path("data/results")
JOBS = Path("data/jobs")
FIRST = "chr21"  # already running when the chain was asked for
ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
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
    out = subprocess.run(
        ["pgrep", "-f", f"enhancer_targets_all.py --chrom {chrom}( |$)"], capture_output=True, text=True
    )
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


def run_one(chrom: str) -> bool:
    """Score one chromosome to completion, retrying runs that end early. False stops the chain."""
    for attempt in range(1, MAX_RETRIES + 1):
        if complete(chrom):
            return True
        if running(chrom):
            wait_for(chrom, "already running elsewhere")
            continue
        log(f"{chrom}: starting (attempt {attempt})")
        with open(JOBS / f"enhancer_targets_all_{chrom}.log", "a") as out:
            proc = subprocess.Popen(
                [sys.executable, "scripts/enhancer_targets_all.py", "--chrom", chrom],
                stdout=out,
                stderr=subprocess.STDOUT,
            )
            while proc.poll() is None:
                heartbeat(JOB)
                time.sleep(POLL)
        if proc.returncode == 2:
            log(f"{chrom}: AlphaGenome is disabled (exit 2); the chain stops here")
            return False
        if complete(chrom):
            sc = scored(chrom)
            log(f"{chrom}: complete, {sc[0]:,} elements" if sc else f"{chrom}: complete")
            return True
        pause = RETRY_PAUSE // 60
        log(f"{chrom}: run ended with code {proc.returncode} before completing; retry in {pause} min")
        time.sleep(RETRY_PAUSE)
    log(f"{chrom}: not complete after {MAX_RETRIES} attempts; the chain stops here")
    return False


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
