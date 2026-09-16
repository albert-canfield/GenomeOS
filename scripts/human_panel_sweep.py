# SPDX-License-Identifier: AGPL-3.0-or-later
"""Read the chromosomes the human panel has not read, one at a time, unattended.

    uv run python scripts/human_panel_sweep.py            # every unread chromosome, largest first
    uv run python scripts/human_panel_sweep.py --dry-run  # say what it would do

The panel (scripts/human_panel.py) reads one chromosome per invocation, and for three days it was
driven by hand: a session started a chromosome, watched the disk, committed the result and started the
next. Twice a session ended mid-chromosome and nothing picked the sweep up, and once the running job
was invisible in the Progress tab because nothing in the registry knew about it. This sequences it the
way scripts/enhancer_targets_all_chain.py sequences the deletion scoring, so the Progress tab can start
it, show it and heal it.

Rules it keeps:
- one chromosome at a time: the panel holds a whole chromosome's alignment store in memory (chr1 and
  chr2 need 5 to 6 GB), so two at once is how the machine runs out;
- a chromosome whose result exists is skipped, so the sweep resumes after any restart;
- it refuses to start a chromosome with less than FLOOR_GB free and stops rather than filling the
  disk, because the disk floor stopped this lane four times on 2026-09-15;
- it streams no model requests and takes no key: the panel calls no model at all;
- it writes a heartbeat, so `jobs.last_activity` can tell a stalled run from a slow one.

Progress goes to data/jobs/human_panel_sweep.log. Results are NOT committed here: a panel result is
read before it is published (seven of 23 chromosomes fail a control), and that reading is a person's.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path

from genomeos.jobs import heartbeat

JOB = "human_panel_sweep"
RESULTS = Path("data/results")
JOBS = Path("data/jobs")
FLOOR_GB = 20.0

# Largest first: the large chromosomes are the ones whose absence a genome-wide claim is judged on,
# and they are the ones a disk stop keeps taking out. chrM has no panel; chrY is read like the rest.
ORDER = [f"chr{n}" for n in range(1, 23)] + ["chrX", "chrY"]


def say(message: str) -> None:
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {message}"
    print(line, flush=True)
    JOBS.mkdir(parents=True, exist_ok=True)
    with open(JOBS / f"{JOB}.log", "a") as fh:
        fh.write(line + "\n")
    heartbeat(JOB)


def free_gb() -> float:
    return shutil.disk_usage(".").free / 1e9


def unread() -> list[str]:
    """The chromosomes with no committed panel result, in the order they are to be read."""
    return [c for c in ORDER if not (RESULTS / f"human_panel_{c}.json").exists()]


def read_one(chrom: str) -> int:
    """Run the panel on one chromosome, its output going to that chromosome's own log."""
    log = JOBS / f"human_panel_{chrom}.log"
    with open(log, "a") as fh:
        proc = subprocess.Popen(  # noqa: S603
            [sys.executable, "scripts/human_panel.py", "--chrom", chrom], stdout=fh, stderr=subprocess.STDOUT
        )
        while proc.poll() is None:
            time.sleep(30)
            heartbeat(JOB)
            if free_gb() < FLOOR_GB / 2:
                proc.terminate()
                say(f"{chrom}: stopped, free space fell to {free_gb():.1f} GB")
                return 1
    return proc.returncode or 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--floor-gb", type=float, default=FLOOR_GB, help="refuse to start below this")
    args = ap.parse_args(argv)

    todo = unread()
    if not todo:
        say(f"every chromosome read: {len(ORDER)} results in {RESULTS}")
        return 0
    say(f"unread: {' '.join(todo)} ({free_gb():.0f} GB free)")
    if args.dry_run:
        return 0

    for chrom in todo:
        if free_gb() < args.floor_gb:
            say(f"stopping before {chrom}: {free_gb():.1f} GB free, floor {args.floor_gb:.0f} GB")
            return 0
        say(f"{chrom}: starting ({free_gb():.0f} GB free)")
        code = read_one(chrom)
        result = RESULTS / f"human_panel_{chrom}.json"
        if code == 0 and result.exists():
            say(f"{chrom}: read, {result.stat().st_size / 1e6:.1f} MB, uncommitted and waiting to be read")
        else:
            say(f"{chrom}: ended with code {code} and no result; stopping rather than running on")
            return code or 1
    say("sweep done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
