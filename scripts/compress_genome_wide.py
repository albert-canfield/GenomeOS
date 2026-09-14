# SPDX-License-Identifier: AGPL-3.0-or-later
"""Run the compression verdict pass over every chromosome, one at a time, then sum the genome.

    uv run python scripts/compress_genome_wide.py            # smallest first; resumes after a stop
    uv run python scripts/compress_genome_wide.py --dry-run  # say what it would do

Rules it keeps (Albert's standing policy, and the disk rule of 2026-09-14):
- one chromosome at a time, each in its own process, so memory is returned between chromosomes;
- a chromosome with a `compress_pass_<chrom>` result is skipped, so the chain resumes after a stop;
- free disk is checked before each chromosome against its expected temporaries (26 bytes a base)
  and every 30 s while it runs; if free space would fall below the floor (15 GB) the chromosome is
  interrupted, its temporaries are deleted and the chain stops and says so, rather than filling the
  disk under the other lanes. The pass checks the same floor itself before every step;
- a chromosome that fails for another reason is retried once, then left for the next start;
- when the loop ends, stopped or finished, the genome is summed from what exists
  (`compress_genome_wide`, with the missing chromosomes listed).
No model API is called. Progress goes to data/jobs/compress_genome_wide.log.
"""

from __future__ import annotations

import argparse
import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

from genomeos.attribution.compress import DISK_FLOOR_GB, KARYOTYPE, DiskGuard, rollup
from genomeos.jobs import heartbeat
from genomeos.results import save_result

JOB = "compress_genome_wide"
JOBS = Path("data/jobs")
RESULTS = Path("data/results")
SPILL = Path("data/cache/compress")
REFERENCE = Path("data/reference")
ORDER = [
    "chr21", "chr22", "chrY", "chr19", "chr20", "chr18", "chr17", "chr16", "chr15", "chr14", "chr13",
    "chr12", "chr11", "chr10", "chr9", "chr8", "chrX", "chr7", "chr6", "chr5", "chr4", "chr3", "chr2",
    "chr1",
]  # fmt: skip
BYTES_PER_BASE = 26  # adaptive and copy maps (10) and count rows waiting on disk (16)
POLL = 30


def log(msg: str) -> None:
    JOBS.mkdir(parents=True, exist_ok=True)
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    with open(JOBS / f"{JOB}.log", "a") as fh:
        fh.write(line + "\n")
    print(line, flush=True)


def alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    return True


def need_bytes(chrom: str) -> int:
    p = REFERENCE / f"{chrom}.fa"
    return BYTES_PER_BASE * (p.stat().st_size if p.exists() else 0)


def clean_spill() -> None:
    """Temporaries left by a pass that is no longer running."""
    if not SPILL.exists():
        return
    for d in SPILL.iterdir():
        pid = d.name.rsplit("_", 1)[-1]
        if d.is_dir() and not (pid.isdigit() and alive(int(pid))):
            shutil.rmtree(d, ignore_errors=True)


def run_one(chrom: str, train: str, guard: DiskGuard) -> int:
    """0 done, 3 stopped at the disk floor, anything else failed."""
    out = open(JOBS / f"compress_pass_{chrom}.log", "a")  # noqa: SIM115  (the child owns it)
    argv = [sys.executable, "scripts/compress.py", "--pass", f"{chrom}:{train}", "--heartbeat", JOB]
    proc = subprocess.Popen(argv, stdout=out, stderr=subprocess.STDOUT)  # noqa: S603
    try:
        while proc.poll() is None:
            time.sleep(POLL)
            heartbeat(JOB)
            try:
                guard.check()
            except Exception as e:  # DiskFloorError: stop the pass through its cleanup
                log(f"{chrom}: {e}; interrupting")
                proc.send_signal(signal.SIGTERM)
                proc.wait(timeout=600)
                return 3
    finally:
        out.close()
    return proc.returncode


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--floor-gb", type=float, default=DISK_FLOOR_GB)
    args = ap.parse_args()
    pidfile = JOBS / f"{JOB}.pid"
    if pidfile.exists() and pidfile.read_text().strip().isdigit() and alive(int(pidfile.read_text())):
        print(f"{JOB} is already running (pid {pidfile.read_text().strip()})")
        return
    guard = DiskGuard(SPILL, args.floor_gb)
    todo = [c for c in ORDER if not (RESULTS / f"compress_pass_{c}.json").exists()]
    if args.dry_run:
        for c in todo:
            print(
                c, f"needs about {need_bytes(c) / 1e9:.1f} GB of temporaries; {guard.free_gb():.1f} GB free"
            )
        return
    JOBS.mkdir(parents=True, exist_ok=True)
    pidfile.write_text(str(os.getpid()))
    clean_spill()
    log(f"start: {len(todo)} chromosomes to do, {guard.free_gb():.1f} GB free, floor {args.floor_gb:.0f} GB")
    stopped = None
    try:
        for chrom in todo:
            train = "chr21" if chrom == "chr22" else "chr22"
            try:
                free = guard.check(need_bytes(chrom))
            except Exception as e:
                stopped = f"{chrom} not started: {e}"
                log(stopped)
                break
            log(f"{chrom}: start (trained on {train}), {free:.1f} GB free")
            code = run_one(chrom, train, guard)
            if code not in (0, 3):
                log(f"{chrom}: exit {code}, retrying once")
                clean_spill()
                code = run_one(chrom, train, guard)
            clean_spill()
            if code == 3:
                stopped = f"{chrom}: stopped at the disk floor"
                log(stopped)
                break
            log(f"{chrom}: {'done' if code == 0 else f'failed (exit {code}), left for the next start'}")
    finally:
        genome = rollup(list(KARYOTYPE), RESULTS)
        genome["stopped"] = stopped
        path = save_result(JOB, genome, RESULTS)
        log(f"rollup of {len(genome['chromosomes'])} chromosomes: {path}; missing {genome['missing']}")
        pidfile.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
