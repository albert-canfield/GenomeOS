# SPDX-License-Identifier: AGPL-3.0-or-later
"""Predicted enhancer targets on every chromosome (AlphaGenome feature b), one after another.

Runs scripts/enhancer_targets.py per chromosome, smallest first, skipping the ones whose sample is
complete; a chromosome that fails (quota, network) is retried after a wait and the loop only ends when
every chromosome has its sample. Per-element answers are cached, so nothing is ever scored twice.

    uv run python scripts/enhancer_targets_genome_wide.py [--sample 200]
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

from genomeos.jobs import heartbeat
from genomeos.results import load_result

# smallest first, so the first results arrive early; chrM has no enhancers
ORDER = [
    "chr21",
    "chr22",
    "chrY",
    "chr19",
    "chr20",
    "chr18",
    "chr17",
    "chr16",
    "chr15",
    "chr14",
    "chr13",
    "chr12",
]
ORDER += ["chr11", "chr10", "chr9", "chr8", "chrX", "chr7", "chr6", "chr5", "chr4", "chr3", "chr2", "chr1"]


def scored(chrom: str) -> int:
    r = load_result(f"enhancer_targets_{chrom}")
    return int((r or {}).get("summary", {}).get("elements_scored", 0))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", type=int, default=200)
    args = ap.parse_args()
    waits = 0
    while True:
        todo = [c for c in ORDER if scored(c) < args.sample and (Path("data/reference") / f"{c}.fa").exists()]
        if not todo:
            done = [c for c in ORDER if scored(c) >= args.sample]
            print(f"done: {len(done)} chromosomes with {args.sample} elements each", flush=True)
            return 0
        chrom = todo[0]
        heartbeat("enhancer_targets_genome_wide")
        print(f"{chrom}: {scored(chrom)}/{args.sample} scored; {len(todo)} chromosomes to go", flush=True)
        rc = subprocess.call(  # noqa: S603
            [sys.executable, "scripts/enhancer_targets.py", "--chrom", chrom, "--sample", str(args.sample)]
        )
        if rc != 0 or scored(chrom) < args.sample:
            waits += 1
            wait = min(900, 60 * waits)
            print(f"{chrom}: not complete (exit {rc}); waiting {wait}s before the next attempt", flush=True)
            time.sleep(wait)
        else:
            waits = 0


if __name__ == "__main__":
    sys.exit(main())
