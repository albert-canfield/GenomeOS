# SPDX-License-Identifier: AGPL-3.0-or-later
"""The composition budget of every human chromosome: constraint over each UNKNOWN block.

Per chromosome, smallest first: read Zoonomia phyloP over the blocks with range
requests (about 3 bytes per base, nothing stored), fetch the 100-vertebrate
conserved elements once into the knowledge cache, save `budget_<chrom>.json`.
Then distil the genome-wide summary. Resumable: a chromosome with a saved budget
is skipped.

    uv run python scripts/budget_genome_wide.py [--chroms chr21 chr22 ...] [--threshold 2.27]
"""

from __future__ import annotations

import argparse
import time

from genomeos.attribution.budget import distil, run_and_save
from genomeos.jobs import heartbeat
from genomeos.results import load_result, save_result

JOB = "budget_genome_wide"
# smallest first so the table fills from the first hour
ORDER = ["chr21", "chr22", "chrY", "chr19", "chr20", "chr18", "chr17", "chr16", "chr15", "chr14", "chr13"]
ORDER += [
    "chr12",
    "chr11",
    "chr10",
    "chr9",
    "chr8",
    "chrX",
    "chr7",
    "chr6",
    "chr5",
    "chr4",
    "chr3",
    "chr2",
    "chr1",
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", nargs="*", default=ORDER)
    ap.add_argument("--threshold", type=float, default=2.27)
    args = ap.parse_args()
    todo = [c for c in args.chroms if not load_result(f"budget_{c}")]
    print(f"{len(todo)} chromosomes to budget: {' '.join(todo)}", flush=True)
    for chrom in todo:
        if not load_result(f"unknown_{chrom}"):
            print(f"{chrom}: no UNKNOWN result, skipped", flush=True)
            continue
        t0 = time.time()
        last = [0.0]

        def progress(done: int, total: int, mb: int, chrom=chrom, t0=t0, last=last) -> None:
            heartbeat(JOB)
            if time.time() - last[0] > 30:
                last[0] = time.time()
                print(
                    f"  {chrom}: {done}/{total} data blocks, {mb / 1e6:.0f} MB, {time.time() - t0:.0f} s",
                    flush=True,
                )

        out = run_and_save(chrom, threshold=args.threshold, progress=progress)
        heartbeat(JOB)
        tiers = ", ".join(f"{t} {v['fraction_of_unknown']:.1%}" for t, v in out["by_tier"].items())
        print(
            f"{chrom}: {len(out['blocks'])} blocks, constrained {out['constrained_fraction']:.2%} "
            f"of measured bases; {tiers}; "
            f"{out['cost'].get('phylop', {}).get('mb_fetched', 0)} MB in {out['cost']['seconds']} s",
            flush=True,
        )
    s = distil()
    save_result(JOB, s)
    print(
        f"genome-wide: {s['chromosomes']} chromosomes, {s['unknown_bp'] / 1e6:.0f} Mb of UNKNOWN blocks, "
        f"constrained {s['constrained_fraction']:.2%} of measured bases, "
        f"guessed at >= 0.5: {s['guessed_fraction']:.1%}",
        flush=True,
    )


if __name__ == "__main__":
    main()
