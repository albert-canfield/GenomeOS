# SPDX-License-Identifier: AGPL-3.0-or-later
"""JASPAR motifs over every chromosome's promoters and scored elements, then the operators per library.

    uv run python scripts/motifs_genome_wide.py [--chroms chr21 chr22 ...]

Resumable: a chromosome with a saved motifs_<chrom> result is skipped; the genome summary is redone.
"""

from __future__ import annotations

import argparse
import time

from genomeos.genome.motifs import distil, load_motifs, run_and_save
from genomeos.jobs import heartbeat
from genomeos.results import load_result, save_result

JOB = "motifs_genome_wide"
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
    args = ap.parse_args()
    motifs = load_motifs()
    todo = [c for c in args.chroms if not load_result(f"motifs_{c}")]
    print(f"{len(motifs)} profiles; {len(todo)} chromosomes to scan: {' '.join(todo)}", flush=True)
    for chrom in todo:
        t0 = time.time()
        try:
            out = run_and_save(chrom, motifs=motifs)
        except FileNotFoundError as e:
            print(f"{chrom}: skipped ({e})", flush=True)
            continue
        heartbeat(JOB)
        top = ", ".join(f"{x['factor']} {x['enrichment']}" for x in out["most_enriched"][:5])
        print(
            f"{chrom}: {out['promoters']} promoters, {len(out['elements'])} elements, "
            f"{out['requires_per_promoter']} requires per promoter; most enriched {top}; "
            f"{time.time() - t0:.0f} s",
            flush=True,
        )
    s = distil()
    save_result(JOB, s)
    print(
        f"genome-wide: {s['promoters']:,} promoters over {s['chromosomes']} chromosomes; "
        f"{len(s['libraries'])} libraries with factors; "
        f"{sum(1 for v in s['libraries'].values() if v['operators'])} with operator pairs",
        flush=True,
    )


if __name__ == "__main__":
    main()
