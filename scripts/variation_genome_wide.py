# SPDX-License-Identifier: AGPL-3.0-or-later
"""The human constraint axis over every chromosome's UNKNOWN blocks and scored elements.

Per chromosome, smallest first: gnomAD Gnocchi over the budget's blocks, the scored
elements (with Zoonomia phyloP for the uniform ones the constrained run never read)
and the VISTA elements, the coding and intron controls, saved as `variation_<chrom>.json`;
then the genome-wide summary. Resumable: a chromosome with a saved result is skipped;
one without a budget is reported and skipped.

    uv run python scripts/variation_genome_wide.py [--chroms chr21 chr22 ...]
"""

from __future__ import annotations

import argparse
import time

from genomeos.attribution.variation import distil, run_and_save
from genomeos.jobs import heartbeat
from genomeos.results import load_result, save_result

JOB = "variation_genome_wide"
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
    todo = [c for c in args.chroms if not load_result(f"variation_{c}")]
    print(f"{len(todo)} chromosomes to read: {' '.join(todo)}", flush=True)
    for chrom in todo:
        if not load_result(f"budget_{chrom}"):
            print(f"{chrom}: no budget result, skipped", flush=True)
            continue
        t0 = time.time()

        def progress(done: int, total: int, mb: int, chrom=chrom) -> None:
            heartbeat(JOB)

        out = run_and_save(chrom, progress=progress)
        heartbeat(JOB)
        ctrl = out["controls"]
        cds = (ctrl.get("canonical_cds") or {}).get("fraction_above")
        tiers = ", ".join(
            f"{t} {v['human_constrained_fraction']:.1%}"
            for t, v in out["by_tier"].items()
            if v["human_constrained_fraction"] is not None
        )
        ec = out["by_element_case"]["by_case"]
        print(
            f"{chrom}: {len(out['blocks'])} blocks, {len(out['elements'])} elements; coding {cds:.1%}"
            if cds is not None
            else f"{chrom}: {len(out['blocks'])} blocks; coding n/a",
            end="",
            flush=True,
        )
        print(
            f"; {tiers}; syntax elements name a gene "
            f"{ec['syntax']['name_a_gene_share'] if ec['syntax']['elements'] else 'n/a'} "
            f"({ec['syntax']['elements']}); {time.time() - t0:.0f} s",
            flush=True,
        )
    s = distil()
    save_result(JOB, s)
    print(
        f"genome-wide: {s['chromosomes']} chromosomes; coding "
        f"{s['controls']['canonical_cds']['fraction_above']}; "
        + "; ".join(f"{t} {v['human_constrained_fraction']}" for t, v in s["by_tier"].items())
        + f"; syntax elements name a gene {s['by_element_case']['syntax']['name_a_gene_share']} "
        f"({s['by_element_case']['syntax']['elements']}); {s['gnocchi_mb_fetched']} MB",
        flush=True,
    )


if __name__ == "__main__":
    main()
