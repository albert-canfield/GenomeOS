# SPDX-License-Identifier: AGPL-3.0-or-later
"""Origin per gene, age per library, paralogues: Ensembl Compara for human, streamed once.

Two passes, both resumable. First every Ensembl vertebrate species is placed on the clade
ladder shared with human through the taxonomy endpoint (one call per species, cached under
data/knowledge/homology). Then the human homology dump (109 MB gzipped) is streamed over
HTTP, never stored, and distilled to `origin_genome_wide.json`: origin and paralogues per
protein-coding gene, counts per stratum, the age distribution per BioLib library.

    uv run python scripts/origin_genome_wide.py [--classify-only] [--release 116]
"""

from __future__ import annotations

import argparse
import time

from genomeos.jobs import heartbeat
from genomeos.knowledge.homology import DUMP_URL, classify_species, ensembl_species, run_and_save

JOB = "origin_genome_wide"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--classify-only", action="store_true")
    ap.add_argument("--release", type=int, default=116)
    args = ap.parse_args()
    species = ensembl_species()
    print(f"{len(species)} Ensembl vertebrate species to place on the ladder", flush=True)

    def sp_progress(i: int, n: int, sp: str, stratum: str | None) -> None:
        heartbeat(JOB)
        print(f"  {i}/{n} {sp}: {stratum}", flush=True)

    placed = classify_species(species, progress=sp_progress)
    missing = sorted(sp for sp, v in placed.items() if not v.get("stratum"))
    print(f"placed {len(placed) - len(missing)} species; unplaced: {' '.join(missing) or 'none'}", flush=True)
    if args.classify_only:
        return
    t0 = time.time()

    def progress(rows: int, genes: int, seconds: float) -> None:
        heartbeat(JOB)
        print(f"  {rows / 1e6:.1f} M rows, {genes:,} genes, {seconds:.0f} s", flush=True)

    out = run_and_save(DUMP_URL.format(release=args.release), progress=progress)
    heartbeat(JOB)
    rows = sum(c.get("rows", 0) for c in out["cost"].values())
    print(
        f"origin: {out['coding_genes']:,} coding genes; by ladder {out['by_ladder']}; "
        f"{len(out['libraries'])} libraries aged; {rows / 1e6:.1f} M rows over {len(out['cost'])} dumps in "
        f"{time.time() - t0:.0f} s",
        flush=True,
    )


if __name__ == "__main__":
    main()
