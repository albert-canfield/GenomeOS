# SPDX-License-Identifier: AGPL-3.0-or-later
"""Read one chromosome across the HPRC human assemblies (attribution/human_panel.py).

Streams the chromosome's Cactus alignment of 90 human assemblies from UCSC (3.3 GB for chr21, about
four minutes) into a small store under data/knowledge/human_panel (git-ignored), then the regions
around the calibration loci and the 69 syntax candidates, and writes
data/results/human_panel_<chrom>.json: every UNKNOWN block's presence, identity, structure and
recurring variation against a GC- and replication-timing-matched background, with coding exons and
the neutral tier as the controls; the fixed, storage and cannot-place catalogues; value domains
with gnomAD, TRExplorer, GTEx and MPRA beside them; the comparison with Gnocchi; the calibration
loci; the candidates; and what the genome would cost. No model is called.

    uv run python scripts/human_panel.py [--chrom chr21] [--skip-regions]
"""

from __future__ import annotations

import argparse
import json
import sys
import time

from genomeos.attribution import human_panel as hp
from genomeos.results import load_result, save_result


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--chrom", default="chr21")
    ap.add_argument(
        "--skip-regions", action="store_true", help="do not stream the calibration and candidate regions"
    )
    args = ap.parse_args(argv)
    t0 = time.time()

    def say(*a) -> None:
        print(f"[{time.time() - t0:6.0f} s]", *a, flush=True)

    store = hp.CACHE / args.chrom
    if not (store / "meta.json").exists():
        say(f"streaming {args.chrom}'s alignment")
        hp.distil_chromosome(
            args.chrom,
            progress=lambda d, t, b, s: (
                say(f"{d / 1e6:.0f} of {t / 1e6:.0f} MB, {b} blocks") if b % 2000 == 0 else None
            ),
        )
    if not args.skip_regions:
        wanted: dict[str, list[tuple[int, int]]] = {}
        for loc in hp.calibration_units():
            for s, e in loc["intervals"]:
                wanted.setdefault(loc["chrom"], []).append(
                    (max(0, s - hp.CALIBRATION_FLANK), e + hp.CALIBRATION_FLANK)
                )
        for c in (load_result("syntax_candidates_genome_wide") or {}).get("candidates", []):
            wanted.setdefault(c["chrom"], []).append(
                (max(0, c["start"] - hp.CALIBRATION_FLANK), c["end"] + hp.CALIBRATION_FLANK)
            )
        for chrom, ivs in sorted(wanted.items()):
            if not (hp.CACHE / f"{chrom}_regions" / "meta.json").exists():
                say(f"streaming {len(ivs)} regions of {chrom}")
                hp.distil_chromosome(chrom, wanted=ivs, name=f"{chrom}_regions")

    say("blocks, controls and units")
    result = hp.build(args.chrom, progress=say)
    edges = result["_edges"]
    say("storage catalogue")
    catalogue = hp.storage_catalogue(result, args.chrom)
    hp.save_catalogue(catalogue, args.chrom)
    say("against Gnocchi")
    gnocchi = hp.against_gnocchi(result, args.chrom)
    gnocchi["per_kilobase"]["poisson_null"] = hp.depletion_null(result)
    say("sensitivity")
    sens = hp.sensitivity(result)
    say("calibration loci")
    calibration = hp.calibrate(edges=edges)
    say("candidates")
    candidates = hp.read_candidates(edges=edges)
    say("genome-wide cost")
    cost = hp.genome_cost(result, catalogue)
    payload = hp.summarise(result, catalogue, gnocchi, sens, calibration, candidates, cost)
    payload["seconds"] = round(time.time() - t0, 1)
    path = save_result(f"human_panel_{args.chrom}", payload)
    say(f"saved {path}")
    json.dump(payload["control"], sys.stdout, indent=1)
    print()


if __name__ == "__main__":
    main()
