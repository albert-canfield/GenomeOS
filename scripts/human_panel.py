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
        "--fill-pending", action="store_true", help="rerun only the stages a saved result records as pending"
    )
    ap.add_argument("--skip-calibration", action="store_true", help="refer to chr21's calibration loci")
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

    stages: dict[str, float] = {}

    def stage(name: str, fn, pending=None):
        say(name)
        t = time.time()
        try:
            return fn()
        except (OSError, TimeoutError) as e:  # a stalled stream is recorded, not waited on
            say(f"{name}: pending ({e})")
            return pending if pending is not None else {"pending": str(e)}
        finally:
            stages[name] = round(time.time() - t, 1)

    result = stage("blocks, controls and units", lambda: hp.build(args.chrom, progress=say))
    if args.fill_pending:
        fill_pending(args.chrom, result, stage, say)
        return
    edges = result["_edges"]
    catalogue = stage("storage catalogue", lambda: hp.storage_catalogue(result, args.chrom))
    if "_catalogue" in catalogue:
        hp.save_catalogue(catalogue, args.chrom)
    else:
        catalogue = {"summary": catalogue, "showcase": {}, "cost": {}}
    gnocchi = stage("against Gnocchi", lambda: hp.against_gnocchi(result, args.chrom))
    if "per_kilobase" in gnocchi:
        gnocchi["per_kilobase"]["poisson_null"] = hp.depletion_null(result)
    sens = stage("sensitivity", lambda: hp.sensitivity(result))
    if args.skip_calibration:
        calibration = [
            {"see": "human_panel_chr21", "note": "the calibration loci do not depend on the chromosome read"}
        ]
    else:
        calibration = stage("calibration loci", lambda: hp.calibrate(edges=edges), pending=[])
    candidates = stage("candidates", lambda: hp.read_candidates(edges=edges))
    local = stage("this chromosome's real unknown", lambda: hp.local_candidates(result, args.chrom))
    cost = stage("genome-wide cost", lambda: hp.genome_cost(result, catalogue))
    payload = hp.summarise(result, catalogue, gnocchi, sens, calibration, candidates, cost)
    payload["local_candidates"] = local
    payload["stage_seconds"] = stages
    payload["seconds"] = round(time.time() - t0, 1)
    path = save_result(f"human_panel_{args.chrom}", payload)
    say(f"saved {path}")
    json.dump(payload["control"], sys.stdout, indent=1)
    print()


def fill_pending(chrom: str, result: dict, stage, say) -> None:
    """Rerun the network stages a saved result left pending and merge them in, nothing else."""
    payload = load_result(f"human_panel_{chrom}")
    filled = []
    if "pending" in (payload.get("storage") or {}):
        catalogue = stage("storage catalogue", lambda: hp.storage_catalogue(result, chrom))
        if "_catalogue" in catalogue:
            hp.save_catalogue(catalogue, chrom)
            payload["storage"] = catalogue["summary"]
            payload["storage_showcase"] = catalogue["showcase"]
            payload["cost"]["storage_tracks"] = catalogue["cost"]
            filled.append("storage")
    if "pending" in (payload.get("against_gnocchi") or {}):
        gnocchi = stage("against Gnocchi", lambda: hp.against_gnocchi(result, chrom))
        if "per_kilobase" in gnocchi:
            gnocchi["per_kilobase"]["poisson_null"] = hp.depletion_null(result)
            payload["against_gnocchi"] = {k: v for k, v in gnocchi.items() if not k.startswith("_")}
            payload["cost"]["gnocchi"] = gnocchi.get("cost")
            filled.append("against_gnocchi")
    payload.setdefault("filled_later", []).extend(filled)
    save_result(f"human_panel_{chrom}", {k: v for k, v in payload.items() if k not in ("result", "date")})
    say(f"filled {filled or 'nothing'}")


if __name__ == "__main__":
    main()
