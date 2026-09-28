# SPDX-License-Identifier: AGPL-3.0-or-later
"""Score the attribution layer once on the independent IGVF MHC CRISPRi screen (attribution/indep.py).

    uv run python scripts/indep_benchmark.py            # the registered test, once
    uv run python scripts/indep_benchmark.py coverage   # label-blind coverage only, no result written

Inputs (git-ignored, data/cache/indep): IGVFFI4093WUVB.tsv.gz (the screen) and scE2G_K562_chr6_MHC.tsv.gz,
the chr6 27-36 Mb rows of IGVFFI1706PNVV. Both are fetched by `fetch()` and checked against the portal's
md5. No AlphaGenome request is made: the deletion is read from the stored sweep and its per-element cache,
one chromosome (chr6) at a time.
"""

from __future__ import annotations

import gzip
import hashlib
import subprocess
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri, indep  # noqa: E402
from genomeos.results import save_result  # noqa: E402

NAME = "indep_mhc_crispri"


def _md5(p: Path) -> str:
    h = hashlib.md5()
    with p.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def fetch() -> None:
    indep.CACHE.mkdir(parents=True, exist_ok=True)
    screen = indep.CACHE / indep.SCREEN_FILE
    if not screen.exists():
        subprocess.run(["curl", "-sfL", "-o", str(screen), indep.SOURCE["url"]], check=True)
    if _md5(screen) != indep.SOURCE["md5"]:
        raise SystemExit(f"{screen}: md5 does not match the portal's")
    distilled = indep.CACHE / indep.SCE2G_FILE
    if not distilled.exists():
        full = indep.CACHE / "IGVFFI1706PNVV.tsv.gz"
        if not full.exists():
            subprocess.run(["curl", "-sfL", "-o", str(full), indep.PUBLISHED_MODEL["url"]], check=True)
        if _md5(full) != indep.PUBLISHED_MODEL["md5"]:
            raise SystemExit(f"{full}: md5 does not match the portal's")
        with gzip.open(full, "rt") as src, gzip.open(distilled, "wt") as dst:
            for line in src:
                if line.startswith("#") or line.startswith("ElementChr"):
                    dst.write(line)
                    continue
                f = line.split("\t", 3)
                if f[0] == "chr6" and int(f[1]) > 27_000_000 and int(f[2]) < 36_000_000:
                    dst.write(line)


def _ap(pairs: list[dict[str, Any]], key: str) -> float | None:
    v = indep.ap_of(pairs, key)
    return None if v is None else round(v, 4)


def _base(pairs: list[dict[str, Any]]) -> float | None:
    return round(sum(p["label"] == "decrease" for p in pairs) / len(pairs), 4) if pairs else None


def _block(pairs: list[dict[str, Any]], keys: list[str], unit: str) -> dict[str, Any]:
    draws = indep.bootstrap(pairs, keys, unit=unit)
    return {
        "pairs": len(pairs),
        "positives": sum(p["label"] == "decrease" for p in pairs),
        "units": len({p[unit] for p in pairs}),
        "unit": unit,
        "base_rate": _base(pairs),
        "ap": {k: _ap(pairs, k) for k in keys},
        "ci95": {k: indep.interval(v) for k, v in draws.items()},
        "draws_kept": len(draws[keys[0]]),
    }


def main(mode: str = "score") -> None:
    fetch()
    rows = indep.load_screen()
    genes = indep.gencode()
    reused = indep.benchmark_pairs()
    sce2g = indep.load_sce2g()
    pairs, excluded = indep.prepare(
        rows, genes, reused, crispri.DeletionTable(), crispri.ElementCache(), sce2g
    )
    reach = [p for p in pairs if p["in_reach"]]
    coverage = {
        "rows": len(rows),
        "excluded": excluded,
        "eligible_pairs": len(pairs),
        "out_of_reach_pairs": len(pairs) - len(reach),
        "in_reach_pairs": len(reach),
        "in_reach_overlapping_a_sweep_element": sum(p["overlaps_sweep"] for p in reach),
        "in_reach_covered": sum(p["covered"] for p in reach),
        "in_reach_scored_by_sce2g": sum(p["sce2g"] is not None for p in reach),
    }
    if mode == "coverage":
        print(coverage)
        return

    for p in pairs:
        p["label"] = indep.outcome(p)
    increases = [p for p in reach if p["label"] == "increase"]
    test = [p for p in reach if p["label"] != "increase"]
    covered = [p for p in test if p["covered"]]
    uncovered = [p for p in test if not p["covered"]]
    positives = sum(p["label"] == "decrease" for p in covered)
    share = len(covered) / len(test) if test else 0.0
    readable = positives >= indep.MIN_POSITIVES and share >= indep.MIN_COVERAGE

    keys = ["deletion", "distance", "annotation", "node"]
    primary = _block(covered, keys, "element") if readable else None
    by_gene = _block(covered, keys, "gene") if readable else None
    both = [p for p in covered if p["sce2g"] is not None]
    published = _block(both, ["deletion", "sce2g", "distance"], "element") if readable and both else None
    ver = indep.verdict(
        primary["ci95"]["deletion-distance"] if primary else None,
        primary["ci95"]["deletion"] if primary else None,
        primary["base_rate"] if primary else 1.0,
        readable,
    )
    kinds = ("decrease", "increase", "not_detected")
    payload = {
        "partition": "independent",
        "source_accession": indep.SOURCE["accession"],
        "question": (
            "Does the stored AlphaGenome deletion sweep rank measured CRISPRi targets above distance on a "
            "screen no script in this repository had read (IGVF K562 MHC CRISPRi Perturb-seq)?"
        ),
        "verdict": ver,
        "preregistered": indep.PREREGISTERED,
        "alphagenome_requests": 0,
        "coverage": {**coverage, "in_reach_test_pairs": len(test), "covered_share": round(share, 4)},
        "outcomes_eligible": {k: sum(p["label"] == k for p in pairs) for k in kinds},
        "outcomes_in_reach": {k: sum(p["label"] == k for p in reach) for k in kinds},
        "primary_element_bootstrap": primary,
        "sensitivity_gene_bootstrap": by_gene,
        "against_published_model_on_shared_pairs": published,
        "selection_check_distance_ap": {
            "covered": {"pairs": len(covered), "ap": _ap(covered, "distance"), "base_rate": _base(covered)},
            "uncovered_in_reach": {
                "pairs": len(uncovered),
                "ap": _ap(uncovered, "distance"),
                "base_rate": _base(uncovered),
            },
        },
        "increases_in_reach": {
            "pairs": len(increases),
            "covered": sum(p["covered"] for p in increases),
            "covered_with_no_predicted_fall": sum(
                1 for p in increases if p["covered"] and p["deletion"] == 0.0
            ),
        },
        "evidence": (
            "experimental: IGVF K562 CRISPRi Perturb-seq of the MHC (IGVFFI4093WUVB); predicted: AlphaGenome "
            "deletion sweep, stored; baseline: scE2G v1.2 K562 (IGVFFI1706PNVV)"
        ),
    }
    manifest = {
        "sources": [
            {**indep.SOURCE, "version": f"released {indep.SOURCE['released']}"},
            {**indep.PUBLISHED_MODEL, "version": "scE2G v1.2.0"},
            {
                "accession": "EngreitzLab/CRISPR_comparison EPCrisprBenchmark",
                "version": "main (sha256 pinned)",
            },
            {"accession": "GENCODE", "version": "v50, chr6"},
            {"accession": "AlphaGenome all-element deletion sweep", "version": "stored, 2026-09-13..16"},
        ],
        "inputs": [
            mf.input_entry(indep.CACHE / indep.SCREEN_FILE, partition="independent"),
            mf.input_entry(indep.CACHE / indep.SCE2G_FILE, partition="n/a: published baseline"),
            mf.input_entry(crispri.KNOWLEDGE / crispri.TRAINING, partition="training (overlap check only)"),
            mf.input_entry(crispri.KNOWLEDGE / crispri.HELDOUT, partition="heldout (overlap check only)"),
            mf.input_entry(indep.REFERENCE / "gencode_v50_chr6.gff3.gz", partition=None),
            mf.input_entry(crispri.ELEMENTS / "chr6.json", partition=None),
            mf.input_entry(crispri.ELEMENT_CACHE / "chr6.json.gz", partition=None),
        ],
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "cell": indep.CELL,
            "fdr": indep.FDR,
            "promoter_window": indep.PROMOTER_WINDOW,
            "reach": indep.REACH,
            "min_positives": indep.MIN_POSITIVES,
            "min_coverage": indep.MIN_COVERAGE,
            "bootstraps": indep.BOOTSTRAPS,
            "seed": indep.SEED,
            "min_distance": crispri.MIN_DISTANCE,
        },
        "exclusions": [{"reason": k, "pairs": v} for k, v in excluded.items()]
        + [
            {"reason": "beyond the sweep's 500 kb reach", "pairs": len(pairs) - len(reach)},
            {"reason": "significant increase, reported apart", "pairs": len(increases)},
            {"reason": "in reach, not covered by the sweep (missing, not imputed)", "pairs": len(uncovered)},
        ],
        "partitions": {
            "independent": "IGVF MHC CRISPRi screen, read by no earlier script; scored once",
            "training": "ENCODE K562 training pairs, read only to exclude overlapping pairs",
            "heldout": "ENCODE held-out pairs, read only to exclude overlapping pairs",
        },
    }
    print(save_result(NAME, payload, manifest=manifest))
    print("verdict:", ver)
    print({"primary": primary, "published": published})


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "score")
