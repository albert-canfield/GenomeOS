# SPDX-License-Identifier: AGPL-3.0-or-later
"""The CRISPRi result on the published pair sets, beside the published figures, and the coverage arms.

    uv run python scripts/crispri_published.py

Scores the frozen model (crispri.FEATURES, weights fitted on the covered K562 training pairs) by the
ENCODE benchmark's own estimator on all 10,356 training pairs (hold-one-chromosome-out) and all
4,378 held-out pairs (pooled, weighted by direct-effect probability), sets the result beside
Gschwind et al. 2026's Supplementary Table 3, and runs the three coverage arms on held-out K562.
The registration is crispri.PREREGISTERED_PUBLISHED, committed before this script existed.
No AlphaGenome request: the deletion values are the sweep's own cache.
"""

from __future__ import annotations

import csv
import gzip
import time
from collections import Counter

from genomeos import manifest as mf
from genomeos.attribution import crispri
from genomeos.attribution.measured import CRISPRI_SPLIT_OF
from genomeos.results import save_result


def invalid_connections(name: str) -> dict[str, int]:
    """The rows crispri.parse drops, counted by the benchmark's own ValidConnection reason."""
    with gzip.open(crispri.KNOWLEDGE / name, "rt") as fh:
        c = Counter(r["ValidConnection"] for r in csv.DictReader(fh, delimiter="\t"))
    c.pop("TRUE", None)
    return dict(sorted(c.items()))


def manifest(training: list[crispri.Pair], heldout: list[crispri.Pair]) -> dict:
    """The provenance contract (review item R9): what this result read, which release and which bytes."""
    inputs = [
        mf.input_entry(crispri.KNOWLEDGE / name, partition=CRISPRI_SPLIT_OF[name], pairs=len(pairs))
        for name, pairs in ((crispri.TRAINING, training), (crispri.HELDOUT, heldout))
    ]
    inputs.append(mf.input_entry(crispri.ELEMENTS, partition=None, role="DeletionTable: elements by overlap"))
    inputs.append(
        mf.input_entry(
            crispri.ELEMENT_CACHE, partition=None, role="ElementCache: every gene's deletion value"
        )
    )
    return {
        "sources": [
            {
                "accession": "EngreitzLab/CRISPR_comparison resources/crispr_data, EPCrisprBenchmark "
                "training_K562 and heldout_5_cell_types (Gschwind et al.)",
                "version": "main branch, unpinned upstream; fetched 2026-09-16; pinned here by sha256",
                "url": crispri.BASE_URL,
            },
            {
                "accession": "Gschwind et al. 2026, Nature, doi:10.1038/s41586-026-10781-4, "
                "Supplementary Table 3",
                "version": "as published (PMC13471189), read 2026-09-27; its AUPRCs and intervals are "
                "transcribed in genomeos/attribution/crispri.PUBLISHED, not read from a file",
                "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC13471189/",
            },
            {
                "accession": "ENCODE SCREEN cCREs scored by AlphaGenome deletion (all-element table and "
                "per-element response cache)",
                "version": "AlphaGenome as served during the 2026-09 all-element sweep (unpinned); "
                "pinned here by sha256",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "features": {k: list(v) for k, v in crispri.FEATURES.items()},
            "fit": "crispri.logistic_fit, lam 1e-3, 25 rounds, on the covered training pairs",
            "training_scoring": "hold-one-chromosome-out over all training pairs; uncovered pairs carry "
            "deletion features of zero",
            "heldout_scoring": "frozen weights, all held-out pairs pooled",
            "estimator": "crispri.benchmark_auprc (yardstick pr_curve, first and last points dropped, "
            "trapezoid)",
            "heldout_weight": "direct_vs_indirect_negative",
            "model_cells": list(crispri.MODEL_CELLS),
            "min_distance": crispri.MIN_DISTANCE,
            "reach": crispri.REACH,
            "bootstraps": crispri.BOOTSTRAPS,
            "bootstrap_seed": 0,
            "match_draws": crispri.MATCH_DRAWS,
            "request_budget": crispri.REQUEST_BUDGET,
            "alphagenome_requests": 0,
        },
        "exclusions": [
            {
                "file": name,
                "dropped_by_valid_connection": invalid_connections(name),
                "why": "the benchmark marks these pairs as no test of an enhancer-gene link",
            }
            for name in (crispri.TRAINING, crispri.HELDOUT)
        ]
        + [
            "held-out HCT116, Jurkat and WTC11 pairs carry no deletion value (the sweep kept four cell "
            "lines); they are scored with deletion features of zero, not dropped",
            "HCT116 as a second cell type: costed, not run (second_cell_type)",
            "per-cell held-out blocks for a cell with no regulated pair are omitted",
        ],
        "partitions": {
            crispri.TRAINING: CRISPRI_SPLIT_OF[crispri.TRAINING],
            crispri.HELDOUT: CRISPRI_SPLIT_OF[crispri.HELDOUT],
            "training": "fitted and scored by hold-one-chromosome-out (training_published_split)",
            "heldout": "evaluation only, frozen weights (heldout_published_pairs, "
            "coverage_arms_k562_heldout, post_hoc_positive_filter)",
        },
    }


def main() -> int:
    t0 = time.time()
    for name in (crispri.TRAINING, crispri.HELDOUT):
        crispri.fetch(name)
    training, heldout = crispri.load(crispri.TRAINING), crispri.load(crispri.HELDOUT)
    result = crispri.score_published(training, heldout, crispri.DeletionTable(), crispri.ElementCache())
    result[mf.KEY] = manifest(training, heldout)  # save_result takes it from the payload
    path = save_result("crispri_published", result)
    print(f"estimator check: {result['estimator_check_raw_distance']}")
    tr = result["training_published_split"]
    for name, m in tr["models"].items():
        print(f"  training LOCO, all pairs  {name:32s} AUPRC {m['auprc']}")
    print(f"  training gain {tr['deletion_gain']}; vs ENCODE-rE2G: {tr['against_encode_re2g']}")
    ho = result["heldout_published_pairs"]
    for name, m in ho["models"].items():
        print(f"  held out pooled weighted  {name:32s} AUPRC {m['auprc']}")
    print(f"  held-out gain {ho['deletion_gain']}; vs ENCODE-rE2G: {ho['against_encode_re2g']}")
    print(f"  baseline vs ABC: {ho['baseline_against_abc']}")
    for cell, c in ho["per_cell_type_weighted"].items():
        m = c["models"]
        print(
            f"  {cell:8s} weighted: a+d {m['activity + distance']['auprc']}  "
            f"+del {m['activity + distance + deletion']['auprc']}  gain {c['deletion_gain']}"
        )
    arms = result["coverage_arms_k562_heldout"]
    for k in ("arm1_all_pairs", "arm2_coverage_matched", "arm3_coverage_indicator_control"):
        print(f"  {k}: {arms[k]}")
    print(f"  invariant to coverage: {arms['invariant_to_coverage']}")
    print(f"second cell type: {result['second_cell_type']['cost']}")
    print(f"({time.time() - t0:.0f} s) -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
