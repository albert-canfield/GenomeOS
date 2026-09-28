# SPDX-License-Identifier: AGPL-3.0-or-later
"""Run R8's registered coupling pretest once and save `joint_pretest`.

    uv run python scripts/joint_pretest.py

Reads the CRISPRi training file only (the held-out file is never opened), the all-element deletion
table and the per-element response cache one chromosome at a time. 0 AlphaGenome requests.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri, measured  # noqa: E402
from genomeos.attribution import joint_pretest as jp  # noqa: E402
from genomeos.results import save_result  # noqa: E402


def manifest() -> dict:
    train = [n for n in measured.CRISPRI_FILES if measured.CRISPRI_SPLIT_OF[n] == measured.TRAINING]
    inputs = [mf.input_entry(measured.CRISPRI_KNOWLEDGE / n, partition=measured.TRAINING) for n in train]
    inputs.append(mf.input_entry(jp.ALL_ELEMENTS, partition=None))
    inputs.append(mf.input_entry(jp.ELEMENT_CACHE, partition=None))
    return {
        "sources": [
            {
                "accession": (
                    "EngreitzLab/CRISPR_comparison EPCrisprBenchmark (Gschwind et al. 2025), training file"
                ),
                "version": "main, as fetched; pinned by sha256",
            },
            {
                "accession": (
                    "this repository, the all-enhancer AlphaGenome deletion sweep and its per-element cache"
                ),
                "version": "pinned by sha256",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "registered": jp.REGISTERED,
            "tau": jp.TAU,
            "variants": list(jp.VARIANTS),
            "control": jp.CONTROL,
            "bootstraps": jp.BOOTSTRAPS,
            "seed": jp.SEED,
            "interval_quantiles": list(jp.INTERVAL),
            "ridge": jp.RIDGE,
            "folds": "leave-one-chromosome-out within the training split",
            "cells_scored": list(jp.CELLS),
            "overlap_reach": crispri.REACH,
            "well_powered": measured.WELL_POWERED,
            "pass_rule": jp.PASS_RULE,
            "falsifier": jp.FALSIFIER,
        },
        "exclusions": [
            "the held-out CRISPRi file: never opened in this lane (R5)",
            "significant increases and underpowered or power-unknown nulls: counted, never scored (R2)",
            "pairs no cached element answers: missing, never imputed",
            "pairs the benchmark marks as not a valid connection: never read",
        ],
        "partitions": {measured.TRAINING: "the CRISPRi benchmark's training file (K562), the only one read"},
    }


def main() -> None:
    t0 = time.time()
    pairs = jp.load_training()
    rows = jp.score(pairs)
    pooled = jp.evaluate(rows)
    out = {
        "question": (
            "does a coupling term (an element's candidate genes competing, or a gene's share of its "
            "regulators) beat independent per-element deletion scoring at naming CRISPRi-measured targets?"
        ),
        "registration": {
            "date": jp.REGISTERED,
            "pass_rule": jp.PASS_RULE,
            "falsifier": jp.FALSIFIER,
            "readings": jp.READINGS,
        },
        "counts": jp.counts(rows),
        "pooled": pooled,
        "per_cell": jp.per_cell(rows, pooled),
        "heldout_overlap_note": (
            "lane-split's split_overlap found 249 held-out pairs sharing bases with training intervals; "
            "irrelevant here, since the held-out file is never read"
        ),
        "alphagenome_requests": 0,
        "seconds": round(time.time() - t0, 1),
    }
    p = save_result("joint_pretest", out, manifest=manifest())
    print(p)
    print(
        {
            k: pooled[k]
            for k in ("pairs_scored", "positives", "negatives", "components", "auprc", "reading", "r8")
        }
    )
    for name, g in pooled["gain_over_independent"].items():
        print(name, g["by_component"], g["by_chromosome"])


if __name__ == "__main__":
    main()
