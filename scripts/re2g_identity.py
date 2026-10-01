# SPDX-License-Identifier: AGPL-3.0-or-later
"""Is the model this lane scored the frozen headline model? (data/results/re2g_identity.json)

    uv run --frozen python scripts/re2g_identity.py

Three checks, because the question has three parts and only the third can be done pair by pair against
the committed file:

1. **Pair by pair, two independent derivations.** The per-pair deletion value and model score from this
   lane's own path (`re2g_paired.our_scores`) against a second derivation that calls the library's steps
   directly in the order `crispri.score_published` calls them. Reports pairs compared and pairs
   differing, on the 1,918 K562 primary and on all 4,378.
2. **Against the committed headline.** `data/results/crispri_published.json` carries aggregate AUPRCs
   and no per-pair column, so a pair-level diff against the file is not possible; that is stated rather
   than worked around. What is compared is the weighted AUPRC of each model on each population, at the
   four decimals the file records.
3. **The two coverage rules reconciled**, by name, with the count each gives and which one the frozen
   model actually applied.

Reads only files already on disk. 0 model requests.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "re2g_identity"
HEADLINE = Path("data/results/crispri_published.json")
PAIRED = Path("data/results/re2g_paired.json")
TOLERANCE = 1e-12  # two derivations of the same arithmetic, so exact equality is the expectation


def _paired_module() -> Any:
    spec = importlib.util.spec_from_file_location(
        "re2g_paired", Path(__file__).resolve().parent / "re2g_paired.py"
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def second_derivation(
    training: list[crispri.Pair], heldout: list[crispri.Pair]
) -> tuple[dict[str, list[float]], list[tuple[float, float]]]:
    """The headline's own steps, called directly: annotate, fit on covered training pairs, score.

    This is `crispri.score_published`'s recipe for the held-out block, run here so the comparison is
    between two derivations rather than between a derivation and itself.
    """
    crispri.annotate(training + heldout, crispri.DeletionTable(), crispri.ElementCache())
    for p in training + heldout:
        p.features["covered"] = float(p.covered)
    train = [p for p in training if p.covered]
    weights = {
        name: crispri.logistic_fit(crispri.matrix(train, cols), [p.regulated for p in train])
        for name, cols in crispri.FEATURES.items()
    }
    scores = {
        name: crispri.logistic_score(weights[name], crispri.matrix(heldout, cols))
        for name, cols in crispri.FEATURES.items()
    }
    features = [(p.features["top_target"], p.features["deletion_drop"]) for p in heldout]
    return scores, features


def main() -> None:
    paired = _paired_module()
    headline = json.loads(HEADLINE.read_text())
    held_block = headline["heldout_published_pairs"]

    # derivation A: this lane's path
    training_a, heldout_a = crispri.load(crispri.TRAINING), crispri.load(crispri.HELDOUT)
    scores_a, _ = paired.our_scores(training_a, heldout_a)
    features_a = [(p.features["top_target"], p.features["deletion_drop"]) for p in heldout_a]

    # derivation B: the library's steps in the headline's order, on freshly loaded pairs
    training_b, heldout_b = crispri.load(crispri.TRAINING), crispri.load(crispri.HELDOUT)
    scores_b, features_b = second_derivation(training_b, heldout_b)

    model = "activity + distance + deletion"
    k562 = [i for i, p in enumerate(heldout_a) if p.cell == "K562"]

    def differing(idx: list[int]) -> dict[str, int]:
        feat = sum(
            1
            for i in idx
            if abs(features_a[i][0] - features_b[i][0]) > TOLERANCE
            or abs(features_a[i][1] - features_b[i][1]) > TOLERANCE
        )
        score = sum(1 for i in idx if abs(scores_a[model][i] - scores_b[model][i]) > TOLERANCE)
        return {"pairs_compared": len(idx), "deletion_values_differing": feat, "scores_differing": score}

    pair_level = {
        "primary_k562": differing(k562),
        "all_heldout": differing(list(range(len(heldout_a)))),
        "what_was_compared": (
            "the per-pair top_target and deletion_drop columns and the per-pair model score, between "
            "this lane's derivation and a second one that calls the library's steps in the order "
            f"crispri.score_published calls them; equal within {TOLERANCE}"
        ),
    }

    def auprc(idx: list[int], name: str) -> float | None:
        rows = [heldout_a[i] for i in idx]
        return crispri.bench_metrics([scores_a[name][i] for i in idx], rows, weighted=True)["auprc"]

    against_headline = {
        "note": (
            "data/results/crispri_published.json records aggregate AUPRCs and no per-pair column, so a "
            "pair-level diff against the committed file is not possible. This compares the weighted "
            "AUPRC of each model on each population at the four decimals the file records"
        ),
        "pooled_4378": {
            name: {
                "headline": held_block["models"][name]["auprc"],
                "ours": auprc(list(range(len(heldout_a))), name),
            }
            for name in crispri.FEATURES
        },
        "k562_1918": {
            name: {
                "headline": held_block["per_cell_type_weighted"]["K562"]["models"][name]["auprc"],
                "ours": auprc(k562, name),
            }
            for name in crispri.FEATURES
        },
    }
    for block in ("pooled_4378", "k562_1918"):
        for name, v in against_headline[block].items():
            v["equal"] = v["headline"] == v["ours"]
    against_headline["all_equal"] = all(
        v["equal"] for block in ("pooled_4378", "k562_1918") for v in against_headline[block].values()
    )

    answered = {
        cell: sum(
            1
            for p in heldout_a
            if p.cell == cell and p.features.get("deletion_answered") == 1.0
        )
        for cell in sorted({p.cell for p in heldout_a})
    }
    in_model_cells = [p for p in heldout_a if p.cell in crispri.MODEL_CELLS]
    reconciliation = {
        "rule_1_cell_membership": {
            "name": "cell membership (crispri.annotate's `cells` argument, default crispri.MODEL_CELLS)",
            "what_it_decides": (
                "whether the two deletion columns are computed for a pair at all. A pair whose own cell "
                "line is not in MODEL_CELLS takes the branch (0.0, []) and both columns are zero"
            ),
            "pairs_with": len(in_model_cells),
            "pairs_without": len(heldout_a) - len(in_model_cells),
            "headline_reports_without": held_block["pairs_without_a_deletion_value"],
            "agrees_with_the_headline": (
                len(heldout_a) - len(in_model_cells) == held_block["pairs_without_a_deletion_value"]
            ),
        },
        "rule_2_answered": {
            "name": "answered (the non-model column `deletion_answered`, set when the sweep returned "
            "any value for this gene in this pair's own cell)",
            "what_it_decides": (
                "whether the sweep actually had something to say about the pair, so that a drop of zero "
                "can be told apart from a table that was silent"
            ),
            "pairs_answered": sum(answered.values()),
            "by_cell": answered,
        },
        "the_gap": {
            "pairs": len(in_model_cells) - sum(answered.values()),
            "what_they_are": (
                "pairs in a model cell line for which the deletion columns were computed but the sweep "
                "returned no value, so top_target and deletion_drop are both 0.0. The headline's "
                "cell-membership count includes them among the pairs 'with a deletion value'; the "
                "answered count does not"
            ),
            "by_cell": {
                cell: sum(1 for p in heldout_a if p.cell == cell) - answered.get(cell, 0)
                for cell in crispri.MODEL_CELLS
                if any(p.cell == cell for p in heldout_a)
            },
        },
        "which_rule_the_frozen_model_applied": (
            "cell membership. The frozen model's two columns exist for every pair in a model cell line "
            "and are zero elsewhere; `deletion_answered` is not a model column and no scorer reads it. "
            "So this lane's 'ours' IS the frozen headline model, not a different feature join, which "
            "the AUPRC comparison above checks at four decimals on both populations rather than asserts"
        ),
        "where_the_feature_is_actually_non_zero": {
            "pairs_with_a_predicted_drop_above_zero": sum(
                1 for p in heldout_a if p.features["deletion_drop"] > 0
            ),
            "pairs_with_the_top_target_flag": sum(1 for p in heldout_a if p.features["top_target"] > 0),
            "pairs_where_both_columns_are_zero": sum(
                1
                for p in heldout_a
                if p.features["deletion_drop"] == 0 and p.features["top_target"] == 0
            ),
            "of_pairs": len(heldout_a),
        },
    }

    payload: dict[str, Any] = {
        "status": "identity check of the model this lane scored against the frozen headline model",
        "lane": "lane-re2g",
        "pair_level_two_derivations": pair_level,
        "against_the_committed_headline": against_headline,
        "coverage_rules_reconciled": reconciliation,
        "verdict": (
            "the model scored by scripts/re2g_paired.py is the frozen headline model: no pair differs "
            "in either deletion column or in its score between two independent derivations, and every "
            "model's weighted AUPRC equals the committed headline's on both populations"
            if pair_level["all_heldout"]["scores_differing"] == 0
            and pair_level["all_heldout"]["deletion_values_differing"] == 0
            and against_headline["all_equal"]
            else "DIFFERENCES FOUND: see pair_level_two_derivations and against_the_committed_headline"
        ),
        "alphagenome_requests": 0,
    }
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "EngreitzLab/CRISPR_comparison training and held-out benchmark tables",
                "version": "v1.0.0",
            },
            {
                "accession": "data/results/crispri_published.json (this project's committed headline)",
                "version": headline.get("date", "committed 2026-09-28"),
            },
        ],
        "inputs": [
            mf.input_entry(
                crispri.KNOWLEDGE / crispri.HELDOUT,
                partition="held-out",
                url=crispri.BASE_URL + crispri.HELDOUT,
            ),
            mf.input_entry(
                crispri.KNOWLEDGE / crispri.TRAINING,
                partition="training",
                url=crispri.BASE_URL + crispri.TRAINING,
            ),
            mf.input_entry(HEADLINE, partition=None),
        ],
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {"tolerance": TOLERANCE, "model": model},
        "exclusions": ["nothing excluded: both derivations run on the published pair sets whole"],
        "partitions": {"primary_k562": 1918, "all_heldout": len(heldout_a)},
    }
    path = save_result(RESULT, payload)
    print(f"{RESULT}: {path}")
    print(json.dumps(pair_level, indent=1))
    print("all AUPRCs equal to the headline:", against_headline["all_equal"])
    print(json.dumps(against_headline["k562_1918"], indent=1))
    print("gap between the two coverage rules:", reconciliation["the_gap"]["pairs"], "pairs")
    print("VERDICT:", payload["verdict"])


if __name__ == "__main__":
    main()
