# SPDX-License-Identifier: AGPL-3.0-or-later
"""Step 1 of a second cell type for the CRISPRi deletion result: a blinded eligibility count.

Writes data/results/cell2_eligibility.json.

    uv run --frozen python scripts/cell2_eligibility.py

Counts, per held-out cell type, the pairs, the measured positives and the independent loci among
those positives, by the convention registered in genomeos/attribution/cell2.py and committed before
this script was first run. Blinded: it reads the benchmark's labels, coordinates, measured genes and
cell types through the project's own loader and nothing else. No prediction, no deletion value, no
model score, no AUPRC; no gain is computed, no model is fitted and no model request is made.

The decision the count serves is whether any second cell type can support a stratified test at all.
Pooled independent loci below cell2.POOLED_LOCUS_FLOOR is a no-go, recorded and stopped there. The
request counts a frozen deletion feature would need per cell type are quoted from the committed
data/results/crispri_published.json (second_cell_type.cost), not re-derived here.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import cell2, crispri  # noqa: E402
from genomeos.attribution.measured import CRISPRI_SPLIT_OF  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "cell2_eligibility"
PUBLISHED = Path("data/results/crispri_published.json")


def quoted_cost() -> dict:
    """The per-cell-type request counts already in the committed CRISPRi result, quoted as they stand."""
    payload = json.loads(PUBLISHED.read_text())
    cost = payload["second_cell_type"]["cost"]
    return {
        "source": f"{PUBLISHED} second_cell_type.cost, quoted not re-derived",
        "per_cell_type": cost,
        "budget_in_that_result": payload["second_cell_type"]["budget"],
        "why_that_result_did_not_run": payload["second_cell_type"]["why"],
        "total_requests_for_all_three": sum(v["requests_frozen_feature"] for v in cost.values()),
    }


def registration_draft(counted: dict, verdict: dict, cost: dict) -> dict:
    """The registration this lane would propose, drafted and not run. The request count goes to Albert
    for authorisation; this lane makes no request."""
    candidates = counted["candidate_second_cell_types"]
    with_feature = [c for c in candidates if counted["per_cell_type"][c]["deletion_value_available"]]
    return {
        "status": "draft only, not registered and not run; the coordinator reviews before anything else",
        "question": (
            "does the deletion feature's contribution to the CRISPRi element-to-gene benchmark hold "
            "outside K562, in a second held-out cell type measured on its own positives"
        ),
        "stratum": (
            "cell type: each candidate second cell type is its own stratum, scored and reported alone, "
            "never pooled with K562 into one number"
        ),
        "unit_of_independence": counted["convention"],
        "unit_of_independence_is_operational": counted["convention_is_operational"],
        "preset_floor": {
            "pooled_independent_loci": counted["floor"],
            "verdict_of_step_1": verdict,
            "note": "the floor was committed in genomeos/attribution/cell2.py before this count was taken",
        },
        "cell_types_that_could_carry_a_deletion_gain_at_all": with_feature,
        "cell_types_that_could_not": [c for c in candidates if c not in with_feature],
        "loci_where_a_gain_could_be_measured_at_all": counted["pooled_over_candidates_with_a_deletion_value"],
        "the_constraint_this_puts_on_the_draft": (
            "the pooled count the floor is read against is carried mostly by cell types that have no "
            "deletion value, where no deletion gain can be measured. Only the cell types listed under "
            "cell_types_that_could_carry_a_deletion_gain_at_all can be a stratum of this test as the "
            "feature stands, and their own locus count is the one that limits it"
        ),
        "cost_table_does_not_cover_every_candidate": (
            "the quoted second_cell_type.cost covers HCT116, WTC11 and Jurkat only. It carries no figure "
            "for a candidate outside those three, so no request count for such a candidate is stated "
            "here: this lane quotes the committed figures and derives none"
        ),
        "availability_rule": (
            "a cell type outside crispri.MODEL_CELLS has no deletion value, so no deletion gain is "
            "measured there: crispri.gain_where_available refuses a number rather than report the "
            "difference between two model forms as a measurement (crispri.UNAVAILABLE_GAIN)"
        ),
        "falsifier": (
            "the claim is that the deletion feature adds to the frozen baseline in a second cell type's "
            "own positives. It is falsified if, in every candidate stratum that has enough independent "
            "loci, the measured gain's 95% interval over loci includes zero or lies below it. A gain "
            "that holds only when the strata are pooled with K562, or only in a stratum whose largest "
            "locus holds most of its positives, does not count as support"
        ),
        "outcomes": {
            "interval_above_zero_in_a_stratum": (
                "the feature's contribution is reproduced outside K562 in that cell type's own "
                "positives, on the loci named in this count"
            ),
            "interval_includes_zero": (
                "no contribution is demonstrated in that cell type at this number of loci; the stratum "
                "is reported as uninformative, not as evidence of no effect"
            ),
            "interval_below_zero": (
                "the feature costs performance outside K562, which contradicts the single-cell-type "
                "reading and must be reported as such"
            ),
        },
        "model_requests_needed_per_cell_type": cost,
        "authorisation": (
            "the request count goes to Albert for authorisation; this lane made no request and holds none"
        ),
        "what_this_lane_did_not_do": [
            "computed no gain and fitted no model",
            "read no prediction, deletion value, model score or AUPRC",
            "made no model request",
        ],
    }


def manifest(heldout: list[crispri.Pair]) -> dict:
    return {
        "sources": [
            {
                "accession": "EngreitzLab/CRISPR_comparison resources/crispr_data, EPCrisprBenchmark "
                "heldout_5_cell_types (Gschwind et al.)",
                "version": "main branch, unpinned upstream; fetched 2026-09-16; pinned here by sha256",
                "url": crispri.BASE_URL,
            },
            {
                "accession": "data/results/crispri_published.json, second_cell_type.cost",
                "version": "the committed result of 2026-09-27; quoted, not re-derived",
            },
        ],
        "inputs": [
            mf.input_entry(
                crispri.KNOWLEDGE / crispri.HELDOUT,
                partition=CRISPRI_SPLIT_OF[crispri.HELDOUT],
                pairs=len(heldout),
                role="labels, coordinates, measured genes and cell types only",
            ),
            mf.input_entry(
                PUBLISHED,
                partition=None,
                role="second_cell_type.cost: the request counts quoted in the draft",
            ),
        ],
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "independent_locus_span_bp": cell2.INDEPENDENT_LOCUS_SPAN,
            "pooled_locus_floor": cell2.POOLED_LOCUS_FLOOR,
            "primary_cell": cell2.PRIMARY_CELL,
            "blinded_fields": list(cell2.BLINDED_FIELDS),
            "model_cells": list(crispri.MODEL_CELLS),
            "alphagenome_requests": 0,
            "performance_metrics_computed": 0,
            "gains_computed": 0,
        },
        "exclusions": [
            "pairs the benchmark marks ValidConnection=FALSE are dropped by crispri.parse, as in every "
            "other use of this table",
            "K562 is not a candidate second cell type: the committed result already covers it; its "
            "counts are reported for comparison only",
            "no training pair is read: the question is about the held-out cell types",
        ],
        "partitions": {
            "heldout": "the 4,378 valid held-out pairs across five cell types; the only pairs read here"
        },
    }


def main() -> int:
    heldout = crispri.load(crispri.HELDOUT)
    counted = cell2.eligibility(heldout, crispri.MODEL_CELLS)
    verdict = cell2.verdict(counted)
    cost = quoted_cost()

    payload = {
        "status": (
            "step 1 only: a blinded eligibility count of independent loci among the measured positives "
            "of each held-out cell type. No gain, no model, no request"
        ),
        "lane": "lane-cell2",
        "date": "2026-10-01",
        "question": (
            "can any second held-out cell type support a stratified test of the CRISPRi deletion "
            "result, judged by independent loci among its measured positives rather than by pair counts"
        ),
        "blinding": (
            "labels, coordinates, measured genes and cell types only, through genomeos.attribution."
            "crispri.load; no prediction, deletion value, model score or AUPRC was read or computed"
        ),
        "code_stamp_note": (
            "the stamp's dirty_code_paths may name files of other lanes working in the same checkout. "
            "The only module of theirs this script reads is genomeos/attribution/crispri.py, and the "
            "parts it reads (Pair, parse, load, MODEL_CELLS) are identical to the committed copy: the "
            "uncommitted change there is the gain-reporting fix, which this script never calls"
        ),
        "counts": counted,
        "verdict": verdict,
        "cost_quoted": cost,
        "alphagenome_requests": 0,
    }
    if verdict["meets_floor"]:
        payload["registration_draft"] = registration_draft(counted, verdict, cost)
    else:
        payload["no_go"] = {
            "decision": "no-go: the pooled independent loci are below the pre-set floor; the lane stops",
            "floor": verdict["floor"],
            "counted": verdict["pooled_independent_loci_genome_wide"],
            "note": "no registration is drafted, no gain is computed and no model request is made",
        }
    payload[mf.KEY] = manifest(heldout)
    path = save_result(RESULT, payload)

    for cell, row in counted["per_cell_type"].items():
        shapes = row["locus_shapes"]
        print(
            f"{cell:8s} pairs={row['pairs']:5d} positives={row['positives']:4d} "
            f"loci={row['independent_loci']:4d} largest={shapes['largest_locus_positives']} "
            f"({shapes['largest_locus_share_of_positives']}) deletion={row['deletion_value_available']}"
            f"{'' if row['candidate_second_cell_type'] else '  [primary, not a candidate]'}"
        )
    pooled = counted["pooled_over_candidates"]
    print(
        f"pooled over {', '.join(counted['candidate_second_cell_types'])}: "
        f"pairs={pooled['pairs']} positives={pooled['positives']} "
        f"loci_genome_wide={pooled['independent_loci_genome_wide']} "
        f"loci_summed={pooled['independent_loci_summed_per_cell_type']}"
    )
    subset = counted["pooled_over_candidates_with_a_deletion_value"]
    print(
        f"of those, in candidates with a deletion value ({', '.join(subset['cell_types']) or 'none'}): "
        f"positives={subset['positives']} loci_genome_wide={subset['independent_loci_genome_wide']} "
        "(reported beside the floor, not read against it)"
    )
    print(f"verdict: {verdict['decision']} (floor {verdict['floor']})")
    print(f"-> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
