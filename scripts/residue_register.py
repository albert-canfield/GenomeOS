# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-residue registration: the census questions, the cap ladder and the decision rule, BEFORE any count.

    uv run python scripts/residue_register.py [--no-save]

Writes `data/results/residue_questions.json`. It declares, in advance: the population (re-derived, not
copied), the eight statistics the census will report, the per-element cap levels the recovery is counted
at and why those levels, the thresholds that map the measurement onto the implication, and the prior
exposure that makes one branch of that rule the expected one. NO COUNT IS COMPUTED HERE.

Nothing is fitted, tuned or selected. No RNA value is opened, no model request is sent, no key is read.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import combine2 as c2  # noqa: E402
from genomeos.attribution import residue as rs  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "residue_questions"
SCALEFIT_CENSUS = Path("data/results/scalefit_census.json")
CENSUS_7D = Path("data/results/combine2_census.json")
MAGCAL = Path("data/results/magcal_result.json")
OWN_CODE = (
    "genomeos/attribution/residue.py",
    "scripts/residue_register.py",
    "scripts/residue_census.py",
    "tests/test_residue.py",
)


def to_be_reconciled_against() -> dict[str, Any]:
    """scalefit's committed counts, recorded as the targets of the reconciliation, NOT as inputs.

    These numbers are read here so that the record states in advance what the census's own re-derivation
    will be checked against. The census computes its counts from the compiled programs and never reads
    this field; a disagreement is a finding and not something to repair.
    """
    m = json.loads(SCALEFIT_CENSUS.read_text())["per_form"]["mean"]
    return {
        "read_from": str(SCALEFIT_CENSUS),
        "used_as": (
            "the target of a reconciliation performed AFTER this lane's own count. The census re-derives"
            " every number from the 24 compiled programs with combine2's own code and never reads these"
        ),
        "eligible_for_a_fit": m["eligible_for_a_fit"],
        "excluded_basal_not_positive": m["excluded_basal_not_positive"],
        "recovered_at_the_scale": m["recovered"],
        "still_basal_not_positive_the_residue": m["still_basal_not_positive"],
        "rule_count_quantiles_of_the_residue_as_committed": m["rule_count_quantiles_of_the_not_recovered"],
        "rule_count_quantiles_of_the_recovered_as_committed": m["rule_count_quantiles_of_the_recovered"],
    }


def record() -> dict[str, Any]:
    return {
        "lane": rs.LANE,
        "date": time.strftime("%Y-%m-%d"),
        "kind": "registration: the census questions, the cap ladder and the decision rule, before any count",
        "question": rs.QUESTION,
        "population": rs.POPULATION,
        "drops": rs.DROPS,
        "questions": rs.QUESTIONS,
        "quantiles_reported": list(rs.QUANTILES),
        "cap_rule": rs.CAP_RULE,
        "cap_levels": list(rs.CAP_LEVELS),
        "cap_primary": rs.CAP_PRIMARY,
        "cap_levels_reason": rs.CAP_LEVELS_REASON,
        "cap_bound": rs.CAP_BOUND,
        "decision_rule": rs.DECISION_RULE,
        "decision_thresholds": {
            "concentrated_median_top1_at_or_above": rs.CONCENTRATED_TOP1,
            "few_median_top3_at_or_above": rs.FEW_TOP3,
            "few_median_k_min_at_or_below": rs.FEW_KMIN,
            "diffuse_median_top1_below": rs.DIFFUSE_TOP1,
        },
        "readings": rs.READINGS,
        "to_be_reconciled_against": to_be_reconciled_against(),
        "prior_exposure": rs.PRIOR_EXPOSURE,
        "what_this_lane_does_not_do": (
            "it proposes NO combination form, registers NO RNA test, and makes NO claim about expression."
            " A gene a fit can accept is not a gene the model predicts: lane-recovered (f775e9b) had to"
            " establish that separately for the genes that recovered at the scale, and nothing here"
            " inherits its reading. Q7 counts what a per-element cap could do to the arithmetic and stops"
            " there; whether a cap belongs in the engine is another lane's registration"
        ),
        "limitations": list(rs.LIMITATIONS),
        "conventions": rs.CONVENTIONS,
        "stop_list_check": (
            "docs/ROADMAP.md §5, 'The stop list of 2026-09-27, reinstated': this lane touches no stop"
            " item. No AlphaGenome request and no deletion request, no new reader biosample, nothing in"
            " area J, no node-model work, no BioForge calibration"
        ),
        "rna_read": "none",
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
    }


def manifest() -> dict[str, Any]:
    inputs = [
        mf.input_entry(SCALEFIT_CENSUS, partition=None),
        mf.input_entry(CENSUS_7D, partition=None),
        mf.input_entry(MAGCAL, partition=None),
    ]
    return {
        "sources": [
            {"accession": "data/results/scalefit_census.json (lane-scalefit)", "version": "7f65544"},
            {"accession": "data/results/combine2_census.json (7d/7e, lane-combine2)", "version": "948cfc9"},
            {"accession": "data/results/magcal_result.json (M2h, lane-magcal)", "version": "ac5580b"},
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "scale": rs.SCALE,
            "scale_is_refitted": False,
            "cap_levels": list(rs.CAP_LEVELS),
            "cap_primary": rs.CAP_PRIMARY,
            "selection": c2.SELECTION,
            **rs.CONVENTIONS,
        },
        "exclusions": ["no gene is excluded here: this record declares the questions and computes no count"],
        "partitions": "n/a: a registration for a predictor-side census; no outcome read and no RNA opened",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args()
    out = record()
    print("questions:", len(out["questions"]))
    for k in out["questions"]:
        print("  ", k)
    print("cap levels:", out["cap_levels"], "primary", out["cap_primary"])
    print("cap bound: a cap at L cannot fit a gene with n_activating >= 1/L")
    print("thresholds:", out["decision_thresholds"])
    reconcile = {k: v for k, v in out["to_be_reconciled_against"].items() if isinstance(v, int)}
    print("to be reconciled against:", reconcile)
    if not args.no_save:
        print("saved", save_result(RESULT, out, manifest=manifest()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
