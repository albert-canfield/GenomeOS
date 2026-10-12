# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-scalefit registration: the census rule and the floor, committed BEFORE the answer is computed.

    uv run python scripts/scalefit_register.py [--no-save]

Writes `data/results/scalefit_rule.json`: 7d's `basal_not_positive` population read first-hand from its
committed census, the scaling by M2h's fixed constant, what counts as recovered, the floor below which
the answer is "the magnitude error does not explain the exclusions", and the prior exposure that makes
that floor not blind. NO COUNT IS COMPUTED HERE. No RNA, no model request, nothing fetched.
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
from genomeos.attribution import scalefit as sf  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "scalefit_rule"
CENSUS_7D = Path("data/results/combine2_census.json")
MAGCAL = Path("data/results/magcal_result.json")
OWN_CODE = (
    "genomeos/attribution/scalefit.py",
    "scripts/scalefit_register.py",
    "scripts/scalefit_census.py",
    "tests/test_scalefit.py",
)


def population() -> dict[str, Any]:
    """7d's excluded counts, read from its committed census (not from how it is described anywhere)."""
    d = json.loads(CENSUS_7D.read_text())
    rs = d["reader_selection"]
    out: dict[str, Any] = {
        "read_from": str(CENSUS_7D),
        "read_from_field": "reader_selection.forms.<form>.not_runnable.basal_not_positive",
        "lane": d["lane"],
        "date": d["date"],
        "eligible_for_a_fit": rs["eligible_for_a_fit"],
        "two_or_more_rules": rs["two_or_more_rules"],
        "per_form": {
            f: {
                "basal_not_positive": rs["forms"][f]["not_runnable"]["basal_not_positive"],
                "inhibitor_strength_above_1": rs["forms"][f]["not_runnable"]["inhibitor_strength_above_1"],
                "runnable": rs["forms"][f]["runnable"],
            }
            for f in sf.CONVENTIONS["forms"]
        },
        "summed_single_removal_drop_quantiles_as_committed_by_7d": rs["summed_single_removal_drop_quantiles"],
        "exclusion_code": (
            "genomeos/attribution/combine2.py:fit — `if not q < 1.0: return Fit(form, False,"
            " 'basal_not_positive', q=q, z=z)`. Read from the code, not from a description of it"
        ),
        "roughly_13000_is": (
            "the mean form's 12,500 of 17,983; ROADMAP §5 item 7d writes 'the roughly 13,000 genes no"
            " form can fit' and 7h's census gives 11,802 informative under neither 7d form. This lane's"
            " primary population is the mean form's exact 12,500, with competition's 14,935 beside it"
        ),
    }
    return out


def record() -> dict[str, Any]:
    pop = population()
    return {
        "lane": sf.LANE,
        "date": time.strftime("%Y-%m-%d"),
        "kind": "registration: the census rule and the floor, before any count",
        "question": (
            "7d/7e excluded 12,500 of 17,983 multi-element genes as basal_not_positive: their summed"
            " single-removal drops reach or exceed the gene's whole intact level, so no basal rate is"
            " left and the gene cannot be fitted at all. M2h measured that the compiled rules run about"
            " five times too large (one fixed scale, +0.2146, held out by sha256 of the symbol). If the"
            " removals are too large by that factor, is their sum exceeding the gene's level arithmetic"
            " rather than biology? How many of the excluded genes have a positive basal rate again when"
            " every drop is scaled by M2h's constant, taken fixed and NOT refitted?"
        ),
        "population": pop,
        "scale": sf.SCALE,
        "scale_rounded": sf.SCALE_ROUNDED,
        "scale_ci": list(sf.SCALE_CI),
        "scale_source": sf.SCALE_SOURCE,
        "nothing_is_fitted": sf.RULE["nothing_is_fitted"],
        "recovery_threshold_on_the_unscaled_ratio": sf.RECOVERY_THRESHOLD,
        "rule": sf.RULE,
        "floor": sf.FLOOR,
        "floor_share_of_the_excluded": sf.FLOOR_SHARE,
        "floor_reason": sf.FLOOR_REASON,
        "strong_bar": sf.STRONG_BAR,
        "strong_bar_reason": sf.STRONG_BAR_REASON,
        "prior_exposure": sf.PRIOR_EXPOSURE,
        "verdict_wordings": {
            "below_the_floor": sf.verdict(sf.FLOOR - 1, pop["per_form"]["mean"]["basal_not_positive"]),
            "floor_cleared": sf.verdict(sf.FLOOR, pop["per_form"]["mean"]["basal_not_positive"]),
            "strong_bar_cleared": sf.verdict(sf.STRONG_BAR, pop["per_form"]["mean"]["basal_not_positive"]),
        },
        "what_this_lane_does_not_claim": (
            "nothing about expression. A recovered gene is a gene the fit can now accept, nothing more."
            " If the count clears the floor, an RNA registration follows as a SEPARATE lane and not this"
            " one; this lane opens no RNA value at any point"
        ),
        "conventions": sf.CONVENTIONS,
        "limitations": list(sf.LIMITATIONS),
        "rna_read": "none",
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
    }


def manifest() -> dict[str, Any]:
    inputs = [mf.input_entry(CENSUS_7D, partition=None), mf.input_entry(MAGCAL, partition=None)]
    return {
        "sources": [
            {"accession": "data/results/combine2_census.json (7d/7e, lane-combine2)", "version": "948cfc9"},
            {"accession": "data/results/magcal_result.json (M2h, lane-magcal)", "version": "ac5580b"},
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "scale": sf.SCALE,
            "scale_is_refitted": False,
            "floor": sf.FLOOR,
            "strong_bar": sf.STRONG_BAR,
            **sf.CONVENTIONS,
        },
        "exclusions": [
            "no gene is excluded here: this record declares the rule and computes no count",
        ],
        "partitions": "n/a: a registration; no outcome read",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args()
    out = record()
    print("excluded, mean form:", out["population"]["per_form"]["mean"]["basal_not_positive"])
    print("excluded, competition:", out["population"]["per_form"]["competition"]["basal_not_positive"])
    print("floor:", out["floor"], "strong bar:", out["strong_bar"])
    print("threshold on the unscaled ratio:", round(out["recovery_threshold_on_the_unscaled_ratio"], 4))
    if not args.no_save:
        print("saved", save_result(RESULT, out, manifest=manifest()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
