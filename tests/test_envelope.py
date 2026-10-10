# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-envelope: the additive envelope's two ends ARE two rules the engine already carries.

The whole finding of this lane rests on one claim about a definition, so it is pinned here rather
than asserted in prose: lane-strengthform's `additive_envelope` is an interval whose upper end, with
the cap it names, is `grn`'s `sum_capped` and whose lower end is `grn`'s `max`. If an edit to either
lane's code or to the engine breaks that, these tests fail and the lane's NOT_TESTABLE reading no
longer follows.

Nothing here reads a measured RNA value, and nothing fetches.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from genomeos.attribution import combine2 as c2
from genomeos.attribution import envelope as ev
from genomeos.attribution import strengthform as sf
from genomeos.runtime import grn

ROOT = Path(__file__).resolve().parent.parent


def test_the_upper_end_with_its_cap_is_the_engines_sum_capped() -> None:
    """A = min(1, W) is exactly what `grn.combine_activators` computes under `sum_capped`."""
    assert ev.FORM == "sum_capped"
    assert ev.FORM in grn.ACTIVATOR_RULES
    for terms in ([0.1, 0.2], [0.4, 0.4, 0.4], [0.9, 0.9], [1.0, 1.0, 1.0], [0.3]):
        assert grn.combine_activators(list(terms), ev.FORM) == min(1.0, sum(terms))


def test_the_lower_end_is_the_engines_max() -> None:
    assert ev.LOWER_FORM == "max"
    assert ev.LOWER_FORM in grn.ACTIVATOR_RULES
    for terms in ([0.1, 0.2], [0.4, 0.9, 0.4], [1.0, 0.0]):
        assert grn.combine_activators(list(terms), ev.LOWER_FORM) == max(terms)


def test_the_envelope_is_an_interval_with_zero_error_inside_it() -> None:
    """The comparator that lane-strengthform did not beat scores nothing between its two ends."""
    assert "additive_envelope" in sf.INTERVAL_RULES
    lower, upper = sf.predicted_boost("additive_envelope", 1.0, 0.5, 1.0)
    assert lower < upper
    assert sf.squared_error((lower + upper) / 2.0, (lower, upper)) == 0.0
    assert sf.squared_error(upper + 1.0, (lower, upper)) == pytest.approx(1.0)
    assert sf.squared_error(lower - 1.0, (lower, upper)) == pytest.approx(1.0)


def test_both_identities_hold_numerically() -> None:
    out = ev.identity_checks(draws=2000)
    assert out["holds"], out
    assert out["worst"] <= ev.IDENTITY_TOLERANCE
    for key in (
        "upper_end_is_sum_capped_in_drive_units",
        "upper_end_is_the_uncapped_additive_sum_in_boost_units",
        "lower_end_is_the_stronger_single_in_boost_units",
        "lower_end_is_grn_max",
    ):
        assert out[key] <= ev.IDENTITY_TOLERANCE, (key, out[key])


def test_the_upper_end_is_additivity_above_the_baseline() -> None:
    """P = (S1 - B) + (S2 - B) + B, i.e. boost = log2(2**b1 + 2**b2 - 1)."""
    b1, b2 = 1.5, 0.75
    _, upper = sf.predicted_boost("additive_envelope", b1, b2, max(b1, b2))
    assert upper == pytest.approx(math.log2(2.0**b1 + 2.0**b2 - 1.0))


def test_the_form_is_the_one_lane_capsumrna_ran_against_rna() -> None:
    """The identity is what licenses reading that lane's census as this form's census."""
    census = json.loads((ROOT / ev.CAPSUM_CENSUS).read_text())
    assert census["form"] == ev.FORM
    assert census["size_check"]["rna_opened"] is False


def test_the_decision_applies_7ds_minimum_and_is_not_invented_here() -> None:
    assert ev.MIN_GENES == c2.MIN_GENES == 100
    holds = {"holds": True}
    assert ev.decide(holds, 3) == "not_testable_below_7d_minimum"
    assert ev.decide(holds, 99) == "not_testable_below_7d_minimum"
    assert ev.decide(holds, 100) == "register_on_the_genes_that_differ"
    assert ev.decide({"holds": False}, 10_000) == "identity_failed"
    assert set(ev.DECISIONS) == {
        "not_testable_below_7d_minimum",
        "register_on_the_genes_that_differ",
        "identity_failed",
    }


def test_the_census_counts_come_from_the_committed_file_field_by_field() -> None:
    counts = ev.capsum_counts(ROOT / ev.CAPSUM_CENSUS)
    assert counts["form"] == ev.FORM
    assert counts["rna_opened"] is False
    assert counts["unique_fit"]["differs_from_the_mean"]["genes"] < ev.MIN_GENES
    assert counts["two_roots"]["genes"] > 0  # lane-rootpick's question, not picked here


def test_the_paired_bootstrap_is_strengthforms_own_with_the_reference_a_parameter() -> None:
    """Two clusters, hand-made errors: the delta is MSE(form) - MSE(reference) and excludes it."""
    clusters = [("a", "b"), ("a", "b"), ("c", "d")]
    errors = {
        ev.REFERENCE_FORM: [1.0, 1.0, 1.0],
        "worse": [2.0, 2.0, 2.0],
        "better": [0.5, 0.5, 0.5],
    }
    out = ev.paired_bootstrap(clusters, errors, draws=50)
    assert ev.REFERENCE_FORM not in out
    assert out["worse"]["delta_mse"] == pytest.approx(1.0)
    assert out["better"]["delta_mse"] == pytest.approx(-0.5)
    assert sf.word(out["worse"]) == "better"  # the reference is closer than `worse`
    assert sf.word(out["better"]) == "worse"


def test_the_bootstrap_must_reproduce_strengthforms_committed_intervals() -> None:
    """`check_against_committed` is what makes the pair side that lane's procedure, not a new one."""
    result = ROOT / ev.STRENGTHFORM_RESULT
    committed = json.loads(result.read_text())["intervals"]
    good = ev.check_against_committed(committed, result)
    assert good["problems"] == []
    assert good["worst"] == 0.0
    nudged = {name: dict(row) for name, row in committed.items()}
    nudged["mean"]["delta_mse"] += 1e-6
    bad = ev.check_against_committed(nudged, result)
    assert any("mean.delta_mse" in problem for problem in bad["problems"])
    missing = {name: row for name, row in committed.items() if name != "max"}
    assert any(
        "max: not recomputed" in problem
        for problem in ev.check_against_committed(missing, result)["problems"]
    )


def test_point_errors_uses_strengthforms_own_prediction_for_both_ends() -> None:
    combos = [
        sf.Combination(
            promoter="P1",
            frag1="e1",
            strand1="+",
            frag2="e2",
            strand2="+",
            baseline=1.0,
            b1=0.5,
            b2=0.5,
            activity=2.0**0.5,
        )
    ]
    out = ev.point_errors(combos, {"P1": 1.0})
    lower, upper = sf.predicted_boost("additive_envelope", 0.5, 0.5, 1.0)
    observed = combos[0].observed_boost
    assert observed == pytest.approx(0.5)
    assert out["envelope_upper_point"][0] == pytest.approx((observed - upper) ** 2)
    assert out["envelope_lower_point"][0] == pytest.approx((observed - lower) ** 2)
    assert out["inside"] == [1.0 if lower <= observed <= upper else 0.0]


def test_this_lane_declares_no_money_no_model_and_no_new_biosample() -> None:
    text = (ROOT / "genomeos" / "attribution" / "envelope.py").read_text()
    assert "No model request, no key, no money." in text
    assert "thirteen reader cells unchanged" in text
    assert "SPDX-License-Identifier: AGPL-3.0-or-later" in text.splitlines()[0]
