# SPDX-License-Identifier: AGPL-3.0-or-later
"""The power simulation behind clause 2's measured experiment: the parts that must be right for the
sample-size range to mean anything, and the guards that keep a reporting floor from being read as a
power result."""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


dp = _load("clause2_design_power")


# ---- the registration itself -----------------------------------------------------------------------


def test_the_model_effect_is_not_a_number_in_this_script():
    """The correction this lane exists for: 0.2725 and 0.4324 are model output. They may be NAMED in
    prose as what is excluded; neither may appear as a numeric literal the code computes with."""
    import ast

    tree = ast.parse((SCRIPTS / "clause2_design_power.py").read_text())
    numbers = {
        node.value
        for node in ast.walk(tree)
        if isinstance(node, ast.Constant) and isinstance(node.value, int | float)
    }
    assert 0.2725 not in numbers
    assert 0.4324 not in numbers


def test_every_grid_carries_a_justification():
    for name, grid in dp.PRE_REGISTRATION["grids"].items():
        assert grid.get("justification"), name
        assert len(grid["justification"]) > 80, name


def test_the_floor_is_declared_as_a_reporting_rule_and_not_power():
    text = dp.PRE_REGISTRATION["the_floor_is_not_power"]
    assert "minimum-reporting" in text
    assert dp.MIN_BLOCKS == 20
    assert min(dp.N_GRID) >= dp.MIN_BLOCKS  # no answer may be printed below the floor


# ---- the statistics --------------------------------------------------------------------------------


def test_sigma_for_icc_inverts_the_latent_scale_identity():
    for icc in (0.05, 0.21, 0.3, 0.35, 0.6):
        s = dp.sigma_for_icc(icc)
        assert dp.sigma_for_icc(0.0) == 0.0
        assert math.isclose(s**2 / (s**2 + math.pi**2 / 3), icc, rel_tol=1e-9)


def test_power_rises_with_n_and_with_the_effect():
    sd = 0.4
    a = [dp.power_at(0.1, sd, n, 0.0)["correct_direction"] for n in (20, 50, 200, 1000)]
    assert a == sorted(a)
    b = [dp.power_at(e, sd, 200, 0.0)["correct_direction"] for e in (0.02, 0.05, 0.1)]
    assert b == sorted(b)


def test_a_zero_effect_gives_the_nominal_false_positive_rate():
    got = dp.power_at(0.0, 0.4, 500, 0.0)["any_direction"]
    assert got == pytest.approx(0.05, abs=1e-3)


def test_chromosome_clustering_only_ever_costs_power():
    base = dp.power_at(0.05, 0.4, 500, 0.0)["correct_direction"]
    for rho in (0.005, 0.014, 0.03):
        assert dp.power_at(0.05, 0.4, 500, rho)["correct_direction"] <= base
    # and the design effect is 1 when every block is its own chromosome or the ICC is zero
    assert dp.power_at(0.05, 0.4, 24, 0.03)["design_effect"] == pytest.approx(1.0)


def test_the_design_effect_matches_the_registered_formula():
    n, rho = 480, 0.02
    want = 1 + (n / dp.CHROMOSOMES - 1) * rho
    assert dp.power_at(0.05, 0.4, n, rho)["design_effect"] == pytest.approx(round(want, 4))


# ---- the simulation --------------------------------------------------------------------------------


def test_the_window_arm_reproduces_the_observed_anchor():
    """The true rate is set so that the OBSERVED window rate comes back at the anchor."""
    rng = np.random.default_rng(7)
    m = dp.moments(0.1345, 1.0, 0.6674, 0.0, 1, 1, 0.0, 200_000, rng)
    assert m["window_endpoint_rate"] == pytest.approx(0.1345, abs=0.005)
    assert m["window_true_rate"] == pytest.approx(0.1345 / 0.6674, abs=1e-3)


def test_equal_rates_give_a_mean_difference_of_zero():
    rng = np.random.default_rng(11)
    m = dp.moments(0.1345, 1.0, 0.6674, 0.0, 6, 3, 0.3, 200_000, rng)
    assert abs(m["mean_difference_points"]) < 0.5


def test_a_lower_block_rate_gives_a_negative_difference_that_grows():
    rng = np.random.default_rng(13)
    got = [
        dp.moments(0.1345, r, 0.6674, 0.0, 6, 3, 0.3, 100_000, rng)["mean_difference_points"]
        for r in (1.0, 0.5, 0.0)
    ]
    assert got[0] > got[1] > got[2]
    assert got[2] < 0


def test_clustering_inside_a_block_raises_the_dispersion_of_the_difference():
    rng = np.random.default_rng(17)
    flat = dp.moments(0.1345, 0.5, 0.6674, 0.0, 6, 3, 0.0, 100_000, rng)
    clustered = dp.moments(0.1345, 0.5, 0.6674, 0.0, 6, 3, 0.35, 100_000, rng)
    assert clustered["sd_of_per_block_difference"] > flat["sd_of_per_block_difference"]


def test_more_tested_windows_shrink_the_dispersion():
    rng = np.random.default_rng(19)
    one = dp.moments(0.1345, 0.5, 0.6674, 0.0, 6, 1, 0.3, 100_000, rng)
    ten = dp.moments(0.1345, 0.5, 0.6674, 0.0, 6, 10, 0.3, 100_000, rng)
    assert ten["sd_of_per_block_difference"] < one["sd_of_per_block_difference"]


def test_a_sensitivity_of_zero_detects_nothing():
    rng = np.random.default_rng(23)
    m = dp.moments(0.1345, 0.0, 0.0, 0.0, 6, 3, 0.3, 10_000, rng)
    assert m["block_endpoint_rate"] == 0.0
    assert m["window_endpoint_rate"] == 0.0


# ---- the answer's shape ----------------------------------------------------------------------------


def test_the_answer_is_a_range_and_never_a_single_number():
    rows = [
        {"ratio": 0.5, "configuration": "a", "n_for_power": {"0.8": 100}},
        {"ratio": 0.5, "configuration": "b", "n_for_power": {"0.8": 500}},
        {"ratio": 0.5, "configuration": "c", "n_for_power": {"0.8": None}},
    ]
    got = dp.ranges(rows)["0.5"]
    assert got["smallest"] == 100 and got["largest"] == 500
    assert got["smallest_under"] == "a" and got["largest_under"] == "b"
    assert got["did_not_reach_80_percent_by_n"] == ["c"]


def test_exact_binomial_interval_brackets_the_anchor():
    lo, hi = dp._exact_binomial_interval(30, 223)
    assert lo == pytest.approx(0.0926, abs=1e-3)
    assert hi == pytest.approx(0.1864, abs=1e-3)
    assert dp._exact_binomial_interval(0, 10)[0] == 0.0
    assert dp._exact_binomial_interval(10, 10)[1] == 1.0


# ---- equivalence -----------------------------------------------------------------------------------


def test_no_equivalence_margin_is_claimed():
    got = dp.equivalence()
    assert got["admissible_margin_found"] is False
    assert got["test_run"] is None
    assert "NO DIFFERENCE DETECTED" in got["finding"]
    assert "chance band" in got["margin_search"]


def test_the_inverted_margin_is_labelled_as_not_a_justified_margin():
    got = dp.equivalence()
    assert "NOT a margin anyone has justified" in got["smallest_margin_note"]
    if got["contrasts"]:
        c = got["contrasts"][0]
        assert c["smallest_margin_a_tost_would_pass_at_points"] == pytest.approx(
            max(abs(c["ci95"][0]), abs(c["ci95"][1])), abs=0.01
        )


# ---- the denominator table -------------------------------------------------------------------------


def test_every_denominator_row_names_its_denominator_and_its_rule():
    got = dp.denominators()
    if not got["rows"]:
        pytest.skip("the committed clause 2 results are not present")
    for r in got["rows"]:
        assert r["denominator"] and r["rule"] and r["source"], r
        assert r["denominator_count"] is None or r["count"] <= r["denominator_count"], r


def test_the_two_quoted_shares_are_the_same_numerator_over_different_denominators():
    got = dp.denominators()
    if not got["rows"]:
        pytest.skip("the committed clause 2 results are not present")
    assert "59 of 882" in got["the_two_shares_that_are_quoted"]["6.7%"]
    assert "59 of the 531" in got["the_two_shares_that_are_quoted"]["11.1%"]
    assert got["the_overlap_rule"]["0.25"].startswith("72")
    assert got["the_overlap_rule"]["0.75"].startswith("10")


# ---- the additive correction in the measured arm ----------------------------------------------------


def test_the_measured_arm_keeps_its_record_and_gains_the_separation():
    """`coverage_needed` must stay exactly as it was; the new function is beside it, not instead."""
    arm = _load("clause2_measured_arm")
    assert hasattr(arm, "coverage_needed")
    assert hasattr(arm, "reporting_floor_and_power_assumptions")
    records = [
        {
            "elements": 1,
            "carrying": 4,
            "drawn": 4,
            "yes": {arm.MODEL_QUESTION: i % 2 == 0},
            "windows_yes": {arm.MODEL_QUESTION: i % 3},
        }
        for i in range(40)
    ]
    old = arm.coverage_needed(records, 1)
    new = arm.reporting_floor_and_power_assumptions(records, 1)
    assert new["record_of_what_was_computed"] == old  # the record is carried, unchanged
    assert new["reporting_floor"]["blocks_short_of_the_floor"] == arm.MIN_BLOCKS - 1
    assert "reporting" in new["reporting_floor"]["what_it_is"]
    assert new["provisional_calculation"]["is_a_power_result_for_the_measured_experiment"] is False
    for part in new["model_derived_assumptions"].values():
        if isinstance(part, dict):
            assert part["is"].startswith("AN ASSUMPTION")
    assert new["where_the_power_question_is_answered"].endswith("clause2_design_power.py")


# ---- the corrected keys on the committed results ----------------------------------------------------


def test_the_annotation_is_additive_and_refuses_otherwise(tmp_path):
    payload = {"result": "clause2_measured_arm", "primary": {"n_blocks_compared": 1}, "keep": [1, 2]}
    (tmp_path / "clause2_measured_arm.json").write_text(__import__("json").dumps(payload))
    done = dp.annotate_committed_results(tmp_path)
    assert done["clause2_measured_arm"]
    assert done["clause2_matched_control"] == "absent"
    got = __import__("json").loads((tmp_path / "clause2_measured_arm.json").read_text())
    assert {k: v for k, v in got.items() if k != dp.CORRECTION_KEY} == payload
    assert got[dp.CORRECTION_KEY]["additive"] is True


def test_every_correction_names_what_it_corrects():
    for name, entries in dp.CORRECTIONS.items():
        assert entries, name
        for field, text in entries.items():
            assert field and len(text) > 60, (name, field)


def test_the_five_review_points_are_each_corrected_somewhere():
    text = " ".join(t for e in dp.CORRECTIONS.values() for t in e.values()).lower()
    assert "target-naming frequency" in text  # (1) what 9.18% against 44.04% measures
    assert "persists after adjustment" in text  # (2) adjustment does not identify a cause
    assert "minimum-reporting rule" in text  # (3) the floor is not a power result
    assert "reciprocal overlap" in text and "denominator" in text  # (4) coverage rule and denominator
    assert "no difference detected" in text  # (5) not equivalence
