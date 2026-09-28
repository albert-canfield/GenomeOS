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


def test_the_annotation_is_its_own_result_and_edits_no_committed_file(tmp_path):
    payload = {"result": "clause2_measured_arm", "primary": {"n_blocks_compared": 1}, "keep": [1, 2]}
    committed = tmp_path / "clause2_measured_arm.json"
    committed.write_text(__import__("json").dumps(payload))
    before = committed.read_bytes()
    done = dp.annotate_committed_results(tmp_path)
    assert done["clause2_measured_arm"]
    assert done["clause2_matched_control"] == "absent"
    assert committed.read_bytes() == before  # item 12 S6 follow-up: no in-place edit
    got = __import__("json").loads((tmp_path / f"{dp.RESULT_CORRECTIONS}.json").read_text())
    block = got["corrections"]["clause2_measured_arm"]
    assert block == dp.correction_block(dp.CORRECTIONS["clause2_measured_arm"]) and block["additive"] is True
    assert got["result_manifest"]["complete"] and "clause2_matched_control" in got["absent"]


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


# ---- item 12 S2: the committed reference against its own anchor (lane-s2) ------------------------------


def test_the_committed_reference_produces_the_reviews_two_figures_and_not_its_anchor():
    """The first draw `sweep` makes, at the reference configuration on the committed seed, gives the
    review's 16.80% element detection and 56.96% positive windows. The anchor it was set from is 13.45%."""
    rng = np.random.default_rng(20260928)
    m = dp.moments(0.1345, 1.0, 0.6674, 0.0, 6, 3, 0.30, 200_000, rng)
    assert m["realised_element_detection_probability"] == 0.16802
    assert m["window_endpoint_rate"] == 0.5696
    assert round(m["realised_element_detection_probability"], 4) == dp.REVIEW_REPORTED["element_detection"]
    assert m["window_endpoint_rate"] - 0.1345 > 0.4


def test_the_committed_result_row_is_the_one_the_review_read():
    import json

    p = dp.RESULTS_DIR / f"{dp.RESULT}.json"
    if not p.exists():
        pytest.skip("committed result not present")
    row = next(
        r
        for r in json.loads(p.read_text())["sweep"]["rows"]
        if r["configuration"] == "reference" and r["ratio"] == 1.0
    )
    assert row["realised_element_detection_probability"] == 0.16802
    assert row["window_endpoint_rate"] == 0.5696


def test_the_closed_form_traces_both_figures_to_their_causes():
    got = dp.quadrature_of_the_committed_reference(0.1345, 0.6674, 0.30, 6)
    # the element rate: the intercept's location is at 0.1345 / 0.6674, its mean is above it
    assert got["true_rate_at_the_intercept_location"] == pytest.approx(0.2015, abs=1e-4)
    assert got["mean_true_rate"] > got["true_rate_at_the_intercept_location"] + 0.04
    assert got["observed_element_rate"] == pytest.approx(0.1680, abs=0.001)
    # the window rate: any of six elements
    assert got["positive_unit_rate_from_any_of_the_elements_alone"] == pytest.approx(0.5797, abs=1e-4)
    assert got["observed_positive_unit_rate"] == pytest.approx(0.5696, abs=0.002)
    # an observed-scale 0.30 entered as latent comes back far below 0.30 on the observed scale
    assert got["observed_scale_icc_it_implies"] < 0.2


def test_the_closed_form_agrees_with_moments_where_the_old_test_looked():
    """At 1 element and ICC 0 the three rates coincide, which is why the old anchor test passed."""
    got = dp.quadrature_of_the_committed_reference(0.1345, 0.6674, 0.0, 1)
    assert got["observed_element_rate"] == pytest.approx(0.1345, abs=1e-4)
    assert got["observed_positive_unit_rate"] == pytest.approx(0.1345, abs=1e-4)


# ---- item 12 S2: the calibrated model's machinery, on toy inputs (lane-s2) ------------------------------


def test_the_model_arm_is_excluded_from_the_calibrated_code_too():
    import ast

    tree = ast.parse((SCRIPTS / "clause2_design_power.py").read_text())
    numbers = {
        node.value
        for node in ast.walk(tree)
        if isinstance(node, ast.Constant) and isinstance(node.value, int | float)
    }
    assert not {27.25, -27.25, 0.2725, 0.4324} & numbers


def test_calibration_matches_the_anchor_at_the_window_level_and_the_icc_on_the_observed_scale():
    toy = {1: 120, 2: 60, 4: 20}
    cal = dp.calibrate(0.2, toy, 0.8, 0.02, 0.25, 0.02)
    assert cal["feasible"]
    assert cal["observed_window_rate_achieved"] == pytest.approx(0.2, abs=1e-4)
    assert cal["observed_icc_achieved"] == pytest.approx(0.25, abs=1e-3)
    assert cal["observed_chromosome_icc_achieved"] == pytest.approx(0.02, abs=1e-3)
    # observed and latent are different numbers, and the latent one is the larger
    assert cal["latent_icc_unit_plus_chromosome"] > cal["observed_icc_achieved"] + 0.05
    # windows with several tested elements are positive more often than one element is called
    assert cal["observed_element_rate_implied"] < 0.2


def test_an_unreachable_correlation_is_reported_infeasible_not_forced():
    cal = dp.calibrate(0.2, {1: 100}, 0.3, 0.0, 0.6, 0.0)
    assert cal["feasible"] is False and "largest reachable" in cal["why"]


def test_the_ratio_applies_to_the_mean_true_rate_not_the_location():
    cal = dp.calibrate(0.2, {1: 100}, 0.8, 0.0, 0.25, 0.0)
    arm = dp.block_arm(cal, 0.5)
    got = dp.marginal_true_rate(arm["mu_block"], cal["sigma_total"])
    assert got == pytest.approx(0.5 * cal["window_true_rate_latent_mean"], rel=1e-4)
    assert dp.block_arm(cal, 0.0)["mu_block"] is None


def test_unit_probability_without_spread_is_the_plain_formula():
    s, f, k, p = 0.7, 0.01, 3, 0.2
    mu = math.log(p / (1 - p))
    want = 1 - (1 - (s * p + f * (1 - p))) ** k
    assert dp.unit_positive_probability(mu, 0.0, s, f, k) == pytest.approx(want, rel=1e-9)
    assert dp.unit_positive_probability(None, 1.0, s, f, k) == pytest.approx(1 - (1 - f) ** k)


def test_the_exact_bootstrap_agrees_with_resampling():
    rng = np.random.default_rng(3)
    n, m = 40, 3
    y = rng.random((1, n)) < 0.3
    wbar = rng.integers(0, m + 1, (1, n)) / m
    t = dp.exact_block_bootstrap_tails(y, wbar, m, 0.0)
    # in integer lattice units, so the atom at 0 is not lost to floating-point sums
    d = y[0].astype(np.int64) * m - np.rint(wbar[0] * m).astype(np.int64)
    boots = d[rng.integers(0, n, (200_000, n))].sum(axis=1)
    assert t["p_below"][0] == pytest.approx((boots < 0).mean(), abs=0.005)
    assert t["p_at_or_below"][0] == pytest.approx((boots <= 0).mean(), abs=0.005)


def test_no_size_above_the_eligible_population_is_searched_or_simulated():
    assert dp.n_grid_for(45) == [20, 30, 45]
    assert dp.n_grid_for(500) == list(dp.S2_N_GRID)
    assert dp.n_grid_for(12) == []
    assert all(n <= 531 for n in dp.n_grid_for(531)) and 531 in dp.n_grid_for(531)
    cal = dp.calibrate(0.2, {1: 100}, 0.8, 0.0, 0.25, 0.0)
    with pytest.raises(ValueError):
        dp.simulate_design(
            cal,
            cal["mu_window"],
            (1, 1, 1),
            30,
            np.zeros(25, dtype=np.int64),
            10,
            np.random.default_rng(0),
            0.0,
        )


def test_cost_is_counted_in_tested_elements():
    assert dp.design_cost_per_block((1, 1, 1)) == 2
    assert dp.design_cost_per_block((6, 3, 1)) == 24
    assert dp.design_cost_per_block((2, 10, 4)) == pytest.approx(7.0)


def test_the_null_is_near_nominal_without_shared_controls_on_a_toy():
    cal = dp.calibrate(0.2, {1: 100}, 0.8, 0.0, 0.2, 0.0)
    chroms = np.repeat(np.arange(24), 20)
    got = dp.simulate_design(
        cal, cal["mu_window"], (1, 3, 1), 200, chroms, 600, np.random.default_rng(5), 0.0
    )
    assert 0.02 <= got["interval_excludes_zero"] <= 0.09
    assert got["expected_difference_points"] == 0.0


def test_a_real_difference_is_detected_more_often_than_the_null_on_a_toy():
    cal = dp.calibrate(0.2, {1: 100}, 0.8, 0.0, 0.2, 0.0)
    arm = dp.block_arm(cal, 0.25)
    delta = dp.expected_difference(cal, arm["mu_block"], 2)
    assert delta < 0
    chroms = np.repeat(np.arange(24), 20)
    got = dp.simulate_design(
        cal, arm["mu_block"], (2, 3, 4), 100, chroms, 300, np.random.default_rng(6), delta
    )
    assert got["interval_below_zero"] > 0.5
    assert got["assay_cost_elements_mean"] < 100 * 2 * (1 + 3)  # shared windows cost less than own ones


def test_the_recorder_leaves_the_committed_draw_unchanged():
    mc = _load("clause2_matched_control")
    rows = [{"start": s, "end": s + 200, "_hit": s % 3 == 0} for s in range(0, 200_000, 700)]
    targets = [
        {"start": 50_000, "end": 58_000, "length": 8_000},
        {"start": 120_000, "end": 131_000, "length": 11_000},
    ]
    preds = {"hit": lambda e: e["_hit"]}
    plain = mc.matched_windows(targets, targets, rows, preds, None, seed=11, draws=20, max_tries=400)
    rec = dp._UnitRecorder({(b["start"] + b["end"]) // 2 for b in targets})
    recorded = mc.matched_windows(
        targets,
        targets,
        rows,
        {**preds, "_record": rec.element},
        None,
        seed=11,
        draws=20,
        max_tries=400,
        raw_fn=rec.raw,
    )
    for a, b in zip(plain, recorded, strict=True):
        assert a["drawn"] == b["drawn"] and a["carrying"] == b["carrying"]
        assert a["windows_yes"]["hit"] == b["windows_yes"]["hit"]
    blocks = [u for u in rec.units if u["kind"] == "block"]
    windows = [u for u in rec.units if u["kind"] == "window"]
    assert len(blocks) == 2
    assert len(windows) == sum(r["drawn"] for r in plain)
    assert sum(1 for u in windows if u["els"]) == sum(r["carrying"] for r in plain)
    assert sum(any(e["_hit"] for e in u["els"]) for u in windows) == sum(
        r["windows_yes"]["hit"] for r in plain
    )


# ---- item 12 S2: the committed calibrated result keeps its own rules (lane-s2) ---------------------------


def _calibrated():
    import json

    p = dp.RESULTS_DIR / f"{dp.RESULT_CALIBRATED}.json"
    if not p.exists():
        pytest.skip("calibrated result not present")
    return json.loads(p.read_text())


def test_the_calibrated_result_prints_no_size_above_its_eligible_population():
    got = _calibrated()
    for row in got["sample_size_at_equal_cost"]:
        for d in row["designs"]:
            assert d["n_blocks"] == "infeasible" or d["n_blocks"] <= d["eligible_population"]
    caps = {
        int(k): v["blocks"] for k, v in got["eligible_population"]["by_elements_tested_per_block"].items()
    }
    for row in got["power_at_fixed_budgets"]:
        if row["status"] == "simulated":
            assert row["n_blocks"] <= caps[row["k"]]
        elif row["status"].startswith("infeasible"):
            assert row["n_blocks"] == "infeasible"
    for r in got["by_design"]:
        assert all(c["n_blocks"] <= r["eligible_population"] for c in r["cells"])


def test_the_calibrated_result_is_gated_and_checked_against_its_anchor():
    got = _calibrated()
    assert got["anchor_redraw_gate"]["passed"] is True
    assert got["anchor"]["windows"] == 223 and got["anchor"]["positive_windows"] == 30
    check = got["anchor_reproduced"]
    assert set(check) >= {"simulated_window_rate", "simulated_observed_icc", "reproduced"}
    cal = got["calibration_reference"]
    assert cal["latent_icc_unit_plus_chromosome"] != cal["observed_icc_achieved"]


def test_the_derived_view_only_names_designs_whose_null_is_calibrated():
    got = _calibrated()
    null = {(c["k"], c["m"], c["g"], c["n_blocks"]): c for c in got["null_calibration"]["cells"]}
    for row in dp.cheapest_with_a_calibrated_null(got):
        best = row["cheapest_with_a_calibrated_null"]
        if best is not None:
            assert null[(best["k"], best["m"], best["g"], best["n_blocks"])]["calibrated"] is True
            assert best["n_blocks"] <= best["eligible_population"]


# ---- the measured arm's reading rule, scored (lane-rule, 2026-09-29) ---------------------------------


def test_the_multi_line_tails_equal_the_single_line_function():
    rng = np.random.default_rng(11)
    y = rng.random((5, 40)) < 0.2
    wbar = (rng.random((5, 40, 3)) < 0.25).mean(axis=2)
    deltas = [0.0, -0.05, -0.1234, 0.03]
    many = dp.exact_block_bootstrap_tails_at(y, wbar, 3, deltas)
    for delta, got in zip(deltas, many, strict=True):
        one = dp.exact_block_bootstrap_tails(y, wbar, 3, delta)
        assert np.allclose(got["p_below"], one["p_below"])
        assert np.allclose(got["p_at_or_below"], one["p_at_or_below"])


def test_reading_lines_change_no_simulated_experiment():
    cal = dp.calibrate(0.2, {1: 100}, 0.8, 0.0, 0.2, 0.0)
    chroms = np.repeat(np.arange(24), 20)
    plain = dp.simulate_design(
        cal, cal["mu_window"], (1, 3, 4), 60, chroms, 300, np.random.default_rng(3), 0.0
    )
    lined = dp.simulate_design(
        cal, cal["mu_window"], (1, 3, 4), 60, chroms, 300, np.random.default_rng(3), 0.0, lines={"z": 0.0}
    )
    lined.pop("_per_experiment")
    assert lined.pop("interval_wholly_below")["z"] == plain["interval_below_zero"]
    assert lined.pop("interval_wholly_above")["z"] == plain["interval_above_zero"]
    assert lined == plain


def _toy_lines(cal, k):
    return {
        dp._line_name(r): dp.expected_difference(cal, dp.block_arm(cal, r)["mu_block"], k)
        for r in dp.S2_RATIOS
    }


def test_a_line_at_the_margin_stops_an_undecided_interval_being_read_as_a_model_failure():
    cal = dp.calibrate(0.2, {1: 100}, 0.8, 0.01, 0.2, 0.0)
    chroms = np.repeat(np.arange(24), 20)
    lines = _toy_lines(cal, 1)
    arm = dp.block_arm(cal, 0.75)
    delta = dp.expected_difference(cal, arm["mu_block"], 1)
    got = dp.simulate_design(
        cal, arm["mu_block"], (1, 3, 1), 100, chroms, 400, np.random.default_rng(8), delta, lines=lines
    )
    committed_model_failed = 1 - got["interval_below_zero"]
    revised_model_failed = got["interval_wholly_above"][dp._line_name(0.75)]
    assert committed_model_failed > 0.5
    assert revised_model_failed < 0.08


def test_reading_probabilities_add_up_and_keep_the_committed_rule_beside():
    names = {dp._line_name(r) for r in dp.S2_RATIOS}
    cell = {
        "interval_below_zero": 0.3,
        "interval_wholly_below": {n: 0.1 for n in names},
        "interval_wholly_above": {n: 0.2 for n in names},
    }
    got = dp.reading_probabilities(cell, 1)
    assert got["committed_c17eedc"] == {"model_failed": 0.7, "wording_wrong": 0.3, "cannot_decide": 0.0}
    rev = got["revised_before_the_cell_check"]
    assert rev == {"model_failed": 0.2, "wording_wrong": 0.1, "cannot_decide": 0.7}
    assert got["interval_excludes_ratio"]["0.5"] == pytest.approx(0.3)


def test_the_best_rule_bounds_its_error_on_every_wrong_ratio():
    rng = np.random.default_rng(4)
    values = {r: rng.normal(r, 0.1, 4000) for r in dp.S2_RATIOS}
    mf = dp._best_rule_on(values, 0.75, "above", 1.0)
    ww = dp._best_rule_on(values, 0.1, "below", 0.0)
    assert mf["largest_error_on_the_wrong_side"] <= 0.05 and ww["largest_error_on_the_wrong_side"] <= 0.05
    # one-sided 5% at 2.5 sd apart: Phi(2.5 - 1.645) = 0.80; at 1 sd apart: Phi(1 - 1.645) = 0.26
    assert mf["probability_at_the_true_pole"] == pytest.approx(0.80, abs=0.03)
    assert ww["probability_at_the_true_pole"] == pytest.approx(0.26, abs=0.03)


def _scored(rev_by_ratio, experiments=2000):
    cells = []
    for r in dp.S2_RATIOS:
        mf, ww = rev_by_ratio[r]
        cells.append(
            {
                "k": 1,
                "m": 3,
                "g": 1,
                "n_blocks": 100,
                "ratio": r,
                "readings": {
                    "committed_c17eedc": {"model_failed": 1.0, "wording_wrong": 0.0, "cannot_decide": 0.0},
                    "revised_before_the_cell_check": {
                        "model_failed": mf,
                        "wording_wrong": ww,
                        "cannot_decide": round(1 - mf - ww, 4),
                    },
                },
            }
        )
    best = {(1, 3, 1, 100): {"best_rule_on_the_estimator": {}, "best_rule_on_the_block_rate_alone": {}}}
    return dp.score_the_cells(cells, best, experiments)[0]


def test_a_cell_failing_a_bound_reads_cannot_decide_whatever_its_interval():
    ok = {
        1.0: (0.9, 0.0),
        0.75: (0.02, 0.0),
        0.5: (0.0, 0.0),
        0.25: (0.0, 0.01),
        0.1: (0.0, 0.03),
        0.0: (0.0, 0.9),
    }
    got = _scored(ok)
    assert got["rule_may_read_here"] and got["decides_clause_2"]
    bad = {**ok, 0.75: (0.2, 0.0)}
    got = _scored(bad)
    assert not got["rule_may_read_here"] and not got["decides_clause_2"]
    assert got["model_failed_error"] == {"largest": 0.2, "at_ratio": 0.75, "meets_the_bound": False}
    assert all(
        v == {"model_failed": 0.0, "wording_wrong": 0.0, "cannot_decide": 1.0}
        for v in got["revised_rule"].values()
    )


def test_deciding_needs_both_poles():
    one_pole = {
        1.0: (0.9, 0.0),
        0.75: (0.02, 0.0),
        0.5: (0.0, 0.0),
        0.25: (0.0, 0.0),
        0.1: (0.0, 0.02),
        0.0: (0.0, 0.3),
    }
    got = _scored(one_pole)
    assert got["rule_may_read_here"] and not got["decides_clause_2"]
    assert got["model_failed_at_ratio_1"] == 0.9 and got["wording_wrong_at_ratio_0"] == 0.3
