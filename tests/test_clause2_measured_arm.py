# SPDX-License-Identifier: AGPL-3.0-or-later
"""The measured arm of milestone 1.3's clause 2: the verdict rules, the denominator and the readings.

These check the parts a wrong answer would hide: that a never-measured element is never counted as a
measured negative, that lentiMPRA's per-tile rule is the one R6 fixed rather than the strongest tile,
that a significant increase counts as a measured effect (R2), that the primary's denominator is the
elements a screen actually tested, and that the three registered readings are chosen by the rules the
registration states and not by the numbers.
"""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import pytest
import tracked_paths as tp

from genomeos.attribution import measured as ms

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "clause2_measured_arm.py"
_spec = importlib.util.spec_from_file_location("clause2_measured_arm", SCRIPT)
arm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(arm)


CODING = {"AAA", "BBB"}


def pair(gene: str, start: int = 100, end: int = 200, **kw):
    fields = {
        "chrom": "chr1",
        "start": start,
        "end": end,
        "gene": gene,
        "cell": "K562",
        "dataset": "d",
        "reference": "r",
        "regulated": False,
        "significant": False,
        "effect_size": 0.0,
        "p_adjusted": 1.0,
        "power_at_effect_size_20": 0.95,
    }
    fields.update(kw)
    return ms.CrispriPair(**fields)


def layer_with(**kw) -> ms.Layer:
    return ms.Layer(chrom="chr1", **kw)


# --- the registration is a registration ------------------------------------------------------------


def test_registration_fixes_all_three_readings_before_the_run():
    assert set(arm.PRE_REGISTRATION["readings"]) == {"model_failed", "wording_wrong", "cannot_decide"}
    assert arm.PRE_REGISTRATION["registered"] == "2026-09-28"
    # the absent-measurement case is registered as a finished lane, not as a failure
    assert "FINISHED LANE" in arm.PRE_REGISTRATION["readings"]["cannot_decide"]


def test_power_floor_is_the_measured_layers_own_bar_not_a_chosen_one():
    assert arm.MIN_BLOCKS == ms.MIN_FOR_A_COMPARISON


def test_the_gate_is_the_committed_matched_primary():
    assert arm.REPRODUCE["matched_difference_points"] == -27.25
    assert arm.REPRODUCE["n_blocks_compared"] == 531


def test_the_draw_is_imported_and_not_restated():
    assert arm.SEED == arm.mc.SEED and arm.DRAWS == arm.mc.DRAWS
    assert arm.mc.matched_windows.__module__ == "clause2_matched_control"


# --- the verdict of one element ---------------------------------------------------------------------


def test_an_element_no_assay_touched_is_not_a_measured_negative():
    v = arm.verdict(layer_with(), 0, 1000, CODING)
    assert v["measured_by_any_assay"] is False
    assert v[arm.DENOMINATOR_QUESTION] is False
    assert v[arm.PRIMARY_QUESTION] is False
    assert v["measured_null_well_powered_on_every_coding_gene"] is False


def test_a_significant_decrease_on_a_coding_gene_is_the_primary():
    pairs = [pair("AAA", regulated=True, significant=True, effect_size=-0.4)]
    v = arm.verdict(layer_with(crispri=pairs), 100, 200, CODING)
    assert v[arm.PRIMARY_QUESTION] is True
    assert v["measured_regulated_decrease_only"] is True
    assert v[arm.DENOMINATOR_QUESTION] is True


def test_a_significant_increase_counts_as_a_measured_effect_and_not_as_a_null():
    pairs = [pair("AAA", significant=True, effect_size=+0.4)]
    v = arm.verdict(layer_with(crispri=pairs), 100, 200, CODING)
    assert v[arm.PRIMARY_QUESTION] is True  # R2: an increase is an effect
    assert v["measured_regulated_decrease_only"] is False  # the benchmark's Regulated is a decrease
    assert v["measured_increase_only"] is True
    assert v["measured_null_well_powered_on_every_coding_gene"] is False


def test_a_non_coding_gene_never_answers_a_clause_that_asks_for_a_coding_gene():
    pairs = [pair("LINC1", regulated=True, significant=True, effect_size=-0.4)]
    v = arm.verdict(layer_with(crispri=pairs), 100, 200, CODING)
    assert v["crispri"] is True  # the element was measured
    assert v[arm.DENOMINATOR_QUESTION] is False  # but never against a coding gene
    assert v[arm.PRIMARY_QUESTION] is False


def test_a_well_powered_null_and_an_underpowered_null_are_different_answers():
    strong = arm.verdict(layer_with(crispri=[pair("AAA", power_at_effect_size_20=0.95)]), 100, 200, CODING)
    weak = arm.verdict(layer_with(crispri=[pair("AAA", power_at_effect_size_20=0.1)]), 100, 200, CODING)
    assert strong["measured_null_well_powered_on_every_coding_gene"] is True
    assert weak["measured_null_well_powered_on_every_coding_gene"] is False
    assert weak["measured_null_underpowered_only"] is True


def test_the_overlap_rule_decides_what_counts_as_measured():
    pairs = [pair("AAA", start=100, end=120, regulated=True, significant=True, effect_size=-0.4)]
    layer = layer_with(crispri=pairs)
    assert arm.verdict(layer, 100, 1000, CODING)["crispri"] is False  # 20 bp of a 900 bp element
    assert arm.verdict(layer, 100, 120, CODING)["crispri"] is True


def test_lentimpra_follows_the_per_tile_rule_and_never_the_strongest_tile():
    from genomeos.attribution import mpra

    tiles = [
        mpra.Element(chrom="chr1", start=100, end=300, name="a", activity={"K562": 5.0}),
        mpra.Element(chrom="chr1", start=100, end=300, name="b", activity={"K562": -1.0}),
        mpra.Element(chrom="chr1", start=100, end=300, name="c", activity={"K562": -1.0}),
    ]
    v = arm.verdict(layer_with(lentimpra=tiles), 100, 300, CODING)
    assert v["lentimpra"] is True
    assert v["lentimpra_active"] is False  # two silent tiles out of three, though one is strong
    assert v["measured_active_but_no_gene_named"] is False


def test_reporter_activity_names_no_gene_so_it_never_enters_the_primary():
    from genomeos.attribution import mpra

    tiles = [mpra.Element(chrom="chr1", start=100, end=300, name="a", activity={"K562": 5.0})]
    v = arm.verdict(layer_with(lentimpra=tiles), 100, 300, CODING)
    assert v["lentimpra_active"] is True
    assert v["measured_active_but_no_gene_named"] is True
    assert v[arm.PRIMARY_QUESTION] is False
    assert v[arm.DENOMINATOR_QUESTION] is False


def test_a_vista_positive_is_activity_and_a_vista_negative_is_still_a_measurement():
    from genomeos.attribution import vista

    pos = vista.VistaElement(
        id="hs1", chrom="chr1", start=100, end=300, status="positive", tissues=("heart",), experiments=1
    )
    neg = vista.VistaElement(
        id="hs2", chrom="chr1", start=100, end=300, status="negative", tissues=(), experiments=1
    )
    vp = arm.verdict(layer_with(vista=[pos]), 100, 300, CODING)
    vn = arm.verdict(layer_with(vista=[neg]), 100, 300, CODING)
    assert vp["vista_positive"] is True and vp["measured_active_but_no_gene_named"] is True
    assert vn["measured_by_any_assay"] is True and vn["vista_negative"] is True
    assert vn["measured_active_but_no_gene_named"] is False
    assert vp[arm.PRIMARY_QUESTION] is False  # VISTA names no gene either


# --- the denominator and the statistic --------------------------------------------------------------


def record(yes_denom: bool, yes_q: bool, win_denom: int, win_q: int) -> dict:
    return {
        "yes": {arm.DENOMINATOR_QUESTION: yes_denom, arm.PRIMARY_QUESTION: yes_q},
        "windows_yes": {arm.DENOMINATOR_QUESTION: win_denom, arm.PRIMARY_QUESTION: win_q},
    }


def test_a_block_no_screen_reached_leaves_the_primary_rather_than_scoring_zero():
    recs = [record(False, False, 10, 5), record(True, True, 10, 5)]
    assert len(arm.restricted(recs)) == 1


def test_a_block_whose_windows_no_screen_reached_also_leaves_the_primary():
    assert arm.restricted([record(True, True, 0, 0)]) == []


def test_below_the_power_floor_no_interval_is_read():
    by_chrom = {"chr1": [record(True, True, 10, 5) for _ in range(arm.MIN_BLOCKS - 1)]}
    out = arm.summarise_measured(by_chrom, arm.PRIMARY_QUESTION)
    assert out["below_the_power_floor"] is True
    assert out["matched_difference_points"] is None
    assert str(arm.MIN_BLOCKS) in out["why"]


def test_at_the_floor_the_interval_is_read_and_the_window_rate_is_over_reached_windows():
    by_chrom = {"chr1": [record(True, True, 10, 5) for _ in range(arm.MIN_BLOCKS)]}
    out = arm.summarise_measured(by_chrom, arm.PRIMARY_QUESTION)
    assert out["n_blocks_compared"] == arm.MIN_BLOCKS
    assert out["matched_difference_points"] == pytest.approx(50.0)  # 1 - 5/10
    assert out["pooled"]["window_rate"] == pytest.approx(0.5)


def test_the_readings_are_chosen_by_the_registered_rule_and_not_by_the_number():
    assert arm.measured_reading({"matched_difference_points": None})["outcome"] == "cannot_decide"
    below = {"matched_difference_points": -30.0, "ci95_over_blocks": [-40.0, -20.0]}
    assert arm.measured_reading(below)["outcome"] == "wording_wrong"
    across = {"matched_difference_points": -5.0, "ci95_over_blocks": [-20.0, +10.0]}
    assert arm.measured_reading(across)["outcome"] == "model_failed"
    above = {"matched_difference_points": +10.0, "ci95_over_blocks": [+2.0, +18.0]}
    assert arm.measured_reading(above)["outcome"] == "model_failed"
    assert arm.measured_reading({"matched_difference_points": None})["decides_clause_2"] is False


def test_coverage_needed_uses_the_registered_formula():
    recs = [
        {
            "elements": 1,
            "carrying": 1,
            "yes": {arm.MODEL_QUESTION: i % 2 == 0},
            "windows_yes": {arm.MODEL_QUESTION: 0},
        }
        for i in range(40)
    ]
    got = arm.coverage_needed(recs, compared_now=3)
    s = got["model_arm_sd_of_per_block_differences"]
    assert got["n_for_80_percent_power"] == math.ceil(arm.POWER_CONSTANT * s**2 / arm.EFFECT_TO_DETECT**2)
    assert got["blocks_short_of_the_floor"] == arm.MIN_BLOCKS - 3


# --- coverage --------------------------------------------------------------------------------------


def blank() -> dict:
    v = dict.fromkeys(tuple(arm.QUESTIONS[1:]) + arm.COVERAGE_FLAGS, False)
    v["crispri_coding_genes_tested"] = 0
    return v


def test_coverage_counts_blocks_and_elements_apart_and_names_the_untouched_blocks():
    els = [{"start": 0, "end": 10, "_measured": blank()}, {"start": 20, "end": 30, "_measured": blank()}]
    els[0]["_measured"]["measured_by_any_assay"] = True
    els[0]["_measured"]["crispri"] = True
    els[0]["_measured"]["crispri_coding_genes_tested"] = 2
    rows = [{"block": "chr1:0-100", "elements": 2}, {"block": "chr1:100-200", "elements": 0}]
    got = arm.coverage(rows, {"chr1:0-100": els, "chr1:100-200": []})
    assert got["blocks"] == 2
    assert got["blocks_carrying_a_scored_element"] == 1
    assert got["scored_elements"] == 2
    assert got["blocks_with"]["measured_by_any_assay"] == 1
    assert got["elements_with"]["measured_by_any_assay"] == 1
    assert got["blocks_with_no_measured_element"] == 1  # the empty block counts as untouched
    assert got["crispri_coding_pairs_tested"] == 2


def test_pooled_coverage_adds_the_chromosomes_up():
    one = arm.coverage([{"block": "b", "elements": 0}], {"b": []})
    pooled = arm.pool_coverage([one, one])
    assert pooled["blocks"] == 2
    assert pooled["blocks_with"]["measured_by_any_assay"] == 0


# --- the reading rule revisited (lane-rule, registered 2026-09-29) ------------------------------------


def test_the_committed_rule_stays_the_record_beside_the_revised_one():
    assert set(arm.PRE_REGISTRATION["readings"]) == {"model_failed", "wording_wrong", "cannot_decide"}
    assert arm.PRE_REGISTRATION["registered"] == "2026-09-28"
    reg = arm.READING_RULE_2026_09_29
    assert reg["registered"] == "2026-09-29" and reg["lane"] == "lane-rule"
    assert set(reg["readings"]) == set(arm.PRE_REGISTRATION["readings"])
    across = {"matched_difference_points": -1.0, "ci95_over_blocks": [-3.5, 1.5]}
    got = arm.revised_reading(across, 1)
    assert got["committed_rule_c17eedc"] == arm.measured_reading(across)["outcome"] == "model_failed"


def test_the_error_bounds_are_registered_before_the_rule_is_scored():
    reg = arm.READING_RULE_2026_09_29["error_rates_a_rule_must_meet"]
    assert "0.75 or below" in reg["model_failed"] and "0.1 or above" in reg["wording_wrong"]
    assert arm.REVISED_RULE_ERROR_BOUND == 0.05 and arm.REVISED_RULE_DECIDES_AT == 0.8
    assert arm.REVISED_RULE_MARGINS == {"model_failed": 0.75, "wording_wrong": 0.1}
    assert "cannot_decide_is_never_an_error" in reg


def test_an_interval_that_merely_reaches_zero_is_not_a_model_failure():
    across = {"matched_difference_points": -1.0, "ci95_over_blocks": [-3.5, 1.5]}
    assert arm.revised_reading(across, 1)["outcome"] == "cannot_decide"
    assert arm.revised_reading(across, 1)["decides_clause_2"] is False


def test_model_failed_needs_the_interval_above_the_three_quarters_line():
    line = arm.REVISED_RULE_LINES_POINTS[1][0.75]
    above = {"matched_difference_points": 0.0, "ci95_over_blocks": [line + 0.01, 2.5]}
    at = {"matched_difference_points": 0.0, "ci95_over_blocks": [line, 2.5]}
    assert arm.revised_reading(above, 1)["outcome"] == "model_failed"
    assert arm.revised_reading(at, 1)["outcome"] == "cannot_decide"


def test_wording_wrong_needs_the_interval_below_the_one_tenth_line_not_below_zero():
    line = arm.REVISED_RULE_LINES_POINTS[1][0.1]
    below_zero_only = {"matched_difference_points": -5.0, "ci95_over_blocks": [-8.0, -2.0]}
    assert arm.measured_reading(below_zero_only)["outcome"] == "wording_wrong"
    assert arm.revised_reading(below_zero_only, 1)["outcome"] == "cannot_decide"
    below = {"matched_difference_points": -12.0, "ci95_over_blocks": [-14.0, line - 0.01]}
    assert arm.revised_reading(below, 1)["outcome"] == "wording_wrong"


def test_no_line_no_floor_or_a_failed_cell_reads_cannot_decide():
    clear = {"matched_difference_points": 0.0, "ci95_over_blocks": [-1.0, 1.0]}
    assert arm.revised_reading(clear, None)["outcome"] == "cannot_decide"  # varying elements per unit
    assert arm.revised_reading(clear, 4)["outcome"] == "cannot_decide"  # no registered design
    assert arm.revised_reading(clear, 1, admissible=False)["outcome"] == "cannot_decide"
    assert arm.revised_reading({"matched_difference_points": None}, 1)["outcome"] == "cannot_decide"


def test_the_ratios_an_interval_excludes_are_reported_and_are_not_a_reading():
    got = arm.revised_reading({"matched_difference_points": -4.5, "ci95_over_blocks": [-6.0, -3.0]}, 1)
    assert got["outcome"] == "cannot_decide"
    assert got["ratios_the_interval_excludes"] == [1.0, 0.75, 0.25, 0.1, 0.0]


def test_the_lines_are_the_committed_calibrations_expected_differences():
    import json

    p = Path(__file__).resolve().parents[1] / "data" / "results" / "clause2_design_power_calibrated.json"
    # the calibrated result is TRACKED, so "not in this checkout" was never a reachable state
    tp.must_be_present(p, was="the calibrated result is not in this checkout")
    seen: dict[int, dict[float, float]] = {}
    for r in json.loads(p.read_text())["by_design"]:
        seen.setdefault(r["design"]["k"], {})[r["ratio"]] = r["expected_difference_points"]
    assert seen == arm.REVISED_RULE_LINES_POINTS
