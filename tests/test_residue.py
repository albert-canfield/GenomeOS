# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-residue: the concentration statistics, the cap's arithmetic bound, and the guards on the record.

Every test here is predictor-side. None opens an RNA value, none fits anything, and none asserts a count
of real genes except through the committed record itself.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest
import tracked_paths as tp

from genomeos.attribution import combine2
from genomeos.attribution import residue as rs
from genomeos.attribution import scalefit as sf

QUESTIONS_FILE = Path("data/results/residue_questions.json")
CENSUS_FILE = Path("data/results/residue_census.json")
GENES_FILE = Path("data/results/residue_genes.json")

# ==== the scale is inherited, not re-declared, and nothing here fits =====================================


def test_the_scale_is_scalefit_s_constant_and_is_not_re_declared_in_this_module() -> None:
    assert rs.SCALE == sf.SCALE
    assert rs.SCALE_ROUNDED == sf.SCALE_ROUNDED
    text = Path("genomeos/attribution/residue.py").read_text()
    assert "0.2145911381525326" not in text, "the scale must come from scalefit, not a second copy"


def test_no_module_or_script_of_this_lane_fits_optimises_or_re_estimates_anything() -> None:
    for path in ("genomeos/attribution/residue.py", "scripts/residue_register.py"):
        text = Path(path).read_text()
        for banned in ("brentq", "minimize", "curve_fit", "lstsq", "polyfit", "least_squares"):
            assert banned not in text, f"{banned} in {path} would make a declared constant a fit"


def test_the_residue_threshold_is_the_reciprocal_of_the_scale() -> None:
    assert pytest.approx(4.660025, abs=1e-5) == rs.RESIDUE_THRESHOLD
    assert pytest.approx(1.0, abs=1e-12) == rs.RESIDUE_THRESHOLD * rs.SCALE


# ==== the per-gene statistics, on planted vectors ========================================================


def test_the_drops_are_the_activating_ones_only_and_come_back_largest_first() -> None:
    # 0.4 and 0.1 are activating drops; 2.0 is an inhibiting lift and enters nothing
    assert rs.activating_drops({"a": 0.6, "b": 0.9, "c": 2.0}) == pytest.approx([0.4, 0.1])


def test_top_share_is_the_largest_over_the_sum() -> None:
    assert rs.top_share([0.8, 0.1, 0.1], 1) == pytest.approx(0.8)
    assert rs.top_share([0.8, 0.1, 0.1], 3) == pytest.approx(1.0)
    assert rs.top_share([0.25] * 4, 1) == pytest.approx(0.25)
    assert rs.top_share([0.25] * 4, 3) == pytest.approx(0.75)


def test_top3_is_forced_to_1_when_the_gene_has_three_or_fewer_activating_elements() -> None:
    assert rs.top_share([0.9, 0.05], 3) == pytest.approx(1.0)
    assert rs.top_share([0.5, 0.3, 0.2], 3) == pytest.approx(1.0)


def test_n_effective_is_1_for_one_element_and_n_for_n_equal_ones() -> None:
    assert rs.n_effective([0.7]) == pytest.approx(1.0)
    for n in (2, 5, 58):
        assert rs.n_effective([0.3] * n) == pytest.approx(float(n))
    # and strictly between for anything in between
    assert 1.0 < rs.n_effective([0.9, 0.1, 0.1]) < 3.0


def test_k_min_is_0_when_the_scaled_sum_already_fits_and_1_when_one_drop_carries_the_excess() -> None:
    # a sum already under 1 at the scale: no element need go
    assert rs.k_min_to_fit([1.0, 1.0]) == 0.0
    # 31 at the scale is 6.65; dropping 25 leaves 6 -> 1.29, still over; dropping 25 and 5 leaves 0.21
    assert rs.k_min_to_fit([25.0, 5.0, 1.0]) == 2.0
    # one dominant drop whose removal is enough
    assert rs.k_min_to_fit([20.0, 1.0, 1.0]) == 1.0
    # k must be BELOW n: a two-element gene that needs both gone is NaN, not 2
    assert math.isnan(rs.k_min_to_fit([25.0, 5.0]))


def test_k_min_is_nan_when_no_k_below_n_brings_the_sum_under_the_level() -> None:
    # every drop is 5, scaled 1.073 > 1: even one element alone exceeds the gene, so no k < n works
    assert math.isnan(rs.k_min_to_fit([5.0] * 4))


def test_a_cap_cannot_fit_a_gene_with_at_least_1_over_L_activating_elements_whatever_its_drops() -> None:
    """CAP_BOUND, the arithmetic stated before the count: capped sum <= n * L."""
    for cap in rs.CAP_LEVELS:
        n_at_bound = math.ceil(1.0 / cap)
        huge = [1e6] * n_at_bound  # every drop is far above the cap, so each contributes exactly L
        assert not rs.fits_under_cap(huge, cap), f"cap {cap} must not fit {n_at_bound} saturated elements"
        if n_at_bound > 1:
            assert rs.fits_under_cap([1e6] * (n_at_bound - 1), cap), f"cap {cap}: n < 1/L must always fit"


def test_the_capped_sum_is_the_scaled_sum_when_no_drop_reaches_the_cap() -> None:
    drops = [0.4, 0.3, 0.2]
    assert rs.capped_sum(drops, 1.0) == pytest.approx(rs.SCALE * sum(drops))


def test_the_cap_agrees_with_combine2_s_own_fit_on_the_capped_observations() -> None:
    """The cross-check CAP_RULE declares: the capped sum's test is 7d's `q < 1`, computed by 7d's code."""
    rho = {f"e{i}": 1.0 - d for i, d in enumerate([3.0, 2.0, 1.5, 0.9, 0.4, 0.2, 0.1])}
    rho["inh"] = 1.7  # an inhibiting element, which the cap leaves alone
    drops = rs.activating_drops(rho)
    for cap in rs.CAP_LEVELS:
        f = combine2.fit("mean", rs.capped_rho(rho, cap))
        assert f.q == pytest.approx(rs.capped_sum(drops, cap))
        assert rs.fits_under_cap(drops, cap) == (f.q < 1.0)


def test_the_cap_leaves_an_inhibiting_lift_at_the_scale_and_never_caps_it() -> None:
    out = rs.capped_rho({"inh": 11.0}, 0.01)
    assert out["inh"] - 1.0 == pytest.approx(rs.SCALE * 10.0)  # scaled, not capped at 0.01


# ==== the statistic that decides is NOT one scalefit already published ===================================


def test_PLANTED_two_genes_with_the_same_element_count_and_the_same_summed_drop_differ_in_top1() -> None:
    """The guard on this whole lane: scalefit published the residue's RULE COUNTS and its summed drops.
    Neither determines the answer to this lane's question, so a census that only re-derived those would
    preserve the quantity it watches and say nothing. These two genes are indistinguishable in both of
    scalefit's statistics and opposite in every statistic of this one."""
    concentrated = [19.0] + [0.1] * 9  # sum 19.9, ten elements
    diffuse = [1.99] * 10  # sum 19.9, ten elements
    assert sum(concentrated) == pytest.approx(sum(diffuse))
    assert len(concentrated) == len(diffuse)
    assert rs.SCALE * sum(concentrated) == pytest.approx(rs.SCALE * sum(diffuse))
    assert rs.top_share(concentrated, 1) > 0.9 > 0.2 > rs.top_share(diffuse, 1)
    assert rs.k_min_to_fit(concentrated) == 1.0
    assert rs.k_min_to_fit(diffuse) > 1.0
    assert rs.n_effective(concentrated) < 2.0 < rs.n_effective(diffuse)


# ==== the decision rule, every branch reachable before the count =========================================


def test_each_branch_of_the_decision_rule_is_reachable_and_the_thresholds_are_the_declared_ones() -> None:
    assert rs.decide(0.6, 0.9, 1.0) == "concentrated"  # median top1 >= 0.5
    assert rs.decide(0.5, 0.5, 9.0) == "concentrated"  # exactly at the threshold is concentrated
    assert rs.decide(0.3, 0.8, 9.0) == "few"  # top3 >= 0.75
    assert rs.decide(0.3, 0.5, 3.0) == "few"  # k_min <= 3
    assert rs.decide(0.1, 0.3, 9.0) == "diffuse"  # top1 < 0.25 and k_min > 3
    assert rs.decide(0.3, 0.5, 9.0) == "mixed"  # between every threshold
    assert rs.decide(0.1, 0.3, math.nan) == "diffuse"  # a residue no k < n can fit is not less diffuse


def test_every_branch_has_a_registered_reading_and_only_the_diffuse_one_rejects_a_cap() -> None:
    assert set(rs.READINGS) == {"concentrated", "few", "diffuse", "mixed"}
    assert "PER-ELEMENT CAP is the indicated next" in rs.READINGS["concentrated"]
    assert "NOT the indicated next form" in rs.READINGS["diffuse"]
    for branch, text in rs.READINGS.items():
        assert rs.reading(branch).startswith(text)
        assert "predict" not in text and "RNA" not in text, f"{branch} must make no expression claim"


def test_the_median_the_decision_reads_is_the_0_point_5_quantile_the_record_reports() -> None:
    xs = [float(i) for i in range(100)]
    assert rs.median(xs) == rs.distribution(xs)["quantiles"]["0.5"]


def test_the_distribution_reports_the_whole_shape_and_counts_the_nans_it_dropped() -> None:
    d = rs.distribution([1.0, 2.0, math.nan, 3.0])
    assert d["n"] == 3 and d["n_nan"] == 1
    assert d["min"] == 1.0 and d["max"] == 3.0 and d["mean"] == pytest.approx(2.0)
    assert set(d["quantiles"]) == {str(q) for q in rs.QUANTILES}


# ==== the per-gene row ===================================================================================


def test_a_gene_row_carries_every_declared_statistic_and_its_cap_verdict_at_each_level() -> None:
    rho = {"a": 0.2, "b": 0.5, "c": 1.6}
    row = rs.gene_residue("SYM", "chr21", rho, in_residue=False).as_row()
    assert row["n_rules"] == 3 and row["n_activating"] == 2 and row["n_inhibiting"] == 1
    assert row["sum_d"] == pytest.approx(1.3)
    assert row["q_after"] == pytest.approx(rs.SCALE * 1.3)
    assert row["top1"] == pytest.approx(0.8 / 1.3)
    assert set(row["fits_under"]) == {str(c) for c in rs.CAP_LEVELS}


# ==== the ladder itself =================================================================================


def test_the_cap_ladder_is_round_levels_only_so_no_level_can_have_been_tuned() -> None:
    for level in rs.CAP_LEVELS:
        assert level in (1.0, 0.5, 0.25, 0.1, 0.05, 0.02, 0.01), f"{level} is not a round ladder rung"
    assert len(set(rs.CAP_LEVELS)) == len(rs.CAP_LEVELS)
    assert list(rs.CAP_LEVELS) == sorted(rs.CAP_LEVELS, reverse=True)


def test_the_module_discloses_the_prior_exposure_and_claims_nothing_about_expression() -> None:
    assert "58" in rs.PRIOR_EXPOSURE, "the residue's published median rule count must be disclosed"
    assert "PROPOSES NO FORM" in rs.__doc__
    assert "NO RNA VALUE IS OPENED" in rs.__doc__
    assert "PROPOSES NO FORM" in rs.CAP_BOUND or "PROPOSES NO FORM" in rs.CAP_BOUND.upper()


# ==== the committed registration ========================================================================


def test_the_registration_is_committed_and_declares_the_cap_ladder_the_code_uses() -> None:
    tp.must_be_present(
        QUESTIONS_FILE,
        was="the questions, the cap ladder and the decision rule, committed before any count",
    )
    d = json.loads(QUESTIONS_FILE.read_text())
    assert d["lane"] == rs.LANE
    assert d["cap_levels"] == list(rs.CAP_LEVELS)
    assert d["cap_primary"] == rs.CAP_PRIMARY
    assert d["cap_primary"] in d["cap_levels"]
    assert d["quantiles_reported"] == list(rs.QUANTILES)
    assert d["questions"] == rs.QUESTIONS
    assert len(d["questions"]) == 8
    assert d["decision_thresholds"] == {
        "concentrated_median_top1_at_or_above": rs.CONCENTRATED_TOP1,
        "few_median_top3_at_or_above": rs.FEW_TOP3,
        "few_median_k_min_at_or_below": rs.FEW_KMIN,
        "diffuse_median_top1_below": rs.DIFFUSE_TOP1,
    }
    assert d["readings"] == rs.READINGS


def test_the_registration_holds_no_count_of_its_own_and_names_scalefit_s_only_as_a_target() -> None:
    """The record that declares the questions must not already contain an answer to one of them."""
    d = json.loads(QUESTIONS_FILE.read_text())
    assert set(d) & {"Q1_counts", "residue", "recovered", "branch", "verdict"} == set()
    target = d["to_be_reconciled_against"]
    assert "AFTER this lane's own count" in target["used_as"]
    assert target["still_basal_not_positive_the_residue"] > 0  # recorded as the target of a check
    assert d["rna_read"] == "none" and d["alphagenome_requests"] == 0
    assert "NO claim about expression" in d["what_this_lane_does_not_do"]
    assert "proposes NO combination form" in d["what_this_lane_does_not_do"]


def test_the_registration_discloses_its_prior_exposure_word_for_word_from_the_module() -> None:
    d = json.loads(QUESTIONS_FILE.read_text())
    assert d["prior_exposure"] == rs.PRIOR_EXPOSURE
    assert "58" in d["prior_exposure"]
    assert d["cap_bound"] == rs.CAP_BOUND
    assert d["decision_rule"] == rs.DECISION_RULE


# ==== the committed census: the outcome is RECOMPUTED, never read as written ==============================


def _census() -> dict:
    tp.must_be_present(
        CENSUS_FILE,
        GENES_FILE,
        was="the one census run of this lane and its per-gene table, written at the committed questions",
    )
    return json.loads(CENSUS_FILE.read_text())


def test_the_branch_and_the_verdict_recompute_from_the_reported_medians_by_the_decision_function() -> None:
    """The lesson of 2026-10-12: a committed reading must be recomputable from the numbers beside it."""
    d = _census()
    inp = d["decision_inputs"]
    branch = rs.decide(
        inp["median_top1_of_the_residue"],
        inp["median_top3_of_the_residue"],
        inp["median_k_min_of_the_residue"],
    )
    assert branch == d["branch"]
    assert d["verdict"].startswith(rs.READINGS[branch])


def test_the_medians_the_decision_read_are_the_0_point_5_quantiles_the_record_publishes() -> None:
    d = _census()
    res, inp = d["residue"], d["decision_inputs"]
    assert res["Q3_top1_share"]["quantiles"]["0.5"] == pytest.approx(inp["median_top1_of_the_residue"])
    assert res["Q4_top3_share"]["quantiles"]["0.5"] == pytest.approx(inp["median_top3_of_the_residue"])
    assert res["Q6_k_min"]["quantiles"]["0.5"] == pytest.approx(inp["median_k_min_of_the_residue"])


def test_the_counts_close_and_the_residue_is_what_the_per_gene_table_says_it_is() -> None:
    d = _census()
    g = json.loads(GENES_FILE.read_text())
    counts = d["Q1_counts"]
    assert counts["recovered_at_the_scale"] + counts["the_residue"] == counts["excluded_basal_not_positive"]
    assert d["residue"]["genes"] == counts["the_residue"]
    assert d["recovered"]["genes"] == counts["recovered_at_the_scale"]
    assert len(g["rows"]) == counts["excluded_basal_not_positive"]
    i = g["columns"].index("in_residue")
    assert sum(1 for r in g["rows"] if r[i]) == counts["the_residue"]


def test_the_residue_is_recomputable_from_the_table_and_every_member_exceeds_the_level_at_the_scale() -> None:
    """A gene is in the residue exactly when its scaled summed drop is not below the gene's level."""
    g = json.loads(GENES_FILE.read_text())
    i_res, i_q, i_sum = (g["columns"].index(c) for c in ("in_residue", "q_after", "sum_d"))
    for row in g["rows"]:
        assert bool(row[i_res]) == (not row[i_q] < 1.0)
        assert row[i_q] == pytest.approx(rs.SCALE * row[i_sum])


def test_every_declared_question_is_answered_for_both_populations_with_its_whole_distribution() -> None:
    d = _census()
    for name in ("residue", "recovered"):
        block = d[name]
        for key in ("Q2_n_activating", "Q3_top1_share", "Q4_top3_share", "Q5_q_after_scaling", "Q8_n_eff"):
            dist = block[key]
            assert set(dist["quantiles"]) == {str(q) for q in rs.QUANTILES}, f"{name}.{key}"
            assert dist["n"] + dist["n_nan"] == block["genes"]


def test_the_cap_curve_covers_the_whole_declared_ladder_and_no_other_level() -> None:
    d = _census()
    per = d["Q7_per_element_cap"]["per_level"]
    assert list(per) == [str(c) for c in rs.CAP_LEVELS]
    assert d["Q7_per_element_cap"]["primary_level"] == rs.CAP_PRIMARY
    for cap, v in per.items():
        assert v["recovered"] + v["still_not_fittable"] == d["residue"]["genes"]
        # CAP_BOUND: a cap fits every gene with fewer than 1/L activating elements, so it is a floor
        assert v["recovered"] >= v["genes_below_the_bound_1_over_L"], f"L={cap} breaks its own bound"


def test_the_cap_counts_recompute_from_the_per_gene_table() -> None:
    d = _census()
    g = json.loads(GENES_FILE.read_text())
    i_res = g["columns"].index("in_residue")
    for cap in rs.CAP_LEVELS:
        i = g["columns"].index(f"fits_under_{cap}")
        n = sum(1 for r in g["rows"] if r[i_res] and r[i])
        assert n == d["Q7_per_element_cap"]["per_level"][str(cap)]["recovered"]


def test_the_census_re_derived_its_population_and_reports_whether_it_agrees_with_what_came_before() -> None:
    d = _census()
    rec = d["reconciliation"]
    assert rec["eligible_for_a_fit"]["re_derived"] == d["Q1_counts"]["eligible_for_a_fit"]
    assert rec["the_residue"]["re_derived"] == d["Q1_counts"]["the_residue"]
    assert rec["agrees"] is True, f"this run disagrees with the committed counts: {rec}"


def test_the_identity_and_the_cap_cross_check_held_on_every_gene_of_the_run() -> None:
    d = _census()
    checks = d["identity_checks"]
    assert checks["max_abs_q_after_minus_scale_times_q_before"] < 1e-9
    assert checks["max_abs_capped_sum_minus_combine2_q"] < 1e-9
    assert checks["cap_verdict_disagreements_with_combine2"] == 0


def test_the_census_opened_no_rna_and_claims_nothing_about_expression_or_any_form() -> None:
    d = _census()
    assert d["rna_read"] == "none" and d["alphagenome_requests"] == 0
    assert d["scale_is_refitted"] is False
    assert "not a claim about expression and not a proposed form" in d["what_this_is_not"]
    assert d["form"] == "mean"
