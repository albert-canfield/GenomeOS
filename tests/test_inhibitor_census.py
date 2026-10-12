# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-inhibitor: the census's definitions, its derived thresholds and the guards on the record.

Every test here is predictor-side and runs on PLANTED removal folds or on the committed record itself.
None opens an RNA value, none fits anything, and none asserts a count of real genes except through the
record. The tests that read the committed result recompute its verdict from its own numbers.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest
import tracked_paths

from genomeos.attribution import bridge
from genomeos.attribution import combine2 as c2
from genomeos.attribution import inhibitor_census as ic
from genomeos.attribution import residuemult as rm
from genomeos.attribution import saturate as sa
from genomeos.attribution import scalefit as sf

MODULE = Path("genomeos/attribution/inhibitor_census.py")
SCRIPT = Path("scripts/inhibitor_census_run.py")
RESULT_FILE = Path("data/results/inhibitor_census_result.json")
RESIDUEMULT_FILE = Path("data/results/residuemult_lookup.json")


def planted(rho: dict[str, float], symbol: str = "PLANTED") -> c2.GeneRules:
    """One gene with the folds a test chooses. No real gene and no RNA is involved."""
    els = list(rho)
    return c2.GeneRules(
        gene=symbol,
        symbol=symbol,
        chrom="chrTEST",
        elements=els,
        actions=["activates" if rho[e] < 1.0 else "inhibits" for e in els],
        rho=dict(rho),
        unreadable=[],
        compiled={e: "planted_cell" for e in els},
        states={e: ["open"] for e in els},
    )


# ==== nothing is declared twice and nothing is fitted ====================================================


def test_the_scale_the_occupancy_and_the_tolerance_are_inherited_not_re_declared() -> None:
    assert ic.SCALE == sf.SCALE
    assert ic.SCALE_ROUNDED == sf.SCALE_ROUNDED
    assert ic.SCALE == rm.SCALE
    assert ic.H == bridge.H_INTACT
    assert ic.STRENGTH_TOLERANCE == bridge.STRENGTH_TOLERANCE
    body = MODULE.read_text()
    assert "0.2145911381525326" not in body, "the scale must be inherited, never re-typed"
    assert ic.CONVENTIONS["scale_is_refitted"] is False


def test_no_module_or_script_of_this_lane_fits_optimises_or_re_estimates_anything() -> None:
    for path in (MODULE, SCRIPT):
        body = path.read_text()
        for word in ("curve_fit", "minimize", "least_squares", "fsolve", "brentq", "LinearRegression"):
            assert word not in body, f"{path} must not fit anything, found {word}"


def test_the_lane_declares_it_opens_no_rna_and_the_script_reads_only_the_symbol_column() -> None:
    body = SCRIPT.read_text()
    assert "NO RNA VALUE IS OPENED" in body
    assert 'index("symbol")' in body
    assert '"D"' not in body and '"P"' not in body
    assert ic.CONVENTIONS["rna_read"] == "none"


def test_the_record_proposes_no_change_to_any_fit_and_claims_nothing_about_expression() -> None:
    assert "PROPOSES NO CHANGE TO ANY FIT" in ic.WHAT_THIS_IS_NOT
    assert "no new form" in ic.WHAT_THIS_IS_NOT
    assert "not a claim about expression" in ic.WHAT_THIS_IS_NOT
    assert "makes no claim about expression" in MODULE.read_text()


# ==== the questions are the committed ones ===============================================================


def test_the_seven_questions_are_in_the_module_and_name_what_they_must() -> None:
    assert list(ic.QUESTIONS) == [
        "Q1_re_derivation",
        "Q2_strength_distribution",
        "Q3_how_many_inhibitors",
        "Q4_one_inhibitor_or_several",
        "Q5_where_the_rules_came_from",
        "Q6_outside_the_39",
        "Q7_threshold_or_tail",
    ]
    assert "QUANTILES" in ic.QUESTIONS["Q2_strength_distribution"]
    assert "never as a median alone" in ic.QUESTIONS["Q2_strength_distribution"]
    assert "OUTSIDE THE 39" in ic.QUESTIONS["Q6_outside_the_39"]
    assert "not a distinguished set" in ic.QUESTIONS["Q6_outside_the_39"]
    assert "the more important answer" in ic.QUESTIONS["Q6_outside_the_39"]
    assert "evidence KIND" in ic.QUESTIONS["Q5_where_the_rules_came_from"]
    assert "ONE INHIBITOR OR SEVERAL" in ic.QUESTIONS["Q4_one_inhibitor_or_several"]


def test_the_script_takes_every_question_from_the_module_and_declares_none_of_its_own() -> None:
    body = SCRIPT.read_text()
    for key in ic.QUESTIONS:
        assert f'ic.QUESTIONS["{key}"]' in body, f"{key} must be read from the committed record"
    assert "QUESTIONS = {" not in body, "the questions live in the record, not in the entry script"


# ==== IMPLIED_STRENGTH is saturate.fit's own expression, on planted folds ================================


@pytest.mark.parametrize("rho", [1.1, 1.5, 2.0, 3.0, 9.0])
def test_implied_strength_is_the_number_saturate_fit_itself_puts_on_an_inhibitor(rho: float) -> None:
    gr = planted({"act": 0.5, "inh": rho})
    f = sa.fit(gr.rho)
    s = ic.implied_strength(rho)
    if s <= 1.0:
        assert f.ok, f.reason
        assert f.s["inh"] == pytest.approx(s, abs=1e-15)
    else:
        assert not f.ok and f.reason == ic.REASON


def test_implied_strength_is_nan_for_an_activator_so_it_can_never_be_summarised_as_one() -> None:
    assert math.isnan(ic.implied_strength(0.5))
    assert math.isnan(ic.implied_strength(1.0))


def test_the_mean_form_strength_is_combine2_s_own_and_a_DIFFERENT_quantity_from_the_one_that_refused() -> (
    None
):
    rho = {"a1": 0.6, "a2": 0.7, "inh": 1.2}
    f = c2.fit("mean", rho)
    assert f.ok, f.reason
    assert ic.mean_form_strength(1.2, f.q) == pytest.approx(f.s["inh"], rel=1e-12)
    assert ic.mean_form_strength(1.2, f.q) != pytest.approx(ic.implied_strength(1.2))


def test_above_bound_is_exactly_what_combine2_strength_ok_refuses() -> None:
    for s in (0.5, 1.0, 1.0 + ic.STRENGTH_TOLERANCE / 2, 1.0 + ic.STRENGTH_TOLERANCE * 10, 2.0):
        assert ic.above_bound(s) is (c2._strength_ok(s) is None)


# ==== the two derived thresholds, against the code's own refusal =========================================


def test_the_unscaled_threshold_is_a_removal_fold_of_2_and_the_code_refuses_just_past_it() -> None:
    t = ic.derived_thresholds()
    assert t["fold_above_which_an_inhibitor_is_refused_unscaled"] == pytest.approx(2.0)
    assert sa.fit({"act": 0.5, "inh": 2.0}).ok
    refused = sa.fit({"act": 0.5, "inh": 2.0 * (1 + 1e-6)})
    assert not refused.ok and refused.reason == ic.REASON


def test_the_scaled_threshold_is_1_plus_1_over_the_scale_and_the_code_refuses_just_past_it() -> None:
    t = ic.derived_thresholds()
    edge = t["fold_above_which_an_inhibitor_is_refused_at_the_scale"]
    assert edge == pytest.approx(1.0 + 1.0 / ic.SCALE)
    assert sa.fit(sf.scaled_rho({"act": 0.5, "inh": edge * (1 - 1e-6)}, ic.SCALE)).ok
    refused = sa.fit(sf.scaled_rho({"act": 0.5, "inh": edge * (1 + 1e-6)}, ic.SCALE))
    assert not refused.ok and refused.reason == ic.REASON


def test_the_scaled_threshold_is_the_larger_so_scaling_can_only_relieve_this_refusal() -> None:
    t = ic.derived_thresholds()
    assert (
        t["fold_above_which_an_inhibitor_is_refused_at_the_scale"]
        > t["fold_above_which_an_inhibitor_is_refused_unscaled"]
    )
    for rho in (1.01, 1.5, 2.0, 4.0, 40.0):
        scaled = sf.scaled_rho({"inh": rho}, ic.SCALE)["inh"]
        assert ic.implied_strength(scaled) <= ic.implied_strength(rho) + 1e-15


# ==== the bound checks can fail: planted violations are REPORTED =========================================


def test_a_planted_gene_refused_for_the_reason_has_an_offending_element_and_the_check_passes() -> None:
    gr = planted({"act": 0.5, "inh": 40.0})
    els = ic.element_rows(gr, {}, c2.fit("mean", gr.rho).q)
    row = ic.gene_row(gr, els, in_residue=True, excluded_by_the_mean=True)
    row.elements = els
    assert row.reason_scaled == ic.REASON
    assert row.n_offending_scaled == 1
    assert ic.gene_bound_violations(row) == []


def test_PLANTED_a_row_whose_offender_count_disagrees_with_the_code_s_refusal_is_REPORTED() -> None:
    gr = planted({"act": 0.5, "inh": 40.0})
    els = ic.element_rows(gr, {}, c2.fit("mean", gr.rho).q)
    row = ic.gene_row(gr, els, in_residue=True, excluded_by_the_mean=True)
    row.elements = els
    row.n_offending_scaled = 0  # the census quantity, planted to disagree with saturate.fit's refusal
    bad = ic.gene_bound_violations(row)
    assert any("disagree" in b for b in bad), bad


def test_PLANTED_an_element_offending_at_the_scale_but_not_unscaled_is_REPORTED() -> None:
    gr = planted({"act": 0.5, "inh": 40.0})
    els = ic.element_rows(gr, {}, c2.fit("mean", gr.rho).q)
    row = ic.gene_row(gr, els, in_residue=True, excluded_by_the_mean=True)
    row.elements = els
    els[0].offends_unscaled = False  # planted: the subset bound (d) violated
    bad = ic.gene_bound_violations(row)
    assert any("offends at the scale but not unscaled (d)" in b for b in bad), bad


def test_element_rows_holds_only_the_inhibiting_elements_and_carries_the_rule_s_provenance() -> None:
    gr = planted({"a1": 0.4, "a2": 0.8, "inh": 3.0})
    els = ic.element_rows(gr, {"inh": ("predicted", "evidence note, log2 fold change", "src")}, 0.8)
    assert [e.element for e in els] == ["inh"]
    assert els[0].kind == "predicted" and els[0].source == "src"
    assert els[0].how == "evidence note, log2 fold change"
    assert els[0].compiled_cell == "planted_cell"
    assert els[0].rho_scaled == pytest.approx(1.0 + ic.SCALE * 2.0)


# ==== the branch rule is the committed one and the numbers do not choose it ==============================


def test_the_branch_is_decided_by_the_near_band_counts_alone() -> None:
    assert ic.branch(10, 10) == "threshold_is_doing_the_work"
    assert ic.branch(11, 10) == "threshold_is_doing_the_work"
    assert ic.branch(9, 10) == "a_separated_tail"
    assert ic.NEAR_BAND == 1.2
    assert "NEAR-BAND COUNTS ALONE" in ic.IMPLICATION_RULE
    assert "does not choose it" in ic.IMPLICATION_RULE


@pytest.mark.parametrize("which", ["threshold_is_doing_the_work", "a_separated_tail"])
def test_both_branch_wordings_fill_from_the_same_numbers_so_neither_is_written_later(which: str) -> None:
    text = ic.reading(which, near_below=1, near_band=1.2, above=2, outside_genes=3, residue=4)
    assert "{" not in text and "}" not in text
    assert ic.BRANCHES[which].split("{")[0][:20] in text


def test_the_threshold_wording_says_the_threshold_and_not_the_biology() -> None:
    assert "THE THRESHOLD IS DOING THE WORK, NOT THE BIOLOGY" in ic.BRANCHES["threshold_is_doing_the_work"]


# ==== distributions are reported as distributions ========================================================


def test_distribution_reports_min_max_mean_and_every_committed_quantile() -> None:
    d = ic.distribution([float(x) for x in range(1, 101)])
    assert d["n"] == 100 and d["min"] == 1.0 and d["max"] == 100.0
    assert d["mean"] == pytest.approx(50.5)
    assert set(d["quantiles"]) == {f"{q:g}" for q in ic.QUANTILES}
    assert d["quantiles"]["0"] == 1.0 and d["quantiles"]["1"] == 100.0
    assert ic.distribution([]) == {"n": 0}


def test_the_quantile_set_includes_both_ends_so_no_distribution_is_a_median_alone() -> None:
    assert 0.0 in ic.QUANTILES and 1.0 in ic.QUANTILES
    assert len(ic.QUANTILES) >= 5


def test_counts_of_gives_the_exact_value_counts_of_a_small_integer_distribution() -> None:
    assert ic.counts_of([1, 1, 2, 5]) == {"1": 2, "2": 1, "5": 1}
    assert ic.counts_of([]) == {}


# ==== the limitations and the identity convention are in the record ======================================


def test_the_record_states_that_the_compiled_rules_are_predicted_and_that_strength_is_not_observable() -> (
    None
):
    joined = " ".join(ic.LIMITATIONS)
    assert "PREDICTED" in joined
    assert "not observable" in joined
    assert "NEAR_BAND of 1.2 is a choice" in joined


def test_the_match_field_is_the_same_identity_lane_residuemult_fixed() -> None:
    tracked_paths.must_be_present(RESIDUEMULT_FILE)
    assert ic.MATCH_FIELD == "symbol"
    assert json.loads(RESIDUEMULT_FILE.read_text())["match_field"] == ic.MATCH_FIELD


def test_the_committed_residuemult_result_puts_every_one_of_the_neither_genes_on_this_one_reason() -> None:
    tracked_paths.must_be_present(RESIDUEMULT_FILE)
    d = json.loads(RESIDUEMULT_FILE.read_text())
    assert d["informative_under_neither_by_reason"] == {ic.REASON: d["counts"]["informative_under_neither"]}
