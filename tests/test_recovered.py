# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-recovered: the population, the fixed scale, and the guards that keep this a registered reading.

Every guard on a record of this lane is `tracked_paths.must_be_present`, never a `pytest.skip`: a skip
keyed on a file the lane itself commits goes inert the moment git tracks it
(`tests/test_inert_skip_guards.py` docstring at `4f2652c`), and converting it afterwards needs a force
the harness denies. So the raising form is written here first, before the record exists.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import tracked_paths as tp

from genomeos.attribution import combine2 as c2
from genomeos.attribution import magcal as mc
from genomeos.attribution import multrna as mr
from genomeos.attribution import recovered as rc
from genomeos.attribution import scalefit as sf

POPULATION_FILE = Path("data/results/recovered_population.json")


# ==== the scale is M2h's constant, and this lane does not refit it ======================================


def test_the_scale_is_m2h_s_number_at_full_precision_and_is_not_refitted() -> None:
    d = json.loads(Path("data/results/magcal_result.json").read_text())
    assert d["scales"]["multiplicative"] == sf.SCALE
    assert rc.CONVENTIONS["scale"] == sf.SCALE
    assert rc.CONVENTIONS["scale_is_refitted"] is False
    assert round(sf.SCALE, 4) == 0.2146


def test_this_lane_s_module_neither_fits_nor_optimises_anything() -> None:
    text = Path("genomeos/attribution/recovered.py").read_text()
    for banned in ("brentq", "minimize", "curve_fit", "lstsq", "polyfit", "least_squares", "fit_scale"):
        assert banned not in text, f"{banned} would make the scale a fit"


# ==== the two forms are different, deliberately =========================================================


def test_the_population_form_is_the_mean_and_the_prediction_form_is_the_multiplicative() -> None:
    assert rc.POPULATION_FORM == "mean"
    assert rc.POPULATION_FORM in c2.FORMS
    assert rc.FORM == mr.FORM
    assert rc.FORM != rc.POPULATION_FORM


# ==== the split is M2h's, carried unchanged =============================================================


def test_held_out_split_is_magcal_s_split_gene_for_gene() -> None:
    symbols = ["OCA2", "HERC2", "TP53", "MYC", "SHH", "GATA1", "FOXP2", "BRCA1"]
    train, held = rc.held_out_split(symbols)
    assert sorted(train + held) == sorted(symbols)
    assert all(mc.is_held_out(s) for s in held)
    assert not any(mc.is_held_out(s) for s in train)


def test_the_split_is_not_vacuous_on_these_symbols() -> None:
    _train, held = rc.held_out_split(["OCA2", "HERC2", "TP53", "MYC", "SHH", "GATA1", "FOXP2", "BRCA1"])
    assert 0 < len(held) < 8, "a split that puts every gene on one side would test nothing"


# ==== recovery is scalefit's, on planted genes ===========================================================


def test_a_planted_gene_just_inside_the_threshold_is_recovered_and_one_outside_is_not() -> None:
    # two activating elements, each dropping the level by `d`: the mean form's q is 2d exactly.
    inside = {"e1": 1.0 - 1.2, "e2": 1.0 - 1.2}  # q = 2.4, under 1/SCALE = 4.6600
    outside = {"e1": 1.0 - 3.0, "e2": 1.0 - 3.0}  # q = 6.0, over it
    ex_in, rec_in, _run_in = rc.is_recovered(inside)
    ex_out, rec_out, _run_out = rc.is_recovered(outside)
    assert (ex_in, rec_in) == (True, True)
    assert (ex_out, rec_out) == (True, False)


def test_a_planted_gene_7d_already_fitted_is_not_in_this_population() -> None:
    fittable = {"e1": 0.8, "e2": 0.9}  # q = 0.3 < 1: 7d fits it, so it was never excluded
    excluded, recovered, runnable = rc.is_recovered(fittable)
    assert (excluded, recovered, runnable) == (False, False, False)


def test_recovery_uses_scalefit_s_own_recovery_and_agrees_with_it() -> None:
    rho = {"e1": 1.0 - 1.2, "e2": 1.0 - 1.2}
    r = sf.recovery(rc.POPULATION_FORM, rho)
    assert rc.is_recovered(rho) == (r.before_reason == "basal_not_positive", r.recovered, r.runnable)


# ==== the overlap is measured, and its words follow the measurement =====================================


def test_overlap_reports_disjoint_subset_and_partial_in_the_right_words() -> None:
    disjoint = rc.overlap({"A", "B"}, {"C"})
    subset = rc.overlap({"A", "B"}, {"A", "B", "C"})
    partial = rc.overlap({"A", "B"}, {"B", "C"})
    assert disjoint["disjoint"] and "DISJOINT" in rc.overlap_sentence(disjoint, "theirs")
    assert subset["mine_is_a_subset_of_theirs"]
    assert "WHOLLY CONTAINED" in rc.overlap_sentence(subset, "theirs")
    assert not partial["disjoint"] and not partial["mine_is_a_subset_of_theirs"]
    assert "OVERLAPS" in rc.overlap_sentence(partial, "theirs")
    assert partial["in_both"] == 1 and partial["mine_only"] == 1


def test_the_brief_s_claim_that_these_genes_were_in_no_previous_reading_is_measured_not_assumed() -> None:
    assert "MEASURED" in rc.OVERLAP_QUESTION
    assert "NOT ASSUMED" in rc.OVERLAP_QUESTION


# ==== the population record =============================================================================


def test_the_population_record_is_there_and_its_harness_agrees() -> None:
    tp.must_be_present(POPULATION_FILE, was="a skip keyed on the population record's absence")
    d = json.loads(POPULATION_FILE.read_text())
    assert d["outcome"] == "POPULATION_ESTABLISHED"
    assert d["population_harness"]["agrees"] is True
    for name, c in d["population_harness"]["counts"].items():
        assert c["same"] is True, f"{name}: re-derived {c['here']}, scalefit {c['scalefit']}"


def test_the_population_record_re_derives_scalefit_s_four_counts() -> None:
    tp.must_be_present(POPULATION_FILE, was="a skip keyed on the population record's absence")
    d = json.loads(POPULATION_FILE.read_text())["re_derived_counts"]
    want = rc.POPULATION_SOURCE_COUNTS
    for key in ("eligible_for_a_fit", "excluded_basal_not_positive", "recovered", "runnable_after"):
        assert d[key] == want[key]


def test_the_population_record_s_symbol_lists_are_consistent_with_its_counts() -> None:
    tp.must_be_present(POPULATION_FILE, was="a skip keyed on the population record's absence")
    d = json.loads(POPULATION_FILE.read_text())
    s = d["symbols"]
    assert len(s["all"]) == d["re_derived_counts"]["runnable_with_a_row_in_7h_s_multiplicative_table"]
    assert len(s["held_out"]) == d["re_derived_counts"]["held_out_of_this_population"]
    assert len(s["training"]) == d["re_derived_counts"]["training_of_this_population"]
    assert sorted(s["held_out"] + s["training"]) == sorted(s["all"])
    assert set(s["held_out"]).isdisjoint(s["training"])
    assert all(mc.is_held_out(x) for x in s["held_out"])
    assert rc.sha256_of_symbols(s["all"]) == d["population_sha256"]
    assert rc.sha256_of_symbols(s["held_out"]) == d["held_out_sha256"]


def test_the_population_record_opened_no_rna_value() -> None:
    tp.must_be_present(POPULATION_FILE, was="a skip keyed on the population record's absence")
    d = json.loads(POPULATION_FILE.read_text())
    assert d["rna_read"].startswith("none")
    assert d["alphagenome_requests"] == 0
    assert d["scale_is_refitted"] is False


def test_the_population_step_reads_only_the_symbol_column_of_7h_s_table() -> None:
    text = Path("scripts/recovered_population.py").read_text()
    body = text.split("def multrna_symbols")[1].split("\ndef ")[0]
    assert 'index("symbol")' in body
    for banned in ('index("P")', 'index("D")', '"P"]', '"D"]'):
        assert banned not in body, f"{banned} would be reading a value in the population step"


def test_the_population_is_NOT_assumed_disjoint_from_m2h_s(  # the brief said it was; it is not
) -> None:
    tp.must_be_present(POPULATION_FILE, was="a skip keyed on the population record's absence")
    d = json.loads(POPULATION_FILE.read_text())["overlap_with_m2h"]
    # whichever way it came out, the record must state it as a measurement with both counts present
    for key in ("all_16509", "m2h_training_8181", "m2h_held_out_8328"):
        ov = d[key]
        assert ov["mine"] > 0 and ov["theirs"] > 0
        assert ov["in_both"] + ov["mine_only"] == ov["mine"]
        assert isinstance(ov["disjoint"], bool)
    assert len(d["sentences"]) == 3


@pytest.mark.parametrize("key", ["population", "population_harness_words", "overlap_question"])
def test_the_population_record_carries_the_registered_words_verbatim(key: str) -> None:
    tp.must_be_present(POPULATION_FILE, was="a skip keyed on the population record's absence")
    d = json.loads(POPULATION_FILE.read_text())
    want = {
        "population": rc.POPULATION,
        "population_harness_words": rc.POPULATION_HARNESS,
        "overlap_question": rc.OVERLAP_QUESTION,
    }[key]
    assert d[key] == want


# ==== step 2: the registration ==========================================================================

REGISTRATION_FILE = Path("data/results/recovered_registration.json")


def test_the_registration_carries_7d_7h_and_m2h_s_endpoint_baseline_control_and_bootstrap_verbatim() -> None:
    assert rc.MEASURED == c2.MEASURED
    assert c2.PREDICTIONS["no_rule_baseline"] == rc.NO_RULE_BASELINE
    assert rc.RNA_SELECTION == c2.RNA_SELECTION
    assert rc.PSEUDOCOUNT == c2.PSEUDOCOUNT
    assert rc.MIN_GENES == c2.MIN_GENES
    assert rc.CONTROL_ALPHA == c2.CONTROL_ALPHA
    assert rc.PERMUTATIONS == c2.PERMUTATIONS
    assert list(rc.CELLS) == list(c2.CELLS)
    assert rc.PREDICTIONS["form"] == mc.PREDICTIONS["form"]
    assert rc.BOOTSTRAP["resamples"] == mc.BOOTSTRAP["resamples"]
    assert rc.BOOTSTRAP["seed"] == mc.BOOTSTRAP["seed"]
    assert rc.BOOTSTRAP["unit"] == mc.BOOTSTRAP["unit"]
    assert mr.FORM_WORDS[mr.FORM] == rc.FORM_WORDS
    assert rc.CARRIED["split"] == mc.SPLIT
    assert rc.TOLERANCE == mc.TOLERANCE


def test_the_interval_is_m2h_s_width_and_the_reason_is_recorded() -> None:
    assert rc.INTERVAL == mc.FORM_INTERVAL
    assert rc.INTERVAL_WORDS == mc.FORM_INTERVAL_WORDS == "98.33% CI"
    assert "conservative choice" in rc.INTERVAL_REASON
    assert "read at one width" in rc.INTERVAL_REASON
    assert "DECIDED NOW" in rc.INTERVAL_REASON


def test_the_four_way_table_is_the_supervisor_s_verbatim_and_names_the_scale_interval_it_uses() -> None:
    assert rc.WHICH_GOVERNS == mc.WHICH_GOVERNS
    assert "calibration removes the rules" in rc.WHICH_GOVERNS
    assert "UNREACHABLE HERE" in rc.WHICH_SCALE_INTERVAL
    # and the reason it is unreachable is a fact about M2h's committed interval, not an opinion
    assert rc.SCALE_CI[0] > 0.0, "a scale interval including 0 would make CALIBRATION_REMOVES reachable"


def test_the_decision_is_magcal_s_decide_and_a_straddle_reads_inconclusive_at_this_scale_interval() -> None:
    lo, hi = rc.SCALE_CI
    assert mc.decide(0.01, 0.05, lo, hi) == "BETTER"
    assert mc.decide(-0.05, -0.01, lo, hi) == "WORSE"
    assert mc.decide(-0.01, 0.05, lo, hi) == "INCONCLUSIVE"
    # the branch this lane declared unreachable IS reachable for a lane that fits its own scale
    assert mc.decide(-0.01, 0.05, -0.1, 0.1) == "CALIBRATION_REMOVES"


def test_the_scale_is_not_refitted_and_the_choice_is_stated_rather_than_omitted() -> None:
    assert rc.PREDICTIONS["calibrated"].count("NOTHING IS FITTED") == 1
    assert "NOT REFITTED ON THIS POPULATION" in rc.NO_REFIT
    assert "CANNOT distinguish" in rc.NO_REFIT


def test_the_exposure_says_plainly_that_the_brief_s_claim_is_false_and_gives_the_measured_counts() -> None:
    for want in ("7,999", "3,990", "4,009", "SUBPOPULATION READING", "THAT IS\n FALSE".replace("\n ", " ")):
        assert want in rc.EXPOSURE, want
    assert "NOT a fresh test of the endpoint" in rc.EXPOSURE
    for want in ("7,999", "4,009", "SUBPOPULATION", "not independent confirmation"):
        assert want in rc.EXPOSURE_CLAUSE, want


def test_the_falsifier_names_worse_as_the_sharper_finding_and_says_it_is_reported_first() -> None:
    assert "upper bound < 0, outcome WORSE" in rc.FALSIFIER
    assert "recovers genes the fit can accept but not genes the model predicts" in rc.FALSIFIER
    assert "REPORTED FIRST AND UNSOFTENED" in rc.DECISION


def test_the_variance_clause_is_filled_from_the_numbers_and_not_written_by_hand() -> None:
    s = rc.variance_clause(0.0434, 2.1524)
    assert "about 2.0% of the variance" in s
    assert "2.1524" in s
    assert rc.variance_clause(-0.0434, 2.1524).startswith("about -2.0%")


def test_the_full_reading_is_magcal_s_sentence_with_this_lane_s_three_clauses_after_it() -> None:
    s = rc.full_reading(
        "BETTER",
        delta=0.0434,
        mean_d2=2.1524,
        held=4009,
        total=8535,
        form=rc.FORM_WORDS,
        d="+0.0434",
        ci=rc.INTERVAL_WORDS,
        lo="+0.03",
        hi="+0.05",
        n=52117,
        g=4009,
        p=0.000999,
        scale="the constant +0.2146",
    )
    assert s.startswith(mc.OUTCOMES["BETTER"].split("{")[0])
    assert mc.PRIOR_EXPOSURE_CLAUSE in s
    assert "4009 genes of the 8535" in s
    assert "SUBPOPULATION" in s
    assert "about 2.0% of the variance" in s


def test_every_bound_field_resolves_and_registered_fields_returns_them_all() -> None:
    fields = rc.registered_fields()
    assert set(fields) == set(rc.BOUND)
    assert fields["which_governs"] == mc.WHICH_GOVERNS
    assert fields["exposure"] == rc.EXPOSURE
    assert fields["scale"] == sf.SCALE


def test_the_registration_record_is_there_and_reproduces_every_bound_field() -> None:
    tp.must_be_present(REGISTRATION_FILE, was="a skip keyed on the registration record's absence")
    d = json.loads(REGISTRATION_FILE.read_text())
    for field, value in rc.registered_fields().items():
        assert d[field] == value, field
    assert sorted(d["bound_fields"]) == sorted(rc.BOUND)


def test_the_registration_record_pins_the_population_of_step_1_by_its_sha256() -> None:
    tp.must_be_present(REGISTRATION_FILE, was="a skip keyed on the registration record's absence")
    reg = json.loads(REGISTRATION_FILE.read_text())["population_pin"]
    pop = json.loads(POPULATION_FILE.read_text())
    assert reg["population_sha256"] == pop["population_sha256"]
    assert reg["held_out_sha256"] == pop["held_out_sha256"]
    assert reg["genes"] == len(pop["symbols"]["all"])
    assert reg["held_out_genes"] == len(pop["symbols"]["held_out"])
    assert reg["held_out_genes"] >= rc.MIN_GENES


def test_the_registration_record_reproduces_m2h_s_scale_and_interval_exactly() -> None:
    tp.must_be_present(REGISTRATION_FILE, was="a skip keyed on the registration record's absence")
    pin = json.loads(REGISTRATION_FILE.read_text())["scale_pin"]
    m2h = json.loads(Path("data/results/magcal_result.json").read_text())
    assert pin["scale"] == m2h["scales"]["multiplicative"]
    assert pin["scale_ci"] == m2h["scale_cis"]["multiplicative"]
    assert pin["scale_interval_excludes_zero"] is True
    assert pin["is_refitted_here"] is False


def test_the_registration_record_opened_no_rna_value() -> None:
    tp.must_be_present(REGISTRATION_FILE, was="a skip keyed on the registration record's absence")
    d = json.loads(REGISTRATION_FILE.read_text())
    assert d["rna_read"].startswith("none")
    assert d["alphagenome_requests"] == 0


def test_the_registration_record_quotes_m2h_and_m2i_word_for_word() -> None:
    tp.must_be_present(REGISTRATION_FILE, was="a skip keyed on the registration record's absence")
    ctx = json.loads(REGISTRATION_FILE.read_text())["context_verbatim"]
    m2h = json.loads(Path("data/results/magcal_result.json").read_text())
    m2i = json.loads(Path("data/results/freshrna_result.json").read_text())
    assert ctx["m2h_multiplicative_reading"] == m2h["per_form"]["multiplicative"]["reading"]
    assert ctx["m2i_verdict_sentence"] == m2i["verdict_sentence"]
    assert (
        ctx["scalefit_verdict"]
        == json.loads(Path("data/results/scalefit_census.json").read_text())["verdict"]
    )


# ==== step 3: the one run ===============================================================================

RESULT_FILE = Path("data/results/recovered_result.json")


def test_the_run_never_fits_a_scale() -> None:
    text = Path("scripts/recovered_run.py").read_text()
    for banned in ("fit_scale", "bootstrap_scale", "brentq", "minimize", "curve_fit", "lstsq", "polyfit"):
        assert banned not in text, f"{banned} would make this lane fit the scale it was told to fix"


def test_the_run_reuses_magcal_s_functions_rather_than_reimplementing_them() -> None:
    text = Path("scripts/recovered_run.py").read_text()
    for needed in ("mc.delta_mse(", "mc.bootstrap_delta(", "mc.permuted_control(", "mc.decide("):
        assert needed in text, needed
    assert "rc.full_reading(" in text


def test_the_result_is_there_and_its_harness_agreed_before_anything_was_read() -> None:
    tp.must_be_present(RESULT_FILE, was="a skip keyed on the result's absence")
    d = json.loads(RESULT_FILE.read_text())
    assert d["harness_check"]["agrees"] is True
    assert d["outcome"] in {"BETTER", "WORSE", "INCONCLUSIVE", "NOT_TESTABLE", "HARNESS_DEFECT"}
    assert d["scale"] == sf.SCALE
    assert d["scale_is_refitted"] is False
    assert d["alphagenome_requests"] == 0


def test_the_result_s_outcome_is_magcal_s_decide_on_its_own_two_intervals() -> None:
    tp.must_be_present(RESULT_FILE, was="a skip keyed on the result's absence")
    d = json.loads(RESULT_FILE.read_text())
    di = d["decision_inputs"]
    lo, hi = di["held_out_ci"]
    slo, shi = di["scale_ci"]
    assert mc.decide(lo, hi, slo, shi) == d["outcome"]
    assert [slo, shi] == list(sf.SCALE_CI)
    assert di["calibration_removes_was_unreachable"] is True
    assert di["which_governs"] == mc.WHICH_GOVERNS


def test_the_result_s_reading_is_the_registered_sentence_for_its_own_outcome() -> None:
    tp.must_be_present(RESULT_FILE, was="a skip keyed on the result's absence")
    d = json.loads(RESULT_FILE.read_text())
    assert d["reading"].startswith(mc.OUTCOMES[d["outcome"]].split("{")[0])
    assert mc.PRIOR_EXPOSURE_CLAUSE in d["reading"]
    assert "SUBPOPULATION" in d["reading"]
    assert "% of the variance" in d["reading"]


def test_the_result_carries_the_variance_share_beside_the_reading_as_m2h_and_m2i_do() -> None:
    tp.must_be_present(RESULT_FILE, was="a skip keyed on the result's absence")
    d = json.loads(RESULT_FILE.read_text())
    p = d["primary"]
    assert p["variance_share"] == pytest.approx(p["delta_mse"] / p["mean_D2"], rel=1e-9)
    assert rc.variance_clause(p["delta_mse"], p["mean_D2"]) in d["reading"]


def test_the_result_reads_the_held_out_half_and_only_the_held_out_half() -> None:
    tp.must_be_present(RESULT_FILE, was="a skip keyed on the result's absence")
    d = json.loads(RESULT_FILE.read_text())
    pop = json.loads(POPULATION_FILE.read_text())
    assert d["primary"]["genes"] == len(pop["symbols"]["held_out"])
    assert d["primary"]["pairs"] == 13 * d["primary"]["genes"]
    assert d["primary"]["genes"] >= rc.MIN_GENES
    assert "descriptive" in d["descriptives_not_verdicts"]["training_side"]["side"]


def test_the_result_s_harness_reproduced_7h_s_delta_mse_and_m2h_s_split() -> None:
    tp.must_be_present(RESULT_FILE, was="a skip keyed on the result's absence")
    h = json.loads(RESULT_FILE.read_text())["harness_check"]
    s4 = h["4_delta_mse_at_scale_1_reproduces_7h"]
    assert abs(s4["recomputed"] - s4["committed"]) <= rc.TOLERANCE
    s6 = h["6_m2h_split_reproduces"]
    assert s6["held_out_genes"] == s6["held_out_genes_committed"] == 8328
    assert s6["held_out_pairs"] == s6["held_out_pairs_committed"] == 108264
    assert h["7_this_population_s_held_out_half_is_held_out_under_m2h_s_split"] is True


def test_the_result_carries_the_exposure_and_quotes_m2h_m2i_and_scalefit_unchanged() -> None:
    tp.must_be_present(RESULT_FILE, was="a skip keyed on the result's absence")
    d = json.loads(RESULT_FILE.read_text())
    assert d["exposure"] == rc.EXPOSURE
    assert d["falsifier"] == rc.FALSIFIER
    assert d["limitations"] == list(rc.LIMITATIONS)
    ctx = d["context_verbatim"]
    m2h = json.loads(Path("data/results/magcal_result.json").read_text())
    assert ctx["m2h_multiplicative_reading"] == m2h["per_form"]["multiplicative"]["reading"]
    assert ctx["m2h_held_out_delta_mse"] == m2h["per_form"]["multiplicative"]["held_out_delta_mse"]
    assert (
        ctx["scalefit_verdict"]
        == json.loads(Path("data/results/scalefit_census.json").read_text())["verdict"]
    )


def test_the_result_s_control_p_is_the_registered_formula() -> None:
    tp.must_be_present(RESULT_FILE, was="a skip keyed on the result's absence")
    d = json.loads(RESULT_FILE.read_text())
    assert 0.0 < d["control_p"] <= 1.0
    assert d["control_p"] >= 1.0 / (1 + rc.PERMUTATIONS["count"])
    assert "REFITS the scale" in d["control_difference"]
