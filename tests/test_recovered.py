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
