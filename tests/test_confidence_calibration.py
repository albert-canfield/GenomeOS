# SPDX-License-Identifier: AGPL-3.0-or-later
"""The reliability table for the project's own confidences: the rule, the two zeros, the refusals.

Everything here is synthetic. What is pinned is not a number but the discipline: that a band below
the bar is never judged, that the two agreement denominators cannot collapse into one, that an assay
which cannot disagree gets no table, and that the failure patterns were named before any rate existed.
"""

from __future__ import annotations

import itertools

import pytest

from genomeos.attribution import confidence_calibration as cc
from genomeos.attribution.measured import AGREES, DISAGREES, NOT_TESTED

CHROM = "chrT"


COUNTER = itertools.count()


def element(conf, verdicts=None, eligible=("crispri",), gc=0.5, tss=10_000.0, length=300.0):
    return cc.Element(
        id=f"e{next(COUNTER)}",
        chrom=CHROM,
        stated_confidence=conf,
        length=length,
        gc=gc,
        nearest_coding_tss=tss,
        eligible=frozenset(eligible),
        verdicts=dict(verdicts or {}),
    )


def population(spec):
    """(confidence, crispri verdict or None, how many) triples, flattened into elements."""
    out = []
    for conf, verdict, count in spec:
        for _ in range(count):
            out.append(element(conf, {"crispri": verdict} if verdict else None))
    return out


# --- the number under test ---------------------------------------------------------------------------
def test_stated_confidence_is_the_compilers_formula_and_not_a_rereading():
    """The confidence is recomputed from the prediction exactly as compile.py writes it."""
    assert cc.stated_confidence({"predicted_coding": {"gene": "A", "confidence": 0.4321}}) == 0.43
    # the cap bites, and it is the compiler's cap
    assert cc.stated_confidence({"predicted_coding": {"gene": "A", "confidence": 0.99}}) == cc.PREDICTED_CAP
    # the floor bites too
    assert cc.stated_confidence({"predicted_coding": {"gene": "A", "confidence": 0.0}}) == 0.05
    # no confidence: the effect size, as the compiler falls back to it
    assert cc.stated_confidence({"predicted_coding": {"gene": "A", "log2_fold_change": -0.33}}) == 0.33
    # an element with no named target states nothing, and None is not zero
    assert cc.stated_confidence({"predicted_coding": {"log2_fold_change": -0.33}}) is None
    assert cc.stated_confidence({}) is None


def test_the_bands_are_not_this_lanes_to_choose():
    """They are target_calibration's, imported unchanged. A band tuned to a rate is not a band."""
    from genomeos.attribution.target_calibration import CONFIDENCE_BANDS

    assert cc.BANDS is CONFIDENCE_BANDS
    assert cc.band_of(0.1) == "0.1-0.25"  # half-open, [lo, hi), as that module reads its own
    assert cc.band_of(0.25) == "0.25-0.5"
    assert cc.band_of(0.7) == "0.5-0.75"


# --- the two denominators ----------------------------------------------------------------------------
def test_the_narrow_denominator_excludes_the_unasked_question():
    """`predicted_gene_not_tested` is not a wrong answer, and it must not be counted as one."""
    pop = population([(0.2, AGREES, 5), (0.2, DISAGREES, 3), (0.2, NOT_TESTED, 100)])
    narrow = cc.table(pop, "crispri", narrow=True)
    wide = cc.table(pop, "crispri", narrow=False)
    assert narrow["elements"] == 8 and narrow["agrees"] == 5
    assert wide["elements"] == 108 and wide["agrees"] == 5
    assert narrow["observed_overall"] == 0.625
    assert wide["observed_overall"] == round(5 / 108, 4)
    # two names, never one number
    assert narrow["denominator_name"] != wide["denominator_name"]


def test_an_assay_report_keeps_the_wide_key_even_when_it_coincides():
    """For lentiMPRA the two denominators are the same rows; the key stays so the narrow rate is
    never the only rate on the page."""
    pop = [element(0.2, {"lentimpra": AGREES}, eligible=("lentimpra",)) for _ in range(40)]
    rep = cc.assay_report(pop, "lentimpra")
    assert rep["denominators_coincide"] is True
    assert "identical_to" in rep["over_all_matched_elements"]


# --- the bar, and that it is not moved ------------------------------------------------------------------
def test_a_band_below_the_bar_is_printed_and_never_judged():
    short = cc.MIN_FOR_A_BAND - 1
    pop = population([(0.2, AGREES, short), (0.6, AGREES, cc.MIN_FOR_A_BAND)])
    tbl = cc.table(pop, "crispri", narrow=True)
    rows = {r["band"]: r for r in tbl["rows"] if r["measured_by_this_assay"]}
    assert rows["0.1-0.25"]["judged"] is False and rows["0.1-0.25"]["measured_by_this_assay"] == short
    assert rows["0.5-0.75"]["judged"] is True
    assert tbl["bands_populated"] == 2 and tbl["bands_judged"] == 1


def test_too_few_judged_bands_refuses_a_verdict_and_names_the_shortfall():
    """A bar missed by one is still a bar. The refusal says by how much rather than moving it."""
    pop = population([(0.2, AGREES, cc.MIN_FOR_A_BAND - 1), (0.6, AGREES, cc.MIN_FOR_A_BAND)])
    v = cc.verdict(cc.table(pop, "crispri", narrow=True), level_interpretable=True)
    assert v["pattern"] == "refused"
    assert v["bands_populated_but_below_the_bar"]["0.1-0.25"]["short_of_the_bar_by"] == 1
    # the offset is still described, because describing is not judging
    assert v["described_not_judged"]["median_gap_over_populated_bands"] is not None


# --- the pre-stated patterns ------------------------------------------------------------------------------
def test_the_patterns_were_named_before_any_rate():
    assert set(cc.PATTERNS) == {
        "uninformative",
        "ordered_but_miscalibrated",
        "calibrated",
        "level_withheld_ordering_only",
    }
    for text in cc.PATTERNS.values():
        assert len(text) > 100  # each says what it would look like, not just what it is called


def _three_bands(rate_low, rate_mid, rate_high, n=400):
    spec = []
    for conf, rate in ((0.15, rate_low), (0.35, rate_mid), (0.65, rate_high)):
        hits = round(n * rate)
        spec += [(conf, AGREES, hits), (conf, DISAGREES, n - hits)]
    return population(spec)


def test_a_flat_table_is_uninformative():
    v = cc.verdict(cc.table(_three_bands(0.3, 0.3, 0.3), "crispri", True), level_interpretable=True)
    assert v["pattern"] == "uninformative"
    assert v["intervals_disjoint"] is False


def test_a_rise_whose_level_is_wrong_is_ordered_but_miscalibrated():
    # the rates rise far apart, and none of them is near the band's stated confidence
    v = cc.verdict(cc.table(_three_bands(0.55, 0.75, 0.95), "crispri", True), level_interpretable=True)
    assert v["pattern"] == "ordered_but_miscalibrated"
    assert v["ordering"] == "rises with the stated confidence"
    assert v["median_signed_offset_observed_minus_stated"] > 0


def test_a_rise_that_lands_on_its_own_bands_is_calibrated():
    v = cc.verdict(cc.table(_three_bands(0.15, 0.35, 0.65), "crispri", True), level_interpretable=True)
    assert v["pattern"] == "calibrated"
    assert v["share_consistent"] >= cc.MIN_SHARE_CONSISTENT


def test_the_level_is_withheld_where_the_assay_is_not_asking_the_compiled_claim():
    """The same table that would be `calibrated` for CRISPRi earns only `ordering` for a reporter."""
    tbl = cc.table(_three_bands(0.15, 0.35, 0.65), "crispri", True)
    assert cc.verdict(tbl, level_interpretable=False)["pattern"] == "level_withheld_ordering_only"


def test_point_estimates_that_rise_inside_overlapping_intervals_do_not_count_as_a_rise():
    """The ordering is judged on the intervals, and the wording says which happened."""
    v = cc.verdict(cc.table(_three_bands(0.30, 0.32, 0.34), "crispri", True), level_interpretable=True)
    assert v["point_estimates_rise"] is True and v["intervals_disjoint"] is False
    assert v["pattern"] == "uninformative"


# --- the assay that cannot disagree ---------------------------------------------------------------------
def test_the_base_level_assay_gets_no_table_at_all():
    pop = [element(0.3, {"satmut": AGREES}, eligible=("satmut",)) for _ in range(50)]
    rep = cc.assay_report(pop, "satmut")
    assert "refused" in rep and "where_the_predicted_gene_was_tested" not in rep
    assert rep["agrees"] == 50  # 1.0 by construction, which is why there is no table


def test_every_assay_says_whether_its_level_may_be_read():
    from genomeos.attribution.measured import ASSAYS

    assert set(cc.WHAT_AGREEMENT_MEANS) == set(ASSAYS)
    interpretable = [a for a, d in cc.WHAT_AGREEMENT_MEANS.items() if d["level_interpretable"]]
    assert interpretable == ["crispri"]  # only the perturbation asks the compiled claim's question


# --- every group prints its covariates and its coverage ---------------------------------------------------
def test_every_band_prints_length_gc_tss_and_three_denominators_in_order():
    pop = population([(0.2, AGREES, 40), (0.2, None, 200), (0.6, AGREES, 40)])
    tbl = cc.table(pop, "crispri", narrow=True)
    for r in tbl["rows"]:
        for key in (
            "compiled_elements_in_the_band",
            "in_this_assays_footprint",
            "measured_by_this_assay",
            "coverage_of_the_band",
            "median_length",
            "median_gc",
            "median_nearest_coding_tss",
        ):
            assert key in r
    low = next(r for r in tbl["rows"] if r["band"] == "0.1-0.25")
    # compiled >= eligible >= measured, and the coverage is the third over the first
    assert low["compiled_elements_in_the_band"] == 240
    assert low["in_this_assays_footprint"] == 240
    assert low["measured_by_this_assay"] == 40
    assert low["coverage_of_the_band"] == pytest.approx(40 / 240, abs=1e-5)


def test_input_presence_says_bought_or_free_before_any_comparison():
    """The stated confidence is on every element; a measurement is on some. One is free, one bought."""
    pop = population([(0.2, AGREES, 30), (0.6, AGREES, 30), (0.2, None, 500)])
    got = cc.carries_information(pop, "crispri", narrow=True)
    presence = got["input_presence"]
    assert presence["a stated confidence on the compiled fact"]["kind"] == "free"
    assert presence["a verdict from crispri over this element"]["kind"] == "bought"


def test_the_comparison_is_refused_below_the_bar_rather_than_run():
    pop = population([(0.2, AGREES, 5), (0.6, AGREES, 5)])
    got = cc.carries_information(pop, "crispri", narrow=True)
    assert "standardised" not in got and "comparison_refused" in got
    assert str(cc.MIN_FOR_A_COMPARISON) in got["comparison_refused"]


def test_the_measured_slice_prints_the_imbalance_and_the_coverage_per_band():
    pop = population([(0.2, AGREES, 30), (0.2, None, 970), (0.6, AGREES, 60), (0.6, None, 40)])
    got = cc.measured_slice(pop)
    assert got["measured_by_any_assay"] == 90 and got["never_measured"] == 1010
    assert "stated_confidence" in got["imbalance"]  # the confidence itself is a covariate here
    assert got["coverage_by_band"]["0.1-0.25"] == pytest.approx(0.03, abs=1e-5)
    assert got["coverage_by_band"]["0.5-0.75"] == pytest.approx(0.6, abs=1e-5)


def test_the_report_holds_no_pooled_rate():
    """Four assays, four events. A single number over them would be the reporter's wearing the
    perturbation's name, and the result must not contain one."""
    pop = population([(0.2, AGREES, 40), (0.6, DISAGREES, 40)])
    got = cc.report(pop, {CHROM: {"compiled_elements_with_a_stated_confidence": len(pop)}})
    assert "no_pooled_rate" in got and got["why_pooling_is_refused"]
    assert "pooled" not in got["per_assay"]
    assert set(got["per_assay"]) == {"crispri", "lentimpra", "vista", "satmut"}
    assert got["scope"] and got["falsifies_transfer"]
