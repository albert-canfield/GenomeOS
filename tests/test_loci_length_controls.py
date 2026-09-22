# SPDX-License-Identifier: AGPL-3.0-or-later
"""The length-matched control arm: synthetic rows only, no network and no model request.

The gate a test like this has to hold is not that the arithmetic runs but that the registration
cannot drift: the decision rule, the strata, the tolerance and the combination rule are asserted
against the text that was committed before the first number existed.
"""

from __future__ import annotations

from genomeos.benchmark import loci_length

L = loci_length  # a short handle: every line here is an assertion about this one module


def row(locus, length, storage, *, values=1, scored=1, gc=0.45, distance=10_000, constrained=0.1):
    in_syntax = 1 if storage else 0
    return L.row_of(
        locus,
        "chr1",
        1_000,
        1_000 + length,
        gc,
        distance,
        constrained,
        {"values": values, "values_in_syntax": in_syntax},
        {"elements_scored": scored},
    )


def test_registration_states_every_number_the_rule_uses():
    text = L.PREREGISTRATION
    for token in ("9 of 17", "15 of 85", ">= 0.20", "<= 0.10", "undecided at this n", "NEVER POOLED"):
        assert token in text, token
    # the revision is written out rather than folded in, and the failed blind comes before the tables
    assert "THE REVISION" in text
    assert text.index("THIS LANE IS NOT BLINDED") < text.index("THE OUTCOMES")
    assert "arm C decides on S1_covariates" in text
    assert L.LENGTH_TOLERANCE_BP == 0
    assert L.CONTROLS_PER_LOCUS == 5
    assert L.PRIMARY_BY_ARM == {"B_existing": "P1_length", "C_length_only": "S1_covariates"}
    assert L.PRIMARY_WITH_COVERAGE_BY_ARM == {
        "B_existing": "P2_length_coverage",
        "C_length_only": "S2_covariates_coverage",
    }
    assert set(L.STRATA) == {
        "P1_length",
        "P2_length_coverage",
        "P3_length_has_values",
        "S1_covariates",
        "S2_covariates_coverage",
    }
    assert L.STRATA["P2_length_coverage"].covariates == ("length", "deletion_scored")


def test_the_two_zeros_are_different_categories():
    assert L.coverage_of(None, None) == "never_looked_syntax"
    assert L.coverage_of(0, 0) == "no_variable_position"
    assert L.coverage_of(7, 0) == "values_none_constrained"
    assert L.coverage_of(7, 2) == "values_in_syntax"
    arm = [row("a", 400, False, values=0), row("a", 400, True), row("a", 400, False, values=4)]
    s = L.summarise_arm("x", arm, "synthetic")
    assert s["coverage"] == {
        "never_looked_syntax": 0,
        "no_variable_position": 1,
        "values_none_constrained": 1,
        "values_in_syntax": 1,
    }
    # a window with nothing to say cannot make the claim, so the conditioned rate has denominator 2
    assert s["storage"] == {"k": 1, "n": 3, "rate": 0.3333}
    assert s["storage_given_a_value_exists"] == {"k": 1, "n": 2, "rate": 0.5}


def test_length_matching_is_checked_from_the_coordinates_not_from_the_prose():
    panel = [row("a", 400, True), row("b", 9_000, False)]
    good = [row("a", 400, False), row("b", 9_000, False)]
    bad = [row("a", 400, False), row("b", 9_400, False)]
    assert L.length_match_check(panel, good)["matched"] is True
    off = L.length_match_check(panel, bad)
    assert off["matched"] is False and off["worst_difference_bp"] == 400


def test_length_only_draw_is_exactly_the_element_length_and_reproducible(monkeypatch):
    class FakeChrom:
        chrom = "chr1"
        length = 50_000_000

        def gc(self, s, e):
            return 0.4

        def nearest_coding(self, pos):
            return ("GENE", 12_345)

    monkeypatch.setattr(L, "keep_out", lambda ch, results_dir=None: [])
    e = L.PANEL[0]
    want = e.element[1] - e.element[0]
    first = L.length_only_windows(FakeChrom(), e)
    again = L.length_only_windows(FakeChrom(), e)
    assert len(first) == L.CONTROLS_PER_LOCUS
    assert all(w.end - w.start == want for w in first)
    assert [(w.start, w.end) for w in first] == [(w.start, w.end) for w in again]


def test_keep_out_blocks_a_window_that_touches_an_excluded_interval():
    ivs = [(100, 200), (5_000, 6_000)]
    starts = [a for a, _b in ivs]
    assert L.blocked(ivs, starts, 150, 160) is True
    assert L.blocked(ivs, starts, 5_900, 7_000) is True
    assert L.blocked(ivs, starts, 300, 400) is False


def test_verdict_calls_an_artefact_only_when_the_excess_collapses():
    panel = [row(f"L{i}", 400, i < 9) for i in range(17)]
    # a control arm that reproduces the panel inside every length bin: no excess anywhere
    same = [row(f"L{i // 5}", 400, (i // 5) < 9) for i in range(85)]
    comps = {"B_existing": L.compare_arms(panel, same), "C_length_only": L.compare_arms(panel, same)}
    v = L.verdict(comps)
    assert v["verdict"] in {"the claim is an artefact of length", "undecided at this n"}
    assert v["expected_before_the_run"] == "undecided at this n"


def test_each_arm_is_judged_on_the_stratification_registered_for_it():
    panel = [row(f"L{i}", 400, i < 12) for i in range(17)]
    controls = [row(f"L{i // 5}", 400, i < 5) for i in range(85)]
    comps = {
        "B_existing": L.compare_arms(panel, controls),
        "C_length_only": L.compare_arms(panel, controls),
    }
    v = L.verdict(comps)
    assert v["decides_on"] == {"B_existing": "P1_length", "C_length_only": "S1_covariates"}
    assert "decides on P1_length" in v["reasons"][0]
    assert "decides on S1_covariates" in v["reasons"][1]


def test_verdict_can_say_survives_when_the_separation_holds_in_every_stratum():
    # not 17/17 against 0/85: at a difference of exactly 1.0 the normal approximation has no
    # standard error left and `compare.difference` returns p 0.5, which is the right refusal
    panel = [row(f"L{i}", 400, i < 15) for i in range(17)]
    controls = [row(f"L{i // 5}", 400, i < 5) for i in range(85)]
    comps = {
        "B_existing": L.compare_arms(panel, controls),
        "C_length_only": L.compare_arms(panel, controls),
    }
    v = L.verdict(comps)
    assert v["verdict"] == "the claim survives"
    assert comps["B_existing"]["P1_length"]["dropped_for_want_of_a_control"] == 0


def test_undecided_when_the_strata_eat_the_targets():
    panel = [row(f"L{i}", 400 * (i + 1), i < 9) for i in range(17)]
    controls = [row("L0", 400, False) for _ in range(85)]
    comps = {
        "B_existing": L.compare_arms(panel, controls),
        "C_length_only": L.compare_arms(panel, controls),
    }
    v = L.verdict(comps)
    assert v["verdict"] == "undecided at this n"


def test_standardisation_is_compare_py_and_carries_the_matched_n():
    panel = [row(f"L{i}", 400, i < 9) for i in range(17)]
    controls = [row(f"L{i // 5}", 400, i < 15) for i in range(85)]
    out = L.compare_arms(panel, controls)["P1_length"]
    assert out["targets"] == 17 and out["controls"] == 85
    assert out["matched"]["a"] is not None
    assert "imbalance" in out and "length" in out["imbalance"]
