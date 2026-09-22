# SPDX-License-Identifier: AGPL-3.0-or-later
"""The per-screen offset, checked on pairs whose truth is constructed rather than measured."""

from __future__ import annotations

import math

from genomeos.attribution import crispri
from genomeos.attribution import target_calibration as tc


def pair(dataset: str, cell: str, regulated: bool, drop: float = 0.0, dist: float = 10000.0):
    p = crispri.Pair(
        chrom="chr1",
        start=0,
        end=100,
        gene="G",
        cell=cell,
        dataset=dataset,
        distance=dist,
        dhs=0.0,
        h3k27ac=0.0,
        regulated=regulated,
    )
    p.covered = True
    p.features = {
        "deletion_drop": drop,
        "top_target": 0.0,
        "log_tss_distance": math.log(dist),
        "node_target": 0.0,
        "scored": 1.0,
        **tc.class_features("dELS"),
    }
    return p


def screen(dataset: str, n: int, k: int, cell: str = "K562", drop: float = 0.0):
    """n pairs of one screen, k of them regulated, all with the same features."""
    return [pair(dataset, cell, i < k, drop) for i in range(n)]


def test_a_screen_is_a_dataset_in_a_cell_line():
    a, b = pair("Nasser2021", "K562", True), pair("Nasser2021", "GM12878", True)
    assert tc.screen_key(a) != tc.screen_key(b)


def test_the_size_rule_is_read_off_pairs_and_positives_and_nothing_else():
    rows = screen("big", 200, 30) + screen("few_positives", 200, 5) + screen("tiny", 20, 15)
    census = tc.screen_census(rows)
    assert census["big|K562"]["own_term"] is True
    assert census["few_positives|K562"]["own_term"] is False  # 5 positives, under 10
    assert census["tiny|K562"]["own_term"] is False  # 20 pairs, under 100
    strata = tc.screen_strata(rows)
    assert strata["big|K562"] == "big|K562"
    assert strata["tiny|K562"] == tc.SMALL_SCREEN
    assert strata["few_positives|K562"] == tc.SMALL_SCREEN


def test_the_small_screens_pool_into_one_stratum_with_their_pooled_rate():
    rows = screen("big", 200, 30) + screen("a", 20, 20) + screen("b", 20, 0)
    strata = tc.screen_strata(rows)
    rates = tc.screen_base_rates(rows, strata)
    assert rates["big|K562"] == 30 / 200
    assert rates[tc.SMALL_SCREEN] == 20 / 40  # 20 of 40, not 1.0 and not 0.0


def test_the_offset_spends_no_parameter_so_a_separated_screen_stays_finite():
    """Six of six regulated would send a free intercept to infinity; the offset is a constant."""
    rows = screen("big", 200, 20) + screen("sep", 6, 6)
    strata = tc.screen_strata(rows)
    rates = tc.screen_base_rates(rows, strata)
    off = tc.screen_offsets(rows, strata, rates, 0.5)
    w = tc.fit_with_offset(rows, ("deletion_drop",), off)
    assert all(math.isfinite(v) for v in w)
    p = tc.predict_with_offset(w, rows, ("deletion_drop",), off)
    assert all(0.0 < v < 1.0 for v in p)


def test_the_offset_carries_the_level_and_the_shape_carries_nothing_of_it():
    """Two screens, identical features, different base rates: only the offset may separate them."""
    rows = screen("low", 400, 20, drop=0.3) + screen("high", 400, 200, drop=0.3)
    strata = tc.screen_strata(rows)
    rates = tc.screen_base_rates(rows, strata)
    off = tc.screen_offsets(rows, strata, rates, 0.5)
    w = tc.fit_with_offset(rows, tc.SWEEP_FEATURES, off)
    p = tc.predict_with_offset(w, rows, tc.SWEEP_FEATURES, off)
    low = sum(p[:400]) / 400
    high = sum(p[400:]) / 400
    assert abs(low - 0.05) < 0.02, low
    assert abs(high - 0.50) < 0.05, high


def test_a_residual_shift_of_zero_means_the_base_rate_was_the_whole_story():
    labels = [i < 50 for i in range(200)]
    p = [0.25] * 200
    assert abs(tc.residual_shift(p, labels)) < 1e-3
    assert tc.residual_shift([0.10] * 200, labels) > 1.0


def test_the_bin_count_follows_the_screen_size_and_not_the_result():
    assert tc.screen_bins(tc.SCREEN_BINS_FULL) == tc.BINS
    assert tc.screen_bins(tc.SCREEN_BINS_FULL - 1) == tc.SMALL_BINS


def test_smoothing_keeps_a_screen_of_seven_of_seven_off_the_asymptote():
    assert 0.9 < tc.smoothed(7, 7) < 1.0
    assert 0.0 < tc.smoothed(0, 7) < 0.1
    assert math.isfinite(tc.logit(tc.smoothed(1, 1)))


def test_leave_one_screen_out_never_lets_a_screen_see_its_own_labels_in_the_shape():
    """The held-back screen's only contribution to its own prediction is its base rate."""
    rows = screen("a", 400, 40, drop=0.4) + screen("b", 400, 20, drop=0.2) + screen("c", 400, 100, drop=0.6)
    out = tc.leave_one_screen_out(rows, ("deletion_drop",))
    assert sorted(out["screens_with_their_own_term"]) == ["a|K562", "b|K562", "c|K562"]
    for s, v in out["per_screen"].items():
        assert v["pooled_model"]["pairs"] == 400, s
        assert v["own_base_rate"] != v["base_rate_of_the_other_screens"], s


def test_the_verdict_is_the_conjunction_of_the_four_registered_clauses():
    rows = screen("a", 400, 40, drop=0.4) + screen("b", 400, 20, drop=0.2)
    out = tc.leave_one_screen_out(rows, ("deletion_drop",))
    assert out["screens_improving_required"] == tc.SCREENS_IMPROVING_REQUIRED
    assert out["weighted_ece_fall_required"] == tc.ECE_FALL_REQUIRED
    assert out["residual_shift_bar"] == tc.RESIDUAL_SHIFT_BAR
    # two screens cannot reach the registered four, so the verdict must be failed whatever it reads
    assert out["verdict"] == "failed"


def test_the_registration_says_a_band_cannot_be_quoted_for_an_untested_target():
    text = tc.PREREGISTERED_PREVALENCE
    assert "IT CANNOT BE QUOTED" in text
    assert "NOT ONE SCREEN NAME IS IN BOTH TABLES" in text
    assert "WHAT SHOULD BE QUOTED INSTEAD" in text


def test_the_offset_solver_reproduces_the_shipped_one_when_every_offset_is_zero():
    """The comparator must differ from the offset model in the offset and in nothing else."""
    rows = screen("a", 300, 40, drop=0.4) + screen("b", 300, 20, drop=0.2)
    shipped = tc.fit(rows, tc.SWEEP_FEATURES)
    same = tc.fit_with_offset(rows, tc.SWEEP_FEATURES, [0.0] * len(rows))
    assert all(abs(a - b) < 1e-6 for a, b in zip(shipped, same, strict=True)), (shipped, same)


def test_the_offset_fit_is_calibrated_in_sample_on_a_near_separable_problem():
    """An undamped Newton step diverges here; the line search is what keeps the fit finite."""
    rows = screen("a", 400, 380, drop=0.9) + screen("b", 400, 20, drop=0.0)
    strata = tc.screen_strata(rows)
    rates = tc.screen_base_rates(rows, strata)
    off = tc.screen_offsets(rows, strata, rates, 0.5)
    w = tc.fit_with_offset(rows, tc.SWEEP_FEATURES, off)
    assert all(abs(v) < 1e4 for v in w), w
    p = tc.predict_with_offset(w, rows, tc.SWEEP_FEATURES, off)
    assert abs(sum(p) / len(p) - 400 / 800) < 0.02
