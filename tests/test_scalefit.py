# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-scalefit: the scaling, the recovery test, and the guards that keep the scale a constant."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import tracked_paths as tp

from genomeos.attribution import combine2
from genomeos.attribution import scalefit as sf

RULE_FILE = Path("data/results/scalefit_rule.json")
CENSUS_FILE = Path("data/results/scalefit_census.json")


# ==== the scale is a constant, not a fit =================================================================


def test_the_scale_is_m2h_s_number_at_full_precision_from_the_committed_result() -> None:
    d = json.loads(Path("data/results/magcal_result.json").read_text())
    assert d["scales"]["multiplicative"] == sf.SCALE
    assert list(sf.SCALE_CI) == d["scale_cis"]["multiplicative"]
    assert round(sf.SCALE, 4) == sf.SCALE_ROUNDED


def test_no_module_of_this_lane_fits_optimises_or_re_estimates_anything() -> None:
    text = Path("genomeos/attribution/scalefit.py").read_text()
    for banned in ("brentq", "minimize", "curve_fit", "lstsq", "polyfit", "least_squares"):
        assert banned not in text, f"{banned} would make the scale a fit"


def test_the_recovery_threshold_is_the_reciprocal_of_the_scale() -> None:
    assert pytest.approx(4.660025, abs=1e-5) == sf.RECOVERY_THRESHOLD
    assert pytest.approx(1.0, abs=1e-12) == sf.RECOVERY_THRESHOLD * sf.SCALE


# ==== the scaling ========================================================================================


def test_scaling_multiplies_an_activating_drop_and_an_inhibiting_lift_alike() -> None:
    out = sf.scaled_rho({"a": 0.5, "b": 3.0, "c": 1.0})
    assert 1.0 - out["a"] == pytest.approx(sf.SCALE * 0.5)  # drop 0.5 -> 0.5 * SCALE
    assert out["b"] - 1.0 == pytest.approx(sf.SCALE * 2.0)  # lift 2.0 -> 2.0 * SCALE
    assert out["c"] == 1.0  # a null effect stays null and stays excluded


def test_the_ratio_is_the_mean_form_s_q_exactly() -> None:
    rho = {"a": 0.3, "b": 0.4, "c": 2.0}
    assert sf.ratio(rho) == pytest.approx(combine2.fit("mean", rho).q)


# ==== q is homogeneous of degree 1 in the drops, under BOTH forms ========================================


@pytest.mark.parametrize("form", ["mean", "competition"])
@pytest.mark.parametrize(
    "rho",
    [
        {"a": 0.4, "b": 0.5},
        {"a": 0.01, "b": 0.02, "c": 0.03, "d": 0.9},
        {"a": 0.2, "b": 0.3, "c": 2.5},
        {f"e{i}": 0.8 for i in range(12)},
    ],
)
def test_q_scales_linearly_so_the_threshold_is_exact(form: str, rho: dict[str, float]) -> None:
    assert sf.q_is_homogeneous(form, rho) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("form", ["mean", "competition"])
def test_a_gene_recovers_exactly_when_its_unscaled_q_is_below_the_threshold(form: str) -> None:
    for n in range(2, 40):
        rho = {f"e{i}": 1.0 - 0.25 for i in range(n)}  # q_mean = 0.25 n
        r = sf.recovery(form, rho)
        if r.before_reason == "basal_not_positive":
            assert r.recovered == (r.q_before < sf.RECOVERY_THRESHOLD)


def test_a_gene_that_was_never_excluded_is_not_counted_as_recovered() -> None:
    r = sf.recovery("mean", {"a": 0.7, "b": 0.8})  # q = 0.5 < 1: 7d fits it already
    assert r.before_reason == "" and not r.recovered and not r.runnable


def test_a_gene_excluded_for_another_reason_is_not_counted_as_recovered() -> None:
    r = sf.recovery("mean", {"a": 0.9, "b": 1.5})  # inhibitor strength above 1, not basal
    assert r.before_reason == "inhibitor_strength_above_1"
    assert not r.recovered and not r.runnable


def test_a_gene_whose_drops_are_far_too_large_does_not_recover() -> None:
    r = sf.recovery("mean", {f"e{i}": 0.5 for i in range(20)})  # q = 10 > 4.66
    assert r.before_reason == "basal_not_positive" and not r.recovered
    assert r.after_reason == "basal_not_positive"


# ==== the strength tests are scale-invariant: the secondary can only lose, never gain ====================


def test_inhibitor_and_activator_strengths_do_not_move_with_the_scale() -> None:
    rho = {"a": 0.8, "b": 0.9, "c": 1.2}  # q = 0.3, inhibitor lift 0.2 <= q
    a = combine2.fit("mean", rho)
    b = combine2.fit("mean", sf.scaled_rho(rho))
    assert a.ok and b.ok
    for el in rho:
        assert a.s[el] == pytest.approx(b.s[el], abs=1e-12)


# ==== the floor and the verdict ==========================================================================


def test_the_floor_is_more_than_half_of_the_excluded_population() -> None:
    assert sf.FLOOR == 6_251 and sf.FLOOR_SHARE == 0.5
    assert sf.FLOOR > 12_500 * sf.FLOOR_SHARE
    assert sf.STRONG_BAR == 10_000 and sf.STRONG_BAR > sf.FLOOR


def test_the_verdict_is_decided_by_the_floor_alone() -> None:
    assert sf.verdict(sf.FLOOR - 1, 12_500).startswith("the magnitude error does not explain")
    assert sf.verdict(sf.FLOOR, 12_500).startswith("the magnitude error partly explains")
    assert sf.verdict(sf.STRONG_BAR - 1, 12_500).startswith("the magnitude error partly explains")
    assert sf.verdict(sf.STRONG_BAR, 12_500).startswith("the magnitude error explains")


def test_the_prior_exposure_is_declared_and_not_hidden() -> None:
    assert "not blind" in sf.PRIOR_EXPOSURE
    assert "1.8518" in sf.PRIOR_EXPOSURE  # 7d's median, the number that was read


def test_the_transfer_across_forms_is_named_a_limitation() -> None:
    assert any("TRANSFER ACROSS FORMS" in x for x in sf.LIMITATIONS)
    assert any("NO EXPRESSION CLAIM IS MADE HERE" in x for x in sf.LIMITATIONS)


# ==== the committed records ==============================================================================


def test_the_rule_was_committed_and_holds_the_floor_and_no_count() -> None:
    if not RULE_FILE.exists():
        pytest.skip("the registration has not been written yet")
    d = json.loads(RULE_FILE.read_text())
    assert d["floor"] == sf.FLOOR and d["strong_bar"] == sf.STRONG_BAR
    assert d["scale"] == sf.SCALE
    assert "copied from ac5580b" in d["rule"]["nothing_is_fitted"]
    assert d["population"]["per_form"]["mean"]["basal_not_positive"] == 12_500
    assert d["rna_read"] == "none" and d["alphagenome_requests"] == 0
    for key in ("recovered", "runnable", "recovered_count"):
        assert key not in d, "the registration must carry no count"


def test_the_census_result_reports_the_population_the_rule_declared() -> None:
    if not CENSUS_FILE.exists():
        pytest.skip("the census has not been computed yet")
    d = json.loads(CENSUS_FILE.read_text())
    assert d["scale"] == sf.SCALE
    assert d["rna_read"] == "none" and d["alphagenome_requests"] == 0
    assert d["floor"] == sf.FLOOR
    mean = d["per_form"]["mean"]
    assert mean["excluded_basal_not_positive"] == 12_500
    assert 0 <= mean["recovered"] <= mean["excluded_basal_not_positive"]
    assert mean["runnable_after"] <= mean["recovered"]
    assert d["verdict"] == sf.verdict(mean["recovered"], mean["excluded_basal_not_positive"])


def test_the_committed_records_are_present_and_their_absence_is_never_a_skip() -> None:
    """The ARMED guard for this lane's committed records (the brief's rule (c), d427f7d).

    It is added here rather than by editing the `exists()` skips above because committing a result is
    what disarms such a skip, and a lane does not --force past its own earlier commit in this checkout
    (scripts/check_staged.py's own note: "The lanes never force in this project -- the coordinator
    does"). So the skips above are now redundant rather than load-bearing: this test raises on a broken
    checkout whatever they do, and deleting them is a two-line coordinator commit.
    """
    tp.must_be_present(RULE_FILE, was="the registration has not been written yet")
