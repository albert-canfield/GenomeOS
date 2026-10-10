# SPDX-License-Identifier: Apache-2.0
"""The opt-in `competition` activator rule of the gene-network runtime (lane-combine2, 2026-10-10).

`competition` is A = W / (K + W), W the sum of a gene's activator terms and K its
`promoter_half_saturation` attribute: Bothma et al. 2015 (eLife 4:e07956) Methods eqs. 7-9, one promoter
engaged by at most one enhancer at a time, read for n enhancers of one transcriptional efficiency. It is
OPT-IN: a run gets it only by `NetworkRuntime(..., activator_combination="competition")`, the default
stays the mean, and the four rules of the 2026-09-28 choice test are untouched. Engine-only: this file
imports nothing from the application.
"""

from __future__ import annotations

import math

import pytest

from genomeos.lang import parse
from genomeos.runtime import grn
from genomeos.runtime.grn import NetworkRuntime

TWO = """
module net.competition.two
gene g { max: 10; basal: 1.0 }
protein S1 {}
protein S2 {}
rule S1 activates g { strength: 1.0; threshold: 1.0; hill: 2 }
rule S2 activates g { strength: 0.5; threshold: 1.0; hill: 2 }
"""
HOURS, DT = 40.0, 0.05


def _level(text: str, clamp: dict[str, float], k: float | None = None, rule: str | None = None) -> float:
    m = parse(text)
    if k is not None:
        m.entities["g"].attrs[grn.PROMOTER_HALF_SATURATION] = k
    rt = NetworkRuntime(m, strict=True, activator_combination=rule)
    return rt.run(HOURS, DT, clamp=clamp, record_every=10_000).final()["g.mRNA"]


def test_competition_is_opt_in_and_the_choice_test_is_untouched():
    assert grn.ACTIVATOR_COMBINATION == "mean"
    assert grn.ACTIVATOR_RULES == ("mean", "sum_capped", "max", "or")
    assert "competition" in grn.OPT_IN_ACTIVATOR_RULES
    assert "competition" not in grn.ACTIVATOR_RULES
    assert set(grn.COMBINATION_TEST["cases"]["bothma2015_kni_early"]["contradicts"]) == set(
        grn.ACTIVATOR_RULES
    )
    assert NetworkRuntime(parse(TWO)).activator_combination is None


def test_combine_activators_competition_values():
    assert grn.combine_activators([0.5, 0.5], "competition") == pytest.approx(0.5)
    assert grn.combine_activators([0.3], "competition", 2.0) == pytest.approx(0.3 / 2.3)
    assert grn.combine_activators([], "competition") == 1.0  # no activator: full rate, as every rule
    for bad in (0.0, -1.0):
        with pytest.raises(ValueError, match="half-saturation"):
            grn.combine_activators([0.2], "competition", bad)
    with pytest.raises(ValueError, match="competition"):
        grn.combine_activators([0.2], "median")  # the refusal names the opt-in rule too


def test_weak_engagement_adds_and_strong_engagement_saturates():
    """Bothma 2015's reading: weak pairs add, strong pairs fall short of adding. Under one K the
    pair against the sum of its singles approaches 1 as W/K -> 0 and falls below it as W/K grows."""

    def ratio(t: float, k: float) -> float:
        pair = grn.combine_activators([t, t], "competition", k)
        single = grn.combine_activators([t], "competition", k)
        return pair / (2 * single)

    assert ratio(0.01, 100.0) == pytest.approx(1.0, abs=1e-3)
    assert ratio(0.5, 0.1) < 0.7
    # monotone: adding an activator never lowers the drive (so it cannot give snail's P < max)
    for k in (0.1, 1.0, 10.0):
        for t in (0.0, 0.2, 1.0):
            assert grn.combine_activators([0.4, t], "competition", k) >= grn.combine_activators(
                [0.4], "competition", k
            )


def test_the_default_run_is_the_mean_exactly_as_before():
    present = {"S1": 1.0, "S2": 1.0}
    h = 0.5  # H(1) at K = 1, n = 2
    expect = 1.0 + 10 * (1.0 * h + 0.5 * h) / 2  # basal + max * mean(t); delta_m = 1 / h at the default
    assert _level(TWO, present) == pytest.approx(expect, rel=1e-9)
    assert _level(TWO, present, rule="mean") == pytest.approx(expect, rel=1e-9)
    # an attribute the default rule does not read changes nothing
    assert _level(TWO, present, k=7.0) == pytest.approx(expect, rel=1e-9)


@pytest.mark.parametrize("k", [0.25, 1.0, 4.0])
def test_a_competition_run_reaches_its_closed_form(k):
    h = 0.5
    for clamp, w in (
        ({"S1": 1.0, "S2": 1.0}, 1.0 * h + 0.5 * h),
        ({"S1": 1.0, "S2": 0.0}, 1.0 * h),
        ({"S1": 0.0, "S2": 1.0}, 0.5 * h),
        ({"S1": 0.0, "S2": 0.0}, 0.0),
    ):
        expect = 1.0 + 10 * w / (k + w)
        assert _level(TWO, clamp, k=k, rule="competition") == pytest.approx(expect, rel=1e-9)


def test_a_gene_without_the_attribute_takes_the_default_half_saturation():
    expect = 1.0 + 10 * 0.75 / (grn.DEFAULT_PROMOTER_HALF_SATURATION + 0.75)
    assert _level(TWO, {"S1": 1.0, "S2": 1.0}, rule="competition") == pytest.approx(expect, rel=1e-9)


def test_refusals():
    with pytest.raises(ValueError, match="unknown activator combination"):
        NetworkRuntime(parse(TWO), activator_combination="median")
    m = parse(TWO)
    m.entities["g"].attrs[grn.PROMOTER_HALF_SATURATION] = 0.0
    with pytest.raises(ValueError, match="must be > 0"):
        NetworkRuntime(m, activator_combination="competition")
    m.entities["g"].attrs[grn.PROMOTER_HALF_SATURATION] = math.nan
    with pytest.raises(ValueError, match="must be > 0"):
        NetworkRuntime(m, activator_combination="competition")
