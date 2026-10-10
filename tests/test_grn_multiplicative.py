# SPDX-License-Identifier: Apache-2.0
"""The opt-in `multiplicative` activator rule of the gene-network runtime (lane-saturate, 2026-10-10).

`multiplicative` is A = exp(G * W), W the sum of a gene's activator terms and G its
`multiplicative_gain` attribute: the log-link model of Zhou et al. 2024 (Cell Genomics 100672), in which
enhancer effects in K562 combined multiplicatively. Removing or holding one activator divides the drive
by that activator's own factor whatever the others do, so a gene's single-removal folds may sum past its
whole level without contradiction. It is OPT-IN: a run gets it only by
`NetworkRuntime(..., activator_combination="multiplicative")`; the default stays the mean, competition
is unchanged, and the four rules of the 2026-09-28 choice test are untouched. Engine-only: this file
imports nothing from the application.
"""

from __future__ import annotations

import math

import pytest

from genomeos.lang import parse
from genomeos.runtime import grn
from genomeos.runtime.grn import NetworkRuntime

THREE = """
module net.multiplicative.three
gene g { max: 10; basal: 0.0 }
protein S1 {}
protein S2 {}
protein S3 {}
protein R1 {}
rule S1 activates g { strength: 1.0; threshold: 1.0; hill: 2 }
rule S2 activates g { strength: 0.5; threshold: 1.0; hill: 2 }
rule S3 activates g { strength: 0.25; threshold: 1.0; hill: 2 }
rule R1 inhibits g { strength: 0.4; threshold: 1.0; hill: 2 }
"""
HOURS, DT = 40.0, 0.05
H = 0.5  # H(1) at K = 1, n = 2
PRESENT = {"S1": 1.0, "S2": 1.0, "S3": 1.0, "R1": 1.0}


def _level(clamp: dict[str, float], gain: float | None = None, rule: str | None = "multiplicative") -> float:
    m = parse(THREE)
    if gain is not None:
        m.entities["g"].attrs[grn.MULTIPLICATIVE_GAIN] = gain
    rt = NetworkRuntime(m, strict=True, activator_combination=rule)
    return rt.run(HOURS, DT, clamp=clamp, record_every=10_000).final()["g.mRNA"]


def test_multiplicative_is_opt_in_and_nothing_else_moves():
    assert grn.ACTIVATOR_COMBINATION == "mean"
    assert grn.ACTIVATOR_RULES == ("mean", "sum_capped", "max", "or")
    assert grn.OPT_IN_ACTIVATOR_RULES == ("competition", "multiplicative")
    assert "multiplicative" not in grn.ACTIVATOR_RULES
    assert NetworkRuntime(parse(THREE)).activator_combination is None
    assert NetworkRuntime(parse(THREE)).gain == {}
    # the default run is the mean exactly as before, and the gain attribute is not read by it
    rep = 1.0 - 0.4 * H
    mean = 10 * (1.0 * H + 0.5 * H + 0.25 * H) / 3 * rep
    assert _level(PRESENT, rule=None) == pytest.approx(mean, rel=1e-9)
    assert _level(PRESENT, gain=7.0, rule=None) == pytest.approx(mean, rel=1e-9)
    # competition still reads its own attribute and nothing else
    assert grn.combine_activators([0.3], "competition", 2.0) == pytest.approx(0.3 / 2.3)


def test_combine_activators_multiplicative_values():
    assert grn.combine_activators([0.5, 0.25], "multiplicative") == pytest.approx(math.exp(0.75))
    assert grn.combine_activators([0.5, 0.25], "multiplicative", gain=4.0) == pytest.approx(math.exp(3.0))
    assert grn.combine_activators([], "multiplicative", gain=4.0) == 1.0  # no activator: as every rule
    for bad in (0.0, -1.0, math.inf, math.nan):
        with pytest.raises(ValueError, match="finite gain"):
            grn.combine_activators([0.2], "multiplicative", gain=bad)
    with pytest.raises(ValueError, match="multiplicative"):
        grn.combine_activators([0.2], "median")  # the refusal names the opt-in rule too


@pytest.mark.parametrize("gain", [0.5, 1.0, 6.0])
def test_a_multiplicative_run_reaches_its_closed_form(gain):
    rep = {True: 1.0 - 0.4 * H, False: 1.0}
    s = {"S1": 1.0, "S2": 0.5, "S3": 0.25}
    for held in ((), ("S1",), ("S2", "S3"), ("S1", "S2", "S3"), ("R1",), ("S1", "R1")):
        clamp = {**PRESENT, **dict.fromkeys(held, 0.0)}
        w = sum(v * H for k, v in s.items() if k not in held)
        expect = 10 * math.exp(gain * w) * rep["R1" not in held]
        assert _level(clamp, gain=gain) == pytest.approx(expect, rel=1e-9)


def test_each_removal_is_its_own_factor_so_single_removals_can_sum_past_the_whole():
    """The property the rule exists for: a fold per element, whatever the others' states, so the
    single-removal drops 1 - RHO_i are not bounded by the whole level as under an additive form."""
    gain = 6.0
    full = _level(PRESENT, gain=gain)
    rho = {el: _level({**PRESENT, el: 0.0}, gain=gain) / full for el in ("S1", "S2", "S3")}
    for el, s in (("S1", 1.0), ("S2", 0.5), ("S3", 0.25)):
        assert rho[el] == pytest.approx(math.exp(-gain * s * H), rel=1e-9)
        # the same fold with another element already held
        other = {**PRESENT, "S1" if el != "S1" else "S2": 0.0}
        assert _level({**other, el: 0.0}, gain=gain) / _level(other, gain=gain) == pytest.approx(
            rho[el], rel=1e-9
        )
    assert sum(1 - r for r in rho.values()) > 1.0  # 0.95 + 0.78 + 0.53: past the whole


def test_a_gene_without_the_attribute_takes_the_default_gain_and_bad_gains_are_refused():
    expect = 10 * math.exp(grn.DEFAULT_MULTIPLICATIVE_GAIN * 1.75 * H) * (1.0 - 0.4 * H)
    assert _level(PRESENT) == pytest.approx(expect, rel=1e-9)
    for bad in (0.0, -2.0, math.nan, math.inf):
        m = parse(THREE)
        m.entities["g"].attrs[grn.MULTIPLICATIVE_GAIN] = bad
        with pytest.raises(ValueError, match="must be a finite number > 0"):
            NetworkRuntime(m, activator_combination="multiplicative")
