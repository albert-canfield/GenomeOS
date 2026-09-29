# SPDX-License-Identifier: AGPL-3.0-or-later
"""Review R3, measured rates (registered in ff2062a before the rate table was built).

Negatives first: a gene with no measured rate gets a marker and never a number, a response the
measured rate cannot carry with a positive basal rate is reported and not clipped, and a borrowed
median is labelled as borrowed. Then the acceptance fixture: a gene with a measured transcription
rate and a measured mRNA half-life reproduces its registered observation AND lands on the level those
two measurements imply, so that the simulation reads in molecules per cell rather than in a.u.
"""

from __future__ import annotations

import math

import pytest

from genomeos.attribution import bridge
from genomeos.lang.parser import parse
from genomeos.runtime.grn import LN2, NetworkRuntime

K562 = {"cell_type": "K562"}
#: one measured gene: Schwanhausser's mouse NIH 3T3 numbers for the human symbol they are carried to
T_MEASURED = 12.0  # molecules/(cell*h)
HALF_LIFE = 4.0  # h


def _program(lfc: float = -1.0, half_life: float = HALF_LIFE) -> object:
    return parse(
        f"param mrna_half_life = {half_life} h {{ evidence: experimental"
        ' "Schwanhausser et al. 2011, mouse NIH 3T3" }\n'
        "element E1 { class: enhancer }\n"
        "gene GENE1 { symbol: GENE1 }\n"
        f"rule E1 {'activates' if lfc < 0 else 'inhibits'} GENE1 {{ strength: {abs(lfc)};"
        f' when: cell_type = K562; evidence: predicted "AlphaGenome deletion, K562" effect {lfc:+}'
        " log2 fold change on deletion, probability unavailable }\n",
        "t",
    )


def _steady(module, clamp, hours=200.0):
    return NetworkRuntime(module, context=K562, strict=True).run(
        hours=hours, dt=0.01, clamp=clamp, record_every=1000
    )


# ---- what a measured number is, and is not --------------------------------------------------------


def test_the_split_returns_the_measured_total_and_the_observed_fold():
    """b + V*s*h is the measured rate and b/(b + V*s*h) is the fold: no fraction was invented."""
    for rho in (0.5, 0.25, 0.9):
        b, v = bridge.split_measured_rate(rho, T_MEASURED)
        assert b + v * bridge.DECLARED_STRENGTH * bridge.H_INTACT == pytest.approx(T_MEASURED)
        assert b / (b + v * bridge.DECLARED_STRENGTH * bridge.H_INTACT) == pytest.approx(rho)
    for rho in (1.2, 1.9):
        b, v = bridge.split_measured_rate(rho, T_MEASURED)
        intact = b + v * (1.0 - bridge.DECLARED_STRENGTH * bridge.H_INTACT)
        assert intact == pytest.approx(T_MEASURED) and (b + v) / intact == pytest.approx(rho)


def test_a_measured_rate_is_not_a_max_rate():
    """The registration's first claim, as arithmetic: T is the total, the ceiling is above it."""
    _, v = bridge.split_measured_rate(0.25, T_MEASURED)
    assert v > T_MEASURED and bridge.MEASURED_RATE_UNITS[
        "transcription rate (vsr) [molecules/(cell*h)]"
    ].startswith("becomes T, the gene's TOTAL transcription")


def test_a_gene_with_no_measured_rate_gets_a_marker_and_no_number():
    b = bridge.parameterize(_program(), K562, rates={"OTHER": T_MEASURED})
    assert [(i.reason, i.subject) for i in b.unresolved] == [("gene_rate_unmeasured", "GENE1")]
    assert b.strengths == {} and b.rates_used == {} and b.tier == {}
    assert "mouse" in b.unresolved[0].detail


def test_an_inhibition_that_more_than_doubles_the_gene_leaves_no_basal_rate():
    """The registered inhibitor bound: reported as its own reason, never clipped to something positive."""
    b = bridge.parameterize(_program(lfc=1.5), K562, rates={"GENE1": T_MEASURED})  # RHO = 2.83
    assert [(i.reason, i.subject) for i in b.unresolved] == [("rate_split_infeasible", "E1")]
    assert b.strengths == {} and "leaves basal" in b.unresolved[0].detail


def test_a_borrowed_median_is_labelled_and_changes_no_strength():
    """Scale invariance, stated in the registration: a fold is dimensionless and the split is linear."""
    own = bridge.parameterize(_program(), K562, rates={"GENE1": T_MEASURED})
    lent = bridge.parameterize(_program(), K562, rates={}, borrowed_rate=T_MEASURED * 37.0)
    assert own.tier == {"GENE1": "measured"} and lent.tier == {"GENE1": "borrowed_median"}
    ((_, s_own),) = own.strengths.items()
    ((_, s_lent),) = lent.strengths.items()
    assert s_own == pytest.approx(s_lent)  # the rate buys the scale, not the regulation
    assert lent.rates_used["GENE1"][0] != own.rates_used["GENE1"][0]


def test_an_explicit_gene_entry_still_wins_over_a_measured_rate():
    b = bridge.parameterize(
        _program(), K562, genes={"GENE1": {"basal_rate": 1.0, "max_rate": 10.0}}, rates={"GENE1": T_MEASURED}
    )
    assert b.rates_used == {} and b.strengths


# ---- acceptance -----------------------------------------------------------------------------------


@pytest.mark.parametrize("lfc", [-1.0, -0.4, 0.5])
def test_a_measured_gene_reproduces_its_observation_and_lands_on_its_measured_level(lfc):
    """RATE_ACCEPTANCE (1) and (2) in one run: the fold within 1%, and T / delta_m within 1%."""
    b = bridge.parameterize(_program(lfc=lfc), K562, rates={"GENE1": T_MEASURED})
    assert not b.unresolved and b.tier == {"GENE1": "measured"}
    intact = _steady(b.module, {"E1": bridge.INTACT}).final()["GENE1.mRNA"]
    removed = _steady(b.module, {"E1": bridge.REMOVED}).final()["GENE1.mRNA"]
    assert removed / intact == pytest.approx(2.0**lfc, rel=0.01)
    assert intact == pytest.approx(T_MEASURED * HALF_LIFE / LN2, rel=0.01)  # molecules per cell


def test_the_level_follows_the_measured_half_life_and_not_the_default():
    """A half-life is a degradation constant: doubling it doubles the level and moves no fold."""
    a = bridge.parameterize(_program(half_life=4.0), K562, rates={"GENE1": T_MEASURED})
    d = bridge.parameterize(_program(half_life=8.0), K562, rates={"GENE1": T_MEASURED})
    la = _steady(a.module, {"E1": bridge.INTACT}).final()["GENE1.mRNA"]
    ld = _steady(d.module, {"E1": bridge.INTACT}).final()["GENE1.mRNA"]
    assert ld == pytest.approx(2.0 * la, rel=0.01)
    assert next(iter(a.strengths.values())) == pytest.approx(next(iter(d.strengths.values())))


def test_the_fitted_strength_is_the_declared_one_and_not_the_compiled_number():
    b = bridge.parameterize(_program(lfc=-1.0), K562, rates={"GENE1": T_MEASURED})
    ((_, s),) = b.strengths.items()
    assert s == pytest.approx(bridge.DECLARED_STRENGTH) and not math.isclose(s, 1.0 * 0.3)
