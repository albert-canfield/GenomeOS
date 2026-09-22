# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Stage 2's gate (b): the derivation that closed its primary, tested rather than asserted.

The amendment in `genomeos/runtime/abundance_gate.py` claims the pool layer hands every gene one
global factor, so its prediction differs from the one-to-one baseline by a constant and the
registered primary cannot move. Three of these tests are the claim; two are the CONTRAST that keeps
the claim from being vacuous — a policy that does differ per gene, and the saturable form the
amendment names as the instrument, both of which this same machinery detects.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genomeos.ir import Allocation, Cost, Pool
from genomeos.runtime.abundance_gate import AMENDMENT, DERIVATION_TOLERANCE_LOG10, PRE_REGISTRATION
from genomeos.runtime.economy import Economy

RESULT = Path(__file__).resolve().parents[1] / "data" / "results" / "abundance_gate.json"


def economy(policy: str, order: list[str] | None = None, n: int = 200) -> tuple[Economy, dict]:
    """n genes, all drawing the same two pools, with demand spanning four decades."""
    genes = [f"g{i:04d}" for i in range(n)]
    econ = Economy(
        pools={
            "Ribosomes": Pool(id="Ribosomes", kind="pool", size=1.0e4),
            "ATP": Pool(id="ATP", kind="pool", size=1.0e7),
        },
        costs={
            g: [
                Cost(pool="Ribosomes", amount=1.0, per="chain"),
                Cost(pool="ATP", amount=4.0, per="residue"),
            ]
            for g in genes
        },
        allocations={
            "Ribosomes": Allocation(
                id="a1", kind="allocation", pool="Ribosomes", policy=policy, order=order or []
            ),
            "ATP": Allocation(id="a2", kind="allocation", pool="ATP", policy=policy, order=order or []),
        },
        lengths={g: 100.0 + 5.0 * i for i, g in enumerate(genes)},
    )
    work = {g: 10.0 ** (i / (n / 4.0)) for i, g in enumerate(genes)}
    return econ, work


@pytest.mark.parametrize("policy", ["proportional", "competitive"])
def test_a_shared_pool_hands_every_gene_the_same_factor(policy: str) -> None:
    """The amendment's whole claim, and it holds for `competitive` as much as for `proportional`."""
    econ, work = economy(policy)
    factors = econ.scale(work)
    assert min(factors.values()) < 1.0, "the pool must be short or the test proves nothing"
    assert len({round(v, 12) for v in factors.values()}) == 1


def test_competitive_differs_from_proportional_only_in_the_constant() -> None:
    """Reported as a runtime finding: `competitive` saturates per POOL, not per demander."""
    prop, work = economy("proportional")
    comp, _ = economy("competitive")
    fp = prop.scale(work)
    fc = comp.scale(work)
    ratios = {round(fc[g] / fp[g], 9) for g in work}
    assert len(ratios) == 1, "two gene-blind policies must differ by one constant, not per gene"
    assert ratios != {1.0}, "and they must not be the same constant, or the contrast is empty"


def test_priority_is_the_one_policy_that_is_not_gene_blind() -> None:
    """The contrast. Without it, the test above could pass on machinery that cannot see a difference."""
    genes = [f"g{i:04d}" for i in range(200)]
    econ, work = economy("priority", order=genes)
    factors = econ.scale(work)
    assert len({round(v, 12) for v in factors.values()}) > 1


def test_a_global_factor_cannot_change_a_fitted_absolute_error() -> None:
    """The algebra the primary died of, in the smallest form that shows it."""
    from scripts.abundance_gate import fitted_shift, mae  # noqa: PLC0415 - script, not a package

    y = [3.0, 1.0, 4.0, 1.5, 5.0, 9.0, 2.0, 6.0]
    x = [2.0, 1.2, 3.1, 2.0, 4.4, 7.0, 2.9, 5.5]
    base = mae(y, [v + fitted_shift(y, x) for v in x])
    for log_f in (-3.0, -0.5, 0.7, 2.0):
        moved = [v + log_f for v in x]
        got = mae(y, [v + fitted_shift(y, moved) for v in moved])
        assert abs(got - base) <= DERIVATION_TOLERANCE_LOG10


def test_a_per_demander_share_would_move_it_which_is_why_the_instrument_is_named() -> None:
    """The amendment names capacity / (capacity + demand_i) as the form that would work. Check it."""
    import math  # noqa: PLC0415

    from scripts.abundance_gate import fitted_shift, mae  # noqa: PLC0415

    capacity = 1.0e4
    demand = [10.0**i for i in (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5)]
    y = [math.log10(d) * 0.5 + 2.0 for d in demand]  # a compressed truth, slope 0.5
    x = [math.log10(d) for d in demand]
    base = mae(y, [v + fitted_shift(y, x) for v in x])
    saturable = [math.log10(d * capacity / (capacity + d)) for d in demand]
    got = mae(y, [v + fitted_shift(y, saturable) for v in saturable])
    assert got < base, "a per-demander share must be able to beat the one-to-one arm"


def test_the_registration_says_the_primary_cannot_move_and_names_the_instrument() -> None:
    assert "CANNOT BE MET" in AMENDMENT["consequence_for_the_primary"]
    assert "per-demander" in AMENDMENT["instrument_that_would_make_the_primary_able_to_move"]
    assert "slope" in AMENDMENT["new_registered_quantity"]
    # §5.3's consequence was agreed before any of this, and the amendment does not soften it
    assert "stays optional" in PRE_REGISTRATION["outcomes"]["fail"]
    assert "stays optional" in AMENDMENT["verdict_rule"]


@pytest.mark.skipif(not RESULT.exists(), reason="gate (b) has not been run in this checkout")
def test_the_recorded_run_fails_the_gate_and_the_derivation_holds() -> None:
    got = json.loads(RESULT.read_text())
    assert got["verdict"]["passed"] is False
    assert got["derivation_check"]["holds"] is True
    assert got["primary"]["improvement"] == 0.0
    for label, row in got["runtime_factors_over_the_joined_genes"].items():
        if isinstance(row, dict):
            assert row["distinct_factors"] == 1, label
    short = got["runtime_factors_over_the_joined_genes"]["proportional_with_the_pools_short"]
    assert short["constrained"] is True, "one arm must actually run the pools short"
