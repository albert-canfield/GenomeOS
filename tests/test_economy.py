# SPDX-License-Identifier: AGPL-3.0-or-later
"""Stage 2's arithmetic: demand, capacity, and who gets cut (BIOLANG-v0.4-ECONOMY.md §5).

The burden gate's falsifier is the test that matters here and it is the last one: the reduction must
disappear when the pool is not the constraint. A model that reduces an unrelated gene's output
whatever the capacity is not reproducing a burden, it is reproducing a bug.
"""

from __future__ import annotations

import pytest

from genomeos.ir import UNKNOWN, Allocation, Cost, Pool
from genomeos.runtime.economy import Economy, EconomyError


def pool(pid: str, size: float | object = 1_000.0, **kw) -> Pool:
    p = Pool(id=pid, kind="pool", **kw)
    p.size = size
    return p


def econ(size: float | object = 1_000.0, policy: str = "proportional", **kw) -> Economy:
    return Economy(
        pools={"R": pool("R", size)},
        costs={"A": [Cost("R", 1.0, "chain")], "B": [Cost("R", 1.0, "chain")]},
        allocations={"R": Allocation(id="a", kind="allocation", pool="R", policy=policy, **kw)},
    )


def test_within_capacity_nobody_is_scaled() -> None:
    assert econ().scale({"A": 100.0, "B": 100.0}) == {"A": 1.0, "B": 1.0}


def test_over_capacity_proportional_scales_everyone_by_the_same_factor() -> None:
    """One factor for all is exactly why the pool layer cannot move a rank correlation."""
    got = econ().scale({"A": 1_500.0, "B": 500.0})

    assert got["A"] == got["B"] == pytest.approx(0.5)


def test_priority_serves_its_order_and_starves_the_rest() -> None:
    e = econ(policy="priority", order=["B"])

    got = e.scale({"A": 800.0, "B": 800.0})

    assert got["B"] == 1.0, "the named order is served first"
    assert got["A"] == pytest.approx(0.25), "what is left over is 200 of the 800 it wanted"


def test_competitive_serves_nobody_in_full() -> None:
    got = econ(policy="competitive").scale({"A": 1_500.0, "B": 500.0})

    assert 0.0 < got["A"] < 1.0 and got["A"] == got["B"]


def test_optimise_is_refused_at_run_time_because_it_is_a_modelling_device() -> None:
    with pytest.raises(EconomyError, match="does not solve an objective"):
        econ(policy="optimise", objective="growth").scale({"A": 1_500.0, "B": 500.0})


def test_an_entity_is_limited_by_its_tightest_pool() -> None:
    """A protein with no ribosome is not made whatever the ATP: the work needs all of its inputs."""
    e = Economy(
        pools={"R": pool("R", 100.0), "ATP": pool("ATP", 1e9)},
        costs={"A": [Cost("R", 1.0, "chain"), Cost("ATP", 1.0, "chain")]},
    )

    assert e.scale({"A": 1_000.0})["A"] == pytest.approx(0.1)


def test_a_cost_charging_an_undeclared_pool_is_refused() -> None:
    """Charging to nothing would make the burden vanish silently."""
    with pytest.raises(EconomyError, match="which no pool declares"):
        Economy(pools={}, costs={"A": [Cost("R", 1.0, "chain")]}).check()


def test_a_cost_on_a_pool_of_unknown_size_is_refused() -> None:
    """No ceiling means no burden can ever appear, so "no burden found" would mean nothing."""
    with pytest.raises(EconomyError, match="no ceiling"):
        Economy(pools={"R": pool("R", UNKNOWN)}, costs={"A": [Cost("R", 1.0, "chain")]}).check()


def test_a_per_length_cost_with_no_length_is_refused() -> None:
    """Charging zero would make a long molecule as cheap as a short one."""
    with pytest.raises(EconomyError, match="no length is known"):
        Economy(pools={"R": pool("R")}, costs={"A": [Cost("R", 4.0, "residue")]}).check()


def test_a_per_length_cost_multiplies_by_the_length() -> None:
    e = Economy(
        pools={"R": pool("R", 1e9)},
        costs={"A": [Cost("R", 4.0, "residue")]},
        lengths={"A": 300.0},
    )

    assert e.demand_by_pool({"A": 10.0}) == {"R": 4.0 * 300.0 * 10.0}


def test_expressing_one_protein_harder_reduces_an_unrelated_one_sharing_the_pool() -> None:
    """Gate (a)'s documented direction (Ceroni 2015, Frei 2020), as arithmetic."""
    e = econ()
    housekeeper = [e.scale({"A": burden, "B": 10.0})["B"] for burden in (10.0, 1_000.0, 10_000.0)]

    assert housekeeper[0] == 1.0
    assert housekeeper[0] > housekeeper[1] > housekeeper[2], "it must fall, and keep falling"


def test_the_reduction_disappears_when_the_pool_is_not_the_constraint() -> None:
    """The falsifier, and the reason the gate can fail.

    §5.3: "Falsified if no reduction appears, or if it appears with the pool switched off." A model
    that cuts an unrelated gene whatever the capacity is not reproducing a burden; it has a bug. So
    the same demands against a pool large enough to serve them must leave every factor at 1.
    """
    huge = econ(size=1e12)

    assert huge.scale({"A": 10_000.0, "B": 10.0}) == {"A": 1.0, "B": 1.0}
    assert huge.report({"A": 10_000.0, "B": 10.0})["any_constrained"] is False
