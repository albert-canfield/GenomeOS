# SPDX-License-Identifier: Apache-2.0
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


# ---- `saturable`: the per-demander policy (2026-10-10, lane-saturable) ------------------------------


def sat(size: float = 1e12, k: float | object = 1_000.0) -> Economy:
    return econ(size=size, policy="saturable", attrs={"half_saturation": k})


def test_saturable_gives_each_demander_its_own_factor_falling_in_its_own_demand() -> None:
    """The gene index the other three policies lack: K / (K + d) per demander."""
    got = sat().scale({"A": 1_500.0, "B": 500.0})

    assert got["A"] == pytest.approx(1_000.0 / 2_500.0)
    assert got["B"] == pytest.approx(1_000.0 / 1_500.0)
    assert got["A"] < got["B"], "the larger demand is cut more"


def test_saturable_is_not_a_burden_while_the_pool_is_not_short() -> None:
    """Gate (a)'s falsifier, kept: with the pool not the constraint, another demander's demand never
    enters B's factor. B is cut by its own saturation only, and that cut does not move."""
    e = sat(size=1e12)
    b = [e.scale({"A": burden, "B": 10.0})["B"] for burden in (10.0, 1_000.0, 1e6)]

    assert b[0] == b[1] == b[2] == pytest.approx(1_000.0 / 1_010.0)


def test_saturable_enforces_the_pool_on_what_is_served_and_then_burdens() -> None:
    """When the served total exceeds the capacity, every served amount is cut by one factor, so the
    pool is used exactly; raising A then lowers B, which is the documented burden direction."""
    e = sat(size=500.0)
    work = {"A": 1_500.0, "B": 500.0}
    got = e.scale(work)
    served = sum(work[x] * got[x] for x in work)

    assert served == pytest.approx(500.0)
    assert got["A"] < got["B"]
    b = [e.scale({"A": burden, "B": 500.0})["B"] for burden in (500.0, 5_000.0, 50_000.0)]
    assert b[0] > b[1] > b[2]


@pytest.mark.parametrize("k", [UNKNOWN, 0.0, -5.0, float("nan"), float("inf")])
def test_saturable_without_a_usable_k_is_refused_not_defaulted(k: object) -> None:
    e = sat(k=k)
    with pytest.raises(EconomyError, match="half_saturation"):
        e.check()
    with pytest.raises(EconomyError, match="does not choose K"):
        e.scale({"A": 1.0, "B": 1.0})


def test_the_three_existing_policies_are_unchanged_by_saturable() -> None:
    """Pinned values, computed by the policies as they stood at 5c620d3, before saturable existed."""
    work = {"A": 1_500.0, "B": 500.0}
    assert econ().scale(work) == {"A": 0.5, "B": 0.5}
    assert econ(policy="competitive").scale(work) == {"A": 1_000.0 / 3_000.0, "B": 1_000.0 / 3_000.0}
    assert econ(policy="priority", order=["B"]).scale(work) == {"A": 500.0 / 1_500.0, "B": 1.0}
    assert econ(size=1e12, policy="competitive").scale(work) == {"A": 1.0, "B": 1.0}


def test_a_program_declares_saturable_with_its_k_and_the_parser_refuses_the_rest() -> None:
    from genomeos.lang import parse
    from genomeos.lang.parser import BioLangError

    head = "module t\n\ncompartment Cytosol { volume: 1.0; translation: yes }\n"
    head += "pool Ribosomes { location: Cytosol; size: 5e6 }\n"
    m = parse(head + "allocation share { pool: Ribosomes; policy: saturable; half_saturation: 40 }")
    al = m.entities["share"]
    assert (al.policy, al.attrs["half_saturation"]) == ("saturable", 40.0)
    Economy.from_module(m)  # runs: K is stated

    for body, why in (
        ("allocation s { pool: Ribosomes; policy: saturable }", "needs a finite half_saturation"),
        ("allocation s { pool: Ribosomes; policy: saturable; half_saturation: 0 }", "above 0"),
        ("allocation s { pool: Ribosomes; policy: saturable; half_saturation: unknown }", "expects a number"),
        ("allocation s { pool: Ribosomes; policy: proportional; half_saturation: 40 }", "read only by"),
        ("regime r { allocation: saturable }", "regime allocation must be one of"),
    ):
        with pytest.raises(BioLangError, match=why):
            parse(head + body)
