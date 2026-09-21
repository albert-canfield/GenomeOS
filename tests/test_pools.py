# SPDX-License-Identifier: AGPL-3.0-or-later
"""Stage 2's language surface: pool, cost and allocation (BIOLANG-v0.4-ECONOMY.md §5).

Stage 2 was held for two years of project-days on open decision 2, which was resolved on 2026-09-19.
These tests cover what the language now accepts and — more to the point — what it refuses, because a
pool whose capacity nobody measured cannot support a burden claim, and a cost without a unit is
ambiguous between two different questions about the same gene.
"""

from __future__ import annotations

import pytest

from genomeos.ir import UNKNOWN, Allocation, Cost, Pool
from genomeos.lang import parse
from genomeos.lang.parser import BioLangError

HEAD = "module t\n\ncompartment Cytosol { volume: 1.0; translation: yes }\n"


def one(body: str):
    return parse(HEAD + body)


def test_a_pool_carries_its_capacity_location_and_return() -> None:
    m = one("pool Ribosomes { location: Cytosol; size: 5e6; regenerates: 1e5 /h; returns_as: Ribosomes }")
    p = m.entities["Ribosomes"]

    assert isinstance(p, Pool)
    assert (p.location, p.size, p.regeneration, p.returns_as) == ("Cytosol", 5e6, 1e5, "Ribosomes")


def test_a_pool_may_be_refilled_by_processes_rather_than_at_a_rate() -> None:
    """ATP is the reason this shape exists: it is not restored at a constant rate."""
    m = one(
        "pool ATP { location: Cytosol; size: 3e9; regenerates: from Glycolysis, OXPHOS; returns_as: ADP }"
    )
    p = m.entities["ATP"]

    assert p.regenerates_from == ["Glycolysis", "OXPHOS"]
    assert p.regeneration is UNKNOWN, "a process list is not a rate and must not be read as one"


def test_an_unmeasured_capacity_is_unknown_rather_than_a_guess() -> None:
    """`unknown` is a statement: no burden claim can rest on an invented ceiling."""
    m = one("pool Ribosomes { location: Cytosol; size: unknown }")

    assert m.entities["Ribosomes"].size is UNKNOWN


def test_a_capacity_of_zero_or_less_is_refused() -> None:
    with pytest.raises(BioLangError, match="capacity in molecules"):
        one("pool Ribosomes { location: Cytosol; size: 0 }")


def test_a_cost_states_its_pool_amount_and_unit() -> None:
    m = one("gene G { location: Cytosol; cost: RNAPII 1 per transcript ; cost: ATP 2 per nucleotide }")

    assert m.entities["G"].costs == [Cost("RNAPII", 1.0, "transcript"), Cost("ATP", 2.0, "nucleotide")]


def test_a_cost_without_a_unit_is_refused_because_it_is_ambiguous() -> None:
    """ "ATP 4" could be per residue or per chain, and those are different claims."""
    with pytest.raises(BioLangError, match="per <unit>"):
        one("protein P { location: Cytosol; cost: ATP 4 }")


def test_a_cost_per_something_outside_the_closed_set_is_refused() -> None:
    with pytest.raises(BioLangError, match="expected one of"):
        one("protein P { location: Cytosol; cost: ATP 4 per fortnight }")


def test_an_allocation_names_the_pool_it_divides() -> None:
    m = one("pool R { location: Cytosol; size: 1e6 }\nallocation share { pool: R; policy: proportional }")
    a = m.entities["share"]

    assert isinstance(a, Allocation)
    assert (a.pool, a.policy) == ("R", "proportional")


def test_an_allocation_with_no_pool_is_refused() -> None:
    with pytest.raises(BioLangError, match="must name the pool"):
        one("allocation share { policy: proportional }")


def test_priority_needs_an_order_and_optimise_needs_an_objective() -> None:
    """A policy whose parameter is missing is a claim nobody stated."""
    with pytest.raises(BioLangError, match="needs an order"):
        one("allocation s { pool: R; policy: priority }")
    with pytest.raises(BioLangError, match="declare its objective"):
        one("allocation s { pool: R; policy: optimise }")


def test_an_unknown_policy_is_refused_rather_than_defaulted() -> None:
    with pytest.raises(BioLangError, match="expected one of"):
        one("allocation s { pool: R; policy: whatever_seems_best }")


def test_the_proportional_trap_is_recorded_where_the_policy_is_defined() -> None:
    """The abundance gate's pre-stated trap must live with the construct, not only in the doc.

    A shared pool under `proportional` rescales every gene by the same factor, so it cannot change a
    rank correlation. Anything crediting allocation with a rank improvement is crediting something
    else, which is why the gate reads absolute error in log copy numbers.
    """
    assert "cannot change a rank correlation" in (Allocation.__doc__ or "")
