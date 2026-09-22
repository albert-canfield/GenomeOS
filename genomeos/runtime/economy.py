# SPDX-License-Identifier: Apache-2.0
# Part of the BioLang engine (language, IR, VM, standard library); see LICENSING.md.
"""Stage 2's arithmetic: what the work costs, what the pools can supply, and who gets cut.

Stage 1 could express a cost and could not charge it. This module is what makes a cost mean
something: it turns each entity's declared draws into a demand per pool, compares that with the
pool's capacity, and returns the factor by which each entity's rate must be scaled when demand
exceeds supply. The scaling is where a burden comes from — expressing one protein hard leaves less
ribosome for every other protein sharing the pool, which is the documented effect stage 2's gate (a)
has to reproduce (Ceroni et al. 2015, Nat Methods 12:415; Frei et al. 2020, Nat Commun 11:4641).

**Two refusals rather than two guesses**, and both are the point rather than defensiveness:

- a cost naming a pool the program never declared is a typo or a missing declaration, and charging it
  to nothing would make the burden vanish silently;
- a cost drawing on a pool whose `size` is UNKNOWN has no ceiling, so it can never be exceeded and no
  burden can ever appear. Running such a program would return "no burden found" for a reason that has
  nothing to do with biology. `Pool`'s own docstring promises this refusal, so it is enforced here.

**What a unit multiplies, stated once.** `per transcript` and `per chain` are charged once per
molecule made. `per nucleotide` and `per residue` are charged per molecule times its length, so they
need a length; a program that asks for one without stating it is refused rather than charged zero,
which would quietly make a long protein as cheap as a short one. `per event` and `per division` are
charged once when the event fires.

Nothing here decides a policy. The policy is the program's claim (§5.3), and `Allocation` carries the
warning that matters for gate (b): under `proportional` every demander is scaled by the same factor,
so the whole pool layer cannot change a rank correlation and a rank improvement credited to it would
be credited to something else.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from genomeos.ir import UNKNOWN, Allocation, Cost, Module, Pool

#: units charged once per molecule of work, against those charged per base or residue of it
PER_MOLECULE = ("transcript", "chain", "event", "division")
PER_LENGTH = ("nucleotide", "residue")


class EconomyError(ValueError):
    """A program whose costs cannot be charged as written."""


@dataclass
class Draw:
    """One entity's demand on one pool, with the work that produced it kept beside the number."""

    entity: str
    pool: str
    per: str
    amount: float
    units_of_work: float
    demand: float


@dataclass
class Economy:
    """The pools a module declares, the costs drawn on them, and the policy that divides each.

    `lengths` maps an entity id to the length its per-base or per-residue costs multiply — transcript
    nucleotides for a gene, residues for a protein. It is supplied by the caller because the runtime
    knows it and the IR does not always: a gene may carry a locus, a sequence, both or neither.
    """

    pools: dict[str, Pool]
    costs: dict[str, list[Cost]]
    allocations: dict[str, Allocation] = field(default_factory=dict)
    lengths: dict[str, float] = field(default_factory=dict)

    @classmethod
    def from_module(cls, module: Module, lengths: dict[str, float] | None = None) -> Economy:
        pools = {e.id: e for e in module.entities.values() if isinstance(e, Pool)}
        allocations = {}
        for e in module.entities.values():
            if isinstance(e, Allocation):
                allocations[e.pool] = e
        costs: dict[str, list[Cost]] = {}
        for e in module.entities.values():
            got = list(getattr(e, "costs", None) or [])
            if got:
                costs[e.id] = got
        for ev in getattr(module, "events", None) or []:  # events are a list, not a mapping
            got = list(getattr(ev, "costs", None) or [])
            if got:
                costs[ev.id] = got
        econ = cls(pools=pools, costs=costs, allocations=allocations, lengths=dict(lengths or {}))
        econ.check()
        return econ

    # ---- refusals --------------------------------------------------------------------------
    def check(self) -> None:
        """Refuse a program whose costs cannot be charged, naming every reason at once."""
        problems: list[str] = []
        for entity, costs in sorted(self.costs.items()):
            for c in costs:
                if c.pool not in self.pools:
                    problems.append(f"{entity} charges {c.pool}, which no pool declares")
                    continue
                if self.pools[c.pool].size is UNKNOWN:
                    problems.append(
                        f"{entity} charges {c.pool}, whose size is unknown: a pool with no ceiling"
                        " cannot be exceeded, so no burden could ever appear"
                    )
                if c.per in PER_LENGTH and entity not in self.lengths:
                    problems.append(
                        f"{entity} charges {c.pool} per {c.per} and no length is known for it:"
                        " charging zero would make a long molecule as cheap as a short one"
                    )
        if problems:
            raise EconomyError("; ".join(problems))

    # ---- demand and supply -----------------------------------------------------------------
    def draws(self, work: dict[str, float]) -> list[Draw]:
        """Every demand implied by `work`, which maps an entity id to molecules made per unit time."""
        out: list[Draw] = []
        for entity, rate in work.items():
            for c in self.costs.get(entity, []):
                units = rate * (self.lengths[entity] if c.per in PER_LENGTH else 1.0)
                out.append(
                    Draw(
                        entity=entity,
                        pool=c.pool,
                        per=c.per,
                        amount=c.amount,
                        units_of_work=units,
                        demand=c.amount * units,
                    )
                )
        return out

    def demand_by_pool(self, work: dict[str, float]) -> dict[str, float]:
        total: dict[str, float] = {}
        for d in self.draws(work):
            total[d.pool] = total.get(d.pool, 0.0) + d.demand
        return total

    def scale(self, work: dict[str, float]) -> dict[str, float]:
        """The factor each entity's rate must be multiplied by, given what the pools can supply.

        An entity drawing on several pools is limited by the tightest of them, because the work needs
        all of its inputs: a protein with no ribosome does not get made whatever the ATP.
        """
        demand = self.demand_by_pool(work)
        share: dict[str, dict[str, float]] = {}
        for pool, wanted in demand.items():
            capacity = float(self.pools[pool].size)
            policy = (self.allocations.get(pool) or Allocation(id="", kind="allocation")).policy
            share[pool] = self._share(pool, policy, capacity, wanted, work)
        factors = dict.fromkeys(work, 1.0)
        for per_entity in share.values():
            for entity, f in per_entity.items():
                factors[entity] = min(factors[entity], f)
        return factors

    def _share(
        self, pool: str, policy: str, capacity: float, wanted: float, work: dict[str, float]
    ) -> dict[str, float]:
        demanders = {d.entity: d.demand for d in self.draws(work) if d.pool == pool and d.demand > 0}
        if wanted <= capacity or not demanders:
            return dict.fromkeys(demanders, 1.0)
        if policy == "proportional":
            # one factor for everybody: this is exactly why the pool layer cannot move a rank
            # correlation, which gate (b) is told in advance
            return dict.fromkeys(demanders, capacity / wanted)
        if policy == "priority":
            order = (self.allocations[pool].order or []) + sorted(demanders)
            left, out = capacity, {}
            for entity in dict.fromkeys(order):
                if entity not in demanders:
                    continue
                want = demanders[entity]
                out[entity] = 1.0 if want <= left else max(0.0, left / want)
                left = max(0.0, left - want)
            return out
        if policy == "competitive":
            # saturable, as a transport is: each demander gets a share that falls off as the pool is
            # oversubscribed, and nobody is served in full
            return {e: capacity / (capacity + wanted) for e in demanders}
        if policy == "optimise":
            raise EconomyError(
                f"pool {pool} is allocated by `optimise`, a declared modelling device: the runtime"
                " does not solve an objective, so the program must state a policy that can be run"
            )
        raise EconomyError(f"pool {pool}: unknown allocation policy {policy!r}")

    def report(self, work: dict[str, float]) -> dict[str, Any]:
        """Demand, capacity and the resulting factor per pool — what a burden claim has to show."""
        demand = self.demand_by_pool(work)
        factors = self.scale(work)
        return {
            "pools": {
                pool: {
                    "capacity": float(p.size),
                    "demand": round(demand.get(pool, 0.0), 6),
                    "oversubscribed": demand.get(pool, 0.0) > float(p.size),
                    "policy": (self.allocations.get(pool) or Allocation(id="", kind="allocation")).policy,
                }
                for pool, p in sorted(self.pools.items())
            },
            "factor": {k: round(v, 6) for k, v in sorted(factors.items())},
            "any_constrained": any(v < 1.0 for v in factors.values()),
        }
