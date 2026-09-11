# SPDX-License-Identifier: Apache-2.0
"""Body runtime: grows an organism from one cell by discrete events.

Every cell (or population: a cell with `count` > 1) reads its context, the
program's decisions say what it does, timers say when, and signals from
other cells change what it reads. The loop per cell:

    read     context = name, lineage, generation, cell type, stage, factors, environment
    decide   die | differentiate | quiesce | divide     (first matching decision per action)
    wait     the matching timer (tempo-scaled, lengthening per generation)
    write    daughters inherit the factors; `asymmetric` keeps a factor in one daughter

A cell for which no decision or no timer applies stops and is reported as
UNKNOWN: the program did not say. Nothing is invented.
"""

from __future__ import annotations

import bisect
import heapq
import math
import random
import re
from collections import Counter
from dataclasses import dataclass, field

from genomeos.ir import Decision, EvidenceKind, Module, Timer, matches, to_minutes

from .spatial import Field2D
from .uncertainty import UncertaintyReport

SUFFIXES = (("a", "p"), ("l", "r"))


@dataclass(slots=True)
class Cell:
    name: str
    lineage: str
    generation: int
    born: float
    cell_type: str
    factors: dict[str, str] = field(default_factory=dict)
    parent: str = ""
    children: list[str] = field(default_factory=list)
    count: float = 1.0  # > 1: a population
    population: bool = False  # counted node: division grows the count in place, so it does not end the node
    divides_at: float | None = None
    dies_at: float | None = None
    terminal_name: str = ""
    quiescent: bool = False
    fired: list[str] = field(default_factory=list)
    unknown: str = ""  # why the program left this cell without a next step
    cull_at: float | None = None  # next fractional death of a population, if a chain is running
    flow_at: dict[str, float] = field(default_factory=dict)  # recurring differentiation flows by decision id
    x: int | None = None  # grid position when the organism has a space
    y: int | None = None
    move_at: float | None = None  # next step of a recurring migration
    divide_seq: int = -1  # event token of the pending division: a rescheduled division runs exactly once
    measured: set[str] = field(default_factory=set)  # factors set by express decisions (the reader)
    stated: set[str] = field(default_factory=set)  # factors stated by mechanism: maternal load, asymmetry

    @property
    def end(self) -> float:
        if self.population:
            return self.dies_at if self.dies_at is not None else math.inf
        ends = [x for x in (self.divides_at, self.dies_at) if x is not None]
        return min(ends) if ends else math.inf

    @property
    def terminal(self) -> bool:
        return not self.children and self.dies_at is None


class Body:
    def __init__(
        self,
        module: Module,
        seed: int | None = None,
        max_cells: int = 200_000,
        means: bool = False,
        knockouts: set[str] | frozenset[str] = frozenset(),
        adds: set[str] | frozenset[str] = frozenset(),
        environment: dict[str, str] | None = None,
    ):
        """`seed` overrides the organism's declared seed; `means=True` runs every timer at its mean;
        `knockouts` are factors never present and signals (by id, ligand or receptor) never sent;
        `adds` are extra factors in the zygote; `environment` overrides the organism's."""
        if module.organism is None:
            raise ValueError(f"module {module.name!r} declares no organism block")
        self.module = module
        self.organism = module.organism
        self.knockouts = set(knockouts)
        self.adds = set(adds)
        self.environment = {**module.organism.environment, **(environment or {})}
        # space: fields on the organism's grid, stepped between events; sites held by positioned cells
        self.fields: dict[str, Field2D] = {}
        self.sources: list[tuple[str, int, int, float]] = []
        self.occupied: dict[tuple[int, int], str] = {}
        self.field_time = 0.0
        self.blocked_divisions = 0  # divisions that found no free site (contact inhibition)
        if module.organism.width > 0 and module.organism.height > 0:
            for f in module.fields:
                self.fields[f.name] = Field2D(
                    f.name, module.organism.width, module.organism.height, f.diffusion, f.decay
                )
                self.sources.extend((f.name, x, y, rate) for x, y, rate in f.sources)
        self._gradients = [sg for sg in module.signals() if sg.mode == "gradient" and sg.sets]
        if seed is None:
            seed = self.organism.seed
        self.rng = random.Random(seed) if seed is not None and not means else None
        self.max_cells = max_cells
        self.cells: dict[str, Cell] = {}
        self.time = 0.0
        self.fired: Counter[str] = Counter()
        self.culled: float = 0.0  # population cells removed by fractional deaths
        self.history: list[tuple[float, float]] = []  # (time, total count) whenever a population changes
        self.unknown: Counter[str] = Counter()
        self.timers_used: Counter[str] = Counter()
        self._queue: list[tuple[float, int, str, str]] = []
        self._seq = 0
        self._decisions_by_action: dict[str, list[Decision]] = {}
        for d in module.decisions:
            self._decisions_by_action.setdefault(d.action, []).append(d)
        self._by_cell: dict[str, list[Decision]] = {}
        for d in module.decisions:
            if "cell" in d.when and "|" not in d.when["cell"] and d.when["cell"] != "any":
                self._by_cell.setdefault(d.when["cell"], []).append(d)
        named = {x.id for v in self._by_cell.values() for x in v}
        self._general = [d for d in module.decisions if d.id not in named]
        self._divide_ids = {d.id for d in module.decisions if d.action == "divide"}
        self._expressed_names = {f for d in module.decisions if d.action == "express" for f in d.sets}
        self._bootstrap()

    # ---- setup ---------------------------------------------------------

    def _bootstrap(self) -> None:
        o = self.organism
        for st in self.module.stages:
            self._push(to_minutes(st.start, st.unit), "", "stage")
        factors = {f: "present" for f in list(o.factors) + sorted(self.adds) if f not in self.knockouts}
        root = Cell(o.root, "", 0, 0.0, o.cell_type, factors)
        root.stated = set(factors)
        root.population = o.resolution == "populations"
        if self.spatial:
            root.x, root.y = o.origin
            self.occupied[o.origin] = root.name
            if self._gradients:
                self._push(o.sense, "", "sense")
        self._add(root, None)

    @property
    def spatial(self) -> bool:
        return self.organism.width > 0 and self.organism.height > 0

    def set_factor(self, name: str, value: str | None) -> None:
        """An organism-wide factor set from outside (a coupled process, an environment change): every
        population and resting cell reads it and decides again; cells born later carry it in their context."""
        current = self.environment.get(name)
        if value is None:
            self.environment.pop(name, None)
        else:
            self.environment[name] = value
        if current != value:
            self._stage_change()

    def _stage_change(self) -> None:
        """Populations and resting cells read the new stage and decide again."""
        for c in self.alive_at(self.time):
            if self.population(c) or c.quiescent:
                self._resolve(c, born=False)

    def _push(self, t: float, name: str, kind: str) -> int:
        self._seq += 1
        heapq.heappush(self._queue, (t, self._seq, name, kind))
        return self._seq

    def population(self, c: Cell) -> bool:
        return c.population

    # ---- context and decisions ------------------------------------------

    def context(self, c: Cell, t: float | None = None) -> dict[str, str]:
        t = self.time if t is None else t
        ctx = {
            "cell": c.name,
            "lineage": c.lineage,
            "generation": str(c.generation),
            "cell_type": c.cell_type,
            "stage": self.module.stage_at(t),
            "count": str(c.count),
        }
        ctx.update(self.environment)
        ctx.update(c.factors)
        if c.x is not None:
            ctx["x"], ctx["y"] = str(c.x), str(c.y)
        return ctx

    def _candidates(self, c: Cell) -> list[Decision]:
        return self._by_cell.get(c.name, []) + self._general

    def _first(self, c: Cell, action: str, ctx: dict[str, str], unfired: bool = False) -> Decision | None:
        for d in self._candidates(c):
            if d.action == action and d.applies(ctx) and not (unfired and d.id in c.fired):
                return d
        return None

    def _timer_for(self, c: Cell, d: Decision, ctx: dict[str, str]) -> Timer | None:
        if d.timer:
            return self.module.timer(d.timer)
        return next((t for t in self.module.timers if t.applies(ctx)), None)

    def _duration(self, timer: Timer, c: Cell) -> float:
        dur = timer.minutes() * (timer.lengthening ** max(0, c.generation - 1)) * self.organism.tempo
        if self.rng is not None and timer.sd > 0:
            dur = max(0.1, self.rng.gauss(dur, timer.sd * self.organism.tempo))
        return dur

    def _resolve(self, c: Cell, born: bool = True) -> None:
        """Decide what a cell does next; called at birth and after a signal changes its factors."""
        ctx = self.context(c)
        if c.dies_at is None and (d := self._first(c, "die", ctx)):
            if self.population(c) and d.fraction < 1.0:
                if (
                    c.cull_at is None
                ):  # start the recurring loss; _cull keeps it going while a decision applies
                    c.cull_at = self.time + (d.after or 0.0) * self.organism.tempo
                    self._push(c.cull_at, c.name, "cull")
            else:
                c.dies_at = c.born + (d.after or 0.0) * self.organism.tempo
                c.fired.append(d.id)
                self.fired[d.id] += 1
                self._push(c.dies_at, c.name, "die")
        self._express(c, ctx)
        ctx = self.context(c)
        for _ in range(32):  # differentiation chains and sequential population splits
            d = self._first(c, "differentiate", ctx, unfired=True)
            if d is None:
                break
            if d.fraction < 1.0 and self.population(c) and d.after is not None:
                # a recurring flow: every `after`, this share moves into the target pool
                if d.id not in c.flow_at:
                    c.flow_at[d.id] = self.time + d.after * self.organism.tempo
                    self._push(c.flow_at[d.id], c.name, f"flow:{d.id}")
            elif d.fraction < 1.0 and self.population(c):
                self._split(c, d)
            else:
                c.cell_type = d.to
                if d.name:
                    c.terminal_name = d.name
            c.fired.append(d.id)
            self.fired[d.id] += 1
            ctx = self.context(c)
        d = self._first(c, "quiesce", ctx)
        c.quiescent = d is not None  # re-evaluated at every decision point
        if d is not None and d.id not in c.fired:
            c.fired.append(d.id)
            self.fired[d.id] += 1
        if c.quiescent and c.divides_at is not None and self.population(c):
            c.divides_at = None
        if (d := self._first(c, "migrate", ctx)) and d.id not in c.fired:
            c.fired.append(d.id)
            self.fired[d.id] += 1
            if c.x is not None:
                self._move(c, d)
                if d.after is not None and c.move_at is None:
                    c.move_at = self.time + d.after * self.organism.tempo
                    self._push(c.move_at, c.name, f"move:{d.id}")
        swapped = False
        if c.divides_at is not None and not born and not c.quiescent:
            # a signal changed what the cell reads: an earlier-precedence division may now apply
            d = self._first(c, "divide", ctx)
            current = next((x for x in c.fired if x in self._divide_ids), None)
            if d is not None and current is not None and d.id != current:
                c.fired.remove(current)
                self.fired[current] -= 1
                if d.after is not None and self.population(c):
                    # same cadence: keep the pending step's time, only the decision changes (a re-timed
                    # step would shift the phase between growth and outflow, which the balance depends on)
                    c.fired.append(d.id)
                    self.fired[d.id] += 1
                else:
                    c.divides_at = None
                    swapped = True  # a cell whose division decision changed is rescheduled from now
        if c.divides_at is None and not c.quiescent and (born or swapped or self.population(c)):
            d = self._first(c, "divide", ctx)
            if d is None:
                if not any(x for x in c.fired if x.startswith(("fate_", "die_"))) and c.cell_type == "":
                    c.unknown = "no decision"
                    self.unknown["no decision"] += 1
                return
            timer = self._timer_for(c, d, ctx)
            if timer is None and d.after is None:
                c.unknown = "no timer"
                self.unknown["no timer"] += 1
                return
            if d.after is not None:
                wait = d.after * self.organism.tempo
            else:
                wait = self._duration(timer, c)
                self.timers_used[timer.name] += 1
            # cells wait from birth (or from now when they decide again later, e.g. after a blocked division);
            # populations wait from now (they decide again after every step)
            c.divides_at = (self.time if self.population(c) else max(c.born, self.time)) + wait
            c.fired.append(d.id)
            self.fired[d.id] += 1
            c.divide_seq = self._push(c.divides_at, c.name, "divide")

    def _record(self) -> None:
        total = sum(x.count for x in self.cells.values() if x.born <= self.time < x.end)
        if self.history and self.history[-1][0] == self.time:
            self.history[-1] = (self.time, total)
        else:
            self.history.append((self.time, total))

    def _express(self, c: Cell, ctx: dict[str, str]) -> None:
        """Every applying `express` decision sets its factors; a cell the reader covers carries exactly the
        measured set, so measured factors inherited from the parent but not listed here are dropped."""
        if not self._expressed_names:
            return
        applying = [d for d in self._candidates(c) if d.action == "express" and d.applies(ctx)]
        if not applying:
            return
        wanted: set[str] = set()
        for d in applying:
            wanted.update(f for f in d.sets if f not in self.knockouts)
            if d.id not in c.fired:
                c.fired.append(d.id)
                self.fired[d.id] += 1
        # only factors that came from the reader are replaced; maternal factors, asymmetric inheritance and
        # signals are stated mechanism and keep their say even when a reporter does not see the protein
        for f in c.measured - wanted - c.stated:
            if c.factors.get(f) == "present":
                del c.factors[f]
        for f in wanted:
            c.factors.setdefault(f, "present")
        c.measured = set(wanted) - c.stated

    def _split(self, c: Cell, d: Decision) -> None:
        """A fraction of a population differentiates into the target pool (created on first use)."""
        part = c.count * d.fraction
        c.count -= part
        self._pour(c, d.to, part)
        self._record()

    def _pour(self, c: Cell, cell_type: str, amount: float) -> None:
        name = f"{c.name}>{cell_type}"
        pool = self.cells.get(name)
        if pool is not None:
            pool.count += amount
            if pool.born <= self.time < pool.end:
                self._resolve(pool, born=False)  # the receiving population decides again too
            return
        child = Cell(
            name, c.lineage, c.generation, self.time, cell_type, dict(c.factors), parent=c.name, count=amount
        )
        child.population = True
        child.measured = set(c.measured)
        child.stated = set(c.stated)
        self._add(child, c)

    def _flow(self, c: Cell, did: str) -> None:
        """One step of a recurring flow; the chain continues while the decision applies."""
        c.flow_at.pop(did, None)
        d = next((x for x in self.module.decisions if x.id == did), None)
        alive = c.born <= self.time < c.end
        if d is None or not alive or not d.applies(self.context(c)):
            if did in c.fired:
                c.fired.remove(did)  # so that a later decision point can start the flow again
            if alive:
                self._resolve(c, born=False)
            return
        part = c.count * d.fraction
        c.count -= part
        self.fired[did] += 1
        self._pour(c, d.to, part)
        c.flow_at[did] = self.time + (d.after or 0.0) * self.organism.tempo
        self._push(c.flow_at[did], c.name, f"flow:{did}")
        self._record()
        self._resolve(c, born=False)  # a changed population decides again (caps and quiescence are re-read)

    # ---- events ----------------------------------------------------------

    def _add(self, c: Cell, parent: Cell | None) -> None:
        if c.name in self.cells:
            raise ValueError(f"cell {c.name!r} already exists (division naming collision)")
        self.cells[c.name] = c
        if parent is not None:
            parent.children.append(c.name)
        self._apply_signals(c)
        self._resolve(c)

    def _read_gradients(self, c: Cell) -> bool:
        """Set or clear the factors of gradient signals from the fields at the cell's site."""
        if c.x is None or not self._gradients:
            return False
        changed = False
        ctx = None
        for sg in self._gradients:
            if {sg.id, sg.field_name, sg.sets} & self.knockouts or sg.field_name not in self.fields:
                continue
            if sg.receiver:
                ctx = ctx or self.context(c)
                if not matches(sg.receiver, ctx):
                    continue
            on = self.fields[sg.field_name].at(c.x, c.y) >= sg.threshold
            if on and c.factors.get(sg.sets) != sg.value:
                c.factors[sg.sets] = sg.value
                self.fired[sg.id] += 1
                changed = True
            elif not on and sg.sets in c.factors:
                del c.factors[sg.sets]
                changed = True
        return changed

    def _sense(self) -> None:
        """Every `sense` minutes, positioned cells re-read the gradients and decide again if one changed."""
        for c in self.alive_at(self.time):
            if self._read_gradients(c):
                self._resolve(c, born=False)
        self._push(self.time + self.organism.sense, "", "sense")

    def _advance_fields(self, t: float) -> None:
        """Step the fields (hours) from the last field time to t (minutes), within the stable step."""
        if not self.fields or t <= self.field_time:
            return
        hours = (t - self.field_time) / 60.0
        self.field_time = t
        max_d = max(f.diffusion for f in self.fields.values()) or 1.0
        h = min(0.25 / max_d, 0.5)
        while hours > 1e-9:
            dt = min(h, hours)
            for name, x, y, rate in self.sources:
                self.fields[name].add(x, y, rate * dt)
            for f in self.fields.values():
                f.step(dt)
            hours -= dt

    def _free_site(self, c: Cell, direction: str = "") -> tuple[int, int] | None:
        """A free site for a daughter: along `direction` first, then the nearest free site."""
        assert c.x is not None and c.y is not None
        w, hgt = self.organism.width, self.organism.height
        steps = {"+x": (1, 0), "-x": (-1, 0), "+y": (0, 1), "-y": (0, -1)}

        def free(x: int, y: int) -> bool:
            return 0 <= x < w and 0 <= y < hgt and (x, y) not in self.occupied

        if direction in steps:
            dx, dy = steps[direction]
            x, y = c.x + dx, c.y + dy
            while 0 <= x < w and 0 <= y < hgt:  # the nearest free site in that direction
                if (x, y) not in self.occupied:
                    return (x, y)
                x, y = x + dx, y + dy
        seen = {(c.x, c.y)}
        frontier = [(c.x, c.y)]
        while frontier:
            nxt = []
            for x, y in frontier:
                for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    site = (x + dx, y + dy)
                    if site in seen or not (0 <= site[0] < w and 0 <= site[1] < hgt):
                        continue
                    if site not in self.occupied:
                        return site
                    seen.add(site)
                    nxt.append(site)
            frontier = nxt
        return None

    def _move(self, c: Cell, d: Decision) -> bool:
        """Migrate by `steps` sites, along `direction` or up the gradient of `toward`; free sites only."""
        assert c.x is not None and c.y is not None
        w, hgt = self.organism.width, self.organism.height
        steps = {"+x": (1, 0), "-x": (-1, 0), "+y": (0, 1), "-y": (0, -1)}
        moved = False
        for _ in range(max(1, d.steps)):
            target = None
            if d.direction in steps:
                dx, dy = steps[d.direction]
                target = (c.x + dx, c.y + dy)
            elif d.toward in self.fields:
                fld = self.fields[d.toward]
                best = fld.at(c.x, c.y)
                for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    x, y = c.x + dx, c.y + dy
                    if 0 <= x < w and 0 <= y < hgt and fld.at(x, y) > best:
                        best, target = fld.at(x, y), (x, y)
            if target is None or not (0 <= target[0] < w and 0 <= target[1] < hgt) or target in self.occupied:
                break
            del self.occupied[(c.x, c.y)]
            c.x, c.y = target
            self.occupied[target] = c.name
            moved = True
        return moved

    def _step_move(self, c: Cell, did: str) -> None:
        c.move_at = None
        d = next((x for x in self.module.decisions if x.id == did), None)
        alive = c.born <= self.time < c.end
        if d is None or not alive or not d.applies(self.context(c)):
            if did in c.fired:
                c.fired.remove(did)
            if alive:
                self._resolve(c, born=False)  # the walk is over: other decisions get their turn
            return
        if self._move(c, d):
            self.fired[did] += 1
            if self._read_gradients(c):
                self._resolve(c, born=False)  # a new reading at the new site may change the decision
                if did not in c.fired:
                    return  # the walk decision no longer applies after re-deciding
        c.move_at = self.time + (d.after or 0.0) * self.organism.tempo
        self._push(c.move_at, c.name, f"move:{did}")

    def _apply_signals(self, newborn: Cell) -> None:
        self._read_gradients(newborn)
        signals = [sg for sg in self.module.signals() if sg.mode != "gradient"]
        if not signals:
            return
        alive = [x for x in self.cells.values() if x.born <= self.time < x.end and x is not newborn]
        ctx_new = self.context(newborn)
        for sg in signals:
            if not sg.sets or {sg.id, sg.ligand, sg.receptor, sg.sets} & self.knockouts:
                continue
            if _match(sg.receiver, ctx_new) and any(_match(sg.sender, self.context(x)) for x in alive):
                newborn.factors[sg.sets] = sg.value
                self.fired[sg.id] += 1
            if _match(sg.sender, ctx_new):
                for x in alive:
                    if _match(sg.receiver, self.context(x)) and x.factors.get(sg.sets) != sg.value:
                        x.factors[sg.sets] = sg.value
                        self.fired[sg.id] += 1
                        self._resolve(x, born=False)

    def _divide(self, c: Cell) -> None:
        d = next((x for x in self._candidates(c) if x.id in c.fired and x.action == "divide"), None)
        if d is None:
            return
        if c.quiescent:
            return
        if self.population(c) and not d.daughters:
            # a population: grows in place by `fraction` (1.0 = doubling), then decides again
            c.count *= 1.0 + d.fraction
            c.divides_at = None
            c.fired = [x for x in c.fired if x not in self._divide_ids]  # the next step is decided afresh
            self._record()
            self._resolve(c)
            return
        if d.daughters:
            names = list(d.daughters)
        else:
            a, b = SUFFIXES[c.generation % 2]
            names = [c.name + a, c.name + b]
        sites: list[tuple[int, int] | None] = [None, None]
        if c.x is not None:
            second = self._free_site(c, d.direction)
            if second is None:  # contact inhibition: no room, the cell waits another cycle
                self.blocked_divisions += 1
                c.divides_at = None
                c.fired = [x for x in c.fired if x not in self._divide_ids]
                self._resolve(c)
                return
            sites = [(c.x, c.y), second]
            del self.occupied[(c.x, c.y)]
        for i, name in enumerate(names):
            factors = dict(c.factors)
            stated = set(c.stated)
            for factor, keeper in d.asymmetric.items():
                keep = (
                    keeper == name or (name.endswith(keeper) and len(keeper) == 1) or keeper == ("a", "p")[i]
                )
                if keep and factor not in self.knockouts:
                    factors[factor] = factors.get(factor, "present")
                    stated.add(factor)
                else:
                    factors.pop(factor, None)
                    stated.discard(factor)
            lineage = d.lineages.get(name, c.lineage)
            generation = 0 if name in d.lineages else c.generation + 1
            child = Cell(
                name, lineage, generation, self.time, c.cell_type, factors, parent=c.name, count=c.count
            )
            child.population = c.population
            child.measured = {f for f in c.measured if f in factors}
            child.stated = {f for f in stated if f in factors}
            if sites[i] is not None:
                child.x, child.y = sites[i]
                self.occupied[sites[i]] = name
            self._add(child, c)

    def _cull(self, c: Cell) -> None:
        """A share of a population dies (turnover); recurs while a fractional death decision applies."""
        c.cull_at = None
        d = self._first(c, "die", self.context(c))
        if d is None or d.fraction >= 1.0 or not (c.born <= self.time < c.end):
            return
        lost = c.count * d.fraction
        c.count -= lost
        self.culled += lost
        self.fired[d.id] += 1
        if d.id not in c.fired:
            c.fired.append(d.id)
        c.cull_at = self.time + (d.after or 0.0) * self.organism.tempo
        self._push(c.cull_at, c.name, "cull")
        self._record()
        self._resolve(c, born=False)  # a smaller population may grow again (quiescence is re-read)

    def run(self, until: float = math.inf) -> Body:
        while self._queue and len(self.cells) < self.max_cells:
            t, seq, name, kind = self._queue[0]
            if t > until:
                break
            heapq.heappop(self._queue)
            self.time = t
            self._advance_fields(t)
            if kind == "stage":
                self._stage_change()
                continue
            if kind == "sense":
                self._sense()
                continue
            c = self.cells[name]
            if (
                kind == "divide"
                and seq == c.divide_seq
                and c.divides_at == t
                and (c.dies_at is None or c.dies_at > t)
            ):
                self._divide(c)
            elif kind == "cull" and c.dies_at is None:
                self._cull(c)
            elif kind.startswith("flow:"):
                self._flow(c, kind[5:])
            elif kind.startswith("move:"):
                self._step_move(c, kind[5:])
            elif kind == "die" and c.dies_at == t and c.x is not None:
                self.occupied.pop((c.x, c.y), None)  # the site is free again
            # deaths need no action: dies_at already ends the cell
        self.time = until if until != math.inf else self.time
        self._advance_fields(self.time)
        return self

    # ---- space -----------------------------------------------------------

    def positions(self) -> list[tuple[str, int, int, str]]:
        return [(c.name, c.x, c.y, c.cell_type) for c in self.alive_at(self.time) if c.x is not None]

    def type_map(self) -> list[str]:
        """One character per site, row by row: first letter of the cell type, '.' for an empty site."""
        if not self.spatial:
            return []
        grid = [["."] * self.organism.width for _ in range(self.organism.height)]
        for _, x, y, cell_type in self.positions():
            grid[y][x] = (cell_type or "?")[0]
        return ["".join(row) for row in grid]

    def bands_along_x(self, y: int | None = None) -> list[tuple[str, int, int]]:
        """Contiguous runs of cell type along x in one row: (type, start, end)."""
        if not self.spatial:
            return []
        row = self.organism.height // 2 if y is None else y
        types = ["" for _ in range(self.organism.width)]
        for _, x, yy, cell_type in self.positions():
            if yy == row:
                types[x] = cell_type
        bands: list[tuple[str, int, int]] = []
        start = 0
        for x in range(1, self.organism.width + 1):
            if x == self.organism.width or types[x] != types[start]:
                bands.append((types[start], start, x))
                start = x
        return bands

    # ---- queries ---------------------------------------------------------

    def alive_at(self, t: float) -> list[Cell]:
        return [c for c in self.cells.values() if c.born <= t < c.end]

    def count_at(self, t: float) -> float:
        """Cells alive at t. Population counts are mutable, so for past times the recorded total is used."""
        if self.history and t < self.time:
            i = bisect.bisect_right(self.history, (t, math.inf)) - 1
            if i >= 0:
                return self.history[i][1]
        return sum(c.count for c in self.alive_at(t))

    def lineage_count_at(self, lineage: str, t: float) -> float:
        return sum(c.count for c in self.alive_at(t) if c.lineage == lineage)

    def deaths_by(self, t: float) -> int:
        return sum(1 for c in self.cells.values() if c.dies_at is not None and c.dies_at <= t)

    def turnover_per_day(self) -> float:
        """Cells lost per day by the fractional death decisions currently applying to populations."""
        total = 0.0
        for c in self.alive_at(self.time):
            if not self.population(c):
                continue
            d = self._first(c, "die", self.context(c))
            if d is not None and d.fraction < 1.0 and d.after:
                total += c.count * d.fraction / (d.after / 1440.0)
        return total

    def fates_at(self, t: float) -> dict[str, float]:
        out: dict[str, float] = {}
        for c in self.alive_at(t):
            if (c.divides_at is None and c.dies_at is None) or c.quiescent:  # terminal, quiescent or UNKNOWN
                key = c.cell_type or "UNKNOWN"
                out[key] = out.get(key, 0) + c.count
        return dict(sorted(out.items(), key=lambda kv: -kv[1]))

    def tree(self, root: str | None = None, depth: int = 3) -> list[str]:
        out: list[str] = []

        def walk(name: str, d: int) -> None:
            c = self.cells[name]
            what = ""
            if c.divides_at is not None:
                what = f"  divides {c.divides_at:.0f} min"
            elif c.dies_at is not None:
                what = f"  dies {c.dies_at:.0f} min"
            elif c.cell_type:
                what = f"  -> {c.cell_type}" + (f" ({c.terminal_name})" if c.terminal_name else "")
            if c.unknown:
                what += f"  UNKNOWN: {c.unknown}"
            n = f" ×{c.count:.4g} [{c.cell_type}]" if self.population(c) else ""
            out.append(f"{'  ' * d}{c.name}{n} (born {c.born:.0f} min){what}")
            if d < depth:
                for ch in c.children:
                    walk(ch, d + 1)

        walk(root or self.organism.root, 0)
        return out

    def uncertainty(self) -> UncertaintyReport:
        rep = UncertaintyReport()
        by_id = {d.id: d for d in self.module.decisions}
        for sid, n in self.fired.items():
            d = by_id.get(sid)
            if d is None:
                continue
            for _ in range(n):
                rep.add("organism", d.confidence, d.evidence.kind, d.id)
        by_name = {t.name: t for t in self.module.timers}
        for tid, n in self.timers_used.items():
            t = by_name[tid]
            for _ in range(n):
                rep.add("cellular", t.confidence, t.evidence.kind, t.name)
        for r in self.module.rules:
            rep.add_rule(r)
        for reason, n in self.unknown.items():
            for _ in range(n):
                rep.add("organism", 0.0, EvidenceKind.NONE, f"UNKNOWN: {reason}")
        return rep

    def check_asserts(self) -> list[dict]:
        return [evaluate_assert(self, a) for a in self.organism.asserts]

    def summary(self) -> dict:
        t = self.time
        return {
            "organism": self.organism.name,
            "species": self.organism.species,
            "time_min": t,
            "cells_born": len(self.cells),
            "alive": self.count_at(t),
            "deaths": self.deaths_by(t),
            "culled": self.culled,
            "turnover_per_day": self.turnover_per_day(),
            "populations": sum(1 for c in self.alive_at(t) if self.population(c)),
            "fates": self.fates_at(t),
            "decisions_fired": sum(self.fired.values()),
            "unknown": dict(self.unknown),
        }


def _match(when: dict[str, str], ctx: dict[str, str]) -> bool:
    return bool(when) and matches(when, ctx)


# ---- asserts ---------------------------------------------------------------


def evaluate_assert(body: Body, text: str) -> dict:
    """`count at 100 min in 24..34` · `deaths at 800 min = 113` · `fate Neuron at 800 min >= 200`
    · `type EPrecursor at 100 min = 2` (alive cells of that type, dividing or not)
    · `lineage AB at 120 min >= 8`."""
    m = re.match(
        r"^(count|deaths|fate|type|lineage)\s*(\S+)?\s+at\s+([-0-9.]+)\s*(\w+)?\s*(in|=|>=|<=|>|<)\s*(.+)$",
        text,
    )
    if not m:
        return {"assert": text, "ok": False, "error": "cannot parse"}
    metric, arg, t, unit, op, rhs = m.groups()
    t_min = to_minutes(float(t), unit or "min")
    if metric == "count":
        value = body.count_at(t_min)
    elif metric == "deaths":
        value = body.deaths_by(t_min)
    elif metric == "fate":
        value = body.fates_at(t_min).get(arg or "", 0)
    elif metric == "type":
        value = sum(c.count for c in body.alive_at(t_min) if c.cell_type == (arg or ""))
    else:
        value = body.lineage_count_at(arg or "", t_min)
    if op == "in":
        lo, hi = (float(x) for x in rhs.split(".."))
        ok = lo <= value <= hi
    else:
        n = float(rhs)
        ok = {"=": value == n, ">=": value >= n, "<=": value <= n, ">": value > n, "<": value < n}[op]
    return {"assert": text, "ok": bool(ok), "value": value}
