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
        self._bootstrap()

    # ---- setup ---------------------------------------------------------

    def _bootstrap(self) -> None:
        o = self.organism
        for st in self.module.stages:
            self._push(to_minutes(st.start, st.unit), "", "stage")
        factors = {f: "present" for f in list(o.factors) + sorted(self.adds) if f not in self.knockouts}
        root = Cell(o.root, "", 0, 0.0, o.cell_type, factors)
        root.population = o.resolution == "populations"
        self._add(root, None)

    def _stage_change(self) -> None:
        """Populations and resting cells read the new stage and decide again."""
        for c in self.alive_at(self.time):
            if self.population(c) or c.quiescent:
                self._resolve(c, born=False)

    def _push(self, t: float, name: str, kind: str) -> None:
        self._seq += 1
        heapq.heappush(self._queue, (t, self._seq, name, kind))

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
        for _ in range(32):  # differentiation chains and sequential population splits
            d = self._first(c, "differentiate", ctx, unfired=True)
            if d is None:
                break
            if d.fraction < 1.0 and self.population(c):
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
        if c.divides_at is not None and not born and not c.quiescent:
            # a signal changed what the cell reads: an earlier-precedence division may now apply
            d = self._first(c, "divide", ctx)
            current = next((x for x in c.fired if x in self._divide_ids), None)
            if d is not None and current is not None and d.id != current:
                c.fired.remove(current)
                self.fired[current] -= 1
                c.divides_at = None
        if c.divides_at is None and not c.quiescent and (born or self.population(c)):
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
            # cells wait from birth; populations wait from now (they decide again after every step)
            c.divides_at = (self.time if self.population(c) else c.born) + wait
            c.fired.append(d.id)
            self.fired[d.id] += 1
            self._push(c.divides_at, c.name, "divide")

    def _record(self) -> None:
        total = sum(x.count for x in self.cells.values() if x.born <= self.time < x.end)
        if self.history and self.history[-1][0] == self.time:
            self.history[-1] = (self.time, total)
        else:
            self.history.append((self.time, total))

    def _split(self, c: Cell, d: Decision) -> None:
        """A fraction of a population differentiates into a new node."""
        part = c.count * d.fraction
        c.count -= part
        name = f"{c.name}>{d.to}"
        child = Cell(
            name, c.lineage, c.generation, self.time, d.to, dict(c.factors), parent=c.name, count=part
        )
        child.population = True
        self._add(child, c)
        self._record()

    # ---- events ----------------------------------------------------------

    def _add(self, c: Cell, parent: Cell | None) -> None:
        if c.name in self.cells:
            raise ValueError(f"cell {c.name!r} already exists (division naming collision)")
        self.cells[c.name] = c
        if parent is not None:
            parent.children.append(c.name)
        self._apply_signals(c)
        self._resolve(c)

    def _apply_signals(self, newborn: Cell) -> None:
        signals = self.module.signals()
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
        for i, name in enumerate(names):
            factors = dict(c.factors)
            for factor, keeper in d.asymmetric.items():
                keep = (
                    keeper == name or (name.endswith(keeper) and len(keeper) == 1) or keeper == ("a", "p")[i]
                )
                if keep and factor not in self.knockouts:
                    factors[factor] = factors.get(factor, "present")
                else:
                    factors.pop(factor, None)
            lineage = d.lineages.get(name, c.lineage)
            generation = 0 if name in d.lineages else c.generation + 1
            child = Cell(
                name, lineage, generation, self.time, c.cell_type, factors, parent=c.name, count=c.count
            )
            child.population = c.population
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
            t, _, name, kind = self._queue[0]
            if t > until:
                break
            heapq.heappop(self._queue)
            self.time = t
            if kind == "stage":
                self._stage_change()
                continue
            c = self.cells[name]
            if kind == "divide" and c.divides_at == t and (c.dies_at is None or c.dies_at > t):
                self._divide(c)
            elif kind == "cull" and c.dies_at is None:
                self._cull(c)
            # deaths need no action: dies_at already ends the cell
        self.time = until if until != math.inf else self.time
        return self

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
