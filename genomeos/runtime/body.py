# SPDX-License-Identifier: Apache-2.0
"""Body runtime: grows an organism from one cell by discrete events.

Every cell (or population: a cell with `count` > 1) reads its context, the
program's decisions say what it does, timers say when, and signals from
other cells change what it reads. The loop per cell:

    read     context = name, lineage, generation, cell type, stage, factors, environment
    decide   die | differentiate | quiesce | divide     (by `priority`, then module order; for differentiate
             see `regime fates`: first (default) takes one fate per decision point, last is the legacy mode)
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
from pathlib import Path

from genomeos.ir import Decision, EvidenceKind, Module, Timer, matches, to_minutes

from .spatial import Field2D
from .uncertainty import UncertaintyReport

SUFFIXES = (("a", "p"), ("l", "r"))
FORCE_LAST = 1 << 62  # event tie-break: a forced factor arrives after everything else at its instant
# an integrated read of a factor over a named window (v0.4 §7.2a): `ELT-2.exposure(lineage)`
READ = re.compile(r"^(.+)\.(exposure|mean)\((cell|lineage)\)$")
# a threshold this close is reached, and a wake-up this soon is now: without both, a crossing computed
# to land exactly on its threshold can miss it by a float's width and reschedule itself for ever
REACHED = 1e-9
SOONEST = 1e-6  # minutes


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
    committed: str = ""  # the programme this cell is committed to, once a `commitment` established (§7.2a)
    fate_priority: int | None = None  # precedence of the decision that settled this cell's fate (§7.3)
    fate_ctx: dict[str, str] | None = None  # what the cell read then, so a loser cannot try again on it
    refused: set[str] = field(default_factory=set)  # fate decisions refused for this cell, counted once
    # integrated reads (v0.4 §7.2a): exposure banked in this cell's own life, what the path to it carried
    # before it was born, and the values in force since `exp_t`
    own: dict[str, float] = field(default_factory=dict)
    inherited: dict[str, float] = field(default_factory=dict)
    exp_t: float = 0.0
    exp_val: dict[str, float] = field(default_factory=dict)
    recheck_at: float | None = None  # when a read this cell is waiting on reaches a threshold (§7.5)
    recheck_seq: int = -1  # event token: only the latest scheduled crossing is honoured
    measured: set[str] = field(default_factory=set)  # factors set by express decisions (the reader)
    stated: set[str] = field(default_factory=set)  # factors stated by mechanism: maternal load, asymmetry
    levels: dict[str, float] = field(default_factory=dict)  # v0.4: this cell's own network state

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
        fates: str | None = None,
        add_at: dict[str, float] | None = None,
        recheck: str | None = None,
    ):
        """`seed` overrides the organism's declared seed; `means=True` runs every timer at its mean;
        `knockouts` are factors never present and signals (by id, ligand or receptor) never sent;
        `adds` are extra factors in the zygote; `add_at` forces a factor at a stated time instead
        (every cell alive then, and every cell born after); `environment` overrides the organism's."""
        if module.organism is None:
            raise ValueError(f"module {module.name!r} declares no organism block")
        self.module = module
        self.organism = module.organism
        self.knockouts = set(knockouts)
        self.add_at = {k: v for k, v in (add_at or {}).items() if k not in self.knockouts}
        self.adds = {f for f in adds if f not in self.add_at}
        self.forced: set[str] = set()  # factors already forced, carried by every cell born afterwards
        self.forced_cells: Counter[str] = Counter()  # cells the forcing reached, per factor
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
        self.ambiguous: Counter[str] = Counter()  # decision points where equal-precedence fates disagreed
        self.revised: Counter[str] = Counter()  # terminal fates changed again at a later decision point
        # fate changes the program refused rather than applied (v0.4 §7.2a), by decision id
        self.outside_competence: Counter[str] = Counter()
        self.refused_committed: Counter[str] = Counter()
        # fate changes refused because a higher-precedence decision had already settled the fate (§7.3)
        self.overruled: Counter[str] = Counter()
        self.windows = {c.name: c for c in module.competences}
        # every integrated read the program asks for, collected once: (key, factor, kind, window)
        guards = [d.when for d in module.decisions] + [ct.establish for ct in module.commitments]
        guards += [sg.receiver for sg in module.signals()] + [sg.sender for sg in module.signals()]
        self._reads = sorted(
            {
                (key, m.group(1), m.group(2), m.group(3))
                for when in guards
                for key in when
                if (m := READ.match(key))
            }
        )
        # §7.5: the thresholds each decision waits for, so a cell can be woken at the instant one is
        # reached. Only lower bounds schedule anything: an integral that has not yet reached a floor
        # will reach it at a computable time, while `< T` is true until it stops being true and needs
        # no wake-up (a fate is never withdrawn).
        self._terms: dict[int, list[tuple[str, str, str, float]]] = {}
        for d in module.decisions:
            terms = []
            for key, want in d.when.items():
                m = READ.match(key)
                op = want[:2] if want[:2] == ">=" else want[:1]
                if m and op in (">=", ">"):
                    try:
                        terms.append((m.group(1), m.group(2), m.group(3), float(want[len(op) :])))
                    except ValueError:
                        continue
            if terms:
                self._terms[id(d)] = terms
        self.recheck = recheck or (module.regime.recheck if module.regime is not None else "crossings")
        if self.recheck not in ("crossings", "none"):
            raise ValueError(f"recheck must be crossings or none, not {self.recheck!r}")
        self.rechecks = 0  # re-decisions taken at a crossing
        self.committed_cells: Counter[str] = Counter()  # cells committed, by programme
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
        # by identity, never by id: two decisions may share an id (mechanism stated before the generated
        # lookup, as the worm's founders do), and indexing by id dropped the one without a `cell` clause
        # out of every candidate list, so it never fired and nothing said so.
        named = {id(x) for v in self._by_cell.values() for x in v}
        self._general = [d for d in module.decisions if id(d) not in named]
        # a shared id is legal but not free: a cell's `fired` list records ids, so two decisions with one
        # name cannot be told apart there. Counted here and reported, like ambiguous_fates.
        self.duplicate_ids = {did: n for did, n in Counter(d.id for d in module.decisions).items() if n > 1}
        # decisions on one cell type are indexed by it; the rest apply to any type. Candidate lists are
        # rebuilt in module order per (cell name, cell type) so precedence is unchanged, and cached.
        self._order = {id(d): i for i, d in enumerate(module.decisions)}
        self._by_type: dict[str, list[Decision]] = {}
        self._untyped: list[Decision] = []
        for d in self._general:
            ct = d.when.get("cell_type")
            if ct and ct != "any" and "|" not in ct and ct[:1] not in ("<", ">") and ct != "absent":
                self._by_type.setdefault(ct, []).append(d)
            else:
                self._untyped.append(d)
        self._cand_cache: dict[tuple[str, str], list[Decision]] = {}
        self._stage_cache: tuple[float, str] = (-1.0, "")
        self._divide_ids = {d.id for d in module.decisions if d.action == "divide"}
        self._expressed_names = {f for d in module.decisions if d.action == "express" for f in d.sets}
        # how a decision point takes a fate, declared in the regime (docs/BIOLANG-v0.4-ECONOMY.md §7.3)
        # `fates` overrides the program; default first since 2026-09-14 (spec §10 decision 8)
        self.fates = fates or (module.regime.fates if module.regime is not None else "first")
        if self.fates not in ("first", "last"):
            raise ValueError(f"fates must be first or last, not {self.fates!r}")
        self._init_contacts(seed, means)
        if self.fates == "last" and any(d.priority for d in module.decisions):
            raise ValueError(
                "decision `priority` needs `regime { fates: first }`: with fates: last the final match wins"
            )
        self._bootstrap()

    # ---- setup ---------------------------------------------------------

    def _bootstrap(self) -> None:
        o = self.organism
        for st in self.module.stages:
            self._push(to_minutes(st.start, st.unit), "", "stage")
        for when in sorted(set(self.add_at.values())):
            # last at its instant: a cell dividing exactly then has already been replaced by its
            # daughters, so a perturbation timed on a division boundary misses no lineage
            heapq.heappush(self._queue, (when, FORCE_LAST, "", "force"))
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

    def _force(self, t: float) -> None:
        """A factor forced at a stated time (`add: HLH-1 at 350 min`): every cell alive now gains it as
        stated mechanism, which the reader cannot drop, and decides again; cells born later inherit it.
        This is the perturbation a competence window is measured with, so it must be able to arrive at
        any time, not only in the zygote."""
        names = sorted(f for f, when in self.add_at.items() if when == t)
        if not names:
            return
        self.forced.update(names)
        for c in self.alive_at(t):
            for f in names:
                c.factors[f] = "present"
                c.stated.add(f)
                c.measured.discard(f)
                self.forced_cells[f] += 1
            self._resolve(c, born=False)

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

    # ---- integrated reads (v0.4 §7.2a) -----------------------------------

    def _sync(self, c: Cell, t: float) -> None:
        """Bank the exposure a cell has accumulated since it was last read, then take the values now in
        force. Called before every read and at the end of every decision point, which is where factors
        change, so an interval is banked with the values that were actually in force during it."""
        t = min(t, c.end)  # a cell stops accumulating when it divides or dies, whenever it is read
        dt = t - c.exp_t
        if dt > 0:
            for f, v in c.exp_val.items():
                if v:
                    c.own[f] = c.own.get(f, 0.0) + v * dt
        c.exp_t = t
        vals: dict[str, float] = dict.fromkeys(c.factors, 1.0)  # a factor the cell carries counts as 1
        vals.update(c.levels)  # a species with a numeric level is integrated at that level instead
        c.exp_val = vals

    def _read(self, c: Cell, factor: str, kind: str, window: str, t: float) -> float:
        """`F.exposure(cell)` minutes carrying F in this cell's own life; `F.exposure(lineage)` the same
        summed along the path from the first cell; `F.mean(...)` that divided by the window's length, so
        a plain factor reads as the fraction of the window it was carried and a species reads as its
        mean level. Area E measured that the instantaneous read is the worst of the three on 555 terminal
        cells (62 errors against 38 and 29), which is why the language has all of them and none is a
        default: the window is declared, because the three do not agree."""
        t = min(t, c.end)
        total = c.own.get(factor, 0.0)
        span = max(t - c.born, 0.0)
        if window == "lineage":
            total += c.inherited.get(factor, 0.0)
            span = max(t, 0.0)  # the path starts with the first cell, at time zero
        if kind == "exposure":
            return total
        # over a window of no length the mean is its own limit, the value in force: a cell reading
        # `mean(cell)` at the instant it is born is taking the instantaneous read, and says so
        return total / span if span > 0 else c.exp_val.get(factor, 0.0)

    def _crossing(self, c: Cell, factor: str, kind: str, window: str, threshold: float, t: float) -> float:
        """How long until this read reaches this threshold at the rate now in force, or infinity if it
        never does. The runtime integrates each read at a rate held constant between decision points
        (§7.2a), so the answer is a closed form rather than a search: exposure grows at `r`, and a mean
        of a growing total over a growing window rises only while `r` exceeds the threshold itself."""
        total = c.own.get(factor, 0.0) + (c.inherited.get(factor, 0.0) if window == "lineage" else 0.0)
        rate = c.exp_val.get(factor, 0.0)
        if kind == "exposure":
            if total >= threshold - REACHED:
                return math.inf  # already met: this decision point has considered it
            return (threshold - total) / rate if rate > 0 else math.inf
        span = max(t if window == "lineage" else t - c.born, 0.0)
        if (total / span if span > 0 else rate) >= threshold - REACHED:
            return math.inf
        if rate <= threshold:
            return math.inf  # the window grows at least as fast as the total, so the mean cannot rise to it
        return max((threshold * span - total) / (rate - threshold), 0.0)

    def _schedule_recheck(self, c: Cell) -> None:
        """One pending re-decision per cell: the earliest instant at which any threshold a decision it
        could still take is waiting for would be reached. A crossing is a re-reading, not a new right -
        precedence and commitment refuse at a crossing exactly what they refuse anywhere else."""
        t = self.time
        if not (c.born <= t < c.end):
            return
        self._sync(c, t)
        soonest = math.inf
        for d in self._candidates(c):
            terms = self._terms.get(id(d))
            if terms is None or d.id in c.fired:
                continue
            for factor, kind, window, threshold in terms:
                dt = self._crossing(c, factor, kind, window, threshold, t)
                if SOONEST < dt < soonest:
                    soonest = dt
        if soonest == math.inf:
            c.recheck_at = None
            return
        c.recheck_at = t + soonest
        c.recheck_seq = self._push(c.recheck_at, c.name, "recheck")

    # ---- context and decisions ------------------------------------------

    def context(self, c: Cell, t: float | None = None, amounts: bool = True) -> dict[str, str]:
        t = self.time if t is None else t
        if self._stage_cache[0] != t:
            self._stage_cache = (t, self.module.stage_at(t))
        ctx = {
            "cell": c.name,
            "lineage": c.lineage,
            "generation": str(c.generation),
            "cell_type": c.cell_type,
            "stage": self._stage_cache[1],
            "count": str(c.count),
        }
        ctx.update(self.environment)
        if c.levels:
            ctx.update({k: f"{v:.6g}" for k, v in c.levels.items()})
        ctx.update(c.factors)
        if self._reads:  # only programs that ask for an integrated read pay for one
            self._sync(c, t)
            for key, factor, kind, window in self._reads:
                ctx[key] = f"{self._read(c, factor, kind, window, t):.6g}"
        if amounts and self._amount_signals and self.network is None:
            for sg in self._amount_signals:  # without networks an amount is the count of touching senders
                if not {sg.id, sg.ligand, sg.receptor, sg.sets} & self.knockouts:
                    ctx[sg.sets] = f"{self.received(c, sg):.6g}"
        if c.x is not None:
            ctx["x"], ctx["y"] = str(c.x), str(c.y)
        return ctx

    def _candidates(self, c: Cell) -> list[Decision]:
        key = (c.name if c.name in self._by_cell else "", c.cell_type)
        lst = self._cand_cache.get(key)
        if lst is None:
            lst = self._by_cell.get(key[0], []) + self._by_type.get(c.cell_type, []) + self._untyped
            lst.sort(
                key=lambda d: (-d.priority, self._order[id(d)])
            )  # explicit precedence, then module order
            self._cand_cache[key] = lst
        return lst

    def _first(self, c: Cell, action: str, ctx: dict[str, str], unfired: bool = False) -> Decision | None:
        for d in self._candidates(c):
            if d.action == action and d.applies(ctx) and not (unfired and d.id in c.fired):
                return d
        return None

    def _pick_fate(self, c: Cell, ctx: dict[str, str], before: dict[str, str] | None) -> Decision | None:
        """The one `differentiate` a decision point takes: the first matching decision by precedence
        (`priority`, then module order). After a change of type the chain continues only through a decision
        the new type enables; one that already applied before the change is a competitor, not a successor,
        so a later rule can no longer silently overwrite an earlier one. Equal-precedence matches that
        disagree on the target are counted as ambiguous."""
        chosen: Decision | None = None
        for d in self._candidates(c):
            if d.action != "differentiate" or d.id in c.fired or not d.applies(ctx):
                continue
            if before is not None and d.applies(before):
                continue
            if d.competence and not self._competent(c, d, ctx):
                self._refuse(c, d, self.outside_competence)  # the program said no, which is not UNKNOWN
                continue
            if (
                self.fates == "first"
                and c.fate_priority is not None
                and d.priority < c.fate_priority
                and d.to != c.cell_type
                and c.fate_ctx is not None
                and d.applies(c.fate_ctx)
                and not (self.population(c) and d.fraction < 1.0)
            ):
                # the same rule as inside a decision point, carried across them: a decision that could
                # already have applied when the fate was settled, and lost on precedence, is a competitor
                # and not a successor, so it cannot take the fate back when the cell decides again (a
                # `cell_network` step, a signal, a stage). A decision whose guard was false then is a new
                # reading and still applies. Populations' shares are splits, not competitors.
                self._refuse(c, d, self.overruled)
                continue
            if c.committed and d.to != c.committed:
                self._refuse(c, d, self.refused_committed)
                continue
            if chosen is None:
                chosen = d
                continue
            share = self.population(c) and (d.fraction < 1.0 or chosen.fraction < 1.0)
            if not share and d.priority == chosen.priority and d.to != chosen.to:
                self.ambiguous[chosen.id] += 1  # a population's shares are meant to split, not compete
            break
        return chosen

    def _refuse(self, c: Cell, d: Decision, counter: Counter[str]) -> None:
        """A fate the program refused, counted once per cell and decision so a cell that decides again
        is not counted again. Never silently applied and never silently dropped (v0.4 §7.2a)."""
        if d.id not in c.refused:
            c.refused.add(d.id)
            counter[d.id] += 1

    def _competent(self, c: Cell, d: Decision, ctx: dict[str, str]) -> bool:
        """Is the window this fate change needs still open for this cell? A window governs the cells its
        `when` matches; it closes at a time, after a generation or on commitment, and `closed_by` names
        the machinery that closes it, so a cell without that machinery stays competent."""
        w = self.windows.get(d.competence)
        if w is None or (w.when and not matches(w.when, ctx)):
            return True
        if w.closed_by and ctx.get(w.closed_by) != "present":
            return True  # the closing machinery is not there, so the window does not close
        if w.closes_at is not None and self.time >= w.closes_at:
            return False
        if w.closes_generation is not None and c.generation > w.closes_generation:
            return False
        return not (w.closes_on_commitment and c.committed)

    def _commit(self, c: Cell, ctx: dict[str, str]) -> None:
        """A cell whose `establish` clause matches is committed to its programme; from then on a fate
        outside it is refused. The state is inherited by the daughters when `inherit: daughters`."""
        if c.committed:
            return
        for ct in self.module.commitments:
            if matches(ct.establish, ctx):
                c.committed = ct.programme or c.cell_type
                self.committed_cells[c.committed] += 1
                return

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
        """Decide what a cell does next; called at birth and after a signal changes its factors, and
        then (§7.5) put the cell down for a re-decision at the instant a read it is waiting on arrives."""
        self._decide(c, born)
        if self._terms and self.recheck == "crossings":
            self._schedule_recheck(c)

    def _decide(self, c: Cell, born: bool = True) -> None:
        """One decision point: die, express, differentiate, quiesce, migrate, divide."""
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
        before: dict[str, str] | None = None  # the context before the last change of type
        had_fate, changed = c.terminal_name, False
        for _ in range(32):  # differentiation chains and sequential population splits
            d = self._pick_fate(c, ctx, before if self.fates == "first" else None)
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
                # a decision that lands on the type the cell already holds revises nothing: it must not
                # be counted as a revision and must not take the cell's terminal name away, or a cell
                # re-read at a crossing (§7.5) would quietly stop answering to its name
                revision = bool(had_fate) and not changed and d.to != c.cell_type
                if revision:
                    self.revised[d.id] += 1  # visible until `commitment` can refuse it (v0.4 §7.2)
                before, changed = ctx, True
                c.cell_type = d.to
                # what a later decision point must beat, and what the losers of this one read
                c.fate_priority, c.fate_ctx = d.priority, ctx
                if d.name or revision:
                    c.terminal_name = d.name  # a revised cell no longer answers to its old terminal name
            c.fired.append(d.id)
            self.fired[d.id] += 1
            ctx = self.context(c)
        if self.module.commitments:
            self._commit(c, ctx)
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
            if d is None and current is not None:
                # no division applies to what the cell has become: the pending one is cancelled
                c.fired.remove(current)
                self.fired[current] -= 1
                c.divides_at = None
            elif d is not None and current is not None and d.id != current:
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
            c.levels = dict(parent.levels)
        elif self.network is not None:
            c.levels = {p.id: p.initial for p in self.module.proteins() if p.initial}
        self._apply_signals(c)
        self._resolve(c)
        if self._amount_signals and self.network is None:
            for name in self.neighbours(c):
                self._resolve(self.cells[name], born=False)

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
        signals = [sg for sg in self.module.signals() if sg.mode != "gradient" and sg.reads != "amount"]
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
            second = self._free_site(c, self._placement(d, names))
            if second is None:  # contact inhibition: no room, the cell waits another cycle
                self.blocked_divisions += 1
                c.divides_at = None
                c.fired = [x for x in c.fired if x not in self._divide_ids]
                self._resolve(c)
                return
            sites = [(c.x, c.y), second]
            del self.occupied[(c.x, c.y)]
        path: dict[str, float] = {}
        if self._reads:  # what the path to these daughters carried, banked up to this division
            self._sync(c, self.time)
            path = dict(c.inherited)
            for f, v in c.own.items():
                path[f] = path.get(f, 0.0) + v
        for i, name in enumerate(names):
            factors = dict(c.factors)
            stated = set(c.stated)
            for factor, keeper in d.asymmetric.items():
                keep = (
                    keeper == name or (name.endswith(keeper) and len(keeper) == 1) or keeper == ("a", "p")[i]
                )
                if keep and factor in factors and factor not in self.knockouts:
                    # a keeper keeps what the mother had; a mother without the factor passes nothing on,
                    # so `asymmetric` segregates a factor and never creates one out of nothing
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
            child.exp_t = self.time
            if path:
                child.inherited = dict(path)
            if c.committed and any(ct.inherit for ct in self.module.commitments):
                child.committed = c.committed  # `inherit: daughters`: the lock is mechanism, not bookkeeping
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

    # ---- contacts and per-cell networks (v0.4 §7.4) -----------------------------------------

    def _init_contacts(self, seed: int | None, means: bool) -> None:
        """Neighbours come from a measured contact table when the organism names one, otherwise from grid
        adjacency. With `cell_network`, every cell carries its own copy of the module's network, all copies
        are stepped together every `cell_network` minutes, and contact signals that read an amount couple
        them: a receiver's input is the ligand summed over its current neighbours, read before anyone moves
        (a synchronous update), so no cell is advantaged by the order cells are visited in."""
        o = self.organism
        self.contacts = read_contacts(o.contacts) if o.contacts else []
        self._contact_times = [t for t, _ in self.contacts]
        self._amount_signals = [
            sg for sg in self.module.signals() if sg.mode == "contact" and sg.reads == "amount" and sg.sets
        ]
        self.network = None
        self.network_steps = 0
        self._rngs: dict[str, random.Random] = {}
        self.network_noise = 0.0
        self.network_seed = None if means else seed
        if o.cell_network <= 0:
            return
        from .grn import NetworkRuntime

        self.network = NetworkRuntime(self.module)
        self.network_noise = float(self.network.params.get("noise", 0.0))
        if self.network_noise > 0 and self.network_seed is None:
            raise ValueError(
                "network noise needs a seed: declare `seed:` or pass one (noise is part of the mechanism)"
            )
        self.network.params["noise"] = 0.0  # noise is applied per cell below, from each cell's own stream
        self._push(o.cell_network, "", "network")

    def neighbours(self, c: Cell, t: float | None = None) -> dict[str, float]:
        """The cells touching `c` now, with contact area: the latest contact snapshot at or before t, or the
        four grid neighbours when there is no table. A cell not alive at t has no neighbours."""
        t = self.time if t is None else t
        if self.contacts:
            i = bisect.bisect_right(self._contact_times, t) - 1
            if i < 0:
                return {}
            snap = self.contacts[i][1].get(c.name, {})
            return {
                n: a
                for n, a in snap.items()
                if n in self.cells and self.cells[n].born <= t < self.cells[n].end
            }
        if c.x is None:
            return {}
        out = {}
        for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            name = self.occupied.get((c.x + dx, c.y + dy))
            if name is not None:
                out[name] = 1.0
        return out

    def received(self, c: Cell, sg) -> float:
        """The amount of a contact signal a cell receives: the ligand level in each neighbour that matches the
        sender condition, weighted by contact area and summed. Without per-cell networks the ligand amount is
        1 per matching neighbour, so the reading is a count of senders touching the cell."""
        total = 0.0
        for name, area in self.neighbours(c).items():
            other = self.cells[name]
            if sg.sender and not _match(sg.sender, self.context(other, amounts=False)):
                continue
            level = other.levels.get(sg.ligand, 0.0) if self.network is not None else 1.0
            total += area * level
        return total

    def _network_step(self) -> None:
        """Advance every cell's network by one `cell_network` interval, all together.

        Inputs are read from every cell before any cell changes, so the update is synchronous; each cell
        draws its noise from its own stream (seeded by the run's seed and the cell's name), so which cell
        wins a symmetric competition is decided by the seed and never by the order cells are stored in."""
        vm = self.network
        assert vm is not None
        hours = self.organism.cell_network / 60.0
        steps = max(1, math.ceil(hours / 0.01))
        dt = hours / steps
        live = [c for c in self.alive_at(self.time) if not self.population(c)]
        signals = [
            sg for sg in self._amount_signals if not {sg.id, sg.ligand, sg.receptor, sg.sets} & self.knockouts
        ]
        inputs = {c.name: {sg.sets: self.received(c, sg) for sg in signals} for c in live}
        for c in live:
            state = dict.fromkeys(vm.species, 0.0)
            state.update(c.levels)
            clamp = inputs[c.name]
            state.update(clamp)
            rng = self._rngs.get(c.name)
            if rng is None and self.network_noise > 0:
                rng = self._rngs[c.name] = random.Random(f"{self.network_seed}:{c.name}")
            for _ in range(steps):
                d = vm._derivatives(state)
                for s in vm.species:
                    if s in clamp:
                        continue
                    kick = 0.0
                    if rng is not None:
                        kick = (
                            self.network_noise
                            * math.sqrt(dt)
                            * rng.gauss(0.0, 1.0)
                            * math.sqrt(max(state[s], 1e-9))
                        )
                    state[s] = max(0.0, state[s] + d[s] * dt + kick)
            c.levels = state
        self.network_steps += 1
        for c in live:
            self._resolve(c, born=False)  # a level crossing a threshold is a new reading
        self._push(self.time + self.organism.cell_network, "", "network")

    def _placement(self, d: Decision, names: list[str]) -> str:
        """Where the second daughter goes: the decision's `direction`, else, with `placement: names`, the axis
        its name implies (a/p along x, l/r and d/v along y), else the nearest free site."""
        if d.direction or self.organism.placement != "names":
            return d.direction
        return {"a": "-x", "p": "+x", "l": "-y", "r": "+y", "d": "-y", "v": "+y"}.get(names[1][-1:], "")

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
            if kind == "force":
                self._force(t)
                continue
            if kind == "sense":
                self._sense()
                continue
            if kind == "network":
                self._network_step()
                continue
            c = self.cells[name]
            if (
                kind == "divide"
                and seq == c.divide_seq
                and c.divides_at == t
                and (c.dies_at is None or c.dies_at > t)
            ):
                self._divide(c)
            elif kind == "recheck" and seq == c.recheck_seq and c.born <= t < c.end:
                self.rechecks += 1  # a read this cell was waiting on has arrived (§7.5)
                self._resolve(c, born=False)
            elif kind == "cull" and c.dies_at is None:
                self._cull(c)
            elif kind.startswith("flow:"):
                self._flow(c, kind[5:])
            elif kind.startswith("move:"):
                self._step_move(c, kind[5:])
            elif kind == "die" and c.dies_at == t:
                # a death changes what its neighbours touch, exactly as a birth does (`_add`), so they
                # read their contacts again; the worm kills 110 cells, so a stale reading is not a corner
                around = list(self.neighbours(c)) if self._amount_signals and self.network is None else []
                if c.x is not None:
                    self.occupied.pop((c.x, c.y), None)  # the site is free again
                for other in around:
                    n = self.cells[other]
                    if n.born <= t < n.end:
                        self._resolve(n, born=False)
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

    def orders(self) -> list[dict]:
        """What a run did with each `order` along time (§7.6): the sequence it actually took, the
        inversions against the declared one, and the members that never happened. An order drives
        nothing, so this changes no result; it reports whether the program was right."""
        out = []
        for od in self.module.orders():
            if od.axis != "time":
                continue
            when = {m: self.cells[m].born for m in od.members if m in self.cells}
            taken = sorted(when, key=lambda m: (when[m], od.members.index(m)))
            # members sharing a step are not ordered against each other, and neither are members the
            # run put at the same instant: an inversion needs both a stated order and a real gap
            rank = {m: i for i, g in enumerate(od.groups) for m in g}
            inversions = sum(
                1
                for i, a in enumerate(taken)
                for b in taken[i + 1 :]
                if rank[a] > rank[b] and when[a] != when[b]
            )
            out.append(
                {
                    "order": od.id,
                    "expected": ["|".join(g) for g in od.groups],
                    "taken": taken,
                    "inversions": inversions,
                    "never happened": [m for m in od.members if m not in when],
                    "ok": inversions == 0 and len(when) == len(od.members),
                }
            )
        return out

    def check_asserts(self) -> list[dict]:
        return [evaluate_assert(self, a) for a in self.organism.asserts if not is_replicate_assert(a)]

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
            "fates_mode": self.fates,
            "network_steps": self.network_steps,
            "network_seed": self.network_seed if self.network is not None else None,
            "neighbours_from": "table" if self.contacts else ("grid" if self.spatial else "none"),
            "ambiguous_fates": sum(self.ambiguous.values()),
            "revised_fates": sum(self.revised.values()),
            "overruled_fates": sum(self.overruled.values()),
            "recheck": self.recheck,
            "rechecks": self.rechecks,
            "orders": self.orders(),
            "duplicate_decision_ids": dict(sorted(self.duplicate_ids.items())),
            "forced": {f: self.forced_cells[f] for f in sorted(self.forced)},
            "outside_competence": sum(self.outside_competence.values()),
            "refused_committed": sum(self.refused_committed.values()),
            "committed": dict(sorted(self.committed_cells.items())),
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


# ---- contacts, neighbours and per-cell networks (v0.4 §7.4) -----------------------------------


def read_contacts(path: str) -> list[tuple[float, dict[str, dict[str, float]]]]:
    """A time-resolved contact table: rows of `time cell cell [area]`, whitespace or tab separated.

    Returns snapshots sorted by time; the neighbours at t are the latest snapshot at or before t,
    each with the contact area (1.0 when the table does not say). A measured table replaces grid
    adjacency, which is why the worm can have neighbours without having a grid."""
    snapshots: dict[float, dict[str, dict[str, float]]] = {}
    for line in Path(path).read_text().splitlines():
        row = line.split("#", 1)[0].split()
        if len(row) < 3 or not row[0].replace(".", "", 1).replace("-", "", 1).isdigit():
            continue  # a header or a comment
        t, a, b = float(row[0]), row[1], row[2]
        area = float(row[3]) if len(row) > 3 else 1.0
        snap = snapshots.setdefault(t, {})
        snap.setdefault(a, {})[b] = area
        snap.setdefault(b, {})[a] = area
    return sorted(snapshots.items())


# ---- replicate runs (v0.4 §7.4): outcomes decided by noise are scored across seeds, never in one run ---

_EXACTLY = re.compile(
    r"^exactly one of\s+(.+?)\s+is\s+(\S+)(?:\s+at\s+([-0-9.]+)\s*(\w+)?)?\s+in\s*>=\s*([0-9.]+)\s*%$"
)
_SHARE = re.compile(r"^(\S+)\s+is\s+(\S+)(?:\s+at\s+([-0-9.]+)\s*(\w+)?)?\s+in\s+([0-9.]+)\.\.([0-9.]+)\s*%$")


def is_replicate_assert(text: str) -> bool:
    return bool(_EXACTLY.match(text.strip()) or _SHARE.match(text.strip()))


def _is(body: Body, name: str, fate: str, t: float | None) -> bool:
    c = body.cells.get(name)
    if c is None or (t is not None and not (c.born <= t)):
        return False
    return fate in (c.cell_type, c.terminal_name)


def replicate(module: Module, until: float, seeds: range | list[int], **kwargs) -> list[Body]:
    """The same program grown once per seed; every other argument is passed to the Body unchanged."""
    return [Body(module, seed=s, **kwargs).run(until=until) for s in seeds]


def evaluate_replicate_asserts(bodies: list[Body], asserts: list[str]) -> list[dict]:
    """`exactly one of A, B is T [at N min] in >= P%`: in at least P% of the runs exactly one of the named
    cells took the fate. `A is T [at N min] in lo..hi%`: the share of runs in which A took it lies in the
    range, which is how a coin-flip decision is tested: a runtime that picks the same winner every run has
    smuggled in an order, and fails it."""
    out = []
    n = len(bodies)
    for text in asserts:
        m = _EXACTLY.match(text.strip())
        if m:
            names = [x.strip() for x in re.split(r",|\bor\b", m.group(1)) if x.strip()]
            t = to_minutes(float(m.group(3)), m.group(4) or "min") if m.group(3) else None
            winners = Counter()
            hits = 0
            for b in bodies:
                took = [x for x in names if _is(b, x, m.group(2), t)]
                if len(took) == 1:
                    hits += 1
                    winners[took[0]] += 1
            share = 100.0 * hits / n if n else 0.0
            out.append(
                {
                    "assert": text,
                    "ok": n > 0 and share >= float(m.group(5)),
                    "value": round(share, 2),
                    "runs": n,
                    "winners": dict(winners),
                }
            )
            continue
        m = _SHARE.match(text.strip())
        if m:
            t = to_minutes(float(m.group(3)), m.group(4) or "min") if m.group(3) else None
            share = 100.0 * sum(_is(b, m.group(1), m.group(2), t) for b in bodies) / n if n else 0.0
            ok = n > 0 and float(m.group(5)) <= share <= float(m.group(6))
            out.append({"assert": text, "ok": ok, "value": round(share, 2), "runs": n})
            continue
        out.append({"assert": text, "ok": False, "error": "not a replicate assert"})
    return out
