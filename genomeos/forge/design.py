"""BioForge (Phase 5): design under constraints, in silico.

Given a BioLang module, a set of tunable knobs (gene max/basal rates, rule
strength/threshold/hill, module parameters), an objective on the simulated
trajectory and hard constraints, search for a design that meets the target.
The search is a simple random-restart, log-space local search: dependency
free and good enough for tens of knobs. Every changed value is emitted with
`predicted` evidence ("BioForge search") and low confidence: a design is a
hypothesis until an experiment confirms it.
"""

from __future__ import annotations

import copy
import math
import random
from collections.abc import Callable
from dataclasses import dataclass, field

from genomeos.ir import Evidence, EvidenceKind, Module, Parameter, Rule
from genomeos.runtime.grn import NetworkRuntime, Trajectory

Objective = Callable[[Trajectory], float]  # lower is better (a loss)
Constraint = Callable[[Trajectory], bool]
FORGE_EVIDENCE = Evidence(EvidenceKind.PREDICTED, "BioForge search (in silico design, not validated)")


@dataclass(slots=True)
class Knob:
    path: str  # "gene:tetR.max", "rule:<id>.strength", "param:translation_rate"
    lo: float
    hi: float

    def _parts(self) -> tuple[str, str, str]:
        kind, _, rest = self.path.partition(":")
        if kind == "param":
            return kind, rest, ""
        name, _, attr = rest.rpartition(".")
        return kind, name, attr

    def get(self, m: Module) -> float:
        kind, name, attr = self._parts()
        if kind == "gene":
            g = m.entities[name]
            return g.basal_rate if attr == "basal" else float(g.attrs.get("max_rate", 0.0))
        if kind == "rule":
            return float(getattr(next(r for r in m.rules if r.id == name), attr))
        return m.parameters[name].value

    def set(self, m: Module, value: float) -> None:
        kind, name, attr = self._parts()
        value = min(self.hi, max(self.lo, value))
        if kind == "gene":
            g = m.entities[name]
            if attr == "basal":
                g.basal_rate = value
            else:
                g.attrs["max_rate"] = value
            g.evidence, g.confidence = FORGE_EVIDENCE, min(g.confidence, 0.3)
        elif kind == "rule":
            r: Rule = next(r for r in m.rules if r.id == name)
            setattr(r, attr, value)
            r.evidence, r.confidence = FORGE_EVIDENCE, min(r.confidence, 0.3)
        else:
            p = m.parameters[name]
            m.parameters[name] = Parameter(p.name, value, p.unit, FORGE_EVIDENCE, min(p.confidence, 0.3))


@dataclass(slots=True)
class Design:
    module: Module
    values: dict[str, float]
    loss: float
    feasible: bool


@dataclass(slots=True)
class DesignResult:
    best: Design
    original: dict[str, float]
    evaluations: int
    history: list[float] = field(default_factory=list)

    def changes(self) -> list[tuple[str, float, float]]:
        return [
            (k, self.original[k], v) for k, v in self.best.values.items() if abs(v - self.original[k]) > 1e-12
        ]


def period_of(traj: Trajectory, species: str) -> float | None:
    """Mean interval between successive peaks; None if fewer than three peaks."""
    xs = traj.levels[species]
    peaks = [traj.times[i] for i in range(1, len(xs) - 1) if xs[i - 1] < xs[i] >= xs[i + 1] and xs[i] > 1e-6]
    if len(peaks) < 3:
        return None
    gaps = [b - a for a, b in zip(peaks[1:], peaks[2:], strict=False)]  # skip the transient first gap
    return sum(gaps) / len(gaps) if gaps else None


def forge(
    module: Module,
    knobs: list[Knob],
    objective: Objective,
    constraints: list[Constraint] | None = None,
    hours: float = 60.0,
    dt: float = 0.02,
    initial: dict[str, float] | None = None,
    iterations: int = 60,
    restarts: int = 3,
    seed: int = 0,
) -> DesignResult:
    rng = random.Random(seed)
    constraints = constraints or []
    original = {k.path: k.get(module) for k in knobs}
    history: list[float] = []
    evaluations = 0

    def simulate(m: Module) -> Trajectory:
        return NetworkRuntime(m).run(hours=hours, dt=dt, initial=initial, record_every=max(1, int(0.1 / dt)))

    def evaluate(values: dict[str, float]) -> Design:
        nonlocal evaluations
        m = copy.deepcopy(module)
        for k in knobs:
            k.set(m, values[k.path])
        traj = simulate(m)
        evaluations += 1
        feasible = all(c(traj) for c in constraints)
        loss = objective(traj) if feasible else math.inf
        return Design(m, dict(values), loss, feasible)

    best = evaluate(original)
    history.append(best.loss)
    for r in range(restarts):
        current = (
            best
            if r == 0
            else evaluate(
                {
                    k.path: math.exp(rng.uniform(math.log(k.lo), math.log(k.hi)))
                    if k.lo > 0
                    else rng.uniform(k.lo, k.hi)
                    for k in knobs
                }
            )
        )
        step = 0.5
        for _ in range(iterations):
            k = rng.choice(knobs)
            values = dict(current.values)
            v = values[k.path]
            values[k.path] = v * math.exp(rng.gauss(0, step)) if v > 0 else v + rng.gauss(0, step)
            cand = evaluate(values)
            if cand.loss < current.loss:
                current = cand
            else:
                step = max(0.05, step * 0.97)
            history.append(current.loss)
        if current.loss < best.loss:
            best = current
    return DesignResult(best, original, evaluations, history)
