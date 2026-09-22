# SPDX-License-Identifier: AGPL-3.0-or-later
"""BioForge over organisms: a `design` block asks which perturbation reaches a goal.

A design names the perturbations BioForge may try (factors to knock out or add, up
to `at_most` at once, and continuous knobs on timers and decisions), the targets to
reach (assert grammar; the loss is the normalised distance from holding) and the
constraints that must keep holding. Every candidate is an experiment run against
the same program, seed and horizon; the answer is the feasible candidate with the
smallest loss and the fewest perturbations, emitted as an `experiment` block with
`predicted` evidence: a design is a hypothesis until the bench confirms it.
"""

from __future__ import annotations

import copy
import itertools
import math
import random
import re
from dataclasses import dataclass, field

from genomeos.ir import Design, Evidence, EvidenceKind, Experiment, Module, to_minutes
from genomeos.runtime.body import Body, evaluate_assert

FORGE_EVIDENCE = Evidence(EvidenceKind.PREDICTED, "BioForge search over the organism program (not validated)")
_ASSERT = re.compile(
    r"^(count|deaths|fate|type|lineage)\s*(\S+)?\s+at\s+([-0-9.]+)\s*(\w+)?\s*(in|=|>=|<=|>|<)\s*(.+)$"
)
MAX_COMBINATIONS = 300


@dataclass(slots=True)
class Knob:
    kind: str  # timer | decision
    name: str
    attr: str  # duration | fraction | after
    lo: float
    hi: float

    @classmethod
    def parse(cls, spec: str) -> Knob:
        kind, name, attr, rng = spec.split()
        lo, _, hi = rng.partition("..")
        return cls(kind, name, attr, float(lo), float(hi))

    @property
    def path(self) -> str:
        return f"{self.kind}:{self.name}.{self.attr}"

    def get(self, m: Module) -> float:
        if self.kind == "timer":
            t = m.timer(self.name)
            if t is None:
                raise ValueError(f"design varies unknown timer {self.name!r}")
            return t.duration
        d = next((x for x in m.decisions if x.id == self.name), None)
        if d is None:
            raise ValueError(f"design varies unknown decision {self.name!r}")
        return float(getattr(d, self.attr) or 0.0)

    def set(self, m: Module, value: float) -> None:
        value = min(self.hi, max(self.lo, value))
        if self.kind == "timer":
            t = m.timer(self.name)
            t.duration = value
            t.evidence, t.confidence = FORGE_EVIDENCE, min(t.confidence, 0.3)
        else:
            d = next(x for x in m.decisions if x.id == self.name)
            setattr(d, self.attr, value)
            d.evidence, d.confidence = FORGE_EVIDENCE, min(d.confidence, 0.3)


def distance(check: dict, text: str) -> float:
    """How far an assert is from holding, normalised: 0 when it holds."""
    if check.get("ok"):
        return 0.0
    m = _ASSERT.match(text)
    if not m or "value" not in check:
        return 10.0
    v, op, rhs = float(check["value"]), m.group(5), m.group(6)
    if op == "in":
        lo, hi = (float(x) for x in rhs.split(".."))
        span = max(hi - lo, 1.0)
        return (lo - v) / span if v < lo else (v - hi) / span
    n = float(rhs)
    scale = max(abs(n), 1.0)
    if op in ("=",):
        return abs(v - n) / scale
    if op in (">=", ">"):
        return max(0.0, n - v) / scale + (1e-3 if op == ">" else 0.0)
    return max(0.0, v - n) / scale + (1e-3 if op == "<" else 0.0)


@dataclass(slots=True)
class Candidate:
    knockouts: list[str]
    adds: list[str]
    values: dict[str, float]
    loss: float
    feasible: bool
    targets: list[dict]
    keeps: list[dict]
    summary: dict

    @property
    def perturbations(self) -> int:
        return len(self.knockouts) + len(self.adds) + len(self.values)

    def label(self) -> str:
        parts = ["-" + k for k in self.knockouts] + ["+" + a for a in self.adds]
        parts += [f"{k}={v:.4g}" for k, v in self.values.items()]
        return ", ".join(parts) or "wild type"


@dataclass(slots=True)
class DesignResult:
    design: Design
    until: float
    candidates: list[Candidate] = field(default_factory=list)
    evaluations: int = 0

    @property
    def best(self) -> Candidate | None:
        return self.candidates[0] if self.candidates else None

    @property
    def solved(self) -> bool:
        b = self.best
        return b is not None and b.feasible and b.loss < 1e-9

    def to_experiment(self) -> Experiment | None:
        b = self.best
        if b is None or not b.feasible:
            return None
        return Experiment(
            name=f"design_{self.design.name}",
            knockouts=list(b.knockouts),
            adds=list(b.adds),
            until=self.until,
            asserts=list(self.design.targets) + list(self.design.keeps),
            expect=f"BioForge: {b.label()} reaches {'; '.join(self.design.targets)}",
            evidence=FORGE_EVIDENCE,
            confidence=0.3,
        )

    def to_bio(self) -> str:
        ex = self.to_experiment()
        if ex is None:
            return f"# design {self.design.name}: no feasible candidate found\n"
        parts = [f"experiment {ex.name} {{"]
        if ex.knockouts:
            parts.append(f"  knockout: {', '.join(ex.knockouts)}")
        if ex.adds:
            parts.append(f"  add: {', '.join(ex.adds)}")
        if self.best and self.best.values:
            parts.append("  # knobs: " + ", ".join(f"{k} = {v:.4g}" for k, v in self.best.values.items()))
        parts.append(f"  until: {self.until:g} min")
        parts.append(f'  expect: "{ex.expect}"')
        parts += [f"  assert: {a}" for a in ex.asserts]
        parts.append(f'  evidence: predicted "{ex.evidence.source}"; confidence: {ex.confidence}')
        parts.append("}")
        return "\n".join(parts) + "\n"

    def to_dict(self) -> dict:
        return {
            "design": self.design.name,
            "until_min": self.until,
            "solved": self.solved,
            "evaluations": self.evaluations,
            "candidates": [
                {
                    "perturbation": c.label(),
                    "knockouts": c.knockouts,
                    "adds": c.adds,
                    "values": c.values,
                    "loss": c.loss,
                    "feasible": c.feasible,
                    "targets": c.targets,
                    "keeps": c.keeps,
                }
                for c in self.candidates
            ],
            "experiment": self.to_bio(),
        }

    def format(self, top: int = 6) -> str:
        d = self.design
        out = [
            f"design {d.name}: {self.evaluations} organism runs to {self.until:g} min; "
            + (
                "solved"
                if self.solved
                else "best effort"
                if self.best and self.best.feasible
                else "no feasible candidate"
            )
        ]
        for t in d.targets:
            out.append(f"  target  {t}")
        for k in d.keeps:
            out.append(f"  keep    {k}")
        for c in self.candidates[:top]:
            got = ", ".join(f"{a.get('value', a.get('error'))!s:.10}" for a in c.targets)
            flag = "ok  " if c.feasible and c.loss < 1e-9 else "feas" if c.feasible else "FAIL"
            out.append(f"  {flag} {c.label():<40} loss {c.loss:8.3g}  targets got {got}")
        if len(self.candidates) > top:
            out.append(f"  ... {len(self.candidates) - top} more candidates")
        if self.best and self.best.feasible:
            out.append("  answer as an experiment block (predicted):")
            out.extend("    " + line for line in self.to_bio().rstrip().splitlines())
        return "\n".join(out)


def _horizon(module: Module, design: Design) -> float:
    if design.until is not None:
        return design.until
    if module.stages:
        return max(to_minutes(s.start, s.unit) for s in module.stages)
    return 800.0


def run_design(
    module: Module,
    design: Design,
    seed: int | None = None,
    means: bool = True,
    iterations: int = 30,
    restarts: int = 1,
) -> DesignResult:
    until = _horizon(module, design)
    knobs = [Knob.parse(spec) for spec in design.vary]
    rng = random.Random(0 if seed is None else seed)
    res = DesignResult(design, until)

    def evaluate(knockouts: list[str], adds: list[str], values: dict[str, float]) -> Candidate:
        m = module
        if knobs:
            m = copy.deepcopy(module)
            for k in knobs:
                k.set(m, values[k.path])
        body = Body(m, seed=seed, means=means, knockouts=set(knockouts), adds=set(adds)).run(until=until)
        res.evaluations += 1
        targets = [evaluate_assert(body, t) for t in design.targets]
        keeps = [evaluate_assert(body, k) for k in design.keeps]
        loss = sum(distance(c, t) for c, t in zip(targets, design.targets, strict=True))
        feasible = all(c["ok"] for c in keeps)
        return Candidate(
            list(knockouts), list(adds), dict(values), loss, feasible, targets, keeps, body.summary()
        )

    options = [("ko", f) for f in design.knockout_any_of] + [("add", f) for f in design.add_any_of]
    combos: list[tuple[list[str], list[str]]] = [([], [])]
    for n in range(1, design.at_most + 1):
        for combo in itertools.combinations(options, n):
            combos.append(([f for k, f in combo if k == "ko"], [f for k, f in combo if k == "add"]))
            if len(combos) >= MAX_COMBINATIONS:
                break
        if len(combos) >= MAX_COMBINATIONS:
            break
    original = {k.path: k.get(module) for k in knobs}
    for knockouts, adds in combos:
        best = evaluate(knockouts, adds, original)
        if knobs:  # log-space local search on the knobs, as in the network BioForge
            for r in range(restarts):
                current = (
                    best
                    if r == 0
                    else evaluate(
                        knockouts,
                        adds,
                        {
                            k.path: math.exp(
                                rng.uniform(math.log(max(k.lo, 1e-9)), math.log(max(k.hi, 1e-9)))
                            )
                            for k in knobs
                        },
                    )
                )
                step = 0.5
                for _ in range(iterations):
                    k = rng.choice(knobs)
                    values = dict(current.values)
                    values[k.path] = values[k.path] * math.exp(rng.gauss(0, step))
                    cand = evaluate(knockouts, adds, values)
                    if (cand.feasible, -cand.loss) > (current.feasible, -current.loss):
                        current = cand
                    else:
                        step = max(0.05, step * 0.97)
                if (current.feasible, -current.loss) > (best.feasible, -best.loss):
                    best = current
        res.candidates.append(best)
    res.candidates.sort(key=lambda c: (not c.feasible, c.loss, c.perturbations, c.label()))
    return res


def run_designs(
    module: Module, names: list[str] | None = None, seed: int | None = None
) -> list[DesignResult]:
    return [run_design(module, d, seed=seed) for d in module.designs if not names or d.name in names]
