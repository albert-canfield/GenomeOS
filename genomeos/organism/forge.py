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

from genomeos.certainty import Certainty
from genomeos.ir import Design, Evidence, EvidenceKind, Experiment, Module, to_minutes
from genomeos.runtime.body import Body, evaluate_assert

FORGE_EVIDENCE = Evidence(EvidenceKind.PREDICTED, "BioForge search over the organism program (not validated)")
_ASSERT = re.compile(
    r"^(count|deaths|fate|type|lineage)\s*(\S+)?\s+at\s+([-0-9.]+)\s*(\w+)?\s*(in|=|>=|<=|>|<)\s*(.+)$"
)
MAX_COMBINATIONS = 300

# Budget constraints (area H, registered 2026-09-27 before the code that honours them was written).
# The unit is the perturbation: each knockout, each addition and each continuous knob moved off its
# value in the program costs 1 unless the caller prices it (`Budget.costs`, e.g. bases or edits). The
# budget is a hard limit on the answer, not a loss term: the answer never spends more than the limit,
# and the result states what the spend bought and which evaluated candidates the limit excluded.
BUDGET_UNIT = "perturbations"  # knockouts + additions + knobs moved; per-name costs override 1 each
BUDGET_SPEC = (
    "A design run with a budget returns an answer whose spend is <= the limit, for every design and "
    "every limit; candidates over the limit are still run once at the program's own knob values so the "
    "result can name what the limit excluded and whether an excluded candidate had a lower loss than "
    "the answer. No budget keeps the previous behaviour byte for byte."
)
BUDGET_FALSIFIER = (
    "Any shipped design, any limit in 0..3: the reported answer spends more than the limit, or a "
    "candidate the limit excluded is missing from the excluded list, or an excluded candidate with "
    "a lower loss than the answer is not flagged."
)
# the known case, from an existing fixture: data/organisms/celegans/designs.bio, design
# two_intestinal_founders, cites Lin et al. 1995 ("without POP-1, MS takes the E fate"), and
# tests/test_design.py pins POP-1 as its unbudgeted answer at loss 0.
BUDGET_KNOWN_CASE = {
    "design": "data/organisms/celegans/designs.bio two_intestinal_founders",
    "source": "Lin et al. 1995, cited in the design block; answer pinned by tests/test_design.py",
    "limit 1": "answer -POP-1, solved, spend 1",
    "limit 0": "answer wild type, not solved, spend 0; -POP-1 listed as excluded with loss 0 < the answer's",
}

# Review item R4 (2026-09-28): the emitted experiment used to carry a fixed confidence of 0.3, a constant
# a reader takes for a probability of being right. No calibration record exists for any outcome of a
# design answer (docs/BIOFORGE-CONFIDENCE.md), so the answer now carries a `Certainty` whose probability
# is None with the reason, and the experiment's numeric confidence is left unstated: 0.0 is what the
# BioLang parser gives an `experiment` block with no `confidence:` key, so the emitted text and the
# returned Experiment agree. It is the absence of a stated value, not a probability of zero.
UNSTATED_CONFIDENCE = 0.0
EFFECT_UNIT = "normalised target distance removed (wild-type loss minus answer loss; 0 = no change)"
MODEL_SCORE_NAME = "loss: normalised distance from every target holding; 0 = all hold; lower ranks first"
UNCERTAINTY_NOTE = (
    "not measured: one simulation per candidate at a fixed seed; the search's spread is not a measurement"
)
PROBABILITY_UNAVAILABLE = (
    "no calibration record: no population of design answers with a published outcome and a negative case"
    " exists, and the 2-4 designs with a published answer were authored from that answer"
    " (docs/BIOFORGE-CONFIDENCE.md)"
)


def design_certainty(answer_loss: float, baseline_loss: float | None) -> Certainty:
    """The certainty record of a design answer. The effect is how much target distance the perturbation
    removed; nothing here turns it, or the loss, into a probability, so a larger effect cannot raise it."""
    effect = None if baseline_loss is None else baseline_loss - answer_loss
    return Certainty(
        evidence_category=FORGE_EVIDENCE.kind.value,
        effect_estimate=effect,
        effect_unit=EFFECT_UNIT,
        measurement_uncertainty=None,
        uncertainty_note=UNCERTAINTY_NOTE,
        model_score=answer_loss,
        model_score_name=MODEL_SCORE_NAME,
        probability=None,
        probability_unavailable=PROBABILITY_UNAVAILABLE,
    )


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
        # the two min(x, 0.3) below cap the IR's evidence-quality score of a value the search set: they
        # only lower it and never read `value`, so no knob position raises it. Not a probability (R4).
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
    baseline_loss: float | None = None  # loss of the unperturbed program at its own knob values

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
            confidence=UNSTATED_CONFIDENCE,  # was a fixed 0.3 until review item R4; see design_certainty
        )

    def certainty(self) -> Certainty | None:
        """Evidence, effect, uncertainty, score and probability of the answer, kept apart (review R4)."""
        b = self.best
        if b is None or not b.feasible:
            return None
        return design_certainty(b.loss, self.baseline_loss)

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
        parts.append(f'  evidence: predicted "{ex.evidence.source}"')
        cert = self.certainty()
        if cert is not None:
            parts += ["  " + line for line in cert.comment_lines()]
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
            "certainty": None if (c := self.certainty()) is None else c.to_dict(),
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


# The budget lives below DesignResult on purpose: tests/test_forge_calibration.py pins lines 93, 97
# and 170 of this file by number, so the budgeted types extend the classes above instead of editing them.
# (Those pins were retired on 2026-09-28 with the literal on line 170, review item R4; the property that
# replaces them is in tests/test_forge_certainty.py.)
@dataclass(slots=True)
class Budget:
    """A hard limit on what a design may spend, in `unit`; each perturbation costs 1 unless priced.

    Names priced in `costs` are factor names (a knockout or an addition of that factor) or knob paths
    (`timer:NAME.duration`); a knockout and an addition of the same factor cost the same."""

    limit: float
    costs: dict[str, float] = field(default_factory=dict)
    unit: str = BUDGET_UNIT

    def __post_init__(self) -> None:
        if self.limit < 0 or any(c < 0 for c in self.costs.values()):
            raise ValueError("a budget limit and its costs are non-negative")

    def cost(self, name: str) -> float:
        return float(self.costs.get(name, 1.0))

    def spend(self, knockouts: list[str], adds: list[str], moved: list[str]) -> float:
        return sum(self.cost(n) for n in [*knockouts, *adds, *moved])

    def to_dict(self) -> dict:
        return {"limit": self.limit, "unit": self.unit, "costs": dict(self.costs)}


@dataclass(slots=True)
class BudgetedCandidate(Candidate):
    spend: float = 0.0  # in the budget's unit
    within_budget: bool = True
    moved: list[str] = field(default_factory=list)  # knob paths off the program's own value


@dataclass(slots=True)
class BudgetedDesignResult(DesignResult):
    """A design run under a budget: the answer is the best candidate that spends no more than the limit."""

    budget: Budget = field(default_factory=lambda: Budget(0))

    @property
    def best(self) -> Candidate | None:
        b = self.candidates[0] if self.candidates else None
        return b if b is not None and b.within_budget else None  # never an answer over the budget

    def budget_report(self) -> dict:
        """What the spend bought, what the limit excluded, and whether the limit cost the target."""
        b = self.best
        excluded = [c for c in self.candidates if not c.within_budget]
        better = [
            c for c in excluded if c.feasible and (b is None or not b.feasible or c.loss < b.loss - 1e-12)
        ]
        return {
            **self.budget.to_dict(),
            "bought": None
            if b is None
            else {
                "perturbation": b.label(),
                "spend": b.spend,
                "loss": b.loss,
                "feasible": b.feasible,
                "solved": self.solved,
            },
            "excluded": [
                {"perturbation": c.label(), "spend": c.spend, "loss": c.loss, "feasible": c.feasible}
                for c in excluded
            ],
            "excluded_better": [c.label() for c in better],
            "limit_cost_the_target": bool(better),
        }

    def to_dict(self) -> dict:
        out = DesignResult.to_dict(self)
        for row, c in zip(out["candidates"], self.candidates, strict=True):
            row.update(spend=c.spend, within_budget=c.within_budget, moved=list(c.moved))
        out["budget"] = self.budget_report()
        return out

    def format(self, top: int = 6) -> str:
        out = [DesignResult.format(self, top)]
        rep = self.budget_report()
        bought = rep["bought"]
        out.append(
            f"  budget  {rep['limit']:g} {rep['unit']}: bought "
            + (f"{bought['perturbation']} for {bought['spend']:g}" if bought else "nothing")
        )
        ex = rep["excluded"]
        out.append(
            f"  budget  excluded {len(ex)} candidate(s) over the limit"
            + (": " + ", ".join(f"{e['perturbation']} ({e['spend']:g})" for e in ex[:6]) if ex else "")
        )
        if rep["excluded_better"]:
            out.append(
                "  budget  the limit cost the target: "
                + ", ".join(rep["excluded_better"][:6])
                + " would have come closer"
            )
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
    budget: Budget | None = None,
) -> DesignResult:
    """Search the design's perturbations. With `budget`, the answer never spends more than
    `budget.limit`; candidates over it are run once at the program's knob values and reported as
    excluded (`BudgetedDesignResult.budget_report`). Without one, behaviour is unchanged."""
    until = _horizon(module, design)
    knobs = [Knob.parse(spec) for spec in design.vary]
    rng = random.Random(0 if seed is None else seed)
    res = (
        DesignResult(design, until) if budget is None else BudgetedDesignResult(design, until, budget=budget)
    )
    original = {k.path: k.get(module) for k in knobs}

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
        args = (list(knockouts), list(adds), dict(values), loss, feasible, targets, keeps, body.summary())
        if budget is None:
            return Candidate(*args)
        moved = [p for p, v in values.items() if abs(v - original[p]) > 1e-12]
        spend = budget.spend(list(knockouts), list(adds), moved)
        return BudgetedCandidate(*args, spend=spend, within_budget=spend <= budget.limit + 1e-12, moved=moved)

    def better(a: Candidate, b: Candidate) -> bool:
        """a replaces b in the knob search; under a budget, never by going over it."""
        if budget is not None and not a.within_budget:
            return False
        return (a.feasible, -a.loss) > (b.feasible, -b.loss)

    def draw(k: Knob) -> float:
        return math.exp(rng.uniform(math.log(max(k.lo, 1e-9)), math.log(max(k.hi, 1e-9))))

    def restart_values(knockouts: list[str], adds: list[str]) -> dict[str, float]:
        left = math.inf if budget is None else budget.limit - budget.spend(knockouts, adds, [])
        if sum(1.0 if budget is None else budget.cost(k.path) for k in knobs) <= left + 1e-12:
            return {k.path: draw(k) for k in knobs}  # every knob affordable: the unbudgeted draw
        # a restart under a budget moves only the knobs the remaining budget can pay for
        values = dict(original)
        for k in rng.sample(knobs, len(knobs)):
            if budget.cost(k.path) <= left + 1e-12:
                values[k.path] = draw(k)
                left -= budget.cost(k.path)
        return values

    options = [("ko", f) for f in design.knockout_any_of] + [("add", f) for f in design.add_any_of]
    combos: list[tuple[list[str], list[str]]] = [([], [])]
    for n in range(1, design.at_most + 1):
        for combo in itertools.combinations(options, n):
            combos.append(([f for k, f in combo if k == "ko"], [f for k, f in combo if k == "add"]))
            if len(combos) >= MAX_COMBINATIONS:
                break
        if len(combos) >= MAX_COMBINATIONS:
            break
    for knockouts, adds in combos:
        best = evaluate(knockouts, adds, original)
        if not knockouts and not adds and res.baseline_loss is None:
            res.baseline_loss = best.loss  # the unperturbed program: the reference the effect is read against
        if budget is not None and not best.within_budget:
            res.candidates.append(best)  # over the limit before any knob moves: run once, report as excluded
            continue
        affordable = budget is None or any(
            budget.cost(k.path) <= budget.limit - best.spend + 1e-12 for k in knobs
        )  # no knob the remaining budget can pay for: every move would be refused, so none is run
        if knobs and affordable:  # log-space local search on the knobs, as in the network BioForge
            for r in range(restarts):
                current = best if r == 0 else evaluate(knockouts, adds, restart_values(knockouts, adds))
                step = 0.5
                for _ in range(iterations):
                    k = rng.choice(knobs)
                    values = dict(current.values)
                    values[k.path] = values[k.path] * math.exp(rng.gauss(0, step))
                    cand = evaluate(knockouts, adds, values)
                    if better(cand, current):
                        current = cand
                    else:
                        step = max(0.05, step * 0.97)
                if better(current, best):
                    best = current
        res.candidates.append(best)
    if budget is None:
        res.candidates.sort(key=lambda c: (not c.feasible, c.loss, c.perturbations, c.label()))
    else:
        res.candidates.sort(
            key=lambda c: (not c.within_budget, not c.feasible, c.loss, c.perturbations, c.label())
        )
    return res


def run_designs(
    module: Module, names: list[str] | None = None, seed: int | None = None, budget: Budget | None = None
) -> list[DesignResult]:
    return [
        run_design(module, d, seed=seed, budget=budget)
        for d in module.designs
        if not names or d.name in names
    ]
