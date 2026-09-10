"""Uncertainty report per biological level (task 3.5).

Every simulation output carries a report saying how grounded it is at each
level. Confidence is aggregated from the evidence attached to the rules,
parameters and events that actually fired, weighted by evidence kind:

    experimental 1.0 · curated 0.9 · predicted 0.6 · inferred 0.4 · none 0.0

A level with nothing behind it is reported as UNKNOWN rather than 0, so that
"no information" is distinguishable from "measured to be unreliable".
"""

from __future__ import annotations

from dataclasses import dataclass, field

from genomeos.ir import Event, EvidenceKind, Module, Parameter, Rule

LEVELS = ("molecular", "cellular", "tissue", "organism")
KIND_WEIGHT = {
    EvidenceKind.EXPERIMENTAL: 1.0,
    EvidenceKind.CURATED: 0.9,
    EvidenceKind.PREDICTED: 0.6,
    EvidenceKind.INFERRED: 0.4,
    EvidenceKind.NONE: 0.0,
}


@dataclass(slots=True)
class LevelReport:
    level: str
    items: int = 0
    confidence: float | None = None  # None = UNKNOWN
    weakest: str = ""

    @property
    def label(self) -> str:
        if self.confidence is None:
            return "UNKNOWN"
        return "high" if self.confidence >= 0.75 else "medium" if self.confidence >= 0.45 else "low"


@dataclass(slots=True)
class UncertaintyReport:
    levels: dict[str, LevelReport] = field(default_factory=lambda: {lv: LevelReport(lv) for lv in LEVELS})

    def add(self, level: str, confidence: float, kind: EvidenceKind, name: str) -> None:
        lr = self.levels[level]
        score = confidence * KIND_WEIGHT[kind]
        if lr.confidence is None:
            lr.confidence = score
            lr.weakest = name
        else:
            lr.confidence = (lr.confidence * lr.items + score) / (lr.items + 1)
            if score < lr.confidence and (not lr.weakest or score < self._score_of(lr)):
                lr.weakest = name
        lr.items += 1

    @staticmethod
    def _score_of(lr: LevelReport) -> float:
        return lr.confidence if lr.confidence is not None else 0.0

    def add_rule(self, rule: Rule, level: str = "molecular") -> None:
        self.add(level, rule.confidence, rule.evidence.kind, rule.id)

    def add_parameter(self, p: Parameter, level: str = "cellular") -> None:
        self.add(level, p.confidence, p.evidence.kind, p.name)

    def add_event(self, e: Event, level: str = "cellular") -> None:
        self.add(level, e.confidence, e.evidence.kind, e.id)

    def to_dict(self) -> dict:
        return {
            lv: {
                "items": r.items,
                "confidence": None if r.confidence is None else round(r.confidence, 3),
                "label": r.label,
                "weakest": r.weakest,
            }
            for lv, r in self.levels.items()
        }

    def format(self) -> str:
        out = ["uncertainty by level:"]
        for lv, r in self.levels.items():
            if r.confidence is None:
                out.append(f"  {lv:<10} UNKNOWN   (nothing simulated at this level)")
            else:
                bar = "█" * int(r.confidence * 20)
                weakest = f"  weakest: {r.weakest}" if r.weakest else ""
                out.append(f"  {lv:<10} {bar:<20} {r.confidence:.2f} {r.label:<6} {r.items} items{weakest}")
        return "\n".join(out)


def report_for_network(module: Module, active_rules: list[Rule]) -> UncertaintyReport:
    rep = UncertaintyReport()
    for r in active_rules:
        rep.add_rule(r, "molecular")
    for p in module.parameters.values():
        rep.add_parameter(p, "molecular")
    for e in module.events:
        rep.add_event(e, "cellular")
    return rep


def report_for_ageing(parameters: list[Parameter]) -> UncertaintyReport:
    rep = UncertaintyReport()
    for p in parameters:
        rep.add_parameter(p, "cellular")
    # a population of one cell type says something about tissue, at reduced confidence
    for p in parameters:
        rep.add(  # tissue-level claims inherit cell-level evidence discounted
            "tissue", p.confidence * 0.6, p.evidence.kind, p.name
        )
    return rep
