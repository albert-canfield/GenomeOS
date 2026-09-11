"""Evidence: what is known, how it is known, and how sure we are.

Every conclusion in the therapeutics pipeline points back to an Evidence
record. The record keeps the source, what kind of source it is, the claim it
supports, the strength of the observation behind it, and a confidence. A
curated database entry, a clinical trial, a patient measurement and a
computational guess are never mixed: `level` separates them and the report
groups by it.

`Missing` is a first-class value. Where a fact cannot be established the
pipeline records why, not a default.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any, Literal

SourceType = Literal[
    "database",
    "publication",
    "clinical_trial",
    "patient_data",
    "prediction",
    "derived",
]

#: Strength of the observation behind a claim, strongest first.
EvidenceLevel = Literal[
    "clinical",  # observed in treated patients
    "human",  # observed in human tissue or human subjects
    "preclinical",  # animal or disease-model observation
    "in_vitro",  # cell or biochemical observation
    "computational",  # produced by a predictor
    "inferred",  # produced by this pipeline's own reasoning
]

LEVEL_ORDER: tuple[str, ...] = (
    "clinical",
    "human",
    "preclinical",
    "in_vitro",
    "computational",
    "inferred",
)

#: How much a level is allowed to contribute to an evidence-strength score.
LEVEL_WEIGHT: dict[str, float] = {
    "clinical": 1.0,
    "human": 0.85,
    "preclinical": 0.6,
    "in_vitro": 0.5,
    "computational": 0.3,
    "inferred": 0.2,
}


@dataclass(slots=True)
class Evidence:
    """One traceable statement and where it came from."""

    source: str
    source_type: SourceType
    claim: str
    level: EvidenceLevel
    confidence: float
    identifier: str | None = None
    url: str | None = None
    retrieved_at: str | None = None
    version: str | None = None
    patient_derived: bool = False

    def to_dict(self) -> dict[str, Any]:
        d: dict[str, Any] = {
            "source": self.source,
            "source_type": self.source_type,
            "claim": self.claim,
            "level": self.level,
            "confidence": round(self.confidence, 3),
            "patient_derived": self.patient_derived,
        }
        for k in ("identifier", "url", "retrieved_at", "version"):
            v = getattr(self, k)
            if v:
                d[k] = v
        return d


def today() -> str:
    return time.strftime("%Y-%m-%d")


def database(
    source: str,
    claim: str,
    confidence: float,
    level: EvidenceLevel = "human",
    identifier: str | None = None,
    url: str | None = None,
    version: str | None = None,
) -> Evidence:
    """A curated database record."""
    return Evidence(source, "database", claim, level, confidence, identifier, url, today(), version)


def derived(source: str, claim: str, confidence: float) -> Evidence:
    """A conclusion this pipeline drew from other evidence."""
    return Evidence(source, "derived", claim, "inferred", confidence, retrieved_at=today())


def prediction(source: str, claim: str, confidence: float) -> Evidence:
    """Output of a predictor. Never a measurement."""
    return Evidence(source, "prediction", claim, "computational", confidence, retrieved_at=today())


def patient(source: str, claim: str, confidence: float, level: EvidenceLevel = "human") -> Evidence:
    """A measurement from this patient's own sample."""
    return Evidence(
        source, "patient_data", claim, level, confidence, retrieved_at=today(), patient_derived=True
    )


def clinical(
    source: str, claim: str, confidence: float, identifier: str | None = None, url: str | None = None
) -> Evidence:
    """An observation from treated patients (approved therapy, trial)."""
    return Evidence(source, "clinical_trial", claim, "clinical", confidence, identifier, url, today())


@dataclass(slots=True)
class Missing:
    """A fact that could not be established, and what would establish it."""

    what: str
    reason: str
    would_be_resolved_by: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "what": self.what,
            "reason": self.reason,
            "would_be_resolved_by": list(self.would_be_resolved_by),
        }


@dataclass(slots=True)
class EvidenceLedger:
    """The evidence collected for one candidate, with the missing facts beside it."""

    items: list[Evidence] = field(default_factory=list)
    missing: list[Missing] = field(default_factory=list)

    def add(self, *evidence: Evidence) -> None:
        self.items.extend(e for e in evidence if e is not None)

    def lack(self, what: str, reason: str, resolved_by: tuple[str, ...] = ()) -> Missing:
        m = Missing(what, reason, resolved_by)
        self.missing.append(m)
        return m

    def by_level(self) -> dict[str, list[Evidence]]:
        out: dict[str, list[Evidence]] = {}
        for e in self.items:
            out.setdefault(e.level, []).append(e)
        return {k: out[k] for k in LEVEL_ORDER if k in out}

    def strength(self) -> float:
        """0-1: how strongly the best evidence supports this candidate at all."""
        if not self.items:
            return 0.0
        return max(LEVEL_WEIGHT.get(e.level, 0.2) * e.confidence for e in self.items)

    def to_dict(self) -> dict[str, Any]:
        return {
            "by_level": {k: [e.to_dict() for e in v] for k, v in self.by_level().items()},
            "strength": round(self.strength(), 3),
            "missing": [m.to_dict() for m in self.missing],
        }


def strongest_level(items: list[Evidence]) -> str:
    """The strongest evidence level present, or 'none'."""
    for level in LEVEL_ORDER:
        if any(e.level == level for e in items):
            return level
    return "none"
