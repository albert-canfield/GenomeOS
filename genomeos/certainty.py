# SPDX-License-Identifier: AGPL-3.0-or-later
"""Certainty as five separate quantities, and a probability that cannot exist without its record.

Review item R4 (2026-09-28): a single `confidence` number had been carrying the kind of evidence, the
size of an effect, the spread of a measurement, a model's score and a probability of being right, and
could be derived from effect magnitude or set as a constant. None of those establishes a probability
of correctness. This module keeps them apart:

- `evidence_category`: what kind of support a claim has (experimental, curated, predicted, ...);
- `effect_estimate` with `effect_unit`: how large the effect is, in a stated unit;
- `measurement_uncertainty`: the spread of that estimate in the same unit, or None with the reason;
- `model_score` with `model_score_name`: whatever the producing model ranks by, named, not a probability;
- `probability`: None unless calibrated. A `Probability` can only be built with a `Calibration` that
  names the outcome it is a probability of, the population it was calibrated on and the evaluation
  method, so every probability identifies all three.

An unavailable probability stays unavailable: `Certainty` refuses `probability=None` without a stated
`probability_unavailable` reason, and nothing here converts an effect, a score or a category into one.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True, slots=True)
class Calibration:
    """The record a probability is quoted against: of what, on whom, by which method."""

    outcome: str  # the event the probability is a probability of
    population: str  # the population it was calibrated on (and so the one it may be quoted for)
    method: str  # how calibration was evaluated
    n: int | None = None  # size of that population, when known

    def __post_init__(self) -> None:
        for name in ("outcome", "population", "method"):
            if not str(getattr(self, name)).strip():
                raise ValueError(f"a calibration record must name its {name}")


@dataclass(frozen=True, slots=True)
class Probability:
    """A probability of one named outcome, calibrated on one named population by one named method."""

    value: float
    calibration: Calibration

    def __post_init__(self) -> None:
        if not 0.0 <= self.value <= 1.0:
            raise ValueError(f"a probability lies in 0..1, not {self.value}")
        if not isinstance(self.calibration, Calibration):
            raise ValueError("a probability needs a Calibration record")

    def to_dict(self) -> dict[str, Any]:
        c = self.calibration
        return {
            "value": self.value,
            "outcome": c.outcome,
            "calibration_population": c.population,
            "method": c.method,
            "n": c.n,
        }


@dataclass(frozen=True, slots=True)
class Certainty:
    """Evidence category, effect, measurement uncertainty, model score and probability, kept apart."""

    evidence_category: str
    effect_estimate: float | None = None
    effect_unit: str = ""
    measurement_uncertainty: float | None = None  # same unit as the effect
    uncertainty_note: str = ""  # what the uncertainty is, or why there is none
    model_score: float | None = None
    model_score_name: str = ""
    probability: Probability | None = None
    probability_unavailable: str = ""  # required when probability is None

    def __post_init__(self) -> None:
        if not self.evidence_category.strip():
            raise ValueError("a certainty record names its evidence category")
        if self.probability is None and not self.probability_unavailable.strip():
            raise ValueError("an unavailable probability must say why it is unavailable")
        if self.probability is not None and self.probability_unavailable.strip():
            raise ValueError("a probability is either present or unavailable, not both")
        if self.effect_estimate is not None and not self.effect_unit.strip():
            raise ValueError("an effect estimate carries its unit")
        if self.model_score is not None and not self.model_score_name.strip():
            raise ValueError("a model score names what it scores")

    def to_dict(self) -> dict[str, Any]:
        return {
            "evidence_category": self.evidence_category,
            "effect_estimate": self.effect_estimate,
            "effect_unit": self.effect_unit,
            "measurement_uncertainty": self.measurement_uncertainty,
            "uncertainty_note": self.uncertainty_note,
            "model_score": self.model_score,
            "model_score_name": self.model_score_name,
            "probability": None if self.probability is None else self.probability.to_dict(),
            "probability_unavailable": self.probability_unavailable or None,
        }

    def comment_lines(self) -> list[str]:
        """The record as BioLang comments: the grammar has no key for it, and a comment invents nothing."""

        def num(x: float | None) -> str:
            return "unavailable" if x is None else f"{x:.4g}"

        out = [
            f"# certainty: evidence {self.evidence_category}",
            f"#   effect: {num(self.effect_estimate)} {self.effect_unit}".rstrip(),
            f"#   measurement uncertainty: {num(self.measurement_uncertainty)}"
            + (f" ({self.uncertainty_note})" if self.uncertainty_note else ""),
            f"#   model score: {num(self.model_score)}"
            + (f" ({self.model_score_name})" if self.model_score_name else ""),
        ]
        if self.probability is None:
            out.append(f"#   probability: unavailable ({self.probability_unavailable})")
        else:
            p = self.probability.to_dict()
            out.append(
                f"#   probability: {p['value']:.4g} of {p['outcome']}; calibrated on"
                f" {p['calibration_population']} by {p['method']}"
            )
        return out
