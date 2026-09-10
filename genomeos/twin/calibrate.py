"""Calibration: turn a measurement into an engine parameter (rules from data).

`calibrate_attrition` fits the one low-confidence knob of the ageing model,
divisions_per_year for a cell type, so that the simulated mean telomere at
the subject's age matches a measured value. The fitted parameter is stored as
a result with `inferred` evidence naming the measurement, and the twin uses it.
"""

from __future__ import annotations

from dataclasses import dataclass

from genomeos.ir import Evidence, EvidenceKind, Parameter
from genomeos.results import save_result
from genomeos.runtime.cell import CELL_TYPES, GLOBAL_PARAMS


@dataclass(slots=True)
class Calibration:
    cell_type: str
    age_years: float
    measured_bp: float
    fitted_divisions_per_year: float
    predicted_bp: float
    net_bp_per_year: float
    parameter: Parameter


def expected_telomere(
    cell_type: str, age_years: float, divisions_per_year: float | None = None, birth_bp: float | None = None
) -> float:
    """Closed-form expectation of the ageing model's mean telomere (no senescence, no noise)."""
    p = CELL_TYPES[cell_type]
    div = p.divisions_per_year.value if divisions_per_year is None else divisions_per_year
    keep = 1.0 - p.telomerase_compensation.value
    start = GLOBAL_PARAMS["telomere_at_birth_bp"].value if birth_bp is None else birth_bp
    return start - div * p.telomere_loss_per_division_bp.value * keep * age_years


def calibrate_attrition(
    cell_type: str, age_years: float, measured_bp: float, source: str, save: bool = True
) -> Calibration:
    """Solve divisions_per_year so that expected_telomere(age) == measured_bp."""
    p = CELL_TYPES[cell_type]
    keep = 1.0 - p.telomerase_compensation.value
    start = GLOBAL_PARAMS["telomere_at_birth_bp"].value
    loss = p.telomere_loss_per_division_bp.value
    if age_years <= 0 or keep <= 0 or loss <= 0:
        raise ValueError(
            "cannot calibrate: age, telomerase compensation and loss per division must allow attrition"
        )
    div = max(0.0, (start - measured_bp) / (loss * keep * age_years))
    predicted = expected_telomere(cell_type, age_years, div)
    param = Parameter(
        "divisions_per_year",
        round(div, 3),
        "1/yr",
        Evidence(
            EvidenceKind.INFERRED,
            f"calibrated to {source}",
            note=(
                f"{measured_bp:.0f} bp at {age_years:.0f} y; birth {start:.0f} bp; "
                f"{loss:.0f} bp/division × {keep:.2f} kept"
            ),
        ),
        0.4,
    )
    cal = Calibration(
        cell_type, age_years, measured_bp, div, predicted, (start - measured_bp) / age_years, param
    )
    if save:
        save_result(
            f"calibration_{cell_type}_attrition",
            {
                "cell_type": cell_type,
                "age_years": age_years,
                "measured_bp": measured_bp,
                "source": source,
                "fitted_divisions_per_year": round(div, 4),
                "net_bp_per_year": round(cal.net_bp_per_year, 2),
                "evidence": {
                    "kind": "inferred",
                    "source": param.evidence.source,
                    "note": param.evidence.note,
                },
                "confidence": 0.4,
            },
        )
    return cal
