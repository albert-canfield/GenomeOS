"""Accessibility is not specificity.

A protein can be beautifully exposed on the tumour cell and still be a bad
target, because the same protein sits on cardiomyocytes. This stage keeps the
two questions apart: can a binder reach it, and does reaching it hit anything
that matters.

Three rules hold throughout:

  * A DNA mutation never implies that the protein is made. Without tumour
    RNA-seq or proteomics, tumour expression is unavailable, and says so.
  * Population healthy-tissue levels (Human Protein Atlas consensus nTPM) are
    labelled population-level and never presented as this patient's tumour.
  * Patient data, when supplied, overrides population assumptions and is
    marked patient-derived.
"""

from __future__ import annotations

import math
from typing import Any

from .evidence import Evidence, database, derived
from .model import Measure, NormalTissueProfile, TissueExpression, TumourState
from .providers import CRITICAL_TISSUES, Answer

#: nTPM above which a healthy tissue counts as expressing the gene at all.
DETECTION_NTPM = 1.0
#: nTPM above which on-target binding in a critical tissue is a real concern.
CONCERN_NTPM = 10.0

SPECIFICITY_HINT = {
    "Tissue enriched": 0.9,
    "Group enriched": 0.75,
    "Tissue enhanced": 0.6,
    "Low tissue specificity": 0.2,
    "Not detected": 1.0,
}


def _tissue_rows(row: dict[str, Any]) -> list[TissueExpression]:
    out: list[TissueExpression] = []
    for tissue, weight in CRITICAL_TISSUES.values():
        key = f"Tissue RNA - {tissue} [nTPM]"
        raw = row.get(key)
        if raw is None:
            # HPA drops null columns; keep the tissue with an explicit absence
            out.append(TissueExpression(tissue, None, "nTPM", None, weight >= 0.8, weight))
            continue
        try:
            value = float(raw)
        except (TypeError, ValueError):
            value = None
        out.append(TissueExpression(tissue, value, "nTPM", None, weight >= 0.8, weight))
    return out


def normal_profile(gene: str, answer: Answer) -> NormalTissueProfile:
    """Healthy-tissue expression and the on-target/off-tumour risk it implies."""
    profile = NormalTissueProfile()
    if not answer.available:
        profile.summary = "unavailable"
        profile.risk_basis = answer.reason
        profile.on_target_off_tumour_risk = "unknown"
        return profile
    row = answer.data
    profile.tissues = _tissue_rows(row)
    profile.specificity = row.get("RNA tissue specificity")
    profile.distribution = row.get("RNA tissue distribution")
    profile.evidence = list(answer.evidence)

    measured = [t for t in profile.tissues if t.value is not None]
    if not measured:
        profile.summary = "no per-tissue values returned"
        profile.on_target_off_tumour_risk = "unknown"
        profile.risk_basis = "Human Protein Atlas returned no consensus nTPM for the queried tissues"
        return profile

    at_risk = sorted(
        (t for t in measured if t.value is not None and t.value >= CONCERN_NTPM),
        key=lambda t: -(t.value or 0) * t.weight,
    )
    profile.tissues_at_risk = [f"{t.tissue} ({t.value:g} nTPM)" for t in at_risk[:8]]
    top = max(measured, key=lambda t: (t.value or 0) * t.weight)
    detected = [t for t in measured if (t.value or 0) >= DETECTION_NTPM]
    profile.summary = (
        f"detected in {len(detected)}/{len(measured)} queried healthy tissues; highest weighted "
        f"{top.tissue} at {top.value:g} nTPM"
    )
    profile.on_target_off_tumour_risk, profile.risk_basis = _risk(at_risk, top)
    profile.evidence.append(
        derived(
            "GenomeOS normal-tissue safety",
            f"{gene}: {profile.summary}; on-target/off-tumour risk {profile.on_target_off_tumour_risk} "
            f"({profile.risk_basis})",
            0.6,
        )
    )
    return profile


def _risk(at_risk: list[TissueExpression], top: TissueExpression) -> tuple[str, str]:
    essential = [t for t in at_risk if t.critical]
    if essential:
        worst = essential[0]
        level = "high" if (worst.value or 0) >= 5 * CONCERN_NTPM else "moderate"
        return level, (
            f"{worst.tissue} is an organ without a spare and carries {worst.value:g} nTPM "
            f"(threshold {CONCERN_NTPM:g})"
        )
    if at_risk:
        return "moderate", (
            f"{len(at_risk)} healthy tissues above {CONCERN_NTPM:g} nTPM, none of them critical organs; "
            f"highest {at_risk[0].tissue} at {at_risk[0].value:g}"
        )
    if (top.value or 0) >= DETECTION_NTPM:
        return "low", (
            f"no queried healthy tissue reaches {CONCERN_NTPM:g} nTPM; highest is {top.tissue} at "
            f"{top.value:g}"
        )
    return "low", "not detected above background in any queried healthy tissue"


def normal_tissue_safety(profile: NormalTissueProfile) -> tuple[float | None, str]:
    """0-1 safety score with the arithmetic stated, or None when unknown.

    The penalty is the largest weighted healthy-tissue level on a log scale,
    so a hundredfold difference in expression matters and a tenfold wobble at
    the bottom does not.
    """
    measured = [t for t in profile.tissues if t.value is not None]
    if not measured:
        return None, "no healthy-tissue measurements available"
    worst = max(measured, key=lambda t: (t.value or 0) * t.weight)
    load = (worst.value or 0.0) * worst.weight
    penalty = min(1.0, math.log10(1 + load) / math.log10(1 + 100.0))
    score = round(1.0 - penalty, 3)
    return score, (
        f"1 - log10(1 + {load:.1f}) / log10(101), where {load:.1f} is {worst.tissue} at "
        f"{worst.value:g} nTPM weighted {worst.weight:g} (Human Protein Atlas consensus, population-level)"
    )


def tumour_state(
    gene: str,
    alterations: list[str],
    patient_rna: Answer,
    normal: NormalTissueProfile,
    hpa_row: dict[str, Any] | None = None,
) -> TumourState:
    """What the tumour is known to do with this gene, and what is not known."""
    state = TumourState(altered=bool(alterations), alterations=list(alterations))
    if patient_rna.available:
        d = patient_rna.data
        state.expression = Measure(
            value=float(d["value"]),
            unit=d["unit"],
            source=f"patient tumour RNA-seq ({d.get('sample') or 'supplied'})",
            level="human",
            patient_specific=True,
            confidence=0.9,
        )
        state.evidence.extend(patient_rna.evidence)
    else:
        state.expression = Measure.unavailable(patient_rna.reason or "no tumour RNA-seq supplied")
        state.evidence.append(
            derived(
                "GenomeOS evidence layering",
                f"{gene}: DNA evidence present, RNA evidence unavailable, protein evidence unavailable; "
                "surface expression in this tumour cannot be established",
                0.9,
            )
        )
    if hpa_row:
        cancer = hpa_row.get("RNA cancer specificity")
        if cancer:
            state.evidence.append(
                database(
                    "Human Protein Atlas (cancer RNA)",
                    f"{gene} across tumour cohorts: {cancer} (population-level, not this patient)",
                    0.6,
                    "human",
                )
            )
    if normal.specificity in SPECIFICITY_HINT:
        state.evidence.append(
            database(
                "Human Protein Atlas (tissue specificity class)",
                f"{gene} is classed '{normal.specificity}' across healthy tissues "
                f"({normal.distribution or 'distribution unstated'})",
                0.7,
                "human",
            )
        )
    return state


def tumour_selectivity(
    gene: str, tumour: TumourState, normal: NormalTissueProfile
) -> tuple[float | None, str, list[Evidence]]:
    """How much more of this target the tumour carries than healthy tissue.

    With patient tumour RNA the ratio is computed and marked patient-derived.
    Without it, only the population specificity class is available, which is a
    hint about the gene, not a measurement of this tumour; it is capped low and
    labelled inferred.
    """
    measured = [t for t in normal.tissues if t.value is not None]
    if tumour.expression.known and tumour.expression.value is not None and measured:
        highest = max(measured, key=lambda t: t.value or 0)
        ratio = (tumour.expression.value + 1.0) / ((highest.value or 0.0) + 1.0)
        score = min(1.0, math.log10(max(ratio, 1e-3) + 1) / math.log10(21))  # 20x ratio saturates
        basis = (
            f"log-ratio of tumour {tumour.expression.value:g} {tumour.expression.unit} to the highest "
            f"healthy tissue ({highest.tissue} {highest.value:g} nTPM); 20x saturates the score"
        )
        ev = [
            derived(
                "GenomeOS tumour selectivity",
                f"{gene} is {ratio:.1f}x the highest queried healthy tissue in this tumour's RNA",
                0.75,
            )
        ]
        return round(score, 3), basis, ev
    if normal.specificity in SPECIFICITY_HINT:
        capped = min(0.5, SPECIFICITY_HINT[normal.specificity])
        basis = (
            f"no tumour RNA-seq; population tissue-specificity class '{normal.specificity}' used as a "
            "weak prior and capped at 0.5"
        )
        ev = [
            derived(
                "GenomeOS tumour selectivity",
                f"{gene} selectivity is a population-level prior, not a measurement of this tumour",
                0.3,
            )
        ]
        return capped, basis, ev
    return None, "no tumour expression and no healthy-tissue specificity class available", []


def target_density(gene: str, tumour: TumourState, normal: NormalTissueProfile) -> dict[str, Any]:
    """Relative abundance, kept in layers so RNA is never read as receptor count."""
    measured = [t for t in normal.tissues if t.value is not None]
    return {
        "absolute_surface_abundance": tumour.surface_abundance.to_dict(),
        "tumour_protein_abundance": tumour.protein_abundance.to_dict(),
        "tumour_rna": tumour.expression.to_dict(),
        "normal_tissue_rna_max": (
            {
                "tissue": max(measured, key=lambda t: t.value or 0).tissue,
                "value": max(t.value or 0 for t in measured),
                "unit": "nTPM",
                "level": "human",
                "patient_specific": False,
            }
            if measured
            else None
        ),
        "copy_number": tumour.copy_number.to_dict(),
        "note": (
            "RNA level is not a receptor count. Converting nTPM into receptors per cell requires "
            "quantitative surface proteomics, which is not supplied; the layers are kept separate."
        ),
    }
