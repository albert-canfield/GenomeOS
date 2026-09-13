"""Normal-cell protection as a first-class result, not a footnote.

A candidate must not rank highly merely because it is abundant on cancer.
This module turns the healthy-tissue profile into the thing a therapy
designer actually needs: a named list of tissues to stay away from, the
weighted worst case, and the explicit on-target/off-tumour risk that a binder
would carry into the patient.

The weighting is a stated policy of this pipeline, not a measurement: tissue
whose loss cannot be replaced or transplanted is weighted hardest.
"""

from __future__ import annotations

from typing import Any

from .evidence import Evidence, derived
from .expression import CONCERN_NTPM, DETECTION_NTPM
from .model import NormalTissueProfile

ESSENTIAL = (
    "heart muscle",
    "cerebral cortex",
    "cerebellum",
    "hypothalamus",
    "bone marrow",
    "lung",
    "liver",
    "kidney",
)


def exclusions(gene: str, profile: NormalTissueProfile) -> dict[str, Any]:
    """What a binder must avoid, and how firmly that is established."""
    measured = [t for t in profile.tissues if t.value is not None]
    unmeasured = [t.tissue for t in profile.tissues if t.value is None]
    avoid = sorted(
        (t for t in measured if (t.value or 0) >= CONCERN_NTPM),
        key=lambda t: -(t.value or 0) * t.weight,
    )
    watch = sorted(
        (t for t in measured if DETECTION_NTPM <= (t.value or 0) < CONCERN_NTPM),
        key=lambda t: -(t.value or 0) * t.weight,
    )
    return {
        "gene": gene,
        "on_target_off_tumour_risk": profile.on_target_off_tumour_risk,
        "basis": profile.risk_basis,
        "avoid": [
            {
                "tissue": t.tissue,
                "value": t.value,
                "unit": t.unit,
                "essential_organ": t.tissue in ESSENTIAL,
                "constraint": f"binding must not exceed a level that damages {t.tissue}, which carries "
                f"{t.value:g} {t.unit} of this transcript in healthy donors",
            }
            for t in avoid[:10]
        ],
        "watch": [{"tissue": t.tissue, "value": t.value, "unit": t.unit} for t in watch[:10]],
        "not_measured": unmeasured,
        "thresholds": {
            "detection_nTPM": DETECTION_NTPM,
            "concern_nTPM": CONCERN_NTPM,
            "note": "thresholds are a stated policy of this pipeline applied to Human Protein Atlas "
            "consensus RNA, not a validated clinical cutoff",
        },
        "layers_not_established": [
            "protein-level (not RNA) abundance in healthy tissue",
            "surface density per healthy cell",
            "isoform differences between tumour and healthy tissue",
            "tumour-specific post-translational state, for example altered glycosylation",
        ],
    }


def risk_evidence(gene: str, profile: NormalTissueProfile) -> Evidence:
    return derived(
        "GenomeOS normal-tissue safety",
        f"{gene} on-target/off-tumour risk is {profile.on_target_off_tumour_risk}: {profile.risk_basis}",
        0.6 if profile.known else 0.0,
    )
