"""Component scores, and one overall score whose arithmetic is on the page.

A single opaque number is useless for research: it cannot be argued with. So
every dimension is scored on its own, each carries the sentence that explains
how it was obtained, and the overall score is a stated weighted mean of the
dimensions that were actually available.

Two rules keep it honest:

  * an unknown dimension is dropped from the mean, never defaulted to a
    middling value, and the number of dropped dimensions is reported as
    component coverage;
  * safety can only ever lower the result. A candidate cannot rank highly
    because the normal-tissue data is missing, and it cannot rank highly with
    a known bad normal-tissue profile either.
"""

from __future__ import annotations

import math
from typing import Any

from .evidence import Evidence, derived
from .model import ScoreComponent, ScoreSet, TumourState

#: Weight of each dimension in the overall research-priority score.
WEIGHTS: dict[str, float] = {
    "surface_accessibility": 1.0,
    "tumour_selectivity": 1.3,
    "tumour_expression": 0.9,
    "normal_tissue_safety": 1.5,
    "clonality": 0.8,
    "internalisation": 0.5,
    "lysosomal_trafficking": 0.3,
    "structural_bindability": 0.6,
    "neoantigen_strength": 0.9,
    "presentation_confidence": 0.7,
    "evidence_strength": 0.8,
    "clinical_precedent": 0.7,
}

#: Dimensions that describe the mechanism's inputs rather than the target's
#: priority; scored and published, but kept out of the overall mean.
DIAGNOSTIC = ("shedding",)

SAFETY_UNKNOWN_CAP = 0.6
SAFETY_POOR_CAP = 0.5
SAFETY_POOR_THRESHOLD = 0.35


def component(
    key: str,
    value: float | None,
    basis: str,
    evidence: list[Evidence] | None = None,
    unknown_reason: str = "",
) -> ScoreComponent:
    return ScoreComponent(
        key=key,
        value=None if value is None else round(max(0.0, min(1.0, value)), 3),
        basis=basis if value is not None else "",
        weight=WEIGHTS.get(key, 0.0),
        unknown_reason=unknown_reason or (basis if value is None else ""),
        evidence=evidence or [],
    )


def surface_accessibility(localisation: Any) -> tuple[float | None, str]:
    """How reachable the protein is from outside, from curated localisation."""
    if localisation.plasma_membrane is False:
        return 0.0, (
            f"curated localisation places {localisation.primary} inside the cell; a circulating binder "
            "cannot reach it"
        )
    best = max(
        (
            localisation.compartments.get(c, 0.0)
            for c in ("plasma_membrane", "cell_surface", "transmembrane", "gpi_anchored")
        ),
        default=0.0,
    )
    if best == 0.0:
        return None, "no localisation evidence places this protein at or outside the membrane"
    if localisation.extracellular_regions:
        longest = max(localisation.extracellular_regions, key=lambda r: (r.end or 0) - (r.start or 0))
        span = (longest.end or 0) - (longest.start or 0) + 1
        bonus = min(0.1, span / 5000)
        return min(1.0, best + bonus), (
            f"curated membrane localisation at confidence {best:.2f}, plus {span} residues of curated "
            f"extracellular topology ({longest.start}-{longest.end}) worth +{bonus:.2f}"
        )
    return best * 0.8, (
        f"membrane localisation at confidence {best:.2f}, reduced by 20% because no extracellular "
        "topological domain is curated"
    )


def clonality(tumour: TumourState) -> tuple[float | None, str]:
    """Present on every tumour cell, or only a subclone?"""
    if tumour.clonality == "clonal":
        return 1.0, f"variant allele fraction {tumour.vaf:.2f} indicates a clonal alteration"
    if tumour.clonality == "subclonal" and tumour.vaf is not None:
        return round(min(1.0, tumour.vaf * 2), 3), (
            f"variant allele fraction {tumour.vaf:.2f}; an antigen on part of the tumour leaves the rest "
            "untouched, so the fraction scales the score"
        )
    return None, tumour.clonality_reason


def tumour_expression_score(tumour: TumourState) -> tuple[float | None, str]:
    """Only patient measurements count. DNA alone never establishes expression."""
    m = tumour.expression
    if m.value is None:
        return None, m.reason or "no tumour RNA-seq or proteomics supplied"
    score = min(1.0, math.log10(1 + m.value) / math.log10(1 + 100.0))
    return round(score, 3), (
        f"log10(1 + {m.value:g} {m.unit}) / log10(101) from this patient's tumour RNA-seq; saturates at 100"
    )


def clinical_precedent(precedent: dict[str, Any] | None, available: bool) -> tuple[float | None, str]:
    """Has anything ever been aimed at this target, and how far did it get?"""
    if not available:
        return None, "therapeutic precedent could not be retrieved"
    if not precedent:
        return 0.0, "no drug or clinical candidate is recorded against this target"
    rows = ((precedent.get("drugAndClinicalCandidates") or {}).get("rows")) or []
    approved = [r for r in rows if r.get("maxClinicalStage") == "APPROVAL"]
    trials = [r for r in rows if str(r.get("maxClinicalStage", "")).startswith("PHASE")]
    tract = [t for t in precedent.get("tractability") or [] if t.get("modality") == "AB" and t.get("value")]
    if approved:
        return 1.0, f"{len(approved)} approved therapies already engage this target"
    if trials:
        return 0.7, f"{len(trials)} agents against this target have entered clinical trials"
    if tract:
        return 0.4, (
            "no clinical agent, but Open Targets places this target in antibody-tractability buckets: "
            + ", ".join(t["label"] for t in tract[:3])
        )
    return 0.0, "no drug, trial or antibody-tractability evidence recorded against this target"


def shedding_score(trafficking: Any) -> tuple[float | None, str]:
    """1 means antigen is shed, which every surface mechanism dislikes."""
    if trafficking.shedding is None:
        return None, "no evidence either way on ectodomain shedding"
    if trafficking.shedding:
        return 1.0, "curated annotation indicates a released or secreted form circulates"
    return 0.0, "no evidence of ectodomain shedding"


def assemble(components: list[ScoreComponent], evidence_confidence: float) -> ScoreSet:
    """Weighted mean over available components, with every adjustment named."""
    s = ScoreSet(components={c.key: c for c in components})
    scored = [c for c in components if c.value is not None and c.key not in DIAGNOSTIC and c.weight > 0]
    declared = sum(WEIGHTS[k] for k in WEIGHTS)
    used = sum(c.weight for c in scored)
    s.coverage = round(used / declared, 3) if declared else 0.0
    if not scored:
        s.overall = None
        s.formula = "no scoreable dimension was available"
        s.confidence = 0.0
        return s
    raw = sum(c.weight * (c.value or 0.0) for c in scored) / used
    s.formula = (
        "overall = sum(weight_i * score_i) / sum(weight_i) over the "
        f"{len(scored)} available dimensions ({', '.join(sorted(c.key for c in scored))}), then the "
        "safety adjustments below; dimensions with no data are excluded rather than defaulted"
    )
    value = raw
    safety = s.components.get("normal_tissue_safety")
    if safety is None or safety.value is None:
        value = min(value, SAFETY_UNKNOWN_CAP)
        s.adjustments.append(
            {
                "adjustment": "safety cap (unknown)",
                "from": round(raw, 3),
                "to": round(value, 3),
                "reason": "normal-tissue expression is not established; a target cannot rank highly on "
                "missing safety data",
            }
        )
    elif safety.value < SAFETY_POOR_THRESHOLD:
        capped = min(value, SAFETY_POOR_CAP)
        if capped < value:
            s.adjustments.append(
                {
                    "adjustment": "safety cap (poor)",
                    "from": round(value, 3),
                    "to": round(capped, 3),
                    "reason": f"normal-tissue safety scores {safety.value:.2f}, below the "
                    f"{SAFETY_POOR_THRESHOLD:g} threshold",
                }
            )
        value = capped
    s.overall = round(value, 3)
    s.confidence = round(s.coverage * max(0.2, evidence_confidence), 3)
    return s


def summary_evidence(s: ScoreSet, gene: str) -> Evidence:
    missing = sorted(k for k, c in s.components.items() if c.value is None and c.weight > 0)
    return derived(
        "GenomeOS scoring",
        f"{gene} scores {s.overall if s.overall is not None else 'n/a'} over "
        f"{s.coverage:.0%} of the declared dimensions"
        + (f"; unavailable: {', '.join(missing)}" if missing else ""),
        s.confidence,
    )
