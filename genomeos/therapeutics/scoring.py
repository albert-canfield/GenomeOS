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

#: Dimensions that describe the candidate's inputs rather than the target's
#: priority; scored and published, but kept out of the overall mean.
#:
#: `alteration_evidence` is here for a reason worth writing down, because it was
#: first written as a weighted dimension and measured before it was believed. The
#: mean answers how good a target this protein is on the dimensions that could be
#: measured; the tier answers whether there is evidence that this tumour involves
#: the gene at all, which is a precondition and not a dimension of goodness.
#: Averaged in, a precondition compensates for weak biology and weak biology
#: compensates for a missing precondition — and it did worse than that in
#: practice: because the mean is taken over the dimensions that were available, a
#: new dimension lifts a candidate with fewer of them further, and because the
#: poor-safety cap clips at a constant, a capped candidate can be overtaken by
#: adding any dimension at all. Both inverted a tumour with twelve copies of
#: ERBB2, putting a mutated PIK3CA first. So the tier gates instead: it orders
#: the candidates in `pipeline.analyse` and it is published per candidate with
#: the sentence behind it, in the same family as safety, which can only ever cap
#: a score and never raise one.
DIAGNOSTIC = ("shedding", "alteration_evidence")

#: Copy number bounds what a cell could display; it is not a measurement of
#: what it does, so it cannot reach the score a measurement can.
COPY_NUMBER_CAP = 0.5

#: A discrete copy-number call carries a direction and no amount, so it scores
#: below what the weakest count producing that call would score. An
#: amplification is at least 4 copies, which the formula above would put at
#: 0.25; a gain is at least 3, which it would put at 0.125.
DISCRETE_COPY_NUMBER: dict[str, float] = {"amplification": 0.2, "gain": 0.1}

#: The count at which a copy number is an alteration rather than a measurement
#: that the gene is normal; the same 4 copies the note above reasons from.
RAISED_COPIES = 4.0

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


#: What is known about *this gene in this patient*, which is not the same
#: question as what is known about the protein. Curated localisation describes
#: the gene whether or not the tumour touched it, so without this dimension a
#: gene carried at twelve copies and a gene named only for neighbouring a
#: mutated one are scored from the same annotation and come out equal.
#:
#: The order is argued rather than assumed. The top two tiers are both
#: measurements of this gene in this patient and differ in degree: a somatic
#: alteration is attributable to the tumour and is what approved indications are
#: written on, while the expression route establishes a ratio against a queried
#: healthy-tissue panel rather than this patient's own normal tissue, and raised
#: transcript is not protein on the surface. The bottom tier differs in kind: it
#: is not a weak measurement of this gene but a measurement of a different gene
#: plus a database association. It is 0.2 and not 0.0 because an association
#: with a disrupted driver is evidence of something and the pipeline is entitled
#: to propose it; it is not entitled to prefer it.
EVIDENCE_TIERS: dict[str, float] = {
    "observed_alteration": 1.0,
    "patient_measurement": 0.6,
    "association_hypothesis": 0.2,
}

#: The tier that is not preferred over the target: no measurement of this gene
#: in this patient at all. `pipeline.analyse` orders this tier below the other
#: two before it looks at the score.
UNMEASURED_TIER = "association_hypothesis"


# --- two rules on the order, registered 2026-09-27 ------------------------------------
#
# The evidence tier asks whether this tumour involves the gene at all. Two questions it
# does not ask were named when it landed. First, nothing asks whether any modelled
# modality reaches the candidate: in the ERBB2-amplified tumour the second place goes to
# a gene whose best mechanism is absent rather than weak. Second, inside the top tier
# twelve copies of a gene and one missense read the same, because the tier is a label and
# carries no amount.
#
# Both are answered on the ORDER and neither enters WEIGHTS, for the reason the tier's
# own first shape established by being measured: the mean is taken over the dimensions
# that were available, so any new dimension lifts the candidate with fewer of them
# further, and the poor-safety cap clips at a constant that an uncapped candidate walks
# past. A precondition averaged against biology also lets each compensate for the other.
# So one is a gate and one is a tiebreak, in the same family as safety, which can only
# ever cap a score and never raise one.

#: The mechanism gate. Whether some modelled modality reaches this candidate with its
#: hard requirements *answered* — `TherapeuticTargetCandidate.best_mechanism` — or not.
#:
#: Two classes and not three, which is the argued part. A mechanism with an unanswered
#: requirement has not been shown to apply; it has only failed to be ruled out, and that
#: is this file's existing rule for which mechanism may head a list. So "provisional
#: only" cannot be preferred over "nothing at all": in the EML4-ALK tumour that would put
#: a mutated PIK3CA, whose nearest mechanism is adcp at 0.25 with `surface_accessible`
#: unanswered, above the fusion the tumour actually carries.
REACH_ESTABLISHED = "established_mechanism"
REACH_NONE = "no_established_mechanism"

#: The gate is read after the evidence tier and before the score, never before the tier:
#: a candidate measured nowhere in this patient does not rise by having a reachable
#: surface, which is the defect the tier exists to prevent.
MECHANISM_GATE_ORDER = (REACH_ESTABLISHED, REACH_NONE)

#: The magnitude tiebreak. The quantities that say *how much*, each read from this
#: patient's own tumour data and never imputed: the copy count in the patient's
#: copy-number table, the variant allele fraction of the observed variant, and whether
#: the observed position is a recorded hotspot. A quantity that is absent is reported as
#: absent; no default fraction, no count inferred from a discrete call, no hotspot
#: inferred from a gene's driver frequency.
#:
#: They are compared like with like and in this order — both candidates carry a count, so
#: the larger count; else both carry a fraction, so the larger; else exactly one is a
#: hotspot. Copies, fractions and a yes/no share no unit, and inventing an exchange rate
#: between them would be the guess this rule exists to refuse, so a pair that has no
#: quantity in common leaves the tie unbroken and the gene-name fallback stands.
MAGNITUDE_QUANTITIES = ("copies", "vaf", "hotspot")

#: What a copy count is read against, so that a magnitude is a distance from normal
#: rather than a raw number.
DIPLOID_COPIES = 2.0


def alteration_evidence(origins: Any, tumour: TumourState) -> tuple[str, float, str]:
    """Which tier of evidence about this gene in this patient reached it here.

    Read from the evidence and never from the name of the route that proposed
    the candidate. That is deliberate: CD19 is reached from the patient's RNA
    and carries no alteration of any kind, so "has a DNA origin" would demote
    the target of four approved therapies. And a gene first proposed as a
    neighbour rises on its own the moment the patient's RNA is supplied, without
    the tier being edited.
    """
    alterations = [o for o in (origins or [])]
    if alterations:
        labels = [
            o.alteration_label or f"{o.variant_type} {o.protein_change or ''}".strip() or "alteration"
            for o in alterations
        ]
        return (
            "observed_alteration",
            EVIDENCE_TIERS["observed_alteration"],
            (
                "this tumour's own DNA carries "
                + ", ".join(dict.fromkeys(labels))
                + "; the alteration is observed in this patient and not inferred"
            ),
        )
    measured = []
    if tumour.expression.value is not None and tumour.expression.patient_specific:
        measured.append(
            f"tumour RNA at {tumour.expression.value:g} {tumour.expression.unit} ({tumour.expression.source})"
        )
    if tumour.protein_abundance.known and tumour.protein_abundance.patient_specific:
        measured.append(f"tumour proteomics ({tumour.protein_abundance.source})")
    if tumour.surface_abundance.known and tumour.surface_abundance.patient_specific:
        measured.append(f"tumour surface proteomics ({tumour.surface_abundance.source})")
    cn = tumour.copy_number
    if cn.patient_specific and (
        (cn.value is not None and cn.value >= RAISED_COPIES) or cn.qualitative in DISCRETE_COPY_NUMBER
    ):
        # A raised count in the patient's copy-number table is an observed
        # alteration whether or not a caller also emitted an event record for
        # it, because the tier reads the evidence and not the record-keeping.
        basis = (
            f"this patient's copy-number table places this gene at {cn.value:g} copies against the diploid 2"
            if cn.value is not None
            else f"this patient's copy-number call for this gene is '{cn.qualitative}' ({cn.source})"
        )
        return "observed_alteration", EVIDENCE_TIERS["observed_alteration"], basis
    if measured:
        return (
            "patient_measurement",
            EVIDENCE_TIERS["patient_measurement"],
            (
                "not altered in this tumour's DNA, and measured in this patient: "
                + "; ".join(measured)
                + ". A measurement of this gene in this patient, scored below an observed alteration "
                "because the comparator is a healthy-tissue panel rather than this patient's own normal "
                "tissue and because raised transcript is not protein on the surface"
            ),
        )
    return (
        "association_hypothesis",
        EVIDENCE_TIERS["association_hypothesis"],
        (
            "nothing about this gene was measured in this patient: it is on the list because a database "
            "associates it with a gene that was altered, which is evidence about the other gene"
        ),
    )


def mechanism_reach(candidate: Any) -> tuple[str, str]:
    """Does any modelled modality reach this candidate, requirements answered?

    Registered 2026-09-27 as a gate on the order and not a dimension of the mean.
    `best_mechanism` is the question: a mechanism whose hard requirement is merely
    unanswered has not been shown to apply, so a candidate carrying only such a
    mechanism is in the same class as one carrying none. Ranking an unanswered
    question above a gene the tumour altered is the defect this rule must not
    introduce while closing the one it was written for.
    """
    best = candidate.best_mechanism
    if best is not None:
        return REACH_ESTABLISHED, (
            f"{best.mechanism} reaches this target with every hard requirement answered, at "
            f"compatibility {best.compatibility:.2f}"
        )
    nearest = candidate.best_provisional_mechanism
    if nearest is not None:
        open_ = ", ".join(nearest.provisional_requirements) or "an unnamed requirement"
        return REACH_NONE, (
            f"no modelled modality is established against this target; the nearest is "
            f"{nearest.mechanism} at compatibility {nearest.compatibility:.2f} with {open_} "
            "unanswered, which is a question and not an option"
        )
    refused = sorted({g for m in candidate.therapeutic_mechanisms for g in m.gates_failed})
    if candidate.therapeutic_mechanisms:
        return REACH_NONE, (
            "every modelled modality was refused against this target"
            + (f" ({', '.join(refused[:4])})" if refused else "")
        )
    return REACH_NONE, "no modelled modality applies to this target at all"


def alteration_magnitude(origins: Any, tumour: TumourState) -> tuple[dict[str, Any], str]:
    """How much this tumour altered the gene, from measurements only.

    Registered 2026-09-27. The evidence tier says an alteration was observed; it
    carries no amount, so twelve copies and a single missense read the same inside
    the top tier. The quantities here are the ones this patient's own data state,
    and an absent quantity is reported as absent: no default allele fraction, no
    count inferred from a discrete `amplification` call, no hotspot inferred from a
    gene's driver frequency. `magnitude_prefers` compares them like with like, and
    only between candidates a score could not separate.
    """
    rows = list(origins or [])
    counts = [o.copy_number for o in rows if o.copy_number is not None]
    cn = tumour.copy_number
    if not counts and cn.patient_specific and cn.value is not None:
        counts = [cn.value]
    fractions = [o.vaf for o in rows if o.vaf is not None]
    if not fractions and tumour.vaf is not None:
        fractions = [tumour.vaf]
    positioned = [o for o in rows if o.position is not None or o.protein_change]
    copies = max(counts) if counts else None
    vaf = max(fractions) if fractions else None
    hotspot = any(o.hotspot for o in positioned) if positioned else None
    magnitude: dict[str, Any] = {
        "copies": copies,
        "copies_above_diploid": None if copies is None else round(copies - DIPLOID_COPIES, 3),
        "vaf": vaf,
        "hotspot": hotspot,
        "quantities": [
            q for q in MAGNITUDE_QUANTITIES if magnitude_value(q, copies, vaf, hotspot) is not None
        ],
    }
    said = []
    if copies is not None:
        said.append(f"{copies:g} copies against the diploid {DIPLOID_COPIES:g} in this patient's table")
    else:
        said.append("no copy count is present in this patient's data")
    if vaf is not None:
        said.append(f"variant allele fraction {vaf:.2f} at the observed variant")
    else:
        said.append("no variant allele fraction was reported")
    if hotspot is None:
        said.append("no observed position to look up, so hotspot status does not apply")
    else:
        said.append(
            "the observed position is a recorded hotspot"
            if hotspot
            else "the observed position is not a recorded hotspot"
        )
    return magnitude, "; ".join(said) + ". Absent quantities are absent and are not imputed"


def magnitude_value(quantity: str, copies: float | None, vaf: float | None, hotspot: bool | None) -> Any:
    return {"copies": copies, "vaf": vaf, "hotspot": hotspot}[quantity]


def magnitude_prefers(a: dict[str, Any], b: dict[str, Any]) -> int:
    """-1 if `a` carries the larger measured magnitude, 1 if `b`, 0 if neither.

    Like with like, in the registered order, and never across kinds: copies, a
    fraction and a yes/no share no unit, so a pair with no quantity in common
    leaves the tie unbroken rather than being ordered by an invented exchange rate.
    """
    for quantity in MAGNITUDE_QUANTITIES:
        av, bv = (a or {}).get(quantity), (b or {}).get(quantity)
        if av is None or bv is None or av == bv:
            continue
        if quantity == "hotspot":
            return -1 if av else 1
        return -1 if av > bv else 1
    return 0


def surface_accessibility(localisation: Any, ectodomain_lost: bool = False) -> tuple[float | None, str]:
    """How reachable the protein is from outside, from curated localisation.

    Curated localisation describes the full-length protein. When the tumour
    makes a fusion that drops this gene's N-terminal ectodomain, the curated
    figure is a description of a protein the tumour does not have, so it is
    not the number to score.
    """
    if ectodomain_lost:
        return 0.0, (
            "this gene is the 3' partner of the fusion, so the product begins with the partner's "
            "sequence and carries neither its signal peptide nor its N-terminal extracellular "
            "domain; the curated topology describes the full-length protein, which this tumour "
            "does not make"
        )
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
    if m.value is not None:
        score = min(1.0, math.log10(1 + m.value) / math.log10(1 + 100.0))
        return round(score, 3), (
            f"log10(1 + {m.value:g} {m.unit}) / log10(101) from this patient's tumour RNA-seq; "
            "saturates at 100"
        )
    cn = tumour.copy_number
    if cn.value is not None:
        # Copy number is not expression. It bounds how much a cell could make,
        # which is weaker than a measurement and is capped to say so.
        score = min(COPY_NUMBER_CAP, max(0.0, (cn.value - 2.0) / 8.0))
        return round(score, 3), (
            f"no tumour RNA-seq; {cn.value:g} copies against the diploid 2, scored as "
            f"({cn.value:g} - 2) / 8 and capped at {COPY_NUMBER_CAP:g} because copy number bounds what "
            "a cell could make and never shows that it does"
        )
    if cn.qualitative in DISCRETE_COPY_NUMBER:
        # A discrete call says which direction the copy number moved and by how
        # much it does not say. It therefore scores below the weakest count that
        # would produce the same call, rather than being given one.
        score = DISCRETE_COPY_NUMBER[cn.qualitative]
        return score, (
            f"no tumour RNA-seq and no copy count; the call is '{cn.qualitative}' ({cn.source}), scored "
            f"at {score:g}, below what an actual count of copies could reach, because a discrete call "
            "gives a direction and no amount"
        )
    return None, m.reason or "no tumour RNA-seq or proteomics supplied"


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
