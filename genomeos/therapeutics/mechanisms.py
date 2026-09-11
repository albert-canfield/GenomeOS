"""Therapeutic mechanism selection.

Not "which payload", because most therapeutic mechanisms carry no payload at
all. An antibody that recruits NK cells has cargo `none` and effector `NK
cells`, and it is no less a mechanism for that. The five concepts stay
separate everywhere: target, binder, mechanism, cargo, effector.

Each mechanism declares what it needs, in two strengths:

    gates    biology without which the mechanism cannot work at all.
             A failed gate sets compatibility to zero and says why.
    weights  properties that make it work better or worse, combined into a
             compatibility score with every contribution listed.

Unknown inputs are never filled in. A weighted term whose input is unknown is
dropped from the average and recorded as a blocking unknown, so a mechanism
cannot score well on ignorance.

Compatibility is a computational compatibility score. It is not a probability
that a patient responds, and it is labelled as such wherever it is emitted.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

from .evidence import Evidence, derived
from .model import Binder, Cargo, Effector, Grade, MechanismFit, MechanismStatus

#: How mature the mechanism class is, independent of any particular target.
#: Applied as a visible multiplier, never hidden inside the weights.
STATUS_FACTOR: dict[str, float] = {
    "established": 1.0,
    "clinical": 0.95,
    "preclinical": 0.75,
    "experimental": 0.6,
}

#: A mechanism evaluated on few of its declared inputs is damped to this
#: fraction at zero coverage, so ignorance never reads as compatibility.
COVERAGE_FLOOR = 0.4

REPAIR_CAVEAT = (
    "Correcting or replacing one gene does not restore a normal cell. An established malignancy "
    "carries several cooperating drivers, passenger mutations, copy-number change, aneuploidy and "
    "epigenetic dysregulation; reverting one of them removes one dependency, not the disease."
)


@dataclass(frozen=True, slots=True)
class Gate:
    """A hard requirement, with the reason it exists."""

    key: str
    test: Callable[[Any], bool | None]
    unmet: str
    unknown_is_failure: bool = False


@dataclass(frozen=True, slots=True)
class MechanismSpec:
    """What a mechanism is and what biology it needs."""

    mechanism: str
    summary: str
    binder_kind: str
    binder_recognises: str
    cargo_kind: str
    payload_required: bool
    effector_kind: str
    immune_components: tuple[str, ...]
    status: MechanismStatus
    gates: tuple[Gate, ...] = ()
    weights: tuple[tuple[str, float], ...] = ()  # score key, weight (higher is better)
    penalties: tuple[tuple[str, float], ...] = ()  # score key, weight (higher is worse)
    immunogenic_cell_death: Grade = "unknown"
    antigen_release: Grade = "unknown"
    notes: tuple[str, ...] = ()
    precedent_drug_types: tuple[str, ...] = ()


# --- gate predicates ------------------------------------------------------------------


def _surface(c: Any) -> bool | None:
    if c.localization.plasma_membrane is False:
        return False
    if c.localization.reachable:
        return True
    return None


def _internalises(c: Any) -> bool | None:
    return c.trafficking.internalises


def _peptide_hla(c: Any) -> bool | None:
    if c.neoantigen is None:
        return None
    if c.neoantigen.novel_peptide_sequence is False:
        return False
    return bool(c.neoantigen.peptides) or None


def _hla_typed(c: Any) -> bool | None:
    return bool(c.neoantigen and c.neoantigen.hla_alleles)


def _protein_altering(c: Any) -> bool | None:
    return any(o.protein_change for o in c.origins) or None


LOF_TYPES = {
    "stop_gained",
    "frameshift_variant",
    "start_lost",
    "splice_acceptor_variant",
    "splice_donor_variant",
}


def _loss_of_function(c: Any) -> bool | None:
    return any(o.variant_type in LOF_TYPES for o in c.origins)


SURFACE_GATE = Gate(
    "surface_accessible",
    _surface,
    "the protein is not reachable from outside the cell, so a circulating binder cannot engage it",
    unknown_is_failure=False,
)
INTERNALISATION_GATE = Gate(
    "internalisation",
    _internalises,
    "the receptor is annotated as not internalising, so nothing can be delivered inside through it",
)
DELIVERY_GATE = Gate(
    "intracellular_delivery",
    _internalises,
    "intracellular cargo delivery requires the target to carry the binder inside; no internalisation "
    "evidence exists",
    unknown_is_failure=True,
)
PEPTIDE_GATE = Gate(
    "peptide_hla_target",
    _peptide_hla,
    "no mutation-derived peptide could be derived, so there is no peptide/HLA complex to recognise",
    unknown_is_failure=True,
)
HLA_GATE = Gate(
    "patient_hla",
    _hla_typed,
    "the patient's HLA genotype is not available, and peptide presentation is allele-specific",
    unknown_is_failure=True,
)
MUTATION_GATE = Gate(
    "mutation_rationale",
    _protein_altering,
    "no protein-altering alteration in this gene, so there is no mutation-specific rationale",
    unknown_is_failure=True,
)
LOF_GATE = Gate(
    "loss_of_function",
    _loss_of_function,
    "restoration only applies to an inactivated gene; no loss-of-function alteration is present",
    unknown_is_failure=True,
)


# --- the ontology ---------------------------------------------------------------------

MECHANISMS: dict[str, MechanismSpec] = {}


def _add(spec: MechanismSpec) -> None:
    MECHANISMS[spec.mechanism] = spec


_add(
    MechanismSpec(
        "blocking_antibody",
        "a binder occupies the target and stops its signalling or ligand binding",
        "antibody-like",
        "an accessible extracellular epitope",
        "none",
        False,
        "the tumour cell's own blocked pathway",
        (),
        "established",
        (SURFACE_GATE,),
        (
            ("surface_accessibility", 1.0),
            ("tumour_selectivity", 0.8),
            ("structural_bindability", 0.7),
            ("normal_tissue_safety", 0.8),
        ),
        (("shedding", 0.4),),
        "hypothetical",
        "unknown",
        (
            "only useful where the target actively drives the tumour; blocking a passenger antigen "
            "achieves nothing",
        ),
        ("Antibody",),
    )
)
_add(
    MechanismSpec(
        "agonist_antibody",
        "a binder deliberately switches the target on, for example a death receptor",
        "antibody-like",
        "an accessible extracellular epitope",
        "none",
        False,
        "the target receptor's own signalling",
        (),
        "clinical",
        (SURFACE_GATE,),
        (
            ("surface_accessibility", 1.0),
            ("tumour_selectivity", 0.9),
            ("structural_bindability", 0.6),
            ("normal_tissue_safety", 0.9),
        ),
        (),
        "probable",
        "probable",
        ("requires the target to have an activating function worth triggering",),
        ("Antibody",),
    )
)
_add(
    MechanismSpec(
        "adcc",
        "an antibody coats the tumour cell and NK cells and macrophages destroy it",
        "antibody-like",
        "an accessible extracellular epitope",
        "none",
        False,
        "endogenous immune system",
        ("NK cells", "macrophages"),
        "established",
        (SURFACE_GATE,),
        (
            ("surface_accessibility", 1.0),
            ("tumour_expression", 0.9),
            ("tumour_selectivity", 0.9),
            ("normal_tissue_safety", 0.9),
            ("structural_bindability", 0.5),
        ),
        (("shedding", 0.5), ("internalisation", 0.3)),
        "probable",
        "probable",
        (
            "needs no payload and no intracellular delivery",
            "rapid internalisation works against it: the antibody must stay on the surface to be seen",
            "adequate surface density matters more than absolute affinity",
        ),
        ("Antibody",),
    )
)
_add(
    MechanismSpec(
        "adcp",
        "an antibody marks the cell and macrophages engulf it",
        "antibody-like",
        "an accessible extracellular epitope",
        "none",
        False,
        "endogenous immune system",
        ("macrophages",),
        "clinical",
        (SURFACE_GATE,),
        (
            ("surface_accessibility", 1.0),
            ("tumour_expression", 0.9),
            ("tumour_selectivity", 0.9),
            ("normal_tissue_safety", 0.9),
        ),
        (("shedding", 0.5),),
        "probable",
        "established",
        ("phagocytosis delivers tumour antigen to antigen-presenting cells",),
        ("Antibody",),
    )
)
_add(
    MechanismSpec(
        "complement_recruitment",
        "an antibody fixes complement on the tumour cell surface",
        "antibody-like",
        "a densely expressed extracellular epitope",
        "none",
        False,
        "complement cascade",
        ("complement C1q", "membrane attack complex"),
        "clinical",
        (SURFACE_GATE,),
        (
            ("surface_accessibility", 1.0),
            ("tumour_expression", 1.0),
            ("tumour_selectivity", 0.8),
            ("normal_tissue_safety", 0.9),
        ),
        (("shedding", 0.5),),
        "probable",
        "probable",
        ("depends strongly on antigen density and on the tumour's complement regulators",),
        ("Antibody",),
    )
)
_add(
    MechanismSpec(
        "t_cell_engager",
        "a bispecific binder holds a T cell against the tumour cell",
        "bispecific antibody-like",
        "an accessible tumour antigen on one arm, CD3 on the other",
        "none",
        False,
        "T cell",
        ("CD3+ T cells",),
        "established",
        (SURFACE_GATE,),
        (
            ("surface_accessibility", 1.0),
            ("tumour_selectivity", 1.2),
            ("normal_tissue_safety", 1.4),
            ("tumour_expression", 0.8),
        ),
        (("shedding", 0.4),),
        "probable",
        "established",
        (
            "needs no payload",
            "the most potent mechanism here and therefore the least forgiving: any expression in "
            "healthy tissue becomes on-target off-tumour toxicity, and low antigen density can still "
            "suffice for killing",
        ),
        ("Antibody",),
    )
)
_add(
    MechanismSpec(
        "nk_cell_engager",
        "a bispecific binder holds an NK cell against the tumour cell",
        "bispecific antibody-like",
        "an accessible tumour antigen on one arm, an NK activating receptor on the other",
        "none",
        False,
        "NK cell",
        ("NK cells",),
        "clinical",
        (SURFACE_GATE,),
        (
            ("surface_accessibility", 1.0),
            ("tumour_selectivity", 1.0),
            ("normal_tissue_safety", 1.2),
            ("tumour_expression", 0.8),
        ),
        (("shedding", 0.4),),
        "probable",
        "probable",
        ("narrower toxicity profile than a T-cell engager but the same on-target risk",),
        ("Antibody",),
    )
)
_add(
    MechanismSpec(
        "adc",
        "an antibody carries a cytotoxic payload that is released after internalisation",
        "antibody-like",
        "an accessible extracellular epitope on an internalising receptor",
        "cytotoxic payload",
        True,
        "released payload inside the tumour cell",
        (),
        "established",
        (SURFACE_GATE, INTERNALISATION_GATE),
        (
            ("surface_accessibility", 1.0),
            ("internalisation", 1.3),
            ("lysosomal_trafficking", 0.9),
            ("tumour_selectivity", 1.0),
            ("normal_tissue_safety", 1.1),
            ("tumour_expression", 0.7),
        ),
        (("shedding", 0.8),),
        "probable",
        "probable",
        (
            "the only mechanism here that genuinely needs the receptor to be swallowed",
            "shed antigen absorbs the conjugate before it reaches the cell",
            "payload chemistry is out of scope: GenomeOS states the trafficking requirement, not a linker",
        ),
        ("Antibody drug conjugate",),
    )
)
_add(
    MechanismSpec(
        "targeted_radionuclide",
        "a binder concentrates a radioisotope at the tumour; the radiation does the work",
        "antibody-like or peptide-like",
        "an accessible extracellular epitope",
        "radionuclide",
        True,
        "emitted radiation, including neighbouring cells",
        (),
        "clinical",
        (SURFACE_GATE,),
        (
            ("surface_accessibility", 1.0),
            ("tumour_selectivity", 1.0),
            ("normal_tissue_safety", 1.3),
        ),
        (),
        "probable",
        "probable",
        (
            "does not require internalisation: crossfire reaches antigen-negative neighbours, which "
            "tolerates heterogeneous antigen expression",
            "bone marrow and kidney dose limits dominate the safety analysis",
        ),
        ("Antibody", "Small molecule"),
    )
)
_add(
    MechanismSpec(
        "immune_marker_delivery",
        "tumour-selective delivery makes the cell display a strongly immunogenic surface marker, which "
        "pre-existing or adaptive antibodies then attack",
        "targeted carrier",
        "an internalising, tumour-selective surface target",
        "instructions for an immunogenic surface marker",
        True,
        "endogenous immune system",
        ("pre-existing antibodies", "complement", "NK cells", "macrophages"),
        "experimental",
        (SURFACE_GATE, DELIVERY_GATE),
        (
            ("surface_accessibility", 1.0),
            ("internalisation", 1.2),
            ("tumour_selectivity", 1.5),
            ("normal_tissue_safety", 1.5),
        ),
        (("shedding", 0.6),),
        "hypothetical",
        "hypothetical",
        (
            "experimental: the whole chain from selective delivery to marker display to effective immune "
            "attack has no clinical precedent that GenomeOS can cite",
            "selectivity carries the entire safety argument, because any healthy cell that takes up the "
            "carrier is marked for destruction too",
            "GenomeOS states the biological requirements only; it does not design markers, carriers or "
            "constructs",
        ),
    )
)
_add(
    MechanismSpec(
        "rna_delivery",
        "a tumour-selective carrier delivers RNA that is translated inside the tumour cell",
        "targeted carrier",
        "an internalising, tumour-selective surface target",
        "mRNA",
        True,
        "the translated biological program",
        (),
        "experimental",
        (SURFACE_GATE, DELIVERY_GATE),
        (
            ("surface_accessibility", 1.0),
            ("internalisation", 1.3),
            ("lysosomal_trafficking", 0.6),
            ("tumour_selectivity", 1.4),
            ("normal_tissue_safety", 1.3),
        ),
        (),
        "unknown",
        "unknown",
        (
            "endosomal escape is the unsolved step and is not modelled here",
            "lysosomal routing helps a cleavable payload and hurts RNA, which must escape before it",
        ),
    )
)
_add(
    MechanismSpec(
        "tumour_suppressor_restoration",
        "the inactivated gene's product is restored inside the tumour cell",
        "targeted carrier",
        "an internalising, tumour-selective surface target for delivery",
        "coding RNA or gene",
        True,
        "the restored protein",
        (),
        "experimental",
        (LOF_GATE, DELIVERY_GATE),
        (
            ("internalisation", 1.0),
            ("tumour_selectivity", 1.5),
            ("normal_tissue_safety", 1.2),
        ),
        (),
        "unknown",
        "unknown",
        (
            REPAIR_CAVEAT,
            "restoration must reach every malignant cell to matter, which is a delivery "
            "problem no current carrier solves",
        ),
    )
)
_add(
    MechanismSpec(
        "genome_editing",
        "the tumour allele is corrected or disrupted in place",
        "targeted carrier",
        "an internalising, tumour-selective surface target for delivery",
        "editing machinery",
        True,
        "the edited genome",
        (),
        "experimental",
        (MUTATION_GATE, DELIVERY_GATE),
        (
            ("internalisation", 1.0),
            ("tumour_selectivity", 1.5),
            ("normal_tissue_safety", 1.3),
        ),
        (),
        "unknown",
        "unknown",
        (
            REPAIR_CAVEAT,
            "editing must also reach the nucleus of every malignant cell; surface access is only the "
            "first of several barriers",
            "GenomeOS describes the rationale and the delivery requirement only; it produces no guide "
            "sequences, constructs or vectors",
        ),
    )
)
_add(
    MechanismSpec(
        "tcr_based",
        "a T-cell receptor recognises the mutant peptide presented by the patient's HLA",
        "TCR-like",
        "a mutant peptide in complex with a specific HLA allele",
        "none",
        False,
        "T cell",
        ("CD8+ T cells",),
        "clinical",
        (PEPTIDE_GATE, HLA_GATE),
        (
            ("neoantigen_strength", 1.3),
            ("presentation_confidence", 1.3),
            ("tumour_selectivity", 0.6),
            ("normal_tissue_safety", 0.8),
        ),
        (),
        "probable",
        "probable",
        (
            "reaches intracellular and nuclear proteins, which no surface binder can",
            "must discriminate the mutant peptide from the wild-type peptide in the same groove",
            "cross-reactivity against similar self peptide/HLA complexes is the documented safety risk",
        ),
        ("Cell therapy", "Protein"),
    )
)
_add(
    MechanismSpec(
        "tcr_mimic",
        "an antibody-like binder recognises the mutant peptide/HLA complex",
        "antibody-like",
        "a mutant peptide in complex with a specific HLA allele",
        "none",
        False,
        "endogenous immune system or a recruited effector",
        ("NK cells", "T cells"),
        "preclinical",
        (PEPTIDE_GATE, HLA_GATE),
        (
            ("neoantigen_strength", 1.2),
            ("presentation_confidence", 1.4),
            ("normal_tissue_safety", 0.8),
        ),
        (),
        "probable",
        "probable",
        (
            "peptide/HLA copy number per cell is far lower than a surface antigen's, so potency is the "
            "limiting factor",
        ),
    )
)


SURFACE_MECHANISMS = tuple(m for m, s in MECHANISMS.items() if SURFACE_GATE in s.gates)
PEPTIDE_MECHANISMS = tuple(m for m, s in MECHANISMS.items() if PEPTIDE_GATE in s.gates)
EXPERIMENTAL_MECHANISMS = tuple(m for m, s in MECHANISMS.items() if s.status == "experimental")


# --- the compatibility engine ---------------------------------------------------------


@dataclass(slots=True)
class Inputs:
    """The scores a mechanism is scored against, plus what is unknown."""

    values: dict[str, float | None] = field(default_factory=dict)
    bases: dict[str, str] = field(default_factory=dict)

    def get(self, key: str) -> float | None:
        return self.values.get(key)


def evaluate(candidate: Any, spec: MechanismSpec, inputs: Inputs) -> MechanismFit:
    """Score one mechanism against one candidate, showing every step."""
    fit = MechanismFit(
        mechanism=spec.mechanism,
        binder=Binder(spec.binder_kind, spec.binder_recognises),
        cargo=Cargo(spec.cargo_kind, spec.payload_required),
        effector=Effector(spec.effector_kind, list(spec.immune_components)),
        payload_required=spec.payload_required,
        compatibility=0.0,
        status=spec.status,
        immunogenic_cell_death=spec.immunogenic_cell_death,
        antigen_release=spec.antigen_release,
        notes=list(spec.notes),
    )
    for gate in spec.gates:
        result = gate.test(candidate)
        if result is False or (result is None and gate.unknown_is_failure):
            state = "not established" if result is None else "not met"
            fit.gates_failed.append(f"{gate.key} {state}: {gate.unmet}")
        elif result is None:
            fit.blocking_unknowns.append(f"{gate.key} is unknown; the mechanism is scored provisionally")

    numerator = 0.0
    denominator = 0.0
    for key, weight in spec.weights:
        v = inputs.get(key)
        if v is None:
            fit.blocking_unknowns.append(f"{key} unknown: {inputs.bases.get(key, 'not established')}")
            continue
        numerator += weight * v
        denominator += weight
        fit.contributions.append(
            {"factor": key, "value": round(v, 3), "weight": weight, "direction": "supports"}
        )
    for key, weight in spec.penalties:
        v = inputs.get(key)
        if v is None:
            continue
        numerator += weight * (1.0 - v)
        denominator += weight
        fit.contributions.append(
            {
                "factor": key,
                "value": round(v, 3),
                "weight": weight,
                "direction": "penalises",
                "contributes": round(1.0 - v, 3),
            }
        )

    if fit.gates_failed:
        fit.compatibility = 0.0
        fit.contributions.append(
            {
                "factor": "hard requirement",
                "value": 0.0,
                "weight": 1.0,
                "direction": "gate",
                "contributes": 0.0,
            }
        )
        return fit
    if denominator == 0.0:
        fit.compatibility = 0.0
        fit.blocking_unknowns.append(
            "every weighted input for this mechanism is unknown; compatibility cannot be computed"
        )
        return fit

    raw = numerator / denominator
    factor = STATUS_FACTOR[spec.status]
    declared = sum(w for _, w in spec.weights) + sum(w for _, w in spec.penalties)
    coverage = denominator / declared if declared else 0.0
    damping = COVERAGE_FLOOR + (1.0 - COVERAGE_FLOOR) * coverage
    available = len(spec.weights) + len(spec.penalties) - len(fit.blocking_unknowns)
    fit.compatibility = round(raw * factor * damping, 3)
    fit.contributions.append(
        {
            "factor": "mechanism maturity",
            "value": factor,
            "weight": 0.0,
            "direction": "multiplier",
            "note": f"{spec.status} mechanism class",
        }
    )
    fit.contributions.append(
        {
            "factor": "input coverage",
            "value": round(coverage, 3),
            "weight": 0.0,
            "direction": "multiplier",
            "note": f"{max(available, 0)} of {len(spec.weights) + len(spec.penalties)} declared inputs "
            f"were available; the score is multiplied by {damping:.2f} so a mechanism cannot rank "
            "highly on missing data",
        }
    )
    fit.evidence.append(
        derived(
            "GenomeOS mechanism compatibility",
            f"{spec.mechanism}: weighted mean {raw:.3f} over {len(spec.weights) + len(spec.penalties)} "
            f"declared factors ({coverage:.0%} available), times {factor:.2f} for a {spec.status} "
            f"mechanism class and {damping:.2f} for input coverage",
            round(coverage * 0.8, 3),
        )
    )
    return fit


def precedent_for(
    spec: MechanismSpec, precedent: dict[str, Any] | None, externally_reachable: bool | None = None
) -> tuple[dict[str, Any] | None, list[Evidence]]:
    """Existing therapies of this mechanism's modality against this target.

    Real precedent is clinical evidence that a target is reachable that way.
    It is never evidence that the therapy suits this patient's tumour, and the
    note says so every time.

    One subtlety decides whether the note is evidence or noise. A mechanism
    that needs an extracellular epitope can cite a small-molecule precedent
    only when the target is *known* to be reachable from outside: radioligands
    against a surface enzyme are real, so the modality is not wrong in
    general. Against an intracellular protein the same citation is a
    non-sequitur, because a kinase inhibitor crossing the membrane says
    nothing about whether a binder can reach the protein from the outside.
    Without positive evidence of an outward-facing part, small-molecule
    precedent is dropped for such mechanisms (docs/DECISIONS.md D29).
    """
    if not precedent or not spec.precedent_drug_types:
        return None, []
    types = spec.precedent_drug_types
    needs_epitope = "extracellular" in spec.binder_recognises
    if needs_epitope and externally_reachable is not True:
        types = tuple(t for t in types if t != "Small molecule")
        if not types:
            return None, []
    rows = ((precedent.get("drugAndClinicalCandidates") or {}).get("rows")) or []
    hits = [r for r in rows if (r.get("drug") or {}).get("drugType") in types]
    if not hits:
        return None, []
    approved = [r for r in hits if r.get("maxClinicalStage") == "APPROVAL"]
    trials = [r for r in hits if str(r.get("maxClinicalStage", "")).startswith("PHASE")]
    names = sorted({(r["drug"]["name"] or "").title() for r in approved})[:5]
    payload = {
        "modalities": list(types),
        "approved": len(approved),
        "in_trials": len(trials),
        "examples": names,
        "caveat": (
            "an existing drug against this gene shows the target is reachable by this modality; it does "
            "not mean the drug is appropriate for this patient's tumour, whose alteration, expression "
            "and context may differ entirely"
        ),
    }
    ev: list[Evidence] = []
    if approved:
        from .evidence import clinical

        ev.append(
            clinical(
                "Open Targets Platform / ChEMBL",
                f"{len(approved)} approved {'/'.join(types).lower()} therapies "
                f"already engage this target: {', '.join(names)}",
                0.95,
            )
        )
    elif trials:
        ev.append(
            Evidence(
                "Open Targets Platform / ChEMBL",
                "clinical_trial",
                f"{len(trials)} agents of this modality against this target have entered trials",
                "clinical",
                0.8,
            )
        )
    return payload, ev
