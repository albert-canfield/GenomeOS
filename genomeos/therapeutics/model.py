"""The therapeutic-target domain model.

Five things are kept apart on purpose, because collapsing them is what makes
a pipeline assume every target needs a toxin:

    Target      what distinguishes the cancer cell
    Binder      what recognises the target
    Mechanism   how therapeutic action follows from recognition
    Cargo       optional material the binder carries (often none)
    Effector    what ultimately acts on the cell

A TherapeuticTargetCandidate is one molecular address on a tumour cell with
everything known about reaching it: where the protein sits, whether the
tumour differs from normal tissue there, what happens after a binder attaches,
what structure exists, whether a mutant peptide could be presented instead,
and which mechanisms the biology supports. Anything unestablished is None
with a stated reason; nothing is defaulted into existence.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

from .evidence import Evidence, EvidenceLedger, Missing

# --- classification vocabularies ------------------------------------------------------

TargetClass = Literal[
    "direct_surface",
    "neoantigen_hla",
    "pathway_induced_surface",
    "secreted",
    "intracellular_only",
    "unsuitable",
    "unknown",
]

#: Where a protein sits. Deliberately finer than "membrane or not".
Compartment = Literal[
    "plasma_membrane",
    "cell_surface",
    "transmembrane",
    "gpi_anchored",
    "secreted",
    "extracellular",
    "cytoplasmic",
    "nuclear",
    "mitochondrial",
    "other_intracellular",
    "unknown",
]

COMPARTMENTS: tuple[str, ...] = (
    "plasma_membrane",
    "cell_surface",
    "transmembrane",
    "gpi_anchored",
    "secreted",
    "extracellular",
    "cytoplasmic",
    "nuclear",
    "mitochondrial",
    "other_intracellular",
    "unknown",
)

#: Compartments from which a circulating binder can physically reach the protein.
REACHABLE = ("plasma_membrane", "cell_surface", "transmembrane", "gpi_anchored")

#: Confidence a reachable compartment needs before accessibility is claimed.
REACHABLE_THRESHOLD = 0.6

TherapeuticMechanism = Literal[
    "blocking_antibody",
    "agonist_antibody",
    "adcc",
    "adcp",
    "complement_recruitment",
    "t_cell_engager",
    "nk_cell_engager",
    "adc",
    "targeted_radionuclide",
    "immune_marker_delivery",
    "rna_delivery",
    "tumour_suppressor_restoration",
    "genome_editing",
    "tcr_based",
    "tcr_mimic",
    "unknown",
]

RecognitionMode = Literal[
    "mutation_specific",
    "expression_differential",
    "isoform_specific",
    "peptide_hla",
    "multi_antigen",
    "pathway_induced",
    "not_addressable",
]

#: How mature a mechanism is as a class of therapy, independent of this target.
MechanismStatus = Literal["established", "clinical", "preclinical", "experimental"]

#: Evidence grade for a mechanistic property such as immunogenic cell death.
Grade = Literal["established", "probable", "preclinical", "hypothetical", "unknown"]

RiskLevel = Literal["low", "moderate", "high", "unknown"]

Clonality = Literal["clonal", "subclonal", "unknown"]


# --- measurements ---------------------------------------------------------------------


@dataclass(slots=True)
class Measure:
    """A number with its provenance, or an explicit absence with its reason.

    A qualitative band (high/medium/low/not_detected) is used where the source
    only supports one; inventing a number from a category is not allowed.
    """

    value: float | None = None
    unit: str = ""
    qualitative: str | None = None
    source: str = ""
    level: str = "inferred"
    patient_specific: bool = False
    confidence: float = 0.0
    reason: str = ""

    @property
    def known(self) -> bool:
        return self.value is not None or self.qualitative is not None

    def to_dict(self) -> dict[str, Any]:
        d: dict[str, Any] = {
            "value": self.value,
            "unit": self.unit,
            "qualitative": self.qualitative,
            "source": self.source,
            "level": self.level,
            "patient_specific": self.patient_specific,
            "confidence": round(self.confidence, 3),
        }
        if not self.known:
            d["reason"] = self.reason or "not established"
        return d

    @classmethod
    def unavailable(cls, reason: str) -> Measure:
        return cls(reason=reason)


@dataclass(slots=True)
class Region:
    """A stretch of a protein: a domain, a topological segment, an epitope window."""

    kind: str
    start: int | None = None
    end: int | None = None
    description: str = ""

    def contains(self, residue: int | None) -> bool:
        if residue is None or self.start is None or self.end is None:
            return False
        return self.start <= residue <= self.end

    def to_dict(self) -> dict[str, Any]:
        return {"kind": self.kind, "start": self.start, "end": self.end, "description": self.description}


# --- the pipeline stages' outputs -----------------------------------------------------


@dataclass(slots=True)
class Localisation:
    """Where the protein is, and which parts of it face outwards."""

    compartments: dict[str, float] = field(default_factory=dict)  # compartment -> confidence
    plasma_membrane: bool | None = None
    extracellular_regions: list[Region] = field(default_factory=list)
    transmembrane_regions: list[Region] = field(default_factory=list)
    intracellular_regions: list[Region] = field(default_factory=list)
    signal_peptide: Region | None = None
    topology: str = "unknown"
    orientation: str = "unknown"
    confidence: float = 0.0
    evidence: list[Evidence] = field(default_factory=list)

    @property
    def primary(self) -> str:
        """The compartment that decides reachability, ties broken towards the outside."""
        if not self.compartments:
            return "unknown"
        rank = {c: len(COMPARTMENTS) - i for i, c in enumerate(COMPARTMENTS)}
        return max(self.compartments.items(), key=lambda kv: (kv[1], rank.get(kv[0], 0)))[0]

    @property
    def reachable(self) -> bool:
        """Can a circulating binder engage it? Needs evidence, not a membrane word."""
        return any(self.compartments.get(c, 0.0) >= REACHABLE_THRESHOLD for c in REACHABLE)

    def extracellular_residue(self, residue: int | None) -> bool | None:
        """Is this residue on the outside? None when topology is unknown."""
        if residue is None:
            return None
        if any(r.contains(residue) for r in self.extracellular_regions):
            return True
        if any(r.contains(residue) for r in self.intracellular_regions):
            return False
        if any(r.contains(residue) for r in self.transmembrane_regions):
            return False
        if not (self.extracellular_regions or self.intracellular_regions):
            return None
        return None

    def to_dict(self) -> dict[str, Any]:
        return {
            "primary": self.primary,
            "compartments": {k: round(v, 3) for k, v in self.compartments.items()},
            "plasma_membrane": self.plasma_membrane,
            "extracellular_regions": [r.to_dict() for r in self.extracellular_regions],
            "transmembrane_regions": [r.to_dict() for r in self.transmembrane_regions],
            "intracellular_regions": [r.to_dict() for r in self.intracellular_regions],
            "signal_peptide": self.signal_peptide.to_dict() if self.signal_peptide else None,
            "topology": self.topology,
            "orientation": self.orientation,
            "confidence": round(self.confidence, 3),
            "evidence": [e.to_dict() for e in self.evidence],
        }


@dataclass(slots=True)
class TumourState:
    """What the tumour does with this gene, as far as the supplied data shows."""

    altered: bool = False
    alterations: list[str] = field(default_factory=list)
    expression: Measure = field(default_factory=lambda: Measure.unavailable("no tumour RNA-seq supplied"))
    protein_abundance: Measure = field(
        default_factory=lambda: Measure.unavailable("no tumour proteomics supplied")
    )
    surface_abundance: Measure = field(
        default_factory=lambda: Measure.unavailable("no surface proteomics supplied")
    )
    copy_number: Measure = field(default_factory=lambda: Measure.unavailable("no copy-number data supplied"))
    vaf: float | None = None
    purity: float | None = None
    clonality: Clonality = "unknown"
    clonality_reason: str = "no variant allele fraction or tumour purity supplied"
    evidence: list[Evidence] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "altered": self.altered,
            "alterations": self.alterations,
            "expression": self.expression.to_dict(),
            "protein_abundance": self.protein_abundance.to_dict(),
            "surface_abundance": self.surface_abundance.to_dict(),
            "copy_number": self.copy_number.to_dict(),
            "vaf": self.vaf,
            "purity": self.purity,
            "clonality": self.clonality,
            "clonality_reason": self.clonality_reason if self.clonality == "unknown" else "",
            "evidence": [e.to_dict() for e in self.evidence],
        }


@dataclass(slots=True)
class TissueExpression:
    tissue: str
    value: float | None
    unit: str = "nTPM"
    qualitative: str | None = None
    critical: bool = False
    weight: float = 0.0

    def to_dict(self) -> dict[str, Any]:
        return {
            "tissue": self.tissue,
            "value": self.value,
            "unit": self.unit,
            "qualitative": self.qualitative,
            "critical": self.critical,
            "weight": self.weight,
        }


@dataclass(slots=True)
class NormalTissueProfile:
    """Healthy-tissue expression, and what it costs a therapy to hit it."""

    tissues: list[TissueExpression] = field(default_factory=list)
    tissues_at_risk: list[str] = field(default_factory=list)
    specificity: str | None = None  # population-level category, e.g. "Tissue enhanced"
    distribution: str | None = None
    summary: str = "unknown"
    on_target_off_tumour_risk: RiskLevel = "unknown"
    risk_basis: str = ""
    evidence: list[Evidence] = field(default_factory=list)

    @property
    def known(self) -> bool:
        return any(t.value is not None for t in self.tissues)

    def to_dict(self) -> dict[str, Any]:
        return {
            "expression_summary": self.summary,
            "specificity": self.specificity,
            "distribution": self.distribution,
            "tissues": [t.to_dict() for t in self.tissues],
            "critical_tissues": self.tissues_at_risk,
            "on_target_off_tumour_risk": self.on_target_off_tumour_risk,
            "risk_basis": self.risk_basis,
            "evidence": [e.to_dict() for e in self.evidence],
        }


@dataclass(slots=True)
class Trafficking:
    """What happens after a binder attaches. Unknown stays unknown."""

    internalises: bool | None = None
    internalisation_rate: Measure = field(
        default_factory=lambda: Measure.unavailable("no internalisation kinetics in any wired source")
    )
    endosomal: bool | None = None
    lysosomal: bool | None = None
    recycling: bool | None = None
    membrane_residence: str | None = None
    shedding: bool | None = None
    soluble_antigen: bool | None = None
    reason: str = ""
    evidence: list[Evidence] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "internalises": self.internalises,
            "internalisation_rate": self.internalisation_rate.to_dict(),
            "endosomal": self.endosomal,
            "lysosomal": self.lysosomal,
            "recycling": self.recycling,
            "membrane_residence": self.membrane_residence,
            "shedding": self.shedding,
            "soluble_antigen": self.soluble_antigen,
            "reason": self.reason,
            "evidence": [e.to_dict() for e in self.evidence],
        }


@dataclass(slots=True)
class Epitope:
    """A candidate recognition site, with its accessibility stated separately.

    Extracellular is not the same as surface-accessible, which is not the same
    as structurally exposed, which is not the same as bindable. Each is its own
    field so none of them can be silently promoted to another.
    """

    epitope_type: str
    protein: str = ""
    protein_position: int | None = None
    region: Region | None = None
    extracellular: bool | None = None
    surface_accessible: bool | None = None
    structurally_resolved: bool | None = None
    predicted_bindable: bool | None = None
    experimentally_validated: bool | None = None
    tumour_specific_change: str | None = None
    reference_context: str | None = None
    mutant_context: str | None = None
    structural_accessibility: Measure = field(
        default_factory=lambda: Measure.unavailable("no solvent-accessibility computation performed")
    )
    normal_human_match_risk: str = "unknown"
    evidence: list[Evidence] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "epitope_type": self.epitope_type,
            "protein": self.protein,
            "protein_position": self.protein_position,
            "region": self.region.to_dict() if self.region else None,
            "extracellular": self.extracellular,
            "surface_accessible": self.surface_accessible,
            "structurally_resolved": self.structurally_resolved,
            "predicted_bindable": self.predicted_bindable,
            "experimentally_validated": self.experimentally_validated,
            "tumour_specific_change": self.tumour_specific_change,
            "reference_context": self.reference_context,
            "mutant_context": self.mutant_context,
            "structural_accessibility": self.structural_accessibility.to_dict(),
            "normal_human_match_risk": self.normal_human_match_risk,
            "evidence": [e.to_dict() for e in self.evidence],
        }


@dataclass(slots=True)
class StructureContext:
    uniprot_id: str | None = None
    pdb_ids: list[dict[str, Any]] = field(default_factory=list)
    predicted_structure_id: str | None = None
    extracellular_structure_available: bool | None = None
    target_domain: str | None = None
    residue_range: list[int] | None = None
    mutation_position: int | None = None
    mutation_resolved_in: int | None = None  # deposited structures covering the altered residue
    surface_accessibility: Measure = field(
        default_factory=lambda: Measure.unavailable("not computed; no coordinates loaded")
    )
    candidate_epitopes: list[Epitope] = field(default_factory=list)
    evidence: list[Evidence] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "uniprot_id": self.uniprot_id,
            "pdb_ids": self.pdb_ids,
            "predicted_structure_id": self.predicted_structure_id,
            "extracellular_structure_available": self.extracellular_structure_available,
            "target_domain": self.target_domain,
            "residue_range": self.residue_range,
            "mutation_position": self.mutation_position,
            "mutation_resolved_in_structures": self.mutation_resolved_in,
            "surface_accessibility": self.surface_accessibility.to_dict(),
            "candidate_epitopes": [e.to_dict() for e in self.candidate_epitopes],
            "evidence": [e.to_dict() for e in self.evidence],
        }


@dataclass(slots=True)
class PeptideCandidate:
    sequence: str
    wild_type_sequence: str
    length: int
    start: int
    end: int
    mutation_offset: int  # 1-based position of the changed residue inside the peptide
    source: str = "mutation_derived"

    def to_dict(self) -> dict[str, Any]:
        return {
            "sequence": self.sequence,
            "wild_type_sequence": self.wild_type_sequence,
            "length": self.length,
            "window": [self.start, self.end],
            "mutation_offset": self.mutation_offset,
            "source": self.source,
        }


@dataclass(slots=True)
class NeoantigenAssessment:
    """The peptide/HLA route for a variant, kept strictly layered.

    A predicted binder is not a presented peptide; a presented peptide is not
    a recognised one. Each layer has its own field and an unavailable layer
    stays unavailable.
    """

    applicable: bool = False
    novel_peptide_sequence: bool | None = None
    reason: str = ""
    peptides: list[PeptideCandidate] = field(default_factory=list)
    hla_alleles: list[str] = field(default_factory=list)
    predicted_binding: dict[str, Any] | None = None
    predicted_processing: dict[str, Any] | None = None
    observed_immunopeptidomics: dict[str, Any] | None = None
    clinical_evidence: dict[str, Any] | None = None
    wild_type_discrimination: str = "unknown"
    evidence: list[Evidence] = field(default_factory=list)
    missing: list[Missing] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "applicable": self.applicable,
            "novel_peptide_sequence": self.novel_peptide_sequence,
            "reason": self.reason,
            "peptide_candidates": [p.to_dict() for p in self.peptides],
            "hla_alleles": self.hla_alleles,
            "predicted_binding": self.predicted_binding,
            "predicted_processing": self.predicted_processing,
            "observed_immunopeptidomics": self.observed_immunopeptidomics,
            "clinical_evidence": self.clinical_evidence,
            "wild_type_discrimination": self.wild_type_discrimination,
            "evidence": [e.to_dict() for e in self.evidence],
            "missing": [m.to_dict() for m in self.missing],
        }


# --- target / binder / mechanism / cargo / effector ------------------------------------


@dataclass(slots=True)
class Target:
    """What distinguishes the cancer cell. Never the therapy itself."""

    id: str
    gene: str
    protein: str = ""
    uniprot: str | None = None
    isoform: str | None = None
    mutation: str | None = None
    recognition_mode: RecognitionMode = "expression_differential"
    region: Region | None = None
    description: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "gene": self.gene,
            "protein": self.protein,
            "uniprot": self.uniprot,
            "isoform": self.isoform,
            "mutation": self.mutation,
            "recognition_mode": self.recognition_mode,
            "region": self.region.to_dict() if self.region else None,
            "description": self.description,
        }


@dataclass(slots=True)
class Binder:
    """What recognises the target. A class of molecule, not a sequence."""

    kind: str  # "antibody-like", "TCR-like", "targeted carrier", "none"
    recognises: str = ""
    note: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {"kind": self.kind, "recognises": self.recognises, "note": self.note}


@dataclass(slots=True)
class Cargo:
    """Optional material carried by the binder. `none` is the common case."""

    kind: str = "none"
    required: bool = False
    note: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {"kind": self.kind, "required": self.required, "note": self.note}


@dataclass(slots=True)
class Effector:
    """What ultimately acts on the cell."""

    kind: str
    components: list[str] = field(default_factory=list)
    note: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {"kind": self.kind, "components": self.components, "note": self.note}


@dataclass(slots=True)
class MechanismFit:
    """One complete therapeutic strategy scored against one candidate."""

    mechanism: str
    binder: Binder
    cargo: Cargo
    effector: Effector
    payload_required: bool
    compatibility: float
    status: MechanismStatus
    contributions: list[dict[str, Any]] = field(default_factory=list)
    gates_failed: list[str] = field(default_factory=list)
    blocking_unknowns: list[str] = field(default_factory=list)
    immunogenic_cell_death: Grade = "unknown"
    antigen_release: Grade = "unknown"
    precedent: dict[str, Any] | None = None
    notes: list[str] = field(default_factory=list)
    evidence: list[Evidence] = field(default_factory=list)

    @property
    def viable(self) -> bool:
        return not self.gates_failed

    def to_dict(self) -> dict[str, Any]:
        return {
            "mechanism": self.mechanism,
            "compatibility": round(self.compatibility, 3),
            "label": "computational compatibility score, not a response probability",
            "status": self.status,
            "binder": self.binder.to_dict(),
            "cargo": self.cargo.to_dict(),
            "payload_required": self.payload_required,
            "effector": self.effector.to_dict(),
            "immune_components": self.effector.components,
            "immunogenic_cell_death": self.immunogenic_cell_death,
            "antigen_release": self.antigen_release,
            "contributions": self.contributions,
            "gates_failed": self.gates_failed,
            "blocking_unknowns": self.blocking_unknowns,
            "precedent": self.precedent,
            "notes": self.notes,
            "evidence": [e.to_dict() for e in self.evidence],
        }


# --- combination logic ----------------------------------------------------------------


@dataclass(slots=True)
class TargetLogic:
    """Boolean recognition over more than one antigen.

    `A AND B` can be more selective than either alone when healthy tissues
    carry one but never both.
    """

    operator: str  # "SINGLE" | "AND" | "OR" | "AND_NOT"
    targets: list[str]
    rationale: str = ""
    normal_tissues_excluded: list[str] = field(default_factory=list)
    normal_tissues_still_at_risk: list[str] = field(default_factory=list)
    selectivity_gain: float | None = None
    evidence: list[Evidence] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "operator": self.operator,
            "targets": self.targets,
            "rationale": self.rationale,
            "normal_tissues_excluded": self.normal_tissues_excluded,
            "normal_tissues_still_at_risk": self.normal_tissues_still_at_risk,
            "selectivity_gain": self.selectivity_gain,
            "evidence": [e.to_dict() for e in self.evidence],
        }


# --- scoring --------------------------------------------------------------------------


@dataclass(slots=True)
class ScoreComponent:
    key: str
    value: float | None
    basis: str
    weight: float = 1.0
    unknown_reason: str = ""
    evidence: list[Evidence] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        d: dict[str, Any] = {
            "value": None if self.value is None else round(self.value, 3),
            "weight": self.weight,
            "basis": self.basis,
        }
        if self.value is None:
            d["unknown_reason"] = self.unknown_reason or "not established"
        return d


@dataclass(slots=True)
class ScoreSet:
    components: dict[str, ScoreComponent] = field(default_factory=dict)
    overall: float | None = None
    formula: str = ""
    adjustments: list[dict[str, Any]] = field(default_factory=list)
    coverage: float = 0.0
    confidence: float = 0.0

    def value(self, key: str) -> float | None:
        c = self.components.get(key)
        return c.value if c else None

    def to_dict(self) -> dict[str, Any]:
        return {
            "overall": None if self.overall is None else round(self.overall, 3),
            "components": {k: v.to_dict() for k, v in self.components.items()},
            "formula": self.formula,
            "adjustments": self.adjustments,
            "component_coverage": round(self.coverage, 3),
            "confidence": round(self.confidence, 3),
            "label": "research-priority score, not a clinical recommendation",
        }


# --- variant origin -------------------------------------------------------------------


@dataclass(slots=True)
class VariantOrigin:
    """The tumour DNA a candidate came from. Never dropped on the way up."""

    sample_id: str = ""
    chromosome: str = ""
    position: int | None = None
    reference: str = ""
    alternate: str = ""
    gene: str = ""
    transcript: str | None = None
    protein_change: str | None = None
    hgvsc: str | None = None
    hgvsp: str | None = None
    residue: int | None = None
    variant_type: str = ""
    somatic_status: str = "unknown"
    somatic_basis: str = ""
    vaf: float | None = None
    copy_number: float | None = None
    clonality: Clonality = "unknown"
    cosmic: list[str] = field(default_factory=list)
    gnomad_af: float | None = None
    driver_frequency: float | None = None
    hotspot: bool = False

    def to_dict(self) -> dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "chromosome": self.chromosome,
            "position": self.position,
            "reference": self.reference,
            "alternate": self.alternate,
            "gene": self.gene,
            "transcript": self.transcript,
            "protein_change": self.protein_change,
            "hgvsc": self.hgvsc,
            "hgvsp": self.hgvsp,
            "residue": self.residue,
            "variant_type": self.variant_type,
            "somatic_status": self.somatic_status,
            "somatic_basis": self.somatic_basis,
            "vaf": self.vaf,
            "copy_number": self.copy_number,
            "clonality": self.clonality,
            "cosmic": self.cosmic,
            "gnomad_af": self.gnomad_af,
            "driver_frequency": self.driver_frequency,
            "hotspot": self.hotspot,
        }


# --- the candidate --------------------------------------------------------------------


@dataclass(slots=True)
class TherapeuticTargetCandidate:
    """One molecular address on the tumour cell, with everything known about it."""

    gene: str
    protein: str = ""
    uniprot: str | None = None
    target_class: TargetClass = "unknown"
    class_reason: str = ""
    origins: list[VariantOrigin] = field(default_factory=list)
    localization: Localisation = field(default_factory=Localisation)
    tumour: TumourState = field(default_factory=TumourState)
    normal_tissue: NormalTissueProfile = field(default_factory=NormalTissueProfile)
    trafficking: Trafficking = field(default_factory=Trafficking)
    structure: StructureContext = field(default_factory=StructureContext)
    neoantigen: NeoantigenAssessment | None = None
    target_logic: TargetLogic | None = None
    therapeutic_mechanisms: list[MechanismFit] = field(default_factory=list)
    scores: ScoreSet = field(default_factory=ScoreSet)
    precedent: dict[str, Any] | None = None
    ledger: EvidenceLedger = field(default_factory=EvidenceLedger)
    limitations: list[str] = field(default_factory=list)
    why_interesting: str = ""

    @property
    def best_mechanism(self) -> MechanismFit | None:
        viable = [m for m in self.therapeutic_mechanisms if m.viable]
        return max(viable, key=lambda m: m.compatibility) if viable else None

    def to_dict(self) -> dict[str, Any]:
        return {
            "gene": self.gene,
            "protein": self.protein,
            "uniprot": self.uniprot,
            "target_class": self.target_class,
            "class_reason": self.class_reason,
            "why_interesting": self.why_interesting,
            "origin_variants": [o.to_dict() for o in self.origins],
            "localization": self.localization.to_dict(),
            "tumour": self.tumour.to_dict(),
            "normal_tissue": self.normal_tissue.to_dict(),
            "trafficking": self.trafficking.to_dict(),
            "structure": self.structure.to_dict(),
            "neoantigen": self.neoantigen.to_dict() if self.neoantigen else None,
            "target_logic": self.target_logic.to_dict() if self.target_logic else None,
            "therapeutic_mechanisms": [m.to_dict() for m in self.therapeutic_mechanisms],
            "scores": self.scores.to_dict(),
            "therapeutic_precedent": self.precedent,
            "evidence": self.ledger.to_dict(),
            "limitations": self.limitations,
        }
