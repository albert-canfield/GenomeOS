"""The Therapeutic Design Dataset.

The output a downstream binder-design system needs, so it never has to redo
the tumour interpretation. It says what must be recognised, and with equal
weight what must not be:

    BIND      the tumour state of this target
    DO NOT    the wild-type form, the closest human paralogues, the same
              protein on healthy tissue, similar exposed motifs

The negative set matters as much as the positive one, because designing a
binder is a selectivity problem, not an affinity problem. Optimising affinity
against a positive target alone produces something that also sticks to the
patient's heart.

What this module does not do, deliberately: it produces no nucleotide or
amino-acid sequence for a therapeutic, no construct, no vector, no expression
cassette, no delivery formulation and no protocol. It describes the molecular
recognition that is required. Deciding what molecule to build is a separate,
downstream, computational stage.
"""

from __future__ import annotations

import time
from typing import Any

import genomeos

from .evidence import derived
from .expression import CONCERN_NTPM, target_density
from .mechanisms import MECHANISMS
from .model import TherapeuticTargetCandidate
from .providers import Providers
from .safety import exclusions

SCHEMA = "GenomeOS.TherapeuticDesignDataset"
SCHEMA_VERSION = "1.0"
PIPELINE_VERSION = "therapeutics/1.0"

NO_CONSTRUCTS = (
    "This dataset defines molecular recognition requirements only. It contains no therapeutic "
    "nucleotide or protein sequence, no construct, no vector, no expression cassette, no delivery "
    "formulation and no laboratory protocol. Producing those is a separate downstream stage with its "
    "own review."
)

#: Modality classification only: whether a candidate is conceptually compatible
#: with a therapeutic that a cell could express, never how to encode one.
ENCODING_STRATEGIES = (
    "protein_administration",
    "dna_or_rna_encoded_binder",
    "cell_expressed_binder",
    "secreted_binder",
    "membrane_bound_therapeutic_receptor",
    "intracellularly_expressed_therapeutic_protein",
)

#: A target nothing can currently be aimed at is not design-ready, whatever
#: else is known about it.
NOT_ADDRESSABLE_CAP = 0.15

READINESS_WEIGHTS = {
    "target_identity": 1.0,
    "tumour_specificity": 1.2,
    "surface_confirmation": 1.2,
    "structural_information": 0.8,
    "normal_tissue_characterisation": 1.2,
    "trafficking_information": 0.6,
}


# --- recognition specification --------------------------------------------------------


def recognition_mode(c: TherapeuticTargetCandidate) -> tuple[str, str]:
    """What actually distinguishes the tumour cell at this target."""
    exposed = [
        e
        for e in c.structure.candidate_epitopes
        if e.extracellular and e.epitope_type == "mutation_specific_surface"
    ]
    if exposed:
        return "mutation_specific", (
            f"the alteration changes residue {exposed[0].protein_position}, which curated topology "
            "places outside the cell: the tumour carries a surface feature healthy cells do not"
        )
    if c.target_class == "neoantigen_hla" or (
        c.neoantigen and c.neoantigen.peptides and not c.localization.reachable
    ):
        return "peptide_hla", (
            "the protein is not reachable from outside, but the alteration yields peptides that HLA "
            "could present; the target is the peptide/HLA complex, not the protein"
        )
    if c.target_class == "pathway_induced_surface":
        return "pathway_induced", (
            "the target is not altered; it is proposed as a surface consequence of a disrupted driver "
            "and requires expression evidence before it is real"
        )
    if c.target_logic and c.target_logic.operator != "SINGLE":
        return "multi_antigen", c.target_logic.rationale
    if not c.localization.reachable:
        return "not_addressable", (
            "the protein is not reachable from outside the cell and this alteration yields no "
            "mutation-derived peptide, so there is no molecular address for a binder to recognise; "
            "the routes left are indirect ones through the pathway this gene disrupts"
        )
    return "expression_differential", (
        "the tumour and healthy cells carry the same protein sequence here; selectivity has to come "
        "from how much of it each displays, which is why tumour and normal expression data dominate "
        "this candidate's uncertainty"
    )


def target_region(c: TherapeuticTargetCandidate) -> dict[str, Any] | None:
    """The part of the protein a binder has to engage."""
    if not c.localization.extracellular_regions:
        return None
    r = max(c.localization.extracellular_regions, key=lambda x: (x.end or 0) - (x.start or 0))
    return {
        "domain": c.structure.target_domain,
        "start": r.start,
        "end": r.end,
        "extracellular": True,
        "topology": c.localization.topology,
        "orientation": c.localization.orientation,
        "structural_accessibility": c.structure.surface_accessibility.to_dict(),
        "structurally_resolved": c.structure.extracellular_structure_available,
    }


def recognition_specification(c: TherapeuticTargetCandidate) -> dict[str, Any]:
    mode, why = recognition_mode(c)
    lead = c.origins[0] if c.origins else None
    return {
        "target_id": f"{c.gene}:{c.uniprot or 'no-accession'}",
        "molecular_target": {
            "gene": c.gene,
            "protein": c.protein,
            "uniprot": c.uniprot,
            "isoform": None,
            "mutation": lead.protein_change if lead else None,
        },
        "recognition_mode": mode,
        "recognition_rationale": why,
        "target_region": target_region(c),
        "desired_specificity": {
            "tumour": _positive_states(c),
            "exclude": [n["label"] for n in negative_targets(c, None)],
        },
        # the strongest few, so this block reads on its own; the full ledger,
        # grouped by evidence level, is the specification's `evidence` section
        "evidence": [e.to_dict() for e in sorted(c.ledger.items, key=_by_level)[:10]],
    }


def _by_level(e: Any) -> tuple[int, float]:
    from .evidence import LEVEL_ORDER

    order = LEVEL_ORDER.index(e.level) if e.level in LEVEL_ORDER else len(LEVEL_ORDER)
    return (order, -e.confidence)


def _positive_states(c: TherapeuticTargetCandidate) -> list[str]:
    out = []
    mode, why = recognition_mode(c)
    lead = c.origins[0] if c.origins else None
    if mode == "not_addressable":
        return [f"none: {why}"]
    if mode == "mutation_specific" and lead:
        out.append(f"{c.gene} carrying {lead.protein_change} on the tumour cell surface")
    elif mode == "peptide_hla" and c.neoantigen:
        alleles = ", ".join(c.neoantigen.hla_alleles) or "the patient's class I alleles"
        out.append(f"mutant {c.gene} peptide presented by {alleles}")
    elif mode == "pathway_induced":
        out.append(f"{c.gene} displayed at elevated density on this tumour, if expression confirms it")
    else:
        out.append(f"{c.gene} displayed at tumour-level density")
    if c.target_logic and c.target_logic.operator == "AND":
        out.append("co-display with " + " and ".join(c.target_logic.targets[1:]))
    return out


# --- the negative set -----------------------------------------------------------------


def negative_targets(c: TherapeuticTargetCandidate, providers: Providers | None) -> list[dict[str, Any]]:
    """What the binder must not recognise, with the reason for each entry."""
    out: list[dict[str, Any]] = []
    mode, _ = recognition_mode(c)
    lead = c.origins[0] if c.origins else None

    if mode == "mutation_specific" and lead:
        out.append(
            {
                "label": f"wild-type {c.gene}",
                "kind": "wild_type_protein",
                "protein": c.gene,
                "uniprot": c.uniprot,
                "state": "wild_type",
                "reason": (
                    f"every healthy cell carries {c.gene} with the reference residue at position "
                    f"{lead.residue}; the binder must tell a single-residue difference apart"
                ),
                "discrimination_required": "single residue",
            }
        )
    elif mode == "peptide_hla" and c.neoantigen:
        out.append(
            {
                "label": f"wild-type {c.gene} peptide/HLA complex",
                "kind": "wild_type_peptide_hla",
                "protein": c.gene,
                "state": "wild_type peptide in the same HLA groove",
                "reason": "healthy cells present the wild-type counterpart of the same peptide",
                "discrimination_required": "single residue inside the HLA groove",
            }
        )
        out.append(
            {
                "label": "similar human peptide/HLA complexes",
                "kind": "self_peptide_hla",
                "reason": (
                    "T-cell receptor cross-reactivity against an unrelated self peptide presented by the "
                    "same allele is the documented cause of severe off-target toxicity"
                ),
                "discrimination_required": "proteome-wide similarity screen, not performed here",
            }
        )
    else:
        out.append(
            {
                "label": f"{c.gene} on healthy cells",
                "kind": "same_protein_healthy_tissue",
                "protein": c.gene,
                "uniprot": c.uniprot,
                "state": "identical protein, lower density",
                "reason": (
                    "the tumour and healthy cells carry the same sequence here, so the only available "
                    "discriminator is display density and no sequence-level discrimination exists"
                ),
                "discrimination_required": "density threshold, which requires quantitative surface data",
            }
        )

    for t in c.normal_tissue.tissues:
        if t.value is not None and t.value >= CONCERN_NTPM:
            out.append(
                {
                    "label": f"{c.gene} in healthy {t.tissue}",
                    "kind": "healthy_tissue_expression",
                    "tissue": t.tissue,
                    "value": t.value,
                    "unit": t.unit,
                    "essential_organ": t.critical,
                    "reason": (
                        f"healthy {t.tissue} carries {t.value:g} {t.unit} of this transcript; binding "
                        "there is on-target off-tumour toxicity"
                    ),
                }
            )

    paras, status = _paralogues(c, providers)
    for i, p in enumerate(paras[:6]):
        rank = "the closest human relative of" if i == 0 else "a human paralogue of"
        out.append(
            {
                "label": p["symbol"],
                "kind": "human_paralogue",
                "ensembl_gene": p["ensembl_gene"],
                "sequence_identity_percent": round(p["identity"], 1) if p.get("identity") else None,
                "reason": (
                    f"{p['symbol']} is {rank} {c.gene} at {p['identity']:.0f}% identity; a binder "
                    "raised against a conserved surface may engage it too"
                )
                if p.get("identity")
                else f"{p['symbol']} is a human paralogue of {c.gene}",
            }
        )
    if not paras:
        out.append(
            {
                "label": "human paralogues",
                "kind": "human_paralogue",
                "reason": (
                    f"no human paralogue of {c.gene} reaches 20% sequence identity, so paralogue "
                    "cross-reactivity is a smaller concern for this target than usual"
                    if status == "looked_up"
                    else "the paralogue lookup did not run, so cross-reactivity against close relatives "
                    "is uncharacterised"
                ),
                "status": "none_found" if status == "looked_up" else status,
            }
        )
    out.append(
        {
            "label": "common polymorphic variants of the target",
            "kind": "polymorphism",
            "reason": (
                "population variation at or near the epitope can abolish binding in some patients or "
                "create unintended recognition; not screened here"
            ),
            "status": "not_assessed",
        }
    )
    return out


def _paralogues(
    c: TherapeuticTargetCandidate, providers: Providers | None
) -> tuple[list[dict[str, Any]], str]:
    """Close human relatives, and whether the lookup ran at all.

    An empty list because the protein has no close relatives is a useful fact.
    An empty list because the lookup failed is the opposite of one, so the two
    are returned separately.
    """
    if providers is None:
        return [], "not_looked_up"
    ann = providers.protein.annotation(c.gene)
    ens = None
    if ann.available:
        ens = ((ann.data.get("sections", {}).get("genomic_origin") or {}).get("items") or {}).get("gene_id")
    answer = providers.homology.paralogues(c.gene, ens)
    if not answer.available:
        return [], "unavailable"
    return [p for p in answer.data if (p.get("identity") or 0) >= 20][:8], "looked_up"


# --- discrimination, density, structure -----------------------------------------------


def discrimination(c: TherapeuticTargetCandidate) -> dict[str, Any]:
    """Tumour versus normal, side by side, with unknowns left unknown."""
    return {
        "target": c.gene,
        "tumour": {
            "altered": c.tumour.altered,
            "alterations": c.tumour.alterations,
            "expression": c.tumour.expression.to_dict(),
            "surface_presence": c.tumour.surface_abundance.to_dict(),
            "clonality": c.tumour.clonality,
            "vaf": c.tumour.vaf,
        },
        "normal_tissues": [
            {"tissue": t.tissue, "expression": t.value, "unit": t.unit, "essential_organ": t.critical}
            for t in c.normal_tissue.tissues
        ],
        "objective": ("maximise tumour binding relative to healthy-tissue binding, not absolute affinity"),
        "density": target_density(c.gene, c.tumour, c.normal_tissue),
    }


def structural_dataset(c: TherapeuticTargetCandidate) -> dict[str, Any]:
    lead = c.origins[0] if c.origins else None
    return {
        "uniprot_id": c.uniprot,
        "pdb_ids": [p["id"] for p in c.structure.pdb_ids],
        "pdb_entries": c.structure.pdb_ids,
        "predicted_structure_id": c.structure.predicted_structure_id,
        "target_domain": c.structure.target_domain,
        "residue_range": c.structure.residue_range,
        "mutation_position": lead.residue if lead else None,
        "extracellular_structure_available": c.structure.extracellular_structure_available,
        "surface_accessibility": c.structure.surface_accessibility.to_dict(),
        "candidate_epitopes": [e.to_dict() for e in c.structure.candidate_epitopes],
        "note": "identifiers only; GenomeOS deposits no coordinates and computes no new structure",
    }


# --- binder requirements and desired behaviour ----------------------------------------


def binder_requirements(c: TherapeuticTargetCandidate, providers: Providers | None) -> dict[str, Any]:
    best = c.best_mechanism
    mode, _ = recognition_mode(c)
    spec = MECHANISMS.get(best.mechanism) if best else None
    constraints = [e["constraint"] for e in exclusions(c.gene, c.normal_tissue)["avoid"]]
    if mode == "mutation_specific":
        constraints.append("must not engage the wild-type protein present on every healthy cell")
    if mode == "peptide_hla":
        constraints.append("must not engage the wild-type peptide in the same HLA molecule")
    if c.trafficking.shedding:
        constraints.append(
            "a soluble form of this antigen circulates and will compete for the binder before it "
            "reaches the cell"
        )
    return {
        "recognition": {
            "positive": _positive_states(c),
            "negative": [n["label"] for n in negative_targets(c, providers)],
        },
        "location": "extracellular" if c.localization.reachable else "not accessible from outside",
        "desired_properties": {
            "high_specificity": True,
            "tumour_selectivity": True,
            "wild_type_discrimination": mode in ("mutation_specific", "peptide_hla"),
            "density_discrimination": mode == "expression_differential",
        },
        "mechanistic_requirements": {
            "mechanism": best.mechanism if best else None,
            "internalisation_required": bool(spec and spec.payload_required),
            "immune_effector_recruitment": bool(spec and spec.immune_components),
            "payload_required": bool(spec and spec.payload_required),
        },
        "constraints": constraints,
        "not_provided": NO_CONSTRUCTS,
    }


def desired_action(c: TherapeuticTargetCandidate) -> dict[str, Any]:
    best = c.best_mechanism
    if best is None:
        return {
            "recognize_tumour": True,
            "recognize_normal": False,
            "after_binding": None,
            "effector": None,
            "therapeutic_goal": "undetermined",
            "reason": "no mechanism passed its hard requirements for this candidate",
        }
    spec = MECHANISMS[best.mechanism]
    internalise = spec.payload_required
    return {
        "recognize_tumour": True,
        "recognize_normal": False,
        "after_binding": {
            "remain_surface_bound": not internalise,
            "internalize": internalise,
            "traffic_to_endosome": internalise,
            "traffic_to_lysosome": internalise and spec.mechanism == "adc",
            "recruit_immune_effector": bool(spec.immune_components),
        },
        "effector": spec.effector_kind,
        "immune_components": list(spec.immune_components),
        "cargo": spec.cargo_kind,
        "therapeutic_goal": ("intracellular_delivery" if internalise else "tumour_cell_elimination")
        if spec.mechanism != "blocking_antibody"
        else "pathway_blockade",
    }


def encoding_compatibility(c: TherapeuticTargetCandidate) -> dict[str, Any]:
    """Modality classification only. No encoding is produced."""
    reachable = c.localization.reachable
    mode, _ = recognition_mode(c)
    options = {
        "protein_administration": (
            reachable or mode == "peptide_hla",
            "an administered protein binder can reach any target displayed on the cell surface",
        ),
        "dna_or_rna_encoded_binder": (
            reachable or mode == "peptide_hla",
            "a binder against an extracellular or peptide/HLA target can in principle be produced from "
            "an encoded template rather than administered as protein",
        ),
        "cell_expressed_binder": (
            reachable or mode == "peptide_hla",
            "an engineered cell can display or secrete a binder against this target class",
        ),
        "secreted_binder": (
            reachable or mode == "peptide_hla",
            "the target is outside the cell, so a secreted binder can reach it",
        ),
        "membrane_bound_therapeutic_receptor": (
            reachable or mode == "peptide_hla",
            "surface and peptide/HLA targets are the two classes a membrane-bound therapeutic receptor "
            "can engage",
        ),
        "intracellularly_expressed_therapeutic_protein": (
            not reachable,
            "only relevant where the target is inside the cell and the therapeutic would have to be "
            "produced there, which also requires a delivery route that is not established here",
        ),
    }
    return {
        "compatible": [k for k, (ok, _why) in options.items() if ok],
        "assessment": {k: {"compatible": ok, "why": why} for k, (ok, why) in options.items()},
        "classification_only": NO_CONSTRUCTS,
    }


# --- design readiness -----------------------------------------------------------------


def design_readiness(c: TherapeuticTargetCandidate) -> dict[str, Any]:
    """How ready this target is for binder design, and what blocks it."""
    blocking: list[str] = []
    comp: dict[str, float | None] = {}

    comp["target_identity"] = 1.0 if c.uniprot else 0.0
    if not c.uniprot:
        blocking.append("no reviewed protein entry, so the target itself is not pinned down")

    comp["tumour_specificity"] = c.scores.value("tumour_selectivity")
    if comp["tumour_specificity"] is None:
        blocking.append("tumour-versus-normal specificity is not established")

    surf = c.scores.value("surface_accessibility")
    if surf is None:
        comp["surface_confirmation"] = None
        blocking.append("surface localisation is not established")
    elif c.tumour.surface_abundance.known:
        comp["surface_confirmation"] = surf
    else:
        comp["surface_confirmation"] = min(0.6, surf)
        blocking.append(
            "surface protein abundance has not been measured in this patient; localisation is curated "
            "annotation, not a measurement of this tumour"
        )

    comp["structural_information"] = c.scores.value("structural_bindability")
    if comp["structural_information"] is None:
        blocking.append("no structure covers the region a binder would engage")

    comp["normal_tissue_characterisation"] = 1.0 if c.normal_tissue.known else None
    if not c.normal_tissue.known:
        blocking.append("healthy-tissue expression is unavailable, so safety cannot be characterised")

    if c.trafficking.internalises is None:
        comp["trafficking_information"] = None
        blocking.append("internalisation kinetics unknown; no trafficking data for this receptor")
    else:
        comp["trafficking_information"] = 1.0 if c.trafficking.lysosomal is not None else 0.6
        if c.trafficking.lysosomal is None:
            blocking.append("post-endocytic routing is not characterised")

    scored = {k: v for k, v in comp.items() if v is not None}
    if scored:
        used = sum(READINESS_WEIGHTS[k] for k in scored)
        overall = sum(READINESS_WEIGHTS[k] * v for k, v in scored.items()) / used
        coverage = used / sum(READINESS_WEIGHTS.values())
    else:
        overall, coverage = 0.0, 0.0
    adjustments = []
    if comp["normal_tissue_characterisation"] is None and overall > 0.5:
        adjustments.append(
            {
                "adjustment": "safety information missing",
                "from": round(overall, 3),
                "to": 0.5,
                "reason": "a target cannot be design-ready while its normal-tissue profile is unknown",
            }
        )
        overall = 0.5
    if not c.tumour.expression.known and overall > 0.7:
        adjustments.append(
            {
                "adjustment": "tumour expression missing",
                "from": round(overall, 3),
                "to": 0.7,
                "reason": "no measurement shows this tumour makes the protein",
            }
        )
        overall = 0.7
    viable = [m for m in c.therapeutic_mechanisms if m.viable and m.compatibility > 0]
    if not viable and overall > NOT_ADDRESSABLE_CAP:
        adjustments.append(
            {
                "adjustment": "no viable mechanism",
                "from": round(overall, 3),
                "to": NOT_ADDRESSABLE_CAP,
                "reason": "no therapeutic mechanism's hard requirements are met for this target, so "
                "there is nothing for a binder to be designed against yet",
            }
        )
        overall = NOT_ADDRESSABLE_CAP
        blocking.insert(0, "no therapeutic mechanism passes its hard requirements for this target")
    return {
        "overall": round(overall, 3),
        "components": {k: (None if v is None else round(v, 3)) for k, v in comp.items()},
        "weights": READINESS_WEIGHTS,
        "component_coverage": round(coverage, 3),
        "adjustments": adjustments,
        "blocking_unknowns": blocking,
        "formula": (
            "weighted mean over the components that were available, then the adjustments above; "
            "components with no data are excluded rather than defaulted"
        ),
    }


# --- the dataset ----------------------------------------------------------------------


def specification(c: TherapeuticTargetCandidate, providers: Providers | None = None) -> dict[str, Any]:
    """One complete target specification: the interface to molecular design."""
    negatives = negative_targets(c, providers)
    return {
        "target_id": f"{c.gene}:{c.uniprot or 'no-accession'}",
        "gene": c.gene,
        "protein": c.protein,
        "target_class": c.target_class,
        "target_class_reason": c.class_reason,
        "origin": [o.to_dict() for o in c.origins],
        "recognition_specification": recognition_specification(c),
        "positive_targets": {
            "states": _positive_states(c),
            "region": target_region(c),
            "logic": c.target_logic.to_dict() if c.target_logic else None,
        },
        "negative_targets": negatives,
        "tumour_versus_normal": discrimination(c),
        "normal_tissue_exclusions": exclusions(c.gene, c.normal_tissue),
        "structure": structural_dataset(c),
        "trafficking": c.trafficking.to_dict(),
        "neoantigen": c.neoantigen.to_dict() if c.neoantigen else None,
        "binder_requirements": binder_requirements(c, providers),
        "desired_action": desired_action(c),
        "encoding_strategy_compatibility": encoding_compatibility(c),
        "therapeutic_mechanisms": [m.to_dict() for m in c.therapeutic_mechanisms if m.viable][:6],
        "scores": c.scores.to_dict(),
        "design_readiness": design_readiness(c),
        "evidence": c.ledger.to_dict(),
        "limitations": c.limitations,
    }


def dataset(analysis: dict[str, Any], providers: Providers | None = None) -> dict[str, Any]:
    """The whole design dataset for one tumour."""
    candidates: list[TherapeuticTargetCandidate] = analysis["candidates"]
    specs = [specification(c, providers) for c in candidates]
    return {
        "schema": SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "provenance": {
            "genomeos_version": genomeos.__version__,
            "pipeline_version": PIPELINE_VERSION,
            "generated": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "sample": analysis["sample"],
            "sample_id": analysis["sample_id"],
            "input_hashes": analysis["input_hashes"],
            "providers": analysis["providers"],
            "data_level": analysis["data_level"]["level_reached"],
            "determinism": (
                "given the same inputs, the same cached database versions and the same pipeline "
                "version, this dataset is reproducible; provider caches live under "
                "data/knowledge/therapeutics/"
            ),
        },
        "disclaimer": analysis["disclaimer"],
        "not_provided": NO_CONSTRUCTS,
        "target_specifications": specs,
        "combination_logic": [x.to_dict() for x in analysis.get("combinations", [])],
        "missing_data": analysis["missing_data"],
    }


def positive_set(dataset_obj: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema": f"{SCHEMA}.PositiveTargets",
        "schema_version": SCHEMA_VERSION,
        "sample_id": dataset_obj["provenance"]["sample_id"],
        "targets": [
            {
                "target_id": s["target_id"],
                "gene": s["gene"],
                "target_class": s["target_class"],
                "recognition_mode": s["recognition_specification"]["recognition_mode"],
                "states": s["positive_targets"]["states"],
                "region": s["positive_targets"]["region"],
                "logic": s["positive_targets"]["logic"],
                "origin": s["origin"],
                "design_readiness": s["design_readiness"]["overall"],
            }
            for s in dataset_obj["target_specifications"]
        ],
    }


def negative_set(dataset_obj: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema": f"{SCHEMA}.NegativeTargets",
        "schema_version": SCHEMA_VERSION,
        "sample_id": dataset_obj["provenance"]["sample_id"],
        "purpose": (
            "what a designed binder must not recognise. Optimising affinity against the positive target "
            "without this set produces cross-reactivity, which is the usual reason a selective binder "
            "fails."
        ),
        "targets": [
            {"target_id": s["target_id"], "gene": s["gene"], "avoid": s["negative_targets"]}
            for s in dataset_obj["target_specifications"]
        ],
    }


def evidence_graph(analysis: dict[str, Any]) -> dict[str, Any]:
    """variant -> transcript -> protein -> structure/expression -> target -> mechanism."""
    nodes: list[dict[str, Any]] = []
    edges: list[dict[str, Any]] = []
    seen: set[str] = set()

    def node(nid: str, kind: str, **attrs: Any) -> str:
        if nid not in seen:
            seen.add(nid)
            nodes.append({"id": nid, "kind": kind, **attrs})
        return nid

    def edge(a: str, b: str, rel: str, evidence: str = "", confidence: float | None = None) -> None:
        edges.append({"from": a, "to": b, "relation": rel, "evidence": evidence, "confidence": confidence})

    sample = node(f"sample:{analysis['sample_id']}", "sample", path=analysis["sample"])
    for c in analysis["candidates"]:
        gene = node(f"gene:{c.gene}", "gene", symbol=c.gene)
        protein = node(f"protein:{c.uniprot or c.gene}", "protein", name=c.protein, uniprot=c.uniprot)
        edge(gene, protein, "encodes", "Ensembl / UniProt", 0.95)
        for o in c.origins:
            vid = node(
                f"variant:{o.chromosome}:{o.position}:{o.reference}>{o.alternate}",
                "variant",
                gene=o.gene,
                consequence=o.variant_type,
                protein_change=o.protein_change,
                somatic_status=o.somatic_status,
            )
            edge(sample, vid, "observed_in", "tumour VCF", 0.9)
            edge(vid, gene, "alters", "Ensembl VEP", 0.95)
            if o.protein_change:
                edge(vid, protein, "changes_residue", f"VEP {o.hgvsp or o.protein_change}", 0.9)
        loc = node(f"localisation:{c.gene}", "localisation", compartment=c.localization.primary)
        edge(protein, loc, "localises_to", "UniProt / HPA / GO", c.localization.confidence)
        if c.normal_tissue.known:
            exp = node(f"expression:{c.gene}", "expression", summary=c.normal_tissue.summary)
            edge(protein, exp, "expressed_in_normal_tissue", "Human Protein Atlas", 0.8)
        if c.structure.pdb_ids or c.structure.predicted_structure_id:
            st = node(
                f"structure:{c.gene}",
                "structure",
                pdb=[p["id"] for p in c.structure.pdb_ids][:5],
                predicted=c.structure.predicted_structure_id,
            )
            edge(protein, st, "has_structure", "PDB / AlphaFold", 0.9)
        if c.neoantigen and c.neoantigen.peptides:
            pep = node(f"peptides:{c.gene}", "peptide_set", count=len(c.neoantigen.peptides))
            edge(protein, pep, "yields_mutant_peptides", "derived from UniProt sequence", 0.7)
        target = node(f"target:{c.gene}", "target", target_class=c.target_class)
        edge(loc, target, "supports_target_class", c.class_reason, c.localization.confidence)
        for m in c.therapeutic_mechanisms:
            if not m.viable or m.compatibility <= 0:
                continue
            mid = node(f"mechanism:{m.mechanism}", "mechanism", status=m.status)
            edge(target, mid, "compatible_with", "GenomeOS mechanism engine", m.compatibility)
    return {
        "schema": f"{SCHEMA}.EvidenceGraph",
        "schema_version": SCHEMA_VERSION,
        "sample_id": analysis["sample_id"],
        "nodes": nodes,
        "edges": edges,
        "traversal": (
            "every target is reachable back to the tumour DNA it came from, and forward to the "
            "mechanisms its biology supports"
        ),
        "evidence": derived("GenomeOS evidence graph", f"{len(nodes)} nodes, {len(edges)} edges", 0.0).claim,
    }
