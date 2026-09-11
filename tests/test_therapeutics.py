"""The therapeutic target pipeline, offline.

Every provider is stubbed so the tests state exactly what the reasoning is
given, and assert what it must and must not conclude from it. The cases are
the ones that decide whether the pipeline is honest: a nuclear protein must
lose the surface route and keep the peptide route, a missing measurement must
stay missing, and a mechanism must not score well on ignorance.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from genomeos.cancer.tumour import TumourVariant
from genomeos.therapeutics import (
    MECHANISMS,
    PatientProfile,
    build_candidate,
    dataset,
    design_readiness,
    evidence_graph,
    machine_report,
    negative_set,
    negative_targets,
    positive_set,
    specification,
    text_report,
)
from genomeos.therapeutics import logic as combination_logic
from genomeos.therapeutics.design import NO_CONSTRUCTS
from genomeos.therapeutics.evidence import patient
from genomeos.therapeutics.localisation import localise
from genomeos.therapeutics.neoantigen import normalise_alleles, peptides_around
from genomeos.therapeutics.pipeline import analyse, vaf_table
from genomeos.therapeutics.providers import Answer, NoNeoantigenPredictor

# --- fixtures: compiled protein definitions, the shape the compiler produces ----------


def defn(
    gene: str,
    accession: str,
    sequence: str,
    locations: list[str],
    keywords: list[str],
    features: list[dict] | None = None,
    pdb: list[dict] | None = None,
    partners: list[dict] | None = None,
    ensembl: str = "ENSG00000000000",
) -> dict:
    return {
        "id": f"UniProt:{accession}",
        "gene": gene,
        "sections": {
            "identity": {
                "items": {
                    "accession": accession,
                    "name": f"{gene} test protein",
                    "length": len(sequence),
                    "sequence": sequence,
                    "keywords": keywords,
                },
                "evidence": "curated",
                "confidence": 0.95,
            },
            "function": {"items": {"summary": [], "location": locations}, "confidence": 0.95},
            "domains": {"items": {"interpro": [], "features": features or []}, "confidence": 0.9},
            "modifications": {"items": []},
            "processing": {"items": []},
            "structures_experimental": {"items": pdb or []},
            "structures_predicted": {"items": [{"id": accession, "method": "PREDICTED"}]},
            "pathways": {"items": []},
            "diseases": {"items": []},
            "genomic_origin": {"items": {"gene_id": ensembl, "transcripts": []}},
            "interactions": {"items": partners or []},
            "expression": {"items": None},
        },
    }


NUCLEAR = defn(
    "RUNX1",
    "Q01196",
    "M" * 200 + "Y" + "K" * 252,
    ["Nucleus"],
    ["DNA-binding", "Nucleus", "Proto-oncogene"],
)
RECEPTOR = defn(
    "ERBB2",
    "P04626",
    "A" * 1255,
    ["Cell membrane", "Early endosome"],
    ["Cell membrane", "Transmembrane", "Receptor", "Endocytosis", "Coated pit"],
    features=[
        {"type": "Topological domain", "description": "Extracellular", "start": 23, "end": 652},
        {"type": "Transmembrane", "description": "Helical", "start": 653, "end": 675},
        {"type": "Topological domain", "description": "Cytoplasmic", "start": 676, "end": 1255},
        {"type": "Domain", "description": "Protein kinase", "start": 720, "end": 987},
    ],
    pdb=[{"source": "PDB", "id": "1N8Z", "method": "X_RAY", "resolution_A": 2.5, "chains": "C=23-624"}],
    partners=[{"partner": "EGFR", "score": 0.99, "physical_evidence": True}],
    ensembl="ENSG00000141736",
)
STATIC_RECEPTOR = defn(
    "FLAT1",
    "P00001",
    "A" * 500,
    ["Cell membrane"],
    ["Cell membrane", "Transmembrane"],
    features=[
        {"type": "Topological domain", "description": "Extracellular", "start": 1, "end": 300},
        {"type": "Transmembrane", "description": "Helical", "start": 301, "end": 321},
        {"type": "Topological domain", "description": "Cytoplasmic", "start": 322, "end": 500},
    ],
)
PARTNER = defn(
    "PART1",
    "P00002",
    "A" * 400,
    ["Cell membrane"],
    ["Cell membrane", "Transmembrane"],
    features=[
        {"type": "Topological domain", "description": "Extracellular", "start": 1, "end": 250},
        {"type": "Transmembrane", "description": "Helical", "start": 251, "end": 271},
    ],
)

DEFINITIONS = {d["gene"]: d for d in (NUCLEAR, RECEPTOR, STATIC_RECEPTOR, PARTNER)}


def hpa_row(gene: str, levels: dict[str, float], specificity: str = "Low tissue specificity") -> dict:
    row = {
        "Gene": gene,
        "Ensembl": "ENSG0",
        "RNA tissue specificity": specificity,
        "RNA tissue distribution": "Detected in many",
    }
    for tissue, value in levels.items():
        row[f"Tissue RNA - {tissue} [nTPM]"] = str(value)
    return row


LOW = {"heart muscle": 0.1, "liver": 0.4, "kidney": 0.2, "lung": 1.1, "bone marrow": 0.3}
HIGH_HEART = {"heart muscle": 180.0, "liver": 2.0, "kidney": 1.0, "lung": 3.0, "bone marrow": 0.4}
PARTNER_LEVELS = {"heart muscle": 0.2, "liver": 40.0, "kidney": 0.5, "lung": 0.8, "bone marrow": 0.1}


# --- stub providers -------------------------------------------------------------------


class StubProtein:
    def __init__(self, definitions: dict[str, dict]) -> None:
        self.definitions = definitions

    def annotation(self, gene: str) -> Answer:
        d = self.definitions.get(gene.upper())
        if d is None:
            return Answer.none(f"no reviewed UniProt entry for {gene}")
        return Answer(d, [], True, "", "2026-01-01")


class StubExpression:
    def __init__(self, rows: dict[str, dict]) -> None:
        self.rows = rows

    def normal_expression(self, gene: str) -> Answer:
        row = self.rows.get(gene.upper())
        if row is None:
            return Answer.none("no Human Protein Atlas record for this gene in the stub")
        return Answer(row, [], True, "", "2026-01-01")


class StubStructure:
    def __init__(self, protein: StubProtein) -> None:
        self.protein = protein

    def structures(self, gene: str, annotation=None) -> Answer:
        d = annotation or self.protein.annotation(gene).data
        if d is None:
            return Answer.none("no annotation")
        exp = (d["sections"]["structures_experimental"] or {}).get("items") or []
        pred = (d["sections"]["structures_predicted"] or {}).get("items") or []
        if not exp and not pred:
            return Answer.none("no structure")
        return Answer({"experimental": exp, "predicted": pred}, [], True, "", "2026-01-01")


class StubTrafficking:
    knowledge_base = None

    def __init__(self, protein: StubProtein) -> None:
        self.protein = protein

    def trafficking(self, gene: str, annotation=None) -> Answer:
        d = annotation or self.protein.annotation(gene).data
        if d is None:
            return Answer.none("no annotation")
        ident = d["sections"]["identity"]["items"]
        fn = d["sections"]["function"]["items"]
        return Answer(
            {
                "keywords": [k.lower() for k in ident["keywords"]],
                "locations": [x.lower() for x in fn["location"]],
                "processing": [],
                "modifications": [],
                "go_terms": [],
            },
            [],
            True,
            "",
            "2026-01-01",
        )


class StubPrecedent:
    def __init__(self, payloads: dict[str, dict] | None = None) -> None:
        self.payloads = payloads or {}

    def precedent(self, gene: str, ensembl_gene=None) -> Answer:
        p = self.payloads.get(gene.upper())
        if p is None:
            return Answer.none("no Open Targets record in the stub")
        return Answer(p, [], True, "", "2026-01-01")


class StubHomology:
    def __init__(self, rows: dict[str, list[dict]] | None = None) -> None:
        self.rows = rows or {}

    def paralogues(self, gene: str, ensembl_gene=None) -> Answer:
        r = self.rows.get(gene.upper())
        if r is None:
            return Answer.none("no paralogue record in the stub")
        return Answer(r, [], True, "", "2026-01-01")


class StubRna:
    def __init__(self, values: dict[str, float] | None = None) -> None:
        self.values = values or {}
        self.path = "stub" if values else ""

    def tumour_expression(self, gene: str) -> Answer:
        if not self.values:
            return Answer.none("no tumour RNA-seq supplied")
        v = self.values.get(gene.upper())
        if v is None:
            return Answer.none(f"{gene} absent from the supplied tumour RNA-seq table")
        ev = patient("tumour RNA-seq (stub)", f"{gene} measured at {v:g} TPM in this tumour", 0.9)
        return Answer({"value": v, "unit": "TPM", "sample": "stub"}, [ev], True, "", "2026-01-01")


class StubProviders:
    def __init__(
        self,
        definitions=None,
        expression=None,
        precedent=None,
        paralogues=None,
        rna=None,
        neoantigen=None,
    ) -> None:
        self.protein = StubProtein(definitions or DEFINITIONS)
        self.expression = StubExpression(expression or {})
        self.structure = StubStructure(self.protein)
        self.trafficking = StubTrafficking(self.protein)
        self.precedent = StubPrecedent(precedent)
        self.homology = StubHomology(paralogues)
        self.patient_rna = StubRna(rna)
        self.neoantigen = neoantigen or NoNeoantigenPredictor()

    def versions(self) -> dict[str, str]:
        return {"stub": "offline test providers"}


def _nucleotide_runs(obj, found=None) -> list[str]:
    """Any string that looks like a nucleotide sequence anywhere in a payload."""
    found = [] if found is None else found
    if isinstance(obj, str):
        if re.fullmatch(r"[ACGTUacgtu]{20,}", obj):
            found.append(obj)
    elif isinstance(obj, dict):
        for v in obj.values():
            _nucleotide_runs(v, found)
    elif isinstance(obj, list):
        for v in obj:
            _nucleotide_runs(v, found)
    return found


def variant(
    gene: str,
    consequence: str,
    change: str = "",
    residue: int | None = None,
    pos: int = 1000,
    driver: float | None = 0.05,
    ref: str = "A",
    alt: str = "T",
) -> TumourVariant:
    v = TumourVariant("chr1", pos, ref, alt, gene=gene, consequence=consequence)
    v.protein_change = change
    v.residue = residue
    v.driver_frequency = driver
    v.score = 1.0
    return v


PROFILE = PatientProfile("test", "tumour.vcf")


def candidate(gene, variants, providers=None, profile=PROFILE):
    return build_candidate(gene, variants, profile, providers or StubProviders(), {})


# --- the cases ------------------------------------------------------------------------


def test_nuclear_protein_loses_the_surface_route_and_keeps_the_peptide_route():
    c = candidate("RUNX1", [variant("RUNX1", "missense_variant", "Y201K", 201)])
    assert c.target_class == "neoantigen_hla"
    assert c.localization.primary == "nuclear"
    assert c.scores.value("surface_accessibility") == 0.0
    surface = {m.mechanism: m for m in c.therapeutic_mechanisms}
    for name in ("adcc", "adc", "t_cell_engager", "blocking_antibody"):
        assert surface[name].compatibility == 0.0
        assert surface[name].gates_failed, name
    # the peptide route survives: peptides derived, and the reason is stated
    assert c.neoantigen is not None and c.neoantigen.applicable
    assert c.neoantigen.peptides
    failed = surface["adcc"].gates_failed[0]
    assert "not reachable" in failed or "cannot engage" in failed


def test_a_truncation_is_not_pretended_to_make_a_neoepitope():
    c = candidate("RUNX1", [variant("RUNX1", "stop_gained", "Y201*", 201)])
    assert c.neoantigen.novel_peptide_sequence is False
    assert c.neoantigen.peptides == []
    assert "identical to wild" in c.neoantigen.reason
    assert c.target_class == "intracellular_only"
    tcr = next(m for m in c.therapeutic_mechanisms if m.mechanism == "tcr_based")
    assert tcr.compatibility == 0.0 and tcr.gates_failed


def test_surface_receptor_is_a_direct_candidate_with_its_topology_captured():
    c = candidate("ERBB2", [variant("ERBB2", "missense_variant", "A310V", 310)])
    assert c.target_class == "direct_surface"
    assert c.localization.plasma_membrane is True
    assert c.localization.topology == "single-pass"
    assert [(r.start, r.end) for r in c.localization.extracellular_regions] == [(23, 652)]
    assert [(r.start, r.end) for r in c.localization.transmembrane_regions] == [(653, 675)]
    assert [(r.start, r.end) for r in c.localization.intracellular_regions] == [(676, 1255)]
    assert c.localization.extracellular_residue(310) is True
    assert c.localization.extracellular_residue(900) is False
    assert c.scores.value("surface_accessibility") > 0.8


def test_a_mutation_in_the_cytoplasmic_tail_is_not_a_surface_epitope():
    c = candidate("ERBB2", [variant("ERBB2", "missense_variant", "A900V", 900)])
    ep = c.structure.candidate_epitopes[0]
    assert ep.extracellular is False
    assert ep.epitope_type == "mutation_specific_internal"
    assert ep.surface_accessible is False
    spec = specification(c, None)
    assert spec["recognition_specification"]["recognition_mode"] == "expression_differential"


def test_a_mutation_in_the_extracellular_domain_is_a_mutation_specific_epitope():
    c = candidate("ERBB2", [variant("ERBB2", "missense_variant", "A310V", 310)])
    ep = c.structure.candidate_epitopes[0]
    assert ep.extracellular is True
    assert ep.epitope_type == "mutation_specific_surface"
    # extracellular is not promoted to bindable
    assert ep.predicted_bindable is None
    assert ep.experimentally_validated is None
    assert not ep.structural_accessibility.known
    spec = specification(c, None)
    assert spec["recognition_specification"]["recognition_mode"] == "mutation_specific"


def test_high_expression_in_an_essential_organ_is_penalised():
    safe = candidate(
        "ERBB2",
        [variant("ERBB2", "missense_variant", "A310V", 310)],
        StubProviders(expression={"ERBB2": hpa_row("ERBB2", LOW)}),
    )
    risky = candidate(
        "ERBB2",
        [variant("ERBB2", "missense_variant", "A310V", 310)],
        StubProviders(expression={"ERBB2": hpa_row("ERBB2", HIGH_HEART)}),
    )
    assert safe.normal_tissue.on_target_off_tumour_risk == "low"
    assert risky.normal_tissue.on_target_off_tumour_risk == "high"
    assert risky.scores.value("normal_tissue_safety") < safe.scores.value("normal_tissue_safety")
    assert "heart muscle" in " ".join(risky.normal_tissue.tissues_at_risk)
    assert risky.scores.overall < safe.scores.overall
    engager_safe = next(m for m in safe.therapeutic_mechanisms if m.mechanism == "t_cell_engager")
    engager_risky = next(m for m in risky.therapeutic_mechanisms if m.mechanism == "t_cell_engager")
    assert engager_risky.compatibility < engager_safe.compatibility


def test_missing_expression_is_unknown_and_never_fabricated():
    c = candidate("ERBB2", [variant("ERBB2", "missense_variant", "A310V", 310)])
    assert c.tumour.expression.value is None
    assert c.tumour.expression.qualitative is None
    assert "RNA-seq" in c.tumour.expression.reason
    assert c.scores.components["tumour_expression"].value is None
    assert c.scores.components["normal_tissue_safety"].value is None
    assert any("cannot be established" in x for x in c.limitations)
    # a DNA mutation must never imply the protein is made
    layering = [e.claim for e in c.ledger.items if "RNA evidence unavailable" in e.claim]
    assert layering, "the DNA/RNA/protein evidence layering must be stated"


def test_the_overall_score_cannot_be_high_without_safety_data():
    c = candidate("ERBB2", [variant("ERBB2", "missense_variant", "A310V", 310)])
    assert c.scores.overall <= 0.6
    assert any(a["adjustment"] == "safety cap (unknown)" for a in c.scores.adjustments)
    assert "excluded rather than defaulted" in c.scores.formula


def test_internalising_receptor_supports_a_payload_mechanism():
    inter = candidate(
        "ERBB2",
        [variant("ERBB2", "missense_variant", "A310V", 310)],
        StubProviders(expression={"ERBB2": hpa_row("ERBB2", LOW)}),
    )
    flat = candidate(
        "FLAT1",
        [variant("FLAT1", "missense_variant", "A100V", 100)],
        StubProviders(expression={"FLAT1": hpa_row("FLAT1", LOW)}),
    )
    assert inter.trafficking.internalises is True
    assert inter.trafficking.endosomal is True
    assert flat.trafficking.internalises is None
    adc_inter = next(m for m in inter.therapeutic_mechanisms if m.mechanism == "adc")
    adc_flat = next(m for m in flat.therapeutic_mechanisms if m.mechanism == "adc")
    assert adc_inter.compatibility > adc_flat.compatibility
    assert any("internalisation" in b for b in adc_flat.blocking_unknowns)


def test_a_non_internalising_receptor_keeps_adcc_and_loses_adc():
    flat = candidate(
        "FLAT1",
        [variant("FLAT1", "missense_variant", "A100V", 100)],
        StubProviders(expression={"FLAT1": hpa_row("FLAT1", LOW)}),
    )
    adcc = next(m for m in flat.therapeutic_mechanisms if m.mechanism == "adcc")
    adc = next(m for m in flat.therapeutic_mechanisms if m.mechanism == "adc")
    assert adcc.compatibility > 0 and adcc.viable
    assert adcc.payload_required is False and adcc.cargo.kind == "none"
    assert adcc.effector.components == ["NK cells", "macrophages"]
    assert adc.compatibility < adcc.compatibility
    # delivery mechanisms need internalisation established, not merely not-denied
    rna = next(m for m in flat.therapeutic_mechanisms if m.mechanism == "rna_delivery")
    assert rna.compatibility == 0.0 and rna.gates_failed


def test_a_mechanism_cannot_score_well_on_missing_inputs():
    """Same biology, more unknowns, lower compatibility."""
    known = candidate(
        "FLAT1",
        [variant("FLAT1", "missense_variant", "A100V", 100)],
        StubProviders(expression={"FLAT1": hpa_row("FLAT1", LOW)}, rna={"FLAT1": 120.0}),
    )
    unknown = candidate("FLAT1", [variant("FLAT1", "missense_variant", "A100V", 100)])
    a = next(m for m in known.therapeutic_mechanisms if m.mechanism == "adcc")
    b = next(m for m in unknown.therapeutic_mechanisms if m.mechanism == "adcc")
    assert a.compatibility > b.compatibility
    assert len(b.blocking_unknowns) > len(a.blocking_unknowns)


def test_without_hla_the_neoantigen_route_is_limited_not_asserted():
    c = candidate("RUNX1", [variant("RUNX1", "missense_variant", "Y201K", 201)])
    assert c.neoantigen.hla_alleles == []
    assert c.neoantigen.predicted_binding is None
    assert c.neoantigen.observed_immunopeptidomics is None
    assert any("HLA" in m.what for m in c.neoantigen.missing)
    for name in ("tcr_based", "tcr_mimic"):
        m = next(x for x in c.therapeutic_mechanisms if x.mechanism == name)
        assert m.compatibility == 0.0
        assert any("HLA" in g or "hla" in g for g in m.gates_failed)
    assert any("HLA genotype" in x for x in c.limitations)


def test_with_hla_the_peptides_are_prioritisable_but_presentation_stays_unproven():
    profile = PatientProfile("test", "tumour.vcf", hla_alleles=["A*02:01", "HLA-B0702"])
    c = build_candidate(
        "RUNX1", [variant("RUNX1", "missense_variant", "Y201K", 201)], profile, StubProviders(), {}
    )
    assert c.neoantigen.hla_alleles == ["HLA-A*02:01", "HLA-B*07:02"]
    assert c.neoantigen.predicted_binding is None  # no predictor is wired in, so none is invented
    assert c.scores.value("presentation_confidence") == 0.0
    tcr = next(x for x in c.therapeutic_mechanisms if x.mechanism == "tcr_based")
    assert tcr.viable
    assert tcr.compatibility > 0.0


def test_combination_of_two_antigens_is_supported_and_its_gain_computed():
    providers = StubProviders(
        expression={
            "FLAT1": hpa_row("FLAT1", {"heart muscle": 0.1, "liver": 60.0, "kidney": 0.2}),
            "PART1": hpa_row("PART1", {"heart muscle": 40.0, "liver": 0.3, "kidney": 0.2}),
        }
    )
    a = candidate("FLAT1", [variant("FLAT1", "missense_variant", "A100V", 100)], providers)
    b = candidate("PART1", [variant("PART1", "missense_variant", "A100V", 100)], providers)
    combos = combination_logic.pairs([a, b])
    assert combos, "two surface candidates with healthy-tissue data must yield a combination"
    x = combos[0]
    assert x.operator == "AND" and set(x.targets) == {"FLAT1", "PART1"}
    assert set(x.normal_tissues_excluded) == {"heart muscle", "liver"}
    assert x.selectivity_gain and x.selectivity_gain > 1
    assert any("same cell" in e.claim for e in x.evidence)
    assert combination_logic.describe(x) in ("FLAT1 AND PART1", "PART1 AND FLAT1")


# --- the design dataset ---------------------------------------------------------------


def test_the_negative_set_names_what_must_not_be_bound():
    providers = StubProviders(
        expression={"ERBB2": hpa_row("ERBB2", HIGH_HEART)},
        paralogues={"ERBB2": [{"symbol": "EGFR", "ensembl_gene": "ENSG00000146648", "identity": 50.7}]},
    )
    c = candidate("ERBB2", [variant("ERBB2", "missense_variant", "A310V", 310)], providers)
    negatives = negative_targets(c, providers)
    kinds = {n["kind"] for n in negatives}
    assert "wild_type_protein" in kinds
    assert "human_paralogue" in kinds
    assert "healthy_tissue_expression" in kinds
    assert "polymorphism" in kinds
    assert any(n.get("tissue") == "heart muscle" for n in negatives)
    para = next(n for n in negatives if n["kind"] == "human_paralogue")
    assert para["label"] == "EGFR" and para["sequence_identity_percent"] == 50.7


def test_a_peptide_target_must_avoid_the_wild_type_complex():
    profile = PatientProfile("test", "tumour.vcf", hla_alleles=["HLA-A*02:01"])
    c = build_candidate(
        "RUNX1", [variant("RUNX1", "missense_variant", "Y201K", 201)], profile, StubProviders(), {}
    )
    negatives = negative_targets(c, None)
    kinds = {n["kind"] for n in negatives}
    assert "wild_type_peptide_hla" in kinds
    assert "self_peptide_hla" in kinds


def test_the_design_dataset_keeps_the_dna_origin_and_refuses_to_build_anything():
    providers = StubProviders(expression={"ERBB2": hpa_row("ERBB2", LOW)})
    c = candidate("ERBB2", [variant("ERBB2", "missense_variant", "A310V", 310, pos=39724740)], providers)
    a = analyse([], PROFILE, providers)
    a["candidates"] = [c]
    d = dataset(a, providers)
    assert d["schema"] == "GenomeOS.TherapeuticDesignDataset"
    assert d["schema_version"] == "1.0"
    assert d["provenance"]["pipeline_version"]
    assert d["provenance"]["genomeos_version"]
    spec = d["target_specifications"][0]
    origin = spec["origin"][0]
    assert origin["chromosome"] == "chr1" and origin["position"] == 39724740
    assert origin["reference"] == "A" and origin["alternate"] == "T"
    assert origin["protein_change"] == "A310V" and origin["gene"] == "ERBB2"
    # nothing laboratory-ready anywhere in the dataset, the disclaimers aside
    assert "no construct" in d["not_provided"].lower()
    text = json.dumps(d).lower().replace(NO_CONSTRUCTS.lower(), "")
    for forbidden in ("plasmid", "expression cassette", "viral vector", "transfection", "guide rna"):
        assert forbidden not in text, forbidden
    assert not _nucleotide_runs(d), "the dataset must contain no therapeutic nucleotide sequence"
    p, n = positive_set(d), negative_set(d)
    assert p["targets"][0]["gene"] == "ERBB2"
    assert n["targets"][0]["avoid"]


def test_design_readiness_cannot_be_high_without_safety_or_expression_data():
    bare = candidate("ERBB2", [variant("ERBB2", "missense_variant", "A310V", 310)])
    r = design_readiness(bare)
    assert r["overall"] <= 0.5
    assert any("healthy-tissue" in b for b in r["blocking_unknowns"])
    assert any("surface protein abundance" in b for b in r["blocking_unknowns"])
    assert r["components"]["normal_tissue_characterisation"] is None
    with_safety = candidate(
        "ERBB2",
        [variant("ERBB2", "missense_variant", "A310V", 310)],
        StubProviders(expression={"ERBB2": hpa_row("ERBB2", LOW)}),
    )
    r2 = design_readiness(with_safety)
    assert r2["overall"] > r["overall"]
    assert r2["overall"] <= 0.7  # tumour expression is still unmeasured


def test_the_evidence_graph_links_a_mechanism_back_to_the_dna():
    providers = StubProviders(expression={"ERBB2": hpa_row("ERBB2", LOW)})
    c = candidate("ERBB2", [variant("ERBB2", "missense_variant", "A310V", 310)], providers)
    a = analyse([], PROFILE, providers)
    a["candidates"] = [c]
    g = evidence_graph(a)
    kinds = {n["kind"] for n in g["nodes"]}
    assert {"variant", "gene", "protein", "localisation", "target", "mechanism"} <= kinds
    relations = {e["relation"] for e in g["edges"]}
    assert {"observed_in", "alters", "encodes", "localises_to", "compatible_with"} <= relations


def test_evidence_never_mixes_computation_with_clinical_observation():
    c = candidate("ERBB2", [variant("ERBB2", "missense_variant", "A310V", 310)])
    by_level = c.ledger.by_level()
    assert by_level, "a candidate must carry evidence"
    for level, items in by_level.items():
        for e in items:
            assert e.level == level
            assert e.source_type in (
                "database",
                "publication",
                "clinical_trial",
                "patient_data",
                "prediction",
                "derived",
            )
    assert "clinical" not in by_level  # nothing here is a clinical observation


def test_patient_rna_overrides_the_population_assumption():
    providers = StubProviders(expression={"ERBB2": hpa_row("ERBB2", LOW)}, rna={"ERBB2": 240.0})
    c = candidate("ERBB2", [variant("ERBB2", "missense_variant", "A310V", 310)], providers)
    assert c.tumour.expression.value == 240.0
    assert c.tumour.expression.patient_specific is True
    assert c.scores.value("tumour_expression") is not None
    assert c.scores.value("tumour_selectivity") is not None
    assert "240 TPM" in c.scores.components["tumour_selectivity"].basis
    assert any(e.patient_derived for e in c.tumour.evidence)


def test_clonality_comes_from_the_vcf_or_stays_unknown(tmp_path):
    vcf = tmp_path / "t.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        "chr1\t1000\t.\tA\tT\t.\tPASS\tAF=0.48\tGT\t0/1\n"
        "chr1\t2000\t.\tC\tG\t.\tPASS\t.\tGT:AD\t0/1:90,10\n"
    )
    table = vaf_table(str(vcf))
    assert table["chr1:1000:A:T"] == 0.48
    assert table["chr1:2000:C:G"] == pytest.approx(0.1)
    profile = PatientProfile("t", str(vcf))
    providers = StubProviders(expression={"ERBB2": hpa_row("ERBB2", LOW)})
    clonal = build_candidate(
        "ERBB2", [variant("ERBB2", "missense_variant", "A310V", 310, pos=1000)], profile, providers, table
    )
    sub = build_candidate(
        "FLAT1",
        [variant("FLAT1", "missense_variant", "A100V", 100, pos=2000, ref="C", alt="G")],
        profile,
        providers,
        table,
    )
    none = candidate("ERBB2", [variant("ERBB2", "missense_variant", "A310V", 310)], providers)
    assert clonal.tumour.clonality == "clonal" and clonal.scores.value("clonality") == 1.0
    assert sub.tumour.clonality == "subclonal" and sub.scores.value("clonality") < 1.0
    assert none.tumour.clonality == "unknown" and none.scores.value("clonality") is None


def test_clinical_precedent_is_evidence_about_the_target_not_about_this_patient():
    precedent = {
        "tractability": [{"label": "Approved Drug", "modality": "AB", "value": True}],
        "drugAndClinicalCandidates": {
            "count": 2,
            "rows": [
                {"maxClinicalStage": "APPROVAL", "drug": {"name": "TRASTUZUMAB", "drugType": "Antibody"}},
                {
                    "maxClinicalStage": "APPROVAL",
                    "drug": {"name": "TRASTUZUMAB EMTANSINE", "drugType": "Antibody drug conjugate"},
                },
            ],
        },
    }
    providers = StubProviders(expression={"ERBB2": hpa_row("ERBB2", LOW)}, precedent={"ERBB2": precedent})
    c = candidate("ERBB2", [variant("ERBB2", "missense_variant", "A310V", 310)], providers)
    assert c.scores.value("clinical_precedent") == 1.0
    adc = next(m for m in c.therapeutic_mechanisms if m.mechanism == "adc")
    assert adc.precedent and adc.precedent["approved"] == 1
    assert "does not mean the drug is appropriate for this patient" in adc.precedent["caveat"]
    assert any(e.level == "clinical" for e in adc.evidence)


def test_repair_mechanisms_carry_the_multi_driver_caveat():
    c = candidate("RUNX1", [variant("RUNX1", "stop_gained", "Y201*", 201)])
    for name in ("genome_editing", "tumour_suppressor_restoration"):
        m = next(x for x in c.therapeutic_mechanisms if x.mechanism == name)
        assert m.status == "experimental"
        assert any("does not restore a normal cell" in n for n in m.notes)
    assert any("does not restore a normal cell" in x for x in c.limitations)


def test_immune_marker_delivery_exists_and_is_experimental():
    spec = MECHANISMS["immune_marker_delivery"]
    assert spec.status == "experimental"
    assert spec.payload_required is True
    assert "endogenous immune system" in spec.effector_kind
    assert any("no clinical precedent" in n for n in spec.notes)
    assert any("does not design markers, carriers or constructs" in n for n in spec.notes)


def test_payloadless_mechanisms_name_their_effector():
    for name in ("adcc", "t_cell_engager", "nk_cell_engager", "adcp"):
        spec = MECHANISMS[name]
        assert spec.payload_required is False
        assert spec.cargo_kind == "none"
        assert spec.immune_components or spec.effector_kind


def test_peptides_span_the_mutation_and_keep_their_wild_type_counterpart():
    seq = "ABCDEFGHIJKLMNOPQRST"
    peptides = peptides_around(seq, 10, "X", lengths=(9,))
    assert all("X" in p.sequence for p in peptides)
    assert all(p.length == 9 for p in peptides)
    assert all(p.sequence[p.mutation_offset - 1] == "X" for p in peptides)
    assert all(p.wild_type_sequence[p.mutation_offset - 1] == "J" for p in peptides)


def test_hla_allele_parsing():
    good, bad = normalise_alleles(["A*02:01", "HLA-B0702", "hla-c*07:01", "nonsense"])
    assert good == ["HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:01"]
    assert bad == ["HLA-NONSENSE"]


def test_localisation_prefers_curated_topology_over_a_weak_keyword():
    loc = localise("ERBB2", RECEPTOR)
    assert loc.compartments["transmembrane"] >= 0.9
    assert loc.confidence >= 0.9
    claims = " ".join(e.claim for e in loc.evidence)
    assert "outside the cell" in claims
    nuclear = localise("RUNX1", NUCLEAR)
    assert nuclear.reachable is False
    assert nuclear.extracellular_residue(100) is None  # no topology, so no claim either way


def test_the_report_renders_and_carries_the_disclaimer():
    providers = StubProviders(expression={"ERBB2": hpa_row("ERBB2", LOW)})
    c = candidate("ERBB2", [variant("ERBB2", "missense_variant", "A310V", 310)], providers)
    n = candidate("RUNX1", [variant("RUNX1", "stop_gained", "Y201*", 201)], providers)
    a = analyse([], PROFILE, providers)
    a["candidates"] = [c, n]
    text = text_report(a)
    assert "THERAPEUTIC TARGET ANALYSIS" in text
    assert "not clinical recommendations" in text
    assert "computational compatibility score" in text.lower()
    for section in (
        "TUMOUR SUMMARY",
        "TOP DIRECT SURFACE TARGETS",
        "TOP NEOANTIGEN / HLA CANDIDATES",
        "NORMAL-TISSUE SAFETY CONCERNS",
        "TRAFFICKING / INTERNALISATION",
        "THERAPEUTIC MECHANISM RECOMMENDATIONS",
        "EXPERIMENTAL MECHANISMS",
        "REQUIRED MISSING DATA",
        "LIMITATIONS",
    ):
        assert section in text, section
    m = machine_report(a)
    assert m["therapeutic_candidates"][0]["gene"] in ("ERBB2", "RUNX1")
    assert "not clinical recommendations" in m["disclaimer"]
    assert m["data_level"] == 1


def test_the_pipeline_never_claims_a_cure():
    providers = StubProviders(expression={"ERBB2": hpa_row("ERBB2", LOW)})
    c = candidate("ERBB2", [variant("ERBB2", "missense_variant", "A310V", 310)], providers)
    a = analyse([], PROFILE, providers)
    a["candidates"] = [c]
    blob = (text_report(a) + json.dumps(dataset(a, providers))).lower()
    for phrase in ("will cure", "cures the", "guaranteed", "clinically recommended"):
        assert phrase not in blob


# --- accessibility: the distinction the whole pipeline rests on ------------------------

INNER_LEAFLET = defn(
    "SRCLIKE",
    "P00003",
    "A" * 536,
    ["Cell membrane", "Nucleus", "Cytoplasm, cytoskeleton"],
    ["Cell membrane", "Membrane", "Myristate", "Lipoprotein", "Nucleus", "Cytoplasm"],
)
INNER_LEAFLET["sections"]["modifications"]["items"] = [
    {"type": "Lipidation", "description": "N-myristoyl glycine", "start": 2, "end": 2}
]
JUNCTION = defn(
    "TIGHT1",
    "P00004",
    "A" * 500,
    ["Cytoplasm", "Nucleus", "Cell junction, tight junction", "Cell membrane"],
    ["Cell membrane", "Membrane", "Nucleus", "Cytoplasm"],
)
NUCLEAR_RECEPTOR = defn(
    "NR1",
    "P00005",
    "A" * 462,
    ["Nucleus", "Cytoplasm"],
    ["Receptor", "Nucleus", "DNA-binding"],
)
GPI = defn(
    "GPI1",
    "P00006",
    "A" * 300,
    ["Cell membrane"],
    ["Cell membrane", "GPI-anchor", "Lipoprotein"],
)
GPI["sections"]["modifications"]["items"] = [
    {"type": "Lipidation", "description": "GPI-anchor amidated serine", "start": 280, "end": 280}
]


def test_a_lipid_anchored_protein_is_not_a_surface_target():
    """SRC is annotated 'Cell membrane' and faces the cytoplasm. The two must not be confused."""
    loc = localise("SRCLIKE", INNER_LEAFLET)
    assert loc.reachable is False
    assert loc.plasma_membrane is False
    assert loc.orientation == "cytoplasmic face (lipid anchor)"
    assert any("cytoplasmic face" in e.claim for e in loc.evidence)
    c = candidate(
        "SRCLIKE",
        [variant("SRCLIKE", "missense_variant", "A100V", 100)],
        StubProviders(definitions={**DEFINITIONS, "SRCLIKE": INNER_LEAFLET}),
    )
    assert c.target_class != "direct_surface"
    assert next(m for m in c.therapeutic_mechanisms if m.mechanism == "adcc").compatibility == 0.0


def test_a_membrane_annotation_without_a_stated_side_does_not_claim_accessibility():
    loc = localise("TIGHT1", JUNCTION)
    assert loc.reachable is False
    assert loc.orientation == "membrane side not established"
    assert any("not established" in e.claim for e in loc.evidence)


def test_the_receptor_keyword_alone_does_not_make_a_nuclear_receptor_accessible():
    loc = localise("NR1", NUCLEAR_RECEPTOR)
    assert loc.primary == "nuclear"
    assert loc.reachable is False


def test_a_gpi_anchored_protein_is_accessible():
    loc = localise("GPI1", GPI)
    assert loc.reachable is True
    assert loc.compartments["gpi_anchored"] >= 0.9
    assert any("GPI anchor" in e.claim for e in loc.evidence)


def test_an_unreachable_target_is_not_design_ready():
    providers = StubProviders(expression={"RUNX1": hpa_row("RUNX1", LOW)})
    c = candidate("RUNX1", [variant("RUNX1", "stop_gained", "Y201*", 201)], providers)
    assert c.best_mechanism is None
    r = design_readiness(c)
    assert r["overall"] <= 0.15
    assert any(a["adjustment"] == "no viable mechanism" for a in r["adjustments"])
    spec = specification(c, None)
    assert spec["recognition_specification"]["recognition_mode"] == "not_addressable"
    assert "no molecular address" in spec["recognition_specification"]["recognition_rationale"]


def test_a_truncation_flags_nonsense_mediated_decay():
    c = candidate("RUNX1", [variant("RUNX1", "stop_gained", "Y201*", 201)])
    assert any("nonsense-mediated decay" in x for x in c.limitations)


def test_no_paralogues_is_reported_differently_from_a_failed_lookup():
    with_none = StubProviders(paralogues={"FLAT1": []})
    c = candidate("FLAT1", [variant("FLAT1", "missense_variant", "A100V", 100)], with_none)
    entry = next(n for n in negative_targets(c, with_none) if n["kind"] == "human_paralogue")
    assert entry["status"] == "none_found"
    assert "smaller concern" in entry["reason"]

    failed = StubProviders()
    c2 = candidate("FLAT1", [variant("FLAT1", "missense_variant", "A100V", 100)], failed)
    entry2 = next(n for n in negative_targets(c2, failed) if n["kind"] == "human_paralogue")
    assert entry2["status"] == "unavailable"
    assert "uncharacterised" in entry2["reason"]


def test_an_unreachable_target_declares_no_positive_recognition_state():
    c = candidate("RUNX1", [variant("RUNX1", "stop_gained", "Y201*", 201)])
    spec = specification(c, None)
    states = spec["positive_targets"]["states"]
    assert len(states) == 1 and states[0].startswith("none:")
    assert "no molecular address" in states[0]


# --- the demo tumour, end to end ------------------------------------------------------

ROOT = Path(__file__).resolve().parent.parent
DEMO_VCF = ROOT / "data/demo/cancer_tumour.vcf"
VEP_CACHE = ROOT / "data/knowledge/vep/vep_cache.jsonl"
PROTEINS = ROOT / "data/knowledge/proteins"


@pytest.mark.skipif(
    not (DEMO_VCF.exists() and VEP_CACHE.exists() and (PROTEINS / "RUNX1.json").exists()),
    reason="needs the cached VEP annotations and compiled proteins",
)
def test_the_demo_tumour_runs_offline_and_keeps_every_route_open():
    from genomeos.therapeutics import analyse_vcf

    a = analyse_vcf(str(DEMO_VCF), hla=["HLA-A*02:01"], top_genes=6, net=False, indirect=False)
    genes = {c.gene: c for c in a["candidates"]}
    assert {"RUNX1", "APP", "SOD1"} <= set(genes)

    # the nuclear driver keeps its DNA origin and loses only the surface route
    runx1 = genes["RUNX1"]
    assert runx1.target_class == "intracellular_only"
    assert runx1.origins[0].chromosome == "chr21" and runx1.origins[0].position == 34792138
    assert runx1.origins[0].protein_change == "Y480*"
    assert runx1.neoantigen is not None  # the route was evaluated, not skipped
    assert runx1.neoantigen.novel_peptide_sequence is False

    # the membrane protein is reachable, but its mutation faces the cytoplasm
    app = genes["APP"]
    assert app.target_class == "direct_surface"
    assert app.localization.topology == "single-pass"
    assert app.localization.extracellular_residue(770) is False
    epitope = app.structure.candidate_epitopes[0]
    assert epitope.epitope_type == "mutation_specific_internal"
    assert app.neoantigen.peptides  # the missense still yields mutant peptides

    assert a["data_level"]["level_reached"] == 1
    assert any("RNA-seq" in m["input"] for m in a["missing_data"])


@pytest.mark.skipif(
    not (DEMO_VCF.exists() and VEP_CACHE.exists() and (PROTEINS / "RUNX1.json").exists()),
    reason="needs the cached VEP annotations and compiled proteins",
)
def test_the_web_endpoint_returns_the_analysis(tmp_path):
    from genomeos.web.server import Api

    out = Api(ROOT).therapeutics("data/demo/cancer_tumour.vcf", hla="HLA-A*02:01", top=4, offline=True)
    assert out["therapeutic_candidates"]
    assert "not clinical recommendations" in out["disclaimer"]
    assert "THERAPEUTIC TARGET ANALYSIS" in out["report"]
    assert out["schema"] == "GenomeOS.TherapeuticDesignDataset"
    gene = out["therapeutic_candidates"][0]["gene"]
    assert "readiness" in out["design"][gene]
    assert "negative_targets" in out["design"][gene]
