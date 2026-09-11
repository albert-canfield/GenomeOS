"""The therapeutic target reasoning pipeline.

    somatic variants
        -> gene and protein consequence
        -> protein localisation
        -> tumour accessibility
        -> expression and selectivity
        -> internalisation and trafficking
        -> neoantigen / HLA
        -> pathway-induced surface phenotype
        -> normal-tissue exclusion
        -> structural and epitope evidence
        -> therapeutic mechanism selection
        -> candidate ranking

Each stage is a function over explicit inputs, so any one of them can be
replaced without touching the rest. The pipeline never discards a variant for
being intracellular: a nuclear driver loses the direct-surface route and keeps
the peptide/HLA and pathway routes, and the report says which is which.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from genomeos.cancer.tumour import CODING, TumourVariant

from . import logic as combination_logic
from . import mechanisms as mech
from .evidence import derived, strongest_level
from .expression import (
    normal_profile,
    normal_tissue_safety,
    tumour_selectivity,
    tumour_state,
)
from .localisation import localise, mutation_topology
from .model import (
    Localisation,
    ScoreComponent,
    TherapeuticTargetCandidate,
    VariantOrigin,
)
from .neoantigen import assess as assess_neoantigen
from .neoantigen import neoantigen_strength, presentation_confidence
from .providers import Providers
from .scoring import (
    assemble,
    clinical_precedent,
    clonality,
    component,
    shedding_score,
    summary_evidence,
    surface_accessibility,
    tumour_expression_score,
)
from .structure import context as structure_context
from .structure import mutation_epitope, structural_bindability
from .trafficking import assess as assess_trafficking
from .trafficking import internalisation_score, lysosomal_score

#: What a patient profile can contain, poorest first. The report states which
#: level was actually available so a conclusion is never read above its data.
DATA_LEVELS: tuple[tuple[int, str, str], ...] = (
    (1, "tumour_vcf", "tumour DNA variants"),
    (2, "matched_normal", "matched normal genome, so somatic calls are certain"),
    (3, "tumour_rna", "tumour RNA-seq, so expression is measured rather than assumed"),
    (4, "copy_number", "copy number and structural variation"),
    (5, "proteomics", "tumour proteomics"),
    (6, "surface_proteomics", "surface proteomics, the only direct measure of a surface target"),
    (7, "hla", "HLA genotype, without which no peptide can be prioritised"),
    (8, "immunopeptidomics", "immunopeptidomics, the only proof a peptide is presented"),
    (9, "single_cell", "single-cell tumour data, for heterogeneity and co-expression"),
)

CLONAL_VAF = 0.35
AF_INFO = re.compile(r"(?:^|;)AF=([0-9.]+)")

DISCLAIMER = (
    "GenomeOS generates computational research hypotheses. Therapeutic target rankings and mechanism "
    "suggestions are not clinical recommendations and require experimental and clinical validation."
)


@dataclass(slots=True)
class PatientProfile:
    """Everything known about this patient, and everything that is not."""

    sample_id: str
    tumour_vcf: str
    normal_vcf: str | None = None
    rna_path: str | None = None
    hla_alleles: list[str] = field(default_factory=list)
    purity: float | None = None
    copy_number: dict[str, float] = field(default_factory=dict)

    def levels(self) -> dict[str, Any]:
        have = {
            "tumour_vcf": bool(self.tumour_vcf),
            "matched_normal": bool(self.normal_vcf),
            "tumour_rna": bool(self.rna_path),
            "copy_number": bool(self.copy_number),
            "proteomics": False,
            "surface_proteomics": False,
            "hla": bool(self.hla_alleles),
            "immunopeptidomics": False,
            "single_cell": False,
        }
        reached = 0
        for n, key, _desc in DATA_LEVELS:
            if have[key]:
                reached = max(reached, n)
        contiguous = 0
        for n, key, _desc in DATA_LEVELS:
            if not have[key]:
                break
            contiguous = n
        return {
            "level_reached": contiguous,
            "highest_input_present": reached,
            "available": [
                {"level": n, "input": key, "purpose": desc, "present": have[key]}
                for n, key, desc in DATA_LEVELS
            ],
            "note": (
                f"analysis ran at data level {contiguous}; conclusions that require a higher level are "
                "reported as unavailable rather than estimated"
            ),
        }


def vaf_table(vcf: str) -> dict[str, float]:
    """Variant allele fractions from a VCF, from INFO AF or FORMAT AD/AF."""
    out: dict[str, float] = {}
    path = Path(vcf)
    if not path.exists():
        return out
    opener = __import__("gzip").open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            key = f"{f[0]}:{f[1]}:{f[3]}:{f[4].split(',')[0]}"
            value: float | None = None
            m = AF_INFO.search(f[7])
            if m:
                try:
                    value = float(m.group(1))
                except ValueError:
                    value = None
            if value is None and len(f) >= 10:
                keys, vals = f[8].split(":"), f[9].split(":")
                sample = dict(zip(keys, vals, strict=False))
                if "AF" in sample:
                    try:
                        value = float(sample["AF"].split(",")[0])
                    except ValueError:
                        value = None
                elif "AD" in sample:
                    try:
                        counts = [int(x) for x in sample["AD"].split(",")]
                        total = sum(counts)
                        value = counts[1] / total if total and len(counts) > 1 else None
                    except ValueError:
                        value = None
            if value is not None:
                out[key] = value
    return out


def origin_of(
    v: TumourVariant,
    sample_id: str,
    vafs: dict[str, float],
    purity: float | None,
    copy_number: dict[str, float],
) -> VariantOrigin:
    """The tumour DNA record a candidate must always be traceable back to."""
    key = f"{v.chrom}:{v.pos}:{v.ref}:{v.alt}"
    vaf = vafs.get(key)
    cn = copy_number.get(v.gene)
    o = VariantOrigin(
        sample_id=sample_id,
        chromosome=v.chrom,
        position=v.pos,
        reference=v.ref,
        alternate=v.alt,
        gene=v.gene,
        protein_change=v.protein_change or None,
        hgvsc=v.hgvsc or None,
        hgvsp=v.hgvsp or None,
        residue=v.residue,
        variant_type=v.consequence,
        somatic_status="somatic_candidate" if not v.likely_germline else "likely_germline",
        somatic_basis=(
            "no matched normal supplied; population frequency used to set inherited variants aside"
        ),
        vaf=vaf,
        copy_number=cn,
        cosmic=list(v.cosmic),
        gnomad_af=v.gnomad_af,
        driver_frequency=v.driver_frequency,
        hotspot=v.hotspot,
    )
    if vaf is not None:
        adjusted = vaf / purity if purity else vaf
        o.clonality = "clonal" if adjusted >= CLONAL_VAF else "subclonal"
    return o


def _go_terms(providers: Providers, gene: str) -> frozenset[str] | None:
    kb = providers.trafficking.knowledge_base
    if kb is None:
        return None
    try:
        return kb.annotations.terms(gene)
    except Exception:  # noqa: BLE001 - GO files are optional
        return None


def build_candidate(
    gene: str,
    variants: list[TumourVariant],
    profile: PatientProfile,
    providers: Providers,
    vafs: dict[str, float],
    origin_class: str = "altered_gene",
) -> TherapeuticTargetCandidate:
    """Run every stage for one gene and return the assembled candidate."""
    c = TherapeuticTargetCandidate(gene=gene)
    c.origins = [origin_of(v, profile.sample_id, vafs, profile.purity, profile.copy_number) for v in variants]

    # --- stage 1: protein annotation
    ann = providers.protein.annotation(gene)
    if not ann.available:
        c.target_class = "unknown"
        c.class_reason = ann.reason
        c.ledger.lack("protein annotation", ann.reason, ("a reviewed UniProt entry for this gene",))
        c.limitations.append(f"no protein record: {ann.reason}")
        c.scores = assemble([], 0.0)
        return c
    defn = ann.data
    c.ledger.add(*ann.evidence)
    ident = (defn.get("sections", {}).get("identity") or {}).get("items") or {}
    c.protein = ident.get("name") or ""
    c.uniprot = ident.get("accession")
    sequence = ident.get("sequence")
    ensembl_gene = ((defn.get("sections", {}).get("genomic_origin") or {}).get("items") or {}).get("gene_id")

    # --- stage 2: localisation and tumour accessibility
    c.localization = localise(gene, defn, _go_terms(providers, gene))
    c.ledger.add(*c.localization.evidence)

    # --- stage 3: expression and selectivity
    hpa = providers.expression.normal_expression(gene)
    c.normal_tissue = normal_profile(gene, hpa)
    c.ledger.add(*c.normal_tissue.evidence)
    if not hpa.available:
        c.ledger.lack(
            "healthy-tissue expression",
            hpa.reason,
            ("Human Protein Atlas consensus RNA, or a local GTEx extract",),
        )
    rna = providers.patient_rna.tumour_expression(gene)
    c.tumour = tumour_state(
        gene,
        [f"{o.variant_type} {o.protein_change or ''}".strip() for o in c.origins],
        rna,
        c.normal_tissue,
        hpa.data if hpa.available else None,
    )
    c.ledger.add(*c.tumour.evidence)
    if not rna.available:
        c.ledger.lack(
            "tumour expression of this gene",
            rna.reason,
            ("tumour RNA-seq", "tumour proteomics", "surface proteomics"),
        )
    clonal = next((o for o in c.origins if o.clonality != "unknown"), None)
    if clonal:
        c.tumour.clonality = clonal.clonality
        c.tumour.vaf = clonal.vaf
        c.tumour.purity = profile.purity
    else:
        c.ledger.lack(
            "clonality",
            c.tumour.clonality_reason,
            ("variant allele fraction in the VCF", "tumour purity", "copy number"),
        )

    # --- stage 4: trafficking
    c.trafficking = assess_trafficking(gene, providers.trafficking.trafficking(gene, defn), c.localization)
    c.ledger.add(*c.trafficking.evidence)
    if c.trafficking.internalises is None:
        c.ledger.lack(
            "internalisation",
            c.trafficking.reason or "no internalisation evidence found",
            ("receptor internalisation assay", "curated endocytosis annotation"),
        )

    # --- stage 5: structure and the mutation-created epitope
    structures = providers.structure.structures(gene, defn)
    lead = _lead_variant(variants)
    residue = lead.residue if lead else None
    c.structure = structure_context(
        gene, defn, c.localization, residue, structures.data if structures.available else None
    )
    c.ledger.add(*c.structure.evidence)
    if not structures.available:
        c.ledger.lack("structure", structures.reason, ("a PDB entry", "an AlphaFold model"))
    for v in variants:
        if not v.protein_change or v.residue is None:
            continue
        _outside, tev = mutation_topology(c.localization, v.residue, gene, v.protein_change)
        ep = mutation_epitope(gene, v.protein_change, v.residue, sequence, c.localization, c.structure, tev)
        if ep:
            c.structure.candidate_epitopes.append(ep)
            c.ledger.add(*ep.evidence)

    # --- stage 6: the peptide/HLA route, run for every variant class
    if lead is not None:
        c.neoantigen = assess_neoantigen(
            gene,
            lead.consequence,
            lead.protein_change,
            lead.residue,
            sequence,
            profile.hla_alleles,
            providers.neoantigen,
        )
        c.ledger.add(*c.neoantigen.evidence)
        c.ledger.missing.extend(c.neoantigen.missing)

    # --- stage 7: therapeutic precedent
    prec = providers.precedent.precedent(gene, ensembl_gene)
    c.precedent = prec.data if prec.available else None
    c.ledger.add(*prec.evidence)
    if not prec.available:
        c.ledger.lack("therapeutic precedent", prec.reason, ("Open Targets Platform access",))

    # --- stage 8: classification
    c.target_class, c.class_reason = classify(c, origin_class)

    # --- stage 9: scores
    c.scores = score_candidate(c, prec.available)
    c.ledger.add(summary_evidence(c.scores, gene))

    # --- stage 10: mechanisms
    c.therapeutic_mechanisms = select_mechanisms(c)
    c.target_logic = combination_logic.single(c)
    c.limitations = limitations_for(c)
    c.why_interesting = why_interesting(c)
    return c


def _lead_variant(variants: list[TumourVariant]) -> TumourVariant | None:
    """The alteration that decides this gene's story: the highest-scoring coding one."""
    coding = [v for v in variants if v.consequence in CODING]
    pool = coding or variants
    return max(pool, key=lambda v: v.score) if pool else None


def classify(c: TherapeuticTargetCandidate, origin_class: str) -> tuple[str, str]:
    """Assign the target class from localisation, never from the mutation alone."""
    loc: Localisation = c.localization
    if origin_class == "pathway_induced":
        return "pathway_induced_surface", (
            "not altered in this tumour; reached as a surface protein associated with a disrupted "
            "driver, so it is a hypothesis about an induced phenotype, not an observed alteration"
        )
    if loc.primary == "unknown" or not loc.compartments or max(loc.compartments.values()) == 0.0:
        return "unknown", "no curated localisation evidence, so accessibility cannot be established"
    if loc.reachable and loc.plasma_membrane is not False:
        outside = [r for r in loc.extracellular_regions]
        if outside:
            longest = max(outside, key=lambda r: (r.end or 0) - (r.start or 0))
            span = (longest.end or 0) - (longest.start or 0) + 1
            return "direct_surface", (
                f"curated topology places {span} residues ({longest.start}-{longest.end}) outside the "
                "cell, within reach of a circulating binder"
            )
        return "direct_surface", (
            f"curated localisation places {c.gene} at the {loc.primary.replace('_', ' ')}; no "
            "extracellular topological domain is curated, so the exposed region is not delimited"
        )
    if loc.compartments.get("secreted", 0.0) >= 0.6:
        return "secreted", (
            f"{c.gene} is secreted: it is reachable in the circulation but not attached to the tumour "
            "cell, so cell-directed mechanisms do not apply"
        )
    if c.neoantigen is not None and c.neoantigen.peptides:
        return "neoantigen_hla", (
            f"{c.gene} is {loc.primary.replace('_', ' ')} and out of reach of a surface binder, but the "
            "alteration yields mutation-derived peptides that HLA could present"
        )
    return "intracellular_only", (
        f"{c.gene} is {loc.primary.replace('_', ' ')}; no circulating binder reaches it and this "
        "alteration yields no mutation-derived peptide"
    )


def score_candidate(c: TherapeuticTargetCandidate, precedent_available: bool) -> Any:
    """Every dimension, each with the sentence that explains it."""
    sa, sa_basis = surface_accessibility(c.localization)
    sel, sel_basis, sel_ev = tumour_selectivity(c.gene, c.tumour, c.normal_tissue)
    exp, exp_basis = tumour_expression_score(c.tumour)
    nts, nts_basis = normal_tissue_safety(c.normal_tissue)
    clo, clo_basis = clonality(c.tumour)
    inter, inter_basis = internalisation_score(c.trafficking)
    lyso, lyso_basis = lysosomal_score(c.trafficking)
    struct, struct_basis = structural_bindability(c.structure, c.localization)
    neo, neo_basis = neoantigen_strength(c.neoantigen)
    pres, pres_basis = presentation_confidence(c.neoantigen)
    prec, prec_basis = clinical_precedent(c.precedent, precedent_available)
    shed, shed_basis = shedding_score(c.trafficking)
    components: list[ScoreComponent] = [
        component("surface_accessibility", sa, sa_basis),
        component("tumour_selectivity", sel, sel_basis, sel_ev),
        component("tumour_expression", exp, exp_basis),
        component("normal_tissue_safety", nts, nts_basis),
        component("clonality", clo, clo_basis),
        component("internalisation", inter, inter_basis),
        component("lysosomal_trafficking", lyso, lyso_basis),
        component("structural_bindability", struct, struct_basis),
        component("neoantigen_strength", neo, neo_basis),
        component("presentation_confidence", pres, pres_basis),
        component("clinical_precedent", prec, prec_basis),
        component(
            "evidence_strength",
            round(c.ledger.strength(), 3),
            f"strongest evidence level available is {strongest_level(c.ledger.items)}, over "
            f"{len(c.ledger.items)} evidence records",
        ),
        component("shedding", shed, shed_basis),
    ]
    return assemble(components, c.ledger.strength())


def select_mechanisms(c: TherapeuticTargetCandidate) -> list[Any]:
    """Score every mechanism in the ontology against this candidate."""
    inputs = mech.Inputs(
        values={k: v.value for k, v in c.scores.components.items()},
        bases={k: (v.unknown_reason or v.basis) for k, v in c.scores.components.items()},
    )
    out = []
    for spec in mech.MECHANISMS.values():
        fit = mech.evaluate(c, spec, inputs)
        payload, ev = mech.precedent_for(spec, c.precedent)
        fit.precedent = payload
        fit.evidence.extend(ev)
        if payload and payload["approved"] and fit.viable:
            fit.notes.append(
                f"{payload['approved']} approved therapies of this modality already engage this target"
            )
        out.append(fit)
    out.sort(key=lambda f: (-f.compatibility, f.mechanism))
    return out


def limitations_for(c: TherapeuticTargetCandidate) -> list[str]:
    out: list[str] = []
    if not c.tumour.expression.known:
        out.append(
            "surface expression in this tumour cannot be established: DNA evidence present, RNA "
            "evidence unavailable, protein evidence unavailable"
        )
    if c.tumour.clonality == "unknown":
        out.append(
            "clonality is unknown, so it is not known whether this alteration is present in every "
            "tumour cell or in a minority subclone"
        )
    if not c.normal_tissue.known:
        out.append("healthy-tissue expression is unavailable, so on-target/off-tumour risk is unassessed")
    if c.trafficking.internalises is None and c.localization.reachable:
        out.append("internalisation is unknown, so payload-delivery mechanisms cannot be assessed")
    if c.neoantigen and c.neoantigen.applicable and not c.neoantigen.hla_alleles:
        out.append("no HLA genotype supplied, so no mutation-derived peptide can be prioritised")
    if any(m.status == "experimental" and m.compatibility > 0 for m in c.therapeutic_mechanisms):
        out.append(
            "experimental mechanisms are listed for completeness; their compatibility scores describe "
            "biological fit, not evidence that the approach works"
        )
    if any(o.variant_type in ("stop_gained", "frameshift_variant") for o in c.origins):
        out.append(
            "a premature stop or frameshift often triggers nonsense-mediated decay, which removes the "
            "transcript rather than producing a truncated protein; whether any product is made here is "
            "unmeasured and would need tumour RNA-seq"
        )
        out.append(mech.REPAIR_CAVEAT)
    return out


def why_interesting(c: TherapeuticTargetCandidate) -> str:
    bits = []
    drivers = [o for o in c.origins if o.driver_frequency]
    if drivers:
        d = max(drivers, key=lambda o: o.driver_frequency or 0)
        bits.append(f"{c.gene} is a recurrent driver ({(d.driver_frequency or 0):.1%} of tumours)")
    if c.target_class == "direct_surface":
        bits.append("its product is physically reachable from outside the cell")
    elif c.target_class == "neoantigen_hla":
        bits.append("its alteration yields peptides that HLA could present, though the protein is internal")
    elif c.target_class == "pathway_induced_surface":
        bits.append("it is a surface protein associated with a disrupted driver")
    if c.precedent:
        rows = ((c.precedent.get("drugAndClinicalCandidates") or {}).get("rows")) or []
        approved = [r for r in rows if r.get("maxClinicalStage") == "APPROVAL"]
        if approved:
            bits.append(f"{len(approved)} approved therapies already engage it")
    if c.normal_tissue.on_target_off_tumour_risk in ("low",):
        bits.append("healthy-tissue expression is low across the queried organs")
    return "; ".join(bits) or f"{c.gene} is altered in this tumour"


# --- pathway-induced candidates -------------------------------------------------------


def pathway_induced(
    disrupted: list[TherapeuticTargetCandidate],
    providers: Providers,
    profile: PatientProfile,
    limit: int = 4,
) -> list[TherapeuticTargetCandidate]:
    """Surface proteins associated with a disrupted driver.

    A nuclear driver that cannot be bound directly may still change what the
    cell displays. GenomeOS cannot observe that change without tumour surface
    data, so it does the part it can: name the surface-localised proteins most
    strongly associated with the disrupted driver, and mark every one of them
    as a hypothesis that needs expression evidence before it means anything.
    """
    out: list[TherapeuticTargetCandidate] = []
    seen = {c.gene for c in disrupted}
    for driver in disrupted:
        ann = providers.protein.annotation(driver.gene)
        if not ann.available:
            continue
        partners = ((ann.data.get("sections", {}).get("interactions") or {}).get("items")) or []
        ranked = sorted((p for p in partners if p.get("physical_evidence")), key=lambda p: -p.get("score", 0))
        for p in ranked:
            gene = p["partner"]
            if gene in seen or len(out) >= limit:
                continue
            sub = providers.protein.annotation(gene)
            if not sub.available:
                continue
            loc = localise(gene, sub.data, _go_terms(providers, gene))
            if not loc.reachable or loc.plasma_membrane is False:
                continue
            seen.add(gene)
            c = build_candidate(gene, [], profile, providers, {}, origin_class="pathway_induced")
            c.ledger.add(
                derived(
                    "GenomeOS pathway-induced reasoning",
                    f"{gene} is a surface protein with physical-evidence association to the disrupted "
                    f"driver {driver.gene} (STRING combined score {p.get('score'):.2f})",
                    0.3,
                )
            )
            c.ledger.lack(
                "induced surface phenotype",
                f"association with {driver.gene} does not show that this tumour changes {gene} on its "
                "surface; that needs measurement",
                ("tumour RNA-seq", "tumour surface proteomics", "single-cell tumour data"),
            )
            c.limitations.insert(
                0,
                f"indirect candidate: {gene} is not altered in this tumour. It is proposed because it is "
                f"a surface protein associated with the disrupted driver {driver.gene}. Whether the "
                "tumour actually displays more of it is unmeasured.",
            )
            c.why_interesting = (
                f"surface protein associated with disrupted {driver.gene}; a candidate indirect route "
                "when the driver itself cannot be reached"
            )
            out.append(c)
        if len(out) >= limit:
            break
    return out


# --- the analysis ---------------------------------------------------------------------


def input_hash(path: str) -> str | None:
    p = Path(path)
    if not p.exists():
        return None
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return f"sha256:{h.hexdigest()[:32]}"


def missing_data_report(candidates: list[TherapeuticTargetCandidate], profile: PatientProfile) -> list[dict]:
    """What would most improve confidence, in the order it would help."""
    wanted = [
        (
            "matched normal genome",
            not profile.normal_vcf,
            "somatic status is estimated from population frequency; rare inherited variants cannot be told "
            "apart without it",
        ),
        (
            "tumour RNA-seq",
            not profile.rna_path,
            "no target's expression in this tumour can be established from DNA alone",
        ),
        (
            "HLA genotype",
            not profile.hla_alleles,
            "presentation is allele-specific; without it no mutant peptide can be prioritised",
        ),
        (
            "tumour purity and variant allele fraction",
            not profile.purity,
            "clonal and subclonal alterations cannot be separated, so an antigen on 10% of cells looks like "
            "one on 100%",
        ),
        (
            "copy number",
            not profile.copy_number,
            "amplification is one of the strongest reasons a surface target is over-displayed",
        ),
        ("tumour proteomics", True, "RNA is not protein; abundance at the protein level is unmeasured"),
        (
            "surface proteomics",
            True,
            "the only direct measurement of what is actually displayed on the tumour cell surface",
        ),
        (
            "immunopeptidomics",
            True,
            "the only evidence that a mutant peptide is genuinely presented rather than predicted",
        ),
        (
            "single-cell tumour data",
            True,
            "heterogeneity and co-expression of two antigens in the same cell are otherwise unknowable",
        ),
    ]
    gaps: dict[str, int] = {}
    for c in candidates:
        for m in c.ledger.missing:
            gaps[m.what] = gaps.get(m.what, 0) + 1
    out = [
        {"input": name, "missing": True, "why": why, "candidates_affected": len(candidates)}
        for name, absent, why in wanted
        if absent
    ]
    for what, n in sorted(gaps.items(), key=lambda kv: -kv[1]):
        if not any(o["input"] == what for o in out):
            out.append(
                {
                    "input": what,
                    "missing": True,
                    "why": "unresolved for this analysis",
                    "candidates_affected": n,
                }
            )
    return out


def analyse(
    ranked: list[TumourVariant],
    profile: PatientProfile,
    providers: Providers | None = None,
    top_genes: int = 12,
    indirect: bool = True,
    log: Any = None,
) -> dict[str, Any]:
    """Run the pipeline over one tumour's ranked variants."""
    providers = providers or Providers.default()
    vafs = vaf_table(profile.tumour_vcf)
    somatic = [v for v in ranked if not v.likely_germline and v.gene]
    by_gene: dict[str, list[TumourVariant]] = {}
    for v in somatic:
        by_gene.setdefault(v.gene, []).append(v)
    order = sorted(by_gene, key=lambda g: -max(v.score for v in by_gene[g]))[:top_genes]

    candidates: list[TherapeuticTargetCandidate] = []
    for i, gene in enumerate(order, 1):
        if log:
            print(f"therapeutics: {i}/{len(order)} {gene}", file=log, flush=True)
        candidates.append(build_candidate(gene, by_gene[gene], profile, providers, vafs))

    if indirect:
        blocked = [
            c
            for c in candidates
            if c.target_class in ("intracellular_only", "neoantigen_hla")
            and any(o.driver_frequency for o in c.origins)
        ]
        if blocked:
            if log:
                print(f"therapeutics: indirect routes for {[c.gene for c in blocked]}", file=log, flush=True)
            candidates.extend(pathway_induced(blocked, providers, profile))

    for c in candidates:
        c.scores = score_candidate(c, c.precedent is not None)
        c.therapeutic_mechanisms = select_mechanisms(c)

    combos = combination_logic.pairs(candidates)
    candidates.sort(key=lambda c: (-(c.scores.overall or 0.0), c.gene))
    return {
        # the live provider bundle, so the dataset stage can still ask questions;
        # it is never serialised (machine_report and dataset take what they need)
        "provider_bundle": providers,
        "sample": profile.tumour_vcf,
        "sample_id": profile.sample_id,
        "generated": derived("GenomeOS", "therapeutic target analysis", 0.0).retrieved_at,
        "disclaimer": DISCLAIMER,
        "data_level": profile.levels(),
        "candidates": candidates,
        "combinations": combos,
        "missing_data": missing_data_report(candidates, profile),
        "providers": providers.versions(),
        "input_hashes": {
            k: v
            for k, v in (
                ("tumour_vcf", input_hash(profile.tumour_vcf)),
                ("normal_vcf", input_hash(profile.normal_vcf) if profile.normal_vcf else None),
                ("tumour_rna", input_hash(profile.rna_path) if profile.rna_path else None),
            )
            if v
        },
    }


def rank_variants(
    tumour_vcf: str,
    normal_vcf: str | None = None,
    chroms: set[str] | None = None,
    knowledge: dict[str, Any] | None = None,
    log: Any = None,
) -> list[TumourVariant]:
    """Annotate and grade a tumour's variants with the existing cancer pipeline.

    With a matched normal the shared variants are subtracted first, which is
    what makes the remainder genuinely somatic rather than estimated.
    """
    from genomeos.cancer.tumour import annotate_vep, grade
    from genomeos.genome.variants import iter_vcf

    if normal_vcf:
        from genomeos.cancer.compare import somatic

        variants = somatic(normal_vcf, tumour_vcf, chroms)
    else:
        variants = list(iter_vcf(tumour_vcf, chroms, pass_only=False))
    vep = annotate_vep(variants, log=log)
    ranked = grade(variants, vep, knowledge)
    if normal_vcf:
        for v in ranked:
            v.likely_germline = False
    return ranked


def analyse_vcf(
    tumour_vcf: str,
    normal_vcf: str | None = None,
    hla: list[str] | None = None,
    rna: str | None = None,
    chroms: set[str] | None = None,
    knowledge: dict[str, Any] | None = None,
    top_genes: int = 12,
    purity: float | None = None,
    net: bool = True,
    indirect: bool = True,
    log: Any = None,
) -> dict[str, Any]:
    """End to end: a tumour VCF in, a therapeutic target analysis out."""
    from genomeos.results import load_result

    from .providers import PatientRnaProvider, Providers

    knowledge = knowledge if knowledge is not None else load_result("cancer_msk_impact_2017")
    ranked = rank_variants(tumour_vcf, normal_vcf, chroms, knowledge, log)
    profile = PatientProfile(
        sample_id=Path(tumour_vcf).stem,
        tumour_vcf=tumour_vcf,
        normal_vcf=normal_vcf,
        rna_path=rna,
        hla_alleles=list(hla or []),
        purity=purity,
    )
    kb = None
    try:
        from genomeos.lib import KnowledgeBase

        if KnowledgeBase.available():
            kb = KnowledgeBase()
    except Exception:  # noqa: BLE001 - the GO files are optional
        kb = None
    patient_rna = PatientRnaProvider.from_file(rna) if rna else PatientRnaProvider()
    providers = Providers.default(net=net, knowledge_base=kb, patient_rna=patient_rna)
    analysis = analyse(ranked, profile, providers, top_genes=top_genes, indirect=indirect, log=log)
    analysis["variants_total"] = len(ranked)
    analysis["somatic_candidates"] = sum(1 for v in ranked if not v.likely_germline)
    return analysis


def write_outputs(analysis: dict[str, Any], out_dir: str | Path, providers: Any = None) -> dict[str, str]:
    """The five dataset layers plus the report, written side by side."""
    import json

    from .design import dataset, evidence_graph, negative_set, positive_set
    from .report import machine_report, text_report

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    design = dataset(analysis, providers or analysis.get("provider_bundle"))
    files = {
        "therapeutic_analysis.json": machine_report(analysis),
        "therapeutic_design_dataset.json": design,
        "target_positive_set.json": positive_set(design),
        "target_negative_set.json": negative_set(design),
        "evidence_graph.json": evidence_graph(analysis),
    }
    written: dict[str, str] = {}
    for name, payload in files.items():
        p = out / name
        p.write_text(json.dumps(payload, indent=1, default=str))
        written[name] = str(p)
    report = out / "therapeutic_report.txt"
    report.write_text(text_report(analysis))
    written["therapeutic_report.txt"] = str(report)
    return written
