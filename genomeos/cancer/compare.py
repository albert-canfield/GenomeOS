"""Healthy versus tumour: somatic differences, drivers, cancer type, targets.

Inputs are two variant lists (VCFs) from the same person: a healthy sample
(blood or normal tissue) and a tumour sample. Somatic variants are those in
the tumour and not in the normal. Each is annotated with its coding
consequence (GenomeOS variant engine), matched against the distilled
cBioPortal knowledge (driver frequencies, hotspots per cancer type), and
ranked. The cancer type is suggested by comparing the sample's driver profile
with the per-type frequencies. Targets are the altered genes whose products
sit on the cell surface, which is where a designed binder could reach them.
Everything is graded: cBioPortal-derived facts are `experimental`, our
ranking is `inferred`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

from genomeos.genome.variants import Variant, iter_vcf
from genomeos.results import load_result

SURFACE_TERMS = ("GO:0005886", "GO:0009986")  # plasma membrane, cell surface


@dataclass(slots=True)
class SomaticVariant:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str = ""
    consequence: str = ""
    protein_change: str = ""
    driver_frequency: float | None = None
    hotspot: bool = False
    score: float = 0.0
    evidence: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "chrom": self.chrom,
            "pos": self.pos,
            "ref": self.ref,
            "alt": self.alt,
            "gene": self.gene,
            "consequence": self.consequence,
            "protein_change": self.protein_change,
            "driver_frequency": self.driver_frequency,
            "hotspot": self.hotspot,
            "score": round(self.score, 3),
            "evidence": self.evidence,
        }


def somatic(normal_vcf: str, tumour_vcf: str, chroms: set[str] | None = None) -> list[Variant]:
    """Variants present in the tumour and absent from the normal (same person)."""
    normal = {(v.chrom, v.pos, v.ref, v.alts) for v in iter_vcf(normal_vcf, chroms, pass_only=False)}
    return [
        v
        for v in iter_vcf(tumour_vcf, chroms, pass_only=False)
        if (v.chrom, v.pos, v.ref, v.alts) not in normal
    ]


def annotate(
    variants: list[Variant], annotation, genome, knowledge: dict | None = None
) -> list[SomaticVariant]:
    """Consequence per variant via the gene models, plus cBioPortal driver knowledge."""
    from genomeos.runtime import classify_all

    knowledge = knowledge or load_result("cancer_msk_impact_2017") or {}
    kgenes = knowledge.get("genes", {})
    module = annotation.to_module("cmp")
    out: list[SomaticVariant] = []
    for v in variants:
        alt = v.alts[0] if v.alts else ""
        sv = SomaticVariant(v.chrom, v.pos, v.ref, alt)
        hits = []
        for g in annotation.protein_coding():
            if g.locus.chrom == v.chrom and g.locus.start - 2000 <= v.pos < g.locus.end + 2000:
                txs = [t for t in module.entities[g.id].transcripts if t.cds_segments]
                if txs:
                    e = classify_all(genome, txs, Variant(v.chrom, v.pos, v.ref, v.alts, gt=(1, 1)))[0]
                    hits.append((g.symbol, e))
        if hits:
            from genomeos.runtime import severity

            sym, e = min(hits, key=lambda h: severity(h[1].consequence))
            sv.gene, sv.consequence, sv.protein_change = sym, e.consequence, e.protein_change
        k = kgenes.get(sv.gene)
        if k:
            sv.driver_frequency = k["frequency"]
            sv.evidence.append(
                f"{sv.gene} mutated in {k['frequency']:.1%} of {knowledge.get('samples')} tumours "
                f"(cBioPortal {knowledge.get('study')})"
            )
            hot = {h[0] for h in k.get("hotspots", [])}
            if sv.protein_change and sv.protein_change in hot:
                sv.hotspot = True
                sv.evidence.append(f"{sv.protein_change} is a recurrent hotspot")
        sv.score = _score(sv)
        out.append(sv)
    out.sort(key=lambda s: -s.score)
    return out


def _score(sv: SomaticVariant) -> float:
    sev = {
        "frameshift_variant": 1.0,
        "nonsense": 1.0,
        "start_lost": 0.9,
        "stop_lost": 0.8,
        "splice_site": 0.8,
        "inframe_deletion": 0.6,
        "inframe_insertion": 0.6,
        "missense_variant": 0.5,
        "synonymous_variant": 0.05,
    }
    s = sev.get(sv.consequence, 0.1)
    if sv.driver_frequency:
        s += min(1.0, math.log10(1 + 100 * sv.driver_frequency))
    if sv.hotspot:
        s += 1.0
    return s


def suggest_cancer_type(
    somatic_genes: set[str], knowledge: dict | None = None, top: int = 5
) -> list[dict[str, Any]]:
    """Rank cancer types by how well the sample's mutated driver genes match per-type frequencies."""
    knowledge = knowledge or load_result("cancer_msk_impact_2017") or {}
    genes = knowledge.get("genes", {})
    types = knowledge.get("cancer_types", {})
    scores: dict[str, float] = dict.fromkeys(types, 0.0)
    for g in somatic_genes:
        k = genes.get(g)
        if not k:
            continue
        base = k["frequency"]
        for ct, f in k.get("by_cancer_type", {}).items():
            if ct in scores:
                scores[ct] += math.log((f + 0.005) / (base + 0.005))  # enrichment of this gene in this type
    ranked = sorted(scores.items(), key=lambda kv: -kv[1])[:top]
    return [
        {
            "cancer_type": ct,
            "log_enrichment": round(s, 3),
            "samples_in_study": types[ct],
            "evidence": "inferred from per-type driver frequencies (cBioPortal); suggestion, not diagnosis",
        }
        for ct, s in ranked
    ]


def surface_targets(genes: set[str], kb=None) -> list[dict[str, Any]]:
    """Altered genes whose products are on the plasma membrane / cell surface (reachable by a binder)."""
    from genomeos.lib import KnowledgeBase

    kb = kb or KnowledgeBase()
    out = []
    for g in sorted(genes):
        terms = kb.annotations.terms(g)
        if any(t in terms for t in SURFACE_TERMS):
            out.append(
                {
                    "gene": g,
                    "location": "plasma membrane / cell surface",
                    "evidence": "GO cellular component annotation (goa_human)",
                    "confidence": 0.7,
                }
            )
    return out


def agent_packet(
    sample: str, ranked: list[SomaticVariant], types: list[dict], targets: list[dict]
) -> dict[str, Any]:
    """The dataset and instructions an AI agent needs for one focused task: propose a binder target."""
    return {
        "task": (
            "From the somatic alterations of one tumour, choose the best cell-surface target for a designed "
            "binder and state what the binder would carry, with the evidence for each step."
        ),
        "sample": sample,
        "inputs": {
            "somatic_variants_ranked": [s.to_dict() for s in ranked[:50]],
            "suggested_cancer_types": types,
            "surface_targets": targets,
        },
        "constraints": [
            "Use only the evidence in inputs; mark any external knowledge as 'assumption'.",
            "Prefer targets altered in the tumour and absent from the normal sample.",
            "Report confidence per claim on a 0-1 scale; never present a prediction as a measurement.",
            "This is a research analysis, not clinical advice.",
        ],
        "expected_output": {
            "target_gene": "string",
            "why": "string",
            "confidence": "0-1",
            "payload_options": ["string"],
            "risks": ["string"],
            "next_experiments": ["string"],
        },
        "minimal_prompt": (
            "You are given a ranked list of somatic variants of a tumour, suggested cancer types, and the "
            "altered genes whose products sit on the cell surface. Choose one target for a designed protein "
            "binder that would reach tumour cells but not healthy ones, justify it from the inputs, say what "
            "payload it should carry (toxin, immune recruiter, radionuclide, or a gene/RNA cargo) and list "
            "risks and the next experiments. Answer in the expected_output JSON."
        ),
    }
