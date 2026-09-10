"""Tumour-only analysis: a tumour's DNA without a matched normal.

Without the normal sample the germline cannot be subtracted, so it is
estimated: a variant common in the population (gnomAD frequency at or above
GERMLINE_AF) is treated as inherited and set aside. Every remaining variant
is annotated once through Ensembl VEP (consequence, gene, HGVS, SIFT and
PolyPhen, COSMIC ids, gnomAD frequencies) on any chromosome, graded with the
cBioPortal driver knowledge, and the top hits are traced further: the
mutant peptide around a missense change (neoantigen candidates) and the
pathway reactions lost when a truncated driver is treated as absent.

All of it is research analysis with the evidence per claim; none of it is
clinical advice.
"""

from __future__ import annotations

import json
import math
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from genomeos.genome.variants import Variant, iter_vcf
from genomeos.results import load_result

VEP_URL = "https://rest.ensembl.org/vep/homo_sapiens/region"
VEP_OPTS = "af_gnomadg=1&af_gnomade=1&canonical=1&hgvs=1&pick=1"
CACHE = Path("data/knowledge/vep")
GERMLINE_AF = 0.01
BATCH = 200
EXOME_MB = 35.0  # size of the coding exome used for the mutation-burden estimate
SEVERITY = {
    "frameshift_variant": 1.0,
    "stop_gained": 1.0,
    "start_lost": 0.9,
    "stop_lost": 0.8,
    "splice_acceptor_variant": 0.8,
    "splice_donor_variant": 0.8,
    "inframe_deletion": 0.6,
    "inframe_insertion": 0.6,
    "missense_variant": 0.5,
    "protein_altering_variant": 0.5,
    "splice_region_variant": 0.2,
    "synonymous_variant": 0.05,
}
CODING = {k for k, v in SEVERITY.items() if v >= 0.5}
AA1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Ter": "*",
}


@dataclass(slots=True)
class TumourVariant:
    chrom: str
    pos: int  # 1-based
    ref: str
    alt: str
    gene: str = ""
    consequence: str = ""
    hgvsp: str = ""
    hgvsc: str = ""
    protein_change: str = ""  # short form, R175H
    residue: int | None = None
    sift: str = ""
    polyphen: str = ""
    cosmic: list[str] = field(default_factory=list)
    gnomad_af: float | None = None
    likely_germline: bool = False
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
            "hgvsp": self.hgvsp,
            "hgvsc": self.hgvsc,
            "protein_change": self.protein_change,
            "sift": self.sift,
            "polyphen": self.polyphen,
            "cosmic": self.cosmic,
            "gnomad_af": self.gnomad_af,
            "likely_germline": self.likely_germline,
            "driver_frequency": self.driver_frequency,
            "hotspot": self.hotspot,
            "score": round(self.score, 3),
            "evidence": self.evidence,
        }


# ---- VEP adapter ---------------------------------------------------------------------


def _key(v: Variant) -> str:
    return f"{v.chrom}:{v.pos}:{v.ref}:{v.alts[0] if v.alts else ''}"


def _vep_line(v: Variant) -> str:
    # Variant.pos is 0-based inside GenomeOS; VEP and VCF are 1-based
    return f"{v.chrom.removeprefix('chr')} {v.pos + 1} . {v.ref} {v.alts[0]} . . ."


def vep_post(lines: list[str], timeout: int = 180, retries: int = 3) -> list[dict[str, Any]]:
    body = json.dumps({"variants": lines}).encode()
    delay = 2.0
    for attempt in range(retries + 1):
        req = urllib.request.Request(
            f"{VEP_URL}?{VEP_OPTS}",
            data=body,
            headers={
                "Content-Type": "application/json",
                "Accept": "application/json",
                "User-Agent": "GenomeOS/0.1",
            },
        )
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:  # noqa: S310
                return json.load(r)
        except urllib.error.HTTPError as e:
            if attempt == retries or e.code not in (429, 500, 502, 503, 504):
                raise
            time.sleep(float(e.headers.get("Retry-After") or delay))
            delay *= 2
    raise RuntimeError("unreachable")


def normalise_vep(rec: dict[str, Any]) -> dict[str, Any]:
    """The fields GenomeOS keeps from one VEP record."""
    tc = (rec.get("transcript_consequences") or [{}])[0]
    freqs: dict[str, float] = {}
    ids: list[str] = []
    for c in rec.get("colocated_variants") or []:
        if c.get("id"):
            ids.append(c["id"])
        for _allele, f in (c.get("frequencies") or {}).items():
            freqs.update({k: v for k, v in f.items() if isinstance(v, int | float)})
    hgvsp = (tc.get("hgvsp") or "").split(":")[-1]
    short, residue = "", None
    m = __import__("re").match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=|\*)?", hgvsp)
    if m:
        a, n, b = m.groups()
        residue = int(n)
        short = f"{AA1.get(a, a)}{n}{AA1.get(b, b) if b and b != '=' else ('=' if b == '=' else '')}"
    return {
        "gene": tc.get("gene_symbol") or "",
        "consequence": rec.get("most_severe_consequence") or "",
        "hgvsp": hgvsp,
        "hgvsc": (tc.get("hgvsc") or "").split(":")[-1],
        "protein_change": short,
        "residue": residue,
        "sift": tc.get("sift_prediction") or "",
        "polyphen": tc.get("polyphen_prediction") or "",
        "cosmic": [i for i in ids if i.startswith("COS")],
        "gnomad_af": max(freqs.values()) if freqs else None,
    }


def annotate_vep(variants: list[Variant], cache_dir: Path = CACHE, log=None) -> dict[str, dict[str, Any]]:
    """Annotate through VEP, batch by batch, keeping a small per-variant cache."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / "vep_cache.jsonl"
    cache: dict[str, dict[str, Any]] = {}
    if cache_file.exists():
        for line in cache_file.read_text().splitlines():
            if line.strip():
                k, _, rest = line.partition("\t")
                cache[k] = json.loads(rest)
    todo = [v for v in variants if v.alts and _key(v) not in cache]
    with cache_file.open("a") as fh:
        for i in range(0, len(todo), BATCH):
            batch = todo[i : i + BATCH]
            recs = vep_post([_vep_line(v) for v in batch])
            by_input = {r.get("input"): r for r in recs}
            for v in batch:
                rec = by_input.get(_vep_line(v))
                norm = normalise_vep(rec) if rec else {"gene": "", "consequence": "unannotated"}
                cache[_key(v)] = norm
                fh.write(f"{_key(v)}\t{json.dumps(norm)}\n")
            if log:
                print(f"VEP: {min(i + BATCH, len(todo))}/{len(todo)} annotated", file=log, flush=True)
    return {_key(v): cache[_key(v)] for v in variants if v.alts}


# ---- grading -------------------------------------------------------------------------


def grade(
    variants: list[Variant], vep: dict[str, dict[str, Any]], knowledge: dict | None = None
) -> list[TumourVariant]:
    knowledge = knowledge or load_result("cancer_msk_impact_2017") or {}
    kgenes = knowledge.get("genes", {})
    out = []
    for v in variants:
        if not v.alts:
            continue
        a = vep.get(_key(v), {})
        t = TumourVariant(
            v.chrom, v.pos + 1, v.ref, v.alts[0], **{k: a[k] for k in a if k in TumourVariant.__slots__}
        )
        if t.gnomad_af is not None and t.gnomad_af >= GERMLINE_AF:
            t.likely_germline = True
            t.evidence.append(
                f"population frequency {t.gnomad_af:.3g} (gnomAD): likely inherited, not somatic"
            )
        k = kgenes.get(t.gene)
        if k:
            t.driver_frequency = k["frequency"]
            t.evidence.append(
                f"{t.gene} mutated in {k['frequency']:.1%} of {knowledge.get('samples')} tumours "
                f"(cBioPortal {knowledge.get('study')})"
            )
            if t.protein_change and t.protein_change in {h[0] for h in k.get("hotspots", [])}:
                t.hotspot = True
                t.evidence.append(f"{t.protein_change} is a recurrent hotspot")
        if t.cosmic:
            t.evidence.append(f"seen in COSMIC ({', '.join(t.cosmic[:3])})")
        if t.sift == "deleterious" or t.polyphen.startswith("probably"):
            t.evidence.append(f"predicted damaging (SIFT {t.sift or '-'}, PolyPhen {t.polyphen or '-'})")
        t.score = _score(t)
        out.append(t)
    out.sort(key=lambda s: (s.likely_germline, -s.score))
    return out


def _score(t: TumourVariant) -> float:
    s = SEVERITY.get(t.consequence, 0.02)
    if t.driver_frequency:
        s += min(1.0, math.log10(1 + 100 * t.driver_frequency))
    if t.hotspot:
        s += 1.0
    if t.cosmic:
        s += 0.3
    if t.sift == "deleterious" or t.polyphen.startswith("probably"):
        s += 0.2
    if t.likely_germline:
        s *= 0.1
    return s


def mutation_burden(ranked: list[TumourVariant]) -> dict[str, Any]:
    coding = [t for t in ranked if not t.likely_germline and t.consequence in CODING]
    return {
        "coding_somatic_variants": len(coding),
        "per_megabase": round(len(coding) / EXOME_MB, 2),
        "assumption": f"whole-exome input over {EXOME_MB:.0f} Mb; a panel or a single chromosome gives a "
        "different denominator",
        "evidence": "inferred",
        "note": "≥10/Mb is the usual threshold for a high burden (immune-checkpoint response signal)",
    }


# ---- deeper: peptides and pathways ------------------------------------------------------


def mutant_peptides(t: TumourVariant, sequence: str, flank: int = 8) -> dict[str, Any] | None:
    """The peptide window around a missense change: wild type and mutant, for neoantigen work."""
    if not t.residue or len(t.protein_change) < 3 or t.consequence != "missense_variant":
        return None
    i = t.residue - 1
    if i >= len(sequence) or sequence[i] != t.protein_change[0]:
        have = sequence[i] if i < len(sequence) else "?"
        return {"error": f"reference residue mismatch at {t.residue}: sequence has {have}"}
    mut = sequence[:i] + t.protein_change[-1] + sequence[i + 1 :]
    a, b = max(0, i - flank), min(len(sequence), i + flank + 1)
    return {
        "gene": t.gene,
        "change": t.protein_change,
        "wild_type": sequence[a:b],
        "mutant": mut[a:b],
        "window": [a + 1, b],
        "evidence": "derived: UniProt canonical sequence with the substitution applied",
        "note": "9-mers containing the changed residue are the candidates; HLA binding is not predicted here",
    }


def analyse(
    vcf: str,
    chroms: set[str] | None = None,
    knowledge: dict | None = None,
    deep: int = 5,
    pathways: bool = True,
    log=None,
) -> dict[str, Any]:
    """The whole tumour-only pipeline; returns the analysis dict (the packet wraps it)."""
    variants = list(iter_vcf(vcf, chroms, pass_only=False))
    vep = annotate_vep(variants, log=log)
    ranked = grade(variants, vep, knowledge)
    somatic = [t for t in ranked if not t.likely_germline]
    germline = [t for t in ranked if t.likely_germline]
    peptides, effects = [], []
    if deep:
        from genomeos.molecules import compile_protein

        for t in [t for t in somatic if t.gene and t.consequence in CODING][:deep]:
            try:
                d = compile_protein(t.gene, sources={"uniprot"})
            except Exception as e:  # noqa: BLE001
                effects.append({"gene": t.gene, "error": str(e)[:100]})
                continue
            ident = d["sections"].get("identity", {}).get("items") or {}
            if ident.get("sequence"):
                pep = mutant_peptides(t, ident["sequence"])
                if pep:
                    peptides.append(pep)
            if (
                pathways
                and t.consequence in ("stop_gained", "frameshift_variant", "start_lost")
                and ident.get("accession")
            ):
                from genomeos.molecules.reactome import PathwayModel, fetch_pathway

                pws = [x["id"] for x in (d["sections"].get("pathways", {}).get("items") or [])][:10]
                lost = 0
                worst = None
                for pid in pws:
                    try:
                        k = PathwayModel.from_sbml(fetch_pathway(pid)).knockout(ident["accession"])
                    except Exception:  # noqa: BLE001
                        continue
                    lost += len(k["reactions_lost"])
                    if worst is None or k["fraction_lost"] > worst["fraction_lost"]:
                        worst = {"pathway": pid, "name": k["name"], "fraction_lost": k["fraction_lost"]}
                effects.append(
                    {
                        "gene": t.gene,
                        "change": t.protein_change,
                        "treated_as": "absent (truncating)",
                        "pathways_checked": len(pws),
                        "reactions_lost": lost,
                        "most_affected": worst,
                        "evidence": "curated: Reactome; logic inferred (reachability)",
                    }
                )
    genes = {t.gene for t in somatic if t.gene and t.consequence in CODING}
    from genomeos.cancer.compare import suggest_cancer_type, surface_targets

    return {
        "sample": vcf,
        "variants_total": len(variants),
        "annotated": sum(1 for t in ranked if t.consequence and t.consequence != "unannotated"),
        "likely_germline": len(germline),
        "somatic_candidates": len(somatic),
        "coding_somatic": sum(1 for t in somatic if t.consequence in CODING),
        "mutation_burden": mutation_burden(ranked),
        "ranked": [t.to_dict() for t in somatic[:50]],
        "germline_set_aside": [t.to_dict() for t in germline[:20]],
        "mutant_peptides": peptides,
        "pathway_effects": effects,
        "suggested_cancer_types": suggest_cancer_type(genes, knowledge),
        "surface_targets": surface_targets(genes),
        "evidence": {
            "consequence": "curated: Ensembl VEP (canonical transcript, pick one)",
            "germline_filter": f"inferred: gnomAD frequency ≥ {GERMLINE_AF} means inherited; rare germline "
            "variants cannot be told apart without the normal sample",
            "drivers": "curated: cBioPortal study frequencies and hotspots",
            "damage": "predicted: SIFT / PolyPhen",
        },
    }


def tumour_packet(analysis: dict[str, Any]) -> dict[str, Any]:
    """The AI-agent packet for a tumour-only sample."""
    return {
        "task": (
            "From the alterations of one tumour (no matched normal), name the likely driver events, the "
            "likely cancer type, the best cell-surface target for a designed binder, and neoantigen "
            "candidates, with evidence and confidence for each step."
        ),
        "sample": analysis["sample"],
        "inputs": {
            k: analysis[k]
            for k in (
                "variants_total",
                "likely_germline",
                "somatic_candidates",
                "coding_somatic",
                "mutation_burden",
                "ranked",
                "germline_set_aside",
                "mutant_peptides",
                "pathway_effects",
                "suggested_cancer_types",
                "surface_targets",
                "evidence",
            )
        },
        "constraints": [
            "Use only the evidence in inputs; mark any external knowledge as 'assumption'.",
            "Variants set aside as likely germline are not tumour-specific; do not build on them.",
            "A predicted damage score is not a measurement; a pathway loss is a reachability inference.",
            "Report confidence per claim on a 0-1 scale.",
            "This is a research analysis, not clinical advice.",
        ],
        "expected_output": {
            "driver_events": [{"gene": "string", "change": "string", "why": "string", "confidence": "0-1"}],
            "cancer_type": {"name": "string", "confidence": "0-1"},
            "target_gene": "string",
            "target_why": "string",
            "neoantigen_candidates": ["string"],
            "payload_options": ["string"],
            "risks": ["string"],
            "next_experiments": ["string"],
        },
        "minimal_prompt": (
            "You are given the annotated variants of one tumour without a matched normal: a ranked somatic "
            "candidate list, the variants set aside as likely germline, the mutation burden, mutant peptide "
            "windows, pathway reactions lost by truncated drivers, suggested cancer types and cell-surface "
            "products among the altered genes. Name the driver events and the likely cancer type, choose one "
            "surface target for a designed binder and its payload, list neoantigen candidates, risks and the "
            "next experiments. Answer in the expected_output JSON."
        ),
    }
