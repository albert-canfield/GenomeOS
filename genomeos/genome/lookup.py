"""Look one variant up through every layer GenomeOS has.

`lookup("chr21", 25897620, "C", "T")` returns, with the evidence of each:

  consequence   Ensembl VEP on any chromosome (canonical transcript): effect,
                gene, HGVS c./p., SIFT and PolyPhen, gnomAD frequency
  known         ClinVar significance, dbSNP id, COSMIC ids, PubMed ids (VEP's
                colocated variants)
  local trace   when the chromosome's models and sequence are local, the
                CentralDogmaTrace consequence, checked against VEP
  protein       the compiled definition's features at the residue (domains,
                regions, sites, modifications), structures available
  pathways      for a truncating change, the Reactome reactions lost when the
                protein is treated as absent (reachability)

One VEP call per variant, cached; nothing else leaves the machine except
the first compile of the protein.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

CACHE = Path("data/knowledge/vep/lookup_cache.jsonl")
TRUNCATING = {
    "stop_gained",
    "frameshift_variant",
    "start_lost",
    "splice_acceptor_variant",
    "splice_donor_variant",
}
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


def parse_variant(text: str) -> tuple[str, int, str, str]:
    """'chr21:25897620 C>T', 'chr21:25897620C>T', '21 25897620 C T' → (chrom, pos 1-based, ref, alt)."""
    t = text.strip().replace(",", "")
    m = re.match(r"^(chr)?([0-9XYM]+|MT)[:\s]+(\d+)[\s:]*([ACGT]+)\s*[>/\s]\s*([ACGT]+)$", t, re.I)
    if not m:
        raise ValueError("variant format: chr21:25897620 C>T (1-based, plus strand)")
    chrom = m.group(2).upper().replace("MT", "M")
    return f"chr{chrom}", int(m.group(3)), m.group(4).upper(), m.group(5).upper()


def normalise_record(rec: dict[str, Any]) -> dict[str, Any]:
    """The fields kept from one raw VEP record, including what is known about the variant."""
    tc = (rec.get("transcript_consequences") or [{}])[0]
    freqs: dict[str, float] = {}
    ids: list[str] = []
    clin: list[str] = []
    pubmed: list[int] = []
    for c in rec.get("colocated_variants") or []:
        if c.get("id"):
            ids.append(c["id"])
        for _allele, f in (c.get("frequencies") or {}).items():
            freqs.update({k: v for k, v in f.items() if isinstance(v, int | float)})
        clin += [x for x in (c.get("clin_sig") or []) if x not in clin]
        pubmed += [p for p in (c.get("pubmed") or []) if p not in pubmed]
    hgvsp = (tc.get("hgvsp") or "").split(":")[-1]
    m = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=|\*)?", hgvsp)
    residue = int(m.group(2)) if m else None
    aa_from = AA1.get(m.group(1), m.group(1)) if m else None
    aa_to = (
        (AA1.get(m.group(3), m.group(3)) if m and m.group(3) and m.group(3) != "=" else aa_from)
        if m
        else None
    )
    return {
        "consequence": rec.get("most_severe_consequence") or "",
        "gene": tc.get("gene_symbol") or "",
        "transcript": tc.get("transcript_id") or "",
        "hgvsc": (tc.get("hgvsc") or "").split(":")[-1],
        "hgvsp": hgvsp,
        "residue": residue,
        "aa_from": aa_from,
        "aa_to": aa_to,
        "sift": tc.get("sift_prediction") or "",
        "polyphen": tc.get("polyphen_prediction") or "",
        "gnomad_af": max(freqs.values()) if freqs else None,
        "dbsnp": [i for i in ids if i.startswith("rs")],
        "cosmic": [i for i in ids if i.startswith("COS")],
        "clinvar": clin,
        "pubmed": pubmed[:20],
        "pubmed_count": len(pubmed),
    }


def _vep(chrom: str, pos: int, ref: str, alt: str) -> dict[str, Any]:
    key = f"{chrom}:{pos}:{ref}:{alt}"
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    if CACHE.exists():
        for line in CACHE.read_text().splitlines():
            k, _, rest = line.partition("\t")
            if k == key:
                return json.loads(rest)
    from genomeos.cancer.tumour import vep_post

    line = f"{chrom.removeprefix('chr')} {pos} . {ref} {alt} . . ."
    recs = vep_post([line])
    norm = normalise_record(recs[0]) if recs else {"consequence": "unannotated"}
    with CACHE.open("a") as fh:
        fh.write(f"{key}\t{json.dumps(norm)}\n")
    return norm


def features_at(defn: dict[str, Any], residue: int) -> list[dict[str, Any]]:
    """UniProt features (domains, regions, sites, modifications) that contain a residue."""
    s = defn.get("sections", {})
    out = []
    for sec in ("domains", "modifications", "processing"):
        items = (s.get(sec) or {}).get("items")
        feats = items.get("features", []) if isinstance(items, dict) else (items or [])
        for f in feats:
            a, b = f.get("start"), f.get("end")
            if a is not None and b is not None and a <= residue <= b:
                out.append({"type": f["type"], "description": f.get("description", ""), "start": a, "end": b})
    return out


INDIVIDUALS = {
    "HG002": ("data/results/HG002_{chrom}.vcf", "GIAB v4.2.1 benchmark calls (a real person, open consent)"),
}
INDIVIDUAL_FALLBACK = "data/reference/{sample}_{chrom}.vcf"  # chromosomes split by `data fetch --individual`


def carriers(chrom: str, pos: int, ref: str, alt: str) -> list[dict[str, Any]]:
    """Which local individuals carry this exact variant, and with what genotype."""
    out = []
    scan: list[tuple[str, Path, str]] = []
    for name, (pattern, source) in INDIVIDUALS.items():
        path = Path(pattern.format(chrom=chrom))
        if not path.exists():
            path = Path(INDIVIDUAL_FALLBACK.format(sample=name, chrom=chrom))
        if path.exists():
            scan.append((name, path, f"measured: {source}"))
    from genomeos.genome.individuals import sources

    scan += [(n, p, e) for n, p, e in sources(chrom) if n not in INDIVIDUALS]  # imported genomes
    for name, path, source in scan:
        found = None
        with path.open() as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 10:
                    continue
                p_ = int(f[1])
                if p_ > pos:
                    break
                if p_ == pos and f[3] == ref and alt in f[4].split(","):
                    found = f[9].split(":")[0]
                    break
        out.append(
            {
                "individual": name,
                "carries": found is not None,
                "genotype": found,
                "evidence": source,
            }
        )
    return out


def lookup(
    chrom: str, pos: int, ref: str, alt: str, annotation=None, genome=None, pathways: int = 10
) -> dict[str, Any]:
    out: dict[str, Any] = {
        "variant": f"{chrom}:{pos:,} {ref}>{alt}",
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
    }
    v = _vep(chrom, pos, ref, alt)
    out["vep"] = {
        **v,
        "evidence": "curated: Ensembl VEP (canonical transcript); ClinVar, dbSNP, COSMIC, gnomAD via VEP",
    }
    known = []
    if v.get("clinvar"):
        known.append(f"ClinVar: {', '.join(v['clinvar'])}")
    if v.get("gnomad_af") is not None:
        known.append(f"population frequency {v['gnomad_af']:.3g} (gnomAD)")
    if v.get("cosmic"):
        known.append(f"in COSMIC ({len(v['cosmic'])} ids)")
    if v.get("pubmed_count"):
        known.append(f"{v['pubmed_count']} PubMed citations")
    out["known"] = known or ["no prior record in ClinVar, dbSNP, COSMIC or gnomAD: a novel variant"]
    out["individuals"] = carriers(chrom, pos, ref, alt)
    # local trace
    if annotation is not None and genome is not None and v.get("gene"):
        try:
            from genomeos.flow import trace_gene

            g = annotation.gene(v["gene"])
            module = annotation.to_module("lookup")
            tr = trace_gene(genome, g, module.entities[g.id].transcripts)
            if tr is not None:
                sub = tr.substitute(pos - 1, ref, alt)
                out["local_trace"] = {
                    "transcript": tr.transcript,
                    "consequence": sub.get("consequence"),
                    "hgvs_c": sub.get("hgvs_c"),
                    "hgvs_p": sub.get("hgvs_p"),
                    "region": sub.get("region"),
                    "agrees_with_vep": bool(sub.get("hgvs_p")) and sub.get("hgvs_p") == v.get("hgvsp"),
                    "evidence": "derived: GenomeOS trace of the canonical transcript on the local sequence",
                }
        except KeyError:
            pass
    # protein context
    if v.get("gene") and v.get("residue"):
        try:
            from genomeos.molecules import compile_protein

            d = compile_protein(v["gene"], sources={"uniprot"})
            ident = d["sections"].get("identity", {}).get("items") or {}
            if ident:
                seq = ident.get("sequence", "")
                r = v["residue"]
                out["protein"] = {
                    "accession": ident.get("accession"),
                    "name": ident.get("name"),
                    "length": ident.get("length"),
                    "residue": r,
                    "reference_residue_matches": bool(seq)
                    and r <= len(seq)
                    and seq[r - 1] == v.get("aa_from"),
                    "features_at_residue": features_at(d, r),
                    "structures_experimental": len(
                        (d["sections"].get("structures_experimental") or {}).get("items") or []
                    ),
                    "alphafold": bool((d["sections"].get("structures_predicted") or {}).get("items")),
                    "evidence": "curated: UniProt features; PDB/AlphaFold availability",
                }
                if v["consequence"] in TRUNCATING and pathways:
                    from genomeos.molecules.reactome import PathwayModel, fetch_pathway

                    pws = [x["id"] for x in (d["sections"].get("pathways", {}).get("items") or [])][:pathways]
                    lost, worst = 0, None
                    for pid in pws:
                        try:
                            k = PathwayModel.from_sbml(fetch_pathway(pid)).knockout(ident["accession"])
                        except Exception:  # noqa: BLE001
                            continue
                        lost += len(k["reactions_lost"])
                        if worst is None or k["fraction_lost"] > worst["fraction_lost"]:
                            worst = {"pathway": pid, "name": k["name"], "fraction_lost": k["fraction_lost"]}
                    out["pathways"] = {
                        "treated_as": "absent (truncating change)",
                        "checked": len(pws),
                        "reactions_lost": lost,
                        "most_affected": worst,
                        "evidence": "curated: Reactome; inferred reachability, 0.6",
                    }
        except Exception as e:  # noqa: BLE001
            out["protein"] = {"error": str(e)[:120]}
    return out
