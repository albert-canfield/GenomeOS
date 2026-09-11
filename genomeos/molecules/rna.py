"""The RNA layer of a gene: its transcripts and where they are expressed.

Transcripts come from the local gene models (GENCODE): isoforms with their
biotype, exon count, spliced length, coding length and the tags that say
which one is canonical or MANE. Expression comes from GTEx (public API, no
key): the median TPM of the gene in each of 54 tissues, measured on bulk
post-mortem tissue, cached once per gene. Between them they answer "which
RNAs does this gene make, and where".
"""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

from genomeos.runtime.central_dogma import coding_sequence, splice

CACHE = Path("data/knowledge/expression")
GTEX = "https://gtexportal.org/api/v2"
GTEX_EVIDENCE = "measured: GTEx v8 median TPM per tissue (bulk RNA-seq, post-mortem donors)"
NONCODING = {
    "lncRNA": "long non-coding RNA",
    "miRNA": "microRNA",
    "snRNA": "small nuclear RNA",
    "snoRNA": "small nucleolar RNA",
    "misc_RNA": "other small RNA",
    "rRNA": "ribosomal RNA",
    "rRNA_pseudogene": "rRNA pseudogene",
    "scaRNA": "small Cajal-body RNA",
    "retained_intron": "retained intron (non-coding isoform)",
    "nonsense_mediated_decay": "NMD target (non-coding isoform)",
    "processed_transcript": "processed transcript (non-coding)",
    "protein_coding_CDS_not_defined": "coding gene, isoform without a defined CDS",
}


def _get(url: str, timeout: int = 60) -> Any:
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "GenomeOS/0.1"})
    with urllib.request.urlopen(req, timeout=timeout) as r:  # noqa: S310
        return json.load(r)


def transcripts_of(symbol: str, annotation, genome=None) -> dict[str, Any]:
    """Every transcript of a gene from the local models, with lengths when the sequence is local."""
    g = annotation.gene(symbol)
    module = annotation.to_module("rna")
    rows = []
    biotypes: dict[str, int] = {}
    for t in module.entities[g.id].transcripts:
        bt = t.attrs.get("transcript_type", "")
        biotypes[bt] = biotypes.get(bt, 0) + 1
        row: dict[str, Any] = {
            "id": t.id,
            "name": t.attrs.get("name", t.id),
            "biotype": bt,
            "biotype_label": NONCODING.get(bt, "protein-coding isoform" if bt == "protein_coding" else bt),
            "exons": len(t.exons),
            "coding": bool(t.cds_segments),
            "canonical": "Ensembl_canonical" in t.tags,
            "mane": any(x.startswith("MANE") for x in t.tags),
            "tags": [
                x
                for x in t.tags
                if x in ("Ensembl_canonical", "MANE_Select", "basic", "CCDS", "GENCODE_Primary")
            ],
            "start": t.exons[0].start if t.exons else None,
            "end": t.exons[-1].end if t.exons else None,
        }
        if genome is not None and t.exons:
            mrna = splice(genome, t)
            cds = coding_sequence(genome, t) if t.cds_segments else ""
            row["spliced_nt"] = len(mrna)
            row["cds_nt"] = len(cds)
            row["protein_aa"] = (
                len(cds) // 3 - 1 if cds and len(cds) % 3 == 0 else (len(cds) // 3 if cds else 0)
            )
        rows.append(row)
    rows.sort(key=lambda r: (not r["canonical"], not r["mane"], not r["coding"], r["name"]))
    return {
        "gene": g.symbol,
        "gene_type": g.type,
        "transcripts": rows,
        "count": len(rows),
        "coding_isoforms": sum(1 for r in rows if r["coding"]),
        "by_biotype": dict(sorted(biotypes.items(), key=lambda kv: -kv[1])),
        "evidence": "curated: GENCODE gene models; lengths derived by splicing the reference sequence",
    }


def gtex_expression(symbol: str, cache_dir: Path = CACHE) -> dict[str, Any]:
    """Median TPM per tissue for a gene, from GTEx; cached."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    p = cache_dir / f"gtex_{symbol.upper()}.json"
    if p.exists():
        return json.loads(p.read_text())
    ref = _get(
        f"{GTEX}/reference/gene?geneId={urllib.parse.quote(symbol)}&gencodeVersion=v26&genomeBuild=GRCh38%2Fhg38"
    )
    hits = [g for g in ref.get("data", []) if g.get("geneSymbol", "").upper() == symbol.upper()]
    if not hits:
        out = {
            "gene": symbol.upper(),
            "gencode_id": None,
            "tissues": {},
            "evidence": "none",
            "error": "not in GTEx",
        }
    else:
        gid = hits[0]["gencodeId"]
        d = _get(
            f"{GTEX}/expression/medianGeneExpression?gencodeId={urllib.parse.quote(gid)}&datasetId=gtex_v8"
        )
        tissues = {x["tissueSiteDetailId"]: round(x["median"], 2) for x in d.get("data", [])}
        out = summarise_tissues(symbol.upper(), gid, tissues)
    p.write_text(json.dumps(out))
    return out


def summarise_tissues(symbol: str, gencode_id: str, tissues: dict[str, float]) -> dict[str, Any]:
    vals = sorted(tissues.values())
    n = len(vals)
    median = vals[n // 2] if n else 0.0
    top = sorted(tissues.items(), key=lambda kv: -kv[1])[:8]
    expressed = [t for t, v in tissues.items() if v >= 1.0]
    if not n or not vals[-1]:
        pattern = "not detected"
    elif len(expressed) <= 0.2 * n:
        pattern = "tissue restricted"
    elif median and vals[-1] / median >= 5:
        pattern = "tissue enhanced"
    elif len(expressed) >= 0.9 * n:
        pattern = "expressed in all tissues"
    else:
        pattern = "mixed"
    return {
        "gene": symbol,
        "gencode_id": gencode_id,
        "tissues": dict(sorted(tissues.items(), key=lambda kv: -kv[1])),
        "tissues_measured": n,
        "tissues_expressed": len(expressed),
        "median_tpm": round(median, 2),
        "max_tpm": round(vals[-1], 2) if n else 0.0,
        "top": top,
        "pattern": pattern,
        "evidence": GTEX_EVIDENCE,
        "confidence": 0.9,
    }


def rna_report(symbol: str, annotation=None, genome=None, expression: bool = True) -> dict[str, Any]:
    out: dict[str, Any] = {"gene": symbol.upper()}
    if annotation is not None:
        try:
            out["transcripts"] = transcripts_of(symbol, annotation, genome)
        except KeyError:
            out["transcripts"] = None
    if expression:
        try:
            out["expression"] = gtex_expression(symbol)
        except Exception as e:  # noqa: BLE001
            out["expression"] = {
                "gene": symbol.upper(),
                "tissues": {},
                "evidence": "none",
                "error": str(e)[:120],
            }
    return out
