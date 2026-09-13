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


GTEX_ISOFORM_EVIDENCE = "measured: GTEx v8 median transcript TPM per tissue (RSEM isoform quantification)"


def gtex_isoforms(
    symbol: str, cache_dir: Path = CACHE, transcripts: list[dict[str, Any]] | None = None
) -> dict:
    """Median TPM per transcript per tissue for a gene, from GTEx; cached. `transcripts` (name, canonical,
    protein_length per Ensembl id, as the compiled definition lists them) label the isoforms."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    p = cache_dir / f"gtex_isoforms_{symbol.upper()}.json"
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
            "isoforms": {},
            "evidence": "none",
            "error": "not in GTEx",
        }
    else:
        gid = hits[0]["gencodeId"]
        rows: list[dict[str, Any]] = []
        page = 0
        while True:
            d = _get(
                f"{GTEX}/expression/medianTranscriptExpression?gencodeId={urllib.parse.quote(gid)}"
                f"&datasetId=gtex_v8&itemsPerPage=1000&page={page}"
            )
            rows.extend(d.get("data", []))
            info = d.get("paging_info", {})
            page += 1
            if page >= int(info.get("numberOfPages", 1)):
                break
        out = summarise_isoforms(symbol.upper(), gid, rows, transcripts or [])
    p.write_text(json.dumps(out))
    return out


def summarise_isoforms(symbol: str, gencode_id: str, rows: list[dict[str, Any]], transcripts: list) -> dict:
    """Per transcript: TPM per tissue and share of the gene; per tissue: the dominant isoform; overall: how
    often the canonical transcript is the one the tissue actually makes."""
    label = {t["transcript"]: t for t in transcripts}
    per_tx: dict[str, dict[str, float]] = {}
    for x in rows:
        tid = x["transcriptId"].split(".")[0]
        per_tx.setdefault(tid, {})[x["tissueSiteDetailId"]] = round(float(x["median"]), 3)
    tissues = sorted({t for v in per_tx.values() for t in v})
    dominant: dict[str, tuple[str, float, float]] = {}
    for tissue in tissues:
        total = sum(v.get(tissue, 0.0) for v in per_tx.values())
        if total <= 0:
            continue
        best = max(per_tx, key=lambda tid: per_tx[tid].get(tissue, 0.0))
        dominant[tissue] = (
            best,
            round(per_tx[best].get(tissue, 0.0), 2),
            round(per_tx[best].get(tissue, 0.0) / total, 3),
        )
    canonical = next((tid for tid, t in label.items() if t.get("canonical")), None)
    iso = {}
    for tid, v in per_tx.items():
        vals = list(v.values())
        iso[tid] = {
            "name": label.get(tid, {}).get("name"),
            "canonical": label.get(tid, {}).get("canonical", False),
            "protein_length": label.get(tid, {}).get("protein_length"),
            "biotype": label.get(tid, {}).get("biotype"),
            "max_tpm": round(max(vals), 2) if vals else 0.0,
            "median_tpm": round(sorted(vals)[len(vals) // 2], 2) if vals else 0.0,
            "tissues_dominant": sum(1 for t, (b, _, _) in dominant.items() if b == tid),
            "tissues": v,
        }
    ranked = sorted(iso.items(), key=lambda kv: -kv[1]["tissues_dominant"])
    non_canonical = {t: d for t, d in dominant.items() if canonical and d[0] != canonical}
    return {
        "gene": symbol,
        "gencode_id": gencode_id,
        "transcripts_measured": len(per_tx),
        "tissues_measured": len(tissues),
        "canonical": canonical,
        "canonical_dominant_in": sum(1 for d in dominant.values() if canonical and d[0] == canonical),
        "tissues_where_another_isoform_dominates": {
            t: {"transcript": d[0], "name": iso[d[0]]["name"], "tpm": d[1], "share": d[2]}
            for t, d in sorted(non_canonical.items(), key=lambda kv: -kv[1][2])
        },
        "dominant_by_tissue": {
            t: {"transcript": d[0], "tpm": d[1], "share": d[2]} for t, d in dominant.items()
        },
        "isoforms": dict(ranked),
        "evidence": GTEX_ISOFORM_EVIDENCE,
        "confidence": 0.9,
    }
