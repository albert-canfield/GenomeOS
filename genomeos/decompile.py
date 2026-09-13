# SPDX-License-Identifier: AGPL-3.0-or-later
"""The decompiled locus: one gene, every layer GenomeOS has read for it, and what stays unknown.

The conversation behind ROADMAP area J imagined a "decompiled view" per locus: the block as
BioLang with its promoter's requirements, its elements and their targets, its node, where it
is read, how conserved it is on both axes, its orthologues and paralogues, its protein, and the
evidence under each line. Every one of those exists as a separate command and result file;
this module only assembles them. It infers nothing, fetches nothing, and names the layers it
could not find as the locus's UNKNOWN residue, so the view is as complete as the results on
disk and says so.
"""

from __future__ import annotations

import glob
import gzip
import json
from pathlib import Path
from typing import Any

from genomeos.results import RESULTS_DIR, load_result

PROTEOME = Path(__file__).resolve().parent / "lib" / "data" / "proteome.json.gz"
LAYERS = (
    "gene",
    "protein",
    "origin",
    "paralogues",
    "requires",
    "node",
    "reader",
    "elements",
    "human_axis",
    "duplication",
    "expression",
)


def _protein(symbol: str) -> dict | None:
    if not PROTEOME.exists():
        return None
    with gzip.open(PROTEOME, "rt") as fh:
        p = json.load(fh).get("proteins", {}).get(symbol)
    return p


def _gene(symbol: str, chrom: str) -> dict | None:
    from genomeos.coords import Strand
    from genomeos.genome.annotation import Annotation, default_gencode

    gff = default_gencode({chrom})
    if gff is None:
        return None
    ann = Annotation.from_gff3(gff, {chrom})
    try:
        g = ann.gene(symbol)
    except KeyError:
        return None
    ts = list(g.transcripts.values())
    canon = [t for t in ts if "Ensembl_canonical" in t.tags] or sorted(
        ts, key=lambda t: -(t.locus.end - t.locus.start)
    )
    t = canon[0] if canon else None
    return {
        "symbol": g.symbol,
        "id": g.id,
        "type": g.type,
        "start": g.locus.start,
        "end": g.locus.end,
        "strand": "+" if g.locus.strand == Strand.PLUS else "-",
        "transcripts": len(ts),
        "canonical": t.id if t else None,
        "exons": len(t.exons) if t else None,
        "cds_segments": len(t.cds) if t else None,
    }


def _node(symbol: str, gene: dict | None, chrom: str, results_dir: Path) -> dict | None:
    d = load_result(f"domains_{chrom}", results_dir) or {}
    for n in d.get("domains", []):
        if symbol in n.get("genes", []) or (gene and n["start"] <= gene["start"] < n["end"]):
            return {
                k: n.get(k)
                for k in (
                    "id",
                    "start",
                    "end",
                    "length",
                    "coding_genes",
                    "promoters",
                    "enhancers",
                    "evidence",
                )
            }
    return None


def _reader(symbol: str, chrom: str, results_dir: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for f in sorted(glob.glob(str(results_dir / f"reader_*_{chrom}.json"))):
        r = load_result(Path(f).stem, results_dir) or {}
        if "silent_genes" not in r or "_vs_" in Path(f).stem:
            continue
        out[r["cell_type"]] = "silent" if symbol in r["silent_genes"] else "read"
    return out


def _elements(symbol: str, chrom: str, results_dir: Path) -> list[dict]:
    seen: set[str] = set()
    var = {e["id"]: e for e in (load_result(f"variation_{chrom}", results_dir) or {}).get("elements", [])}
    dup_ids = set(
        ((load_result(f"duplication_{chrom}", results_dir) or {}).get("elements") or {}).get("ids", [])
    )
    mot = {e["id"]: e for e in (load_result(f"motifs_{chrom}", results_dir) or {}).get("elements", [])}
    out = []
    for name, origin in (("constrained_targets", "constrained"), ("enhancer_targets", "uniform")):
        r = load_result(f"{name}_{chrom}", results_dir) or {}
        for e in r.get("elements", []):
            pc = e.get("predicted_coding") or {}
            if pc.get("gene") != symbol or e["id"] in seen:
                continue
            seen.add(e["id"])
            v = var.get(e["id"], {})
            out.append(
                {
                    "id": e["id"],
                    "start": e["start"],
                    "end": e["end"],
                    "origin": origin,
                    "action": pc.get("action"),
                    "log2_fold_change": pc.get("log2_fold_change"),
                    "tissue": pc.get("tissue"),
                    "confidence": pc.get("confidence"),
                    "mammal_fraction": v.get("mammal_fraction", e.get("constrained_fraction")),
                    "human_fraction": (v.get("gnocchi") or {}).get("fraction_above"),
                    "case": (v.get("case") or {}).get("case"),
                    "duplicated": e["id"] in dup_ids,
                    "requires": [x["factor"] for x in (mot.get(e["id"]) or {}).get("requires", [])[:5]],
                }
            )
    out.sort(key=lambda e: e["start"])
    return out


def _expression(symbol: str) -> dict | None:
    p = Path("data/knowledge/expression") / f"gtex_{symbol}.json"
    if not p.exists():
        return None
    d = json.loads(p.read_text())
    tissues = d.get("tissues") or d.get("medianTranscriptExpression") or d
    if isinstance(tissues, dict) and all(isinstance(v, (int, float)) for v in tissues.values()):
        top = sorted(tissues.items(), key=lambda kv: -kv[1])[:6]
        return {"tissues": len(tissues), "top": [{"tissue": t, "tpm": round(v, 1)} for t, v in top]}
    return {"cached": True}


def decompile(symbol: str, chrom: str, results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    gene = _gene(symbol, chrom)
    protein = _protein(symbol)
    origin = ((load_result("origin_genome_wide", results_dir) or {}).get("genes") or {}).get(symbol)
    motifs = ((load_result(f"motifs_{chrom}", results_dir) or {}).get("genes") or {}).get(symbol)
    elements = _elements(symbol, chrom, results_dir)
    layers = {
        "gene": gene,
        "protein": (
            {
                k: protein.get(k)
                for k in (
                    "accession",
                    "length",
                    "function",
                    "location",
                    "domains",
                    "pathways",
                    "tissue_pattern",
                )
            }
            if protein
            else None
        ),
        "origin": {k: origin[k] for k in ("origin", "ladder", "species", "one2one")} if origin else None,
        "paralogues": origin["paralogues"] if origin else None,
        "requires": motifs["requires"] if motifs else None,
        "node": _node(symbol, gene, chrom, results_dir),
        "reader": _reader(symbol, chrom, results_dir) or None,
        "elements": elements or None,
        "human_axis": [e for e in elements if e["case"]] or None,
        "duplication": [e["id"] for e in elements if e["duplicated"]] if elements else None,
        "expression": _expression(symbol),
    }
    unknown = [k for k in LAYERS if layers.get(k) is None]
    return {"symbol": symbol, "chrom": chrom, "layers": layers, "unknown": unknown}


def _ident(s: str) -> str:
    import re

    out = re.sub(r"[^A-Za-z0-9_]", "_", s)
    return out if not out[:1].isdigit() else f"g_{out}"


def render(d: dict[str, Any]) -> str:
    """The BioLang-flavoured view with an evidence note per line, and the UNKNOWN residue last."""
    lay = d["layers"]
    sym, chrom = d["symbol"], d["chrom"]
    name = _ident(sym)
    lines = [f"# {sym} on {chrom}, decompiled from what GenomeOS has read; this view infers nothing", ""]
    g = lay["gene"]
    if g:
        lines += [
            f"gene {name} {{",
            f"  locus: {chrom}:{g['start']}-{g['end']}({g['strand']})",
            f"  # {g['type']}, {g['transcripts']} transcripts; canonical {g['canonical']}: "
            f"{g['exons']} exons, {g['cds_segments']} coding segments",
            '  evidence: curated "GENCODE 50"',
        ]
        o = lay["origin"]
        if o:
            lines.append(
                f"  # origin {o['origin']} ({o['ladder']}); orthologues in {o['species']} species, "
                f"{o['one2one']} one-to-one   [curated: Ensembl Compara]"
            )
        if lay["paralogues"] is not None:
            shown = ", ".join(lay["paralogues"][:8]) or "none"
            more = " …" if len(lay["paralogues"]) > 8 else ""
            lines.append(f"  # paralogues: {shown}{more}   [curated: Ensembl Compara]")
        if lay["requires"]:
            factors = ", ".join(_ident(r["factor"]) for r in lay["requires"])
            lines.append(f"  requires: {factors}   # promoter motifs against a shuffle [predicted: JASPAR]")
        if lay["reader"]:
            read = ", ".join(sorted(c for c, st in lay["reader"].items() if st == "read")) or "no cell type"
            silent = ", ".join(sorted(c for c, st in lay["reader"].items() if st == "silent")) or "none"
            lines.append(
                f"  # read in {read}; silent in {silent}   [inferred: promoter open in ENCODE DNase]"
            )
        n = lay["node"]
        if n:
            lines.append(
                f"  domain: {_ident(n['id'])}   # {(n['length'] or 0) / 1e3:.0f} kb, "
                f"{n['coding_genes']} coding genes, {n['enhancers']} enhancer-like elements "
                f"[{(n['evidence'] or '')[:60]}]"
            )
        if lay["protein"]:
            lines.append(f"  produces: {name}_protein")
        lines.append("}")
        pr = lay["protein"]
        if pr:
            doms = pr.get("domains") or []
            more = " …" if len(doms) > 8 else ""
            lines += [
                "",
                f"protein {name}_protein {{",
                f"  accession: {pr.get('accession')}",
                f"  # {pr.get('length')} aa; {pr.get('function') or 'function not recorded'}",
                f"  # location {', '.join(pr.get('location') or []) or 'not recorded'}; "
                f"{pr.get('tissue_pattern') or ''}",
                f"  domains: {', '.join(doms[:8])}{more}",
                f"  pathways: {', '.join((pr.get('pathways') or [])[:6])}",
                '  evidence: curated "UniProt, InterPro, Reactome"',
                "}",
            ]
    for e in lay["elements"] or []:
        axes = []
        if e["mammal_fraction"] is not None:
            axes.append(f"mammals {e['mammal_fraction'] * 100:.0f}%")
        if e["human_fraction"] is not None:
            axes.append(f"people {e['human_fraction'] * 100:.0f}%")
        if e["case"]:
            axes.append(f"case {e['case']}")
        if e["duplicated"]:
            axes.append("inside a segmental duplication")
        lines += [
            "",
            f"element {e['id']} {{",
            f"  locus: {chrom}:{e['start']}-{e['end']}",
            f"  targets: {name}",
            f"  # {e['action']} by {abs(e['log2_fold_change'] or 0):.2f} log2 in {e['tissue']} when deleted "
            f"[predicted: AlphaGenome, {e['confidence']}]",
        ]
        if axes:
            lines.append(f"  # {'; '.join(axes)}   [curated: Zoonomia, gnomAD Gnocchi, UCSC superdups]")
        if e["requires"]:
            lines.append(f"  # motifs: {', '.join(e['requires'])}   [predicted: JASPAR]")
        lines.append("}")
    ex = lay["expression"]
    if ex and ex.get("top"):
        top = ", ".join(f"{t['tissue']} {t['tpm']}" for t in ex["top"])
        lines += ["", f"# expression (GTEx median TPM): {top}"]
    lines += ["", "unknown {"]
    for k in d["unknown"]:
        lines.append(f"  {k}: not read for this locus")
    if not d["unknown"]:
        lines.append("  # every layer GenomeOS reads is present for this locus")
    lines.append("}")
    return "\n".join(lines) + "\n"
