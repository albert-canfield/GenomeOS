# SPDX-License-Identifier: AGPL-3.0-or-later
"""Human trait associations (the GWAS Catalog) inside the attributed elements: does a consequence land here?

An eQTL says an element moves a gene; a reporter says it drives transcription; VISTA says it is an
enhancer in an embryo. None says the element matters to a person. The GWAS Catalog's lead variants
do, weakly and statistically: a lead variant is the marker most associated with a trait in a locus,
and it sits inside the causal element only sometimes (linkage carries the signal a kilobase or more).
So the question here is one of enrichment: do the elements the attribution singles out (a named
target, a strong effect, constraint, a VISTA positive) hold a lead variant more often than the rest,
and which traits? The catalog's association table (600 MB unpacked, 1.2 million rows) is streamed
from its zip once; only the rows inside or within MARGIN of an element are kept under
data/knowledge/gwas, local, with the rows on the same elements shifted SHIFT bases along the
chromosome, which is the chance level the enrichment is read against. The catalog's mapped gene is its
own nearest-gene call, recorded next to the predicted and inferred targets for the reader's eye, not as a
measurement. Evidence: `curated`.
"""

from __future__ import annotations

import io
import json
import time
import urllib.request
import zipfile
from pathlib import Path
from typing import Any

from genomeos.attribution.eqtl import Intervals

GWAS_ZIP = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations-full.zip"
KNOWLEDGE = Path("data/knowledge/gwas")
MARGIN = 1_000  # a lead variant this close counts as landing on the element
SHIFT = 100_000  # the control: the same elements moved this far along the chromosome
SHIFTED = "|shifted"
EVIDENCE = "curated: NHGRI-EBI GWAS Catalog lead variants (GRCh38 positions), all associations"
HIT_COLUMNS = ("chrom", "pos", "rs", "trait", "mapped_gene", "context", "pvalue_mlog", "pubmed", "elements")


def parse_rows(fh):
    """(chrom, pos, rs, trait, mapped_gene, context, -log10 p, pubmed) for rows with one GRCh38 position."""
    header = None
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if header is None:
            header = {k: i for i, k in enumerate(f)}
            continue
        try:
            chrom, pos = f[header["CHR_ID"]], f[header["CHR_POS"]]
        except IndexError:
            continue
        if not chrom or not pos or ";" in pos or "x" in pos or not pos.isdigit():
            continue
        yield (
            f"chr{chrom}" if not chrom.startswith("chr") else chrom,
            int(pos),
            f[header["SNPS"]],
            f[header["DISEASE/TRAIT"]],
            f[header["MAPPED_GENE"]],
            f[header["CONTEXT"]],
            f[header["PVALUE_MLOG"]],
            f[header["PUBMEDID"]],
        )


def distil(
    intervals: Intervals,
    knowledge: Path = KNOWLEDGE,
    source: str | Path = GWAS_ZIP,
    progress=None,
) -> dict[str, Any]:
    """Stream the catalog once; keep the associations that land on an element (± MARGIN)."""
    knowledge.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if isinstance(source, str) and source.startswith(("http://", "https://")):
        req = urllib.request.Request(source, headers={"User-Agent": "GenomeOS/0.1 (stream)"})
        with urllib.request.urlopen(req, timeout=900) as resp:  # noqa: S310
            blob = io.BytesIO(resp.read())
    else:
        blob = io.BytesIO(Path(source).read_bytes())
    scanned = kept = 0
    dest = knowledge / "hits.tsv"
    with zipfile.ZipFile(blob) as zf, dest.open("w") as out:
        out.write("\t".join(HIT_COLUMNS) + "\n")
        name = next(n for n in zf.namelist() if n.endswith(".tsv"))
        with zf.open(name) as raw, io.TextIOWrapper(raw, encoding="utf-8", errors="replace") as fh:
            for chrom, pos, rs, trait, gene, ctx, mlog, pmid in parse_rows(fh):
                scanned += 1
                if progress and scanned % 200_000 == 0:
                    progress(f"GWAS Catalog: {scanned:,} associations read, {kept:,} on elements")
                ids = intervals.at(chrom, pos)
                if ids:
                    kept += 1
                    out.write(
                        f"{chrom}\t{pos}\t{rs}\t{trait}\t{gene}\t{ctx}\t{mlog}\t{pmid}\t{','.join(ids)}\n"
                    )
    summary = {
        "source": str(source),
        "associations_read": scanned,
        "hits": kept,
        "seconds": round(time.time() - t0),
        "elements_indexed": intervals.n,
        "margin_bp": MARGIN,
        "evidence": EVIDENCE,
    }
    (knowledge / "distil_summary.json").write_text(json.dumps(summary, indent=1))
    return summary


def load_hits(knowledge: Path = KNOWLEDGE) -> dict[str, list[dict[str, Any]]]:
    out: dict[str, list[dict[str, Any]]] = {}
    p = knowledge / "hits.tsv"
    if not p.exists():
        return out
    with p.open() as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            row = {
                "chrom": f[0],
                "pos": int(f[1]),
                "rs": f[2],
                "trait": f[3],
                "mapped_gene": f[4],
                "context": f[5],
                "pvalue_mlog": float(f[6]) if f[6] else None,
                "pubmed": f[7],
            }
            for eid in f[8].split(","):
                out.setdefault(eid, []).append(row)
    return out


def index(elements: list[dict[str, Any]], margin: int = MARGIN, shift: int = SHIFT) -> Intervals:
    """Every element with its margin, and a copy shifted along the chromosome as the chance control."""
    iv = Intervals()
    for e in elements:
        iv.add(e["chrom"], max(0, e["start"] - margin), e["end"] + margin, e["key"])
        iv.add(e["chrom"], max(0, e["start"] - margin + shift), e["end"] + margin + shift, e["key"] + SHIFTED)
    return iv.freeze()


def _frac(rows, pred) -> float | None:
    return round(sum(1 for r in rows if pred(r)) / len(rows), 3) if rows else None


def summarise(elements: list[dict[str, Any]], hits: dict[str, list[dict[str, Any]]]) -> dict[str, Any]:
    """Enrichment of lead variants over the attribution's own partitions of the elements."""
    for e in elements:
        e["gwas"] = hits.get(e["key"], [])
        e["gwas_shifted"] = len(hits.get(e["key"] + SHIFTED, []))
    has = lambda e: bool(e["gwas"])  # noqa: E731
    shifted = _frac(elements, lambda e: e["gwas_shifted"] > 0)
    named = [e for e in elements if (e.get("predicted") or {}).get("gene")]
    unnamed = [e for e in elements if not (e.get("predicted") or {}).get("gene") and "predicted" in e]
    strong = [e for e in named if e["predicted"].get("strength") == "strong"]
    cons = [e for e in elements if (e.get("constrained_fraction") or 0) >= 0.2]
    uncons = [
        e for e in elements if e.get("constrained_fraction") is not None and e["constrained_fraction"] < 0.2
    ]
    with_hit = [e for e in elements if e["gwas"]]
    agree_pred = [
        e
        for e in with_hit
        if (e.get("predicted") or {}).get("gene")
        and any(
            (e["predicted"]["gene"] in (h["mapped_gene"] or "").replace(" - ", ",").split(","))
            for h in e["gwas"]
        )
    ]
    agree_inf = [
        e
        for e in with_hit
        if (e.get("inferred") or {}).get("gene")
        and any(
            (e["inferred"]["gene"] in (h["mapped_gene"] or "").replace(" - ", ",").split(","))
            for h in e["gwas"]
        )
    ]
    traits: dict[str, int] = {}
    for e in with_hit:
        for t in {h["trait"] for h in e["gwas"]}:
            traits[t] = traits.get(t, 0) + 1
    out = {
        "elements": len(elements),
        "with_lead_variant": len(with_hit),
        "fraction_with_lead_variant": _frac(elements, has),
        "shifted_control_fraction": shifted,
        "enrichment_over_shifted": (
            round(_frac(elements, has) / shifted, 2) if shifted and _frac(elements, has) is not None else None
        ),
        "associations": sum(len(e["gwas"]) for e in elements),
        "by_partition": {
            "named_target": {"elements": len(named), "with_lead_variant": _frac(named, has)},
            "no_target": {"elements": len(unnamed), "with_lead_variant": _frac(unnamed, has)},
            "strong_effect": {"elements": len(strong), "with_lead_variant": _frac(strong, has)},
            "constrained": {"elements": len(cons), "with_lead_variant": _frac(cons, has)},
            "not_constrained": {"elements": len(uncons), "with_lead_variant": _frac(uncons, has)},
        },
        "catalog_gene_is_the_predicted_target": round(len(agree_pred) / len(with_hit), 3)
        if with_hit
        else None,
        "catalog_gene_is_the_inferred_target": round(len(agree_inf) / len(with_hit), 3) if with_hit else None,
        "top_traits": dict(sorted(traits.items(), key=lambda kv: -kv[1])[:15]),
        "margin_bp": MARGIN,
        "evidence": EVIDENCE,
    }
    if any("status" in e for e in elements):
        pos = [e for e in elements if e.get("status") == "positive"]
        neg = [e for e in elements if e.get("status") == "negative"]
        out["by_partition"]["vista_positive"] = {"elements": len(pos), "with_lead_variant": _frac(pos, has)}
        out["by_partition"]["vista_negative"] = {"elements": len(neg), "with_lead_variant": _frac(neg, has)}
    return out
