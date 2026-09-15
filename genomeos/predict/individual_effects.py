# SPDX-License-Identifier: AGPL-3.0-or-later
"""A person's regulatory variants, predicted on one gene (AlphaGenome feature d).

The twin question is not "what does this variant do" but "what do this person's
variants do to this gene". For one gene the regulation layer already names the
elements that reach it (its promoter elements, the enhancers inside its CTCF
node, the ones AlphaGenome itself said move the gene). This module collects the
individual's variants inside those elements, asks the variant scorer what each
one does to the gene's expression per tissue, and sums the effects per
haplotype where the genotypes are phased. One request per variant, about eight
seconds, cached per variant under data/knowledge/alphagenome/variants; a gene
with sixty regulatory variants is eight minutes once and instant after.

The model's own caveats apply and are repeated in the output: it reads one
unphased sequence at a time, so a haplotype sum is a sum of single-variant
effects, not a prediction on the whole haplotype; personal genomes are a known
weakness; everything is `predicted`, capped at 0.7.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

CACHE = Path("data/knowledge/alphagenome/variants")
MIN_EFFECT = 0.1
MAX_INDEL = 50
Scorer = Callable[
    [str, int, str, str], list[tuple[str, str, float]]
]  # chrom, pos1, ref, alt -> (gene, tissue, lfc)


def intervals_for(reg: dict[str, Any]) -> list[dict[str, Any]]:
    """The elements that reach the gene, most likely first: promoters, enhancers the model already
    tied to this gene, near enhancers, the rest of the node."""
    out = []
    for e in reg.get("promoters", []):
        out.append({**e, "priority": 0})
    for e in reg.get("enhancers", []):
        if e.get("predicted_this_gene"):
            pr = 1
        elif e["distance"] <= 20_000:
            pr = 2
        else:
            pr = 3
        out.append({**e, "priority": pr})
    out.sort(key=lambda e: (e["priority"], e["distance"]))
    return out


def collect(vcf: Path, intervals: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """The person's variants inside the elements (1-based positions, first alt of each row)."""
    from genomeos.genome.individuals import rows_in

    seen: set[tuple[int, str, str]] = set()
    out = []
    for e in intervals:
        for f in rows_in(vcf, e["start"] + 1, e["end"]):
            ref, alt = f[3], f[4].split(",")[0]
            if len(ref) > MAX_INDEL or len(alt) > MAX_INDEL or not set(ref + alt) <= set("ACGTN"):
                continue
            key = (int(f[1]), ref, alt)
            if key in seen:
                continue
            seen.add(key)
            out.append(
                {
                    "pos": int(f[1]),
                    "ref": ref,
                    "alt": alt,
                    "genotype": f[9].split(":")[0] if len(f) > 9 else "",
                    "element": e["id"],
                    "element_class": e["class"],
                    "distance_to_tss": e["distance"],
                    "priority": e["priority"],
                }
            )
    return out


def cache_path(chrom: str, pos: int, ref: str, alt: str, cache: Path = CACHE) -> Path:
    return cache / chrom / f"{pos}_{ref[:20]}_{alt[:20]}.json"


def score_variant(
    scorer: Scorer, chrom: str, pos: int, ref: str, alt: str, cache: Path = CACHE
) -> list[list[Any]]:
    """Every (gene, tissue, log2FC) the scorer returns for the variant, cached."""
    p = cache_path(chrom, pos, ref, alt, cache)
    if p.exists():
        return json.loads(p.read_text())
    effects = [[g, t, round(float(v), 4)] for g, t, v in scorer(chrom, pos, ref, alt)]
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(effects))
    return effects


def _on_gene(effects: list[list[Any]], gene: str) -> dict[str, Any]:
    mine = [(t, v) for g, t, v in effects if g == gene]
    if not mine:
        return {"log2_fold_change": 0.0, "tissue": None, "tracks_moved": 0}
    t, v = max(mine, key=lambda x: abs(x[1]))
    return {
        "log2_fold_change": v,
        "tissue": t,
        "tracks_moved": sum(1 for _, x in mine if abs(x) >= MIN_EFFECT),
        "mean_log2_fold_change": round(sum(x for _, x in mine) / len(mine), 4),
    }


def haplotype_sums(rows: list[dict[str, Any]], min_effect: float = MIN_EFFECT) -> dict[str, Any]:
    """Sum of single-variant effects per haplotype (phased), with unphased heterozygotes listed apart."""
    h1 = h2 = 0.0
    unphased: list[float] = []
    for r in rows:
        v = r["effect"]["log2_fold_change"]
        if abs(v) < min_effect:
            continue
        gt = r["genotype"].replace("/", "|") if r["genotype"].count("1") == 2 else r["genotype"]
        if "|" in gt:
            a, _, b = gt.partition("|")
            h1 += v if a == "1" else 0.0
            h2 += v if b == "1" else 0.0
        elif gt.count("1") == 1:
            unphased.append(v)
        else:
            h1 += v
            h2 += v
    return {
        "hap1_log2_fold_change": round(h1, 3),
        "hap2_log2_fold_change": round(h2, 3),
        "unphased_heterozygous_effects": [round(x, 3) for x in unphased],
        "note": "sum of single-variant predictions; the model reads one variant at a time, so this is not "
        "a prediction on the whole haplotype",
    }


def predict_gene(
    name: str,
    gene: str,
    chrom: str,
    reg: dict[str, Any],
    scorer: Scorer,
    vcf: Path,
    max_variants: int = 60,
    cache: Path = CACHE,
    progress=None,
    min_effect: float = MIN_EFFECT,
) -> dict[str, Any]:
    """Collect the person's variants in the elements that reach the gene and score them one by one."""
    intervals = intervals_for(reg)
    found = collect(vcf, intervals)
    todo = found[:max_variants]
    rows = []
    for i, v in enumerate(todo, 1):
        effects = score_variant(scorer, chrom, v["pos"], v["ref"], v["alt"], cache)
        on = _on_gene(effects, gene)
        others = [(g, t, x) for g, t, x in effects if g != gene and abs(x) >= min_effect]
        strongest_other = max(others, key=lambda e: abs(e[2])) if others else None
        rows.append(
            {
                **v,
                "effect": on,
                "moves_gene": abs(on["log2_fold_change"]) >= min_effect,
                "strongest_other_gene": (
                    {
                        "gene": strongest_other[0],
                        "tissue": strongest_other[1],
                        "log2_fold_change": strongest_other[2],
                    }
                    if strongest_other
                    else None
                ),
            }
        )
        if progress:
            progress(f"{gene}: {i}/{len(todo)} variants scored ({v['element_class']} {v['pos']:,})")
    rows.sort(key=lambda r: -abs(r["effect"]["log2_fold_change"]))
    moving = [r for r in rows if r["moves_gene"]]
    return {
        "individual": name,
        "gene": gene,
        "chrom": chrom,
        "elements_considered": len(intervals),
        "variants_in_elements": len(found),
        "variants_scored": len(rows),
        "variants_moving_gene": len(moving),
        "strongest": rows[0] if rows else None,
        "haplotypes": haplotype_sums(rows, min_effect),
        "variants": rows,
        "min_effect_log2fc": min_effect,
        "evidence": "measured genotypes; predicted effect per variant (AlphaGenome RNA-seq gene scorer), "
        "confidence capped at 0.7",
    }
