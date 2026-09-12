# SPDX-License-Identifier: AGPL-3.0-or-later
"""Genome diff as code review: two people, or one person against the reference, gene by gene.

A VCF diff lists positions. A geneticist wants to know which proteins differ between two people and
what each difference does: the variant's consequence on the protein, whether it is homozygous, what
AlphaMissense predicts for it when the person has been scored, whether ClinVar lists it. This module
takes two people's protein-changing variants (the coding inventory keeps every one under the person as
`coding_variants.json`), normalises the alleles so one indel written two ways compares, and writes the
difference the way a code review reads: per gene, the lines only the first person has (+), only the
second (-), and the count they share. Against the reference, the second side is empty and the page is
the person's own edits. Nothing here leaves the people's directories; the Markdown is printed or
written where the user says.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from genomeos.genome.individuals import ROOT, coding_inventory, normalise_variant

REFERENCE = "reference"
CLASS_RANK = {"likely_pathogenic": 0, "ambiguous": 1, "likely_benign": 3, None: 2}
CONSEQUENCE_RANK = {"nonsense": 0, "frameshift": 0, "start_lost": 0, "stop_lost": 1, "missense": 2}


def _load(root: Path, name: str, fname: str):
    p = root / name / fname
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text())
    except (OSError, json.JSONDecodeError):
        return None


def variants_of(name: str, root: Path | None = None, progress=None) -> list[dict[str, Any]]:
    """The person's protein-changing variants, computed once when the inventory has not been run."""
    root = root or ROOT
    if name == REFERENCE:
        return []
    rows = _load(root, name, "coding_variants.json")
    if rows is None:
        if progress:
            progress(f"{name}: tracing every coding gene once (a few minutes)")
        coding_inventory(name, None, root, progress=progress)
        rows = _load(root, name, "coding_variants.json") or []
    return rows


def _key(v: dict[str, Any]) -> tuple[str, int, str, str]:
    pos, ref, alt = normalise_variant(v["pos"], v["ref"], v["alt"])
    return (v["chrom"], pos, ref, alt)


def _predictions(name: str, root: Path) -> dict[str, dict[str, Any]]:
    """AlphaMissense scores the person has, keyed chrom:pos:ref:alt (the highest across transcripts)."""
    cache = _load(root, name, "alphamissense.json") or {}
    out = {}
    for key, entries in (cache.get("scores") or {}).items():
        if entries:
            out[key] = max(entries, key=lambda e: e["score"])
    return out


def _clinvar(name: str, root: Path) -> dict[str, dict[str, Any]]:
    screen = _load(root, name, "clinvar_screen.json") or {}
    return {f"{h['chrom']}:{h['pos']}:{h['ref']}:{h['alt']}": h for h in screen.get("hits", [])}


def _decorate(v: dict[str, Any], preds, clin) -> dict[str, Any]:
    key = f"{v['chrom']}:{v['pos']}:{v['ref']}:{v['alt']}"
    p = preds.get(key)
    c = clin.get(key)
    return {
        **v,
        "predicted": {"score": p["score"], "class": p["class"]} if p else None,
        "clinvar": {"significance": c["significance"], "conditions": c["conditions"]} if c else None,
    }


def _rank(v: dict[str, Any]) -> tuple:
    pc = (v.get("predicted") or {}).get("class")
    return (
        0 if v.get("clinvar") else 1,
        CONSEQUENCE_RANK.get(v["consequence"], 5),
        CLASS_RANK.get(pc, 2),
        -((v.get("predicted") or {}).get("score") or 0.0),
        v["zygosity"] != "homozygous",
    )


def genome_diff(a: str, b: str = REFERENCE, root: Path | None = None, progress=None) -> dict[str, Any]:
    """Protein-changing variants only one of the two carries, per gene, with consequence and evidence."""
    root = root or ROOT
    va, vb = variants_of(a, root, progress), variants_of(b, root, progress)
    ka = {_key(v): v for v in va}
    kb = {_key(v): v for v in vb}
    pa, pb = _predictions(a, root), (_predictions(b, root) if b != REFERENCE else {})
    ca, cb = _clinvar(a, root), (_clinvar(b, root) if b != REFERENCE else {})
    genes: dict[str, dict[str, Any]] = {}
    for k, v in ka.items():
        g = genes.setdefault(
            v["gene"], {"gene": v["gene"], "chrom": v["chrom"], "only_a": [], "only_b": [], "shared": 0}
        )
        if k in kb:
            g["shared"] += 1
        else:
            g["only_a"].append(_decorate(v, pa, ca))
    for k, v in kb.items():
        if k not in ka:
            g = genes.setdefault(
                v["gene"], {"gene": v["gene"], "chrom": v["chrom"], "only_a": [], "only_b": [], "shared": 0}
            )
            g["only_b"].append(_decorate(v, pb, cb))
    for g in genes.values():
        g["only_a"].sort(key=_rank)
        g["only_b"].sort(key=_rank)
    differing = [g for g in genes.values() if g["only_a"] or g["only_b"]]
    differing.sort(key=lambda g: min([_rank(v) for v in g["only_a"] + g["only_b"]] or [(9,)]))

    def count(side: str, pred) -> int:
        return sum(1 for g in differing for v in g[side] if pred(v))

    lp = lambda v: (v.get("predicted") or {}).get("class") == "likely_pathogenic"  # noqa: E731
    trunc = lambda v: v["consequence"] in ("nonsense", "frameshift", "start_lost", "stop_lost")  # noqa: E731
    return {
        "a": a,
        "b": b,
        "variants_a": len(ka),
        "variants_b": len(kb),
        "shared": sum(g["shared"] for g in genes.values()),
        "only_a": count("only_a", lambda v: True),
        "only_b": count("only_b", lambda v: True),
        "genes_differing": len(differing),
        "genes_with_a_change_in_either": len(genes),
        "likely_pathogenic_only_a": count("only_a", lp),
        "likely_pathogenic_only_b": count("only_b", lp),
        "truncating_only_a": count("only_a", trunc),
        "truncating_only_b": count("only_b", trunc),
        "clinvar_only_a": count("only_a", lambda v: bool(v.get("clinvar"))),
        "clinvar_only_b": count("only_b", lambda v: bool(v.get("clinvar"))),
        "predictions_available": {a: bool(pa), b: bool(pb)},
        "genes": differing,
        "evidence": "measured genotypes; consequences derived on the canonical transcript; predicted: "
        "AlphaMissense where the person has been scored; curated: ClinVar where the person has been screened",
        "note": "alleles normalised before comparison (trimmed, indels left-aligned); a variant only one "
        "side carries is a difference in the call sets as much as in the people when the files come from "
        "different callers",
    }


def _line(sign: str, v: dict[str, Any]) -> str:
    p = v.get("predicted")
    c = v.get("clinvar")
    tags = []
    if p:
        tags.append(f"{p['class'].replace('_', ' ')} {p['score']:.2f}")
    if c:
        tags.append(f"ClinVar {c['significance']}" + (f": {c['conditions'][:50]}" if c["conditions"] else ""))
    change = v["hgvs_p"] or v["consequence"]
    return f"{sign} {change} ({v['consequence'].replace('_', ' ')}, {v['zygosity']}, {v['genotype']})" + (
        "  [" + "; ".join(tags) + "]" if tags else ""
    )


def render(d: dict[str, Any], top: int = 40) -> str:
    """The diff as Markdown, the way a review reads: a summary, then the genes with their + and - lines."""
    a, b = d["a"], d["b"]
    against = "the reference" if b == REFERENCE else b
    lines = [
        f"# {a} against {against}: protein-changing variants",
        "",
        f"{d['variants_a']:,} protein-changing variants in {a}"
        + (f", {d['variants_b']:,} in {b}, {d['shared']:,} shared" if b != REFERENCE else "")
        + f"; {d['genes_differing']:,} genes differ. "
        + f"Only {a}: {d['only_a']:,} variants ({d['truncating_only_a']} truncating, "
        f"{d['likely_pathogenic_only_a']} predicted likely pathogenic, {d['clinvar_only_a']} in ClinVar)"
        + (
            f"; only {b}: {d['only_b']:,} ({d['truncating_only_b']} truncating, "
            f"{d['likely_pathogenic_only_b']} likely pathogenic, {d['clinvar_only_b']} in ClinVar)."
            if b != REFERENCE
            else "."
        ),
        "",
        f"_{d['evidence']}._ {d['note']}",
        "",
    ]
    for g in d["genes"][:top]:
        lines.append(f"## {g['gene']} ({g['chrom']}, {g['shared']} shared)")
        lines.append("")
        lines.append("```diff")
        for v in g["only_a"]:
            lines.append(_line("+", v))
        for v in g["only_b"]:
            lines.append(_line("-", v))
        lines.append("```")
        lines.append("")
    if len(d["genes"]) > top:
        lines.append(f"… and {len(d['genes']) - top:,} more genes.")
    return "\n".join(lines)
