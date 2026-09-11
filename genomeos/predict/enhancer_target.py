# SPDX-License-Identifier: AGPL-3.0-or-later
"""Enhancer to gene, predicted (AlphaGenome feature b).

ENCODE says where an enhancer-like element is; nothing in the registry says
which gene it reaches. `genomeos regulation` infers the target as the nearest
transcription start inside the CTCF domain (inferred, 0.4). This module asks
AlphaGenome a sharper question: delete the element from its 1 Mb window and
read the predicted expression change of every gene in the window, per tissue
track. The gene whose expression drops most is the predicted target, the track
where it drops most is the tissue, and the drop is the magnitude. A rise means
the element behaves as a silencer for that gene.

Everything here is `predicted` evidence, capped at confidence 0.7 (docs/
ALPHAGENOME.md). One request per element, about eight seconds, so a chromosome
is a sampled job and the per-element answers are cached under
data/knowledge/alphagenome/elements (local, not committed); the summary that is
committed says how often the prediction agrees with the domain inference, which
is a test of the node model itself.
"""

from __future__ import annotations

import json
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

from genomeos.coords import Locus

CACHE = Path("data/knowledge/alphagenome/elements")
MIN_EFFECT = 0.1  # smallest |log2 fold change| that names a target
STRONG_EFFECT = 0.3
CONFIDENCE_CAP = 0.7

Scorer = Callable[[str, int, str, str], list[tuple[str, str, float]]]


def cache_path(chrom: str, element_id: str, cache: Path = CACHE) -> Path:
    return cache / chrom / f"{element_id}.json"


def load_cached(chrom: str, element_id: str, cache: Path = CACHE) -> dict[str, Any] | None:
    p = cache_path(chrom, element_id, cache)
    return json.loads(p.read_text()) if p.exists() else None


def aggregate(effects: list[tuple[str, str, float]]) -> list[dict[str, Any]]:
    """Per gene over every track: mean, the largest drop and the largest rise with their tissues."""
    by_gene: dict[str, dict[str, Any]] = {}
    for gene, tissue, val in effects:
        g = by_gene.setdefault(
            gene,
            {
                "gene": gene,
                "n_tracks": 0,
                "sum": 0.0,
                "min": 0.0,
                "min_tissue": "",
                "max": 0.0,
                "max_tissue": "",
            },
        )
        g["n_tracks"] += 1
        g["sum"] += val
        if val < g["min"]:
            g["min"], g["min_tissue"] = val, tissue
        if val > g["max"]:
            g["max"], g["max_tissue"] = val, tissue
    rows = []
    for g in by_gene.values():
        rows.append(
            {
                "gene": g["gene"],
                "n_tracks": g["n_tracks"],
                "mean_log2fc": round(g["sum"] / g["n_tracks"], 4),
                "max_drop_log2fc": round(g["min"], 4),
                "max_drop_tissue": g["min_tissue"],
                "max_rise_log2fc": round(g["max"], 4),
                "max_rise_tissue": g["max_tissue"],
            }
        )
    rows.sort(key=lambda r: r["max_drop_log2fc"])
    return rows


def predict_target(rows: list[dict[str, Any]], min_effect: float = MIN_EFFECT) -> dict[str, Any] | None:
    """The gene the deleted element reaches, or None when no gene moves by `min_effect`.

    Loss of expression on deletion means the element activates the gene; a rise means it
    represses it. The larger of the two wins; ties go to activation, the usual case.
    """
    best = None
    for r in rows:
        drop, rise = -r["max_drop_log2fc"], r["max_rise_log2fc"]
        if drop >= rise:
            cand = (drop, "activates", r["max_drop_tissue"], -drop)
        else:
            cand = (rise, "represses", r["max_rise_tissue"], rise)
        if cand[0] >= min_effect and (best is None or cand[0] > best[0]):
            best = (*cand, r["gene"])
    if best is None:
        return None
    size, action, tissue, log2fc, gene = best
    return {
        "gene": gene,
        "action": action,
        "log2_fold_change": round(log2fc, 4),
        "tissue": tissue,
        "strength": "strong" if size >= STRONG_EFFECT else "weak",
        "confidence": round(min(CONFIDENCE_CAP, size), 3),
        "basis": "predicted: expression change on deleting the element (AlphaGenome RNA-seq gene scorer)",
    }


def compare(
    prediction: dict[str, Any] | None,
    inferred_targets: list[dict[str, Any]],
    domain_genes: list[str] | None = None,
) -> str:
    """How the prediction relates to the domain inference: one of four verdicts."""
    if prediction is None:
        return "no predicted effect"
    genes = [t["gene"] for t in inferred_targets]
    if genes and prediction["gene"] == genes[0]:
        return "agrees with nearest TSS in domain"
    if prediction["gene"] in genes or prediction["gene"] in (domain_genes or []):
        return "another gene in the same domain"
    return "gene outside the domain"


def score_element(
    scorer: Scorer,
    fetch: Callable[[Locus], str],
    chrom: str,
    element_id: str,
    start: int,
    end: int,
    inferred_targets: list[dict[str, Any]] | None = None,
    domain_genes: list[str] | None = None,
    coding: set[str] | None = None,
    cache: Path = CACHE,
    min_effect: float = MIN_EFFECT,
) -> dict[str, Any]:
    """Delete one element (0-based half-open start..end) as a variant and read every gene in the window.

    Cached per element id; a cached answer is returned without a request. `predicted` names the gene
    that moves most, whatever its type; `predicted_coding` restricts the choice to `coding` symbols,
    which is the comparison the domain inference (coding genes only) can be held to.
    """
    hit = load_cached(chrom, element_id, cache)
    if hit is None:
        seq = str(fetch(Locus(chrom, start - 1, end))).upper()  # VCF-style: anchor base then the element
        t0 = time.time()
        effects = scorer(chrom, start, seq, seq[0])
        rows = aggregate(effects)
        hit = {
            "id": element_id,
            "chrom": chrom,
            "start": start,
            "end": end,
            "length": end - start,
            "genes_in_window": len(rows),
            "tracks": rows[0]["n_tracks"] if rows else 0,
            "seconds": round(time.time() - t0, 1),
            "genes": rows,
        }
        p = cache_path(chrom, element_id, cache)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(hit))
    pred = predict_target(hit["genes"], min_effect)
    pred_coding = (
        predict_target([g for g in hit["genes"] if g["gene"] in coding], min_effect) if coding else pred
    )
    out = {k: v for k, v in hit.items() if k != "genes"}
    out["predicted"] = pred
    out["predicted_coding"] = pred_coding
    out["top_genes"] = hit["genes"][:5]
    if inferred_targets is not None:
        out["inferred"] = inferred_targets[0] if inferred_targets else None
        out["verdict"] = compare(pred, inferred_targets, domain_genes)
        out["verdict_coding"] = compare(pred_coding, inferred_targets, domain_genes)
    return out


def summarise(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """The committed summary of a chromosome's sampled elements."""
    verdicts: dict[str, int] = {}
    verdicts_coding: dict[str, int] = {}
    tissues: dict[str, int] = {}
    strong = weak = repress = 0
    distances = []
    for r in rows:
        verdicts[r.get("verdict", "?")] = verdicts.get(r.get("verdict", "?"), 0) + 1
        vc = r.get("verdict_coding", r.get("verdict", "?"))
        verdicts_coding[vc] = verdicts_coding.get(vc, 0) + 1
        p = r.get("predicted")
        if p:
            tissues[p["tissue"]] = tissues.get(p["tissue"], 0) + 1
            strong += p["strength"] == "strong"
            weak += p["strength"] == "weak"
            repress += p["action"] == "represses"
        if r.get("inferred") and p and r.get("verdict") == "agrees with nearest TSS in domain":
            distances.append(r["inferred"]["distance"])
    n = len(rows)
    named = strong + weak
    agree = verdicts.get("agrees with nearest TSS in domain", 0)
    in_domain = agree + verdicts.get("another gene in the same domain", 0)
    return {
        "elements_scored": n,
        "with_predicted_target": named,
        "fraction_with_target": round(named / n, 3) if n else 0.0,
        "strong": strong,
        "weak": weak,
        "silencer_like": repress,
        "verdicts": dict(sorted(verdicts.items(), key=lambda kv: -kv[1])),
        "verdicts_coding": dict(sorted(verdicts_coding.items(), key=lambda kv: -kv[1])),
        "coding_target_agrees_with_nearest": round(
            verdicts_coding.get("agrees with nearest TSS in domain", 0)
            / max(1, n - verdicts_coding.get("no predicted effect", 0)),
            3,
        ),
        "fraction_agreeing_with_nearest": round(agree / named, 3) if named else None,
        "fraction_inside_domain": round(in_domain / named, 3) if named else None,
        "top_tissues": dict(sorted(tissues.items(), key=lambda kv: -kv[1])[:12]),
        "median_distance_when_agreeing": sorted(distances)[len(distances) // 2] if distances else None,
    }


class Context:
    """One chromosome's elements, domains and genes, loaded once for many predictions."""

    def __init__(self, chrom: str, reference: Path = Path("data/reference")) -> None:
        from genomeos.genome import Annotation, IndexedGenome, default_gencode
        from genomeos.genome.domains import infer_domains
        from genomeos.genome.regulation import assign_targets
        from genomeos.genome.regulatory import load_ccres

        self.chrom = chrom
        self.ccres = load_ccres(chrom)
        fasta = reference / f"{chrom}.fa"
        if not self.ccres or not fasta.exists():
            raise FileNotFoundError(
                f"{chrom} is not fetched; run `genomeos data fetch --chrom {chrom}` first"
            )
        self.genome = IndexedGenome(str(fasta))
        self.annotation = Annotation.from_gff3(default_gencode({chrom}), {chrom})
        self.domains = infer_domains(chrom, self.genome.lengths[chrom], self.ccres, self.annotation)
        self.domain_by_id = {d.id: d for d in self.domains}
        self.elements = assign_targets(chrom, self.ccres, self.annotation, self.domains)
        self.element_by_id = {e.id: e for e in self.elements}
        self.coding = {g.symbol for g in self.annotation.genes.values() if g.type == "protein_coding"}
        self._distal_ids = {c.id for c in self.ccres if c.cls == "dELS"}

    def close(self) -> None:
        self.genome.close()

    def distal_enhancers(self, min_distance: int = 2_000, max_distance: int = 400_000) -> list:
        """Distal enhancer-like elements whose inferred target is clear of the TSS and inside the window."""
        out = [
            e
            for e in self.elements
            if e.cls == "enhancer"
            and e.id in self._distal_ids
            and e.targets
            and min_distance <= e.targets[0]["distance"] <= max_distance
        ]
        out.sort(key=lambda e: e.locus.start)
        return out

    def element_at(self, start: int, end: int):
        """The ENCODE element overlapping start..end, or None."""
        for e in self.elements:
            if e.locus.start < end and start < e.locus.end:
                return e
        return None

    def score(
        self, scorer: Scorer, element, cache: Path = CACHE, min_effect: float = MIN_EFFECT
    ) -> dict[str, Any]:
        dom = self.domain_by_id.get(element.domain)
        r = score_element(
            scorer,
            self.genome.fetch,
            self.chrom,
            element.id,
            element.locus.start,
            element.locus.end,
            element.targets,
            list(dom.genes) if dom is not None else [],
            coding=self.coding,
            cache=cache,
            min_effect=min_effect,
        )
        r["domain"] = element.domain
        return r

    def score_region(
        self, scorer: Scorer, start: int, end: int, cache: Path = CACHE, min_effect: float = MIN_EFFECT
    ) -> dict[str, Any]:
        """Score an ENCODE element if one overlaps the region, else the region itself as an ad-hoc element."""
        e = self.element_at(start, end)
        if e is not None:
            return self.score(scorer, e, cache, min_effect)
        r = score_element(
            scorer,
            self.genome.fetch,
            self.chrom,
            f"{self.chrom}_{start}_{end}",
            start,
            end,
            None,
            None,
            coding=self.coding,
            cache=cache,
            min_effect=min_effect,
        )
        r["domain"] = ""
        return r


def cached_prediction(chrom: str, element_id: str, cache: Path = CACHE) -> dict[str, Any] | None:
    """The predicted target of an element already scored, for layers that never call the API."""
    hit = load_cached(chrom, element_id, cache)
    return predict_target(hit["genes"]) if hit else None
