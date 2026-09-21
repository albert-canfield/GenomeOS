# SPDX-License-Identifier: AGPL-3.0-or-later
"""What scoring every element of a chromosome bought the 98%, decided on chromosome 21.

Chromosome 21 is the one chromosome whose 12,139 registry elements have all been deleted in
AlphaGenome with the effect read per gene and per cell line (5,409 requests). Before more quota
goes to the sweep, this module asks what that bought for the UNKNOWN space, from the tables
already on disk and with no model call:

1. the ablation, the only test that can falsify: the gene-level closure recomputed from the
   committed `closure_<chrom>` table with the elements that overlap an UNKNOWN block dropped from
   every gene's input, and again with only those elements, each against random removals of the
   same size per gene, so "the closure fell" cannot be the amount of input removed;
2. whether elements over unknown blocks differ from the rest at all, on movement, magnitude,
   silencer share, how many cell lines they act in and the distance to their target, raw and
   inside strata of length, GC and distance to the nearest coding TSS, because an unknown block
   is gene-poor sequence by construction and that alone moves every one of those numbers;
3. whether the syntax reading gains: the block's case on the two axes (area J) and its mammalian
   constraint against the model's effects, and what the scoring now says about the 69 syntax
   candidates of `syntax_candidates_genome_wide` that lie on a scored chromosome;
4. validation against measurement where chr21 has any: lentiMPRA activity per cell, VISTA;
5. the coverage the scoring added to the unknown space, block by block and base by base, against
   what the two samples had already named, with the request cost beside it.

Every comparison states its control and its null. A difference that does not survive its strata
is reported as not surviving.
"""

from __future__ import annotations

import json
import random
import time
from collections import Counter
from pathlib import Path
from typing import Any

from genomeos.attribution.closure import attributed_elements, judge

# imported rather than restated: this scoring writes it into its result as `min_effect_log2`,
# so two homes would publish two bars under one name
from genomeos.predict.enhancer_target import MIN_EFFECT
from genomeos.results import RESULTS_DIR, load_result, save_result

LENGTH_BINS = 3
GC_BINS = 4
TSS_BINS = 40  # distance is the covariate an unknown block shifts most: 8 bins left it 260 kb against 173 kb
PERMUTATIONS = 2000
RANDOM_DRAWS = 50  # random removals of the same size, the control for the ablation
SEED = 20260913
REFERENCE = Path("data/reference")
EVIDENCE = {
    "elements": "predicted: AlphaGenome deletion per element, per gene and per cell line (already scored)",
    "blocks": "inferred: genomeos unknown classes and the budget's tiers; curated sequence classes",
    "axes": "inferred: Zoonomia phyloP against gnomAD Gnocchi per block (area J)",
    "closure": "experimental: ENCODE total RNA-seq and DNase (the committed closure table)",
    "measured": "experimental: ENCODE4 lentiMPRA, VISTA e11.5",
}


# ---- elements and their blocks ----------------------------------------------------------


def scored_elements(chrom: str, results_dir: Path = RESULTS_DIR) -> list[dict[str, Any]]:
    """Every element of the whole-chromosome run, which must be complete."""
    from genomeos.attribution.targets import run_elements

    summary = load_result(f"enhancer_targets_all_{chrom}", results_dir) or {}
    if not summary.get("complete"):
        raise FileNotFoundError(
            f"enhancer_targets_all_{chrom} is not complete; this reading is for a fully scored chromosome"
        )
    return run_elements("enhancer_targets_all", chrom, results_dir)


def unknown_blocks(chrom: str, results_dir: Path = RESULTS_DIR) -> list[dict[str, Any]]:
    """The UNKNOWN blocks with their tier, their mammalian constraint and their case on the two axes."""
    unknown = load_result(f"unknown_{chrom}", results_dir) or {}
    budget = {
        (b["start"], b["end"]): b
        for b in (load_result(f"budget_{chrom}", results_dir) or {}).get("blocks", [])
    }
    var = {
        (b["start"], b["end"]): b
        for b in (load_result(f"variation_{chrom}", results_dir) or {}).get("blocks", [])
    }
    out = []
    for b in unknown.get("blocks", []):
        key = (b["start"], b["end"])
        bu, v = budget.get(key) or {}, var.get(key) or {}
        out.append(
            {
                "start": b["start"],
                "end": b["end"],
                "length": b["length"],
                "class": b["class"],
                "tier": (bu.get("guess") or {}).get("tier"),
                "mammal_fraction": (bu.get("phylop") or {}).get("fraction_above"),
                "case": (v.get("case") or {}).get("case"),
                "human_fraction": (v.get("gnocchi") or {}).get("fraction_above"),
            }
        )
    out.sort(key=lambda b: b["start"])
    return out


def gc_fraction(seq: str) -> float | None:
    acgt = sum(seq.count(b) for b in "ACGT")
    return round((seq.count("G") + seq.count("C")) / acgt, 4) if acgt else None


def annotate(
    chrom: str, elements: list[dict[str, Any]], results_dir: Path = RESULTS_DIR
) -> list[dict[str, Any]]:
    """Per element: the UNKNOWN block it overlaps (edges count), length, GC, distance to the nearest
    coding TSS and to its predicted target, and what the deletion said."""
    from genomeos.coords import Locus
    from genomeos.genome import Annotation
    from genomeos.genome.annotation import default_gencode
    from genomeos.genome.index import IndexedGenome

    blocks = unknown_blocks(chrom, results_dir)
    starts = [b["start"] for b in blocks]
    ann = Annotation.from_gff3(default_gencode({chrom}), {chrom})
    tss_of = {
        g.symbol: (g.locus.end - 1 if g.locus.strand.value == "-" else g.locus.start)
        for g in ann.genes.values()
        if g.locus.chrom == chrom
    }
    coding_tss = sorted(
        (g.locus.end - 1 if g.locus.strand.value == "-" else g.locus.start)
        for g in ann.protein_coding()
        if g.locus.chrom == chrom
    )
    genome = IndexedGenome(REFERENCE / f"{chrom}.fa")
    import bisect

    bodies: list[list[int]] = []
    for lo, hi in sorted((g.locus.start, g.locus.end) for g in ann.genes.values() if g.locus.chrom == chrom):
        if bodies and lo <= bodies[-1][1]:
            bodies[-1][1] = max(bodies[-1][1], hi)
        else:
            bodies.append([lo, hi])
    body_starts = [b[0] for b in bodies]

    rows = []
    for e in elements:
        s, t = e["start"], e["end"]
        i = max(0, bisect.bisect_right(starts, s) - 1)
        best, best_bp = None, 0
        while i < len(blocks) and blocks[i]["start"] < t:
            ov = min(blocks[i]["end"], t) - max(blocks[i]["start"], s)
            if ov > best_bp:
                best, best_bp = blocks[i], ov
            i += 1
        mid = (s + t) // 2
        j = bisect.bisect_left(coding_tss, mid)
        near = min(
            (abs(coding_tss[k] - mid) for k in (j - 1, j) if 0 <= k < len(coding_tss)),
            default=None,
        )
        pc = e.get("predicted_coding") or {}
        pred = e.get("predicted") or {}
        by_cell = e.get("predicted_by_cell") or {}
        rows.append(
            {
                "id": e["id"],
                "start": s,
                "end": t,
                "length": t - s,
                "gc": gc_fraction(str(genome.fetch(Locus(chrom, s, t))).upper()),
                "nearest_coding_tss": near,
                "in_gene_body": any(
                    bodies[k][0] < t and bodies[k][1] > s
                    for k in range(
                        max(0, bisect.bisect_right(body_starts, s) - 1), bisect.bisect_left(body_starts, t)
                    )
                ),
                "in_unknown": best is not None,
                "unknown_bp": best_bp,
                "block": [best["start"], best["end"]] if best else None,
                "tier": best["tier"] if best else None,
                "block_class": best["class"] if best else None,
                "case": best["case"] if best else None,
                "block_mammal_fraction": best["mammal_fraction"] if best else None,
                "moves_gene": bool(pred),
                "names_coding": bool(pc),
                "abs_log2": round(abs(pred["log2_fold_change"]), 4) if pred else None,
                "silencer": (pred.get("action") == "represses") if pred else None,
                # among elements that move a gene: in how many of the four lines the deletion reaches the bar
                "cells_acting": (
                    sum(1 for v in by_cell.values() if v is not None and abs(v) >= MIN_EFFECT)
                    if pred
                    else None
                ),
                "target": pc.get("gene"),
                "target_tss_distance": (abs(tss_of[pc["gene"]] - mid) if pc.get("gene") in tss_of else None),
            }
        )
    genome.close()
    return rows


# ---- strata and the comparison ------------------------------------------------------------


def _bins(values: list[float], n: int) -> list[float]:
    xs = sorted(v for v in values if v is not None)
    return [xs[int(round(k * (len(xs) - 1) / n))] for k in range(1, n)] if xs else []


def _bin_of(value: float | None, edges: list[float]) -> int:
    if value is None:
        return -1
    return sum(1 for e in edges if value >= e)


def add_strata(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Coarsened exact matching: length, GC and distance to the nearest coding TSS in bins."""
    le = _bins([r["length"] for r in rows], LENGTH_BINS)
    ge = _bins([r["gc"] for r in rows], GC_BINS)
    te = _bins([r["nearest_coding_tss"] for r in rows], TSS_BINS)
    for r in rows:
        r["stratum"] = (
            _bin_of(r["length"], le),
            _bin_of(r["gc"], ge),
            _bin_of(r["nearest_coding_tss"], te),
        )
    return rows


def _mean(values: list[float]) -> float | None:
    vs = [v for v in values if v is not None]
    return sum(vs) / len(vs) if vs else None


def stratified(
    rows: list[dict[str, Any]], metric: str, flag: str = "in_unknown", seed: int = SEED
) -> dict[str, Any]:
    """One metric, unknown against the rest: raw means, means inside shared strata, and the
    permutation p of the stratified difference from shuffling the flag within each stratum."""
    usable = [r for r in rows if r.get(metric) is not None]
    a = [r for r in usable if r[flag]]
    b = [r for r in usable if not r[flag]]
    by_stratum: dict[tuple, dict[str, list[float]]] = {}
    for r in usable:
        d = by_stratum.setdefault(r["stratum"], {"a": [], "b": []})
        d["a" if r[flag] else "b"].append(r[metric])
    shared = {k: v for k, v in by_stratum.items() if v["a"] and v["b"]}
    weight = {k: len(v["a"]) for k, v in shared.items()}
    total = sum(weight.values()) or 1

    def matched(side: str) -> float | None:
        if not shared:
            return None
        return sum(weight[k] * (_mean(v[side]) or 0.0) for k, v in shared.items()) / total

    obs_a, obs_b = matched("a"), matched("b")
    diff = (obs_a - obs_b) if obs_a is not None and obs_b is not None else None
    rng = random.Random(seed)
    null = []
    if diff is not None:
        pooled = {k: v["a"] + v["b"] for k, v in shared.items()}
        for _ in range(PERMUTATIONS):
            da = db = 0.0
            for k, v in pooled.items():
                rng.shuffle(v)
                na = len(shared[k]["a"])
                da += weight[k] * (_mean(v[:na]) or 0.0)
                db += weight[k] * (_mean(v[na:]) or 0.0)
            null.append(da / total - db / total)
    p = (sum(1 for x in null if abs(x) >= abs(diff)) + 1) / (len(null) + 1) if null else None
    return {
        "metric": metric,
        "unknown_n": len(a),
        "rest_n": len(b),
        "unknown_raw": round(_mean([r[metric] for r in a]), 4) if a else None,
        "rest_raw": round(_mean([r[metric] for r in b]), 4) if b else None,
        "unknown_matched": round(obs_a, 4) if obs_a is not None else None,
        "rest_matched": round(obs_b, 4) if obs_b is not None else None,
        "difference_matched": round(diff, 4) if diff is not None else None,
        "p_permuted_within_strata": round(p, 4) if p is not None else None,
        "strata_shared": len(shared),
        "unknown_in_shared_strata": total,
    }


METRICS = ("moves_gene", "names_coding", "abs_log2", "silencer", "cells_acting", "target_tss_distance")


def comparison(rows: list[dict[str, Any]]) -> dict[str, Any]:
    out = {m: stratified(rows, m) for m in METRICS}
    out["balance_inside_strata"] = {
        c: {k: v[k] for k in ("unknown_matched", "rest_matched")}
        for c, v in ((c, stratified(rows, c)) for c in ("length", "gc", "nearest_coding_tss"))
    }
    out["covariates"] = {
        side: {
            "n": len(sel),
            "length": round(_mean([r["length"] for r in sel]), 1),
            "gc": round(_mean([r["gc"] for r in sel]), 4),
            "nearest_coding_tss": round(_mean([r["nearest_coding_tss"] for r in sel]), 1),
            "overlapping_a_gene_body": round(_mean([r["in_gene_body"] for r in sel]), 4),
        }
        for side, sel in (
            ("unknown", [r for r in rows if r["in_unknown"]]),
            ("rest", [r for r in rows if not r["in_unknown"]]),
        )
    }
    return out


# ---- the ablation ------------------------------------------------------------------------


def _rebuild(closure_result: dict, by_gene: dict[str, list[dict]], indexes: dict, keep) -> list[dict]:
    """The closure's gene by cell table rebuilt with only the elements `keep` admits."""
    rows = []
    for r in closure_result["genes"]:
        els = [e for e in by_gene.get(r["gene"], []) if keep(e)]
        row = {"gene": r["gene"], "elements": len(els), "cells": {}}
        for cell in closure_result["cells"]:
            x = r["cells"][cell]
            idx = indexes[cell]
            active = [e for e in els if idx.overlapping(e["start"], e["end"])]
            scored = [e for e in els if cell in e["by_cell"]]
            open_scored = [e for e in active if cell in e["by_cell"]]
            row["cells"][cell] = {
                "expression": x["expression"],
                "expressed": x["expressed"],
                "promoter_open": x["promoter_open"],
                "promoter_percentile": x["promoter_percentile"],
                "active_elements": len(active),
                "input_both": round(sum(e["by_cell"][cell] for e in open_scored), 3) if scored else None,
                "scored_elements": len(scored),
            }
        rows.append(row)
    return rows


def _across(rows: list[dict], cells: list[str], permutations: int = 1000, seed: int = 0) -> dict[str, Any]:
    j = judge(rows, cells, seed, "input_both", permutations)
    a = j["across_cells"]
    return {
        "genes_tested": a["genes_tested"],
        "most_active_cell_is_most_expressed": a["most_active_cell_is_most_expressed"],
        "shuffled": a["same_under_shuffled_cells"],
        "p_value": a["p_value"],
        "mean_rho": a["mean_rho_input_vs_expression"],
        "rejected": j["rejected_count"],
        "unexplained": j["unexplained_count"],
    }


def ablation(chrom: str, results_dir: Path = RESULTS_DIR, rows: list[dict] | None = None) -> dict[str, Any]:
    """The closure with and without the elements that overlap an UNKNOWN block, against random
    removals of the same size per gene. The full-input row must reproduce the committed closure."""
    from genomeos.genome.reader import PeakIndex, load_peaks

    c = load_result(f"closure_{chrom}", results_dir)
    if not c or not c.get("genes"):
        raise FileNotFoundError(f"no closure_{chrom} with a gene table; run genomeos closure --chrom {chrom}")
    flags = {r["id"]: r["in_unknown"] for r in (rows or [])}
    by_gene = attributed_elements(chrom)
    indexes = {cell: PeakIndex(load_peaks(cell, chrom)) for cell in c["cells"]}
    unknown_count = {g: sum(1 for e in els if flags.get(e["id"])) for g, els in by_gene.items()}
    total = sum(len(v) for v in by_gene.values())
    in_unknown = sum(unknown_count.values())

    def run(keep, permutations=1000):
        return _across(_rebuild(c, by_gene, indexes, keep), c["cells"], permutations)

    full = run(lambda e: True)
    without = run(lambda e: not flags.get(e["id"]))
    only = run(lambda e: bool(flags.get(e["id"])))
    rng = random.Random(SEED)

    def draws(drop: bool, observed: float | None) -> dict[str, Any]:
        hits, genes = [], []
        for _ in range(RANDOM_DRAWS):
            chosen: set[str] = set()
            for g, els in by_gene.items():
                k = unknown_count.get(g, 0)
                if drop:
                    chosen.update(e["id"] for e in rng.sample(els, k))
                else:
                    chosen.update(e["id"] for e in rng.sample(els, len(els) - k))
            out = run(lambda e, chosen=chosen: e["id"] not in chosen, permutations=0)
            if out["most_active_cell_is_most_expressed"] is not None:
                hits.append(out["most_active_cell_is_most_expressed"])
                genes.append(out["genes_tested"])
        return {
            "draws": len(hits),
            "mean_hit": round(sum(hits) / len(hits), 4) if hits else None,
            "min_hit": min(hits) if hits else None,
            "max_hit": max(hits) if hits else None,
            "mean_genes_tested": round(sum(genes) / len(genes), 1) if genes else None,
            "draws_at_or_below_observed": (
                round(sum(1 for h in hits if h <= observed) / len(hits), 3)
                if hits and observed is not None
                else None
            ),
        }

    return {
        "elements_in_input": total,
        "elements_over_unknown_blocks": in_unknown,
        "genes_losing_every_element": sum(
            1 for g, els in by_gene.items() if els and unknown_count.get(g, 0) == len(els)
        ),
        "full_input": full,
        "reproduces_committed_closure": (
            full["most_active_cell_is_most_expressed"]
            == (c["modes"]["both"]["across_cells"]["most_active_cell_is_most_expressed"])
        ),
        "without_unknown_block_elements": without,
        "random_removal_of_the_same_size": draws(True, without["most_active_cell_is_most_expressed"]),
        "only_unknown_block_elements": only,
        "random_subset_of_the_same_size": draws(False, only["most_active_cell_is_most_expressed"]),
    }


# ---- the syntax reading -------------------------------------------------------------------


def by_case(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Do the elements over blocks constrained on both axes carry larger model effects? The null
    permutes the case labels over blocks, not elements, because elements inside a block share it."""
    sel = [r for r in rows if r["in_unknown"]]
    groups: dict[str, list[dict]] = {}
    for r in sel:
        groups.setdefault(r["case"] or "unmeasured", []).append(r)
    per_case = {
        k: {
            "elements": len(v),
            "blocks": len({tuple(r["block"]) for r in v}),
            "moves_gene": round(_mean([r["moves_gene"] for r in v]), 4),
            "abs_log2_when_moving": round(_mean([r["abs_log2"] for r in v]), 4)
            if any(r["abs_log2"] for r in v)
            else None,
            "cells_acting": round(_mean([r["cells_acting"] for r in v]), 3),
        }
        for k, v in sorted(groups.items())
    }
    # the mammalian axis as a two-group test inside strata: blocks at or above the budget's 5% bar
    for r in sel:
        r["_mammal"] = bool((r["block_mammal_fraction"] or 0) >= 0.05)
    mammal = {m: stratified(sel, m, flag="_mammal") for m in ("moves_gene", "abs_log2", "cells_acting")}
    # elements inside one block share its label, so the element-level p is optimistic: the same test with
    # blocks as the unit, block means compared and the label permuted over blocks
    groups_by_block: dict[tuple, list[dict]] = {}
    for r in sel:
        groups_by_block.setdefault(tuple(r["block"]), []).append(r)
    for m in mammal:
        mammal[m]["block_level"] = block_level(groups_by_block, m)
    return {
        "per_case": per_case,
        "mammal_constrained_block": mammal,
        "cases_present": sorted(groups),
        "syntax_blocks": sum(
            1 for b in {tuple(r["block"]): r["case"] for r in sel}.values() if b == "syntax"
        ),
    }


def block_level(blocks: dict[tuple, list[dict]], metric: str, seed: int = SEED) -> dict[str, Any] | None:
    units = []
    for v in blocks.values():
        m = _mean([r[metric] for r in v])
        if m is not None:
            units.append((bool(v[0]["_mammal"]), m))
    yes = [x for f, x in units if f]
    no = [x for f, x in units if not f]
    if not yes or not no:
        return None
    obs = _mean(yes) - _mean(no)
    labels = [f for f, _ in units]
    vals = [x for _, x in units]
    rng = random.Random(seed)
    null = []
    for _ in range(PERMUTATIONS):
        rng.shuffle(labels)
        a = [x for f, x in zip(labels, vals, strict=True) if f]
        b = [x for f, x in zip(labels, vals, strict=True) if not f]
        null.append(_mean(a) - _mean(b))
    p = (sum(1 for x in null if abs(x) >= abs(obs)) + 1) / (len(null) + 1)
    return {
        "constrained_blocks": len(yes),
        "other_blocks": len(no),
        "difference": round(obs, 4),
        "p_permuted_over_blocks": round(p, 4),
    }


def candidates_on_scored_chromosomes(results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """What the whole-chromosome scoring now says about each syntax candidate that lies on a chromosome
    whose run is complete: a target, a cell and a magnitude, or still nothing."""
    from genomeos.attribution.targets import run_elements

    cand = (load_result("syntax_candidates_genome_wide", results_dir) or {}).get("candidates", [])
    complete = {
        c
        for c in {x["chrom"] for x in cand}
        if (load_result(f"enhancer_targets_all_{c}", results_dir) or {}).get("complete")
    }
    rows = []
    for chrom in sorted(complete):
        elements = run_elements("enhancer_targets_all", chrom, results_dir)
        for c in (x for x in cand if x["chrom"] == chrom):
            inside = [e for e in elements if e["start"] < c["end"] and e["end"] > c["start"]]
            movers = [e for e in inside if e.get("predicted")]
            best = max(movers, key=lambda e: abs(e["predicted"]["log2_fold_change"]), default=None)
            by_cell = (
                (best or {}).get("predicted_coding_by_cell") or (best or {}).get("predicted_by_cell") or {}
            )
            cells = sorted(k for k, v in by_cell.items() if v is not None and abs(v) >= MIN_EFFECT)
            pc = (best or {}).get("predicted_coding") or {}
            rows.append(
                {
                    "block": f"{chrom}:{c['start']}-{c['end']}",
                    "reading_before": c["reading"]["class"],
                    "elements_scored": len(inside),
                    "elements_moving_a_gene": len(movers),
                    "target": pc.get("gene") or (best or {}).get("predicted", {}).get("gene"),
                    "target_is_coding": bool(pc.get("gene")),
                    "log2": (pc or (best or {}).get("predicted") or {}).get("log2_fold_change"),
                    "tissue": (pc or (best or {}).get("predicted") or {}).get("tissue"),
                    "cells_at_bar": cells,
                    "gain": (
                        "target, cell and magnitude"
                        if best and cells
                        else ("target and magnitude, no cell line at the bar" if best else "nothing")
                    ),
                }
            )
    return {
        "candidates": len(cand),
        "on_completely_scored_chromosomes": len(rows),
        "chromosomes": sorted(complete),
        "rows": rows,
    }


# ---- validation against measurement --------------------------------------------------------


def spearman(x: list[float], y: list[float]) -> float | None:
    from genomeos.attribution.closure import spearman as rho

    return rho(x, y)


def rho_p(x: list[float], y: list[float], permutations: int = PERMUTATIONS, seed: int = SEED) -> float | None:
    """Two-sided permutation p of a rank correlation."""
    if len(x) < 3:
        return None
    obs = spearman(x, y)
    if obs is None:
        return None
    rng = random.Random(seed)
    yy = list(y)
    hits = 0
    for _ in range(permutations):
        rng.shuffle(yy)
        r = spearman(x, yy)
        if r is not None and abs(r) >= abs(obs):
            hits += 1
    return round((hits + 1) / (permutations + 1), 4)


def against_mpra(chrom: str, rows: list[dict[str, Any]], elements: list[dict[str, Any]]) -> dict[str, Any]:
    """Predicted deletion magnitude against measured lentiMPRA activity, per cell, for the elements
    both cover; and the same split by whether the element overlaps an unknown block."""
    from genomeos.attribution import mpra

    measured = mpra.load_rows(chrom)
    by_id = {e["id"]: e for e in elements}
    flag = {r["id"]: r for r in rows}
    pairs: list[dict[str, Any]] = []
    for m in measured:
        over = [r for r in rows if r["start"] < m["end"] and r["end"] > m["start"]]
        for r in over:
            e = by_id.get(r["id"]) or {}
            by_cell = e.get("predicted_by_cell") or {}
            pairs.append(
                {
                    "id": r["id"],
                    "in_unknown": flag[r["id"]]["in_unknown"],
                    "activity": m["activity"],
                    "predicted": {c: by_cell.get(c) for c in ("K562", "HepG2")},
                    "abs_log2": r["abs_log2"] or 0.0,
                }
            )
    out: dict[str, Any] = {"elements_covered_by_both": len(pairs)}
    for cell in ("K562", "HepG2"):
        have = [
            p for p in pairs if p["predicted"].get(cell) is not None and p["activity"].get(cell) is not None
        ]
        out[cell] = {
            "n": len(have),
            "rho_abs_effect_vs_activity": (
                spearman([abs(p["predicted"][cell]) for p in have], [p["activity"][cell] for p in have])
                if len(have) >= 3
                else None
            ),
            "rho_signed_effect_vs_activity": (
                spearman([-p["predicted"][cell] for p in have], [p["activity"][cell] for p in have])
                if len(have) >= 3
                else None
            ),
            "p_signed_permuted": rho_p(
                [-p["predicted"][cell] for p in have], [p["activity"][cell] for p in have]
            ),
            "p_abs_permuted": rho_p(
                [abs(p["predicted"][cell]) for p in have], [p["activity"][cell] for p in have]
            ),
            "unknown": {
                "n": sum(1 for p in have if p["in_unknown"]),
                "rho": (
                    spearman(
                        [abs(p["predicted"][cell]) for p in have if p["in_unknown"]],
                        [p["activity"][cell] for p in have if p["in_unknown"]],
                    )
                    if sum(1 for p in have if p["in_unknown"]) >= 3
                    else None
                ),
            },
            "rest": {
                "n": sum(1 for p in have if not p["in_unknown"]),
                "rho": (
                    spearman(
                        [abs(p["predicted"][cell]) for p in have if not p["in_unknown"]],
                        [p["activity"][cell] for p in have if not p["in_unknown"]],
                    )
                    if sum(1 for p in have if not p["in_unknown"]) >= 3
                    else None
                ),
            },
        }
    return out


def against_vista(chrom: str, rows: list[dict[str, Any]]) -> dict[str, Any]:
    p = Path("data/knowledge/vista") / f"rows_{chrom}.json"
    if not p.exists():
        return {"elements": 0}
    vista = json.loads(p.read_text())
    out = {"positive": [], "negative": []}
    for v in vista:
        over = [r for r in rows if r["start"] < v["end"] and r["end"] > v["start"]]
        if v["status"] in out:
            out[v["status"]].append(
                {
                    "id": v["id"],
                    "elements": len(over),
                    "moving": sum(1 for r in over if r["moves_gene"]),
                    "max_abs_log2": max((r["abs_log2"] or 0.0 for r in over), default=None),
                }
            )
    return {
        "elements": len(vista),
        "positives": len(out["positive"]),
        "negatives": len(out["negative"]),
        "positives_with_a_moving_element": sum(1 for x in out["positive"] if x["moving"]),
        "negatives_with_a_moving_element": sum(1 for x in out["negative"] if x["moving"]),
        "mean_max_abs_log2_positive": round(_mean([x["max_abs_log2"] for x in out["positive"]]), 4)
        if out["positive"]
        else None,
        "mean_max_abs_log2_negative": round(_mean([x["max_abs_log2"] for x in out["negative"]]), 4)
        if out["negative"]
        else None,
        "note": "too few VISTA elements here to conclude from; reported, not read",
    }


# ---- coverage and cost --------------------------------------------------------------------


def coverage(chrom: str, rows: list[dict[str, Any]], results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """How much of the unknown space now carries a target, a cell and a magnitude, against what the
    two samples had named before the sweep."""
    from genomeos.attribution.targets import run_elements

    before_ids = {
        e["id"]
        for name in ("constrained_targets", "enhancer_targets")
        for e in run_elements(name, chrom, results_dir)
        if (e.get("predicted_coding") or {}).get("gene")
    }
    blocks = unknown_blocks(chrom, results_dir)
    by_block: dict[tuple, list[dict]] = {}
    for r in rows:
        if r["block"]:
            by_block.setdefault(tuple(r["block"]), []).append(r)
    tiers = Counter(b["tier"] for b in blocks)
    out: dict[str, Any] = {"unknown_blocks": len(blocks), "unknown_bp": sum(b["length"] for b in blocks)}
    per_tier: dict[str, dict[str, Any]] = {}
    for tier in sorted(tiers):
        sel = [b for b in blocks if b["tier"] == tier]
        named_after = [
            b for b in sel if any(r["names_coding"] for r in by_block.get((b["start"], b["end"]), []))
        ]
        named_before = [
            b
            for b in sel
            if any(
                r["names_coding"] and r["id"] in before_ids for r in by_block.get((b["start"], b["end"]), [])
            )
        ]
        with_cell = [
            b
            for b in sel
            if any(r["names_coding"] and r["cells_acting"] for r in by_block.get((b["start"], b["end"]), []))
        ]
        per_tier[tier] = {
            "blocks": len(sel),
            "bp": sum(b["length"] for b in sel),
            "elements_scored": sum(len(by_block.get((b["start"], b["end"]), [])) for b in sel),
            "blocks_with_a_named_target_before": len(named_before),
            "blocks_with_a_named_target_after": len(named_after),
            "blocks_with_a_target_and_a_cell": len(with_cell),
            "bp_with_a_named_target_after": sum(b["length"] for b in named_after),
        }
    out["by_tier"] = per_tier
    out["blocks_with_a_named_target_before"] = sum(
        v["blocks_with_a_named_target_before"] for v in per_tier.values()
    )
    out["blocks_with_a_named_target_after"] = sum(
        v["blocks_with_a_named_target_after"] for v in per_tier.values()
    )
    out["blocks_with_a_target_and_a_cell"] = sum(
        v["blocks_with_a_target_and_a_cell"] for v in per_tier.values()
    )
    out["bp_with_a_named_target_after"] = sum(v["bp_with_a_named_target_after"] for v in per_tier.values())
    named = [r for r in rows if r["in_unknown"] and r["names_coding"]]
    out["element_bp_with_a_named_target_inside_unknown"] = sum(r["unknown_bp"] for r in named)
    out["element_bp_with_a_named_target_inside_unknown_before"] = sum(
        r["unknown_bp"] for r in named if r["id"] in before_ids
    )
    out["elements_over_unknown_blocks"] = sum(1 for r in rows if r["in_unknown"])
    out["elements_over_unknown_naming_a_gene"] = sum(1 for r in rows if r["in_unknown"] and r["names_coding"])
    out["elements_over_unknown_moving_any_gene"] = sum(1 for r in rows if r["in_unknown"] and r["moves_gene"])
    return out


def cost(chrom: str, results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """The requests each scored chromosome took, and what the alternatives would cost per element."""
    per: dict[str, dict[str, Any]] = {}
    for p in sorted(results_dir.glob("enhancer_targets_all_chr*.json")):
        r = json.loads(p.read_text())
        if r.get("chrom") in ("genome_wide", None):
            continue
        per[r["chrom"]] = {
            "elements_total": r.get("elements_total"),
            "scored": r.get("scored"),
            "requests_this_run": r.get("requests_this_run"),
            "complete": bool(r.get("complete")),
        }
    done = [v for v in per.values() if v["complete"] and v["requests_this_run"]]
    ratio = sum(v["requests_this_run"] for v in done) / sum(v["scored"] for v in done) if done else None
    registry = 0
    for c in [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]:
        registry += (load_result(f"enhancer_targets_{c}", results_dir) or {}).get(
            "elements_on_chromosome"
        ) or 0
    chr21 = per.get(chrom) or {}
    in_node = (chr21.get("elements_total") or 0) / max(
        1, (load_result(f"enhancer_targets_{chrom}", results_dir) or {}).get("elements_on_chromosome") or 1
    )
    scored = sum(v["scored"] or 0 for v in per.values())
    remaining = registry * in_node - scored
    return {
        "per_chromosome": per,
        "requests_per_element_on_complete_runs": round(ratio, 3) if ratio else None,
        "registry_elements_genome": registry,
        "share_inside_a_node_on_this_chromosome": round(in_node, 3),
        "elements_scored_so_far": scored,
        "elements_remaining_estimate": round(remaining),
        "requests_remaining_estimate": round(remaining * ratio) if ratio else None,
    }


# ---- the run ------------------------------------------------------------------------------


def run(chrom: str = "chr21", results_dir: Path = RESULTS_DIR, progress=None) -> dict[str, Any]:
    t0 = time.time()

    def say(m: str) -> None:
        if progress:
            progress(m)

    elements = scored_elements(chrom, results_dir)
    say(f"{len(elements)} scored elements")
    rows = add_strata(annotate(chrom, elements, results_dir))
    say(f"{sum(1 for r in rows if r['in_unknown'])} overlap an UNKNOWN block")
    abl = ablation(chrom, results_dir, rows)
    say("ablation done")
    comp = comparison(rows)
    say("comparison done")
    cases = by_case(rows)
    cand = candidates_on_scored_chromosomes(results_dir)
    mpra_out = against_mpra(chrom, rows, elements)
    vista_out = against_vista(chrom, rows)
    cov = coverage(chrom, rows, results_dir)
    return {
        "chrom": chrom,
        "elements_scored": len(elements),
        "elements_over_unknown_blocks": sum(1 for r in rows if r["in_unknown"]),
        "by_tier": Counter(r["tier"] for r in rows if r["in_unknown"]),
        "blocks_touched_by_tier": {
            tier: len({tuple(r["block"]) for r in rows if r["tier"] == tier})
            for tier in sorted({r["tier"] for r in rows if r["in_unknown"]})
        },
        "ablation": abl,
        "comparison": comp,
        "syntax": cases,
        "syntax_candidates": cand,
        "mpra": mpra_out,
        "vista": vista_out,
        "coverage": cov,
        "cost": cost(chrom, results_dir),
        "thresholds": {
            "min_effect_log2": MIN_EFFECT,
            "strata": f"{LENGTH_BINS} length x {GC_BINS} GC x {TSS_BINS} distance-to-nearest-coding-TSS bins",
            "permutations": PERMUTATIONS,
            "random_draws": RANDOM_DRAWS,
        },
        "evidence": EVIDENCE,
        "seconds": round(time.time() - t0, 1),
    }


def run_and_save(chrom: str = "chr21", results_dir: Path = RESULTS_DIR, progress=None) -> dict[str, Any]:
    out = run(chrom, results_dir, progress)
    save_result(f"unknown_scoring_{chrom}", out, results_dir)
    return out


def main() -> None:  # pragma: no cover - the job entry point
    out = run_and_save(progress=lambda m: print("  " + m, flush=True))
    print(json.dumps({k: out[k] for k in ("ablation", "comparison", "syntax", "coverage")}, indent=1)[:8000])
