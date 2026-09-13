# SPDX-License-Identifier: AGPL-3.0-or-later
"""The gene-level closure test: do the attributions reproduce what the cell makes?

The house argument, applied where it has enough bandwidth to falsify. For every coding
gene of a chromosome and every cell type that has both a DNase reader and a measured
RNA-seq track, three things are read independently:

- the promoter's openness (a DNase peak within 1 kb of the TSS), the reader's claim;
- the regulatory input from the elements attributed to the gene by the AlphaGenome
  deletions: an element is active in the cell when it overlaps a DNase peak there, and
  the input is the signed sum of the predicted magnitudes (activators positive,
  repressors negative) over the active elements;
- the measured expression, the fraction of the canonical exons' bases carrying
  ENCODE total RNA-seq signal on the gene's strand.

Then the closure: within a cell, do open promoters and positive input go with
expression; across cells, is the cell where a gene's elements are most active the cell
where the gene is most expressed, against a control that shuffles the cells. Genes that
are expressed with a closed promoter and no active element, and genes whose promoter
and elements are active but which are silent, are listed by name: those are the
attributions the cell rejects, and they are the point of the exercise.
"""

from __future__ import annotations

import random
import time
from typing import Any

from genomeos.coords import Strand
from genomeos.genome.annotation import Annotation, default_gencode
from genomeos.genome.reader import PROMOTER_WINDOW, PeakIndex, load_peaks
from genomeos.genome.rna_measured import MeasuredRna
from genomeos.results import save_result

CELLS = ("K562", "HepG2", "GM12878", "IMR-90")  # a reader and a total RNA-seq track exist for each
EXPRESSED = 0.3  # mean covered fraction of the canonical exons from which a gene counts as made
EVIDENCE = {
    "expression": "experimental: ENCODE total RNA-seq signal over the canonical exons, per strand",
    "promoter": "experimental: ENCODE DNase-seq peak within 1 kb of the TSS (the reader)",
    "elements": "predicted: AlphaGenome deletion target and magnitude, active where a DNase peak overlaps",
}


def canonical_exons(g) -> list[tuple[int, int]]:
    """The Ensembl-canonical transcript's exons, else the longest transcript's."""
    ts = list(g.transcripts.values())
    if not ts:
        return []
    canon = [t for t in ts if "Ensembl_canonical" in t.tags] or sorted(
        ts, key=lambda t: -(t.locus.end - t.locus.start)
    )
    return sorted((e.start, e.end) for e in canon[0].exons)


def merge(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def attributed_elements(chrom: str) -> dict[str, list[dict]]:
    """Elements with a predicted coding target, by target symbol; the whole-chromosome run first."""
    from genomeos.attribution.targets import attributed

    by: dict[str, list[dict]] = {}
    for e in attributed(chrom):
        pc = e["predicted_coding"]
        sign = 1.0 if pc.get("action") == "activates" else -1.0
        # per-cell: the deletion scored on the cell line's own track; a fall in expression on
        # deletion (negative log2) means the element activates the gene in that cell
        by_cell = {
            c: -float(v) for c, v in (e.get("predicted_coding_by_cell") or {}).items() if v is not None
        }
        by.setdefault(pc["gene"], []).append(
            {
                "id": e["id"],
                "start": e["start"],
                "end": e["end"],
                "effect": sign * abs(pc["log2_fold_change"]),
                "by_cell": by_cell,
            }
        )
    return by


def spearman(x: list[float], y: list[float]) -> float | None:
    """Rank correlation with average ranks for ties; None when either side is constant."""
    n = len(x)
    if n < 3:
        return None

    def ranks(v: list[float]) -> list[float]:
        order = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]:
                j += 1
            for k in range(i, j + 1):
                r[order[k]] = (i + j) / 2 + 1
            i = j + 1
        return r

    rx, ry = ranks(x), ranks(y)
    mx, my = sum(rx) / n, sum(ry) / n
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    if sxx == 0 or syy == 0:
        return None
    return sum((a - mx) * (b - my) for a, b in zip(rx, ry, strict=True)) / (sxx * syy) ** 0.5


def _covered(rna: MeasuredRna, merged: dict, strand: str, a: int, b: int) -> float:
    # the merged interval that contains the exon carries its covered fraction
    for s, e in merged[strand]:
        if s <= a and b <= e:
            return rna.covered_fraction(strand, s, e)
    return 0.0


def measure(chrom: str, genes: list, cells: tuple[str, ...], signal: float, progress=None) -> dict:
    """Per cell: expression per gene (mean covered fraction of canonical exons), promoter openness,
    and the peak index for element activity."""
    exons_by_strand: dict[str, list[tuple[int, int]]] = {"+": [], "-": []}
    for g in genes:
        exons_by_strand[g.locus.strand.value].extend(canonical_exons(g))
    merged = {s: merge(v) for s, v in exons_by_strand.items()}
    out: dict[str, dict[str, Any]] = {}
    for cell in cells:
        peaks = load_peaks(cell, chrom)
        if not peaks:
            continue
        idx = PeakIndex(peaks)
        rna, orientation = oriented(cell, chrom, signal, merged)
        if rna is None:
            continue
        rna.prepare(merged, progress)
        expr: dict[str, float] = {}
        prom: dict[str, bool] = {}
        for g in genes:
            ex = canonical_exons(g)
            fr = [_covered(rna, merged, g.locus.strand.value, a, b) for a, b in ex]
            expr[g.symbol] = round(sum(fr) / len(fr), 4) if fr else 0.0
            tss = g.locus.end - 1 if g.locus.strand is Strand.MINUS else g.locus.start
            prom[g.symbol] = bool(idx.overlapping(tss - PROMOTER_WINDOW, tss + PROMOTER_WINDOW))
        out[cell] = {
            "expression": expr,
            "promoter_open": prom,
            "index": idx,
            "mb": rna.bytes_fetched / 1e6,
            "strand_orientation": orientation,
        }
    return out


def oriented(
    cell: str, chrom: str, signal: float, merged: dict[str, list[tuple[int, int]]], probe: int = 40
) -> tuple[MeasuredRna | None, str]:
    """The cell's RNA-seq tracks in the orientation that puts signal on the genes' own strands.

    ENCODE labels a track "plus strand signal" by the read, and for some library protocols the read
    is antisense to the transcript: IMR-90's total RNA-seq reads APP (minus strand) on the plus track.
    A probe over a few dozen exons per strand decides; the choice is recorded with the result.
    """
    base = MeasuredRna(cell, chrom, signal=signal)
    if base.missing or not base.tracks.get("tracks"):
        return None, "no track"
    tracks = base.tracks
    swapped = {**tracks, "tracks": {"+": tracks["tracks"].get("-"), "-": tracks["tracks"].get("+")}}
    swapped["tracks"] = {k: v for k, v in swapped["tracks"].items() if v}
    sample = {s: v[:probe] for s, v in merged.items()}
    scores = {}
    for name, tr in (("as labelled", tracks), ("swapped", swapped)):
        m = MeasuredRna(cell, chrom, signal=signal, tracks=tr)
        m.prepare(sample)
        scores[name] = sum(m.covered_fraction(s, a, b) for s, v in sample.items() for a, b in v)
    if scores["swapped"] > 2 * scores["as labelled"]:
        return MeasuredRna(cell, chrom, signal=signal, tracks=swapped), "swapped"
    return base, "as labelled"


def closure(
    chrom: str, cells: tuple[str, ...] = CELLS, signal: float = 0.05, progress=None, seed: int = 0
) -> dict:
    t0 = time.time()
    gff = default_gencode({chrom})
    if gff is None:
        raise FileNotFoundError(f"no GENCODE models for {chrom}; run genomeos data fetch --chrom {chrom}")
    ann = Annotation.from_gff3(gff, {chrom})
    genes = [g for g in ann.protein_coding() if g.locus.chrom == chrom]
    by_gene = attributed_elements(chrom)
    measured = measure(chrom, genes, cells, signal, progress)
    used = list(measured)
    rows: list[dict] = []
    for g in genes:
        row: dict[str, Any] = {"gene": g.symbol, "elements": len(by_gene.get(g.symbol, [])), "cells": {}}
        for cell in used:
            m = measured[cell]
            els = by_gene.get(g.symbol, [])
            active = [e for e in els if m["index"].overlapping(e["start"], e["end"])]
            scored = [e for e in els if cell in e["by_cell"]]
            row["cells"][cell] = {
                "expression": m["expression"][g.symbol],
                "expressed": m["expression"][g.symbol] >= EXPRESSED,
                "promoter_open": m["promoter_open"][g.symbol],
                "active_elements": len(active),
                # dnase: tissue-agnostic magnitude, gated by a DNase peak on the element
                "input": round(sum(e["effect"] for e in active), 3),
                # cell: the deletion scored on this cell's own track, no gating
                "input_cell": round(sum(e["by_cell"][cell] for e in scored), 3) if scored else None,
                # both: the cell's own magnitude, only where a DNase peak also sits on the element
                "input_both": (
                    round(sum(e["by_cell"][cell] for e in scored if e in active), 3) if scored else None
                ),
                "scored_elements": len(scored),
            }
        rows.append(row)
    has_cell = any(x["scored_elements"] for r in rows for x in r["cells"].values())
    modes = {"dnase": judge(rows, used, seed, "input")}
    if has_cell:
        modes["cell"] = judge(rows, used, seed, "input_cell")
        modes["both"] = judge(rows, used, seed, "input_both")
    return {
        "chrom": chrom,
        "cells": used,
        "strand_orientation": {c: measured[c]["strand_orientation"] for c in used},
        "coding_genes": len(genes),
        "genes_with_elements": sum(1 for r in rows if r["elements"]),
        "signal_threshold": signal,
        "expressed_threshold": EXPRESSED,
        **modes["dnase"],
        "modes": modes,
        "per_cell_scores": has_cell,
        "genes": rows,
        "cost": {
            "seconds": round(time.time() - t0, 1),
            "mb_fetched": round(sum(m["mb"] for m in measured.values()), 1),
        },
        "evidence": EVIDENCE,
    }


def tie_fair_hit(inp: list[float], exp: list[float]) -> float:
    """Probability that a random cell among the most active and a random cell among the most expressed
    are the same cell."""
    a = {i for i, v in enumerate(inp) if v == max(inp)}
    b = {i for i, v in enumerate(exp) if v == max(exp)}
    return len(a & b) / (len(a) * len(b))


def judge(
    rows: list[dict], cells: list[str], seed: int = 0, key: str = "input", permutations: int = 1000
) -> dict:
    """The closure tests over the assembled gene by cell table, with `key` as the element input:
    `input` (DNase-gated, tissue-agnostic magnitude), `input_cell` (the cell's own deletion score) or
    `input_both`. A None input means the model has no score for that cell; the pair is left out."""
    within: dict[str, dict] = {}
    for cell in cells:
        c = [dict(r["cells"][cell], input=r["cells"][cell].get(key)) for r in rows]
        c = [x for x in c if x["input"] is not None]

        def frac(sel):
            sel = list(sel)
            return (
                (round(sum(1 for x in sel if x["expressed"]) / len(sel), 3), len(sel)) if sel else (None, 0)
            )

        within[cell] = {
            "expressed_promoter_open": frac(x for x in c if x["promoter_open"]),
            "expressed_promoter_closed": frac(x for x in c if not x["promoter_open"]),
            "expressed_open_with_activating_input": frac(
                x for x in c if x["promoter_open"] and x["input"] > 0
            ),
            "expressed_open_without_input": frac(
                x for x in c if x["promoter_open"] and x["active_elements"] == 0
            ),
            "expressed_open_with_repressing_input": frac(
                x for x in c if x["promoter_open"] and x["input"] < 0
            ),
        }
    # across cells: genes with elements whose input and expression both vary. Ties are broken fairly:
    # a gene scores the probability that a random pick among its most-active cells and a random pick
    # among its most-expressed cells land on the same cell, so no cell order can inflate the rate.
    # The null permutes each gene's inputs across cells, the same gene set and the same tie rule.
    rng = random.Random(seed)
    rhos, prom_rhos, pairs = [], [], []
    for r in rows:
        if not r["elements"] or len(cells) < 3:
            continue
        if any(r["cells"][c].get(key) is None for c in cells):
            continue
        inp = [r["cells"][c][key] for c in cells]
        exp = [r["cells"][c]["expression"] for c in cells]
        prom = [1.0 if r["cells"][c]["promoter_open"] else 0.0 for c in cells]
        if len(set(inp)) > 1 and len(set(exp)) > 1:
            rho = spearman(inp, exp)
            if rho is not None:
                rhos.append(rho)
            prho = spearman(prom, exp)
            if prho is not None:
                prom_rhos.append(prho)
            pairs.append((inp, exp))
    argmax_n = len(pairs)
    observed = sum(tie_fair_hit(i, e) for i, e in pairs) / argmax_n if argmax_n else None
    null = []
    for _ in range(permutations if argmax_n else 0):
        h = 0.0
        for i, e in pairs:
            y = i[:]
            rng.shuffle(y)
            h += tie_fair_hit(y, e)
        null.append(h / argmax_n)
    p_value = (sum(1 for h in null if h >= observed) + 1) / (len(null) + 1) if null else None
    rejected = []
    unexplained = []
    for r in rows:
        for c in cells:
            x = r["cells"][c]
            inp = x.get(key)
            if inp is not None and x["promoter_open"] and inp > 0 and not x["expressed"]:
                rejected.append({"gene": r["gene"], "cell": c, "input": inp, "expression": x["expression"]})
            if x["expressed"] and not x["promoter_open"] and x["active_elements"] == 0:
                unexplained.append({"gene": r["gene"], "cell": c, "expression": x["expression"]})
    rejected.sort(key=lambda d: -d["input"])
    unexplained.sort(key=lambda d: -d["expression"])
    return {
        "input": key,
        "within_cell": within,
        "across_cells": {
            "genes_tested": argmax_n,
            "mean_rho_input_vs_expression": round(sum(rhos) / len(rhos), 3) if rhos else None,
            "mean_rho_promoter_vs_expression": round(sum(prom_rhos) / len(prom_rhos), 3)
            if prom_rhos
            else None,
            "most_active_cell_is_most_expressed": round(observed, 3) if observed is not None else None,
            "same_under_shuffled_cells": round(sum(null) / len(null), 3) if null else None,
            "p_value": round(p_value, 4) if p_value is not None else None,
            "permutations": len(null),
            "chance": round(1 / len(cells), 3) if cells else None,
        },
        "rejected_attributions": rejected[:40],
        "rejected_count": len(rejected),
        "unexplained_expression": unexplained[:40],
        "unexplained_count": len(unexplained),
    }


def run_and_save(chrom: str, **kw) -> dict:
    out = closure(chrom, **kw)
    save_result(f"closure_{chrom}", {k: v for k, v in out.items()})
    return out
