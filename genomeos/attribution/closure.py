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
from genomeos.results import load_result, save_result

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
    """Elements with a predicted coding target, by target symbol; the constrained run first."""
    seen: set[str] = set()
    by: dict[str, list[dict]] = {}
    for name in ("constrained_targets", "enhancer_targets"):
        for e in (load_result(f"{name}_{chrom}") or {}).get("elements", []):
            pc = e.get("predicted_coding") or {}
            if not pc.get("gene") or e["id"] in seen:
                continue
            seen.add(e["id"])
            sign = 1.0 if pc.get("action") == "activates" else -1.0
            by.setdefault(pc["gene"], []).append(
                {
                    "id": e["id"],
                    "start": e["start"],
                    "end": e["end"],
                    "effect": sign * abs(pc["log2_fold_change"]),
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
            active = [e for e in by_gene.get(g.symbol, []) if m["index"].overlapping(e["start"], e["end"])]
            row["cells"][cell] = {
                "expression": m["expression"][g.symbol],
                "expressed": m["expression"][g.symbol] >= EXPRESSED,
                "promoter_open": m["promoter_open"][g.symbol],
                "active_elements": len(active),
                "input": round(sum(e["effect"] for e in active), 3),
            }
        rows.append(row)
    return {
        "chrom": chrom,
        "cells": used,
        "strand_orientation": {c: measured[c]["strand_orientation"] for c in used},
        "coding_genes": len(genes),
        "genes_with_elements": sum(1 for r in rows if r["elements"]),
        "signal_threshold": signal,
        "expressed_threshold": EXPRESSED,
        **judge(rows, used, seed),
        "genes": rows,
        "cost": {
            "seconds": round(time.time() - t0, 1),
            "mb_fetched": round(sum(m["mb"] for m in measured.values()), 1),
        },
        "evidence": EVIDENCE,
    }


def judge(rows: list[dict], cells: list[str], seed: int = 0) -> dict:
    """The closure tests over the assembled gene by cell table."""
    within: dict[str, dict] = {}
    for cell in cells:
        c = [r["cells"][cell] for r in rows]

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
    # across cells: genes with elements whose input and expression both vary
    rng = random.Random(seed)
    rhos, prom_rhos, argmax_hits, argmax_n = [], [], 0, 0
    shuffled_hits = 0
    for r in rows:
        if not r["elements"] or len(cells) < 3:
            continue
        inp = [r["cells"][c]["input"] for c in cells]
        exp = [r["cells"][c]["expression"] for c in cells]
        prom = [1.0 if r["cells"][c]["promoter_open"] else 0.0 for c in cells]
        if len(set(inp)) > 1 and len(set(exp)) > 1:
            rho = spearman(inp, exp)
            if rho is not None:
                rhos.append(rho)
            prho = spearman(prom, exp)
            if prho is not None:
                prom_rhos.append(prho)
            best_in = max(range(len(cells)), key=lambda i: inp[i])
            best_ex = max(range(len(cells)), key=lambda i: exp[i])
            argmax_n += 1
            argmax_hits += int(best_in == best_ex)
            sh = list(range(len(cells)))
            rng.shuffle(sh)
            shuffled_hits += int(sh[best_in] == best_ex)
    rejected = []
    unexplained = []
    for r in rows:
        for c in cells:
            x = r["cells"][c]
            if x["promoter_open"] and x["input"] > 0 and not x["expressed"]:
                rejected.append(
                    {"gene": r["gene"], "cell": c, "input": x["input"], "expression": x["expression"]}
                )
            if x["expressed"] and not x["promoter_open"] and x["active_elements"] == 0:
                unexplained.append({"gene": r["gene"], "cell": c, "expression": x["expression"]})
    rejected.sort(key=lambda d: -d["input"])
    unexplained.sort(key=lambda d: -d["expression"])
    return {
        "within_cell": within,
        "across_cells": {
            "genes_tested": argmax_n,
            "mean_rho_input_vs_expression": round(sum(rhos) / len(rhos), 3) if rhos else None,
            "mean_rho_promoter_vs_expression": round(sum(prom_rhos) / len(prom_rhos), 3)
            if prom_rhos
            else None,
            "most_active_cell_is_most_expressed": round(argmax_hits / argmax_n, 3) if argmax_n else None,
            "same_under_shuffled_cells": round(shuffled_hits / argmax_n, 3) if argmax_n else None,
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
