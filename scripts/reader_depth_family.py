"""Score the reader's per-biosample family against assay depth, on the registration
committed in docs/NODES-READER-WRITER.md ("Is the rest of the reader's family assay
depth too?", 41bc7d6).

Reads only what is already on disk: data/results/reader_genome_wide.json for the
readings, data/results/epigenome_chr*.json for the mark peak counts, and
data/results/dnase_*_chr*.bed.gz for DNase peak width. No network, no model request.

Writes data/results/reader_depth_family.json.
"""

from __future__ import annotations

import glob
import gzip
import json
import os
import random
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

RESULTS = os.path.join(ROOT, "data", "results")
ODD = [f"chr{i}" for i in range(1, 23, 2)]
EVEN = [f"chr{i}" for i in range(2, 23, 2)]
PERMUTATIONS = 10_000
SEED = 20260922


def _rank(v: list[float]) -> list[float]:
    order = sorted(range(len(v)), key=lambda i: v[i])
    out = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            out[order[k]] = avg
        i = j + 1
    return out


def _pearson(x: list[float], y: list[float]) -> float:
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    num = sum((a - mx) * (b - my) for a, b in zip(x, y, strict=True))
    den = (sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y)) ** 0.5
    return round(num / den, 4) if den else 0.0


def spearman(x: list[float], y: list[float]) -> float:
    return _pearson(_rank(x), _rank(y))


def fit_residuals(x: list[float], y: list[float]) -> tuple[float, float, list[float]]:
    """Least-squares y ~ a + b*x; returns (a, b, residuals)."""
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y, strict=True))
    b = sxy / sxx if sxx else 0.0
    a = my - b * mx
    return a, b, [yy - (a + b * xx) for xx, yy in zip(x, y, strict=True)]


def permutation_p(x: list[float], y: list[float], observed: float, n: int, seed: int) -> float:
    """One-sided p for Spearman rho >= observed, by shuffling y's labels."""
    rng = random.Random(seed)
    shuffled = list(y)
    hits = 0
    for _ in range(n):
        rng.shuffle(shuffled)
        if spearman(x, shuffled) >= observed:
            hits += 1
    return round((hits + 1) / (n + 1), 5)


def load_readings() -> tuple[list[str], dict, dict]:
    with open(os.path.join(RESULTS, "reader_genome_wide.json")) as f:
        d = json.load(f)
    chroms = d["chromosomes"]
    cells = [c for c in d["totals"] if c != "seconds"]
    return cells, chroms, d["totals"]


def summed(chroms: dict, cell: str, field: str, only: list[str] | None = None) -> int:
    keys = only if only is not None else list(chroms)
    return sum((chroms[c][cell].get(field) or 0) for c in keys if cell in chroms.get(c, {}))


def mark_peaks() -> dict[str, dict[str, int]]:
    """Genome-wide peak count per (biosample, mark), summed over the per-chromosome
    epigenome summaries already on disk."""
    tot: dict[str, dict[str, int]] = {}
    for path in glob.glob(os.path.join(RESULTS, "epigenome_chr*.json")):
        if "genome_wide" in path:
            continue
        with open(path) as f:
            d = json.load(f)
        for cell, row in d.get("cell_types", {}).items():
            for mark, p in (row.get("peaks") or {}).items():
                tot.setdefault(cell, {}).setdefault(mark, 0)
                tot[cell][mark] += p["n"]
    return tot


def dnase_peak_width(cells: list[str]) -> dict[str, dict[str, float]]:
    """Mean DNase peak width and total open bp per biosample, over the autosomes, read
    from the peak files on disk. The file stem is the biosample with non-word characters
    replaced, as `slug` writes it."""
    out: dict[str, dict[str, float]] = {}
    for cell in cells:
        slug = re.sub(r"[^A-Za-z0-9]+", "_", cell).strip("_")
        n, bp = 0, 0
        for chrom in ODD + EVEN:
            path = os.path.join(RESULTS, f"dnase_{slug}_{chrom}.bed.gz")
            if not os.path.exists(path):
                continue
            with gzip.open(path, "rt") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.split("\t")
                    if len(parts) < 2:
                        continue
                    n += 1
                    bp += int(parts[1]) - int(parts[0])
        out[cell] = {
            "peaks": n,
            "open_bp": bp,
            "mean_width": round(bp / n, 2) if n else 0.0,
        }
    return out


def main() -> None:
    cells, chroms, totals = load_readings()
    depth = {c: summed(chroms, c, "peaks") for c in cells}

    # 1. the family, banded against the right covariate
    marks = mark_peaks()
    readings = {
        "enhancers_active": ("DNase peaks", [summed(chroms, c, "enhancers_active") for c in cells]),
        "genes_read_by_marks": ("DNase peaks", [totals[c]["genes_read_by_marks"] for c in cells]),
        "nodes_silent": ("DNase peaks", [totals[c]["nodes_silent"] for c in cells]),
        "genes_read_open": ("DNase peaks", [totals[c]["genes_read_open"] for c in cells]),
        "genes_read": ("DNase peaks", [totals[c]["genes_read"] for c in cells]),
        "nodes_open": ("DNase peaks", [summed(chroms, c, "nodes_open") for c in cells]),
    }
    family = {}
    for name, (cov, vals) in readings.items():
        family[name] = {
            "covariate": cov,
            "min": min(vals),
            "max": max(vals),
            "span": round(max(vals) / max(1, min(vals)), 2),
            "rho": spearman([depth[c] for c in cells], vals),
        }
    # genes_poised gets the mark that calls it, not DNase
    poised = [totals[c]["genes_poised"] for c in cells]
    k27me3 = [marks.get(c, {}).get("H3K27me3", 0) for c in cells]
    k27ac = [marks.get(c, {}).get("H3K27ac", 0) for c in cells]
    family["genes_poised"] = {
        "covariate": "H3K27me3 peaks",
        "min": min(poised),
        "max": max(poised),
        "span": round(max(poised) / max(1, min(poised)), 2),
        "rho": spearman(k27me3, poised),
        "rho_vs_dnase": spearman([depth[c] for c in cells], poised),
        "rho_vs_h3k27ac": spearman(k27ac, poised),
    }

    # 2. kill test 1: does the depth residual of enhancers_active replicate across
    #    disjoint halves of the same genome?
    halves = {}
    resid = {}
    for label, keys in (("odd", ODD), ("even", EVEN)):
        x = [summed(chroms, c, "peaks", keys) for c in cells]
        y = [summed(chroms, c, "enhancers_active", keys) for c in cells]
        a, b, r = fit_residuals([float(v) for v in x], [float(v) for v in y])
        halves[label] = {"intercept": round(a, 1), "slope": round(b, 4), "chroms": len(keys)}
        resid[label] = r
    rho_halves = spearman(resid["odd"], resid["even"])
    p_halves = permutation_p(resid["odd"], resid["even"], rho_halves, PERMUTATIONS, SEED)
    replicates = rho_halves >= 0.5 and p_halves < 0.05

    # 3. kill test 2: is whatever survives depth just peak-calling shape?
    width = dnase_peak_width(cells)
    x_all = [float(depth[c]) for c in cells]
    y_all = [float(summed(chroms, c, "enhancers_active")) for c in cells]
    a_all, b_all, r_all = fit_residuals(x_all, y_all)
    rho_width = spearman([width[c]["mean_width"] for c in cells], r_all)
    is_shape = abs(rho_width) >= 0.7

    # 4. kill test 3: does the rate over-correct?
    rate = [summed(chroms, c, "enhancers_active") / max(1, depth[c]) * 100_000 for c in cells]
    rho_rate = spearman(x_all, rate)
    over_corrects = rho_rate <= -0.7
    removes_depth = abs(rho_rate) < 0.4

    # 5. does the rate move a published conclusion?
    by_count = sorted(cells, key=lambda c: -summed(chroms, c, "enhancers_active"))
    by_rate = sorted(cells, key=lambda c: -rate[cells.index(c)])
    k562 = cells.index("K562")
    hepg2 = cells.index("HepG2")
    reversal = {
        "K562_count": summed(chroms, "K562", "enhancers_active"),
        "HepG2_count": summed(chroms, "HepG2", "enhancers_active"),
        "K562_rate": round(rate[k562], 1),
        "HepG2_rate": round(rate[hepg2], 1),
        "count_says": "K562 > HepG2"
        if rate and by_count.index("K562") < by_count.index("HepG2")
        else "HepG2 > K562",
        "rate_says": "K562 > HepG2" if rate[k562] > rate[hepg2] else "HepG2 > K562",
        "reverses": (rate[k562] > rate[hepg2]) != (y_all[k562] > y_all[hepg2]),
    }

    verdict = (
        "relabel only: the depth residual does not replicate across disjoint halves of the genome"
        if not replicates
        else (
            "relabel only: the residual replicates but tracks mean DNase peak width, so it is a "
            "second assay property"
            if is_shape
            else "ship the rate beside the count"
        )
    )

    out = {
        "result": "reader_depth_family",
        "registration": "docs/NODES-READER-WRITER.md, committed 41bc7d6 before this ran",
        "biosamples": cells,
        "dnase_peaks": depth,
        "family": family,
        "residual_replication": {
            "halves": halves,
            "rho": rho_halves,
            "p_permutation": p_halves,
            "permutations": PERMUTATIONS,
            "threshold": "rho >= 0.5 and p < 0.05",
            "replicates": replicates,
        },
        "residual_vs_peak_width": {
            "full_fit": {"intercept": round(a_all, 1), "slope": round(b_all, 4)},
            "mean_width_bp": {c: width[c]["mean_width"] for c in cells},
            "rho": rho_width,
            "threshold": "|rho| >= 0.7 means peak-calling shape, still the assay",
            "is_shape": is_shape,
        },
        "rate": {
            "enhancers_active_per_100k_peaks": {c: round(rate[cells.index(c)], 1) for c in cells},
            "rho_vs_dnase_peaks": rho_rate,
            "over_corrects": over_corrects,
            "removes_depth": removes_depth,
            "order_by_count": by_count,
            "order_by_rate": by_rate,
            "published_comparison": reversal,
        },
        "verdict": verdict,
        "evidence": (
            "derived: ENCODE DNase-seq and Histone ChIP-seq peak counts already on disk; no request made"
        ),
    }
    path = os.path.join(RESULTS, "reader_depth_family.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: out[k] for k in ("family", "residual_replication", "verdict")}, indent=1))
    print(json.dumps(out["residual_vs_peak_width"], indent=1))
    print(json.dumps(out["rate"], indent=1))
    print(f"written: {path}", flush=True)


if __name__ == "__main__":
    main()
