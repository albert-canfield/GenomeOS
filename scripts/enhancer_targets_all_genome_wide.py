# SPDX-License-Identifier: AGPL-3.0-or-later
"""One table over every chromosome's all-elements summary, as the chain lands them.

    uv run python scripts/enhancer_targets_all_genome_wide.py

Reads data/results/enhancer_targets_all_chr*.json (summaries only; the element tables stay local) and
writes enhancer_targets_all_genome_wide.json: per chromosome the counts and rates, the totals over the
chromosomes complete so far with the requests spent, the spread of each rate across chromosomes rather
than one pooled number, the controls each rate has to be read against, and whether the rates drift with
chromosome size (chr21, chr22 and chrY are acrocentric or nearly gene-free, and a reading taken on them
alone has already misled one lane).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from statistics import median

from genomeos.molecules.proteome import CHROM_LENGTHS
from genomeos.results import save_result

# Chromosomes whose architecture is not the genome's: two acrocentric arms and the male-specific
# chromosome, all three small and gene-poor. Rates are reported with and without them.
UNUSUAL = ("chr21", "chr22", "chrY")

# What each rate has to be read against. Both come from measurements this project already made, and
# both qualify a rate that reads as a result on its own.
CONTROLS = {
    "fraction_with_target": {
        "control": 0.87,
        "what": "a derived layer names some target at 52 of 60 matched random windows (87%), "
        "against 11 of 12 published loci",
        "source": "docs/LOCI-BENCHMARK.md, data/results/loci_benchmark.json (commits 303d0b0, ef80075)",
        "reading": "a naming rate is not evidence by itself; what separates a real element from a "
        "matched window is a direction and a variable position on constrained sequence",
    },
    "coding_target_inside_domain": {
        "control": 0.791,
        "measured": 0.817,
        "excess": 0.026,
        "what": "over the whole deletion archive of 113,399 elements the same caller keeps a coding "
        "target inside the element's own CTCF node for 81.7% of elements, against 79.1% for the same "
        "number of boundaries placed at random: 2.6 points, not nine times out of ten",
        "source": "data/results/domains_oriented_comparison.json (commit dbad593)",
        "reading": "the node's advantage over a random partition at the same resolution is small, and "
        "it is measured on a model's reading rather than on observed enhancer-gene pairs",
    },
    "coding_target_agrees_with_nearest": {
        "control": None,
        "what": "no matched-window control has been measured for nearest-TSS agreement; the nearest "
        "TSS heuristic names the published target at 8 of 12 loci, the same count the derived layers "
        "reach on a different 8",
        "source": "docs/LOCI-BENCHMARK.md (commit 303d0b0)",
        "reading": "quote it as agreement between two readings of the same element, not as accuracy",
    },
}


def _spearman(xs: list[float], ys: list[float]) -> float | None:
    """Rank correlation, no ties expected among chromosome lengths or rates."""
    n = len(xs)
    if n < 4:
        return None

    def ranks(vs: list[float]) -> list[float]:
        order = sorted(range(n), key=lambda i: vs[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and vs[order[j + 1]] == vs[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    rx, ry = ranks(xs), ranks(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry, strict=True))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return round(num / den, 3) if den else None


def _rates(counts: dict) -> dict:
    """The three headline rates from one chromosome's or one pool's counts."""
    agree = counts["verdicts"].get("agrees with nearest TSS in domain", 0)
    inside = agree + counts["verdicts"].get("another gene in the same domain", 0)
    coding = counts["coding_named"]
    return {
        "fraction_with_target": round(counts["named"] / counts["scored"], 3) if counts["scored"] else None,
        "coding_target_agrees_with_nearest": round(agree / coding, 3) if coding else None,
        "coding_target_inside_domain": round(inside / coding, 3) if coding else None,
    }


def _pool(rows: list[dict]) -> dict:
    tot = {"elements": 0, "scored": 0, "named": 0, "coding_named": 0, "strong": 0, "silencer_like": 0}
    verdicts: dict[str, int] = {}
    for r in rows:
        for k in ("elements", "scored", "named", "coding_named", "strong", "silencer_like"):
            tot[k] += r["counts"][k]
        for k, n in r["counts"]["verdicts"].items():
            verdicts[k] = verdicts.get(k, 0) + n
    tot["verdicts"] = verdicts
    return tot


def _spread(rows: list[dict], key: str) -> dict | None:
    vals = [(r["rates"][key], r["chrom"]) for r in rows if r["rates"][key] is not None]
    if not vals:
        return None
    vals.sort()
    return {
        "n": len(vals),
        "min": vals[0][0],
        "lowest": vals[0][1],
        "median": round(median(v for v, _ in vals), 3),
        "max": vals[-1][0],
        "highest": vals[-1][1],
        "range": round(vals[-1][0] - vals[0][0], 3),
        "per_chromosome": {c: v for v, c in vals},
    }


def aggregate(results_dir: Path = Path("data/results")) -> dict:
    rows: dict[str, dict] = {}
    complete: list[dict] = []
    requests = 0
    for f in sorted(results_dir.glob("enhancer_targets_all_chr*.json")):
        d = json.loads(f.read_text())
        s = d.get("summary") or {}
        vc = s.get("verdicts_coding") or {}
        coding_named = sum(n for k, n in vc.items() if k != "no predicted effect")
        counts = {
            "elements": d.get("elements_total") or 0,
            "scored": d.get("scored") or 0,
            "named": s.get("with_predicted_target") or 0,
            "coding_named": coding_named,
            "strong": s.get("strong") or 0,
            "silencer_like": s.get("silencer_like") or 0,
            "verdicts": dict(vc),
        }
        rates = _rates(counts)
        rows[d["chrom"]] = {
            "elements": d.get("elements_total"),
            "scored": d.get("scored", 0),
            "complete": bool(d.get("complete")),
            "length": CHROM_LENGTHS.get(d["chrom"]),
            "fraction_with_target": s.get("fraction_with_target"),
            "coding_target_agrees_with_nearest": s.get("coding_target_agrees_with_nearest"),
            "coding_target_inside_domain": rates["coding_target_inside_domain"],
            "fraction_inside_domain": s.get("fraction_inside_domain"),
            "strong": s.get("strong"),
            "silencer_like": s.get("silencer_like"),
            "requests_this_run": d.get("requests_this_run"),
        }
        requests += d.get("requests_this_run") or 0
        if d.get("complete"):
            complete.append({"chrom": d["chrom"], "counts": counts, "rates": rates})

    keys = ("fraction_with_target", "coding_target_agrees_with_nearest", "coding_target_inside_domain")
    pooled = _pool(complete)
    genome = {
        **{k: v for k, v in pooled.items() if k != "verdicts"},
        **_rates(pooled),
        "verdicts_coding": dict(sorted(pooled["verdicts"].items(), key=lambda kv: -kv[1])),
    }
    usual = [r for r in complete if r["chrom"] not in UNUSUAL]
    unusual = [r for r in complete if r["chrom"] in UNUSUAL]

    sized = [r for r in complete if r["chrom"] in CHROM_LENGTHS]
    lengths = [CHROM_LENGTHS[r["chrom"]] for r in sized]
    density = [r["counts"]["elements"] / CHROM_LENGTHS[r["chrom"]] * 1e6 for r in sized]
    trend = {
        "axis": "hg38 chromosome length against the chromosome's own rate, over the complete chromosomes",
        "chromosomes": len(lengths),
        **{k: _spearman(lengths, [r["rates"][k] for r in sized]) for k in keys},
        "against_element_density": {
            "axis": "registry elements per Mb, the stand-in for gene density this run already holds",
            **{k: _spearman(density, [r["rates"][k] for r in sized]) for k in keys},
            "per_chromosome": dict(
                sorted(
                    ((r["chrom"], round(d, 1)) for r, d in zip(sized, density, strict=True)),
                    key=lambda p: -p[1],
                )
            ),
        },
        "excluding_acrocentric_and_chrY": {
            "chromosomes": len(usual),
            "pooled": _rates(_pool(usual)) if usual else None,
        },
        "acrocentric_and_chrY_only": {
            "chromosomes": len(unusual),
            "chroms": [r["chrom"] for r in unusual],
            "pooled": _rates(_pool(unusual)) if unusual else None,
        },
        "reading": "a rate that moves with chromosome size is a rate the small chromosomes bought; "
        "chr21 and chr22 are acrocentric and chrY is nearly gene-free, and the compression lane's "
        "per-chromosome layers already looked positive there and went negative genome-wide",
    }

    return {
        "chromosomes": rows,
        "complete_chromosomes": len(complete),
        "chromosomes_complete": sorted(r["chrom"] for r in complete),
        "genome": genome,
        "spread": {k: _spread(complete, k) for k in keys},
        "controls": CONTROLS,
        "size_trend": trend,
        "requests_total": requests,
        "cells": ["K562", "HepG2", "GM12878", "IMR-90"],
        "elements_where": "data/knowledge/alphagenome/all_elements/<chrom>.json (local)",
        "evidence": "predicted: AlphaGenome deletion effect per element, per gene and per cell line; "
        "inferred: CTCF-only nodes",
    }


def main() -> int:
    out = aggregate()
    save_result("enhancer_targets_all_genome_wide", out)
    g, sp = out["genome"], out["spread"]
    print(
        f"{out['complete_chromosomes']} chromosomes complete: {g['scored']:,} elements, "
        f"{g['fraction_with_target']} name a gene (87% at matched random windows), "
        f"nearest TSS {g['coding_target_agrees_with_nearest']}, "
        f"inside the node {g['coding_target_inside_domain']} (random boundaries 0.791, +2.6 points); "
        f"{out['requests_total']:,} requests so far"
    )
    for k, s in sp.items():
        if s:
            print(f"  {k}: {s['min']} ({s['lowest']}) to {s['max']} ({s['highest']}), median {s['median']}")
    t = out["size_trend"]
    dens = t["against_element_density"]
    print(f"  over {t['chromosomes']} chromosomes, rho with length / with elements per Mb:")
    for k in sp:
        print(f"    {k}: {t[k]} / {dens[k]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
