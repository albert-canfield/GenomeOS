"""Score the reader's per-biosample family against assay depth, on the registration
committed in docs/NODES-READER-WRITER.md ("Is the rest of the reader's family assay
depth too?", 41bc7d6).

Reads only what is already on disk: data/results/reader_genome_wide.json for the
readings, data/results/epigenome_chr*.json for the mark peak counts, and
data/results/dnase_*_chr*.bed.gz for DNase peak width. No network, no model request.

Writes data/results/reader_depth_family.json.

With --normalised, scores the three readings registered on 2026-09-27 ("The reader family
normalised, to close milestone 1.1": nodes_open, genes_poised, genes_read and genes_read_open)
from the same files plus the per-chromosome node tables in data/results/reader_*_chr*.json, and
writes data/results/reader_normalised.json. Also no network and no model request.
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

RESULTS = os.path.join(ROOT, "data", "results")  # read here; written through save_result (item 12 S6)
ODD = [f"chr{i}" for i in range(1, 23, 2)]
EVEN = [f"chr{i}" for i in range(2, 23, 2)]
PERMUTATIONS = 10_000
SEED = 20260922


def _slug(cell: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", cell).strip("_")


def _in_results(names: list[str]) -> list[str]:
    return [p for p in (os.path.join(RESULTS, n) for n in names) if os.path.exists(p)]


def manifest(
    parameters: dict,
    cells: list[str],
    chroms: list[str],
    *,
    readings: list[str] = ("reader_genome_wide.json",),
    epigenome: bool = True,
    node_tables: list[str] = (),
    gencode: bool = False,
    partitions: dict | str = "n/a: no evaluation split",
) -> dict:
    """What a reader-family result read: the committed readings, the per-chromosome epigenome summaries,
    the DNase peak files of the biosamples scored, and (normalised, rarefied) the node tables."""
    from genomeos import manifest as mf

    inputs = [mf.input_entry(p) for p in _in_results(list(readings))]
    if epigenome:
        epi = sorted(
            p for p in glob.glob(os.path.join(RESULTS, "epigenome_chr*.json")) if "genome_wide" not in p
        )
        inputs.append(mf.files_entry("data/results/epigenome_chr*.json (per-chromosome peak counts)", epi))
    peaks = _in_results([f"dnase_{_slug(c)}_{ch}.bed.gz" for c in cells for ch in chroms])
    inputs.append(
        mf.files_entry(f"data/results/dnase_<biosample>_<chrom>.bed.gz, {len(cells)} biosamples", peaks)
    )
    if node_tables:
        tables = _in_results(list(node_tables))
        inputs.append(mf.files_entry("data/results/reader_<biosample>_<chrom>.json (node tables)", tables))
    if gencode:
        from genomeos.genome import default_gencode

        gffs = sorted({str(default_gencode({ch})) for ch in chroms})
        inputs.append(mf.files_entry("GENCODE v50 GFF3 per chromosome", gffs))
    return {
        "sources": [
            {
                "accession": "ENCODE DNase-seq peak files (one per biosample) and Histone ChIP-seq peak "
                "counts, as fetched by genomeos.genome.reader and genomeos.genome.epigenome",
                "version": "as on disk at the run; the files' sha256 is in inputs",
            },
            {
                "accession": "the reader's own committed readings (data/results/reader_*.json)",
                "version": "on disk",
            },
        ]
        + (
            [{"accession": "GENCODE v50 comprehensive annotation (GFF3)", "version": "v50"}]
            if gencode
            else []
        ),
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": parameters,
        "exclusions": ["epigenome_genome_wide*.json (a sum of the per-chromosome files, not read twice)"],
        "partitions": partitions,
    }


SPLIT_HALF = {
    "odd_chromosomes": "chr1, chr3, ... chr21: the first half of the registered split-half replication",
    "even_chromosomes": "chr2, chr4, ... chr22: the second half",
}


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
    from genomeos.results import save_result

    params = {"permutations": PERMUTATIONS, "seed": SEED}
    path = save_result(
        "reader_depth_family", out, manifest=manifest(params, cells, ODD + EVEN, partitions=SPLIT_HALF)
    )
    print(json.dumps({k: out[k] for k in ("family", "residual_replication", "verdict")}, indent=1))
    print(json.dumps(out["residual_vs_peak_width"], indent=1))
    print(json.dumps(out["rate"], indent=1))
    print(f"written: {path}", flush=True)


NORMALISED_SEED = 20260927
AUTOSOMES = ODD + EVEN


def per_chrom_nodes_at_reference(cells: list[str], chroms: list[str]) -> dict[str, dict[str, dict]]:
    """Per (biosample, chromosome): nodes at or above the frozen reference density, and nodes."""
    from genomeos.genome.reader import NODE_OPEN_REFERENCE_DENSITY

    out: dict[str, dict[str, dict]] = {}
    for cell in cells:
        slug = re.sub(r"[^A-Za-z0-9]+", "_", cell).strip("_")
        for chrom in chroms:
            with open(os.path.join(RESULTS, f"reader_{slug}_{chrom}.json")) as f:
                table = json.load(f)["node_table"]
            out.setdefault(cell, {})[chrom] = {
                "open": sum(1 for n in table if n["peaks_per_100kb"] >= NODE_OPEN_REFERENCE_DENSITY),
                "nodes": len(table),
            }
    return out


def per_chrom_mark(cells: list[str], mark: str) -> dict[str, dict[str, int]]:
    out: dict[str, dict[str, int]] = {c: {} for c in cells}
    for path in glob.glob(os.path.join(RESULTS, "epigenome_chr*.json")):
        if "genome_wide" in path:
            continue
        with open(path) as f:
            d = json.load(f)
        for cell in cells:
            row = (d.get("cell_types", {}).get(cell) or {}).get("peaks") or {}
            if mark in row:
                out[cell][d["chrom"]] = row[mark]["n"]
    return out


def score_reading(cells, x_of, y_of, width, seed, second=None):
    """The three registered bars for one reading. x_of/y_of map a chromosome list to per-cell
    values; `second` is the covariate of the second-assay test (default: mean DNase peak width)."""
    from genomeos.genome.reader import depth_residuals, loo_residuals

    full_chroms = None  # all 24, as the genome-wide totals are
    xs, ys = x_of(full_chroms), y_of(full_chroms)
    z = depth_residuals(xs, ys)
    loo = loo_residuals(xs, ys)
    halves = {}
    for label, keys in (("odd", ODD), ("even", EVEN)):
        halves[label] = depth_residuals(x_of(keys), y_of(keys))
    rho_h = spearman(halves["odd"], halves["even"])
    p_h = permutation_p(halves["odd"], halves["even"], rho_h, PERMUTATIONS, seed)
    rho_in = spearman(xs, z)
    rho_loo = spearman(xs, loo)
    if second is None:
        second_name, second_vals = "mean DNase peak width", [width[c]["mean_width"] for c in cells]
    else:
        second_name, second_vals = second
    rho_second = spearman(second_vals, z)
    replicates = rho_h >= 0.5 and p_h < 0.05
    removes = abs(rho_loo) < 0.4
    is_second = abs(rho_second) >= 0.7
    return {
        "raw": dict(zip(cells, ys, strict=True)),
        "covariate": dict(zip(cells, xs, strict=True)),
        "raw_rho_vs_covariate": spearman(xs, ys),
        "z": {c: round(v, 4) for c, v in zip(cells, z, strict=True)},
        "z_loo": {c: round(v, 4) for c, v in zip(cells, loo, strict=True)},
        "z_odd": {c: round(v, 4) for c, v in zip(cells, halves["odd"], strict=True)},
        "z_even": {c: round(v, 4) for c, v in zip(cells, halves["even"], strict=True)},
        "bar1_split_half": {"rho": rho_h, "p": p_h, "passes": replicates},
        "bar2_depth_removed": {
            "rho_in_sample": rho_in,
            "rho_leave_one_out": rho_loo,
            "passes": removes,
        },
        "bar3_second_assay": {"covariate": second_name, "rho": rho_second, "fires": is_second},
        "kept": replicates and removes and not is_second,
        "order_raw": sorted(cells, key=lambda c: -ys[cells.index(c)]),
        "order_normalised": sorted(cells, key=lambda c: -z[cells.index(c)]),
    }


def k562_hepg2(score: dict) -> dict:
    k, h = "K562", "HepG2"
    signs = [score[f][k] > score[f][h] for f in ("z", "z_odd", "z_even")]
    verdict = "survives" if all(signs) else ("reverses" if not any(signs) else "unclaimable")
    return {
        "raw": [score["raw"][k], score["raw"][h]],
        "z": [score["z"][k], score["z"][h]],
        "z_odd": [score["z_odd"][k], score["z_odd"][h]],
        "z_even": [score["z_even"][k], score["z_even"][h]],
        "verdict": verdict,
    }


def main_normalised() -> None:
    from genomeos.genome.reader import normalise_family

    with open(os.path.join(RESULTS, "reader_genome_wide.json")) as f:
        gw = json.load(f)
    chroms = gw["chromosomes"]
    cells = [c for c in gw["totals"] if c != "seconds"]
    all_chroms = list(chroms)
    nodes = per_chrom_nodes_at_reference(cells, all_chroms)
    k27 = per_chrom_mark(cells, "H3K27me3")
    width = dnase_peak_width(cells)

    def keys(only):
        return only if only is not None else all_chroms

    def dnase(only):
        return [float(summed(chroms, c, "peaks", keys(only))) for c in cells]

    def field(name):
        return lambda only: [float(summed(chroms, c, name, keys(only))) for c in cells]

    def mark(only):
        return [float(sum(k27[c].get(ch, 0) for ch in keys(only))) for c in cells]

    def open_ref(only):
        return [float(sum(nodes[c][ch]["open"] for ch in keys(only))) for c in cells]

    scores = {
        "nodes_open_depth_residual": score_reading(cells, dnase, open_ref, width, NORMALISED_SEED),
        "genes_poised_mark_residual": score_reading(
            cells, mark, field("genes_poised"), width, NORMALISED_SEED, ("DNase peaks", dnase(None))
        ),
        "genes_read_depth_residual": score_reading(cells, dnase, field("genes_read"), width, NORMALISED_SEED),
        "genes_read_open_depth_residual": score_reading(
            cells, dnase, field("genes_read_open"), width, NORMALISED_SEED
        ),
    }
    # nodes: floored or saturated biosamples are flagged, not ranked
    total_nodes = {c: sum(nodes[c][ch]["nodes"] for ch in all_chroms) for c in cells}
    frac = {c: round(scores["nodes_open_depth_residual"]["raw"][c] / total_nodes[c], 4) for c in cells}
    scores["nodes_open_depth_residual"]["fraction_at_reference"] = frac
    scores["nodes_open_depth_residual"]["floored_or_saturated"] = [
        c for c in cells if frac[c] < 0.05 or frac[c] > 0.95
    ]
    scores["nodes_open_depth_residual"]["median_split_raw"] = {
        c: summed(chroms, c, "nodes_open") for c in cells
    }
    # genes_poised without testis, the known under-called H3K27me3 experiment
    rest = [c for c in cells if c != "testis"]
    idx = [cells.index(c) for c in rest]

    def sub(fn):
        return lambda only: [fn(only)[i] for i in idx]

    no_testis = score_reading(
        rest, sub(mark), sub(field("genes_poised")), width, NORMALISED_SEED, ("DNase peaks", sub(dnase)(None))
    )
    gp = scores["genes_poised_mark_residual"]
    gp["without_testis"] = {
        k: no_testis[k] for k in ("bar1_split_half", "bar2_depth_removed", "bar3_second_assay", "kept", "z")
    }
    gp["testis_dependent"] = no_testis["kept"] != gp["kept"]
    if gp["testis_dependent"]:
        gp["kept"] = False

    # the published comparison, reported first
    first = {
        f: k562_hepg2(scores[f]) for f in ("genes_read_open_depth_residual", "genes_read_depth_residual")
    }
    # P3: both tissues above the cultured biosamples' median residual
    cultured = [c for c in cells if c not in ("testis", "ovary")]
    p3 = {}
    for f in ("genes_read_depth_residual", "genes_read_open_depth_residual"):
        z = scores[f]["z"]
        med = sorted(z[c] for c in cultured)[len(cultured) // 2]
        p3[f] = {
            "cultured_median_z": med,
            "testis": z["testis"],
            "ovary": z["ovary"],
            "holds": z["testis"] > med and z["ovary"] > med,
        }

    # the shipped function must agree with the scoring
    totals = {
        c: {
            "peaks": summed(chroms, c, "peaks"),
            "nodes_open_at_reference": scores["nodes_open_depth_residual"]["raw"][c],
            "genes_poised": gw["totals"][c]["genes_poised"],
            "h3k27me3_peaks": scores["genes_poised_mark_residual"]["covariate"][c],
            "genes_read": gw["totals"][c]["genes_read"],
            "genes_read_open": gw["totals"][c]["genes_read_open"],
        }
        for c in cells
    }
    shipped = normalise_family(totals)
    for f, sc in scores.items():
        assert all(abs(shipped[c][f] - sc["z"][c]) < 1e-3 for c in cells), f

    out = {
        "result": "reader_normalised",
        "registration": "docs/NODES-READER-WRITER.md, 'The reader family normalised, to close "
        "milestone 1.1: the registration', committed 63929e1 before this ran",
        "biosamples": cells,
        "k562_vs_hepg2_first": first,
        "p3_tissues_above_cultured": p3,
        "readings": scores,
        "per_biosample": shipped,
        "kept": {f: sc["kept"] for f, sc in scores.items()},
        "evidence": "derived: ENCODE DNase-seq and Histone ChIP-seq peak counts already on disk; "
        "no request made",
    }
    from genomeos.results import save_result

    path = save_result(
        "reader_normalised",
        out,
        manifest=manifest(
            {"permutations": PERMUTATIONS, "seed": NORMALISED_SEED},
            cells,
            AUTOSOMES,
            node_tables=[f"reader_{_slug(c)}_{ch}.json" for c in cells for ch in all_chroms],
            partitions=SPLIT_HALF,
        ),
    )
    print(json.dumps(first, indent=1))
    for f, sc in scores.items():
        print(
            f,
            {
                k: sc[k]
                for k in (
                    "raw_rho_vs_covariate",
                    "bar1_split_half",
                    "bar2_depth_removed",
                    "bar3_second_assay",
                    "kept",
                )
            },
        )
    print(f"written: {path}", flush=True)


def rarefied_peaks(cell: str, n: int, chroms: list[str]) -> dict[str, list[tuple[int, int, float]]]:
    """The biosample's n strongest DNase peaks genome-wide (signalValue descending, ties broken by
    genome order): the peaks a shallower experiment on the same cells would most likely still call."""
    from genomeos.genome.reader import load_peaks

    allp = []
    for ci, chrom in enumerate(chroms):
        allp += [(-v, ci, s0, e0, v) for s0, e0, v in load_peaks(cell, chrom)]
    allp.sort()
    out: dict[str, list[tuple[int, int, float]]] = {c: [] for c in chroms}
    for _, ci, s0, e0, v in allp[:n]:
        out[chroms[ci]].append((s0, e0, v))
    return out


def main_rarefied(n: int | None = None) -> None:
    """Not registered in advance: added after the registered residual reversed the K562-against-
    HepG2 read share with K562 alone at the top of the depth axis. Every biosample is cut to the
    same number of peaks (the shallowest's, GM12878's) by keeping its strongest, and the reader is
    re-run on the cached files, so depth is matched by construction rather than modelled. The node
    tables are the ones on disk; no domain is re-inferred."""
    from types import SimpleNamespace

    from genomeos.genome import Annotation, default_gencode
    from genomeos.genome import reader as rd

    with open(os.path.join(RESULTS, "reader_normalised.json")) as f:
        norm = json.load(f)
    cells = norm["biosamples"]
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    depth = {c: sum(len(rd.load_peaks(c, ch)) for ch in chroms) for c in cells}
    # default: the shallowest biosample's depth, so all thirteen take part; a larger n drops every
    # biosample that has fewer peaks than n, because a sample cannot be rarefied upwards
    n = n or min(depth.values())
    cells = [c for c in cells if depth[c] >= n]
    kept = {c: rarefied_peaks(c, n, chroms) for c in cells}
    fields = ("peaks", "genes_read", "genes_read_open", "nodes_open_at_reference")
    tot = {c: dict.fromkeys(fields, 0) for c in cells}
    original = rd.load_peaks
    try:
        for chrom in chroms:
            ann = Annotation.from_gff3(default_gencode({chrom}), {chrom})
            with open(os.path.join(RESULTS, f"reader_K562_{chrom}.json")) as f:
                table = json.load(f)["node_table"]
            doms = [
                SimpleNamespace(
                    id=t["id"],
                    start=t["start"],
                    end=t["end"],
                    length=t["end"] - t["start"],
                    coding_genes=t["coding_genes"],
                )
                for t in table
            ]
            for cell in cells:
                if not kept[cell].get(chrom):
                    # every strong peak sits elsewhere (chrY in a female or a shallow sample):
                    # nothing read, which is what the reader would say of an empty file
                    continue
                rd.load_peaks = lambda _c, ch, _k=kept[cell]: _k.get(ch, [])
                r = rd.read_chromosome(cell, chrom, ann, doms, [])
                for k in fields:
                    tot[cell][k] += r[k] or 0
            print(chrom, flush=True)
    finally:
        rd.load_peaks = original
    compare_to = {
        "genes_read": "genes_read_depth_residual",
        "genes_read_open": "genes_read_open_depth_residual",
        "nodes_open_at_reference": "nodes_open_depth_residual",
    }
    agreement = {}
    for raw, zname in compare_to.items():
        z = [norm["readings"][zname]["z"][c] for c in cells]
        r_full = [norm["readings"][zname]["raw"][c] for c in cells]
        rar = [tot[c][raw] for c in cells]
        agreement[raw] = {
            "rho_rarefied_vs_residual": spearman(rar, z),
            "rho_rarefied_vs_raw": spearman(rar, r_full),
            "k562": tot["K562"][raw],
            "hepg2": tot["HepG2"][raw],
            "k562_vs_hepg2": "K562 > HepG2" if tot["K562"][raw] > tot["HepG2"][raw] else "HepG2 >= K562",
            "min": min(rar),
            "max": max(rar),
            "order_rarefied": sorted(cells, key=lambda c: -tot[c][raw]),
        }
    out = {
        "result": "reader_rarefied",
        "registered": False,
        "why": main_rarefied.__doc__.strip(),
        "peaks_per_biosample": n,
        "biosamples": cells,
        "full_depth": depth,
        "rarefied": tot,
        "agreement": agreement,
        "evidence": "derived: ENCODE DNase-seq peaks already on disk, strongest n kept; no request made",
    }
    from genomeos.results import load_result, save_result

    existing = load_result("reader_rarefied") or {}
    runs = existing.get("runs", {})
    runs[str(n)] = out
    path = save_result(
        "reader_rarefied",
        {"result": "reader_rarefied", "registered": False, "runs": runs},
        manifest=manifest(
            {"peaks_kept_per_biosample": n, "runs_in_file": sorted(runs, key=int)},
            cells,
            chroms,
            readings=["reader_normalised.json"],
            epigenome=False,
            node_tables=[f"reader_K562_{ch}.json" for ch in chroms],
            gencode=True,
        ),
    )
    print(json.dumps(agreement, indent=1))
    print(f"written: {path}", flush=True)


if __name__ == "__main__":
    if "--rarefied" in sys.argv:
        _rest = [a for a in sys.argv[1:] if a.isdigit()]
        main_rarefied(int(_rest[0]) if _rest else None)
    elif "--normalised" in sys.argv:
        main_normalised()
    else:
        main()
