# SPDX-License-Identifier: AGPL-3.0-or-later
"""The epigenome layer on real data: fetch, per-chromosome summaries, and three questions.

    uv run python scripts/epigenome.py fetch [--signal chr21,chr22] [--marks H3K27me3,...]
        [--methylation all] [--workers 6]
    uv run python scripts/epigenome.py summary chr21 chr22 ...
    uv run python scripts/epigenome.py genome         # the summaries rolled up
    uv run python scripts/epigenome.py direction      # marks against registry class, chr21 and chr22
    uv run python scripts/epigenome.py direction-transfer  # fit on chr21+chr22, score held-out chromosomes
    uv run python scripts/epigenome.py hox            # H3K27me3 over the HOX clusters, matched promoters
    uv run python scripts/epigenome.py fossil         # methylation of the fossil tier, GC and CpG matched
    uv run python scripts/epigenome.py reader-check   # read, poised and silent genes against measured RNA
    uv run python scripts/epigenome.py alu            # Alu against L1 remains, age, region and WCGW matched

No model is called: the deletion effects are read from the per-element tables the
all-enhancer chain already wrote (data/knowledge/alphagenome/all_elements).
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import random
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos.genome import epigenome as ep  # noqa: E402
from genomeos.genome import reader  # noqa: E402

HOX_CHROMS = ("chr7", "chr17", "chr12", "chr2")


def _task(kind: str, cell: str, mark: str | None, chrom: str | None) -> str:
    m = ep.load_manifest()
    row = m["cell_types"][cell]
    if kind == "peaks":
        e = row["marks"][mark]
        if all(ep.peaks_path(cell, mark, c).exists() for c in ep.CHROMS):
            return f"{cell} {mark} peaks cached"
        kept = ep.fetch_mark_peaks(e)
        return f"{cell} {mark} peaks {e['peaks']['accession']}: {sum(kept.values()):,}"
    if kind == "signal":
        e = row["marks"][mark]
        if ep.load_signal_profile(cell, mark, chrom) is not None:
            return f"{cell} {mark} {chrom} signal cached"
        c = ep.fetch_signal_profile(e, chrom)
        return f"{cell} {mark} {chrom} signal {c['mb']} MB {c['seconds']} s"
    e = row["methylation"]
    if ep.load_methylation_profile(cell, chrom) is not None:
        return f"{cell} {chrom} WGBS cached"
    c = ep.fetch_methylation_profile(e, chrom)
    return f"{cell} {chrom} WGBS {c['calls']:,} calls {c['mb']} MB {c['seconds']} s"


def fetch(argv: list[str]) -> None:
    opts = dict(zip(argv[::2], argv[1::2], strict=False))
    workers = int(opts.get("--workers", 6))
    # all 24 by default: a cached profile costs nothing to skip, and the connection-reusing range
    # reader took the whole genome from "too slow to attempt" to a couple of hours at eight workers
    sig = opts.get("--signal", "all")
    sig_chroms = list(ep.CHROMS) if sig == "all" else [c for c in sig.split(",") if c]
    meth = opts.get("--methylation", "all")
    meth_chroms = list(ep.CHROMS) if meth == "all" else [c for c in meth.split(",") if c]
    m = ep.load_manifest() or ep.build_manifest()
    tasks = []
    for cell, row in m["cell_types"].items():
        for mark, e in row["marks"].items():
            if e.get("peaks"):
                tasks.append(("peaks", cell, mark, None))
    marks = opts.get("--marks", ",".join(ep.MARKS)).split(",")
    for chrom in sig_chroms:
        for cell, row in m["cell_types"].items():
            for mark, e in row["marks"].items():
                if e.get("signal") and mark in marks:
                    tasks.append(("signal", cell, mark, chrom))
    for chrom in meth_chroms:
        for cell, row in m["cell_types"].items():
            if "cpg" in row["methylation"]:
                tasks.append(("methylation", cell, None, chrom))
    print(f"{len(tasks)} tasks, {workers} workers", flush=True)
    failed = 0
    with ProcessPoolExecutor(workers) as pool:
        futs = {pool.submit(_task, *t): t for t in tasks}
        for i, f in enumerate(as_completed(futs), 1):
            try:
                print(f"[{i}/{len(tasks)}] {f.result()}", flush=True)
            except Exception as exc:  # noqa: BLE001 - one file failing must not stop the rest
                failed += 1
                print(f"[{i}/{len(tasks)}] FAILED {futs[f]}: {exc!r}", flush=True)
    print(f"done, {failed} failed", flush=True)


# ------------------------------------------------------------------------------------------
# Question 1: do the measured marks predict the model's deletion direction better than the
# registry class?
# ------------------------------------------------------------------------------------------

DIRECTION_CELLS = ("K562", "HepG2", "GM12878", "IMR-90")  # the four lines the chain scored per cell
ACTS = 0.1  # |log2 fold change| the project calls an effect
FLANK = 500  # marks sit on the nucleosomes around an element, not on its open core


def _element_rows(chrom: str) -> list[dict]:
    table = json.loads(Path(f"data/knowledge/alphagenome/all_elements/{chrom}.json").read_text())
    cls = {}
    with gzip.open(f"data/results/ccres_{chrom}.bed.gz", "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            cls[f[3]] = (f[4], f[5] == "1")
    rows = []
    for e in table:
        by = e.get("predicted_by_cell") or {}
        if e["id"] not in cls or not by:
            continue
        rows.append({**e, "cls": cls[e["id"]][0], "ctcf": cls[e["id"]][1], "by": by})
    return rows


def _features(layer: ep.Layer, cell: str, start: int, end: int) -> dict:
    mid = (start + end) // 2
    lo, hi = mid - FLANK, mid + FLANK
    out = {}
    dn = layer.dnase.get(cell)
    hits = dn.overlapping(start, end) if dn else []
    out["dnase_peak"] = 1.0 if hits else 0.0
    out["dnase_signal"] = math.log1p(max((v for _, _, v in hits), default=0.0))
    for m in ep.MARKS:
        idx = layer.peaks.get((cell, m))
        out[f"{m}_peak"] = 1.0 if idx and idx.overlapping(lo, hi) else 0.0
        prof = layer.signal.get((cell, m))
        fc = ep.profile_mean(prof, lo, hi) if prof is not None else None
        out[f"{m}_fc"] = math.log1p(fc) if fc is not None else None
    cols = layer.meth.get(cell)
    mm = ep.methylation_over(cols, lo, hi) if cols is not None else {"fraction": None}
    out["methylation"] = mm["fraction"]
    out["state"] = ep.chromatin_state({m: bool(out[f"{m}_peak"]) for m in ep.MARKS})
    return out


def _impute(rows: list[dict], keys: list[str]) -> None:
    for k in keys:
        vals = sorted(r[k] for r in rows if r[k] is not None)
        med = vals[len(vals) // 2] if vals else 0.0
        for r in rows:
            r[k + "_missing"] = 1.0 if r[k] is None else 0.0
            if r[k] is None:
                r[k] = med


def direction(argv: list[str]) -> None:
    from genomeos.results import save_result

    chroms = argv or ["chr21", "chr22"]
    units: list[dict] = []
    for chrom in chroms:
        layer = ep.Layer(chrom, DIRECTION_CELLS)
        for e in _element_rows(chrom):
            gid = int(hashlib.md5(e["id"].encode()).hexdigest()[:8], 16)  # noqa: S324 - fold id only
            for cell in DIRECTION_CELLS:
                if cell not in e["by"]:
                    continue
                f = _features(layer, cell, e["start"], e["end"])
                units.append(
                    {
                        "id": e["id"],
                        "group": gid,
                        "chrom": chrom,
                        "cell": cell,
                        "cls": e["cls"],
                        "ctcf": e["ctcf"],
                        "lfc": e["by"][cell],
                        **f,
                    }
                )
        print(f"{chrom}: {len(units):,} element-cell units so far", flush=True)
    mark_keys = [f"{m}_fc" for m in ep.MARKS] + ["methylation"]
    _impute(units, mark_keys)
    # controls that keep every marginal: the same element's marks in another of the four lines, and
    # marks shuffled among units of the same chromosome, cell, registry class and DNase call
    by_key = {(u["id"], u["cell"]): u for u in units}
    rot = dict(zip(DIRECTION_CELLS, DIRECTION_CELLS[1:] + DIRECTION_CELLS[:1], strict=True))
    mark_cols = [f"{m}_peak" for m in ep.MARKS] + mark_keys + ["methylation_missing"]
    rng = random.Random(21)
    strata: dict[tuple, list[dict]] = {}
    for u in units:
        strata.setdefault((u["chrom"], u["cell"], u["cls"], u["dnase_peak"]), []).append(u)
    for group in strata.values():
        donors = group[:]
        rng.shuffle(donors)
        for u, d in zip(group, donors, strict=True):
            u["perm"] = {c: d[c] for c in mark_cols}
    for u in units:
        other = by_key.get((u["id"], rot[u["cell"]]))
        u["other"] = {c: other[c] for c in mark_cols} if other else None
    units = [u for u in units if u["other"] is not None]

    def x_of(u: dict, model: str) -> list[float]:
        base = [1.0 if u["cls"] == "pELS" else 0.0, 1.0 if u["ctcf"] else 0.0]
        base += [1.0 if u["cell"] == c else 0.0 for c in DIRECTION_CELLS[1:]]
        base += [1.0 if u["chrom"] == c else 0.0 for c in chroms[1:]]
        openness = [u["dnase_peak"], u["dnase_signal"]]
        src = {"marks": u, "other_cell": u["other"], "permuted": u["perm"]}
        out = list(base)
        if model != "registry":
            out += openness
        if model in ("marks", "other_cell", "permuted", "marks_only"):
            d = src.get(model, u)
            vals = [d[c] for c in mark_cols]
            out = (base[2:] + vals) if model == "marks_only" else out + vals
        return out

    models = ("registry", "registry+dnase", "marks", "other_cell", "permuted", "marks_only")
    labels = {
        "registry": "registry class (pELS/dELS, CTCF-bound)",
        "registry+dnase": "registry class + DNase in the cell (reader v1)",
        "marks": "registry + DNase + measured marks and methylation in the cell",
        "other_cell": "control: same, marks from another of the four lines",
        "permuted": "control: same, marks shuffled within chromosome x cell x class x DNase",
        "marks_only": "measured marks and methylation alone (with cell and chromosome)",
    }
    tasks = {
        "acts": ([1.0 if abs(u["lfc"]) >= ACTS else 0.0 for u in units], list(range(len(units)))),
    }
    acting = [i for i, u in enumerate(units) if abs(units[i]["lfc"]) >= ACTS]
    tasks["silencer_among_acting"] = ([1.0 if units[i]["lfc"] > 0 else 0.0 for i in acting], acting)
    tasks["magnitude"] = ([math.log10(abs(u["lfc"]) + 1e-3) for u in units], list(range(len(units))))
    results: dict = {}
    preds: dict = {}
    for task, (y, idx) in tasks.items():
        groups = [units[i]["group"] for i in idx]
        results[task] = {"n": len(idx)}
        for model in models:
            x = ep.standardise([x_of(units[i], model) for i in idx])
            p = ep.cross_validated(x, y, groups, folds=5, lam=1.0)
            preds[(task, model)] = p
            score = ep.spearman(p, y) if task == "magnitude" else ep.auc(p, [int(v) for v in y])
            results[task][model] = round(score, 4) if score is not None else None
            print(f"{task:24s} {model:16s} {results[task][model]}", flush=True)
        # element bootstrap of the two differences that answer the question
        rng2 = random.Random(7)
        by_group: dict[int, list[int]] = {}
        for j, g in enumerate(groups):
            by_group.setdefault(g, []).append(j)
        keys = list(by_group)
        diffs: dict[str, list[float]] = {
            "marks-registry+dnase": [],
            "marks-other_cell": [],
            "marks-permuted": [],
        }
        for _ in range(200):
            sample = [j for _g in (rng2.choice(keys) for _ in keys) for j in by_group[_g]]
            ys = [y[j] for j in sample]
            for name in diffs:
                a, b = name.split("-", 1)
                if task == "magnitude":
                    sa = ep.spearman([preds[(task, a)][j] for j in sample], ys)
                    sb = ep.spearman([preds[(task, b)][j] for j in sample], ys)
                else:
                    yi = [int(v) for v in ys]
                    sa = ep.auc([preds[(task, a)][j] for j in sample], yi)
                    sb = ep.auc([preds[(task, b)][j] for j in sample], yi)
                diffs[name].append(sa - sb)
        results[task]["bootstrap_95"] = {
            k: [
                round(sorted(v)[4], 4),
                round(sorted(v)[194], 4),
                round(sum(1 for d in v if d <= 0) / len(v), 3),
            ]
            for k, v in diffs.items()
        }
    # the same question read as a table: silencer share by registry class and by measured state
    table: dict = {"by_class": {}, "by_state": {}, "by_h3k27me3": {}}
    for u in units:
        if abs(u["lfc"]) < ACTS:
            continue
        for key, val in (
            ("by_class", u["cls"]),
            ("by_state", u["state"]),
            ("by_h3k27me3", "H3K27me3 peak" if u["H3K27me3_peak"] else "no H3K27me3 peak"),
        ):
            t = table[key].setdefault(val, {"acting": 0, "silencer": 0, "abs_lfc": []})
            t["acting"] += 1
            t["silencer"] += u["lfc"] > 0
            t["abs_lfc"].append(abs(u["lfc"]))
    for key in table.values():
        for t in key.values():
            a = sorted(t.pop("abs_lfc"))
            t["silencer_share"] = round(t["silencer"] / t["acting"], 4)
            t["median_abs_lfc"] = round(a[len(a) // 2], 4)
    state_all: dict[str, dict] = {}
    for u in units:
        t = state_all.setdefault(u["state"], {"units": 0, "acting": 0})
        t["units"] += 1
        t["acting"] += abs(u["lfc"]) >= ACTS
    for t in state_all.values():
        t["acting_share"] = round(t["acting"] / t["units"], 4)
    table["acting_by_state"] = state_all
    out = {
        "chromosomes": chroms,
        "cell_types": list(DIRECTION_CELLS),
        "units": len(units),
        "elements": len({u["id"] for u in units}),
        "acts_threshold_abs_log2fc": ACTS,
        "window": f"element midpoint +/- {FLANK} bp for marks and methylation; the element for DNase",
        "outcome": "AlphaGenome predicted_by_cell: signed log2 fold change of the gene the deletion moves "
        "most, per cell line; positive = expression rises on deletion (silencer-like)",
        "metrics": "acts and silencer: out-of-fold AUC; magnitude: Spearman of out-of-fold prediction with "
        "log10 |lfc|; 5 folds grouped by element; ridge linear model on standardised features",
        "models": labels,
        "results": results,
        "bootstrap": "95% interval of the difference over 200 element resamples, and the share of resamples "
        "with no gain",
        "tables": table,
        "caveat": "the outcome is a model's prediction, and AlphaGenome was trained on ENCODE histone "
        "ChIP-seq, DNase and expression for these same cell lines: marks predicting its direction partly "
        "measures what the model learned from those tracks, not the biology of the element",
        "evidence": {
            "marks": ep.EVIDENCE_MARK,
            "methylation": ep.EVIDENCE_METHYLATION,
            "outcome": "predicted: AlphaGenome deletion scoring, chain tables, no new model call",
        },
    }
    save_result(
        "epigenome_direction_chr21_chr22" if chroms == ["chr21", "chr22"] else "epigenome_direction", out
    )
    print("saved", flush=True)


# ------------------------------------------------------------------------------------------
# Shared: annotation, promoters and sequence composition
# ------------------------------------------------------------------------------------------


def _annotation(chrom: str):
    from genomeos.genome import Annotation, default_gencode

    gff = default_gencode({chrom})
    if gff is None:
        raise SystemExit(f"{chrom}: no local GENCODE models (genomeos data fetch --chrom {chrom})")
    return Annotation.from_gff3(gff, {chrom})


def _promoters(chrom: str) -> list[dict]:
    """Coding promoters (TSS +/- 1 kb) of one chromosome with GC and CpG density."""
    from genomeos.coords import Locus
    from genomeos.genome import IndexedGenome

    ann = _annotation(chrom)
    g = IndexedGenome(f"data/reference/{chrom}.fa.gz")
    out = []
    try:
        for gene in ann.genes.values():
            if gene.locus.chrom != chrom or gene.type != "protein_coding":
                continue
            t = ep._tss(gene)
            lo, hi = max(0, t - ep.PROMOTER_WINDOW), t + ep.PROMOTER_WINDOW
            gc, cpg = ep.gc_cpg(str(g.fetch(Locus(chrom, lo, min(hi, g.lengths[chrom])))))
            if gc != gc:
                continue
            out.append(
                {"gene": gene.symbol, "chrom": chrom, "tss": t, "lo": lo, "hi": hi, "gc": gc, "cpg": cpg}
            )
    finally:
        g.close()
    return out


def _repression(layer: ep.Layer, cell: str, lo: int, hi: int) -> dict:
    called = layer._marks_called(cell, lo, hi)
    dn = layer.dnase.get(cell)
    cols = layer.meth.get(cell)
    meth = ep.methylation_over(cols, lo, hi) if cols is not None else None
    frac = meth["fraction"] if meth and meth["cpg_calls_covered"] >= 4 else None
    prof = layer.signal.get((cell, "H3K27me3"))
    return {
        "open": None if dn is None else bool(dn.overlapping(lo, hi)),
        "called": called,
        "state": ep.chromatin_state(called),
        "methylation": frac,
        "h3k27me3_fc": ep.profile_mean(prof, lo, hi) if prof is not None else None,
    }


def _silence_cause(r: dict) -> str:
    """The first measured mechanism that accounts for a closed promoter, in a fixed order."""
    c = r["called"]
    if c.get("H3K27me3"):
        return "H3K27me3"
    if c.get("H3K9me3"):
        return "H3K9me3"
    if r["methylation"] is not None and r["methylation"] >= 0.5:
        return "CpG methylation >= 0.5"
    if r["methylation"] is None:
        return "none measured (no methylation data)"
    return "none of the three"


# ------------------------------------------------------------------------------------------
# Question 2: does H3K27me3 explain the loci that behave as repressed, the HOX clusters first?
# ------------------------------------------------------------------------------------------

HOX = re.compile(r"^HOX[ABCD]\d+$")


def hox(argv: list[str]) -> None:
    from genomeos.results import save_result

    chroms = argv or list(HOX_CHROMS)
    cells = list(ep.CELL_TYPES)
    proms: list[dict] = []
    for chrom in chroms:
        proms += _promoters(chrom)
    per_cell: dict = {}
    clusters: dict = {}
    for chrom in chroms:
        layer = ep.Layer(chrom, cells)
        mine = [p for p in proms if p["chrom"] == chrom]
        for p in mine:
            p["by_cell"] = {c: _repression(layer, c, p["lo"], p["hi"]) for c in cells}
        hx = [p for p in mine if HOX.match(p["gene"])]
        if not hx:
            continue
        name = "HOX" + hx[0]["gene"][3]
        lo, hi = min(p["tss"] for p in hx), max(p["tss"] for p in hx)
        cl = {"chrom": chrom, "start": lo, "end": hi, "genes": sorted(p["gene"] for p in hx), "cells": {}}
        for c in cells:
            idx = layer.peaks.get((c, "H3K27me3"))
            covered = idx.covered_bp(lo, hi) if idx else None
            prof = layer.signal.get((c, "H3K27me3"))
            fc_cluster = ep.profile_mean(prof, lo, hi) if prof is not None else None
            chrom_vals = [v for v in prof if v == v] if prof is not None else []
            fc_chrom = sum(chrom_vals) / len(chrom_vals) if chrom_vals else None
            k4 = layer.peaks.get((c, "H3K4me3"))
            cl["cells"][c] = {
                "h3k27me3_bp_share": round(covered / (hi - lo), 4) if covered is not None else ep.UNKNOWN,
                "h3k27me3_fc_cluster": round(fc_cluster, 3) if fc_cluster is not None else ep.UNKNOWN,
                "h3k27me3_fc_chromosome_mean": round(fc_chrom, 3) if fc_chrom is not None else ep.UNKNOWN,
                "h3k4me3_bp_share": round(k4.covered_bp(lo, hi) / (hi - lo), 4) if k4 else ep.UNKNOWN,
                "promoters": {
                    p["gene"]: {
                        "open": p["by_cell"][c]["open"],
                        "state": p["by_cell"][c]["state"],
                        "methylation": p["by_cell"][c]["methylation"],
                    }
                    for p in sorted(hx, key=lambda q: q["tss"])
                },
            }
        clusters[name] = cl
    # matched control: non-HOX coding promoters in the same GC x CpG-density stratum
    for c in cells:
        rows = []
        for group, silent_only in (("all", False), ("silent", True)):
            hx = [
                p
                for p in proms
                if HOX.match(p["gene"]) and (not silent_only or p["by_cell"][c]["open"] is False)
            ]
            ctl = [
                p
                for p in proms
                if not HOX.match(p["gene"]) and (not silent_only or p["by_cell"][c]["open"] is False)
            ]
            by_s: dict[str, list[dict]] = {}
            for p in ctl:
                by_s.setdefault(ep.stratum(p["gc"], p["cpg"]), []).append(p)
            n_hox = n_used = 0
            k27_hox = k27_exp = 0.0
            biv_hox = biv_exp = 0.0
            for p in hx:
                pool = by_s.get(ep.stratum(p["gc"], p["cpg"]), [])
                n_hox += 1
                if len(pool) < 5:
                    continue
                n_used += 1
                k27_hox += bool(p["by_cell"][c]["called"].get("H3K27me3"))
                biv_hox += p["by_cell"][c]["state"] == "bivalent"
                k27_exp += sum(bool(q["by_cell"][c]["called"].get("H3K27me3")) for q in pool) / len(pool)
                biv_exp += sum(q["by_cell"][c]["state"] == "bivalent" for q in pool) / len(pool)
            rows.append(
                {
                    "promoters": group,
                    "hox": n_hox,
                    "hox_matched": n_used,
                    "h3k27me3_share_hox": round(k27_hox / n_used, 4) if n_used else None,
                    "h3k27me3_share_matched": round(k27_exp / n_used, 4) if n_used else None,
                    "bivalent_share_hox": round(biv_hox / n_used, 4) if n_used else None,
                    "bivalent_share_matched": round(biv_exp / n_used, 4) if n_used else None,
                }
            )
        causes: dict[str, dict[str, int]] = {"silent": {}, "read": {}}
        for p in proms:
            r = p["by_cell"][c]
            if r["open"] is None:
                continue
            key = "read" if r["open"] else "silent"
            cause = _silence_cause(r)
            causes[key][cause] = causes[key].get(cause, 0) + 1
        per_cell[c] = {"matched": rows, "promoter_mechanism": causes}
    out = {
        "chromosomes": chroms,
        "coding_promoters": len(proms),
        "hox_promoters": sum(1 for p in proms if HOX.match(p["gene"])),
        "clusters": clusters,
        "cell_types": per_cell,
        "control": "non-HOX coding promoters on the same chromosomes in the same GC (0.05) x CpG per 100 bp "
        "(0.5, 1, 2, 4, 8) stratum; expectation averaged over each HOX promoter's stratum; strata with "
        "fewer than 5 controls are dropped",
        "mechanism_order": "H3K27me3 peak, else H3K9me3 peak, else CpG methylation >= 0.5 (4+ covered "
        "calls), else none; for promoters the DNase reader calls silent, beside the same count for read ones",
        "evidence": {
            "marks": ep.EVIDENCE_MARK,
            "methylation": ep.EVIDENCE_METHYLATION,
            "open": "DNase reader",
        },
    }
    save_result("epigenome_hox", out)
    for c in cells:
        a, sil = per_cell[c]["matched"]
        print(
            f"{c:24s} HOX H3K27me3 {a['h3k27me3_share_hox']} vs matched {a['h3k27me3_share_matched']}; "
            f"silent HOX {sil['hox']} -> {sil['h3k27me3_share_hox']} vs {sil['h3k27me3_share_matched']}",
            flush=True,
        )


# ------------------------------------------------------------------------------------------
# Question 3: is the fossil tier more methylated than the other tiers, matched for GC and CpG?
# ------------------------------------------------------------------------------------------

WINDOW = 1_000
MIN_WINDOW_CALLS = 2


def fossil(argv: list[str]) -> None:
    from genomeos.coords import Locus
    from genomeos.genome import IndexedGenome
    from genomeos.results import load_result, save_result

    m = ep.load_manifest()
    cells = [c for c, r in m["cell_types"].items() if "cpg" in r["methylation"]]
    wanted = argv or list(ep.CHROMS)
    chroms = [
        c
        for c in wanted
        if load_result(f"budget_{c}")
        and all(ep.load_methylation_profile(x, c) is not None for x in cells[:1])
        and all(ep.methylation_path(x, c).exists() for x in cells)
    ]
    tiers = ("fossil", "neutral", "regulatory", "constrained_unknown", "structural")
    # acc[cell][tier][block_key][stratum] = [sum_fraction, windows_measured]
    acc = {c: {t: {} for t in tiers} for c in cells}
    k9 = {c: {t: {} for t in tiers} for c in ep.CELL_TYPES}  # H3K9me3 peak overlap, same shape
    tally = {t: {"blocks": 0, "bp": 0, "windows": 0, "ref_cpg": 0} for t in tiers}
    measured = {c: {t: {"windows": 0, "covered_calls": 0} for t in tiers} for c in cells}
    families: dict = {c: {} for c in cells}
    genome_tier_bp = (load_result("budget_genome_wide") or {}).get("by_tier", {})
    for chrom in chroms:
        budget = load_result(f"budget_{chrom}")
        cols = {c: ep.load_methylation_profile(c, chrom) for c in cells}
        k9idx = {}
        for c in ep.CELL_TYPES:
            pk = ep.load_mark_peaks(c, "H3K9me3", chrom)
            k9idx[c] = reader.PeakIndex(pk) if pk else None
        g = IndexedGenome(f"data/reference/{chrom}.fa.gz")
        length = g.lengths[chrom]
        for bi, blk in enumerate(budget["blocks"]):
            tier = blk["guess"]["tier"]
            if tier not in tally or blk["class"] == "gap":
                continue
            key = f"{chrom}:{bi}"
            tally[tier]["blocks"] += 1
            tally[tier]["bp"] += blk["end"] - blk["start"]
            w0 = -(-blk["start"] // WINDOW) * WINDOW
            for ws in range(w0, min(blk["end"], length) - WINDOW + 1, WINDOW):
                seq = str(g.fetch(Locus(chrom, ws, ws + WINDOW)))
                if seq.upper().count("N") > WINDOW // 10:
                    continue
                gc, cpg = ep.gc_cpg(seq)
                st = ep.stratum(gc, cpg)
                tally[tier]["windows"] += 1
                tally[tier]["ref_cpg"] += ep.cpg_sites(seq)
                for c in cells:
                    mo = ep.methylation_over(cols[c], ws, ws + WINDOW)
                    measured[c][tier]["covered_calls"] += mo["cpg_calls_covered"]
                    if mo["cpg_calls_covered"] < MIN_WINDOW_CALLS:
                        continue
                    measured[c][tier]["windows"] += 1
                    cell_blk = acc[c][tier].setdefault(key, {})
                    cs = cell_blk.setdefault(st, [0.0, 0])
                    cs[0] += mo["fraction"]
                    cs[1] += 1
                    if tier == "fossil":
                        fam = blk["class"].replace("interspersed_repeat_", "") or "unspecified"
                        fs = families[c].setdefault(fam, {}).setdefault(key, {}).setdefault(st, [0.0, 0])
                        fs[0] += mo["fraction"]
                        fs[1] += 1
                for c in ep.CELL_TYPES:
                    if k9idx[c] is None:
                        continue
                    ks = k9[c][tier].setdefault(key, {}).setdefault(st, [0.0, 0])
                    ks[0] += 1.0 if k9idx[c].overlapping(ws, ws + WINDOW) else 0.0
                    ks[1] += 1
        g.close()
        print(f"{chrom}: windows so far {sum(t['windows'] for t in tally.values()):,}", flush=True)

    def pooled(blocks: dict) -> dict[str, list[float]]:
        out: dict[str, list[float]] = {}
        for strata in blocks.values():
            for st, (sm, n) in strata.items():
                o = out.setdefault(st, [0.0, 0])
                o[0] += sm
                o[1] += n
        return out

    def raw_mean(blocks: dict) -> float | None:
        p = pooled(blocks)
        n = sum(v[1] for v in p.values())
        return sum(v[0] for v in p.values()) / n if n else None

    def matched(a_blocks: dict, b_blocks: dict) -> tuple[float | None, float]:
        """Mean of a minus mean of b inside strata, weighted by a's windows; share of a's windows used."""
        pa, pb = pooled(a_blocks), pooled(b_blocks)
        tot = sum(v[1] for v in pa.values())
        num = den = 0.0
        for st, (sa, na) in pa.items():
            sb, nb = pb.get(st, (0.0, 0))
            if na < 5 or nb < 5:
                continue
            num += na * (sa / na - sb / nb)
            den += na
        return (num / den if den else None), (den / tot if tot else 0.0)

    def boot(a_blocks: dict, b_blocks: dict, reps: int = 200, seed: int = 3) -> list[float] | None:
        rng = random.Random(seed)
        ka, kb = list(a_blocks), list(b_blocks)
        if not ka or not kb:
            return None
        vals = []
        for _ in range(reps):
            ra = {i: a_blocks[k] for i, k in enumerate(rng.choice(ka) for _ in ka)}
            rb = {i: b_blocks[k] for i, k in enumerate(rng.choice(kb) for _ in kb)}
            d, _ = matched(ra, rb)
            if d is not None:
                vals.append(d)
        if len(vals) < reps // 2:
            return None
        vals.sort()
        lo_i, hi_i = int(0.025 * len(vals)), int(0.975 * len(vals)) - 1
        return [
            round(vals[lo_i], 4),
            round(vals[hi_i], 4),
            round(sum(1 for v in vals if v <= 0) / len(vals), 3),
        ]

    results: dict = {}
    for c in cells:
        row: dict = {"raw_mean": {}, "matched_fossil_minus": {}, "coverage": {}}
        for t in tiers:
            rm = raw_mean(acc[c][t])
            row["raw_mean"][t] = round(rm, 4) if rm is not None else None
            w = tally[t]["windows"]
            row["coverage"][t] = {
                "windows_measured_share": round(measured[c][t]["windows"] / w, 4) if w else None,
                "cpg_strand_calls_covered_share": round(
                    measured[c][t]["covered_calls"] / (2 * tally[t]["ref_cpg"]), 4
                )
                if tally[t]["ref_cpg"]
                else None,
            }
        for t in tiers[1:]:
            d, used = matched(acc[c]["fossil"], acc[c][t])
            row["matched_fossil_minus"][t] = {
                "difference": round(d, 4) if d is not None else None,
                "fossil_windows_in_shared_strata": round(used, 4),
                "bootstrap_95_and_share_no_excess": boot(acc[c]["fossil"], acc[c][t]),
            }
        fam_rows = {}
        for fam, blocks in families[c].items():
            d, used = matched(blocks, acc[c]["neutral"])
            rm = raw_mean(blocks)
            fam_rows[fam] = {
                "raw_mean": round(rm, 4) if rm is not None else None,
                "matched_minus_neutral": round(d, 4) if d is not None else None,
                "windows": sum(v[1] for v in pooled(blocks).values()),
            }
        row["fossil_families"] = fam_rows
        results[c] = row
        print(
            c,
            row["raw_mean"],
            {k: v["difference"] for k, v in row["matched_fossil_minus"].items()},
            flush=True,
        )
    k9res = {}
    for c in ep.CELL_TYPES:
        if not any(k9[c]["fossil"].values()):
            continue
        k9res[c] = {
            "share_windows_with_peak": {
                t: (round(raw_mean(k9[c][t]), 4) if raw_mean(k9[c][t]) is not None else None) for t in tiers
            },
            "matched_fossil_minus": {
                t: round(matched(k9[c]["fossil"], k9[c][t])[0] or 0.0, 4) for t in tiers[1:]
            },
        }
    fossil_bp_read = tally["fossil"]["bp"]
    fossil_bp_genome = (genome_tier_bp.get("fossil") or {}).get("bp")
    out = {
        "chromosomes": chroms,
        "cell_types_with_wgbs": cells,
        "cell_types_without": [c for c in ep.CELL_TYPES if c not in cells],
        "window": f"{WINDOW} bp windows on a fixed grid inside each block, over 10% N dropped; a window is "
        f"measured with {MIN_WINDOW_CALLS}+ CpG strand calls of {ep.MIN_COVERAGE}+ reads",
        "match": "GC (0.05 bins) x CpG per 100 bp (0.5, 1, 2, 4, 8) strata; the fossil-minus-tier "
        "difference is the per-stratum difference of mean window methylation weighted by fossil windows, "
        "strata with 5+ windows on both sides; bootstrap resamples blocks, not windows",
        "tiers": tally,
        "fossil_bp_read": fossil_bp_read,
        "fossil_bp_genome": fossil_bp_genome,
        "fossil_bp_share_of_genome_tier": round(fossil_bp_read / fossil_bp_genome, 4)
        if fossil_bp_genome
        else None,
        "methylation": results,
        "h3k9me3": k9res,
        "evidence": {
            "methylation": ep.EVIDENCE_METHYLATION,
            "tiers": "inferred: budget tiers (attribution/budget.py), sequence class and Zoonomia constraint",
            "h3k9me3": ep.EVIDENCE_MARK + ", replicated peaks",
        },
    }
    save_result("epigenome_fossil" if not argv else "epigenome_fossil_" + "_".join(chroms), out)
    print("saved", flush=True)


# ------------------------------------------------------------------------------------------
# Per-chromosome summaries
# ------------------------------------------------------------------------------------------


def summary(argv: list[str]) -> None:
    from genomeos.genome.regulatory import load_ccres
    from genomeos.results import load_result, save_result

    for chrom in argv or ["chr21", "chr22"]:
        ann = _annotation(chrom)
        s = ep.summarise_chromosome(chrom, ann, load_ccres(chrom), budget=load_result(f"budget_{chrom}"))
        s["manifest"] = "epigenome_manifest"
        p = save_result(f"epigenome_{chrom}", s)
        print(f"{chrom}: {p} {p.stat().st_size / 1e3:.0f} KB", flush=True)


def genome(argv: list[str]) -> None:
    """Roll the per-chromosome summaries up into one genome-wide table per cell type."""
    from genomeos.results import load_result, save_result

    chroms = [c for c in ep.CHROMS if load_result(f"epigenome_{c}")]
    cells: dict = {}
    for chrom in chroms:
        s = load_result(f"epigenome_{chrom}")
        for cell, row in s["cell_types"].items():
            out = cells.setdefault(
                cell,
                {
                    "marks_measured": row["marks_measured"],
                    "methylation": row["methylation"],
                    "promoters": {},
                    "registry_classes": {},
                    "methylation_by_tier": {},
                    "peaks": {},
                },
            )
            for m, pk in row.get("peaks", {}).items():
                o = out["peaks"].setdefault(m, {"n": 0, "bp": 0})
                o["n"] += pk["n"]
                o["bp"] += pk["bp"]
            for group, b in row["promoters"].items():
                o = out["promoters"].setdefault(group, {"n": 0, "states": {}})
                o["n"] += b["n"]
                for m in ep.MARKS:
                    if isinstance(b.get(m), dict):
                        o[m] = o.get(m, 0) + b[m]["n"]
                for st, k in b.get("states", {}).items():
                    o["states"][st] = o["states"].get(st, 0) + k
            for cls, b in row["registry_classes"].items():
                o = out["registry_classes"].setdefault(cls, {"n": 0, "open": 0, "states": {}})
                o["n"] += b["n"]
                o["open"] += b["open"]
                for m in ep.MARKS:
                    if isinstance(b.get(m), dict):
                        o[m] = o.get(m, 0) + b[m]["n"]
                for st, k in b["states"].items():
                    o["states"][st] = o["states"].get(st, 0) + k
            for tier, t in (row.get("methylation_by_tier") or {}).items():
                o = out["methylation_by_tier"].setdefault(
                    tier, {"bp": 0, "calls": 0, "covered": 0, "_f": 0.0}
                )
                o["bp"] += t["bp"]
                o["calls"] += t["calls"]
                o["covered"] += t["covered"]
                o["_f"] += (t["fraction"] or 0.0) * t["covered"]
    for out in cells.values():
        for group in (*out["promoters"].values(), *out["registry_classes"].values()):
            for m in ep.MARKS:
                if m in group:
                    group[m] = {
                        "n": group[m],
                        "share": round(group[m] / group["n"], 4) if group["n"] else None,
                    }
        for t in out["methylation_by_tier"].values():
            f = t.pop("_f")
            t["fraction"] = round(f / t["covered"], 4) if t["covered"] else None
            t["covered_share_of_calls"] = round(t["covered"] / t["calls"], 4) if t["calls"] else None
    sig = ep.signal_coverage()
    signal_chroms = sorted(
        {c for marks in sig.values() for cs in marks.values() for c in cs}, key=lambda c: ep.CHROMS.index(c)
    )
    complete = [c for c in signal_chroms if all(c in cs for marks in sig.values() for cs in marks.values())]
    save_result(
        "epigenome_genome_wide",
        {
            "chromosomes": chroms,
            "cell_types": cells,
            "fold_change_profiles": {
                "chromosomes_all_cells_all_marks": complete,
                "note": "fold change over control is read on these chromosomes only; every other chromosome "
                "carries replicated peaks for all five marks and WGBS where it exists, and a record's "
                "fold_change field is UNKNOWN there with that reason",
            },
        },
    )
    for cell, out in cells.items():
        rd, sl = out["promoters"].get("read", {}), out["promoters"].get("silent", {})
        k4 = (rd.get("H3K4me3") or {}).get("share")
        k27 = (sl.get("H3K27me3") or {}).get("share")
        fos = (out["methylation_by_tier"].get("fossil") or {}).get("fraction")
        print(
            f"{cell:24s} read {rd.get('n')} H3K4me3 {k4}; silent {sl.get('n')} H3K27me3 {k27}; "
            f"fossil methylation {fos}"
        )


# ------------------------------------------------------------------------------------------
# The Alu lead: is the methylation Alu remains keep in hypomethylated lines a family effect?
# ------------------------------------------------------------------------------------------

HYPOMETHYLATED = ("K562", "HepG2", "GM12878", "SK-N-SH")
REGION = 50_000


def _element_methylation(cols, start: int, end: int) -> tuple[float, int] | None:
    """Mean fraction over the bins at least half inside [start, end), and their covered calls."""
    lo = -(-start // ep.BIN) if start % ep.BIN > ep.BIN // 2 else start // ep.BIN
    hi = (end // ep.BIN) + (1 if end % ep.BIN >= ep.BIN // 2 else 0)
    covered = sum(cols[1][lo:hi])
    if covered < 2:
        return None
    return sum(cols[2][lo:hi]) / covered / 1000, covered


def _wcgw_share(seq: str) -> tuple[float, int]:
    """Share of an element's CpGs flanked by A or T on both sides (solo-WCGW, Zhou et al. 2018)."""
    s = seq.upper()
    n = w = 0
    for i in range(1, len(s) - 2):
        if s[i] == "C" and s[i + 1] == "G":
            n += 1
            w += s[i - 1] in "AT" and s[i + 2] in "AT"
    return (w / n if n else math.nan), n


def alu(argv: list[str]) -> None:
    from genomeos.coords import Locus
    from genomeos.genome import IndexedGenome
    from genomeos.results import load_result, save_result

    m = ep.load_manifest()
    cells = [c for c, r in m["cell_types"].items() if "cpg" in r["methylation"]]
    chroms = argv or [c for c in ep.CHROMS if load_result(f"budget_{c}")]
    units: list[dict] = []
    for chrom in chroms:
        rp = Path(f"data/results/rmsk_{chrom}.bed.gz")
        budget = load_result(f"budget_{chrom}")
        if not rp.exists() or not budget:
            continue
        blocks = sorted((b["start"], b["end"]) for b in budget["blocks"] if b["guess"]["tier"] == "fossil")
        starts = [b[0] for b in blocks]
        cols = {c: ep.load_methylation_profile(c, chrom) for c in cells}
        if any(v is None for v in cols.values()):
            continue
        g = IndexedGenome(f"data/reference/{chrom}.fa.gz")
        import bisect

        with gzip.open(rp, "rt") as fh:
            for line in fh:
                if line[0] == "#":
                    continue
                f = line.rstrip("\n").split("\t")
                s0, e0, cls, fam = int(f[0]), int(f[1]), f[2], f[3]
                family = (
                    "Alu" if (cls, fam) == ("SINE", "Alu") else "L1" if (cls, fam) == ("LINE", "L1") else None
                )
                if family is None or e0 - s0 < 200:
                    continue
                i = bisect.bisect_right(starts, s0) - 1
                if i < 0 or not (blocks[i][0] <= s0 and e0 <= blocks[i][1]):
                    continue
                seq = str(g.fetch(Locus(chrom, s0, e0)))
                gc, cpg = ep.gc_cpg(seq)
                wcgw, n_cpg = _wcgw_share(seq)
                if gc != gc or n_cpg < 2:
                    continue
                meth = {}
                for c in cells:
                    em = _element_methylation(cols[c], s0, e0)
                    if em:
                        meth[c] = em[0]
                if not meth:
                    continue
                units.append(
                    {
                        "family": family,
                        "region": f"{chrom}:{s0 // REGION}",
                        "region_mb": f"{chrom}:{s0 // 1_000_000}",
                        "gc": gc,
                        "cpg": cpg,
                        "div": float(f[5]),
                        "wcgw": wcgw,
                        "meth": meth,
                    }
                )
        g.close()
        print(f"{chrom}: {len(units):,} Alu and L1 remains in fossil blocks so far", flush=True)

    def div_bin(d: float) -> int:
        return min(int(d / 0.05), 6)

    def wcgw_bin(w: float) -> int:
        return min(int(w / 0.2), 4)

    def matched_difference(cell: str, key) -> tuple[float | None, int, int]:
        """Alu minus L1 inside strata given by `key`, weighted by Alu elements; strata need both."""
        strata: dict = {}
        for u in units:
            if cell not in u["meth"]:
                continue
            st = strata.setdefault(key(u), {"Alu": [0.0, 0], "L1": [0.0, 0]})
            st[u["family"]][0] += u["meth"][cell]
            st[u["family"]][1] += 1
        num = den = 0.0
        used_alu = used_l1 = 0
        for st in strata.values():
            (sa, na), (sl, nl) = st["Alu"], st["L1"]
            if na < 3 or nl < 3:
                continue
            num += na * (sa / na - sl / nl)
            den += na
            used_alu += na
            used_l1 += nl
        return (round(num / den, 4) if den else None), used_alu, used_l1

    def base(u):
        return (ep.stratum(u["gc"], u["cpg"]), div_bin(u["div"]))

    designs = {
        "raw": lambda u: 0,
        "gc_cpg": lambda u: ep.stratum(u["gc"], u["cpg"]),
        "gc_cpg_age": base,
        "gc_cpg_age_wcgw": lambda u: (*base(u), wcgw_bin(u["wcgw"])),
        "same_50kb_region": lambda u: u["region"],
        "same_1mb_gc_cpg_wcgw": lambda u: (
            u["region_mb"],
            ep.stratum(u["gc"], u["cpg"]),
            wcgw_bin(u["wcgw"]),
        ),
    }
    out_cells = {}
    for c in cells:
        row = {}
        for name, key in designs.items():
            d, na, nl = matched_difference(c, key)
            row[name] = {"alu_minus_l1": d, "alu": na, "l1": nl}
        out_cells[c] = row
        print(c, {k: v["alu_minus_l1"] for k, v in row.items()}, flush=True)
    # context of the families: how different are they on the matching variables
    fam_desc = {}
    for fam in ("Alu", "L1"):
        us = [u for u in units if u["family"] == fam]
        if not us:
            continue
        fam_desc[fam] = {
            "elements": len(us),
            "median_gc": round(sorted(u["gc"] for u in us)[len(us) // 2], 3),
            "median_cpg_per_100bp": round(sorted(u["cpg"] for u in us)[len(us) // 2], 2),
            "median_divergence": round(sorted(u["div"] for u in us)[len(us) // 2], 3),
            "median_wcgw_share": round(sorted(u["wcgw"] for u in us)[len(us) // 2], 3),
        }
    # the WCGW gradient itself, in every family: methylation by solo-WCGW share
    gradient = {}
    for c in cells:
        gradient[c] = {}
        for fam in ("Alu", "L1"):
            bins: dict[int, list[float]] = {}
            for u in units:
                if u["family"] == fam and c in u["meth"]:
                    bins.setdefault(wcgw_bin(u["wcgw"]), []).append(u["meth"][c])
            gradient[c][fam] = {
                f"wcgw_{b * 20}_{b * 20 + 20}pct": {"n": len(v), "mean": round(sum(v) / len(v), 4)}
                for b, v in sorted(bins.items())
                if len(v) >= 20
            }
    out = {
        "chromosomes": chroms,
        "elements": len(units),
        "families": fam_desc,
        "cell_types": out_cells,
        "hypomethylated_lines": list(HYPOMETHYLATED),
        "wcgw_gradient": gradient,
        "designs": {
            "raw": "Alu mean minus L1 mean",
            "gc_cpg": "inside GC (0.05) x CpG per 100 bp strata",
            "gc_cpg_age": "and RepeatMasker divergence in 0.05 bins",
            "gc_cpg_age_wcgw": "and the share of solo-WCGW CpGs in 0.2 bins",
            "same_50kb_region": "inside the same 50 kb window (regional domain held)",
            "same_1mb_gc_cpg_wcgw": "inside the same 1 Mb window, GC x CpG stratum and solo-WCGW share",
        },
        "note": "RepeatMasker Alu (SINE/Alu) and L1 (LINE/L1) records of 200+ bp wholly inside fossil-tier "
        "blocks; methylation from 200 bp bins at least half inside the element with 2+ calls of "
        f"{ep.MIN_COVERAGE}+ reads, so a flank can dilute an element's value toward its neighbourhood"
        f"; differences are weighted by Alu "
        "elements in strata with 3+ of each family",
        "evidence": {"methylation": ep.EVIDENCE_METHYLATION, "repeats": "curated: UCSC RepeatMasker (rmsk)"},
    }
    save_result("epigenome_fossil_alu", out)
    print("saved", flush=True)


# ------------------------------------------------------------------------------------------
# Does calling poised promoters make the reader more right? Measured RNA says
# ------------------------------------------------------------------------------------------

RNA_CELLS = ("K562", "HepG2", "GM12878", "IMR-90")
EXPRESSED = 0.3  # mean covered exon fraction the segment filter uses (genome/rna_measured.py)


def _direct_url(href: str) -> str:
    """The S3 URL behind an ENCODE download href, so range requests skip the portal's redirect."""
    acc = href.rstrip("/").split("/")[-1].split(".")[0]
    meta = ep._portal(f"/files/{acc}/?format=json&frame=object")
    return (meta.get("cloud_metadata") or {}).get("url") or href


def _reader_check_chrom(chrom: str) -> tuple[str, dict, dict]:
    import bisect

    from genomeos.genome import rna_measured
    from genomeos.results import load_result

    tallies: dict = {c: {} for c in RNA_CELLS}
    orientation: dict = {}
    try:
        ann = _annotation(chrom)
    except SystemExit:
        return chrom, tallies, orientation
    genes = [g for g in ann.genes.values() if g.locus.chrom == chrom and g.type == "protein_coding"]
    by_strand: dict[str, list[tuple[int, int]]] = {"+": [], "-": []}
    gene_exons: dict[str, tuple[str, list[tuple[int, int]]]] = {}
    for g in genes:
        ex = sorted({(e.start, e.end) for t in g.transcripts.values() for e in t.exons})
        if ex:
            gene_exons[g.symbol] = (g.locus.strand.value, ex)
            by_strand[g.locus.strand.value].extend(ex)
    merged: dict[str, list[tuple[int, int]]] = {}
    for st, ivs in by_strand.items():
        out: list[list[int]] = []
        for a, b in sorted(ivs):
            if out and a <= out[-1][1]:
                out[-1][1] = max(out[-1][1], b)
            else:
                out.append([a, b])
        merged[st] = [(a, b) for a, b in out]
    for cell in RNA_CELLS:
        r = load_result(f"reader_{reader.slug(cell)}_{chrom}")
        if not r or "poised_genes" not in r:
            continue
        rna = rna_measured.MeasuredRna(cell, chrom)
        for tr in rna.tracks_by_cell.values():
            for t in tr.get("tracks", {}).values():
                t["href"] = _direct_url(t["href"])
        rna.prepare(merged)
        orientation[cell] = rna.orientation.get(cell)
        poised = set(r["poised_genes"])
        silent = set(r["silent_genes"]) - poised
        marked = set(r.get("read_by_marks", []))
        for sym, (st, ex) in gene_exons.items():
            starts = [a for a, _ in merged[st]]
            fr = []
            for a, _b in ex:
                i = bisect.bisect_right(starts, a) - 1
                if i >= 0:
                    m0, m1 = merged[st][i]
                    fr.append(rna.covered_fraction(st, m0, m1))
            score = sum(fr) / len(fr) if fr else 0.0
            call = (
                "poised"
                if sym in poised
                else "closed"
                if sym in silent
                else "read_by_marks"
                if sym in marked
                else "read_open"
            )
            groups = [call, *(["read"] if call.startswith("read") else [])]
            for grp in groups:
                t = tallies[cell].setdefault(grp, {"genes": 0, "expressed": 0})
                t["genes"] += 1
                t["expressed"] += score >= EXPRESSED
    return chrom, tallies, orientation


def reader_check(argv: list[str]) -> None:
    from genomeos.genome import rna_measured
    from genomeos.results import save_result

    chroms = argv or [f"chr{i}" for i in range(1, 23)] + ["chrX"]
    tallies: dict = {c: {} for c in RNA_CELLS}
    orientation: dict = {}
    with ProcessPoolExecutor(6) as pool:
        for chrom, t_chrom, orient in pool.map(_reader_check_chrom, chroms):
            for cell, groups in t_chrom.items():
                for grp, v in groups.items():
                    t = tallies[cell].setdefault(grp, {"genes": 0, "expressed": 0})
                    t["genes"] += v["genes"]
                    t["expressed"] += v["expressed"]
            for cell, o in orient.items():
                orientation.setdefault(cell, {})[chrom] = o
            print(chrom, {c: {k: v["genes"] for k, v in t.items()} for c, t in t_chrom.items()}, flush=True)
    for t in tallies.values():
        for v in t.values():
            v["expressed_share"] = round(v["expressed"] / v["genes"], 4) if v["genes"] else None
        if "read_open" in t and "poised" in t:
            open_genes = t["read_open"]["genes"] + t["poised"]["genes"]
            open_expr = t["read_open"]["expressed"] + t["poised"]["expressed"]
            t["read_by_openness_only"] = {
                "genes": open_genes,
                "expressed": open_expr,
                "expressed_share": round(open_expr / open_genes, 4),
            }
    save_result(
        "reader_poised_check",
        {
            "chromosomes": chroms,
            "cell_types": tallies,
            "orientation": orientation,
            "expressed": f"mean covered fraction of the gene's exons at or above {EXPRESSED}, ENCODE total "
            f"RNA-seq on the gene's strand (signal {rna_measured.SIGNAL})",
            "evidence": {
                "rna": rna_measured.EVIDENCE,
                "calls": "reader_<cell>_<chrom> after the poised rule",
            },
        },
    )
    for c, t in tallies.items():
        print(c, {k: (v["genes"], v["expressed_share"]) for k, v in t.items()}, flush=True)


# ------------------------------------------------------------------------------------------
# Question 1 again, on held-out chromosomes: does the marks' gain transfer off chr21 and chr22,
# and does any of it survive beside the model's own predicted direction in the other three lines?
#
# The chr21/chr22 comparison was fitted and scored on the same two chromosomes by cross-validation.
# Here the model is fitted on chr21 and chr22 only and scored on chromosomes it never saw, and the
# marks compete with two things at once: their own controls (the same twelve columns taken from
# another of the four lines, and shuffled within chromosome x cell x class x DNase call) and the
# outcome's sibling lines, which a separate lane measured at 0.865 AUC for direction.
# ------------------------------------------------------------------------------------------

TRAIN = ("chr21", "chr22")  # the two chromosomes the first comparison was fitted and scored on
MARK_COLS = (
    [f"{m}_peak" for m in ep.MARKS] + [f"{m}_fc" for m in ep.MARKS] + ["methylation", "methylation_missing"]
)
BASE_COLS = ["pels", "ctcf", "dnase_peak", "dnase_signal"] + [f"cell_{c}" for c in DIRECTION_CELLS[1:]]
SIB_COLS = ["sib_mean", "sib_net_sign", "sib_abs_mean"]
OTHER_COLS = [f"other_{c}" for c in MARK_COLS]  # control: the marks of another of the four lines
PERM_COLS = [f"shuffled_{c}" for c in MARK_COLS]  # control: marks shuffled within a stratum
ALL_COLS = BASE_COLS + SIB_COLS + MARK_COLS + OTHER_COLS + PERM_COLS
#: an all-element archive counts as held out only when the chain has scored at least this share of
#: the chromosome's cCREs and nothing has written to it lately: a sweep in progress is not a
#: held-out chromosome, it is the first per cent of one.
MIN_ARCHIVE_SHARE = 0.5
ARCHIVE_QUIET_SECONDS = 600
#: chromosomes kept out on purpose, with the reason, rather than by whether a file happens to exist
EXCLUDE = {
    "chrY": "two of the four lines (GM12878, IMR-90) are female, so their chrY marks and predicted "
    "effects are not about a chromosome the cell has"
}


def _profiles_complete(chrom: str) -> bool:
    return all(ep.signal_path(c, m, chrom).exists() for c in DIRECTION_CELLS for m in ep.MARKS) and all(
        ep.methylation_path(c, chrom).exists() for c in DIRECTION_CELLS
    )


def _archive_share(chrom: str) -> float:
    """How much of the chromosome's cCRE set the all-enhancer chain has scored per cell line."""
    p = Path(f"data/knowledge/alphagenome/all_elements/{chrom}.json")
    with gzip.open(f"data/results/ccres_{chrom}.bed.gz", "rt") as fh:
        total = sum(1 for line in fh if line[0] != "#")
    scored = sum(1 for e in json.loads(p.read_text()) if e.get("predicted_by_cell"))
    return scored / total if total else 0.0


def _held_out_reason(chrom: str, share: float) -> str | None:
    """Why a chromosome cannot serve as held out, or None when it can."""
    if chrom in EXCLUDE:
        return EXCLUDE[chrom]
    if not _profiles_complete(chrom):
        return "no fold-change profile for one of the four lines"
    p = Path(f"data/knowledge/alphagenome/all_elements/{chrom}.json")
    if time.time() - p.stat().st_mtime < ARCHIVE_QUIET_SECONDS:
        return "the all-element archive is being written right now by the chromosome sweep"
    if share < MIN_ARCHIVE_SHARE:
        return f"the all-element archive covers {share:.1%} of the chromosome's cCREs, the sweep is partway"
    return None


def _chrom_units(np, chrom: str) -> tuple:
    """One chromosome's element-line units as a float32 matrix over ALL_COLS, with the outcome and
    an element id per row (the rows of one element resample together). A missing mark stays NaN and
    is filled from the training median later, so nothing on a held-out chromosome informs its own
    imputation."""
    layer = ep.Layer(chrom, DIRECTION_CELLS)
    rot = dict(zip(DIRECTION_CELLS, DIRECTION_CELLS[1:] + DIRECTION_CELLS[:1], strict=True))
    rows: list[list[float]] = []
    lfcs: list[float] = []
    groups: list[int] = []
    strata: list[tuple] = []
    blank = [math.nan] * len(MARK_COLS)
    for e in _element_rows(chrom):
        by = e["by"]
        cells = [c for c in DIRECTION_CELLS if c in by]
        if len(cells) < 2:
            continue  # no sibling line, so no fair comparison to make
        feats = {c: _features(layer, c, e["start"], e["end"]) for c in cells}
        marks = {}
        for c, f in feats.items():
            f["methylation_missing"] = 1.0 if f["methylation"] is None else 0.0
            marks[c] = [math.nan if f[k] is None else float(f[k]) for k in MARK_COLS]
        gid = int(hashlib.md5(e["id"].encode()).hexdigest()[:8], 16)  # noqa: S324 - resample key only
        pels = 1.0 if e["cls"] == "pELS" else 0.0
        ctcf = 1.0 if e["ctcf"] else 0.0
        for cell in cells:
            f = feats[cell]
            others = [by[c] for c in cells if c != cell]
            rows.append(
                [pels, ctcf, f["dnase_peak"], f["dnase_signal"]]
                + [1.0 if cell == c else 0.0 for c in DIRECTION_CELLS[1:]]
                # the model's own behaviour in the sibling lines: the signed mean, the net sign
                # (genomeos-t1's two direction columns) and how large the effect is there
                + [
                    sum(others) / len(others),
                    sum((o > 0) - (o < 0) for o in others) / len(others),
                    sum(abs(o) for o in others) / len(others),
                ]
                + marks[cell]
                + marks.get(rot[cell], blank)
                + blank  # filled by the shuffle below
            )
            lfcs.append(by[cell])
            groups.append(gid)
            strata.append((cell, e["cls"], f["dnase_peak"]))
    if not rows:
        return np.zeros((0, len(ALL_COLS)), np.float32), np.zeros(0, np.float32), np.zeros(0, np.int64)
    x = np.array(rows, dtype=np.float32)
    own = len(BASE_COLS) + len(SIB_COLS)
    perm = len(ALL_COLS) - len(PERM_COLS)
    rng = random.Random(21)
    by_stratum: dict[tuple, list[int]] = {}
    for i, k in enumerate(strata):
        by_stratum.setdefault(k, []).append(i)
    for idx in by_stratum.values():
        donors = idx[:]
        rng.shuffle(donors)
        x[idx, perm:] = x[donors, own : own + len(MARK_COLS)]
    return x, np.array(lfcs, dtype=np.float32), np.array(groups, dtype=np.int64)


def _ties(np, v):
    """The sort order of `v` and a tie-group number per sorted position. The rank machinery below
    is weighted, so ties are grouped once rather than broken by position: the registry-only model
    gives thousands of rows exactly the same score, and breaking those by position inflates it."""
    order = np.argsort(v, kind="mergesort")
    s = v[order]
    new = np.empty(len(s), dtype=bool)
    if len(s):
        new[0] = True
        np.not_equal(s[1:], s[:-1], out=new[1:])
    return order, np.cumsum(new) - 1


def _auc_w(np, prep, y, w) -> float | None:
    """AUC with a weight per row, ties averaged. All weights 1 is the plain AUC; a bootstrap over
    elements is the multiplicity of each row's element, which then costs one pass, not a new sort."""
    order, gid = prep
    if not len(order):
        return None
    ys, ws = y[order].astype(bool), w[order]
    negw = np.where(~ys, ws, 0.0)
    posw = np.where(ys, ws, 0.0)
    p, n = posw.sum(), negw.sum()
    if not p or not n:
        return None
    grp = np.bincount(gid, weights=negw, minlength=int(gid[-1]) + 1)
    before = np.concatenate(([0.0], np.cumsum(grp)[:-1]))
    return float((posw * (before[gid] + 0.5 * grp[gid])).sum() / (p * n))


def _ranks_w(np, prep, w):
    """Weighted average ranks, in the original row order."""
    order, gid = prep
    grp = np.bincount(gid, weights=w[order], minlength=int(gid[-1]) + 1)
    before = np.concatenate(([0.0], np.cumsum(grp)[:-1]))
    r = np.empty(len(order), dtype=np.float64)
    r[order] = before[gid] + (grp[gid] + 1.0) / 2.0
    return r


def _spearman_w(np, prep_x, prep_y, w) -> float | None:
    if not len(w):
        return None
    rx, ry = _ranks_w(np, prep_x, w), _ranks_w(np, prep_y, w)
    tw = w.sum()
    mx, my = (w * rx).sum() / tw, (w * ry).sum() / tw
    sx = math.sqrt(float((w * (rx - mx) ** 2).sum() / tw))
    sy = math.sqrt(float((w * (ry - my) ** 2).sum() / tw))
    if not sx or not sy:
        return None
    return float((w * (rx - mx) * (ry - my)).sum() / tw / (sx * sy))


def _score(np, metric: str, prep_x, prep_y, y, w) -> float | None:
    return _auc_w(np, prep_x, y, w) if metric == "auc" else _spearman_w(np, prep_x, prep_y, w)


def direction_transfer(argv: list[str]) -> None:
    import numpy as np

    from genomeos.results import save_result

    candidates = argv or [
        c for c in ep.CHROMS if Path(f"data/knowledge/alphagenome/all_elements/{c}.json").exists()
    ]
    shares, left_out = {}, {}
    for c in candidates:
        shares[c] = round(_archive_share(c), 3)
        why = _held_out_reason(c, shares[c])
        if why:
            left_out[c] = why
    usable = [c for c in candidates if c not in left_out]
    test_chroms = [c for c in usable if c not in TRAIN]
    if not all(c in usable for c in TRAIN) or not test_chroms:
        raise SystemExit(f"need {TRAIN} and one held-out chromosome; usable {usable}, left out {left_out}")
    print(f"train {list(TRAIN)}, held out {test_chroms}", flush=True)
    for c, why in left_out.items():
        print(f"left out {c}: {why}", flush=True)
    parts = {}
    for c in TRAIN + tuple(test_chroms):
        parts[c] = _chrom_units(np, c)
        print(f"{c}: {len(parts[c][1]):,} element-line units", flush=True)
    xtr = np.vstack([parts[c][0] for c in TRAIN])
    ltr = np.concatenate([parts[c][1] for c in TRAIN])
    xte = np.vstack([parts[c][0] for c in test_chroms])
    lte = np.concatenate([parts[c][1] for c in test_chroms])
    gte = np.concatenate([parts[c][2] for c in test_chroms])
    cte = np.concatenate([np.full(len(parts[c][1]), i) for i, c in enumerate(test_chroms)])
    parts.clear()

    models = {
        "registry+dnase": BASE_COLS,
        "registry+dnase+marks": BASE_COLS + MARK_COLS,
        "registry+dnase+marks_other_line": BASE_COLS + OTHER_COLS,
        "registry+dnase+marks_shuffled": BASE_COLS + PERM_COLS,
        "sibling_direction": SIB_COLS + [f"cell_{c}" for c in DIRECTION_CELLS[1:]],
        "sibling+registry+dnase": SIB_COLS + BASE_COLS,
        "sibling+registry+dnase+marks": SIB_COLS + BASE_COLS + MARK_COLS,
        "sibling+registry+dnase+marks_other_line": SIB_COLS + BASE_COLS + OTHER_COLS,
    }
    comparisons = {
        "marks_over_registry_dnase": ("registry+dnase+marks", "registry+dnase"),
        "marks_over_other_line_marks": ("registry+dnase+marks", "registry+dnase+marks_other_line"),
        "marks_over_shuffled_marks": ("registry+dnase+marks", "registry+dnase+marks_shuffled"),
        "marks_over_sibling": ("sibling+registry+dnase+marks", "sibling+registry+dnase"),
        "marks_over_sibling_with_other_line_marks": (
            "sibling+registry+dnase+marks",
            "sibling+registry+dnase+marks_other_line",
        ),
        "sibling_over_the_marks_model": ("sibling+registry+dnase", "registry+dnase+marks"),
    }
    at = {c: i for i, c in enumerate(ALL_COLS)}

    def fit_predict(cols: list[str], keep_tr, keep_te, ytr):
        j = [at[c] for c in cols]
        a = xtr[np.ix_(keep_tr, j)].astype(np.float64)
        b = xte[np.ix_(keep_te, j)].astype(np.float64)
        med = np.nanmedian(a, axis=0)
        med = np.where(np.isnan(med), 0.0, med)
        a = np.where(np.isnan(a), med, a)
        b = np.where(np.isnan(b), med, b)
        mu, sd = a.mean(axis=0), a.std(axis=0)
        sd = np.where(sd > 0, sd, 1.0)
        a, b = (a - mu) / sd, (b - mu) / sd
        a = np.hstack([np.ones((len(a), 1)), a])
        b = np.hstack([np.ones((len(b), 1)), b])
        lam = np.eye(a.shape[1])
        lam[0, 0] = 0.0  # the intercept is not shrunk
        return b @ np.linalg.solve(a.T @ a + lam, a.T @ ytr)

    tasks = {
        "acts": ("auc", None, lambda lfc: (np.abs(lfc) >= ACTS).astype(np.float64)),
        "silencer_among_acting": ("auc", ACTS, lambda lfc: (lfc > 0).astype(np.float64)),
        "magnitude": ("spearman", None, lambda lfc: np.log10(np.abs(lfc) + 1e-3)),
    }
    wanted = sorted({m for pair in comparisons.values() for m in pair})
    results: dict = {}
    for task, (metric, floor, y_of) in tasks.items():
        keep_tr = np.ones(len(ltr), bool) if floor is None else np.abs(ltr) >= floor
        keep_te = np.ones(len(lte), bool) if floor is None else np.abs(lte) >= floor
        ytr, yte = y_of(ltr[keep_tr]).astype(np.float64), y_of(lte[keep_te]).astype(np.float64)
        chrom_of, groups = cte[keep_te], gte[keep_te]
        preds = {name: fit_predict(cols, keep_tr, keep_te, ytr) for name, cols in models.items()}
        prep = {name: _ties(np, p) for name, p in preds.items()}
        prep_y = _ties(np, yte) if metric == "spearman" else None
        ones = np.ones(len(yte))
        per_chrom = {}
        for i, c in enumerate(test_chroms):
            mask = chrom_of == i
            ys = yte[mask]
            py = _ties(np, ys) if metric == "spearman" else None
            row: dict = {"units": int(mask.sum())}
            for name in models:
                s = _score(np, metric, _ties(np, preds[name][mask]), py, ys, np.ones(len(ys)))
                row[name] = round(s, 4) if s is not None else None
            per_chrom[c] = row
            print(task, c, row, flush=True)
        pooled = {}
        for name in models:
            s = _score(np, metric, prep[name], prep_y, yte, ones)
            pooled[name] = round(s, 4) if s is not None else None
        # element bootstrap on the pooled held-out units: a resampled element carries all its rows
        uniq, inv = np.unique(groups, return_inverse=True)
        rng = np.random.default_rng(int(hashlib.md5(task.encode()).hexdigest()[:8], 16))  # noqa: S324
        diffs: dict[str, list[float]] = {k: [] for k in comparisons}
        for _ in range(200):
            counts = np.bincount(rng.integers(0, len(uniq), len(uniq)), minlength=len(uniq))
            w = counts.astype(np.float64)[inv]
            cache = {n: _score(np, metric, prep[n], prep_y, yte, w) for n in wanted}
            for k, (a, b) in comparisons.items():
                diffs[k].append(cache[a] - cache[b])
        results[task] = {
            "train_units": int(keep_tr.sum()),
            "held_out_units": int(keep_te.sum()),
            "per_chromosome": per_chrom,
            "pooled": pooled,
            "bootstrap_95_and_share_with_no_gain": {
                k: [
                    round(float(np.percentile(v, 2.5)), 4),
                    round(float(np.percentile(v, 97.5)), 4),
                    round(float(np.mean(np.array(v) <= 0)), 3),
                ]
                for k, v in diffs.items()
            },
        }
        print(task, "pooled", pooled, results[task]["bootstrap_95_and_share_with_no_gain"], flush=True)
    save_result(
        "epigenome_direction_transfer",
        {
            "train": list(TRAIN),
            "held_out": test_chroms,
            "left_out": left_out,
            "archive_share_of_ccres": shares,
            "cell_types": list(DIRECTION_CELLS),
            "models": models,
            "results": results,
            "sibling": "the same element's predicted signed log2 fold change averaged over the other "
            "three lines, the net sign over them (genomeos-t1's two direction columns) and the mean "
            "absolute change there; the line being predicted never enters its own features",
            "controls": "marks_other_line takes the same twelve columns from another of the four lines "
            "(the rotation K562 -> HepG2 -> GM12878 -> IMR-90), marks_shuffled shuffles them among units "
            "of the same chromosome, cell, registry class and DNase call; both keep every marginal, so a "
            "gain that survives neither is not the cell's own chromatin",
            "method": "ridge (lambda 1) on standardised features, fitted on chr21 and chr22 only and "
            "scored on chromosomes it never saw; the training median fills a missing mark on both sides; "
            "acts and direction by AUC with ties averaged, magnitude by Spearman with log10 |lfc|; 200 "
            "resamples of the held-out elements",
            "caveat": "every outcome and the sibling feature are AlphaGenome predictions, and AlphaGenome "
            "was trained on ENCODE histone ChIP-seq, DNase and RNA-seq of these four lines: a mark "
            "predicting its direction can be the model reading back a track it was trained on",
            "evidence": {
                "marks": ep.EVIDENCE_MARK,
                "methylation": ep.EVIDENCE_METHYLATION,
                "outcome": "predicted: AlphaGenome deletion scoring, chain tables, no new model call",
            },
        },
    )
    print("saved", flush=True)


def main() -> None:
    if len(sys.argv) < 2:
        print(__doc__)
        return
    action, rest = sys.argv[1], sys.argv[2:]
    if action == "fetch":
        fetch(rest)
    elif action == "direction":
        direction(rest)
    elif action == "hox":
        hox(rest)
    elif action == "fossil":
        fossil(rest)
    elif action == "summary":
        summary(rest)
    elif action == "genome":
        genome(rest)
    elif action == "alu":
        alu(rest)
    elif action == "reader-check":
        reader_check(rest)
    elif action == "direction-transfer":
        direction_transfer(rest)
    else:
        raise SystemExit(f"unknown action {action}")


if __name__ == "__main__":
    main()
