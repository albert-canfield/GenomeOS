#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Discovered element types against the registry's classes (area I).

    uv run python scripts/element_types.py
    uv run python scripts/element_types.py --fit chr21,chr22 --transfer chr19,chr20 --describe chr17

Clusters the elements of the completed chromosomes on their AlphaGenome deletion behaviour alone,
says what the clusters correspond to, runs the four controls, and scores discovered types against
the curated ones on the epigenome layer's held-out task. Writes `element_types`.

No model call: everything is read from the deletion archive and the committed results.
"""

from __future__ import annotations

import argparse
import json
import math
import random
import resource
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos.attribution import element_types as etypes  # noqa: E402
from genomeos.genome import epigenome as ep  # noqa: E402
from genomeos.results import save_result  # noqa: E402

CELLS = etypes.CELLS
MARK_PEAKS = [f"{m}_peak" for m in ep.MARKS]
MARK_FC = [f"{m}_fc" for m in ep.MARKS]
SEED = etypes.SEED
score = etypes.score


def peak_mem_mb() -> float:
    usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return round(usage / (1e6 if sys.platform == "darwin" else 1e3), 1)


def raw_profile(vectors: list[list[float]], labels: list[int], k: int) -> dict[str, dict[str, float]]:
    """Per cluster, the mean of every behaviour feature in its own units."""
    out: dict[str, dict[str, float]] = {}
    for j in range(k):
        rows = [v for v, lab in zip(vectors, labels, strict=True) if lab == j]
        if not rows:
            continue
        out[str(j)] = {
            name: round(sum(r[i] for r in rows) / len(rows), 4) for i, name in enumerate(etypes.FEATURES)
        }
    return out


def median(xs: list[float]) -> float | None:
    xs = sorted(x for x in xs if x is not None)
    return round(xs[len(xs) // 2], 4) if xs else None


def mean(xs: list[float]) -> float | None:
    xs = [x for x in xs if x is not None]
    return round(sum(xs) / len(xs), 4) if xs else None


def order_by_strength(fit: dict, k: int) -> dict:
    """Relabel clusters 0..k-1 by rising mean log strongest effect, so a table reads weakest first."""
    order = sorted(range(k), key=lambda j: fit["centroids"][j][0])
    new = {old: i for i, old in enumerate(order)}
    return {
        "centroids": [fit["centroids"][old] for old in order],
        "labels": [new[lab] for lab in fit["labels"]],
        "inertia": fit["inertia"],
    }


def tier_of(row: dict) -> str:
    return row["tier"] or "annotated (outside UNKNOWN blocks)"


def context(rows: list[dict], labels: list[int], k: int, units: list[dict] | None) -> dict:
    """What each cluster is, read from everything the vector was blind to."""
    out: dict = {}
    by_unit: dict[int, list[dict]] = {}
    for u in units or []:
        by_unit.setdefault(u["element"], []).append(u)
    for j in range(k):
        idx = [i for i, lab in enumerate(labels) if lab == j]
        members = [rows[i] for i in idx]
        if not members:
            continue
        us = [u for i in idx for u in by_unit.get(i, [])]
        states = Counter(u["state"] for u in us)
        out[str(j)] = {
            "elements": len(members),
            "share": round(len(members) / len(rows), 4),
            "registry_class": {
                c: round(n / len(members), 4) for c, n in Counter(r["cls"] for r in members).most_common()
            },
            "ctcf_bound": round(sum(r["ctcf"] for r in members) / len(members), 4),
            "budget_tier": {
                t: round(n / len(members), 4) for t, n in Counter(tier_of(r) for r in members).most_common()
            },
            "conserved_fraction_mean": mean([r.get("conserved_fraction") for r in members]),
            "share_with_conserved_bases": round(
                sum(1 for r in members if (r.get("conserved_fraction") or 0) > 0) / len(members), 4
            ),
            "length_median": median([r["length"] for r in members]),
            "gc_mean": mean([r["gc"] for r in members]),
            "nearest_gene_distance_median": median([r["nearest_any_gene_distance"] for r in members]),
            "dnase_open_share": round(sum(u["dnase_peak"] for u in us) / len(us), 4) if us else None,
            "h3k27me3_peak_share": round(sum(u["H3K27me3_peak"] for u in us) / len(us), 4) if us else None,
            "state": {s: round(n / len(us), 4) for s, n in states.most_common()} if us else None,
        }
    return out


def class_adjusted_states(rows: list[dict], labels: list[int], k: int, units: list[dict]) -> dict:
    """Per cluster, the measured-state share minus what its registry-class mix alone would give:
    a cluster the registry cannot see but the marks can shows large departures here."""
    by_class: dict[str, Counter] = {}
    for u in units:
        by_class.setdefault(rows[u["element"]]["cls"], Counter())[u["state"]] += 1
    class_share = {c: {s: n / sum(cnt.values()) for s, n in cnt.items()} for c, cnt in by_class.items()}
    out = {}
    for j in range(k):
        us = [u for u in units if labels[u["element"]] == j]
        if not us:
            continue
        observed = Counter(u["state"] for u in us)
        expected: Counter = Counter()
        for u in us:
            for s, p in class_share[rows[u["element"]]["cls"]].items():
                expected[s] += p
        dev = {s: round((observed[s] - expected[s]) / len(us), 4) for s in set(observed) | set(expected)}
        top = sorted(dev.items(), key=lambda kv: -abs(kv[1]))[:4]
        out[str(j)] = {
            "largest_departures": dict(top),
            "total_variation": round(sum(abs(v) for v in dev.values()) / 2, 4),
        }
    return out


def build_units(
    chroms: tuple[str, ...], rows_by: dict[str, list[dict]], offset_of: dict[str, int]
) -> list[dict]:
    """One unit per element and cell line, with the line's chromatin as the epigenome layer reads it."""
    units: list[dict] = []
    for chrom in chroms:
        layer = ep.Layer(chrom, CELLS)
        for i, r in enumerate(rows_by[chrom]):
            for cell in CELLS:
                if cell not in r["outcome"] or cell not in r["cells"]:
                    continue
                units.append(
                    {
                        "element": offset_of[chrom] + i,
                        "id": r["id"],
                        "chrom": chrom,
                        "cell": cell,
                        "lfc": r["outcome"][cell],
                        **etypes.chromatin(layer, cell, r["start"], r["end"]),
                    }
                )
        del layer
        print(f"  {chrom}: {len(units):,} element-cell units so far ({peak_mem_mb()} MB peak)", flush=True)
    return units


def features_of(
    u: dict, row: dict, model: str, label: int, k: int, vector: list[float], chroms: tuple[str, ...]
) -> list[float]:
    """The feature row of one element-cell unit under one model. `label` is the discovered type for
    this cell line (assigned from the other three lines), `vector` its leave-the-line-out behaviour."""
    cells = [1.0 if u["cell"] == c else 0.0 for c in CELLS[1:]]
    chrom = [1.0 if u["chrom"] == c else 0.0 for c in chroms[1:]]
    registry = [1.0 if row["cls"] == "pELS" else 0.0, 1.0 if row["ctcf"] else 0.0]
    dnase = [u["dnase_peak"], u["dnase_signal"]]
    peaks = [u[c] for c in MARK_PEAKS] + [u["methylation"], u["methylation_missing"]]
    full = peaks + [u[c] for c in MARK_FC] + [u[c + "_missing"] for c in MARK_FC]
    types = etypes.one_hot(label, k)
    return {
        "registry": cells + chrom + registry,
        "registry+dnase": cells + chrom + registry + dnase,
        "registry+dnase+marks": cells + chrom + registry + dnase + full,
        "registry+dnase+marks_peaks": cells + chrom + registry + dnase + peaks,
        "types": cells + chrom + types,
        "types+registry+dnase+marks_peaks": cells + chrom + registry + dnase + peaks + types,
        "types_shuffled": cells + chrom + types,
        "behaviour_vector": cells + chrom + vector,
        "other_line_direction": cells + chrom + [vector[1], vector[2]],
    }[model]


MODELS = (
    "registry",
    "registry+dnase",
    "registry+dnase+marks",
    "registry+dnase+marks_peaks",
    "types",
    "types+registry+dnase+marks_peaks",
    "types_shuffled",
    "behaviour_vector",
    "other_line_direction",
)

LABELS = {
    "registry": "curated: registry class (pELS/dELS, CTCF-bound)",
    "registry+dnase": "curated: class + DNase in the line (reader v1)",
    "registry+dnase+marks": "curated: class + DNase + the line's five marks (peaks and fold change) "
    "and methylation",
    "registry+dnase+marks_peaks": "curated: class + DNase + the line's mark peaks and methylation (no fold "
    "change: no signal profile is cached beyond chr21 and chr22)",
    "types": "discovered: the element's type, assigned from the other three lines' deletion behaviour",
    "types+registry+dnase+marks_peaks": "discovered type on top of the curated features (peaks)",
    "types_shuffled": "control: types taken from another element of the same class and DNase openness",
    "behaviour_vector": "ceiling: the leave-the-line-out behaviour vector itself, uncompressed",
    "other_line_direction": "two columns of that vector: the sign and the net direction of the deletion in "
    "the other three lines, which is what a type mostly carries",
}


def label_for(u: dict, loo: dict, donors: list[int] | None, model: str, key: str) -> int:
    element = u["element"]
    if model == "types_shuffled" and donors is not None:
        element = donors[element]
    return loo[u["cell"]][key][element]


def evaluate_cv(
    units: list[dict], rows: list[dict], loo: dict, donors: list[int], k: int, chroms: tuple[str, ...]
) -> dict:
    """The epigenome layer's own scoring: out-of-fold, five folds grouped by element."""
    import hashlib

    out: dict = {}
    groups = [int(hashlib.md5(u["id"].encode()).hexdigest()[:8], 16) for u in units]  # noqa: S324 - fold id
    for task, (y, idx) in etypes.tasks_of(units).items():
        out[task] = {"n": len(idx)}
        for model in MODELS:
            x = ep.standardise(
                [
                    features_of(
                        units[i],
                        rows[units[i]["element"]],
                        model,
                        label_for(units[i], loo, donors, model, "fit_labels"),
                        k,
                        loo[units[i]["cell"]]["fit_vectors"][units[i]["element"]],
                        chroms,
                    )
                    for i in idx
                ]
            )
            p = ep.cross_validated(x, y, [groups[i] for i in idx], folds=5, lam=1.0)
            out[task][model] = score(task, p, y)
            print(f"    cv {task:18s} {model:34s} {out[task][model]}", flush=True)
    return out


def evaluate_transfer(
    train_units: list[dict],
    train_rows: list[dict],
    test_units: list[dict],
    test_rows: list[dict],
    loo: dict,
    donors: list[int],
    test_donors: list[int],
    k: int,
) -> dict:
    """Fitted on the fit chromosomes and applied to held-out chromosomes: clusters, scales and ridge
    weights all come from the fit set. Chromosome indicators are dropped, since they cannot transfer."""
    one = ("all",)
    out: dict = {}
    preds: dict = {}
    train_tasks = etypes.tasks_of(train_units)
    test_tasks = etypes.tasks_of(test_units)
    for task in train_tasks:
        y_tr, idx_tr = train_tasks[task]
        y_te, idx_te = test_tasks[task]
        out[task] = {"n_train": len(idx_tr), "n_test": len(idx_te)}
        for model in MODELS:
            x_tr = [
                features_of(
                    train_units[i],
                    train_rows[train_units[i]["element"]],
                    model,
                    label_for(train_units[i], loo, donors, model, "fit_labels"),
                    k,
                    loo[train_units[i]["cell"]]["fit_vectors"][train_units[i]["element"]],
                    one,
                )
                for i in idx_tr
            ]
            x_te = [
                features_of(
                    test_units[i],
                    test_rows[test_units[i]["element"]],
                    model,
                    label_for(test_units[i], loo, test_donors, model, "transfer_labels"),
                    k,
                    loo[test_units[i]["cell"]]["transfer_vectors"][test_units[i]["element"]],
                    one,
                )
                for i in idx_te
            ]
            p = etypes.fit_predict(x_tr, y_tr, x_te)
            preds[(task, model)] = p
            out[task][model] = score(task, p, y_te)
            print(f"    transfer {task:18s} {model:34s} {out[task][model]}", flush=True)
        out[task]["bootstrap_95"] = bootstrap(
            task, y_te, [test_units[i]["element"] for i in idx_te], preds, PAIRS
        )
    return out


PAIRS = (
    ("types", "registry+dnase+marks_peaks"),
    ("types+registry+dnase+marks_peaks", "registry+dnase+marks_peaks"),
    ("types", "types_shuffled"),
    ("behaviour_vector", "registry+dnase+marks_peaks"),
    ("types", "other_line_direction"),
)


def bootstrap(
    task: str,
    y: list[float],
    elements: list[int],
    preds: dict,
    pairs: tuple[tuple[str, str], ...],
    draws: int = 200,
    seed: int = SEED,
) -> dict:
    """95% interval of each difference over element resamples of the scored set."""
    rng = random.Random(seed)
    by_element: dict[int, list[int]] = {}
    for j, e in enumerate(elements):
        by_element.setdefault(e, []).append(j)
    keys = list(by_element)
    diffs: dict[str, list[float]] = {f"{a} - {b}": [] for a, b in pairs}
    for _ in range(draws):
        sample = [j for e in (rng.choice(keys) for _ in keys) for j in by_element[e]]
        ys = [y[j] for j in sample]
        for a, b in pairs:
            sa = score(task, [preds[(task, a)][j] for j in sample], ys)
            sb = score(task, [preds[(task, b)][j] for j in sample], ys)
            if sa is not None and sb is not None:
                diffs[f"{a} - {b}"].append(sa - sb)
    out = {}
    for name, vals in diffs.items():
        vals.sort()
        if vals:
            lo, hi = vals[int(0.025 * len(vals))], vals[int(0.975 * len(vals)) - 1]
            out[name] = [round(lo, 4), round(hi, 4), round(sum(1 for v in vals if v <= 0) / len(vals), 3)]
    return out


def geometry_control(rows: list[dict], labels: list[int], vectors: list[list[float]], k: int) -> dict:
    """How much of the structure length, GC and distance to the nearest gene explain on their own."""
    geo = [
        [
            math.log10(max(r["length"], 1)),
            r["gc"] if r["gc"] is not None else 0.5,
            math.log10(1 + (r["nearest_any_gene_distance"] or 0)),
        ]
        for r in rows
    ]
    x = ep.standardise(geo)
    groups = list(range(len(rows)))
    membership = {}
    for j in range(k):
        y = [1.0 if lab == j else 0.0 for lab in labels]
        p = ep.cross_validated(x, y, groups, folds=5)
        membership[str(j)] = score("acts", p, y)
    explained = {}
    for i, name in enumerate(etypes.FEATURES):
        y = [v[i] for v in vectors]
        p = ep.cross_validated(x, y, groups, folds=5)
        m = sum(y) / len(y)
        sse = sum((a - b) ** 2 for a, b in zip(y, p, strict=True))
        sst = sum((a - m) ** 2 for a in y)
        explained[name] = round(1 - sse / sst, 4) if sst else None
    geo_scaled = etypes.scale(geo, *etypes.scaler(geo))
    geo_fit = etypes.kmeans(geo_scaled, k, SEED)
    return {
        "auc_of_cluster_membership_from_length_gc_distance": membership,
        "out_of_fold_r2_of_each_behaviour_feature": explained,
        "adjusted_rand_with_clusters_of_length_gc_distance": round(
            etypes.adjusted_rand(labels, geo_fit["labels"]), 4
        ),
        "reading": "AUC 0.5 and r2 0 would mean geometry explains nothing; 1 would mean the clusters are "
        "geometry renamed",
    }


def assign_chrom(
    rows: list[dict], means: list[float], sds: list[float], centroids: list[list[float]]
) -> list[int]:
    return etypes.assign(etypes.scale([etypes.behaviour(r) for r in rows], means, sds), centroids)


def shares(labels: list[int], k: int) -> dict[str, float]:
    n = len(labels) or 1
    counts = Counter(labels)
    return {str(j): round(counts.get(j, 0) / n, 4) for j in range(k)}


def leave_line_out(fit_rows: list[dict], tr_rows: list[dict], k: int, labels_fit: list[int]) -> dict:
    """For each line, k-means on the other three lines' vectors of the fit set, and the transfer set
    assigned to those centroids."""
    loo: dict = {}
    for cell in CELLS:
        others = tuple(c for c in CELLS if c != cell)
        v_fit = [etypes.behaviour(r, others) for r in fit_rows]
        m3, s3 = etypes.scaler(v_fit)
        f3 = etypes.kmeans(etypes.scale(v_fit, m3, s3), k, SEED)
        v_tr = [etypes.behaviour(r, others) for r in tr_rows]
        loo[cell] = {
            "fit_labels": f3["labels"],
            "fit_vectors": v_fit,
            "transfer_labels": etypes.assign(etypes.scale(v_tr, m3, s3), f3["centroids"]),
            "transfer_vectors": v_tr,
            "adjusted_rand_with_four_line_types": round(etypes.adjusted_rand(f3["labels"], labels_fit), 4),
        }
    return loo


def k_sweep(
    fit_rows: list[dict],
    tr_rows: list[dict],
    units_fit: list[dict],
    units_tr: list[dict],
    labels_fit: list[int],
) -> dict:
    """Sensitivity, not selection: the discovered types alone on the transfer task at every k, so the
    verdict does not hang on the k the criterion picked."""
    out: dict = {}
    train_tasks = etypes.tasks_of(units_fit)
    test_tasks = etypes.tasks_of(units_tr)
    for k in etypes.K_RANGE:
        loo = leave_line_out(fit_rows, tr_rows, k, labels_fit)
        out[str(k)] = {}
        for task in train_tasks:
            y_tr, idx_tr = train_tasks[task]
            y_te, idx_te = test_tasks[task]
            row = []
            for units, idx, key in ((units_fit, idx_tr, "fit_labels"), (units_tr, idx_te, "transfer_labels")):
                row.append(
                    [
                        [1.0 if units[i]["cell"] == c else 0.0 for c in CELLS[1:]]
                        + etypes.one_hot(loo[units[i]["cell"]][key][units[i]["element"]], k)
                        for i in idx
                    ]
                )
            p = etypes.fit_predict(row[0], y_tr, row[1])
            out[str(k)][task] = score(task, p, y_te)
        print(f"    k sweep {k}: {out[str(k)]}", flush=True)
    return out


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fit", default="chr21,chr22")
    ap.add_argument("--transfer", default="chr19,chr20")
    ap.add_argument("--describe", default="chr17,chr18,chrY")
    args = ap.parse_args(argv)
    fit_chroms = tuple(args.fit.split(","))
    transfer_chroms = tuple(args.transfer.split(","))
    describe_chroms = tuple(c for c in args.describe.split(",") if c)
    t0 = time.time()
    stages: dict[str, float] = {}

    def lap(name: str, since: float) -> float:
        stages[name] = round(time.time() - since, 1)
        print(
            f"[{name}] {stages[name]} s, {peak_mem_mb()} MB peak, {etypes.COST.bytes_read / 1e6:.0f} MB read"
        )
        return time.time()

    rows_by: dict[str, list[dict]] = {}
    t = time.time()
    for chrom in fit_chroms + transfer_chroms + describe_chroms:
        rows_by[chrom] = etypes.load_chromosome(chrom)
        print(f"  {chrom}: {len(rows_by[chrom]):,} elements ({peak_mem_mb()} MB peak)", flush=True)
    t = lap("load deletion archive, chain table, registry, GENCODE, sequence, constraint", t)

    def concat(chroms: tuple[str, ...]) -> tuple[list[dict], dict[str, int]]:
        rows: list[dict] = []
        offsets: dict[str, int] = {}
        for c in chroms:
            offsets[c] = len(rows)
            rows.extend(rows_by[c])
        return rows, offsets

    fit_rows, fit_off = concat(fit_chroms)
    tr_rows, tr_off = concat(transfer_chroms)

    # 1. behaviour vectors and the clusters, blind to class and marks
    vec_fit = [etypes.behaviour(r) for r in fit_rows]
    means, sds = etypes.scaler(vec_fit)
    z_fit = etypes.scale(vec_fit, means, sds)
    t = lap("behaviour vectors", t)
    sel = etypes.choose_k(z_fit)
    k = sel["k"]
    fit = order_by_strength(sel["fit"], k)
    labels_fit = fit["labels"]
    print(f"  k = {k}; curve {json.dumps(sel['curve'])}", flush=True)
    t = lap("choose k and cluster (k-means, silhouette over k = 2..10)", t)

    # the chromatin each element-cell unit carries (for the mapping and the held-out task)
    units_fit = build_units(fit_chroms, rows_by, fit_off)
    units_tr = build_units(transfer_chroms, rows_by, tr_off)
    for units in (units_fit, units_tr):
        etypes.impute(units, [*MARK_FC, "methylation"])
    t = lap("chromatin per element and cell line (epigenome layer)", t)

    # 2. what the clusters are
    unit_labels = [labels_fit[u["element"]] for u in units_fit]
    mapping = {
        "centroids_in_feature_units": raw_profile(vec_fit, labels_fit, k),
        "context": context(fit_rows, labels_fit, k, units_fit),
        "states_beyond_the_class_mix": class_adjusted_states(fit_rows, labels_fit, k, units_fit),
        "cramers_v": {
            "cluster_vs_registry_class": round(etypes.cramers_v(labels_fit, [r["cls"] for r in fit_rows]), 4),
            "cluster_vs_ctcf_bound": round(etypes.cramers_v(labels_fit, [r["ctcf"] for r in fit_rows]), 4),
            "cluster_vs_budget_tier": round(etypes.cramers_v(labels_fit, [tier_of(r) for r in fit_rows]), 4),
            "cluster_vs_state_per_line": round(
                etypes.cramers_v(unit_labels, [u["state"] for u in units_fit]), 4
            ),
            "registry_class_vs_state_per_line": round(
                etypes.cramers_v(
                    [fit_rows[u["element"]]["cls"] for u in units_fit], [u["state"] for u in units_fit]
                ),
                4,
            ),
            "cluster_vs_conserved_bases": round(
                etypes.cramers_v(labels_fit, [(r.get("conserved_fraction") or 0) > 0 for r in fit_rows]), 4
            ),
        },
    }
    t = lap("map clusters onto class, state, tier, constraint", t)

    # 3. controls
    col = etypes.column_shuffle(vec_fit, SEED)
    z_col = etypes.scale(col, *etypes.scaler(col))
    col_fit = etypes.kmeans(z_col, k, SEED)
    col_choice = etypes.choose_k(z_col)
    open_lines = [0] * len(fit_rows)
    for u in units_fit:
        open_lines[u["element"]] += int(u["dnase_peak"])
    donors = etypes.stratified_donors(
        [(r["cls"], r["ctcf"], open_lines[i]) for i, r in enumerate(fit_rows)], SEED
    )
    shuffled = [labels_fit[d] for d in donors]
    controls: dict = {
        "features_permuted_independently": {
            "silhouette_real": sel["curve"][k]["silhouette"],
            "silhouette_permuted_same_k": round(etypes.silhouette(z_col, col_fit["labels"]) or 0.0, 4),
            "inertia_real": round(fit["inertia"], 1),
            "inertia_permuted_same_k": round(col_fit["inertia"], 1),
            "k_the_criterion_picks_on_permuted": col_choice["k"],
            "silhouette_curve_permuted": {kk: c["silhouette"] for kk, c in col_choice["curve"].items()},
            "reading": "marginals kept, dependence between features destroyed; the real silhouette must "
            "clear the permuted one for the clusters to be more than the shape of each feature",
        },
    }
    shuffled_units = [shuffled[u["element"]] for u in units_fit]
    controls["vectors_shuffled_within_class_ctcf_and_open_lines"] = {
        "strata": len({(r["cls"], r["ctcf"], open_lines[i]) for i, r in enumerate(fit_rows)}),
        "cramers_v_real": mapping["cramers_v"],
        "cramers_v_shuffled": {
            "cluster_vs_registry_class": round(etypes.cramers_v(shuffled, [r["cls"] for r in fit_rows]), 4),
            "cluster_vs_budget_tier": round(etypes.cramers_v(shuffled, [tier_of(r) for r in fit_rows]), 4),
            "cluster_vs_state_per_line": round(
                etypes.cramers_v(shuffled_units, [u["state"] for u in units_fit]), 4
            ),
            "cluster_vs_conserved_bases": round(
                etypes.cramers_v(shuffled, [(r.get("conserved_fraction") or 0) > 0 for r in fit_rows]), 4
            ),
        },
        "states_beyond_the_class_mix_shuffled": class_adjusted_states(fit_rows, shuffled, k, units_fit),
        "reading": "whole vectors move between elements of one stratum, so the cloud (and its silhouette) is "
        "unchanged by construction; what the shuffle removes is which element carries which behaviour. "
        "Associations with class survive by design; associations with state, tier and constraint that "
        "survive are carried by class and openness, not by the behaviour",
    }
    controls["length_gc_nearest_gene_distance"] = geometry_control(fit_rows, labels_fit, vec_fit, k)
    t = lap("controls: permuted features, stratified shuffle, geometry", t)

    # transfer: the fit set's scale and centroids applied unchanged, against a refit on the new set
    vec_tr = [etypes.behaviour(r) for r in tr_rows]
    labels_tr = etypes.assign(etypes.scale(vec_tr, means, sds), fit["centroids"])
    refit = etypes.kmeans(etypes.scale(vec_tr, *etypes.scaler(vec_tr)), k, SEED)
    refit_choice = etypes.choose_k(etypes.scale(vec_tr, *etypes.scaler(vec_tr)))
    transfer: dict = {
        "fit": list(fit_chroms),
        "applied_to": list(transfer_chroms),
        "adjusted_rand_transferred_vs_refit": round(etypes.adjusted_rand(labels_tr, refit["labels"]), 4),
        "matched_agreement_transferred_vs_refit": round(
            etypes.matched_agreement(labels_tr, refit["labels"]), 4
        ),
        "k_the_criterion_picks_on_the_transfer_set": refit_choice["k"],
        "silhouette_of_transferred_labels": round(
            etypes.silhouette(etypes.scale(vec_tr, means, sds), labels_tr) or 0.0, 4
        ),
        "shares": {"fit": shares(labels_fit, k), "transferred": shares(labels_tr, k)},
        "shares_by_chromosome": {},
        "profile_transferred": raw_profile(vec_tr, labels_tr, k),
        "cramers_v_transferred": {
            "cluster_vs_registry_class": round(etypes.cramers_v(labels_tr, [r["cls"] for r in tr_rows]), 4),
            "cluster_vs_budget_tier": round(etypes.cramers_v(labels_tr, [tier_of(r) for r in tr_rows]), 4),
        },
    }
    for chrom in fit_chroms + transfer_chroms + describe_chroms:
        rows = rows_by[chrom]
        if rows:
            transfer["shares_by_chromosome"][chrom] = {
                "elements": len(rows),
                **shares(assign_chrom(rows, means, sds, fit["centroids"]), k),
            }
    t = lap("transfer to held-out chromosomes", t)
    # 4. the held-out task: types for each line from the other three lines only
    loo = leave_line_out(fit_rows, tr_rows, k, labels_fit)
    open_tr = [0] * len(tr_rows)
    for u in units_tr:
        open_tr[u["element"]] += int(u["dnase_peak"])
    donors_tr = etypes.stratified_donors(
        [(r["cls"], r["ctcf"], open_tr[i]) for i, r in enumerate(tr_rows)], SEED
    )
    t = lap("leave-the-line-out clusterings", t)
    print("  held-out task, cross-validated within the fit chromosomes", flush=True)
    within = evaluate_cv(units_fit, fit_rows, loo, donors, k, fit_chroms)
    t = lap("held-out task: five folds grouped by element (fit chromosomes)", t)
    print("  held-out task, fitted on the fit chromosomes and scored on the transfer chromosomes", flush=True)
    across = evaluate_transfer(units_fit, fit_rows, units_tr, tr_rows, loo, donors, donors_tr, k)
    t = lap("held-out task: fit chromosomes to transfer chromosomes", t)
    sweep = k_sweep(fit_rows, tr_rows, units_fit, units_tr, labels_fit)
    t = lap("held-out task: types alone at every k (sensitivity)", t)
    loaded = [r for c in rows_by for r in rows_by[c]]
    loaded_mb = sum(chrom_length(c) for c in rows_by) / 1e6
    genome_mb = 3_088.3  # GRCh38 primary assembly, chr1..chr22, X, Y
    scale_up = genome_mb / loaded_mb if loaded_mb else None
    classify_s = sum(v for name, v in stages.items() if name.startswith(("load", "behaviour", "choose")))
    model_s = sum(r.get("model_seconds") or 0.0 for r in loaded)
    efficiency = {
        "wall_clock_seconds_total": round(time.time() - t0, 1),
        "stages_seconds": stages,
        "peak_memory_mb": peak_mem_mb(),
        "bytes_read": etypes.COST.bytes_read,
        "files_read": etypes.COST.files,
        "elements": len(loaded),
        "element_cell_units": len(units_fit) + len(units_tr),
        "chromosomes_loaded_mb": round(loaded_mb, 1),
        "classification_seconds": round(classify_s, 1),
        "classification_seconds_per_1000_elements": round(1000 * classify_s / len(loaded), 2),
        "alphagenome_seconds_already_spent_on_these_elements": round(model_s),
        "genome_wide_estimate": {
            "method": "linear in sequence length from the chromosomes loaded here",
            "elements": round(len(loaded) * scale_up) if scale_up else None,
            "bytes_read_gb": round(etypes.COST.bytes_read * scale_up / 1e9, 1) if scale_up else None,
            "classification_minutes": round(classify_s * scale_up / 60, 1) if scale_up else None,
            "alphagenome_hours_needed_first": round(model_s * scale_up / 3600, 1) if scale_up else None,
            "peak_memory": "set by the largest chromosome's archive, parsed whole (chr1 would be the "
            "largest); not by the genome, since chromosomes are read one at a time",
        },
    }
    out = {
        "chromosomes": {
            "fit": list(fit_chroms),
            "transfer": list(transfer_chroms),
            "described": list(describe_chroms),
        },
        "cells": list(CELLS),
        "features": list(etypes.FEATURES),
        "k_criterion": f"highest mean silhouette on a fixed {etypes.SILHOUETTE_SAMPLE}-element subsample "
        f"(seed {SEED}) over k = {etypes.K_RANGE[0]}..{etypes.K_RANGE[-1]}, among the k whose smallest "
        f"cluster holds at least {etypes.MIN_CLUSTER_SHARE:.0%} of the elements; "
        "fixed before any result was read",
        "k": k,
        "silhouette_curve": sel["curve"],
        "elements_fit": len(fit_rows),
        "clusters": mapping,
        "controls": controls,
        "transfer": transfer,
        "held_out_task": {
            "outcome": "AlphaGenome predicted_by_cell: signed log2 fold change of the element's predicted "
            "target on each line's own track; positive = rise on deletion",
            "types_for_a_line": "k-means with the same k on the other three lines' behaviour vectors, fitted "
            "on the fit chromosomes; each element assigned to the nearest centroid",
            "leave_the_line_out_types_vs_four_line_types_adjusted_rand": {
                c: loo[c]["adjusted_rand_with_four_line_types"] for c in CELLS
            },
            "models": LABELS,
            "within_fit_chromosomes_cv": within,
            "fit_to_transfer_chromosomes": across,
            "fit_to_transfer_note": "no mark fold-change profile is cached for chr19 or chr20, so the "
            "registry+dnase+marks model meets constant columns there and is not a fair transfer; "
            "registry+dnase+marks_peaks is the curated comparison on held-out chromosomes",
            "types_alone_at_every_k_fit_to_transfer": sweep,
            "reference": "docs/NODES-READER-WRITER.md: registry 0.532, + DNase 0.594, + marks 0.616 for rise "
            "among acting on chr21 and chr22",
        },
        "caveat": "every behaviour value is AlphaGenome's prediction, and AlphaGenome was trained on ENCODE "
        "histone ChIP-seq, DNase and RNA-seq of these same lines: clusters that agree with ENCODE-derived "
        "marks partly read the model's inputs back, and types assigned from three lines that predict the "
        "fourth measure the model's consistency across lines rather than an element's biology",
        "evidence": "predicted: AlphaGenome deletion scoring (archive, no new model call); experimental: "
        "ENCODE marks, DNase, WGBS; curated: ENCODE cCREs v3, GENCODE 50; UCSC phastConsElements100way",
        "efficiency": efficiency,
    }
    name = (
        "element_types"
        if (fit_chroms, transfer_chroms) == (("chr21", "chr22"), ("chr19", "chr20"))
        else ("element_types_" + "_".join(fit_chroms + transfer_chroms))
    )
    print("saved", save_result(name, out), json.dumps(efficiency, indent=1), flush=True)


def chrom_length(chrom: str) -> int:
    fai = etypes.REFERENCE / f"{chrom}.fa.fai"
    return int(fai.read_text().split("\t")[1]) if fai.exists() else 0


if __name__ == "__main__":
    main()
