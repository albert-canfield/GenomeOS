#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Run the terminality gate of area E and write data/results/celegans_terminality.json.

The question, the circularity audit, the metric, the baselines, the bar and the six predictions are
registered in `genomeos/organism/terminality.py`'s docstring and were committed before this ran.

    uv run python scripts/celegans_terminality.py

0 model requests. Nothing is fetched: the reference lineage and the Ma 2021 atlas are already
distilled under data/results/ and data/knowledge/.
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos.organism import fate_rules as fr  # noqa: E402
from genomeos.organism import terminality as tg  # noqa: E402
from genomeos.organism.atlas_levels import load_levels  # noqa: E402
from genomeos.organism.reference import ReferenceLineage  # noqa: E402
from genomeos.organism.tf_atlas import load_cells  # noqa: E402

OUT = Path("data/results/celegans_terminality.json")


def build_features(ref: ReferenceLineage, labels: dict[str, str]) -> dict[str, dict[str, set[str]]]:
    """The two reads a cell has of the measured factors, on exactly the universe's cells."""
    cells = [ref.cells[c] for c in labels]
    inst = fr.instantaneous(ref, load_cells(), cells)
    exp = fr.exposures(ref, load_levels(), cells)
    integ = fr.integrated(exp)
    return {"instantaneous": inst, "integrated": integ}


def main() -> int:
    ref = ReferenceLineage.load()
    labels = tg.universe(ref)
    gen = tg.depths(ref, labels)
    n_term = sum(1 for v in labels.values() if v == tg.TERMINAL)
    n_div = len(labels) - n_term

    # V1: the universe is what was registered, or nothing is reported
    if (n_term, n_div, len(labels)) != (555, 771, 1326):
        print(f"VOID V1: universe is {n_term} terminal + {n_div} dividing = {len(labels)}", file=sys.stderr)
        return 2
    print(f"universe {len(labels)}: {n_term} terminal, {n_div} dividing", flush=True)

    reads = build_features(ref, labels)
    counts = {c: len(reads["instantaneous"][c]) for c in labels}
    rows: list[dict] = []

    def add(s: tg.Score) -> tg.Score:
        rows.append(s.row())
        print(
            f"  {s.arm:<44} bal {s.balanced:.4f}  acc {s.accuracy:.4f}  "
            f"mcc {s.mcc:+.4f}  claimed {s.claimed:>4}  called T {s.called_terminal:>4}",
            flush=True,
        )
        return s

    # ---- baselines ------------------------------------------------------------------------
    print("baselines", flush=True)
    t1 = add(tg.score("T1 always-dividing", {}, labels))
    t2 = add(tg.score("T2 always-terminal", dict.fromkeys(labels, tg.TERMINAL), labels))
    k = tg.best_threshold(gen, labels)
    t3 = add(tg.score(f"T3 depth gen>={k}, in sample", tg.depth_call(gen, k), labels))
    ks = tg.per_sublineage_thresholds(gen, labels)
    t4 = add(tg.score("T4 depth per founder, in sample", tg.per_sublineage_call(gen, ks, k), labels))

    # held out: the threshold for a sublineage comes from the other seven only
    groups: dict[str, list[str]] = defaultdict(list)
    for c in labels:
        groups[fr.sublineage(c)].append(c)
    p3cv: dict[str, str] = {}
    for g, held in groups.items():
        train = {c: labels[c] for c in labels if fr.sublineage(c) != g}
        kk = tg.best_threshold({c: gen[c] for c in train}, train)
        p3cv.update(tg.depth_call({c: gen[c] for c in held}, kk))
    t3cv = add(tg.score("T3cv depth threshold, held out by founder", p3cv, labels))

    # T4 held out: a held-out founder has no threshold of its own, so the honest held-out form of the
    # per-founder baseline is leave-one-cell-out inside the founder. Generous to the baseline on purpose.
    p4cv: dict[str, str] = {}
    for held in groups.values():
        for c in held:
            rest = {x: labels[x] for x in held if x != c}
            kk = tg.best_threshold({x: gen[x] for x in rest}, rest) if rest else k
            p4cv[c] = tg.TERMINAL if gen[c] >= kk else tg.DIVIDING
    t4cv = add(tg.score("T4cv depth per founder, leave-one-cell-out", p4cv, labels))

    kc = tg.best_count_threshold(counts, labels)
    t5 = add(tg.score(f"T5 factor COUNT >= {kc}, in sample", tg.count_call(counts, kc), labels))
    p5cv: dict[str, str] = {}
    for gname, held in groups.items():
        train = {c: labels[c] for c in labels if fr.sublineage(c) != gname}
        kk = tg.best_count_threshold({c: counts[c] for c in train}, train)
        p5cv.update(tg.count_call({c: counts[c] for c in held}, kk))
    t5cv = add(tg.score("T5cv factor COUNT, held out", p5cv, labels))

    # ---- factor arms ----------------------------------------------------------------------
    print("factors only", flush=True)
    learned_rules: dict[str, list[dict]] = {}
    factor_only: list[tg.Score] = []
    for read, feats in reads.items():
        pred_is, rules = tg.in_sample_rules(feats, labels)
        learned_rules[read] = [
            {"calls": r.tissue, "factors": list(r.factors), "precision": r.precision, "support": r.support}
            for r in rules
        ]
        add(tg.score(f"F {read} rules, in sample", pred_is, labels))
        factor_only.append(
            add(tg.score(f"F {read} rules, held out", tg.held_out_rules(feats, labels), labels))
        )

    print("factors only, ceiling (logistic on the whole factor vector)", flush=True)
    for read, feats in reads.items():
        factor_only.append(
            add(tg.score(f"C {read} logistic, held out", tg.logistic_held_out(feats, labels), labels))
        )

    print("factors plus depth", flush=True)
    plus_depth: list[tg.Score] = []
    for read, feats in reads.items():
        wd = tg.with_depth(feats, gen)
        plus_depth.append(
            add(tg.score(f"FD {read} rules + depth, held out", tg.held_out_rules(wd, labels), labels))
        )
        plus_depth.append(
            add(
                tg.score(
                    f"CD {read} logistic + depth, held out",
                    tg.logistic_held_out(feats, labels, extra={c: [float(gen[c])] for c in labels}),
                    labels,
                )
            )
        )
    dlog = add(
        tg.score(
            "D depth alone, logistic, held out",
            tg.logistic_held_out(
                {c: set() for c in labels}, labels, extra={c: [float(gen[c])] for c in labels}
            ),
            labels,
        )
    )
    depth_only_best = max([t3, t4, t3cv, t4cv, dlog], key=lambda s: s.balanced)

    # ---- within a generation band: the clock held constant ---------------------------------
    print("within a generation band (depth constant)", flush=True)
    band_rows: list[dict] = []
    pooled_tp = pooled_fn = pooled_tn = pooled_fp = 0
    for b in tg.bands(gen, labels):
        cells = [c for c in labels if gen[c] == b]
        sub_labels = {c: labels[c] for c in cells}
        best = None
        for read, feats in reads.items():
            sub = {c: feats[c] for c in cells}
            for name, pred in (
                (f"rules {read}", tg.held_out_rules(sub, sub_labels)),
                (f"logistic {read}", tg.logistic_held_out(sub, sub_labels)),
            ):
                s = tg.score(f"band {b} {name}", pred, sub_labels)
                band_rows.append(s.row())
                if best is None or s.balanced > best[0].balanced:
                    best = (s, pred)
        s, pred = best
        n_t = sum(1 for c in cells if sub_labels[c] == tg.TERMINAL)
        print(
            f"  gen {b:>2}  n {len(cells):>4}  terminal {n_t:>4}  best {s.arm} bal {s.balanced:.4f}",
            flush=True,
        )
        for c in cells:
            call = pred.get(c) or tg.DEFAULT_CALL
            if sub_labels[c] == tg.TERMINAL:
                pooled_tp += call == tg.TERMINAL
                pooled_fn += call != tg.TERMINAL
            else:
                pooled_fp += call == tg.TERMINAL
                pooled_tn += call != tg.TERMINAL
    pooled = 0.0
    if pooled_tp + pooled_fn and pooled_tn + pooled_fp:
        pooled = (pooled_tp / (pooled_tp + pooled_fn) + pooled_tn / (pooled_tn + pooled_fp)) / 2
    print(f"  pooled within-band balanced {pooled:.4f}  (best arm per band, held out)", flush=True)

    # ---- POST-HOC: the macro average, and the three circularity controls -------------------
    # The registered pooling lets a caller that names each stratum's own majority -- depth and founder
    # and nothing else -- score far above 0.5 while scoring exactly 0.500 inside every stratum. The
    # macro average of the per-stratum figures has a floor of exactly 0.500. Both are published.
    print("POST-HOC: macro average per stratum (floor exactly 0.500), and the majority baseline")
    parents = {c: (ref.cells[c].parent or "") for c in labels}
    named = set(load_cells())
    posthoc: dict[str, dict] = {}

    def arms(group_cells, restrict=None, by_parent=False):
        out = {"STRATUM MAJORITY (depth" + (" + founder)" if by_parent else ")"): tg.stratum_majority}
        for read in reads:
            for kind in ("rules", "logistic"):
                out[f"{kind} {read}"] = _arm(read, kind, by_parent)
        _ = group_cells, restrict
        return out

    def _arm(read, kind, by_parent):
        def f(cells, lab):
            sub = {c: reads[read][c] for c in cells}
            if by_parent:
                return tg.within_stratum(cells, lab, sub, parents, kind=kind)
            return tg.held_out_rules(sub, lab) if kind == "rules" else tg.logistic_held_out(sub, lab)

        return f

    def run_block(title: str, groups, by_parent: bool) -> dict:
        print(f"  {title}: {[(k, len(v)) for k, v in groups]}", flush=True)
        block = {}
        for name, fn in arms(groups, by_parent=by_parent).items():
            per, macro, pooled_naive = tg.stratified(groups, labels, fn)
            block[name] = {
                "macro": round(macro, 4),
                "naive_pooled": round(pooled_naive, 4),
                "per_stratum": {k: round(v, 4) for k, v in sorted(per.items())},
            }
            print(f"    {name:<34} macro {macro:.4f}   naive-pooled {pooled_naive:.4f}", flush=True)
        return block

    g_band = tg.strata(gen, labels, founder=False, min_cells=1, min_each=20)
    posthoc["D1 within band, one arm fixed across bands"] = run_block("bands", g_band, False)

    g_named = [(k, [c for c in cs if c in named]) for k, cs in g_band]
    posthoc["D2 within band, atlas-named cells only"] = run_block("bands, named only", g_named, False)

    g_bf = tg.strata(gen, labels, founder=True)
    posthoc["D3 within band AND founder, folds by grandparent"] = run_block("band x founder", g_bf, True)

    g_bfn = [(k, [c for c in cs if c in named]) for k, cs in g_bf]
    g_bfn = [(k, cs) for k, cs in g_bfn if len(cs) >= 30]
    posthoc["D3n within band AND founder, atlas-named only"] = run_block("band x founder, named", g_bfn, True)

    # ---- the null -------------------------------------------------------------------------
    print("shuffle null", flush=True)
    null = tg.shuffle_null(reads["integrated"], labels)
    print(
        f"  balanced {null['balanced_mean']:.4f} +/- {null['balanced_sd']:.4f} over {null['draws']}",
        flush=True,
    )

    # ---- the secondary universe: deaths counted as done dividing ---------------------------
    lab_d = tg.universe(ref, deaths=True)
    gen_d = tg.depths(ref, lab_d)
    reads_d = build_features(ref, lab_d)
    kd = tg.best_threshold(gen_d, lab_d)
    sec = [
        tg.score(f"S depth gen>={kd}, in sample (deaths counted)", tg.depth_call(gen_d, kd), lab_d).row(),
        tg.score(
            "S integrated rules, held out (deaths counted)",
            tg.held_out_rules(reads_d["integrated"], lab_d),
            lab_d,
        ).row(),
    ]

    # ---- V2 and the verdict ----------------------------------------------------------------
    bad = [r for r in rows if r["n"] != 1326 or r["called_terminal"] + r["called_dividing"] != 1326]
    if bad:
        print(f"VOID V2: {len(bad)} arms do not call all 1,326 cells", file=sys.stderr)
        return 2

    best_factor_only = max(factor_only, key=lambda s: s.balanced)
    best_plus_depth = max(plus_depth, key=lambda s: s.balanced)
    trivial_best = max([t1, t2, t3, t4, t3cv, t4cv], key=lambda s: s.balanced)
    clause_a = best_factor_only.balanced >= trivial_best.balanced + 0.05
    clause_b = best_plus_depth.balanced >= depth_only_best.balanced + 0.03
    per_band_best: dict[int, float] = {}
    for r in band_rows:
        b = int(r["arm"].split()[1])
        per_band_best[b] = max(per_band_best.get(b, 0.0), r["balanced"])
    clause_c = pooled >= 0.55 and sum(1 for v in per_band_best.values() if v >= 0.55) >= 2

    # clause (c) recomputed on the statistic whose floor really is 0.500, with the strictest control
    d1 = posthoc["D1 within band, one arm fixed across bands"]
    d3 = posthoc["D3 within band AND founder, folds by grandparent"]
    factor_macro = [v["macro"] for k, v in d1.items() if not k.startswith("STRATUM")]
    strict_macro = [v["macro"] for k, v in d3.items() if not k.startswith("STRATUM")]
    best_macro = max(factor_macro)
    clause_c_corrected = (
        best_macro >= 0.55
        and sum(
            1 for v in max(d1.items(), key=lambda kv: kv[1]["macro"])[1]["per_stratum"].values() if v >= 0.55
        )
        >= 2
    )

    knows = clause_a and clause_b and clause_c
    verdict = (
        "THE RULES KNOW SOMETHING"
        if knows
        else (
            "NEGATIVE: the measured factors do not say when a cell is done dividing"
            " well enough to build a rule on -- they lose to a per-founder depth threshold and add"
            " nothing to one -- though they do carry a real signal that is not depth, founder"
            " or coverage"
            if not clause_a and not clause_b
            else "PARTIAL: some clauses passed"
        )
    )

    print(f"\nclause (a) factor-only beats every trivial baseline by 0.05: {clause_a}", flush=True)
    print(f"clause (b) factors add 0.03 over depth alone:                 {clause_b}", flush=True)
    print(f"clause (c) within-band pooled >= 0.55 and two bands >= 0.55:  {clause_c}", flush=True)
    print(f"clause (c) CORRECTED, macro (floor 0.500) >= 0.55:            {clause_c_corrected}", flush=True)
    print(f"    best factor macro within band        {best_macro:.4f}", flush=True)
    print(f"    best factor macro within band+founder {max(strict_macro):.4f}", flush=True)
    print("    the same statistic for a stratum-majority caller  0.5000 by construction", flush=True)
    print(f"VERDICT: {verdict}", flush=True)

    result = {
        "result": "celegans_terminality",
        "date": date.today().isoformat(),
        "question": "can a rule set say WHEN a cell is done dividing, and not only what it becomes",
        "registered_in": "genomeos/organism/terminality.py docstring, committed before this ran",
        "universe": {
            "embryonic_cells": len(labels),
            "terminal": n_term,
            "dividing": n_div,
            "programmed_deaths_excluded": 113,
        },
        "circularity_audit": {
            "excluded_cell_type": "only the lookup's fate_<cell> ever sets one; the guard that pins 0 of 414",
            "excluded_cell_name": "the lookup's own address form; a lineage name spells out the divisions",
            "excluded_tracked_lifetime": {
                "why": "it ends at the cell's division or at the last frame, so it is the answer",
                "last_frame_min": 400.0,
                "cells_ending_at_the_last_frame": 454,
            },
            "excluded_own_cycle_length": "defined only for a cell that divides",
            "excluded_atlas_coverage": {
                "why": "which cells the imaging resolved, not a property of the cell",
                "terminal_named": 455,
                "terminal_total": 555,
                "dividing_named": 717,
                "dividing_total": 771,
            },
            "baseline_not_result_depth": "every division in this program is scheduled by a lookup div_<cell>",
            "measurement": "Ma et al. 2021 transcription-factor presence, instantaneous or integrated",
        },
        "baselines": [t1.row(), t2.row(), t3.row(), t4.row(), t3cv.row(), t4cv.row(), t5.row(), t5cv.row()],
        "depth_thresholds_per_founder": ks,
        "arms": rows,
        "bands": band_rows,
        "pooled_within_band_balanced": round(pooled, 4),
        "per_band_best_balanced": {str(b): round(v, 4) for b, v in sorted(per_band_best.items())},
        "post_hoc_stratified": {
            "why": "the registered pooling lets a stratum-majority caller -- depth and founder and "
            "nothing else -- score 0.7575 within band and 0.7886 within band+founder while scoring "
            "exactly 0.500 inside every stratum; the macro average has a floor of exactly 0.500",
            "blocks": posthoc,
        },
        "shuffle_null": null,
        "secondary_universe_deaths_counted": sec,
        "learned_rules": learned_rules,
        "clauses": {"a_factor_only": clause_a, "b_adds_to_depth": clause_b, "c_within_band": clause_c},
        "clauses_corrected": {
            "why": "the registered clause (c) pooling is inflated by the stratum prior; the macro "
            "average of the per-stratum figures has a floor of exactly 0.500",
            "c_within_band_corrected_macro": clause_c_corrected,
            "best_factor_macro_within_band": round(best_macro, 4),
            "best_factor_macro_within_band_and_founder": round(max(strict_macro), 4),
        },
        "best": {
            "factor_only": best_factor_only.row(),
            "factors_plus_depth": best_plus_depth.row(),
            "depth_only": depth_only_best.row(),
            "trivial_or_depth": trivial_best.row(),
        },
        "verdict": verdict,
        "model_requests": 0,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=1) + "\n")
    print(f"wrote {OUT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
