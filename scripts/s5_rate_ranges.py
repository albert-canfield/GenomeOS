# SPDX-License-Identifier: AGPL-3.0-or-later
"""Item 12 S5: run the simulable tier across its registered ranges and say which answers survive.

    uv run python scripts/s5_rate_ranges.py [--chroms chr21] [--sample 1000] [--no-save]

The registration is in 0ea0a90 (`attribution/bridge.py`, the `S5_*` constants, and the dated section
of docs/DESIGN-MINIMAL-CELL.md). In order, this script:

1. re-reads the range widths from the rate sources and stops if any differs from the registered
   number in the fourth decimal (the gate: the ranges are the registered ones, not re-chosen);
2. collects every measured-tier simulable gene-cell pair from the 24 compiled programs through
   `bridge.parameterize`, with its transfer label;
3. runs the registered sample (1,000 pairs, seed 20260928) through `NetworkRuntime` at every point
   of R0, R1, R2 and R3, one gene per program because the runtime holds one mRNA half-life per module,
   checks every trajectory's label, and runs the nominal point under all four activator rules;
4. computes the same predictions for every pair from the model's closed form (steady state =
   production / delta_m; half-way time after removal = t_half), and reports that census only if it
   agrees with the runtime on every sampled pair;
5. classifies each prediction as stable or unstable over the primary range (R0 + R1 + R2) and over
   R3, counts within-context rank comparisons, and counts the one human context's matched pairs.

Stable across these ranges is not validated: `bridge.SIMULATION_STATUS` travels on every record.
Reads the compiled programs, data/results/gene_rates.json and data/cache/rates/; no network, no
AlphaGenome request.
"""

from __future__ import annotations

import argparse
import bisect
import json
import math
import random
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from genomeos import manifest as mf
from genomeos.attribution import bridge
from genomeos.evidence import COMPILED_DIR
from genomeos.lang.parser import parse
from genomeos.results import RESULTS_DIR, save_result
from genomeos.runtime import grn

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bridge_audit  # noqa: E402  (the audit's rate loaders, one reading of the table)

CHROMS = bridge_audit.CHROMS
RATE_TABLE = RESULTS_DIR / "gene_rates.json"
CACHE = Path("data/cache/rates")
SCHWAN = CACHE / "schwanhausser2011_TableS3.xls"
SCHOFIELD = CACHE / "schofield2018_TableS2_halflives.xlsx"
LN2 = math.log(2)
ONE_MOLECULE = 1.0  # the on/off threshold registered in S5_PREDICTIONS
TWO_FOLD = 2.0  # the level / response-time criterion registered in S5_PREDICTIONS
GATE_DECIMALS = 4
REPLICATES = {
    "w_T": (
        "transcription rate (vsr) experiment [molecules/(cell*h)]",
        "transcription rate (vsr) replicate [molecules/(cell*h)]",
    ),
    "w_t_half": ("mRNA half-life experiment [h]", "mRNA half-life replicate [h]"),
}


# ---- 1. the ranges, re-read and gated against the registration --------------------------------------


def read_ranges() -> dict:
    """The S5_RANGES statistics recomputed from the sources, exactly as registered."""
    out: dict = {}
    d = pd.read_excel(SCHWAN)
    for key, (e, r) in REPLICATES.items():
        x = d[[e, r]].dropna()
        x = x[(x[e] > 0) & (x[r] > 0)]
        out[key] = float(np.percentile(np.abs(np.log2(x[e] / x[r])), 95))
        out[key + "_rows"] = len(x)
    table = json.loads(RATE_TABLE.read_text())["rates"]
    k562 = pd.read_excel(SCHOFIELD, sheet_name="Table S2_K562")
    human = dict(zip(k562["transcript"].astype(str), k562["mean_half_life"].astype(float), strict=True))
    logs = [
        math.log2(human[g] / r["mrna_half_life_h"])
        for g, r in table.items()
        if g in human and r.get("mrna_half_life_h") and human[g] > 0
    ]
    out["q05"], out["q50"], out["q95"] = (float(v) for v in np.percentile(logs, [5, 50, 95]))
    out["r2_genes"] = len(logs)
    ts = [r["transcription_rate"] for r in table.values()]
    out["R3_T"] = tuple(float(v) for v in np.percentile(ts, [5, 95]))
    out["R3_t_half"] = tuple(float(v) for v in np.percentile(k562["mean_half_life"].astype(float), [5, 95]))
    out["k562_half_life"] = human
    return out


def gate(ranges: dict) -> list[str]:
    """Each registered number against its re-reading; empty means the ranges are the registered ones."""
    reg = bridge.S5_RANGES
    pairs = [
        ("w_T", reg["R1_measurement"]["w_T"]),
        ("w_t_half", reg["R1_measurement"]["w_t_half"]),
        ("q05", reg["R2_transfer_to_a_human_cell"]["q05"]),
        ("q50", reg["R2_transfer_to_a_human_cell"]["q50"]),
        ("q95", reg["R2_transfer_to_a_human_cell"]["q95"]),
    ]
    bad = [
        f"{k}: read {ranges[k]:.6f}, registered {v}" for k, v in pairs if round(ranges[k], GATE_DECIMALS) != v
    ]
    for k, name in (("R3_T", "T"), ("R3_t_half", "t_half")):
        got = tuple(round(v, GATE_DECIMALS) for v in ranges[k])
        if got != tuple(reg["R3_no_transfer"][name]):
            bad.append(f"{k}: read {got}, registered {reg['R3_no_transfer'][name]}")
    return bad


def points(rate: float, t: float, rg: dict) -> dict[str, list[tuple[str, float, float]]]:
    """Every point of each registered range for one pair: (name, T, t_half)."""
    w_rate, wt = rg["w_T"], rg["w_t_half"]
    r1 = [
        (f"R1 T{sr:+d} t{st:+d}", rate * 2.0 ** (sr * w_rate), t * 2.0 ** (st * wt))
        for sr in (-1, 1)
        for st in (-1, 1)
    ]
    r2 = []
    for q in ("q05", "q95"):
        r = 2.0 ** rg[q]
        r2.append((f"R2 {q} rate carries over", rate, t * r))
        r2.append((f"R2 {q} copy number carries over", rate / r, t * r))
    (tlo, thi), (hlo, hhi) = rg["R3_T"], rg["R3_t_half"]
    r3 = [
        (f"R3 T{a} t{b}", x, y) for a, x in (("lo", tlo), ("hi", thi)) for b, y in (("lo", hlo), ("hi", hhi))
    ]
    return {"R0": [("R0 nominal", rate, t)], "R1": r1, "R2": r2, "R3": r3}


# ---- 2. the simulable pairs, with their labels -------------------------------------------------------


def simulable_pairs(chroms: list[str]) -> list[dict]:
    """Every measured-tier simulable (gene, cell) pair, as the audit counts them, with its label."""
    rates, _, _ = bridge_audit.load_rates()
    per, _ = bridge_audit.load_provenance()
    table = json.loads(RATE_TABLE.read_text())["rates"]
    out = []
    for chrom in chroms:
        path = Path(COMPILED_DIR) / f"noncoding_{chrom}.bio"
        m = parse(path.read_text(), path.stem)
        contexts = sorted({r.when.get("cell_type") for r in m.rules if r.when.get("cell_type")})
        for cell in contexts:
            ctx = {"cell_type": cell}
            n_rules = Counter(
                bridge.mechanism_key(r)
                for r in m.active_rules(ctx)
                if r.action in (bridge.Action.ACTIVATE, bridge.Action.INHIBIT)
            )
            b = bridge.parameterize(m, ctx, build=False, rates=rates, rate_provenance=per)
            for key, ob in b.observations.items():
                gene = key[1]
                out.append(
                    {
                        "chrom": chrom,
                        "cell": cell,
                        "element": key[0],
                        "gene": gene,
                        "kind": ob.kind,
                        "rho": ob.rho,
                        "T": b.rates_used[gene][0],
                        "t_half": table[gene]["mrna_half_life_h"],
                        "rules": n_rules[key],
                        "label": b.transferred.get(gene),
                    }
                )
    return out


def sort_key(p: dict) -> tuple:
    return (CHROMS.index(p["chrom"]), p["cell"], p["gene"], p["element"])


# ---- 3. the model, two ways ------------------------------------------------------------------------


def closed_form(rho: float, rate: float, t: float) -> dict:
    """The one-gene model's own answer: production / delta_m, and t_half for the half-way time."""
    intact = rate * t / LN2
    removed = rho * intact
    return {"intact": intact, "removed": removed, "response_time": t}


def one_gene_program(rho: float, t: float, cell: str) -> str:
    """A pair restated as a one-gene program: its fold as the observation, its half-life declared."""
    x = math.log2(rho)
    action = "activates" if rho < 1.0 else "inhibits"
    return (
        f"param mrna_half_life = {t!r} h {{ evidence: experimental"
        ' "Schwanhausser et al. 2011, Mus musculus NIH 3T3 fibroblast (transferred)" }\n'
        "element E { class: enhancer }\n"
        "gene G { symbol: G }\n"
        f"rule E {action} G {{ strength: 0.5; when: cell_type = {cell};"
        f' evidence: predicted "restated fold" effect {x!r} log2 fold change on deletion }}\n'
    )


def half_way(times: list[float], xs: list[float], start: float, end: float) -> float:
    """Time at which a monotone relaxation from `start` to `end` has covered half the distance."""
    target = start + 0.5 * (end - start)
    for i in range(1, len(xs)):
        a, b = xs[i - 1] - target, xs[i] - target
        if a == 0.0:
            return times[i - 1]
        if a * b <= 0.0 and a != b:
            return times[i - 1] + (times[i] - times[i - 1]) * a / (a - b)
    return float("nan")


def runtime(
    rho: float, rate: float, t: float, cell: str, label_source: dict, rule: str | None = None
) -> dict:
    """The runtime's answer for one pair at one point, and the label its trajectories carried."""
    m = parse(one_gene_program(rho, t, cell), "s5")
    ctx = {"cell_type": cell}
    b = bridge.parameterize(m, ctx, rates={"G": rate}, rate_provenance={"G": label_source})
    saved = grn.ACTIVATOR_COMBINATION
    if rule is not None:
        grn.ACTIVATOR_COMBINATION = rule
    try:
        dt, hours = t / 100.0, 20.0 * t
        rt = grn.NetworkRuntime(b.module, context=ctx, strict=True)
        up = rt.run(hours=hours, dt=dt, clamp={"E": bridge.INTACT}, record_every=100)
        intact = up.final()["G.mRNA"]
        down = rt.run(
            hours=hours, dt=dt, initial={"G.mRNA": intact}, clamp={"E": bridge.REMOVED}, record_every=1
        )
    finally:
        grn.ACTIVATOR_COMBINATION = saved
    removed = down.final()["G.mRNA"]
    return {
        "intact": intact,
        "removed": removed,
        "response_time": half_way(down.times, down.levels["G.mRNA"], intact, removed),
        "labels": [up.provenance(), down.provenance()],
    }


# ---- 4. the classification -------------------------------------------------------------------------


def classify(rho: float, vals: list[dict]) -> dict[str, bool]:
    """Stable (True) or not for each per-pair prediction over one set of points (S5_PREDICTIONS)."""
    intact = [v["intact"] for v in vals]
    effect = [abs(v["intact"] - v["removed"]) for v in vals]
    resp = [v["response_time"] for v in vals]
    on = {v["intact"] >= ONE_MOLECULE for v in vals}
    signs = {(v["removed"] > v["intact"]) - (v["removed"] < v["intact"]) for v in vals}
    return {
        "direction": len(signs) == 1 and signs != {0},
        "fold": all(abs(v["removed"] / v["intact"] - rho) <= 0.01 * rho for v in vals),
        "level": max(intact) <= TWO_FOLD * min(intact),
        "effect_size": min(effect) > 0 and max(effect) <= TWO_FOLD * min(effect),
        "on_off": len(on) == 1,
        "switch": switch_stable(rho, min(intact), max(intact)),
        "response_time": max(resp) <= TWO_FOLD * min(resp),
    }


def switch_stable(rho: float, lo: float, hi: float) -> bool:
    """The switch over the whole level interval [lo, hi] the points span (S5_AMENDMENT).

    Removal switches an expressed gene off when 1 <= intact < 1/RHO, a band in the level, so the answer
    is not monotone in it and the corners alone can miss the band's interior. The primary range and R3
    are connected intervals in the level, so the answer is stable iff [lo, hi] lies wholly inside the
    band or wholly outside it. An inhibition (RHO >= 1) raises the gene and can never switch it off.
    """
    if rho >= 1.0:
        return True
    a, b = ONE_MOLECULE, ONE_MOLECULE / rho
    return hi < a or lo >= b or (lo >= a and hi < b)


def answers(vals: list[dict]) -> dict:
    """The qualitative answers at the nominal point, for reporting which way a stable answer points."""
    v = vals[0]
    return {
        "on": v["intact"] >= ONE_MOLECULE,
        "switch": v["intact"] >= ONE_MOLECULE and v["removed"] < ONE_MOLECULE,
    }


def near_threshold(vals: list[dict]) -> bool:
    """A value within 1% of the one-molecule threshold, where formula and runtime may round apart."""
    return any(abs(v[k] - ONE_MOLECULE) <= 0.01 for v in vals for k in ("intact", "removed"))


def stable_comparisons(intervals: dict[str, list[tuple[float, float]]]) -> tuple[int, int]:
    """(stable, total) within-context comparisons: two ranges that do not overlap keep their order."""
    stable = total = 0
    for iv in intervals.values():
        n = len(iv)
        total += n * (n - 1) // 2
        los = sorted(lo for lo, _ in iv)
        for _, hi in iv:
            stable += n - bisect.bisect_right(los, hi)  # pairs whose whole range lies above this one
    return stable, total


def wilson(k: int, n: int, z: float = 1.959964) -> list[float]:
    if n == 0:
        return [float("nan"), float("nan")]
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return [round(c - h, 4), round(c + h, 4)]


KINDS = ("direction", "fold", "level", "effect_size", "on_off", "switch", "response_time")
PRIMARY = ("R0", "R1", "R2")


def pair_values(p: dict, rg: dict, how: str, label_source: dict | None = None) -> dict[str, list[dict]]:
    pts = points(p["T"], p["t_half"], rg)
    out = {}
    for name, plist in pts.items():
        if how == "formula":
            out[name] = [closed_form(p["rho"], x, t) for _, x, t in plist]
        else:
            out[name] = [runtime(p["rho"], x, t, p["cell"], label_source or {}) for _, x, t in plist]
    return out


def summarise(classes: list[dict[str, bool]], weights: list[int], rhos: list[float]) -> dict:
    n = len(classes)
    out = {}
    for k in KINDS:
        s = sum(1 for c in classes if c[k])
        out[k] = {
            "stable_pairs": s,
            "of_pairs": n,
            "stable_rules": sum(w for c, w in zip(classes, weights, strict=True) if c[k]),
            "of_rules": sum(weights),
            "wilson95": wilson(s, n),
        }
    both = sum(1 for c in classes if c["on_off"] and c["switch"])
    out["on_off_and_switch"] = {
        "stable_pairs": both,
        "of_pairs": n,
        "stable_rules": sum(w for c, w in zip(classes, weights, strict=True) if c["on_off"] and c["switch"]),
        "wilson95": wilson(both, n),
    }
    # S5_AMENDMENT: removal can switch off only an activated gene, so the switch and on/off are also
    # counted among the activating pairs, where the registered comparison between them is meaningful
    act = [c for c, r in zip(classes, rhos, strict=True) if r < 1.0]
    out["activating_pairs_only"] = {
        k: {
            "stable_pairs": sum(1 for c in act if c[k]),
            "of_pairs": len(act),
            "wilson95": wilson(sum(1 for c in act if c[k]), len(act)),
        }
        for k in ("on_off", "switch")
    }
    return out


# ---- 5. the run ------------------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--chroms", default="")
    ap.add_argument("--sample", type=int, default=bridge.S5_SAMPLE["size"])
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args(argv)
    chroms = args.chroms.split(",") if args.chroms else CHROMS
    t0 = time.time()

    rg = read_ranges()
    bad = gate(rg)
    if bad:
        print("GATE: the ranges read are not the registered ones:", *bad, sep="\n  ")
        return 2
    print("gate: ranges equal the registration", flush=True)

    pairs = sorted(simulable_pairs(chroms), key=sort_key)
    print(f"pairs {len(pairs)}, rules {sum(p['rules'] for p in pairs)}, {time.time() - t0:.0f}s", flush=True)
    unlabelled = [p for p in pairs if not p["label"]]
    per, _ = bridge_audit.load_provenance()

    # the registered sample, through the runtime
    n = min(args.sample, len(pairs))
    sample = random.Random(20260928).sample(pairs, n)
    s_classes: dict[str, list[dict]] = {"primary": [], "R3": []}
    agreement = Counter()
    label_failures = []
    combination_diffs = 0
    s_answers = Counter()
    for i, p in enumerate(sample):
        rv = pair_values(p, rg, "runtime", per.get(p["gene"]))
        fv = pair_values(p, rg, "formula")
        for name in rv:
            for a, b in zip(rv[name], fv[name], strict=True):
                agreement["points"] += 1
                ok_level = all(abs(a[k] - b[k]) <= 0.01 * abs(b[k]) for k in ("intact", "removed"))
                ok_time = abs(a["response_time"] - b["response_time"]) <= 0.02 * b["response_time"]
                agreement["levels_within_1pct"] += ok_level
                agreement["times_within_2pct"] += ok_time
                for lab in a["labels"]:
                    tr = lab["transferred"]
                    good = (
                        len(tr) == 1
                        and tr[0]["species_cell"] == per[p["gene"]]["species_cell"]
                        and tr[0]["source"] == per[p["gene"]]["source"]
                        and tr[0]["used_in"] == p["cell"]
                        and tr[0]["tier"] == "measured"
                        and lab["defaulted"] == []
                    )
                    if not good:
                        label_failures.append((p["gene"], p["cell"], lab))
        prim_r = [v for k in PRIMARY for v in rv[k]]
        prim_f = [v for k in PRIMARY for v in fv[k]]
        cr, cf = classify(p["rho"], prim_r), classify(p["rho"], prim_f)
        cr3, cf3 = classify(p["rho"], rv["R3"]), classify(p["rho"], fv["R3"])
        if cr == cf and cr3 == cf3:
            agreement["classification_identical"] += 1
        elif near_threshold(prim_f + fv["R3"]):
            agreement["classification_differs_near_threshold"] += 1
        else:
            agreement["classification_differs"] += 1
        s_classes["primary"].append(cr)
        s_classes["R3"].append(cr3)
        if cr["on_off"]:
            s_answers["on_off stable: on" if answers(prim_r)["on"] else "on_off stable: off"] += 1
        # the combination rule, at the nominal point, under every rule
        base = rv["R0"][0]
        for rule in bridge.S5_COMBINATION_RANGE["rules"]:
            alt = runtime(p["rho"], p["T"], p["t_half"], p["cell"], per.get(p["gene"], {}), rule=rule)
            if any(abs(alt[k] - base[k]) > 1e-9 * max(1.0, abs(base[k])) for k in ("intact", "removed")):
                combination_diffs += 1
        if (i + 1) % 100 == 0:
            print(f"sample {i + 1}/{n}, {time.time() - t0:.0f}s", flush=True)
    s_weights = [p["rules"] for p in sample]
    s_rhos = [p["rho"] for p in sample]
    census_ok = (
        agreement["classification_differs"] == 0 and agreement["levels_within_1pct"] == agreement["points"]
    )
    census_ok = census_ok and agreement["times_within_2pct"] == agreement["points"]

    # the census, by the closed form
    c_classes: dict[str, list[dict]] = {"primary": [], "R3": []}
    c_answers = Counter()
    intervals: dict[str, dict[str, dict[str, list]]] = {
        rng: {k: defaultdict(list) for k in ("level", "effect_size", "response_time")}
        for rng in ("primary", "R3")
    }
    for p in pairs:
        fv = pair_values(p, rg, "formula")
        prim = [v for k in PRIMARY for v in fv[k]]
        c = classify(p["rho"], prim)
        c_classes["primary"].append(c)
        c_classes["R3"].append(classify(p["rho"], fv["R3"]))
        a = answers(prim)
        if c["on_off"]:
            c_answers["on_off stable: on" if a["on"] else "on_off stable: off"] += 1
        if c["switch"]:
            c_answers["switch stable: a switch" if a["switch"] else "switch stable: not a switch"] += 1
        for rng, vals in (("primary", prim), ("R3", fv["R3"])):
            for k, f in (
                ("level", lambda v: v["intact"]),
                ("effect_size", lambda v: abs(v["intact"] - v["removed"])),
                ("response_time", lambda v: v["response_time"]),
            ):
                xs = [math.log2(f(v)) for v in vals]
                intervals[rng][k][p["cell"]].append((min(xs), max(xs)))
    ranks = {}
    for rng in ("primary", "R3"):
        ranks[rng] = {}
        for k, by_cell in intervals[rng].items():
            s, tot = stable_comparisons(by_cell)
            ranks[rng][k] = {"stable": s, "comparisons": tot, "share": round(s / tot, 4) if tot else None}
    weights = [p["rules"] for p in pairs]
    rhos = [p["rho"] for p in pairs]

    # the one human context: its pairs and what is matched in it
    by_cell = Counter(p["cell"] for p in pairs)
    k562 = [p for p in pairs if p["cell"] == "K562"]
    human_hl = rg["k562_half_life"]
    k562_with_hl = [p for p in k562 if p["gene"] in human_hl]
    k562_in_r2 = sum(
        1
        for p in k562_with_hl
        if p["t_half"] * 2.0 ** rg["q05"] <= human_hl[p["gene"]] <= p["t_half"] * 2.0 ** rg["q95"]
    )

    out = {
        "result": "s5_rate_ranges",
        "review_item": "item 12 S5 (second external review, 2026-09-28)",
        "registered_in": "0ea0a90",
        "simulation_status": bridge.SIMULATION_STATUS,
        "validated_human_kinetics": bridge.VALIDATED_HUMAN_KINETICS,
        "reading": (
            "stable means the same qualitative answer at every registered point; it is NOT validated."
            " Direction and fold are the observation reproduced, so their stability is not a result"
        ),
        "ranges": {
            "registered": bridge.S5_RANGES,
            "re_read": {k: v for k, v in rg.items() if k != "k562_half_life"},
            "gate": "equal to the registration in the fourth decimal",
        },
        "combination_range": bridge.S5_COMBINATION_RANGE,
        "not_propagated": list(bridge.S5_NOT_PROPAGATED),
        "predictions": bridge.S5_PREDICTIONS,
        "expected": list(bridge.S5_EXPECTED),
        "falsifiers": bridge.S5_FALSIFIERS,
        "runtime_gap": bridge.S5_RUNTIME_GAP,
        "pairs": {
            "simulable_pairs": len(pairs),
            "compiled_rules": sum(weights),
            "unlabelled_pairs": len(unlabelled),
            "labels_by_species_cell": dict(Counter(p["label"]["species_cell"] for p in pairs if p["label"])),
            "by_kind": dict(Counter(p["kind"] for p in pairs)),
            "activating": sum(1 for p in pairs if p["rho"] < 1.0),
            "inhibiting": sum(1 for p in pairs if p["rho"] > 1.0),
            "contexts": len(by_cell),
            "largest_contexts": by_cell.most_common(10),
        },
        "sample": {
            "design": bridge.S5_SAMPLE,
            "size": n,
            "agreement_with_formula": dict(agreement),
            "label_failures": len(label_failures),
            "label_failure_examples": [str(x)[:300] for x in label_failures[:3]],
            "combination_rule_differences": combination_diffs,
            "stable_primary": summarise(s_classes["primary"], s_weights, s_rhos),
            "stable_R3": summarise(s_classes["R3"], s_weights, s_rhos),
            "on_off_answers": dict(s_answers),
        },
        "census": {
            "reported": census_ok,
            "why": "the closed form agreed with the runtime on every sampled pair"
            if census_ok
            else "WITHDRAWN: the closed form disagreed with the runtime on the sample (S5_FALSIFIERS)",
            "stable_primary": summarise(c_classes["primary"], weights, rhos) if census_ok else None,
            "stable_R3": summarise(c_classes["R3"], weights, rhos) if census_ok else None,
            "amendment": bridge.S5_AMENDMENT,
            "answers": dict(c_answers) if census_ok else None,
            "rank_comparisons": ranks if census_ok else None,
        },
        "human_context": {
            "kinds": list(bridge.S5_HUMAN_CONTEXT_KINDS),
            "k562_simulable_pairs": len(k562),
            "k562_pairs_by_kind": dict(Counter(p["kind"] for p in k562)),
            "k562_pairs_with_a_human_k562_half_life": len(k562_with_hl),
            "k562_human_half_life_inside_R2": k562_in_r2,
            "in_sample_note": (
                "R2 was computed from these very genes' half-lives, so the share inside it is in-sample"
                " by construction and is not a test"
            ),
        },
        "seconds": round(time.time() - t0, 1),
    }
    if not args.no_save:
        print("saved", save_result("s5_rate_ranges", out, manifest=manifest(chroms)))
    brief = {k: out[k] for k in ("pairs", "human_context")}
    brief["sample"] = {k: v for k, v in out["sample"].items() if k != "design"}
    brief["census"] = out["census"]
    print(json.dumps(brief, indent=1, default=str))
    return 0


def manifest(chroms: list[str]) -> dict:
    paths = [Path(COMPILED_DIR) / f"noncoding_{c}.bio" for c in chroms]
    return {
        "sources": [
            {
                "accession": "the compiled non-coding programs, scripts/compile_genome_programs.py",
                "version": "as compiled from this revision's attribution/compile.py, pinned by sha256",
            },
            {
                "accession": f"schwanhausser2011: {bridge.MEASURED_RATE_SOURCES['schwanhausser2011']}",
                "version": "as cached 2026-09-28",
            },
            {
                "accession": f"schofield2018: {bridge.MEASURED_RATE_SOURCES['schofield2018']}",
                "version": "as cached 2026-09-28",
            },
        ],
        "inputs": [mf.input_entry(p, partition=None) for p in paths if p.exists()]
        + [mf.input_entry(p, partition=None) for p in (RATE_TABLE, SCHWAN, SCHOFIELD) if p.exists()],
        "assembly": "GRCh38",
        "coordinates": "n/a: gene-cell pairs by symbol and context, no interval is reported",
        "parameters": {
            "ranges": bridge.S5_RANGES,
            "sample": bridge.S5_SAMPLE,
            "predictions": bridge.S5_PREDICTIONS,
            "on_off_threshold_molecules_per_cell": ONE_MOLECULE,
            "stability_fold": TWO_FOLD,
            "combination_rules": list(bridge.S5_COMBINATION_RANGE["rules"]),
            "not_propagated": list(bridge.S5_NOT_PROPAGATED),
            "species": bridge.RATE_SPECIES,
            "status": bridge.SIMULATION_STATUS,
        },
        "exclusions": [
            "pairs outside the measured tier (borrowed median, unmeasured, not identifiable) are not run",
        ],
        "partitions": "n/a: a propagation of declared parameters, not an evaluation against held-out data",
    }


if __name__ == "__main__":
    raise SystemExit(main())
