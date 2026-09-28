# SPDX-License-Identifier: AGPL-3.0-or-later
"""Milestone 1.3's second clause re-tested against a control matched on length and scored-element count.

    uv run python scripts/clause2_element_count_control.py [--chroms chr21] [--no-save]

The clause ("the constrained-unknown blocks attributed to a gene and a tissue") is held as not met.
d717b28 put the real unknown -37.23 points below length-matched random windows on "names a coding
gene"; e1dbcf3/0c8b82d matched those windows on local coding-gene density as well and the gap stayed
-27.25 points (95% over blocks -30.91 to -23.58, n = 531), with the neutral tier equally low. That run
matched its covariate (4.11 coding TSSs against 4.05) but reported one residual imbalance and named it
the last live explanation for the gap besides "the instrument does not work": a carrying window holds
13.27 scored elements against 6.18 per carrying block, so a window has more chances to name a gene.

This script asks the same question once more against a control matched on that quantity: windows of
each block's exact length, drawn exactly as d717b28 drew them, accepted only when the number of scored
elements whose midpoint falls inside them lies in the block's own decile of scored-element count. The
per-element rate -- which does not depend on how many elements anything holds -- is reported beside it
on d717b28's own unmatched windows, because it answers the same imbalance without matching at all.

The registration (`PRE_REGISTRATION`, and its dated section in docs/ATTRIBUTION.md) was committed
before any element-count-matched figure was computed, and it includes the admissibility judgement:
the imbalance is expected by construction, so a difference that stays far below 0 is informative and
one that moves to 0 is not, on its own, a rescue. Model requests: none. It reads the stored
all-element archive one chromosome at a time and never opens the per-element response cache.
"""

from __future__ import annotations

import argparse
import bisect
import importlib.util
import random
import time
from collections.abc import Callable, Hashable
from pathlib import Path
from typing import Any

from genomeos import manifest as mf
from genomeos.attribution import organise
from genomeos.results import save_result


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).resolve().parent / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


mc = _load("clause2_matched_control")  # the density-matched control; its draw, gate and estimator
cut = mc.cut  # constrained_unknown_targets, which drew d717b28's windows

RESULT = "clause2_element_count_control"
DECILES = mc.DECILES
DRAWS = mc.DRAWS  # 50 accepted windows per block, as d717b28
SEED = mc.SEED  # 20260913 per chromosome, as d717b28, so the unmatched pass reproduces it to the digit
MATCHED_TRIES = mc.MATCHED_TRIES  # 20,000, as e1dbcf3
UNMATCHED_TRIES = mc.UNMATCHED_TRIES  # 4,000, d717b28's cap, kept for the reproduction gate
BOOTSTRAP = mc.BOOTSTRAP
CHANCE_BAND = mc.CHANCE_BAND  # 5 points, the band d717b28 registered
PRIMARY_QUESTION = mc.PRIMARY_QUESTION
QUESTIONS = mc.QUESTIONS
REPRODUCE = mc.REPRODUCE  # d717b28's pooled counts, which the unmatched pass must give back
UNMATCHED_GAP = mc.UNMATCHED_GAP  # -37.23
FALSIFIER_GAP = mc.FALSIFIER_GAP  # -18.61, half of it: e1dbcf3's threshold, kept unchanged
DENSITY_MATCHED_GAP = -27.25  # 0c8b82d, the gap this run tries to explain away
HALF_WINDOW = mc.HALF_WINDOW

PRE_REGISTRATION = {
    "registered": "2026-09-28",
    "chosen_after_the_earlier_runs": (
        "this control was chosen after seeing d717b28 and 0c8b82d, which matched on length and on local "
        "coding-gene density and left the real unknown -27.25 points below its windows while reporting "
        "one residual imbalance: 13.27 scored elements per carrying window against 6.18 per carrying "
        "block. That is admissible because those runs only named the imbalance; the matching variable, "
        "its bins, the estimator, the interval, the pass rule, the falsifier and every reading below are "
        "fixed and committed here before the element-count-matched comparison is computed on any "
        "chromosome, and the same windows are asked of the neutral tier so the control is not tuned"
    ),
    "admissibility": (
        "the imbalance is expected by construction, and the registration says so before the numbers. An "
        "organiser block is an intergenic gap: genome/unknown.unknown_blocks carves blocks out of the "
        "sequence left between annotated gene bodies, so no block contains a gene, while a window placed "
        "uniformly usually straddles one; and the budget's classifier sends CpG-island blocks to the "
        "regulatory tier, so the most promoter-like -- hence most cCRE-dense -- intergenic blocks are "
        "removed from constrained_unknown before this comparison starts. Scored elements are ENCODE "
        "SCREEN cCREs, which cluster at promoters and inside gene bodies. A constrained-unknown block is "
        "therefore element-poorer than a random window of its length by definition, not by accident, and "
        "matching on element count conditions on a variable the target set's own definition sets. Two "
        "consequences are registered now, not chosen later. (1) The matched window set is pulled toward "
        "the target set itself: windows as element-poor as a block, at the block's length, are "
        "disproportionately the same intergenic deserts the blocks are, so this control can only "
        "attenuate a difference, never inflate one. A difference that stays far below 0 is therefore "
        "informative and cannot be an element-count artefact. (2) A difference that moves to 0 is weak "
        "evidence: it is equally consistent with the element-count artefact and with the control having "
        "become a copy of the target, and on its own it does not restore clause 2 -- it only replaces "
        "'below chance' with 'unresolved by this instrument' and hands the question to the per-element "
        "estimator and to measurement. Because of (2) the per-element rate is registered as a secondary "
        "and read first when the matched difference moves: it removes the number-of-trials imbalance "
        "arithmetically, on d717b28's own unmatched windows, without matching anything away"
    ),
    "targets": (
        "the real unknown: constrained_unknown blocks with the copies out (882 blocks, 531 carrying a "
        "scored element), from organise.blocks on all 24 chromosomes; the neutral tier (2,632 blocks, "
        "1,181 carrying) as the secondary target set"
    ),
    "control": (
        "for every target block, in start order, windows of the block's exact length placed uniformly "
        "inside the span of the chromosome's scored-element midpoints and rejected on overlap with any "
        "organiser block -- d717b28's draw, same seed per chromosome, one random number per try -- and "
        "additionally rejected unless the number of scored elements whose midpoint falls inside the "
        "window lies in the same decile of scored-element count as the block's own number; 50 accepted "
        "windows per block, 20,000 tries at most, exactly e1dbcf3's cap. A block for which no window is "
        "accepted inside the cap is counted as undrawable and leaves the comparison, and the count of "
        "undrawable blocks is reported; a block that carries an element but gets no accepted window "
        "carrying one also leaves the comparison and is reported separately. Both are e1dbcf3's rules, "
        "unchanged"
    ),
    "banding": (
        "deciles of scored-element count: the nearest-rank 10th..90th percentiles of the target set's own "
        "per-block element counts, pooled genome-wide over every block of the set including the blocks "
        "that hold none, ties collapsed so there may be fewer than ten bins -- the same `decile_edges` "
        "e1dbcf3 used for coding-TSS count. Each target set is banded on its own distribution. A block "
        "holding no element falls in the lowest bin and can only draw windows holding none, so it leaves "
        "the comparison as it already did"
    ),
    "question": (
        "does a block (or window) carrying at least one scored element carry one whose deletion names a "
        "protein-coding gene (predicted_coding, the chr21 instrument's question and the one behind the "
        "milestone's 49.1% against 65.3%)"
    ),
    "primary": (
        "real unknown against its length- and element-count-matched windows, question names_a_coding_gene: "
        "the per-block matched difference, the mean over blocks that carry an element and have at least "
        "one accepted window carrying an element of (block yes, 0 or 1) minus (share of that block's "
        "carrying windows saying yes), in points -- e1dbcf3's statistic unchanged"
    ),
    "interval": (
        "95% percentile bootstrap over blocks, 10,000 resamples, numpy default_rng(20260913); a bootstrap "
        "over the 24 chromosomes (resampled with replacement, blocks pooled) is reported beside it as the "
        "sensitivity for correlation along a chromosome, and the pooled rates as d717b28 printed them"
    ),
    "secondary": [
        "the per-element rate, read first if the matched difference moves toward 0: for each compared "
        "block, the share of its own scored elements that name a coding gene minus the share of the "
        "elements of its accepted carrying windows that do, in points, with the same two bootstraps. "
        "Reported on d717b28's unmatched windows (where it answers the imbalance with no matching at "
        "all) and on the element-count-matched windows",
        "the neutral tier under the same matching (its own deciles), same statistic and interval",
        "real unknown minus neutral on the matched difference, bootstrap over blocks within each set",
        "the question moves_a_gene (any gene, the script's own) on both sets",
    ],
    "sensitivities": [
        "matching jointly on the scored-element decile and the coding-TSS-count decile of e1dbcf3, so "
        "both imbalances are closed at once; blocks with no window in both bins inside the cap are "
        "undrawable and reported",
        "matching on the exact scored-element count rather than its decile",
        "rejecting windows only on overlap with the target set itself rather than every organiser block",
    ],
    "gate": (
        "before any element-count-matched figure is read, the same draw function with the element test "
        "switched off and d717b28's 4,000-try cap must reproduce d717b28's pooled counts for both sets to "
        "the digit (real unknown 44,100 windows drawn, 31,676 carrying, 21,217 naming a coding gene, 531 "
        "and 158 blocks; neutral 131,600 / 88,317 / 57,212, 1,181 and 340). If it does not, no matched "
        "figure is reported"
    ),
    "balance_check": (
        "reported, not gated: scored elements per carrying block against per carrying window, which is "
        "the quantity being matched and must come close to equal for the matching to have worked, and "
        "the mean coding-TSS count of the compared blocks against their windows', which is now free and "
        "may drift back apart"
    ),
    "pass_rule": (
        "clause 2 as the milestone words it -- the blocks attributed to a gene and a tissue, a labelled "
        "lead where no measurement exists (section 5 item 10) -- passes on this instrument only if the "
        "primary's 95% interval lies wholly above 0: a real-unknown block then names a coding gene more "
        "often than sequence of the same length holding as many scored elements. The clause's third part "
        "(scored against measurement) is untouched by any outcome here"
    ),
    "readings": {
        "above": (
            "interval wholly above 0: clause 2 passes as a labelled lead against an element-count-matched "
            "control, with measurement still outstanding"
        ),
        "at_chance": (
            "interval covers 0: the element-count imbalance was the gap. The below-chance readings of "
            "d717b28 and 0c8b82d were an element-count artefact and clause 2's status must be "
            "re-examined: it stops being 'below chance' and becomes 'unresolved by this instrument', not "
            "'met'. Per the admissibility judgement this outcome is weak on its own, because the matched "
            "control is pulled toward the target set; the per-element secondary is read next, and the "
            "clause stays not met until a measured arm answers it"
        ),
        "below": (
            "interval wholly below 0: element count is not the cause of the gap either. With length, "
            "coding-gene density (0c8b82d) and now the number of scored elements all matched and the "
            "difference still far below 0, no covariate named so far explains the reading, and clause 2 "
            "stays not met and below chance on the ground it already stood on"
        ),
    },
    "falsifier": (
        "the element-count account of the gap is falsified if the element-count-matched difference keeps "
        f"more than half of d717b28's unmatched gap, i.e. is below {FALSIFIER_GAP} points (half of "
        f"{UNMATCHED_GAP}) -- e1dbcf3's threshold, unchanged so the two controls are read on one scale. "
        f"The share of the unmatched gap that element count closes, and the share of the density-matched "
        f"{DENSITY_MATCHED_GAP} points it closes, are both reported"
    ),
    "second_route_not_run": (
        "the measured arm stays registered and unrun: it belongs to its own lane and no figure from it is "
        "read here. Nothing measured exists for these blocks (0.45% of unknown space has any assay)"
    ),
    "cost": "0 AlphaGenome requests; the response cache is never opened; one chromosome's archive in memory",
}


# ---- the draw -----------------------------------------------------------------------------


def element_matched_windows(
    targets: list[dict[str, Any]],
    exclude: list[dict[str, Any]],
    rows: list[dict[str, Any]],
    predicates: dict[str, Callable[[dict[str, Any]], bool]],
    key_fn: Callable[[int, int], Hashable] | None = None,
    seed: int = SEED,
    draws: int = DRAWS,
    max_tries: int = MATCHED_TRIES,
    raw_fn: Callable[[int], float] | None = None,
) -> list[dict[str, Any]]:
    """Per target block, its answer and the answers of its accepted windows.

    `constrained_unknown_targets.matched_random_windows` step for step (one randrange per try, then the
    overlap test), with one added rejection: when `key_fn` is given, a window whose key differs from the
    block's is rejected after the overlap test and costs no extra random number. The key is computed
    from the midpoint and the number of scored elements inside, so element-count bands, the joint band
    with coding-TSS count, and the exact count all go through the same function. With `key_fn=None` and
    `max_tries=4000` it draws d717b28's windows exactly, which is the reproduction gate.

    Each record also counts elements, not only blocks: `elements_yes` and `window_elements_yes` carry
    the per-element numerators for the registered per-element secondary.
    """
    if not targets or not rows:
        return []
    mids = sorted(((r["start"] + r["end"]) // 2, i) for i, r in enumerate(rows))
    keys = [m for m, _ in mids]
    lo, hi = keys[0], keys[-1]
    unknown = sorted((b["start"], b["end"]) for b in exclude)
    u_starts = [s for s, _ in unknown]

    def overlaps(s: int, t: int) -> bool:
        i = max(0, bisect.bisect_right(u_starts, s) - 1)
        while i < len(unknown) and unknown[i][0] < t:
            if unknown[i][1] > s:
                return True
            i += 1
        return False

    def inside_window(s: int, t: int) -> list[dict[str, Any]]:
        return [rows[i] for _, i in mids[bisect.bisect_left(keys, s) : bisect.bisect_left(keys, t)]]

    rng = random.Random(seed)
    out = []
    for b in sorted(targets, key=lambda b: b["start"]):
        length = b["length"]
        bmid = (b["start"] + b["end"]) // 2
        sel_b = inside_window(b["start"], b["end"])
        want = key_fn(bmid, len(sel_b)) if key_fn else None
        rec: dict[str, Any] = {
            "start": b["start"],
            "length": length,
            "key": want,
            "raw": raw_fn(bmid) if raw_fn else None,
            "elements": len(sel_b),
            "elements_yes": {q: sum(1 for r in sel_b if t(r)) for q, t in predicates.items()},
            "yes": {q: any(t(r) for r in sel_b) for q, t in predicates.items()},
            "drawn": 0,
            "carrying": 0,
            "window_elements": 0,
            "window_raw": 0.0,
            "windows_yes": dict.fromkeys(predicates, 0),
            "window_elements_yes": dict.fromkeys(predicates, 0),
        }
        got = tries = 0
        while got < draws and tries < max_tries:
            tries += 1
            top = hi - length
            if top <= lo:
                break
            s = rng.randrange(lo, top)
            if overlaps(s, s + length):
                continue
            wmid = s + length // 2
            sel = inside_window(s, s + length)
            if key_fn is not None and key_fn(wmid, len(sel)) != want:
                continue
            got += 1
            if raw_fn:
                rec["window_raw"] += raw_fn(wmid) or 0
            if not sel:
                continue
            rec["carrying"] += 1
            rec["window_elements"] += len(sel)
            for q, t in predicates.items():
                n_yes = sum(1 for r in sel if t(r))
                rec["windows_yes"][q] += n_yes > 0
                rec["window_elements_yes"][q] += n_yes
        rec["drawn"] = got
        rec["tries"] = tries
        out.append(rec)
    return out


# ---- the per-element secondary ------------------------------------------------------------


def per_element_differences(records: list[dict[str, Any]], q: str) -> list[float]:
    """Block share of elements naming a gene minus its windows' share, per compared block."""
    return [
        r["elements_yes"][q] / r["elements"] - r["window_elements_yes"][q] / r["window_elements"]
        for r in mc.compared(records)
        if r["window_elements"]
    ]


def per_element_summary(
    by_chrom: dict[str, list[dict[str, Any]]], q: str, seed: int = SEED, n_boot: int = BOOTSTRAP
) -> dict[str, Any]:
    """The registered per-element secondary: the same two bootstraps on a rate that ignores how many
    elements anything holds."""
    import numpy as np

    chroms = [c for c in by_chrom if per_element_differences(by_chrom[c], q)]
    per = {c: np.array(per_element_differences(by_chrom[c], q)) for c in chroms}
    d = np.concatenate([per[c] for c in chroms]) if chroms else np.array([])
    comp = [r for c in by_chrom for r in mc.compared(by_chrom[c]) if r["window_elements"]]
    out: dict[str, Any] = {"question": q, "n_blocks_compared": len(comp)}
    if not len(d):
        out["per_element_difference_points"] = None
        return out
    rng = np.random.default_rng(seed)
    boots = d[rng.integers(0, len(d), size=(n_boot, len(d)))].mean(axis=1)
    sums = np.array([per[c].sum() for c in chroms])
    ns = np.array([len(per[c]) for c in chroms])
    cidx = rng.integers(0, len(chroms), size=(n_boot, len(chroms)))
    cboots = sums[cidx].sum(axis=1) / ns[cidx].sum(axis=1)
    be, bn = sum(r["elements_yes"][q] for r in comp), sum(r["elements"] for r in comp)
    we, wn = sum(r["window_elements_yes"][q] for r in comp), sum(r["window_elements"] for r in comp)
    out.update(
        per_element_difference_points=round(100 * float(d.mean()), 2),
        ci95_over_blocks=[round(100 * float(x), 2) for x in np.percentile(boots, [2.5, 97.5])],
        ci95_over_chromosomes=[round(100 * float(x), 2) for x in np.percentile(cboots, [2.5, 97.5])],
        pooled={
            "block_elements": bn,
            "block_elements_yes": be,
            "block_rate": round(be / bn, 4) if bn else None,
            "window_elements": wn,
            "window_elements_yes": we,
            "window_rate": round(we / wn, 4) if wn else None,
            "difference_in_points": round(100 * (be / bn - we / wn), 2) if bn and wn else None,
        },
    )
    return out


def reading(ci: list[float], point: float) -> dict[str, Any]:
    lo, hi = ci
    key = "above" if lo > 0 else "below" if hi < 0 else "at_chance"
    return {
        "outcome": key,
        "clause_2": "passes as a labelled lead" if key == "above" else "stays not met",
        "below_the_chance_band": point < -CHANCE_BAND,
        "falsifier_fired": point < FALSIFIER_GAP,
        "share_of_the_unmatched_gap_closed": round(1 - point / UNMATCHED_GAP, 3),
        "share_of_the_density_matched_gap_closed": round(1 - point / DENSITY_MATCHED_GAP, 3),
        "text": PRE_REGISTRATION["readings"][key],
    }


# ---- the run ------------------------------------------------------------------------------


def collect(chroms: list[str]) -> dict[str, Any]:
    t0 = time.time()
    # pass 1: blocks, element counts and TSSs only, so every band exists before a window is drawn
    blocks: dict[str, list[dict[str, Any]]] = {}
    tss: dict[str, list[int]] = {}
    counts: dict[str, list[int]] = {"real_unknown": [], "neutral": []}
    for c in chroms:
        blocks[c] = organise.blocks(c)
        tss[c] = mc.coding_tss(c)
        els = sorted(cut.elements_of(c), key=lambda e: e["start"])
        emids = sorted((e["start"] + e["end"]) // 2 for e in els)
        for name, sel in mc.target_sets(blocks[c]).items():
            for b in sel:
                counts[name].append(
                    bisect.bisect_left(emids, b["end"]) - bisect.bisect_left(emids, b["start"])
                )
        print(f"{c}: {len(blocks[c])} blocks, {len(tss[c])} coding TSSs, {len(els)} elements", flush=True)
        del els, emids
    edges = {
        f"{name}:element_count": mc.decile_edges(counts[name], DECILES)
        for name in ("real_unknown", "neutral")
    }
    edges.update(
        {
            f"{name}:tss_count": mc.decile_edges(
                [
                    mc.tss_count(tss[c], (b["start"] + b["end"]) // 2)
                    for c in chroms
                    for b in mc.target_sets(blocks[c])[name]
                ]
            )
            for name in ("real_unknown", "neutral")
        }
    )

    runs: dict[str, dict[str, list[dict[str, Any]]]] = {
        k: {}
        for k in (
            "unmatched_real_unknown",
            "unmatched_neutral",
            "matched_real_unknown",
            "matched_neutral",
            "sens_joint_real_unknown",
            "sens_exact_count_real_unknown",
            "sens_exclude_targets_only_real_unknown",
        )
    }
    qs = list(QUESTIONS)
    for c in chroms:
        els = cut.elements_of(c)
        els.sort(key=lambda e: e["start"])
        got = cut.read_chromosome(c, els)  # the organiser's blocks as d717b28 saw them
        if not got:
            continue
        sets = mc.target_sets(got)
        t = tss[c]
        for name in ("real_unknown", "neutral"):
            runs[f"unmatched_{name}"][c] = element_matched_windows(
                sets[name], got, els, QUESTIONS, None, max_tries=UNMATCHED_TRIES
            )
            e = edges[f"{name}:element_count"]
            runs[f"matched_{name}"][c] = element_matched_windows(
                sets[name],
                got,
                els,
                QUESTIONS,
                lambda mid, n, e=e: mc.bin_of(n, e),
                raw_fn=lambda mid, t=t: mc.tss_count(t, mid),
            )
        e = edges["real_unknown:element_count"]
        te = edges["real_unknown:tss_count"]
        runs["sens_joint_real_unknown"][c] = element_matched_windows(
            sets["real_unknown"],
            got,
            els,
            QUESTIONS,
            lambda mid, n, e=e, te=te, t=t: (mc.bin_of(n, e), mc.bin_of(mc.tss_count(t, mid), te)),
            raw_fn=lambda mid, t=t: mc.tss_count(t, mid),
        )
        runs["sens_exact_count_real_unknown"][c] = element_matched_windows(
            sets["real_unknown"], got, els, QUESTIONS, lambda mid, n: n
        )
        runs["sens_exclude_targets_only_real_unknown"][c] = element_matched_windows(
            sets["real_unknown"],
            sets["real_unknown"],
            els,
            QUESTIONS,
            lambda mid, n, e=e: mc.bin_of(n, e),
            raw_fn=lambda mid, t=t: mc.tss_count(t, mid),
        )
        del els, got
        print(f"{c}: drawn ({time.time() - t0:.0f} s)", flush=True)

    gate = {}
    for name in ("real_unknown", "neutral"):
        recs = [r for c in chroms for r in runs[f"unmatched_{name}"].get(c, [])]
        gate[name] = {"got": mc.lifted_counts(recs, qs), "want": REPRODUCE[name]}
        gate[name]["same"] = gate[name]["got"] == gate[name]["want"]
    passed = all(g["same"] for g in gate.values())
    out: dict[str, Any] = {
        "result": RESULT,
        "chromosomes": chroms,
        "registration": PRE_REGISTRATION,
        "bands": edges,
        "gate_reproduces_d717b28": {"passed": passed, **gate},
    }
    if not passed:
        out["reading"] = "the reproduction gate failed: no element-count-matched figure is reported"
        out["seconds"] = round(time.time() - t0, 1)
        return out
    res = {run: {q: mc.summarise(by_chrom, q) for q in qs} for run, by_chrom in runs.items()}

    import numpy as np

    contrast = {}
    for q in qs:
        a, b = res["matched_real_unknown"][q], res["matched_neutral"][q]
        if a.get("matched_difference_points") is None or b.get("matched_difference_points") is None:
            continue
        diff = a["_boots"] - b["_boots"]
        contrast[q] = {
            "points": round(a["matched_difference_points"] - b["matched_difference_points"], 2),
            "ci95": [round(100 * float(x), 2) for x in np.percentile(diff, [2.5, 97.5])],
        }
    for run in res.values():
        for s in run.values():
            s.pop("_boots", None)
    prim = res["matched_real_unknown"][PRIMARY_QUESTION]
    out.update(
        primary=prim,
        primary_reading=reading(prim["ci95_over_blocks"], prim["matched_difference_points"]),
        secondary_per_element={
            "on_unmatched_windows": {q: per_element_summary(runs["unmatched_real_unknown"], q) for q in qs},
            "on_element_matched_windows": {
                q: per_element_summary(runs["matched_real_unknown"], q) for q in qs
            },
            "neutral_on_unmatched_windows": {
                q: per_element_summary(runs["unmatched_neutral"], q) for q in qs
            },
        },
        secondary_neutral=res["matched_neutral"],
        secondary_neutral_reading={
            q: reading(s["ci95_over_blocks"], s["matched_difference_points"])["outcome"]
            for q, s in res["matched_neutral"].items()
            if s.get("matched_difference_points") is not None
        },
        secondary_real_minus_neutral=contrast,
        secondary_real_unknown_moves_a_gene=res["matched_real_unknown"]["moves_a_gene"],
        sensitivities={
            "element_and_tss_deciles_jointly": res["sens_joint_real_unknown"],
            "exact_element_count": res["sens_exact_count_real_unknown"],
            "exclude_targets_only": res["sens_exclude_targets_only_real_unknown"],
        },
        unmatched_same_estimator={
            "real_unknown": res["unmatched_real_unknown"],
            "neutral": res["unmatched_neutral"],
        },
        seconds=round(time.time() - t0, 1),
    )
    return out


@mf.depends_on_models("alphagenome")  # the sweep's model, as far as the disk says (R9)
def manifest(chroms: list[str]) -> dict[str, Any]:
    m = mc.manifest(chroms)
    m["sources"] = [
        *m["sources"],
        {
            "accession": "this repository, data/results/clause2_matched_control.json at 0c8b82d",
            "version": "git 0c8b82d, the density-matched control whose residual element-count imbalance "
            "this run matches on",
        },
    ]
    m["parameters"] = {
        "half_window_bp": HALF_WINDOW,
        "deciles": DECILES,
        "draws_per_block": DRAWS,
        "matched_tries": MATCHED_TRIES,
        "unmatched_tries": UNMATCHED_TRIES,
        "seed_per_chromosome": SEED,
        "bootstrap": BOOTSTRAP,
        "chance_band_points": CHANCE_BAND,
        "falsifier_gap_points": FALSIFIER_GAP,
        "density_matched_gap_points": DENSITY_MATCHED_GAP,
        "primary_question": PRIMARY_QUESTION,
        "matched_on": "block length (exact) and scored-element count (decile of the target set's own)",
    }
    m["exclusions"] = [
        "constrained_unknown blocks that are copies are out of the real unknown",
        "an element belongs to the block or window holding its midpoint",
        "windows overlapping any organiser block are rejected (primary); windows in another "
        "scored-element decile than their block are rejected; blocks with no accepted window are "
        "undrawable and leave the comparison",
        "blocks with no scored element, or no accepted window carrying one, leave the comparison",
    ]
    return m


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--chroms", default="", help="comma-separated; default every chromosome with a table")
    ap.add_argument("--no-save", action="store_true", help="print without writing the result")
    args = ap.parse_args(argv)
    chroms = (
        args.chroms.split(",")
        if args.chroms
        else sorted((p.stem for p in cut.ELEMENTS.glob("chr*.json")), key=lambda c: (len(c), c))
    )
    out = collect(chroms)
    if args.no_save or args.chroms:
        print("not saved (a partial run never overwrites the genome's)")
    else:
        print(f"saved {save_result(RESULT, out, manifest=manifest(chroms))}")
    print("gate passed:", out["gate_reproduces_d717b28"]["passed"])
    if "primary" in out:
        p = out["primary"]
        print(
            f"primary: {p['matched_difference_points']:+.2f} points, CI {p['ci95_over_blocks']}, "
            f"n={p['n_blocks_compared']} blocks; balance {p['balance']}"
        )
        print("reading:", out["primary_reading"]["outcome"], out["primary_reading"]["clause_2"])
        pe = out["secondary_per_element"]["on_unmatched_windows"][PRIMARY_QUESTION]
        print(
            f"per element (unmatched windows): {pe['per_element_difference_points']:+.2f} points, "
            f"CI {pe['ci95_over_blocks']}, n={pe['n_blocks_compared']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
