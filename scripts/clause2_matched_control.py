# SPDX-License-Identifier: AGPL-3.0-or-later
"""Milestone 1.3's second clause re-tested against a control matched on length and gene density.

    uv run python scripts/clause2_matched_control.py [--chroms chr21] [--no-save]

The clause ("the constrained-unknown blocks attributed to a gene and a tissue") is held as not met
because the real unknown names a coding gene far below length-matched random windows (d717b28:
29.8% of 531 blocks against 67.0% of the windows, -37.2 points). The same control puts the neutral
tier -36.0 points down, and inside length deciles the real unknown sits +4.6 points above neutral:
the audit named the confounder, which is that organiser blocks lie in gene-poorer sequence than a
uniformly placed window. This script asks the same question once more against a control it can pass:
windows of each block's exact length, drawn exactly as d717b28 drew them, but accepted only when the
number of coding TSSs in the 1 Mb window around their midpoint falls in the block's own decile.

The registration (`PRE_REGISTRATION`, and its dated section in docs/ATTRIBUTION.md) was committed
before any matched figure was computed. Model requests: none. It reads the stored all-element archive
(`data/knowledge/alphagenome/all_elements`, one chromosome at a time) and never opens the response
cache, so no cache reader runs at all.
"""

from __future__ import annotations

import argparse
import bisect
import importlib.util
import random
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

from genomeos import manifest as mf
from genomeos.attribution import organise
from genomeos.results import save_result

_spec = importlib.util.spec_from_file_location(
    "constrained_unknown_targets", Path(__file__).resolve().parent / "constrained_unknown_targets.py"
)
cut = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cut)

RESULT = "clause2_matched_control"
HALF_WINDOW = 524_288  # half of AlphaGenome's 1 Mb input (dna_client.SEQUENCE_LENGTH_1MB = 1,048,576)
DECILES = 10
DRAWS = cut.RANDOM_DRAWS  # 50 accepted windows per block, as d717b28
SEED = cut.SEED  # 20260913 per chromosome, as d717b28, so the unmatched run reproduces it to the digit
MATCHED_TRIES = 20_000  # a density-matched window is rarer than an unmatched one; d717b28 allowed 4,000
UNMATCHED_TRIES = 4_000  # d717b28's cap, kept for the reproduction gate
BOOTSTRAP = 10_000
CHANCE_BAND = cut.CHANCE_BAND  # 5 points, the band d717b28 registered
PRIMARY_QUESTION = "names_a_coding_gene"
QUESTIONS: dict[str, Callable[[dict[str, Any]], bool]] = {
    "names_a_coding_gene": cut.moves_coding,
    "moves_a_gene": cut.moves,
}
# d717b28's pooled counts, which the unmatched pass of this script must give back to the digit
REPRODUCE = {
    "real_unknown": {
        "windows_drawn": 44_100,
        "windows_carrying_an_element": 31_676,
        "blocks_carrying_an_element": 531,
        "names_a_coding_gene": {"blocks_yes": 158, "windows_yes": 21_217},
        "moves_a_gene": {"blocks_yes": 331, "windows_yes": 27_254},
    },
    "neutral": {
        "windows_drawn": 131_600,
        "windows_carrying_an_element": 88_317,
        "blocks_carrying_an_element": 1_181,
        "names_a_coding_gene": {"blocks_yes": 340, "windows_yes": 57_212},
        "moves_a_gene": {"blocks_yes": 768, "windows_yes": 74_466},
    },
}
UNMATCHED_GAP = -37.23  # d717b28, real unknown, names a coding gene, pooled points
FALSIFIER_GAP = round(UNMATCHED_GAP / 2, 2)  # -18.62: matching must close at least half the gap

PRE_REGISTRATION = {
    "registered": "2026-09-28",
    "chosen_after_the_audit": (
        "this control was chosen after seeing the audit (d717b28 and lane-scoring 2026-09-27), which "
        "found the neutral tier as far below its windows as the real unknown and the real unknown +4.6 "
        "points above neutral inside length deciles. That is admissible because the audit only named "
        "the confounder (local gene density); the matching variable, its bins, the estimator, the "
        "interval, the pass rule and the falsifier are fixed and committed here before the matched "
        "comparison is computed on any chromosome, and the same windows are asked of the neutral tier "
        "so the control is not tuned to the target"
    ),
    "targets": (
        "the real unknown: constrained_unknown blocks with the copies out (882 blocks, 531 carrying a "
        "scored element), from organise.blocks on all 24 chromosomes; the neutral tier (2,632 blocks, "
        "1,181 carrying) as the secondary target set"
    ),
    "control": (
        "for every target block, in start order, windows of the block's exact length placed uniformly "
        "inside the span of the chromosome's scored-element midpoints and rejected on overlap with any "
        "organiser block -- d717b28's draw, same seed per chromosome -- and additionally rejected unless "
        "the count of GENCODE protein-coding TSSs within 524,288 bp either side of the window's midpoint "
        "falls in the same decile as the block's own count; 50 accepted windows per block, 20,000 tries "
        "at most; a block with no accepted window is counted as undrawable and leaves the comparison. "
        "Decile cuts are the nearest-rank 10th..90th percentiles of the target set's own block counts "
        "pooled genome-wide (ties collapse, so there may be fewer than ten bins); each target set is "
        "binned on its own distribution"
    ),
    "question": (
        "does a block (or window) carrying at least one scored element carry one whose deletion names a "
        "protein-coding gene (predicted_coding, the chr21 instrument's question and the one behind the "
        "milestone's 49.1% against 65.3%)"
    ),
    "primary": (
        "real unknown against its density- and length-matched windows, question names_a_coding_gene: the "
        "per-block matched difference, the mean over blocks that carry an element and have at least one "
        "accepted window carrying an element of (block yes, 0 or 1) minus (share of that block's "
        "carrying windows saying yes), in points"
    ),
    "interval": (
        "95% percentile bootstrap over blocks, 10,000 resamples, numpy default_rng(20260913); a bootstrap "
        "over the 24 chromosomes (resampled with replacement, blocks pooled) is reported beside it as the "
        "sensitivity for correlation along a chromosome, and the pooled rates as d717b28 printed them"
    ),
    "secondary": [
        "the neutral tier under the same matching (its own deciles), same statistic and interval",
        "real unknown minus neutral on the matched difference, bootstrap over blocks within each set",
        "the question moves_a_gene (any gene, the script's own) on both sets",
    ],
    "sensitivities": [
        "matching on the distance from the midpoint to the nearest coding TSS, in deciles, instead of "
        "the count",
        "rejecting windows only on overlap with the target set itself rather than every organiser block",
    ],
    "gate": (
        "before any matched figure is read, the same function with the density test switched off and "
        "d717b28's 4,000-try cap must reproduce d717b28's pooled counts for both sets to the digit "
        "(real unknown 44,100 windows drawn, 31,676 carrying, 21,217 naming a coding gene, 531 and 158 "
        "blocks; neutral 131,600 / 88,317 / 57,212, 1,181 and 340). If it does not, no matched figure is "
        "reported"
    ),
    "balance_check": (
        "reported, not gated: mean coding-TSS count of the compared blocks against the mean of their "
        "windows' means, and scored elements per carrying block against per carrying window. A matched "
        "window that still holds more elements than its block is a residual confounder the reading names"
    ),
    "pass_rule": (
        "clause 2 as the milestone words it -- the blocks attributed to a gene and a tissue, a labelled "
        "lead where no measurement exists (section 5 item 10) -- passes on this instrument only if the "
        "primary's 95% interval lies wholly above 0: a real-unknown block then names a coding gene more "
        "often than sequence of the same length and gene density, and its lead carries information. The "
        "clause's third part (scored against measurement) is untouched by any outcome here"
    ),
    "readings": {
        "above": "interval wholly above 0: clause 2 passes as a lead; it moves from not met to met as a "
        "labelled lead above a matched control, with measurement still outstanding",
        "at_chance": "interval covers 0: the gap was the confounder and the lead carries no block-level "
        "information; clause 2 stays not met, restated from 'below chance' to 'at chance against a "
        "matched control'",
        "below": "interval wholly below 0: clause 2 stays not met, now against a control it could have "
        f"passed; if the point is also below -{CHANCE_BAND} points it stays 'below chance' as held",
    },
    "falsifier": (
        "the audit's account (most of the gap is gene density) is falsified if the matched difference "
        f"keeps more than half of the unmatched gap, i.e. is below {FALSIFIER_GAP} points "
        f"(half of d717b28's {UNMATCHED_GAP}); then density is not what drives the gap and the "
        "milestone's reading stands on its original ground"
    ),
    "second_route_not_run": (
        "a measured arm per block: the designed experiment (284,001 oligos over 878 blocks, unrun) tiles "
        "each block and a density-matched window set, so an MPRA readout would give measured activity "
        "per block against its own control rather than a model's named gene. Registered as an option "
        "only; nothing measured exists for it (0.45% of unknown space has any assay) and it is not run"
    ),
    "cost": "0 AlphaGenome requests; the response cache is never opened; one chromosome's archive in memory",
}


# ---- the covariate ------------------------------------------------------------------------


def coding_tss(chrom: str) -> list[int]:
    """GENCODE protein-coding TSSs of one chromosome, the way unknown_scoring.annotate places them."""
    from genomeos.genome import Annotation
    from genomeos.genome.annotation import default_gencode

    ann = Annotation.from_gff3(default_gencode({chrom}), {chrom})
    return sorted(
        (g.locus.end - 1 if g.locus.strand.value == "-" else g.locus.start)
        for g in ann.protein_coding()
        if g.locus.chrom == chrom
    )


def tss_count(tss: list[int], mid: int, half: int = HALF_WINDOW) -> int:
    """Coding TSSs within `half` bp either side of `mid` (half-open)."""
    return bisect.bisect_left(tss, mid + half) - bisect.bisect_left(tss, mid - half)


def tss_distance(tss: list[int], mid: int) -> int:
    """Distance from `mid` to the nearest coding TSS; a chromosome with none returns a huge value."""
    j = bisect.bisect_left(tss, mid)
    near = [abs(tss[k] - mid) for k in (j - 1, j) if 0 <= k < len(tss)]
    return min(near) if near else 10**12


COVARIATES: dict[str, Callable[[list[int], int], int]] = {
    "tss_count": tss_count,
    "tss_distance": tss_distance,
}


def decile_edges(values: list[float], n: int = DECILES) -> list[float]:
    """Nearest-rank cut points, ties collapsed -- as `by_length_decile` cuts lengths."""
    pool = sorted(values)
    if not pool:
        return []
    return sorted({pool[int(len(pool) * i / n)] for i in range(1, n)})


def bin_of(value: float, edges: list[float]) -> int:
    return sum(value >= c for c in edges)


# ---- the draw -----------------------------------------------------------------------------


def matched_windows(
    targets: list[dict[str, Any]],
    exclude: list[dict[str, Any]],
    rows: list[dict[str, Any]],
    predicates: dict[str, Callable[[dict[str, Any]], bool]],
    bin_fn: Callable[[int], int] | None = None,
    seed: int = SEED,
    draws: int = DRAWS,
    max_tries: int = MATCHED_TRIES,
    raw_fn: Callable[[int], float] | None = None,
) -> list[dict[str, Any]]:
    """Per target block, its answer and the answers of its accepted windows.

    The draw is `constrained_unknown_targets.matched_random_windows` step for step (one randrange per
    try, the overlap test, then acceptance), with one added rejection: when `bin_fn` is given, a window
    whose midpoint's bin differs from the block's is rejected after the overlap test and costs no extra
    random number. With `bin_fn=None` and `max_tries=4000` it therefore draws d717b28's windows exactly.
    `raw_fn` gives the covariate's raw value for the balance check.
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
        want = bin_fn(bmid) if bin_fn else None
        sel_b = inside_window(b["start"], b["end"])
        rec: dict[str, Any] = {
            "start": b["start"],
            "length": length,
            "bin": want,
            "raw": raw_fn(bmid) if raw_fn else None,
            "elements": len(sel_b),
            "yes": {q: any(t(r) for r in sel_b) for q, t in predicates.items()},
            "drawn": 0,
            "carrying": 0,
            "window_elements": 0,
            "window_raw": 0.0,
            "windows_yes": dict.fromkeys(predicates, 0),
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
            if bin_fn is not None and bin_fn(wmid) != want:
                continue
            got += 1
            sel = inside_window(s, s + length)
            if raw_fn:
                rec["window_raw"] += raw_fn(wmid) or 0
            if not sel:
                continue
            rec["carrying"] += 1
            rec["window_elements"] += len(sel)
            for q, t in predicates.items():
                rec["windows_yes"][q] += any(t(r) for r in sel)
        rec["drawn"] = got
        out.append(rec)
    return out


def lifted_counts(records: list[dict[str, Any]], questions: list[str]) -> dict[str, Any]:
    """The records summed into d717b28's pooled shape, for the reproduction gate."""
    return {
        "windows_drawn": sum(r["drawn"] for r in records),
        "windows_carrying_an_element": sum(r["carrying"] for r in records),
        "blocks_carrying_an_element": sum(1 for r in records if r["elements"]),
        **{
            q: {
                "blocks_yes": sum(1 for r in records if r["elements"] and r["yes"][q]),
                "windows_yes": sum(r["windows_yes"][q] for r in records),
            }
            for q in questions
        },
    }


# ---- the statistic ------------------------------------------------------------------------


def compared(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Blocks that carry an element and have at least one accepted window carrying one."""
    return [r for r in records if r["elements"] and r["carrying"]]


def differences(records: list[dict[str, Any]], q: str) -> list[float]:
    return [float(r["yes"][q]) - r["windows_yes"][q] / r["carrying"] for r in compared(records)]


def summarise(
    by_chrom: dict[str, list[dict[str, Any]]], q: str, seed: int = SEED, n_boot: int = BOOTSTRAP
) -> dict[str, Any]:
    """The per-block matched difference with both bootstraps, the pooled rates and the balance check."""
    import numpy as np

    chroms = [c for c in by_chrom if compared(by_chrom[c])]
    per = {c: np.array(differences(by_chrom[c], q)) for c in chroms}
    allrec = [r for c in by_chrom for r in by_chrom[c]]
    comp = compared(allrec)
    d = np.concatenate([per[c] for c in chroms]) if chroms else np.array([])
    out: dict[str, Any] = {
        "question": q,
        "target_blocks": len(allrec),
        "blocks_carrying_an_element": sum(1 for r in allrec if r["elements"]),
        "undrawable_blocks": sum(1 for r in allrec if not r["drawn"]),
        "carrying_blocks_without_a_carrying_window": sum(
            1 for r in allrec if r["elements"] and not r["carrying"]
        ),
        "n_blocks_compared": len(comp),
        "windows_drawn": sum(r["drawn"] for r in allrec),
        "windows_carrying_an_element_compared": sum(r["carrying"] for r in comp),
    }
    if not len(d):
        out["matched_difference_points"] = None
        return out
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(d), size=(n_boot, len(d)))
    boots = d[idx].mean(axis=1)
    sums = np.array([per[c].sum() for c in chroms])
    ns = np.array([len(per[c]) for c in chroms])
    cidx = rng.integers(0, len(chroms), size=(n_boot, len(chroms)))
    cboots = sums[cidx].sum(axis=1) / ns[cidx].sum(axis=1)
    by = sum(1 for r in comp if r["yes"][q])
    wy = sum(r["windows_yes"][q] for r in comp)
    wc = sum(r["carrying"] for r in comp)
    out.update(
        matched_difference_points=round(100 * float(d.mean()), 2),
        ci95_over_blocks=[round(100 * float(x), 2) for x in np.percentile(boots, [2.5, 97.5])],
        ci95_over_chromosomes=[round(100 * float(x), 2) for x in np.percentile(cboots, [2.5, 97.5])],
        block_rate=round(by / len(comp), 4),
        mean_window_rate_per_block=round(
            float(np.mean([r["windows_yes"][q] / r["carrying"] for r in comp])), 4
        ),
        pooled={
            "blocks_yes": by,
            "windows_yes": wy,
            "windows_carrying": wc,
            "block_rate": round(by / len(comp), 4),
            "window_rate": round(wy / wc, 4),
            "difference_in_points": round(100 * (by / len(comp) - wy / wc), 2),
        },
    )
    if comp[0]["raw"] is not None:
        out["balance"] = {
            "covariate_block_mean": round(float(np.mean([r["raw"] for r in comp])), 3),
            "covariate_window_mean": round(float(np.mean([r["window_raw"] / r["drawn"] for r in comp])), 3),
        }
    else:
        out["balance"] = {}
    out["balance"].update(
        elements_per_carrying_block=round(float(np.mean([r["elements"] for r in comp])), 3),
        elements_per_carrying_window=round(
            float(np.mean([r["window_elements"] / r["carrying"] for r in comp])), 3
        ),
    )
    out["_boots"] = boots  # dropped before saving; kept for the real-minus-neutral contrast
    return out


def reading(ci: list[float], point: float) -> dict[str, Any]:
    lo, hi = ci
    if lo > 0:
        key = "above"
    elif hi < 0:
        key = "below"
    else:
        key = "at_chance"
    return {
        "outcome": key,
        "clause_2": "passes as a labelled lead" if key == "above" else "stays not met",
        "below_the_chance_band": point < -CHANCE_BAND,
        "falsifier_fired": point < FALSIFIER_GAP,
        "text": PRE_REGISTRATION["readings"][key],
    }


# ---- the run ------------------------------------------------------------------------------


def target_sets(blocks: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    return {
        "real_unknown": [b for b in blocks if b["tier"] == "constrained_unknown" and not b["copy"]],
        "neutral": [b for b in blocks if b["tier"] == "neutral"],
    }


def collect(chroms: list[str]) -> dict[str, Any]:
    t0 = time.time()
    # pass 1: blocks and TSSs only, so the decile cuts exist before any window is drawn
    blocks: dict[str, list[dict[str, Any]]] = {}
    tss: dict[str, list[int]] = {}
    for c in chroms:
        blocks[c] = organise.blocks(c)
        tss[c] = coding_tss(c)
        print(f"{c}: {len(blocks[c])} blocks, {len(tss[c])} coding TSSs", flush=True)
    edges = {
        (name, cov): decile_edges(
            [fn(tss[c], (b["start"] + b["end"]) // 2) for c in chroms for b in target_sets(blocks[c])[name]]
        )
        for name in ("real_unknown", "neutral")
        for cov, fn in COVARIATES.items()
    }
    runs = {
        "unmatched_real_unknown": {},
        "unmatched_neutral": {},
        "matched_real_unknown": {},
        "matched_neutral": {},
        "sens_distance_real_unknown": {},
        "sens_exclude_targets_only_real_unknown": {},
    }
    qs = list(QUESTIONS)
    for c in chroms:
        els = cut.elements_of(c)
        els.sort(key=lambda e: e["start"])
        got = cut.read_chromosome(c, els)  # the organiser's blocks as d717b28 saw them
        if not got:
            continue
        sets = target_sets(got)
        t = tss[c]
        for name in ("real_unknown", "neutral"):
            runs[f"unmatched_{name}"][c] = matched_windows(
                sets[name], got, els, QUESTIONS, None, max_tries=UNMATCHED_TRIES
            )
            e = edges[(name, "tss_count")]
            runs[f"matched_{name}"][c] = matched_windows(
                sets[name],
                got,
                els,
                QUESTIONS,
                lambda m, e=e, t=t: bin_of(tss_count(t, m), e),
                raw_fn=lambda m, t=t: tss_count(t, m),
            )
        e = edges[("real_unknown", "tss_distance")]
        runs["sens_distance_real_unknown"][c] = matched_windows(
            sets["real_unknown"],
            got,
            els,
            QUESTIONS,
            lambda m, e=e, t=t: bin_of(tss_distance(t, m), e),
            raw_fn=lambda m, t=t: tss_distance(t, m),
        )
        e = edges[("real_unknown", "tss_count")]
        runs["sens_exclude_targets_only_real_unknown"][c] = matched_windows(
            sets["real_unknown"],
            sets["real_unknown"],
            els,
            QUESTIONS,
            lambda m, e=e, t=t: bin_of(tss_count(t, m), e),
            raw_fn=lambda m, t=t: tss_count(t, m),
        )
        del els, got
        print(f"{c}: drawn ({time.time() - t0:.0f} s)", flush=True)

    gate = {}
    for name in ("real_unknown", "neutral"):
        recs = [r for c in chroms for r in runs[f"unmatched_{name}"].get(c, [])]
        gate[name] = {"got": lifted_counts(recs, qs), "want": REPRODUCE[name]}
        gate[name]["same"] = gate[name]["got"] == gate[name]["want"]
    passed = all(g["same"] for g in gate.values())
    out: dict[str, Any] = {
        "result": RESULT,
        "chromosomes": chroms,
        "registration": PRE_REGISTRATION,
        "decile_edges": {f"{n}:{c}": v for (n, c), v in edges.items()},
        "gate_reproduces_d717b28": {"passed": passed, **gate},
    }
    if not passed:
        out["reading"] = "the reproduction gate failed: no matched figure is reported"
        out["seconds"] = round(time.time() - t0, 1)
        return out
    res: dict[str, Any] = {}
    for run, by_chrom in runs.items():
        res[run] = {q: summarise(by_chrom, q) for q in qs}
    # real unknown minus neutral, matched, by independent bootstraps of the two sets
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
        secondary_neutral=res["matched_neutral"],
        secondary_neutral_reading={
            q: reading(s["ci95_over_blocks"], s["matched_difference_points"])["outcome"]
            for q, s in res["matched_neutral"].items()
            if s.get("matched_difference_points") is not None
        },
        secondary_real_minus_neutral=contrast,
        secondary_real_unknown_moves_a_gene=res["matched_real_unknown"]["moves_a_gene"],
        sensitivities={
            "tss_distance_deciles": res["sens_distance_real_unknown"],
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
    from genomeos.genome.annotation import default_gencode

    inputs = []
    for c in chroms:
        inputs += [mf.input_entry(p, partition=None) for p in organise.inputs(c)]
        p = cut.ELEMENTS / f"{c}.json"
        if not any(i["path"] == str(p) for i in inputs):
            inputs.append(mf.input_entry(p, partition=None))
    gencode = default_gencode({chroms[0]}) if chroms else None
    if gencode is not None and Path(gencode).exists():
        inputs.append(mf.input_entry(Path(gencode), partition=None))
    return {
        "sources": [
            {
                "accession": "ENCODE SCREEN cCREs scored by AlphaGenome deletion (the all-element archive)",
                "version": "AlphaGenome as served during the 2026-09 all-element sweep (unpinned); "
                "pinned here by sha256",
                "path": str(cut.ELEMENTS),
            },
            {
                "accession": "GENCODE protein-coding genes (the TSS density covariate)",
                "version": "the repository's default GENCODE GFF3, pinned by sha256",
            },
            {
                "accession": "Zoonomia cactus241way phyloP, gnomAD Gnocchi, UCSC genomicSuperDups "
                "(the organiser's tiers, axes and copy flag)",
                "version": "as pinned by the organiser's inputs",
            },
            {
                "accession": "this repository, data/results/constrained_unknown_targets.json at d717b28",
                "version": "git d717b28, the unmatched control the gate reproduces",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "half_window_bp": HALF_WINDOW,
            "deciles": DECILES,
            "draws_per_block": DRAWS,
            "matched_tries": MATCHED_TRIES,
            "unmatched_tries": UNMATCHED_TRIES,
            "seed_per_chromosome": SEED,
            "bootstrap": BOOTSTRAP,
            "chance_band_points": CHANCE_BAND,
            "falsifier_gap_points": FALSIFIER_GAP,
            "primary_question": PRIMARY_QUESTION,
        },
        "exclusions": [
            "constrained_unknown blocks that are copies are out of the real unknown",
            "an element belongs to the block or window holding its midpoint",
            "windows overlapping any organiser block are rejected (primary); windows in another coding-TSS "
            "decile than their block are rejected; blocks with no accepted window are undrawable",
            "blocks with no scored element, or no accepted window carrying one, leave the comparison",
        ],
        "partitions": "n/a: arithmetic over model answers already on disk; nothing fitted, nothing held out",
    }


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
            f"n={p['n_blocks_compared']} blocks; {out['primary_reading']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
