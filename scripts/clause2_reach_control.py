# SPDX-License-Identifier: AGPL-3.0-or-later
"""Milestone 1.3 clause 2: the two descriptive controls lane-elemcount named and did not run.

    uv run python scripts/clause2_reach_control.py [--chroms chr21] [--no-save]

The clause is held as not met. The sharpest number behind that is per element and needs no matching:
9.18% of the 3,280 elements inside the real unknown name a coding gene against 44.04% of the 456,573
elements in their windows, a pooled -34.86 points, -28.77 as a mean over blocks (faeb0da). lane-elemcount
named two descriptive controls it did not run and called the second the sharper:

(a) class composition -- a block's cCREs may be dELS or CTCF-only where a window's include PLS and pELS;
(b) reach -- how many protein-coding TSSs lie inside the model's own input window centred on the element,
    because the model can only name a gene it can see. Earlier TSS matching (0c8b82d) was at block
    midpoint level, which is not the same quantity.

This script runs both, descriptively, on the draw that produced the -28.77: d717b28's unmatched windows,
reproduced to the digit by the same gate. It adds no matching of its own, so lane-elemcount's cost note
(chr2 alone 3,811 s, because window contents are evaluated on every try) does not apply: every per-element
quantity is computed once per element before any window is drawn, and a window costs one sum per element
it holds. Model requests: none; the stored all-element archive is read one chromosome at a time and the
per-element response cache is never opened.

The registration (`PRE_REGISTRATION`, and its dated section in docs/ATTRIBUTION.md) was committed before
any reach or class figure was computed.
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import importlib.util
import math
import random
import time
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

from genomeos import manifest as mf
from genomeos.results import save_result


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).resolve().parent / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


mc = _load("clause2_matched_control")  # the density-matched control: its draw, gate, bands and bootstraps
cut = mc.cut  # constrained_unknown_targets, which drew d717b28's windows

RESULT = "clause2_reach_control"


def scorer_window() -> int:
    """The scorer's own input length, read from the scorer rather than assumed.

    `AlphaGenomeAdapter._live_scorer` resizes the variant's reference interval to
    `dna_client.SEQUENCE_LENGTH_1MB` before every deletion score (genomeos/predict/alphagenome_adapter.py),
    and that is the window every element in the archive was read in. When the client is not installed the
    predict layer's own copy of the constant is used; a test pins the two together.
    """
    try:  # pragma: no cover - the client is not installed in every environment
        from alphagenome.models import dna_client  # type: ignore[import-not-found]

        return int(dna_client.SEQUENCE_LENGTH_1MB)
    except Exception:
        from genomeos.predict.splice_sites import WINDOW

        return int(WINDOW)


SCORER_WINDOW = scorer_window()  # 1,048,576 bp
HALF_WINDOW = SCORER_WINDOW // 2  # 524,288 bp either side of the element midpoint
DRAWS = mc.DRAWS  # 50 windows per block, as d717b28
SEED = mc.SEED  # 20260913 per chromosome, as d717b28, so the draw reproduces it to the digit
UNMATCHED_TRIES = mc.UNMATCHED_TRIES  # 4,000, d717b28's cap
MATCHED_TRIES = mc.MATCHED_TRIES  # 20,000, e1dbcf3's cap, used only by the density-matched sensitivity
BOOTSTRAP = mc.BOOTSTRAP  # 10,000
DECILES = mc.DECILES
REPRODUCE = mc.REPRODUCE  # d717b28's pooled counts, which the draw must give back
PRIMARY_QUESTION = mc.PRIMARY_QUESTION
QUESTIONS = mc.QUESTIONS
POOLED_GAP = -34.86  # faeb0da, per element, pooled: 9.18% against 44.04%
PER_BLOCK_GAP = -28.77  # faeb0da, per element, mean over the 531 compared blocks
HALF_POOLED_GAP = round(POOLED_GAP / 2, 2)  # -17.43: the half-the-gap convention of e1dbcf3 and da5764e
MIN_CELL = 30  # elements per arm for a stratum to enter a standardisation
REACH_STRATA: tuple[tuple[int, int | None], ...] = (
    (0, 0),
    (1, 1),
    (2, 2),
    (3, 4),
    (5, 6),
    (7, 9),
    (10, 14),
    (15, 24),
    (25, None),
)
CLASSES = ("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3", "unclassified")
PROMOTER_LIKE = ("PLS", "pELS", "DNase-H3K4me3")

PRE_REGISTRATION = {
    "registered": "2026-09-28",
    "lane": "lane-tssreach",
    "why_now": (
        "lane-elemcount (da5764e/faeb0da) closed the matched-control route and named two descriptive "
        "controls it did not run, calling reach the sharper. They are the last cheap things that could "
        "explain the per-element gap without measurement, so they are registered and run before anything "
        "is concluded about what clause 2's failure means. Nothing here restates, weakens or re-opens any "
        "existing claim: clause 2 stays not met on the ground d717b28, 0c8b82d and faeb0da put it on, "
        "whatever these two describe"
    ),
    "chosen_after_the_earlier_runs": (
        "both controls were chosen after seeing faeb0da, which is admissible because that run only named "
        "them as unmeasured descriptions of the same elements. The quantities, their cut points, the "
        "estimators, the intervals, the standardisation rule, the thresholds and every reading below are "
        "fixed and committed here before either is computed on any chromosome, and both arms of both "
        "comparisons are read from the draw that produced the number they are meant to explain"
    ),
    "targets": (
        "the real unknown: constrained_unknown blocks with the copies out, from organise.blocks on all 24 "
        "chromosomes (882 blocks, 531 carrying a scored element, 3,280 elements); the neutral tier as the "
        "secondary target set. An element belongs to the block or window holding its midpoint, as before"
    ),
    "control": (
        "d717b28's unmatched draw, unchanged and reproduced to the digit by the gate below: windows of the "
        "block's exact length placed uniformly inside the span of the chromosome's scored-element "
        "midpoints, rejected on overlap with any organiser block, 50 per block, 4,000 tries at most, one "
        "random number per try, seed 20260913 per chromosome. No matching of any kind is added, because "
        "the quantity to be described is the same in every window the -34.86 was computed on. The "
        "coding-TSS-density-matched draw of 0c8b82d (20,000 tries) is run once as a registered sensitivity"
    ),
    "quantity_b_reach": (
        "per element, on its midpoint m = (start + end) // 2, all three computed before any window is "
        "drawn and never recomputed: (1) reach_count, the number of GENCODE protein-coding TSSs in "
        f"[m - {HALF_WINDOW}, m + {HALF_WINDOW}), half-open -- {HALF_WINDOW} being half of the input "
        "length the scorer itself uses, read from the scorer and not assumed: "
        "AlphaGenomeAdapter._live_scorer resizes the variant's reference interval to "
        f"dna_client.SEQUENCE_LENGTH_1MB = {SCORER_WINDOW} bp before every deletion score, and every "
        "element in the archive was read in that window; (2) reach_distance, |m - nearest coding TSS| in "
        "bp, uncapped; (3) genes_in_window, the number of protein-coding genes whose annotated span "
        "intersects the same window, registered as the sensitivity that is closer to what the RNA-seq "
        "gene scorer actually enumerated, since a gene whose body reaches the window but whose TSS does "
        "not is still scored. TSSs are placed as unknown_scoring.annotate places them (start, or end - 1 "
        "on the minus strand), the same list 0c8b82d used"
    ),
    "quantity_a_class": (
        "per element, its ENCODE SCREEN cCRE class from data/results/ccres_<chrom>.bed.gz field 5 (PLS, "
        "pELS, dELS, CTCF-only, DNase-H3K4me3), joined by element id; an id absent from the registry file "
        "is counted as 'unclassified' and reported, never dropped"
    ),
    "comparison": (
        "elements inside real-unknown blocks against the elements of that block's own accepted windows, "
        "for every block that carries an element and has at least one accepted window carrying one -- the "
        "531 blocks the -28.77 was computed on, by faeb0da's own `compared` rule"
    ),
    "estimators": {
        "reach_difference": (
            "primary for (b): the per-block difference in mean reach_count, the mean over compared blocks "
            "of (mean reach_count of the block's own elements) minus (mean reach_count of the elements of "
            "its accepted carrying windows), in coding TSSs. The same statistic is reported for "
            "log10(1 + reach_distance) and for genes_in_window, and the pooled means, medians and "
            "reach_count = 0 shares are reported beside them"
        ),
        "reach_standardised_rate": (
            "the statistic that decides between the readings for (b): the per-element rate of "
            "names_a_coding_gene, directly standardised on reach. Elements of both arms are placed in the "
            "fixed strata of reach_count 0, 1, 2, 3-4, 5-6, 7-9, 10-14, 15-24, 25+ -- cut points fixed "
            "here, not read from the data. A stratum enters the standardisation when it holds at least "
            f"{MIN_CELL} elements in each arm; the block elements in every excluded stratum are reported "
            "with their share and their own rate, never silently dropped. The standardised difference is "
            "(block rate over the included strata) minus (sum over included strata of the block's share "
            "of those elements times the window rate in that stratum), in points. It answers: if the "
            "block's elements had been placed in windows of the same reach, how much of the gap would be "
            "left"
        ),
        "class_standardised_rate": (
            "the same direct standardisation for (a), on the six classes instead of the reach strata, "
            "with the same minimum cell and the same reporting of excluded classes; and the per-class "
            "per-element rate in each arm, so it can be seen whether the gap lives between classes or "
            "inside them. The per-block difference in the promoter-like share (PLS + pELS + "
            "DNase-H3K4me3) is reported as the composition statistic"
        ),
        "joint": (
            "class x reach stratum standardised jointly, same minimum cell, as a sensitivity: the two "
            "descriptions together"
        ),
    },
    "interval": (
        "95% percentile bootstrap over the compared blocks, 10,000 resamples, numpy default_rng(20260913), "
        "and a bootstrap over the 24 chromosomes (resampled with replacement, blocks pooled) beside it -- "
        "mc.summarise's two bootstraps, unchanged. For a standardised difference the resample is over "
        "blocks: the block stratum weights and the block rates are recomputed on each resample and the "
        "window per-stratum rates are held at their pooled values, because the window arm is the "
        "reference population being standardised to and has 100 times the elements"
    ),
    "gate": (
        "before any reach or class figure is read, the draw in this script must reproduce d717b28's pooled "
        "counts for both target sets to the digit (real unknown 44,100 windows drawn, 31,676 carrying, "
        "21,217 naming a coding gene, 531 and 158 blocks; neutral 131,600 / 88,317 / 57,212, 1,181 and "
        "340) and the unstandardised per-element pooled rates must come back as faeb0da wrote them "
        "(3,280 block elements, 301 naming, 456,573 window elements, 201,072 naming). If either fails, no "
        "reach or class figure is reported"
    ),
    "thresholds": (
        "the half-the-gap convention of e1dbcf3 and da5764e, on the pooled per-element scale this run "
        f"works on: a description explains the gap when standardising on it leaves the pooled difference "
        f"above {HALF_POOLED_GAP} points, i.e. closes more than half of the {POOLED_GAP} points, and does "
        f"not when the standardised difference stays at or below {HALF_POOLED_GAP}. 'Far fewer' for reach "
        "is fixed here as either the pooled mean reach_count inside the blocks being below half the "
        "window mean, or the share of block elements with reach_count = 0 exceeding the window share by "
        "more than 10 points. 'Comparable' is fixed as the per-block reach difference's 95% interval "
        "lying wholly inside +/-1.0 coding TSS of 0 and the two reach_count = 0 shares differing by less "
        "than 2 points"
    ),
    "readings": {
        "reach_short_and_explains": (
            "reach is far lower inside the blocks by the threshold above and standardising on it closes "
            "more than half the gap: the instrument was asked to name a gene it could not see, and the "
            "-28.77 (pooled -34.86) is partly an artefact of reach. What clause 2's failure means changes: "
            "it stops being a statement about the blocks' sequence and becomes a statement about where "
            "the blocks sit relative to coding genes, which is what the organiser selected them for. The "
            "clause still does not pass -- the standardised difference is the number to quote, and it is "
            "still negative -- but the reading 'these elements name genes far less often than sequence "
            "like them' must be replaced by 'these elements sit where the model can see fewer genes'"
        ),
        "reach_short_but_gap_survives": (
            "reach is far lower inside the blocks yet standardising on it keeps more than half the gap: "
            "the reach imbalance is real and must be stated, but it is not the explanation. The number to "
            "quote becomes the reach-standardised difference and the failure stays about the sequence"
        ),
        "reach_comparable": (
            "reach is comparable by the threshold above: the model could see as many coding genes from "
            "inside the blocks as from their windows, the instrument was not asked to do the impossible, "
            "and clause 2's failure is about the sequence, not the window. Nothing about the clause's "
            "status changes and the -28.77 stands as it was written"
        ),
        "reach_comparable_but_gap_closes": (
            "reach is comparable and yet standardising on it closes more than half the gap: that "
            "combination is incoherent on its face and is to be reported as a likely implementation "
            "error in this script, not as a finding about clause 2, until it is traced"
        ),
        "class_explains": (
            "the class-standardised difference closes more than half the gap: the blocks hold a different "
            "mix of cCRE classes and the per-element gap is largely composition. The classes and their "
            "per-class rates become the number to quote"
        ),
        "class_does_not_explain": (
            "the class-standardised difference keeps at least half the gap: inside every class with "
            "enough elements the block elements still name a coding gene far less often, so class "
            "composition is a real difference between the two element sets and not the cause of the gap"
        ),
    },
    "what_this_cannot_do": (
        "neither control can make clause 2 pass. Both are descriptions of the same elements, and the best "
        "either can do is move the -34.86 toward 0, which replaces 'below chance' with 'explained by "
        "where the blocks sit' -- not with 'attributed'. The clause's third part, scored against "
        "measurement, is untouched by every outcome here, and no measured arm is run"
    ),
    "cost": (
        "0 AlphaGenome requests; the per-element response cache is never opened; one chromosome's archive "
        "in memory at a time. No matching is added, so no window's contents are evaluated on a rejected "
        "try: lane-elemcount's 3,811 s chr2 cost came from exactly that and does not arise here"
    ),
}


# ---- per-element quantities ---------------------------------------------------------------


def coding_loci(chrom: str) -> tuple[list[int], list[int], list[int]]:
    """Coding TSSs, and the sorted starts and ends of the protein-coding gene spans, for one chromosome."""
    from genomeos.genome import Annotation
    from genomeos.genome.annotation import default_gencode

    ann = Annotation.from_gff3(default_gencode({chrom}), {chrom})
    genes = [g for g in ann.protein_coding() if g.locus.chrom == chrom]
    tss = sorted((g.locus.end - 1 if g.locus.strand.value == "-" else g.locus.start) for g in genes)
    starts = sorted(g.locus.start for g in genes)
    ends = sorted(g.locus.end for g in genes)
    return tss, starts, ends


def genes_in_window(starts: list[int], ends: list[int], mid: int, half: int = HALF_WINDOW) -> int:
    """Protein-coding genes whose span intersects [mid - half, mid + half)."""
    return bisect.bisect_left(starts, mid + half) - bisect.bisect_right(ends, mid - half)


def classes_of(chrom: str, results: Path = Path("data/results")) -> dict[str, str]:
    """cCRE class per element id, from the registry file the element layer already reads."""
    out: dict[str, str] = {}
    p = results / f"ccres_{chrom}.bed.gz"
    if not p.exists():
        return out
    with gzip.open(p, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) > 4:
                out[f[3]] = f[4]
    return out


def stratum_of(reach: int) -> str:
    for lo, hi in REACH_STRATA:
        if reach >= lo and (hi is None or reach <= hi):
            return f"{lo}" if hi == lo else (f"{lo}+" if hi is None else f"{lo}-{hi}")
    return "?"


def element_facts(chrom: str, els: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Reach, distance, window gene count, class and both answers, once per element."""
    tss, starts, ends = coding_loci(chrom)
    cls = classes_of(chrom)
    facts = []
    for e in els:
        mid = (e["start"] + e["end"]) // 2
        reach = mc.tss_count(tss, mid, HALF_WINDOW)
        facts.append(
            {
                "reach": reach,
                "stratum": stratum_of(reach),
                "distance": mc.tss_distance(tss, mid),
                "genes": genes_in_window(starts, ends, mid),
                "class": cls.get(e["id"], "unclassified"),
                "yes": {q: t(e) for q, t in QUESTIONS.items()},
            }
        )
    return facts


# ---- accumulation -------------------------------------------------------------------------


def _blank(questions: Sequence[str]) -> dict[str, Any]:
    return {
        "n": 0,
        "reach": 0,
        "logdist": 0.0,
        "genes": 0,
        "zero_reach": 0,
        "promoter_like": 0,
        "yes": dict.fromkeys(questions, 0),
        "by_stratum": {},
        "by_class": {},
        "by_joint": {},
    }


def _add(acc: dict[str, Any], f: dict[str, Any]) -> None:
    acc["n"] += 1
    acc["reach"] += f["reach"]
    acc["logdist"] += math.log10(1 + f["distance"])
    acc["genes"] += f["genes"]
    acc["zero_reach"] += f["reach"] == 0
    acc["promoter_like"] += f["class"] in PROMOTER_LIKE
    for q, v in f["yes"].items():
        acc["yes"][q] += v
    for book, key in (
        ("by_stratum", f["stratum"]),
        ("by_class", f["class"]),
        ("by_joint", f"{f['class']}|{f['stratum']}"),
    ):
        cell = acc[book].setdefault(key, {"n": 0, **dict.fromkeys(f["yes"], 0)})
        cell["n"] += 1
        for q, v in f["yes"].items():
            cell[q] += v


def _merge(into: dict[str, Any], other: dict[str, Any]) -> None:
    for k in ("n", "reach", "logdist", "genes", "zero_reach", "promoter_like"):
        into[k] += other[k]
    for q, v in other["yes"].items():
        into["yes"][q] += v
    for book in ("by_stratum", "by_class", "by_joint"):
        for key, cell in other[book].items():
            tgt = into[book].setdefault(key, {k: 0 for k in cell})
            for k, v in cell.items():
                tgt[k] = tgt.get(k, 0) + v


def reach_windows(
    targets: list[dict[str, Any]],
    exclude: list[dict[str, Any]],
    rows: list[dict[str, Any]],
    facts: list[dict[str, Any]],
    key_fn: Callable[[int], Any] | None = None,
    seed: int = SEED,
    draws: int = DRAWS,
    max_tries: int = UNMATCHED_TRIES,
) -> list[dict[str, Any]]:
    """`clause2_matched_control.matched_windows` step for step, accumulating the per-element facts.

    One `randrange` per try and the same overlap test in the same order, so with `key_fn=None` and
    `max_tries=4000` it draws d717b28's windows exactly; the reproduction gate checks that it did. The
    per-element quantities are looked up, never recomputed, so a window costs one addition per element.
    `key_fn` takes the window midpoint only (the coding-TSS band of 0c8b82d), so no rejected try ever
    touches a window's contents.
    """
    if not targets or not rows:
        return []
    mids = sorted(((r["start"] + r["end"]) // 2, i) for i, r in enumerate(rows))
    keys = [m for m, _ in mids]
    lo, hi = keys[0], keys[-1]
    unknown = sorted((b["start"], b["end"]) for b in exclude)
    u_starts = [s for s, _ in unknown]
    qs = list(QUESTIONS)

    def overlaps(s: int, t: int) -> bool:
        i = max(0, bisect.bisect_right(u_starts, s) - 1)
        while i < len(unknown) and unknown[i][0] < t:
            if unknown[i][1] > s:
                return True
            i += 1
        return False

    def inside(s: int, t: int) -> list[int]:
        return [i for _, i in mids[bisect.bisect_left(keys, s) : bisect.bisect_left(keys, t)]]

    rng = random.Random(seed)
    out = []
    for b in sorted(targets, key=lambda b: b["start"]):
        length = b["length"]
        bmid = (b["start"] + b["end"]) // 2
        sel_b = inside(b["start"], b["end"])
        want = key_fn(bmid) if key_fn else None
        block = _blank(qs)
        for i in sel_b:
            _add(block, facts[i])
        rec: dict[str, Any] = {
            "start": b["start"],
            "length": length,
            "elements": len(sel_b),
            "yes": {q: block["yes"][q] > 0 for q in qs},
            "block": block,
            "drawn": 0,
            "carrying": 0,
            "windows_yes": dict.fromkeys(qs, 0),
            "window": _blank(qs),
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
            if key_fn is not None and key_fn(s + length // 2) != want:
                continue
            got += 1
            sel = inside(s, s + length)
            if not sel:
                continue
            rec["carrying"] += 1
            seen = dict.fromkeys(qs, 0)
            for i in sel:
                f = facts[i]
                _add(rec["window"], f)
                for q, v in f["yes"].items():
                    seen[q] += v
            for q in qs:
                rec["windows_yes"][q] += seen[q] > 0
        rec["drawn"] = got
        out.append(rec)
    return out


def compared(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """faeb0da's rule: blocks that carry an element and have an accepted window carrying one."""
    return [r for r in records if r["elements"] and r["carrying"] and r["window"]["n"]]


# ---- statistics ---------------------------------------------------------------------------


def _boot(per: dict[str, Any], seed: int, n_boot: int) -> dict[str, Any]:
    """The two bootstraps of mc.summarise on one per-block quantity, keyed by chromosome."""
    import numpy as np

    chroms = [c for c in per if len(per[c])]
    if not chroms:
        return {}
    arrs = {c: np.asarray(per[c], dtype=float) for c in chroms}
    d = np.concatenate([arrs[c] for c in chroms])
    rng = np.random.default_rng(seed)
    boots = d[rng.integers(0, len(d), size=(n_boot, len(d)))].mean(axis=1)
    sums = np.array([arrs[c].sum() for c in chroms])
    ns = np.array([len(arrs[c]) for c in chroms])
    cidx = rng.integers(0, len(chroms), size=(n_boot, len(chroms)))
    cboots = sums[cidx].sum(axis=1) / ns[cidx].sum(axis=1)
    return {
        "mean": round(float(d.mean()), 4),
        "ci95_over_blocks": [round(float(x), 4) for x in np.percentile(boots, [2.5, 97.5])],
        "ci95_over_chromosomes": [round(float(x), 4) for x in np.percentile(cboots, [2.5, 97.5])],
        "n_blocks": int(len(d)),
    }


def per_block(
    by_chrom: dict[str, list[dict[str, Any]]],
    f: Callable[[dict[str, Any]], float],
    seed: int = SEED,
    n_boot: int = BOOTSTRAP,
) -> dict[str, Any]:
    """Block value minus its windows' value, per compared block, with both bootstraps."""
    return _boot({c: [f(r) for r in compared(rs)] for c, rs in by_chrom.items()}, seed, n_boot)


def _mean_diff(key: str) -> Callable[[dict[str, Any]], float]:
    return lambda r: r["block"][key] / r["block"]["n"] - r["window"][key] / r["window"]["n"]


def standardise(records: list[dict[str, Any]], book: str, q: str, min_cell: int = MIN_CELL) -> dict[str, Any]:
    """Direct standardisation of the block arm's rate onto the window arm's per-cell rates."""
    blk: dict[str, dict[str, int]] = {}
    win: dict[str, dict[str, int]] = {}
    for r in records:
        for src, tgt in ((r["block"], blk), (r["window"], win)):
            for key, cell in src[book].items():
                t = tgt.setdefault(key, {"n": 0, q: 0})
                t["n"] += cell["n"]
                t[q] += cell[q]
    included = sorted(k for k in blk if blk[k]["n"] >= min_cell and win.get(k, {}).get("n", 0) >= min_cell)
    bn = sum(blk[k]["n"] for k in included)
    excluded = sorted(set(blk) - set(included))
    wtot = sum(c["n"] for c in win.values())
    out: dict[str, Any] = {
        "book": book,
        "min_cell": min_cell,
        "cells_included": included,
        "cells_excluded": {
            k: {
                "block_elements": blk[k]["n"],
                "block_yes": blk[k][q],
                "window_elements": win.get(k, {}).get("n", 0),
            }
            for k in excluded
        },
        "block_elements_standardised": bn,
        "block_elements_total": sum(c["n"] for c in blk.values()),
        "per_cell": {
            k: {
                "block_elements": blk[k]["n"],
                "block_rate": round(blk[k][q] / blk[k]["n"], 4),
                "window_elements": win[k]["n"],
                "window_rate": round(win[k][q] / win[k]["n"], 4),
                "difference_in_points": round(100 * (blk[k][q] / blk[k]["n"] - win[k][q] / win[k]["n"]), 2),
                "block_share": round(blk[k]["n"] / bn, 4) if bn else None,
                "window_share": round(win[k]["n"] / wtot, 4) if wtot else None,
            }
            for k in included
        },
    }
    if not bn:
        out["standardised_difference_points"] = None
        return out
    brate = sum(blk[k][q] for k in included) / bn
    wrate = sum((blk[k]["n"] / bn) * (win[k][q] / win[k]["n"]) for k in included)
    out.update(
        block_rate=round(brate, 4),
        window_rate_standardised_to_the_blocks=round(wrate, 4),
        window_rate_crude=round(sum(c[q] for c in win.values()) / wtot, 4) if wtot else None,
        standardised_difference_points=round(100 * (brate - wrate), 2),
    )
    return out


def standardised_bootstrap(
    by_chrom: dict[str, list[dict[str, Any]]],
    book: str,
    q: str,
    min_cell: int = MIN_CELL,
    seed: int = SEED,
    n_boot: int = BOOTSTRAP,
) -> dict[str, Any]:
    """The block resample of the standardised difference, window per-cell rates held fixed."""
    import numpy as np

    recs = [r for rs in by_chrom.values() for r in compared(rs)]
    base = standardise(recs, book, q, min_cell)
    if base.get("standardised_difference_points") is None:
        return base
    included = base["cells_included"]
    win_rate = {k: base["per_cell"][k]["window_rate"] for k in included}
    chroms = [c for c in by_chrom if compared(by_chrom[c])]
    cells = {
        c: [
            (
                [r["block"][book].get(k, {"n": 0})["n"] for k in included],
                [r["block"][book].get(k, {q: 0}).get(q, 0) for k in included],
            )
            for r in compared(by_chrom[c])
        ]
        for c in chroms
    }
    ns = np.array([n for c in chroms for n, _ in cells[c]], dtype=float)
    ys = np.array([y for c in chroms for _, y in cells[c]], dtype=float)
    w = np.array([win_rate[k] for k in included], dtype=float)
    owner = np.concatenate([np.full(len(cells[c]), i) for i, c in enumerate(chroms)])

    def stat(rows: Any) -> Any:
        n = ns[rows].sum(axis=-2)
        y = ys[rows].sum(axis=-2)
        tot = n.sum(axis=-1)
        out = np.full(np.shape(tot), np.nan)
        ok = tot > 0
        out[ok] = 100 * (y.sum(axis=-1)[ok] / tot[ok] - (n[ok] / tot[ok, None]) @ w)
        return out

    rng = np.random.default_rng(seed)
    boots = stat(rng.integers(0, len(ns), size=(n_boot, len(ns))))
    cidx = rng.integers(0, len(chroms), size=(n_boot, len(chroms)))
    members = [np.where(owner == j)[0] for j in range(len(chroms))]
    cboots = np.array([float(stat(np.concatenate([members[j] for j in row]))) for row in cidx])
    base["ci95_over_blocks"] = [round(float(x), 2) for x in np.nanpercentile(boots, [2.5, 97.5])]
    base["ci95_over_chromosomes"] = [round(float(x), 2) for x in np.nanpercentile(cboots, [2.5, 97.5])]
    base["n_blocks"] = int(len(ns))
    base["closes_more_than_half_the_gap"] = base["standardised_difference_points"] > HALF_POOLED_GAP
    base["share_of_the_pooled_gap_closed"] = round(1 - base["standardised_difference_points"] / POOLED_GAP, 3)
    return base


def pooled_facts(by_chrom: dict[str, list[dict[str, Any]]], q: str) -> dict[str, Any]:
    """The crude per-element rates and the pooled descriptions of both arms."""
    recs = [r for rs in by_chrom.values() for r in compared(rs)]
    qs = list(QUESTIONS)
    tot = {"block": _blank(qs), "window": _blank(qs)}
    for r in recs:
        _merge(tot["block"], r["block"])
        _merge(tot["window"], r["window"])
    out: dict[str, Any] = {"n_blocks_compared": len(recs)}
    for arm, acc in tot.items():
        n = acc["n"] or 1
        out[arm] = {
            "elements": acc["n"],
            "names_a_coding_gene": acc["yes"][q],
            "rate": round(acc["yes"][q] / n, 4),
            "mean_reach": round(acc["reach"] / n, 3),
            "mean_log10_distance": round(acc["logdist"] / n, 3),
            "mean_genes_in_window": round(acc["genes"] / n, 3),
            "zero_reach_share": round(acc["zero_reach"] / n, 4),
            "promoter_like_share": round(acc["promoter_like"] / n, 4),
            "class_shares": {
                c: round(acc["by_class"][c]["n"] / n, 4) for c in CLASSES if acc["by_class"].get(c)
            },
            "reach_stratum_shares": {k: round(v["n"] / n, 4) for k, v in sorted(acc["by_stratum"].items())},
        }
    b, w = out["block"], out["window"]
    out["difference_in_points"] = round(100 * (b["rate"] - w["rate"]), 2)
    out["reach_ratio_block_over_window"] = (
        round(b["mean_reach"] / w["mean_reach"], 3) if w["mean_reach"] else None
    )
    out["zero_reach_share_difference_in_points"] = round(
        100 * (b["zero_reach_share"] - w["zero_reach_share"]), 2
    )
    return out


def reading(pooled: dict[str, Any], reach_diff: dict[str, Any], std: dict[str, Any]) -> dict[str, Any]:
    """The registered 2x2 for (b), decided by the thresholds fixed before the run."""
    ratio = pooled["reach_ratio_block_over_window"]
    zero_gap = pooled["zero_reach_share_difference_in_points"]
    short = (ratio is not None and ratio < 0.5) or zero_gap > 10
    lo, hi = reach_diff["ci95_over_blocks"]
    comparable = lo > -1.0 and hi < 1.0 and abs(zero_gap) < 2
    closes = bool(std.get("closes_more_than_half_the_gap"))
    if comparable and not short:
        key = "reach_comparable_but_gap_closes" if closes else "reach_comparable"
    else:
        key = "reach_short_and_explains" if closes else "reach_short_but_gap_survives"
    return {
        "outcome": key,
        "reach_is_far_lower_inside_the_blocks": short,
        "reach_is_comparable": comparable,
        "neither_threshold_met": not short and not comparable,
        "standardising_on_reach_closes_more_than_half": closes,
        "standardised_difference_points": std.get("standardised_difference_points"),
        "clause_2": "stays not met; no outcome here can make it pass",
        "text": PRE_REGISTRATION["readings"][key],
    }


# ---- the run ------------------------------------------------------------------------------


def collect(chroms: list[str]) -> dict[str, Any]:
    t0 = time.time()
    runs: dict[str, dict[str, list[dict[str, Any]]]] = {
        k: {} for k in ("real_unknown", "neutral", "sens_density_matched_real_unknown")
    }
    for c in chroms:
        els = cut.elements_of(c)
        els.sort(key=lambda e: e["start"])
        got = cut.read_chromosome(c, els)
        if not got:
            continue
        facts = element_facts(c, els)
        sets = mc.target_sets(got)
        for name in ("real_unknown", "neutral"):
            runs[name][c] = reach_windows(sets[name], got, els, facts)
        tss = mc.coding_tss(c)
        edges = mc.decile_edges(
            [mc.tss_count(tss, (b["start"] + b["end"]) // 2) for b in sets["real_unknown"]]
        )
        runs["sens_density_matched_real_unknown"][c] = reach_windows(
            sets["real_unknown"],
            got,
            els,
            facts,
            key_fn=lambda mid, t=tss, e=edges: mc.bin_of(mc.tss_count(t, mid), e),
            max_tries=MATCHED_TRIES,
        )
        print(f"{c}: {len(els)} elements described, windows drawn ({time.time() - t0:.0f} s)", flush=True)
        del els, got, facts

    qs = list(QUESTIONS)
    gate: dict[str, Any] = {}
    for name in ("real_unknown", "neutral"):
        recs = [r for c in chroms for r in runs[name].get(c, [])]
        gate[name] = {"got": mc.lifted_counts(recs, qs), "want": REPRODUCE[name]}
        gate[name]["same"] = gate[name]["got"] == gate[name]["want"]
    crude = pooled_facts(runs["real_unknown"], PRIMARY_QUESTION)
    gate["per_element_rates"] = {
        "got": {
            "block_elements": crude["block"]["elements"],
            "block_elements_yes": crude["block"]["names_a_coding_gene"],
            "window_elements": crude["window"]["elements"],
            "window_elements_yes": crude["window"]["names_a_coding_gene"],
        },
        "want": {
            "block_elements": 3280,
            "block_elements_yes": 301,
            "window_elements": 456573,
            "window_elements_yes": 201072,
        },
    }
    gate["per_element_rates"]["same"] = gate["per_element_rates"]["got"] == gate["per_element_rates"]["want"]
    passed = all(g["same"] for g in gate.values())
    out: dict[str, Any] = {
        "result": RESULT,
        "chromosomes": chroms,
        "registration": PRE_REGISTRATION,
        "scorer_input_window_bp": SCORER_WINDOW,
        "gate_reproduces_d717b28_and_faeb0da": {"passed": passed, **gate},
    }
    if not passed:
        out["reading"] = "the reproduction gate failed: no reach or class figure is reported"
        out["seconds"] = round(time.time() - t0, 1)
        return out

    reach_diff = per_block(runs["real_unknown"], _mean_diff("reach"))
    std_reach = standardised_bootstrap(runs["real_unknown"], "by_stratum", PRIMARY_QUESTION)
    std_class = standardised_bootstrap(runs["real_unknown"], "by_class", PRIMARY_QUESTION)
    out.update(
        pooled=crude,
        reach={
            "per_block_difference_in_coding_tss": reach_diff,
            "per_block_difference_in_log10_distance": per_block(runs["real_unknown"], _mean_diff("logdist")),
            "per_block_difference_in_genes_in_window": per_block(runs["real_unknown"], _mean_diff("genes")),
            "standardised": std_reach,
        },
        reach_reading=reading(crude, reach_diff, std_reach),
        classes={
            "per_block_difference_in_promoter_like_share": per_block(
                runs["real_unknown"], _mean_diff("promoter_like")
            ),
            "standardised": std_class,
            "reading": PRE_REGISTRATION["readings"][
                "class_explains"
                if std_class.get("closes_more_than_half_the_gap")
                else "class_does_not_explain"
            ],
        },
        sensitivities={
            "class_and_reach_jointly": standardised_bootstrap(
                runs["real_unknown"], "by_joint", PRIMARY_QUESTION
            ),
            "reach_on_density_matched_windows": {
                "pooled": pooled_facts(runs["sens_density_matched_real_unknown"], PRIMARY_QUESTION),
                "per_block_difference_in_coding_tss": per_block(
                    runs["sens_density_matched_real_unknown"], _mean_diff("reach")
                ),
                "standardised": standardised_bootstrap(
                    runs["sens_density_matched_real_unknown"], "by_stratum", PRIMARY_QUESTION
                ),
            },
            "neutral_tier": {
                "pooled": pooled_facts(runs["neutral"], PRIMARY_QUESTION),
                "per_block_difference_in_coding_tss": per_block(runs["neutral"], _mean_diff("reach")),
                "standardised": standardised_bootstrap(runs["neutral"], "by_stratum", PRIMARY_QUESTION),
            },
            "moves_a_gene": {
                "standardised_on_reach": standardised_bootstrap(
                    runs["real_unknown"], "by_stratum", "moves_a_gene"
                ),
                "standardised_on_class": standardised_bootstrap(
                    runs["real_unknown"], "by_class", "moves_a_gene"
                ),
            },
        },
        seconds=round(time.time() - t0, 1),
    )
    return out


@mf.depends_on_models("alphagenome")  # the all-element sweep's model, as far as the disk says (R9)
def manifest(chroms: list[str]) -> dict[str, Any]:
    m = mc.manifest(chroms)
    m["sources"] = [
        *m["sources"],
        {
            "accession": "ENCODE SCREEN cCRE registry classes, data/results/ccres_<chrom>.bed.gz",
            "version": "ENCODE cCREs v3 as distilled into this repository",
        },
        {
            "accession": "this repository, data/results/clause2_element_count_control.json at faeb0da",
            "version": "git faeb0da, the per-element gap these two controls describe",
        },
    ]
    m["parameters"] = {
        "scorer_input_window_bp": SCORER_WINDOW,
        "half_window_bp": HALF_WINDOW,
        "draws_per_block": DRAWS,
        "unmatched_tries": UNMATCHED_TRIES,
        "matched_tries": MATCHED_TRIES,
        "seed_per_chromosome": SEED,
        "bootstrap": BOOTSTRAP,
        "deciles": DECILES,
        "min_cell_elements": MIN_CELL,
        "reach_strata": [list(s) for s in REACH_STRATA],
        "primary_question": PRIMARY_QUESTION,
        "pooled_gap_points": POOLED_GAP,
        "half_pooled_gap_points": HALF_POOLED_GAP,
        "matched_on": "nothing: d717b28's unmatched draw, described",
    }
    m["exclusions"] = [
        "constrained_unknown blocks that are copies are out of the real unknown",
        "an element belongs to the block or window holding its midpoint",
        "windows overlapping any organiser block are rejected, as d717b28",
        "blocks with no scored element, or no accepted window carrying one, leave the comparison",
        f"a stratum or class holding fewer than {MIN_CELL} elements in either arm is out of the "
        "standardisation and is reported with its counts",
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
    print("gate passed:", out["gate_reproduces_d717b28_and_faeb0da"]["passed"])
    if "pooled" in out:
        p = out["pooled"]
        print(
            f"reach: block {p['block']['mean_reach']} coding TSSs in the {SCORER_WINDOW} bp window "
            f"against window {p['window']['mean_reach']}; zero-reach "
            f"{p['block']['zero_reach_share']} against {p['window']['zero_reach_share']}"
        )
        print(
            "per-block reach difference:",
            out["reach"]["per_block_difference_in_coding_tss"]["mean"],
            out["reach"]["per_block_difference_in_coding_tss"]["ci95_over_blocks"],
        )
        for name, blk in (
            ("reach", out["reach"]["standardised"]),
            ("class", out["classes"]["standardised"]),
        ):
            print(
                f"{name}-standardised per-element difference: "
                f"{blk['standardised_difference_points']} points, CI {blk.get('ci95_over_blocks')} "
                f"(crude {p['difference_in_points']})"
            )
        print("reading:", out["reach_reading"]["outcome"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
