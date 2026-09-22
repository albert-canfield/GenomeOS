# SPDX-License-Identifier: AGPL-3.0-or-later
"""A measured confidence for every predicted enhancer-to-gene target: the CRISPRi screens as a scale.

Three measurements of chromatin structure have failed to say which element acts on which gene
(CTCF orientation, measured boundary strength, measured Hi-C contact at 5 kb); the predicted
deletion works (docs/LESSONS.md, "Three measurements of chromatin structure, three nulls"). What
the project therefore carries genome-wide is the deletion's output, and the number beside it is
currently the model's own effect size: the sweep writes `confidence = min(1, |log2 fold change|)`.
That is a rescaled effect size, not a probability, and nothing measured says what it means.

This module turns the CRISPRi benchmark into that missing scale. On the K562 training pairs it fits
the probability that a screen calls a pair `Regulated` from the features the genome-wide sweep also
has, so the fit can be applied to all 961,227 deleted elements (612,323 of them with a target):

- `deletion_drop`: the predicted drop in K562 when the pair's gene is the element's top predicted
  target (the negated log2 fold change, floored at zero), zero otherwise;
- `top_target`: whether the pair's gene is that top target at all;
- `log_tss_distance`: the log distance from the element's midpoint to the gene's GENCODE v50 TSS,
  computed the same way on the benchmark and on the sweep, so the feature transfers;
- `node_target`: whether the gene is the nearest TSS inside the element's CTCF node;
- the registry class of the element (dELS is the reference level; pELS, PLS, CTCF-only,
  DNase-H3K4me3 and `unknown` are the levels fitted).

The benchmark's measured activity (DNase x H3K27ac in the screen's own cell) is deliberately left
out of the calibration: the project has cached DNase for one chromosome, so a calibration that used
it could not be applied to the sweep. What it would have bought is reported beside the result.

Two fits, because the population matters more than the fit does:

- the **tested-pair** calibration, on every K562 training pair on a deleted element whose gene
  GENCODE knows (8,796 of 9,237 pairs, 437 regulated, 5.0%). It answers: if a K562 CRISPRi screen
  tested this pair, how often is it called regulated. `top_target` is a feature here.
- the **predicted-target** calibration, on the 245 of those pairs where the tested gene *is* the
  element's top predicted target. That is the sweep's own population, one gene per element, and its
  base rate is nothing like the other's (188 of 245, 77%). Only three features are fitted on it
  (drop, distance, node), and 40 held-out pairs cannot judge it, which is said and not hidden.

`PREREGISTERED_CALIBRATION` was fixed in this file before any held-out pair was scored. The claim
is reliability, not ranking: the observed rate must sit inside the predicted bin's interval in a
stated majority of equal-count bins, and the expected calibration error and the Brier score must
both beat the effect size read as a confidence. AUPRC is reported because the benchmark reports it,
but a better AUPRC with a worse calibration would fail this test, as it should.

The scope is narrow and stated in the result itself: K562, on cCREs inside a CTCF node, near a gene
a screen chose to test. Everything else is extrapolation (`SCOPE`, `FALSIFIES_TRANSFER`).
"""

from __future__ import annotations

import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from genomeos.attribution import crispri, targets
from genomeos.attribution.element_types import gene_starts, registry_classes

ELEMENTS = crispri.ELEMENTS
RESULTS = Path("data/results")
REFERENCE = Path("data/reference")
CELL = "K562"  # the only cell line the calibration is fitted in, and the only one it speaks for
MIN_DISTANCE = crispri.MIN_DISTANCE
CLASS_REFERENCE = "dELS"  # the commonest class; the fitted levels are the others
CLASS_LEVELS = ("pELS", "PLS", "CTCF-only", "DNase-H3K4me3", "unknown")
SWEEP_FEATURES = (
    "deletion_drop",
    "top_target",
    "log_tss_distance",
    "node_target",
    *(f"class_{c}" for c in CLASS_LEVELS),
)
TARGET_FEATURES = ("deletion_drop", "log_tss_distance", "node_target")  # the sweep's own population
ACTIVITY_FEATURES = (*SWEEP_FEATURES, "log_activity")  # reported, never applied: DNase is not cached
BINS = 10  # equal-count bins for the reliability table
SMALL_BINS = 3  # what 40 held-out predicted-target pairs can carry
MIN_BINS_CONSISTENT = 7  # of BINS, the pre-registered majority
DROP_BANDS = ((0.0, 0.0), (0.0, 0.1), (0.1, 0.2), (0.2, 0.5), (0.5, 99.0))
CONFIDENCE_BANDS = (
    (0.0, 0.02),
    (0.02, 0.05),
    (0.05, 0.1),
    (0.1, 0.25),
    (0.25, 0.5),
    (0.5, 0.75),
    (0.75, 0.9),
    (0.9, 1.01),
)
CHROMS = (*(f"chr{i}" for i in range(1, 23)), "chrX", "chrY")
CLIP = 1e-6  # a probability of exactly zero makes the log loss infinite; the clip is reported
EVIDENCE = (
    "experimental: CRISPRi enhancer-gene screens, ENCODE benchmark (EngreitzLab/CRISPR_comparison, "
    "Gschwind et al. 2025), K562; predicted: AlphaGenome deletion per element; curated: ENCODE cCRE "
    "classes and GENCODE v50 gene starts; inferred: CTCF-only nodes"
)
PREREGISTERED_CALIBRATION = (
    "fitted on the K562 training pairs that lie on a deleted element, using only the features the "
    "genome-wide sweep also has (the predicted drop in K562, the top-target flag, the log distance to "
    "the gene's GENCODE TSS, whether the gene is the node's nearest TSS, and the registry class), the "
    "calibrated probability is reliable on the held-out K562 pairs and better calibrated than the "
    "model's own effect size read as a confidence: (a) in at least 7 of the 10 equal-count bins the "
    "bin's mean predicted probability lies inside the 95% Wilson interval of the bin's observed rate, "
    "and (b) both the expected calibration error and the Brier score are lower than those of "
    "min(1, |log2 fold change|) scored on the same pairs. The predicted-target calibration fitted on "
    "the top-target pairs alone is reported beside it and is not part of the claim: 40 held-out pairs "
    "cannot judge it"
)
CRITERION_NOTE = (
    "The Hosmer-Lemeshow statistic is reported for every reliability table but is deliberately not part "
    "of the pre-registered rule. It was read on the training pairs only, before any held-out pair was "
    "scored, and it rejects the fit's own in-sample curve at 8,796 pairs (chi2 22.5 on 8 df, p 0.004) "
    "while the leave-chromosome-out curve, which is the honest estimate, gives p 0.013 at an expected "
    "calibration error of 0.006. At this many pairs the statistic answers a question about counts of "
    "one and two positives per bin, not about whether a probability is usable, so it is a diagnostic "
    "here and the per-bin interval plus the two errors are the test"
)
SCOPE = (
    "The calibration is measured on one cell line and one kind of element: K562 CRISPRi pairs whose "
    "element is an ENCODE cCRE inside a CTCF node and whose gene a screen chose to test, which means "
    "a gene expressed in K562 within about a megabase. The probability is conditional on that test "
    "happening: it says how often a K562 screen calls such a pair regulated, not how often an element "
    "regulates a gene. Applying it to another cell line, to elements outside the registry, to genes "
    "not expressed in the line, or to distances the screens do not reach is extrapolation, and for a "
    "sweep element whose strongest predicted tissue is not K562 the drop entering the calibration is "
    "the K562 drop, so the band is a statement about K562 and about nothing else. The objection that "
    "AlphaGenome was trained on ENCODE tracks of these same cell lines is narrowed and not closed: E1's "
    "2,260 fresh pairs (commit 2a07c80) have units agreeing 637 of 1,096 (0.581) against matched nulls "
    "490 of 968 (0.506), +0.075 with an upper 95% bound of 0.111 at one-sided p 0.00037 — the weak band, "
    "since the pre-registration asked for 0.10 in size, and a small p does not convert a weak effect "
    "into the declared one. One cell line carries it: GM12878 +0.089 on 1,895 pairs, Jurkat -0.001 on "
    "365. So in lymphoblastoid cells the model weakly and detectably tracks a reporter assay's "
    "direction, and outside them that endpoint says nothing. The declared secondary did pass, and it is "
    "the first met prediction of that series: agreement is higher where DAP-G fine-maps the variant, "
    "+0.198 on 132 units against 90 controls (p 0.0018) against +0.057 for the rest, which is a reason "
    "to expect the transfer where the causal variant is known, not evidence that it happens. A "
    "calibration fitted on K562 CRISPRi is therefore a calibration against measurement in K562."
)
FALSIFIES_TRANSFER = (
    "The transfer fails if, on a CRISPRi screen in another cell line with an AlphaGenome line in the "
    "deletion table (GM12878, HepG2 or IMR-90 at the same coverage), the bands do not hold: the "
    "observed rate per band outside the band's interval in most populated bands, or an expected "
    "calibration error above that of the effect size read as a confidence. It also fails if a screen "
    "that tests elements chosen without regard to K562 activity (an unbiased tiling rather than a "
    "candidate list) finds the top-band rate far below the band, which is what the selection on "
    "testability would look like from outside. A third way to fail is narrower and already measurable: "
    "E1's met secondary says agreement is higher where the causal variant is fine-mapped, so if the "
    "bands hold no better on the fine-mapped subset of a new screen than off it, the reason to expect "
    "any transfer at all is gone. The prevalence shift this result found is a fourth: a screen whose "
    "base rate differs from these screens' 5% will need its own intercept, and a band quoted without "
    "the population's rate beside it is not a measurement of anything"
)
PREREGISTERED_GATE = (
    "Registered 2026-09-22, before the calibration was re-fitted through the sweep's own per-element "
    "response cache and before any re-fitted held-out pair was scored.\n"
    "\n"
    "WHAT IS GATED. `top_target` is the flag that the pair's measured gene is the single gene the "
    "compact all_elements table kept for that element. It gates this module in three places. (1) The "
    "magnitude: `deletion_drop` is read from the compact table, so it is a structural zero for every "
    "pair whose gene is not that one gene, and the zero is a property of the projection, not of what "
    "the sweep predicted. On the fitted population it is zero for 8,589 of 8,796 K562 training pairs "
    "(97.6%) and 1,680 of 1,715 held-out K562 pairs (97.9%). (2) The element: `matched_element` picks "
    "the overlapping element whose top predicted target is the pair's gene, and otherwise the element "
    "with the largest predicted magnitude for whatever other gene the table named, which is the "
    "element the class and distance features are then taken from. (3) The population: the "
    "predicted-target calibration is fitted on the 245 training pairs the gate admits (188 regulated, "
    "76.7%) and read on 40 held-out pairs (36 regulated, 90.0%), and it is that fit whose weights band "
    "all 612,323 sweep targets.\n"
    "\n"
    "WHAT THE GATE EXCLUDES. 8,551 of the 8,796 fitted training pairs (249 regulated, 2.91%) and 1,675 "
    "of the 1,715 held-out K562 pairs (76 regulated, 4.54%). Asked through `targets.ElementResponses`, "
    "which carries every gene in the scorer's 1 Mb window with a signed change on the cell's own "
    "track, 5,571 of the excluded training pairs (65.1%) and 1,112 of the excluded held-out pairs "
    "(66.4%) carry a number the sweep did predict; the remaining 2,980 and 563 are the named silence "
    "`the gene is not in the scorer's window at this element` and stay unanswerable. No excluded pair "
    "is a `not on this cell's own track` or a `not cached` silence. So the gate is two thirds a "
    "censoring and one third a real limit, and the calibration can be fitted with a magnitude that "
    "means something on 5,816 training pairs instead of 245.\n"
    "\n"
    "WHAT IS RE-FITTED AND WHAT IS NOT. The published 2026-09-17 curve is annotated, never rewritten: "
    "`score()` keeps its default of no cache and reproduces the published numbers exactly, and the "
    "re-fit is a second arm that passes the cache. `SWEEP_FEATURES` and `TARGET_FEATURES` do not "
    "change, because a feature the sweep does not have cannot be fitted; only what `deletion_drop` "
    "means changes, from `the compact table's entry for the one gene it kept` to `what the sweep "
    "predicted for this pair's own gene at this element`.\n"
    "\n"
    "THE DIRECTION EXPECTED, AND WHY. Reliability is expected NOT to improve. The 2026-09-17 claim "
    "failed on the population and not on the curve: 6 of 10 bins inside their intervals, all four "
    "failures in the same direction, mean predicted 0.0441 against an observed 0.0653, and a single "
    "log-odds shift of +0.679 restoring 9 of 10. An uncensored magnitude is a better feature; it is "
    "not an intercept, and it cannot move a prevalence that is a property of how a screen chose its "
    "pairs. The registered expectation is therefore 6 of 10 bins give or take one, with the prevalence "
    "ratio still near 1.31, and a rise in AUPRC as the only movement a better feature buys.\n"
    "\n"
    "WHAT A CONFIDENCE MEANS FOR A PAIR THE GATE WOULD HAVE EXCLUDED. This is the clause that decides "
    "whether the re-fit is worth anything. The predicted-target curve is fitted on 245 pairs whose "
    "base rate is 76.7% -- the model's most confident calls, the one gene per element it was surest "
    "of -- and the sweep quotes it for 612,323 targets. A user who now asks about any other gene in "
    "the window, which the reader makes askable and which is about fifty genes per element rather "
    "than one, would be quoted a curve fitted on a population whose base rate is twenty-six times "
    "theirs. That is the classic form of a calibration that looks reliable and is not. The measurement "
    "registered here is direct: score the shipped 245-pair curve on the held-out pairs the gate "
    "excluded but the sweep did answer, and report its mean predicted probability against their "
    "observed rate. WHAT WOULD SHOW THE FAILURE IS HAPPENING: the shipped curve's mean predicted "
    "probability on those pairs sits far above their observed rate -- a gap of the order of the 0.767 "
    "against 0.029 base-rate gap rather than of the 1.31 prevalence ratio already found -- and its "
    "bins fall outside their intervals in one direction. WHAT WOULD SHOW IT IS NOT: the shipped curve "
    "lands near their observed rate, which would mean the three features carry the population "
    "difference and the gate was only selecting on them. Either way the number is reported, and a "
    "confidence for an off-gate pair is quoted from a curve fitted on off-gate pairs or it is not "
    "quoted at all.\n"
    "\n"
    "IF RELIABILITY IMPROVES. An improvement is the outcome that would tempt a lane to stop checking, "
    "so the checks are fixed in advance and run whether it improves or not. A curve that predicts the "
    "base rate everywhere is trivially inside every equal-count bin, so a rise in bins_consistent is "
    "reported only beside (a) the width of the predicted range across the ten bins, which must not "
    "shrink, (b) AUPRC on the same held-out pairs, which must not fall, and (c) the prevalence ratio "
    "and the log-odds shift, which must have moved towards 1 and 0 for the improvement to be a "
    "calibration rather than a flattening. If bins_consistent reaches 7 or more while the predicted "
    "range narrows or AUPRC falls, the improvement is recorded as a flattening and the 2026-09-17 "
    "verdict of failed is not upgraded. The verdict is only upgraded if the prevalence gap itself "
    "closes, and nothing in this change acts on the intercept, so that is not expected."
)


# ------------------------------------------------------------------------------------------
# Features: what the benchmark and the sweep both have
# ------------------------------------------------------------------------------------------


def class_features(cls: str) -> dict[str, float]:
    """One-hot of the registry class, with dELS as the reference level."""
    return {f"class_{c}": 1.0 if c == cls else 0.0 for c in CLASS_LEVELS}


def matched_element(elements: list[dict[str, Any]], gene: str) -> dict[str, Any] | None:
    """The one overlapping element a pair's features come from.

    A pair can sit on more than one deleted element. The element whose top predicted target is the
    pair's gene is the one the deletion feature already speaks for; failing that, the one with the
    largest predicted magnitude in K562; failing that, the first.
    """
    if not elements:
        return None
    for e in elements:
        for key in ("predicted_coding", "predicted"):
            p = e.get(key)
            if p and p["gene"] == gene:
                return e
    return max(elements, key=lambda e: abs((e.get("predicted_by_cell") or {}).get(CELL) or 0.0))


def tss_distance(midpoint: int, tss: int | None) -> float | None:
    return None if tss is None else float(max(MIN_DISTANCE, abs(midpoint - tss)))


def add_features(
    pairs: list[crispri.Pair],
    table: crispri.DeletionTable,
    reference: Path = REFERENCE,
    results: Path = RESULTS,
) -> dict[str, Any]:
    """Add the transferable features to already-annotated pairs; count every value that is missing.

    A pair whose gene has no GENCODE v50 gene entry has no distance under the definition the sweep
    uses, and is counted per arm rather than given a distance of zero. `scored` marks the pairs that
    have every feature; nothing else is fitted or scored.
    """
    counts: Counter[str] = Counter()
    by_chrom: dict[str, list[crispri.Pair]] = defaultdict(list)
    for p in pairs:
        by_chrom[p.chrom].append(p)
    for chrom, rows in sorted(by_chrom.items()):
        starts = gene_starts(chrom, reference)
        classes = registry_classes(chrom, results) if (results / f"ccres_{chrom}.bed.gz").exists() else {}
        if not starts:
            counts["chromosomes_without_a_gencode_file"] += 1
        for p in rows:
            els = table.overlapping(p.chrom, p.start, p.end)
            el = matched_element(els, p.gene)
            cls = classes.get((el or {}).get("id", ""), ("unknown", False))[0]
            d = tss_distance(p.midpoint, starts.get(p.gene))
            arm = "regulated" if p.regulated else "not regulated"
            if el is None:
                counts[f"no_deleted_element_{arm}"] += 1
            if cls == "unknown":
                counts[f"no_registry_class_{arm}"] += 1
            if d is None:
                counts[f"no_gencode_tss_{arm}"] += 1
            p.features.update(class_features(cls))
            p.features["log_tss_distance"] = math.log(d) if d is not None else 0.0
            p.features["node_target"] = p.features["node_nearest"]
            p.features["scored"] = 1.0 if (el is not None and d is not None) else 0.0
            p.features["benchmark_distance_error"] = abs(d - p.distance) if d is not None else -1.0
            counts[f"scored_{arm}"] += int(bool(p.features["scored"]))
            counts[f"pairs_{arm}"] += 1
    return dict(counts)


def scored(pairs: list[crispri.Pair]) -> list[crispri.Pair]:
    return [p for p in pairs if p.covered and p.features.get("scored")]


# ------------------------------------------------------------------------------------------
# The gate: which pairs `top_target` admits, and what the window reader says about the rest
# ------------------------------------------------------------------------------------------


def gate_population(rows: list[crispri.Pair]) -> dict[str, Any]:
    """Split a set of pairs by the top-target gate and report each side's rate.

    The two sides are the whole point: the gate admits the model's most confident calls, and a
    calibration fitted on them alone is a curve fitted on a base rate that is not the population's.
    """

    def arm(rs: list[crispri.Pair]) -> dict[str, Any]:
        k = sum(p.regulated for p in rs)
        return {"pairs": len(rs), "regulated": k, "rate": round(k / len(rs), 4) if rs else None}

    return {
        "all": arm(rows),
        "admitted, the gene is the top target": arm([p for p in rows if p.features["top_target"]]),
        "excluded by the gate": arm([p for p in rows if not p.features["top_target"]]),
        "deletion_drop is zero": sum(1 for p in rows if p.features["deletion_drop"] == 0.0),
        "the sweep answered this gene": sum(1 for p in rows if p.features.get("deletion_answered")),
    }


def gate_silences(
    rows: list[crispri.Pair], table: crispri.DeletionTable, responses: targets.ElementResponses
) -> dict[str, Any]:
    """For each pair, what the sweep's per-element cache says about its own gene, or which silence.

    An excluded pair is either a censoring — the sweep predicted a change for this gene and the
    compact table dropped it — or a real limit, the gene lying outside the scorer's window. Only the
    second is a reason the calibration cannot speak. The pairs are walked chromosome by chromosome
    because the reader holds one chromosome's archive at a time.
    """
    out: dict[str, Counter[str]] = {"admitted": Counter(), "excluded": Counter()}
    ranks: list[int] = []
    window: list[int] = []
    by_chrom: dict[str, list[crispri.Pair]] = defaultdict(list)
    for p in rows:
        by_chrom[p.chrom].append(p)
    for chrom in sorted(by_chrom):
        for p in by_chrom[chrom]:
            els = table.overlapping(p.chrom, p.start, p.end)
            best: targets.Response | None = None
            for e in els:
                r = responses.response(p.chrom, e["id"], p.gene, p.cell)
                if best is None or (r.answered and (not best.answered or r.value < best.value)):
                    best = r
            side = "admitted" if p.features["top_target"] else "excluded"
            out[side][best.reason if best is not None else "no overlapping deleted element"] += 1
            if els:
                order = [g for g, _ in responses.ranked(p.chrom, els[0]["id"], p.cell)]
                if order:
                    window.append(len(order))
                    if p.gene in order:
                        ranks.append(order.index(p.gene) + 1)
    return {
        "by_side": {k: dict(sorted(v.items())) for k, v in out.items()},
        "measured_gene_rank_in_the_window": {
            "pairs": len(ranks),
            "median_rank": sorted(ranks)[len(ranks) // 2] if ranks else None,
            "rank_1": sum(1 for r in ranks if r == 1),
            "top_3": sum(1 for r in ranks if r <= 3),
            "top_5": sum(1 for r in ranks if r <= 5),
            "median_genes_in_the_window": sorted(window)[len(window) // 2] if window else None,
        },
    }


def distance_agreement(pairs: list[crispri.Pair]) -> dict[str, Any]:
    """How far the GENCODE gene-level TSS sits from the distance the benchmark itself reports.

    The two disagree wherever a gene has several TSSs and the screen used another one. It is a
    property of the feature, not of the model, so it is measured once and reported.
    """
    errs = sorted(
        p.features["benchmark_distance_error"]
        for p in pairs
        if p.features.get("benchmark_distance_error", -1.0) >= 0
    )
    if not errs:
        return {"pairs": 0}
    return {
        "pairs": len(errs),
        "median_bp": round(errs[len(errs) // 2]),
        "p90_bp": round(errs[int(0.9 * (len(errs) - 1))]),
        "within_1kb": round(sum(e <= 1_000 for e in errs) / len(errs), 4),
        "within_10kb": round(sum(e <= 10_000 for e in errs) / len(errs), 4),
    }


# ------------------------------------------------------------------------------------------
# Calibration: the fit, the reliability table, the errors
# ------------------------------------------------------------------------------------------


def sigmoid(z: float) -> float:
    return 1 / (1 + math.exp(-max(-30.0, min(30.0, z))))


def fit(pairs: list[crispri.Pair], cols: tuple[str, ...], lam: float = 1e-3) -> list[float]:
    return crispri.logistic_fit(crispri.matrix(pairs, cols), [p.regulated for p in pairs], lam=lam)


def predict(w: list[float], pairs: list[crispri.Pair], cols: tuple[str, ...]) -> list[float]:
    return [sigmoid(z) for z in crispri.logistic_score(w, crispri.matrix(pairs, cols))]


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """95% Wilson interval for a rate; the interval the reliability table is judged against."""
    if not n:
        return (0.0, 1.0)
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def stratum(p: crispri.Pair) -> str:
    """A coarse stratum every covered pair has, scored or not: the top-target flag and the drop band.

    A calibration bin is a rate conditioned on the pair having been scored at all, so each bin
    inherits the coverage of both arms. The strata are built from the two features that exist for
    every covered pair — the flag and the predicted drop — and never from the distance or the class,
    which are the features a pair can be missing. A bin can then say which strata it draws from and
    how completely each arm of those strata was scored, before it says anything about its own rate.
    """
    d = p.features["deletion_drop"]
    band = next(
        (f"{lo:g}-{hi:g}" for lo, hi in DROP_BANDS[1:] if lo < d <= hi),
        "0" if d == 0.0 else f">{DROP_BANDS[-1][0]:g}",
    )
    return f"top={p.features['top_target']:g}|drop={band}"


def stratum_coverage(pairs: list[crispri.Pair]) -> dict[str, dict[str, list[int]]]:
    """Per stratum and per arm, how many covered pairs have every feature out of how many there are."""
    out: dict[str, dict[str, list[int]]] = defaultdict(lambda: {"regulated": [0, 0], "not regulated": [0, 0]})
    for p in pairs:
        if not p.covered:
            continue
        cell = out[stratum(p)]["regulated" if p.regulated else "not regulated"]
        cell[0] += int(bool(p.features.get("scored")))
        cell[1] += 1
    return dict(out)


def bin_coverage(strata: list[str], coverage: dict[str, dict[str, list[int]]]) -> dict[str, Any]:
    """The coverage of both arms in the strata a bin draws from, to be read before the bin's rate."""
    out: dict[str, Any] = {"strata": len(set(strata))}
    for arm in ("regulated", "not regulated"):
        ok = sum(coverage.get(s, {}).get(arm, [0, 0])[0] for s in set(strata))
        total = sum(coverage.get(s, {}).get(arm, [0, 0])[1] for s in set(strata))
        out[f"coverage_{arm.replace(' ', '_')}"] = [ok, total, round(ok / total, 4) if total else None]
    return out


def chi2_sf(x: float, df: int) -> float:
    """Upper tail of the chi-square distribution, standard library only (the recurrence in df)."""
    if x <= 0:
        return 1.0
    sf = math.erfc(math.sqrt(x / 2)) if df % 2 else math.exp(-x / 2)
    k = 1 if df % 2 else 2
    while k + 2 <= df:
        sf += (x / 2) ** (k / 2) * math.exp(-x / 2) / math.gamma(k / 2 + 1)
        k += 2
    return min(1.0, max(0.0, sf))


def hosmer_lemeshow(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Goodness of fit over the reliability bins: is the whole curve consistent with the counts.

    The per-bin interval answers one bin at a time and is easily failed by one noisy count; this is
    the textbook statistic over all of them at once, with two degrees of freedom spent on the fit.
    """
    chi = 0.0
    for r in rows:
        n, p = r["pairs"], r["mean_predicted"]
        if 0 < p < 1:
            chi += (r["regulated"] - n * p) ** 2 / (n * p * (1 - p))
    df = max(1, len(rows) - 2)
    return {"chi2": round(chi, 3), "df": df, "p": round(chi2_sf(chi, df), 4)}


def equal_count_bins(p: list[float], bins: int) -> list[list[int]]:
    """Indices grouped into `bins` nearly equal groups by predicted probability."""
    order = sorted(range(len(p)), key=lambda i: p[i])
    if not order:
        return []
    n, out = len(order), []
    for b in range(bins):
        lo, hi = b * n // bins, (b + 1) * n // bins
        if hi > lo:
            out.append(order[lo:hi])
    return out


def reliability(
    p: list[float],
    labels: list[bool],
    bins: int = BINS,
    strata: list[str] | None = None,
    coverage_per_stratum: dict[str, dict[str, list[int]]] | None = None,
) -> dict[str, Any]:
    """Predicted probability against observed rate per bin, with the interval each bin is judged by.

    Each row states the coverage of both arms in the strata it draws from before it states its own
    rate: a bin whose rate rises because only well-covered elements landed in it has to be visible
    as a number, not smoothed into the curve.
    """
    rows, ece, mce = [], 0.0, 0.0
    groups = equal_count_bins(p, bins)
    for g in groups:
        k = sum(labels[i] for i in g)
        mean_p = sum(p[i] for i in g) / len(g)
        obs = k / len(g)
        lo, hi = wilson(k, len(g))
        gap = abs(mean_p - obs)
        ece += gap * len(g) / len(p)
        mce = max(mce, gap)
        cov = (
            bin_coverage([strata[i] for i in g], coverage_per_stratum)
            if strata is not None and coverage_per_stratum is not None
            else {}
        )
        rows.append(
            {
                "pairs": len(g),
                **cov,
                "predicted_range": [round(min(p[i] for i in g), 4), round(max(p[i] for i in g), 4)],
                "mean_predicted": round(mean_p, 4),
                "regulated": k,
                "observed": round(obs, 4),
                "ci95": [round(lo, 4), round(hi, 4)],
                "consistent": bool(lo <= mean_p <= hi),
            }
        )
    return {
        "bins": len(rows),
        "bins_consistent": sum(r["consistent"] for r in rows),
        "ece": round(ece, 4),
        "mce": round(mce, 4),
        "hosmer_lemeshow": hosmer_lemeshow(rows),
        "rows": rows,
    }


def calibration(
    p: list[float],
    labels: list[bool],
    bins: int = BINS,
    strata: list[str] | None = None,
    coverage_per_stratum: dict[str, dict[str, list[int]]] | None = None,
) -> dict[str, Any]:
    """Every number one scorer earns on one set of pairs: coverage first, calibration, ranking last."""
    n = len(labels)
    if not n:
        return {"pairs": 0}
    brier = sum((pi - (1.0 if y else 0.0)) ** 2 for pi, y in zip(p, labels, strict=True)) / n
    clipped = [min(1 - CLIP, max(CLIP, pi)) for pi in p]
    logloss = -sum(math.log(pi) if y else math.log(1 - pi) for pi, y in zip(clipped, labels, strict=True)) / n
    return {
        **crispri.metrics(p, labels),
        "mean_predicted": round(sum(p) / n, 4),
        "brier": round(brier, 5),
        "log_loss": round(logloss, 4),
        "log_loss_clip": CLIP,
        "reliability": reliability(p, labels, bins, strata, coverage_per_stratum),
    }


def logit(p: float) -> float:
    q = min(1 - CLIP, max(CLIP, p))
    return math.log(q / (1 - q))


def prevalence_shift(p: list[float], labels: list[bool], bins: int = BINS) -> dict[str, Any]:
    """What the same probabilities look like once one constant is added to match the set's own rate.

    This is a description of a failure, not a test, and it is fitted on the very labels it is read
    against: the one number it adds is the shift in log odds that makes the mean predicted probability
    equal the observed rate. If the curve is right in shape and wrong in level — the signature of a
    prevalence that moved between the screens the fit saw and the screens it is applied to — the
    shifted curve becomes reliable while the shipped one does not. It never changes the verdict.
    """
    n = len(labels)
    if not n:
        return {"pairs": 0}
    observed = sum(labels) / n
    lo, hi = -12.0, 12.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if sum(sigmoid(logit(x) + mid) for x in p) / n < observed:
            lo = mid
        else:
            hi = mid
    shift = (lo + hi) / 2
    moved = [sigmoid(logit(x) + shift) for x in p]
    out = calibration(moved, labels, bins)
    return {
        "note": (
            "fitted on the held-out labels' own prevalence: a description of the failure, not a test, "
            "and not the calibration that is shipped"
        ),
        "mean_predicted_before": round(sum(p) / n, 4),
        "observed_rate": round(observed, 4),
        "log_odds_shift": round(shift, 4),
        "bins_consistent": out["reliability"]["bins_consistent"],
        "bins": out["reliability"]["bins"],
        "ece": out["reliability"]["ece"],
        "brier": out["brier"],
        "hosmer_lemeshow": out["reliability"]["hosmer_lemeshow"],
    }


def raw_confidence(pairs: list[crispri.Pair]) -> list[float]:
    """The number the sweep writes today: the effect size read as a confidence, min(1, |log2 fc|)."""
    return [min(1.0, p.features["deletion_drop"]) for p in pairs]


def scorers(pairs: list[crispri.Pair], weights: dict[str, list[float]]) -> dict[str, list[float]]:
    """The three numbers compared on the same pairs: the calibration, the drop fitted alone, today's."""
    return {
        "calibrated, sweep features": predict(weights["sweep"], pairs, SWEEP_FEATURES),
        "drop alone, fitted": predict(weights["drop"], pairs, ("deletion_drop",)),
        "effect size as a confidence": raw_confidence(pairs),
    }


def judge(held: dict[str, Any]) -> tuple[bool, dict[str, Any]]:
    """The pre-registered rule, read off the held-out numbers: reliability and two errors."""
    cal = held["calibrated, sweep features"]
    raw = held["effect size as a confidence"]
    checks = {
        "bins_consistent": cal["reliability"]["bins_consistent"],
        "bins_required": MIN_BINS_CONSISTENT,
        "reliable": cal["reliability"]["bins_consistent"] >= MIN_BINS_CONSISTENT,
        "ece": [cal["reliability"]["ece"], raw["reliability"]["ece"]],
        "ece_better": cal["reliability"]["ece"] < raw["reliability"]["ece"],
        "brier": [cal["brier"], raw["brier"]],
        "brier_better": cal["brier"] < raw["brier"],
    }
    return (checks["reliable"] and checks["ece_better"] and checks["brier_better"]), checks


# ------------------------------------------------------------------------------------------
# The measured table: the sweep's own population, with no model between it and the screens
# ------------------------------------------------------------------------------------------


def drop_band_label(drop: float) -> str:
    """The band a predicted drop falls in, spelled the same way in the measured table and in the sweep."""
    for lo, hi in DROP_BANDS:
        if drop == 0.0 if hi == 0.0 else lo < drop <= hi:
            return "= 0" if hi == 0.0 else f"{lo:g} < drop <= {'inf' if hi > 1 else f'{hi:g}'}"
    return "= 0"


DROP_LABELS = tuple(
    drop_band_label(d) for d in (0.0, 0.05, 0.15, 0.3, 1.0)
)  # the bands in order, for the sweep's counters


def by_drop_band(pairs: list[crispri.Pair]) -> list[dict[str, Any]]:
    """The observed rate per band of predicted drop, on every covered pair whose gene is the top target.

    This is the calibration a reader can check by hand: no fit, one interval per band, and the only
    feature it needs is the drop, which every covered pair has — so the band is not a rate conditioned
    on a pair having been scorable. Each band still reports how much of each arm carries the model's
    other features, before its rate, because that is the subset the fitted rows above are read on.
    """
    out = []
    for label in DROP_LABELS:
        rows = [p for p in pairs if drop_band_label(p.features["deletion_drop"]) == label]
        k = sum(p.regulated for p in rows)
        w = wilson(k, len(rows))
        row: dict[str, Any] = {"drop": label, "pairs": len(rows)}
        for arm, flag in (("regulated", True), ("not_regulated", False)):
            arm_rows = [p for p in rows if p.regulated is flag]
            ok = sum(bool(p.features.get("scored")) for p in arm_rows)
            row[f"coverage_{arm}"] = [
                ok,
                len(arm_rows),
                round(ok / len(arm_rows), 4) if arm_rows else None,
            ]
        row.update(
            {
                "regulated": k,
                "rate": round(k / len(rows), 4) if rows else None,
                "ci95": [round(w[0], 4), round(w[1], 4)] if rows else None,
                "elements": len({p.element for p in rows}),
            }
        )
        out.append(row)
    return out


def baselines(pairs: list[crispri.Pair]) -> dict[str, Any]:
    """The two distance baselines, so a subset can be checked for being an easier problem."""
    lab = [p.regulated for p in pairs]
    return {
        "distance": crispri.metrics([-p.features["log_distance"] for p in pairs], lab),
        "activity over distance": crispri.metrics([p.features["activity_over_distance"] for p in pairs], lab),
    }


def coverage(pairs: list[crispri.Pair], counts: dict[str, Any]) -> dict[str, Any]:
    """Coverage of the features per arm, and whether the scored subset is the same problem."""
    covered = [p for p in pairs if p.covered]
    arms = {}
    for arm, flag in (("regulated", True), ("not regulated", False)):
        rows = [p for p in covered if p.regulated is flag]
        ok = sum(bool(p.features.get("scored")) for p in rows)
        arms[arm] = {
            "pairs_on_a_deleted_element": len(rows),
            "with_every_feature": ok,
            "fraction": round(ok / len(rows), 4) if rows else None,
        }
    return {
        "arms": arms,
        "missing_values": counts,
        "baselines": {
            "on a deleted element": baselines(covered),
            "and with every feature": baselines(scored(pairs)),
        },
        "gencode_against_benchmark_distance": distance_agreement(covered),
    }


# ------------------------------------------------------------------------------------------
# The whole measurement
# ------------------------------------------------------------------------------------------


def score(
    training: list[crispri.Pair],
    heldout: list[crispri.Pair],
    table: crispri.DeletionTable,
    reference: Path = REFERENCE,
    results: Path = RESULTS,
) -> dict[str, Any]:
    """Fit on K562 training pairs, score the held-out K562 pairs once, judge the pre-registration."""
    crispri.annotate(training, table)
    crispri.annotate(heldout, table)
    found = {
        "training": add_features(training, table, reference, results),
        "heldout": add_features(heldout, table, reference, results),
    }
    train = scored(training)
    labels = [p.regulated for p in train]
    train_top = [p for p in train if p.features["top_target"]]
    train_cov = stratum_coverage(training)
    held_cov = stratum_coverage([p for p in heldout if p.cell == CELL])
    train_strata = [stratum(p) for p in train]

    weights = {
        "sweep": fit(train, SWEEP_FEATURES),
        "drop": fit(train, ("deletion_drop",)),
        "activity": fit(train, ACTIVITY_FEATURES),
        "target": fit(train_top, TARGET_FEATURES),
    }

    loco: dict[str, list[float]] = {}
    chroms = sorted({p.chrom for p in train})
    for name, cols in (("sweep", SWEEP_FEATURES), ("activity", ACTIVITY_FEATURES)):
        if len(chroms) < 2:  # one chromosome: there is nothing to leave out, and no fold is reported
            continue
        s = [0.0] * len(train)
        for c in chroms:
            rows = [p for p in train if p.chrom != c]
            w = crispri.logistic_fit(crispri.matrix(rows, cols), [p.regulated for p in rows])
            idx = [i for i, p in enumerate(train) if p.chrom == c]
            for i, v in zip(idx, predict(w, [train[i] for i in idx], cols), strict=True):
                s[i] = v
        loco[name] = s

    held_k562 = scored([p for p in heldout if p.cell == CELL])
    held_labels = [p.regulated for p in held_k562]
    held_strata = [stratum(p) for p in held_k562]
    held = {
        name: calibration(s, held_labels, BINS, held_strata, held_cov)
        for name, s in scorers(held_k562, weights).items()
    }
    passed, checks = judge(held)

    held_top = [p for p in held_k562 if p.features["top_target"]]
    train_top_all = [p for p in training if p.covered and p.features.get("top_target")]
    held_top_all = [p for p in heldout if p.cell == CELL and p.covered and p.features.get("top_target")]
    other_cells: dict[str, Any] = {}
    for cell in sorted({p.cell for p in heldout}):
        if cell == CELL:
            continue
        rows = scored([p for p in heldout if p.cell == cell])
        if cell not in crispri.MODEL_CELLS:
            other_cells[cell] = {"refused": "no AlphaGenome line for this cell type in the deletion table"}
        elif len(rows) < 100 or sum(p.regulated for p in rows) < 20:
            other_cells[cell] = {
                "refused": "too few pairs to read a reliability table",
                "pairs": len(rows),
                "regulated": sum(p.regulated for p in rows),
                "calibrated": calibration(
                    predict(weights["sweep"], rows, SWEEP_FEATURES), [p.regulated for p in rows], SMALL_BINS
                )
                if rows
                else {"pairs": 0},
            }
        else:
            other_cells[cell] = calibration(
                predict(weights["sweep"], rows, SWEEP_FEATURES), [p.regulated for p in rows]
            )

    return {
        "evidence": EVIDENCE,
        "preregistered": PREREGISTERED_CALIBRATION,
        "verdict": "passed" if passed else "failed",
        "checks": checks,
        "criterion_note": CRITERION_NOTE,
        "scope": SCOPE,
        "falsifies_the_transfer": FALSIFIES_TRANSFER,
        "populations": {
            "tested pair": {
                "what": "a K562 pair a screen tested, on a cCRE inside a CTCF node, every feature present",
                "training_pairs": len(train),
                "training_regulated": sum(labels),
                "training_rate": round(sum(labels) / len(train), 4) if train else None,
                "heldout_pairs": len(held_k562),
                "heldout_regulated": sum(held_labels),
                "heldout_rate": round(sum(held_labels) / len(held_k562), 4) if held_k562 else None,
            },
            "predicted target": {
                "what": "the same, and the tested gene is the element's top predicted target",
                "training_pairs": len(train_top),
                "training_regulated": sum(p.regulated for p in train_top),
                "training_rate": round(sum(p.regulated for p in train_top) / len(train_top), 4)
                if train_top
                else None,
                "heldout_pairs": len(held_top),
                "heldout_regulated": sum(p.regulated for p in held_top),
                "heldout_rate": round(sum(p.regulated for p in held_top) / len(held_top), 4)
                if held_top
                else None,
                "covered_pairs_before_any_feature_is_required": [len(train_top_all), len(held_top_all)],
            },
        },
        "coverage": {
            "training": coverage(training, found["training"]),
            "heldout, every cell type": coverage(heldout, found["heldout"]),
            "heldout K562, the set the claim is judged on": coverage(
                [p for p in heldout if p.cell == CELL], {}
            ),
        },
        "coverage_by_arm_of_the_deletion": {
            "training": crispri.coverage_by_arm(training),
            "heldout": crispri.coverage_by_arm(heldout),
        },
        "training_leave_chromosome_out": {
            name: calibration(s, labels, BINS, train_strata, train_cov) for name, s in loco.items()
        }
        or {"refused": "fewer than two chromosomes in the training pairs"},
        "training_in_sample": {
            name: calibration(s, labels, BINS, train_strata, train_cov)
            for name, s in scorers(train, weights).items()
        },
        "heldout_k562": held,
        "heldout_prevalence_shift": {
            name: prevalence_shift(s, held_labels) for name, s in scorers(held_k562, weights).items()
        },
        "heldout_other_cells": other_cells,
        "weights": {
            name: dict(
                zip(
                    ("intercept", *cols),
                    (round(v, 4) for v in weights[key]),
                    strict=True,
                )
            )
            for name, key, cols in (
                ("calibrated, sweep features", "sweep", SWEEP_FEATURES),
                ("drop alone, fitted", "drop", ("deletion_drop",)),
                ("with measured activity (not applied)", "activity", ACTIVITY_FEATURES),
                ("predicted-target population", "target", TARGET_FEATURES),
            )
        },
        "measured_by_drop_band": {
            "population": (
                "every covered K562 pair whose tested gene is the element's top predicted target (the "
                "sweep's own pair); no feature but the drop enters, so no pair is dropped for a missing "
                "value, and each band reports the coverage of both arms before its rate"
            ),
            "K562 training": by_drop_band(train_top_all),
            "K562 held out": by_drop_band(held_top_all),
            "K562 pooled (both, after the held-out test)": by_drop_band(train_top_all + held_top_all),
        },
        "predicted_target_calibration": {
            "note": (
                "fitted on the training pairs of the sweep's own population; the held-out reliability "
                "rests on a few dozen pairs and is not part of the pre-registration"
            ),
            "training_in_sample": calibration(
                predict(weights["target"], train_top, TARGET_FEATURES),
                [p.regulated for p in train_top],
                SMALL_BINS,
                [stratum(p) for p in train_top],
                train_cov,
            ),
            "heldout": calibration(
                predict(weights["target"], held_top, TARGET_FEATURES),
                [p.regulated for p in held_top],
                SMALL_BINS,
                [stratum(p) for p in held_top],
                held_cov,
            )
            if held_top
            else {"pairs": 0},
        },
        "what_measured_activity_would_add": {
            "note": (
                "DNase x H3K27ac measured in the screen's own cell is in the benchmark's columns but "
                "cached for one chromosome only, so a calibration using it cannot be applied to the "
                "sweep; this row says what leaving it out costs"
            ),
            "leave_chromosome_out_ece": {
                name: calibration(loco[name], labels)["reliability"]["ece"]
                for name in ("sweep", "activity")
                if name in loco
            },
            "heldout_k562": calibration(
                predict(weights["activity"], held_k562, ACTIVITY_FEATURES), held_labels
            ),
        },
    }


# ------------------------------------------------------------------------------------------
# The genome-wide sweep, read through the calibration
# ------------------------------------------------------------------------------------------


def band_of(p: float) -> str:
    for lo, hi in CONFIDENCE_BANDS:
        if lo <= p < hi:
            return f"{lo:g}-{hi if hi <= 1 else 1:g}"
    return "1"


def element_row(
    el: dict[str, Any], key: str, starts: dict[str, int], classes: dict[str, tuple[str, bool]]
) -> dict[str, Any] | None:
    """The features of one element's predicted target, or None with a reason counted by the caller."""
    t = el.get(key)
    if not t or not t.get("gene"):
        return None
    gene = t["gene"]
    tss = starts.get(gene)
    d = tss_distance((el["start"] + el["end"]) // 2, tss)
    if d is None:
        return {"missing": "no_gencode_tss_for_the_predicted_gene"}
    cls = classes.get(el.get("id", ""), ("unknown", False))[0]
    _, drop = crispri.deletion_for([el], gene, CELL)
    return {
        "features": {
            "deletion_drop": drop,
            "top_target": 1.0,
            "log_tss_distance": math.log(d),
            "node_target": 1.0 if (el.get("inferred") or {}).get("gene") == gene else 0.0,
            **class_features(cls),
        },
        "cls": cls,
        "drop": drop,
    }


def sweep_chromosome(
    chrom: str,
    weights: dict[str, list[float]],
    elements: Path = ELEMENTS,
    reference: Path = REFERENCE,
    results: Path = RESULTS,
) -> dict[str, Any]:
    """Every deleted element of one chromosome, banded by its calibrated probability."""
    p = elements / f"{chrom}.json"
    if not p.exists():
        return {"refused": "no deleted elements cached for this chromosome"}
    els = json.loads(p.read_text())
    starts = gene_starts(chrom, reference)
    classes = registry_classes(chrom, results) if (results / f"ccres_{chrom}.bed.gz").exists() else {}
    out: dict[str, Any] = {"elements": len(els)}
    for key, label in (("predicted", "any gene"), ("predicted_coding", "coding gene")):
        missing: Counter[str] = Counter()
        bands: Counter[str] = Counter()
        target_bands: Counter[str] = Counter()
        raw_bands: Counter[str] = Counter()
        by_class: dict[str, Counter[str]] = defaultdict(Counter)
        total, ps = 0, 0.0
        for el in els:
            row = element_row(el, key, starts, classes)
            if row is None:
                missing["no_predicted_target"] += 1
                continue
            if "missing" in row:
                missing[row["missing"]] += 1
                continue
            x = [[row["features"][n] for n in SWEEP_FEATURES]]
            xt = [[row["features"][n] for n in TARGET_FEATURES]]
            prob = sigmoid(crispri.logistic_score(weights["sweep"], x)[0])
            prob_t = sigmoid(crispri.logistic_score(weights["target"], xt)[0])
            bands[band_of(prob)] += 1
            target_bands[band_of(prob_t)] += 1
            raw_bands[drop_band_label(row["drop"])] += 1
            by_class[row["cls"]][band_of(prob_t)] += 1
            total += 1
            ps += prob_t
        out[label] = {
            "targets": total,
            "mean_predicted_target_probability": round(ps / total, 4) if total else None,
            "tested_pair_bands": dict(sorted(bands.items())),
            "predicted_target_bands": dict(sorted(target_bands.items())),
            "drop_bands": {k: raw_bands[k] for k in DROP_LABELS if k in raw_bands},
            "by_registry_class": {c: dict(sorted(v.items())) for c, v in sorted(by_class.items())},
            "missing": dict(sorted(missing.items())),
        }
    return out


def sweep(
    weights: dict[str, list[float]],
    chroms: tuple[str, ...] = CHROMS,
    elements: Path = ELEMENTS,
    reference: Path = REFERENCE,
    results: Path = RESULTS,
    progress: Any = None,
) -> dict[str, Any]:
    """The calibration applied to every chromosome of the sweep, per chromosome and per class."""
    per: dict[str, Any] = {}
    for c in chroms:
        per[c] = sweep_chromosome(c, weights, elements, reference, results)
        if progress:
            progress(f"{c}: {per[c].get('coding gene', {}).get('targets', 0)} coding targets banded")
    totals: dict[str, Any] = {}
    for label in ("any gene", "coding gene"):
        bands: Counter[str] = Counter()
        tbands: Counter[str] = Counter()
        rbands: Counter[str] = Counter()
        by_class: dict[str, Counter[str]] = defaultdict(Counter)
        missing: Counter[str] = Counter()
        total = 0
        for v in per.values():
            d = v.get(label)
            if not d:
                continue
            total += d["targets"]
            bands.update(d["tested_pair_bands"])
            tbands.update(d["predicted_target_bands"])
            rbands.update(d["drop_bands"])
            missing.update(d["missing"])
            for cls, bb in d["by_registry_class"].items():
                by_class[cls].update(bb)
        totals[label] = {
            "targets": total,
            "tested_pair_bands": dict(sorted(bands.items())),
            "predicted_target_bands": dict(sorted(tbands.items())),
            "drop_bands": {k: rbands[k] for k in DROP_LABELS if k in rbands},
            "by_registry_class": {c: dict(sorted(v.items())) for c, v in sorted(by_class.items())},
            "missing": dict(sorted(missing.items())),
        }
    return {"genome_wide": totals, "per_chromosome": per}
