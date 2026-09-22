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

MIN_POOLED_FOR_A_BAND = 30  # pairs; below this no stratum is given a numeric band, however tight
NOT_CALIBRATED = "not calibrated here"  # the band for a population whose interval names no one band

PREREGISTERED_BANDS = (
    "Registered 2026-09-22, after the gate was counted and removed and before any swept target was "
    "re-banded. The band table this judges is `genome_wide.predicted_target_bands` in the "
    "`target_calibration` result, written by `sweep()`.\n"
    "\n"
    "WHY THE OBVIOUS SPLIT IS THE WRONG ONE, AND HOW THAT WAS ESTABLISHED. The lesson of 2026-09-22 "
    "is that the 245-pair curve quotes 0.4126 against a measured 0.0603 one step off its gate, so the "
    "expected repair was to split the swept targets by whether the target is the element's top "
    "predicted target and band the two sides from two curves. That split does not exist in the sweep. "
    "`crispri.deletion_values` sets `top_target` to 1.0 when the pair's gene matches EITHER "
    "`predicted` OR `predicted_coding` at an overlapping element, and `sweep_chromosome` bands "
    "exactly those two keys and nothing else. Every one of the 612,323 swept targets is therefore "
    "inside the gate, and the x6.85 is not an error in the published table: it is the error the "
    "module would make on the population `targets.ElementResponses` has just made askable and which "
    "the sweep does not yet band. Counted before anything was re-fitted: of the 440,377 coding-arm "
    "targets, 331,209 are the element's window head and 109,168 are the coding head only -- but both "
    "are compact-table entries, so both carry `top_target` = 1 and both were inside the fit's "
    "population. The split is real and it is not the one that separates the sweep from its curve.\n"
    "\n"
    "WHAT THE POPULATION MISMATCH ACTUALLY IS. The fitted 245 and the banded 612,323 are both inside "
    "the gate and are nevertheless different populations, on two measured axes. (a) MAGNITUDE: 34.3% "
    "of the fitted pairs carry a predicted K562 drop above 0.2, against 5.8% of the swept targets, "
    "and 15.5% of the fitted pairs sit at a drop of exactly zero against 46.5% of the swept ones. "
    "(b) CELL: 78 of the 245 fitted pairs (31.8%) sit on an element whose target the sweep named on "
    "K562, the one cell line this curve is measured in, against 32,597 of 612,323 swept targets "
    "(5.3%) -- a factor of 6.0. The remaining 94.7% were named on placenta, CD14-positive monocyte, "
    "HepG2, testis, psoas muscle and 300-odd other tracks, and enter a K562 curve through a K562 "
    "drop that is zero for most of them.\n"
    "\n"
    "THE POPULATIONS, AND HOW A TARGET IS ASSIGNED TO ONE. The axis is the predicted K562 deletion "
    "drop, in the five `DROP_BANDS` strata the module already reports `measured_by_drop_band` on. It "
    "is chosen over the top-target axis because it is the axis the two populations actually differ "
    "on, it is computable for every swept target and for every benchmark pair with no extra data, "
    "and it is the axis the fitted rate moves along (0.579, 0.566, 0.875, 1.000, 1.000 on the 245). "
    "Each swept target is assigned the stratum of its own predicted K562 drop -- the same number "
    "`element_row` already computes -- and each benchmark pair the stratum of its own. The counts "
    "registered in advance, 'any gene' arm: drop = 0, 275,821 swept targets against 38 training and "
    "5 held-out pairs; 0 < drop <= 0.1, 239,487 against 83 and 9; 0.1 < drop <= 0.2, 44,299 against "
    "40 and 5; 0.2 < drop <= 0.5, 24,281 against 35 and 11; 0.5 < drop, 9,877 against 49 and 10. So "
    "86.8% of the genome's bands rest on 121 training and 14 held-out pairs. The 'coding gene' arm "
    "is 205,541 / 173,234 / 37,488 / 18,338 / 5,776 on the same strata.\n"
    "\n"
    "THE RULE. A stratum is given a numeric band only if (a) it holds at least "
    "MIN_POOLED_FOR_A_BAND = 30 pooled training and held-out pairs, and (b) the 95% Wilson interval "
    "of its pooled observed rate lies inside a single one of the eight `CONFIDENCE_BANDS`. Otherwise "
    "the stratum's targets are banded NOT_CALIBRATED = 'not calibrated here' and carry the stratum's "
    "pair count, observed rate and interval in place of a number. Pooling training with held-out "
    "pairs for the rate is deliberate and is declared here: the held-out 40 alone give intervals of "
    "width 0.65 and would mark every stratum uncalibrated, which is true but says nothing, so the "
    "interval is read on all 285 on-gate pairs and the quantity being interval-bounded is a base "
    "rate rather than a fitted curve's error. WHAT IS DONE IF A STRATUM IS TOO THIN TO FIT AT ALL: "
    "it is banded 'not calibrated here' with its count, and no curve is fitted for it. That is the "
    "expected outcome for the drop = 0 stratum, whose 43 pooled pairs give 0.5814 [0.4335, 0.7160], "
    "an interval spanning two published bands.\n"
    "\n"
    "THE DIRECTION EXPECTED, AND BY HOW MUCH. Most of the genome loses its number. The registered "
    "prediction is that the two bottom strata (515,308 of 593,765 'any gene' targets, 86.8%) and the "
    "0.1-0.2 stratum (44,299, 7.5%) all fail clause (b) and become 'not calibrated here', and that "
    "only the two top strata (34,158 targets, 5.8%) keep a numeric band, which will be 0.9-1 for "
    "both. That is between 80% and 95% of the table losing its band, and the direction for every "
    "target that keeps one is unchanged or up, never down -- because within the gate the curve was "
    "measured honest in every stratum (quoted 0.4932 against an observed 0.5789 at drop = 0, 0.6139 "
    "against 0.5663, 0.8746 against 0.8750, 0.9808 against 1.0000, 0.9999 against 1.0000). The "
    "correction this lane makes is therefore NOT that the published levels are wrong. It is that "
    "the published table states eight-way bands to two decimal places for 593,765 targets on the "
    "evidence of 285 pairs, 86.8% of them priced by 135, and a band is a claim about a rate whose "
    "interval has to fit inside it.\n"
    "\n"
    "THE FALSIFIER. The split is the wrong one if the two bottom strata's pooled intervals do fit "
    "inside a single published band, because then the drop axis separates nothing the curve has not "
    "already absorbed and the published table stands as written. It is also the wrong one if the "
    "alternative axis measured here -- the target named on K562 against named on another track -- "
    "yields single-band intervals where the drop axis does not; both are computed and both are "
    "reported, and if the cell axis is the cleaner one it replaces the drop axis and this "
    "registration is recorded as wrong on the axis while right on the direction. A third way to "
    "fail: if re-banding moves fewer than half the targets, the registered magnitude was wrong and "
    "is reported as wrong rather than rounded towards.\n"
    "\n"
    "WHAT IS NOT TOUCHED. The published 2026-09-17 band table is annotated and kept, never "
    "overwritten. `sweep()` keeps its output and its keys, the re-banding is a second block beside "
    "it under a new key, and both can be read side by side from the same result. The eight "
    "`CONFIDENCE_BANDS` and the fitted weights do not change, so a reader who wants the old number "
    "still has it; what is added is the population each number was fitted on and whether that "
    "population can carry it. 0 AlphaGenome requests: every number here is on disk."
)


PREREGISTERED_UNION = (
    "Registered 2026-09-22, after the two axes were crossed and BEFORE the union of them was fitted "
    "to anything or used to band a single target. Every count below was computed with no label "
    "read: `crossing_counts` and `crossing_census` take pairs and elements and return sizes, and "
    "neither touches `regulated`.\n"
    "\n"
    "HOW THIS AXIS WAS CHOSEN, WHICH IS THE WHOLE PROBLEM. On 2026-09-22 the re-banding lane "
    "measured both axes it had registered -- the predicted K562 deletion drop and the track the "
    "sweep named the element's target on -- crossed them on the 285 on-gate pairs, and found that "
    "three of the four cells are at 100%: 70 of 70 with both signals, 35 of 35 with the drop alone, "
    "23 of 23 with the track alone, against 96 of 157 with neither. Their union is 128 of 128, "
    "1.0000 [0.9709, 1.0]. It refused to adopt that, and this registration exists because of the "
    "refusal. The axis was chosen AFTER its result was seen, on the only 285 pairs that can be read, "
    "so the interval it reports is the interval of the subset that won a search, and the search is "
    "not in the interval. Nothing below is allowed to forget this.\n"
    "\n"
    "THE AXIS, AND HOW A TARGET IS ASSIGNED TO A STRATUM. Two binary signals, both already computed "
    "for every swept target and every benchmark pair with no extra data. (a) MAGNITUDE: the "
    "predicted deletion drop exceeds UNION_DROP = 0.2, which is the boundary between the two bottom "
    "and the two top `DROP_BANDS` strata and is the split the 2026-09-22 registration already named. "
    "(b) TRACK: the element's predicted target was named on the K562 track, the one cell line this "
    "curve speaks for. A target is in the `union` stratum if either signal is present and in "
    "`neither` if both are absent; `crossing` keeps all four cells (`both`, `drop only`, "
    "`track only`, `neither`) so the union's own homogeneity can be read. For a benchmark pair the "
    "track is read at the element (`pair_track`): on the gate the pair's gene IS the element's "
    "predicted target so this is exactly the axis of 2026-09-22, and off the gate it is the only "
    "reading that exists and is also precisely what a swept target carries, so the axis means one "
    "thing on both populations instead of collapsing to a constant off the gate. The drop for a "
    "benchmark pair is read in that pair's own screen cell line; off the gate it is the sweep's own "
    "windowed value, because the compact table has nothing to say about a gene it did not name.\n"
    "\n"
    "THE POPULATIONS AND THEIR SIZES, COUNTED BEFORE ANY RATE. On the gate, K562, the 285 pairs "
    "that suggested the axis (245 training + 40 held out): both 70, drop only 35, track only 23, "
    "neither 157; union 128. Their screens: Gasperini2019 207, Nasser2021 28, Xie 19, "
    "Schraivogel2020 10, Morris 8, Klann 7, K562_DC_TAP 5, Reilly 1 -- so 72.6% of the evidence for "
    "this axis is one screen. Off the gate, read through the sweep's own window, 10,226 pairs: both "
    "19, drop only 14, track only 1,419, neither 8,774; union 1,452. (That 10,226 pools both "
    "benchmark arms. The 2026-09-22 driver called `add_features` on the training arm only in its "
    "windowed pass, so `scored` dropped every held-out pair and its windowed off-gate table was "
    "training-only; corrected here.) Held out in the five cell lines that are not K562, on the gate "
    "and scored: 8 pairs -- drop only 2, track only 1, neither 5 -- out of 2,063 non-K562 held-out "
    "pairs (WTC11 1,602, HCT116 335, Jurkat 64, GM12878 62), which fail the gate rather than the "
    "cell. Genome-wide, `any gene`, of 593,765 banded targets: both 8,370, drop only 25,788, track "
    "only 23,468, neither 536,139 -- union 57,626 (9.70%), neither 536,139 (90.30%). `coding gene`, "
    "of 440,377: both 6,593, drop only 17,521, track only 20,661, neither 395,602; union 44,775 "
    "(10.17%). These reconcile with 2026-09-22's 32,597 K562-named and 35,217 above a drop of 0.2, "
    "which counted all 612,323 swept targets; over the 593,765 that carry a GENCODE TSS and are "
    "actually banded the same counts are 31,838 and 34,158, and 34,158 is to the digit the count "
    "that registration predicted for its two top strata. On the elements GTEx has distilled eQTLs "
    "for: 6,619 elements carry a hit, 5,100 of them are in the all-enhancer table, and they hold "
    "3,331 banded `any gene` targets -- both 43, drop only 100, track only 118, neither 3,070; "
    "union 261 -- and 2,758 `coding gene` targets, union 228.\n"
    "\n"
    "THE RULE, UNCHANGED FROM 2026-09-22. A stratum is given a numeric band only if it holds at "
    "least MIN_POOLED_FOR_A_BAND = 30 pooled pairs and the 95% Wilson interval of its pooled "
    "observed rate lies inside a single one of the eight `CONFIDENCE_BANDS`. Otherwise its targets "
    "are banded NOT_CALIBRATED = 'not calibrated here' and carry the stratum's count, rate and "
    "interval in place of a number. Training and held-out pairs are pooled for the rate, declared "
    "here for the same reason as before: the quantity bounded is a base rate, not a fitted curve's "
    "error, and the 40 held-out pairs alone bound nothing.\n"
    "\n"
    "THE CLAUSE THAT MATTERS: WHAT IS A TEST HERE AND WHAT IS NOT. The 128 of 128 runs three "
    "questions together, and they are separated here because only one of them has evidence that is "
    "not the hypothesis restated.\n"
    "  (1) THE LEVEL -- is the union stratum's rate inside 0.9-1. The only pairs on which this can "
    "be read are the 285 that chose the axis. A hypothesis and its confirmation cannot be the same "
    "128 pairs, so THIS IS REGISTERED AS A DESCRIPTION AND NOT AS A TEST. No prediction is made for "
    "it, no falsifier is offered for it, and its interval is reported with the statement that it is "
    "the interval of a subset selected for having exactly that interval. Whatever it reads, it "
    "prices no genome-wide target.\n"
    "  (2) THE SEARCH -- could a crossing this clean come out of two axes that carry nothing. This "
    "is priced by permuting the labels over the 285 with the four cells held at their sizes, and "
    "taking as the statistic the best Wilson lower bound any grouping of the four cells reaches "
    "under the adoption rule, so the whole family the search really ran over is inside the null. It "
    "is a genuine test of 'neither axis carries information'. It is ALSO nearly worthless, and that "
    "is registered in advance: the drop axis was registered and measured to separate on 2026-09-22, "
    "so that null is already known to be false, and a small p here is not evidence for the union "
    "over the drop axis alone. It bounds the search; it does not choose the axis.\n"
    "  (3) THE SEPARATION -- does either arm carry information on a population that took no part in "
    "choosing it. This is the only real test available, and it is run on populations that were not "
    "looked at when the union was picked.\n"
    "\n"
    "WHAT INDEPENDENT EVIDENCE WOULD LOOK LIKE, AND WHAT IS REACHABLE AT 0 REQUESTS. It would be a "
    "set of measured element-gene pairs, in the union stratum, that took no part in this choice, "
    "and enough of them to bound a rate inside one published band. Reachable at 0 requests, in "
    "descending order of what it can settle: (i) THE OFF-GATE PAIRS, 1,452 union against 8,774 "
    "neither, read through the sweep's own window; the drop arm there is already known to separate "
    "(2026-09-22 measured 0.0091 at drop = 0 against 84.6% and 83.3% above 0.1) so it is not news, "
    "but the TRACK arm's 1,419 off-gate pairs were never looked at and are the sharpest thing this "
    "lane can ask. (ii) THE GTEx eQTL ARM, an independent assay: 261 union against 3,070 neither "
    "banded targets on elements with distilled cis-eQTLs, asking how often the element's predicted "
    "target is among its eGenes. Three caveats registered in advance -- the 6,619 elements were "
    "distilled against a pre-selected index of 9,286, not a random draw from the sweep; the label "
    "is association in bulk GTEx tissue, not a CRISPRi call, so its base rate is not this one's; "
    "and GTEx has no K562, so the track arm is being asked a cross-context question. It prices "
    "direction and nothing else. (iii) THE NON-K562 HELD-OUT PAIRS, registered in advance as "
    "UNUSABLE: 8 scored on-gate pairs cannot bound anything, and saying so is the result for that "
    "arm. NOT reachable at 0 requests, and named here as what would actually settle the level: a "
    "CRISPRi screen outside the EngreitzLab benchmark with on-gate union pairs; the benchmark's "
    "other cell lines rescored so their pairs reach the gate at all; or the prospective form, which "
    "is to publish the union stratum's top targets as a list and have them tested. The level of "
    "this axis is not confirmable from anything on this disk, and if the separation tests come back "
    "empty then the honest result of this lane is that the union axis cannot be adopted yet.\n"
    "\n"
    "REGISTERED PREDICTIONS AND WHAT EACH OUTCOME MEANS.\n"
    "  P1 (level): none registered; it is a description.\n"
    "  P2 (search): the permutation prices the search at p < 0.01. If it does not, the crossing is "
    "inside what the search can produce from noise and the axis is dead on the spot.\n"
    "  P3 (off the gate, the track arm): the 1,419 track-only pairs are predicted NOT to separate "
    "from the 8,774 neither pairs with non-overlapping 95% intervals. The reason for expecting a "
    "null is stated in advance: 'named on K562' is a property of the element's top target, and the "
    "on-gate association is as easily an account of which pairs a K562 screen could test as of a "
    "mechanism. If the intervals DO separate, the track arm carries information where it was never "
    "looked at and the union has one arm that is more than a testability marker; if they do not, "
    "the union's 23 track-only pairs are a selection artefact and the union collapses to the drop "
    "axis, which was already registered.\n"
    "  P4 (the union's own falsifier): the `neither` stratum is predicted to be HETEROGENEOUS on the "
    "drop axis that was registered first -- its sub-strata were measured at 0.5789, 0.5663 and "
    "0.8750 on the training 245, which no single band contains. If that holds, then banding "
    "`neither` with one number aggregates populations with different rates over 536,139 targets, "
    "90.30% of the genome, which is precisely the failure the track axis was refused for, and the "
    "union axis does not fix the genome: it fixes 9.70% of it and repeats the mistake on the rest.\n"
    "  P5 (eQTL): the union stratum is predicted to carry its predicted target among the element's "
    "eGenes more often than `neither` does. A separation supports the direction of the axis on an "
    "independent assay. No outcome here licenses the 0.9-1 level.\n"
    "\n"
    "THE ADOPTION RULE, FIXED NOW. The union axis is adopted and the genome banded from it only if "
    "(a) P2 holds, AND (b) at least one population that took no part in the choice separates with "
    "non-overlapping 95% intervals, AND (c) the `neither` stratum is homogeneous enough for one "
    "number -- every drop sub-stratum of it with at least 30 pairs has its interval inside the band "
    "the pooled `neither` rate names. If (c) fails but (a) and (b) hold, only the union stratum is "
    "banded and `neither` is banded 'not calibrated here', which bands 9.70% of the genome and "
    "refuses 90.30% of it; that is the registered fallback and it is not a failure. If (b) fails, "
    "NOTHING is banded from this axis and the lane reports that the union cannot be confirmed from "
    "anything on this disk. Under every outcome the union stratum's 0.9-1 is written down as a "
    "description of the 128 pairs that chose it.\n"
    "\n"
    "WHAT IS NOT TOUCHED. The published 2026-09-17 band table and the 2026-09-22 re-banding are "
    "annotated and kept, never overwritten. `sweep()`, `reband()` on the `drop` and `track` axes, "
    "the eight `CONFIDENCE_BANDS` and the fitted weights do not change; the union is a third axis "
    "beside them under its own key, so all four readings can be held against each other from disk. "
    "0 AlphaGenome requests: every number here is cached."
)


# ------------------------------------------------------------------------------------------
# Features: what the benchmark and the sweep both have
# ------------------------------------------------------------------------------------------


def class_features(cls: str) -> dict[str, float]:
    """One-hot of the registry class, with dELS as the reference level."""
    return {f"class_{c}": 1.0 if c == cls else 0.0 for c in CLASS_LEVELS}


def matched_element(
    elements: list[dict[str, Any]],
    gene: str,
    responses: targets.ElementResponses | None = None,
    chrom: str = "",
    cell: str = CELL,
) -> dict[str, Any] | None:
    """The one overlapping element a pair's features come from.

    A pair can sit on more than one deleted element. The element whose top predicted target is the
    pair's gene is the one the deletion feature already speaks for; failing that, the one with the
    largest predicted magnitude in K562; failing that, the first.

    That fallback is itself conditioned on the compact table's one gene: the largest magnitude it
    ranks on is the magnitude predicted for *whatever other gene the table named*, not for the pair's
    own gene, and the class and distance features are then taken from an element chosen on a number
    about a different gene. With `responses`, the sweep's own per-element cache is asked about this
    pair's gene at each candidate and the element with the strongest predicted fall for it wins; the
    old rule is kept for the pairs the cache cannot answer, so `responses=None` reproduces exactly.
    """
    if not elements:
        return None
    for e in elements:
        for key in ("predicted_coding", "predicted"):
            p = e.get(key)
            if p and p["gene"] == gene:
                return e
    if responses is not None and chrom:
        asked = [(e, responses.value(chrom, e["id"], gene, cell)) for e in elements]
        answered = [(e, v) for e, v in asked if v is not None]
        if answered:
            return min(answered, key=lambda ev: ev[1])[0]
    return max(elements, key=lambda e: abs((e.get("predicted_by_cell") or {}).get(CELL) or 0.0))


def tss_distance(midpoint: int, tss: int | None) -> float | None:
    return None if tss is None else float(max(MIN_DISTANCE, abs(midpoint - tss)))


def add_features(
    pairs: list[crispri.Pair],
    table: crispri.DeletionTable,
    reference: Path = REFERENCE,
    results: Path = RESULTS,
    responses: targets.ElementResponses | None = None,
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
            el = matched_element(els, p.gene, responses, p.chrom, p.cell)
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


# ------------------------------------------------------------------------------------------
# Re-banding: a band drawn from the population it is quoted for, or none at all
# ------------------------------------------------------------------------------------------

BAND_LABELS = tuple(f"{lo:g}-{hi if hi <= 1 else 1:g}" for lo, hi in CONFIDENCE_BANDS)


def bands_touched(lo: float, hi: float) -> list[str]:
    """Every published band the interval [lo, hi] reaches into."""
    return [f"{a:g}-{b if b <= 1 else 1:g}" for a, b in CONFIDENCE_BANDS if b > lo and a <= hi]


def observed_band(k: int, n: int, min_pooled: int = MIN_POOLED_FOR_A_BAND) -> dict[str, Any]:
    """The band a measured population can carry, or the reason it can carry none.

    A band is a claim about a rate, so the claim is allowed only when the rate's 95% Wilson
    interval fits inside one band and rests on enough pairs to be worth stating. Anything else
    returns `NOT_CALIBRATED` and carries the count and the interval in place of a number, which is
    the honest answer for a stratum the benchmark barely reaches.
    """
    if n == 0:
        return {"band": NOT_CALIBRATED, "why": "no measured pair in this population", "pairs": 0}
    rate = k / n
    lo, hi = wilson(k, n)
    touched = bands_touched(lo, hi)
    out = {
        "pairs": n,
        "regulated": k,
        "observed_rate": round(rate, 4),
        "ci95": [round(lo, 4), round(hi, 4)],
        "bands_the_interval_touches": touched,
    }
    if n < min_pooled:
        return {**out, "band": NOT_CALIBRATED, "why": f"fewer than {min_pooled} pooled pairs"}
    if len(touched) != 1:
        return {
            **out,
            "band": NOT_CALIBRATED,
            "why": f"the interval spans {len(touched)} published bands",
        }
    return {**out, "band": touched[0], "why": "the interval fits inside one published band"}


def population_bands(
    pairs: list[crispri.Pair], key: Any, min_pooled: int = MIN_POOLED_FOR_A_BAND
) -> dict[str, Any]:
    """Split pairs by `key` and give each population the band its own measured rate can carry."""
    acc: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    for p in pairs:
        s = key(p)
        acc[s][0] += 1
        acc[s][1] += int(p.regulated)
    return {s: observed_band(k, n, min_pooled) for s, (n, k) in sorted(acc.items())}


def drop_stratum(p: crispri.Pair) -> str:
    return drop_band_label(p.features["deletion_drop"])


def element_tissue(el: dict[str, Any], key: str) -> str:
    """The alternative axis: the track the sweep named this element's target on."""
    return CELL if ((el.get(key) or {}).get("tissue") or "?") == CELL else "another track"


UNION_DROP = 0.2  # the drop above which a target sits in the union axis's magnitude arm
UNION_CELLS = ("both", "drop only", "track only", "neither")
UNION_STRATA = ("union", "neither")


def union_cell(drop: float, track: str) -> str:
    """Which of the four cells of the two measured axes crossed a target or a pair sits in."""
    a, b = drop > UNION_DROP, track == CELL
    if a and b:
        return "both"
    if a:
        return "drop only"
    return "track only" if b else "neither"


def union_stratum_of(drop: float, track: str) -> str:
    """The two-stratum union axis: either signal present, or neither."""
    return "neither" if union_cell(drop, track) == "neither" else "union"


def pair_track(p: crispri.Pair, table: crispri.DeletionTable) -> str:
    """The track the sweep named this pair's element's target on, read the element's way.

    On the gate the pair's gene *is* that target, so this returns exactly what the 2026-09-22
    alternative axis read. Off the gate the pair's gene is not the element's target and only the
    element-level reading exists — which is also precisely what a swept target carries, so the axis
    means the same thing on both populations instead of collapsing to one value off the gate.
    """
    el = matched_element(table.overlapping(p.chrom, p.start, p.end), p.gene)
    if el is None:
        return "another track"
    for k in ("predicted_coding", "predicted"):
        q = el.get(k)
        if q and q.get("gene") == p.gene:
            return CELL if (q.get("tissue") or "?") == CELL else "another track"
    return element_tissue(el, "predicted")


def pair_union_cell(p: crispri.Pair, table: crispri.DeletionTable) -> str:
    return union_cell(p.features["deletion_drop"], pair_track(p, table))


def pair_union_stratum(p: crispri.Pair, table: crispri.DeletionTable) -> str:
    return union_stratum_of(p.features["deletion_drop"], pair_track(p, table))


def crossing_counts(pairs: list[crispri.Pair], table: crispri.DeletionTable) -> dict[str, int]:
    """How many pairs sit in each cell of the crossing. Counts only: no label is read here."""
    acc: Counter[str] = Counter()
    for p in pairs:
        acc[pair_union_cell(p, table)] += 1
    return {c: acc.get(c, 0) for c in UNION_CELLS}


def crossing_census(
    chroms: tuple[str, ...] = CHROMS,
    elements: Path = ELEMENTS,
    reference: Path = REFERENCE,
    results: Path = RESULTS,
    progress: Any = None,
) -> dict[str, Any]:
    """How many swept targets sit in each cell of the crossing, genome-wide.

    No weight is fitted and no label is read, so this can be counted before the registration is
    written -- which is the point: a registration that states its population sizes in advance has
    to be able to count them without touching the outcome.
    """
    per: dict[str, Any] = {}
    for chrom in chroms:
        p = elements / f"{chrom}.json"
        if not p.exists():
            continue
        els = json.loads(p.read_text())
        starts = gene_starts(chrom, reference)
        classes = registry_classes(chrom, results) if (results / f"ccres_{chrom}.bed.gz").exists() else {}
        row: dict[str, Any] = {}
        for key, label in (("predicted", "any gene"), ("predicted_coding", "coding gene")):
            acc: Counter[str] = Counter()
            for el in els:
                r = element_row(el, key, starts, classes)
                if r is None or "missing" in r:
                    continue
                acc[union_cell(r["drop"], element_tissue(el, key))] += 1
            row[label] = {c: acc.get(c, 0) for c in UNION_CELLS}
        per[chrom] = row
        if progress:
            progress(f"{chrom}: {sum(row['any gene'].values())} targets counted")
    totals: dict[str, Any] = {}
    for label in ("any gene", "coding gene"):
        acc = Counter()
        for v in per.values():
            acc.update(v[label])
        cells = {c: acc.get(c, 0) for c in UNION_CELLS}
        totals[label] = {
            "targets": sum(cells.values()),
            "by_crossing_cell": cells,
            "union": sum(v for k, v in cells.items() if k != "neither"),
            "neither": cells["neither"],
        }
    return {"genome_wide": totals, "per_chromosome": per}


def eqtl_separation(
    hits: dict[str, list[dict[str, Any]]],
    symbol_of: dict[str, str],
    chroms: tuple[str, ...] = CHROMS,
    elements: Path = ELEMENTS,
    reference: Path = REFERENCE,
    progress: Any = None,
) -> dict[str, Any]:
    """Does the union stratum's predicted target turn up among the element's measured eGenes?

    An independent assay, and the only one on this disk that speaks about *which gene* rather than
    whether a sequence is an enhancer: a significant cis-eQTL inside an element ties a gene to it in
    hundreds of donors. The label is association in bulk GTEx tissue rather than a CRISPRi call, so
    this prices the direction of the axis and never the level of a band; the registration says so,
    along with the fact that the distilled elements are a pre-selected index rather than a draw.
    """
    per: dict[str, dict[str, list[int]]] = {
        label: {c: [0, 0] for c in UNION_CELLS} for label in ("any gene", "coding gene")
    }
    for chrom in chroms:
        p = elements / f"{chrom}.json"
        if not p.exists():
            continue
        els = [e for e in json.loads(p.read_text()) if f"{chrom}:{e['start']}-{e['end']}" in hits]
        if not els:
            continue
        starts = gene_starts(chrom, reference)
        for el in els:
            egenes = {
                symbol_of.get(h["gene_id"], h["gene_id"]) for h in hits[f"{chrom}:{el['start']}-{el['end']}"]
            }
            for key, label in (("predicted", "any gene"), ("predicted_coding", "coding gene")):
                row = element_row(el, key, starts, {})
                if row is None or "missing" in row:
                    continue
                cell = union_cell(row["drop"], element_tissue(el, key))
                per[label][cell][0] += 1
                per[label][cell][1] += int((el.get(key) or {}).get("gene") in egenes)
        if progress:
            progress(f"{chrom}: {len(els)} elements with an eQTL")
    out: dict[str, Any] = {}
    for label, cells in per.items():
        n = sum(v[0] for c, v in cells.items() if c != "neither")
        k = sum(v[1] for c, v in cells.items() if c != "neither")
        out[label] = {
            "by crossing cell": {c: observed_band(kk, nn, min_pooled=1) for c, (nn, kk) in cells.items()},
            "union": observed_band(k, n, min_pooled=1),
            "neither": observed_band(cells["neither"][1], cells["neither"][0], min_pooled=1),
        }
        out[label]["separated"] = separated(out[label]["union"], out[label]["neither"])
    return out


def separated(a: dict[str, Any], b: dict[str, Any]) -> bool | None:
    """Do two measured rates have non-overlapping 95% intervals? None when either has no pairs."""
    if not a.get("pairs") or not b.get("pairs"):
        return None
    return bool(a["ci95"][0] > b["ci95"][1] or b["ci95"][0] > a["ci95"][1])


def move_of(old: str, new: str) -> str:
    """Which way a target's band went: `lost`, `up`, `down` or `unchanged`."""
    if new == NOT_CALIBRATED:
        return "lost its band"
    if old == new:
        return "unchanged"
    return "up" if BAND_LABELS.index(new) > BAND_LABELS.index(old) else "down"


def reband_chromosome(
    chrom: str,
    weights: dict[str, list[float]],
    stratum_band: dict[str, Any],
    elements: Path = ELEMENTS,
    reference: Path = REFERENCE,
    results: Path = RESULTS,
    axis: str = "drop",
) -> dict[str, Any]:
    """One chromosome's targets, banded from their own stratum, beside the band they carry today.

    `axis` picks which population a target belongs to: `drop`, the registered axis, is the predicted
    K562 deletion drop in its `DROP_BANDS` stratum; `track` is the alternative the registration named
    as the thing that would replace it, the track the sweep named this element's target on.

    Two more since 2026-09-22: `union` is the two-stratum axis registered as `PREREGISTERED_UNION`,
    either signal present or neither, and `crossing` keeps the same two axes as four cells so the
    union stratum's own homogeneity can be read rather than assumed.
    """
    p = elements / f"{chrom}.json"
    if not p.exists():
        return {"refused": "no deleted elements cached for this chromosome"}
    els = json.loads(p.read_text())
    starts = gene_starts(chrom, reference)
    classes = registry_classes(chrom, results) if (results / f"ccres_{chrom}.bed.gz").exists() else {}
    out: dict[str, Any] = {"elements": len(els)}
    for key, label in (("predicted", "any gene"), ("predicted_coding", "coding gene")):
        published: Counter[str] = Counter()
        rebanded: Counter[str] = Counter()
        moves: Counter[str] = Counter()
        transitions: Counter[str] = Counter()
        by_tissue: Counter[str] = Counter()
        by_stratum: Counter[str] = Counter()
        total = 0
        for el in els:
            row = element_row(el, key, starts, classes)
            if row is None or "missing" in row:
                continue
            xt = [[row["features"][n] for n in TARGET_FEATURES]]
            old = band_of(sigmoid(crispri.logistic_score(weights["target"], xt)[0]))
            tissue = element_tissue(el, key)
            s = tissue if axis == "track" else drop_band_label(row["drop"])
            if axis == "union":  # registered 2026-09-22, beside the two axes above, not over them
                s = union_stratum_of(row["drop"], tissue)
            elif axis == "crossing":
                s = union_cell(row["drop"], tissue)
            new = stratum_band.get(s, {}).get("band", NOT_CALIBRATED)
            published[old] += 1
            rebanded[new] += 1
            moves[move_of(old, new)] += 1
            transitions[f"{old} -> {new}"] += 1
            by_tissue[tissue] += 1
            by_stratum[s] += 1
            total += 1
        out[label] = {
            "targets": total,
            "published_bands": dict(sorted(published.items())),
            "rebanded": dict(sorted(rebanded.items())),
            "moves": dict(sorted(moves.items())),
            "transitions": dict(sorted(transitions.items())),
            "named_on_track": dict(sorted(by_tissue.items())),
            "by_stratum": dict(sorted(by_stratum.items())),
        }
    return out


def reband(
    weights: dict[str, list[float]],
    stratum_band: dict[str, Any],
    chroms: tuple[str, ...] = CHROMS,
    elements: Path = ELEMENTS,
    reference: Path = REFERENCE,
    results: Path = RESULTS,
    progress: Any = None,
    axis: str = "drop",
) -> dict[str, Any]:
    """Every chromosome re-banded, with the published band kept beside the new one."""
    per: dict[str, Any] = {}
    for c in chroms:
        per[c] = reband_chromosome(c, weights, stratum_band, elements, reference, results, axis)
        if progress:
            progress(f"{c}: {per[c].get('any gene', {}).get('targets', 0)} targets re-banded")
    totals: dict[str, Any] = {}
    for label in ("any gene", "coding gene"):
        acc = {
            k: Counter() for k in ("published_bands", "rebanded", "moves", "transitions", "named_on_track")
        }
        acc["by_stratum"] = Counter()  # added 2026-09-22: coverage has to be readable whatever the bands are
        total = 0
        for v in per.values():
            d = v.get(label)
            if not d:
                continue
            total += d["targets"]
            for k, c in acc.items():
                c.update(d[k])
        totals[label] = {
            "targets": total,
            **{k: dict(sorted(c.items())) for k, c in acc.items()},
        }
    return {"genome_wide": totals, "per_chromosome": per}


# ------------------------------------------------------------------------------------------
# The gate removed: the same calibration fitted through the sweep's own window
# ------------------------------------------------------------------------------------------


def answered(pairs: list[crispri.Pair]) -> list[crispri.Pair]:
    """The pairs the sweep actually predicted a change for, on this gene, in this cell.

    This is the population the magnitude feature means anything on. Without the cache it is the
    top-target pairs and nothing else, which is the gate; with it, every pair whose gene the scorer
    kept in its window.
    """
    return [p for p in pairs if p.features.get("deletion_answered")]


def flattening_check(before: dict[str, Any], after: dict[str, Any]) -> dict[str, Any]:
    """Did reliability improve, or did the curve merely flatten onto the base rate?

    Registered in `PREREGISTERED_GATE` before the re-fit: a curve that predicts the base rate
    everywhere sits inside every equal-count bin and has learnt nothing, so a rise in the bin count
    is only an improvement if the spread of the predicted probabilities held and the ranking did not
    get worse. `flattened` is True when the bins improved on a narrower range or a lower AUPRC.
    """

    def spread(cal: dict[str, Any]) -> float:
        rows = cal["reliability"]["rows"]
        return round(max(r["mean_predicted"] for r in rows) - min(r["mean_predicted"] for r in rows), 4)

    bins = (before["reliability"]["bins_consistent"], after["reliability"]["bins_consistent"])
    rng = (spread(before), spread(after))
    auprc = (before.get("auprc"), after.get("auprc"))
    improved = bins[1] > bins[0]
    return {
        "bins_consistent": list(bins),
        "predicted_range_across_the_bins": list(rng),
        "auprc": list(auprc),
        "reliability_improved": improved,
        "range_narrowed": rng[1] < rng[0],
        "auprc_fell": bool(auprc[0] is not None and auprc[1] is not None and auprc[1] < auprc[0]),
        "flattened": bool(improved and (rng[1] < rng[0] or (auprc[1] is not None and auprc[1] < auprc[0]))),
    }


def off_gate_quote(
    weights: list[float], cols: tuple[str, ...], rows: list[crispri.Pair], what: str
) -> dict[str, Any]:
    """What a curve says about pairs its own population excludes, against what those pairs did.

    The registered clause: the predicted-target curve is fitted on 245 pairs at a base rate of 76.7%
    and the sweep quotes it for every target. Asked about a pair the gate would have excluded, does
    it land near that pair's own rate or near the rate it was fitted on?
    """
    if not rows:
        return {"pairs": 0, "what": what}
    p = predict(weights, rows, cols)
    k = sum(r.regulated for r in rows)
    observed = k / len(rows)
    mean_p = sum(p) / len(p)
    lo, hi = wilson(k, len(rows))
    return {
        "what": what,
        "pairs": len(rows),
        "regulated": k,
        "observed_rate": round(observed, 4),
        "observed_ci95": [round(lo, 4), round(hi, 4)],
        "mean_predicted": round(mean_p, 4),
        "quoted_over_observed": round(mean_p / observed, 2) if observed else None,
        "log_odds_gap": round(logit(mean_p) - logit(observed), 4),
        "inside_the_observed_interval": bool(lo <= mean_p <= hi),
    }


def window_coverage(
    responses: targets.ElementResponses,
    chroms: tuple[str, ...],
    elements: Path = ELEMENTS,
    cell: str = CELL,
) -> dict[str, Any]:
    """How many (element, gene) questions the sweep can be asked, against the one it used to answer.

    The compact table keeps one gene per element, so a user could ask about one gene per element and
    got a structural zero for every other. The cache holds the scorer's whole window, so the
    denominator of "what fraction of the pairs a user would ask about can be answered" changes by the
    number of genes in a window. Measured on the named chromosomes only: the tree is 775 MB and one
    archive is held at a time.
    """
    out: dict[str, Any] = {}
    for chrom in chroms:
        p = elements / f"{chrom}.json"
        if not p.exists():
            out[chrom] = {"refused": "no deleted elements cached for this chromosome"}
            continue
        els = json.loads(p.read_text())
        named = sum(1 for e in els if (e.get("predicted") or {}).get("gene"))
        genes, cached, sizes = 0, 0, []
        for e in els:
            window = responses.ranked(chrom, e["id"], cell)
            if window:
                cached += 1
                genes += len(window)
                sizes.append(len(window))
        out[chrom] = {
            "elements": len(els),
            "elements_the_compact_table_names_a_gene_for": named,
            "elements_in_the_response_cache": cached,
            "askable_pairs_before": named,
            "askable_pairs_now": genes,
            "genes_per_window_median": sorted(sizes)[len(sizes) // 2] if sizes else None,
            "factor": round(genes / named, 1) if named else None,
        }
    return out


def compare_the_gate(
    training: list[crispri.Pair],
    heldout: list[crispri.Pair],
    table: crispri.DeletionTable,
    cache: crispri.ElementCache,
    responses: targets.ElementResponses,
    reference: Path = REFERENCE,
    results: Path = RESULTS,
    coverage_chroms: tuple[str, ...] = ("chr21", "chr22"),
) -> dict[str, Any]:
    """The published calibration and the same calibration re-fitted through the window, side by side.

    Two arms on the same pairs. The shipped arm annotates from the compact table, exactly as
    `score()` does, and reproduces the published numbers — that is the control. The re-fitted arm
    passes the sweep's per-element cache, so `deletion_drop` is what the sweep predicted for the
    pair's own gene and `matched_element` chooses on the same. Nothing else moves: the feature
    columns, the bins, the Wilson rule and the held-out population are the published ones, so the
    bin counts are comparable.
    """
    # Which held-out pairs the gate excludes but the sweep did answer, fixed once and by index:
    # both arms annotate the same Pair objects in place, so the membership cannot be carried on them.
    crispri.annotate(heldout, table, cache)
    add_features(heldout, table, reference, results, responses)
    k562 = [p for p in heldout if p.cell == CELL]
    off_idx = [
        i
        for i, p in enumerate(k562)
        if p.covered
        and p.features.get("scored")
        and not p.features["top_target"]
        and p.features.get("deletion_answered")
    ]
    top_idx = [
        i for i, p in enumerate(k562) if p.covered and p.features.get("scored") and p.features["top_target"]
    ]

    arms: dict[str, Any] = {}
    kept: dict[str, Any] = {}
    for arm, c, r in (
        ("shipped, the compact table", None, None),
        ("re-fitted, the window", cache, responses),
    ):
        crispri.annotate(training, table, c)
        crispri.annotate(heldout, table, c)
        add_features(training, table, reference, results, r)
        add_features(heldout, table, reference, results, r)
        train = scored(training)
        held = scored([p for p in heldout if p.cell == CELL])
        held_labels = [p.regulated for p in held]
        train_top = [p for p in train if p.features["top_target"]]
        w = {
            "sweep": fit(train, SWEEP_FEATURES),
            "drop": fit(train, ("deletion_drop",)),
            "target": fit(train_top, TARGET_FEATURES),
        }
        held_cov = stratum_coverage([p for p in heldout if p.cell == CELL])
        held_strata = [stratum(p) for p in held]
        cal = {
            name: calibration(s, held_labels, BINS, held_strata, held_cov)
            for name, s in scorers(held, w).items()
        }
        passed, checks = judge(cal)
        train_ans, held_ans = answered(train), answered(held)
        off = [k562[i] for i in off_idx]
        on = [k562[i] for i in top_idx]
        arms[arm] = {
            "populations": {
                "training, every feature present": gate_population(train),
                "held-out K562, every feature present": gate_population(held),
                "training, the sweep answered this gene": gate_population(train_ans),
                "held-out K562, the sweep answered this gene": gate_population(held_ans),
            },
            "heldout_k562": cal,
            "verdict": "passed" if passed else "failed",
            "checks": checks,
            "prevalence_shift": {
                name: prevalence_shift(s, held_labels) for name, s in scorers(held, w).items()
            },
            "weights": dict(
                zip(("intercept", *SWEEP_FEATURES), (round(v, 4) for v in w["sweep"]), strict=True)
            ),
            "predicted_target_weights": dict(
                zip(("intercept", *TARGET_FEATURES), (round(v, 4) for v in w["target"]), strict=True)
            ),
            "measured_by_drop_band, the answered held-out pairs": by_drop_band(held_ans),
            "measured_by_drop_band, the answered pairs the gate excludes": by_drop_band(off),
            "measured_by_drop_band, the answered training pairs the gate excludes": by_drop_band(
                [p for p in train_ans if not p.features["top_target"]]
            ),
            "measured_by_drop_band, both, after the held-out test": by_drop_band(
                [p for p in train_ans if not p.features["top_target"]] + off
            ),
            "what_a_confidence_means_off_the_gate": {
                "the predicted-target curve, on the pairs the gate excludes": off_gate_quote(
                    w["target"],
                    TARGET_FEATURES,
                    off,
                    "fitted on the top-target pairs alone; the curve the sweep's bands come from",
                ),
                "the predicted-target curve, on the pairs it admits": off_gate_quote(
                    w["target"], TARGET_FEATURES, on, "its own population, for comparison"
                ),
                "the tested-pair curve, on the pairs the gate excludes": off_gate_quote(
                    w["sweep"], SWEEP_FEATURES, off, "fitted on every scored training pair"
                ),
            },
        }
        kept[arm] = {"w": w, "held": held}

    ship, refit = kept["shipped, the compact table"], kept["re-fitted, the window"]
    return {
        "preregistered": PREREGISTERED_GATE,
        "arms": arms,
        "off_gate_note": (
            "the registered clause: what a confidence means when it is quoted for a pair the gate "
            "would have excluded. Each arm reports it with its own features, so the shipped arm says "
            "what a user is quoted today and the re-fitted arm says what the same curve says once the "
            "sweep's own answer for that gene reaches it"
        ),
        "the_same_heldout_population": {
            "pairs": len(ship["held"]),
            "regulated": sum(p.regulated for p in ship["held"]),
            "identical_in_both_arms": len(ship["held"]) == len(refit["held"]),
            "pairs_the_gate_excludes_and_the_sweep_answers": len(off_idx),
            "pairs_the_gate_admits": len(top_idx),
        },
        "flattening_check": flattening_check(
            arms["shipped, the compact table"]["heldout_k562"]["calibrated, sweep features"],
            arms["re-fitted, the window"]["heldout_k562"]["calibrated, sweep features"],
        ),
        "coverage": {
            "note": (
                "what fraction of the pairs a user would ask about can be answered at all: on the "
                "benchmark, the pairs whose own gene the sweep predicted a change for; on the sweep, "
                "the (element, gene) questions the window makes askable against the one gene the "
                "compact table kept. The sweep figure is measured on two chromosomes, not genome-wide"
            ),
            "benchmark, training K562 on a deleted element": {
                "pairs": sum(1 for p in training if p.covered),
                "answerable_before": sum(1 for p in training if p.covered and p.features["top_target"]),
                "answerable_now": sum(
                    1 for p in training if p.covered and p.features.get("deletion_answered")
                ),
            },
            "benchmark, held-out K562 on a deleted element": {
                "pairs": sum(1 for p in heldout if p.cell == CELL and p.covered),
                "answerable_before": sum(
                    1 for p in heldout if p.cell == CELL and p.covered and p.features["top_target"]
                ),
                "answerable_now": sum(
                    1 for p in heldout if p.cell == CELL and p.covered and p.features.get("deletion_answered")
                ),
            },
            "sweep": window_coverage(responses, coverage_chroms),
        },
    }


# ------------------------------------------------------------------------------------------
# The prevalence term: a per-screen intercept, registered before it is fitted
# ------------------------------------------------------------------------------------------

SCREEN_MIN_PAIRS = 100  # a base rate read on fewer pairs than this is not a base rate
SCREEN_MIN_POSITIVES = 10
SMALL_SCREEN = "small screens, pooled"
SCREENS_IMPROVING_REQUIRED = 4  # of the 6 screens that clear the size rule
ECE_FALL_REQUIRED = 0.005  # half the gap the single pooled shift closed (0.0231 -> 0.0100)
RESIDUAL_SHIFT_BAR = 0.34  # half the +0.679 / +0.7025 the single pooled correction needed
PREREGISTERED_PREVALENCE = (
    "Registered 2026-09-22 (sixth), after the screen census below was counted and before any offset "
    "was fitted or any reliability table was read.\n"
    "\n"
    "THE CENSUS, COUNTED FIRST. On the scored population (covered, every feature present) the "
    "calibration is fitted on three screens -- Gasperini2019 4,776 pairs / 335 regulated / 7.01%, "
    "Nasser2021 2,931 / 80 / 2.73%, Schraivogel2020 1,089 / 22 / 2.02% -- pooling to 8,796 at 4.97%. "
    "It is read on five -- K562_DC_TAP 1,084 / 11 / 1.01%, Xie 416 / 41 / 9.86%, Morris 177 / 34 / "
    "19.21%, Klann 32 / 20 / 62.50%, Reilly 6 / 6 / 100% -- pooling to 1,715 at 6.53%. The spread "
    "within each table is larger than the gap between them: 3.5x across the fitted screens and 99x "
    "across the read ones, against the 1.31x the 2026-09-17 section named. NOT ONE SCREEN NAME IS IN "
    "BOTH TABLES: the intersection of the two screen sets is empty. On the gate, which is the "
    "population that bands the genome, there are 285 pairs and Gasperini2019 is 207 of them (72.6%); "
    "no other screen reaches 30, so no rule of any kind can fit a per-screen intercept there.\n"
    "\n"
    "THE MODEL. logit p_i = logit(pi_s(i)) + a + b.x_i, where pi_s is screen s's own measured base "
    "rate on the scored population, entering as an OFFSET -- a known constant with its coefficient "
    "fixed at 1 and no parameter spent -- while a single intercept a and a single shape b are shared "
    "by every screen. The comparator is the shipped model, logit p_i = a + b.x_i, fitted and read on "
    "the same pairs.\n"
    "\n"
    "WHY AN OFFSET AND NOT A RE-FIT OR A RECALIBRATION. A free per-screen intercept spends one "
    "parameter per screen, and Reilly's 6 pairs at 6 of 6 separate, so that parameter is infinite; "
    "the offset spends none and a screen of any size can carry one. A recalibration on the held-out "
    "labels (Platt, isotonic) moves slope and level together, so it cannot be checked out of sample "
    "at all and it would dissolve the pre-registered claim rather than correct it. The offset also "
    "states the hypothesis the 2026-09-17 failure actually supports and no more: the shape transfers "
    "and the level does not -- all four failing bins failed in the same direction, none in the other "
    "-- so exactly one number per screen should move, and that number is a property the screen "
    "reports about itself, its own base rate, with no model in it.\n"
    "\n"
    "WHICH SCREENS GET THEIR OWN TERM, DECIDED ON SIZE ALONE. A screen carries its own offset when "
    "its scored population has at least 100 pairs AND at least 10 regulated pairs; below either "
    "threshold it joins one `small screens, pooled` stratum which takes the pooled base rate of the "
    "below-threshold screens in its own table. On the census that admits all three fitted screens and "
    "three of the five read ones (K562_DC_TAP, Xie, Morris), and pools Klann and Reilly into 38 pairs "
    "at 68.42%. Six screens carry a term. This rule was fixed on the census counts above, before any "
    "reliability table was read.\n"
    "\n"
    "WHAT IS EXPECTED, AND WHAT CAN FAIL. (E1) Read on the held-out K562 pairs with each screen's own "
    "held-out base rate, bins inside their Wilson intervals go 6 of 10 to 9 or 10 of 10. This is "
    "registered as a DESCRIPTION THAT CANNOT FAIL and it upgrades nothing: five numbers fitted to the "
    "same labels must beat the one number that already gave 9 of 10. It is reported because the "
    "per-screen shifts are the quantity of interest, not because it tests anything. (E2) The test "
    "that can fail is leave-one-screen-out over the six screens that carry a term: b is fitted on the "
    "other five with their offsets, and the held-back screen is scored with its own base rate as its "
    "offset and nothing else of its own. Each screen is read in 10 equal-count bins if it has at "
    "least 300 scored pairs and in 3 below that, the same count for both models, so the comparison is "
    "like for like; a tie does not count as an improvement. Registered: the offset model puts more "
    "bins inside their Wilson intervals than the pooled model on at least 4 of the 6 screens, and the "
    "pairs-weighted expected calibration error over the six falls by at least 0.005. (E3) The "
    "residual per-screen "
    "log-odds shift still needed AFTER the offset is applied has a median absolute value below 0.34, "
    "half of what the single pooled correction needed; if the base rate is the whole story it is near "
    "zero.\n"
    "\n"
    "THE FALSIFIER. The prevalence term is the wrong correction if fewer than 4 of the 6 screens "
    "improve their bin count under leave-one-screen-out, or if the pairs-weighted expected "
    "calibration error does not fall by 0.005, or if the median residual shift is not below 0.34 -- "
    "any of the three, and the level gap is not a per-screen prevalence and something else moved. A "
    "fourth kills it even if the three pass: if the offset model's AUPRC under leave-one-screen-out "
    "falls below the pooled model's, the bins were bought by flattening the curve onto the base rate, "
    "which `flattening_check` refuses everywhere else in this module and which is refused here.\n"
    "\n"
    "WHAT A BAND MEANS FOR A TARGET IN NO SCREEN AT ALL, WHICH IS ALMOST THE WHOLE GENOME. Registered "
    "before the fit, because the honest answer is negative and no result of this fit can change it: "
    "IT CANNOT BE QUOTED. (a) The offset's input is a screen's own base rate. A sweep target is in no "
    "screen, so that number does not exist for it; the only substitute is a guess at what some future "
    "screen's base rate would be, and across the six K562 screens of this one benchmark that number "
    "runs from 1.01% to 100%, a hundredfold. A correction whose input spans a hundredfold is not a "
    "correction. (b) The two screen sets do not intersect, so the term is not transportable even "
    "between TESTED pairs: Xie's intercept cannot be fitted on Gasperini, Nasser and Schraivogel, it "
    "can only be read off Xie's own labels. A per-screen intercept is a thing that exists after a "
    "screen has been run and never before it. (c) The population that bands the genome is the on-gate "
    "one, 285 pairs, 72.6% of them one screen, where no size rule can fit per-screen terms at all. "
    "WHAT SHOULD BE QUOTED INSTEAD, registered now: the two quantities that need no intercept. The "
    "SHAPE -- what leave-one-screen-out tests -- reported as an ordering or a likelihood ratio and "
    "never as a probability; and the MEASURED BAND ON A STRATUM, which is the 2026-09-22 re-banding's "
    "answer, a rate read directly on a stratum of pairs with its population and its interval beside "
    "it and a refusal where the stratum is thin. The prevalence term does not unlock the 593,765 "
    "targets the 2026-09-17 table priced. It says why that table was never quotable, and it leaves "
    "the re-banding's 5.8% and the union axis's 9.70% as the honest genome-wide coverage"
)


def screen_key(p: crispri.Pair) -> str:
    """A screen is a dataset in a cell line: Nasser2021 in K562 and in GM12878 chose pairs differently."""
    return f"{p.dataset}|{p.cell}"


def screen_census(pairs: list[crispri.Pair]) -> dict[str, Any]:
    """Every screen contributing to a scored population, with its size, its rate and its on-gate share.

    Counted before anything is fitted. `own_term` applies the registered size rule and nothing else.
    """
    n: Counter[str] = Counter()
    k: Counter[str] = Counter()
    gn: Counter[str] = Counter()
    gk: Counter[str] = Counter()
    for p in pairs:
        s = screen_key(p)
        n[s] += 1
        k[s] += int(p.regulated)
        if p.features.get("top_target"):
            gn[s] += 1
            gk[s] += int(p.regulated)
    return {
        s: {
            "pairs": n[s],
            "regulated": k[s],
            "rate": round(k[s] / n[s], 4),
            "on_gate_pairs": gn[s],
            "on_gate_regulated": gk[s],
            "on_gate_rate": round(gk[s] / gn[s], 4) if gn[s] else None,
            "own_term": bool(n[s] >= SCREEN_MIN_PAIRS and k[s] >= SCREEN_MIN_POSITIVES),
        }
        for s in sorted(n, key=lambda s: (-n[s], s))
    }


def screen_strata(pairs: list[crispri.Pair]) -> dict[str, str]:
    """Screen name to the stratum that carries its offset: itself, or the pooled small one."""
    census = screen_census(pairs)
    return {s: (s if c["own_term"] else SMALL_SCREEN) for s, c in census.items()}


def screen_base_rates(pairs: list[crispri.Pair], strata: dict[str, str]) -> dict[str, float]:
    """Each stratum's own measured base rate: the one number the offset carries, with no model in it."""
    n: Counter[str] = Counter()
    k: Counter[str] = Counter()
    for p in pairs:
        s = strata[screen_key(p)]
        n[s] += 1
        k[s] += int(p.regulated)
    return {s: k[s] / n[s] for s in n}


def screen_offsets(
    pairs: list[crispri.Pair], strata: dict[str, str], rates: dict[str, float], fallback: float
) -> list[float]:
    """The log odds of each pair's own screen's base rate, which enters the fit with coefficient 1."""
    return [logit(rates.get(strata.get(screen_key(p), SMALL_SCREEN), fallback)) for p in pairs]


def logistic_fit_offset(
    x: list[list[float]], y: list[bool], off: list[float], lam: float = 1e-3, rounds: int = 25
) -> list[float]:
    """`crispri.logistic_fit` with a per-row offset: a known constant added to the linear predictor.

    The offset spends no parameter, so a screen of six pairs can carry one and nothing separates.
    """
    k = len(x[0]) + 1
    rows = [[1.0, *r] for r in x]
    w = [0.0] * k
    for _ in range(rounds):
        grad = [0.0] * k
        hess = [[0.0] * k for _ in range(k)]
        for r, yi, o in zip(rows, y, off, strict=True):
            z = o + sum(a * b for a, b in zip(w, r, strict=True))
            p = sigmoid(z)
            g, h = (1.0 if yi else 0.0) - p, p * (1 - p)
            for i in range(k):
                grad[i] += g * r[i]
                hi = hess[i]
                for j in range(i, k):
                    hi[j] += h * r[i] * r[j]
        for i in range(k):
            for j in range(i):
                hess[i][j] = hess[j][i]
            if i:
                hess[i][i] += lam
                grad[i] -= lam * w[i]
        step = crispri.solve(hess, grad)
        w = [a + b for a, b in zip(w, step, strict=True)]
        if max(abs(s) for s in step) < 1e-6:
            break
    return w


def fit_with_offset(
    pairs: list[crispri.Pair], cols: tuple[str, ...], off: list[float], lam: float = 1e-3
) -> list[float]:
    return logistic_fit_offset(crispri.matrix(pairs, cols), [p.regulated for p in pairs], off, lam=lam)


def predict_with_offset(
    w: list[float], pairs: list[crispri.Pair], cols: tuple[str, ...], off: list[float]
) -> list[float]:
    z = crispri.logistic_score(w, crispri.matrix(pairs, cols))
    return [sigmoid(a + b) for a, b in zip(z, off, strict=True)]


def residual_shift(p: list[float], labels: list[bool]) -> float:
    """The log-odds constant still needed to match this set's own rate once the offset is in place."""
    n = len(labels)
    if not n:
        return 0.0
    observed = sum(labels) / n
    lo, hi = -12.0, 12.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if sum(sigmoid(logit(x) + mid) for x in p) / n < observed:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


SCREEN_BINS_FULL = 300  # a screen with fewer scored pairs than this is read in SMALL_BINS, not BINS


def screen_bins(n: int) -> int:
    return BINS if n >= SCREEN_BINS_FULL else SMALL_BINS


def one_screen(p: list[float], labels: list[bool]) -> dict[str, Any]:
    """One screen's reliability, at the registered bin count for its size, plus its residual shift."""
    bins = screen_bins(len(labels))
    cal = calibration(p, labels, bins)
    return {
        "pairs": len(labels),
        "regulated": sum(labels),
        "observed_rate": round(sum(labels) / len(labels), 4) if labels else None,
        "mean_predicted": cal["mean_predicted"],
        "bins": bins,
        "bins_consistent": cal["reliability"]["bins_consistent"],
        "ece": cal["reliability"]["ece"],
        "brier": cal["brier"],
        "auprc": cal["auprc"],
        "residual_log_odds_shift": round(residual_shift(p, labels), 4),
    }


def by_screen(rows: list[crispri.Pair], p: list[float], strata: dict[str, str]) -> dict[str, dict[str, Any]]:
    """The same probabilities split by the screen that produced the pair, each judged on its own size."""
    idx: dict[str, list[int]] = defaultdict(list)
    for i, pair in enumerate(rows):
        idx[strata.get(screen_key(pair), SMALL_SCREEN)].append(i)
    return {
        s: one_screen([p[i] for i in ii], [rows[i].regulated for i in ii])
        for s, ii in sorted(idx.items(), key=lambda kv: -len(kv[1]))
    }


def leave_one_screen_out(rows: list[crispri.Pair], cols: tuple[str, ...] = SWEEP_FEATURES) -> dict[str, Any]:
    """The registered test (E2): the shape fitted without a screen, the level taken from the screen.

    For each screen that carries its own term, the shape `b` is fitted on every other screen with
    their offsets in place, and the held-back screen is then scored with its own measured base rate
    as its offset and nothing else of its own. The comparator is the identical split with no offset
    anywhere. Nothing in the held-back screen's labels enters the shape; only the one number the
    screen reports about itself enters the level, which is the whole of the claim being tested.
    """
    census = screen_census(rows)
    own = [s for s, c in census.items() if c["own_term"]]
    per: dict[str, Any] = {}
    for s in own:
        held = [r for r in rows if screen_key(r) == s]
        rest = [r for r in rows if screen_key(r) != s]
        if not held or not rest or len(set(r.regulated for r in rest)) < 2:
            continue
        rest_strata = screen_strata(rest)
        rest_rates = screen_base_rates(rest, rest_strata)
        rest_pooled = sum(r.regulated for r in rest) / len(rest)
        w_off = fit_with_offset(rest, cols, screen_offsets(rest, rest_strata, rest_rates, rest_pooled))
        own_rate = census[s]["rate"]
        p_off = predict_with_offset(w_off, held, cols, [logit(own_rate)] * len(held))
        w_pool = fit(rest, cols)
        p_pool = predict(w_pool, held, cols)
        a = one_screen(p_pool, [r.regulated for r in held])
        b = one_screen(p_off, [r.regulated for r in held])
        per[s] = {
            "own_base_rate": round(own_rate, 4),
            "base_rate_of_the_other_screens": round(rest_pooled, 4),
            "pooled_model": a,
            "offset_model": b,
            "bins_improved": b["bins_consistent"] > a["bins_consistent"],
            "auprc_fell": bool(a["auprc"] is not None and b["auprc"] is not None and b["auprc"] < a["auprc"]),
        }
    total = sum(v["pooled_model"]["pairs"] for v in per.values())
    ece_pool = sum(v["pooled_model"]["ece"] * v["pooled_model"]["pairs"] for v in per.values())
    ece_off = sum(v["offset_model"]["ece"] * v["offset_model"]["pairs"] for v in per.values())
    residuals = sorted(abs(v["offset_model"]["residual_log_odds_shift"]) for v in per.values())
    median = (
        0.0
        if not residuals
        else (
            residuals[len(residuals) // 2]
            if len(residuals) % 2
            else (residuals[len(residuals) // 2 - 1] + residuals[len(residuals) // 2]) / 2
        )
    )
    improved = sum(1 for v in per.values() if v["bins_improved"])
    fell = round((ece_pool - ece_off) / total, 5) if total else 0.0
    flattened = any(v["auprc_fell"] for v in per.values())
    return {
        "what": (
            "leave one screen out: the shape fitted without the screen, the level taken from the "
            "screen's own base rate, against the same split with no offset anywhere"
        ),
        "screens_with_their_own_term": own,
        "per_screen": per,
        "screens_improving": [improved, len(per)],
        "screens_improving_required": SCREENS_IMPROVING_REQUIRED,
        "weighted_ece": [round(ece_pool / total, 5), round(ece_off / total, 5)] if total else None,
        "weighted_ece_fell_by": fell,
        "weighted_ece_fall_required": ECE_FALL_REQUIRED,
        "median_absolute_residual_shift": round(median, 4),
        "residual_shift_bar": RESIDUAL_SHIFT_BAR,
        "any_auprc_fell": flattened,
        "verdict": (
            "passed"
            if (
                improved >= SCREENS_IMPROVING_REQUIRED
                and fell >= ECE_FALL_REQUIRED
                and median < RESIDUAL_SHIFT_BAR
                and not flattened
            )
            else "failed"
        ),
    }


def smoothed(k: int, n: int) -> float:
    """(k + 0.5) / (n + 1): a rate a screen of 7 of 7 can still have a finite log odds for."""
    return (k + 0.5) / (n + 1)


def genome_under_each_screen(
    weights: dict[str, list[float]],
    shifts: dict[str, float],
    chroms: tuple[str, ...] = CHROMS,
    elements: Path = ELEMENTS,
    reference: Path = REFERENCE,
    results: Path = RESULTS,
) -> dict[str, Any]:
    """The published band table recomputed under one log-odds offset per screen, on the same walk.

    The published band is the `published` column, taken at a shift of zero on the same pass, so the
    table is its own control. Every other column is what the sweep would quote if the target it is
    quoting for belonged to that screen's population — which no sweep target does, and which is the
    point the columns are here to make visible.
    """
    counts: dict[str, dict[str, Counter[str]]] = {
        label: defaultdict(Counter) for label in ("any gene", "coding gene")
    }
    totals: Counter[str] = Counter()
    for chrom in chroms:
        p = elements / f"{chrom}.json"
        if not p.exists():
            continue
        els = json.loads(p.read_text())
        starts = gene_starts(chrom, reference)
        classes = registry_classes(chrom, results) if (results / f"ccres_{chrom}.bed.gz").exists() else {}
        for key, label in (("predicted", "any gene"), ("predicted_coding", "coding gene")):
            for el in els:
                row = element_row(el, key, starts, classes)
                if row is None or "missing" in row:
                    continue
                xt = [[row["features"][n] for n in TARGET_FEATURES]]
                z = crispri.logistic_score(weights["target"], xt)[0]
                totals[label] += 1
                counts[label]["published"][band_of(sigmoid(z))] += 1
                for name, d in shifts.items():
                    counts[label][name][band_of(sigmoid(z + d))] += 1
    return {
        label: {
            "targets": totals[label],
            "bands": {name: dict(sorted(c.items())) for name, c in counts[label].items()},
            "top_band_0.9_1": {name: c.get("0.9-1", 0) for name, c in counts[label].items()},
        }
        for label in ("any gene", "coding gene")
    }


def prevalence_terms(
    training: list[crispri.Pair],
    heldout: list[crispri.Pair],
    table: crispri.DeletionTable,
    reference: Path = REFERENCE,
    results: Path = RESULTS,
    chroms: tuple[str, ...] = CHROMS,
) -> dict[str, Any]:
    """The registered prevalence work: the census, the two reads, and what it leaves the genome."""
    crispri.annotate(training, table)
    crispri.annotate(heldout, table)
    add_features(training, table, reference, results)
    add_features(heldout, table, reference, results)
    train = scored(training)
    held = scored([p for p in heldout if p.cell == CELL])
    held_labels = [p.regulated for p in held]

    train_census, held_census = screen_census(train), screen_census(held)
    train_strata, held_strata = screen_strata(train), screen_strata(held)
    train_rates, held_rates = (
        screen_base_rates(train, train_strata),
        screen_base_rates(held, held_strata),
    )
    train_pooled = sum(p.regulated for p in train) / len(train)
    held_pooled = sum(held_labels) / len(held)

    shipped = predict(fit(train, SWEEP_FEATURES), held, SWEEP_FEATURES)
    before = calibration(shipped, held_labels, BINS)
    pooled_shift = prevalence_shift(shipped, held_labels, BINS)
    w_off = fit_with_offset(
        train, SWEEP_FEATURES, screen_offsets(train, train_strata, train_rates, train_pooled)
    )
    after_p = predict_with_offset(
        w_off, held, SWEEP_FEATURES, screen_offsets(held, held_strata, held_rates, held_pooled)
    )
    after = calibration(after_p, held_labels, BINS)

    both = train + held
    loso = leave_one_screen_out(both)

    on_gate = [p for p in both if p.features["top_target"]]
    gate_census = screen_census(on_gate)
    gate_total = len(on_gate)
    shifts = {}
    base = logit(smoothed(sum(p.regulated for p in on_gate), gate_total))
    for s, c in gate_census.items():
        if c["on_gate_pairs"]:
            shifts[s] = logit(smoothed(c["on_gate_regulated"], c["on_gate_pairs"])) - base
    weights = {"target": fit([p for p in train if p.features["top_target"]], TARGET_FEATURES)}
    genome = genome_under_each_screen(weights, shifts, chroms, ELEMENTS, reference, results)

    return {
        "evidence": EVIDENCE,
        "preregistered": PREREGISTERED_PREVALENCE,
        "census": {
            "what": "every screen in each population, counted before anything was fitted",
            "fitted (training K562)": {
                "pairs": len(train),
                "regulated": sum(p.regulated for p in train),
                "rate": round(train_pooled, 4),
                "screens": train_census,
            },
            "read (held-out K562)": {
                "pairs": len(held),
                "regulated": sum(held_labels),
                "rate": round(held_pooled, 4),
                "screens": held_census,
            },
            "screens_in_both_populations": sorted(set(train_census) & set(held_census)),
            "base_rate_strata": {
                "training": {k: round(v, 4) for k, v in sorted(train_rates.items())},
                "held_out": {k: round(v, 4) for k, v in sorted(held_rates.items())},
            },
            "on_the_gate": {
                "what": "the population whose curve bands the genome",
                "pairs": gate_total,
                "regulated": sum(p.regulated for p in on_gate),
                "screens_with_their_own_term_under_the_registered_rule": [
                    s for s, c in gate_census.items() if c["on_gate_pairs"] >= SCREEN_MIN_PAIRS
                ],
                "largest_screen_share": max(
                    (round(c["on_gate_pairs"] / gate_total, 4) for c in gate_census.values()),
                    default=None,
                ),
                "screens": {
                    s: {
                        "pairs": c["on_gate_pairs"],
                        "regulated": c["on_gate_regulated"],
                        "rate": c["on_gate_rate"],
                    }
                    for s, c in gate_census.items()
                    if c["on_gate_pairs"]
                },
            },
        },
        "held_out_before_and_after": {
            "what": (
                "E1, registered as a description that cannot fail: the offset takes each held-out "
                "screen's level from that screen's own held-out labels"
            ),
            "shipped curve": {
                "bins_consistent": before["reliability"]["bins_consistent"],
                "bins": before["reliability"]["bins"],
                "ece": before["reliability"]["ece"],
                "brier": before["brier"],
                "auprc": before["auprc"],
                "mean_predicted": before["mean_predicted"],
                "observed_rate": round(held_pooled, 4),
                "per_screen": by_screen(held, shipped, held_strata),
            },
            "one pooled shift (the published description)": pooled_shift,
            "per-screen offsets": {
                "bins_consistent": after["reliability"]["bins_consistent"],
                "bins": after["reliability"]["bins"],
                "ece": after["reliability"]["ece"],
                "brier": after["brier"],
                "auprc": after["auprc"],
                "mean_predicted": after["mean_predicted"],
                "per_screen": by_screen(held, after_p, held_strata),
            },
        },
        "leave_one_screen_out": loso,
        "verdict": loso["verdict"],
        "genome_under_each_screen": genome,
    }
