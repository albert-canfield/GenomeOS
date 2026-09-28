# SPDX-License-Identifier: AGPL-3.0-or-later
"""How large the measured experiment behind milestone 1.3's clause 2 would have to be, simulated for
the design that would actually be run, and registered before the first simulated block.

    uv run python scripts/clause2_design_power.py [--replicates N] [--no-save]

`scripts/clause2_measured_arm.py` ran the arm the project can afford today and could not decide:
3 of the 882 real-unknown blocks hold an element CRISPRi ever tested against a coding gene, 1 of them
also has a tested matched window, and the registered floor for reading an interval is 20 blocks. That
run also published a field called `coverage_needed` whose `n_for_80_percent_power` is 20 blocks. Two
things are wrong with reading that as the size of the experiment, and this script exists because of
them.

1. Both of that formula's inputs are **model output**. The effect `d = 0.2725` is the model arm's own
   matched difference in target-naming frequency and `s = 0.4324` is the standard deviation of the
   model arm's per-block differences. The experiment they were used to size has a **measured**
   endpoint: whether an assay finds that an element regulates a protein-coding gene. Nothing licenses
   carrying a model's effect size and a model's dispersion across to a different measurement.
2. The "19 blocks short" that was reported from it is `MIN_BLOCKS - compared_now` = 20 - 1, the
   distance to `measured.MIN_FOR_A_COMPARISON`, a **minimum-reporting floor**. The power formula
   happened to return 20 as well, which made a reporting rule look like a power result.

So the power question is asked here from the beginning, for the measured endpoint, over registered
grids, with the one empirical anchor the project holds and with everything else varied rather than
assumed. The output is a **range of sample sizes with its assumptions**, never a single number, and
power is a probability that an experiment detects an effect of an assumed size, never a guarantee
that it decides anything.

Nothing here re-opens, weakens or restates any committed figure. It computes no new descriptive
statistic about the blocks, and it makes 0 AlphaGenome requests: the model is not called, and the
model's -27.25 is deliberately **not** the effect being powered for.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from genomeos import manifest as mf
from genomeos.attribution import measured as ms
from genomeos.results import RESULTS_DIR, save_result

RESULT = "clause2_design_power"

#: the committed results this lane reads and never rewrites
MEASURED_ARM = "clause2_measured_arm"
MATCHED_CONTROL = "clause2_matched_control"
ELEMENT_COUNT_CONTROL = "clause2_element_count_control"
REACH_CONTROL = "clause2_reach_control"

#: the project's minimum-reporting floor, imported and not restated. It is a rule about when a number
#: may be printed, not a statement about power, and this script keeps the two apart everywhere.
MIN_BLOCKS = ms.MIN_FOR_A_COMPARISON

Z95 = 1.959964  # the two-sided 5% normal quantile, the committed estimator's own level
CHROMOSOMES = 24  # the blocks' clusters at the second level, as every clause 2 lane counted them

#: the sample sizes searched, in compared blocks (blocks with a tested element in BOTH arms)
N_GRID = (20, 30, 50, 75, 100, 150, 200, 300, 500, 750, 1000, 1500, 2000, 3000, 5000)
#: the powers reported for every configuration, so no single target hides the shape of the curve
POWER_TARGETS = (0.5, 0.8, 0.9)


PRE_REGISTRATION: dict[str, Any] = {
    "registered": "2026-09-28",
    "lane": "lane-design",
    "registered_before": (
        "any simulated block was drawn. The grids, their justifications, the endpoint, the estimator, "
        "the test, the clustering treatment, the reference configuration, the reported quantities and "
        "the words the answer may be stated in are all fixed in this dict, which was committed before "
        "the script was run"
    ),
    "question": (
        "for the experiment clause 2 actually needs -- a CRISPRi (or equivalently powered) screen that "
        "tests elements inside constrained-unknown blocks against protein-coding genes, against "
        "elements in matched windows -- how many compared blocks would be needed, over what range of "
        "assumptions, for the committed estimator to detect a difference of an assumed size"
    ),
    "endpoint": (
        "MEASURED, not modelled: per element, does the assay find a significant change in a "
        "protein-coding gene's expression on silencing that element (either sign, as R2 settled). Per "
        "block, the block's endpoint is 1 if at least one of its tested elements is found to regulate "
        "a coding gene. This is `clause2_measured_arm.PRIMARY_QUESTION`, unchanged"
    ),
    "estimator": (
        "the committed one, unchanged: d_i = Y_i - Wbar_i over compared blocks, where Y_i is the "
        "block's 0/1 endpoint and Wbar_i is the share of that block's TESTED matched windows whose "
        "endpoint is 1. The reported statistic is the mean of d_i in points"
    ),
    "test": (
        "the committed reading rule: the 95% interval over blocks excludes 0. Power is computed from "
        "the normal approximation to that interval (|mean| > 1.959964 * SE), because the committed "
        "interval is a percentile bootstrap of a mean over blocks and the two agree in the range that "
        "matters. That agreement is not assumed: a registered subset of cells is also run with the "
        "actual percentile bootstrap (2,000 resamples, 1,000 simulated experiments) and the largest "
        "disagreement in power is reported in the result. If it exceeds 0.05 the bootstrap numbers are "
        "the ones to read, and the result says so"
    ),
    "two_powers_reported": (
        "`power_any_direction` (the interval excludes 0 either way) and `power_correct_direction` (it "
        "excludes 0 on the side the simulated truth is on). The second is the one a reader wants and "
        "the one quoted"
    ),
    "the_floor_is_not_power": (
        f"`measured.MIN_FOR_A_COMPARISON` = {MIN_BLOCKS} blocks is a minimum-reporting rule: below it "
        "no interval is printed at all. It is applied here as a hard floor on the answer (an "
        "experiment of fewer than that many compared blocks reports nothing, whatever its power), and "
        "it is never reported as a power result. Where a power calculation returns a number below the "
        "floor, both are printed, separately named"
    ),
    "the_model_effect_is_excluded_by_name": (
        "the model arm's -27.25 points is NOT used as the effect to detect anywhere in this script. It "
        "is a difference in how often a deletion model names a target, measured on the same sequence; "
        "it is not an estimate of a difference in measured regulation, and there is no basis for "
        "treating it as one. The same goes for s = 0.4324, the model arm's per-block dispersion"
    ),
    "grids": {
        "window_true_rate": {
            "values_observed_scale": [0.0926, 0.1345, 0.1864],
            "justification": (
                "the one empirical anchor the project holds: 30 of 223 matched windows carrying a "
                "CRISPRi-tested element hold one measured to move a coding gene, 13.45%, from "
                "`clause2_measured_arm.json`. Its exact binomial 95% interval is 9.26% to 18.64% and "
                "that interval is the grid. Two independent anchors agree and are reported beside it "
                "rather than pooled in: over the whole ENCODE benchmark, 12.15% of the 3,941 distinct "
                "tested elements in the training file and 14.56% of the 1,697 in the held-out file "
                "have at least one significant pair. NOTE, and this is a correction to the write-up "
                "the anchor comes from: the 223 and the 30 are WINDOWS containing such an element, not "
                "tested elements; the element-level anchors are the ones that are element-level. The "
                "rate is carried on the OBSERVED scale, so the simulated true rate is "
                "observed / sensitivity, capped at 1"
            ),
        },
        "block_rate_ratio": {
            "values": [1.0, 0.75, 0.5, 0.25, 0.1, 0.0],
            "justification": (
                "the true per-element regulation rate inside the blocks as a multiple of the windows'. "
                "1.0 is the null the design must be able to sit at (the blocks regulate as often as "
                "their windows); 0.0 is the clause's strongest form (they regulate nothing); the rest "
                "span the range between. No value is taken from the model arm, and no published "
                "estimate of this ratio exists to take one from, which is exactly why a range is "
                "reported and not a number"
            ),
        },
        "assay_sensitivity": {
            "values": "the mean of each of the five power columns of the ENCODE benchmark files",
            "justification": (
                "`measured.POWER_COLUMNS` (PowerAtEffectSize10..50) state, per tested pair, the "
                "probability the screen would have called an effect of that size. Their means over the "
                "14,734 valid pairs are the sensitivities a real screen of this kind achieves, read "
                "from the project's own tables rather than assumed. They are reported by split as "
                "well, because the held-out file is better powered than the training file. The "
                "sequence is not monotone in effect size (the 25% column is the benchmark's own "
                "filter and sits near 1 for every pair); the columns are used as five separate "
                "sensitivity levels and the non-monotonicity is reported, not smoothed"
            ),
        },
        "false_positive_rate": {
            "values": [0.0, 0.01],
            "justification": (
                "the benchmark controls its own false discovery rate, so 0 is the design's intent and "
                "0.01 per tested element is the sensitivity analysis. A false positive rate acts on "
                "both arms and mostly costs power by flattening the contrast"
            ),
        },
        "tested_elements_per_block": {
            "values": [1, 2, 6],
            "justification": (
                "6 is 6.177, the scored elements per carrying block in the committed matched control, "
                "so it is what a screen that tested every scored element in a block would have; 1 is "
                "the minimum that makes a block testable and is what the 3 tested blocks that exist "
                "today actually have; 2 is between them. The same count is used per tested window"
            ),
        },
        "tested_windows_per_block": {
            "values": [1, 3, 10],
            "justification": (
                "the matched control draws 50 windows per block, of which only some would be tested. "
                "Today 223 windows over 134 blocks carry a tested element, about 1.7 each, so 1 is the "
                "realistic low end; 3 and 10 say what buying more control measurements is worth. More "
                "windows shrink the variance of Wbar_i and therefore of d_i"
            ),
        },
        "elements_within_block_icc": {
            "values": [0.0, 0.21, 0.30, 0.35],
            "justification": (
                "measured, not guessed. Whether an element is found to regulate a coding gene is "
                "correlated between nearby elements: the one-way ANOVA estimator on the benchmark's "
                "5,638 distinct tested elements gives an intraclass correlation of 0.375 within 25 kb "
                "clusters, 0.350 within 50 kb, 0.300 within 100 kb, 0.273 within 250 kb and 0.209 "
                "within 1 Mb. A constrained-unknown block is an intergenic gap of that order, so the "
                "grid is 0 (independence, for contrast), 0.21, 0.30 and 0.35. The script recomputes "
                "these from the tables at run time and reports them"
            ),
        },
        "blocks_within_chromosome_icc": {
            "values": [0.0, 0.014, 0.03],
            "justification": (
                "the same estimator with chromosomes as clusters gives 0.0142 over the 23 chromosomes "
                "the benchmark touches. 0 and 0.03 bracket it"
            ),
        },
    },
    "clustering_treatment": {
        "elements_within_a_block": (
            "HIERARCHICAL SIMULATION. Each block carries a random intercept u ~ Normal(0, sigma) on "
            "the logit of its elements' regulation probability, with sigma fixed from the registered "
            "ICC by the latent-scale identity sigma^2 / (sigma^2 + pi^2/3) = ICC. Its elements are "
            "then conditionally independent given u. Each tested window carries its own independent "
            "intercept, because a window is somewhere else in the genome. The realised binary-scale "
            "ICC of the simulated elements is reported beside the latent one, since they differ"
        ),
        "blocks_within_a_chromosome": (
            "DESIGN EFFECT, applied to the variance of the mean and not simulated: "
            "Var(mean d) = sigma_d^2 / N * (1 + (N / 24 - 1) * rho_chrom), the standard equal-cluster "
            "design effect with the N blocks spread over the 24 chromosomes. This is deliberate and it "
            "is the conservative direction to be explicit about: the committed interval is a bootstrap "
            "over BLOCKS, which assumes blocks are independent, so where rho_chrom > 0 the committed "
            "test is anti-conservative and the design effect is what says by how much"
        ),
    },
    "reference_configuration": {
        "window_observed_rate": 0.1345,
        "block_rate_ratio": "swept",
        "assay_sensitivity": "PowerAtEffectSize20, pooled mean",
        "false_positive_rate": 0.0,
        "tested_elements_per_block": 6,
        "tested_windows_per_block": 3,
        "elements_within_block_icc": 0.30,
        "blocks_within_chromosome_icc": 0.014,
        "justification": (
            "each field is the middle or the measured value of its own grid, so the sensitivity "
            "analysis varies one factor at a time from a configuration that is defensible on its own"
        ),
    },
    "sweep": (
        "the primary sweep is the block-rate ratio at the reference configuration, over N_GRID. The "
        "sensitivity analysis is one factor at a time: every other value of every grid, at each ratio, "
        "over the same N_GRID. Cells are not crossed exhaustively, because a fully crossed grid would "
        "report a range whose extremes are combinations no one would defend"
    ),
    "how_the_moments_are_obtained": (
        "for each configuration, 200,000 simulated blocks with their windows are drawn once and the "
        "mean and standard deviation of d_i are taken from them; power at every N then follows in "
        "closed form from the normal approximation and the chromosome design effect. This is a "
        "Monte Carlo estimate of the estimator's own moments, not a Monte Carlo of the test, and its "
        "own error is reported (the standard error of the simulated mean)"
    ),
    "reported": (
        "for every configuration: the smallest N in N_GRID reaching 50%, 80% and 90% correct-direction "
        "power, or 'above the largest N searched'; the floor applied separately; and across the "
        "configurations, the RANGE of the 80% number at each ratio. The headline is a range with its "
        "assumptions attached, and the result states in its own text that power is a probability of "
        "detection under assumptions, not a guarantee that the experiment decides clause 2"
    ),
    "equivalence_decision_registered_before_computing": (
        "for the neutral tier against the real unknown (committed contrast -0.10 points, 95% -4.55 to "
        "+4.30 on `names_a_coding_gene`; -0.71, -5.58 to +4.17 element-count-matched), an equivalence "
        "claim needs a margin justified from OUTSIDE these data. The rule fixed here, before looking: "
        "a margin is admissible only if it comes from a source that is not this project's own analysis "
        "of these blocks -- a published threshold for a meaningful difference in regulatory annotation "
        "rate, a decision-theoretic cost, or a downstream requirement stated elsewhere in the project "
        "before this lane. The project's own 5-point chance band is NOT admissible: it was chosen in "
        "d717b28 for this very comparison. If no admissible margin exists, no equivalence test is run "
        "and the finding is stated as 'no difference detected'; the script then reports, descriptively "
        "and labelled as such, the smallest margin at which a TOST on the committed interval would "
        "pass -- which is a statement about what the data would need, not a margin anyone has "
        "justified"
    ),
    "denominator_table": (
        "every population used across the four clause 2 lanes, each count with its denominator and the "
        "rule that produced it, read from the committed result files rather than retyped, including "
        "the overlap-rule sensitivity at 0.25, 0.5 and 0.75"
    ),
    "falsifier": (
        "if the registered bootstrap check disagrees with the normal approximation by more than 0.05 "
        "in power at any checked cell, the normal-approximation table is reported as unreliable and "
        "the bootstrap numbers replace it in the reading"
    ),
    "cost": "0 AlphaGenome requests; no model is called and no per-element cache is opened",
}


# ---- the empirical anchors, measured from the project's own tables ---------------------------------


def _icc_oneway(groups: list[list[int]]) -> dict[str, Any] | None:
    """One-way random-effects ICC of a 0/1 outcome over unequal clusters (ANOVA estimator)."""
    groups = [g for g in groups if len(g) >= 2]
    k = len(groups)
    n = sum(len(g) for g in groups)
    if k < 2 or n <= k:
        return None
    flat = [y for g in groups for y in g]
    ybar = sum(flat) / n
    msb = sum(len(g) * (sum(g) / len(g) - ybar) ** 2 for g in groups) / (k - 1)
    msw = sum(sum((y - sum(g) / len(g)) ** 2 for y in g) for g in groups) / (n - k)
    m0 = (n - sum(len(g) ** 2 for g in groups) / n) / (k - 1)
    denom = msb + (m0 - 1) * msw
    return {
        "clusters": k,
        "elements": n,
        "mean_cluster_size": round(m0, 3),
        "icc": round((msb - msw) / denom, 4) if denom > 0 else None,
    }


def anchors() -> dict[str, Any]:
    """Every number this simulation is anchored on, recomputed from the tables it came from."""
    pairs, invalid = ms.load_crispri()
    by_element: dict[tuple[str, int, int], list[ms.CrispriPair]] = defaultdict(list)
    for p in pairs:
        by_element[(p.chrom, p.start, p.end)].append(p)
    els = [(c, s, e, int(any(x.significant for x in v))) for (c, s, e), v in by_element.items()]

    power = {}
    for col in ms.POWER_COLUMNS:
        attr = "power_at_effect_size_" + col.removeprefix("PowerAtEffectSize")
        vals = [getattr(p, attr) for p in pairs if getattr(p, attr) is not None]
        by_split = {}
        for split in (ms.TRAINING, ms.HELDOUT):
            sv = [getattr(p, attr) for p in pairs if p.split == split and getattr(p, attr) is not None]
            by_split[split] = round(sum(sv) / len(sv), 4) if sv else None
        power[col] = {
            "pairs": len(vals),
            "mean": round(sum(vals) / len(vals), 4) if vals else None,
            "share_at_or_above_0.8": round(sum(v >= 0.8 for v in vals) / len(vals), 4) if vals else None,
            "mean_by_split": by_split,
        }

    element_rate = {}
    for split in (ms.TRAINING, ms.HELDOUT, ms.ALL):
        sub = [(c, s, e, v) for (c, s, e), v in by_element.items() if split in (ms.ALL, v[0].split)]
        moved = sum(1 for *_, v in sub if any(x.significant for x in v))
        element_rate[split] = {
            "tested_elements": len(sub),
            "with_a_significant_pair": moved,
            "rate": round(moved / len(sub), 4) if sub else None,
        }

    within = {}
    for bp in (25_000, 50_000, 100_000, 250_000, 1_000_000):
        g: dict[tuple[str, int], list[int]] = defaultdict(list)
        for c, s, e, y in els:
            g[(c, (s + e) // 2 // bp)].append(y)
        within[f"{bp // 1000}kb"] = _icc_oneway(list(g.values()))
    per_chrom: dict[str, list[int]] = defaultdict(list)
    for c, _s, _e, y in els:
        per_chrom[c].append(y)

    return {
        "crispri_pairs": len(pairs),
        "crispri_invalid_pairs": invalid,
        "power_columns": power,
        "power_columns_are_not_monotone_in_effect_size": (
            "PowerAtEffectSize25 sits near 1 for essentially every pair because it is the benchmark's "
            "own inclusion filter, so the column means do not increase with effect size. The five "
            "columns are used as five sensitivity levels, in the order they are reported, not as a "
            "curve"
        ),
        "element_level_regulation_rate": element_rate,
        "elements_within_cluster_icc": within,
        "blocks_within_chromosome_icc": _icc_oneway(list(per_chrom.values())),
    }


def window_anchor_from_the_measured_arm(results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """The 30-of-223 anchor, read from the committed result rather than retyped, with what it is."""
    p = results_dir / f"{MEASURED_ARM}.json"
    if not p.exists():
        return {}
    w = json.loads(p.read_text())["coverage"]["windows_matched_real_unknown"]
    tested = w["windows_with"]["crispri_tested_against_a_coding_gene"]
    moved = w["windows_with"]["measured_to_regulate_a_coding_gene"]
    lo, hi = _exact_binomial_interval(moved, tested)
    return {
        "numerator": moved,
        "denominator": tested,
        "rate": round(moved / tested, 4) if tested else None,
        "exact_95_interval": [round(lo, 4), round(hi, 4)],
        "unit": "matched windows containing at least one CRISPRi-tested element",
        "correction": (
            "the lane-measured write-up calls these 'tested elements'. They are windows: the count "
            "comes from `windows_with`, which counts windows satisfying the predicate. Windows of the "
            "same block can contain the same tested element, so the 223 are not 223 independent "
            "measurements and the exact interval above is narrower than the truth. The element-level "
            "rates over the benchmark's own distinct elements are the element-level anchors"
        ),
    }


def _exact_binomial_interval(k: int, n: int) -> tuple[float, float]:
    """Clopper-Pearson, without scipy: the beta quantiles by bisection on the regularised beta."""
    if n == 0:
        return (0.0, 1.0)

    def betainc(a: float, b: float, x: float) -> float:
        # regularised incomplete beta by the continued fraction, enough digits for an interval
        if x <= 0:
            return 0.0
        if x >= 1:
            return 1.0
        lbeta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
        front = math.exp(math.log(x) * a + math.log(1 - x) * b - lbeta) / a
        f, c, d = 1.0, 1.0, 0.0
        for i in range(0, 300):
            m = i // 2
            if i == 0:
                num = 1.0
            elif i % 2 == 0:
                num = (m * (b - m) * x) / ((a + 2 * m - 1) * (a + 2 * m))
            else:
                num = -((a + m) * (a + b + m) * x) / ((a + 2 * m) * (a + 2 * m + 1))
            d = 1.0 + num * d
            d = 1e-30 if abs(d) < 1e-30 else d
            d = 1.0 / d
            c = 1.0 + num / (c if abs(c) > 1e-30 else 1e-30)
            f *= c * d
            if abs(1.0 - c * d) < 1e-12:
                break
        got = front * (f - 1.0)
        return got if a >= (a + b) * x else 1.0 - betainc(b, a, 1 - x)

    def solve(target: float, a: float, b: float) -> float:
        lo, hi = 0.0, 1.0
        for _ in range(200):
            mid = (lo + hi) / 2
            if betainc(a, b, mid) < target:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2

    low = 0.0 if k == 0 else solve(0.025, k, n - k + 1)
    high = 1.0 if k == n else solve(0.975, k + 1, n - k)
    return (low, high)


# ---- the simulation --------------------------------------------------------------------------------


def sigma_for_icc(icc: float) -> float:
    """The logit-normal intercept standard deviation whose latent-scale ICC is `icc`."""
    if icc <= 0:
        return 0.0
    return math.sqrt(icc / (1 - icc) * math.pi**2 / 3)


def _expit(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-x))


def moments(
    window_observed_rate: float,
    ratio: float,
    sensitivity: float,
    false_positive: float,
    elements_per_block: int,
    windows_per_block: int,
    icc: float,
    blocks: int,
    rng: np.random.Generator,
) -> dict[str, Any]:
    """Mean and standard deviation of the per-block difference under one configuration.

    Hierarchical: every block and every window carries its own logit-normal intercept and its
    elements are conditionally independent given it. The window arm's true rate is set so that the
    OBSERVED window rate matches the anchor at this sensitivity.
    """
    p_win_true = min(1.0, window_observed_rate / sensitivity) if sensitivity > 0 else 0.0
    p_block_true = min(1.0, p_win_true * ratio)
    sigma = sigma_for_icc(icc)

    def arm(p_true: float, n_clusters: int) -> tuple[np.ndarray, np.ndarray]:
        """(endpoint per cluster, per-element detection probability) for `n_clusters` clusters."""
        if p_true <= 0:
            q = np.full(n_clusters, false_positive)
        else:
            lo = math.log(p_true / (1 - p_true)) if p_true < 1 else 30.0
            u = rng.normal(0.0, sigma, size=n_clusters) if sigma else np.zeros(n_clusters)
            p = _expit(lo + u)
            q = p * sensitivity + (1 - p) * false_positive
        hits = rng.binomial(elements_per_block, q)
        return (hits > 0).astype(float), q

    y, q_block = arm(p_block_true, blocks)
    w, _ = arm(p_win_true, blocks * windows_per_block)
    wbar = w.reshape(blocks, windows_per_block).mean(axis=1)
    d = y - wbar
    mean = float(d.mean())
    sd = float(d.std(ddof=1))
    return {
        "window_true_rate": round(p_win_true, 4),
        "block_true_rate": round(p_block_true, 4),
        "block_endpoint_rate": round(float(y.mean()), 4),
        "window_endpoint_rate": round(float(w.mean()), 4),
        "mean_difference_points": round(100 * mean, 3),
        "sd_of_per_block_difference": round(sd, 4),
        "monte_carlo_se_of_the_mean_points": round(100 * sd / math.sqrt(blocks), 4),
        "simulated_blocks": blocks,
        "latent_icc": icc,
        "realised_element_detection_probability": round(float(q_block.mean()), 5),
    }


def power_at(mean: float, sd: float, n: int, rho_chrom: float, chroms: int = CHROMOSOMES) -> dict[str, float]:
    """Two-sided and correct-direction power of the committed interval at `n` compared blocks."""
    if sd <= 0:
        return {"any_direction": 1.0 if mean != 0 else 0.0, "correct_direction": 1.0 if mean != 0 else 0.0}
    deff = 1.0 + (max(n / chroms, 1.0) - 1.0) * rho_chrom
    se = sd / math.sqrt(n) * math.sqrt(deff)
    lam = mean / se
    upper = _normal_sf(Z95 - lam)  # P(estimate above the upper critical value)
    lower = _normal_sf(Z95 + lam)  # P(estimate below the lower one)
    correct = upper if mean > 0 else lower
    return {
        "any_direction": round(upper + lower, 4),
        "correct_direction": round(correct, 4),
        "design_effect": round(deff, 4),
    }


def _normal_sf(x: float) -> float:
    return 0.5 * math.erfc(x / math.sqrt(2))


def smallest_n(mean: float, sd: float, rho_chrom: float, target: float) -> int | None:
    """The smallest N in N_GRID whose correct-direction power reaches `target`, or None."""
    for n in N_GRID:
        if power_at(mean, sd, n, rho_chrom)["correct_direction"] >= target:
            return n
    return None


def bootstrap_check(
    mean_cfg: dict[str, Any],
    cfg: dict[str, Any],
    n: int,
    experiments: int,
    resamples: int,
    rng: np.random.Generator,
) -> dict[str, Any]:
    """Run the committed percentile bootstrap on simulated experiments of `n` blocks, and compare.

    The registered check on the normal approximation: this is the actual reading rule of the four
    clause 2 lanes, applied to data simulated under the same configuration, with no approximation at
    all. The comparison is against the same configuration's normal-approximation power.
    """
    truth_positive = mean_cfg["mean_difference_points"] > 0
    hits = correct = 0
    for _ in range(experiments):
        y, w = _draw_arms(cfg, n, rng)
        diff = y - w
        idx = rng.integers(0, n, size=(resamples, n))
        boots = diff[idx].mean(axis=1)
        lo, hi = np.percentile(boots, [2.5, 97.5])
        if lo > 0 or hi < 0:
            hits += 1
            if (lo > 0) == truth_positive:
                correct += 1
    mean, sd = mean_cfg["mean_difference_points"] / 100, mean_cfg["sd_of_per_block_difference"]
    normal = power_at(mean, sd, n, cfg["rho_chrom"])
    out = {
        "n_blocks": n,
        "experiments": experiments,
        "resamples": resamples,
        "bootstrap_power_any_direction": round(hits / experiments, 4),
        "bootstrap_power_correct_direction": round(correct / experiments, 4),
        "normal_power_any_direction": normal["any_direction"],
        "normal_power_correct_direction": normal["correct_direction"],
    }
    out["largest_disagreement"] = round(
        max(
            abs(out["bootstrap_power_any_direction"] - out["normal_power_any_direction"]),
            abs(out["bootstrap_power_correct_direction"] - out["normal_power_correct_direction"]),
        ),
        4,
    )
    return out


def _draw_arms(cfg: dict[str, Any], n: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """One experiment's raw per-block endpoints and window shares."""
    s = cfg["sensitivity"]
    p_win = min(1.0, cfg["window_observed_rate"] / s) if s > 0 else 0.0
    p_blk = min(1.0, p_win * cfg["ratio"])
    sigma = sigma_for_icc(cfg["icc"])

    def arm(p_true: float, clusters: int) -> np.ndarray:
        if p_true <= 0:
            q = np.full(clusters, cfg["false_positive"])
        else:
            lo = math.log(p_true / (1 - p_true)) if p_true < 1 else 30.0
            u = rng.normal(0.0, sigma, size=clusters) if sigma else np.zeros(clusters)
            q = _expit(lo + u) * s + (1 - _expit(lo + u)) * cfg["false_positive"]
        return (rng.binomial(cfg["elements_per_block"], q) > 0).astype(float)

    y = arm(p_blk, n)
    w = arm(p_win, n * cfg["windows_per_block"]).reshape(n, cfg["windows_per_block"]).mean(axis=1)
    return y, w


# ---- the registered sweep --------------------------------------------------------------------------


def reference(anch: dict[str, Any], window_rate: float) -> dict[str, Any]:
    """The reference configuration of PRE_REGISTRATION, with its sensitivity read from the tables."""
    return {
        "window_observed_rate": window_rate,
        "ratio": 1.0,
        "sensitivity": anch["power_columns"]["PowerAtEffectSize20"]["mean"],
        "sensitivity_column": "PowerAtEffectSize20",
        "false_positive": 0.0,
        "elements_per_block": 6,
        "windows_per_block": 3,
        "icc": 0.30,
        "rho_chrom": 0.014,
    }


def variants(anch: dict[str, Any], window_rate: float) -> list[tuple[str, dict[str, Any]]]:
    """One factor at a time from the reference, each with the grid value it changes to."""
    ref = reference(anch, window_rate)
    out: list[tuple[str, dict[str, Any]]] = [("reference", dict(ref))]
    for col in ms.POWER_COLUMNS:
        m = anch["power_columns"][col]["mean"]
        if col != ref["sensitivity_column"]:
            out.append((f"sensitivity={col}", {**ref, "sensitivity": m, "sensitivity_column": col}))
    for k in (1, 2, 6):
        if k != ref["elements_per_block"]:
            out.append((f"elements_per_block={k}", {**ref, "elements_per_block": k}))
    for m in (1, 3, 10):
        if m != ref["windows_per_block"]:
            out.append((f"windows_per_block={m}", {**ref, "windows_per_block": m}))
    for icc in (0.0, 0.21, 0.30, 0.35):
        if icc != ref["icc"]:
            out.append((f"element_icc={icc}", {**ref, "icc": icc}))
    for rho in (0.0, 0.014, 0.03):
        if rho != ref["rho_chrom"]:
            out.append((f"chromosome_icc={rho}", {**ref, "rho_chrom": rho}))
    lo, hi = PRE_REGISTRATION["grids"]["window_true_rate"]["values_observed_scale"][0::2]
    for w in (lo, hi):
        out.append((f"window_observed_rate={w}", {**ref, "window_observed_rate": w}))
    out.append(("false_positive=0.01", {**ref, "false_positive": 0.01}))
    return out


def sweep(anch: dict[str, Any], window_rate: float, blocks: int, seed: int) -> dict[str, Any]:
    """Every registered configuration at every registered ratio, and the N each one needs."""
    rng = np.random.default_rng(seed)
    ratios = PRE_REGISTRATION["grids"]["block_rate_ratio"]["values"]
    rows: list[dict[str, Any]] = []
    for name, cfg in variants(anch, window_rate):
        for ratio in ratios:
            c = {**cfg, "ratio": ratio}
            m = moments(
                c["window_observed_rate"],
                ratio,
                c["sensitivity"],
                c["false_positive"],
                c["elements_per_block"],
                c["windows_per_block"],
                c["icc"],
                blocks,
                rng,
            )
            mean, sd = m["mean_difference_points"] / 100, m["sd_of_per_block_difference"]
            row = {
                "configuration": name,
                "ratio": ratio,
                "settings": {k: v for k, v in c.items() if k != "ratio"},
                **m,
                "n_for_power": {str(t): smallest_n(mean, sd, c["rho_chrom"], t) for t in POWER_TARGETS},
                "power_at_the_reporting_floor": power_at(mean, sd, MIN_BLOCKS, c["rho_chrom"]),
                "power_at_100_blocks": power_at(mean, sd, 100, c["rho_chrom"]),
                "power_at_500_blocks": power_at(mean, sd, 500, c["rho_chrom"]),
            }
            rows.append(row)
    return {"configurations": len(variants(anch, window_rate)), "ratios": list(ratios), "rows": rows}


def ranges(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """The answer the lane owes: a RANGE of compared blocks per ratio, with what drives each end."""
    out: dict[str, Any] = {}
    for ratio in sorted({r["ratio"] for r in rows}, reverse=True):
        at = [r for r in rows if r["ratio"] == ratio]
        got = [(r["n_for_power"]["0.8"], r["configuration"]) for r in at]
        found = [(n, c) for n, c in got if n is not None]
        out[str(ratio)] = {
            "configurations": len(at),
            "reached_80_percent_within_the_searched_range": len(found),
            "smallest": min(found)[0] if found else None,
            "smallest_under": min(found)[1] if found else None,
            "largest": max(found)[0] if found else None,
            "largest_under": max(found)[1] if found else None,
            "did_not_reach_80_percent_by_n": [c for n, c in got if n is None],
            "reference_configuration": next(
                (r["n_for_power"]["0.8"] for r in at if r["configuration"] == "reference"), None
            ),
        }
    return out


# ---- equivalence -----------------------------------------------------------------------------------


EQUIVALENCE_MARGIN_SEARCH = (
    "sources looked at for a margin, before any computation: (1) published thresholds for a "
    "'meaningful' difference in the rate at which sequence is annotated as regulatory -- none states "
    "one; the ENCODE benchmark and the enhancer-gene literature report effect sizes and power, not "
    "equivalence bounds; (2) a decision-theoretic cost, which would need a stated consequence of "
    "treating the neutral tier as the real unknown -- the project states none; (3) an earlier "
    "requirement in this project fixed before this comparison -- the only candidate is the 5-point "
    "chance band of d717b28, which was chosen inside this analysis for this comparison and is ruled "
    "inadmissible by the registration above. No admissible margin was found"
)


def equivalence(results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """The registered decision for 'the neutral tier behaves identically', applied."""
    out: dict[str, Any] = {
        "claim_reviewed": "the neutral tier and the real unknown are indistinguishable / behave identically",
        "margin_search": EQUIVALENCE_MARGIN_SEARCH,
        "admissible_margin_found": False,
        "test_run": None,
        "finding": (
            "NO DIFFERENCE DETECTED. Not equivalence, and not 'identical'. No equivalence test is run, "
            "because no margin can be justified from outside these data, and the registration fixed "
            "that rule before the contrasts were read"
        ),
        "contrasts": [],
    }
    for name, key in (
        (MATCHED_CONTROL, "secondary_real_minus_neutral"),
        (ELEMENT_COUNT_CONTROL, "secondary_real_minus_neutral"),
    ):
        p = results_dir / f"{name}.json"
        if not p.exists():
            continue
        for q, c in json.loads(p.read_text())[key].items():
            lo, hi = c["ci95"]
            out["contrasts"].append(
                {
                    "result": name,
                    "question": q,
                    "points": c["points"],
                    "ci95": [lo, hi],
                    "smallest_margin_a_tost_would_pass_at_points": round(max(abs(lo), abs(hi)), 2),
                    "what_the_interval_permits": (
                        f"a real-unknown rate anywhere from {abs(lo):.2f} points below the neutral "
                        f"tier's to {abs(hi):.2f} points above it. Differences in either direction of "
                        "that size are compatible with the data"
                    ),
                }
            )
    out["smallest_margin_note"] = (
        "`smallest_margin_a_tost_would_pass_at_points` is the margin the DATA would need for a "
        "two-one-sided-test to declare equivalence, computed by inverting the committed 95% interval "
        "(conservative: a TOST at 5% uses the 90% interval, which is narrower). It is a description of "
        "what these data could support, NOT a margin anyone has justified, and it must never be quoted "
        "as one. Reading it as a margin would be choosing the bound after seeing the estimate"
    )
    return out


# ---- the denominator table -------------------------------------------------------------------------


def denominators(results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """Every population the four clause 2 lanes count over, with its denominator and its rule."""

    def read(name: str) -> dict[str, Any]:
        p = results_dir / f"{name}.json"
        return json.loads(p.read_text()) if p.exists() else {}

    arm, matched, elem, reach = (
        read(MEASURED_ARM),
        read(MATCHED_CONTROL),
        read(ELEMENT_COUNT_CONTROL),
        read(REACH_CONTROL),
    )
    rows: list[dict[str, Any]] = []

    def row(population: str, count: Any, denominator: str, of: Any, rule: str, source: str) -> None:
        share = round(count / of, 4) if isinstance(count, int) and isinstance(of, int) and of else None
        rows.append(
            {
                "population": population,
                "count": count,
                "denominator": denominator,
                "denominator_count": of,
                "share": share,
                "rule": rule,
                "source": source,
            }
        )

    if arm:
        ru = arm["coverage"]["real_unknown"]
        nt = arm["coverage"]["neutral"]
        row(
            "constrained-unknown blocks, copies out (the real unknown)",
            ru["blocks"],
            "itself: the tier",
            ru["blocks"],
            "organise.blocks on 24 chromosomes, the `constrained_unknown` tier with duplicate copies removed",
            f"{MEASURED_ARM}.json coverage.real_unknown.blocks",
        )
        row(
            "real-unknown blocks carrying at least one scored element",
            ru["blocks_carrying_a_scored_element"],
            "real-unknown blocks",
            ru["blocks"],
            "at least one ENCODE SCREEN cCRE in the all-element archive has its midpoint inside the block",
            f"{MEASURED_ARM}.json coverage.real_unknown",
        )
        row(
            "scored elements inside real-unknown blocks",
            ru["scored_elements"],
            "elements, not blocks",
            ru["scored_elements"],
            "every scored element of a carrying block; the per-element estimator's denominator",
            f"{MEASURED_ARM}.json coverage.real_unknown.scored_elements",
        )
        for label, key in (
            ("eligible for an assay (shares one base with a tested interval)", "eligible_for_an_assay"),
            ("measured by an element-level assay", "measured_by_an_element_level_assay"),
            ("CRISPRi-tested against a protein-coding gene", "crispri_tested_against_a_coding_gene"),
            ("measured to move a coding gene", "measured_to_regulate_a_coding_gene"),
        ):
            row(
                f"real-unknown blocks {label}",
                ru["blocks_with"][key],
                "real-unknown blocks",
                ru["blocks"],
                "reciprocal overlap 0.5 both ways, except eligibility which is one shared base",
                f"{MEASURED_ARM}.json coverage.real_unknown.blocks_with",
            )
            row(
                f"real-unknown blocks {label}, against the carrying blocks",
                ru["blocks_with"][key],
                "real-unknown blocks carrying a scored element",
                ru["blocks_carrying_a_scored_element"],
                "same counts, the other denominator: this is where 6.7% and 11.1% differ",
                f"{MEASURED_ARM}.json coverage.real_unknown.blocks_with",
            )
        row(
            "real-unknown blocks with no measured element at all",
            ru["blocks_with_no_measured_element"],
            "real-unknown blocks",
            ru["blocks"],
            "no assay measures ANY of the block's scored elements at reciprocal overlap 0.5. It does "
            "NOT mean no part of the block was ever measured, and it does not speak about a block "
            "that carries no scored element",
            f"{MEASURED_ARM}.json coverage.real_unknown.blocks_with_no_measured_element",
        )
        row(
            "neutral-tier blocks",
            nt["blocks"],
            "itself: the tier",
            nt["blocks"],
            "the secondary target set, same construction",
            f"{MEASURED_ARM}.json coverage.neutral.blocks",
        )
        row(
            "neutral-tier blocks carrying a scored element",
            nt["blocks_carrying_a_scored_element"],
            "neutral-tier blocks",
            nt["blocks"],
            "as above",
            f"{MEASURED_ARM}.json coverage.neutral",
        )
        w = arm["coverage"]["windows_matched_real_unknown"]
        row(
            "matched windows drawn",
            w["windows_drawn"],
            "windows",
            w["windows_drawn"],
            "50 accepted windows per block of the block's exact length, coding-TSS decile matched, "
            "20,000 tries, seed 20260913",
            f"{MEASURED_ARM}.json coverage.windows_matched_real_unknown",
        )
        row(
            "matched windows carrying a scored element",
            w["windows_carrying_a_scored_element"],
            "matched windows drawn",
            w["windows_drawn"],
            "at least one scored element's midpoint inside the window",
            f"{MEASURED_ARM}.json coverage.windows_matched_real_unknown",
        )
        row(
            "matched windows carrying a CRISPRi-tested element",
            w["windows_with"]["crispri_tested_against_a_coding_gene"],
            "matched windows drawn",
            w["windows_drawn"],
            "WINDOWS, not elements: several windows of one block can carry the same tested element, "
            "so these are not independent measurements",
            f"{MEASURED_ARM}.json coverage.windows_matched_real_unknown.windows_with",
        )
        row(
            "matched windows holding an element measured to move a coding gene",
            w["windows_with"]["measured_to_regulate_a_coding_gene"],
            "matched windows carrying a CRISPRi-tested element",
            w["windows_with"]["crispri_tested_against_a_coding_gene"],
            "the 13.45% anchor. Its unit is windows; the lane-measured write-up calls it elements",
            f"{MEASURED_ARM}.json coverage.windows_matched_real_unknown.windows_with",
        )
        sens = arm["coverage"].get("sensitivity_by_overlap", {})
        for f, v in sorted(sens.items(), key=lambda kv: float(kv[0])):
            row(
                f"real-unknown blocks with a measured element at reciprocal overlap {f}",
                v["real_unknown"]["blocks_with"]["measured_by_an_element_level_assay"],
                "real-unknown blocks",
                ru["blocks"],
                f"the overlap rule set to {f} instead of 0.5; the count is a property of the rule",
                f"{MEASURED_ARM}.json coverage.sensitivity_by_overlap",
            )
            row(
                f"real-unknown blocks CRISPRi-tested against a coding gene at overlap {f}",
                v["real_unknown"]["blocks_with"]["crispri_tested_against_a_coding_gene"],
                "real-unknown blocks",
                ru["blocks"],
                f"the overlap rule set to {f}; no setting of it produces a usable arm",
                f"{MEASURED_ARM}.json coverage.sensitivity_by_overlap",
            )
        row(
            "blocks the measured primary could compare",
            arm["primary"]["n_blocks_compared"],
            "real-unknown blocks",
            ru["blocks"],
            "block AND at least one of its matched windows both carry a CRISPRi-tested element",
            f"{MEASURED_ARM}.json primary.n_blocks_compared",
        )
    if matched:
        p = matched["primary"]
        row(
            "blocks the model arm compares (density-matched)",
            p["n_blocks_compared"],
            "real-unknown blocks",
            p["target_blocks"],
            "block carries a scored element and has at least one accepted window carrying one",
            f"{MATCHED_CONTROL}.json primary",
        )
        row(
            "window-elements the per-element estimator uses",
            elem["secondary_per_element"]["on_unmatched_windows"]["names_a_coding_gene"]["pooled"][
                "window_elements"
            ]
            if elem
            else None,
            "elements in the unmatched windows of the 531 compared blocks",
            elem["secondary_per_element"]["on_unmatched_windows"]["names_a_coding_gene"]["pooled"][
                "window_elements"
            ]
            if elem
            else None,
            "every scored element inside every drawn window; 9.18% vs 44.04% is over these, and both "
            "are TARGET-NAMING FREQUENCIES of the deletion model, not accuracies",
            f"{ELEMENT_COUNT_CONTROL}.json secondary_per_element",
        )
    if reach:
        row(
            "block elements entering the reach standardisation",
            reach["reach"]["standardised"]["block_elements_standardised"]
            if "standardised" in reach.get("reach", {})
            else None,
            "scored elements inside real-unknown blocks",
            arm["coverage"]["real_unknown"]["scored_elements"] if arm else None,
            "elements in strata holding at least 30 elements in each arm; excluded strata are reported",
            f"{REACH_CONTROL}.json reach",
        )
    return {
        "rows": [r for r in rows if r["count"] is not None],
        "the_two_shares_that_are_quoted": {
            "6.7%": "59 of 882 tier blocks measured at reciprocal overlap 0.5",
            "11.1%": "59 of the 531 blocks that carry a scored element, same numerator",
            "note": "both are correct; neither may be quoted without its denominator",
        },
        "the_overlap_rule": {
            "0.25": "72 measured blocks",
            "0.5": "59 measured blocks (the committed rule)",
            "0.75": "10 measured blocks",
            "note": "the coverage count is a property of the rule as much as of the data",
        },
    }


# ---- corrected keys for the committed results ------------------------------------------------------

#: The top-level key this lane adds to a committed result. It is ADDITIVE: no existing key or value in
#: any of the four files is changed, and each entry names the string it corrects so a reader who has
#: the old wording in hand can find the correction from it.
CORRECTION_KEY = "corrections_after_the_statistical_review_2026_09_28"

CORRECTIONS: dict[str, dict[str, Any]] = {
    MEASURED_ARM: {
        "coverage_needed": (
            "this field is kept as the record of what the run computed. It is NOT a power result for "
            "the measured experiment: `effect_to_detect` and `model_arm_sd_of_per_block_differences` "
            "are both model output, used to size an experiment whose endpoint is a measurement, and "
            "`blocks_short_of_the_floor` is 20 minus the 1 block compared -- the distance to "
            "`measured.MIN_FOR_A_COMPARISON`, a minimum-reporting rule. That the formula also returned "
            "20 is a coincidence, not a confirmation. The power question is answered in "
            "`clause2_design_power`, from the measured endpoint"
        ),
        "primary_reading.text": (
            "the phrase 'the coverage that would be needed ... for 80% power against an effect the "
            "size of the model arm's own -27.25 points' names the mismatch in its own words: the "
            "effect is the model arm's. Read it as a provisional assumption, not as the size of the "
            "experiment"
        ),
        "coverage.windows_matched_real_unknown.windows_with": (
            "these are WINDOWS satisfying each predicate, not elements. The 223 and the 30 behind the "
            "13.45% rate are windows carrying a CRISPRi-tested element; several windows of one block "
            "can carry the same element, so they are not 223 independent measurements. The "
            "element-level rates over the benchmark's own distinct tested elements are 12.15% "
            "(training, 479 of 3,941) and 14.56% (held-out, 247 of 1,697)"
        ),
        "coverage.real_unknown.blocks_with_no_measured_element": (
            "823 means: no assay measures any of the block's scored elements at reciprocal overlap "
            "0.5. It does not mean no part of those blocks was ever measured. At overlap 0.25 the "
            "measured count is 72 blocks and at 0.75 it is 10, against 59 at the committed rule, so "
            "the count is a property of the rule as much as of the data"
        ),
        "the_two_shares": (
            "59 measured blocks is 6.7% of the 882 tier blocks and 11.1% of the 531 that carry a "
            "scored element. Both are correct and neither may be quoted without its denominator"
        ),
        "sample_size_answer": (
            "see `clause2_design_power.sample_size_range_by_ratio`: a RANGE of compared blocks per "
            "assumed effect, not a number"
        ),
    },
    MATCHED_CONTROL: {
        "secondary_real_minus_neutral": (
            "an interval of -4.55 to +4.30 points is NO DIFFERENCE DETECTED, not equivalence and not "
            "'the two tiers behave identically'. It permits a real difference in either direction of "
            "up to about four and a half points. No equivalence margin can be justified from outside "
            "these data, so no equivalence test is run; see `clause2_design_power.equivalence`"
        ),
        "primary": (
            "`names_a_coding_gene` is the deletion model naming a target. The difference is in "
            "TARGET-NAMING FREQUENCY, not in prediction accuracy, and it shows neither that the model "
            "failed nor that these sequences lack regulatory function"
        ),
    },
    ELEMENT_COUNT_CONTROL: {
        "secondary_real_minus_neutral": (
            "an interval of -5.58 to +4.17 points is NO DIFFERENCE DETECTED, not equivalence. See "
            "`clause2_design_power.equivalence` for why no margin can be justified"
        ),
        "secondary_per_element": (
            "9.18% against 44.04% is the share of scored elements for which the deletion model names a "
            "protein-coding gene, in each arm. It is a difference in target-naming frequency, not in "
            "accuracy: no experimental outcome enters it, so it cannot show that the model failed nor "
            "that the blocks' sequence lacks regulatory function"
        ),
    },
    REACH_CONTROL: {
        "reach_reading.text": (
            "'the failure stays about the sequence' is withdrawn. What the run shows is that the gap "
            "persists after standardising on the tested covariates, reach and class. That the tested "
            "covariates do not explain it is not evidence of what does: unmeasured context, selection "
            "effects, model limitations and remaining geometric differences are all still open. The "
            "defensible sentence is: the gap persists after adjustment for the tested reach and class "
            "variables"
        ),
        "classes.reading": (
            "same correction: class composition not explaining the gap does not identify what does"
        ),
    },
}

#: what every correction above shares, written once into each file it touches
CORRECTION_PREAMBLE = (
    "Added by lane-design on 2026-09-28. Every key and value already in this file is unchanged; this "
    "entry is additive and names the field or the string it corrects. The corrections follow a "
    "statistical review of milestone 1.3 clause 2's notes, whose full text is in docs/ROADMAP.md "
    "beneath milestone 1.3 and in the dated lane-design sections of docs/ATTRIBUTION.md."
)


def annotate_committed_results(results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """Add the corrected keys beside the committed ones. Nothing existing is read back out."""
    done: dict[str, Any] = {}
    for name, entries in CORRECTIONS.items():
        p = results_dir / f"{name}.json"
        if not p.exists():
            done[name] = "absent"
            continue
        payload = json.loads(p.read_text())
        before = {k: v for k, v in payload.items() if k != CORRECTION_KEY}
        payload[CORRECTION_KEY] = {
            "added": "2026-09-28",
            "by": "lane-design",
            "preamble": CORRECTION_PREAMBLE,
            "additive": True,
            "corrects": entries,
        }
        after = {k: v for k, v in payload.items() if k != CORRECTION_KEY}
        if after != before:  # never reachable; the guard is the point
            raise AssertionError(f"{p}: an existing value changed; refusing to write")
        p.write_text(json.dumps(payload, indent=2, default=str))
        done[name] = sorted(entries)
    return done


# ---- the run ---------------------------------------------------------------------------------------


def collect(blocks: int, seed: int, check_experiments: int, check_resamples: int) -> dict[str, Any]:
    t0 = time.time()
    anch = anchors()
    win = window_anchor_from_the_measured_arm()
    window_rate = win.get("rate") or 0.1345
    sw = sweep(anch, window_rate, blocks, seed)
    rng = np.random.default_rng(seed + 1)
    checks = []
    ref = reference(anch, window_rate)
    for ratio in (0.5, 0.25):
        for n in (20, 50, 200):
            cfg = {**ref, "ratio": ratio, "rho_chrom": 0.0}
            m = moments(
                cfg["window_observed_rate"],
                ratio,
                cfg["sensitivity"],
                cfg["false_positive"],
                cfg["elements_per_block"],
                cfg["windows_per_block"],
                cfg["icc"],
                blocks,
                rng,
            )
            checks.append(
                {"ratio": ratio, **bootstrap_check(m, cfg, n, check_experiments, check_resamples, rng)}
            )
    worst = max(c["largest_disagreement"] for c in checks) if checks else None
    # Added after the registered check fired its falsifier, and named as such wherever it is read:
    # the registered cells show the disagreement only at the reporting floor, so the same check is
    # run at more sample sizes to say WHERE the approximation becomes safe. It changes no grid, no
    # threshold and no reading; the registered cells above are drawn first and are unaffected by it.
    extended = []
    if worst is not None and worst > 0.05:
        for ratio in (0.5, 0.25):
            for n in (30, 75, 100):
                cfg = {**ref, "ratio": ratio, "rho_chrom": 0.0}
                m = moments(
                    cfg["window_observed_rate"],
                    ratio,
                    cfg["sensitivity"],
                    cfg["false_positive"],
                    cfg["elements_per_block"],
                    cfg["windows_per_block"],
                    cfg["icc"],
                    blocks,
                    rng,
                )
                extended.append(
                    {"ratio": ratio, **bootstrap_check(m, cfg, n, check_experiments, check_resamples, rng)}
                )
    null_rows = [r for r in sw["rows"] if r["ratio"] == 1.0]
    calibration = {
        "what_this_is": (
            "at ratio 1.0 the two arms have the same true rate, so `power_any_direction` is the "
            "test's false-positive rate and nothing else. It is reported because a power table whose "
            "null is not calibrated is not readable"
        ),
        "false_positive_rate_at_the_reporting_floor": [
            {
                "configuration": r["configuration"],
                "n_blocks": MIN_BLOCKS,
                "rate": r["power_at_the_reporting_floor"]["any_direction"],
            }
            for r in null_rows
        ],
        "false_positive_rate_at_500_blocks": [
            {"configuration": r["configuration"], "rate": r["power_at_500_blocks"]["any_direction"]}
            for r in null_rows
        ],
        "monte_carlo_error_of_the_null_mean_points": [
            {
                "configuration": r["configuration"],
                "se": r["monte_carlo_se_of_the_mean_points"],
                "mean": r["mean_difference_points"],
            }
            for r in null_rows
        ],
    }
    return {
        "result": RESULT,
        "registration": PRE_REGISTRATION,
        "anchors": anch,
        "window_rate_anchor": win,
        "sweep": sw,
        "sample_size_range_by_ratio": ranges(sw["rows"]),
        "bootstrap_check": {
            "cells": checks,
            "largest_disagreement_in_power": worst,
            "threshold": 0.05,
            "normal_approximation_usable": (worst is not None and worst <= 0.05),
            "note": (
                "the check cells set the chromosome ICC to 0, because the design effect is applied to "
                "the normal formula and not simulated; the check is of the approximation, not of the "
                "clustering"
            ),
            "extended_cells_added_after_the_falsifier_fired": extended,
            "extended_cells_are_post_hoc": (
                "these sample sizes were not in the registered check. They were added once the "
                "registered check failed, to locate where the normal approximation becomes safe, and "
                "they change no grid, threshold or reading. Read them as a diagnostic of the "
                "approximation, not as part of the registered analysis"
            ),
            "largest_disagreement_above_the_floor": (
                round(
                    max(c["largest_disagreement"] for c in checks + extended if c["n_blocks"] > MIN_BLOCKS),
                    4,
                )
                if any(c["n_blocks"] > MIN_BLOCKS for c in checks + extended)
                else None
            ),
        },
        "null_calibration": calibration,
        "equivalence": equivalence(),
        "denominators": denominators(),
        "reporting_floor_is_not_a_power_result": {
            "floor_blocks": MIN_BLOCKS,
            "what_it_is": "measured.MIN_FOR_A_COMPARISON, a minimum-reporting rule",
            "what_it_is_not": "a power calculation, an effect size, or a sample size",
            "compared_now": 1,
            "blocks_short_of_the_floor": MIN_BLOCKS - 1,
            "the_withdrawn_sentence": (
                "'19 more blocks would decide it' and 'an experiment whose size is known'. 19 is "
                f"{MIN_BLOCKS} minus the 1 block compared today, the distance to a reporting floor. "
                "The power formula beside it returned 20 from a model effect and a model dispersion, "
                "which made the coincidence look like a result"
            ),
        },
        "power_is_a_probability": (
            "every number here is the probability that an experiment of that size would produce an "
            "interval excluding 0, IF the assumptions of its configuration hold. It is not a "
            "guarantee that such an experiment decides clause 2, and an experiment at 80% power fails "
            "to detect a real effect of the assumed size one time in five"
        ),
        "seconds": round(time.time() - t0, 1),
    }


def manifest(blocks: int, seed: int) -> dict[str, Any]:
    inputs = []
    for name in ms.CRISPRI_FILES:
        p = ms.CRISPRI_KNOWLEDGE / name
        if p.exists():
            inputs.append(mf.input_entry(p, partition=ms.CRISPRI_SPLIT_OF[name]))
    for name in (MEASURED_ARM, MATCHED_CONTROL, ELEMENT_COUNT_CONTROL, REACH_CONTROL):
        p = RESULTS_DIR / f"{name}.json"
        if p.exists():
            inputs.append(mf.input_entry(p, partition=None))
    return {
        "sources": [
            {
                "accession": "EngreitzLab/CRISPR_comparison EPCrisprBenchmark (Gschwind et al. 2025)",
                "version": "main, as fetched; pinned by sha256. Read for its power columns and its "
                "element-level regulation rates only",
            },
            {
                "accession": "this repository, data/results/clause2_*.json",
                "version": "0c8b82d, faeb0da, a7f207f, e8aa571: the committed clause 2 results, read "
                "and never rewritten",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "simulated_blocks_per_configuration": blocks,
            "seed": seed,
            "n_grid": list(N_GRID),
            "power_targets": list(POWER_TARGETS),
            "min_blocks_reporting_floor": MIN_BLOCKS,
            "chromosomes": CHROMOSOMES,
            "block_rate_ratio_grid": PRE_REGISTRATION["grids"]["block_rate_ratio"]["values"],
            "element_icc_grid": PRE_REGISTRATION["grids"]["elements_within_block_icc"]["values"],
            "chromosome_icc_grid": PRE_REGISTRATION["grids"]["blocks_within_chromosome_icc"]["values"],
        },
        "exclusions": [
            "the model arm's -27.25 points is not used as an effect size anywhere",
            "the model arm's per-block dispersion 0.4324 is not used as a dispersion anywhere",
            "no AlphaGenome request is made and no per-element response cache is opened",
            "no new descriptive statistic about the blocks is computed; the committed ones are read",
            "configurations are not fully crossed: one factor at a time from a registered reference",
        ],
        "partitions": {
            "training": "ENCODE benchmark training file, read for power columns and element rates",
            "heldout": "ENCODE benchmark held-out file, read for the same and reported apart",
        },
    }


# ==== item 12 S2 (second external review, 2026-09-28, lane-s2): the committed reference against its anchor


#: what the review reports the committed reference configuration produces, read before anything changed
REVIEW_REPORTED = {"element_detection": 0.1680, "positive_windows": 0.5696}

#: Gauss-Hermite nodes for expectations over a standard normal intercept (probabilists' weights, summing
#: to 1). 96 nodes put the quadrature error far below the Monte Carlo error of any figure it is set beside
_GH_X, _GH_W = np.polynomial.hermite_e.hermegauss(96)
_GH_W = _GH_W / _GH_W.sum()


def _observed_probability(eta: np.ndarray, sensitivity: float, false_positive: float) -> np.ndarray:
    """P(an element is observed positive | its logit): the assay applied to the latent state."""
    p = _expit(eta)
    return p * sensitivity + (1 - p) * false_positive


def quadrature_of_the_committed_reference(
    window_rate: float, sensitivity: float, icc: float, elements: int, false_positive: float = 0.0
) -> dict[str, float]:
    """What `moments` computes, in closed form: the committed model's own expectations.

    `moments` puts the logit-normal intercept's LOCATION at logit(rate / sensitivity), so the mean true
    rate is not rate / sensitivity once the intercept has any spread, and it reads a window as positive if
    ANY of its `elements` tested elements is. Both are computed here without simulation, so the two
    figures the review reports can be traced to their two causes.
    """
    p_loc = min(1.0, window_rate / sensitivity)
    lo = math.log(p_loc / (1 - p_loc))
    sigma = sigma_for_icc(icc)
    p = _expit(lo + sigma * _GH_X)
    q = _observed_probability(lo + sigma * _GH_X, sensitivity, false_positive)
    p_bar, q_bar = float(_GH_W @ p), float(_GH_W @ q)
    return {
        "true_rate_at_the_intercept_location": round(p_loc, 4),
        "mean_true_rate": round(p_bar, 4),
        "observed_element_rate": round(q_bar, 5),
        "observed_positive_unit_rate": round(1 - float(_GH_W @ (1 - q) ** elements), 4),
        "positive_unit_rate_from_any_of_the_elements_alone": round(1 - (1 - window_rate) ** elements, 4),
        "latent_icc_given": icc,
        "observed_scale_icc_it_implies": round(float(_GH_W @ (q - q_bar) ** 2) / (q_bar * (1 - q_bar)), 4),
    }


def reproduce_the_review(
    results_dir: Path = RESULTS_DIR, blocks: int = 200_000, seed: int = 20260928
) -> dict[str, Any]:
    """The review's two figures, from the committed code at its reference configuration, before any change.

    The first draw `sweep` makes is `moments` at the reference configuration and ratio 1.0 on
    `default_rng(seed)`, so the same call on a fresh generator is that draw exactly. The committed row is
    read from the committed result beside it, and the closed form says where each figure comes from.
    """
    anch = anchors()
    win = window_anchor_from_the_measured_arm(results_dir)
    rate = win.get("rate") or 0.1345
    ref = reference(anch, rate)
    m = moments(
        ref["window_observed_rate"],
        1.0,
        ref["sensitivity"],
        ref["false_positive"],
        ref["elements_per_block"],
        ref["windows_per_block"],
        ref["icc"],
        blocks,
        np.random.default_rng(seed),
    )
    committed = {}
    p = results_dir / f"{RESULT}.json"
    if p.exists():
        row = next(
            r
            for r in json.loads(p.read_text())["sweep"]["rows"]
            if r["configuration"] == "reference" and r["ratio"] == 1.0
        )
        committed = {k: row[k] for k in ("realised_element_detection_probability", "window_endpoint_rate")}
    closed = quadrature_of_the_committed_reference(
        rate, ref["sensitivity"], ref["icc"], ref["elements_per_block"], ref["false_positive"]
    )
    tier = json.loads((results_dir / f"{MEASURED_ARM}.json").read_text())["coverage"]["real_unknown"]
    return {
        "reference_configuration": ref,
        "anchor": {"rate": rate, "unit": win.get("unit")},
        "rerun_now": {
            "element_detection": m["realised_element_detection_probability"],
            "positive_windows": m["window_endpoint_rate"],
            "positive_blocks": m["block_endpoint_rate"],
            "seed": seed,
            "simulated_blocks": blocks,
        },
        "committed_row": committed,
        "review_reported": REVIEW_REPORTED,
        "reproduced": (
            round(m["realised_element_detection_probability"], 4) == REVIEW_REPORTED["element_detection"]
            and round(m["window_endpoint_rate"], 4) == REVIEW_REPORTED["positive_windows"]
        ),
        "closed_form": closed,
        "why_the_anchor_is_not_reproduced": [
            "the anchor counts WINDOWS carrying a tested element (30 of 223), and `moments` reads a "
            "window as positive if any of its 6 tested elements is, so the quantity it matched to the "
            "anchor is not the quantity the anchor measured",
            "even at the element level it does not come back: the intercept's location is set at "
            "logit(0.1345 / sensitivity) and a logit-normal intercept with spread raises the MEAN rate "
            "above its location, so the observed element rate is sensitivity x mean rate, not 0.1345",
            "the registered ICC (0.30) is an observed-scale correlation measured on 0/1 calls, and it "
            "was entered into the latent-scale identity; the observed-scale correlation the simulation "
            "then produces is the figure under `closed_form.observed_scale_icc_it_implies`. The "
            "registration promised the realised binary-scale ICC beside the latent one; the committed "
            "result carries only the latent one",
            "the committed test that watches the anchor runs `moments` at 1 element, 1 window and ICC "
            "0, the one configuration where window, element and location rates coincide",
        ],
        "sizes_above_the_tier": {
            "largest_n_searched": max(N_GRID),
            "tier_blocks": tier["blocks"],
            "blocks_carrying_a_scored_element": tier["blocks_carrying_a_scored_element"],
            "n_grid_entries_above_the_tier": [n for n in N_GRID if n > tier["blocks"]],
            "n_grid_entries_above_the_carrying_blocks": [
                n for n in N_GRID if n > tier["blocks_carrying_a_scored_element"]
            ],
        },
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--blocks", type=int, default=200_000, help="simulated blocks per configuration")
    ap.add_argument("--seed", type=int, default=20260928)
    ap.add_argument("--check-experiments", type=int, default=1000)
    ap.add_argument("--check-resamples", type=int, default=2000)
    ap.add_argument("--no-save", action="store_true")
    ap.add_argument(
        "--annotate",
        action="store_true",
        help="add the corrected keys beside the committed clause 2 results (additive, never in place)",
    )
    ap.add_argument(
        "--reproduce-review",
        action="store_true",
        help="item 12 S2: rerun the committed reference configuration and print what it produces",
    )
    args = ap.parse_args(argv)
    if args.reproduce_review:
        print(json.dumps(reproduce_the_review(), indent=1))
        return 0
    if args.annotate:
        for name, what in annotate_committed_results().items():
            print(f"{name}: {what}")
        return 0
    out = collect(args.blocks, args.seed, args.check_experiments, args.check_resamples)
    if args.no_save:
        print("not saved")
    else:
        print(f"saved {save_result(RESULT, out, manifest=manifest(args.blocks, args.seed))}")
    for ratio, r in out["sample_size_range_by_ratio"].items():
        print(
            f"ratio {ratio}: 80% power at {r['smallest']} to {r['largest']} compared blocks "
            f"({r['reference_configuration']} at the reference), "
            f"{len(r['did_not_reach_80_percent_by_n'])} configurations never reach it"
        )
    print("bootstrap check, largest disagreement:", out["bootstrap_check"]["largest_disagreement_in_power"])
    print("equivalence:", out["equivalence"]["admissible_margin_found"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
