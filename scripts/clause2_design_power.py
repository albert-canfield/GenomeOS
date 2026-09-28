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
import importlib.util
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


#: Since the item 12 S6 follow-up (lane-contract) `--annotate` writes the corrections as their own result
#: through save_result instead of editing the four results in place. lane-design's first run (70801de)
#: added CORRECTION_KEY to each committed file; those keys stay exactly as committed and are never
#: rewritten. An in-place edit would sit under a manifest that describes another run, and rewriting a
#: file through save_result would replace its committed code stamp.
RESULT_CORRECTIONS = "clause2_statistical_corrections"


def correction_block(entries: dict[str, Any]) -> dict[str, Any]:
    return {
        "added": "2026-09-28",
        "by": "lane-design",
        "preamble": CORRECTION_PREAMBLE,
        "additive": True,
        "corrects": entries,
    }


def annotate_committed_results(results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """The corrections of the committed results, as the result RESULT_CORRECTIONS (one block per result,
    the block lane-design first added beside each file's own keys). No committed file is written."""
    done: dict[str, Any] = {}
    blocks: dict[str, Any] = {}
    read: list[Path] = []
    for name, entries in CORRECTIONS.items():
        p = results_dir / f"{name}.json"
        if not p.exists():
            done[name] = "absent"
            continue
        read.append(p)
        blocks[name] = correction_block(entries)
        done[name] = sorted(entries)
    payload = {
        "correction_key": CORRECTION_KEY,
        "corrections": blocks,
        "absent": sorted(n for n, v in done.items() if v == "absent"),
        "first_written": f"each of the {len(blocks)} results carries its block under {CORRECTION_KEY}, added "
        "by lane-design (70801de) and kept as committed",
    }
    manifest = {
        "sources": [
            {
                "accession": "statistical review of milestone 1.3 clause 2's notes (docs/ROADMAP.md) and "
                "lane-design's audit (f6cb4d7, 70801de)",
                "version": "2026-09-28",
            }
        ],
        "inputs": [mf.input_entry(p) for p in read],
        "assembly": "n/a: corrections of wording and of field meanings, no genomic data",
        "coordinates": "n/a: no genomic intervals",
        "parameters": {"correction_key": CORRECTION_KEY},
        "exclusions": [f"{n}: absent from the results directory" for n in payload["absent"]],
        "partitions": "n/a: no evaluation split",
    }
    save_result(RESULT_CORRECTIONS, payload, results_dir, manifest=manifest)
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


# ==== item 12 S2: the calibrated model, registered by lane-s2 on 2026-09-28 before its first run =========
#
# Everything above is kept as it was computed. What follows is a second model beside it, built to meet the
# review's acceptance: the endpoint as the committed estimator reads it, the anchor matched at the level it
# was measured and with its clustering, the assay's sensitivity and false positives applied at the element
# level on both sides of the calibration, observed and latent correlation kept apart, the planned analysis
# simulated end to end under the null and each alternative with shared controls and chromosome clustering,
# designs compared in assay units, and the eligible block population as a hard cap.

RESULT_CALIBRATED = "clause2_design_power_calibrated"

#: the designs compared at equal cost: tested elements per block AND per window (k), tested windows per
#: block (m), and blocks sharing one set of tested windows (g, the shared-control factor)
S2_K = (1, 2, 3, 6)
S2_M = (1, 3, 10)
S2_G = (1, 4)
#: lane-design's reference design kept for continuity, and the design every carrying block can enter
S2_REFERENCE_DESIGNS = ((6, 3, 1), (1, 3, 1))
#: compared blocks searched; every design's grid is cut at its eligible population and that population
#: itself is added, so no size above it is ever simulated or printed as a number
S2_N_GRID = (20, 30, 50, 75, 100, 150, 200, 300, 400, 500)
#: total assay budgets, in tested elements, for the equal-cost table
S2_BUDGETS = (250, 500, 1000, 2000, 4000, 8000)
S2_RATIOS = (1.0, 0.75, 0.5, 0.25, 0.1, 0.0)
S2_POWER_TARGETS = (0.5, 0.8, 0.9)
S2_FALSE_POSITIVE = (0.0, 0.01, 0.03)
S2_REFERENCE_FALSE_POSITIVE = 0.01
S2_REFERENCE_SENSITIVITY_COLUMN = "PowerAtEffectSize20"
S2_CHROMOSOME_ICC_VARIANTS = (0.0, 0.03)
S2_EXPERIMENTS = 2000  # simulated experiments per cell
S2_CHROMOSOME_RESAMPLES = 1000  # Monte Carlo resamples of the chromosome interval, per experiment
S2_ANCHOR_BOOTSTRAP = 10_000  # cluster-bootstrap resamples of the anchor
S2_ANCHOR_REPLICATES = 4000  # simulated replicates of the anchor's own windows in the reproduction check
S2_ANCHOR_TOLERANCE = 0.005  # the anchor is reproduced if the simulated window rate is within half a point
S2_ICC_TOLERANCE = 0.02  # and the simulated observed-scale ICC within 0.02 of its target
S2_BOOTSTRAP_CHECK_EXPERIMENTS = 200  # the exact block bootstrap against the committed 10,000-resample one
S2_COMMITTED_RESAMPLES = 10_000  # clause2_matched_control.BOOTSTRAP, the committed interval's resamples
S2_SEED = 2026092812
S2_ICC_BINS_KB = (25, 50, 100, 250, 1000)


S2_REGISTRATION: dict[str, Any] = {
    "registered": "2026-09-28",
    "lane": "lane-s2",
    "item": "docs/ROADMAP.md section 5, item 12, row S2 (second external review)",
    "registered_before": (
        "the anchor's windows were redrawn and before any configuration of this model was simulated. The "
        "constants, the equations, the grids, the checks and their tolerances, the reported quantities and "
        "the words the answer may be stated in are fixed in this dict and in the constants beside it"
    ),
    "what_it_replaces_and_what_it_keeps": (
        "lane-design's committed table (70801de) is kept as the record of what that run computed and is read "
        "as exploratory. Its simulation reproduced neither its anchor's quantity nor its measured "
        "correlation "
        "(b903e1e: 0.16802 element detection and 0.5696 positive windows against a 13.45% window anchor; an "
        "observed-scale ICC of 0.115 against the 0.30 measured). This model is written beside it under a new "
        "result name, clause2_design_power_calibrated"
    ),
    "review_acceptance": [
        "calibrate the actual experimental endpoint",
        "incorporate sensitivity and false positives consistently",
        "distinguish observed from latent correlation",
        "simulate the planned analysis under the null and alternatives, including shared controls and "
        "clustering",
        "compare designs at equal total assay cost",
        "flag any sample size that exceeds the eligible block population",
    ],
    "endpoint": (
        "exactly as the committed estimator reads it (clause2_measured_arm.summarise_measured): per tested "
        "element, the assay calls a significant change in a protein-coding gene (either sign); a BLOCK is "
        "positive if any of its tested elements is called; a WINDOW is positive if any of its tested "
        "elements "
        "is called. Per compared block d = Y - (positive tested windows / tested windows). The endpoint is "
        "therefore a property of a unit and of how many elements it had tested, and every rate below says "
        "which unit and how many elements"
    ),
    "estimator_and_reading": (
        "the committed ones, unchanged: the mean of d over compared blocks; the 95% percentile bootstrap "
        "interval over blocks, and the one over chromosomes beside it; the registered readings of "
        "clause2_measured_arm, which with at least 20 compared blocks read 'wording_wrong' when the interval "
        "over blocks lies wholly below 0 and 'model_failed' otherwise. Probability of detection is the share "
        "of simulated experiments whose interval over blocks lies wholly below 0 when the blocks truly "
        "regulate less; the false-positive rate is the share whose interval excludes 0 either way at the "
        "null. The probability of each registered reading is reported too, because 'model_failed' is read "
        "whenever the interval reaches 0, including when an experiment is simply too small"
    ),
    "latent_model": (
        "per element a TRUE state (regulates a coding gene or not) drawn with probability expit(mu_arm + c + "
        "u): c ~ Normal(0, tau) per chromosome, shared by a block and the windows on its chromosome; u ~ "
        "Normal(0, sigma_unit) per unit (a block, or one tested window), independent between units. Given c "
        "and u, elements are independent. The assay then calls a true regulator with probability s "
        "(sensitivity) and a non-regulator with probability f (false-positive rate), independently per "
        "element. Calls, not states, make the endpoint"
    ),
    "anchor_at_its_level": (
        "the anchor is 30 of 223 WINDOWS carrying a CRISPRi-tested element, and it is matched as that: the "
        "observed rate of positive windows over those windows, each with its own number of tested elements, "
        "must come back at 30/223. To know those numbers the anchor's windows are redrawn with the committed "
        "draw (clause2_matched_control.matched_windows, the committed seed and decile edges, the measured "
        "arm's own verdicts), recording which scored elements each unit holds. Gate, to the digit, before "
        "any use: 44,081 windows drawn, 31,554 carrying a scored element, 223 with a tested element, 30 with "
        "a regulating one, 134 and 27 blocks with such windows, 882 blocks, 531 carrying. If the gate fails "
        "nothing is calibrated and the result says so"
    ),
    "anchor_clustering": (
        "223 windows are not 223 independent units: windows of one block, and windows of different blocks on "
        "one chromosome, can hold the same tested element. The anchor's uncertainty is a cluster bootstrap "
        "(10,000 resamples, ratio of summed positives to summed windows) over (a) connected components of "
        "windows linked by a shared tested element and (b) blocks; the sensitivity grid's ends are the wider "
        "of the two intervals, end by end. The exact binomial interval is reported beside them, labelled as "
        "the as-if-independent one"
    ),
    "sensitivity_and_false_positives": (
        "applied at the ELEMENT level on both sides of the calibration: the same s and f that deconvolve the "
        "anchor are the ones the planned experiment is simulated with (a screen of the benchmark's kind). "
        "s is the mean of a benchmark power column (reference PowerAtEffectSize20, 0.6674; the four others "
        "as one-factor variants). f reference 0.01: if calls are made at a false discovery rate between 5% "
        "and 10% -- the usual range for screens of this kind, not a figure read from these tables -- a call "
        "rate near 13% implies about 0.008 to 0.015 per non-regulating element; 0 and 0.03 are the variants"
    ),
    "correlation_observed_and_latent": (
        "the benchmark ICCs (0.375 at 25 kb to 0.209 at 1 Mb; 0.0142 between chromosomes) are one-way ANOVA "
        "ICCs of 0/1 CALLS: observed-scale figures. They are matched on the observed scale, after the assay "
        "model: the within-unit spread is solved so that Var(q) / (qbar (1 - qbar)) of the calls equals the "
        "measured ICC, and the chromosome spread so that the calls' between-chromosome share equals 0.0142. "
        "The latent ICCs that result are reported beside the observed ones and never substituted for them. "
        "The reference scale is the benchmark bin nearest the median length of the carrying blocks (log "
        "scale); the other bins are one-factor variants"
    ),
    "calibration": (
        "three equations in three unknowns (window location mu, total spread sigma, chromosome part tau), "
        "solved by nested bisection over Gauss-Hermite quadrature: (1) the observed positive-window rate "
        "over the anchor's windows equals the anchor; (2) the observed-scale within-unit ICC equals its "
        "target; (3) the observed-scale chromosome ICC equals its target. A configuration with no solution "
        "is "
        "reported as infeasible, with the largest reachable value, and is not simulated"
    ),
    "the_check_s2_exists_for": (
        "after calibration, SIMULATE (not the quadrature it was solved with) 4,000 replicates of the "
        "anchor's "
        "own 223 windows with their own tested-element counts, and 20,000 clusters of four calls. The anchor "
        "is reproduced if the simulated positive-window rate is within 0.5 points of 30/223 AND the "
        "simulated "
        "observed-scale ICC within 0.02 of its target. If the reference calibration fails either, no size "
        "from this run is read as calibrated and the result says so first. The implied observed ELEMENT rate "
        "is reported beside the benchmark's element-level rates (12.15%, 14.56%) as an out-of-sample "
        "comparison, not a gate"
    ),
    "block_arm": (
        "the block arm's MEAN true per-element rate is ratio x the window arm's, at the same spreads, over "
        "the "
        "registered ratios 1.0, 0.75, 0.5, 0.25, 0.1, 0.0. (lane-design applied the ratio to the location, "
        "which is not the mean once there is spread.) The estimand at each design is computed in closed form "
        "and the interval's coverage of it is reported"
    ),
    "designs": (
        "k tested elements per block AND per tested window (1, 2, 3, 6), so both arms' endpoints count the "
        "same number of calls and the null stays a null; m tested windows per block (1, 3, 10); g blocks of "
        "one chromosome sharing one set of m tested windows (1: none; 4: shared controls). 24 designs. "
        "Reference designs for the one-factor sensitivity analysis: (6, 3, 1), lane-design's reference, and "
        "(1, 3, 1), which every carrying block can enter"
    ),
    "shared_controls_and_clustering_simulated": (
        "blocks are drawn without replacement from the eligible population, so their chromosomes are the "
        "real "
        "ones; the chromosome intercept is drawn once per experiment and shared by blocks and windows on "
        "that "
        "chromosome; with g = 4 the blocks of a chromosome are grouped in fours and each group's blocks read "
        "the same windows, so their differences are correlated exactly as shared controls make them. The "
        "committed interval over blocks treats blocks as independent; the simulation measures what that "
        "costs rather than assuming it"
    ),
    "cost_model": (
        "assay units are tested elements. A design with N compared blocks costs N k for the blocks and "
        "(number of distinct window sets) m k for the windows; nominal per block k (1 + m / g), realised "
        "cost "
        "averaged over the simulated experiments. Designs are compared (a) by the realised cost of the "
        "smallest searched N reaching 80% probability of detection, the cheapest feasible design named per "
        "ratio, and (b) by probability of detection at fixed budgets of 250, 500, 1,000, 2,000, 4,000 and "
        "8,000 tested elements"
    ),
    "eligible_population_is_a_hard_cap": (
        "a k-element design can only use real-unknown blocks holding at least k scored elements; that count, "
        "from the redraw, is its eligible population. N is never simulated above it; the searched grid "
        "(20, 30, 50, 75, 100, 150, 200, 300, 400, 500) is cut at it and the population itself is added. A "
        "power target not reached within it is printed as 'infeasible', never as a number, and so is a "
        "budget "
        "that buys more blocks than it. Every entry of lane-design's committed table above its own design's "
        "eligible population is listed as infeasible beside it. The window arm's supply is reported (share "
        "of carrying windows holding at least k scored elements) and not capped, since windows are drawn "
        "from the rest of the chromosome"
    ),
    "simulation": (
        "2,000 simulated experiments per cell. The interval over blocks is computed EXACTLY: each block's "
        "difference is a multiple of 1/m, so the resampled mean's distribution is the n-fold convolution of "
        "the experiment's own empirical distribution (FFT), the quantity the committed 10,000-resample "
        "percentile interval estimates; a registered check runs the committed resampling interval on 200 "
        "simulated experiments at 8 cells and reports how often the two readings differ. The interval over "
        "chromosomes is the committed resampling one, 1,000 resamples per experiment. Alternatives stop at "
        "the "
        "first searched N reaching 90%; the null runs every N"
    ),
    "null_calibration_rule": (
        "a (design, N) cell is calibrated if its false-positive rate is at most 0.05 + 2 Monte Carlo "
        "standard "
        "errors; otherwise it is reported as anti-conservative, and its sizes are read as optimistic by that "
        "much. Expected before running, stated so it can be scored: shared controls (g = 4) are "
        "anti-conservative, because the committed interval treats blocks that share windows as independent"
    ),
    "sensitivity_analysis": (
        "one factor at a time from the reference calibration, at the two reference designs: the anchor at "
        "its clustered interval's ends; the four other power columns; f = 0 and 0.03; the four other ICC "
        "scales; chromosome ICC 0 and 0.03. Each variant is recalibrated, re-checked against the anchor "
        "(1,000 replicates) and searched to 80%"
    ),
    "assumptions_stated_not_tested": [
        "the anchor's tested elements were chosen by the benchmark's designers; the planned experiment's "
        "elements would be scored elements chosen for the design. The window arm's per-element rate is "
        "carried "
        "across that difference",
        "one sensitivity per element, taken from pair-level power; an element tested against several genes "
        "is treated as one call",
        "matched windows holding k scored elements exist for every block; their supply is reported",
        "blocks and windows share the latent spread; only the mean differs between arms",
    ],
    "exclusions": (
        "the model arm's -27.25 points and its per-block dispersion are not used anywhere, as in "
        "lane-design's "
        "registration, and the test that enforces it covers this code too. 0 AlphaGenome requests: the "
        "element "
        "archive is opened only for which elements are scored"
    ),
    "words": (
        "every size is a probability-of-detection statement under the calibrated assumptions -- 'N compared "
        "blocks detect a difference of this size with probability p' -- never a guarantee that an experiment "
        "decides clause 2"
    ),
}


# ---- the closed forms: expectations over the unit's logit-normal intercept -----------------------------


def unit_positive_probability(
    mu: float | None, sigma: float, sensitivity: float, false_positive: float, k: int
) -> float:
    """P(a unit with k tested elements is observed positive): 1 - E[(1 - q)^k] over its intercept.

    `mu` None means the unit's elements never truly regulate (ratio 0), so only false positives remain.
    """
    if mu is None:
        return 1 - (1 - false_positive) ** k
    q = _observed_probability(mu + sigma * _GH_X, sensitivity, false_positive)
    return float(1 - _GH_W @ (1 - q) ** k)


def marginal_true_rate(mu: float | None, sigma: float) -> float:
    """The mean true per-element regulation rate, E[expit(mu + sigma Z)]: the LATENT rate, not a location."""
    if mu is None:
        return 0.0
    return float(_GH_W @ _expit(mu + sigma * _GH_X))


def observed_icc(mu: float, sigma: float, sensitivity: float, false_positive: float) -> float:
    """The OBSERVED-scale intraclass correlation of 0/1 calls of elements sharing an intercept.

    Given the intercept the calls are independent Bernoulli(q), so the between-unit share of the variance
    of a call is Var(q) / (qbar (1 - qbar)). This is the quantity a one-way ANOVA ICC of 0/1 calls
    estimates, which is what the benchmark's 0.30 is.
    """
    q = _observed_probability(mu + sigma * _GH_X, sensitivity, false_positive)
    qbar = float(_GH_W @ q)
    return float(_GH_W @ (q - qbar) ** 2) / (qbar * (1 - qbar)) if 0 < qbar < 1 else 0.0


def observed_chromosome_icc(
    mu: float, sigma_total: float, tau: float, sensitivity: float, false_positive: float
) -> float:
    """The observed-scale ICC of calls sharing only a chromosome intercept (sd tau) of the total sigma."""
    su = math.sqrt(max(sigma_total**2 - tau**2, 0.0))
    eta = mu + tau * _GH_X[:, None] + su * _GH_X[None, :]
    q = _observed_probability(eta, sensitivity, false_positive)
    per_c = q @ _GH_W
    qbar = float(_GH_W @ per_c)
    return float(_GH_W @ (per_c - qbar) ** 2) / (qbar * (1 - qbar)) if 0 < qbar < 1 else 0.0


def latent_icc(sigma: float) -> float:
    """The latent-scale ICC of an intercept of sd sigma: sigma^2 / (sigma^2 + pi^2 / 3)."""
    return sigma**2 / (sigma**2 + math.pi**2 / 3)


def _bisect(fn: Any, target: float, lo: float, hi: float, iters: int = 100) -> float:
    """Root of an INCREASING fn(x) = target on [lo, hi]."""
    for _ in range(iters):
        mid = (lo + hi) / 2
        if fn(mid) < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def window_rate_over(
    mu: float, sigma: float, sensitivity: float, false_positive: float, k_counts: dict[int, int]
) -> float:
    """The observed positive-window rate over windows whose tested-element counts are `k_counts`."""
    total = sum(k_counts.values())
    return (
        sum(
            c * unit_positive_probability(mu, sigma, sensitivity, false_positive, k)
            for k, c in k_counts.items()
        )
        / total
    )


def mu_for_window_rate(
    rate: float, sigma: float, sensitivity: float, false_positive: float, k_counts: dict[int, int]
) -> float | None:
    """The window arm's logit location that makes the observed window rate `rate` over those windows."""
    floor = window_rate_over(-40.0, sigma, sensitivity, false_positive, k_counts)
    ceil = window_rate_over(40.0, sigma, sensitivity, false_positive, k_counts)
    if not floor < rate < ceil:
        return None
    return _bisect(
        lambda m: window_rate_over(m, sigma, sensitivity, false_positive, k_counts), rate, -40.0, 40.0
    )


def mu_for_marginal_rate(pi: float, sigma: float) -> float | None:
    """The logit location whose MEAN true rate is `pi` at intercept spread sigma; None for pi = 0."""
    if pi <= 0:
        return None
    return _bisect(lambda m: marginal_true_rate(m, sigma), pi, -40.0, 40.0)


def calibrate(
    anchor_rate: float,
    k_counts: dict[int, int],
    sensitivity: float,
    false_positive: float,
    icc_target: float,
    chromosome_icc_target: float,
) -> dict[str, Any]:
    """Solve the latent model so that, AFTER the assay, it reproduces the anchor at the window level and
    the measured correlations on the observed scale.

    Unknowns: the window arm's logit location mu, the total intercept spread sigma (unit plus chromosome)
    and the chromosome part tau. Equations: (1) the observed rate of positive windows over the anchor's own
    windows, each with its own number of tested elements, equals the anchor; (2) the observed-scale ICC of
    calls sharing an intercept equals the measured one; (3) the observed-scale ICC of calls sharing only a
    chromosome equals the measured chromosome ICC. The latent ICCs that result are reported beside the
    observed ones and never substituted for them.
    """
    out: dict[str, Any] = {
        "anchor_rate": anchor_rate,
        "sensitivity": sensitivity,
        "false_positive": false_positive,
        "observed_icc_target": icc_target,
        "observed_chromosome_icc_target": chromosome_icc_target,
        "anchor_tested_elements_per_window": {str(k): c for k, c in sorted(k_counts.items())},
    }
    if sensitivity <= false_positive:
        return {
            **out,
            "feasible": False,
            "why": "an assay no more likely to call a regulator than a non-regulator",
        }

    def icc_at(sigma: float) -> float:
        mu = mu_for_window_rate(anchor_rate, sigma, sensitivity, false_positive, k_counts)
        return -1.0 if mu is None else observed_icc(mu, sigma, sensitivity, false_positive)

    top = 12.0
    if icc_target <= 0:
        sigma = 0.0
    elif icc_at(top) < icc_target:
        return {
            **out,
            "feasible": False,
            "why": (
                f"no intercept spread up to {top} on the logit produces an observed-scale ICC of "
                f"{icc_target} "
                f"at this sensitivity and anchor; the largest reachable is {round(icc_at(top), 4)}"
            ),
        }
    else:
        sigma = _bisect(icc_at, icc_target, 0.0, top)
    mu = mu_for_window_rate(anchor_rate, sigma, sensitivity, false_positive, k_counts)
    if mu is None:
        return {**out, "feasible": False, "why": "the anchor is outside what this assay can produce"}
    tau = 0.0
    if chromosome_icc_target > 0 and sigma > 0:
        if observed_chromosome_icc(mu, sigma, sigma, sensitivity, false_positive) < chromosome_icc_target:
            tau = sigma
        else:
            tau = _bisect(
                lambda t: observed_chromosome_icc(mu, sigma, t, sensitivity, false_positive),
                chromosome_icc_target,
                0.0,
                sigma,
            )
    pi_w = marginal_true_rate(mu, sigma)
    return {
        **out,
        "feasible": True,
        "mu_window": round(mu, 6),
        "sigma_total": round(sigma, 6),
        "tau_chromosome": round(tau, 6),
        "sigma_unit": round(math.sqrt(max(sigma**2 - tau**2, 0.0)), 6),
        "window_true_rate_latent_mean": round(pi_w, 5),
        "window_rate_at_the_location": round(float(_expit(np.array(mu))), 5),
        "observed_element_rate_implied": round(sensitivity * pi_w + false_positive * (1 - pi_w), 5),
        "observed_window_rate_achieved": round(
            window_rate_over(mu, sigma, sensitivity, false_positive, k_counts), 5
        ),
        "observed_icc_achieved": round(observed_icc(mu, sigma, sensitivity, false_positive), 5),
        "observed_chromosome_icc_achieved": round(
            observed_chromosome_icc(mu, sigma, tau, sensitivity, false_positive), 5
        ),
        "latent_icc_unit_plus_chromosome": round(latent_icc(sigma), 5),
        "latent_icc_chromosome": round(tau**2 / (sigma**2 + math.pi**2 / 3), 5),
    }


def block_arm(cal: dict[str, Any], ratio: float) -> dict[str, Any]:
    """The block arm at `ratio`: its MEAN true rate is ratio x the window arm's, at the same spreads."""
    pi_b = ratio * cal["window_true_rate_latent_mean"]
    mu_b = mu_for_marginal_rate(pi_b, cal["sigma_total"])
    return {"ratio": ratio, "block_true_rate_latent_mean": round(pi_b, 5), "mu_block": mu_b}


def expected_difference(cal: dict[str, Any], mu_b: float | None, k: int) -> float:
    """E[Y - W] per compared block for a design testing k elements per unit: the estimand, in points/100."""
    s, f, sg = cal["sensitivity"], cal["false_positive"], cal["sigma_total"]
    return unit_positive_probability(mu_b, sg, s, f, k) - unit_positive_probability(
        cal["mu_window"], sg, s, f, k
    )


# ---- the anchor, redrawn at the level it was measured ---------------------------------------------------


def _measured_arm_module() -> Any:
    spec = importlib.util.spec_from_file_location(
        "clause2_measured_arm", Path(__file__).resolve().parent / "clause2_measured_arm.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class _UnitRecorder:
    """Records, for the committed window draw, which scored elements each block and each window holds.

    `matched_windows` calls `raw_fn(mid)` once when it opens a block and once per accepted window, and it
    evaluates every predicate over a unit's elements. A last predicate that always answers False therefore
    sees every element of every unit in order, and `raw_fn` marks where each unit starts. The draw itself
    consumes no random number on either call, so the windows are the committed windows; the gate checks it.
    """

    def __init__(self, block_mids: set[int]):
        self.block_mids = block_mids
        self.units: list[dict[str, Any]] = []

    def raw(self, mid: int) -> int:
        self.units.append({"kind": "block" if mid in self.block_mids else "window", "mid": mid, "els": []})
        return 0

    def element(self, e: dict[str, Any]) -> bool:
        self.units[-1]["els"].append(e)
        return False


def redraw_the_anchor(results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """The 223 windows behind the 13.45% anchor, redrawn with the committed code and seed, with what each
    holds: its tested elements, whether one regulates a coding gene, and which windows share an element.

    Reads the element archive (for which elements are scored, never for a model answer), the organiser's
    blocks, GENCODE and the measured layer. 0 AlphaGenome requests.
    """
    from genomeos.attribution import organise

    arm = _measured_arm_module()
    mcm, cutm = arm.mc, arm.cut
    committed = json.loads((results_dir / f"{MEASURED_ARM}.json").read_text())
    chroms = committed["chromosomes"]
    denom, primary = arm.DENOMINATOR_QUESTION, arm.PRIMARY_QUESTION
    blocks = {c: organise.blocks(c) for c in chroms}
    tss = {c: mcm.coding_tss(c) for c in chroms}
    coding = arm.coding_symbols(chroms)
    edges = mcm.decile_edges(
        [
            mcm.tss_count(tss[c], (b["start"] + b["end"]) // 2)
            for c in chroms
            for b in mcm.target_sets(blocks[c])["real_unknown"]
        ]
    )
    windows: list[dict[str, Any]] = []
    block_rows: list[dict[str, Any]] = []
    totals = {"windows_drawn": 0, "windows_carrying_a_scored_element": 0}
    for c in chroms:
        els = cutm.elements_of(c)
        if not els:
            continue
        els.sort(key=lambda e: e["start"])
        got = cutm.read_chromosome(c, els)
        if not got:
            continue
        layer = ms.Layer.load(c)
        arm.attach(els, layer, coding)
        targets = mcm.target_sets(got)["real_unknown"]
        rec = _UnitRecorder({(b["start"] + b["end"]) // 2 for b in targets})
        preds = {
            denom: lambda e: bool(e["_measured"][denom]),
            primary: lambda e: bool(e["_measured"][primary]),
            "_record": rec.element,
        }
        t = tss[c]
        runs = mcm.matched_windows(
            targets, got, els, preds, lambda m, t=t: mcm.bin_of(mcm.tss_count(t, m), edges), raw_fn=rec.raw
        )
        by_start = {r["start"]: r for r in runs}
        current = None
        for u in rec.units:
            if u["kind"] == "block":
                current = next(b for b in targets if (b["start"] + b["end"]) // 2 == u["mid"])
                r = by_start[current["start"]]
                block_rows.append(
                    {
                        "chrom": c,
                        "block": f"{c}:{current['start']}-{current['end']}",
                        "length": current["length"],
                        "scored_elements": len(u["els"]),
                        "elements_field": r["elements"],
                        "tested": any(e["_measured"][denom] for e in u["els"]),
                        "windows_drawn": r["drawn"],
                        "windows_carrying": r["carrying"],
                        "windows_tested": r["windows_yes"][denom],
                        "windows_regulating": r["windows_yes"][primary],
                    }
                )
                continue
            totals["windows_drawn"] += 1
            if not u["els"]:
                continue
            totals["windows_carrying_a_scored_element"] += 1
            tested = [e for e in u["els"] if e["_measured"][denom]]
            if not tested:
                windows.append({"chrom": c, "tested": 0, "scored": len(u["els"])})
                continue
            windows.append(
                {
                    "chrom": c,
                    "block": f"{c}:{current['start']}-{current['end']}",
                    "window_mid": u["mid"],
                    "scored": len(u["els"]),
                    "tested": len(tested),
                    "positive": any(e["_measured"][primary] for e in u["els"]),
                    "tested_elements": sorted(f"{c}:{e['start']}-{e['end']}" for e in tested),
                    "regulating_elements": sorted(
                        f"{c}:{e['start']}-{e['end']}" for e in tested if e["_measured"][primary]
                    ),
                }
            )
        del els, got, layer, rec
        print(f"{c}: redrawn", flush=True)
    anchor = [w for w in windows if w["tested"]]
    want = committed["coverage"]["windows_matched_real_unknown"]
    got_counts = {
        "windows_drawn": totals["windows_drawn"],
        "windows_carrying_a_scored_element": totals["windows_carrying_a_scored_element"],
        "windows_with_a_tested_element": len(anchor),
        "windows_with_a_regulating_element": sum(w["positive"] for w in anchor),
        "blocks_with_a_tested_window": len({w["block"] for w in anchor}),
        "blocks_with_a_regulating_window": len({w["block"] for w in anchor if w["positive"]}),
        "real_unknown_blocks": len(block_rows),
        "blocks_carrying_a_scored_element": sum(1 for b in block_rows if b["scored_elements"]),
        "per_block_counts_agree_with_the_draw": all(
            b["scored_elements"] == b["elements_field"] for b in block_rows
        )
        and sum(b["windows_tested"] for b in block_rows) == len(anchor),
    }
    want_counts = {
        "windows_drawn": want["windows_drawn"],
        "windows_carrying_a_scored_element": want["windows_carrying_a_scored_element"],
        "windows_with_a_tested_element": want["windows_with"][denom],
        "windows_with_a_regulating_element": want["windows_with"][primary],
        "blocks_with_a_tested_window": want["blocks_with_at_least_one_such_window"][denom],
        "blocks_with_a_regulating_window": want["blocks_with_at_least_one_such_window"][primary],
        "real_unknown_blocks": committed["coverage"]["real_unknown"]["blocks"],
        "blocks_carrying_a_scored_element": committed["coverage"]["real_unknown"][
            "blocks_carrying_a_scored_element"
        ],
        "per_block_counts_agree_with_the_draw": True,
    }
    return {
        "chromosomes": chroms,
        "gate": {"want": want_counts, "got": got_counts, "passed": want_counts == got_counts},
        "anchor_windows": anchor,
        "blocks": block_rows,
        "window_scored_elements": [w["scored"] for w in windows],
    }


def anchor_structure(anchor: list[dict[str, Any]], seed: int = S2_SEED) -> dict[str, Any]:
    """The anchor at the level it was measured: its windows' tested-element counts, the distinct elements
    beneath them, and its uncertainty with the windows clustered (shared elements; blocks) rather than
    counted as independent."""
    n = len(anchor)
    pos = sum(w["positive"] for w in anchor)
    k_counts: dict[int, int] = defaultdict(int)
    for w in anchor:
        k_counts[w["tested"]] += 1
    # components: windows linked when they share a tested element
    parent = list(range(n))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    first: dict[str, int] = {}
    for i, w in enumerate(anchor):
        for e in w["tested_elements"]:
            if e in first:
                a, b = find(i), find(first[e])
                if a != b:
                    parent[a] = b
            else:
                first[e] = i
    comp = [find(i) for i in range(n)]
    distinct = {e for w in anchor for e in w["tested_elements"]}
    regulating = {e for w in anchor for e in w["regulating_elements"]}

    def cluster_interval(keys: list[Any]) -> dict[str, Any]:
        groups: dict[Any, list[int]] = defaultdict(list)
        for key, w in zip(keys, anchor, strict=True):
            groups[key].append(int(w["positive"]))
        ys = np.array([sum(g) for g in groups.values()], dtype=float)
        ns = np.array([len(g) for g in groups.values()], dtype=float)
        rng = np.random.default_rng(seed)
        idx = rng.integers(0, len(ys), size=(S2_ANCHOR_BOOTSTRAP, len(ys)))
        boots = ys[idx].sum(axis=1) / ns[idx].sum(axis=1)
        lo, hi = np.percentile(boots, [2.5, 97.5])
        # the design effect of the clustering on the variance of the rate
        p = ys.sum() / ns.sum()
        var_iid = p * (1 - p) / ns.sum()
        return {
            "clusters": len(ys),
            "largest_cluster_windows": int(ns.max()),
            "interval_95": [round(float(lo), 4), round(float(hi), 4)],
            "design_effect": round(float(boots.var()) / var_iid, 3) if var_iid > 0 else None,
            "effective_independent_windows": round(n / (float(boots.var()) / var_iid), 1)
            if var_iid > 0
            else None,
        }

    by_component = cluster_interval(comp)
    by_block = cluster_interval([w["block"] for w in anchor])
    lo = min(by_component["interval_95"][0], by_block["interval_95"][0])
    hi = max(by_component["interval_95"][1], by_block["interval_95"][1])
    return {
        "windows": n,
        "positive_windows": pos,
        "rate": round(pos / n, 4) if n else None,
        "unit": "matched windows carrying at least one CRISPRi-tested element (tested against a coding gene)",
        "tested_elements_per_window": {str(k): c for k, c in sorted(k_counts.items())},
        "distinct_tested_elements": len(distinct),
        "distinct_tested_elements_regulating": len(regulating),
        "element_level_rate_beneath_the_windows": round(len(regulating) / len(distinct), 4)
        if distinct
        else None,
        "positive_windows_whose_positive_is_an_untested_element": sum(
            1 for w in anchor if w["positive"] and not w["regulating_elements"]
        ),
        "windows_sharing_a_tested_element_with_another": sum(
            1 for i in range(n) if sum(1 for j in comp if j == comp[i]) > 1
        ),
        "cluster_bootstrap_by_shared_element": by_component,
        "cluster_bootstrap_by_block": by_block,
        "exact_binomial_interval_as_if_independent": [round(x, 4) for x in _exact_binomial_interval(pos, n)],
        "grid_ends_used": [round(lo, 4), round(hi, 4)],
        "grid_rule": "the wider of the two clustered intervals, end by end",
    }


# ---- the planned experiment, simulated end to end -------------------------------------------------------


def exact_block_bootstrap_tails(
    y: np.ndarray, wbar: np.ndarray, m: int, delta: float = 0.0
) -> dict[str, np.ndarray]:
    """The committed percentile bootstrap over blocks, computed exactly rather than by resampling.

    Each block's difference y - wbar is a multiple of 1/m, so the resampled mean is a lattice sum of n iid
    draws from the experiment's own empirical distribution; its exact distribution is the n-fold
    convolution of that distribution (by FFT). The committed 10,000-resample percentile interval estimates
    exactly these quantiles. With the inverted-CDF percentile, the interval lies wholly below `delta` iff
    P*(mean < delta) >= 0.975 and wholly above it iff P*(mean <= delta) < 0.025.
    """
    e, n = y.shape
    code = (y.astype(np.int64) * m - np.rint(wbar * m).astype(np.int64)) + m  # 0 .. 2m
    u = 2 * m + 1
    pmf = np.stack([(code == j).sum(axis=1) for j in range(u)], axis=1) / n  # (e, u)
    size = 1 << int(math.ceil(math.log2(2 * m * n + 1)))
    phi = np.fft.rfft(pmf, size, axis=1)
    dist = np.fft.irfft(phi**n, size, axis=1)[:, : 2 * m * n + 1]
    dist = np.clip(dist, 0.0, None)
    dist /= dist.sum(axis=1, keepdims=True)
    # lattice point j is the mean (j - m n) / (m n)
    cut = delta * m * n + m * n
    js = np.arange(2 * m * n + 1)
    below = dist[:, js < cut - 1e-9].sum(axis=1)
    at_or_below = dist[:, js <= cut + 1e-9].sum(axis=1)
    return {"p_below": below, "p_at_or_below": at_or_below}


def _chromosome_interval(
    d: np.ndarray, chrom: np.ndarray, resamples: int, rng: np.random.Generator
) -> tuple[np.ndarray, np.ndarray]:
    """The committed interval over chromosomes, per experiment: resample the chromosomes that hold a
    compared block, with replacement, and take the ratio of summed differences to summed blocks."""
    e, _ = d.shape
    sums = np.zeros((e, CHROMOSOMES))
    ns = np.zeros((e, CHROMOSOMES))
    rows = np.repeat(np.arange(e), d.shape[1])
    np.add.at(sums, (rows, chrom.ravel()), d.ravel())
    np.add.at(ns, (rows, chrom.ravel()), 1)
    present = ns > 0
    order = np.argsort(~present, axis=1, kind="stable")  # present chromosomes first
    counts = present.sum(axis=1)
    s_sorted = np.take_along_axis(sums, order, axis=1)
    n_sorted = np.take_along_axis(ns, order, axis=1)
    lo = np.empty(e)
    hi = np.empty(e)
    for i in range(e):
        held = counts[i]
        idx = rng.integers(0, held, size=(resamples, held))
        b = s_sorted[i, :held][idx].sum(axis=1) / n_sorted[i, :held][idx].sum(axis=1)
        lo[i], hi[i] = np.percentile(b, [2.5, 97.5])
    return lo, hi


def simulate_design(
    cal: dict[str, Any],
    mu_b: float | None,
    design: tuple[int, int, int],
    n: int,
    eligible_chroms: np.ndarray,
    experiments: int,
    rng: np.random.Generator,
    delta: float,
    chromosome_resamples: int = S2_CHROMOSOME_RESAMPLES,
    chunk: int = 250,
) -> dict[str, Any]:
    """`experiments` runs of the planned experiment and its committed analysis at `n` compared blocks.

    Per experiment: n blocks drawn without replacement from the eligible population (so their chromosomes
    are the real ones); one latent intercept per chromosome shared by the blocks and the windows on it; one
    per block and one per window set; k tested elements per unit, each called by the assay with the
    registered sensitivity and false-positive rate; a unit is positive if any of its calls is. Windows are
    shared by g blocks of the same chromosome. The committed estimator and both committed intervals are then
    computed on the simulated data.
    """
    k, m, g = design
    if n > len(eligible_chroms):
        raise ValueError(f"{n} compared blocks exceeds the eligible population of {len(eligible_chroms)}")
    s, f = cal["sensitivity"], cal["false_positive"]
    tau, su, mu_w = cal["tau_chromosome"], cal["sigma_unit"], cal["mu_window"]

    def positive(mu: float | None, shift: np.ndarray) -> np.ndarray:
        q = np.full(shift.shape, f) if mu is None else _observed_probability(mu + shift, s, f)
        return rng.random(shift.shape) < 1 - (1 - q) ** k

    stats = defaultdict(list)
    for start in range(0, experiments, chunk):
        e = min(chunk, experiments - start)
        pick = np.argsort(rng.random((e, len(eligible_chroms))), axis=1)[:, :n]
        chrom = eligible_chroms[pick]
        c = rng.normal(0.0, tau, (e, CHROMOSOMES))
        cb = np.take_along_axis(c, chrom, axis=1)
        y = positive(mu_b, cb + rng.normal(0.0, su, (e, n)))
        if g == 1:
            w = positive(mu_w, cb[:, :, None] + rng.normal(0.0, su, (e, n, m)))
            wbar = w.mean(axis=2)
            window_sets = np.full(e, n)
        else:
            order = np.argsort(chrom, axis=1, kind="stable")
            sc = np.take_along_axis(chrom, order, axis=1)
            posn = np.broadcast_to(np.arange(n), (e, n))
            starts = np.ones((e, n), dtype=bool)
            starts[:, 1:] = sc[:, 1:] != sc[:, :-1]
            rank = posn - np.maximum.accumulate(np.where(starts, posn, 0), axis=1)
            per = n // g + 1
            gid_sorted = sc * per + rank // g
            gid = np.empty_like(gid_sorted)
            np.put_along_axis(gid, order, gid_sorted, axis=1)
            grp_chrom = np.arange(CHROMOSOMES * per) // per
            cg = c[:, grp_chrom]
            w = positive(mu_w, cg[:, :, None] + rng.normal(0.0, su, cg.shape + (m,)))
            wbar = np.take_along_axis(w.mean(axis=2), gid, axis=1)
            window_sets = np.array([len(np.unique(row)) for row in gid])
        d = y - wbar
        tails0 = exact_block_bootstrap_tails(y, wbar, m, 0.0)
        tails_t = exact_block_bootstrap_tails(y, wbar, m, delta)
        clo, chi = _chromosome_interval(d, chrom, chromosome_resamples, rng)
        stats["mean"].append(d.mean(axis=1))
        stats["below"].append(tails0["p_below"] >= 0.975)  # interval wholly below 0
        stats["above"].append(tails0["p_at_or_below"] < 0.025)  # wholly above 0
        stats["covers"].append((tails_t["p_below"] < 0.975) & (tails_t["p_at_or_below"] >= 0.025))
        stats["chrom_below"].append(chi < 0)
        stats["chrom_above"].append(clo > 0)
        stats["cost"].append(n * k + window_sets * m * k)
    a = {key: np.concatenate(v) for key, v in stats.items()}
    below, above = a["below"].mean(), a["above"].mean()

    def se(p: float) -> float:
        return round(math.sqrt(max(p * (1 - p), 1e-12) / experiments), 4)

    return {
        "n_blocks": n,
        "experiments": experiments,
        "expected_difference_points": round(100 * delta, 3),
        "mean_estimate_points": round(100 * float(a["mean"].mean()), 3),
        "sd_of_the_estimate_points": round(100 * float(a["mean"].std(ddof=1)), 3),
        "interval_below_zero": round(float(below), 4),
        "interval_above_zero": round(float(above), 4),
        "interval_excludes_zero": round(float(below + above), 4),
        "monte_carlo_se": se(float(below + above)),
        "coverage_of_the_expected_difference": round(float(a["covers"].mean()), 4),
        "registered_reading_wording_wrong": round(float(below), 4),
        "registered_reading_model_failed": round(1 - float(below), 4),
        "chromosome_interval_below_zero": round(float(a["chrom_below"].mean()), 4),
        "chromosome_interval_excludes_zero": round(
            float(a["chrom_below"].mean() + a["chrom_above"].mean()), 4
        ),
        "assay_cost_elements_mean": round(float(a["cost"].mean()), 1),
    }


def design_cost_per_block(design: tuple[int, int, int]) -> float:
    """Nominal tested elements per compared block: its own k, and its share of m windows of k each."""
    k, m, g = design
    return k * (1 + m / g)


def eligible_population(blocks: list[dict[str, Any]], k: int) -> int:
    """Real-unknown blocks holding at least k scored elements: the most blocks a k-element design can use."""
    return sum(1 for b in blocks if b["scored_elements"] >= k)


def n_grid_for(cap: int) -> list[int]:
    """The searched sizes for a design: the registered grid cut at its cap, and the cap itself."""
    grid = [n for n in S2_N_GRID if n <= cap and n >= MIN_BLOCKS]
    if cap >= MIN_BLOCKS and cap not in grid:
        grid.append(cap)
    return grid


def bootstrap_agreement_check(
    cal: dict[str, Any],
    mu_b: float | None,
    design: tuple[int, int, int],
    n: int,
    eligible_chroms: np.ndarray,
    experiments: int,
    rng: np.random.Generator,
) -> dict[str, Any]:
    """The exact block bootstrap against the committed one (10,000 resamples, percentile, as
    `summarise_measured` runs it), on the same simulated experiments: the share of experiments where the
    two read the interval differently."""
    k, m, g = design
    s, f = cal["sensitivity"], cal["false_positive"]
    tau, su, mu_w = cal["tau_chromosome"], cal["sigma_unit"], cal["mu_window"]
    disagree = 0
    for _ in range(experiments):
        pick = rng.permutation(len(eligible_chroms))[:n]
        chrom = eligible_chroms[pick]
        c = rng.normal(0.0, tau, CHROMOSOMES)
        shift_b = c[chrom] + rng.normal(0.0, su, n)
        qb = np.full(n, f) if mu_b is None else _observed_probability(mu_b + shift_b, s, f)
        y = rng.random(n) < 1 - (1 - qb) ** k
        shift_w = c[chrom][:, None] + rng.normal(0.0, su, (n, m))
        qw = _observed_probability(mu_w + shift_w, s, f)
        wbar = (rng.random((n, m)) < 1 - (1 - qw) ** k).mean(axis=1)
        d = y.astype(float) - wbar
        idx = rng.integers(0, n, size=(S2_COMMITTED_RESAMPLES, n))
        lo, hi = np.percentile(d[idx].mean(axis=1), [2.5, 97.5])
        committed = (hi < 0, lo > 0)
        t = exact_block_bootstrap_tails(y[None, :], wbar[None, :], m, 0.0)
        exact = (bool(t["p_below"][0] >= 0.975), bool(t["p_at_or_below"][0] < 0.025))
        disagree += committed != exact
    return {
        "design": list(design),
        "n_blocks": n,
        "experiments": experiments,
        "committed_resamples": S2_COMMITTED_RESAMPLES,
        "readings_that_differ": disagree,
        "share_that_differ": round(disagree / experiments, 4),
    }


# ---- the run -------------------------------------------------------------------------------------------


def _anova_icc_equal(calls: np.ndarray) -> float:
    """One-way ANOVA ICC for equal cluster sizes, rows = clusters: the estimator the benchmark ICC used."""
    k, n = calls.shape[0], calls.shape[1]
    means = calls.mean(axis=1)
    grand = calls.mean()
    msb = n * ((means - grand) ** 2).sum() / (k - 1)
    msw = ((calls - means[:, None]) ** 2).sum() / (k * (n - 1))
    return float((msb - msw) / (msb + (n - 1) * msw))


def anchor_reproduction_check(
    cal: dict[str, Any],
    k_counts: dict[int, int],
    rng: np.random.Generator,
    replicates: int = S2_ANCHOR_REPLICATES,
) -> dict[str, Any]:
    """The check S2 exists for, by simulation rather than by the closed form the model was solved with:
    simulate the anchor's own windows, each with its own number of tested elements, and compare the
    observed positive-window rate with the anchor; simulate clusters of calls and compare the observed-
    scale ICC with its target."""
    s, f, sg, mu = cal["sensitivity"], cal["false_positive"], cal["sigma_total"], cal["mu_window"]
    ks = np.repeat(np.array(sorted(k_counts)), [k_counts[k] for k in sorted(k_counts)])
    q = _observed_probability(mu + rng.normal(0.0, sg, (replicates, len(ks))), s, f)
    rate = float((rng.random(q.shape) < 1 - (1 - q) ** ks).mean())
    clusters = 20_000
    qq = _observed_probability(mu + rng.normal(0.0, sg, clusters), s, f)
    calls = (rng.random((clusters, 4)) < qq[:, None]).astype(float)
    icc = _anova_icc_equal(calls)
    diff = rate - cal["anchor_rate"]
    return {
        "anchor_rate": cal["anchor_rate"],
        "simulated_window_rate": round(rate, 5),
        "difference_points": round(100 * diff, 3),
        "tolerance_points": 100 * S2_ANCHOR_TOLERANCE,
        "window_rate_reproduced": abs(diff) <= S2_ANCHOR_TOLERANCE,
        "observed_icc_target": cal["observed_icc_target"],
        "simulated_observed_icc": round(icc, 4),
        "icc_tolerance": S2_ICC_TOLERANCE,
        "observed_icc_reproduced": abs(icc - cal["observed_icc_target"]) <= S2_ICC_TOLERANCE,
        "simulated_element_rate": round(float(calls.mean()), 5),
        "reproduced": abs(diff) <= S2_ANCHOR_TOLERANCE
        and abs(icc - cal["observed_icc_target"]) <= S2_ICC_TOLERANCE,
        "replicates_of_the_anchor_windows": replicates,
        "simulated_clusters_of_four": clusters,
    }


def _search(
    cal: dict[str, Any],
    design: tuple[int, int, int],
    ratio: float,
    eligible_chroms: np.ndarray,
    experiments: int,
    rng: np.random.Generator,
    stop_at: float | None,
) -> dict[str, Any]:
    """Every searched size of one design at one ratio, up to its cap; alternatives stop once the highest
    power target is reached, the null runs the whole grid."""
    k = design[0]
    arm = block_arm(cal, ratio)
    delta = expected_difference(cal, arm["mu_block"], k)
    cap = len(eligible_chroms)
    cells = []
    for n in n_grid_for(cap):
        cell = simulate_design(cal, arm["mu_block"], design, n, eligible_chroms, experiments, rng, delta)
        cells.append(cell)
        if stop_at is not None and cell["interval_below_zero"] >= stop_at:
            break
    sizes: dict[str, Any] = {}
    for t in S2_POWER_TARGETS:
        hit = next((c for c in cells if c["interval_below_zero"] >= t), None)
        if hit is not None:
            sizes[str(t)] = {
                "n_blocks": hit["n_blocks"],
                "assay_cost_elements": hit["assay_cost_elements_mean"],
                "power": hit["interval_below_zero"],
            }
        else:
            at_cap = cells[-1] if cells else None
            sizes[str(t)] = {
                "n_blocks": "infeasible",
                "why": (
                    f"{int(100 * t)}% probability of detection is not reached within the eligible "
                    "population of "
                    f"{cap} blocks"
                    + (f"; at all {cap} of them it is {at_cap['interval_below_zero']}" if at_cap else "")
                ),
            }
    return {
        "design": {"k": design[0], "m": design[1], "g": design[2]},
        "ratio": ratio,
        "block_true_rate_latent_mean": arm["block_true_rate_latent_mean"],
        "expected_difference_points": round(100 * delta, 3),
        "eligible_population": cap,
        "cells": cells,
        "smallest_n_for_power": sizes,
    }


def collect_calibrated(experiments: int = S2_EXPERIMENTS, seed: int = S2_SEED) -> dict[str, Any]:
    from statistics import median

    t0 = time.time()
    rng = np.random.default_rng(seed)
    out: dict[str, Any] = {
        "result": RESULT_CALIBRATED,
        "registration": S2_REGISTRATION,
        "reproduction_of_the_review": reproduce_the_review(),
    }
    anch = anchors()
    out["benchmark_anchors"] = {
        "power_column_means": {c: v["mean"] for c, v in anch["power_columns"].items()},
        "element_level_regulation_rate": anch["element_level_regulation_rate"],
        "elements_within_cluster_icc_observed_scale": anch["elements_within_cluster_icc"],
        "blocks_within_chromosome_icc_observed_scale": anch["blocks_within_chromosome_icc"],
    }
    red = redraw_the_anchor()
    out["anchor_redraw_gate"] = red["gate"]
    if not red["gate"]["passed"]:
        out["reading"] = "the redraw did not reproduce the anchor's windows: nothing is calibrated from it"
        out["seconds"] = round(time.time() - t0, 1)
        return out
    struct = anchor_structure(red["anchor_windows"], seed)
    out["anchor"] = struct
    out["anchor_windows"] = red["anchor_windows"]
    k_counts = {int(k): c for k, c in struct["tested_elements_per_window"].items()}
    chroms = red["chromosomes"]
    cidx = {c: i for i, c in enumerate(chroms)}
    eligible = {
        k: np.array([cidx[b["chrom"]] for b in red["blocks"] if b["scored_elements"] >= k], dtype=np.int64)
        for k in S2_K
    }
    carrying_windows = [x for x in red["window_scored_elements"] if x]
    out["eligible_population"] = {
        "tier_blocks": len(red["blocks"]),
        "by_elements_tested_per_block": {
            str(k): {
                "blocks": int(len(eligible[k])),
                "rule": f"real-unknown blocks holding at least {k} scored element(s)",
                "carrying_windows_holding_at_least_k": round(
                    sum(1 for x in carrying_windows if x >= k) / len(carrying_windows), 4
                ),
            }
            for k in S2_K
        },
    }
    lengths = [b["length"] for b in red["blocks"] if b["scored_elements"]]
    med = float(median(lengths))
    icc_bin = min(S2_ICC_BINS_KB, key=lambda kb: abs(math.log(kb * 1000 / med)))
    out["icc_scale"] = {
        "median_length_of_carrying_blocks_bp": med,
        "bin_used_kb": icc_bin,
        "rule": "the benchmark ICC bin nearest the median carrying-block length on a log scale",
    }
    icc_target = anch["elements_within_cluster_icc"][f"{icc_bin}kb"]["icc"]
    chrom_icc = anch["blocks_within_chromosome_icc"]["icc"]
    s_ref = anch["power_columns"][S2_REFERENCE_SENSITIVITY_COLUMN]["mean"]
    ref_settings = {
        "anchor_rate": struct["rate"],
        "k_counts": k_counts,
        "sensitivity": s_ref,
        "false_positive": S2_REFERENCE_FALSE_POSITIVE,
        "icc_target": icc_target,
        "chromosome_icc_target": chrom_icc,
    }
    cal = calibrate(**ref_settings)
    out["calibration_reference"] = cal
    if not cal["feasible"]:
        out["reading"] = "the reference calibration is infeasible: " + cal["why"]
        out["seconds"] = round(time.time() - t0, 1)
        return out
    check = anchor_reproduction_check(cal, k_counts, rng)
    out["anchor_reproduced"] = check
    out["the_committed_reference_under_the_same_check"] = {
        "positive_windows": out["reproduction_of_the_review"]["rerun_now"]["positive_windows"],
        "observed_scale_icc": out["reproduction_of_the_review"]["closed_form"][
            "observed_scale_icc_it_implies"
        ],
        "anchor": struct["rate"],
        "observed_icc_measured": out["reproduction_of_the_review"]["reference_configuration"]["icc"],
    }
    designs = [(k, m, g) for k in S2_K for m in S2_M for g in S2_G]
    top = max(S2_POWER_TARGETS)
    by_design = []
    for d in designs:
        for ratio in S2_RATIOS:
            by_design.append(
                _search(cal, d, ratio, eligible[d[0]], experiments, rng, None if ratio == 1.0 else top)
            )
        print(f"design {d}: done ({time.time() - t0:.0f} s)", flush=True)
    out["by_design"] = by_design
    null_cells = [(r["design"], c) for r in by_design if r["ratio"] == 1.0 for c in r["cells"]]
    thresh = [0.05 + 2 * c["monte_carlo_se"] for _, c in null_cells]
    out["null_calibration"] = {
        "what_this_is": (
            "at ratio 1.0 blocks and windows regulate equally often, so the share of simulated experiments "
            "whose committed interval excludes 0 is the analysis's false-positive rate; nominal 0.05"
        ),
        "cells": [
            {
                **dd,
                "n_blocks": c["n_blocks"],
                "false_positive_rate": c["interval_excludes_zero"],
                "monte_carlo_se": c["monte_carlo_se"],
                "calibrated": c["interval_excludes_zero"] <= t,
                "chromosome_interval_false_positive_rate": c["chromosome_interval_excludes_zero"],
                "false_wording_wrong_reading": c["registered_reading_wording_wrong"],
            }
            for (dd, c), t in zip(null_cells, thresh, strict=True)
        ],
    }
    cells = out["null_calibration"]["cells"]
    out["null_calibration"]["summary"] = {
        "cells": len(cells),
        "calibrated": sum(c["calibrated"] for c in cells),
        "largest_false_positive_rate": max(c["false_positive_rate"] for c in cells),
        "largest_without_shared_controls": max(c["false_positive_rate"] for c in cells if c["g"] == 1),
        "largest_with_shared_controls": max(c["false_positive_rate"] for c in cells if c["g"] > 1),
        "largest_chromosome_interval_rate": max(c["chromosome_interval_false_positive_rate"] for c in cells),
    }
    # the sample-size table at equal cost: every design's 80% size with its cost, the cheapest per ratio
    table = []
    for ratio in S2_RATIOS[1:]:
        rows = []
        for r in (x for x in by_design if x["ratio"] == ratio):
            got = r["smallest_n_for_power"]["0.8"]
            rows.append({**r["design"], "eligible_population": r["eligible_population"], **got})
        feasible = [x for x in rows if x["n_blocks"] != "infeasible"]
        best = min(feasible, key=lambda x: (x["assay_cost_elements"], x["n_blocks"])) if feasible else None
        table.append(
            {
                "ratio": ratio,
                "designs": rows,
                "infeasible_designs": len(rows) - len(feasible),
                "cheapest_feasible_design_at_80_percent": best,
            }
        )
    out["sample_size_at_equal_cost"] = table
    # power at fixed budgets
    budget_rows = []
    for budget in S2_BUDGETS:
        for d in designs:
            k = d[0]
            n = int(budget // design_cost_per_block(d))
            cap = len(eligible[k])
            row: dict[str, Any] = {"budget_elements": budget, "k": d[0], "m": d[1], "g": d[2], "n_blocks": n}
            if n < MIN_BLOCKS:
                row["status"] = f"below the reporting floor: this budget buys {n} compared blocks"
            elif n > cap:
                row["n_blocks"] = "infeasible"
                row["status"] = (
                    f"infeasible: this budget buys {n} compared blocks and only {cap} hold {k} scored "
                    "element(s)"
                )
            else:
                row["status"] = "simulated"
                row["by_ratio"] = {}
                for ratio in S2_RATIOS:
                    arm = block_arm(cal, ratio)
                    delta = expected_difference(cal, arm["mu_block"], k)
                    c = simulate_design(cal, arm["mu_block"], d, n, eligible[k], experiments, rng, delta)
                    row["by_ratio"][str(ratio)] = {
                        "probability_of_detection": c["interval_below_zero"],
                        "false_positive_rate" if ratio == 1.0 else "interval_excludes_zero": c[
                            "interval_excludes_zero"
                        ],
                        "assay_cost_elements": c["assay_cost_elements_mean"],
                    }
            budget_rows.append(row)
        print(f"budget {budget}: done ({time.time() - t0:.0f} s)", flush=True)
    out["power_at_fixed_budgets"] = budget_rows
    # one factor at a time from the reference calibration, at the two reference designs
    lo, hi = struct["grid_ends_used"]
    variants: list[tuple[str, dict[str, Any]]] = [
        (f"anchor={lo}", {**ref_settings, "anchor_rate": lo}),
        (f"anchor={hi}", {**ref_settings, "anchor_rate": hi}),
    ]
    for col in ms.POWER_COLUMNS:
        if col != S2_REFERENCE_SENSITIVITY_COLUMN:
            variants.append(
                (f"sensitivity={col}", {**ref_settings, "sensitivity": anch["power_columns"][col]["mean"]})
            )
    for fp in S2_FALSE_POSITIVE:
        if fp != S2_REFERENCE_FALSE_POSITIVE:
            variants.append((f"false_positive={fp}", {**ref_settings, "false_positive": fp}))
    for kb in S2_ICC_BINS_KB:
        if kb != icc_bin:
            variants.append(
                (
                    f"icc_scale={kb}kb",
                    {**ref_settings, "icc_target": anch["elements_within_cluster_icc"][f"{kb}kb"]["icc"]},
                )
            )
    for rho in S2_CHROMOSOME_ICC_VARIANTS:
        variants.append((f"chromosome_icc={rho}", {**ref_settings, "chromosome_icc_target": rho}))
    sens = []
    for name, settings in variants:
        vcal = calibrate(**settings)
        entry: dict[str, Any] = {"variant": name, "calibration": vcal}
        if vcal["feasible"]:
            entry["anchor_check"] = anchor_reproduction_check(vcal, k_counts, rng, replicates=1000)
            entry["designs"] = []
            for d in S2_REFERENCE_DESIGNS:
                per = {"k": d[0], "m": d[1], "g": d[2], "eligible_population": int(len(eligible[d[0]]))}
                for ratio in S2_RATIOS:
                    r = _search(
                        vcal, d, ratio, eligible[d[0]], experiments, rng, None if ratio == 1.0 else 0.8
                    )
                    if ratio == 1.0:
                        per["null_false_positive_rate_by_n"] = {
                            str(c["n_blocks"]): c["interval_excludes_zero"] for c in r["cells"]
                        }
                    else:
                        per[f"n_for_80_percent_at_ratio_{ratio}"] = r["smallest_n_for_power"]["0.8"][
                            "n_blocks"
                        ]
                entry["designs"].append(per)
        sens.append(entry)
        print(f"variant {name}: done ({time.time() - t0:.0f} s)", flush=True)
    out["sensitivity_one_factor_at_a_time"] = sens
    # the registered check of the exact bootstrap against the committed resampling one
    agree = []
    for d in S2_REFERENCE_DESIGNS:
        for ratio in (1.0, 0.5):
            for n in (20, 100):
                if n <= len(eligible[d[0]]):
                    arm = block_arm(cal, ratio)
                    agree.append(
                        {
                            "ratio": ratio,
                            **bootstrap_agreement_check(
                                cal,
                                arm["mu_block"],
                                d,
                                n,
                                eligible[d[0]],
                                S2_BOOTSTRAP_CHECK_EXPERIMENTS,
                                rng,
                            ),
                        }
                    )
    out["exact_bootstrap_against_the_committed_one"] = agree
    out["committed_table_against_the_eligible_population"] = flag_the_committed_table(
        {k: int(len(v)) for k, v in eligible.items()}
    )
    out["power_is_a_probability"] = (
        "every size here is the smallest searched number of compared blocks at which the simulated committed "
        "analysis detects the assumed difference with the stated probability, IF the calibrated assumptions "
        "hold. It is never a guarantee that an experiment of that size decides clause 2; at 80% one such "
        "experiment in five fails to detect a real difference of the assumed size"
    )
    out["seconds"] = round(time.time() - t0, 1)
    return out


def flag_the_committed_table(
    eligible_by_k: dict[int, int], results_dir: Path = RESULTS_DIR
) -> dict[str, Any]:
    """lane-design's committed sizes, each set against the eligible population of its own design. Nothing in
    the committed file is changed; the flags are written here, beside it."""
    committed = json.loads((results_dir / f"{RESULT}.json").read_text())
    flagged = []
    for r in committed["sweep"]["rows"]:
        k = r["settings"]["elements_per_block"]
        cap = eligible_by_k.get(k)
        if cap is None:
            cap = eligible_by_k[max(x for x in eligible_by_k if x <= k)]
        for t, n in r["n_for_power"].items():
            if n is not None and n > cap:
                flagged.append(
                    {
                        "configuration": r["configuration"],
                        "ratio": r["ratio"],
                        "power": t,
                        "n_blocks": n,
                        "eligible": cap,
                    }
                )
    ranges_flagged = {
        ratio: {"largest": v["largest"], "reference": v["reference_configuration"]}
        for ratio, v in committed["sample_size_range_by_ratio"].items()
        if (v["largest"] or 0) > eligible_by_k[1] or (v["reference_configuration"] or 0) > eligible_by_k[1]
    }
    return {
        "entries_above_their_eligible_population": flagged,
        "count": len(flagged),
        "ranges_whose_ends_exceed_the_carrying_blocks": ranges_flagged,
        "reading": (
            "each of these is a size no experiment on this tier can have: infeasible, not a number to plan "
            "with. "
            "The committed table is kept as the record of what that run computed"
        ),
    }


def cheapest_with_a_calibrated_null(result: dict[str, Any]) -> list[dict[str, Any]]:
    """Derived AFTER the registered run, from its own null rule and nothing else: per ratio, the cheapest
    design reaching 80% whose null cell at the same number of blocks is calibrated. The registered
    `cheapest_feasible_design_at_80_percent` is kept as computed; the registration says a design with an
    anti-conservative null has sizes to be read as optimistic, and this names the cheapest one that has not.
    """
    null = {(c["k"], c["m"], c["g"], c["n_blocks"]): c for c in result["null_calibration"]["cells"]}
    out = []
    for row in result["sample_size_at_equal_cost"]:
        ok = []
        for x in row["designs"]:
            if x["n_blocks"] == "infeasible":
                continue
            cell = null[(x["k"], x["m"], x["g"], x["n_blocks"])]
            if cell["calibrated"]:
                ok.append({**x, "null_false_positive_rate": cell["false_positive_rate"]})
        reg = row["cheapest_feasible_design_at_80_percent"]
        out.append(
            {
                "ratio": row["ratio"],
                "registered_cheapest": reg,
                "registered_cheapest_null_false_positive_rate": (
                    null[(reg["k"], reg["m"], reg["g"], reg["n_blocks"])]["false_positive_rate"]
                    if reg
                    else None
                ),
                "feasible_designs": sum(1 for x in row["designs"] if x["n_blocks"] != "infeasible"),
                "feasible_with_a_calibrated_null": len(ok),
                "cheapest_with_a_calibrated_null": (
                    min(ok, key=lambda x: (x["assay_cost_elements"], x["n_blocks"])) if ok else None
                ),
            }
        )
    return out


@mf.depends_on_models("alphagenome")  # the scored-element set the windows are drawn over is the sweep's
def manifest_calibrated(experiments: int, seed: int) -> dict[str, Any]:
    arm = _measured_arm_module()
    base = arm.manifest(json.loads((RESULTS_DIR / f"{MEASURED_ARM}.json").read_text())["chromosomes"])
    inputs = list(base["inputs"])
    for name in (MEASURED_ARM, RESULT):
        p = RESULTS_DIR / f"{name}.json"
        if p.exists():
            inputs.append(mf.input_entry(p, partition=None))
    base["inputs"] = inputs
    base["sources"] = list(base["sources"]) + [
        {
            "accession": "this repository, data/results/clause2_design_power.json at 70801de",
            "version": "git 70801de: lane-design's committed table, read for the reproduction and the flags",
        }
    ]
    base["parameters"] = {
        **base["parameters"],
        "seed": seed,
        "experiments_per_cell": experiments,
        "designs_k": list(S2_K),
        "designs_m": list(S2_M),
        "designs_g": list(S2_G),
        "n_grid": list(S2_N_GRID),
        "budgets": list(S2_BUDGETS),
        "ratios": list(S2_RATIOS),
        "false_positive_grid": list(S2_FALSE_POSITIVE),
        "reference_false_positive": S2_REFERENCE_FALSE_POSITIVE,
        "reference_sensitivity_column": S2_REFERENCE_SENSITIVITY_COLUMN,
        "chromosome_resamples": S2_CHROMOSOME_RESAMPLES,
        "anchor_bootstrap": S2_ANCHOR_BOOTSTRAP,
        "anchor_tolerance": S2_ANCHOR_TOLERANCE,
        "icc_tolerance": S2_ICC_TOLERANCE,
    }
    base["exclusions"] = list(base["exclusions"]) + [
        "the model arm's -27.25 points and its per-block dispersion are not used anywhere",
        "no model answer is read: the element archive is opened for which elements are scored",
        "no size above a design's eligible population is simulated or printed as a number",
    ]
    return base


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
    ap.add_argument(
        "--calibrated",
        action="store_true",
        help="item 12 S2: run the calibrated model registered in S2_REGISTRATION (writes a new result)",
    )
    ap.add_argument("--experiments", type=int, default=S2_EXPERIMENTS, help="simulated experiments per cell")
    ap.add_argument(
        "--read-calibrated",
        action="store_true",
        help="item 12 S2: print the derived cheapest-with-a-calibrated-null view of the committed result",
    )
    args = ap.parse_args(argv)
    if args.read_calibrated:
        got = json.loads((RESULTS_DIR / f"{RESULT_CALIBRATED}.json").read_text())
        print(json.dumps(cheapest_with_a_calibrated_null(got), indent=1))
        return 0
    if args.calibrated:
        out = collect_calibrated(args.experiments, S2_SEED)
        if args.no_save:
            print("not saved")
        else:
            man = manifest_calibrated(args.experiments, S2_SEED)
            print(f"saved {save_result(RESULT_CALIBRATED, out, manifest=man)}")
        print("anchor reproduced:", (out.get("anchor_reproduced") or {}).get("reproduced"))
        print("null calibration:", (out.get("null_calibration") or {}).get("summary"))
        return 0
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
