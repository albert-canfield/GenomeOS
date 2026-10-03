# SPDX-License-Identifier: AGPL-3.0-or-later
"""Why 44 of N1's 56 candidate factors have no responders, from the committed result alone.

This module contains no measurement. Every function here is pure: it takes the figures that
`data/results/n1_result_amendment_2.json` already records and the text of the committed module that
wrote them, and it returns arithmetic and line numbers. Nothing in it opens a study file, parses a
p-value, or touches a perturbed row. The question it serves is a question about a finished run, not a
new reading of the data, and the distinction is the point: N1b step 2's value is that outcome
blindness has not been spent, so an explanation of N1's coverage must be derivable without spending
any of it.

Three things are worth separating, because one label in the committed result covers only two of them:

* `no_responders` is emitted when a factor's responder class is EMPTY (`stratified_auroc`, `pos == 0`);
* `all_responders` and `no_stratum_with_both_classes` are separate return values of the same function,
  so "every gene responded" and "decile matching left no comparison group" are counted apart and are
  not inside the 44;
* a gene whose p-value is missing leaves the factor's universe and is counted as `genes_missing`.

What `no_responders` therefore still conflates is exactly two causes: a perturbation for which the
producers called nothing differentially expressed anywhere, and a perturbation whose differentially
expressed genes all lie outside the analysis universe. `missing_decomposition_field` names the field
that would separate them and says what reading it would cost.
"""

from __future__ import annotations

import math
import re
from collections.abc import Sequence
from typing import Any

#: The floor the registration froze: defined paired differences needed before an estimate is reported.
FLOOR = 30

#: Exact source strings looked for in the committed analysis module, so that every line number this
#: analysis reports is FOUND rather than typed in. A rename or a move changes the numbers; it cannot
#: leave a stale number looking right. Each key is what the site is, each value the text that marks it.
SITES = {
    "responder_predicate": "labels = [pvals[genes[j]] < AD_LEVEL for j in keep]",
    "ad_level_constant": "AD_LEVEL = 0.05",
    "no_responders_emission": 'return None, "no_responders", 0',
    "all_responders_emission": 'return None, "all_responders", 0',
    "no_stratum_emission": 'return None, "no_stratum_with_both_classes", 0',
    "values_pass_row_skip": "if gene not in wanted_rows:",
    "labels_pass_row_count_discarded": 'seen = Counter(labels["rows"])',
    "responds_t_constant": "RESPONDS_T = 3.0",
    "pooled_t_definition": "def pooled_t(",
    "run_v3_definition": "def run_v3(",
    "factor_result_v3_definition": "def factor_result_v3(",
    "stratified_auroc_definition": "def stratified_auroc(",
}

#: The keys that tell the two endpoints apart in a committed result, without reading any code. Each
#: `factor_result_*` emits a key the other does not, so the result's own bytes name the path that wrote
#: them: the Anderson-Darling path counts the genes whose p-value was missing, the T path counts the
#: genes whose response statistic was not finite.
ENDPOINT_KEYS = {
    "published_anderson_darling": "genes_missing",
    "derived_t_statistic": "genes_nonfinite_response",
}


def source_sites(text: str, sites: dict[str, str] | None = None) -> dict[str, Any]:
    """Where each of `SITES` occurs in a module's text: 1-based line numbers, in file order.

    A site that does not occur returns an empty list rather than an error, so the absence of a site is
    itself reportable: that is how `t_is_absent_from` establishes that the T statistic is not on the
    amendment-2 path instead of asserting it.
    """
    lines = text.splitlines()
    out: dict[str, Any] = {}
    for name, needle in (sites or SITES).items():
        out[name] = [i for i, line in enumerate(lines, 1) if needle in line]
    return out


def function_span(text: str, name: str) -> tuple[int, int] | None:
    """The 1-based line span of a top-level `def name(` in a module's text, or None when it has none.

    The end is the last line before the next top-level statement, which is what a reader needs to check
    whether a name is used inside a function without running anything.
    """
    lines = text.splitlines()
    start = None
    for i, line in enumerate(lines, 1):
        if re.match(rf"^def {re.escape(name)}\(", line):
            start = i
            break
    if start is None:
        return None
    for j in range(start + 1, len(lines) + 1):
        line = lines[j - 1]
        if line and not line[0].isspace() and not line.startswith(")"):
            return start, j - 1
    return start, len(lines)


def t_is_absent_from(text: str, function: str, markers: Sequence[str]) -> dict[str, Any]:
    """Whether a function's own body contains any of `markers`, with the lines where they DO occur.

    Used to show, checkably, that `run_v3` neither thresholds the derived statistic nor pools it: the
    marker lines are all outside its span, and the span is reported beside them so a reader can compare
    the two sets rather than take the conclusion.
    """
    span = function_span(text, function)
    if span is None:
        return {"function": function, "span": None, "absent": None, "marker_lines": {}}
    lo, hi = span
    found = {m: [i for i, line in enumerate(text.splitlines(), 1) if m in line] for m in markers}
    inside = {m: [i for i in lines if lo <= i <= hi] for m, lines in found.items()}
    return {
        "function": function,
        "span": [lo, hi],
        "absent": all(not v for v in inside.values()),
        "marker_lines": found,
        "marker_lines_inside_the_span": {m: v for m, v in inside.items() if v},
    }


def endpoint_from_result_keys(factor_record: dict[str, Any]) -> dict[str, Any]:
    """Which endpoint a committed per-factor record came from, from its key set alone.

    Independent of the code argument: whichever of `ENDPOINT_KEYS` the record carries names the writer,
    and carrying both or neither is returned as `ambiguous` rather than resolved by preference.
    """
    present = {name: key for name, key in ENDPOINT_KEYS.items() if key in factor_record}
    return {
        "keys_present": sorted(present.values()),
        "keys_absent": sorted(k for k in ENDPOINT_KEYS.values() if k not in factor_record),
        "endpoint": next(iter(present)) if len(present) == 1 else "ambiguous",
    }


def concentration(counts: Sequence[int], floor: int = FLOOR) -> dict[str, Any]:
    """The responder counts' first two moments and what an EVEN spread at the same mean would give.

    This is the figure that decides what the shortfall means. Under any model that spreads responders
    evenly across factors -- a homogeneous Poisson at the observed mean -- the expected number of
    factors with at least one responder follows from the mean alone, and is compared here against the
    floor. The observed zero count is placed in standard deviations of that model's own zero count, so
    the model is refuted or not by the committed counts rather than assumed away.
    """
    n = len(counts)
    total = sum(counts)
    mean = total / n
    var = sum((x - mean) ** 2 for x in counts) / (n - 1)
    zeros = sum(1 for x in counts if x == 0)
    p0 = math.exp(-mean)
    expected_zeros = n * p0
    sd = math.sqrt(n * p0 * (1 - p0))
    return {
        "factors": n,
        "responder_cells": total,
        "per_factor_mean": mean,
        "per_factor_variance_ddof1": var,
        "variance_to_mean_ratio": var / mean if mean else None,
        "zero_factors_observed": zeros,
        "defined_observed": n - zeros,
        "floor": floor,
        "even_spread": {
            "model": "homogeneous Poisson at the observed per-factor mean: every factor the same rate",
            "zero_probability": p0,
            "expected_zero_factors": expected_zeros,
            "expected_defined": n - expected_zeros,
            "clears_the_floor": (n - expected_zeros) >= floor,
            "sd_of_the_zero_count": sd,
            "observed_zeros_in_sd": (zeros - expected_zeros) / sd if sd else None,
        },
    }


def gamma_mixed_shape(mean: float, var: float) -> float | None:
    """The gamma mixing shape k of a gamma-mixed Poisson matched to a mean and variance, by moments.

    None when the variance does not exceed the mean, because then no gamma mixture fits and the even
    spread is not refuted.
    """
    return mean * mean / (var - mean) if var > mean else None


def gamma_mixed_zero_probability(mu: float, k: float) -> float:
    """P(count = 0) for a gamma-mixed Poisson of mean `mu` and mixing shape `k`."""
    return (k / (k + mu)) ** k


def gamma_mixed_universe(counts: Sequence[int], genes: int, universes: Sequence[int]) -> dict[str, Any]:
    """What a gamma-mixed Poisson fitted to the observed counts expects at larger gene universes.

    The mean scales with the universe and the mixing shape is held, which is the assumption to state
    rather than bury: it says a factor's own per-gene rate is the same outside the analysed genes as
    inside them. It is a MODEL-BASED extrapolation from twelve counts and it is never to be quoted
    without `distribution_free_bracket`, which gives the two ends it cannot exclude.
    """
    n = len(counts)
    mean = sum(counts) / n
    var = sum((x - mean) ** 2 for x in counts) / (n - 1)
    k = gamma_mixed_shape(mean, var)
    if k is None:
        return {"fits": False, "why": "the variance does not exceed the mean, so no gamma mixture fits"}
    rows = []
    for u in universes:
        mu = mean * u / genes
        p0 = gamma_mixed_zero_probability(mu, k)
        rows.append(
            {
                "universe_genes": u,
                "mean_responders_per_factor": mu,
                "zero_probability": p0,
                "expected_zero_factors": n * p0,
                "expected_defined": n * (1 - p0),
            }
        )
    return {
        "fits": True,
        "model": "gamma-mixed Poisson (negative binomial), shape and mean matched to the observed "
        "counts by moments; the mean scales with the universe and the shape is held",
        "shape_k": k,
        "mixing_coefficient_of_variation": 1 / math.sqrt(k),
        "assumption": "a factor's per-gene responder rate outside the analysed genes equals its rate "
        "inside them",
        "at_universes": rows,
    }


def gamma_mixed_universe_for_defined(
    counts: Sequence[int], genes: int, wanted_defined: int
) -> dict[str, Any]:
    """The gene universe the fitted gamma mixture would need before `wanted_defined` factors are defined.

    Reported to show the scale of the shortfall, not as a route: when the answer exceeds every gene the
    study holds, no rule about which genes to analyse reaches the floor.
    """
    n = len(counts)
    mean = sum(counts) / n
    var = sum((x - mean) ** 2 for x in counts) / (n - 1)
    k = gamma_mixed_shape(mean, var)
    if k is None or not 0 < wanted_defined < n:
        return {"solvable": False}
    target_p0 = (n - wanted_defined) / n
    mu = k * (target_p0 ** (-1 / k) - 1)
    return {
        "solvable": True,
        "wanted_defined": wanted_defined,
        "mean_responders_per_factor_needed": mu,
        "universe_genes_needed": genes * mu / mean,
    }


def exact_upper_rate(genes: int, observed: int = 0, level: float = 0.05) -> float:
    """The one-sided upper confidence bound on a per-gene responder rate, exact and binomial.

    `observed = 0` is the case that matters: a factor with no responder among `genes` genes still
    admits any rate whose probability of giving none is at least `level`. Solved by bisection on the
    binomial tail so the bound is the definition's own and not a normal approximation of it.
    """
    if observed != 0:
        raise NotImplementedError("the zero-count bound is what this analysis needs")
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if (1 - mid) ** genes < level:
            hi = mid
        else:
            lo = mid
    return hi


def distribution_free_bracket(
    counts: Sequence[int], genes: int, added_genes: int, level: float = 0.05
) -> dict[str, Any]:
    """The two ends the committed counts cannot separate, stated together with the model's answer.

    LOW end: each factor's own observed rate, so a factor with none stays at zero and the defined count
    never moves, whatever the universe. HIGH end: a factor with none is at its exact upper confidence
    bound, where the chance that at least one of the added genes responds is large. Both are honest
    readings of the same committed counts; the quantity that would separate them is named by
    `missing_decomposition_field` and reading it spends outcome blindness.
    """
    n = len(counts)
    zeros = sum(1 for x in counts if x == 0)
    rate = exact_upper_rate(genes, 0, level)
    p_any = 1 - (1 - rate) ** added_genes
    return {
        "low_end": {
            "assumption": "each factor's per-gene rate is its own observed rate",
            "defined_at_any_universe": n - zeros,
            "reading": "a factor with no responder among the analysed genes stays undefined however "
            "many genes are added",
        },
        "high_end": {
            "assumption": f"a factor with no responder sits at its exact one-sided {1 - level:.0%} "
            "upper confidence bound",
            "genes_observed": genes,
            "upper_per_gene_rate": rate,
            "responders_per_observed_universe_at_that_rate": rate * genes,
            "added_genes": added_genes,
            "probability_at_least_one_responder_among_the_added_genes": p_any,
            "reading": "a case exists in which most of the undefined factors would become defined, and "
            "the committed counts do not exclude it",
        },
        "what_separates_them": "whether a factor with none among the analysed genes has a rate near "
        "zero or was merely unlucky in them, which only its p-values at the other genes would say",
    }


def missing_decomposition_field() -> dict[str, Any]:
    """The field that would decompose the undefined factors, and what forming it would cost.

    One field, one cost, stated so neither is mistaken for the other: the count exists nowhere in the
    committed record, and forming it parses published p-values at genes the run deliberately skipped.
    """
    return {
        "field": "responders_outside_universe",
        "would_separate": [
            "a perturbation for which the producers called nothing differentially expressed anywhere",
            "a perturbation whose differentially expressed genes all lie outside the analysis universe",
        ],
        "why_it_was_never_formed": "the values pass skips every row outside the wanted rows after "
        "reading its first field, so a p-value outside the universe is never converted and no count "
        "of such calls exists in the result",
        "cost_to_form": "it parses published adjusted p-values at genes the authorised run did not "
        "read, so it spends outcome blindness and was not taken here",
    }
