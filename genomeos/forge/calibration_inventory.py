# SPDX-License-Identifier: AGPL-3.0-or-later
"""The inventory that could become BioForge's calibration record, and the rule that decides if it is one.

Area H, 2026-10-10 (docs/BIOFORGE-CONFIDENCE.md). A design answer carries no probability because no
calibration record exists: no population of design answers with published outcomes. This module states,
before anything is counted, what such a population would have to be: a set of (design question,
BioForge's answer, published outcome) triples whose outcomes the model did not already hold. It fixes the
smallest population at which a calibration is worth computing, from the width of an exact binomial
interval, and it fixes the rules that mark a triple circular. The census that applies them follows the
registration in this file; it never fabricates an outcome, and it computes no match rate on a population
the registration calls circular, because an answer that reproduces its own training phenotype is not a
prediction.
"""

from __future__ import annotations

import math
from typing import Any

REGISTERED_ON = "2026-10-10"
CONFIDENCE_LEVEL = 0.95
MAX_INTERVAL_WIDTH = 1 / 3


def _binom_cdf(k: int, n: int, p: float) -> float:
    """P(X <= k) for X ~ Binomial(n, p), summed exactly in floating point (n is small here)."""
    if k < 0:
        return 0.0
    if k >= n:
        return 1.0
    return min(1.0, sum(math.comb(n, i) * p**i * (1.0 - p) ** (n - i) for i in range(k + 1)))


def _bisect(f, lo: float = 0.0, hi: float = 1.0, steps: int = 80) -> float:
    """The root of a function increasing on [lo, hi]."""
    for _ in range(steps):
        mid = (lo + hi) / 2
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def clopper_pearson(k: int, n: int, level: float = CONFIDENCE_LEVEL) -> tuple[float, float]:
    """The exact (Clopper-Pearson) two-sided interval for a binomial rate, k successes of n."""
    if n <= 0 or not 0 <= k <= n:
        raise ValueError(f"need 0 <= k <= n and n > 0, got k={k}, n={n}")
    a = (1.0 - level) / 2
    # lower: the p at which P(X >= k) = a; P(X >= k) rises with p
    lo = 0.0 if k == 0 else _bisect(lambda p: (1.0 - _binom_cdf(k - 1, n, p)) - a)
    # upper: the p at which P(X <= k) = a; P(X <= k) falls with p, so bisect its negation
    hi = 1.0 if k == n else _bisect(lambda p: a - _binom_cdf(k, n, p))
    return lo, hi


def worst_width(n: int, level: float = CONFIDENCE_LEVEL) -> float:
    """The widest exact interval any outcome count can give at n: the width a record must plan for."""
    return max(hi - lo for lo, hi in (clopper_pearson(k, n, level) for k in range(n + 1)))


def minimum_triples(max_width: float = MAX_INTERVAL_WIDTH, level: float = CONFIDENCE_LEVEL) -> int:
    """The smallest n at which every possible count gives an exact interval narrower than `max_width`."""
    n = 1
    while worst_width(n, level) >= max_width:
        n += 1
    return n


# Fixed here and pinned by tests/test_forge_calibration_inventory.py against `minimum_triples()`.
MINIMUM_INDEPENDENT_TRIPLES = 38

REGISTRATION: dict[str, Any] = {
    "written": (
        "2026-10-10, by lane-forgecal, after an exploratory read of which programs hold design and"
        " experiment blocks and whose citations their modules repeat, and BEFORE the census below was"
        " coded, before any inverse design question was posed to BioForge and before any match was"
        " counted. The exploratory read already suggested that no triple would be independent; the"
        " minimum below is a function of the interval width alone and was not chosen with any count in"
        " view, so it could not be tuned to make the census pass or fail"
    ),
    "the_triple": (
        "(design question, BioForge's answer, published outcome): a question BioForge can be asked as"
        " a `design` block over an organism program (targets and keeps in assert grammar, a candidate"
        " list of knockouts or additions), the perturbation BioForge reports as its best feasible"
        " candidate, and the perturbation a cited publication reports for the same phenotype"
    ),
    "the_match": (
        "exact set equality between the knockouts and additions BioForge reports and those the"
        " publication reports. A tie between equally good candidates that BioForge breaks by label"
        " order counts as what BioForge reported; the tie is recorded beside the triple"
    ),
    "the_minimum": {
        "n": MINIMUM_INDEPENDENT_TRIPLES,
        "rule": (
            "the smallest n at which the exact 95% Clopper-Pearson interval on the match rate is"
            " narrower than 1/3 for EVERY possible count k in 0..n (worst case k = n/2)"
        ),
        "why_one_third": (
            "a probability on a design answer is worth stating only if the record behind it can tell"
            " a reader at least whether the answer is more often right than wrong. An interval"
            " narrower than 1/3 cannot contain both 1/3 and 2/3, so whatever the count, the record"
            " rules out at least one of 'mostly wrong' (rate <= 1/3) and 'mostly right' (rate >= 2/3)."
            " Below that n some count leaves both readings open and the record says nothing a reader"
            " can act on. A narrower width (0.25 needs 67, 0.2 needs 104) would be better and is not"
            " required: this is the floor, not the target"
        ),
        "what_it_buys": (
            "calibration-in-the-large only: one rate for one population of design answers. A"
            " confidence that varies across answers would need at least this many independent"
            " triples PER reliability bin, and discrimination needs both outcomes present"
        ),
    },
    "circular_if_any": {
        "same_source": (
            "the outcome's publication (first author's surname and year) is cited by an element of"
            " the model the question runs over (a decision, timer, rule, cell type or other evidence"
            " bearing element, imports included), or by a script the program names in its own text as"
            " the source of its gates (scripts/intestine_registration.py for the crypt module). The"
            " search then traverses a rule written from the answer"
        ),
        "acceptance_check": (
            "the outcome is shipped as an `experiment` block run against the same model by `bio test`"
            " in CI, so the model is held to reproduce it: an answer that matches it reproduces a"
            " phenotype the model was built and tested to give"
        ),
        "authored_frame": (
            "the question's candidate list was written by hand by an author who held the answer, and"
            " contains it (every shipped `design` block)"
        ),
        "evidence_is_the_model": (
            "the design block's own evidence names the model's requirements, not a publication, as the"
            " source of its expected answer"
        ),
    },
    "independent": (
        "a triple none of the four circularity rules marks. A triple that no rule marks is not proven"
        " independent: the rules catch a shared citation and a shipped check, not an author who read"
        " the paper and cited another. The count of independent triples is therefore an upper bound"
    ),
    "eligible_shapes": (
        "an `experiment` with a published citation (evidence experimental) can be posed as an inverse"
        " design question only when BioForge's grammar can express it: knockouts at the zygote over a"
        " program that grows an organism (stages declared). Additions forced at a later time, ablations"
        " given as an environment, and programs with no organism are counted and set aside"
    ),
    "decision": (
        f"if the independent triples number at least {MINIMUM_INDEPENDENT_TRIPLES}, register a"
        " calibration-in-the-large on them (outcome, population and method named, as genomeos.certainty"
        " .Calibration requires) and compute it. Otherwise record the shortfall and compute no rate:"
        " no inverse design is run and no match is counted on circular triples, because a rate on them"
        " measures retrieval of the model's own rules"
    ),
    "never": (
        "an outcome is never written that is not in the repository with its citation; a design answer"
        " that matches the mutant it was derived from is never counted as an independent outcome"
    ),
}
