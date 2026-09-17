# SPDX-License-Identifier: AGPL-3.0-or-later
"""Comparing two groups without fooling yourself: standardisation, with the imbalance printed.

Four readings had to be corrected or withdrawn on 2026-09-16 and 2026-09-17, and every one of them was
the same mistake in different clothes:

- a tier's rate read against another tier's rate, when the two differ tenfold in distance to a promoter
  (the 882-block reading, withdrawn);
- a rate conditioned on "was this measured", when one arm was measured more than the other (the locus
  benchmark's direction claim, withdrawn; the CRISPRi lane's coverage check, which passed only because
  it was made to look);
- a matched comparison implemented by pooling every control once per target, which inflates the sample
  and manufactures the p-value (+0.284 at p 0.00014 where the answer is +0.100 at p 0.18);
- a calibration whose curve transfers while its level does not, because the population's base rate
  differs (the calibration lane, failed as registered).

Each was found by hand, after publication. This module exists so the next one is found by the code:
`standardised` computes the comparison the right way, and returns the covariate medians of both groups
and the coverage of both arms **in the same result**, so an imbalance is a number in the output rather
than something a reader has to think to ask for.

    from genomeos.compare import standardised, Strata

    strata = Strata(length=(200, 400), gc=(0.35, 0.45, 0.55), tss=(1_000, 5_000, 20_000, 100_000))
    out = standardised(targets, controls, strata, hit="moved")
    out["matched"]["difference"]   # direct standardisation, not a pooled rate
    out["imbalance"]["tss"]        # the medians that decide whether the comparison means anything

Nothing here is specific to the unknown space; it is arithmetic with a memory.
"""

from __future__ import annotations

import bisect
import math
from collections import defaultdict
from dataclasses import dataclass, field
from statistics import median
from typing import Any

Row = dict[str, Any]


@dataclass(frozen=True)
class Strata:
    """Bin edges per covariate. A row is assigned the tuple of its bins, and only rows sharing a tuple
    are compared. Coarse bins are honest as long as the imbalance is reported beside the result."""

    bins: dict[str, tuple[float, ...]] = field(default_factory=dict)

    def __init__(self, **covariates: tuple[float, ...]) -> None:
        object.__setattr__(self, "bins", dict(covariates))

    def key(self, row: Row) -> tuple[int, ...]:
        return tuple(bisect.bisect_right(cuts, row[name]) for name, cuts in self.bins.items())

    @property
    def covariates(self) -> tuple[str, ...]:
        return tuple(self.bins)


def difference(a_hits: float, a_n: float, b_hits: float, b_n: float) -> dict[str, Any]:
    """Difference in two shares with a one-sided p and a 95% upper bound, normal approximation.

    Both arms are given as (hits, n) so a standardised control arm - whose n is the number of matched
    targets, not the number of control rows - can be passed without pretending it has more precision
    than it does.
    """
    if not a_n or not b_n:
        return {"a": None, "b": None, "difference": None, "p_one_sided": None, "upper_95": None}
    pa, pb = a_hits / a_n, b_hits / b_n
    se = math.sqrt(pa * (1 - pa) / a_n + pb * (1 - pb) / b_n)
    z = (pa - pb) / se if se else 0.0
    return {
        "a": round(pa, 4),
        "b": round(pb, 4),
        "difference": round(pa - pb, 4),
        "p_one_sided": round(0.5 * math.erfc(z / math.sqrt(2)), 6),
        "upper_95": round(pa - pb + 1.96 * se, 4),
    }


def stratum_rates(rows: list[Row], strata: Strata, hit: str = "moved") -> dict[tuple[int, ...], float]:
    """Each stratum's hit rate, computed once over the whole pool.

    Once, not once per target: computing it inside the target loop is both quadratic and the source of
    the pooling defect, since a stratum with many rows then counts many times over.
    """
    hits: dict[tuple[int, ...], int] = defaultdict(int)
    seen: dict[tuple[int, ...], int] = defaultdict(int)
    for r in rows:
        key = strata.key(r)
        seen[key] += 1
        hits[key] += bool(r[hit])
    return {k: hits[k] / n for k, n in seen.items() if n}


def imbalance(targets: list[Row], controls: list[Row], strata: Strata) -> dict[str, Any]:
    """The medians of every covariate in both groups, and the ratio between them.

    This is the number that would have caught the 882-block reading a day earlier: the arms differed by
    a factor of ten in distance to a promoter, and nothing in the comparison said so.
    """
    out = {}
    for name in strata.covariates:
        t = [r[name] for r in targets if r.get(name) is not None]
        c = [r[name] for r in controls if r.get(name) is not None]
        mt, mc = (median(t) if t else None), (median(c) if c else None)
        out[name] = {
            "target_median": round(mt, 4) if mt is not None else None,
            "control_median": round(mc, 4) if mc is not None else None,
            "ratio": round(mt / mc, 3) if mt and mc else None,
        }
    return out


def standardised(
    targets: list[Row], controls: list[Row], strata: Strata, hit: str = "moved"
) -> dict[str, Any]:
    """Compare two groups by direct standardisation, with the imbalance and the coverage beside it.

    `raw` is the unstandardised difference, kept because it is what an unmatched claim would have said.
    `matched` averages each target's stratum rate over the targets, which is the comparison that holds
    the covariates fixed. `dropped` counts the targets with no control in their stratum: a comparison
    that keeps a tenth of its targets is not the comparison it looks like, so the number is returned
    rather than logged.
    """
    rates = stratum_rates(controls, strata, hit)
    kept: list[Row] = []
    control_hits = 0.0
    for r in targets:
        rate = rates.get(strata.key(r))
        if rate is None:
            continue
        kept.append(r)
        control_hits += rate
    return {
        "targets": len(targets),
        "controls": len(controls),
        "targets_matched": len(kept),
        "dropped_for_want_of_a_control": len(targets) - len(kept),
        "strata": {name: list(cuts) for name, cuts in strata.bins.items()},
        "raw": difference(
            sum(1 for r in targets if r[hit]), len(targets), sum(1 for r in controls if r[hit]), len(controls)
        ),
        "matched": difference(sum(1 for r in kept if r[hit]), len(kept), round(control_hits, 6), len(kept)),
        "imbalance": imbalance(targets, controls, strata),
    }
