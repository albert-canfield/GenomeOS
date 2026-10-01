# SPDX-License-Identifier: AGPL-3.0-or-later
"""The wiring diagnostic: do the existing element-to-gene deletion assignments beat matched rewired ones?

The CRISPRi result (`crispri.score_published`, data/results/crispri_published.json) adds the sweep's cached
deletion block -- `top_target` and `deletion_drop`, what deleting the registry elements under a tested
element does to the measured gene on K562's own track -- to an activity-and-distance model. On the 10,356
K562 training pairs, hold-one-chromosome-out, the benchmark's AUPRC goes from 0.5068 to 0.7241. That
benchmark has already been examined, and the deletion features were chosen after the single predictors
were read on these same pairs, so anything this module supports is development evidence, not
independent validation.

The question here is narrower. Does the *assignment* of a cached value to the measured gene carry
information beyond what a matched alternative assignment carries? A rewiring keeps every evaluated pair,
its CRISPRi label, its own distance and activity features, the split and the estimator, and changes only
which gene's cached deletion value is attached to the pair. A pair with no admissible alternative keeps
its real value in every rewiring and is counted; no pair is ever added, dropped or moved.

Three rewirings are defined (NULLS); each answers a different question:

- `element_kept` (primary): the element is kept, so its activity and its links are kept, and the value
  attached is the same element's cached value on another gene that the same element was measured
  against, at a matched distance to TSS. It asks whether the gene the value is assigned to matters.
- `gene_kept`: the gene is kept and the value attached is the same gene's cached value from another
  element measured against it, at a matched distance and activity decile. It asks whether the element
  the value comes from matters, with everything about the gene held.
- `strata`: values are permuted among links in the same distance bin and activity decile, always from a
  different element. It is a reference floor: it keeps distance, activity and coverage and nothing about
  the element or the gene, so it answers whether the value carries anything beyond those strata.

The links a rewiring may touch are the population: covered pairs whose gene the sweep's cache answered
(`deletion_answered`), so coverage is held exactly. Within a population, links are matched in strata: the
null's group (the element, the gene, or the activity decile), optionally a K562 expression class
(SCHEMES), and a distance bin of width `width_log10` in log10 distance with a random offset drawn per
rewiring. Inside a stratum the values are permuted so that no link receives a value of its own gene (or,
for the element-breaking nulls, its own element); a permutation keeps every gene's in-degree and every
element's number of links exactly.

This module holds the matching, the diagnostics and the thresholds they are judged by. THRESHOLDS,
SCHEMES and the regulated floors were committed before any diagnostic was computed on the real pairs.
Nothing in this file computes a performance metric of a rewired assignment; the analysis is registered
separately and run once, after the registration has been checked.
"""

from __future__ import annotations

import bisect
import math
import random
from collections import Counter, defaultdict
from collections.abc import Hashable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from genomeos.attribution import crispri

BLOCK = ("top_target", "deletion_drop")  # the cached deletion block a rewiring moves; nothing else moves
SEED = 20261001
FEASIBILITY_DRAWS = 20  # rewirings drawn for the matching diagnostics; the run's first 20 are these draws
DERANGE_TRIES = 1_000  # uniform rejection draws before the constructive fallback
ACTIVITY_BINS = 10  # deciles of crispri's log_activity over the population
EXPRESSION_BINS = 4  # quartiles of log1p K562 expression over the population's genes with a value

NULLS: dict[str, str] = {
    "element_kept": (
        "primary: the element (so its activity and its links) is kept; the value attached is the same "
        "element's cached value on another gene the same element was measured against, at a matched distance"
    ),
    "gene_kept": (
        "secondary: the gene is kept; the value attached is the same gene's cached value from another "
        "element measured against it, at a matched distance and activity decile"
    ),
    "strata": (
        "secondary, a reference floor: values permuted among links in the same distance bin and activity "
        "decile, always from a different element; element and gene are not held"
    ),
}

#: Tried in this order for the primary null on the training pairs; the first that meets every threshold
#: is the one registered. Fixed before any diagnostic was computed. A scheme with expression classes
#: matches only within the same K562 expression quartile (an unknown expression is its own class).
SCHEMES: tuple[dict[str, Any], ...] = (
    {"name": "S1", "width_log10": 0.20, "expression_classes": False},
    {"name": "S2", "width_log10": 0.20, "expression_classes": True},
    {"name": "S3", "width_log10": 0.30, "expression_classes": False},
    {"name": "S4", "width_log10": 0.30, "expression_classes": True},
)

#: Adequacy thresholds, fixed 2026-10-01 before any diagnostic was computed on the real pairs. Shares
#: are of the population (covered pairs with a cached value); balance is measured on moved links,
#: pooled over the draws, within each label class, attached against real.
THRESHOLDS: dict[str, float] = {
    # 1. meaningful rewiring
    "movable_share_min": 0.50,  # links with >= 1 admissible alternative; overall and regulated
    "moved_share_min": 0.50,  # links moved per rewiring, mean over draws; overall and regulated
    "two_alternatives_share_min": 0.50,  # movable links with >= 2 distinct alternatives; overall, regulated
    # 2. balance, attached against real, in each label class
    "smd_max": 0.10,  # |standardised mean difference|: distance, activity, expression, links per gene
    "ks_max": 0.10,  # two-sample Kolmogorov-Smirnov D: distance, expression
    "expression_known_min": 0.80,  # moved links whose genes both have a cached value (or are one gene)
    "expression_unknown_gap_max": 0.05,  # |share without a value, attached - real|
}
#: Regulated links that must move per rewiring (mean over draws): the training pairs, then held-out K562.
REGULATED_FLOOR = {"training": 100, "heldout_k562": 30}
COVARIATES = ("distance", "activity", "expression", "links_per_gene")
KS_COVARIATES = ("distance", "expression")
CLASSES = ("regulated", "not_regulated")


@dataclass(frozen=True)
class Link:
    """One population pair as the matching sees it. `regulated` is read by the diagnostics only."""

    index: int  # the pair's position in the evaluated list, which never changes
    element: tuple[str, int, int]
    gene: str
    chrom: str
    position: float  # log10 distance to the gene's TSS, floored at crispri.MIN_DISTANCE
    activity: float  # crispri's log_activity of the element
    expression: float | None  # log1p K562 expression of the gene; None where nothing is cached
    regulated: bool


# --- the population ---------------------------------------------------------------------------------------


def element_of(p: crispri.Pair) -> tuple[str, int, int]:
    return (p.chrom, p.start, p.end)


def log10_distance(p: crispri.Pair) -> float:
    return math.log10(max(crispri.MIN_DISTANCE, p.distance))


def in_population(p: crispri.Pair) -> bool:
    """A covered pair whose gene the sweep's cache answered: the only pairs a rewiring may touch."""
    return bool(p.covered) and bool(p.features.get("deletion_answered"))


def links_of(pairs: Sequence[crispri.Pair], expression: dict[str, float]) -> list[Link]:
    """The population as links; `expression` maps a gene symbol to log1p expression, absent if unknown."""
    return [
        Link(
            index=i,
            element=element_of(p),
            gene=p.gene,
            chrom=p.chrom,
            position=log10_distance(p),
            activity=p.features["log_activity"],
            expression=expression.get(p.gene),
            regulated=p.regulated,
        )
        for i, p in enumerate(pairs)
        if in_population(p)
    ]


def cuts(values: Sequence[float], k: int) -> list[float]:
    """The k-quantile cut points (k - 1 of them), by linear interpolation."""
    v = sorted(values)
    if not v:
        return []
    out = []
    for i in range(1, k):
        q = i / k * (len(v) - 1)
        lo = math.floor(q)
        out.append(v[lo] + (v[min(lo + 1, len(v) - 1)] - v[lo]) * (q - lo))
    return out


def bin_of(value: float, edges: Sequence[float]) -> int:
    return bisect.bisect_right(edges, value)


def keys_for(
    links: Sequence[Link], null: str, expression_classes: bool
) -> tuple[list[Hashable], list[Hashable]]:
    """(group, distinct) per link: a link may receive a value only from its own group, never from a link
    sharing its distinct key (the gene for `element_kept`, the element otherwise)."""
    act = cuts([lk.activity for lk in links], ACTIVITY_BINS)
    expr = cuts(
        sorted({lk.gene: lk.expression for lk in links if lk.expression is not None}.values()),
        EXPRESSION_BINS,
    )
    groups: list[Hashable] = []
    distinct: list[Hashable] = []
    for lk in links:
        xc: Hashable = None
        if expression_classes:
            xc = "unknown" if lk.expression is None else bin_of(lk.expression, expr)
        ad = bin_of(lk.activity, act)
        if null == "element_kept":
            groups.append((lk.element, xc))
            distinct.append(lk.gene)
        elif null == "gene_kept":  # the expression class is the gene's own, so it changes nothing here
            groups.append((lk.gene, ad))
            distinct.append(lk.element)
        elif null == "strata":
            groups.append((ad, xc))
            distinct.append(lk.element)
        else:
            raise ValueError(f"unknown null {null!r}")
    return groups, distinct


# --- one rewiring -----------------------------------------------------------------------------------------


def derange(labels: Sequence[Hashable], rng: random.Random, tries: int = DERANGE_TRIES) -> list[int] | None:
    """A permutation `src` of range(n) with labels[src[k]] != labels[k] for every k, or None if none exists.

    One exists exactly when no label fills more than half the stratum. Uniform by rejection; if `tries`
    shuffles all fail, a valid one is built by laying the labels out in contiguous blocks (blocks and
    members in random order) and shifting by the largest block, which is valid but not uniform.
    """
    n = len(labels)
    if n < 2:
        return None
    count = Counter(labels)
    top = max(count.values())
    if 2 * top > n:
        return None
    order = list(range(n))
    for _ in range(tries):
        rng.shuffle(order)
        if all(labels[order[k]] != labels[k] for k in range(n)):
            return list(order)
    blocks: dict[Hashable, list[int]] = defaultdict(list)
    for k, lab in enumerate(labels):
        blocks[lab].append(k)
    names = list(blocks)
    rng.shuffle(names)
    seq: list[int] = []
    for name in names:
        members = blocks[name]
        rng.shuffle(members)
        seq.extend(members)
    src = [0] * n
    for pos, k in enumerate(seq):
        src[k] = seq[(pos + top) % n]
    return src


def rewire(
    links: Sequence[Link],
    groups: Sequence[Hashable],
    distinct: Sequence[Hashable],
    width: float,
    rng: random.Random,
) -> dict[int, int]:
    """One matched rewiring as {target pair index: source pair index}; a link not in it keeps its value.

    Strata are (group, distance bin); the bins have width `width` in log10 distance and a random offset
    drawn once per rewiring, so two links closer than `width` share a bin with probability 1 - gap/width.
    Labels are never read.
    """
    offset = rng.random() * width
    strata: dict[Hashable, list[int]] = {}
    for k, lk in enumerate(links):
        strata.setdefault((groups[k], math.floor((lk.position - offset) / width)), []).append(k)
    out: dict[int, int] = {}
    for members in strata.values():
        src = derange([distinct[k] for k in members], rng)
        if src is None:
            continue
        for pos, k in enumerate(members):
            out[links[k].index] = links[members[src[pos]]].index
    return out


def rng_for(pair_set: str, null: str, scheme: str, draw: int) -> random.Random:
    """The generator of one rewiring; string seeds are hashed by sha512, so this is the same everywhere."""
    return random.Random(f"{SEED}:{pair_set}:{null}:{scheme}:{draw}")


def attach(pairs: Sequence[crispri.Pair], rewiring: dict[int, int]) -> list[dict[str, float]]:
    """Every pair's feature row with the rewired block attached: the same pairs, in the same order.

    Only the BLOCK columns of a moved pair change, to its source pair's real values; every other column
    and every unmoved pair (a population link with no admissible alternative, or a pair outside the
    population) keeps its real value. Labels are not touched because they are not here.
    """
    rows = [dict(p.features) for p in pairs]
    for t, s in rewiring.items():
        for col in BLOCK:
            rows[t][col] = pairs[s].features[col]
    return rows


# --- diagnostics --------------------------------------------------------------------------------------------


def alternatives(
    links: Sequence[Link], groups: Sequence[Hashable], distinct: Sequence[Hashable], width: float
) -> list[set[Hashable]]:
    """For each link, the distinct values (genes or elements) it could receive under the matching: links
    of its group whose distance differs by less than `width` and whose distinct key differs."""
    by_group: dict[Hashable, list[int]] = defaultdict(list)
    for k in range(len(links)):
        by_group[groups[k]].append(k)
    alts: list[set[Hashable]] = [set() for _ in links]
    for members in by_group.values():
        members.sort(key=lambda k: links[k].position)
        for a, i in enumerate(members):
            for j in members[a + 1 :]:
                if links[j].position - links[i].position >= width:
                    break
                if distinct[j] != distinct[i]:
                    alts[i].add(distinct[j])
                    alts[j].add(distinct[i])
    return alts


def mean(v: Sequence[float]) -> float:
    return sum(v) / len(v)


def variance(v: Sequence[float]) -> float:
    m = mean(v)
    return sum((x - m) ** 2 for x in v) / (len(v) - 1)


def smd(real: Sequence[float], attached: Sequence[float]) -> float | None:
    """(mean attached - mean real) over the pooled standard deviation sqrt((var real + var attached) / 2)."""
    if len(real) < 2 or len(attached) < 2:
        return None
    sd = math.sqrt((variance(real) + variance(attached)) / 2)
    diff = mean(attached) - mean(real)
    if sd == 0:
        return 0.0 if diff == 0 else math.copysign(math.inf, diff)
    return diff / sd


def ks(a: Sequence[float], b: Sequence[float]) -> float | None:
    """Two-sample Kolmogorov-Smirnov statistic: the largest gap between the two empirical CDFs."""
    if not a or not b:
        return None
    x, y = sorted(a), sorted(b)
    i = j = 0
    d = 0.0
    while i < len(x) and j < len(y):
        v = min(x[i], y[j])
        while i < len(x) and x[i] == v:
            i += 1
        while j < len(y) and y[j] == v:
            j += 1
        d = max(d, abs(i / len(x) - j / len(y)))
    return d


def quantiles(v: Sequence[float], qs: Sequence[float] = (0.1, 0.25, 0.5, 0.75, 0.9)) -> list[float] | None:
    s = sorted(v)
    if not s:
        return None
    out = []
    for q in qs:
        p = q * (len(s) - 1)
        lo = math.floor(p)
        out.append(round(s[lo] + (s[min(lo + 1, len(s) - 1)] - s[lo]) * (p - lo), 4))
    return out


def _r(x: float | None, nd: int = 4) -> float | None:
    return None if x is None or not math.isfinite(x) else round(x, nd)


def _share(a: int, b: int) -> float | None:
    return round(a / b, 4) if b else None


def diagnostics(
    links: Sequence[Link],
    groups: Sequence[Hashable],
    distinct: Sequence[Hashable],
    width: float,
    rewirings: Sequence[dict[int, int]],
    pairs_total: int,
) -> dict[str, Any]:
    """Everything the thresholds judge, for one null under one scheme, from the draws in `rewirings`.

    Labels split the counts into classes and are never used to choose a match. `pairs_total` is the size
    of the evaluated list, so the pairs outside the population (kept real in every rewiring) are counted.
    """
    by_index = {lk.index: lk for lk in links}
    pos_of = {lk.index: k for k, lk in enumerate(links)}
    gene_links = Counter(lk.gene for lk in links)
    n_all = len(links)
    n_reg = sum(lk.regulated for lk in links)

    alts = alternatives(links, groups, distinct, width)
    movable = [k for k in range(n_all) if alts[k]]
    movable_reg = [k for k in movable if links[k].regulated]
    sizes = Counter(min(len(alts[k]), 5) for k in movable)
    used = Counter(a for k in movable for a in alts[k])
    total_assignments = sum(len(alts[k]) for k in movable)
    top_decile = sorted(used.values(), reverse=True)[: max(1, math.ceil(len(used) / 10))] if used else []

    def two_plus(ks_: list[int]) -> float | None:
        return _share(sum(len(alts[k]) >= 2 for k in ks_), len(ks_))

    meaningful = {
        "pairs_evaluated": pairs_total,
        "pairs_outside_population": pairs_total - n_all,
        "population_links": {"all": n_all, "regulated": n_reg},
        "immovable_links": {"all": n_all - len(movable), "regulated": n_reg - len(movable_reg)},
        "movable_share": {"all": _share(len(movable), n_all), "regulated": _share(len(movable_reg), n_reg)},
        "distinct_alternatives_per_movable_link": {
            "1": sizes[1],
            "2": sizes[2],
            "3": sizes[3],
            "4": sizes[4],
            "5 or more": sizes[5],
            "median": quantiles([len(alts[k]) for k in movable], (0.5,))[0] if movable else None,
            "mean": _r(mean([len(alts[k]) for k in movable])) if movable else None,
        },
        "two_alternatives_share": {"all": two_plus(movable), "regulated": two_plus(movable_reg)},
        "links_with_only_one_or_two_alternatives": {
            "all": sum(len(alts[k]) <= 2 for k in movable),
            "regulated": sum(len(alts[k]) <= 2 for k in movable_reg),
        },
        "distinct_alternative_assignments_total": total_assignments,
        "distinct_alternative_values_total": len(used),
        "assignments_from_the_most_used_tenth_of_values": _share(sum(top_decile), total_assignments),
    }

    moved_all, moved_reg = [], []
    exact = {
        "permutation": True,
        "same_group": True,
        "distinct_differs": True,
        "within_width": True,
        "population_only": True,
        "gene_in_degree_preserved": True,
    }
    real: dict[str, dict[str, list[float]]] = {c: defaultdict(list) for c in CLASSES}
    att: dict[str, dict[str, list[float]]] = {c: defaultdict(list) for c in CLASSES}
    expr_seen = {c: {"moved": 0, "both": 0, "real_unknown": 0, "attached_unknown": 0} for c in CLASSES}
    sources_of: dict[int, Counter] = defaultdict(Counter)

    def cov(lk: Link, name: str) -> float | None:
        if name == "distance":
            return lk.position
        if name == "activity":
            return lk.activity
        if name == "expression":
            return lk.expression
        return math.log(gene_links[lk.gene])

    for rw in rewirings:
        if sorted(rw) != sorted(rw.values()):
            exact["permutation"] = False
        if not all(t in by_index and s in by_index for t, s in rw.items()):
            exact["population_only"] = False
            continue
        if Counter(by_index[s].gene for s in rw.values()) != Counter(by_index[t].gene for t in rw):
            exact["gene_in_degree_preserved"] = False
        moved_all.append(len(rw))
        moved_reg.append(sum(by_index[t].regulated for t in rw))
        for t, s in rw.items():
            a, b = by_index[t], by_index[s]
            if groups[pos_of[t]] != groups[pos_of[s]]:
                exact["same_group"] = False
            if distinct[pos_of[t]] == distinct[pos_of[s]]:
                exact["distinct_differs"] = False
            if abs(a.position - b.position) >= width:
                exact["within_width"] = False
            sources_of[t][s] += 1
            c = "regulated" if a.regulated else "not_regulated"
            for name in COVARIATES:
                if name == "expression":
                    # the same gene carries the same expression, known or not: held by identity
                    same = a.gene == b.gene
                    expr_seen[c]["moved"] += 1
                    expr_seen[c]["real_unknown"] += a.expression is None
                    expr_seen[c]["attached_unknown"] += b.expression is None
                    expr_seen[c]["both"] += same or (a.expression is not None and b.expression is not None)
                    if a.expression is None or b.expression is None:
                        continue
                real[c][name].append(cov(a, name))
                att[c][name].append(cov(b, name))

    draws = len(rewirings)
    realised = {
        "draws": draws,
        "moved_share_mean": {
            "all": _r(mean(moved_all) / n_all) if draws and n_all else None,
            "regulated": _r(mean(moved_reg) / n_reg) if draws and n_reg else None,
        },
        "moved_share_range": {
            "all": [_share(min(moved_all), n_all), _share(max(moved_all), n_all)] if draws else None,
            "regulated": [_share(min(moved_reg), n_reg), _share(max(moved_reg), n_reg)] if draws else None,
        },
        "moved_regulated_mean": _r(mean(moved_reg), 1) if draws else None,
        "modal_source_share_mean": {
            cls: _r(mean(vals)) if vals else None
            for cls, vals in (
                ("all", [max(c.values()) / sum(c.values()) for c in sources_of.values()]),
                (
                    "regulated",
                    [
                        max(c.values()) / sum(c.values())
                        for t, c in sources_of.items()
                        if by_index[t].regulated
                    ],
                ),
            )
        },
        "distinct_sources_realised_mean": {
            "all": _r(mean([len(c) for c in sources_of.values()])) if sources_of else None,
            "regulated": _r(mean(v))
            if (v := [len(c) for t, c in sources_of.items() if by_index[t].regulated])
            else None,
        },
    }

    balance: dict[str, Any] = {}
    for name in COVARIATES:
        balance[name] = {}
        for c in CLASSES:
            r, a = real[c][name], att[c][name]
            entry: dict[str, Any] = {"n": len(r), "smd": _r(smd(r, a))}
            if name in KS_COVARIATES:
                entry["ks"] = _r(ks(r, a))
                entry["real_quantiles"] = quantiles(r)
                entry["attached_quantiles"] = quantiles(a)
            balance[name][c] = entry
    for c in CLASSES:
        e = expr_seen[c]
        balance["expression"][c] |= {
            "known_share": _share(e["both"], e["moved"]),
            "real_unknown_share": _share(e["real_unknown"], e["moved"]),
            "attached_unknown_share": _share(e["attached_unknown"], e["moved"]),
        }
    return {
        "meaningful": meaningful | {"realised": realised},
        "balance": balance,
        "exact": exact,
        "coverage": "held by construction: links are covered pairs with a cached value and only exchange "
        "values with each other; pairs outside the population keep their real values",
        "links_per_element": "held by construction: no pair is added, dropped or moved between elements",
    }


def adequacy(diag: dict[str, Any], regulated_floor: int, thresholds: dict[str, float] = THRESHOLDS) -> dict:
    """Every registered check, its value, its bound and whether it is met; adequate only if all are."""
    checks: list[dict[str, Any]] = []

    def check(name: str, value: Any, bound: Any, kind: str) -> None:
        if kind == "min":
            met = value is not None and value >= bound
        elif kind == "max_abs":
            met = value is not None and abs(value) <= bound
        else:  # "true"
            met = value is True
        checks.append({"check": name, "value": value, "bound": bound, "kind": kind, "met": met})

    m, rz, b = diag["meaningful"], diag["meaningful"]["realised"], diag["balance"]
    for cls in ("all", "regulated"):
        check(f"movable share, {cls}", m["movable_share"][cls], thresholds["movable_share_min"], "min")
        check(
            f"moved share per rewiring, {cls}",
            rz["moved_share_mean"][cls],
            thresholds["moved_share_min"],
            "min",
        )
        check(
            f"movable links with two or more distinct alternatives, {cls}",
            m["two_alternatives_share"][cls],
            thresholds["two_alternatives_share_min"],
            "min",
        )
    check("regulated links moved per rewiring", rz["moved_regulated_mean"], regulated_floor, "min")
    for name in COVARIATES:
        for c in CLASSES:
            check(f"{name} |SMD|, {c}", b[name][c]["smd"], thresholds["smd_max"], "max_abs")
            if name in KS_COVARIATES:
                check(f"{name} KS, {c}", b[name][c]["ks"], thresholds["ks_max"], "max_abs")
    for c in CLASSES:
        e = b["expression"][c]
        check(f"expression known share, {c}", e["known_share"], thresholds["expression_known_min"], "min")
        gap = (
            None
            if e["attached_unknown_share"] is None
            else round(e["attached_unknown_share"] - e["real_unknown_share"], 4)
        )
        check(f"expression unknown-share gap, {c}", gap, thresholds["expression_unknown_gap_max"], "max_abs")
    for name, ok in diag["exact"].items():
        check(f"exact: {name}", ok, True, "true")
    return {
        "adequate": all(c["met"] for c in checks),
        "checks": checks,
        "failed": [c for c in checks if not c["met"]],
    }


def assess(
    pairs: Sequence[crispri.Pair],
    expression: dict[str, float],
    pair_set: str,
    null: str,
    scheme: dict[str, Any],
    draws: int = FEASIBILITY_DRAWS,
) -> dict[str, Any]:
    """The matching diagnostics of one null under one scheme on one pair set, and their verdict.

    No model is fitted and no metric is computed: the draws exist only to measure the matching.
    """
    links = links_of(pairs, expression)
    groups, distinct = keys_for(links, null, scheme["expression_classes"])
    width = scheme["width_log10"]
    rewirings = [
        rewire(links, groups, distinct, width, rng_for(pair_set, null, scheme["name"], d))
        for d in range(draws)
    ]
    diag = diagnostics(links, groups, distinct, width, rewirings, len(pairs))
    return {
        "pair_set": pair_set,
        "null": null,
        "scheme": scheme,
        "diagnostics": diag,
        "verdict": adequacy(diag, REGULATED_FLOOR[pair_set]),
    }


def choose(primary: Sequence[dict[str, Any]]) -> dict[str, Any] | None:
    """The first scheme, in SCHEMES order, whose primary assessment is adequate; None is the no-go."""
    for a in primary:
        if a["verdict"]["adequate"]:
            return a["scheme"]
    return None


# --- the expression covariate -------------------------------------------------------------------------------


def control_expression(h5ad: Path) -> dict[str, float]:
    """Mean raw pseudobulk over the non-targeting rows of Replogle et al.'s K562 file, per Ensembl gene ID
    with its version stripped: N1 amendment 2's strata variable (n1_perturb_response.raw_control_expression),
    read here from the cached file. Nothing but obs/gene_transcript, var/gene_id and those rows of X is
    read."""
    import h5py

    from genomeos.attribution import n1_perturb_response as n1

    with h5py.File(h5ad, "r") as h:
        obs = [x.decode() if isinstance(x, bytes) else str(x) for x in h["obs/gene_transcript"][()]]
        var = [x.decode() if isinstance(x, bytes) else str(x) for x in h["var/gene_id"][()]]
        rows = sorted(i for i, lab in enumerate(obs) if n1.parse_row(lab)["non_targeting"])
        x = h["X"][rows]
    means = x.astype("float64").mean(axis=0)
    return {gid.split(".")[0]: float(v) for gid, v in zip(var, means, strict=True)}


def expression_by_symbol(ensembl_of: dict[str, str], by_gene_id: dict[str, float]) -> dict[str, float]:
    """log1p expression per benchmark gene symbol; a symbol with no cached value is absent (unknown)."""
    return {
        sym: math.log1p(max(0.0, by_gene_id[gid])) for sym, gid in ensembl_of.items() if gid in by_gene_id
    }
