# SPDX-License-Identifier: AGPL-3.0-or-later
"""Is commitment, competence or lateral inhibition visible in the worm's measured factors?

Fate as a function of position, signals, state, genome, epigenome and time predicts signatures in a per-cell
factor atlas that a pure lookup does not need. Four tests over the Ma et al. 2021 atlas (266 factors, cells
named by lineage, 1.25 min per frame) and Sulston's lineage, each against a null:

1. states     cells born later occupy fewer, more separated factor states: per birth-time bin (equal-size
              subsamples), the mean pairwise Jaccard distance, the nearest-neighbour distance relative to it
              (discreteness), the silhouette of cells grouped by the tissue their descendants mostly make
              (fate separation), and the discreteness in excess of a curveball randomisation that keeps
              every cell's factor count and every factor's frequency in the bin. Trend: least-squares slope
              over bins; null: birth times shuffled among cells.
2. programmes competing tissue programmes become mutually exclusive: among cells with a textbook programme
              factor on (intestine, body-wall muscle, hypodermis, neuron; PHA-4 and CEH-22 mark an organ whose
              cells are also muscle, neuron or epithelium, so the pharynx is not a competing programme), the
              share with two or more programmes, observed and relative to the curveball expectation.
   ratchet    on mother -> daughter pairs: a programme factor on in the mother stays on (persistence), and a
              daughter of a cell with a programme gains a competing one (switching); by the mother's birth.
              Reporter protein perdures across a division, so persistence is an upper bound; switching is not
              inflated by perdurance.
3. competence the same factor switching on at a different time leads to a different fate: for every de novo
              onset (factor present, mother without it), the mutual information between onset time bin and
              the dominant descendant tissue, given the factor and the founder sublineage; null: onset times
              permuted within factor x sublineage.
4. sisters    equivalent sisters diverge: sister pairs by division axis (a/p sisters differ by POP-1/Wnt, l/r
              sisters are the equivalent case), factor distance against the measurement floor of the same
              factor seen by two reporter strains in the same cell, fate discordance, and the Notch target
              REF-1 (Neves & Priess 2005, Dev Cell 8:867) discordant between sisters against factors of the
              same frequency.

Cells enter only with a complete lifetime (their last frame before the end of imaging), because a cell cut off
by the end of imaging has fewer frames in which a factor can be seen.
"""

from __future__ import annotations

import math
import random
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass, field

from .fate_rules import sublineage
from .reference import ReferenceLineage

BIN_EDGES = (
    0,
    40,
    80,
    120,
    160,
    200,
)  # birth frame (1.25 min each); later births rarely finish their lives in view
END_OF_IMAGING_FRAME = 318
PROGRAMMES = {
    "intestine": ("ELT-2", "ELT-7"),
    "muscle": ("HLH-1", "UNC-120"),
    "hypoderm": ("ELT-1", "LIN-26", "NHR-25"),
    "neuron": ("LIN-32", "NGN-1", "HLH-14", "CND-1"),
}
PROGRAMME_SOURCES = {
    "intestine": "Fukushige et al. 1998; Sommermann et al. 2010",
    "muscle": "Krause et al. 1990; Fukushige et al. 2006",
    "hypoderm": "Page et al. 1997; Labouesse et al. 1994; Gissendanner & Sluder 2000",
    "neuron": "Zhao & Emmons 1995 (LIN-32); Nakano et al. 2010 (NGN-1); Frank et al. 2003 (HLH-14); "
    "Hallam et al. 2000 (CND-1)",
}
NOTCH_TARGETS = ("REF-1",)


@dataclass(slots=True)
class AtlasCell:
    id: str
    born: int  # first frame
    last: int
    factors: frozenset[str]
    bits: int
    sublineage: str
    tissues: Counter = field(default_factory=Counter)  # surviving terminal descendants by tissue

    @property
    def dominant(self) -> str:
        return self.tissues.most_common(1)[0][0] if self.tissues else ""


def _descendant_tissues(ref: ReferenceLineage, cid: str) -> Counter:
    out: Counter = Counter()
    stack = [cid]
    while stack:
        c = ref.cells[stack.pop()]
        if c.children:
            stack.extend(c.children)
        elif c.dies is None and c.tissue:
            out[c.tissue] += 1
    return out


def atlas_cells(
    ref: ReferenceLineage, atlas: dict[str, list[str]], lifetimes: dict[str, list[int]], complete: bool = True
) -> dict[str, AtlasCell]:
    names = sorted({f for fs in atlas.values() for f in fs})
    index = {f: i for i, f in enumerate(names)}
    out = {}
    for cid, fs in atlas.items():
        if cid not in ref.cells or cid not in lifetimes:
            continue
        born, last = lifetimes[cid]
        if complete and (last >= END_OF_IMAGING_FRAME or born <= 0 or born >= BIN_EDGES[-1]):
            continue
        bits = 0
        for f in fs:
            bits |= 1 << index[f]
        out[cid] = AtlasCell(
            cid, born, last, frozenset(fs), bits, sublineage(cid), _descendant_tissues(ref, cid)
        )
    return out


def time_bin(frame: int, edges: tuple[int, ...] = BIN_EDGES) -> int:
    for i in range(len(edges) - 1):
        if edges[i] <= frame < edges[i + 1]:
            return i
    return len(edges) - 2 if frame >= edges[-1] else 0


def _jaccard(a: int, b: int) -> float:
    u = (a | b).bit_count()
    return 1.0 - (a & b).bit_count() / u if u else 0.0


def _slope(ys: list[float]) -> float:
    xs = list(range(len(ys)))
    pts = [(x, y) for x, y in zip(xs, ys, strict=True) if y is not None and not math.isnan(y)]
    if len(pts) < 2:
        return 0.0
    mx = statistics.fmean(x for x, _ in pts)
    my = statistics.fmean(y for _, y in pts)
    sxx = sum((x - mx) ** 2 for x, _ in pts)
    return sum((x - mx) * (y - my) for x, y in pts) / sxx if sxx else 0.0


def curveball(rows: list[int], rng: random.Random, trades: int | None = None) -> list[int]:
    """Strona et al. 2014 curveball on a presence matrix given as row bitsets: keeps every row's count and
    every column's count, randomises which rows share which columns."""
    rows = list(rows)
    n = len(rows)
    if n < 2:
        return rows
    for _ in range(trades if trades is not None else 5 * n):
        i, j = rng.sample(range(n), 2)
        a, b = rows[i], rows[j]
        only_a, only_b = a & ~b, b & ~a
        pool = only_a | only_b
        if not only_a or not only_b:
            continue
        bits = [k for k in range(pool.bit_length()) if pool >> k & 1]
        rng.shuffle(bits)
        na = only_a.bit_count()
        new_a = 0
        for k in bits[:na]:
            new_a |= 1 << k
        new_b = pool & ~new_a
        common = a & b
        rows[i], rows[j] = common | new_a, common | new_b
    return rows


# ---- 1. states ------------------------------------------------------------------


def _silhouette(d: list[list[float]], labels: list[str]) -> float:
    n = len(labels)
    sil = []
    groups = Counter(labels)
    for i in range(n):
        if not labels[i] or groups[labels[i]] < 2 or len(groups) < 2:
            continue
        own = [d[i][j] for j in range(n) if j != i and labels[j] == labels[i]]
        others = defaultdict(list)
        for j in range(n):
            if labels[j] and labels[j] != labels[i]:
                others[labels[j]].append(d[i][j])
        if not others:
            continue
        a = statistics.fmean(own)
        b = min(statistics.fmean(v) for v in others.values())
        sil.append((b - a) / max(a, b) if max(a, b) > 0 else 0.0)
    return statistics.fmean(sil) if sil else float("nan")


def _state_stats(bits: list[int], labels: list[str], lineages: list[str] | None = None) -> dict[str, float]:
    n = len(bits)
    d = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d[i][j] = d[j][i] = _jaccard(bits[i], bits[j])
    pair = [d[i][j] for i in range(n) for j in range(i + 1, n)]
    spread = statistics.fmean(pair) if pair else 0.0
    nn = statistics.fmean(min(d[i][j] for j in range(n) if j != i) for i in range(n)) if n > 1 else 0.0
    return {
        "factors_per_cell": statistics.fmean(b.bit_count() for b in bits),
        "spread": spread,
        "nearest_neighbour": nn,
        "discreteness": nn / spread if spread else float("nan"),
        "fate_silhouette": _silhouette(d, labels),
        "lineage_silhouette": _silhouette(d, lineages) if lineages else float("nan"),
    }


def _binned(cells: list[AtlasCell], births: list[int]) -> list[list[AtlasCell]]:
    bins: list[list[AtlasCell]] = [[] for _ in range(len(BIN_EDGES) - 1)]
    for c, b in zip(cells, births, strict=True):
        bins[time_bin(b)].append(c)
    return bins


def _bin_series(
    bins: list[list[AtlasCell]], rng: random.Random, n_sub: int, repeats: int, excess: bool
) -> dict:
    keys = (
        "factors_per_cell",
        "spread",
        "nearest_neighbour",
        "discreteness",
        "fate_silhouette",
        "lineage_silhouette",
    )
    series: dict[str, list[float]] = {k: [] for k in keys}
    if excess:
        series["discreteness_excess"] = []
    for cs in bins:
        acc: dict[str, list[float]] = defaultdict(list)
        for _ in range(repeats):
            sub = rng.sample(cs, min(n_sub, len(cs)))
            st = _state_stats([c.bits for c in sub], [c.dominant for c in sub], [c.sublineage for c in sub])
            for k, v in st.items():
                acc[k].append(v)
            if excess:
                rnd = curveball([c.bits for c in cs], rng)
                idx = rng.sample(range(len(cs)), min(n_sub, len(cs)))
                null = _state_stats([rnd[i] for i in idx], [""] * len(idx))
                acc["discreteness_excess"].append(st["discreteness"] - null["discreteness"])
        for k in series:
            vals = [v for v in acc[k] if not math.isnan(v)]
            series[k].append(statistics.fmean(vals) if vals else float("nan"))
    return series


def states_test(
    cells: dict[str, AtlasCell], n_sub: int = 36, repeats: int = 12, permutations: int = 100, seed: int = 0
) -> dict:
    rng = random.Random(seed)
    cs = list(cells.values())
    births = [c.born for c in cs]
    observed = _bin_series(_binned(cs, births), rng, n_sub, repeats, excess=True)
    slopes = {k: _slope(v) for k, v in observed.items()}
    null: dict[str, list[float]] = defaultdict(list)
    for _ in range(permutations):
        shuffled = births[:]
        rng.shuffle(shuffled)
        s = _bin_series(_binned(cs, shuffled), rng, n_sub, max(3, repeats // 3), excess=True)
        for k, v in s.items():
            null[k].append(_slope(v))
    out = {
        "bins_birth_frame": [f"{BIN_EDGES[i]}-{BIN_EDGES[i + 1]}" for i in range(len(BIN_EDGES) - 1)],
        "cells_per_bin": [len(b) for b in _binned(cs, births)],
        "subsample": n_sub,
        "series": {k: [round(x, 4) for x in v] for k, v in observed.items()},
        "slope": {k: round(v, 5) for k, v in slopes.items()},
        "p_shuffled_time": {},
        "permutations": permutations,
    }
    for k, v in slopes.items():
        if k in null:
            out["p_shuffled_time"][k] = {
                "lower": round((1 + sum(1 for x in null[k] if x <= v)) / (1 + len(null[k])), 4),
                "higher": round((1 + sum(1 for x in null[k] if x >= v)) / (1 + len(null[k])), 4),
            }
    return out


# ---- 2. programmes and the ratchet ------------------------------------------------


def programmes_on(factors: frozenset[str] | set[str]) -> set[str]:
    return {p for p, fs in PROGRAMMES.items() if any(f in factors for f in fs)}


def _coexpression(sets: list[frozenset[str]]) -> float:
    with_one = [programmes_on(s) for s in sets]
    with_one = [p for p in with_one if p]
    return sum(1 for p in with_one if len(p) >= 2) / len(with_one) if with_one else float("nan")


def _unbits(rows: list[int], names: list[str]) -> list[frozenset[str]]:
    return [frozenset(names[k] for k in range(r.bit_length()) if r >> k & 1) for r in rows]


def programme_test(
    cells: dict[str, AtlasCell], atlas_names: list[str], permutations: int = 200, seed: int = 0
) -> dict:
    rng = random.Random(seed)
    cs = list(cells.values())
    births = [c.born for c in cs]

    def series(bins: list[list[AtlasCell]], expected: bool) -> tuple[list[float], list[float]]:
        obs, ratio = [], []
        for b in bins:
            share = _coexpression([c.factors for c in b])
            obs.append(share)
            if expected:
                exp = statistics.fmean(
                    _coexpression(_unbits(curveball([c.bits for c in b], rng), atlas_names))
                    for _ in range(10)
                )
                ratio.append(share / exp if exp and not math.isnan(share) else float("nan"))
        return obs, ratio

    obs, ratio = series(_binned(cs, births), expected=True)
    null_obs = []
    for _ in range(permutations):
        sh = births[:]
        rng.shuffle(sh)
        null_obs.append(_slope(series(_binned(cs, sh), expected=False)[0]))
    s = _slope(obs)
    return {
        "programmes": {p: list(fs) for p, fs in PROGRAMMES.items()},
        "sources": PROGRAMME_SOURCES,
        "cells_with_a_programme_per_bin": [
            sum(1 for c in b if programmes_on(c.factors)) for b in _binned(cs, births)
        ],
        "share_with_two_or_more": [round(x, 3) for x in obs],
        "relative_to_curveball": [round(x, 3) for x in ratio],
        "slope": round(s, 5),
        "slope_relative": round(_slope(ratio), 5),
        "p_shuffled_time_lower": round((1 + sum(1 for x in null_obs if x <= s)) / (1 + len(null_obs)), 4),
    }


def _pearson(xs: list[float], ys: list[float]) -> float:
    if len(xs) < 3:
        return float("nan")
    mx, my = statistics.fmean(xs), statistics.fmean(ys)
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    sy = math.sqrt(sum((y - my) ** 2 for y in ys))
    return (
        sum((x - mx) * (y - my) for x, y in zip(xs, ys, strict=True)) / (sx * sy)
        if sx and sy
        else float("nan")
    )


def ratchet_test(
    cells: dict[str, AtlasCell], ref: ReferenceLineage, permutations: int = 1000, seed: int = 0
) -> dict:
    """Per mother -> daughter pair: persistence (programme factors on in the mother that stay on in the
    daughter) and switching (a daughter of a programme-carrying mother gains a competing programme), each
    correlated with the mother's birth frame; null: birth frames shuffled among pairs."""
    rng = random.Random(seed)
    persist: list[tuple[int, float]] = []
    switch: list[tuple[int, float]] = []
    for c in cells.values():
        for kid in ref.cells[c.id].children:
            d = cells.get(kid)
            if d is None:
                continue
            mother_f = [f for fs in PROGRAMMES.values() for f in fs if f in c.factors]
            if mother_f:
                persist.append((c.born, sum(1 for f in mother_f if f in d.factors) / len(mother_f)))
            mp = programmes_on(c.factors)
            if mp:
                switch.append((c.born, float(bool(programmes_on(d.factors) - mp))))

    def test(rows: list[tuple[int, float]], sign: int) -> dict:
        if len(rows) < 5:
            return {"pairs": len(rows)}
        xs, ys = [float(b) for b, _ in rows], [y for _, y in rows]
        r = _pearson(xs, ys)
        null = []
        for _ in range(permutations):
            sh = xs[:]
            rng.shuffle(sh)
            null.append(_pearson(sh, ys))
        bins = defaultdict(list)
        for b, y in rows:
            bins[time_bin(b)].append(y)
        return {
            "pairs": len(rows),
            "mean": round(statistics.fmean(ys), 3),
            "by_bin": {str(k): [round(statistics.fmean(v), 3), len(v)] for k, v in sorted(bins.items())},
            "r_with_mother_birth": round(r, 4),
            "p": round((1 + sum(1 for x in null if sign * x >= sign * r)) / (1 + permutations), 4),
        }

    return {
        "persistence_rises": test(persist, +1),
        "switching_falls": test(switch, -1),
        "caveat": "reporter protein perdures across a division, so persistence is an upper bound",
    }


# ---- 3. competence ------------------------------------------------------------------


def _soft_cmi(rows: list[tuple[str, str, dict[str, float]]]) -> float:
    """I(time; fate | stratum) in bits from (stratum, time, fate weights) rows; a row's weights sum to 1, so a
    cell whose descendants make several tissues contributes to each in proportion."""
    by: dict[str, list[tuple[str, dict[str, float]]]] = defaultdict(list)
    for s, t, w in rows:
        by[s].append((t, w))
    n = len(rows)
    total = 0.0
    for v in by.values():
        joint: dict[tuple[str, str], float] = defaultdict(float)
        a: dict[str, float] = defaultdict(float)
        b: dict[str, float] = defaultdict(float)
        for t, w in v:
            for f, x in w.items():
                joint[(t, f)] += x
                a[t] += x
                b[f] += x
        m = len(v)
        mi = sum(k / m * math.log2(k * m / (a[t] * b[f])) for (t, f), k in joint.items() if k > 0)
        total += m / n * mi
    return total if n else 0.0


def competence_test(
    cells: dict[str, AtlasCell],
    ref: ReferenceLineage,
    factors: set[str] | None = None,
    permutations: int = 500,
    min_onsets: int = 4,
    by_generation: bool = False,
    seed: int = 0,
) -> dict:
    """De novo onsets of a factor (present in the cell, absent in its tracked mother), stratified by factor
    and founder sublineage; does the onset time bin carry information about fate beyond the stratum? Two fate
    reads: the dominant descendant tissue, and the full descendant composition. The first rises with time
    for a trivial reason (a later onset is a smaller clone, so its majority is purer and can be a minority
    tissue of the sublineage); the composition read is not biased that way, so it is the test. With
    `by_generation` the stratum also fixes the number of divisions since the founder, so time can only differ
    through cycle speed between branches of the same depth."""
    rng = random.Random(seed)
    rows: list[tuple[str, str, dict[str, float], str]] = []
    for c in cells.values():
        parent = ref.cells[c.id].parent
        p = cells.get(parent) if parent else None
        if p is None or not c.tissues:
            continue
        n = sum(c.tissues.values())
        comp = {t: k / n for t, k in c.tissues.items()}
        for f in c.factors - p.factors:
            if factors is None or f in factors:
                stratum = f"{f}|{c.sublineage}" + (f"|{ref.cells[c.id].generation}" if by_generation else "")
                rows.append((stratum, str(time_bin(c.born)), comp, c.dominant))
    strata = Counter(r[0] for r in rows)
    spans = defaultdict(set)
    for r in rows:
        spans[r[0]].add(r[1])
    rows = [r for r in rows if strata[r[0]] >= min_onsets and len(spans[r[0]]) > 1]
    by = defaultdict(list)
    for i, r in enumerate(rows):
        by[r[0]].append(i)

    def stats(times: list[str]) -> tuple[float, float]:
        soft = _soft_cmi([(r[0], t, r[2]) for r, t in zip(rows, times, strict=True)])
        hard = _soft_cmi([(r[0], t, {r[3]: 1.0}) for r, t in zip(rows, times, strict=True)])
        return soft, hard

    observed = [r[1] for r in rows]
    obs_soft, obs_hard = stats(observed)
    null_soft, null_hard = [], []
    for _ in range(permutations):
        times = observed[:]
        for idx in by.values():
            ts = [times[i] for i in idx]
            rng.shuffle(ts)
            for i, t in zip(idx, ts, strict=True):
                times[i] = t
        a, b = stats(times)
        null_soft.append(a)
        null_hard.append(b)

    def summary(obs: float, null: list[float]) -> dict:
        if not null:
            return {"bits": round(obs, 4)}
        med = statistics.median(null)
        return {
            "bits": round(obs, 4),
            "null_median_bits": round(med, 4),
            "excess_bits": round(obs - med, 4),
            "p": round((1 + sum(1 for x in null if x >= obs)) / (1 + len(null)), 4),
        }

    return {
        "onsets": len(rows),
        "strata": len(by),
        "stratum": "factor x sublineage" + (" x generation" if by_generation else ""),
        "composition": summary(obs_soft, null_soft),
        "dominant_tissue": summary(obs_hard, null_hard),
        "permutations": permutations,
    }


# ---- 4. sisters and equivalent cells ---------------------------------------------------


def _axis(a: str, b: str) -> str:
    return {frozenset("ap"): "a/p", frozenset("lr"): "l/r", frozenset("dv"): "d/v"}.get(
        frozenset(a[-1] + b[-1]), "other"
    )


def _tv(a: Counter, b: Counter) -> float:
    na, nb = sum(a.values()), sum(b.values())
    if not na or not nb:
        return float("nan")
    return 0.5 * sum(abs(a[k] / na - b[k] / nb) for k in set(a) | set(b))


def _pair_summary(pairs: list[tuple[AtlasCell, AtlasCell]]) -> dict:
    dist = [_jaccard(a.bits, b.bits) for a, b in pairs]
    tv = [x for x in (_tv(a.tissues, b.tissues) for a, b in pairs) if not math.isnan(x)]
    differ = [(a, b) for a, b in pairs if a.dominant and b.dominant and a.dominant != b.dominant]
    agree = [(a, b) for a, b in pairs if a.dominant and a.dominant == b.dominant]

    def med(ps: list[tuple[AtlasCell, AtlasCell]]) -> float | None:
        return round(statistics.median(_jaccard(a.bits, b.bits) for a, b in ps), 3) if ps else None

    return {
        "pairs": len(pairs),
        "factor_distance_median": round(statistics.median(dist), 3) if dist else None,
        "fate_total_variation_median": round(statistics.median(tv), 3) if tv else None,
        "dominant_tissue_differs": len(differ),
        "factor_distance_when_fates_differ": med(differ),
        "factor_distance_when_fates_agree": med(agree),
        "examples_fates_differ": [f"{a.id}:{a.dominant} / {b.id}:{b.dominant}" for a, b in differ[:8]],
    }


def sister_test(
    cells: dict[str, AtlasCell],
    ref: ReferenceLineage,
    strain_peaks: dict[str, dict[str, list[float]]],
    levels: dict[str, dict[str, list]],
    fraction: float = 0.2,
    permutations: int = 2000,
    seed: int = 0,
) -> dict:
    rng = random.Random(seed)
    sisters: list[tuple[AtlasCell, AtlasCell, str]] = []
    for c in cells.values():
        kids = ref.cells[c.id].children
        if len(kids) == 2 and all(k in cells for k in kids):
            a, b = cells[kids[0]], cells[kids[1]]
            sisters.append((a, b, _axis(a.id, b.id)))
    by_axis = {
        axis: _pair_summary([(a, b) for a, b, x in sisters if x == axis]) for axis in ("a/p", "l/r", "d/v")
    }
    # bilateral homologues: ABpl... and ABpr... with the same suffix are mirror images (Sulston 1983)
    homologues = [
        (c, cells["ABpr" + c.id[4:]])
        for c in cells.values()
        if c.id.startswith("ABpl") and "ABpr" + c.id[4:] in cells
    ]
    # the measurement floor: a factor two reporter strains follow, called present (each at >= `fraction` of
    # its own strain's maximum) in one and not the other, in the same cell
    floor_disagree = floor_n = 0
    for per_cell in strain_peaks.values():
        for v in per_cell.values():
            calls = [x >= fraction for x in v]
            if any(calls):
                floor_n += 1
                floor_disagree += not all(calls)

    def disagreement(pairs: list[tuple[AtlasCell, AtlasCell]]) -> float | None:
        k = n = 0
        for a, b in pairs:
            for f in strain_peaks:
                ia, ib = f in a.factors, f in b.factors
                if ia or ib:
                    n += 1
                    k += ia != ib
        return round(k / n, 3) if n else None

    # POP-1: anterior (or left, dorsal) sister over the other, mean level over each life
    pop = levels.get("POP-1", {})
    ap = []
    for a, b, axis in sisters:
        if axis != "a/p":
            continue
        first, second = (a, b) if a.id.endswith("a") else (b, a)
        x, y = pop.get(first.id, [0, 0])[1], pop.get(second.id, [0, 0])[1]
        if x > 0 and y > 0:
            ap.append(
                (
                    math.log2(x / y),
                    first.dominant != second.dominant and bool(first.dominant and second.dominant),
                )
            )
    higher = sum(1 for r, _ in ap if r > 0)
    # two-sided sign test by normal approximation
    z = (higher - len(ap) / 2) / math.sqrt(len(ap) / 4) if ap else 0.0
    differ = [abs(r) for r, d in ap if d]
    agree = [abs(r) for r, d in ap if not d]
    gap = statistics.fmean(differ) - statistics.fmean(agree) if differ and agree else float("nan")
    null = []
    flags = [d for _, d in ap]
    mags = [abs(r) for r, _ in ap]
    for _ in range(permutations if differ and agree else 0):
        rng.shuffle(flags)
        dd = [m for m, f in zip(mags, flags, strict=True) if f]
        aa = [m for m, f in zip(mags, flags, strict=True) if not f]
        null.append(statistics.fmean(dd) - statistics.fmean(aa))
    # REF-1 discordance between sisters and fate discordance, against factors of similar frequency
    freq = Counter(f for c in cells.values() for f in c.factors)
    typed = [(a, b) for a, b, _ in sisters if a.dominant and b.dominant]

    def association(f: str) -> float | None:
        disc = [a.dominant != b.dominant for a, b in typed if (f in a.factors) != (f in b.factors)]
        conc = [a.dominant != b.dominant for a, b in typed if (f in a.factors) == (f in b.factors)]
        if len(disc) < 3 or not conc:
            return None
        return sum(disc) / len(disc) - sum(conc) / len(conc)

    notch = {}
    for f in NOTCH_TARGETS:
        if f not in freq:
            notch[f] = {"status": "not in atlas"}
            continue
        obs = association(f)
        similar = [g for g in freq if g != f and 0.5 * freq[f] <= freq[g] <= 2 * freq[f]]
        others = [x for x in (association(g) for g in similar) if x is not None]
        notch[f] = {
            "cells": freq[f],
            "fate_discordance_gain_when_sisters_differ_on_it": None if obs is None else round(obs, 3),
            "same_frequency_factors": len(others),
            "rank_among_them": None if obs is None else 1 + sum(1 for x in others if x > obs),
        }
    return {
        "sister_pairs": len(sisters),
        "by_axis": by_axis,
        "bilateral_homologues_ABpl_ABpr": _pair_summary(homologues),
        "measurement_floor": {
            "factors_with_two_strains": len(strain_peaks),
            "strain_disagreement": round(floor_disagree / floor_n, 3) if floor_n else None,
            "sister_disagreement_same_factors": disagreement([(a, b) for a, b, _ in sisters]),
            "homologue_disagreement_same_factors": disagreement(homologues),
        },
        "pop1_anterior_over_posterior": {
            "pairs": len(ap),
            "anterior_higher": higher,
            "log2_ratio_median": round(statistics.median(r for r, _ in ap), 3) if ap else None,
            "sign_test_z": round(z, 2),
            "abs_log2_ratio_fates_differ_minus_agree": None if math.isnan(gap) else round(gap, 3),
            "p_asymmetry_larger_when_fates_differ": round(
                (1 + sum(1 for x in null if x >= gap)) / (1 + len(null)), 4
            )
            if null
            else None,
        },
        "notch_targets": notch,
    }
