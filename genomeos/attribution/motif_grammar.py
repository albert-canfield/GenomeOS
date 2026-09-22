# SPDX-License-Identifier: AGPL-3.0-or-later
"""Does the *arrangement* of motif sites predict measured enhancer activity better than their presence?

Area J has kept saying "grammar" without a test that could fail. Section 8 read a promoter as a
dependency list of factors; §15 asked whether sites held across species separate enhancers from
conserved negatives (they do not); §16 measured which bases matter and found that JASPAR sites at the
usual 0.85 threshold cover 99.9% of a regulatory element's bases, so a loose site call says nothing,
while a strict one (0.95) is selective and modestly informative. What nobody in this area has tested
is the claim the DNA-as-code reading rests on: that *word order* matters -- that how far apart two
factor sites sit, on which strand, and in which helical phase, carries information beyond which
factors are present and how often.

ENCODE4's joint lentiMPRA library (attribution/mpra.py) makes that testable: 51,376 elements of 200 bp,
each with a measured log2(RNA/DNA) in K562, HepG2 and WTC11. This module builds three nested feature
sets per element and compares them out of sample:

  (a) composition only: GC, GC squared, CpG observed over expected;
  (b) plus strict motif counts: log1p sites per TFClass family unit (JASPAR 2026 CORE vertebrates at
      relative score 0.95, the threshold measurement calibrated in §16), the COUNT_UNITS most frequent
      family units, plus total sites and distinct families;
  (c) plus grammar: for every unordered pair of the PAIR_UNITS most frequent family units (including a
      family with itself), the number of non-overlapping site pairs in each gap bin (GAP_BINS), the
      number on the same strand, and the number whose centre-to-centre distance is in helical phase
      (a multiple of 10.5 bp, within a quarter turn).

A fourth model is the falsifier: the same grammar features computed after the site labels and strands
are permuted *within* each element, which keeps every count and every position but destroys which
family sits where and how it is turned. If (c) beats (b) only because pair counts are a nonlinear
function of counts, the shuffled control beats (b) by as much.

Ridge regression (dependency-free, epigenome.py's solve is the same idea with a Cholesky factor here so
one factorisation serves three cell lines) is fitted on the training chromosomes; the penalty is chosen
by fold-by-chromosome cross-validation inside the training set only. Held-out chromosomes (HELD_OUT)
are scored once, after PREREGISTERED was committed. Intervals are bootstrap over elements.

Evidence: activity is `experimental` (ENCODE4 lentiMPRA); sites are `predicted` (a matrix score, not a
binding event); the comparison between feature sets is `inferred`.
"""

from __future__ import annotations

import json
import math
import random
from array import array
from collections import defaultdict
from pathlib import Path
from typing import Any

from genomeos.results import RESULTS_DIR, save_result

PREREGISTERED = (
    "(c) beats (b) in Spearman with measured activity on held-out chromosomes in each cell line: "
    "the 95% bootstrap interval of rho(c) - rho(b) over held-out elements lies above 0 in K562, HepG2 "
    "and WTC11. The arrangement reading additionally requires rho(c) - rho(c shuffled) above 0 by the "
    "same interval; if only (c) - (b) holds, the gain is pair counts, not arrangement."
)
SITE_THRESHOLD = 0.95  # relative matrix score; below it JASPAR sites cover 97%+ of bases (see §16)
HELD_OUT = ("chr8", "chr9", "chr21", "chr22")  # scored once; chr1-7, 10-20, X, Y train
COUNT_UNITS = 40  # family units with count features, the most frequent on the training chromosomes
PAIR_UNITS = 10  # family units entering pair (grammar) features: PAIR_UNITS * (PAIR_UNITS + 1) / 2 pairs
GAP_BINS = ((0, 10), (10, 20), (20, 50), (50, None))  # edge-to-edge gap in bp between two sites
HELICAL_PERIOD = 10.5  # bp per turn of B-form DNA
PHASE_TOLERANCE = HELICAL_PERIOD / 4  # centre distance within a quarter turn of a whole turn
LAMBDAS = (1.0, 10.0, 100.0, 1000.0, 10000.0)  # ridge penalty in standardised units
INNER_FOLDS = 4  # folds by chromosome inside the training set, for the penalty only
BOOTSTRAP = 1000
SEED = 20260916
CELLS = ("K562", "HepG2", "WTC11")
KNOWLEDGE = Path("data/knowledge/mpra")
CHROMS = tuple(f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"])
EVIDENCE = {
    "activity": "experimental: ENCODE4 lentiMPRA (ENCSR106SZM), log2(RNA/DNA) per 200 bp element",
    "sites": (
        f"predicted: JASPAR 2026 CORE vertebrates at relative score >= {SITE_THRESHOLD}, "
        "collapsed to one site per TFClass family unit per position"
    ),
    "comparison": (
        "inferred: ridge on the training chromosomes, Spearman on the held-out ones, bootstrap by element"
    ),
}


# ---------------------------------------------------------------- sites


def collapse(hits: list[tuple[str, int, int, int, float]]) -> list[tuple[str, int, int, int]]:
    """One site per family unit per position: the best-scoring hit, then nothing overlapping it.

    Matrices of one family bind the same sequence, so a GC box read by KLF1, KLF5 and SP2 is one site,
    not three; overlapping hits of the same unit on either strand are the same site seen twice.
    """
    out: list[tuple[str, int, int, int]] = []
    by_unit: dict[str, list[tuple[str, int, int, int, float]]] = defaultdict(list)
    for h in hits:
        by_unit[h[0]].append(h)
    for unit, group in by_unit.items():
        kept: list[tuple[int, int, int]] = []
        for _u, start, end, strand, _rel in sorted(group, key=lambda h: (-h[4], h[1], h[3])):
            if any(start < e and s < end for s, e, _ in kept):
                continue
            kept.append((start, end, strand))
        out.extend((unit, s, e, sd) for s, e, sd in kept)
    return sorted(out, key=lambda s: (s[1], s[2], s[0]))


def scan_sites(motifs: list, seqs: list[str], families: dict[str, dict[str, str]]) -> list[list[tuple]]:
    """Every strict site of every family unit on each sequence, as (unit, start, end, strand).

    The motifs must be prepared at SITE_THRESHOLD. The k-mer index of genome/motifs.py is built once
    over the whole batch and every feasible core of every profile is looked up in it.
    """
    from genomeos.genome.motifs import MAX_POS, Index, family_unit

    idx = Index(seqs, {m.core_k for m in motifs})
    raw: list[list[tuple[str, int, int, int, float]]] = [[] for _ in seqs]
    for m in motifs:
        unit = family_unit(m.name, families)
        table = idx.tables[m.core_k]
        w, span = m.width, m.maximum - m.minimum
        for code in m.feasible_cores():
            for packed in table.get(code, ()):
                i, strand, pos = packed >> 12, (packed >> 11) & 1, packed & (MAX_POS - 1)
                start = pos - m.core_offset
                codes = idx.rev[i] if strand else idx.fwd[i]
                if start < 0 or start + w > len(codes):
                    continue
                s = m.score(codes, start)
                if s is None or s < m.threshold:
                    continue
                rel = (s - m.minimum) / span if span else 1.0
                fwd = len(codes) - start - w if strand else start
                raw[i].append((unit, fwd, fwd + w, strand, rel))
    return [collapse(h) for h in raw]


def shuffle_sites(sites: list[tuple], rng: random.Random) -> list[tuple]:
    """The same site positions with their family labels and strands permuted within the element.

    Every count is preserved exactly (the multiset of units and of positions is unchanged); what is
    destroyed is which family sits at which position and how each site is turned, which is precisely
    what the grammar features read.
    """
    units = [s[0] for s in sites]
    strands = [s[3] for s in sites]
    rng.shuffle(units)
    rng.shuffle(strands)
    return sorted(
        ((u, s[1], s[2], sd) for u, s, sd in zip(units, sites, strands, strict=True)),
        key=lambda s: (s[1], s[2], s[0]),
    )


# ---------------------------------------------------------------- features


def composition(seq: str) -> list[float]:
    """GC, GC squared and CpG observed over expected: the feature set every model starts from."""
    s = seq.upper()
    c, g = s.count("C"), s.count("G")
    n = c + g + s.count("A") + s.count("T")
    if not n:
        return [0.0, 0.0, 0.0]
    gc = (c + g) / n
    oe = (s.count("CG") * n / (c * g)) if c and g else 0.0
    return [gc, gc * gc, oe]


def count_features(sites: list[tuple], units: list[str]) -> dict[int, float]:
    """log1p sites per listed family unit, then total sites and distinct families (columns 0..len+1)."""
    n: dict[str, int] = defaultdict(int)
    for s in sites:
        n[s[0]] += 1
    out = {i: math.log1p(n[u]) for i, u in enumerate(units) if n[u]}
    if sites:
        out[len(units)] = math.log1p(len(sites))
        out[len(units) + 1] = math.log1p(len(n))
    return out


def gap_bin(gap: int) -> int | None:
    for i, (lo, hi) in enumerate(GAP_BINS):
        if gap >= lo and (hi is None or gap < hi):
            return i
    return None


def in_phase(centre_distance: float) -> bool:
    r = centre_distance % HELICAL_PERIOD
    return min(r, HELICAL_PERIOD - r) <= PHASE_TOLERANCE


def pair_index(units: list[str]) -> dict[tuple[str, str], int]:
    """Column block per unordered pair of pair-units, a family with itself included."""
    out: dict[tuple[str, str], int] = {}
    k = 0
    for i, a in enumerate(units):
        for b in units[i:]:
            out[(a, b)] = k
            k += 1
    return out


PAIR_FEATURES = len(GAP_BINS) + 2  # gap bins, same strand, in helical phase


def grammar_features(sites: list[tuple], pairs: dict[tuple[str, str], int]) -> dict[int, float]:
    """Per pair of families: non-overlapping site pairs per gap bin, on the same strand, in phase.

    Overlapping sites are skipped: two factors cannot occupy the same bases, and an overlap between
    two families is usually one stretch of sequence matching two related matrices.
    """
    raw: dict[int, int] = defaultdict(int)
    wanted = {u for pair in pairs for u in pair}
    rel = [s for s in sites if s[0] in wanted]
    for i, (ua, sa, ea, da) in enumerate(rel):
        for ub, sb, eb, db in rel[i + 1 :]:
            gap = sb - ea
            if gap < 0:
                continue
            block = pairs.get((ua, ub)) or pairs.get((ub, ua))
            if block is None:
                continue
            base = block * PAIR_FEATURES
            b = gap_bin(gap)
            if b is not None:
                raw[base + b] += 1
            if da == db:
                raw[base + len(GAP_BINS)] += 1
            if in_phase(((sb + eb) - (sa + ea)) / 2):
                raw[base + len(GAP_BINS) + 1] += 1
    return {k: math.log1p(v) for k, v in raw.items()}


def assemble(
    comp: list[float], sites: list[tuple], units: list[str], pairs: dict[tuple[str, str], int]
) -> dict[int, float]:
    """One element's sparse row: [1, composition (3), counts, grammar], the nested sets in order."""
    row: dict[int, float] = {0: 1.0}
    for i, v in enumerate(comp):
        row[1 + i] = v
    off = 4
    for i, v in count_features(sites, units).items():
        row[off + i] = v
    off += len(units) + 2
    for i, v in grammar_features(sites, pairs).items():
        row[off + i] = v
    return row


def row_features(
    seq: str, sites: list[tuple], units: list[str], pairs: dict[tuple[str, str], int]
) -> dict[int, float]:
    return assemble(composition(seq), sites, units, pairs)


def block_widths(units: list[str], pairs: dict[tuple[str, str], int]) -> dict[str, int]:
    """Number of leading columns of each nested model: (a), (b), (c)."""
    a = 4
    b = a + len(units) + 2
    return {"a": a, "b": b, "c": b + len(pairs) * PAIR_FEATURES}


# ---------------------------------------------------------------- ridge


def accumulate(rows: list[dict[int, float]], ys: dict[str, list[float]], k: int):
    """Normal equations of one chromosome: X'X (flat, upper filled) and X'y per cell line."""
    xtx = array("d", bytes(8 * k * k))
    xty = {c: array("d", bytes(8 * k)) for c in ys}
    for n, row in enumerate(rows):
        items = sorted(row.items())
        for c, y in ys.items():
            yv = y[n]
            t = xty[c]
            for i, v in items:
                t[i] += v * yv
        for p, (i, vi) in enumerate(items):
            base = i * k
            for j, vj in items[p:]:
                xtx[base + j] += vi * vj
    return xtx, xty


def add(into: array, other: array, sign: float = 1.0) -> None:
    for i, v in enumerate(other):
        if v:
            into[i] += sign * v


def cholesky(a: list[list[float]]) -> list[list[float]] | None:
    """Lower factor of a symmetric positive-definite matrix; None if it is not positive definite."""
    n = len(a)
    low: list[list[float]] = [[0.0] * n for _ in range(n)]
    for j in range(n):
        s = a[j][j] - sum(v * v for v in low[j][:j])
        if s <= 1e-12:
            return None
        d = math.sqrt(s)
        low[j][j] = d
        lj = low[j][:j]
        for i in range(j + 1, n):
            low[i][j] = (a[i][j] - sum(x * y for x, y in zip(low[i][:j], lj, strict=True))) / d
    return low


def cholesky_solve(low: list[list[float]], b: list[float]) -> list[float]:
    n = len(b)
    y = [0.0] * n
    for i in range(n):
        y[i] = (b[i] - sum(low[i][j] * y[j] for j in range(i))) / low[i][i]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(low[j][i] * x[j] for j in range(i + 1, n))) / low[i][i]
    return x


def penalised(xtx: array, k: int, cols: int, lam: float) -> list[list[float]]:
    """The leading `cols` block of X'X as a dense matrix, with lam on the diagonal in standardised units.

    Penalising column j by lam * var(j) is exactly ridge on standardised features with an unpenalised
    intercept, and needs no second pass over the data: the variances are inside X'X already.
    """
    a = [[0.0] * cols for _ in range(cols)]
    for i in range(cols):
        base = i * k
        for j in range(i, cols):
            v = xtx[base + j]
            a[i][j] = v
            a[j][i] = v
    n = a[0][0] or 1.0
    for j in range(1, cols):
        var = max(a[j][j] / n - (a[0][j] / n) ** 2, 0.0)
        a[j][j] += lam * var + 1e-9
    return a


def fit(xtx: array, xty: array, k: int, cols: int, lam: float) -> list[float] | None:
    low = cholesky(penalised(xtx, k, cols, lam))
    return cholesky_solve(low, list(xty[:cols])) if low else None


def predict(w: list[float], rows: list[dict[int, float]]) -> list[float]:
    cols = len(w)
    return [sum(v * w[i] for i, v in row.items() if i < cols) for row in rows]


# ---------------------------------------------------------------- statistics


def ranks(v: list[float]) -> list[float]:
    n = len(v)
    order = sorted(range(n), key=lambda i: v[i])
    out = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and v[order[j + 1]] == v[order[i]]:
            j += 1
        r = (i + j) / 2 + 1
        for t in range(i, j + 1):
            out[order[t]] = r
        i = j + 1
    return out


def spearman(xs: list[float], ys: list[float]) -> float | None:
    """Rank correlation, unrounded (the bootstrap works on differences of these)."""
    n = len(xs)
    if n < 3:
        return None
    rx, ry = ranks(xs), ranks(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(rx, ry, strict=True))
    vx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    vy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return cov / (vx * vy) if vx and vy else None


def interval(values: list[float], lo: float = 2.5, hi: float = 97.5) -> tuple[float, float]:
    s = sorted(values)
    n = len(s)

    def at(p: float) -> float:
        k = min(n - 1, max(0, int(round(p / 100 * (n - 1)))))
        return s[k]

    return at(lo), at(hi)


def bootstrap_models(
    truth: list[float], preds: dict[str, list[float]], draws: int = BOOTSTRAP, seed: int = SEED
) -> dict[str, Any]:
    """Spearman per model and the differences between them, resampling elements with replacement."""
    rng = random.Random(seed)
    n = len(truth)
    names = list(preds)
    obs = {m: spearman(preds[m], truth) for m in names}
    got: dict[str, list[float]] = defaultdict(list)
    for _ in range(draws):
        idx = [rng.randrange(n) for _ in range(n)]
        y = [truth[i] for i in idx]
        rho = {}
        for m in names:
            p = preds[m]
            rho[m] = spearman([p[i] for i in idx], y) or 0.0
        for m in names:
            got[m].append(rho[m])
        for a, b in (("b", "a"), ("c", "b"), ("c", "shuffled"), ("shuffled", "b")):
            if a in rho and b in rho:
                got[f"{a}-{b}"].append(rho[a] - rho[b])
    out: dict[str, Any] = {}
    for key, vals in got.items():
        lo, hi = interval(vals)
        point = obs[key] if key in obs else obs[key.split("-")[0]] - obs[key.split("-")[1]]
        out[key] = {
            "point": round(point, 4),
            "ci95": [round(lo, 4), round(hi, 4)],
            "above_zero": round(sum(1 for v in vals if v > 0) / len(vals), 3),
        }
    return out


# ---------------------------------------------------------------- the run


def sites_path(chrom: str, knowledge: Path = KNOWLEDGE) -> Path:
    return knowledge / f"grammar_sites_{chrom}.json"


def scan_chromosome(chrom: str, knowledge: Path = KNOWLEDGE) -> int:
    """Strict sites of every lentiMPRA element of one chromosome, cached locally (they are large)."""
    from genomeos.attribution.mpra import load_rows
    from genomeos.coords import Locus
    from genomeos.genome.fetch import REFERENCE
    from genomeos.genome.genome import Genome
    from genomeos.genome.motifs import load_families, load_motifs

    rows = load_rows(chrom, knowledge)
    fa = REFERENCE / f"{chrom}.fa.gz"
    if not rows or not fa.exists():
        return 0
    genome = Genome.from_fasta(fa)
    seqs = [str(genome.fetch(Locus(chrom, r["start"], r["end"]))) for r in rows]
    sites = scan_sites(load_motifs(relative=SITE_THRESHOLD), seqs, load_families())
    payload = {
        "chrom": chrom,
        "threshold": SITE_THRESHOLD,
        "elements": [
            {"key": r["key"], "gc": composition(s), "sites": [list(x) for x in st]}
            for r, s, st in zip(rows, seqs, sites, strict=True)
        ],
    }
    p = sites_path(chrom, knowledge)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(payload))
    return len(rows)


def load_sites(chrom: str, knowledge: Path = KNOWLEDGE) -> list[dict[str, Any]]:
    p = sites_path(chrom, knowledge)
    return json.loads(p.read_text())["elements"] if p.exists() else []


def training_chromosomes(chroms=CHROMS) -> list[str]:
    return [c for c in chroms if c not in HELD_OUT]


def vocabulary(
    per_chrom: dict[str, list[dict[str, Any]]], chroms: list[str]
) -> tuple[list[str], dict[tuple[str, str], int]]:
    """Family units ranked by the share of *training* elements carrying a site, and the pair blocks."""
    n = 0
    seen: dict[str, int] = defaultdict(int)
    for c in chroms:
        for e in per_chrom[c]:
            n += 1
            for u in {s[0] for s in e["sites"]}:
                seen[u] += 1
    order = sorted(seen, key=lambda u: (-seen[u], u))
    units = order[:COUNT_UNITS]
    return units, pair_index(units[:PAIR_UNITS])


def totals(blocks: dict[str, tuple[array, dict[str, array]]], chroms: list[str], k: int, cells=CELLS):
    """Normal equations summed over a set of chromosomes."""
    xtx = array("d", bytes(8 * k * k))
    xty = {cell: array("d", bytes(8 * k)) for cell in cells}
    for c in chroms:
        add(xtx, blocks[c][0])
        for cell in cells:
            add(xty[cell], blocks[c][1][cell])
    return xtx, xty


def element_rows(
    elements: list[dict[str, Any]], units: list[str], pairs: dict[tuple[str, str], int], rng: random.Random
) -> tuple[list[dict[int, float]], list[dict[int, float]]]:
    """Real and shuffled-grammar rows of one chromosome's elements (the shuffle keeps every count)."""
    real, shuf = [], []
    for e in elements:
        sites = [tuple(s) for s in e["sites"]]
        real.append(assemble(e["gc"], sites, units, pairs))
        shuf.append(assemble(e["gc"], shuffle_sites(sites, rng), units, pairs))
    return real, shuf


def fold_assignment(sizes: dict[str, int], folds: int = INNER_FOLDS) -> dict[str, int]:
    """Training chromosomes split into balanced folds, largest first: the penalty is chosen inside them."""
    load = [0] * folds
    out: dict[str, int] = {}
    for c in sorted(sizes, key=lambda c: (-sizes[c], c)):
        f = min(range(folds), key=lambda i: load[i])
        out[c] = f
        load[f] += sizes[c]
    return out


def inner_folds(
    blocks: dict[str, tuple[array, dict[str, array]]],
    rows: dict[str, list[dict[int, float]]],
    ys: dict[str, dict[str, list[float]]],
    k: int,
    folds: dict[str, int],
    cells=CELLS,
) -> list[tuple]:
    """Per fold: the other folds' normal equations, and the fold's own rows and activity."""
    train = list(folds)
    total_x, total_y = totals(blocks, train, k, cells)
    out = []
    for f in range(max(folds.values()) + 1):
        te = [c for c in train if folds[c] == f]
        xtx = array("d", total_x)
        xty = {cell: array("d", total_y[cell]) for cell in cells}
        for c in te:
            add(xtx, blocks[c][0], -1.0)
            for cell in cells:
                add(xty[cell], blocks[c][1][cell], -1.0)
        out.append(
            (
                xtx,
                xty,
                [r for c in te for r in rows[c]],
                {cell: [v for c in te for v in ys[c][cell]] for cell in cells},
            )
        )
    return out


def choose_lambda(prepared: list[tuple], k: int, cols: int, cells=CELLS) -> dict[str, dict[str, Any]]:
    """Per cell line: the penalty with the best mean fold-by-chromosome Spearman inside the training set."""
    scores: dict[str, dict[str, list[float]]] = {cell: defaultdict(list) for cell in cells}
    for lam in LAMBDAS:
        for xtx, xty, xs, truth in prepared:
            low = cholesky(penalised(xtx, k, cols, lam))
            if low is None:
                continue
            for cell in cells:
                w = cholesky_solve(low, list(xty[cell][:cols]))
                rho = spearman(predict(w, xs), truth[cell])
                if rho is not None:
                    scores[cell][str(lam)].append(rho)
    best: dict[str, dict[str, Any]] = {}
    for cell, rec in scores.items():
        means = {lam: sum(v) / len(v) for lam, v in rec.items()}
        if not means:
            continue
        lam = max(means, key=lambda x: means[x])
        best[cell] = {
            "lambda": float(lam),
            "inner_spearman": round(means[lam], 4),
            "by_lambda": {x: round(v, 4) for x, v in means.items()},
        }
    return best


def run(
    knowledge: Path = KNOWLEDGE,
    results_dir: Path = RESULTS_DIR,
    score_held_out: bool = True,
    progress=None,
) -> dict[str, Any]:
    """Build the three feature sets, choose the penalty inside the training set, score the held-out folds."""
    say = progress or (lambda _m: None)
    chroms = [c for c in CHROMS if sites_path(c, knowledge).exists()]
    per_chrom = {c: load_sites(c, knowledge) for c in chroms}
    train = [c for c in training_chromosomes(chroms) if per_chrom[c]]
    test = [c for c in chroms if c in HELD_OUT and per_chrom[c]]
    units, pairs = vocabulary(per_chrom, train)
    widths = block_widths(units, pairs)
    k = widths["c"]
    say(f"{sum(len(per_chrom[c]) for c in chroms):,} elements, {len(units)} family units, {k} columns")
    rng = random.Random(SEED)
    activity = {c: {cell: [] for cell in CELLS} for c in chroms}
    rows: dict[str, list[dict[int, float]]] = {}
    srows: dict[str, list[dict[int, float]]] = {}
    keep: dict[str, list[int]] = {}
    from genomeos.attribution.mpra import load_rows

    for c in chroms:
        meas = {r["key"]: r["activity"] for r in load_rows(c, knowledge)}
        idx = [i for i, e in enumerate(per_chrom[c]) if all(cell in meas.get(e["key"], {}) for cell in CELLS)]
        keep[c] = idx
        els = [per_chrom[c][i] for i in idx]
        for cell in CELLS:
            activity[c][cell] = [meas[e["key"]][cell] for e in els]
        rows[c], srows[c] = element_rows(els, units, pairs, rng)
        say(f"  {c}: {len(els):,} elements with activity in all three cell lines")
    blocks = {c: accumulate(rows[c], activity[c], k) for c in chroms}
    sblocks = {c: accumulate(srows[c], activity[c], k) for c in chroms}
    folds = fold_assignment({c: len(rows[c]) for c in train})
    say("choosing the penalty inside the training chromosomes")
    prepared = inner_folds(blocks, rows, activity, k, folds)
    sprepared = inner_folds(sblocks, srows, activity, k, folds)
    chosen = {}
    for model in ("a", "b", "c"):
        chosen[model] = choose_lambda(prepared, k, widths[model])
        say(
            f"  ({model}) inner: " + ", ".join(f"{c} {v['inner_spearman']}" for c, v in chosen[model].items())
        )
    chosen["shuffled"] = choose_lambda(sprepared, k, widths["c"])
    say(
        "  (shuffled) inner: "
        + ", ".join(f"{c} {v['inner_spearman']}" for c, v in chosen["shuffled"].items())
    )
    out: dict[str, Any] = {
        "preregistered": PREREGISTERED,
        "site_threshold": SITE_THRESHOLD,
        "training_chromosomes": train,
        "held_out_chromosomes": test,
        "training_elements": sum(len(rows[c]) for c in train),
        "columns": widths,
        "family_units": units,
        "pair_units": units[:PAIR_UNITS],
        "pair_features_per_pair": PAIR_FEATURES,
        "gap_bins": [list(b) for b in GAP_BINS],
        "sites_per_element_median": _median(
            [len(e["sites"]) for c in chroms for e in per_chrom[c]],
        ),
        "inner_cross_validation": chosen,
        "folds": folds,
        "evidence": EVIDENCE,
    }
    if not score_held_out or not test:
        out["held_out"] = None
        return out
    say(f"scoring the held-out chromosomes: {', '.join(test)}")
    train_x, train_y = totals(blocks, train, k)
    train_sx, train_sy = totals(sblocks, train, k)
    test_rows = [r for c in test for r in rows[c]]
    test_srows = [r for c in test for r in srows[c]]
    cells_out: dict[str, Any] = {}
    for cell in CELLS:
        truth = [v for c in test for v in activity[c][cell]]
        preds = {}
        for model in ("a", "b", "c"):
            lam = chosen[model][cell]["lambda"]
            w = fit(train_x, train_y[cell], k, widths[model], lam)
            preds[model] = predict(w, test_rows) if w else [0.0] * len(truth)
        lam = chosen["shuffled"][cell]["lambda"]
        w = fit(train_sx, train_sy[cell], k, widths["c"], lam)
        preds["shuffled"] = predict(w, test_srows) if w else [0.0] * len(truth)
        cells_out[cell] = {
            "held_out_elements": len(truth),
            "lambda": {m: chosen[m][cell]["lambda"] for m in ("a", "b", "c", "shuffled")},
            **bootstrap_models(truth, preds),
        }
        say(
            f"  {cell}: a {cells_out[cell]['a']['point']:.3f} b {cells_out[cell]['b']['point']:.3f} "
            f"c {cells_out[cell]['c']['point']:.3f} shuffled {cells_out[cell]['shuffled']['point']:.3f}"
        )
    out["held_out"] = cells_out
    out["verdict"] = verdict(cells_out)
    return out


def verdict(cells: dict[str, Any]) -> dict[str, Any]:
    """The pre-registered claim, decided: every cell line's interval above 0, counts and arrangement."""
    beats_counts = {c: cells[c]["c-b"]["ci95"][0] > 0 for c in cells}
    beats_shuffle = {c: cells[c]["c-shuffled"]["ci95"][0] > 0 for c in cells}
    passed = all(beats_counts.values())
    return {
        "grammar_beats_counts": beats_counts,
        "grammar_beats_shuffled_grammar": beats_shuffle,
        "preregistered_claim_holds": passed,
        "arrangement_reading_holds": passed and all(beats_shuffle.values()),
        "statement": (
            "arrangement (spacing, orientation, helical phase) adds measurable signal over counts"
            if passed and all(beats_shuffle.values())
            else "pair features add signal over counts, but their shuffled control adds as much: the gain "
            "is co-occurrence, not arrangement"
            if passed
            else "grammar features add nothing detectable over strict motif counts"
        ),
    }


def _median(values: list[float]) -> float | None:
    s = sorted(values)
    return s[len(s) // 2] if s else None


def run_and_save(
    knowledge: Path = KNOWLEDGE,
    results_dir: Path = RESULTS_DIR,
    score_held_out: bool = True,
    progress=None,
) -> dict[str, Any]:
    out = run(knowledge, results_dir, score_held_out, progress)
    name = "motif_grammar" if score_held_out else "motif_grammar_preregistration"
    save_result(name, out, results_dir)
    return out
