# SPDX-License-Identifier: AGPL-3.0-or-later
"""Element types discovered from deletion behaviour, held against the registry's classes.

The registry names an element enhancer-like or promoter-like from its biochemistry. The project's
deletion scoring says the name is a weak guide to what the element does: the class predicts
whether a deletion lowers or raises its gene at AUC 0.532 (docs/NODES-READER-WRITER.md). This
module asks whether categories read off the behaviour itself do better.

1. A behaviour vector per element, from the AlphaGenome deletion archive only: the signed effect on
   every gene in the 1 Mb window on each of four cell lines' own RNA-seq tracks (K562, HepG2,
   GM12878, IMR-90), plus where the most-moved gene's start sits. Nothing from the registry class,
   the marks, DNase, sequence or constraint enters it.
2. k-means on the standardised vector. k is chosen by a criterion fixed before any result was read
   (`choose_k`): the highest mean silhouette on a fixed subsample over k = 2..10, among the k whose
   smallest cluster holds at least MIN_CLUSTER_SHARE of the elements.
3. What the clusters are: their mix of registry classes, measured chromatin states, budget tiers
   and constraint, beside their behaviour centroids.
4. Controls: the vector with each feature permuted on its own (does the joint structure exceed the
   marginals?); whole vectors shuffled among elements of the same registry class and DNase openness
   (do the clusters' associations and predictions survive losing the element?); centroids fitted
   on two chromosomes and applied unchanged to others; and how much length, GC and distance to the
   nearest gene explain alone.
5. The held-out task of the epigenome layer: per element and cell line, whether the deletion acts,
   whether it raises the gene, and by how much. The type used for a cell line is assigned from the
   other three lines only, so the outcome never enters the vector that names the type.

Everything the vector holds is `predicted` evidence from one model. AlphaGenome was trained on ENCODE
tracks of these same cell lines, so agreement between these types and ENCODE-derived marks is partly
the model reading back its inputs; and a type that predicts the model's output in one line from its
output in three others measures the model's consistency across lines, not an element's biology.

Clustering and silhouettes use numpy, imported inside the functions that need it (as the lexicon
does). The features, ridge models and scores are standard library and reuse the epigenome layer's
statistics, so the numbers sit beside its table.
"""

from __future__ import annotations

import bisect
import gzip
import json
import math
import random
from collections import Counter
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

ARCHIVE = Path("data/knowledge/alphagenome/elements")  # every gene in the window, per cell line
CHAIN = Path("data/knowledge/alphagenome/all_elements")  # the chain's target, nearest gene, outcome
REFERENCE = Path("data/reference")
CONSTRAINT = Path("data/knowledge/constraint")
RESULTS = Path("data/results")
CELLS = ("K562", "HepG2", "GM12878", "IMR-90")
ACTS = 0.1  # |log2 fold change| the project calls an effect
FLOOR = 1e-3  # added to a magnitude before its log
K_RANGE = tuple(range(2, 11))
MIN_CLUSTER_SHARE = 0.02
SILHOUETTE_SAMPLE = 5000
SEED = 21
FLANK = 500  # marks are read over the element midpoint +/- FLANK, as in the epigenome layer's task

FEATURES = (
    "log_strongest",  # log10 of the largest |effect| on any gene in any line
    "sign_strongest",  # +1 when that effect is a rise on deletion, -1 when a drop
    "net_direction",  # sum of signed effects over sum of |effects|, every gene and line
    "log_genes_moving",  # log1p of the genes moved by ACTS or more in any line
    "share_window_moving",  # those genes over the genes in the window
    "lines_acting",  # share of lines whose strongest effect reaches ACTS
    "direction_differs",  # 1 when acting lines disagree in sign
    "log_line_spread",  # log10 of (largest - smallest) of the lines' strongest signed effects
    "acts_on_nearest",  # 1 when the most-moved gene is the nearest TSS in the domain
    "log_target_distance",  # log10 of 1 + bases from the element midpoint to the most-moved gene's start
    "log_mean_effect",  # log10 of the mean |effect| over genes and lines: how diffuse the deletion is
)


class Cost:
    """Bytes read from disk, counted by file size as each file is opened."""

    def __init__(self) -> None:
        self.bytes_read = 0
        self.files = 0

    def add(self, path: Path) -> Path:
        if path.exists():
            self.bytes_read += path.stat().st_size
            self.files += 1
        return path


COST = Cost()


def _numpy():
    import numpy as np

    return np


# ------------------------------------------------------------------------------------------
# Loading: the archive, the chain's table, the registry, GENCODE starts, sequence, constraint
# ------------------------------------------------------------------------------------------


def gene_starts(chrom: str, reference: Path = REFERENCE) -> dict[str, int]:
    """0-based TSS per gene symbol (the archive names unnamed genes by their Ensembl id, as GENCODE does)."""
    from genomeos.genome.annotation import iter_gff3

    out: dict[str, int] = {}
    p = COST.add(reference / f"gencode_v50_{chrom}.gff3.gz")
    if not p.exists():
        return out
    for row in iter_gff3(p, types={"gene"}):
        name = row.attrs.get("gene_name") or row.attrs.get("gene_id", "").split(".")[0]
        out.setdefault(name, row.end - 1 if row.strand == "-" else row.start - 1)
    return out


def registry_classes(chrom: str, results: Path = RESULTS) -> dict[str, tuple[str, bool]]:
    out = {}
    with gzip.open(COST.add(results / f"ccres_{chrom}.bed.gz"), "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            out[f[3]] = (f[4], f[5] == "1")
    return out


def summarise_cells(genes: Sequence[dict], cells: Sequence[str] = CELLS) -> dict[str, dict]:
    """Per line: the strongest signed effect and its gene, the genes moving, and the effect sums."""
    per: dict[str, dict] = {}
    for cell in cells:
        best, gene, moving, total, total_abs, n = 0.0, None, [], 0.0, 0.0, 0
        for g in genes:
            v = (g.get("by_cell") or {}).get(cell)
            if v is None:
                continue
            n += 1
            total += v
            total_abs += abs(v)
            if abs(v) >= ACTS:
                moving.append(g["gene"])
            if abs(v) > abs(best):
                best, gene = v, g["gene"]
        if n:
            per[cell] = {
                "top": best,
                "gene": gene,
                "moving": tuple(moving),
                "sum": total,
                "abs": total_abs,
                "n": n,
            }
    return per


def behaviour(el: dict, cells: Sequence[str] = CELLS) -> list[float] | None:
    """The behaviour vector (FEATURES) of one element over `cells`, or None when no line was scored."""
    per = [(c, el["cells"][c]) for c in cells if c in el["cells"]]
    if not per:
        return None
    tops = [p["top"] for _, p in per]
    strongest_cell, strongest = max(per, key=lambda cp: abs(cp[1]["top"]))
    top = strongest["top"]
    moving = set().union(*(p["moving"] for _, p in per))
    acting = [t for t in tops if abs(t) >= ACTS]
    total_abs = sum(p["abs"] for _, p in per)
    total = sum(p["sum"] for _, p in per)
    n_genes = max(p["n"] for _, p in per)
    distance = el["target_distance"].get(strongest_cell)
    if distance is None:
        distance = el.get("nearest_distance")
    return [
        math.log10(abs(top) + FLOOR),
        1.0 if top > 0 else -1.0,
        total / total_abs if total_abs else 0.0,
        math.log1p(len(moving)),
        len(moving) / n_genes,
        len(acting) / len(per),
        1.0 if any(t > 0 for t in acting) and any(t < 0 for t in acting) else 0.0,
        math.log10(max(tops) - min(tops) + FLOOR),
        1.0 if strongest["gene"] is not None and strongest["gene"] == el.get("nearest_gene") else 0.0,
        math.log10(1 + distance) if distance is not None else 5.0,
        math.log10(total_abs / (len(per) * n_genes) + FLOOR),
    ]


def _fai(path: Path) -> dict[str, tuple[int, int, int, int]]:
    out = {}
    for line in COST.add(path).read_text().splitlines():
        f = line.split("\t")
        out[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
    return out


def gc_of(
    chrom: str, intervals: Sequence[tuple[int, int]], reference: Path = REFERENCE
) -> list[float | None]:
    """GC fraction over called bases of each 0-based half-open interval, read by seek through the .fai.

    The bytes come from whichever form of the chromosome is cached: the flat
    `.fa` seeked directly, or the blocked `.fa.gz` seeked through its `.gzi`.
    """
    from genomeos.genome.index import resolve_fasta

    fa = resolve_fasta(reference / f"{chrom}.fa")
    fai = fa.with_name(fa.name + ".fai")
    if not fa.exists() or not fai.exists():
        return [None] * len(intervals)
    _length, offset, per_line, line_bytes = _fai(fai)[chrom]
    out: list[float | None] = []
    fh = _seekable(fa)
    try:
        for s, e in intervals:
            a = offset + (s // per_line) * line_bytes + s % per_line
            b = offset + (e // per_line) * line_bytes + e % per_line
            fh.seek(a)
            chunk = fh.read(b - a)
            COST.bytes_read += len(chunk)
            seq = chunk.replace(b"\n", b"").upper()
            acgt = sum(seq.count(x) for x in (b"A", b"C", b"G", b"T"))
            out.append((seq.count(b"G") + seq.count(b"C")) / acgt if acgt else None)
    finally:
        fh.close()
    return out


def _seekable(fa: Path):
    """A binary handle over the FASTA bytes that supports seek, compressed or not."""
    from genomeos.genome import bgzf

    if fa.suffix == ".gz" and bgzf.is_bgzf(fa):
        return bgzf.BgzfReader(fa)
    return open(fa, "rb")  # noqa: SIM115  (closed by the caller)


def load_chromosome(chrom: str) -> list[dict]:
    """The elements of one chromosome the chain scored per cell line, with behaviour and context."""
    from genomeos.attribution.constraint import elements_over_blocks, load_elements

    chain = {e["id"]: e for e in json.loads(COST.add(CHAIN / f"{chrom}.json").read_text())}
    classes = registry_classes(chrom)
    starts = gene_starts(chrom)
    tss_sorted = sorted(starts.values())
    with gzip.open(COST.add(ARCHIVE / f"{chrom}.json.gz"), "rt") as fh:
        archive = json.load(fh)
    rows: list[dict] = []
    for eid, hit in archive.items():
        c = chain.get(eid)
        if not c or not c.get("predicted_by_cell") or eid not in classes:
            continue
        mid = (c["start"] + c["end"]) // 2
        cells = summarise_cells(hit.get("genes") or [])
        if not cells:
            continue
        inferred = c.get("inferred") or {}
        j = bisect.bisect_left(tss_sorted, mid)
        near = [abs(tss_sorted[i] - mid) for i in (j - 1, j) if 0 <= i < len(tss_sorted)]
        rows.append(
            {
                "id": eid,
                "chrom": chrom,
                "start": c["start"],
                "end": c["end"],
                "length": c["end"] - c["start"],
                "cls": classes[eid][0],
                "ctcf": classes[eid][1],
                "outcome": c["predicted_by_cell"],
                "model_seconds": hit.get("seconds"),
                "cells": cells,
                "nearest_gene": inferred.get("gene"),
                "nearest_distance": inferred.get("distance"),
                "nearest_any_gene_distance": min(near) if near else None,
                "target_distance": {
                    cell: abs(starts[p["gene"]] - mid)
                    for cell, p in cells.items()
                    if p["gene"] is not None and p["gene"] in starts
                },
            }
        )
    del archive
    rows.sort(key=lambda r: r["start"])
    for r, g in zip(rows, gc_of(chrom, [(r["start"], r["end"]) for r in rows]), strict=True):
        r["gc"] = g
    cons = CONSTRAINT / f"phastConsElements100way_{chrom}.bed.gz"
    if cons.exists() and rows:
        conserved = load_elements(COST.add(cons))
        for r, o in zip(
            rows, elements_over_blocks(conserved, [(r["start"], r["end"]) for r in rows]), strict=True
        ):
            r["conserved_fraction"] = o["fraction"]
    tiers = _block_tiers(chrom)
    for r in rows:
        r["tier"] = _at(tiers, (r["start"] + r["end"]) // 2)
    return rows


def _block_tiers(chrom: str) -> list[tuple[int, int, str]]:
    p = COST.add(RESULTS / f"budget_{chrom}.json")
    if not p.exists():
        return []
    blocks = json.loads(p.read_text())["blocks"]
    return sorted((b["start"], b["end"], (b.get("guess") or {}).get("tier") or "unknown") for b in blocks)


def _at(blocks: list[tuple[int, int, str]], pos: int) -> str | None:
    i = bisect.bisect_right(blocks, (pos, math.inf, "")) - 1
    return blocks[i][2] if 0 <= i < len(blocks) and blocks[i][0] <= pos < blocks[i][1] else None


# ------------------------------------------------------------------------------------------
# Clustering: scaling, k-means, silhouette, the k criterion, agreement between labelings
# ------------------------------------------------------------------------------------------


def scaler(x: Sequence[Sequence[float]]) -> tuple[list[float], list[float]]:
    """Column means and standard deviations (1 where a column is constant)."""
    n, k = len(x), len(x[0])
    means = [sum(r[j] for r in x) / n for j in range(k)]
    sds = [math.sqrt(sum((r[j] - means[j]) ** 2 for r in x) / n) or 1.0 for j in range(k)]
    return means, sds


def scale(x: Sequence[Sequence[float]], means: Sequence[float], sds: Sequence[float]) -> list[list[float]]:
    return [[(v - m) / s for v, m, s in zip(r, means, sds, strict=True)] for r in x]


def _distances(np, a, c):
    return np.maximum((a * a).sum(1)[:, None] - 2 * a @ c.T + (c * c).sum(1)[None, :], 0.0)


def kmeans(
    x: Sequence[Sequence[float]], k: int, seed: int = SEED, restarts: int = 4, iters: int = 100
) -> dict[str, Any]:
    """Lloyd's k-means with k-means++ seeding; the best of `restarts` by within-cluster sum of squares."""
    np = _numpy()
    a = np.asarray(x, dtype=float)
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(restarts):
        centres = [a[rng.integers(len(a))]]
        d2 = ((a - centres[0]) ** 2).sum(1)
        for _ in range(k - 1):
            total = d2.sum()
            i = rng.choice(len(a), p=d2 / total) if total > 0 else rng.integers(len(a))
            centres.append(a[i])
            d2 = np.minimum(d2, ((a - a[i]) ** 2).sum(1))
        c = np.array(centres)
        labels = None
        for _ in range(iters):
            new = _distances(np, a, c).argmin(1)
            if labels is not None and (new == labels).all():
                break
            labels = new
            for j in range(k):
                members = a[labels == j]
                if len(members):
                    c[j] = members.mean(0)
        inertia = float(_distances(np, a, c)[np.arange(len(a)), labels].sum())
        if best is None or inertia < best["inertia"]:
            best = {"centroids": c.tolist(), "labels": labels.tolist(), "inertia": inertia}
    return best


def assign(x: Sequence[Sequence[float]], centroids: Sequence[Sequence[float]]) -> list[int]:
    np = _numpy()
    return _distances(np, np.asarray(x, dtype=float), np.asarray(centroids, dtype=float)).argmin(1).tolist()


def silhouette(
    x: Sequence[Sequence[float]], labels: Sequence[int], sample: int = SILHOUETTE_SAMPLE, seed: int = SEED
) -> float | None:
    """Mean silhouette over a fixed random subsample, distances within the subsample (Euclidean)."""
    np = _numpy()
    a = np.asarray(x, dtype=float)
    lab = np.asarray(labels)
    idx = np.random.default_rng(seed).permutation(len(a))[: min(sample, len(a))]
    a, lab = a[idx], lab[idx]
    ks = np.unique(lab)
    if len(ks) < 2:
        return None
    sums = np.zeros((len(a), len(ks)))
    for s in range(0, len(a), 1000):
        d = np.sqrt(_distances(np, a[s : s + 1000], a))
        for j, kk in enumerate(ks):
            sums[s : s + 1000, j] = d[:, lab == kk].sum(1)
    counts = np.array([(lab == kk).sum() for kk in ks], dtype=float)
    own = np.searchsorted(ks, lab)
    rows = np.arange(len(a))
    own_n = counts[own] - 1
    intra = np.where(own_n > 0, sums[rows, own] / np.maximum(own_n, 1), 0.0)
    other = sums / counts[None, :]
    other[rows, own] = np.inf
    inter = other.min(1)
    s = np.where(own_n > 0, (inter - intra) / np.maximum(np.maximum(inter, intra), 1e-12), 0.0)
    return float(s.mean())


def choose_k(x: Sequence[Sequence[float]], ks: Iterable[int] = K_RANGE, seed: int = SEED) -> dict[str, Any]:
    """The criterion fixed in advance: highest subsample silhouette among the k whose smallest cluster
    holds at least MIN_CLUSTER_SHARE of the elements."""
    curve = {}
    fits = {}
    for k in ks:
        fit = kmeans(x, k, seed)
        sizes = Counter(fit["labels"])
        smallest = min(sizes.get(j, 0) for j in range(k)) / len(x)
        curve[k] = {
            "silhouette": round(silhouette(x, fit["labels"], seed=seed) or 0.0, 4),
            "smallest_share": round(smallest, 4),
            "inertia": round(fit["inertia"], 1),
        }
        fits[k] = fit
    eligible = [k for k, c in curve.items() if c["smallest_share"] >= MIN_CLUSTER_SHARE]
    k = max(eligible or list(curve), key=lambda kk: curve[kk]["silhouette"])
    return {"k": k, "curve": curve, "fit": fits[k]}


def adjusted_rand(a: Sequence[int], b: Sequence[int]) -> float:
    """Adjusted Rand index of two labelings of the same items."""

    def c2(n: float) -> float:
        return n * (n - 1) / 2

    n = len(a)
    pairs = Counter(zip(a, b, strict=True))
    sum_ij = sum(c2(v) for v in pairs.values())
    sum_a = sum(c2(v) for v in Counter(a).values())
    sum_b = sum(c2(v) for v in Counter(b).values())
    expected = sum_a * sum_b / c2(n) if n > 1 else 0.0
    top = (sum_a + sum_b) / 2
    return (sum_ij - expected) / (top - expected) if top != expected else 1.0


def matched_agreement(a: Sequence[int], b: Sequence[int]) -> float:
    """Share of items on which two labelings agree after matching each label of `b` greedily to the
    label of `a` it shares most items with (a readable companion to the Rand index)."""
    pairs = Counter(zip(a, b, strict=True))
    used_a: set[int] = set()
    used_b: set[int] = set()
    agree = 0
    for (la, lb), n in pairs.most_common():
        if la in used_a or lb in used_b:
            continue
        used_a.add(la)
        used_b.add(lb)
        agree += n
    return agree / len(a) if a else 0.0


def cramers_v(a: Sequence[Any], b: Sequence[Any]) -> float | None:
    """Association of two categorical labelings, 0 (none) to 1 (one determines the other)."""
    n = len(a)
    if not n:
        return None
    ca, cb = Counter(a), Counter(b)
    if len(ca) < 2 or len(cb) < 2:
        return None
    obs = Counter(zip(a, b, strict=True))
    chi2 = 0.0
    for x, nx in ca.items():
        for y, ny in cb.items():
            e = nx * ny / n
            chi2 += (obs.get((x, y), 0) - e) ** 2 / e
    return math.sqrt(chi2 / (n * (min(len(ca), len(cb)) - 1)))


# ------------------------------------------------------------------------------------------
# Controls
# ------------------------------------------------------------------------------------------


def column_shuffle(x: Sequence[Sequence[float]], seed: int = SEED) -> list[list[float]]:
    """Each feature permuted on its own: the marginals kept, every dependence between features lost."""
    rng = random.Random(seed)
    cols = [list(c) for c in zip(*x, strict=True)]
    for c in cols:
        rng.shuffle(c)
    return [list(r) for r in zip(*cols, strict=True)]


def stratified_donors(keys: Sequence[Any], seed: int = SEED) -> list[int]:
    """For each item, the index of another item of the same stratum whose data it takes (a permutation
    inside every stratum)."""
    rng = random.Random(seed)
    strata: dict[Any, list[int]] = {}
    for i, key in enumerate(keys):
        strata.setdefault(key, []).append(i)
    donor = list(range(len(keys)))
    for members in strata.values():
        shuffled = members[:]
        rng.shuffle(shuffled)
        for i, d in zip(members, shuffled, strict=True):
            donor[i] = d
    return donor


# ------------------------------------------------------------------------------------------
# The held-out task: chromatin features as the epigenome layer reads them, ridge models, scores
# ------------------------------------------------------------------------------------------


def chromatin(layer, cell: str, start: int, end: int) -> dict[str, Any]:
    """DNase over the element; peak calls, fold change and methylation over midpoint +/- FLANK; the state.

    The same reading as the epigenome layer's direction task. Fold change is None where no signal
    profile is cached for the chromosome."""
    from genomeos.genome import epigenome as ep

    mid = (start + end) // 2
    lo, hi = mid - FLANK, mid + FLANK
    dn = layer.dnase.get(cell)
    hits = dn.overlapping(start, end) if dn else []
    out: dict[str, Any] = {
        "dnase_peak": 1.0 if hits else 0.0,
        "dnase_signal": math.log1p(max((v for _, _, v in hits), default=0.0)),
    }
    for m in ep.MARKS:
        idx = layer.peaks.get((cell, m))
        out[f"{m}_peak"] = 1.0 if idx and idx.overlapping(lo, hi) else 0.0
        prof = layer.signal.get((cell, m))
        fc = ep.profile_mean(prof, lo, hi) if prof is not None else None
        out[f"{m}_fc"] = math.log1p(fc) if fc is not None else None
    cols = layer.meth.get(cell)
    out["methylation"] = ep.methylation_over(cols, lo, hi)["fraction"] if cols is not None else None
    out["state"] = ep.chromatin_state({m: bool(out[f"{m}_peak"]) for m in ep.MARKS})
    return out


def impute(rows: list[dict], keys: Iterable[str]) -> None:
    """Median for a missing value, with a missing flag column."""
    for k in keys:
        vals = sorted(r[k] for r in rows if r[k] is not None)
        med = vals[len(vals) // 2] if vals else 0.0
        for r in rows:
            r[k + "_missing"] = 1.0 if r[k] is None else 0.0
            if r[k] is None:
                r[k] = med


def one_hot(label: int, k: int) -> list[float]:
    """k - 1 indicator columns (the first type is the reference)."""
    return [1.0 if label == j else 0.0 for j in range(1, k)]


def fit_predict(
    train_x: Sequence[Sequence[float]],
    train_y: Sequence[float],
    test_x: Sequence[Sequence[float]],
    lam: float = 1.0,
) -> list[float]:
    """A ridge fit on standardised training columns, applied with the training scale to new rows."""
    from genomeos.genome import epigenome as ep

    means, sds = scaler(train_x)
    w = ep.ridge_fit([[1.0, *r] for r in scale(train_x, means, sds)], list(train_y), lam)
    return ep.predict(w, [[1.0, *r] for r in scale(test_x, means, sds)])


def score(task: str, pred: Sequence[float], y: Sequence[float]) -> float | None:
    from genomeos.genome import epigenome as ep

    s = ep.spearman(list(pred), list(y)) if task == "magnitude" else ep.auc(list(pred), [int(v) for v in y])
    return round(s, 4) if s is not None else None


def tasks_of(units: Sequence[dict]) -> dict[str, tuple[list[float], list[int]]]:
    """The three questions of the epigenome layer's table: acts, rise among acting, magnitude."""
    acting = [i for i, u in enumerate(units) if abs(u["lfc"]) >= ACTS]
    return {
        "acts": ([1.0 if abs(u["lfc"]) >= ACTS else 0.0 for u in units], list(range(len(units)))),
        "rise_among_acting": ([1.0 if units[i]["lfc"] > 0 else 0.0 for i in acting], acting),
        "magnitude": ([math.log10(abs(u["lfc"]) + FLOOR) for u in units], list(range(len(units)))),
    }


def share_table(labels: Sequence[int], values: Sequence[Any]) -> dict[str, dict[str, float]]:
    """Per cluster, the share of each value (None values counted as 'unmeasured')."""
    counts: dict[int, Counter] = {}
    for lab, v in zip(labels, values, strict=True):
        counts.setdefault(lab, Counter())["unmeasured" if v is None else str(v)] += 1
    out = {}
    for lab in sorted(counts):
        n = sum(counts[lab].values())
        out[str(lab)] = {k: round(c / n, 4) for k, c in counts[lab].most_common()}
    return out


__all__ = [
    "CELLS",
    "FEATURES",
    "adjusted_rand",
    "assign",
    "behaviour",
    "choose_k",
    "column_shuffle",
    "cramers_v",
    "kmeans",
    "load_chromosome",
    "matched_agreement",
    "silhouette",
    "stratified_donors",
    "summarise_cells",
]
