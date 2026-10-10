# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Simpler genomes, stage 2: the worm's constraint tiers, derived the way the human ones were.

THE HUMAN TIERS. A human base is "constrained" when its Zoonomia phyloP score is positive and its
q value is at most 0.05 (Christmas et al. 2023, Science 380:eabn3943, Methods: "We computed a
mammalian phyloP threshold by converting the p values corresponding to the phyloP scores into q
values using a FDR correction. We considered any column with a resulting q <= 0.05 to be
significantly evolutionarily constrained or accelerated, as determined by the sign of the score").
That procedure gives 2.27 on the 241-mammal track (`constraint.PHYLOP_THRESHOLD`). A block of the
regulatory class is "under constraint" when at least 5% of its measured bases pass
(`budget.CONSTRAINED_MIN`, `budget._rule`); otherwise selection is weak or not detected. The number
2.27 belongs to that track and does not transfer.

THE WORM TRACKS. UCSC's ce11 phyloP135way (134 genomes aligned to C. elegans: 111 nematodes, 22
flatworms, Ciona) and phyloP26way (26 nematodes). Neither the READMEs nor the track descriptions
give a significance cut-off; they say the scores are -log10 p-values under a null of neutral
evolution, computed per column by a likelihood-ratio test (--method LRT) in CONACC mode, with the
neutral model fitted by phyloFit on 4-fold degenerate sites (checked 2026-10-10). Neither track has
a publication of its own. So the cut is not borrowed and not chosen: it is re-derived here by
Zoonomia's own procedure applied to each track's own scores. p = 10^-|score| for every scored
base, Benjamini-Hochberg at q = 0.05 over all scored bases of both signs, constrained = positive
score at or above the resulting cut.

READ BY RANGE, NEVER DOWNLOADED. The score distribution comes from a uniform sample of the track's
data blocks: the blocks of the six nuclear chromosomes, in file order, are cut into consecutive
groups of SAMPLE_GROUP and SAMPLE_FRACTION of the groups are drawn with a fixed seed. Every scored
base sits in exactly one block, so every scored base has the same chance of being in the sample. The
cut's sampling interval is a bootstrap over 1-Mb windows of the sample. The elements are then read
block by block through `bigwig.BigWig`'s index, and nothing is written to disk but the counts.

THE ELEMENT TIER. Exon bases (any transcript in the WBcel235 GFF3) are removed from each element
first, because the human tiers are measured on unknown space, which excludes exons. Per element:
the scored bases left (`measured`) and how many of them pass the cut (`above`). Constrained when
above / measured >= CONSTRAINED_MIN; not constrained otherwise; not measured when no scored base
is left. Conservation is not function, and the tier is not one either.
"""

from __future__ import annotations

import bisect
import gzip
import math
import zlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np

from genomeos.attribution import simpler
from genomeos.attribution.bigwig import BigWig, LeafItem, coalesce, decode_section
from genomeos.attribution.budget import CONSTRAINED_MIN

UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath/ce11"
TRACKS: dict[str, dict[str, Any]] = {
    "phyloP135way": {
        "url": f"{UCSC}/phyloP135way/ce11.phyloP135way.bw",
        "aligned": "134 genomes to C. elegans: 111 nematodes, 22 flatworms, Ciona intestinalis",
        "note": (
            "guide-tree topology from shared k-mers (README: 'not necessarily biologically correct'); "
            "branch lengths for phyloP fitted by phyloFit on 4-fold degenerate sites (track description)"
        ),
    },
    "phyloP26way": {
        "url": f"{UCSC}/phyloP26way/ce11.phyloP26way.bw",
        "aligned": "26 nematode genomes to C. elegans",
        "note": "tree model fitted by phyloFit on 4-fold degenerate sites of the 26-way alignment",
    },
}
#: The primary track: the most Caenorhabditis genomes, where worm non-coding sequence aligns at all.
PRIMARY = "phyloP135way"
SECONDARY = "phyloP26way"
NUCLEAR = ("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX")

FDR_Q = 0.05
#: Histogram of |score| in 0.01 bins (the human cut is quoted to 0.01); |score| >= 20 in the last bin.
HIST_BIN = 0.01
HIST_LIMIT = 20.0
N_BINS = int(round(HIST_LIMIT / HIST_BIN))
SAMPLE_GROUP = 16
SAMPLE_FRACTION = 0.10
CUT_BOOTSTRAP = 1_000
SEED = simpler.SEED
WORKERS = 3
#: Blocks closer than this in the file are fetched in one request.
ELEMENT_MAX_GAP = 16_384
CLASSES = ("cds", "other_exon", "intron", "intergenic")

#: Below this many elements in either tier (in strata holding both), the run stops unscored.
MIN_PER_TIER = 200
PERMUTATIONS = 1_000
BOOTSTRAP = simpler.BOOTSTRAP
#: Measured-length bins: tertile edges of `measured` over the measured population, fixed in the census.
LENGTH_QUANTILES = (1 / 3, 2 / 3)


# --------------------------------------------------------------------------------------------- intervals


def merge(intervals: list[tuple[int, int]]) -> tuple[np.ndarray, np.ndarray]:
    """Sorted, merged half-open intervals as (starts, ends)."""
    out: list[list[int]] = []
    for s, e in sorted(intervals):
        if e <= s:
            continue
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    a = np.array(out, dtype=np.int64).reshape(-1, 2)
    return a[:, 0].copy(), a[:, 1].copy()


def subtract(start: int, end: int, starts: np.ndarray, ends: np.ndarray) -> list[tuple[int, int]]:
    """[start, end) minus merged intervals (starts, ends)."""
    out = []
    pos = start
    i = max(int(np.searchsorted(ends, start, side="right")), 0)
    while i < len(starts) and starts[i] < end:
        s, e = int(starts[i]), int(ends[i])
        if s > pos:
            out.append((pos, min(s, end)))
        pos = max(pos, e)
        i += 1
    if pos < end:
        out.append((pos, end))
    return out


def inside(pos: np.ndarray, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    """Per position, whether it lies in one of the merged intervals."""
    if len(starts) == 0:
        return np.zeros(len(pos), dtype=bool)
    i = np.searchsorted(starts, pos, side="right") - 1
    ok = i >= 0
    out = np.zeros(len(pos), dtype=bool)
    out[ok] = pos[ok] < ends[i[ok]]
    return out


def read_annotation(path: Path = simpler.GFF3) -> dict[str, dict[str, tuple[np.ndarray, np.ndarray]]]:
    """Per nuclear chromosome: merged CDS, merged exons (every transcript) and merged protein-coding
    gene spans, 0-based half-open (the GFF3 is 1-based closed)."""
    raw: dict[str, dict[str, list[tuple[int, int]]]] = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t", 8)
            if len(p) < 9 or p[0] == "MtDNA":
                continue
            kind = p[2]
            if kind == "CDS":
                key = "cds"
            elif kind == "exon":
                key = "exon"
            elif kind == "gene" and "biotype=protein_coding" in p[8]:
                key = "gene"
            else:
                continue
            per = raw.setdefault("chr" + p[0], {"cds": [], "exon": [], "gene": []})
            per[key].append((int(p[3]) - 1, int(p[4])))
    return {c: {k: merge(v) for k, v in d.items()} for c, d in raw.items()}


def classify(pos: np.ndarray, ann: dict[str, tuple[np.ndarray, np.ndarray]]) -> np.ndarray:
    """0 CDS, 1 other exon (UTR, non-coding RNA), 2 intron (in a protein-coding gene), 3 intergenic."""
    cls = np.full(len(pos), 3, dtype=np.int8)
    cls[inside(pos, *ann["gene"])] = 2
    cls[inside(pos, *ann["exon"])] = 1
    cls[inside(pos, *ann["cds"])] = 0
    return cls


# --------------------------------------------------------------------------------------------- the cut


def abs_bins(vals: np.ndarray) -> np.ndarray:
    """|score| to its 0.01 bin. Scores carry three decimals; the small offset keeps a float32 2.27
    (2.2699999...) in the bin that starts at 2.27, matching `vals >= np.float32(2.27)`."""
    idx = np.floor(np.abs(vals.astype(np.float64)) / HIST_BIN + 1e-4).astype(np.int64)
    return np.minimum(idx, N_BINS - 1)


def bh_cut(abs_hist: np.ndarray, q: float = FDR_Q) -> float | None:
    """Benjamini-Hochberg on p = 10^-|score|, from counts of |score| in 0.01 bins.

    BH rejects every p at or below p_(k*), k* the largest k with p_(k) <= q k / m. In score terms:
    the smallest bin edge c at which 10^-c <= q * #(|score| >= c) / m, m every scored base. The
    lower edge of a bin stands for its values (the conservative side: a value in the bin is >= c,
    so its p is <= 10^-c). None when no edge qualifies."""
    h = np.asarray(abs_hist, dtype=np.float64)
    m = h.sum()
    if m <= 0:
        return None
    tail = np.cumsum(h[::-1])[::-1]  # tail[i] = #(|score| >= i * BIN)
    edges = np.arange(len(h)) * HIST_BIN
    ok = (tail > 0) & (10.0 ** (-edges) <= q * tail / m)
    if not ok.any():
        return None
    return round(float(edges[int(np.argmax(ok))]), 2)


def bh_cut_batch(abs_hists: np.ndarray, q: float = FDR_Q) -> np.ndarray:
    """`bh_cut` for each row of a (reps, bins) matrix; NaN where none qualifies."""
    h = np.asarray(abs_hists, dtype=np.float64)
    m = h.sum(axis=1, keepdims=True)
    tail = np.cumsum(h[:, ::-1], axis=1)[:, ::-1]
    edges = np.arange(h.shape[1]) * HIST_BIN
    with np.errstate(divide="ignore", invalid="ignore"):
        ok = (tail > 0) & (10.0 ** (-edges) <= q * tail / m)
    first = np.argmax(ok, axis=1).astype(np.float64)
    first[~ok.any(axis=1)] = np.nan
    return np.round(first * HIST_BIN, 2)


def tier(measured: int, above: int) -> str:
    """The human rule for the regulatory class (budget._rule), on exon-free scored bases."""
    if measured <= 0:
        return "not_measured"
    return "constrained" if above / measured >= CONSTRAINED_MIN else "not_constrained"


# --------------------------------------------------------------------------------------------- reading


def _runs(raw: bytes) -> list[tuple[int, np.ndarray]]:
    """A decompressed section as (start, per-base float32 values)."""
    out = []
    for start, step, vals in decode_section(raw):
        a = np.frombuffer(vals.tobytes(), dtype=np.float32)
        if step != 1:
            a = np.repeat(a, step)
        out.append((start, a))
    return out


def _fetch_groups(bw: BigWig, groups: list[list[LeafItem]], workers: int = WORKERS):
    """Yield (group, blob) with up to `workers` range requests in flight, in order."""

    def get(g: list[LeafItem]) -> bytes:
        base = g[0].offset
        return bw.src.read(base, g[-1].offset + g[-1].size - base)

    with ThreadPoolExecutor(max_workers=workers) as ex:
        yield from zip(groups, ex.map(get, groups), strict=True)


def _sections(bw: BigWig, g: list[LeafItem], blob: bytes):
    base = g[0].offset
    for it in g:
        raw = blob[it.offset - base : it.offset - base + it.size]
        if bw.header.uncompress_buf_size:
            raw = zlib.decompress(raw)
        yield int.from_bytes(raw[:4], "little"), raw


def all_leaves(bw: BigWig, chroms: tuple[str, ...] = NUCLEAR) -> list[LeafItem]:
    """Every data block of the given chromosomes, in file order."""
    out: list[LeafItem] = []
    for c in chroms:
        if c in bw.chroms:
            out += bw.leaf_items(c, [(0, bw.chroms[c][1])])
    uniq = {(it.offset, it.size): it for it in out}
    return sorted(uniq.values(), key=lambda it: it.offset)


def sample_groups(leaves: list[LeafItem], seed: int = SEED) -> list[list[LeafItem]]:
    """Consecutive groups of SAMPLE_GROUP blocks; SAMPLE_FRACTION of them, drawn with the seed."""
    groups = [leaves[i : i + SAMPLE_GROUP] for i in range(0, len(leaves), SAMPLE_GROUP)]
    k = math.ceil(SAMPLE_FRACTION * len(groups))
    rng = np.random.default_rng(seed)
    pick = sorted(rng.choice(len(groups), size=k, replace=False).tolist())
    return [groups[i] for i in pick]


def genome_sample(
    bw: BigWig, ann: dict[str, dict[str, tuple[np.ndarray, np.ndarray]]], seed: int = SEED, progress=None
) -> dict[str, Any]:
    """The score distribution on a uniform sample of blocks: |score| histograms by sign and by class,
    and per 1-Mb window (both signs, every class) for the cut's interval."""
    names = {cid: c for c, (cid, _size) in bw.chroms.items()}
    leaves = all_leaves(bw)
    picked = sample_groups(leaves, seed)
    pos_h = np.zeros((len(CLASSES), N_BINS), dtype=np.int64)
    neg_h = np.zeros((len(CLASSES), N_BINS), dtype=np.int64)
    windows: dict[tuple[str, int], np.ndarray] = {}
    bases = 0
    data_bytes = 0
    for n, (g, blob) in enumerate(_fetch_groups(bw, picked), 1):
        data_bytes += len(blob)
        for cid, raw in _sections(bw, g, blob):
            chrom = names.get(cid)
            if chrom not in ann:
                continue
            for start, vals in _runs(raw):
                pos = start + np.arange(len(vals), dtype=np.int64)
                cls = classify(pos, ann[chrom])
                b = abs_bins(vals)
                positive = vals > 0
                for k in range(len(CLASSES)):
                    sel = cls == k
                    pos_h[k] += np.bincount(b[sel & positive], minlength=N_BINS)
                    neg_h[k] += np.bincount(b[sel & ~positive], minlength=N_BINS)
                win = (chrom, start // simpler.CLUSTER_BP)
                if win not in windows:
                    windows[win] = np.zeros(N_BINS, dtype=np.int64)
                windows[win] += np.bincount(b, minlength=N_BINS)
                bases += len(vals)
        if progress:
            progress(n, len(picked), data_bytes)
    return {
        "blocks_total": len(leaves),
        "groups_total": math.ceil(len(leaves) / SAMPLE_GROUP),
        "groups_sampled": len(picked),
        "blocks_sampled": sum(len(g) for g in picked),
        "bases_sampled": bases,
        "data_bytes": data_bytes,
        "pos": pos_h,
        "neg": neg_h,
        "windows": windows,
    }


def derive_cut(sample: dict[str, Any], seed: int = SEED) -> dict[str, Any]:
    """The 5% FDR cut, its interval over 1-Mb windows, and what passes it, by class."""
    total = sample["pos"].sum(axis=0) + sample["neg"].sum(axis=0)
    cut = bh_cut(total)
    keys = sorted(sample["windows"])
    mat = np.stack([sample["windows"][k] for k in keys]) if keys else np.zeros((0, N_BINS))
    rng = np.random.default_rng(seed)
    w = np.stack(
        [
            np.bincount(rng.integers(0, len(keys), len(keys)), minlength=len(keys))
            for _ in range(CUT_BOOTSTRAP)
        ]
    )
    cuts = bh_cut_batch(w @ mat)
    finite = cuts[~np.isnan(cuts)]
    if len(finite):
        interval = [round(float(np.quantile(finite, 0.025)), 2), round(float(np.quantile(finite, 0.975)), 2)]
    else:
        interval = None
    out: dict[str, Any] = {
        "cut": cut,
        "cut_interval_95": interval,
        "cut_bootstrap": {
            "reps": CUT_BOOTSTRAP,
            "windows_1mb": len(keys),
            "reps_without_a_cut": int(np.isnan(cuts).sum()),
        },
        "scored_bases_sampled": int(total.sum()),
    }
    if cut is None:
        return out
    i = int(round(cut / HIST_BIN))
    by_class = {}
    for k, name in enumerate(CLASSES):
        n = int(sample["pos"][k].sum() + sample["neg"][k].sum())
        a = int(sample["pos"][k][i:].sum())
        by_class[name] = {
            "scored_bases": n,
            "constrained_bases": a,
            "constrained_share": round(a / n, 5) if n else None,
        }
    n_all = int(total.sum())
    a_all = int(sample["pos"].sum(axis=0)[i:].sum())
    nc = sum(by_class[c]["scored_bases"] for c in ("intron", "intergenic"))
    ac = sum(by_class[c]["constrained_bases"] for c in ("intron", "intergenic"))
    out["constrained_share_all_scored"] = round(a_all / n_all, 5) if n_all else None
    out["constrained_share_non_exon"] = round(ac / nc, 5) if nc else None
    out["by_class"] = by_class
    cds, intron = by_class["cds"]["constrained_share"], by_class["intron"]["constrained_share"]
    out["instrument_check"] = {
        "rule": "the cut must place a larger share of CDS bases than of intron bases above it",
        "cds_share": cds,
        "intron_share": intron,
        "passed": bool(cds is not None and intron is not None and cds > intron),
    }
    return out


def element_counts(
    bw: BigWig,
    elements: list[tuple[str, int, int]],
    ann: dict[str, dict[str, tuple[np.ndarray, np.ndarray]]],
    cut: float,
    progress=None,
) -> tuple[list[tuple[int, int]], int]:
    """Per element (chrom, start, end), exon bases removed: (scored bases, bases at or above the cut).
    Returns the counts in the elements' order and the data bytes fetched."""
    cut32 = np.float32(cut)
    measured = [0] * len(elements)
    above = [0] * len(elements)
    data_bytes = 0
    by_chrom: dict[str, list[tuple[int, int, int]]] = {}
    for idx, (c, s, e) in enumerate(elements):
        if c not in ann:
            continue
        for a, b in subtract(s, e, *ann[c]["exon"]):
            by_chrom.setdefault(c, []).append((a, b, idx))
    done = 0
    for c in sorted(by_chrom):
        if c not in bw.chroms:
            continue
        subs = sorted(by_chrom[c])
        starts = [a for a, _, _ in subs]
        maxlen = max(b - a for a, b, _ in subs)
        ms, me = merge([(a, b) for a, b, _ in subs])
        query = list(zip(ms.tolist(), me.tolist(), strict=True))
        cid = bw.chroms[c][0]
        items = bw.leaf_items(c, query)
        for g, blob in _fetch_groups(bw, coalesce(items, max_gap=ELEMENT_MAX_GAP)):
            data_bytes += len(blob)
            for sec_cid, raw in _sections(bw, g, blob):
                if sec_cid != cid:
                    continue
                for start, vals in _runs(raw):
                    end = start + len(vals)
                    j = bisect.bisect_left(starts, start - maxlen)
                    hi = bisect.bisect_left(starts, end)
                    for a, b, idx in subs[j:hi]:
                        lo_, hi_ = max(a, start), min(b, end)
                        if lo_ < hi_:
                            part = vals[lo_ - start : hi_ - start]
                            measured[idx] += hi_ - lo_
                            above[idx] += int(np.count_nonzero(part >= cut32))
        done += 1
        if progress:
            progress(done, len(by_chrom), data_bytes)
    return list(zip(measured, above, strict=True)), data_bytes


# --------------------------------------------------------------------------------------------- the population


def population(
    chroms: dict[str, simpler.Chromosome], elements: list[simpler.Element]
) -> list[dict[str, Any]]:
    """Stage 1's primary population with its domain-rule gene: single-tissue putative enhancers
    whose chromosome has genes and whose domain rule names one. Positions and strands only."""
    rows = []
    for e in elements:
        if e.kind != simpler.ELEMENT_TYPE:
            continue
        ch = chroms.get(e.chrom)
        if ch is None:
            continue
        d, case = simpler.domain_rule(ch, e.mid)
        if d is None or ch.nearest_tss(e.mid) is None:
            continue
        rows.append(
            {
                "chrom": e.chrom,
                "start": e.start,
                "end": e.end,
                "mid": e.mid,
                "tissue": e.tissue,
                "domain_gene": d.gid,
                "domain_case": case,
                "domain_dist": abs(d.tss - e.mid),
            }
        )
    return rows


def length_edges(measured: list[int]) -> list[float]:
    """Tertile edges of measured bases over the measured elements (the registered length bins)."""
    m = np.array([x for x in measured if x > 0], dtype=np.float64)
    return [float(np.quantile(m, q)) for q in LENGTH_QUANTILES] if len(m) else []


def fraction_edges(rows: list[dict[str, Any]]) -> list[float]:
    """Tertile edges of the constrained fraction (above / measured) over the measured elements."""
    f = np.array([r["above"] / r["measured"] for r in rows if r["measured"] > 0], dtype=np.float64)
    return [float(np.quantile(f, q)) for q in LENGTH_QUANTILES] if len(f) else []


def by_fraction_tertile(rows: list[dict[str, Any]], fedges: list[float]) -> list[dict[str, Any]]:
    """Top tertile of the constrained fraction relabelled constrained, bottom tertile not
    constrained, the middle left out: the dose view of the same track, deciding nothing."""
    out = []
    for r in rows:
        if r["measured"] <= 0 or len(fedges) < 2:
            continue
        f = r["above"] / r["measured"]
        if f >= fedges[1] and f > fedges[0]:
            out.append({**r, "tier": "constrained"})
        elif f < fedges[0]:
            out.append({**r, "tier": "not_constrained"})
    return out


def stratum(row: dict[str, Any], edges: list[float]) -> tuple[str, bool, int, int]:
    """Tissue, inside a gene body or not, distance bin of the domain-rule gene, measured-length bin."""
    return (
        row["tissue"],
        row["domain_case"] == "host",
        simpler.distance_bin(row["domain_dist"]),
        bisect.bisect_right(edges, row["length_measured"]),
    )


#: Columns of a stored per-element row: chrom, start, end, then (measured, above) per track.
ELEMENT_COLUMNS = (
    "chrom",
    "start",
    "end",
    f"{PRIMARY}_measured",
    f"{PRIMARY}_above",
    f"{SECONDARY}_measured",
    f"{SECONDARY}_above",
)


def attach_tiers(
    rows: list[dict[str, Any]], tiers_rows: list[list[Any]], track: str = PRIMARY
) -> list[dict[str, Any]]:
    """A copy of the population with one track's stored counts and tier, joined by (chrom, start, end).
    The length used for matching is always the primary track's measured bases, so both tracks are
    compared within the same strata."""
    col = {PRIMARY: (3, 4), SECONDARY: (5, 6)}[track]
    by_key = {(r[0], r[1], r[2]): (r[col[0]], r[col[1]], r[3]) for r in tiers_rows}
    out = []
    for r in rows:
        m, a, length = by_key[(r["chrom"], r["start"], r["end"])]
        out.append({**r, "measured": m, "above": a, "tier": tier(m, a), "length_measured": length})
    return out


def _count(rows: list[dict[str, Any]], pred, t: str) -> int:
    return sum(1 for r in rows if pred(r) and r["tier"] == t)


def census(rows: list[dict[str, Any]], edges: list[float]) -> dict[str, Any]:
    """Tier counts, common support and power in advance, from positions and the track only."""
    measured = [r for r in rows if r["tier"] != "not_measured"]
    groups: dict[tuple, list[int]] = {}
    for r in measured:
        g = groups.setdefault(stratum(r, edges), [0, 0])
        g[0 if r["tier"] == "constrained" else 1] += 1
    support = {k for k, (nc, nu) in groups.items() if nc and nu}
    n_c = sum(groups[k][0] for k in support)
    n_u = sum(groups[k][1] for k in support)

    def split(pred) -> dict[str, int]:
        return {t: _count(measured, pred, t) for t in ("constrained", "not_constrained")}

    # Worst-case Bernoulli variance 1/4 per element; 2.8 SE is 80% power at two-sided 5%, before any
    # inflation for clustering or loss to stratification.
    mde = 2.8 * math.sqrt(0.25 * (1 / n_c + 1 / n_u)) if n_c and n_u else None
    return {
        "elements": len(rows),
        "not_measured": len(rows) - len(measured),
        "constrained": _count(measured, lambda r: True, "constrained"),
        "not_constrained": _count(measured, lambda r: True, "not_constrained"),
        "by_tissue": {t: split(lambda r, t=t: r["tissue"] == t) for t in simpler.TISSUES},
        "inside_gene_body": split(lambda r: r["domain_case"] == "host"),
        "between_genes": split(lambda r: r["domain_case"] != "host"),
        "length_edges_measured_bases": [round(x, 3) for x in edges],
        "strata": len(groups),
        "strata_with_both_tiers": len(support),
        "constrained_in_support": n_c,
        "not_constrained_in_support": n_u,
        "mde_points_unclustered": round(100 * mde, 3) if mde is not None else None,
        "fraction_tertile_edges": [round(x, 6) for x in fraction_edges(rows)],
        "stop_rule_met": n_c < MIN_PER_TIER or n_u < MIN_PER_TIER,
    }


# --------------------------------------------------------------------------------------------- scoring


def _arrays(rows: list[dict[str, Any]], edges: list[float]):
    keys = sorted({stratum(r, edges) for r in rows}, key=repr)
    sid = {k: i for i, k in enumerate(keys)}
    clusters = sorted({(r["chrom"], r["mid"] // simpler.CLUSTER_BP) for r in rows})
    cid = {k: i for i, k in enumerate(clusters)}
    s = np.array([sid[stratum(r, edges)] for r in rows], dtype=np.int64)
    c = np.array([cid[(r["chrom"], r["mid"] // simpler.CLUSTER_BP)] for r in rows], dtype=np.int64)
    t = np.array([r["tier"] == "constrained" for r in rows], dtype=bool)
    y = np.array([r["c_domain"] for r in rows], dtype=np.float64)
    return s, c, t, y, len(keys), len(clusters)


def _std_diff(nc, kc, nu, ku) -> tuple[float, float, float, float, float]:
    """Standardised to the constrained elements' strata: (difference, constrained rate, matched
    not-constrained rate, constrained n, not-constrained n), strata holding both tiers only."""
    nc, kc, nu, ku = (np.asarray(x, dtype=np.float64) for x in (nc, kc, nu, ku))
    ok = (nc > 0) & (nu > 0)
    if not ok.any():
        return math.nan, math.nan, math.nan, 0.0, 0.0
    w = nc[ok]
    rc = kc[ok] / nc[ok]
    ru = ku[ok] / nu[ok]
    tot = w.sum()
    return (
        float((w * (rc - ru)).sum() / tot),
        float((w * rc).sum() / tot),
        float((w * ru).sum() / tot),
        float(nc[ok].sum()),
        float(nu[ok].sum()),
    )


def _permuted_labels(s: np.ndarray, n_per: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """The constrained label shuffled within each stratum: n_per[k] elements of stratum k drawn."""
    order = np.lexsort((rng.random(len(s)), s))
    first = np.searchsorted(s[order], np.arange(len(n_per)))
    rank = np.empty(len(s), dtype=np.int64)
    rank[order] = np.arange(len(s)) - first[s[order]]
    return rank < n_per[s]


def matched_difference(
    rows: list[dict[str, Any]],
    edges: list[float],
    reps: int = BOOTSTRAP,
    perms: int = PERMUTATIONS,
    seed: int = SEED,
) -> dict[str, Any]:
    """Concordance of the domain-rule gene, constrained minus not constrained, matched by
    stratification; 95% percentile interval over 1-Mb windows; within-stratum permutation p."""
    rows = [r for r in rows if r["tier"] != "not_measured"]
    if not rows:
        return {"n": 0, "delta_points": None, "delta_ci95_points": None}
    s, c, t, y, n_s, n_cl = _arrays(rows, edges)
    flat = c * n_s + s
    size = n_cl * n_s

    def per_cluster(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        n = np.bincount(flat[mask], minlength=size).reshape(n_cl, n_s).astype(np.float64)
        k = np.bincount(flat[mask], weights=y[mask], minlength=size).reshape(n_cl, n_s)
        return n, k

    nc, kc = per_cluster(t)
    nu, ku = per_cluster(~t)
    point, rate_c, rate_u, n_cs, n_us = _std_diff(nc.sum(0), kc.sum(0), nu.sum(0), ku.sum(0))
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(reps):
        w = np.bincount(rng.integers(0, n_cl, n_cl), minlength=n_cl).astype(np.float64)
        d = _std_diff(w @ nc, w @ kc, w @ nu, w @ ku)[0]
        if not math.isnan(d):
            boots.append(d)
    boots.sort()
    n_per = np.bincount(s[t], minlength=n_s)
    prng = np.random.default_rng(seed + 1)
    ge = 0
    for _ in range(perms):
        tp = _permuted_labels(s, n_per, prng)
        d = _std_diff(
            np.bincount(s[tp], minlength=n_s),
            np.bincount(s[tp], weights=y[tp], minlength=n_s),
            np.bincount(s[~tp], minlength=n_s),
            np.bincount(s[~tp], weights=y[~tp], minlength=n_s),
        )[0]
        if d >= point - 1e-12:
            ge += 1
    has = n_cs > 0 and not math.isnan(point)
    raw_c = [r for r in rows if r["tier"] == "constrained"]
    raw_u = [r for r in rows if r["tier"] != "constrained"]
    return {
        "n": len(rows),
        "constrained_in_support": int(n_cs),
        "not_constrained_in_support": int(n_us),
        "constrained_concordant_in_support": int(round(rate_c * n_cs)) if has else None,
        "constrained_rate_pct": round(100 * rate_c, 3) if has else None,
        "matched_not_constrained_rate_pct": round(100 * rate_u, 3) if has else None,
        "delta_points": round(100 * point, 3) if has else None,
        "delta_ci95_points": (
            [
                round(100 * boots[int(0.025 * len(boots))], 3),
                round(100 * boots[int(0.975 * len(boots)) - 1], 3),
            ]
            if boots
            else None
        ),
        "bootstrap": {"reps": reps, "clusters_1mb": n_cl, "reps_without_support": reps - len(boots)},
        "permutation": {
            "perms": perms,
            "seed": seed + 1,
            "at_or_above_observed": ge,
            "p_one_sided": round((1 + ge) / (1 + perms), 5),
        },
        "unmatched": {
            "constrained": len(raw_c),
            "constrained_concordant": int(sum(r["c_domain"] for r in raw_c)),
            "not_constrained": len(raw_u),
            "not_constrained_concordant": int(sum(r["c_domain"] for r in raw_u)),
        },
    }


def outcome_label(lo: float, hi: float) -> str:
    if lo > 0:
        return "CONSTRAINT_ADDS"
    if hi < 0:
        return "CONSTRAINT_BELOW"
    return "NO_DIFFERENCE_DETECTED"


def score(
    rows: list[dict[str, Any]],
    secondary: list[dict[str, Any]],
    expr: dict[str, tuple[float, ...]],
    edges: list[float],
    cen: dict[str, Any],
    instrument_passed: bool,
) -> dict[str, Any]:
    """The registered run. `rows` carry the primary tiers, `secondary` the same elements with the
    secondary track's; `expr` is the endpoint, read only here."""
    if not instrument_passed:
        res: dict[str, Any] = {"outcome": "INSTRUMENT_FAILED"}
        res["reading"] = reading(res)
        return res
    if cen["stop_rule_met"]:
        res = {"outcome": "UNDERPOWERED"}
        res["reading"] = reading(res)
        return res
    for rr in (rows, secondary):
        for r in rr:
            r["c_domain"] = int(simpler.concordant(expr.get(r["domain_gene"]), r["tissue"]))
    primary = matched_difference(rows, edges)
    label = outcome_label(*primary["delta_ci95_points"])
    falsifier = "FALSIFIED" if primary["permutation"]["p_one_sided"] >= simpler.ALPHA else "NOT_FALSIFIED"
    strata = {
        "inside_gene_body": matched_difference([r for r in rows if r["domain_case"] == "host"], edges),
        "between_genes": matched_difference([r for r in rows if r["domain_case"] != "host"], edges),
    }
    for t in simpler.TISSUES:
        strata[t] = matched_difference([r for r in rows if r["tissue"] == t], edges)
    strata[f"track_{SECONDARY}"] = matched_difference(secondary, edges)
    strata["fraction_top_vs_bottom_tertile"] = matched_difference(
        by_fraction_tertile(rows, cen["fraction_tertile_edges"]), edges
    )
    res = {
        "outcome": label,
        "falsifier": falsifier,
        "primary": primary,
        "strata_beside": strata,
        "domain_genes_absent_from_rna_sheet": sum(1 for r in rows if r["domain_gene"] not in expr),
    }
    res["reading"] = reading(res)
    return res


def reading(res: dict[str, Any]) -> str:
    """The registered outcome wording, filled in. Templates: REGISTRATION['outcome_wording']."""
    w = REGISTRATION["outcome_wording"]
    if res["outcome"] in ("INSTRUMENT_FAILED", "UNDERPOWERED"):
        return w[res["outcome"]]
    p = res["primary"]
    lo, hi = p["delta_ci95_points"]
    vals = {
        "n": f"{p['n']:,}",
        "nc": f"{p['constrained_in_support']:,}",
        "kc": f"{p['constrained_concordant_in_support']:,}",
        "nu": f"{p['not_constrained_in_support']:,}",
        "rc": f"{p['constrained_rate_pct']:.1f}%",
        "ru": f"{p['matched_not_constrained_rate_pct']:.1f}%",
        "d": f"{p['delta_points']:+.2f}",
        "lo": f"{lo:+.2f}",
        "hi": f"{hi:+.2f}",
        "pp": f"{p['permutation']['p_one_sided']:.3g}",
        "falsifier": res["falsifier"],
    }
    text = w["common"].format(**vals) + " " + w[res["outcome"]]
    between = res["strata_beside"]["between_genes"]
    if res["outcome"] == "CONSTRAINT_ADDS" and (
        not between.get("delta_ci95_points") or between["delta_ci95_points"][0] <= 0
    ):
        text += " " + w["between_genes_caveat"]
    return text + " " + w["always"]


REGISTRATION: dict[str, Any] = {
    "claim": (
        "Among C. elegans putative enhancers accessible in exactly one of five tissues, the gene the "
        "domain rule attributes is tissue-concordant more often for elements in the constrained tier "
        "than for elements not in it, at the same tissue, position class, distance to that gene and "
        "measured length: constraint carries tissue information beyond the domain rule's position."
    ),
    "constraint_cut": (
        "Per track, Zoonomia's own procedure on the track's own scores: p = 10^-|phyloP| for every "
        "scored base, Benjamini-Hochberg at q = 0.05 over all scored bases of both signs; constrained "
        "base = positive score at or above the cut. Distribution from a uniform sample of data blocks "
        "(consecutive groups of 16 blocks of chrI-V and chrX, 10% of the groups, seed 20261010), "
        "|phyloP| in 0.01 bins; cut interval: 1,000 bootstrap resamples of 1-Mb windows of the sample. "
        "No UCSC page or publication for either ce11 track gives a cut-off."
    ),
    "tracks": (
        "Primary: UCSC ce11 phyloP135way (the most Caenorhabditis genomes, where worm non-coding "
        "sequence aligns at all). Beside, deciding nothing: phyloP26way with its own cut."
    ),
    "element_tier": (
        "Exon bases of any WBcel235 transcript removed first (the human tiers are measured on unknown "
        "space, which excludes exons). Constrained: at least 5% of the remaining scored bases at or "
        "above the cut (budget.CONSTRAINED_MIN, the human rule for the regulatory class). Not "
        "constrained: below 5%. Not measured: no scored base left; excluded and counted."
    ),
    "population": (
        "Stage 1's primary population (simpler_worm_registration): Serizay et al. 2020 putative "
        "enhancers accessible in exactly one of the five tissues by the source's call, with the domain "
        "rule's gene."
    ),
    "endpoint": (
        "Stage 1's endpoint, unchanged: the domain-rule gene's TPM in the element's tissue >= 5 and "
        "strictly higher than in each of the other four; a gene absent from RNA_metrics is not concordant."
    ),
    "matching": (
        "Strata: tissue (5) x inside a gene body or between genes (2) x distance bin of the domain-rule "
        "gene's 5' end (stage 1's 8 bins) x measured-length tertile (edges fixed in the census, from "
        "the primary track's measured bases). Only strata holding both tiers count."
    ),
    "statistic": (
        "delta = concordance(constrained) - concordance(not constrained), standardised to the "
        "constrained elements' strata (each stratum weighted by its constrained count), in points; 95% "
        "percentile interval over 1-Mb windows, 2,000 resamples, seed 20261010."
    ),
    "falsifier": (
        "The claim is falsified if delta is not above its within-stratum permutation null: the "
        "constrained label shuffled within each stratum 1,000 times (seed 20261011), one-sided "
        "p = (1 + permutations with delta >= observed) / 1,001; FALSIFIED if p >= 0.05."
    ),
    "decision": (
        "CONSTRAINT_ADDS if the interval's lower bound > 0; CONSTRAINT_BELOW if its upper bound < 0; "
        "NO_DIFFERENCE_DETECTED otherwise."
    ),
    "instrument_check": (
        "Before any scoring: in the genome sample the primary cut must place a larger share of CDS "
        "bases than of intron bases above it; otherwise the track read is declared broken "
        "(INSTRUMENT_FAILED) and nothing is scored."
    ),
    "strata_beside": (
        "Deciding nothing: inside gene bodies and between genes apart; each tissue; the phyloP26way "
        "tiers on the same elements; the dose view (the primary's constrained fraction, top tertile "
        "against bottom tertile at edges fixed in the census, matched the same way), named because "
        "the 5% rule puts most measured elements in the constrained tier when 7.2% of the worm's "
        "non-exon scored bases pass the cut; the unmatched counts. If the primary is CONSTRAINT_ADDS "
        "and the between-genes lower bound is <= 0, the between-genes caveat is appended."
    ),
    "confound": (
        "Stage 1 named transcription: an element inside a gene body can be open where its host is "
        "transcribed for that reason alone, and most worm transcripts are trans-spliced, so an element "
        "upstream of an annotated 5' end can lie in a transcribed outron. Constraint shares that "
        "confound: transcribed sequence can be constrained for splicing, trans-splicing or RNA-level "
        "reasons. Exon bases are removed and the position class is a stratum, but introns and outrons "
        "are not removed (outrons are not annotated). A difference therefore does not separate "
        "regulatory constraint from transcription-linked constraint."
    ),
    "stop_rule": (
        "One run. Fewer than 200 elements of either tier in strata holding both: stop unscored "
        "(UNDERPOWERED). Nothing (cut, tier rule, strata, endpoint) changes after the RNA sheet is "
        "read; another definition is a new registration committed before it is run."
    ),
    "outcome_wording": {
        "common": (
            "On {n} C. elegans putative enhancers accessible in exactly one of five tissues (Serizay et "
            "al. 2020) with constraint measured on their exon-free bases (UCSC ce11 phyloP135way, cut "
            "re-derived at 5% FDR by Zoonomia's procedure), the domain-rule gene is expressed at 5 TPM "
            "or more and highest in the element's own tissue for {kc} of {nc} constrained elements ({rc}), "
            "against {ru} for the {nu} not constrained at the same tissue, position class, distance "
            "and measured length; difference {d} points, 95% {lo} to {hi} over 1-Mb clusters "
            "(within-stratum permutation p {pp}; {falsifier})."
        ),
        "CONSTRAINT_ADDS": (
            "On this endpoint constraint carries tissue information beyond the domain rule's position."
        ),
        "NO_DIFFERENCE_DETECTED": (
            "No difference was detected: on this endpoint constraint adds no detectable tissue "
            "information to the domain rule's position. This is not evidence that constrained "
            "elements lack function."
        ),
        "CONSTRAINT_BELOW": (
            "Constrained elements' genes are less often tissue-concordant than matched unconstrained "
            "ones' on this endpoint; one reading named in advance is that constrained elements sit "
            "near broadly expressed genes, which an endpoint asking for the highest tissue counts "
            "against."
        ),
        "between_genes_caveat": (
            "The difference is not carried by the elements between genes, where transcription by a "
            "host gene is not the explanation."
        ),
        "always": (
            "Conservation and accessibility are not function, and this run does not separate "
            "regulatory constraint from transcription-linked constraint."
        ),
        "INSTRUMENT_FAILED": (
            "The cut placed no larger share of CDS than of intron bases above it, so the track read is "
            "declared broken and nothing was scored."
        ),
        "UNDERPOWERED": "Fewer than 200 elements of one tier in strata holding both tiers: stopped unscored.",
    },
}
