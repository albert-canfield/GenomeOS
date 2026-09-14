# SPDX-License-Identifier: AGPL-3.0-or-later
"""The lexicon's two axes at base resolution, and the JASPAR scan at the calibrated threshold.

The first lexicon (attribution/lexicon.py, 2026-09-13) asked whether its units separate into syntax
(held across mammals and among people) and value slots (held across mammals, free among people),
and found nothing. It named why the answer might be soft: conservation was membership in UCSC's
100-vertebrate conserved elements, the human axis was one Gnocchi score per UNKNOWN block, and the
motif sites were the 0.85 hits the motif runs had recorded inside selected elements. This module
replaces all three and asks again.

*Mammals.* Zoonomia phyloP over 241 placental mammals, every base of the chromosome, streamed from
UCSC's bigWig by range requests through the reader the budget uses (attribution/bigwig.py) and
distilled to one signed byte per base, placed so that the 5% FDR threshold of 2.27 falls on a bin
edge. The byte track is the only thing kept (data/knowledge/lexicon, git-ignored); the float
sections are decoded and dropped.

*People.* Three measures exist and they disagree (genomeos-h1, 2026-09-14: Gnocchi against the
89-assembly panel is rank correlation -0.08, and per-kilobase panel counts are 15.6 times
overdispersed). Per base, a unit can only be read by a per-base measure, so the axis here is the
HPRC panel's recurring substitutions (a column where two or more of 89 assemblies carry a base other
than the commonest), read from h1's distilled store through its own reader, with the base's
trinucleotide as the mutation-rate control. Gnocchi per kilobase is joined too, for the record and
the comparison, not as the axis. The overdispersion is handled where it arises: every unit's test
uses a variance built from its own occurrences (a sandwich estimate), not a Poisson one, and every
family of tests is calibrated on shifted copies of the units, which keep their clustering.

*Motif sites.* Every JASPAR 2026 CORE vertebrate profile scanned over the whole chromosome at 0.95
of the matrix range on both strands (motifs.py's calibration: at 0.85 some profile covers 99.9% of
an element's bases; at 0.95, 65%). The scan is indexed: one sort of the chromosome's 8-mers, each
profile's feasible core k-mers looked up as ranges of that order, candidates scored in blocks.
"""

from __future__ import annotations

import gzip
import json
import math
import time
import zlib
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np

KNOWLEDGE = Path("data/knowledge/lexicon")
PHYLOP_THRESHOLD = 2.27
PHYLOP_MISSING = -128
# a byte per base: q = floor((p - 2.27) * 10) + 23, so q >= 23 is exactly p >= 2.27
PHYLOP_SCALE = 10
PHYLOP_ANCHOR = 23
PHYLOP_BINS = (-128, -9, 0, 10, 23, 40, 128)  # edges on q: missing | < -1.1 | < 1.3 | < 2.27 | < 4 | >= 4
MOTIF_RELATIVE = 0.95
SCAN_BLOCK = 4_000_000  # candidate windows scored per numpy block
CODE = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    CODE[_b] = _i
    CODE[_b + 32] = _i


# ---- the mammal axis --------------------------------------------------------------------------


def quantise_phylop(values: np.ndarray) -> np.ndarray:
    """phyloP floats to signed bytes with 2.27 on a bin edge; NaN becomes the missing code."""
    q = np.floor((values - PHYLOP_THRESHOLD) * PHYLOP_SCALE) + PHYLOP_ANCHOR
    q = np.clip(q, -127, 127)
    q[np.isnan(values)] = PHYLOP_MISSING
    return q.astype(np.int8)


def dequantise_phylop(q: np.ndarray) -> np.ndarray:
    """The bin midpoint back as a float; missing bases are NaN."""
    out = (q.astype(np.float32) - PHYLOP_ANCHOR + 0.5) / PHYLOP_SCALE + PHYLOP_THRESHOLD
    out[q == PHYLOP_MISSING] = np.nan
    return out


def track_bytes(url: str, chrom: str, length: int, progress=None) -> tuple[np.ndarray, dict[str, Any]]:
    """A whole chromosome of a bigWig as float32 per base (NaN where the track has no value),
    streamed in coalesced range requests; returns the array and what it cost."""
    from genomeos.attribution.bigwig import BigWig, coalesce, decode_section

    t0 = time.time()
    bw = BigWig(url)
    out = np.full(length, np.nan, dtype=np.float32)
    try:
        cid = bw.chroms[chrom][0]
        items = bw.leaf_items(chrom, [(0, length)])
        if progress:
            progress(
                f"{chrom}: {len(items)} data blocks, {sum(i.size for i in items) / 1e6:.0f} MB to stream"
            )
        for n, g in enumerate(coalesce(items)):
            base = g[0].offset
            blob = bw.src.read(base, g[-1].offset + g[-1].size - base)
            for it in g:
                raw = blob[it.offset - base : it.offset - base + it.size]
                if bw.header.uncompress_buf_size:
                    raw = zlib.decompress(raw)
                if int.from_bytes(raw[:4], "little") != cid:
                    continue
                for start, step, vals in decode_section(raw):
                    v = np.frombuffer(vals, dtype=np.float32)
                    if step == 1:
                        a, b = max(0, start), min(length, start + len(v))
                        if b > a:
                            out[a:b] = v[a - start : b - start]
                    else:
                        for k, x in enumerate(v):
                            a = start + k * step
                            out[max(0, a) : min(length, a + step)] = x
            if progress and n % 20 == 0:
                progress(f"  {bw.src.bytes_fetched / 1e6:.0f} MB fetched")
        cost = {
            "requests": bw.src.requests,
            "mb_fetched": round(bw.src.bytes_fetched / 1e6, 1),
            "seconds": round(time.time() - t0, 1),
        }
    finally:
        bw.close()
    return out, cost


def phylop_track(chrom: str, length: int, progress=None, directory: Path = KNOWLEDGE) -> np.ndarray:
    """The distilled phyloP byte track of one chromosome, streamed once and kept as bytes."""
    from genomeos.attribution.constraint import PHYLOP_241_URL

    path = directory / f"phylop241_{chrom}.i8.gz"
    if path.exists():
        with gzip.open(path, "rb") as fh:
            return np.frombuffer(fh.read(), dtype=np.int8).copy()
    values, cost = track_bytes(PHYLOP_241_URL, chrom, length, progress)
    q = quantise_phylop(values)
    del values
    directory.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wb", compresslevel=6) as fh:
        fh.write(q.tobytes())
    (directory / f"phylop241_{chrom}.cost.txt").write_text(f"{cost}\n")
    return q


def gnocchi_track(chrom: str, length: int, progress=None, directory: Path = KNOWLEDGE) -> np.ndarray:
    """Gnocchi Z per base (its kilobase value repeated), kept as a float16 track."""
    from genomeos.attribution.variation import GNOCCHI_URL

    path = directory / f"gnocchi_{chrom}.f16.gz"
    if path.exists():
        with gzip.open(path, "rb") as fh:
            return np.frombuffer(fh.read(), dtype=np.float16).astype(np.float32)
    values, _ = track_bytes(GNOCCHI_URL, chrom, length, progress)
    directory.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wb", compresslevel=6) as fh:
        fh.write(values.astype(np.float16).tobytes())
    return values


# ---- the human axis ---------------------------------------------------------------------------


def panel_tracks(chrom: str, length: int, min_present: int = 80, min_minor: int = 2):
    """From the HPRC panel store: per base, whether the base is informative (inside an alignment
    block where at least `min_present` of the assemblies are present) and whether it is a recurring
    substitution (at least `min_minor` assemblies away from the commonest state, one of them a base).
    Returns (informative, variable) as uint8 arrays, or None when the chromosome has no store."""
    from genomeos.attribution.human_panel import CACHE, Panel

    directory = CACHE / chrom
    if not (directory / "meta.json").exists():
        return None
    panel = Panel(directory)
    informative = np.zeros(length, dtype=np.uint8)
    for s, e, st in zip(panel.block_start, panel.block_end, panel.block_status, strict=True):
        if st.count(b"s") >= min_present:
            informative[s:e] = 1
    variable = np.zeros(length, dtype=np.uint8)
    for k in range(len(panel.site_pos)):
        c = panel.counts(k)
        if c.minor >= min_minor and c.alt:
            variable[c.pos] = 1
    variable &= informative
    return informative, variable


def trinucleotide_class(codes: np.ndarray) -> np.ndarray:
    """Per base the strand-collapsed trinucleotide (32 classes, the middle base A or C after
    reverse-complementing G and T centres); 255 where a neighbour is not A, C, G or T."""
    n = len(codes)
    out = np.full(n, 255, dtype=np.uint8)
    if n < 3:
        return out
    left, mid, right = codes[:-2].astype(np.int32), codes[1:-1].astype(np.int32), codes[2:].astype(np.int32)
    ok = (left < 4) & (mid < 4) & (right < 4)
    tri = left * 16 + mid * 4 + right
    rc = (3 - right) * 16 + (3 - mid) * 4 + (3 - left)
    canon = np.where(mid < 2, tri, rc)  # centre A or C kept, G or T flipped
    # 32 classes: centre (0/1) * 16 + left * 4 + right
    cls = (canon // 4 % 4) * 16 + (canon // 16) * 4 + canon % 4
    out[1:-1] = np.where(ok, cls, 255).astype(np.uint8)
    return out


# ---- the motif scan -----------------------------------------------------------------------------


def encode(seq: str) -> np.ndarray:
    return CODE[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]


def reverse_complement_motif(m):
    """The same profile read on the other strand, as a motif scored on forward coordinates."""
    from genomeos.genome.motifs import Motif

    counts = [list(reversed(m.counts[3 - b])) for b in range(4)]
    return Motif(m.id, m.name, counts).prepare(MOTIF_RELATIVE)


def kmer_order(codes: np.ndarray, k: int = 8) -> tuple[np.ndarray, np.ndarray]:
    """Positions sorted by their k-mer code, and the offsets of each code in that order. Windows
    touching a non-ACGT base are left out."""
    n = len(codes) - k + 1
    kcode = np.zeros(n, dtype=np.int32)
    bad = np.zeros(n, dtype=bool)
    for j in range(k):
        c = codes[j : j + n]
        kcode = (kcode << 2) | np.minimum(c, 3).astype(np.int32)
        bad |= c > 3
    kcode[bad] = 1 << (2 * k)  # past every real code
    order = np.argsort(kcode, kind="stable").astype(np.int32)
    counts = np.bincount(kcode, minlength=(1 << (2 * k)) + 1)
    offsets = np.concatenate(([0], np.cumsum(counts)))
    return order, offsets


def scan_profile(m, codes: np.ndarray, order: np.ndarray, offsets: np.ndarray, k: int = 8) -> np.ndarray:
    """Forward start positions where one strand-specific profile scores at or above its threshold."""
    w = m.width
    pwm = np.full((w, 5), -1e9, dtype=np.float32)
    pwm[:, :4] = np.asarray(m.pwm, dtype=np.float32)
    cores = m.feasible_cores()
    shift = 2 * (k - m.core_k)
    chunks = [order[offsets[c << shift] : offsets[(c + 1) << shift]] for c in cores]
    if not chunks:
        return np.zeros(0, dtype=np.int64)
    cand = np.concatenate(chunks).astype(np.int64) - m.core_offset
    cand = cand[(cand >= 0) & (cand + w <= len(codes))]
    hits = []
    cols = np.arange(w)
    for i in range(0, len(cand), max(1, SCAN_BLOCK // w)):
        starts = cand[i : i + SCAN_BLOCK // w]
        win = codes[starts[:, None] + cols]
        score = pwm[cols, win].sum(axis=1)
        hits.append(starts[score >= m.threshold])
    return np.sort(np.concatenate(hits)) if hits else np.zeros(0, dtype=np.int64)


def decoy_motifs(motifs: list, seed: int = 20260914) -> list:
    """Every profile with its columns permuted: the same base composition and information per column,
    none of the motif's order. A site effect that decoys share belongs to the letters, not the motif."""
    import random

    from genomeos.genome.motifs import Motif

    rng = random.Random(seed)
    out = []
    for m in motifs:
        order = list(range(m.width))
        for _ in range(10):
            rng.shuffle(order)
            if order != sorted(order):
                break
        counts = [[m.counts[b][j] for j in order] for b in range(4)]
        out.append(Motif(m.id, m.name, counts).prepare(MOTIF_RELATIVE))
    return out


def dedupe_sites(
    sites: dict[str, list[tuple[int, int]]],
) -> tuple[dict[str, list[tuple[int, int]]], dict[str, list[str]]]:
    """Factors whose profiles give the identical site set are one unit (JASPAR lists several names for
    one matrix family); the kept name carries the others as aliases."""
    seen: dict[tuple, str] = {}
    kept: dict[str, list[tuple[int, int]]] = {}
    aliases: dict[str, list[str]] = {}
    for name in sorted(sites):
        key = (len(sites[name]), tuple(sites[name][:50]), tuple(sites[name][-50:]))
        if key in seen:
            aliases[seen[key]].append(name)
            continue
        seen[key] = name
        kept[name] = sites[name]
        aliases[name] = []
    return kept, aliases


def scan_chromosome(
    seq: str, progress=None, motifs: list | None = None, index=None
) -> dict[str, list[tuple[int, int]]]:
    """Every factor's sites on both strands at MOTIF_RELATIVE, overlapping sites of one factor merged."""
    from genomeos.genome.motifs import load_motifs

    t0 = time.time()
    codes = encode(seq)
    order, offsets = index if index is not None else kmer_order(codes)
    motifs = motifs if motifs is not None else load_motifs(relative=MOTIF_RELATIVE)
    by_factor: dict[str, list[np.ndarray]] = {}
    widths: dict[str, list[int]] = {}
    for n, m in enumerate(motifs):
        name = m.name.upper()
        for strand_motif in (m, reverse_complement_motif(m)):
            starts = scan_profile(strand_motif, codes, order, offsets)
            if len(starts):
                by_factor.setdefault(name, []).append(starts)
                widths.setdefault(name, []).append(strand_motif.width)
        if progress and n % 200 == 0:
            progress(f"  {n} of {len(motifs)} profiles scanned, {time.time() - t0:.0f} s")
    out: dict[str, list[tuple[int, int]]] = {}
    for name, arrays in by_factor.items():
        ivs = sorted((int(s), int(s) + w) for arr, w in zip(arrays, widths[name], strict=True) for s in arr)
        merged: list[tuple[int, int]] = []
        for s, e in ivs:
            if merged and s <= merged[-1][1]:
                if e > merged[-1][1]:
                    merged[-1] = (merged[-1][0], e)
            else:
                merged.append((s, e))
        out[name] = merged
    return out


# ---- per-base arrays and the expectation tables --------------------------------------------------

CONTEXT_GROUP = {
    "cds": 0,
    "utr_or_exon": 0,
    "exon_noncoding": 0,
    "promoter_like": 1,
    "enhancer_like": 1,
    "ctcf_only": 1,
    "unknown_regulatory": 1,
    "unknown_structural": 2,
    "unknown_fossil": 2,
    "unknown_constrained": 2,
    "unknown_neutral": 2,
    "intron": 3,
    "other": 3,
}
PSEUDO_BASES = 500  # a thin expectation cell borrows this many bases' worth of its marginal
MIN_OCC = 20
FLANK_GAP = 10  # bases between an occurrence and its flanks
FLANK_MAX = 1_000  # a flank is as long as its occurrence, up to this
MIN_EXPECTED_VARIANTS = 5.0


def phylop_bin(q: np.ndarray) -> np.ndarray:
    """0 missing, then < -0.93, < 0, < 1, < 2.27, < 4, >= 4 (edges on the byte scale)."""
    out = np.digitize(q, np.array(PHYLOP_BINS[1:-1])).astype(np.uint8) + 1
    out[q == PHYLOP_MISSING] = 0
    return out


def _rates(
    cells: np.ndarray, mask: np.ndarray, hits: np.ndarray, size: int, marginal_cells=None, marginal_size=0
):
    """Share of `hits` among masked bases per cell, shrunk toward a marginal cell when thin; the
    last slot (index `size`) is the invalid cell with rate 0."""
    n = np.bincount(cells[mask], minlength=size + 1)[: size + 1].astype(np.float64)
    k = np.bincount(cells[mask], weights=hits[mask].astype(np.float64), minlength=size + 1)[: size + 1]
    if marginal_cells is None:
        rate = np.divide(k, n, out=np.zeros_like(k), where=n > 0)
    else:
        mn = np.bincount(marginal_cells[mask], minlength=marginal_size + 1)[: marginal_size + 1].astype(float)
        mk = np.bincount(
            marginal_cells[mask], weights=hits[mask].astype(np.float64), minlength=marginal_size + 1
        )[: marginal_size + 1]
        mrate = np.divide(mk, mn, out=np.zeros_like(mk), where=mn > 0)
        # every cell's marginal: the first cell index that maps to it
        first = np.full(size + 1, marginal_size, dtype=np.int64)
        seen = np.zeros(size + 1, dtype=bool)
        valid = mask & (cells < size)
        c, m = cells[valid], marginal_cells[valid]
        uniq, where = np.unique(c, return_index=True)
        first[uniq] = m[where]
        seen[uniq] = True
        prior = mrate[first]
        rate = (k + PSEUDO_BASES * prior) / (n + PSEUDO_BASES)
        rate[~seen] = 0.0
    rate[size] = 0.0
    return rate.astype(np.float32)


def base_arrays(seq: str, maps, windows, q: np.ndarray, human, gnocchi: np.ndarray | None) -> dict[str, Any]:
    """Every per-base array and expectation table the unit tests read."""
    from genomeos.attribution.lexicon import CONTEXTS, REPEAT_BINS

    length = maps.length
    ctx = np.frombuffer(maps.context, dtype=np.uint8)
    stratum = np.full(length, 255, dtype=np.uint8)
    for r in windows.rows:
        stratum[r["start"] : r["end"]] = r["stratum"][0] * REPEAT_BINS + r["stratum"][1]
    n_strata = 5 * REPEAT_BINS
    ok_st = stratum < n_strata
    gcbin = np.where(ok_st, stratum // REPEAT_BINS, 0).astype(np.int64)
    tri = trinucleotide_class(encode(seq))
    group = np.array([CONTEXT_GROUP[c] for c in CONTEXTS], dtype=np.int64)[ctx]
    pb = phylop_bin(q).astype(np.int64)
    measured = q != PHYLOP_MISSING
    constrained = q >= PHYLOP_ANCHOR
    n_ctx = len(CONTEXTS)
    # mammals: the constrained share of the same context in the same stratum, and of the stratum alone
    size_m = n_ctx * n_strata
    cell_m = np.where(ok_st, ctx.astype(np.int64) * n_strata + stratum, size_m)
    rate_m = _rates(cell_m, measured, constrained, size_m)
    cell_ms = np.where(ok_st, stratum.astype(np.int64), n_strata)
    rate_ms = _rates(cell_ms, measured, constrained, n_strata)
    out: dict[str, Any] = {
        "length": length,
        "ctx": ctx,
        "q": q,
        "measured": measured,
        "constrained": constrained,
        "cell_m": cell_m.astype(np.uint16),
        "rate_m": rate_m,
        "cell_ms": cell_ms.astype(np.uint8),
        "rate_ms": rate_ms,
        "gnocchi": gnocchi,
        "human": human is not None,
    }
    if human is not None:
        informative, variable = human
        ok_tri = tri < 32
        ok = (informative == 1) & ok_st & ok_tri
        t = tri.astype(np.int64)
        # E0: trinucleotide x context x GC bin; the marginal is trinucleotide x GC bin
        size0, size_t = 32 * n_ctx * 5, 32 * 5
        cell_t = np.where(ok, t * 5 + gcbin, size_t)
        cell0 = np.where(ok, (t * n_ctx + ctx) * 5 + gcbin, size0)
        var = variable == 1
        rate0 = _rates(cell0, ok, var, size0, cell_t, size_t)
        rate_t = _rates(cell_t, ok, var, size_t)
        # E1: trinucleotide x phyloP bin x context group x GC bin; the marginal is trinucleotide x phyloP bin
        size1, size_tp = 32 * 7 * 4 * 5, 32 * 7
        cell_tp = np.where(ok, t * 7 + pb, size_tp)
        cell1 = np.where(ok, ((t * 7 + pb) * 4 + group) * 5 + gcbin, size1)
        rate1 = _rates(cell1, ok, var, size1, cell_tp, size_tp)
        out.update(
            {
                "informative": ok,
                "variable": var,
                "cell0": cell0.astype(np.uint16),
                "rate0": rate0,
                "cell_t": cell_t.astype(np.uint16),
                "rate_t": rate_t,
                "cell1": cell1.astype(np.uint16),
                "rate1": rate1,
                "variant_rate": float(var[ok].mean()) if ok.any() else None,
                "informative_bases": int(ok.sum()),
            }
        )
    return out


# ---- per-occurrence sums and the unit tests ----------------------------------------------------

OCC_FIELDS = (
    "bases",
    "m_bases",
    "m_obs",
    "m_exp",
    "m_exp_free",
    "q_sum",
    "h_bases",
    "h_obs",
    "h_e0",
    "h_e1",
    "h_free",
    "g_sum",
    "g_bases",
)


def occurrence_sums(starts: np.ndarray, ends: np.ndarray, arr: dict[str, Any]) -> dict[str, np.ndarray]:
    """For each occurrence, the sums every test needs, gathered over its bases at once."""
    size = arr["length"]
    s = np.clip(starts, 0, size)
    e = np.clip(ends, 0, size)
    keep = e > s
    s, e = s[keep], e[keep]
    lens = e - s
    n = len(s)
    out = {"keep": keep}
    if n == 0:
        for f in OCC_FIELDS:
            out[f] = np.zeros(0)
        return out
    offs = np.concatenate(([0], np.cumsum(lens)[:-1]))
    idx = np.arange(int(lens.sum()), dtype=np.int64) - np.repeat(offs, lens) + np.repeat(s, lens)

    def red(values: np.ndarray) -> np.ndarray:
        return np.add.reduceat(values.astype(np.float64), offs)

    meas = arr["measured"][idx]
    out["bases"] = lens.astype(np.float64)
    out["m_bases"] = red(meas)
    out["m_obs"] = red(arr["constrained"][idx])
    out["m_exp"] = red(arr["rate_m"][arr["cell_m"][idx]] * meas)
    out["m_exp_free"] = red(arr["rate_ms"][arr["cell_ms"][idx]] * meas)
    qv = arr["q"][idx].astype(np.float64)
    out["q_sum"] = red(np.where(meas, (qv - PHYLOP_ANCHOR + 0.5) / PHYLOP_SCALE + PHYLOP_THRESHOLD, 0.0))
    if arr["human"]:
        inf = arr["informative"][idx]
        out["h_bases"] = red(inf)
        out["h_obs"] = red(arr["variable"][idx] & inf)
        out["h_e0"] = red(arr["rate0"][arr["cell0"][idx]] * inf)
        out["h_e1"] = red(arr["rate1"][arr["cell1"][idx]] * inf)
        out["h_free"] = red(arr["rate_t"][arr["cell_t"][idx]] * inf)
    else:
        for f in ("h_bases", "h_obs", "h_e0", "h_e1", "h_free"):
            out[f] = np.zeros(n)
    if arr["gnocchi"] is not None:
        g = arr["gnocchi"][idx]
        ok = ~np.isnan(g)
        out["g_sum"] = red(np.where(ok, g, 0.0))
        out["g_bases"] = red(ok)
    else:
        out["g_sum"] = np.zeros(n)
        out["g_bases"] = np.zeros(n)
    out["mid_ctx"] = arr["ctx"][np.minimum((s + e) // 2, size - 1)]
    return out


def normal_sf(z: float) -> float:
    return 0.5 * math.erfc(z / math.sqrt(2))


def sandwich_z(obs: np.ndarray, exp: np.ndarray) -> float | None:
    """(O - E) over a variance built from the occurrences' own residuals, floored at Poisson."""
    observed, expected = float(obs.sum()), float(exp.sum())
    if expected <= 0:
        return None
    var = max(float(((obs - exp) ** 2).sum()), expected)
    return (observed - expected) / math.sqrt(var)


def row_stats(
    sums: dict[str, np.ndarray], sel: np.ndarray, flanks: dict[str, np.ndarray] | None = None
) -> dict[str, Any]:
    """The tests and the descriptive numbers for one set of occurrences: against the strata, against
    the strata and phyloP, and against the occurrences' own flanks."""
    g = {f: sums[f][sel] for f in OCC_FIELDS}
    occ = int(sel.sum())
    mb = float(g["m_bases"].sum())
    row: dict[str, Any] = {
        "occurrences": occ,
        "bases": int(g["bases"].sum()),
        "phylop_measured_fraction": round(mb / max(1.0, float(g["bases"].sum())), 3),
        "phylop_mean": round(float(g["q_sum"].sum()) / mb, 3) if mb else None,
        "constrained_fraction": round(float(g["m_obs"].sum()) / mb, 4) if mb else None,
        "constrained_expected": round(float(g["m_exp"].sum()) / mb, 4) if mb else None,
        "constrained_expected_stratum_only": round(float(g["m_exp_free"].sum()) / mb, 4) if mb else None,
    }
    zm = sandwich_z(g["m_obs"], g["m_exp"]) if mb >= 200 else None
    row["z_mammals"] = round(zm, 3) if zm is not None else None
    row["p_mammals"] = normal_sf(zm) if zm is not None else None
    hb = float(g["h_bases"].sum())
    e0 = float(g["h_e0"].sum())
    row["human_informative_fraction"] = round(hb / max(1.0, float(g["bases"].sum())), 3)
    row["human_variants"] = int(g["h_obs"].sum())
    row["human_expected_trinucleotide"] = round(e0, 2)
    row["human_expected_given_phylop"] = round(float(g["h_e1"].sum()), 2)
    row["human_oe"] = round(float(g["h_obs"].sum()) / e0, 3) if e0 > 0 else None
    row["human_oe_given_phylop"] = (
        round(float(g["h_obs"].sum()) / float(g["h_e1"].sum()), 3) if g["h_e1"].sum() > 0 else None
    )
    row["human_oe_free"] = (
        round(float(g["h_obs"].sum()) / float(g["h_free"].sum()), 3) if g["h_free"].sum() > 0 else None
    )
    if e0 >= MIN_EXPECTED_VARIANTS:
        z0 = sandwich_z(g["h_obs"], g["h_e0"])
        z1 = sandwich_z(g["h_obs"], g["h_e1"])
        row["z_human"] = round(z0, 3)
        row["p_human_depleted"] = 1.0 - normal_sf(z0)
        row["z_human_given_phylop"] = round(z1, 3) if z1 is not None else None
        row["p_human_given_phylop"] = min(1.0, 2 * normal_sf(abs(z1))) if z1 is not None else None
    else:
        row["z_human"] = row["p_human_depleted"] = row["z_human_given_phylop"] = row[
            "p_human_given_phylop"
        ] = None
    if flanks is not None:
        fl = {f: flanks[f][sel] for f in ("m_bases", "m_obs", "h_obs", "h_e0")}
        fmb = float(fl["m_bases"].sum())
        share = float(fl["m_obs"].sum()) / fmb if fmb else None
        row["constrained_flank_fraction"] = round(share, 4) if share is not None else None
        zf = sandwich_z(g["m_obs"], g["m_bases"] * share) if share and mb >= 200 else None
        row["z_mammals_flank"] = round(zf, 3) if zf is not None else None
        row["p_mammals_flank"] = normal_sf(zf) if zf is not None else None
        fe = float(fl["h_e0"].sum())
        ratio = float(fl["h_obs"].sum()) / fe if fe > 0 else None
        row["human_oe_flank"] = round(ratio, 3) if ratio is not None else None
        if ratio and e0 * ratio >= MIN_EXPECTED_VARIANTS:
            zh = sandwich_z(g["h_obs"], g["h_e0"] * ratio)
            row["z_human_flank"] = round(zh, 3)
            row["p_human_flank_depleted"] = 1.0 - normal_sf(zh)
            row["p_human_flank_two_sided"] = min(1.0, 2 * normal_sf(abs(zh)))
        else:
            row["z_human_flank"] = row["p_human_flank_depleted"] = row["p_human_flank_two_sided"] = None
    gb = float(g["g_bases"].sum())
    row["gnocchi_mean"] = round(float(g["g_sum"].sum()) / gb, 3) if gb else None
    return row


def unit_table(
    name: str, kind: str, positions: list[tuple[int, int]], arr: dict[str, Any]
) -> list[dict[str, Any]]:
    """One row for the unit across its contexts and one per context it occupies at least MIN_OCC times."""
    from genomeos.attribution.lexicon import CONTEXTS

    if not positions:
        return []
    p = np.asarray(positions, dtype=np.int64)
    w = np.minimum(p[:, 1] - p[:, 0], FLANK_MAX)
    # only occurrences whose two flanks lie on the chromosome, so sites and flanks stay aligned
    ok = (p[:, 1] > p[:, 0]) & (p[:, 0] - FLANK_GAP - w >= 0) & (p[:, 1] + FLANK_GAP + w <= arr["length"])
    p, w = p[ok], w[ok]
    if len(p) < MIN_OCC:
        return []
    sums = occurrence_sums(p[:, 0], p[:, 1], arr)
    left = occurrence_sums(p[:, 0] - FLANK_GAP - w, p[:, 0] - FLANK_GAP, arr)
    right = occurrence_sums(p[:, 1] + FLANK_GAP, p[:, 1] + FLANK_GAP + w, arr)
    flanks = {f: left[f] + right[f] for f in ("m_bases", "m_obs", "h_obs", "h_e0")}
    n = len(sums["bases"])
    rows = [{"unit": name, "kind": kind, "context": "*", **row_stats(sums, np.ones(n, dtype=bool), flanks)}]
    counts = np.bincount(sums["mid_ctx"], minlength=len(CONTEXTS))
    for code in np.nonzero(counts >= MIN_OCC)[0]:
        sel = sums["mid_ctx"] == code
        rows.append({"unit": name, "kind": kind, "context": CONTEXTS[code], **row_stats(sums, sel, flanks)})
    return rows


def shifted(positions: list[tuple[int, int]], length: int, rng) -> list[tuple[int, int]]:
    """The unit moved as one piece by a random offset around the chromosome: its occurrences keep
    their spacing and clustering and lose their places."""
    off = rng.randrange(1_000_000, max(1_000_001, length - 1_000_000))
    return [((a + off) % length, (a + off) % length + (b - a)) for a, b in positions]


# ---- the questions ---------------------------------------------------------------------------------

FDR = 0.05


def _bh(rows: list[dict[str, Any]], key: str) -> tuple[float, int, int]:
    from genomeos.attribution.lexicon import benjamini_hochberg

    ps = [r[key] for r in rows if r.get(key) is not None]
    thr, passing = benjamini_hochberg(ps, FDR)
    return thr, passing, len(ps)


FAMILIES = {
    "mammals": "p_mammals",
    "people_depleted": "p_human_depleted",
    "people_given_phylop": "p_human_given_phylop",
    "mammals_vs_flanks": "p_mammals_flank",
    "people_depleted_vs_flanks": "p_human_flank_depleted",
    "people_two_sided_vs_flanks": "p_human_flank_two_sided",
}
CLASSES = (
    "held_mammals",
    "held_people",
    "syntax",
    "people_tighter_than_mammals_predict",
    "value_slot",
    "held_mammals_vs_flanks",
    "held_people_vs_flanks",
    "syntax_vs_flanks",
    "value_slot_vs_flanks",
)


def classify(rows: list[dict[str, Any]], fixed: dict[str, dict] | None = None) -> dict[str, Any]:
    """Syntax, value slot and people-tighter per (unit, context) row, each test against its own
    Benjamini-Hochberg family; with `fixed`, the thresholds of another family are applied instead
    (the decoys and the shifted copies are read with the real units' thresholds)."""
    ctx_rows = [r for r in rows if r["context"] != "*"]
    tests: dict[str, dict] = {}
    for fam, key in FAMILIES.items():
        thr, passing, n = _bh(ctx_rows, key)
        tests[fam] = {"tests": n, "passing": passing, "bh_threshold": thr}
    use = fixed or tests

    def ok(r: dict, fam: str) -> bool:
        p = r.get(FAMILIES[fam])
        return p is not None and use[fam]["passing"] > 0 and p <= use[fam]["bh_threshold"]

    counts: Counter = Counter()
    for r in ctx_rows:
        z1 = r.get("z_human_given_phylop") or 0.0
        r["held_mammals"] = ok(r, "mammals")
        r["held_people"] = ok(r, "people_depleted")
        r["syntax"] = r["held_mammals"] and r["held_people"]
        r["people_tighter_than_mammals_predict"] = ok(r, "people_given_phylop") and z1 < 0
        r["value_slot"] = r["held_mammals"] and ok(r, "people_given_phylop") and z1 > 0
        r["human_measured"] = r.get("p_human_depleted") is not None
        r["held_mammals_vs_flanks"] = ok(r, "mammals_vs_flanks")
        r["held_people_vs_flanks"] = ok(r, "people_depleted_vs_flanks")
        r["syntax_vs_flanks"] = r["held_mammals_vs_flanks"] and r["held_people_vs_flanks"]
        r["value_slot_vs_flanks"] = (
            r["held_mammals_vs_flanks"]
            and ok(r, "people_two_sided_vs_flanks")
            and (r.get("z_human_flank") or 0) > 0
        )
        counts["rows"] += 1
        counts["human_measured"] += r["human_measured"]
        for k in CLASSES:
            counts[k] += r[k]
    measured = [r for r in ctx_rows if r["human_measured"] and r.get("p_mammals") is not None]
    a = sum(1 for r in measured if r["held_mammals"] and r["held_people"])
    b = sum(1 for r in measured if r["held_mammals"] and not r["held_people"])
    c = sum(1 for r in measured if not r["held_mammals"] and r["held_people"])
    d = sum(1 for r in measured if not r["held_mammals"] and not r["held_people"])
    return {
        "tests": tests,
        "thresholds_applied": "own" if fixed is None else "the real units'",
        "counts": dict(counts),
        "two_by_two": {
            "held_both": a,
            "mammals_only": b,
            "people_only": c,
            "neither": d,
            "odds_ratio": round((a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5)), 2),
        },
    }


PAIR_CONTEXTS = (
    "*",
    "intron",
    "enhancer_like",
    "promoter_like",
    "ctcf_only",
    "unknown_regulatory",
    "unknown_fossil",
    "unknown_neutral",
)
PAIR_MIN_EXPECTED = 20.0


def motif_against_decoy(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Each factor against its own decoy in the same context, both read against their own flanks:
    the site-over-flank ratio of the motif minus that of the decoy, on each axis. A positive
    mammal difference and a negative people difference are what a motif, rather than its letters,
    would add. Sign tests over factors; factors with similar matrices share sites, so the p values
    are optimistic."""
    idx = {(r["kind"], r["unit"].split(":", 1)[-1], r["context"]): r for r in rows}
    out: dict[str, Any] = {}
    for ctx in PAIR_CONTEXTS:
        people: list[float] = []
        mammals: list[float] = []
        for (kind, factor, c), x in idx.items():
            if kind != "jaspar95" or c != ctx:
                continue
            y = idx.get(("decoy95", factor, c))
            if y is None:
                continue
            if (
                all(v.get("human_oe") and v.get("human_oe_flank") for v in (x, y))
                and min(x["human_expected_trinucleotide"], y["human_expected_trinucleotide"])
                >= PAIR_MIN_EXPECTED
            ):
                people.append(x["human_oe"] / x["human_oe_flank"] - y["human_oe"] / y["human_oe_flank"])
            if x.get("constrained_flank_fraction") and y.get("constrained_flank_fraction"):
                mammals.append(
                    x["constrained_fraction"] / x["constrained_flank_fraction"]
                    - y["constrained_fraction"] / y["constrained_flank_fraction"]
                )
        lower = sum(1 for d in people if d < 0)
        higher = sum(1 for d in mammals if d > 0)
        out[ctx] = {
            "people_pairs": len(people),
            "people_motif_less_variable": lower,
            "people_sign_p": _sign_p(lower, len(people)),
            "people_median_difference": round(float(np.median(people)), 3) if people else None,
            "mammal_pairs": len(mammals),
            "mammal_motif_more_constrained": higher,
            "mammal_sign_p": _sign_p(higher, len(mammals)),
            "mammal_median_difference": round(float(np.median(mammals)), 3) if mammals else None,
        }
    return out


def conditional_people(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Among rows held across mammals against their flanks, the share also depleted among people
    against their flanks, for motifs and for decoys: the people axis adds something specific to the
    motif only if that share is higher for motifs."""
    out = {}
    for kind in ("jaspar95", "decoy95"):
        held = [
            r for r in rows if r["kind"] == kind and r["context"] != "*" and r.get("held_mammals_vs_flanks")
        ]
        both = sum(1 for r in held if r.get("held_people_vs_flanks"))
        out[kind] = {
            "held_mammals_vs_flanks": len(held),
            "also_people": both,
            "share": round(both / len(held), 4) if held else None,
        }
    return out


def rank_correlation(x: list[float], y: list[float]) -> float | None:
    if len(x) < 10:
        return None
    rx = np.argsort(np.argsort(np.asarray(x)))
    ry = np.argsort(np.argsort(np.asarray(y)))
    return round(float(np.corrcoef(rx, ry)[0, 1]), 3)


def own_element(
    factor_sites: dict[str, list[tuple[int, int]]], ccres: list[tuple[int, int]], arr: dict[str, Any]
) -> dict[str, Any]:
    """Each factor's sites inside a cCRE against the rest of the same cCRE: on phyloP, the sites'
    constrained bases against the rest's constrained share; among people, the sites' recurring
    substitutions against their trinucleotide expectation scaled by the rest's own observed over
    expected. A site that only borrows its element's conservation or its element's variation does
    not pass."""
    from genomeos.attribution.lexicon import benjamini_hochberg

    hosts = np.asarray(sorted(ccres), dtype=np.int64)
    host_sums = occurrence_sums(hosts[:, 0], hosts[:, 1], arr)
    hs, he = hosts[host_sums["keep"], 0], hosts[host_sums["keep"], 1]
    rows = []
    pool = Counter()
    for name, sites in factor_sites.items():
        if len(sites) < MIN_OCC:
            continue
        p = np.asarray(sites, dtype=np.int64)
        i = np.searchsorted(hs, p[:, 0], side="right") - 1
        inside = (i >= 0) & (p[:, 1] <= he[np.maximum(i, 0)])
        if inside.sum() < MIN_OCC:
            continue
        p, i = p[inside], i[inside]
        s = occurrence_sums(p[:, 0], p[:, 1], arr)
        i = i[s["keep"]]
        rest_mb = host_sums["m_bases"][i] - s["m_bases"]
        rest_mo = host_sums["m_obs"][i] - s["m_obs"]
        share = np.divide(rest_mo, rest_mb, out=np.zeros_like(rest_mo), where=rest_mb > 0)
        m_exp = s["m_bases"] * share
        zm = sandwich_z(s["m_obs"], m_exp)
        row = {
            "factor": name,
            "sites_in_ccres": int(len(i)),
            "sites_total": len(sites),
            "constrained_fraction": round(float(s["m_obs"].sum()) / max(1.0, float(s["m_bases"].sum())), 4),
            "rest_of_element_fraction": round(float(m_exp.sum()) / max(1.0, float(s["m_bases"].sum())), 4),
            "z_mammals": round(zm, 3) if zm is not None else None,
            "p_mammals": normal_sf(zm) if zm is not None else None,
        }
        pool["m_obs"] += float(s["m_obs"].sum())
        pool["m_exp"] += float(m_exp.sum())
        pool["m_bases"] += float(s["m_bases"].sum())
        if arr["human"]:
            rest_o = host_sums["h_obs"][i] - s["h_obs"]
            rest_e = host_sums["h_e0"][i] - s["h_e0"]
            ratio = np.divide(rest_o, rest_e, out=np.ones_like(rest_o), where=rest_e > 0.5)
            h_exp = s["h_e0"] * ratio
            if h_exp.sum() >= MIN_EXPECTED_VARIANTS:
                zh = sandwich_z(s["h_obs"], h_exp)
                row.update(
                    {
                        "human_variants": int(s["h_obs"].sum()),
                        "human_expected_from_rest": round(float(h_exp.sum()), 2),
                        "z_human": round(zh, 3),
                        "p_human_depleted": 1.0 - normal_sf(zh),
                    }
                )
                pool["h_obs"] += float(s["h_obs"].sum())
                pool["h_exp"] += float(h_exp.sum())
        rows.append(row)
    thr_m, pass_m = benjamini_hochberg([r["p_mammals"] for r in rows if r["p_mammals"] is not None], FDR)
    hp = [r["p_human_depleted"] for r in rows if r.get("p_human_depleted") is not None]
    thr_h, pass_h = benjamini_hochberg(hp, FDR)
    both = [
        r
        for r in rows
        if r["p_mammals"] is not None
        and r["p_mammals"] <= thr_m
        and pass_m
        and r.get("p_human_depleted") is not None
        and r["p_human_depleted"] <= thr_h
        and pass_h
    ]
    rows.sort(key=lambda r: -(r["z_mammals"] or 0))
    return {
        "factors_tested": len(rows),
        "mammals": {
            "tests": sum(1 for r in rows if r["p_mammals"] is not None),
            "passing": pass_m,
            "bh_threshold": thr_m,
        },
        "people": {"tests": len(hp), "passing": pass_h, "bh_threshold": thr_h},
        "passing_both": len(both),
        "pooled": {
            "constrained_fraction": round(pool["m_obs"] / max(1.0, pool["m_bases"]), 4),
            "rest_of_element_fraction": round(pool["m_exp"] / max(1.0, pool["m_bases"]), 4),
            "human_observed_over_expected_from_rest": (
                round(pool["h_obs"] / pool["h_exp"], 3) if pool["h_exp"] else None
            ),
            "note": "sites of different factors overlap, so pooled numbers count shared bases more than once",
        },
        "both_axes": [
            {
                k: r.get(k)
                for k in (
                    "factor",
                    "sites_in_ccres",
                    "constrained_fraction",
                    "rest_of_element_fraction",
                    "z_mammals",
                    "z_human",
                )
            }
            for r in both[:25]
        ],
        "top_mammals": [
            {
                k: r.get(k)
                for k in (
                    "factor",
                    "sites_in_ccres",
                    "constrained_fraction",
                    "rest_of_element_fraction",
                    "z_mammals",
                    "z_human",
                )
            }
            for r in rows[:15]
        ],
    }


def fossil_tier(rows: list[dict[str, Any]], tier_rows: list[dict[str, Any]]) -> dict[str, Any]:
    """The fossil tier read on both axes: the same repeat subfamily inside the tier against its copies
    in introns, with the context-free expectations (the tier and the intron differ by definition,
    so a context-matched expectation would absorb the very difference asked about). phyloP selected
    the tier (low block constraint), so the mammal comparison is circular by construction; the
    people axis is not."""
    by_unit: dict[str, dict[str, dict]] = defaultdict(dict)
    for r in rows:
        if r["kind"] == "repeat":
            by_unit[r["unit"]][r["context"]] = r
    pairs = []
    for unit, ctxs in by_unit.items():
        f, i = ctxs.get("unknown_fossil"), ctxs.get("intron")
        if not f or not i:
            continue
        if f["human_oe_free"] is None or i["human_oe_free"] is None:
            continue
        if (
            f["human_expected_trinucleotide"] < MIN_EXPECTED_VARIANTS
            or i["human_expected_trinucleotide"] < MIN_EXPECTED_VARIANTS
        ):
            continue
        pairs.append(
            {
                "subfamily": unit.split(":", 1)[1],
                "fossil_occurrences": f["occurrences"],
                "intron_occurrences": i["occurrences"],
                "fossil_human_oe": f["human_oe_free"],
                "intron_human_oe": i["human_oe_free"],
                "fossil_constrained": f["constrained_fraction"],
                "intron_constrained": i["constrained_fraction"],
                "fossil_phylop_mean": f["phylop_mean"],
                "intron_phylop_mean": i["phylop_mean"],
            }
        )
    lower_h = sum(1 for p in pairs if p["fossil_human_oe"] < p["intron_human_oe"])
    lower_m = sum(1 for p in pairs if (p["fossil_constrained"] or 0) < (p["intron_constrained"] or 0))

    sign_p = _sign_p

    hs = sorted(p["fossil_human_oe"] - p["intron_human_oe"] for p in pairs)
    return {
        "subfamilies_paired": len(pairs),
        "people": {
            "fossil_copies_less_variable": lower_h,
            "sign_test_p": sign_p(lower_h, len(pairs)),
            "median_difference_oe": round(hs[len(hs) // 2], 3) if hs else None,
        },
        "mammals": {
            "fossil_copies_less_constrained": lower_m,
            "sign_test_p": sign_p(lower_m, len(pairs)),
            "note": "circular: the tier is defined by low block-level phyloP",
        },
        "tiers": tier_rows,
        "pairs": sorted(pairs, key=lambda p: -p["fossil_occurrences"])[:20],
    }


def replication_strata(chrom: str, strata: int = 3) -> dict[int, int] | None:
    """Kilobase -> replication-timing stratum (terciles of the Repli-seq signal), from the panel
    store's own timing table where genomeos-h1 has measured it."""
    from genomeos.attribution.human_panel import CACHE

    path = CACHE / chrom / "replication_timing.json"
    if not path.exists():
        return None
    rt = {int(k): float(v) for k, v in json.loads(path.read_text()).items()}
    values = sorted(rt.values())
    edges = [values[int(len(values) * k / strata)] for k in range(1, strata)]
    return {kb: sum(1 for e in edges if v >= e) for kb, v in rt.items()}


REFERENCE_CONTEXTS = ("intron", "other", "unknown_neutral", "unknown_regulatory")


def fossil_matched(
    units: dict[str, list[tuple[int, int]]], arr: dict[str, Any], rt: dict[int, int] | None
) -> dict[str, Any]:
    """The same repeat subfamily inside the fossil tier against its copies in each reference context,
    on the people axis, matched within replication-timing terciles (chr21) and with the
    trinucleotide x GC expectation of each base. For every subfamily the fossil copies' variants are
    compared with what the reference copies' observed-over-expected, tercile by tercile, predicts."""
    from genomeos.attribution.lexicon import CONTEXT_CODE

    fossil_code = CONTEXT_CODE["unknown_fossil"]
    out: dict[str, Any] = {"replication_timing_matched": rt is not None}
    for ref in REFERENCE_CONTEXTS:
        ref_code = CONTEXT_CODE[ref]
        per = []
        pool_obs = pool_exp = 0.0
        resid: list[float] = []
        for name, positions in units.items():
            if not name.startswith("repeat:") or len(positions) < 2 * MIN_OCC:
                continue
            p = np.asarray(positions, dtype=np.int64)
            sm = occurrence_sums(p[:, 0], p[:, 1], arr)
            kept = p[sm["keep"]]
            mids = (kept[:, 0] + kept[:, 1]) // 2
            if rt is not None:
                tercile = np.array([rt.get(int(m) // 1000, -1) for m in mids])
            else:
                tercile = np.zeros(len(mids), dtype=int)
            f = (sm["mid_ctx"] == fossil_code) & (tercile >= 0)
            r = (sm["mid_ctx"] == ref_code) & (tercile >= 0)
            if f.sum() < MIN_OCC or r.sum() < MIN_OCC:
                continue
            exp = np.zeros(len(mids))
            usable = np.zeros(len(mids), dtype=bool)
            for t in np.unique(tercile[f]):
                rt_sel = r & (tercile == t)
                re = float(sm["h_free"][rt_sel].sum())
                if re < MIN_EXPECTED_VARIANTS:
                    continue
                ratio = float(sm["h_obs"][rt_sel].sum()) / re
                sel = f & (tercile == t)
                exp[sel] = sm["h_free"][sel] * ratio
                usable |= sel
            e = float(exp[usable].sum())
            if e < MIN_EXPECTED_VARIANTS:
                continue
            o = float(sm["h_obs"][usable].sum())
            per.append(
                {
                    "subfamily": name.split(":", 1)[1],
                    "fossil_copies": int(usable.sum()),
                    "observed": int(o),
                    "expected": round(e, 1),
                    "ratio": round(o / e, 3),
                }
            )
            pool_obs += o
            pool_exp += e
            resid.extend((sm["h_obs"][usable] - exp[usable]).tolist())
        more = sum(1 for x in per if x["ratio"] > 1)
        var = max(sum(x * x for x in resid), pool_exp) if resid else 0.0
        z = (pool_obs - pool_exp) / math.sqrt(var) if var else None
        out[ref] = {
            "subfamilies": len(per),
            "fossil_copies_more_variable": more,
            "sign_test_p": _sign_p(more, len(per)),
            "pooled_ratio": round(pool_obs / pool_exp, 3) if pool_exp else None,
            "z": round(z, 2) if z is not None else None,
            "largest": sorted(per, key=lambda x: -x["expected"])[:10],
        }
    return out


def _sign_p(k: int, n: int) -> float | None:
    if not n:
        return None
    lo = _binom_cdf(k, n)
    hi = 1.0 - _binom_cdf(k - 1, n)
    return float(min(1.0, 2 * min(lo, hi)))


def _binom_cdf(k: int, n: int) -> float:
    if k < 0:
        return 0.0
    return sum(math.comb(n, j) for j in range(0, min(k, n) + 1)) / 2**n


def run(chrom: str, results_dir=None, progress=None) -> dict[str, Any]:
    """The re-ask on one chromosome."""
    import random

    from genomeos.attribution import lexicon as L  # noqa: N812 - the first lexicon, read as a namespace
    from genomeos.coords import Locus
    from genomeos.genome.index import IndexedGenome
    from genomeos.results import RESULTS_DIR, load_result

    results_dir = results_dir or RESULTS_DIR
    t0 = time.time()

    def say(m: str) -> None:
        if progress:
            progress(f"{time.time() - t0:6.0f} s  {m}")

    genome = IndexedGenome(L.REFERENCE / f"{chrom}.fa")
    length = genome.lengths[chrom]
    seq = str(genome.fetch(Locus(chrom, 0, length))).upper()
    genome.close()
    maps = L.Maps(chrom, length, results_dir)
    windows = L.Windows(seq, maps)
    q = phylop_track(chrom, length, progress=say)
    say(
        f"phyloP: {(q != PHYLOP_MISSING).mean():.3f} of bases measured, "
        f"{(q >= PHYLOP_ANCHOR).mean():.4f} constrained"
    )
    gnocchi = gnocchi_track(chrom, length, progress=say)
    human = panel_tracks(chrom, length)
    say(
        "panel store read"
        if human is not None
        else "no panel store for this chromosome: the people axis is absent"
    )
    arr = base_arrays(seq, maps, windows, q, human, gnocchi)
    say("per-base arrays and expectation tables built")
    from genomeos.genome.motifs import load_motifs

    codes = encode(seq)
    index = kmer_order(codes)
    del codes
    profiles = load_motifs(relative=MOTIF_RELATIVE)
    sites, aliases = dedupe_sites(scan_chromosome(seq, motifs=profiles, index=index))
    decoys, _ = dedupe_sites(scan_chromosome(seq, motifs=decoy_motifs(profiles), index=index))
    del index
    say(
        f"JASPAR at {MOTIF_RELATIVE}: {len(sites)} distinct site sets, "
        f"{sum(len(v) for v in sites.values())} merged sites"
    )
    units = {
        k: v for k, v in L.curated_units(chrom, maps, results_dir).items() if not k.startswith("jaspar:")
    }
    for name, v in sites.items():
        units[f"jaspar95:{name}"] = v
    for name, v in decoys.items():
        units[f"decoy95:{name}"] = v
    found = L.seeds(seq)
    for w, pos in found.items():
        units[f"seed:{w}"] = [(p, p + L.SEED_K) for p in pos[:-1]]
    say(f"{len(units)} units")
    rows: list[dict[str, Any]] = []
    for name, positions in units.items():
        rows.extend(unit_table(name, name.split(":")[0], positions, arr))
    say(f"{len(rows)} unit rows")
    rng = random.Random(L.SEED)
    null_rows: list[dict[str, Any]] = []
    for name, positions in units.items():
        null_rows.extend(unit_table(name, name.split(":")[0], shifted(positions, length, rng), arr))
    say(f"{len(null_rows)} shifted rows")
    real_rows = [r for r in rows if r["kind"] != "decoy95"]
    decoy_rows = [r for r in rows if r["kind"] == "decoy95"]
    real = classify(real_rows)
    null = classify([r for r in null_rows if r["kind"] != "decoy95"], fixed=real["tests"])
    decoy = classify(decoy_rows, fixed=real["tests"])
    paired = motif_against_decoy(rows)
    conditional = conditional_people(rows)
    motif_vs_decoy = {}
    for cls in CLASSES:
        jr = [r for r in real_rows if r["kind"] == "jaspar95" and r["context"] != "*"]
        motif_vs_decoy[cls] = {
            "motif_rows": sum(1 for r in jr if r[cls]),
            "motif_tested": len(jr),
            "decoy_rows": sum(1 for r in decoy_rows if r["context"] != "*" and r[cls]),
            "decoy_tested": sum(1 for r in decoy_rows if r["context"] != "*"),
        }
    # the two people measures against each other, per context row
    both = [
        r
        for r in real_rows
        if r["context"] != "*" and r["z_human"] is not None and r["gnocchi_mean"] is not None
    ]
    agreement = {
        "rows": len(both),
        "rank_correlation_gnocchi_vs_panel_depletion": rank_correlation(
            [r["gnocchi_mean"] for r in both], [-r["z_human"] for r in both]
        ),
        "rank_correlation_phylop_vs_panel_depletion": rank_correlation(
            [r["z_mammals"] or 0.0 for r in both], [-r["z_human"] for r in both]
        ),
        "rank_correlation_phylop_vs_gnocchi": rank_correlation(
            [r["z_mammals"] or 0.0 for r in both], [r["gnocchi_mean"] for r in both]
        ),
    }
    ccres = [(int(f[1]), int(f[2])) for f in L._bed(results_dir / f"ccres_{chrom}.bed.gz")]
    element = own_element(sites, ccres, arr)
    element_decoy = own_element(decoys, ccres, arr)
    say("motif sites against their own elements")
    tier_rows = []
    for tier in (
        "unknown_fossil",
        "unknown_neutral",
        "unknown_regulatory",
        "unknown_constrained",
        "intron",
        "cds",
    ):
        code = L.CONTEXT_CODE[tier]
        mask = arr["ctx"] == code
        idx = np.nonzero(mask)[0]
        if not len(idx):
            continue
        # contiguous runs of the context as its occurrences
        cut = np.nonzero(np.diff(idx) > 1)[0]
        starts = np.concatenate(([idx[0]], idx[cut + 1]))
        ends = np.concatenate((idx[cut] + 1, [idx[-1] + 1]))
        s = occurrence_sums(starts, ends, arr)
        mb = float(s["m_bases"].sum())
        tier_rows.append(
            {
                "context": tier,
                "bases": int(s["bases"].sum()),
                "constrained_fraction": round(float(s["m_obs"].sum()) / mb, 4) if mb else None,
                "phylop_mean": round(float(s["q_sum"].sum()) / mb, 3) if mb else None,
                "human_oe_free": (
                    round(float(s["h_obs"].sum()) / float(s["h_free"].sum()), 3)
                    if s["h_free"].sum()
                    else None
                ),
                "human_informative_fraction": round(
                    float(s["h_bases"].sum()) / max(1.0, float(s["bases"].sum())), 3
                ),
            }
        )
    fossil = fossil_tier(rows, tier_rows)
    fossil["matched_on_people"] = fossil_matched(units, arr, replication_strata(chrom))
    say("fossil tier matched")
    candidates = [
        (c["start"], c["end"])
        for c in (load_result("syntax_candidates_genome_wide", results_dir) or {}).get("candidates", [])
        if c["chrom"] == chrom
    ]
    cand_rows = []
    for a, b in candidates:
        s = occurrence_sums(np.array([a]), np.array([b]), arr)
        cand_rows.append(
            {
                "start": a,
                "end": b,
                "constrained_fraction": round(float(s["m_obs"][0] / s["m_bases"][0]), 4)
                if s["m_bases"][0]
                else None,
                "human_variants": int(s["h_obs"][0]),
                "human_expected": round(float(s["h_e0"][0]), 2),
                "human_expected_given_phylop": round(float(s["h_e1"][0]), 2),
            }
        )
    ctx_rows = [r for r in real_rows if r["context"] != "*"]
    keys = (
        "unit",
        "context",
        "occurrences",
        "constrained_fraction",
        "constrained_expected",
        "z_mammals",
        "human_variants",
        "human_expected_trinucleotide",
        "human_oe",
        "human_oe_given_phylop",
        "z_human",
        "z_human_given_phylop",
        "gnocchi_mean",
    )

    def pick(flag: str) -> list[dict[str, Any]]:
        def order(r: dict[str, Any]) -> float:
            return (r["z_human_given_phylop"] or 0) if flag == "value_slot" else -(r["z_mammals"] or 0)

        return sorted((r for r in ctx_rows if r.get(flag)), key=order)

    out = {
        "chrom": chrom,
        "length": length,
        "axes": {
            "mammals": (
                "Zoonomia phyloP 241 mammals per base, constrained at >= 2.27 (byte track, bin edge at 2.27)"
            ),
            "people": (
                "HPRC 89-assembly panel per base (genomeos-h1's store): recurring substitution columns "
                "among informative bases, expected from trinucleotide x context x GC bin, "
                "and again with the phyloP bin"
                if human is not None
                else "absent: no panel store for this chromosome"
            ),
            "people_for_comparison": "gnomAD Gnocchi Z per kilobase, joined per base, not used as the axis",
            "motifs": (
                f"JASPAR 2026 CORE vertebrates at {MOTIF_RELATIVE} of the matrix range, whole chromosome, "
                "both strands; column-permuted decoys scanned the same way"
            ),
        },
        "phylop": {
            "measured_fraction": round(float(arr["measured"].mean()), 4),
            "constrained_fraction": round(float(arr["constrained"].sum() / max(1, arr["measured"].sum())), 4),
        },
        "panel": {"informative_bases": arr.get("informative_bases"), "variant_rate": arr.get("variant_rate")},
        "units": len(units),
        "unit_rows": len(rows),
        "context_rows": len(ctx_rows),
        "real": real,
        "shifted_calibration": null,
        "decoy_motifs": decoy,
        "motif_against_decoy": motif_vs_decoy,
        "motif_against_own_decoy_paired": paired,
        "people_given_mammals_motif_and_decoy": conditional,
        "motif_aliases": {k: v for k, v in aliases.items() if v},
        "measures_agreement": agreement,
        "motif_sites_against_own_element": element,
        "decoy_sites_against_own_element": {
            k: v for k, v in element_decoy.items() if k not in ("both_axes", "top_mammals")
        },
        "fossil_tier": fossil,
        "syntax_candidates": cand_rows,
        "syntax_rows": [{k: r.get(k) for k in keys} for r in pick("syntax")[:30]],
        "value_slot_rows": [{k: r.get(k) for k in keys} for r in pick("value_slot")[:30]],
        "people_tighter_rows": [
            {k: r.get(k) for k in keys} for r in pick("people_tighter_than_mammals_predict")[:20]
        ],
        "syntax_vs_flanks_rows": [
            {
                k: r.get(k)
                for k in (
                    *keys,
                    "constrained_flank_fraction",
                    "z_mammals_flank",
                    "human_oe_flank",
                    "z_human_flank",
                )
            }
            for r in sorted(
                (r for r in ctx_rows if r.get("syntax_vs_flanks")), key=lambda r: -(r["z_mammals_flank"] or 0)
            )[:30]
        ],
        "value_slot_vs_flanks_rows": [
            {
                k: r.get(k)
                for k in (
                    *keys,
                    "constrained_flank_fraction",
                    "z_mammals_flank",
                    "human_oe_flank",
                    "z_human_flank",
                )
            }
            for r in sorted(
                (r for r in ctx_rows if r.get("value_slot_vs_flanks")),
                key=lambda r: -(r["z_human_flank"] or 0),
            )[:30]
        ],
        "syntax_vs_flanks_by_kind": dict(Counter(r["kind"] for r in ctx_rows if r.get("syntax_vs_flanks"))),
        "syntax_vs_flanks_by_context": dict(
            Counter(r["context"] for r in ctx_rows if r.get("syntax_vs_flanks"))
        ),
        "syntax_by_kind": dict(Counter(r["kind"] for r in ctx_rows if r.get("syntax"))),
        "syntax_by_context": dict(Counter(r["context"] for r in ctx_rows if r.get("syntax"))),
        "value_slots_by_kind": dict(Counter(r["kind"] for r in ctx_rows if r.get("value_slot"))),
        "seconds": round(time.time() - t0, 1),
    }
    directory = KNOWLEDGE
    directory.mkdir(parents=True, exist_ok=True)
    with gzip.open(directory / f"lexicon_axes_rows_{chrom}.json.gz", "wt") as fh:
        json.dump(rows, fh)
    return out


def run_and_save(chrom: str, progress=None) -> dict[str, Any]:
    from genomeos.results import save_result

    out = run(chrom, progress=progress)
    save_result(f"lexicon_axes_{chrom}", out)
    return out


def main() -> None:  # pragma: no cover - the job entry point
    import sys

    for chrom in sys.argv[1:] or ["chr21", "chr22"]:
        out = run_and_save(chrom, progress=lambda m: print("  " + m, flush=True))
        print(
            json.dumps(
                {
                    k: out[k]
                    for k in ("chrom", "real", "shifted_calibration", "measures_agreement", "seconds")
                },
                indent=1,
            )
        )
