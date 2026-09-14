# SPDX-License-Identifier: AGPL-3.0-or-later
"""Compression as understanding: bits per base of a chromosome under generic models and under the
biology the project knows.

Albert's probe of 2026-09-14: find which way of classifying the genome buys the most understanding
per unit of effort. This is the information-theoretic one. The best model of a sequence is the one
that describes it in the fewest bits (minimum description length), so every piece of knowledge the
project holds is turned into a predictive model of the next base, the chromosome is coded under
it, and the gain is measured in bits against generic baselines. A layer that does not shorten the
code is not knowledge about the sequence, however it reads in a table. It needs no labels to score
against, which is why it can speak about the 98%.

The code lengths are ideal (-log2 of the probability the model gave the base that came). An
arithmetic coder reaches them within two bits per chromosome, so nothing here is left to a coder.

How the models are kept honest:

- **Held-out fitting.** Every static model, generic or conditioned on an annotation, is fitted on
  a different chromosome (chr22 for chr21 and chr21 for chr22 by default), both strands, and applied
  to the chromosome being coded. Its parameters are not charged to that chromosome; they are
  knowledge learned elsewhere, the way a repeat library is.
- **Adaptive models pay for themselves.** The models fitted on the chromosome itself learn as they
  code (Krichevsky-Trofimov counts updated after each base, the reverse strand's contexts added as
  soon as their bases are known), so the decoder can rebuild them and no parameter is free.
- **Annotations are paid for.** A layer's annotation (where the repeats are and which subfamily,
  where the coding exons start and their phase, the CpG islands, the cCREs, the motif sites, the
  blocks and their tiers) is side information the decoder needs, so its description length
  (Elias delta for positions and lengths, adaptive categorical codes for labels, the label names
  as text) is subtracted from the layer's gain. An annotation found by reading the sequence
  (RepeatMasker, a CpG island, a motif hit) is not cheating for that reason: its cost is paid.
- **Mixing is causal.** Models are combined by a windowed Bayesian mixture: each model's weight at
  a base is its likelihood over the previous W bases, so the mixture uses only bases already
  coded. W and the mixing temperature are chosen from a grid of twelve on the generic stack and
  the log2(12) bits of that choice are charged. A layer's models are inactive outside the region
  its annotation claims (the decoder knows it), so a layer is never charged for diluting the
  mixture where it has nothing to say.
- **A control prices the duplication alignment.** The generic adaptive models are also run primed
  with the partners' sequence, so the alignment is measured against a coder that holds the same
  bases.
- **CpG islands are modelled before anything reads them as structure.** genomeos-i1's lexicon
  found CG-containing words everywhere because an order-2 chain fitted over a context cannot
  absorb the clustering of CpG islands. Here a causal GC and CpG observed/expected state over the
  previous kilobase (and the island annotation) is a layer of its own, and the regulatory, motif
  and tier layers are also read on top of everything else (leave-one-out), so what they add is not
  GC composition by another name.

The layers: RepeatMasker subfamilies with divergence, segmental duplications as copies of a locus
the decoder already holds, GC and CpG structure, GENCODE coding sequence with its codon phase,
ENCODE cCREs, the JASPAR sites the motif runs recorded, and the UNKNOWN blocks with their budget
tiers. Each is reported alone over the naive and the generic stacks, over the repeat-aware stack,
cumulatively and left out of the full stack, both per base of the region it claims and per
megabase of the chromosome, net of its annotation. Bits per base are then read per budget tier,
on all bases and on unique bases (no repeat, no duplication) standardised to the neutral tier's
GC, with a block bootstrap for constrained_unknown against neutral.

Nothing is fetched but the UCSC CpG islands (one request per chromosome, cached); the sequence,
RepeatMasker, duplications, GENCODE, cCREs, motif runs and budgets are already on disk.
"""

from __future__ import annotations

import bz2
import gzip
import json
import lzma
import math
import resource
import sys
import time
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from genomeos.results import RESULTS_DIR, load_result, save_result

REFERENCE = Path("data/reference")
CACHE = Path("data/knowledge/compress")
UCSC_API = "https://api.genome.ucsc.edu/getData/track"
KARYOTYPE = tuple([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])
TIERS = ("structural", "fossil", "regulatory", "constrained_unknown", "neutral")

SCALE = 4096  # model code lengths are kept as uint16 in 1/4096 of a bit
MAX_CODE = 65534  # a model probability below 2^-16 is floored there; the mixture floor sits far above
INACTIVE = 65535  # a model's code where it does not apply: it takes no weight in the mixture
NAIVE_ORDERS = (1, 2, 3, 4, 6, 8, 10, 12)  # static orders fitted on the other chromosome
ADAPTIVE_ORDERS = (12, 16, 20, 24)  # orders fitted on the chromosome itself as it is coded: blind copies
MARKOV_ORDERS = tuple(range(17))  # the baseline table
HIGH_ORDER = 8  # above this order the estimator's pseudocount drops, since long contexts are deterministic
ALPHA_LOW, ALPHA_HIGH = 0.5, 1 / 16
MIX_WINDOWS = (8, 16, 32, 64)
MIX_BETAS = (0.25, 0.5, 1.0)
MIX_FLOOR = 1e-3  # share of the uniform distribution in every mixture: no base costs more than 12 bits
MIX_CHUNK = 1 << 19
DIRECT_BITS = 26  # contexts, state and base that fit this many bits are counted in a direct table
PART_ELEMENTS = 1 << 25  # a count is taken in parts of the hash space, each at most this many rows
GC_WINDOW = 1_000  # the causal composition window
GC_EDGES = (0.35, 0.40, 0.45, 0.50, 0.55)
OE_EDGES = (0.25, 0.60)  # CpG observed/expected: depleted, intermediate, island-like
DIV_EDGES = (0.05, 0.10, 0.20, 0.30)  # RepeatMasker divergence bins
COPY_K = 12  # anchor width of the duplication copy model
TANDEM_CLASSES = {"Simple_repeat", "Low_complexity", "Satellite"}
BOOTSTRAP = 1_000
SEED = 20260914
GOLD = np.uint64(0x9E3779B97F4A7C15)
LN2 = math.log(2.0)

EVIDENCE = {
    "sequence": "curated: GRCh38 primary assembly, A/C/G/T only (other letters are assembly gaps)",
    "repeats": "curated: UCSC RepeatMasker (rmsk_<chrom>)",
    "duplications": "curated: UCSC genomicSuperDups (superdups_<chrom>)",
    "cpg_islands": "curated: UCSC cpgIslandExt, fetched once per chromosome through the REST API",
    "coding": "curated: GENCODE v50 canonical CDS of protein-coding genes, with phase",
    "registry": "curated: ENCODE SCREEN cCREs v3",
    "motifs": "predicted: JASPAR 2026 sites recorded by the motif runs (motifs_<chrom>)",
    "tiers": "inferred: budget tiers (budget_<chrom>) and UNKNOWN block classes (unknown_<chrom>)",
    "code_length": "measured: ideal code length, -log2 of the mixture's probability for each base",
}

_LUT = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _LUT[_b] = _i
    _LUT[_b + 32] = _i


class Ledger:
    """What a run cost: seconds per stage, bytes read from disk, requests and bytes over the network."""

    def __init__(self, progress=None) -> None:
        self.t0 = time.perf_counter()
        self.last = self.t0
        self.stages: dict[str, float] = {}
        self.disk_bytes = 0
        self.requests = 0
        self.network_bytes = 0
        self.progress = progress

    def read(self, path: Path) -> bytes:
        data = Path(path).read_bytes()
        self.disk_bytes += len(data)
        return data

    def lap(self, stage: str) -> None:
        now = time.perf_counter()
        self.stages[stage] = round(self.stages.get(stage, 0.0) + now - self.last, 2)
        self.last = now
        if self.progress:
            self.progress(f"{now - self.t0:7.1f} s  {stage}")

    def summary(self) -> dict[str, Any]:
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        rss_bytes = rss if sys.platform == "darwin" else rss * 1024
        return {
            "seconds": round(time.perf_counter() - self.t0, 1),
            "stages_seconds": self.stages,
            "peak_memory_gb": round(rss_bytes / 1e9, 2),
            "disk_bytes_read": self.disk_bytes,
            "requests": self.requests,
            "network_bytes": self.network_bytes,
        }


# ---- the sequence ---------------------------------------------------------------------------


def parse_fasta_codes(raw: bytes) -> np.ndarray:
    """One FASTA record as base codes over genome coordinates: A0 C1 G2 T3, anything else 4."""
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    if raw.startswith(b">"):
        raw = raw[raw.index(b"\n") + 1 :]
    raw = raw.replace(b"\n", b"").replace(b"\r", b"")
    return _LUT[np.frombuffer(raw, dtype=np.uint8)]


@dataclass
class Chromosome:
    """A chromosome as the stream of its A/C/G/T bases, with the genome position of each."""

    name: str
    genome: np.ndarray  # uint8 codes over genome coordinates
    s: np.ndarray = field(init=False)  # uint8 codes of the A/C/G/T bases only
    gpos: np.ndarray = field(init=False)  # int64 genome position of each stream base

    parts: list[Chromosome] = field(init=False)  # the chromosomes a training set is made of

    def __post_init__(self) -> None:
        keep = self.genome < 4
        self.s = self.genome[keep]
        self.gpos = np.flatnonzero(keep)
        self.parts = [self]

    @classmethod
    def training_set(cls, chroms: list[Chromosome]) -> Chromosome:
        """Several chromosomes as one training stream; annotations are painted part by part."""
        if len(chroms) == 1:
            return chroms[0]
        out = cls.__new__(cls)
        out.name = ",".join(c.name for c in chroms)
        out.genome = np.zeros(0, dtype=np.uint8)
        out.s = np.concatenate([c.s for c in chroms])
        out.gpos = np.zeros(0, dtype=np.int64)
        out.parts = list(chroms)
        return out

    @classmethod
    def load(cls, name: str, reference: Path = REFERENCE, ledger: Ledger | None = None) -> Chromosome:
        p = reference / f"{name}.fa"
        if not p.exists():
            p = reference / f"{name}.fa.gz"
        raw = ledger.read(p) if ledger else p.read_bytes()
        return cls(name, parse_fasta_codes(raw))

    @classmethod
    def from_text(cls, name: str, text: str) -> Chromosome:
        return cls(name, _LUT[np.frombuffer(text.encode(), dtype=np.uint8)])

    @property
    def n(self) -> int:
        return len(self.s)

    def stream_range(self, starts, ends) -> tuple[np.ndarray, np.ndarray]:
        """Genome half-open intervals as stream index ranges."""
        return (
            np.searchsorted(self.gpos, np.asarray(starts, dtype=np.int64)),
            np.searchsorted(self.gpos, np.asarray(ends, dtype=np.int64)),
        )

    def paint(self, starts, ends, values, dtype=np.uint16) -> np.ndarray:
        """A per-base state: each interval painted with its value, later intervals over earlier ones."""
        out = np.zeros(self.n, dtype=dtype)
        a, b = self.stream_range(starts, ends)
        for x, y, v in zip(a.tolist(), b.tolist(), list(values), strict=True):
            if y > x:
                out[x:y] = v
        return out


# ---- contexts and counts ---------------------------------------------------------------------


def contexts(s: np.ndarray, k: int) -> np.ndarray:
    """The k bases before each position, most recent in the lowest bits; zeros before the start."""
    n = len(s)
    x = np.zeros(n, dtype=np.uint64)
    for j in range(1, min(k, n) + 1):
        x[j:] |= s[:-j].astype(np.uint64) << np.uint64(2 * (j - 1))
    return x


def rc_contexts(s: np.ndarray, k: int) -> np.ndarray:
    """The context of the reverse strand at each position: the complements of the k bases after it,
    nearest in the lowest bits. Positions within k of the end are incomplete and never used."""
    n = len(s)
    x = np.zeros(n, dtype=np.uint64)
    comp = (3 - s.astype(np.int16)).astype(np.uint64)
    for m in range(1, min(k, n - 1) + 1):
        x[: n - m] |= comp[m:] << np.uint64(2 * (m - 1))
    return x


def _hash(ctx: np.ndarray, state: np.ndarray | None, k: int, state_bits: int, bits: int) -> np.ndarray:
    """A context and its state as a key of `bits` bits: exact when it fits, multiplicative hash else.
    A collision only blurs two contexts together; the code stays decodable."""
    x = ctx if state is None else ctx | (state.astype(np.uint64) << np.uint64(2 * k))
    if 2 * k + state_bits <= bits:
        return x
    return (x * GOLD) >> np.uint64(64 - bits)


def _time_bits(n: int) -> int:
    return (2 * n + 4).bit_length()


def _group_offsets(g: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    start = np.empty(len(g), dtype=bool)
    start[0] = True
    np.not_equal(g[1:], g[:-1], out=start[1:])
    idx = np.flatnonzero(start)
    return idx, np.diff(np.append(idx, len(g)))


def _count_part(q, events, tb: int, adaptive: bool):
    tbu = np.uint64(tb)
    tmask = np.uint64((1 << tb) - 1)

    def composite(with_sym: bool) -> np.ndarray:
        rows = []
        for h, s, t in [q, *events]:
            key = h << np.uint64(2)
            if with_sym:
                key = key | s.astype(np.uint64)
            rows.append((key << tbu) | (np.uint64(1) if t is None else t))
        c = np.concatenate(rows) if rows else np.zeros(0, dtype=np.uint64)
        c.sort()
        return c

    def prior(c: np.ndarray, shift: int) -> tuple[np.ndarray, np.ndarray]:
        t = c & tmask
        ev = (t & np.uint64(1)).astype(np.int64)
        cum = np.cumsum(ev)
        idx, lens = _group_offsets(c >> np.uint64(tb + shift))
        within = cum - np.repeat(cum[idx] - ev[idx], lens)
        isq = ev == 0
        return ((t[isq] >> np.uint64(1)) - np.uint64(1)).astype(np.int64), within[isq]

    c = composite(True)
    qi, n_ca = prior(c, 0)
    if adaptive:
        c = composite(False)
        qi2, n_c = prior(c, 2)
    else:
        t = c & tmask
        ev = (t & np.uint64(1)).astype(np.int64)
        idx, lens = _group_offsets(c >> np.uint64(tb + 2))
        totals = np.repeat(np.add.reduceat(ev, idx), lens)
        isq = ev == 0
        qi2, n_c = ((t[isq] >> np.uint64(1)) - np.uint64(1)).astype(np.int64), totals[isq]
    return qi, n_ca, qi2, n_c


def _count(q_hash, q_sym, events, tb: int, adaptive: bool) -> tuple[np.ndarray, np.ndarray]:
    """For every query position i, the events with its context and base (n_ca) and with its context
    (n_c). Static events precede every query; adaptive events carry odd times and count only when
    they precede the query's time 2(i+1)."""
    n = len(q_hash)
    total = n + sum(len(e[0]) for e in events)
    pb = max(0, math.ceil(math.log2(max(1.0, total / PART_ELEMENTS))))
    hb = 62 - tb
    n_ca = np.zeros(n, dtype=np.int64)
    n_c = np.zeros(n, dtype=np.int64)
    q_time = np.arange(1, n + 1, dtype=np.uint64) << np.uint64(1)
    shift = np.uint64(max(0, hb - pb))
    for part in range(1 << pb):
        if pb:
            m = (q_hash >> shift) == np.uint64(part)
            q = (q_hash[m], q_sym[m], q_time[m])
            ev = []
            for h, s, t in events:
                em = (h >> shift) == np.uint64(part)
                ev.append((h[em], s[em], None if t is None else t[em]))
        else:
            q, ev = (q_hash, q_sym, q_time), events
        if len(q[0]) == 0:
            continue
        qi, ca, qi2, c = _count_part(q, ev, tb, adaptive)
        n_ca[qi] = ca
        n_c[qi2] = c
    return n_ca, n_c


def quantise(bits: np.ndarray) -> np.ndarray:
    return np.minimum(np.rint(np.asarray(bits, dtype=np.float64) * SCALE), MAX_CODE).astype(np.uint16)


def _codes(n_ca: np.ndarray, n_c: np.ndarray, alpha: float) -> np.ndarray:
    return quantise(-np.log2((n_ca + alpha) / (n_c + 4 * alpha)))


def alpha_for(k: int) -> float:
    return ALPHA_LOW if k <= HIGH_ORDER else ALPHA_HIGH


def static_counts(
    test_s: np.ndarray,
    train_s: np.ndarray,
    k: int,
    test_state: np.ndarray | None = None,
    train_state: np.ndarray | None = None,
    train_rc_state: np.ndarray | None = None,
    state_bits: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Counts of an order-k model (optionally conditioned on a per-base state) fitted on both strands
    of another chromosome, read at every base of the test chromosome."""
    if 2 * k + state_bits + 2 <= DIRECT_BITS:  # a direct table: no sort needed
        return _direct_counts(test_s, train_s, k, test_state, train_state, train_rc_state)
    tb = _time_bits(len(test_s))
    hb = 62 - tb
    qh = _hash(contexts(test_s, k), test_state, k, state_bits, hb)
    fh = _hash(contexts(train_s, k), train_state, k, state_bits, hb)
    events = [(fh, train_s, None)]
    m = len(train_s) - k
    if m > 0:
        rs = None if train_rc_state is None else train_rc_state[:m]
        rh = _hash(rc_contexts(train_s, k)[:m], rs, k, state_bits, hb)
        events.append((rh, (3 - train_s[:m].astype(np.int16)).astype(np.uint8), None))
    return _count(qh, test_s, events, tb, adaptive=False)


def _direct_counts(test_s, train_s, k, test_state, train_state, train_rc_state):
    def key(ctx, state):
        x = ctx.astype(np.int64)
        return x if state is None else x | (state.astype(np.int64) << (2 * k))

    fwd = key(contexts(train_s, k), train_state) * 4 + train_s
    table = np.bincount(fwd, minlength=int(fwd.max()) + 1 if len(fwd) else 1)
    m = len(train_s) - k
    if m > 0:
        rs = None if train_rc_state is None else train_rc_state[:m]
        rc = key(rc_contexts(train_s, k)[:m], rs) * 4 + (3 - train_s[:m].astype(np.int64))
        extra = np.bincount(rc)
        size = max(len(table), len(extra))
        table = np.pad(table, (0, size - len(table))) + np.pad(extra, (0, size - len(extra)))
    q = key(contexts(test_s, k), test_state)
    size = int(max(q.max() + 1, (len(table) + 3) // 4)) * 4
    table = np.pad(table, (0, size - len(table)))
    return table[q * 4 + test_s], table.reshape(-1, 4).sum(axis=1)[q]


def static_model(test_s, train_s, k, alpha=None, **state) -> np.ndarray:
    n_ca, n_c = static_counts(test_s, train_s, k, **state)
    return _codes(n_ca, n_c, alpha_for(k) if alpha is None else alpha)


def adaptive_counts(
    s: np.ndarray, k: int, revcomp: bool = True, prime: np.ndarray | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """Counts of an order-k model that learns as it codes: a base's own context and base are added
    after it is coded, and the reverse strand's context at j once base j+k is known. `prime` is
    sequence the decoder already holds (both strands counted before the first base)."""
    n = len(s)
    tb = _time_bits(n)
    hb = 62 - tb
    qh = _hash(contexts(s, k), None, k, 0, hb)
    t_fwd = (np.arange(1, n + 1, dtype=np.uint64) << np.uint64(1)) | np.uint64(1)
    events = [(qh, s, t_fwd)]
    m = n - k - 1
    if revcomp and m > 0:
        j = np.arange(m, dtype=np.uint64)
        rh = _hash(rc_contexts(s, k)[:m], None, k, 0, hb)
        rt = ((j + np.uint64(k + 1)) << np.uint64(1)) | np.uint64(1)
        events.append((rh, (3 - s[:m].astype(np.int16)).astype(np.uint8), rt))
    if prime is not None and len(prime) > k:
        events.append((_hash(contexts(prime, k), None, k, 0, hb), prime, None))
        pm = len(prime) - k
        rp = (3 - prime[:pm].astype(np.int16)).astype(np.uint8)
        events.append((_hash(rc_contexts(prime, k)[:pm], None, k, 0, hb), rp, None))
    return _count(qh, s, events, tb, adaptive=True)


def adaptive_model(s, k, alpha=None, revcomp=True, prime=None) -> np.ndarray:
    n_ca, n_c = adaptive_counts(s, k, revcomp, prime)
    return _codes(n_ca, n_c, alpha_for(k) if alpha is None else alpha)


# ---- the baselines ---------------------------------------------------------------------------


def _lgamma_fast(x: np.ndarray) -> np.ndarray:
    """log Gamma for positive arrays: exact table below 256, Stirling's series above (error < 1e-12)."""
    x = np.asarray(x, dtype=np.float64)
    out = np.empty_like(x)
    small = x < 256
    if small.any():
        out[small] = np.vectorize(math.lgamma, otypes=[np.float64])(x[small])
    big = ~small
    if big.any():
        z = x[big]
        out[big] = (z - 0.5) * np.log(z) - z + 0.5 * math.log(2 * math.pi) + 1 / (12 * z) - 1 / (360 * z**3)
    return out


def kt_bits(counts: np.ndarray, alpha: float) -> float:
    """Code length of an adaptive Dirichlet(alpha) code over rows of symbol counts (contexts x 4):
    independent of the order the bases came in, so it needs only the final counts."""
    counts = np.asarray(counts, dtype=np.float64)
    if counts.ndim == 1:
        counts = counts[None, :]
    a = counts.shape[1] * alpha
    total = counts.sum(axis=1)
    nats = (_lgamma_fast(total + a) - math.lgamma(a)).sum() - (
        _lgamma_fast(counts + alpha) - math.lgamma(alpha)
    ).sum()
    return float(nats / LN2)


def markov_row(test: Chromosome, train: Chromosome, k: int) -> dict[str, Any]:
    """One order of the baseline table, in bits per base of the test chromosome."""
    n = test.n
    key = contexts(test.s, k) * np.uint64(4) + test.s.astype(np.uint64)
    if 2 * k + 2 <= DIRECT_BITS:
        table = np.bincount(key.astype(np.int64), minlength=4 ** (k + 1)).reshape(-1, 4).astype(np.float64)
        table = table[table.sum(axis=1) > 0]
    else:
        uniq, counts = np.unique(key, return_counts=True)
        ctx = uniq >> np.uint64(2)
        sym = (uniq & np.uint64(3)).astype(np.int64)
        cidx, _ = _group_offsets(ctx)
        row_of = np.repeat(np.arange(len(cidx)), np.diff(np.append(cidx, len(ctx))))
        table = np.zeros((len(cidx), 4), dtype=np.float64)
        table[row_of, sym] = counts
        del uniq, ctx, sym, row_of
    del key
    tot = table.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        ent = -np.nansum(np.where(table > 0, table * np.log2(table / tot[:, None]), 0.0))
    adaptive = kt_bits(table, ALPHA_LOW)
    held = static_model(test.s, train.s, k).astype(np.float64).sum() / SCALE
    return {
        "order": k,
        "contexts_seen": len(table),
        "in_sample_bits_per_base": round(ent / n, 4),
        "adaptive_bits_per_base": round(adaptive / n, 4),
        "held_out_bits_per_base": round(held / n, 4),
    }


def pack2(s: np.ndarray) -> bytes:
    """Four bases per byte."""
    pad = (-len(s)) % 4
    x = np.concatenate([s, np.zeros(pad, dtype=np.uint8)]).reshape(-1, 4)
    return ((x[:, 0] << 6) | (x[:, 1] << 4) | (x[:, 2] << 2) | x[:, 3]).astype(np.uint8).tobytes()


def general_compressors(s: np.ndarray, ascii_xz: bool = True) -> dict[str, Any]:
    """xz and bzip2 on the 2-bit packed bases, and xz on the letters, in bits per base."""
    n = len(s)
    packed = pack2(s)
    out: dict[str, Any] = {}
    t = time.perf_counter()
    out["xz_packed"] = {
        "bits_per_base": round(len(lzma.compress(packed, preset=9 | lzma.PRESET_EXTREME)) * 8 / n, 4)
    }
    out["xz_packed"]["seconds"] = round(time.perf_counter() - t, 1)
    t = time.perf_counter()
    out["bzip2_packed"] = {"bits_per_base": round(len(bz2.compress(packed, 9)) * 8 / n, 4)}
    out["bzip2_packed"]["seconds"] = round(time.perf_counter() - t, 1)
    if ascii_xz:
        letters = np.frombuffer(b"ACGT", dtype=np.uint8)[s].tobytes()
        t = time.perf_counter()
        out["xz_letters"] = {"bits_per_base": round(len(lzma.compress(letters, preset=9)) * 8 / n, 4)}
        out["xz_letters"]["seconds"] = round(time.perf_counter() - t, 1)
    return out


# ---- the mixture -----------------------------------------------------------------------------


def mix(
    models: list[np.ndarray],
    window: int,
    beta: float,
    floor: float = MIX_FLOOR,
    chunk: int = MIX_CHUNK,
) -> np.ndarray:
    """Per-base code length (bits, float32) of a causal mixture: each model weighted by
    2^(-beta * its code length over the previous `window` bases), with a floor of the uniform.

    A model whose code is INACTIVE at a base (a layer outside the region its annotation claims,
    which the decoder knows because the annotation is transmitted) takes no weight there, and
    counts as two bits in the window sums that later bases read."""
    m = len(models)
    n = len(models[0])
    out = np.empty(n, dtype=np.float32)
    for a in range(0, n, chunk):
        b = min(n, a + chunk)
        a0 = max(0, a - window)
        raw = np.stack([x[a0:b] for x in models])
        inactive = raw == INACTIVE
        lengths = np.where(inactive, 2.0, raw / SCALE)
        cum = np.zeros((m, b - a0 + 1), dtype=np.float64)
        np.cumsum(lengths, axis=1, out=cum[:, 1:])
        r = np.arange(a - a0, b - a0)
        logw = cum[:, np.maximum(0, r - window)] - cum[:, r]
        logw *= beta * LN2
        now = inactive[:, a - a0 :]
        logw[now] = -np.inf
        top = logw.max(axis=0)
        top[~np.isfinite(top)] = 0.0
        w = np.exp(logw - top)
        p = (w * np.power(2.0, -lengths[:, a - a0 :])).sum(axis=0)
        total = w.sum(axis=0)
        p = np.divide(p, total, out=np.full_like(p, 0.25), where=total > 0)
        out[a:b] = -np.log2((1 - floor) * p + floor / 4)
    return out


def restrict(codes: np.ndarray, active: np.ndarray) -> np.ndarray:
    """A model made inactive outside the bases it claims."""
    codes[~active] = INACTIVE
    return codes


# ---- the costs of an annotation --------------------------------------------------------------


def delta_bits(values) -> float:
    """Elias delta code lengths of positive integers, summed."""
    v = np.maximum(1, np.asarray(values, dtype=np.float64))
    lg = np.floor(np.log2(v))
    return float((lg + 2 * np.floor(np.log2(lg + 1)) + 1).sum())


def integer_bits(values) -> float:
    """Positive integers under an adaptive code: each value's bit length from an adaptive categorical
    code over lengths, then the bits below the leading one written out."""
    v = np.maximum(1, np.asarray(values, dtype=np.int64))
    if len(v) == 0:
        return 0.0
    lengths = np.floor(np.log2(v)).astype(np.int64) + 1
    return categorical_bits(lengths, 64) + float((lengths - 1).sum())


def categorical_bits(ids, alphabet: int, alpha: float = ALPHA_LOW) -> float:
    """An adaptive code for a sequence of labels drawn from `alphabet` values."""
    ids = np.asarray(ids, dtype=np.int64)
    if len(ids) == 0 or alphabet <= 1:
        return 0.0
    return kt_bits(np.bincount(ids, minlength=alphabet)[None, :], alpha)


def interval_bits(starts, ends) -> float:
    """Intervals sorted by start: each start as the gap from the previous start, each length."""
    starts = np.asarray(starts, dtype=np.int64)
    ends = np.asarray(ends, dtype=np.int64)
    if len(starts) == 0:
        return 0.0
    order = np.argsort(starts, kind="stable")
    st = starts[order]
    gaps = np.diff(np.concatenate([[0], st])) + 1
    return integer_bits(gaps) + integer_bits(ends[order] - st + 1)


def names_bits(names) -> float:
    """A label dictionary written out as text, a byte per character and a terminator. Reported, not
    charged: the dictionary is part of the models fitted on the other chromosome."""
    return float(sum(8 * (len(x) + 1) for x in names))


# ---- the layers ------------------------------------------------------------------------------


@dataclass
class Layer:
    name: str
    title: str
    models: dict[str, np.ndarray]
    claim: np.ndarray
    annotation_bits: float
    items: int
    notes: dict[str, Any] = field(default_factory=dict)


def _bed_rows(path: Path, ledger: Ledger | None):
    if not path.exists():
        return []
    raw = ledger.read(path) if ledger else path.read_bytes()
    text = gzip.decompress(raw).decode() if path.suffix == ".gz" else raw.decode()
    return [ln.split("\t") for ln in text.splitlines() if ln and not ln.startswith("#")]


def _states(test: Chromosome, train: Chromosome, fn) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """A per-base state on the test chromosome, and on the training set with its reverse-strand view.
    `fn(chromosome)` returns (state, reverse-strand state) for one chromosome."""
    t = fn(test)[0]
    parts = [fn(p) for p in train.parts]
    return t, np.concatenate([p[0] for p in parts]), np.concatenate([p[1] for p in parts])


def _state_models(test, train, name, states, bits, orders, active=None) -> dict[str, np.ndarray]:
    """Static models over one state, for several orders, inactive outside `active`."""
    t_state, r_state, r_rc_state = states
    out = {}
    for k in orders:
        codes = static_model(
            test.s,
            train.s,
            k,
            test_state=t_state,
            train_state=r_state,
            train_rc_state=r_rc_state,
            state_bits=bits,
        )
        out[f"{name}{k}"] = codes if active is None else restrict(codes, active)
    return out


def _rmsk(c: Chromosome, results_dir: Path, ledger) -> list[list[str]]:
    return _bed_rows(results_dir / f"rmsk_{c.name}.bed.gz", ledger)


def repeats_layer(test: Chromosome, train: Chromosome, results_dir: Path, ledger=None) -> Layer:
    """Every RepeatMasker copy as a copy of its subfamily, and of its family at its divergence."""
    rows = {c.name: _rmsk(c, results_dir, ledger) for c in (test, *train.parts)}
    subs = sorted({r[4] for rr in rows.values() for r in rr})
    fams = sorted({f"{r[2]}/{r[3]}" for rr in rows.values() for r in rr})
    sub_id = {x: i for i, x in enumerate(subs)}
    fam_id = {x: i for i, x in enumerate(fams)}

    def div_bin(x: str) -> int:
        return int(np.searchsorted(DIV_EDGES, float(x), side="right"))

    def sub_state(c):
        rr = rows[c.name]
        st = c.paint([int(r[0]) for r in rr], [int(r[1]) for r in rr], [sub_id[r[4]] + 1 for r in rr])
        return st, st

    def fd_state(c):
        rr = rows[c.name]
        values = [fam_id[f"{r[2]}/{r[3]}"] * 5 + div_bin(r[5]) + 1 for r in rr]
        st = c.paint([int(r[0]) for r in rr], [int(r[1]) for r in rr], values)
        return st, st

    subs_st = _states(test, train, sub_state)
    claim = subs_st[0] > 0
    models = _state_models(
        test, train, "subfamily", subs_st, (len(subs) + 1).bit_length(), (8, 12, 16), claim
    )
    models |= _state_models(
        test,
        train,
        "family_divergence",
        _states(test, train, fd_state),
        (5 * len(fams) + 1).bit_length(),
        (11,),
        claim,
    )
    rr = rows[test.name]
    cost = (
        interval_bits([int(r[0]) for r in rr], [int(r[1]) for r in rr])
        + categorical_bits([sub_id[r[4]] for r in rr], len(subs))
        + categorical_bits([div_bin(r[5]) for r in rr], 5)
    )
    tandem = [r for r in rr if r[2] in TANDEM_CLASSES]
    return Layer(
        "repeats",
        "RepeatMasker subfamilies with divergence",
        models,
        claim,
        cost,
        len(rr),
        {
            "subfamilies": len(subs),
            "families": len(fams),
            "dictionary_bits_not_charged": round(names_bits(subs) + names_bits(fams)),
            "tandem": test.paint(
                [int(r[0]) for r in tandem], [int(r[1]) for r in tandem], [1] * len(tandem), np.uint8
            ),
        },
    )


def _window_keys(a: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    """For each position x, the k bases before it as a key, and whether all k are bases."""
    n = len(a)
    key = np.zeros(n, dtype=np.uint64)
    bad = np.zeros(n, dtype=bool)
    bad[: min(k, n)] = True
    for j in range(1, min(k, n) + 1):
        key[j:] |= np.minimum(a[:-j], 3).astype(np.uint64) << np.uint64(2 * (j - 1))
        bad[j:] |= a[:-j] > 3
    return key, ~bad


def copy_offsets(t: np.ndarray, p: np.ndarray, k: int = COPY_K) -> np.ndarray:
    """For each position x of `t`, the position of `p` it is predicted to copy (-1 where none).

    Causal: the prediction for t[x] comes from the last position x' <= x whose preceding k bases
    t[x'-k:x'] occur in p, at the occurrence nearest the linear expectation; the diagonal found
    there is carried forward over mismatches until the next anchor."""
    lt, lp = len(t), len(p)
    out = np.full(lt, -1, dtype=np.int64)
    if lt <= k or lp <= k:
        return out
    kt, vt = _window_keys(t, k)
    kp, vp = _window_keys(p, k)
    ys = np.flatnonzero(vp)
    xs = np.flatnonzero(vt)
    if len(ys) == 0 or len(xs) == 0:
        return out
    table = np.sort((kp[ys] << np.uint64(32)) | ys.astype(np.uint64))
    y_exp = np.minimum(lp - 1, (xs * (lp / lt)).astype(np.int64))
    i = np.searchsorted(table, (kt[xs] << np.uint64(32)) | y_exp.astype(np.uint64))
    best = np.full(len(xs), -1, dtype=np.int64)
    best_d = np.full(len(xs), np.iinfo(np.int64).max, dtype=np.int64)
    for j in (i - 1, i):
        ok = (j >= 0) & (j < len(table))
        cand = table[np.clip(j, 0, len(table) - 1)]
        ok &= (cand >> np.uint64(32)) == kt[xs]
        y = (cand & np.uint64(0xFFFFFFFF)).astype(np.int64)
        d = np.abs(y - y_exp)
        better = ok & (d < best_d)
        best[better] = y[better]
        best_d[better] = d[better]
    found = (best >= 0) & (best_d <= max(500, int(0.2 * max(lt, lp))))
    diag = np.zeros(lt, dtype=np.int64)
    last = np.full(lt, -1, dtype=np.int64)
    diag[xs[found]] = best[found] - xs[found]
    last[xs[found]] = xs[found]
    np.maximum.accumulate(last, out=last)
    x = np.flatnonzero(last >= 0)
    y = x + diag[last[x]]
    inside = (y >= 0) & (y < lp)
    out[x[inside]] = y[inside]
    return out


def copy_codes(t: np.ndarray, p: np.ndarray, offsets: np.ndarray, usable: np.ndarray) -> tuple:
    """Code lengths of the copied bases: the probability of a copy is its pair's running agreement
    (Krichevsky-Trofimov over the predictions made so far), the rest shared by the other three."""
    x = np.flatnonzero(usable)
    pred = p[offsets[x]]
    ok = pred < 4
    x, pred = x[ok], pred[ok]
    hit = (pred == t[x]).astype(np.int64)
    before = np.cumsum(hit) - hit
    q = (before + 0.5) / (np.arange(len(x)) + 1.0)
    bits = np.where(hit == 1, -np.log2(q), -np.log2((1 - q) / 3))
    return x, bits, int(hit.sum())


def segdup_layer(
    test: Chromosome, results_dir: Path, reference: Path = REFERENCE, ledger=None, partners=None
) -> Layer:
    """Segmental duplications as copies of a locus the decoder already holds: a partner on an
    earlier chromosome in karyotype order, or earlier on the same chromosome."""
    rows = _bed_rows(results_dir / f"superdups_{test.name}.bed.gz", ledger)
    order = {c: i for i, c in enumerate(KARYOTYPE)}
    here = order.get(test.name, len(KARYOTYPE))
    allowed = [
        r
        for r in rows
        if r[2] in order and (order[r[2]] < here or (r[2] == test.name and int(r[3]) < int(r[0])))
    ]
    allowed.sort(key=lambda r: float(r[5]))  # the most similar partner is painted last and wins
    codes = np.full(test.n, INACTIVE, dtype=np.uint16)
    claim = np.zeros(test.n, dtype=bool)
    partners = partners if partners is not None else {}
    lengths: dict[str, int] = {}
    by_chrom: dict[str, list] = {}
    for r in allowed:
        by_chrom.setdefault(r[2], []).append(r)
    predicted = hits = 0
    used = []
    primes = []  # partner sequence on other chromosomes: what a genome-wide generic coder would hold
    same = np.zeros(test.n, dtype=bool)
    for pc in sorted(by_chrom, key=lambda c: order[c]):
        genome = test.genome if pc == test.name else partners.get(pc)
        if genome is None:
            p = reference / f"{pc}.fa"
            if not p.exists():
                continue
            genome = parse_fasta_codes(ledger.read(p) if ledger else p.read_bytes())
        lengths[pc] = len(genome)
        for r in by_chrom[pc]:
            ta, tb_, pa, pb_ = int(r[0]), int(r[1]), int(r[3]), int(r[4])
            minus = r[6] == "-"
            t = test.genome[ta:tb_]
            p = genome[pa:pb_]
            if minus:
                p = np.where(p < 4, 3 - p, 4).astype(np.uint8)[::-1]
            off = copy_offsets(t, p)
            usable = (off >= 0) & (t < 4)
            if pc == test.name:  # the copied base must already be decoded
                src = np.where(minus, pb_ - 1 - off, pa + off)
                usable &= src < ta + np.arange(len(t))
            x, bits, h = copy_codes(t, p, off, usable)
            if len(x) == 0:
                continue
            codes[np.searchsorted(test.gpos, ta + x)] = quantise(bits)
            a, b = test.stream_range([ta], [tb_])
            claim[a[0] : b[0]] = True
            if pc == test.name:
                same[a[0] : b[0]] = True
            else:
                primes.append(p[p < 4])
            predicted += len(x)
            hits += h
            used.append(r)
        del genome
    cost = 0.0
    if used:
        cost = interval_bits([int(r[0]) for r in used], [int(r[1]) for r in used])
        cost += len(used) * (math.log2(len(KARYOTYPE)) + 1)  # partner chromosome and strand
        cost += sum(math.log2(max(2, lengths[r[2]])) for r in used)  # partner start
        cost += integer_bits([int(r[4]) - int(r[3]) + 1 for r in used])
    return Layer(
        "duplications",
        "segmental duplications as copies of a decoded locus",
        {"copy": codes},
        claim,
        cost,
        len(used),
        {
            "pairs": len(rows),
            "pairs_with_decoded_partner": len(allowed),
            "pairs_used": len(used),
            "claimed_bases_with_a_partner_on_the_same_chromosome": int((same & claim).sum()),
            "partner_bases_on_other_chromosomes": int(sum(len(x) for x in primes)),
            "prime": np.concatenate(primes) if primes else None,
            "predictions_made": predicted,
            "prediction_agreement": round(hits / predicted, 4) if predicted else None,
            "any": test.paint(
                [int(r[0]) for r in rows], [int(r[1]) for r in rows], [1] * len(rows), np.uint8
            ),
        },
    )


def gc_state(s: np.ndarray, window: int = GC_WINDOW, causal: bool = True) -> np.ndarray:
    """GC share and CpG observed/expected over the `window` bases before each base (causal), or after
    it (the reverse strand's view, for training events): 18 states, never 0."""
    n = len(s)

    def cum(a: np.ndarray) -> np.ndarray:
        out = np.zeros(n + 1, dtype=np.int64)
        np.cumsum(a, out=out[1:])
        return out

    c, g = s == 1, s == 2
    cg = np.zeros(n, dtype=bool)
    cg[1:] = c[:-1] & g[1:]  # the dinucleotide ending at each base
    gc_c, c_c, g_c, cg_c = cum(c | g), cum(c), cum(g), cum(cg)
    i = np.arange(n, dtype=np.int64)
    if causal:
        lo, hi = np.maximum(0, i - window), i
    else:
        lo, hi = np.minimum(n, i + 1), np.minimum(n, i + 1 + window)
    length = hi - lo
    gc = (gc_c[hi] - gc_c[lo]) / np.maximum(1, length)
    cpg = cg_c[hi] - cg_c[np.minimum(lo + 1, hi)]
    expected = (c_c[hi] - c_c[lo]) * (g_c[hi] - g_c[lo]) / np.maximum(1, length)
    oe = np.where(expected > 0, cpg / np.maximum(expected, 1e-9), 0.0)
    gb = np.searchsorted(np.asarray(GC_EDGES), gc, side="right")
    ob = np.searchsorted(np.asarray(OE_EDGES), oe, side="right")
    return (gb * 3 + ob + 1).astype(np.uint8)


def cpg_islands(
    chrom: str, cache: Path = CACHE, ledger: Ledger | None = None
) -> list[tuple[int, int]] | None:
    """UCSC's CpG islands for one chromosome: one REST request, cached as a BED."""
    p = cache / f"cpgIslandExt_{chrom}.bed.gz"
    if not p.exists():
        url = f"{UCSC_API}?genome=hg38;track=cpgIslandExt;chrom={chrom}"
        req = urllib.request.Request(url, headers={"User-Agent": "GenomeOS/0.9 (compression probe)"})
        try:
            with urllib.request.urlopen(req, timeout=180) as r:  # noqa: S310
                raw = r.read()
        except OSError:
            return None
        if ledger:
            ledger.requests += 1
            ledger.network_bytes += len(raw)
        items = json.loads(raw).get("cpgIslandExt", [])
        if isinstance(items, dict):
            items = items.get(chrom, [])
        cache.mkdir(parents=True, exist_ok=True)
        with gzip.open(p, "wt") as fh:
            fh.write(f"# UCSC hg38 cpgIslandExt, {chrom}, fetched by GenomeOS through {UCSC_API}\n")
            for it in sorted(items, key=lambda x: x["chromStart"]):
                fh.write(f"{chrom}\t{it['chromStart']}\t{it['chromEnd']}\n")
    return [(int(r[1]), int(r[2])) for r in _bed_rows(p, ledger)]


def gc_layer(test: Chromosome, train: Chromosome, ledger=None, cache: Path = CACHE) -> Layer:
    """GC and CpG structure: the composition of the kilobase just coded (no annotation, active
    everywhere), and the CpG islands (annotation, active inside them)."""
    models = _state_models(
        test,
        train,
        "gc_cpg",
        _states(test, train, lambda c: (gc_state(c.s), gc_state(c.s, causal=False))),
        5,
        (3, 6, 10),
    )
    islands = {c.name: cpg_islands(c.name, cache, ledger) for c in (test, *train.parts)}
    cost = 0.0
    notes: dict[str, Any] = {"cpg_islands": None}
    if all(v is not None for v in islands.values()):

        def island_state(c):
            iv = islands[c.name]
            st = c.paint([a for a, _ in iv], [b for _, b in iv], [1] * len(iv), np.uint8)
            return st, st

        st = _states(test, train, island_state)
        models |= _state_models(test, train, "cpg_island", st, 1, (4, 8), st[0] > 0)
        iv = islands[test.name]
        cost = interval_bits([a for a, _ in iv], [b for _, b in iv])
        notes = {"cpg_islands": len(iv), "island_mask": st[0] > 0}
    return Layer(
        "gc_cpg",
        "GC and CpG composition, CpG islands",
        models,
        np.ones(test.n, dtype=bool),
        cost,
        len(islands[test.name] or []),
        notes,
    )


def canonical_cds(chrom: str, ledger=None) -> list[tuple[bool, list[tuple[int, int, int]]]]:
    """(minus strand, [(start, end, phase)]) of the canonical coding transcript of every coding gene."""
    from genomeos.genome.annotation import Annotation, default_gencode

    gff = default_gencode({chrom})
    if gff is None:
        return []
    if ledger:
        ledger.disk_bytes += gff.stat().st_size
    ann = Annotation.from_gff3(gff, {chrom})
    out = []
    for g in ann.genes.values():
        if g.locus.chrom != chrom or g.type != "protein_coding":
            continue
        ts = [t for t in g.transcripts.values() if t.cds]
        if not ts:
            continue
        canon = [t for t in ts if "Ensembl_canonical" in t.tags] or sorted(
            ts, key=lambda t: -sum(c.end - c.start for c, _ in t.cds)
        )
        t = canon[0]
        out.append((t.locus.strand.value == "-", sorted((c.start, c.end, int(ph)) for c, ph in t.cds)))
    out.sort(key=lambda x: x[1][0][0])
    return out


def coding_states(c: Chromosome, transcripts) -> tuple[np.ndarray, np.ndarray]:
    """Per base: 1 + 3 * minus + codon position in the transcript's direction, and the same state as
    the reverse strand sees it (strand flipped, codon position kept)."""
    g = np.zeros(len(c.genome), dtype=np.uint8)
    for minus, segs in transcripts:
        for a, b, ph in segs:
            o = np.arange(b - a)
            if minus:
                o = o[::-1]
            g[a:b] = 1 + 3 * int(minus) + (o - ph) % 3
    fwd = g[c.gpos]
    rc = np.where(fwd > 0, np.where(fwd > 3, fwd - 3, fwd + 3), 0).astype(np.uint8)
    return fwd, rc


def coding_layer(test: Chromosome, train: Chromosome, ledger=None) -> Layer:
    """GENCODE coding sequence with its reading frame: an order-k model per strand and codon position."""
    txs = {c.name: canonical_cds(c.name, ledger) for c in (test, *train.parts)}
    st = _states(test, train, lambda c: coding_states(c, txs[c.name]))
    claim = st[0] > 0
    models = _state_models(test, train, "codon", st, 3, (2, 5, 8), claim)
    tt = txs[test.name]
    cost = 0.0
    if tt:
        cost = integer_bits(np.diff([0] + [segs[0][0] for _, segs in tt]) + 1)  # transcript starts
        cost += len(tt) * 3  # strand and the first segment's phase
        cost += integer_bits([len(segs) for _, segs in tt])
        introns, exons = [], []
        for _, segs in tt:
            introns += [a - e + 1 for (a, _, _), (_, e, _) in zip(segs[1:], segs[:-1], strict=True)]
            exons += [b - a for a, b, _ in segs]
        cost += integer_bits(introns) + integer_bits(exons)
    plus_starts = [segs[0][0] for minus, segs in tt if not minus]
    atg = sum(1 for a in plus_starts if test.genome[a : a + 3].tolist() == [0, 3, 2])
    return Layer(
        "coding",
        "GENCODE canonical CDS with codon phase",
        models,
        claim,
        cost,
        len(tt),
        {
            "transcripts": len(tt),
            "plus_strand_starting_atg": round(atg / len(plus_starts), 3) if plus_starts else None,
        },
    )


def ccre_layer(test: Chromosome, train: Chromosome, results_dir: Path, ledger=None) -> Layer:
    """ENCODE cCREs: an order-k model per registry class."""
    rows = {c.name: _bed_rows(results_dir / f"ccres_{c.name}.bed.gz", ledger) for c in (test, *train.parts)}
    classes = sorted({r[4] for rr in rows.values() for r in rr})
    cid = {x: i + 1 for i, x in enumerate(classes)}

    def state(c):
        rr = rows[c.name]
        st = c.paint([int(r[1]) for r in rr], [int(r[2]) for r in rr], [cid[r[4]] for r in rr], np.uint8)
        return st, st

    st = _states(test, train, state)
    claim = st[0] > 0
    models = _state_models(test, train, "ccre", st, (len(classes) + 1).bit_length(), (4, 8), claim)
    rr = rows[test.name]
    cost = interval_bits([int(r[1]) for r in rr], [int(r[2]) for r in rr])
    cost += categorical_bits([cid[r[4]] - 1 for r in rr], len(classes))
    return Layer("registry", "ENCODE cCREs by class", models, claim, cost, len(rr), {"classes": len(classes)})


def motif_layer(test: Chromosome, results_dir: Path, ledger=None) -> Layer | None:
    """The JASPAR sites the motif runs recorded, each base predicted by its matrix column."""
    from genomeos.genome.motifs import JASPAR_PATH, load_motifs, promoter_loci

    d = load_result(f"motifs_{test.name}", results_dir)
    if not d or not JASPAR_PATH.exists():
        return None
    if ledger:
        ledger.disk_bytes += JASPAR_PATH.stat().st_size
    by_name: dict[str, Any] = {}
    for m in load_motifs():
        by_name.setdefault(m.name, m)
    loci = {sym: a for sym, a, _ in promoter_loci(test.name, d.get("flank", 1000)) or []}
    sites = []
    for sym, g in d.get("genes", {}).items():
        if sym in loci:
            sites += [(loci[sym] + r["position"], r["factor"], r["strand"]) for r in g.get("requires", [])]
    for e in d.get("elements", []):
        sites += [(e["start"] + r["position"], r["factor"], r["strand"]) for r in e.get("requires", [])]
    sites = sorted({s for s in sites if s[1] in by_name})
    codes = np.full(test.n, INACTIVE, dtype=np.uint16)
    claim = np.zeros(test.n, dtype=bool)
    recorded_better = 0
    for f, name, strand in sites:
        m = by_name[name]
        w = m.width
        counts = np.asarray(m.counts, dtype=np.float64)
        prob = (counts + 0.25) / (counts.sum(axis=0) + 1.0)
        bases = test.genome[f : f + w]
        if len(bases) < w or (bases > 3).any():
            continue
        cols = np.arange(w)
        plus = prob[bases, cols]
        minus = prob[3 - bases[::-1], cols][::-1]
        recorded, other = (minus, plus) if strand in ("-", 1) else (plus, minus)
        recorded_better += int(np.log2(recorded).sum() >= np.log2(other).sum())
        idx = np.searchsorted(test.gpos, f + cols)
        codes[idx] = quantise(-np.log2(recorded))
        claim[idx] = True
    factors = sorted({s[1] for s in sites})
    fid = {x: i for i, x in enumerate(factors)}
    cost = 0.0
    if sites:
        cost = integer_bits(np.diff([0] + [s[0] for s in sites]) + 1) + len(sites)  # position, strand
        cost += categorical_bits([fid[s[1]] for s in sites], len(factors))
    return Layer(
        "motifs",
        "JASPAR sites recorded by the motif runs",
        {"pwm": codes},
        claim,
        cost,
        len(sites),
        {
            "sites": len(sites),
            "factors": len(factors),
            "recorded_strand_scores_at_least_the_other": round(recorded_better / len(sites), 3)
            if sites
            else None,
        },
    )


def tier_blocks(chrom: str, results_dir: Path) -> list[dict[str, Any]]:
    budget = load_result(f"budget_{chrom}", results_dir) or {}
    unknown = {
        (b["start"], b["end"]): b["class"]
        for b in (load_result(f"unknown_{chrom}", results_dir) or {}).get("blocks", [])
    }
    out = []
    for b in budget.get("blocks", []):
        tier = (b.get("guess") or {}).get("tier")
        if tier in TIERS:
            out.append(
                {
                    "start": b["start"],
                    "end": b["end"],
                    "tier": tier,
                    "class": unknown.get((b["start"], b["end"]), b.get("class")),
                }
            )
    return sorted(out, key=lambda b: b["start"])


def tiers_layer(test: Chromosome, train: Chromosome, results_dir: Path) -> Layer:
    """The UNKNOWN blocks: an order-k model per budget tier and per sequence class."""
    blocks = {c.name: tier_blocks(c.name, results_dir) for c in (test, *train.parts)}
    classes = sorted({b["class"] for bb in blocks.values() for b in bb})
    cid = {x: i + 1 for i, x in enumerate(classes)}

    def tier_state(c):
        bb = blocks[c.name]
        values = [TIERS.index(b["tier"]) + 1 for b in bb]
        st = c.paint([b["start"] for b in bb], [b["end"] for b in bb], values, np.uint8)
        return st, st

    def class_state(c):
        bb = blocks[c.name]
        st = c.paint(
            [b["start"] for b in bb], [b["end"] for b in bb], [cid[b["class"]] for b in bb], np.uint8
        )
        return st, st

    tiers = _states(test, train, tier_state)
    claim = tiers[0] > 0
    models = _state_models(test, train, "tier", tiers, 3, (4, 8, 12), claim)
    cbits = (len(classes) + 1).bit_length()
    models |= _state_models(test, train, "block_class", _states(test, train, class_state), cbits, (8,), claim)
    bb = blocks[test.name]
    cost = interval_bits([b["start"] for b in bb], [b["end"] for b in bb])
    cost += categorical_bits([TIERS.index(b["tier"]) for b in bb], len(TIERS))
    cost += categorical_bits([cid[b["class"]] - 1 for b in bb], len(classes))
    block_id = test.paint(
        [b["start"] for b in bb], [b["end"] for b in bb], list(range(1, len(bb) + 1)), np.int32
    )
    return Layer(
        "tiers",
        "UNKNOWN blocks by budget tier and class",
        models,
        claim,
        cost,
        len(bb),
        {"blocks": bb, "tier_state": tiers[0], "block_id": block_id},
    )


# ---- stacks, regions and the summary ---------------------------------------------------------

COMBO_SHAPE = (len(TIERS) + 1, 2, 3, 2, len(GC_EDGES) + 1)  # tier, coding, repeat kind, duplicated, GC bin


def window_gc_bins(c: Chromosome, window: int = 1_000) -> np.ndarray:
    """The GC bin of the (non-causal) kilobase window each base sits in: the matching variable."""
    w = c.gpos // window
    gc = np.bincount(w, weights=((c.s == 1) | (c.s == 2)).astype(np.float64)) / np.maximum(1, np.bincount(w))
    return np.searchsorted(np.asarray(GC_EDGES), gc, side="right")[w].astype(np.int64)


def stack_name(key: tuple[str, ...]) -> str:
    return "+".join(key)


def plan_stacks(layers: list[str], control: bool = False) -> dict[str, list[tuple[str, ...]]]:
    """The stacks to mix, by the question each answers. The control gives the generic models the
    duplication partners' sequence without the alignment, to price the alignment alone."""
    naive, generic = ("naive",), ("naive", "adaptive")
    aware = generic + ("repeats",)
    full = generic + tuple(layers)
    plan = {
        "baselines": [naive, generic],
        "over_naive": [naive + (x,) for x in layers],
        "over_generic": [generic + (x,) for x in layers],
        "over_repeat_aware": [aware + (x,) for x in layers if x != "repeats"],
        "cumulative": [generic + tuple(layers[: i + 1]) for i in range(len(layers))],
        "leave_one_out": [tuple(g for g in full if g != x) for x in layers],
    }
    if control:
        primed = ("naive", "adaptive_primed")
        plan["control"] = [primed, primed + ("duplications",)]
    return plan


def summarise_regions(aggs, bases_combo, bases_block, blocks, n, rng) -> dict[str, Any]:
    """Bits per base by tier, all bases and unique bases, GC-standardised, with the block bootstrap."""
    tier, cds, rep, dup, gcb = np.unravel_index(np.arange(np.prod(COMBO_SHAPE)), COMBO_SHAPE)
    unique = (rep == 0) & (dup == 0)

    def bpb(key, sel) -> float | None:
        bp = bases_combo[sel].sum()
        return round(float(aggs[key]["combo"][sel].sum() / bp), 4) if bp else None

    keys = list(aggs)
    neutral_unique = (tier == TIERS.index("neutral") + 1) & unique
    weights = np.array(
        [bases_combo[neutral_unique & (gcb == b)].sum() for b in range(COMBO_SHAPE[4])], dtype=float
    )

    def standardised(key, sel) -> float | None:
        num = den = 0.0
        for b in range(COMBO_SHAPE[4]):
            s = sel & (gcb == b)
            bp = bases_combo[s].sum()
            if bp and weights[b]:
                num += weights[b] * aggs[key]["combo"][s].sum() / bp
                den += weights[b]
        return round(num / den, 4) if den else None

    regions = {
        "chromosome": np.ones_like(tier, dtype=bool),
        "coding_cds": cds == 1,
        "outside_unknown_non_coding": (tier == 0) & (cds == 0),
        **{t: tier == i + 1 for i, t in enumerate(TIERS)},
    }
    rows = []
    for name, sel in regions.items():
        bp = int(bases_combo[sel].sum())
        if not bp:
            continue
        rows.append(
            {
                "region": name,
                "bases": bp,
                "interspersed_share": round(float(bases_combo[sel & (rep == 1)].sum() / bp), 4),
                "tandem_share": round(float(bases_combo[sel & (rep == 2)].sum() / bp), 4),
                "duplicated_share": round(float(bases_combo[sel & (dup == 1)].sum() / bp), 4),
                "unique_bases": int(bases_combo[sel & unique].sum()),
                "bits_per_base": {k: bpb(k, sel) for k in keys},
                "unique_bits_per_base": {k: bpb(k, sel & unique) for k in keys},
                "interspersed_bits_per_base": {k: bpb(k, sel & (rep == 1)) for k in keys},
                "tandem_bits_per_base": {k: bpb(k, sel & (rep == 2)) for k in keys},
                "duplicated_bits_per_base": {k: bpb(k, sel & (dup == 1)) for k in keys},
                "unique_bits_per_base_gc_standardised": {k: standardised(k, sel & unique) for k in keys},
            }
        )
    # constrained_unknown against neutral, unique bases, resampling blocks
    boot = {}
    ids = {t: [i + 1 for i, b in enumerate(blocks) if b["tier"] == t] for t in TIERS}
    for key in keys:
        per = aggs[key]["block"]

        def tier_arrays(t, per=per):
            idx = np.asarray(ids[t], dtype=np.int64)
            return per[idx * 2 + 1], bases_block[idx * 2 + 1]

        (cb, cn), (nb, nn) = tier_arrays("constrained_unknown"), tier_arrays("neutral")
        (fb, fn) = tier_arrays("fossil")
        if cn.sum() == 0 or nn.sum() == 0:
            continue
        diffs, fdiffs = [], []
        for _ in range(BOOTSTRAP):
            a = rng.integers(0, len(cb), len(cb))
            b = rng.integers(0, len(nb), len(nb))
            if cn[a].sum() and nn[b].sum():
                diffs.append(cb[a].sum() / cn[a].sum() - nb[b].sum() / nn[b].sum())
            if len(fb) and fn.sum():
                f = rng.integers(0, len(fb), len(fb))
                if fn[f].sum() and nn[b].sum():
                    fdiffs.append(fb[f].sum() / fn[f].sum() - nb[b].sum() / nn[b].sum())
        boot[key] = {
            "constrained_minus_neutral": round(float(cb.sum() / cn.sum() - nb.sum() / nn.sum()), 4),
            "ci95": [round(float(np.percentile(diffs, 2.5)), 4), round(float(np.percentile(diffs, 97.5)), 4)],
            "fossil_minus_neutral": round(float(fb.sum() / fn.sum() - nb.sum() / nn.sum()), 4)
            if fn.sum()
            else None,
            "fossil_ci95": [
                round(float(np.percentile(fdiffs, 2.5)), 4),
                round(float(np.percentile(fdiffs, 97.5)), 4),
            ]
            if fdiffs
            else None,
            "blocks": {
                "constrained_unknown": int((cn > 0).sum()),
                "neutral": int((nn > 0).sum()),
                "fossil": int((fn > 0).sum()),
            },
        }
    return {"rows": rows, "unique_bootstrap": boot}


def run(
    chrom: str,
    train: str,
    results_dir: Path = RESULTS_DIR,
    reference: Path = REFERENCE,
    progress=None,
    markov: bool = True,
    compressors: bool = True,
    spill: Path | None = None,
) -> dict[str, Any]:
    """Code `chrom` under every stack. `spill` keeps model codes in files there instead of memory
    (2 bytes a base a model), for chromosomes whose models do not fit."""
    ledger = Ledger(progress)
    test = Chromosome.load(chrom, reference, ledger)
    trn = Chromosome.training_set([Chromosome.load(c, reference, ledger) for c in train.split(",")])
    n = test.n
    gaps = np.flatnonzero(np.diff(np.concatenate([[0], (test.genome > 3).astype(np.int8), [0]])))
    gap_bits = interval_bits(gaps[0::2], gaps[1::2]) if len(gaps) else 0.0
    ledger.lap("sequence")
    out: dict[str, Any] = {
        "chrom": chrom,
        "trained_on": train,
        "bases_coded": n,
        "other_letters": int(len(test.genome) - n),
        "gap_description_bits": round(gap_bits, 1),
        "evidence": EVIDENCE,
    }
    if markov:
        table = [markov_row(test, trn, k) for k in MARKOV_ORDERS]
        out["markov_note"] = (
            "in_sample is the empirical entropy with every parameter free: not a code, shown to say how "
            "far fitting a sequence on itself cheats; adaptive pays for its parameters as it learns; "
            f"held_out is fitted on {train}"
        )
        out["markov"] = {
            "rows": table,
            "best_order_adaptive": min(table, key=lambda r: r["adaptive_bits_per_base"])["order"],
            "best_order_held_out": min(table, key=lambda r: r["held_out_bits_per_base"])["order"],
        }
        ledger.lap("markov table")
    if compressors:
        out["general_compressors"] = general_compressors(test.s)
        ledger.lap("general compressors")

    def keep(label: str, codes: np.ndarray) -> np.ndarray:
        if spill is None:
            return codes
        spill.mkdir(parents=True, exist_ok=True)
        mm = np.lib.format.open_memmap(
            spill / f"{chrom}_{label}.npy", mode="w+", dtype=codes.dtype, shape=codes.shape
        )
        mm[:] = codes
        mm.flush()
        return mm

    groups: dict[str, dict[str, np.ndarray]] = {
        "naive": {f"static{k}": keep(f"static{k}", static_model(test.s, trn.s, k)) for k in NAIVE_ORDERS}
    }
    ledger.lap("naive models")
    groups["adaptive"] = {
        f"adaptive{k}": keep(f"adaptive{k}", adaptive_model(test.s, k)) for k in ADAPTIVE_ORDERS
    }
    ledger.lap("adaptive models")

    layers: list[Layer] = []
    builders = [
        ("repeats", lambda: repeats_layer(test, trn, results_dir, ledger)),
        ("duplications", lambda: segdup_layer(test, results_dir, reference, ledger)),
        ("gc_cpg", lambda: gc_layer(test, trn, ledger)),
        ("coding", lambda: coding_layer(test, trn, ledger)),
        ("registry", lambda: ccre_layer(test, trn, results_dir, ledger)),
        ("motifs", lambda: motif_layer(test, results_dir, ledger)),
        ("tiers", lambda: tiers_layer(test, trn, results_dir)),
    ]
    for name, build in builders:
        layer = build()
        if layer is not None:
            layers.append(layer)
            groups[layer.name] = {k: keep(f"{layer.name}_{k}", v) for k, v in layer.models.items()}
            layer.models = groups[layer.name]
        ledger.lap(f"layer {name}")
    by = {x.name: x for x in layers}
    prime = by["duplications"].notes.pop("prime", None) if "duplications" in by else None
    if prime is not None:
        groups["adaptive_primed"] = {
            f"adaptive_primed{k}": keep(f"adaptive_primed{k}", adaptive_model(test.s, k, prime=prime))
            for k in ADAPTIVE_ORDERS
        }
        del prime
        ledger.lap("control: adaptive models primed with the partners")

    # the regions every stack is read over
    tiers = by["tiers"].notes
    rep = np.where(by["repeats"].notes["tandem"] > 0, 2, np.where(by["repeats"].claim, 1, 0)).astype(np.int64)
    dup = by["duplications"].notes["any"].astype(np.int64)
    cds = by["coding"].claim.astype(np.int64) if "coding" in by else np.zeros(n, dtype=np.int64)
    combo = np.ravel_multi_index(
        (tiers["tier_state"].astype(np.int64), cds, rep, dup, window_gc_bins(test)), COMBO_SHAPE
    )
    unique = (rep == 0) & (dup == 0)
    block = tiers["block_id"].astype(np.int64) * 2 + unique
    size_combo, size_block = int(np.prod(COMBO_SHAPE)), 2 * (len(tiers["blocks"]) + 1)
    bases_combo = np.bincount(combo, minlength=size_combo)
    bases_block = np.bincount(block, minlength=size_block)
    claims = {x.name: x.claim for x in layers}
    if "island_mask" in by["gc_cpg"].notes:
        claims["cpg_islands"] = by["gc_cpg"].notes["island_mask"]
    ledger.lap("regions")

    names = [x.name for x in layers]
    canonical = ["naive", "adaptive", "adaptive_primed", *names]

    def models_of(key: tuple[str, ...]) -> list[np.ndarray]:
        return [a for g in key for a in groups[g].values()]

    # the mixing window and temperature, chosen on the generic stack and charged
    grid = []
    generic = ("naive", "adaptive")
    for w in MIX_WINDOWS:
        for beta in MIX_BETAS:
            bits = mix(models_of(generic), w, beta)
            grid.append(
                {"window": w, "beta": beta, "bits_per_base": round(float(bits.sum(dtype=np.float64)) / n, 5)}
            )
    chosen = min(grid, key=lambda r: r["bits_per_base"])
    hyper_bits = math.log2(len(grid))
    ledger.lap("mixing grid")

    plan = plan_stacks(names, control="adaptive_primed" in groups)
    aggs: dict[str, dict[str, Any]] = {}
    for keys in plan.values():
        for key in keys:
            key = tuple(sorted(key, key=canonical.index))
            name = stack_name(key)
            if name in aggs:
                continue
            bits = mix(models_of(key), chosen["window"], chosen["beta"])
            aggs[name] = {
                "key": key,
                "total": float(bits.sum(dtype=np.float64)) + hyper_bits,
                "combo": np.bincount(combo, weights=bits, minlength=size_combo),
                "block": np.bincount(block, weights=bits, minlength=size_block),
                "claims": {c: float(bits[m].sum(dtype=np.float64)) for c, m in claims.items()},
            }
            del bits
            if progress:
                progress(f"  mixed {name}")
    ledger.lap("mixtures")

    ann = {x.name: x.annotation_bits for x in layers}
    stacks = []
    for name, a in aggs.items():
        paid = sum(ann.get(g, 0.0) for g in a["key"])
        stacks.append(
            {
                "stack": name,
                "models": len(models_of(a["key"])),
                "bits_per_base": round(a["total"] / n, 4),
                "annotation_bits_per_base": round(paid / n, 5),
                "bits_per_base_with_annotation": round((a["total"] + paid) / n, 4),
            }
        )

    def gain(base: tuple[str, ...], with_: tuple[str, ...], layer: Layer) -> dict[str, Any]:
        a = aggs[stack_name(tuple(sorted(base, key=canonical.index)))]
        b = aggs[stack_name(tuple(sorted(with_, key=canonical.index)))]
        claimed = int(layer.claim.sum())
        whole = a["total"] - b["total"]
        inside = a["claims"][layer.name] - b["claims"][layer.name]
        net = whole - layer.annotation_bits
        return {
            "gross_bits_saved": round(whole),
            "saved_inside_claim_bits": round(inside),
            "saved_outside_claim_bits": round(whole - inside),
            "net_bits_saved": round(net),
            "gross_bits_per_claimed_base": round(inside / claimed, 4) if claimed else None,
            "net_bits_per_claimed_base": round(net / claimed, 4) if claimed else None,
            "net_bits_per_mb_of_chromosome": round(net / (n / 1e6)),
            "pays_for_itself": bool(net > 0),
        }

    full = ("naive", "adaptive", *names)
    layer_rows = []
    for i, x in enumerate(layers):
        claimed = int(x.claim.sum())
        row = {
            "layer": x.name,
            "title": x.title,
            "items": x.items,
            "claimed_bases": claimed,
            "claimed_share": round(claimed / n, 4),
            "annotation_bits": round(x.annotation_bits),
            "annotation_bits_per_claimed_base": round(x.annotation_bits / claimed, 4) if claimed else None,
            "over_naive": gain(("naive",), ("naive", x.name), x),
            "over_generic": gain(generic, (*generic, x.name), x),
            "over_repeat_aware": gain((*generic, "repeats"), (*generic, "repeats", x.name), x)
            if x.name != "repeats"
            else None,
            "cumulative_step": gain((*generic, *names[:i]), (*generic, *names[: i + 1]), x),
            "leave_one_out": gain(tuple(g for g in full if g != x.name), full, x),
            "over_generic_primed_with_partner_sequence": gain(
                ("naive", "adaptive_primed"), ("naive", "adaptive_primed", x.name), x
            )
            if x.name == "duplications" and "adaptive_primed" in groups
            else None,
            "claim_bits_per_base": {
                s: round(
                    aggs[stack_name(tuple(sorted(s_key, key=canonical.index)))]["claims"][x.name] / claimed, 4
                )
                for s, s_key in (
                    ("naive", ("naive",)),
                    ("generic", generic),
                    ("generic_plus_layer", (*generic, x.name)),
                    ("full", full),
                )
            }
            if claimed
            else None,
            "notes": {k: v for k, v in x.notes.items() if not isinstance(v, (np.ndarray, list))},
        }
        layer_rows.append(row)
    if "cpg_islands" in claims:
        m = int(claims["cpg_islands"].sum())
        g0 = aggs[stack_name(generic)]["claims"]["cpg_islands"]
        g1 = aggs[stack_name((*generic, "gc_cpg"))]["claims"]["cpg_islands"]
        l0 = aggs[stack_name(tuple(g for g in full if g != "gc_cpg"))]["claims"]["cpg_islands"]
        l1 = aggs[stack_name(full)]["claims"]["cpg_islands"]
        out["cpg_island_bases"] = {
            "bases": m,
            "generic_bits_per_base": round(g0 / m, 4),
            "generic_plus_gc_cpg_bits_per_base": round(g1 / m, 4),
            "full_without_gc_cpg_bits_per_base": round(l0 / m, 4),
            "full_bits_per_base": round(l1 / m, 4),
        }
    regions = summarise_regions(
        aggs, bases_combo, bases_block, tiers["blocks"], n, np.random.default_rng(SEED)
    )
    constrained = []
    for i, b in enumerate(tiers["blocks"]):
        if b["tier"] != "constrained_unknown":
            continue
        j = i + 1
        bp_all = int(bases_block[2 * j] + bases_block[2 * j + 1])
        if not bp_all:
            continue
        constrained.append(
            {
                "start": b["start"],
                "end": b["end"],
                "class": b["class"],
                "bases": bp_all,
                "unique_share": round(float(bases_block[2 * j + 1] / bp_all), 3),
                **{
                    f"{lab}_bits_per_base": round(
                        float(
                            (aggs[stack_name(k)]["block"][2 * j] + aggs[stack_name(k)]["block"][2 * j + 1])
                            / bp_all
                        ),
                        3,
                    )
                    for lab, k in (("naive", ("naive",)), ("generic", generic), ("full", full))
                },
            }
        )
    cost = ledger.summary()
    mb = n / 1e6
    cost |= {
        "megabases_coded": round(mb, 2),
        "seconds_per_megabase": round(cost["seconds"] / mb, 2),
        "models": sum(len(g) for g in groups.values()),
        "mixtures": len(aggs) + len(grid),
        "model_memory_gb": round(sum(len(g) for g in groups.values()) * 2 * n / 1e9, 2),
    }
    out |= {
        "method": {
            "fitting": f"static models fitted on both strands of {train}; adaptive models fitted as "
            "they code; annotation costs charged; mixture window and temperature chosen from a grid "
            "and charged",
            "naive_orders": NAIVE_ORDERS,
            "adaptive_orders": ADAPTIVE_ORDERS,
            "pseudocounts": {"order_up_to_8": ALPHA_LOW, "above": ALPHA_HIGH},
            "mixing_grid": grid,
            "mixing_chosen": chosen,
            "mixing_choice_bits": round(hyper_bits, 2),
            "mixing_floor": MIX_FLOOR,
            "duplication_partners": "an earlier chromosome in karyotype order, or earlier on the same one",
            "gc_window": GC_WINDOW,
        },
        "stacks": sorted(stacks, key=lambda r: r["bits_per_base"], reverse=True),
        "layers": layer_rows,
        "regions": regions,
        "constrained_unknown_blocks": constrained,
        "cost": cost,
    }
    return out


def run_and_save(
    chrom: str, train: str, results_dir: Path = RESULTS_DIR, name: str | None = None, **kw
) -> Path:
    """compress_<chrom>, or the given name for a sensitivity run."""
    return save_result(name or f"compress_{chrom}", run(chrom, train, results_dir, **kw), results_dir)
