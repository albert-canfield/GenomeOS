# SPDX-License-Identifier: AGPL-3.0-or-later
"""The lexicon: one index of the genome's units at every scale, with the nulls that make counts real.

Albert's framing of 2026-09-13: treat the genome as a maths problem, take the pieces of the code
from small to big, index them, and use the index as the base for comparison, statistics and the
attempt at syntax. The levels, small to big:

1. k-mers, counted per context rather than globally, against an order-2 Markov expectation fitted
   inside the same context (so GC and dinucleotide composition are already accounted for);
2. recurring segments found from the sequence itself: exact seeds of SEED_K bases that occur at
   least SEED_MIN times, extended by consensus over their occurrences. The data-driven vocabulary,
   not a curated one;
3. curated units the project already holds: RepeatMasker families and subfamilies, ENCODE cCRE
   classes, GENCODE coding and non-coding exons, segmental duplications, JASPAR sites where the
   motif runs have recorded them;
4. UNKNOWN blocks with their class and budget tier;
5. CTCF nodes;
6. libraries (`genomeos/lib`) and the chromosome.

Every unit carries its occurrences, the context they sit in (coding, exonic, regulatory, intronic,
or an unknown tier), the composition around them (GC, repeat share), its conservation (the
100-vertebrate conserved elements per base) and the human axis of the blocks it falls in (gnomAD
Gnocchi at a kilobase, from `variation_<chrom>`).

The part that decides whether any of it is real is the background model, and it is inside the
index rather than on top of it. Every count is read against a stratified expectation: the
chromosome is cut into windows, each window carries a GC bin and a repeat-share bin, and a unit's
expected occurrences in a context are its own density inside each stratum spread over that
context's bases in the same strata. The p is a Poisson tail on that expectation and the family of
tests is corrected by Benjamini-Hochberg, with the number of tests reported. The project has been
burned twice for want of this (the promoter operator pairs that vanished under a GC-by-repeat null,
and the VISTA panel where a held motif site said only that the sequence was conserved), so nothing
here is reported as enrichment without its null beside it.

Nothing is streamed: the index is built from the sequence and the committed results already on
disk. Per-base conservation is the local 100-vertebrate element cache and the human axis is the
block-level Gnocchi of `variation_<chrom>`; the per-base phyloP and per-kilobase Gnocchi upgrades
would cost range requests, and the cost is stated in the summary rather than spent here.
"""

from __future__ import annotations

import bisect
import gzip
import json
import math
import random
import time
from array import array
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from genomeos.results import RESULTS_DIR, load_result, save_result

KNOWLEDGE = Path("data/knowledge/lexicon")
REFERENCE = Path("data/reference")
KMER_K = 6  # k-mers counted per context; 4,096 units, an order-2 Markov null inside each context
SEED_K = 16  # the data-driven vocabulary: exact seeds of this width
SEED_MIN = 50  # a seed must occur at least this often on the chromosome to be a unit
SEED_TOP = 2_000  # the most frequent seeds kept as units (the tail is a long one)
SEED_OCCURRENCE_CAP = 400  # occurrences kept per seed, sampled beyond this
SEED_STORE = 50_000  # positions stored per seed before sampling; a satellite beyond it keeps its first
SEED_EXTEND = 120  # bases of consensus extension either side of a seed
SEED_CONSENSUS = 0.7  # share of occurrences that must agree on a base for the extension to go on
WINDOW = 10_000  # the stratification window
GC_BINS = 5
REPEAT_BINS = 4
PERMUTATIONS = 1_000
FDR = 0.05
MIN_OCCURRENCES = 20  # a unit is tested for context enrichment only above this many occurrences
CONSERVED_MIN = 0.10  # share of a unit's bases inside a 100-vertebrate element: held across species
HUMAN_MIN = 0.25  # share of its blocks' kilobases in Gnocchi's top decile (variation.py's bar)
SEED = 20260913

CONTEXTS = (
    "cds",
    "exon_noncoding",
    "utr_or_exon",
    "promoter_like",
    "enhancer_like",
    "ctcf_only",
    "intron",
    "unknown_structural",
    "unknown_fossil",
    "unknown_regulatory",
    "unknown_constrained",
    "unknown_neutral",
    "other",
)
CONTEXT_CODE = {c: i for i, c in enumerate(CONTEXTS)}
# a unit that paints a context is not tested for enrichment in it: the answer is the definition
SELF_DEFINED = {
    "ccre:PLS": ("promoter_like",),
    "ccre:DNase-H3K4me3": ("promoter_like",),
    "ccre:pELS": ("enhancer_like",),
    "ccre:dELS": ("enhancer_like",),
    "ccre:CTCF-only": ("ctcf_only",),
    "gencode:cds": ("cds",),
    "gencode:exon": ("cds", "utr_or_exon", "exon_noncoding"),
}
TIER_CONTEXT = {
    "structural": "unknown_structural",
    "fossil": "unknown_fossil",
    "regulatory": "unknown_regulatory",
    "constrained_unknown": "unknown_constrained",
    "neutral": "unknown_neutral",
}
CCRE_CONTEXT = {
    "PLS": "promoter_like",
    "DNase-H3K4me3": "promoter_like",
    "pELS": "enhancer_like",
    "dELS": "enhancer_like",
    "CTCF-only": "ctcf_only",
}
EVIDENCE = {
    "sequence": "curated: GRCh38 primary assembly",
    "repeats": "curated: UCSC RepeatMasker (rmsk_<chrom>)",
    "annotation": "curated: GENCODE v50 (CDS, exons)",
    "registry": "curated: ENCODE SCREEN cCREs v3",
    "blocks": "inferred: genomeos unknown classes with the budget's tiers",
    "nodes": "inferred: CTCF-only elements as node boundaries",
    "conservation": "experimental-derived: UCSC phastConsElements100way, per base (local cache)",
    "human_axis": "inferred: gnomAD Gnocchi per block at a kilobase (variation_<chrom>)",
    "null": (
        "stratified: a unit's own density inside GC x repeat-share strata spread over the context's "
        "bases; Poisson tail; Benjamini-Hochberg over the family of tests, whose size is reported"
    ),
}


# ---- statistics --------------------------------------------------------------------------


def poisson_tail(observed: int, expected: float) -> float:
    """P(X >= observed) for X ~ Poisson(expected); 1.0 when nothing is expected or seen."""
    if observed <= 0:
        return 1.0
    if expected <= 0:
        return 0.0

    def log_term(k: int) -> float:
        return -expected + k * math.log(expected) - math.lgamma(k + 1)

    # anchored at the observed count in logs, so a large expectation cannot underflow exp(-expected)
    if observed > expected:
        # upper tail summed outward from the observed count; the terms fall geometrically
        term = math.exp(log_term(observed))
        total, k = 0.0, observed
        while term > 0 and (term > 1e-15 * total or total == 0):
            total += term
            k += 1
            term *= expected / k
        return max(0.0, min(1.0, total))
    # lower tail summed downward from observed - 1, then subtracted
    k = observed - 1
    term = math.exp(log_term(k))
    total = 0.0
    while k >= 0 and term > 1e-17 * max(total, 1e-300):
        total += term
        term *= k / expected
        k -= 1
    return max(0.0, min(1.0, 1.0 - total))


def benjamini_hochberg(pvalues: list[float], fdr: float = FDR) -> tuple[float, int]:
    """The largest p that passes Benjamini-Hochberg at `fdr`, and how many pass."""
    if not pvalues:
        return 0.0, 0
    ordered = sorted(pvalues)
    n = len(ordered)
    threshold, passing = 0.0, 0
    for i, p in enumerate(ordered, start=1):
        if p <= fdr * i / n:
            threshold, passing = p, i
    return threshold, passing


def cosine(a: dict[str, float], b: dict[str, float]) -> float:
    keys = set(a) & set(b)
    if not keys:
        return 0.0
    num = sum(a[k] * b[k] for k in keys)
    na = math.sqrt(sum(v * v for v in a.values()))
    nb = math.sqrt(sum(v * v for v in b.values()))
    return num / (na * nb) if na and nb else 0.0


# ---- the per-base maps -------------------------------------------------------------------


def _bed(path: Path):
    if not path.exists():
        return
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            yield line.rstrip("\n").split("\t")


class Maps:
    """Per-base arrays for one chromosome: context, repeat flag, conserved flag, human-axis flag.

    A byte per base and three bits elsewhere: chromosome 21 costs about 100 MB of memory and
    nothing on disk. Contexts are assigned by priority, coding first and "other" last.
    """

    def __init__(self, chrom: str, length: int, results_dir: Path = RESULTS_DIR) -> None:
        self.chrom = chrom
        self.length = length
        self.context = array("B", bytes([CONTEXT_CODE["other"]])) * length
        self.repeat = bytearray(length)
        self.conserved = bytearray(length)
        self.human = bytearray(length)
        self.human_measured = bytearray(length)
        self.cds: list[tuple[int, int]] = []
        self.exons: list[tuple[int, int]] = []
        self.blocks: list[dict[str, Any]] = []
        self._fill(results_dir)

    def _paint(self, start: int, end: int, ctx: str) -> None:
        a, b = max(0, start), min(self.length, end)
        if b > a:
            self.context[a:b] = array("B", bytes([CONTEXT_CODE[ctx]])) * (b - a)

    @staticmethod
    def _flag(arr: bytearray, start: int, end: int) -> None:
        a, b = max(0, start), min(len(arr), end)
        if b > a:
            arr[a:b] = b"\x01" * (b - a)

    def _fill(self, results_dir: Path) -> None:
        from genomeos.genome import Annotation
        from genomeos.genome.annotation import default_gencode

        # painted in rising priority: blocks, gene bodies, the registry, exons, coding sequence
        budget = {
            (b["start"], b["end"]): (b["guess"] or {}).get("tier")
            for b in (load_result(f"budget_{self.chrom}", results_dir) or {}).get("blocks", [])
        }
        var = {
            (b["start"], b["end"]): b
            for b in (load_result(f"variation_{self.chrom}", results_dir) or {}).get("blocks", [])
        }
        for b in (load_result(f"unknown_{self.chrom}", results_dir) or {}).get("blocks", []):
            key = (b["start"], b["end"])
            tier = budget.get(key)
            v = var.get(key) or {}
            human = ((v.get("gnocchi") or {}).get("fraction_above")) or 0.0
            self.blocks.append(
                {
                    "start": b["start"],
                    "end": b["end"],
                    "length": b["length"],
                    "class": b["class"],
                    "tier": tier,
                    "case": (v.get("case") or {}).get("case"),
                    "human_fraction": human,
                }
            )
            if tier:
                self._paint(b["start"], b["end"], TIER_CONTEXT[tier])
            if (v.get("gnocchi") or {}).get("bases"):
                self._flag(self.human_measured, b["start"], b["end"])
            if human >= HUMAN_MIN:
                self._flag(self.human, b["start"], b["end"])
        self.blocks.sort(key=lambda b: b["start"])
        ann = Annotation.from_gff3(default_gencode({self.chrom}), {self.chrom})
        genes = [g for g in ann.genes.values() if g.locus.chrom == self.chrom]
        self.coding_tss = {
            g.symbol: (g.locus.end - 1 if g.locus.strand.value == "-" else g.locus.start)
            for g in genes
            if g.type == "protein_coding"
        }
        for g in genes:
            self._paint(g.locus.start, g.locus.end, "intron")
        for f in _bed(results_dir / f"ccres_{self.chrom}.bed.gz"):
            ctx = CCRE_CONTEXT.get(f[4])
            if ctx:
                self._paint(int(f[1]), int(f[2]), ctx)
        for g in genes:
            coding = g.type == "protein_coding"
            for t in g.transcripts.values():
                for e in t.exons:
                    self._paint(e.start, e.end, "utr_or_exon" if coding else "exon_noncoding")
                    self.exons.append((e.start, e.end))
        for g in genes:
            for t in g.transcripts.values():
                for c, _ in t.cds:
                    self._paint(c.start, c.end, "cds")
                    self.cds.append((c.start, c.end))
        for f in _bed(results_dir / f"rmsk_{self.chrom}.bed.gz"):
            if f[2] not in ("Simple_repeat", "Low_complexity"):
                self._flag(self.repeat, int(f[0]), int(f[1]))
        for f in _bed(Path("data/knowledge/constraint") / f"phastConsElements100way_{self.chrom}.bed.gz"):
            self._flag(self.conserved, int(f[1]), int(f[2]))

    # summaries over an interval ------------------------------------------------------------
    def context_at(self, pos: int) -> str:
        return CONTEXTS[self.context[pos]]

    def share(self, kind: str, start: int, end: int) -> float:
        arr = getattr(self, kind)
        a, b = max(0, start), min(self.length, end)
        return arr[a:b].count(1) / max(1, b - a)


# ---- windows and strata ------------------------------------------------------------------


class Windows:
    """The stratification: WINDOW-base windows with a GC bin, a repeat bin and their context mix."""

    def __init__(self, seq: str, maps: Maps, window: int = WINDOW) -> None:
        import numpy as np

        self.window = window
        self.rows: list[dict[str, Any]] = []
        context = np.frombuffer(maps.context, dtype=np.uint8)
        conserved = np.frombuffer(maps.conserved, dtype=np.uint8)
        for start in range(0, maps.length, window):
            end = min(maps.length, start + window)
            s = seq[start:end]
            acgt = sum(s.count(b) for b in "ACGT")
            if acgt < 0.5 * (end - start):
                continue  # an assembly gap: no window
            gc = (s.count("G") + s.count("C")) / acgt
            rep = maps.share("repeat", start, end)
            chunk = maps.context[start:end].tobytes()
            mix = {k: chunk.count(bytes([k])) for k in range(len(CONTEXTS))}
            mix = {k: v for k, v in mix.items() if v}
            held = np.bincount(context[start:end][conserved[start:end] == 1], minlength=len(CONTEXTS))
            self.rows.append(
                {
                    "start": start,
                    "end": end,
                    "bases": end - start,
                    "gc": gc,
                    "repeat": rep,
                    "context_bp": {CONTEXTS[k]: v for k, v in mix.items()},
                    "context_conserved": {CONTEXTS[k]: int(v) for k, v in enumerate(held) if v},
                    "conserved": maps.conserved[start:end].count(1),
                    "human": maps.human[start:end].count(1),
                    "human_measured": maps.human_measured[start:end].count(1),
                }
            )
        self._bin()

    def _bin(self) -> None:
        gcs = sorted(r["gc"] for r in self.rows)
        reps = sorted(r["repeat"] for r in self.rows)

        def edges(xs: list[float], n: int) -> list[float]:
            return [xs[int(k * (len(xs) - 1) / n)] for k in range(1, n)] if xs else []

        self.gc_edges = edges(gcs, GC_BINS)
        self.repeat_edges = edges(reps, REPEAT_BINS)
        for r in self.rows:
            r["stratum"] = (
                sum(1 for e in self.gc_edges if r["gc"] >= e),
                sum(1 for e in self.repeat_edges if r["repeat"] >= e),
            )
        self.starts = [r["start"] for r in self.rows]
        # bases per (stratum, context), the denominator of every expectation
        self.bp: dict[tuple, dict[str, int]] = defaultdict(lambda: defaultdict(int))
        for r in self.rows:
            for ctx, n in r["context_bp"].items():
                self.bp[r["stratum"]][ctx] += n
        self.stratum_bp = {k: sum(v.values()) for k, v in self.bp.items()}
        self.conserved_bp: dict[tuple, int] = defaultdict(int)
        self.human_bp: dict[tuple, int] = defaultdict(int)
        for r in self.rows:
            self.conserved_bp[r["stratum"]] += r["conserved"]
            self.human_bp[r["stratum"]] += r["human"]
        self.conserved_share = {k: self.conserved_bp[k] / v for k, v in self.stratum_bp.items() if v}
        held: dict[tuple, dict[str, int]] = defaultdict(lambda: defaultdict(int))
        for r in self.rows:
            for ctx, n in r["context_conserved"].items():
                held[r["stratum"]][ctx] += n
        # the conserved share of one context inside one stratum: the null for a unit's conservation
        self.conserved_share_in_context = {
            (st, ctx): held[st].get(ctx, 0) / n for st, per in self.bp.items() for ctx, n in per.items() if n
        }
        measured: dict[tuple, int] = defaultdict(int)
        for r in self.rows:
            measured[r["stratum"]] += r["human_measured"]
        # the human axis is a share of the bases where it was measured, never of all bases
        self.human_share = {k: self.human_bp[k] / v for k, v in measured.items() if v}

    def stratum_at(self, pos: int) -> tuple | None:
        i = bisect.bisect_right(self.starts, pos) - 1
        if 0 <= i < len(self.rows) and self.rows[i]["start"] <= pos < self.rows[i]["end"]:
            return self.rows[i]["stratum"]
        return None


def context_enrichment(
    occurrences: list[tuple[int, str, tuple | None]], windows: Windows
) -> dict[str, dict[str, Any]]:
    """One unit's occurrences against the stratified expectation, context by context.

    `occurrences` are (position, context, stratum). Inside a stratum the unit's own density is
    spread over that stratum's bases; a context's expectation is the sum over strata of that
    density times the context's bases in the stratum. No context effect, no enrichment.
    """
    per_stratum: Counter = Counter()
    for _, _, st in occurrences:
        if st is not None:
            per_stratum[st] += 1
    density = {st: n / windows.stratum_bp[st] for st, n in per_stratum.items() if windows.stratum_bp.get(st)}
    observed: Counter = Counter(ctx for _, ctx, st in occurrences if st is not None)
    out: dict[str, dict[str, Any]] = {}
    for ctx in CONTEXTS:
        exp = sum(d * windows.bp[st].get(ctx, 0) for st, d in density.items())
        obs = observed.get(ctx, 0)
        if obs == 0 and exp < 1:
            continue
        out[ctx] = {
            "observed": obs,
            "expected": round(exp, 2),
            "ratio": round(obs / exp, 3) if exp else None,
            "p": poisson_tail(obs, exp),
        }
    return out


# ---- level 1: k-mers per context ----------------------------------------------------------


def kmer_table(seq: str, maps: Maps, k: int = KMER_K) -> dict[str, Any]:
    """Every k-mer counted inside every context, with an order-2 Markov expectation fitted in the
    same context: the null keeps that context's GC and dinucleotide composition. The context of a
    k-mer is the context of its middle base."""
    counts: dict[int, Counter] = defaultdict(Counter)
    ctx = maps.context
    half = k // 2
    for i in range(len(seq) - k + 1):
        counts[ctx[i + half]][seq[i : i + k]] += 1
    out: dict[str, Any] = {}
    for code, table in counts.items():
        c = CONTEXTS[code]
        table = Counter({w: n for w, n in table.items() if "N" not in w})
        total = sum(table.values())
        if total < 10_000:
            continue
        di: Counter = Counter()
        tri: Counter = Counter()
        for w, n in table.items():
            for x in range(k - 1):
                di[w[x : x + 2]] += n
            for x in range(k - 2):
                tri[w[x : x + 3]] += n
        rows = []
        for w, obs in table.items():
            exp = _markov2_expectation(w, di, tri, total)
            rows.append({"kmer": w, "observed": obs, "expected": round(exp, 2), "p": poisson_tail(obs, exp)})
        thr, passing = benjamini_hochberg([r["p"] for r in rows])
        twofold = sum(1 for r in rows if r["p"] <= thr and r["observed"] >= 2 * r["expected"])
        rows.sort(key=lambda r: -(r["observed"] / max(0.5, r["expected"])))
        out[c] = {
            "kmers_seen": len(rows),
            "positions": total,
            "tests": len(rows),
            "bh_threshold": thr,
            "passing_bh": passing,
            "passing_bh_twofold_over_markov": twofold,
            "top_enriched": [
                {**r, "ratio": round(r["observed"] / max(0.5, r["expected"]), 2)} for r in rows[:15]
            ],
        }
    return out


def _markov2_expectation(w: str, di: Counter, tri: Counter, total: int) -> float:
    """The expected count of a word under an order-2 Markov chain fitted on the same context:
    P(w) = P(w0 w1) * prod P(w_i | w_{i-2} w_{i-1})."""
    tri_total = sum(tri.values()) or 1
    di_total = sum(di.values()) or 1
    p = (di.get(w[:2], 0) + 0.5) / di_total
    for i in range(2, len(w)):
        num = (tri.get(w[i - 2 : i + 1], 0) + 0.5) / tri_total
        den = (di.get(w[i - 2 : i], 0) + 0.5) / di_total
        p *= num / den
    return max(1e-9, p) * total


# ---- level 2: recurring segments from the sequence itself -----------------------------------


def seeds(seq: str, k: int = SEED_K, minimum: int = SEED_MIN, top: int = SEED_TOP) -> dict[str, list[int]]:
    """The recurring vocabulary: exact k-mers occurring at least `minimum` times, the `top` most
    frequent kept with up to SEED_OCCURRENCE_CAP positions each (the true count under `count:`).
    Two passes keep memory bounded: a hashed counter wide enough that collisions cannot make a
    bucket hot, then the exact strings only at the positions of hot buckets."""
    buckets = 1 << 26
    counts = bytearray(buckets)
    n = len(seq)
    for i in range(0, n - k + 1):
        w = seq[i : i + k]
        if "N" in w:
            continue
        h = hash(w) & (buckets - 1)
        if counts[h] < 255:
            counts[h] += 1
    floor = min(minimum, 255)
    total: Counter = Counter()
    positions: dict[str, array] = defaultdict(lambda: array("I"))
    for i in range(0, n - k + 1):
        w = seq[i : i + k]
        if counts[hash(w) & (buckets - 1)] < floor or "N" in w:
            continue
        total[w] += 1
        if total[w] <= SEED_STORE:
            positions[w].append(i)
    del counts
    kept = [w for w, c in total.most_common(top) if c >= minimum]
    positions = {w: sample_loci(positions[w]) for w in kept}
    return {w: positions[w] + [-total[w]] for w in kept}


def is_tandem(word: str, longest: int = 6) -> bool:
    """A word that is one short period repeated (a microsatellite), which RepeatMasker's
    interspersed classes do not cover and the repeat flag leaves out."""
    return any(word == (word[:p] * (len(word) // p + 1))[: len(word)] for p in range(1, longest + 1))


def _locus_key(pos: int) -> int:
    """A fixed pseudo-random key per 256-base bin (a splitmix step), shared by every seed."""
    z = ((pos >> 8) + SEED + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
    return z ^ (z >> 31)


def sample_loci(positions, cap: int | None = None) -> list[int]:
    """Occurrences kept to the cap by the smallest locus keys: uniform along the chromosome and
    consistent across seeds, so two seeds of one longer segment keep the same loci (a sample drawn
    per seed would break the co-occurrence the families are found by)."""
    cap = SEED_OCCURRENCE_CAP if cap is None else cap
    if len(positions) <= cap:
        return list(positions)
    return sorted(sorted(positions, key=_locus_key)[:cap])


def extend(seq: str, positions: list[int], k: int = SEED_K, reach: int = SEED_EXTEND) -> dict[str, Any]:
    """Consensus extension of a seed's occurrences: how far the flanks agree, and on what."""
    sample = positions[:SEED_OCCURRENCE_CAP]
    left = right = 0
    for step in range(1, reach + 1):
        bases = [seq[p - step] for p in sample if p - step >= 0]
        if not bases:
            break
        base, count = Counter(bases).most_common(1)[0]
        if base == "N" or count / len(bases) < SEED_CONSENSUS:
            break
        left = step
    for step in range(reach):
        bases = [seq[p + k + step] for p in sample if p + k + step < len(seq)]
        if not bases:
            break
        base, count = Counter(bases).most_common(1)[0]
        if base == "N" or count / len(bases) < SEED_CONSENSUS:
            break
        right = step + 1
    return {"left": left, "right": right, "consensus_length": k + left + right}


def seed_families(found: dict[str, list[int]], k: int = SEED_K, share: float = 0.3) -> list[dict[str, Any]]:
    """Seeds that are pieces of one longer recurring segment: two seeds join when at least `share` of
    the rarer one's occurrences sit at one fixed offset (under SEED_EXTEND) from the other's. The
    families are the vocabulary at segment scale; their span is the offsets they cover plus k."""
    placed = sorted((p, w) for w, pos in found.items() for p in pos)
    pair: Counter = Counter()
    for i, (p, w) in enumerate(placed):
        j = i + 1
        while j < len(placed) and placed[j][0] - p < SEED_EXTEND:
            q, v = placed[j]
            if v != w and q > p:
                pair[(w, v, q - p)] += 1
            j += 1
    parent = {w: w for w in found}

    def root(w: str) -> str:
        while parent[w] != w:
            parent[w] = parent[parent[w]]
            w = parent[w]
        return w

    offsets: dict[str, int] = {}
    for (w, v, d), n in pair.items():
        if n >= share * min(len(found[w]), len(found[v])):
            parent[root(v)] = root(w)
            offsets[v] = max(offsets.get(v, 0), d)
    groups: dict[str, list[str]] = defaultdict(list)
    for w in found:
        groups[root(w)].append(w)
    out = []
    for members in groups.values():
        members.sort(key=lambda w: -len(found[w]))
        out.append(
            {
                "seeds": members,
                "occurrences_of_leading_seed": len(found[members[0]]),
                "span": k + max((offsets.get(w, 0) for w in members), default=0),
            }
        )
    out.sort(key=lambda f: -f["occurrences_of_leading_seed"])
    return out


# ---- the index ---------------------------------------------------------------------------


def unit_rows(
    name: str,
    kind: str,
    positions: list[tuple[int, int]],
    maps: Maps,
    windows: Windows,
    block_at,
    candidates: list[tuple[int, int]] | None = None,
) -> dict[str, Any]:
    """One unit: its occurrences, their contexts, conservation, human axis and composition."""
    occ = [(a, maps.context_at((a + b) // 2), windows.stratum_at((a + b) // 2)) for a, b in positions]
    bp = sum(b - a for a, b in positions)
    cons = sum(maps.conserved[a:b].count(1) for a, b in positions)
    human = sum(maps.human[a:b].count(1) for a, b in positions)
    measured = sum(maps.human_measured[a:b].count(1) for a, b in positions)
    rep = sum(maps.repeat[a:b].count(1) for a, b in positions)
    tiers: Counter = Counter()
    cases: Counter = Counter()
    for a, b in positions:
        blk = block_at((a + b) // 2)
        if blk:
            tiers[blk["tier"] or "untiered"] += 1
            cases[blk["case"] or "unmeasured"] += 1
    bp_by_stratum: Counter = Counter()
    bp_by_cell: Counter = Counter()
    measured_by_stratum: Counter = Counter()
    for (a, b), (_, ctx, st) in zip(positions, occ, strict=True):
        if st is not None:
            bp_by_stratum[st] += b - a
            bp_by_cell[(st, ctx)] += b - a
            measured_by_stratum[st] += maps.human_measured[a:b].count(1)
    exp_cons_strata = sum(n * windows.conserved_share.get(st, 0.0) for st, n in bp_by_stratum.items())
    exp_cons = sum(n * windows.conserved_share_in_context.get(cell, 0.0) for cell, n in bp_by_cell.items())
    exp_human = sum(n * windows.human_share.get(st, 0.0) for st, n in measured_by_stratum.items())
    # bases inside one occurrence are not independent draws (a conserved element covers a run of
    # them), so the tail is taken on occurrence-equivalents: bases over the mean occurrence length
    mean_len = bp / len(positions) if positions else 1.0
    return {
        "unit": name,
        "kind": kind,
        "occurrences": len(positions),
        "bp": bp,
        "contexts": dict(Counter(c for _, c, _ in occ)),
        "conserved_fraction": round(cons / bp, 4) if bp else None,
        "human_measured_fraction": round(measured / bp, 4) if bp else None,
        "human_constrained_fraction": round(human / measured, 4) if measured else None,
        "repeat_fraction": round(rep / bp, 4) if bp else None,
        "conserved_expected_fraction_strata": (
            round(exp_cons_strata / bp, 4) if bp and bp_by_stratum else None
        ),
        "conserved_expected_fraction": (round(exp_cons / bp, 4) if bp and bp_by_cell else None),
        "conserved_p": (poisson_tail(round(cons / mean_len), exp_cons / mean_len) if bp_by_cell else None),
        "human_expected_fraction": (round(exp_human / measured, 4) if measured else None),
        "human_p": poisson_tail(round(human / mean_len), exp_human / mean_len) if measured else None,
        "in_syntax_candidates": (
            sum(1 for a, b in positions for c0, c1 in candidates if a < c1 and b > c0) if candidates else 0
        ),
        "block_tiers": dict(tiers),
        "block_cases": dict(cases),
        "enrichment": (
            {c: e for c, e in context_enrichment(occ, windows).items() if c not in SELF_DEFINED.get(name, ())}
            if len(positions) >= MIN_OCCURRENCES
            else {}
        ),
    }


def curated_units(
    chrom: str, maps: Maps, results_dir: Path = RESULTS_DIR
) -> dict[str, list[tuple[int, int]]]:
    """Level 3: the units the project already curates, as occurrence lists per unit name."""
    out: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for f in _bed(results_dir / f"rmsk_{chrom}.bed.gz"):
        cls, fam, sub = f[2], f[3], f[4]
        if cls in ("Simple_repeat", "Low_complexity"):
            continue
        out[f"repeat:{fam}:{sub}"].append((int(f[0]), int(f[1])))
        out[f"repeat_family:{fam}"].append((int(f[0]), int(f[1])))
    for f in _bed(results_dir / f"ccres_{chrom}.bed.gz"):
        out[f"ccre:{f[4]}"].append((int(f[1]), int(f[2])))
    out["gencode:cds"] = sorted(set(maps.cds))
    out["gencode:exon"] = sorted(set(maps.exons))
    for f in _bed(results_dir / f"superdups_{chrom}.bed.gz"):
        out["segmental_duplication"].append((int(f[0]), int(f[1])))
    for name, (sites, _) in jaspar_sites(chrom, results_dir).items():
        out[name] = sites
    return out


def jaspar_sites(
    chrom: str, results_dir: Path = RESULTS_DIR
) -> dict[str, tuple[list[tuple[int, int]], list[tuple[int, int]]]]:
    """JASPAR sites as the motif runs recorded them: per factor, the 8-base sites and, aligned, the
    element each sits in. The elements are the constrained and enhancer targets, so a site's
    conservation has to be read against its own element, not against the chromosome."""
    out: dict[str, tuple[list, list]] = defaultdict(lambda: ([], []))
    motifs = load_result(f"motifs_{chrom}", results_dir) or {}
    for el in motifs.get("elements", []):
        for r in el.get("requires", []):
            sites, hosts = out[f"jaspar:{r['factor']}"]
            a = el["start"] + r["position"]
            sites.append((a, a + 8))
            hosts.append((el["start"], el["end"]))
    return dict(out)


def within_host(
    unit: dict[str, Any], sites: list[tuple[int, int]], hosts: list[tuple[int, int]], maps: Maps
) -> None:
    """The site's conserved bases against the conserved share of the rest of its own element."""
    cons = sum(maps.conserved[a:b].count(1) for a, b in sites)
    expected = 0.0
    for (a, b), (h0, h1) in zip(sites, hosts, strict=True):
        rest = (h1 - h0) - (b - a)
        held = maps.conserved[h0:h1].count(1) - maps.conserved[a:b].count(1)
        expected += (b - a) * (held / rest if rest > 0 else 0.0)
    bp = sum(b - a for a, b in sites)
    unit["host_conserved_expected_fraction"] = round(expected / bp, 4) if bp else None
    unit["host_conserved_p"] = poisson_tail(round(cons / 8), expected / 8) if bp else None


def build(chrom: str, results_dir: Path = RESULTS_DIR, progress=None, kmers: bool = True):
    """The index for one chromosome: every level, every unit with its statistics."""
    from genomeos.coords import Locus
    from genomeos.genome.index import IndexedGenome

    t0 = time.time()

    def say(m: str) -> None:
        if progress:
            progress(f"{time.time() - t0:6.0f} s  {m}")

    genome = IndexedGenome(REFERENCE / f"{chrom}.fa")
    length = genome.lengths[chrom]
    seq = str(genome.fetch(Locus(chrom, 0, length))).upper()
    genome.close()
    say(f"{chrom}: {length / 1e6:.1f} Mb read")
    maps = Maps(chrom, length, results_dir)
    say("per-base maps built")
    windows = Windows(seq, maps)
    say(f"{len(windows.rows)} windows, {len(windows.stratum_bp)} strata")
    block_starts = [b["start"] for b in maps.blocks]

    def block_at(pos: int) -> dict[str, Any] | None:
        i = bisect.bisect_right(block_starts, pos) - 1
        if 0 <= i < len(maps.blocks) and maps.blocks[i]["start"] <= pos < maps.blocks[i]["end"]:
            return maps.blocks[i]
        return None

    candidates = [
        (c["start"], c["end"])
        for c in (load_result("syntax_candidates_genome_wide", results_dir) or {}).get("candidates", [])
        if c["chrom"] == chrom
    ]
    units: list[dict[str, Any]] = []
    curated = curated_units(chrom, maps, results_dir)
    hosts = jaspar_sites(chrom, results_dir)
    for name, positions in curated.items():
        kind = name.split(":")[0]
        row = unit_rows(name, kind, positions, maps, windows, block_at, candidates)
        if name in hosts:
            within_host(row, hosts[name][0], hosts[name][1], maps)
        units.append(row)
    say(f"{len(curated)} curated units")
    found = seeds(seq)
    for w, pos in found.items():
        count = -pos[-1]
        pos = pos[:-1]
        row = unit_rows(
            f"seed:{w}",
            "seed",
            [(p, p + SEED_K) for p in pos[:SEED_OCCURRENCE_CAP]],
            maps,
            windows,
            block_at,
            candidates,
        )
        row["extension"] = extend(seq, pos)
        row["occurrences_total"] = count
        units.append(row)
    families = seed_families({w: p[:-1] for w, p in found.items()})
    family_of = {w: i for i, fam in enumerate(families) for w in fam["seeds"]}
    for u in units:
        if u["kind"] == "seed":
            u["seed_family"] = family_of.get(u["unit"][5:])
    for fam in families:
        members = [
            u for u in units if u["kind"] == "seed" and u.get("seed_family") == family_of[fam["seeds"][0]]
        ]
        fam["repeat_fraction"] = round(sum(u["repeat_fraction"] or 0 for u in members) / len(members), 3)
        fam["conserved_fraction"] = round(
            sum(u["conserved_fraction"] or 0 for u in members) / len(members), 4
        )
        fam["dominant_context"] = Counter(
            c for u in members for c, n in u["contexts"].items() for _ in range(n)
        ).most_common(1)[0][0]
    seed_rows = [u for u in units if u["kind"] == "seed"]
    tandem = [u for u in seed_rows if (u["repeat_fraction"] or 0) < 0.5 and is_tandem(u["unit"][5:])]
    unexplained = [u for u in seed_rows if (u["repeat_fraction"] or 0) < 0.5 and not is_tandem(u["unit"][5:])]
    seed_makeup = {
        "mostly_repeatmasker": sum(1 for u in seed_rows if (u["repeat_fraction"] or 0) >= 0.5),
        "tandem_period_up_to_6": len(tandem),
        "neither": len(unexplained),
        "neither_examples": [
            {k: u.get(k) for k in ("unit", "occurrences_total", "repeat_fraction", "conserved_fraction")}
            for u in sorted(unexplained, key=lambda u: -u["occurrences_total"])[:10]
        ],
    }
    say(f"{len(found)} recurring seeds in {len(families)} families")
    kmer = kmer_table(seq, maps) if kmers else {}
    say("k-mers done" if kmers else "k-mers skipped")
    # levels 4 to 6: blocks, nodes, libraries, chromosome
    nodes = (load_result(f"domains_{chrom}", results_dir) or {}).get("domains", [])
    libraries = library_units(chrom, results_dir)
    tests = sum(len(u["enrichment"]) for u in units)
    all_p = [e["p"] for u in units for e in u["enrichment"].values()]
    thr, passing = benjamini_hochberg(all_p)
    for u in units:
        for ctx, e in u["enrichment"].items():
            e["passes_bh"] = bool(e["p"] <= thr and e["observed"] > e["expected"])
            del ctx
    index = {
        "chrom": chrom,
        "length": length,
        "levels": {
            "kmers": {"k": KMER_K, "per_context": kmer},
            "seeds": {
                "k": SEED_K,
                "minimum": SEED_MIN,
                "kept": len(found),
                "families": len(families),
                "largest_families": sorted(families, key=lambda f: -len(f["seeds"]))[:15],
                "makeup": seed_makeup,
            },
            "curated": len(curated),
            "blocks": len(maps.blocks),
            "nodes": len(nodes),
            "libraries": len(libraries),
        },
        "units": units,
        "windows": {
            "window": WINDOW,
            "count": len(windows.rows),
            "gc_edges": [round(x, 4) for x in windows.gc_edges],
            "repeat_edges": [round(x, 4) for x in windows.repeat_edges],
            "strata": len(windows.stratum_bp),
            "bases_per_context": {c: sum(v.get(c, 0) for v in windows.bp.values()) for c in CONTEXTS},
        },
        "library_units": libraries,
        "tests": {"context_enrichment": tests, "bh_threshold": thr, "passing": passing, "fdr": FDR},
        "evidence": EVIDENCE,
        "seconds": round(time.time() - t0, 1),
    }
    return index, {"seq": seq, "maps": maps, "windows": windows, "block_at": block_at, "nodes": nodes}


def library_units(chrom: str, results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """Level 6: the libraries whose genes sit on this chromosome, with their genes' nodes."""
    from genomeos.genome import Annotation
    from genomeos.genome.annotation import default_gencode
    from genomeos.lib.membership import KnowledgeBase

    distilled = KnowledgeBase.distilled() or {}
    members = distilled.get("members") or {}
    ann = Annotation.from_gff3(default_gencode({chrom}), {chrom})
    symbols = {g.symbol for g in ann.protein_coding() if g.locus.chrom == chrom}
    out = {}
    for lib, genes in members.items():
        here = sorted(symbols & set(genes))
        if len(here) >= 5:
            out[lib] = here
    return out


# ---- question 1: is the fossil tier a library for its node? -------------------------------


def fossil_composition(chrom: str, maps: Maps, results_dir: Path = RESULTS_DIR) -> list[dict[str, Any]]:
    """Per fossil-tier block: the repeat subfamily bases inside it, its GC-free length and its node."""
    blocks = [b for b in maps.blocks if b["tier"] == "fossil"]
    starts = [b["start"] for b in blocks]
    comp: dict[int, Counter] = {i: Counter() for i in range(len(blocks))}
    for f in _bed(results_dir / f"rmsk_{chrom}.bed.gz"):
        if f[2] in ("Simple_repeat", "Low_complexity"):
            continue
        a, b = int(f[0]), int(f[1])
        i = bisect.bisect_right(starts, a) - 1
        if 0 <= i < len(blocks) and blocks[i]["end"] > a:
            comp[i][f"{f[3]}:{f[4]}"] += min(blocks[i]["end"], b) - max(blocks[i]["start"], a)
    out = []
    for i, b in enumerate(blocks):
        if sum(comp[i].values()) < 500:
            continue
        out.append({**b, "composition": dict(comp[i]), "repeat_bp": sum(comp[i].values())})
    return out


def assign_nodes(items: list[dict[str, Any]], nodes: list[dict[str, Any]]) -> None:
    starts = [n["start"] for n in nodes]
    for it in items:
        mid = (it["start"] + it["end"]) // 2
        i = bisect.bisect_right(starts, mid) - 1
        it["node"] = nodes[i]["id"] if 0 <= i < len(nodes) and nodes[i]["end"] > mid else None


def fossil_node_similarity(
    blocks: list[dict[str, Any]], nodes: list[dict[str, Any]], seed: int = SEED
) -> dict[str, Any]:
    """Are the fossils of one node more alike than fossils of other nodes?

    Two nulls. The first permutes the node label over blocks inside strata of GC-free length and
    repeat share, so like is compared with like. The second shifts every node boundary by one
    random offset, which keeps the nodes' sizes and their contiguity but moves the boundaries: it
    asks whether the CTCF partition carries anything beyond the fact that neighbours are alike.
    """
    assign_nodes(blocks, nodes)
    sel = [b for b in blocks if b["node"]]
    if len(sel) < 10:
        return {"blocks": len(sel), "note": "too few fossil blocks with a node to test"}
    lengths = sorted(b["repeat_bp"] for b in sel)
    edges = [lengths[int(k * (len(lengths) - 1) / 4)] for k in range(1, 4)]
    for b in sel:
        b["_stratum"] = sum(1 for e in edges if b["repeat_bp"] >= e)
    pairs = [(i, j) for i in range(len(sel)) for j in range(i + 1, len(sel))]
    cos = {(i, j): cosine(sel[i]["composition"], sel[j]["composition"]) for i, j in pairs}

    def within(labels: list[str]) -> float:
        vals = [cos[(i, j)] for i, j in pairs if labels[i] == labels[j]]
        return sum(vals) / len(vals) if vals else 0.0

    labels = [b["node"] for b in sel]
    observed = within(labels)
    between = sum(cos[(i, j)] for i, j in pairs if labels[i] != labels[j]) / max(
        1, sum(1 for i, j in pairs if labels[i] != labels[j])
    )
    rng = random.Random(seed)
    by_stratum: dict[int, list[int]] = defaultdict(list)
    for i, b in enumerate(sel):
        by_stratum[b["_stratum"]].append(i)
    null = []
    for _ in range(PERMUTATIONS):
        shuffled = list(labels)
        for idx in by_stratum.values():
            picked = [labels[i] for i in idx]
            rng.shuffle(picked)
            for i, lab in zip(idx, picked, strict=True):
                shuffled[i] = lab
        null.append(within(shuffled))
    p_perm = (sum(1 for x in null if x >= observed) + 1) / (len(null) + 1)
    # the boundary-shift null
    span = max(n["end"] for n in nodes)
    shift_null = []
    node_starts = [n["start"] for n in nodes]
    for _ in range(PERMUTATIONS // 5):
        off = rng.randrange(span)
        moved = []
        for b in sel:
            mid = ((b["start"] + b["end"]) // 2 + off) % span
            i = bisect.bisect_right(node_starts, mid) - 1
            moved.append(nodes[i]["id"] if 0 <= i < len(nodes) else "none")
        shift_null.append(within(moved))
    p_shift = (sum(1 for x in shift_null if x >= observed) + 1) / (len(shift_null) + 1)
    return {
        "blocks": len(sel),
        "nodes_with_two_or_more": sum(1 for _, n in Counter(labels).items() if n >= 2),
        "within_node_mean_cosine": round(observed, 4),
        "between_node_mean_cosine": round(between, 4),
        "null_label_permutation": round(sum(null) / len(null), 4),
        "p_label_permutation": round(p_perm, 4),
        "null_boundary_shift": round(sum(shift_null) / len(shift_null), 4),
        "p_boundary_shift": round(p_shift, 4),
        "permutations": len(null),
        "shifts": len(shift_null),
    }


def fossil_library_link(
    chrom: str,
    maps: Maps,
    windows: Windows,
    nodes: list[dict[str, Any]],
    libraries: dict[str, list[str]],
    results_dir: Path = RESULTS_DIR,
    seed: int = SEED,
) -> dict[str, Any]:
    """Does a node's fossil composition relate to the libraries its genes belong to?

    For every library and every repeat subfamily: the subfamily's occurrences inside the fossil
    tier of the nodes that hold the library's genes, against the stratified expectation (the
    subfamily's own density inside each GC x repeat stratum spread over those nodes' fossil bases
    in the same stratum). Benjamini-Hochberg over all pairs tested, and the same test on random
    gene sets of the same sizes as the calibration, because occurrences inside a node are not
    independent and a Poisson tail will over-call.
    """
    node_starts = [n["start"] for n in nodes]

    def node_at(pos: int) -> str | None:
        i = bisect.bisect_right(node_starts, pos) - 1
        return nodes[i]["id"] if 0 <= i < len(nodes) and nodes[i]["end"] > pos else None

    # fossil bases per (node, stratum) and the subfamily occurrences inside them
    fossil_bp: dict[tuple[str, tuple], int] = defaultdict(int)
    for r in windows.rows:
        n = node_at((r["start"] + r["end"]) // 2)
        if n:
            fossil_bp[(n, r["stratum"])] += r["context_bp"].get("unknown_fossil", 0)
    occ: dict[str, list[tuple[str, tuple]]] = defaultdict(list)
    for f in _bed(results_dir / f"rmsk_{chrom}.bed.gz"):
        if f[2] in ("Simple_repeat", "Low_complexity"):
            continue
        mid = (int(f[0]) + int(f[1])) // 2
        if maps.context_at(mid) != "unknown_fossil":
            continue
        n, st = node_at(mid), windows.stratum_at(mid)
        if n and st:
            occ[f"{f[3]}:{f[4]}"].append((n, st))
    totals: dict[tuple, int] = defaultdict(int)
    for (_, st), bp in fossil_bp.items():
        totals[st] += bp
    families = {k: v for k, v in occ.items() if len(v) >= MIN_OCCURRENCES}
    genes_of = {lib: g for lib, g in libraries.items()}
    all_genes = sorted({g for gs in libraries.values() for g in gs})
    rng = random.Random(seed)

    def test(gene_set: list[str]) -> list[dict[str, Any]]:
        their_nodes = {node_at(maps.coding_tss[g]) for g in gene_set if g in maps.coding_tss}
        their_nodes.discard(None)
        target_bp: dict[tuple, int] = defaultdict(int)
        for (n, st), bp in fossil_bp.items():
            if n in their_nodes:
                target_bp[st] += bp
        rows = []
        for fam, places in families.items():
            density = {st: 0.0 for st in totals}
            counts: Counter = Counter(st for _, st in places)
            for st in totals:
                if totals[st]:
                    density[st] = counts.get(st, 0) / totals[st]
            expected = sum(density[st] * target_bp.get(st, 0) for st in totals)
            observed = sum(1 for n, _ in places if n in their_nodes)
            if expected <= 0 and observed == 0:
                continue
            rows.append(
                {
                    "family": fam,
                    "observed": observed,
                    "expected": round(expected, 2),
                    "ratio": round(observed / expected, 3) if expected else None,
                    "p": poisson_tail(observed, expected),
                }
            )
        return rows

    results = {}
    every_p = []
    for lib, gs in genes_of.items():
        rows = test(gs)
        results[lib] = rows
        every_p.extend(r["p"] for r in rows)
    thr, passing = benjamini_hochberg(every_p)
    top = []
    for lib, rows in results.items():
        for r in rows:
            if r["p"] <= thr and (r["ratio"] or 0) > 1:
                top.append({"library": lib, **r})
    top.sort(key=lambda r: r["p"])
    # the calibration: random gene sets of the same sizes
    random_passing = []
    for _ in range(10):
        ps = []
        for gs in genes_of.values():
            draw = rng.sample(all_genes, min(len(gs), len(all_genes)))
            ps.extend(r["p"] for r in test(draw))
        random_passing.append(benjamini_hochberg(ps)[1])
    return {
        "libraries_tested": len(genes_of),
        "families_tested": len(families),
        "tests": len(every_p),
        "bh_threshold": thr,
        "passing_bh": passing,
        "passing_with_ratio_above_one": len(top),
        "top": top[:15],
        "random_gene_sets_passing": random_passing,
        "reading": (
            "the random gene sets say how many pairs a Poisson tail passes when the only structure is "
            "that repeats cluster inside nodes; a positive claim needs the real libraries to pass more"
        ),
    }


# ---- question 2: syntax against value slots ------------------------------------------------


def syntax_classes(index: dict[str, Any], results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """Units invariant across mammals and among people are syntax candidates; units everywhere but
    variable are value slots. Each unit is also read against its own neighbourhood: the conserved
    and human-constrained shares expected from the strata its occurrences sit in."""
    rows = [u for u in index["units"] if u["occurrences"] >= MIN_OCCURRENCES and u["bp"]]
    cons_p = [u["conserved_p"] for u in rows if u["conserved_p"] is not None]
    thr_cons, pass_cons = benjamini_hochberg(cons_p)
    human_p = [u["human_p"] for u in rows if u["human_p"] is not None]
    thr_human, pass_human = benjamini_hochberg(human_p)
    host_p = [u["host_conserved_p"] for u in rows if u.get("host_conserved_p") is not None]
    thr_host, pass_host = benjamini_hochberg(host_p)
    classes: Counter = Counter()
    per_kind: dict[str, Counter] = defaultdict(Counter)
    syntax_units = []
    slot_units = []
    survival: Counter = Counter()
    for u in rows:
        held = (u["conserved_fraction"] or 0) >= CONSERVED_MIN
        if u["human_constrained_fraction"] is None:
            u["axis_case"] = "unmeasured"
            classes["unmeasured"] += 1
            per_kind[u["kind"]]["unmeasured"] += 1
            continue
        human = u["human_constrained_fraction"] >= HUMAN_MIN
        case = "syntax" if held and human else "relaxed" if held else "recent" if human else "tolerant"
        u["axis_case"] = case
        u["conserved_above_neighbourhood"] = bool(
            u["conserved_p"] is not None
            and u["conserved_p"] <= thr_cons
            and (u["conserved_fraction"] or 0) > (u["conserved_expected_fraction"] or 0)
        )
        if u.get("host_conserved_p") is not None:
            u["conserved_above_own_element"] = bool(
                u["host_conserved_p"] <= thr_host
                and (u["conserved_fraction"] or 0) > (u["host_conserved_expected_fraction"] or 0)
            )
        classes[case] += 1
        per_kind[u["kind"]][case] += 1
        wide = len(u["contexts"]) >= 4 and u["occurrences"] >= 100
        if case == "syntax":
            survival["syntax_by_thresholds"] += 1
            if u["conserved_above_neighbourhood"]:
                survival["and_above_context_matched_neighbourhood"] += 1
                if u.get("conserved_above_own_element") is False:
                    survival["jaspar_lost_to_own_element_null"] += 1
                else:
                    survival["surviving_every_null"] += 1
                    syntax_units.append(u)
        if case in ("relaxed", "tolerant") and wide and u["conserved_above_neighbourhood"]:
            slot_units.append(u)
    jaspar = [u for u in index["units"] if u.get("host_conserved_expected_fraction") is not None]
    j_obs = sum((u["conserved_fraction"] or 0) * u["bp"] for u in jaspar)
    j_exp = sum(u["host_conserved_expected_fraction"] * u["bp"] for u in jaspar)
    jaspar_pooled = {
        "factors": len(jaspar),
        "site_bases": sum(u["bp"] for u in jaspar),
        "conserved_fraction": round(j_obs / max(1, sum(u["bp"] for u in jaspar)), 4),
        "own_element_expected_fraction": round(j_exp / max(1, sum(u["bp"] for u in jaspar)), 4),
        "p_site_equivalents": poisson_tail(round(j_obs / 8), j_exp / 8) if jaspar else None,
        "note": "sites of different factors overlap, so the pooled p is optimistic",
    }
    by_tier: dict[str, Counter] = defaultdict(Counter)
    by_block_case: dict[str, Counter] = defaultdict(Counter)
    for u in rows:
        if u["block_tiers"]:
            by_tier[max(u["block_tiers"], key=u["block_tiers"].get)][u["axis_case"]] += 1
        if u["block_cases"]:
            by_block_case[max(u["block_cases"], key=u["block_cases"].get)][u["axis_case"]] += 1
    in_candidates = [u for u in rows if u.get("in_syntax_candidates")]
    syntax_units.sort(key=lambda u: -(u["conserved_fraction"] or 0))
    slot_units.sort(key=lambda u: -u["occurrences"])
    return {
        "units_tested": len(rows),
        "by_case": dict(classes),
        "by_kind": {k: dict(v) for k, v in per_kind.items()},
        "conservation_tests": len(cons_p),
        "conservation_bh_threshold": thr_cons,
        "conservation_passing": pass_cons,
        "human_tests": len(human_p),
        "human_bh_threshold": thr_human,
        "human_passing": pass_human,
        "own_element_tests": len(host_p),
        "own_element_bh_threshold": thr_host,
        "own_element_passing": pass_host,
        "jaspar_pooled_against_own_element": jaspar_pooled,
        "jaspar_above_own_element": sum(1 for u in rows if u.get("conserved_above_own_element")),
        "syntax_survival": dict(survival),
        "value_slots_found": len(slot_units),
        "syntax_candidates": [
            {
                k: u[k]
                for k in (
                    "unit",
                    "kind",
                    "occurrences",
                    "conserved_fraction",
                    "conserved_expected_fraction",
                    "human_constrained_fraction",
                    "block_tiers",
                    "contexts",
                )
            }
            for u in syntax_units[:20]
        ],
        "value_slot_candidates": [
            {
                k: u[k]
                for k in (
                    "unit",
                    "kind",
                    "occurrences",
                    "conserved_fraction",
                    "conserved_expected_fraction",
                    "human_constrained_fraction",
                    "contexts",
                )
            }
            for u in slot_units[:20]
        ],
        "by_dominant_tier": {k: dict(v) for k, v in by_tier.items()},
        "by_dominant_block_case": {k: dict(v) for k, v in by_block_case.items()},
        "units_in_the_syntax_candidates": [
            {
                k: u.get(k)
                for k in (
                    "unit",
                    "kind",
                    "in_syntax_candidates",
                    "occurrences",
                    "axis_case",
                    "conserved_fraction",
                )
            }
            for u in sorted(in_candidates, key=lambda u: -u["in_syntax_candidates"])[:20]
        ],
        "thresholds": {"conserved_min": CONSERVED_MIN, "human_min": HUMAN_MIN, "fdr": FDR},
    }


# ---- the run ------------------------------------------------------------------------------


def summary(index: dict[str, Any], questions: dict[str, Any], index_path: Path) -> dict[str, Any]:
    """What is committed: the counts, the tests, the two questions and the cost, never the units."""
    units = index["units"]
    by_kind: Counter = Counter(u["kind"] for u in units)
    enriched = [
        {"unit": u["unit"], "context": ctx, **e}
        for u in units
        for ctx, e in u["enrichment"].items()
        if e.get("passes_bh") and e["observed"] >= 10
    ]
    enriched.sort(key=lambda r: -(r["ratio"] or 0))
    return {
        "chrom": index["chrom"],
        "length": index["length"],
        "levels": {
            k: (v if not isinstance(v, dict) else {kk: vv for kk, vv in v.items() if kk != "per_context"})
            for k, v in index["levels"].items()
        }
        | {
            "seeds": {
                **index["levels"]["seeds"],
                "largest_families": [
                    {**f, "seeds": f["seeds"][:8], "seed_count": len(f["seeds"])}
                    for f in index["levels"]["seeds"]["largest_families"]
                ],
            }
        },
        "kmers_per_context": {
            c: {k: v for k, v in d.items() if k != "top_enriched"}
            for c, d in index["levels"]["kmers"]["per_context"].items()
        },
        "kmer_examples": {
            c: d["top_enriched"][:5] for c, d in index["levels"]["kmers"]["per_context"].items()
        },
        "units": len(units),
        "units_by_kind": dict(by_kind),
        "windows": index["windows"],
        "tests": index["tests"],
        "context_enrichments_passing": len(enriched),
        "most_enriched": enriched[:25],
        "fossil_as_node_library": questions["fossil"],
        "fossil_library_link": questions["libraries"],
        "syntax_against_slots": questions["syntax"],
        "index_file": str(index_path),
        "index_bytes": index_path.stat().st_size if index_path.exists() else None,
        "evidence": index["evidence"],
        "seconds": index["seconds"],
    }


def run_and_save(
    chrom: str, results_dir: Path = RESULTS_DIR, progress=None, kmers: bool = True
) -> dict[str, Any]:
    t0 = time.time()
    index, ctx = build(chrom, results_dir, progress, kmers)
    maps, windows, nodes = ctx["maps"], ctx["windows"], ctx["nodes"]
    blocks = fossil_composition(chrom, maps, results_dir)
    fossil = fossil_node_similarity(blocks, nodes)
    libraries = fossil_library_link(chrom, maps, windows, nodes, index["library_units"], results_dir)
    syn = syntax_classes(index, results_dir)
    questions = {"fossil": fossil, "libraries": libraries, "syntax": syn}
    KNOWLEDGE.mkdir(parents=True, exist_ok=True)
    path = KNOWLEDGE / f"lexicon_{chrom}.json.gz"
    with gzip.open(path, "wt") as fh:
        json.dump(index, fh)
    out = summary(index, questions, path)
    out["seconds_with_questions"] = round(time.time() - t0, 1)
    out["cost"] = {
        "seconds_per_mb": round(out["seconds_with_questions"] / (index["length"] / 1e6), 2),
        "index_bytes_per_mb": (
            round(out["index_bytes"] / (index["length"] / 1e6)) if out["index_bytes"] else None
        ),
    }
    save_result(f"lexicon_{chrom}", out, results_dir)
    return out


def main() -> None:  # pragma: no cover - the job entry point
    import sys

    chroms = sys.argv[1:] or ["chr21", "chr22"]
    for c in chroms:
        out = run_and_save(c, progress=lambda m: print("  " + m, flush=True))
        print(
            json.dumps(
                {
                    k: out[k]
                    for k in (
                        "chrom",
                        "units",
                        "units_by_kind",
                        "tests",
                        "context_enrichments_passing",
                        "fossil_as_node_library",
                        "fossil_library_link",
                        "syntax_against_slots",
                        "index_bytes",
                        "seconds",
                    )
                },
                indent=1,
            )[:6000]
        )
