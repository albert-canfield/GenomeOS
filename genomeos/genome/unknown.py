"""Investigate UNKNOWN blocks: classify the space between genes from the
sequence itself, largest block first.

Classes, each with the evidence that produced it:

    gap                    runs of N (assembly gap)                          curated
    telomere               TTAGGG arrays                                     curated
    centromere             CENP-B boxes / ~171 bp periodicity                curated/inferred
    tandem_repeat          simple tandem repeats dominate the block          predicted
    low_complexity         low k-mer entropy dominates the block             predicted
    interspersed_repeat    high-copy 16-mers (learned from the chromosome)   predicted
    long_orf               an open reading frame of 250-2,500 aa: unannotated
                           gene, pseudogene or transposon ORF (L1 ORF2 ~1,275 aa)   predicted
    satellite_array        a stop-free frame > 2,500 aa: a repeat array         inferred
    promoter_like          CpG island present                                predicted
    gene_desert            very large block with nothing above               inferred
    user pattern classes   from a patterns file (regex, class, evidence)     as declared
    unclassified           none of the above dominated                       none

The interspersed-repeat signal is learned, not assumed: every 16-mer of the
chromosome is counted; positions whose 16-mer occurs many times belong to a
repeated element. Named families come only from the pattern file (Alu core,
CENP-B box) so that names are never guessed.
"""

from __future__ import annotations

import math
import re
from array import array
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from genomeos.genome.sequence import Sequence

K = 16
BUCKETS = 1 << 24
HIGH_COPY = 40  # a 16-mer seen this often in one chromosome is part of a repeated element
INTERSPERSED_FRACTION = 0.5  # human intergenic space is about half interspersed repeats
ORF_MIN_AA, ORF_MAX_AA = 250, 2500  # a real coding remnant; longer frames without stops are repeats
DEFAULT_PATTERNS = Path("data/patterns/known_motifs.tsv")


@dataclass(slots=True)
class Pattern:
    name: str
    cls: str
    regex: re.Pattern
    evidence: str
    confidence: float
    note: str = ""


def load_patterns(path: str | Path = DEFAULT_PATTERNS) -> list[Pattern]:
    out: list[Pattern] = []
    p = Path(path)
    if not p.exists():
        return out
    for line in p.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 5:
            continue
        out.append(Pattern(f[0], f[1], re.compile(f[2]), f[3], float(f[4]), f[5] if len(f) > 5 else ""))
    return out


# ---------------------------------------------------------------- k-mers ----

_ENC = {"A": 0, "C": 1, "G": 2, "T": 3}


def kmer_table(seq: str, k: int = K) -> array:
    """Count every k-mer (hashed into 2^24 buckets) over a sequence. Pure Python, ~1 min for 50 Mb."""
    counts = array("H", [0]) * BUCKETS
    mask = (1 << (2 * k)) - 1
    code = 0
    valid = 0
    enc = _ENC
    for ch in seq:
        v = enc.get(ch)
        if v is None:
            valid = 0
            code = 0
            continue
        code = ((code << 2) | v) & mask
        valid += 1
        if valid >= k:
            b = (code * 0x9E3779B1) & (BUCKETS - 1)
            if counts[b] < 65535:
                counts[b] += 1
    return counts


def high_copy_fraction(seq: str, counts: array, k: int = K, threshold: int = HIGH_COPY) -> float:
    mask = (1 << (2 * k)) - 1
    code = valid = 0
    hi = n = 0
    for ch in seq:
        v = _ENC.get(ch)
        if v is None:
            valid = code = 0
            continue
        code = ((code << 2) | v) & mask
        valid += 1
        if valid >= k:
            n += 1
            if counts[(code * 0x9E3779B1) & (BUCKETS - 1)] >= threshold:
                hi += 1
    return hi / n if n else 0.0


ALU_CORE = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCA"  # 5' start of the Alu consensus (36 nt)
ALU_CORE_RC = str(Sequence(ALU_CORE).reverse_complement())
ALU_MAX_MISMATCH = 7  # ~20% divergence: old Alu copies still match; random sequence does not
ALU_PER_KB_HUMAN = 1 / 3  # about one Alu per 3 kb genome-wide


def fuzzy_count(seq: str, core: str, max_mismatch: int, step: int = 1) -> int:
    """Occurrences of `core` with at most `max_mismatch` substitutions (no indels)."""
    n, m = len(seq), len(core)
    hits = 0
    i = 0
    while i <= n - m:
        mm = 0
        for j in range(m):
            if seq[i + j] != core[j]:
                mm += 1
                if mm > max_mismatch:
                    break
        if mm <= max_mismatch:
            hits += 1
            i += m  # do not count the same copy twice
        else:
            i += step
    return hits


def alu_density(seq: str) -> float:
    """Alu-core matches (both strands, <= 7 mismatches) per kb."""
    if not seq:
        return 0.0
    return (
        (fuzzy_count(seq, ALU_CORE, ALU_MAX_MISMATCH) + fuzzy_count(seq, ALU_CORE_RC, ALU_MAX_MISMATCH))
        / len(seq)
        * 1000
    )


def kmer_set(seq: str, k: int = 20, step: int = 7) -> set[int]:
    return {hash(seq[i : i + k]) for i in range(0, max(0, len(seq) - k), step) if "N" not in seq[i : i + k]}


# ------------------------------------------------------------- features ----


def entropy_low_fraction(seq: str, window: int = 64) -> float:
    """Fraction of windows whose 3-mer entropy is low (< 2.5 bits): low-complexity sequence."""
    low = n = 0
    for i in range(0, len(seq) - window + 1, window):
        w = seq[i : i + window]
        c = Counter(w[j : j + 3] for j in range(window - 2))
        tot = sum(c.values())
        h = -sum(v / tot * math.log2(v / tot) for v in c.values())
        n += 1
        if h < 2.5:
            low += 1
    return low / n if n else 0.0


_TANDEM = re.compile(r"(?:A{10,}|C{10,}|G{10,}|T{10,})|(?:([ACGT]{2,6}))\1{4,}")


def tandem_fraction(seq: str) -> float:
    return sum(m.end() - m.start() for m in _TANDEM.finditer(seq)) / len(seq) if seq else 0.0


def periodicity_171(seq: str, max_len: int = 20000) -> float:
    """Autocorrelation of base identity at lag 171 (alpha-satellite monomer) minus background."""
    s = seq[:max_len]
    if len(s) < 600:
        return 0.0

    def match(lag: int) -> float:
        n = len(s) - lag
        return sum(1 for i in range(0, n, 3) if s[i] == s[i + lag]) / (n / 3)

    return match(171) - (match(100) + match(250)) / 2


def cpg_islands(seq: str, window: int = 200) -> int:
    n = 0
    run = 0
    for i in range(0, len(seq) - window + 1, window):
        w = seq[i : i + window]
        c, g = w.count("C"), w.count("G")
        if c and g and (c + g) / window > 0.5 and w.count("CG") * window / (c * g) > 0.6:
            run += 1
        else:
            if run >= 2:
                n += 1
            run = 0
    return n + (1 if run >= 2 else 0)


def longest_orf_aa(seq: str, min_aa: int = 100) -> int:
    from genomeos.runtime import find_orfs

    best = 0
    for o in find_orfs(seq, "u", min_aa=min_aa):
        best = max(best, o.length_aa)
    return best


# ------------------------------------------------------------ classify ----


@dataclass(slots=True)
class UnknownBlock:
    start: int
    end: int
    cls: str = "unclassified"
    evidence: str = "none"
    confidence: float = 0.0
    features: dict[str, Any] = field(default_factory=dict)
    patterns: dict[str, int] = field(default_factory=dict)
    similar_to: list[str] = field(default_factory=list)

    @property
    def length(self) -> int:
        return self.end - self.start

    def to_dict(self) -> dict[str, Any]:
        return {
            "start": self.start,
            "end": self.end,
            "length": self.length,
            "class": self.cls,
            "evidence": self.evidence,
            "confidence": self.confidence,
            "features": self.features,
            "patterns": self.patterns,
            "similar_to": self.similar_to,
        }


def composition(seq: str, counts: array | None, window: int = 5000) -> dict[str, float]:
    """Fraction of 5 kb windows dominated by each signal: the mix inside a large block."""
    kinds = {"interspersed": 0, "tandem": 0, "low_complexity": 0, "unique": 0}
    n = 0
    for i in range(0, max(1, len(seq) - window + 1), window):
        w = seq[i : i + window]
        if w.count("N") > window * 0.5:
            continue
        n += 1
        if tandem_fraction(w) > 0.5:
            kinds["tandem"] += 1
        elif entropy_low_fraction(w) > 0.5:
            kinds["low_complexity"] += 1
        elif (counts is not None and high_copy_fraction(w, counts) > 0.5) or alu_density(
            w
        ) >= ALU_PER_KB_HUMAN:
            kinds["interspersed"] += 1
        else:
            kinds["unique"] += 1
    return {k: round(v / n, 3) for k, v in kinds.items()} if n else {}


def classify_block(
    seq: str, counts: array | None, patterns: list[Pattern], ccre: dict[str, int] | None = None
) -> tuple[str, str, float, dict, dict]:
    n = len(seq)
    f: dict[str, Any] = {}
    if ccre:
        f["ccre"] = dict(ccre)
        f["ccre_per_10kb"] = round(sum(ccre.values()) / n * 10_000, 3)
    f["n_fraction"] = round(seq.count("N") / n, 3) if n else 0.0
    if f["n_fraction"] > 0.5:
        return "gap", "curated", 1.0, f, {}
    hits: dict[str, int] = {}
    cover: dict[str, float] = {}
    for p in patterns:
        ms = list(p.regex.finditer(seq))
        if ms:
            hits[p.name] = len(ms)
            cover[p.name] = sum(m.end() - m.start() for m in ms) / n
    acgt = seq.replace("N", "")
    f["gc"] = round((acgt.count("G") + acgt.count("C")) / max(1, len(acgt)), 3)
    f["tandem_fraction"] = round(tandem_fraction(seq), 3)
    f["low_complexity_fraction"] = round(entropy_low_fraction(seq), 3)
    f["high_copy_fraction"] = round(high_copy_fraction(seq, counts), 3) if counts is not None else None
    f["cpg_islands"] = cpg_islands(seq)
    f["longest_orf_aa"] = longest_orf_aa(seq, 100) if n <= 2_000_000 else longest_orf_aa(seq[:2_000_000], 100)
    f["periodicity_171"] = round(periodicity_171(seq), 3)
    f["alu_per_kb"] = round(alu_density(seq), 3)
    f["composition"] = composition(seq, counts) if n >= 20_000 else {}
    # user/declared patterns first: they carry their own evidence
    by_class: dict[str, tuple[int, str, float]] = {}
    for p in patterns:
        c = hits.get(p.name, 0)
        if not c:
            continue
        per_kb = c / n * 1000
        strong = (
            (p.cls == "telomere" and cover.get(p.name, 0) > 0.3)
            or (p.cls == "centromere" and per_kb > 1.0)
            or (p.cls.startswith("interspersed_repeat") and per_kb > 0.5)
        )
        if strong:
            by_class[p.cls] = (c, p.evidence, p.confidence)
    if "telomere" in by_class:
        return "telomere", by_class["telomere"][1], by_class["telomere"][2], f, hits
    # experimental chromatin evidence beats every sequence-only rule below
    if (
        ccre
        and f.get("ccre_per_10kb", 0) >= 2
        and (ccre.get("PLS", 0) + ccre.get("pELS", 0) + ccre.get("dELS", 0)) >= 2
    ):
        return "regulatory", "curated", 0.7, f, hits
    if "centromere" in by_class:
        return "centromere", by_class["centromere"][1], by_class["centromere"][2], f, hits
    if f["tandem_fraction"] > 0.5:
        return "tandem_repeat", "predicted", 0.7, f, hits
    if f["periodicity_171"] > 0.15 and f["tandem_fraction"] < 0.3:
        return "centromere", "inferred", 0.5, f, hits
    if f["low_complexity_fraction"] > 0.5:
        return "low_complexity", "predicted", 0.6, f, hits
    named = [k for k in by_class if k.startswith("interspersed_repeat")]
    if f["alu_per_kb"] >= 3 * ALU_PER_KB_HUMAN and n < 20_000:
        return "interspersed_repeat_SINE", "predicted", 0.6, f, hits
    if f["high_copy_fraction"] is not None and f["high_copy_fraction"] > INTERSPERSED_FRACTION:
        cls = named[0] if named else "interspersed_repeat"
        return cls, "predicted", 0.6 if named else 0.5, f, hits
    orf_ok = (
        ORF_MIN_AA <= f["longest_orf_aa"] <= ORF_MAX_AA
        and f["tandem_fraction"] < 0.2
        and (f["high_copy_fraction"] is None or f["high_copy_fraction"] < 0.5)
    )
    if f["longest_orf_aa"] > ORF_MAX_AA:
        # a stop-free frame tens of kb long is an array of a repeat unit with no stop in one frame
        return "satellite_array", "inferred", 0.4, f, hits
    if orf_ok:
        return "long_orf", "predicted", 0.4, f, hits
    islands_per_100kb = f["cpg_islands"] / n * 100_000
    if (n <= 20_000 and f["cpg_islands"] >= 1) or islands_per_100kb >= 2:
        return "promoter_like", "predicted", 0.5, f, hits
    comp = f["composition"]
    if comp:
        if comp["interspersed"] >= 0.5:
            return "interspersed_repeat", "predicted", 0.5, f, hits
        if comp["unique"] >= 0.7:
            # mostly unique, non-repetitive sequence between genes: where regulation is expected to live
            return "unique_intergenic", "inferred", 0.3, f, hits
        return "mixed_intergenic", "predicted", 0.4, f, hits
    if n >= 500_000:
        return "gene_desert", "inferred", 0.4, f, hits
    return "unclassified", "none", 0.0, f, hits


def unknown_blocks(annotation, chrom: str, length: int, min_size: int = 1000) -> list[UnknownBlock]:
    genes = sorted(
        (g for g in annotation.genes.values() if g.locus.chrom == chrom), key=lambda g: g.locus.start
    )
    out: list[UnknownBlock] = []
    cursor = 0
    for g in genes:
        if g.locus.start - cursor >= min_size:
            out.append(UnknownBlock(cursor, g.locus.start))
        cursor = max(cursor, g.locus.end)
    if length - cursor >= min_size:
        out.append(UnknownBlock(cursor, length))
    out.sort(key=lambda b: -b.length)
    return out


def investigate(
    seq: Sequence | str,
    annotation,
    chrom: str,
    patterns: list[Pattern] | None = None,
    min_size: int = 1000,
    progress=None,
    max_blocks: int | None = None,
) -> dict[str, Any]:
    from genomeos.genome.regulatory import ccre_index, count_in, load_ccres

    s = str(Sequence(str(seq)))
    patterns = patterns if patterns is not None else load_patterns()
    ccres = ccre_index(load_ccres(chrom))
    if progress:
        progress("counting 16-mers over the chromosome")
    counts = kmer_table(s)
    blocks = unknown_blocks(annotation, chrom, len(s), min_size)
    if max_blocks:
        blocks = blocks[:max_blocks]
    for i, b in enumerate(blocks):
        cc = count_in(ccres, b.start, b.end) if ccres else None
        cls, ev, conf, f, hits = classify_block(s[b.start : b.end], counts, patterns, cc)
        b.cls, b.evidence, b.confidence, b.features, b.patterns = cls, ev, conf, f, hits
        if progress and (i % 25 == 0 or i == len(blocks) - 1):
            progress(f"{i + 1}/{len(blocks)} blocks ({b.length:,} bp -> {cls})")
    # duplication candidates among the 60 largest blocks (shared 20-mers)
    big = [b for b in blocks[:60] if b.cls not in ("gap",)]
    sets = {id(b): kmer_set(s[b.start : b.end]) for b in big}
    for i, a in enumerate(big):
        for c in big[i + 1 :]:
            sa, sc = sets[id(a)], sets[id(c)]
            if not sa or not sc:
                continue
            j = len(sa & sc) / len(sa | sc)
            if j > 0.05:
                a.similar_to.append(f"{chrom}:{c.start}-{c.end} (jaccard {j:.2f})")
                c.similar_to.append(f"{chrom}:{a.start}-{a.end} (jaccard {j:.2f})")
    by_class: dict[str, dict[str, int]] = {}
    for b in blocks:
        d = by_class.setdefault(b.cls, {"blocks": 0, "bp": 0})
        d["blocks"] += 1
        d["bp"] += b.length
    total = sum(b.length for b in blocks)
    return {
        "chrom": chrom,
        "unknown_blocks": len(blocks),
        "unknown_bp": total,
        "by_class": dict(sorted(by_class.items(), key=lambda kv: -kv[1]["bp"])),
        "classified_fraction": round(1 - by_class.get("unclassified", {"bp": 0})["bp"] / total, 4)
        if total
        else None,
        "high_copy_threshold": HIGH_COPY,
        "patterns": [p.name for p in patterns],
        "ccres_used": len(ccres),
        "blocks": [b.to_dict() for b in blocks],
    }
