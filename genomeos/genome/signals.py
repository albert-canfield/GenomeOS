"""Sequence signals learned from data: the real "tags" of the genome.

A gene is not delimited by exact tokens the way HTML is. It is delimited by
statistical signals that the cell's machinery recognises with some
probability: a start codon in Kozak context, GT...AG splice sites with their
neighbourhoods, a stop codon, a polyadenylation signal, a CpG-rich promoter.
This module learns those signals from an annotation as position weight
matrices (PWMs) and scores raw sequence with them, so that GenomeOS finds
structure the way gene finders do, with a log-odds score and a confidence,
rather than by fitting DNA into categories we invented.
"""

from __future__ import annotations

import json
import math
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from genomeos.coords import Locus, Strand
from genomeos.genome.sequence import Sequence

BASES = "ACGT"


@dataclass(slots=True)
class Pwm:
    """Position weight matrix with a background; scores are log2 odds."""

    name: str
    offset: int  # position of the anchor (e.g. the G of GT) inside the window
    length: int
    counts: list[dict[str, int]]
    background: dict[str, float] = field(default_factory=lambda: {b: 0.25 for b in BASES})
    examples: int = 0
    consensus: str = ""

    @classmethod
    def learn(cls, name: str, seqs: list[str], offset: int) -> Pwm:
        seqs = [s.upper() for s in seqs if s and all(c in BASES for c in s)]
        length = len(seqs[0])
        counts = [dict.fromkeys(BASES, 0) for _ in range(length)]
        for s in seqs:
            if len(s) != length:
                continue
            for i, c in enumerate(s):
                counts[i][c] += 1
        bg = Counter(c for s in seqs for c in s)
        total = sum(bg.values()) or 1
        pwm = cls(name, offset, length, counts, {b: bg[b] / total for b in BASES}, len(seqs))
        pwm.consensus = "".join(
            max(BASES, key=lambda b: c[b]) if max(c.values()) / max(1, sum(c.values())) > 0.5 else "n"
            for c in counts
        )
        return pwm

    def matrix(self) -> list[dict[str, float]]:
        out = []
        for c in self.counts:
            tot = sum(c.values()) + 4 * 0.5
            out.append({b: math.log2(((c[b] + 0.5) / tot) / self.background[b]) for b in BASES})
        return out

    def score(self, window: str) -> float:
        if len(window) != self.length:
            return -math.inf
        m = self.matrix()
        s = 0.0
        for i, ch in enumerate(window.upper()):
            if ch not in BASES:
                return -math.inf
            s += m[i][ch]
        return s

    def max_score(self) -> float:
        return sum(max(col.values()) for col in self.matrix())

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "offset": self.offset,
            "length": self.length,
            "counts": self.counts,
            "background": self.background,
            "examples": self.examples,
            "consensus": self.consensus,
        }

    @classmethod
    def from_dict(cls, d: dict) -> Pwm:
        return cls(
            d["name"], d["offset"], d["length"], d["counts"], d["background"], d["examples"], d["consensus"]
        )


@dataclass(slots=True)
class SignalSet:
    """The learned delimiters of one annotation: donor, acceptor, start (Kozak), stop, polyA."""

    pwms: dict[str, Pwm]
    stats: dict[str, object] = field(default_factory=dict)

    def save(self, path: str | Path) -> None:
        Path(path).write_text(
            json.dumps({"pwms": {k: v.to_dict() for k, v in self.pwms.items()}, "stats": self.stats})
        )

    @classmethod
    def load(cls, path: str | Path) -> SignalSet:
        d = json.loads(Path(path).read_text())
        return cls({k: Pwm.from_dict(v) for k, v in d["pwms"].items()}, d.get("stats", {}))


def learn_signals(annotation, genome, chrom: str) -> SignalSet:
    """Learn PWMs from the canonical protein-coding transcripts of one chromosome."""
    module = annotation.to_module("signals")

    def fetch(s: int, e: int, strand: Strand) -> str:
        if s < 0:
            return ""
        return str(genome.fetch(Locus(chrom, s, e, strand)))

    donors, acceptors, starts, stops, tails = [], [], [], [], []
    donor2 = Counter()
    acceptor2 = Counter()
    for gene in annotation.protein_coding():
        for tx in module.entities[gene.id].transcripts:
            if not tx.cds_segments or "Ensembl_canonical" not in tx.tags:
                continue
            strand = tx.exons[0].strand
            minus = strand is Strand.MINUS
            exons = sorted(tx.exons, key=lambda l: l.start, reverse=minus)
            for a, b in zip(exons, exons[1:], strict=False):
                i_s, i_e = (b.end, a.start) if minus else (a.end, b.start)
                if i_e - i_s < 20:
                    continue
                if minus:
                    donors.append(fetch(i_e - 6, i_e + 3, strand))
                    acceptors.append(fetch(i_s - 1, i_s + 14, strand))
                    donor2[fetch(i_e - 2, i_e, strand)] += 1
                    acceptor2[fetch(i_s, i_s + 2, strand)] += 1
                else:
                    donors.append(fetch(i_s - 3, i_s + 6, strand))
                    acceptors.append(fetch(i_e - 14, i_e + 1, strand))
                    donor2[fetch(i_s, i_s + 2, strand)] += 1
                    acceptor2[fetch(i_e - 2, i_e, strand)] += 1
            cds = sorted(tx.cds_segments, key=lambda l: l.start, reverse=minus)
            first, last = cds[0], cds[-1]
            if minus:
                starts.append(fetch(first.end - 4, first.end + 6, strand))
                stops.append(fetch(last.start, last.start + 3, strand))
                tails.append(fetch(exons[-1].start, exons[-1].start + 40, strand))
            else:
                starts.append(fetch(first.start - 6, first.start + 4, strand))
                stops.append(fetch(last.end - 3, last.end, strand))
                tails.append(fetch(exons[-1].end - 40, exons[-1].end, strand))
    n_intron = sum(donor2.values())
    pwms = {
        "splice_donor": Pwm.learn(
            "splice_donor", donors, offset=3
        ),  # window -3..+6, anchor = first intron base
        "splice_acceptor": Pwm.learn(
            "splice_acceptor", acceptors, offset=14
        ),  # window -14..+1, anchor = first exon base
        "start_kozak": Pwm.learn("start_kozak", starts, offset=6),  # window -6..+4, anchor = A of ATG
    }
    # how separable each signal is: relative scores of the true sites vs random windows
    import random

    rng = random.Random(0)
    chrom_len = genome.lengths[chrom] if hasattr(genome, "lengths") else genome.chromosomes[chrom].length
    separability: dict[str, dict[str, float]] = {}
    for name, examples in (("splice_donor", donors), ("splice_acceptor", acceptors), ("start_kozak", starts)):
        pwm = pwms[name]
        mx = pwm.max_score()
        rels = sorted(
            pwm.score(w) / mx for w in examples if len(w) == pwm.length and pwm.score(w) > -math.inf
        )
        if not rels:
            continue
        p10 = rels[len(rels) // 10]
        bg = []
        while len(bg) < 5000:
            st = rng.randrange(0, chrom_len - pwm.length)
            w = str(genome.fetch(Locus(chrom, st, st + pwm.length)))
            sc = pwm.score(w)
            if sc > -math.inf:
                bg.append(sc / mx)
        separability[name] = {
            "true_p10": round(p10, 3),
            "true_median": round(rels[len(rels) // 2], 3),
            "background_hits_per_kb_at_p10": round(sum(1 for x in bg if x >= p10) / len(bg) * 1000, 1),
        }
    stop_counts = Counter(s for s in stops if len(s) == 3)
    polya = sum(1 for t in tails if "AATAAA" in t or "ATTAAA" in t)
    stats = {
        "chromosome": chrom,
        "transcripts": len(starts),
        "introns": n_intron,
        "donor_dinucleotide": dict(donor2.most_common(4)),
        "acceptor_dinucleotide": dict(acceptor2.most_common(4)),
        "donor_GT_fraction": round(donor2["GT"] / n_intron, 4) if n_intron else None,
        "acceptor_AG_fraction": round(acceptor2["AG"] / n_intron, 4) if n_intron else None,
        "start_codons": dict(Counter(s[6:9] for s in starts if len(s) == 10).most_common(3)),
        "stop_codons": dict(stop_counts.most_common(4)),
        "polyA_signal_in_last_40nt": round(polya / max(1, len(tails)), 4),
        "separability": separability,
    }
    return SignalSet(pwms, stats)


@dataclass(slots=True)
class Hit:
    signal: str
    pos: int  # anchor position (0-based, forward coordinates)
    strand: Strand
    score: float
    relative: float  # score / max score, 0..1


def scan(seq: Sequence | str, signals: SignalSet, min_relative: float = 0.6) -> list[Hit]:
    """Score every position of a sequence with every PWM (both strands); return hits."""
    fwd = str(Sequence(str(seq)))
    rc = str(Sequence(fwd).reverse_complement())
    n = len(fwd)
    hits: list[Hit] = []
    for name, pwm in signals.pwms.items():
        mx = pwm.max_score()
        for strand, s in ((Strand.PLUS, fwd), (Strand.MINUS, rc)):
            for i in range(0, n - pwm.length + 1):
                sc = pwm.score(s[i : i + pwm.length])
                if sc == -math.inf or sc / mx < min_relative:
                    continue
                anchor = i + pwm.offset
                pos = anchor if strand is Strand.PLUS else n - anchor - 1
                hits.append(Hit(name, pos, strand, round(sc, 2), round(sc / mx, 3)))
    hits.sort(key=lambda h: (h.pos, h.signal))
    return hits
