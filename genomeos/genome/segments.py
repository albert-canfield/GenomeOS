"""Segment parser v1: chain learned signals into candidate genes, never asserting.

The sequence grammar (docs/SEQUENCE-GRAMMAR.md) reads: start codon in a
Kozak context → coding exon → donor (GT) → intron → acceptor (AG) → coding
exon … → stop codon. `genomeos signals learn` gives position weight
matrices for the donor, the acceptor and the start; this module adds a
coding-potential table (codon log-odds, coding vs chromosome background)
learned from the same annotation, and runs a Viterbi-style dynamic
programme over candidate signal positions to find the best-scoring gene
structures on each strand.

The point is not to replace GENCODE. It is to measure how far sequence
signals alone carry: every prediction is scored against the annotation
(exon, splice-site and gene level sensitivity and precision) and saved as
`predicted` evidence with that measurement beside it.

Complexity is kept linear-ish with prefix sums (exon coding score in O(1)),
a next-in-frame-stop table, and a sliding-window maximum for introns (a
flat intron prior, no length model in v1).
"""

from __future__ import annotations

import bisect
import math
from collections import deque
from dataclasses import dataclass
from typing import Any

from genomeos.coords import Strand
from genomeos.genome.sequence import Sequence
from genomeos.genome.signals import SignalSet

STOPS = ("TAA", "TAG", "TGA")
MIN_EXON, MAX_EXON = 30, 5_000
MIN_INTRON, MAX_INTRON = 60, 50_000
START_PRIOR = -10.0  # log2 prior of a gene starting here (against random ATGs)
INTRON_PRIOR = -16.0  # log2 prior per intron: a splice pair must score well above chance to pay for itself
STOP_BONUS = 0.0
MIN_CDS = 300  # total coding bases a candidate must have (100 residues)
MIN_GENE_SCORE = 10.0  # bits over background a candidate must reach
EVIDENCE = "predicted: segment parser v1 (learned PWMs + codon log-odds, Viterbi over candidate signals)"
EVIDENCE_PREDICTED_SITES = (
    "predicted: segment parser v1 with AlphaGenome splice sites as the donor/acceptor candidates "
    "(learned Kozak matrix for starts, codon log-odds, Viterbi)"
)


@dataclass(slots=True)
class Prediction:
    chrom: str
    strand: Strand
    start: int  # first base of the start codon (forward coordinates, 0-based)
    end: int  # one past the last base of the stop codon
    exons: list[tuple[int, int]]  # CDS segments, forward coordinates, sorted
    score: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "chrom": self.chrom,
            "strand": self.strand.value,
            "start": self.start,
            "end": self.end,
            "exons": self.exons,
            "score": round(self.score, 2),
        }


# ---- coding potential --------------------------------------------------------------


def codon_log_odds(coding_dna: list[str], background_seq: str, sample: int = 3_000_000) -> dict[str, float]:
    """log2(P(codon | coding, in frame) / P(codon | chromosome, any frame))."""
    cod: dict[str, int] = {}
    for cds in coding_dna:
        s = cds.upper()
        for i in range(0, len(s) - 2, 3):
            c = s[i : i + 3]
            if "N" not in c:
                cod[c] = cod.get(c, 0) + 1
    bg: dict[str, int] = {}
    b = background_seq.upper()[:sample]
    for i in range(0, len(b) - 2):
        c = b[i : i + 3]
        if "N" not in c:
            bg[c] = bg.get(c, 0) + 1
    tc, tb = sum(cod.values()) + 64, sum(bg.values()) + 64
    return {
        a + x + y: math.log2(((cod.get(a + x + y, 0) + 1) / tc) / ((bg.get(a + x + y, 0) + 1) / tb))
        for a in "ACGT"
        for x in "ACGT"
        for y in "ACGT"
    }


# ---- the parser ----------------------------------------------------------------------


@dataclass(slots=True)
class _Beginning:
    pos: int  # first coding base of the exon
    phase: int  # bases already used of the current codon (0..2) at pos
    score: float
    prev: int  # index of the donor state this exon follows, or -1 for a start
    kind: str  # "start" | "acceptor"


@dataclass(slots=True)
class _Donor:
    pos: int  # first intronic base (the G of GT)
    phase: int  # bases of the current codon used at the exon end
    score: float
    begin: int  # index of the beginning state of the exon that ended here


class SegmentParser:
    """Viterbi over candidate signals. `sites`, when given, replaces the donor/acceptor matrices with an
    external oracle: sites(strand, genomic_offset, n) -> (donors, acceptors) in the string's own coordinates
    (feature c: AlphaGenome splice-site tracks). Start codons always come from the learned Kozak matrix."""

    def __init__(
        self,
        signals: SignalSet,
        codon_lo: dict[str, float],
        min_relative: float = 0.6,
        sites=None,
        start_windows: dict[str, list[tuple[int, int]]] | None = None,
    ) -> None:
        self.signals = signals
        self.codon_lo = codon_lo
        self.min_relative = min_relative
        self.sites = sites
        # optional: a gene may only begin inside these genomic windows per strand (e.g. downstream of
        # a promoter-like element), sorted and non-overlapping
        self.start_windows = start_windows

    def _prefix_and_stops(self, s: str) -> tuple[list[list[float]], list[list[int]]]:
        n = len(s)
        pre = [[0.0] * (n + 1) for _ in range(3)]
        nxt = [[n] * (n + 1) for _ in range(3)]
        lo = self.codon_lo
        for f in range(3):
            acc = 0.0
            row = pre[f]
            for i in range(f, n - 2, 3):
                acc += lo.get(s[i : i + 3], 0.0)
                row[i + 3] = acc
            # fill gaps so that pre[f][p] = coding score of codons in frame f ending at or before p
            last = 0.0
            for p in range(n + 1):
                if p >= 3 and (p - f) % 3 == 0 and p - 3 >= f:
                    last = row[p]
                row[p] = last
            nrow = nxt[f]
            nstop = n
            for i in range(n - 3, -1, -1):
                if (i - f) % 3 == 0 and s[i : i + 3] in STOPS:
                    nstop = i
                nrow[i] = nstop if (i - f) % 3 == 0 else nrow[i + 1] if i + 1 <= n else n
        return pre, nxt

    def _coding(self, pre: list[list[float]], b: int, e: int) -> float:
        """Coding score of whole codons in [b, e) read in the frame that starts at b."""
        f = b % 3
        return pre[f][e] - pre[f][b]

    def _candidates(
        self, s: str, offset: int = 0, strand: str = "+"
    ) -> tuple[list[tuple[int, float]], list[tuple[int, float]], list[tuple[int, float]]]:
        pw = self.signals.pwms
        donor, acceptor, start = pw["splice_donor"], pw["splice_acceptor"], pw["start_kozak"]
        dmax, amax, smax = donor.max_score(), acceptor.max_score(), start.max_score()
        dm, am, sm = donor.matrix(), acceptor.matrix(), start.matrix()
        n = len(s)

        def score_at(m: list[dict[str, float]], anchor: int, offset: int, length: int) -> float:
            i = anchor - offset
            if i < 0 or i + length > n:
                return -math.inf
            tot = 0.0
            for k in range(length):
                ch = s[i + k]
                if ch not in "ACGT":
                    return -math.inf
                tot += m[k][ch]
            return tot

        donors, acceptors, starts = [], [], []
        external = self.sites is not None
        for p in range(1, n - 2):
            two = s[p : p + 2]
            if not external and two == "GT":
                sc = score_at(dm, p, donor.offset, donor.length)
                if sc / dmax >= self.min_relative:
                    donors.append((p, sc))
            if not external and s[p - 1 : p + 1] == "AG":
                sc = score_at(am, p + 1, acceptor.offset, acceptor.length)  # anchor: first exonic base
                if sc / amax >= self.min_relative:
                    acceptors.append((p + 1, sc))  # first exonic base after AG
            if s[p : p + 3] == "ATG":
                sc = score_at(sm, p, start.offset, start.length)
                if sc / smax >= self.min_relative:
                    starts.append((p, sc))
        if external:
            donors, acceptors = self.sites(strand, offset, n)
            donors = [(p, sc) for p, sc in donors if 1 <= p < n - 2]
            acceptors = [(p, sc) for p, sc in acceptors if 1 <= p < n]
        if self.start_windows is not None:
            wins = self.start_windows.get(strand, [])
            lo_ends = [b for _, b in wins]

            def allowed(p: int) -> bool:
                g = offset + p if strand == "+" else offset + (n - 1 - p)
                i = bisect.bisect_right(lo_ends, g)
                return i < len(wins) and wins[i][0] <= g < wins[i][1]

            starts = [(p, sc) for p, sc in starts if allowed(p)]
        return donors, acceptors, starts

    def parse_strand(
        self, s: str, offset: int = 0, strand: str = "+"
    ) -> list[tuple[int, int, list[tuple[int, int]], float]]:
        """Best non-overlapping gene structures on one strand of `s` (coordinates in `s`)."""
        n = len(s)
        pre, nxt = self._prefix_and_stops(s)
        donors, acceptors, starts = self._candidates(s, offset, strand)
        beginnings: list[_Beginning] = []
        donor_states: list[_Donor] = []
        # events in order of position; a donor closes exons that began before it, an acceptor opens one
        events = sorted(
            [(p, 0, sc) for p, sc in starts]
            + [(p, 1, sc) for p, sc in acceptors]
            + [(p, 2, sc) for p, sc in donors]
        )
        # sliding window of donor states per phase: donors with pos in [a - MAX_INTRON, a - MIN_INTRON]
        windows: list[deque[int]] = [deque() for _ in range(3)]
        win_next = 0  # index into donor_states not yet admitted into windows
        genes: list[tuple[int, int, list[tuple[int, int]], float]] = []
        begin_by_pos_from = 0  # beginnings older than MAX_EXON are never needed again
        for pos, kind, sc in events:
            if kind == 0:  # start codon: a new beginning with phase 0
                beginnings.append(_Beginning(pos, 0, sc + START_PRIOR, -1, "start"))
                continue
            if kind == 1:  # acceptor: best donor within intron range, per phase
                while win_next < len(donor_states) and donor_states[win_next].pos <= pos - MIN_INTRON:
                    d = donor_states[win_next]
                    w = windows[d.phase]
                    while w and donor_states[w[-1]].score <= d.score:
                        w.pop()
                    w.append(win_next)
                    win_next += 1
                for ph in range(3):
                    w = windows[ph]
                    while w and donor_states[w[0]].pos < pos - MAX_INTRON:
                        w.popleft()
                    if w:
                        best = donor_states[w[0]]
                        beginnings.append(
                            _Beginning(pos, ph, best.score + sc + INTRON_PRIOR, w[0], "acceptor")
                        )
                continue
            # kind == 2: donor closes an exon that began within [pos - MAX_EXON, pos - MIN_EXON]
            while begin_by_pos_from < len(beginnings) and beginnings[begin_by_pos_from].pos < pos - MAX_EXON:
                begin_by_pos_from += 1
            best_score, best_i, best_phase = -math.inf, -1, 0
            for i in range(begin_by_pos_from, len(beginnings)):
                b = beginnings[i]
                if b.pos > pos - MIN_EXON:
                    break
                first_full = b.pos + (3 - b.phase) % 3
                if first_full >= pos:
                    continue
                # no in-frame stop inside the exon
                if nxt[first_full % 3][first_full] < pos - 2:
                    continue
                total = b.score + self._coding(pre, first_full, pos) + sc
                if total > best_score:
                    best_score, best_i, best_phase = total, i, (b.phase + (pos - b.pos)) % 3
            if best_i >= 0:
                donor_states.append(_Donor(pos, best_phase, best_score, best_i))
        # gene ends: every beginning whose frame runs into a stop within MAX_EXON closes a gene there
        for b in beginnings:
            first_full = b.pos + (3 - b.phase) % 3
            stop = nxt[first_full % 3][first_full] if first_full < n else n
            if stop >= n or stop - b.pos > MAX_EXON or stop - b.pos < MIN_EXON:
                continue
            total = b.score + self._coding(pre, first_full, stop) + STOP_BONUS
            if total < MIN_GENE_SCORE:
                continue
            exons = [(b.pos, stop + 3)]
            j = b.prev
            while j >= 0:
                d = donor_states[j]
                bb = beginnings[d.begin]
                exons.append((bb.pos, d.pos))
                j = bb.prev
            exons.reverse()
            if sum(e - a for a, e in exons) < MIN_CDS:
                continue
            genes.append((exons[0][0], stop + 3, exons, total))
        # keep the best non-overlapping structures with a positive score
        genes.sort(key=lambda g: -g[3])
        taken: list[tuple[int, int]] = []
        out = []
        for g in genes:
            if g[3] <= 0:
                break
            if any(g[0] < e and g[1] > s0 for s0, e in taken):
                continue
            taken.append((g[0], g[1]))
            out.append(g)
        out.sort()
        return out

    def parse(self, chrom: str, seq: str, offset: int = 0) -> list[Prediction]:
        """Both strands of `seq`, which starts at genomic `offset` (used only by an external site oracle)."""
        fwd = seq.upper()
        n = len(fwd)
        preds = [
            Prediction(chrom, Strand.PLUS, s, e, ex, sc)
            for s, e, ex, sc in self.parse_strand(fwd, offset, "+")
        ]
        rc = str(Sequence(fwd).reverse_complement())
        for s, e, ex, sc in self.parse_strand(rc, offset, "-"):
            exons = sorted((n - b, n - a) for a, b in ex)
            preds.append(Prediction(chrom, Strand.MINUS, n - e, n - s, exons, sc))
        preds.sort(key=lambda p: p.start)
        return preds


# ---- evaluation against the annotation -----------------------------------------------


def evaluate(preds: list[Prediction], annotation, chrom: str) -> dict[str, Any]:
    """Exon-, splice-site- and gene-level sensitivity and precision against every coding transcript."""
    module = annotation.to_module("segments")
    true_exons: set[tuple[int, int, str]] = set()
    true_sites: set[tuple[int, str]] = set()
    canon_exons: set[tuple[int, int, str]] = set()  # the canonical transcript's CDS segments only
    true_genes: list[tuple[int, int, str]] = []
    for g in annotation.protein_coding():
        if g.locus.chrom != chrom:
            continue
        strand = g.locus.strand.value
        true_genes.append((g.locus.start, g.locus.end, strand))
        for t in module.entities[g.id].transcripts:
            canonical = "Ensembl_canonical" in t.tags
            for seg in t.cds_segments:
                true_exons.add((seg.start, seg.end, strand))
                true_sites.add((seg.start, strand))
                true_sites.add((seg.end, strand))
                if canonical:
                    canon_exons.add((seg.start, seg.end, strand))
    pred_exons = {(a, b, p.strand.value) for p in preds for a, b in p.exons}
    pred_sites = {x for p in preds for a, b in p.exons for x in ((a, p.strand.value), (b, p.strand.value))}
    exon_tp = len(pred_exons & true_exons)
    canon_tp = len(pred_exons & canon_exons)
    site_tp = len(pred_sites & true_sites)
    gene_hit = sum(
        1 for s, e, st in true_genes if any(p.strand.value == st and p.start < e and p.end > s for p in preds)
    )
    pred_hit = sum(
        1 for p in preds if any(p.strand.value == st and p.start < e and p.end > s for s, e, st in true_genes)
    )

    def frac(a: int, b: int) -> float | None:
        return round(a / b, 4) if b else None

    return {
        "predictions": len(preds),
        "predicted_exons": len(pred_exons),
        "true_exons": len(true_exons),
        "exon_sensitivity": frac(exon_tp, len(true_exons)),
        "exon_precision": frac(exon_tp, len(pred_exons)),
        "canonical_exons": len(canon_exons),
        "canonical_exon_sensitivity": frac(canon_tp, len(canon_exons)),
        "site_sensitivity": frac(site_tp, len(true_sites)),
        "site_precision": frac(site_tp, len(pred_sites)),
        "gene_sensitivity": frac(gene_hit, len(true_genes)),
        "gene_precision": frac(pred_hit, len(preds)),
        "true_genes": len(true_genes),
    }


def parse_chromosome(
    chrom: str,
    seq: str,
    signals: SignalSet,
    annotation,
    min_relative: float = 0.5,
    window: int = 2_000_000,
    overlap: int = 200_000,
    sites=None,
    start_windows: dict[str, list[tuple[int, int]]] | None = None,
) -> dict[str, Any]:
    """Learn coding potential from the annotation, parse the chromosome in windows, evaluate.
    `sites` is an optional external splice-site oracle, `start_windows` an optional restriction of
    where a gene may begin (see SegmentParser)."""
    from genomeos.genome import Genome
    from genomeos.runtime.central_dogma import coding_sequence

    genome = Genome(chrom)
    from genomeos.genome import Chromosome

    genome.chromosomes[chrom] = Chromosome(chrom, Sequence(seq))
    module = annotation.to_module("segments")
    cds = []
    for g in annotation.protein_coding():
        if g.locus.chrom != chrom:
            continue
        for t in module.entities[g.id].transcripts:
            if t.cds_segments and "Ensembl_canonical" in t.tags:
                cds.append(coding_sequence(genome, t).replace("U", "T"))
    lo = codon_log_odds(cds, seq)
    parser = SegmentParser(signals, lo, min_relative, sites, start_windows)
    preds: list[Prediction] = []
    n = len(seq)
    pos = 0
    while pos < n:
        end = min(n, pos + window)
        for p in parser.parse(chrom, seq[pos:end], pos):
            shifted = Prediction(
                chrom, p.strand, p.start + pos, p.end + pos, [(a + pos, b + pos) for a, b in p.exons], p.score
            )
            # keep a structure from the window it lies wholly inside (overlap handled by dedupe below)
            preds.append(shifted)
        pos = end - overlap if end < n else n
    # dedupe overlapping windows: keep the best-scoring of overlapping predictions
    preds.sort(key=lambda p: -p.score)
    kept: list[Prediction] = []
    for p in preds:
        if any(q.strand is p.strand and q.start < p.end and q.end > p.start for q in kept):
            continue
        kept.append(p)
    kept.sort(key=lambda p: p.start)
    ev = evaluate(kept, annotation, chrom)
    return {
        "chrom": chrom,
        "min_relative": min_relative,
        "codon_table_size": len(lo),
        "training_cds": len(cds),
        **ev,
        "predictions_top": [p.to_dict() for p in sorted(kept, key=lambda p: -p.score)[:50]],
        "evidence": EVIDENCE if sites is None else EVIDENCE_PREDICTED_SITES,
        "note": "scored against every coding transcript's CDS; a prediction is a candidate, not an assertion",
    }


def promoter_start_windows(
    ccres, reach: int = 5_000, upstream: int = 500
) -> dict[str, list[tuple[int, int]]]:
    """Where a gene may begin if it begins near an ENCODE promoter-like element: from `upstream` bases
    before the element to `reach` bases after it, in the direction of each strand. Curated evidence,
    no expression data: an element marks a candidate promoter, not which strand it serves, so both
    strands get a window."""
    plus, minus = [], []
    for c in ccres:
        if c.cls != "PLS":
            continue
        plus.append((max(0, c.start - upstream), c.end + reach))
        minus.append((max(0, c.start - reach), c.end + upstream))

    def merge(ws: list[tuple[int, int]]) -> list[tuple[int, int]]:
        ws.sort()
        out: list[tuple[int, int]] = []
        for a, b in ws:
            if out and a <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], b))
            else:
                out.append((a, b))
        return out

    return {"+": merge(plus), "-": merge(minus)}
