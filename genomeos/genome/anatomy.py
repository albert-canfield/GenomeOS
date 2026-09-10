"""Genome anatomy: count the blocks and elements of a real chromosome.

Given a sequence and (optionally) its annotation, produce an inventory:

    composition   what the bases are spent on: CDS, UTR, intron, intergenic, gaps
    genes         counts by type, strand balance, transcripts per gene
    structure     exons per transcript, exon/intron/gene lengths, single-exon genes
    layout        gene density, spacing, orientation of neighbours, overlaps, nesting
    elements      CpG islands, homopolymer runs, microsatellites, telomeric blocks
                  (from the sequence alone)

Comparing anatomies across genomes (mitochondrion, human chr21, worm III) is
how the organisation principles show up, and what a genome designed from
zero would have to reproduce.
"""

from __future__ import annotations

import re
import statistics
from dataclasses import dataclass, field
from typing import Any

from genomeos.coords import Strand
from genomeos.genome.sequence import Sequence

CODES = {"intergenic": 0, "intron": 1, "utr": 2, "cds": 3, "gap": 4}
NAMES = {v: k for k, v in CODES.items()}


def _median(xs: list[float]) -> float | None:
    return round(statistics.median(xs), 1) if xs else None


@dataclass(slots=True)
class Anatomy:
    name: str
    length: int
    composition_bp: dict[str, int] = field(default_factory=dict)
    sequence: dict[str, Any] = field(default_factory=dict)
    elements: dict[str, Any] = field(default_factory=dict)
    genes: dict[str, Any] = field(default_factory=dict)
    structure: dict[str, Any] = field(default_factory=dict)
    layout: dict[str, Any] = field(default_factory=dict)

    def composition_fraction(self) -> dict[str, float]:
        return {k: round(v / self.length, 4) for k, v in self.composition_bp.items()} if self.length else {}

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "length": self.length,
            "composition_bp": self.composition_bp,
            "composition_fraction": self.composition_fraction(),
            "sequence": self.sequence,
            "elements": self.elements,
            "genes": self.genes,
            "structure": self.structure,
            "layout": self.layout,
        }


# ---------------------------------------------------------------- sequence ---


def sequence_elements(seq: Sequence | str, window: int = 200) -> tuple[dict[str, Any], dict[str, Any]]:
    """Composition and annotation-free elements of a sequence."""
    s = str(Sequence(str(seq)))
    n = len(s)
    counts = {b: s.count(b) for b in "ACGTN"}
    gc = (counts["G"] + counts["C"]) / max(1, n - counts["N"])
    cpg = s.count("CG")
    cpg_oe = cpg * n / max(1, counts["C"] * counts["G"])
    # GC in 100 kb bins
    bins = [s[i : i + 100_000] for i in range(0, n, 100_000)]
    gcs = []
    for b in bins:
        acgt = len(b) - b.count("N")
        if acgt > 1000:
            gcs.append((b.count("G") + b.count("C")) / acgt)
    # CpG islands: 200 bp windows with GC > 0.5 and obs/exp > 0.6, merged when adjacent
    islands: list[tuple[int, int]] = []
    for i in range(0, n - window + 1, window):
        w = s[i : i + window]
        c, g = w.count("C"), w.count("G")
        if c and g and (c + g) / window > 0.5 and w.count("CG") * window / (c * g) > 0.6:
            if islands and islands[-1][1] == i:
                islands[-1] = (islands[-1][0], i + window)
            else:
                islands.append((i, i + window))
    islands = [(a, b) for a, b in islands if b - a >= 400]
    homopolymers = len(re.findall(r"A{12,}|C{12,}|G{12,}|T{12,}", s))
    micro = len(
        re.findall(
            r"(?:AC){10,}|(?:AG){10,}|(?:AT){10,}|(?:CA){10,}|(?:CT){10,}|(?:GA){10,}|(?:GT){10,}|(?:TA){10,}|(?:TC){10,}|(?:TG){10,}|(?:CG){10,}|(?:GC){10,}",
            s,
        )
    )
    telomeric = re.findall(r"(?:TTAGGG){3,}|(?:CCCTAA){3,}", s)
    gaps = [m.span() for m in re.finditer(r"N{100,}", s)]
    sequence = {
        "gc": round(gc, 4),
        "n_fraction": round(counts["N"] / max(1, n), 4),
        "cpg_obs_exp": round(cpg_oe, 3),
        "gc_100kb_min": round(min(gcs), 3) if gcs else None,
        "gc_100kb_max": round(max(gcs), 3) if gcs else None,
        "gaps": len(gaps),
        "gap_bp": sum(b - a for a, b in gaps),
    }
    elements = {
        "cpg_islands": len(islands),
        "cpg_island_bp": sum(b - a for a, b in islands),
        "homopolymer_runs_ge12": homopolymers,
        "microsatellites_dinuc_ge20bp": micro,
        "telomeric_blocks": len(telomeric),
        "telomeric_bp": sum(len(t) for t in telomeric),
    }
    return sequence, elements


# -------------------------------------------------------------- annotation ---


def annotate_anatomy(anatomy: Anatomy, annotation, chrom: str, seq: Sequence | str) -> None:
    """Fill composition, genes, structure and layout from an Annotation for one chromosome."""
    from collections import Counter

    s = str(seq)
    n = len(s)
    mask = bytearray(n)
    for m in re.finditer(r"N{100,}", s):
        a, b = m.span()
        mask[a:b] = bytes([CODES["gap"]]) * (b - a)

    genes = [g for g in annotation.genes.values() if g.locus.chrom == chrom]
    by_type = Counter(g.type for g in genes)
    strands = Counter(g.locus.strand.value for g in genes)
    pc = [g for g in genes if g.type == "protein_coding"]

    def paint(a: int, b: int, code: int) -> None:
        a, b = max(0, a), min(n, b)
        for i in range(a, b):
            if mask[i] < code and mask[i] != CODES["gap"]:
                mask[i] = code

    tx_per_gene, exons_per_tx, exon_len, intron_len, gene_len, cds_len_aa = [], [], [], [], [], []
    single_exon = 0
    canonical_tx = 0
    for g in pc:
        gene_len.append(g.locus.length)
        txs = [t for t in g.transcripts.values() if t.type == "protein_coding"]
        tx_per_gene.append(len(txs))
        paint(g.locus.start, g.locus.end, CODES["intron"])
        for t in txs:
            for e in t.exons:
                paint(e.start, e.end, CODES["utr"])
            for c, _ in t.cds:
                paint(c.start, c.end, CODES["cds"])
            if "Ensembl_canonical" in t.tags or len(txs) == 1:
                canonical_tx += 1
                ex = sorted(t.exons, key=lambda l: l.start)
                exons_per_tx.append(len(ex))
                if len(ex) == 1:
                    single_exon += 1
                exon_len.extend(e.length for e in ex)
                intron_len.extend(
                    b.start - a.end for a, b in zip(ex, ex[1:], strict=False) if b.start > a.end
                )
                if t.cds:
                    cds_len_aa.append(sum(c.length for c, _ in t.cds) // 3)
    for g in genes:
        if g.type != "protein_coding":
            paint(
                g.locus.start, g.locus.end, CODES["intron"]
            )  # non-coding gene bodies count as "intron" bucket
            for t in g.transcripts.values():
                for e in t.exons:
                    paint(e.start, e.end, CODES["utr"])
    comp = Counter(mask)
    anatomy.composition_bp = {NAMES[k]: comp.get(k, 0) for k in NAMES}
    anatomy.composition_bp["noncoding_exon_or_utr"] = anatomy.composition_bp.pop("utr")

    anatomy.genes = {
        "total": len(genes),
        "by_type": dict(by_type.most_common(12)),
        "strand_plus": strands.get("+", 0),
        "strand_minus": strands.get("-", 0),
        "protein_coding": len(pc),
        "transcripts_per_coding_gene_mean": round(statistics.mean(tx_per_gene), 2) if tx_per_gene else None,
        "transcripts_per_coding_gene_max": max(tx_per_gene) if tx_per_gene else None,
    }
    anatomy.structure = {
        "canonical_transcripts": canonical_tx,
        "exons_per_transcript_median": _median(exons_per_tx),
        "exons_per_transcript_max": max(exons_per_tx) if exons_per_tx else None,
        "single_exon_fraction": round(single_exon / canonical_tx, 3) if canonical_tx else None,
        "exon_length_median": _median(exon_len),
        "intron_length_median": _median(intron_len),
        "intron_length_max": max(intron_len) if intron_len else None,
        "gene_length_median": _median(gene_len),
        "gene_length_max": max(gene_len) if gene_len else None,
        "protein_length_aa_median": _median(cds_len_aa),
        "introns_total": len(intron_len),
    }
    # layout: density, spacing, orientation, overlap, nesting (protein-coding genes)
    ordered = sorted(pc, key=lambda g: g.locus.start)
    spacing, orient, overlaps, nested = [], Counter(), 0, 0
    for a, b in zip(ordered, ordered[1:], strict=False):
        gap = b.locus.start - a.locus.end
        if gap < 0:
            overlaps += 1
            if b.locus.end <= a.locus.end:
                nested += 1
        else:
            spacing.append(gap)
        sa, sb = a.locus.strand, b.locus.strand
        if sa == sb:
            orient["tandem"] += 1
        elif sa is Strand.PLUS and sb is Strand.MINUS:
            orient["convergent"] += 1
        else:
            orient["divergent"] += 1
    non_gap = n - anatomy.composition_bp.get("gap", 0)
    anatomy.layout = {
        "coding_genes_per_mb": round(len(pc) / max(1, non_gap) * 1e6, 2),
        "intergenic_spacing_median": _median(spacing),
        "intergenic_spacing_max": max(spacing) if spacing else None,
        "neighbour_orientation": dict(orient),
        "overlapping_coding_pairs": overlaps,
        "nested_coding_genes": nested,
    }


def anatomy_of(name: str, seq: Sequence | str, annotation=None, chrom: str | None = None) -> Anatomy:
    an = Anatomy(name=name, length=len(seq))
    an.sequence, an.elements = sequence_elements(seq)
    if annotation is not None and chrom is not None:
        annotate_anatomy(an, annotation, chrom, seq)
    else:
        gaps = an.sequence["gap_bp"]
        an.composition_bp = {"unannotated": an.length - gaps, "gap": gaps}
    return an


def design_lessons(anatomies: list[Anatomy]) -> list[str]:
    """Plain statements a genome designer can act on, derived by comparing inventories."""
    out = []
    for a in anatomies:
        f = a.composition_fraction()
        cds = f.get("cds", 0)
        introns = a.structure.get("intron_length_median")
        density = a.layout.get("coding_genes_per_mb")
        exons = a.structure.get("exons_per_transcript_median")
        if not a.genes:
            continue
        out.append(
            f"{a.name}: {a.genes['protein_coding']} coding genes in {a.length / 1e6:.2f} Mb "
            f"({density} per Mb); {cds:.1%} of bases code for protein; "
            f"median transcript has {exons} exons"
            + (f" and introns of {introns:.0f} bp" if introns else " and no introns")
            + f"; neighbours {a.layout.get('neighbour_orientation')}."
        )
    return out
