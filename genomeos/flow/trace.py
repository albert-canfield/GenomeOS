"""Base ↔ transcript ↔ codon ↔ residue mapping for one transcript.

Coordinates: genomic positions are 0-based on the chromosome. mRNA positions
are 0-based in the spliced transcript (5' to 3', so on the minus strand they
run against genomic coordinates). CDS positions are 0-based from the first
base of the start codon after phase trimming. Residues are 1-based, as in
UniProt and in variant notation (p.Asp120Asn).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from genomeos.coords import Locus, Strand
from genomeos.ir.model import Gene, Transcript
from genomeos.runtime.central_dogma import (
    SELENOCYSTEINE_CODE,
    STANDARD_CODE,
    VERTEBRATE_MITOCHONDRIAL_CODE,
    translate_cds,
)

AA3 = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "Q": "Gln",
    "E": "Glu",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
    "*": "Ter",
    "X": "Xaa",
}

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


@dataclass(slots=True)
class ExonMap:
    index: int  # 1-based, in transcript order (5' → 3')
    locus: Locus
    mrna_start: int  # offset of this exon's first base in the mRNA
    mrna_end: int


@dataclass(slots=True)
class CentralDogmaTrace:
    gene: str
    transcript: str
    chrom: str
    strand: Strand
    exons: list[ExonMap]
    mrna: str  # spliced, as RNA
    cds_start: int  # mRNA offset of the first CDS base (after phase trimming); -1 if non-coding
    cds_end: int  # mRNA offset one past the last CDS base
    protein: str
    table_name: str = "standard"
    genomic_to_mrna: dict[int, int] = field(default_factory=dict)
    mrna_to_genomic: list[int] = field(default_factory=list)
    tags: list[str] = field(default_factory=list)

    # ---- sizes
    @property
    def utr5(self) -> int:
        return self.cds_start if self.cds_start >= 0 else 0

    @property
    def utr3(self) -> int:
        return len(self.mrna) - self.cds_end if self.cds_start >= 0 else 0

    @property
    def cds(self) -> str:
        return self.mrna[self.cds_start : self.cds_end] if self.cds_start >= 0 else ""

    # ---- upward: base → residue
    def mrna_position(self, genomic_pos: int) -> int | None:
        return self.genomic_to_mrna.get(genomic_pos)

    def residue_of(self, genomic_pos: int) -> dict[str, Any] | None:
        """Where a genomic base ends up: mRNA offset, region, codon, residue."""
        m = self.mrna_position(genomic_pos)
        if m is None:
            return {"genomic": genomic_pos, "region": "intron" if self._inside(genomic_pos) else "outside"}
        out: dict[str, Any] = {"genomic": genomic_pos, "mrna": m, "base": self.mrna[m]}
        if self.cds_start < 0 or m < self.cds_start:
            out["region"] = "5'UTR" if self.cds_start >= 0 else "non-coding"
            return out
        if m >= self.cds_end:
            out["region"] = "3'UTR"
            return out
        c = m - self.cds_start
        codon_i = c // 3
        codon = self.cds[codon_i * 3 : codon_i * 3 + 3]
        aa = self.protein[codon_i] if codon_i < len(self.protein) else "*"
        out.update(
            {
                "region": "CDS",
                "cds": c,
                "codon_index": codon_i + 1,
                "codon": codon,
                "codon_position": c % 3 + 1,
                "residue": codon_i + 1,
                "aa": aa,
                "aa3": AA3.get(aa, aa),
            }
        )
        return out

    # ---- downward: residue → bases
    def genomic_of_residue(self, residue: int) -> list[int]:
        """The three genomic positions (5'→3' in the mRNA) coding a residue (1-based)."""
        if self.cds_start < 0 or residue < 1:
            return []
        c0 = self.cds_start + (residue - 1) * 3
        return [self.mrna_to_genomic[i] for i in range(c0, min(c0 + 3, len(self.mrna_to_genomic)))]

    # ---- a single-base change, traced upward
    def substitute(self, genomic_pos: int, ref: str, alt: str) -> dict[str, Any]:
        """Consequence of a single-nucleotide change on the plus strand of the genome."""
        where = self.residue_of(genomic_pos) or {}
        region = where.get("region", "outside")
        if region in ("intron", "outside", "5'UTR", "3'UTR", "non-coding"):
            return {"consequence": region, "hgvs_p": None, **where}
        # transcribe the change: on the minus strand the mRNA carries the complement
        r, a = ref.upper(), alt.upper()
        if self.strand is Strand.MINUS:
            r, a = COMPLEMENT.get(r, "N"), COMPLEMENT.get(a, "N")
        r_rna, a_rna = r.replace("T", "U"), a.replace("T", "U")
        m = where["mrna"]
        if self.mrna[m] != r_rna:
            return {"consequence": "reference mismatch", "expected": self.mrna[m], **where}
        c = where["cds"]
        codon_i = c // 3
        codon = list(self.cds[codon_i * 3 : codon_i * 3 + 3])
        codon[c % 3] = a_rna
        table = VERTEBRATE_MITOCHONDRIAL_CODE if self.table_name == "mito" else STANDARD_CODE
        new_aa = table.get("".join(codon), "X")
        old_aa = where["aa"]
        if new_aa == old_aa:
            kind = "synonymous"
        elif new_aa == "*":
            kind = "nonsense"
        elif old_aa == "*":
            kind = "stop_lost"
        elif codon_i == 0:
            kind = "start_lost"
        else:
            kind = "missense"
        return {
            "consequence": kind,
            "codon_before": "".join(self.cds[codon_i * 3 : codon_i * 3 + 3]),
            "codon_after": "".join(codon),
            "aa_before": old_aa,
            "aa_after": new_aa,
            "hgvs_p": f"p.{AA3.get(old_aa, old_aa)}{codon_i + 1}{AA3.get(new_aa, new_aa)}",
            "hgvs_c": f"c.{c + 1}{r}>{a}",
            **where,
        }

    def _inside(self, pos: int) -> bool:
        return any(e.locus.start <= pos < e.locus.end for e in self.exons) or (
            self.exons
            and min(e.locus.start for e in self.exons) <= pos < max(e.locus.end for e in self.exons)
        )

    def to_dict(self, sequences: bool = True) -> dict[str, Any]:
        d: dict[str, Any] = {
            "gene": self.gene,
            "transcript": self.transcript,
            "chrom": self.chrom,
            "strand": self.strand.value,
            "exons": [
                {
                    "index": e.index,
                    "start": e.locus.start,
                    "end": e.locus.end,
                    "mrna_start": e.mrna_start,
                    "mrna_end": e.mrna_end,
                }
                for e in self.exons
            ],
            "mrna_length": len(self.mrna),
            "utr5": self.utr5,
            "cds_start": self.cds_start,
            "cds_end": self.cds_end,
            "utr3": self.utr3,
            "cds_length": len(self.cds),
            "protein_length": len(self.protein),
            "codon_table": self.table_name,
            "tags": self.tags,
            "evidence": {
                "exons": "curated: GENCODE",
                "mrna": "derived: splice of the reference sequence",
                "protein": "derived: translation with the named codon table",
            },
        }
        if sequences:
            d["mrna"] = self.mrna
            d["protein"] = self.protein
            d["mrna_to_genomic"] = self.mrna_to_genomic
        return d


def trace(
    genome, transcript: Transcript, gene_symbol: str = "", table_name: str | None = None
) -> CentralDogmaTrace:
    """Build the base ↔ residue map of one transcript against a genome."""
    if not transcript.exons:
        raise ValueError(f"transcript {transcript.id} has no exons")
    strand = transcript.exons[0].strand
    chrom = transcript.exons[0].chrom
    exons_sorted = sorted(transcript.exons, key=lambda l: l.start, reverse=(strand is Strand.MINUS))
    mrna_parts: list[str] = []
    m2g: list[int] = []
    exon_maps: list[ExonMap] = []
    for i, e in enumerate(exons_sorted, 1):
        start = len(m2g)
        seq = str(genome.fetch(e)).upper().replace("T", "U")
        mrna_parts.append(seq)
        positions = range(e.end - 1, e.start - 1, -1) if strand is Strand.MINUS else range(e.start, e.end)
        m2g.extend(positions)
        exon_maps.append(ExonMap(i, e, start, len(m2g)))
    mrna = "".join(mrna_parts)
    g2m = {g: m for m, g in enumerate(m2g)}
    cds_start = cds_end = -1
    protein = ""
    if transcript.cds_segments:
        segs = sorted(transcript.cds_segments, key=lambda l: l.start, reverse=(strand is Strand.MINUS))
        first = segs[0].end - 1 if strand is Strand.MINUS else segs[0].start
        last = segs[-1].start if strand is Strand.MINUS else segs[-1].end - 1
        if first in g2m and last in g2m:
            cds_start = g2m[first] + transcript.cds_phase
            cds_end = g2m[last] + 1
            if table_name is None:
                table_name = (
                    "mito"
                    if chrom in ("chrM", "MT")
                    else ("seleno" if "seleno" in transcript.tags else "standard")
                )
            table = {"mito": VERTEBRATE_MITOCHONDRIAL_CODE, "seleno": SELENOCYSTEINE_CODE}.get(
                table_name, STANDARD_CODE
            )
            protein = translate_cds(mrna[cds_start:cds_end], table)
    return CentralDogmaTrace(
        gene=gene_symbol or transcript.gene_id,
        transcript=transcript.attrs.get("name", transcript.id),
        chrom=chrom,
        strand=strand,
        exons=exon_maps,
        mrna=mrna,
        cds_start=cds_start,
        cds_end=cds_end,
        protein=protein,
        table_name=table_name or "standard",
        genomic_to_mrna=g2m,
        mrna_to_genomic=m2g,
        tags=list(transcript.tags),
    )


def trace_gene(genome, gene: Gene, transcripts: list[Transcript]) -> CentralDogmaTrace | None:
    """Trace a gene's canonical coding transcript (Ensembl_canonical tag, else the longest CDS)."""
    coding = [t for t in transcripts if t.cds_segments]
    pool = coding or transcripts
    if not pool:
        return None
    canon = next((t for t in pool if "Ensembl_canonical" in t.tags), None)
    if canon is None:
        canon = max(pool, key=lambda t: sum(s.length for s in (t.cds_segments or t.exons)))
    return trace(genome, canon, gene.symbol)
