"""GFF3 annotation compiler: GENCODE (or any Ensembl-style GFF3) -> BioIR.

Streams the file, keeps only the chromosomes asked for, and produces Gene and
Transcript entities with exons and CDS segments. GENCODE coordinates are
1-based inclusive; BioIR loci are 0-based half-open, converted here once.
"""

from __future__ import annotations

import gzip
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import unquote

from genomeos.ir import Evidence, EvidenceKind, Gene, Module, Transcript

from .sequence import Locus, Strand

GENCODE_FULL = Path("data/reference/gencode.v50.annotation.gff3.gz")
GENCODE_SUBSET = Path("data/results/gencode_v50_chr21_chrM.gff3.gz")  # distilled: chr21 + chrM only


def default_gencode(chroms: set[str] | None = None) -> Path | None:
    """Genome-wide GENCODE if present, else the distilled chr21+chrM subset when it covers `chroms`."""
    if GENCODE_FULL.exists():
        return GENCODE_FULL
    if GENCODE_SUBSET.exists() and (chroms is None or chroms <= {"chr21", "chrM"}):
        return GENCODE_SUBSET
    return None


GENCODE_EVIDENCE = Evidence(EvidenceKind.CURATED, "GENCODE 50 (GRCh38.p14)", note="Ensembl/HAVANA annotation")
ENSEMBL_EVIDENCE = Evidence(EvidenceKind.CURATED, "Ensembl GFF3", note="Ensembl gene build")

# Ensembl-style GFF3 uses SO feature types for transcripts and `biotype`/`Name` attributes
GENE_TYPES = {"gene", "ncRNA_gene", "pseudogene"}
TRANSCRIPT_TYPES = {
    "transcript",
    "mRNA",
    "lnc_RNA",
    "ncRNA",
    "miRNA",
    "tRNA",
    "rRNA",
    "snRNA",
    "snoRNA",
    "scRNA",
    "piRNA",
    "pseudogenic_transcript",
    "unconfirmed_transcript",
    "primary_transcript",
    "V_gene_segment",
    "J_gene_segment",
    "C_gene_segment",
    "D_gene_segment",
    "antisense_RNA",
    "pre_miRNA",
    "SRP_RNA",
    "RNase_MRP_RNA",
    "Y_RNA",
}


def _strip_prefix(value: str) -> str:
    """Ensembl ids look like gene:WBGene00000001 / transcript:Y74C9A.3.1."""
    return (
        value.split(":", 1)[1]
        if ":" in value and value.split(":", 1)[0] in ("gene", "transcript", "CDS")
        else value
    )


@dataclass(slots=True)
class Gff3Row:
    seqid: str
    source: str
    type: str
    start: int  # 1-based inclusive as in the file
    end: int
    strand: str
    phase: int | None
    attrs: dict[str, str]

    @property
    def locus(self) -> Locus:
        return Locus(
            self.seqid, self.start - 1, self.end, Strand.MINUS if self.strand == "-" else Strand.PLUS
        )


def _open(path: str | Path):
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def iter_gff3(
    path: str | Path, chroms: set[str] | None = None, types: set[str] | None = None
) -> Iterator[Gff3Row]:
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if chroms is not None and parts[0] not in chroms:
                continue
            if types is not None and parts[2] not in types:
                continue
            attrs: dict[str, str] = {}
            for kv in parts[8].split(";"):
                k, sep, v = kv.partition("=")
                if sep:
                    attrs[k] = unquote(v)
            yield Gff3Row(
                parts[0],
                parts[1],
                parts[2],
                int(parts[3]),
                int(parts[4]),
                parts[6],
                None if parts[7] == "." else int(parts[7]),
                attrs,
            )


@dataclass(slots=True)
class TranscriptRecord:
    id: str
    gene_id: str
    type: str
    name: str
    locus: Locus
    exons: list[Locus] = field(default_factory=list)
    cds: list[tuple[Locus, int]] = field(default_factory=list)  # (locus, phase)
    tags: list[str] = field(default_factory=list)


@dataclass(slots=True)
class GeneRecord:
    id: str
    symbol: str
    type: str
    locus: Locus
    transcripts: dict[str, TranscriptRecord] = field(default_factory=dict)


class Annotation:
    """Parsed gene models for a set of chromosomes."""

    def __init__(self) -> None:
        self.genes: dict[str, GeneRecord] = {}
        self.by_symbol: dict[str, str] = {}
        self._tx_index: dict[str, str] = {}  # transcript id -> gene id

    @classmethod
    def from_gff3(cls, path: str | Path, chroms: set[str] | None = None) -> Annotation:
        """Parse GENCODE or Ensembl GFF3. Both dialects are handled: GENCODE uses
        `transcript` rows with gene_type/gene_name/transcript_type; Ensembl uses SO
        types (mRNA, ncRNA, ...) with biotype/Name and `gene:`/`transcript:` id prefixes."""
        ann = cls()
        wanted = GENE_TYPES | TRANSCRIPT_TYPES | {"exon", "CDS"}
        for row in iter_gff3(path, chroms, wanted):
            a = row.attrs
            if row.type in GENE_TYPES:
                gid = _strip_prefix(a["ID"])
                symbol = a.get("gene_name") or a.get("Name") or gid
                gtype = a.get("gene_type") or a.get("biotype") or row.type
                ann.genes[gid] = GeneRecord(gid, symbol, gtype, row.locus)
                ann.by_symbol.setdefault(symbol, gid)
            elif row.type in TRANSCRIPT_TYPES:
                gid = _strip_prefix(a.get("Parent", ""))
                if gid not in ann.genes:
                    continue
                tid = _strip_prefix(a["ID"])
                tags = [t for t in a.get("tag", "").split(",") if t]
                ttype = a.get("transcript_type") or a.get("biotype") or row.type
                if ttype == "mRNA":
                    ttype = "protein_coding"
                ann.genes[gid].transcripts[tid] = TranscriptRecord(
                    tid, gid, ttype, a.get("transcript_name") or a.get("Name") or tid, row.locus, tags=tags
                )
                ann._tx_index[tid] = gid
            elif row.type in ("exon", "CDS"):
                tid = _strip_prefix(a.get("Parent", ""))
                gid = ann._tx_index.get(tid)
                if gid is None:
                    continue
                tx = ann.genes[gid].transcripts[tid]
                if row.type == "exon":
                    tx.exons.append(row.locus)
                else:
                    tx.cds.append((row.locus, row.phase or 0))
        for g in ann.genes.values():
            for tx in g.transcripts.values():
                tx.exons.sort(key=lambda l: l.start)
                tx.cds.sort(key=lambda c: c[0].start)
        return ann

    # ---- queries ---------------------------------------------------------

    def gene(self, symbol_or_id: str) -> GeneRecord:
        gid = self.by_symbol.get(symbol_or_id, symbol_or_id)
        return self.genes[gid]

    def protein_coding(self) -> list[GeneRecord]:
        return [g for g in self.genes.values() if g.type == "protein_coding"]

    def summary(self) -> dict[str, int]:
        out: dict[str, int] = {}
        for g in self.genes.values():
            out[g.type] = out.get(g.type, 0) + 1
        return dict(sorted(out.items(), key=lambda kv: -kv[1]))

    # ---- compile to BioIR ------------------------------------------------

    def to_module(self, name: str) -> Module:
        m = Module(name=name)
        for g in self.genes.values():
            gene = Gene(
                id=g.id,
                kind="gene",
                symbol=g.symbol,
                locus=g.locus,
                attrs={"gene_type": g.type},
                evidence=GENCODE_EVIDENCE,
                confidence=0.95 if g.type == "protein_coding" else 0.8,
            )
            m.add(gene)
            for t in g.transcripts.values():
                first_phase = 0
                if t.cds:
                    # phase of the segment where translation starts (strand-aware)
                    first = t.cds[0] if t.locus.strand is Strand.PLUS else t.cds[-1]
                    first_phase = first[1]
                span = None
                if t.cds:
                    span = Locus(
                        t.locus.chrom,
                        min(c.start for c, _ in t.cds),
                        max(c.end for c, _ in t.cds),
                        t.locus.strand,
                    )
                tx = Transcript(
                    id=t.id,
                    kind="transcript",
                    gene_id=g.id,
                    exons=list(t.exons),
                    cds=span,
                    cds_segments=[c for c, _ in t.cds],
                    cds_phase=first_phase,
                    tags=list(t.tags),
                    attrs={"transcript_type": t.type, "name": t.name},
                    evidence=GENCODE_EVIDENCE,
                    confidence=0.9,
                )
                m.add(tx)
                gene.transcripts.append(tx)
        return m
