"""Genome container.

Design note: v0.1 stores one linear sequence per chromosome name. The design
target (see docs/ARCHITECTURE.md) is a pangenome graph in which an individual's
chromosome is a path. The API below is written so that fetch(locus) can later be
served by a graph-backed store without changing callers.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from .fasta import iter_fasta
from .sequence import Locus, Sequence, Strand


@dataclass(slots=True)
class Chromosome:
    name: str
    sequence: Sequence
    description: str = ""

    @property
    def length(self) -> int:
        return len(self.sequence)


@dataclass(slots=True)
class Genome:
    """A set of chromosomes for one individual or one reference assembly.

    `ploidy_label` records which haplotype this object represents when a diploid
    genome is stored as two Genome objects (e.g. "maternal", "paternal").
    """

    name: str
    chromosomes: dict[str, Chromosome] = field(default_factory=dict)
    ploidy_label: str = "reference"

    @classmethod
    def from_fasta(cls, path: str | Path, name: str | None = None) -> Genome:
        g = cls(name=name or Path(path).name)
        for chrom_name, desc, seq in iter_fasta(path):
            g.chromosomes[chrom_name] = Chromosome(chrom_name, seq, desc)
        return g

    def __len__(self) -> int:
        return sum(c.length for c in self.chromosomes.values())

    def fetch(self, locus: Locus) -> Sequence:
        """Return the sequence at a locus, reverse-complemented for the minus strand."""
        chrom = self.chromosomes[locus.chrom]
        if locus.end > chrom.length:
            raise IndexError(f"{locus} exceeds {locus.chrom} length {chrom.length}")
        seq = chrom.sequence[locus.start : locus.end]
        return seq.reverse_complement() if locus.strand is Strand.MINUS else seq

    def summary(self) -> list[dict]:
        rows = []
        for c in self.chromosomes.values():
            rows.append(
                {
                    "chromosome": c.name,
                    "length_bp": c.length,
                    "gc": round(c.sequence.gc_content(), 4),
                    "n_fraction": round(c.sequence.n_fraction(), 4),
                }
            )
        return rows
