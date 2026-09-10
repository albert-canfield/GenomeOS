"""Genomic coordinates: strand and locus.

Lives at the top level so that both the genome engine and BioIR can import it
without a circular dependency. Coordinates are 0-based, half-open.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum


class Strand(StrEnum):
    PLUS = "+"
    MINUS = "-"


@dataclass(frozen=True, slots=True)
class Locus:
    chrom: str
    start: int
    end: int
    strand: Strand = Strand.PLUS

    def __post_init__(self) -> None:
        if self.start < 0 or self.end < self.start:
            raise ValueError(f"invalid locus {self.chrom}:{self.start}-{self.end}")

    @property
    def length(self) -> int:
        return self.end - self.start

    def __str__(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}({self.strand.value})"

    @classmethod
    def parse(cls, text: str) -> Locus:
        """Parse 'chr7:1000000-1000500' or 'chr7:1000000-1000500(-)'."""
        strand = Strand.PLUS
        if text.endswith("(-)"):
            strand, text = Strand.MINUS, text[:-3]
        elif text.endswith("(+)"):
            text = text[:-3]
        chrom, _, span = text.partition(":")
        start, _, end = span.partition("-")
        return cls(chrom, int(start.replace(",", "")), int(end.replace(",", "")), strand)
