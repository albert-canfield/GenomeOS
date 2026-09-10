"""Nucleotide sequences and genomic coordinates.

Coordinates are 0-based, half-open, like BED and Python slices:
    Locus("chr1", 0, 10) is the first ten bases of chr1.
"""

from __future__ import annotations

from genomeos.coords import Locus, Strand  # noqa: F401  (re-exported)

# IUPAC nucleotide alphabet. N and the ambiguity codes are legal because real
# assemblies contain them; a lossless genome engine must not reject them.
DNA_ALPHABET = frozenset("ACGTNRYKMSWBDHV")
_COMPLEMENT = str.maketrans("ACGTNRYKMSWBDHVacgtnrykmswbdhv", "TGCANYRMKSWVHDBtgcanyrmkswvhdb")


class Sequence:
    """An immutable DNA string with biological operations.

    Stored uppercase. Validation is opt-in because validating a 250 Mb
    chromosome on load is wasteful; call validate() when you need the guarantee.
    """

    __slots__ = ("_s",)

    def __init__(self, s: str) -> None:
        self._s = s.upper()

    def __str__(self) -> str:
        return self._s

    def __len__(self) -> int:
        return len(self._s)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Sequence):
            return self._s == other._s
        if isinstance(other, str):
            return self._s == other.upper()
        return NotImplemented

    def __hash__(self) -> int:
        return hash(self._s)

    def __getitem__(self, item: int | slice) -> Sequence:
        return Sequence(self._s[item])

    def __repr__(self) -> str:
        head = self._s[:30] + ("..." if len(self._s) > 30 else "")
        return f"Sequence({head!r}, len={len(self._s)})"

    def validate(self) -> None:
        bad = set(self._s) - DNA_ALPHABET
        if bad:
            raise ValueError(f"non-IUPAC characters in sequence: {sorted(bad)}")

    def complement(self) -> Sequence:
        return Sequence(self._s.translate(_COMPLEMENT))

    def reverse_complement(self) -> Sequence:
        return Sequence(self._s.translate(_COMPLEMENT)[::-1])

    def gc_content(self) -> float:
        """Fraction of G+C among unambiguous bases (N and ambiguity codes excluded)."""
        acgt = sum(self._s.count(b) for b in "ACGT")
        if acgt == 0:
            return 0.0
        return (self._s.count("G") + self._s.count("C")) / acgt

    def n_fraction(self) -> float:
        return self._s.count("N") / len(self._s) if self._s else 0.0

    def count(self, motif: str) -> int:
        """Count non-overlapping occurrences of a motif on this strand."""
        return self._s.count(motif.upper())

    def telomeric_repeats(self, unit: str = "TTAGGG") -> tuple[int, int]:
        """Count telomeric repeat units on the forward and reverse strands.

        Human telomeres are tandem TTAGGG repeats (5'->3' on the G-rich strand).
        Returns (forward_count, reverse_count).
        """
        rc = Sequence(unit).reverse_complement()._s
        return self._s.count(unit.upper()), self._s.count(rc)
