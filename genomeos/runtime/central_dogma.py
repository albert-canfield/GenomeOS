"""The central dogma as executable code: DNA -> RNA -> protein.

This is the lowest level of BioVM. It is deterministic; regulation of *when*
and *how much* transcription happens lives in the network runtime.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass

from genomeos.genome.sequence import Locus, Sequence, Strand
from genomeos.ir import Transcript

_BASES = "UCAG"
_AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

# Standard genetic code (NCBI translation table 1), built from the canonical
# UCAG ordering so the table is derived rather than typed by hand.
STANDARD_CODE: dict[str, str] = {
    a + b + c: _AMINO[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}

# Vertebrate mitochondrial code (NCBI table 2) differs at four codons.
VERTEBRATE_MITOCHONDRIAL_CODE = {**STANDARD_CODE, "UGA": "W", "AUA": "M", "AGA": "*", "AGG": "*"}

# Selenoproteins (25 human genes, GENCODE tag "seleno"): an in-frame UGA is read as
# selenocysteine (U) when a SECIS element sits in the 3' UTR; the real stop is UAA or UAG.
SELENOCYSTEINE_CODE = {**STANDARD_CODE, "UGA": "U"}

START_CODON = "AUG"
STOP = "*"

# Initiator codons per NCBI translation table. Whatever the codon, the first
# amino acid of a protein is methionine (formyl-Met in mitochondria).
STANDARD_START_CODONS = frozenset({"AUG"})
VERTEBRATE_MITOCHONDRIAL_START_CODONS = frozenset({"AUU", "AUA", "AUC", "AUG", "GUG"})


def start_codons_for(table: dict[str, str]) -> frozenset[str]:
    return (
        VERTEBRATE_MITOCHONDRIAL_START_CODONS
        if table is VERTEBRATE_MITOCHONDRIAL_CODE
        else STANDARD_START_CODONS
    )


def transcribe(dna: Sequence | str, strand: Strand = Strand.PLUS) -> str:
    """Return the RNA transcribed from the given template orientation.

    The input is the *coding* (sense) sequence as written in the genome for the
    plus strand; for the minus strand we reverse-complement first, which is
    what the polymerase effectively reads.
    """
    seq = Sequence(str(dna))
    if strand is Strand.MINUS:
        seq = seq.reverse_complement()
    return str(seq).replace("T", "U")


def translate(
    rna: str,
    table: dict[str, str] = STANDARD_CODE,
    require_start: bool = False,
    to_stop: bool = True,
    initiator: bool = False,
) -> str:
    """Translate RNA (or DNA; T is accepted) in frame 0.

    Unknown or ambiguous codons (containing N etc.) translate to 'X'.
    With initiator=True the first codon is read as methionine if it is one of
    the table's start codons (needed for mitochondrial AUU/AUA starts).
    """
    rna = rna.upper().replace("T", "U")
    if require_start and not rna.startswith(START_CODON):
        raise ValueError("no start codon at position 0")
    n = len(rna) - len(rna) % 3
    get = table.get
    aas = [get(rna[i : i + 3], "X") for i in range(0, n, 3)]
    if aas and initiator and rna[:3] in start_codons_for(table):
        aas[0] = "M"
    protein = "".join(aas)
    if to_stop:
        stop = protein.find(STOP)
        if stop >= 0:
            protein = protein[:stop]
    return protein


@dataclass(frozen=True, slots=True)
class Orf:
    locus: Locus
    protein: str

    @property
    def length_aa(self) -> int:
        return len(self.protein)


def find_orfs(
    seq: Sequence | str,
    chrom: str = "seq",
    min_aa: int = 30,
    table: dict[str, str] = STANDARD_CODE,
    start_codons: frozenset[str] | None = None,
) -> Iterator[Orf]:
    """Yield open reading frames (start ... stop) on both strands, all three frames.

    Start codons default to the table's initiators (AUG for the standard code;
    AUU/AUA/AUC/AUG/GUG for vertebrate mitochondria). Loci are reported in
    forward-strand coordinates. An ORF is emitted only if it is closed by a
    stop codon inside the sequence, and the initiator is emitted as M.
    """
    starts = start_codons if start_codons is not None else start_codons_for(table)
    forward = str(Sequence(str(seq)))
    n = len(forward)
    for strand in (Strand.PLUS, Strand.MINUS):
        s = forward if strand is Strand.PLUS else str(Sequence(forward).reverse_complement())
        rna = s.replace("T", "U")
        for frame in range(3):
            i = frame
            while i + 3 <= n:
                if rna[i : i + 3] in starts:
                    j = i
                    protein: list[str] = ["M"]
                    closed = False
                    j += 3
                    while j + 3 <= n:
                        aa = table.get(rna[j : j + 3], "X")
                        if aa == STOP:
                            closed = True
                            break
                        protein.append(aa)
                        j += 3
                    if closed and len(protein) >= min_aa:
                        start, end = i, j + 3
                        if strand is Strand.MINUS:
                            start, end = n - end, n - start
                        yield Orf(Locus(chrom, start, end, strand), "".join(protein))
                    i = j + 3 if closed else n
                else:
                    i += 3


def splice(genome, transcript: Transcript) -> str:
    """Join a transcript's exons from the genome into mature mRNA (as RNA)."""
    if not transcript.exons:
        raise ValueError(f"transcript {transcript.id} has no exons")
    strand = transcript.exons[0].strand
    exons = sorted(transcript.exons, key=lambda l: l.start, reverse=(strand is Strand.MINUS))
    parts = [str(genome.fetch(e)) for e in exons]
    return "".join(parts).replace("T", "U")


def coding_sequence(genome, transcript: Transcript) -> str:
    """Join a transcript's CDS segments (strand-aware, phase-trimmed) into the
    coding sequence as RNA. Empty string if the transcript has no CDS."""
    if not transcript.cds_segments:
        return ""
    strand = transcript.cds_segments[0].strand
    segs = sorted(transcript.cds_segments, key=lambda l: l.start, reverse=(strand is Strand.MINUS))
    rna = "".join(str(genome.fetch(seg)) for seg in segs).replace("T", "U")
    return rna[transcript.cds_phase :]


def table_for(transcript: Transcript, chrom: str | None = None) -> dict[str, str]:
    """The codon table a transcript is read with: mitochondrial on chrM, selenocysteine when tagged."""
    if chrom in ("chrM", "MT"):
        return VERTEBRATE_MITOCHONDRIAL_CODE
    if "seleno" in getattr(transcript, "tags", ()):
        return SELENOCYSTEINE_CODE
    return STANDARD_CODE


def translate_cds(rna: str, table: dict[str, str] = STANDARD_CODE) -> str:
    """Translate a coding sequence; with the selenocysteine table a final UGA is still the stop."""
    protein = translate(rna, table=table, initiator=True)
    rna = rna.upper().replace("T", "U")
    # with the selenocysteine table a UGA in the last codon of the CDS is the terminator, not Sec
    last = rna[len(protein) * 3 - 3 : len(protein) * 3]
    if (
        table is SELENOCYSTEINE_CODE
        and protein.endswith("U")
        and last == "UGA"
        and len(rna) - len(protein) * 3 < 3
    ):
        protein = protein[:-1]
    return protein


def translate_transcript(genome, transcript: Transcript, table: dict[str, str] | None = None) -> str:
    """Protein encoded by a transcript, stopping at the first stop codon."""
    chrom = transcript.exons[0].chrom if transcript.exons else None
    return translate_cds(coding_sequence(genome, transcript), table or table_for(transcript, chrom))
