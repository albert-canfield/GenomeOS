"""Variant effect on a transcript (task 1.6).

Applies a variant to the transcript's coding sequence in genomic coordinates,
translates reference and alternate, and classifies the consequence with the
Sequence Ontology vocabulary ClinVar uses (synonymous_variant, missense_variant,
nonsense, frameshift_variant, inframe_insertion/deletion, start_lost, stop_lost,
splice_site, intron_variant, 5/3_prime_UTR_variant, non_coding).
"""

from __future__ import annotations

from dataclasses import dataclass

from genomeos.genome.sequence import Locus, Strand
from genomeos.genome.variants import Variant
from genomeos.ir import Transcript

from .central_dogma import STANDARD_CODE, translate

_REF_CACHE: dict[tuple[int, str], list[str]] = {}
_TRANSLATION_CACHE: dict[tuple[tuple[int, str], int], str] = {}


@dataclass(frozen=True, slots=True)
class Effect:
    transcript_id: str
    consequence: str
    protein_change: str = ""  # e.g. R123*, A45V, M1?, fs
    ref_protein_len: int = 0
    alt_protein_len: int = 0
    note: str = ""


def _region(transcript: Transcript, variant: Variant) -> str | None:
    """Non-coding classification, or None if the variant touches the CDS."""
    v_start, v_end = variant.pos, variant.pos + max(1, len(variant.ref))
    if not transcript.exons:
        return "non_coding"
    if not transcript.cds_segments:
        return (
            "non_coding_transcript_exon"
            if any(e.start < v_end and v_start < e.end for e in transcript.exons)
            else "intron_variant"
        )
    if any(c.start < v_end and v_start < c.end for c in transcript.cds_segments):
        return None
    for e in transcript.exons:
        if e.start < v_end and v_start < e.end:
            cds_start = min(c.start for c in transcript.cds_segments)
            cds_end = max(c.end for c in transcript.cds_segments)
            strand = transcript.exons[0].strand
            before = v_end <= cds_start
            if (before and strand is Strand.PLUS) or (
                not before and strand is Strand.MINUS and v_start >= cds_end
            ):
                return "5_prime_UTR_variant"
            return "3_prime_UTR_variant"
    for e in transcript.exons:
        if (
            abs(v_start - e.start) <= 2
            or abs(v_end - e.end) <= 2
            or abs(v_start - e.end) <= 2
            or abs(v_end - e.start) <= 2
        ):
            return "splice_site"
    span_start = min(e.start for e in transcript.exons)
    span_end = max(e.end for e in transcript.exons)
    return "intron_variant" if span_start <= v_start < span_end else "intergenic"


def classify(
    genome, transcript: Transcript, variant: Variant, table: dict[str, str] = STANDARD_CODE
) -> Effect:
    region = _region(transcript, variant)
    if region is not None:
        return Effect(transcript.id, region)

    strand = transcript.cds_segments[0].strand
    segs = sorted(transcript.cds_segments, key=lambda l: l.start)
    v_start, v_end = variant.pos, variant.pos + len(variant.ref)
    ref_allele = variant.ref
    alt_allele = variant.allele(1) or variant.allele(0) or (variant.alts[0] if variant.alts else "")
    delta = len(alt_allele) - len(ref_allele)

    # fetch only the coding segments (plus strand), cached per transcript; edit the one with the variant
    key = (id(genome), transcript.id)
    if key not in _REF_CACHE:
        _REF_CACHE[key] = [str(genome.fetch(Locus(s.chrom, s.start, s.end))) for s in segs]
        if len(_REF_CACHE) > 4096:
            _REF_CACHE.pop(next(iter(_REF_CACHE)))
    ref_pieces = _REF_CACHE[key]
    alt_pieces: list[str] = []
    for s, piece in zip(segs, ref_pieces, strict=True):
        if v_end <= s.start or v_start >= s.end:
            alt_pieces.append(piece)
        elif s.start <= v_start and v_end <= s.end:
            rel = v_start - s.start
            if piece[rel : rel + len(ref_allele)].upper() != ref_allele.upper():
                return Effect(transcript.id, "reference_mismatch", note="VCF REF does not match the genome")
            alt_pieces.append(piece[:rel] + alt_allele + piece[rel + len(ref_allele) :])
        else:
            return Effect(transcript.id, "complex_indel", note="variant crosses a CDS boundary")

    def cds(pieces: list[str]) -> str:
        joined = "".join(pieces)
        if strand is Strand.MINUS:
            from genomeos.genome.sequence import Sequence

            joined = str(Sequence(joined).reverse_complement())
        return joined[transcript.cds_phase :]

    alt_cds = cds(alt_pieces)
    tkey = (key, id(table))
    if tkey not in _TRANSLATION_CACHE:
        _TRANSLATION_CACHE[tkey] = translate(cds(ref_pieces), table=table, to_stop=False, initiator=True)
        if len(_TRANSLATION_CACHE) > 4096:
            _TRANSLATION_CACHE.pop(next(iter(_TRANSLATION_CACHE)))
    ref_full = _TRANSLATION_CACHE[tkey]
    # translate the alternate only from the first changed codon onward
    ref_cds = cds(ref_pieces)
    d = next(
        (k for k, (a, b) in enumerate(zip(ref_cds, alt_cds, strict=False)) if a != b),
        min(len(ref_cds), len(alt_cds)),
    )
    codon0 = d // 3
    if codon0 == 0:
        alt_full = translate(alt_cds, table=table, to_stop=False, initiator=True)
    else:
        alt_full = ref_full[:codon0] + translate(alt_cds[codon0 * 3 :], table=table, to_stop=False)
    ref_p, alt_p = ref_full.split("*")[0], alt_full.split("*")[0]

    def first_diff() -> int:
        for i in range(codon0, min(len(ref_full), len(alt_full))):
            if ref_full[i] != alt_full[i]:
                return i
        return min(len(ref_full), len(alt_full))

    i = first_diff()
    if delta % 3 != 0:
        return Effect(
            transcript.id,
            "frameshift_variant",
            f"{ref_full[i] if i < len(ref_full) else '?'}{i + 1}fs",
            len(ref_p),
            len(alt_p),
        )
    if ref_full == alt_full:
        return Effect(transcript.id, "synonymous_variant", "=", len(ref_p), len(alt_p))
    if i == 0 and ref_full.startswith("M") and not alt_full.startswith("M"):
        return Effect(transcript.id, "start_lost", "M1?", len(ref_p), len(alt_p))
    if delta == 0:
        if i < len(alt_full) and alt_full[i] == "*" and i < len(ref_full) and ref_full[i] != "*":
            return Effect(transcript.id, "nonsense", f"{ref_full[i]}{i + 1}*", len(ref_p), len(alt_p))
        if i < len(ref_full) and ref_full[i] == "*" and i < len(alt_full) and alt_full[i] != "*":
            return Effect(transcript.id, "stop_lost", f"*{i + 1}{alt_full[i]}ext", len(ref_p), len(alt_p))
        if alt_p == ref_p:
            return Effect(transcript.id, "synonymous_variant", "=", len(ref_p), len(alt_p))
        return Effect(
            transcript.id, "missense_variant", f"{ref_full[i]}{i + 1}{alt_full[i]}", len(ref_p), len(alt_p)
        )
    if "*" in alt_full[: len(ref_p)] and len(alt_p) < len(ref_p) - abs(delta) // 3:
        return Effect(
            transcript.id,
            "nonsense",
            f"{ref_full[i]}{i + 1}*",
            len(ref_p),
            len(alt_p),
            note="in-frame indel introduces a stop",
        )
    kind = "inframe_insertion" if delta > 0 else "inframe_deletion"
    return Effect(
        transcript.id,
        kind,
        f"{ref_full[i] if i < len(ref_full) else '?'}{i + 1}{'ins' if delta > 0 else 'del'}",
        len(ref_p),
        len(alt_p),
    )


SEVERITY = [
    "frameshift_variant",
    "nonsense",
    "start_lost",
    "stop_lost",
    "splice_site",
    "inframe_deletion",
    "inframe_insertion",
    "missense_variant",
    "complex_indel",
    "synonymous_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon",
    "intron_variant",
    "non_coding",
    "intergenic",
    "reference_mismatch",
]


def severity(consequence: str) -> int:
    return SEVERITY.index(consequence) if consequence in SEVERITY else len(SEVERITY)


def classify_all(
    genome, transcripts: list[Transcript], variant: Variant, table: dict[str, str] = STANDARD_CODE
) -> list[Effect]:
    """Effects on every transcript, most severe first."""
    effects = [classify(genome, t, variant, table) for t in transcripts]
    effects.sort(key=lambda e: severity(e.consequence))
    return effects
