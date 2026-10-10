# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""What the alteration leaves of the protein, read against the curated topology.

Registered in docs/THERAPEUTICS.md, "The surface score reads the alteration,
pre-registered 2026-10-10". The surface dimension used to answer from the gene's
curated annotation whatever the tumour did to the gene. `ectodomain_lost` put
the product question into the score for a fusion's 3' partner and
`product_extent_unknown` for its 5' partner; this module answers the same
question for everything else the record carries: a premature stop, a frameshift,
an in-frame indel, a splice-site change, a copy-number event.

Two things are kept strictly apart here.

  * The *extent*: which reference residues the tumour's product still contains.
    It comes from the reconstructed altered sequence where one exists, and from
    the HGVS protein notation otherwise, and it says nothing about topology.
  * The *verdict*: the extent held against the curated transmembrane,
    extracellular and signal-peptide features. It has five values, and three of
    them leave the score exactly as it reads today.

`unread` is as important as the rest: where the record curates no membrane
anchor there is nothing for an alteration to remove, the surface figure is
already the discounted one, and replacing a low number with a dropped dimension
would *raise* a candidate's mean. An unknown may not improve a candidate.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any

from genomeos.cancer.tumour import AA1

#: Three-letter to one-letter, and the one-letter codes mapped to themselves, so
#: a notation in either form can be checked against the curated sequence.
_ONE_LETTER = {**AA1, **{v: v for v in AA1.values()}, "*": "*", "X": "X"}

#: Consequence terms whose product may differ in sequence from the reference.
#: Anything outside this set leaves the product unchanged as far as the record
#: says, which is an answer and not a gap: a synonymous or intronic variant
#: makes the same protein.
PRODUCT_CHANGING = frozenset(
    {
        "stop_gained",
        "frameshift_variant",
        "inframe_deletion",
        "inframe_insertion",
        "missense_variant",
        "protein_altering_variant",
        "start_lost",
        "stop_lost",
        "splice_acceptor_variant",
        "splice_donor_variant",
    }
)

#: Terms whose effect on the product this module declines to read, with the
#: reason it declines. A splice-site change alters which exons the transcript
#: carries, and `neoantigen.py` already states that the retained-intron or
#: exon-skipped sequence is not reconstructed; nothing downstream of that can
#: say whether the signal peptide, the outward face or the anchor survives.
DECLINED: dict[str, str] = {
    "splice_acceptor_variant": (
        "a splice-acceptor change alters which exons the transcript carries, and the exon-skipped or "
        "retained-intron sequence is not reconstructed, so whether the product keeps this gene's "
        "signal peptide, outward face or membrane anchor cannot be read"
    ),
    "splice_donor_variant": (
        "a splice-donor change alters which exons the transcript carries, and the exon-skipped or "
        "retained-intron sequence is not reconstructed, so whether the product keeps this gene's "
        "signal peptide, outward face or membrane anchor cannot be read"
    ),
    "start_lost": (
        "translation initiation is lost, and whether a downstream start produces a protein that still "
        "reaches the membrane is not on the record; the curated figure describes the full-length "
        "protein made from the lost start"
    ),
}

#: Copy-number and structural kinds, and what each does to the product's
#: sequence. A deep deletion removes the gene rather than editing a product, and
#: `mechanisms._product_present` already refuses every mechanism for a candidate
#: whose every origin removes the product, so it is not re-answered here.
ALTERATION_KINDS: dict[str, str] = {
    "amplification": (
        "an amplification changes how many copies are transcribed and not the sequence of the product"
    ),
    "gain": "a copy gain changes how many copies are transcribed and not the sequence of the product",
    "shallow_deletion": (
        "a shallow deletion changes how many copies are transcribed and not the sequence of the product"
    ),
    "deep_deletion": (
        "a deep deletion removes the gene rather than editing a product; that the cell makes none of "
        "the protein is answered by the product-present gate and not by this dimension"
    ),
}

#: `p.` notation, in either the three-letter or the one-letter form, with the
#: optional second position of a range and whatever follows it. `hgvsp` is the
#: field to parse and the short `protein_change` form is not: the short form is
#: built by one regex over the first residue only, so `p.Gln1430ArgfsTer12`
#: comes out as `Q1430R` and a frameshift reads as a substitution.
HGVSP = re.compile(
    r"^p\.\(?(?P<a1>[A-Z][a-z]{2}|[A-Z*])(?P<n1>\d+)"
    r"(?:_(?P<a2>[A-Z][a-z]{2}|[A-Z*])(?P<n2>\d+))?(?P<rest>.*?)\)?$"
)

#: Which extent kinds each consequence term can produce. A notation that
#: contradicts its own consequence term is declined rather than preferred: the
#: two disagreeing is a reason to stop, not a tie to break.
EXPECTED_KINDS: dict[str, frozenset[str]] = {
    "stop_gained": frozenset({"truncating"}),
    "frameshift_variant": frozenset({"frameshift"}),
    "inframe_deletion": frozenset({"deletion", "insertion"}),
    "inframe_insertion": frozenset({"deletion", "insertion"}),
    "missense_variant": frozenset({"substitution"}),
    "protein_altering_variant": frozenset({"substitution", "deletion", "insertion"}),
    "stop_lost": frozenset({"insertion", "substitution"}),
}


@dataclass(frozen=True, slots=True)
class Extent:
    """Which reference residues the tumour's product still contains.

    `last_residue` is the last reference residue in the product, inclusive;
    `None` means the product runs to the reference C-terminus. `removed` holds
    reference spans the product does not contain. `novel_tail` says the product
    continues in sequence that has no reference counterpart, which is a
    frameshift's tail and is not characterised here.
    """

    kind: str
    last_residue: int | None = None
    removed: tuple[tuple[int, int], ...] = ()
    novel_tail: bool = False
    first_altered: int | None = None
    basis: str = ""
    source: str = ""


@dataclass(frozen=True, slots=True)
class ProductSurface:
    """What the alteration does to the outward face, and the sentence for it.

    `retained` and `unread` both leave the score reading the annotation, and
    they are different statements: the first says the product keeps the curated
    arrangement, the second says this module has nothing to say about it.
    """

    verdict: str  # retained | reduced | not_anchored | unknown | unread
    reason: str
    extracellular_span: int | None = None
    source: str = ""
    kind: str = ""


#: The sentence carried by every verdict that rests on a premature stop or a
#: frameshift. Whether the transcript is degraded turns on the exon structure
#: and the position of the last exon-exon junction, and `VariantOrigin` carries
#: no exon coordinates, so this change does not predict it.
NMD_UNREAD = (
    "whether the premature stop triggers nonsense-mediated decay, which would leave no product at "
    "all, is not read here: the exon structure and the last exon-exon junction are not carried"
)


# --- the extent, from a sequence or from a notation ------------------------------------


def extent_from_sequence(rebuilt: dict[str, Any]) -> Extent:
    """Read the extent off `neoantigen.altered_protein`'s comparison.

    The orders of business are forced by how that comparison works. A premature
    stop leaves every retained residue identical, so the first difference is
    found only where the mutant has ended: `first_altered_residue` is then
    `mutant_length + 1` and no residue is novel. An in-frame deletion differs
    earlier than that, and the replaced reference span is
    `(reference - mutant) + novel` residues wide from the first difference,
    which covers a plain deletion, a delins and an insertion in one expression.
    """
    ref = int(rebuilt["reference_length"])
    mut = int(rebuilt["mutant_length"])
    first = int(rebuilt["first_altered_residue"])
    novel = int(rebuilt["novel_residues"])
    source = "the reconstructed altered sequence"
    if rebuilt.get("transcript"):
        source = f"the altered {rebuilt['transcript']} translated and compared with the reference"
    if rebuilt.get("frame_shifted"):
        return Extent(
            "frameshift",
            last_residue=max(0, first - 1),
            novel_tail=novel > 0,
            first_altered=first,
            basis=(
                f"read from {source}: the reading frame shifts at residue {first}, so the product "
                f"carries reference "
                f"residues 1-{first - 1} and then {novel} residues that have no reference counterpart"
            ),
            source=source,
        )
    if novel == 0 and first > mut:
        return Extent(
            "truncating",
            last_residue=mut,
            first_altered=first,
            basis=(
                f"read from {source}: the altered sequence translates to {mut} residues against {ref} "
                f"with no residue "
                f"differing, so the product is reference residues 1-{mut}"
            ),
            source=source,
        )
    if ref == mut:
        return Extent(
            "substitution",
            first_altered=first,
            basis=(
                f"read from {source}: {novel or 1} residue(s) replaced at residue {first}; the product "
                "keeps its length"
            ),
            source=source,
        )
    replaced = (ref - mut) + novel
    if replaced > 0:
        span = (first, first + replaced - 1)
        return Extent(
            "deletion",
            removed=(span,),
            first_altered=first,
            basis=(
                f"read from {source}: reference residues {span[0]}-{span[1]} are replaced by {novel} "
                f"residue(s), so the "
                f"product does not contain them"
            ),
            source=source,
        )
    return Extent(
        "insertion",
        first_altered=first,
        basis=(
            f"read from {source}: {novel} residue(s) are inserted after residue {first - 1}; every "
            "curated feature is "
            "still in the product, with the ones after the insertion shifted downstream"
        ),
        source=source,
    )


def extent_from_notation(variant_type: str, hgvsp: str, residue: int | None) -> Extent | None:
    """Read the extent off the HGVS protein notation, or decline.

    `None` means the notation does not establish the extent. The consequence
    term is used as a cross-check and never as a substitute: a notation whose
    own form contradicts its term is declined.
    """
    kind, last, removed, novel_tail, basis = "", None, (), False, ""
    m = HGVSP.match((hgvsp or "").strip())
    if m:
        n1 = int(m.group("n1"))
        n2 = int(m.group("n2")) if m.group("n2") else n1
        rest = m.group("rest") or ""
        if "fs" in rest:
            kind, last, novel_tail = "frameshift", max(0, n1 - 1), True
            basis = (
                f"{hgvsp}: the frame shifts at residue {n1}, so the product carries reference residues "
                f"1-{n1 - 1} and then sequence with no reference counterpart"
            )
        elif rest.startswith("delins"):
            kind, removed = "deletion", ((n1, n2),)
            basis = f"{hgvsp}: reference residues {n1}-{n2} are replaced, so the product lacks them"
        elif rest.startswith("del"):
            kind, removed = "deletion", ((n1, n2),)
            basis = f"{hgvsp}: reference residues {n1}-{n2} are removed from the product"
        elif rest.startswith(("dup", "ins")) or "ext" in rest:
            kind = "insertion"
            basis = (
                f"{hgvsp}: residues are added and none removed, so every curated feature is in the product"
            )
        elif rest in ("Ter", "*", "X"):
            kind, last = "truncating", max(0, n1 - 1)
            basis = f"{hgvsp}: the product ends at reference residue {n1 - 1}"
        elif rest == "=":
            kind = "substitution"
            basis = f"{hgvsp}: the protein sequence is unchanged"
        elif rest in _ONE_LETTER:
            kind = "substitution"
            basis = f"{hgvsp}: one residue is replaced at {n1} and none removed"
        if kind == "":
            return None
        expected = EXPECTED_KINDS.get(variant_type)
        if expected is not None and kind not in expected:
            return None
        return Extent(
            kind,
            last_residue=last,
            removed=removed,
            novel_tail=novel_tail,
            first_altered=n1,
            basis=basis,
            source="the HGVS protein notation",
        )
    # No parsable notation. The consequence term plus the residue index still
    # establishes the extent for a stop and for a frameshift, because both are
    # "the product ends here"; it does not establish an indel's span.
    if residue is None:
        return None
    if variant_type == "stop_gained":
        return Extent(
            "truncating",
            last_residue=max(0, residue - 1),
            first_altered=residue,
            basis=(
                f"a premature stop at residue {residue} from the consequence term and the residue "
                f"index, so the product ends at reference residue {residue - 1}"
            ),
            source="the consequence term and the residue index",
        )
    if variant_type == "frameshift_variant":
        return Extent(
            "frameshift",
            last_residue=max(0, residue - 1),
            novel_tail=True,
            first_altered=residue,
            basis=(
                f"a frameshift at residue {residue} from the consequence term and the residue index, "
                f"so the product carries reference residues 1-{residue - 1} and then sequence with no "
                "reference counterpart"
            ),
            source="the consequence term and the residue index",
        )
    if variant_type in ("missense_variant", "protein_altering_variant"):
        return Extent(
            "substitution",
            first_altered=residue,
            basis=f"one residue is replaced at {residue} and none removed",
            source="the consequence term and the residue index",
        )
    return None


def extent_of(origin: Any, rebuilt: dict[str, Any] | None) -> tuple[Extent | None, ProductSurface | None]:
    """The extent for one origin record, or the verdict that stands in for it.

    The second element is what the candidate's answer uses when the extent
    cannot be read, and the difference between its two values carries weight.
    `unknown` is an open question about a product this module does read: a
    splice-site change, a lost start, an indel whose span the record does not
    give. `unread` says the record is outside what this module reads at all —
    an alteration kind that is not one of the recorded five — and leaves the
    curated annotation to answer as it did before.
    """
    kind = getattr(origin, "alteration_kind", "") or ""
    if kind == "fusion":
        return None, ProductSurface(
            "unknown",
            "which part of this gene a fusion product contains is not on the record; where every "
            "origin is a fusion the three fusion states answer it, and here one origin is not",
            kind="fusion",
        )
    if kind:
        said = ALTERATION_KINDS.get(kind)
        if said is None:
            return None, ProductSurface(
                "unread",
                f"alteration kind '{kind}' is not one this module reads against the curated topology",
                kind=kind,
            )
        return Extent("unchanged", basis=said, source="the alteration record"), None
    variant_type = getattr(origin, "variant_type", "") or ""
    if variant_type in DECLINED:
        return None, ProductSurface("unknown", DECLINED[variant_type], kind=variant_type)
    if variant_type not in PRODUCT_CHANGING:
        return (
            Extent(
                "unchanged",
                basis=(
                    f"consequence '{variant_type or 'none recorded'}' does not change the protein's "
                    "sequence as the record carries it"
                ),
                source="the consequence term",
            ),
            None,
        )
    if rebuilt is not None and _matches(origin, rebuilt):
        return extent_from_sequence(rebuilt), None
    extent = extent_from_notation(
        variant_type, getattr(origin, "hgvsp", "") or "", getattr(origin, "residue", None)
    )
    if extent is None:
        return None, ProductSurface(
            "unknown",
            f"the extent of a '{variant_type}' is not established by this record: "
            f"'{getattr(origin, 'hgvsp', '') or 'no protein notation'}' does not give the span the "
            "product loses, and no reconstructed sequence for it is carried",
            kind=variant_type,
        )
    return extent, None


def _matches(origin: Any, rebuilt: dict[str, Any]) -> bool:
    """Is this reconstruction the one for this origin?

    Keyed on the coding change and the transcript it was applied to, because a
    candidate's reconstruction is built for its lead variant and a candidate can
    carry several.
    """
    hgvsc = rebuilt.get("hgvsc")
    transcript = rebuilt.get("transcript")
    return bool(
        hgvsc
        and transcript
        and (getattr(origin, "hgvsc", "") or "") == hgvsc
        and (getattr(origin, "transcript", "") or "") == transcript
    )


# --- the extent held against the curated topology --------------------------------------


def _coordinates(regions: Any) -> list[tuple[int, int]]:
    return sorted((r.start, r.end) for r in (regions or []) if r.start and r.end)


def _overlap(a: tuple[int, int], b: tuple[int, int]) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)


def surface_from_extent(localisation: Any, extent: Extent, numbering: str) -> ProductSurface:
    """Hold one extent against the curated topology map.

    `numbering` is empty when the extent's residue numbering could not be
    confirmed against the curated sequence, and a definite answer is then capped
    at `unknown`: `hgvsp` numbers the VEP transcript, the topology numbers the
    UniProt canonical sequence, and the repository holds a case where the two
    differ (RUNX1 Y480* on a 481-residue transcript against a 453-residue
    canonical entry).
    """
    tms = _coordinates(localisation.transmembrane_regions)
    ecto = _coordinates(localisation.extracellular_regions)
    if not tms:
        return ProductSurface(
            "unread",
            "no transmembrane segment with coordinates is curated, so the record names no membrane "
            "anchor for this alteration to remove",
            kind=extent.kind,
        )
    if extent.kind in ("unchanged", "substitution", "insertion"):
        return ProductSurface(
            "retained",
            f"the product keeps every curated feature: {extent.basis}",
            source=extent.source,
            kind=extent.kind,
        )
    anchor_lo, anchor_hi = tms[0][0], tms[-1][1]
    outer_hi = max([anchor_hi] + [e[1] for e in ecto])
    nmd = f"; {NMD_UNREAD}" if extent.kind in ("truncating", "frameshift") else ""

    cut = extent.last_residue
    if cut is not None:
        if cut >= outer_hi:
            return ProductSurface(
                "retained",
                f"the product keeps every curated feature, the last of which ends at residue "
                f"{outer_hi}: {extent.basis}{nmd}",
                source=extent.source,
                kind=extent.kind,
            )
        if cut < anchor_lo:
            if extent.kind == "frameshift":
                return ProductSurface(
                    "unknown",
                    f"the product carries reference residues 1-{cut} and then sequence with no "
                    f"reference counterpart, so the curated transmembrane segment at {anchor_lo}-"
                    f"{tms[0][1]} is not in it, and whether the novel stretch anchors the product in "
                    f"the membrane is not characterised: {extent.basis}{nmd}",
                    source=extent.source,
                    kind=extent.kind,
                )
            if not numbering:
                return ProductSurface(
                    "unknown",
                    f"the product would end at reference residue {cut}, before the curated "
                    f"transmembrane segment at {anchor_lo}-{tms[0][1]}, but the residue numbering of "
                    f"'{extent.basis}' could not be confirmed against the curated sequence, so the "
                    f"answer stays open rather than definite{nmd}",
                    source=extent.source,
                    kind=extent.kind,
                )
            return ProductSurface(
                "not_anchored",
                f"the product ends at reference residue {cut}, before the curated transmembrane "
                f"segment at {anchor_lo}-{tms[0][1]}, so it carries no membrane anchor and is not "
                "retained at the cell surface; whether an unanchored outward face is released into "
                f"solution is a different question and is not answered here: {extent.basis}, and "
                f"{numbering}{nmd}",
                source=extent.source,
                kind=extent.kind,
            )
        return ProductSurface(
            "unknown",
            f"the product ends at reference residue {cut}, inside the curated arrangement rather "
            f"than beyond it (transmembrane {anchor_lo}-{anchor_hi}"
            + (f", outward face to {outer_hi}" if ecto else "")
            + "), and whether what survives is still held in the membrane with its outward face "
            f"exposed is not something the record answers: {extent.basis}{nmd}",
            source=extent.source,
            kind=extent.kind,
        )

    if not extent.removed:
        return ProductSurface(
            "retained",
            f"the product keeps every curated feature: {extent.basis}",
            source=extent.source,
            kind=extent.kind,
        )
    whole = [tm for tm in tms if any(_overlap(tm, span) == tm[1] - tm[0] + 1 for span in extent.removed)]
    partial = [tm for tm in tms if tm not in whole and any(_overlap(tm, span) for span in extent.removed)]
    if whole and not partial and len(whole) == len(tms):
        if not numbering:
            return ProductSurface(
                "unknown",
                f"the removed span would contain every curated transmembrane segment, but the "
                f"residue numbering of '{extent.basis}' could not be confirmed against the curated "
                "sequence, so the answer stays open rather than definite",
                source=extent.source,
                kind=extent.kind,
            )
        return ProductSurface(
            "not_anchored",
            f"the alteration removes the whole of the curated transmembrane segment at "
            f"{tms[0][0]}-{tms[0][1]}, so the product carries no membrane anchor and is not retained "
            f"at the cell surface: {extent.basis}, and {numbering}",
            source=extent.source,
            kind=extent.kind,
        )
    if whole or partial:
        return ProductSurface(
            "unknown",
            f"the alteration removes part of the curated transmembrane arrangement, and whether what "
            f"survives still integrates into the membrane is not something the record answers: "
            f"{extent.basis}",
            source=extent.source,
            kind=extent.kind,
        )
    if not ecto:
        return ProductSurface(
            "retained",
            f"the removed span contains no curated transmembrane segment and no extracellular "
            f"segment is curated: {extent.basis}",
            source=extent.source,
            kind=extent.kind,
        )
    longest = max(ecto, key=lambda e: e[1] - e[0])
    gone = sum(_overlap(longest, span) for span in extent.removed)
    if not gone:
        return ProductSurface(
            "retained",
            f"the removed span falls outside the curated outward face and the membrane anchor: "
            f"{extent.basis}",
            source=extent.source,
            kind=extent.kind,
        )
    span = longest[1] - longest[0] + 1
    return ProductSurface(
        "reduced",
        f"the product keeps its membrane anchor and {span - gone} of the {span} curated "
        f"extracellular residues at {longest[0]}-{longest[1]}, {gone} of them being removed by the "
        f"alteration, so the outward face a binder could engage is shorter than the curated one: "
        f"{extent.basis}",
        extracellular_span=max(0, span - gone),
        source=extent.source,
        kind=extent.kind,
    )


# --- the candidate's answer ------------------------------------------------------------

#: Strongest first. The composition takes the first verdict present, so a
#: candidate whose record holds one alteration that keeps the product whole
#: reads the annotation: the missed-target direction, and the behaviour
#: 3fe7c27's 5'-fusion-plus-missense control depends on.
ORDER = ("unread", "retained", "reduced", "unknown", "not_anchored")


def _numbering(origin: Any, extent: Extent, sequence: str) -> str:
    """Is the extent's numbering the curated sequence's? The sentence, or empty."""
    if not sequence:
        return ""
    residue = extent.first_altered
    if residue is None or residue < 1 or residue > len(sequence):
        return ""
    m = HGVSP.match((getattr(origin, "hgvsp", "") or "").strip())
    letter = _ONE_LETTER.get(m.group("a1")) if m else None
    if letter is None:
        letter = (getattr(origin, "protein_change", "") or "")[:1] or None
    if letter is None or letter != sequence[residue - 1]:
        return ""
    return (
        f"the notation's wild-type {letter} at residue {residue} is the residue the curated sequence "
        "carries there, so the two number the same protein"
    )


def _sequence_numbering(rebuilt: dict[str, Any], sequence: str) -> str:
    """The rebuilt protein and the curated one are the same length, or they are not."""
    if not sequence or not rebuilt:
        return ""
    if int(rebuilt.get("reference_length") or 0) != len(sequence):
        return ""
    return (
        f"the rebuilt reference protein is {len(sequence)} residues, the length of the curated "
        "sequence, so the two number the same protein"
    )


def reconstruction(candidate: Any) -> dict[str, Any] | None:
    """The altered-sequence comparison the peptide route already built, if any."""
    neo = getattr(candidate, "neoantigen", None)
    return getattr(neo, "reconstruction", None) if neo is not None else None


def altered_surface(candidate: Any) -> ProductSurface | None:
    """What this tumour's alterations do to the candidate's outward face.

    `None` means there is nothing for this module to say: no origin record at
    all (the expression route reaches a gene the tumour never altered), or every
    origin is a fusion, which the three fusion states own.
    """
    origins = list(getattr(candidate, "origins", None) or [])
    if not origins:
        return None
    if all((getattr(o, "alteration_kind", "") or "") == "fusion" for o in origins):
        return None
    sequence = getattr(candidate, "reference_sequence", "") or ""
    rebuilt = reconstruction(candidate)
    verdicts: list[ProductSurface] = []
    for origin in origins:
        extent, declined = extent_of(origin, rebuilt)
        if extent is None:
            verdicts.append(declined)
            continue
        if extent.source.startswith(("the altered ", "the reconstructed")):
            numbering = _sequence_numbering(rebuilt or {}, sequence)
        else:
            numbering = _numbering(origin, extent, sequence)
        verdicts.append(surface_from_extent(candidate.localization, extent, numbering))
    for verdict in ORDER:
        for v in verdicts:
            if v.verdict == verdict:
                return v
    return None
