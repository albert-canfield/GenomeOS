# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The surface score read from the alteration, registered 2026-10-10.

`surface_accessibility` answered from the gene's curated annotation whatever the
tumour did to the gene. The fusion sections of 2026-10-03 closed that for a
fusion's two ends; these are the plants for every other class the record
carries. The fixtures are constructed topologies rather than benchmark cases on
purpose: no benchmark case carries a truncating, splice-site or in-frame-indel
alteration, which is why the registration predicts that no committed number
moves and why this file is where the change is actually measured.

Each test names the mistake it would catch. The three that matter most:

  * a truncation before the anchor must answer 0.0 and a truncation in the
    cytoplasmic tail must answer exactly what it answered before, or this is a
    blanket penalty on truncating variants rather than a reading of the product;
  * a frameshift must answer `unknown` where a clean stop answers 0.0, because
    the novel tail is not characterised and the two are not the same evidence;
  * nothing may be answered for a protein with no curated anchor, because
    dropping a dimension *raises* a mean when the dropped value is below it, and
    an unknown may not improve a candidate.
"""

from __future__ import annotations

import genomeos.therapeutics.mechanisms as mech
from genomeos.therapeutics import product as prod
from genomeos.therapeutics.model import (
    Localisation,
    NeoantigenAssessment,
    Region,
    TherapeuticTargetCandidate,
    VariantOrigin,
)
from genomeos.therapeutics.pipeline import score_candidate
from genomeos.therapeutics.scoring import surface_accessibility

SIGNAL = Region("signal_peptide", 1, 18, "Signal peptide")


def type_i() -> Localisation:
    """ALK's curated arrangement, the one the fusion sections are built on.

    Signal peptide 1-18, ectodomain to 1038, one pass, cytoplasmic tail to 1620.
    """
    return Localisation(
        compartments={"plasma_membrane": 0.95, "cell_surface": 0.9},
        plasma_membrane=True,
        extracellular_regions=[Region("topo_dom", 19, 1038, "Extracellular")],
        transmembrane_regions=[Region("transmem", 1039, 1059, "Helical")],
        intracellular_regions=[Region("topo_dom", 1060, 1620, "Cytoplasmic")],
        signal_peptide=SIGNAL,
        topology="single-pass type I",
        orientation="N-terminus outside",
    )


def cytoplasmic() -> Localisation:
    """BRAF's shape: a membrane compartment under the threshold, no topology.

    The scope line's fixture. The curated figure here is already the discounted
    one, 0.45 x 0.8 = 0.36, and it sits below a candidate's mean, so answering
    `unknown` would raise the overall score of a gene this tumour truncated.
    """
    return Localisation(
        compartments={"plasma_membrane": 0.45},
        plasma_membrane=None,
        topology="unknown",
    )


def short_ectodomain() -> Localisation:
    """A curated outward face short enough for the span bonus to be under its cap.

    19-400 is 382 residues, worth +0.0764 of the +0.1 maximum, so an in-frame
    deletion inside it can move the number at all.
    """
    return Localisation(
        compartments={"plasma_membrane": 0.9},
        plasma_membrane=True,
        extracellular_regions=[Region("topo_dom", 19, 400, "Extracellular")],
        transmembrane_regions=[Region("transmem", 401, 421, "Helical")],
        intracellular_regions=[Region("topo_dom", 422, 700, "Cytoplasmic")],
        signal_peptide=SIGNAL,
        topology="single-pass type I",
    )


def sequence_for(positions: dict[int, str], length: int = 1620) -> str:
    """A curated sequence of the length `type_i`'s topology numbers.

    The residue at each named position is the one that fixture's notation
    claims, so that the numbering check has something real to agree with.
    """
    seq = ["A"] * length
    for residue, letter in positions.items():
        seq[residue - 1] = letter
    return "".join(seq)


def mutation(
    variant_type: str,
    hgvsp: str = "",
    residue: int | None = None,
    protein_change: str = "",
) -> VariantOrigin:
    return VariantOrigin(
        gene="ALK",
        variant_type=variant_type,
        hgvsp=hgvsp,
        residue=residue,
        protein_change=protein_change,
    )


def fusion(orientation: str = "3'") -> VariantOrigin:
    o = VariantOrigin(gene="ALK", alteration_kind="fusion")
    o.fusion_orientation = orientation
    return o


def candidate(
    loc: Localisation,
    origins: list[VariantOrigin],
    sequence: str = "",
    reconstruction: dict | None = None,
) -> TherapeuticTargetCandidate:
    c = TherapeuticTargetCandidate(gene="ALK", localization=loc, origins=origins)
    c.reference_sequence = sequence
    if reconstruction is not None:
        c.neoantigen = NeoantigenAssessment(reconstruction=reconstruction)
    return c


def score(c: TherapeuticTargetCandidate):
    """The dimension as the pipeline assembles it, all four states passed."""
    return surface_accessibility(
        c.localization,
        ectodomain_lost=mech.ectodomain_lost(c),
        ectodomain_uncertain=mech.ectodomain_uncertain(c),
        product_extent_unknown=mech.product_extent_unknown(c),
        product=prod.altered_surface(c),
    )


def annotation_only(c: TherapeuticTargetCandidate):
    """What the score answered before this change: the gene's annotation."""
    return surface_accessibility(
        c.localization,
        ectodomain_lost=mech.ectodomain_lost(c),
        ectodomain_uncertain=mech.ectodomain_uncertain(c),
        product_extent_unknown=mech.product_extent_unknown(c),
    )


STOP_BEFORE_ANCHOR = sequence_for({500: "Y"})
STOP_IN_TAIL = sequence_for({1401: "Y"})

#: Every class the registration answers, for the two grid tests below.
ALL_CLASSES = [
    mutation("stop_gained", "p.Tyr500Ter", 500),
    mutation("stop_gained", "p.Tyr1401Ter", 1401),
    mutation("frameshift_variant", "p.Tyr500ArgfsTer12", 500),
    mutation("inframe_deletion", "p.Gly1030_Arg1065del", 1030, protein_change="G1030"),
    mutation("inframe_deletion", "p.Gly100_Arg299del", 100, protein_change="G100"),
    mutation("inframe_insertion", "p.Tyr772_Ala775dup", 772, protein_change="Y772"),
    mutation("missense_variant", "p.Leu858Arg", 858, protein_change="L858R"),
    mutation("splice_donor_variant", "", None),
    mutation("start_lost", "p.Met1?", 1),
    mutation("synonymous_variant", "p.Leu858=", 858),
    VariantOrigin(gene="ALK", alteration_kind="amplification"),
    VariantOrigin(gene="ALK", alteration_kind="deep_deletion"),
]

GRID_SEQUENCE = sequence_for({500: "Y", 1401: "Y", 1030: "G", 100: "G", 858: "L", 772: "Y"})


# --- the truncation, where the whole change turns ---------------------------------------


def test_a_stop_before_the_anchor_answers_zero_and_not_the_annotation():
    """The plant. The same gene and the same topology as the fusion fixtures.

    ALK truncated at residue 499 has its signal peptide and most of its outward
    face and none of its transmembrane segment, so the tumour's product is not
    held in the membrane. The curated figure, 1.0, describes a protein this
    tumour does not make.
    """
    c = candidate(type_i(), [mutation("stop_gained", "p.Tyr500Ter", 500)], STOP_BEFORE_ANCHOR)
    assert annotation_only(c)[0] == 1.0, "the fixture is one the old reading scored at the maximum"
    value, basis = score(c)
    assert value == 0.0, f"the score read the annotation at {value} for a product with no anchor"
    assert "before the curated transmembrane segment at 1039-1059" in basis
    assert "carries no membrane anchor" in basis
    assert prod.NMD_UNREAD in basis, "decay is named rather than assumed either way"


def test_a_stop_in_the_cytoplasmic_tail_changes_nothing_at_all():
    """The control that separates a reading from a penalty on truncating variants.

    The same consequence term, the same gene, the same topology: only the
    position moves. A stop at 1401 leaves the signal peptide, the whole outward
    face and the anchor in the product, so the curated figure still describes
    what the tumour makes and the number may not move.
    """
    c = candidate(type_i(), [mutation("stop_gained", "p.Tyr1401Ter", 1401)], STOP_IN_TAIL)
    assert prod.altered_surface(c).verdict == "retained"
    assert score(c) == annotation_only(c), "a tail truncation is not a surface defect"
    assert score(c)[0] == 1.0
    assert prod.NMD_UNREAD in prod.altered_surface(c).reason


def test_a_stop_inside_the_anchor_is_an_open_question_and_not_an_answer():
    """Half a membrane segment is not an answer either way.

    0.0 would claim the shortened segment fails to integrate and 1.0 would
    claim it holds; the record says neither, so the dimension is dropped.
    """
    c = candidate(type_i(), [mutation("stop_gained", "p.Tyr1050Ter", 1050)], sequence_for({1050: "Y"}))
    value, basis = score(c)
    assert value is None, f"the score answered {value} where the cut falls inside the anchor"
    assert "inside the curated arrangement" in basis


def test_a_frameshift_before_the_anchor_is_unknown_where_a_clean_stop_is_zero():
    """The evidence is not the same, so the answer is not the same.

    A premature stop leaves exactly the reference residues before it. A
    frameshift leaves them and then a stretch with no reference counterpart,
    and whether that stretch anchors the product in the membrane is not
    characterised anywhere in this repository. The two must not collapse into
    one answer, in either direction.
    """
    shifted = candidate(
        type_i(), [mutation("frameshift_variant", "p.Tyr500ArgfsTer12", 500)], STOP_BEFORE_ANCHOR
    )
    stopped = candidate(type_i(), [mutation("stop_gained", "p.Tyr500Ter", 500)], STOP_BEFORE_ANCHOR)
    assert score(shifted)[0] is None, "the novel stretch is not characterised"
    assert "not characterised" in score(shifted)[1]
    assert score(stopped)[0] == 0.0, "a clean stop is an answer"


def test_the_frameshift_is_read_from_hgvsp_and_not_from_the_short_form():
    """The trap in the record, which a reader of `protein_change` falls into.

    `cancer.tumour` builds the short form with one regex over the first residue
    only, so `p.Gln1430ArgfsTer12` becomes `Q1430R`: read from the short form a
    frameshift is a substitution, which keeps the product whole. Read from
    `hgvsp` it is a frameshift at 1430, with the reference part ending at 1429
    and a novel tail.
    """
    o = mutation("frameshift_variant", "p.Gln1430ArgfsTer12", 1430, protein_change="Q1430R")
    extent, declined = prod.extent_of(o, None)
    assert declined is None and extent is not None
    assert extent.kind == "frameshift", "the short form would have read this as a substitution"
    assert extent.last_residue == 1429 and extent.novel_tail


# --- the in-frame indel -----------------------------------------------------------------


def test_an_in_frame_deletion_of_the_whole_anchor_answers_zero():
    c = candidate(
        type_i(),
        [mutation("inframe_deletion", "p.Gly1030_Arg1065del", 1030, protein_change="G1030")],
        sequence_for({1030: "G"}),
    )
    value, basis = score(c)
    assert value == 0.0, f"the score answered {value} for a product with its only anchor deleted"
    assert "removes the whole of the curated transmembrane segment at 1039-1059" in basis


def test_an_in_frame_deletion_in_the_outward_face_scores_the_residues_that_survive():
    """The epitope-region case: the anchor holds and the outward face is shorter.

    The bonus is the only part of this dimension that reads a length, so the
    surviving residues are what it reads. For a 382-residue curated face with
    200 removed, the bonus falls from +0.0764 to +0.0364. Nothing about the
    compartment confidence moves, because nothing about the alteration says the
    product left the membrane.
    """
    c = candidate(
        short_ectodomain(),
        [mutation("inframe_deletion", "p.Gly100_Arg299del", 100, protein_change="G100")],
        sequence_for({100: "G"}, length=700),
    )
    surface = prod.altered_surface(c)
    assert surface.verdict == "reduced"
    assert surface.extracellular_span == 182, "382 curated residues less the 200 the deletion removes"
    assert round(annotation_only(c)[0], 4) == 0.9764
    assert round(score(c)[0], 4) == 0.9364
    assert "182 residues of curated extracellular topology that this tumour's product keeps" in score(c)[1]


def test_an_in_frame_deletion_in_the_tail_changes_nothing():
    """The control for the one above: same class, outside the curated face."""
    c = candidate(
        type_i(),
        [mutation("inframe_deletion", "p.Gly1200_Arg1240del", 1200, protein_change="G1200")],
        sequence_for({1200: "G"}),
    )
    assert prod.altered_surface(c).verdict == "retained"
    assert score(c) == annotation_only(c)


def test_an_insertion_and_a_substitution_keep_the_product_whole():
    """Six of the ten benchmark cases are substitutions; none of them may move."""
    for origin in (
        mutation("missense_variant", "p.Leu858Arg", 858, protein_change="L858R"),
        mutation("inframe_insertion", "p.Tyr772_Ala775dup", 772, protein_change="Y772"),
        mutation("stop_lost", "p.Ter1621LeuextTer5", 1621),
        mutation("synonymous_variant", "p.Leu858=", 858),
    ):
        c = candidate(type_i(), [origin], sequence_for({858: "L", 772: "Y"}))
        assert prod.altered_surface(c).verdict == "retained", origin.variant_type
        assert score(c) == annotation_only(c), origin.variant_type


# --- the classes this change declines to read -------------------------------------------


def test_a_splice_site_change_is_unknown_and_says_why():
    """The signal-peptide case. Which exons the transcript keeps is not read."""
    for term in ("splice_donor_variant", "splice_acceptor_variant"):
        c = candidate(type_i(), [mutation(term, "", None)], STOP_BEFORE_ANCHOR)
        value, basis = score(c)
        assert value is None, f"{term} scored {value} from an annotation of the full-length protein"
        assert "exon-skipped or retained-intron sequence is not reconstructed" in basis
        assert "signal peptide" in basis


def test_a_lost_start_is_unknown():
    c = candidate(type_i(), [mutation("start_lost", "p.Met1?", 1)], STOP_BEFORE_ANCHOR)
    assert score(c)[0] is None
    assert "translation initiation is lost" in score(c)[1]


def test_an_indel_whose_span_the_record_does_not_give_is_unknown():
    """A consequence term is not an extent. `inframe_deletion` with no notation."""
    c = candidate(type_i(), [mutation("inframe_deletion", "", 500)], STOP_BEFORE_ANCHOR)
    assert score(c)[0] is None
    assert "does not give the span the product loses" in score(c)[1]


def test_a_notation_that_contradicts_its_consequence_term_is_declined():
    """Two fields disagreeing is a reason to stop, not a tie to break."""
    extent, declined = prod.extent_of(mutation("stop_gained", "p.Leu858Arg", 858), None)
    assert extent is None and declined.verdict == "unknown"
    assert "not established by this record" in declined.reason


def test_a_copy_number_event_leaves_the_product_whole():
    """Two of the ten benchmark cases are copy-number calls."""
    for kind in ("amplification", "gain", "shallow_deletion", "deep_deletion"):
        c = candidate(type_i(), [VariantOrigin(gene="ERBB2", alteration_kind=kind)])
        assert prod.altered_surface(c).verdict == "retained", kind
        assert score(c) == annotation_only(c), kind


# --- the scope line, and the rule that an unknown may not improve a candidate -----------


def test_nothing_is_answered_for_a_protein_with_no_curated_anchor():
    """BRAF's shape, and the reason the scope line exists.

    The curated figure is already discounted to 0.36 and sits below the mean, so
    replacing it with a dropped dimension would *raise* the overall score of a
    gene this tumour truncated. There is also nothing to read: with no curated
    transmembrane segment the record names no anchor for the alteration to
    remove.
    """
    for origin in (
        mutation("stop_gained", "p.Tyr500Ter", 500),
        mutation("frameshift_variant", "p.Tyr500ArgfsTer12", 500),
        mutation("inframe_deletion", "p.Gly100_Arg299del", 100, protein_change="G100"),
    ):
        c = candidate(cytoplasmic(), [origin], STOP_BEFORE_ANCHOR)
        assert prod.altered_surface(c).verdict == "unread", origin.variant_type
        assert score(c) == annotation_only(c), origin.variant_type
        assert round(score(c)[0], 2) == 0.36


def test_no_candidate_s_overall_score_rises_anywhere_on_the_grid():
    """Registered regression condition 5, measured rather than argued.

    This change can remove a value or lower one and never add evidence, so no
    candidate's overall may rise. The grid crosses every alteration class with
    the three arrangements, and each candidate is compared with the same
    candidate whose surface dimension reads the annotation, which is the old
    reading exactly.
    """
    rose = []
    for loc in (type_i(), short_ectodomain(), cytoplasmic()):
        for origin in ALL_CLASSES:
            c = candidate(loc, [origin], GRID_SEQUENCE)
            after = score_candidate(c, precedent_available=False)
            before = score_candidate(candidate(loc, [origin], GRID_SEQUENCE), precedent_available=False)
            # the before-figure: the same assembled score with the dimension put
            # back to what the annotation says, which is what moved here
            value, basis = annotation_only(c)
            before.components["surface_accessibility"].value = value
            before.components["surface_accessibility"].basis = basis
            baseline = _mean(before)
            if after.overall is not None and baseline is not None and after.overall > baseline + 1e-9:
                rose.append((loc.topology, origin.variant_type, baseline, after.overall))
    assert not rose, f"an unknown improved a candidate: {rose}"


def _mean(scores) -> float | None:
    """The weighted mean of a score set, recomputed from its components."""
    from genomeos.therapeutics.scoring import WEIGHTS

    pairs = [
        (WEIGHTS.get(name, 0.0), c.value)
        for name, c in scores.components.items()
        if c.value is not None and name in WEIGHTS
    ]
    total = sum(w for w, _ in pairs)
    return sum(w * v for w, v in pairs) / total if total else None


# --- numbering, which a definite answer may not take on trust ---------------------------


def test_a_definite_answer_requires_the_numbering_confirmed():
    """The RUNX1 shape: the notation numbers a transcript, the topology a protein.

    With no curated sequence to check against, and with a sequence whose residue
    at that position is not the one the notation claims, the same alteration that
    answers 0.0 on confirmed numbering answers `unknown`. The answer may be
    open; it may not be definite on trust.
    """
    origin = mutation("stop_gained", "p.Tyr500Ter", 500)
    confirmed = candidate(type_i(), [origin], STOP_BEFORE_ANCHOR)
    assert score(confirmed)[0] == 0.0
    assert "number the same protein" in score(confirmed)[1]

    for sequence in ("", sequence_for({500: "W"})):
        unconfirmed = candidate(type_i(), [origin], sequence)
        value, basis = score(unconfirmed)
        assert value is None, f"a definite 0.0 was given on numbering that is not confirmed ({value})"
        assert "could not be confirmed against the curated sequence" in basis


def test_the_reconstructed_sequence_is_preferred_and_confirmed_by_length():
    """The measurement beats the notation, and says which transcript it read.

    `neoantigen.altered_protein`'s comparison is already built for the peptide
    route. Here it is the truncation of the same fixture, carried on the
    assessment and keyed on the coding change and the transcript, and the reason
    names the transcript rather than the notation.
    """
    rebuilt = {
        "reference_length": 1620,
        "mutant_length": 499,
        "first_altered_residue": 500,
        "novel_residues": 0,
        "frame_shifted": False,
        "hgvsc": "c.1500T>A",
        "transcript": "ENST00000389048",
    }
    origin = mutation("stop_gained", "p.Tyr500Ter", 500)
    origin.hgvsc, origin.transcript = "c.1500T>A", "ENST00000389048"
    c = candidate(type_i(), [origin], STOP_BEFORE_ANCHOR, reconstruction=rebuilt)
    value, basis = score(c)
    assert value == 0.0
    assert "ENST00000389048" in basis, "the reason names the measurement it read"
    assert "the length of the curated sequence" in basis

    # A reconstruction of another length is another protein's numbering.
    shorter = candidate(type_i(), [origin], sequence_for({}, length=453), reconstruction=rebuilt)
    assert score(shorter)[0] is None
    assert "could not be confirmed" in score(shorter)[1]


def test_a_reconstruction_for_another_variant_is_not_applied_to_this_one():
    """Keyed on the coding change and the transcript, because a candidate has many."""
    rebuilt = {
        "reference_length": 1620,
        "mutant_length": 499,
        "first_altered_residue": 500,
        "novel_residues": 0,
        "frame_shifted": False,
        "hgvsc": "c.1500T>A",
        "transcript": "ENST00000389048",
    }
    origin = mutation("stop_gained", "p.Tyr1401Ter", 1401)
    origin.hgvsc, origin.transcript = "c.4203T>A", "ENST00000389048"
    c = candidate(type_i(), [origin], STOP_IN_TAIL, reconstruction=rebuilt)
    assert prod.altered_surface(c).verdict == "retained", "the other variant's rebuild was applied"
    assert "1401" in prod.altered_surface(c).reason


def test_the_sequence_reader_tells_a_truncation_from_an_in_frame_deletion():
    """Both are shorter with no novel residue, and they are not the same product.

    The comparison distinguishes them by where the first difference falls: a
    premature stop differs only past the end of the mutant protein, an in-frame
    deletion differs at the deletion itself.
    """
    stop = prod.extent_from_sequence(
        {
            "reference_length": 1620,
            "mutant_length": 499,
            "first_altered_residue": 500,
            "novel_residues": 0,
            "frame_shifted": False,
        }
    )
    deletion = prod.extent_from_sequence(
        {
            "reference_length": 1620,
            "mutant_length": 1615,
            "first_altered_residue": 746,
            "novel_residues": 0,
            "frame_shifted": False,
        }
    )
    assert stop.kind == "truncating" and stop.last_residue == 499
    assert deletion.kind == "deletion" and deletion.removed == ((746, 750),)
    assert deletion.last_residue is None, "a deletion does not end the product"


# --- composition over several origins ---------------------------------------------------


def test_one_alteration_that_keeps_the_product_whole_reads_the_annotation():
    """The missed-target direction, and 3fe7c27's control.

    A tumour whose record holds a truncation *and* a missense in the same gene
    makes a product that keeps the outward face from at least one of them, so
    the annotation is read. The opposite composition would lose a real target on
    the strength of one allele.
    """
    c = candidate(
        type_i(),
        [
            mutation("stop_gained", "p.Tyr500Ter", 500),
            mutation("missense_variant", "p.Leu858Arg", 858, protein_change="L858R"),
        ],
        sequence_for({500: "Y", 858: "L"}),
    )
    assert prod.altered_surface(c).verdict == "retained"
    assert score(c)[0] == 1.0

    both = candidate(
        type_i(),
        [mutation("stop_gained", "p.Tyr500Ter", 500), mutation("stop_gained", "p.Tyr600Ter", 600)],
        sequence_for({500: "Y", 600: "Y"}),
    )
    assert score(both)[0] == 0.0, "every origin removes the anchor"


def test_a_fusion_beside_a_mutation_leaves_the_fusion_question_open():
    """The mixed record. The fusion's extent is not on the record either way."""
    c = candidate(
        type_i(),
        [fusion("5'"), mutation("splice_donor_variant", "", None)],
        STOP_BEFORE_ANCHOR,
    )
    assert not mech.product_extent_unknown(c), "the fusion states need every origin to be a fusion"
    assert score(c)[0] is None


def test_the_fusion_region_is_untouched_by_this_change():
    """3fe7c27 and its predecessor keep every published answer they had.

    `altered_surface` returns None for a fusion-only candidate, so the three
    fusion states answer alone: 0.0 for the curated 3' partner, which is the
    benchmark's published ALK row, and None for the 5' partner and the
    unoriented fusion.
    """
    assert score(candidate(type_i(), [fusion("3'")]))[0] == 0.0
    for orientation in ("5'", ""):
        c = candidate(type_i(), [fusion(orientation)])
        assert prod.altered_surface(c) is None, "a fusion-only candidate is not this module's"
        assert score(c)[0] is None

    control = candidate(
        type_i(),
        [fusion("5'"), VariantOrigin(gene="ALK", alteration_kind="point_mutation")],
    )
    assert prod.altered_surface(control).verdict == "unread", "an unrecorded kind is not read"
    assert score(control)[0] == 1.0, "3fe7c27's published control"


def test_a_candidate_the_tumour_never_altered_is_not_read():
    """CD19's route: no origin record, so there is no alteration to read."""
    c = candidate(type_i(), [])
    assert prod.altered_surface(c) is None
    assert score(c) == annotation_only(c)


# --- the gate, which moves with the score -----------------------------------------------


def test_the_gate_and_the_score_agree_across_every_alteration_class():
    """The disagreement the fusion sections closed, not reopened with sides swapped.

    A gate that answered `True` while the score answered 0.0 would offer a
    blocking antibody as established against a product with no membrane anchor.
    """
    for loc in (type_i(), short_ectodomain()):
        for origin in ALL_CLASSES:
            c = candidate(loc, [origin], GRID_SEQUENCE)
            gate, (value, _basis) = mech._surface(c), score(c)
            if value == 0.0:
                assert gate is False, f"gate {gate} with the score at 0.0 ({origin.variant_type})"
            elif value is None:
                assert gate is None, f"gate {gate} with the score unknown ({origin.variant_type})"
            else:
                assert gate is True, f"gate {gate} with the score at {value} ({origin.variant_type})"


def test_the_pipeline_passes_the_reading_and_does_not_default_it():
    """A flag nothing passes is a flag that never fires.

    The dimension has to arrive at `score_candidate` as 0.0 with its reason, and
    the coverage has to fall where it arrives as unknown rather than being
    averaged in at the annotation's value.
    """
    truncated = candidate(type_i(), [mutation("stop_gained", "p.Tyr500Ter", 500)], STOP_BEFORE_ANCHOR)
    component = score_candidate(truncated, precedent_available=False).components["surface_accessibility"]
    assert component.value == 0.0
    assert "carries no membrane anchor" in component.basis

    spliced = candidate(type_i(), [mutation("splice_donor_variant", "", None)], STOP_BEFORE_ANCHOR)
    whole = candidate(type_i(), [], STOP_BEFORE_ANCHOR)
    open_ = score_candidate(spliced, precedent_available=False)
    assert open_.components["surface_accessibility"].value is None
    assert open_.coverage < score_candidate(whole, precedent_available=False).coverage
