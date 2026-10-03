# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The surface score's last annotation read, registered 2026-10-03.

`ectodomain_lost` is the conjunction of three clauses, and its third one asks
for a *curated* type-I arrangement on purpose, so that a multi-pass protein, a
C-terminal ectodomain or an uncurated topology answers no and the question
stays open. The mechanism gate honours that and returns `None`. The score was
handed a bool, could not tell "open" from "the ectodomain survives", and read
the full-length annotation at full value — defaulting an unknown to the
maximum, which is the rule at the top of `scoring.py` turned inside out.

These tests are built from constructed topologies rather than from the
benchmark, because no benchmark case reaches this subcase: ALK's signal peptide
is curated, so it answers 0.0 and always did. That is stated in the
registration and is the reason no committed number moves here.
"""

from __future__ import annotations

import genomeos.therapeutics.mechanisms as mech
from genomeos.therapeutics.model import Localisation, Region, TherapeuticTargetCandidate, VariantOrigin
from genomeos.therapeutics.pipeline import score_candidate
from genomeos.therapeutics.scoring import surface_accessibility

SIGNAL = Region("signal_peptide", 1, 18, "Signal peptide")


def type_i(signal_peptide: Region | None = SIGNAL) -> Localisation:
    """ALK's curated arrangement: signal peptide, ectodomain, one pass, cytoplasm."""
    return Localisation(
        compartments={"plasma_membrane": 0.95, "cell_surface": 0.9},
        plasma_membrane=True,
        extracellular_regions=[Region("topo_dom", 19, 1038, "Extracellular")],
        transmembrane_regions=[Region("transmem", 1039, 1059, "Helical")],
        intracellular_regions=[Region("topo_dom", 1060, 1620, "Cytoplasmic")],
        signal_peptide=signal_peptide,
        topology="single-pass type I",
        orientation="N-terminus outside",
    )


def fusion(orientation: str = "3'") -> VariantOrigin:
    o = VariantOrigin(gene="ALK", alteration_kind="fusion")
    o.fusion_orientation = orientation
    return o


def candidate(loc: Localisation, origins: list[VariantOrigin]) -> TherapeuticTargetCandidate:
    return TherapeuticTargetCandidate(gene="ALK", localization=loc, origins=origins)


def score(c: TherapeuticTargetCandidate):
    return surface_accessibility(
        c.localization,
        ectodomain_lost=mech.ectodomain_lost(c),
        ectodomain_uncertain=mech.ectodomain_uncertain(c),
    )


def test_an_uncurated_arrangement_scores_unknown_and_not_the_annotation():
    """The plant. The same gene, the same fusion, one curation field missing.

    Nothing about the product changed between this candidate and the one
    below: it is the 3' partner either way and its outward face is curated
    either way. Only the signal peptide is absent, which is what makes the
    arrangement an open question rather than an answer. The score must say so
    instead of returning the full-length figure at full value.
    """
    c = candidate(type_i(signal_peptide=None), [fusion("3'")])
    assert mech.ectodomain_uncertain(c), "the arrangement is uncurated, so the question is open"
    assert not mech.ectodomain_lost(c), "an uncurated arrangement is not an answer of yes"
    value, basis = score(c)
    assert value is None, (
        f"the surface score read the annotation at {value} for a product whose outward face is "
        "an open question; an unknown dimension is dropped and never defaulted, least of all to "
        "the maximum"
    )
    assert "does not settle" in basis and "3' partner" in basis


def test_a_curated_arrangement_still_answers_zero_and_not_unknown():
    """`ectodomain_lost` keeps precedence. A repair must not soften an answer."""
    c = candidate(type_i(), [fusion("3'")])
    assert mech.ectodomain_lost(c)
    assert not mech.ectodomain_uncertain(c)
    assert score(c)[0] == 0.0


def test_uncertain_and_lost_are_never_both_true():
    """A registered falsifier: the arrangement cannot be curated and uncurated at once."""
    arrangements = [
        type_i(),
        type_i(signal_peptide=None),
        Localisation(plasma_membrane=True, compartments={"plasma_membrane": 0.9}),
        Localisation(
            compartments={"plasma_membrane": 0.9},
            plasma_membrane=True,
            extracellular_regions=[Region("topo_dom", 19, 1038)],
            transmembrane_regions=[Region("transmem", 1039, 1059), Region("transmem", 1200, 1220)],
            signal_peptide=SIGNAL,
        ),
        Localisation(
            compartments={"plasma_membrane": 0.9},
            plasma_membrane=True,
            extracellular_regions=[Region("topo_dom", 19, 1500)],
            transmembrane_regions=[Region("transmem", 1039, 1059)],
            signal_peptide=SIGNAL,
        ),
    ]
    for loc in arrangements:
        for orientation in ("3'", "5'", None):
            c = candidate(loc, [fusion(orientation)] if orientation else [])
            assert not (mech.ectodomain_lost(c) and mech.ectodomain_uncertain(c)), loc.topology


def test_a_curated_multipass_or_c_terminal_ectodomain_is_an_answer_and_not_a_gap():
    """The scope line, because the opposite mistake is the one made in 2026-09-21.

    A protein curated with two passes, or with an extracellular segment running
    past its first one, has been described, and the description says its
    outward face is not only N-terminal: a 3' partner may still contribute part
    of it. That is out of scope of this repair, and leaving the annotation
    reading in place is how the repair stays the width it was registered at.

    Each arrangement appears twice, with and without a curated signal peptide,
    and the second of each pair is the one that makes the exclusion do any
    work. With a signal peptide present the predicate's last clause already
    answers no, so removing either exclusion left this test green — both were
    unreached guards until these two fixtures were added.
    """

    def multipass(signal_peptide):
        return Localisation(
            compartments={"plasma_membrane": 0.9},
            plasma_membrane=True,
            extracellular_regions=[Region("topo_dom", 19, 1038)],
            transmembrane_regions=[Region("transmem", 1039, 1059), Region("transmem", 1200, 1220)],
            signal_peptide=signal_peptide,
        )

    def c_terminal(signal_peptide):
        return Localisation(
            compartments={"plasma_membrane": 0.9},
            plasma_membrane=True,
            extracellular_regions=[Region("topo_dom", 19, 1500)],
            transmembrane_regions=[Region("transmem", 1039, 1059)],
            signal_peptide=signal_peptide,
        )

    for build in (multipass, c_terminal):
        for signal_peptide in (SIGNAL, None):
            loc = build(signal_peptide)
            c = candidate(loc, [fusion("3'")])
            assert not mech.ectodomain_uncertain(c), (
                f"a curated arrangement read as a gap (signal peptide {signal_peptide})"
            )
            assert score(c)[0] is not None, "a curated answer must not be turned into an unknown"


def test_a_five_prime_partner_and_a_multiply_reached_gene_keep_the_annotation():
    """The rule is about which end the gene contributes, and about whether any
    full-length product exists at all.

    A 5' partner keeps the N-terminus it contributes. A gene reached by a
    fusion *and* by another route makes a full-length product too, so the
    curated figure is the right reading for that allele. Without both of these
    the repair would decay into "a fusion is never a surface target", which is
    false.
    """
    five_prime = candidate(type_i(signal_peptide=None), [fusion("5'")])
    assert not mech.ectodomain_uncertain(five_prime)
    assert score(five_prime)[0] == 1.0

    also_mutated = candidate(
        type_i(signal_peptide=None),
        [fusion("3'"), VariantOrigin(gene="ALK", alteration_kind="point_mutation")],
    )
    assert not mech.ectodomain_uncertain(also_mutated)
    assert score(also_mutated)[0] == 1.0

    no_fusion = candidate(type_i(signal_peptide=None), [])
    assert not mech.ectodomain_uncertain(no_fusion)
    assert score(no_fusion)[0] == 1.0


def test_a_protein_curated_inside_the_cell_answers_zero_before_the_question_is_asked():
    """Order matters, and it mirrors the gate's. `plasma_membrane is False` is an
    answer about the whole protein, so it is read before the fusion question and
    gives 0.0 rather than an unknown."""
    inside = Localisation(
        compartments={"nucleus": 0.9}, plasma_membrane=False, signal_peptide=None, topology="nuclear"
    )
    c = candidate(inside, [fusion("3'")])
    assert mech.ectodomain_uncertain(c), "the fusion question is open on its own terms"
    assert score(c)[0] == 0.0, "a protein curated inside the cell is not an open question"
    assert mech._surface(c) is False


def test_the_score_and_the_gate_agree_on_the_question_this_change_is_about():
    """The registered invariant, and it is narrower than "the two always agree".

    This test was first written as the wider claim and it failed, on the case
    the registration had already excluded: a **5' partner** with a fully
    curated arrangement. The gate answers `None` there too, but about a
    different question — it declines to call any fusion product reachable even
    when this gene contributes its own N-terminus, because the junction and the
    product's trafficking are unestablished. The score reads the annotation
    there by registered intent. Closing that gap would turn every 5'-partner
    fusion's surface score into an unknown, which moves numbers and is a
    second decision; it is recorded as open rather than taken here.

    So the invariant asserted is the registered one: for a gene contributed as
    the 3' partner of every fusion that reached it, the score answers exactly
    when the gate answers.
    """
    for loc in (type_i(), type_i(signal_peptide=None)):
        c = candidate(loc, [fusion("3'")])
        gate, value = mech._surface(c), score(c)[0]
        if gate is None:
            assert value is None, f"gate open, score {value}"
        else:
            assert value is not None, f"gate answered {gate}, score unknown"


def test_the_five_prime_divergence_is_named_rather_than_left_uncovered():
    """A disagreement left standing on purpose has to be countable.

    For a 5' partner the gate says the surface question is open and the score
    says the protein is maximally accessible. That is outside the registration
    of 2026-10-03 and is asserted here as today's behaviour, so that closing it
    is a deliberate change with this test as its plant, and so that it cannot
    be mistaken for something this change already handled.
    """
    c = candidate(type_i(), [fusion("5'")])
    assert mech._surface(c) is None, "the gate does not call a fusion product reachable"
    assert not mech.ectodomain_uncertain(c), "a 5' partner keeps the N-terminus it contributes"
    assert score(c)[0] == 1.0, "today's behaviour: the score reads the full-length annotation"


def test_the_pipeline_passes_the_third_state_and_does_not_default_it():
    """The call site, because a flag nothing passes is a flag that never fires.

    `score_candidate` is where the dimension is assembled, so the unknown has
    to arrive there and be counted as missing coverage rather than averaged in
    at 1.0.
    """
    c = candidate(type_i(signal_peptide=None), [fusion("3'")])
    scores = score_candidate(c, precedent_available=False)
    component = scores.components["surface_accessibility"]
    assert component.value is None, (
        f"the pipeline scored surface accessibility at {component.value} for a product whose "
        "outward face it cannot establish"
    )
    assert "does not settle" in (component.unknown_reason or component.basis)

    curated = candidate(type_i(), [fusion("3'")])
    assert (
        score_candidate(curated, precedent_available=False).components["surface_accessibility"].value == 0.0
    )
