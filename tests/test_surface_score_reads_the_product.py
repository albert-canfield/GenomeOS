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

import pytest

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
        product_extent_unknown=mech.product_extent_unknown(c),
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


@pytest.mark.xfail(
    strict=True,
    reason=(
        "superseded 2026-10-03 by the 5' partner's registration: a 5' partner keeps the N-terminus "
        "it contributes, but no field records whether the junction falls beyond its own "
        "transmembrane segment, so the curated figure is no longer read at full value for it. The "
        "two claims this test makes that survive - a multiply-reached gene and a gene with no "
        "fusion origin keep the annotation reading - are carried forward live in "
        "test_a_route_other_than_a_fusion_keeps_the_annotation_reading."
    ),
)
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


@pytest.mark.xfail(
    strict=True,
    reason=(
        "this is the plant it was written to be, and it has now fired. The divergence it recorded "
        "was decided on 2026-10-03: the score answers None for a 5' partner, so its third "
        "assertion - today's behaviour of reading the full-length annotation - is superseded. Its "
        "other two claims are carried forward live in "
        "test_a_five_prime_partner_scores_unknown_because_the_junction_is_not_recorded."
    ),
)
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


# --- the 5' partner and the unoriented fusion, registered 2026-10-03 --------------------
#
# The section above closed the 3' partner's uncurated arrangement and left this half
# standing on purpose. These are its plants. Every fixture below varies the *origins*,
# because that is the axis `product_extent_unknown` turns on; the topology is varied
# across the same grid so that nothing here can pass because of a curation field.


def mixed_origins() -> list[VariantOrigin]:
    """The 3' partner of one fusion and the 5' partner of another."""
    return [fusion("3'"), fusion("5'")]


def multipass(signal_peptide: Region | None = SIGNAL) -> Localisation:
    return Localisation(
        compartments={"plasma_membrane": 0.9},
        plasma_membrane=True,
        extracellular_regions=[Region("topo_dom", 19, 1038)],
        transmembrane_regions=[Region("transmem", 1039, 1059), Region("transmem", 1200, 1220)],
        signal_peptide=signal_peptide,
        topology="multi-pass",
    )


def c_terminal(signal_peptide: Region | None = SIGNAL) -> Localisation:
    return Localisation(
        compartments={"plasma_membrane": 0.9},
        plasma_membrane=True,
        extracellular_regions=[Region("topo_dom", 19, 1500)],
        transmembrane_regions=[Region("transmem", 1039, 1059)],
        signal_peptide=signal_peptide,
        topology="ectodomain past the first pass",
    )


def inside() -> Localisation:
    return Localisation(
        compartments={"nucleus": 0.9}, plasma_membrane=False, signal_peptide=None, topology="nuclear"
    )


#: Every curated arrangement in this file, so that an origins-axis guard cannot be
#: satisfied by a topology that answered before it was reached.
ARRANGEMENTS = (
    type_i(),
    type_i(signal_peptide=None),
    multipass(),
    multipass(signal_peptide=None),
    c_terminal(),
    c_terminal(signal_peptide=None),
    Localisation(plasma_membrane=True, compartments={"plasma_membrane": 0.9}),
)


def test_a_five_prime_partner_scores_unknown_because_the_junction_is_not_recorded():
    """The plant. A 5' partner keeps its N-terminus and no field says where the break is.

    The fixture is ALK's own curated arrangement, and ALK is the worked case: the
    recorded annotation of the EML4-ALK fixture is EML4 exons 1-20 with ALK exons
    20-29, and ALK's ectodomain is exons 1-19. A fusion contributing ALK as the 5'
    partner at that junction carries ALK's signal peptide and its whole ectodomain
    and none of its transmembrane segment - an unanchored outward face, which is not
    a surface target and was scored at the maximum.

    The two claims carried forward from the superseded divergence test are asserted
    here too, live: the gate answers None, and a 5' partner is not the 3' partner's
    uncertainty.
    """
    c = candidate(type_i(), [fusion("5'")])
    assert mech.product_extent_unknown(c), "a 5' partner's junction is not on the record"
    assert not mech.ectodomain_uncertain(c), "a 5' partner keeps the N-terminus it contributes"
    assert not mech.ectodomain_lost(c), "a 5' partner's outward face is not established as absent"
    assert mech._surface(c) is None, "the gate does not call a fusion product reachable"
    value, basis = score(c)
    assert value is None, (
        f"the surface score read the annotation at {value} for a product that may be an anchored "
        "receptor and may be an unanchored ectodomain; an unknown is dropped and never defaulted, "
        "least of all to the maximum"
    )
    assert "which part of this gene" in basis and "membrane anchor" in basis


def test_an_unreported_orientation_scores_unknown_and_not_the_annotation():
    """Strictly less is known here than for a 5' partner, so it cannot score higher.

    `read_sv_table` refuses to infer an orientation from the gene order in a
    fusion's name, calling that a naming convention and not a measurement. Reading
    the annotation at full value here while the better-characterised 5' case
    answered None would make the score fall as the record improves, which no
    scoring rule may do.
    """
    for loc in ARRANGEMENTS:
        c = candidate(loc, [fusion("")])
        assert mech.product_extent_unknown(c), f"unreported orientation read as settled ({loc.topology})"
        value, basis = score(c)
        assert value is None, f"annotation read at {value} for an unoriented fusion ({loc.topology})"
        assert "which part of this gene" in basis


def test_a_gene_on_both_ends_of_two_fusions_is_unknown_rather_than_either_answer():
    """The 3' partner of one fusion and the 5' of another. Neither answer is the record's."""
    for loc in ARRANGEMENTS:
        c = candidate(loc, mixed_origins())
        assert mech.product_extent_unknown(c), f"mixed orientations read as settled ({loc.topology})"
        assert not mech.ectodomain_lost(c), "one fusion contributes this gene's N-terminus"
        assert not mech.ectodomain_uncertain(c), "that predicate is about the 3' partner"
        assert score(c)[0] is None


def test_the_three_prime_partners_answers_are_not_swallowed_by_the_new_state():
    """The plant for the `not _n_terminus_lost` clause, and the scope line.

    Removing that clause makes `product_extent_unknown` true for a 3' partner too.
    The score would not move for the first two arrangements, because
    `ectodomain_lost` and `ectodomain_uncertain` are read first - which is exactly
    how a guard passes for the wrong reason. What catches it is the third group: the
    curated multi-pass and C-terminal arrangements, which the registration of the
    section above deliberately left reading the annotation because a described
    arrangement is an answer and not a gap. Those must keep a number.
    """
    curated_type_i = candidate(type_i(), [fusion("3'")])
    assert not mech.product_extent_unknown(curated_type_i)
    assert score(curated_type_i)[0] == 0.0, "ectodomain_lost keeps precedence"

    uncurated_type_i = candidate(type_i(signal_peptide=None), [fusion("3'")])
    assert not mech.product_extent_unknown(uncurated_type_i)
    assert score(uncurated_type_i)[0] is None
    assert "does not settle" in score(uncurated_type_i)[1], "the 3' partner's own reason, not this one"

    for build in (multipass, c_terminal):
        for signal_peptide in (SIGNAL, None):
            loc = build(signal_peptide)
            c = candidate(loc, [fusion("3'")])
            assert not mech.product_extent_unknown(c), (
                f"a described arrangement read as an unrecorded extent ({loc.topology}, signal "
                f"peptide {signal_peptide})"
            )
            assert score(c)[0] is not None, (
                "the section above left this reading the annotation on purpose; this registration "
                "does not reopen it"
            )


def test_a_route_other_than_a_fusion_keeps_the_annotation_reading():
    """The plant for the `_fusion_only` clause, and the carried-forward falsifier.

    Removing that clause leaves `not _n_terminus_lost`, which is true of a candidate
    with no origins at all and of a gene reached by a 5' fusion *and* by a point
    mutation. Both make a full-length product, so both must keep the annotation. The
    3'-plus-mutation case alone would not catch the removal, because its orientation
    set is still exactly {"3'"} - that is the fixture variation this guard needs.
    """
    no_fusion = candidate(type_i(signal_peptide=None), [])
    assert not mech.product_extent_unknown(no_fusion), "no fusion reached this gene"
    assert score(no_fusion)[0] == 1.0

    five_prime_and_mutated = candidate(
        type_i(),
        [fusion("5'"), VariantOrigin(gene="ALK", alteration_kind="point_mutation")],
    )
    assert not mech.product_extent_unknown(five_prime_and_mutated), (
        "a point mutation establishes a full-length product, so the annotation is the right reading"
    )
    assert score(five_prime_and_mutated)[0] == 1.0

    three_prime_and_mutated = candidate(
        type_i(signal_peptide=None),
        [fusion("3'"), VariantOrigin(gene="ALK", alteration_kind="point_mutation")],
    )
    assert not mech.product_extent_unknown(three_prime_and_mutated)
    assert score(three_prime_and_mutated)[0] == 1.0

    unoriented_and_mutated = candidate(
        type_i(),
        [fusion(""), VariantOrigin(gene="ALK", alteration_kind="amplification")],
    )
    assert not mech.product_extent_unknown(unoriented_and_mutated)
    assert score(unoriented_and_mutated)[0] == 1.0


def test_the_four_states_are_pairwise_exclusive():
    """A registered falsifier, over the whole origins-by-topology grid.

    `product_extent_unknown` requires that this gene is *not* established as the 3'
    partner of every fusion, and both of the others require that it is, so no two of
    the three can hold together for any input.
    """
    origin_sets = {
        "3' only": [fusion("3'")],
        "5' only": [fusion("5'")],
        "unreported": [fusion("")],
        "mixed": mixed_origins(),
        "3' plus a mutation": [fusion("3'"), VariantOrigin(gene="ALK", alteration_kind="point_mutation")],
        "5' plus a mutation": [fusion("5'"), VariantOrigin(gene="ALK", alteration_kind="point_mutation")],
        "no origin": [],
    }
    for name, origins in origin_sets.items():
        for loc in (*ARRANGEMENTS, inside()):
            c = candidate(loc, origins)
            states = [
                mech.ectodomain_lost(c),
                mech.ectodomain_uncertain(c),
                mech.product_extent_unknown(c),
            ]
            assert sum(bool(x) for x in states) <= 1, f"{name} / {loc.topology}: {states}"


def test_a_five_prime_partner_curated_inside_the_cell_answers_zero_before_the_fusion_question():
    """The plant for the order of the branches, which mirrors the gate's.

    `plasma_membrane is False` is an answer about the whole protein, so it is read
    before the fusion question. Moving the new branch above it turns a protein
    curated in the nucleus into an open question.
    """
    for origins in ([fusion("5'")], [fusion("")], mixed_origins()):
        c = candidate(inside(), origins)
        assert mech.product_extent_unknown(c), "the fusion question is open on its own terms"
        assert score(c)[0] == 0.0, "a protein curated inside the cell is not an open question"
        assert mech._surface(c) is False


def test_the_pipeline_passes_the_fourth_state_and_does_not_default_it():
    """The call site, because a flag nothing passes is a flag that never fires.

    The dimension has to arrive at `score_candidate` as missing coverage rather than
    averaged in at 1.0, and the overall score may only fall or hold: the value
    removed is the maximum, so dropping it from a mean cannot raise it.
    """
    before = score_candidate(
        candidate(type_i(), [fusion("5'"), VariantOrigin(gene="ALK", alteration_kind="point_mutation")]),
        precedent_available=False,
    )
    c = candidate(type_i(), [fusion("5'")])
    after = score_candidate(c, precedent_available=False)
    component = after.components["surface_accessibility"]
    assert component.value is None, (
        f"the pipeline scored surface accessibility at {component.value} for a product whose extent "
        "it cannot establish"
    )
    assert "which part of this gene" in (component.unknown_reason or component.basis)
    assert after.coverage < before.coverage, "a dropped dimension has to be reported as missing"
    if before.overall is not None and after.overall is not None:
        assert after.overall <= before.overall + 1e-9, (
            f"the overall score rose from {before.overall} to {after.overall} when a dimension worth "
            "the maximum was dropped"
        )


def test_the_gate_and_the_score_now_agree_across_the_fusion_region_with_one_named_exception():
    """The invariant widened to what this registration closed, exception named.

    `_surface` answers None for every fusion-only candidate except the curated
    type-I 3' partner, where it answers False. The score now answers across that
    whole region too, with exactly one disagreement left standing: the 3' partner
    whose curated arrangement is multi-pass or whose extracellular segment runs past
    its first pass. The gate calls that open; the section above registered the score
    reading the annotation there, because a described arrangement is an answer about
    which part of the outward face a 3' partner could still contribute. That is a
    different question from this one and it is not reopened - it is asserted here so
    that it stays countable rather than becoming invisible.
    """
    exceptions = 0
    for origins in ([fusion("3'")], [fusion("5'")], [fusion("")], mixed_origins()):
        for loc in ARRANGEMENTS:
            c = candidate(loc, origins)
            gate, value = mech._surface(c), score(c)[0]
            described_three_prime = mech._n_terminus_lost(c) and not (
                mech.ectodomain_lost(c) or mech.ectodomain_uncertain(c)
            )
            if described_three_prime:
                assert gate is None and value is not None, "the named exception changed shape"
                exceptions += 1
                continue
            if gate is None:
                assert value is None, f"gate open, score {value} ({loc.topology})"
            else:
                assert value is not None, f"gate answered {gate}, score unknown ({loc.topology})"
    assert exceptions == 4, (
        f"the one remaining disagreement is countable and there are {exceptions} of it: the curated "
        "and uncurated multi-pass and C-terminal arrangements, 3' partner only"
    )
