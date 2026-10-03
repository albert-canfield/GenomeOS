# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Where the ranking reads the weighted mean, and what the mean does with a
known-absent dimension against an unestablished one. Registered 2026-10-03 in
docs/THERAPEUTICS.md under "Does the known-absent surface actually invert a
ranking".

The measurement these guards hold in place: over the ten committed benchmark
cases no candidate carries an unestablished surface, so no pair of the
registered form exists there and no ranking inverts. That is a fact about those
ten tumours and not about the ordering rule, and the difference is the whole
point of this file. The rule does let the mean decide: the evidence tier enters
`_rank_key` as a two-class partition rather than as its three tiers, the
mechanism gate puts a known-absent surface and an unestablished one in the SAME
class, and inside one such cell nothing but the mean and the gene name is left.

Built from constructed candidates rather than from the benchmark, because the
benchmark does not reach the pair: its 42 candidates hold seven known-absent
surfaces and no unestablished one.
"""

from __future__ import annotations

import genomeos.therapeutics.mechanisms as mech
from genomeos.therapeutics.model import Localisation, Region, TherapeuticTargetCandidate, VariantOrigin
from genomeos.therapeutics.pipeline import _rank_key, score_candidate, select_mechanisms
from genomeos.therapeutics.scoring import (
    EVIDENCE_TIERS,
    MECHANISM_GATE_ORDER,
    REACH_NONE,
    UNMEASURED_TIER,
    WEIGHTS,
    assemble,
    component,
    mechanism_reach,
    surface_accessibility,
)

SIGNAL = Region("signal_peptide", 1, 18, "Signal peptide")


def type_i() -> Localisation:
    """A curated single-pass type-I receptor: signal peptide, ectodomain, one pass."""
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


def fusion(orientation: str) -> VariantOrigin:
    o = VariantOrigin(gene="FUS1", alteration_kind="fusion")
    o.fusion_orientation = orientation
    return o


def scored(gene: str, orientation: str) -> TherapeuticTargetCandidate:
    """One candidate, carried through the stages the ranking reads.

    The only thing varied between the two candidates this file compares is the
    recorded fusion orientation, which is the axis the surface dimension turns
    on: `3'` is an answer of no and `` is no record at all.
    """
    c = TherapeuticTargetCandidate(gene=gene, localization=type_i(), origins=[fusion(orientation)])
    c.scores = score_candidate(c, False)
    c.therapeutic_mechanisms = select_mechanisms(c)
    c.mechanism_reach, c.mechanism_reach_reason = mechanism_reach(c)
    return c


def test_the_tier_enters_the_rank_key_as_two_classes_and_not_as_three_tiers():
    """The plant. Three tiers exist; the key carries two classes of them.

    If the key ordered on the tier itself, a `patient_measurement` candidate
    could never be ranked by the mean against an `observed_alteration` one, and
    the mean would reach far fewer rankings than it does. It carries the coarse
    cut instead, which is registered and argued, and this guard is what keeps
    the reading of it honest.
    """
    assert len(EVIDENCE_TIERS) == 3, "three tiers are declared"
    firsts = {}
    for tier in EVIDENCE_TIERS:
        c = TherapeuticTargetCandidate(gene="G", localization=type_i())
        c.scores = score_candidate(c, False)
        c.evidence_tier = tier
        c.mechanism_reach = REACH_NONE
        firsts[tier] = _rank_key(c)[0]
    assert firsts["observed_alteration"] == firsts["patient_measurement"], (
        "two tiers that are both measurements of this gene in this patient must share the rank "
        "key's first element, or the mean never compares across them"
    )
    assert firsts[UNMEASURED_TIER] != firsts["observed_alteration"], (
        "the unmeasured tier must be the one the first element separates"
    )
    assert len(set(firsts.values())) == 2, (
        f"the rank key's first element took {len(set(firsts.values()))} values over three tiers; "
        "the registered partition is two classes"
    )


def test_a_known_absent_surface_and_an_unestablished_one_share_the_gate_class():
    """The two classes the gate does NOT separate, which is why the mean decides.

    Varied on the orientation record and on nothing else. A recorded 3' end is
    an answer of no, so the surface requirement is a failed hard requirement; no
    record at all leaves it unanswered, so the requirement is provisional. The
    registered rule puts a provisional requirement and a failed one in one
    class, so both candidates reach the mean together.
    """
    absent = scored("AAA", "3'")
    unestablished = scored("ZZZ", "")
    assert absent.scores.components["surface_accessibility"].value == 0.0
    assert unestablished.scores.components["surface_accessibility"].value is None
    assert absent.mechanism_reach == REACH_NONE
    assert unestablished.mechanism_reach == REACH_NONE, (
        "an unanswered surface requirement must not be an established mechanism; three gate "
        "classes were argued against and refused when the gate was registered"
    )
    assert _rank_key(absent)[1] == _rank_key(unestablished)[1] == MECHANISM_GATE_ORDER.index(REACH_NONE), (
        "the gate class must be equal for the two, or this pair never reaches the mean at all"
    )


def cell_pair():
    """Two candidates differing in one component and in nothing else.

    The component lists are identical except for the surface dimension, which
    one candidate has as a measured `0.0` and the other has as unestablished.
    Normal-tissue safety is supplied on both at 0.8 so that neither safety cap
    fires: a cap clips at a constant and would hide the very gap being measured,
    which is itself worth knowing and is why it is excluded here on purpose.
    """
    shared = [
        component("tumour_selectivity", 1.0, "basis"),
        component("normal_tissue_safety", 0.8, "basis"),
    ]
    out = []
    for gene, surface in (("AAA", 0.0), ("ZZZ", None)):
        c = TherapeuticTargetCandidate(gene=gene, localization=type_i())
        basis = "basis" if surface is not None else "unknown"
        c.scores = assemble([component("surface_accessibility", surface, basis), *shared], 0.5)
        c.evidence_tier = "observed_alteration"
        c.mechanism_reach = REACH_NONE
        out.append(c)
    return out


def test_a_provisional_only_mechanism_is_still_not_an_established_one():
    """The registered branch of `mechanism_reach`, reached on purpose.

    The fixture above cannot reach it: with no provider data every mechanism
    scores a compatibility of 0.0, so `best_provisional_mechanism` is None and
    `mechanism_reach` answers from the refused-everything branch instead. A
    guard on a branch no fixture reaches is a guard that cannot fail, so this
    one supplies one positive factor beside the unanswered surface requirement
    and lands the candidate in the provisional branch itself.
    """
    c = TherapeuticTargetCandidate(gene="ZZZ", localization=type_i(), origins=[fusion("")])
    inputs = mech.Inputs(
        values={"surface_accessibility": None, "tumour_selectivity": 0.9, "tumour_expression": 0.8},
        bases={"surface_accessibility": "the record does not establish the product's extent"},
    )
    c.therapeutic_mechanisms = [mech.evaluate(c, mech.MECHANISMS["adcc"], inputs)]
    nearest = c.best_provisional_mechanism
    assert nearest is not None and nearest.compatibility > 0.0, (
        "this fixture exists to reach the provisional branch; a compatibility of zero means it "
        "did not and the guard below would be unreachable"
    )
    assert nearest.provisional_requirements == ["surface_accessible"]
    assert c.best_mechanism is None
    reach, reason = mechanism_reach(c)
    assert reach == REACH_NONE, (
        "a mechanism whose hard requirement is merely unanswered has not been shown to apply, so "
        "it may not be classed as reaching the target"
    )
    assert "a question and not an option" in reason


def test_inside_one_cell_the_known_absent_surface_sorts_below_the_unestablished_one():
    """The measured asymmetry, with both absolute levels stated.

    Same tier class, same gate class, every other component identical. The gene
    names are chosen so the fallback would put the known-absent candidate first:
    `AAA` sorts before `ZZZ`. It does not come first, so the mean decided, and
    what the mean did was keep a measured zero in the average and drop an
    unmeasured dimension out of it.
    """
    absent, unestablished = cell_pair()
    assert _rank_key(absent)[:2] == _rank_key(unestablished)[:2], "the pair must share the cell"
    order = sorted([absent, unestablished], key=_rank_key)
    assert [c.gene for c in order] == ["ZZZ", "AAA"], (
        "the alphabetical fallback would have ranked the known-absent candidate first, so this "
        "order is the weighted mean's and nothing else's"
    )
    low, high = absent.scores.overall, unestablished.scores.overall
    assert (low, high) == (0.658, 0.893), (
        f"on the published 0.0-to-1.0 overall scale the known-absent candidate scores {low} and "
        f"the unestablished one {high}, with no cap firing on either"
    )
    assert absent.scores.adjustments == [] and unestablished.scores.adjustments == [], (
        "no safety cap may fire here, or the gap being measured is a cap and not the mean"
    )
    assert absent.scores.coverage > unestablished.scores.coverage, (
        "the known-absent candidate is the one with MORE of the declared dimensions covered, which "
        "is the shape of the asymmetry: being measured is what costs it the places"
    )


def test_assemble_keeps_a_measured_zero_in_the_mean_and_drops_an_unknown_out_of_it():
    """The arithmetic underneath, on explicit numbers and both absolute levels.

    Surface at weight 1.0, selectivity at 1.3 held at 1.0, safety at 1.5 held at
    0.8, so no cap fires. Known zero: (1.0*0.0 + 1.3*1.0 + 1.5*0.8) / 3.8 =
    2.5/3.8 = 0.658. Unestablished: 2.5/2.8 = 0.893. The difference is 0.235 of
    the 0.0-to-1.0 scale and 35.7% of the lower absolute 0.658.
    """
    weights = (
        WEIGHTS["surface_accessibility"],
        WEIGHTS["tumour_selectivity"],
        WEIGHTS["normal_tissue_safety"],
    )
    assert weights == (1.0, 1.3, 1.5)
    absent, unestablished = cell_pair()
    assert absent.scores.overall == 0.658, (
        f"a measured zero at weight 1.0 inside a 3.8 denominator gives 0.658, not {absent.scores.overall}"
    )
    assert unestablished.scores.overall == 0.893, (
        f"dropping that weight leaves a 2.8 denominator and 0.893, not {unestablished.scores.overall}"
    )
    assert round(unestablished.scores.overall - absent.scores.overall, 3) == 0.235, (
        "the dimension that was measured and found absent is the one that lowers the mean; the "
        "dimension that was never established lowers nothing"
    )
    assert absent.scores.coverage > unestablished.scores.coverage, (
        "the lower-scoring candidate is the better-covered one"
    )


def test_the_four_doors_into_the_two_surface_classes_stay_where_they_are():
    """Which inputs answer 0.0 and which answer None, so a count can be read.

    Two doors into the known-absent class and two into the unestablished one.
    This guard exists because the measurement counts candidates by which class
    their surface fell into, and a count is only as good as the classification.
    """
    inside = Localisation(compartments={"cytoplasm": 0.9}, plasma_membrane=False)
    no_evidence = Localisation(compartments={"cytoplasm": 0.9})
    loc = type_i()
    assert surface_accessibility(loc, ectodomain_lost=True)[0] == 0.0
    assert surface_accessibility(inside)[0] == 0.0
    assert surface_accessibility(loc, ectodomain_uncertain=True)[0] is None
    assert surface_accessibility(loc, product_extent_unknown=True)[0] is None
    assert surface_accessibility(no_evidence)[0] is None, (
        "a protein with no membrane evidence at all answers unestablished, which is the door the "
        "benchmark's own candidates never walked through"
    )
    assert surface_accessibility(loc)[0] == 1.0, "the full-length annotation still reads at its own value"


def test_the_fusion_only_predicate_is_what_moves_a_candidate_between_the_classes():
    """The axis each fixture above is varied on, asserted as an axis.

    `ectodomain_lost` and `product_extent_unknown` must never hold together,
    and the orientation record is the only thing that decides which holds.
    """
    recorded = TherapeuticTargetCandidate(gene="AAA", localization=type_i(), origins=[fusion("3'")])
    unrecorded = TherapeuticTargetCandidate(gene="ZZZ", localization=type_i(), origins=[fusion("")])
    assert mech.ectodomain_lost(recorded) and not mech.product_extent_unknown(recorded)
    assert mech.product_extent_unknown(unrecorded) and not mech.ectodomain_lost(unrecorded), (
        "with no orientation on the record the product's extent is the unknown, and the known-loss "
        "clause must not fire on a record that says nothing"
    )
