# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The therapeutic pipeline against targets whose answer is already known.

The benchmark itself needs the network and several minutes
(`scripts/therapeutic_benchmark.py`); these tests read its committed result,
so a regression that loses a known target fails in CI within a second.

Three questions are kept apart on purpose. Recovering the target and calling
the route correctly are asserted, because we believe them. Whether the
mechanism ranked first is defensible was *recorded* rather than asserted for as
long as the pipeline got it wrong, because pretending otherwise would have made
the benchmark decorative. The count is pinned: it may fall, never rise, and it
now reads 0. A pin at 0 is only worth having if the thing it counts can still
be counted, so `test_a_preferred_mechanism_is_never_merely_unrefused` checks the
property the last two defects violated rather than trusting the zero.
"""

from __future__ import annotations

import pytest

from genomeos.results import load_result

RESULT = load_result("therapeutic_benchmark")
pytestmark = pytest.mark.skipif(not RESULT, reason="benchmark result not present")

#: Known defects, counted today. Lower this when one is fixed; a rise is a regression.
#: 3 at the start: an agonist antibody offered against an activating driver, and
#: two mechanisms that headed a list on a requirement nobody had answered —
#: blocking_antibody for BRAF and adcp for PIK3CA, both at 0.25 against a
#: cytoplasmic protein whose membrane compartment is curated at 0.45, under the
#: 0.6 reachability threshold. The agonist went on 2026-09-11, the other two on
#: 2026-09-15 when a preferred mechanism was required to be established and not
#: merely unrefused. Both are still scored and still listed, now as the nearest
#: provisional mechanism with `surface_accessible` named as the open question.
KNOWN_MECHANISM_DEFECTS = 0  # was 2; the two unanswered-requirement cases closed 2026-09-15


def rows():
    return RESULT["rows"]


def test_every_known_target_is_recovered():
    lost = [r["gene"] for r in rows() if not r["recovered"]]
    assert not lost, f"the pipeline lost known targets: {lost}"


def test_the_approved_antibody_targets_are_found_as_surface_targets():
    """EGFR, ERBB2 and CD19 have approved antibodies. If these are not surface targets, nothing is."""
    surface = [r for r in rows() if r["expected"] == "surface"]
    assert len(surface) >= 3
    for r in surface:
        assert r["pass"], f"{r['gene']} has an approved antibody but was not called a surface target"
        assert (r["surface_accessibility"] or 0) > 0.5


def test_no_surface_route_is_claimed_for_an_intracellular_target():
    """The four small-molecule cases pass by *not* claiming a route that does not exist."""
    out = [r for r in rows() if r["expected"] == "out_of_scope_expected"]
    assert len(out) >= 4
    for r in out:
        assert r["pass"], f"{r['gene']}: {r['verdict']}"
        assert (r["surface_accessibility"] or 0) <= 0.5


def test_an_intracellular_driver_still_gets_the_peptide_route():
    """GenomeOS's distinctive claim: a nuclear or cytoplasmic mutation is never 'no target'.

    Scoped to point mutations, which is what it always meant and could leave
    implicit while every case was one. A peptide comes from a changed sequence,
    and for a fusion the changed sequence is the junction, which GenomeOS does
    not reconstruct: the EML4-ALK case therefore offers no peptide route and
    must not be read as this claim failing.
    """
    missense = [
        r for r in rows() if r["expected"] == "out_of_scope_expected" and r["driver_call"] == "point mutation"
    ]
    assert len(missense) >= 4
    assert all(r.get("peptide_route") for r in missense), (
        "an intracellular missense driver must still offer the peptide/HLA route"
    )


def test_a_fusion_is_judged_on_the_product_and_not_the_curated_gene():
    """The case where every database says surface and the tumour says otherwise.

    ALK is curated as a single-pass receptor with a signal peptide and a
    1,020-residue ectodomain. As the fusion's 3' partner it keeps none of it,
    and the approved drugs are small molecules acting inside the cell. The
    benchmark passes this case by *not* claiming a surface route — the same way
    the four intracellular missense cases pass — and it could not have done so
    before the orientation was read.
    """
    alk = next(r for r in rows() if r["driver_call"] == "structural variant")
    assert alk["gene"] == "ALK"
    assert alk["recovered"] and alk["pass"]
    assert alk["target_class"] == "intracellular_only"
    assert alk["surface_accessibility"] == 0.0, "scored from the product the fusion makes"
    assert not alk["best_mechanism"], "no antibody-like route exists against a cytoplasmic kinase"


def test_the_known_mechanism_defects_do_not_grow():
    flagged = [r["gene"] for r in rows() if not r.get("mechanism_sane", True)]
    assert len(flagged) <= KNOWN_MECHANISM_DEFECTS, (
        f"the top-ranked mechanism became indefensible for more targets: {flagged}"
    )
    for r in rows():
        assert "mechanism_sane" in r, "the benchmark stopped scoring the mechanism question"


def test_a_preferred_mechanism_is_never_merely_unrefused():
    """The property the two closed defects violated, checked rather than counted.

    A gate has three answers and the pipeline read two of them. Requiring the
    count of defects to be zero only helps if the pipeline can still report one,
    so this asserts the rule directly: a mechanism that heads a list has had its
    hard requirements answered. And a row with no preferred mechanism has to say
    which question is in the way, because the failure mode of the fix is a
    pipeline that goes quiet instead of one that says "not established".
    """
    for r in rows():
        if not r["recovered"]:
            continue
        assert "best_mechanism_established" in r, "the benchmark stopped recording the distinction"
        if r["best_mechanism"]:
            assert r["best_mechanism_established"], (
                f"{r['gene']}: {r['best_mechanism']} heads the list with "
                f"{r['requirements_unanswered']} unanswered"
            )
        elif r["nearest_provisional_mechanism"]:
            assert r["requirements_unanswered"], (
                f"{r['gene']}: a provisional mechanism must name the requirement that is open"
            )


def test_a_target_is_reached_by_something_other_than_a_point_mutation():
    """Six point mutations was a narrower benchmark than it read as.

    Copy number, structural variants and expression all reach the candidate
    list, and until CD19 was added no scored case turned on any of them, so the
    end-to-end evidence for three of the four routes was zero. CD19 carries no
    alteration in any tumour and is the target of four approved therapies, so
    the only measurement that can reach it is the patient's own RNA.
    """
    calls = RESULT["cases_by_driver_call"]
    assert sum(calls.values()) == RESULT["cases"]
    assert set(calls) - {"point mutation"}, "every case is still driven by a point mutation"
    cd19 = next(r for r in rows() if r["gene"] == "CD19")
    assert cd19["driver_call"] == "expression"
    assert cd19["recovered"] and cd19["pass"]
    assert cd19["best_mechanism_established"], "the approved CD19 drugs are antibody-like"


def test_the_copy_number_route_is_scored_and_not_merely_unit_tested():
    """The same target twice, so that the route is the only thing that differs.

    ERBB2 is in this benchmark under a point mutation and again under nothing
    but its amplification. Varying one thing is the whole reason the second
    case exists: CD19 changes the target and the route together, so a failure
    there cannot be attributed. Amplification is also the clinically honest
    call, since it is what trastuzumab is prescribed on.
    """
    amp = next(r for r in rows() if r["driver_call"] == "copy number")
    assert amp["gene"] == "ERBB2"
    assert amp["recovered"] and amp["pass"]
    assert amp["target_class"] == "direct_surface", "reached as itself, not as a pathway hypothesis"
    assert amp["best_mechanism_established"], "an approved antibody must not be provisional here"
    assert amp["evidence_tier"] == "observed_alteration"
    assert amp["alteration_evidence"] == 1.0, "twelve copies are an observed alteration"
    mutated = next(r for r in rows() if r["gene"] == "ERBB2" and r["driver_call"] == "point mutation")
    assert mutated["rank"] == 1
    assert not amp["outranked_by_hypotheses"], (
        "the finding this case existed to record, now closed: nothing with no measurement in this "
        "tumour outranks the target the approved therapy is prescribed on"
    )
    assert amp["rank"] <= AMPLIFIED_TARGET_RANK, (
        f"the amplified target fell below the pinned rank {AMPLIFIED_TARGET_RANK}"
    )
    # Recorded, not repaired: both routes to the same gene still produce the same
    # number, because the tier orders the candidates and is published per
    # candidate rather than averaged into the mean. What the amplification buys
    # is the place, which is what the clinic reads. The assertion allows the
    # number to rise and not to fall.
    assert amp["score"] >= mutated["score"]


#: Where the amplified ERBB2 case ranks: the case the evidence tier was built
#: for. It may not rise. It was 4, behind KDR, EGFR and PDGFRB, none of them
#: altered in that tumour and all three named only for neighbouring the mutated
#: PIK3CA. It is 1 since the evidence tier of 2026-09-27, with every score in the
#: benchmark unchanged: the tier orders the candidates and is published per
#: candidate, and it is deliberately kept out of the weighted mean. Averaging it
#: in was implemented first and measured: it put PIK3CA above ERBB2 in this very
#: tumour, because the mean is taken over the dimensions that were available and
#: so lifts the candidate with fewer of them further, and because ERBB2's raw
#: score is clipped by the poor-normal-tissue-safety cap that PIK3CA's mediocre
#: safety walks past.
AMPLIFIED_TARGET_RANK = 1  # was 4

#: Targets outranked by candidates with nothing measured about them in this
#: patient, counted today. It may fall, never rise. It read 2 while the score had
#: no dimension for the alteration — the amplified gene and the same gene as a
#: guess scored identically at 0.494 — and it reads 0 since the evidence tier,
#: with the metric also corrected: it counted "candidates with no origins
#: record", under which CD19's own route read as burial.
BURIED_SURFACE_TARGETS = 0  # was 2; the evidence tier, 2026-09-27


def test_a_recovered_target_is_not_quietly_buried_under_hypotheses():
    """Recovery was never the hard question, and three scored questions miss this one.

    A pathway-induced candidate has no origin: it is named for neighbouring
    something altered, and needs expression evidence before it means anything.
    When such a candidate outranks a target whose approved antibody exists, the
    pipeline has recovered the target without preferring it, and every verdict
    above still reads as a pass. The count is pinned rather than asserted to
    zero, because lowering it is a scoring change and this test exists to keep
    it visible until someone makes that change deliberately.
    """
    buried = [
        r for r in rows() if r["expected"] == "surface" and r["recovered"] and r["outranked_by_hypotheses"]
    ]
    assert len(buried) <= BURIED_SURFACE_TARGETS, (
        f"more approved-antibody targets are buried under evidence-free candidates than the "
        f"pinned {BURIED_SURFACE_TARGETS}: {[(r['gene'], r['outranked_by_hypotheses']) for r in buried]}"
    )
    for r in rows():
        assert r["gene"] not in r["outranked_by_hypotheses"], (
            f"{r['gene']} outranks itself: the same gene reached twice, once by its alteration and "
            "once as a hypothesis about a neighbour"
        )


def test_a_target_is_preferred_for_its_alteration_and_not_only_recovered():
    """The question the nine cases could not ask, and the trap inside asking it.

    Recovering a target and preferring it are different results. Until
    2026-09-27 nothing in the score read the alteration: `surface_accessibility`
    reads curated localisation, which describes the gene whether or not the
    tumour touched it, so twelve copies of ERBB2 and a STRING neighbour of a
    mutated gene were worth the same and three unaltered genes outranked the
    target trastuzumab is prescribed on.

    The trap is that the obvious rule — penalise a candidate with no DNA origin
    — demotes CD19, which has no origin either and is the target of four
    approved therapies. So the tier is read from the evidence about the gene in
    this patient, and CD19 sits in its middle tier rather than its bottom one.
    This test asserts both halves: no target is under an unmeasured candidate,
    and the case that carries no alteration of any kind did not pay for it.
    """
    for r in rows():
        if not r["recovered"]:
            continue
        assert r["evidence_tier"] in ("observed_alteration", "patient_measurement"), (
            f"{r['gene']}: an approved target is not a hypothesis about a neighbour"
        )
        assert not r["outranked_by_hypotheses"], (
            f"{r['gene']} is outranked by candidates measured nowhere in this tumour: "
            f"{r['outranked_by_hypotheses']}"
        )
        if r["driver_call"] != "expression":
            assert r["alteration_evidence"] == 1.0
    cd19 = next(r for r in rows() if r["gene"] == "CD19")
    assert cd19["evidence_tier"] == "patient_measurement", (
        "CD19 carries no alteration; it is reached from this patient's own RNA, which is a "
        "measurement of this gene and not an association with another one"
    )
    assert cd19["alteration_evidence"] == 0.6
    assert cd19["rank"] <= 2, "the rule aimed at hypotheses must not demote the expression route"


def test_every_route_into_the_candidate_list_is_scored():
    """The gap this test used to describe is closed, so it now asserts the closure.

    It read "the two routes no scored case turns on have to be named as
    uncovered" while copy number and structural variants were measured by unit
    tests alone. Both are scored now, and the assertion is inverted rather than
    deleted: an end-to-end case per route, so that a route cannot quietly stop
    working while its unit test keeps passing.
    """
    assert set(RESULT["cases_by_driver_call"]) == {
        "point mutation",
        "copy number",
        "structural variant",
        "expression",
    }
    assert RESULT["cases"] == len(rows()) >= 9


def test_the_benchmark_states_what_it_does_not_cover():
    """What remains out of scope is the modality, not a route: no small molecules."""
    note = RESULT["note"].lower()
    assert "small molecule" in note
    out = [r for r in rows() if r["expected"] == "out_of_scope_expected"]
    assert len(out) >= 5, "the small-molecule cases are the ones that state the boundary"
    assert all(r["approved_modality"] == "small_molecule" for r in out)
