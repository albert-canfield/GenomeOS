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
    """GenomeOS's distinctive claim: a nuclear or cytoplasmic mutation is never 'no target'."""
    missense = [r for r in rows() if r["expected"] == "out_of_scope_expected"]
    assert all(r.get("peptide_route") for r in missense), (
        "an intracellular missense driver must still offer the peptide/HLA route"
    )


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
    mutated = next(r for r in rows() if r["gene"] == "ERBB2" and r["driver_call"] == "point mutation")
    assert mutated["rank"] == 1 and amp["rank"] > mutated["rank"], (
        "the finding this case exists to record: the same target, reached by the measurement the "
        "approved therapy is actually prescribed on, ranks below where a coding change puts it"
    )


#: Targets outranked by candidates carrying no alteration in the tumour, counted
#: today. It may fall, never rise. ERBB2 reached by its amplification sits behind
#: three such hypotheses and CD19 behind one; the four intracellular cases are not
#: expected to head their lists at all. The cause is one thing: for a surface
#: target the score reads the gene's curated annotation and not the alteration, so
#: the amplified gene and the same gene as a guess scored identically at 0.494.
BURIED_SURFACE_TARGETS = 2


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


def test_the_benchmark_states_what_it_does_not_cover():
    note = RESULT["note"].lower()
    assert "small molecule" in note
    assert "structural variant" in note, "the route no scored case turns on has to be named as uncovered"
    assert "copy-number call" in note, "the route a scored case now turns on has to be named as covered"
    assert "structural variant" not in RESULT["cases_by_driver_call"], (
        "no scored case is driven by a structural variant yet; when one is, say so here"
    )
    assert RESULT["cases"] == len(rows()) >= 8
