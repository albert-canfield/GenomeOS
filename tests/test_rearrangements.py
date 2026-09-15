# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The node model against published rearrangements, free half (docs/LOCI-BENCHMARK.md 11b).

The coordinate transforms and the expectation record are tested offline. The committed result
(`data/results/loci_rearrangements.json`) is read where it exists, so a verdict that changes fails
CI in a second. No test makes a model request; the run does not either.
"""

from __future__ import annotations

import pytest

from genomeos.benchmark.rearrangements import (
    ANALYTIC_NOTE,
    CASES,
    KINDS,
    NOT_YET_RUNNABLE,
    PREDICTION,
    rearranged_length,
    transform,
)
from genomeos.results import load_result

RESULT = load_result("loci_rearrangements")
needs_result = pytest.mark.skipif(not RESULT, reason="rearrangement result not present")


# ------------------------------------------------------------------- the transforms, offline
def test_a_deletion_removes_the_span_and_shifts_what_is_above_it():
    span = (100, 200)
    assert transform(50, span, "deletion") == 50
    assert transform(150, span, "deletion") is None
    assert transform(200, span, "deletion") == 100
    assert transform(1000, span, "deletion") == 900
    assert rearranged_length(1_000, span, "deletion") == 900


def test_an_inversion_mirrors_inside_the_span_and_leaves_the_rest():
    span = (100, 200)
    assert transform(50, span, "inversion") == 50
    assert transform(100, span, "inversion") == 199
    assert transform(199, span, "inversion") == 100
    assert transform(1000, span, "inversion") == 1000
    assert rearranged_length(1_000, span, "inversion") == 1_000


def test_a_duplication_pushes_everything_above_the_span_up():
    span = (100, 200)
    assert transform(50, span, "duplication") == 50
    assert transform(150, span, "duplication") == 150
    assert transform(200, span, "duplication") == 300
    assert rearranged_length(1_000, span, "duplication") == 1_100


def test_neither_endpoint_of_a_case_is_swallowed_by_its_own_deletion():
    """A position inside the interval being removed goes nowhere.

    The first control construction put the donor's transcription start at the span's lower edge, so
    it resolved to None and every count downstream was nonsense. The span must lie strictly between
    the two genes.
    """
    span = (100, 200)
    assert transform(99, span, "deletion") is not None
    assert transform(200, span, "deletion") is not None


# ------------------------------------------------------------------- the expectations, offline
def test_every_case_is_written_down_with_its_sources_and_its_span_provenance():
    for c in CASES:
        assert c.kind in KINDS, c.locus
        assert c.donor and c.recipient and c.donor != c.recipient, c.locus
        assert c.citations and c.answer_from and c.phenotype, c.locus
        assert c.span_source in ("published", "stated"), c.locus
        assert len(c.span_citation) > 80, (
            f"{c.locus}: a span the panel states needs to say on whose authority"
        )


def test_the_cases_that_cannot_run_say_what_each_is_waiting_for():
    assert NOT_YET_RUNNABLE, "dropping a case silently is what the panel exists to prevent"
    for name, why in NOT_YET_RUNNABLE.items():
        assert len(why) > 60, name
        assert "liftover" not in why.lower(), (
            f"{name}: the blocker is that no breakpoint coordinate is published, not a liftover -"
            " saying otherwise sends the next session after a chain file it does not need"
        )


def test_no_case_is_placed_approximately():
    """The rule: if a breakpoint cannot be had cleanly, drop the case rather than guess it.

    Applied to all four, including the EPHA4-to-PAX3 deletion an earlier run placed approximately.
    """
    assert CASES == (), "a case is being run without published breakpoints"
    assert "EPHA4_PAX3_deletion" in NOT_YET_RUNNABLE


def test_the_analytic_note_says_why_two_of_the_three_claims_cannot_discriminate():
    assert "guaranteed by construction" in ANALYTIC_NOTE
    assert "arithmetic, not biology" in ANALYTIC_NOTE


# ------------------------------------------------------------------------- the gate, on the result
@needs_result
def test_the_prediction_is_kept_verbatim_in_the_result():
    """Registered before the first run; it is only worth anything if it cannot be edited afterwards."""
    assert RESULT["prediction_registered_before_the_run"] == PREDICTION
    assert "node model has failed its sharpest test" in PREDICTION


@needs_result
def test_the_node_model_is_never_counted_as_a_derived_hit():
    assert "inferred" in RESULT["evidence"]["node_model"]
    assert "NOT RUN" in RESULT["evidence"]["derived"]


@needs_result
def test_the_result_says_plainly_that_no_case_is_runnable():
    assert RESULT["cases"] == []
    assert "nothing to lift over" in RESULT["no_case_is_runnable"]


@needs_result
def test_a_deletion_of_the_published_size_loses_a_boundary_wherever_it_is_put():
    """The measurement that survives: claims one and two fire everywhere, so they discriminate nothing."""
    b = RESULT["boundary_behaviour"]
    assert b["rates"]["separated_before"]["n"] >= 1
    r = b["rates"]
    assert r["boundary_lost"]["k"] == r["boundary_lost"]["n"], (
        "a size-matched deletion that did not lose a boundary would contradict the analytic note"
    )
    assert r["new_adjacency"]["k"] == r["new_adjacency"]["n"]


@needs_result
def test_naming_the_recipient_is_not_quoted_without_coverage():
    """The one claim that could discriminate, and it cannot be scored until chr2 is swept."""
    b = RESULT["boundary_behaviour"]
    if not b["scored_near_the_recipient"]:
        assert b["rates"]["recipient_named"]["k"] == 0, (
            "a recipient was named with no deletion scored near it: that cannot happen"
        )
