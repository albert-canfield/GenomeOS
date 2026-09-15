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
def test_every_claim_is_scored_at_the_controls_too():
    a = RESULT["aggregate"]
    for claim in ("separated_before", "boundary_lost", "new_adjacency"):
        assert a["published"][claim]["n"] >= 1, claim
        assert a["controls"][claim]["n"] >= 1, f"{claim} has no matched control"


@needs_result
def test_naming_the_recipient_is_not_quoted_without_coverage_on_both_sides():
    """The claim that could separate the cases is the one most easily faked by where requests went."""
    n = RESULT["aggregate"]["recipient_named_where_a_deletion_was_scored"]
    if not n["judgeable"]:
        assert n["controls"]["n"] == 0, "not judgeable but the controls do have data"
    assert n["published"]["n"] + n["published"]["not_scored"] == n["published"]["of"]


@needs_result
def test_the_boundary_census_is_asked_before_anything_else():
    """Whether a CTCF-only element sits between donor and recipient at all: at HOXD there was none."""
    for c in RESULT["cases"]:
        b = c["boundary_census"]
        assert b["ctcf_only_between_donor_and_recipient"] >= 0
        assert b["there_is_a_boundary_at_all"] == bool(b["ctcf_only_between_donor_and_recipient"])


@needs_result
def test_a_case_whose_pair_the_node_model_does_not_separate_says_nothing():
    """`separated_before` is the precondition, not a result: no boundary, no experiment."""
    for c in RESULT["cases"]:
        claims = c["published"]["claims"]
        if not claims["separated_before"]:
            assert not claims["new_adjacency"], (
                f"{c['locus']}: a new adjacency where the pair was never separated is a bug"
            )
