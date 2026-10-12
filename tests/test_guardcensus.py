# SPDX-License-Identifier: AGPL-3.0-or-later
"""guardcensus: the rule, pinned before anything is counted."""

from __future__ import annotations

from genomeos.attribution import guardcensus as gcn


def test_condition_4_is_measurement_and_the_first_three_are_description() -> None:
    """The whole point of the rule, and the thing globcensus stage 2 lacked."""
    assert "A VIOLATION WAS ACTUALLY PLANTED AND THE GUARD" in gcn.CANNOT_FIRE_TEST
    assert "(1) to (3) are description; (4) is measurement" in gcn.CANNOT_FIRE_TEST
    assert "was not FOUND, never that none EXISTS" in gcn.CANNOT_FIRE_TEST


def test_the_rule_names_where_it_came_from() -> None:
    assert "globcensus stage 2 met its own three conditions" in gcn.WHY_CONDITION_4
    assert "labelling artefact" in gcn.WHY_CONDITION_4
    assert "open the thing" in gcn.WHY_CONDITION_4


def test_an_unplantable_violation_is_NOT_ESTABLISHED_and_never_cannot_fire() -> None:
    """The branch that would otherwise repeat stage 2's error."""
    assert "NOT_ESTABLISHED" in gcn.IF_NO_VIOLATION_CAN_BE_PLANTED
    assert "a gap in this lane, not a verdict about the guard" in gcn.IF_NO_VIOLATION_CAN_BE_PLANTED


def test_the_generating_mechanism_is_named_as_a_mechanism() -> None:
    assert "COMMITTING A RESULT DISARMS THE SKIP GUARD KEYED ON ITS ABSENCE" in gcn.GENERATING_MECHANISM
    assert "4f2652c" in gcn.GENERATING_MECHANISM, "the prevention must be citable"
    assert "Eleven guards" in gcn.GENERATING_MECHANISM


def test_all_five_cases_are_present_and_each_says_how_it_surfaced() -> None:
    assert len(gcn.CASES) == 5, len(gcn.CASES)
    for c in gcn.CASES:
        assert c["case"] and c["kind"] and c["what"] and c["surfaced_by"], c["case"]
    kinds = {c["kind"] for c in gcn.CASES}
    assert "could never PASS" in kinds, "the always-fires mirror image must be in scope"
    assert any("UNDIAGNOSED" in c["case"] for c in gcn.CASES), "my own instrument must be a case"


def test_the_census_refuses_to_repair_and_says_why() -> None:
    assert "A LIST AND A COUNT, no repairs" in gcn.SCOPE
    assert "doing it here would hide how many were unexercised" in gcn.SCOPE


def test_the_limits_are_fixed_before_the_count() -> None:
    for phrase in (
        "not that an unexercised guard is wrong",
        "NOT_ESTABLISHED will be the",
        "a plant can be",
    ):
        assert phrase in gcn.DOES_NOT_CLAIM, phrase
