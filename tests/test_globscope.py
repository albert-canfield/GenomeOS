# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Guards on globscope's RULE, committed with the rule and before any call site is read.

Every test here checks the rule's shape, not a finding: that the classes are mutually exclusive, that
the decision order is the one written down, that each class can actually be reached by a planted set of
facts (a class no fact can reach is a class that will never accuse anything), that UNDECIDABLE cannot
be squeezed out by a residual, that the inherited count is the same object and not a retyped literal,
and that the inherited limit and reading are carried word for word.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genomeos.attribution import globscope

ROOT = Path(__file__).resolve().parent.parent


def test_the_five_classes_are_named_and_each_carries_a_decidable_test():
    assert set(globscope.CLASSES) == {
        "REFUSES",
        "CANNOT_SHRINK",
        "SCOPE_IS_DECLARED_ELSEWHERE",
        "CAN_SHRINK_SILENTLY",
        "UNDECIDABLE_FROM_THE_CODE",
    }
    for name, text in globscope.CLASSES.items():
        assert "DECIDABLE TEST" in text or "DECLARED IN ADVANCE" in text, name


def test_decision_order_covers_every_class_exactly_once():
    assert sorted(globscope.DECISION_ORDER) == sorted(globscope.CLASSES)
    assert len(set(globscope.DECISION_ORDER)) == len(globscope.DECISION_ORDER)


def test_undecidable_is_second_and_the_reason_is_written_down():
    assert globscope.DECISION_ORDER[1] == "UNDECIDABLE_FROM_THE_CODE"
    assert len(globscope.WHY_UNDECIDABLE_IS_SECOND) > 200


# ------------------------------------------------------------- every class is reachable by some fact

_FULL = dict(
    refusal_in_enclosing_scope=False,
    glob_root_resolved=True,
    matched_files_on_disk=3,
    all_matched_files_tracked=False,
    independent_scope_recorded=False,
)


@pytest.mark.parametrize(
    "override,expected",
    [
        ({"refusal_in_enclosing_scope": True}, "REFUSES"),
        ({"glob_root_resolved": False, "all_matched_files_tracked": None}, "UNDECIDABLE_FROM_THE_CODE"),
        ({"matched_files_on_disk": 0, "all_matched_files_tracked": None}, "UNDECIDABLE_FROM_THE_CODE"),
        ({"all_matched_files_tracked": None}, "UNDECIDABLE_FROM_THE_CODE"),
        ({"all_matched_files_tracked": True}, "CANNOT_SHRINK"),
        ({"independent_scope_recorded": True}, "SCOPE_IS_DECLARED_ELSEWHERE"),
        ({}, "CAN_SHRINK_SILENTLY"),
    ],
)
def test_PLANTED_each_class_is_reachable(override, expected):
    assert globscope.classify_site(**{**_FULL, **override}) == expected


def test_PLANTED_a_refusal_wins_over_every_other_fact():
    """REFUSES is first: a site that raises on a short set is not accused of shrinking, whatever its
    store looks like."""
    assert (
        globscope.classify_site(
            **{**_FULL, "refusal_in_enclosing_scope": True, "all_matched_files_tracked": True}
        )
        == "REFUSES"
    )


def test_PLANTED_an_unresolved_glob_does_NOT_fall_into_can_shrink_silently():
    """The ordering decision that matters: an unreadable directory says so rather than accusing."""
    got = globscope.classify_site(**{**_FULL, "glob_root_resolved": False, "all_matched_files_tracked": None})
    assert got == "UNDECIDABLE_FROM_THE_CODE"
    assert got != "CAN_SHRINK_SILENTLY"


def test_PLANTED_a_tracked_directory_is_not_assumed_from_a_resolved_root_alone():
    """CANNOT_SHRINK needs tracked-ness MEASURED true, not merely a root that resolved."""
    assert globscope.classify_site(**_FULL) == "CAN_SHRINK_SILENTLY"
    assert (
        globscope.classify_site(**{**_FULL, "all_matched_files_tracked": None}) == "UNDECIDABLE_FROM_THE_CODE"
    )


# ----------------------------------------------------------------- the varying-group rule

_VAR = dict(
    artefact_states_a_scope=True,
    declared_count=22,
    scope_count=24,
    missing_members_present_in_store_today=True,
    missing_members_named_in_exclusions=False,
    missing_members_carry_eligible_data=True,
)


@pytest.mark.parametrize(
    "override,expected",
    [
        ({}, "SILENT_SHRINK"),
        ({"artefact_states_a_scope": False}, "UNDECIDABLE"),
        ({"scope_count": None}, "UNDECIDABLE"),
        ({"declared_count": 24}, "LEGITIMATE_SCOPE"),
        ({"missing_members_named_in_exclusions": True}, "LEGITIMATE_SCOPE"),
        ({"missing_members_carry_eligible_data": False}, "LEGITIMATE_SCOPE"),
        ({"missing_members_present_in_store_today": False}, "LEGITIMATE_SCOPE"),
        ({"missing_members_carry_eligible_data": None}, "UNDECIDABLE"),
    ],
)
def test_PLANTED_the_variation_rule_reaches_each_of_its_three_classes(override, expected):
    assert globscope.classify_variation(**{**_VAR, **override}) == expected


def test_PLANTED_no_stated_scope_is_UNDECIDABLE_and_never_legitimate():
    """The clause the brief insisted on: an artefact silent about its scope is not 'probably fine'."""
    got = globscope.classify_variation(**{**_VAR, "artefact_states_a_scope": False})
    assert got == "UNDECIDABLE"
    assert got != "LEGITIMATE_SCOPE"


def test_PLANTED_the_fourth_condition_can_block_a_shrink_verdict():
    """lane-shrinkimpact's condition, carried whole: without it the first three conditions alone would
    have called all six of stage 2 a shrink, which was refuted."""
    assert globscope.classify_variation(**_VAR) == "SILENT_SHRINK"
    assert (
        globscope.classify_variation(**{**_VAR, "missing_members_carry_eligible_data": False})
        == "LEGITIMATE_SCOPE"
    )
    assert "is_shrink_v2" in globscope.VARIATION_CLASSES["SILENT_SHRINK"]


# ----------------------------------------------------------------- the pinned inherited count


def test_the_inherited_18_is_READ_from_globcensus_and_not_retyped():
    assert globscope.inherited_site_count() == 18
    d = json.loads((ROOT / "data/results/globcensus.json").read_text())
    assert globscope.inherited_site_count() == d["call_sites"]["exists_filtered"]


def test_PLANTED_the_pinned_count_follows_its_source_rather_than_a_literal(tmp_path):
    """If globcensus.json said 19, this module would say 19. A literal 18 in the source would make the
    two disagree silently, which is the exact failure mode the pin exists to prevent."""
    d = json.loads((ROOT / "data/results/globcensus.json").read_text())
    d["call_sites"]["exists_filtered"] = 19
    (tmp_path / "data/results").mkdir(parents=True)
    (tmp_path / globscope.INHERITED_COUNT_SOURCE).write_text(json.dumps(d))
    assert globscope.inherited_site_count(root=tmp_path) == 19


def test_source_text_carries_no_bare_literal_count_of_the_sites():
    src = (ROOT / "genomeos/attribution/globscope.py").read_text()
    assert "INHERITED_COUNT_FIELD" in src
    assert "d = d[key]" in src


# ----------------------------------------------------------------- the carried wording


def test_the_inherited_limit_is_carried_word_for_word():
    assert "NOT shown safe" in globscope.INHERITED_LIMIT
    assert "CONSTANT at a shrunken value is invisible" in globscope.INHERITED_LIMIT
    assert "82 constant groups" in globscope.INHERITED_LIMIT
    assert "the 13 that all say 312" in globscope.INHERITED_LIMIT


def test_the_inherited_reading_is_carried_and_not_strengthened():
    assert globscope.INHERITED_READING == (
        "the defect is real as a mechanism in 18 places and so far unrealised where it was found"
    )
    for banned in ("proven", "shown safe", "widespread", "systemic"):
        assert banned not in globscope.INHERITED_READING


def test_every_negative_carries_its_own_boundary():
    for part in ("18 sites", "AST scan", "6 varying", "nothing wider"):
        assert part in globscope.NEGATIVE_BOUNDARY


def test_the_outcomes_were_written_before_the_scan_and_cover_all_three_cases():
    assert set(globscope.PREWRITTEN_OUTCOMES) == {
        "all_18_cannot_shrink",
        "some_can_shrink_silently",
        "most_undecidable",
        "the_six_varying",
    }
    assert "COMPLETE ANSWER" in globscope.PREWRITTEN_OUTCOMES["most_undecidable"].upper()


def test_does_not_claim_keeps_the_constant_groups_unclaimed():
    assert "NOT SHOWN SAFE" in globscope.DOES_NOT_CLAIM.upper()
    assert "0 of 6" in globscope.DOES_NOT_CLAIM


def test_the_stop_list_was_read_and_nothing_on_it_is_touched():
    assert "2026-09-27" in globscope.STOP_LIST_READ
    assert "0 model requests" in globscope.STOP_LIST_READ
