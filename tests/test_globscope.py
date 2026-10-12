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


# ====================================================================== the result, read back
#
# These guards READ data/results/globscope.json. They are in the same file as the rule's guards above
# but were added in their own commit, AFTER the record: OWN_CODE declares this test file, so
# check_staged condition (a) refuses a record written while its test was edited. No pytest.skip is
# keyed on the record -- tracked_paths.must_be_present raises instead, from this version on.

from tests import tracked_paths  # noqa: E402

RESULT = ROOT / "data/results/globscope.json"


def _result() -> dict:
    tracked_paths.must_be_present("data/results/globscope.json")
    return json.loads(RESULT.read_text())


def test_the_classes_sum_to_the_inherited_eighteen():
    d = _result()
    assert sum(d["class_counts"].values()) == globscope.inherited_site_count() == 18
    assert len(d["sites"]) == 18


def test_every_site_carries_exactly_one_class_from_the_registered_set():
    d = _result()
    for s in d["sites"]:
        assert s["class"] in globscope.CLASSES, s["path"]


def test_the_class_of_every_site_is_RECOMPUTED_from_its_own_recorded_facts():
    """The record cannot carry a class its own facts do not produce under the committed rule."""
    d = _result()
    for s in d["sites"]:
        assert s["class"] == globscope.classify_site(
            refusal_in_enclosing_scope=s["refusal_in_enclosing_scope"],
            glob_root_resolved=s["glob_root_resolved"],
            matched_files_on_disk=s["matched_files_on_disk"],
            all_matched_files_tracked=s["all_matched_files_tracked"],
            independent_scope_recorded=s["independent_scope_recorded"],
        ), f"{s['path']}:{s['line']}"


def test_the_counts_agree_with_the_sites_they_count():
    d = _result()
    for name, n in d["class_counts"].items():
        assert sum(1 for s in d["sites"] if s["class"] == name) == n, name


def test_the_four_refusals_are_the_check_registration_functions_and_each_cites_its_raise():
    d = _result()
    refusers = [s for s in d["sites"] if s["class"] == "REFUSES"]
    assert len(refusers) == 4
    assert {s["enclosing_function"] for s in refusers} == {"check_registration"}
    for s in refusers:
        ev = " ".join(s.get("read_evidence", []))
        assert "pinned_sha256" in ev and "SystemExit" in ev, s["path"]
        # the boundary travels with the fact, not in a separate paragraph
        assert "REGISTRATION'S PIN" in ev, s["path"]


def test_the_one_live_site_is_named_and_its_untracked_store_is_cited():
    d = _result()
    live = [s for s in d["sites"] if s["class"] == "CAN_SHRINK_SILENTLY"]
    assert len(live) == 1
    assert live[0]["path"] == "scripts/reader_genome_wide.py"
    assert live[0]["line"] == 88
    ev = " ".join(live[0]["read_evidence"])
    assert ".gitignore:16" in ev and "mutually exclusive" in ev.lower()


def test_no_site_is_CANNOT_SHRINK_and_the_ignored_stores_are_why():
    """The uncomfortable half of the answer: not one glob reads a wholly tracked store."""
    d = _result()
    assert d["class_counts"]["CANNOT_SHRINK"] == 0
    assert all(s["all_matched_files_tracked"] is not True for s in d["sites"])
    evidence = " ".join(" ".join(s.get("read_evidence", [])) for s in d["sites"])
    assert ".gitignore:50" in evidence and "299" in evidence


def test_no_read_fact_names_a_class_and_every_one_cites_a_path_line():
    d = _result()
    for a in d["read_facts_applied"]:
        assert "class" not in a["sets"] and "classification" not in a["sets"]
        assert ":" in a["evidence"]
        assert set(a["sets"]) <= {
            "refusal_in_enclosing_scope",
            "glob_root_resolved",
            "matched_files_on_disk",
            "all_matched_files_tracked",
            "independent_scope_recorded",
        }


def test_a_read_fact_never_DENIES_a_refusal():
    d = _result()
    for a in d["read_facts_applied"]:
        if "refusal_in_enclosing_scope" in a["sets"]:
            assert a["sets"]["refusal_in_enclosing_scope"] is True


# ----------------------------------------------------------------- the six varying groups


def test_the_inherited_six_are_six_and_are_the_ones_globcensus_named():
    d = _result()
    rows = d["varying_groups_the_inherited_six_CORRECTED"]
    assert len(rows) == 6
    inherited = json.loads((ROOT / "data/results/globcensus.json").read_text())
    assert {r["group"] for r in rows} == set(inherited["groups_whose_count_varies"])


def test_every_varying_class_is_RECOMPUTED_from_the_committed_rule():
    d = _result()
    for r in d["varying_groups_the_inherited_six_CORRECTED"]:
        classes = set()
        for x in r["per_result"]:
            assert x["class"] in globscope.VARIATION_CLASSES
            classes.add(x["class"])
        expected = (
            next(iter(classes))
            if len(classes) == 1
            else "MIXED_ACROSS_RESULTS: " + ", ".join(sorted(classes))
        )
        assert r["group_class"] == expected, r["group"]


def test_no_varying_group_is_called_a_SILENT_SHRINK_by_this_rule():
    d = _result()
    for r in d["varying_groups_the_inherited_six_CORRECTED"]:
        for x in r["per_result"]:
            assert x["class"] != "SILENT_SHRINK", (r["group"], x["result"])


def test_every_LEGITIMATE_SCOPE_names_the_field_that_decided_it():
    d = _result()
    for r in d["varying_groups_the_inherited_six_CORRECTED"]:
        for x in r["per_result"]:
            if x["class"] == "LEGITIMATE_SCOPE":
                assert x["commensurable"] is True
                assert "parameters." in str(x["the_field_that_decided_it"])
                assert x["declared"] <= x["scope_count"]


def test_PLANTED_an_incommensurable_scope_is_UNDECIDABLE_and_not_legitimate():
    """The first pass's error, held as a test: 73 files against a 24-chromosome list decides nothing."""
    d = _result()
    sweep = next(
        r
        for r in d["varying_groups_the_inherited_six_CORRECTED"]
        if r["group"].startswith("the all-enhancer sweep")
    )
    assert sweep["group_class"] == "UNDECIDABLE"
    for x in sweep["per_result"]:
        assert x["commensurable"] is False
        assert "no commensurable scope" in x["the_field_that_decided_it"]


def test_the_result_keeps_BOTH_passes_and_says_which_is_the_answer():
    d = _result()
    assert d["varying_groups_the_inherited_six"], "the first pass is kept"
    assert "CORRECTED" in d["which_pass_is_the_answer"]
    assert "flatter" in d["why_a_second_pass"] or "flatter" in d["which_pass_is_the_answer"]


def test_the_refinement_not_taken_is_recorded_as_a_decision():
    d = _result()
    note = d["the_refinement_not_taken"]
    assert "NOT DO IT" in note.upper()
    assert "flattering direction" in note
    assert "0 of 6" in note


# ------------------------------------------------- the carried wording, inside the record itself


def test_the_record_carries_the_inherited_reading_word_for_word():
    d = _result()
    assert d["inherited"]["reading_carried_verbatim"] == globscope.INHERITED_READING
    assert d["inherited"]["limit_carried_verbatim"] == globscope.INHERITED_LIMIT


def test_the_record_carries_shrinkimpacts_refutation_and_does_not_restate_it_as_its_own():
    from genomeos.attribution import globcensus

    d = _result()
    assert d["inherited"]["stage_2_refuted"] == globcensus.STAGE_2_REFUTED
    assert d["inherited"]["stage_2_summary"] == globcensus.STAGE_2_SUMMARY
    assert "realised in 0 of the 6 checked" in d["inherited"]["stage_2_summary"]


def test_the_population_move_is_recorded_rather_than_absorbed():
    d = _result()
    moved = d["the_population_moved_while_this_lane_worked"]
    assert moved["inherited_sites_classified"] == 18
    assert moved["exists_filtered_sites_in_the_tree_now"] == 23
    assert len(moved["new_since_globcensus"]) == 5
    assert len(d["varying_groups_new_since_globcensus"]) == 2
    for named in moved["new_since_globcensus"]:
        assert not any(s["path"] in named and str(s["line"]) in named for s in d["sites"])


def test_the_record_claims_no_money_and_no_model_request():
    d = _result()
    assert d["model_requests"] == 0
    assert d["money_spent"] == "$0.00"
    assert "2026-09-27" in d["stop_list"]


def test_the_manifest_is_complete_and_its_inputs_are_repo_relative():
    d = _result()
    m = d["result_manifest"]
    assert m["complete"] is True
    assert m["code_cleanliness"]["own_code_is_committed"] is True
    for entry in m["inputs"]:
        assert not str(entry["path"]).startswith("/"), entry["path"]
