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


# ------------------------------------------------------------------------ cases added after the rule


def test_the_ordering_mechanism_is_named_with_its_measured_instance() -> None:
    for phrase in ("036eec6", "40b0a8d", "3 failed", "7,258 passed", "passes 57 at HEAD"):
        assert phrase in gcn.ORDERING_MECHANISM, phrase
    assert "IN THE SAME COMMIT as that record, or after it, never before" in gcn.ORDERING_MECHANISM


def test_the_ordering_case_explains_why_it_is_in_a_cannot_fire_census() -> None:
    """It is the one case where nothing is broken, so its membership needs a reason."""
    assert "nothing in it is broken" in gcn.WHY_ORDERING_IS_THIS_FAMILY
    assert "two instruments that are each right and jointly wrong" in gcn.WHY_ORDERING_IS_THIS_FAMILY


def test_the_extras_lock_blindness_is_quoted_from_the_source() -> None:
    """The claim rests on one line, so the line is checked against the file."""
    from pathlib import Path

    src = (Path(__file__).resolve().parent / "extras_lock.py").read_text()
    assert "names = [node.module] if node.level == 0 and node.module else []" in src
    case = [c for c in gcn.CASES_ADDED if "extras_lock" in c["case"]][0]
    assert "names = [node.module] if" in case["what"]
    assert "INVISIBLE" in case["what"]


def test_the_requires_extra_cost_is_counted_as_a_cost_not_a_fault() -> None:
    case = [c for c in gcn.CASES_ADDED if "requires_extra" in c["case"]][0]
    assert "36 real tests" in case["what"]
    assert "DELIBERATE and its reason is good" in case["note"]
    assert "not a defect in it" in case["note"]


def test_the_inherited_fourth_condition_is_recorded_before_counting() -> None:
    assert "CARRY DATA the result's own scope would have used" in gcn.FOURTH_CONDITION_INHERITED
    assert "globcensus stage 2 was refuted" in gcn.FOURTH_CONDITION_INHERITED


# ------------------------------------------------------ the vacuity leg: condition (3), measured


def test_the_instrument_separates_a_tautology_from_the_nan_idiom() -> None:
    """The discriminator is SHAPE, not the pool: both are invariant over non-NaN values."""
    import ast

    taut = ast.parse("total == compared + (total - compared)", mode="eval").body
    idiom = ast.parse("v == v", mode="eval").body
    assert gcn.classify_comparison(taut)[0] == "tautology"
    assert gcn.classify_comparison(idiom)[0] == "nan_idiom"


def test_the_pool_keeps_every_constant_whatever_its_type() -> None:
    """Version 4's defect: bytes fell between the int and str filters and six sites were lost."""
    import ast

    for src in ("magic == b'BAM\\x01'", "k == 's'", "value == -32768", "n == 512 * 1024 * 1024"):
        node = ast.parse(src, mode="eval").body
        # the value that SATISFIES the comparison must be reachable, or the probe measures itself
        assert gcn.classify_comparison(node)[0] == "varies", src


def test_a_trial_that_raised_is_not_a_trial_that_agreed() -> None:
    """Version 3's second defect: the usable floor, or ~99 of 400 trials decided a verdict."""
    import ast

    node = ast.parse("ne - ns == e - s", mode="eval").body
    verdict, _value, usable = gcn.classify_comparison(node, trials=40, min_usable=10_000)
    assert verdict == "not_measured"
    assert usable < 10_000


def test_the_four_wrong_versions_are_kept_with_their_counts() -> None:
    counts = [v["count"] for v in gcn.VACUITY_INSTRUMENT_WAS_WRONG_FOUR_TIMES]
    assert counts == [1600, 20, 7, 1], counts
    assert "measuring its own pool" in gcn.VACUITY_INSTRUMENT_WAS_WRONG_FOUR_TIMES[0]["defect"]
    assert "DELETED THE ONE TRUE POSITIVE" in gcn.VACUITY_INSTRUMENT_WAS_WRONG_FOUR_TIMES[2]["defect"]


def test_the_vacuity_result_is_stated_as_the_negative_it_is() -> None:
    assert "ALREADY IN CASES" in gcn.VACUITY_RESULT
    assert "found nothing new" in gcn.VACUITY_RESULT
    assert "very nearly EMPTY" in gcn.VACUITY_RESULT
    assert len(gcn.NAN_IDIOM_SET_ASIDE) == 16


def test_the_known_positive_is_still_a_tautology_in_the_tree() -> None:
    """The validation that rejected version 3: this site must appear, or the leg is suspect."""
    import ast
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    src = (root / "scripts" / "manifest_rebuild.py").read_text()
    found = [
        n
        for n in ast.walk(ast.parse(src))
        if isinstance(n, ast.Compare) and gcn.classify_comparison(n)[0] == "tautology"
    ]
    assert found, "manifest_rebuild.py no longer contains the known tautology"
