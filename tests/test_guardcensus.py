# SPDX-License-Identifier: AGPL-3.0-or-later
"""guardcensus: the rule, pinned before anything is counted."""

from __future__ import annotations

import pytest

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


@pytest.mark.xfail(
    strict=False,
    reason="SUPERSEDED, and this lane's own defect: `ne - ns == e - s` genuinely VARIES, and the"
    " varies branch returns as soon as two verdicts are seen -- before the usable floor is ever"
    " consulted. So this assertion was wrong when 271912f committed it, and FLAKY rather than"
    " simply red, because classify_comparison draws random values: it passed once and has mostly"
    " failed since. The correct test is"
    " test_the_usable_floor_is_reached_when_every_evaluation_raises below. Left in place with"
    " xfail(strict=False) rather than deleted, because removing it needs --force inside the undo"
    " window and the rule is: never --force; a red you caused yourself is still a report, not an"
    " exception. strict=False and not strict=True precisely because it is flaky -- a strict xfail"
    " would then fail on the runs where it happens to pass.",
)
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


# --------------------------------- condition (4) satisfied: the extras_lock closure-walk blindness


def test_the_blindness_is_in_the_first_party_leg_and_the_census_says_so() -> None:
    """The third-party leg is CORRECT, and the census's correction must match the source.

    `Import(name.split(".")[0], ...)` is the right answer to "which distribution does this file
    need", so `from scipy import stats` recording scipy is not a defect. Only the closure walk
    loses the edge.
    """
    from pathlib import Path

    src = (Path(__file__).resolve().parent / "extras_lock.py").read_text()
    assert 'imports.append(Import(name.split(".")[0], path, node.lineno))' in src
    assert "first_party.add(name)" in src
    assert "THE THIRD-PARTY LEG IS CORRECT" in gcn.BLINDNESS_CORRECTION
    assert "first-party closure walk" in gcn.BLINDNESS_CORRECTION


def test_the_overstated_case_is_left_standing_not_edited() -> None:
    """As globcensus's refuted headline was: the record shows what was claimed and when."""
    case = [c for c in gcn.CASES_ADDED if "extras_lock" in c["case"]][0]
    assert "INVISIBLE" in case["what"], "the original wording must remain, for the record"
    assert "too broad and the entry is left standing" in gcn.BLINDNESS_CORRECTION


@pytest.mark.requires_extra("predict")
def test_the_closure_walk_still_drops_the_submodule_edge() -> None:
    """The mechanism, on a synthetic file, so no live result record is read.

    WHEN extras_lock IS FIXED, THIS IS THE TEST TO UPDATE, in the same commit as the fix: invert
    the assertion and amend BLINDNESS_IS_REALISED to record that it was repaired and when. A
    census must not become a barrier to its own repair, so this pins the measured state and says
    how to move it rather than forbidding the move.
    """
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import extras_lock

    tmp = Path(__file__).resolve().parent / "extras_lock.py"  # any real first-party file
    assert extras_lock._first_party_file("genomeos.attribution.pbrun") is not None, (
        "the fixture for this test is gone: genomeos/attribution/pbrun.py no longer exists"
    )
    assert tmp.exists()
    # the shape the walk cannot follow
    resolved = extras_lock._first_party_file("genomeos.attribution")
    assert resolved is not None and resolved.name == "__init__.py", resolved


def test_condition_4_is_recorded_as_SATISFIED_here_and_not_merely_described() -> None:
    assert "AS COMMITTED: 0 violations" in gcn.BLINDNESS_IS_REALISED
    assert "WITH THE EDGE RECORDED: 1 violation" in gcn.BLINDNESS_IS_REALISED
    assert "condition (4) is SATISFIED" in gcn.BLINDNESS_IS_REALISED
    assert "established CANNOT-FIRE" in gcn.BLINDNESS_IS_REALISED


def test_the_scale_is_stated_with_the_caveat_that_keeps_it_honest() -> None:
    assert "1,942" in gcn.BLINDNESS_SCALE
    assert "NOT 1,942" in gcn.BLINDNESS_SCALE
    assert "exactly ONE" in gcn.BLINDNESS_SCALE


def test_the_instance_records_why_it_is_not_red_yet_and_why_the_local_gate_misses_it() -> None:
    inst = gcn.BLINDNESS_INSTANCE
    assert inst["test"] == "tests/test_pbrun.py"
    assert "h5py" in inst["through"]
    assert "zero requires_extra calls" in inst["declared"]
    assert "on no remote branch" in inst["why_not_red_yet"]
    assert "installed in the shared .venv" in inst["why_the_local_gate_misses_it"]
    assert "openpyxl" in inst["same_class_as"]


def test_no_repair_is_made_by_this_lane() -> None:
    assert "reported and not fixed" in gcn.BLINDNESS_NO_REPAIR
    assert "requires_extra('predict')" in gcn.BLINDNESS_NO_REPAIR


# ----------------------------------- the stopgap measured in CI's condition, and procedure records


def test_the_instance_is_superseded_by_the_measurement_not_edited() -> None:
    """f4a9faa made 'test-bare fails at COLLECTION' false; the old wording stays for the record."""
    assert "fails at COLLECTION" in gcn.BLINDNESS_INSTANCE["why_not_red_yet"]
    # `... or True` stood here until ruff's SIM222 refused it: a guard that CANNOT FAIL, written
    # into the census of guards that cannot fail. Recorded in CASES_SELF rather than hidden.
    assert "SO THE STOPGAP WORKS" in gcn.INSTANCE_SUPERSEDED
    assert "7,500 tests" in gcn.INSTANCE_SUPERSEDED
    assert "ZERO" in gcn.INSTANCE_SUPERSEDED


def test_absence_was_simulated_by_veiling_finders_not_by_monkeypatching_a_return() -> None:
    """The distinction is the whole value of the measurement."""
    assert "WRAPS every finder" in gcn.HOW_ABSENCE_WAS_SIMULATED
    assert "was NOT monkeypatched" in gcn.HOW_ABSENCE_WAS_SIMULATED
    assert "No" in gcn.HOW_ABSENCE_WAS_SIMULATED and "uv sync" in gcn.HOW_ABSENCE_WAS_SIMULATED


def test_the_stopgaps_inert_condition_is_recorded_as_a_caveat_not_a_fault() -> None:
    assert "CANNOT ARISE IN CI" in gcn.STOPGAP_INERT_CONDITION
    assert "anndata requires h5py" in gcn.STOPGAP_INERT_CONDITION
    assert "not mistaken for a fault" in gcn.STOPGAP_INERT_CONDITION


def test_the_announcement_case_names_all_three_measurements() -> None:
    case = gcn.ANNOUNCEMENT_CANNOT_BE_HEARD
    assert "without -q the line is printed" in case["measured"]
    assert "-q -rs --no-fold-skipped" in case["measured"]
    assert "never a silent absence" in case["against_its_own_words"]
    assert "the instrument is RIGHT and its reading never reaches anyone" in case["family"]


def test_the_cli_flags_that_suppress_the_header_are_still_the_ones_ci_uses() -> None:
    """If CI stops passing -q, this case is void and must be withdrawn, not left standing."""
    from pathlib import Path

    ci = (Path(__file__).resolve().parents[1] / ".github" / "workflows" / "ci.yml").read_text()
    assert "pytest -q -rs --no-fold-skipped" in ci, (
        "CI's flags have changed; ANNOUNCEMENT_CANNOT_BE_HEARD must be re-measured"
    )


def test_the_collect_ignore_confound_is_recorded_for_the_next_lane() -> None:
    assert "governs DIRECTORY collection" in gcn.COLLECT_IGNORE_AND_EXPLICIT_PATHS
    assert "may read it as the stopgap failing" in gcn.COLLECT_IGNORE_AND_EXPLICIT_PATHS


def test_the_no_force_rule_is_carried_verbatim() -> None:
    assert gcn.NO_FORCE_RULE == (
        "never `--force`; a red you caused yourself is still a report, not an exception"
    )
    assert "every bad force also says" in gcn.WHY_THE_WORDING_WAS_TIGHTENED


def test_the_procedure_deviation_lists_every_fact_that_made_it_acceptable() -> None:
    dev = gcn.PROCEDURE_DEVIATIONS[0]
    assert dev["commit"] == "b2bf090"
    assert len(dev["facts"]) == 6, dev["facts"]
    assert "2 added, 1 removed, one file" in dev["facts"][2]
    assert "procedure, not instrument" in dev["class"]


def test_the_census_records_its_own_cannot_fire_guard() -> None:
    """Found by ruff SIM222 in the commit that recorded the stopgap; repaired, and recorded."""
    case = gcn.CASES_SELF[0]
    assert "or True" in case["what"]
    assert "ruff SIM222" in case["found_by"]
    assert "condition (4) is not pedantry" in case["why_recorded"]
    assert case["fixed"].startswith("yes")


def test_the_usable_floor_is_reached_when_every_evaluation_raises() -> None:
    """The corrected form of the test above: an expression where NO trial can succeed.

    `n % 0 == 1` raises ZeroDivisionError on every substitution, so usable stays 0, the varies
    branch is never reached, and the floor decides -- which is what version 3's second defect was
    about.
    """
    import ast

    node = ast.parse("n % 0 == 1", mode="eval").body
    verdict, _value, usable = gcn.classify_comparison(node, trials=40)
    assert verdict == "not_measured", verdict
    assert usable == 0, usable


def test_the_census_records_its_own_flaky_guard_and_how_it_was_caught() -> None:
    """The undo guard caught it, not a test and not review."""
    case = gcn.CASES_SELF_FLAKY
    assert "3 xfailed, 2 xpassed" in case["worse_than_wrong"]
    assert "intermittent" in case["worse_than_wrong"]
    assert "UNDO GUARD refusing the commit" in case["how_it_was_caught"]
    assert "ZERO removed" in case["why_not_deleted"]
    assert "guard that fires at random" in case["why_strict_false"]


def test_the_reading_of_the_count_is_registered_before_the_count() -> None:
    """Fixed before fired.json existed, so it cannot be fitted to the numbers."""
    assert "1,736 raise sites" in gcn.COUNT_POPULATION
    assert "1,550 guards" in gcn.COUNT_POPULATION
    assert "not doing the work" in gcn.COUNT_POPULATION


def test_never_fired_is_not_claimed_to_be_a_defect() -> None:
    assert "NOT a defect and NOT a cannot-fire" in gcn.COUNT_READING
    assert "only condition (4) can tell the two apart" in gcn.COUNT_READING
    assert "EXPECTED and nearly uninformative" in gcn.COUNT_READING


def test_the_population_warning_cites_the_lanes_own_mistake() -> None:
    assert "CASES_SELF_RED is the reason" in gcn.COUNT_POPULATION_WARNING
    assert "two readers before the guard disproved it" in gcn.COUNT_POPULATION_WARNING


# ------------------------------------------------- the abandoned suite leg, and the static count


def test_the_abandoned_leg_states_its_cost_and_that_nothing_was_salvaged() -> None:
    assert "44% in 160 minutes" in gcn.SUITE_LEG_ABANDONED
    assert "MACHINE-WIDE suite lock" in gcn.SUITE_LEG_ABANDONED
    assert "atexit does not run on SIGTERM" in gcn.SUITE_LEG_ABANDONED


def test_the_cause_is_named_as_one_the_callback_could_not_fix() -> None:
    """The dedupe was real and irrelevant; saying so is the point."""
    assert "deoptimises the interpreter" in gcn.WHY_IT_WAS_SLOW
    assert "AT THE SAME 15x" in gcn.WHY_IT_WAS_SLOW
    assert "real and irrelevant" in gcn.WHY_IT_WAS_SLOW


def test_the_three_instances_are_recorded_as_ONE_finding() -> None:
    """The vacuity pool, the narrow glob and the unrepresentative benchmark are one mistake."""
    w = gcn.WHY_THE_BENCHMARK_LIED
    assert "THIRD TIME THIS LANE MADE THE SAME MISTAKE" in w
    assert "one finding, not three" in w
    assert "EXERCISED ON THE POPULATION IT WILL JUDGE" in w


def test_a_sample_inverts_the_reading_and_the_inversion_is_stated() -> None:
    alt = gcn.SUITE_LEG_ALTERNATIVE
    assert "only the POSITIVE is sound" in alt
    assert "NOT MEASURED" in alt
    assert "can never establish the count of guards that are not" in alt


def test_the_static_count_claims_only_the_population() -> None:
    c = gcn.STATIC_COUNT
    assert c["guards"] == 1556
    assert c["raise_sites"] == 1742
    assert c["guards"] + c["reraise"] + c["unconditional"] + c["loop"] == c["raise_sites"]
    assert sum(c["guards_by_directory"].values()) == c["guards"]
    assert "nothing about any guard's state" in c["what_it_establishes"]
    assert "NOT_ESTABLISHED will be the largest class" in c["what_it_establishes"]


def test_the_count_records_that_the_tree_moved_under_it() -> None:
    assert "1,736 raise sites at 25db524 and 1,742" in gcn.STATIC_COUNT["tree_moves_under_the_count"]
    assert "NOT instrument drift" in gcn.STATIC_COUNT["tree_moves_under_the_count"]


def test_the_instrumented_run_rule_is_carried_verbatim() -> None:
    """Supervisor, 2026-10-12 19:25, on this lane's own 160-minute lock hold."""
    r = gcn.INSTRUMENTED_RUNS_RULE
    assert r.startswith("an INSTRUMENTED or traced run never takes the machine-wide suite lock")
    assert "scoped to the files it needs" in r
    assert "past 45 minutes is reported to the coordinator" in r


def test_the_rule_records_that_this_lane_was_the_cause_not_a_bystander() -> None:
    assert "NO COMMIT" in gcn.WHY_THAT_RULE_EXISTS
    assert "cause rather than a bystander" in gcn.WHY_THAT_RULE_EXISTS
    assert "The lock did" in gcn.WHY_THAT_RULE_EXISTS


def test_duration_estimates_must_be_measured_over_a_meaningful_sample() -> None:
    d = gcn.DURATION_ESTIMATES_ARE_MEASURED
    assert "90-second window" in d
    assert "four and a half hours" in d
    assert "reported as 'rate unknown' rather than as an estimate" in d


def test_the_veil_is_recorded_as_instrumentation_not_as_an_exemption() -> None:
    """This lane filed an instrumented run under the exemption by calling it uninstrumented."""
    v = gcn.VEIL_IS_INSTRUMENTATION
    assert "That is WRONG and the wording" in v
    assert "70% slower" in v
    assert "sidestepped by vocabulary rather than by argument" in v


def test_the_replacement_method_tests_the_condition_instead_of_simulating_it() -> None:
    m = gcn.BARE_CONDITION_METHOD
    assert "git clone --no-local" in m
    assert "13 base packages" in m
    assert "OUTSIDE the suite lock" in m
    assert "right constraint and the wrong solution" in m
