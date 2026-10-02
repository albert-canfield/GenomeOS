"""The amended rebuild rule: the science must reproduce at a leaf, the environment may not excuse it.

Adopted 2026-10-02. As first written the rule asked for 0 differences, which a result honest about a shared
checkout could never meet: it records the other lanes' outstanding files, and a clean worktree has none. So
one class of field is exempt by exact path, and two fields must hold on both sides so the exemption can
never cover a result whose own code was uncommitted.
"""

import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location("mr", Path("scripts/manifest_rebuild.py"))
mr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mr)


def test_the_exempt_fields_are_exact_paths_not_patterns():
    assert "/code_cleanliness/dirty" in mr.ENVIRONMENT_FIELDS
    assert "/result_manifest/code_cleanliness/foreign_uncommitted_code" in mr.ENVIRONMENT_FIELDS
    # a pattern would have swept these in; an exact-path list must not
    assert not any(f.endswith("/gain") or f.endswith("/positives") for f in mr.ENVIRONMENT_FIELDS)


def test_an_environment_difference_is_set_aside_and_reported():
    real, env = mr.environment_differences(
        ["/code_cleanliness/dirty: True vs False", "/code_cleanliness/foreign_uncommitted_code: 6 items vs 0"]
    )
    assert real == []
    assert len(env) == 2


def test_a_scientific_difference_is_never_set_aside():
    """The defect the rule exists to catch: a number that does not reproduce."""
    real, env = mr.environment_differences(
        ["/verdict/pooled_independent_loci_genome_wide: 22 vs 21", "/code_cleanliness/dirty: True vs False"]
    )
    assert real == ["/verdict/pooled_independent_loci_genome_wide: 22 vs 21"]
    assert len(env) == 1


def test_a_difference_inside_an_unexempt_part_of_the_cleanliness_record_is_real():
    real, _ = mr.environment_differences(["/code_cleanliness/counting_path: 45 items vs 44"])
    assert real == ["/code_cleanliness/counting_path: 45 items vs 44"]


def test_an_element_of_an_exempt_list_is_exempt():
    _, env = mr.environment_differences(["/code_cleanliness/foreign_uncommitted_code[2]: 'a' vs 'b'"])
    assert len(env) == 1


def test_uncommitted_own_code_fails_on_either_side():
    """The exemption must never excuse a result no commit reproduces."""
    bad = {"code_cleanliness": {"own_code_is_committed": False}}
    good = {"code_cleanliness": {"own_code_is_committed": True}}
    assert mr.must_hold_failures(bad, good), "an uncommitted original must fail"
    assert mr.must_hold_failures(good, bad), "an uncommitted rebuild must fail"
    assert mr.must_hold_failures(good, good) == []


def test_a_foreign_file_on_the_counting_path_fails():
    bad = {"code_cleanliness": {"foreign_uncommitted_code_on_the_counting_path": ["x.py"]}}
    good = {"code_cleanliness": {"foreign_uncommitted_code_on_the_counting_path": []}}
    assert mr.must_hold_failures(bad, good)
    assert mr.must_hold_failures(good, good) == []


def test_must_hold_is_checked_at_any_depth():
    nested = {"result_manifest": {"code_cleanliness": {"own_code_is_committed": False}}}
    assert mr.must_hold_failures(nested, {})


def test_leaves_counts_values_not_top_level_keys():
    """A count of top-level keys can hide a nested difference, which is why leaves are reported."""
    payload = {"a": 1, "b": {"c": 2, "d": [3, 4, 5]}}
    assert len(payload) == 2
    assert mr.leaves(payload) == 5


def test_the_leaf_counts_reconcile_with_their_denominator():
    """A count needs its denominator: compared plus not-compared must equal the whole file."""
    payload = {
        "date": "2026-10-02",
        "seconds": 1.5,
        "a": {"b": 1, "c": [2, 3]},
        "result_manifest": {"code": {"git_sha": "x", "dirty": False}, "complete": True},
    }
    r = mr.leaf_reconciliation(payload)
    assert r["total"] == mr.leaves(payload)
    assert r["compared"] + r["not_compared"] == r["total"]
    assert r["reconciles"] is True
    assert r["not_compared"] > 0, "the date, the timing key and the code block must be accounted for"


def test_the_dropped_leaves_are_named_by_cause_not_only_counted():
    r = mr.leaf_reconciliation({"date": "x", "a": 1})
    why = r["not_compared_because"]
    assert "date" in why["ignored_keys"]
    assert "timing" in why["timing_keys"].lower()
    assert "code" in why["manifest_code_block"]
