"""The amended rebuild rule: the science must reproduce at a leaf, the environment may not excuse it.

Adopted 2026-10-02. As first written the rule asked for 0 differences, which a result honest about a shared
checkout could never meet: it records the other lanes' outstanding files, and a clean worktree has none. So
one class of field is exempt by exact path, and two fields must hold on both sides so the exemption can
never cover a result whose own code was uncommitted.
"""

import hashlib
import importlib.util
from pathlib import Path

from genomeos import manifest as mf

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


# ---------------------------------------------------------------------------------------------------
# A grouped input is opened and hashed file by file.
#
# Until 2026-10-02 the rebuild resolved every declared input by `entry["path"]`, where
# `manifest.files_entry` records a label rather than a path. So a result whose inputs were grouped had
# every group reported "absent" with not one file opened or hashed, and the entry never reached the
# report's `inputs` list, which carried no denominator and so read as though every input had matched.
# The error ran in the direction that flattered the result: such a result was excused as resting on an
# unavailable dependency while its files sat on disk, unchecked and uncompared.
# ---------------------------------------------------------------------------------------------------

BODIES = {"one.txt": b"the first", "two.txt": b"the second", "three.txt": b"the third"}


def _a_group(tmp_path, monkeypatch, bodies=None):
    """A files_entry over three files in `tmp_path`, recorded the way a writer records it."""
    (tmp_path / "data").mkdir(exist_ok=True)
    for name, body in (bodies or BODIES).items():
        (tmp_path / "data" / name).write_bytes(body)
    monkeypatch.chdir(tmp_path)
    return mf.files_entry("reader_v1_dnase_peak_sets", [f"data/{n}" for n in (bodies or BODIES)])


def test_a_grouped_input_is_opened_and_hashed_member_by_member(tmp_path, monkeypatch):
    """The old code could not pass this: it resolved the label as a path, found nothing there, and
    reported an input absent whose three files were on disk."""
    entry = _a_group(tmp_path, monkeypatch)
    got = mr.check_inputs(tmp_path, [entry], tmp_path)
    assert got["unavailable"] == []
    assert not any("absent" in u for u in got["unavailable"])
    assert (got["declared"], got["checked"]) == (1, 1)
    assert got["files_opened"] == 3, "each file in the group must be opened and hashed on its own"
    assert got["entries"][0]["members_opened_and_hashed"] == 3
    assert got["entries"][0]["members_declared"] == 3
    assert got["entries"][0]["sha256_matches"] is True


def test_a_group_with_one_missing_member_cannot_pass(tmp_path, monkeypatch):
    entry = _a_group(tmp_path, monkeypatch)
    (tmp_path / "data" / "two.txt").unlink()
    got = mr.check_inputs(tmp_path, [entry], tmp_path)
    assert got["checked"] == 0, "a group missing a member must not count as checked"
    assert any("data/two.txt" in u for u in got["unavailable"]), got["unavailable"]
    assert got["entries"][0]["sha256_matches"] is False
    assert len(got["entries"]) == got["declared"], "a failed input must not leave the list"


def test_a_member_whose_bytes_changed_is_named_even_at_the_same_length(tmp_path, monkeypatch):
    """A byte count cannot catch this; only hashing each member can."""
    entry = _a_group(tmp_path, monkeypatch)
    (tmp_path / "data" / "two.txt").write_bytes(b"the SECOND")
    got = mr.check_inputs(tmp_path, [entry], tmp_path)
    assert got["checked"] == 0
    assert got["entries"][0]["members_with_different_bytes"] == ["data/two.txt"]
    assert got["entries"][0]["sha256_matches"] is False


def test_a_group_recorded_before_its_members_were_named_is_a_loud_failure(tmp_path):
    """A manifest written before 2026-10-02 names no members, so nothing can be opened for it. That is
    a reason to stop and say which input could not be resolved, never a reason to report a verdict."""
    legacy = {"path": "reader_v1_dnase_peak_sets", "sha256": "a" * 64, "bytes": 27697843, "files": 312}
    got = mr.check_inputs(tmp_path, [legacy], tmp_path)
    assert got["checked"] == 0
    assert got["files_opened"] == 0
    assert any("312" in u and "label" in u for u in got["unavailable"]), got["unavailable"]
    assert len(got["entries"]) == 1


def test_every_declared_input_is_accounted_for(tmp_path, monkeypatch):
    """compared plus not compared equals declared, for inputs as for leaves: a list an unchecked input
    can leave invites the reader to assume the list is the denominator."""
    entry = _a_group(tmp_path, monkeypatch)
    (tmp_path / "data" / "alone.txt").write_bytes(b"on its own")
    alone = mf.input_entry("data/alone.txt")
    missing = {"path": "data/never_written.txt", "sha256": "b" * 64, "partition": None}
    got = mr.check_inputs(tmp_path, [entry, alone, missing], tmp_path)
    assert got["declared"] == 3
    assert got["checked"] == 2
    assert len(got["entries"]) == 3
    assert got["unavailable"], "the third input must be named, not dropped"


def test_a_manifest_declaring_no_inputs_cannot_report_a_clean_check(tmp_path):
    for inputs in ([], None, "data/x.json"):
        got = mr.check_inputs(tmp_path, inputs, tmp_path)
        assert got["unavailable"], f"{inputs!r} must stop the rebuild"
        assert (got["declared"], got["checked"], got["files_opened"]) == (0, 0, 0)


def test_the_writer_and_the_checker_share_one_group_digest(tmp_path, monkeypatch):
    """Two implementations of the same digest can answer differently about which bytes it covers."""
    entry = _a_group(tmp_path, monkeypatch)
    names = sorted(m["path"] for m in entry["members"])
    digest, total, members = mf.group_digest(names, root=tmp_path)
    assert digest == entry["sha256"]
    assert total == entry["bytes"]
    assert [m["sha256"] for m in members] == [m["sha256"] for m in entry["members"]]


def test_the_group_digest_is_unchanged_so_a_digest_recorded_earlier_still_verifies(tmp_path, monkeypatch):
    """Naming the members must not move the group's own digest, or every committed group would read as
    changed bytes."""
    entry = _a_group(tmp_path, monkeypatch)
    h = hashlib.sha256()  # the formula as it stood before the members were named
    for name in sorted(f"data/{n}" for n in BODIES):
        h.update(name.encode() + b"\0")
        h.update((tmp_path / name).read_bytes())
    assert entry["sha256"] == h.hexdigest()


def test_an_absolute_input_outside_the_linked_stores_is_not_reported_as_checked(tmp_path):
    """`worktree / path` leaves the worktree when path is absolute, so the bytes hashed would not be the
    bytes the rebuild reads. Under the linked stores they are the same file; anywhere else they are not."""
    (tmp_path / "elsewhere.json").write_bytes(b"{}")
    entry = {"path": str(tmp_path / "elsewhere.json"), "sha256": "c" * 64, "partition": None}
    got = mr.check_inputs(tmp_path / "worktree", [entry], tmp_path / "root")
    assert got["checked"] == 0
    assert any("absolute path outside the linked data stores" in u for u in got["unavailable"])


def test_a_recorded_group_carries_enough_to_open_its_files_again(tmp_path, monkeypatch):
    """The recorder's half of the defect: the entry's `path` is a label, so a reader that resolves an
    input by `path` opens nothing. What the group stood for has to be in the entry."""
    entry = _a_group(tmp_path, monkeypatch)
    assert not (tmp_path / entry["path"]).exists(), "the label is not a path: nothing opens at it"
    assert mr.is_group(entry)
    assert [m["path"] for m in entry["members"]] == sorted(f"data/{n}" for n in BODIES)
    for m in entry["members"]:
        assert (tmp_path / m["path"]).exists(), f"{m['path']} must be openable from the entry alone"
