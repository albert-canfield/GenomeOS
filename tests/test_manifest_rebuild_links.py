"""The machine-local inputs a rebuild must read, linked read-only, and the two ways it must refuse.

A result computed from other results declares inputs under data/results, which is git-ignored: reader v1's
cached DNase peak sets are on one machine and in no commit. Until 2026-10-02 the rebuild linked only
data/reference, data/knowledge and data/cache, so those inputs were named absent one by one and no
comparison was attempted at all -- response_map_increment2 reached 111 of 410 inputs opened, context_evidence
the same wall. They are linked here instead, and the point of these tests is the other half of that change:
linking must not become a way to pass. So the two failing cases are constructed rather than described -- an
input that is absent, and an input whose bytes differ -- and each is required to stop the rebuild by name.
"""

import importlib.util
import json
import os
import subprocess
from pathlib import Path

import pytest

from genomeos import manifest as mf

spec = importlib.util.spec_from_file_location("mrl", Path("scripts/manifest_rebuild.py"))
mrl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mrl)


def _machine(tmp_path: Path, bodies: dict[str, bytes]) -> tuple[Path, Path]:
    """A machine checkout whose data/results is git-ignored, and an empty worktree beside it."""
    root, wt = tmp_path / "machine", tmp_path / "worktree"
    (root / "data" / "results").mkdir(parents=True)
    (wt / "data" / "results").mkdir(parents=True)
    # the real .gitignore names machine-local results by pattern (data/results/dnase_*_chr*.bed.gz and
    # the rest), never the whole directory, so a tracked result can sit beside them
    (root / ".gitignore").write_text("data/results/*.json\ndata/results/*.bed.gz\n")
    subprocess.run(["git", "-C", str(root), "init", "-q"], check=True)
    for name, body in bodies.items():
        (root / "data" / "results" / name).write_bytes(body)
    return root, wt


def _entry(root: Path, name: str, sha: str | None = None) -> dict:
    p = f"data/results/{name}"
    digest, size, _ = mf.sha256_of(root / p)
    return {"path": p, "sha256": sha or digest, "bytes": size, "partition": None}


# --- the two required failing cases -------------------------------------------------------------------


def test_a_declared_input_that_is_absent_stops_the_rebuild_by_name(tmp_path):
    """Case one: the manifest declares an input nothing on this machine holds."""
    root, wt = _machine(tmp_path, {"present.json": b"{}"})
    gone = {"path": "data/results/gone.json", "sha256": "a" * 64, "bytes": 2}
    inputs = [_entry(root, "present.json"), gone]

    links = mrl.link_machine_local_inputs(wt, root, inputs)
    assert links["absent_on_this_machine"] == ["data/results/gone.json"]
    assert links["linked"] == ["data/results/present.json"]

    got = mrl.check_inputs(wt, inputs, root, machine_local=set(links["linked"]))
    assert got["declared"] == 2
    assert got["checked"] == 1, "an input nobody opened cannot be counted as checked"
    assert any("gone.json is absent" in u for u in got["unavailable"])
    # and the verdict is refused, not merely annotated
    assert mrl.no_verdict_reason(got["declared"], got["checked"]) is not None


def test_a_declared_input_whose_bytes_differ_is_reported_as_a_difference(tmp_path):
    """Case two: the file is there and linked, and its bytes are not the bytes the manifest declares.

    This is the failure that would flatter: the input resolves, so a tool that only asked whether the path
    existed would count it and report a clean check over bytes it had not compared.
    """
    root, wt = _machine(tmp_path, {"drifted.json": b'{"x": 1}'})
    declared = _entry(root, "drifted.json")
    (root / "data" / "results" / "drifted.json").write_bytes(b'{"x": 2}')  # same length, different bytes

    links = mrl.link_machine_local_inputs(wt, root, [declared])
    assert links["linked"] == ["data/results/drifted.json"], "it is present, so it is linked"

    got = mrl.check_inputs(wt, [declared], root, machine_local=set(links["linked"]))
    assert got["checked"] == 0
    assert got["entries"][0]["sha256_matches"] is False
    assert any("different bytes" in u for u in got["unavailable"])
    assert mrl.no_verdict_reason(got["declared"], got["checked"]) is not None


# --- read-only is the point ----------------------------------------------------------------------------


def test_a_linked_input_cannot_be_written_and_the_machines_copy_keeps_its_mode(tmp_path):
    """A reader that opens an input for writing must fail loudly, not mutate the only copy there is."""
    root, wt = _machine(tmp_path, {"only_copy.json": b'{"kept": true}'})
    src = root / "data" / "results" / "only_copy.json"
    before = oct(src.stat().st_mode & 0o777)

    mrl.link_machine_local_inputs(wt, root, [_entry(root, "only_copy.json")])
    link = wt / "data" / "results" / "only_copy.json"

    assert link.read_bytes() == b'{"kept": true}', "the rebuild must be able to read it"
    assert not link.stat().st_mode & 0o222, "no write bit on the worktree's copy"
    with pytest.raises(PermissionError):
        link.open("w")
    with pytest.raises(PermissionError):
        link.open("ab")
    # the machine's one copy is untouched in both bytes and mode: the clone's chmod is its own inode's
    assert src.read_bytes() == b'{"kept": true}'
    assert oct(src.stat().st_mode & 0o777) == before
    assert not link.is_symlink(), "a symlink would carry the machine's mode and could not refuse a write"
    assert src.stat().st_ino != link.stat().st_ino, "a hard link's chmod would disarm the machine's copy"


def test_the_write_bits_come_back_so_the_worktree_can_be_removed(tmp_path):
    root, wt = _machine(tmp_path, {"a.json": b"1"})
    links = mrl.link_machine_local_inputs(wt, root, [_entry(root, "a.json")])
    mrl.unlock_links(wt, links["linked"])
    assert (wt / "data" / "results" / "a.json").stat().st_mode & 0o200
    assert os.stat(root / "data" / "results" / "a.json").st_mode & 0o777 == 0o644


# --- what may be linked, and what may not --------------------------------------------------------------


def test_only_a_git_ignored_path_is_linked_from_this_machine(tmp_path):
    """A path git does not ignore is read from the commit, whatever is sitting at it on this machine."""
    root, wt = _machine(tmp_path, {"local.json": b"1"})
    (root / "data" / "results" / "tracked.json").write_bytes(b"2")
    (root / ".gitignore").write_text("data/results/*.json\n!data/results/tracked.json\n")
    inputs = [_entry(root, "local.json"), _entry(root, "tracked.json")]

    links = mrl.link_machine_local_inputs(wt, root, inputs)
    assert links["linked"] == ["data/results/local.json"]
    assert links["not_linked_because_git_does_not_ignore_it"] == ["data/results/tracked.json"]
    assert not (wt / "data" / "results" / "tracked.json").exists()
    # and the one not linked is then named absent, never passed over
    got = mrl.check_inputs(wt, inputs, root, machine_local=set(links["linked"]))
    assert any("tracked.json is absent" in u for u in got["unavailable"])


def test_a_file_the_commit_carries_is_read_from_the_commit(tmp_path):
    """The worktree's own copy wins: linking must not put this machine's bytes over a committed file."""
    root, wt = _machine(tmp_path, {"both.json": b"machine"})
    (wt / "data" / "results" / "both.json").write_bytes(b"committed")
    links = mrl.link_machine_local_inputs(wt, root, [_entry(root, "both.json")])
    assert links["linked"] == []
    assert (wt / "data" / "results" / "both.json").read_bytes() == b"committed"


def test_the_result_being_rebuilt_is_an_output_and_is_never_linked(tmp_path):
    """Linking the output would hand the run a read-only file to write, and a false difference to report."""
    root, wt = _machine(tmp_path, {"out.json": b"{}"})
    links = mrl.link_machine_local_inputs(
        wt, root, [_entry(root, "out.json")], output="data/results/out.json"
    )
    assert links["not_linked_because_it_is_the_output"] == ["data/results/out.json"]
    assert links["linked"] == []


def test_nothing_outside_data_results_and_no_path_that_escapes_is_linked(tmp_path):
    root, wt = _machine(tmp_path, {"ok.json": b"1"})
    inputs = [
        _entry(root, "ok.json"),
        {"path": "data/knowledge/x.json", "sha256": "b" * 64},
        {"path": "data/results/../../escape.json", "sha256": "c" * 64},
        {"path": "/etc/hosts", "sha256": "d" * 64},
    ]
    links = mrl.link_machine_local_inputs(wt, root, inputs)
    assert links["linked"] == ["data/results/ok.json"]
    assert links["absent_on_this_machine"] == []


# --- the count that says how much rested on this machine -----------------------------------------------


def test_the_machine_local_count_names_what_came_from_this_machine(tmp_path):
    root, wt = _machine(tmp_path, {"local1.json": b"1", "local2.json": b"2"})
    (wt / "data" / "results" / "committed.json").write_bytes(b"3")
    inputs = [_entry(root, "local1.json"), _entry(root, "local2.json")]
    digest, size, _ = mf.sha256_of(wt / "data" / "results" / "committed.json")
    inputs.append({"path": "data/results/committed.json", "sha256": digest, "bytes": size})

    links = mrl.link_machine_local_inputs(wt, root, inputs)
    got = mrl.check_inputs(wt, inputs, root, machine_local=set(links["linked"]))
    assert (got["declared"], got["checked"]) == (3, 3)
    assert got["machine_local_opened"] == 2, "two of the three came from a git-ignored path"
    assert [e["path"] for e in got["entries"] if e.get("machine_local")] == [
        "data/results/local1.json",
        "data/results/local2.json",
    ]


def test_a_machine_local_input_is_counted_even_though_its_bytes_differ(tmp_path):
    """The count says what the rebuild rested on, not what passed, so it cannot be read as a pass rate."""
    root, wt = _machine(tmp_path, {"d.json": b"aaa"})
    declared = _entry(root, "d.json")
    (root / "data" / "results" / "d.json").write_bytes(b"bbb")
    links = mrl.link_machine_local_inputs(wt, root, [declared])
    got = mrl.check_inputs(wt, [declared], root, machine_local=set(links["linked"]))
    assert got["machine_local_opened"] == 1
    assert got["checked"] == 0
    assert got["unavailable"]


# --- the rule that may not be broken -------------------------------------------------------------------


def test_no_verdict_is_allowed_over_inputs_that_were_not_opened():
    assert mrl.no_verdict_reason(410, 111) is not None, "the figure of 2026-10-02 must refuse a verdict"
    assert mrl.no_verdict_reason(13, 1) is not None
    assert mrl.no_verdict_reason(0, 0) is not None, "no inputs declared is not a clean check"
    assert mrl.no_verdict_reason(410, 410) is None, "every input opened and hashed: a verdict may be given"


def test_the_verdict_keys_are_written_only_after_that_gate():
    """A later edit must not be able to move a verdict key above the gate without this failing."""
    src = Path("scripts/manifest_rebuild.py").read_text()
    gate = src.index("why = no_verdict_reason(")
    for key in ("identical_bytes_except_date_and_run", 'report["identical_bytes"]'):
        assert src.index(key) > gate, f"{key} is written before the no-verdict gate"


def test_declared_paths_names_a_groups_members_and_not_its_label():
    """A label is not a path, so a group contributes the files it names and its label is not looked for.

    An entry that is not marked a group contributes its own path, even one carrying a file count: at this
    point nothing distinguishes a label written before files_entry named its members from an input recorded
    on a directory, and it is _check_one_path that tells the two apart and names the first as unresolved.
    """
    group = {
        "path": "compiler_inputs:budget_chr*.json",
        "group": True,
        "members": [{"path": "data/results/budget_chr1.json", "sha256": "b" * 64}],
        "sha256": "c" * 64,
    }
    assert mrl.declared_paths([{"path": "data/results/one.json", "sha256": "a" * 64}, group]) == [
        "data/results/one.json",
        "data/results/budget_chr1.json",
    ]
    assert "compiler_inputs:budget_chr*.json" not in mrl.declared_paths([group])


def test_a_group_whose_members_were_never_named_is_not_linked_into_a_pass(tmp_path):
    """The twelve groups context_evidence declares are labels with no member list, written before
    files_entry named its members. Linking cannot help them: there is no path to link. The tool must keep
    saying so rather than let the label resolve to nothing and count."""
    root, wt = _machine(tmp_path, {"budget_chr1.json": b"1"})
    inputs = [{"path": "compiler_inputs:budget_chr*.json", "sha256": "e" * 64, "bytes": 1, "files": 24}]
    links = mrl.link_machine_local_inputs(wt, root, inputs)
    assert links["linked"] == []
    got = mrl.check_inputs(wt, inputs, root, machine_local=set(links["linked"]))
    assert got["checked"] == 0
    assert any("files_entry named its members" in u for u in got["unavailable"])


def test_the_report_carries_the_machine_local_count_under_a_name_that_says_so(tmp_path):
    """The reader of a rebuild must see at a glance how much of it depended on one machine."""
    src = Path("scripts/manifest_rebuild.py").read_text()
    assert '"inputs_satisfied_from_machine_local_paths"' in src
    assert json.dumps("machine_local_inputs")[1:-1] in src


# --- run-resource readings, exempt by exact key --------------------------------------------------------


def test_a_peak_memory_reading_is_ignored_like_a_wall_clock_second():
    """Supervisor's ruling of 2026-10-02: peak memory is a run-resource reading, the class of the timing
    keys, and not code cleanliness -- so it is exempt here and not through ENVIRONMENT_FIELDS.

    These two keys were the one leaf standing between two otherwise identical rebuilds and 0 differences:
    `/cost/peak_memory_mb` 1142.5 against 1148.8 in placement_audit, `/compute/peak_rss_mb` 1300.5 against
    1298.2 in context_contrast_feasibility.
    """
    assert mrl.is_resource("peak_rss_mb")
    assert mrl.is_resource("peak_memory_mb")
    a = {"cost": {"peak_memory_mb": 1142.5}, "compute": {"peak_rss_mb": 1300.5}}
    b = {"cost": {"peak_memory_mb": 1148.8}, "compute": {"peak_rss_mb": 1298.2}}
    assert mrl.diff(a, b), "the raw payloads do differ"
    assert mrl.diff(mrl.comparable(a), mrl.comparable(b)) == [], "comparable() drops them"
    assert sorted(mrl.resource_paths(a)) == ["/compute/peak_rss_mb", "/cost/peak_memory_mb"]
    # and the exemption is not smuggled into the code-cleanliness list, which is a different kind of thing
    assert not any("peak" in f for f in mrl.ENVIRONMENT_FIELDS)
    for key in mrl.RESOURCE_KEYS:
        assert not any(f.endswith(key) for f in mrl.ENVIRONMENT_FIELDS)


def test_a_scientific_leaf_with_a_similar_name_is_still_compared():
    """What an exact-key list buys. Each of these would fall to a loose pattern -- `*_mb`, `peak_*`,
    `*memory*` -- and each is a number a run computed, so each must still be compared."""
    similar = {
        "peak_signal_mb": 3.5,  # a measured signal, which a `*_mb` pattern would swallow
        "peak_rss_mb_per_cell": 12.0,  # a per-cell figure, which a `peak_rss_mb*` prefix would swallow
        "memory_mb_budget": 2048,  # a declared parameter, which a `*memory*` pattern would swallow
        "peak_openness": 0.91,  # a measured peak, which a `peak_*` pattern would swallow
        "peak_memory_mb_of_the_cache": 900,  # the exact key is not a prefix either
    }
    for key, value in similar.items():
        assert not mrl.is_resource(key), f"{key} is not a run-resource reading"
        a, b = {"verdict": {key: value}}, {"verdict": {key: value + 1}}
        assert mrl.diff(mrl.comparable(a), mrl.comparable(b)), f"a change in {key} must be a difference"


def test_what_was_set_aside_is_reported_and_not_merely_dropped():
    src = Path("scripts/manifest_rebuild.py").read_text()
    assert '"resource_fields_ignored"' in src
    assert '"timing_fields_ignored"' in src
    assert "resource_keys" in src, "the leaf reconciliation must name the cause too"


# --- the inputs-only mode may never report a verdict ---------------------------------------------------


def test_inputs_only_stops_before_the_command_and_reports_no_verdict():
    """The mode item 10's capacity rule forces on a lane that may not read the per-element cache. It must
    buy the input figure and nothing else: no comparison, no verdict, rebuilt false."""
    src = Path("scripts/manifest_rebuild.py").read_text()
    branch = src.index("if inputs_only:")
    assert branch < src.index('subprocess.run(["uv", "run", "python"'), "it must stop before the command"
    assert branch < src.index('report["differences"]'), "it must stop before any comparison"
    tail = src[branch : src.index("env = {k: v for k, v in os.environ.items()")]
    assert '"rebuilt": False' in tail
    assert "no verdict is " in tail
