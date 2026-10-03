"""The three holes `lane-silencerebuild` measured in scripts/manifest_rebuild.py, each with its planted case.

It ran the committed tool against two committed results and got NO VERDICT on either, for reasons that
were defects in the tool rather than differences in a value (commit fbc4c88):

(A) a run that read back ITS OWN OUTPUT was refused a verdict before a single leaf was compared, though
    `link_machine_local_inputs` already knows that concept on the input side
    (`not_linked_because_it_is_the_output`) and the output path held nothing before the command ran;
(B) every RELATIVE read was resolved with `os.path.realpath` in the TOOL'S working directory rather than
    the child's, so it resolved into the wrong checkout and was dropped as "not under data/" with
    nothing in the report saying a path had been discarded -- 26 of silenceragree's 3,340 recorded reads,
    two of which the manifest DECLARES and which the reconciliation then reported as never opened;
(C) an input declared as an absolute path was refused even when it resolved inside the repository root,
    which left 9 committed results unrebuildable at step 3.

Two of the three fixes make the tool refuse LESS, so each one's narrowing is planted BOTH WAYS here: the
case that must now pass, and beside it the case that must still refuse. A fix that made the tool blind
would show up as the second case passing.
"""

import importlib.util
import json
import os
import subprocess
from pathlib import Path

import pytest

from genomeos import manifest as mf

spec = importlib.util.spec_from_file_location("mr_reads", Path("scripts/manifest_rebuild.py"))
mr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mr)


def _declared(paths: dict[str, str]) -> dict:
    return {"inputs": [{"path": p, "sha256": s, "partition": None} for p, s in paths.items()]}


# --- (A) a run that read its own output ---------------------------------------------------------------


OUTPUT = "data/results/silenceragree_registration.json"


def test_a_read_of_this_run_s_own_output_is_not_a_read_of_an_undeclared_input():
    """The planted case: the only undeclared read under data/ is the file this run writes, and the tool
    measured that the worktree held nothing there. `read_but_not_declared` must come back empty, which
    is what lets `rebuild` go on to compare the leaves instead of returning before the first one."""
    got = mr.reads_against_declared(
        Path(),
        _declared({"data/cache/x.tsv": "a" * 64}),
        ["data/cache/x.tsv", OUTPUT],
        [],
        output=OUTPUT,
        output_held_bytes_before_the_run=False,
    )
    assert got["read_but_not_declared"] == []
    assert got["read_but_not_declared_because_it_is_this_run_s_own_output"] == [OUTPUT]
    assert got["every_read_is_declared_or_is_this_run_s_own_output"] is True
    assert OUTPUT in got["own_output_exempt_because"]


def test_the_exemption_does_not_print_as_the_stronger_claim():
    """The claim "every read is declared" and the claim that allows this run's own output are
    different claims, and the second is never printable as the first: the output WAS read and is NOT
    declared."""
    got = mr.reads_against_declared(
        Path(), _declared({}), [OUTPUT], [], output=OUTPUT, output_held_bytes_before_the_run=False
    )
    assert got["every_read_is_declared"] is False
    assert got["every_read_is_declared_or_is_this_run_s_own_output"] is True


def test_an_undeclared_read_that_is_not_the_output_is_still_refused():
    """The narrowing must be exactly "its own declared output" and nothing wider. A second undeclared
    file read by the same run keeps the refusal, and the output does not shelter it."""
    other = "data/results/unknown_chr21.json"
    got = mr.reads_against_declared(
        Path(),
        _declared({}),
        [OUTPUT, other],
        [],
        output=OUTPUT,
        output_held_bytes_before_the_run=False,
    )
    assert got["read_but_not_declared"] == [other]
    assert got["every_read_is_declared_or_is_this_run_s_own_output"] is False


@pytest.mark.parametrize(
    "candidate",
    [
        "data/results/silenceragree_registration.json.gz",
        "data/results/silenceragree_registration_v2.json",
        "data/results/silenceragree.json",
        "data/results",
    ],
)
def test_the_exemption_is_one_exact_path_and_not_a_prefix_or_a_directory(candidate):
    """Not a prefix, not the directory, not the name with another suffix: a directory or prefix match
    would have swept in every committed result, which is the shape the framework-read exemption was
    forced into an exact path over."""
    got = mr.reads_against_declared(
        Path(),
        _declared({}),
        [candidate],
        [],
        output=OUTPUT,
        output_held_bytes_before_the_run=False,
    )
    assert got["read_but_not_declared"] == [candidate]
    assert got["read_but_not_declared_because_it_is_this_run_s_own_output"] == []


@pytest.mark.parametrize("held", [True, None])
def test_the_output_is_not_exempt_when_the_path_was_not_measured_empty(held):
    """The condition is PROVEN, not assumed. If the worktree already held bytes at that path, the run
    read bytes no sha256 pins and the refusal stands; and an unmeasured None is never read as False."""
    got = mr.reads_against_declared(
        Path(), _declared({}), [OUTPUT], [], output=OUTPUT, output_held_bytes_before_the_run=held
    )
    assert got["read_but_not_declared"] == [OUTPUT]
    assert got["read_but_not_declared_because_it_is_this_run_s_own_output"] == []


def test_the_absence_is_measured_before_the_command_runs_and_not_after():
    """The command writes that path, so after it runs no observation can tell what was there first.
    The measurement must sit above the line that starts the child."""
    src = Path("scripts/manifest_rebuild.py").read_text()
    measured = src.index("output_held_bytes_before_the_run = ")
    started = src.index('run = subprocess.run(["uv", "run", "python", *argv]')
    assert measured < started, "the output's absence is measured after the command could have written it"


# --- (B) a relative read -------------------------------------------------------------------------------


def test_a_relative_read_is_resolved_against_the_child_s_directory_and_not_this_tool_s(tmp_path):
    """The measured defect: `os.path.realpath("data/results/x.json")` resolves in the TOOL'S cwd. The
    child ran in the worktree, so the same spelling must resolve there -- and here the two directories
    are different, which is exactly the condition the old code could not tell apart."""
    wt = tmp_path / "worktree"
    (wt / "data" / "results").mkdir(parents=True)
    (wt / "data" / "results" / "x.json").write_text("{}")
    got = mr.reads_against_declared(
        wt, _declared({"data/results/x.json": "a" * 64}), ["data/results/x.json"], []
    )
    assert got["reads_under_data"] == ["data/results/x.json"]
    assert got["relative_reads_resolved_under_data"] == ["data/results/x.json"]
    assert got["relative_reads_recorded_count"] == 1
    assert got["read_but_not_declared"] == []
    assert os.path.realpath(str(wt)) == got["relative_reads_resolved_against"]


def test_a_relative_read_the_manifest_does_not_declare_is_now_caught(tmp_path):
    """The fix is not only a rescue: a relative read of an UNDECLARED file was invisible too, so the
    tool's "every read is declared" was weaker than it read as for every result in the registry."""
    wt = tmp_path / "worktree"
    (wt / "data" / "results").mkdir(parents=True)
    got = mr.reads_against_declared(wt, _declared({}), ["data/results/unknown_chr21.json"], [])
    assert got["read_but_not_declared"] == ["data/results/unknown_chr21.json"]
    assert got["every_read_is_declared"] is False


def test_a_relative_read_through_a_linked_store_is_matched_too(tmp_path):
    """A store read resolves outside the worktree through a symlink. Relative or absolute, it is the
    same file, and dropping it would make the reconciliation vacuous for every writer that reads one."""
    wt, store = tmp_path / "worktree", tmp_path / "machine" / "data" / "cache"
    (store / "silenceragree").mkdir(parents=True)
    (store / "silenceragree" / "rese_records.json").write_text("{}")
    (wt / "data").mkdir(parents=True)
    (wt / "data" / "cache").symlink_to(store)
    got = mr.reads_against_declared(
        wt,
        _declared({"data/cache/silenceragree/rese_records.json": "a" * 64}),
        ["data/cache/silenceragree/rese_records.json"],
        [str(store.resolve())],
    )
    assert got["reads_under_data"] == ["data/cache/silenceragree/rese_records.json"]
    assert got["read_but_not_declared"] == []


def test_a_relative_read_that_resolves_to_nothing_is_named_rather_than_dropped(tmp_path):
    """What the resolution still cannot do: a child that chdir'd opens a relative path somewhere this
    tool does not know. Such a path must land in a named list, because silence is the defect itself."""
    wt = tmp_path / "worktree"
    (wt / "data").mkdir(parents=True)
    got = mr.reads_against_declared(wt, _declared({}), ["somewhere/else/table.tsv"], [])
    assert got["relative_reads_that_are_not_there"] == ["somewhere/else/table.tsv"]
    assert got["reads_under_data"] == []


def test_an_absolute_read_outside_data_is_still_not_an_input(tmp_path):
    """The interpreter, the virtual environment and the standard library are not inputs, and the
    relative-path fix must not have widened what counts as one."""
    got = mr.reads_against_declared(tmp_path, _declared({}), ["/usr/lib/python3.12/json/__init__.py"], [])
    assert got["reads_under_data"] == []
    assert got["relative_reads_recorded"] == []


# --- (C) an absolute input that resolves inside the repository -----------------------------------------


def _machine(tmp_path: Path, bodies: dict[str, bytes]) -> tuple[Path, Path]:
    """A machine checkout whose data/results is git-ignored, and an empty worktree beside it.

    The same fixture as tests/test_manifest_rebuild_links.py, because the case under test is one of
    those inputs declared the other way round.
    """
    root, wt = tmp_path / "machine", tmp_path / "worktree"
    (root / "data" / "results").mkdir(parents=True)
    (wt / "data" / "results").mkdir(parents=True)
    (root / ".gitignore").write_text("data/results/*.json\n")
    subprocess.run(["git", "-C", str(root), "init", "-q"], check=True)
    for name, body in bodies.items():
        (root / "data" / "results" / name).write_bytes(body)
    return root, wt


def _abs_entry(root: Path, name: str) -> dict:
    """An input declared the way `input_entry` records `ROOT / REGISTRATION`: the absolute path verbatim."""
    p = root / "data" / "results" / name
    digest, size, _ = mf.sha256_of(p)
    return {"path": str(p), "sha256": digest, "bytes": size, "partition": None}


def test_an_in_repo_absolute_input_is_relativised_linked_and_checked(tmp_path):
    """The planted case that must now pass: silenceragree's fourth input, in miniature. It is declared
    absolutely, it is git-ignored, it is absent from the worktree, and its bytes are on this machine."""
    root, wt = _machine(tmp_path, {"registration.json": b'{"value": 7}'})
    declared = _abs_entry(root, "registration.json")

    links = mr.link_machine_local_inputs(wt, root, [declared])
    assert links["linked"] == ["data/results/registration.json"]

    got = mr.check_inputs(wt, [declared], root, machine_local=set(links["linked"]))
    assert got["unavailable"] == []
    assert got["checked"] == 1 and got["declared"] == 1
    assert got["machine_local_opened"] == 1, "counted by the relative form, or the figure silently drops"
    assert got["relativised_absolute_declarations"][0]["declared"] == declared["path"]
    assert got["relativised_absolute_declarations"][0]["read_as"] == "data/results/registration.json"
    assert len(got["relativised_absolute_declarations"]) == 1
    assert got["entries"][0]["read_as"] == "data/results/registration.json"


def test_an_out_of_repo_absolute_input_is_still_refused_by_name(tmp_path):
    """The other way round, and this is the narrowing: a path outside the repository root is not a path
    this worktree reads, and the refusal must stand with the root it was tested against named."""
    root, wt = _machine(tmp_path, {"registration.json": b"{}"})
    elsewhere = tmp_path / "another_machine" / "registration.json"
    elsewhere.parent.mkdir(parents=True)
    elsewhere.write_text("{}")
    digest, size, _ = mf.sha256_of(elsewhere)
    declared = {"path": str(elsewhere), "sha256": digest, "bytes": size, "partition": None}

    links = mr.link_machine_local_inputs(wt, root, [declared])
    assert links["linked"] == []
    got = mr.check_inputs(wt, [declared], root, machine_local=set())
    assert got["checked"] == 0
    assert len(got["unavailable"]) == 1
    assert "is an absolute path outside the linked data stores" in got["unavailable"][0]
    assert os.path.realpath(str(root)) in got["unavailable"][0]
    assert got["relativised_absolute_declarations"] == []


def test_an_in_repo_absolute_input_whose_bytes_differ_is_a_named_failure(tmp_path):
    """Relativising must not become a way to pass: the bytes linked are hashed against the declared
    sha256 afterwards, so a machine-local input that has since changed stops the rebuild by name."""
    root, wt = _machine(tmp_path, {"registration.json": b'{"value": 7}'})
    declared = _abs_entry(root, "registration.json")
    (root / "data" / "results" / "registration.json").write_bytes(b'{"value": 8}')

    links = mr.link_machine_local_inputs(wt, root, [declared])
    got = mr.check_inputs(wt, [declared], root, machine_local=set(links["linked"]))
    assert got["checked"] == 0
    assert any("has different bytes" in u for u in got["unavailable"])


def test_the_boundary_is_a_resolved_path_and_a_separator_not_a_string_prefix(tmp_path):
    """`..`, symlinks and a sibling whose name merely starts with the root's. A string prefix test would
    have let `machine-old/` in and kept a `..` path out of nothing."""
    root, wt = _machine(tmp_path, {"registration.json": b"{}"})
    sibling = tmp_path / "machine-old" / "data" / "results"
    sibling.mkdir(parents=True)
    (sibling / "registration.json").write_text("{}")
    assert str(sibling).startswith(str(root))

    assert mr.repo_relative_declaration(str(sibling / "registration.json"), root) is None
    assert mr.repo_relative_declaration(str(root / "data" / "results"), root) == "data/results"
    # a path that climbs out and back in resolves inside, and is read by where it LANDS
    assert (
        mr.repo_relative_declaration(str(root / "data" / ".." / "data" / "results"), root) == "data/results"
    )
    # a relative declaration that climbs out of the repository is outside by the same rule
    assert mr.repo_relative_declaration("../elsewhere/x.json", root) is None
    assert mr.repo_relative_declaration("data/results/x.json", root) == "data/results/x.json"
    assert wt.exists()


def test_the_nine_committed_results_that_declare_an_absolute_input_are_measured_not_asserted():
    """The nine are a measurement over the committed registry, so the figure is re-measured here rather
    than quoted: an absolute declaration under the repository root is what (C) now reads.

    This is a statement about the files as they stand. It does not re-verdict any of them: whether a
    rebuild of one succeeds depends on its inputs' bytes and its command, which this does not run.
    """
    results = Path("data/results")
    if not results.is_dir():
        pytest.skip("data/results is git-ignored: this measurement needs the registry on this machine")
    absolute: dict[str, list[str]] = {}
    for f in sorted(results.glob("*.json")):
        try:
            payload = json.loads(f.read_text())
        except (OSError, ValueError):
            continue
        m = payload.get(mf.KEY) if isinstance(payload, dict) else None
        inputs = m.get("inputs") if isinstance(m, dict) else None
        if not isinstance(inputs, list):
            continue
        abs_paths = [
            i["path"]
            for i in inputs
            if isinstance(i, dict) and isinstance(i.get("path"), str) and os.path.isabs(i["path"])
        ]
        if abs_paths:
            absolute[f.stem] = abs_paths
    # every absolutely declared input in the registry points inside this repository, which is why (C)
    # covers all of them: not one of them names another machine
    root = Path(__file__).resolve().parent.parent
    outside = [
        p for paths in absolute.values() for p in paths if mr.repo_relative_declaration(p, root) is None
    ]
    assert outside == [], f"an absolute declaration outside this repository: {outside[:3]}"
    assert len(absolute) >= 9, f"{len(absolute)} results declare an absolute input, 9 were measured"
