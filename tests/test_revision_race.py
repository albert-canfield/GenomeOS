# SPDX-License-Identifier: AGPL-3.0-or-later
"""The revision race: recorded at the write, named and decided at the rebuild (2026-10-02).

A result records HEAD twice. `code_cleanliness.git_sha` is sampled while the writer computes;
`code.git_sha` is stamped when `save_result` writes. In the checkout several lanes share, a peer can
commit between the two, and the result then names two commits. A rebuild builds its worktree at
`code.git_sha`, so the cleanliness sha it recomputes there cannot equal the recorded one.

An exemption for that leaf was proposed and refused: silencing it would hide that HEAD moved during the
write, which is the one thing it says that nothing else does. So the leaf stays compared and the race is
made decidable -- named, then settled by `git diff A B` over the result's own counting path.

Both cases are planted here in a real repository, from real commits, and the verdict comes back from real
git. Neither is a fixture that cannot fail:

* the BENIGN case commits a change OFF the counting path, and the test asserts that the unrestricted
  diff between those same two revisions is NON-empty -- so "benign" is the pathspec doing its work and
  not "nothing changed";
* the REAL case commits a change ON the counting path and the test asserts the changed file is named.

Drop the pathspec and the benign case reads real; stop running git and the real case reads benign.
"""

import importlib.util
import json
import os
import subprocess
from pathlib import Path
from typing import Any

import planted_cleanliness
import pytest

from genomeos import manifest as mf
from genomeos import results
from genomeos.results import save_result

BENIGN = "benign race: the counting path is identical at both revisions."

ENTRY = "scripts/planted_entry.py"
MODULE = "scripts/planted_module.py"
OFF_PATH = "docs/planted_notes.md"

#: The entry puts its own directory on sys.path exactly as `counting_path` requires, so the import of
#: the planted module resolves to a file and the counting path really holds it.
ENTRY_SOURCE = """\
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import planted_module

print(planted_module.VALUE)
"""

#: Variables that would point git at another repository's index, directory or work tree. The functions
#: under test shell out with the inherited environment, and `scripts/commit_own.sh` runs the suite with
#: GIT_INDEX_FILE set to a private index, so without this the planted repository would be read against a
#: different one and every verdict below would be meaningless.
GIT_ENV_TO_DROP = ("GIT_INDEX_FILE", "GIT_DIR", "GIT_WORK_TREE", "GIT_OBJECT_DIRECTORY")


@pytest.fixture(autouse=True)
def _no_inherited_git_environment(monkeypatch: pytest.MonkeyPatch) -> None:
    for name in GIT_ENV_TO_DROP:
        monkeypatch.delenv(name, raising=False)


def _git(repo: Path, *args: str) -> str:
    r = subprocess.run(
        ["git", "-C", str(repo), *args], check=True, capture_output=True, text=True, env=dict(os.environ)
    )
    return r.stdout.strip()


@pytest.fixture
def planted(tmp_path: Path) -> dict[str, Any]:
    """A repository with three commits: the counting path changes between B and C and not between A and B.

    A -- everything as written.
    B -- `docs/planted_notes.md` edited. OFF the counting path.
    C -- `scripts/planted_module.py` edited. ON the counting path.
    """
    repo = tmp_path / "planted"
    (repo / "scripts").mkdir(parents=True)
    (repo / "docs").mkdir(parents=True)
    (repo / MODULE).write_text("VALUE = 1\n")
    (repo / ENTRY).write_text(ENTRY_SOURCE)
    (repo / OFF_PATH).write_text("first\n")
    _git(repo, "init", "-b", "main")
    _git(repo, "config", "user.email", "lane@example.invalid")
    _git(repo, "config", "user.name", "planted")
    _git(repo, "config", "commit.gpgsign", "false")
    _git(repo, "add", "--", ENTRY, MODULE, OFF_PATH)
    _git(repo, "commit", "-m", "A")
    a = _git(repo, "rev-parse", "HEAD")

    (repo / OFF_PATH).write_text("second\n")
    _git(repo, "add", "--", OFF_PATH)
    _git(repo, "commit", "-m", "B")
    b = _git(repo, "rev-parse", "HEAD")

    (repo / MODULE).write_text("VALUE = 2\n")
    _git(repo, "add", "--", MODULE)
    _git(repo, "commit", "-m", "C")
    c = _git(repo, "rev-parse", "HEAD")

    # the counting path is computed by the shared function from the planted imports, not written here
    block = mf.code_cleanliness(ENTRY, (), repo)
    return {"repo": repo, "a": a, "b": b, "c": c, "block": block}


def _result(planted: dict[str, Any], sampled: str, stamped: str) -> dict[str, Any]:
    """A result whose cleanliness block was sampled at `sampled` and whose write stamped `stamped`."""
    block = {**planted["block"], "git_sha": sampled}
    return {
        "result": "planted",
        mf.KEY: {"code": {"git_sha": stamped, "dirty": False}, "code_cleanliness": block},
    }


_spec = importlib.util.spec_from_file_location("mr", Path("scripts/manifest_rebuild.py"))
mr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mr)


# --- the planting itself, so neither case below can be vacuous -------------------------------------


def test_the_counting_path_holds_the_module_and_not_the_notes(planted: dict[str, Any]) -> None:
    path = planted["block"]["counting_path"]
    assert MODULE in path, "the planted module must be on the counting path or the real case proves nothing"
    assert ENTRY in path
    assert OFF_PATH not in path, "the off-path file must be off it or the benign case proves nothing"


def test_the_benign_pair_really_differs_outside_the_counting_path(planted: dict[str, Any]) -> None:
    """The counterfactual that keeps "benign" from meaning "the two revisions are the same commit"."""
    everything = _git(planted["repo"], "diff", "--name-only", planted["a"], planted["b"])
    assert everything.splitlines() == [OFF_PATH], (
        "A and B must differ, and only off the counting path: otherwise an implementation that ignored "
        "the pathspec would still read benign"
    )


# --- case one: a race whose counting path is identical at both revisions ---------------------------


def test_an_identical_counting_path_reads_benign(planted: dict[str, Any]) -> None:
    race = mr.revision_race(_result(planted, planted["a"], planted["b"]), planted["repo"])
    assert race["found"] is True, "a race between two different revisions is found, benign or not"
    assert race["counting_path_files_changed"] == []
    assert BENIGN in race["verdict"]


def test_the_benign_race_is_still_named_and_printed(planted: dict[str, Any]) -> None:
    """Benign is not silent: the name carries both revisions, in the direction HEAD moved."""
    race = mr.revision_race(_result(planted, planted["a"], planted["b"]), planted["repo"])
    assert race["race"] == f"revision race ({planted['a'][:7]} → {planted['b'][:7]})"
    assert race["race"] in race["verdict"]
    assert race["stamps"]["revision_race"] is True


def test_the_benign_verdict_rests_on_a_pathspec_that_was_actually_passed(planted: dict[str, Any]) -> None:
    argv = mr.revision_race(_result(planted, planted["a"], planted["b"]), planted["repo"])["git_diff_argv"]
    assert argv[:5] == ["git", "-C", str(planted["repo"]), "diff", "--name-only"]
    assert "--" in argv and MODULE in argv[argv.index("--") + 1 :], (
        "the diff must be restricted to the counting path, and the counting path must reach it"
    )


# --- case two: a race where one counting-path file changed -----------------------------------------


def test_a_changed_counting_path_file_reads_real_and_is_named(planted: dict[str, Any]) -> None:
    race = mr.revision_race(_result(planted, planted["b"], planted["c"]), planted["repo"])
    assert race["found"] is True
    assert race["counting_path_files_changed"] == [MODULE]
    assert "a REAL difference" in race["verdict"]
    assert MODULE in race["verdict"], "a real difference names the files that changed"
    assert BENIGN not in race["verdict"]


def test_only_the_counting_path_file_is_named_when_both_kinds_changed(planted: dict[str, Any]) -> None:
    """A to C changed the notes AND the module; only the counting-path file may be reported."""
    race = mr.revision_race(_result(planted, planted["a"], planted["c"]), planted["repo"])
    assert race["counting_path_files_changed"] == [MODULE]
    assert OFF_PATH not in race["verdict"]


# --- no race, and the states that cannot be established --------------------------------------------


def test_two_equal_stamps_are_not_a_race(planted: dict[str, Any]) -> None:
    race = mr.revision_race(_result(planted, planted["c"], planted["c"]), planted["repo"])
    assert race["found"] is False
    assert "no revision race" in race["verdict"]
    assert BENIGN not in race["verdict"], "no race and a benign race are different claims"


def test_a_missing_cleanliness_sha_cannot_say_and_is_not_agreement(planted: dict[str, Any]) -> None:
    bare = {"result": "planted", mf.KEY: {"code": {"git_sha": planted["c"], "dirty": False}}}
    race = mr.revision_race(bare, planted["repo"])
    assert race["found"] is None
    assert race["stamps"]["agree"] is None and race["stamps"]["revision_race"] is None
    assert "cannot say" in race["verdict"]


def test_an_unknown_revision_is_undecided_and_never_benign(planted: dict[str, Any]) -> None:
    unknown = "0" * 40
    race = mr.revision_race(_result(planted, unknown, planted["c"]), planted["repo"])
    assert race["found"] is True
    assert "undecided" in race["verdict"]
    assert BENIGN not in race["verdict"]


def test_an_absent_counting_path_is_refused_rather_than_called_benign(planted: dict[str, Any]) -> None:
    r = _result(planted, planted["a"], planted["b"])
    r[mf.KEY]["code_cleanliness"]["counting_path"] = []
    race = mr.revision_race(r, planted["repo"])
    assert race["found"] is True
    assert BENIGN not in race["verdict"]
    assert "nothing to diff" in race["verdict"]


# --- the leaf stays compared: the race explains it, it does not excuse it ---------------------------


def test_the_race_leaves_are_not_an_exemption_list() -> None:
    for field in mr.RACE_FIELDS:
        assert field not in mr.ENVIRONMENT_FIELDS, f"{field} must not be set aside as environment"
        assert field not in mr.RUN_FIELDS
    real, env = mr.environment_differences([f"{mr.RACE_FIELDS[1]}: 'aaa' vs 'bbb'"])
    assert real and not env, "a revision-race leaf is reported as a difference like any other"


# --- the write records both stamps and flags the disagreement ---------------------------------------


def _full(tmp_path: Path) -> dict[str, Any]:
    src = tmp_path / "in.tsv"
    src.write_text("chrom\tstart\tend\nchr21\t0\t10\n")
    return {
        "sources": [{"accession": "ENCSR000XXX", "version": "v1"}],
        "inputs": [mf.input_entry(src, partition="heldout")],
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {"threshold": 0.5},
        "exclusions": [],
        "partitions": {"heldout": "the test"},
        # 2026-10-02: the entry script is counted now (genomeos.manifest.entry_script), so under pytest
        # `sys.argv[0]` is a test runner and no in-process call can return a passing block. This
        # scaffolding block is built in a planted repository whose entry script really is committed and
        # clean, which is what it always meant to assert, instead of naming a file of this shared
        # checkout and reddening when a peer edits it (tests/planted_cleanliness.py).
        "code_cleanliness": planted_cleanliness.planted_clean_block(tmp_path, "race_clean"),
    }


def _saved(tmp_path: Path, cleanliness_sha: str) -> dict[str, Any]:
    """A result written by `save_result` whose cleanliness block was sampled at `cleanliness_sha`.

    The disagreement is made the way a real one is made: a real commit of this repository in the
    cleanliness block, and this repository's real HEAD in the stamp `save_result` takes. Nothing is
    patched, so a `save_result` that stopped recording the two stamps fails here.
    """
    m = _full(tmp_path)
    m["code_cleanliness"] = {**m["code_cleanliness"], "git_sha": cleanliness_sha}
    out = tmp_path / "out"
    out.mkdir(exist_ok=True)
    p = save_result("planted", {"value": 7}, out, manifest=m)
    return json.loads(p.read_text())[mf.KEY][mf.REVISION_STAMPS]


def test_a_write_whose_head_moved_records_both_stamps_and_the_flag(tmp_path: Path) -> None:
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    peer = subprocess.check_output(["git", "rev-parse", "HEAD~1"], text=True).strip()
    stamps = _saved(tmp_path, peer)
    assert stamps["code_cleanliness_git_sha"] == peer, "the sha the writer sampled is kept"
    assert stamps["code_git_sha"] == head, "the sha the write stamped is kept"
    assert stamps["revision_race"] is True and stamps["agree"] is False
    assert peer[:7] in stamps["reading"] and head[:7] in stamps["reading"]


def test_a_write_whose_head_held_still_records_agreement(tmp_path: Path) -> None:
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    stamps = _saved(tmp_path, head)
    assert stamps["revision_race"] is False and stamps["agree"] is True
    assert "no revision race" in stamps["reading"]


def test_the_flag_does_not_refuse_the_write(tmp_path: Path) -> None:
    """A race is a fact to record, not an incompleteness: nothing is quarantined over one."""
    peer = subprocess.check_output(["git", "rev-parse", "HEAD~1"], text=True).strip()
    m = _full(tmp_path)
    m["code_cleanliness"] = {**m["code_cleanliness"], "git_sha": peer}
    out = tmp_path / "out"
    out.mkdir()
    p = save_result("planted", {"value": 7}, out, manifest=m, strict=True)
    body = json.loads(p.read_text())
    assert body[mf.KEY]["complete"] is True
    assert body[mf.KEY][mf.REVISION_STAMPS]["revision_race"] is True
    assert not results.quarantine_dir(out).exists()


# --- the real-world case the ruling was raised on ---------------------------------------------------


def test_the_recorded_astrocyte_race_is_the_benign_one() -> None:
    """`data/results/astroreg2_astrocyte_activity.json`: 53f3b33 sampled, 675e54a stamped, 130 files.

    Checked against the committed result and this repository's real history, not a copy of it, so the
    reference point the ruling named is a test and not a note. Skipped only if the result is absent.
    """
    p = Path("data/results/astroreg2_astrocyte_activity.json")
    if not p.exists():
        pytest.skip(f"{p} is not in this checkout")
    race = mr.revision_race(json.loads(p.read_text()), Path.cwd())
    assert race["found"] is True
    assert race["race"] == "revision race (53f3b33 → 675e54a)"
    assert race["counting_path_count"] == 130
    assert race["counting_path_files_changed"] == []
    assert BENIGN in race["verdict"]
    # The sentence that makes this a measurement and not a coincidence: the two revisions DO differ, and
    # the one file they differ in is off the counting path. So the empty diff above is the pathspec doing
    # its work, not nothing having changed. Without this line the assertion could pass on two revisions
    # that were the same tree, which is the vacuous form of the same test.
    everything = subprocess.check_output(
        ["git", "diff", "--name-only", "53f3b33", "675e54a"], text=True
    ).split()
    assert everything == ["docs/ROADMAP.md"], everything
    counting = json.loads(p.read_text())[mf.KEY]["code_cleanliness"]["counting_path"]
    assert "docs/ROADMAP.md" not in counting
