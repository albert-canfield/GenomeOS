# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The fast pre-push's window, its fallbacks, and the two assertions that make its store-free leg mean
something.

WHAT WENT WRONG AND WHAT THESE TESTS ARE FOR. On 2026-10-03 three CI runs went red -- 3 failures, then
6 -- on one class: a test that assumes this machine's git-ignored data stores. The local full suite
cannot see that class, because this checkout HAS the stores. `scripts/fast_prepush.py` catches it by
running the targeted tests a second time with the stores removed. Three of the tests below are the
plants for its three window cases and one of them, (c), is the one that is usually got wrong.

PLANT (c) ASSERTS THE FILE COUNT, NOT THE MESSAGE. When the last-green-CI record is stale, the rule is
to USE ITS SHA ANYWAY and widen, because a stale record still names a sha CI once went green on. A
hook that printed the staleness warning and then narrowed to the since-origin window would pass any
test that only looked for the warning text -- so `test_c_...` asserts the NUMBER OF FILES the window
selects, requires it to be strictly greater than the since-origin window's, and requires the wider set
to CONTAIN the narrower one. `test_c_the_count_assertion_is_not_vacuous` then shows the planted
repository really does distinguish the two, so a green (c) cannot come from both windows happening to
be the same size.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from datetime import UTC, datetime, timedelta
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import fast_prepush as fp  # noqa: E402
import last_green_ci as lg  # noqa: E402

# --------------------------------------------------------------------------- a planted repository


def _git(repo: Path, *args: str) -> str:
    out = subprocess.run(["git", *args], cwd=repo, capture_output=True, text=True, check=True)
    return out.stdout.strip()


def _commit(repo: Path, rel: str, body: str, message: str) -> str:
    (repo / rel).parent.mkdir(parents=True, exist_ok=True)
    (repo / rel).write_text(body)
    _git(repo, "add", "--", rel)
    _git(repo, "commit", "-q", "-m", message)
    return _git(repo, "rev-parse", "HEAD")


@pytest.fixture
def planted(tmp_path: Path) -> dict:
    """A four-commit repository whose two windows are DIFFERENT SIZES, which is what plant (c) needs.

    The history is the shape of the incident: the green CI sha is two commits behind the remote tip,
    so the files changed since green (3) are strictly more than the files changed since origin (1).
    That gap is the whole reason the window is measured from green and not from origin -- a file
    changed in a commit that reached origin under a red run is invisible to the narrow window.
    """
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(repo, "init", "-q", "-b", "dev")
    _git(repo, "config", "user.email", "planted@example.invalid")
    _git(repo, "config", "user.name", "planted")
    _git(repo, "config", "commit.gpgsign", "false")
    green = _commit(repo, "tests/test_one.py", "def test_one():\n    assert True\n", "one")
    _commit(repo, "tests/test_two.py", "def test_two():\n    assert True\n", "two")
    at_origin = _commit(repo, "tests/test_three.py", "def test_three():\n    assert True\n", "three")
    head = _commit(repo, "tests/test_four.py", "def test_four():\n    assert True\n", "four")
    _git(repo, "update-ref", "refs/remotes/origin/dev", at_origin)
    return {"repo": repo, "green": green, "at_origin": at_origin, "head": head}


def _record(planted: dict, written: datetime, sha: str | None = None) -> dict:
    return lg.write(
        planted["repo"],
        sha or planted["green"],
        run_id="18200000000",
        date="2026-10-01T09:14:00Z",
        now=written,
    )


# --------------------------------------------------------------------------- plant (a): FRESH


def test_a_a_fresh_record_takes_the_window_to_the_green_sha(planted: dict) -> None:
    now = datetime.now(UTC)
    _record(planted, now)
    w = fp.resolve_window(planted["repo"], now=now + timedelta(hours=1))
    assert w.case == "fresh", w.as_dict()
    assert w.base == planted["green"]
    assert not w.block
    files = fp.changed_test_files(planted["repo"], w.base)
    assert files == ["tests/test_four.py", "tests/test_three.py", "tests/test_two.py"]
    assert len(files) == 3


# --------------------------------------------------------------------------- plant (b): MISSING


def test_b_a_missing_record_falls_back_to_since_origin_and_does_not_block(planted: dict) -> None:
    assert not lg.last_green_path(planted["repo"]).exists()
    w = fp.resolve_window(planted["repo"])
    assert w.case == "missing", w.as_dict()
    assert w.base == planted["at_origin"], "the fallback base is the merge base with origin/dev"
    assert w.block is False, "a fresh clone looks exactly like this and nothing is wrong"
    assert any("no last-green-CI record" in n for n in w.notices), w.notices
    assert fp.changed_test_files(planted["repo"], w.base) == ["tests/test_four.py"]


def test_b_an_unreadable_record_is_reported_as_unreadable_and_not_as_missing(planted: dict) -> None:
    """Two different facts get two different words, so nobody is told a file is absent when it is there."""
    path = lg.last_green_path(planted["repo"])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("{not json")
    w = fp.resolve_window(planted["repo"])
    assert w.case == "missing" and w.base == planted["at_origin"]
    assert any("unreadable" in n for n in w.notices), w.notices
    assert w.block is False


def test_b_a_record_naming_a_commit_this_checkout_lacks_is_its_own_case(planted: dict) -> None:
    path = lg.last_green_path(planted["repo"])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"sha": "0" * 40, "written": datetime.now(UTC).isoformat()}))
    w = fp.resolve_window(planted["repo"])
    assert w.case == "unknown-sha", w.as_dict()
    assert w.base == planted["at_origin"]
    assert w.block is False


def test_b_no_remote_ref_either_leaves_always_run_and_still_does_not_block(tmp_path: Path) -> None:
    repo = tmp_path / "bare"
    repo.mkdir()
    _git(repo, "init", "-q", "-b", "dev")
    _git(repo, "config", "user.email", "planted@example.invalid")
    _git(repo, "config", "user.name", "planted")
    _commit(repo, "tests/test_one.py", "def test_one():\n    assert True\n", "one")
    w = fp.resolve_window(repo)
    assert w.case == "no-base" and w.base is None
    assert w.block is False
    assert fp.changed_test_files(repo, w.base) == []


# --------------------------------------------------------------------------- plant (c): STALE


def test_c_a_stale_record_uses_its_sha_anyway_and_the_window_is_wider_BY_FILE_COUNT(planted: dict) -> None:
    """THE PLANT THAT MATTERS. Stale must WIDEN, and the evidence is the count, not the warning.

    A hook that printed the staleness warning and then fell back to the since-origin window would
    satisfy a test that only read the message. So the message is checked AND the selection is
    measured: strictly more files, and every file of the narrow window among them.
    """
    written = datetime.now(UTC) - timedelta(hours=72)
    _record(planted, written)
    w = fp.resolve_window(planted["repo"], now=datetime.now(UTC))

    assert w.case == "stale", w.as_dict()

    # THE COUNT FIRST, deliberately. A mutation of scripts/fast_prepush.py that keeps the staleness
    # warning word for word and narrows the base to the since-origin window was run against this test
    # on 2026-10-03: with the sha assertion first, the failure arrived on the sha and the count
    # assertion was never reached, so the count could not be shown to be doing any work. In this order
    # the mutation fails on the COUNT line, which is the assertion this plant exists to make.
    stale_files = fp.changed_test_files(planted["repo"], w.base)
    origin_files = fp.changed_test_files(planted["repo"], fp.since_origin_base(planted["repo"]))

    assert len(origin_files) == 1
    assert len(stale_files) == 3, (
        f"the stale window selected {len(stale_files)} file(s) where the record's own sha gives 3: "
        f"staleness changed the window instead of only warning about it"
    )
    assert len(stale_files) > len(origin_files), (
        f"the stale window selected {len(stale_files)} file(s) and the since-origin window "
        f"{len(origin_files)}: staleness NARROWED the window, which is 'we know less, so we check "
        f"less' -- the defect this test exists to catch"
    )
    assert set(origin_files) < set(stale_files), "the wider window must CONTAIN the narrower one"

    assert w.base == planted["green"], "the stale record's OWN sha, not the remote's"
    assert w.block is False
    assert w.age_seconds is not None and w.age_seconds > fp.STALE_SECONDS

    assert any("48 h" in n or "48 h)" in n for n in w.notices), w.notices
    assert any("USING ITS SHA ANYWAY" in n for n in w.notices), w.notices


def test_c_the_count_assertion_is_not_vacuous(planted: dict) -> None:
    """The planted repository really does make the two windows different sizes.

    Without this, a (c) that passed because both windows happened to select the same files would read
    as a demonstration that staleness widens, and it would be no such thing.
    """
    from_green = fp.changed_test_files(planted["repo"], planted["green"])
    from_origin = fp.changed_test_files(planted["repo"], planted["at_origin"])
    assert len(from_green) == 3 and len(from_origin) == 1
    assert from_green != from_origin


def test_c_a_record_just_under_the_threshold_is_fresh(planted: dict) -> None:
    """The boundary is checked from both sides, so "stale" cannot be a state nothing ever leaves."""
    now = datetime.now(UTC)
    _record(planted, now - timedelta(hours=47))
    assert fp.resolve_window(planted["repo"], now=now).case == "fresh"
    _record(planted, now - timedelta(hours=49))
    assert fp.resolve_window(planted["repo"], now=now).case == "stale"


# --------------------------------------------------------------------------- the path, in a worktree


def test_the_record_lives_under_the_git_COMMON_dir_because_a_worktree_dot_git_is_a_FILE(
    planted: dict,
) -> None:
    """The platform assumption that cost three CI runs on 2026-10-03, checked rather than assumed.

    In a LINKED worktree `.git` is a regular file holding a `gitdir:` pointer, so a path built as
    `<root>/.git/genomeos-check/last-green-ci` is a path under a FILE and can never exist. The reader
    would then see "missing" forever and narrow its window on every push, silently. Both the fast
    pre-push and the existing push check run inside linked worktrees, so this is the live case.
    """
    repo = planted["repo"]
    wt = repo.parent / "linked"
    _git(repo, "worktree", "add", "-q", "--detach", str(wt), planted["head"])
    try:
        dot_git = wt / ".git"
        assert dot_git.is_file(), "a linked worktree's .git is a FILE; the whole point of this test"
        assert not dot_git.is_dir()
        assert dot_git.read_text().startswith("gitdir:")

        naive = wt / ".git" / "genomeos-check" / "last-green-ci"
        assert not naive.exists(), "a path under a file cannot exist, which is why --git-dir is wrong"

        from_main = lg.last_green_path(repo)
        from_worktree = lg.last_green_path(wt)
        assert from_worktree.resolve() == from_main.resolve(), (
            "the record must be ONE file for the checkout and all of its worktrees"
        )

        written = _record(planted, datetime.now(UTC))
        assert json.loads(from_worktree.read_text())["sha"] == written["sha"], (
            "written from the main checkout and read from inside the linked worktree"
        )
        w = fp.resolve_window(wt, now=datetime.now(UTC))
        assert w.case == "fresh" and w.base == planted["green"]
    finally:
        _git(repo, "worktree", "remove", "--force", str(wt))


def test_the_record_sits_beside_the_verdict_status_files(planted: dict) -> None:
    """Pinned to the convention `genomeos/verdict.py` already set, so the two cannot drift apart."""
    from genomeos import verdict

    assert (
        lg.last_green_path(planted["repo"]).parent.resolve() == verdict.status_dir(planted["repo"]).resolve()
    )
    assert lg.DIRNAME == "genomeos-check"


def test_the_writer_refuses_a_sha_this_checkout_does_not_have(planted: dict) -> None:
    with pytest.raises(SystemExit) as e:
        lg.write(planted["repo"], "0" * 40, run_id=None, date=None)
    assert "not a commit in this checkout" in str(e.value)
    assert not lg.last_green_path(planted["repo"]).exists()


# --------------------------------------------------------------------------- the store-free leg


def test_assert_stores_absent_passes_on_an_empty_tree_and_REFUSES_on_a_present_store(tmp_path: Path) -> None:
    """The assertion that makes the leg's green mean something, checked in both directions.

    A leg that silently found the stores present would pass everything and report nothing wrong. So
    its failure must be a refusal, and the refusal must name the store it found.
    """
    tree = tmp_path / "tree"
    (tree / "data").mkdir(parents=True)
    clean = fp.assert_stores_absent(tree)
    assert clean["ok"] is True
    assert clean["present"] == []
    assert set(clean["checked"]) == {"data/knowledge", "data/cache", "data/reference"}

    (tree / "data" / "cache").mkdir()
    dirty = fp.assert_stores_absent(tree)
    assert dirty["ok"] is False
    assert dirty["present"] == ["data/cache"]
    assert "REFUSED" in dirty["why"] and "data/cache" in dirty["why"]


def test_assert_stores_absent_sees_a_store_reached_through_a_symlink(tmp_path: Path) -> None:
    """A linked store is present as far as every test is concerned, and linking is how they get in."""
    real = tmp_path / "real_store"
    real.mkdir()
    tree = tmp_path / "tree"
    (tree / "data").mkdir(parents=True)
    os.symlink(real, tree / "data" / "reference")
    got = fp.assert_stores_absent(tree)
    assert got["ok"] is False and got["present"] == ["data/reference"]


def test_the_import_path_assertion_REFUSES_a_tree_that_carries_no_source(tmp_path: Path) -> None:
    """Non-vacuity for the second required assertion: it is known to fire, not merely known to pass.

    Against a directory with no `genomeos/` in it, `import genomeos` resolves to the MAIN CHECKOUT's
    source, which this venv can reach -- the leg would then be testing the very tree it exists to
    exclude -- and the assertion must refuse by name.
    """
    got = fp.assert_import_path(ROOT, tmp_path)
    assert got["ok"] is False
    assert "REFUSED" in got["why"]
    assert str(ROOT) in got["imported"]
    assert str(tmp_path) in got["why"]


def test_run_pytest_calls_a_non_zero_exit_red_and_an_empty_set_nothing_to_run(tmp_path: Path) -> None:
    empty = fp.run_pytest(ROOT, None, [], "probe")
    assert empty["ok"] is True and empty["skipped"] == "nothing to run" and empty["exit"] is None


def test_leg_env_puts_the_worktree_first_and_never_syncs(tmp_path: Path) -> None:
    """`uv sync` mutates the one environment every session on this machine shares; it is never run.

    And no leg resolves an environment of its own: every child is `sys.executable`, so the hook
    cannot behave one way here and another in CI over where a venv happens to live -- which would be
    the very defect class it exists to catch.
    """
    env = fp.leg_env(ROOT, tmp_path)
    assert env["UV_NO_SYNC"] == "1"
    assert "UV_PROJECT_ENVIRONMENT" not in env
    assert env["PYTHONPATH"] == str(tmp_path)
    assert "PYTHONPATH" not in fp.leg_env(ROOT, None)
    assert sys.executable == fp.PYTHON
    source = (ROOT / "scripts" / "fast_prepush.py").read_text()
    assert '"uv", "run"' not in source, "the legs run sys.executable; uv is not in the inner loop"


# --------------------------------------------------------------------------- both directions of the marker

#: a test that reads a git-ignored store. WITHOUT the marker it must make the store-free leg red.
UNMARKED = """\
from pathlib import Path


def test_reads_a_git_ignored_store_without_declaring_it():
    assert Path("data/cache/fast_prepush_plant.txt").read_text().strip() == "planted"
"""

#: the same read, declared. It must SKIP BY NAME where the store is absent, and the leg stays green.
MARKED = """\
from pathlib import Path

import pytest


@pytest.mark.needs_local_data(
    "data/cache",
    how="nothing fetches it: this store is planted by this test file and removed with the worktree",
)
def test_reads_a_git_ignored_store_and_declares_it():
    assert Path("data/cache/fast_prepush_plant.txt").read_text().strip() == "planted"
"""


def test_an_undeclared_read_of_a_git_ignored_store_makes_the_store_free_leg_RED(tmp_path: Path) -> None:
    """The defect class of 2026-10-03, planted, and the leg must block on it.

    Planted inside the worktree's own `tests/` and not in a scratch directory, because pytest finds a
    conftest by walking up from the test file: outside `tests/` the project's `needs_local_data`
    wiring is not loaded at all and the plant would prove nothing about this project.
    """
    out = fp.store_free_leg(
        ROOT,
        ["tests/test_fast_prepush_plant_unmarked.py"],
        plant={"tests/test_fast_prepush_plant_unmarked.py": UNMARKED},
    )
    assert out["worktree_added"] is True
    assert out["stores"]["ok"] is True, out["stores"]
    assert out["import_path"]["ok"] is True, out["import_path"]
    assert out["ok"] is False, "an undeclared read of an absent store must make this leg red"
    assert out["pytest"]["exit"] not in (0, 5), out["pytest"]
    assert "FileNotFoundError" in out["pytest"]["output_tail"] or "1 failed" in out["pytest"]["summary"]
    assert out["removal"]["ok"] is True, out["removal"]


def test_the_SAME_read_with_the_marker_passes_the_store_free_leg(tmp_path: Path) -> None:
    """The other direction: the marker moves where a test runs, so the leg is green and says why."""
    out = fp.store_free_leg(
        ROOT,
        ["tests/test_fast_prepush_plant_marked.py"],
        plant={"tests/test_fast_prepush_plant_marked.py": MARKED},
    )
    assert out["stores"]["ok"] is True, out["stores"]
    assert out["import_path"]["ok"] is True, out["import_path"]
    assert out["ok"] is True, out.get("why") or out["pytest"]
    assert out["pytest"]["exit"] == 0, out["pytest"]
    assert "skipped" in out["pytest"]["summary"], out["pytest"]["summary"]
    assert out["removal"]["ok"] is True, out["removal"]


def test_the_import_path_assertion_reads_the_worktrees_own_source() -> None:
    """The second required assertion: PYTHONPATH precedence is checked, not assumed."""
    out = fp.store_free_leg(ROOT, [], plant=None)
    assert out["import_path"]["ok"] is True, out["import_path"]
    imported = Path(out["import_path"]["imported"])
    assert str(imported).startswith(out["import_path"]["worktree"] + os.sep), out["import_path"]
    assert imported.name == "__init__.py" and imported.parent.name == "genomeos"
    assert out["removal"]["ok"] is True, out["removal"]


# --------------------------------------------------------------------------- the selection


def test_the_selection_always_includes_ALWAYS_RUN_and_names_it_from_one_place() -> None:
    members = fp.always_run(ROOT)
    assert members and all(m.startswith("tests/test_") for m in members)
    assert set(fp.select(ROOT, []).files) == set(members)
    assert "ALWAYS_RUN" not in fp.select(ROOT, []).notices


def test_a_changed_helper_with_no_tests_brings_in_the_modules_that_import_it() -> None:
    """The seven imported-not-collected files under tests/, handled by what they are."""
    sel = fp.select(ROOT, ["tests/local_data.py"])
    importers = fp.importers_of(ROOT, "local_data")
    assert importers, "tests/local_data.py is imported by at least one test module"
    assert set(importers) <= set(sel.files)
    assert "tests/local_data.py" not in sel.files, "it carries no tests; naming it would collect nothing"
    assert sel.helpers["tests/local_data.py"] == importers


def test_a_changed_helper_that_does_carry_tests_is_named_directly() -> None:
    """`tests/always_run.py` has four tests and a name `python_files` never matches; pytest collects
    an explicitly named file anyway, so it is run rather than resolved to importers."""
    assert fp.has_tests(ROOT / "tests" / "always_run.py") is True
    sel = fp.select(ROOT, ["tests/always_run.py"])
    assert "tests/always_run.py" in sel.files
    assert "tests/always_run.py" not in sel.helpers


def test_a_changed_conftest_is_a_notice_and_not_an_expansion() -> None:
    sel = fp.select(ROOT, ["tests/conftest.py"])
    assert sel.helpers["tests/conftest.py"] == []
    assert any("LOADED BY EVERY pytest run" in n for n in sel.notices), sel.notices
    assert set(sel.files) == set(fp.always_run(ROOT))


def test_every_one_of_the_seven_helpers_is_classified_and_none_is_silently_dropped() -> None:
    seven = (
        "always_run",
        "committed_data",
        "conftest",
        "extras_lock",
        "guard_injection",
        "local_data",
        "planted_cleanliness",
    )
    rels = [f"tests/{s}.py" for s in seven]
    assert all((ROOT / r).exists() for r in rels)
    sel = fp.select(ROOT, rels)
    for rel in rels:
        assert rel in sel.files or rel in sel.helpers, f"{rel} was neither run nor resolved"
    assert fp.has_tests(ROOT / "tests" / "always_run.py")
    assert not any(fp.has_tests(ROOT / f"tests/{s}.py") for s in seven if s != "always_run")


def test_a_file_deleted_since_the_window_is_named_and_skipped() -> None:
    sel = fp.select(ROOT, ["tests/test_a_file_that_was_removed.py"])
    assert "tests/test_a_file_that_was_removed.py" not in sel.files
    assert any("is not in the tree now" in n for n in sel.notices), sel.notices


def test_nothing_here_holds_an_allowlist_or_a_skip_for_the_window() -> None:
    """The rule `tests/always_run.py` sets for itself: this adds a path, it never excuses one.

    Asked of the CODE and not of the prose. The first version of this test grepped the file for
    "exempt" and failed on the paragraph that says there is no exemption -- a word test cannot tell a
    mechanism from a sentence about one. So: no binding whose name offers an allowlist or an
    exemption, no skip marker, and -- the claim that actually matters -- ALWAYS_RUN survives every
    input the selection takes, including one that names each of its own members as a deletion.
    """
    import ast as _ast

    tree = _ast.parse((ROOT / "scripts" / "fast_prepush.py").read_text())
    named = []
    for node in _ast.walk(tree):
        targets = (
            [node.target]
            if isinstance(node, _ast.AnnAssign)
            else (node.targets if isinstance(node, _ast.Assign) else [])
        )
        named += [
            _ast.unparse(t)
            for t in targets
            if any(w in _ast.unparse(t).lower() for w in ("allow", "exempt", "skiplist"))
        ]
        if isinstance(node, _ast.FunctionDef):
            named += [_ast.unparse(d) for d in node.decorator_list if "skip" in _ast.unparse(d)]
    assert named == [], named

    members = set(fp.always_run(ROOT))
    for probe in ([], list(members), ["tests/conftest.py"], ["tests/test_a_file_that_was_removed.py"]):
        assert members <= set(fp.select(ROOT, probe).files), probe
