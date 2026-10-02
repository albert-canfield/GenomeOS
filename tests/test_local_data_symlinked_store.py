# SPDX-License-Identifier: AGPL-3.0-or-later
"""A store reached through a symbolic link must SKIP BY NAME, and never ERROR.

WHAT THIS PLANTS. On 2026-10-02 twenty tests ERRORED and two RED status files were written for trees
that were nothing of the kind -- `status-7446015a503fcc6f281173590591709d1f6e7aab` (4,418 passed, 0
failed, 20 errors) and `status-579b4949ea024ef98722a9ea14f9825478818724` (4,436 passed, 0 failed, 20
errors, which refused a push). Every one of the forty was the same line:

    local_data.NotMachineLocalError: git could not say whether
    data/cache/entex/alphagenome_track_metadata_copy.csv is ignored
    (fatal: pathspec '...' is beyond a symbolic link)

THE CAUSE, measured rather than reasoned, because two diagnoses were offered and both were wrong.
`.git/hooks/pre-push` was a STALE COPY taken at 04:34 which still did
`ln -s "$root/data/$d" "$tmp/data/$d"`, while `scripts/pre-push.sh` had moved that day to read-only
APFS clones (`scripts/link_stores_read_only.py`). The hook is a copy and not a symlink -- see
`scripts/install-hooks.sh` -- so editing the script in the tree changed nothing the push ran. The
push's worktree therefore had `data/cache` as a symbolic link into the main checkout, and git stops
at a symbolic link.

NOT the cause, and both were checked before this test was written. `$TMPDIR` lies behind
`/var -> private/var` on macOS, but a worktree under `$TMPDIR` with the stores as real directories
answers 0: git resolves a linked worktree's root physically and a relative pathspec is taken from the
process's physical cwd. And `realpath` before `check-ignore` does not fix it -- it asks about a
different NAME, and in the exact failing shape the resolved path lies in another working tree of the
same repository, which git refuses with 128 again. Both readings are in `local_data.ignored_by_name`.

WHAT THESE TESTS ASSERT, which is the difference between a fix and a cover-up. A test that needs a
store it cannot find must SKIP WITH THE STORE NAMED AND THE FETCH COMMAND GIVEN. It must not pass, it
must not be quietly deselected, and it must not error. Twenty silent passes would have been worse than
the red, because then nobody could tell "not run here" from "passed".
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import local_data
import pytest

#: The real store path from the forty errors, used verbatim so the planted case is the measured one.
FAILING_PATH = "data/cache/entex/alphagenome_track_metadata_copy.csv"


def _git(cwd: Path, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *args], cwd=str(cwd), capture_output=True, text=True)


def _repo(root: Path) -> Path:
    """A repository that ignores `data/cache`, with one commit, made without porcelain `git commit`."""
    root.mkdir(parents=True, exist_ok=True)
    _git(root, "init", "-q", "-b", "main", ".")
    (root / ".gitignore").write_text("data/cache\n")
    _git(root, "add", ".gitignore")
    tree = _git(root, "write-tree").stdout.strip()
    sha = _git(
        root,
        "-c",
        "user.email=t@example.invalid",
        "-c",
        "user.name=t",
        "commit-tree",
        tree,
        "-m",
        "init",
    ).stdout.strip()
    _git(root, "update-ref", "refs/heads/main", sha, "")
    return root


def test_git_itself_refuses_a_path_beyond_a_symlink_and_answers_for_a_real_directory(tmp_path) -> None:
    """The measurement the fix rests on, taken here so a git change cannot silently retire the fix."""
    root = _repo(tmp_path / "repo")
    (root / "store" / "entex").mkdir(parents=True)
    (root / "store" / "entex" / "x.csv").write_text("x\n")
    (root / "data").mkdir()
    (root / "data" / "cache").symlink_to(root / "store")

    linked = _git(root, "check-ignore", "-q", "--no-index", "data/cache/entex/x.csv")
    assert linked.returncode == 128, linked
    assert local_data.BEYOND_A_SYMLINK in linked.stderr

    (root / "data" / "cache").unlink()
    (root / "data" / "cache" / "entex").mkdir(parents=True)
    (root / "data" / "cache" / "entex" / "x.csv").write_text("x\n")
    real = _git(root, "check-ignore", "-q", "--no-index", "data/cache/entex/x.csv")
    assert real.returncode == 0, real


def test_a_symlinked_store_is_answered_by_name_and_does_not_raise(tmp_path, monkeypatch) -> None:
    """The twenty errors, planted: the same path, the same symlink, and now an ANSWER."""
    root = _repo(tmp_path / "repo")
    (root / "store" / "entex").mkdir(parents=True)
    (root / "store" / "entex" / "alphagenome_track_metadata_copy.csv").write_text("x\n")
    (root / "data").mkdir()
    (root / "data" / "cache").symlink_to(root / "store")
    monkeypatch.chdir(root)
    local_data.ignored_by_name.cache_clear()

    assert local_data.is_git_ignored(FAILING_PATH) is True
    ignored, answered_about = local_data.ignored_by_name(FAILING_PATH)
    assert (ignored, answered_about) == (True, "data/cache"), (ignored, answered_about)
    # and the marker check, which is what pytest's setup calls, no longer raises
    local_data.check_is_machine_local(FAILING_PATH)
    local_data.ignored_by_name.cache_clear()


def test_a_worktree_under_a_symlinked_tmpdir_was_never_the_cause(tmp_path, monkeypatch) -> None:
    """The second diagnosis, planted as the negative it is: a leading symlink above the repo is fine.

    `/var -> private/var` was proposed as the cause. Here the whole repository sits behind a symlinked
    parent, exactly as a `$TMPDIR` worktree does, and the stores are real directories: the question is
    answerable and the answer is "ignored". A fix aimed at this would have left the push red.
    """
    real = tmp_path / "physical"
    _repo(real / "repo")
    (real / "repo" / "data" / "cache" / "entex").mkdir(parents=True)
    (real / "repo" / "data" / "cache" / "entex" / "x.csv").write_text("x\n")
    (tmp_path / "behind").symlink_to(real)
    monkeypatch.chdir(tmp_path / "behind" / "repo")
    local_data.ignored_by_name.cache_clear()

    assert local_data.symlink_boundary("data/cache/entex/x.csv") is None
    assert local_data.ignored_by_name("data/cache/entex/x.csv") == (True, "data/cache/entex/x.csv")
    local_data.ignored_by_name.cache_clear()


def test_realpath_before_check_ignore_would_not_have_worked(tmp_path) -> None:
    """The first proposed fix, planted as the negative: resolving the path asks a DIFFERENT question.

    The failing shape exactly: a linked worktree whose `data/cache` is a symlink into the main
    checkout. The resolved path lies in another working tree of the same repository, and git refuses
    it with 128 a second time. Nothing here depends on git's wording beyond the code being 128.
    """
    root = _repo(tmp_path / "main")
    (root / "data" / "cache" / "entex").mkdir(parents=True)
    (root / "data" / "cache" / "entex" / "x.csv").write_text("x\n")
    worktree = tmp_path / "wt"
    added = _git(root, "worktree", "add", "-q", "--detach", str(worktree), "main")
    assert added.returncode == 0, added.stderr
    (worktree / "data").mkdir(exist_ok=True)
    (worktree / "data" / "cache").symlink_to(root / "data" / "cache")

    as_given = _git(worktree, "check-ignore", "-q", "--no-index", "data/cache/entex/x.csv")
    assert as_given.returncode == 128, as_given
    resolved = (worktree / "data" / "cache" / "entex" / "x.csv").resolve()
    as_resolved = _git(worktree, "check-ignore", "-q", "--no-index", str(resolved))
    assert as_resolved.returncode == 128, as_resolved
    # ... while the fix's question, asked at the symlink, is answerable
    at_boundary = _git(worktree, "check-ignore", "-q", "--no-index", "data/cache")
    assert at_boundary.returncode == 0, at_boundary


def test_a_store_missing_behind_a_symlink_skips_by_name_with_the_fetch_command(tmp_path, monkeypatch) -> None:
    """Not a pass and not an error: a SKIP whose reason names the store and how to get it."""
    root = _repo(tmp_path / "repo")
    (root / "data").mkdir()
    (root / "data" / "cache").symlink_to(root / "gone")  # a dangling link: the store is not here
    monkeypatch.chdir(root)
    local_data.ignored_by_name.cache_clear()

    absent = local_data.missing((FAILING_PATH,))
    assert absent == [FAILING_PATH]
    reason = local_data.skip_reason(absent, "scripts/entex_feasibility.py writes it (AG_METADATA)")
    assert "NOT RUN HERE, not passed" in reason
    assert FAILING_PATH in reason
    assert "scripts/entex_feasibility.py writes it (AG_METADATA)" in reason
    local_data.ignored_by_name.cache_clear()


def test_a_path_git_does_not_ignore_still_raises(tmp_path, monkeypatch) -> None:
    """The rule the fix must not weaken: a marker naming a committed path is wrong and says so."""
    root = _repo(tmp_path / "repo")
    monkeypatch.chdir(root)
    local_data.ignored_by_name.cache_clear()
    with pytest.raises(local_data.NotMachineLocalError, match="git does not ignore"):
        local_data.check_is_machine_local(".gitignore")
    local_data.ignored_by_name.cache_clear()


def test_a_symlink_that_is_not_an_ignored_name_is_refused_and_not_guessed(tmp_path, monkeypatch) -> None:
    """Truncating to the symlink answers only when that name IS ignored. Otherwise it refuses."""
    root = _repo(tmp_path / "repo")
    (root / "elsewhere").mkdir()
    (root / "notignored").symlink_to(root / "elsewhere")
    monkeypatch.chdir(root)
    local_data.ignored_by_name.cache_clear()
    with pytest.raises(local_data.NotMachineLocalError, match="stops at the symbolic link"):
        local_data.check_is_machine_local("notignored/deep/file.csv")
    local_data.ignored_by_name.cache_clear()


#: Every `needs_local_data` call site in the suite, read as text. Collecting the real markers would
#: mean importing every test module; the fetch command is written at the call site, so the call site
#: is what this reads.
_MARKER_CALL = re.compile(r"needs_local_data\(\s*(.*?)\)\s*\n(?:def |@|needs_|\s*def )", re.S)


def test_every_needs_local_data_marker_carries_the_command_that_fetches_the_store() -> None:
    """A skip reason without a `how` names the store but not the way back. None may be without one."""
    without = []
    for path in sorted(Path("tests").glob("test_*.py")):
        text = path.read_text()
        for match in re.finditer(r"needs_local_data\(", text):
            start = match.end()
            depth, i = 1, start
            while i < len(text) and depth:
                depth += (text[i] == "(") - (text[i] == ")")
                i += 1
            if "how=" not in text[start : i - 1]:
                without.append(f"{path}:{text.count(chr(10), 0, match.start()) + 1}")
    assert without == [], (
        f"these needs_local_data markers give no `how=`, so their skip reason cannot name the command "
        f"that fetches the store: {without}"
    )
