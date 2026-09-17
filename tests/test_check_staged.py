# SPDX-License-Identifier: AGPL-3.0-or-later
"""The stale-base guard, against a repository built to contain the accident it is for."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "check_staged.py"

BASE = """# Grammar

| field | type | meaning |
|---|---|---|
| `fraction` | `number` | a share of what is left when it runs |
"""

PEER_ADDITION = """# Grammar

| field | type | meaning |
|---|---|---|
| `fraction` | `number` | a share of what is left when it runs |
| `share` | `number` | a share of the whole at this decision point, order-free |
"""

OWN_EDIT = """# Grammar

| field | type | meaning |
|---|---|---|
| `fraction` | `number` | a share of what is left when it runs |
| `share` | `number` | a share of the whole at this decision point, order-free |
| `order` | `list` | the sequence a cluster is read in |
"""


def git(repo: Path, *args: str, index: str | None = None) -> str:
    env = dict(os.environ)
    env.pop("GIT_INDEX_FILE", None)
    if index:
        env["GIT_INDEX_FILE"] = index
    done = subprocess.run(["git", *args], cwd=repo, capture_output=True, text=True, env=env, check=True)
    return done.stdout


def run_check(repo: Path, index: str, *extra: str) -> subprocess.CompletedProcess[str]:
    env = dict(os.environ)
    env["GIT_INDEX_FILE"] = index
    return subprocess.run(
        [sys.executable, str(SCRIPT), *extra],
        cwd=repo,
        capture_output=True,
        text=True,
        env=env,
    )


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    """A shared file, a peer's commit adding a feature row to it, and a stale copy to stage."""
    repo = tmp_path / "checkout"
    repo.mkdir()
    git(repo, "init", "-q", "-b", "main")
    git(repo, "config", "user.email", "lane@example.test")
    git(repo, "config", "user.name", "A Lane")
    shared = repo / "GRAMMAR.md"
    shared.write_text(BASE)
    git(repo, "add", "GRAMMAR.md")
    git(repo, "commit", "-qm", "the grammar table")
    shared.write_text(PEER_ADDITION)
    git(repo, "add", "GRAMMAR.md")
    git(repo, "commit", "-qm", "share: the same partition, stated as its numbers")
    return repo


def private_index(repo: Path, tmp_path: Path) -> str:
    index = str(tmp_path / "idx")
    git(repo, "read-tree", "HEAD", index=index)
    return index


def stage_blob(repo: Path, index: str, path: str, content: str) -> None:
    """Put content into the index without touching the working copy, as a staging tool would."""
    tmp = repo / ".stage-tmp"
    tmp.write_text(content)
    sha = git(repo, "hash-object", "-w", str(tmp), index=index).strip()
    tmp.unlink()
    git(repo, "update-index", "--add", "--cacheinfo", f"100644,{sha},{path}", index=index)


def test_stale_copy_is_refused_and_names_the_peer(repo: Path, tmp_path: Path) -> None:
    """The accident of 2026-09-15: staging the copy you read before the peer's commit landed."""
    index = private_index(repo, tmp_path)
    stage_blob(repo, index, "GRAMMAR.md", BASE)  # the pre-share copy, staged whole

    done = run_check(repo, index)

    assert done.returncode == 2, done.stdout
    assert "REFUSED" in done.stdout
    assert "share: the same partition" in done.stdout
    assert "`share`" in done.stdout


def test_an_addition_on_top_of_the_peer_passes(repo: Path, tmp_path: Path) -> None:
    """Adding your own row to the current copy removes nobody's line and must not be refused."""
    index = private_index(repo, tmp_path)
    stage_blob(repo, index, "GRAMMAR.md", OWN_EDIT)

    done = run_check(repo, index)

    assert done.returncode == 0, done.stdout
    assert "takes nothing back out" in done.stdout


def test_force_reports_without_refusing(repo: Path, tmp_path: Path) -> None:
    """A deliberate removal is sometimes right; it costs a flag and still prints what it takes out."""
    index = private_index(repo, tmp_path)
    stage_blob(repo, index, "GRAMMAR.md", BASE)

    done = run_check(repo, index, "--force")

    assert done.returncode == 0, done.stdout
    assert "REFUSED" in done.stdout  # seen, not hidden


def test_a_staged_deletion_of_a_tracked_file_is_reported(repo: Path, tmp_path: Path) -> None:
    """Losing a peer's whole file is the unrecoverable case, so it is never silent."""
    index = private_index(repo, tmp_path)
    git(repo, "update-index", "--force-remove", "GRAMMAR.md", index=index)

    done = run_check(repo, index)

    assert done.returncode == 2
    assert "staged as DELETED" in done.stdout


def test_an_empty_commit_is_refused(repo: Path, tmp_path: Path) -> None:
    """A push that looks like it failed may have landed; committing again produces an empty commit."""
    index = private_index(repo, tmp_path)

    done = run_check(repo, index)

    assert done.returncode == 2
    assert "would be empty" in done.stdout


def test_a_new_file_reverts_nothing(repo: Path, tmp_path: Path) -> None:
    """A path HEAD does not have cannot take anything back out, whatever it contains."""
    index = private_index(repo, tmp_path)
    stage_blob(repo, index, "NOTES.md", BASE)

    done = run_check(repo, index)

    assert done.returncode == 0, done.stdout


def test_old_removals_pass_and_the_window_decides(repo: Path, tmp_path: Path) -> None:
    """The guard is a stale-base detector, not a merge policeman: outside the window it says nothing."""
    index = private_index(repo, tmp_path)
    stage_blob(repo, index, "GRAMMAR.md", BASE)

    refused = run_check(repo, index)
    passed = run_check(repo, index, "--since", "2099-01-01")  # a window with no commit in it

    assert refused.returncode == 2
    assert passed.returncode == 0, passed.stdout
