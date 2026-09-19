# SPDX-License-Identifier: AGPL-3.0-or-later
"""The lost-commit audit, against a repository built to contain the loss it looks for."""

from __future__ import annotations

import importlib.util
import os
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("check_history", ROOT / "scripts" / "check_history.py")
assert spec and spec.loader
ch = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ch)


def git(repo: Path, *args: str) -> str:
    env = dict(os.environ)
    env.pop("GIT_INDEX_FILE", None)
    done = subprocess.run(["git", *args], cwd=repo, capture_output=True, text=True, env=env, check=True)
    return done.stdout.strip()


@pytest.fixture
def repo(tmp_path: Path, monkeypatch) -> Path:
    """A branch with one ordinary commit, then a ref moved sideways over a peer's commit."""
    repo = tmp_path / "checkout"
    repo.mkdir()
    git(repo, "init", "-q", "-b", "dev")
    git(repo, "config", "user.email", "lane@example.test")
    git(repo, "config", "user.name", "A Lane")
    (repo / "a.txt").write_text("one\n")
    git(repo, "add", "a.txt")
    git(repo, "commit", "-qm", "base")
    base = git(repo, "rev-parse", "HEAD")

    # a peer's commit lands
    (repo / "peer.txt").write_text("peer work\n")
    git(repo, "add", "peer.txt")
    git(repo, "commit", "-qm", "the peer's commit")

    # another session, which read HEAD at `base`, moves the ref over it with no old value
    mine = git(repo, "commit-tree", f"{base}^{{tree}}", "-p", base, "-m", "my commit, built on a stale read")
    git(repo, "update-ref", "refs/heads/dev", mine)

    monkeypatch.chdir(repo)
    return repo


def test_a_ref_moved_over_a_peers_commit_is_found(repo: Path) -> None:
    """The failure with no conflict and no error: the peer's commit is simply not in the history."""
    rows = ch.sideways("dev")

    assert len(rows) == 1
    assert rows[0]["looks_like"] == "a ref moved over a commit"
    assert rows[0]["old_subject"] == "the peer's commit"
    assert rows[0]["new_subject"] == "my commit, built on a stale read"


def test_the_lost_commit_is_reported_as_not_re_landed(repo: Path) -> None:
    """`relanded_as` is what separates "this happened" from "this was recovered"."""
    assert ch.sideways("dev")[0]["relanded_as"] is None


def test_re_landing_the_same_subject_is_recognised(repo: Path) -> None:
    """All six on this project's dev re-landed; the audit must be able to say so."""
    (repo / "peer.txt").write_text("peer work\n")
    git(repo, "add", "peer.txt")
    git(repo, "commit", "-qm", "the peer's commit")

    assert ch.sideways("dev")[0]["relanded_as"] is not None


def test_a_fast_forward_history_reports_nothing(tmp_path: Path, monkeypatch) -> None:
    """The common case must be silent, or nobody will run it twice."""
    repo = tmp_path / "clean"
    repo.mkdir()
    git(repo, "init", "-q", "-b", "dev")
    git(repo, "config", "user.email", "l@example.test")
    git(repo, "config", "user.name", "L")
    for i in range(3):
        (repo / f"{i}.txt").write_text(f"{i}\n")
        git(repo, "add", ".")
        git(repo, "commit", "-qm", f"commit {i}")
    monkeypatch.chdir(repo)

    assert ch.sideways("dev") == []


def test_an_amend_is_labelled_apart_from_a_loss(tmp_path: Path, monkeypatch) -> None:
    """An amend moves the ref sideways too, and calling it a lost commit would cry wolf."""
    repo = tmp_path / "amended"
    repo.mkdir()
    git(repo, "init", "-q", "-b", "dev")
    git(repo, "config", "user.email", "l@example.test")
    git(repo, "config", "user.name", "L")
    (repo / "a.txt").write_text("one\n")
    git(repo, "add", "a.txt")
    git(repo, "commit", "-qm", "the only subject")
    (repo / "a.txt").write_text("two\n")
    git(repo, "add", "a.txt")
    git(repo, "commit", "-q", "--amend", "-m", "the only subject")
    monkeypatch.chdir(repo)

    rows = ch.sideways("dev")

    assert len(rows) == 1
    assert rows[0]["looks_like"] == "amend or rewrite"
