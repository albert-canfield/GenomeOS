# SPDX-License-Identifier: AGPL-3.0-or-later
"""A verdict worktree gets the git-ignored stores READ-ONLY, and still hashes to the committed tree.

Two claims, and both were load-bearing guesses until they were run.

(1) A write into a linked store is REFUSED. Until 2026-10-02 `scripts/pre-push.sh` symlinked the three
    stores into the worktree every push is verified in, and a symlink carries its TARGET's mode: the
    stores' directories are 0755, so every push verdict this project has taken was produced in a
    worktree whose input stores were writable. Measured on the real tree before the fix:
    `open("<worktree>/data/reference/HG002_chr1.vcf.gz", "r+b")` through the link SUCCEEDED.

(2) TREE IDENTITY SURVIVES. A verdict must name the committed tree, `tree_begin == tree_end ==` it, and
    a clone is a real DIRECTORY where there used to be a symlink -- so this is not self-evident.
    `.gitignore` spells the three stores without a trailing slash, which matches a directory as well as
    a symlink, and the comment there records that the spelling was chosen for this. Measured on the real
    tree: e49a3467eeae3aa183c2b9426cdfadcb4e2c21e5 before linking and the same after, equal to
    `git rev-parse HEAD^{tree}`.

The tests below are HERMETIC -- a planted repository with planted stores -- so they run anywhere and in
a second rather than cloning 11.2 GB. The real-tree figures above are the measurement; these are the
mechanism, and they fail if either claim stops holding.
"""

from __future__ import annotations

import importlib.util
import os
import subprocess
from pathlib import Path

import pytest

_spec = importlib.util.spec_from_file_location("lsro", Path("scripts/link_stores_read_only.py"))
lsro = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(lsro)


@pytest.fixture(autouse=True)
def _no_inherited_git_environment(monkeypatch: pytest.MonkeyPatch) -> None:
    for name in ("GIT_INDEX_FILE", "GIT_DIR", "GIT_WORK_TREE", "GIT_OBJECT_DIRECTORY"):
        monkeypatch.delenv(name, raising=False)


def _git(repo: Path, *args: str) -> str:
    r = subprocess.run(
        ["git", "-C", str(repo), *args], check=True, capture_output=True, text=True, env=dict(os.environ)
    )
    return r.stdout.strip()


@pytest.fixture
def planted(tmp_path: Path) -> dict[str, Path]:
    """A repository with the three store names git-ignored WITHOUT a trailing slash, as this one spells
    them, holding a file in each; and a worktree of its one commit with no stores in it."""
    repo = tmp_path / "planted"
    (repo / "data").mkdir(parents=True)
    (repo / "tracked.txt").write_text("in the commit\n")
    # the spelling under test: no trailing slash, so it matches a real directory as well as a symlink
    (repo / ".gitignore").write_text("data/reference\ndata/knowledge\ndata/cache\n")
    for name in lsro.STORES:
        d = repo / "data" / name
        d.mkdir()
        (d / f"{name}_table.tsv").write_text(f"{name}\n")
        (d / "nested").mkdir()
        (d / "nested" / "deep.tsv").write_text("deep\n")
    _git(repo, "init", "-b", "main")
    _git(repo, "config", "user.email", "lane@example.invalid")
    _git(repo, "config", "user.name", "planted")
    _git(repo, "config", "commit.gpgsign", "false")
    _git(repo, "add", "--", "tracked.txt", ".gitignore")
    _git(repo, "commit", "-m", "planted")
    wt = tmp_path / "worktree"
    _git(repo, "worktree", "add", "-q", "--detach", str(wt), "HEAD")
    return {"repo": repo, "wt": wt, "tree": _git(repo, "rev-parse", "HEAD^{tree}")}


@pytest.fixture
def linked(planted: dict[str, Path]):
    """The stores linked read-only, and unlocked again however the test ends, so tmp_path can be removed."""
    out = lsro.link(planted["repo"], planted["wt"])
    try:
        yield {**planted, "out": out}
    finally:
        lsro.unlink(planted["wt"])


# --- the stores arrive, and arrive read-only --------------------------------------------------------


def test_every_store_is_linked_and_its_bytes_are_there(linked: dict) -> None:
    for name in lsro.STORES:
        d = linked["wt"] / "data" / name
        assert d.is_dir() and not d.is_symlink(), f"{name} must be a clone, not a symlink"
        assert (d / f"{name}_table.tsv").read_text() == f"{name}\n"
        assert (d / "nested" / "deep.tsv").read_text() == "deep\n"
    assert linked["out"]["ok"] is True and linked["out"]["failures"] == []


def test_a_write_into_a_linked_store_is_refused_by_the_kernel(linked: dict) -> None:
    """Stronger than the audit hook in scripts/rebuild_write_guard.py, and complementary to it: a mode
    is enforced against a C library -- htslib, pysam -- and an audit hook is not."""
    d = linked["wt"] / "data" / "reference"
    with pytest.raises(OSError), open(d / "new_file.tmp", "w"):
        pass
    with pytest.raises(OSError), open(d / "reference_table.tsv", "r+b"):
        pass
    assert not (linked["repo"] / "data" / "reference" / "new_file.tmp").exists()


def test_the_probe_is_part_of_linking_so_a_writable_store_is_a_named_failure(linked: dict) -> None:
    """Linking does not assume the lock took: it tries a write and reports whether it was refused."""
    for store in linked["out"]["stores"]:
        assert store["a_write_into_it_is_refused"] is True, store


def test_the_machine_s_own_copy_keeps_its_write_bits(linked: dict) -> None:
    """The reason for a clone rather than a hard link or a chmod of the original: other sessions share
    the real store, and taking its write bits off would be taking them off for everyone."""
    real = linked["repo"] / "data" / "reference" / "reference_table.tsv"
    assert real.stat().st_mode & 0o200, "the real store must stay writable"
    real.write_text("the machine can still write its own copy\n")
    assert (linked["wt"] / "data" / "reference" / "reference_table.tsv").read_text() == "reference\n", (
        "a clone shares bytes until one side writes; the worktree must keep what it was given"
    )


# --- tree identity, which is what a verdict names ---------------------------------------------------


def test_the_worktree_still_hashes_to_the_committed_tree(planted: dict[str, Path]) -> None:
    """The claim pre-push.sh now rests on. A clone is a real directory where a symlink was, so this is
    checked rather than assumed; if .gitignore ever grew a trailing slash, this fails."""
    before = lsro.tree_hash(planted["wt"])
    assert before == planted["tree"], f"{before} != {planted['tree']}"
    out = lsro.link(planted["repo"], planted["wt"], check_tree=True)
    try:
        assert out["tree_before"] == planted["tree"]
        assert out["tree_after"] == planted["tree"]
        assert out["tree_unchanged"] is True
    finally:
        lsro.unlink(planted["wt"])


def test_a_trailing_slash_in_the_ignore_rule_would_break_it(planted: dict[str, Path]) -> None:
    """The counterfactual, run rather than argued. `data/reference/` matches a DIRECTORY only in the
    sense git means, and the .gitignore comment in this repository records that the three stores are
    spelled without a slash for the pre-push worktree's sake. Here the rules are removed entirely, which
    is the limit of that case: the cloned stores then count as untracked content and the hash MOVES."""
    (planted["repo"] / ".gitignore").write_text("# nothing ignored\n")
    _git(planted["repo"], "add", "--", ".gitignore")
    _git(planted["repo"], "commit", "-m", "stop ignoring the stores")
    wt2 = planted["wt"].parent / "worktree2"
    _git(planted["repo"], "worktree", "add", "-q", "--detach", str(wt2), "HEAD")
    try:
        before = lsro.tree_hash(wt2)
        lsro.link(planted["repo"], wt2)
        after = lsro.tree_hash(wt2)
        assert before != after, (
            "with the stores unignored the hash MUST move; if it does not, the identity check in "
            "pre-push.sh proves nothing"
        )
    finally:
        lsro.unlink(wt2)
        _git(planted["repo"], "worktree", "remove", "--force", str(wt2))


# --- unlocking, so a worktree can be taken away -----------------------------------------------------


def test_unlocking_puts_the_write_bits_back_so_the_worktree_can_be_removed(
    planted: dict[str, Path],
) -> None:
    """Without this a read-only directory refuses `rm -rf` and an interrupted push leaves a tree nothing
    can take away. pre-push.sh does the same with a shell chmod, which needs no venv inside a trap."""
    lsro.link(planted["repo"], planted["wt"])
    out = lsro.unlink(planted["wt"])
    assert {u["store"] for u in out["unlocked"]} == set(lsro.STORES)
    (planted["wt"] / "data" / "reference" / "now_writable.tmp").write_text("x")
    _git(planted["repo"], "worktree", "remove", "--force", str(planted["wt"]))
    assert not planted["wt"].exists()


# --- the clone is reported as a clone, and a real copy would be reported as one ----------------------


def test_the_report_says_whether_apfs_cloned_or_the_bytes_were_copied(linked: dict) -> None:
    """A clone costs nothing and a copy costs the store's size; the two read the same in a report
    otherwise, and a rebuild that filled a disk would have had no way to say why."""
    for store in linked["out"]["stores"]:
        assert store["cloned"] is not store["copied_not_cloned"]
        assert isinstance(store["seconds"], float)
        assert store["entries_locked"] >= 3, store
    assert "free_bytes_spent" in linked["out"]


def test_linking_twice_leaves_the_first_link_alone(linked: dict) -> None:
    """Idempotent, because pre-push may be re-run over a worktree it already prepared."""
    again = lsro.link(linked["repo"], linked["wt"])
    assert all(s["linked"] is False for s in again["stores"])
    assert all(s["why"] == "already in the worktree" for s in again["stores"])
    assert again["ok"] is True


def test_a_store_not_on_this_machine_is_named_rather_than_failing(tmp_path: Path) -> None:
    out = lsro.link(tmp_path / "no-such-root", tmp_path / "wt")
    assert [s["why"] for s in out["stores"]] == ["not on this machine"] * len(lsro.STORES)
    assert out["ok"] is True, "a store that is simply absent is not a failure to link"
