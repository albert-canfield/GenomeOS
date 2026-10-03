# SPDX-License-Identifier: AGPL-3.0-or-later
"""A worktree holding read-only store clones is removed by git, never by `rm`, and a worktree holding
the REAL store is refused with the path named.

WHY THIS FILE EXISTS. `scripts/link_stores_read_only.py` clones the three git-ignored stores into a
verdict worktree at 0444/0555 so a verification run cannot modify the inputs it verifies against. On
2026-10-02 a lane then had to remove such a worktree, found `git worktree remove --force` refused, and
fell back to `chmod -R u+w` over the whole worktree followed by `rm -rf`. The rule in force since is
that no worktree is removed by `rm`, and `scripts/remove_worktree.py` is the tool behind that rule.

Every test here is HERMETIC -- a planted repository with planted stores under `tmp_path`, which pytest
cleans up -- so nothing in this file can reach the machine's real 11.2 GB of stores, and nothing in it
calls `rm`. The two stale worktrees in the real checkout are not touched by any test here.

THE FOUR PLANTS, and each must FAIL if the guard it tests is taken out:

(a) a worktree whose recorded data path IS the real store, by two routes -- a record pointing straight
    at the store, and a symlink to it, which is what `pre-push.sh` did until 2026-10-02 and which is
    dangerous precisely because a symlink carries its TARGET's mode. Both must be REFUSED with the path
    named, and the store's inode AND mode are asserted before and after.
(b) a worktree with genuine `cp -c` clones at 0444 is removed successfully, and the real stores' inodes
    and modes are asserted unchanged afterwards.
(c) a path absent from `git worktree list` is refused. An orphan left by a half-done removal lands here.
(d) the MEASUREMENT the tool's existence rests on: plain `git worktree remove`, with and without
    `--force`, on a worktree whose clones are at 0444. Measured on 2026-10-03 with git 2.50.0:
    exit 255, `error: failed to delete '<wt>': Permission denied`, and the worktree survives. The
    second finding, which was not expected: the failed removal is NOT atomic -- it had already deleted
    the tracked files and the `.git` file and DEREGISTERED the worktree, so `--force` then answers
    `fatal: '<wt>' is not a working tree` and the directory on disk is an orphan git can no longer
    touch by any flag. That is the dead end that leaves a lane with no git route at all, and it is why
    `remove_worktree.py` unlocks BEFORE it calls git rather than after a first attempt fails.
"""

from __future__ import annotations

import importlib.util
import json
import os
import stat
import subprocess
from pathlib import Path

import pytest

_rspec = importlib.util.spec_from_file_location("rmwt", Path("scripts/remove_worktree.py"))
rmwt = importlib.util.module_from_spec(_rspec)
_rspec.loader.exec_module(rmwt)
lsro = rmwt.lsro


@pytest.fixture(autouse=True)
def _no_inherited_git_environment(monkeypatch: pytest.MonkeyPatch) -> None:
    """A planted repository must not be read through this checkout's shared index."""
    for name in ("GIT_INDEX_FILE", "GIT_DIR", "GIT_WORK_TREE", "GIT_OBJECT_DIRECTORY"):
        monkeypatch.delenv(name, raising=False)


def _git(repo: Path, *args: str, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["git", "-C", str(repo), *args], check=check, capture_output=True, text=True, env=dict(os.environ)
    )


def _state(p: Path) -> tuple[int, int, str]:
    """(device, inode, mode) -- the reading every plant asserts before and after."""
    st = p.stat()
    return (st.st_dev, st.st_ino, oct(stat.S_IMODE(st.st_mode)))


@pytest.fixture
def planted(tmp_path: Path) -> dict[str, Path]:
    """A repository with the three store names git-ignored as this one spells them, each holding a file
    in a subdirectory, and ONE registered worktree of its single commit with no stores in it yet."""
    repo = tmp_path / "planted"
    (repo / "data").mkdir(parents=True)
    (repo / "tracked.txt").write_text("in the commit\n")
    (repo / ".gitignore").write_text("".join(f"data/{n}\n" for n in lsro.STORES))
    for name in lsro.STORES:
        (repo / "data" / name / "sub").mkdir(parents=True)
        (repo / "data" / name / "sub" / "f.bin").write_bytes(b"store bytes\n")
    _git(repo, "init", "-q", ".")
    _git(repo, "config", "user.email", "lane@example.invalid")
    _git(repo, "config", "user.name", "lane")
    _git(repo, "add", "tracked.txt", ".gitignore")
    _git(repo, "commit", "-qm", "the one commit")
    worktree = tmp_path / "wt"
    _git(repo, "worktree", "add", "-q", "--detach", str(worktree), "HEAD")
    return {"repo": repo, "worktree": worktree}


@pytest.fixture
def linked(planted: dict[str, Path]) -> dict[str, Path]:
    """`planted`, with the stores cloned in read-only by the real linker and RECORDED by it."""
    r = lsro.link(planted["repo"], planted["worktree"])
    assert r["ok"], r
    assert Path(r["manifest"]).is_file(), r
    return planted


# ---------------------------------------------------------------- the linker's record


def test_the_linker_records_what_it_cloned(linked: dict[str, Path]) -> None:
    """The record is the remover's only source of what may be made writable, so its shape is a test.
    Before 2026-10-03 the linker printed what it cloned and recorded nothing, and a remover with no
    record has only a recursive chmod available -- which is how the incident began."""
    record = lsro.read_manifest(linked["worktree"].resolve(), linked["repo"])
    assert sorted(c["store"] for c in record["clones"]) == sorted(lsro.STORES)
    for c in record["clones"]:
        assert Path(c["path"]).is_dir()
        assert c["inode"] != c["source_inode"], c
        assert Path(c["path"]).stat().st_ino == c["inode"]


def test_the_record_is_outside_the_worktree_so_tree_identity_survives(linked: dict[str, Path]) -> None:
    """A record written INSIDE the worktree would be untracked content, and `tree_before ==
    tree_after` is what `link_stores_read_only` exists to preserve."""
    manifest = lsro.manifest_path(linked["worktree"].resolve(), linked["repo"])
    assert not rmwt._inside(manifest, linked["worktree"])
    assert _git(linked["worktree"], "status", "--porcelain").stdout.strip() == ""


def test_the_record_is_keyed_per_worktree(planted: dict[str, Path], tmp_path: Path) -> None:
    """Two worktrees must never share a record, or one removal reads the other's paths."""
    other = tmp_path / "wt2"
    _git(planted["repo"], "worktree", "add", "-q", "--detach", str(other), "HEAD")
    a = lsro.manifest_path(planted["worktree"].resolve(), planted["repo"])
    b = lsro.manifest_path(other.resolve(), planted["repo"])
    assert a != b


def test_a_second_link_run_does_not_drop_the_first_runs_record(linked: dict[str, Path]) -> None:
    """`link` leaves a store it finds already there alone, so a second run records nothing new; if it
    overwrote the record, the remover would then unlock nothing and git would refuse."""
    again = lsro.link(linked["repo"], linked["worktree"].resolve())
    assert all(not s["linked"] for s in again["stores"]), again
    record = lsro.read_manifest(linked["worktree"].resolve(), linked["repo"])
    assert sorted(c["store"] for c in record["clones"]) == sorted(lsro.STORES)


# ---------------------------------------------------------------- plant (a): the real store


def test_a_record_pointing_at_the_real_store_is_refused_and_the_store_is_untouched(
    linked: dict[str, Path],
) -> None:
    """The plant: the worktree's recorded data path IS the real store, same inode. Removing it would
    remove the machine's only copy, so it must be refused with the path named."""
    repo, worktree = linked["repo"], linked["worktree"].resolve()
    real = repo / "data" / "reference"
    before = {n: _state(repo / "data" / n) for n in lsro.STORES}

    manifest = lsro.manifest_path(worktree, repo)
    record = json.loads(manifest.read_text())
    for c in record["clones"]:
        if c["store"] == "reference":
            st = real.stat()
            c["path"], c["inode"], c["device"] = str(real), st.st_ino, st.st_dev
    manifest.write_text(json.dumps(record))

    with pytest.raises(rmwt.RefusedError) as e:
        rmwt.remove(repo, worktree)
    assert str(real) in str(e.value)
    # The reading is required to be the INODE one specifically. This assertion read
    # `"SAME INODE" in ... or "inside" in ...` until a mutation run showed that taking the inode check
    # out left the file green: the containment check behind it fired instead and the disjunction
    # accepted it. A test that any one of two guards can satisfy measures neither.
    assert "SAME INODE" in str(e.value), str(e.value)
    assert f"inode {real.stat().st_ino}" in str(e.value)

    assert {n: _state(repo / "data" / n) for n in lsro.STORES} == before
    assert (real / "sub" / "f.bin").read_bytes() == b"store bytes\n"
    assert worktree.is_dir(), "a refusal must not have removed anything"


def test_a_path_inside_the_worktree_that_is_the_real_store_is_refused(planted: dict[str, Path]) -> None:
    """The same hazard with the record telling the truth: the recorded path really is inside the
    worktree, and the directory sitting there really is the machine's store, because `data/reference`
    in the checkout is a symlink to it. Nothing about the record is forged here, so this plant reaches
    the inode check through the ordinary path and not past the locality checks."""
    repo, worktree = planted["repo"], planted["worktree"].resolve()
    # the store itself now LIVES in the worktree, and the checkout reaches it through a link
    (worktree / "data").mkdir(parents=True, exist_ok=True)
    (repo / "data" / "reference" / "sub" / "f.bin").rename(worktree / "data" / "f.bin")
    (repo / "data" / "reference" / "sub").rmdir()
    (repo / "data" / "reference").rmdir()
    store_now = worktree / "data" / "reference"
    store_now.mkdir()
    (worktree / "data" / "f.bin").rename(store_now / "f.bin")
    (repo / "data" / "reference").symlink_to(store_now)

    st = store_now.stat()
    lsro.write_manifest(
        worktree,
        repo,
        [
            {
                "store": "reference",
                "path": str(store_now),
                "source": str(repo / "data" / "reference"),
                "inode": st.st_ino,
                "device": st.st_dev,
                "source_inode": st.st_ino,
                "source_device": st.st_dev,
                "cloned": False,
                "entries_locked": 0,
            }
        ],
    )
    before = _state(repo / "data" / "reference")

    with pytest.raises(rmwt.RefusedError) as e:
        rmwt.remove(repo, worktree)
    assert "SAME INODE" in str(e.value), str(e.value)
    assert str(store_now) in str(e.value)
    assert _state(repo / "data" / "reference") == before
    assert (store_now / "f.bin").read_bytes() == b"store bytes\n"
    assert worktree.is_dir()


def test_a_symlink_to_the_real_store_is_refused_and_the_store_is_untouched(
    planted: dict[str, Path],
) -> None:
    """The other route to the same hazard, and the arrangement `pre-push.sh` really used until
    2026-10-02: the worktree's store is a SYMLINK to the real one. A symlink carries its target's
    mode, so restoring write through it reaches the store itself."""
    repo, worktree = planted["repo"], planted["worktree"].resolve()
    before = {n: _state(repo / "data" / n) for n in lsro.STORES}

    (worktree / "data").mkdir(parents=True, exist_ok=True)
    (worktree / "data" / "reference").symlink_to(repo / "data" / "reference")
    st = (repo / "data" / "reference").stat()
    lsro.write_manifest(
        worktree,
        repo,
        [
            {
                "store": "reference",
                "path": str(worktree / "data" / "reference"),
                "source": str(repo / "data" / "reference"),
                "inode": st.st_ino,
                "device": st.st_dev,
                "source_inode": st.st_ino,
                "source_device": st.st_dev,
                "cloned": False,
                "entries_locked": 0,
            }
        ],
    )

    with pytest.raises(rmwt.RefusedError) as e:
        rmwt.remove(repo, worktree)
    assert str(worktree / "data" / "reference") in str(e.value)
    assert "SYMLINK" in str(e.value)

    assert {n: _state(repo / "data" / n) for n in lsro.STORES} == before
    assert (repo / "data" / "reference" / "sub" / "f.bin").read_bytes() == b"store bytes\n"


def test_the_refusal_comes_before_anything_is_made_writable(linked: dict[str, Path]) -> None:
    """A refusal that happened after the chmod would be no protection at all: the clones must still be
    at 0444 after the tool has said no."""
    repo, worktree = linked["repo"], linked["worktree"].resolve()
    manifest = lsro.manifest_path(worktree, repo)
    record = json.loads(manifest.read_text())
    # the LAST recorded store is made the dangerous one, so the refusal has to arrive after the tool
    # has already looked at the safe ones -- the order is what makes the "before anything" claim real
    st = (repo / "data" / "reference").stat()
    record["clones"][-1] |= {
        "path": str(repo / "data" / "reference"),
        "inode": st.st_ino,
        "device": st.st_dev,
    }
    manifest.write_text(json.dumps(record))
    with pytest.raises(rmwt.RefusedError):
        rmwt.remove(repo, worktree)
    for c in record["clones"][:-1]:
        f = Path(c["path"]) / "sub" / "f.bin"
        assert oct(stat.S_IMODE(f.stat().st_mode)) == "0o444", f


# ---------------------------------------------------------------- plant (b): genuine clones


def test_a_worktree_of_genuine_read_only_clones_is_removed(linked: dict[str, Path]) -> None:
    """The success case, and the whole point: 0444 clones, removed by git, stores unchanged."""
    repo, worktree = linked["repo"], linked["worktree"].resolve()
    before = {n: _state(repo / "data" / n) for n in lsro.STORES}
    assert oct(stat.S_IMODE((worktree / "data" / "reference" / "sub" / "f.bin").stat().st_mode)) == "0o444"

    r = rmwt.remove(repo, worktree)

    assert r["ok"], r
    assert r["removed"] is True, r
    assert r["git_exit"] == 0, r
    assert not worktree.exists()
    assert worktree not in rmwt.registered_worktrees(repo)
    after = {n: _state(repo / "data" / n) for n in lsro.STORES}
    assert after == before, (before, after)
    assert (repo / "data" / "reference" / "sub" / "f.bin").read_bytes() == b"store bytes\n"
    assert r["stores_unchanged"] is True
    assert not lsro.manifest_path(worktree, repo).exists(), "the record outlived the worktree"


def test_the_tool_unlocks_only_the_recorded_paths(linked: dict[str, Path]) -> None:
    """A file in the worktree that NO record mentions keeps its mode, while the recorded clones become
    writable: the tool restores write on the record and not on a tree, and a recursive chmod is what
    reaches a store nobody recorded.

    The unlock has to be OBSERVED after it ran, which means the worktree has to survive the call. So
    git is made to refuse: a modified tracked file and no `--force` is refused with
    `contains modified or untracked files`, measured 2026-10-03, and git refuses that one BEFORE
    deleting anything (unlike the 0444 refusal in plant (d), which deletes half the worktree first).
    An earlier version of this test used `check_only=True` and therefore never ran the unlock at all:
    a mutation that widened the unlock to the whole worktree left it green.
    """
    repo, worktree = linked["repo"], linked["worktree"].resolve()
    unrecorded = worktree / "data" / "unrecorded.bin"
    unrecorded.write_bytes(b"not in the record\n")
    os.chmod(unrecorded, 0o444)
    (worktree / "tracked.txt").write_text("a lane changed this\n")

    r = rmwt.remove(repo, worktree, force=False)

    assert r["removed"] is False and r["ok"] is False, r
    assert {i["store"] for i in r["plan"]} == set(lsro.STORES)
    assert {u["path"] for u in r["unlocked"]} == {str(worktree / "data" / n) for n in lsro.STORES}, r[
        "unlocked"
    ]
    # the recorded clones are writable now, and the file nobody recorded is still 0444
    for name in lsro.STORES:
        f = worktree / "data" / name / "sub" / "f.bin"
        assert oct(stat.S_IMODE(f.stat().st_mode)) == "0o644", f
    assert oct(stat.S_IMODE(unrecorded.stat().st_mode)) == "0o444"


def test_a_git_refusal_is_reported_verbatim_and_not_escalated(linked: dict[str, Path]) -> None:
    """Requirement 3: if git still refuses, report its exact message and STOP. The escalation on offer
    from here is `rm -rf`, which is the one that can destroy a store, so the tool must have no path to
    it -- and the worktree must still be standing afterwards for its owner to deal with."""
    repo, worktree = linked["repo"], linked["worktree"].resolve()
    before = {n: _state(repo / "data" / n) for n in lsro.STORES}
    (worktree / "tracked.txt").write_text("a lane changed this\n")

    r = rmwt.remove(repo, worktree, force=False)

    assert r["git_exit"] != 0
    assert r["refused_by_git"] == r["git_stderr"]
    assert "contains modified or untracked files" in r["refused_by_git"], r
    assert worktree.is_dir(), "a git refusal must leave the worktree for its owner"
    assert lsro.manifest_path(worktree, repo).is_file(), "the record outlives a refused removal"
    assert {n: _state(repo / "data" / n) for n in lsro.STORES} == before


def test_force_is_passed_to_git_and_clears_untracked_scratch(linked: dict[str, Path]) -> None:
    """Requirement 3's other half: `--force` means git's `--force`, for untracked scratch files, and
    never a stronger tool."""
    repo, worktree = linked["repo"], linked["worktree"].resolve()
    before = {n: _state(repo / "data" / n) for n in lsro.STORES}
    (worktree / "scratch.log").write_text("untracked scratch\n")
    r = rmwt.remove(repo, worktree, force=True)
    assert r["git_command"].split()[:4] == ["git", "worktree", "remove", "--force"]
    assert r["removed"] and r["ok"], r
    assert not worktree.exists()
    assert {n: _state(repo / "data" / n) for n in lsro.STORES} == before


def test_a_record_whose_inode_no_longer_matches_is_refused(linked: dict[str, Path]) -> None:
    """The record has to describe what is actually there. If the path was replaced since the linker
    wrote the record, the record is not evidence about it and is not acted on."""
    repo, worktree = linked["repo"], linked["worktree"].resolve()
    manifest = lsro.manifest_path(worktree, repo)
    record = json.loads(manifest.read_text())
    record["clones"][0]["inode"] = record["clones"][0]["inode"] + 1
    manifest.write_text(json.dumps(record))
    with pytest.raises(rmwt.RefusedError) as e:
        rmwt.remove(repo, worktree)
    assert "no longer describes this path" in str(e.value)
    assert record["clones"][0]["path"] in str(e.value)
    assert worktree.is_dir()


def test_check_only_changes_nothing(linked: dict[str, Path]) -> None:
    repo, worktree = linked["repo"], linked["worktree"].resolve()
    before = {n: _state(repo / "data" / n) for n in lsro.STORES}
    r = rmwt.remove(repo, worktree, check_only=True)
    assert r["ok"] and r["checked_only"]
    assert worktree.is_dir()
    assert oct(stat.S_IMODE((worktree / "data" / "reference" / "sub" / "f.bin").stat().st_mode)) == "0o444"
    assert {n: _state(repo / "data" / n) for n in lsro.STORES} == before


# ---------------------------------------------------------------- plant (c): not git's worktree


def test_a_path_git_does_not_list_is_refused(planted: dict[str, Path], tmp_path: Path) -> None:
    stranger = tmp_path / "not-a-worktree"
    stranger.mkdir()
    (stranger / "data").mkdir()
    with pytest.raises(rmwt.RefusedError) as e:
        rmwt.remove(planted["repo"], stranger)
    assert "not in `git worktree list`" in str(e.value)
    assert str(stranger) in str(e.value)
    assert stranger.is_dir()


def test_the_main_worktree_is_refused(planted: dict[str, Path]) -> None:
    with pytest.raises(rmwt.RefusedError) as e:
        rmwt.remove(planted["repo"], planted["repo"])
    assert "MAIN worktree" in str(e.value)
    assert planted["repo"].is_dir()


def test_an_orphan_left_by_a_half_done_removal_is_refused_not_removed(linked: dict[str, Path]) -> None:
    """The state plant (d) measures into existence: a directory on disk that git has deregistered.
    It is the state in which a lane has no git route left, and the answer is a report, not an `rm`."""
    repo, worktree = linked["repo"], linked["worktree"].resolve()
    plain = _git(repo, "worktree", "remove", str(worktree), check=False)
    assert plain.returncode != 0
    assert worktree.is_dir()
    assert worktree not in rmwt.registered_worktrees(repo)  # git deregistered it as it failed
    with pytest.raises(rmwt.RefusedError) as e:
        rmwt.remove(repo, worktree, force=True)
    assert "orphan" in str(e.value)
    # left for its owner, and the stores it was cloned from are intact
    assert (repo / "data" / "reference" / "sub" / "f.bin").read_bytes() == b"store bytes\n"
    lsro.unlock(worktree / "data")  # so pytest's tmp_path cleanup can do its job


# ---------------------------------------------------------------- plant (d): the measurement


def test_d_plain_git_worktree_remove_fails_on_a_read_only_clone(linked: dict[str, Path]) -> None:
    """The tool's reason for existing, MEASURED rather than asserted. If this ever starts passing,
    `remove_worktree.py`'s premise is wrong and that is a finding to report, not a test to relax."""
    repo, worktree = linked["repo"], linked["worktree"].resolve()
    before = {n: _state(repo / "data" / n) for n in lsro.STORES}

    plain = _git(repo, "worktree", "remove", str(worktree), check=False)
    assert plain.returncode != 0, (
        "plain `git worktree remove` SUCCEEDED on a 0444 clone. The premise of "
        "scripts/remove_worktree.py is then wrong and the 2026-10-02 incident had another cause."
    )
    assert "Permission denied" in plain.stderr, plain.stderr
    assert worktree.is_dir(), "the worktree survived the failed removal"

    forced = _git(repo, "worktree", "remove", "--force", str(worktree), check=False)
    assert forced.returncode != 0, forced
    assert worktree.is_dir()

    # the second finding: the failed removal is not atomic. It deleted the tracked file and the .git
    # file and deregistered the worktree before it reached the read-only data, so --force cannot find a
    # working tree to remove and nothing git offers can clear what is left.
    assert not (worktree / "tracked.txt").exists()
    assert not (worktree / ".git").exists()
    assert "is not a working tree" in forced.stderr, forced.stderr

    # and through all of that the real stores are untouched
    assert {n: _state(repo / "data" / n) for n in lsro.STORES} == before
    lsro.unlock(worktree / "data")  # so pytest's tmp_path cleanup can do its job


def test_d_the_tool_succeeds_where_plain_removal_failed(planted: dict[str, Path]) -> None:
    """The same worktree, the same clones, the same git: the difference is unlocking the RECORDED
    paths first. Without that step this is the failure measured above."""
    repo, worktree = planted["repo"], planted["worktree"].resolve()
    assert lsro.link(repo, worktree)["ok"]
    before = {n: _state(repo / "data" / n) for n in lsro.STORES}
    r = rmwt.remove(repo, worktree)
    assert r["removed"] and r["git_exit"] == 0, r
    assert {n: _state(repo / "data" / n) for n in lsro.STORES} == before


def test_a_locked_clone_the_linker_could_not_record_is_a_FAILURE(
    planted: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """The record is written last, so a run that cannot write it has already locked the clones. That
    leaves a 0444 clone nothing records: the remover unlocks only recorded paths, so nobody can unlock
    it and nobody is told why -- which is the exact state that leaves a lane reaching for `rm -rf`.
    It must therefore be reported as a failure and not as a linked store."""
    repo, worktree = planted["repo"], planted["worktree"].resolve()
    monkeypatch.setattr(lsro, "write_manifest", lambda *a, **k: (_ for _ in ()).throw(OSError("no room")))
    out = lsro.link(repo, worktree)
    assert out["ok"] is False, out
    assert out["manifest"] is None
    assert any("could not be written" in f for f in out["failures"]), out["failures"]
    lsro.unlock(worktree / "data")  # so pytest's tmp_path cleanup can do its job
