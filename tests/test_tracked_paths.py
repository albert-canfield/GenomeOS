# SPDX-License-Identifier: AGPL-3.0-or-later
"""`tests/tracked_paths.py`: the raise that replaces a skip which could never fire.

Both raises are exercised on real paths of this repository, and the second one -- the refusal to be
used on a path that can genuinely be absent -- is the one that matters most. Converting a skip into a
raise on an untracked path turns a quiet skip into a red push, which is worse than the skip it
replaced, so the helper refuses rather than trusting the lane that calls it.
"""

from __future__ import annotations

import pytest
import tracked_paths as tp

#: committed, so present in every checkout of this commit
A_TRACKED_RESULT = "data/results/clause1.json"
#: in a git-ignored store, so genuinely absent on a machine that has not fetched it
A_LOCAL_STORE_FILE = "data/cache/entex/alphagenome_track_metadata_copy.csv"


def test_a_tracked_path_that_is_present_passes_and_says_nothing() -> None:
    tp.must_be_present(A_TRACKED_RESULT, "data/results", tp.ROOT / A_TRACKED_RESULT)


def test_a_tracked_path_that_is_absent_RAISES_and_names_it_as_a_broken_checkout(monkeypatch) -> None:
    """The conversion's whole point: absence of a committed file is red, by name, never a skip.

    A file that is tracked AND absent cannot be made by deleting one -- that would damage the shared
    checkout -- so git's answer is the thing substituted: the tracked set is told it carries a path
    that is not on disk, which is exactly the state a broken checkout is in.
    """
    gone = "data/results/a_result_this_commit_does_not_carry.json"
    with pytest.raises(tp.NotTrackedError):  # untracked, so refused before existence is asked
        tp.must_be_present(gone)
    monkeypatch.setattr(tp, "tracked_files", lambda: frozenset({gone}))
    with pytest.raises(tp.BrokenCheckoutError, match="TRACKED by git") as raised:
        tp.must_be_present(gone, was="the result is not on this machine")
    assert gone in str(raised.value)
    assert "the result is not on this machine" in str(raised.value), (
        "the words of the skip it replaced are carried into the failure, so nothing said is lost"
    )


def test_the_helper_REFUSES_a_path_git_does_not_track() -> None:
    """A conversion aimed at a machine-local store would be a red push on every clean machine."""
    with pytest.raises(tp.NotTrackedError, match="needs_local_data"):
        tp.must_be_present(A_LOCAL_STORE_FILE)


def test_the_two_paths_this_test_stands_on_are_what_it_says_they_are() -> None:
    """Non-vacuity: if `clause1.json` left the commit, or `data/cache` entered it, both tests above
    would still pass while meaning the opposite."""
    assert tp.is_tracked(A_TRACKED_RESULT)
    assert not tp.is_tracked(A_LOCAL_STORE_FILE)
    assert tp.as_relative(tp.ROOT / "data" / "results") == "data/results"
    assert tp.is_tracked("data/results"), "a directory holding tracked files counts as tracked"
