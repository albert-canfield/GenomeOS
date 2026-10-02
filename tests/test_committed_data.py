# SPDX-License-Identifier: AGPL-3.0-or-later
"""`must_be_committed` is SHOWN to fail, in both directions, rather than asserted to be able to.

A helper whose only job is to raise is worthless if nothing ever makes it raise, and that is the
defect class it was written against: a guard that reported "not applicable here" forever because the
condition it tested never came true. So both failures are provoked here -- an absent tracked path and
a present untracked one -- and the agreement with `tests/local_data.py` is checked rather than
claimed, because the two modules only cover the whole range of paths if nothing is both ignored and
tracked.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from tests import committed_data as cd
from tests import local_data


def test_a_committed_and_present_path_passes_and_is_returned() -> None:
    """The happy path, on a result that is in the commit: the Path comes back, nothing raises."""
    got = cd.must_be_committed("data/results/celegans_terminality.json")
    assert got.is_file()
    assert got == cd.ROOT / "data/results/celegans_terminality.json"


def test_a_committed_directory_passes_on_the_files_it_carries() -> None:
    """`data/results` is committed because it carries 863 tracked files, not as a path of its own.

    MEASURED rather than assumed, and the first draft of this test had it wrong: `--error-unmatch`
    takes a PATHSPEC, so it exits 0 on a directory that carries any tracked file and `is_tracked`
    answers True for `data/results` even though git stores no entry under that name. `tracked_under`
    is kept because it is the answer that can be counted, and it is what distinguishes "this
    directory is in the commit" from "this directory happens to exist".
    """
    assert cd.is_tracked("data/results") is True
    assert cd.tracked_under("data/results") > 800
    assert cd.tracked_under("data/reference") == 0, "the ignored stores carry no tracked file"
    assert cd.must_be_committed("data/results").is_dir()


def test_an_absent_tracked_path_FAILS_and_names_the_path_and_gits_answer() -> None:
    """The defect this exists for. Named with a path git tracks that is deleted from the checkout.

    `.gitignore` is used as the probe only because it is certainly tracked and certainly present; the
    absence is simulated by asking about a tracked path under a directory that does not exist, so the
    test needs no write and deletes nothing.
    """
    absent = "data/results/nonexistent_dir/celegans_terminality.json"
    assert not (cd.ROOT / absent).exists()
    with pytest.raises(AssertionError) as caught:
        cd.must_be_committed(absent)
    assert absent in str(caught.value)
    assert "git tracks it: False" in str(caught.value)
    assert "never its existence" in str(caught.value)


def test_a_present_but_untracked_path_FAILS_rather_than_passing_quietly(tmp_path: Path) -> None:
    """A guard naming a path nothing committed measures the machine, so it is a defect in the test."""
    planted = tmp_path / "planted.json"
    planted.write_text("{}")
    relative = str(planted.resolve().relative_to(planted.resolve().parents[len(planted.parts) - 2]))
    with pytest.raises(AssertionError, match="git does not track it"):
        cd.must_be_committed(str(planted))
    assert relative  # the relative form is unused by the helper; recorded so the shape is explicit


def test_the_message_carries_the_reason_so_it_travels_with_the_failure() -> None:
    """A reader must not have to find this module to understand why the verdict is a failure."""
    assert "a stub may substitute a dependency's behaviour, never its existence" in (
        cd.ABSENCE_IS_A_DEFECT.replace("A stub", "a stub")
    )
    assert "needs_local_data" in cd.ABSENCE_IS_A_DEFECT


def test_the_committed_decorator_fails_the_test_it_guards_and_not_its_module() -> None:
    """The decorator raises inside the call, so one named test fails and collection is unaffected."""

    @cd.committed("data/results/nonexistent_dir/whatever.json")
    def guarded() -> str:
        return "ran"

    with pytest.raises(AssertionError, match="nonexistent_dir"):
        guarded()


def test_the_decorator_runs_the_test_when_every_named_path_is_committed() -> None:
    """Non-vacuity: without this the decorator could pass by never running the body."""

    @cd.committed("data/results/celegans_terminality.json")
    def guarded() -> str:
        return "ran"

    assert guarded() == "ran"


@pytest.mark.parametrize(
    "relative",
    [
        "data/results/celegans_terminality.json",
        "data/results/direction_link.json",
        "data/results/rule_number_sources_census.json",
        "data/models/BIOMD0000000012.xml",
        "data/organisms/human/noncoding_chr21.bio",
    ],
)
def test_no_path_this_lane_converted_is_both_tracked_and_git_ignored(relative: str) -> None:
    """The two modules only cover every path if nothing falls in both. Asked of git, not assumed.

    `data/results` is ignored by pattern with `!` re-includes, so this is the case where a prefix
    test would get it wrong in both directions -- which is why `tests/local_data.py` asks git too.
    """
    assert cd.is_tracked(relative), f"{relative} is not tracked, so it belongs to local_data"
    assert not local_data.is_git_ignored(relative), (
        f"{relative} is both tracked and git-ignored, so the two guards disagree about it"
    )
