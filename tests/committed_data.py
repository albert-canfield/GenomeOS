# SPDX-License-Identifier: AGPL-3.0-or-later
"""A committed artefact's absence is a DEFECT, not a machine difference, so it fails and never skips.

The mirror of `tests/local_data.py`, and the same question asked the other way round. That module
exists for the git-IGNORED stores -- `data/cache`, `data/reference`, `data/knowledge` -- which are
machine-local by definition, so a test needing one SKIPS BY NAME where it is absent; and it RAISES
on a path git does not ignore, because a marker naming a path that is in the commit would otherwise
skip for no stated reason. This module is the other half: a path git DOES track is present in every
checkout and every worktree of the commit, so its absence is something broken rather than something
missing, and a guard on it must fail.

THE RULE IT ENFORCES is the project's own, from docs/LESSONS.md: **a stub may substitute a
dependency's behaviour, never its existence.** A `pytest.skip` on a committed result substitutes its
existence -- a deleted, truncated or corrupted published artefact reads as "not applicable here"
instead of as the defect it is, and a reader cannot tell "not run" from "passed".

WHY IT MATTERS HERE RATHER THAN IN GENERAL. 863 files under `data/results` are committed, while
`data/results` is itself ignored by PATTERN with `!` re-includes, so whether a path under it is
machine-local cannot be decided from its prefix and is asked of git -- exactly as
`tests/local_data.ignored_by_name` asks. The two modules therefore agree by construction: a path is
either ignored (skip by name, with the store named) or tracked (fail, with the path named), and
nothing is both.

MEASURED, 2026-10-02 (lane-noskip): 55 skip guards across 30 test files stood on paths git tracks,
with 8 further guards mixing a tracked path and an ignored store in one condition. The three in
`tests/test_respmap_v2.py` were converted first, at `4f44dbf`, where `must_be_committed` was written;
it is moved here rather than copied so there is one helper and not two. All 63 stood on paths that
are PRESENT in this checkout, so the conversion changed no verdict -- which is the point: the guards
were never reporting a machine difference, they were covering for one that had not happened yet.
"""

from __future__ import annotations

import subprocess
from collections.abc import Callable
from functools import cache, wraps
from pathlib import Path
from typing import Any, TypeVar

#: The repository root, derived from this file so it does not depend on the working directory.
ROOT = Path(__file__).resolve().parents[1]

#: Why a guard here fails instead of skipping. Carried into every message so the reason travels with
#: the failure and a reader never has to find this module to understand the verdict.
ABSENCE_IS_A_DEFECT = (
    "a committed result's absence FAILS and does not skip. A stub may substitute a dependency's "
    "behaviour, never its existence, and 863 files under data/results are in the commit, so a "
    "deleted or corrupted published artefact is a defect rather than a machine difference. "
    "tests/local_data.py's needs_local_data marker is for the git-IGNORED stores and raises on a "
    "path git does not ignore, which is this same distinction asked the other way round"
)

F = TypeVar("F", bound=Callable[..., Any])


@cache
def is_tracked(relative: str) -> bool:
    """Whether git tracks `relative`, asked of git. Cached: asked once per path per session.

    `--error-unmatch` is what makes this an answer about the INDEX rather than about the filesystem,
    which is the whole point: the question is whether the commit carries this path, and that stays
    answerable when the file on disk has been deleted.
    """
    return (
        subprocess.run(
            ["git", "ls-files", "--error-unmatch", "--", relative],
            cwd=ROOT,
            capture_output=True,
        ).returncode
        == 0
    )


@cache
def tracked_under(relative: str) -> int:
    """How many paths git tracks under `relative`. A directory is committed if it carries files."""
    out = subprocess.run(["git", "ls-files", "--", relative], cwd=ROOT, capture_output=True, text=True)
    return len([line for line in out.stdout.splitlines() if line])


def must_be_committed(relative: str) -> Path:
    """The path, or an AssertionError. A committed artefact's absence is a defect, not a difference.

    `AssertionError` and not a custom exception so pytest reports FAILED rather than ERROR: this is a
    claim about the tree that has come out false, which is what a failure is.

    Git's answer is REPORTED rather than assumed, in both directions. An absent path says whether git
    tracks it, so "the result was deleted" is distinguishable from "this test names a path that was
    never committed" -- and the second is a defect in the test, which is why a present-but-untracked
    path fails too: nothing in the commit fixes it, so the guard would have been measuring the
    machine either way.
    """
    path = ROOT / relative
    tracked = is_tracked(relative) or tracked_under(relative) > 0
    if not path.exists():
        raise AssertionError(
            f"{relative} is absent from this checkout and git tracks it: {tracked}. {ABSENCE_IS_A_DEFECT}"
        )
    if not tracked:
        raise AssertionError(
            f"{relative} is present but git does not track it, so nothing in the commit fixes it "
            f"and a guard on it would measure the machine. {ABSENCE_IS_A_DEFECT}"
        )
    return path


def committed(*relatives: str) -> Callable[[F], F]:
    """Decorator: every named path must be committed and present, or the test FAILS naming it.

    The one-line replacement for `@pytest.mark.skipif(not RESULT.exists(), reason=...)` on a tracked
    path. It is a decorator and not a pytest marker on purpose: a marker is read at collection, and a
    marker whose job is to fail is a marker that fails the collection of its whole module, which
    tells a reader less than a named failure in the one test that depends on the artefact.
    """

    def decorate(test: F) -> F:
        @wraps(test)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            for relative in relatives:
                must_be_committed(relative)
            return test(*args, **kwargs)

        return wrapper  # type: ignore[return-value]

    return decorate
