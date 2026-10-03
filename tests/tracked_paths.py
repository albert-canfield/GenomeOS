# SPDX-License-Identifier: AGPL-3.0-or-later
"""The mirror of `tests/local_data.py`: an input git TRACKS, whose absence is a broken checkout.

`local_data.needs_local_data` covers the git-ignored stores: a test that needs one SKIPS BY NAME
where the store is absent and RUNS where it is present. This module covers the other half of the same
rule, which until now had no helper and was therefore written as a skip 31 times.

A path git tracks is in the commit, so it is present in this checkout, in CI, in a verdict worktree
and in the store-free pre-push leg. A `skipif` or a `pytest.skip()` waiting for such a file to be
absent can never fire. It is not a guard; it is a sentence that reads to every later reviewer as "this
test is conditional" while the condition is a constant. Worse, the identical shape keyed on a path
that is SOMETIMES absent hides a real failure as a skip, and the two are indistinguishable by eye --
which is why `scripts/inert_skip_guards.py` resolves the path and asks git instead of reading reasons.

So the rule, as the project adopted it on 2026-10-03: a skip keyed on a TRACKED path must RAISE.

`must_be_present` RAISES ON ITS OWN MISUSE TOO, for the reason `local_data.check_is_machine_local`
does. If the path it names is not tracked -- a typo, or a file that has since been git-ignored or
deleted from the commit -- then calling it here is a claim about the tree that is no longer true, and
the fix is to name a real input or to use `needs_local_data`. Converting a skip into a raise on a path
that can genuinely be absent would turn a quiet skip into a red push, and that is worse than the skip
it replaced; this is the check that stops the conversion being done on the wrong path.
"""

from __future__ import annotations

import subprocess
from functools import cache
from pathlib import Path

#: The repository root, as every test is run from it and as `tests/local_data.py` spells it.
ROOT = Path(__file__).resolve().parents[1]


class NotTrackedError(Exception):
    """`must_be_present` was called on a path git does not track, so the claim it makes is wrong."""


class BrokenCheckoutError(AssertionError):
    """A tracked input is missing. Not a machine without a store: a checkout that is not the commit."""


@cache
def tracked_files() -> frozenset[str]:
    """Every path git tracks, as git spells it. Asked once per session, and asked of git.

    Asked of `git ls-files` rather than matched against store prefixes because `data/results` is
    git-ignored by PATTERN with `!` re-includes and hundreds of its files are committed, so a prefix
    test would get it wrong in both directions.
    """
    out = subprocess.run(["git", "ls-files"], capture_output=True, text=True, check=True, cwd=ROOT)
    return frozenset(out.stdout.splitlines())


def as_relative(path: str | Path) -> str:
    """`path` as git spells it: repository-relative, forward slashes. Absolute paths are relativised."""
    p = Path(path)
    if p.is_absolute():
        try:
            p = p.resolve().relative_to(ROOT.resolve())
        except ValueError:
            return str(p)
    return p.as_posix()


def is_tracked(path: str | Path) -> bool:
    """Whether git tracks `path`, counting a directory as tracked when it holds a tracked file."""
    rel = as_relative(path)
    if rel in tracked_files():
        return True
    prefix = rel.rstrip("/") + "/"
    return any(t.startswith(prefix) for t in tracked_files())


def must_be_present(*paths: str | Path, was: str | None = None) -> None:
    """Raise unless every path is there. `was` carries the words of the skip this replaced.

    Two raises, and they are different findings. `NotTrackedError` says this call is wrong about the
    tree. `BrokenCheckoutError` says the tree is wrong: a committed file is missing from a checkout,
    which no amount of fetching data fixes and which must never read as "not run here".
    """
    for path in paths:
        rel = as_relative(path)
        if not is_tracked(path):
            raise NotTrackedError(
                f"must_be_present({rel!r}): git does not track it, so it CAN be absent and raising "
                f"here would turn a skip into a red push. Name a tracked input, or mark the test "
                f"with needs_local_data naming the git-ignored store it actually reads"
            )
        if not Path(path).exists():
            said = f' The skip this replaced said: "{was}".' if was else ""
            raise BrokenCheckoutError(
                f"{rel} is TRACKED by git, so it is present in every checkout of this commit -- "
                f"here, in CI, in a verdict worktree and in the store-free pre-push leg. Its "
                f"absence is a broken checkout, not a machine without a store, and it is not a "
                f"reason to skip.{said}"
            )
