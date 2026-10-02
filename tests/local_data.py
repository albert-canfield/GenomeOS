# SPDX-License-Identifier: AGPL-3.0-or-later
"""Which machine-local data stores a test needs, and whether this checkout has them.

The half of `needs_local_data` that is not pytest wiring (that is `tests/conftest.py`). It exists as a
module beside the conftest for the same reason `tests/extras_lock.py` does: the decision is testable on
its own, and a decision only a conftest hook can make is a decision nothing checks.

The rule this implements, as the supervisor amended it on 2026-10-02. The project's acceptance rule was
"a status-file verdict from a worktree of the committed tree", and that rule was unachievable as stated:
`data/cache`, `data/reference` and `data/knowledge` are git-ignored machine-local data, absent from a
worktree by definition, and 21 tests therefore FAILED in one -- 16 of them in `tests/test_context_evidence.py`
alone, every one on `data/cache/entex/alphagenome_track_metadata_copy.csv`. The module was behaving
correctly: it refuses to guess a cell-to-biosample mapping from names. A verdict determined by the
environment rather than by the code is the defect class the project spent 2026-10-02 removing, and it had
reappeared inside the rule adopted to escape it.

The amended rule has two halves and this is the second. The verdict is a worktree of the committed tree
WITH the git-ignored stores linked read-only, so a test that needs local data RUNS there; and where the
stores are genuinely absent -- CI -- such a test SKIPS BY NAME and never fails. By name means the skip
reason identifies the missing store, so a reader can tell "not run here" from "passed".

TWO THINGS IT MUST NOT DO, both of them lessons rather than preferences.

A marker never weakens an assertion. It changes WHERE a test runs, never what it claims. A test whose
only content is that local data exists is a readiness check filed as a unit test, and
docs/LESSONS.md ("The verdict moved twice without the mechanism moving once") is about exactly that
class: such a test is left alone and reported, not marked.

A path that is NOT git-ignored RAISES rather than skipping. `requires_extra` set that rule for an extra
pyproject does not provide, and the reason carries over exactly: a marker naming a path that is in the
commit would otherwise skip, or run, for no stated reason, and a typo would read as a clean skip
forever. `data/results` is the case that makes this necessary -- it is ignored by PATTERN with `!`
re-includes, 848 of its files are committed, so whether a path under it is machine-local cannot be
decided by its prefix and is asked of git.
"""

from __future__ import annotations

import glob
import subprocess
from functools import cache
from pathlib import Path

#: The repository root, as every test is run from it.
ROOT = Path()


class NotMachineLocalError(Exception):
    """A `needs_local_data` path that git does not ignore: the marker is wrong, so it raises."""


def literal_prefix(pattern: str) -> str:
    """The part of a path before its first wildcard, which is what git can be asked about.

    `data/knowledge/human_panel/chr*/storage_catalogue.json.gz` is a real marker argument (24 per-
    chromosome catalogues, none committed), and `git check-ignore` takes a path and not a glob.
    """
    out: list[str] = []
    for part in Path(pattern).parts:
        if any(c in part for c in "*?["):
            break
        out.append(part)
    return str(Path(*out)) if out else pattern


@cache
def is_git_ignored(path: str) -> bool:
    """Whether git ignores `path`, asked of git. Cached: a marked test asks once per session.

    Asked of git rather than matched against a list of store prefixes because `data/results` is ignored
    by pattern with re-includes and 848 of its files are committed, so a prefix test would get it wrong
    in both directions.
    """
    r = subprocess.run(
        ["git", "check-ignore", "-q", "--no-index", path],
        capture_output=True,
        cwd=ROOT or None,
    )
    # 0: ignored. 1: not ignored. anything else: git could not answer, which is not a licence to skip.
    if r.returncode not in (0, 1):
        raise NotMachineLocalError(
            f"git could not say whether {path} is ignored ({r.stderr.decode()[-200:].strip()}), so "
            f"whether this test needs machine-local data is unknown and it is not skipped"
        )
    return r.returncode == 0


def check_is_machine_local(pattern: str) -> None:
    """Raise unless `pattern` names a git-ignored path. The marker's own correctness, checked."""
    if not is_git_ignored(literal_prefix(pattern)):
        raise NotMachineLocalError(
            f"needs_local_data({pattern!r}): git does not ignore this path, so it is in the commit and "
            f"present in any worktree. A marker here would skip for no stated reason; name the "
            f"git-ignored store the test actually needs, or drop the marker"
        )


def present(pattern: str) -> bool:
    """Whether this checkout holds the store the pattern names. A glob needs at least one match."""
    if any(c in pattern for c in "*?["):
        return bool(glob.glob(pattern))
    return Path(pattern).exists()


def missing(patterns: tuple[str, ...]) -> list[str]:
    """Those of `patterns` this checkout does not hold, each checked for being machine-local first."""
    out = []
    for pattern in patterns:
        check_is_machine_local(pattern)
        if not present(pattern):
            out.append(pattern)
    return out


def skip_reason(patterns: list[str], how: str | None = None) -> str:
    """The sentence a skipped test prints. It NAMES the store, so "not run here" cannot read as "passed".

    `how` is the recovery instruction, in the convention the rest of the suite's skip reasons use
    ("fetch chr21 first", "run scripts/celegans_fates.py"): where it is known, it is said.
    """
    tail = f"; {how}" if how else ""
    return (
        f"NOT RUN HERE, not passed: this test needs machine-local data that is git-ignored and absent "
        f"from this checkout: {', '.join(patterns)}. It runs in a checkout or a verdict worktree that "
        f"has the stores{tail}"
    )
