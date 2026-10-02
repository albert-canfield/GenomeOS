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

A STORE REACHED THROUGH A SYMLINK SKIPS BY NAME. `git check-ignore` stops at a symbolic link: it exits
128 with "fatal: pathspec ... is beyond a symbolic link" when the store is linked and 0 when it is a
real directory, so a worktree with linked stores made this module raise and 20 tests ERROR -- a RED
status file for a tree that was nothing of the kind, twice. The name is therefore truncated to the
symlink and git is asked about that, because an ignore rule is about a name and a directory's rule
covers its subtree. `ignored_by_name` records what was measured and the two fixes that do not work.

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


#: git's own words when a pathspec passes through a symbolic link, recorded for the reader. The retry
#: below is NOT decided by matching this string -- it is decided by asking the filesystem which
#: component is a symlink -- because a message is a version detail and a symlink is a fact.
BEYOND_A_SYMLINK = "beyond a symbolic link"


def _ask_git(path: str) -> tuple[bool | None, str]:
    """(whether git ignores `path`, git's stderr). None when git refused to answer the question at all.

    0 is ignored, 1 is not ignored, and anything else -- 128, in practice -- is a refusal and never a
    licence to skip.
    """
    r = subprocess.run(
        ["git", "check-ignore", "-q", "--no-index", path],
        capture_output=True,
        cwd=ROOT or None,
    )
    if r.returncode in (0, 1):
        return r.returncode == 0, ""
    return None, r.stderr.decode(errors="replace")[-300:].strip()


def symlink_boundary(path: str) -> str | None:
    """The shortest prefix of `path` that is itself a symlink, or None when no component is one.

    This is the deepest name `git check-ignore` can be asked about, because git stops at a symlink:
    `data/cache/entex/x.csv` with `data/cache` a symlink is refused with exit 128 and
    "fatal: pathspec ... is beyond a symbolic link", while `data/cache` itself answers 0. Measured
    2026-10-02; the same probe gets 0 against a real directory, which is why nothing noticed until a
    worktree had the stores linked rather than cloned.
    """
    parts = Path(path).parts
    for i in range(1, len(parts) + 1):
        candidate = Path(*parts[:i])
        if candidate.is_symlink():
            return str(candidate)
    return None


@cache
def ignored_by_name(path: str) -> tuple[bool, str]:
    """(whether git ignores it, the NAME git actually answered about). Cached: asked once per session.

    Asked of git rather than matched against a list of store prefixes because `data/results` is ignored
    by pattern with re-includes and 848 of its files are committed, so a prefix test would get it wrong
    in both directions.

    A STORE REACHED THROUGH A SYMLINK IS ANSWERED BY NAME, which is the whole of the 2026-10-02 fix.
    Twenty tests ERRORED and two RED status files were written for trees that were nothing of the kind,
    because `.git/hooks/pre-push` was a stale copy from 04:34 that still did
    `ln -s "$root/data/$d" "$tmp/data/$d"` while `scripts/pre-push.sh` had moved to read-only clones:
    the push's worktree had `data/cache` as a symlink, every path under it was "beyond a symbolic link"
    to git, and `is_git_ignored` raised. So the name is truncated to the symlink and git is asked about
    THAT: a gitignore entry that ignores a directory ignores its whole subtree, so "data/cache is
    ignored" answers "data/cache/entex/... is machine-local" -- by name, which is what the marker is for.

    TWO THINGS MEASURED AND REJECTED, written down so they are not tried again.

    `realpath` BEFORE `check-ignore` DOES NOT WORK, and it was the first fix proposed. The ignore
    decision is about a NAME, so resolving the path asks a different question. Measured 2026-10-02 in
    the exact failing configuration -- a worktree whose `data/cache` is a symlink into the main checkout
    -- the resolved path is in a DIFFERENT working tree of the same repository and git exits 128 again:
    "fatal: /.../R/data/cache/entex/x.csv is outside repository at '/.../W'". Where the target happens to
    sit inside the same worktree under another name it is worse than an error: git exits 1, "not
    ignored", and the marker is reported as wrong instead of unanswerable. `--no-index` makes no
    difference to either reading; it was checked rather than assumed.

    `$TMPDIR` BEING BEHIND `/var -> private/var` IS NOT THE CAUSE, which was the second diagnosis. A
    worktree created under `$TMPDIR` with the stores as REAL directories answers 0, measured twice: with
    the main repository also under `/var`, and with it on an unsymlinked path. git resolves a linked
    worktree's root physically (`rev-parse --show-toplevel` prints `/private/var/...`) and a relative
    pathspec is taken from the process's physical cwd, so the leading symlink never enters the question.
    """
    answer, why = _ask_git(path)
    if answer is not None:
        return answer, path
    boundary = symlink_boundary(path)
    if boundary is None:
        raise NotMachineLocalError(
            f"git could not say whether {path} is ignored ({why}), and no component of it is a "
            f"symbolic link, so whether this test needs machine-local data is unknown and it is "
            f"not skipped"
        )
    if boundary != path:
        at_boundary, boundary_why = _ask_git(boundary)
        if at_boundary:
            return True, boundary
        if at_boundary is None:
            raise NotMachineLocalError(
                f"git could not say whether {path} is ignored ({why}), and asking it about the "
                f"symbolic link {boundary} that it stops at was refused too ({boundary_why}), so "
                f"whether this test needs machine-local data is unknown and it is not skipped"
            )
    raise NotMachineLocalError(
        f"needs_local_data({path!r}): git stops at the symbolic link {boundary} and does not ignore "
        f"that name, so whether the path beyond it is machine-local cannot be asked and has not been "
        f"guessed. Name the git-ignored store itself, or drop the marker"
    )


def is_git_ignored(path: str) -> bool:
    """Whether git ignores `path`, asked of git. The answer alone; `ignored_by_name` carries the name."""
    return ignored_by_name(path)[0]


def check_is_machine_local(pattern: str) -> None:
    """Raise unless `pattern` names a git-ignored path. The marker's own correctness, checked."""
    path = literal_prefix(pattern)
    ignored, answered_about = ignored_by_name(path)
    if not ignored:
        raise NotMachineLocalError(
            f"needs_local_data({pattern!r}): git does not ignore {answered_about}, so it is in the "
            f"commit and present in any worktree. A marker here would skip for no stated reason; name "
            f"the git-ignored store the test actually needs, or drop the marker"
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
