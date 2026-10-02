# SPDX-License-Identifier: AGPL-3.0-or-later
"""`code_cleanliness` against a planted git repository, so the verdict does not depend on this tree.

Until 2026-10-02 the only test that two of the fields worked at all was a pair of assertions in
`tests/test_crispri_benchmark_v2.py` and `tests/test_crispri_published_v2.py` reading
`own_code_is_committed is True` and `foreign_uncommitted_code_on_the_counting_path == []` against the
LIVE checkout. Those assert the state of a shared working tree, not the behaviour of any code: four
lanes work in this checkout, so at least one usually has uncommitted code on some counting path, and
the assertions went red for a peer's edit and green in a quiet moment while the mechanism behind them
was never exercised either way. Two consecutive full runs had that assertion as the only red in about
4,160 tests, naming a different lane's file each time.

What they claimed to check is checked here instead, hermetically: a temporary repository is built with
an entry script that imports a module, and the three states are made by what is committed in it, not
by what any lane happens to be editing. The invariants that do hold in every tree state stayed in the
two suites, in the style of `tests/test_context_evidence.py`.

Nothing here touches enforcement, and the refusal that matters never lived in a unit test.
`scripts/check_staged.py` refuses the commit of any staged `data/results/*.json` whose block has
`own_code_is_committed` false or `foreign_uncommitted_code_on_the_counting_path` non-empty, which is
where a result gets published; `genomeos/results.py` quarantines a block that did not come from the
shared function at all, which `mf.cleanliness_problems` decides by KEY SET and not by either of these
two values. So this file moves a readiness check out of the unit suite, where it was testing the tree,
and takes nothing away from either refusal.
"""

from __future__ import annotations

import inspect
import os
import subprocess
import textwrap
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from genomeos import manifest as mf

#: The planted entry script and the module it imports. Both sit under `scripts/`, because `mf.is_code`
#: counts a path as code only under CODE_ROOTS: a planted `entry.py` at a repository root would be
#: ignored by the revision stamp and all three cases below would pass for the wrong reason.
ENTRY = "scripts/planted_entry.py"
MODULE = "scripts/planted_helper.py"

#: The entry is the writing lane's own code; the module it imports belongs to another lane.
OWN: tuple[str, ...] = (ENTRY,)

ENTRY_SOURCE = """\
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import planted_helper

print(planted_helper.VALUE)
"""


#: Variables that would point git at another repository's index, directory or work tree.
#: `mf.code_revision` shells out with the inherited environment, and `scripts/commit_own.sh` runs the
#: suite with GIT_INDEX_FILE set to a private index, so without this the planted repository's status
#: could be read against a different repository and all three cases below would be meaningless.
GIT_ENV_TO_DROP = ("GIT_INDEX_FILE", "GIT_DIR", "GIT_WORK_TREE", "GIT_OBJECT_DIRECTORY")


@pytest.fixture(autouse=True)
def _no_inherited_git_environment(monkeypatch: pytest.MonkeyPatch) -> None:
    for name in GIT_ENV_TO_DROP:
        monkeypatch.delenv(name, raising=False)


def _git(repo: Path, *args: str) -> None:
    env = dict(os.environ)
    subprocess.run(["git", "-C", str(repo), *args], check=True, capture_output=True, text=True, env=env)


def _planted_repo(tmp_path: Path, committed: tuple[str, ...]) -> Path:
    """A repository holding ENTRY and MODULE, with exactly `committed` committed and the rest untracked."""
    repo = tmp_path / "planted"
    (repo / "scripts").mkdir(parents=True)
    (repo / MODULE).write_text("VALUE = 1\n")
    (repo / ENTRY).write_text(ENTRY_SOURCE)
    _git(repo, "init", "-b", "main")
    _git(repo, "config", "user.email", "lane@example.invalid")
    _git(repo, "config", "user.name", "planted")
    _git(repo, "config", "commit.gpgsign", "false")
    _git(repo, "add", "--", *committed)
    _git(repo, "commit", "-m", "planted")
    return repo


def _the_planted_module_is_named(block: dict[str, Any]) -> None:
    """Case (a)'s claim, as one callable so the counterfactual can hold the same claim against a copy."""
    assert MODULE in block["foreign_uncommitted_code_on_the_counting_path"]
    for p in block["foreign_uncommitted_code_on_the_counting_path"]:
        assert p in block["counting_path"] and p in block["foreign_uncommitted_code"]


def test_the_counting_path_of_the_planted_repository_holds_both_files(tmp_path: Path) -> None:
    """The closure is read from the imports, so it is the same in all three cases below."""
    repo = _planted_repo(tmp_path, (ENTRY, MODULE))
    assert mf.counting_path(ENTRY, repo) == [ENTRY, MODULE]


# --- (a) a foreign file left uncommitted ----------------------------------------------------------


def test_an_uncommitted_foreign_module_on_the_path_is_named(tmp_path: Path) -> None:
    repo = _planted_repo(tmp_path, (ENTRY,))
    block = mf.code_cleanliness(ENTRY, OWN, repo)
    _the_planted_module_is_named(block)
    assert block["own_uncommitted_code"] == []
    assert block["own_code_is_committed"] is True
    assert block["dirty"] is True


def test_an_uncommitted_foreign_module_off_the_path_is_not_named(tmp_path: Path) -> None:
    """The `on_the_counting_path` field is the filter it says it is, not a copy of the dirty list."""
    repo = _planted_repo(tmp_path, (ENTRY, MODULE))
    (repo / "scripts" / "planted_elsewhere.py").write_text("VALUE = 2\n")
    block = mf.code_cleanliness(ENTRY, OWN, repo)
    assert block["foreign_uncommitted_code"] == ["scripts/planted_elsewhere.py"]
    assert block["foreign_uncommitted_code_on_the_counting_path"] == []


# --- (b) everything committed ---------------------------------------------------------------------


def test_with_everything_committed_the_block_is_clean(tmp_path: Path) -> None:
    repo = _planted_repo(tmp_path, (ENTRY, MODULE))
    block = mf.code_cleanliness(ENTRY, OWN, repo)
    assert block["foreign_uncommitted_code_on_the_counting_path"] == []
    assert block["own_code_is_committed"] is True
    assert block["own_uncommitted_code"] == []
    assert block["foreign_uncommitted_code"] == []
    assert block["dirty"] is False


# --- (c) the entry script itself left uncommitted -------------------------------------------------


def test_an_uncommitted_entry_script_is_the_writing_lane_s_own(tmp_path: Path) -> None:
    repo = _planted_repo(tmp_path, (MODULE,))
    block = mf.code_cleanliness(ENTRY, OWN, repo)
    assert block["own_code_is_committed"] is False
    assert block["own_uncommitted_code"] == [ENTRY]
    assert block["foreign_uncommitted_code"] == []
    assert block["foreign_uncommitted_code_on_the_counting_path"] == []


def test_a_modified_committed_entry_script_counts_the_same_as_an_untracked_one(tmp_path: Path) -> None:
    """Case (c) by the other route: committed and then edited, which is how a lane actually gets there."""
    repo = _planted_repo(tmp_path, (ENTRY, MODULE))
    (repo / ENTRY).write_text(ENTRY_SOURCE + "print('edited')\n")
    block = mf.code_cleanliness(ENTRY, OWN, repo)
    assert block["own_code_is_committed"] is False
    assert block["own_uncommitted_code"] == [ENTRY]


# --- the counterfactual ---------------------------------------------------------------------------

#: The one line of `mf.code_cleanliness` that computes the field. The counterfactual below strips it
#: from a copy of the function's own source rather than from a paraphrase, so the copy cannot drift
#: from the implementation: if this text stops matching, the test fails instead of passing vacuously.
PATH_FILTER = '"foreign_uncommitted_code_on_the_counting_path": [p for p in foreign if p in path],'


def _cleanliness_without_the_path_filter() -> Callable[..., dict[str, Any]]:
    source = inspect.getsource(mf.code_cleanliness)
    assert source.count(PATH_FILTER) == 1, (
        "the line this counterfactual strips is no longer in genomeos/manifest.py as written here, so "
        "the counterfactual would pass without testing anything: update PATH_FILTER"
    )
    stripped = source.replace(PATH_FILTER, '"foreign_uncommitted_code_on_the_counting_path": [],')
    namespace = dict(vars(mf))
    exec(textwrap.dedent(stripped), namespace)  # noqa: S102 - a copy of this repository's own source
    return namespace["code_cleanliness"]


def test_case_a_fails_against_a_copy_with_the_foreign_on_path_computation_stripped(tmp_path: Path) -> None:
    """Without this, the invariants in the two suites are unverified: they would hold of a block that
    never names anything. The copy keeps every other field and loses only the filter."""
    repo = _planted_repo(tmp_path, (ENTRY,))
    stripped = _cleanliness_without_the_path_filter()
    block = stripped(ENTRY, OWN, repo)
    assert block["foreign_uncommitted_code"] == [MODULE]  # the copy still sees the uncommitted file
    assert block["counting_path"] == [ENTRY, MODULE]  # and still computes the closure
    with pytest.raises(AssertionError):
        _the_planted_module_is_named(block)
