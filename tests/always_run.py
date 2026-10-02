# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""ALWAYS_RUN: the cross-cutting invariants every lane acceptance runs, on top of its own tests.

WHY IT EXISTS, and it is a procedure change rather than a new check. The amended light-lane rule says
a lane's acceptance is its TARGETED tests plus lint, in a worktree of its committed tree, and the
full suite is the push's job. That is right for cost and it has one hole: an invariant that no lane
owns is in no lane's targeted set, so it is skipped BY CONSTRUCTION and not by anyone's choice. On
2026-10-02 `genomeos/lang/rule_cell_provenance.py` was committed importing
`genomeos/attribution/measured.py`, which put the Apache-2.0 engine in the position of importing the
AGPL-3.0 application against LICENSING.md decision D40, and reddened every sha after it. The lane's
own 60 tests passed. `tests/test_engine_boundary.py` would have caught it in under two seconds and
was not in the lane's set.

So from 2026-10-02 every lane acceptance runs its targeted tests AND this list. The members are
chosen for one property: each is a claim about the WHOLE tree that no single lane owns, and each is
cheap.

    uv run pytest -q $(uv run python tests/always_run.py) tests/test_<your module>.py

WHAT THIS FILE IS NOT. It is not an exemption mechanism and it holds no allowlist. It adds a second
place a cross-cutting test is named; it never excuses one, shortens one or skips one.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

#: The members, each with the category it covers and why it belongs. Paths relative to the repo root.
MEMBERS: tuple[dict[str, str], ...] = (
    {
        "path": "tests/test_engine_boundary.py",
        "category": "engine/licensing boundary",
        "why": (
            "the Apache-2.0 engine may not import the AGPL-3.0 application (LICENSING.md D40, "
            "Albert's split of 2026-09-11). A legal boundary, not a style rule. It is a STATIC "
            "analysis: it walks `ast.Import`/`ast.ImportFrom` over the source on disk and imports "
            "nothing, so its verdict is about the bytes of the files in the tree"
        ),
    },
    {
        "path": "tests/test_engine_package.py",
        "category": "engine/licensing boundary",
        "why": (
            "the same boundary from the packaging side: it builds the engine package and refuses a "
            "docstring that still names an application module, because in a package of its own that "
            "reference points at something that is not there"
        ),
    },
    {
        "path": "tests/test_headline_registry_matches_readme.py",
        "category": "README-registry match",
        "why": "a number quoted in the README must be the number the registry holds",
    },
    {
        "path": "tests/test_committed_data.py",
        "category": "committed-artefact presence",
        "why": (
            "a committed artefact's absence is a defect and must FAIL rather than skip "
            "(`must_be_committed`, `4f44dbf`). A guard that skips on it measures the machine"
        ),
    },
    {
        "path": "tests/test_results_writers_guard.py",
        "category": "the write guard",
        "why": (
            "nothing writes into data/results/ except through `save_result`. The coordinator named "
            "this category and not a file; this lane resolved it to this path. The other candidate "
            "in the tree is `tests/test_rebuild_write_guard.py`, which guards a narrower thing - a "
            "rebuild writing to the stores it verifies against - and is named here so a reader can "
            "correct the choice rather than guess at it"
        ),
    },
)

#: Just the paths, in order.
ALWAYS_RUN: tuple[str, ...] = tuple(m["path"] for m in MEMBERS)

#: MEASURED on 2026-10-02 in a fresh interpreter per module, by importing the module and reading
#: `sys.modules`, NOT by grepping its import lines. The two answers disagree: at file level none of
#: these names `genomeos.manifest` or `genomeos.results`, and transitively a module that imports any
#: part of the application pulls both in. Only the transitive answer is true, and it is the one a
#: declared main-tree deviation has to quote, because `genomeos/manifest.py` and `genomeos/results.py`
#: are held uncommitted by a peer.
#:
#: Every ALWAYS_RUN member loads ZERO genomeos modules, so no member's verdict can be affected by an
#: uncommitted application file. `tests/test_rule_cell_provenance.py` loads 48, including both.
TRANSITIVE_GENOMEOS_MODULES_LOADED: dict[str, int] = {
    "tests/test_engine_boundary.py": 0,
    "tests/test_engine_package.py": 0,
    "tests/test_headline_registry_matches_readme.py": 0,
    "tests/test_committed_data.py": 0,
    "tests/test_results_writers_guard.py": 0,
}

HOW_THE_COUNT_WAS_TAKEN = (
    "one fresh interpreter per module: `importlib.import_module(name)` then every key of "
    "`sys.modules` under `genomeos`. A grep of the file's own import lines gives a different and "
    "wrong answer, because what matters is what the imports pull in"
)


def test_every_member_exists_and_is_a_test_module() -> None:
    for member in MEMBERS:
        path = ROOT / member["path"]
        assert path.exists(), f"{member['path']} is named in ALWAYS_RUN and is not in the tree"
        assert "def test_" in path.read_text(), member["path"]
        assert member["category"].strip() and member["why"].strip(), member["path"]


def test_the_four_categories_the_coordinator_named_are_all_covered() -> None:
    categories = {m["category"] for m in MEMBERS}
    assert categories == {
        "engine/licensing boundary",
        "README-registry match",
        "committed-artefact presence",
        "the write guard",
    }


def test_this_list_carries_no_allowlist_and_no_skip() -> None:
    """It names tests; it never excuses one. A member that could be skipped would defeat its purpose."""
    import ast

    tree = ast.parse(Path(__file__).read_text())
    offenders = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = [node.target] if isinstance(node, ast.AnnAssign) else node.targets
            offenders += [
                ast.unparse(t)
                for t in targets
                if "allow" in ast.unparse(t).lower() or "exempt" in ast.unparse(t).lower()
            ]
        if isinstance(node, ast.FunctionDef):
            offenders += [ast.unparse(d) for d in node.decorator_list if "skip" in ast.unparse(d)]
        if isinstance(node, ast.Call) and "skip" in ast.unparse(node.func):
            offenders.append(ast.unparse(node.func))
    assert offenders == [], offenders
    assert all((ROOT / member).exists() for member in ALWAYS_RUN)


def test_the_transitive_measurement_covers_every_member() -> None:
    assert set(TRANSITIVE_GENOMEOS_MODULES_LOADED) == set(ALWAYS_RUN)
    assert HOW_THE_COUNT_WAS_TAKEN.strip()


if __name__ == "__main__":
    print(" ".join(ALWAYS_RUN))
