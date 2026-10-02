# SPDX-License-Identifier: AGPL-3.0-or-later
"""One `counting_path` and one `code_cleanliness`, and the copies cannot come back.

Before 2026-10-02 `code_cleanliness` was copy-pasted into ten files and the import closure behind it
into eight. Nothing checked that the copies agreed, and they did not. Measured against the same entry
script on the same tree, the copies fell into three groups:

    41 paths   six copies that followed only the `genomeos` package
    68 paths   `scripts/crispri_published_v2.py`, which followed `scripts.*` and parent `__init__.py`
    44 paths   `scripts/cell2_eligibility.py`, which had never taken the cd263bc fix and so admitted
               `genomeos/genome/Genome.py`, `genomeos/genome/Annotation.py` and
               `genomeos/knowledge/Reactome.py` -- class names that no file spells -- because this
               filesystem is case-insensitive

Every published result's honesty about which code was uncommitted when it was written rests on that
block. These tests are what keep the copies from coming back: the copies were the defect.
"""

from __future__ import annotations

import ast
import importlib
import sys
from pathlib import Path

import pytest

from genomeos import manifest as mf

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

#: The ten files that each held a copy, as the writers that must now all read the same implementation.
#: `genomeos/attribution/context_evidence.py` is reached as a module; the other nine as scripts.
CALLERS = (
    "scripts.astroreg_register",
    "scripts.response_map_coverage",
    "scripts.crispri_published_v2",
    "scripts.placement_cause_198",
    "scripts.fresh_crispri_eligibility",
    "scripts.crispri_benchmark_v2",
    "scripts.re2g_paired",
    "scripts.context_evidence_census",
    "scripts.cell2_eligibility",
    "genomeos.attribution.context_evidence",
)

#: The one file allowed to define them. Anything else is a copy.
SHARED_HOME = "genomeos/manifest.py"
FORBIDDEN = ("counting_path", "code_cleanliness")


# (d) the copies cannot come back ---------------------------------------------------------------------


def test_no_file_outside_the_shared_home_defines_its_own_counting_path_or_cleanliness():
    """The supervisor's rule, 2026-10-02: the copies are the defect, this test is what keeps them out.

    A re-export (`code_cleanliness = _mf.code_cleanliness`) is not a definition and is how a module
    that writers already reach through keeps its name; a `def` is a second implementation.
    """
    offenders = []
    for path in sorted([*ROOT.glob("genomeos/**/*.py"), *ROOT.glob("scripts/**/*.py")]):
        rel = path.relative_to(ROOT).as_posix()
        if rel == SHARED_HOME:
            continue
        try:
            tree = ast.parse(path.read_text())
        except (OSError, SyntaxError):
            continue
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef | ast.AsyncFunctionDef) and node.name in FORBIDDEN:
                offenders.append(f"{rel}:{node.lineno} def {node.name}")
    assert offenders == [], (
        f"a second implementation is a drift waiting to happen; import genomeos.manifest instead: {offenders}"
    )


def test_the_shared_home_does_define_them_so_the_test_above_can_fail():
    """A guard that passes because it looked in the wrong place is worse than no guard."""
    tree = ast.parse((ROOT / SHARED_HOME).read_text())
    defined = {n.name for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name in FORBIDDEN}
    assert defined == set(FORBIDDEN)


# (c) the ten callers all report the same shape -------------------------------------------------------


@pytest.fixture(scope="module")
def blocks() -> dict[str, dict]:
    """Each former copy-holder's cleanliness block, computed through the shared function."""
    out = {}
    for name in CALLERS:
        m = importlib.import_module(name)
        own = getattr(m, "OWN_CODE", ())
        entry = getattr(m, "ENTRY", None)
        if entry is None:  # the module, which is a library and not an entry script
            entry = "scripts/context_evidence_census.py"
        out[name] = mf.code_cleanliness(entry, own, ROOT)
    return out


def test_every_caller_returns_the_same_keys(blocks):
    """What the copies drifted on first. One implementation means one key set."""
    assert set(blocks) == set(CALLERS), "every caller must be reachable"
    for name, block in blocks.items():
        assert set(block) == set(mf.CLEANLINESS_KEYS), name


def test_every_caller_returns_the_same_wording_for_the_fields_that_do_not_vary(blocks):
    """The copies had four different notes and one had none, so a reader of two results was told two
    different things about what the block meant."""
    for field in ("counting_path_is_computed", "note"):
        said = {block[field] for block in blocks.values()}
        assert len(said) == 1, f"{field} is worded {len(said)} ways"


def test_the_callers_that_write_a_result_name_their_own_entry(blocks):
    """`ENTRY` is the script whose closure is the counting path. A writer that did not name one would
    be guessing, which is what the shared function refuses to do."""
    for name in CALLERS:
        m = importlib.import_module(name)
        if name == "genomeos.attribution.context_evidence":
            continue  # a library: its callers pass their own entry
        assert name.replace(".", "/") + ".py" == m.ENTRY, name
        assert m.ENTRY in mf.counting_path(m.ENTRY, ROOT), f"{name} must be on its own counting path"


# what the unified closure is, and the two defects it closes ------------------------------------------


def test_a_class_name_can_never_enter_the_counting_path_as_a_file():
    """Finding 2, the worst of the three: `data/results/cell2_eligibility.json` published
    `genomeos/genome/Genome.py`, `genomeos/genome/Annotation.py` and `genomeos/knowledge/Reactome.py`.
    The files are genome.py, annotation.py and reactome.py; those three are the classes they export,
    admitted by a plain `Path.is_file()` on a case-insensitive filesystem. The result's counting path
    is three entries shorter where the filesystem is case-sensitive, and a counting-path difference is
    a real difference, not an environment field (tests/test_manifest_rebuild_environment.py).
    """
    path = mf.counting_path("scripts/cell2_eligibility.py", ROOT)
    for not_a_file in (
        "genomeos/genome/Genome.py",
        "genomeos/genome/Annotation.py",
        "genomeos/knowledge/Reactome.py",
    ):
        assert not_a_file not in path
    assert len({p.lower() for p in path}) == len(path), "no two entries differ only by case"
    for rel in path:
        assert mf._is_file_exactly(ROOT, rel.split("/")), f"{rel} is not spelled as its directory spells it"


def test_a_peer_reached_through_sys_path_is_on_the_counting_path():
    """Finding 3, turned into a test at the supervisor's direction. `scripts/placement_cause_198.py`
    does `sys.path.insert(0, str(ROOT / "scripts"))` and then `import response_map_coverage as rmc`,
    whose `attaches_to` defines the 198 the result is about. Its published counting path of 45 entries
    did not name that file, so the block's own note -- a foreign file off the path cannot have entered
    the count -- did not cover a file whose code does enter it. No copy caught this, including the
    broad one; resolving bare imports against the directories the source puts on sys.path is what does.
    """
    path = mf.counting_path("scripts/placement_cause_198.py", ROOT)
    assert "scripts/response_map_coverage.py" in path
    assert "scripts/placement_cause_198.py" in path


def test_the_closure_resolves_a_literal_sys_path_directory(tmp_path):
    """The mechanism of the test above, on a tree of three files rather than on the repository."""
    (tmp_path / "scripts").mkdir()
    (tmp_path / "entry.py").write_text(
        "import sys\nfrom pathlib import Path\n"
        "ROOT = Path(__file__).resolve().parent\n"
        'sys.path.insert(0, str(ROOT / "scripts"))\n'
        "import peer\n"
    )
    (tmp_path / "scripts" / "peer.py").write_text("import helper\n")
    (tmp_path / "scripts" / "helper.py").write_text("x = 1\n")
    assert mf.counting_path("entry.py", tmp_path) == [
        "entry.py",
        "scripts/helper.py",
        "scripts/peer.py",
    ]


def test_a_sys_path_target_that_is_not_literal_is_skipped(tmp_path):
    """`sys.path.insert(0, sys.argv[1])` names no directory at write time, so it adds nothing."""
    (tmp_path / "entry.py").write_text("import sys\nsys.path.insert(0, sys.argv[1])\nimport peer\n")
    (tmp_path / "peer.py").write_text("x = 1\n")
    assert mf.counting_path("entry.py", tmp_path) == ["entry.py", "peer.py"]


def test_parent_and_parents_are_told_apart(tmp_path):
    """`Path(__file__).resolve().parent` is the file's own directory; `.parents[1]` is one above it.
    Reading the expression rather than its text is what distinguishes them, and it has to be read: the
    two differ by one directory and a script that inserts the wrong one imports nothing.

    `only_in_pkg` exists only below `pkg/`, so it is reachable from the first entry and not the second.
    The repository root is always searched, since that is where a writer is run from, so a name that
    resolves in both places resolves to both: a counting path is a claim about what could have entered
    a number, and where two files could have, the honest answer names both.
    """
    (tmp_path / "pkg").mkdir()
    (tmp_path / "pkg" / "only_in_pkg.py").write_text("x = 1\n")
    (tmp_path / "pkg" / "own_dir.py").write_text(
        "import sys\nfrom pathlib import Path\n"
        "sys.path.insert(0, str(Path(__file__).resolve().parent))\nimport only_in_pkg\n"
    )
    (tmp_path / "pkg" / "one_above.py").write_text(
        "import sys\nfrom pathlib import Path\n"
        "sys.path.insert(0, str(Path(__file__).resolve().parents[1]))\nimport only_in_pkg\n"
    )
    assert mf.counting_path("pkg/own_dir.py", tmp_path) == ["pkg/only_in_pkg.py", "pkg/own_dir.py"]
    assert mf.counting_path("pkg/one_above.py", tmp_path) == ["pkg/one_above.py"]


def test_the_closure_counts_the_init_of_every_package_above_a_module(tmp_path):
    (tmp_path / "pkg" / "sub").mkdir(parents=True)
    for rel in ("pkg/__init__.py", "pkg/sub/__init__.py"):
        (tmp_path / rel).write_text("")
    (tmp_path / "pkg" / "sub" / "thing.py").write_text("x = 1\n")
    (tmp_path / "entry.py").write_text("import pkg.sub.thing\n")
    assert mf.counting_path("entry.py", tmp_path) == [
        "entry.py",
        "pkg/__init__.py",
        "pkg/sub/__init__.py",
        "pkg/sub/thing.py",
    ]


def test_a_module_outside_the_repository_stays_off_the_counting_path(tmp_path):
    (tmp_path / "entry.py").write_text("import json\nimport numpy\nfrom typing import Any\n")
    assert mf.counting_path("entry.py", tmp_path) == ["entry.py"]


def test_the_cleanliness_block_sets_the_uncommitted_files_against_the_counting_path():
    block = mf.code_cleanliness("scripts/context_evidence_census.py", ("a.py",), ROOT)
    assert block["own_uncommitted_code"] == []  # "a.py" is not a file this repository has
    assert block["counting_path_count"] == len(block["counting_path"])
    for p in block["foreign_uncommitted_code_on_the_counting_path"]:
        assert p in block["counting_path"] and p in block["foreign_uncommitted_code"]


def test_the_entry_and_the_own_code_are_arguments_and_neither_is_guessed():
    """The defect the copies shared: each closed over its own module globals, so the same function name
    meant a different thing in every file."""
    import inspect

    params = list(inspect.signature(mf.code_cleanliness).parameters)
    assert params[:2] == ["entry", "own_code"]
    a = mf.code_cleanliness("scripts/cell2_eligibility.py", (), ROOT)
    b = mf.code_cleanliness("scripts/re2g_paired.py", (), ROOT)
    assert a["counting_path"] != b["counting_path"], "the entry must decide the counting path"
