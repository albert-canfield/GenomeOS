# SPDX-License-Identifier: AGPL-3.0-or-later
"""`scripts/inert_skip_guards.py`: the finder that resolves what a skip guard is keyed on.

WHAT THIS CHECKS AND WHY IT IS NOT THE ONE ALREADY IN THE SUITE.
`tests/test_local_data_symlinked_store.py` carries `test_no_skip_guard_is_keyed_on_a_path_that_git_
TRACKS`, which reads `skipif` DECORATORS and nothing else. It was right about its own population and
blind to the larger one: this suite has 88 `skipif` lines and 78 `pytest.skip()` calls, and a
`pytest.skip()` under an `if not P.exists()` inside a test body is the identical defect in a shape an
`args[0]` walk never reaches. Measured on 2026-10-03 against the real suite: that invariant passes
with 31 guards keyed on a path git TRACKS still in the tree, 7 of them `skipif` decorators it does
read, because its resolver cannot follow `Path(__file__).resolve().parent.parent`, a constant imported
by name, or a dict entry. The finder here resolves all three, so it sees them.

NOTHING HERE WEAKENS THAT TEST. It keeps its own resolver, its own plants and its own non-vacuity
pins. This is the wider instrument beside it, and the two disagreeing about a guard is a finding
either way.

THE PLANTS ARE THE EVIDENCE. Each one is written as source and read, so the finder is shown to fail on
a case it must catch rather than asserted to be able to:
  - a guard on a TRACKED path is named INERT, in both shapes (decorator and body);
  - the same guard on a path in a git-ignored store is NOT called inert -- it is LIVE;
  - a path the finder cannot resolve is UNRESOLVED and is NOT folded into either answer, which is the
    regression this finder actually had: `<?>/chr21.json.gz` is untracked only because nothing can
    match a hole, and reading that as "the guard can fire" made an unknown guard look live.
"""

from __future__ import annotations

import ast
import importlib.util
import subprocess
import sys
from functools import cache
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
FINDER = ROOT / "scripts" / "inert_skip_guards.py"


@cache
def _finder():
    """The script, loaded as a module. It is a script and not a package, as its siblings are.

    Registered in `sys.modules` before it is executed, because its `@dataclass` resolves a string
    annotation through `sys.modules[cls.__module__]` and an unregistered module makes that `None`.
    """
    spec = importlib.util.spec_from_file_location("inert_skip_guards", FINDER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


#: a path git tracks, so a guard on its absence can never fire
A_TRACKED_RESULT = "data/results/clause1.json"
#: a path in a git-ignored store, so a guard on its absence is doing real work
A_LOCAL_STORE_FILE = "data/cache/entex/alphagenome_track_metadata_copy.csv"


def _classify(source: str, name: str = "planted.py"):
    """{where: (class, terms)} for one planted module, read against the real tracked set."""
    finder = _finder()
    tracked = finder.tracked_paths()
    where = f"tests/{name}"
    with finder.reading(where):
        guards = finder._guards(ast.parse(source, filename=where), where, tracked)
    return {g.where: (g.verdict, g.terms) for g in guards}


def test_the_finder_names_an_inert_skipif_planted_on_a_tracked_path() -> None:
    """Mutation leg one, the decorator shape. The plant must be NAMED, not merely counted."""
    source = (
        "from pathlib import Path\n"
        "import pytest\n"
        f'NEEDED = Path("{A_TRACKED_RESULT}")\n'
        '@pytest.mark.skipif(not NEEDED.exists(), reason="not written yet")\n'
        "def test_reads_a_committed_result():\n"
        "    assert NEEDED.read_text()\n"
    )
    found = _classify(source)
    assert found["tests/planted.py:4"] == ("INERT", (A_TRACKED_RESULT,)), found


def test_the_finder_names_an_inert_skip_INSIDE_A_BODY_the_shape_the_other_check_misses() -> None:
    """Mutation leg one, the body shape: `pytest.skip()` under an `if`, keyed through a local name.

    The condition is on the `if` and not on the call, and the path is bound one statement earlier from
    a `Path(__file__).resolve().parent.parent` root. Reading the call alone -- which is what a walk
    over `skipif(...)` arguments does -- classifies this as carrying no path at all.
    """
    source = (
        "from pathlib import Path\n"
        "import pytest\n"
        "ROOT = Path(__file__).resolve().parent.parent\n"
        "def test_reads_a_committed_result():\n"
        f'    p = ROOT / "{A_TRACKED_RESULT}"\n'
        "    if not p.exists():\n"
        '        pytest.skip("the result is not on this machine")\n'
        "    assert p.read_text()\n"
    )
    found = _classify(source)
    assert found["tests/planted.py:7"] == ("INERT", (A_TRACKED_RESULT,)), found


def test_the_finder_does_NOT_call_a_guard_on_a_git_ignored_store_inert() -> None:
    """Mutation leg two. The same two shapes keyed on `data/cache` are LIVE, and a finder that called
    them inert would have had a lane delete the guards the suite genuinely needs."""
    source = (
        "from pathlib import Path\n"
        "import pytest\n"
        "ROOT = Path(__file__).resolve().parent.parent\n"
        f'LOCAL = ROOT / "{A_LOCAL_STORE_FILE}"\n'
        '@pytest.mark.skipif(not LOCAL.exists(), reason="not on this machine")\n'
        "def test_reads_a_store():\n"
        "    assert LOCAL.read_text()\n"
        "def test_reads_the_store_in_a_body():\n"
        f'    q = ROOT / "{A_LOCAL_STORE_FILE}"\n'
        "    if not q.exists():\n"
        '        pytest.skip("not on this machine")\n'
        "    assert q.read_text()\n"
    )
    found = _classify(source)
    assert found["tests/planted.py:5"] == ("LIVE", (A_LOCAL_STORE_FILE,)), found
    assert found["tests/planted.py:11"] == ("LIVE", (A_LOCAL_STORE_FILE,)), found


def test_a_guard_the_finder_cannot_resolve_is_UNRESOLVED_and_not_folded_into_either_answer() -> None:
    """The regression this finder HAD, and the reason UNRESOLVED is a class and not a footnote.

    `data/results/<?>.json` is not tracked -- nothing can match a hole -- so a classifier that asked
    only "is it tracked?" called it LIVE and reported a guard of unknown class as one doing real work.
    Both plants here are unknown, and both must say so: one with an opaque leaf under a directory that
    holds committed and machine-local files alike, one with no resolvable path at all.
    """
    source = (
        "from pathlib import Path\n"
        "import pytest\n"
        "ROOT = Path(__file__).resolve().parent.parent\n"
        "def test_leaf_from_a_fixture(mod):\n"
        '    p = ROOT / "data" / "results" / f"{mod.NAME}.json"\n'
        "    if not p.exists():\n"
        '        pytest.skip("not here")\n'
        "def test_nothing_resolves(mod):\n"
        "    if not mod.whatever().exists():\n"
        '        pytest.skip("not here")\n'
    )
    found = _classify(source)
    assert found["tests/planted.py:7"] == ("UNRESOLVED", ("data/results/<?>.json",)), found
    assert found["tests/planted.py:10"] == ("UNRESOLVED", (None,)), found


def test_a_guard_whose_second_term_can_be_absent_is_LIVE_and_is_not_converted() -> None:
    """One term that can fire is enough. A compound guard is LIVE even with a tracked term in it."""
    source = (
        "from pathlib import Path\n"
        "import pytest\n"
        "ROOT = Path(__file__).resolve().parent.parent\n"
        f'A = ROOT / "{A_TRACKED_RESULT}"\n'
        f'B = ROOT / "{A_LOCAL_STORE_FILE}"\n'
        '@pytest.mark.skipif(not (A.exists() and B.exists()), reason="needs both")\n'
        "def test_needs_both():\n"
        "    assert A.read_text() and B.read_text()\n"
    )
    found = _classify(source)
    assert found["tests/planted.py:6"][0] == "LIVE", found
    assert set(found["tests/planted.py:6"][1]) == {A_TRACKED_RESULT, A_LOCAL_STORE_FILE}, found


def test_the_project_marker_is_never_inert_because_it_polices_itself() -> None:
    """`needs_local_data` reaches `local_data.check_is_machine_local`, which asks git and RAISES on a
    tracked path. It cannot be an inert guard, and the finder says so by its own class rather than by
    resolving the path and then guessing what the marker would have done with it."""
    source = (
        "import pytest\n"
        f'@pytest.mark.needs_local_data("{A_LOCAL_STORE_FILE}", how="fetch it")\n'
        "def test_reads_a_store():\n"
        "    pass\n"
    )
    found = _classify(source)
    assert found["tests/planted.py:2"] == ("LIVE_BY_RULE", (A_LOCAL_STORE_FILE,)), found


def test_a_tracked_path_tested_for_PRESENCE_is_its_own_class_and_not_called_inert() -> None:
    """`if P.exists(): pytest.skip(...)` on a tracked path skips EVERY run, so the test never runs at
    all. That is worse than a guard that never fires and it is not reported as the same thing. The
    real suite has none; this plant is why the class exists and the only place it is exercised."""
    source = (
        "from pathlib import Path\n"
        "import pytest\n"
        f'NEEDED = Path("{A_TRACKED_RESULT}")\n'
        "def test_skips_whenever_the_result_is_there():\n"
        "    if NEEDED.exists():\n"
        '        pytest.skip("there is a result, so this is not the case under test")\n'
    )
    found = _classify(source)
    assert found["tests/planted.py:6"] == ("ALWAYS_FIRES", (A_TRACKED_RESULT,)), found


def test_the_two_paths_the_plants_stand_on_are_what_this_test_says_they_are() -> None:
    """Non-vacuity of the plants themselves: if `clause1.json` stopped being tracked, or `data/cache`
    stopped being ignored, every plant above would still pass while meaning the opposite."""
    finder = _finder()
    tracked = finder.tracked_paths()
    assert A_TRACKED_RESULT in tracked
    assert A_LOCAL_STORE_FILE not in tracked
    assert finder.under_ignored_store(A_LOCAL_STORE_FILE)
    ignored = subprocess.run(
        ["git", "check-ignore", "-q", "--no-index", "data/cache"], capture_output=True, cwd=ROOT
    )
    assert ignored.returncode == 0, "data/cache is no longer git-ignored, so leg two proves nothing"


def test_the_finder_resolves_paths_over_the_REAL_suite_and_is_not_vacuous() -> None:
    """A resolver that answered None everywhere would pass every plant that asserts a class ON a
    plant. Over the real suite it must still resolve real paths in real numbers."""
    finder = _finder()
    guards = finder.all_guards(ROOT / "tests")
    resolved = [g for g in guards if any(t is not None and not finder.has_hole(t) for t in g.terms)]
    assert len(guards) >= 150, f"the finder has stopped finding guards: {len(guards)}"
    assert len(resolved) >= 60, f"the finder has stopped resolving paths: {len(resolved)}"
    assert {g.verdict for g in guards} >= {"LIVE", "LIVE_BY_RULE", "NOT_PATH"}


def test_the_finder_runs_as_a_command_and_prints_its_counts() -> None:
    """It is a tool a lane runs, not only a library this test imports."""
    out = subprocess.run(["python", str(FINDER)], capture_output=True, text=True, cwd=ROOT, check=True)
    assert "LIVE=" in out.stdout and "UNRESOLVED=" in out.stdout, out.stdout


@pytest.mark.parametrize("shape", ["needs_local_data", "importorskip", "mark.skip"])
def test_each_shape_the_finder_knows_is_present_in_the_real_suite(shape: str) -> None:
    """A corpus pin: the shapes are read off the suite, so a shape going unread is visible."""
    finder = _finder()
    assert any(g.shape == shape for g in finder.all_guards(ROOT / "tests")), shape


def test_a_condition_that_CALLS_a_helper_is_read_through_the_helper_not_filed_as_NOT_PATH() -> None:
    """The one way NOT_PATH could hide an inert guard, found on 2026-10-03 and closed.

    `if _benchmark_present() is None: pytest.skip(...)` -- a real guard in
    `tests/test_increase_links.py` -- carries no existence call and no path, so a reading that stops
    at the condition files it as "the condition tests no path" while the whole question lives in the
    three lines of the helper. It is read through the helper now, and the real case it was found on
    comes back LIVE on `data/knowledge/crispri/<?>`, which is what it is.
    """
    source = (
        "from pathlib import Path\n"
        "import pytest\n"
        "ROOT = Path(__file__).resolve().parent.parent\n"
        "def _present():\n"
        f'    if not (ROOT / "{A_TRACKED_RESULT}").exists():\n'
        "        return None\n"
        "    return True\n"
        "def test_reads_it():\n"
        "    if _present() is None:\n"
        '        pytest.skip("not here")\n'
    )
    found = _classify(source)
    assert found["tests/planted.py:10"] == ("INERT", (A_TRACKED_RESULT,)), found
    finder = _finder()
    live = {
        g.where: g.terms
        for g in finder.all_guards(ROOT / "tests")
        if "test_increase_links.py" in g.where and g.verdict == "LIVE"
    }
    assert live, "the real case this was found on no longer reads as LIVE"
    assert any("data/knowledge" in t for terms in live.values() for t in terms if t), live


def test_NO_skip_guard_IN_THE_SUITE_is_keyed_on_a_path_that_git_TRACKS() -> None:
    """The deliverable, over the real suite and over all four shapes a guard comes in.

    31 guards were in this class on 2026-10-03 and every one of them now raises instead. A tracked
    file is in every checkout of its commit, so a skip waiting on its absence is not a guard: it is a
    sentence that reads as protection while the condition is a constant. The companion check in
    `tests/test_local_data_symlinked_store.py` makes the same claim over `skipif` decorators alone
    and keeps its own resolver; this one covers `pytest.skip()` in a body and a fixture as well,
    which is where 24 of the 31 lived.

    SCOPED TO TEST FILES GIT TRACKS, and that is the claim correctly aimed rather than a weaker one.
    Several sessions edit this one checkout at once, so an untracked `tests/test_*.py` is a peer's
    work in progress. Reading it would turn one lane's half-written file into every other lane's red
    push, which happened within the hour of this being written: `tests/test_manifest_rebuild_reads.py`
    appeared untracked carrying `if not Path("data/results").is_dir(): pytest.skip("data/results is
    git-ignored: ...")`, whose premise is the one `tests/local_data.py` exists to correct --
    `data/results` is ignored by PATTERN with `!` re-includes and hundreds of its files are committed.
    It was reported to the coordinator rather than edited, and this check will refuse it the moment it
    is committed, which is exactly when it becomes part of the suite.
    """
    finder = _finder()
    inert = [
        g for g in finder.all_guards(ROOT / "tests", root=ROOT, tracked_only=True) if g.verdict == "INERT"
    ]
    assert not inert, (
        "these skip guards test for a file git TRACKS, so they are present in every checkout and the "
        "skip can never fire. Convert each to tracked_paths.must_be_present, which RAISES by name: "
        f"{[(g.where, g.terms) for g in inert]}"
    )


def test_no_test_is_SKIPPED_EVERY_RUN_by_a_guard_on_a_tracked_path() -> None:
    """The mirror class, over the real suite: a tracked path tested for PRESENCE never runs at all.
    There are none, and this is the check that says so rather than the absence of a report."""
    finder = _finder()
    always = [
        g
        for g in finder.all_guards(ROOT / "tests", root=ROOT, tracked_only=True)
        if g.verdict == "ALWAYS_FIRES"
    ]
    assert not always, f"these tests skip on every run, so they never run: {always}"


#: THE ELEVEN UNRESOLVED GUARDS, read by hand on 2026-10-03 and left alone, with what each one is.
#: Reported by name rather than folded into a class, because an UNRESOLVED count of zero would have
#: to be earned and this finder has not earned it. Every one of them is LIVE or needs a judgement no
#: resolver makes, and NOT ONE is keyed on a tracked path:
#:   tests/conftest.py:47                      `local_data.missing(tuple(mark.args))` -- the
#:                                             needs_local_data machinery itself, whose args are a
#:                                             runtime marker. It cannot be inert: the helper raises.
#:   tests/test_clause2_design_power.py:293/475/656   `dp.RESULTS_DIR / f"{dp.RESULT}.json"` reached
#:                                             through a module loaded by importlib from a path.
#:   tests/test_crispri_split.py:66/242        `measured.CRISPRI_KNOWLEDGE` under data/knowledge.
#:   tests/test_joint_pretest.py:115           the same store, through `dict(measured.CRISPRI_SPLIT_OF)`.
#:   tests/test_kd_semantics_census.py:50      `m.CACHE / m.DOCS[...]["file"]`, m a dynamic import.
#:   tests/test_loci_gene_input.py:220         a result name from a loop over a runtime frame.
#:   tests/test_pre_push_hook_is_the_wrapper.py:49   `.git/hooks/pre-push`, which is not in the tree
#:                                             at all, so no `ls-files` answer applies to it.
#:   tests/test_response_map_increment4.py:59  `ROOT / i4.INCREMENT_3_COUNT`, i4 a fixture.
