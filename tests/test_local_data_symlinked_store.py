# SPDX-License-Identifier: AGPL-3.0-or-later
"""A store reached through a symbolic link must SKIP BY NAME, and never ERROR.

WHAT THIS PLANTS. On 2026-10-02 twenty tests ERRORED and two RED status files were written for trees
that were nothing of the kind -- `status-7446015a503fcc6f281173590591709d1f6e7aab` (4,418 passed, 0
failed, 20 errors) and `status-579b4949ea024ef98722a9ea14f9825478818724` (4,436 passed, 0 failed, 20
errors, which refused a push). Every one of the forty was the same line:

    local_data.NotMachineLocalError: git could not say whether
    data/cache/entex/alphagenome_track_metadata_copy.csv is ignored
    (fatal: pathspec '...' is beyond a symbolic link)

THE CAUSE, measured rather than reasoned, because two diagnoses were offered and both were wrong.
`.git/hooks/pre-push` was a STALE COPY taken at 04:34 which still did
`ln -s "$root/data/$d" "$tmp/data/$d"`, while `scripts/pre-push.sh` had moved that day to read-only
APFS clones (`scripts/link_stores_read_only.py`). The hook is a copy and not a symlink -- see
`scripts/install-hooks.sh` -- so editing the script in the tree changed nothing the push ran. The
push's worktree therefore had `data/cache` as a symbolic link into the main checkout, and git stops
at a symbolic link.

NOT the cause, and both were checked before this test was written. `$TMPDIR` lies behind
`/var -> private/var` on macOS, but a worktree under `$TMPDIR` with the stores as real directories
answers 0: git resolves a linked worktree's root physically and a relative pathspec is taken from the
process's physical cwd. And `realpath` before `check-ignore` does not fix it -- it asks about a
different NAME, and in the exact failing shape the resolved path lies in another working tree of the
same repository, which git refuses with 128 again. Both readings are in `local_data.ignored_by_name`.

WHAT THESE TESTS ASSERT, which is the difference between a fix and a cover-up. A test that needs a
store it cannot find must SKIP WITH THE STORE NAMED AND THE FETCH COMMAND GIVEN. It must not pass, it
must not be quietly deselected, and it must not error. Twenty silent passes would have been worse than
the red, because then nobody could tell "not run here" from "passed".
"""

from __future__ import annotations

import ast
import subprocess
from pathlib import Path

import local_data
import pytest

#: The real store path from the forty errors, used verbatim so the planted case is the measured one.
FAILING_PATH = "data/cache/entex/alphagenome_track_metadata_copy.csv"


def _git(cwd: Path, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *args], cwd=str(cwd), capture_output=True, text=True)


def _repo(root: Path) -> Path:
    """A repository that ignores `data/cache`, with one commit, made without porcelain `git commit`."""
    root.mkdir(parents=True, exist_ok=True)
    _git(root, "init", "-q", "-b", "main", ".")
    (root / ".gitignore").write_text("data/cache\n")
    _git(root, "add", ".gitignore")
    tree = _git(root, "write-tree").stdout.strip()
    sha = _git(
        root,
        "-c",
        "user.email=t@example.invalid",
        "-c",
        "user.name=t",
        "commit-tree",
        tree,
        "-m",
        "init",
    ).stdout.strip()
    _git(root, "update-ref", "refs/heads/main", sha, "")
    return root


def test_git_itself_refuses_a_path_beyond_a_symlink_and_answers_for_a_real_directory(tmp_path) -> None:
    """The measurement the fix rests on, taken here so a git change cannot silently retire the fix."""
    root = _repo(tmp_path / "repo")
    (root / "store" / "entex").mkdir(parents=True)
    (root / "store" / "entex" / "x.csv").write_text("x\n")
    (root / "data").mkdir()
    (root / "data" / "cache").symlink_to(root / "store")

    linked = _git(root, "check-ignore", "-q", "--no-index", "data/cache/entex/x.csv")
    assert linked.returncode == 128, linked
    assert local_data.BEYOND_A_SYMLINK in linked.stderr

    (root / "data" / "cache").unlink()
    (root / "data" / "cache" / "entex").mkdir(parents=True)
    (root / "data" / "cache" / "entex" / "x.csv").write_text("x\n")
    real = _git(root, "check-ignore", "-q", "--no-index", "data/cache/entex/x.csv")
    assert real.returncode == 0, real


def test_a_symlinked_store_is_answered_by_name_and_does_not_raise(tmp_path, monkeypatch) -> None:
    """The twenty errors, planted: the same path, the same symlink, and now an ANSWER."""
    root = _repo(tmp_path / "repo")
    (root / "store" / "entex").mkdir(parents=True)
    (root / "store" / "entex" / "alphagenome_track_metadata_copy.csv").write_text("x\n")
    (root / "data").mkdir()
    (root / "data" / "cache").symlink_to(root / "store")
    monkeypatch.chdir(root)
    local_data.ignored_by_name.cache_clear()

    assert local_data.is_git_ignored(FAILING_PATH) is True
    ignored, answered_about = local_data.ignored_by_name(FAILING_PATH)
    assert (ignored, answered_about) == (True, "data/cache"), (ignored, answered_about)
    # and the marker check, which is what pytest's setup calls, no longer raises
    local_data.check_is_machine_local(FAILING_PATH)
    local_data.ignored_by_name.cache_clear()


def test_a_worktree_under_a_symlinked_tmpdir_was_never_the_cause(tmp_path, monkeypatch) -> None:
    """The second diagnosis, planted as the negative it is: a leading symlink above the repo is fine.

    `/var -> private/var` was proposed as the cause. Here the whole repository sits behind a symlinked
    parent, exactly as a `$TMPDIR` worktree does, and the stores are real directories: the question is
    answerable and the answer is "ignored". A fix aimed at this would have left the push red.
    """
    real = tmp_path / "physical"
    _repo(real / "repo")
    (real / "repo" / "data" / "cache" / "entex").mkdir(parents=True)
    (real / "repo" / "data" / "cache" / "entex" / "x.csv").write_text("x\n")
    (tmp_path / "behind").symlink_to(real)
    monkeypatch.chdir(tmp_path / "behind" / "repo")
    local_data.ignored_by_name.cache_clear()

    assert local_data.symlink_boundary("data/cache/entex/x.csv") is None
    assert local_data.ignored_by_name("data/cache/entex/x.csv") == (True, "data/cache/entex/x.csv")
    local_data.ignored_by_name.cache_clear()


def test_realpath_before_check_ignore_would_not_have_worked(tmp_path) -> None:
    """The first proposed fix, planted as the negative: resolving the path asks a DIFFERENT question.

    The failing shape exactly: a linked worktree whose `data/cache` is a symlink into the main
    checkout. The resolved path lies in another working tree of the same repository, and git refuses
    it with 128 a second time. Nothing here depends on git's wording beyond the code being 128.
    """
    root = _repo(tmp_path / "main")
    (root / "data" / "cache" / "entex").mkdir(parents=True)
    (root / "data" / "cache" / "entex" / "x.csv").write_text("x\n")
    worktree = tmp_path / "wt"
    added = _git(root, "worktree", "add", "-q", "--detach", str(worktree), "main")
    assert added.returncode == 0, added.stderr
    (worktree / "data").mkdir(exist_ok=True)
    (worktree / "data" / "cache").symlink_to(root / "data" / "cache")

    as_given = _git(worktree, "check-ignore", "-q", "--no-index", "data/cache/entex/x.csv")
    assert as_given.returncode == 128, as_given
    resolved = (worktree / "data" / "cache" / "entex" / "x.csv").resolve()
    as_resolved = _git(worktree, "check-ignore", "-q", "--no-index", str(resolved))
    assert as_resolved.returncode == 128, as_resolved
    # ... while the fix's question, asked at the symlink, is answerable
    at_boundary = _git(worktree, "check-ignore", "-q", "--no-index", "data/cache")
    assert at_boundary.returncode == 0, at_boundary


def test_a_store_missing_behind_a_symlink_skips_by_name_with_the_fetch_command(tmp_path, monkeypatch) -> None:
    """Not a pass and not an error: a SKIP whose reason names the store and how to get it."""
    root = _repo(tmp_path / "repo")
    (root / "data").mkdir()
    (root / "data" / "cache").symlink_to(root / "gone")  # a dangling link: the store is not here
    monkeypatch.chdir(root)
    local_data.ignored_by_name.cache_clear()

    absent = local_data.missing((FAILING_PATH,))
    assert absent == [FAILING_PATH]
    reason = local_data.skip_reason(absent, "scripts/entex_feasibility.py writes it (AG_METADATA)")
    assert "NOT RUN HERE, not passed" in reason
    assert FAILING_PATH in reason
    assert "scripts/entex_feasibility.py writes it (AG_METADATA)" in reason
    local_data.ignored_by_name.cache_clear()


def test_a_path_git_does_not_ignore_still_raises(tmp_path, monkeypatch) -> None:
    """The rule the fix must not weaken: a marker naming a committed path is wrong and says so."""
    root = _repo(tmp_path / "repo")
    monkeypatch.chdir(root)
    local_data.ignored_by_name.cache_clear()
    with pytest.raises(local_data.NotMachineLocalError, match="git does not ignore"):
        local_data.check_is_machine_local(".gitignore")
    local_data.ignored_by_name.cache_clear()


def test_a_symlink_that_is_not_an_ignored_name_is_refused_and_not_guessed(tmp_path, monkeypatch) -> None:
    """Truncating to the symlink answers only when that name IS ignored. Otherwise it refuses."""
    root = _repo(tmp_path / "repo")
    (root / "elsewhere").mkdir()
    (root / "notignored").symlink_to(root / "elsewhere")
    monkeypatch.chdir(root)
    local_data.ignored_by_name.cache_clear()
    with pytest.raises(local_data.NotMachineLocalError, match="stops at the symbolic link"):
        local_data.check_is_machine_local("notignored/deep/file.csv")
    local_data.ignored_by_name.cache_clear()


# ---- the marker invariant, read by SYNTAX and not by text ---------------------------------------
#
# WHY AN AST AND NOT A REGEX. Until 2026-10-03 this invariant was `re.finditer(r"needs_local_data\(")`
# over every `tests/test_*.py`, balancing parentheses by hand. A scan over source TEXT cannot tell a
# marker from a MENTION of one, and on 2026-10-03 it went red on a marker that exists only inside a
# STRING LITERAL -- `tests/test_fast_prepush.py`'s `MARKED` template, a test file that another test
# writes to disk. `c9509e1` patched that one literal as the smallest honest fix and recorded that the
# instrument was still wrong; the next plant that writes a marker into a string would have tripped it
# again. A syntax tree cannot see inside a string literal, which is exactly the property wanted.
#
# THE CLAIM IS UNCHANGED, and is the reason this exists: a skip reason without a `how=` names the
# store but not the way back to it, so a reader learns what is missing and not how to get it.
#
# COUNTS, measured at 0199fbe so a later divergence can be read against them: the text scan found 24
# call sites, this walk finds 23, and the one difference is the string literal above. The walk finds
# no site the text scan missed.


def markers_without_how(source: str, filename: str) -> list[str]:
    """Every `needs_local_data(...)` CALL in one module's syntax tree that gives no `how=` keyword.

    Resolves the callee by its final name, so all three forms in use are found: a direct
    `@pytest.mark.needs_local_data(...)` decorator, a module-level alias
    (`needs_x = pytest.mark.needs_local_data(...)`) and either of those stacked under another
    decorator -- a decorator is an expression like any other and `ast.walk` reaches it.
    """
    without = []
    for node in ast.walk(ast.parse(source, filename=filename)):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        if isinstance(func, ast.Attribute):
            name = func.attr
        elif isinstance(func, ast.Name):
            name = func.id
        else:
            continue
        if name != "needs_local_data":
            continue
        if not any(kw.arg == "how" for kw in node.keywords):
            without.append(f"{filename}:{node.lineno}")
    return without


def marker_call_sites(source: str, filename: str) -> list[int]:
    """The line of every `needs_local_data(...)` call in one module, found the same way."""
    lines = []
    for node in ast.walk(ast.parse(source, filename=filename)):
        if isinstance(node, ast.Call):
            func = node.func
            name = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", None)
            if name == "needs_local_data":
                lines.append(node.lineno)
    return sorted(lines)


def test_every_needs_local_data_marker_carries_the_command_that_fetches_the_store() -> None:
    """A skip reason without a `how` names the store but not the way back. None may be without one."""
    without: list[str] = []
    for path in sorted(Path("tests").glob("test_*.py")):
        without.extend(markers_without_how(path.read_text(), str(path)))
    assert without == [], (
        f"these needs_local_data markers give no `how=`, so their skip reason cannot name the command "
        f"that fetches the store: {without}"
    )


def test_the_marker_invariant_has_a_corpus_to_read_and_finds_all_three_decorator_forms() -> None:
    """Non-vacuity over the REAL suite: a walk that resolved nothing would pass the test above.

    The three forms are pinned by file, not by count, so adding a marker does not move this test:
    `tests/test_astroargmax.py` carries direct `@pytest.mark.needs_local_data(...)` decorators,
    `tests/test_compiled_defaults.py` a module-level alias, and in that same file the alias is
    STACKED with a second decorator.
    """
    found = {
        str(path): marker_call_sites(path.read_text(), str(path))
        for path in sorted(Path("tests").glob("test_*.py"))
    }
    total = sum(len(v) for v in found.values())
    assert total >= 20, f"the marker corpus has all but vanished, so the invariant reads nothing: {total}"
    assert len(found["tests/test_astroargmax.py"]) >= 8, found["tests/test_astroargmax.py"]
    assert found["tests/test_compiled_defaults.py"], "the module-level alias form was not found"
    stacked = Path("tests/test_compiled_defaults.py").read_text()
    assert "@needs_the_chr21_element_table\n@needs_chr21\n" in stacked, "the stacked form is gone"


def test_the_marker_invariant_REFUSES_a_real_marker_that_carries_no_how() -> None:
    """The plant. A check that cannot be shown to refuse is a list, not a check."""
    planted = (
        "import pytest\n"
        "\n"
        "\n"
        '@pytest.mark.needs_local_data("data/cache/planted")\n'
        "def test_reads_a_store_without_saying_how_to_get_it():\n"
        "    pass\n"
    )
    assert markers_without_how(planted, "planted.py") == ["planted.py:4"]
    # and the same marker, with a `how=`, is accepted: it is the `how=` being read, not the call
    with_how = planted.replace('"data/cache/planted")', '"data/cache/planted", how="fetch it")')
    assert markers_without_how(with_how, "planted.py") == []
    assert marker_call_sites(with_how, "planted.py") == [4]


def test_a_marker_written_inside_a_STRING_LITERAL_is_ignored() -> None:
    """The case that caused this, planted so nobody re-tightens the instrument back to a text scan.

    A test file that WRITES a test file carries markers in its string literals. Those are not call
    sites of this suite: the generated file has its own run, with its own conftest, and this
    invariant does not reach it. The planted literal here deliberately carries NO `how=`, so a text
    scan would refuse it; the AST sees a string and nothing else.
    """
    planted = (
        'TEMPLATE = """\\\n'
        "import pytest\n"
        "\n"
        "\n"
        '@pytest.mark.needs_local_data("data/cache/planted")\n'
        "def test_written_to_disk_by_another_test():\n"
        "    pass\n"
        '"""\n'
    )
    assert 'needs_local_data("data/cache/planted")' in planted  # the text a regex would have matched
    assert markers_without_how(planted, "planted.py") == []
    assert marker_call_sites(planted, "planted.py") == []


# ---- a skip guard keyed on a path git TRACKS can never fire -------------------------------------
#
# THE CLASS, found on 2026-10-03. Twelve `skipif` guards were keyed on the existence of a file that
# git TRACKS. A tracked file is in every checkout -- this one, CI, a verdict worktree, the store-free
# pre-push leg, which keeps `data/results` because it is committed -- so those skips had never fired
# anywhere and could not. They were not failures: the tests ran and passed. They read as protection
# that did not exist, and a reader who believed them would think a test was conditional when it was
# not. The twelve were settled (eleven conditions removed, one turned into a fixture that FAILS) and
# this is the check that stops the class coming back.
#
# WHAT IS NOT REFUSED, and why the rule is "every term", not "any term". A guard whose condition
# names several paths can fire as soon as ONE of them can be absent, and a term this cannot resolve
# statically -- an imported module constant, a function call -- might be anything. So a guard is
# called inert only when EVERY path term it carries resolves AND every one of them is tracked. The
# two compound guards in `tests/test_context_evidence.py` are of exactly that shape: inert on their
# first term (`data/results/budget_chr21.json`, tracked) and live on their second, and they pass.
#
# WHY NOT `needs_local_data`. The marker is for the git-IGNORED machine-local stores and RAISES on a
# path git tracks (`local_data.ignored_by_name`). Converting one of these guards to it would be a
# lie about the path and would raise at once.


def _tracked_paths() -> frozenset[str]:
    """Every path git tracks in this checkout, as git spells it: repository-relative, forward slashes."""
    out = subprocess.run(["git", "ls-files"], capture_output=True, text=True, check=True)
    return frozenset(out.stdout.splitlines())


def _is_repo_root(node: ast.AST) -> bool:
    """`Path(__file__).resolve().parents[1]` and its kin: the repository root spelled from a file."""
    return any(isinstance(n, ast.Attribute) and n.attr == "parents" for n in ast.walk(node))


def _as_path(node: ast.AST, consts: dict[str, str]) -> str | None:
    """A repository-relative path for an expression, or None when it cannot be resolved by syntax.

    None is not a failure: it is the honest answer for `ce.TRACK_METADATA` or `default_gencode(...)`,
    and it is what keeps such a guard out of the inert class.
    """
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    if isinstance(node, ast.Name):
        return consts.get(node.id)
    if isinstance(node, ast.Call):
        func = node.func
        name = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", None)
        if name in {"Path", "str"} and len(node.args) == 1 and not node.keywords:
            return _as_path(node.args[0], consts)
        return None
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
        left = _as_path(node.left, consts)
        if left is None and _is_repo_root(node.left):
            left = ""  # ROOT / "data" / ... is the same path as "data/..." asked of git
        right = _as_path(node.right, consts)
        if left is None or right is None:
            return None
        return f"{left}/{right}".lstrip("/")
    return None


def _existence_terms(node: ast.AST, consts: dict[str, str], predicates: dict[str, list]) -> list:
    """Every `<something>.exists()` term inside an expression, each resolved to a path or to None.

    Follows a module-level boolean through its name, which is the form `needs_chr21` used:
    `HAS_CHR21 = (ROOT / "data" / "results" / "budget_chr21.json").exists()` and then
    `skipif(not HAS_CHR21, ...)`. Without that the condition is a bare `Name` and carries no path.
    """
    terms: list = []
    for inner in ast.walk(node):
        if (
            isinstance(inner, ast.Call)
            and isinstance(inner.func, ast.Attribute)
            and inner.func.attr == "exists"
            and not inner.args
        ):
            terms.append(_as_path(inner.func.value, consts))
        elif isinstance(inner, ast.Name) and inner.id in predicates:
            terms.extend(predicates[inner.id])
    return terms


def _module_level(tree: ast.Module) -> tuple[dict[str, str], dict[str, list]]:
    """Module-level names: those that resolve to a path, and those that are an existence predicate."""
    consts: dict[str, str] = {}
    predicates: dict[str, list] = {}
    for stmt in tree.body:
        if not (isinstance(stmt, ast.Assign) and len(stmt.targets) == 1):
            continue
        target = stmt.targets[0]
        if not isinstance(target, ast.Name):
            continue
        resolved = _as_path(stmt.value, consts)
        if resolved is not None:
            consts[target.id] = resolved
        elif _is_repo_root(stmt.value):
            consts[target.id] = ""
        terms = _existence_terms(stmt.value, consts, predicates)
        if terms:
            predicates[target.id] = terms
    return consts, predicates


def skip_guards(source: str, filename: str) -> list[tuple[str, list]]:
    """Every `skipif(...)` call site in one module, with the path terms its condition tests."""
    tree = ast.parse(source, filename=filename)
    consts, predicates = _module_level(tree)
    guards = []
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call) and node.args):
            continue
        func = node.func
        name = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", None)
        if name != "skipif":
            continue
        terms = _existence_terms(node.args[0], consts, predicates)
        if terms:
            guards.append((f"{filename}:{node.lineno}", terms))
    return guards


def inert_skip_guards(source: str, filename: str, tracked: frozenset[str]) -> list[tuple[str, list]]:
    """The guards that can never fire: every path their condition tests is tracked by git."""
    return [
        (where, terms)
        for where, terms in skip_guards(source, filename)
        if all(path is not None and path in tracked for path in terms)
    ]


def test_no_skip_guard_is_keyed_on_a_path_that_git_TRACKS() -> None:
    """A tracked file is in every checkout, so a skip waiting on its absence is not a guard at all."""
    tracked = _tracked_paths()
    inert = []
    for path in sorted(Path("tests").glob("test_*.py")):
        inert.extend(inert_skip_guards(path.read_text(), str(path), tracked))
    assert inert == [], (
        f"these skipif guards test for a file git TRACKS, so they are present in every checkout and "
        f"the skip can never fire -- remove the condition and let the test run, or make the absence "
        f"of a tracked input a failure: {inert}"
    )


def test_the_tracked_path_check_has_a_corpus_and_resolves_the_two_guards_that_DO_fire() -> None:
    """Non-vacuity over the real suite. A resolver that answered None everywhere would pass above.

    The two pinned guards are the ones measured on 2026-10-03 to key on genuinely UNTRACKED paths,
    and they are named by path rather than by line so an edit above them does not move this test.
    """
    tracked = _tracked_paths()
    resolved = {}
    for path in sorted(Path("tests").glob("test_*.py")):
        for where, terms in skip_guards(path.read_text(), str(path)):
            resolved[where] = terms
    with_a_path = [t for t in resolved.values() if any(p is not None for p in t)]
    assert len(with_a_path) >= 20, f"the resolver has stopped resolving paths: {len(with_a_path)}"

    gencode = [t for w, t in resolved.items() if w.startswith("tests/test_gene_identity.py:")]
    assert ["data/reference/gencode_v50_chr21.gff3.gz"] in gencode, gencode
    assert "data/reference/gencode_v50_chr21.gff3.gz" not in tracked

    compound = [t for w, t in resolved.items() if w.startswith("tests/test_context_evidence.py:")]
    assert compound, "the compound guards are gone"
    for terms in compound:
        assert "data/results/budget_chr21.json" in terms and None in terms, terms
    assert (
        inert_skip_guards(
            Path("tests/test_context_evidence.py").read_text(), "tests/test_context_evidence.py", tracked
        )
        == []
    ), "a compound guard that can fire on its second term must not be refused"


def test_the_tracked_path_check_REFUSES_a_guard_on_a_tracked_path() -> None:
    """The plant, both ways: the same guard refuses on a tracked path and passes on an untracked one."""
    planted = (
        "from pathlib import Path\n"
        "\n"
        "import pytest\n"
        "\n"
        'NEEDED = Path("data/results/planted.json")\n'
        "\n"
        "\n"
        '@pytest.mark.skipif(not NEEDED.exists(), reason="not written yet")\n'
        "def test_reads_a_committed_result():\n"
        "    pass\n"
    )
    where = [w for w, _ in skip_guards(planted, "planted.py")]
    assert where == ["planted.py:8"], where

    refused = inert_skip_guards(planted, "planted.py", frozenset({"data/results/planted.json"}))
    assert refused == [("planted.py:8", ["data/results/planted.json"])], refused
    # the identical guard on a path git does not track is left alone
    assert inert_skip_guards(planted, "planted.py", frozenset({"data/results/something_else.json"})) == []


def test_the_tracked_path_check_passes_a_guard_whose_second_term_can_be_absent() -> None:
    """The two shapes that must survive: one untracked term, and one term it cannot resolve."""
    both_tracked = (
        "from pathlib import Path\n"
        "\n"
        "import pytest\n"
        "\n"
        "import somewhere\n"
        "\n"
        "ROOT = Path(__file__).resolve().parents[1]\n"
        "\n"
        "\n"
        "@pytest.mark.skipif(\n"
        '    not ((ROOT / "data" / "results" / "a.json").exists() and Path("data/cache/b").exists()),\n'
        '    reason="one of the two is not here",\n'
        ")\n"
        "def test_needs_both():\n"
        "    pass\n"
        "\n"
        "\n"
        "@pytest.mark.skipif(\n"
        '    not ((ROOT / "data" / "results" / "a.json").exists() and somewhere.TABLE.exists()),\n'
        '    reason="one of the two is not here",\n'
        ")\n"
        "def test_needs_both_one_of_them_opaque():\n"
        "    pass\n"
    )
    guards = dict(skip_guards(both_tracked, "planted.py"))
    assert guards["planted.py:10"] == ["data/results/a.json", "data/cache/b"], guards
    assert guards["planted.py:18"] == ["data/results/a.json", None], guards
    only_a = frozenset({"data/results/a.json"})
    assert inert_skip_guards(both_tracked, "planted.py", only_a) == []
    # and with BOTH terms tracked the same compound guard IS refused, so the rule is not toothless
    both = frozenset({"data/results/a.json", "data/cache/b"})
    assert [w for w, _ in inert_skip_guards(both_tracked, "planted.py", both)] == ["planted.py:10"]
