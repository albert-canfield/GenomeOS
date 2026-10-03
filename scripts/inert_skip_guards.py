# SPDX-License-Identifier: AGPL-3.0-or-later
"""Every skip guard in `tests/`, with the path it is keyed on resolved and classified.

WHY A RESOLVER AND NOT A GREP. A skip guard keyed on a path git TRACKS can never fire: the file is in
the commit, so it is present in this checkout, in CI, in a verdict worktree and in the store-free
pre-push leg. The `skipif` reads to every later reviewer as "this test is conditional" when it is not,
and the identical shape keyed on a path that is SOMETIMES absent hides a real failure as a skip. The
rule the project adopted on 2026-10-03: `needs_local_data` skips BY NAME where a git-ignored store is
absent and RUNS where present, and a skip keyed on a TRACKED path must RAISE instead.

A text scan cannot find these. Two `skipif` lines in the whole of `tests/` name `data/results` as a
literal string; the rest reach it through a module constant, a `ROOT / "data" / "results" / ...`
expression, a local variable inside the test body, or `tests/local_data.present`. So the path is
resolved by syntax, and where syntax cannot resolve it the guard is reported UNRESOLVED by name rather
than folded into either answer.

FOUR SHAPES OF GUARD, and the fourth is the one the earlier instrument could not see.

  1. `@pytest.mark.skipif(not P.exists(), ...)` -- a decorator, condition in `args[0]`.
  2. `@pytest.mark.needs_local_data("data/cache/...", how=...)` -- the project marker. It CANNOT be
     inert: `tests/local_data.check_is_machine_local` asks git and RAISES on a path git does not
     ignore, so the marker polices itself. Reported as LIVE_BY_RULE and never converted.
  3. `pytest.importorskip("module")` -- keyed on an import, not on a path. Out of the population.
  4. `pytest.skip(...)` inside a test body or a fixture, under an enclosing `if`. This is where the
     population actually lives: 78 call sites against 88 `skipif` lines, and the committed invariant
     in `tests/test_local_data_symlinked_store.py` reads only shape 1.

CLASSES. `INERT` -- every path term resolved and every one of them tracked, so the condition is
constant. `LIVE` -- at least one resolved term is NOT tracked, so the guard can fire and is doing real
work. `UNRESOLVED` -- a term this cannot resolve, so the guard's class is unknown and is said to be.
`ALWAYS_FIRES` -- a tracked path tested for PRESENCE, so the skip is taken every time and the test
never runs; worse than inert and reported apart from it. `NOT_PATH` -- the condition tests no path at
all (a missing `node` binary, a shallow clone with no history, an errno this machine does not produce),
earned by finding no existence call and no resolvable path in the condition rather than assumed.

Run it: `uv run python scripts/inert_skip_guards.py` for the summary, `--class INERT` for one class,
`--json` for the rows.
"""

from __future__ import annotations

import argparse
import ast
import fnmatch
import json
import subprocess
import sys
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from functools import cache
from pathlib import Path

#: The git-ignored machine-local stores. Here for the reader and for `--stores`; the tracked/ignored
#: decision is asked of `git ls-files` and never taken from this list, because `data/results` is
#: ignored by PATTERN with `!` re-includes and hundreds of its files are committed, so a prefix test
#: would get it wrong in both directions.
IGNORED_STORES = ("data/cache", "data/knowledge", "data/individuals", "data/reference")

#: The calls that ask whether a path is there. `present` and `missing` are `tests/local_data`'s.
EXISTENCE_CALLS = frozenset({"exists", "is_file", "is_dir", "present", "missing", "glob", "iglob"})

INERT = "INERT"
LIVE = "LIVE"
LIVE_BY_RULE = "LIVE_BY_RULE"
UNRESOLVED = "UNRESOLVED"
ALWAYS_FIRES = "ALWAYS_FIRES"
NOT_PATH = "NOT_PATH"


@dataclass(frozen=True)
class Guard:
    """One skip guard: where it is, what shape it has, the paths it tests, and its class."""

    where: str
    shape: str
    verdict: str
    terms: tuple[str | None, ...]
    reason: str

    def __str__(self) -> str:
        shown = ", ".join("?" if t is None else t for t in self.terms) or "-"
        return f"{self.verdict:12} {self.shape:14} {self.where}  [{shown}]  {self.reason}"


def tracked_paths(root: Path | None = None) -> frozenset[str]:
    """Every path git tracks, as git spells it: repository-relative with forward slashes."""
    out = subprocess.run(["git", "ls-files"], capture_output=True, text=True, check=True, cwd=root or None)
    return frozenset(out.stdout.splitlines())


def is_tracked(path: str, tracked: frozenset[str]) -> bool:
    """Whether git tracks `path`, or -- for a glob -- whether it tracks anything the glob matches.

    A glob that matches a tracked file always finds that file, so such a guard is as constant as one
    on an exact tracked path. A directory counts as tracked when git tracks anything inside it.
    """
    if HOLE in path:
        return False  # a path with an unresolved fragment is never CALLED tracked
    if any(c in path for c in "*?["):
        return any(fnmatch.fnmatch(t, path) for t in tracked)
    if path in tracked:
        return True
    prefix = path.rstrip("/") + "/"
    return any(t.startswith(prefix) for t in tracked)


#: The file whose AST is being read, as a stack, so `Path(__file__).resolve().parent.parent` can be
#: turned into the repository-relative directory it actually names. It is a stack and not an argument
#: because `as_path` recurses through eight node kinds and an imported module is read in the middle of
#: reading a test; the alternative was threading one unchanging value through every one of them.
_ORIGIN: list[str] = ["tests/_nothing.py"]


@contextmanager
def reading(filename: str):
    """Read `filename`'s AST with `__file__` expressions resolved against its own place in the tree."""
    _ORIGIN.append(filename)
    try:
        yield
    finally:
        _ORIGIN.pop()


def _ancestor(levels: int) -> str:
    """The directory `levels` above the file being read, repository-relative. "" is the root."""
    parts = Path(_ORIGIN[-1]).parts
    return "/".join(parts[: max(len(parts) - levels, 0)])


def levels_above_file(node: ast.AST) -> int | None:
    """How far above `__file__` an expression climbs, or None when it is not such an expression.

    `Path(__file__).resolve().parent.parent` is two, `Path(__file__).resolve().parents[1]` is two,
    and `Path(__file__).parent` is one. COUNTED rather than treated as "a repository root somewhere",
    because the two are not the same answer: from a file in `tests/` one `parent` is `tests` and two
    is the root, and reading the first as the root silently moves every path the guard names.
    """
    if isinstance(node, ast.Call):
        func = node.func
        name = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", None)
        if (
            name == "Path"
            and len(node.args) == 1
            and isinstance(node.args[0], ast.Name)
            and node.args[0].id == "__file__"
        ):
            return 0
        if name in {"resolve", "absolute"} and isinstance(func, ast.Attribute):
            return levels_above_file(func.value)
        return None
    if isinstance(node, ast.Attribute):
        if node.attr == "parent":
            inner = levels_above_file(node.value)
            return None if inner is None else inner + 1
        return None
    if isinstance(node, ast.Subscript):
        holder = node.value
        if isinstance(holder, ast.Attribute) and holder.attr == "parents":
            inner = levels_above_file(holder.value)
            index = node.slice.value if isinstance(node.slice, ast.Constant) else None
            if inner is None or not isinstance(index, int):
                return None
            return inner + index + 1
    return None


def _is_repo_root(node: ast.AST) -> bool:
    """`Path(__file__).resolve().parents[1]` and its kin: the repository root spelled from a file."""
    return any(isinstance(n, ast.Attribute) and n.attr == "parents" for n in ast.walk(node))


#: The sentinel an unresolvable fragment of a path leaves behind, so the LITERAL part of the path is
#: still readable. `REFERENCE / (name + ".gz")` resolves to `data/reference/<?>`: the leaf is unknown
#: and the term can never be called tracked, but the prefix is enough to see it lies in a git-ignored
#: store, which is what makes the guard LIVE. Without this the term is a bare None and the guard's
#: class is unknown for no reason -- the prefix was there to be read.
HOLE = "<?>"


def has_hole(path: str | None) -> bool:
    """Whether a resolved path still carries an unresolved fragment."""
    return path is None or HOLE in path


def under_ignored_store(path: str | None) -> bool:
    """Whether the LITERAL prefix of `path` lies inside a git-ignored store.

    A gitignore entry for a directory ignores its whole subtree, so `data/reference/<?>` is
    machine-local whatever the leaf turns out to be. This is the one place a store prefix is used
    instead of asking git, and it is sound for exactly that reason: the question is about the prefix,
    which is literal, and not about the file, which is not.
    """
    if path is None:
        return False
    head = path.split(HOLE)[0]
    return any(head == s or head.startswith(s + "/") for s in IGNORED_STORES)


@cache
def _module_consts(dotted: str, depth: int = 0) -> dict[str, str]:
    """The path-valued module-level names of an imported module, read from its source.

    `mo.JASPAR_PATH`, `measured.CRISPRI_KNOWLEDGE`, `dp.RESULTS_DIR` and `dp.RESULT` are all real
    guard terms, and every one of them is a constant in another module. Reading the module's AST is
    how the finder resolves the path instead of reporting the guard unresolved -- and it is read, not
    imported, so finding out what a test is keyed on never runs that module's code.
    """
    if depth > 2:
        return {}
    source = Path(dotted.replace(".", "/") + ".py")
    if not source.is_file():
        source = Path(dotted.replace(".", "/")) / "__init__.py"
    if not source.is_file():
        return {}
    try:
        tree = ast.parse(source.read_text(), filename=str(source))
    except (OSError, SyntaxError):
        return {}
    consts: dict[str, str] = {}
    _ORIGIN.append(str(source))
    for alias, target in _imported_modules(tree).items():
        consts.update(
            {f"{alias}.{k}": v for k, v in _module_consts(target, depth + 1).items() if "." not in k}
        )
    for stmt in tree.body:
        if isinstance(stmt, ast.Assign) and len(stmt.targets) == 1:
            target_node = stmt.targets[0]
            if not isinstance(target_node, ast.Name):
                continue
            consts.update(_dict_entries(target_node.id, stmt.value, consts))
            resolved = as_path(stmt.value, consts)
            if resolved is not None:
                consts[target_node.id] = resolved
            elif _is_repo_root(stmt.value):
                consts[target_node.id] = ""
    _ORIGIN.pop()
    return consts


def _dict_entries(name: str, value: ast.AST, consts: dict[str, str]) -> dict[str, str]:
    """`NAME[key]` -> path, for a dict literal of literal keys and path values.

    `REFERENCES = {"celegans": ROOT / "data" / ...}` and then `Path(REFERENCES["celegans"])` is a real
    guard in `tests/test_body.py`, and without this the guard carries no path at all.
    """
    out: dict[str, str] = {}
    if not isinstance(value, ast.Dict):
        return out
    for key, item in zip(value.keys, value.values, strict=True):
        if isinstance(key, ast.Constant) and isinstance(key.value, str):
            resolved = as_path(item, consts)
            if resolved is not None:
                out[f"{name}[{key.value}]"] = resolved
    return out


def _imported_names(tree: ast.Module) -> dict[str, tuple[str, str]]:
    """Local name -> (module, attribute), for `from genomeos.x import CONSTANT`.

    A constant imported by name is the commonest unresolved term in this suite -- `PACKAGED`,
    `CELLS_FILE`, `KNOWLEDGE`, `REFERENCES` -- and each one is a path written out in the module it
    comes from. Resolving it there is the difference between a named class and a shrug.
    """
    out: dict[str, tuple[str, str]] = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module and not node.level:
            for alias in node.names:
                out[alias.asname or alias.name] = (node.module, alias.name)
    return out


def _imported_modules(tree: ast.Module) -> dict[str, str]:
    """Local name -> dotted module, for every import anywhere in the file, including inside a test.

    `from genomeos.genome import motifs as mo` sits inside the body of the test it guards, so the
    imports are taken from the whole tree and not from its top level.
    """
    out: dict[str, str] = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                out[alias.asname or alias.name.split(".")[0]] = (
                    alias.name if alias.asname else alias.name.split(".")[0]
                )
        elif isinstance(node, ast.ImportFrom) and node.module and not node.level:
            for alias in node.names:
                out[alias.asname or alias.name] = f"{node.module}.{alias.name}"
    return out


def as_path(node: ast.AST, consts: dict[str, str]) -> str | None:
    """A repository-relative path for an expression, or None when syntax cannot resolve it at all.

    `HOLE` where one fragment is unknown but the rest is literal; None where nothing is. None is the
    honest answer for `default_gencode(...)`, and it is what keeps such a guard out of INERT.
    """
    levels = levels_above_file(node)
    if levels is not None:
        return _ancestor(levels)
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    if isinstance(node, ast.Name):
        return consts.get(node.id)
    if isinstance(node, ast.Attribute):
        base = getattr(node.value, "id", None)
        return consts.get(f"{base}.{node.attr}") if base else None
    if isinstance(node, ast.Subscript) and isinstance(node.slice, ast.Constant):
        # `REFERENCES["celegans"]`: a literal key into a mapping of literal paths resolves as a name.
        holder = as_path(node.value, consts)
        key = node.slice.value
        if holder is None and isinstance(key, str):
            name = (
                node.value.id
                if isinstance(node.value, ast.Name)
                else f"{getattr(node.value.value, 'id', '?')}.{node.value.attr}"
                if isinstance(node.value, ast.Attribute)
                else None
            )
            if name is not None:
                return consts.get(f"{name}[{key}]")
        return None
    if isinstance(node, ast.JoinedStr):  # f-string: literal parts kept, the rest a hole
        parts = []
        for value in node.values:
            if isinstance(value, ast.Constant) and isinstance(value.value, str):
                parts.append(value.value)
            elif isinstance(value, ast.FormattedValue):
                parts.append(as_path(value.value, consts) or HOLE)
            else:
                parts.append(HOLE)
        return "".join(parts)
    if isinstance(node, ast.Call):
        func = node.func
        name = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", None)
        if name in {"Path", "str"} and len(node.args) == 1 and not node.keywords:
            return as_path(node.args[0], consts)
        if name == "format":
            return None
        return None
    if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Div, ast.Add)):
        left = as_path(node.left, consts)
        if left is None and _is_repo_root(node.left):
            left = ""  # ROOT / "data" / ... is the same path as "data/..." asked of git
        right = as_path(node.right, consts)
        if left is None and right is None:
            return None
        joiner = "/" if isinstance(node.op, ast.Div) else ""
        return f"{left if left is not None else HOLE}{joiner}{right if right is not None else HOLE}".lstrip(
            "/"
        )
    return None


def _existence_calls(node: ast.AST, into_comps: bool = False) -> list[ast.Call]:
    """Every call inside `node` that asks whether something is on disk.

    `into_comps` is False by default because a comprehension is read by SUBSTITUTION instead -- its
    loop variable is bound to each literal element, which turns one hole into the real paths.
    """
    found = []
    for inner in _walk(node, into_comps=into_comps):
        if not isinstance(inner, ast.Call):
            continue
        func = inner.func
        name = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", None)
        if name in EXISTENCE_CALLS:
            found.append(inner)
    return found


COMPREHENSIONS = (ast.ListComp, ast.SetComp, ast.GeneratorExp, ast.DictComp)


def _walk(node: ast.AST, into_comps: bool = True):
    """`ast.walk`, optionally not descending into comprehensions."""
    todo = [node]
    while todo:
        current = todo.pop()
        yield current
        for child in ast.iter_child_nodes(current):
            if into_comps or not isinstance(child, COMPREHENSIONS):
                todo.append(child)


def _comprehensions(node: ast.AST) -> list[ast.AST]:
    """The outermost comprehensions in `node`, counting `node` itself when it is one.

    `absent = [s for s in PLANT_OVERSHOOT if not (RESULTS / f"{s}.json").exists()]` binds the
    comprehension DIRECTLY to the name, so the expression being read IS the comprehension; missing
    that case left the loop variable unbound and reported `data/results/<?>.json` for a guard whose
    two paths are both written out in the tuple three lines above it.
    """
    if isinstance(node, COMPREHENSIONS):
        return [node]
    return [n for n in ast.iter_child_nodes(node) if isinstance(n, COMPREHENSIONS)] + [
        c for n in ast.iter_child_nodes(node) if not isinstance(n, COMPREHENSIONS) for c in _comprehensions(n)
    ]


def literal_sequence(
    node: ast.AST, consts: dict[str, str], sequences: dict[str, tuple[str, ...]]
) -> tuple[str, ...] | None:
    """A tuple/list/set of string constants, or a name bound to one. None when it is not literal."""
    if isinstance(node, (ast.Tuple, ast.List, ast.Set)):
        out = []
        for elt in node.elts:
            resolved = as_path(elt, consts)  # a bare name, a literal, or a `ROOT / ...` expression
            if resolved is None:
                return None
            out.append(resolved)
        return tuple(out)
    if isinstance(node, ast.Name):
        return sequences.get(node.id)
    if isinstance(node, ast.Attribute):
        base = getattr(node.value, "id", None)
        return sequences.get(f"{base}.{node.attr}") if base else None
    if isinstance(node, ast.Call):
        func = node.func
        name = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", None)
        if name == "sorted" and len(node.args) == 1:
            return literal_sequence(node.args[0], consts, sequences)
    return None


def _direct_terms(
    node: ast.AST,
    consts: dict[str, str],
    predicates: dict[str, list[str | None]],
    into_comps: bool = False,
) -> list[str | None]:
    """The path terms of one expression, not expanding comprehensions."""
    terms: list[str | None] = []
    for call in _existence_calls(node, into_comps=into_comps):
        func = call.func
        if isinstance(func, ast.Attribute) and func.attr in {"exists", "is_file", "is_dir", "glob"}:
            if func.attr == "glob" and call.args:
                base = as_path(func.value, consts)
                pat = as_path(call.args[0], consts)
                terms.append(None if base is None or pat is None else f"{base}/{pat}".lstrip("/"))
            else:
                terms.append(as_path(func.value, consts))
        else:  # present(...), missing(...), glob.glob(...): the path is an argument
            for arg in call.args:
                elts = arg.elts if isinstance(arg, (ast.Tuple, ast.List)) else [arg]
                terms.extend(as_path(e, consts) for e in elts)
    for inner in _walk(node, into_comps=into_comps):
        if isinstance(inner, ast.Name) and inner.id in predicates:
            terms.extend(predicates[inner.id])
    return terms


def helper_terms(
    node: ast.AST,
    helpers: dict[str, ast.AST],
    consts: dict[str, str],
    sequences: dict[str, tuple[str, ...]],
) -> list[str | None]:
    """The path terms of any module-local helper the condition CALLS, read from the helper's body.

    `if _benchmark_present() is None: pytest.skip(...)` carries no existence call and no path, so a
    reading that stops at the condition puts it in NOT_PATH -- "the condition tests no path" -- when
    the whole question is inside the three lines of `_benchmark_present`. That is the one way the
    NOT_PATH class could hide an inert guard, and the fix is to read the helper. Found on 2026-10-03
    while checking that a peer's claim about a cached-table skip was not being undone.
    """
    terms: list[str | None] = []
    for inner in ast.walk(node):
        if not isinstance(inner, ast.Call):
            continue
        name = getattr(inner.func, "id", None)
        if name is None or name not in helpers:
            continue
        body = helpers[name]
        inner_consts, inner_predicates = dict(consts), {}
        inner_sequences = dict(sequences)
        _bindings(body.body, inner_consts, inner_predicates, inner_sequences, nested=True)
        for statement in _walk(ast.Module(body=body.body, type_ignores=[])):
            if isinstance(statement, ast.If):
                terms.extend(existence_terms(statement.test, inner_consts, inner_predicates, inner_sequences))
            elif isinstance(statement, ast.Return) and statement.value is not None:
                terms.extend(
                    existence_terms(statement.value, inner_consts, inner_predicates, inner_sequences)
                )
    return terms


def module_helpers(tree: ast.Module) -> dict[str, ast.AST]:
    """The module's own non-test functions, by name, so a condition that calls one can be read."""
    return {
        node.name: node
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and not node.name.startswith("test")
    }


def existence_terms(
    node: ast.AST,
    consts: dict[str, str],
    predicates: dict[str, list[str | None]],
    sequences: dict[str, tuple[str, ...]] | None = None,
) -> list[str | None]:
    """Every path whose presence `node` tests, each resolved or None. Follows names and expands
    comprehensions.

    A predicate name is `HAS_CHR21 = (ROOT / "data" / "results" / "x.json").exists()` followed by
    `skipif(not HAS_CHR21, ...)`; without following the name the condition is a bare `Name` carrying
    no path. A comprehension is `absent = [s for s in PLANT_OVERSHOOT if not (RESULTS /
    f"{s}.json").exists()]`, and expanding it over the literal tuple is the difference between four
    named tracked results and one unresolved hole.
    """
    sequences = sequences or {}
    terms: list[str | None] = []
    for comp in _comprehensions(node):
        generator = comp.generators[0]
        elements = literal_sequence(generator.iter, consts, sequences)
        target = generator.target
        if elements is None or not isinstance(target, ast.Name):
            terms.extend(_direct_terms(comp, consts, predicates, into_comps=True))
            continue
        for value in elements:
            terms.extend(_direct_terms(comp, {**consts, target.id: value}, predicates, True))
    if not isinstance(node, COMPREHENSIONS):
        terms.extend(_direct_terms(node, consts, predicates))
    return terms


def _bindings(
    body: list[ast.stmt],
    consts: dict[str, str],
    predicates: dict[str, list[str | None]],
    sequences: dict[str, tuple[str, ...]],
    nested: bool = False,
) -> None:
    """Fold `name = ...` and `for name in ...` into `consts`, `predicates` and `sequences`, in place.

    `nested` reads the whole body including loops and `with` blocks, which is where a test binds the
    path it then guards: `for name in measured.CRISPRI_FILES: path = measured.CRISPRI_KNOWLEDGE /
    name`. A loop variable is bound to `HOLE` and NOT to one of its values, because which iteration
    the skip is taken on is not a thing syntax knows; the literal prefix is still read.
    """
    statements = list(_walk(ast.Module(body=body, type_ignores=[]))) if nested else list(body)
    for stmt in statements:
        if isinstance(stmt, ast.For) and isinstance(stmt.target, ast.Name):
            consts[stmt.target.id] = HOLE
            continue
        if not (isinstance(stmt, ast.Assign) and len(stmt.targets) == 1):
            continue
        target = stmt.targets[0]
        if not isinstance(target, ast.Name):
            continue
        consts.update(_dict_entries(target.id, stmt.value, consts))
        as_sequence = literal_sequence(stmt.value, consts, sequences)
        if as_sequence is not None:
            sequences[target.id] = as_sequence
        resolved = as_path(stmt.value, consts)
        if resolved is not None:
            consts[target.id] = resolved
        elif _is_repo_root(stmt.value):
            consts[target.id] = ""
        terms = existence_terms(stmt.value, consts, predicates, sequences)
        if terms:
            predicates[target.id] = terms


def module_level(
    tree: ast.Module,
) -> tuple[dict[str, str], dict[str, list[str | None]], dict[str, tuple[str, ...]]]:
    """The module's own names: paths, presence predicates and literal string sequences.

    The constants and sequences of every module it imports come in first, under the local alias, so
    `mo.JASPAR_PATH` and `measured.CRISPRI_FILES` resolve instead of being reported unresolved.
    """
    consts: dict[str, str] = {}
    predicates: dict[str, list[str | None]] = {}
    sequences: dict[str, tuple[str, ...]] = {}
    for alias, dotted in _imported_modules(tree).items():
        for key, value in _module_consts(dotted).items():
            if "." not in key:
                consts[f"{alias}.{key}"] = value
        for key, value in _module_sequences(dotted).items():
            sequences[f"{alias}.{key}"] = value
    for local, (dotted, attribute) in _imported_names(tree).items():
        from_module = _module_consts(dotted)
        if attribute in from_module:
            consts[local] = from_module[attribute]
        for key, value in from_module.items():
            if key.startswith(f"{attribute}["):
                consts[local + key[len(attribute) :]] = value
        from_sequences = _module_sequences(dotted)
        if attribute in from_sequences:
            sequences[local] = from_sequences[attribute]
    _bindings(tree.body, consts, predicates, sequences)
    for stmt in tree.body:  # a path bound inside `if TYPE_CHECKING:` or a try/except is still a path
        if isinstance(stmt, (ast.If, ast.Try)):
            _bindings(stmt.body, consts, predicates, sequences)
    return consts, predicates, sequences


@cache
def _module_sequences(dotted: str) -> dict[str, tuple[str, ...]]:
    """The literal string sequences a module defines at its top level, read from its source."""
    source = Path(dotted.replace(".", "/") + ".py")
    if not source.is_file():
        return {}
    try:
        tree = ast.parse(source.read_text(), filename=str(source))
    except (OSError, SyntaxError):
        return {}
    out: dict[str, tuple[str, ...]] = {}
    for stmt in tree.body:
        if isinstance(stmt, ast.Assign) and len(stmt.targets) == 1:
            target = stmt.targets[0]
            if isinstance(target, ast.Name):
                found = literal_sequence(stmt.value, {}, out)
                if found is not None:
                    out[target.id] = found
    return out


def _enclosing_function(tree: ast.Module, call: ast.Call) -> ast.AST | None:
    """The innermost function whose body holds `call`."""
    best, chosen = None, None
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.lineno <= call.lineno <= (
            node.end_lineno or node.lineno
        ):
            span = (node.end_lineno or node.lineno) - node.lineno
            if best is None or span <= best:
                best, chosen = span, node
    return chosen


def _enclosing_if(tree: ast.Module, call: ast.Call) -> tuple[ast.If | None, bool]:
    """The innermost `if` whose branch holds `call`, and whether that branch is the `else`.

    The guard of a `pytest.skip()` in a body is that `if`'s condition: the call itself carries only a
    reason string, so reading the call alone classifies every one of them as carrying no path. This
    is the shape the committed `skipif` invariant cannot see, and it is where the population lives.
    """
    best, chosen = None, (None, False)
    for node in ast.walk(tree):
        if not isinstance(node, ast.If):
            continue
        for in_else, branch in ((False, node.body), (True, node.orelse)):
            for stmt in branch:
                if stmt.lineno <= call.lineno <= (stmt.end_lineno or stmt.lineno):
                    span = (stmt.end_lineno or stmt.lineno) - stmt.lineno
                    if best is None or span <= best:
                        best, chosen = span, (node, in_else)
    return chosen


def parameter_values(tree: ast.Module, function: ast.AST) -> dict[str, tuple[str, ...]]:
    """For a helper, the literal strings every call site in the module passes for each parameter.

    `def _result(name): p = RESULTS / f"{name}.json"; if not p.exists(): pytest.skip(...)` is the
    commonest guard shape in this suite, and read on its own the path is `data/results/<?>.json`.
    Where EVERY call site passes a string literal the parameter is known to take exactly those
    values, so the guard is keyed on exactly those paths -- and where one call site passes anything
    else, nothing is bound and the term stays a hole. The test functions are excluded because their
    arguments come from `parametrize` and not from a call.
    """
    name = getattr(function, "name", "")
    parameters = [a.arg for a in getattr(function, "args", ast.arguments()).args]
    if not parameters or name.startswith("test"):
        return {}
    seen: dict[str, list[str]] = {p: [] for p in parameters}
    calls = 0
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call) and getattr(node.func, "id", None) == name):
            continue
        calls += 1
        for parameter, argument in zip(parameters, node.args, strict=False):
            if isinstance(argument, ast.Constant) and isinstance(argument.value, str):
                seen[parameter].append(argument.value)
    if not calls:
        return {}
    return {p: tuple(dict.fromkeys(v)) for p, v in seen.items() if len(v) == calls}


def _tests_for_presence(condition: ast.AST, in_else: bool) -> bool:
    """Whether the branch holding the skip is the one taken when the path IS there.

    `if P.exists(): pytest.skip(...)` on a tracked path skips EVERY run, so the test never runs at
    all: a different and worse defect than a skip that never fires, and not reported as the same one.
    Read as: an existence call under a `not`, or a `missing(...)`/`absent` shape, means "when gone".
    """
    negated = any(isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.Not) for n in ast.walk(condition))
    inverted = any(
        isinstance(n, ast.Call)
        and (n.func.attr if isinstance(n.func, ast.Attribute) else getattr(n.func, "id", None)) == "missing"
        for n in ast.walk(condition)
    )
    names = {n.id.lower() for n in ast.walk(condition) if isinstance(n, ast.Name)}
    worded = any("absent" in n or "missing" in n for n in names)
    if not _existence_calls(condition, into_comps=True):
        # The direction is NOT readable here. `if _present() is None:` has no `not` and no existence
        # call, and the inversion lives inside the helper's `return None`, so reading the surface
        # syntax called a perfectly ordinary absent-guard ALWAYS_FIRES. A constant condition is still
        # constant -- that is what INERT says -- but which constant it is, is not claimed on syntax
        # that does not carry it. Measured against a plant on 2026-10-03.
        return False
    return (negated or inverted or worded) == in_else


def _verdict(terms: list[str | None], tracked: frozenset[str], for_presence: bool) -> tuple[str, str]:
    """The class of a guard from its path terms, and the sentence that says why."""
    if not terms:
        return NOT_PATH, "the condition tests no path"
    # A term with a hole is LIVE only on its LITERAL prefix lying in an ignored store. "not tracked"
    # is NOT enough: `<?>/chr21.json.gz` is untracked only because nothing can match it, and reading
    # that as "the guard can fire" is the very fold -- unresolved into live -- the classes exist to
    # prevent. It was doing exactly that until it was measured.
    live = [
        t
        for t in terms
        if t is not None and (under_ignored_store(t) or (not has_hole(t) and not is_tracked(t, tracked)))
    ]
    if live:
        where = live[0]
        why = "lies in a git-ignored store" if under_ignored_store(where) else "is not tracked by git"
        return LIVE, f"{where} {why}, so the guard can fire"
    unknown = [t for t in terms if has_hole(t)]
    if unknown:
        shown = ", ".join(t if t is not None else "an expression syntax cannot resolve" for t in unknown)
        return UNRESOLVED, f"{len(unknown)} of {len(terms)} terms do not resolve: {shown}"
    if for_presence:
        return ALWAYS_FIRES, "every path is tracked and the skip is taken when it IS there"
    return INERT, "git tracks every path the condition tests, so it is never absent"


def guards_in_module(source: str, filename: str, tracked: frozenset[str]) -> list[Guard]:
    """Every skip guard in one test module, classified."""
    tree = ast.parse(source, filename=filename)
    with reading(filename):
        return _guards(tree, filename, tracked)


def _guards(tree: ast.Module, filename: str, tracked: frozenset[str]) -> list[Guard]:
    """The body of `guards_in_module`, with `filename` on the `_ORIGIN` stack."""
    consts, predicates, sequences = module_level(tree)
    helpers = module_helpers(tree)
    out: list[Guard] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        name = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", None)
        where = f"{filename}:{node.lineno}"
        on_mark = isinstance(func, ast.Attribute) and getattr(func.value, "attr", None) == "mark"
        if name == "needs_local_data":
            out.append(
                Guard(
                    where,
                    "needs_local_data",
                    LIVE_BY_RULE,
                    tuple(as_path(a, consts) for a in node.args),
                    "local_data.check_is_machine_local asks git and RAISES on a tracked path",
                )
            )
        elif name == "importorskip":
            out.append(
                Guard(
                    where,
                    "importorskip",
                    NOT_PATH,
                    tuple(as_path(a, consts) for a in node.args[:1]),
                    "keyed on an import, not on a path",
                )
            )
        elif name == "skip" and on_mark:
            out.append(Guard(where, "mark.skip", NOT_PATH, (), "an unconditional skip, keyed on nothing"))
        elif name == "skipif" and node.args:
            terms = existence_terms(node.args[0], consts, predicates, sequences)
            terms += helper_terms(node.args[0], helpers, consts, sequences)
            verdict, why = _verdict(terms, tracked, for_presence=False)
            out.append(Guard(where, "skipif", verdict, tuple(terms), why))
        elif name == "skip" and isinstance(func, ast.Attribute):
            branch, in_else = _enclosing_if(tree, node)
            if branch is None:
                out.append(Guard(where, "skip()", NOT_PATH, (), "unconditional: it is under no `if` at all"))
                continue
            scope = _enclosing_function(tree, node)
            local_consts, local_predicates = dict(consts), dict(predicates)
            local_sequences = dict(sequences)
            variants: list[dict[str, str]] = [{}]
            if scope is not None:
                _bindings(scope.body, local_consts, local_predicates, local_sequences, nested=True)
                for parameter, values in parameter_values(tree, scope).items():
                    variants = [{**v, parameter: value} for v in variants for value in values]
            terms: list[str | None] = []
            for variant in variants:
                if not variant:
                    bound = (local_consts, local_predicates, local_sequences)
                else:  # the parameter is bound FIRST, then the body's own `p = ...` is folded again
                    bound = ({**local_consts, **variant}, {}, dict(local_sequences))
                    if scope is not None:
                        _bindings(scope.body, *bound, nested=True)
                terms.extend(existence_terms(branch.test, *bound))
            terms.extend(helper_terms(branch.test, helpers, local_consts, local_sequences))
            verdict, why = _verdict(terms, tracked, for_presence=_tests_for_presence(branch.test, in_else))
            out.append(Guard(where, "skip()", verdict, tuple(dict.fromkeys(terms)), why))
    return out


def all_guards(
    tests: Path = Path("tests"), root: Path | None = None, tracked_only: bool = False
) -> list[Guard]:
    """Every skip guard in every test module, classified. The whole population, in file order.

    `tracked_only` leaves out test files git does not track. The CLI shows everything, because a lane
    wants its own new file checked before it commits it; the suite-wide INVARIANT asks for tracked
    only, and the reason is this shared checkout. Several sessions edit this tree at once, so an
    untracked `tests/test_*.py` is a peer's work in progress and not part of the suite. An invariant
    that read it would turn one lane's half-written file into every other lane's red push -- and it
    did, within the hour: `tests/test_manifest_rebuild_reads.py` appeared untracked with a guard on
    `data/results` while this was being written. Scoping to the commit is not a weaker claim; it is
    the claim correctly aimed, and the guard is still caught the moment the file is committed.
    """
    tracked = tracked_paths(root)
    out: list[Guard] = []
    for path in sorted(tests.glob("*.py")):
        if tracked_only:
            try:
                relative = path.resolve().relative_to((root or Path()).resolve()).as_posix()
            except ValueError:
                relative = path.as_posix()
            if relative not in tracked:
                continue
        out.extend(guards_in_module(path.read_text(), str(path), tracked))
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--class", dest="only", default=None, help="show one class only")
    parser.add_argument("--json", action="store_true", help="the rows, as JSON")
    parser.add_argument("--tests", default="tests", help="the directory to read")
    args = parser.parse_args(argv)

    guards = all_guards(Path(args.tests))
    if args.only:
        guards = [g for g in guards if g.verdict == args.only]
    if args.json:
        print(json.dumps([asdict(g) for g in guards], indent=2))
        return 0
    for guard in guards:
        print(guard)
    counts: dict[str, int] = {}
    for guard in all_guards(Path(args.tests)):
        counts[guard.verdict] = counts.get(guard.verdict, 0) + 1
    print("\n  " + "  ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    return 0


if __name__ == "__main__":
    sys.exit(main())
