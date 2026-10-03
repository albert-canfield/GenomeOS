# SPDX-License-Identifier: AGPL-3.0-or-later
"""What the project declares, read from pyproject's lock rather than from a list anybody maintains.

Ten tests passed, for as long as they existed, only because somebody had once installed pandas and
openpyxl into this checkout by hand (2026-10-02). Seven of them passed *transitively*: CI's `test` job
runs `uv sync --extra compose`, and eight locked packages pull pandas in behind that one flag, so the
seven were green for a reason unrelated to this project declaring pandas. The other three needed
openpyxl, which arrives only with biolearn in the `clocks` extra, which no CI job installs: they were
never in CI at all. This module exists so the next one is a failure and not a surprise.

Three sets, all derived from `uv.lock`:

    bare_distributions()   what `uv sync --frozen` installs: the project plus its dev group, resolved
                           through the lock's own dependency edges. This is the environment the
                           `test-bare` CI job builds, and the only one that shows what the project
                           declares rather than what somebody once installed.
    extra_distributions()  per extra, the distributions it adds BEYOND the bare set. numpy is in the
                           `analysis` extra and in the dev group, so `analysis` adds nothing and numpy
                           is not extra-only: the subtraction is what makes the set honest.
    bare_modules()         the top-level importable names the bare set provides.

The one thing the lock does not hold is which module name a distribution installs: nothing in a lock
file says that scikit-learn is imported as `sklearn`. Two sources are used, in order:

    1. `importlib.metadata.packages_distributions()`, which reads the metadata of what is installed in
       the running interpreter. For the bare set this is authoritative wherever the suite can run at
       all, because the bare set is by definition installed.
    2. the normalised distribution name (`process-bigraph` -> `process_bigraph`), for anything not
       installed here -- which, in a bare environment, is every extra-only distribution.

So `extras_for_module` is exact in a full checkout and name-shaped in a bare one. That blind spot is
named in `scan`'s docstring, and `scan` is built so the blind spot cannot hide an import: a module that
no source accounts for is reported too, rather than assumed innocent.
"""

from __future__ import annotations

import ast
import errno
import sys
import tomllib
from collections.abc import Callable
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOCK = ROOT / "uv.lock"

#: the lock's name for this project, as `[project] name` in pyproject.toml spells it
PROJECT = "genomeos"

#: the dependency group `uv sync` installs when nobody asks for one (uv's default-groups default)
DEFAULT_GROUP = "dev"

#: directories whose Python files are this repository's own code, searched in this order for a module
FIRST_PARTY_ROOTS = (ROOT, ROOT / "scripts", ROOT / "tests")

#: the top-level names those directories provide, so an import of one is never a third-party import
FIRST_PARTY_TOP = ("genomeos", "scripts", "tests", "conftest")


def _norm(name: str) -> str:
    """A distribution name as a module name would spell it; PEP 503 normalisation, then PEP 8."""
    return name.lower().replace("-", "_").replace(".", "_")


@lru_cache(maxsize=1)
def _lock() -> dict:
    return tomllib.loads(LOCK.read_text())


def reset_caches() -> None:
    """Forget every set derived from the lock.

    A test points `LOCK` at another file -- the counterfactual lock that does not declare pandas -- and
    needs an answer about that one; without this, whichever reader came first would have cached the real
    lock for the rest of the session. Call it after changing `LOCK`, and again when putting it back.
    """
    for name, value in list(globals().items()):
        if hasattr(value, "cache_clear") and name != "reset_caches":
            value.cache_clear()


@lru_cache(maxsize=1)
def _by_name() -> dict[str, dict]:
    return {package["name"]: package for package in _lock()["package"]}


def _requirements(package: dict, extra: str | None) -> list[dict]:
    if extra is None:
        return package.get("dependencies", [])
    return package.get("optional-dependencies", {}).get(extra, [])


def _closure(roots: list[tuple[str, str | None]]) -> frozenset[str]:
    """Every distribution reachable from `roots`, following the lock's edges and nothing else.

    A root is a (distribution, extra) pair, because an edge in a lock file may ask for a dependency's
    extra (`{name = "x", extra = ["y"]}`), and that extra's requirements are then part of the closure
    while the same distribution's other extras are not. Markers are deliberately not evaluated: see
    `scan`.
    """
    seen: set[tuple[str, str | None]] = set()
    stack = list(roots)
    while stack:
        name, extra = stack.pop()
        if (name, extra) in seen:
            continue
        seen.add((name, extra))
        package = _by_name().get(name)
        if package is None:
            continue
        for requirement in _requirements(package, extra):
            stack.append((requirement["name"], None))
            for dependency_extra in requirement.get("extra") or ():
                stack.append((requirement["name"], dependency_extra))
    return frozenset(name for name, _ in seen)


@lru_cache(maxsize=1)
def bare_distributions() -> frozenset[str]:
    """What `uv sync --frozen` installs: this project, its dependencies, and its default dev group."""
    project = _by_name()[PROJECT]
    roots: list[tuple[str, str | None]] = [(PROJECT, None)]
    for requirement in project.get("dev-dependencies", {}).get(DEFAULT_GROUP, []):
        roots.append((requirement["name"], None))
        for extra in requirement.get("extra") or ():
            roots.append((requirement["name"], extra))
    return _closure(roots)


@lru_cache(maxsize=1)
def extra_distributions() -> dict[str, frozenset[str]]:
    """Per extra, the distributions it adds that the bare environment does not already hold."""
    bare = bare_distributions()
    project = _by_name()[PROJECT]
    return {extra: _closure([(PROJECT, extra)]) - bare for extra in project.get("optional-dependencies", {})}


@lru_cache(maxsize=1)
def known_extras() -> tuple[str, ...]:
    """The extras pyproject provides, as the lock records them."""
    return tuple(sorted(_by_name()[PROJECT].get("optional-dependencies", {})))


@lru_cache(maxsize=1)
def _installed_modules() -> dict[str, frozenset[str]]:
    """Distribution -> the top-level modules its installed metadata says it provides."""
    from importlib.metadata import packages_distributions

    found: dict[str, set[str]] = {}
    for module, distributions in packages_distributions().items():
        for distribution in distributions:
            found.setdefault(_norm(distribution), set()).add(module)
    return {name: frozenset(modules) for name, modules in found.items()}


def modules_of(distribution: str) -> frozenset[str]:
    """The module names a distribution provides: its metadata if installed, else its own name."""
    return _installed_modules().get(_norm(distribution)) or frozenset({_norm(distribution)})


@lru_cache(maxsize=1)
def bare_modules() -> frozenset[str]:
    """Every top-level module an import may reach in a bare environment."""
    modules: set[str] = set(FIRST_PARTY_TOP)
    for distribution in bare_distributions():
        modules |= modules_of(distribution)
    return frozenset(modules)


@lru_cache(maxsize=1)
def _module_extras() -> dict[str, tuple[str, ...]]:
    found: dict[str, set[str]] = {}
    for extra, distributions in extra_distributions().items():
        for distribution in distributions:
            for module in modules_of(distribution):
                found.setdefault(module, set()).add(extra)
    return {module: tuple(sorted(extras)) for module, extras in found.items()}


def extras_for_module(module: str) -> tuple[str, ...]:
    """The extras that would install a top-level module, or () if no extra reaches it."""
    return _module_extras().get(module, ())


def _importable(module: str) -> bool:
    """Whether THIS ENVIRONMENT provides a top-level module, asked of the finders and not of sys.modules.

    `importlib.util.find_spec`, which this used until 2026-10-03, answers from `sys.modules` BEFORE it
    asks the filesystem, and a pytest session is one process. So the answer moved with whatever an
    earlier test had put there, in both directions, and both were measured on this machine:

      * a `types.ModuleType("alphagenome")` -- exactly what tests/test_adapter_gene_axis.py and
        tests/test_model_version_pin.py build to stand in for the uninstalled client -- has
        `__spec__` None, so `find_spec` raised `ValueError: alphagenome.__spec__ is None`, the
        `except (ImportError, ValueError)` below read that as absent, and `missing_modules_for_extra`
        reported the INSTALLED `predict` extra missing. Five tests in tests/test_astrorun.py then
        skipped and the suite stayed green: measured 173 passed / 6 skipped against 178 passed /
        1 skipped. A skip is how a test disappears quietly, which is the failure this whole module
        was written against.
      * a stub carrying a real spec went the other way: `compose` read as PRESENT with
        process_bigraph not installed, so a test that must skip would run against the stand-in --
        the 2026-10-02 defect itself, a pass for a reason unrelated to the project declaring
        anything.

    Both stub sites use `monkeypatch.setitem` today, so they restore and the gate (which runs in
    `pytest_runtest_setup`, before a test's fixtures) never sees them. The dependence was therefore
    LATENT and not active -- the same word the 900 MB and 4 GiB `ru_maxrss` ceilings earned the same
    day, which were also safe only by the order the files happened to run in.

    The finders are asked in `sys.meta_path` order, which is what an actual `import` would consult
    after the `sys.modules` shortcut, so a module provided by an editable install's own finder still
    answers True. The residual blind spot, named rather than closed: a test that installs a
    MetaPathFinder of its own and leaves it there moves this answer too. Nothing in this suite does;
    `importlib.machinery.PathFinder` alone would be narrower but would miss editable installs.
    """
    for finder in sys.meta_path:
        try:
            spec = finder.find_spec(module, None)
        except (ImportError, AttributeError, TypeError, ValueError):
            continue
        if spec is not None:
            return True
    return False


def missing_modules_for_extra(extra: str, present: Callable[[str], bool] | None = None) -> tuple[str, ...]:
    """The modules of an extra's own requirements that this environment does not provide.

    Only the extra's direct requirements are checked, not its whole closure: those are what pyproject
    names, and an extra is here when what it names can be imported. Empty means the extra is installed.

    `present` is the reading, INJECTED so that a test can prove this answers about what it was handed
    rather than about the process it happens to run in. A real call passes nothing and gets
    `_importable`, which is the live environment: tests/conftest.py's skip gate calls this with the
    extra alone, and tests/test_extra_only_imports.py pins that it does. A guard keyed to what its
    caller supplies would be no guard, so the default is the only reading production ever gets.
    """
    reading = _importable if present is None else present
    if extra not in known_extras():
        raise LookupError(f"no extra named {extra!r} in {LOCK.name}; pyproject provides {known_extras()}")
    missing = []
    for requirement in _by_name()[PROJECT]["optional-dependencies"][extra]:
        for module in sorted(modules_of(requirement["name"])):
            if not reading(module):
                missing.append(module)
    return tuple(missing)


# --------------------------------------------------------------------------------------------------
# reading the tests


@dataclass(frozen=True)
class Import:
    """One top-level module name an import statement reaches, and where it is written."""

    module: str
    path: Path
    line: int


@dataclass(frozen=True)
class Violation:
    """A test module that reaches a package the bare environment does not hold, without saying so."""

    test: Path
    module: str
    through: Path
    line: int
    extras: tuple[str, ...]

    def __str__(self) -> str:
        where = f"{self.through.name}:{self.line}"
        if self.extras:
            how = "mark it @pytest.mark.requires_extra(" + " or ".join(repr(e) for e in self.extras) + ")"
        else:
            how = "and no extra installs it either: the project declares it nowhere"
        return f"{self.test.name}: imports {self.module} at {where} -- {how}"


def _exists(path: Path) -> bool:
    """`is_file`, but case-sensitively: this project is developed on a case-insensitive filesystem.

    Without the listing check, `genomeos/genome/Genome.py` answers True on macOS and the walk follows a
    file no name in the tree spells (the same trap cd263bc fixed in the cleanliness counter).

    A candidate is a module name turned into a path, and a module name is not always a legal one. On
    2026-10-02 this raised `OSError` 63, ENAMETOOLONG, on `scripts/push_own/…`: a long dotted-looking
    STRING inside tests/test_push_own.py was read as a dotted module and statted. The scanner then
    died, so no tree could earn a green verdict and every push gated on one was blocked — a scanner
    that cannot answer "no" for an impossible name takes the whole suite down with it. An unaskable
    question is answered False, not raised: a name the filesystem cannot hold is not a file in it.
    The errors are named rather than swallowed as a bare `except OSError`, so a permission or I/O
    fault on a path that COULD exist still surfaces instead of reading as a clean miss.
    """
    try:
        return path.is_file() and path.name in {entry.name for entry in path.parent.iterdir()}
    except OSError as e:
        if e.errno in (errno.ENAMETOOLONG, errno.ENOENT, errno.ENOTDIR, errno.EINVAL, errno.ELOOP):
            return False
        raise


@lru_cache(maxsize=4096)
def _first_party_file(module: str) -> Path | None:
    parts = module.split(".")
    for base in FIRST_PARTY_ROOTS:
        candidate = base.joinpath(*parts)
        for path in (candidate.with_suffix(".py"), candidate / "__init__.py"):
            if _exists(path):
                return path
    return None


_IMPORT_ERRORS = frozenset({"ImportError", "ModuleNotFoundError", "Exception", "BaseException"})


def _dotted(node: ast.AST | None) -> str:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        return f"{_dotted(node.value)}.{node.attr}"
    return ""


def _catches_import_error(handler: ast.ExceptHandler) -> bool:
    if handler.type is None:
        return True
    caught = handler.type.elts if isinstance(handler.type, ast.Tuple) else [handler.type]
    return any(_dotted(element) in _IMPORT_ERRORS for element in caught)


def _guarded(stack: list[ast.AST]) -> bool:
    """True if an import at this position cannot raise: a try that catches it, or TYPE_CHECKING.

    A handler body is not a guard -- `except ImportError: import other` runs `other` unprotected -- so
    only the try body counts, which is why the child is compared against the parent's `body`.
    """
    for parent, child in zip(stack, stack[1:], strict=False):
        if not any(child is statement for statement in getattr(parent, "body", ())):
            continue
        if isinstance(parent, ast.Try) and any(_catches_import_error(h) for h in parent.handlers):
            return True
        if isinstance(parent, ast.If) and _dotted(parent.test).endswith("TYPE_CHECKING"):
            return True
    return False


def _deferred(stack: list[ast.AST]) -> bool:
    """True if an import does not run when the file is imported, because a function holds it."""
    return any(isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda)) for node in stack)


def _parse(path: Path) -> ast.Module | None:
    try:
        return ast.parse(path.read_text(), filename=str(path))
    except (OSError, SyntaxError, UnicodeDecodeError):
        return None


@lru_cache(maxsize=4096)
def read_module(path: Path) -> tuple[tuple[Import, ...], frozenset[str], frozenset[str]]:
    """The import-time imports, the first-party modules, and the extras declared, in one file.

    Only imports that run when the file is imported are collected: not the ones written inside a
    function, and not the ones inside a try that catches ImportError or an `if TYPE_CHECKING`. That is
    the line pyproject already draws -- "everywhere else numpy is imported inside the function that
    needs it, which is the pattern to keep" -- and it is the line the defect crossed. A module-level
    import of a package nobody declared does not fail one test: the file does not COLLECT, which reads
    as "no such tests" rather than as a missing dependency.

    First-party edges are taken from import statements and from string literals, because the closure
    that hid openpyxl was a string: `tests/test_code_cleanliness_shared.py` names
    `"scripts.astroreg_register"` in a tuple and imports it with importlib, and
    `scripts/astroreg_register.py` is where `import openpyxl` is written. A bare name counts as
    first-party when a file of that name exists in one of `FIRST_PARTY_ROOTS`, which is how the tests
    that insert `scripts/` on `sys.path` reach `import capacity_gate`.

    Cached, and so returning only immutable containers: the 300-odd test modules of this suite share one
    closure of project files, and parsing each file once rather than once per test turns the whole-suite
    check from three quarters of a minute into a second.
    """
    tree = _parse(path)
    imports: list[Import] = []
    first_party: set[str] = set()
    extras: set[str] = set()
    if tree is None:
        return (), frozenset(), frozenset()

    stack: list[ast.AST] = []

    def visit(node: ast.AST) -> None:
        stack.append(node)
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            else:
                names = [node.module] if node.level == 0 and node.module else []
            for name in names:
                if name.split(".")[0] in FIRST_PARTY_TOP or _first_party_file(name) is not None:
                    first_party.add(name)
                elif not _guarded(stack) and not _deferred(stack):
                    imports.append(Import(name.split(".")[0], path, node.lineno))
        elif isinstance(node, ast.Constant) and isinstance(node.value, str):
            text = node.value.removesuffix(".py").replace("/", ".")
            if text.split(".")[0] in FIRST_PARTY_TOP and _first_party_file(text) is not None:
                first_party.add(text)
        elif isinstance(node, ast.Call) and _dotted(node.func).endswith("requires_extra"):
            # `update`, not `|=`: an augmented assignment would make `extras` local to `visit` and
            # every marked module would raise UnboundLocalError instead of being cleared.
            extras.update(
                argument.value
                for argument in node.args
                if isinstance(argument, ast.Constant) and isinstance(argument.value, str)
            )
        for child in ast.iter_child_nodes(node):
            visit(child)
        stack.pop()

    visit(tree)
    return tuple(imports), frozenset(first_party), frozenset(extras)


def scan(test_files: list[Path]) -> list[Violation]:
    """Every test module that can reach a package a bare environment does not hold, without saying so.

    For each test file the first-party closure is walked -- the project's own modules it imports, and
    theirs, by import statement or by module-name string -- and every unguarded third-party import in
    that closure is checked against `bare_modules()`. An import the bare set does not provide is a
    violation unless the test file carries `requires_extra` naming an extra that reaches it.

    What it cannot see, and each is a reason a clean run is not a proof:

      * an import written as a string the walk does not recognise (`__import__` of a computed name, a
        plugin named in a config file, an entry point);
      * a module name that neither the installed metadata nor normalisation gets right, in an
        environment where that distribution is not installed: such an import is still reported, because
        it is not in `bare_modules()`, but the message cannot say which extra installs it;
      * a package the lock reaches only under a marker this walk ignores (platform, python version), so
        a distribution installed on no runner still counts as bare. Ignoring markers is the safe
        direction for the extras sets and the unsafe one here;
      * whether a guarded import's fallback works, or whether a test asserts anything once it skips. A
        guard and a skip are both taken at their word.
    """
    bare = bare_modules()
    violations: list[Violation] = []
    for test in sorted(test_files):
        _, _, declared = read_module(test)
        seen: set[Path] = set()
        queue = [test]
        while queue:
            path = queue.pop()
            if path in seen:
                continue
            seen.add(path)
            imports, first_party, _ = read_module(path)
            for module in first_party:
                found = _first_party_file(module)
                if found is not None and found not in seen:
                    queue.append(found)
            for item in imports:
                if item.module in bare or item.module in sys.stdlib_module_names:
                    continue
                extras = extras_for_module(item.module)
                if declared.intersection(extras):
                    continue
                violations.append(Violation(test, item.module, item.path, item.line, extras))
    return violations
