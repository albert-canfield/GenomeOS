# SPDX-License-Identifier: AGPL-3.0-or-later
"""Every place the code writes into data/results/ other than through `save_result` (item 12 S6 follow-up).

    uv run python scripts/results_writer_census.py            # the table, one line per site
    uv run python scripts/results_writer_census.py --json     # the sites, for tooling

Item 12 S6 (docs/DATA.md, "The result registry") made `genomeos.results.save_result` refuse a new
result without a complete manifest. lane-s6's census found 13 writers that put results into
data/results/ without calling it, by a pattern (a name assigned a path into data/results, then
`.write_text` or `open(name, "w")` on it), so a lower bound. This census reads every tracked Python
file under genomeos/ and scripts/ with a static taint analysis, so the guard test
(tests/test_results_writers_guard.py) can hold the line afterwards:

- a *results path* is a string naming data/results, the registry constant `RESULTS_DIR`, a path
  joined from "data" and "results", a parameter named `results_dir`, a `--flag` whose default is a
  results path *in any form, a literal or a name* (read back as `args.flag`), or anything built from
  one of those by path operations
  (`/`, `Path(...)`, `str(...)`, `os.path.join`, f-strings, `.with_suffix`, `.glob`, a `for` over
  one, a container holding one, a local or imported function that returns one);
- a *sink* is a call that writes to its target: `.write_text`, `.write_bytes`, `open(..., "w"|"a"|
  "x"|"+")` and `Path.open`, `gzip/bz2/lzma/io/codecs/tarfile.open` in a write mode, `shutil.copy*`
  and `shutil.move` (destination), `os.replace` and `os.rename` (destination), `Path.rename` and
  `Path.replace` with one argument, `.touch`, `.symlink_to`, `.hardlink_to`, `np.save*`, and a
  DataFrame's `.to_csv/.to_json/.to_parquet/.to_pickle` or a figure's `.savefig`;
- a *writer helper* is a function (in the same file or imported from another analysed file) one of
  whose parameters reaches a sink; a call that passes a results path to that parameter is a site;
- a write mode that is not a constant, and a subprocess command that names a results path, are
  reported as `possible` for a reader to decide.

`save_result` itself (genomeos/results.py) is the sanctioned sink and is not reported. Removals
(`.unlink`, `os.remove`, `shutil.rmtree`) are reported apart: they change the registry without
writing a result. Flow-insensitive within a function, so it over-reports rather than under-reports;
every site is read by hand in the census.
"""

from __future__ import annotations

import argparse
import ast
import json
import re
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOTS = ("genomeos", "scripts")
#: The one module allowed to write into the registry: it is the contract.
SANCTIONED = ("genomeos/results.py",)

RESULTS_TEXT = re.compile(r"(?:^|/)data/results(?:/|$)")
DATA_TEXT = re.compile(r"(?:^|/)data/?$")
RESULTS_CONSTANTS = {"RESULTS_DIR"}
RESULTS_PARAMS = {"results_dir"}
WRITE_MODE = re.compile(r"[wax+]")
PATH_CALLS = {"Path", "PurePath", "PosixPath", "str", "fspath"}
OS_PATH_CALLS = {"join", "abspath", "realpath", "expanduser", "normpath"}
PATH_METHODS = {
    "joinpath",
    "with_suffix",
    "with_name",
    "with_stem",
    "resolve",
    "absolute",
    "expanduser",
    "glob",
    "rglob",
    "iterdir",
    "relative_to",
    "__truediv__",
}
FILE_OPENERS = {"gzip", "bz2", "lzma", "io", "codecs", "tarfile", "bgzf", "xopen"}
SUBPROCESS = {"run", "call", "check_call", "check_output", "Popen"}

R = "R"  # the results label; parameter labels are "P:<function qualname>:<parameter>"


@dataclass
class Site:
    file: str
    line: int
    function: str
    kind: str  # write | possible | remove
    sink: str
    target: str
    via: str = ""  # the writer helper, for a helper call
    reads_target: bool = False  # the same function also reads a results path (an in-place edit)
    notes: list[str] = field(default_factory=list)


class _Scope:
    def __init__(self, name: str, parent: _Scope | None, params: list[str]):
        self.name = name
        self.parent = parent
        self.params = params
        self.local: set[str] = set(params)


def _const(n: ast.AST) -> str | None:
    return n.value if isinstance(n, ast.Constant) and isinstance(n.value, str) else None


def _div_parts(n: ast.AST) -> list[ast.AST]:
    if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Div):
        return _div_parts(n.left) + _div_parts(n.right)
    return [n]


def _mode_of(c: ast.Call, pos: int) -> tuple[str | None, bool]:
    """(mode, known): the mode argument at `pos` or `mode=`; known False when it is not a constant."""
    node = c.args[pos] if len(c.args) > pos else next((k.value for k in c.keywords if k.arg == "mode"), None)
    if node is None:
        return "r", True
    s = _const(node)
    return (s, True) if s is not None else (None, False)


def _writes(mode: str | None) -> bool:
    return bool(mode and WRITE_MODE.search(mode))


def _module_of(path: str) -> str:
    return path[:-3].replace("/", ".").removesuffix(".__init__")


class FileAnalysis:
    """Taint labels per (scope, name) and per attribute key, to a fixpoint, then the sinks."""

    def __init__(self, path: str, tree: ast.Module, world: World):
        self.path = path
        self.module = _module_of(path)
        self.tree = tree
        self.world = world
        self.labels: dict[tuple[str, str], set[str]] = {}
        self.attr_labels: dict[str, set[str]] = {}
        self.arg_dests: set[str] = set()
        self.arg_defaults: list[tuple[_Scope, ast.AST, str]] = []  # scope, `default=` expression, dest
        self.scopes: dict[ast.AST, _Scope] = {}
        self.assigns: list[tuple[_Scope, ast.AST, ast.AST, str]] = []  # scope, target, value, how
        self.calls: list[tuple[_Scope, ast.Call]] = []
        self.returns: dict[str, list[tuple[_Scope, ast.AST]]] = {}
        self.functions: dict[str, ast.AST] = {}
        self.imports: dict[str, str] = {}  # local alias -> dotted module or module.attr
        self._data_names: set[str] | None = None
        self._collect()

    # --- collection -----------------------------------------------------------------------------
    def _collect(self) -> None:
        module_scope = _Scope("<module>", None, [])
        self.scopes[self.tree] = module_scope
        self._visit(self.tree, module_scope, "")

    def _visit(self, node: ast.AST, scope: _Scope, prefix: str) -> None:
        for child in ast.iter_child_nodes(node):
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda)):
                a = child.args
                params = [x.arg for x in a.posonlyargs + a.args + a.kwonlyargs]
                name = f"{prefix}{child.name}" if not isinstance(child, ast.Lambda) else f"{prefix}<lambda>"
                s = _Scope(name, scope, params)
                self.scopes[child] = s
                if not isinstance(child, ast.Lambda):
                    self.functions.setdefault(name, child)
                    self.functions.setdefault(child.name, child)
                    scope.local.add(child.name)
                pos = a.posonlyargs + a.args
                for arg, default in zip(pos[len(pos) - len(a.defaults) :], a.defaults, strict=True):
                    self.assigns.append((s, ast.Name(id=arg.arg), default, "default"))
                for arg, default in zip(a.kwonlyargs, a.kw_defaults, strict=True):
                    if default is not None:
                        self.assigns.append((s, ast.Name(id=arg.arg), default, "default"))
                for p in params:
                    self.labels.setdefault((name, p), set()).add(f"P:{name}:{p}")
                    if p in RESULTS_PARAMS:
                        self.labels[(name, p)].add(R)
                self._visit(child, s, name + ".")
                continue
            if isinstance(child, ast.ClassDef):
                scope.local.add(child.name)
                self._visit(child, scope, f"{prefix}{child.name}.")
                init = self.functions.get(f"{prefix}{child.name}.__init__")
                if init is not None:
                    self.functions.setdefault(child.name, init)  # a constructor call reaches __init__
                continue
            if isinstance(child, ast.Import):
                for al in child.names:
                    self.imports[al.asname or al.name.split(".")[0]] = (
                        al.name if al.asname else al.name.split(".")[0]
                    )
            elif isinstance(child, ast.ImportFrom) and child.module:
                for al in child.names:
                    self.imports[al.asname or al.name] = f"{child.module}.{al.name}"
            elif isinstance(child, ast.Assign):
                for t in child.targets:
                    self._assign(scope, t, child.value, "assign")
            elif (
                isinstance(child, (ast.AnnAssign, ast.AugAssign)) and child.value is not None
            ) or isinstance(child, ast.NamedExpr):
                self._assign(scope, child.target, child.value, "assign")
            elif isinstance(child, (ast.For, ast.AsyncFor, ast.comprehension)):
                self._assign(scope, child.target, child.iter, "iter")
            elif isinstance(child, ast.withitem) and child.optional_vars is not None:
                self._assign(scope, child.optional_vars, child.context_expr, "with")
            elif isinstance(child, ast.Return) and child.value is not None:
                self.returns.setdefault(scope.name, []).append((scope, child.value))
            elif isinstance(child, ast.Call):
                self.calls.append((scope, child))
                self._argparse(child, scope)
            self._visit(child, scope, prefix)

    def _assign(self, scope: _Scope, target: ast.AST, value: ast.AST, how: str) -> None:
        if isinstance(target, (ast.Tuple, ast.List)):
            if isinstance(value, (ast.Tuple, ast.List)) and len(value.elts) == len(target.elts):
                for t, v in zip(target.elts, value.elts, strict=True):
                    self._assign(scope, t, v, how)
            else:
                for t in target.elts:
                    self._assign(scope, t, value, how)
            return
        if isinstance(target, ast.Name):
            scope.local.add(target.id)
        self.assigns.append((scope, target, value, how))

    def _argparse(self, c: ast.Call, scope: _Scope) -> None:
        """Record a flag's `default=` expression with its destination and its scope. Whether that
        expression is a results path is NOT decided here: collection runs before the propagation pass,
        so a `default=` that is a bare NAME has no label yet, and resolving a name needs the enclosing
        scope. `_resolve_arg_dests` asks in `propagate` instead, where every name has its label."""
        f = c.func
        if not (isinstance(f, ast.Attribute) and f.attr == "add_argument"):
            return
        default = next((k.value for k in c.keywords if k.arg == "default"), None)
        if default is None:
            return
        dest = next((_const(k.value) for k in c.keywords if k.arg == "dest"), None)
        if dest is None:
            flags = [s for s in (_const(a) for a in c.args) if s]
            long = next((s for s in flags if s.startswith("--")), flags[0] if flags else "")
            dest = long.lstrip("-").replace("-", "_")
        if dest:
            self.arg_defaults.append((scope, default, dest))

    # --- labels -----------------------------------------------------------------------------------
    def _lookup(self, scope: _Scope | None, name: str) -> set[str]:
        s = scope
        while s is not None:
            if name in s.local:
                return self.labels.get((s.name, name), set())
            s = s.parent
        if ("<module>", name) in self.labels:
            return self.labels[("<module>", name)]
        return self.world.imported(self, name)

    def _results(self, n: ast.AST, scope: _Scope | None) -> bool:
        return R in self.label(n, scope)

    def label(self, n: ast.AST | None, scope: _Scope | None) -> set[str]:
        """The taint labels an expression carries: R when it is a results path, P:x when it is built
        from parameter x of the enclosing function."""
        if n is None:
            return set()
        if isinstance(n, ast.Constant):
            s = _const(n)
            return {R} if s is not None and RESULTS_TEXT.search(s) else set()
        if isinstance(n, ast.JoinedStr):
            out: set[str] = set()
            for v in n.values:
                if isinstance(v, ast.Constant) and _const(v) and RESULTS_TEXT.search(_const(v) or ""):
                    out.add(R)
                elif isinstance(v, ast.FormattedValue):
                    out |= self.label(v.value, scope)
            return out
        if isinstance(n, ast.Name):
            if n.id in RESULTS_CONSTANTS:
                return {R}
            return set(self._lookup(scope, n.id))
        if isinstance(n, ast.Attribute):
            if n.attr in RESULTS_CONSTANTS:
                return {R}
            out = set(self.attr_labels.get(ast.unparse(n), set()))
            if isinstance(n.value, ast.Name) and n.value.id in self.imports:
                out |= self.world.module_constant(self.imports[n.value.id], n.attr)
            if n.attr in self.arg_dests and isinstance(n.value, ast.Name):
                out.add(R)
            if n.attr in ("name", "stem", "suffix", "suffixes"):  # a string, not a path
                return out
            return out | self.label(n.value, scope)
        if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Div):
            parts = _div_parts(n)
            out = set()
            for i, p in enumerate(parts):
                out |= self.label(p, scope)
                if _const(p) == "results" and i > 0:
                    prev = parts[i - 1]
                    if (_const(prev) is not None and DATA_TEXT.search(_const(prev) or "")) or self._data_dir(
                        prev, scope
                    ):
                        out.add(R)
            return out
        if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Add):
            return self.label(n.left, scope) | self.label(n.right, scope)
        if isinstance(n, ast.Call):
            return self._call_label(n, scope)
        if isinstance(n, ast.IfExp):
            return self.label(n.body, scope) | self.label(n.orelse, scope)
        if isinstance(n, ast.BoolOp):
            out = set()
            for v in n.values:
                out |= self.label(v, scope)
            return out
        if isinstance(n, (ast.List, ast.Tuple, ast.Set)):
            out = set()
            for v in n.elts:
                out |= self.label(v, scope)
            return out
        if isinstance(n, ast.Dict):
            out = set()
            for v in n.values:
                out |= self.label(v, scope)
            return out
        if isinstance(n, (ast.ListComp, ast.SetComp, ast.GeneratorExp)):
            return self.label(n.elt, self.scopes.get(n, scope))
        if isinstance(n, ast.DictComp):
            return self.label(n.value, scope)
        if isinstance(n, ast.Subscript):
            return self.label(n.value, scope)
        if isinstance(n, ast.Starred):
            return self.label(n.value, scope)
        if isinstance(n, ast.NamedExpr):
            return self.label(n.value, scope)
        return set()

    def _data_dir(self, n: ast.AST, scope: _Scope | None, names: set[str] | None = None) -> bool:
        """Whether an expression names a data/ directory (so `x / "results"` is the registry)."""
        names = self.data_names() if names is None else names
        s = _const(n)
        if s is not None:
            return bool(DATA_TEXT.search(s))
        if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Div):
            return self._data_dir(_div_parts(n)[-1], scope, names)
        if isinstance(n, ast.Name):
            return n.id in names
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id in PATH_CALLS and n.args:
            return self._data_dir(n.args[-1], scope, names)
        return False

    def data_names(self) -> set[str]:
        """Names assigned a data/ directory anywhere in the file (flow-insensitive), to a fixpoint."""
        if self._data_names is not None:
            return self._data_names
        out: set[str] = set()
        while True:
            new = {
                t.id
                for _s, t, v, _how in self.assigns
                if isinstance(t, ast.Name) and t.id not in out and self._data_dir(v, None, out)
            }
            if not new:
                break
            out |= new
        self._data_names = out
        return out

    def _call_label(self, c: ast.Call, scope: _Scope | None) -> set[str]:
        f = c.func
        args = list(c.args) + [k.value for k in c.keywords]
        if isinstance(f, ast.Name) and f.id in PATH_CALLS:
            out = set()
            for a in args:
                out |= self.label(a, scope)
            consts = [_const(a) for a in c.args]
            if any(
                consts[i] == "results" and consts[i - 1] and DATA_TEXT.search(consts[i - 1] or "")
                for i in range(1, len(consts))
            ):
                out.add(R)
            return out
        if isinstance(f, ast.Attribute):
            owner = ast.unparse(f.value)
            if f.attr in OS_PATH_CALLS and owner in ("os.path", "path", "osp"):
                out = set()
                for a in args:
                    out |= self.label(a, scope)
                consts = [_const(a) for a in c.args]
                if any(
                    consts[i] == "results" and DATA_TEXT.search(consts[i - 1] or "")
                    for i in range(1, len(consts))
                ):
                    out.add(R)
                return out
            if f.attr == "fspath" and owner == "os":
                return self.label(c.args[0], scope) if c.args else set()
            if f.attr in PATH_METHODS:
                out = self.label(f.value, scope)
                if f.attr == "joinpath":
                    for a in args:
                        out |= self.label(a, scope)
                return out
        ret = self.world.returns_results(self, f)
        return {R} if ret else set()

    def resolve(self, f: ast.AST) -> tuple[str, str] | None:
        """(module, function) a call's callee names, when it is a function of an analysed file."""
        if isinstance(f, ast.Name):
            if f.id in self.functions:
                return (self.module, f.id)
            dotted = self.imports.get(f.id)
            if dotted and "." in dotted:
                mod, fn = dotted.rsplit(".", 1)
                return (mod, fn)
            return None
        if isinstance(f, ast.Attribute):
            owner = ast.unparse(f.value)
            if owner in ("self", "cls"):
                return (self.module, f.attr) if f.attr in self.functions else None
            dotted = self.imports.get(owner.split(".")[0])
            if dotted:
                rest = owner.split(".")[1:]
                return (".".join([dotted, *rest]), f.attr)
        return None

    def _resolve_arg_dests(self) -> bool:
        """Which argparse destinations read back a results path, decided with the labels of the pass
        this runs in. A `default=` is an ordinary expression and gets no special treatment: a name is
        resolved here exactly as a name on the right of an assignment is, in its own scope, and a dest
        found in one round taints `args.<dest>` in the next because `propagate` runs to a fixpoint.

        WHAT EACH PART STOPS, measured by removing it alone (tests/test_results_writers_guard.py):

        - running in the propagation pass instead of during collection. Remove it (decide in
          `_argparse` again) and three positives plus the reconstructed writer fail: a `default=` that
          is a NAME has no label yet while the tree is being walked, so it resolved to nothing. This is
          the part that stops the breach that happened -- data/results/astroreg_request_plan.json at
          723d802, written outside save_result through `default=DEFAULT_OUT`, with the census finding 13
          sites and none of them that one.
        - passing the enclosing scope instead of None. Remove it alone and exactly one positive fails,
          a default that is a name LOCAL to the enclosing function; the module constant is still caught,
          because `_lookup` falls back to ("<module>", name) with no scope at all. So this part is NOT
          what caught the astroreg shape and must not be credited with it: it is a guard against a
          different harm, a function-local default, and it stops nothing without the part above.
        """
        changed = False
        for scope, default, dest in self.arg_defaults:
            if dest not in self.arg_dests and self._results(default, scope):
                self.arg_dests.add(dest)
                changed = True
        return changed

    def propagate(self) -> bool:
        changed = self._resolve_arg_dests()
        for scope, target, value, _how in self.assigns:
            lab = self.label(value, scope)
            if not lab:
                continue
            if isinstance(target, ast.Subscript) and isinstance(target.value, ast.Name):
                target = target.value  # a container that holds a results path
            if isinstance(target, ast.Name):
                owner = scope
                while owner is not None and target.id not in owner.local:
                    owner = owner.parent
                key = ((owner or scope).name, target.id)
                cur = self.labels.setdefault(key, set())
            elif isinstance(target, (ast.Attribute, ast.Subscript)):
                base = target if isinstance(target, ast.Attribute) else target.value
                cur = self.attr_labels.setdefault(ast.unparse(base), set())
            else:
                continue
            if not lab <= cur:
                cur |= lab
                changed = True
        return changed

    # --- sinks ------------------------------------------------------------------------------------
    def sinks(self) -> list[tuple[_Scope, ast.Call, str, ast.AST | None, str]]:
        """(scope, call, sink, target, kind) for every call that writes or removes."""
        out = []
        for scope, c in self.calls:
            f = c.func
            name = f.attr if isinstance(f, ast.Attribute) else f.id if isinstance(f, ast.Name) else ""
            owner = ast.unparse(f.value) if isinstance(f, ast.Attribute) else ""
            a0 = c.args[0] if c.args else None
            a1 = c.args[1] if len(c.args) > 1 else next((k.value for k in c.keywords if k.arg == "dst"), None)
            if isinstance(f, ast.Attribute) and name in (
                "write_text",
                "write_bytes",
                "touch",
                "symlink_to",
                "hardlink_to",
            ):
                out.append((scope, c, name, f.value, "write"))
            elif isinstance(f, ast.Attribute) and name == "open" and owner.split(".")[-1] in FILE_OPENERS:
                mode, known = _mode_of(c, 1)
                if _writes(mode) or not known:
                    out.append((scope, c, f"{owner}.open", a0, "write" if known else "possible"))
            elif isinstance(f, ast.Attribute) and name == "open":
                mode, known = _mode_of(c, 0)
                if _writes(mode) or not known:
                    out.append((scope, c, ".open", f.value, "write" if known else "possible"))
            elif isinstance(f, ast.Name) and name == "open":
                mode, known = _mode_of(c, 1)
                if _writes(mode) or not known:
                    out.append((scope, c, "open", a0, "write" if known else "possible"))
            elif owner == "shutil" and name in ("copy", "copy2", "copyfile", "copytree", "move"):
                out.append((scope, c, f"shutil.{name}", a1, "write"))
            elif owner == "os" and name in ("replace", "rename", "link", "symlink"):
                out.append((scope, c, f"os.{name}", a1, "write"))
            elif (
                isinstance(f, ast.Attribute)
                and name in ("rename", "replace")
                and len(c.args) == 1
                and not c.keywords
            ):
                out.append((scope, c, f".{name}", a0, "write"))
            elif owner in ("np", "numpy") and name in ("save", "savez", "savez_compressed", "savetxt"):
                out.append((scope, c, f"np.{name}", a0, "write"))
            elif isinstance(f, ast.Attribute) and name in (
                "to_csv",
                "to_json",
                "to_parquet",
                "to_pickle",
                "savefig",
            ):
                out.append((scope, c, f".{name}", a0, "write"))
            elif isinstance(f, ast.Attribute) and name in ("unlink", "rmdir"):
                out.append((scope, c, f".{name}", f.value, "remove"))
            elif (owner == "os" and name in ("remove", "unlink")) or (owner == "shutil" and name == "rmtree"):
                out.append((scope, c, f"{owner}.{name}", a0, "remove"))
            elif owner == "subprocess" and name in SUBPROCESS:
                out.append((scope, c, f"subprocess.{name}", a0, "possible"))
        return out


class World:
    """Every analysed file, and the writer helpers and path-returning functions across them."""

    def __init__(self, root: Path = ROOT, files: list[str] | None = None):
        self.root = root
        if files is None:
            files = _git("ls-files", "--", *(f"{d}/*.py" for d in SOURCE_ROOTS), root=root).splitlines()
            files += _git(
                "ls-files",
                "--others",
                "--exclude-standard",
                "--",
                *(f"{d}/*.py" for d in SOURCE_ROOTS),
                root=root,
            ).splitlines()
        self.files: dict[str, FileAnalysis] = {}
        self.by_module: dict[str, FileAnalysis] = {}
        self.helpers: dict[tuple[str, str], set[str]] = {}  # (module, fn) -> parameters that reach a sink
        self.methods: dict[str, set[tuple[str, str]]] = {}  # method name -> helper methods of that name
        self.path_returns: set[tuple[str, str]] = set()
        for f in sorted(set(files)):
            p = root / f
            if not p.is_file():
                continue
            try:
                tree = ast.parse(p.read_text(), filename=f)
            except (SyntaxError, UnicodeDecodeError):
                continue
            fa = FileAnalysis(f, tree, self)
            self.files[f] = fa
            self.by_module[fa.module] = fa
        self._fixpoint()

    def imported(self, fa: FileAnalysis, name: str) -> set[str]:
        """An imported module-level name (`from m import NAME`) carries R when NAME is a results path in m."""
        dotted = fa.imports.get(name)
        if not dotted or "." not in dotted:
            return set()
        mod, attr = dotted.rsplit(".", 1)
        return self.module_constant(mod, attr)

    def module_constant(self, mod: str, attr: str) -> set[str]:
        src = self.by_module.get(mod)
        if src is None or src.path in SANCTIONED:
            return set()
        return {R} if R in src.labels.get(("<module>", attr), set()) else set()

    def returns_results(self, fa: FileAnalysis, f: ast.AST) -> bool:
        key = fa.resolve(f)
        return key in self.path_returns if key else False

    def _fixpoint(self) -> None:
        for _ in range(50):
            changed = False
            for fa in self.files.values():
                while fa.propagate():
                    changed = True
            for path, fa in self.files.items():
                if path in SANCTIONED:
                    continue
                for fname, rets in fa.returns.items():
                    key = (fa.module, fname.rsplit(".", 1)[-1])
                    if key not in self.path_returns and any(R in fa.label(v, s) for s, v in rets):
                        self.path_returns.add(key)
                        changed = True
                for scope, _c, sink, target, kind in fa.sinks():
                    if kind == "remove" or sink.startswith("subprocess"):
                        continue
                    labels = {x for x in fa.label(target, scope) if x.startswith("P:")}
                    changed |= self._register(fa, labels)
                for scope, c in fa.calls:
                    changed |= self._register(fa, self._helper_params(fa, scope, c))
            if not changed:
                return

    def _register(self, fa: FileAnalysis, labels: set[str]) -> bool:
        """Record that parameter x of function q reaches a sink, for each label P:q:x. A parameter of
        __init__ also makes the class a helper (a constructor call reaches it); a method is indexed by
        its name too, for calls on an object whose class is not resolved."""
        changed = False
        for lab in labels:
            _p, qual, param = lab.split(":", 2)
            if qual == "<module>":
                continue
            parts = qual.split(".")
            keys = [(fa.module, parts[-1])]
            if parts[-1] == "__init__" and len(parts) >= 2:
                keys.append((fa.module, parts[-2]))
            for key in keys:
                cur = self.helpers.setdefault(key, set())
                if param not in cur:
                    cur.add(param)
                    changed = True
            if len(parts) >= 2 and parts[-1] != "__init__":
                self.methods.setdefault(parts[-1], set()).add((fa.module, parts[-1]))
        return changed

    def _callee_args(self, fa: FileAnalysis, c: ast.Call) -> tuple[str, list[tuple[str, ast.AST]]] | None:
        """(callee, [(parameter, argument)]) for a call to a writer helper whose writing parameters it
        fills, or None. An attribute call on an object of an unresolved class is matched by method name."""
        key = fa.resolve(c.func)
        keys = [key] if key is not None and key in self.helpers else []
        if not keys and key is None and isinstance(c.func, ast.Attribute):
            keys = sorted(self.methods.get(c.func.attr, ()))
        found: list[str] = []
        hits: list[tuple[str, ast.AST]] = []
        for k in keys:
            callee_fa = self.by_module.get(k[0])
            fn = callee_fa.functions.get(k[1]) if callee_fa else None
            if fn is None or isinstance(fn, ast.Lambda):
                continue
            names = [x.arg for x in fn.args.posonlyargs + fn.args.args]
            if names and names[0] in ("self", "cls"):
                names = names[1:]
            pairs = [(names[i], a) for i, a in enumerate(c.args) if i < len(names)]
            pairs += [(kw.arg, kw.value) for kw in c.keywords if kw.arg]
            want = self.helpers.get(k, set())
            hit = [(p, a) for p, a in pairs if p in want]
            if hit:
                found.append(f"{k[0]}.{k[1]}")
                hits += [h for h in hit if h not in hits]
        if not found:
            return None
        return (found[0] if keys == [key] else "one of " + ", ".join(found) + " (by method name)"), hits

    def _helper_params(self, fa: FileAnalysis, scope: _Scope, c: ast.Call) -> set[str]:
        got = self._callee_args(fa, c)
        out: set[str] = set()
        for _p, a in got[1] if got else []:
            out |= {x for x in fa.label(a, scope) if x.startswith("P:")}
        return out

    def sites(self) -> list[Site]:
        out: list[Site] = []
        for path, fa in self.files.items():
            if path in SANCTIONED:
                continue
            reads = self._reading_scopes(fa)
            for scope, c, sink, target, kind in fa.sinks():
                if kind == "possible" and sink.startswith("subprocess"):
                    hit = not _git_read(target) and any(
                        R in fa.label(a, scope) for a in [*c.args, *(k.value for k in c.keywords)]
                    )
                elif kind == "possible":
                    hit = R in fa.label(target, scope)
                else:
                    hit = R in fa.label(target, scope)
                if hit:
                    out.append(
                        Site(
                            path,
                            c.lineno,
                            scope.name,
                            kind,
                            sink,
                            ast.unparse(target) if target else "",
                            reads_target=scope.name in reads,
                        )
                    )
            for scope, c in fa.calls:
                got = self._callee_args(fa, c)
                if not got:
                    continue
                callee, pairs = got
                hit = [(p, a) for p, a in pairs if R in fa.label(a, scope)]
                if hit:
                    out.append(
                        Site(
                            path,
                            c.lineno,
                            scope.name,
                            "write",
                            "helper",
                            ast.unparse(hit[0][1])[:120],
                            via=f"{callee}({', '.join(p for p, _a in hit)})",
                            reads_target=scope.name in reads,
                        )
                    )
        return sorted(out, key=lambda s: (s.file, s.line))

    def _reading_scopes(self, fa: FileAnalysis) -> set[str]:
        """Functions that also read a results path (read_text, json.load of open(...), load_result)."""
        out = set()
        for scope, c in fa.calls:
            f = c.func
            name = f.attr if isinstance(f, ast.Attribute) else f.id if isinstance(f, ast.Name) else ""
            if (
                name in ("read_text", "read_bytes")
                and isinstance(f, ast.Attribute)
                and R in fa.label(f.value, scope)
            ):
                out.add(scope.name)
            elif name == "open" and isinstance(f, ast.Name) and c.args and R in fa.label(c.args[0], scope):
                mode, _ = _mode_of(c, 1)
                if not _writes(mode):
                    out.add(scope.name)
            elif name == "load_result":
                out.add(scope.name)
        return out


#: git subcommands that read the repository and never write the working tree.
GIT_READS = {"show", "log", "diff", "cat-file", "ls-files", "ls-tree", "rev-parse", "status", "blame"}


def _git_read(cmd: ast.AST | None) -> bool:
    """A subprocess command list that runs a git subcommand which only reads (a `git show` of a
    committed result is a read, not a write)."""
    if not isinstance(cmd, (ast.List, ast.Tuple)) or len(cmd.elts) < 2:
        return False
    first, second = _const(cmd.elts[0]), _const(cmd.elts[1])
    if first != "git":
        return False
    if second == "-C" and len(cmd.elts) > 3:
        second = _const(cmd.elts[3])
    return second in GIT_READS


def _git(*args: str, root: Path = ROOT) -> str:
    return subprocess.run(["git", "-C", str(root), *args], capture_output=True, text=True, check=True).stdout


def find_sites(root: Path = ROOT, files: list[str] | None = None) -> list[Site]:
    """Every site under genomeos/ and scripts/ that writes (or may write) into data/results/ other
    than through save_result."""
    return World(root, files).sites()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    ap.add_argument("--json", action="store_true", help="print the sites as JSON")
    args = ap.parse_args(argv)
    sites = find_sites()
    if args.json:
        print(json.dumps([asdict(s) for s in sites], indent=1))
        return 0
    kinds: dict[str, int] = {}
    for s in sites:
        kinds[s.kind] = kinds.get(s.kind, 0) + 1
    print(f"sites writing into data/results without save_result: {len(sites)} {kinds}")
    print(f"files: {len({s.file for s in sites})}")
    for s in sites:
        extra = f" via {s.via}" if s.via else ""
        rd = " [reads a results path too]" if s.reads_target else ""
        print(f"  {s.kind:8s} {s.file}:{s.line} {s.function} {s.sink}({s.target[:70]}){extra}{rd}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
