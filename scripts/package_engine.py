# SPDX-License-Identifier: AGPL-3.0-or-later
"""Package the BioLang engine on its own and prove it runs with GenomeOS absent (ROADMAP area A, 2.0).

`tests/test_engine_boundary.py` shows that no engine file *imports* the application. That is a
different claim from the one the milestone makes, which is that the engine subtree is **complete**:
that nothing it needs is left behind when the rest of the project is not there. Only building it
proves that, so this builds it.

    uv run python scripts/package_engine.py [--out DIR] [--name biolang]

It copies the Apache-2.0 half — lang, ir, runtime, std and the three loose modules — into a package
of its own, rewrites its imports to the new name, gives it a pyproject of its own declaring Apache-2.0
and no dependencies, and then runs it **in a Python with no site-packages and no PYTHONPATH**, from a
directory that holds nothing else. In that interpreter GenomeOS cannot be imported at all, so anything
the engine still needs from it fails loudly rather than being quietly satisfied by the venv.

What it then runs is the toolchain's own test set, in the project's idiom: the standard library's
`.bio` modules testing themselves, plus `check`, `compile`, `run` and `test` over a program written
here, so all four verbs are exercised without the application.

It also carries the engine's **own pytest suite**: every test file under `tests/` whose GenomeOS
imports all fall inside the engine (and which loads nothing from `scripts/`) is copied into the
package with the same import rewrite the source gets, together with the `.bio`, `.bnet`, SBML and
grammar fixtures those tests read, each followed through its `import X.bio` lines. The ninth check
runs that suite in the same isolated interpreter; a missing fixture stops the session rather than
skipping a test. And it gets a version of its own (ENGINE_VERSION), because the engine and the
application are released on different cadences and cannot share one string.
"""

import argparse
import ast
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PACKAGES = ("lang", "ir", "runtime", "std")  # the same list the boundary test enforces
MODULES = ("coords", "version", "bio")

# The engine's own version, separate from the application's (genomeos/version.py, 1.0.0). It starts
# fresh at 0.1.0 rather than tracking the language spec (v0.4): the spec versions the language, which
# BioIR already records as `bioir_version`, and a package release has to move for a runtime fix that
# changes no syntax, so tying the two would force a spec bump for every bug fix or leave them lying.
ENGINE_VERSION = "0.1.0"

# what pytest needs in an interpreter with no site-packages: pytest and its runtime dependencies
# only, copied beside the build so the isolated run can import them and still cannot import GenomeOS
PYTEST_ONLY = ("pytest", "_pytest", "pluggy", "iniconfig", "packaging", "pygments", "py")

VERSION_PY = '''# SPDX-License-Identifier: Apache-2.0
"""The BioLang engine's own version.

The engine is released on its own cadence, so it does not share the GenomeOS application's version
string. It started fresh at 0.1.0 rather than tracking the language spec (v0.4): the spec is recorded
in every compiled module as `bioir_version`, and a package release has to be able to move for a
runtime fix that changes no syntax.
"""

__version__ = "{version}"
'''

CONFTEST = '''# SPDX-License-Identifier: Apache-2.0
"""Fixtures for the engine's own test suite: every file these tests read must be here.

A test whose fixture is missing is not a passing test, so a missing fixture stops the session with a
non-zero exit instead of letting a `skipif(not path.exists())` pass it over.
"""

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
FIXTURES = {fixtures!r}


def pytest_sessionstart(session):
    missing = [f for f in FIXTURES if not (ROOT / f).is_file()]
    if missing:
        pytest.exit(f"fixtures missing, the suite cannot run: {{missing}}", returncode=4)
'''

INIT = '''"""BioLang: a language for biology, its IR and the engines that run it.

    {name}.lang     parser: BioLang source to BioIR
    {name}.ir       BioIR: the types every engine reads
    {name}.runtime  the engines: networks, Boolean, located cells, the Body
    {name}.std      the standard library, written in BioLang

Apache-2.0. GenomeOS (AGPL-3.0-or-later) is an application built on this and is
not required to run it.
"""
# SPDX-License-Identifier: Apache-2.0

from {name}.version import __version__

__all__ = ["__version__"]
'''

PYPROJECT = """[project]
name = "{name}"
version = "{version}"
description = "BioLang: a language for biology, its intermediate representation, and the engines that run it"
readme = "README.md"
requires-python = ">=3.11"
license = "Apache-2.0"
license-files = ["LICENSE"]
dependencies = []

[project.scripts]
bio = "{name}.bio:main"

[dependency-groups]
dev = ["pytest>=8"]

[tool.pytest.ini_options]
testpaths = ["tests"]

[project.optional-dependencies]
compose = ["process-bigraph>=1.8.4"]
sbml = ["python-libsbml>=5.20"]

[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[tool.hatch.build.targets.wheel]
packages = ["{name}"]
"""

README = """# BioLang

The language, its IR and the engines that run it, packaged on its own.

    bio check FILE     compile, resolve references, report evidence and confidence
    bio compile FILE   BioLang to BioIR JSON
    bio run FILE       run it
    bio test PATH...   run every .bio file and evaluate its `# test:` lines
    bio repl           type BioLang, run it, inspect it

Apache-2.0, no dependencies. Built from the GenomeOS tree by
`scripts/package_engine.py`; GenomeOS is an application on this engine and is
not needed to run it.

## Tests

    pip install pytest    # or: uv sync --group dev
    pytest

`tests/` is the engine's own suite: {test_files} files and {test_count} tests,
the GenomeOS tests that touch nothing but the engine, with their imports
rewritten to `{name}`. The fixtures they read are under `data/` and `docs/`;
`tests/conftest.py` lists them and stops the run if one is missing rather than
letting a test skip. `bio test {name}/std` runs the standard library's own
`.bio` tests.

## Version

`{name}` is versioned on its own, starting at {version}; the GenomeOS
application has its own version and releases on its own cadence. It starts
fresh rather than tracking the language spec (v0.4): the spec is recorded in
every compiled module as `bioir_version`, and a package release has to be able
to move for a runtime fix that changes no syntax.
"""

# a program that exercises the language without any of the application's data
SMOKE = """module smoke
# test: rules >= 2
# test: unknowns == 0
gene A { locus: chr1:100-400(+); basal_rate: 1.0; evidence: curated "smoke"; confidence: 0.9 }
gene B { locus: chr1:900-1200(+); basal_rate: 0.2; evidence: curated "smoke"; confidence: 0.9 }
protein Ap { evidence: curated "smoke"; confidence: 0.9 }
protein Bp { evidence: curated "smoke"; confidence: 0.9 }
rule A produces Ap { rate: 1.0; evidence: curated "smoke"; confidence: 0.9 }
rule B produces Bp { rate: 1.0; evidence: curated "smoke"; confidence: 0.9 }
rule Ap inhibits Bp { rate: 0.5; evidence: curated "smoke"; confidence: 0.8 }
order along { members: A, B; axis: position; direction: increasing }
"""

# every module in the package, imported: `bio` touches most of the engine but not all of it, and a
# file that is missing something only fails when it is imported
EVERY = """
import importlib, pkgutil, sys
name = sys.argv[1]
package = importlib.import_module(name)
failed, extras = [], []
for info in pkgutil.walk_packages(package.__path__, name + "."):
    try:
        importlib.import_module(info.name)
    except ImportError as e:
        # an optional extra (process-bigraph, libsbml) is declared in the pyproject and is not the
        # engine reaching into the application; anything else is a hole in the package
        (extras if getattr(e, "name", "") and not e.name.startswith((name, "genomeos")) else failed).append(
            f"{info.name}: {e}"
        )
print("EXTRAS", "; ".join(extras) or "none")
assert not failed, failed
print("IMPORTED every module")
"""

ORGANISM = """module smoke.organism
# test: cells >= 6
# test: unknowns == 0
import bio.std.development
cell_type Skin { parent: PostMitotic; evidence: curated "smoke"; confidence: 0.8 }
organism Smoke { root: P0; cell_type: Blastomere; observe: count, fates
  assert: count at 0.5 min = 1
  evidence: curated "smoke"; confidence: 0.8 }
stage Early { from: 0 min; to: 40 min }
stage Late { from: 40 min }
timer cycle { duration: 10 min; when: generation = <2 }
decision grow { action: divide; when: cell_type = Blastomere, generation = <2 }
decision done { action: differentiate; when: stage = Late, generation = >=2; to: Skin }
order births { members: P0, P0a|P0p; axis: time; observe: birth }
"""

CHECKS = """
import importlib.util, sys
from pathlib import Path
name = sys.argv[1]
assert importlib.util.find_spec("genomeos") is None, "GenomeOS is importable: this is not an isolated run"
mod = importlib.import_module(name + ".bio")
loaded = sorted(m for m in sys.modules if m.split(".")[0] in (name, "genomeos"))
assert all(m.startswith(name) for m in loaded), loaded
version = importlib.import_module(name).__version__
assert version == sys.argv[2], f"the package says {version}, the build set {sys.argv[2]}"
import contextlib, io
out = io.StringIO()
with contextlib.redirect_stdout(out):
    try:
        mod.main(["--version"])
    except SystemExit:
        pass
assert out.getvalue().strip() == f"bio {version}", out.getvalue()
print("VERSION", version, "|", out.getvalue().strip())
print("LOADED", " ".join(loaded))
"""

# the engine's pytest suite, in the isolated interpreter: pytest and its dependencies are on the
# path from a directory that holds only them, and GenomeOS is checked to be absent before it starts
SUITE = """
import importlib.util, sys
sys.path.insert(0, sys.argv[1])
assert importlib.util.find_spec("genomeos") is None, "GenomeOS is importable: this is not an isolated run"
import pytest
sys.exit(pytest.main(["-q", "-rs", "-p", "no:cacheprovider", "tests"]))
"""


def engine_files() -> list[Path]:
    files = [p for pkg in PACKAGES for p in (ROOT / "genomeos" / pkg).rglob("*") if p.is_file()]
    files += [ROOT / "genomeos" / f"{m}.py" for m in MODULES]
    return [p for p in files if "__pycache__" not in p.parts]


def _rewrite(text: str, name: str) -> str:
    new = re.sub(r"\bgenomeos\b(?=\.[a-z_]+)", name, text)
    return new.replace("from genomeos import", f"from {name} import")


def _imported(path: Path) -> tuple[set[str], set[str]]:
    """(GenomeOS subpackages or modules a file imports, the other top-level names it imports)."""
    own, other = set(), set()
    for node in ast.walk(ast.parse(path.read_text())):
        if isinstance(node, ast.Import):
            names = [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom) and not node.level:
            names = (
                [f"genomeos.{a.name}" for a in node.names]
                if node.module == "genomeos"
                else [node.module or ""]
            )
        else:
            continue
        for n in names:
            parts = n.split(".")
            (own.add(parts[1] if len(parts) > 1 else "") if parts[0] == "genomeos" else other.add(parts[0]))
    return own, other


def engine_tests() -> list[Path]:
    """Test files that import GenomeOS only through the engine and load nothing from `scripts/`: the
    ones that test the engine and can travel with it. Found by reading the imports, not listed, so a
    new engine test joins the suite without anyone remembering to add it."""
    found = []
    for path in sorted((ROOT / "tests").glob("test_*.py")):
        own, other = _imported(path)
        if (
            own
            and own <= {*PACKAGES, *MODULES}
            and "scripts" not in other
            and not re.search(r"\bscripts\b[/.\"]", path.read_text())
        ):
            found.append(path)
    return found


def _bio_imports(path: Path) -> list[Path]:
    """The files a `.bio` fixture pulls in through `import`, resolved the way the parser does for
    anything outside the standard library (which travels with the package anyway)."""
    found = []
    for line in path.read_text().splitlines():
        m = re.match(r"\s*import\s+(\S+)", line)
        if not m or m.group(1).startswith("bio.std."):
            continue
        imp = m.group(1)
        rel = Path(imp) if imp.endswith(".bio") else Path(*imp.split(".")).with_suffix(".bio")
        if (path.parent / rel).is_file():
            found.append(path.parent / rel)
    return found


def fixtures(tests: list[Path]) -> list[str]:
    """Every file under data/ or docs/ the tests read, named as a path literal (`"data/x.bio"`) or
    built with `/` (`ROOT / "data" / "x.bio"`), plus the defaults engine functions they call read,
    closed over `.bio` imports."""
    texts = [t.read_text() for t in tests] + [p.read_text() for p in engine_files() if p.suffix == ".py"]
    names = set()
    for text in texts:
        names |= set(re.findall(r"[\"'](?:data|docs)/[\w./-]+\.\w+[\"']", text))
        for joined in re.findall(r'"(?:data|docs)"(?:\s*/\s*"[^"]+")+', text):
            names.add("/".join(re.findall(r'"([^"]+)"', joined)))
    todo = [ROOT / n.strip("\"'") for n in names]
    todo = [p for p in todo if p.is_file()]
    seen: set[Path] = set()
    while todo:
        p = todo.pop()
        if p in seen:
            continue
        seen.add(p)
        if p.suffix == ".bio":
            todo += _bio_imports(p)
    return sorted(str(p.relative_to(ROOT)) for p in seen)


def pytest_only(out: Path) -> Path:
    """pytest and its dependencies, copied beside the build (not into it) so the isolated run can
    import them without the venv, where GenomeOS is installed, ever being on the path."""
    import importlib.util

    where = out.with_name(out.name + "-pytest")
    shutil.rmtree(where, ignore_errors=True)
    where.mkdir(parents=True)
    for mod in PYTEST_ONLY:
        spec = importlib.util.find_spec(mod)
        if spec is None or spec.origin is None:
            raise SystemExit(f"{mod} is not installed; the engine's test suite needs pytest (uv sync)")
        src = Path(spec.origin)
        if src.name == "__init__.py":
            shutil.copytree(src.parent, where / src.parent.name, ignore=shutil.ignore_patterns("__pycache__"))
        else:
            shutil.copy2(src, where / src.name)
    return where


def build(out: Path, name: str) -> dict:
    """Copy the engine into `out/name`, rewrite its imports, and give it a pyproject of its own."""
    shutil.rmtree(out, ignore_errors=True)
    package = out / name
    package.mkdir(parents=True)
    copied, rewritten = [], 0
    for path in engine_files():
        target = package / path.relative_to(ROOT / "genomeos")
        target.parent.mkdir(parents=True, exist_ok=True)
        if path.suffix == ".py":
            text = path.read_text()
            new = _rewrite(text, name)
            rewritten += int(new != text)
            target.write_text(new)
        else:
            shutil.copy2(path, target)
        copied.append(str(target.relative_to(out)))
    (package / "__init__.py").write_text(INIT.format(name=name))
    # the engine's own version: genomeos/version.py stays the application's
    (package / "version.py").write_text(VERSION_PY.format(version=ENGINE_VERSION))
    (out / "pyproject.toml").write_text(PYPROJECT.format(name=name, version=ENGINE_VERSION))
    # the engine's own pytest suite, and every file it reads
    tests = engine_tests()
    (out / "tests").mkdir()
    for path in tests:
        (out / "tests" / path.name).write_text(_rewrite(path.read_text(), name))
    needed = fixtures(tests)
    for rel in needed:
        target = out / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        if rel.startswith("docs/"):  # generated from the engine's source, so it names the engine's modules
            target.write_text(_rewrite((ROOT / rel).read_text(), name))
        else:
            shutil.copy2(ROOT / rel, target)
    (out / "tests" / "conftest.py").write_text(CONFTEST.format(fixtures=needed))
    test_count = sum(
        len(re.findall(r"^def test_", p.read_text(), re.M)) for p in tests
    )  # functions; parametrised cases are counted when the suite runs
    (out / "README.md").write_text(
        README.format(name=name, version=ENGINE_VERSION, test_files=len(tests), test_count=test_count)
    )
    shutil.copy2(ROOT / "LICENSE-APACHE", out / "LICENSE")
    shutil.copy2(ROOT / "NOTICE", out / "NOTICE")
    (out / "smoke.bio").write_text(SMOKE)
    (out / "organism.bio").write_text(ORGANISM)
    # references to the other half that survive the rewrite: they import nothing, so they cannot break
    # a run, but in a package of its own they point at a module that is not there
    dangling = sorted(
        {
            m
            for path in package.rglob("*.py")
            for m in re.findall(rf"{name}\.[a-z_]+(?:\.[a-z_]+)*", path.read_text())
            if m.split(".")[1] not in (*PACKAGES, *MODULES)
        }
    )
    return {
        "files": len(copied),
        "rewritten": rewritten,
        "package": str(package),
        "dangling references": dangling,
        "version": ENGINE_VERSION,
        "test files": [p.name for p in tests],
        "test functions": test_count,
        "fixtures": needed,
    }


def isolated(out: Path, args: list[str]) -> subprocess.CompletedProcess:
    """Run in a Python with no site-packages and no PYTHONPATH, from the package directory: `-S` drops
    the venv where GenomeOS is installed, so nothing outside the standard library and this directory
    can satisfy an import. The current directory stays on the path, which is how the package is found."""
    env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "PYTHONHOME")}
    return subprocess.run([sys.executable, "-S", *args], capture_output=True, text=True, cwd=out, env=env)


def verify(out: Path, name: str) -> list[dict]:
    results = []

    def check(label: str, proc: subprocess.CompletedProcess, want: str = "") -> None:
        ok = proc.returncode == 0 and (want in proc.stdout if want else True)
        results.append(
            {
                "check": label,
                "ok": ok,
                "output": " | ".join((proc.stdout.strip().splitlines() or [""])[-3:])[:400],
                "error": proc.stderr.strip().splitlines()[-1][:200] if proc.stderr.strip() else "",
            }
        )

    check(
        "GenomeOS is not importable, the engine imports only itself, and bio --version is its own",
        isolated(out, ["-c", CHECKS, name, ENGINE_VERSION]),
        f"bio {ENGINE_VERSION}",
    )
    cli = ["-m", f"{name}.bio"]
    check("bio check", isolated(out, [*cli, "check", "smoke.bio"]), "entities")
    check("bio compile", isolated(out, [*cli, "compile", "smoke.bio"]), "bioir_version")
    check("bio run", isolated(out, [*cli, "run", "smoke.bio", "--hours", "4"]), "Ap")
    check("every engine module imports", isolated(out, ["-c", EVERY, name]), "IMPORTED every module")
    check("bio test (the program written here)", isolated(out, [*cli, "test", "smoke.bio"]), "checks passed")
    check(
        "bio test (an organism, so the Body runs too)",
        isolated(out, [*cli, "test", "organism.bio"]),
        "checks passed",
    )
    check(
        "bio test (the standard library testing itself)",
        isolated(out, [*cli, "test", f"{name}/std"]),
        "checks passed",
    )
    results.append(suite(out))
    return results


def suite(out: Path) -> dict:
    """The ninth check: the engine's own pytest suite, in the interpreter where GenomeOS cannot be
    imported. It passes only if nothing failed or errored, at least one test passed, and every skip is
    an optional extra that is not installed (process-bigraph, libsbml), never a missing file."""
    proc = isolated(out, ["-c", SUITE, str(pytest_only(out))])
    lines = proc.stdout.strip().splitlines()
    tail = lines[-1] if lines else ""
    counts = {
        k: int(v) for v, k in re.findall(r"(\d+) (passed|failed|skipped|errors?|xfailed|xpassed)", tail)
    }
    reasons = [ln for ln in lines if ln.startswith("SKIPPED")]
    bad_skips = [r for r in reasons if "not installed" not in r]
    ok = (
        proc.returncode == 0
        and counts.get("passed", 0) > 0
        and not counts.get("failed")
        and not counts.get("error")
        and not counts.get("errors")
        and not bad_skips
    )
    return {
        "check": "the engine's own pytest suite",
        "ok": ok,
        "output": tail[:400],
        "error": "; ".join(bad_skips)[:400]
        or ("" if ok else " | ".join(lines[-5:])[:400] or proc.stderr.strip()[-400:]),
        "passed": counts.get("passed", 0),
        "failed": counts.get("failed", 0),
        "skipped": counts.get("skipped", 0),
        "skip reasons": reasons,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default="build/biolang", help="where to build the package")
    ap.add_argument("--name", default="biolang", help="the package name to build it under")
    args = ap.parse_args()
    out = Path(args.out).resolve()
    built = build(out, args.name)
    print(
        f"built {args.name} {built['version']} from {built['files']} files"
        f" ({built['rewritten']} with imports rewritten)"
    )
    print(
        f"  its own suite: {len(built['test files'])} test files, {built['test functions']} test functions,"
        f" {len(built['fixtures'])} fixtures"
    )
    for ref in built["dangling references"]:
        print(f"  note: {ref} is named in a docstring but is not part of the engine")
    results = verify(out, args.name)
    for r in results:
        print(f"  {'ok  ' if r['ok'] else 'FAIL'} {r['check']}")
        if not r["ok"]:
            print(f"       {r['error'] or r['output']}")
        elif "passed" in r:
            print(f"       {r['output']}")
    ok = all(r["ok"] for r in results)
    print(f"\nthe engine runs with the application absent: {ok}")
    try:
        from genomeos.results import save_result

        print(
            "wrote",
            save_result(
                "engine_package",
                {"name": args.name, "built": built, "checks": results, "passed": ok},
            ),
        )
    except Exception as e:  # noqa: BLE001  (the report is a convenience, not the measurement)
        print(f"(no result file: {e})", file=sys.stderr)
    print(json.dumps({"passed": ok}))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
