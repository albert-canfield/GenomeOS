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
"""

import argparse
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
print("VERSION", importlib.import_module(name).__version__)
print("LOADED", " ".join(loaded))
"""


def engine_files() -> list[Path]:
    files = [p for pkg in PACKAGES for p in (ROOT / "genomeos" / pkg).rglob("*") if p.is_file()]
    files += [ROOT / "genomeos" / f"{m}.py" for m in MODULES]
    return [p for p in files if "__pycache__" not in p.parts]


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
            new = re.sub(r"\bgenomeos\b(?=\.[a-z_]+)", name, text)
            new = new.replace("from genomeos import", f"from {name} import")
            rewritten += int(new != text)
            target.write_text(new)
        else:
            shutil.copy2(path, target)
        copied.append(str(target.relative_to(out)))
    (package / "__init__.py").write_text(INIT.format(name=name))
    version = re.search(r'__version__\s*=\s*"([^"]+)"', (ROOT / "genomeos" / "version.py").read_text())
    (out / "pyproject.toml").write_text(
        PYPROJECT.format(name=name, version=version.group(1) if version else "0.0.0")
    )
    (out / "README.md").write_text(README)
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
        "GenomeOS is not importable, and the engine imports only itself",
        isolated(out, ["-c", CHECKS, name]),
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
    return results


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default="build/biolang", help="where to build the package")
    ap.add_argument("--name", default="biolang", help="the package name to build it under")
    args = ap.parse_args()
    out = Path(args.out).resolve()
    built = build(out, args.name)
    print(f"built {args.name} from {built['files']} files ({built['rewritten']} with imports rewritten)")
    for ref in built["dangling references"]:
        print(f"  note: {ref} is named in a docstring but is not part of the engine")
    results = verify(out, args.name)
    for r in results:
        print(f"  {'ok  ' if r['ok'] else 'FAIL'} {r['check']}")
        if not r["ok"]:
            print(f"       {r['error'] or r['output']}")
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
