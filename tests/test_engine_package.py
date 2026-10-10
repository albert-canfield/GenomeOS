# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The engine packaged on its own, and run with GenomeOS absent (ROADMAP area A, milestone 2.0).

`test_engine_boundary.py` shows that no engine file *imports* the application. This is the other
half of the claim, and the one the milestone actually makes: that the engine subtree is **complete**,
so nothing it needs is left behind when the rest of the project is not there. The only way to know
that is to build it and run it somewhere GenomeOS does not exist, which is what
`scripts/package_engine.py` does and what this test runs.

It is a licence check as much as a design one: the engine is Apache-2.0 and the application
AGPL-3.0-or-later (LICENSING.md D40), and a package that could not be built without the application
could not honestly be distributed under Apache.
"""

from __future__ import annotations

import json
import sys
import tomllib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

package_engine = pytest.importorskip("package_engine")


@pytest.fixture(scope="module")
def built(tmp_path_factory) -> tuple[Path, dict, list[dict]]:
    out = tmp_path_factory.mktemp("engine")
    info = package_engine.build(out, "biolang")
    return out, info, package_engine.verify(out, "biolang")


def test_the_engine_builds_into_a_package_of_its_own(built):
    _, info, _ = built
    assert info["files"] > 20, "the engine file list looks wrong"
    assert info["rewritten"] > 10, "nothing had its imports rewritten; the copy is suspect"


def test_nothing_in_the_engine_still_names_the_application(built):
    """A docstring that names an application module imports nothing and breaks no run, but in a
    package of its own it points at something that is not there."""
    _, info, _ = built
    assert info["dangling references"] == []


def test_it_runs_with_genomeos_not_importable(built):
    """Every check the packaging script makes, in a Python with no site-packages: the venv where
    GenomeOS is installed is not on the path, so an engine that still needed it would fail here."""
    _, _, checks = built
    failed = [c for c in checks if not c["ok"]]
    assert not failed, "\n".join(f"{c['check']}: {c['error'] or c['output']}" for c in failed)


def test_the_four_verbs_and_both_runtimes_are_covered(built):
    """The packaged test set is the toolchain's own: all four verbs, the standard library testing
    itself, and an organism so the Body runs and not only the network."""
    _, _, checks = built
    labels = " ".join(c["check"] for c in checks)
    for verb in ("bio check", "bio compile", "bio run", "bio test"):
        assert verb in labels
    assert "the Body runs" in labels and "standard library" in labels
    assert "every engine module imports" in labels


def test_the_package_declares_apache_and_no_dependencies(built):
    out, _, _ = built
    pyproject = (out / "pyproject.toml").read_text()
    assert 'license = "Apache-2.0"' in pyproject
    assert "dependencies = []" in pyproject  # the engine core takes nothing; the heavy engines are extras
    assert "Apache License" in (out / "LICENSE").read_text()


#: Every third-party module an engine file imports, with the extra that installs it. The engine core
#: takes nothing, so each is imported inside the function that needs it.
ENGINE_EXTRAS = {"process_bigraph": "compose", "roadrunner": "roadrunner"}


def test_every_extra_the_engine_imports_is_declared_as_the_application_declares_it(built):
    """runtime/compose.py and runtime/roadrunner.py import process-bigraph and libroadrunner lazily, so
    the package runs without them. Its own pyproject still has to offer each as an extra, under the
    name and the requirement the application's pyproject gives it, or `pip install biolang[roadrunner]`
    installs nothing. A new third-party import in the engine fails here until it is declared."""
    out, _, _ = built
    imported = set()
    for path in (out / "biolang").rglob("*.py"):
        _, other = package_engine._imported(path)
        imported |= {m for m in other if m and m != "biolang" and m not in sys.stdlib_module_names}
    assert imported == set(ENGINE_EXTRAS), f"the engine imports {sorted(imported)}"
    packaged = tomllib.loads((out / "pyproject.toml").read_text())["project"]["optional-dependencies"]
    application = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]["optional-dependencies"]
    for module, extra in ENGINE_EXTRAS.items():
        assert packaged.get(extra) == application[extra], (
            f"{module} is installed by the application's `{extra}` extra {application[extra]}, "
            f"and the package declares {packaged.get(extra)} for it"
        )


def test_the_engine_carries_its_own_pytest_suite_and_it_passes_in_isolation(built):
    """Milestone 2.0 asks for two test suites. The engine's is the GenomeOS tests that touch nothing
    but the engine, rewritten to `biolang`, run as the ninth check in the interpreter where GenomeOS
    cannot be imported. Nothing may fail, and a skip counts only for an optional extra that is not
    installed, never for a missing fixture."""
    out, info, checks = built
    suite = next(c for c in checks if c["check"] == "the engine's own pytest suite")
    assert suite["ok"], suite["error"] or suite["output"]
    assert suite["failed"] == 0 and suite["passed"] >= 140, suite["output"]
    assert all("not installed" in r for r in suite["skip reasons"]), suite["skip reasons"]
    assert len(info["test files"]) >= 24 and "test_lang_v04.py" in info["test files"]
    assert len(checks) == 9
    assert "genomeos" not in (out / "tests" / "test_lang.py").read_text()
    assert 'dev = ["pytest>=8"]' in (out / "pyproject.toml").read_text()


def test_every_fixture_the_suite_reads_travels_with_it(built):
    """The fixtures are found from the tests and followed through `import X.bio`, and the conftest
    that ships with them stops the run if one is missing."""
    out, info, _ = built
    for rel in (
        "data/demo/repressilator.bio",
        "data/organisms/celegans/embryo.bio",
        "data/organisms/human/body.bio",
    ):
        assert rel in info["fixtures"]
    assert "data/organisms/celegans/lineage_larva.bio" in info["fixtures"]  # reached only by import
    assert all((out / f).is_file() for f in info["fixtures"])
    assert "pytest.exit" in (out / "tests" / "conftest.py").read_text()


def test_the_engine_has_a_version_of_its_own(built):
    """Two packages on two cadences cannot share one string: `bio --version` in the package prints
    the engine's version, and genomeos/version.py stays the application's."""
    out, info, _ = built
    from genomeos.version import __version__ as application

    assert info["version"] == package_engine.ENGINE_VERSION != application
    assert f'version = "{info["version"]}"' in (out / "pyproject.toml").read_text()
    assert f'__version__ = "{info["version"]}"' in (out / "biolang" / "version.py").read_text()
    assert info["version"] in (out / "README.md").read_text()


@pytest.fixture(scope="module")
def described(built) -> list[dict]:
    out, _, _ = built
    return package_engine.describe(out, "biolang")


def test_the_package_describes_its_own_language(built, described):
    """The grammar and every specification travel inside the package, so the wheel carries them; the
    grammar is re-rendered by the packaged parser in isolation and must match what was copied."""
    out, info, _ = built
    checks = described
    assert len(checks) == 3
    docs = out / "biolang" / "docs"
    assert set(info["docs"]) == {p.name for p in (ROOT / "docs").glob("BIOLANG-*.md")}
    assert {"BIOLANG-GRAMMAR.md", "BIOLANG-v0.1.md", "BIOLANG-v0.3.md", "BIOLANG-v0.4-ECONOMY.md"} <= set(
        info["docs"]
    )
    assert all((docs / d).is_file() for d in info["docs"])
    assert "genomeos" not in (docs / "BIOLANG-GRAMMAR.md").read_text()
    assert "Creative Commons Attribution 4.0" in (docs / "README.md").read_text()
    assert 'packages = ["biolang"]' in (out / "pyproject.toml").read_text()  # docs/ is inside it
    c = next(c for c in checks if c["check"].startswith("the grammar and the specs"))
    assert c["ok"], c["error"] or c["output"]


def test_bio_repl_runs_in_isolation(described):
    c = next(c for c in described if c["check"].startswith("bio repl"))
    assert c["ok"], c["error"] or c["output"]


def test_import_protein_without_genomeos_is_a_named_error(described):
    """The engine ships no proteome: the import must stop, exit 2, name the missing resolver and how to
    supply one, with no traceback; not an empty module, not a network fetch."""
    c = next(c for c in described if c["check"].startswith("import protein:"))
    assert c["ok"], c["error"] or c["output"]
    assert "the engine ships none" in c["error"]


# ---------------------------------------------------------------- what the writer records, and where


def test_the_recorded_package_path_names_no_machine(built):
    """A build outside the repository (here a pytest temporary directory) is recorded relative to its
    `--out` with the flag set: the committed result named `/Users/...` until 2026-10-10, and a scratch
    rebuild named a `/private/tmp` scratchpad that no longer exists."""
    out, info, _ = built
    assert info["package"] == "biolang"
    assert info["package outside repository"] is True
    assert str(out) not in json.dumps(info) and str(out.resolve()) not in json.dumps(info)


def test_a_package_inside_the_repository_is_recorded_relative_to_it():
    out = package_engine.DEFAULT_OUT
    assert package_engine.recorded_path(out / "biolang", out) == {
        "package": "build/biolang/biolang",
        "package outside repository": False,
    }


def _stub_run(monkeypatch) -> list[dict]:
    """main() with the build and the checks stubbed out, and save_result recording its calls, so the
    test asks only where a result would go and writes nothing anywhere."""
    calls: list[dict] = []
    built = {
        "files": 1,
        "rewritten": 1,
        "package": "x",
        "package outside repository": True,
        "dangling references": [],
        "version": package_engine.ENGINE_VERSION,
        "test files": [],
        "test functions": 0,
        "fixtures": [],
        "docs": [],
    }
    monkeypatch.setattr(package_engine, "build", lambda out, name: built)
    monkeypatch.setattr(package_engine, "verify", lambda out, name: [{"check": "c", "ok": True}])
    monkeypatch.setattr(package_engine, "describe", lambda out, name: [])
    import genomeos.results

    def save_result(name, payload, results_dir=genomeos.results.RESULTS_DIR, **kw):
        calls.append({"name": name, "payload": payload, "results_dir": Path(results_dir)})
        return Path(results_dir) / f"{name}.json"

    monkeypatch.setattr(genomeos.results, "save_result", save_result)
    return calls


@pytest.mark.parametrize("extra", [["--out", "SCRATCH"], ["--name", "other"]])
def test_a_scratch_run_never_writes_the_canonical_result(monkeypatch, tmp_path, capsys, extra):
    """Any --out or --name but the default is a scratch run: its result goes beside its own build and
    never to data/results, which is how the uncommitted 2026-10-02 scratch rebuild replaced the
    canonical result in the shared tree."""
    calls = _stub_run(monkeypatch)
    argv = [str(tmp_path / "scratch") if a == "SCRATCH" else a for a in extra]
    assert package_engine.main(argv) == 0
    [call] = calls
    assert call["results_dir"].resolve() != (ROOT / "data" / "results").resolve()
    assert call["results_dir"] != Path("data/results")
    assert call["payload"]["canonical run"] is False
    assert "scratch run" in capsys.readouterr().out


def test_the_default_run_writes_the_canonical_result(monkeypatch):
    calls = _stub_run(monkeypatch)
    assert package_engine.main([]) == 0
    [call] = calls
    assert call["results_dir"] == Path("data/results")
    assert call["payload"]["canonical run"] is True
