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

import sys
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
