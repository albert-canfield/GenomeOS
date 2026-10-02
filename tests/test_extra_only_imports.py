# SPDX-License-Identifier: AGPL-3.0-or-later
"""No test may need a package that only an optional extra installs without saying which extra.

The defect this closes (2026-10-02): ten tests passed, for as long as they existed, only because
somebody had once installed pandas and openpyxl into this checkout by hand. Seven passed transitively --
CI's `test` job installs `--extra compose`, and eight locked packages pull pandas in behind it -- and
three needed openpyxl, which arrives only with biolearn in the `clocks` extra, which no CI job installs,
so those three were never in CI at all. Both quantities the gates watched were preserved the whole time:
the tests passed, and the count of passing tests went up.

Three mechanisms, in the order they fire:

    the `test-bare` CI job   installs the dev group and no extras, so a package nobody declared is a red
                             run and not a silent pass
    `requires_extra`         a test that genuinely needs an extra declares which one and skips by that
                             declaration, naming it (tests/conftest.py)
    this file                `test_no_test_module_needs_an_undeclared_package` fails on an import-time
                             import of a package the bare environment does not hold, before anyone has
                             to notice a skip

The last one is the one that would have caught openpyxl, and it is checked in both directions below: it
is shown firing on a planted module, and shown clearing the same module once the marker is added. A
detector shown only to fire is half-tested.
"""

from __future__ import annotations

import subprocess
import sys
import tomllib
from pathlib import Path

import extras_lock
import pytest

ROOT = Path(__file__).resolve().parents[1]
TESTS = ROOT / "tests"

#: the real site the project reaches openpyxl through, as a test module reaches it: by name, in a string
OPENPYXL_SITE = ("tests/test_code_cleanliness_shared.py", "scripts/astroreg_register.py")

#: the real site the project reaches pandas through
PANDAS_SITE = ("tests/test_s5_rate_ranges.py", "scripts/s5_rate_ranges.py")


@pytest.fixture
def counterfactual_lock(tmp_path):
    """A lock that is this project's, minus a declaration: `extras_lock` reads it instead of uv.lock.

    `build(...)` writes a copy of uv.lock whose dev group no longer names the given distributions, which
    is the tree as it stood before 8c4dbe0 declared pandas and openpyxl. Nothing else changes: the
    packages are still in the file, still reachable through `clocks` and seven other paths, exactly as
    they were when the ten tests were passing for the wrong reason.
    """
    text = (ROOT / "uv.lock").read_text()
    original = extras_lock.LOCK

    def build(*distributions: str) -> Path:
        document = tomllib.loads(text)
        lines = text.splitlines(keepends=True)
        # only inside this project's own package block: `{ name = "pandas" },` is also how eight other
        # packages declare their dependency on it, and deleting those would remove pandas from the
        # extras too -- a lock in which nothing needs pandas, which is not the lock that was there.
        first = next(i for i, line in enumerate(lines) if line.strip() == f'name = "{extras_lock.PROJECT}"')
        last = next(i for i, line in enumerate(lines[first:], first) if line.startswith("[[package]]"))
        for distribution in distributions:
            assert distribution in extras_lock.bare_distributions(), (
                f"{distribution} is not declared today, so removing it proves nothing"
            )
            target = f'{{ name = "{distribution}" }},'
            kept = [line for i, line in enumerate(lines) if not (first < i < last and line.strip() == target)]
            assert len(kept) == len(lines) - 1, (
                f"one dev-group line for {distribution}, not {len(lines) - len(kept)}"
            )
            lines = kept
            last -= 1
        written = tmp_path / "uv.lock"
        written.write_text("".join(lines))
        after = tomllib.loads(written.read_text())
        assert {package["name"] for package in after["package"]} == {
            package["name"] for package in document["package"]
        }, "the counterfactual must drop a declaration, not a package"
        return written

    def use(lock: Path) -> None:
        extras_lock.LOCK = lock
        extras_lock.reset_caches()

    yield build, use
    extras_lock.LOCK = original
    extras_lock.reset_caches()


def _plant(tmp_path: Path, body: str) -> Path:
    path = tmp_path / "test_planted.py"
    path.write_text(body)
    return path


# --------------------------------------------------------------------------------------------------
# what "reachable only through an extra" is derived from


def test_an_extra_adds_nothing_when_the_dev_group_already_declares_it():
    """`analysis` is numpy and the dev group is too, so `analysis` is extra-only for no package.

    This is the subtraction that makes the rest honest. Without it numpy would read as extra-only, every
    test that imports it would want a marker, and the marker would be a lie: numpy is declared.
    """
    assert extras_lock.extra_distributions()["analysis"] == frozenset()
    assert "numpy" in extras_lock.bare_distributions()
    assert extras_lock.extras_for_module("numpy") == ()


def test_the_two_packages_the_ten_tests_needed_are_declared_now():
    assert {"pandas", "openpyxl"} <= extras_lock.bare_distributions()
    assert extras_lock.extras_for_module("pandas") == ()
    assert extras_lock.extras_for_module("openpyxl") == ()


def test_an_extras_closure_follows_the_locks_own_edges():
    """`compose` is one declared requirement, process-bigraph, and it drags in 23 more packages.

    Eight of those was the whole defect: the `test` job asks for one extra and gets a closure, and a
    closure is not a declaration. The sets are read from the lock, so a new extra needs nobody to
    remember this file.
    """
    compose = extras_lock.extra_distributions()["compose"]
    assert "process-bigraph" in compose
    assert "bigraph-schema" in compose, "a transitive requirement of the extra is part of the extra"
    assert extras_lock.extras_for_module("process_bigraph") == ("compose",)
    assert set(extras_lock.known_extras()) == set(extras_lock.extra_distributions())


def test_a_package_several_extras_reach_names_them_all():
    """scipy arrives with clocks, compose, hla and predict; the message must not pick one."""
    assert len(extras_lock.extras_for_module("scipy")) > 1


# --------------------------------------------------------------------------------------------------
# the detector, both directions


PLANTED_IMPORT = '''\
"""A planted test module that needs an extra and does not say so."""

import process_bigraph


def test_nothing():
    assert process_bigraph
'''

PLANTED_MARKED = '''\
"""The same planted module, declaring the extra it needs."""

import pytest

pytestmark = pytest.mark.requires_extra("compose")

import process_bigraph  # noqa: E402


def test_nothing():
    assert process_bigraph
'''


def test_a_planted_module_that_needs_an_extra_without_saying_so_is_named(tmp_path):
    (violation,) = extras_lock.scan([_plant(tmp_path, PLANTED_IMPORT)])
    assert violation.module == "process_bigraph"
    assert violation.extras == ("compose",)
    assert "requires_extra('compose')" in str(violation)


def test_the_marker_clears_the_same_planted_module(tmp_path):
    """The other direction: the detector must stop firing once the declaration is there.

    Same file, same import, one marker added. A detector that cannot be satisfied is not a check but a
    wall, and the first person who met it would delete it.
    """
    assert extras_lock.scan([_plant(tmp_path, PLANTED_MARKED)]) == []


def test_a_marker_naming_the_wrong_extra_does_not_clear_it(tmp_path):
    body = PLANTED_MARKED.replace('requires_extra("compose")', 'requires_extra("clocks")')
    (violation,) = extras_lock.scan([_plant(tmp_path, body)])
    assert violation.module == "process_bigraph"


def test_an_import_inside_a_function_is_not_named(tmp_path):
    """pyproject's own rule: "imported inside the function that needs it ... the pattern to keep".

    Such an import cannot stop a module from collecting, and this project uses it everywhere: 730 sites
    across 66 test closures reach alphagenome or mhcflurry that way. Calling one of those functions in a
    bare environment still raises; what this check is about is import time, where the failure reads as
    "no such tests" rather than as a missing dependency.
    """
    body = "def test_nothing():\n    import process_bigraph\n\n    assert process_bigraph\n"
    assert extras_lock.scan([_plant(tmp_path, body)]) == []


def test_a_guarded_import_is_not_named(tmp_path):
    body = (
        "try:\n    import process_bigraph\n\n    AVAILABLE = True\n"
        "except ImportError:  # pragma: no cover\n    AVAILABLE = False\n\n\n"
        "def test_nothing():\n    assert AVAILABLE in (True, False)\n"
    )
    assert extras_lock.scan([_plant(tmp_path, body)]) == []


def test_a_fallback_import_in_the_handler_is_named(tmp_path):
    """`except ImportError: import other` runs `other` unprotected; only the try body is a guard."""
    body = (
        "try:\n    import json\nexcept ImportError:\n    import process_bigraph as json\n\n\n"
        "def test_nothing():\n    assert json\n"
    )
    (violation,) = extras_lock.scan([_plant(tmp_path, body)])
    assert violation.module == "process_bigraph"


def test_a_package_no_extra_installs_is_named_as_declared_nowhere(tmp_path):
    """The other half: an import the project declares in no extra either is still a finding.

    The module-name map is the weak part of the derivation -- nothing in a lock says scikit-learn is
    imported as `sklearn` -- so an import nothing accounts for is reported rather than assumed innocent.
    It is then a question for a person: declare it, or do not import it.
    """
    body = "import nobody_declares_this\n\n\ndef test_nothing():\n    assert nobody_declares_this\n"
    (violation,) = extras_lock.scan([_plant(tmp_path, body)])
    assert violation.extras == ()
    assert "declares it nowhere" in str(violation)


# --------------------------------------------------------------------------------------------------
# the defect itself, as a counterfactual


def test_the_detector_names_the_openpyxl_and_pandas_tests_if_their_declaration_goes(counterfactual_lock):
    """The ten tests, on the lock as it stood before 8c4dbe0: both sites named, through real files.

    openpyxl is the hard one, and the reason the closure follows strings as well as imports.
    `tests/test_code_cleanliness_shared.py` writes `"scripts.astroreg_register"` in a tuple and imports
    it with importlib; `scripts/astroreg_register.py:28` is where `import openpyxl` is written. No import
    statement in any test module mentions openpyxl at all.

    This is also the guard against the declaration being removed again: take pandas out of the dev group
    and seven tests are named here, before anyone reads a green run that owed its colour to a
    hand-installed package.
    """
    build, use = counterfactual_lock
    use(build("pandas", "openpyxl"))

    openpyxl_test, openpyxl_source = OPENPYXL_SITE
    pandas_test, pandas_source = PANDAS_SITE
    violations = extras_lock.scan([ROOT / openpyxl_test, ROOT / pandas_test])
    found = {(v.test.name, v.module, v.through.relative_to(ROOT).as_posix()) for v in violations}

    assert (Path(openpyxl_test).name, "openpyxl", openpyxl_source) in found
    assert (Path(pandas_test).name, "pandas", pandas_source) in found
    assert all(v.extras for v in violations if v.module == "openpyxl"), (
        "openpyxl is reachable through the clocks extra, so the message must name it"
    )


def test_the_counterfactual_drops_a_declaration_and_not_a_package(counterfactual_lock):
    """Or the test above would prove nothing: a lock missing the package proves a different claim."""
    build, use = counterfactual_lock
    use(build("pandas"))
    assert "pandas" not in extras_lock.bare_distributions()
    assert "openpyxl" in extras_lock.bare_distributions(), "only what was asked for is removed"
    # and this is why the `test` job could not see it: undeclared, pandas still arrives with any of four
    # extras, so one `--extra compose` was enough to make seven tests pass for an unrelated reason.
    assert extras_lock.extras_for_module("pandas") == ("clocks", "compose", "hla", "predict")


# --------------------------------------------------------------------------------------------------
# the gate


def test_no_test_module_needs_an_undeclared_package():
    """Every test module in this suite, against what `uv sync --frozen` installs.

    A failure here is not a licence to add a dependency so the test passes. It is one of three things:
    the package belongs in the dev group (then declare it there, with its licence read from the
    installed metadata), or the test genuinely needs an extra (then mark it `requires_extra`), or the
    import belongs inside the function that needs it.
    """
    violations = sorted(str(violation) for violation in extras_lock.scan(sorted(TESTS.glob("*.py"))))
    assert violations == [], "\n".join(
        ["tests that need a package the project does not declare:", *violations]
    )


def test_the_bare_environment_is_what_ci_builds():
    """The derivation and the `test-bare` job must mean the same thing by "bare", or neither holds."""
    workflow = (ROOT / ".github" / "workflows" / "ci.yml").read_text()
    assert "uv sync --frozen\n" in workflow, "test-bare must install from the lock, with no extra"
    assert "--extra" not in workflow.split("test-bare:")[1], "an extra in test-bare defeats the job"


# --------------------------------------------------------------------------------------------------
# the marker itself


def test_a_marker_naming_an_extra_pyproject_does_not_provide_raises():
    """A typo in an extra name must not read as "skipped", which is how a test disappears quietly."""
    with pytest.raises(LookupError):
        extras_lock.missing_modules_for_extra("compsoe")


def test_the_marker_is_registered_so_pytest_does_not_warn_about_it(tmp_path):
    """In a subprocess, with the unknown-mark warning made an error: this session has it registered."""
    marked = tmp_path / "test_marked.py"
    marked.write_text(
        'import pytest\n\n\n@pytest.mark.requires_extra("compose")\ndef test_nothing():\n    pass\n'
    )
    done = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "-rs",
            "-W",
            "error::pytest.PytestUnknownMarkWarning",
            "-p",
            "no:cacheprovider",
            "-c",
            str(ROOT / "pyproject.toml"),
            "--rootdir",
            str(ROOT),
            str(marked),
        ],
        capture_output=True,
        text=True,
        cwd=ROOT,
        check=False,
    )
    assert "PytestUnknownMarkWarning" not in done.stdout + done.stderr, done.stdout
    assert done.returncode == 0, done.stdout + done.stderr


@pytest.mark.requires_extra("compose")
def test_the_compose_runtime_reports_itself_available_when_the_extra_is_installed():
    """A test that genuinely needs an extra, declaring it: in a bare environment this skips by marker.

    What it checks is worth checking on its own. `genomeos/runtime/compose.py` decides AVAILABLE by a
    guarded import, and every test in tests/test_compose.py is skipped unless that flag is True, so a
    flag stuck False would retire the whole file silently -- the same failure shape as the one this file
    is about. Here the extra is installed by declaration, so the flag must be True.
    """
    from genomeos.runtime import compose

    assert compose.AVAILABLE is True
    assert extras_lock.missing_modules_for_extra("compose") == ()
