# SPDX-License-Identifier: AGPL-3.0-or-later
"""`scripts/check.sh` must not lint a non-Python argument as Python, and must not die over one.

WHAT THIS PLANTS. On 2026-10-02 a lane ran `scripts/check.sh ... data/results/manifest_headlines.json`.
ruff read the result file as Python and reported 96 errors in it, `set -e` stopped the run at the
ruff-check leg with `counts: null` -- no test ran at all -- and a RED status file was written for the
tree. Twice: `status-b17ef4d0b2800f21f64e34c4dbee62ea707cc180` (scope_files named
constrained_unknown_targets_v2.json and manifest_headlines.json) and
`status-f8a9cf884bdf2f665f9dde009400d718749dc560` (manifest_headlines.json). Two false reds from
nothing but an argument, in a record other sessions read as real.

AND THE ERROR CLASS, which is the other half. Until that day a red file could not say whether the tree
was at fault. Five reds that day were tooling and environment faults and each read exactly like a test
failure. The last two tests here are the demonstration asked for: a genuine test failure and a tooling
refusal, told apart FROM THE STATUS FILE ALONE, with no log opened and no exit code consulted.

HOW IT IS RUN. The real `scripts/check.sh`, copied into a throwaway git repository, with `uv` and
`shellcheck` stubbed on PATH so the legs cost milliseconds and so the stub can RECORD WHAT RUFF WAS
GIVEN. Only `genomeos.verdict` is real, because the status file is the thing under test. A stub is
used rather than the project's own suite for the reason `tests/test_verdict.py` gives for its own:
the shape is the real one and the cost is not.
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import pytest

from genomeos import verdict

CHECK_SH = Path("scripts/check.sh").resolve()

#: scripts/check.sh SOURCES this to take the machine-wide suite lock, so a harness that copies the
#: one script and not the other gets `No such file or directory` and `set -e` kills the run before
#: any leg: all twelve tests below failed that way the moment the lock was added. Copied, not
#: stubbed, so these runs queue by the real rule -- with a TMPDIR of their own, see `harness`.
SUITE_LOCK_SH = Path("scripts/suite_lock.sh").resolve()
REPO_ROOT = Path.cwd().resolve()

#: The exact argument that produced the two false reds, used verbatim so the planted case is the
#: measured one rather than a stand-in.
THE_JSON = "data/results/manifest_headlines.json"

#: What the stub `uv` writes one line per invocation into, so a test can assert on what ruff was
#: handed and not merely on the verdict that came out.
CALLS = "uv_calls.txt"

_UV_STUB = r"""#!/usr/bin/env bash
# Stands in for `uv` so the legs cost nothing. Every invocation is recorded with its arguments, which
# is the point: the question "was the .json handed to ruff" is answered from this file.
set -u
printf '%s\n' "$*" >> "$RECORD"
shift || true                      # drop `run`
while [ "$#" -gt 0 ]; do           # drop uv's own flags before the tool name
  case "$1" in --quiet | --no-project | --with | --with=*) [ "$1" = "--with" ] && shift; shift ;;
  --) shift; break ;;
  *) break ;;
  esac
done
tool=${1:-}
case "$tool" in
  ruff) exit "${STUB_RUFF_EXIT:-0}" ;;
  pytest) printf '%s\n' "${STUB_PYTEST_SUMMARY:-3 passed in 0.10s}"; exit "${STUB_PYTEST_EXIT:-0}" ;;
  shellcheck) exit 0 ;;
  bio) exit "${STUB_BIO_EXIT:-0}" ;;
  *) exit 0 ;;
esac
"""

_SHELLCHECK_STUB = "#!/usr/bin/env bash\nexit 0\n"


def _git(cwd: Path, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *args], cwd=str(cwd), capture_output=True, text=True)


@pytest.fixture
def harness(tmp_path: Path) -> dict:
    """A throwaway repo holding the real check.sh, with stubs on PATH and its own status directory."""
    repo = tmp_path / "repo"
    (repo / "scripts").mkdir(parents=True)
    for src in (CHECK_SH, SUITE_LOCK_SH):
        (repo / "scripts" / src.name).write_bytes(src.read_bytes())
        (repo / "scripts" / src.name).chmod(0o755)
    (repo / "data" / "results").mkdir(parents=True)
    (repo / "data" / "results" / "manifest_headlines.json").write_text('{"ok": true}\n')
    (repo / "keep.py").write_text("x = 1\n")
    _git(repo, "init", "-q", "-b", "main", ".")

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "uv").write_text(_UV_STUB)
    (bin_dir / "uv").chmod(0o755)
    (bin_dir / "shellcheck").write_text(_SHELLCHECK_STUB)
    (bin_dir / "shellcheck").chmod(0o755)

    status_dir = tmp_path / "status"
    record = tmp_path / CALLS
    env = dict(os.environ)
    env.update(
        PATH=f"{bin_dir}:{env['PATH']}",
        PYTHONPATH=str(REPO_ROOT),
        GENOMEOS_CHECK_STATUS_DIR=str(status_dir),
        RECORD=str(record),
        # Since 2026-10-02 scripts/check.sh takes "$TMPDIR/genomeos-suite.lock" so that one full
        # suite runs on the machine at a time (scripts/suite_lock.sh). A TMPDIR of its own is what
        # keeps these runs off the REAL lock: inheriting it, this file would block for twelve
        # minutes behind a peer's check and the lock would be a test dependency instead of a lock.
        # tests/test_push_lock.py isolates TMPDIR for the push lock for the same reason.
        TMPDIR=str(tmp_path / "tmp"),
    )
    (tmp_path / "tmp").mkdir()
    return {"repo": repo, "env": env, "status_dir": status_dir, "record": record}


def _run(harness: dict, *args: str, **stub_env: str) -> subprocess.CompletedProcess:
    env = dict(harness["env"])
    env.update(stub_env)
    return subprocess.run(
        ["bash", "scripts/check.sh", *args],
        cwd=str(harness["repo"]),
        env=env,
        capture_output=True,
        text=True,
    )


def _verdict_file(harness: dict) -> dict:
    """THE verdict, read the way the project requires: from the artefact, never from an exit code."""
    latest = harness["status_dir"] / "latest.json"
    assert latest.is_file(), f"no verdict file was written; {sorted(harness['status_dir'].glob('*'))}"
    return json.loads(latest.read_text())


def _ruff_calls(harness: dict) -> list[str]:
    if not harness["record"].is_file():
        return []
    return [line for line in harness["record"].read_text().splitlines() if " ruff " in f" {line} "]


# ---- item 2: a non-Python argument is refused by name and does not kill the run ----------------


def test_a_json_argument_is_never_handed_to_ruff(harness: dict) -> None:
    """The first half of the false red: 96 ruff errors in a result file, because ruff was given it."""
    r = _run(harness, THE_JSON, "keep.py")
    calls = _ruff_calls(harness)
    assert calls, "ruff was not run at all; the test is not exercising the lint legs"
    for call in calls:
        assert THE_JSON not in call, f"ruff was handed the .json: {call}"
    assert any("keep.py" in call for call in calls), calls
    assert r.returncode == 0, r.stderr[-800:]


def test_a_json_argument_produces_a_named_refusal_and_not_a_red_verdict(harness: dict) -> None:
    """The second half: the run died at ruff-check with counts null. It must now be green and say why."""
    _run(harness, THE_JSON)
    status = _verdict_file(harness)
    assert status["verdict"] == "green", status
    assert status["error_class"] == "", status
    assert THE_JSON in status["note"], status["note"]
    assert status["counts"] is not None, "the tests must have run; a refused argument may not skip them"


def test_the_refusal_names_the_argument_on_stderr_and_offers_the_right_tool(harness: dict) -> None:
    r = _run(harness, THE_JSON)
    assert "REFUSED these lint arguments by name" in r.stderr
    assert THE_JSON in r.stderr
    assert "json.load" in r.stderr, "a refusal that does not say what WOULD check a .json is half a refusal"


def test_an_argument_list_with_no_python_in_it_does_not_become_a_project_lint(harness: dict) -> None:
    """`ruff check` with no paths lints everything. A file-scoped run must skip the leg, not widen."""
    _run(harness, THE_JSON)
    assert _ruff_calls(harness) == [], _ruff_calls(harness)
    status = _verdict_file(harness)
    assert status["scope"] == "files", status
    assert status["scope_files"] == [THE_JSON], status


def test_a_real_lint_failure_on_a_python_file_is_still_red(harness: dict) -> None:
    """Not weakened: refusing a .json must not make ruff's own verdict on a .py file optional."""
    _run(harness, "keep.py", STUB_RUFF_EXIT="1")
    status = _verdict_file(harness)
    assert status["verdict"] == "red", status
    assert status["failed_leg"] == "ruff-check", status


# ---- item 3: a tooling red and a test red, told apart from the status file alone ----------------


def test_a_genuine_test_failure_is_classified_test(harness: dict) -> None:
    _run(harness, "keep.py", STUB_PYTEST_EXIT="1", STUB_PYTEST_SUMMARY="1 failed, 2 passed in 0.10s")
    status = _verdict_file(harness)
    assert (status["verdict"], status["error_class"], status["failed_leg"]) == ("red", "test", "pytest")
    assert status["counts"]["failed"] == 1


def test_a_tooling_refusal_is_classified_tooling(harness: dict) -> None:
    _run(harness, "keep.py", STUB_RUFF_EXIT="1")
    status = _verdict_file(harness)
    assert (status["verdict"], status["error_class"], status["failed_leg"]) == (
        "red",
        "tooling",
        "ruff-check",
    )


def test_the_two_reds_are_distinguishable_from_the_status_file_alone(harness: dict) -> None:
    """The demonstration asked for: one field, read without opening a log or trusting an exit code."""
    _run(harness, "keep.py", STUB_PYTEST_EXIT="1", STUB_PYTEST_SUMMARY="1 failed, 2 passed in 0.10s")
    test_red = _verdict_file(harness)
    _run(harness, "keep.py", STUB_RUFF_EXIT="1")
    tooling_red = _verdict_file(harness)

    assert test_red["verdict"] == tooling_red["verdict"] == "red"
    assert test_red["exit_code"] == tooling_red["exit_code"] == 1
    assert test_red["reason"] == tooling_red["reason"]  # "the check exited 1": identical, as before
    assert test_red["error_class"] != tooling_red["error_class"]
    assert {test_red["error_class"], tooling_red["error_class"]} == {"test", "tooling"}


def test_twenty_setup_errors_and_no_failures_are_test_setup_and_not_test(harness: dict) -> None:
    """The symlink false red, by its counts: 20 errors, 0 failures. Nothing was asserted wrongly."""
    _run(
        harness,
        "keep.py",
        STUB_PYTEST_EXIT="1",
        STUB_PYTEST_SUMMARY="4436 passed, 11 skipped, 7 xfailed, 20 errors in 900.00s",
    )
    status = _verdict_file(harness)
    assert status["error_class"] == "test_setup", status
    assert (status["counts"]["errors"], status["counts"]["failed"]) == (20, 0)


def test_a_green_run_carries_no_error_class(harness: dict) -> None:
    _run(harness, "keep.py")
    status = _verdict_file(harness)
    assert (status["verdict"], status["error_class"]) == ("green", "")


def test_a_tree_that_moved_is_classified_tree_moved_whatever_the_counts_say() -> None:
    assert verdict.classify("tree_moved", "pytest", {"failed": 3, "passed": 1}) == "tree_moved"


def test_every_class_written_is_one_of_the_classes_the_file_lists(harness: dict) -> None:
    """The list travels with the value, so a reader needs neither this test nor the docstring."""
    _run(harness, "keep.py", STUB_RUFF_EXIT="1")
    status = _verdict_file(harness)
    assert status["error_classes"] == list(verdict.ERROR_CLASSES)
    assert status["error_class"] in status["error_classes"]


def test_a_red_of_any_class_is_still_refused_for_a_push(harness: dict) -> None:
    """A reading aid, not an excuse: `require` must not start accepting a merely-tooling red."""
    _run(harness, "keep.py", STUB_RUFF_EXIT="1")
    status = _verdict_file(harness)
    code = subprocess.run(
        ["python3", "-m", "genomeos.verdict", "require", "--tree", status["tree_begin"]],
        cwd=str(harness["repo"]),
        env=harness["env"],
        capture_output=True,
        text=True,
    )
    assert code.returncode == verdict.EXIT_REFUSED, code
    assert "error_class 'tooling'" in code.stderr, code.stderr
    assert "NOT a test failure" in code.stderr, code.stderr
