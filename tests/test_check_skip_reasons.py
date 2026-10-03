# SPDX-License-Identifier: AGPL-3.0-or-later
"""`scripts/check.sh`'s log must name every skip BY TEST ID, so a verdict can be read without a re-run.

THE DEFECT THIS PLANTS AGAINST. Until 2026-10-03 the pytest leg ran bare `pytest -q`. A green status
file therefore recorded a skip COUNT and nothing else: push19's verdict for tree c78041b1 said
"5,053 passed, 12 skipped, 8 xfailed" and could not show, for any named test, that it RAN. Every
requirement of the form "the signed tree's verdict must show test X PASSED, not skipped" then needs
somebody to re-run those tests by hand, which is what happened that night -- a separate 260-test run
to establish what the kept file could not.

AND WHY `-rs` ALONE DOES NOT FIX IT, which is the finding worth keeping. pytest FOLDS its skip
summary by (file, line, reason): `-rs` prints `SKIPPED [3] tests/test_x.py:12: reason`, a count and a
location, no test id -- and three parametrised cases of one test collapse into that single line.
`--no-fold-skipped` prints one line per test, `SKIPPED tests/test_x.py::test_y - Skipped: reason`.
Both flags are needed and this file pins both, by consequence rather than by spelling.

HOW IT IS RUN, AND WHY IT IS NOT A GREP FOR `-rs`. A test that searched scripts/check.sh for the
string `-rs` would prove the flag is WRITTEN, not that the log carries skip ids; this project has a
lesson from the same night about exactly that distinction (0cfe1e4: a guard that fired only because a
test patched the source to raise proved the branch, not its reachability). So this PLANTS SKIPS and
reads the real artefact: a throwaway git repository holds the real scripts/check.sh and
scripts/suite_lock.sh, `ruff`, `shellcheck` and `bio` are stubbed to cost milliseconds, and THE
PYTEST LEG IS A REAL PYTEST over planted tests that skip with a reason chosen here. The assertions
then read the log that the status file itself names, and look for the planted ids in it.

NO LOCAL DATA. Everything is built under tmp_path, so nothing here needs data/cache, data/knowledge
or data/reference and nothing here can skip for want of them -- which matters, since a planted skip
that became a CI failure would be this file failing at its own subject.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

CHECK_SH = Path("scripts/check.sh").resolve()
SUITE_LOCK_SH = Path("scripts/suite_lock.sh").resolve()
REPO_ROOT = Path.cwd().resolve()

#: The reason the planted tests skip with. Long and unmistakable on purpose: a reader of the log must
#: be able to tell this apart from any real skip, and a substring search for it must not match by
#: accident anywhere else in a suite's output.
PLANTED_REASON = "planted-skip-reason-for-check-sh-log-by-test-id"

#: The planted file, as the log will spell it (relative to the harness repo, which is pytest's
#: rootdir there because of the pyproject.toml the harness writes).
PLANTED_FILE = "tests/test_planted_skips.py"

#: One test that skips on its own, and one parametrised test whose three cases SHARE a file, a line
#: and a reason. The parametrised three are the discriminating plant: folded they are one line
#: reading `SKIPPED [3] ...`, unfolded they are three lines each naming its own id.
_PLANTED_TESTS = f'''import pytest

REASON = "{PLANTED_REASON}"


@pytest.mark.skip(reason=REASON)
def test_planted_skip_alone():
    raise AssertionError("this body must never run")


@pytest.mark.parametrize("n", [1, 2, 3])
@pytest.mark.skip(reason=REASON)
def test_planted_skip_folded(n):
    raise AssertionError("this body must never run")


def test_planted_pass():
    assert True
'''

#: The ids the log must carry, one line each.
PLANTED_IDS = (
    f"{PLANTED_FILE}::test_planted_skip_alone",
    f"{PLANTED_FILE}::test_planted_skip_folded[1]",
    f"{PLANTED_FILE}::test_planted_skip_folded[2]",
    f"{PLANTED_FILE}::test_planted_skip_folded[3]",
)

#: What the stub `uv` writes one line per invocation into, so a test can read what each tool was
#: handed rather than only the verdict that came out.
CALLS = "uv_calls.txt"

#: `-p no:cacheprovider` is added by the HARNESS and not by check.sh: the inner pytest would
#: otherwise write .pytest_cache into the harness repo mid-run, the working tree would differ between
#: tree_begin and tree_end, and every verdict here would be `tree_moved` instead of green. The
#: harness repo also carries a .gitignore for the same reason, since working_tree_hash respects it.
_UV_STUB = r"""#!/usr/bin/env bash
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
shift || true
case "$tool" in
  ruff | shellcheck | bio) exit 0 ;;
  pytest)
    if [ "${STUB_PYTEST_KNOWS_NO_FOLD:-1}" != "1" ]; then
      # Stands in for a pytest older than 8.3: the flag is absent from --help, and passing it is a
      # USAGE ERROR with no test run -- which is how an unprobed flag would turn a green tree red.
      for a in "$@"; do
        if [ "$a" = "--help" ]; then
          "$REAL_PYTHON" -m pytest --help | sed '/--no-fold-skipped/d'
          exit 0
        fi
      done
      for a in "$@"; do
        if [ "$a" = "--no-fold-skipped" ]; then
          echo "ERROR: unrecognized arguments: --no-fold-skipped" >&2
          exit 4
        fi
      done
    fi
    exec "$REAL_PYTHON" -m pytest -p no:cacheprovider "$@"
    ;;
  *) exit 0 ;;
esac
"""

_SHELLCHECK_STUB = "#!/usr/bin/env bash\nexit 0\n"


def _git(cwd: Path, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *args], cwd=str(cwd), capture_output=True, text=True)


def _build(tmp_path: Path) -> dict:
    repo = tmp_path / "repo"
    (repo / "scripts").mkdir(parents=True)
    for src in (CHECK_SH, SUITE_LOCK_SH):
        (repo / "scripts" / src.name).write_bytes(src.read_bytes())
        (repo / "scripts" / src.name).chmod(0o755)
    # testpaths, so that the bare `uv run pytest -q ...` check.sh issues collects the planted file
    # the same way it collects this project's own tests, and so pytest's rootdir is this repo and the
    # ids in the log are spelled `tests/test_planted_skips.py::...`.
    (repo / "pyproject.toml").write_text('[tool.pytest.ini_options]\ntestpaths = ["tests"]\n')
    (repo / ".gitignore").write_text("__pycache__/\n.pytest_cache/\n")
    (repo / "tests").mkdir()
    (repo / "tests" / "test_planted_skips.py").write_text(_PLANTED_TESTS)
    (repo / "keep.py").write_text("x = 1\n")
    _git(repo, "init", "-q", "-b", "main", ".")

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "uv").write_text(_UV_STUB)
    (bin_dir / "uv").chmod(0o755)
    (bin_dir / "shellcheck").write_text(_SHELLCHECK_STUB)
    (bin_dir / "shellcheck").chmod(0o755)

    status_dir = tmp_path / "status"
    env = dict(os.environ)
    env.pop("PYTEST_CURRENT_TEST", None)
    env.pop("PYTEST_ADDOPTS", None)
    env.update(
        PATH=f"{bin_dir}:{env['PATH']}",
        PYTHONPATH=str(REPO_ROOT),
        GENOMEOS_CHECK_STATUS_DIR=str(status_dir),
        RECORD=str(tmp_path / CALLS),
        REAL_PYTHON=sys.executable,
        # A TMPDIR of its own, so these runs queue on a lock of their own and never behind a peer's
        # twelve-minute suite on the machine's real one (tests/test_check_sh_arguments.py, same note).
        TMPDIR=str(tmp_path / "tmp"),
    )
    (tmp_path / "tmp").mkdir()
    return {"repo": repo, "env": env, "status_dir": status_dir, "record": tmp_path / CALLS}


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


def _verdict(harness: dict) -> dict:
    latest = harness["status_dir"] / "latest.json"
    assert latest.is_file(), f"no verdict file was written; {sorted(harness['status_dir'].glob('*'))}"
    return json.loads(latest.read_text())


def _log_from_the_verdict(status: dict) -> str:
    """The log THE STATUS FILE NAMES, not one this test went looking for.

    That is the whole point of the lane: a reader holding a signed verdict must be able to reach the
    skip ids from the verdict alone, with no knowledge of how the run was invoked.
    """
    assert status["pytest_log"], f"the verdict names no pytest log: {status}"
    log = Path(status["pytest_log"])
    assert log.is_file(), f"the verdict names {log}, which does not exist"
    return log.read_text(errors="replace")


@pytest.fixture(scope="module")
def ran(tmp_path_factory: pytest.TempPathFactory) -> dict:
    """One real check.sh run over the planted skips, shared by the assertions below."""
    harness = _build(tmp_path_factory.mktemp("skipids"))
    proc = _run(harness)
    status = _verdict(harness)
    return {"harness": harness, "proc": proc, "status": status, "log": _log_from_the_verdict(status)}


@pytest.fixture(scope="module")
def ran_without_the_flag(tmp_path_factory: pytest.TempPathFactory) -> dict:
    """The same run against a pytest that does not know --no-fold-skipped."""
    harness = _build(tmp_path_factory.mktemp("skipids-oldpytest"))
    proc = _run(harness, STUB_PYTEST_KNOWS_NO_FOLD="0")
    status = _verdict(harness)
    return {"harness": harness, "proc": proc, "status": status, "log": _log_from_the_verdict(status)}


# ---- the planted skip, named by id, with its reason ---------------------------------------------


def test_the_log_names_the_planted_skip_by_test_id_and_quotes_its_reason(ran: dict) -> None:
    """The lane's requirement, read out of the artefact a verdict points at."""
    lines = [ln for ln in ran["log"].splitlines() if PLANTED_IDS[0] in ln]
    assert len(lines) == 1, f"expected one line naming {PLANTED_IDS[0]}; got {lines}"
    line = lines[0]
    assert line.startswith("SKIPPED"), line
    assert PLANTED_REASON in line, f"the line names the test but not its reason: {line}"


def test_every_planted_skip_has_a_line_of_its_own_so_folding_cannot_hide_one(ran: dict) -> None:
    """The three parametrised cases share a file, a line and a reason: folded they are ONE line."""
    log = ran["log"]
    for test_id in PLANTED_IDS:
        found = [ln for ln in log.splitlines() if test_id in ln and ln.startswith("SKIPPED")]
        assert len(found) == 1, f"expected exactly one SKIPPED line for {test_id}; got {found}"
    folded = [ln for ln in log.splitlines() if ln.startswith("SKIPPED [")]
    assert folded == [], f"skips are still folded, so an id is unrecoverable: {folded}"


def test_a_test_that_ran_can_be_told_from_one_that_skipped_without_a_re_run(ran: dict) -> None:
    """The consequence Albert asked for: no named skip line means the test was not skipped."""
    log = ran["log"]
    skipped_ids = {
        ln.split(" - ")[0].removeprefix("SKIPPED").strip()
        for ln in log.splitlines()
        if ln.startswith("SKIPPED")
    }
    assert set(PLANTED_IDS) <= skipped_ids, skipped_ids
    assert f"{PLANTED_FILE}::test_planted_pass" not in skipped_ids, skipped_ids


def test_the_skip_ids_are_reachable_from_the_verdict_alone(ran: dict) -> None:
    """`pytest_log` is an existing file, inside the status directory, holding the ids."""
    status = ran["status"]
    log = Path(status["pytest_log"])
    assert log.is_file(), status["pytest_log"]
    assert log.parent.resolve() == Path(status["pytest_log"]).parent.resolve()
    assert log.parent.name == ran["harness"]["status_dir"].name, log
    assert PLANTED_IDS[0] in log.read_text(), "the log the verdict names does not carry the ids"


# ---- the defect itself, measured rather than asserted -------------------------------------------


def test_without_the_flags_the_same_planted_skip_is_unfindable_by_id(
    tmp_path: Path,
) -> None:
    """What push19's log could not answer, reproduced: bare `-q`, and `-rs` without unfolding.

    Not a test of check.sh -- a measurement of the two outputs this lane replaced, so the plant above
    is known to discriminate and is not passing for some other reason.
    """
    harness = _build(tmp_path)
    repo = harness["repo"]

    def real_pytest(*flags: str) -> str:
        return subprocess.run(
            [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", *flags],
            cwd=str(repo),
            env=harness["env"],
            capture_output=True,
            text=True,
        ).stdout

    bare = real_pytest()
    assert "4 skipped" in bare, bare[-400:]
    assert PLANTED_IDS[0] not in bare, "bare -q already names ids; the premise of this lane is wrong"
    assert PLANTED_REASON not in bare, "bare -q already names reasons"

    folded = real_pytest("-rs")
    assert PLANTED_REASON in folded, folded[-400:]
    assert PLANTED_IDS[0] not in folded, "-rs alone names the test id; --no-fold-skipped is then moot"
    assert "SKIPPED [3]" in folded, f"-rs did not fold the parametrised three: {folded[-400:]}"


# ---- nothing else about the verdict moves -------------------------------------------------------


def test_the_run_is_green_and_carries_no_error_class(ran: dict) -> None:
    status = ran["status"]
    assert (status["verdict"], status["error_class"]) == ("green", ""), status
    assert status["exit_code"] == 0, status
    assert ran["proc"].returncode == 0, ran["proc"].stderr[-800:]


def test_the_counts_are_still_read_from_the_summary_line_unchanged(ran: dict) -> None:
    """The flags add lines BEFORE the summary; the line the counts come from is untouched."""
    status = ran["status"]
    assert status["counts_read"] is True, status
    assert status["counts"]["passed"] == 1, status["counts"]
    assert status["counts"]["skipped"] == 4, status["counts"]
    assert status["counts"]["failed"] == 0 and status["counts"]["errors"] == 0, status["counts"]


def test_the_note_is_empty_when_the_ids_were_recorded(ran: dict) -> None:
    """A note is for what a run carries LESS of than its reader thinks. This run carries them all."""
    assert "FOLDED" not in ran["status"]["note"], ran["status"]["note"]


def test_a_real_test_failure_is_still_red_and_still_classified_test(tmp_path: Path) -> None:
    """Not weakened: naming skips must not make a failure any softer."""
    harness = _build(tmp_path)
    (harness["repo"] / "tests" / "test_planted_fail.py").write_text(
        "def test_planted_failure():\n    assert False\n"
    )
    _run(harness)
    status = _verdict(harness)
    assert (status["verdict"], status["error_class"], status["failed_leg"]) == ("red", "test", "pytest")
    assert status["counts"]["failed"] == 1, status["counts"]


# ---- an older pytest degrades loudly and is never a red ------------------------------------------


def test_a_pytest_without_the_flag_is_green_and_not_a_usage_error(ran_without_the_flag: dict) -> None:
    """The flag is probed, not assumed: an unknown flag would be exit 4 with no test run at all."""
    status = ran_without_the_flag["status"]
    assert (status["verdict"], status["error_class"]) == ("green", ""), status
    assert status["exit_code"] == 0, status
    assert status["counts"]["skipped"] == 4, status["counts"]


def test_a_pytest_without_the_flag_says_what_the_log_then_lacks(ran_without_the_flag: dict) -> None:
    """Loudly on stderr AND in the verdict, the shellcheck leg's rule."""
    status = ran_without_the_flag["status"]
    assert "FOLDED" in status["note"], status["note"]
    assert "test id" in status["note"], status["note"]
    assert "--no-fold-skipped" in ran_without_the_flag["proc"].stderr


def test_a_pytest_without_the_flag_still_records_every_reason(ran_without_the_flag: dict) -> None:
    """Degraded is not nothing: `-rs` alone still names each skip's file, line and reason."""
    log = ran_without_the_flag["log"]
    assert PLANTED_REASON in log, log[-600:]
    assert "SKIPPED [3]" in log, f"the folded form is what this case is about: {log[-600:]}"
