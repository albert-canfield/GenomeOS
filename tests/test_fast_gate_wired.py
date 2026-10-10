# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The fast gate as WIRED into `scripts/pre-push.sh` on 2026-10-10 (ROADMAP section 5, E1).

THE RULING (the supervisor, 2026-10-03; ROADMAP-HISTORY beside `1b2e768`; CONTRIBUTING step 4). A push
to dev takes the FAST gate -- `ruff check`, `ruff format --check`, ALWAYS_RUN and the test files changed
since the last green CI, in both legs, normal and store-free -- and its status file is scope "fast" and
is never read as a project verdict. The FULL suite stays where it decides something: when the daily CI
is red or unread, on the promotion to main, and on any paid path. The stated trade: a red in a test
file that did not change does not block a dev push.

HOW THESE TESTS ARE BUILT. The plants run THE REAL `scripts/pre-push.sh` end to end, in a throwaway
repository that carries real copies of the scripts it calls (`fast_prepush.py`, `last_green_ci.py`,
`link_stores_read_only.py`, `remove_worktree.py`, `push_own.sh`, `genomeos/verdict.py`), with two
stand-ins and nothing else: `uv`, which only drops `run` and execs this interpreter, and
`scripts/check.sh`, which prints a marker and writes a green project verdict -- the full suite itself is
not what is under test here, WHICH GATE RUNS is. So a test that says "the full gate was taken" sees
the marker, and a test that says "the fast gate blocked" sees the fast leg go red and the hook exit 1.
The ruling's four plants are all here: a format-only defect refuses, an ALWAYS_RUN failure refuses, a
red in a CHANGED test file refuses, and a red in an UNCHANGED file does not block.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from datetime import UTC, datetime, timedelta
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import fast_prepush as fp  # noqa: E402
import last_green_ci as lg  # noqa: E402

from genomeos import verdict  # noqa: E402

COPIED = (
    "scripts/pre-push.sh",
    "scripts/fast_prepush.py",
    "scripts/last_green_ci.py",
    "scripts/link_stores_read_only.py",
    "scripts/remove_worktree.py",
    "scripts/push_own.sh",
    "genomeos/verdict.py",
)

FULL_MARKER = "STUB-FULL-CHECK-RAN"
CHECK_STUB = f"""#!/usr/bin/env bash
# stand-in for the full suite: says it ran, and writes a green PROJECT verdict for the tree it is in
set -euo pipefail
cd "$(dirname "$0")/.."
echo "{FULL_MARKER} in $(pwd)"
tree=$(python3 -m genomeos.verdict tree)
log=$(mktemp "${{TMPDIR:-/tmp}}/stub-check.XXXXXX")
echo "1 passed in 0.01s" > "$log"
python3 -m genomeos.verdict write --exit-code 0 --tree-begin "$tree" --scope project \\
  --started-at "$(date +%s)" --pytest-log "$log"
"""

PASS = "def test_ok():\n    assert True\n"
FAIL = 'def test_red():\n    assert False, "planted red"\n'
LINT = "import os\n\n\ndef test_ok():\n    assert True\n"  # F401, a lint defect only
FORMAT = "def test_ok():\n    assert  True\n"  # two spaces: a format defect only, lint clean

ALWAYS_RUN = (
    'MEMBERS = ({"path": "tests/test_always.py", "category": "plant", '
    '"why": "the planted ALWAYS_RUN member"},)\n'
)


def _git(repo: Path, *args: str) -> str:
    out = subprocess.run(["git", *args], cwd=repo, capture_output=True, text=True, check=True, env=_env(repo))
    return out.stdout.strip()


def _env(repo: Path, **extra: str) -> dict[str, str]:
    """A clean environment: no GIT_* from whatever ran this suite, no user git config, own TMPDIR."""
    home = repo.parent / "home"
    tmp = repo.parent / "tmpdir"
    home.mkdir(exist_ok=True)
    tmp.mkdir(exist_ok=True)
    bindir = repo.parent / "bin"
    env = {
        "PATH": os.pathsep.join(
            [str(bindir), str(Path(sys.executable).parent), "/usr/bin", "/bin", "/usr/sbin", "/sbin"]
        ),
        "HOME": str(home),
        "TMPDIR": str(tmp),
        "GIT_CONFIG_NOSYSTEM": "1",
        "LANG": "C.UTF-8",
    }
    env.update(extra)
    return env


def _write(repo: Path, rel: str, body: str, mode: int | None = None) -> None:
    p = repo / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(body)
    if mode is not None:
        p.chmod(mode)


def _commit(repo: Path, files: dict[str, str], message: str) -> str:
    for rel, body in files.items():
        _write(repo, rel, body)
        _git(repo, "add", "--", rel)
    _git(repo, "commit", "-q", "-m", message)
    return _git(repo, "rev-parse", "HEAD")


def build(tmp_path: Path, always: str = PASS) -> dict:
    """A repository whose green-CI sha carries a FAILING test file that nothing changes afterwards."""
    repo = tmp_path / "repo"
    repo.mkdir()
    bindir = tmp_path / "bin"
    bindir.mkdir()
    # `uv run [--flags] python ARGS` -> this interpreter; anything else is executed as named
    _write(
        bindir,
        "uv",
        "#!/bin/sh\n"
        '[ "$1" = run ] && shift\n'
        'while [ $# -gt 0 ]; do case "$1" in --*) shift ;; *) break ;; esac; done\n'
        f'case "$1" in python|python3) shift; exec "{sys.executable}" "$@" ;; esac\n'
        'exec "$@"\n',
        0o755,
    )
    _git(repo, "init", "-q", "-b", "dev")
    _git(repo, "config", "user.email", "planted@example.invalid")
    _git(repo, "config", "user.name", "planted")
    _git(repo, "config", "commit.gpgsign", "false")
    base = {rel: (ROOT / rel).read_text() for rel in COPIED}
    base.update(
        {
            "genomeos/__init__.py": '"""the planted package: the import-path assertions need one here."""\n',
            "scripts/check.sh": CHECK_STUB,
            ".gitignore": "data/knowledge\ndata/cache\ndata/reference\n__pycache__/\n.pytest_cache/\n",
            # the copied scripts follow the project's ruff settings, which are not copied, so ruff here
            # reads the planted tests only: the defects under test are planted there
            "ruff.toml": 'line-length = 110\nextend-exclude = ["scripts", "genomeos"]\n',
            "tests/always_run.py": ALWAYS_RUN,
            "tests/test_always.py": always,
            "tests/test_unchanged.py": FAIL,
        }
    )
    for rel, body in base.items():
        _write(repo, rel, body, 0o755 if rel.endswith(".sh") else None)
        _git(repo, "add", "--", rel)
    _git(repo, "commit", "-q", "-m", "the green CI sha")
    green = _git(repo, "rev-parse", "HEAD")
    _git(repo, "update-ref", "refs/remotes/origin/dev", green)
    return {"repo": repo, "green": green}


def read_daily(h: dict, conclusion: str = "success", hours_ago: float = 1.0) -> None:
    """The coordinator's morning reading of CI, recorded the way the procedure records it."""
    when = datetime.now(UTC) - timedelta(hours=hours_ago)
    lg.write_daily(h["repo"], conclusion, h["green"], "18200000000", "2026-10-10T06:00:00Z", read_at=when)


def push(h: dict, sha: str, remote_ref: str = "refs/heads/dev", **env: str) -> subprocess.CompletedProcess:
    repo = h["repo"]
    return subprocess.run(
        ["bash", "scripts/pre-push.sh", "origin", "https://example.invalid/x.git"],
        cwd=repo,
        input=f"refs/heads/dev {sha} {remote_ref} {'0' * 40}\n",
        capture_output=True,
        text=True,
        timeout=600,
        env=_env(repo, **env),
    )


def fast_file(h: dict, sha: str) -> dict:
    tree = _git(h["repo"], "rev-parse", f"{sha}^{{tree}}")
    path = h["repo"] / ".git" / "genomeos-check" / f"fast-{tree}.json"
    assert path.exists(), f"no fast verdict at {path}"
    return json.loads(path.read_text())


def out(r: subprocess.CompletedProcess) -> str:
    return r.stdout + r.stderr


@pytest.fixture
def h(tmp_path: Path) -> dict:
    return build(tmp_path)


# --------------------------------------------------------------------------- the ruling's plants


def test_a_red_in_a_CHANGED_test_file_blocks_the_push(h: dict) -> None:
    read_daily(h)
    sha = _commit(h["repo"], {"tests/test_changed.py": FAIL}, "a changed test file that fails")
    r = push(h, sha)
    assert r.returncode != 0, out(r)
    assert "FAST gate" in r.stdout and FULL_MARKER not in out(r), out(r)
    status = fast_file(h, sha)
    assert status["scope"] == "fast" and status["verdict"] == "red", status
    assert status["failed_leg"] == "normal", status
    assert "tests/test_changed.py" in status["plan"]["selection"]["files"]


def test_a_red_in_an_UNCHANGED_test_file_does_not_block_the_push_THE_STATED_TRADE(h: dict) -> None:
    read_daily(h)
    sha = _commit(h["repo"], {"tests/test_changed.py": PASS}, "a changed test file that passes")
    # the trade is only a trade if the unchanged file really is red: shown, not assumed
    control = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", "tests/test_unchanged.py"],
        cwd=h["repo"],
        capture_output=True,
        text=True,
        env=_env(h["repo"]),
    )
    assert control.returncode == 1, control.stdout
    r = push(h, sha)
    assert r.returncode == 0, out(r)
    assert "FAST gate" in r.stdout and FULL_MARKER not in out(r), out(r)
    assert "GENOMEOS_FAST_OK" in r.stdout, r.stdout
    status = fast_file(h, sha)
    assert status["verdict"] == "green" and status["scope"] == "fast", status
    files = status["plan"]["selection"]["files"]
    assert "tests/test_changed.py" in files and "tests/test_always.py" in files, files
    assert "tests/test_unchanged.py" not in files, files
    legs = {leg["leg"]: leg for leg in status["legs"]}
    assert set(legs) == {"ruff-check", "ruff-format", "normal", "store-free"}, legs
    assert legs["store-free"]["stores"]["ok"] and legs["store-free"]["import_path"]["ok"], legs["store-free"]


def test_an_ALWAYS_RUN_failure_blocks_the_push_although_its_file_did_not_change(tmp_path: Path) -> None:
    h = build(tmp_path, always=FAIL)
    read_daily(h)
    sha = _commit(h["repo"], {"tests/test_changed.py": PASS}, "a passing change")
    r = push(h, sha)
    assert r.returncode != 0, out(r)
    status = fast_file(h, sha)
    assert status["verdict"] == "red" and status["failed_leg"] == "normal", status
    assert "tests/test_always.py" in status["plan"]["selection"]["files"]


def test_a_FORMAT_only_defect_blocks_the_push(h: dict) -> None:
    read_daily(h)
    sha = _commit(h["repo"], {"tests/test_changed.py": FORMAT}, "passes and lints, not formatted")
    r = push(h, sha)
    assert r.returncode != 0, out(r)
    status = fast_file(h, sha)
    assert status["verdict"] == "red" and status["failed_leg"] == "ruff-format", status
    assert [leg["leg"] for leg in status["legs"]] == ["ruff-check", "ruff-format"], status["legs"]


def test_a_LINT_defect_blocks_the_push(h: dict) -> None:
    read_daily(h)
    sha = _commit(h["repo"], {"tests/test_changed.py": LINT}, "passes, an unused import")
    r = push(h, sha)
    assert r.returncode != 0, out(r)
    status = fast_file(h, sha)
    assert status["verdict"] == "red" and status["failed_leg"] == "ruff-check", status


# --------------------------------------------------------------------------- where the full suite stays


@pytest.mark.parametrize(
    "reading",
    [
        pytest.param({"conclusion": "failure"}, id="red"),
        pytest.param({"conclusion": "cancelled"}, id="cancelled-is-not-green"),
        pytest.param(None, id="unread-no-reading"),
        pytest.param({"conclusion": "success", "hours_ago": 30.0}, id="unread-today-stale-reading"),
    ],
)
def test_a_red_or_unread_daily_CI_sends_the_push_to_the_FULL_gate(h: dict, reading: dict | None) -> None:
    if reading is not None:
        read_daily(h, **reading)
    # a changed red test: the fast gate would refuse it, the stand-in full suite does not, so a green
    # push here can only mean the full gate was taken
    sha = _commit(h["repo"], {"tests/test_changed.py": FAIL}, "a red the fast gate would refuse")
    r = push(h, sha)
    assert FULL_MARKER in r.stdout, out(r)
    assert "GATE full" in r.stdout and "FAST gate" not in r.stdout, out(r)
    assert r.returncode == 0, out(r)
    tree = _git(h["repo"], "rev-parse", f"{sha}^{{tree}}")
    assert not (h["repo"] / ".git" / "genomeos-check" / f"fast-{tree}.json").exists()


@pytest.mark.parametrize(
    ("remote_ref", "env"),
    [
        pytest.param("refs/heads/promote-abc1234", {}, id="the-promotion-branch"),
        pytest.param("refs/heads/main", {}, id="main"),
        pytest.param("refs/heads/dev", {"GENOMEOS_FULL_GATE": "1"}, id="the-paid-path-switch"),
    ],
)
def test_the_promotion_and_the_paid_path_take_the_FULL_gate_on_a_green_daily_CI(
    h: dict, remote_ref: str, env: dict
) -> None:
    read_daily(h)
    sha = _commit(h["repo"], {"tests/test_changed.py": PASS}, "a change")
    r = push(h, sha, remote_ref=remote_ref, **env)
    assert FULL_MARKER in r.stdout, out(r)
    assert "FAST gate" not in r.stdout, out(r)
    assert r.returncode == 0, out(r)


def test_a_ref_the_hook_does_not_check_is_still_not_checked(h: dict) -> None:
    read_daily(h)
    sha = _commit(h["repo"], {"tests/test_changed.py": FAIL}, "a feature branch")
    r = push(h, sha, remote_ref="refs/heads/feature")
    assert r.returncode == 0 and FULL_MARKER not in r.stdout and "FAST gate" not in r.stdout, out(r)


# --------------------------------------------------------------------------- never a project verdict


def test_the_fast_verdict_is_NEVER_read_as_a_project_verdict_by_any_reader(h: dict) -> None:
    read_daily(h)
    sha = _commit(h["repo"], {"tests/test_changed.py": PASS}, "a change")
    repo = h["repo"]
    tree = _git(repo, "rev-parse", f"{sha}^{{tree}}")
    r = push(h, sha)
    assert r.returncode == 0 and fast_file(h, sha)["verdict"] == "green", out(r)
    status_dir = repo / ".git" / "genomeos-check"
    assert not (status_dir / f"status-{tree}.json").exists() and not (status_dir / "latest.json").exists()

    def require(*scope: str) -> subprocess.CompletedProcess:
        return subprocess.run(
            [sys.executable, "-m", "genomeos.verdict", "require", "--tree", tree, *scope],
            cwd=repo,
            capture_output=True,
            text=True,
            env=_env(repo),
        )

    # genomeos/verdict.py, with the scope the hooks pass and with none at all
    for scope in (("--scope", "project"), ()):
        got = require(*scope)
        assert got.returncode == verdict.EXIT_REFUSED and "GENOMEOS_VERDICT_OK" not in got.stdout, got
    # and even COPIED to a project verdict's own name it is refused, by its schema before its scope
    shutil.copy(status_dir / f"fast-{tree}.json", status_dir / f"status-{tree}.json")
    for scope in (("--scope", "project"), ()):
        got = require(*scope)
        assert got.returncode == verdict.EXIT_REFUSED and "is not a" in got.stderr, got
    (status_dir / f"status-{tree}.json").unlink()
    # scripts/push_own.sh: with the check skipped it requires a project verdict and refuses this one
    got = subprocess.run(
        ["bash", "scripts/push_own.sh", "--sha", sha],
        cwd=repo,
        capture_output=True,
        text=True,
        env=_env(repo, GENOMEOS_SKIP_CHECK="1"),
    )
    assert got.returncode != 0 and "no green project-wide verdict" in got.stderr, got


def test_the_fast_verdict_never_overwrites_a_project_verdict_for_the_same_tree(h: dict) -> None:
    """The status files are keyed by TREE ALONE, so a fast file under that name would destroy one."""
    read_daily(h)
    repo = h["repo"]
    sha = _commit(repo, {"tests/test_changed.py": PASS}, "a change")
    tree = _git(repo, "rev-parse", f"{sha}^{{tree}}")
    subprocess.run(
        [sys.executable, "-m", "genomeos.verdict", "write", "--exit-code", "0", "--tree-begin", tree]
        + ["--scope", "project", "--started-at", "0"],
        cwd=repo,
        check=True,
        capture_output=True,
        env=_env(repo),
    )
    project = repo / ".git" / "genomeos-check" / f"status-{tree}.json"
    before = project.read_bytes()
    assert push(h, sha).returncode == 0
    assert project.read_bytes() == before


def test_no_reader_of_status_files_names_the_fast_file() -> None:
    """The readers: genomeos/verdict.py, scripts/push_own.sh, scripts/promote_main.sh and check.sh.

    None of them may mention the fast verdict's name, schema or token, and verdict.py's own reader and
    pruner look only at `status-*`; the one shell that reads the fast file, pre-push.sh, reads it only
    through `fast_prepush.py --require-fast`.
    """
    for rel in ("genomeos/verdict.py", "scripts/push_own.sh", "scripts/promote_main.sh", "scripts/check.sh"):
        text = (ROOT / rel).read_text()
        # "fast-forward" is in two of them and is not the file; the file is named fast-<tree>.json
        assert not re.search(r"fast-[{$<*]", text), f"{rel} names the fast verdict file"
        for name in (fp.FAST_SCHEMA, fp.FAST_OK_TOKEN, "fast_prepush"):
            assert name not in text, f"{rel} mentions {name!r}"
    assert verdict.status_path(ROOT, "0" * 40).name.startswith("status-")
    assert fp.fast_status_path(ROOT, "0" * 40).name == f"fast-{'0' * 40}.json"
    assert fp.FAST_SCHEMA != verdict.SCHEMA and fp.FAST_SCOPE not in ("project", "files")
    hook = (ROOT / "scripts" / "pre-push.sh").read_text()
    assert "--require-fast" in hook and "--scope project" in hook


# --------------------------------------------------------------------------- the decision and the record


def test_the_decision_in_each_case(h: dict) -> None:
    repo = h["repo"]
    dev = fp.DEV_REF
    assert fp.gate_decision(repo, dev, env={})["mode"] == "full"  # unread
    read_daily(h)
    assert fp.gate_decision(repo, dev, env={})["mode"] == "fast"
    assert fp.gate_decision(repo, dev, env={"GENOMEOS_FULL_GATE": "0"})["mode"] == "fast"
    assert fp.gate_decision(repo, dev, env={"GENOMEOS_FULL_GATE": "1"})["mode"] == "full"
    assert fp.gate_decision(repo, "refs/heads/promote-1234567", env={})["mode"] == "full"
    later = datetime.now(UTC) + timedelta(seconds=fp.DAILY_READ_MAX_SECONDS + 3600)
    assert fp.gate_decision(repo, dev, env={}, now=later)["mode"] == "full"
    read_daily(h, "failure")
    d = fp.gate_decision(repo, dev, env={"GENOMEOS_FULL_GATE": "1"})
    assert d["mode"] == "full" and len(d["reasons"]) == 2, d  # every reason, not the first
    lg.daily_read_path(repo).write_text("{not json")
    assert fp.gate_decision(repo, dev, env={})["mode"] == "full"


def test_a_green_reading_moves_the_last_green_record_and_a_red_one_does_not(h: dict) -> None:
    repo = h["repo"]
    lg.write(repo, h["green"], "1", "2026-10-01", now=datetime.now(UTC))
    red = _commit(repo, {"tests/test_changed.py": PASS}, "a later sha")
    lg.write_daily(repo, "failure", red, "2", None)
    assert json.loads(lg.last_green_path(repo).read_text())["sha"] == h["green"]
    lg.write_daily(repo, "success", red, "3", None)
    assert json.loads(lg.last_green_path(repo).read_text())["sha"] == red
    assert json.loads(lg.daily_read_path(repo).read_text())["green"] is True


def test_a_reading_of_a_sha_this_checkout_lacks_writes_nothing(h: dict) -> None:
    with pytest.raises(SystemExit):
        lg.write_daily(h["repo"], "success", "f" * 40, None, None)
    assert not lg.daily_read_path(h["repo"]).exists()


def test_the_cache_reading_takes_read_at_from_the_cache_and_reads_its_real_spelling(
    h: dict, tmp_path: Path
) -> None:
    """`ci_status_cache.py` writes sha/id/created, not headSha/databaseId; and a cache is as old as it is."""
    cache = tmp_path / "ci_status.json"
    row = {"status": "completed", "branch": "dev", "sha": h["green"], "id": 7}
    cache.write_text(
        json.dumps(
            {
                "written": "2026-09-30T04:57:59Z",
                "runs": [
                    {**row, "conclusion": "failure", "created": "2026-09-30T04:00:00Z", "id": 8},
                    {**row, "conclusion": "success", "created": "2026-09-29T04:00:00Z"},
                    {**row, "conclusion": "success", "created": "2026-10-01T04:00:00Z", "branch": "main"},
                ],
                "latest_success": {**row, "conclusion": "success", "created": "2026-09-29T04:00:00Z"},
            }
        )
    )
    conclusion, sha, rid, _, read_at = lg.daily_from_cache(cache)
    assert (conclusion, sha, rid) == ("failure", h["green"], "8")
    assert read_at == datetime(2026, 9, 30, 4, 57, 59, tzinfo=UTC)
    assert lg.from_cache(cache)[:2] == (h["green"], "7")
    lg.write_daily(h["repo"], conclusion, sha, rid, None, read_at=read_at)
    assert fp.gate_decision(h["repo"], fp.DEV_REF, env={})["mode"] == "full"


def test_require_fast_refuses_everything_but_a_green_fast_verdict_for_this_tree(h: dict) -> None:
    repo = h["repo"]
    tree = _git(repo, "rev-parse", "HEAD^{tree}")
    assert fp.require_fast(repo, tree) == fp.EXIT_REFUSED  # none
    good = {"schema": fp.FAST_SCHEMA, "scope": "fast", "tree": tree, "tree_begin": tree, "tree_end": tree}
    for bad in (
        {**good, "verdict": "red"},
        {**good, "verdict": "green", "scope": "project"},
        {**good, "verdict": "green", "schema": verdict.SCHEMA},
        {**good, "verdict": "green", "tree_end": "0" * 40},
    ):
        fp.write_fast_status(repo, bad)
        assert fp.require_fast(repo, tree) == fp.EXIT_REFUSED, bad
    fp.write_fast_status(repo, {**good, "verdict": "green"})
    assert fp.require_fast(repo, tree) == 0
