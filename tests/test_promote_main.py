# SPDX-License-Identifier: AGPL-3.0-or-later
"""The promotion gate, against a throwaway remote and a stand-in for `gh`.

Every run here is a dry run or a refusal: nothing in this file lets the script push, and the
remote it could reach is a bare repository under the test's temporary directory. What is checked
is the decision: which shas the gate lets through and which it refuses, read from the check runs
a stand-in `gh` reports.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "promote_main.sh"

pytestmark = pytest.mark.skipif(
    not (shutil.which("bash") and shutil.which("git") and shutil.which("python3")),
    reason="needs bash, git and python3",
)

FAKE_GH = """\
import json, os, sys
state = json.load(open(os.environ["FAKE_GH_STATE"]))
args = sys.argv[1:]
with open(os.environ["FAKE_GH_STATE"] + ".log", "a") as log:
    log.write(" ".join(args) + "\\n")
if args[:2] == ["repo", "view"]:
    print("owner/repo")
elif args[:1] == ["api"] and "/check-runs" in args[1]:
    sha = args[1].split("/commits/")[1].split("/")[0]
    runs = [
        {
            "status": "completed",
            "conclusion": c["conclusion"],
            "completed_at": c["completed_at"],
            "html_url": "https://example.invalid/" + c["conclusion"],
            "app": {"slug": c.get("app", "github-actions")},
        }
        for c in state["checks"].get(sha, [])
    ]
    print(json.dumps({"total_count": len(runs), "check_runs": runs}))
elif args[:1] == ["api"] and "/actions/workflows/ci.yml/runs" in args[1]:
    print(json.dumps({"workflow_runs": [{"head_sha": s} for s in state["runs"]]}))
else:
    sys.exit("fake gh: unexpected call " + " ".join(args))
"""


def _git(cwd: Path, *args: str) -> str:
    env = {
        **os.environ,
        "GIT_AUTHOR_NAME": "t",
        "GIT_AUTHOR_EMAIL": "t@example.invalid",
        "GIT_COMMITTER_NAME": "t",
        "GIT_COMMITTER_EMAIL": "t@example.invalid",
    }
    out = subprocess.run(
        ["git", "-c", "core.hooksPath=/dev/null", *args],
        cwd=cwd,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    return out.stdout.strip()


class Repo:
    """main at A; dev at A-B-C; a commit S off A that is on no remote branch."""

    def __init__(self, tmp: Path):
        self.tmp = tmp
        self.bare = tmp / "remote.git"
        self.work = tmp / "work"
        _git(tmp, "init", "-q", "--bare", str(self.bare))
        self.work.mkdir()
        _git(self.work, "init", "-q", "-b", "dev")
        _git(self.work, "remote", "add", "origin", str(self.bare))
        self.a = self._commit("a")
        _git(self.work, "push", "-q", "origin", "dev:dev", "dev:main")
        self.b = self._commit("b")
        self.c = self._commit("c")
        _git(self.work, "push", "-q", "origin", "dev:dev")
        _git(self.work, "checkout", "-q", "-b", "side", self.a)
        self.s = self._commit("s")
        _git(self.work, "checkout", "-q", "dev")
        self.state = {"checks": {}, "runs": []}
        bindir = tmp / "bin"
        bindir.mkdir()
        (tmp / "fake_gh.py").write_text(FAKE_GH)
        gh = bindir / "gh"
        # a shell wrapper, because a shebang cannot carry a path with a space in it
        gh.write_text(f'#!/bin/sh\nexec "{sys.executable}" "{tmp / "fake_gh.py"}" "$@"\n')
        gh.chmod(0o755)
        self.path = f"{bindir}{os.pathsep}{os.environ['PATH']}"

    def _commit(self, name: str) -> str:
        (self.work / name).write_text(name)
        _git(self.work, "add", name)
        _git(self.work, "commit", "-q", "-m", name)
        return _git(self.work, "rev-parse", "HEAD")

    def check(self, sha: str, conclusion: str, when: str, app: str = "github-actions") -> None:
        entry = {"conclusion": conclusion, "completed_at": when, "app": app}
        self.state["checks"].setdefault(sha, []).append(entry)

    def remote_main(self) -> str:
        return _git(self.bare, "rev-parse", "refs/heads/main")

    def run(self, *args: str) -> subprocess.CompletedProcess:
        # the only remote this script can reach in a test is the bare repository made above
        assert Path(_git(self.work, "remote", "get-url", "origin")) == self.bare
        state = self.tmp / "gh_state.json"
        state.write_text(json.dumps(self.state))
        env = {**os.environ, "PATH": self.path, "FAKE_GH_STATE": str(state)}
        env.pop("PROMOTE_REPO", None)
        return subprocess.run(
            ["bash", str(SCRIPT), *args], cwd=self.work, env=env, capture_output=True, text=True
        )


@pytest.fixture
def repo(tmp_path: Path) -> Repo:
    return Repo(tmp_path)


def test_a_green_sha_is_let_through_as_a_dry_run_by_default(repo):
    repo.check(repo.c, "success", "2026-09-28T10:00:00Z")
    for args in ([repo.c], ["--dry-run", repo.c]):
        r = repo.run(*args)
        assert r.returncode == 0, r.stderr
        assert f"would run: git push origin {repo.c}:refs/heads/main" in r.stdout
        assert "2 commits" in r.stdout
    assert repo.remote_main() == repo.a, "a dry run moves nothing"


def test_the_latest_run_decides_and_a_later_red_one_refuses(repo):
    repo.check(repo.c, "success", "2026-09-27T10:00:00Z")
    repo.check(repo.c, "failure", "2026-09-28T10:00:00Z")
    r = repo.run(repo.c)
    assert r.returncode == 1
    assert "concluded failure" in r.stderr and "2 on record" in r.stderr
    assert repo.remote_main() == repo.a


def test_a_red_run_followed_by_a_green_rerun_is_let_through(repo):
    repo.check(repo.c, "failure", "2026-09-27T10:00:00Z")
    repo.check(repo.c, "success", "2026-09-28T10:00:00Z")
    assert repo.run(repo.c).returncode == 0


@pytest.mark.parametrize("conclusion", ["cancelled", "skipped", "neutral", "timed_out"])
def test_only_success_counts(repo, conclusion):
    """Main's branch rule accepts skipped and neutral; the gate does not."""
    repo.check(repo.c, conclusion, "2026-09-28T10:00:00Z")
    r = repo.run(repo.c)
    assert r.returncode == 1 and f"concluded {conclusion}" in r.stderr


def test_a_sha_nothing_has_tested_is_refused(repo):
    r = repo.run(repo.b)
    assert r.returncode == 1 and "has no completed `test` run" in r.stderr


def test_a_green_check_named_test_from_another_app_does_not_count(repo):
    repo.check(repo.c, "success", "2026-09-28T10:00:00Z", app="someone-elses-ci")
    r = repo.run(repo.c)
    assert r.returncode == 1 and "has no completed `test` run" in r.stderr


def test_a_green_sha_that_is_not_on_dev_is_refused(repo):
    repo.check(repo.s, "success", "2026-09-28T10:00:00Z")
    r = repo.run(repo.s)
    assert r.returncode == 1 and "is not on origin/dev" in r.stderr


def test_a_move_that_is_not_a_fast_forward_is_refused(repo):
    _git(repo.work, "push", "-q", "origin", f"{repo.c}:refs/heads/main")
    repo.check(repo.b, "success", "2026-09-28T10:00:00Z")
    r = repo.run(repo.b)
    assert r.returncode == 1 and "not a fast-forward" in r.stderr
    assert repo.remote_main() == repo.c


def test_main_already_there_is_nothing_to_do(repo):
    r = repo.run(repo.a)
    assert r.returncode == 0 and "already at" in r.stdout


def test_latest_green_takes_the_newest_green_sha_between_main_and_dev(repo):
    # newest first, as GitHub lists runs: one sha this clone lacks, one off dev, then B
    repo.state["runs"] = ["f" * 40, repo.s, repo.b, repo.a]
    repo.check(repo.b, "success", "2026-09-28T10:00:00Z")
    r = repo.run("--latest-green")
    assert r.returncode == 0, r.stderr
    assert f"would run: git push origin {repo.b}:refs/heads/main" in r.stdout
    assert repo.remote_main() == repo.a


def test_latest_green_refuses_when_no_green_sha_is_between(repo):
    repo.state["runs"] = [repo.s]
    r = repo.run("--latest-green")
    assert r.returncode == 1 and "no sha between origin/main and origin/dev" in r.stderr


@pytest.mark.parametrize("args", [[], ["--latest-green", "HEAD"], ["--bogus", "HEAD"], ["HEAD", "HEAD~1"]])
def test_usage_errors_exit_2_before_any_check(repo, args):
    r = repo.run(*args)
    assert r.returncode == 2
    assert not (repo.tmp / "gh_state.json.log").exists(), "nothing was asked of GitHub"
