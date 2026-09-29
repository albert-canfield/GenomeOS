# SPDX-License-Identifier: AGPL-3.0-or-later
"""The promotion gate, against a throwaway remote and a stand-in for `gh`.

Every run here is a dry run or a refusal: nothing in this file lets the script push, and the
remote it could reach is a bare repository under the test's temporary directory. What is checked
is the decision: which shas the gate lets through and which it refuses, read from the check runs
a stand-in `gh` reports.

2026-09-29, beside the above: the `--prepare` tests do let the script push, one branch named
`promote-<first 7 of the sha>`, to that bare repository and nowhere else. These tests check the
script's decisions and git operations against a local remote. GitHub's enforcement of main's branch
rule is not simulated, and nothing here shows what GitHub would accept: the stand-in `gh` reports
the check runs and events a test gives it, and decides nothing.
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
            "check_suite": {"id": c.get("suite")},
        }
        for c in state["checks"].get(sha, [])
    ]
    print(json.dumps({"total_count": len(runs), "check_runs": runs}))
elif args[:1] == ["api"] and "/actions/workflows/ci.yml/runs" in args[1]:
    print(json.dumps({"workflow_runs": [{"head_sha": s} for s in state["runs"]]}))
elif args[:1] == ["api"] and "/actions/runs?head_sha=" in args[1]:
    sha = args[1].split("head_sha=")[1].split("&")[0]
    runs = [
        {"head_sha": sha, "check_suite_id": c["suite"], "event": c["event"]}
        for c in state["checks"].get(sha, [])
        if c.get("event")
    ]
    print(json.dumps({"total_count": len(runs), "workflow_runs": runs}))
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

    def check_event(self, sha: str, conclusion: str, when: str, event: str | None) -> None:
        """A `test` check run whose workflow run was triggered by `event`; None: no workflow run found."""
        self.check(sha, conclusion, when)
        self.state["checks"][sha][-1]["event"] = event

    def remote_main(self) -> str:
        return _git(self.bare, "rev-parse", "refs/heads/main")

    def remote_refs(self) -> dict[str, str]:
        out = _git(self.bare, "for-each-ref", "--format=%(refname) %(objectname)")
        return dict(line.split() for line in out.splitlines())

    def run(self, *args: str) -> subprocess.CompletedProcess:
        # the only remote this script can reach in a test is the bare repository made above
        assert Path(_git(self.work, "remote", "get-url", "origin")) == self.bare
        state = self.tmp / "gh_state.json"
        # 2026-09-29: each check run gets its own check suite, and each suite a workflow run whose
        # event is `schedule` unless the test said otherwise (check_event)
        entries = [e for es in self.state["checks"].values() for e in es]
        for entry in entries:
            if "suite" not in entry:
                entry["suite"] = 1 + max((e.get("suite", 0) for e in entries), default=0)
            entry.setdefault("event", "schedule")
        state.write_text(json.dumps(self.state))
        env = {**os.environ, "PATH": self.path, "FAKE_GH_STATE": str(state)}
        env.pop("PROMOTE_REPO", None)
        # 2026-09-29: --prepare pushes to the bare repository; no hook of this machine runs on it
        env.update(GIT_CONFIG_COUNT="1", GIT_CONFIG_KEY_0="core.hooksPath", GIT_CONFIG_VALUE_0="/dev/null")
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


# 2026-09-29: the event of the run the gate relies on, --prepare, --verify, and the retired scripts.
# Every sentence below is about the script's own decisions against a local remote; none of them is
# about what GitHub's branch rule would do.

PROMOTE_SH = SCRIPT.parent / "promote.sh"


def _gh_calls(repo) -> str:
    log = repo.tmp / "gh_state.json.log"
    return log.read_text() if log.exists() else ""


@pytest.mark.parametrize("event", ["workflow_dispatch", "schedule"])
def test_a_green_manual_or_scheduled_run_is_a_pre_check_only(repo, event):
    repo.check_event(repo.c, "success", "2026-09-28T10:00:00Z", event)
    before = repo.remote_refs()
    r = repo.run(repo.c)
    assert r.returncode == 0, r.stderr
    assert f"triggered by {event}" in r.stdout
    assert "CI pre-check passed" in r.stdout
    assert "a pre-check only; it does not make the sha eligible" in r.stdout
    assert "eligible for protected promotion:" not in r.stdout
    assert "is an inference, not established" in r.stdout
    assert f"actions/runs?head_sha={repo.c}" in _gh_calls(repo)
    assert "that push is no longer made" in r.stdout
    assert repo.remote_refs() == before, "a dry run moves nothing"


@pytest.mark.parametrize("event", ["pull_request", "push"])
def test_a_green_pull_request_or_push_run_is_eligible(repo, event):
    repo.check_event(repo.c, "success", "2026-09-28T10:00:00Z", event)
    r = repo.run(repo.c)
    assert r.returncode == 0, r.stderr
    assert f"eligible for protected promotion: the `test` run relied on was triggered by {event}" in r.stdout
    assert "pre-check only" not in r.stdout


def test_the_event_is_that_of_the_run_the_gate_relies_on(repo):
    """An earlier green pull_request run does not lend its event to a later green manual run."""
    repo.check_event(repo.c, "success", "2026-09-27T10:00:00Z", "pull_request")
    repo.check_event(repo.c, "success", "2026-09-28T10:00:00Z", "workflow_dispatch")
    r = repo.run(repo.c)
    assert r.returncode == 0, r.stderr
    assert "triggered by workflow_dispatch" in r.stdout and "pre-check only" in r.stdout


def test_a_run_whose_workflow_run_is_not_found_is_a_pre_check_only(repo):
    repo.check_event(repo.c, "success", "2026-09-28T10:00:00Z", None)
    r = repo.run(repo.c)
    assert r.returncode == 0, r.stderr
    assert "triggered by unknown" in r.stdout and "pre-check only" in r.stdout


def test_push_is_refused_before_anything_is_read_or_moved(repo):
    repo.check_event(repo.c, "success", "2026-09-28T10:00:00Z", "pull_request")
    before = repo.remote_refs()
    for args in (["--push", repo.c], [repo.c, "--push"], ["--dry-run", "--push", repo.c]):
        r = repo.run(*args)
        assert r.returncode == 2
        assert "--push is retired" in r.stderr and "--prepare SHA" in r.stderr
    assert _gh_calls(repo) == "", "nothing was asked of GitHub"
    assert repo.remote_refs() == before


def test_prepare_pushes_only_the_frozen_branch(repo):
    repo.check_event(repo.c, "success", "2026-09-28T10:00:00Z", "workflow_dispatch")
    before = repo.remote_refs()
    r = repo.run("--prepare", repo.c)
    assert r.returncode == 0, r.stderr
    branch = f"refs/heads/promote-{repo.c[:7]}"
    assert repo.remote_refs() == {**before, branch: repo.c}, "one new branch, at exactly that sha"
    assert repo.remote_main() == repo.a
    assert f"https://github.com/owner/repo/compare/main...promote-{repo.c[:7]}?expand=1" in r.stdout
    assert f"scripts/promote_main.sh --verify {repo.c}" in r.stdout
    assert "pre-check only" in r.stdout
    reads = ("repo view", "api repos/owner/repo/commits/", "api repos/owner/repo/actions/runs?")
    assert all(line.startswith(reads) for line in _gh_calls(repo).splitlines()), "only reads of GitHub"
    # again: the branch is already there at that sha, so nothing more is pushed
    r = repo.run("--prepare", repo.c)
    assert r.returncode == 0, r.stderr
    assert "nothing pushed" in r.stdout
    assert repo.remote_refs() == {**before, branch: repo.c}


def test_prepare_refuses_a_branch_that_exists_at_another_sha(repo):
    repo.check_event(repo.c, "success", "2026-09-28T10:00:00Z", "workflow_dispatch")
    branch = f"refs/heads/promote-{repo.c[:7]}"
    _git(repo.work, "push", "-q", "origin", f"{repo.b}:{branch}")
    before = repo.remote_refs()
    r = repo.run("--prepare", repo.c)
    assert r.returncode == 1
    assert f"already has promote-{repo.c[:7]} at {repo.b[:7]}" in r.stderr
    assert repo.remote_refs() == before, "neither the branch nor main moved"


def test_prepare_makes_the_same_checks_first(repo):
    repo.check_event(repo.c, "failure", "2026-09-28T10:00:00Z", "workflow_dispatch")
    repo.check_event(repo.s, "success", "2026-09-28T10:00:00Z", "workflow_dispatch")
    before = repo.remote_refs()
    r = repo.run("--prepare", repo.c)
    assert r.returncode == 1 and "concluded failure" in r.stderr
    r = repo.run("--prepare", repo.s)
    assert r.returncode == 1 and "is not on origin/dev" in r.stderr
    assert repo.remote_refs() == before


def _merge_on_main(repo, *parents: str) -> str:
    """Stand in for the owner's merge: a commit with the chosen sha's tree, pushed to main."""
    args = ["commit-tree", f"{repo.c}^{{tree}}", "-m", "merge"]
    for parent in parents:
        args += ["-p", parent]
    merge = _git(repo.work, *args)
    _git(repo.work, "push", "-q", "origin", f"{merge}:refs/heads/main")
    return merge


def test_verify_accepts_main_with_the_sha_s_tree_and_the_sha_as_ancestor(repo):
    _merge_on_main(repo, repo.a, repo.c)
    r = repo.run("--verify", repo.c)
    assert r.returncode == 0, r.stderr
    assert "verified" in r.stdout
    assert _gh_calls(repo) == "", "--verify reads git only"


def test_verify_refuses_main_before_the_merge(repo):
    r = repo.run("--verify", repo.c)
    assert r.returncode == 1 and "does not have" in r.stderr and "file tree" in r.stderr


def test_verify_refuses_the_same_tree_without_the_sha_as_ancestor(repo):
    """A squash merge carries the files but not the sha."""
    _merge_on_main(repo, repo.a)
    r = repo.run("--verify", repo.c)
    assert r.returncode == 1 and "is not an ancestor" in r.stderr


def test_verify_refuses_the_sha_as_ancestor_with_another_tree(repo):
    _git(repo.work, "push", "-q", "origin", f"{repo.c}:refs/heads/main")
    r = repo.run("--verify", repo.b)
    assert r.returncode == 1 and "does not have" in r.stderr


@pytest.mark.parametrize(
    "args",
    [
        ["--prepare", "--latest-green"],
        ["--verify", "--latest-green"],
        ["--prepare", "--verify", "HEAD"],
        ["--dry-run", "--prepare", "HEAD"],
        ["--prepare"],
        ["--verify"],
    ],
)
def test_new_usage_errors_exit_2_before_any_check(repo, args):
    before = repo.remote_refs()
    r = repo.run(*args)
    assert r.returncode == 2
    assert _gh_calls(repo) == "", "nothing was asked of GitHub"
    assert repo.remote_refs() == before


def test_help_prints_the_new_usage_only(repo):
    r = repo.run("--help")
    assert r.returncode == 0
    assert "--prepare" in r.stderr and "--verify" in r.stderr and "--push" not in r.stderr


def test_promote_sh_is_retired_and_runs_nothing(tmp_path):
    bindir = tmp_path / "bin"
    bindir.mkdir()
    log = tmp_path / "calls.log"
    for tool in ("git", "gh", "uv"):
        stub = bindir / tool
        stub.write_text(f'#!/bin/sh\necho "{tool} $*" >> "{log}"\n')
        stub.chmod(0o755)
    env = {**os.environ, "PATH": f"{bindir}{os.pathsep}{os.environ['PATH']}"}
    r = subprocess.run(["bash", str(PROMOTE_SH)], cwd=tmp_path, env=env, capture_output=True, text=True)
    assert r.returncode != 0
    assert "retired" in r.stderr and "scripts/promote_main.sh" in r.stderr
    assert not log.exists(), "no git, gh or uv call"
