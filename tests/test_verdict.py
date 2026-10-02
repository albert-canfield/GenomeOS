# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The verdict file: what it records, how it is written, and what reading it refuses.

Each test here stands for one of the five faults that produced genomeos/verdict.py. The two that
matter most are the near-misses at the bottom: a matching green verdict must be ACCEPTED, or the
guard is useless and will be removed by whoever it blocks next; and a stale green verdict from a
DIFFERENT tree must be REFUSED, which is the nine-minute hazard and the reason the file names its
subject at all. The counterfactual at the very bottom strips the tree comparison from a copy and
shows the stale-green case passing, so the guard is credited with something it demonstrably stops.
"""

import json
import os
import subprocess
import time
from pathlib import Path

import pytest

from genomeos import verdict

GIT_ENV = {
    "GIT_AUTHOR_NAME": "t",
    "GIT_AUTHOR_EMAIL": "t@example.invalid",
    "GIT_COMMITTER_NAME": "t",
    "GIT_COMMITTER_EMAIL": "t@example.invalid",
}


def _git(repo: Path, *args: str) -> str:
    env = dict(os.environ)
    env.update(GIT_ENV)
    out = subprocess.run(["git", *args], cwd=str(repo), env=env, check=True, capture_output=True, text=True)
    return out.stdout.strip()


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    """A one-commit git repository, small enough that hashing it is instant."""
    r = tmp_path / "repo"
    r.mkdir()
    _git(r, "init", "-q", "-b", "main")
    (r / "a.py").write_text("x = 1\n")
    (r / ".gitignore").write_text("ignored_dir\n")
    _git(r, "add", "a.py", ".gitignore")
    _git(r, "commit", "-qm", "one")
    return r


@pytest.fixture
def status_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    d = tmp_path / "verdicts"
    d.mkdir()
    monkeypatch.setenv("GENOMEOS_CHECK_STATUS_DIR", str(d))
    return d


# ---------------------------------------------------------------- which tree hash gets recorded


def test_clean_checkout_hashes_to_the_commit_tree(repo: Path) -> None:
    """The claim the consumers depend on: in a clean checkout the working-tree hash IS the commit's.

    scripts/pre-push.sh compares the verdict against `git rev-parse <sha>^{tree}`, so if this ever
    stopped holding the comparison would be a tolerance rather than an identity.
    """
    assert verdict.working_tree_hash(repo) == _git(repo, "rev-parse", "HEAD^{tree}")


def test_an_untracked_file_is_part_of_the_tree_that_was_judged(repo: Path) -> None:
    """Untracked files count, because pytest collects them and ruff reads them.

    This is the fifth fault in miniature: a test file nobody staged was collected and run. The hash
    has to move when it appears, or a verdict could be credited to content that was never on disk.
    """
    before = verdict.working_tree_hash(repo)
    (repo / "test_nobody_staged_me.py").write_text("def test_x():\n    assert True\n")
    assert verdict.working_tree_hash(repo) != before


def test_an_ignored_file_is_not_part_of_the_tree(repo: Path) -> None:
    """Ignored paths do not, or every fetched dataset and cache write would invalidate the verdict."""
    before = verdict.working_tree_hash(repo)
    (repo / "ignored_dir").mkdir()
    (repo / "ignored_dir" / "big.bin").write_bytes(b"\0" * 1024)
    assert verdict.working_tree_hash(repo) == before


def test_hashing_does_not_stage_anything(repo: Path) -> None:
    """The hash is built in a throwaway index: in a shared checkout it must not stage a peer's work."""
    (repo / "peer_half_written.py").write_text("def f(:\n")
    index_before = (repo / ".git" / "index").read_bytes()
    verdict.working_tree_hash(repo)
    assert (repo / ".git" / "index").read_bytes() == index_before
    assert _git(repo, "status", "--porcelain") == "?? peer_half_written.py"


# ---------------------------------------------------------------- the counts, and what green means


@pytest.mark.parametrize(
    ("line", "expected"),
    [
        ("3258 passed, 4 skipped in 412.33s", {"passed": 3258, "failed": 0, "errors": 0, "skipped": 4}),
        ("19 failed, 3239 passed in 400.11s", {"passed": 3239, "failed": 19, "errors": 0, "skipped": 0}),
        ("2 errors in 1.20s", {"passed": 0, "failed": 0, "errors": 2, "skipped": 0}),
        ("1 failed, 1 error in 0.50s", {"passed": 0, "failed": 1, "errors": 1, "skipped": 0}),
    ],
)
def test_counts_are_read_from_the_summary_line(line: str, expected: dict) -> None:
    counts = verdict.parse_pytest_counts(f"some output\n{line}\n")
    assert counts is not None
    for key, value in expected.items():
        assert counts[key] == value


def test_no_summary_line_reads_as_unknown_and_never_as_zero() -> None:
    """A killed run has no summary. Unknown counts must not be reported as nothing having failed."""
    assert verdict.parse_pytest_counts("collecting ...\n") is None
    assert verdict.decide(0, "a" * 40, "a" * 40, None)[0] == "unreadable"


def test_exit_zero_over_failing_counts_is_red() -> None:
    """The backgrounded check that returned 0 while its log said `19 failed`."""
    counts = {"passed": 3239, "failed": 19, "errors": 0, "skipped": 0}
    state, reason = verdict.decide(0, "a" * 40, "a" * 40, counts)
    assert state == "red"
    assert "19 failed" in reason


def test_a_tree_that_moved_mid_run_is_not_a_verdict() -> None:
    """The nine-minute hazard: four failures that were in neither the starting nor the ending tree."""
    counts = {"passed": 10, "failed": 0, "errors": 0, "skipped": 0}
    state, reason = verdict.decide(0, "a" * 40, "b" * 40, counts)
    assert state == "tree_moved"
    assert "changed during the run" in reason


# ---------------------------------------------------------------- the write is all-or-nothing


def test_a_failed_write_leaves_the_previous_verdict_whole(status_dir: Path, monkeypatch) -> None:
    """Either the old file or the complete new one, never half of one, and no stray temp left behind."""
    target = status_dir / "status-x.json"
    verdict._atomic_write_json(target, {"verdict": "green", "mark": "first"})
    first = target.read_text()

    def explode(*args, **kwargs):
        raise OSError("disk full halfway through")

    monkeypatch.setattr(verdict.json, "dump", explode)
    with pytest.raises(OSError):
        verdict._atomic_write_json(target, {"verdict": "red", "mark": "second"})

    assert target.read_text() == first
    assert json.loads(target.read_text())["mark"] == "first"
    assert list(status_dir.glob(".tmp-status-*")) == []


def test_a_reader_never_sees_a_partial_file(status_dir: Path) -> None:
    """Every intermediate state of the target path parses, because only `os.replace` ever touches it."""
    target = status_dir / "status-y.json"
    for mark in ("a", "bb", "ccc" * 5000):
        verdict._atomic_write_json(target, {"verdict": "green", "mark": mark})
        assert json.loads(target.read_text())["mark"] == mark


# ---------------------------------------------------------------- the planted case: $? 0, verdict red


def test_a_background_launch_exits_zero_while_the_verdict_says_red(repo: Path, status_dir: Path) -> None:
    """The planted case, run for real: the launcher's `$?` is 0 and the check it launched went red.

    A stub stands in for scripts/check.sh so the test costs milliseconds, but the shape is the one
    that cost a night: `cmd &` returns the status of the LAUNCH, which succeeded, and the shell that
    reads it is told 0 while the thing it launched failed nineteen tests.
    """
    tree = verdict.working_tree_hash(repo)
    # The log goes OUTSIDE the working tree, as scripts/check.sh keeps its own under .git. A log
    # written into the repo would be untracked content and would move the tree mid-run, which the
    # first draft of this test did and which `tree_moved` duly caught.
    log = status_dir / "stub-pytest.log"
    stub = status_dir / "stub_check.sh"
    stub.write_text(
        "#!/usr/bin/env bash\n"
        f"cd {repo}\n"
        f"printf '19 failed, 3239 passed in 400.11s\\n' > {log}\n"
        "python3 -m genomeos.verdict write --exit-code 19 "
        f"--tree-begin {tree} --scope project --started-at {int(time.time())} "
        f"--pytest-log {log} --failed-leg pytest\n"
        "exit 19\n"
    )
    env = dict(os.environ)
    env["PYTHONPATH"] = str(Path(verdict.__file__).resolve().parents[1])
    launcher = subprocess.run(
        [
            "bash",
            "-c",
            f'bash "{stub}" >/dev/null 2>&1 & launched=$!; echo "launcher_status=$?"; wait $launched',
        ],
        cwd=str(repo),
        env=env,
        capture_output=True,
        text=True,
    )
    assert "launcher_status=0" in launcher.stdout, launcher.stderr

    # ... and the only reading that matters refuses.
    status = verdict.read_status(repo, tree)
    assert status is not None, "the backgrounded run wrote no verdict at all"
    assert status["verdict"] == "red"
    assert status["counts"]["failed"] == 19
    assert verdict.require(repo, tree, scope="project") == verdict.EXIT_REFUSED


def test_no_verdict_at_all_is_refused(repo: Path, status_dir: Path) -> None:
    """A missing verdict is not a pass: a killed run leaves nothing, and nothing must refuse."""
    assert verdict.require(repo, verdict.working_tree_hash(repo), scope="project") == verdict.EXIT_REFUSED


def test_a_file_scoped_verdict_does_not_answer_for_the_project(repo: Path, status_dir: Path) -> None:
    """`check.sh FILE...` lints only those files, so its green is not the green a push needs."""
    tree = verdict.working_tree_hash(repo)
    log = status_dir / "p.log"
    log.write_text("3258 passed, 4 skipped in 412.33s\n")
    verdict.write_status(
        repo,
        exit_code=0,
        tree_begin=tree,
        scope="files",
        scope_files=["a.py"],
        started_at=time.time(),
        pytest_log=log,
    )
    assert verdict.require(repo, tree, scope=None) == 0
    assert verdict.require(repo, tree, scope="project") == verdict.EXIT_REFUSED


# ---------------------------------------------------------------- the two near-misses


def _write_green(repo: Path, tree: str) -> dict:
    log = verdict.status_dir(repo) / "green.log"
    log.write_text("3258 passed, 4 skipped in 412.33s\n")
    return verdict.write_status(
        repo,
        exit_code=0,
        tree_begin=tree,
        scope="project",
        scope_files=[],
        started_at=time.time(),
        pytest_log=log,
    )


def test_a_matching_green_verdict_is_accepted(repo: Path, status_dir: Path, capsys) -> None:
    """The near-miss that matters as much as the refusals: the guard must let the real case through."""
    tree = verdict.working_tree_hash(repo)
    payload = _write_green(repo, tree)
    assert payload["verdict"] == "green"
    assert verdict.require(repo, tree, scope="project") == 0
    printed = capsys.readouterr().out
    assert printed.startswith(f"{verdict.OK_TOKEN} {tree} ")
    assert "passed=3258" in printed


def test_a_stale_green_verdict_from_a_different_tree_is_refused(repo: Path, status_dir: Path) -> None:
    """The nine-minute hazard, as a test rather than a hope.

    The check passed on the tree it started with; nine minutes later a peer's commit moved dev, and
    the sha being pushed carries a different tree. The old verdict is green, complete and honest --
    and it is about something else, so it is not an answer.
    """
    judged = verdict.working_tree_hash(repo)
    _write_green(repo, judged)
    (repo / "a_peer_committed_this.py").write_text("y = 2\n")
    now_being_pushed = verdict.working_tree_hash(repo)
    assert now_being_pushed != judged

    assert verdict.require(repo, judged, scope="project") == 0, "the old tree is still green, as it was"
    assert verdict.require(repo, now_being_pushed, scope="project") == verdict.EXIT_REFUSED


def test_a_verdict_whose_two_readings_disagree_is_refused(repo: Path, status_dir: Path) -> None:
    """A tampered or mid-run-moved file names two trees; neither reading alone may be trusted."""
    tree = verdict.working_tree_hash(repo)
    payload = _write_green(repo, tree)
    payload["tree_end"] = "f" * 40
    payload["verdict"] = "green"
    verdict._atomic_write_json(verdict.status_path(repo, tree), payload)
    assert verdict.require(repo, tree, scope="project") == verdict.EXIT_REFUSED


# ---------------------------------------------------------------- the counterfactual, by removal


def _require_with_the_tree_comparison_stripped(repo: Path) -> int:
    """scripts/pre-push.sh's reading with the tree comparison taken out, and nothing else changed.

    What is left is what the consumers did before this existed: find a verdict, see that it is green,
    and go. It does not ask what the verdict was about, which is the whole of the removed guard.
    """
    try:
        status = json.loads(verdict.latest_path(repo).read_text())
    except (OSError, json.JSONDecodeError):
        return verdict.EXIT_REFUSED
    if status.get("verdict") == "green" and status.get("scope") == "project":
        return 0
    return verdict.EXIT_REFUSED


def test_removing_the_tree_comparison_lets_the_stale_green_through(repo: Path, status_dir: Path) -> None:
    """The guard, credited only with what it is shown to stop.

    Same repository, same verdict file, same stale-green situation as the test above. With the tree
    comparison in place the push is refused; with it stripped out the push is allowed, carrying a
    tree on which nothing was ever run. The harm happens when the guard is removed, so this is a
    protection and not a diagnosis.
    """
    judged = verdict.working_tree_hash(repo)
    _write_green(repo, judged)
    (repo / "a_peer_committed_this.py").write_text("y = 2\n")
    now_being_pushed = verdict.working_tree_hash(repo)

    assert verdict.require(repo, now_being_pushed, scope="project") == verdict.EXIT_REFUSED
    assert _require_with_the_tree_comparison_stripped(repo) == 0


def test_removing_the_comparison_also_admits_any_green_file_lying_about(repo: Path, status_dir: Path) -> None:
    """And the same removal admits any green file in the directory, whatever has happened since."""
    judged = verdict.working_tree_hash(repo)
    _write_green(repo, judged)
    for _ in range(3):
        (repo / f"more_{time.time_ns()}.py").write_text("z = 3\n")
    assert _require_with_the_tree_comparison_stripped(repo) == 0
    assert verdict.require(repo, verdict.working_tree_hash(repo), scope="project") == verdict.EXIT_REFUSED
