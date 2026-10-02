# SPDX-License-Identifier: AGPL-3.0-or-later
"""One full `scripts/check.sh` at a time on the machine, with the takeover branch pinned both ways.

WHAT THIS PLANTS. On 2026-10-02 the acceptance rule in force was "a status-file verdict of the
committed tree", and the only way to earn one was the whole ~4,500-test suite. Lanes ran it three at
once. The machine paged -- swapouts rising 90,000 to 114,000 per sample, 14 GB of swap on disk, free
disk under the 10 GB floor -- and the capacity gate then refused the heavy jobs real work needed, so
a hygiene rule finished by holding a release. `scripts/suite_lock.sh` makes those runs queue.

WHY IT IS THE PUSH LOCK'S ALGORITHM. `scripts/pre-push.sh` already serialises one push's check at a
time and `tests/test_push_lock.py` pins its two cases. A second lock in the same project behaving
differently would be worse than no second lock, so the branches here are the same three, tested the
same way:

  - a DEAD holder is taken over, because a crashed run must not block every later run for the cap;
  - a LIVE holder is NEVER taken -- on the same day a stale log line naming a holder that had since
    exited was read as the takeover branch failing, and the proposed remedy, deleting the lock by
    hand, would have taken a live peer's lock;
  - a waiter RE-PRINTS THE HOLDER WHEN IT CHANGES, so somebody watching learns who they are waiting
    for instead of reading one line and then silence for half an hour.

HOW IT IS RUN. Through the REAL `scripts/check.sh`, which is what has to take the lock, copied with
the real `scripts/suite_lock.sh` into a throwaway git repository with `uv` and `shellcheck` stubbed
on PATH, after the pattern `tests/test_check_sh_arguments.py` set: the shape is the real one and the
legs cost milliseconds instead of twelve minutes. Only `genomeos.verdict` is real.

The lock lives at `$TMPDIR/genomeos-suite.lock`, so every test here sets TMPDIR to its own directory
and clears `GENOMEOS_SUITE_LOCK_HELD`. Without the first these tests would contend with whatever
check is running in this checkout; without the second, a run of this file from INSIDE a project-wide
`check.sh` would inherit that run's lock marker, take the re-entrant path and test nothing.
"""

from __future__ import annotations

import json
import os
import subprocess
import time
from pathlib import Path

import pytest

TESTS_DIR = Path(__file__).resolve().parent
REPO_ROOT = TESTS_DIR.parent
CHECK_SH = REPO_ROOT / "scripts" / "check.sh"
SUITE_LOCK_SH = REPO_ROOT / "scripts" / "suite_lock.sh"

#: The takeover line, quoted from scripts/suite_lock.sh and shared word for word with
#: scripts/pre-push.sh, which is the point: one lock behaviour, one sentence for it.
TAKEN = "which is gone; taking it"

_UV_STUB = """#!/usr/bin/env bash
# Stands in for `uv` so the legs cost nothing; this file is about the lock, not about the legs.
set -u
shift || true
while [ "$#" -gt 0 ]; do
  case "$1" in --quiet | --no-project | --with | --with=*) [ "$1" = "--with" ] && shift; shift ;;
  --) shift; break ;;
  *) break ;;
  esac
done
case "${1:-}" in
  pytest) printf '%s\\n' "3 passed in 0.10s"; exit 0 ;;
  *) exit 0 ;;
esac
"""

_SHELLCHECK_STUB = "#!/usr/bin/env bash\nexit 0\n"

_WAITING_FOR = "waiting for the suite running in process"


@pytest.fixture()
def harness(tmp_path: Path) -> dict:
    """A throwaway repo with the real check.sh and suite_lock.sh, its own TMPDIR and status dir."""
    repo = tmp_path / "repo"
    (repo / "scripts").mkdir(parents=True)
    for src in (CHECK_SH, SUITE_LOCK_SH):
        dst = repo / "scripts" / src.name
        dst.write_bytes(src.read_bytes())
        dst.chmod(0o755)
    (repo / "keep.py").write_text("x = 1\n")
    subprocess.run(["git", "init", "-q", "-b", "main", "."], cwd=str(repo), capture_output=True)

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for name, body in (("uv", _UV_STUB), ("shellcheck", _SHELLCHECK_STUB)):
        (bin_dir / name).write_text(body)
        (bin_dir / name).chmod(0o755)

    # TMPDIR of its own: the lock path is derived from it, so this is what keeps the test off the
    # real lock. The status dir likewise, so a verdict written here never lands in .git.
    lock_home = tmp_path / "tmp"
    lock_home.mkdir()
    status_dir = tmp_path / "status"
    env = dict(os.environ)
    env.pop("GENOMEOS_SUITE_LOCK_HELD", None)
    env.update(
        PATH=f"{bin_dir}:{env['PATH']}",
        PYTHONPATH=str(REPO_ROOT),
        GENOMEOS_CHECK_STATUS_DIR=str(status_dir),
        TMPDIR=str(lock_home),
    )
    return {
        "repo": repo,
        "env": env,
        "status_dir": status_dir,
        "lock": lock_home / "genomeos-suite.lock",
    }


def _spawn(harness: dict) -> subprocess.Popen:
    return subprocess.Popen(
        ["bash", "scripts/check.sh", "keep.py"],
        cwd=str(harness["repo"]),
        env=harness["env"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


def _held_lock(harness: dict, pid: int) -> Path:
    lock = harness["lock"]
    lock.mkdir()
    (lock / "pid").write_text(f"{pid}\n")
    return lock


def _a_dead_pid() -> int:
    p = subprocess.Popen(["true"])
    p.wait()
    # the pid of a reaped child is not running; a high invented number might be in use
    return p.pid


def _verdict(harness: dict) -> dict:
    latest = harness["status_dir"] / "latest.json"
    assert latest.is_file(), f"no verdict file was written; {sorted(harness['status_dir'].glob('*'))}"
    return json.loads(latest.read_text())


def _still_waiting_after(proc: subprocess.Popen, seconds: float) -> str:
    """Require the run to be STILL WAITING after `seconds`, then kill it and return its output.

    `Popen.communicate` raises TimeoutExpired with no output attached, so the output is collected
    from a second `communicate` after the kill. Reading it from the exception would read an empty
    string and every assertion below would pass on nothing.
    """
    try:
        proc.communicate(timeout=seconds)
        raise AssertionError("the waiter exited; it was supposed to still be waiting")
    except subprocess.TimeoutExpired:
        pass
    finally:
        proc.kill()
    out, err = proc.communicate()
    return f"{out}{err}"


def test_a_suite_lock_whose_holder_is_dead_is_taken_over(harness: dict) -> None:
    """A run that crashed holding the lock must not stop every later check for the cap."""
    lock = _held_lock(harness, _a_dead_pid())
    proc = _spawn(harness)
    out, err = proc.communicate(timeout=90)
    assert TAKEN in out, (out, err)
    assert proc.returncode == 0, (out, err)
    # The verdict is read from the artefact, never from the exit code.
    verdict = _verdict(harness)
    assert verdict["exit_code"] == 0, verdict
    # The run WAS serialised, so the note must not claim otherwise: the "not serialised" note exists
    # only for a lock that could not be had, and a green that quietly skipped its own serialisation
    # is the class of green this project keeps writing up.
    assert "not serialised" not in (verdict.get("note") or ""), verdict
    # And it gave the lock back, so the next run in the queue does not wait out the stale-lock timer.
    assert not lock.exists(), sorted(p.name for p in lock.parent.iterdir())


def test_a_suite_lock_whose_holder_is_alive_is_never_taken(harness: dict) -> None:
    """The dangerous direction: a lock removed because it looked stale is a live peer's lock."""
    holder = subprocess.Popen(["sleep", "30"])
    try:
        lock = _held_lock(harness, holder.pid)
        started = time.monotonic()
        out = _still_waiting_after(_spawn(harness), 8)
        assert time.monotonic() - started >= 8
        assert TAKEN not in out, out
        assert f"{_WAITING_FOR} {holder.pid}" in out, out
        assert (lock / "pid").read_text().strip() == str(holder.pid)
        # It never reached a leg, so it must not have written a verdict for this tree either.
        assert not (harness["status_dir"] / "latest.json").exists()
    finally:
        holder.kill()
        holder.wait()


def test_a_nested_check_does_not_wait_for_the_lock_its_own_parent_holds(harness: dict) -> None:
    """The re-entrant path, which is a deadlock if it is wrong rather than a missing optimisation.

    `tests/test_check_sh_arguments.py` runs the real `scripts/check.sh` from inside the suite, so a
    project-wide run reaches the lock from a process whose ancestor already holds it. Waiting there
    would stop the whole suite until the cap. The marker must name a LIVE process, and the holder
    here is a DIFFERENT live process, so this cannot be satisfied by the dead-holder branch.
    """
    outer = subprocess.Popen(["sleep", "40"])
    other = subprocess.Popen(["sleep", "40"])
    try:
        lock = _held_lock(harness, other.pid)
        harness["env"]["GENOMEOS_SUITE_LOCK_HELD"] = str(outer.pid)
        proc = _spawn(harness)
        out, err = proc.communicate(timeout=90)
        assert proc.returncode == 0, (out, err)
        assert _WAITING_FOR not in out, out
        assert TAKEN not in out, out
        assert f"process {outer.pid} holds it and started this run" in out, out
        # It held no lock of its own, so it must not have removed the one that was there.
        assert (lock / "pid").read_text().strip() == str(other.pid)
        assert _verdict(harness)["exit_code"] == 0
    finally:
        for p in (outer, other):
            p.kill()
            p.wait()


def test_a_stale_marker_naming_a_dead_process_does_not_skip_the_lock(harness: dict) -> None:
    """The re-entrant path is not a bypass: a leftover variable cannot stand in for the lock."""
    harness["env"]["GENOMEOS_SUITE_LOCK_HELD"] = str(_a_dead_pid())
    holder = subprocess.Popen(["sleep", "30"])
    try:
        lock = _held_lock(harness, holder.pid)
        out = _still_waiting_after(_spawn(harness), 8)
        assert f"{_WAITING_FOR} {holder.pid}" in out, out
        assert TAKEN not in out, out
        assert (lock / "pid").read_text().strip() == str(holder.pid)
    finally:
        holder.kill()
        holder.wait()


def test_a_waiter_reprints_the_holder_when_it_changes(harness: dict) -> None:
    """Silence is what gets a check backgrounded, and a backgrounded check is where exit codes lied.

    Both holders are alive throughout, so nothing here can be satisfied by the takeover branch.
    """
    first = subprocess.Popen(["sleep", "40"])
    second = subprocess.Popen(["sleep", "40"])
    try:
        lock = _held_lock(harness, first.pid)
        proc = _spawn(harness)
        # Long enough for the first announcement, well inside the five-second poll.
        time.sleep(2)
        (lock / "pid").write_text(f"{second.pid}\n")
        out = _still_waiting_after(proc, 12)
        assert f"{_WAITING_FOR} {first.pid}" in out, out
        assert f"{_WAITING_FOR} {second.pid}" in out, out
        assert TAKEN not in out, out
        assert (lock / "pid").read_text().strip() == str(second.pid)
    finally:
        for p in (first, second):
            p.kill()
            p.wait()
