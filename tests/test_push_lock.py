# SPDX-License-Identifier: AGPL-3.0-or-later
"""The push lock's takeover branch, which had never been shown to fire.

One push's check runs at a time (scripts/pre-push.sh), and a push whose process dies holding
the lock must not stop every later push for an hour. The branch that takes over a dead
holder was written on 2026-10-02 and nothing exercised it. On the same day a stale log line
naming a holder that had since exited was read as that branch failing, and the proposed
remedy was to delete the lock by hand -- which would have taken a live peer's lock. So the
branch is pinned both ways: a dead holder is taken over, and a live one is never touched.

The lock lives at $TMPDIR/genomeos-push.lock, so every test here sets TMPDIR to its own
directory. Without that these tests would contend for the real lock with whatever push is
running in this checkout.
"""

import os
import subprocess
import time

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

TAKEN = "which is gone; taking it"
# a ref that is neither dev nor main: the hook acquires the lock, finds nothing to check and
# exits, so these tests cost a lock acquisition and not a nine-minute suite
OTHER_REF = "refs/heads/not-a-push-target"
SHA = "0" * 39 + "1"
STDIN = f"{OTHER_REF} {SHA} {OTHER_REF} {SHA}\n"


def run_hook(tmpdir, timeout):
    return subprocess.run(
        ["bash", "scripts/pre-push.sh", "origin", "https://example.invalid/x.git"],
        cwd=REPO_ROOT,
        input=STDIN,
        capture_output=True,
        text=True,
        timeout=timeout,
        env={"PATH": "/usr/bin:/bin:/usr/sbin:/sbin", "TMPDIR": str(tmpdir), "HOME": str(tmpdir)},
    )


@pytest.fixture()
def lock_dir(tmp_path):
    d = tmp_path / "genomeos-push.lock"
    d.mkdir()
    return d


def a_dead_pid():
    p = subprocess.Popen(["true"])
    p.wait()
    # the pid of a reaped child is not running; a high invented number might be in use
    return p.pid


def test_a_lock_whose_holder_is_dead_is_taken_over(tmp_path, lock_dir):
    (lock_dir / "pid").write_text(f"{a_dead_pid()}\n")
    r = run_hook(tmp_path, timeout=30)
    assert TAKEN in r.stdout, (r.stdout, r.stderr)
    assert r.returncode == 0, (r.stdout, r.stderr)
    assert not lock_dir.exists() or (lock_dir / "pid").read_text().strip() != "", r.stdout


def test_a_lock_whose_holder_is_alive_is_never_taken(tmp_path, lock_dir):
    holder = subprocess.Popen(["sleep", "30"])
    try:
        (lock_dir / "pid").write_text(f"{holder.pid}\n")
        started = time.monotonic()
        with pytest.raises(subprocess.TimeoutExpired) as caught:
            run_hook(tmp_path, timeout=8)
        waited = time.monotonic() - started
        out = (
            (caught.value.stdout or b"").decode()
            if isinstance(caught.value.stdout, bytes)
            else (caught.value.stdout or "")
        )
        assert TAKEN not in out, out
        assert waited >= 8, waited
        assert (lock_dir / "pid").read_text().strip() == str(holder.pid)
    finally:
        holder.kill()
        holder.wait()
