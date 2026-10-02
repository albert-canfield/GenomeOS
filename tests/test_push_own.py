# SPDX-License-Identifier: AGPL-3.0-or-later
"""scripts/push_own.sh: the checked sha is what is pushed, and the ref-lock race is recovered safely.

On 2026-10-02 a lane and the coordinator each lost the same race within an hour. `git push` reads the
remote's refs when it opens the connection, THEN runs the pre-push hook, which here runs a
nine-minute check. A peer pushing inside that window makes the old value git promised stale, and the
remote refuses: "cannot lock ref 'refs/heads/dev': is at X but expected Y". Both diagnosed it by
hand, and the instinct the message invites -- retry with --force -- would discard the peer's work.

So the race is reproduced exactly rather than approximated: the test repository's own pre-push hook
moves the remote while the push is in flight, which is the only way the stale old value arises. Both
branches are planted. A remote tip that is an ANCESTOR of the pushed sha is a fast-forward and is
retried without re-running the check; a tip that is NOT an ancestor is refused.

The refusal is a diagnosis, not a safety property, and the difference was measured rather than
assumed: with the ancestor test stripped out, the diverged case still fails and the peer's commit is
still on the remote, because the retry carries no --force and git refuses a non-fast-forward on its
own. What the refusal buys is an actionable message, one fewer round trip, and the rule written where
someone would otherwise add --force.
"""

import os
import subprocess

import pytest

SCRIPT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts", "push_own.sh")
CANNOT_LOCK = "cannot lock ref"


def git(cwd, *args, check=True):
    return subprocess.run(("git",) + args, cwd=cwd, check=check, capture_output=True, text=True)


def commit(cwd, name):
    (cwd / name).write_text(name)
    git(cwd, "add", name)
    git(cwd, "commit", "-q", "-m", f"add {name}")
    return git(cwd, "rev-parse", "HEAD").stdout.strip()


@pytest.fixture()
def world(tmp_path):
    """A bare remote, our clone `ours`, and a peer clone `peer` that will push mid-flight."""
    bare = tmp_path / "bare.git"
    git(tmp_path, "init", "-q", "--bare", "-b", "dev", str(bare))
    ours = tmp_path / "ours"
    git(tmp_path, "clone", "-q", str(bare), str(ours))
    for d in (ours,):
        git(d, "config", "user.email", "lane@example.invalid")
        git(d, "config", "user.name", "A Lane")
    commit(ours, "base.txt")
    git(ours, "push", "-q", "origin", "HEAD:refs/heads/dev")
    peer = tmp_path / "peer"
    git(tmp_path, "clone", "-q", str(bare), str(peer))
    git(peer, "config", "user.email", "peer@example.invalid")
    git(peer, "config", "user.name", "A Peer")
    git(peer, "remote", "add", "ours", str(ours))
    return bare, ours, peer


def arm_hook(ours, peer, ref_to_push):
    """Our pre-push hook moves the remote while our own push is in flight.

    That is what makes the old value git already promised stale. It honours GENOMEOS_SKIP_CHECK so
    the script's retry does not move the remote a second time -- the real hook does the same.
    """
    hook = ours / ".git" / "hooks" / "pre-push"
    hook.write_text(
        "#!/bin/sh\n"
        '[ "${GENOMEOS_SKIP_CHECK:-0}" = "1" ] && exit 0\n'
        f'git -C "{peer}" push -q origin {ref_to_push}:refs/heads/dev\n'
        "exit 0\n"
    )
    hook.chmod(0o755)


def run_push(ours, sha):
    return subprocess.run(
        ["bash", SCRIPT, "--sha", sha], cwd=ours, capture_output=True, text=True, timeout=120
    )


def test_a_remote_that_moved_to_an_ancestor_is_retried_as_a_fast_forward(world):
    bare, ours, peer = world
    mid = commit(ours, "mid.txt")  # the commit the peer will push from under us
    tip = commit(ours, "tip.txt")  # what we are pushing
    git(peer, "fetch", "-q", "ours", "dev")
    git(peer, "update-ref", "refs/heads/dev", mid)
    arm_hook(ours, peer, "refs/heads/dev")
    r = run_push(ours, tip)
    assert CANNOT_LOCK in r.stdout, (r.stdout, r.stderr)
    assert "an ancestor of" in r.stdout, (r.stdout, r.stderr)
    assert r.returncode == 0, (r.stdout, r.stderr)
    assert git(bare, "rev-parse", "refs/heads/dev").stdout.strip() == tip
    # the peer's commit is still in history: nothing was discarded to get here
    git(bare, "merge-base", "--is-ancestor", mid, tip)


def test_a_remote_that_diverged_is_refused_and_never_forced(world):
    bare, ours, peer = world
    tip = commit(ours, "tip.txt")
    their = commit(peer, "theirs.txt")  # not in our history at all
    arm_hook(ours, peer, their)
    r = run_push(ours, tip)
    assert CANNOT_LOCK in r.stdout, (r.stdout, r.stderr)
    assert "REFUSED to retry" in r.stderr, (r.stdout, r.stderr)
    assert "Never --force here" in r.stderr, (r.stdout, r.stderr)
    assert r.returncode == 1, (r.stdout, r.stderr)
    # the peer's commit is untouched on the remote, which is the whole point
    assert git(bare, "rev-parse", "refs/heads/dev").stdout.strip() == their


def test_the_sha_is_pushed_not_the_branch_name(world):
    """A commit made after the check starts must not ride along on the push."""
    bare, ours, peer = world
    checked = commit(ours, "checked.txt")
    later = commit(ours, "later.txt")
    r = run_push(ours, checked)
    assert r.returncode == 0, (r.stdout, r.stderr)
    assert git(bare, "rev-parse", "refs/heads/dev").stdout.strip() == checked
    assert git(bare, "rev-parse", "refs/heads/dev").stdout.strip() != later
