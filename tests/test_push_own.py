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
import time

import pytest

from genomeos import verdict

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO_ROOT, "scripts", "push_own.sh")
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


def write_green_verdict(ours, sha, scope="project"):
    """Give push_own.sh the verdict it now requires, for the tree of one sha.

    Since 2026-10-02 the retry no longer rests on having watched an exit code: it reads the status
    file scripts/check.sh writes and requires it to name the tree of the sha being pushed. So a test
    of the retry has to supply one, and supplying it for the RIGHT tree is part of what is tested.

    Call this while the working tree IS at `sha`, which is how a real check comes to be green: the
    verdict's two tree readings must agree, and they are taken from the working tree. The log goes
    under .git so writing it does not itself move the tree.
    """
    tree = git(ours, "rev-parse", f"{sha}^{{tree}}").stdout.strip()
    log = ours / ".git" / "verdict-pytest.log"
    log.write_text("3758 passed, 1 skipped in 400.11s\n")
    payload = verdict.write_status(
        ours,
        exit_code=0,
        tree_begin=tree,
        scope=scope,
        scope_files=[],
        started_at=time.time(),
        pytest_log=log,
    )
    assert payload["verdict"] == "green", payload
    return tree


def run_push(ours, sha):
    # PYTHONPATH, because push_own.sh reads the verdict with `python3 -m genomeos.verdict` and these
    # repositories are synthetic: in the real checkout the package is simply there beside the script.
    env = dict(os.environ)
    env["PYTHONPATH"] = REPO_ROOT + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run(
        ["bash", SCRIPT, "--sha", sha], cwd=ours, env=env, capture_output=True, text=True, timeout=120
    )


def test_a_remote_that_moved_to_an_ancestor_is_retried_as_a_fast_forward(world):
    bare, ours, peer = world
    mid = commit(ours, "mid.txt")  # the commit the peer will push from under us
    tip = commit(ours, "tip.txt")  # what we are pushing
    write_green_verdict(ours, tip)  # the check passed on exactly this tree, in writing
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


def test_a_green_verdict_for_another_tree_does_not_let_the_retry_through(world):
    """The nine-minute hazard, at the place it would actually do harm.

    The check passes, and then the tree moves: a peer commits, or this session commits again, and the
    sha now being pushed carries different content. The verdict sitting in the directory is green,
    complete and honest -- and it is about something else. The retry skips the check, so if that
    verdict were accepted the push would carry a tree nothing was ever run on.

    Reproduced the way it happens rather than by writing a file by hand: the verdict is taken while
    the working tree is at `mid`, and only then is `tip` committed on top.
    """
    bare, ours, peer = world
    mid = commit(ours, "mid.txt")
    stale_tree = write_green_verdict(ours, mid)  # green, and about mid, which is about to be stale
    tip = commit(ours, "tip.txt")  # the tree moves, as it does in a checkout several sessions share
    assert git(ours, "rev-parse", f"{tip}^{{tree}}").stdout.strip() != stale_tree
    git(peer, "fetch", "-q", "ours", "dev")
    git(peer, "update-ref", "refs/heads/dev", mid)
    arm_hook(ours, peer, "refs/heads/dev")

    r = run_push(ours, tip)
    # the race really happened and the fast-forward really was available, so this cannot pass vacuously
    assert CANNOT_LOCK in r.stdout, (r.stdout, r.stderr)
    assert "an ancestor of" in r.stdout, (r.stdout, r.stderr)
    # and the retry was refused on the verdict alone
    assert r.returncode == 1, (r.stdout, r.stderr)
    assert "REFUSED" in r.stderr, (r.stdout, r.stderr)
    assert stale_tree in r.stderr, "the refusal must name the tree the verdict actually judged"
    # nothing reached the remote: it is still at the peer's commit
    assert git(bare, "rev-parse", "refs/heads/dev").stdout.strip() == mid


def test_no_verdict_at_all_stops_the_retry(world):
    """A missing verdict is not a pass, and the retry is the one path with no check behind it."""
    bare, ours, peer = world
    mid = commit(ours, "mid.txt")
    tip = commit(ours, "tip.txt")
    git(peer, "fetch", "-q", "ours", "dev")
    git(peer, "update-ref", "refs/heads/dev", mid)
    arm_hook(ours, peer, "refs/heads/dev")
    r = run_push(ours, tip)
    assert CANNOT_LOCK in r.stdout, (r.stdout, r.stderr)
    assert "an ancestor of" in r.stdout, (r.stdout, r.stderr)
    assert r.returncode == 1, (r.stdout, r.stderr)
    assert "no check verdict exists" in r.stderr, (r.stdout, r.stderr)
    assert git(bare, "rev-parse", "refs/heads/dev").stdout.strip() == mid


def test_a_skipped_check_is_refused_before_the_push_when_there_is_no_verdict(world):
    """GENOMEOS_SKIP_CHECK means nothing checks the push, so the verdict must already exist."""
    bare, ours, peer = world
    base = git(ours, "rev-parse", "refs/heads/dev").stdout.strip()
    tip = commit(ours, "tip.txt")
    env = dict(os.environ)
    env["PYTHONPATH"] = REPO_ROOT + os.pathsep + env.get("PYTHONPATH", "")
    env["GENOMEOS_SKIP_CHECK"] = "1"
    r = subprocess.run(
        ["bash", SCRIPT, "--sha", tip], cwd=ours, env=env, capture_output=True, text=True, timeout=120
    )
    assert r.returncode == 1, (r.stdout, r.stderr)
    assert "no check verdict exists" in r.stderr, (r.stdout, r.stderr)
    # refused BEFORE the push, not after it: the remote never moved
    assert git(bare, "rev-parse", "refs/heads/dev").stdout.strip() == base

    # ... and with the verdict in writing, the same push goes through
    write_green_verdict(ours, tip)
    r = subprocess.run(
        ["bash", SCRIPT, "--sha", tip], cwd=ours, env=env, capture_output=True, text=True, timeout=120
    )
    assert r.returncode == 0, (r.stdout, r.stderr)
    assert "the status file's verdict is for this exact tree" in r.stdout, (r.stdout, r.stderr)
    assert git(bare, "rev-parse", "refs/heads/dev").stdout.strip() == tip


# --- the outcome, not the exit code of the thing that was supposed to produce it ------------------
#
# On 2026-10-03 push_own.sh was reported as EXITING 0 ON A FAILED PUSH, twice at the terminal: git
# printed `error: failed to push some refs`, the script printed its own correct diagnosis, and the
# caller was handed 0. Reproduced against a local throwaway bare remote, the script exits 1 on that
# path and always has -- `test_a_push_refused_by_the_check_exits_non_zero` below is that path, and it
# passed before anything was changed. The 0 came from the INVOCATION: `push_own.sh 2>&1 | tail -N`
# hands back tail's status and `push_own.sh ... &` hands back the launcher's, both measured at 0 over
# the same red, with the remote unmoved. No script can defend its own status against either, which is
# why the rest of these tests do not try to: they make the CLAIM a fact about the remote, so a
# swallowed status is no longer the only thing standing between a reader and the truth.
#
# None of these skips, so the standing invariant about skip guards keyed on tracked paths does not
# arise; each one is planted, and each was shown to fail against the script as it stood.


def arm_red_check(ours):
    """A pre-push hook that refuses exactly as the real one did in the reported transcript.

    A tree with no green verdict is the cheapest faithful reproduction there is: it needs no network
    and none of the real check's nine minutes, and the lines it prints are the lines that were read.
    """
    hook = ours / ".git" / "hooks" / "pre-push"
    hook.write_text(
        "#!/bin/sh\n"
        '[ "${GENOMEOS_SKIP_CHECK:-0}" = "1" ] && exit 0\n'
        'echo "pre-push: REFUSED: no green verdict exists for the tree being pushed." >&2\n'
        'echo "pre-push: red; fix and push again (see CONTRIBUTING.md)"\n'
        "exit 1\n"
    )
    hook.chmod(0o755)


def test_a_push_refused_by_the_check_exits_non_zero(world):
    """The reported finding, as its own planted case. Nothing covered this path before."""
    bare, ours, peer = world
    base = git(bare, "rev-parse", "refs/heads/dev").stdout.strip()
    tip = commit(ours, "tip.txt")
    arm_red_check(ours)
    r = run_push(ours, tip)
    assert "failed to push some refs" in r.stdout, (r.stdout, r.stderr)
    assert "not the ref-lock race" in r.stderr, (r.stdout, r.stderr)
    assert r.returncode != 0, ("a failed push must not hand back 0", r.stdout, r.stderr)
    # and no claim was made either: the status is not the only thing a reader has
    assert "is on origin/dev" not in r.stdout.lower(), (r.stdout, r.stderr)
    assert git(bare, "rev-parse", "refs/heads/dev").stdout.strip() == base


def test_a_push_that_exits_zero_without_moving_the_branch_is_refused(world):
    """`git push` returning 0 is not the fact this script exists to establish.

    The push's exit status says the server accepted a ref update when it was asked. It does not say
    the branch carries that sha afterwards, and the two really can part: here a server-side hook moves
    `dev` elsewhere once the push is accepted, which is the shape a peer force-push also has. Before
    2026-10-03 this printed `push_own: <sha> is on origin/dev` and exited 0 with the remote at another
    commit entirely -- measured, not reasoned.
    """
    bare, ours, peer = world
    tip = commit(ours, "tip.txt")
    elsewhere = commit(ours, "elsewhere.txt")
    git(ours, "push", "-q", "origin", f"{elsewhere}:refs/heads/parked")  # the bare repo needs the object
    hook = bare / "hooks" / "post-receive"
    hook.write_text(f"#!/bin/sh\ngit update-ref refs/heads/dev {elsewhere}\nexit 0\n")
    hook.chmod(0o755)

    r = run_push(ours, tip)
    # the push really did succeed, so this cannot pass vacuously
    assert "-> dev" in r.stdout, (r.stdout, r.stderr)
    assert r.returncode != 0, (r.stdout, r.stderr)
    assert "REFUSED to claim the push landed" in r.stderr, (r.stdout, r.stderr)
    assert elsewhere in r.stderr, "the refusal must name what is actually on the remote"
    assert "is on origin/dev" not in r.stdout.lower(), (r.stdout, r.stderr)
    assert git(bare, "rev-parse", "refs/heads/dev").stdout.strip() == elsewhere


def test_a_successful_push_still_exits_zero_and_says_it_read_the_remote_back(world):
    """The other direction, which is not optional: a fix that makes everything fail is not a fix."""
    bare, ours, peer = world
    tip = commit(ours, "tip.txt")
    r = run_push(ours, tip)
    assert r.returncode == 0, (r.stdout, r.stderr)
    assert f"IS on origin/dev: read back from the remote as {tip}" in r.stdout, (r.stdout, r.stderr)
    assert git(bare, "rev-parse", "refs/heads/dev").stdout.strip() == tip


def test_a_retry_that_fails_is_not_reported_as_a_push(world):
    """The retry path used to end on an `echo`, so its claim rested on that echo.

    The race is real here -- a peer moves the remote to an ancestor mid-flight, the verdict names the
    pushed tree, the retry is the fast-forward the script says it is -- and only the retry's own push
    is refused.

    WHAT THIS DOES AND DOES NOT PROVE, because the measurement was taken rather than assumed. Run
    against the version this replaced, the STATUS was already non-zero: `set -e` aborted on the failed
    retry before the trailing `echo` could run, so no credit for a status fix is owed here. What that
    version did NOT do is say anything -- it died silently, with its last printed line still promising
    a fast-forward retry, and a reader of the output alone had no sentence to go on. The two
    assertions below are therefore about the claim and not about `set -e`: a refusal is printed, and
    no line says the sha landed.
    """
    bare, ours, peer = world
    mid = commit(ours, "mid.txt")
    tip = commit(ours, "tip.txt")
    write_green_verdict(ours, tip)
    git(peer, "fetch", "-q", "ours", "dev")
    git(peer, "update-ref", "refs/heads/dev", mid)
    arm_hook(ours, peer, "refs/heads/dev")
    # ... and now make the retry itself fail, which the armed hook cannot do: it exits 0 on the skip
    hook = ours / ".git" / "hooks" / "pre-push"
    hook.write_text(
        hook.read_text().replace(
            '[ "${GENOMEOS_SKIP_CHECK:-0}" = "1" ] && exit 0',
            '[ "${GENOMEOS_SKIP_CHECK:-0}" = "1" ] && { echo "the retry is refused too" >&2; exit 1; }',
        )
    )
    hook.chmod(0o755)

    r = run_push(ours, tip)
    assert CANNOT_LOCK in r.stdout, (r.stdout, r.stderr)
    assert "an ancestor of" in r.stdout, "the retry must really have been attempted"
    assert r.returncode != 0, (r.stdout, r.stderr)
    assert "the retry failed" in r.stderr, (r.stdout, r.stderr)
    assert "is on origin/dev" not in r.stdout.lower(), (r.stdout, r.stderr)
    assert git(bare, "rev-parse", "refs/heads/dev").stdout.strip() == mid
