# SPDX-License-Identifier: AGPL-3.0-or-later
"""The three refusals commit_own.sh makes before it writes a commit.

Five commits in this checkout have carried another lane's message, because the scratchpad is
shared and `-F msg1.txt` reads whatever was there. The check that existed compared the whole
message to HEAD's, and caught none of the four after the first: the reused text was older than
HEAD every time. One test per refusal, and one that a message written for its own commit passes
all three.

The harness builds a throwaway repository with three commits, so a duplicate subject can be
further back than HEAD. Nothing here runs scripts/check_staged.py: these three checks come before
it, and the passing case asserts only that the message guard let it through.
"""

import os
import shutil
import subprocess
import time

import pytest

SCRIPT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts", "commit_own.sh")
NAME_REFUSED = 67
STALE_MESSAGE_REFUSED = 68
DUPLICATE_SUBJECT_REFUSED = 66


def git(repo, *args):
    subprocess.run(("git",) + args, cwd=repo, check=True, capture_output=True)


@pytest.fixture()
def repo(tmp_path):
    r = tmp_path / "checkout"
    r.mkdir()
    git(r, "init", "-q", "-b", "dev")
    git(r, "config", "user.email", "lane@example.invalid")
    git(r, "config", "user.name", "A Lane")
    for n, subject in enumerate(["the first lane's work", "the second lane's work", "the third lane's work"]):
        (r / f"f{n}.txt").write_text(f"{n}\n")
        git(r, "add", f"f{n}.txt")
        git(r, "commit", "-q", "-m", subject)
    (r / "mine.txt").write_text("the work this commit is for\n")
    shutil.copy(SCRIPT, r / "commit_own.sh")
    return r


def run(repo, msg_name, text, *extra, age_seconds=0, paths=("mine.txt",)):
    msg = repo / msg_name
    msg.write_text(text)
    if age_seconds:
        then = time.time() - age_seconds
        os.utime(msg, (then, then))
    cmd = ["bash", "./commit_own.sh", "-F", msg_name, *extra, *paths]
    return subprocess.run(cmd, cwd=repo, capture_output=True, text=True)


def fresh_name(lane="mylane"):
    return f"msg-{lane}-purpose-{int(time.time())}.txt"


def test_a_name_any_lane_would_pick_is_refused(repo):
    r = run(repo, "msg1.txt", "a message for this commit alone\n")
    assert r.returncode == NAME_REFUSED, r.stderr
    assert "not a lane-unique message name" in r.stderr


def test_a_name_without_the_committing_lane_is_refused(repo):
    r = run(repo, fresh_name("otherlane"), "a message for this commit alone\n", "-L", "mylane")
    assert r.returncode == NAME_REFUSED, r.stderr
    assert "does not carry the lane name 'mylane'" in r.stderr


def test_a_message_older_than_the_work_is_refused(repo):
    r = run(repo, fresh_name(), "a message for this commit alone\n", age_seconds=3 * 60 * 60)
    assert r.returncode == STALE_MESSAGE_REFUSED, r.stderr
    assert "before the oldest change it describes" in r.stderr


def test_a_message_half_an_hour_older_than_the_work_is_allowed(repo):
    # writing the message and then rerunning the generator that rewrites a result is normal order
    r = run(repo, fresh_name(), "a message for this commit alone\n", age_seconds=20 * 60)
    assert "REFUSED" not in r.stderr, r.stderr


def test_a_first_line_already_in_history_is_refused_even_when_it_is_not_heads(repo):
    r = run(repo, fresh_name(), "the first lane's work\n\nand a body of its own\n")
    assert r.returncode == DUPLICATE_SUBJECT_REFUSED, r.stderr
    assert "already in the last 200 commits" in r.stderr
    third = subprocess.run(
        ("git", "log", "-1", "--skip=2", "--pretty=%h"), cwd=repo, capture_output=True, text=True, check=True
    ).stdout.strip()
    assert f"{third} already begins with this line" in r.stderr


def test_heads_own_subject_is_still_refused(repo):
    r = run(repo, fresh_name(), "the third lane's work\n")
    assert r.returncode == DUPLICATE_SUBJECT_REFUSED, r.stderr


def test_a_message_written_for_its_own_commit_passes_all_three(repo):
    r = run(repo, fresh_name(), "a message written for this commit and no other\n", "-L", "mylane")
    assert "REFUSED" not in r.stderr, r.stderr
    assert r.returncode not in (NAME_REFUSED, STALE_MESSAGE_REFUSED, DUPLICATE_SUBJECT_REFUSED), r.stderr


@pytest.mark.parametrize(
    "name,text,age",
    [
        ("msg1.txt", "a message for this commit alone\n", 0),
        (None, "a message for this commit alone\n", 3 * 60 * 60),
        (None, "the first lane's work\n", 0),
    ],
)
def test_force_carries_past_each_refusal(repo, name, text, age):
    r = run(repo, name or fresh_name(), text, "--force", age_seconds=age)
    assert r.returncode not in (NAME_REFUSED, STALE_MESSAGE_REFUSED, DUPLICATE_SUBJECT_REFUSED), r.stderr
