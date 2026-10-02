# SPDX-License-Identifier: AGPL-3.0-or-later
"""The stale-base guard, against a repository built to contain the accident it is for."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "check_staged.py"

BASE = """# Grammar

| field | type | meaning |
|---|---|---|
| `fraction` | `number` | a share of what is left when it runs |
"""

PEER_ADDITION = """# Grammar

| field | type | meaning |
|---|---|---|
| `fraction` | `number` | a share of what is left when it runs |
| `share` | `number` | a share of the whole at this decision point, order-free |
"""

OWN_EDIT = """# Grammar

| field | type | meaning |
|---|---|---|
| `fraction` | `number` | a share of what is left when it runs |
| `share` | `number` | a share of the whole at this decision point, order-free |
| `order` | `list` | the sequence a cluster is read in |
"""


def git(repo: Path, *args: str, index: str | None = None) -> str:
    env = dict(os.environ)
    env.pop("GIT_INDEX_FILE", None)
    if index:
        env["GIT_INDEX_FILE"] = index
    done = subprocess.run(["git", *args], cwd=repo, capture_output=True, text=True, env=env, check=True)
    return done.stdout


def run_check(repo: Path, index: str, *extra: str) -> subprocess.CompletedProcess[str]:
    env = dict(os.environ)
    env["GIT_INDEX_FILE"] = index
    return subprocess.run(
        [sys.executable, str(SCRIPT), *extra],
        cwd=repo,
        capture_output=True,
        text=True,
        env=env,
    )


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    """A shared file, a peer's commit adding a feature row to it, and a stale copy to stage."""
    repo = tmp_path / "checkout"
    repo.mkdir()
    git(repo, "init", "-q", "-b", "main")
    git(repo, "config", "user.email", "lane@example.test")
    git(repo, "config", "user.name", "A Lane")
    shared = repo / "GRAMMAR.md"
    shared.write_text(BASE)
    git(repo, "add", "GRAMMAR.md")
    git(repo, "commit", "-qm", "the grammar table")
    shared.write_text(PEER_ADDITION)
    git(repo, "add", "GRAMMAR.md")
    git(repo, "commit", "-qm", "share: the same partition, stated as its numbers")
    return repo


def private_index(repo: Path, tmp_path: Path) -> str:
    index = str(tmp_path / "idx")
    git(repo, "read-tree", "HEAD", index=index)
    return index


def stage_blob(repo: Path, index: str, path: str, content: str) -> None:
    """Put content into the index without touching the working copy, as a staging tool would."""
    tmp = repo / ".stage-tmp"
    tmp.write_text(content)
    sha = git(repo, "hash-object", "-w", str(tmp), index=index).strip()
    tmp.unlink()
    git(repo, "update-index", "--add", "--cacheinfo", f"100644,{sha},{path}", index=index)


def test_stale_copy_is_refused_and_names_the_peer(repo: Path, tmp_path: Path) -> None:
    """The accident of 2026-09-15: staging the copy you read before the peer's commit landed."""
    index = private_index(repo, tmp_path)
    stage_blob(repo, index, "GRAMMAR.md", BASE)  # the pre-share copy, staged whole

    done = run_check(repo, index)

    assert done.returncode == 2, done.stdout
    assert "REFUSED" in done.stdout
    assert "share: the same partition" in done.stdout
    assert "`share`" in done.stdout


def test_an_addition_on_top_of_the_peer_passes(repo: Path, tmp_path: Path) -> None:
    """Adding your own row to the current copy removes nobody's line and must not be refused."""
    index = private_index(repo, tmp_path)
    stage_blob(repo, index, "GRAMMAR.md", OWN_EDIT)

    done = run_check(repo, index)

    assert done.returncode == 0, done.stdout
    assert "takes nothing back out" in done.stdout


def test_force_reports_without_refusing(repo: Path, tmp_path: Path) -> None:
    """A deliberate removal is sometimes right; it costs a flag and still prints what it takes out."""
    index = private_index(repo, tmp_path)
    stage_blob(repo, index, "GRAMMAR.md", BASE)

    done = run_check(repo, index, "--force")

    assert done.returncode == 0, done.stdout
    assert "REFUSED" in done.stdout  # seen, not hidden


def test_a_staged_deletion_of_a_tracked_file_is_reported(repo: Path, tmp_path: Path) -> None:
    """Losing a peer's whole file is the unrecoverable case, so it is never silent."""
    index = private_index(repo, tmp_path)
    git(repo, "update-index", "--force-remove", "GRAMMAR.md", index=index)

    done = run_check(repo, index)

    assert done.returncode == 2
    assert "staged as DELETED" in done.stdout


def test_an_empty_commit_is_refused(repo: Path, tmp_path: Path) -> None:
    """A push that looks like it failed may have landed; committing again produces an empty commit."""
    index = private_index(repo, tmp_path)

    done = run_check(repo, index)

    assert done.returncode == 2
    assert "would be empty" in done.stdout


def test_a_new_file_reverts_nothing(repo: Path, tmp_path: Path) -> None:
    """A path HEAD does not have cannot take anything back out, whatever it contains."""
    index = private_index(repo, tmp_path)
    stage_blob(repo, index, "NOTES.md", BASE)

    done = run_check(repo, index)

    assert done.returncode == 0, done.stdout


def test_old_removals_pass_and_the_window_decides(repo: Path, tmp_path: Path) -> None:
    """The guard is a stale-base detector, not a merge policeman: outside the window it says nothing."""
    index = private_index(repo, tmp_path)
    stage_blob(repo, index, "GRAMMAR.md", BASE)

    refused = run_check(repo, index)
    passed = run_check(repo, index, "--since", "2099-01-01")  # a window with no commit in it

    assert refused.returncode == 2
    assert passed.returncode == 0, passed.stdout


def test_an_import_of_an_untracked_module_is_refused(repo: Path, tmp_path: Path) -> None:
    """38087a1: a test committed without the module it imports passes in the shared tree only."""
    index = private_index(repo, tmp_path)
    stage_blob(repo, index, "tests/test_new.py", "from genomeos.certainty import Certainty\n")

    done = run_check(repo, index)

    assert done.returncode == 2, done.stdout
    assert "genomeos.certainty" in done.stdout
    assert "genomeos/certainty.py" in done.stdout


def test_force_does_not_pass_a_missing_module(repo: Path, tmp_path: Path) -> None:
    """The fix is to stage the module, so the flag that waves a revert through does not apply."""
    index = private_index(repo, tmp_path)
    stage_blob(repo, index, "tests/test_new.py", "import genomeos.certainty\n")

    done = run_check(repo, index, "--force")

    assert done.returncode == 2, done.stdout


def test_an_import_staged_with_its_module_passes(repo: Path, tmp_path: Path) -> None:
    """Module and importer in one commit, including a package's __init__ and a nested name."""
    index = private_index(repo, tmp_path)
    stage_blob(repo, index, "genomeos/__init__.py", "")
    stage_blob(repo, index, "genomeos/certainty.py", "class Certainty: ...\n")
    stage_blob(repo, index, "genomeos/lang/__init__.py", "")
    body = "from genomeos.certainty import Certainty\nfrom genomeos import lang\nimport genomeos.lang\n"
    stage_blob(repo, index, "tests/test_new.py", body)

    done = run_check(repo, index)

    assert done.returncode == 0, done.stdout


# The false positive of 2026-10-02, and the hole its fix must not open.
#
# Git reports a line put inside an `if`/`else` as a removal plus an addition of the same stripped
# content. The guard read only the removal, so the honest shape -- keep the peer's line, branch around
# it -- was refused as a revert of the very line it preserves. That is not merely annoying: every false
# refusal is pressure to reach for `--force`, and `--force` is the flag that can really lose a peer's
# work. The fix subtracts what the same diff adds back, which makes the guard more precise rather than
# weaker, so the tests below have to show BOTH halves: the false refusal gone, and every real removal
# still refused -- including a replacement line that merely looks like the held one.

CODE_BASE = """def register(payload):
    path = RESULT
    return path
"""

# the peer commit, inside the window: it adds exactly the two lines the tests below move or remove
CODE_PEER = """def register(payload):
    path = RESULT
    p = save_result(RESULT, payload)
    stamp = record_stamp(p)
    return path
"""

HELD = "p = save_result(RESULT, payload)"
PEER_SUBJECT = "save_result keeps the name it wrote"


@pytest.fixture
def code_repo(tmp_path: Path) -> Path:
    """A python file whose last commit, minutes old, added the two lines the tests act on."""
    repo = tmp_path / "code"
    repo.mkdir()
    git(repo, "init", "-q", "-b", "main")
    git(repo, "config", "user.email", "peer@example.test")
    git(repo, "config", "user.name", "A Peer")
    source = repo / "reg.py"
    source.write_text(CODE_BASE)
    git(repo, "add", "reg.py")
    git(repo, "commit", "-qm", "the register skeleton")
    source.write_text(CODE_PEER)
    git(repo, "add", "reg.py")
    git(repo, "commit", "-qm", PEER_SUBJECT)
    return repo


def test_a_genuine_removal_of_a_peer_line_is_still_refused(code_repo: Path, tmp_path: Path) -> None:
    """Shown to fire: the line is gone, nothing puts it back, and the refusal reads as it always did."""
    index = private_index(code_repo, tmp_path)
    stage_blob(
        code_repo,
        index,
        "reg.py",
        "def register(payload):\n    path = RESULT\n    stamp = record_stamp(path)\n    return path\n",
    )

    done = run_check(code_repo, index)

    assert done.returncode == 2, done.stdout
    assert "REFUSED: this commit would undo work that is already in HEAD" in done.stdout
    assert PEER_SUBJECT in done.stdout
    assert HELD in done.stdout


def test_force_still_reports_a_genuine_removal_without_refusing(code_repo: Path, tmp_path: Path) -> None:
    """--force is unchanged by the subtraction: it still prints the finding and still exits 0."""
    index = private_index(code_repo, tmp_path)
    stage_blob(
        code_repo,
        index,
        "reg.py",
        "def register(payload):\n    path = RESULT\n    stamp = record_stamp(path)\n    return path\n",
    )

    done = run_check(code_repo, index, "--force")

    assert done.returncode == 0, done.stdout
    assert "REFUSED" in done.stdout
    assert HELD in done.stdout


def test_a_reindented_peer_line_passes(code_repo: Path, tmp_path: Path) -> None:
    """The false positive: moving the held line into a branch indents it and preserves it entirely."""
    index = private_index(code_repo, tmp_path)
    stage_blob(
        code_repo,
        index,
        "reg.py",
        "def register(payload):\n"
        "    path = RESULT\n"
        '    if payload.get("name"):\n'
        f"        {HELD}\n"
        "    else:\n"
        f"        {HELD}\n"
        "    stamp = record_stamp(p)\n"
        "    return path\n",
    )

    done = run_check(code_repo, index)

    assert done.returncode == 0, done.stdout
    assert "takes nothing back out" in done.stdout


def test_a_peer_line_moved_elsewhere_at_a_new_indentation_passes(code_repo: Path, tmp_path: Path) -> None:
    """Content, not position: the line leaves one function for another and is still in the file."""
    index = private_index(code_repo, tmp_path)
    stage_blob(
        code_repo,
        index,
        "reg.py",
        "def stamped(p):\n"
        "    if p:\n"
        "        stamp = record_stamp(p)\n"
        "    return None\n"
        "\n"
        "\n"
        "def register(payload):\n"
        "    path = RESULT\n"
        f"    {HELD}\n"
        "    return path\n",
    )

    done = run_check(code_repo, index)

    assert done.returncode == 0, done.stdout
    assert "takes nothing back out" in done.stdout


def test_a_replacement_differing_only_in_a_variable_name_is_still_refused(
    code_repo: Path, tmp_path: Path
) -> None:
    """The near-miss: if the subtraction matched loosely this would pass, and the guard would be gone."""
    index = private_index(code_repo, tmp_path)
    stage_blob(
        code_repo,
        index,
        "reg.py",
        "def register(payload):\n"
        "    path = RESULT\n"
        "    p = save_result(RESULT, data)\n"
        "    stamp = record_stamp(p)\n"
        "    return path\n",
    )

    done = run_check(code_repo, index)

    assert done.returncode == 2, done.stdout
    assert PEER_SUBJECT in done.stdout
    assert HELD in done.stdout


def test_a_replacement_with_a_comment_appended_is_still_refused(code_repo: Path, tmp_path: Path) -> None:
    """The other near-miss: the held text is a prefix of the new line, and a prefix is not the line."""
    index = private_index(code_repo, tmp_path)
    stage_blob(
        code_repo,
        index,
        "reg.py",
        "def register(payload):\n"
        "    path = RESULT\n"
        f"    {HELD}  # the chosen name\n"
        "    stamp = record_stamp(p)\n"
        "    return path\n",
    )

    done = run_check(code_repo, index)

    assert done.returncode == 2, done.stdout
    assert PEER_SUBJECT in done.stdout


def test_a_reindentation_does_not_cover_a_real_removal_in_the_same_diff(
    code_repo: Path, tmp_path: Path
) -> None:
    """Per line, not per diff: one held line is re-indented, the other is deleted, and only that one
    is named. A subtraction that waved the whole path through would report nothing here."""
    index = private_index(code_repo, tmp_path)
    stage_blob(
        code_repo,
        index,
        "reg.py",
        "def register(payload):\n"
        "    path = RESULT\n"
        '    if payload.get("name"):\n'
        f"        {HELD}\n"
        "    else:\n"
        f"        {HELD}\n"
        "    return path\n",
    )

    done = run_check(code_repo, index)

    assert done.returncode == 2, done.stdout
    assert "removes 1 lines" in done.stdout
    assert "stamp = record_stamp(p)" in done.stdout
    assert HELD not in done.stdout  # preserved, so never named as taken out


def test_the_subtraction_is_per_path_so_a_line_moved_to_another_file_is_still_refused(
    code_repo: Path, tmp_path: Path
) -> None:
    """The scope of the change, pinned. Each path is judged by its own diff, so content that leaves
    this file for a different one is still a removal from this file and still costs a look. That is
    the conservative reading and it is deliberate: a cross-file move is exactly the shape where the
    author should confirm the destination is really the same content, and `--force` says so in one
    word. The subtraction excuses a removal only where the same path puts the same text back."""
    index = private_index(code_repo, tmp_path)
    stage_blob(
        code_repo,
        index,
        "reg.py",
        "def register(payload):\n    path = RESULT\n    stamp = record_stamp(p)\n    return path\n",
    )
    stage_blob(code_repo, index, "helper.py", f"def helped(payload):\n    {HELD}\n    return p\n")

    done = run_check(code_repo, index)

    assert done.returncode == 2, done.stdout
    assert "reg.py: removes" in done.stdout
    assert PEER_SUBJECT in done.stdout
