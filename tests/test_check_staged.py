# SPDX-License-Identifier: AGPL-3.0-or-later
"""The stale-base guard, against a repository built to contain the accident it is for."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "check_staged.py"
# the counterfactual below cuts the guard between the marker lines the guard itself names,
# so the markers cannot drift from the block they delimit
sys.path.insert(0, str(SCRIPT.parent))

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


# --------------------------------------------------------------------------------------------------
# The result-reproducibility conditions (2026-10-02). A result stamps the revision it was written at;
# when the code that wrote it was uncommitted, that sha names a tree without the code and NO COMMIT
# REPRODUCES THE FILE. Two results went in that way the same night, and a third went in with no
# manifest at all and was caught fifteen minutes into a push, about a tree that no longer existed.
# --------------------------------------------------------------------------------------------------

#: The shape of data/results/response_map_increment3_count.json at 3360c49, field for field: the sha it
#: stamped, the dirty tree, and its own script named as uncommitted.
AT_3360C49 = {
    "git_sha": "727c9363dfcd23428e85f00da66fa99f3a818569",
    "dirty": True,
    "own_uncommitted_code": ["scripts/response_map_increment3_count.py"],
    "own_code_is_committed": False,
    "foreign_uncommitted_code": [],
    "foreign_uncommitted_code_on_the_counting_path": [],
}
#: A tree with other lanes' work uncommitted but none of it on the counting path: the committed shape
#: of the same result after its repair, and the near-miss that must still commit.
CLEAN_BUT_DIRTY_TREE = {
    "git_sha": "7509b1ada433ed33712fa66ffeeb2a95bc10f620",
    "dirty": True,
    "own_uncommitted_code": [],
    "own_code_is_committed": True,
    "foreign_uncommitted_code": ["genomeos/attribution/astrorun.py", "tests/test_astrorun.py"],
    "foreign_uncommitted_code_on_the_counting_path": [],
}
#: The seven-result shape: a peer's module on the import closure that computed the number.
FOREIGN_ON_THE_PATH = dict(
    CLEAN_BUT_DIRTY_TREE,
    foreign_uncommitted_code=["genomeos/attribution/astrorun.py"],
    foreign_uncommitted_code_on_the_counting_path=["genomeos/attribution/astrorun.py"],
)


def result_json(cleanliness: dict | None, *, manifest: bool = True) -> str:
    """A result as the registry holds one: a number, and a manifest unless the caller wants none."""
    body: dict = {"loci": 110, "delivered": 1322}
    if manifest:
        body["result_manifest"] = {
            "sources": [{"accession": "ENCSR000XXX", "version": "v1"}],
            "complete": True,
            **({"code_cleanliness": dict(cleanliness)} if cleanliness is not None else {}),
        }
    return json.dumps(body, indent=1)


@pytest.fixture
def registry_repo(tmp_path: Path) -> Path:
    """A checkout with a committed legacy allowlist, as the real one has beside data/results."""
    repo = tmp_path / "registry"
    repo.mkdir()
    git(repo, "init", "-q", "-b", "main")
    git(repo, "config", "user.email", "lane@example.test")
    git(repo, "config", "user.name", "A Lane")
    (repo / "data" / "results").mkdir(parents=True)
    (repo / "data" / "results_legacy.txt").write_text(
        "# name<TAB>source\nold_reader_chr1\ttracked\nabundance_gate\tignored-local\n"
    )
    (repo / "README.md").write_text("# a checkout\n")
    git(repo, "add", "data/results_legacy.txt", "README.md")
    git(repo, "commit", "-qm", "the legacy allowlist")
    return repo


def test_a_result_whose_own_code_was_uncommitted_is_refused(registry_repo: Path, tmp_path: Path) -> None:
    """3360c49, exactly: the stamp names 727c936, a tree without the script named beside it. This is
    the only fault here that makes a number permanently unverifiable instead of merely awkward."""
    index = private_index(registry_repo, tmp_path)
    stage_blob(registry_repo, index, "data/results/increment3_count.json", result_json(AT_3360C49))

    done = run_check(registry_repo, index)

    assert done.returncode == 2, done.stdout
    assert "no commit reproduces" in done.stdout
    assert "data/results/increment3_count.json" in done.stdout  # the file
    assert "own_code_is_committed" in done.stdout  # the field
    assert "727c936" in done.stdout
    assert "scripts/response_map_increment3_count.py" in done.stdout


def test_a_peers_uncommitted_module_on_the_counting_path_is_refused(
    registry_repo: Path, tmp_path: Path
) -> None:
    """The seven-result shape: the number was computed through a file the stamped tree does not have."""
    index = private_index(registry_repo, tmp_path)
    stage_blob(registry_repo, index, "data/results/astro_count.json", result_json(FOREIGN_ON_THE_PATH))

    done = run_check(registry_repo, index)

    assert done.returncode == 2, done.stdout
    assert "data/results/astro_count.json" in done.stdout
    assert "foreign_uncommitted_code_on_the_counting_path" in done.stdout
    assert "genomeos/attribution/astrorun.py" in done.stdout


def test_a_result_with_no_manifest_is_refused_at_the_commit(registry_repo: Path, tmp_path: Path) -> None:
    """1f880b2: astroreg_request_plan.json, committed with no manifest. The push-time test did catch
    it -- fifteen minutes later, in a clean worktree, about a tree the repair had already replaced.
    Here the same question is asked of the index, where the answer can still stop something."""
    index = private_index(registry_repo, tmp_path)
    stage_blob(
        registry_repo,
        index,
        "data/results/astroreg_request_plan.json",
        result_json(None, manifest=False),
    )

    done = run_check(registry_repo, index)

    assert done.returncode == 2, done.stdout
    assert "data/results/astroreg_request_plan.json" in done.stdout
    assert "result_manifest" in done.stdout  # the field
    assert "legacy allowlist" in done.stdout


def test_a_clean_result_commits(registry_repo: Path, tmp_path: Path) -> None:
    """The near-miss that matters most: a guard that refuses the honest case is one people route round."""
    index = private_index(registry_repo, tmp_path)
    clean = dict(CLEAN_BUT_DIRTY_TREE, dirty=False, foreign_uncommitted_code=[])
    stage_blob(registry_repo, index, "data/results/good.json", result_json(clean))

    done = run_check(registry_repo, index)

    assert done.returncode == 0, done.stdout
    assert "takes nothing back out" in done.stdout


def test_an_empty_foreign_list_in_a_dirty_tree_commits(registry_repo: Path, tmp_path: Path) -> None:
    """Several lanes share this checkout, so a tree is nearly always dirty with somebody's work. The
    condition is the counting path, not the tree: uncommitted peer code that no import reached cannot
    have entered the number, and refusing it would refuse almost every honest result."""
    index = private_index(registry_repo, tmp_path)
    stage_blob(registry_repo, index, "data/results/good.json", result_json(CLEAN_BUT_DIRTY_TREE))

    done = run_check(registry_repo, index)

    assert done.returncode == 0, done.stdout
    assert "takes nothing back out" in done.stdout


def test_a_legacy_allowlisted_name_with_no_manifest_commits(registry_repo: Path, tmp_path: Path) -> None:
    """955 results were written before the contract. The exemption is the committed list that
    `save_result` reads, and this reads the same list through the same function."""
    index = private_index(registry_repo, tmp_path)
    stage_blob(registry_repo, index, "data/results/old_reader_chr1.json", result_json(None, manifest=False))

    done = run_check(registry_repo, index)

    assert done.returncode == 0, done.stdout


def test_the_allowlist_is_the_one_save_result_reads(registry_repo: Path, tmp_path: Path) -> None:
    """Not a copy of it: a name absent from the list is refused and the same name added to the list is
    not, so the file decides. Two copies of 955 names would drift the first time anybody touched one."""
    index = private_index(registry_repo, tmp_path)
    stage_blob(registry_repo, index, "data/results/only_in_the_list.json", result_json(None, manifest=False))
    assert run_check(registry_repo, index).returncode == 2

    listed = registry_repo / "data" / "results_legacy.txt"
    listed.write_text(listed.read_text() + "only_in_the_list\ttracked\n")

    assert run_check(registry_repo, index).returncode == 0


def test_a_json_outside_the_registry_is_untouched(registry_repo: Path, tmp_path: Path) -> None:
    """The contract is the result registry's. A config, a fixture or a nested file is not a result."""
    index = private_index(registry_repo, tmp_path)
    body = result_json(None, manifest=False)
    stage_blob(registry_repo, index, "data/knowledge/panel.json", body)
    stage_blob(registry_repo, index, "tests/data/case.json", body)
    stage_blob(registry_repo, index, "data/results/held/inner.json", body)

    done = run_check(registry_repo, index)

    assert done.returncode == 0, done.stdout


def test_force_reports_an_unreproducible_result_without_refusing(registry_repo: Path, tmp_path: Path) -> None:
    """`--force` behaves here as it does for every other condition: it prints and does not refuse."""
    index = private_index(registry_repo, tmp_path)
    stage_blob(registry_repo, index, "data/results/increment3_count.json", result_json(AT_3360C49))

    done = run_check(registry_repo, index, "--force")

    assert done.returncode == 0, done.stdout
    assert "own_code_is_committed" in done.stdout  # seen, not hidden


def test_the_staged_blob_decides_and_not_the_working_copy(registry_repo: Path, tmp_path: Path) -> None:
    """What this commit publishes is the index. A repaired file on disk does not excuse a staged one."""
    index = private_index(registry_repo, tmp_path)
    stage_blob(registry_repo, index, "data/results/increment3_count.json", result_json(AT_3360C49))
    (registry_repo / "data" / "results" / "increment3_count.json").write_text(
        result_json(CLEAN_BUT_DIRTY_TREE)
    )

    done = run_check(registry_repo, index)

    assert done.returncode == 2, done.stdout
    assert "own_code_is_committed" in done.stdout


def without_the_conditions(tmp_path: Path) -> Path:
    """A copy of the guard with the three conditions cut out, between the markers the script names."""
    import check_staged as cs

    lines = SCRIPT.read_text().splitlines(keepends=True)
    begin = next(i for i, line in enumerate(lines) if line.rstrip("\n") == cs.CONDITIONS_BEGIN)
    end = next(i for i, line in enumerate(lines) if line.rstrip("\n") == cs.CONDITIONS_END)
    assert begin < end
    cut = tmp_path / "check_staged_without.py"
    cut.write_text("".join(lines[:begin] + lines[end + 1 :]))
    assert "own_code_is_committed` is false" not in cut.read_text()
    return cut


@pytest.mark.parametrize(
    ("name", "cleanliness", "manifest", "shape"),
    [
        ("increment3_count", AT_3360C49, True, "3360c49"),
        ("astro_count", FOREIGN_ON_THE_PATH, True, "the seven-result shape"),
        ("astroreg_request_plan", None, False, "1f880b2"),
    ],
)
def test_without_the_conditions_the_harm_happens(
    registry_repo: Path, tmp_path: Path, name: str, cleanliness: dict | None, manifest: bool, shape: str
) -> None:
    """Credited with preventing something, shown to prevent it by removal.

    The rule here is that a guard is demoted to a diagnosis unless the harm happens when it is taken
    away -- on 2026-10-02 that test demoted a guard in scripts/push_own.sh, because git's own refusal
    was doing the work. These three are protections: with the block cut out, each of the two shapes
    that went in that night, and the peer-code shape, is committed without a word.
    """
    index = private_index(registry_repo, tmp_path)
    stage_blob(registry_repo, index, f"data/results/{name}.json", result_json(cleanliness, manifest=manifest))

    refused = run_check(registry_repo, index)
    assert refused.returncode == 2, f"{shape} must be refused: {refused.stdout}"

    env = dict(os.environ)
    env["GIT_INDEX_FILE"] = index
    allowed = subprocess.run(
        [sys.executable, str(without_the_conditions(tmp_path))],
        cwd=registry_repo,
        capture_output=True,
        text=True,
        env=env,
    )

    assert allowed.returncode == 0, allowed.stdout + allowed.stderr
    assert "takes nothing back out" in allowed.stdout  # through, silently
    assert name not in allowed.stdout


# --------------------------------------------------------------- a file another lane is working on
#
# 2026-10-03: a lane staged a whole test file carrying another lane's in-flight tests and committed
# it under its own message. The guard was clean, exit 0, no --force, and HEAD went red on NameError
# because the tests went in without their imports. The holder had declared the file on the work
# board, and the edit-time hook says so; the commit never asked.

PEER_TASK = "the kinetics rows and their fixture, half written"


def hold(repo: Path, who: str, files: list[str], *, minutes: float = 3, state: str = "working") -> None:
    """One work-board entry, in the shape `genomeos work start` writes."""
    d = repo / "data" / "work"
    d.mkdir(parents=True, exist_ok=True)
    stamp = time.time() - minutes * 60
    (d / f"{who}.json").write_text(
        json.dumps(
            {
                "who": who,
                "task": PEER_TASK,
                "area": None,
                "files": files,
                "next": None,
                "state": state,
                "note": None,
                "started": stamp,
                "updated": stamp,
                "finished": stamp if state == "done" else None,
            }
        )
    )


def without_block(tmp_path: Path, begin: str, end: str, name: str) -> Path:
    """A copy of the guard with one marked block cut out, between the markers the script names."""
    lines = SCRIPT.read_text().splitlines(keepends=True)
    first = next(i for i, line in enumerate(lines) if line.rstrip("\n") == begin)
    last = next(i for i, line in enumerate(lines) if line.rstrip("\n") == end)
    assert first < last
    cut = tmp_path / name
    cut.write_text("".join(lines[:first] + lines[last + 1 :]))
    return cut


def run_cut(repo: Path, index: str, cut: Path, *extra: str) -> subprocess.CompletedProcess[str]:
    env = dict(os.environ)
    env["GIT_INDEX_FILE"] = index
    return subprocess.run(
        [sys.executable, str(cut), *extra], cwd=repo, capture_output=True, text=True, env=env
    )


def test_a_file_a_live_lane_holds_is_refused(repo: Path, tmp_path: Path) -> None:
    """The commit removes nothing, so only the hold can refuse it -- and it must."""
    hold(repo, "lane-kinetics", ["GRAMMAR.md"])
    index = private_index(repo, tmp_path)
    stage_blob(repo, index, "GRAMMAR.md", OWN_EDIT)

    done = run_check(repo, index, "--lane", "lane-mine")

    assert done.returncode == 2, done.stdout
    assert "another lane is working on" in done.stdout
    assert "lane-kinetics" in done.stdout
    assert PEER_TASK in done.stdout
    assert "3 min ago" in done.stdout
    assert "hand the file to the coordinator" in done.stdout
    assert "stage only your own hunk" in done.stdout


def test_a_space_separated_files_field_still_holds_the_file(repo: Path, tmp_path: Path) -> None:
    """Entries on this board store four files in one string, because --files split only on commas."""
    hold(repo, "lane-kinetics", ["docs/x.md GRAMMAR.md scripts/y.py"])
    index = private_index(repo, tmp_path)
    stage_blob(repo, index, "GRAMMAR.md", OWN_EDIT)

    done = run_check(repo, index, "--lane", "lane-mine")

    assert done.returncode == 2, done.stdout
    assert "lane-kinetics" in done.stdout


def test_force_does_not_pass_another_lanes_held_file(repo: Path, tmp_path: Path) -> None:
    """Lanes never force in this project: the holder is a person to ask, not a flag."""
    hold(repo, "lane-kinetics", ["GRAMMAR.md"])
    index = private_index(repo, tmp_path)
    stage_blob(repo, index, "GRAMMAR.md", OWN_EDIT)

    done = run_check(repo, index, "--force", "--lane", "lane-mine")

    assert done.returncode == 2, done.stdout
    assert "lane-kinetics" in done.stdout


def test_your_own_hold_passes_when_the_lane_is_named(repo: Path, tmp_path: Path) -> None:
    hold(repo, "lane-mine", ["GRAMMAR.md"])
    index = private_index(repo, tmp_path)
    stage_blob(repo, index, "GRAMMAR.md", OWN_EDIT)

    done = run_check(repo, index, "--lane", "lane-mine")

    assert done.returncode == 0, done.stdout
    assert "takes nothing back out" in done.stdout


def test_your_own_hold_passes_through_the_message_file_name(repo: Path, tmp_path: Path) -> None:
    """commit_own.sh already requires the message name to carry the lane, so it identifies you."""
    hold(repo, "lane-mine", ["GRAMMAR.md"])
    index = private_index(repo, tmp_path)
    stage_blob(repo, index, "GRAMMAR.md", OWN_EDIT)

    done = run_check(repo, index, "--message-name", "msg-lane-mine-grammar-1791040000.txt")

    assert done.returncode == 0, done.stdout


@pytest.mark.parametrize(
    ("minutes", "state", "why"),
    [
        (7 * 60, "working", "stale: six hours untouched and the board stops calling it current"),
        (3 * 24 * 60, "working", "abandoned: the session that wrote it is gone"),
        (5, "done", "finished, and a finished entry holds nothing"),
    ],
)
def test_a_hold_that_is_not_live_passes(
    repo: Path, tmp_path: Path, minutes: float, state: str, why: str
) -> None:
    hold(repo, "lane-kinetics", ["GRAMMAR.md"], minutes=minutes, state=state)
    index = private_index(repo, tmp_path)
    stage_blob(repo, index, "GRAMMAR.md", OWN_EDIT)

    done = run_check(repo, index, "--lane", "lane-mine")

    assert done.returncode == 0, f"{why}: {done.stdout}"
    assert "lane-kinetics" not in done.stdout


def test_a_directory_hold_does_not_refuse(repo: Path, tmp_path: Path) -> None:
    """Deliberate, and the hole this leaves: entries hold `tests/` and `scripts/` for a day at a
    time, so refusing on a region would stop every other lane committing any test or any script."""
    hold(repo, "lane-kinetics", ["docs/"])
    index = private_index(repo, tmp_path)
    stage_blob(repo, index, "docs/note.md", "a new file of my own\n")

    done = run_check(repo, index, "--lane", "lane-mine")

    assert done.returncode == 0, done.stdout


def test_without_the_live_hold_condition_the_peers_file_goes_in_silently(repo: Path, tmp_path: Path) -> None:
    """Credited with preventing something, shown to prevent it by removal."""
    import check_staged as cs

    hold(repo, "lane-kinetics", ["GRAMMAR.md"])
    index = private_index(repo, tmp_path)
    stage_blob(repo, index, "GRAMMAR.md", OWN_EDIT)

    refused = run_check(repo, index, "--lane", "lane-mine")
    assert refused.returncode == 2, refused.stdout

    cut = without_block(tmp_path, cs.HELD_BEGIN, cs.HELD_END, "check_staged_without_hold.py")
    allowed = run_cut(repo, index, cut, "--lane", "lane-mine")

    assert allowed.returncode == 0, allowed.stdout + allowed.stderr
    assert "takes nothing back out" in allowed.stdout
    assert "lane-kinetics" not in allowed.stdout


# ------------------------------------------------------------- the whole-copy notice, two sources
#
# The hook's notice fired only for a hardcoded list of eight paths, and a test file is not on it,
# which is why nothing was said on 2026-10-03. Being held by a live lane is now the second way a
# file counts as shared; the list stays, because a file can be shared with nobody holding it.

README = "# GenomeOS\n\nthe engine and the application\n"


@pytest.fixture
def shared_repo(tmp_path: Path) -> Path:
    """A listed shared file and an unlisted one, both committed, both dirty in the worktree."""
    repo = tmp_path / "shared"
    repo.mkdir()
    git(repo, "init", "-q", "-b", "main")
    git(repo, "config", "user.email", "lane@example.test")
    git(repo, "config", "user.name", "A Lane")
    (repo / "README.md").write_text(README)
    (repo / "tests").mkdir()
    (repo / "tests" / "test_kinetics.py").write_text("def test_one():\n    assert True\n")
    git(repo, "add", "README.md", "tests/test_kinetics.py")
    git(repo, "commit", "-qm", "the readme and a test")
    return repo


def stage_worktree(repo: Path, index: str, path: str, content: str) -> None:
    """Write the working copy and stage it whole, which is what commit_own.sh does with a path."""
    f = repo / path
    f.parent.mkdir(parents=True, exist_ok=True)
    f.write_text(content)
    git(repo, "add", "--", path, index=index)


def test_a_listed_shared_file_staged_whole_is_noticed(shared_repo: Path, tmp_path: Path) -> None:
    """The static list keeps working: it is a second source now, not a deleted one."""
    index = private_index(shared_repo, tmp_path)
    stage_worktree(shared_repo, index, "README.md", README + "\na line of my own\n")

    done = run_check(shared_repo, index, "--lane", "lane-mine")

    assert done.returncode == 0, done.stdout
    assert "NOTICE" in done.stdout
    assert "README.md is a shared file" in done.stdout
    assert "whole working copy" in done.stdout


def test_a_held_file_that_is_on_no_list_is_noticed(shared_repo: Path, tmp_path: Path) -> None:
    """The file that went in silently on 2026-10-03 was a test file, on nobody's list."""
    hold(shared_repo, "lane-kinetics", ["tests/test_kinetics.py"])
    index = private_index(shared_repo, tmp_path)
    stage_worktree(
        shared_repo,
        index,
        "tests/test_kinetics.py",
        "def test_one():\n    assert True\n\n\ndef test_two():\n    assert True\n",
    )

    done = run_check(shared_repo, index, "--lane", "lane-mine")

    assert "NOTICE" in done.stdout
    assert "tests/test_kinetics.py is held on the work board by lane-kinetics" in done.stdout
    assert "uncommitted hunk(s)" in done.stdout
    assert done.returncode == 2, done.stdout  # and the refusal of finding 1 as well


def test_a_directory_hold_is_noticed_even_though_it_does_not_refuse(
    shared_repo: Path, tmp_path: Path
) -> None:
    """A region claim is the case the refusal deliberately leaves alone, so the notice carries it."""
    hold(shared_repo, "lane-kinetics", ["tests/"])
    index = private_index(shared_repo, tmp_path)
    stage_worktree(
        shared_repo,
        index,
        "tests/test_kinetics.py",
        "def test_one():\n    assert True\n\n\ndef test_two():\n    assert True\n",
    )

    done = run_check(shared_repo, index, "--lane", "lane-mine")

    assert done.returncode == 0, done.stdout
    assert "inside tests/, held on the work board by lane-kinetics" in done.stdout


def test_a_file_neither_listed_nor_held_is_not_noticed(shared_repo: Path, tmp_path: Path) -> None:
    index = private_index(shared_repo, tmp_path)
    stage_worktree(shared_repo, index, "tests/test_kinetics.py", "def test_one():\n    assert 1\n")

    done = run_check(shared_repo, index, "--lane", "lane-mine")

    assert done.returncode == 0, done.stdout
    assert "NOTICE" not in done.stdout


def test_staging_your_own_hunk_instead_says_nothing(shared_repo: Path, tmp_path: Path) -> None:
    """The notice is about taking the WHOLE working copy, so a partial staging silences it."""
    hold(shared_repo, "lane-mine", ["README.md"])
    index = private_index(shared_repo, tmp_path)
    stage_blob(shared_repo, index, "README.md", README + "\nmy own line\n")
    (shared_repo / "README.md").write_text(README + "\nmy own line\na peer's line\n")

    done = run_check(shared_repo, index, "--lane", "lane-mine")

    assert done.returncode == 0, done.stdout
    assert "NOTICE" not in done.stdout


def test_without_the_whole_copy_condition_a_held_file_is_not_noticed(
    shared_repo: Path, tmp_path: Path
) -> None:
    """Credited with preventing something, shown by removal: with the block cut, nothing is said."""
    import check_staged as cs

    hold(shared_repo, "lane-kinetics", ["tests/"])
    index = private_index(shared_repo, tmp_path)
    stage_worktree(
        shared_repo,
        index,
        "tests/test_kinetics.py",
        "def test_one():\n    assert True\n\n\ndef test_two():\n    assert True\n",
    )

    said = run_check(shared_repo, index, "--lane", "lane-mine")
    assert "NOTICE" in said.stdout

    cut = without_block(tmp_path, cs.WHOLE_COPY_BEGIN, cs.WHOLE_COPY_END, "check_staged_without_whole.py")
    quiet = run_cut(shared_repo, index, cut, "--lane", "lane-mine")

    assert quiet.returncode == 0, quiet.stdout + quiet.stderr
    assert "NOTICE" not in quiet.stdout
    assert "lane-kinetics" not in quiet.stdout


# ------------------------------------------------- a removal cancelled only inside a comment
#
# The subtraction that stops a moved line being read as a revert takes content back from anywhere in
# the file, on stripped text. So a removed line pasted verbatim into a docstring stopped counting as
# removed. Quoting superseded PROSE as history is this project's practice and keeps working; pasting
# a removed line of CODE into a comment must not silence the guard.

# The quotation is re-indented, as a quotation of code inside prose is. A line pasted at its own
# original indentation is not reported by git as removed at all -- it matches the line it replaces --
# so there is nothing for any guard to subtract and nothing to refuse; that case is named in
# `live_put_back`. This is the case that was reaching the subtraction and being cancelled by it.
QUOTED_INTO_A_DOCSTRING = """def register(payload):
    \"\"\"Superseded on 2026-10-03, and this is what the line said:

        p = save_result(RESULT, payload)

    The stamp below takes the name from the path instead.
    \"\"\"
    path = RESULT
    stamp = record_stamp(p)
    return path
"""

PROSE_QUOTED_AS_HISTORY = """# Grammar

| field | type | meaning |
|---|---|---|
| `fraction` | `number` | a share of what is left when it runs |

## Superseded

Until today the table also carried, word for word:

| `share` | `number` | a share of the whole at this decision point, order-free |
"""

PROSE_BASE = '''def register(payload):
    """Register the payload."""
    path = RESULT
    return path
'''

# the peer commit adds two lines of PROSE to a python file: a docstring sentence and a comment
PROSE_PEER = '''def register(payload):
    """Register the payload.

    The name save_result wrote is not always the one that was asked for.
    """
    # keep the asked-for name, not the written one
    path = RESULT
    return path
'''

PEER_SENTENCE = "The name save_result wrote is not always the one that was asked for."
PEER_COMMENT = "# keep the asked-for name, not the written one"


def test_a_removed_code_line_quoted_into_a_docstring_is_still_refused(
    code_repo: Path, tmp_path: Path
) -> None:
    """The only removal here is the peer's line, and the only thing putting it back is a docstring."""
    index = private_index(code_repo, tmp_path)
    stage_blob(code_repo, index, "reg.py", QUOTED_INTO_A_DOCSTRING)

    done = run_check(code_repo, index)

    assert done.returncode == 2, done.stdout
    assert PEER_SUBJECT in done.stdout
    assert HELD in done.stdout


def test_prose_quoted_as_history_still_passes(repo: Path, tmp_path: Path) -> None:
    """The established additive practice: the superseded row is kept, word for word, as history."""
    index = private_index(repo, tmp_path)
    stage_blob(repo, index, "GRAMMAR.md", PROSE_QUOTED_AS_HISTORY)

    done = run_check(repo, index)

    assert done.returncode == 0, done.stdout
    assert "takes nothing back out" in done.stdout


@pytest.fixture
def prose_repo(tmp_path: Path) -> Path:
    """A python file whose recent peer commit added PROSE: a docstring sentence and a comment."""
    repo = tmp_path / "prose"
    repo.mkdir()
    git(repo, "init", "-q", "-b", "main")
    git(repo, "config", "user.email", "peer@example.test")
    git(repo, "config", "user.name", "A Peer")
    source = repo / "reg.py"
    source.write_text(PROSE_BASE)
    git(repo, "add", "reg.py")
    git(repo, "commit", "-qm", "the register skeleton")
    source.write_text(PROSE_PEER)
    git(repo, "add", "reg.py")
    git(repo, "commit", "-qm", PEER_SUBJECT)
    return repo


def test_a_peers_prose_rewrapped_in_a_python_file_still_passes(prose_repo: Path, tmp_path: Path) -> None:
    """No false refusal inside python either: a sentence and a comment re-indented word for word.

    Both lines move to a new indentation, so git reports each as a removal plus an addition -- the
    same shape as the quotation above. Neither was live code in HEAD, so neither is restricted, and
    re-wrapping a docstring stays ordinary work.
    """
    index = private_index(prose_repo, tmp_path)
    stage_blob(
        prose_repo,
        index,
        "reg.py",
        "def register(payload):\n"
        '    """Register the payload.\n'
        "\n"
        f"        {PEER_SENTENCE}\n"
        '    """\n'
        f"        {PEER_COMMENT}\n"
        "    path = RESULT\n"
        "    return path\n",
    )

    done = run_check(prose_repo, index)

    assert done.returncode == 0, done.stdout
    assert "takes nothing back out" in done.stdout


def test_a_file_that_cannot_be_tokenized_keeps_the_old_behaviour(code_repo: Path, tmp_path: Path) -> None:
    """A refusal resting on a file this script could not read would be a guess, so it does not."""
    index = private_index(code_repo, tmp_path)
    stage_blob(
        code_repo,
        index,
        "reg.py",
        "def register(payload:\n    path = RESULT\n" + f"    # {HELD}\n" + "    stamp = record_stamp(p)\n",
    )

    done = run_check(code_repo, index)

    assert done.returncode == 2, done.stdout  # the removal is still a removal; nothing cancelled it


def test_without_the_quoted_cancellation_condition_the_quotation_silences_the_guard(
    code_repo: Path, tmp_path: Path
) -> None:
    """Credited with preventing something, shown by removal: with the block cut, it goes through."""
    import check_staged as cs

    index = private_index(code_repo, tmp_path)
    stage_blob(code_repo, index, "reg.py", QUOTED_INTO_A_DOCSTRING)

    refused = run_check(code_repo, index)
    assert refused.returncode == 2, refused.stdout

    cut = without_block(tmp_path, cs.QUOTED_BEGIN, cs.QUOTED_END, "check_staged_without_quoted.py")
    allowed = run_cut(code_repo, index, cut)

    assert allowed.returncode == 0, allowed.stdout + allowed.stderr
    assert "takes nothing back out" in allowed.stdout
    assert HELD not in allowed.stdout


def test_the_tokenizer_marks_a_docstring_interior_as_not_live(tmp_path: Path) -> None:
    """The distinction itself, directly: the same text as code and as a quotation of code."""
    import check_staged as cs

    live = cs.live_code_text('def f():\n    """\n    x = 1\n    """\n    y = 2\n    # z = 3\n')

    assert live is not None
    assert "y = 2" in live
    assert "x = 1" not in live
    assert "# z = 3" not in live
    assert cs.live_code_text("def f(:\n") is None
