# SPDX-License-Identifier: AGPL-3.0-or-later
"""The two things a re-record under a new name rests on.

1. `--result` changes the name written and **nothing else**. Four results had to be re-recorded under
   new names on 2026-10-02 (their manifests named no member of any input group, and their committed
   bytes are a declared input of ten other committed results, so they could not be rewritten in place),
   and the whole claim "no figure moved" is worthless if the writer's numbers can depend on the name it
   writes under. That is checked from the source of every writer involved rather than from one run: the
   value parsed from `--result` may reach `save_result`, a `print`, or a comparison that refuses the run
   -- never anything a figure is computed from.
2. `scripts/rerecord_compare.py` sorts the differences between the two files into the ones the rename
   accounts for and the figures that moved, and it has to fail in the right direction: a figure that
   moved must never land in the forgiven pile. A verification tool that errs towards flattering is this
   project's own documented failure mode (docs/LESSONS.md, 2026-10-01), so the classifier is tested on
   a difference it must not forgive as well as on the ones it must.
"""

from __future__ import annotations

import ast
import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import rerecord_compare as rc  # noqa: E402

from genomeos import manifest as mf  # noqa: E402

#: The writers of the four results re-recorded under new names that can be asked for a name today.
WRITERS = (
    "scripts/context_evidence_census.py",
    "scripts/context_evidence_baserate.py",
    "scripts/placement_audit.py",
    "scripts/n1_register.py",
)

#: The fourth, which cannot be asked yet, and why. `p = save_result(RESULT, payload)` has to become
#: `p = save_result(result_name, payload)`, and that line was added by b359867 inside the two days
#: scripts/check_staged.py protects, so the change cannot be committed without --force, which is not
#: available to the lane that wrote this. The coordinator has ruled the removal correct. When it lands,
#: move this path into WRITERS above and delete this constant: a list of three with no record of the
#: fourth would read as though all four were done. It is not asserted on, because these checks read the
#: working tree and the change is already sitting in it, uncommittable.


def _parents(tree: ast.AST) -> dict[ast.AST, ast.AST]:
    out: dict[ast.AST, ast.AST] = {}
    for node in ast.walk(tree):
        for child in ast.iter_child_nodes(node):
            out[child] = node
    return out


def _call_name(node: ast.AST) -> str | None:
    if not isinstance(node, ast.Call):
        return None
    f = node.func
    return f.id if isinstance(f, ast.Name) else (f.attr if isinstance(f, ast.Attribute) else None)


def _holders(tree: ast.AST) -> set[str]:
    """The names that hold the value of `--result`: `args.result` and anything assigned from it.

    Spelling, not scope: a writer that binds the result name to a name it also uses for something else
    (`for name in ...`) cannot be followed this way, so `_reused` refuses such a name rather than
    reading on and reporting nothing.
    """
    holders = {"args.result"}
    for _ in range(4):  # to a fixpoint; four passes is more than any writer needs
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign) and _spelling(node.value) in holders:
                holders |= {_spelling(t) for t in node.targets if _spelling(t)}
    return {h for h in holders if h}


def _reused(tree: ast.AST, holders: set[str]) -> list[str]:
    """Every place a name that holds the result name is bound to something else, which would make the
    rest of this check read a different variable and say nothing. A writer must not reuse the name."""
    out = []
    for node in ast.walk(tree):
        bound: list[ast.AST] = []
        if isinstance(node, ast.Assign) and _spelling(node.value) not in holders:
            bound = list(node.targets)
        elif isinstance(node, (ast.For, ast.comprehension)):
            bound = [node.target]
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            bound = [ast.Name(id=a.arg) for a in node.args.args + node.args.kwonlyargs]
        for t in bound:
            for sub in ast.walk(t):
                if (s := _spelling(sub)) in holders and s != "args.result":
                    out.append(f"line {getattr(node, 'lineno', 0)}: {s} is bound to something else")
    return out


def _spelling(node: ast.AST | None) -> str | None:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        base = _spelling(node.value)
        return f"{base}.{node.attr}" if base else None
    return None


def _how_it_is_used(node: ast.AST, parents: dict[ast.AST, ast.AST]) -> str:
    """Why this use of the result name is harmless, or `"a value the run computes"` when it is not."""
    cur: ast.AST | None = node
    while cur is not None and not isinstance(cur, ast.stmt):
        up = parents.get(cur)
        if isinstance(up, ast.Compare):
            return "a comparison that can refuse the run"
        if (name := _call_name(up)) in {"save_result", "print"}:
            return f"an argument of {name}()"
        cur = up
    if isinstance(cur, ast.Assign) and _spelling(cur.value) is not None:
        return "assigned to another name, followed below"
    return "a value the run computes"


@pytest.mark.parametrize("writer", WRITERS)
def test_the_result_option_changes_only_the_name_that_is_written(writer):
    """Every use of `--result`'s value in the writer is the output name, a message, or a refusal."""
    tree = ast.parse((ROOT / writer).read_text())
    parents, holders = _parents(tree), _holders(tree)
    assert "args.result" in holders, f"{writer} does not parse a --result option"
    assert _reused(tree, holders) == [], f"{writer} reuses the name that holds --result"
    problems = []
    for node in ast.walk(tree):
        if not isinstance(node, (ast.Name, ast.Attribute)) or _spelling(node) not in holders:
            continue
        if isinstance(node.ctx, (ast.Store, ast.Del)):
            continue
        how = _how_it_is_used(node, parents)
        if how == "a value the run computes":
            problems.append(f"{writer}:{node.lineno} {_spelling(node)} is {how}")
    assert problems == []


def test_every_writer_of_the_four_offers_the_option_at_all():
    """A writer that cannot be asked for another name can only re-record in place, which is what the
    freeze on these four results forbids."""
    missing = [w for w in WRITERS if '"--result"' not in (ROOT / w).read_text()]
    assert missing == []


def _two_files(tmp_path: Path, figure: int, members: bool) -> tuple[Path, Path]:
    old = {
        "result": "r",
        "date": "2026-10-01",
        "rules": 440589,
        "open_share_of_assessable_rules": 0.3252,
        "per_state": {"open_in_reader": figure},
        "result_manifest": {
            "inputs": [{"path": "group_label", "sha256": "ab", "bytes": 1, "files": 3, "partition": None}],
            "code": {"git_sha": "old", "dirty": False},
            "code_cleanliness": {
                "git_sha": "old",
                "own_code_is_committed": True,
                "foreign_uncommitted_code_on_the_counting_path": [],
            },
        },
    }
    new = json.loads(json.dumps(old))
    new["result"] = "r_v2"
    new["date"] = "2026-10-02"
    new["per_state"]["open_in_reader"] = 26551
    new["result_manifest"]["code"] = {"git_sha": "new", "dirty": False}
    new["result_manifest"]["code_cleanliness"]["git_sha"] = "new"
    new["result_manifest"]["supersedes"] = {"file": "data/results/r.json", "sha256": "ab"}
    if members:
        new["result_manifest"]["inputs"][0]["group"] = True
        new["result_manifest"]["inputs"][0]["members"] = [{"path": "a", "sha256": "c", "bytes": 1}]
    a, b = tmp_path / "old.json", tmp_path / "new.json"
    a.write_text(json.dumps(old))
    b.write_text(json.dumps(new))
    return b, a


def test_the_rename_accounts_for_the_name_the_supersedes_note_and_the_members(tmp_path):
    new, old = _two_files(tmp_path, figure=26551, members=True)
    r = rc.compare(new, old)
    assert r["no_figure_moved"] is True
    assert r["figures_that_moved"] == []
    assert sorted(d.split(":")[0] for d in r["accounted_for_by_the_rename"]) == [
        "/result",
        "/result_manifest/code_cleanliness/git_sha",
        "/result_manifest/inputs[0]/group",
        "/result_manifest/inputs[0]/members",
        "/result_manifest/supersedes",
    ]
    assert r["committed"]["inputs"]["groups_naming_their_members"] == 0
    assert r["committed"]["inputs"]["declared_paths"] == 1  # the label, resolving to nothing
    assert r["re_recorded"]["inputs"]["groups_naming_their_members"] == 1
    assert r["re_recorded"]["inputs"]["member_files_with_their_own_sha256"] == 1
    assert r["must_hold_failures"] == []


def test_a_figure_that_moved_is_never_forgiven(tmp_path):
    """The direction this must fail in. One leaf differs and nothing else."""
    new, old = _two_files(tmp_path, figure=26550, members=True)
    r = rc.compare(new, old)
    assert r["figures_that_moved"] == ["/per_state/open_in_reader: 26550 vs 26551"]
    assert r["no_figure_moved"] is False
    assert rc.main([str(new), str(old)]) == 1


def test_uncommitted_own_code_on_either_side_is_reported_although_cleanliness_is_exempt(tmp_path):
    new, old = _two_files(tmp_path, figure=26551, members=True)
    payload = json.loads(new.read_text())
    payload["result_manifest"]["code_cleanliness"]["own_code_is_committed"] = False
    new.write_text(json.dumps(payload))
    r = rc.compare(new, old)
    assert r["no_figure_moved"] is True  # the block itself is exempt
    assert r["must_hold_failures"] == [
        "rebuilt/result_manifest/code_cleanliness/own_code_is_committed: False, must be True"
    ]


def test_the_classifier_forgives_no_other_manifest_key(tmp_path):
    new, old = _two_files(tmp_path, figure=26551, members=False)
    payload = json.loads(new.read_text())
    payload["result_manifest"]["parameters"] = {"openness_call": "changed"}
    payload["result_manifest"]["inputs"][0]["sha256"] = "ff"
    new.write_text(json.dumps(payload))
    r = rc.compare(new, old)
    assert sorted(d.split(":")[0] for d in r["figures_that_moved"]) == [
        "/result_manifest/inputs[0]/sha256",
        "/result_manifest/parameters",
    ]


def _writer_module(writer: str):
    spec = importlib.util.spec_from_file_location(Path(writer).stem, ROOT / writer)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("writer", WRITERS)
def test_the_supersedes_note_carries_the_committed_bytes_it_is_written_beside(writer):
    """`--supersedes` is the other half: a re-record under a new name is only readable if the new file
    says which committed file it is written beside, at which sha256, and why that one is kept."""
    note = _writer_module(writer).supersedes()
    committed = ROOT / note["file"]
    assert committed.is_file(), note["file"]
    assert note["sha256"] == mf.input_entry(committed)["sha256"]
    assert note["bytes"] == committed.stat().st_size
    assert note["kept"].startswith("unchanged")
    assert "new name" in note["why"]
    assert note["declared_as_an_input_by"], "the reason the committed bytes cannot be replaced"
    for other in note["declared_as_an_input_by"]:
        declared = json.loads((ROOT / other).read_text())["result_manifest"]["inputs"]
        assert note["file"] in [i.get("path") for i in declared], f"{other} does not declare {note['file']}"
