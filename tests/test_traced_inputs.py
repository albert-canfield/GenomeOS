# SPDX-License-Identifier: AGPL-3.0-or-later
"""Inputs are traced, not declared (2026-10-02, lane-tracer).

The defect these tests reproduce. The 24 `enhancer_targets_all_chr*` summaries declare
`compiler_inputs:enhancer_targets_all_chr*.json` -- 24 files, 39,683 bytes in all, about 1.6 KB each,
because those files are POINTERS. The bytes the run reads are the roughly 1.4 GB of
data/knowledge/alphagenome/all_elements/chr*.json, declared nowhere and hashed nowhere.
`genomeos/attribution/targets.py:run_elements` takes a path out of an `elements_where` field of a
result file and opens whatever that field names, so the path lives in DATA rather than in code and no
reading of the source can find it.

A missing input would be bad enough. What makes this worse is that a check does not stop on it:
data/knowledge is linked read-only into every checkout, so a rebuild from a clean worktree finds the
tables, the run completes, and `scripts/manifest_rebuild.py` reports success WITHOUT having hashed the
largest thing the run read. docs/LESSONS.md calls that "A verification tool can fail in the direction
that flatters". So the reconciliation is a refusal, and `test_counterfactual_*` below removes it from a
copy of `genomeos/results.py` and shows the pointer case then writes -- the standing rule here being
that a guard credited with preventing something must be shown to prevent it, by removing it.

The pointer case drives the real `run_elements`, not a reimplementation of it: a test that opened the
table itself would prove only that the hook sees an `open` call, not that the indirection this lane
exists for is caught.
"""

import json
import os
import subprocess
import tempfile
import types
from pathlib import Path

import pytest

from genomeos import manifest as mf
from genomeos import results
from genomeos.attribution import targets

REPO_ROOT = Path(__file__).resolve().parent.parent

#: Where the pointer sends the writer, as the repository names it. The real one is
#: data/knowledge/alphagenome/all_elements/<chrom>.json; data/knowledge is a read-only link.
TABLE = "data/knowledge/alphagenome/all_elements/chrT.json"


# --- the fixtures -------------------------------------------------------------------------------------


@pytest.fixture
def run(tmp_path, monkeypatch):
    """A checkout-shaped tree with a registry, a legacy allowlist and a pointer result in it.

    The test stands in it (`chdir`) and the tracer's window is reopened from there, because what the
    hook canonicalises against is data/ where the process stands. It is a git repository because the
    revision stamp every manifest carries is read from the tree the writer stands in. Returns the
    registry.
    """
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    subprocess.run(
        [
            "git",
            "-C",
            str(tmp_path),
            "-c",
            "user.email=t@t",
            "-c",
            "user.name=t",
            "commit",
            "-q",
            "--allow-empty",
            "-m",
            "the tree a writer stands in",
        ],
        check=True,
    )
    reg = tmp_path / "data" / "results"
    reg.mkdir(parents=True)
    (tmp_path / "data" / "results_legacy.txt").write_text("# a header\nold_result\ttracked\n")
    monkeypatch.setattr(results, "RESULTS_DIR", reg)
    monkeypatch.setattr(results, "LEGACY_ALLOWLIST", tmp_path / "data" / "results_legacy.txt")
    monkeypatch.chdir(tmp_path)
    table = tmp_path / TABLE
    table.parent.mkdir(parents=True)
    table.write_text(json.dumps([{"id": "e1", "predicted_coding": {"gene": "SOD1"}}] * 4))
    # the pointer: 1.6 KB of result that names the table in a FIELD, which is how the real one does it
    (reg / "enhancer_targets_all_chrT.json").write_text(
        json.dumps({"result": "enhancer_targets_all_chrT", "elements_where": TABLE})
    )
    mf.trace_begin()
    return reg


def contract(reg, declare_table=False, extra_inputs=()):
    """The manifest a writer of this shape hands `save_result`: the pointer declared as a group, the
    way the real summaries declare theirs, and the table declared only when asked for."""
    inputs = [
        mf.files_entry(
            "compiler_inputs:enhancer_targets_all_chrT.json", [reg / "enhancer_targets_all_chrT.json"]
        )
    ]
    if declare_table:
        inputs.append(mf.input_entry(TABLE, partition=None))
    inputs.extend(extra_inputs)
    return {
        "sources": [{"accession": "ENCSR000XXX", "version": "v1"}],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {"min_effect": 0.1},
        "exclusions": [],
        "partitions": "n/a: a compilation, not an evaluation",
        "code_cleanliness": mf.code_cleanliness("tests/test_traced_inputs.py", (), root=REPO_ROOT),
    }


def writer(reg, name="compiled_chrT", declare_table=False, save=results.save_result, extra_inputs=()):
    """The writer under test: it reads the elements through the real indirection and writes a result.

    `run_elements` loads the pointer result, reads its `elements_where` field and opens the path that
    field names. Nothing in this function names the table.
    """
    elements = targets.run_elements("enhancer_targets_all", "chrT", results_dir=reg)
    payload = {
        "elements_compiled": len(elements),
        "genes": sorted({e["predicted_coding"]["gene"] for e in elements}),
    }
    return save(name, payload, reg, manifest=contract(reg, declare_table, extra_inputs))


# --- the pointer case: the defect, reproduced ---------------------------------------------------------


def test_the_pointer_case_is_refused_and_the_message_names_the_table(run):
    """The point of the lane. A writer whose declared inputs are the pointer files, and which reads a
    table named by a field inside one of them, cannot enter the registry."""
    assert targets.run_elements("enhancer_targets_all", "chrT", results_dir=run), "the indirection runs"
    mf.trace_begin()

    with pytest.raises(mf.ManifestError) as e:
        writer(run)
    assert TABLE in str(e.value), "the refusal must name the file that was read and not declared"
    assert "not declared and hashed in inputs" in str(e.value)
    assert not (run / "compiled_chrT.json").exists(), "nothing may enter the registry"


def test_the_pointer_case_quarantines_the_payload_with_its_reason(run):
    """The existing quarantine path carries it, so hours of compute are kept and the run still fails."""
    with pytest.raises(mf.ManifestError, match="quarantined at"):
        writer(run)
    kept = json.loads((results.quarantine_dir(run) / "compiled_chrT.json").read_text())
    assert kept["elements_compiled"] == 4, "the payload is kept"
    assert kept[mf.KEY]["complete"] is False
    assert TABLE in kept[mf.KEY]["traced_inputs"]["undeclared"]
    assert any(TABLE in p for p in kept["quarantine"]["problems"])


def test_the_declared_pointer_group_is_1_6_kb_against_the_table_it_points_at(run):
    """Why a byte count cannot stand in for the reconciliation: the declared group is a few kilobytes
    whatever the table weighs, so no total and no count of inputs distinguishes the two cases."""
    group = mf.files_entry(
        "compiler_inputs:enhancer_targets_all_chrT.json", [run / "enhancer_targets_all_chrT.json"]
    )
    assert group["bytes"] < (Path(TABLE).stat().st_size), "the pointer is smaller than what it points at"
    assert group["files"] == 1 and group["group"] is True


# --- the counterfactual: the guard is shown to prevent it, by removing it -----------------------------

_BEGIN = "    # --- traced-input reconciliation"
_END = "    # --- end traced-input reconciliation"


def _results_without_reconciliation() -> types.ModuleType:
    """A copy of `genomeos/results.py` with the reconciliation block cut out of `save_result`.

    A real removal, from the shipped source, not a flag and not a monkeypatch: the excision is checked
    to have removed the call it names, so a block that had been renamed or moved could not pass here as
    a block that had been removed.
    """
    src = Path(results.__file__).read_text()
    start, end = src.index(_BEGIN), src.index(_END)
    assert start < end
    cut = src[:start] + src[end:]
    cut = cut[: cut.index(_END)] + cut[cut.index("\n", cut.index(_END)) + 1 :]
    assert "reconciliation_problems" not in cut, "the excision must remove the guard, not merely its comment"
    assert "traced_inputs" not in cut.split('"""')[-1], "nor leave the block attached"
    mod = types.ModuleType("results_without_reconciliation")
    mod.__file__ = results.__file__
    exec(compile(cut, "results_without_reconciliation", "exec"), mod.__dict__)  # noqa: S102
    return mod


def test_counterfactual_without_the_reconciliation_the_pointer_case_writes(run, monkeypatch):
    """Remove the guard and the harm happens: the same writer, the same undeclared 1.4 GB-shaped read,
    and the result enters the registry reporting a complete manifest. So the guard is a protection, not
    a diagnosis of something another check was already stopping."""
    mod = _results_without_reconciliation()
    mod.RESULTS_DIR = run
    mod.LEGACY_ALLOWLIST = results.LEGACY_ALLOWLIST

    p = writer(run, save=mod.save_result)

    assert p.exists(), "without the reconciliation the pointer case is admitted"
    written = json.loads(p.read_text())
    assert written[mf.KEY]["complete"] is True, "and it reports a complete manifest"
    assert "traced_inputs" not in written[mf.KEY], "with nothing recording what it read"
    assert written["elements_compiled"] == 4


def test_counterfactual_the_copy_keeps_every_other_refusal(run):
    """The excision removes this lane's block and nothing else: the copy still refuses a manifest with
    no cleanliness block, so the counterfactual above is not passing because the copy refuses nothing."""
    mod = _results_without_reconciliation()
    mod.RESULTS_DIR = run
    mod.LEGACY_ALLOWLIST = results.LEGACY_ALLOWLIST
    m = contract(run, declare_table=True)
    del m["code_cleanliness"]
    with pytest.raises(mf.ManifestError, match="code_cleanliness"):
        mod.save_result("compiled_chrT", {"n": 1}, run, manifest=m)


# --- the near-misses: what must NOT trip it -----------------------------------------------------------


def test_near_miss_declaring_the_table_writes(run):
    """The fix the refusal asks for, shown to work: the same writer, the table declared and hashed."""
    p = writer(run, declare_table=True)
    assert p.exists()
    written = json.loads(p.read_text())
    assert written[mf.KEY]["complete"] is True
    assert written[mf.KEY]["traced_inputs"]["undeclared"] == []
    assert written[mf.KEY]["traced_inputs"]["active"] is True
    declared = {i["path"] for i in written[mf.KEY]["inputs"]}
    assert TABLE in declared


def test_near_miss_a_file_written_under_data_is_output_not_input(run):
    """A writer that writes a side file under data/ is not answerable for it as an input."""
    (run.parent / "side").mkdir()
    (run.parent / "side" / "scratch.json").write_text("{}")  # a write, through the same hook
    p = writer(run, declare_table=True)
    assert p.exists()
    block = json.loads(p.read_text())[mf.KEY]["traced_inputs"]
    assert block["undeclared"] == []
    assert block["files_written"] >= 1, "the write was seen, and kept apart from the reads"


def test_near_miss_a_file_read_outside_data_is_not_an_input(run, monkeypatch):
    """Code, the lock file and anything else outside data/ is not a data input: the contract's `inputs`
    are the bytes of data a result rests on, and the revision stamp covers the code."""
    (run.parent.parent / "elsewhere.tsv").write_text("a\tb\n")
    Path(run.parent.parent / "elsewhere.tsv").read_text()
    p = writer(run, declare_table=True)
    assert json.loads(p.read_text())[mf.KEY]["traced_inputs"]["undeclared"] == []


def test_near_miss_a_file_under_a_declared_directory_is_declared(run):
    """A declared directory is hashed as a tree (`sha256_of`), so a file under it is covered by prefix
    and does not have to be named on its own."""
    d = Path("data/knowledge/alphagenome/all_elements")
    (d / "extra.json").write_text("[]")
    (d / "extra.json").read_text()
    p = writer(run, declare_table=True, extra_inputs=[mf.input_entry(str(d), partition=None)])
    assert json.loads(p.read_text())[mf.KEY]["traced_inputs"]["undeclared"] == []


def test_near_miss_the_window_does_not_carry_an_earlier_results_reads(run):
    """Each result is answerable for its own reads: a chain writer's second result is not refused for a
    table its first result read and declared."""
    writer(run, name="first", declare_table=True)
    p = results.save_result("second", {"n": 1}, run, manifest=contract(run, declare_table=False))
    assert p.exists(), "the second result reads only the pointer group it declares"
    block = json.loads(p.read_text())[mf.KEY]["traced_inputs"]
    assert block["undeclared"] == [], (
        "nor is the second result charged with data/results_legacy.txt, which the enforcement "
        "machinery itself opened while writing the first"
    )


def test_near_miss_two_spellings_of_one_file_are_one_input(run, tmp_path):
    """data/reference, data/knowledge and data/cache are links, so one file has two paths. Reading it
    through the link and through its target must be one input, and declaring it by either spelling must
    satisfy the reconciliation."""
    store = tmp_path / "store"
    (store / "tables").mkdir(parents=True)
    (store / "tables" / "t.json").write_text("[]")
    link = Path("data") / "linked"
    link.symlink_to(store)
    mf.trace_begin()

    (link / "tables" / "t.json").read_text()  # through the link
    (store / "tables" / "t.json").read_text()  # through its target
    held = mf.trace_close()
    assert held["reads"] == ["data/linked/tables/t.json"], "one file, one spelling"

    declared = {"inputs": [{"path": str(store / "tables" / "t.json"), "sha256": "x", "partition": None}]}
    assert mf.undeclared_reads(declared, held["reads"]) == [], "declared by its target's spelling"
    declared = {"inputs": [{"path": "data/linked/tables/t.json", "sha256": "x", "partition": None}]}
    assert mf.undeclared_reads(declared, held["reads"]) == [], "declared through the link"


# --- the legacy allowlist keeps working ---------------------------------------------------------------


def test_a_legacy_allowlist_name_still_writes_with_an_undeclared_read(run):
    """The names written before the contract must keep working. A legacy name warns and writes, exactly
    as it does for the rest of the contract, so the historical writers keep regenerating."""
    with pytest.warns(results.ManifestWarning, match="not declared and hashed"):
        p = writer(run, name="old_result")
    assert p.exists(), "a legacy name is not refused"
    written = json.loads(p.read_text())
    assert written["elements_compiled"] == 4
    assert TABLE in written[mf.KEY]["traced_inputs"]["undeclared"], "but the finding is on the record"
    assert not (results.quarantine_dir(run) / "old_result.json").exists()


def test_the_rule_does_not_reach_outside_the_registry(run, tmp_path):
    """A scratch directory is not the registry: the rule is about what enters the registry. The reading
    is still on the record, so nothing is hidden by not being enforced."""
    out = tmp_path / "scratch"
    out.mkdir()
    targets.run_elements("enhancer_targets_all", "chrT", results_dir=run)  # the undeclared read
    p = results.save_result("anything", {"n": 1}, out, manifest=contract(run))
    assert p.exists()
    assert TABLE in json.loads(p.read_text())[mf.KEY]["traced_inputs"]["undeclared"]


def test_strict_enforces_outside_the_registry(run, tmp_path):
    """`strict=True` enforces in any directory, which is what the docstring has always promised."""
    out = tmp_path / "scratch"
    out.mkdir()
    targets.run_elements("enhancer_targets_all", "chrT", results_dir=run)
    with pytest.raises(mf.ManifestError, match="not declared and hashed"):
        results.save_result("anything", {"n": 1}, out, manifest=contract(run), strict=True)


# --- the other direction: an over-declaration is reported, not refused --------------------------------
#
# The second manifest fault, opposite in direction (lane-rename, 2026-10-02): the committed
# `context_evidence` names 30 domains_chr*.json, 25 unknown_chr*.json and 25 ccres_chr*.bed.gz where the
# run reads 24 of each, because a later commit narrowed the globs after the result was written. The
# choice made here is deliberate, not an accident of how the set difference was written.


def test_an_over_declaration_is_reported_and_not_refused(run):
    """Declaring more than you read is a wrong record, not a verification that passes by not looking:
    the files it names are pinned, they were simply not needed. Refusing it would fire on results that
    are honest, so it is reported and the result writes."""
    never = Path("data/knowledge/alphagenome/all_elements/chrNEVER.json")
    never.write_text("[]")
    entry = mf.input_entry(str(never), partition=None)  # hashed here, before the window opens
    mf.trace_begin()

    p = writer(run, declare_table=True, extra_inputs=[entry])

    assert p.exists(), "an over-declaration does not keep a result out of the registry"
    block = json.loads(p.read_text())[mf.KEY]["traced_inputs"]
    assert block["undeclared"] == []
    assert str(never) in block["declared_not_read"], "but it is on the record"


def test_a_declared_directory_with_a_read_under_it_is_not_an_over_declaration(run):
    """A directory is declared by name and read by its members, so it must not report itself unread."""
    d = Path("data/knowledge/alphagenome/all_elements")
    p = writer(run, declare_table=True, extra_inputs=[mf.input_entry(str(d), partition=None)])
    block = json.loads(p.read_text())[mf.KEY]["traced_inputs"]
    assert str(d) not in block["declared_not_read"]
    assert block["undeclared"] == []


def test_the_limit_a_file_opened_only_to_hash_it_counts_as_read():
    """What this cannot establish, stated as a test rather than left to be assumed. A manifest is built
    by hashing what it declares, and hashing opens the file, so an input declared but never used is
    recorded as read. The context_evidence over-declaration is of exactly that shape and the hook cannot
    see it; only reading the writer's globs against its reads can."""
    with tempfile.TemporaryDirectory() as d:
        root = Path(d)
        (root / "data").mkdir()
        (root / "data" / "unused.json").write_text("[]")
        here = os.getcwd()
        try:
            os.chdir(root)
            mf.trace_begin()
            entry = mf.input_entry("data/unused.json", partition=None)  # hashed, and so opened
            reads = mf.trace_close()["reads"]
        finally:
            os.chdir(here)
    assert "data/unused.json" in reads
    assert mf.declared_not_read({"inputs": [entry]}, reads) == [], "hashing it counts as reading it"


# --- an unwatched result does not read like a clean one -----------------------------------------------


def test_an_inactive_tracer_is_itself_a_problem(run, monkeypatch):
    """A result written with nothing watching cannot say its declared inputs are what it read, so the
    absence of a record is a refusal rather than a pass. A check that passes because it did not look is
    the failure docs/LESSONS.md records."""
    unwatched = dict(mf.trace_close(), active=False)
    monkeypatch.setattr(mf, "trace_close", lambda: unwatched)
    with pytest.raises(mf.ManifestError, match="the input tracer was not active"):
        writer(run, declare_table=True)


def test_reconciliation_problems_on_a_block_that_is_not_a_block():
    for bad in (None, "a string", 42, [], {}):
        assert mf.reconciliation_problems(bad) != []


# --- the pieces --------------------------------------------------------------------------------------


def test_declared_inputs_reads_a_group_member_by_member_and_not_by_its_label():
    """A group's label is not a file; its members are, each with its own sha256. That is what makes a
    large store hashable member by member rather than as a directory."""
    m = {
        "inputs": [
            {
                "path": "compiler_inputs:x",
                "group": True,
                "sha256": "g",
                "members": [
                    {"path": "data/a.json", "sha256": "aa", "bytes": 1},
                    {"path": "data/b.json", "sha256": "bb", "bytes": 2},
                ],
            }
        ]
    }
    files, dirs = mf.declared_inputs(m)
    assert files == {"data/a.json", "data/b.json"}
    assert dirs == set()
    assert "compiler_inputs:x" not in files
    assert mf.undeclared_reads(m, ["data/a.json", "data/c.json"]) == ["data/c.json"]


def test_a_group_member_with_no_digest_is_not_declared():
    """Declared is not enough: the member must have been hashed."""
    m = {"inputs": [{"path": "g", "group": True, "members": [{"path": "data/a.json"}]}]}
    assert mf.undeclared_reads(m, ["data/a.json"]) == ["data/a.json"]


def test_group_digest_hashes_each_member_once_in_one_pass(tmp_path, monkeypatch):
    """Requirement of the lane: a large store is hashed member by member, once, at write time. The
    group digest and every member digest come out of one read of each file."""
    monkeypatch.chdir(tmp_path)
    for n, text in (("a.json", "aaa"), ("b.json", "bbbb")):
        Path(n).write_text(text)
    mf.trace_begin()
    digest, total, members = mf.group_digest(["a.json", "b.json"])
    assert total == 7
    assert [m["bytes"] for m in members] == [3, 4]
    assert len({m["sha256"] for m in members}) == 2, "each member has its own digest, not the group's"
    assert digest != members[0]["sha256"]
    opens = mf.trace_close()["opens"]
    assert opens == 2, f"one read per member, not two passes: {opens}"


@pytest.mark.parametrize(
    ("mode", "flags", "write"),
    [("r", 0, False), ("rb", 0, False), ("w", 0, True), ("a", 0, True), ("r+", 0, True), ("x", 0, True)],
)
def test_the_hook_tells_a_read_from_a_write(mode, flags, write):
    assert mf._is_write(mode, flags) is write


def test_the_hook_tells_a_read_from_a_write_for_os_open():
    assert mf._is_write(None, os.O_RDONLY | os.O_CLOEXEC) is False
    assert mf._is_write(None, os.O_WRONLY | os.O_CREAT) is True


def test_a_file_descriptor_is_not_a_path(tmp_path, monkeypatch):
    """`open(fd)` raises its own event with an int where a path would be; the open that produced the
    descriptor was already recorded, so the int contributes nothing rather than an exception."""
    assert mf._TRACE.canonical(7) is None
    assert mf._TRACE.canonical(None) is None
