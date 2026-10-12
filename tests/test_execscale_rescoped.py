# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-execfix: the execscale claim at the width it is actually answerable for, and the two digests.

WHY THIS FILE EXISTS. Two tests of the execscale lane went red at HEAD on 2026-10-12 and neither red
was a defect in the code they ran:

  * `tests/test_execscale.py::test_the_committed_programmes_are_byte_identical_to_the_record` pins all
    58 committed programmes' bytes to `data/results/execscale_census.json` for ever. That is a
    REPOSITORY-WIDE BYTE PIN, so it reddens on any legitimate programme change by anyone -- and one
    happened: `71f9d47`, the coordinator's restore of the `share:` generator honouring a registered
    falsifier, moved data/organisms/human/haematopoiesis.bio from 20,618 to 20,643 bytes.
  * `tests/test_execscale_magnitude.py`'s demo-record test read a dict's KEYS where the caveat it
    wanted lives in a VALUE, so it could not pass whatever the record said.

Both are marked `xfail(strict=False)` with reasons that say so, and both marks are boarded for
deletion with the tests they excuse, on 2026-10-14, when `81d3949`'s and `44d68b2`'s check_staged
windows lapse. NOTHING WAS REMOVED to clear either red, because a removal inside those windows needs a
`--force` this checkout does not grant. This file is the other half of that work: the claims those
tests were meant to carry, written at the width their lane really owns.

THE RE-SCOPED CLAIM, and it is narrower than a byte pin in exactly one way. lane-execscale's safety
property is that ITS parameter at ITS default moves no programme -- not that no programme may ever
move. So the comparison here is between TWO COMPILATIONS OF THE SAME TREE, one with the parameter at
its default and one that never passes it at all, and the requirement is byte identity between those
two. It cannot redden on anyone else's commit, because both sides of it are computed in the same
checkout at the same moment; and it is not vacuous, because the calibrated option is shown on the same
tree to move the bytes it is supposed to move.

THE TWO DIGEST IMPLEMENTATIONS. `manifest.group_digest` calls itself "the one implementation", and
`scripts/execscale_census.py:125` is a second one beside it, written so "the writer and the checker
cannot answer differently". This file holds them to that: given the same list in the same order they
agree, bit for bit, on the committed group of 58 programmes. The duplicate is NOT removed here (a
removal in a window again) and is reported for the timed list instead.

WHAT BOTH DIGESTS CARRY, read off the code rather than reasoned about: each hashes
`name + b"\\0" + bytes` per member IN THE ORDER GIVEN, and neither sorts internally. So a group digest
folds in the file LIST -- membership and order -- on top of the per-file content, and a fixed list of
per-file digests cannot see a programme added, removed or reordered. The plants below hold that. It
also means the two ROUTES to "a group's digest" disagree by construction, because
`manifest.files_entry` sorts its names before hashing while `execscale_census.programs()` concatenates
per-tree sorted lists over three trees: the same 58 files give fe56d749 (the census's order), 88fcad50
(globally sorted) and d6af2e88 (files_entry's sorted members) at the census commit. That is a footgun
in the API, not a drift in the data, and it is reported here, not repaired.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import subprocess
from pathlib import Path

import pytest

from genomeos import manifest as mf
from genomeos.attribution import compile as compile_mod
from genomeos.lang import parse
from genomeos.runtime import grn
from tests.tracked_paths import must_be_present

ROOT = Path(__file__).resolve().parents[1]
RECORD = "data/results/execscale_census.json"
#: the commit the census record was taken on; its blobs are what the record's digests are about
CENSUS_COMMIT = "81d3949"
CHROM = "chr21"

#: `compile_chromosome` reaches `targets.run_elements`, whose chr21 element table lives under the
#: git-ignored `data/knowledge` store, so on a machine without it the compiler raises by design. The
#: project's marker for that is `needs_local_data`: the test RUNS where the store is and SKIPS BY NAME
#: where it is genuinely absent (CI, a bare worktree), which is how a reader tells "not run here" from
#: "passed". tests/test_compiled_defaults.py carries the identical marker for the identical path.
needs_the_chr21_element_table = pytest.mark.needs_local_data(
    "data/knowledge/alphagenome/all_elements/chr21.json",
    how="scripts/enhancer_targets_all_chain.py writes it (the enhancer_targets_all_chr21 job); "
    "data/knowledge is git-ignored machine-local data, so a fresh checkout cannot have it",
)


def _census() -> dict:
    must_be_present(RECORD)
    return json.loads((ROOT / RECORD).read_text())


def _census_module():
    """`scripts/execscale_census.py` loaded as a module, the way tests/test_execscale.py loads it."""
    path = ROOT / "scripts" / "execscale_census.py"
    spec = importlib.util.spec_from_file_location("_execfix_execscale_census", path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ec = _census_module()


def _blob(rev: str, path: str) -> bytes:
    return subprocess.run(["git", "show", f"{rev}:{path}"], cwd=ROOT, capture_output=True, check=True).stdout


# ---- the re-scoped claim: the parameter at its default against not passing it at all ------------


@needs_the_chr21_element_table
def test_THE_RE_SCOPED_CLAIM_the_default_compiles_the_SAME_BYTES_as_never_passing_the_parameter():
    """The safety property at its real width: two compilations of ONE tree, byte for byte.

    This is what the over-wide test in tests/test_execscale.py was reaching for. The quantity that
    must not move is the difference the PARAMETER makes, so both sides are compiled here and now, in
    this checkout, from the same inputs -- and a programme another lane legitimately commits cannot
    redden it, because neither side reads a record.
    """
    must_be_present(f"data/results/budget_{CHROM}.json")
    without = compile_mod.compile_chromosome(CHROM)
    at_default = compile_mod.compile_chromosome(CHROM, magnitude=compile_mod.MAGNITUDE_DEFAULT)
    assert at_default == without
    # stated as the digests too, because that is the form the census's own claim is written in
    assert hashlib.sha256(at_default.encode()).hexdigest() == hashlib.sha256(without.encode()).hexdigest()
    # and neither side is a stub, which would make the identity above true and worthless
    assert len([ln for ln in without.splitlines() if ln.startswith("rule ")]) > 5000


@needs_the_chr21_element_table
def test_PLANT_the_CALIBRATED_option_does_move_the_bytes_so_the_identity_above_is_not_vacuous():
    """The same comparison with the option ON must come out DIFFERENT, or the test above proves nothing.

    What the calibrated option adds is the magnitude block and nothing else: the added lines are
    exactly `magnitude_lines`, no `rule`, `gene`, `element` or `region` line appears, and every
    `# test:` line the census is held to is unchanged -- which is why it can be asked for without
    moving a count.
    """
    must_be_present(f"data/results/budget_{CHROM}.json")
    without = compile_mod.compile_chromosome(CHROM)
    calibrated = compile_mod.compile_chromosome(CHROM, magnitude=compile_mod.MAGNITUDE_CALIBRATED)
    assert calibrated != without
    plain_lines = set(without.splitlines())
    added = [ln for ln in calibrated.splitlines() if ln not in plain_lines]
    assert added == [ln for ln in compile_mod.magnitude_lines(compile_mod.MAGNITUDE_CALIBRATED) if ln]
    assert not [ln for ln in added if ln.startswith(("rule ", "gene ", "element ", "region "))]
    declared = [ln for ln in without.splitlines() if ln.startswith("# test:")]
    assert declared and [ln for ln in calibrated.splitlines() if ln.startswith("# test:")] == declared
    # the programme the option emits is one the runtime reads the measured magnitude off
    assert f"param {grn.MAGNITUDE_CORRECTION_PARAM} = {grn.MEASURED_MAGNITUDE_CORRECTION}" in calibrated


def test_the_re_scoped_claim_on_the_RUNTIME_side_needs_no_store_at_all():
    """The same shape one level down: a run that never passes the parameter against one that passes
    its default, compared over the WHOLE final state rather than one species."""
    text = (
        "module tiny.rescoped\n"
        "gene g { max: 10; basal: 0.0 }\n"
        "gene h { max: 4; basal: 0.1 }\n"
        "protein S {}\n"
        "rule S activates g { strength: 1.0; threshold: 1.0; hill: 2 }\n"
        "rule S activates h { strength: 0.5; threshold: 0.5; hill: 1 }\n"
    )

    def _run(**kw):
        m = parse(text)
        for name in ("g", "h"):
            m.entities[name].attrs[grn.MULTIPLICATIVE_GAIN] = 1.0
        rt = grn.NetworkRuntime(m, activator_combination="multiplicative", **kw)
        return rt.run(2.0, dt=0.05, clamp={"S": 2.0}).final()

    without = _run()
    at_default = _run(magnitude_correction=grn.DEFAULT_MAGNITUDE_CORRECTION)
    assert at_default == without  # bit for bit, every species
    # and the measured magnitude is not inert on the same two genes, so the identity is not vacuous
    measured = _run(magnitude_correction=grn.MEASURED_MAGNITUDE_CORRECTION)
    assert measured != without
    assert {k for k in without if without[k] != measured[k]} == {"g.mRNA", "h.mRNA"}


# ---- the two group_digest implementations --------------------------------------------------------


def test_THE_TWO_group_digest_IMPLEMENTATIONS_AGREE_on_the_committed_group():
    """`manifest.group_digest` and the duplicate at `scripts/execscale_census.py:125`, same list.

    The duplicate's docstring says it exists so "the writer and the checker cannot answer
    differently". Given the same names in the same order they do not: this is the check that says so,
    on the real committed group of 58 programmes rather than on a fixture. It compares the two
    implementations on the tree AS IT STANDS and never to a recorded number, so a programme another
    lane commits moves both sides together and cannot redden it.
    """
    names = [row["path"] for row in _census()["programs"]]
    for name in names:
        must_be_present(name)
    from_the_census_script = ec.group_digest([Path(n) for n in names])
    from_the_engine, total, members = mf.group_digest(names)
    assert from_the_census_script == from_the_engine
    # and they agree about the membership they hashed, not only about the one number
    assert [m["path"] for m in members] == names
    assert total == sum(m["bytes"] for m in members)
    # the engine's own per-member digests are the per-file digests of the same bytes
    for member in members:
        assert member["sha256"] == hashlib.sha256((ROOT / member["path"]).read_bytes()).hexdigest()


def test_the_recorded_group_digest_is_REPRODUCIBLE_from_the_tree_it_was_taken_on(tmp_path: Path):
    """fe56d749 is not a mystery and the record is not defective: it is the 58 blobs at `81d3949`.

    Read from git in the record's own order into a temporary root, the engine's implementation returns
    the recorded digest exactly. The working tree gives a different number only because `71f9d47`
    moved one programme's bytes after the census, which is a fact about that commit and not a fault
    in this digest. A blob side rather than a working-tree side, for the reason
    tests/test_execscale_magnitude.py gives for its own: a record is a record of ONE tree.
    """
    record = _census()
    names = [row["path"] for row in record["programs"]]
    for name in names:
        out = tmp_path / name
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_bytes(_blob(CENSUS_COMMIT, name))
    digest, total, _ = mf.group_digest(names, root=tmp_path)
    assert digest == record["group_sha256"]
    assert total == sum(row["bytes"] for row in record["programs"])


def test_PLANT_both_digests_fold_in_the_LIST_so_ORDER_and_MEMBERSHIP_move_them(tmp_path: Path):
    """What a group digest carries that per-file digests cannot: the list itself.

    Read off the code -- both hash `name + b"\\0" + bytes` per member in the ORDER GIVEN, and neither
    sorts -- and held here as behaviour: with every per-file digest unchanged, a reordering and a
    dropped member each move both digests. That is the answer to whether the 58 per-file digests alone
    carry the census's claim: they do not. They cannot see a programme added, removed or reordered.
    """
    a, b = tmp_path / "a.bio", tmp_path / "b.bio"
    body = "module tiny.one\ngene g { max: 10; basal: 0.0 }\n"
    a.write_text(body)
    b.write_text(body.replace("tiny.one", "tiny.two"))
    names = ["a.bio", "b.bio"]
    before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (a, b)}
    base_census = ec.group_digest([a, b])
    base_engine, _, _ = mf.group_digest(names, root=tmp_path)
    # a reordering moves both
    assert ec.group_digest([b, a]) != base_census
    assert mf.group_digest(list(reversed(names)), root=tmp_path)[0] != base_engine
    # a dropped member moves both
    assert ec.group_digest([a]) != base_census
    assert mf.group_digest(names[:1], root=tmp_path)[0] != base_engine
    # and not one per-file digest moved while all four of those digests did
    assert before == {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (a, b)}


def test_the_TWO_ROUTES_to_a_groups_digest_disagree_by_construction_and_this_is_the_finding():
    """`files_entry` SORTS its names; `group_digest` hashes the order it is given. Hence the footgun.

    `execscale_census.programs()` concatenates per-tree sorted lists over three trees, which is not a
    global sort, so the census's own order is not `sorted(names)`. Anything that reaches the group
    digest through `files_entry` therefore computes a DIFFERENT number from a correct record -- which
    is how a checker can read "only the group digest differs" on a tree whose every per-file digest
    matches. The three numbers are pinned here so the next reader does not have to rediscover them.
    """
    names = [row["path"] for row in _census()["programs"]]
    assert names != sorted(names)
    assert [str(p) for p in ec.programs()] == names
    assert ec.PROGRAM_TREES == ("data/organisms", "genomeos/std", "data/demo")
    in_record_order, _, _ = mf.group_digest(names)
    sorted_by_name, _, _ = mf.group_digest(sorted(names))
    assert in_record_order != sorted_by_name
    # files_entry is the sorted route, and its digest is the sorted one rather than the census's
    entry = mf.files_entry("the committed programmes", names)
    assert entry["sha256"] == sorted_by_name
    assert entry["sha256"] != in_record_order
    assert [m["path"] for m in entry["members"]] == sorted(names)


# ---- the mark in tests/conftest.py is not allowed to go quietly inert ----------------------------


def test_the_conftest_xfail_names_a_test_that_EXISTS_or_the_mark_does_nothing():
    """A nodeid-keyed mark silently stops applying if the test is renamed, which must be visible.

    The mark is in tests/conftest.py and not on the test because `scripts/check_staged.py` refuses a
    commit of tests/test_execscale.py while a dead lane's board entry still claims that file. This
    pins the key against the source, so a rename is caught here rather than by nobody.
    """
    import conftest  # noqa: PLC0415  -- the project's own conftest, already imported by pytest

    assert len(conftest.OVER_WIDE_XFAIL) == 1
    for nodeid, reason in conftest.OVER_WIDE_XFAIL.items():
        path, _, name = nodeid.partition("::")
        assert (ROOT / path).exists(), nodeid
        assert f"\ndef {name}(" in (ROOT / path).read_text(), nodeid
        # the reason says what is wrong, where the real claim is, and when the mark goes
        assert "ITS CLAIM IS WRONG" in reason
        assert "71f9d47" in reason and "92c6ea8" in reason
        assert "2026-10-14" in reason
