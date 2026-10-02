# SPDX-License-Identifier: AGPL-3.0-or-later
"""The entry-script census: its verdicts, its boundary split, and the wording it carries.

The census is READ-ONLY with respect to every result it examines, and one test holds it to that by
hashing a result's bytes either side of a read. The rest hold the thing the lane was opened for: the
918 results that record no `code.argv` are not one finding, and a test that let the split fold an
unplaceable result into whichever side was smaller would let the artefact read better than the record.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
HEADLINES = ROOT / "data/results/manifest_headlines.json"


def _load():
    spec = importlib.util.spec_from_file_location(
        "entry_script_census", ROOT / "scripts/entry_script_census.py"
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


census_mod = _load()


def _git_history_is_there() -> bool:
    """Whether this checkout can answer about its own history: a shallow clone cannot."""
    out = census_mod._git("log", "-1", "--format=%H", "--", "genomeos/manifest.py")
    return bool(out and out.strip())


# --- the verdict for one recorded argv[0] -------------------------------------------------------


def test_verdict_for_reads_each_form_from_the_record():
    assert census_mod.verdict_for(None) == "unrecorded"
    assert census_mod.verdict_for("-c") == "inline_source"
    assert census_mod.verdict_for("-") == "stdin"
    assert census_mod.verdict_for("scripts/entry_script_census.py") == "inside_repository"
    assert census_mod.verdict_for("/private/tmp/scratchpad/write.py") == "outside_repository"
    assert census_mod.verdict_for(str(census_mod.ROOT / "scripts/x.py")) == "inside_recorded_absolute"


def test_every_verdict_the_census_can_return_is_described():
    """A verdict with no entry in VERDICTS would be a column a reader has to guess at."""
    returned = {census_mod.verdict_for(v) for v in (None, "-c", "-", "a.py", "/tmp/a.py", "pytest")}
    assert returned <= set(census_mod.VERDICTS)
    assert "absent" in census_mod.VERDICTS


# --- read-only ----------------------------------------------------------------------------------


def test_reading_a_result_does_not_touch_its_bytes():
    """Several results' sha256 are registrations; a census that rewrote one would break a pin."""
    if not HEADLINES.is_file():
        pytest.skip(f"{HEADLINES} is not in this checkout, so there is nothing to read")
    before = hashlib.sha256(HEADLINES.read_bytes()).hexdigest()
    row = census_mod.read_one("manifest_headlines")
    assert row["result"] == "manifest_headlines"
    assert hashlib.sha256(HEADLINES.read_bytes()).hexdigest() == before


# --- the boundary, read from git and not from a date ---------------------------------------------


def test_the_boundary_is_read_from_git_two_ways_that_agree():
    if not _git_history_is_there():
        pytest.skip("this checkout has no history for genomeos/manifest.py (a shallow clone)")
    b = census_mod.argv_boundary()
    assert b["field"].startswith("result_manifest.code.argv")
    assert b["the_two_queries_agree"] is True, b["found_by"]
    assert b["sha"] and len(b["sha"]) == 40
    assert b["date"]
    # the two queries are asked of git, so the record names them rather than a date in prose
    assert len(b["found_by"]) == 2
    assert all(v == b["sha"] for v in b["found_by"].values())


def test_the_boundary_commit_really_added_the_function():
    if not _git_history_is_there():
        pytest.skip("this checkout has no history for genomeos/manifest.py (a shallow clone)")
    sha = census_mod.argv_boundary()["sha"]
    assert sha
    diff = subprocess.run(
        ["git", "show", sha, "--", "genomeos/manifest.py"],
        capture_output=True,
        text=True,
        cwd=ROOT,
        check=True,
    ).stdout
    assert any(line.startswith("+def _argv(") for line in diff.splitlines())


# --- the split ----------------------------------------------------------------------------------


def test_a_result_no_commit_holds_is_named_unplaceable_and_not_folded():
    """The one rule the lane turns on: an unplaceable result gets its own class."""
    split = census_mod.split_by_boundary(["a_result_that_is_in_no_commit_at_all"])
    assert split["counts"] == {"never_committed_not_ignored": 1}
    assert split["the_finding"]["count"] == 0
    assert split["the_history"]["count"] == 0
    assert split["the_unplaceable"]["count"] == 1
    assert split["the_unplaceable"]["by_class"]["never_committed_not_ignored"] == [
        "a_result_that_is_in_no_commit_at_all"
    ]


def test_an_untracked_result_a_committed_rule_keeps_local_is_named_apart():
    """A committed .gitignore rule saying a result stays local is a decision, not an omission.

    `reader_*_chr*.json` is the rule, committed in .gitignore; `check-ignore` answers for the path
    whether or not the file is on disk, so this holds in a fresh checkout too.
    """
    ignored = "reader_a_cell_type_that_does_not_exist_chr1"
    split = census_mod.split_by_boundary([ignored, "a_result_that_is_in_no_commit_at_all"])
    assert split["counts"] == {
        "never_committed_and_git_ignored": 1,
        "never_committed_not_ignored": 1,
    }, split["counts"]
    by_class = split["the_unplaceable"]["by_class"]
    assert by_class["never_committed_and_git_ignored"] == [ignored]
    assert split["the_unplaceable"]["count"] == 2
    assert split["the_finding"]["count"] == 0 and split["the_history"]["count"] == 0


def test_the_finding_and_the_history_are_named_apart():
    split = census_mod.split_by_boundary([])
    assert split["the_finding"]["class"] == "after_the_field_existed"
    assert "FINDING" in census_mod.SIDES["after_the_field_existed"]
    assert split["the_history"]["class"] == "before_the_field_existed"
    assert "HISTORY" in census_mod.SIDES["before_the_field_existed"]


def test_the_side_is_not_read_from_the_date_field():
    assert "`date` field" in split_source()
    assert "never its `date` field" in census_mod.split_by_boundary([])["the_side_is_read_from"]


def split_source() -> str:
    return (ROOT / "scripts/entry_script_census.py").read_text()


def test_every_side_is_described_and_every_unplaceable_class_is_a_side():
    assert set(census_mod.UNPLACEABLE) <= set(census_mod.SIDES)
    assert "before_the_field_existed" not in census_mod.UNPLACEABLE
    assert "after_the_field_existed" not in census_mod.UNPLACEABLE
    # the two placeable sides plus every unplaceable class is the whole of SIDES: a class the split
    # can return but that is in neither group would be counted by nobody
    assert set(census_mod.SIDES) == {
        "before_the_field_existed",
        "after_the_field_existed",
        *census_mod.UNPLACEABLE,
    }


# --- the two stdin-written results --------------------------------------------------------------


def test_stdin_written_counts_them_and_claims_nothing_for_them():
    rows = [
        {
            "result": "manifest_headlines",
            "verdict": "stdin",
            "argv0": "-",
            "claims_own_code_is_committed": None,
            "bytes_state": "committed",
        },
        {
            "result": "crispri_split_audit",
            "verdict": "stdin",
            "argv0": "-",
            "claims_own_code_is_committed": None,
            "bytes_state": "committed",
        },
        {
            "result": "other",
            "verdict": "inside_repository",
            "argv0": "scripts/x.py",
            "claims_own_code_is_committed": True,
            "bytes_state": "committed",
        },
    ]
    block = census_mod.stdin_written(rows)
    assert block["count"] == 2
    assert sorted(r["result"] for r in block["results"]) == [
        "crispri_split_audit",
        "manifest_headlines",
    ]
    assert "neither is falsely certified" in block["reading"]
    assert "not rebuildable" in block["reading"] or "neither is \nrebuildable" in block["reading"].replace(
        "rebuildable", "\nrebuildable"
    )


# --- the population, counted rather than described ----------------------------------------------


def test_population_does_not_call_the_glob_committed():
    pop = census_mod.population(["manifest_headlines", "a_result_that_is_in_no_commit_at_all"])
    assert pop["files"] == 2
    assert pop["tracked_by_git"] + pop["untracked"] == 2
    assert "every committed result" in pop["it_is_not"]


# --- the wording that must not drift ------------------------------------------------------------


def test_the_registered_readings_are_carried_word_for_word():
    r = census_mod.READINGS
    assert r["falsely_certified"] == "0 of 1,109 committed results falsely certified"
    assert "It does NOT mean every result is rebuildable." in r["falsely_certified_means"]
    assert "NEVER PUBLISHED in that state" in " ".join(r["organised_chr"].split())
    assert "no `result_manifest` at all" in r["organised_chr"]
    assert "all 24 as `modified`, none `committed`" in r["organised_chr"]
    assert "do not repeat the stronger claim" in r["organised_chr_do_not_say"]
    assert "superseded at e1def2e" in r["organised_chr_do_not_say"]


def test_the_census_does_not_restate_the_stronger_organised_chr_claim():
    """The source may quote the claim only inside the sentence that forbids repeating it."""
    text = split_source()
    for line in text.splitlines():
        if "asserted own_code_is_committed" in line:
            assert "do not repeat" in text


# --- the manifest the census writes for itself ---------------------------------------------------


def test_the_manifest_meets_the_contract_and_counts_this_script_as_the_entry():
    if not HEADLINES.is_file():
        pytest.skip(f"{HEADLINES} is not in this checkout, so there is nothing to hash")
    m = census_mod.manifest_for(["manifest_headlines"])
    from genomeos import manifest as mf

    # `code` is the one field the writer never supplies: `save_result` stamps it, so the contract is
    # read on the stamped manifest and not on the half the caller builds.
    assert mf.validate(m) == ["missing code"]
    assert mf.validate(mf.stamp(m)) == []
    clean = m["code_cleanliness"]
    assert set(mf.CLEANLINESS_KEYS) <= set(clean)
    assert isinstance(clean["entry_script"], dict)
    assert clean["counting_path"], "the counting path must be computed, not empty"
    assert "scripts/entry_script_census.py" in clean["counting_path"]
    # one group entry naming every examined file with its own sha256
    (entry,) = m["inputs"]
    assert entry["group"] is True
    assert [x["path"] for x in entry["members"]] == ["data/results/manifest_headlines.json"]
    assert len(entry["sha256"]) == 64


def test_the_artefact_carries_the_split_the_readings_and_the_population():
    if not HEADLINES.is_file():
        pytest.skip(f"{HEADLINES} is not in this checkout, so there is nothing to examine")
    a = census_mod.artefact(["manifest_headlines"])
    assert a["population"]["files"] == 1
    assert a["all_results"]["examined"] == 1
    assert a["unrecorded_split"]["boundary"]["field"].startswith("result_manifest.code.argv")
    assert a["readings_carried_word_for_word"] == census_mod.READINGS
    assert a["examined_files"] == ["data/results/manifest_headlines.json"]
    assert "writes none of them" in a["nothing_was_written"]
    json.dumps(a)  # the payload a result is written from must be JSON
