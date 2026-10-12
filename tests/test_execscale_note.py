# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-execfix: the additive note beside the execscale census, and what it is allowed to say.

`data/results/execscale_census.json` (`81d3949`) is a record of ONE tree and is NOT rewritten when
the tree moves -- a record that is edited when its subject changes is not a record. So the fact that
`71f9d47` legitimately moved data/organisms/human/haematopoiesis.bio after the census was taken is
written in a SIBLING record, `data/results/execscale_census_note.json`, by
`scripts/execscale_census_note.py` (`2f8830c`, committed before it ran). These are that record's
tests.

WHY ITS OWN FILE, and it is not a preference. `scripts/execscale_census_note.py` declares
tests/test_execscale_rescoped.py in its `OWN_CODE`, and `scripts/check_staged.py` refuses a result
whose own code was uncommitted when it was written. Appending these tests to that file therefore made
the pair uncommittable in either order: the record cannot be written while its own declared code is
dirty, and a test that reads a record must never be committed BEFORE the record (the ordering rule in
tests/test_inert_skip_guards.py, which cost two pushes on 2026-10-12). Editing the `OWN_CODE` line
would be a removal inside `2f8830c`'s window, which needs a `--force` this checkout does not grant.
A new file is the additive answer: it is not on the writer's own-code list, so the record stamps a
revision that does contain the code that wrote it, and record and tests still land together.

`tracked_paths.must_be_present` and not a skip, so while the record is absent these raise BY NAME.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

from genomeos import manifest as mf
from tests.tracked_paths import must_be_present

ROOT = Path(__file__).resolve().parents[1]
RECORD = "data/results/execscale_census.json"
#: the commit the census record was taken on, and the tree its digests are about
CENSUS_COMMIT = "81d3949"


def _blob(rev: str, path: str) -> bytes:
    return subprocess.run(["git", "show", f"{rev}:{path}"], cwd=ROOT, capture_output=True, check=True).stdout


def _census() -> dict:
    must_be_present(RECORD)
    return json.loads((ROOT / RECORD).read_text())


# ---- the additive note beside the census ---------------------------------------------------------
# These read data/results/execscale_census_note.json, which is committed in the same commit as they
# are (tests/test_inert_skip_guards.py's ordering rule), and `must_be_present` raises by name rather
# than skipping while it is absent, so nothing here can go quietly inert.

NOTE = "data/results/execscale_census_note.json"


def _note() -> dict:
    must_be_present(NOTE)
    return json.loads((ROOT / NOTE).read_text())


def test_the_note_is_a_SIBLING_and_the_census_it_notes_is_untouched():
    """The note's whole value is that the census was not rewritten, so that is asserted first."""
    note = _note()
    assert note["result_manifest"]["complete"] is True
    assert "not rewritten" in note["claim"]
    assert "nothing about RNA" in note["does_not_claim"]
    # the census record's own bytes are still the bytes of the commit that registered it
    census_blob = _blob(CENSUS_COMMIT, RECORD)
    assert hashlib.sha256((ROOT / RECORD).read_bytes()).hexdigest() == hashlib.sha256(census_blob).hexdigest()
    # and the note declares the census as an input rather than standing in for it
    paths = [entry.get("path") for entry in note["result_manifest"]["inputs"]]
    assert RECORD in paths


def test_the_note_names_the_ONE_programme_that_moved_and_the_COMMIT_that_moved_it():
    """What the over-wide byte pin could not say: the difference is a commit, not a slip."""
    note = _note()
    assert note["programs_in_the_census"] == 58
    moved = note["programs_that_moved"]
    assert [row["path"] for row in moved] == ["data/organisms/human/haematopoiesis.bio"]
    only = moved[0]
    assert (only["census_bytes"], only["now_bytes"]) == (20618, 20643)
    assert only["is_committed"] is True
    assert any("71f9d47" in line for line in only["commits_after_the_census"])
    # and no programme carries an uncommitted difference, which is the property this lane owes
    assert note["programs_with_an_uncommitted_difference"] == []
    assert note["programs_identical_to_the_census"] == 57


def test_the_note_records_the_group_digests_a_reader_may_be_HOLDING():
    """Three numbers, so the next reader does not rediscover the ordering footgun."""
    note = _note()
    names = [row["path"] for row in _census()["programs"]]
    assert note["group_sha256_in_the_census"] == _census()["group_sha256"]
    assert note["group_sha256_of_the_working_tree_in_the_census_order"] == mf.group_digest(names)[0]
    assert note["group_sha256_of_the_working_tree_sorted_by_name"] == mf.group_digest(sorted(names))[0]
    # the census's number is not the working tree's, and the note says which commit explains that
    assert note["group_sha256_in_the_census"] != note["group_sha256_of_the_working_tree_in_the_census_order"]


def test_the_notes_three_READINGS_say_what_this_lane_found_and_name_their_checks():
    note = _note()
    r = note["readings"]
    assert "reproducible" in r["the_census_record_is_sound"]
    assert "81d3949" in r["the_census_record_is_sound"]
    assert "same digest" in r["the_two_group_digest_implementations_agree"]
    assert "NOT removed" in r["the_two_group_digest_implementations_agree"]
    carries = r["a_group_digest_carries_the_LIST_as_well_as_the_CONTENT"]
    assert "ORDER GIVEN" in carries and "do not carry" in carries
    assert "files_entry" in carries
