# SPDX-License-Identifier: AGPL-3.0-or-later
"""A regenerated headline is renamed in place AND appended to `superseded`. Both, or neither counts.

THE CONVENTION, set by the supervisor on 2026-10-02. The ACTIVE rows of
`data/results/manifest_headlines.json` are exactly what README and the docs rest on, so when a result
is regenerated the row is RENAMED IN PLACE and the pointer moves with it in the same commit: a reader
of the active set must never be able to reach a figure the project no longer stands behind. AND every
regeneration ALSO APPENDS to `superseded`, with the old row name, the new one, the sha256 of each
file, the date of each and the reason -- so the rename leaves a trail instead of erasing one.

WHY IT IS TESTED AND NOT ONLY WRITTEN DOWN. `therapeutic_benchmark` and `constrained_unknown_targets`
were both renamed in place on 2026-10-02 (6d195a1 and bba1994) and neither left any trace: the old
names were simply gone, and the only record of what they had been was in git history and in the
`result_manifest.supersedes` block of the file that replaced each. A convention nothing checks is a
convention the next lane will not know about, and the lane after that will read the active set as
though it had always said this.

THE GATE NEVER READS `superseded`, and that is the point of keeping it a separate list.
`tests/test_headline_registry_matches_readme.py` and `tests/test_manifest_headlines.py` build their
sets from `rebuilt` and `pending` only, which is what makes the active set a statement about the
present rather than a pile of everything ever registered. These tests check that the history is
complete and consistent; nothing here can make a superseded figure count as a live one.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

RESULTS = Path("data/results")
RECORD = json.loads((RESULTS / "manifest_headlines.json").read_text())
ACTIVE = {e["result"] for e in RECORD["rebuilt"] + RECORD.get("pending", [])}
SUPERSEDED = RECORD.get("superseded", [])


def test_the_convention_is_written_in_the_file_a_lane_opens() -> None:
    """Not in a doc nobody opens: in the registry itself, next to the list it governs."""
    text = RECORD.get("superseded_convention", "")
    assert text, "manifest_headlines.json carries no superseded_convention"
    assert "RENAMED IN PLACE" in text
    assert "THE GATE NEVER READS" in text
    assert "NEVER restamped" in text or "NEVER touched" in text


def test_every_active_row_whose_name_says_it_replaced_one_has_a_superseded_entry() -> None:
    """The failure that prompted this: two rows renamed in place with no trail at all."""
    replaced = {name for name in ACTIVE if re.search(r"_v\d+$", name)}
    recorded = {e["registry_row_now"] for e in SUPERSEDED}
    # crispri_published and crispri_published_v2 are BOTH active -- two registered results, not a
    # supersession -- so a v-suffixed name is only required to have an entry when the name it
    # replaced is no longer active.
    missing = sorted(name for name in replaced - recorded if re.sub(r"_v\d+$", "", name) not in ACTIVE)
    assert missing == [], (
        f"these active rows look like regenerations whose predecessor is gone from the active set, "
        f"and `superseded` does not record them: {missing}. Renaming in place without appending "
        f"erases the trail the rename was supposed to leave"
    )


def test_no_superseded_entry_names_a_row_that_is_still_active() -> None:
    """A superseded name that is still live would mean the pointer did not actually move."""
    still_active = sorted(e["registry_row_was"] for e in SUPERSEDED if e["registry_row_was"] in ACTIVE)
    assert still_active == [], (
        f"these rows are recorded as superseded and are still in the active set: {still_active}. "
        f"The active rows are what README rests on, so a superseded figure may not sit among them"
    )


def test_every_superseded_entry_carries_both_names_both_shas_both_dates_and_a_reason() -> None:
    for e in SUPERSEDED:
        where = e.get("registry_row_was", "?")
        for key in ("registry_row_was", "registry_row_now", "superseded_on", "reason"):
            assert e.get(key), f"{where}: {key} is missing"
        assert "predecessor_was_a_registry_row" in e, (
            f"{where}: it is not stated whether the predecessor was ever a registry row, so a "
            f"reader cannot tell a moved pointer from a name that was new when it was coined"
        )
        for side in ("old", "new"):
            for key in ("file", "sha256", "date"):
                assert e[side].get(key), f"{where}: {side}.{key} is missing"
            assert re.fullmatch(r"[0-9a-f]{64}", e[side]["sha256"]), f"{where}: {side}.sha256"


def test_the_recorded_sha256_of_each_superseded_file_still_matches_its_bytes() -> None:
    """The sha256 is what makes the kept file a checkable record rather than a claim about one.

    It is also the guard against the one thing the convention forbids: a superseded file's bytes are
    never touched and it is never restamped. If one were, this fails and names it.
    """
    wrong = []
    for e in SUPERSEDED:
        for side in ("old", "new"):
            path = Path(e[side]["file"])
            if not path.is_file():
                wrong.append(f"{e[side]['file']}: absent")
                continue
            actual = hashlib.sha256(path.read_bytes()).hexdigest()
            if actual != e[side]["sha256"]:
                wrong.append(f"{e[side]['file']}: recorded {e[side]['sha256']}, on disk {actual}")
    assert wrong == [], (
        f"a superseded entry's sha256 no longer matches the file: {wrong}. A kept historical record "
        f"whose bytes moved is not a record; nothing may restamp or rewrite one"
    )


def test_a_supersession_does_not_by_itself_claim_the_figure_moved() -> None:
    """Two of the three 2026-10-02 regenerations were about what the manifest DECLARES, not the value.

    An entry that does not say so could be read as "the number was wrong", which would be a claim
    neither run makes. Every entry states whether the figure changed, and says it in words too.
    """
    for e in SUPERSEDED:
        assert "figure_changed" in e, f"{e['registry_row_was']}: figure_changed is not stated"
        assert e.get("figure_note"), f"{e['registry_row_was']}: figure_note is missing"
        if not e["figure_changed"]:
            assert e["old"].get("quoted_as") == e["new"].get("quoted_as"), (
                f"{e['registry_row_was']}: figure_changed is false but the quoted figures differ"
            )
        elif e["old"].get("quoted_as") is None:
            # A predecessor that was never a registry row has no registered wording, and one must not
            # be written after the fact. The entry says so in words rather than leaving a blank.
            assert e["old"].get("quoted_as_note"), (
                f"{e['registry_row_was']}: old.quoted_as is null and nothing says why"
            )


def test_the_node_containment_audit_row_records_its_2026_10_02_pass_and_its_limitation() -> None:
    """A gate pass that lives only on the work board does not exist: the board retires after 48 hours.

    The limitation is carried word for word. The result has no leaf in any set-aside class, so
    "nothing was set aside" is true of it by construction and says nothing about the cwd fix.
    """
    row = next(e for e in RECORD["rebuilt"] if e["result"] == "node_containment_audit")
    pass_ = row.get("fourth_pass_2026_10_02")
    assert pass_, "node_containment_audit's 2026-10-02 rebuild is not in the registry"
    assert pass_["inputs_declared_and_checked"] == "96 of 96"
    assert pass_["leaves"] == 3416
    assert re.fullmatch(r"[0-9a-f]{40}", pass_["verification_tool_at"])
    assert pass_["limitation"] == "vacuous: no set-aside class present"
    assert "says nothing" in pass_["limitation_means"]
