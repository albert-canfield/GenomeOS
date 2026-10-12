# SPDX-License-Identifier: AGPL-3.0-or-later
"""globcensus guards. The one that matters is that the AST scan sees what a one-line grep cannot."""

from __future__ import annotations

import subprocess
from pathlib import Path

from genomeos.attribution import globcensus as gc

ROOT = Path(__file__).resolve().parent.parent


def test_the_ast_scan_finds_the_call_a_one_line_grep_cannot_see() -> None:
    """The brief's count was 17 from a one-line grep; the AST sees 18. This is why."""
    calls = gc.files_entry_calls(ROOT)
    filtered = {(c["path"], c["line"]) for c in calls if c["exists_filtered"]}
    out = subprocess.run(
        ["git", "grep", "-n", "files_entry", "--", "genomeos/", "scripts/"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    one_line = set()
    for line in out.splitlines():
        if "exists()" in line:
            f, n, _ = line.split(":", 2)
            one_line.add((f, int(n)))
    missed = filtered - one_line
    assert missed, "if a one-line grep now sees every call, this test no longer measures anything"
    assert ("scripts/trio_q100_phase.py", 72) in missed, missed
    spans = {c["spans_lines"] for c in calls if (c["path"], c["line"]) in missed}
    assert max(spans) > 1, "the missed call must be the multi-line one, which is the whole reason"


def test_the_exists_filtered_sites_are_a_minority_of_all_files_entry_calls() -> None:
    calls = gc.files_entry_calls(ROOT)
    filtered = [c for c in calls if c["exists_filtered"]]
    assert len(filtered) >= 18, len(filtered)
    assert len(filtered) < len(calls), "the census would be vacuous if every call were filtered"


def test_the_reader_dnase_group_declares_the_same_count_in_every_result() -> None:
    """The census's NEGATIVE, pinned: the shrink has not bitten this group in any committed result."""
    groups = gc.declared_groups()
    reader = groups.get("reader v1 DNase narrowPeak files")
    assert reader, "the group this census was opened on is no longer declared anywhere"
    counts = {n for _, n in reader}
    assert counts == {312}, counts
    assert "reader v1 DNase narrowPeak files" not in gc.varying(groups)


def test_varying_names_only_groups_whose_count_really_differs() -> None:
    groups = gc.declared_groups()
    var = gc.varying(groups)
    for name, counts in var.items():
        seen = {n for _, n in groups[name]}
        assert len(seen) > 1, name
        assert sorted(counts, key=lambda x: (x is None, x)) == counts
    for name, rows in groups.items():
        if len({n for _, n in rows}) == 1:
            assert name not in var, name


def test_PLANTED_a_group_constant_at_a_shrunken_value_is_INVISIBLE_here() -> None:
    """The limit DOES_NOT_CLAIM states, as a test rather than a sentence.

    A group that is the same wrong size in every result looks exactly like one that is the same right
    size. If this ever fails, `varying` has gained a power it does not have and the caveat is stale.
    """
    planted = {"a group shrunken identically everywhere": [("r1", 13), ("r2", 13), ("r3", 13)]}
    assert gc.varying(planted) == {}, "a uniformly shrunken group must be invisible to a count test"


def test_the_census_says_what_it_does_not_claim() -> None:
    assert "not that any committed result IS wrong" in gc.DOES_NOT_CLAIM
    assert "constant at a shrunken value" in gc.DOES_NOT_CLAIM


# ------------------------------------------------------------------------------------------ stage 2


def _group(result: str, name_fragment: str) -> dict:
    import json

    rm = json.loads((ROOT / f"data/results/{result}.json").read_text())["result_manifest"]
    hits = [
        i
        for i in rm["inputs"]
        if isinstance(i, dict) and i.get("group") is True and name_fragment in str(i.get("path"))
    ]
    assert len(hits) == 1, f"{result}: {len(hits)} groups match {name_fragment!r}"
    return {
        "group": hits[0],
        "scope": rm["parameters"].get("chromosomes"),
        "exclusions": rm.get("exclusions"),
    }


def test_every_recorded_shrink_is_re_measured_from_its_own_result() -> None:
    """The finding is pinned to the artefacts, not to the list in the module."""
    import re

    for row in gc.SHRINKS_FOUND:
        got = _group(row["result"], row["group"].split(" per-chromosome")[0])
        assert got["group"]["files"] == row["declared"], row["result"]
        assert len(got["scope"]) == row["scope"], row["result"]
        have = {c for m in got["group"]["members"] for c in re.findall(r"chr[0-9XYM]+", str(m["path"]))}
        missing = {str(c) for c in got["scope"]} - have
        assert missing == set(row["missing"]), f"{row['result']}: {sorted(missing)} vs {row['missing']}"


def test_no_shrunken_result_names_the_missing_chromosome_in_its_exclusions() -> None:
    """Condition (3), and the part that makes the defect SILENT rather than merely narrow."""
    for row in gc.SHRINKS_FOUND:
        got = _group(row["result"], row["group"].split(" per-chromosome")[0])
        text = " ".join(str(x) for x in (got["exclusions"] or []))
        named = tuple(c for c in row["missing"] if c in text)
        assert named == row.get("named_in_exclusions", ()), (
            f"{row['result']}: exclusions name {named}, the record says {row.get('named_in_exclusions', ())}"
        )


def test_the_shrink_test_needs_all_three_conditions() -> None:
    assert gc.is_shrink(22, 24, True, False)
    assert not gc.is_shrink(24, 24, True, False), "a full declaration is never a shrink"
    assert not gc.is_shrink(22, 24, False, False), "a store that merely grew is not a shrink"
    assert not gc.is_shrink(22, 24, True, True), "a named exclusion is a narrower analysis, not silence"


def test_stage_2_separates_what_it_shows_from_what_it_does_not() -> None:
    assert "DEFECT IS REALISED" in gc.STAGE_2_READING
    assert "two stores" in gc.STAGE_2_READING and "eight days apart" in gc.STAGE_2_READING
    assert "not that any of the six results' NUMBERS are wrong" in gc.STAGE_2_DOES_NOT_CLAIM
    assert "chrY is routinely empty" in gc.STAGE_2_DOES_NOT_CLAIM
    assert "declare no chromosome scope at all" in gc.STAGE_2_DOES_NOT_CLAIM


# -------------------------------------------------------------------- stage 2's refutation, 2026-10-12


def test_the_realised_claim_is_refuted_and_the_summary_says_so() -> None:
    assert "DEFECT IS REALISED" in gc.STAGE_2_READING, "the wrong claim is kept standing, not edited away"
    assert "REALISED IN 0 OF 6" in gc.STAGE_2_REFUTED
    assert gc.STAGE_2_SUMMARY == "the mechanism is real in 18 places and realised in 0 of the 6 checked"
    module_text = (ROOT / "genomeos/attribution/globcensus.py").read_text()
    assert "9a68bf2" in module_text, "the fold that refuted it must be named"


def test_the_fourth_condition_is_what_the_three_could_not_see() -> None:
    for phrase in ("LABELLING ARTEFACT", "EMPTY CHROMOSOME", "Only opening the member settles it"):
        assert phrase in gc.FOURTH_CONDITION, phrase


def test_PLANTED_the_three_condition_test_PASSES_every_case_that_was_not_a_shrink() -> None:
    """Why the old test was insufficient, as a measurement rather than a concession.

    A labelling artefact and an empty chromosome both satisfy all three original conditions. If this
    ever fails, SHRINK_TEST has gained a discrimination it did not have.
    """
    assert gc.is_shrink(22, 24, True, False), "a labelling artefact passes the 3-condition test"
    assert gc.is_shrink(23, 24, True, False), "an empty chromosome passes it too"
    assert not gc.is_shrink_v2(22, 24, True, False, False), "the 4th condition refuses both"
    assert gc.is_shrink_v2(22, 24, True, False, True), "and still admits a real shrink"


def test_chr21_really_is_present_in_another_group_of_the_same_result() -> None:
    """The prefix-grouping artefact, verified from the artefact rather than accepted."""
    import json
    import re

    rm = json.loads((ROOT / "data/results/closure_cell_inventory.json").read_text())["result_manifest"]
    where = {}
    for i in rm["inputs"]:
        if isinstance(i, dict) and i.get("group") is True:
            chroms = {
                c for m in (i.get("members") or []) for c in re.findall(r"chr[0-9XYM]+", str(m.get("path")))
            }
            where[str(i.get("path"))[:40]] = "chr21" in chroms
    assert any(where.values()), "chr21 is in no group at all, so it was not merely mislabelled"
    assert not all(where.values()), "chr21 is in every group, so there was nothing to explain"
