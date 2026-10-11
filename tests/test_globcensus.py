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
