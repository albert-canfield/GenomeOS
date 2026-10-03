"""The roadmap truth pass asserts only what cannot go stale, and these tests pin that.

Two invariants hold over docs/ROADMAP.md and must keep holding: every commit it cites
is an ancestor of HEAD or carries a recorded reason why not, and every result it names
in a Data bullet exists on disk. Both are about existence, so a failure is a defect and
never merely an old figure.

The counts the script prints are deliberately *not* asserted: areas B, C and J state
figures that were correct on their own dates and have since been overtaken, and a test
that failed on those would force a lane to edit a row it does not own, which is how a
true-as-of-its-date sentence gets overwritten instead of corrected additively. The
tests below therefore check that the counter itself is right — on a fixture whose
answer is known — rather than that the document agrees with it.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import roadmap_truth_pass as tp  # noqa: E402


@pytest.fixture(scope="module")
def found() -> dict:
    return tp.run()


def test_every_cited_commit_is_an_ancestor_or_has_a_recorded_reason(found: dict) -> None:
    """A cited commit the history does not carry is a claim with nothing behind it."""
    commits = found["commits"]
    assert commits["cited"] > 400, "the sha extractor stopped finding citations"
    unexplained = commits["unexplained"]
    assert not unexplained, (
        "cited but not an ancestor of HEAD and not in NOT_ANCESTOR_EXPECTED — either "
        "the work never landed on dev or the row cites a relabelled commit: "
        f"{[u['commit'] for u in unexplained]}"
    )


def test_the_reasons_list_stays_pruned(found: dict) -> None:
    """A commit that has since merged must leave the allowlist, or it hides the next one."""
    stray = {s["commit"] for s in found["commits"]["not ancestors"]}
    stale_entries = sorted(set(tp.NOT_ANCESTOR_EXPECTED) - stray)
    assert not stale_entries, (
        "these are ancestors of HEAD now, so their exemption is dead weight and would "
        f"mask a real stray: {stale_entries}"
    )


def test_every_result_named_in_a_data_bullet_exists(found: dict) -> None:
    named = [r for r in found["data bullets"] if r["kind"] == "file"]
    assert named, "the Data-bullet parser found no named results"
    missing = [(r["line"], r["claim"]) for r in named if not r["exists"]]
    assert not missing, f"named in a Data bullet but absent from data/results: {missing}"


def test_milestone_states_have_one_reading(found: dict) -> None:
    """The ✅/◑ markers are parsed by the same function the web status card uses."""
    ms = found["milestones"]
    assert ms["total"] == 7, f"the milestone table changed shape: {ms['milestones']}"
    states = {m["milestone"]: m["state"] for m in ms["milestones"]}
    assert states["1.3 the 98%"] == "partial", (
        "1.3 is held not met by area I's unknown-scoring measurement; a ✅ here would "
        "need that measurement answered, not a marker change"
    )


def test_the_chromosome_counter_counts_chromosomes_not_files(tmp_path, monkeypatch) -> None:
    """The counter's own failure mode: rmsk ships two files per chromosome.

    Written because the first version of this script matched ``<stem>chrN.json`` against
    a stem that already ended in ``chr``, and so reported 0 for every family while
    looking like a clean pass. A counter that can only read low is worse than none.
    """
    for name in (
        "rmsk_chr1.json",
        "rmsk_chr1.bed.gz",
        "rmsk_chr2.json",
        "rmsk_chr2.bed.gz",
        "rmsk_chrX.json",
        "rmsk_chr1_contact.json",  # a different family under the same prefix
        "unrelated.json",
    ):
        (tmp_path / name).write_text("{}")
    monkeypatch.setattr(tp, "RESULTS", tmp_path)
    chromosomes, files = tp._on_disk("rmsk_chr")
    assert chromosomes == 3, "chr1, chr2 and chrX, counted once each"
    assert files == 6, "every file carrying the prefix, variants included"


def test_a_family_with_no_files_reads_zero(tmp_path, monkeypatch) -> None:
    (tmp_path / "other_chr1.json").write_text("{}")
    monkeypatch.setattr(tp, "RESULTS", tmp_path)
    assert tp._on_disk("absent_chr") == (0, 0)


@pytest.mark.parametrize(
    ("after", "expected"),
    [
        (" (24 of 24), `budget_genome_wide`", 24),
        (" **25 of 25** plus `proteome_genome_wide`", 25),
        (" (18)", 18),
        (" rows (18)", None),  # the count is not this name's; "rows" intervenes
        (", and the rest", None),
    ],
)
def test_the_stated_count_is_read_in_both_forms(after: str, expected: int | None) -> None:
    """Both forms the document uses, and no count where the prose does not give one."""
    match = tp.STATED_COUNT.search(after)
    got = int(match.group(1)) if match else None
    assert got == expected
