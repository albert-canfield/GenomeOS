# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""A third set of published enhancer-gene loci (docs/LOCI-BENCHMARK.md section 21).

The expectation and registration tests need nothing but the code, so they pin the set before it is
ever scored. The rest read the committed result, so a locus that used to pass and stops passing
fails CI within a second.
"""

from __future__ import annotations

import pytest

from genomeos.benchmark import loci
from genomeos.benchmark.loci import MODEL_WINDOW, PANEL, READER_CELLS
from genomeos.benchmark.loci_candidates import CANDIDATES
from genomeos.benchmark.loci_third import (
    INTERVALS as INTERVALS_NAME,
)
from genomeos.benchmark.loci_third import (
    NEEDED_LOCI_HUNK,
    POSITIVE_CONTROL,
    PREREGISTRATION,
    REPRESSORS,
    THIRD,
    TRAPS,
    keep_out_all_three_sets,
    register_stated_intervals,
)
from genomeos.results import load_result

RESULT = load_result("loci_third")
INTERVALS = load_result(INTERVALS_NAME)
needs_result = pytest.mark.skipif(not RESULT, reason="third-set result not present")


# ------------------------------------------------------------------- the expectations, offline
def test_every_locus_is_written_down_with_its_sources():
    for e in THIRD:
        assert e.window[0] <= e.element[0] < e.element[1] <= e.window[1], e.locus
        assert e.targets and e.citations and e.answer_from, e.locus
        assert set(e.cells) <= set(READER_CELLS), e.locus
        assert e.direction in ("activates", "represses", "coding", "none"), e.locus
        if e.element_source == "stated":
            assert e.element_citation, f"{e.locus}: a stated interval must carry its own citation"
        for v in e.variants:
            if v.get("pos"):
                assert e.window[0] <= v["pos"] <= e.window[1], f"{e.locus} {v.get('rsid')}"


def test_no_locus_repeats_an_element_either_of_the_other_two_sets_holds():
    """A third FRAME, not a third reading. Targets may be shared on purpose; elements may not."""
    others = (*PANEL, *CANDIDATES)
    names = {e.locus for e in others}
    for e in THIRD:
        assert e.locus not in names, e.locus
        clash = [
            p
            for p in others
            if p.chrom == e.chrom and p.element[0] < e.element[1] and p.element[1] > e.element[0]
        ]
        assert not clash, f"{e.locus} overlaps {[p.locus for p in clash]}"


def test_the_set_carries_the_two_properties_it_exists_for():
    """Section 19 could test neither: a published repressive direction, and a range of lengths."""
    repressors = {e.locus for e in THIRD if e.direction == "represses"}
    assert repressors == set(REPRESSORS), repressors
    assert len(repressors) >= 2, "the set exists to ask the direction question the other way"
    assert all(e.direction == "activates" for e in CANDIDATES), (
        "if the second set ever gains a repressor this set's reason for existing has changed"
    )
    lengths = sorted(e.element[1] - e.element[0] for e in THIRD)
    assert lengths[0] <= 500 and lengths[-1] >= 10_000, (
        f"section 20 asked for a set spanning the panel's length range, got {lengths}"
    )


def test_a_deliberate_target_overlap_is_declared_rather_than_discovered():
    """Three targets are already panel targets. That is the design, so it is written down."""
    panel_targets = {t for e in PANEL for t in e.targets}
    shared = sorted({e.locus for e in THIRD if set(e.targets) & panel_targets})
    assert shared == sorted(PREREGISTRATION["deliberate_target_overlaps"]["loci"]), shared


def test_the_traps_are_geometry_written_down_before_scoring():
    """A nearest-TSS trap is arithmetic over GENCODE, so it belongs in the file, not in the result."""
    named = {e.locus: e.nearest_gene_trap for e in THIRD if e.nearest_gene_trap}
    assert set(named) == set(TRAPS), (named, set(TRAPS))
    for locus, trap in named.items():
        assert TRAPS[locus]["trap"] == trap
        assert TRAPS[locus]["trap_distance"] < TRAPS[locus]["target_distance"], locus


def test_the_registration_fixes_the_rules_before_the_run():
    for key in (
        "hit_rules",
        "denominators",
        "negative_controls",
        "predictions",
        "reach_filter",
        "aims_that_failed_before_the_run",
    ):
        assert PREREGISTRATION[key], key
    assert PREREGISTRATION["loci"] == len(THIRD)
    # the repressor arm's n is declared too small to decide anything BEFORE the number exists
    assert "n = 2" in PREREGISTRATION["predictions"]["direction_repressors"]
    assert POSITIVE_CONTROL in PREREGISTRATION["what_would_falsify_the_run"]


def test_the_positive_control_is_the_gene_s_own_promoter():
    e = next(x for x in THIRD if x.locus == POSITIVE_CONTROL)
    assert e.element[0] <= 1_295_113 <= e.element[1], "the C228T hotspot must be inside the interval"
    assert e.distance is not None and e.distance < 1_000


# ------------------------------------------------------------- the two shims, and what they promise
def test_the_stated_interval_registration_is_idempotent_and_only_appends():
    before = loci.STATED_INTERVAL_RESULTS
    try:
        register_stated_intervals()
        once = loci.STATED_INTERVAL_RESULTS
        register_stated_intervals()
        assert once == loci.STATED_INTERVAL_RESULTS
        assert set(before) <= set(once) and INTERVALS_NAME in once
        assert NEEDED_LOCI_HUNK["to"] == once[: len(NEEDED_LOCI_HUNK["to"])] or INTERVALS_NAME in once
    finally:
        loci.STATED_INTERVAL_RESULTS = before


def test_the_keep_out_extension_covers_all_three_sets_and_is_restored():
    before = loci.PANEL
    with keep_out_all_three_sets() as extended:
        names = {e.locus for e in extended}
        assert {e.locus for e in PANEL} <= names
        assert {e.locus for e in CANDIDATES} <= names
        assert {e.locus for e in THIRD} <= names
    assert loci.PANEL is before, "the global must be put back, or the next panel draws wrong"


def test_the_keep_out_extension_is_restored_after_an_exception():
    before = loci.PANEL
    with pytest.raises(RuntimeError), keep_out_all_three_sets():
        raise RuntimeError("boom")
    assert loci.PANEL is before


# ------------------------------------------------------------------------- the committed result
@needs_result
def test_the_result_grades_what_the_reach_filter_kept():
    names = [r["locus"] for r in RESULT["loci"]]
    assert names, "the run graded nothing"
    assert set(names) <= {e.locus for e in THIRD}


@needs_result
def test_the_reach_filter_is_counted_rather_than_asserted():
    reach = RESULT["reach"]
    assert reach["loci"] == len(THIRD)
    assert reach["died_at_the_reach_filter"] == len(THIRD) - len(RESULT["loci"])
    for row in reach["which"]:
        assert row["distance_kb"] and row["distance_kb"] > MODEL_WINDOW // 2 // 1000


@needs_result
def test_every_request_is_accounted_for_per_locus():
    """An unaskable locus is never asked, an already-scored element is never re-asked, and every
    request that was spent landed on a graded locus and is visible as a scored element there."""
    rows = {r["locus"]: r for r in RESULT["plan"]}
    for row in rows.values():
        if not row["askable"]:
            assert row["requests"] == 0, f"{row['locus']}: a request on an unaskable locus"
        if row["already_scored_over_element"]:
            assert row["requests"] == 0, f"{row['locus']}: asked for what is already scored"
    for got in (INTERVALS or {}).get("read", []):
        row = rows[got["locus"]]
        assert row["graded"], f"{got['locus']}: a request was spent on a locus that is not graded"
        assert row["already_scored_over_element"], f"{got['locus']}: the request did not land"


@needs_result
def test_the_positive_control_passed_or_the_run_says_it_did_not():
    """A miss here is not a miss, it is a broken reading, and the result has to say which."""
    pc = RESULT["positive_control"]
    assert pc["locus"] == POSITIVE_CONTROL and pc["ran"]
    assert pc["passed"] is not None


@needs_result
def test_the_two_directions_are_reported_apart_and_never_pooled():
    d = RESULT["directions"]
    assert set(d["activators"]) == {"k", "n"} and set(d["repressors"]) == {"k", "n"}
    assert "pooled" not in {k.lower() for k in d} - {"not_pooled"}
    per = {r["locus"]: r for r in d["per_locus"]}
    for name in REPRESSORS:
        assert per[name]["published"] == "represses"
        assert per[name]["repressor"] is True


@needs_result
def test_no_hit_is_annotation_lookup_alone():
    assert RESULT["aggregate"]["target_only_looked_up"] == []


@needs_result
def test_three_frames_are_named_and_no_rate_pools_them():
    a = RESULT["aggregate"]
    for key in ("target_derived", "target_derived_where_the_model_could_answer"):
        assert a[key]["n"] and "misses" in a[key]
    beside = RESULT["beside_the_other_two"]
    assert set(beside["frames"]) == {
        "the_seventeen_panel",
        "the_nine_candidates",
        "the_third_set",
    }
    assert set(beside["target_derived"]) == set(beside["frames"])


@needs_result
def test_presence_is_counted_before_coverage_is_made_a_stratum():
    """Section 19's trap: a null under a coverage stratum on an input nobody ever had to buy is
    arithmetic, not evidence. `input_presence` says which it is, and it runs first."""
    ac = RESULT["against_controls"]
    assert "input_presence" in ac and "direction_given_coverage" in ac
    kinds = {v["kind"] for v in ac["input_presence"].values()}
    assert kinds <= {"free", "bought"} and kinds


@needs_result
def test_the_registration_travels_with_the_result():
    assert RESULT["preregistration"]["written"] == PREREGISTRATION["written"]
    assert RESULT["preregistration"]["predictions"] == PREREGISTRATION["predictions"]
