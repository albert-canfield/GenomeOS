# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Eleven more published enhancer-gene loci (docs/LOCI-BENCHMARK.md section 19).

The expectation and filter tests need nothing but the code. The rest read the committed result, so a
candidate that used to pass and stops passing fails CI within a second.
"""

from __future__ import annotations

import pytest

from genomeos.benchmark.loci import MODEL_WINDOW, PANEL, READER_CELLS, Chromosome, read_reach
from genomeos.benchmark.loci_candidates import (
    CANDIDATES,
    EXHIBITS_ONLY,
    PREREGISTRATION,
)
from genomeos.results import load_result

RESULT = load_result("loci_candidates")
INTERVALS = load_result("loci_candidate_intervals")
needs_result = pytest.mark.skipif(not RESULT, reason="candidate result not present")


class _Locus:
    def __init__(self, start: int, end: int) -> None:
        self.start, self.end = start, end


class _Gene:
    def __init__(self, symbol: str, start: int, end: int) -> None:
        self.symbol, self.locus = symbol, _Locus(start, end)


class _Annotation:
    def __init__(self, genes: dict) -> None:
        self.genes = genes


class _Chromosome:
    """Only what `read_reach` touches: GENCODE gene bodies."""

    def __init__(self, genes: dict) -> None:
        self.annotation = _Annotation(genes)


# ------------------------------------------------------- the filter, and the way it must not be written
def test_the_reach_filter_is_gene_body_overlap_and_not_tss_distance():
    """RUNX1 is 1.22 Mb long, so an element 600 kb from its start still sits inside it.

    On the finished chr21 sweep, 79 of 5,174 elements with a coding target name a gene whose TSS is
    524 kb or more away - up to 986 kb - and every one is a gene whose body overlaps the window.
    A filter written in TSS distance would drop them all, and would fail in the direction that looks
    safe. This pins the semantics so a rewrite cannot pass.
    """
    from genomeos.benchmark.loci import Expect

    element = (35_000_000, 35_000_500)
    mid = sum(element) // 2
    runx1 = _Gene("RUNX1", mid - 600_000, mid + 620_000)  # TSS 600 kb away, body over the element
    far = _Gene("FAR", mid + 900_000, mid + 910_000)  # wholly outside the input
    ch = _Chromosome({"RUNX1": runx1, "FAR": far})
    tss_distance = abs(runx1.locus.start - mid)
    assert tss_distance > MODEL_WINDOW // 2, "the case is only a case if the TSS is out of range"

    def expect(target: str):
        return Expect(
            locus="t",
            chrom="chr21",
            element=element,
            window=element,
            classes=("program",),
            targets=(target,),
            direction="activates",
            tissues=(),
        )

    assert read_reach(ch, expect("RUNX1"))["askable"], "a body over the window is in reach"
    assert not read_reach(ch, expect("FAR"))["askable"], "a body outside the window is not"


def test_the_registry_lookup_finds_an_element_that_starts_before_the_window():
    """`load_ccres` is not in start order, so bisecting the raw list silently misses elements.

    It reported 1 of the 12 cCREs over the MC1R window and 0 of the 1 over the ZRS, which is where
    section 10's `no cCRE over the published element` count came from. The lookup sorts its own copy.
    """
    ch = Chromosome.__new__(Chromosome)
    els = [_Gene("late", 9_000, 9_500), _Gene("early", 1_000, 5_100), _Gene("middle", 4_900, 5_000)]
    for e in els:  # the shape ccres_in reads: .start and .end
        e.start, e.end = e.locus.start, e.locus.end
    ch.ccres = els  # unsorted, as it comes
    ch._ccre_sorted = sorted(els, key=lambda c: c.start)
    ch._ccre_starts = [c.start for c in ch._ccre_sorted]
    got = {c.symbol for c in ch.ccres_in(5_000, 6_000)}
    assert got == {"early"}, f"an element starting before the window must be found, got {got}"


# ------------------------------------------------------------------- the expectations, offline
def test_every_candidate_is_written_down_with_its_sources():
    for e in CANDIDATES:
        assert e.window[0] <= e.element[0] < e.element[1] <= e.window[1], e.locus
        assert e.targets and e.citations and e.answer_from, e.locus
        assert set(e.cells) <= set(READER_CELLS), e.locus
        assert e.direction in ("activates", "represses", "coding", "none"), e.locus
        if e.element_source == "stated":
            assert e.element_citation, f"{e.locus}: a stated interval must carry its own citation"
        for v in e.variants:
            if v.get("pos"):
                assert e.window[0] <= v["pos"] <= e.window[1], f"{e.locus} {v.get('rsid')}"


def test_no_candidate_repeats_a_locus_the_panel_already_holds():
    names = {e.locus for e in PANEL}
    for e in CANDIDATES:
        assert e.locus not in names, e.locus
        clash = [
            p
            for p in PANEL
            if p.chrom == e.chrom and p.element[0] < e.element[1] and p.element[1] > e.element[0]
        ]
        assert not clash, f"{e.locus} overlaps the panel's {[p.locus for p in clash]}"


def test_an_approximately_placed_candidate_is_never_graded():
    """Section 13 withdrew a case for being placed on an offset rather than a coordinate.

    The two kept here are kept for the reach question alone, which an offset can still answer.
    """
    assert {e.locus for e in CANDIDATES} >= EXHIBITS_ONLY
    for name in EXHIBITS_ONLY:
        e = next(x for x in CANDIDATES if x.locus == name)
        assert e.element_source == "stated" and "EXHIBIT ONLY" in (e.element_citation or "")


def test_the_registration_fixes_the_rules_before_the_run():
    for key in ("hit_rules", "denominators", "negative_controls", "predictions", "reach_filter"):
        assert PREREGISTRATION[key], key
    assert "body" in PREREGISTRATION["reach_filter"]["why_body_not_tss"].lower()


# ------------------------------------------------------------------------- the committed result
@needs_result
def test_the_result_is_the_graded_set_and_the_exhibits_are_not_in_it():
    names = [r["locus"] for r in RESULT["loci"]]
    assert names, "the run graded nothing"
    assert not (set(names) & EXHIBITS_ONLY), "an exhibit was graded"
    assert set(names) <= {e.locus for e in CANDIDATES}


@needs_result
def test_the_reach_filter_is_counted_rather_than_asserted():
    reach = RESULT["reach"]
    assert reach["candidates"] == len(CANDIDATES)
    assert reach["died_at_the_reach_filter"] == len(CANDIDATES) - len(RESULT["loci"])
    for row in reach["loci"]:
        assert row["distance_kb"] and row["distance_kb"] > MODEL_WINDOW // 2 // 1000


@needs_result
def test_every_request_is_accounted_for_per_locus():
    """The plan counts the requests before any exist, and a locus out of reach is never asked.

    The plan is re-derived on every run, so once a stated interval has been scored its row shows
    nothing left to spend. What must hold is the accounting: an unaskable locus is never asked, a
    locus whose element the sweep already covers is never asked, and every request that was spent
    landed on a graded locus and is now visible as a scored element there.
    """
    rows = {r["locus"]: r for r in RESULT["plan"]}
    for row in rows.values():
        if not row["askable"] or row["exhibit_only"]:
            assert row["requests"] == 0, f"{row['locus']}: a request on a locus that is not graded"
        if row["already_scored_over_element"]:
            assert row["requests"] == 0, f"{row['locus']}: asked for what is already scored"
    for got in (INTERVALS or {}).get("read", []):
        row = rows[got["locus"]]
        assert row["graded"], f"{got['locus']}: a request was spent on a locus that is not graded"
        assert row["already_scored_over_element"], f"{got['locus']}: the request did not land"
        assert not got["substituted_annotated_element"], (
            f"{got['locus']}: the scorer answered about an annotated element, not the stated one"
        )


@needs_result
def test_no_candidate_hit_is_annotation_lookup_alone():
    assert RESULT["aggregate"]["target_only_looked_up"] == []


@needs_result
def test_both_denominators_are_reported_and_no_third_one_is_invented():
    a = RESULT["aggregate"]
    for key in ("target_derived", "target_derived_where_the_model_could_answer"):
        assert a[key]["n"] and "misses" in a[key]
    beside = RESULT["beside_the_panel"]
    assert set(beside["target_derived"]) == {"candidates", "panel"}


@needs_result
def test_the_direction_claim_is_reported_against_its_coverage():
    """Section 8's retraction, applied to this panel: a direction is a statement about being scored.

    Standardising on GC, distance and constraint leaves the artefact in place, so the run must carry
    the reading that holds coverage fixed as well.
    """
    ac = RESULT["against_controls"]
    assert "direction_given_coverage" in ac
    conditioned = RESULT["negative_controls"]["given_deletion_data"]
    assert conditioned["negatives"]["without_deletion_data"]["direction"]["k"] == 0, (
        "a window with no scored element produced a direction, which is impossible"
    )


@needs_result
def test_the_registration_travels_with_the_result():
    assert RESULT["preregistration"]["written"] == PREREGISTRATION["written"]
    assert RESULT["preregistration"]["predictions"] == PREREGISTRATION["predictions"]
