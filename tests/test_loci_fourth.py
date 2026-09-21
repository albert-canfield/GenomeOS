# SPDX-License-Identifier: AGPL-3.0-or-later
"""The fourth frame's rule, tested before the rule is run.

Every test here names the failure it prevents. The frame's whole claim is that it was fixed before
anybody saw which loci it returns, so the tests are about the rule's arithmetic and about the ways a
frame can quietly become a choice: a locus kept for the wrong reason, a cap that acts as a quality
filter, a control that scores into a rate, a keep-out that leaks an earlier frame's locus in.
"""

from __future__ import annotations

import pytest

from genomeos.attribution.crispri import Pair
from genomeos.benchmark import loci
from genomeos.benchmark import loci_fourth as fourth


def pair(chrom="chr1", start=1_000_000, end=1_000_500, gene="AAA", cell="K562", regulated=True, dist=50_000):
    return Pair(
        chrom=chrom,
        start=start,
        end=end,
        gene=gene,
        cell=cell,
        dataset="Gasperini2019",
        distance=dist,
        dhs=1.0,
        h3k27ac=1.0,
        regulated=regulated,
    )


def element(chrom="chr1", start=1_000_000, end=1_000_500, targets=("AAA",), distance=50_000, cell="K562"):
    return {
        "chrom": chrom,
        "start": start,
        "end": end,
        "cell": cell,
        "mid": (start + end) // 2,
        "targets": list(targets),
        "distances": {t: distance for t in targets},
        "datasets": ["Gasperini2019"],
        "tested_genes": len(targets),
        "distance": distance,
    }


# ------------------------------------------------------------------ the registration itself
def test_the_registration_states_the_rule_the_budget_and_what_a_failure_would_be():
    """A registration missing any of these is a formality, which is the artefact it exists to stop."""
    p = fourth.PREREGISTRATION
    assert "before" in p["written"]
    assert p["the_rule"]["steps_in_order"] and len(p["the_rule"]["steps_in_order"]) >= 8
    assert p["requests"]["budget"] == fourth.REQUEST_BUDGET
    for outcome in ("derived_rate_holds_near_the_other_three", "derived_rate_falls_to_the_chance_floor"):
        assert p["what_each_outcome_would_mean"][outcome]
    assert "positive control" in p["what_would_count_as_a_failure_of_the_RUN_rather_than_a_result"]
    # the reach prediction must be a NUMBER registered in advance, and it must not be the third set's
    assert "1 to 4" in p["predictions"]["reach_fatalities"]
    assert "0 of 9" in p["predictions"]["reach_fatalities"]
    # and the direction axis must be declared unaskable before the draw rather than after it
    assert "NOT ASKABLE" in p["predictions"]["direction"]


def test_the_frame_is_a_rule_and_not_a_list():
    """A hand-written tuple of loci here would mean the frame was picked after looking, whatever the
    docstring said. The module must not contain one."""
    assert not any(
        isinstance(v, tuple) and v and isinstance(v[0], loci.Expect)
        for k, v in vars(fourth).items()
        if k.isupper()
    )


# ------------------------------------------------------------------------------ step 3: grouping
def test_group_elements_groups_by_element_and_cell_and_keeps_only_regulated_genes():
    rows = fourth.group_elements(
        [
            pair(gene="AAA", regulated=True, dist=80_000),
            pair(gene="BBB", regulated=True, dist=40_000),
            pair(gene="CCC", regulated=False),
            pair(gene="AAA", cell="GM12878", regulated=True, dist=90_000),
            pair(chrom="chr2", gene="DDD", regulated=False),
        ]
    )
    assert [(r["cell"], r["targets"]) for r in rows] == [("GM12878", ["AAA"]), ("K562", ["AAA", "BBB"])]
    # the published distance is the CLOSEST of the element's targets, so a second far target can
    # never make a locus look further from its target than it is
    assert [r["distance"] for r in rows] == [90_000, 40_000]
    # an element with no regulated gene is not a locus
    assert all(r["chrom"] == "chr1" for r in rows)


# ------------------------------------------------------- steps 4 to 6: the rule, one element at a time
def test_the_shortcut_being_right_is_what_rejects_an_element():
    v = fourth.assess(element(targets=("AAA",)), ("AAA", 50_000), {"AAA", "ZZZ"}, [])
    assert v["locus"] is None
    assert v["reason"].startswith("the nearest coding TSS IS")


def test_an_element_whose_target_is_not_the_nearest_coding_tss_is_kept_with_its_trap():
    v = fourth.assess(element(targets=("AAA",)), ("ZZZ", 12_000), {"AAA", "ZZZ"}, [])
    assert v["reason"] is None
    assert v["locus"] == "AAA_K562_chr1_1000k"
    assert (v["nearest_coding"], v["nearest_coding_distance"]) == ("ZZZ", 12_000)


def test_every_target_must_clear_the_shortcut_not_just_one_of_them():
    """An element regulating two genes, one of them the nearest, is an element where the cheap answer
    works. Keeping it because its OTHER target is far would put the shortcut's own successes in the
    frame built to defeat it."""
    v = fourth.assess(element(targets=("AAA", "BBB")), ("BBB", 3_000), {"AAA", "BBB"}, [])
    assert v["reason"].startswith("the nearest coding TSS IS")


def test_a_non_coding_target_is_dropped_and_counted_apart():
    """It tests the H19 `predicted_coding`-first defect, not the geometry question, and mixing the two
    would confound both."""
    v = fourth.assess(element(targets=("LINC1",)), ("ZZZ", 9_000), {"ZZZ"}, [])
    assert v["reason"] == "a regulated target is not protein coding in GENCODE"
    assert v["non_coding_targets"] == ["LINC1"]


def test_the_non_coding_check_runs_before_the_geometry_check():
    """Otherwise a non-coding target passes the geometry rule for free - it can never be the nearest
    CODING TSS - and the frame fills up with the wrong defect."""
    v = fourth.assess(element(targets=("LINC1",)), ("ZZZ", 9_000), set(), [])
    assert v["reason"] == "a regulated target is not protein coding in GENCODE"


def test_an_element_beside_a_locus_of_an_earlier_frame_is_rejected():
    row = element()
    v = fourth.assess(row, ("ZZZ", 9_000), {"AAA", "ZZZ"}, [("chr1", row["mid"] + 50_000)])
    assert v["reason"] == f"within {fourth.KEEP_OUT} bp of a locus in an earlier frame"
    far = fourth.assess(row, ("ZZZ", 9_000), {"AAA", "ZZZ"}, [("chr1", row["mid"] + 250_000)])
    assert far["reason"] is None
    other = fourth.assess(row, ("ZZZ", 9_000), {"AAA", "ZZZ"}, [("chr2", row["mid"])])
    assert other["reason"] is None


def test_the_earlier_frames_are_all_three_of_them():
    names = {e.locus for e in fourth.earlier_frames()}
    for one in ("SHH_ZRS", "MYC_BENC", "TERT_promoter"):
        assert one in names
    assert len(fourth.earlier_frames()) == 17 + 11 + 9


# ---------------------------------------------------------------------------- step 7: order and cap
def test_genome_order_is_chromosome_then_start_and_not_lexicographic():
    rows = [element(chrom=c, start=s) for c, s in (("chr2", 5), ("chr10", 1), ("chr1", 9), ("chr1", 2))]
    assert [(r["chrom"], r["start"]) for r in sorted(rows, key=fourth.genome_key)] == [
        ("chr1", 2),
        ("chr1", 9),
        ("chr2", 5),
        ("chr10", 1),
    ]


# ------------------------------------------------------------------------- the budget, and the control
def test_the_budget_drops_from_the_end_of_genome_order_and_never_the_control():
    rows = [
        {"locus": "a", "requests": 1, "graded": True, "control": False},
        {"locus": "b", "requests": 1, "graded": True, "control": False},
        {"locus": "CONTROL_c", "requests": 1, "graded": True, "control": True},
        {"locus": "d", "requests": 0, "graded": True, "control": False},
    ]
    out = fourth.within_budget([dict(r) for r in rows], budget=2)
    assert [r["locus"] for r in out if r["graded"]] == ["a", "CONTROL_c", "d"]
    assert sum(r["requests"] for r in out) == 2
    assert out[1]["dropped_for_budget"] is True
    # the control's own request is inside the budget, so a budget of one buys the control and
    # nothing else. It is never the thing dropped, because without it the run is uninterpretable.
    tight = fourth.within_budget([dict(r) for r in rows], budget=1)
    assert [r["locus"] for r in tight if r["graded"]] == ["CONTROL_c", "d"]


def test_the_budget_leaves_a_run_that_already_fits_alone():
    rows = [{"locus": "a", "requests": 1, "graded": True, "control": False}]
    assert fourth.within_budget([dict(r) for r in rows], budget=5)[0]["graded"] is True


def test_the_control_is_excluded_from_the_reach_count_and_from_the_baselines():
    """It is the easiest case the data can offer, so scoring it into a rate would lift every rate."""
    rows = [
        {"locus": "a", "control": False, "askable": False, "targets_out_of_reach": [{"target": "X"}]},
        {"locus": "CONTROL_b", "control": True, "askable": True, "targets_out_of_reach": []},
    ]
    reach = fourth.reach_fatalities(rows)
    assert (reach["loci"], reach["died_at_the_reach_filter"]) == (1, 1)

    result = {
        "loci": [
            {"locus": "a", "expected": {"targets": ["AAA"]}, "readings": {"node": {"target": "AAA"}}},
            {"locus": "CONTROL_b", "expected": {"targets": ["BBB"]}, "readings": {"node": {"target": "BBB"}}},
        ]
    }
    out = fourth.baselines(result)
    assert (out["nearest_tss_in_node"]["k"], out["nearest_tss_in_node"]["n"]) == (1, 1)


def test_the_control_is_kept_out_of_the_aggregate_the_build_computed_with_it_in():
    """`loci.build` aggregates every row it is handed and this frame hands it the control, so the
    frame's own rate has to be recomputed. Reporting build's figure would lift every rate by the one
    locus chosen for being the easiest case in the data."""
    result = {
        "loci": [
            {
                "locus": "a",
                "expected": {
                    "targets": ["AAA"],
                    "direction": "activates",
                    "cells": [],
                    "gtex_tissues": [],
                    "nearest_gene_trap": "ZZZ",
                },
                "coding_genes_in_window": 10,
                "score": {
                    "target_hit_derived": False,
                    "target_hit_heuristic": False,
                    "target_hit_looked_up": False,
                    "cell_hit_derived": False,
                    "direction_hit_derived": False,
                    "class_hit_derived": False,
                    "heuristic_fell_in_trap": True,
                    "reachable_by_a_derived_target_layer": True,
                    "deletion_unaskable": None,
                    "scored": {"direction": {"judged": False}},
                },
                "readings": {},
            },
            {
                "locus": "CONTROL_b",
                "expected": {
                    "targets": ["BBB"],
                    "direction": "activates",
                    "cells": [],
                    "gtex_tissues": [],
                    "nearest_gene_trap": None,
                },
                "coding_genes_in_window": 10,
                "score": {
                    "target_hit_derived": True,
                    "target_hit_heuristic": True,
                    "target_hit_looked_up": False,
                    "cell_hit_derived": False,
                    "direction_hit_derived": False,
                    "class_hit_derived": False,
                    "heuristic_fell_in_trap": False,
                    "reachable_by_a_derived_target_layer": True,
                    "deletion_unaskable": None,
                    "scored": {"direction": {"judged": False}},
                },
                "readings": {},
            },
        ]
    }
    assert loci.aggregate(result["loci"])["target_derived"] == {
        "k": 1,
        "n": 2,
        "rate": 0.5,
        "hits": ["CONTROL_b"],
        "misses": ["a"],
    }
    without = fourth.without_the_control(result)
    assert (without["target_derived"]["k"], without["target_derived"]["n"]) == (0, 1)
    assert without["target_derived"]["hits"] == []


def test_what_the_layers_named_separates_the_target_from_the_trap():
    """The count between `named the target` and `named the nearest coding TSS` is the whole finding
    this frame exists to produce; a layer that names neither must not fall into either bucket."""

    def row(locus, deletion, trap, targets):
        return {
            "locus": locus,
            "expected": {"targets": targets, "nearest_gene_trap": trap},
            "readings": {"deletion": {"target": deletion}, "node": {"target": None}},
            "score": {
                "scored": {
                    "target": {
                        "by_layer": {"deletion": {"provenance": "derived", "hit": deletion in targets}}
                    }
                }
            },
        }

    result = {
        "loci": [
            row("a", "AAA", "ZZZ", ["AAA"]),
            row("b", "ZZZ", "ZZZ", ["AAA"]),
            row("c", "QQQ", "ZZZ", ["AAA"]),
            row("d", None, "ZZZ", ["AAA"]),
            row("CONTROL_e", "EEE", "ZZZ", ["EEE"]),
        ]
    }
    out = fourth.what_the_layers_named(result)
    assert out["what_it_named_instead"]["deletion"] == {
        "the published target": 1,
        "the nearest coding TSS": 1,
        "another gene": 1,
        "nothing": 1,
    }
    assert out["by_layer"]["deletion"] == {"provenance": "derived", "hit": 1, "n": 4}


def test_the_two_nearest_gene_rules_are_reported_apart():
    """The rule that DREW the frame is 0 by construction; the node rule is a real reading. Reporting
    one as the other is how a constructed frame turns into a beaten baseline."""
    result = {
        "loci": [
            {
                "locus": "a",
                "expected": {"targets": ["AAA"]},
                "readings": {"node": {"target": "AAA", "nearest_coding_anywhere": "ZZZ"}},
            },
            {
                "locus": "b",
                "expected": {"targets": ["BBB"]},
                "readings": {"node": {"target": None, "nearest_coding_anywhere": "YYY"}},
            },
        ]
    }
    out = fourth.baselines(result)
    assert out["nearest_coding_anywhere"]["k"] == 0
    assert (out["nearest_tss_in_node"]["k"], out["nearest_tss_in_node"]["n"]) == (1, 2)
    assert out["silences"] == 1


def test_the_direction_axis_reports_an_empty_repressor_arm_rather_than_pooling():
    result = {
        "loci": [
            {
                "locus": "a",
                "expected": {"direction": "activates"},
                "score": {"scored": {"direction": {"judged": True, "hit_derived": True}}},
            }
        ]
    }
    out = fourth.direction_readings(result)
    assert out["activators"] == {"k": 1, "n": 1}
    assert out["repressors"] == {"k": 0, "n": 0}


# ------------------------------------------------------------------------------- the shared globals
def test_the_keep_out_restores_the_shared_panel_even_when_the_draw_raises():
    """`loci.PANEL` is module-level state several sessions read. A draw that leaves it extended would
    silently change every other panel's negative controls."""
    original = loci.PANEL
    with pytest.raises(RuntimeError), fourth.keep_out_all_four_sets(()) as extended:
        assert len(extended) > len(original)
        raise RuntimeError("the draw failed")
    assert loci.PANEL is original


def test_register_stated_intervals_is_idempotent_and_only_ever_appends():
    before = loci.STATED_INTERVAL_RESULTS
    try:
        fourth.register_stated_intervals()
        fourth.register_stated_intervals()
        assert loci.STATED_INTERVAL_RESULTS.count(fourth.INTERVALS) == 1
        assert loci.STATED_INTERVAL_RESULTS[: len(before)] == before
    finally:
        loci.STATED_INTERVAL_RESULTS = before


def test_the_cell_counterpart_claims_a_reader_cell_only_where_there_is_one():
    """Handing HCT116 a colon GTEx tissue would win the cell axis a hit the reader never earned."""
    for cell, counterpart in fourth.CELL_COUNTERPART.items():
        for c in counterpart["cells"]:
            assert c in loci.READER_CELLS
        if cell in ("HCT116", "Jurkat", "WTC11"):
            assert counterpart["cells"] == () and counterpart["gtex_tissues"] == ()
