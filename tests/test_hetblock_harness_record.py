# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-hetblock step 2's record, guarded by the test that ships WITH it and never before it.

What is guarded here: that data/results/hetblock_harness.json says lane-scalehet's committed figures
came back BIT-IDENTICALLY -- not within a tolerance, equal -- through scalehet's own functions at
scalehet's own count and seed, because the whole lane is conditional on that; that the address-level
block structure it commits depends on no measured value and selects the governing unit by the
registered RULE and not by preference; that the singleton-block null reproduces scalehet's null draw
for draw, which is what shows the exchangeable unit is the only change; that every draw of both block
nulls holds every chromosome's gene count exactly; that the governing unit is not degenerate; and that
NO Q, NO R AND NO p OF ANY BLOCK NULL OF THIS LANE appears in this record, since it is committed
BEFORE the registration.

This file reads data/results/hetblock_harness.json and is committed in the same commit as it
(`22e12de`: a test that reads a record never ships before the record). There is no pytest.skip keyed
on that record either (`4f2652c`): the record is tracked at this commit, so the test must simply pass.
"""

from __future__ import annotations

import json
from pathlib import Path

import tracked_paths as tp

from genomeos.attribution import hetblock as hb
from genomeos.attribution import scalehet as sh

ROOT = Path(__file__).resolve().parents[1]
RECORD = ROOT / "data/results/hetblock_harness.json"
SCALEHET_RESULT = ROOT / "data/results/scalehet_result.json"


def record() -> dict:
    tp.must_be_present(RECORD)
    return json.loads(RECORD.read_text())


def test_the_record_is_tracked_and_is_this_lanes_step_two() -> None:
    d = record()
    assert d["result"] == "hetblock_harness"
    assert d["lane"] == hb.LANE
    assert d["step"].startswith("2")
    assert d["outcome"] == "HARNESS_REPRODUCES"
    assert d["problems"] == []


def test_every_one_of_scalehets_committed_figures_came_back_BIT_IDENTICALLY() -> None:
    rep = record()["scalehet_reproduction"]
    assert rep["all_bit_identical"] is True
    assert set(rep["bit_identical"]) >= {
        "q_observed",
        "r_observed",
        "global_scale",
        "p_q",
        "p_r",
        "null_q_median",
        "null_q_98_33",
        "null_r_median",
        "null_r_98_33",
        "admitted_count",
        "genes_admitted",
        "pairs_admitted",
    }
    assert all(rep["bit_identical"].values())
    # equal, not close: the recomputed block IS the committed block
    assert rep["recomputed"] == rep["committed"]


def test_the_reproduced_figures_are_the_ones_scalehets_own_result_carries() -> None:
    rep = record()["scalehet_reproduction"]
    tp.must_be_present(SCALEHET_RESULT)
    committed = hb.SCALEHET_COMMITTED
    assert rep["recomputed"]["q_observed"] == committed["q_observed"]
    assert rep["recomputed"]["r_observed"] == committed["r_observed"]
    assert rep["recomputed"]["p_q"] == committed["p_q"]
    assert rep["recomputed"]["p_r"] == committed["p_r"]
    assert rep["recomputed"]["admitted_count"] == committed["admitted_count"]
    assert rep["recomputed"]["genes_admitted"] == committed["genes_admitted"]


def test_the_count_that_picks_the_unit_is_address_level_and_the_graph_has_no_edges() -> None:
    g = record()["element_sharing_graph"]
    assert g["elements_driving_two_or_more_genes_of_this_table"] == 0
    assert g["every_component_is_one_gene"] is True
    assert g["largest_component"] == 1
    assert g["components_over_this_tables_genes"] == g["single_gene_components"]
    assert g["min_shared_elements"] == hb.MIN_SHARED_ELEMENTS
    assert "no measured value" in g["note"].lower()
    assert "no rna is read" in g["note"].lower()


def test_the_governing_unit_is_fixed_by_the_registered_rule_and_is_the_1mb_locus_block() -> None:
    u = record()["governing_unit"]
    assert u["rule"] == hb.BLOCK_SELECTION_RULE
    assert u["min_shared_elements"] == hb.MIN_SHARED_ELEMENTS
    assert u["elements_driving_two_or_more_genes_of_this_table"] < hb.MIN_SHARED_ELEMENTS
    assert u["governing_unit"] == "locus_block_1mb"


def test_the_block_structure_is_the_registered_spans_and_the_gene_join_is_complete() -> None:
    s = record()["locus_block_structure"]
    assert s["genes_without_a_start"] == 0
    assert s["genes_on_a_different_chromosome"] == 0
    assert s["this_lanes_starts_equal_freshrnas"] is True
    assert s["gene_table"] == hb.GENE_TABLE
    assert s["one_mb"]["span"] == hb.BLOCK_SPAN
    assert s["five_mb"]["span"] == hb.ROBUST_SPAN
    # strictly coarser: fewer and larger blocks over the same genes
    assert s["five_mb"]["blocks"] < s["one_mb"]["blocks"]
    assert s["five_mb"]["largest"] >= s["one_mb"]["largest"]
    assert s["five_mb"]["genes"] == s["one_mb"]["genes"]


def test_one_gene_per_block_reproduces_scalehets_null_DRAW_FOR_DRAW() -> None:
    e = record()["singleton_block_equivalence"]
    assert e["bit_identical"] is True
    assert e["q_draws_equal"] == e["draws"] == sh.PERMUTATIONS["count"]
    assert e["r_draws_equal"] == e["draws"]
    assert e["split_blocks_per_draw"]["max"] == 0


def test_every_draw_of_both_block_nulls_holds_every_chromosomes_gene_count_exactly() -> None:
    c = record()["exact_gene_counts"]
    assert len(c["observed_sizes"]) == hb.SCALEHET_COMMITTED["admitted_count"]
    assert sum(c["observed_sizes"]) == hb.SCALEHET_COMMITTED["genes_admitted"]
    for span in ("one_mb", "five_mb"):
        got = c["checked"][span]
        assert got["draws"] == sh.PERMUTATIONS["count"]
        assert got["draws_holding_every_count_exactly"] == got["draws"]
        assert got["genes"] == hb.SCALEHET_COMMITTED["genes_admitted"]


def test_the_governing_unit_is_not_degenerate_by_the_registered_rule() -> None:
    d = record()["degeneracy_check"]
    assert d["rule"] == hb.DEGENERACY_RULE
    assert d["degenerate"] is None
    assert d["admitted_chromosomes"] >= hb.MIN_CHROMOSOMES
    assert d["blocks_of_the_governing_unit"] >= hb.MIN_BLOCKS
    assert d["largest_block_share"] <= hb.MAX_BLOCK_SHARE


def test_NO_BLOCK_STATISTIC_OF_THIS_LANE_IS_IN_THIS_RECORD() -> None:
    """The harness is committed before the registration, so it may hold no Q, R or p of a block null."""
    d = record()
    text = json.dumps(
        {
            k: v
            for k, v in d.items()
            if k not in ("question", "why_not_already_answered", "control", "harness")
        }
    )
    assert "p_q_block" not in text and "q_block" not in text
    for key in ("block_null_q", "block_null_r", "p_q_1mb", "p_q_5mb", "outcome_block"):
        assert key not in text
    assert d["nothing_of_this_lanes_statistic_was_read"]
    assert d["alphagenome_requests"] == 0


def test_the_record_declares_its_readers_and_its_inputs_relatively() -> None:
    d = record()
    assert d["readers"], "every reader is declared"
    inputs = d["result_manifest"]["inputs"]
    assert inputs
    for entry in inputs:
        path = entry["path"]
        assert not path.startswith("/"), f"input_entry paths are RELATIVE (d72f959): {path}"
