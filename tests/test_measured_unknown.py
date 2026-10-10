# SPDX-License-Identifier: AGPL-3.0-or-later
"""M5: the census of public measurement over the real-unknown blocks (attribution/measured_unknown.py)."""

from __future__ import annotations

import json
from dataclasses import fields
from pathlib import Path

import pytest
import tracked_paths as tp

from genomeos.attribution import measured_unknown as mu

RESULT = Path("data/results/measured_unknown_inventory.json")
ELEMENTS_HOW = "the AlphaGenome all-element archive (data/knowledge/alphagenome/all_elements), from the sweep"


def block(chrom="chr1", start=1000, end=2000, elements=((1100, 1300),), case="relaxed", genes=(), coding=()):
    return {
        "chrom": chrom,
        "start": start,
        "end": end,
        "length": end - start,
        "case": case,
        "elements": list(elements),
        "named_genes": list(genes),
        "named_cells": [],
        "named_coding_genes": list(coding),
        "named_coding_cells": [],
    }


# ---- parsers: coordinates and the gene, never a value ---------------------------------------------------


def test_a_tested_row_has_no_field_for_a_readout():
    assert [f.name for f in fields(mu.Tested)] == ["chrom", "start", "end", "gene", "cell"]


def test_region_reads_coordinates_and_ignores_a_strand():
    assert mu.region("chr1:100-200") == ("chr1", 100, 200)
    assert mu.region("chr10:80712806-80713306:.") == ("chr10", 80712806, 80713306)
    assert mu.region("*") is None
    assert mu.region("negative_control") is None


def test_bed_rows_skip_headers_tracks_and_malformed_rows():
    text = [
        "track name=x\n",
        "#c\n",
        "chrom\tstart\tend\n",
        "chr2\t5\t9\t42.0\n",
        "chr2\t9\t5\n",
        "x\t1\t2\n",
    ]
    assert list(mu.bed_rows(text, cell="K562")) == [mu.Tested("chr2", 5, 9, None, "K562")]


def test_table_rows_select_named_columns_and_drop_the_rest():
    text = [
        "chrPerturbationTarget\tstartPerturbationTarget\tendPerturbationTarget\tEffectSize\tmeasuredGeneSymbol\n",
        "chr3\t10\t20\t-0.9\tGATA1\n",
        "chr3\t30\t40\t0.1\tNA\n",
        "NTC\tNTC\tNTC\t0\tMKX\n",
    ]
    cols = ("chrPerturbationTarget", "startPerturbationTarget", "endPerturbationTarget")
    got = list(mu.table_rows(text, *cols, gene="measuredGeneSymbol"))
    assert got == [mu.Tested("chr3", 10, 20, "GATA1"), mu.Tested("chr3", 30, 40, None)]
    with pytest.raises(KeyError):
        list(mu.table_rows(text, "chr", "start", "end"))


def test_table_or_bed_rows_reads_a_headerless_table_by_its_first_three_columns():
    headerless = ["chr8\t127367301\t127367304\tMYC|chr8:x\t676\t+\n"]
    assert mu.table_or_bed_rows(headerless, "chrPerturbationTarget", "a", "b") == [
        mu.Tested("chr8", 127367301, 127367304)
    ]
    with pytest.raises(KeyError):
        mu.table_or_bed_rows(["name\ttype\n", "a\tb\n"], "chr", "start", "stop")


def test_region_column_and_fasta_name_rows():
    assert list(mu.region_column_rows(["OligoID\ttarget\n", "o1\tchrX:5-9\n", "o2\t*\n"], "target")) == [
        mu.Tested("chrX", 5, 9)
    ]
    fa = [">chr1_540951_540970_+_chr1.1.1_GGAG\n", "GGAGAGCTG\n", ">control_AATG_NTC\n", "AATG\n"]
    assert list(mu.fasta_name_rows(fa)) == [(mu.Tested("chr1", 540951, 540971), "+", "GGAGAGCTG")]


# ---- the lift ------------------------------------------------------------------------------------------


CHAIN = """chain 1000 chrA 1000 + 100 400 chrB 2000 + 500 800 1
100 50 50
150
chain 900 chrA 1000 + 600 700 chrC 500 - 0 100 2
100

"""


def test_chain_lifts_inside_a_block_and_refuses_a_gap_or_a_short_match():
    c = mu.Chain.from_lines(CHAIN.splitlines())
    assert c.lift("chrA", 110, 120) == ("chrB", 510, 520)
    # the second block of chain 1 starts at t=250, q=650
    assert c.lift("chrA", 260, 270) == ("chrB", 660, 670)
    assert c.lift("chrA", 210, 220) is None  # wholly in the gap
    assert c.lift("chrA", 190, 260) is None  # 20 of 70 bases mapped: under minMatch 0.95
    assert c.lift("chrZ", 1, 2) is None


def test_chain_lifts_onto_the_minus_strand():
    c = mu.Chain.from_lines(CHAIN.splitlines())
    # t 600..700 maps to q- 0..100 of a 500 bp query, which is forward 400..500
    assert c.lift("chrA", 600, 610) == ("chrC", 490, 500)
    lifted, lost = mu.lift_all(c, [mu.Tested("chrA", 600, 610, "G", "K562"), mu.Tested("chrA", 210, 220)])
    assert lifted == [mu.Tested("chrC", 490, 500, "G", "K562")] and lost == 1


# ---- the population and the census -----------------------------------------------------------------------


def test_population_of_keeps_real_unknown_and_assigns_elements_by_midpoint():
    blocks = [
        {"start": 0, "end": 100, "tier": "constrained_unknown", "copy": False, "case": "syntax"},
        {"start": 200, "end": 300, "tier": "constrained_unknown", "copy": True, "case": "syntax"},
        {"start": 400, "end": 500, "tier": "neutral", "copy": False, "case": None},
    ]
    elements = [
        {"start": 90, "end": 130, "predicted": {"gene": "A", "tissue": "liver", "log2_fold_change": 0.5}},
        {"start": 10, "end": 20, "predicted": {"gene": "B", "tissue": "lung", "log2_fold_change": 0.05}},
        {
            "start": 40,
            "end": 60,
            "predicted": {"gene": "C", "log2_fold_change": -0.2},
            "predicted_coding": {"gene": "C", "tissue": "gut"},
        },
    ]
    (row,) = mu.population_of("chr9", blocks, elements)
    assert (row["start"], row["end"], row["case"]) == (0, 100, "syntax")
    assert row["elements"] == [(10, 20), (40, 60)]  # 90-130 has its midpoint at 110, outside
    assert row["named_genes"] == ["C"]  # B moves under the 0.1 bar
    assert row["named_coding_genes"] == ["C"] and row["named_coding_cells"] == ["gut"]
    assert mu.gate([row])["ok"] is False


def test_census_counts_touch_inside_reciprocal_match_and_genes():
    rows = [
        block(start=1000, end=2000, elements=[(1100, 1300)], coding=["GATA1"]),
        block(start=5000, end=6000, elements=[], case="tolerant"),
    ]
    tested = [
        mu.Tested("chr1", 1120, 1290, "GATA1"),  # inside, matches the element both ways, attributed gene
        mu.Tested("chr1", 1950, 2100, "LINC1"),  # touches, midpoint outside, no match
        mu.Tested("chr1", 5500, 5600, "MYC"),
        mu.Tested("chr1", 3000, 3100, "MYC"),  # between blocks
        mu.Tested("chr2", 1100, 1300, "GATA1"),  # another chromosome
    ]
    c = mu.census(rows, tested, coding={"GATA1", "MYC"})
    assert c["tested_intervals"] == 5
    assert c["blocks_touched"] == 2 and c["blocks_touched_carrying"] == 1
    assert c["blocks_element_matched"] == 1 and c["elements_matched"] == 1
    assert c["tested_inside"] == 2
    assert c["blocks_with_tested_inside"] == {"at_least_1": 2, "at_least_2": 0, "at_least_4": 0}
    assert c["blocks_tested_against_a_coding_gene"] == 2
    assert c["blocks_tested_against_an_attributed_gene"] == 1
    assert c["blocks_touched_by_case"] == {"relaxed": 1, "tolerant": 1}


def test_census_without_genes_reports_no_gene_fields():
    c = mu.census([block()], [mu.Tested("chr1", 1500, 1501)], coding=None)
    assert c["blocks_touched"] == 1
    assert "blocks_tested_against_a_coding_gene" not in c


def test_shortfall_sizes_only_a_gene_level_endpoint():
    assert mu.shortfall(13, "gene") == {
        "gene_level_endpoint": True,
        "blocks_usable_for_the_clause": 13,
        "short_of_the_floor": 62,
        "short_of_a_half": 518,
        "designs_within_reach": [],
    }
    assert mu.shortfall(600, "activity")["blocks_usable_for_the_clause"] == 0
    assert mu.shortfall(531, "variant_gene")["designs_within_reach"] == ["half", "quarter", "tenth", "none"]
    assert min(v[0] for v in mu.SIZING.values() if v) == mu.FLOOR_BLOCKS


@pytest.mark.needs_local_data("data/knowledge/alphagenome/all_elements/chr21.json", how=ELEMENTS_HOW)
def test_population_of_chr21_is_the_organisers_real_unknown():
    from genomeos.attribution import organise

    want = [b for b in organise.blocks("chr21") if b["tier"] == "constrained_unknown" and not b["copy"]]
    got = mu.population(["chr21"])
    assert [(r["start"], r["end"]) for r in got] == [(b["start"], b["end"]) for b in want]


# ---- the committed result ---------------------------------------------------------------------------------


def test_committed_inventory_reads_no_value_and_its_shortfall_follows_its_counts():
    # the inventory is committed (7906142), so its absence is a broken checkout, not "not run here"
    tp.must_be_present(RESULT, was="the inventory result is not written yet")
    d = json.loads(RESULT.read_text())
    assert d["population"]["gate"]["ok"] is True
    assert d["readout_values_read"] == 0
    u = d["gene_level_union"]
    coding = u["gene_resolved_sources"]["blocks_tested_against_a_coding_gene"]
    assert u["shortfall"]["short_of_the_floor"] == max(0, mu.FLOOR_BLOCKS - coding)
    if coding < mu.FLOOR_BLOCKS:
        assert d["verdict"]["element_perturbation_suffices"] is False
        assert d["verdict"]["decision"]["registered_test"] is False
    upper = u["with_unresolved_libraries_as_upper_bound"]["blocks_touched"]
    assert upper >= u["gene_resolved_sources"]["blocks_touched"]
    keys = {r["source"] for r in d["table"]}
    for need in (
        "epcrispr_training_k562",
        "gasperini2019_at_scale",
        "encode_mpra_element_references",
        "vista",
        "gtex_v8_finemapped_eqtl",
    ):
        assert need in keys
    for row in d["table"]:
        if row["endpoint"] not in mu.GENE_ENDPOINTS:
            assert row["blocks_tested_against_a_coding_gene"] is None
