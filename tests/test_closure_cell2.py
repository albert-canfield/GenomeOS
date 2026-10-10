# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The search for M3's successor: the gate, the arm counting and the record's integrity."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from genomeos.attribution import closure_cell as cc
from genomeos.attribution import closure_cell2 as cc2

RESULT = Path("data/results/closure_cell2_inventory.json")


def test_the_gate_reproduces_m3s_own_detectable_separation():
    """M3's registered power number, from this module's formula: 35 against a floor of 5,427."""
    assert cc2.smallest_detectable(35, 5427) == cc2.M3_DETECTABLE
    assert cc2.smallest_detectable(145, 5427) == 0.5604
    assert cc2.smallest_detectable(0, 5427) is None
    assert cc2.smallest_detectable(35, 0) is None


def test_the_gate_needs_both_size_and_power():
    assert cc2.qualifies(cc2.M3_ARM_A, cc2.M3_DETECTABLE) is True
    assert cc2.qualifies(cc2.M3_ARM_A - 1, cc2.M3_DETECTABLE) is False  # smaller than M3
    assert cc2.qualifies(400, 0.6218) is False  # bigger but worse powered
    assert cc2.qualifies(400, 0.55) is True
    assert cc2.qualifies(35, None) is False


def test_only_cells_with_a_track_and_a_dependency_column_are_reachable():
    """The chain needs the sweep's track and a DepMap column; the record says which cells have both."""
    assert set(cc2.CELLS) == {"K562", "HepG2"}
    assert set(cc2.CELLS_WITHOUT_A_DEPENDENCY_COLUMN) == {"GM12878", "IMR-90"}
    for candidate in cc2.CANDIDATES:
        if candidate["reachable"]:
            assert candidate["cell"] in cc2.CELLS, candidate["key"]
        else:
            assert candidate["cell"] not in cc2.CELLS, candidate["key"]
            assert candidate["not_reachable_because"], candidate["key"]


def test_every_candidate_names_a_library_a_screen_and_a_readout_it_does_not_fetch():
    keys = [c["key"] for c in cc2.CANDIDATES]
    assert len(keys) == len(set(keys))
    for candidate in cc2.CANDIDATES:
        assert candidate["library"] in cc2.LIBRARIES, candidate["key"]
        assert candidate["screen"] and candidate["assay"] and candidate["lab"], candidate["key"]
        readout = candidate["readout_files"]
        assert readout["screen"] and readout["bytes"] > 0, candidate["key"]
        named = {candidate["screen"], *(candidate.get("control") or "").split()}
        assert cc.SCREEN["experiment"] not in named, "M3's own screen is not a candidate"
        assert cc.CONTROL["experiment"] not in named, "M3's own control is not a candidate"


def test_effect_in_is_per_cell_and_signed_like_m3s():
    element = {"predicted_coding_by_cell": {"K562": -0.4, "HepG2": 0.25}}
    assert cc2.effect_in(element, "K562") == 0.4
    assert cc2.effect_in(element, "HepG2") == -0.25
    assert cc2.effect_in(element, "GM12878") is None
    assert cc2.effect_in({}, "K562") is None
    assert cc2.effect_in(element, "K562") == cc.element_effect(element)


def _rows():
    """One element per arm of M3's rule, plus a promoter, a weak one and one with no value."""
    base = {"chrom": "chr1", "nearest_coding_distance": 50_000}
    return [
        {
            **base,
            "id": "A",
            "start": 100,
            "end": 200,
            "gene": "ESS",
            "nearest_coding_gene": "NOT",
            "effects": {"K562": 0.5},
        },
        {
            **base,
            "id": "D",
            "start": 300,
            "end": 400,
            "gene": "NOT",
            "nearest_coding_gene": "NOT",
            "effects": {"K562": 0.5},
        },
        {
            **base,
            "id": "B",
            "start": 500,
            "end": 600,
            "gene": "ESS",
            "nearest_coding_gene": "ESS",
            "effects": {"K562": 0.5},
        },
        {
            **base,
            "id": "C",
            "start": 700,
            "end": 800,
            "gene": "NOT",
            "nearest_coding_gene": "ESS",
            "effects": {"K562": 0.5},
        },
        {
            **base,
            "id": "P",
            "start": 900,
            "end": 1000,
            "gene": "ESS",
            "nearest_coding_gene": "ESS",
            "effects": {"K562": 0.5},
            "nearest_coding_distance": 500,
        },
        {
            **base,
            "id": "W",
            "start": 1100,
            "end": 1200,
            "gene": "ESS",
            "nearest_coding_gene": "NOT",
            "effects": {"K562": 0.02},
        },
        {
            **base,
            "id": "U",
            "start": 1300,
            "end": 1400,
            "gene": "ESS",
            "nearest_coding_gene": "NOT",
            "effects": {"K562": None},
        },
    ]


def _design_covering(rows, per_element: int):
    """`per_element` single-base intervals inside every row's interval."""
    out: dict[str, list[tuple[int, int, None]]] = {}
    for row in rows:
        for i in range(per_element):
            out.setdefault(row["chrom"], []).append((row["start"] + i, row["start"] + i + 1, None))
    return out


def test_arms_over_divides_the_covered_elements_the_way_m3s_rule_does():
    rows = _rows()
    counted = cc2.arms_over(rows, _design_covering(rows, cc.MIN_GUIDES), {"ESS": -0.9, "NOT": 0.0}, "K562")
    assert counted["counts"]["elements_with_a_guide"] == len(rows)
    assert counted["counts"]["elements_with_min_guides"] == len(rows)
    assert counted["counts"]["no_value_on_the_cells_own_track"] == 1  # U
    assert counted["counts"]["promoter_elements"] == 1  # P
    assert counted["sizes"] == {
        "attributed_essential_neighbourhood_not": 1,
        "attributed_not_neighbourhood_not": 1,
        "attributed_essential_neighbourhood_essential": 1,
        "attributed_not_neighbourhood_essential": 1,
        cc.POSITIVE_CONTROL: 1,
    }
    assert counted["in_no_arm"] == {"prediction_under_the_magnitude_bar": 1}  # W
    assert counted["arm_a"] == 1
    assert counted["floor_arm"] == 1
    assert counted["qualifies"] is False


def test_an_element_under_min_guides_is_not_in_the_population():
    rows = _rows()
    counted = cc2.arms_over(
        rows, _design_covering(rows, cc.MIN_GUIDES - 1), {"ESS": -0.9, "NOT": 0.0}, "K562"
    )
    assert counted["counts"]["elements_with_a_guide"] == len(rows)
    assert counted["counts"]["elements_with_min_guides"] == 0
    assert counted["sizes"] == {}
    assert counted["smallest_detectable_separation"] is None


def test_arms_over_reads_the_cell_it_is_asked_about():
    rows = [
        {
            "chrom": "chr1",
            "id": "A",
            "start": 100,
            "end": 200,
            "gene": "ESS",
            "nearest_coding_gene": "NOT",
            "nearest_coding_distance": 50_000,
            "effects": {"K562": 0.5, "HepG2": 0.02},
        }
    ]
    design = _design_covering(rows, cc.MIN_GUIDES)
    effects = {"ESS": -0.9, "NOT": 0.0}
    assert cc2.arms_over(rows, design, effects, "K562")["arm_a"] == 1
    hepg2 = cc2.arms_over(rows, design, effects, "HepG2")
    assert hepg2["arm_a"] == 0
    assert hepg2["in_no_arm"] == {"prediction_under_the_magnitude_bar": 1}


def test_read_design_counts_targeting_rows_only(tmp_path):
    accession = "ENCFF135FSS"  # a csv with a guideType column
    path = tmp_path / f"{accession}.csv"
    path.write_text(
        "chrom,chromStart,chromEnd,guideType\n"
        "chr1,100,103,targeting\n"
        "chr1,200,203,negative_control\n"
        "chr2,300,303,targeting\n"
        "notachrom,1,2,targeting\n"
    )
    intervals, census = cc2.read_design(accession, store=tmp_path)
    assert census == {"rows": 4, "counted": 2, "not_targeting": 1, "malformed": 1}
    assert intervals == {"chr1": [(100, 103, None)], "chr2": [(300, 303, None)]}
    assert cc2.median_width(intervals) == 3


def test_read_design_of_a_bed_counts_every_placed_row(tmp_path):
    accession = "ENCFF526YZH"  # a bed3 of deletion intervals
    path = tmp_path / f"{accession}.bed.gz"
    with gzip.open(path, "wt") as fh:
        fh.write("track name=x\nchr1\t100\t1100\nchr1\t200\t1300\nbad\n")
    intervals, census = cc2.read_design(accession, store=tmp_path)
    assert census["counted"] == 2
    assert census["malformed"] == 1
    assert cc2.median_width(intervals) == 1050


def test_fetch_design_refuses_a_file_whose_md5_is_not_encodes(tmp_path):
    accession = "ENCFF232FVA"
    path = tmp_path / f"{accession}.bed.gz"
    with gzip.open(path, "wt") as fh:
        fh.write("chr1\t1\t4\n")
    with pytest.raises(ValueError, match="is not ENCODE's"):
        cc2.fetch_design(accession, store=tmp_path, net=False)


@pytest.mark.skipif(not RESULT.exists(), reason="the search has not been run on this machine")
def test_the_committed_search_reports_every_candidate_against_the_same_gate():
    record = json.loads(RESULT.read_text())
    assert record["gate_fixed_before_the_counting"]["arm_a_at_least"] == cc2.M3_ARM_A
    assert len(record["candidates"]) == len(cc2.CANDIDATES)
    assert record["alphagenome_requests"] == 0
    qualifying = [c["key"] for c in record["candidates"] if c.get("qualifies")]
    assert record["qualifying"] == qualifying
    assert record["outcome"] == ("A_CANDIDATE_QUALIFIES" if qualifying else "NO_CANDIDATE_QUALIFIES")
    for candidate in record["candidates"]:
        if not candidate["reachable"]:
            assert candidate["qualifies"] is False
            continue
        coverage = candidate["coverage"]
        assert candidate["qualifies"] == cc2.qualifies(
            coverage["arm_a"], coverage["smallest_detectable_separation"]
        )
        assert coverage["arm_a"] <= coverage["counts"]["elements_with_min_guides"]


def test_the_positive_controls_key_says_k562_and_the_record_says_whose_count_it_is():
    """M3's control arm carries K562 in its key; in HepG2 the count under it is HepG2's.

    The keys are M3's so a candidate's table reads against M3's result field for field, which
    leaves one key naming a cell it may not be counting. The sentence the record carries is what
    keeps that from being read as a K562 number.
    """
    assert cc.POSITIVE_CONTROL == "promoter_of_a_gene_k562_depends_on"
    assert "k562" in cc2.POSITIVE_CONTROL_READS.lower()
    assert "hepg2's own" in cc2.POSITIVE_CONTROL_READS.lower()
    rows = [
        {
            "chrom": "chr1",
            "id": "P",
            "start": 100,
            "end": 200,
            "gene": "ESS",
            "nearest_coding_gene": "ESS",
            "nearest_coding_distance": 500,
            "effects": {"K562": 0.5, "HepG2": 0.5},
        }
    ]
    design = _design_covering(rows, cc.MIN_GUIDES)
    in_hepg2 = cc2.arms_over(rows, design, {"ESS": -0.9}, "HepG2")
    assert in_hepg2["sizes"] == {cc.POSITIVE_CONTROL: 1}  # counted from HepG2's dependency
    assert cc2.arms_over(rows, design, {"ESS": 0.0}, "HepG2")["sizes"] == {}
