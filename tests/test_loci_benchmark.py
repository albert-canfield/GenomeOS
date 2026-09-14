# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The machinery against loci whose answer is already published (docs/LOCI-BENCHMARK.md).

The run needs local reference data and the network (scripts/loci_benchmark.py). These tests read
its committed result, so a locus that used to pass and stops passing fails CI within a second.
The expectation and scoring tests below need nothing but the code.
"""

from __future__ import annotations

import pytest

from genomeos.benchmark.loci import (
    CLASSES,
    PANEL,
    READER_CELLS,
    Expect,
    score_class,
    score_target,
)
from genomeos.results import load_result

RESULT = load_result("loci_benchmark")
needs_result = pytest.mark.skipif(not RESULT, reason="benchmark result not present")


# ------------------------------------------------------------------ the expectations, offline
def test_the_panel_is_small_and_covers_every_class():
    assert 10 <= len(PANEL) <= 14, "the panel must stay small enough to run and to read"
    seen = {c for e in PANEL for c in e.classes}
    assert seen == set(CLASSES)
    names = {e.locus for e in PANEL}
    assert {"HBB_LCR", "HOXD"} <= names, "ordered in time and ordered in space are both required"


def test_every_expectation_is_written_down_with_its_sources():
    for e in PANEL:
        assert e.window[0] <= e.element[0] < e.element[1] <= e.window[1], e.locus
        assert e.targets and e.citations and e.answer_from, e.locus
        assert set(e.classes) <= set(CLASSES), e.locus
        assert set(e.cells) <= set(READER_CELLS), e.locus
        assert e.direction in ("activates", "represses", "coding", "none"), e.locus
        for v in e.variants:
            if v.get("pos"):
                assert e.window[0] <= v["pos"] <= e.window[1], f"{e.locus} {v.get('rsid')}"


def _toy(classes=("program",)) -> Expect:
    return Expect(
        locus="toy",
        chrom="chr1",
        element=(100, 200),
        window=(0, 1000),
        classes=classes,
        targets=("GOOD",),
        direction="activates",
        tissues=("x",),
        citations=("none",),
        answer_from=("none",),
    )


def test_a_layer_hits_only_when_it_ranks_the_published_target_first():
    readings = {
        "deletion": {"provenance": "derived", "targets": [{"gene": "BAD"}, {"gene": "GOOD"}]},
        "node": {"provenance": "heuristic", "target": "GOOD"},
    }
    s = score_target(_toy(), readings)
    assert not s["hit_derived"] and s["hit_derived_among"]
    assert s["hit_heuristic"] and not s["hit_looked_up"]


def test_a_lookup_never_counts_as_a_derived_hit():
    readings = {"lookups": {"provenance": "looked_up", "gwas_genes": ["GOOD"], "clinvar_genes": []}}
    s = score_target(_toy(), readings)
    assert s["hit_looked_up"] and not s["hit_derived"]


def test_the_class_reading_does_not_look_at_the_expectation():
    readings = {
        "constraint": {"mammals": {"fraction_above": 0.5}, "humans": None},
        "syntax_values": {"values_in_syntax": 0},
        "deletion": {"target": "GOOD"},
    }
    a = score_class(_toy(("program",)), readings)
    b = score_class(_toy(("storage",)), readings)
    assert a["read"] == b["read"] == ["syntax", "program"]
    assert a["hit_derived"] and not b["hit_derived"]


# ------------------------------------------------------------------------- the gate, on the result
#: loci whose target a derived layer named on the first run (2026-09-14). A locus may join this set;
#: one that leaves it is a regression.
DERIVED_TARGET_PASSES = {"SHH_ZRS", "HBB_LCR", "BCL11A_enhancer", "ABO", "HLA_DRB1", "APP", "TP53", "HOXD"}
#: loci whose published causal variant constraint ranks first among the window's variable positions
CAUSAL_VARIANT_FIRST = {"HERC2_OCA2": "rs12913832", "BCL11A_enhancer": "rs1427407"}


def loci() -> dict:
    return {r["locus"]: r for r in RESULT["loci"]}


@needs_result
def test_no_locus_that_passed_stops_passing():
    now = {name for name, r in loci().items() if r["score"]["target_hit_derived"]}
    lost = DERIVED_TARGET_PASSES - now
    assert not lost, f"a derived layer stopped naming the published target at {sorted(lost)}"


@needs_result
def test_the_result_is_the_panel_in_the_code():
    assert [r["locus"] for r in RESULT["loci"]] == [e.locus for e in PANEL]


@needs_result
def test_no_hit_is_annotation_lookup_alone():
    assert RESULT["aggregate"]["target_only_looked_up"] == []


@needs_result
def test_derived_and_looked_up_rates_are_reported_apart():
    a = RESULT["aggregate"]
    for key in ("target_derived", "target_heuristic", "target_looked_up"):
        assert a[key]["n"] == len(PANEL) and "misses" in a[key]
    for r in RESULT["loci"]:
        for layer in r["score"]["scored"]["target"]["by_layer"].values():
            assert layer["provenance"] in ("derived", "targeted", "heuristic", "looked_up")


@needs_result
def test_negative_controls_are_scored_through_the_same_claims():
    n = RESULT["negative_controls"]
    assert n["windows"] >= 5 * len(PANEL)
    for claim in ("target", "cell", "direction", "storage"):
        assert n["positives"][claim]["n"] == len(PANEL)
        assert n["negatives"][claim]["n"] == n["windows"]


@needs_result
def test_a_value_on_syntax_still_separates_the_panel_from_its_controls():
    """The one claim the controls did not erase: it must stay at least twice as common at the panel."""
    n = RESULT["negative_controls"]
    assert n["positives"]["storage"]["rate"] >= 2 * n["negatives"]["storage"]["rate"]


@needs_result
def test_the_causal_variants_constraint_found_stay_found():
    for name, rsid in CAUSAL_VARIANT_FIRST.items():
        rank = loci()[name]["score"]["scored"]["variant"]["constraint_rank"][rsid]
        assert rank["rank"] == 1, f"{name}: {rsid} fell to {rank}"


@needs_result
def test_the_coding_calibration_derives_the_published_protein_changes():
    for name, rsid, change in (("APP", "rs63750847", "p.Ala673Thr"), ("TP53", "rs1042522", "p.Pro72Arg")):
        got = loci()[name]["score"]["scored"]["variant"]["protein_change"][rsid]
        assert got == change, f"{name}: the engine derived {got}"


@needs_result
def test_the_hypersensitive_ladder_of_the_globin_locus_is_found_from_dnase_alone():
    ladder = loci()["HBB_LCR"]["readings"]["ladder"]
    assert ladder["erythroid_specific"] >= 4
    assert 2_000 <= ladder["median_spacing"] <= 5_000


@needs_result
def test_every_unreachable_reading_says_why():
    for r in RESULT["loci"]:
        for layer, reading in r["readings"].items():
            if isinstance(reading, dict) and reading.get("pending") is not None:
                assert len(reading["pending"]) > 20, (
                    f"{r['locus']} {layer}: a pending reading without a reason"
                )
    order = loci()["HOXD"]["readings"]["order"]
    assert "order" in order["pending"], "the untestable order of HOXD must stay recorded as pending"


@needs_result
def test_the_value_domain_tells_the_many_valued_slot_from_the_two_valued_ones():
    read = {n: r["score"]["scored"]["class"]["value_domain_read"] for n, r in loci().items()}
    assert read["HLA_DRB1"] == "many"
    assert read["HERC2_OCA2"] == "two" and read["APP"] == "two"
    assert read["ABO"] == "few"


@needs_result
def test_the_named_common_values_come_back_with_their_population():
    expected = {
        "HERC2_OCA2": ("rs12913832", "European"),
        "MCM6_LCT": ("rs4988235", "European"),
        "FTO_IRX3": ("rs1421085", "European"),
        "TP53": ("rs1042522", "European"),
    }
    for name, (rsid, population) in expected.items():
        f = loci()[name]["readings"]["frequencies"]["named"][rsid]
        assert f["af"] >= 0.05 and population in (f["commonest_in"] or ""), f"{name}: {f}"


@needs_result
def test_the_target_rate_is_reported_against_its_chance_floor():
    a = RESULT["aggregate"]
    chance = a["target_by_chance"]
    assert chance["of"] == len(PANEL) and 0 < chance["expected"] < len(PANEL)
