# SPDX-License-Identifier: AGPL-3.0-or-later
"""The direction measurement's rules, above all the one that stops a one-sign subset reading as a pass."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

from genomeos.attribution import crispri_direction as cd

COLUMNS = [
    "chrom",
    "chromStart",
    "chromEnd",
    "EffectSize",
    "startTSS",
    "measuredGeneSymbol",
    "Significant",
    "ValidConnection",
    "CellType",
    "Regulated",
    "Dataset",
    "distanceToTSS",
]


def _row(**kw):
    base = {
        "chrom": "chr1",
        "chromStart": "1000",
        "chromEnd": "1500",
        "EffectSize": "-0.5",
        "startTSS": "50000",
        "measuredGeneSymbol": "GENEA",
        "Significant": "TRUE",
        "ValidConnection": "TRUE",
        "CellType": "K562",
        "Regulated": "TRUE",
        "Dataset": "Test2025",
        "distanceToTSS": "48750",
    }
    return {**base, **{k: str(v) for k, v in kw.items()}}


def _write_table(path: Path, rows):
    with gzip.open(path, "wt") as f:
        f.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(r[c] for c in COLUMNS) + "\n")


def _element(eid, start, end, coding=None, predicted=None, coding_by=None, by=None):
    e = {"id": eid, "start": start, "end": end}
    if coding:
        e["predicted_coding"] = {"gene": coding}
        e["predicted_coding_by_cell"] = coding_by or {}
    if predicted:
        e["predicted"] = {"gene": predicted}
        e["predicted_by_cell"] = by or {}
    return e


def test_wilson_is_finite_at_a_perfect_rate():
    """The normal interval collapses to zero width at k == n; Wilson must not."""
    w = cd.wilson(44, 44)
    assert w["rate"] == 1.0
    assert w["ci95"][0] < 1.0
    assert w["ci95"][1] == 1.0
    assert cd.wilson(0, 0)["rate"] is None


def test_signed_keeps_upward_significant_pairs_that_regulated_would_drop():
    table = [
        _row(EffectSize=-0.5, Significant="TRUE", Regulated="TRUE"),
        _row(EffectSize=0.4, Significant="TRUE", Regulated="FALSE"),
        _row(EffectSize=-0.2, Significant="FALSE", Regulated="FALSE"),
        _row(EffectSize=0.0, Significant="TRUE", Regulated="FALSE"),
    ]
    kept = cd.signed(table)
    assert len(kept) == 2
    assert sorted(float(r["EffectSize"]) for r in kept) == [-0.5, 0.4]


def test_in_reach_is_the_scorers_half_window():
    assert cd.in_reach(_row(distanceToTSS=cd.HALF_WINDOW))
    assert not cd.in_reach(_row(distanceToTSS=cd.HALF_WINDOW + 1))
    assert cd.in_reach(_row(distanceToTSS=-cd.HALF_WINDOW))


def test_model_value_prefers_the_strongest_statement_and_flags_disagreement():
    els = [
        _element("A", 900, 1600, coding="GENEA", coding_by={"K562": -0.1}),
        _element("B", 950, 1550, predicted="GENEA", by={"K562": 0.8}),
    ]
    m = cd.model_value(els, "GENEA", "K562")
    assert m["value"] == 0.8
    assert m["matches"] == 2
    assert m["matches_disagree"] is True
    assert cd.model_value(els, "GENEA", "HepG2") is None
    assert cd.model_value(els, "GENEB", "K562") is None


def test_agreement_counts_signs_and_excludes_a_zero_prediction():
    answered = [
        {"cell": "K562", "measured_sign": "down", "predicted_sign": "down"},
        {"cell": "K562", "measured_sign": "down", "predicted_sign": "up"},
        {"cell": "K562", "measured_sign": "up", "predicted_sign": "up"},
        {"cell": "K562", "measured_sign": "down", "predicted_sign": "zero"},
    ]
    a = cd.agreement(answered)
    assert (a["k"], a["n"]) == (2, 3)
    assert a["predicted_zero_excluded"] == 1
    assert a["both_measured_signs_present"] is True


def test_one_sign_subset_is_undecidable_however_high_the_agreement():
    """The rule this module exists for: a constant-down model on a measured-down subset is not a pass."""
    answered = [{"cell": "K562", "measured_sign": "down", "predicted_sign": "down"} for _ in range(44)]
    head = cd.agreement(answered)
    assert head["rate"] == 1.0
    assert head["both_measured_signs_present"] is False
    v = cd.judge(head, marginal=0.97)
    assert v["verdict"] == "undecidable"
    assert "constant-sign caller" in v["why"]


def test_judge_passes_only_when_both_margins_clear():
    head = {"rate": 0.9, "ci95": [0.8, 0.95], "both_measured_signs_present": True}
    assert cd.judge(head, marginal=0.55)["verdict"] == "passed"
    assert cd.judge(head, marginal=0.85)["verdict"] == "undecidable"  # margin of 0.10 not met
    assert cd.judge(head, marginal=0.82)["verdict"] == "undecidable"  # interval not clear of it
    assert (
        cd.judge({"rate": 0.6, "ci95": [0.4, 0.7], "both_measured_signs_present": True}, 0.1)["verdict"]
        == "failed"
    )
    assert (
        cd.judge({"rate": None, "ci95": None, "both_measured_signs_present": False}, 0.5)["verdict"]
        == "refused"
    )


def test_collect_walks_every_stratum(tmp_path):
    els = tmp_path / "elements"
    els.mkdir()
    (els / "chr1.json").write_text(
        json.dumps(
            [
                _element("A", 900, 1600, coding="GENEA", coding_by={"K562": -0.3}),
                _element("C", 5000, 5600, coding="OTHER", coding_by={"K562": -0.2}),
            ]
        )
    )
    table = [
        _row(),  # answerable
        _row(chromStart=5000, chromEnd=5600, measuredGeneSymbol="GENEA"),  # covered, not the top target
        _row(chromStart=90000, chromEnd=90500),  # element not in the sweep
        _row(distanceToTSS=cd.HALF_WINDOW + 10),  # out of reach
        _row(CellType="WTC11"),  # no model line
    ]
    got = cd.collect(table, ("K562", "GM12878"), elements=els)
    s = got["strata"]
    assert s["answerable::K562"] == 1
    assert s["covered_but_gene_not_the_top_target::K562"] == 1
    assert s["element_not_in_the_sweep::K562"] == 1
    assert s["out_of_reach::K562"] == 1
    assert s["refused_no_model_line::WTC11"] == 1
    assert got["answered"][0]["predicted_sign"] == "down"
    assert got["not_top_target"][0]["model_named_instead"] == ["OTHER"]


def test_marginal_down_rate_counts_the_whole_sweep(tmp_path):
    els = tmp_path / "elements"
    els.mkdir()
    (els / "chr1.json").write_text(
        json.dumps(
            [
                _element("A", 1, 2, coding="G", coding_by={"K562": -0.3}),
                _element("B", 3, 4, coding="G", coding_by={"K562": -0.1}),
                _element("C", 5, 6, coding="G", coding_by={"K562": 0.2}),
                _element("D", 7, 8, predicted="G", by={"K562": -0.4}),
                _element("E", 9, 10, coding="G", coding_by={}),
            ]
        )
    )
    m = cd.marginal_down_rate(("K562",), elements=els)
    assert m["K562"]["k"] == 3
    assert m["K562"]["n"] == 4
    assert m["K562"]["rate"] == 0.75
    assert m["K562"]["no_value"] == 1
    assert m["elements_in_the_sweep"] == 5


def test_table_round_trip(tmp_path):
    p = tmp_path / "t.tsv.gz"
    _write_table(p, [_row(), _row(ValidConnection="overlaps target gene exon")])
    got = cd.rows("t.tsv.gz", knowledge=tmp_path)
    assert len(got) == 1
