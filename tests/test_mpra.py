# SPDX-License-Identifier: AGPL-3.0-or-later
"""lentiMPRA elements parsed and merged, read blind per cell, held against activity in and across cells."""

from types import SimpleNamespace

from genomeos.attribution import mpra
from genomeos.predict.chromatin_tracks import WINDOW, window_for

ROWS_K562 = [
    "chr21\t100\t300\tHepG2_peak1\t212\t+\t1.5\t0.3\t0.9\t-1\t-1\n",
    "chr21\t100\t300\tHepG2_peak1_Reversed:\t212\t-\t2.5\t0.3\t0.9\t-1\t-1\n",
    "chr21\t500\t700\tK562_peak2\t236\t+\t-0.5\t0.2\t0.1\t-1\t-1\n",
    "chr22\t900\t1100\tK562_peak3\t236\t+\t0.2\t0.2\t0.1\t-1\t-1\n",
]
ROWS_HEPG2 = [
    "chr21\t100\t300\tHepG2_peak1\t212\t+\t-1.0\t0.3\t0.9\t-1\t-1\n",
    "chr21\t500\t700\tK562_peak2\t236\t+\t1.2\t0.2\t0.1\t-1\t-1\n",
]


def test_parse_merges_strands_and_cells():
    into = {}
    mpra.parse(ROWS_K562, "K562", "chr21", into)
    mpra.parse(ROWS_HEPG2, "HepG2", "chr21", into)
    els = sorted(into.values(), key=lambda e: e.start)
    assert [e.key for e in els] == ["chr21:100-300", "chr21:500-700"]
    assert els[0].name == "HepG2_peak1" and els[0].activity == {"K562": 2.0, "HepG2": -1.0}
    assert els[1].activity == {"K562": -0.5, "HepG2": 1.2}


def test_window_for_holds_the_element_whole():
    length = 5 * WINDOW
    assert window_for(10, 210, length) == 0
    assert window_for(WINDOW + 10, WINDOW + 210, length) == WINDOW
    ws = window_for(WINDOW - 50, WINDOW + 150, length)
    assert ws <= WINDOW - 50 and ws + WINDOW >= WINDOW + 150
    assert window_for(length - 100, length - 20, length) == length - WINDOW


def test_spearman():
    assert mpra.spearman([1, 2, 3, 4], [10, 20, 30, 40]) == 1.0
    assert mpra.spearman([1, 2, 3, 4], [40, 30, 20, 10]) == -1.0
    assert mpra.spearman([1, 2], [1, 2]) is None
    assert mpra.spearman([1, 1, 2, 3], [1, 2, 3, 4]) is not None


def _row(key, cls, act_k, act_h, open_k, open_h, cf, pk=None, ph=None):
    r = {
        "key": key,
        "ccre_class": cls,
        "activity": {"K562": act_k, "HepG2": act_h},
        "active": {"K562": act_k >= mpra.ACTIVE, "HepG2": act_h >= mpra.ACTIVE},
        "open": {"K562": open_k, "HepG2": open_h},
        "constrained_fraction": cf,
    }
    if pk is not None:
        r["predicted_dnase"] = {"K562": pk, "HepG2": ph}
    return r


def test_summarise_same_cell_and_across():
    rows = [
        _row("a", "dELS", 2.0, -1.0, True, False, 0.5, 0.9, 0.1),
        _row("b", "dELS", 1.5, -0.5, True, True, 0.1, 0.8, 0.2),
        _row("c", "PLS", -0.5, 1.2, False, True, 0.3, 0.2, 0.9),
        _row("d", "none", -1.0, -1.0, False, False, 0.0, 0.1, 0.1),
    ]
    s = mpra.summarise(rows)
    k = s["cells"]["K562"]
    assert k["elements"] == 4 and k["active"] == 0.5
    assert k["active_by_ccre_class"]["dELS"] == {"elements": 2, "active": 1.0}
    assert k["active_when_constrained"] == 0.5 and k["active_when_not_constrained"] == 0.5
    same = k["reader_open_predicts_active"]["same_cell"]
    assert same["precision"] == 1.0 and same["recall"] == 1.0 and same["base_rate"] == 0.5
    cross = k["reader_open_predicts_active"]["other_cells"]["HepG2"]
    assert cross["precision"] == 0.5
    pred = k["predicted_dnase_vs_activity"]
    assert pred["same_cell"]["spearman"] == 1.0 and pred["same_cell"]["active_in_top_quartile"] == 1.0
    assert pred["other_cells"]["HepG2"]["spearman"] < 0
    spec = s["specific_elements"]
    assert spec["elements"] == 3 and spec["active_in"] == {"K562": 2, "HepG2": 1}
    assert spec["predicted_dnase_higher_in_the_active_cell"] == 1.0
    assert spec["reader_open_in_the_active_cell_only"] == 0.667


def test_annotate_reads_class_reader_and_constraint():
    els = [
        mpra.Element("chr21", 100, 300, "x", {"K562": 1.5}),
        mpra.Element("chr21", 900, 1100, "y", {"K562": 0.1}),
    ]
    ccres = [SimpleNamespace(start=250, end=400, cls="dELS"), SimpleNamespace(start=260, end=290, cls="PLS")]
    idx = SimpleNamespace(covered_bp=lambda s, e: 50 if s == 100 else 0)
    stats = [SimpleNamespace(bases=200, fraction_above=0.3), SimpleNamespace(bases=0, fraction_above=None)]
    rows = mpra.annotate(els, ccres, {"K562": idx}, stats)
    assert (
        rows[0]["ccre_class"] == "PLS"
        and rows[0]["open"] == {"K562": True}
        and rows[0]["active"] == {"K562": True}
    )
    assert rows[0]["constrained_fraction"] == 0.3
    assert (
        rows[1]["ccre_class"] == "none"
        and rows[1]["open"] == {"K562": False}
        and rows[1]["constrained_fraction"] is None
    )
