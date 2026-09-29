# SPDX-License-Identifier: AGPL-3.0-or-later
"""The pilot's prior-only lead, tested on its own (lane-prior, 2026-09-29): the admissibility record.
Synthetic rows only: no benchmark file, no cache, no request."""

from __future__ import annotations

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("pot", ROOT / "scripts" / "prior_only_test.py")
pot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pot)

HARNESS_STUDIES = (
    "Gasperini2019",
    "HCT116",
    "K562_DC_TAP",
    "Klann",
    "Morris",
    "Nasser2021",
    "Reilly",
    "Schraivogel2020",
    "WTC11_DC_TAP",
    "Xie",
)


# --- negatives first -------------------------------------------------------------------------------


def test_every_held_out_file_study_is_inadmissible():
    # the benchmark kept a held-out positive only with H3K27ac at the element; the prior reads H3K27ac
    for src, v in pot.STUDIES.items():
        if v["file"] == pot.HELDOUT_FILE or pot.HELDOUT_FILE in v["file"]:
            assert v["class"] == pot.INADMISSIBLE, src


def test_the_other_half_of_the_lead_is_not_admissible():
    assert pot.STUDIES["crispri:Xie"]["class"] == pot.INADMISSIBLE
    assert pot.STUDIES["crispri:Morris"]["class"] == pot.INADMISSIBLE
    assert "crispri:Xie" not in pot.ADMISSIBLE_SOURCES


def test_no_admissible_study_is_unqualified():
    # both admissible screens chose their tested candidates on chromatin: the selection is stated
    assert set(pot.ADMISSIBLE_SOURCES) == {"crispri:Gasperini2019", "crispri:Schraivogel2020"}
    for src in pot.ADMISSIBLE_SOURCES:
        assert pot.STUDIES[src]["class"] == pot.ADMISSIBLE_STATED


# --- the record ------------------------------------------------------------------------------------


def test_every_harness_study_has_a_record_with_a_quoted_source():
    assert set(pot.STUDIES) == {f"crispri:{s}" for s in HARNESS_STUDIES}
    for src, v in pot.STUDIES.items():
        assert v["class"] in (pot.INADMISSIBLE, pot.ADMISSIBLE_STATED, pot.ADMISSIBLE), src
        assert v["quotes"], src
        for q in v["quotes"]:
            assert q["text"] and q["url"].startswith("https://"), src


def test_the_benchmark_evidence_names_the_filter_and_its_reproduction():
    held = pot.BENCHMARK_EVIDENCE["held_out_file"]
    assert any("no H3K27ac categories" in q["text"] for q in held["quotes"])
    assert "4,378 of 4,378" in held["reproduced_by_this_lane"]
    assert "no chromatin filter" in pot.BENCHMARK_EVIDENCE["training_file"]["filters"]


def test_local_counts_count_positives_outside_the_h3k27ac_categories():
    rows = {
        "crispri:A": [
            {"_file": "f", "Regulated": "TRUE", "elementChromatinCategory": "High H3K27ac"},
            {"_file": "f", "Regulated": "TRUE", "elementChromatinCategory": "CTCF element"},
            {"_file": "f", "Regulated": "FALSE", "elementChromatinCategory": "No H3K27ac"},
            {"_file": "f", "Regulated": "FALSE", "elementChromatinCategory": "H3K27ac"},
        ]
    }
    c = pot.local_counts(rows)["crispri:A"]["f"]
    assert (c["pairs"], c["positives"]) == (4, 2)
    assert c["positives_outside_h3k27ac_categories"] == 1
    assert c["pairs_outside_h3k27ac_categories"] == 2


def test_partition_counts_apply_the_harness_floors():
    ho = pot.ho
    ms = pot.ms
    units = []
    for i in range(30):
        outcome = ms.DECREASE if i < 21 else ms.NULL_INFORMATIVE
        units.append(ho.Unit("crispri:A", "chr2", i * 2_000_000, i * 2_000_000 + 500, outcome, gene=f"G{i}"))
    units.append(ho.Unit("crispri:A", "chr1", 0, 500, ms.DECREASE, gene="X"))
    got = pot.partition_counts({"crispri:A": tuple(units)}, ("chr2",))["crispri:A"]
    assert (got["units"], got["positives"], got["negatives"]) == (30, 21, 9)
    assert got["clears_floors"] is False  # 9 negatives, below the floor of 20
