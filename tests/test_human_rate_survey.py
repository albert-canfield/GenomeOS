# SPDX-License-Identifier: AGPL-3.0-or-later
"""Review R3, the second rate source: the survey for a HUMAN absolute transcription rate.

The lane's result is a negative, so these tests guard a negative. They check that the survey states,
for every candidate, what would let a reader disagree (species, cell, units, licence, URL and the
reason it cannot supply T); that no candidate is both human and an absolute per-gene rate, which is
the finding itself and the test that will fail the day one exists; that the rate table says mouse on
every ROW rather than once in a header; and that the audit counts zero pairs standing on a human rate.

The registered invariant is that a survey moves no number: the counts below are the ones the measured
tier already had, and a change in any of them is a defect, not a finding.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from genomeos.attribution import bridge

RESULTS = Path("data/results")
FIELDS = (
    "citation",
    "species",
    "cell",
    "quantity",
    "units",
    "absolute",
    "n_genes",
    "licence",
    "url",
    "why_not",
)
#: the registered invariant, from data/results/bridge_audit.json as lane-rates left it
MEASURED_PAIRS = 22576
BORROWED_PAIRS = 160392
GENE_RATE_UNMEASURED = 138543
NOT_IDENTIFIABLE = 76469


def _result(name: str) -> dict:
    p = RESULTS / f"{name}.json"
    if not p.exists():
        pytest.skip(f"{p} not built")
    return json.loads(p.read_text())


# ---- the survey states enough to be argued with ---------------------------------------------------


def test_every_candidate_states_species_units_licence_and_a_reason():
    assert len(bridge.HUMAN_RATE_SOURCE_SURVEY) >= 8
    for key, rec in bridge.HUMAN_RATE_SOURCE_SURVEY.items():
        for f in FIELDS:
            assert f in rec, f"{key} does not state {f}"
            assert str(rec[f]).strip(), f"{key}.{f} is empty"
        assert rec["url"].startswith("https://"), key


def test_no_candidate_is_both_human_and_an_absolute_per_gene_rate():
    """The finding, as a test. It fails the day a human absolute source is added, which is the point."""
    for key, rec in bridge.HUMAN_RATE_SOURCE_SURVEY.items():
        if str(rec["absolute"]).strip().lower() == "yes":
            assert "Homo" not in rec["species"], (
                f"{key} claims an absolute rate for a human cell: the survey's verdict is refuted and"
                " the lane reopens"
            )


def test_the_one_absolute_source_found_is_mouse_and_undistributed():
    absolute = [k for k, r in bridge.HUMAN_RATE_SOURCE_SURVEY.items() if str(r["absolute"]).lower() == "yes"]
    assert absolute == ["shao2022"]
    rec = bridge.HUMAN_RATE_SOURCE_SURVEY["shao2022"]
    assert rec["species"] == "Mus musculus"
    assert "cell^-1 min^-1" in rec["units"]


def test_a_half_life_source_is_recorded_as_unable_to_supply_a_rate():
    rec = bridge.HUMAN_RATE_SOURCE_SURVEY["schofield2018"]
    assert "A HALF-LIFE CANNOT SUPPLY A SYNTHESIS RATE" in rec["why_not"]
    assert "construction" in bridge.HUMAN_RATE_CONSTRUCTION_NOT_TAKEN


# ---- the rate table says mouse on every row -------------------------------------------------------


def test_every_rate_row_carries_its_source_and_species_and_none_is_human():
    d = _result("gene_rates")
    assert d["rates"], "no rates in the table"
    for gene, row in d["rates"].items():
        assert row["source"] in bridge.MEASURED_RATE_SOURCES, gene
        assert row["species_cell"].startswith("Mus musculus"), gene
    assert d["human_measured_rates"] == 0
    assert d["rates_by_source"] == {"schwanhausser2011": len(d["rates"])}


def test_no_rate_in_the_table_is_the_copy_number_over_the_half_life_construction():
    """The construction the lane refused, checked in the data: T is not N * ln2 / t_half."""
    d = _result("gene_rates")
    diffs = [
        abs(r["transcription_rate"] - r["mrna_copies"] * math.log(2) / r["mrna_half_life_h"])
        / r["transcription_rate"]
        for r in d["rates"].values()
        if r["mrna_copies"] and r["mrna_half_life_h"] and r["transcription_rate"]
    ]
    assert len(diffs) > 1000
    diffs.sort()
    assert diffs[len(diffs) // 2] > 0.1, "the rates would be indistinguishable from the construction"


# ---- the audit reports the species of the rate each pair stands on --------------------------------


def test_the_audit_counts_no_pair_on_a_human_rate_and_every_pair_on_a_borrowed_one():
    d = _result("bridge_audit")
    assert d["genes_with_a_human_measured_rate"] == 0
    for tier, counts in d["pooled_by_tier"].items():
        assert counts["pairs_on_a_human_measured_rate"] == 0, tier
        assert counts["pairs_on_a_borrowed_species_rate"] == counts["genes_on_their_own_measured_rate"], tier


def test_the_survey_moved_no_count():
    d = _result("bridge_audit")
    assert d["pooled_by_tier"]["measured"]["simulable_gene_cell_pairs"] == MEASURED_PAIRS
    assert d["pooled_by_tier"]["borrowed"]["simulable_gene_cell_pairs"] == BORROWED_PAIRS
    assert d["unresolved_by_tier"]["measured"]["gene_rate_unmeasured"] == GENE_RATE_UNMEASURED
    for tier in ("none", "measured", "borrowed"):
        assert d["unresolved_by_tier"][tier]["not_identifiable"] == NOT_IDENTIFIABLE, tier
