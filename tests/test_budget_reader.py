# SPDX-License-Identifier: AGPL-3.0-or-later
"""R7 follow-up, consumers: every budget reader goes through budget.read_axes, and only names move."""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from genomeos.attribution import budget

BANNED = re.compile(r"\b(neutral|fossil|dead|best guess)\b")


def _stored(tmp: Path) -> dict:
    rows = [
        {
            "start": 0,
            "end": 100,
            "length": 100,
            "class": "interspersed_repeat_LINE",
            "class_evidence": "curated",
            "class_confidence": 0.8,
            "phylop": {"bases": 100, "mean": -0.2, "max": 1.0, "above": 1, "fraction_above": 0.01},
            "elements": {"n": 2, "bp": 30, "max_lod": 20, "fraction": 0.3},
            "guess": {"tier": "fossil", "label": "old label", "confidence": 0.8},
        },
        {
            "start": 200,
            "end": 300,
            "length": 100,
            "class": "unique_intergenic",
            "class_evidence": "predicted",
            "class_confidence": 0.5,
            "phylop": {"bases": 100, "mean": 0.1, "max": 2.0, "above": 0, "fraction_above": 0.0},
            "elements": {"n": 0, "bp": 0, "max_lod": 0, "fraction": 0.0},
            "guess": {"tier": "neutral", "label": "old label", "confidence": 0.6},
        },
    ]
    rec = {
        "chrom": "chrT",
        "chromosome_length": 1000,
        "unknown_bp": 200,
        "composition_bp": {"x": 1},
        "cost": {"seconds": 1.0},
        "blocks": rows,
        "by_tier": {"fossil": {"blocks": 1, "bp": 100}, "neutral": {"blocks": 1, "bp": 100}},
    }
    (tmp / "budget_chrT.json").write_text(json.dumps(rec))
    return rec


def _axes(tmp: Path, stored: dict, **change) -> None:
    rows = []
    for b in stored["blocks"]:
        g = b["guess"]
        rows.append(
            {
                "start": b["start"],
                "end": b["end"],
                "length": b["length"],
                "class": b["class"],
                "phylop": {k: b["phylop"][k] for k in ("bases", "above", "fraction_above")},
                "interspersed_coverage": 0.9,
                "guess": {
                    "tier": budget.tier_name(g["tier"]),
                    "legacy_tier": change.get("legacy", g["tier"]),
                    "label": "new label",
                    "evidence_status": "selection_not_detected",
                    "origin": "repeat_derived/LINE",
                    "confidence": g["confidence"],
                },
            }
        )
    rec = {
        "chrom": "chrT",
        "chromosome_length": 1000,
        "unknown_bp": 200,
        "blocks": rows,
        "by_tier": {budget.tier_name(t): v for t, v in stored["by_tier"].items()},
    }
    (tmp / "budget_axes_chrT.json").write_text(json.dumps(rec))


def test_the_reader_returns_the_axes_record_with_the_stored_measurements(tmp_path):
    stored = _stored(tmp_path)
    _axes(tmp_path, stored)
    r = budget.read_axes("chrT", tmp_path)
    assert r["axes_source"] == "budget_axes_chrT"
    assert set(r["by_tier"]) == {"repeat_unconstrained", "unconstrained_unknown"}
    assert r["composition_bp"] == {"x": 1} and r["cost"] == {"seconds": 1.0}
    for new, old in zip(r["blocks"], stored["blocks"], strict=True):
        assert new["guess"]["label"] == "new label" and new["guess"]["legacy_tier"] == old["guess"]["tier"]
        # every measurement a consumer reads is the stored one
        for k in ("phylop", "elements", "class_evidence", "class_confidence", "start", "end", "length"):
            assert new[k] == old[k]
        assert new["guess"]["confidence"] == old["guess"]["confidence"]


def test_the_reader_refuses_records_that_disagree(tmp_path):
    stored = _stored(tmp_path)
    _axes(tmp_path, stored, legacy="regulatory")
    with pytest.raises(ValueError, match="disagree"):
        budget.read_axes("chrT", tmp_path)


def test_without_an_axes_record_the_stored_budget_reads_under_the_new_names(tmp_path):
    stored = _stored(tmp_path)
    r = budget.read_axes("chrT", tmp_path)
    assert "not restated" in r["axes_source"]
    assert [b["guess"]["tier"] for b in r["blocks"]] == ["repeat_unconstrained", "unconstrained_unknown"]
    assert [b["guess"]["legacy_tier"] for b in r["blocks"]] == ["fossil", "neutral"]
    assert r["by_tier"]["repeat_unconstrained"] == stored["by_tier"]["fossil"]
    assert budget.read_axes("chrNone", tmp_path) is None


def test_tier_name_maps_only_the_two_legacy_keys():
    assert budget.tier_name("fossil") == "repeat_unconstrained"
    assert budget.tier_name("neutral") == "unconstrained_unknown"
    for t in ("structural", "regulatory", "constrained_unknown", "cds", None):
        assert budget.tier_name(t) == t


@pytest.mark.parametrize("chrom", ["chr21", "chrY"])
def test_on_the_committed_results_no_number_moves(chrom):
    """legacy_tier is the stored tier on every block, and the tier totals are the stored ones."""
    from genomeos.results import load_result

    stored = load_result(f"budget_{chrom}")
    r = budget.read_axes(chrom)
    if not stored or not r or r["axes_source"] != f"budget_axes_{chrom}":
        pytest.skip("committed budget results not present")
    assert [b["guess"]["legacy_tier"] for b in r["blocks"]] == [b["guess"]["tier"] for b in stored["blocks"]]
    for t, v in stored["by_tier"].items():
        assert r["by_tier"][budget.tier_name(t)]["bp"] == v["bp"]
        assert r["by_tier"][budget.tier_name(t)]["blocks"] == v["blocks"]


def test_the_budget_text_a_person_reads_names_no_function_by_absence():
    html = Path("genomeos/web/static/index.html").read_text()
    start = html.index("The 98%: composition budget")
    card = html[start : html.index("</h2>", start)]
    js = html[html.index("const col = {structural") : html.index("$('#j-budget-start')")]
    assert not BANNED.search(card) and not BANNED.search(js)
    assert not BANNED.search(budget.CONSUMERS_REGISTERED["invariant"])
