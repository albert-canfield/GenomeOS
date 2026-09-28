# SPDX-License-Identifier: AGPL-3.0-or-later
"""The HCT116 arm: HCT116 is read only when asked for, and a re-scored result is its own result, written
through save_result, never an edit of crispri_published.json (item 12 S6 follow-up)."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

from genomeos.attribution import crispri

SPEC = importlib.util.spec_from_file_location("crispri_hct116", Path("scripts/crispri_hct116.py"))
hct = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(hct)


class _Table:
    def overlapping(self, chrom, start, end):
        return [{"id": "E1", "start": start, "end": end}]


class _Cache:
    def value(self, chrom, element_id, gene, cell):
        return -0.5 if cell == "HCT116" else None


def _pair(cell: str) -> crispri.Pair:
    return crispri.Pair("chr1", 100, 200, "G", cell, "d", 5_000.0, 1.0, 1.0, True)


def test_hct116_has_no_deletion_value_unless_its_cell_is_kept():
    default, kept = _pair("HCT116"), _pair("HCT116")
    crispri.annotate([default], _Table(), _Cache())
    crispri.annotate([kept], _Table(), _Cache(), cells=(*crispri.MODEL_CELLS, "HCT116"))
    assert default.features["deletion_drop"] == 0.0 and not default.features["deletion_answered"]
    assert kept.features["deletion_drop"] == 0.5 and kept.features["deletion_answered"]


def test_kept_cells_are_the_sweeps_four_plus_hct116():
    assert hct.KEPT_CELLS[:4] == crispri.MODEL_CELLS and hct.KEPT_CELLS[4] == "HCT116"
    assert hct.ELEMENT_CACHE_HCT116 != crispri.ELEMENT_CACHE  # never written over the sweep's answers
    assert hct.CAP == 760


def test_a_rescore_is_its_own_result_and_leaves_the_headline_file_alone(tmp_path, monkeypatch):
    headline = tmp_path / "crispri_published.json"
    headline.write_text(json.dumps({"a": 1, hct.KEY: {"verdict": "first run"}}, indent=2) + "\n")
    before = headline.read_bytes()
    monkeypatch.setattr(hct, "RESULT", headline)
    manifest = {
        "sources": [{"accession": "a", "version": "1"}],
        "inputs": [{"path": "p", "sha256": "0" * 64, "partition": None}],
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {},
        "exclusions": [],
        "partitions": "n/a: test",
    }
    p = hct.write_result({"verdict": "v", "n": [1, 2]}, manifest, results_dir=tmp_path)
    assert p == tmp_path / f"{hct.NAME}.json" and headline.read_bytes() == before
    data = json.loads(p.read_text())
    assert data["verdict"] == "v" and data["n"] == [1, 2] and data["result"] == hct.NAME
    assert data["result_manifest"]["complete"] and hct.KEY in data["first_run"]
