# SPDX-License-Identifier: AGPL-3.0-or-later
"""The HCT116 arm: HCT116 is read only when asked for, and the result is written without changing a line."""

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


def test_result_is_inserted_without_changing_an_existing_line(tmp_path, monkeypatch):
    f = tmp_path / "r.json"
    original = {"a": 1, "second_cell_type": {"run": False}, "post_hoc_positive_filter": {"x": 2}, "z": 0}
    f.write_text(json.dumps(original, indent=2) + "\n")
    monkeypatch.setattr(hct, "RESULT", f)
    before = f.read_text().splitlines()
    hct.write_additively({"verdict": "v", "n": [1, 2]})
    after = f.read_text().splitlines()
    assert [line for line in before if line not in after] == []
    data = json.loads(f.read_text())
    assert data[hct.KEY] == {"verdict": "v", "n": [1, 2]} and data["z"] == 0
