# SPDX-License-Identifier: AGPL-3.0-or-later
"""The Cells view's data: eleven cell types reading one genome."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genomeos.web.server import Api, ApiError


def _write(root: Path, payload: dict) -> None:
    (root / "data" / "results").mkdir(parents=True, exist_ok=True)
    (root / "data" / "results" / "reader_genome_wide.json").write_text(json.dumps(payload))


def test_every_cell_gets_its_own_read_fraction_and_chromosomes(tmp_path: Path) -> None:
    _write(
        tmp_path,
        {
            "date": "2026-09-14",
            "evidence": "experimental: ENCODE DNase-seq",
            "cell_types": ["K562", "hepatocyte"],
            "totals": {
                "K562": {"coding_genes": 100, "genes_read": 60, "genes_read_open": 70},
                "hepatocyte": {"coding_genes": 100, "genes_read": 80, "genes_read_open": 90},
            },
            "chromosomes": {
                "chr1": {"K562": {"genes_read": 6}, "hepatocyte": {"genes_read": 8}},
                "chr2": {"K562": {"genes_read": 5}},
            },
        },
    )

    out = Api(tmp_path).cells()

    assert [c["cell_type"] for c in out["cells"]] == ["hepatocyte", "K562"]  # most read first
    assert out["cells"][0]["read_fraction"] == 0.8
    assert out["cells"][0]["open_fraction"] == 0.9
    assert sorted(out["cells"][1]["chromosomes"]) == ["chr1", "chr2"]
    assert sorted(out["cells"][0]["chromosomes"]) == ["chr1"]  # absent from chr2, not invented as zero
    assert out["chromosomes"] == ["chr1", "chr2"]


def test_a_cell_with_no_coding_genes_reports_no_fraction_rather_than_zero(tmp_path: Path) -> None:
    """A missing denominator is not a read fraction of 0: that would rank it as the least read."""
    _write(
        tmp_path,
        {"cell_types": ["X"], "totals": {"X": {"genes_read": 5}}, "chromosomes": {}},
    )

    out = Api(tmp_path).cells()

    assert out["cells"][0]["read_fraction"] is None


def test_it_refuses_rather_than_returning_an_empty_page(tmp_path: Path) -> None:
    (tmp_path / "data" / "results").mkdir(parents=True)

    with pytest.raises(ApiError, match="reader_genome_wide"):
        Api(tmp_path).cells()
