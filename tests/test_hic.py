# SPDX-License-Identifier: AGPL-3.0-or-later
"""4DN boundary files as measured node edges (no network: files written by hand)."""

from genomeos.genome import hic


def test_load_and_compare(tmp_path, monkeypatch):
    d = tmp_path / "GM12878"
    d.mkdir()
    (d / "chr21.bed").write_text("chr21\t100000\t110000\nchr21\t500000\t510000\nchr21\t900000\t910000\n")
    (d / "manifest.json").write_text('{"accession": "x", "chromosomes": {"chr21": 3}}')
    bounds = hic.load_boundaries("GM12878", "chr21", tmp_path)
    assert bounds == [105000, 505000, 905000]
    c = hic.compare_with_inferred([104000, 300000, 907000], bounds, 1_000_000)
    assert c["inferred_on_a_predicted_boundary"] == 2 and c["fraction_inferred_supported"] == 0.667
    assert 0.0 <= c["random_control"] <= 1.0
    monkeypatch.delenv("FOURDN_KEY", raising=False)
    monkeypatch.delenv("FOURDN_SECRET", raising=False)
    monkeypatch.setattr(hic, "credentials", lambda: None)
    assert hic.status()["enabled"] is False and "4DN" in hic.status()["how"]
    try:
        hic.fetch_boundaries("K562", tmp_path)
    except PermissionError as ex:
        assert "account key" in str(ex)
    else:  # pragma: no cover
        raise AssertionError("a download without a key must be refused before any request")
