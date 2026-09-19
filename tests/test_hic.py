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


def test_authorization_is_dropped_when_the_redirect_leaves_the_portal():
    import urllib.request

    from genomeos.genome.hic import _DropAuthOnRedirect

    h = _DropAuthOnRedirect()
    req = urllib.request.Request(
        "https://data.4dnucleome.org/files/x", headers={"Authorization": "Basic abc"}
    )
    to_s3 = h.redirect_request(
        req, None, 307, "Temporary Redirect", {}, "https://s3.amazonaws.com/b/x.bed.gz"
    )
    assert to_s3 is not None and not to_s3.has_header("Authorization")
    same = h.redirect_request(req, None, 307, "Temporary Redirect", {}, "https://data.4dnucleome.org/files/y")
    assert same is not None and same.has_header("Authorization")


def test_a_boundary_keeps_the_strength_its_call_was_thresholded_on(tmp_path):
    """4DN's files carry a label and a strength; reading the call without the quantity loses the part
    that separates a strong boundary from a marginal one, which is what the rearrangement test needed.
    A file written before the strength was kept must report None rather than zero, because zero would
    read as 'measured and weak' when the truth is 'not recorded'."""
    from genomeos.genome import hic

    d = tmp_path / "GM12878"
    d.mkdir(parents=True)
    (d / "chr2.bed").write_text(
        "chr2\t100\t200\tStrong\t1.90\n"  # scored, as files written since the fix are
        "chr2\t400\t500\tWeak\t0.57\n"
        "chr2\t700\t800\n"  # written before the fix: interval only
    )
    rows = hic.load_scored_boundaries("GM12878", "chr2", tmp_path)
    assert [r["centre"] for r in rows] == [150, 450, 750]
    assert [r["label"] for r in rows] == ["Strong", "Weak", None]
    assert [r["strength"] for r in rows] == [1.90, 0.57, None]
    # the positions the rest of the project reads are unchanged by carrying the extra columns
    assert hic.load_boundaries("GM12878", "chr2", tmp_path) == [150, 450, 750]
