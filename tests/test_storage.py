"""Stream, distil, discard: results summaries stand in for raw data."""

from pathlib import Path

import pytest

from genomeos import storage
from genomeos.results import list_results, load_result, save_result


def test_results_roundtrip(tmp_path):
    p = save_result("demo", {"value": 1}, results_dir=tmp_path)
    assert p.exists() and load_result("demo", tmp_path)["value"] == 1
    assert list_results(tmp_path)[0]["name"] == "demo"


def test_status_and_manifest():
    st = storage.status()
    assert "disk" in st and st["disk"]["free"] > 0
    m = storage.manifest()
    assert {"clock_GSE41169", "clinvar_chr21_agreement", "library_members", "hg002_chr21"} <= set(m)


def test_clean_is_a_dry_run_by_default():
    rows = storage.clean(dry_run=True)
    assert all(action == "would delete" for _, _, action in rows)
    for path, _, _ in rows:
        assert Path(path).exists()


@pytest.mark.skipif(load_result("hg002_telomere_stream") is None, reason="no streamed telomere result yet")
def test_streamed_hg002_telomere_is_plausible():
    r = load_result("hg002_telomere_stream")
    assert r["disk_used_bytes"] == 0 and r["reads"] >= 1_000_000
    # TelSeq-scale lengths for adults are ~2-5 kb (about half of Southern-blot lengths)
    assert 1500 < r["telomere_bp"] < 8000, r["telomere_bp"]
    assert r["southern_equivalent_bp_inferred"] == 2 * r["telomere_bp"]


@pytest.mark.skipif(load_result("clock_GSE41169") is None, reason="distil first")
def test_distilled_clock_result_holds():
    r = load_result("clock_GSE41169")
    assert r["horvath"]["r"] > 0.8 and r["horvath"]["mae_years"] < 10 and r["samples"] >= 50


@pytest.mark.skipif(load_result("clinvar_chr21_agreement") is None, reason="distil first")
def test_distilled_clinvar_result_holds():
    r = load_result("clinvar_chr21_agreement")
    for k, v in r["agreement"].items():
        assert v["rate"] > 0.9, (k, v)


@pytest.mark.skipif(load_result("library_members") is None, reason="distil first")
def test_distilled_library_result_holds():
    r = load_result("library_members")
    assert r["agreement"] > 0.94 and "BRCA1" in r["members"]["core.dna_repair"]
