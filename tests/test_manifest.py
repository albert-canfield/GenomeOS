"""Result manifests (review item R9): the contract, its enforcement for new results, and the
tolerant reader that keeps the 955 historical results readable."""

import json
import warnings
from pathlib import Path

import pytest

from genomeos import manifest as mf
from genomeos import results
from genomeos.results import ManifestWarning, load_manifest, save_result


def full(tmp_path: Path) -> dict:
    src = tmp_path / "in.tsv"
    src.write_text("chrom\tstart\tend\nchr21\t0\t10\n")
    return {
        "sources": [{"accession": "ENCSR000XXX", "version": "v1"}],
        "inputs": [mf.input_entry(src, partition="heldout")],
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {"threshold": 0.5},
        "exclusions": [],
        "partitions": {"heldout": "the test"},
    }


def test_a_complete_manifest_is_stamped_with_the_code_revision(tmp_path):
    p = save_result("demo", {"value": 1}, tmp_path, manifest=full(tmp_path), strict=True)
    m = json.loads(p.read_text())[mf.KEY]
    assert m["complete"] and "problems" not in m
    assert len(m["code"]["git_sha"]) == 40 and isinstance(m["code"]["dirty"], bool)
    assert m["inputs"][0]["partition"] == "heldout" and len(m["inputs"][0]["sha256"]) == 64


def test_the_caller_cannot_state_its_own_code_revision(tmp_path):
    given = {**full(tmp_path), "code": {"git_sha": "0" * 40, "dirty": False}}
    p = save_result("demo", {}, tmp_path, manifest=given, strict=True)
    assert json.loads(p.read_text())[mf.KEY]["code"]["git_sha"] != "0" * 40


def test_a_new_result_without_the_contract_is_quarantined_then_refused(tmp_path):
    """Item 12 S6: until 2026-09-28 this test asserted the defect, that the failed result was written
    where it was meant to go. It is kept, but in the quarantine."""
    reg = tmp_path / "results"
    with pytest.raises(mf.ManifestError, match="missing sources"):
        save_result("demo", {"value": 1}, reg, strict=True)
    assert not (reg / "demo.json").exists()
    kept = json.loads((results.quarantine_dir(reg) / "demo.json").read_text())
    assert kept["value"] == 1 and kept[mf.KEY]["complete"] is False  # the compute is not lost


def test_strict_is_the_default_for_a_name_off_the_allowlist_in_the_registry_only(tmp_path, monkeypatch):
    reg = tmp_path / "results"
    monkeypatch.setattr(results, "RESULTS_DIR", reg)
    monkeypatch.setattr(results, "LEGACY_ALLOWLIST", tmp_path / "legacy.txt")
    (tmp_path / "legacy.txt").write_text("# header\nold_one\ttracked\n")
    with pytest.raises(mf.ManifestError):
        save_result("brand_new", {"value": 1}, reg)
    # a retry is refused again: before item 12 S6 it found its own file and only warned
    with pytest.raises(mf.ManifestError):
        save_result("brand_new", {"value": 2}, reg)
    assert not (reg / "brand_new.json").exists()
    # a name on the allowlist warns: historical results keep regenerating
    with pytest.warns(ManifestWarning, match="manifest incomplete"):
        save_result("old_one", {"value": 2}, reg)
    # outside the registry (tests, scratch) a write warns rather than fails
    other = tmp_path / "elsewhere"
    with pytest.warns(ManifestWarning):
        save_result("x", {}, other)


def test_the_manifest_can_travel_in_the_payload(tmp_path):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        p = save_result("demo", {"value": 1, mf.KEY: full(tmp_path)}, tmp_path, strict=True)
    d = json.loads(p.read_text())
    assert d[mf.KEY]["complete"] and list(d)[-1] == mf.KEY


def test_an_older_manifest_key_is_left_alone(tmp_path):
    """26 historical results use `manifest` for a file name; the contract lives under its own key."""
    with pytest.warns(ManifestWarning):
        p = save_result("demo", {"manifest": "epigenome_manifest"}, tmp_path)
    assert json.loads(p.read_text())["manifest"] == "epigenome_manifest"


@pytest.mark.parametrize(
    "change, problem",
    [
        ({"coordinates": {"base": 2, "interval": "half-open"}}, "coordinates"),
        ({"coordinates": "n/a"}, "coordinates"),
        ({"sources": [{"accession": "X"}]}, "accession and version"),
        ({"inputs": [{"path": "a", "sha256": "b"}]}, "partition"),
        ({"exclusions": None}, "exclusions"),
        ({"partitions": "none"}, "partitions"),
    ],
)
def test_validate_names_the_fault(tmp_path, change, problem):
    m = {**full(tmp_path), "code": {"git_sha": "a", "dirty": False}, **change}
    assert any(problem in p for p in mf.validate(m))


def test_not_applicable_must_say_why(tmp_path):
    m = {**full(tmp_path), "code": {"git_sha": "a", "dirty": False}}
    m.update(coordinates="n/a: a count per cell type has no coordinates", partitions="n/a: no evaluation")
    assert mf.validate(m) == []


def test_a_directory_input_hashes_the_same_wherever_it_sits(tmp_path):
    for d in ("a", "b"):
        (tmp_path / d / "sub").mkdir(parents=True)
        (tmp_path / d / "sub" / "x.json").write_text("[1]")
        (tmp_path / d / "y.json").write_text("[2]")
    assert mf.sha256_of(tmp_path / "a") == mf.sha256_of(tmp_path / "b")
    (tmp_path / "b" / "y.json").write_text("[3]")
    assert mf.sha256_of(tmp_path / "a")[0] != mf.sha256_of(tmp_path / "b")[0]


def test_the_reader_tolerates_every_historical_shape(tmp_path):
    assert mf.read([1, 2])["declared"] is False
    legacy = mf.read({"result": "x", "genome": "data/reference/chr21.fa.gz", "thresholds": {}})
    assert not legacy["declared"] and not legacy["complete"]
    assert legacy["fields"]["assembly"] == "legacy" and legacy["fields"]["parameters"] == "legacy"
    assert legacy["fields"]["code"] is None
    save_result("demo", {}, tmp_path, manifest=full(tmp_path), strict=True)
    r = load_manifest("demo", tmp_path)
    assert r["declared"] and r["complete"] and set(r["fields"].values()) == {"declared"}
    assert load_manifest("absent", tmp_path) is None


def test_the_census_reads_the_real_registry_without_rewriting_it():
    import importlib.util

    root = Path(__file__).resolve().parent.parent
    spec = importlib.util.spec_from_file_location("manifest_census", root / "scripts" / "manifest_census.py")
    census = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(census)
    before = {p: p.stat().st_mtime_ns for p in (root / "data" / "results").glob("*.json")}
    c = census.census(root / "data" / "results")
    assert c["results"] == len(before) and not c["unreadable"]
    assert set(c["fields"]) == set(mf.REQUIRED)
    for n in c["fields"].values():
        assert n["declared"] + n["legacy"] + n["none"] == c["results"]
    assert {p: p.stat().st_mtime_ns for p in before} == before


def test_the_rebuild_compares_field_by_field():
    import importlib.util

    root = Path(__file__).resolve().parent.parent
    spec = importlib.util.spec_from_file_location(
        "manifest_rebuild", root / "scripts" / "manifest_rebuild.py"
    )
    rb = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rb)
    a = {"date": "2026-01-01", "x": [1, {"y": 2}], mf.KEY: {"code": {"git_sha": "a"}, "assembly": "GRCh38"}}
    b = {"date": "2026-09-28", "x": [1, {"y": 3}], mf.KEY: {"code": {"git_sha": "b"}, "assembly": "GRCh38"}}
    assert rb.diff(rb.comparable(a), rb.comparable(b)) == ["/x[1]/y: 2 vs 3"]
