# SPDX-License-Identifier: AGPL-3.0-or-later
"""The `save_result` cleanliness rule (the supervisor's review order, 2026-10-02).

A new result name cannot enter the registry without `result_manifest.code_cleanliness` from the one
shared function. The rule exists because the block was copy-pasted into ten files whose answers had
drifted apart: the same tree gave 41, 68 or 44 counting-path entries depending on which copy ran, and
`scripts/cell2_eligibility.py` published three paths that name no file. Nothing checked that the ten
agreed, and every published result's honesty about which code was uncommitted rests on them.

The legacy allowlist is respected: the names written before the contract keep working.
"""

import json

import planted_cleanliness
import pytest

from genomeos import manifest as mf
from genomeos import results
from genomeos.results import save_result


def full(tmp_path):
    """A manifest that meets every part of the contract except the cleanliness block."""
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


@pytest.fixture
def registry(tmp_path, monkeypatch):
    """A registry at tmp/data/results with a one-name legacy allowlist beside it."""
    reg = tmp_path / "data" / "results"
    reg.mkdir(parents=True)
    allow = tmp_path / "data" / "results_legacy.txt"
    allow.write_text("# a header line\nold_result\ttracked\n")
    monkeypatch.setattr(results, "RESULTS_DIR", reg)
    monkeypatch.setattr(results, "LEGACY_ALLOWLIST", allow)
    return reg


# (a) the refusal ---------------------------------------------------------------------------------------


def test_a_registry_write_without_cleanliness_is_refused_and_the_quarantine_holds_the_payload(
    registry, tmp_path
):
    with pytest.raises(mf.ManifestError, match="quarantined at") as e:
        save_result("brand_new", {"value": 7}, registry, manifest=full(tmp_path))
    assert "missing code_cleanliness" in str(e.value)
    assert not (registry / "brand_new.json").exists(), "nothing may be written to the registry"

    q = results.quarantine_dir(registry) / "brand_new.json"
    kept = json.loads(q.read_text())
    assert kept["value"] == 7, "the payload is kept, so hours of compute are not lost"
    assert kept[mf.KEY]["complete"] is False
    assert any("missing code_cleanliness" in p for p in kept["quarantine"]["problems"])
    assert kept["quarantine"]["meant_for"] == str(registry / "brand_new.json")


def test_the_refusal_names_the_shared_function_so_the_writer_knows_what_to_call(registry, tmp_path):
    with pytest.raises(mf.ManifestError, match=r"genomeos\.manifest\.code_cleanliness"):
        save_result("brand_new", {"value": 1}, registry, manifest=full(tmp_path))


def test_a_cleanliness_block_built_by_hand_is_refused_by_key_set(registry, tmp_path):
    """The point of the rule: a block that did not come from the shared function cannot be trusted to
    have counted the same things, so a near-miss is refused rather than published."""
    m = full(tmp_path)
    m["code_cleanliness"] = {"git_sha": "abc123", "dirty": False, "counting_path": ["a.py"]}
    with pytest.raises(mf.ManifestError, match="code_cleanliness lacks") as e:
        save_result("brand_new", {"value": 1}, registry, manifest=m)
    assert "own_code_is_committed" in str(e.value)
    assert not (registry / "brand_new.json").exists()


def test_a_complete_manifest_with_the_shared_block_writes(registry, tmp_path):
    m = full(tmp_path)
    # 2026-10-02: the entry script is counted now (genomeos.manifest.entry_script), so under pytest
    # `sys.argv[0]` is a test runner and no in-process call can return a passing block. This
    # scaffolding block is built in a planted repository whose entry script really is committed and
    # clean, which is what it always meant to assert, instead of naming a file of this shared
    # checkout and reddening when a peer edits it (tests/planted_cleanliness.py).
    m["code_cleanliness"] = planted_cleanliness.planted_clean_block(tmp_path, "refusal_clean")
    p = save_result("brand_new", {"value": 1}, registry, manifest=m)
    assert p.exists() and json.loads(p.read_text())[mf.KEY]["complete"] is True


# (b) the legacy allowlist still writes -----------------------------------------------------------------


def test_a_legacy_allowlist_name_still_writes_without_cleanliness(registry, tmp_path):
    """The names written before the contract must keep working: the rule applies to a name that is not
    on the allowlist, so a legacy name writes with no cleanliness block and is not quarantined."""
    p = save_result("old_result", {"value": 1}, registry, manifest=full(tmp_path))
    assert p.exists() and json.loads(p.read_text())["value"] == 1
    assert not (results.quarantine_dir(registry) / "old_result.json").exists()


def test_a_legacy_name_with_an_incomplete_manifest_warns_rather_than_fails(registry):
    """A legacy result warns, never fails, so the historical writers keep regenerating."""
    with pytest.warns(results.ManifestWarning):
        p = save_result("old_result", {"value": 1}, registry)
    assert p.exists()


def test_the_rule_does_not_reach_outside_the_registry(tmp_path):
    """A scratch directory is not the registry: the rule is about what enters the registry."""
    out = tmp_path / "scratch"
    out.mkdir()
    with pytest.warns(results.ManifestWarning):
        p = save_result("anything", {"value": 1}, out)
    assert p.exists()


# the contract itself ----------------------------------------------------------------------------------


def test_cleanliness_problems_accepts_what_the_shared_function_returns():
    block = mf.code_cleanliness("tests/test_results_cleanliness_refusal.py", ())
    assert set(block) >= set(mf.CLEANLINESS_KEYS)
    assert mf.cleanliness_problems({"code_cleanliness": block}) == []


@pytest.mark.parametrize("bad", [None, "a string", 42, []])
def test_cleanliness_problems_refuses_a_block_that_is_not_a_dict(bad):
    assert mf.cleanliness_problems({"code_cleanliness": bad}) != []
