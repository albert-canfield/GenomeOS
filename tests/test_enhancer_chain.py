"""The chromosome chain decides completion from the scoring script's own result flag."""

import importlib.util
import json
from pathlib import Path


def _chain():
    spec = importlib.util.spec_from_file_location("chain", Path("scripts/enhancer_targets_all_chain.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_complete_and_scored_read_the_result_flag(tmp_path, monkeypatch):
    chain = _chain()
    monkeypatch.setattr(chain, "RESULTS", tmp_path)
    assert chain.complete("chr1") is False and chain.scored("chr1") is None
    (tmp_path / "enhancer_targets_all_chr1.json").write_text(
        json.dumps({"scored": 200, "elements_total": 900, "complete": False})
    )
    assert chain.complete("chr1") is False and chain.scored("chr1") == (200, 900)
    (tmp_path / "enhancer_targets_all_chr1.json").write_text(
        json.dumps({"scored": 900, "elements_total": 900, "complete": True})
    )
    assert chain.complete("chr1") is True
    (tmp_path / "enhancer_targets_all_chr2.json").write_text("{not json")
    assert chain.complete("chr2") is False and chain.scored("chr2") is None
    assert chain.ORDER[0] == "chr1" and chain.ORDER[-2:] == ["chrX", "chrY"] and "chr21" in chain.ORDER
