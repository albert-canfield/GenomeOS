# SPDX-License-Identifier: AGPL-3.0-or-later
"""The compiled-genome evidence counts call a fact weak only when it STATES a low confidence (R4d).

Before, scripts/compile_genome_programs.py and scripts/measured_layer.py compared `confidence <= 0.5`
on every row, and since R4 the parser reads every predicted fact's missing confidence as 0.0, so
facts nobody judged were counted weak. evidence.py's rows carry `stated`; both scripts now use it.
"""

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def _script(name: str):
    spec = importlib.util.spec_from_file_location(f"_{name}", ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SCRIPTS = ("compile_genome_programs", "measured_layer")


def _rows(compiled: str) -> list[dict]:
    p = f"{compiled}/noncoding_chrT.bio"
    return [
        {"path": p, "evidence": "predicted", "confidence": 0.0, "stated": False},  # unstated
        {"path": p, "evidence": "predicted", "confidence": 0.0, "stated": False},  # unstated
        {"path": p, "evidence": "curated", "confidence": 0.0, "stated": True},  # a stated 0.0 is weak
        {"path": p, "evidence": "inferred", "confidence": 0.4, "stated": True},  # weak
        {"path": p, "evidence": "experimental", "confidence": 0.9, "stated": True},  # strong
    ]


@pytest.mark.parametrize("name", SCRIPTS)
def test_an_unstated_confidence_is_never_counted_weak(name, monkeypatch):
    mod = _script(name)
    rows = _rows(str(mod.COMPILED_DIR))
    monkeypatch.setattr(mod, "collect", lambda root, compiled: {"rows": rows})
    got = mod.read_evidence()
    for tally in (got["pooled"], got["per_chromosome"]["chrT"]):
        assert tally["unstated"] == 2 and tally["stated"] == 3
        assert tally["weak"] == 2 and tally["strong"] == 1
        assert tally["weak"] + tally["strong"] + tally["unstated"] == tally["facts"] == 5
        assert tally["weak_before_r4_counting_unstated"] == 4  # the old rule, kept for comparison
    assert "STATED" in got["weak_means"]


@pytest.mark.parametrize("name", SCRIPTS)
def test_a_row_without_the_flag_reads_as_stated(name):
    """Rows from before the flag existed were all stated; the default keeps their count."""
    mod = _script(name)
    tally = dict.fromkeys(mod.STATED_COUNTS, 0)
    mod.count_stated(tally, {"confidence": 0.3})
    assert tally["weak"] == 1 and tally["unstated"] == 0
