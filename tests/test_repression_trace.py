# SPDX-License-Identifier: AGPL-3.0-or-later
"""The repression trace's registration, pinned before the trace was written (lane-repress, 2026-09-29)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

_spec = importlib.util.spec_from_file_location("repression_trace", ROOT / "scripts/repression_trace.py")
rt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rt)


def test_the_brief_classes_are_registered_and_only_f_was_added_after_reading() -> None:
    assert list(rt.CLASSES) == ["a", "b", "c", "d", "e", "f"]
    late = [k for k, c in rt.CLASSES.items() if c["added_after_first_reading"]]
    assert late == ["f"]
    for c in rt.CLASSES.values():
        assert c["name"] and c["definition"] and c["decided_by"]
    reg = rt.registration()
    assert reg["registered"] == "2026-09-29"
    assert "not a validation" in reg["status"]
    assert "not blind" in reg["written_after"]


def test_the_four_cases_and_the_comparison_are_the_ones_s4_judged_wrong() -> None:
    assert [c["gene"] for c in rt.CASES] == ["CD83", "HEMGN", "ID1", "BEX4"]
    assert all(c["value"] == "represses_target" for c in rt.CASES)
    assert all(c["cell"] != c["deciding"]["cell"] for c in rt.CASES)
    assert rt.COMPARISON["gene"] == "CCND1" and rt.COMPARISON["value"] == "activates_target"


def test_primary_rule_order() -> None:
    yes, no, unk = rt.PRESENT, rt.ABSENT, rt.UNKNOWN
    base = {"a": no, "b": yes, "c": no, "d": yes, "e": no, "f": yes}
    assert rt.primary(base) == "b"
    assert rt.primary({**base, "c": unk}) == "b"
    assert rt.primary({**base, "c": yes}) == "c"
    assert rt.primary({**base, "c": yes, "e": yes}) == "e"
    assert rt.primary({**base, "c": yes, "e": yes, "a": yes}) == "a"
    with pytest.raises(ValueError):
        rt.primary({"a": no, "b": no, "c": no, "d": yes, "e": no, "f": yes})  # (d) and (f) are never primary
