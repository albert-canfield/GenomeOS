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


# --- the traced causes, pinned from the result (lane-repress, 2026-09-29) -------------------------------
RESULT = ROOT / "data/results/repression_trace.json"


def _result() -> dict:
    import json

    return json.loads(RESULT.read_text())


def test_the_four_traced_causes() -> None:
    r = _result()
    got = {t["gene"]: (t["primary"], t["classes"]) for t in r["four"]}
    absent, present, unknown = rt.ABSENT, rt.PRESENT, rt.UNKNOWN
    assert got == {
        "CD83": ("b", {"a": absent, "b": present, "c": absent, "d": absent, "e": absent, "f": present}),
        "HEMGN": ("b", {"a": absent, "b": present, "c": absent, "d": present, "e": absent, "f": present}),
        "ID1": ("c", {"a": absent, "b": present, "c": present, "d": present, "e": absent, "f": present}),
        "BEX4": ("b", {"a": absent, "b": present, "c": absent, "d": present, "e": absent, "f": present}),
    }
    c = r["comparison"]
    assert (c["gene"], c["primary"]) == ("CCND1", "b")
    assert c["classes"]["c"] == unknown and c["classes"]["d"] == absent
    assert r["comparison_shares_a_cause_with_the_four"] == {
        "a": False,
        "b": True,
        "c": False,
        "d": False,
        "e": False,
        "f": True,
    }
    assert r["labels_or_verdicts_changed"] == 0 and r["bugs_fixed"] == []
    assert r["alphagenome_requests"] == 0


def test_every_sign_step_holds_and_the_verdict_reproduces() -> None:
    for t in [*_result()["four"], _result()["comparison"]]:
        steps = t["label_side"]["steps"] + t["measurement_side"]["steps"]
        assert [s["step"][:2] for s in steps] == ["L1", "L2", "L3", "L4", "M1", "M2", "M3", "M4"]
        assert all(s["consistent"] for s in steps), t["gene"]
        assert t["s4_rerun"]["activity"]["verdict"] == "incorrect"
        assert t["s4_rerun"]["target"]["verdict"] == "correct"
        assert t["s4_rerun"]["context"]["verdict"] == "not_judged"


def test_the_raw_rows_and_the_model_values_behind_each_class() -> None:
    by = {t["gene"]: t for t in [*_result()["four"], _result()["comparison"]]}
    rows = {g: t["measurement_side"]["deciding_row"] for g, t in by.items()}
    assert {g: (r["split"], r["line"]) for g, r in rows.items()} == {
        "CD83": ("heldout", 67),
        "HEMGN": ("training", 9323),
        "ID1": ("training", 6322),
        "BEX4": ("training", 10306),
        "CCND1": ("heldout", 104),
    }
    assert all(float(r["row"]["EffectSize"]) < 0 for g, r in rows.items() if g != "CCND1")
    assert float(rows["CCND1"]["row"]["EffectSize"]) > 0
    # the model's own value for the measured cell: agrees with the screen for three, not for ID1
    cells = {g: t["label_side"]["cached_gene_row"]["by_cell"] for g, t in by.items()}
    assert cells["CD83"]["GM12878"] < 0 and cells["HEMGN"]["K562"] < 0 and cells["BEX4"]["K562"] < 0
    assert cells["ID1"]["K562"] > 0
    assert "HCT116" not in cells["CCND1"]
    # the direction came from the single most extreme of 371 tracks, the rise beating the drop
    margins = {g: t["label_side"]["cached_gene_row"] for g, t in by.items()}
    for g in ("CD83", "HEMGN", "ID1", "BEX4"):
        m = margins[g]
        assert (
            m["max_drop_log2fc"] < 0 < m["max_rise_log2fc"] and m["max_rise_log2fc"] > -m["max_drop_log2fc"]
        )
    # indirect: MIR3193 is the any-gene prediction at ID1's element; ANP32B and TCEAL8 fell in the same screen
    assert by["ID1"]["label_side"]["run_row"]["predicted"]["gene"] == "MIR3193"
    fell = {
        g: {
            o["gene"]
            for o in by[g]["measurement_side"]["other_genes_same_screen"]
            if o["outcome"] == "significant_decrease"
        }
        for g in ("HEMGN", "BEX4", "CD83", "ID1")
    }
    assert fell == {"HEMGN": {"ANP32B"}, "BEX4": {"TCEAL8"}, "CD83": set(), "ID1": set()}


def test_the_recomputed_direction_is_the_cached_answers() -> None:
    from genomeos.predict import enhancer_target as et

    for t in [*_result()["four"], _result()["comparison"]]:
        g = t["label_side"]["cached_gene_row"]
        p = et.predict_target([g])
        assert p["action"] == t["label_side"]["run_row"]["predicted_coding"]["action"]
        assert p["tissue"] == t["label_side"]["run_row"]["predicted_coding"]["tissue"]


def test_the_judge_counts_and_what_refutable_means() -> None:
    j = _result()["judge"]
    assert j["committed_s4_counts_by_axis"] == {
        "activity": {"correct, refutable false": 93, "incorrect, refutable false": 5},
        "context": {"correct, refutable true": 39},
        "target": {"correct, refutable false": 59, "correct, refutable true": 39},
    }
    assert j["activity_verdicts_decided_in_the_stated_cell"] == {"correct": 39}
    assert j["activity_verdicts_decided_only_in_another_cell"] == {"correct": 54, "incorrect": 5}
    assert "not computed" in j["refutable_means"]


def test_the_judge_refutes_a_direction_from_another_cell_and_never_a_target() -> None:
    """The asymmetry the trace found, as the registered rules give it (no change is made)."""
    from genomeos.attribution import correctness as co

    def claim(axis: str, value: str) -> co.Claim:
        return co.Claim("E1", "chr1", 100, 400, axis, value, "G1", "placenta")

    rose_elsewhere = [co.Observation("crispri_decrease", "crispri:X", "G1", "K562", "training")]
    a = co.verdict_of(claim(co.ACTIVITY, "represses_target"), rose_elsewhere)
    assert (a.verdict, a.refutable) == (co.INCORRECT, False)
    null_elsewhere = [co.Observation("crispri_null_well_powered", "crispri:X", "G1", "K562", "training")]
    t = co.verdict_of(claim(co.TARGET, "G1"), null_elsewhere)
    assert (t.verdict, t.detail) == (co.NOT_JUDGED, co.NULL_ELSEWHERE)


@pytest.mark.skipif(
    not (ROOT / "data/knowledge/crispri").exists(), reason="the CRISPRi benchmark is local and untracked"
)
def test_the_id1_row_is_read_from_the_benchmark_as_committed(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(ROOT)
    rows = rt.benchmark_rows("chr20", 31608511, 31608860)
    id1 = [r for r in rows if r["row"]["measuredGeneSymbol"] == "ID1"]
    assert [(r["line"], r["row"]["EffectSize"], r["outcome"]) for r in id1] == [
        (6322, "-0.158249932905106", "significant_decrease")
    ]
