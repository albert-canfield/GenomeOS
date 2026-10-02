# SPDX-License-Identifier: AGPL-3.0-or-later
"""The reading of direction v2's disagreeing rules (lane-contradict, 2026-10-02).

Each branch of the derivation is planted on a row written here, so the claim is checked against
`direction_v2` as committed rather than argued in prose. Nothing is read from disk but the two
scripts' own source, and no cache is opened.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomeos.attribution import direction_v2 as dv  # noqa: E402
from genomeos.predict import enhancer_target as et  # noqa: E402


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / f"scripts/{name}.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


cr = _load("contradiction_reading")
cc = cr.cc

#: the diagnostic's condition, as scripts/direction_v2.chromosome_counts codes it inline. The test
#: below asserts those exact fragments are in that file, so this cannot drift from it.
_DIAGNOSTIC_SOURCE = (
    'other = [v for src, v in dv.retained_values(row, v1.cell) if src == "by_cell"]',
    'selected = float(pc["log2_fold_change"])',
    "if other and other[0] != 0.0 and (other[0] > 0) != (selected > 0):",
)


def _diagnostic_fires(pc: dict[str, object], row: dict[str, object]) -> bool:
    v1 = dv.direction_call(pc, dv.V1)
    other = [v for src, v in dv.retained_values(row, v1.cell) if src == "by_cell"]
    selected = float(pc["log2_fold_change"])  # type: ignore[arg-type]
    return bool(other and other[0] != 0.0 and (other[0] > 0) != (selected > 0))


def test_the_diagnostic_condition_is_still_coded_the_way_the_derivation_reads_it() -> None:
    src = (ROOT / "scripts/direction_v2.py").read_text(encoding="utf-8")
    for fragment in _DIAGNOSTIC_SOURCE:
        assert fragment in src, fragment


# ---- the derivation, branch by branch ------------------------------------------------------------


def test_the_ordinary_branch_a_diagnostic_hit_is_also_signs_disagree_in_cell() -> None:
    """SUBSET_DERIVATION's main claim: the two tallies overlap, so their sum counts rules twice."""
    pc = {"tissue": "K562", "log2_fold_change": -1.2, "action": "activates"}
    row = {"max_drop_tissue": "K562", "max_drop_log2fc": -1.2, "by_cell": {"K562": 0.4}}
    assert _diagnostic_fires(pc, row) is True
    call = dv.direction_call(pc, dv.V2, row=row)
    assert call.resolved is False and call.reason == "signs_disagree_in_cell"
    assert "the two tallies overlap" in cr.SUBSET_DERIVATION
    assert "their sum is not a count of distinct rules" in cr.SUBSET_DERIVATION


def test_escape_branch_one_a_zero_third_value_reaches_zero_value_in_cell_first() -> None:
    pc = {"tissue": "K562", "log2_fold_change": -1.2, "action": "activates"}
    row = {
        "max_drop_tissue": "K562",
        "max_drop_log2fc": -1.2,
        "max_rise_tissue": "K562",
        "max_rise_log2fc": 0.0,
        "by_cell": {"K562": 0.4},
    }
    assert _diagnostic_fires(pc, row) is True
    assert dv.direction_call(pc, dv.V2, row=row).reason == "zero_value_in_cell"
    assert "zero_value_in_cell" in cr.SUBSET_DERIVATION


def test_escape_branch_two_a_row_naming_another_tissue_reaches_one_value_only() -> None:
    """The branch that makes the ROADMAP row's one_value_only gloss too strong."""
    pc = {"tissue": "K562", "log2_fold_change": -1.2, "action": "activates"}
    row = {
        "max_drop_tissue": "HepG2",
        "max_drop_log2fc": -1.2,
        "max_rise_tissue": "HepG2",
        "max_rise_log2fc": 0.3,
        "by_cell": {"K562": 0.4},
    }
    assert _diagnostic_fires(pc, row) is True
    call = dv.direction_call(pc, dv.V2, row=row)
    assert call.reason == "one_value_only"
    assert call.values == (("by_cell", 0.4),)  # the one value is NOT the selected extreme
    assert "is by_cell[cell] and not the selected extreme" in cr.ONE_VALUE_ONLY_GLOSS_IS_TOO_STRONG


def test_signs_disagree_in_cell_can_occur_with_no_diagnostic_hit() -> None:
    """So the two tallies are not equal either: neither contains the other by force."""
    pc = {"tissue": "K562", "log2_fold_change": 0.9, "action": "represses"}
    row = {
        "max_drop_tissue": "K562",
        "max_drop_log2fc": -0.5,
        "max_rise_tissue": "K562",
        "max_rise_log2fc": 0.9,
    }
    assert _diagnostic_fires(pc, row) is False
    assert dv.direction_call(pc, dv.V2, row=row).reason == "signs_disagree_in_cell"
    assert "not an identity either" in cr.NEITHER_TALLY_CONTAINS_THE_OTHER


def test_predict_target_takes_the_extreme_and_its_tissue_from_one_row() -> None:
    """The step SUBSET_DERIVATION rests on, read from enhancer_target rather than asserted."""
    rows = [
        {
            "gene": "G",
            "max_drop_log2fc": -1.5,
            "max_drop_tissue": "K562",
            "max_rise_log2fc": 0.2,
            "max_rise_tissue": "HepG2",
        }
    ]
    pc = et.predict_target(rows)
    assert pc is not None
    assert pc["log2_fold_change"] == -1.5 and pc["tissue"] == "K562"
    assert dv.direction_call(pc, dv.V1).cell == "K562"


# ---- what could not be read ----------------------------------------------------------------------


def test_the_counting_script_writes_no_file_and_records_no_identifier() -> None:
    src = (ROOT / "scripts/direction_v2.py").read_text(encoding="utf-8")
    assert "json.dump" not in src and "open(" not in src and "write_text" not in src
    # the tallies are Counter keys; no element id, locus or gene is put in either Counter
    assert "tally[" in src and "reasons[" in src
    assert 'tally[f"unresolved_from_{axis}"]' in src
    assert "writes no file at all" in cr.ENUMERATION_ABSENT
    assert "the cache read this lane is forbidden" in cr.ENUMERATION_ABSENT


def test_the_reading_names_no_rule_and_says_why() -> None:
    r = cr.reading()
    assert r["rules_named"] == []
    assert r["rules_named_why_empty"] is cr.ENUMERATION_ABSENT
    assert r["registration_committed_at"] == "c6183c2"
    assert "no rule read, no class assigned" in r["read_after"]


# ---- the class assignment ------------------------------------------------------------------------


def test_the_reading_covers_every_registered_class_exactly_once() -> None:
    assert set(cr.POPULATION_CLASSES) == set(cc.CLASSES)
    for key, (r, why) in cr.POPULATION_CLASSES.items():
        assert r in (cc.PRESENT, cc.ABSENT, cc.UNKNOWN), key
        assert why


def test_every_rankable_class_is_unknown_and_only_g_is_present() -> None:
    present = [k for k, (r, _) in cr.POPULATION_CLASSES.items() if r == cc.PRESENT]
    assert present == [cc.WITHIN_CELL_SPLIT]
    for key in ("a", "b", "c", "d", "e", "f"):
        assert cr.POPULATION_CLASSES[key][0] == cc.UNKNOWN, key
    assert not [k for k, (r, _) in cr.POPULATION_CLASSES.items() if r == cc.ABSENT]


def test_the_primary_rule_cannot_name_a_class_for_this_population() -> None:
    readings = {k: r for k, (r, _) in cr.POPULATION_CLASSES.items()}
    with pytest.raises(ValueError):
        cc.primary(readings)
    assert "names no primary class for this population" in cr.PRIMARY_CLASS_UNREACHABLE
    assert "a model disagreeing with itself" in cr.PRIMARY_CLASS_UNREACHABLE


def test_g_is_present_because_it_is_the_same_condition_not_a_second_reading() -> None:
    why = cr.POPULATION_CLASSES[cc.WITHIN_CELL_SPLIT][1]
    assert "IS _v2's condition for that reason" in why
    assert "the same fact" in why
    assert "signs_disagree_in_cell" in cc.CLASSES[cc.WITHIN_CELL_SPLIT]["decided_by"]


# ---- the verdict, and what it may not be called ---------------------------------------------------


def test_no_genuine_contradiction_survives_and_the_reason_is_not_overstated() -> None:
    g = cr.GENUINE_CONTRADICTIONS
    assert "(c) was never reachable" in g
    assert "NOT a finding that the disagreements are explained away" in g
    assert "NOT a finding that no contradiction exists in the compiled genome" in g
    assert "cannot be answered in the vocabulary it asked for" in g
    assert any("clean bill of health" in s for s in cr.MAY_NOT_BE_CALLED)
    assert any("not because candidates were eliminated by evidence" in s for s in cr.MAY_NOT_BE_CALLED)


def test_both_available_widenings_are_refused() -> None:
    assert "both are refused" in cr.NO_WIDENING
    assert "(b) was not widened" in cr.NO_WIDENING
    assert "NOT (b)" in cr.POPULATION_CLASSES[cc.CELL][1]


def test_nothing_is_edited_or_proposed_for_reversal() -> None:
    assert "no rule is edited, proposed for edit or named for reversal here" in cr.WHAT_WOULD_SETTLE_IT
    assert "neither of them this lane's to do" in cr.WHAT_WOULD_SETTLE_IT


def test_the_registration_is_imported_unchanged() -> None:
    assert cr.cc.NAME == "contradiction_classes"
    assert list(cc.CLASSES) == ["a", "b", "c", "d", "e", "f", "g"]
    assert cc.INHERITED_CLASSES is cc.rt.CLASSES


def test_main_names_no_rule_and_says_why(capsys: pytest.CaptureFixture[str]) -> None:
    assert cr.main([]) == 0
    out = capsys.readouterr().out
    assert "no rule is named, because no committed artefact names one" in out
    assert "no rule is changed" in out
    assert cr.main(["--json"]) == 0
