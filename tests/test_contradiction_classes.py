# SPDX-License-Identifier: AGPL-3.0-or-later
"""The cause-class registration for direction v2's disagreeing rules (lane-contradict, 2026-10-02).

Every test here pins the registration: that the vocabulary is the one already registered and not a
copy of it, that exactly one class was added and is marked as added, that the imported primary
order is unchanged and the added class is never primary, and that class (c) cannot be reached on a
population with no measurement.
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


cc = _load("contradiction_classes")
rt = _load("repression_trace")


# ---- the vocabulary is the registered one, imported and not restated -----------------------------


def test_the_five_plus_f_come_from_the_repression_trace_and_are_not_copied() -> None:
    assert cc.INHERITED_CLASSES is cc.rt.CLASSES
    assert list(cc.INHERITED_CLASSES) == ["a", "b", "c", "d", "e", "f"]
    for key in ("a", "b", "c", "d", "e", "f"):
        assert cc.CLASSES[key] is cc.rt.CLASSES[key]
    assert cc.INHERITED_PRIMARY_RULE is cc.rt.PRIMARY_RULE


def test_the_class_letters_match_the_repression_trace_constants() -> None:
    assert (cc.SIGN, cc.CELL, cc.CONTRADICTION, cc.INDIRECT, cc.JUDGE, cc.SPLIT) == (
        rt.SIGN,
        rt.CELL,
        rt.CONTRADICTION,
        rt.INDIRECT,
        rt.JUDGE,
        rt.SPLIT,
    )
    assert (cc.PRESENT, cc.ABSENT, cc.UNKNOWN) == (rt.PRESENT, rt.ABSENT, rt.UNKNOWN)
    assert cc.WITHIN_CELL_SPLIT not in cc.INHERITED_CLASSES


def test_exactly_one_class_is_added_and_it_says_so() -> None:
    added = [k for k, c in cc.CLASSES.items() if c.get("added_for_this_population")]
    assert added == [cc.WITHIN_CELL_SPLIT]
    assert list(cc.CLASSES) == ["a", "b", "c", "d", "e", "f", "g"]
    g = cc.CLASSES[cc.WITHIN_CELL_SPLIT]
    assert g["name"] and g["definition"] and g["decided_by"] and g["why_added"]
    # it is a narrowing of (f), and it says which class it narrows
    assert "(f)" in g["definition"] and "(f)" in g["why_added"]


def test_no_inherited_class_was_edited_here() -> None:
    """`rt` is a second, independent load of the same file: equality proves nothing was mutated."""
    assert rt is not cc.rt
    for key, c in cc.INHERITED_CLASSES.items():
        assert "added_for_this_population" not in c, key
        assert c == rt.CLASSES[key], key
    assert cc.INHERITED_PRIMARY_RULE == rt.PRIMARY_RULE


# ---- the assignment rule -------------------------------------------------------------------------


def test_the_primary_order_is_the_inherited_one_and_g_is_never_primary() -> None:
    yes, no, unk = cc.PRESENT, cc.ABSENT, cc.UNKNOWN
    base = {"a": no, "b": yes, "c": unk, "d": no, "e": no, "f": no, "g": yes}
    assert cc.primary(base) == "b"
    assert cc.primary({**base, "c": yes}) == "c"
    assert cc.primary({**base, "c": yes, "e": yes}) == "e"
    assert cc.primary({**base, "a": yes, "c": yes, "e": yes}) == "a"
    # (g) present with nothing else is not a primary class: no class is primary, and that raises
    with pytest.raises(ValueError):
        cc.primary({"a": no, "b": no, "c": unk, "d": no, "e": no, "f": no, "g": yes})
    with pytest.raises(ValueError):
        cc.primary({"a": no, "b": no, "c": unk, "d": yes, "e": no, "f": yes, "g": yes})


def test_a_reading_that_is_not_present_absent_or_unknown_is_refused() -> None:
    base = {"a": cc.ABSENT, "b": cc.PRESENT, "c": cc.UNKNOWN, "d": cc.ABSENT, "e": cc.ABSENT, "f": cc.ABSENT}
    with pytest.raises(ValueError):
        cc.primary({**base, "g": "yes"})
    with pytest.raises(ValueError):
        cc.primary({**base, "g": True})


def test_contradiction_is_the_residual_class_and_needs_a_measurement() -> None:
    assert "only when no other class covers it" in cc.CONTRADICTION_IS_RESIDUAL
    assert "never the first reading" in cc.CONTRADICTION_IS_RESIDUAL
    # the clause is anchored in (c)'s own registered words, so it cannot drift from them
    c = cc.INHERITED_CLASSES["c"]
    assert "opposite sign to the measured effect" in c["definition"]
    assert "opposite in sign to the measurement" in c["decided_by"]
    assert "UNKNOWN by construction" in cc.CONTRADICTION_NEEDS_A_MEASUREMENT


def test_an_unevaluable_class_is_unknown_and_never_absent() -> None:
    assert "never recorded absent" in cc.ASSIGNMENT_RULE
    assert "absent is a reading and unknown is not" in cc.ASSIGNMENT_RULE


# ---- the population, named by the committed identifiers that produce it ---------------------------


def test_the_population_names_only_tallies_that_exist_in_committed_code() -> None:
    tallies = [p["tally"] for p in cc.POPULATION]
    assert tallies == ["signs_disagree_in_cell", "by_cell_value_opposite_to_the_selected_extreme"]
    assert "signs_disagree_in_cell" in dv.UNRESOLVED_REASONS
    counter = (ROOT / "scripts/direction_v2.py").read_text(encoding="utf-8")
    assert "by_cell_value_opposite_to_the_selected_extreme" in counter
    for p in cc.POPULATION:
        assert p["where"] and p["means"]


def test_no_count_is_restated_and_the_two_are_not_assumed_disjoint() -> None:
    reg = cc.registration()
    text = " ".join(str(v) for v in reg.values())
    for figure in ("14", "10", "24"):
        assert f" {figure} rules" not in text
    assert "sum is NOT" in cc.POPULATION_NOT_ASSUMED_DISJOINT
    assert "never by adding" in cc.POPULATION_NOT_ASSUMED_DISJOINT
    assert "records no element identifier" in cc.NO_ENUMERATION_ASSUMED


def test_the_diagnostic_and_the_no_reversal_clause_are_carried_word_for_word() -> None:
    assert cc.DIAGNOSTIC_READING is dv.SELECTION_EXCLUDED_DIAGNOSTIC
    assert cc.V2_CANNOT_REVERSE is dv.FLIPS_ARE_STRUCTURALLY_IMPOSSIBLE
    assert "never as a count of wrong rules" in cc.DIAGNOSTIC_READING
    assert "can therefore withhold a direction and can never reverse one" in cc.V2_CANNOT_REVERSE


def test_g_decided_by_is_the_condition_v2_already_codes() -> None:
    """(g)'s test must be v2's own `signs_disagree_in_cell`, not a second reading of it."""
    assert et._cell_summary([-1.0, 0.5])["signs_disagree"] is True
    assert et._cell_summary([-1.0, -0.5])["signs_disagree"] is False
    row = {"max_drop_tissue": "K562", "max_drop_log2fc": -1.0, "by_cell": {"K562": 0.5}}
    call = dv.direction_call({"tissue": "K562", "log2_fold_change": -1.0}, dv.V2, row=row)
    assert call.resolved is False and call.reason == "signs_disagree_in_cell"
    g = cc.CLASSES[cc.WITHIN_CELL_SPLIT]
    assert "signs_disagree_in_cell" in g["decided_by"]
    assert "retained_values" in g["decided_by"]


# ---- what this registration may and may not be called --------------------------------------------


def test_the_registration_claims_no_validation_and_changes_nothing() -> None:
    reg = cc.registration()
    assert reg["registered"] == "2026-10-02"
    assert "not a validation" in reg["status"]
    assert "no rule, label, direction or verdict is changed" in cc.NO_CHANGE
    assert "separate decision" in cc.NO_CHANGE
    assert any("validation" in s for s in cc.MAY_NOT_BE_CALLED)
    assert any("no reading here licences a reversal" in s for s in cc.MAY_NOT_BE_CALLED)
    assert "NOT established biological independence" in cc.CELL_GROUPING_CAVEAT


def test_main_classifies_nothing_and_says_so(capsys: pytest.CaptureFixture[str]) -> None:
    assert cc.main([]) == 0
    out = capsys.readouterr().out
    assert "no rule has been classified" in out
    assert "the registration and nothing else" in out
    assert cc.main(["--json"]) == 0
