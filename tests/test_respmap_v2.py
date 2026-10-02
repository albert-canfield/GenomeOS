# SPDX-License-Identifier: AGPL-3.0-or-later
"""The per-assertion direction-v2 emitter: the class comes from the rule, and the rule is unchanged.

Every test here is on the classifier and on the wording it must carry. None of them reads the
committed result, so a test passes or fails on the code and not on a number.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genomeos.attribution import direction_v2 as dv
from genomeos.attribution import respmap_v2 as rv
from genomeos.predict import enhancer_target as et

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/results/response_map_increment2.json"


def assertion(
    *,
    element: str = "EH38E0000001",
    chrom: str = "chr21",
    gene: str = "GENE1",
    cell: str = "K562",
    action: str = "activates",
    value: float = -0.8,
) -> dict:
    """One response-map assertion, shaped as `response_map2._compiled_rule` writes it."""
    pc = {
        "gene": gene,
        "action": action,
        "log2_fold_change": value,
        "tissue": cell,
        "basis": "predicted: expression change on deleting the element",
    }
    return {
        "id": f"r2|{chrom}:100-200|{element}|{gene}",
        "status": "predicted",
        "measurement": {"direction": action, "value": value, "is_a_measurement": False},
        "quoted": {"id": element, "predicted_coding": pc},
    }


def row(**fields) -> dict:
    """A cached gene row, with only the fields `direction_v2.retained_values` reads."""
    base = {
        "gene": "GENE1",
        "max_drop_tissue": "",
        "max_drop_log2fc": 0.0,
        "max_rise_tissue": "",
        "max_rise_log2fc": 0.0,
        "by_cell": {},
    }
    base.update(fields)
    return base


# ---- the premise this lane had to check before counting anything ------------------------------


def test_direction_v2_script_writes_no_file() -> None:
    """The premise check, held as a test: the counting script records no per-assertion class."""
    src = (ROOT / "scripts/direction_v2.py").read_text()
    for forbidden in ("json.dump", "write_text", "save_result", "open("):
        assert forbidden not in src, f"{forbidden} now appears; the emitter may be redundant"


def test_no_committed_result_carried_a_v2_class() -> None:
    """Had any result carried a reason token, the class could have been read instead of emitted."""
    results = ROOT / "data/results"
    if not results.is_dir():
        pytest.skip("no result registry on this machine")
    hits = [
        p.name
        for p in results.glob("*.json")
        if p.name != "respmap_direction_v2.json" and "one_track_seen_twice" in p.read_text(errors="ignore")
    ]
    assert hits == [], f"a result now carries v2 reasons: {hits}"


# ---- the class comes from the rule -------------------------------------------------------------


def test_unanimous_two_values_resolve_and_agree() -> None:
    a = assertion(action="activates", value=-0.8)
    r = row(max_drop_tissue="K562", max_drop_log2fc=-0.8, by_cell={"K562": -0.5})
    out = rv.classify(a, row=r)
    assert out["v2_class"] == "resolved_agrees_with_published"
    assert out["v2_action"] == dv.ACTIVATES
    assert out["v2_unresolved_reason"] is None
    assert out["published_value_is_among_retained"] is True


def test_disagreeing_signs_are_unresolved_with_their_reason() -> None:
    a = assertion(value=-0.8)
    r = row(max_drop_tissue="K562", max_drop_log2fc=-0.8, by_cell={"K562": 0.6})
    out = rv.classify(a, row=r)
    assert out["v2_class"] == "unresolved"
    assert out["v2_unresolved_reason"] == "signs_disagree_in_cell"
    assert out["v2_action"] is None


def test_one_value_only_is_unresolved_and_absent_row_names_its_own_reason() -> None:
    a = assertion(value=-0.8)
    one = rv.classify(a, row=row(max_drop_tissue="K562", max_drop_log2fc=-0.8))
    assert one["v2_class"] == "unresolved"
    assert one["v2_unresolved_reason"] == "one_value_only"
    none = rv.classify(a, row=None)
    assert none["v2_class"] == "unresolved"
    assert none["v2_unresolved_reason"] == "no_row_for_target"
    assert none["cached_row_found"] is False


def test_two_fields_one_track_is_unresolved_under_amendment_1() -> None:
    a = assertion(value=-0.8)
    r = row(max_drop_tissue="K562", max_drop_log2fc=-0.8, by_cell={"K562": -0.8})
    assert rv.classify(a, row=r)["v2_unresolved_reason"] == "one_track_seen_twice"


def test_below_the_imported_floor_is_unresolved_and_the_floor_is_not_restated() -> None:
    small = dv.MAGNITUDE_FLOOR / 10
    a = assertion(value=-small)
    r = row(max_drop_tissue="K562", max_drop_log2fc=-small, by_cell={"K562": -small / 2})
    assert rv.classify(a, row=r)["v2_unresolved_reason"] == "below_magnitude_floor"
    assert dv.MAGNITUDE_FLOOR == et.MIN_EFFECT


def test_absence_gloss_correction_is_counted_and_not_glossed() -> None:
    """The escape branch: the one retained value is by_cell, not the published value, opposite sign."""
    a = assertion(value=-0.8)
    r = row(by_cell={"K562": 0.6})
    out = rv.classify(a, row=r)
    assert out["v2_unresolved_reason"] == "one_value_only"
    assert out["published_value_is_among_retained"] is False
    assert out["only_retained_is_by_cell_opposite_to_published"] is True
    # and the same-sign case is not counted as the branch
    same = rv.classify(a, row=row(by_cell={"K562": -0.6}))
    assert same["only_retained_is_by_cell_opposite_to_published"] is False


def test_a_flip_would_be_recorded_as_a_flip_if_the_rule_ever_produced_one() -> None:
    """The flip class is reachable in code, so a measured 0 is a measurement and not an impossibility.

    This plants a row whose retained values are unanimous and opposite to the published action,
    which the registered reading says the committed cache cannot present; the point is only that
    the classifier would not swallow it.
    """
    a = assertion(action="activates", value=-0.8)
    r = row(max_drop_tissue="K562", max_drop_log2fc=0.9, by_cell={"K562": 0.4})
    out = rv.classify(a, row=r)
    assert out["v2_class"] == "resolved_opposite_to_published"
    assert out["v2_action"] == dv.INHIBITS
    assert out["published_value_is_among_retained"] is False


def test_represses_and_inhibits_are_one_axis() -> None:
    a = assertion(action="represses", value=0.8)
    r = row(max_rise_tissue="K562", max_rise_log2fc=0.8, by_cell={"K562": 0.4})
    out = rv.classify(a, row=r)
    assert out["published_axis"] == "represses_target"
    assert out["v2_action"] == dv.INHIBITS
    assert out["v2_class"] == "resolved_agrees_with_published"


# ---- the shape of the result -------------------------------------------------------------------


def test_classes_are_exhaustive_and_counts_add_to_the_denominator() -> None:
    rows = [
        rv.classify(
            assertion(), row=row(max_drop_tissue="K562", max_drop_log2fc=-0.8, by_cell={"K562": -0.4})
        ),
        rv.classify(
            assertion(), row=row(max_drop_tissue="K562", max_drop_log2fc=-0.8, by_cell={"K562": 0.4})
        ),
        rv.classify(assertion(), row=None),
    ]
    c = rv.counts(rows)
    assert sum(c["by_class"].values()) == c["denominator"] == 3
    assert set(c["by_class"]) == set(rv.CLASSES)
    assert set(c["unresolved_by_reason"]) <= set(dv.UNRESOLVED_REASONS)
    absent = set(c["unresolved_reasons_with_no_assertion"])
    assert absent | set(c["unresolved_by_reason"]) == set(dv.UNRESOLVED_REASONS)


def test_a_reshaped_id_fails_loudly() -> None:
    a = assertion()
    a["id"] = "chr21:100-200"
    with pytest.raises(ValueError):
        rv.locus_of(a)
    b = assertion(element="EH38E0000001")
    b["quoted"]["id"] = "EH38E0000002"
    with pytest.raises(ValueError):
        rv.locus_of(b)


def test_only_predicted_assertions_are_in_the_population() -> None:
    payload = {
        "assertions": [
            assertion(),
            {"id": "x", "status": "observed"},
            {"id": "y", "status": "inferred"},
        ]
    }
    got = rv.predicted_assertions(payload)
    assert [a["status"] for a in got] == ["predicted"]


def test_the_published_payload_fixes_the_denominator() -> None:
    if not SOURCE.exists():
        pytest.skip("the published response map is not on this machine")
    payload = json.loads(SOURCE.read_text())
    assert payload["counts"]["by_status"]["predicted"] == 166
    assert payload["counts"]["by_status"]["observed"] == 361
    assert payload["counts"]["assertions"] == 527
    assert len(rv.predicted_assertions(payload)) == 166


def test_increment_3_is_a_different_population() -> None:
    """The brief attached increment 2's counts to increment 3's name; the code records which is which."""
    p = ROOT / "data/results/response_map_increment3.json"
    if not p.exists():
        pytest.skip("increment 3 is not on this machine")
    counts = json.loads(p.read_text())["counts"]
    assert counts["assertions"] == 574
    assert counts["by_status"]["predicted"] == 73
    assert counts["chains"] == 93
    assert "574 assertions over 337 entities in 93 chains" in rv.INCREMENT_3_IS_NOT_THIS_POPULATION


# ---- the readings this lane may not weaken -----------------------------------------------------


def test_the_registered_readings_are_carried_and_not_restated() -> None:
    assert dv.FLIPS_ARE_STRUCTURALLY_IMPOSSIBLE in rv.FLIP_IS_REPORTED_AS_MEASURED
    assert rv.DIAGNOSTIC_IS_NOT_A_V2_OUTPUT == dv.SELECTION_EXCLUDED_DIAGNOSTIC
    assert rv.UNRESOLVED_MEANS == dv.UNRESOLVED_MEANS
    assert dv.COUNTS_VALIDATE_NOTHING in rv.VALIDATES_NOTHING
    assert dv.NOT_VALIDATION in rv.VALIDATES_NOTHING


def test_the_too_strong_absence_gloss_is_not_repeated() -> None:
    text = " ".join(v for v in vars(rv).values() if isinstance(v, str)) + " ".join(
        s for t in vars(rv).values() if isinstance(t, tuple) for s in t if isinstance(s, str)
    )
    assert "may NOT be glossed" in rv.ABSENCE_GLOSS_CORRECTED
    assert "one value retained for that cell" not in text
    assert "388,997" not in text


def test_the_two_overlapping_tallies_are_not_added_anywhere() -> None:
    assert "may not be added" in rv.EARLIER_TALLIES_NOT_TOUCHED
    assert "14 `signs_disagree_in_cell`" in rv.EARLIER_TALLIES_NOT_TOUCHED
    # no module constant performs arithmetic on them
    src = (ROOT / "genomeos/attribution/respmap_v2.py").read_text()
    assert "14 +" not in src and "+ 10" not in src


def test_the_lane_makes_no_recommendation_about_the_grammar_token() -> None:
    assert "makes no recommendation" in rv.NO_RECOMMENDATION
    for word in ("should", "recommend that", "we recommend"):
        assert word not in rv.NO_RECOMMENDATION.replace("makes no recommendation", "")


def test_nothing_is_rescored_and_no_threshold_is_defined_here() -> None:
    src = (ROOT / "genomeos/attribution/respmap_v2.py").read_text()
    assert "MIN_EFFECT =" not in src
    assert "MAGNITUDE_FLOOR =" not in src
    assert "no threshold of its own" in rv.NOTHING_RESCORED


def test_grouping_is_not_independence() -> None:
    joined = " ".join(rv.CANNOT_ESTABLISH)
    assert "NOT established biological independence" in joined
