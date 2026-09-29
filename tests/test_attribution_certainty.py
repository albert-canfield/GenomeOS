# SPDX-License-Identifier: AGPL-3.0-or-later
"""Review item R4b, attribution half outside the compiler: effect magnitude alone never raises certainty.

The census is genomeos/attribution/confidence_census.py. Two producers derived a stated confidence
from how far a measured quantity lies past a bar: human_panel.block_class (all four called classes)
and variation.case_of. budget.guess and candidates hand-set one constant per rule or per evidence
layer. These tests sweep each effect with everything else fixed and require that nothing reading as
certainty moves: not a `confidence`, not the evidence category, not a model score, not a
probability. The effect itself may, and should, move in `certainty.effect_estimate`.

A probability, where one appears, names its outcome, calibration population and method; where none
exists it is None with its reason.

Committed first as strict xfails against the old code (ef31028). The build makes every one of them
pass; the markers stay as that record and are lifted at the end of this file.
"""

from __future__ import annotations

import pytest

from genomeos.attribution import human_panel as hp
from genomeos.attribution.budget import guess
from genomeos.attribution.confidence_census import CENSUS
from genomeos.attribution.variation import case_of

BASE = {"aligned_share": 1.0, "missing_share": 0.0, "presence": 0.999, "touched_recurring_share": 0.01}


def reported(out: dict) -> tuple:
    """Everything in a result that a reader could take as certainty, as one comparable value."""
    c = out.get("certainty") or {}
    p = c.get("probability")
    return (
        out.get("confidence"),
        c.get("evidence_category"),
        None if p is None else p["value"],
        c.get("model_score"),
    )


def assert_certainty_record(out: dict) -> None:
    c = out["certainty"]
    assert c["evidence_category"]
    if c["probability"] is None:
        assert c["probability_unavailable"], "an unavailable probability says why"
    else:
        for k in ("outcome", "calibration_population", "method"):
            assert c["probability"][k]
    if c["effect_estimate"] is not None:
        assert c["effect_unit"]


def test_the_census_names_every_magnitude_site() -> None:
    magnitude = [e["site"] for e in CENSUS if e["kind"] == "magnitude"]
    assert len(magnitude) == 5
    assert all("human_panel.py" in s or "variation.py" in s for s in magnitude)
    assert all(e["consumers"] for e in CENSUS)


@pytest.mark.xfail(strict=True, reason="R4b census: block_class confidence rises with the absence")
def test_lineage_restricted_certainty_does_not_rise_as_presence_falls() -> None:
    outs = [hp.block_class({**BASE, "presence": p}, 18, 20.0) for p in (0.49, 0.4, 0.25, 0.1, 0.0)]
    assert {o["class"] for o in outs} == {"lineage_restricted"}
    assert len({reported(o) for o in outs}) == 1


@pytest.mark.xfail(strict=True, reason="R4b census: 0.3 + min(0.3, gap * 3)")
def test_polymorphic_certainty_does_not_rise_with_the_gap() -> None:
    outs = [hp.block_class({**BASE, "presence": p}, 18, 20.0) for p in (0.94, 0.9, 0.8, 0.6, 0.51)]
    outs += [hp.block_class({**BASE, "touched_recurring_share": t}, 18, 20.0) for t in (0.11, 0.2, 0.5)]
    assert {o["class"] for o in outs} == {"polymorphic"}
    assert len({reported(o) for o in outs}) == 1


@pytest.mark.xfail(strict=True, reason="R4b census: core confidence rises with the depletion")
def test_core_certainty_does_not_rise_with_the_depletion() -> None:
    outs = [hp.block_class(BASE, k, 40.0) for k in (19, 15, 10, 5, 2, 0)]
    assert {o["class"] for o in outs} == {"core"}
    assert len({reported(o) for o in outs}) == 1


@pytest.mark.xfail(strict=True, reason="R4b census: variable confidence rises with |ratio - 0.5|")
def test_variable_certainty_does_not_rise_with_the_distance_from_the_bar() -> None:
    outs = [hp.block_class(BASE, k, 20.0) for k in (12, 18, 25, 40, 80)]
    assert {o["class"] for o in outs} == {"variable"}
    assert len({reported(o) for o in outs}) == 1


@pytest.mark.xfail(strict=True, reason="R4b census: block_class carries no certainty record yet")
def test_every_panel_class_carries_a_certainty_record_with_its_effect() -> None:
    cases = [
        (BASE, 2, 20.0),
        (BASE, 18, 20.0),
        (BASE, 0, 2.0),
        ({**BASE, "presence": 0.9}, 18, 20.0),
        ({**BASE, "presence": 0.3}, 18, 20.0),
        ({**BASE, "aligned_share": 0.2}, 18, 20.0),
    ]
    for m, obs, exp in cases:
        out = hp.block_class(m, obs, exp)
        assert "confidence" not in out, "a magnitude-derived confidence is not reported"
        assert_certainty_record(out)
        if out["class"] != "unplaced":
            assert out["certainty"]["effect_estimate"] is not None


@pytest.mark.xfail(strict=True, reason="R4b census: case_of confidence rises with both distances")
def test_variation_case_certainty_does_not_rise_with_the_distance_from_the_bars() -> None:
    outs = [case_of(m, h) for m, h in ((0.06, 0.26), (0.1, 0.4), (0.3, 0.7), (1.0, 1.0))]
    assert {o["case"] for o in outs} == {"syntax"}
    assert len({reported(o) for o in outs}) == 1
    for o in outs:
        assert "confidence" not in o
        assert_certainty_record(o)


BUDGET_LABELS = [
    ("interspersed_repeat_LINE", (0.0, 0.01, 0.029)),  # fossil, unconstrained
    ("interspersed_repeat_LINE", (0.03, 0.04, 0.049)),  # fossil, weak constraint
    ("interspersed_repeat_LINE", (0.05, 0.2, 0.9)),  # constrained_unknown, exapted
    ("regulatory", (0.05, 0.3, 1.0)),
    ("regulatory", (0.0, 0.02, 0.049)),
    ("long_orf", (0.05, 0.5)),
    ("unique_intergenic", (0.05, 0.1, 0.8)),
    ("unique_intergenic", (0.0, 0.01, 0.029)),
    ("unique_intergenic", (0.03, 0.045)),
]


@pytest.mark.parametrize(("cls", "fractions"), BUDGET_LABELS)
def test_budget_score_does_not_move_with_the_constrained_fraction_inside_one_label(cls, fractions) -> None:
    outs = [guess(cls, 0.6, {"fraction_above": fa}, {"n": 0}) for fa in fractions]
    assert len({(o["tier"], o["label"]) for o in outs}) == 1
    assert len({o["confidence"] for o in outs}) == 1


@pytest.mark.xfail(strict=True, reason="R4b census: budget.guess carries no certainty record yet")
def test_budget_guess_carries_the_constrained_fraction_as_its_effect_and_no_probability() -> None:
    measured = guess("unique_intergenic", 0.3, {"fraction_above": 0.07}, {"n": 0})
    assert_certainty_record(measured)
    assert measured["certainty"]["effect_estimate"] == 0.07
    assert measured["certainty"]["probability"] is None
    unmeasured = guess("unique_intergenic", 0.3, None, None)
    assert_certainty_record(unmeasured)
    assert unmeasured["certainty"]["effect_estimate"] is None
    assert unmeasured["confidence"] <= 0.3


# The build after ef31028 makes each strict xfail above pass, so the tests run as plain tests now. The
# markers are kept as the record of the census commit and lifted here rather than deleted.
for _test in (
    test_lineage_restricted_certainty_does_not_rise_as_presence_falls,
    test_polymorphic_certainty_does_not_rise_with_the_gap,
    test_core_certainty_does_not_rise_with_the_depletion,
    test_variable_certainty_does_not_rise_with_the_distance_from_the_bar,
    test_every_panel_class_carries_a_certainty_record_with_its_effect,
    test_variation_case_certainty_does_not_rise_with_the_distance_from_the_bars,
    test_budget_guess_carries_the_constrained_fraction_as_its_effect_and_no_probability,
):
    _test.pytestmark = [m for m in getattr(_test, "pytestmark", []) if m.name != "xfail"]


def test_no_xfail_is_left_on_the_certainty_properties() -> None:
    marked = [
        n
        for n, f in globals().items()
        if n.startswith("test_") and any(m.name == "xfail" for m in getattr(f, "pytestmark", []))
    ]
    assert marked == []
