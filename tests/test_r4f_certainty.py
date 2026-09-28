# SPDX-License-Identifier: AGPL-3.0-or-later
"""R4f: the last magnitude-derived confidences are gone, and effect magnitude alone never raises
certainty (AlphaGenome variant effects, AlphaMissense scores, the compiler's enhancer-link record)."""

from __future__ import annotations

import random

import pytest

from genomeos.attribution.compile import element_certainty
from genomeos.genome import missense as am
from genomeos.ir.model import confidence_stated
from genomeos.predict import AlphaGenomeAdapter, PredictedEffect
from genomeos.predict.enhancer_target import link_certainty

RNG = random.Random(20260928)
EFFECTS = [0.0, 0.01, 0.05, 0.3, 0.7, 1.0, 5.0] + [RNG.uniform(-6, 6) for _ in range(300)]


def _without_magnitude(d: dict) -> dict:
    d = dict(d)
    d.pop("effect_estimate"), d.pop("model_score")
    return d


def test_predicted_effect_states_no_confidence():
    assert not hasattr(PredictedEffect("G", "liver", 0.8), "confidence")


@pytest.mark.parametrize("lfc", EFFECTS[:7])
def test_variant_effect_record_carries_effect_and_unit_and_no_probability(lfc):
    c = PredictedEffect("G", "liver", lfc).certainty
    assert c.effect_estimate == lfc and "log2 fold change" in c.effect_unit
    assert c.model_score == abs(lfc) and c.model_score_name
    assert c.probability is None and c.probability_unavailable


def test_variant_effect_magnitude_alone_never_raises_certainty():
    records = [_without_magnitude(PredictedEffect("G", "t", x).certainty.to_dict()) for x in EFFECTS]
    assert all(r == records[0] for r in records)


def test_adapter_rules_state_no_confidence_at_any_effect_size():
    scorer = lambda *_: [("G", f"t{i}", x) for i, x in enumerate(EFFECTS)]  # noqa: E731
    ad = AlphaGenomeAdapter(scorer=scorer)
    m = ad.to_module("chr21", 1, "A", "G", ad.predict("chr21", 1, "A", "G"))
    assert len(m.rules) == len(EFFECTS)
    assert {r.confidence for r in m.rules} == {0.0}
    assert not any(confidence_stated(r.confidence) for r in m.rules)


def test_missense_score_record_moves_only_in_its_score():
    records = [am.score_certainty(RNG.random()).to_dict() for _ in range(200)]
    for r in records:
        assert r["probability"] is None and r["probability_unavailable"]
        assert r["effect_estimate"] is None  # a pathogenicity score is not an effect size
        r.pop("model_score")
    assert all(r == records[0] for r in records)


def test_compiler_link_record_reads_both_link_forms_alike():
    for x in EFFECTS:
        lfc = round(x, 4)
        old = {"gene": "G", "log2_fold_change": lfc, "confidence": round(min(0.7, abs(lfc)), 3)}
        new = {"gene": "G", "log2_fold_change": lfc, "certainty": link_certainty(lfc).to_dict()}
        a, b = element_certainty(old), element_certainty(new)
        assert a == b and a.model_score == abs(lfc)  # the new form no longer loses its score
        assert a.probability is None


def test_compiler_link_record_magnitude_alone_never_raises_certainty():
    records = [
        _without_magnitude(element_certainty({"gene": "G", "log2_fold_change": round(x, 4)}).to_dict())
        for x in EFFECTS
    ]
    assert all(r == records[0] for r in records)
