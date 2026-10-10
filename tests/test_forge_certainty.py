# SPDX-License-Identifier: AGPL-3.0-or-later
"""Review item R4, BioForge half: certainty is five separate things, and a probability needs a record.

The 2026-09-28 review found BioForge stamping a fixed 0.3 on every design answer, where a reader takes
it for a probability of correctness. These tests hold the replacement to the review's acceptance:

- evidence category, effect estimate, measurement uncertainty, model score and calibrated probability
  are separate fields;
- a probability cannot be constructed without its outcome, calibration population and evaluation
  method, and an unavailable probability stays unavailable (None, with the reason) rather than
  becoming a number;
- increasing the effect magnitude alone never increases the reported certainty.

They replace the line-number pins of tests/test_forge_calibration.py, which only preserved the old
constant.
"""

from __future__ import annotations

from dataclasses import replace

import pytest

from genomeos.certainty import Calibration, Certainty, Probability
from genomeos.forge.design_answers import results_for
from genomeos.ir import confidence_stated
from genomeos.lang import parse
from genomeos.organism.forge import UNSTATED_CONFIDENCE, Knob, design_certainty, run_designs

DESIGN_FILES = ("data/organisms/celegans/designs.bio", "data/organisms/human/haematopoiesis_designs.bio")


def _shipped() -> list:
    # one run per process, shared with tests/test_forge_design_answers.py and test_forge_calibration.py
    return [r for path in DESIGN_FILES for r in results_for(path)]


def test_a_probability_cannot_exist_without_its_calibration_record() -> None:
    record = Calibration(
        outcome="the returned knockout set equals the published one",
        population="design blocks with a held-out published answer, n = 40",
        method="Clopper-Pearson interval on the observed rate, one bin",
    )
    p = Probability(0.25, record)
    assert p.to_dict()["outcome"] and p.to_dict()["calibration_population"] and p.to_dict()["method"]
    for missing in ("outcome", "population", "method"):
        with pytest.raises(ValueError):
            replace(record, **{missing: " "})
    with pytest.raises(ValueError):
        Probability(1.3, record)
    with pytest.raises(TypeError):
        Probability(0.3)  # type: ignore[call-arg]  - a bare number is not a probability


def test_an_unavailable_probability_stays_unavailable_and_says_why() -> None:
    with pytest.raises(ValueError):
        Certainty(evidence_category="predicted", probability=None, probability_unavailable="")
    c = Certainty(evidence_category="predicted", probability_unavailable="no calibration record")
    d = c.to_dict()
    assert d["probability"] is None and d["probability_unavailable"] == "no calibration record"
    for key in (
        "evidence_category",
        "effect_estimate",
        "effect_unit",
        "measurement_uncertainty",
        "model_score",
    ):
        assert key in d


def test_bioforge_emits_no_probability_and_no_invented_confidence() -> None:
    """Since 2026-10-10 (area H, BIOFORGE-CONFIDENCE §6) the one number on the answer is the design
    author's stated confidence, carried unchanged and labelled as the author's; BioForge adds none."""
    results = _shipped()
    assert len(results) == 4
    for r in results:
        c = r.certainty()
        assert c is not None and c.probability is None and c.probability_unavailable
        assert c.evidence_category == "predicted"
        ex = r.to_experiment()
        assert ex is not None and confidence_stated(r.design.confidence)
        assert ex.confidence == r.design.confidence and ex.confidence != UNSTATED_CONFIDENCE
        assert (
            f"author of design {r.design.name}" in ex.evidence.note and "not BioForge's" in ex.evidence.note
        )
        bio = r.to_bio()
        assert f"confidence: {r.design.confidence:g}" in bio and bio.count("confidence: ") == 2  # note + key
        assert "0.3" not in bio  # no old constant
        assert "probability: unavailable" in bio
        d = r.to_dict()
        assert d["certainty"]["probability"] is None and d["certainty"]["evidence_category"] == "predicted"
        assert (
            d["author_confidence"] == r.design.confidence and "not BioForge's" in d["author_confidence_note"]
        )


def test_a_design_that_states_no_confidence_emits_none() -> None:
    """The author's number survives only where the author wrote one: no `confidence:` on the design, no
    `confidence:` key on the answer, UNSTATED on the object and null in the record."""
    src = """
module toy.design
import bio.std.development
organism T { root: R; cell_type: Zygote; resolution: populations }
timer cycle { duration: 60 min; when: cell_type = Zygote }
decision grow { action: divide; when: cell_type = Zygote; fraction: 1.0 }
design pace { vary: timer cycle duration 10..120; target: count at 300 min in 30..34; until: 300 min }
"""
    m = parse(src)
    r = run_designs(m, None, seed=1)[0]
    ex = r.to_experiment()
    assert ex is not None and not confidence_stated(ex.confidence) and ex.evidence.note == ""
    assert "confidence: " not in r.to_bio() and r.to_dict()["author_confidence"] is None
    assert r.certainty().probability is None


def test_effect_magnitude_alone_never_raises_the_reported_certainty() -> None:
    """Across the shipped designs the effect differs; the reported certainty does not move with it."""
    results = _shipped()
    effects = {r.certainty().effect_estimate for r in results}
    assert len(effects) > 1, "the shipped designs no longer differ in effect; this check is vacuous"
    # the probability is None for every effect, and the one number on the experiment is the author's,
    # read off the design block, so it equals what the author wrote whatever the effect
    assert {r.certainty().probability for r in results} == {None}
    assert all(r.to_experiment().confidence == r.design.confidence for r in results)
    # and directly: sweep the effect over nine orders of magnitude at a fixed answer loss
    sweep = [design_certainty(0.0, baseline) for baseline in (0.0, 1e-6, 0.1, 1.0, 10.0, 1e3)]
    assert [c.effect_estimate for c in sweep] == sorted(c.effect_estimate for c in sweep)
    assert {(c.probability, c.probability_unavailable, c.evidence_category) for c in sweep} == {
        (sweep[0].probability, sweep[0].probability_unavailable, sweep[0].evidence_category)
    }
    assert sweep[0].probability is None


def test_a_knob_the_search_moves_gets_a_ceiling_that_ignores_the_value() -> None:
    """The min(x, 0.3) sites cap an evidence-quality score; the cap never depends on how far a knob moved."""
    src = """
module toy.design
import bio.std.development
organism T { root: R; cell_type: Zygote; resolution: populations }
timer cycle { duration: 60 min; when: cell_type = Zygote; evidence: curated "toy"; confidence: 0.9 }
decision grow { action: divide; when: cell_type = Zygote; fraction: 1.0; confidence: 0.2 }
"""
    seen = set()
    for value in (10.0, 60.0, 120.0):
        m = parse(src)
        Knob.parse("timer cycle duration 10..120").set(m, value)
        Knob.parse("decision grow fraction 0..1").set(m, value / 120)
        seen.add((m.timer("cycle").confidence, next(d for d in m.decisions if d.id == "grow").confidence))
    assert seen == {(0.3, 0.2)}  # lowered to the ceiling, never raised, the same at every value
