import os

import pytest

from genomeos.ir import Action, EvidenceKind
from genomeos.predict import AlphaGenomeAdapter


def fake_scorer(chrom, pos, ref, alt):
    return [("APP", "brain", -0.8), ("APP", "liver", 0.1), ("CYYR1", "brain", 0.4)]


def test_adapter_emits_predicted_rules_only():
    ad = AlphaGenomeAdapter(scorer=fake_scorer)
    assert ad.available
    effects = ad.predict("chr21", 25897620, "C", "T")
    m = ad.to_module("chr21", 25897620, "C", "T", effects)
    assert len(m.rules) == 3
    assert all(r.evidence.kind is EvidenceKind.PREDICTED for r in m.rules)
    assert all(r.confidence <= 0.7 for r in m.rules)
    brain = next(r for r in m.rules if r.when == {"tissue": "brain"} and r.target == "APP")
    assert brain.action is Action.INHIBIT and brain.confidence > 0.35
    assert m.confidence_report()["rule"] < 0.7


@pytest.mark.skipif(
    not os.environ.get("ALPHAGENOME_API_KEY"), reason="set ALPHAGENOME_API_KEY for a live call"
)
def test_live_prediction():
    ad = AlphaGenomeAdapter()
    effects = ad.predict("chr21", 25897620, "C", "T")
    assert isinstance(effects, list)
