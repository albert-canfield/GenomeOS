# SPDX-License-Identifier: AGPL-3.0-or-later
"""Review R4, compiler half (2026-09-28): what the compiled program states about certainty.

Before this, `attribution/compile.py` wrote `min(0.7, max(0.05, pc.confidence or |log2 fc|))` as the
`confidence:` of every predicted element and rule, so the size of a predicted effect was read as how
sure the program is; every gene stub said 0.9 and every domain 0.4 whatever supported them. These are
properties, not pinned numbers: no effect size may move a stated certainty, and any probability the
program states must name its outcome, calibration population and method.
"""

from __future__ import annotations

import json
import re

import pytest

from genomeos.attribution.compile import compile_chromosome, element_certainty
from genomeos.attribution.measured import Layer
from genomeos.lang.parser import parse

EFFECTS = (0.01, 0.05, 0.2, 0.5, 0.69, 0.9, 1.5, 5.0)


def _fixture(tmp_path, effects):
    block = {
        "start": 0,
        "end": 1000,
        "length": 1000,
        "class": "regulatory",
        "phylop": {"bases": 1000, "fraction_above": 0.2, "above": 200},
        "elements": {"n": 2},
        "guess": {"tier": "regulatory", "label": "regulatory elements", "confidence": 0.7},
    }
    (tmp_path / "budget_chrT.json").write_text(
        json.dumps({"chrom": "chrT", "unknown_bp": 1000, "constrained_fraction": 0.2, "blocks": [block]})
    )
    (tmp_path / "domains_chrT.json").write_text(
        json.dumps({"domains": [{"id": "chrT:D1", "start": 0, "end": 100000, "confidence": 0.4}]})
    )
    els = []
    for i, lfc in enumerate(effects):
        for sign in (1, -1):
            v = sign * lfc
            els.append(
                {
                    "id": f"E{i}_{'p' if sign > 0 else 'm'}",
                    "start": 100 + 40 * len(els),
                    "end": 120 + 40 * len(els),
                    "domain": "chrT:D1",
                    "predicted_coding": {
                        "gene": f"G{i}",
                        "action": "activates" if v < 0 else "represses",
                        "log2_fold_change": v,
                        "tissue": "K562",
                        "strength": "strong" if abs(v) > 0.5 else "weak",
                        "confidence": round(abs(v), 3),
                    },
                }
            )
    (tmp_path / "enhancer_targets_chrT.json").write_text(json.dumps({"elements": els}))
    return els


@pytest.fixture
def program(tmp_path):
    els = _fixture(tmp_path, EFFECTS)
    text = compile_chromosome("chrT", tmp_path, layer=Layer(chrom="chrT"))
    return els, text, parse(text, "chrT")


def test_effect_magnitude_never_raises_a_stated_certainty(program):
    """Across effects from 0.01 to 5 log2, every element and rule states the same (unstated) certainty."""
    els, _, m = program
    by_effect = sorted(
        (abs(e["predicted_coding"]["log2_fold_change"]), m.entities[e["id"]].confidence) for e in els
    )
    assert len({c for _, c in by_effect}) == 1, by_effect
    rules = [r for r in m.rules]
    assert len(rules) == len(els) and len({r.confidence for r in rules}) == 1
    assert rules[0].confidence == 0.0  # the parser's default: nothing stated


def test_the_record_carries_the_effect_in_its_unit_and_no_probability(program):
    els, text, m = program
    for e in els:
        pc = e["predicted_coding"]
        c = element_certainty(pc)
        assert c.effect_estimate == pc["log2_fold_change"] and "log2" in c.effect_unit
        assert c.probability is None and c.probability_unavailable
        assert c.model_score_name  # the score the target run ranked by is named, not called a probability
        note = m.entities[e["id"]].evidence.note
        assert f"{pc['log2_fold_change']:+.3g} log2" in note and "probability unavailable" in note
    # the full record once in the header, as comments
    assert "# certainty: evidence predicted" in text and "#   probability: unavailable (" in text


def test_the_record_does_not_move_with_magnitude_except_its_effect():
    a = element_certainty({"gene": "G", "log2_fold_change": 0.01, "confidence": 0.01}).to_dict()
    b = element_certainty({"gene": "G", "log2_fold_change": 4.0, "confidence": 4.0}).to_dict()
    for d in (a, b):
        d.pop("effect_estimate"), d.pop("model_score")
    assert a == b


def test_every_stated_probability_names_outcome_population_and_method(program):
    _, text, _ = program
    mentions = 0
    for line in text.splitlines():
        for m in re.finditer(r"probability\W+(\S+)", line):
            mentions += 1
            said = m.group(1)
            if re.match(r"[0-9.]", said):  # a number stated as a probability must carry its record
                assert " of " in line and "calibrated on" in line and " by " in line, line
            else:
                assert any(w in line for w in ("unavailable", "no probability", "not a probability")), line
    assert mentions > 0


def test_stub_genes_and_domains_state_no_constant_certainty(program):
    _, text, m = program
    genes = [e for e in m.entities.values() if e.kind == "gene"]
    domains = [e for e in m.entities.values() if e.kind == "domain"]
    assert genes and domains
    assert {g.confidence for g in genes} == {0.0} and {d.confidence for d in domains} == {0.0}
    assert "confidence: 0.9" not in text and "confidence: 0.40" not in text


def test_region_keeps_the_budget_score_labelled_as_hand_set(program):
    _, _, m = program
    (region,) = [e for e in m.entities.values() if e.kind == "region"]
    assert region.confidence == 0.7
    assert "hand-set" in region.evidence.note and "not a probability" in region.evidence.note
