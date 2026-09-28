# SPDX-License-Identifier: AGPL-3.0-or-later
"""R7 follow-up: the classifier records origin beside the class, and the budget's new writes read
constraint as selection only, never missing constraint as no function."""

from __future__ import annotations

import json
import re

import pytest

from genomeos.attribution import budget
from genomeos.attribution import compile as comp
from genomeos.genome.repeats import Repeat, RepeatIndex
from genomeos.genome.unknown import REPEAT_DERIVED_MIN, classify_block, load_patterns, repeat_origin

SEQ = "".join("ACGTTGCA"[i % 8] for i in range(8000))
CCRE = {"PLS": 1, "dELS": 2}
BANNED = re.compile(r"neutral|fossil|dead|best guess|lineage-specific|exapted|no evidence of function")


def test_a_regulatory_block_that_is_mostly_line_says_both():
    pats = load_patterns()
    cls, ev, conf, f, _ = classify_block(SEQ, None, pats, CCRE, {"LINE": 6000, "SINE": 400})
    assert (cls, ev, conf) == ("regulatory", "curated", 0.7)
    assert f["origin"] == "repeat_derived/LINE" and f["interspersed_coverage"] == 0.8
    # the class, evidence and score are what they were without RepeatMasker
    assert classify_block(SEQ, None, pats, CCRE)[:3] == (cls, ev, conf)


def test_origin_is_recorded_when_rmsk_was_consulted_and_found_nothing():
    pats = load_patterns()
    with_empty = classify_block(SEQ, None, pats, None, {})
    without = classify_block(SEQ, None, pats, None, None)
    assert with_empty[:3] == without[:3]
    assert with_empty[3]["origin"] == "unique" and "origin" not in without[3]


@pytest.mark.parametrize(
    ("cov", "n"), [({"LINE": 50, "SINE": 10}, 100), ({"LTR": 20, "Simple_repeat": 70}, 100), ({}, 100)]
)
def test_the_classifier_origin_is_the_compiler_origin(cov, n):
    o = repeat_origin(cov, n)
    inter = {k: v for k, v in cov.items() if k in comp._INTERSPERSED}
    top = max(inter, key=lambda k: inter[k]) if inter else ""
    assert o["origin"] == comp._origin(sum(inter.values()) / n, top)
    assert REPEAT_DERIVED_MIN == comp.REPEAT_DERIVED_MIN
    assert repeat_origin(None, n) == {"origin": "unknown"}


def test_budget_tables_match_the_compiler():
    assert budget.CLASS_ORIGIN == comp.CLASS_ORIGIN
    for fa in (None, 0.0, 0.029, 0.03, 0.049, 0.05, 0.9):
        assert budget.selection(fa) == comp.selection(fa)


CLASSES = (
    "interspersed_repeat_LINE",
    "regulatory",
    "promoter_like",
    "long_orf",
    "unique_intergenic",
    "mixed_intergenic",
    "unclassified",
    "gene_desert",
    "centromere",
    "satellite_array",
    "gap",
)


@pytest.mark.parametrize("cls", CLASSES)
def test_no_new_label_reads_missing_constraint_as_no_function(cls):
    for fa in (None, 0.0, 0.01, 0.035, 0.06, 0.5):
        for n in (0, 4):
            g = budget.guess(cls, 0.5, None if fa is None else {"fraction_above": fa}, {"n": n})
            assert not BANNED.search(g["label"]), g["label"]
            assert g["tier"] in budget.AXES_TIERS and g["legacy_tier"] in budget.TIERS
            assert budget.LEGACY_TIER.get(g["tier"], g["tier"]) == g["legacy_tier"]
            assert g["evidence_status"] == (None if cls == "gap" else budget.selection(fa))


def test_guess_origin_from_class_or_the_classifier():
    reg = budget.guess("regulatory", 0.7, {"fraction_above": 0.01}, None, {"origin": "repeat_derived/LINE"})
    assert reg["tier"] == "regulatory" and reg["origin"] == "repeat_derived/LINE"
    assert budget.guess("interspersed_repeat_SINE", 0.9, {"fraction_above": 0.0}, None)["origin"] == (
        "repeat_derived/SINE"
    )
    assert budget.guess("gap", 1.0, None, None)["origin"] == "assembly_gap"
    assert budget.guess("unique_intergenic", 0.3, {"fraction_above": 0.0}, None)["origin"] == "unknown"


def _stored(tmp_path):
    blocks = [
        {
            "start": 0,
            "end": 100,
            "length": 100,
            "class": "gap",
            "class_confidence": 1.0,
            "phylop": None,
            "elements": None,
            "guess": {
                "tier": "structural",
                "label": "assembly gap: no sequence to attribute",
                "confidence": 1.0,
            },
        },
        {
            "start": 100,
            "end": 300,
            "length": 200,
            "class": "regulatory",
            "class_confidence": 0.7,
            "phylop": {"bases": 200, "above": 2, "fraction_above": 0.01},
            "elements": {"n": 0},
            "guess": {"tier": "regulatory", "label": "x", "confidence": 0.5},
        },
        {
            "start": 300,
            "end": 1000,
            "length": 700,
            "class": "unique_intergenic",
            "class_confidence": 0.3,
            "phylop": {"bases": 700, "above": 7, "fraction_above": 0.01},
            "elements": {"n": 0},
            "guess": {"tier": "neutral", "label": "y", "confidence": 0.6},
        },
    ]
    (tmp_path / "budget_chrT.json").write_text(
        json.dumps({"chromosome_length": 2000, "unknown_bp": 1000, "threshold": 2.27, "blocks": blocks})
    )
    (tmp_path / "unknown_chrT.json").write_text(json.dumps({"blocks": [{"start": 100, "features": {}}]}))
    return tmp_path


def test_restate_writes_a_new_name_and_leaves_the_stored_budget(tmp_path):
    rd = _stored(tmp_path)
    before = (rd / "budget_chrT.json").read_bytes()
    idx = RepeatIndex(
        [Repeat(100, 280, "LINE", "L1", "L1PA2", 0.1), Repeat(400, 450, "SINE", "Alu", "AluY", 0.0)]
    )
    out = budget.restate("chrT", rd, rindex=idx)
    tiers = [b["guess"]["tier"] for b in out["blocks"]]
    assert tiers == ["structural", "regulatory", "unconstrained_unknown"]
    assert [b["guess"]["legacy_tier"] for b in out["blocks"]] == ["structural", "regulatory", "neutral"]
    assert out["blocks"][1]["guess"]["origin"] == "repeat_derived/LINE"
    assert out["blocks"][2]["guess"]["origin"] == "partly_repeat_derived/SINE"
    assert out["blocks"][2]["guess"]["evidence_status"] == "selection_not_detected"
    assert out["by_tier"]["unconstrained_unknown"]["bp"] == 700 and "neutral" not in out["by_tier"]
    assert (rd / "budget_chrT.json").read_bytes() == before
