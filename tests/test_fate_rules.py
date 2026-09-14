"""Terminal fates from measured factors: what a terminal cell reads, the rules, their score, and their
precedence over the observed lineage lookup in the Body."""

import json
from pathlib import Path

import pytest

from genomeos.lang import parse
from genomeos.organism import fate_rules as fr
from genomeos.organism.reference import ReferenceLineage
from genomeos.runtime.body import Body

RESULT = Path("data/results/celegans_fate_rules.json")


def _toy_ref() -> ReferenceLineage:
    rows = [
        ("P0", None, 0, 10, ""),
        ("E", "P0", 10, 20, ""),
        ("MS", "P0", 10, 20, ""),
        ("Ea", "E", 20, None, "intestine"),
        ("Ep", "E", 20, None, "intestine"),
        ("MSa", "MS", 20, None, "muscle"),
        ("MSp", "MS", 20, None, "neuron"),
    ]
    return ReferenceLineage.from_dict(
        {
            "cells": [
                {"id": i, "name": i, "parent": p, "born": b, "divides": d, "dies": None, "tissue": t}
                | {"founder": i[:2] if i != "P0" else "P", "generation": len(i) - 2 if i != "P0" else 0}
                for i, p, b, d, t in rows
            ]
        }
    )


def test_reads_rules_and_score():
    ref = _toy_ref()
    term = fr.embryonic_terminal(ref)
    labels = {c.id: c.tissue for c in term}
    atlas = {"E": ["END-1"], "Ea": ["ELT-2"], "MS": ["HLH-1"], "MSp": ["CND-1", "HLH-1"]}
    inst = fr.instantaneous(ref, atlas, term)
    assert inst["Ep"] == {"END-1"} and inst["MSa"] == {"HLH-1"}  # not tracked: the nearest tracked ancestor
    pred = fr.apply_rules(fr.textbook_rules(), inst)
    assert pred == {"Ea": "intestine", "MSa": "muscle", "MSp": "muscle"}
    s = fr.score(pred, labels)
    assert s["cells"] == 4 and s["decided_by_factors"] == 3 and s["wrong"] == 1
    assert s["after"] == 0.75 and s["factor_only"] == 0.5 and s["tissues"]["neuron"]["lost"] == 1
    assert s["tissues"]["intestine"]["after"] == 1.0 and s["confusion"] == {"neuron->muscle": 1}
    # exposure along the lineage path: the mother's exposure counts for her daughters
    levels = {"ELT-2": {"E": [0.1, 0.1, 10.0, -1, 8], "Ea": [1.0, 0.5, 30.0, 3, 8]}}
    integ = fr.integrated(fr.exposures(ref, levels, term), fraction=0.3)
    assert integ["Ea"] == {"ELT-2"} and integ["Ep"] == set()  # Ep: 10 of the largest 40
    assert fr.integrated(fr.exposures(ref, levels, term), fraction=0.2)["Ep"] == {"ELT-2"}
    # presence at a fraction of the factor's max, from the level table
    assert fr.presence(levels, 0.5) == {"Ea": ["ELT-2"]}
    # a learned list finds the separating factor and the lineage-majority baseline uses none
    rules = fr.learn(
        {"a": {"X"}, "b": {"X"}, "c": {"X"}, "d": {"Y"}, "e": {"Y"}, "f": {"Y"}},
        {"a": "gut", "b": "gut", "c": "gut", "d": "skin", "e": "skin", "f": "skin"},
        min_support=2,
    )
    assert {(r.tissue, r.factors) for r in rules} == {("gut", ("X",)), ("skin", ("Y",))}
    assert fr.lineage_majority({"MSa": "muscle", "MSp": "muscle", "MSq": "neuron"}) == {
        "MSa": "muscle",
        "MSp": "muscle",
        "MSq": "muscle",
    }


def test_factor_rules_take_precedence_over_the_lookup_in_the_body():
    rules = [
        fr.Rule("intestine", ("ELT-2",), "Fukushige et al. 1998"),
        fr.Rule("muscle", ("HLH-1",), "Krause"),
    ]
    src = (
        "module toy.fates\nimport bio.std.development\n"
        "organism O { root: Z; cell_type: Zygote }\n"
        "stage Morphogenesis { from: 0 min }\n"
        "cell_type Neuron { parent: PostMitotic; }\ncell_type Intestine { parent: PostMitotic; }\n"
        "cell_type Muscle { parent: PostMitotic; }\n"
        "decision div { action: divide; when: cell = Z; daughters: A, B; after: 1 min }\n"
        "decision rd_A { action: express; when: cell = A; sets: ELT-2_integrated, HLH-1_integrated }\n"
        "decision fate_A { action: differentiate; when: cell = A; to: Neuron }\n"
        "decision fate_B { action: differentiate; when: cell = B; to: Neuron }\n"
    ) + "\n".join(
        line
        for line in fr.to_bio_fates(rules, integrated_read=True).splitlines()
        if line.startswith("decision")
    )
    src = src.replace(f"cell_type = {fr.TERMINAL_TYPES}", "cell_type = Neuron|Intestine|Muscle")
    b = Body(parse(src)).run(until=5)
    # both rules apply to A; the first in precedence (intestine) wins; B keeps the observed fate
    assert b.cells["A"].cell_type == "Intestine" and b.cells["B"].cell_type == "Neuron"
    text = fr.to_bio_fates(rules, integrated_read=False)
    assert text.index("factor_fate_01") < text.index("factor_fate_00")  # precedence runs bottom to top
    exposure = fr.to_bio_exposure(
        {"A": {"ELT-2"}, "C": {"HLH-1"}}, {"A": {"ELT-2"}, "C": {"UNC-3"}}, {"A": []}, ["ELT-2", "HLH-1"]
    )
    assert "sets: ELT-2_integrated;" in exposure and "sets: UNC-3, HLH-1_integrated;" in exposure


@pytest.mark.skipif(not RESULT.exists(), reason="run scripts/celegans_fates.py")
def test_distilled_fate_rules():
    r = json.loads(RESULT.read_text())
    tb = r["textbook"]
    assert tb["integrated"]["cells"] == 555 and tb["instantaneous"]["cells"] == 555
    # reading the factors over time beats the instantaneous threshold, at the same threshold
    assert tb["integrated"]["wrong"] < tb["instantaneous"]["wrong"]
    assert tb["integrated"]["factor_correct"] > tb["instantaneous"]["factor_correct"]
    assert tb["integrated"]["tissues"]["intestine"]["factor_correct"] == 34
    # the factors know more than the lineage composition: shuffled tissues within a sublineage score lower
    learned = r["learned"]["instantaneous"]
    assert learned["cross_validated"]["factor_only"] > max(
        learned["null_tissues_shuffled_within_sublineage"]["factor_only"]
    )
    body = r["body"]
    assert body["lookup"]["fate_accuracy"] == 1.0 and body["lookup"]["fates_checked"] == 555
    assert 0.8 < body["factors"]["fate_accuracy"] < 1.0
    assert body["factors_to_adult"]["fates_checked"] == 961
    assert (
        Path("data/organisms/celegans/fates.bio").exists()
        and Path("data/organisms/celegans/exposure.bio").exists()
    )
