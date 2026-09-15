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
READS_RESULT = Path("data/results/celegans_fate_reads.json")


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
    # A carries both factors for 20 minutes and then divides; its daughter Aa reads that carriage as
    # `F.exposure(lineage)`, which is the whole point of the rewrite: nothing is precomputed and the
    # deciding cell reads only what the path to it actually carried. B's branch carries nothing.
    src = (
        "module toy.fates\nimport bio.std.development\n"
        "organism O { root: Z; cell_type: Zygote }\n"
        "stage Morphogenesis { from: 0 min }\n"
        "cell_type Neuron { parent: PostMitotic; }\ncell_type Intestine { parent: PostMitotic; }\n"
        "cell_type Muscle { parent: PostMitotic; }\n"
        "decision div_Z { action: divide; when: cell = Z; daughters: A, B; after: 1 min }\n"
        "decision rd_A { action: express; when: cell = A; sets: ELT-2, HLH-1 }\n"
        "decision div_A { action: divide; when: cell = A; daughters: Aa, Ab; after: 20 min }\n"
        "decision div_B { action: divide; when: cell = B; daughters: Ba, Bb; after: 20 min }\n"
        "decision fate_Aa { action: differentiate; when: cell = Aa; to: Neuron }\n"
        "decision fate_Ba { action: differentiate; when: cell = Ba; to: Neuron }\n"
        "regime r { fates: first }\n"
    ) + "\n".join(line for line in fr.to_bio_fates(rules).splitlines() if line.startswith("decision"))
    src = src.replace(f"cell_type = {fr.TERMINAL_TYPES}", "cell_type = Neuron|Intestine|Muscle")
    b = Body(parse(src)).run(until=40)
    # both rules apply to Aa; the first in precedence (intestine) wins; Ba keeps the observed fate
    assert b.cells["Aa"].cell_type == "Intestine" and b.cells["Ba"].cell_type == "Neuron"
    text = fr.to_bio_fates(rules)
    # the read names its window and its threshold, in minutes, with nothing precomputed
    assert "ELT-2.exposure(lineage) >= 15" in text and fr.INTEGRATED_SUFFIX not in text
    # precedence is explicit: the first rule has the highest priority, and line order carries none
    assert "factor_fate_00" in text and "priority: 2;" in text.split("factor_fate_00")[1].split("\n")[0]
    assert "priority: 1;" in text.split("factor_fate_01")[1].split("\n")[0]
    assert b.summary()["ambiguous_fates"] == 0
    # a rule written against the instantaneous read still names the factor the reader gives the cell
    glia = fr.to_bio_fates(fr.glia_rules())
    assert "PROS-1 = present, NHR-25 = present, SOX-2 = present" in glia
    # the superseded precomputed emitter is kept only so the old read can still be scored
    exposure = fr.to_bio_exposure(
        {"A": {"ELT-2"}, "C": {"HLH-1"}}, {"A": {"ELT-2"}, "C": {"UNC-3"}}, {"A": []}, ["ELT-2", "HLH-1"]
    )
    assert "sets: ELT-2_integrated;" in exposure and "sets: UNC-3, HLH-1_integrated;" in exposure


def test_runtime_features_threshold_the_read_the_body_gives_the_cell():
    reads = {
        "Ea": {"ELT-2.exposure(lineage)": 104.0, "ELT-2.mean(lineage)": 0.29, "ELT-2.exposure(cell)": 0.0},
        "MSa": {"ELT-2.exposure(lineage)": 0.0, "ELT-2.mean(lineage)": 0.0, "ELT-2.exposure(cell)": 0.0},
    }
    assert fr.runtime_features(reads) == {"Ea": {"ELT-2"}, "MSa": set()}
    assert fr.runtime_features(reads, "mean", "lineage", 0.5) == {"Ea": set(), "MSa": set()}
    # a terminal cell decides at birth, so its own window is empty whatever the threshold
    assert fr.runtime_features(reads, "exposure", "cell", 1.0) == {"Ea": set(), "MSa": set()}


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
    from genomeos.lang import parse_file

    program = parse_file("data/organisms/celegans/embryo_factors.bio")
    assert program.regime is not None and program.regime.fates == "first"
    assert all(d.priority > 0 for d in program.decisions if d.id.startswith("factor_fate_"))
    assert Path("data/organisms/celegans/fates.bio").exists()
    # the precomputed lookup is gone: the rules read what the runtime integrates (§7.2a)
    assert not Path("data/organisms/celegans/exposure.bio").exists()


@pytest.mark.skipif(not READS_RESULT.exists(), reason="run scripts/celegans_fate_reads.py")
def test_fate_rules_read_what_the_runtime_integrates():
    from genomeos.lang import parse_file

    text = Path("data/organisms/celegans/fates.bio").read_text()
    assert ".exposure(lineage) >= 15" in text and fr.INTEGRATED_SUFFIX not in text
    program = parse_file("data/organisms/celegans/embryo_factors.bio")
    reads = {k for d in program.decisions for k in d.when if ".exposure(" in k or ".mean(" in k}
    assert reads and all(k.endswith("(lineage)") for k in reads)  # every read names its window

    r = json.loads(READS_RESULT.read_text())
    # a terminal cell decides once, at birth, so its own window is empty when it decides: the `cell`
    # reads are zero for every one of the 555, and rules written on them fire for nothing
    assert r["decision_points_per_terminal_cell"] == {"1": 555}
    assert r["terminal_cells_whose_cell_window_is_empty_when_they_decide"] == 555
    for name in ("exposure(cell) >= 1", "mean(cell) >= 0.01"):
        assert r["arms"][name]["in_sample"]["decided_by_factors"] < 10
    # the ordering the precomputed read showed survives a real runtime: on the cited rules alone, with
    # no fitting anywhere, the instantaneous read makes far more errors than the integrated one
    rows = {(x["read"], x["threshold"]): x for x in r["threshold_plateau"]}
    assert rows[("instantaneous (F = present)", None)]["wrong"] > rows[("exposure(lineage)", 15)]["wrong"]
    # and the threshold is a plateau, not a fitted knob
    assert {rows[("exposure(lineage)", t)]["wrong"] for t in (1, 15, 30, 45, 60)} == {
        rows[("exposure(lineage)", 15)]["wrong"]
    }
    arms = r["arms"]
    base = arms["baseline: precomputed _integrated lookup"]
    new = arms[f"exposure(lineage) >= {fr.THRESHOLD_MIN:g}"]
    assert base["in_sample"]["fates_checked"] == new["in_sample"]["fates_checked"] == 555
    # the rewrite is ahead of the lookup it replaces, in sample and held out, and claims fewer cells
    assert new["in_sample"]["fates_correct"] > base["in_sample"]["fates_correct"]
    assert new["held_out"]["fates_correct"] > base["held_out"]["fates_correct"]
    assert new["in_sample"]["decided_by_factors"] < base["in_sample"]["decided_by_factors"]
