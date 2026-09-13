"""The C. elegans reference lineage: parse, distil, query, generate BioLang."""

import json
from pathlib import Path

import pytest

from genomeos.organism.reference import (
    KNOWLEDGE,
    ReferenceLineage,
    _blast_class,
    distil,
    parse_wormweb,
)


def _node(name, did, level, total, death=0, tissue="", children=()):
    return {
        "id": did.lower(),
        "name": name,
        "did": did,
        "data": {"levelDistance": level, "totalDistance": total, "deathDistance": death, "type": tissue},
        "children": list(children),
    }


TOY_TREE = _node(
    "P0",
    "P0",
    0,
    0,
    children=[
        _node(
            "AB",
            "P0a",
            28,
            28,
            children=[_node("N1", "ABa", 100, 128, tissue="neuron"), _node("ABp", "ABp", 100, 128, death=40)],
        ),
        _node(
            "P1",
            "P1",
            30,
            30,
            children=[
                _node(
                    "EMS",
                    "EMS",
                    40,
                    70,
                    children=[
                        _node(
                            "MS",
                            "MS",
                            900,
                            970,
                            children=[
                                _node("MSa", "MSa", 10, 980, tissue="muscle"),
                                _node("MSp", "MSp", 10, 980, tissue="muscle"),
                            ],
                        ),
                        _node("E", "E", 200, 270, tissue="intestine"),
                    ],
                ),
                _node("P2", "P2", 50, 80, tissue="repro"),
            ],
        ),
    ],
)
TOY = "function getJson() { var json; json = " + json.dumps(TOY_TREE) + "; return json; }"


def test_parse_and_distil_toy():
    ref = ReferenceLineage.from_dict(distil(parse_wormweb(TOY)))
    assert set(ref.cells) == {"P0", "AB", "ABa", "ABp", "P1", "EMS", "MS", "MSa", "MSp", "E", "P2"}
    ab = ref.cells["AB"]
    assert ab.founder == "AB" and ab.generation == 0 and ab.divides == 28 and ab.children == ["ABa", "ABp"]
    assert ref.cells["ABa"].name == "N1" and ref.cells["ABa"].cell_type == "Neuron"
    assert ref.cells["ABp"].dies == 68 and ref.cells["ABp"].terminal
    assert ref.cells["E"].founder == "E" and ref.cells["EMS"].founder == "EMS"
    assert ref.cells["MS"].founder == "MS"  # a blast: born in the embryo, divides after hatching
    assert ref.cells["MSa"].founder == "MS"
    assert ref.count_at(0.5) == 2 and ref.count_at(29) == 3 and ref.count_at(100) == 4
    assert [c.id for c in ref.deaths()] == ["ABp"]
    assert ref.fates() == {"Muscle": 2, "Intestine": 1, "Neuron": 1, "Reproductive": 1}
    stats = ref.cycle_stats()
    assert stats[("AB", 0)]["mean"] == 28.0 and stats[("EMS", 0)]["n"] == 1
    timers = ref.to_bio_timers()
    assert "timer cycle_AB_0 { duration: 28.0 min" in timers and "generation = 0" in timers
    prog = ref.to_bio_program()
    assert "decision div_P0 { action: divide; when: cell = P0; daughters: AB, P1;" in prog
    assert "decision die_ABp { action: die; when: cell = ABp; after: 40 min;" in prog
    assert "decision fate_ABa { action: differentiate; when: cell = ABa; to: Neuron; name: N1;" in prog


def test_blast_classes():
    assert [_blast_class(n) for n in ("V1L", "P1", "TL", "Md", "QR", "Z1", "H2Rp", "W", "G2")] == [
        "V",
        "Pn",
        "T",
        "M",
        "Q",
        "Z",
        "H",
        "W",
        "G",
    ]


@pytest.mark.skipif(not Path(KNOWLEDGE).exists(), reason="run `genomeos data distil --only celegans_lineage`")
def test_distilled_reference_holds():
    ref = ReferenceLineage.load()
    s = ref.summary()
    assert s["cells"] == 2183 and s["deaths"] == 131 and s["alive_adult"] == 961  # 959 somatic + Z2/Z3
    assert s["alive_at"]["0.5"] == 2 and s["fates"]["Neuron"] == 302 and s["fates"]["Intestine"] == 34
    assert ref.cells["E"].parent == "EMS" and ref.cells["P4"].parent == "P3"
    assert ref.cells["AB"].parent == "P0"
    assert {c.founder for c in ref.cells.values()} >= {"AB", "MS", "E", "C", "D", "P", "P4", "V", "Pn", "Z"}
    assert 60 <= len(ref.cycle_stats()) <= 120
