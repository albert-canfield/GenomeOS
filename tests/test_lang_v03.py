"""BioLang v0.3: organism, stage, timer, signal, decision, domain blocks; diamond imports."""

import pytest

from genomeos.ir import Module
from genomeos.lang import BioLangError, parse, parse_file

SRC = """
module test.v03
import bio.std.development

organism Toy {
  species: Toyus toyus; genome: none; tempo: 2.0; root: Z0
  cell_type: Zygote
  factors: SKN-1, PIE-1
  environment: temperature = 20
  observe: count, fates
  assert: count at 100 min in 2..4
  assert: deaths at 200 min = 0
  reference: toy
  evidence: curated "toy"; confidence: 0.5
}
stage Cleavage { from: 0 min; to: 2 h }
timer ab_cycle {
  duration: 20 min; sd: 2; lengthening: 1.1; when: lineage = AB
  evidence: experimental "x"; confidence: 0.7
}
signal Notch1 {
  mode: contact; ligand: APX-1; receptor: GLP-1; from: cell = P2; to: cell = ABp; sets: Notch = received
}
decision first {
  action: divide; when: cell = Z0; daughters: AB, P1; lineages: AB = AB, P1 = P; asymmetric: PIE-1 -> P1
}
decision ab { action: divide; when: lineage = AB, generation = >=1; timer: ab_cycle }
decision germ { action: differentiate; when: cell = P1; to: GermCell; name: P1 }
decision gone { action: die; when: cell = ABp; after: 1 h }
domain node1 { locus: chr21:1000-2000; genes: APP; boundaries: E1, E2 }
"""


def test_v03_blocks_compile_and_round_trip():
    m = parse(SRC)
    o = m.organism
    assert o.name == "Toy" and o.root == "Z0" and o.tempo == 2.0 and o.factors == ["SKN-1", "PIE-1"]
    assert o.environment == {"temperature": "20"} and o.observe == ["count", "fates"] and len(o.asserts) == 2
    assert m.stages[0].end == 120.0 and m.stages[0].contains(119) and not m.stages[0].contains(120)
    assert m.stage_at(10) == "Cleavage" and m.stage_at(500) == ""
    t = m.timers[0]
    assert t.sd == 2 and t.lengthening == 1.1 and t.applies({"lineage": "AB"}) and not t.applies({})
    sg = m.signals()[0]
    assert sg.sender == {"cell": "P2"} and sg.receiver == {"cell": "ABp"} and sg.sets == "Notch"
    assert sg.value == "received"
    first, ab, germ, gone = m.decisions
    assert first.daughters == ["AB", "P1"] and first.lineages == {"AB": "AB", "P1": "P"}
    assert first.asymmetric == {"PIE-1": "P1"} and ab.timer == "ab_cycle" and germ.to == "GermCell"
    assert gone.after == 60.0 and ab.applies({"lineage": "AB", "generation": "3"})
    assert not ab.applies({"lineage": "AB", "generation": "0"})
    dm = m.domains()[0]
    assert dm.genes == ["APP"] and dm.boundaries == ["E1", "E2"] and dm.locus.length == 1000
    m2 = Module.from_dict(m.to_dict())
    assert (
        m2.organism.factors == o.factors and len(m2.decisions) == 4 and m2.timers[0].when == {"lineage": "AB"}
    )
    assert m2.signals()[0].sets == "Notch" and m2.domains()[0].genes == ["APP"] and m2.stages[0].end == 120.0


@pytest.mark.parametrize(
    "bad, message",
    [
        ("decision x { action: fly; when: a = b }", "decision action"),
        ("decision x { action: differentiate; to: Nope }", "undeclared cell_type"),
        ("decision x { action: divide; timer: nope }", "undeclared timer"),
        ("decision x { action: divide; daughters: A }", "two daughters"),
        ("timer t { duration: 3 parsecs }", "unknown time unit"),
        ("signal s { mode: telepathy }", "signal mode"),
        ("organism A { } organism B { }", "one organism"),
    ],
)
def test_v03_errors(bad, message):
    with pytest.raises((BioLangError, ValueError), match=message):
        parse("module t\n" + bad)


def test_diamond_imports_merge_once():
    m = parse_file("data/organisms/celegans/embryo.bio")
    assert "Zygote" in m.entities and "Neuron" in m.entities
    assert m.organism.name == "Celegans" and len(m.decisions) > 2000 and len(m.timers) > 50
    # mechanism first, lookup after: the founders' EMS decision precedes the generated one
    ids = [d.id for d in m.decisions if d.id == "div_EMS"]
    first = m.decisions[[d.id for d in m.decisions].index("div_EMS")]
    assert (
        len(ids) == 2 and first.asymmetric == {"POP-1": "MS"} and first.evidence.source.startswith("Sulston")
    )


def test_v02_files_still_compile():
    assert len(parse_file("data/demo/cell_context.bio").rules) > 0
    assert parse_file("data/demo/repressilator.bio").name == "synthetic.repressilator"
