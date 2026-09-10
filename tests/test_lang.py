import json

import pytest

from genomeos.ir import UNKNOWN, Action, Module
from genomeos.lang import BioLangError, parse, parse_file


def test_parse_repressilator():
    m = parse_file("data/demo/repressilator.bio")
    assert m.name == "synthetic.repressilator"
    assert {g.id for g in m.genes()} == {"tetR", "lacI", "cI"}
    assert {p.id for p in m.proteins()} == {"TetR", "LacI", "CI"}
    inhibits = [r for r in m.rules if r.action is Action.INHIBIT]
    assert len(inhibits) == 3
    assert [r for r in m.rules if r.action is Action.PRODUCE][0].target == "TetR"
    assert m.parameters["translation_rate"].value == 5
    assert m.unknowns()[0].role is UNKNOWN


def test_roundtrip_json():
    m = parse_file("data/demo/repressilator.bio")
    data = json.loads(json.dumps(m.to_dict()))
    m2 = Module.from_dict(data)
    assert set(m2.entities) == set(m.entities)
    assert len(m2.rules) == len(m.rules)
    assert m2.rules[0].evidence == m.rules[0].evidence
    assert m2.unknowns()[0].role is UNKNOWN


def test_undeclared_reference_is_compile_error():
    src = """
    module bad
    gene a { produces: P }
    """
    with pytest.raises(BioLangError, match="undeclared entity 'P'"):
        parse(src)


def test_bad_evidence_kind():
    src = 'module bad\ngene a { evidence: rumour "x" }\n'
    with pytest.raises(BioLangError, match="unknown evidence kind"):
        parse(src)


def test_rule_context_gate():
    src = """
    module ctx
    gene a { max: 10 }
    protein X { }
    rule X inhibits a { when: cell_type = neuron }
    """
    m = parse(src)
    r = m.rules[0]
    assert r.applies({"cell_type": "neuron"})
    assert not r.applies({"cell_type": "hepatocyte"})
    assert not r.applies({})
