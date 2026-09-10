"""Task 2.1: BioLang v0.2 — nested transcripts, cell types, events, imports."""

import pytest

from genomeos.ir import Action, CellType
from genomeos.lang import BioLangError, parse, parse_file


def test_v01_files_still_compile():
    m = parse_file("data/demo/repressilator.bio")
    assert len(m.rules) == 6 and m.name == "synthetic.repressilator"


def test_nested_transcript_cell_types_events_and_imports():
    m = parse_file("data/demo/cell_context.bio")
    assert "NKX2-5-201" in m.entities
    tx = m.entities["NKX2-5-201"]
    assert len(tx.exons) == 2 and len(tx.cds_segments) == 2 and tx.cds.length > 0
    # imported prelude content is merged
    assert isinstance(m.entities["Cardiomyocyte"], CellType)
    assert any(e.id == "divide" for e in m.events)
    div = next(e for e in m.events if e.id == "divide")
    assert div.rate == 0.5 and div.rate_unit == "/yr"
    assert [(f.target, f.op, f.value, f.unit) for f in div.effects] == [
        ("telomere_bp", "-=", 70.0, "bp"),
        ("divisions", "+=", 1.0, ""),
    ]
    assert "telomere_at_birth_bp" in m.parameters


def test_cell_type_context_silences_unexpressed_genes():
    m = parse_file("data/demo/cell_context.bio")
    assert m.silenced_genes({"cell_type": "Cardiomyocyte2"}) == {"SYN1"}
    assert m.silenced_genes({"cell_type": "Neuron2"}) == {"NKX2-5", "MYH6"}
    heart = m.active_rules({"cell_type": "Cardiomyocyte2"})
    neuron = m.active_rules({"cell_type": "Neuron2"})
    assert any(r.action is Action.ACTIVATE for r in heart)
    assert not any(r.action is Action.ACTIVATE for r in neuron)
    assert m.silenced_genes({}) == set()


def test_json_roundtrip_keeps_events_and_cell_types():
    from genomeos.ir import Module

    m = parse_file("data/demo/cell_context.bio")
    m2 = Module.from_dict(m.to_dict())
    assert len(m2.events) == len(m.events) and m2.events[0].effects[0].op == "-="
    assert isinstance(m2.entities["Neuron2"], CellType) and m2.entities["Neuron2"].expresses == ["SYN1"]


def test_errors():
    with pytest.raises(BioLangError, match="undeclared parent"):
        parse("module x\ncell_type A { parent: Nope }")
    with pytest.raises(BioLangError, match="cannot resolve import"):
        parse("module x\nimport bio.std.nothing_here")
    with pytest.raises(BioLangError, match="nested inside a gene"):
        parse("module x\ntranscript t { exons: c:1-2 }")
    with pytest.raises(BioLangError, match="bad effect"):
        parse("module x\nevent e { effect: telomere ~ 5 }")
