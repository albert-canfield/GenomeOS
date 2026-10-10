# SPDX-License-Identifier: Apache-2.0
"""The written grammar is generated from the parser and the IR, and cannot drift from them."""

from pathlib import Path

from genomeos.lang import grammar


def test_grammar_table_matches_the_parser_and_the_ir():
    assert grammar.check() == []
    text = grammar.render()
    assert "### `decision`" in text and "### `design`" in text and "- **Module**:" in text
    assert "`action` | `divide | differentiate | migrate | quiesce | die | express`" in text


def test_grammar_document_is_current():
    written = Path("docs/BIOLANG-GRAMMAR.md").read_text()
    assert written == grammar.render(), "run `python -m genomeos.lang.grammar > docs/BIOLANG-GRAMMAR.md`"


def test_ir_type_list_is_derived_and_complete():
    """2026-10-10 (lane-engdebts): the hand-kept list stood at 26 of 32 and check() could not see it."""
    import dataclasses

    from genomeos.ir import model

    defined = [
        n
        for n, v in vars(model).items()
        if isinstance(v, type) and dataclasses.is_dataclass(v) and v.__module__ == model.__name__
    ]
    assert defined == grammar.IR_TYPES
    assert len(grammar.IR_TYPES) == 32  # exact, 2026-10-10 (lane-engdebts)
    for name in ("Pool", "Allocation", "Cost", "Order", "Competence", "Commitment"):
        assert name in grammar.IR_TYPES
        assert f"- **{name}**:" in grammar.render()
    assert "Locus" not in grammar.IR_TYPES  # imported from genomeos.coords, not a BioIR type


def test_check_names_an_ir_dataclass_the_package_does_not_export(monkeypatch):
    import dataclasses

    from genomeos.ir import model

    @dataclasses.dataclass
    class Planted:
        x: int = 0

    Planted.__module__ = model.__name__
    monkeypatch.setattr(model, "Planted", Planted, raising=False)
    monkeypatch.setattr(grammar, "IR_TYPES", grammar._ir_types())
    assert "Planted" in grammar.IR_TYPES
    assert grammar.check() == ["IR dataclass 'Planted' is not exported in genomeos.ir.model.__all__"]
