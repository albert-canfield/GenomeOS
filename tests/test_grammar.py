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
