"""The response map's section of the Evidence tab reads its endpoint and keeps the evidence statuses apart."""

import re
from pathlib import Path

PAGE = Path(__file__).resolve().parent.parent / "genomeos" / "web" / "static" / "index.html"


def section() -> str:
    text = PAGE.read_text()
    start = text.index("// ---- response map:")
    return text[start : text.index("setTimeout(loadResponseMap);", start)]


def test_the_page_reads_the_response_map_endpoint_and_computes_nothing():
    s = section()
    assert "/api/evidence/response-map?example=globin_k562" in s
    assert "fetch(" not in s  # through the page's own api() helper only
    assert 'id="rmap"' in PAGE.read_text()


def test_the_four_statuses_are_shown_apart_with_their_meanings():
    s = section()
    for status, meaning in (
        ("observed", "what an experiment measured; never a mechanism"),
        ("predicted", "model or compiled predictions"),
        ("inferred", "interpretations inherited from existing sources"),
        ("unknown", "questions or missing relationships, not facts"),
    ):
        assert re.search(rf"{status}: '{re.escape(meaning)}'", s), status


def test_an_unresolved_relation_names_its_candidates_and_a_claim_is_not_given_a_relation():
    s = section()
    assert "relation unresolved between" in s
    assert "compiled claim (no relation type among the eight)" in s


def test_the_section_claims_no_validation_or_cause():
    s = section().lower()
    for word in ("causal", "identified", "validated", "rescued"):
        assert word not in s, word
