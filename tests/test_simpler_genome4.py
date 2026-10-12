# SPDX-License-Identifier: AGPL-3.0-or-later
"""The fourth genome's inventory criteria, fixed before any source is opened."""

from __future__ import annotations

from pathlib import Path

from genomeos.attribution import simpler
from genomeos.attribution import simpler_genome4 as g4

ROOT = Path(__file__).resolve().parents[1]


def test_the_stop_rule_is_CARRIED_and_not_restated() -> None:
    """One definition, so the fourth genome cannot be held to a different 200 from the first three."""
    assert g4.MIN_DISCORDANT is simpler.MIN_DISCORDANT
    assert g4.MIN_DISCORDANT == 200


def test_the_three_species_numbers_match_the_records_they_came_from() -> None:
    """Carried wording must be checkable against its source, or it is a paraphrase."""
    doc = (ROOT / "docs" / "SIMPLER-GENOMES.md").read_text()
    src = (ROOT / "genomeos" / "attribution" / "simpler_genome4.py").read_text()
    for figure in ("+0.26 points [-0.01, +0.53]", "0.74 points"):
        assert figure in doc, f"{figure} is no longer in the record"
        assert figure in src, f"{figure} is not carried into the module"
    # the worm's own reading, which this lane must not restate differently
    assert "+4.47 [+3.16, +5.78]" in src
    assert "1,745" in src


def test_the_question_names_why_three_genomes_cannot_answer_it() -> None:
    assert "the only genome with the mechanism is also the only genome with the lead" in g4.QUESTION


def test_each_candidate_family_says_what_it_CANNOT_settle() -> None:
    """A candidate admitted without that is a candidate that will be over-read."""
    assert len(g4.CANDIDATE_FAMILIES) == 2
    nematode, other = g4.CANDIDATE_FAMILIES
    assert "C. elegans specifically from nematodes generally" in nematode["separates"]
    assert "the mechanism" in nematode["cannot_settle"]
    assert "operons and trans-splicing as a mechanism" in other["separates"]
    assert "a single fourth genome cannot" in other["cannot_settle"]
    for family in g4.CANDIDATE_FAMILIES:
        assert family["cannot_settle"], family["family"]


def test_the_four_criteria_are_all_present_and_the_same_source_clause_is_explicit() -> None:
    q = g4.QUALIFIES
    assert "ALL FOUR" in q
    assert "FROM THE SAME SOURCE" in q
    assert "between-study difference" in q
    assert "MIN_DISCORDANT" in q


def test_what_disqualifies_is_named_so_the_inventory_cannot_relax() -> None:
    d = g4.DISQUALIFIES
    assert "DIFFERENT studies" in d
    assert "no tissue resolution" in d
    assert "without strands" in d
    assert "UNDERPOWERED" in d


def test_no_source_qualifies_is_declared_a_complete_answer_in_advance() -> None:
    n = g4.NO_SOURCE_QUALIFIES_IS_COMPLETE
    assert "COMPLETE ANSWER and the likely one" in n
    assert "NOT CURRENTLY DECIDABLE on public data" in n


def test_the_wording_constraint_is_carried_word_for_word() -> None:
    w = g4.WORDING_CONSTRAINT
    assert w.startswith("say 'on these sources', nothing wider")
    assert "does not license any statement about the domain" in w


def test_the_stop_list_is_named_and_the_biosample_clause_is_addressed() -> None:
    """ROADMAP section 5: a lane brief names any stop item it touches."""
    src = (ROOT / "genomeos" / "attribution" / "simpler_genome4.py").read_text()
    assert "THE STOP LIST OF 2026-09-27 WAS READ, NOT IGNORED" in src
    assert "touches none of its six items" in src
    assert "READER_TERMS" in src
    assert "no AlphaGenome requests on the 882" in src
    assert "makes no model request" in src


def test_it_does_not_claim_a_fourth_genome_would_settle_the_question() -> None:
    assert "one fourth genome separates two" in g4.DOES_NOT_CLAIM
    assert "leaves every other difference" in g4.DOES_NOT_CLAIM


# ------------------------------- the limit, enforced as a limit rather than written as a caveat


def test_a_result_WITHOUT_the_limit_is_REFUSED_by_its_own_check() -> None:
    """The shape globcensus's invisible-shrink test has: the limit cannot be dropped silently."""
    assert g4.result_carries_its_limit({g4.LIMIT_KEY: g4.SECOND_NEMATODE_LIMIT})
    assert not g4.result_carries_its_limit({}), "a result with no limit must be refused"
    assert not g4.result_carries_its_limit({g4.LIMIT_KEY: None})
    # and a WEAKENED restatement must also be refused, which is the case that matters
    weakened = (
        "if the source is another nematode the result is a little less general, which is"
        " worth bearing in mind when reading the comparison"
    )
    assert not g4.result_carries_its_limit({g4.LIMIT_KEY: weakened}), (
        "a softened restatement of the limit must not pass as the limit"
    )


def test_the_limit_says_which_comparison_reaches_the_question() -> None:
    limit = g4.SECOND_NEMATODE_LIMIT
    assert limit.startswith("IF THE QUALIFYING SOURCE IS A SECOND NEMATODE")
    assert "separates neither" in limit
    assert "only that comparison reaches the question as asked" in limit


def test_the_search_protocol_is_declared_before_searching() -> None:
    """A negative's boundary IS the negative, so the archives and species are fixed in advance."""
    p = g4.SEARCH_PROTOCOL
    assert len(p["species"]) == 5
    assert any("Ciona" in s for s in p["species"]), "the non-nematode arm must be in the protocol"
    assert any("ANISEED" in a for a in p["archives"])
    assert "at most 4" in p["queries_per_species"]
    assert "never as 'no such data exists'" in p["what_a_search_can_and_cannot_show"]


def test_a_negative_must_travel_with_its_boundary() -> None:
    n = g4.NEGATIVE_CARRIES_ITS_BOUNDARY
    assert "IN THE SAME SENTENCE" in n
    assert "which archives" in n and "how many queries" in n and "which species" in n
