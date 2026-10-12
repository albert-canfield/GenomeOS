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
