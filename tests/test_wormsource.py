# SPDX-License-Identifier: AGPL-3.0-or-later
"""The wormsource registration, fixed before any value is read and before anything is fetched."""

from __future__ import annotations

from pathlib import Path

from genomeos.attribution import simpler
from genomeos.attribution import wormsource as ws

ROOT = Path(__file__).resolve().parents[1]


def test_the_stop_rule_is_CARRIED_as_the_same_object() -> None:
    assert ws.MIN_DISCORDANT is simpler.MIN_DISCORDANT
    assert ws.stated_min_discordant() == 200


def test_both_target_figures_are_still_in_the_record_they_were_read_from() -> None:
    """Pinned to docs/SIMPLER-GENOMES.md rather than retyped, as genome4 did."""
    text = (ROOT / ws.SOURCE_OF_RECORD).read_text()
    assert ws.targets_in_the_record(text)
    assert ws.PRIMARY_TARGET == "+3.73 [+2.79, +4.74]"
    assert ws.PRIMARY_ELEMENTS == 2335
    assert ws.OUTRON_SURVIVING_TARGET == "+4.47 [+3.16, +5.78]"
    assert ws.OUTRON_SURVIVING_ELEMENTS == 1745


def test_the_primary_is_named_primary_so_the_larger_number_cannot_be_promoted() -> None:
    w = ws.WHICH_IS_PRIMARY
    assert w.startswith("the PRIMARY is +3.73")
    assert "is secondary" in w
    assert "did not call it a strengthening" in w


# ------------------------------------- the compound confound: source, stage and assay together


def test_the_three_differences_are_FIELDS_not_prose() -> None:
    """A later reader must see all three without reading a sentence."""
    d = ws.THREE_DIFFERENCES
    assert set(d) == {"stage", "assay", "source"}
    assert d["stage"]["durham_2021"] == "L2 larvae"
    assert d["stage"]["serizay_2020"] == "young adults"
    assert "single-cell" in d["assay"]["durham_2021"]
    assert "bulk" in d["assay"]["serizay_2020"]
    for axis in d.values():
        assert set(axis) == {"serizay_2020", "durham_2021"}


def test_every_outcome_speaks_of_source_stage_and_assay_and_none_of_source_alone() -> None:
    for name, outcome in ws.OUTCOMES.items():
        assert not ws.result_attributes_only_to_source(outcome["concludes"]), name
        assert not ws.result_attributes_only_to_source(outcome["falsifies"]), name


def test_the_NOT_REPLICATED_wording_is_carried_word_for_word() -> None:
    assert ws.NOT_REPLICATED_WORDING == (
        "the lead does not hold on an L2 single-cell source; whether that is source, stage or"
        " assay is not separable here"
    )
    assert ws.OUTCOMES["NOT_REPLICATED"]["concludes"] == ws.NOT_REPLICATED_WORDING


def test_a_reading_that_blames_SOURCE_ALONE_is_REFUSED() -> None:
    """The refusal that matters. Here the weakening a summary produces is a SHARPENING."""
    for wrong in (
        "the lead is source-dependent",
        "The worm's lead is source-dependent and must be qualified downstream.",
        "it depends on the source",
    ):
        assert ws.result_attributes_only_to_source(wrong), wrong
    # and the registered sentence, which names all three, must pass
    assert not ws.result_attributes_only_to_source(ws.NOT_REPLICATED_WORDING)
    assert not ws.result_attributes_only_to_source(
        "source-dependent, though stage and assay differ too and are not separable"
    )


def test_REPLICATED_is_stronger_AND_still_does_not_isolate() -> None:
    """So the positive cannot be over-read either."""
    c = ws.OUTCOMES["REPLICATED"]["concludes"]
    assert "MORE ROBUST" in c
    assert "survived three changes at once" in c
    assert "STILL DOES NOT ISOLATE WHICH FACTOR" in c


def test_the_tissue_intersection_is_necessary_and_NOT_sufficient() -> None:
    r = ws.TISSUE_INTERSECTION_RULE
    assert "RE-READ" in r and "SERIZAY ON THAT SAME INTERSECTION" in r
    assert "NECESSARY AND" in r and "NOT SUFFICIENT" in r
    assert "neither stage nor assay" in r
    assert "removes a FOURTH difference" in ws.WHAT_IS_COMPARED


# ---------------------------------------------------------------- the limit, and the rest of the rules


def test_a_result_that_drops_OR_SOFTENS_the_limit_is_refused() -> None:
    assert ws.result_carries_its_limit({ws.LIMIT_KEY: ws.NOT_THE_FOURTH_GENOME_QUESTION})
    assert not ws.result_carries_its_limit({})
    weakened = (
        "this is the same species, so it is a slightly narrower comparison than a fourth"
        " genome would have been"
    )
    assert not ws.result_carries_its_limit({ws.LIMIT_KEY: weakened})


def test_the_limit_says_it_cannot_touch_the_fourth_genome_question() -> None:
    limit = ws.NOT_THE_FOURTH_GENOME_QUESTION
    assert "CANNOT ADDRESS THE FOURTH-GENOME QUESTION AT ALL" in limit
    assert "same operons" in limit and "SL trans-splicing" in limit


def test_two_of_the_four_outcomes_conclude_nothing_and_are_named() -> None:
    assert set(ws.WHICH_OUTCOMES_CONCLUDE_NOTHING) == {"INCONCLUSIVE", "COMPARISON_CANNOT_BE_MADE"}
    for name in ws.WHICH_OUTCOMES_CONCLUDE_NOTHING:
        assert "nothing" in ws.OUTCOMES[name]["concludes"].lower()


def test_cannot_be_made_is_a_complete_answer_and_not_a_licence_to_widen() -> None:
    c = ws.COMPARISON_CANNOT_BE_MADE_IS_COMPLETE
    assert "COMPLETE ANSWER" in c
    assert "not a licence to widen the tissue set" in c


def test_the_endpoint_is_not_adjusted_to_fit_the_new_source() -> None:
    e = ws.ENDPOINT_IS_NOT_ADJUSTED
    assert "A FINDING TO REPORT BEFORE" in e
    assert "NOT A CHANGE TO MAKE" in e


def test_one_comparison_one_endpoint() -> None:
    assert "ONE pre-specified comparison on ONE pre-specified endpoint" in ws.MULTIPLICITY
    assert "nobody may later read a chosen one" in ws.MULTIPLICITY


def test_the_order_rule_and_the_data_terms_are_registered() -> None:
    assert "A record" in ws.ORDER and "never one commit" in ws.ORDER
    assert "3 failed and 7,258 passed" in ws.ORDER
    assert "Nothing is fetched before the registration is committed" in ws.ORDER
    assert "RELATIVE input_entry" in ws.DATA_TERMS
    assert "STOP AND REPORT" in ws.DATA_TERMS
    assert "0 model" in ws.DATA_TERMS


def test_the_stop_list_is_read_and_named_and_the_biosample_clause_addressed() -> None:
    src = (ROOT / "genomeos" / "attribution" / "wormsource.py").read_text()
    assert "2026-09-27 STOP LIST was read and this lane touches none of its six items" in src
    assert "READER_TERMS" in src
    assert "adds none, names none and widens none" in src


def test_the_refusal_checker_states_its_own_blind_spot() -> None:
    """A string match cannot tell a claim from a prohibition of that claim; it says so."""
    doc = ws.result_attributes_only_to_source.__doc__ or ""
    # the halves, not the phrase: the docstring wraps as "a CLAIM from a\n    PROHIBITION",
    # and a flat substring cannot match across the break. Second time today, same trap.
    assert "cannot tell a CLAIM from a" in doc
    assert "PROHIBITION of that claim" in doc
    assert "reworded rather than the check loosened" in doc
    # and the blind spot is real, demonstrated rather than asserted
    prohibition = "this must never be called source-dependent"
    assert ws.result_attributes_only_to_source(prohibition), (
        "the blind spot the docstring describes must actually be present, or the warning is stale"
    )
