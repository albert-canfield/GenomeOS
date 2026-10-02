"""The sourced-SOX17-repressor registration, checked as committed code.

The list is closed, every entry carries a citation that records what was FETCHED and what the
source does and does not claim, every exclusion is named with its reason, the bound and the number
convention are the prior lane's objects rather than copies of its words, and no entry touches the
sign correction at 5cbce26.
"""

from __future__ import annotations

import pytest

from genomeos.lang import parse_file
from scripts import mesoderm_diagnosis as prior
from scripts.sox17_repressor_search import BOUND_LAST_FRACTION, BOUND_SHARE, REGISTRATION


@pytest.fixture(scope="module")
def reg() -> dict:
    return REGISTRATION


def test_the_list_is_closed_and_its_count_is_registered(reg: dict) -> None:
    ids = [e["id"] for e in reg["closed_list"]]
    assert ids == ["C1", "C2", "C3"]
    assert len(ids) == len(set(ids)) == reg["closed_list_n"] == 3
    assert "set-one-in only" in reg["no_pairwise"]


def test_every_candidate_cites_a_source_that_was_actually_fetched(reg: dict) -> None:
    for e in reg["closed_list"]:
        assert any(c.isdigit() for c in e["citation"]), e["id"]
        fetched = e["fetched"]
        assert fetched["level"] in {"abstract", "full text"}, e["id"]
        assert fetched["url"].startswith("https://www.ebi.ac.uk/europepmc/"), e["id"]
        assert e["source_claims"].strip(), e["id"]
        # the entry must say what the source does NOT claim, not only what it does
        assert e["source_does_not_claim"].strip(), e["id"]


def test_every_candidate_represses_sox17_and_none_touches_the_corrected_sign(reg: dict) -> None:
    for e in reg["closed_list"]:
        assert e.get("remove") is None, e["id"]  # nothing is removed: R1 is a locator, not a remedy
        targets = {(r["src"], r["dst"], r["action"]) for r in e["rules"]}
        assert any(dst == "SOX17" and action == "inhibits" for _, dst, action in targets), e["id"]
        for src, dst, _ in targets:
            assert not (src == "Tbxt" and dst == "SOX17"), e["id"]


def test_every_added_rule_carries_its_own_citation(reg: dict) -> None:
    for e in reg["closed_list"]:
        for rule in e["rules"]:
            assert any(c.isdigit() for c in rule["cite"]), (e["id"], rule)


def test_exclusions_are_named_with_reasons(reg: dict) -> None:
    assert [x["id"] for x in reg["excluded"]] == ["Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Y7", "Y8"]
    for x in reg["excluded"]:
        assert "EXCLUDED" in x["reason"], x["id"]
        assert x["citation"].strip(), x["id"]


def test_the_bound_and_the_convention_are_the_prior_lanes_objects(reg: dict) -> None:
    """Not the same words: the same objects, so 'the same weak bound' cannot drift."""
    assert reg["bound"] is prior.REGISTRATION["bound"]
    assert reg["run"] is prior.REGISTRATION["run"]
    assert reg["number_convention"] is prior.REGISTRATION["number_convention"]
    assert BOUND_SHARE is prior.BOUND_SHARE == 0.01
    assert BOUND_LAST_FRACTION is prior.BOUND_LAST_FRACTION == 0.1
    assert reg["bound"]["share"] < 0.694  # the CS7 mesoderm range is not the bound


def test_the_three_readings_are_carried_without_strengthening(reg: dict) -> None:
    m1, r1, structural = reg["builds_on"]["readings_carried_unchanged"]
    assert "CLEARS THE WEAK BOUND" in m1 and "0.0167" in m1
    assert "does not restore mesoderm" in m1
    assert "LOCATES WHERE THE BAND RESTS" in r1 and "not a fix" in r1
    assert "the rule stays" in r1
    assert "names no biological factor" in structural


def test_the_reading_is_one_of_exactly_two(reg: dict) -> None:
    assert "one of exactly two" in reg["question"]
    assert "restores a mesoderm band" in reg["question"]
    assert "no sourced repressor on the list does" in reg["question"]


def test_the_refusal_to_tune_is_inherited_verbatim(reg: dict) -> None:
    inherited = prior.REGISTRATION["what_this_is_not"]
    assert reg["what_this_is_not"][: len(inherited)] == inherited
    assert any("8bb9123" in s and "19423bd" in s for s in reg["what_this_is_not"])
    assert reg["number_convention"]["not_done"].startswith("no strength, threshold, hill")


def test_the_false_header_is_reported_and_not_edited(reg: dict) -> None:
    hc = reg["header_consequence"]
    assert "Mutual repression makes the choice sharp" in hc["line"]
    assert hc["status_before_this_lane"].startswith("FALSE")
    assert "no run in this file edits data/demo/gastrulation.bio" in hc["this_lane_does_not_change_it"]


def test_the_prior_lanes_seven_each_get_a_verdict(reg: dict) -> None:
    ver = reg["prior_lane_verification"]
    assert [e["id"] for e in ver["entries"]] == ["W1", "W2", "F1", "F2", "M1", "E1", "A1"]
    for e in ver["entries"]:
        assert e["verdict"].startswith(("CONFIRMED", "MIS-CITED", "NOT FOUND")), e["id"]
        assert "PMID" in e["record"], e["id"]
        assert e["supports_the_rule"].strip(), e["id"]
    # the one verdict that qualifies the prior row, stated and not softened
    m1 = next(e for e in ver["entries"] if e["id"] == "M1")
    assert "ITS TWO RULES ARE NOT CONFIRMED BY IT" in m1["verdict"]
    assert "CLEARS THE WEAK BOUND" in m1["supports_the_rule"]


def test_the_module_still_states_the_corrected_sign() -> None:
    module = parse_file("data/demo/gastrulation.bio")
    signs = {(r.source, r.target): r.action.name for r in module.rules}
    assert signs[("Tbxt", "SOX17")] == "ACTIVATE"
    assert signs[("Sox17", "TBXT")] == "INHIBIT"
    # the structural fact this lane exists to answer: SOX17 is repressed by nothing
    assert not [r for r in module.rules if r.target == "SOX17" and r.action.name == "INHIBIT"]
