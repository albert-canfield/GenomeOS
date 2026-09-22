"""Why the worm still reads `exposure(lineage)`, and why 414 cells claimed is not 414 cells decided.

`scripts/celegans_fate_read_choice.py` measured two things on 2026-09-22: that the 332-to-414 move in
the cells a factor rule fires on is the crossing re-decisions of §7.5 and nothing in area E's program,
and that `mean(lineage) >= 0.25` wins every published number while losing on the one measure that is
not credit for agreeing with a lock. Both are pinned here, because both are the kind of fact a later
change would move silently: the first is a property of the runtime, the second is a standing invitation
to switch the read for +11.
"""

import json
from pathlib import Path

import pytest

from genomeos.lang import parse_file
from genomeos.organism import fate_rules as fr

RESULT = Path("data/results/celegans_fate_read_choice.json")
ORG = Path("data/organisms/celegans")


def test_the_shipped_program_still_reads_exposure_over_the_lineage():
    """The read the worm ships, checked against the program rather than against a doc. If a lane changes
    it, this test is where the measurement that refused the change has to be answered."""
    text = (ORG / "fates.bio").read_text()
    assert f".exposure(lineage) >= {fr.THRESHOLD_MIN:g}" in text
    assert ".mean(lineage)" not in text
    rules = [d for d in parse_file(ORG / "embryo_factors.bio").decisions if d.id.startswith("factor_fate_")]
    assert rules, "no factor rules in the shipped program"


def test_the_commitment_that_makes_a_crossing_unable_to_say_anything_else():
    """The mechanism the 332-to-414 explanation rests on: the worm declares a commitment that locks
    `cell_type` over every terminal type, so a rule firing at a later crossing can only name the type
    the cell already has. Without this block the added cells would be revisions, not agreements."""
    module = parse_file(ORG / "embryo_factors.bio")
    lock = next((c for c in module.commitments if c.locks == "cell_type"), None)
    assert lock is not None, "the terminal-fate commitment is gone"
    established = set(lock.establish["cell_type"].split("|"))
    assert set(fr.TERMINAL_TYPES.split("|")) <= established


@pytest.mark.skipif(not RESULT.exists(), reason="run scripts/celegans_fate_read_choice.py")
def test_the_move_from_332_to_414_is_the_crossing_re_decisions_and_costs_nothing():
    """82 cells joined, none left, every one of them claimed-and-right, and the published score is 522
    either way -- so the score cannot see the move because there is nothing in it for the score to see."""
    m = json.loads(RESULT.read_text())["the_332_to_414_move"]
    assert m["cells_claimed_with_the_crossing_re_decisions"] == 414
    assert m["cells_claimed_without_them"] == 332
    assert m["cells_only_the_crossings_arm_claims"] == 82
    assert m["cells_only_the_recheck_none_arm_claims"] == 0
    assert m["of_the_added_cells_claimed_and_right"] == 82
    assert m["of_the_added_cells_claimed_and_wrong"] == 0
    assert m["claimed_and_wrong"] == [33, 33]
    assert m["score_with_the_lookup_as_fallback"] == [522, 522]
    assert m["crossings_taken"][1] == 0 and m["crossings_taken"][0] > 0


@pytest.mark.skipif(not RESULT.exists(), reason="run scripts/celegans_fate_read_choice.py")
def test_neither_read_decides_a_fate_the_lineage_had_not_written():
    """The count the bar was on, for every arm of both reads. It is zero and the guard makes it zero, so
    no choice of read can move it: this is what stops a higher score from being a better result."""
    r = json.loads(RESULT.read_text())
    for name, arm in r["arms"].items():
        assert arm["decided_a_fate_the_lookup_had_not"] == 0, name
        assert arm["claimed"] == arm["lookup_wrote_it_first"], name
        assert arm["abstained_wrong"] == 0, name


@pytest.mark.skipif(not RESULT.exists(), reason="run scripts/celegans_fate_read_choice.py")
def test_the_mean_read_wins_the_headline_and_loses_the_measure_the_bar_was_on():
    """The finding, in the shape that makes it a finding: `mean` is ahead on everything the earlier
    measurement published and behind on free credit, because nearly all of its claims are made at a
    crossing, after the fate is written and locked."""
    r = json.loads(RESULT.read_text())
    bar = r["the_registered_bar"]
    assert bar["honest_in_sample_and_held_out"]["mean"] > bar["honest_in_sample_and_held_out"]["exposure"]
    assert (
        bar["published_score_with_the_lookup_as_fallback"]["mean"]
        > bar["published_score_with_the_lookup_as_fallback"]["exposure"]
    )
    assert bar["free_credit_in_sample"]["mean"] < bar["free_credit_in_sample"]["exposure"]
    assert (
        not bar["clause_a_independence_rises"] and not bar["clause_b_free_credit_rises_by_more_than_5_fates"]
    )
    assert bar["registered_bar_met"] is False
    assert "stays as it is" in bar["decision"]
    # the generous reading of "free" is the recheck=none arm, and it orders them the same way
    arms = r["arms"]
    assert (
        arms["mean(lineage) >= 0.25 | recheck=none"]["honest"]
        < arms["exposure(lineage) >= 15 (shipped) | recheck=none"]["honest"]
    )
    # and the reason: where each read's claims are made
    shipped = arms["exposure(lineage) >= 15 (shipped) | recheck=crossings"]["when_the_rule_fired"]
    mean = arms["mean(lineage) >= 0.25 | recheck=crossings"]["when_the_rule_fired"]
    for f in (shipped, mean):
        assert f["fired_again_at_a_crossing"] == f["already_committed_to_the_type_the_rule_named"]
    assert mean["fired_only_at_the_cell_s_birth"] < shipped["fired_only_at_the_cell_s_birth"]


@pytest.mark.skipif(not RESULT.exists(), reason="run scripts/celegans_fate_read_choice.py")
def test_the_registration_was_not_rewritten_to_match_the_run():
    """All six predictions held, so there is nothing to diagnose -- but the discipline is the same one
    542d613 had to use when three of its seven broke: the docstring keeps what it was committed with."""
    text = Path("scripts/celegans_fate_read_choice.py").read_text()
    assert "Point prediction 332; anything within 332 +/- 3 is the same finding." in text
    assert "A rise in H alone is NOT an improvement" in text
    pre = json.loads(RESULT.read_text())["preregistered"]
    assert len(pre) == 6
    assert all(claim["held"] for claim in pre.values())
