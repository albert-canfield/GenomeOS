"""What the measured factors decide on their own, and why the published fate score cannot show it.

These tests exist because `522 of 555` is a number about a program and not about the factors: the
observed lineage lookup supplies every fate no factor rule claims, and every fate it supplies is right by
construction. The metric that removes it (`scripts/celegans_fate_credit.py`) is therefore the instrument,
and an instrument that quietly credits the lookup would report a large number and be wrong, so the
accounting is tested here rather than trusted.
"""

import json
from dataclasses import dataclass, field
from pathlib import Path

import pytest

from genomeos.lang import parse_file
from genomeos.organism import fate_rules as fr
from genomeos.organism.reference import ReferenceLineage
from genomeos.runtime.body import Body
from scripts.celegans_fate_credit import credit

RESULT = Path("data/results/celegans_fate_credit.json")
ORG = Path("data/organisms/celegans")


@dataclass
class _Cell:
    cell_type: str
    fired: list[str] = field(default_factory=list)


class _Body:
    def __init__(self, cells: dict[str, _Cell]) -> None:
        self.cells = cells


def _ref() -> ReferenceLineage:
    rows = [("a", "intestine"), ("b", "muscle"), ("c", "neuron"), ("d", "neuron")]
    return ReferenceLineage.from_dict(
        {
            "cells": [
                {
                    "id": i,
                    "name": i,
                    "parent": None,
                    "born": 0,
                    "divides": None,
                    "dies": None,
                    "tissue": t,
                    "founder": "AB",
                    "generation": 1,
                }
                for i, t in rows
            ]
        }
    )


def test_the_lookup_gets_no_credit_and_the_rule_that_set_the_type_gets_it():
    """The whole point of the metric: a cell the factors said nothing about is not a fate the factors
    decided, even though the program has it right. `a` and `b` were claimed and are right, `c` was claimed
    and is wrong, `d` was left to the lookup and counts for nothing."""
    ref = _ref()
    targets = {"factor_fate_00": "Intestine", "factor_fate_01": "Muscle", "factor_fate_02": "Muscle"}
    body = _Body(
        {
            "a": _Cell("Intestine", ["fate_a", "factor_fate_00"]),
            "b": _Cell("Muscle", ["fate_b", "factor_fate_01"]),
            "c": _Cell("Muscle", ["fate_c", "factor_fate_02"]),
            "d": _Cell("Neuron", ["fate_d"]),
        }
    )
    c = credit(body, ref, ["a", "b", "c", "d"], targets)
    assert c["with_lookup_fallback"] == 3  # what the published metric reports
    assert c["honest"] == 2 and c["claimed"] == 3 and c["abstained"] == 1  # what the factors decided
    assert c["net_effect_on_the_score"] == -1  # consulting them cost `c`
    # every claimed cell here had the lookup's own decision fire first, as in the worm
    assert c["lookup_wrote_it_first"] == 3 and c["decided_a_fate_the_lookup_had_not"] == 0


def test_credit_refuses_to_be_attributed_to_a_rule_that_did_not_set_the_type():
    """A factor rule can appear in a cell's `fired` list without being the decision the cell's type came
    from. If that ever happens the count is crediting the wrong decision, so it is an assertion and not a
    silently-lower number."""
    ref = _ref()
    body = _Body({"a": _Cell("Neuron", ["fate_a", "factor_fate_00"])})
    with pytest.raises(AssertionError, match="misattributed"):
        credit(body, ref, ["a"], {"factor_fate_00": "Intestine"})


def test_the_factor_rules_cannot_fire_before_the_lookup_has_written_a_fate():
    """Why the honest number has the floor it has. Every generated rule is guarded on a terminal
    `cell_type`, and in this program only the lineage lookup ever puts a terminal type on a cell, so a
    factor rule is a successor to the lookup in the same differentiation chain and never a competitor.
    It follows that the factors can agree with the answer sheet or be wrong, never correct it -- the
    program's 522 is 555 minus the errors the factors introduce -- and it is checked on the real run."""
    types = set(fr.TERMINAL_TYPES.split("|"))
    program = parse_file(ORG / "embryo_factors.bio")
    rules = [d for d in program.decisions if d.id.startswith("factor_fate_")]
    assert rules and all(set(d.when["cell_type"].split("|")) <= types for d in rules)

    ref = ReferenceLineage.load()
    terminal = {c.id: c.cell_type for c in fr.embryonic_terminal(ref)}
    body = Body(program, means=True).run(until=800.0)
    wrong = [cid for cid, want in terminal.items() if body.cells[cid].cell_type != want]
    for cid in terminal:
        cell = body.cells[cid]
        assert f"fate_{cid}" in cell.fired  # the lookup fired on all 555, whatever the factors then did
        if cid in wrong:
            assert any(d.startswith("factor_fate_") for d in cell.fired)  # every error is the factors'


@pytest.mark.skipif(not RESULT.exists(), reason="run scripts/celegans_fate_credit.py")
def test_the_ladder_and_its_denominators():
    r = json.loads(RESULT.read_text())
    ladder = r["answer"]["ladder"]
    assert r["answer"]["denominator"] == 555
    # the published number first, then one borrowing removed per rung, and each rung is lower
    assert [x["of_555"] for x in ladder] == sorted((x["of_555"] for x in ladder), reverse=True)
    assert ladder[0]["of_555"] == 522 and ladder[-1]["of_555"] == 0
    # nothing fitted anywhere is the honest headline, and it is below the factor-free baseline
    cited = r["arms"]["cited only, nothing fitted"]
    assert cited["honest"] < r["baselines_that_use_no_factors"]["majority_tissue"]["right"]
    # holding out rules that were never fitted must change nothing at all
    assert r["arms"]["cited only, held out (must be identical)"]["honest"] == cited["honest"]
    # the ceiling is stated: four tissues of seventeen have a cited rule
    assert (
        r["denominators"]["ceiling_from_the_rule_vocabulary"][
            "cited rules (intestine, muscle, hypoderm, sheath)"
        ]
        == 247
    )
    # F2: the accounting closes on every arm and the lookup is a clean fallback
    for a in r["arms"].values():
        assert a["claimed_right"] + a["claimed_wrong"] + a["abstained"] == 555
        assert a["abstained_wrong"] == 0
        assert a["claimed"] == a["lookup_wrote_it_first"]  # never a fate the lookup had not written
    # F4: the cited rules are reading the cell's factors, not the composition of the classes
    null = r["shuffle_null_for_the_cited_rules"]
    assert cited["honest"] > null["right_mean"] + 2 * null["right_sd"]


@pytest.mark.skipif(not RESULT.exists(), reason="run scripts/celegans_fate_credit.py")
def test_the_terminality_gate_is_measured_and_not_asserted():
    """The deepest borrowing: with the lookup's 555 fate decisions deleted the rules are never even asked,
    and when the guard is opened so that they can be, they call terminality wrong more often than right."""
    r = json.loads(RESULT.read_text())
    g = r["terminality_gate"]
    shipped, opened = g["guard as shipped"], g["guard opened to any undifferentiated cell"]
    assert shipped["lookup_fate_decisions_deleted"] == opened["lookup_fate_decisions_deleted"] == 555
    assert shipped["cells_a_factor_rule_fired_on"] == 0 and shipped["honest_of_555"] == 0
    assert opened["cells_a_factor_rule_fired_on"] > 0  # the arm is measuring something
    assert opened["of_those_not_terminal"] > opened["of_those_terminal"]
    # and the fates they do get right that way are a fraction of what they get with the gate in place
    assert opened["honest_of_555"] < r["arms"]["cited only, nothing fitted"]["honest"] / 3


@pytest.mark.skipif(not RESULT.exists(), reason="run scripts/celegans_fate_credit.py")
def test_the_pre_registration_was_not_rewritten_to_match_the_run():
    """Three of the seven predictions broke. The discipline is that the docstring keeps the numbers it was
    committed with (e2defff, before the instrument had been run) and the result file records what happened,
    with a diagnosis on each break -- not that the prediction is edited afterwards."""
    text = Path("scripts/celegans_fate_credit.py").read_text()
    assert "P2 (derived)  the whole shipped rule set, in sample, decides correctly 299 of 555" in text
    assert "P1 (derived)" in text and "130 to 160 of 555" in text
    pre = json.loads(RESULT.read_text())["preregistered"]
    assert len(pre) == 10  # seven predictions and three of the four void conditions
    for name, claim in pre.items():
        assert "observed" in claim and isinstance(claim["held"], bool)
        assert claim["held"] or claim.get("why_it_broke"), f"{name} broke with no diagnosis"
    # the claim the measurement turns on held, and the derived ones broke in the generous direction
    assert pre["P3 the honest number is below the factor-free majority-tissue baseline"]["held"]
    assert not pre["P1 cited rules alone decide 130 to 160 of 555 correctly"]["held"]
    assert pre["P1 cited rules alone decide 130 to 160 of 555 correctly"]["observed"] > 160
