# SPDX-License-Identifier: Apache-2.0
"""One `when` matcher for every block that has a `when` (lane-when, 2026-09-28).

Decisions and timers read `when` with `genomeos.ir.model.matches`; rules and events compared by
plain equality, so `a|b`, `absent` and the comparisons the grammar allows compiled on a rule or an
event and never matched: `when: cell_type = K562|HepG2` gave a rule that ran in no cell and said
nothing. These tests pin each construct for rules and events, with R1's `unknown` (matches no
context) kept, and check the four blocks agree clause for clause.
"""

from __future__ import annotations

import pytest

from genomeos.ir.model import Action, Decision, Event, Rule, Timer, matches
from genomeos.lang import parse


def rule(when):
    return Rule(id="r", source="A", action=Action.ACTIVATE, target="B", when=when)


def event(when):
    return Event(id="e", rate=1.0, when=when)


MAKERS = {"rule": rule, "event": event}

# (when, context, expected) for every construct docs/BIOLANG-GRAMMAR.md allows in a `when`
CASES = [
    # alternatives
    ({"cell_type": "K562|HepG2"}, {"cell_type": "K562"}, True),
    ({"cell_type": "K562|HepG2"}, {"cell_type": "HepG2"}, True),
    ({"cell_type": "K562|HepG2"}, {"cell_type": "WTC11"}, False),
    ({"cell_type": "K562|HepG2"}, {}, False),
    # absent
    ({"host": "absent"}, {}, True),
    ({"host": "absent"}, {"cell_type": "K562"}, True),
    ({"host": "absent"}, {"host": "ecoli"}, False),
    # the four comparisons, both sides of each boundary
    ({"generation": ">=3"}, {"generation": "3"}, True),
    ({"generation": ">=3"}, {"generation": "2"}, False),
    ({"generation": "<=3"}, {"generation": "3"}, True),
    ({"generation": "<=3"}, {"generation": "4"}, False),
    ({"generation": ">3"}, {"generation": "4"}, True),
    ({"generation": ">3"}, {"generation": "3"}, False),
    ({"generation": "<3"}, {"generation": "2.5"}, True),
    ({"generation": "<3"}, {"generation": "3"}, False),
    ({"generation": ">=3"}, {}, False),  # a comparison on a missing key does not hold
    ({"generation": ">=3"}, {"generation": "early"}, False),  # nor on a value that is not a number
    # unknown: R1, matches no context, including one that says "unknown" or "any"
    ({"cell_type": "unknown"}, {}, False),
    ({"cell_type": "unknown"}, {"cell_type": "K562"}, False),
    ({"cell_type": "unknown"}, {"cell_type": "unknown"}, False),
    ({"cell_type": "unknown"}, {"cell_type": "any"}, False),
    # any and plain equality, unchanged
    ({"cell_type": "any"}, {}, True),
    ({"cell_type": "any"}, {"cell_type": "K562"}, True),
    ({"cell_type": "K562"}, {"cell_type": "K562"}, True),
    ({"cell_type": "K562"}, {"cell_type": "HepG2"}, False),
    ({"cell_type": "K562"}, {}, False),
    ({}, {}, True),
    # clauses conjoin
    ({"cell_type": "K562|HepG2", "generation": ">=3"}, {"cell_type": "K562", "generation": "4"}, True),
    ({"cell_type": "K562|HepG2", "generation": ">=3"}, {"cell_type": "K562", "generation": "1"}, False),
]


@pytest.mark.parametrize("kind", sorted(MAKERS))
@pytest.mark.parametrize("when,context,expected", CASES)
def test_rules_and_events_read_every_when_construct(kind, when, context, expected):
    assert MAKERS[kind](when).applies(context) is expected


@pytest.mark.parametrize("when,context,expected", CASES)
def test_rule_event_decision_and_timer_agree_clause_for_clause(when, context, expected):
    got = {
        "matches": matches(when, context),
        "rule": rule(when).applies(context),
        "event": event(when).applies(context),
        "decision": Decision(id="d", when=when).applies(context),
        "timer": Timer(name="t", duration=1.0, when=when).applies(context),
    }
    assert got == dict.fromkeys(got, expected)


def test_a_rule_written_with_alternatives_fires_through_the_parser_and_active_rules():
    """The defect as it was found: this rule compiled, and ran in no cell."""
    module = parse(
        "module t\n"
        'gene A { symbol: A; evidence: curated "x"; confidence: 0.9 }\n'
        'gene B { symbol: B; evidence: curated "x"; confidence: 0.9 }\n'
        "rule A activates B {\n"
        "  when: cell_type = K562|HepG2, generation >= 3\n"
        '  evidence: curated "x"; confidence: 0.5\n'
        "}\n"
    )
    (r,) = module.rules
    assert r.when == {"cell_type": "K562|HepG2", "generation": ">=3"}
    assert module.active_rules({"cell_type": "HepG2", "generation": "3"}) == [r]
    assert module.active_rules({"cell_type": "WTC11", "generation": "3"}) == []
    assert module.active_rules({"cell_type": "K562", "generation": "2"}) == []


def test_an_event_written_with_absent_holds_only_while_the_key_is_missing():
    module = parse(
        "module t\nevent repair {\n  rate: 2 /yr\n  when: senescent = absent\n  effect: damage -= 1\n}\n"
    )
    (e,) = module.events
    assert e.applies({}) and not e.applies({"senescent": "1"})
