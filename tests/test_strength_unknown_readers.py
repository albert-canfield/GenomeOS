# SPDX-License-Identifier: AGPL-3.0-or-later
"""ROADMAP §5 item 13, A: before a compiler emits `strength: unknown`, its readers carry it.

lane-engdebts (be2d00a) made `strength: unknown` a legal statement: `Module.active_rules` holds such a
rule at strength 0 and `strength_gate` names it with why. Two readers took `Rule.strength` as a number:
the Cell view's rule table (`Api.cell`), which put the bare UNKNOWN singleton into a JSON payload, and
the bridge (`attribution/bridge.py`), which compared it with floats and, through `active_rules`, read
the held copy's 0 as if the program had stated it. Each must show or carry the unknown with its reason,
never a number in its place and never a crash (lane-strengthread, 2026-10-10).
"""

from __future__ import annotations

import json

from genomeos.attribution import bridge
from genomeos.ir import STRENGTH_UNKNOWN_WHY, UNKNOWN
from genomeos.lang.parser import parse
from genomeos.web.server import Api

# ---- the web rule table ---------------------------------------------------------------------------

PROGRAM = """module t.unknown_strength
gene A { max: 10; basal: 1; produces: PA }
gene T { max: 10; basal: 1 }
gene U { max: 10; basal: 1 }
protein PA { evidence: curated "test"; confidence: 0.9 }
cell_type Here { expresses: A, T, U }
cell_type There { expresses: A, T, U }
rule PA activates T { strength: unknown; threshold: 5; when: cell_type = Here }
rule PA activates U { strength: 0.4; threshold: 5 }
"""


def _cell(tmp_path, cell_type: str) -> dict:
    (tmp_path / "data").mkdir(exist_ok=True)
    (tmp_path / "data" / "u.bio").write_text(PROGRAM)
    return Api(tmp_path).cell("data/u.bio", cell_type)


def test_the_rule_table_shows_an_unknown_strength_as_unknown_with_why(tmp_path):
    d = _cell(tmp_path, "Here")
    rows = {(r["source"], r["target"]): r for r in d["rules"]}
    held, known = rows[("PA", "T")], rows[("PA", "U")]

    assert held["strength"] == "unknown", "the word the program states, not 0.0 and not 1.0"
    assert (held["status"], held["reason"]) == ("held", STRENGTH_UNKNOWN_WHY)
    assert (known["status"], known["reason"], known["strength"]) == ("active", "", 0.4)
    assert d["counts"] == {"total": 3, "active": 2, "gated_out": 0, "silenced_out": 0, "held": 1}
    assert any("`strength: unknown`" in n for n in d["notes"])
    (edge,) = [e for e in d["graph"]["edges"] if e["target"] == "T"]
    assert (edge["status"], edge["strength"]) == ("held", "unknown")
    json.dumps(d)  # what Api._json does with it: the payload must serialise, not raise TypeError


def test_a_gated_out_unknown_rule_is_still_shown_unknown_and_not_counted_held(tmp_path):
    d = _cell(tmp_path, "There")
    (row,) = [r for r in d["rules"] if r["target"] == "T"]
    assert (row["status"], row["strength"]) == ("gated out", "unknown")
    assert "held" not in d["counts"], "a program with no held rule keeps the counts it always had"
    json.dumps(d)


# ---- the bridge -----------------------------------------------------------------------------------

K562 = {"cell_type": "K562"}
GENES = {"GENE1": {"basal_rate": 1.0, "max_rate": 10.0}}
PREDICTED = (
    "element E1 { class: enhancer }\ngene GENE1 { symbol: GENE1 }\n"
    'rule E1 activates GENE1 { strength: STRENGTH; when: cell_type = K562; evidence: predicted "AlphaGenome'
    ' deletion, K562" effect -1 log2 fold change on deletion, probability unavailable }\n'
)
MEASURED = (
    "element E1_measured { class: enhancer; evidence: experimental"
    ' "CRISPRi enhancer-gene screens" links GENE1 in K562 from 1 pairs }\n'
    "gene GENE1 { symbol: GENE1 }\n"
    "rule E1_measured activates GENE1 { strength: STRENGTH; when: cell_type = K562; evidence: experimental"
    ' "CRISPRi enhancer-gene screens, silenced in K562"; confidence: 0.9 }\n'
)


def _unknown_reported(b: bridge.Bridged, subject: str) -> None:
    assert [(i.reason, i.subject) for i in b.unresolved] == [("strength_unknown", subject)]
    assert STRENGTH_UNKNOWN_WHY in b.unresolved[0].detail
    assert b.strengths == {} and b.observations == {}
    assert "strength_unknown" in bridge.UNRESOLVED


def test_the_bridge_fits_the_same_programs_when_the_strength_is_stated():
    """The control: each fixture below is resolvable once a number is stated, so the unknown is the cause."""
    for program in (PREDICTED, MEASURED):
        b = bridge.parameterize(parse(program.replace("STRENGTH", "0.3"), "t"), K562, genes=GENES)
        assert b.unresolved == [] and len(b.strengths) == 1


def test_an_unknown_predicted_strength_is_not_overruled_by_its_note():
    """Through active_rules the bridge saw the held copy at 0 and fitted the note's -1 as if the rule ran."""
    m = parse(PREDICTED.replace("STRENGTH", "unknown"), "t")
    _unknown_reported(bridge.parameterize(m, K562, genes=GENES), "E1")


def test_an_unknown_measured_strength_is_not_read_as_a_zero_effect():
    """Read as 0 it became RHO = 1 + 0, a mechanism that 'moves nothing': a finding nobody measured."""
    m = parse(MEASURED.replace("STRENGTH", "unknown"), "t")
    _unknown_reported(bridge.parameterize(m, K562, genes=GENES), "E1_measured")
    assert not any(
        r.action.value == "activates" for r in bridge.parameterize(m, K562, genes=GENES).module.rules
    )


def test_observation_of_the_programs_own_rule_does_not_crash():
    """A caller passing module.rules directly gets the UNKNOWN itself: `0.0 < UNKNOWN` raised TypeError."""
    for program, subject in ((PREDICTED, "E1"), (MEASURED, "E1_measured")):
        m = parse(program.replace("STRENGTH", "unknown"), "t")
        (rule,) = [r for r in m.rules if r.action.value == "activates"]
        assert rule.strength is UNKNOWN
        o = bridge.observation_of(rule, m)
        assert (o.reason, o.subject) == ("strength_unknown", subject)
        assert rule.strength is UNKNOWN, "the program's statement is left as it was"


def test_an_unknown_mechanism_still_counts_against_identifiability():
    """Dropping it would fit the other element alone, as though the unknown one were absent."""
    other = (
        'rule E9 activates GENE1 { strength: 0.3; when: cell_type = K562; evidence: predicted "AlphaGenome'
        ' deletion, K562" effect -1 log2 fold change on deletion, probability unavailable }\n'
        "element E9 { class: enhancer }\n"
    )
    m = parse(PREDICTED.replace("STRENGTH", "unknown") + other, "t")
    b = bridge.parameterize(m, K562, genes=GENES)
    assert [(i.reason, i.subject) for i in b.unresolved] == [("not_identifiable", "GENE1")]
