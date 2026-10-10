# SPDX-License-Identifier: Apache-2.0
"""ROADMAP §5 item 13, A: the engine's debts (lane-engdebts, 2026-10-10).

`strength: unknown` holds a rule instead of running it at 1.0."""

import pytest

from genomeos.bio import evaluate, parse_tests
from genomeos.ir import STRENGTH_UNKNOWN_WHY, UNKNOWN, Module
from genomeos.lang import parse
from genomeos.lang.parser import BioLangError
from genomeos.lang.tools import check_module, run_module
from genomeos.runtime import NetworkRuntime
from genomeos.runtime.located import LocatedRuntime

# ---- debt 1: `strength: unknown` -------------------------------------------------------------------

NET = """module t
gene ACT { max: 40; basal: 1; produces: Act }
gene TGT { max: 60; basal: 0.05; produces: Tgt }
protein Act { half_life: 2 }
protein Tgt { half_life: 2 }
"""


def _final(module: Module, species: str = "TGT.mRNA", ctx: dict | None = None) -> float:
    return NetworkRuntime(module, ctx or {}, seed=0).run(24, dt=0.05).final()[species]


def _reg(m: Module) -> list:
    """The regulatory rules only: each `gene { produces: P }` adds a produces rule of its own."""
    return [r for r in m.rules if r.action.value in ("activates", "inhibits")]


def _with(strength_line: str) -> Module:
    return parse(NET + f"rule Act activates TGT {{ threshold: 5{strength_line} }}\n")


def test_strength_unknown_parses_to_unknown_in_either_case():
    for word in ("unknown", "UNKNOWN", "Unknown"):
        (r,) = _reg(_with(f"; strength: {word}"))
        assert r.strength is UNKNOWN


def test_no_strength_line_still_runs_at_one_and_a_stated_strength_is_unchanged():
    """The decision for a rule with NO strength line: 1.0, as before. Inline programs rely on it."""
    assert _reg(_with(""))[0].strength == 1.0
    assert _reg(_with("; strength: 0.4"))[0].strength == 0.4
    for m in (_with(""), _with("; strength: 0.4")):
        assert m.strength_gate({}) == [] and m.held_rules({}) == []
        assert all(a is b for a, b in zip(m.active_rules({}), m.rules, strict=True))


def test_an_unknown_strength_is_held_at_zero_and_named_with_why():
    m = _with("; strength: unknown")
    (rule,) = _reg(m)
    (held,) = _reg(Module(name="x", rules=m.active_rules({})))
    assert held.strength == 0.0 and held.id == rule.id
    assert rule.strength is UNKNOWN, "the hold is a copy; the program's statement stays"
    assert m.strength_gate({}) == [(rule, STRENGTH_UNKNOWN_WHY)]
    assert m.held_rules({}) == [rule]


def test_a_held_unknown_activator_behaves_as_strength_zero_not_as_one_and_not_as_dropped():
    unknown, zero, unstated = _with("; strength: unknown"), _with("; strength: 0"), _with("")
    assert _final(unknown) == pytest.approx(_final(zero))
    assert _final(unknown) < 1, "held: the target stays at basal"
    assert _final(unstated) > 50, "a rule with no strength line still runs at 1.0"
    dropped = parse(NET)
    assert _final(dropped) > 50, "dropping the rule would drive the target to max: why it is held"


def test_a_held_unknown_inhibitor_lets_its_target_through():
    m = parse(NET + "rule Act inhibits TGT { threshold: 5; strength: unknown }\n")
    free = parse(NET + "rule Act inhibits TGT { threshold: 5; strength: 0 }\n")
    assert _final(m) == pytest.approx(_final(free))


def test_unknown_strength_is_refused_where_a_strength_carries_no_effect():
    with pytest.raises(BioLangError, match="cannot state `strength: unknown`"):
        parse(NET + "rule ACT produces Act { strength: unknown }\n")
    with pytest.raises(BioLangError, match=r"strength \(a number or `unknown`\) expects a number"):
        parse(NET + "rule Act activates TGT { strength: high }\n")


def test_unknown_strength_survives_the_json_round_trip():
    m = _with("; strength: unknown")
    d = m.to_dict()
    assert [r["strength"] for r in d["rules"] if r["action"] == "activates"] == ["UNKNOWN"]
    back = Module.from_dict(d)
    assert _reg(back)[0].strength is UNKNOWN
    assert _reg(Module(name="x", rules=back.active_rules({})))[0].strength == 0.0
    assert _reg(Module.from_dict(_with("; strength: 0.4").to_dict()))[0].strength == 0.4


def test_the_reader_and_an_unknown_strength_hold_a_rule_once():
    src = (
        NET
        + "element ENH1 { class: enhancer; targets: TGT }\ncell_type Closed { }\n"
        + "reader Closed { closed: ENH1 }\n"
        + "rule Act activates TGT { through: ENH1; threshold: 5; strength: unknown }\n"
        + "rule Act inhibits TGT { threshold: 5; strength: unknown; id: second }\n"
    )
    m = parse(src)
    ctx = {"cell_type": "Closed"}
    assert len(m.held_rules(ctx)) == 2
    assert [r.strength for r in _reg(Module(name="x", rules=m.active_rules(ctx)))] == [0.0, 0.0]
    (line,) = evaluate(m, parse_tests("# test: held_rules == 2 in cell_type=Closed\n"))
    assert line["ok"], line


def test_check_and_run_name_the_held_rules_with_the_reason():
    m = _with("; strength: unknown")
    text = check_module(m)
    assert "strength unknown: 1 rules held at strength 0 in every cell" in text
    assert STRENGTH_UNKNOWN_WHY in text and "held Act activates TGT: strength unknown" in text
    assert "strength unknown: 1 rules held at strength 0;" in check_module(m, {"cell_type": "X"})
    assert "held Act activates TGT: strength unknown" in run_module(m, hours=2)[0]
    assert "strength unknown" not in check_module(_with(""))


def test_the_located_runtime_holds_an_unknown_strength_too():
    cell = """module loc
compartment Cytosol { translation: yes; volume: 1.0 }
protein Rp { location: Cytosol; initial: 20000; half_life: 1e6 h }
gene T { location: Cytosol; basal: 1; max: 1000; produces: Tp }
protein Tp { location: Cytosol }
rule Rp inhibits T { threshold: 100; hill: 8; strength: STRENGTH }
"""

    def final(strength: str) -> float:
        return (
            LocatedRuntime(parse(cell.replace("STRENGTH", strength)))
            .run(hours=24.0, dt=0.05)
            .final("Tp@Cytosol")
        )

    assert final("unknown") == pytest.approx(final("0"))
    assert final("unknown") > 10 * final("1.0"), "held, not run at full strength"
