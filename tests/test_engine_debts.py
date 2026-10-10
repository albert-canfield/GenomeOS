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


# ---- debt 4: an event's rate unit and a param's unit are checked ---------------------------------------


@pytest.mark.parametrize("unit", ["/yr", "/day", "1/yr", "/h", "/min", "1/wk", "/years"])
def test_an_event_rate_per_time_is_accepted(unit):
    (e,) = parse(f"module t\nevent divide {{ rate: 0.5 {unit} }}\n").events
    assert e.rate == 0.5 and e.rate_unit == unit


def test_an_event_rate_with_no_unit_keeps_its_default():
    (e,) = parse("module t\nevent divide { rate: 0.5 }\n").events
    assert e.rate_unit == "1/yr"


@pytest.mark.parametrize("unit", ["yr", "per_yr", "/fortnight", "/", "1/", "2/yr", "bp"])
def test_an_event_rate_that_is_not_per_time_is_refused_naming_the_unit_and_the_set(unit):
    with pytest.raises(BioLangError) as err:
        parse(f"module t\nevent divide {{ rate: 0.5 {unit} }}\n")
    msg = str(err.value)
    assert f"event 'divide': rate unit {unit!r} is not a rate" in msg
    assert "/min, /h, /d, /wk, /yr" in msg


@pytest.mark.parametrize(
    "unit",
    [
        "",
        "h",
        "min",
        "bp",
        "CpGs",
        "sites",
        "divisions",
        "fraction",
        "a.u.",
        "/division",
        "cells/h",
        "nM",
        "fL",
    ]
    + ["bp/division", "1/h", "/yr", "molecules/min"],
)
def test_the_param_units_the_programs_use_and_their_kin_are_accepted(unit):
    m = parse(f"module t\nparam x = 2 {unit} {{ }}\n")
    assert m.parameters["x"].unit == unit


@pytest.mark.parametrize("unit", ["hrs", "bps", "per_division", "cells/fortnight", "a.u", "%", "kg"])
def test_a_param_unit_outside_the_set_is_refused_naming_the_unit_and_the_set(unit):
    with pytest.raises(BioLangError) as err:
        parse(f"module t\nparam x = 2 {unit} {{ }}\n")
    msg = str(err.value)
    assert f"param 'x': unit {unit!r} is not one a param may state" in msg
    assert "bp, nt, aa, CpGs" in msg and "division" in msg


def test_an_unknown_param_value_still_takes_a_checked_unit():
    assert parse("module t\nparam x = unknown /division { }\n").parameters["x"].unit == "/division"
    with pytest.raises(BioLangError, match="unit '/divisions' is not one"):
        parse("module t\nparam x = unknown /divisions { }\n")


def test_a_half_life_the_runtimes_read_in_hours_refuses_another_unit():
    """grn.py and located.py read these two in hours and convert nothing, so `30 min` ran as 30 h."""
    for name in ("mrna_half_life", "protein_half_life"):
        for ok in ("", "h", "hours"):
            assert parse(f"module t\nparam {name} = 0.5 {ok} {{ }}\n").parameters[name].value == 0.5
        with pytest.raises(BioLangError, match=f"param '{name}' has unit 'min'; .* in hours"):
            parse(f"module t\nparam {name} = 30 min {{ }}\n")


# ---- debt 5: v0.4 §2's confidence cap is reported, not applied -------------------------------------------


def test_the_cap_names_blocks_stated_above_it_on_a_construct_not_past_its_gate():
    from genomeos.lang.tools import CONFIDENCE_CAP, UNGATED, cap_warnings

    assert CONFIDENCE_CAP == 0.3
    assert set(UNGATED) == {"pool", "allocation", "cost", "partition", "integrated read"}
    src = """module t
compartment Cell { volume: 1.0; absolute_volume: 1 pL }
pool ATP { location: Cell; size: 100; confidence: 0.6 }
pool GTP { location: Cell; size: 100; confidence: 0.3 }
pool UTP { location: Cell; size: 100 }
allocation share { pool: ATP; policy: proportional; confidence: 0.5 }
protein P { location: Cell; cost: ATP 1 per chain; confidence: 0.9 }
event divide { rate: 1 /day; partition: contents; confidence: 0.4 }
"""
    warn = cap_warnings(parse(src))
    assert warn == [
        ("pool ATP", "pool", 0.6),
        ("allocation share", "allocation", 0.5),
        ("protein P", "cost", 0.9),
        ("event divide", "partition", 0.4),
    ], "exactly at the cap, or stating no confidence, is not above it"


def test_an_integrated_read_in_any_clause_is_found():
    from genomeos.lang.tools import cap_warnings

    src = """module t
cell_type A { }
decision d1 { action: differentiate; when: ELT-2.exposure(lineage) >= 15; to: A; confidence: 0.8 }
decision d2 { action: differentiate; when: ELT-2 >= 15; to: A; confidence: 0.8 }
"""
    assert cap_warnings(parse(src)) == [("decision d1", "integrated read", 0.8)]


def test_check_reports_the_cap_and_changes_no_confidence():
    from genomeos.lang import parse_file

    m = parse_file("data/demo/stage2_pools.bio")
    before = {k: e.confidence for k, e in m.entities.items()}
    text = check_module(m)
    assert (
        "confidence above v0.4 §2's cap of 0.3 on a construct not past its gate (reported, not applied): 5"
        in text
    )
    assert "pool RNAPII: 0.6 (pool)" in text
    assert {k: e.confidence for k, e in m.entities.items()} == before
    assert "cap of 0.3" not in check_module(parse_file("data/demo/repressilator.bio"))


# ---- debt 3: version labels ---------------------------------------------------------------------------


def test_the_grammar_and_bioir_carry_one_label_until_albert_moves_it():
    """0.4 with v0.5's reader in: the number is Albert's (BIOLANG-v0.5-READER.md, open decision 3).
    Held equal so that neither label moves without the other (lane-engdebts, 2026-10-10)."""
    from genomeos.ir import BIOIR_VERSION
    from genomeos.lang import grammar
    from genomeos.lang import parser as parser_module

    assert grammar.VERSION == BIOIR_VERSION == "0.4"
    assert parse("module t\ngene A { }\n").to_dict()["bioir_version"] == BIOIR_VERSION
    assert "v0.3 parser" not in (parser_module.__doc__ or "")
    assert "Version label 0.4" in parser_module.__doc__ and "open decision 3" in parser_module.__doc__
    assert "Albert's decision" in grammar.render()
