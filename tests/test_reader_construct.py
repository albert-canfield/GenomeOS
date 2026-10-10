# SPDX-License-Identifier: Apache-2.0
"""BioLang v0.5, the reader: `reader <CellType> { open/poised/closed/unknown }` and a rule's
`through: E` (docs/BIOLANG-v0.5-READER.md). A rule acting through an element runs only in a cell whose
reader states the element open; anywhere else it is held at strength 0, `unknown` included."""

import dataclasses
import json
from pathlib import Path

import pytest

from genomeos.bio import CheckLine, evaluate, parse_tests
from genomeos.ir import READER_STATES, Module, Reader, Rule
from genomeos.lang import parse, parse_file
from genomeos.lang.parser import BioLangError
from genomeos.lang.tools import check_module
from genomeos.runtime import NetworkRuntime

DEMO = Path("data/demo/reader_gating.bio")

BASE = """module t
gene ACT { max: 40; basal: 1; produces: Act }
gene TGT { max: 60; basal: 0.05; produces: Tgt }
protein Act { half_life: 2 }
protein Tgt { half_life: 2 }
element ENH1 { class: enhancer; targets: TGT }
region R1 { locus: chr1:100-200 }
cell_type Open { name: open cell }
cell_type Closed { }
cell_type Poised { }
cell_type Stated { }
cell_type Silent { }
cell_type Unread { }
reader Open { open: ENH1; evidence: experimental "DNase-seq"; confidence: 0.8 }
reader Closed { closed: ENH1 }
reader Poised { poised: ENH1 }
reader Stated { unknown: ENH1 }
reader Silent { open: R1 }
"""
GATED = BASE + "rule Act activates TGT { through: ENH1; threshold: 5 }\n"


def _gated(module: Module, cell: str) -> Rule:
    ctx = {"cell_type": cell} if cell else {}
    return next(r for r in module.active_rules(ctx) if r.through)


def _final(module: Module, cell: str, species: str) -> float:
    return NetworkRuntime(module, {"cell_type": cell}, seed=0).run(24, dt=0.05).final()[species]


# ---- the IR ----------------------------------------------------------------------------------------


def test_a_reader_and_a_through_rule_compile_to_the_ir():
    m = parse(GATED)
    assert [rd.cell_type for rd in m.readers] == ["Open", "Closed", "Poised", "Stated", "Silent"]
    assert m.readers[0].states == {"ENH1": "open"}
    assert m.readers[0].evidence.source == "DNase-seq" and m.readers[0].confidence == 0.8
    assert m.readers[4].states == {"R1": "open"}, "a region is readable as an element is"
    (rule,) = [r for r in m.rules if r.through]
    assert rule.through == "ENH1" and rule.id == "Act activates TGT"
    assert set(READER_STATES) == {"open", "poised", "closed", "unknown"}


def test_readers_and_through_survive_the_json_round_trip():
    m = parse(GATED)
    back = Module.from_dict(json.loads(json.dumps(m.to_dict())))
    assert back.readers == m.readers
    assert [r.through for r in back.rules] == [r.through for r in m.rules]
    for cell in ("Open", "Closed", "Poised", "Stated", "Silent", "Unread", ""):
        ctx = {"cell_type": cell} if cell else {}
        assert [(r.id, r.strength) for r in back.active_rules(ctx)] == [
            (r.id, r.strength) for r in m.active_rules(ctx)
        ]


def test_bioir_written_before_the_reader_reads_back_with_none():
    d = parse(BASE).to_dict()
    d.pop("readers")
    for r in d["rules"]:
        r.pop("through")
    back = Module.from_dict(d)
    assert back.readers == [] and all(r.through == "" for r in back.rules)


def test_module_merge_carries_the_readers():
    m = Module(name="host")
    m.merge(parse(GATED))
    assert len(m.readers) == 5


# ---- the gate ---------------------------------------------------------------------------------------


def test_open_runs_the_rule_as_written():
    m = parse(GATED)
    rule = _gated(m, "Open")
    assert rule.strength == 1.0
    assert rule is next(r for r in m.rules if r.through), "an open rule is the program's own rule"
    ((_, state, why),) = m.reader_gate({"cell_type": "Open"})
    assert state == "open" and "experimental DNase-seq" in why


@pytest.mark.parametrize(
    ("cell", "state", "reason"),
    [
        ("Closed", "closed", "states ENH1 closed"),
        ("Poised", "poised", "states ENH1 poised"),
        ("Stated", "unknown", "states ENH1 unknown"),
        ("Silent", "unknown", "the reader of Silent does not state ENH1"),
        ("Unread", "unknown", "no reader for cell_type Unread"),
        ("", "unknown", "the context names no cell_type"),
    ],
)
def test_every_state_but_open_holds_the_rule_at_strength_zero_and_says_why(cell, state, reason):
    m = parse(GATED)
    ctx = {"cell_type": cell} if cell else {}
    assert _gated(m, cell).strength == 0.0
    ((rule, got, why),) = m.reader_gate(ctx)
    assert (rule.id, got) == ("Act activates TGT", state)
    assert reason in why
    assert next(r for r in m.rules if r.through).strength == 1.0, "holding copies; the program is unchanged"


def test_a_context_may_name_the_cell_type_by_its_name():
    m = parse(GATED)
    assert _gated(m, "open cell").strength == 1.0


def test_a_held_activator_leaves_its_target_at_basal_and_dropping_it_would_drive_it_to_max():
    """Why held is strength 0 and not removal: the network runtime transcribes a gene with no
    activator rule at its full max_rate, so a dropped rule would read as the open cell."""
    m = parse(GATED)
    assert _final(m, "Open", "TGT.mRNA") > 50
    assert _final(m, "Closed", "TGT.mRNA") < 1
    assert _final(m, "Unread", "TGT.mRNA") < 1, "unknown is never silently open"
    dropped = dataclasses.replace(m, rules=[r for r in m.rules if not r.through])
    assert _final(dropped, "Closed", "TGT.mRNA") > 50, "the trap the held copy avoids"


def test_a_held_inhibitor_lets_its_target_through():
    m = parse(BASE + "rule Act inhibits TGT { through: ENH1; threshold: 5; strength: 0.9 }\n")
    assert _final(m, "Closed", "TGT.mRNA") > 50
    assert _final(m, "Open", "TGT.mRNA") < 10


def test_the_reader_gates_only_after_when_and_expresses():
    m = parse(GATED.replace("threshold: 5 }", "threshold: 5; when: cell_type = Closed }"))
    assert m.reader_gate({"cell_type": "Open"}) == [], "a rule `when` drops is not the reader's to hold"
    assert [s for _, s, _ in m.reader_gate({"cell_type": "Closed"})] == ["closed"]


def test_rules_without_through_are_selected_exactly_as_before():
    for path in ("data/demo/cell_context.bio", "data/demo/repressilator.bio"):
        m = parse_file(path)
        assert not m.readers and not any(r.through for r in m.rules)
        for ctx in ({}, {"cell_type": "Cardiomyocyte2"}, {"cell_type": "Neuron2"}):
            silenced = m.silenced_genes(ctx)
            before = [
                r for r in m.rules if r.applies(ctx) and r.source not in silenced and r.target not in silenced
            ]
            after = m.active_rules(ctx)
            assert len(after) == len(before) and all(a is b for a, b in zip(after, before, strict=True))


def test_check_lists_held_rules_by_state_and_never_folds_unknown_into_closed():
    m = parse(GATED)
    text = check_module(m, {"cell_type": "Unread"})
    assert "held at strength 0: poised 0, closed 0, unknown 1" in text
    assert "no reader for cell_type Unread" in text
    assert "readers: 5 for 5 cell types" in check_module(m)


# ---- what the parser refuses -------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("extra", "message"),
    [
        ("reader Open Closed { open: ENH1 }", "reader header must be one cell_type id"),
        ("reader Unread { }", "states no element"),
        ("reader Unread { open: ENH1, ENH1 }", "states 'ENH1' twice"),
        ("reader Unread { open: ENH1; closed: ENH1 }", "both open and closed"),
        ("reader Unread { silent: ENH1 }", "has no key 'silent'"),
        ("reader Nowhere { open: ENH1 }", "undeclared cell_type 'Nowhere'"),
        ("reader Unread { open: TGT }", "'TGT', which is not a declared element or region"),
        ("reader Open { closed: ENH1 }", "readers disagree"),
        ("rule Act activates TGT { through: ENH1, R1 }", "names exactly one element"),
        ("rule Act activates TGT { through: NOPE }", "acts through 'NOPE'"),
        ("rule ACT produces Tgt { through: ENH1 }", "cannot act `through`"),
    ],
)
def test_malformed_readers_and_through_are_parse_errors(extra, message):
    with pytest.raises(BioLangError, match=message):
        parse(BASE + extra + "\n")


def test_two_readers_of_one_cell_type_that_agree_are_one_reading():
    m = parse(BASE + 'reader Open { open: ENH1, R1; evidence: curated "a second assay" }\n')
    assert {e: s for e, (s, _) in m.reader_states("Open").items()} == {"ENH1": "open", "R1": "open"}


def test_a_located_program_cannot_use_through_yet():
    src = "module loc\ncompartment Cell { volume: 1.0 }\n" + GATED.split("\n", 1)[1]
    with pytest.raises(BioLangError, match="located program cannot use `through`"):
        parse(src)


# ---- bio test in a context, and the demo ---------------------------------------------------------------


def test_test_lines_carry_a_context_and_stay_four_tuples():
    lines = parse_tests("# test: TGT.mRNA final > 50 in cell_type=Open, stage=late\n# test: rules == 3\n")
    assert lines[0] == ("TGT.mRNA", "final", ">", 50.0) and isinstance(lines[0], CheckLine)
    assert lines[0].context == {"cell_type": "Open", "stage": "late"}
    assert lines[1].context == {}


def test_the_demo_reads_the_same_rules_differently_in_five_cells():
    text = DEMO.read_text()
    results = evaluate(parse(text), parse_tests(text), hours=48)
    assert results and all(r["ok"] for r in results), [r for r in results if not r["ok"]]
    got = {r["test"]: r["got"] for r in results}
    assert got["TGT.mRNA final > 50 in cell_type=OpenCell"] > 50
    assert got["TGT.mRNA final < 1 in cell_type=ClosedCell"] < 1
    assert got["TGT.mRNA final < 1 in cell_type=UnreadCell"] < 1
    assert sum(1 for r in results if "in cell_type=" in r["test"]) >= 12


def test_a_failing_context_line_fails():
    m = parse(DEMO.read_text())
    (bad,) = evaluate(m, parse_tests("# test: held_rules == 0 in cell_type=ClosedCell\n"))
    assert not bad["ok"] and bad["got"] == 2


def test_the_reader_dataclass_is_exported():
    assert Reader(cell_type="X").states == {}
