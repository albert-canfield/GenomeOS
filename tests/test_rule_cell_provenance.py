# SPDX-License-Identifier: AGPL-3.0-or-later
"""The cell-provenance proposal tested against the code and the committed corpus.

The proposal and its eight predictions were registered at `ea3ef2e`, in a commit of its own, before
any of this existed. Nothing here edits a `.bio` file, the compiler, the IR, the grammar or the
parser: where the proposal would need such an edit, the test demonstrates the REFUSAL the code gives
today rather than removing it.

Every test that reads the committed corpus goes through `tests.committed_data`, so a missing
committed program FAILS and names itself. A `skipif` on a tracked path would be the defect `4f44dbf`
was written against.
"""

from __future__ import annotations

import ast
import re
from functools import cache
from pathlib import Path

import pytest

from genomeos.attribution.measured import CONTEXT_UNKNOWN, SOURCES
from genomeos.ir.model import Rule
from genomeos.lang import rule_cell_provenance as cp
from genomeos.lang import rule_evidence_tier as tier
from genomeos.lang.parser import BioLangError, parse, parse_file
from genomeos.runtime.grn import NetworkRuntime
from tests.committed_data import ROOT, committed, must_be_committed

COMPILED = "data/organisms/human/noncoding_chr21.bio"
RUNTIME_FIXTURE = "data/demo/cell_context.bio"
COMPILER = "genomeos/attribution/compile.py"

#: The hand-authored programs git tracks that carry at least one `rule`. Named here rather than
#: globbed so a test cannot quietly lose its population to a renamed file.
HAND_AUTHORED = (
    "data/demo/cell_context.bio",
    "data/demo/concentration_threshold.bio",
    "data/demo/gastrulation.bio",
    "data/demo/lateral_inhibition.bio",
    "data/demo/repressilator.bio",
    "data/demo/segmentation_clock.bio",
    "data/organisms/celegans/acvu.bio",
    "data/organisms/human/erythrocyte.bio",
    "data/organisms/human/oxphos.bio",
)

#: The ONE adjudication this lane makes by hand, declared as an adjudication. `provenance_of` raises
#: on this rule because no signature the compiler writes matches it; the cell is asserted on the rule
#: line beside a citation of the factor, and nothing in the program ties that citation to the cell.
ADJUDICATED = {"Nkx2-5 activates MYH6": cp.AUTHOR_DECLARED}


@cache
def _cost():
    """The cost report as a module. Imported by path so the test does not depend on `scripts` being
    an importable package."""
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "rule_cell_provenance_cost", ROOT / "scripts/rule_cell_provenance_cost.py"
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@cache
def _parsed(relative: str):
    return parse_file(must_be_committed(relative))


def _cell_gated_rule(relative: str):
    """The program's one rule that names a cell. `rules[0]` is a rule the parser synthesised from a
    `gene ... { produces: ... }` clause and it names none, which is prediction P8."""
    return next(r for r in _parsed(relative).rules if cp.names_a_cell(r))


def _written_rule_headers(relative: str) -> set[str]:
    """The rule headers a program's TEXT writes, so a synthesised rule can be told from a written one."""
    out = set()
    for line in must_be_committed(relative).read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("rule "):
            out.add(stripped[len("rule ") :].split("{")[0].strip())
    return out


# --- the axis, its default and its refusal to be ordered --------------------------------------------


def test_the_absent_default_is_the_siblings_object_and_not_a_second_spelling() -> None:
    assert cp.NOT_ASSESSED is tier.NOT_ASSESSED
    assert cp.registration()["not_assessed"] is tier.NOT_ASSESSED


def test_no_assignment_in_the_module_binds_the_default_as_a_literal() -> None:
    """`is` on an identifier-shaped string is satisfied by interning, so it is not the evidence.

    The evidence is that the import is the ONLY way the value arrives: no assignment in the module
    binds the string itself. The word appears in prose throughout the file, and prose is not a
    second declaration, so this walks the syntax tree instead of grepping.
    """
    assert cp.binds_the_default_as_a_literal() == []


def test_a_planted_literal_assignment_is_found_so_the_walk_is_not_vacuous() -> None:
    source = f'X = {cp.NOT_ASSESSED!r}\nY: str = "something else"\n'
    found = [
        node.lineno
        for node in ast.walk(ast.parse(source))
        if isinstance(node, (ast.Assign, ast.AnnAssign))
        and isinstance(node.value, ast.Constant)
        and node.value.value == cp.NOT_ASSESSED
    ]
    assert found == [1]


def test_the_default_is_not_a_class_and_the_axis_is_exactly_one_wider() -> None:
    assert not cp.is_class(cp.NOT_ASSESSED)
    assert all(cp.is_class(c) for c in cp.CLASSES)
    assert len(cp.AXIS_VALUES) == len(cp.CLASSES) + 1
    assert cp.NOT_ASSESSED not in cp.CLASSES


@pytest.mark.parametrize("value", [*cp.CLASSES, cp.NOT_ASSESSED, "nonsense"])
def test_rank_raises_for_every_value_including_the_classes(value: str) -> None:
    """The sibling ranks its tiers and raises only on its default. This axis has no order at all."""
    with pytest.raises(TypeError, match="not ordered"):
        cp.rank(value)


def test_there_is_no_floor_over_this_axis() -> None:
    assert not hasattr(cp, "floor")
    assert tier.floor([tier.TIERS[0]]) == tier.TIERS[0]  # the sibling has one; this axis does not


def test_the_set_a_population_rests_on_keeps_the_default_in_it() -> None:
    got = cp.provenances_relied_on(
        [cp.ARGMAX_OF_PREDICTED_EFFECT, cp.NOT_ASSESSED, cp.ARGMAX_OF_PREDICTED_EFFECT]
    )
    assert cp.NOT_ASSESSED in got
    assert got == sorted({cp.ARGMAX_OF_PREDICTED_EFFECT, cp.NOT_ASSESSED})


def test_a_value_off_the_axis_cannot_enter_the_set() -> None:
    with pytest.raises(ValueError, match="not on the axis"):
        cp.provenances_relied_on([cp.ARGMAX_OF_PREDICTED_EFFECT, "measured"])


# --- the phrase, and the refusal that keeps it honest ------------------------------------------------


@pytest.mark.parametrize("bad", cp.FORBIDDEN_OF_THE_ARGMAX_CLASS)
def test_an_argmax_record_that_claims_an_observation_is_refused(bad: str) -> None:
    rec = {"cell_provenance": cp.ARGMAX_OF_PREDICTED_EFFECT, "note": f"K562, {bad} this cell line"}
    with pytest.raises(ValueError, match="turns a selection"):
        cp.check_no_observation_claim([rec])


def test_the_refusal_is_class_specific_so_a_screen_may_say_it_measured_something() -> None:
    cp.check_no_observation_claim(
        [{"cell_provenance": cp.MEASURED_PERTURBATION_IN_THAT_CELL, "note": "measured in K562"}]
    )


def test_the_registered_argmax_phrase_makes_no_observation_claim_of_its_own() -> None:
    cp.check_no_observation_claim(
        [
            {
                "cell_provenance": cp.ARGMAX_OF_PREDICTED_EFFECT,
                "phrase": cp.PHRASES[cp.ARGMAX_OF_PREDICTED_EFFECT],
            }
        ]
    )
    phrase = cp.PHRASES[cp.ARGMAX_OF_PREDICTED_EFFECT]
    assert "biosample name" in phrase  # the metadata copy carries no accession, so a NAME is the most
    assert "selection" in phrase


# --- P7: both refusals the storage decision rests on are real ---------------------------------------


def _bare_rule() -> Rule:
    from genomeos.ir.model import Action

    return Rule(id="r", source="A", action=Action.ACTIVATE, target="B")


def test_the_mark_cannot_be_set_on_a_rule_because_the_dataclass_has_slots() -> None:
    """P7a. `ir.model.Rule` is `@dataclass(slots=True)`, so the proposed field cannot be bolted on."""
    rule = _bare_rule()
    assert not hasattr(rule, "__dict__")
    with pytest.raises(AttributeError):
        rule.cell_provenance = cp.ARGMAX_OF_PREDICTED_EFFECT  # type: ignore[attr-defined]


def test_a_program_that_writes_the_property_on_a_rule_does_not_parse_today() -> None:
    """P7b. `parser._check_keys` refuses a rule key the grammar table does not name."""
    source = (
        "module t\n"
        "gene A { max: 10; produces: Pa; evidence: curated 'x' }\n"
        "gene B { max: 10; produces: Pb; evidence: curated 'x' }\n"
        "protein Pa { evidence: curated 'x' }\n"
        "protein Pb { evidence: curated 'x' }\n"
        "rule Pa activates B { strength: 1.0; "
        f"cell_provenance: {cp.ARGMAX_OF_PREDICTED_EFFECT}; evidence: curated 'x' }}\n"
    )
    with pytest.raises(BioLangError, match="no key 'cell_provenance'"):
        parse(source)


def test_the_property_is_absent_from_the_grammar_table_so_adoption_needs_that_edit() -> None:
    from genomeos.lang import grammar

    assert "cell_provenance" not in grammar.BLOCKS["rule"]["props"]
    assert "context_evidence" in grammar.BLOCKS["rule"]["props"]  # the precedent IS there


# --- P5: the runtime does nothing with the mark -------------------------------------------------------


@committed(RUNTIME_FIXTURE)
def test_the_mark_can_be_carried_without_editing_the_ir_and_is_still_a_rule() -> None:
    rule = _cell_gated_rule(RUNTIME_FIXTURE)
    marked = cp.mark(rule, cp.ARGMAX_OF_PREDICTED_EFFECT)
    assert isinstance(marked, Rule)
    assert cp.provenance_mark(marked) == cp.ARGMAX_OF_PREDICTED_EFFECT
    assert cp.provenance_mark(rule) is cp.NOT_ASSESSED  # an unmarked rule reads as the default
    for field in ("id", "source", "action", "target", "strength", "threshold", "hill", "when"):
        assert getattr(marked, field) == getattr(rule, field)


def test_mark_refuses_a_value_off_the_axis() -> None:
    with pytest.raises(ValueError, match="not on the axis"):
        cp.mark(_bare_rule(), "measured")


@pytest.mark.parametrize("value", cp.AXIS_VALUES)
@committed(RUNTIME_FIXTURE)
def test_the_active_rule_set_does_not_move_when_every_rule_is_marked(value: str) -> None:
    for context in ({"cell_type": "Cardiomyocyte2"}, {"cell_type": "Neuron2"}, {}):
        plain = _parsed(RUNTIME_FIXTURE)
        before = [(r.id, r.source, r.target) for r in plain.active_rules(context)]
        marked_module = parse_file(ROOT / RUNTIME_FIXTURE)
        marked_module.rules = [cp.mark(r, value) for r in marked_module.rules]
        after = [(r.id, r.source, r.target) for r in marked_module.active_rules(context)]
        assert before == after


@pytest.mark.parametrize("value", cp.AXIS_VALUES)
@committed(RUNTIME_FIXTURE)
def test_the_trajectory_does_not_move_when_every_rule_is_marked(value: str) -> None:
    context = {"cell_type": "Cardiomyocyte2"}
    base = NetworkRuntime(parse_file(ROOT / RUNTIME_FIXTURE), context=context, seed=0).run(
        hours=8, dt=0.05, record_every=4
    )
    m = parse_file(ROOT / RUNTIME_FIXTURE)
    m.rules = [cp.mark(r, value) for r in m.rules]
    got = NetworkRuntime(m, context=context, seed=0).run(hours=8, dt=0.05, record_every=4)
    assert got.times == base.times
    assert got.species == base.species
    assert got.levels == base.levels


@committed(RUNTIME_FIXTURE)
def test_applies_reads_the_when_clause_and_not_the_mark() -> None:
    rule = _cell_gated_rule(RUNTIME_FIXTURE)
    wanted = rule.when["cell_type"]
    for value in cp.AXIS_VALUES:
        marked = cp.mark(rule, value)
        assert marked.applies({"cell_type": wanted}) is True
        assert marked.applies({"cell_type": "Neuron2"}) is False


# --- P1: the compiler has exactly two sites that put a cell on a rule --------------------------------


@committed(COMPILER)
def test_the_compiler_emits_a_rule_line_from_exactly_two_sites() -> None:
    tree = ast.parse(must_be_committed(COMPILER).read_text())
    sites = []
    for func in ast.walk(tree):
        if not isinstance(func, ast.FunctionDef):
            continue
        for node in ast.walk(func):
            if not isinstance(node, ast.JoinedStr) or not node.values:
                continue
            first = node.values[0]
            if isinstance(first, ast.Constant) and str(first.value).startswith("rule "):
                sites.append(func.name)
    assert sorted(sites) == ["_measured_blocks", "compile_chromosome"], sites


def test_the_two_signatures_are_read_from_the_code_and_not_copied_as_literals() -> None:
    assert cp._crispri_source_prefix() == SOURCES["crispri"]
    source = Path(cp.__file__).read_text()
    assert SOURCES["crispri"] not in source


# --- P2 and P6: the committed compiled program --------------------------------------------------------


@cache
def _compiled_classes() -> dict[str, int]:
    rules = _parsed(COMPILED).rules
    counts: dict[str, int] = {}
    for r in rules:
        counts[cp.provenance_of(r)] = counts.get(cp.provenance_of(r), 0) + 1
    counts["_rules"] = len(rules)
    return counts


@committed(COMPILED)
def test_every_rule_of_the_committed_compiled_program_classifies_from_the_rule_alone() -> None:
    """P2. No rule needs a human, and the two classes partition the program."""
    counts = _compiled_classes()
    assert counts["_rules"] > 0
    assert set(counts) - {"_rules"} == {cp.ARGMAX_OF_PREDICTED_EFFECT, cp.MEASURED_PERTURBATION_IN_THAT_CELL}
    assert (
        counts[cp.ARGMAX_OF_PREDICTED_EFFECT] + counts[cp.MEASURED_PERTURBATION_IN_THAT_CELL]
        == counts["_rules"]
    )


@committed(COMPILED)
def test_the_argmax_class_is_almost_the_whole_compiled_program() -> None:
    counts = _compiled_classes()
    assert counts[cp.ARGMAX_OF_PREDICTED_EFFECT] / counts["_rules"] > 0.99


@committed(COMPILED)
def test_a_gate_at_the_measured_class_leaves_under_one_per_cent_integrable() -> None:
    """P6. A gate is not a stricter corpus but an unrunnable one; counted on this axis's population.

    The sibling (`6150aa9`) counted this on its own axis and its own program; nothing of its figure
    is reused here.
    """
    rules = _parsed(COMPILED).rules
    marked = [cp.mark(r, cp.provenance_of(r)) for r in rules]
    survivors = cp.gate_survivors(marked, cp.MEASURED_PERTURBATION_IN_THAT_CELL)
    assert 0 < len(survivors) / len(rules) < 0.01


def test_a_gate_can_only_be_set_at_a_class_and_never_at_the_default() -> None:
    with pytest.raises(ValueError, match="only be set at a class"):
        cp.gate_survivors([], cp.NOT_ASSESSED)


# --- P4: the code-reachable class that does not occur -------------------------------------------------


@committed(COMPILED)
def test_no_rule_of_the_committed_compiled_program_is_gated_on_an_unrecorded_cell() -> None:
    """P4. `compile.context` can write it; the committed program has none of it."""
    rules = _parsed(COMPILED).rules
    assert [r.id for r in rules if r.when.get("cell_type") == CONTEXT_UNKNOWN] == []
    assert all(cp.names_a_cell(r) for r in rules)


def test_a_rule_with_no_recorded_cell_is_outside_the_population_and_already_runs_nowhere() -> None:
    rule = _bare_rule()
    rule.when = {"cell_type": CONTEXT_UNKNOWN}
    assert not cp.names_a_cell(rule)
    assert rule.applies({"cell_type": "K562"}) is False  # the runtime's existing refusal, untouched
    with pytest.raises(cp.CellProvenanceUndecidableError, match="outside the population"):
        cp.provenance_of(rule)


# --- P3: the hand-authored corpus, and the one adjudication ------------------------------------------


@committed(*HAND_AUTHORED)
def test_the_hand_authored_corpus_carries_cell_gated_rules_that_no_signature_classifies() -> None:
    """P3. These are the marks a human would have to decide, and the count is the adoption cost."""
    undecidable = []
    for relative in HAND_AUTHORED:
        for rule in _parsed(relative).rules:
            if not cp.names_a_cell(rule):
                continue
            with pytest.raises(cp.CellProvenanceUndecidableError, match="must adjudicate"):
                cp.provenance_of(rule)
            undecidable.append(rule.id)
    assert undecidable == sorted(ADJUDICATED)
    assert set(ADJUDICATED.values()) == {cp.AUTHOR_DECLARED}


@committed(RUNTIME_FIXTURE)
def test_the_adjudicated_rule_carries_a_citation_that_names_no_cell() -> None:
    """Why the adjudication is `author_declared` and not the measured class: the program cites a
    source for the FACTOR and says nothing that ties it to the cell on the rule line."""
    rule = next(r for r in _parsed(RUNTIME_FIXTURE).rules if r.id in ADJUDICATED)
    assert rule.when["cell_type"] not in (rule.evidence.source or "")
    assert not (rule.evidence.source or "").startswith(SOURCES["crispri"])


# --- P8: the parser synthesises rules, and not one of them names a cell -------------------------------


@committed(COMPILED, *HAND_AUTHORED)
def test_every_rule_the_parser_synthesises_names_no_cell() -> None:
    """P8. The synthesis is real and large, and it leaves this proposal's population untouched -
    the opposite of what it did to the sibling's (`6150aa9`, cited, not restated)."""
    synthesised = 0
    for relative in (COMPILED, *HAND_AUTHORED):
        written = _written_rule_headers(relative)
        for rule in _parsed(relative).rules:
            if rule.id in written:
                continue
            synthesised += 1
            assert not cp.names_a_cell(rule), (relative, rule.id, rule.when)
    assert synthesised > 0


@committed(COMPILED)
def test_the_compiled_program_synthesises_nothing_so_its_rules_are_all_written() -> None:
    written = _written_rule_headers(COMPILED)
    rules = _parsed(COMPILED).rules
    assert len(rules) == len(written)
    assert {r.id for r in rules} == written


# --- what the committed program does and does not say about its own cells -----------------------------


@committed(COMPILED)
def test_the_compiled_program_nowhere_states_that_a_cell_is_a_selection() -> None:
    """The gap the proposal addresses, measured on the artefact rather than argued."""
    text = must_be_committed(COMPILED).read_text()
    for phrase in ("argmax", "moves most", "strongest track", "largest predicted"):
        assert phrase not in text


@committed(COMPILED)
def test_the_compiled_header_already_asserts_the_reading_the_proposal_says_is_false() -> None:
    """An adoption cost nobody asked for: the compiler's OWN header says the cell is one the rule
    was "measured or predicted in". On the argmax rules that is the misreading itself, written into
    the artefact by the tool, so adopting the proposal means changing this line too."""
    text = must_be_committed(COMPILED).read_text()
    assert "gated on the cell it was measured or predicted in" in text
    compiler = must_be_committed(COMPILER).read_text()
    assert "gated on the cell it was measured or predicted in" in compiler


@committed(COMPILED)
def test_the_class_is_separable_today_so_the_mark_adds_words_and_not_a_distinction() -> None:
    """Reported as a cost to the proposal, not hidden: a reader with the code in hand can already
    tell the two classes apart from `evidence:`. What the file does not carry is the SENTENCE."""
    rules = _parsed(COMPILED).rules
    kinds = {r.evidence.kind.value for r in rules}
    assert kinds == {"predicted", "experimental"}
    by_kind = {k: set() for k in kinds}
    for r in rules:
        by_kind[r.evidence.kind.value].add(cp.provenance_of(r))
    assert by_kind == {
        "predicted": {cp.ARGMAX_OF_PREDICTED_EFFECT},
        "experimental": {cp.MEASURED_PERTURBATION_IN_THAT_CELL},
    }


@committed(COMPILED)
def test_the_precedent_field_is_present_in_the_grammar_and_absent_on_every_committed_rule() -> None:
    """`context_evidence` is the shape this proposal copies, and the only compiled program git
    tracks carries it on NO rule: the mapping table was unavailable when it was generated, so no
    reading was taken. A field that is empty everywhere it occurs is the precedent's own warning."""
    assert all(r.context_evidence == "" for r in _parsed(COMPILED).rules)
    assert "context_evidence" not in must_be_committed(COMPILED).read_text()


# --- the registration itself --------------------------------------------------------------------------


def test_the_registration_recommends_nothing() -> None:
    reg = cp.registration()
    assert "no recommendation" in reg["no_recommendation"].lower()
    text = Path(cp.__file__).read_text().lower()
    for word in ("we recommend", "should be adopted", "i recommend"):
        assert word not in text


def test_the_registration_left_the_cost_and_the_outcomes_empty() -> None:
    """Pinned so the commit that fills them is visibly later than the one that registered them."""
    reg = cp.registration()
    assert reg["adoption_cost"] == {}
    assert reg["prediction_outcomes"] == {}


def test_eight_predictions_each_with_what_would_falsify_it() -> None:
    ids = [p["id"] for p in cp.PREDICTIONS]
    assert ids == [f"P{i}" for i in range(1, 9)]
    assert all(p["falsified_by"].strip() for p in cp.PREDICTIONS)


def test_the_carried_readings_are_present_and_unstrengthened() -> None:
    carried = cp.CARRIED_READINGS
    assert "can never reverse one" in carried["withholds_never_reverses"]
    assert "NOT absence of regulation" in carried["unresolved_is_not_absence"]
    assert '"not known to be biological replicates"' in carried["two_track_cell"]
    assert "ONE BIOSAMPLE NAME is the most" in carried["one_biosample_name"]


def test_a_committed_path_that_git_does_not_track_fails_rather_than_skipping() -> None:
    """The guard shape `4f44dbf` requires: absence of a committed artefact is a defect, not a skip."""
    with pytest.raises(AssertionError, match="git tracks it: False"):
        must_be_committed("data/organisms/human/noncoding_chr99.bio")


def test_no_test_in_this_file_skips_on_a_committed_path() -> None:
    tree = ast.parse(Path(__file__).read_text())
    guards = []
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            guards += [ast.unparse(d) for d in node.decorator_list if "skip" in ast.unparse(d)]
        if isinstance(node, ast.Call) and "skip" in ast.unparse(node.func):
            guards.append(ast.unparse(node.func))
    assert guards == [], guards
    assert re.search(r"must_be_committed|@committed", Path(__file__).read_text())


# --- the adoption cost, counted rather than estimated -------------------------------------------------


@committed(COMPILED, *HAND_AUTHORED)
def test_the_counted_adoption_cost_is_one_human_decision_and_the_rest_mechanical() -> None:
    """The proposal's own scale. A per-rule mark over every tracked program, split by who pays."""
    cost = _cost().count()
    assert cost["marks_a_human_must_decide"] == 1
    assert cost["marks_a_human_must_decide_rules"] == [
        {"program": RUNTIME_FIXTURE, "rule": next(iter(ADJUDICATED))}
    ]
    assert cost["marks_a_machine_can_derive"] == cost["rules_in_the_population"] - 1
    assert cost["marks_a_machine_can_derive_by_class"] == {
        cp.ARGMAX_OF_PREDICTED_EFFECT: len(_parsed(COMPILED).rules) - 2,
        cp.MEASURED_PERTURBATION_IN_THAT_CELL: 2,
    }
    assert cost["synthesised_rules_naming_a_cell"] == 0
    assert cost["default_for_every_one_of_them_before_a_decision"] is cp.NOT_ASSESSED


@committed(COMPILED, *HAND_AUTHORED)
def test_the_cost_scan_raises_rather_than_skipping_when_a_tracked_program_is_absent(monkeypatch) -> None:
    """A guard that skips on an absent COMMITTED artefact is the defect `4f44dbf` was written
    against, so the scan's own absence path raises. Shown by naming a tracked path that is not there."""
    cost = _cost()
    monkeypatch.setattr(cost, "_tracked_bio", lambda: ["data/organisms/human/noncoding_chr99.bio"])
    with pytest.raises(AssertionError, match="absence is a defect"):
        cost.count()


def test_every_prediction_has_an_outcome_whose_test_exists_in_this_file() -> None:
    cost = _cost()
    assert sorted(cost.OUTCOMES) == [p["id"] for p in cp.PREDICTIONS]
    here = Path(__file__).read_text()
    for pid, outcome in cost.OUTCOMES.items():
        assert outcome["verdict"].startswith("held") or outcome["verdict"].startswith("false"), pid
        assert outcome["figure"].strip(), pid
        assert f"def {outcome['test']}(" in here, (pid, outcome["test"])


def test_what_was_found_without_a_prediction_is_reported_and_not_empty() -> None:
    assert len(_cost().FOUND_WITHOUT_A_PREDICTION) >= 3
    assert all(line.strip() for line in _cost().FOUND_WITHOUT_A_PREDICTION)


def test_the_cost_report_makes_no_recommendation() -> None:
    source = Path(_cost().__file__).read_text().lower()
    assert "no recommendation" in source
    for word in ("we recommend", "should be adopted", "i recommend"):
        assert word not in source
