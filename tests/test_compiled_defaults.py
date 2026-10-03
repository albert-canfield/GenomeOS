# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The defaults a compiled program leans on, and the per-slot evidence tier, as the program says them.

Adopted 2026-10-03 (Albert, item 7). Three claims are tested, and every guard behind them is PLANTED
against before it is believed: a guard that cannot be shown to REFUSE is a list and not a check.

1. the per-slot mark is the registered proposal's own vocabulary, keyed by field name and never by
   rule, and `not_assessed` is the absent value - never `U_unsourced`, which asserts a search;
2. the numeric defaults the engine applies where a rule line is silent are PRINTED in the program,
   read off `genomeos/ir/model.py`'s own `Rule` dataclass rather than transcribed;
3. an absent `context_evidence` comes out as `not_assessed`, which is not `not_assessable`.

NOTHING HERE WRITES A `.bio` FILE. The published programs are read only, and two tests pin their
bytes: v1 against the sha256 `data/results/label_gene.json` names, v2 against its committed blob.
"""

from __future__ import annotations

import dataclasses
import hashlib
import subprocess
from pathlib import Path

import pytest

from genomeos.attribution import compile as cp
from genomeos.ir.model import Rule
from genomeos.provenance import compiled_defaults as cd
from genomeos.provenance import rule_evidence_tier as ret

ROOT = Path(__file__).resolve().parents[1]
COMPILER = ROOT / "genomeos" / "attribution" / "compile.py"
V1 = ROOT / "data" / "organisms" / "human" / "noncoding_chr21.bio"
V2 = ROOT / "data" / "organisms" / "human" / "noncoding_chr21_v2.bio"

#: the sha256 `data/results/label_gene.json` and `label_gene_registration.json` name for v1, so v1's
#: bytes are a pin: a cut that moved them would break a rebuild of those results
V1_PINNED_SHA256 = "778615512beaa46ddbb02f3cbee8502ceb1ff8d07dd43ce553e086b29e035192"

#: a real cut needs chr21's attribution results on this machine
HAS_CHR21 = (ROOT / "data" / "results" / "budget_chr21.json").exists()
needs_chr21 = pytest.mark.skipif(not HAS_CHR21, reason="chr21's budget is not on this machine")


@dataclasses.dataclass
class _DriftedRule:
    """A `Rule` whose threshold default has moved. Injected, so the live `Rule` every parse uses is
    not touched: the question is whether the guard notices a drift, not whether a drift can be made."""

    strength: float = 1.0
    threshold: float = 1.5
    hill: float = 2.0
    context_evidence: str = ""


@dataclasses.dataclass
class _RuleWithoutHill:
    """A `Rule` the engine has stopped carrying a `hill` on at all."""

    strength: float = 1.0
    threshold: float = 1.0
    context_evidence: str = ""


class _MarkedRule:
    """A rule carrying per-slot marks, which today's `Rule` does not. Duck-typed on purpose."""

    def __init__(self, marks: dict[str, str]) -> None:
        self.evidence_tier = marks


def _block() -> str:
    return "\n".join(cd.declaration_lines())


# --- 1. the vocabulary is the registered proposal's, as the same objects ---------------------------


def test_the_absent_value_and_the_slots_are_the_proposals_own_objects_not_a_second_copy():
    """A second spelling of `not_assessed` in the language is how the two levels start disagreeing."""
    assert cd.ABSENT is ret.NOT_ASSESSED
    assert cd.TIERED_FIELDS is ret.TIERED_FIELDS
    assert cd.TIERS is ret.TIERS
    assert cd.AXIS_VALUES is ret.AXIS_VALUES


def test_the_absent_value_is_not_a_tier_and_cannot_be_ordered_against_one():
    assert cd.ABSENT not in cd.TIERS
    assert not ret.is_tier(cd.ABSENT)
    with pytest.raises(ValueError, match="no rank"):
        ret.rank(cd.ABSENT)


def test_the_absent_value_is_not_spelled_unsourced_anywhere_it_is_written():
    """`U_unsourced` is a finding - a search was made and came back empty. The compiler made none."""
    assert cd.FORBIDDEN_ABSENT_SPELLING not in cd.ABSENT
    for field in cd.TIERED_FIELDS:
        for tier in cd.TIERS:
            assert f"{field} = {tier}" not in _block()


# --- 2. the mark is per slot and not per rule -------------------------------------------------------


def test_a_rule_with_no_marks_reads_not_assessed_on_every_slot_and_on_no_other_field():
    marks = cd.slot_tiers(Rule(id="r", source="A", action=None, target="B"))
    assert marks == dict.fromkeys(cd.TIERED_FIELDS, cd.ABSENT)
    assert len(marks) == len(ret.TIERED_FIELDS) == 3


def test_three_slots_of_one_rule_can_hold_three_different_marks():
    """The census's own finding (`4d4003d`): tiers DIFFER within one rule, which a per-rule mark
    would discard. This is the whole reason the mark is keyed by field name."""
    marks = cd.slot_tiers(_MarkedRule({"strength": "E_existence_only", "hill": "F_fitted_or_modelled"}))
    assert marks == {
        "strength": "E_existence_only",
        "threshold": cd.ABSENT,
        "hill": "F_fitted_or_modelled",
    }
    assert len(set(marks.values())) == 3


def test_PLANT_a_mark_off_the_axis_is_refused_and_never_counted_as_unassessed():
    with pytest.raises(ValueError, match="not on the axis"):
        cd.slot_tiers(_MarkedRule({"strength": "probably_fine"}))


def test_PLANT_a_mark_on_a_field_that_carries_no_number_is_refused():
    with pytest.raises(ValueError, match="carries no number"):
        cd.slot_tiers(_MarkedRule({"strength": cd.ABSENT, "cell_type": "M_measured_quoted"}))


def test_PLANT_evidence_tier_that_is_not_a_mapping_is_refused():
    rule = _MarkedRule({})
    rule.evidence_tier = ["M_measured_quoted"]
    with pytest.raises(ValueError, match="must be a mapping"):
        cd.slot_tiers(rule)


# --- 3. the defaults are read off the dataclass, not transcribed -----------------------------------


def test_the_declared_defaults_are_the_live_dataclass_defaults():
    """Recomputed here from `dataclasses.fields` independently of the module, so the two readings
    have to agree rather than share one source of error."""
    live = {f.name: f.default for f in dataclasses.fields(Rule)}
    assert cd.rule_defaults() == {n: live[n] for n in cd.DECLARED}
    assert cd.rule_defaults()["threshold"] == 1.0
    assert cd.rule_defaults()["hill"] == 2.0


def test_the_block_prints_every_numeric_default_with_its_live_value():
    block = _block()
    for name in cd.DECLARED_NUMBERS:
        assert f"`{name}` is {cd._number(cd.rule_defaults()[name])}" in block


def test_the_block_says_what_an_absent_context_evidence_reads():
    assert f"`context_evidence:` reads `{cd.ABSENT}`" in _block()
    assert "not `not_assessable`" in _block()


def test_PLANT_rule_defaults_refuses_a_dataclass_that_has_dropped_a_declared_field():
    with pytest.raises(ValueError, match=r"has no field \['hill'\]"):
        cd.rule_defaults(_RuleWithoutHill)


def test_PLANT_rule_defaults_refuses_something_that_is_not_a_dataclass():
    with pytest.raises(ValueError, match="not a dataclass"):
        cd.rule_defaults(int)


# --- 4. the declaration guard REFUSES, shown four ways ---------------------------------------------


def test_the_guard_accepts_the_block_the_writer_writes():
    cd.check_declaration(_block())


def test_PLANT_the_guard_refuses_a_program_with_no_declaration_at_all():
    with pytest.raises(ValueError, match="declares none of the defaults"):
        cd.check_declaration("module m\n\n# nothing about defaults here\n")


def test_PLANT_the_guard_refuses_a_declaration_that_has_drifted_from_the_dataclass():
    """Injection rather than mutation: the engine's `Rule` is left alone and a drifted one is
    handed in, so this shows the guard comparing against the dataclass and not against itself."""
    with pytest.raises(ValueError, match="threshold=1.5"):
        cd.check_declaration(_block(), rule_type=_DriftedRule)


def test_PLANT_the_guard_refuses_a_block_that_marks_no_tier_for_a_number_the_census_tiers():
    """What happens when a fourth tiered number is added upstream: the block must grow a mark for
    it, and until it does the guard refuses rather than passing three marks for four numbers."""
    with pytest.raises(ValueError, match="no tier for 'decay'"):
        cd.check_declaration(_block(), tiered_fields=(*cd.TIERED_FIELDS, "decay"))


@pytest.mark.parametrize("tier", list(ret.TIERS))
def test_PLANT_the_guard_refuses_a_slot_marked_with_a_tier_the_compiler_never_adjudicated(tier):
    planted = _block().replace(f"hill = {cd.ABSENT}", f"hill = {tier}")
    assert f"hill = {tier}" in planted
    with pytest.raises(ValueError, match="a tier: a compiled program asserts no adjudication"):
        cd.check_declaration(planted)


def test_PLANT_the_guard_refuses_a_block_that_drops_the_absent_context_evidence_sentence():
    planted = _block().replace(f"`context_evidence:` reads `{cd.ABSENT}`", "`context_evidence:` is blank")
    with pytest.raises(ValueError, match="does not say what an absent"):
        cd.check_declaration(planted)


# --- 5. not_assessed carries no reason, and that guard REFUSES too ---------------------------------


def test_an_absent_context_evidence_reads_not_assessed_and_a_present_one_passes_through():
    assert cd.context_evidence_value("") == cd.ABSENT
    assert cd.context_evidence_value("   ") == cd.ABSENT
    assert cd.context_evidence_value(None) == cd.ABSENT  # type: ignore[arg-type]
    assert cd.context_evidence_value("open_in_reader, K562") == "open_in_reader, K562"


def test_a_reading_that_was_attempted_and_not_takeable_keeps_its_reason():
    """`not_assessable` is the one that carries a reason, and this guard must not touch it."""
    cd.check_not_assessed_carries_no_reason(
        "rule e1 ACTIVATE g1 { context_evidence: not_assessable, no_reader_for_this_cell, K562; }"
    )


def test_PLANT_not_assessed_with_a_reason_is_refused_because_a_reason_means_an_attempt():
    with pytest.raises(ValueError, match="carries a reason"):
        cd.check_not_assessed_carries_no_reason(
            f"rule e1 ACTIVATE g1 {{ context_evidence: {cd.ABSENT}, no_reader_for_this_cell; }}"
        )


# --- 6. the compiler writes it, and REFUSES when it is wrong: reachability, not a branch -----------


def test_the_omitted_field_reads_not_assessed_even_though_the_line_stays_silent():
    """WITHDRAWN AND SAID SO: the rule LINE keeps the behaviour another lane's test pins - no reader
    means no field - and what that omission reads is `not_assessed`, which the program's header now
    states. The reading is the adoption; writing it onto the line is that lane's decision.

    The pin: `tests/test_context_evidence.py::test_the_state_is_added_beside_the_rule_and_deletes_
    nothing`, "the only difference the reading makes to a program is the field it adds". The second
    reason: `attribution.context_evidence.parse_value` refuses a value that does not begin with one
    of its three states, and `not_assessed` is not one of them.
    """
    assert cp.context_evidence_field(None, "K562", 1, 2) == ""
    assert cd.context_evidence_value("") == cd.ABSENT
    assert f"`context_evidence:` reads `{cd.ABSENT}`" in _block()


def test_the_field_passes_a_real_reading_through_unchanged():
    def state(label, start, end):
        return f"open_in_reader, {label}"

    assert cp.context_evidence_field(state, "K562", 1, 2) == "context_evidence: open_in_reader, K562; "


def test_every_place_the_compiler_writes_the_field_goes_through_the_one_helper():
    """Two emission sites wrote the field inline, so a change to what an absent reading means had two
    places to be made and one to be forgotten. There is now one. Planted below."""
    src = COMPILER.read_text()
    assert src.count("context_evidence: {context_state(") == 1
    assert src.count("def context_evidence_field(") == 1
    assert src.count("ctx = context_evidence_field(") == 2


def test_PLANT_a_second_inline_emission_site_in_a_copy_of_the_compiler_is_caught(tmp_path):
    copy = tmp_path / "compile.py"
    copy.write_text(
        COMPILER.read_text().replace(
            '            ctx = context_evidence_field(context_state, label, e["start"], e["end"])',
            '            ctx = ""\n'
            "            if context_state is not None:\n"
            "                ctx = f\"context_evidence: {context_state(label, e['start'], e['end'])}; \"",
        )
    )
    src = copy.read_text()
    assert src.count("context_evidence: {context_state(") == 2
    assert src.count("ctx = context_evidence_field(") == 1


@needs_chr21
def test_a_real_cut_carries_the_declaration_and_the_guards_pass_on_it():
    text = cp.compile_chromosome("chr21")
    assert cd.BLOCK_MARKER in text
    cd.check_declaration(text)
    cd.check_not_assessed_carries_no_reason(text)


@needs_chr21
def test_the_declaration_is_in_the_header_and_on_no_rule_line_so_a_cut_moves_no_rule(monkeypatch):
    """A mark on every rule LINE would move all 5,176 of chr21's on the next cut, which is a re-cut
    and Albert's decision. The declaration is program-level, so the rule lines are untouched by it."""
    text = cp.compile_chromosome("chr21")
    rules = [ln for ln in text.splitlines() if ln.startswith("rule ")]
    assert len(rules) > 5000
    for line in rules:
        assert f"= {cd.ABSENT}" not in line
        assert cd.BLOCK_MARKER not in line
        assert "context_evidence: " in line


@needs_chr21
def test_PLANT_the_compiler_itself_refuses_a_cut_whose_declaration_drifted(monkeypatch):
    """Reachability, not a branch: the guard is reached through `compile_chromosome`, the only entry
    point that writes a program, with a drifted block injected where the real one is written."""
    drifted = [ln.replace("`threshold` is 1.0", "`threshold` is 9.9") for ln in cd.declaration_lines()]
    assert drifted != cd.declaration_lines()
    monkeypatch.setattr(cd, "declaration_lines", lambda *a, **k: drifted)
    with pytest.raises(ValueError, match="the declaration does not say"):
        cp.compile_chromosome("chr21")


@needs_chr21
def test_PLANT_the_compiler_itself_refuses_a_cut_that_hands_not_assessed_a_reason(monkeypatch):
    """Reachability again: injected at the one place the field is written, the refusal comes out of
    `compile_chromosome` and not out of a direct call to the guard."""
    monkeypatch.setattr(
        cp,
        "context_evidence_field",
        lambda *a, **k: f"context_evidence: {cd.ABSENT}, no_reader_for_this_cell; ",
    )
    with pytest.raises(ValueError, match="carries a reason"):
        cp.compile_chromosome("chr21")


# --- 7. the published programs' bytes did not move -------------------------------------------------


def test_v1s_bytes_are_still_the_sha256_the_committed_results_pin():
    assert hashlib.sha256(V1.read_bytes()).hexdigest() == V1_PINNED_SHA256


def test_v2s_bytes_are_still_its_committed_blob():
    """v2 is another lane's artefact and this lane does not re-cut it: an adopted default that
    changes a published program's bytes is a re-cut, and a re-cut is Albert's decision."""
    blob = subprocess.run(
        ["git", "show", "HEAD:data/organisms/human/noncoding_chr21_v2.bio"],
        cwd=ROOT,
        capture_output=True,
    )
    assert blob.returncode == 0, blob.stderr.decode()[:200]
    assert V2.read_bytes() == blob.stdout


def _bio_digests() -> dict[str, str]:
    return {
        p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted((ROOT / "data").rglob("*.bio"))
    }


@needs_chr21
def test_a_cut_and_every_guard_in_this_suite_leave_every_bio_file_byte_identical():
    """The adoption changes what a cut WOULD write; it writes nothing. Digested before and after a
    real cut and both guards, over the whole corpus and not only the two programs read above."""
    before = _bio_digests()
    assert before, "the corpus under test must not be empty or this test proves nothing"
    text = cp.compile_chromosome("chr21")
    cd.check_declaration(text)
    cd.check_not_assessed_carries_no_reason(text)
    assert _bio_digests() == before
