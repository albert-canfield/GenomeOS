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

#: `needs_local_data`, the project's marker for a test that reads a git-ignored machine-local store
#: (`tests/local_data.py`, the supervisor's AMENDED acceptance rule of 2026-10-02): such a test RUNS in a
#: checkout or a verdict worktree that has the store, and where the store is genuinely absent -- CI, or a
#: bare worktree -- it SKIPS BY NAME, so a reader can tell "not run here" from "passed".
#:
#: Why the five tests below need it, and why `needs_chr21` above was not enough. `compile_chromosome`
#: reaches `targets.run_elements`, and the committed summary `enhancer_targets_all_chr21` carries no
#: inline elements: it points at `data/knowledge/alphagenome/all_elements/chr21.json`, and
#: `data/knowledge` is git-ignored (.gitignore:17), machine-local by the data boundary. `run_elements`
#: then raises FileNotFoundError -- "points at ..., which is not on this machine" -- which is the module
#: behaving correctly: a summary whose table is missing must not read as a run with no elements.
#: `needs_chr21` guards on `data/results/budget_chr21.json`, which git TRACKS, so it is present in every
#: checkout and that skipif never fires; it is kept because it names the other input a real cut needs,
#: and the store the cut actually lacks in CI is named here. CI run 37124753522 on 5c4589d failed all
#: five, in the `test` job, while every one of them was green in this checkout -- so a local green proved
#: nothing about them. tests/test_header_cell_sentence.py carries the same marker for the same defect in
#: the same compiler path (`b04f48a`), and this names the same store in the same form.
#:
#: The marker moves WHERE these tests run, never what they claim. The two reachability PLANTS still
#: reach the guard through `compile_chromosome`, the only entry point that writes a program, rather than
#: supplying their own trigger -- that is the property they exist for. The table is NOT stubbed and NOT
#: mocked, deliberately: the counts below (> 5000 rule lines, and the whole-corpus digest) are readings
#: of the real table, so a planted table would make the plants assert against a fiction.
needs_the_chr21_element_table = pytest.mark.needs_local_data(
    "data/knowledge/alphagenome/all_elements/chr21.json",
    how="scripts/enhancer_targets_all_chain.py writes it (the enhancer_targets_all_chr21 job); "
    "data/knowledge is git-ignored machine-local data, so a fresh checkout cannot have it",
)


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


# SUPERSEDED 2026-10-03, hours after it was written, and KEPT rather than edited: the clause it
# records as another lane's decision has now been decided. Albert's item (7) names `not_assessed` in
# his own words, so `context_evidence_field(None, ...)` no longer returns `""` and the first assertion
# below is FALSE. `strict=True` means it must keep failing. Both named blockers were moved without
# being weakened: `parse_value` was TAUGHT the fourth value (the three readings byte-identical), and
# `test_the_state_is_added_beside_the_rule_and_deletes_nothing` was superseded additively in its own
# file with `xfail(strict=True)` and a corrected replacement. The corrected assertion here is
# `test_the_line_now_states_not_assessed_instead_of_staying_silent`, immediately below, which keeps
# the two assertions of this test that are still true.
@pytest.mark.xfail(
    strict=True,
    reason=(
        "superseded 2026-10-03 (Albert, item 7): the withdrawn piece this test records was adopted, "
        "so context_evidence_field(None, ...) now returns 'context_evidence: not_assessed; ' and not "
        "''. Replaced by test_the_line_now_states_not_assessed_instead_of_staying_silent. Kept so the "
        "record shows what was withdrawn, why, and that neither blocker was weakened to move it."
    ),
)
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


def test_the_line_now_states_not_assessed_instead_of_staying_silent():
    """The adopted form of the test xfailed above: the line SAYS it, it is not inferred from silence.

    Two of that test's three assertions are still true and are kept here. The READING of an empty
    stored value does not move, and it must not: v1 and v2 were cut before today and carry no field at
    all, so `context_evidence_value("")` is what makes them readable, and the header still has to say
    so for a program that stands on its own.
    """
    assert cp.context_evidence_field(None, "K562", 1, 2) == f"context_evidence: {cd.ABSENT}; "
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


@needs_the_chr21_element_table
@needs_chr21
def test_a_real_cut_carries_the_declaration_and_the_guards_pass_on_it():
    text = cp.compile_chromosome("chr21")
    assert cd.BLOCK_MARKER in text
    cd.check_declaration(text)
    cd.check_not_assessed_carries_no_reason(text)


@needs_the_chr21_element_table
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


@needs_the_chr21_element_table
@needs_chr21
def test_PLANT_the_compiler_itself_refuses_a_cut_whose_declaration_drifted(monkeypatch):
    """Reachability, not a branch: the guard is reached through `compile_chromosome`, the only entry
    point that writes a program, with a drifted block injected where the real one is written."""
    drifted = [ln.replace("`threshold` is 1.0", "`threshold` is 9.9") for ln in cd.declaration_lines()]
    assert drifted != cd.declaration_lines()
    monkeypatch.setattr(cd, "declaration_lines", lambda *a, **k: drifted)
    with pytest.raises(ValueError, match="the declaration does not say"):
        cp.compile_chromosome("chr21")


@needs_the_chr21_element_table
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


@needs_chr21
def test_a_cut_that_READ_NOTHING_states_not_assessed_on_every_rule_line(monkeypatch):
    """Reachability, not a branch: the None path is reached through `compile_chromosome` itself, with
    the mapping table taken away where the compiler asks for it, which is what a machine without that
    git-ignored store is. 5,176 rule lines then carry the fourth value, bare, and both guards pass."""
    from genomeos.attribution import context_evidence as ce

    def no_table(*a, **k):
        raise FileNotFoundError("the cell-to-biosample mapping table is not on this machine")

    monkeypatch.setattr(ce, "mapping", no_table)
    text = cp.compile_chromosome("chr21")
    rules = [ln for ln in text.splitlines() if ln.startswith("rule ")]
    assert len(rules) > 5000
    for line in rules:
        assert f"context_evidence: {cd.ABSENT}; " in line, line[:120]
        written = line.split("context_evidence: ", 1)[1].split(";", 1)[0]
        assert ce.parse_value(written) == (cd.ABSENT, ()), line[:120]
    # a bare `not_assessed` carries no reason, so the guard that refuses one passes on a real cut
    cd.check_declaration(text)
    cd.check_not_assessed_carries_no_reason(text)


@needs_chr21
def test_the_superseded_section_opening_is_kept_in_the_source_and_no_longer_written(monkeypatch):
    """Zero deletions, shown rather than claimed: the two sentences the no-reader section used to open
    with are preserved verbatim in the compiler and are absent from what it now writes.

    One of them said "no rule states the field", which is false of every cut made since Albert's item
    (7). Rewriting it in place would have removed a committed line; keeping it where a reader can
    compare it with its replacement is what `dbb5d4a` established as the route.
    """
    from genomeos.attribution import context_evidence as ce

    kept = cp.SUPERSEDED_ABSENT_SECTION_OPENING
    assert kept == (
        "# ---- context evidence: NONE was read for this cut, so no rule states the field and",
        "# every rule's reading is `not_assessed`. The mapping table from cell label to ontology",
    )
    monkeypatch.setattr(ce, "mapping", lambda *a, **k: (_ for _ in ()).throw(FileNotFoundError("none")))
    text = cp.compile_chromosome("chr21")
    for line in kept:
        assert line.lstrip("# -").strip() not in text, line
    # and what IS written says the opposite, in the program and not only in this test
    assert "NONE was read for this cut, so every rule line STATES the" in text


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


#: The re-cut made under this lane, under a NEW NAME and beside v2 rather than instead of it. Whether
#: the PUBLISHED program moves to v3 is the re-cut question still on Albert's list and is not decided
#: here: v1 stays what `data/results/label_gene.json` pins and v2 stays what the README points at.
V3 = ROOT / "data" / "organisms" / "human" / "noncoding_chr21_v3.bio"


def test_v3_differs_from_v2_in_the_context_evidence_FIELD_AND_IN_NOTHING_ELSE():
    """Measured field by field over all 5,176 rules, not claimed from how it was generated.

    `Rule` carries twelve fields. The comparison walks every one of them on every rule, in order, and
    collects the names that ever disagree; the answer has to be the single-element set. Asserting the
    set rather than `!=` on one field is deliberate: a comparison that only looks at the field it
    expects to move cannot notice the one it did not expect.
    """
    m2, m3 = parse(V2.read_text()), parse(V3.read_text())
    assert len(m2.rules) == len(m3.rules) == 5176
    names = [f.name for f in dataclasses.fields(Rule)]
    assert "context_evidence" in names and len(names) >= 12, names
    moved = {n for a, b in zip(m2.rules, m3.rules) for n in names if getattr(a, n) != getattr(b, n)}
    assert moved == {"context_evidence"}, sorted(moved)
    # and the field moved on EVERY rule, so the set above is not one rule's accident
    assert all(r.context_evidence == "" for r in m2.rules)
    assert all(r.context_evidence != "" for r in m3.rules)
    from genomeos.attribution import context_evidence as ce

    assert {ce.parse_value(r.context_evidence)[0] for r in m3.rules} <= set(ce.STATES)


def test_v3_and_v2_are_the_same_program_once_the_field_is_dropped():
    """The whole BioIR module and not only the rules: regions, gene stubs, domains, measured blocks,
    parameters, timers, every top-level key. `to_dict` is the form S1 and R1 check a cut in."""
    d2, d3 = parse(V2.read_text()).to_dict(), parse(V3.read_text()).to_dict()
    assert set(d2) == set(d3)
    assert d2 != d3, "identical before the field is dropped would mean v3 carries no reading at all"
    for r in (*d2["rules"], *d3["rules"]):
        r.pop("context_evidence", None)
    assert d2 == d3


def test_every_line_that_differs_between_v2_and_v3_outside_the_rules_is_a_COMMENT():
    """The text claim beside the structural one: v3's extra lines are the generation date, the
    declaration block and the context-evidence summary - all comments. Counted as a multiset, so a
    line that merely moved is not reported as a difference."""
    n2 = [ln for ln in V2.read_text().splitlines() if not ln.startswith("rule ")]
    n3 = [ln for ln in V3.read_text().splitlines() if not ln.startswith("rule ")]
    c2, c3 = collections.Counter(n2), collections.Counter(n3)
    only = [*(c2 - c3).elements(), *(c3 - c2).elements()]
    assert only, "no difference at all would mean v3 is v2 and this lane cut nothing"
    assert [ln for ln in only if ln.strip() and not ln.lstrip().startswith("#")] == []
    # the rule lines are the same lines once the field is taken off, which is the other half
    r2 = [ln for ln in V2.read_text().splitlines() if ln.startswith("rule ")]
    r3 = [ln for ln in V3.read_text().splitlines() if ln.startswith("rule ")]
    assert len(r2) == len(r3) == 5176
    stripped = [re.sub(r"context_evidence: [^;]*; ", "", ln) for ln in r3]
    assert stripped == r2


def _bio_digests() -> dict[str, str]:
    return {
        p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted((ROOT / "data").rglob("*.bio"))
    }


@needs_the_chr21_element_table
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
