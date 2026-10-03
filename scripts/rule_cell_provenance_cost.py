#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""The adoption cost of the cell-provenance proposal, counted over the committed corpus.

The proposal is registered at `ea3ef2e` (`genomeos.provenance.rule_cell_provenance`). It recommends
nothing and this script recommends nothing: it counts, because a proposal that hides its own scale
is not a proposal.

Two numbers matter and they are counted separately, because they cost different things:

  * marks a MACHINE can derive. A compiled rule's class is readable from the rule alone, on the
    evidence kind together with the quoted source, so re-marking an existing compiled program is
    mechanical. `rule_cell_provenance.provenance_of` does it and raises rather than guessing.
  * marks a HUMAN must decide. Every rule on which that function raises. In the committed corpus
    these are the hand-authored rules that carry a cell: the program asserts a cell and states no
    measurement that fixes it, so only a person can say what the cell rests on.

WHAT IS COUNTED AND WHAT IS NOT. Only what git tracks. `data/organisms/human/noncoding_chr21.bio`
is the one compiled program in the tree; the genome-wide compiled corpus is not committed, so no
figure for it is produced or guessed here. A count for it exists in another lane's note and is not
this lane's to restate.

0 `.bio` files written, 0 cells changed, 0 model requests, no money, no network.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from genomeos.lang.parser import parse_file  # noqa: E402
from genomeos.provenance import rule_cell_provenance as cp  # noqa: E402

#: The generated program git tracks. Named, not globbed, so a rename is a failure and not a silent
#: change of population.
COMPILED = "data/organisms/human/noncoding_chr21.bio"


def vocabulary() -> cp.CorpusVocabulary:
    """The application-side strings the engine module takes as a parameter. THE ONLY PLACE they cross.

    `rule_cell_provenance` used to live at `genomeos/lang/rule_cell_provenance.py`, part of the
    Apache-2.0 engine, which may not import the AGPL-3.0 application (LICENSING.md decision D40,
    Albert's split of 2026-09-11). It did import it, in order to honour the other rule this project
    holds - import a value, never copy its literal - and `tests/test_engine_boundary.py` failed on it,
    correctly. The two rules were reconciled here and not traded against each other: this file is
    application-side, so the import is legal, and it hands the values to the census as data. Nothing is
    copied.

    Since 2026-10-03 the census is application-side too (`genomeos/provenance/`), so D40 no longer
    forbids it the import. The injection stays: it is kept for the reason given in `CorpusVocabulary`'s
    own docstring - a corpus with a different vocabulary is then a different argument - and this
    function remains the one place the two values cross.

    `tests/test_rule_cell_provenance.py` asserts that each value IS the object this reads, by identity
    against `measured.SOURCES["crispri"]` and `measured.CONTEXT_UNKNOWN`, and - because one of those
    is an identifier-shaped string that CPython interns, where identity would pass for the wrong
    reason - also that no assignment in the engine module binds either string as a literal.
    """
    from genomeos.attribution import measured

    return cp.CorpusVocabulary(
        unrecorded_cell=measured.CONTEXT_UNKNOWN,
        measured_source_prefix=measured.SOURCES["crispri"],
    )


#: Every registered prediction's outcome, with the figure that settled it and the test that measured
#: it. They live here and NOT in `rule_cell_provenance.registration()`, which was committed with
#: these fields empty so that the commit filling them is visibly later than the commit registering
#: them. A prediction that had come out false would say so here in the same words.
OUTCOMES: dict[str, dict[str, str]] = {
    "P1": {
        "verdict": "held",
        "figure": (
            "exactly two functions in `attribution/compile.py` emit a `rule` line: "
            "`compile_chromosome` (the predicted element rule) and `_measured_blocks` (the "
            "`<id>_measured` rule). Counted by walking the compiler's syntax tree, not by grep"
        ),
        "test": "test_the_compiler_emits_a_rule_line_from_exactly_two_sites",
    },
    "P2": {
        "verdict": "held",
        "figure": (
            "all 5,176 rules of the one committed compiled program classify from the rule alone, "
            "5,174 argmax and 2 measured, and the evidence kind partitions them exactly: predicted "
            "maps only to the argmax class and experimental only to the measured one"
        ),
        "test": "test_every_rule_of_the_committed_compiled_program_classifies_from_the_rule_alone",
    },
    "P3": {
        "verdict": "held, and narrowed by its own smallness",
        "figure": (
            "the whole tracked hand-authored corpus carries ONE rule that names a cell, "
            "`Nkx2-5 activates MYH6` in `data/demo/cell_context.bio`. `provenance_of` raises on it, "
            "and this lane adjudicates it `author_declared` because the program cites a source for "
            "the factor and nothing in it ties that source to the cell. So the class occurs and is "
            "not invented, but it rests on a single adjudication and is reported as resting on one"
        ),
        "test": "test_the_hand_authored_corpus_carries_cell_gated_rules_that_no_signature_classifies",
    },
    "P4": {
        "verdict": "held",
        "figure": (
            "0 of 5,176 committed compiled rules are gated on `measured.CONTEXT_UNKNOWN`, so the "
            "code-reachable fourth kind is unoccupied and stays off the axis as decision 2 says"
        ),
        "test": "test_no_rule_of_the_committed_compiled_program_is_gated_on_an_unrecorded_cell",
    },
    "P5": {
        "verdict": "held",
        "figure": (
            "for each of the 4 axis values in turn, marking every rule of a committed program left "
            "`Module.active_rules` identical in 3 contexts and left the GRN trajectory identical "
            "level for level. `applies` matched on the `when` clause and ignored the mark"
        ),
        "test": "test_the_trajectory_does_not_move_when_every_rule_is_marked",
    },
    "P6": {
        "verdict": "held",
        "figure": (
            "a gate at `measured_perturbation_in_that_cell` leaves 2 of 5,176 rules integrable, "
            "0.04 per cent. The sibling's own gate figure is cited (`6150aa9`) and not reused"
        ),
        "test": "test_a_gate_at_the_measured_class_leaves_under_one_per_cent_integrable",
    },
    "P7": {
        "verdict": "held",
        "figure": (
            "both refusals are real. A parsed `Rule` has no `__dict__` and `setattr` of the mark "
            "raises AttributeError; a program writing `cell_provenance:` on a rule raises "
            "BioLangError `has no key 'cell_provenance'`, and the key is absent from "
            "`grammar.BLOCKS['rule']['props']` while `context_evidence` is present"
        ),
        "test": "test_a_program_that_writes_the_property_on_a_rule_does_not_parse_today",
    },
    "P8": {
        "verdict": "held",
        "figure": (
            "the parser synthesises 119 rules no tracked program writes, and 0 of the 119 name a "
            "cell. The compiled program synthesises none of its 5,176. So the synthesis leaves this "
            "proposal's population untouched, which is the opposite of what it did to the sibling's"
        ),
        "test": "test_every_rule_the_parser_synthesises_names_no_cell",
    },
}

#: Found while testing, covered by no prediction, and reported rather than worked around.
FOUND_WITHOUT_A_PREDICTION: tuple[str, ...] = (
    "the compiler writes the misreading itself. Its emitted header says every rule is `gated on the "
    "cell it was measured or predicted in`, and on the 5,174 argmax rules the cell is neither: it is "
    "where an extreme fell. Adopting the proposal therefore also means changing that sentence, which "
    "no part of the registration had counted",
    "the class is already separable from the rule today, on `evidence.kind` plus the quoted source. "
    "So the mark adds the SENTENCE and not a distinction the file lacked, and that is a cost to the "
    "proposal rather than a point for it",
    "`context_evidence`, the precedent this proposal copies its shape from, is the empty string on "
    "every one of the 5,176 committed compiled rules and the word does not appear in the file at "
    "all: the mapping table was unavailable when the program was generated, so no reading was taken. "
    "A field that is absent everywhere it occurs is the precedent's own warning about this one",
    "the absent-default cannot be assured by `is` alone. The value is an identifier-shaped string "
    "that CPython interns, so a plain literal in this module would satisfy `is` for the wrong "
    "reason; the identity is therefore backed by a syntax-tree walk requiring that no assignment in "
    "the module binds the string",
)


def _tracked_bio() -> list[str]:
    out = subprocess.run(["git", "ls-files", "*.bio"], cwd=ROOT, capture_output=True, text=True, check=True)
    return [line for line in out.stdout.splitlines() if line]


def _written_headers(path: Path) -> set[str]:
    out = set()
    for line in path.read_text().splitlines():
        s = line.strip()
        if s.startswith("rule "):
            out.add(s[len("rule ") :].split("{")[0].strip())
    return out


def count() -> dict[str, object]:
    """Every figure this script reports, as plain data."""
    machine: dict[str, int] = {}
    human: list[dict[str, str]] = []
    synthesised = 0
    synthesised_naming_a_cell = 0
    outside_population = 0
    rules_total = 0
    programs_with_a_cell_gated_rule: list[str] = []
    vocab = vocabulary()

    for relative in _tracked_bio():
        path = ROOT / relative
        if not path.exists():
            raise AssertionError(
                f"{relative} is tracked by git and absent from this checkout. A committed "
                f"artefact's absence is a defect, not a machine difference"
            )
        module = parse_file(path)
        written = _written_headers(path)
        has_cell = False
        for rule in module.rules:
            rules_total += 1
            if rule.id not in written:
                synthesised += 1
                if cp.names_a_cell(rule, vocab):
                    synthesised_naming_a_cell += 1
            if not cp.names_a_cell(rule, vocab):
                outside_population += 1
                continue
            has_cell = True
            try:
                klass = cp.provenance_of(rule, vocab)
            except cp.CellProvenanceUndecidableError:
                human.append({"program": relative, "rule": rule.id})
                continue
            machine[klass] = machine.get(klass, 0) + 1
        if has_cell:
            programs_with_a_cell_gated_rule.append(relative)

    in_population = sum(machine.values()) + len(human)
    return {
        "programs_scanned": len(_tracked_bio()),
        "rules_parsed": rules_total,
        "rules_outside_the_population": outside_population,
        "rules_in_the_population": in_population,
        "programs_with_a_cell_gated_rule": programs_with_a_cell_gated_rule,
        "marks_a_machine_can_derive": sum(machine.values()),
        "marks_a_machine_can_derive_by_class": dict(sorted(machine.items())),
        "marks_a_human_must_decide": len(human),
        "marks_a_human_must_decide_rules": human,
        "default_for_every_one_of_them_before_a_decision": cp.NOT_ASSESSED,
        "synthesised_rules": synthesised,
        "synthesised_rules_naming_a_cell": synthesised_naming_a_cell,
        "compiled_program": COMPILED,
        "the_genome_wide_corpus": ("not committed, so no figure for it is produced here and none is guessed"),
        "an_adoption_cost_outside_the_rules": (
            "the compiler's own emitted header says every rule is `gated on the cell it was measured "
            "or predicted in`. On the argmax rules that sentence IS the misreading, so adopting the "
            "proposal means changing that line in `genomeos/attribution/compile.py` as well as "
            "adding the field"
        ),
    }


def main(argv: list[str] | None = None) -> int:
    argparse.ArgumentParser(description=__doc__).parse_args(argv)
    c = count()
    print(f"programs scanned (git-tracked .bio): {c['programs_scanned']}")
    print(f"rules parsed: {c['rules_parsed']}")
    print(f"  outside the population (name no cell): {c['rules_outside_the_population']}")
    print(f"  in the population: {c['rules_in_the_population']}")
    print(f"marks a machine can derive: {c['marks_a_machine_can_derive']}")
    for klass, n in c["marks_a_machine_can_derive_by_class"].items():  # type: ignore[union-attr]
        print(f"  {klass}: {n}")
    print(f"marks a human must decide: {c['marks_a_human_must_decide']}")
    for row in c["marks_a_human_must_decide_rules"]:  # type: ignore[union-attr]
        print(f"  {row['program']}: {row['rule']}")
    default = c["default_for_every_one_of_them_before_a_decision"]
    print(f"default for every one of them before a decision: {default}")
    print(f"rules the parser synthesises: {c['synthesised_rules']}")
    print(f"  of those, naming a cell: {c['synthesised_rules_naming_a_cell']}")
    print(f"genome-wide: {c['the_genome_wide_corpus']}")
    print(f"beyond the rules: {c['an_adoption_cost_outside_the_rules']}")
    print("prediction outcomes:")
    for pid, o in OUTCOMES.items():
        print(f"  {pid} {o['verdict']}: {o['figure']}")
    print("found without a prediction:")
    for line in FOUND_WITHOUT_A_PREDICTION:
        print(f"  - {line}")
    print("0 .bio files written, 0 cells changed, 0 model requests, no money, no recommendation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
