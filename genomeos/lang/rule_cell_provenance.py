# SPDX-License-Identifier: AGPL-3.0-or-later
"""PROPOSAL, registered before it is built: carry a compiled rule's CELL PROVENANCE.

A compiled rule says `when: cell_type = K562`. A reader takes that as "where this rule acts". What
it records, on almost every rule in the compiled corpus, is "the biosample name of the track that
held the largest predicted effect of deleting this element on this gene". Those are different
claims and only the second is true of such a rule.

WHY, and the figure is NOT this lane's. `lane-clause1` measured, over its own population of
resolved v2 calls, that every one of those rules' cells is an ARGMAX over the model's track axis
(registration `ba8c41f`, result `b4622ad`, `data/results/clause1.json`). Read its counts there;
they are not restated here and this lane re-derives none of them. Its registration also states the
code fact this proposal rests on, which is checkable without any sweep: `direction_v2.ASSIGNED_CELL`
is `predicted_coding['tissue']`, and `genomeos/predict/enhancer_target.py's `predict_target`` sets that
field to `max_drop_tissue` or `max_rise_tissue` - the tissue of the larger of the biggest drop and
the biggest rise over all tracks, chosen jointly with the gene. `genomeos/attribution/compile.py` then
writes `context(pc.get("tissue"))` into `when: cell_type` on the rule line. So the cell on such a
rule is the place an extreme fell, and the compiled text says nothing that tells a reader so.

WHAT THIS PROPOSAL IS NOT. It recommends nothing, and whether to adopt it is Albert's decision.
It edits no `.bio` file, no compiler, no IR and no grammar; it changes no cell, no action, no
strength and no number, and it is not a route to changing one. A provenance records how a cell got
onto a rule and never which cell should be there. It is also not a gate: see decision 4.

PAIRED WITH, NOT DUPLICATING, THE SIBLING AXIS. `genomeos.lang.rule_evidence_tier` (registration
`709aebc`, proposal `6150aa9`) marks what a rule's NUMBERS rest on. This marks how a rule's CELL
was chosen. They are orthogonal: a rule can carry a measured strength and an argmax cell, or an
author's strength and a screen's cell. The two axes share exactly one thing, the absent-default,
and this module IMPORTS that one object rather than re-spelling it (section 1).
"""

from __future__ import annotations

from typing import Any

from genomeos.lang.rule_evidence_tier import NOT_ASSESSED as _SIBLING_NOT_ASSESSED

# --- 1. THE AXIS ------------------------------------------------------------------------------------
#
# Three classes. Each is here because it OCCURS in the committed corpus and was found by following
# the code that writes `when: cell_type`, not by imagining a taxonomy. The counts land in the
# implementation and in `scripts/rule_cell_provenance_cost.py`; nothing in this section asserts one.

#: The cell is the biosample name of the track on which the largest predicted effect fell. Written by
#: `compile.compile_chromosome` from `predict_target`'s `tissue`.
ARGMAX_OF_PREDICTED_EFFECT = "argmax_of_predicted_effect"

#: The cell is the cell line a CRISPRi screen did the perturbation in. Written by
#: `compile._measured_blocks` from `measured.rule_links`, whose own docstring says "The returned cell
#: is the one the compiler gates the rule on (`when: cell_type = <cell>`)".
MEASURED_PERTURBATION_IN_THAT_CELL = "measured_perturbation_in_that_cell"

#: The cell is asserted by the program's author. The program cites a source beside the rule, but
#: nothing in the program ties that source to the cell, so the cell rests on the author's assertion.
AUTHOR_DECLARED = "author_declared"

#: The three classes, in the order they are defined above and with no order implied between them.
CLASSES: tuple[str, ...] = (
    ARGMAX_OF_PREDICTED_EFFECT,
    MEASURED_PERTURBATION_IN_THAT_CELL,
    AUTHOR_DECLARED,
)

#: The absent-default, and it is NOT a class on the axis: it is the absence of one. "Nobody has said
#: how this rule's cell was chosen" and "the cell rests on the author's word" are different states of
#: the world, and only the second is a finding. Decided by `genomeos.lang.rule_evidence_tier`
#: (`709aebc`) for its own axis and taken from it UNCHANGED and AS THE SAME OBJECT: a second spelling
#: of the word in this file would be a second default that could drift from the first.
#:
#: `is` ALONE IS NOT EVIDENCE HERE and the suite does not rest on it. The value is an identifier-
#: shaped string, so CPython interns it: a plain literal in this file would satisfy `is` for the
#: wrong reason and the assertion would pass while the two defaults were separate declarations. So
#: the suite asserts the identity AND walks this module's own syntax tree to require that no
#: assignment in it binds that string as a literal - the import is the only way the value can arrive.
NOT_ASSESSED = _SIBLING_NOT_ASSESSED

#: Every value the axis can hold: the three classes plus the absence of one. The difference between
#: `len(CLASSES)` and `len(AXIS_VALUES)` is load-bearing.
AXIS_VALUES: tuple[str, ...] = (*CLASSES, NOT_ASSESSED)

#: What each class means, in words a compiled program can carry verbatim. The argmax phrase NEVER
#: says "observed in", "measured in", "acts in" or "active in"; `check_no_observation_claim` refuses
#: a record that makes it say so, and the suite plants each forbidden phrase to show the refusal is
#: real.
PHRASES: dict[str, str] = {
    ARGMAX_OF_PREDICTED_EFFECT: (
        "argmax of predicted deletion effect: the cell is the biosample name of the track on which "
        "the largest predicted effect of deleting this element on this gene fell, selected over the "
        "model's whole track axis jointly with the gene. It is a selection, not an observation in "
        "this cell, and it is not evidence that the element acts here"
    ),
    MEASURED_PERTURBATION_IN_THAT_CELL: (
        "measured perturbation in this cell: a CRISPRi screen silenced this element in this cell line "
        "and measured this gene. The cell is the screen's own, chosen before any value was read"
    ),
    AUTHOR_DECLARED: (
        "declared by the program's author: the cell is asserted on the rule line and the program "
        "states no measurement that fixes it. Searched and nothing found, which is not the same as "
        "never searched"
    ),
}

#: Phrases a provenance record may not contain for the argmax class, each because it converts a
#: selection into an observation. Checked, not trusted.
FORBIDDEN_OF_THE_ARGMAX_CLASS: tuple[str, ...] = (
    "observed in",
    "measured in",
    "acts in",
    "active in",
    "expressed in",
)

#: Precedents in this codebase for marking the absence of a reading rather than defaulting it. The
#: decision in `NOT_ASSESSED` is theirs and is not new here; neither is changed by this proposal.
PRECEDENTS_FOR_AN_ABSENT_READING: dict[str, str] = {
    "rule_evidence_tier.NOT_ASSESSED": (
        "the sibling axis (`709aebc`, `6150aa9`) makes its default a value that is not a tier and "
        "whose `rank()` raises, so that 'nobody looked' cannot be ordered against a real class. This "
        "module imports that object; it does not copy the decision and it does not re-spell the word"
    ),
    "ir.model.Rule.context_evidence": (
        'its own comment: "" where no reading was taken at all, which is not the same as '
        "`not_assessable`, a reading that was attempted and could not be taken. It is the closest "
        "structural precedent for this mark: a field on the rule, written by the compiler, recording "
        "a reading beside the rule, which `applies` does not read and which gates nothing"
    ),
    "measured.CONTEXT_UNKNOWN": (
        "a rule whose cell was never recorded says `cell_type = unknown`, which matches no cell, so "
        "it runs nowhere rather than everywhere. An absent cell is already refused rather than "
        "defaulted, and this proposal leaves that untouched"
    ),
}

# --- 2. THE FOUR DECISIONS, FIXED HERE BEFORE ANYTHING IS BUILT -------------------------------------

DECISIONS: dict[str, str] = {
    "1_the_values": (
        "three classes, each occurring in the committed corpus and each found by following the code "
        "that writes `when: cell_type`: the two rule-emitting sites in "
        "`genomeos/attribution/compile.py` (the predicted element rule and the `<id>_measured` rule) "
        "and the hand-authored programs git tracks. No class is on the axis because it would be "
        "tidy. A fourth kind is CODE-REACHABLE and is deliberately NOT a class: see decision 2."
    ),
    "2_the_population_and_the_absent_default": (
        "the mark is per rule and the population is the rules whose `when` names a cell other than "
        "`measured.CONTEXT_UNKNOWN`. A rule gated on `unknown` names no cell, so there is no cell "
        "whose provenance could be stated; it is excluded from the population rather than given a "
        "value, in the same way the sibling axis narrowed to the rules whose numbers are read. The "
        "default for a rule IN the population is `NOT_ASSESSED`, which is the sibling's object and "
        "is not a class on this axis."
    ),
    "3_what_the_runtime_does": (
        "NOTHING, and that is tested rather than assumed. `Rule.applies` does not read the mark, "
        "`Module.active_rules` does not read it, and the GRN trajectory does not move when every "
        "rule in a program is marked. What a RESULT may do is state the set of provenances it rests "
        "on, so that a cell name cannot travel out of a result as though it were an observation."
    ),
    "4_no_gate": (
        "no gate at any class, and the reason is counted rather than argued. The sibling proposal "
        "(`6150aa9`) found that a gate at its measured tier left none of its rules integrable; the "
        "same count is run here over this axis's own population and reported. A gate that empties "
        "the corpus is not a stricter corpus, it is an unrunnable one."
    ),
}

# --- 3. THE AXIS IS NOT ORDERED ---------------------------------------------------------------------
#
# The sibling's four tiers form a cascade and it gives them a `rank()`, with `not_assessed` raising
# so that absence cannot be ordered against a class. THIS AXIS HAS NO ORDER AT ALL. An argmax cell
# and a screen's cell are not two points on one scale: one is a selection over tracks and the other
# is a choice of cell line made before any value was read, and no number says which is "higher". So
# `rank()` here raises for EVERY value, the three classes included, and there is no `floor()`.
#
# The only aggregate is the SET of provenances a population rests on. `NOT_ASSESSED` in that set is
# reported in it and never dropped: a result that rests partly on rules nobody classified says so.

NO_ORDER = (
    "this axis has no order. `rank` raises for every value including the three classes, and there is "
    "no `floor`: the aggregate over a population is the SET of provenances it rests on, with "
    "`not_assessed` kept in that set rather than dropped from it"
)

# --- 4. WHERE THE MARK IS STORED --------------------------------------------------------------------

STORAGE_OPTIONS: dict[str, str] = {
    "on_the_rule_line_in_the_program": (
        "a `cell_provenance:` property on the rule, parsed into a field on `ir.model.Rule` beside "
        "`context_evidence`, defaulting to `NOT_ASSESSED`. A reader of the program sees the mark "
        "next to the cell. Costs an entry in `genomeos.lang.grammar.BLOCKS['rule']['props']`, a "
        "field on `ir.model.Rule`, two lines in `genomeos/attribution/compile.py`, and a regeneration "
        "of docs/BIOLANG-GRAMMAR.md."
    ),
    "a_side_car_record_keyed_by_rule_id": (
        "a separate file mapping rule id to provenance. Costs no grammar change and no IR change, "
        "and this is the only option that could be adopted with no edit to the language at all."
    ),
    "inside_the_existing_evidence_string": (
        "extend the quoted `evidence:` text the compiler already writes. No grammar change, but the "
        "mark would then be prose inside a free-text field that no code can read as a value, and the "
        "field is already where the misreading happens."
    ),
}

STORAGE_DECISION = (
    "on the rule line in the program. The harm this proposal addresses is a READER of a compiled "
    "program taking `when: cell_type = K562` as an observation, and a side-car fixes nothing about "
    "that: the program text is unchanged and still reads as an observation. The mark has to sit "
    "beside the cell it qualifies, in the artefact the reader reads. `Rule.context_evidence` is the "
    "precedent for exactly this shape and it is already in the grammar."
)

STORAGE_LIMIT_OF_THE_CHOSEN_OPTION = (
    "it cannot be adopted without editing the language, and this lane may edit none of it. Two "
    "refusals in the code as it stands make that concrete and both are demonstrated in the suite "
    "rather than asserted: `ir.model.Rule` is a `@dataclass(slots=True)`, so the mark cannot be "
    "attached to a Rule instance at all; and `lang.parser._check_keys` refuses a rule property the "
    "grammar table does not name, so `cell_provenance:` in a program is a syntax error today. "
    "Adopting this therefore means a commit to `genomeos/ir/model.py`, `genomeos/lang/grammar.py` "
    "and `genomeos/attribution/compile.py`, and that commit is Albert's."
)

# --- 5. PREDICTIONS, REGISTERED BEFORE THEY ARE TESTED ----------------------------------------------
#
# Each is falsifiable from the code or from a committed program. A prediction that comes out false is
# reported as false and the decision it supported is revisited in the open.

PREDICTIONS: tuple[dict[str, str], ...] = (
    {
        "id": "P1",
        "claim": (
            "`genomeos/attribution/compile.py` has exactly two sites that emit a `rule` line, and the "
            "class of the cell each writes is fixed by the site: the predicted element rule writes "
            "`predict_target`'s argmax tissue, and the `<id>_measured` rule writes the cell "
            "`measured.rule_links` took from the screen. No third route puts a cell on a compiled "
            "rule."
        ),
        "falsified_by": "a third emitting site, or either site taking its cell from somewhere else",
    },
    {
        "id": "P2",
        "claim": (
            "in the one compiled program git tracks, every rule's class is separable from the rule "
            "ALONE by `Rule.evidence.kind`: PREDICTED for the argmax class, EXPERIMENTAL for the "
            "measured class, and no other kind occurs. This is what makes re-marking an existing "
            "program mechanical rather than a human judgement, and it also means the mark adds "
            "WORDS and not a distinction the file lacked - which is a cost to the proposal and is "
            "reported as one."
        ),
        "falsified_by": "a compiled rule with a third evidence kind, or a kind that maps to both classes",
    },
    {
        "id": "P3",
        "claim": (
            "`author_declared` occurs in the hand-authored programs git tracks, with a count of at "
            "least one and far below the compiled corpus's. It is a real class and not a synonym for "
            "`not_assessed`, by the sibling's own distinction: searched and nothing found is a "
            "finding, never searched is not."
        ),
        "falsified_by": "no hand-authored rule carrying a cell, which would make the class invented",
    },
    {
        "id": "P4",
        "claim": (
            "`cell_type = unknown` is reachable from `compile.context` and occurs ZERO times in the "
            "committed compiled program. The class is therefore excluded from the axis as decision 2 "
            "says, and the exclusion is reported with its count rather than left implicit."
        ),
        "falsified_by": "a committed compiled rule gated on `unknown`, which would need a decision",
    },
    {
        "id": "P5",
        "claim": (
            "the runtime does nothing with the mark. Marking every rule of a program with every value "
            "in turn leaves `Module.active_rules` identical for each context and leaves the GRN "
            "trajectory identical to the unmarked run, value by value."
        ),
        "falsified_by": "any context whose active-rule set moves, or any trajectory that moves",
    },
    {
        "id": "P6",
        "claim": (
            "a gate at `measured_perturbation_in_that_cell` leaves under one per cent of the "
            "committed compiled corpus's rules integrable. The sibling found its own gate left none "
            "of its population; this is the same question asked of this axis and its own population, "
            "and the sibling's figure is cited and not reused."
        ),
        "falsified_by": "a measured share at or above one per cent",
    },
    {
        "id": "P7",
        "claim": (
            "both refusals in STORAGE_LIMIT_OF_THE_CHOSEN_OPTION are real: `setattr` of the mark on a "
            "parsed `Rule` raises, and a program carrying `cell_provenance:` on a rule fails to parse."
        ),
        "falsified_by": "either one succeeding, which would make adoption cheaper than stated",
    },
    {
        "id": "P8",
        "claim": (
            "the parser synthesises rules no program writes (the sibling measured this on its own "
            "program, `6150aa9`; it is cited, not restated), and EVERY synthesised rule names no "
            "cell. The synthesis therefore does not change this proposal's population, which is the "
            "opposite of what it did to the sibling's."
        ),
        "falsified_by": "one synthesised rule carrying a `when` that names a cell",
    },
)

# --- 6. WHAT THIS LANE HAD ALREADY SEEN, DECLARED BEFORE IT COUNTS ANYTHING -------------------------

NOT_BLIND = (
    "before this registration was written this lane had read `genomeos/attribution/compile.py`'s two "
    "rule-emitting sites, `predict.enhancer_target.predict_target`, `measured.rule_links`, "
    "`ir.model.Rule`, `lang.parser._check_keys`, `genomeos.lang.rule_evidence_tier` and "
    "`data/results/clause1_registration.json`; it had counted `rule` lines and `when: cell_type` "
    "labels in `data/organisms/human/noncoding_chr21.bio` with grep, and it had parsed that program "
    "and the tracked hand-authored programs once each to see which rules carry a cell. So the three "
    "classes and the predictions above were chosen KNOWING roughly what the counts would be, and "
    "every prediction is written to be falsifiable by a re-count rather than confirmed by one."
)

EXCLUSIONS: tuple[str, ...] = (
    "no `.bio` file is edited, and no compiler, IR, grammar or parser file is edited",
    "no cell, action, strength, threshold or confidence on any rule is changed or proposed changed",
    "no gate is added and none is proposed; decision 4 refuses one and counts the reason",
    "no figure from `b4622ad` or `6150aa9` is restated as this lane's; both are cited",
    "no re-derivation of any rule's argmax: that would need the deletion sweep, which is a git-"
    "ignored store and a heavier read than this proposal needs. The CLASS of a rule's cell is "
    "established from the code that wrote it and from the rule's own evidence kind; whether a given "
    "cell is in fact the extreme over the track axis is `lane-clause1`'s measurement over its own "
    "population and is cited as that",
    "no model request, no money, no network",
)

#: Readings this lane carries from peers WORD FOR WORD and does not strengthen.
CARRIED_READINGS: dict[str, str] = {
    "withholds_never_reverses": (
        "v2 can WITHHOLD a direction and can never reverse one. It is cited here only to say what it "
        "does NOT touch: v2 decides a rule's ACTION and never its cell, so no v2 outcome adds, "
        "removes or changes a cell provenance"
    ),
    "unresolved_is_not_absence": (
        "an unresolved class is NOT absence of regulation. It never means the element has no action "
        "and it is not a measurement of absence"
    ),
    "two_track_cell": (
        'a two-track cell\'s tracks are "not known to be biological replicates", in '
        "`enhancer_target._cell_summary`'s own words, so it is not a validated cell-level effect. "
        "Nothing in this proposal makes them replicates and no mark here is evidence either way"
    ),
    "one_biosample_name": (
        "the track metadata copy carries no biosample accession, so ONE BIOSAMPLE NAME is the most "
        "that can be said. The argmax phrase in PHRASES therefore says `biosample name` and never "
        "names a biosample"
    ),
}

REPORTS: tuple[str, ...] = (
    "the classes found occurring, with the count of each in the committed corpus and how each was "
    "established",
    "the code-reachable class that occurs zero times, with its count",
    "what the runtime was TESTED to do with the mark, named test by test",
    "the gate count on this axis's own population",
    "the adoption cost, split into marks a human must decide and marks a machine can derive",
    "every prediction's outcome, false ones included, and anything the code contradicts",
)

WHICHEVER_WAY_IT_FALLS = (
    "if the classes come out as predicted, the proposal is a per-rule mark whose compiled part costs "
    "no human judgement. If P2 holds, the same finding says the mark adds WORDS and not a "
    "discrimination the compiled file lacked, which a reader may count against adopting it, and it "
    "is reported in those terms. If P3 is falsified the third class is struck. If P5 is falsified "
    "the proposal is no longer inert and this registration is wrong about its own central claim. No "
    "outcome here produces a recommendation."
)

NO_RECOMMENDATION = (
    "this module makes no recommendation about adopting the proposal, and neither will the result. "
    "Implementing it IS the adoption, and that is Albert's."
)


def registration() -> dict[str, Any]:
    """The whole registration as data, so a result can carry it and be compared against it."""
    return {
        "proposal": "carry a compiled rule's cell provenance",
        "lane": "lane-cellprov",
        "why_cited_not_restated": {
            "argmax_measured_by_a_peer": "ba8c41f (registration), b4622ad (result)",
            "sibling_axis": "709aebc (registration), 6150aa9 (proposal)",
            "committed_absence_is_a_defect": "4f44dbf",
        },
        "classes": list(CLASSES),
        "axis_values": list(AXIS_VALUES),
        "not_assessed": NOT_ASSESSED,
        "not_assessed_is_the_default_and_is_not_a_class": True,
        "not_assessed_is_the_sibling_object": True,
        "phrases": dict(PHRASES),
        "forbidden_of_the_argmax_class": list(FORBIDDEN_OF_THE_ARGMAX_CLASS),
        "precedents_for_an_absent_reading": dict(PRECEDENTS_FOR_AN_ABSENT_READING),
        "decisions": dict(DECISIONS),
        "no_order": NO_ORDER,
        "storage_options": dict(STORAGE_OPTIONS),
        "storage_decision": STORAGE_DECISION,
        "storage_limit_of_the_chosen_option": STORAGE_LIMIT_OF_THE_CHOSEN_OPTION,
        "predictions": [dict(p) for p in PREDICTIONS],
        "not_blind": NOT_BLIND,
        "exclusions": list(EXCLUSIONS),
        "carried_readings": dict(CARRIED_READINGS),
        "reports": list(REPORTS),
        "whichever_way_it_falls": WHICHEVER_WAY_IT_FALLS,
        "no_recommendation": NO_RECOMMENDATION,
        # written in as empty, so that the commit which fills them is visibly later than this one
        "adoption_cost": {},
        "prediction_outcomes": {},
    }


# ====================================================================================================
# IMPLEMENTATION. Everything above this line was committed at `ea3ef2e`, in a commit of its own,
# before any of the following existed and before a single rule of any program was classified. Nothing
# below changes a decision above; where a prediction came out false it is recorded as false in
# `scripts/rule_cell_provenance_cost.py`'s output and reported, not quietly amended.
# ====================================================================================================

#: THIS MODULE IS AN UNADOPTED PROPOSAL, and it sits inside the Apache-2.0 engine because
#: `genomeos/lang/` is where a language proposal about rules belongs. Nothing in the engine calls it:
#: its only callers are `tests/test_rule_cell_provenance.py` and
#: `scripts/rule_cell_provenance_cost.py`, both application-side. Adopting the proposal means a
#: commit to `genomeos/ir/model.py`, `genomeos/lang/grammar.py` and `genomeos/attribution/compile.py`,
#: and that commit is Albert's.
#:
#: WHETHER AN UNADOPTED PROPOSAL SHOULD SHIP INSIDE THE PACKAGED ENGINE at the 2.0 milestone is NOT
#: decided here and is not this module's to decide: it is a question for Albert's package. What is
#: recorded here is only that it is unadopted, so the question can be asked of a file that says so.
#: One measured fact belongs with it: this file carries an `AGPL-3.0-or-later` SPDX header and
#: `scripts/package_engine.py` copies `genomeos/lang/` into the Apache package, as do the two
#: neighbours it is built on. That is reported, not resolved.
UNADOPTED = (
    "an unadopted proposal. No engine module calls it; adopting it is a commit to the IR, the grammar "
    "and the compiler, and that commit is Albert's. Whether it should ship inside the packaged engine "
    "is a question for the package and is not decided here"
)

import ast  # noqa: E402
import dataclasses  # noqa: E402
import inspect  # noqa: E402
from collections.abc import Iterable  # noqa: E402

#: The quoted `evidence:` source a compiled PREDICTED rule carries, as `attribution.compile` writes
#: it: `predicted "AlphaGenome deletion, <tissue>"`. Matched on the prefix, because the tissue name
#: follows it on every rule and is the very thing whose provenance is in question.
PREDICTED_SOURCE_PREFIX = "AlphaGenome deletion"


@dataclasses.dataclass(frozen=True)
class CorpusVocabulary:
    """The two application-side strings this module needs, INJECTED and never imported.

    THE LICENSING BOUNDARY IS WHY THIS EXISTS, and it is worth writing down because two of this
    project's rules collided here. "Import a value, never copy its literal" said to read the CRISPRi
    citation from the application module that owns it. `genomeos/lang/` is part of the Apache-2.0
    engine, so doing that made the engine import the AGPL-3.0 application, against LICENSING.md
    decision D40 (Albert's split of 2026-09-11) - a legal boundary, not a style rule, and
    `tests/test_engine_boundary.py` failed on it exactly as it should have. Both rules are right, so
    the resolution is structural rather than a choice between them: the engine takes the values as a
    PARAMETER and the application-side caller supplies them from the module that owns them. Nothing
    is copied and nothing crosses the boundary. The boundary test is untouched: no exemption, no
    allowlist, no per-file skip.

    `unrecorded_cell` is the application's `measured.CONTEXT_UNKNOWN` and `measured_source_prefix` is
    its `measured.SOURCES["crispri"]`. `scripts/rule_cell_provenance_cost.py` is the one place that
    reads them, and a test there asserts each value IS that object, so "import, never copy" still
    holds where the import is legal.
    """

    unrecorded_cell: str
    measured_source_prefix: str

    def __post_init__(self) -> None:
        for name in ("unrecorded_cell", "measured_source_prefix"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(
                    f"{name} must be a non-empty string supplied by the application side; got "
                    f"{value!r}. This module may not import it: see the class docstring and D40"
                )


class CellProvenanceUndecidableError(ValueError):
    """The rule's cell provenance cannot be read off the rule. Raised, never guessed.

    This exception IS the adoption cost. Every rule that raises it is a rule a human would have to
    adjudicate; every rule that does not is one a machine can mark. `scripts/rule_cell_provenance_cost.py`
    counts the two populations by running this function over the committed corpus, so the cost is
    measured by the code rather than estimated.
    """


def is_class(value: str) -> bool:
    """True for the three provenance classes, False for the absent-default and for anything off the axis."""
    return value in CLASSES


def rank(value: str) -> int:
    """ALWAYS raises, for every value including the three classes. See NO_ORDER.

    The sibling axis (`6150aa9`) is a cascade and ranks its tiers, with its default raising. This
    axis is not a cascade: an argmax over tracks and a screen's choice of cell line are not two
    points on one scale, and nothing in the project says which is "higher". A rank would invite a
    comparison the axis cannot support, so the function exists and refuses rather than being absent
    and leaving a reader to write their own.
    """
    raise TypeError(
        f"cell provenance is not ordered, so {value!r} has no rank. {NO_ORDER}. "
        f"Use `provenances_relied_on` for the set."
    )


def provenances_relied_on(marks: Iterable[str]) -> list[str]:
    """The SET of provenances a population rests on, sorted, with the absent-default KEPT in it.

    The absent-default is not dropped and is not collapsed into a class: a result that rests partly
    on rules nobody classified says so in the same field that names the classes. An unknown value
    raises; a silent pass-through would let a typo travel as a class.
    """
    out = set()
    for m in marks:
        if m not in AXIS_VALUES:
            raise ValueError(f"{m!r} is not on the axis; the values are {list(AXIS_VALUES)}")
        out.add(m)
    return sorted(out)


def check_no_observation_claim(records: Iterable[dict[str, Any]]) -> None:
    """Refuse a provenance record whose argmax phrasing claims an observation.

    Modelled on `attribution.increase_links.check_no_mechanism_claim`, for the same reason: the one
    thing this mark exists to prevent is the sentence it could so easily be written as. A record for
    the argmax class may not contain "observed in", "measured in", "acts in", "active in" or
    "expressed in" anywhere in its text, whatever field carries it.
    """
    for rec in records:
        if rec.get("cell_provenance") != ARGMAX_OF_PREDICTED_EFFECT:
            continue
        text = " ".join(str(v) for v in rec.values()).lower()
        for bad in FORBIDDEN_OF_THE_ARGMAX_CLASS:
            if bad in text:
                raise ValueError(
                    f"a record of class {ARGMAX_OF_PREDICTED_EFFECT} says {bad!r}, which turns a "
                    f"selection over the model's track axis into an observation in the cell. "
                    f"{PHRASES[ARGMAX_OF_PREDICTED_EFFECT]}"
                )


def names_a_cell(rule: Any, vocabulary: CorpusVocabulary) -> bool:
    """Whether the rule is IN the population: its `when` names a cell other than the unrecorded one.

    A rule gated on the application's `measured.CONTEXT_UNKNOWN` names no cell, so there is no cell
    whose provenance could be stated (decision 2). The runtime already refuses such a rule in every
    context - `Rule.applies` returns False on it - and this proposal leaves that untouched. The value
    arrives in `vocabulary` because this module may not import it; see `CorpusVocabulary`.
    """
    cell = (rule.when or {}).get("cell_type")
    return bool(cell) and cell != vocabulary.unrecorded_cell


def provenance_of(rule: Any, vocabulary: CorpusVocabulary) -> str:
    """The rule's cell provenance read off the rule itself, or `CellProvenanceUndecidableError`.

    Only the two signatures `attribution.compile` writes are decidable from a rule alone, and they
    are decided on the evidence KIND together with the quoted SOURCE - not on the kind alone, because
    a hand-authored rule can carry `experimental` with a citation of its own and must not be mistaken
    for a screen's measurement. Everything else raises; the caller adjudicates and says that it did.
    """
    if not names_a_cell(rule, vocabulary):
        raise CellProvenanceUndecidableError(
            f"rule {rule.id!r} names no cell (when={rule.when!r}), so it is outside the population"
        )
    kind = getattr(rule.evidence.kind, "value", rule.evidence.kind)
    source = rule.evidence.source or ""
    if kind == "predicted" and source.startswith(PREDICTED_SOURCE_PREFIX):
        return ARGMAX_OF_PREDICTED_EFFECT
    if kind == "experimental" and source.startswith(vocabulary.measured_source_prefix):
        return MEASURED_PERTURBATION_IN_THAT_CELL
    raise CellProvenanceUndecidableError(
        f"rule {rule.id!r} carries evidence {kind!r} {source[:60]!r}, which matches neither "
        f"signature the compiler writes. A human must adjudicate this cell's provenance"
    )


class MarkedRule:
    """A rule carrying `cell_provenance`, built WITHOUT editing `ir.model.Rule`.

    `ir.model.Rule` is a `@dataclass(slots=True)`, so the mark cannot be set on an instance of it:
    that is prediction P7 and the suite demonstrates the refusal. A subclass declaring no `__slots__`
    of its own gets a `__dict__` and accepts the attribute, which is how the suite can run the real
    runtime over marked rules without a line of the IR changing. It is a TEST INSTRUMENT and not the
    proposed storage; the proposed storage is a field on `Rule` itself, which is Albert's commit.
    """

    def __init_subclass__(cls, **kw: Any) -> None:  # pragma: no cover - guard only
        raise TypeError("MarkedRule is built by `mark` and is not a base class")


def _marked_class() -> type:
    """The `Rule` subclass that accepts the attribute, built once and reused."""
    from genomeos.ir.model import Rule

    global _MARKED
    try:
        return _MARKED
    except NameError:
        _MARKED = type("Rule_with_cell_provenance", (Rule,), {"__doc__": MarkedRule.__doc__})
        return _MARKED


def mark(rule: Any, provenance: str) -> Any:
    """A copy of `rule` carrying `cell_provenance`, field for field identical otherwise."""
    if provenance not in AXIS_VALUES:
        raise ValueError(f"{provenance!r} is not on the axis; the values are {list(AXIS_VALUES)}")
    fields = {f.name: getattr(rule, f.name) for f in dataclasses.fields(rule)}
    out = _marked_class()(**fields)
    out.cell_provenance = provenance
    return out


def provenance_mark(rule: Any) -> str:
    """The mark a rule carries, or the absent-default. An unmarked rule reads as `not_assessed`."""
    return getattr(rule, "cell_provenance", NOT_ASSESSED)


def gate_survivors(rules: Iterable[Any], at: str) -> list[Any]:
    """The rules a hypothetical gate at one class would leave integrable.

    No gate is proposed (decision 4). This exists so the cost of one can be COUNTED on this axis's
    own population rather than argued about, the way `6150aa9` counted it on the sibling's.
    """
    if at not in CLASSES:
        raise ValueError(f"a gate can only be set at a class, not at {at!r}; classes are {list(CLASSES)}")
    return [r for r in rules if provenance_mark(r) == at]


def binds_as_a_literal(value: str) -> list[int]:
    """The line numbers in THIS module where an assignment binds `value` as a string literal.

    Walking the syntax tree rather than grepping, because these words appear in prose throughout this
    file and prose is not a declaration. Only an assignment whose value IS the literal counts.

    It polices two different rules with one instrument. For the absent-default the import is the only
    way the value may arrive (see NOT_ASSESSED), and `is` cannot show that on its own because the
    value is an identifier-shaped string CPython interns. For an injected application string the
    module must hold no copy at all, and there `is` IS evidence - but the walk is kept anyway, since a
    stale copy sitting unused in the engine would still be a copy of AGPL-side text.
    """
    tree = ast.parse(inspect.getsource(inspect.getmodule(provenance_of)))
    out = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            bound = node.value
            if isinstance(bound, ast.Constant) and bound.value == value:
                out.append(node.lineno)
    return out


def binds_the_default_as_a_literal() -> list[int]:
    """`binds_as_a_literal` for the absent-default. Must be empty."""
    return binds_as_a_literal(NOT_ASSESSED)
