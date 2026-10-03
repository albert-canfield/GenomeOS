# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""PROPOSAL, registered before it is built: mark rule evidence tier in BioLang programs.

This module is the registration. It fixes the vocabulary, where the mark lives, what the runtime
does with it, the default, and the predictions the proposal will be tested against - all BEFORE any
of it is implemented and before any corpus count is written down. The implementation and the counts
land in later commits; nothing here asserts a figure.

WHY, and it is not this lane's finding: the source census of `data/demo/gastrulation.bio`
(`data/results/rule_number_sources_census.json`, registration `4392d24`, result `4d4003d`) reports
that none of that program's rule numbers has a measured value with a quoted source. Read its own
counts there; they are not restated here. What follows from them is the problem this proposal
addresses: the distinction between a rule number resting on a measurement, one resting on a citation
that supports only the interaction, and one resting on nothing lives TODAY ONLY IN A CENSUS RESULT.
In the program it is invisible, and so it is invisible in every dynamic computed from the program.

WHAT THIS PROPOSAL IS NOT. It recommends nothing. Adopting it changes every hand-authored program,
so whether to adopt it is Albert's decision and this module does not make it or argue for it. It
edits no `.bio` file, it changes no number, no threshold, no strength and no Hill coefficient, and
it is not a route to tuning one: a tier records what a number rests on and never what it should be.
This area has already refused to tune two free numbers toward unsourced targets (`8bb9123`,
`19423bd`) and that refusal is inherited here unchanged.
"""

from __future__ import annotations

from typing import Any

from genomeos.provenance import rule_number_sources as rns

# --- 1. THE AXIS -----------------------------------------------------------------------------------
#
# The census's four classes, REUSED AS THE SAME OBJECTS rather than re-spelled. `TIERS is
# rns.CASCADE` and `TIER_DEFINITIONS is rns.CLASSES` are true, and the test suite asserts both with
# `is`, so a change to the census's vocabulary cannot leave a second divergent copy behind in the
# language. There is no second axis of classes and no renaming.

#: The four tiers, lowest-information last, in the census's own cascade order. Same tuple object.
TIERS: tuple[str, ...] = rns.CASCADE

#: What each tier means. Same dict object as the census's `CLASSES`.
TIER_DEFINITIONS: dict[str, str] = rns.CLASSES

#: The three numbers a `rule` declaration carries, one slot each. Same tuple object.
TIERED_FIELDS: tuple[str, ...] = rns.RULE_FIELDS

#: The census's SECOND axis, what the program's own citation does for the declaration it sits on.
#: Same dict object. This axis carries the rule-level mark; see decision 2.
INTERACTION_STANDING: dict[str, str] = rns.ATTRIBUTION

#: The value an unmarked slot reads, and the DEFAULT. It is not a tier: it is the absence of one.
#: "Searched and nothing was found" (`U_unsourced`) and "never searched" are different states of the
#: world, and only the first is a finding. The census covers ONE program; every other program in the
#: corpus has never been searched, so a default of `U_unsourced` would make every unsearched slot in
#: the corpus assert a census it never had, and a default of `M_measured_quoted` would be a lie by
#: omission. Neither is admissible, so the default is a fifth VALUE on the axis that is not a tier.
NOT_ASSESSED = "not_assessed"

#: Every value the axis can hold: the four tiers plus the absence of one. `len(TIERS)` is the number
#: of tiers; `len(AXIS_VALUES)` is one more, and the difference is load-bearing.
AXIS_VALUES: tuple[str, ...] = (*TIERS, NOT_ASSESSED)

#: Two precedents in this codebase for marking the absence of a reading rather than defaulting it,
#: both cited because the decision above is theirs and not new. Neither is changed by this proposal.
PRECEDENTS_FOR_AN_ABSENT_READING = {
    "UnstatedConfidence": (
        "`genomeos/ir/model.py`: a confidence the program did not state equals 0.0 in every "
        "calculation but is not a stated 0.0, because `confidence: 0.0` also exists in the corpus "
        "and the number alone cannot tell 'stated as low' from 'not stated'. `not_assessed` is the "
        "same distinction one level up: 'nobody looked' against 'someone looked and found nothing'."
    ),
    "context_evidence": (
        "`genomeos/ir/model.py` Rule.context_evidence: the empty string where no chromatin reading "
        "was taken at all, which the field's own docstring says 'is not the same as "
        "`not_assessable`, a reading that was attempted and could not be taken'. It also states "
        "that it 'records a measurement beside the rule and changes nothing the rule does: "
        "`applies` does not read it' - which is precisely the runtime behaviour decision 3 takes."
    ),
}

# --- 2. THE FOUR DECISIONS, FIXED HERE BEFORE ANYTHING IS BUILT --------------------------------------

DECISIONS: dict[str, str] = {
    "1_what_the_tiers_are": (
        "The census's four classes, imported as the same objects, with NO fifth class and no "
        "renaming, plus `not_assessed` which is the absence of a tier rather than a tier. One "
        "mismatch is declared rather than papered over: the census defines `M_measured_quoted` as "
        "'a source that was FETCHED AND READ IN THIS LANE', which is indexed to an act of fetching "
        "by one lane and cannot be asserted by a program. A tier written in a program therefore "
        "asserts the class AS ADJUDICATED BY A NAMED ADJUDICATION, and must carry a reference to "
        "it. That is a second FIELD beside the tier, not a second axis, so the axis stays identical."
    ),
    "2_where_the_tier_lives": (
        "BOTH, on the two axes the census already keeps apart, because a single tier per rule would "
        "discard information the census has already measured: it found tiers DIFFERING WITHIN one "
        "rule across `strength`, `threshold` and `hill` (read `4d4003d`). So: the four-tier axis "
        "lives PER SLOT - one mark per number, keyed by field name, three per rule. The rule itself "
        "carries a mark too, but on the census's SECOND axis (`INTERACTION_STANDING`), because the "
        "rule-level claim is the EXISTENCE of the interaction and classes defined for the value of "
        "a quantity cannot class an existence claim. A per-rule FLOOR over the slot tiers is "
        "computed and never written, so it cannot disagree with the slots it summarises. The two "
        "levels are not independent - the census's own definition of `U_unsourced` says the slot's "
        "citation 'does not even establish the interaction' - so storing both requires a coherence "
        "check between them, which is built and tested rather than assumed."
    ),
    "3_what_the_runtime_does": (
        "NOTHING. The tier gates nothing: no rule is refused, skipped, down-weighted or reordered "
        "because of it, and no integration result moves by a digit when tiers are present. The "
        "instinct is registered as a PREDICTION to be tested and not as an assumption - see P1, P4 "
        "and P5. What is added is downstream of the integration and not inside it: a result "
        "computed from a program can STATE THE LOWEST TIER IT RESTS ON, over the rules that "
        "actually fired, so a figure cannot travel without it. The aggregation is a FLOOR and "
        "explicitly not a mean, and that choice is made against an existing instrument rather than "
        "in the abstract: `genomeos/runtime/uncertainty.py` already attaches a report to every "
        "simulation output and aggregates by running MEAN over evidence kind and confidence. P2 and "
        "P3 test whether that existing instrument already answers this question."
    ),
    "4_the_default_and_the_cost": (
        "The default is `NOT_ASSESSED`, for the reason recorded on that constant: `M_measured_quoted` "
        "would be a lie by omission, and `U_unsourced` would be a false certificate - it would have "
        "every unsearched slot in the corpus assert that a search was made and came back empty. "
        "Adopting this therefore marks almost the whole corpus `not_assessed` on day one, including "
        "`data/demo/gastrulation.bio` itself until the census at `4d4003d` is transcribed into it. "
        "That is the honest state rather than a defect, and the proposal must SAY THE NUMBER: the "
        "adoption cost is counted by `scripts/rule_evidence_tier_cost.py` over the corpus and is "
        "not asserted here."
    ),
}

# --- 3. WHERE THE MARK IS STORED -------------------------------------------------------------------

STORAGE_OPTIONS: dict[str, str] = {
    "A_in_the_program_text": (
        "a new `rule` block property, so the mark is readable FROM THE PROGRAM ALONE without "
        "consulting any census result. Costs a grammar entry, a parser branch, an IR field and a "
        "BioIR JSON round trip, and every hand-authored rule must be re-authored to say anything "
        "other than `not_assessed`."
    ),
    "B_a_side_car_adjudication": (
        "a record keyed by (program sha256, rule id, field), which is the shape the census result "
        "already has. Costs no grammar change and touches no program, and a program edit "
        "invalidates the adjudication automatically because the sha256 moves."
    ),
}

STORAGE_DECISION = (
    "A, in the program text, and the reason is a requirement and not a preference: the mark must "
    "distinguish `not_assessed` from `U_unsourced` FROM THE PROGRAM ALONE, and a side-car record is "
    "by construction not the program alone. B cannot meet that requirement, so it is recorded as "
    "the rejected option with its advantages intact rather than dropped."
)

STORAGE_LIMIT_OF_THE_CHOSEN_OPTION = (
    "A tier is not a property of the program. It is a property of THE LITERATURE AT A TIME, "
    "adjudicated by someone. Written into the program text it can go stale while the file does not "
    "change - a measurement published tomorrow does not move a byte of the `.bio` file, and no "
    "registration travels with a hand-edited line. This is the cost of meeting the "
    "read-from-the-program-alone requirement and it is not solved by this proposal. It is only "
    "MITIGATED, by requiring each mark to name its adjudication (decision 1), so that a reader can "
    "see which adjudication a mark came from and ask whether it is current, instead of reading a "
    "bare class with no date and no author."
)

# --- 4. PREDICTIONS, REGISTERED BEFORE THEY ARE TESTED ---------------------------------------------
#
# Each is falsifiable from the code and will be tested rather than asserted. A prediction that comes
# out false is reported as false and the decision it supported is revisited in the open.

PREDICTIONS: tuple[dict[str, str], ...] = (
    {
        "id": "P1",
        "claim": (
            "A gate at `M_measured_quoted` leaves NO rule of `data/demo/gastrulation.bio` "
            "integrable: every one of its rules is refused, so the program cannot be run at all. "
            "This is why the tier gates nothing; the alternative is not a stricter corpus, it is an "
            "unrunnable one."
        ),
        "how": (
            "construct the gate over the census's own per-slot classes and count the surviving "
            "rules. No `.bio` file is edited and no run of the program is required."
        ),
    },
    {
        "id": "P2",
        "claim": (
            "The existing instrument does NOT already answer this question: "
            "`uncertainty.report_for_network` reports gastrulation's molecular level with a label "
            "other than the lowest one available, while the census classes none of that program's "
            "rule numbers as measured. If so, the existing report is not merely silent on the "
            "tier - a reader can mistake its label for grounding OF THE NUMBERS, when what it "
            "aggregates is the declaration's evidence kind and confidence, which are about the "
            "INTERACTION. THIS LANE HAS NOT COMPUTED THAT LABEL. The prediction is registered "
            "before it is computed and the computed value is reported whichever way it falls."
        ),
        "how": "call the existing reporter on the parsed program and read its label and mean.",
    },
    {
        "id": "P3",
        "claim": (
            "`UncertaintyReport` cannot express a floor, for two separate reasons, and both are "
            "properties of the committed code rather than opinions about it. (a) Its aggregation is "
            "a running mean, so one unsourced item is diluted by sourced neighbours. (b) Its "
            "`weakest` field is NOT the minimum: the guard is `score < lr.confidence and (not "
            "lr.weakest or score < self._score_of(lr))` and `_score_of` returns `lr.confidence`, "
            "the running MEAN, so `weakest` holds the last item that fell below a moving average. "
            "Predicted: an input can be constructed on which `weakest` names an item that is not "
            "the lowest-scoring one."
        ),
        "how": (
            "feed a constructed sequence of scores to `UncertaintyReport.add` and compare `weakest` "
            "with the true argmin. Nothing in `uncertainty.py` is edited: if the prediction holds, "
            "it is reported as a defect of an existing module for its owner to decide on, and this "
            "proposal's floor is computed in its own code."
        ),
    },
    {
        "id": "P4",
        "claim": (
            "No module under `genomeos/runtime/` imports this module, directly or transitively, so "
            "the tier CANNOT gate an integration however it is later filled in. An import-closure "
            "check is a stronger statement than a passing simulation, because it holds for inputs "
            "nobody ran."
        ),
        "how": "walk the import closure of the runtime package and assert this module is absent.",
    },
    {
        "id": "P5",
        "claim": (
            "The floor over a set of slot marks is `not_assessed` - UNDETERMINED - whenever any "
            "slot in the set is unmarked, and is never silently the lowest tier present. A result "
            "resting on one unassessed rule does not get to report `U_unsourced` as its floor, "
            "because that would claim a search that was not made."
        ),
        "how": "exhaustive over the five axis values for small sets.",
    },
    {
        "id": "P6",
        "claim": (
            "The coherence rule between the two levels is satisfied by the census's own 21 "
            "adjudications, so the rule is validated against a real adjudication rather than "
            "asserted. Specifically: a slot classed `E_existence_only` requires a fetched source "
            "supporting the interaction, so it may not sit on a rule whose interaction standing is "
            "`no_source_cited` or `names_neither_factor`; and a slot classed `U_unsourced` may not "
            "sit on a rule whose standing is `supports_declaration`."
        ),
        "how": (
            "run the coherence check over `rule_number_sources_findings.FINDINGS` and its "
            "`ATTRIBUTION_FINDINGS` and require zero problems, plus a counterfactual that plants an "
            "incoherent pair and requires the check to fire."
        ),
    },
)

# --- 5. WHAT THIS LANE HAD ALREADY SEEN, DECLARED BEFORE IT COUNTS ANYTHING -------------------------

NOT_BLIND = (
    "THE ADOPTION COUNT IS NOT BLIND AND SAYING SO IS PART OF THE REGISTRATION. Before writing this "
    "registration this lane ran `grep -c` for rule-declaration lines over the `.bio` corpus by "
    "directory and saw approximate per-directory figures, and it had read the header comments of "
    "`data/organisms/human/*.bio`. So the committed cost script re-derives the counts, and its "
    "figures are a re-derivation of numbers this lane has already glanced at, not a first look. The "
    "line-scan it replaces is also named as a weaker instrument than a parse, which is why the cost "
    "script parses the hand-authored programs rather than scanning them."
)

EXCLUSIONS: tuple[str, ...] = (
    "No `.bio` file is edited, created or written to by this proposal or by any of its code or "
    "tests - not even to demonstrate the mark. Fixtures and temporary copies only.",
    "No number is changed anywhere: no strength, no threshold, no Hill coefficient, no confidence, "
    "no parameter. The proposal has no path that writes a value into a program.",
    "No recommendation is made about whether to adopt. The cost and the limits are produced; the "
    "decision is Albert's.",
    "The grammar, the parser, the IR and the BioIR JSON representation are NOT changed. Storage "
    "option A is costed, not implemented, because implementing it is the adoption this proposal "
    "does not pre-empt.",
    "`genomeos/runtime/uncertainty.py` is not edited. P3 may find a defect in it; a defect found is "
    "reported to its owner, not fixed here under cover of a proposal.",
    "The census's figures are cited to `4d4003d` and are not restated as this lane's findings.",
    "`data/demo/gastrulation.bio` line 3 is reported and not edited; the census already recorded it "
    "at `4d4003d` with `edited: false` and it is still unedited in the program.",
    "0 model requests. No money. No network.",
)

REPORTS: tuple[str, ...] = (
    "the axis, its five values, and which four of them are tiers",
    "where the mark lives at each level, and the coherence rule between the levels",
    "what the runtime does, and which predictions were tested rather than assumed",
    "the adoption cost in programs and in rules, split hand-authored against machine-emitted",
    "the default, and why it is neither `M_measured_quoted` nor `U_unsourced`",
    "every prediction that came out false, with the decision it had supported",
)

WHICHEVER_WAY_IT_FALLS = (
    "If P2 comes out false - if the existing uncertainty report already labels gastrulation at its "
    "lowest - then the gap this proposal fills is smaller than the brief states, and that is "
    "reported in those words rather than worked around. If P3 comes out false, the floor could be "
    "computed by the existing module and this proposal's own floor code is redundant; that is "
    "reported too. If P6 comes out false, the two-level design is unsound as registered and the "
    "coherence rule is withdrawn rather than weakened until it passes."
)


def registration() -> dict[str, Any]:
    """The registered proposal, with nothing filled in: no count, no verdict, no corpus figure."""
    return {
        "proposal": "rule_evidence_tier",
        "lane": "lane-evtier",
        "question": (
            "Should a BioLang program mark, per rule number, what that number rests on - a "
            "measurement with a quoted source, a fitted or modelled value, a citation that supports "
            "only the interaction, nothing, or no assessment at all - so that dynamics computed "
            "from the program cannot be read as mechanism?"
        ),
        "why_cited_not_restated": (
            "data/results/rule_number_sources_census.json, registration 4392d24, result 4d4003d"
        ),
        "axis_values": list(AXIS_VALUES),
        "tiers": list(TIERS),
        "tier_definitions": dict(TIER_DEFINITIONS),
        "tiers_are_the_census_objects": (
            "TIERS is rule_number_sources.CASCADE and TIER_DEFINITIONS is "
            "rule_number_sources.CLASSES, the same objects, asserted by `is` in the suite"
        ),
        "not_assessed": NOT_ASSESSED,
        "not_assessed_is_the_default_and_is_not_a_tier": True,
        "tiered_fields": list(TIERED_FIELDS),
        "interaction_standing_values": dict(INTERACTION_STANDING),
        "precedents_for_an_absent_reading": dict(PRECEDENTS_FOR_AN_ABSENT_READING),
        "decisions": dict(DECISIONS),
        "storage_options": dict(STORAGE_OPTIONS),
        "storage_decision": STORAGE_DECISION,
        "storage_limit_of_the_chosen_option": STORAGE_LIMIT_OF_THE_CHOSEN_OPTION,
        "predictions": [dict(p) for p in PREDICTIONS],
        "not_blind": NOT_BLIND,
        "exclusions": list(EXCLUSIONS),
        "reports": list(REPORTS),
        "whichever_way_it_falls": WHICHEVER_WAY_IT_FALLS,
        # Registered empty, and these four counts prove the emptiness rather than claim it.
        "adoption_cost": {},
        "adoption_cost_entries": 0,
        "prediction_outcomes": {},
        "prediction_outcomes_entries": 0,
        "recommendation": "",
        "recommendation_length": 0,
    }


# ===================================================================================================
# IMPLEMENTATION. Everything above this line was committed at `0d4247b` and written into
# `data/results/rule_evidence_tier_proposal_registration.json` at `709aebc` before any of the
# following existed. Nothing below changes a decision above; where a prediction came out false it is
# recorded in `scripts/rule_evidence_tier_cost.py`'s output and reported, not quietly amended.
# ===================================================================================================

from collections.abc import Iterable  # noqa: E402
from dataclasses import dataclass  # noqa: E402


def is_tier(value: str) -> bool:
    """True for the four tiers, False for `not_assessed` and for anything off the axis."""
    return value in TIERS


def rank(value: str) -> int:
    """Position on the census's cascade: 0 is the most-grounded tier, `len(TIERS) - 1` the least.

    `not_assessed` HAS NO RANK and raises, deliberately. It is not a worse tier than `U_unsourced`;
    it is not a tier, and giving it a rank is exactly how it would start being compared with one.
    """
    if value == NOT_ASSESSED:
        raise ValueError(
            f"{NOT_ASSESSED!r} has no rank: it is the absence of a tier, not the lowest one. "
            "A caller that needs to order it is comparing a finding with the lack of one."
        )
    if value not in TIERS:
        raise ValueError(f"{value!r} is not on the axis; the axis is {list(AXIS_VALUES)}")
    return TIERS.index(value)


def floor(values: Iterable[str]) -> str:
    """The lowest tier a set of marks rests on, or `NOT_ASSESSED` when that is undetermined.

    Undetermined in two cases, and both return `NOT_ASSESSED` rather than a tier:

    * any mark in the set is `NOT_ASSESSED` - one unsearched slot means the set's floor is unknown,
      and reporting `U_unsourced` instead would claim a search nobody made;
    * the set is EMPTY - nothing was marked, so nothing is known. An empty set returning the most
      grounded tier is the vacuous-truth bug that would let an unmarked result read as measured.
    """
    vals = list(values)
    off_axis = sorted({v for v in vals if v not in AXIS_VALUES})
    if off_axis:
        raise ValueError(f"not on the axis: {off_axis}; the axis is {list(AXIS_VALUES)}")
    if not vals or NOT_ASSESSED in vals:
        return NOT_ASSESSED
    return max(vals, key=rank)


@dataclass(frozen=True, slots=True)
class SlotMark:
    """One mark on one rule number. `adjudication` names WHO decided and under what registration.

    The tier alone is not self-supporting: the census defines `M_measured_quoted` by an act of
    fetching, so a mark asserts a class as adjudicated by a named adjudication (decision 1). An
    empty `adjudication` on a mark that is not `NOT_ASSESSED` is a problem the coherence check
    reports.
    """

    rule: str
    field: str
    tier: str = NOT_ASSESSED
    adjudication: str = ""


@dataclass(frozen=True, slots=True)
class RuleMark:
    """The rule-level mark: what the rule's own citation does for the interaction it states.

    On the census's SECOND axis, not the tier axis, because the claim is the EXISTENCE of the
    interaction and the four tiers are defined for the value of a quantity.
    """

    rule: str
    interaction_standing: str = ""
    adjudication: str = ""


#: The coherence rules between the two levels, as registered in prediction P6. Each is derived from
#: the census's OWN class definitions and not from the data it produced; both are validated against
#: the census's 21 adjudications and each has a counterfactual in the suite.
COHERENCE_RULES = {
    "E_needs_a_source_for_the_interaction": (
        "`E_existence_only` is defined as a fetched source supporting the EXISTENCE of the "
        "interaction while stating no value. A slot may therefore not carry it on a rule whose "
        "interaction standing is `no_source_cited` or `names_neither_factor`: there would be no "
        "source for the existence for it to rest on."
    ),
    "U_needs_the_citation_not_to_establish_the_interaction": (
        "`U_unsourced` is defined to include that the slot's own citation 'does not even establish "
        "the interaction'. A slot may therefore not carry it on a rule whose standing is "
        "`supports_declaration`; such a slot is `E_existence_only`."
    ),
    "a_tier_names_its_adjudication": (
        "A mark that is not `not_assessed` asserts a class that was decided by someone, so it must "
        "name the adjudication that decided it. A bare tier with no adjudication is the kind of "
        "claim no artefact supports."
    ),
}


def coherence_problems(
    slot_marks: Iterable[SlotMark], rule_marks: Iterable[RuleMark]
) -> list[dict[str, str]]:
    """Every place the per-slot tier and the per-rule standing cannot both be true. Empty is good."""
    standing = {m.rule: m.interaction_standing for m in rule_marks}
    problems: list[dict[str, str]] = []
    for m in slot_marks:
        if m.tier not in AXIS_VALUES:
            problems.append({"rule": m.rule, "field": m.field, "rule_broken": "off_the_axis", "tier": m.tier})
            continue
        st = standing.get(m.rule, "")
        if m.tier == "E_existence_only" and st in ("no_source_cited", "names_neither_factor"):
            problems.append(
                {
                    "rule": m.rule,
                    "field": m.field,
                    "rule_broken": "E_needs_a_source_for_the_interaction",
                    "tier": m.tier,
                    "interaction_standing": st,
                }
            )
        if m.tier == "U_unsourced" and st == "supports_declaration":
            problems.append(
                {
                    "rule": m.rule,
                    "field": m.field,
                    "rule_broken": "U_needs_the_citation_not_to_establish_the_interaction",
                    "tier": m.tier,
                    "interaction_standing": st,
                }
            )
        if is_tier(m.tier) and not m.adjudication:
            problems.append(
                {
                    "rule": m.rule,
                    "field": m.field,
                    "rule_broken": "a_tier_names_its_adjudication",
                    "tier": m.tier,
                    "interaction_standing": st,
                }
            )
    return problems


def result_floor(active_rule_ids: Iterable[str], slot_marks: Iterable[SlotMark]) -> dict[str, Any]:
    """What a result computed from a program states: the lowest tier it rests on, over what fired.

    `active_rule_ids` are the rules the computation actually used, so a rule gated out by `when`
    does not drag a figure down and a rule that fired cannot be left out of the floor. Every
    `TIERED_FIELDS` slot of every active rule is counted, and a slot with no mark counts as
    `NOT_ASSESSED` - so a result over an unmarked program is UNDETERMINED, which is the state the
    whole corpus is in on the day this is adopted.
    """
    ids = list(active_rule_ids)
    marks = {(m.rule, m.field): m.tier for m in slot_marks}
    per_slot = {(rid, f): marks.get((rid, f), NOT_ASSESSED) for rid in ids for f in TIERED_FIELDS}
    counts = {v: 0 for v in AXIS_VALUES}
    for v in per_slot.values():
        counts[v] += 1
    fl = floor(per_slot.values())
    determined = fl != NOT_ASSESSED
    return {
        "rules_that_fired": len(ids),
        "slots": len(per_slot),
        "counts": counts,
        "lowest_tier_it_rests_on": fl,
        "determined": determined,
        "statement": (
            f"this figure rests on {counts[fl]} of {len(per_slot)} rule numbers at {fl}"
            if determined
            else (
                f"the lowest tier this figure rests on is UNDETERMINED: {counts[NOT_ASSESSED]} of "
                f"{len(per_slot)} rule numbers behind it have never been assessed"
            )
        ),
    }


def gate_survivors(active_rule_ids: Iterable[str], slot_marks: Iterable[SlotMark], at: str) -> list[str]:
    """Which rules a hypothetical gate at tier `at` would leave integrable. Hypothetical: NOTHING in
    this repository calls it on a real run, and the runtime does not import this module at all. It
    exists so the cost of gating can be COUNTED instead of argued about (prediction P1)."""
    keep = {t for t in TIERS if rank(t) <= rank(at)}
    marks = {(m.rule, m.field): m.tier for m in slot_marks}
    return [
        rid
        for rid in active_rule_ids
        if all(marks.get((rid, f), NOT_ASSESSED) in keep for f in TIERED_FIELDS)
    ]


def census_slot_marks(adjudication: str = "4d4003d") -> tuple[list[SlotMark], list[RuleMark]]:
    """The census's own 21 adjudications, read out of its findings module, as marks.

    This is the only adjudication that exists. It is READ, never written back into any program, and
    it is what the coherence rule is validated against: a rule derived from the census's class
    definitions must hold on the census's own output or it is wrong (prediction P6).
    """
    from genomeos.provenance import rule_number_sources_findings as rnf

    slots = [
        SlotMark(rule=rule, field=field, tier=rec["class"], adjudication=adjudication)
        for (rule, field), rec in rnf.findings().items()
    ]
    rules = [
        RuleMark(rule=f["rule"], interaction_standing=f["attribution"], adjudication=adjudication)
        for f in rnf.attribution_findings()
    ]
    return slots, rules


# --- NARROWED AFTER REGISTRATION, and named as narrowed with the reason ----------------------------
#
# The registration says a result states the lowest tier over "the rules that actually fired". The
# code contradicts the obvious reading of that phrase, so it is narrowed here rather than left to be
# read generously later.
#
# `genomeos/lang/parser.py` SYNTHESISES a rule from every `gene ... { produces: ... }` clause:
# `data/demo/gastrulation.bio` writes 7 rules and compiles to 11, and the four extra carry
# `strength` 1.0, `threshold` 1.0 and `hill` 2.0 from the `genomeos/ir/model.py` dataclass defaults -
# twelve numbers no line of any program states. All four are in `Module.active_rules` and all four
# therefore "fire". But `genomeos/runtime/grn.py` NEVER READS their numbers: for a PRODUCE rule it
# keeps `r.source` alone (`self.produces[r.target].append(r.source)`), and translation is driven by
# `translation_rate`, not by the rule's Hill function.
#
# Taking "the rules that fired" literally would therefore leave every floor in the corpus permanently
# UNDETERMINED on twelve numbers the dynamics never touch - a false alarm rather than a false
# certificate, but false either way, and a floor that is always undetermined says nothing. The
# population is narrowed to the rules whose tiered fields are read. A test perturbs the synthesised
# rules' three numbers and requires the trajectory not to move, so this rests on the runtime's
# behaviour rather than on a reading of its source.
#
# The narrowing CUTS BOTH WAYS and is not a convenience: the hypothetical gate of P1 is applied over
# the same narrowed population, so it does not get to pass four rules it never judged.

#: The actions whose `strength`, `threshold` and `hill` the GRN runtime reads. Named by string so
#: this module does not import the runtime - P4 requires the closure to stay one-way.
TIERED_ACTIONS: tuple[str, ...] = ("ACTIVATE", "INHIBIT")

TIERED_ACTIONS_WHY = (
    "`genomeos/runtime/grn.py` reads a rule's strength, threshold and hill only for ACTIVATE and "
    "INHIBIT. A PRODUCE rule contributes its source and nothing else, so its three numbers are "
    "inert and giving them a tier would mark a number nothing reads. MODIFY, BIND and DEGRADE are "
    "excluded for the same reason and not because they matter less: no path in the GRN runtime "
    "reads their tiered fields either. If a runtime later reads them, this tuple is what has to "
    "change, and it is one place."
)


def tiered_rules(rules: Iterable[Any]) -> list[Any]:
    """The rules whose tiered numbers a GRN result rests on. Takes IR rules; duck-typed on `action`
    so this module imports neither the runtime nor the IR."""
    return [r for r in rules if getattr(r.action, "name", str(r.action)) in TIERED_ACTIONS]
