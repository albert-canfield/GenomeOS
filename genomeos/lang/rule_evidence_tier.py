# SPDX-License-Identifier: AGPL-3.0-or-later
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

from genomeos.lang import rule_number_sources as rns

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
