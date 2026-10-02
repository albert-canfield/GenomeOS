# SPDX-License-Identifier: AGPL-3.0-or-later
"""The registration of a second measured-layer extractor, fixed and committed BEFORE the extractor
changes (`data/results/increase_links_registration.json`).

**The finding this acts on, established by lane-increase and not re-established here.**
`measured.rule_links` emits a link only from a pair the benchmark calls `Regulated`, which it awards a
significant *decrease* only, so every CRISPRi link's action is `activates` and a significant *increase*
raises no rule at all. 143 of 195 assessable measured-layer rules sit on a pure `activates` axis and 0
on a pure `represses` axis, which is why an earlier lane's repression gate returned 0 eligible links.
The `inhibits` branch of `rule_links` already exists and is unreachable only because the candidate set
is gated on `regulated`: this lane makes a written-but-dead path reachable and adds no new direction
rule. Verified against the two cached benchmark tables before this registration was written: of 14,734
valid pairs, 661 are `Regulated TRUE` and **every one of them has an effect size below zero**, so the
`inhibits` branch is reached by nothing; 820 pairs are significant, 661 decreases and 159 increases.

**What this module is.** The eligibility predicate; what an increase-derived link's action, evidence
status and outcome are, and what its molecular role is not; the conflict rule; both floors, imported
from where the project fixed them and chosen by nobody here; and the by-construction property the new
links share with the old ones. The extractor that applies all of it is `measured.rule_links` under
`EXTRACTOR_V2`, changed only after this registration is committed. Nothing here moves afterwards for
any reason.

**The wording that is binding (review item R2), and it is binding in code and not only in prose.** An
increase-derived link records the MEASURED NET DIRECTION and nothing more. A measured increase on
knockdown is *an increase on knockdown*. It is never a silencer, never a repressor and never evidence
of a repression mechanism, because an increase under CRISPRi can be indirect - KRAB spread onto a
neighbouring promoter, competition between promoters for a shared enhancer - and a rule that said
otherwise would state a mechanism no measurement here distinguishes. `FORBIDDEN_OF_AN_INCREASE_LINK`
names the words that may not appear on such a link, and `check_no_mechanism_claim` refuses them, so the
prohibition is enforced by a function and not by a reader's care.

**What is not done, and this is a condition of the approval.** No direction is read, no repression call
is tested, no rate is reported and no interval is taken: building the population is this lane, testing
it is a later registration with its own floors. `rule_links` v1 stays the default and byte-identical,
so no pinned count moves and no downstream result shifts under anyone; v2 is selected explicitly. No
model request is made, nothing is downloaded, no money is spent.
"""

from __future__ import annotations

from typing import Any

from genomeos.attribution import cell2, fresh
from genomeos.attribution import direction_link as dl
from genomeos.attribution import increases as inc
from genomeos.attribution import repress2 as rp

REGISTERED = "2026-10-02"

# ---- what is inherited, imported rather than restated so the two cannot drift ---------------------

#: lane-increase's registered GO reading, carried word for word and never strengthened.
INHERITED_INCREASE_GO = inc.GO
#: lane-repress2's registered no-go, carried through lane-increase's own constant.
INHERITED_REPRESS2_NO_GO = rp.GATE_NO_GO
#: the structural cause, imported from where it was first recorded.
STRUCTURAL_CAUSE = rp.MEASURED_AXIS_IS_NOT_THE_LINK_DIRECTION
#: where that cause is in the code, named by file and expression, as lane-increase recorded it.
CAUSE_IN_CODE = inc.CAUSE_IN_CODE

INHERITED_FILES = (
    "data/results/increase_population.json",
    "data/results/increase_registration.json",
    "data/results/repress2_population.json",
)
INHERITED_COMMITS = {
    "increase_registration, committed before any count": "e1c2671",
    "genomeos/attribution/increases.py": "6ca8095",
}
#: Every figure of lane-increase's that bounds this lane. None is recomputed here, and none may be
#: quoted under a noun its own breakdown does not support.
INHERITED_FIGURES = {
    "valid_benchmark_pairs": 14_734,
    "significant_pairs": 820,
    "significant_decreases": 661,
    "significant_increases": 159,
    "increases_on_an_attributed_element": 48,
    "increase_links_over_independent_loci": 33,
    "decreases_on_an_attributed_element": 212,
    "increase_links_contested_by_a_regulated_pair": 0,
    "assessable_measured_layer_rules": 195,
    "on_a_pure_activates_axis": 143,
    "on_a_pure_represses_axis": 0,
}

# ---- the versioned extractor ---------------------------------------------------------------------

#: The version names, following the shape the judge's target rule set in
#: `genomeos.attribution.correctness` (RULE_V1/RULE_V2/RULE_V3, RULES, DEFAULT_RULE): the earlier
#: version keeps the values it registered, the new one is selected by name, and a result names the
#: version it was made under. `measured` owns the two strings the code branches on; this module
#: registers what they MEAN, and `test_increase_links` ties the two together so they cannot drift.
EXTRACTOR_V1, EXTRACTOR_V2 = "v1", "v2"
EXTRACTORS = (EXTRACTOR_V1, EXTRACTOR_V2)
#: v1 stays the default, which is the whole reason no pinned count moves. This differs on purpose from
#: the judge, where the new rule became the default: there the old reading was wrong, here the old
#: reading is right and merely incomplete, and every committed caller of `rule_links` is entitled to it.
DEFAULT_EXTRACTOR = EXTRACTOR_V1

EXTRACTOR_RULE = {
    EXTRACTOR_V1: (
        "v1 (the committed extractor, unchanged and byte-identical): one link per (gene, cell) carrying "
        "a pair the benchmark calls `Regulated`, which is a significant decrease only. Its action is "
        "`activates`, its strength is the largest |EffectSize| among the regulated pairs of that (gene, "
        "cell) preferring training pairs, and its split is `training` if any regulated training pair "
        "exists and `heldout` otherwise. A significant increase raises nothing. This is the version a "
        "call that names none uses, and every count pinned on it therefore stands"
    ),
    EXTRACTOR_V2: (
        "v2 (from 2026-10-02, selected explicitly and never the default): v1's links, plus one link per "
        "(gene, cell) carrying a pair whose outcome is a significant increase, with action `inhibits`, "
        "strength and split taken by v1's own rules applied to the increase pairs of that (gene, cell), "
        "MINUS every (gene, cell) the conflict rule refuses. v2 adds a population and changes no rule "
        "v1 states: on a row with no conflict its output is v1's output followed by the new links, "
        "identical to v1's in count and content on the shared part. v2 reaches the `inhibits` branch "
        "that has been written in `rule_links` all along and introduces no new direction rule"
    ),
}

#: Downstream names, fixed here so that no result made under v2 can ever be mistaken for one made under
#: v1. This lane writes none of them: it builds the population and stops.
V2_CARRIES_NEW_NAMES = (
    "a compiled program or downstream result made under v2 carries a NEW NAME and never overwrites a "
    "v1 one: the suffix `_v2` on the result name, `extractor: v2` in the manifest parameters, and "
    "`extractor_v2` on the compiled evidence source. No such program or result is built by this lane, "
    "so the names are registered before anything can claim them rather than after"
)

# ---- the eligibility predicate -------------------------------------------------------------------

#: The one label the new links are built from, read off the cached row's own `outcome` field and never
#: recomputed. Imported from lane-increase so the two lanes cannot come to mean different things by it.
INCREASE = inc.INCREASE

ELIGIBILITY_CALL = (
    "a (gene, cell) on a measured row is an increase-derived link candidate when at least one of that "
    f"row's cached CRISPRi pairs for that gene and that cell has `outcome` exactly {INCREASE!r}, read "
    "off the cached row's own `outcome` field and never recomputed, exactly as lane-increase's "
    "predicate does for a pair. That is the whole eligibility predicate and it is fixed here before the "
    "extractor changes. The overlap rule, the reach and the cell gating are the committed ones and are "
    "not touched: the candidate set is widened at the sign and nowhere else. A link is one (gene, cell) "
    "as review item R1 requires, so two cells that measured the same gene are two links and never the "
    "stronger of the two"
)

STRENGTH_AND_SPLIT_CALL = (
    "an increase-derived link takes its strength and split by v1's own rules applied to its own pairs: "
    "strength is min(1.0, |EffectSize|) of the increase pair with the largest |EffectSize| among the "
    "training increase pairs of that (gene, cell) if any exist and among all of them otherwise, rounded "
    "to three places; split is `training` if any increase training pair exists and `heldout` otherwise, "
    "so a held-out measurement never sets a compiled number. No new magnitude rule is introduced, and "
    "the strength is a magnitude with no sign in it: the sign lives in the action"
)

# ---- the conflict rule ---------------------------------------------------------------------------

#: Decided by the supervisor on approval, registered now and NEVER resolved by choosing a sign. It is 0
#: of 48 on today's data, which is exactly why it is registered now rather than when it first bites.
CONFLICT_RULE = (
    "a (gene, cell) carrying BOTH a significant decrease and a significant increase emits NO RULE under "
    "v2 - not an `activates` rule, not an `inhibits` rule, nothing. Both observations go into a "
    "separate conflicts list with their element, their effect sizes, their splits and their datasets, "
    "and the list is counted and reported. The conflict is NEVER resolved by choosing a sign, by taking "
    "the larger magnitude, by preferring the training pair, by preferring the earlier dataset or by any "
    "other tie-break: two significant measurements of opposite sign on one gene in one cell are a "
    "disagreement between measurements, and a rule emitted from either would state as observed a "
    "direction the data does not agree on. Suppressing the v1 `activates` link of a conflicted (gene, "
    "cell) is the one place v2 can remove a link v1 emitted; it removes 0 today, the conflicts list "
    "makes any future removal visible and counted rather than silent, and the additivity check in "
    "`test_increase_links` asserts the conflicts list is empty at the same time as it asserts v2 is v1 "
    "plus the new links, so the day the rule first bites is the day that test says so"
)
CONFLICT_IS_NOT_A_FINDING = (
    "a conflict count of 0 is not a finding that the measurements agree about direction: it is a count "
    "of (gene, cell) pairs where the benchmark holds significant observations of both signs, over the "
    "small population that reaches an attributed element at all. It is reported as the count it is"
)

# ---- what an increase-derived link records, and what it does not ---------------------------------

ACTION = "inhibits"
EVIDENCE_STATUS = "observed"
OUTCOME_TEXT = "increase on knockdown"
#: The R7 axes keep the molecular role unresolved on an increase-derived link, and `None` is the value
#: that says so: not a role, not an empty string standing in for one, and never inferred from a sign.
MOLECULAR_ROLE = None
MOLECULAR_ROLE_TEXT = "unresolved"

R2_WORDING = (
    "a measured increase on knockdown is 'an increase on knockdown'. It is never a silencer, never a "
    "repressor and never evidence of a repression mechanism. This test is about sign agreement, not "
    "mechanism"
)
RECORDS_THE_MEASURED_NET_DIRECTION = (
    f"an increase-derived link records the MEASURED NET DIRECTION: action {ACTION!r}, evidence status "
    f"{EVIDENCE_STATUS!r}, outcome {OUTCOME_TEXT!r}, molecular role {MOLECULAR_ROLE_TEXT} (the R7 axes "
    "leave it unresolved and this lane does not resolve it). The reason is substantive and not a "
    "formality: an increase on CRISPRi can be indirect - KRAB spread onto a neighbouring promoter, "
    "competition between promoters for a shared enhancer - so the net direction is what was measured "
    "and the mechanism is not. The link says the gene went up when the element was silenced; it does "
    "not say the element represses the gene, and nothing downstream may read it as though it did"
)
#: The words that may not appear on an increase-derived link, in its role or in any text it carries.
#: `check_no_mechanism_claim` refuses them, and a test plants one to show the refusal is real.
FORBIDDEN_OF_AN_INCREASE_LINK = ("silencer", "repressor", "represses", "repression", "represser")
PERMITTED_OUTCOME_PHRASE = OUTCOME_TEXT

# ---- the floors, imported and not chosen ---------------------------------------------------------

#: Both floors come from where the project fixed them before this lane existed. Neither is chosen here,
#: neither is applied here - this lane builds the population and takes no gate - and neither may move
#: after an outcome is seen. They are recorded so the later test cannot pick its own.
POSITIVE_FLOOR = fresh.POSITIVE_FLOOR
LOCUS_FLOOR = fresh.LOCUS_FLOOR
FLOORS = {
    "link_floor": POSITIVE_FLOOR,
    "link_floor_imported_from": "genomeos.attribution.fresh.POSITIVE_FLOOR",
    "locus_floor": LOCUS_FLOOR,
    "locus_floor_imported_from": "genomeos.attribution.fresh.LOCUS_FLOOR, which is "
    "genomeos.attribution.cell2.POOLED_LOCUS_FLOOR",
    "not_applied_here": (
        "no gate is taken by this lane and no floor is compared against anything. The floors are "
        "recorded because the lane that tests this population must use these and not choose its own"
    ),
}
LOCUS_RULE = cell2.INDEPENDENT_LOCUS_RULE
NOT_BIOLOGICAL_INDEPENDENCE = inc.NOT_BIOLOGICAL_INDEPENDENCE

# ---- the by-construction property, stated by the links themselves -------------------------------

#: The trap, imported in the wording lane-direction registered it in: the decrease arm's 212 of 212 is
#: the extractor's own construction and not model coverage.
DECREASE_ARM_IS_NOT_MODEL_COVERAGE = dl.DECREASE_ARM_IS_NOT_MODEL_COVERAGE

#: Carried on EVERY increase-derived link, so the result states the property of itself rather than
#: leaving a reader to infer it from a count that looks like coverage.
BY_CONSTRUCTION = (
    "this link exists because a significant increase was measured on an element the compiled program "
    "had already attributed, and the link then takes its gene and its cell FROM THAT MEASUREMENT. Its "
    "presence is the extractor's own construction and is NOT evidence that the model covered this gene, "
    "predicted this cell or said anything about this link: a completeness like the decrease arm's 212 of "
    "212 follows from how the extractor is built and never from model coverage. Counting these links as "
    "coverage would repeat the error docs/LESSONS.md records under 'A count is not a measurement of the "
    "thing you want to count'"
)

WHAT_THIS_IS_NOT = (
    "no direction is read, no repression call is tested, no rate is reported and no interval is taken "
    "by this lane. The population is built and the lane stops. A later registration with its own floors "
    "tests it, and the gap it must watch is concrete: the direction test already on the record had 21 "
    "ANSWERABLE links of 48 MEASURED, so a count of what exists is not a count of what can be answered",
)
OUTCOME_BREAKDOWN_RULE = (
    "every count this lane reports carries the exhaustive outcome breakdown at its own denominator. The "
    "reason is on the record twice: 198 loci where the model says nothing were once quoted as 198 places "
    "it is wrong when 66 held a measured decrease, and 48 measured links cleared a floor while only 21 "
    "were answerable. A count is not a measurement of the thing you want"
)


# ---- the prohibition, enforced ------------------------------------------------------------------


def check_no_mechanism_claim(records: list[dict[str, Any]]) -> None:
    """Raise if any increase-derived link claims a mechanism. R2, binding in code and not only in prose.

    An increase-derived link may say the gene went up when the element was silenced. It may not carry a
    molecular role of silencer or repressor, and no text on it may name a repression mechanism. The one
    permitted phrase is `OUTCOME_TEXT`, which describes the observation and not a cause.
    """
    for r in records:
        if r.get("derived_from") != INCREASE:
            continue
        role = r.get("molecular_role")
        if role is not None:
            raise ValueError(
                f"increase-derived link {r.get('gene')}/{r.get('cell')} carries molecular_role {role!r}; "
                f"the R7 axes leave it {MOLECULAR_ROLE_TEXT} and this lane does not resolve it. {R2_WORDING}"
            )
        for key, value in r.items():
            if not isinstance(value, str) or key == "by_construction":
                continue
            text = value.lower()
            for word in FORBIDDEN_OF_AN_INCREASE_LINK:
                if word in text and PERMITTED_OUTCOME_PHRASE not in text:
                    raise ValueError(
                        f"increase-derived link {r.get('gene')}/{r.get('cell')} carries {word!r} in "
                        f"{key!r}: {value!r}. {R2_WORDING}"
                    )


def registration() -> dict[str, Any]:
    """Everything registered here, as the result file states it. Nothing is computed."""
    return {
        "lane": "lane-inhibits",
        "registered": REGISTERED,
        "what_is_registered": (
            "a second version of the measured-layer extractor, which raises a link from a significant "
            "increase so that the project can later test its own repression calls"
        ),
        "acts_on": {
            "lane": "lane-increase",
            "files": list(INHERITED_FILES),
            "commits": dict(INHERITED_COMMITS),
            "registered_go_carried_verbatim": INHERITED_INCREASE_GO,
            "and_behind_it": {
                "lane": "lane-repress2",
                "registered_no_go_carried_verbatim": INHERITED_REPRESS2_NO_GO,
            },
            "figures": dict(INHERITED_FIGURES),
            "structural_cause": STRUCTURAL_CAUSE,
            "cause_in_code": CAUSE_IN_CODE,
        },
        "extractor_is_versioned_not_replaced": {
            "versions": list(EXTRACTORS),
            "default": DEFAULT_EXTRACTOR,
            "rule": dict(EXTRACTOR_RULE),
            "precedent": (
                "the judge's versioned target rule in genomeos.attribution.correctness (v1, v2, v3 with "
                "RULES and DEFAULT_RULE, each version keeping the values it registered), whose shape "
                "this follows. It differs in one respect, on purpose: there the new rule became the "
                "default because the old reading was wrong, here v1 stays the default because the old "
                "reading is right and merely incomplete"
            ),
            "why_v1_stays_the_default": (
                "every pinned count rests on v1 and none of them may move silently. chr21 carries "
                "`# test: rules == 5176` and a BioLang program asserts it, and every result built on the "
                "measured layer would shift under everyone if the default moved. v1 is byte-identical "
                "and selected by every call that names no version, which is all of them"
            ),
            "downstream_names": V2_CARRIES_NEW_NAMES,
        },
        "eligibility_predicate": ELIGIBILITY_CALL,
        "strength_and_split": STRENGTH_AND_SPLIT_CALL,
        "increase_derived_link_records": RECORDS_THE_MEASURED_NET_DIRECTION,
        "binding_wording_r2": R2_WORDING,
        "forbidden_of_an_increase_link": list(FORBIDDEN_OF_AN_INCREASE_LINK),
        "prohibition_is_enforced_by": (
            "genomeos.attribution.increase_links.check_no_mechanism_claim, which raises on a molecular "
            "role or on any text naming a repression mechanism. tests/test_increase_links.py plants a "
            "silencer label and a repressor label and asserts the refusal, so the test fails if the "
            "prohibition is ever weakened"
        ),
        "conflict_rule": CONFLICT_RULE,
        "conflict_is_not_a_finding": CONFLICT_IS_NOT_A_FINDING,
        "no_existing_link_may_change": (
            "v1's output is byte-identical, and on data holding no conflict v2's output is v1's output "
            "plus the new links, equal in count and in content on the shared part. Both are asserted by "
            "test against the genome-wide v1 output captured before the extractor changed. If either "
            "fails the lane stops: the repression links are additive or they are wrong"
        ),
        "floors": dict(FLOORS),
        "independent_locus_rule": LOCUS_RULE,
        "not_biological_independence": NOT_BIOLOGICAL_INDEPENDENCE,
        "by_construction": BY_CONSTRUCTION,
        "decrease_arm_is_not_model_coverage": DECREASE_ARM_IS_NOT_MODEL_COVERAGE,
        "outcome_breakdown_rule": OUTCOME_BREAKDOWN_RULE,
        "what_this_lane_does_not_do": list(WHAT_THIS_IS_NOT),
        "budget": "0 model requests, no money, no downloads, no network",
        "order": [
            "1. this registration, committed before genomeos/attribution/measured.py changes",
            "2. the extractor versioned: v1 byte-identical and the default, v2 added and explicit",
            "3. the proofs: v1 unchanged genome-wide, v2 = v1 + the new links, the conflicts list and "
            "its count, and the R2 refusal shown on a planted silencer and repressor label",
            "4. stop. No direction is read and no repression call is tested here",
        ],
    }
