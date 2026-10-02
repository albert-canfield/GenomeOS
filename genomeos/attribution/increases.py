# SPDX-License-Identifier: AGPL-3.0-or-later
"""Could a measured repression population be built at all: the registration, and the ladder it fixes.

**What this follows from.** lane-repress2 took its gate 1 blind and recorded a no-go in these words:
*"the measured layer holds no repression call to test"* - 0 eligible links and 0 independent loci for
`represses_target` in the rule's own cell, against the same two imported floors this module imports.
Its blind ladder located the zero at the cell: **156,925** repression-axis rules, **357** on an element
a CRISPRi pair covers, **11** with a pair on the rule's own gene, **0** in the cell the rule is gated
on (`data/results/repress2_population.json`, genome-wide). And it recorded the structural cause, which
is this lane's question: `measured.rule_links` emits a link only from a pair the benchmark calls
`Regulated`, which is `Significant AND EffectSize < 0`, so every CRISPRi link's action is `activates`
and a significant *increase* raises no rule at all. It measured the consequence: of **195** assessable
measured-layer rules, **143** are on a pure `activates_target` axis and **0** on a pure
`represses_target` axis (same file).

**The question.** The compiled program makes repression calls that cannot currently be tested, because
the extractor that builds measured links structurally discards the measurements that would test them.
So: how many significant increases exist in the measured data, where are they, and would they form a
population that clears the registered floors? If they would, the project can test repression for the
first time. If they would not, the no-go becomes structural and final rather than an accident of one
arm, and that is worth knowing with the same confidence.

**What this module is.** The registration, fixed before any count: one eligibility predicate, the
locus convention imported from `cell2`, both floors imported from `fresh`, the ladder's steps, the two
populations the floors are applied to, the reading of each outcome written before the outcome was seen,
and the map from the step that kills a population to whether that step is a property of the data or of
this repository's code. The counting code lives in `scripts/increase_population.py` and runs after the
registration is committed.

**What this module is not, and what this lane does not do.** No extractor is changed. `rule_links` and
every other extractor is read and left exactly as committed: this lane counts what a changed extractor
could reach, it does not change one, and if the count says the change is worth making the deliverable
is a proposal and not an edit. No floor is set here, no floor or predicate may move after a count, no
direction of a measured effect enters a ladder step, no rule is deleted, no compiled label changes, no
verdict moves, nothing is refitted, nothing is downloaded and no model request is made.

**Why the ladder carries an outcome breakdown at every step.** docs/LESSONS.md, "A count is not a
measurement of the thing you want to count" (2026-10-02): twice in one night a correct count was quoted
under the wrong noun, once a count of zero-valued columns read as a feature's contribution and once a
count of loci where the model says nothing read as a count of places it is wrong (198 claimed, 66
actual). So every step of this ladder reports the full outcome breakdown at its own denominator, and
the step tests themselves are blind to direction: a pair clears a position step on its gene, its cell
and its element and never on the sign of its effect.
"""

from __future__ import annotations

from typing import Any

from genomeos.attribution import cell2, fresh
from genomeos.attribution import measured as ms
from genomeos.attribution import not_open_profile as nop
from genomeos.attribution import repress2 as rp

# ---- what this lane inherits, carried word for word with its file -------------------------------

#: lane-repress2's registered reading of its own gate 1, carried verbatim. It is not restated in this
#: lane's own words anywhere, and it is never strengthened.
INHERITED_NO_GO = "the measured layer holds no repression call to test"
INHERITED_FROM = "data/results/repress2_population.json"
INHERITED_LADDER = {
    "scope": "genome-wide, all 24 chromosomes",
    "file": INHERITED_FROM,
    "activity_axis": rp.REPRESSES,
    "rules_of_this_axis": 156_925,
    "element_carries_a_crispri_pair": 357,
    "and_a_pair_on_the_rules_own_gene": 11,
    "and_in_the_rules_own_cell": 0,
    "and_that_pair_is_significant": 0,
    "gate_1_reading": rp.GATE_NO_GO,
}
INHERITED_MEASURED_LAYER = {
    "file": INHERITED_FROM,
    "assessable_measured_layer_rules": 195,
    "on_a_pure_activates_target_axis": 143,
    "on_a_pure_represses_target_axis": 0,
    "the_rest": (
        "52 rules whose activity axis holds more than one term: 25 "
        "'activates_target, active_in_reporter', 19 'activates_target, inactive_in_reporter', 4 "
        "'activates_target, represses_target' and 4 "
        "'activates_target, represses_target, inactive_in_reporter'"
    ),
}

#: The structural cause, imported from `repress2` rather than restated so the two cannot drift.
STRUCTURAL_CAUSE = rp.MEASURED_AXIS_IS_NOT_THE_LINK_DIRECTION

#: Where that cause is in the code, named by file and line so the claim can be checked without a
#: search. Read and left as committed; this lane changes neither line.
CAUSE_IN_CODE = (
    "genomeos/attribution/measured.py: `rule_links` builds its candidate set as "
    "`{(p['gene'], p['cell']) for p in c['pairs'] if p['regulated']}`, and `regulated` is set in "
    "`parse_crispri` from `r.get('Regulated') == 'TRUE'`, which the benchmark awards a significant "
    "decrease only (`DECREASE`/`INCREASE` on `CrispriPair.outcome`). The `action` line below it reads "
    "`'activates' if strongest['effect_size'] < 0 else 'inhibits'`, so the `inhibits` branch exists "
    "but is unreachable: every pair that reaches it is regulated and therefore has a negative effect "
    "size. The increase is not lost from the program - `Layer.for_element` records it as "
    "`genes_increased` on the measured block, kept apart from `genes_not_regulated` since review item "
    "R2 - it is only never raised to a rule"
)

#: This lane does not touch either line, and says so where a reader of the result will look.
NO_EXTRACTOR_CHANGED = (
    "no extractor is changed by this lane. `measured.rule_links`, `measured.parse_crispri`, "
    "`measured.Layer.for_element` and `measured.rows` are imported and called exactly as committed, "
    "and the overlap rule, the eligibility rule and the compiled rule enumeration are the committed "
    "ones. The count below is of what a changed extractor could reach; if it says such a change is "
    "worth making, the deliverable is a written proposal and not an edit"
)

# ---- the eligibility predicate -------------------------------------------------------------------

#: The one label this lane's population is built from, read off the cached row's own `outcome` field
#: and never recomputed. `measured.CrispriPair.outcome` gives it to a pair that is `Significant` with a
#: non-negative effect size, which is the pair the benchmark calls `Regulated` FALSE and which
#: `rule_links` therefore discards.
INCREASE = ms.INCREASE
#: Both labels a significant pair takes. The ladder's denominators are taken over these, so no step
#: test reads a direction; the split between them is reported as a breakdown, never as a filter.
SIGNIFICANT_OUTCOMES = rp.SIGNIFICANT_OUTCOMES
#: Every label a pair can take, so a breakdown at a denominator is exhaustive and sums to it.
ALL_OUTCOMES = ms.OUTCOMES

ELIGIBILITY_CALL = (
    "a cached CRISPRi pair is an increase when measured.CrispriPair.outcome is exactly "
    f"{INCREASE!r}, read off the cached row's own outcome field and never recomputed. That is the "
    "whole eligibility predicate and it is fixed here before any count. The ladder's position steps - "
    "the element, the gene and the cell - are applied to every significant pair without reading which "
    f"of {SIGNIFICANT_OUTCOMES} it is, and the direction enters only as the breakdown reported at each "
    "step's own denominator"
)

BREAKDOWN_CALL = (
    "every step of the ladder reports the count of each of "
    f"{list(ALL_OUTCOMES)} at its own denominator, so no count of this lane's can be quoted under a "
    "noun the breakdown does not support. This is docs/LESSONS.md, 'A count is not a measurement of "
    "the thing you want to count' (2026-10-02), applied in advance rather than after a review"
)

# ---- the locus convention and the floors, both imported ------------------------------------------

#: The locus convention, imported from `cell2` rather than invented here, so a locus count in this
#: lane means exactly what it means there and in lane-repress2 and the three can be compared.
LOCUS_RULE = cell2.INDEPENDENT_LOCUS_RULE
LOCUS_SPAN = cell2.INDEPENDENT_LOCUS_SPAN
#: Carried with every locus count this lane reports. It is an operational grouping for deciding whether
#: a population has enough distinct places in the genome to measure, and it is **not** established
#: biological independence.
NOT_BIOLOGICAL_INDEPENDENCE = rp.NOT_BIOLOGICAL_INDEPENDENCE

#: Both floors imported, neither chosen by this lane, neither movable after a count.
LOCUS_FLOOR = fresh.LOCUS_FLOOR
POSITIVE_FLOOR = fresh.POSITIVE_FLOOR
FLOORS_CALL = (
    f"independent loci at least {LOCUS_FLOOR} (fresh.LOCUS_FLOOR, which is "
    f"cell2.POOLED_LOCUS_FLOOR) and links at least {POSITIVE_FLOOR} (fresh.POSITIVE_FLOOR). Both are "
    "imported from where the project fixed them before this lane existed, both are the floors "
    "lane-repress2's gate 1 was taken against, neither is set here and neither may be moved after a "
    "count"
)

# ---- the ladder ----------------------------------------------------------------------------------

#: The steps a cached pair has to clear, in order. Every test is blind to direction: a gene, a cell, an
#: element and whether an outcome is one of the two significant labels.
LADDER_STEPS = (
    "pairs_in_the_benchmark",
    "and_significant",
    "and_on_an_attributed_element",
    "and_a_rule_on_that_element_names_the_pairs_own_gene",
    "and_a_rule_on_that_element_and_gene_is_gated_on_the_pairs_own_cell",
)

#: How each step of this ladder lines up with lane-repress2's, so the two can be read step for step.
#: Theirs is anchored on a rule and walks towards a measurement; this one is anchored on a measurement
#: and walks towards a rule. The two meet at the gene step and at the cell step, which is where their
#: zero fell.
LADDER_MIRRORS = {
    "pairs_in_the_benchmark": (
        "no counterpart: lane-repress2's ladder starts at the rules of one axis, this one at the "
        "benchmark rows, and the two denominators are different populations by construction"
    ),
    "and_significant": (
        "lane-repress2's last step, `and_that_pair_is_significant`, applied first here: the same "
        "predicate (repress2.significant, membership in the two significant labels) on the same cached "
        "outcome field"
    ),
    "and_on_an_attributed_element": (
        "the other side of their `element_carries_a_crispri_pair` (357 for the repression axis): the "
        "same join, the same committed reciprocal-overlap rule, read from the pair instead of from the "
        "rule"
    ),
    "and_a_rule_on_that_element_names_the_pairs_own_gene": (
        "the other side of their `and_a_pair_on_the_rules_own_gene` (11 for the repression axis)"
    ),
    "and_a_rule_on_that_element_and_gene_is_gated_on_the_pairs_own_cell": (
        "the other side of their `and_in_the_rules_own_cell` (0 for the repression axis), which is "
        "where their population died"
    ),
}

LADDER_CALL = (
    "pairs, not links: how many cached CRISPRi pairs clear each step in turn, so that a zero or a "
    "short count at the end can be read as the step it fell at. A pair is counted once at each step "
    "however many elements or rules satisfy it, and it clears a step when any one of them does. Every "
    "step is blind - a gene, a cell, an element and membership in the two significant labels - and "
    "none of them reads the sign of an effect; the outcome breakdown is reported at each step's own "
    "denominator beside the count"
)

OVERLAP_CALL = (
    "a pair is on an attributed element when the element and the pair's tested interval reciprocally "
    f"overlap at measured.RECIPROCAL_OVERLAP = {ms.RECIPROCAL_OVERLAP} or more, computed by "
    "measured.reciprocal_overlap over the candidates measured.Layer.near returns at "
    f"measured.REACH = {ms.REACH}. This is the committed rule measured.rows applies and it is not "
    "varied here; it is the same join lane-repress2's `element_carries_a_crispri_pair` step took"
)

RULES_CALL = (
    "the compiled rules are genomeos.attribution.not_open_profile.rules unchanged, which mirrors "
    "compile.compile_chromosome block for block, so this run and lane-repress2's cannot describe "
    "different rules. A rule's cell is the context its `when: cell_type` names and a rule's gene is "
    "the gene it targets; both are read as the compiler wrote them"
)

# ---- the two populations the floors are applied to ----------------------------------------------

P1 = "measured_repression_links_a_changed_extractor_would_emit"
P2 = "increases_that_reach_an_existing_represses_target_rule"
POPULATIONS = (P1, P2)

P1_CALL = (
    "one link per (attributed element, gene, cell) for which that element's overlapping cached pairs "
    f"hold a pair whose outcome is {INCREASE!r} on that gene in that cell. This mirrors "
    "measured.rule_links exactly: that function emits one link per (gene, cell) carrying a pair the "
    "benchmark calls Regulated, and this counts one link per (gene, cell) carrying a significant "
    "increase. It is the population a changed extractor could reach, counted without changing one"
)
P1_CONTEST_RULE = (
    "a (gene, cell) of the same element that carries a regulated pair as well as a significant "
    "increase is **not** in P1. The committed extractor already emits an `activates` link there, so "
    "the direction at that (gene, cell) is contested by the measurement itself and a repression "
    "population built on it would be testing a disagreement rather than a repression. Such a "
    "(gene, cell) is counted in its own row, `contested_by_a_regulated_pair`, and is in neither "
    "population. This rule is fixed here, before any count, and is not relaxed if P1 falls short"
)
P1_SPLIT_RULE = (
    "a pair is counted whatever its evaluation partition, and the training/held-out breakdown is "
    "reported beside the count. The committed extractor already admits a link found only in held-out "
    "pairs and marks it with measured.HELDOUT_MARK, so excluding them here would make the "
    "counterfactual narrower than the thing it is a counterfactual of. No feature, fit or label is "
    "built on anything in this lane, so the split audit's bar is not engaged"
)
P2_CALL = (
    "one link per (attributed element, gene, cell) at the last ladder step for which at least one "
    f"compiled rule on that element, gene and cell has an activity axis of exactly {rp.REPRESSES!r}. "
    "This is lane-repress2's gate 1 population with the direction of the measurement now admitted: "
    "the increases that would test a repression call the compiled program already makes"
)
P2_RECONCILES_WITH = (
    "P2's link count must be at or below lane-repress2's `and_in_the_rules_own_cell` for "
    f"{rp.REPRESSES!r}, which is {INHERITED_LADDER['and_in_the_rules_own_cell']} genome-wide in "
    f"{INHERITED_FROM}, because P2 adds a direction filter to a population their step already counted. "
    "A larger number means the two enumerations disagree, and it is reported as a discrepancy to be "
    "resolved and never as a population. This reconciliation is registered before the count"
)
NEVER_POOLED = (
    "P1 and P2 are gated separately and are never pooled, with each other or with lane-repress2's "
    "activation populations, to reach a floor. They answer different questions: P1 whether a measured "
    "repression population could be built at all, P2 whether a measured increase can reach a "
    "repression call the compiled program already makes"
)

# ---- the readings, written before the outcome was seen ------------------------------------------

GO = (
    "at or above both registered floors: a measured repression population could be built from the "
    "increases the committed extractor discards. This lane stops here. It reads no direction, reports "
    "no rate, tests no repression call and draws no conclusion about whether the compiled program's "
    "repression calls are right: a feasibility count that rolls straight into the test it enables is "
    "how a floor gets moved after the outcome is seen. The test is registered as its own lane, and "
    "this lane's remaining deliverable is a written proposal naming the extractor change the count "
    "shows is worth making"
)
NO_GO = (
    "below a registered floor: a no-go. The compiled program's repression calls cannot be tested "
    "against a measured repression population, and this count is the reason they cannot, not a step "
    "towards a population that could. The count goes on the record, no floor is moved, no population "
    "is pooled with another to reach a floor, no direction is read and no rate is reported. A count "
    "short of a floor is reported as a no-go and never as close, promising, nearly enough, a good "
    "start or enough for a pilot; a short count has no encouraging reading here, and the margin it "
    "falls short by is reported because a reader is owed it and not because a small margin is better "
    "than a large one"
)
READINGS = (GO, NO_GO)

#: What follows each reading, fixed before the count so that neither outcome can be answered with a
#: new population, a new predicate or a smaller floor.
WHAT_FOLLOWS = {
    "go": (
        "stop at the verdict and report. The repression test is a separate lane with its own "
        "registration; the proposal naming the extractor change is this lane's"
    ),
    "no_go": (
        "stop at the verdict and report. The no-go lane-repress2 recorded, "
        f"{INHERITED_NO_GO!r}, is then structural: it is not an accident of the one arm they took, "
        "because the measurements a changed extractor could reach do not make a testable population "
        "either. Nothing is re-counted under a different predicate to find a population, and the "
        "floors stand"
    ),
}

# ---- the data-or-code attribution, fixed before the count ---------------------------------------

#: Which step a short population is attributed to, and whether that step is a property of the data or
#: of this repository's code. Fixed here, before the count, because deciding it afterwards is how an
#: inconvenient attribution gets re-read. "code" never means the code is wrong; it means the step is
#: decided by a choice this repository made and could in principle revisit.
DATA = "a property of the data: the measurements that exist"
CODE = "a property of this repository's code: a choice this repository made and could revisit"
ATTRIBUTION_OF_STEPS = {
    "and_significant": (
        DATA,
        "how many pairs the screens called significant, and how many of those were increases, is the "
        "benchmark's own content. No code in this repository decides it",
    ),
    "and_on_an_attributed_element": (
        CODE,
        "which elements the attribution layer emits, and the reciprocal-overlap rule "
        f"(measured.RECIPROCAL_OVERLAP = {ms.RECIPROCAL_OVERLAP}) that decides whether a tested "
        "interval is a measurement *of* one of them. A different element set or a weaker overlap bar "
        "(measured.ELIGIBILITY_RULE is the weaker one the project already defines) would give a "
        "different count",
    ),
    "and_a_rule_on_that_element_names_the_pairs_own_gene": (
        CODE,
        "which gene a compiled rule targets: for a predicted rule the one gene the deletion sweep's "
        "`predicted_coding` block names, out of every gene the screen measured against that element",
    ),
    "and_a_rule_on_that_element_and_gene_is_gated_on_the_pairs_own_cell": (
        CODE,
        "the `when: cell_type` gate the compiler writes on a rule, and the context "
        "compile.context derives for a predicted rule from the one most extreme track's tissue. This "
        "is the step lane-repress2's repression population died at",
    ),
    "the_contest_rule": (
        CODE,
        "P1_CONTEST_RULE, this lane's own pre-registered exclusion of a (gene, cell) whose "
        "measurement is contested by a regulated pair on the same element",
    ),
    "the_extractor": (
        CODE,
        "CAUSE_IN_CODE: that `measured.rule_links` raises no rule from a significant increase. This is "
        "the step that makes P1 a counterfactual at all, so it is attributed to code by construction "
        "and is not the finding of this count",
    ),
}
ATTRIBUTION_CALL = (
    "the step a population falls short at is attributed to the data or to this repository's code by "
    "ATTRIBUTION_OF_STEPS, which is fixed before the count. A population short at `and_significant` "
    "is a property of the data and the no-go is final for this benchmark; a population short at any "
    "later step is a property of this repository's code and names what would have to change. The "
    "distinction is the value of the lane and it is not decided after the number is seen"
)

# ---- what the result cannot establish ------------------------------------------------------------

CANNOT_ESTABLISH = (
    "a count of the measurements a changed extractor could reach is not a measurement of any "
    "repression call. Whatever it says, nothing here establishes that any compiled repression call is "
    "right or wrong, that any of these increases is a repression mechanism, or that an element whose "
    "silencing raises a gene represses it: an increase on silencing is one screen's signed effect in "
    "one cell and this lane reads no mechanism from it. A locus count is cell2's operational grouping "
    "and never established biological independence. P1 is a counterfactual: the links it counts do not "
    "exist in the compiled program and this lane does not create them. Clearing a floor would say "
    "that a population exists to test, never that a test of it would come out one way; falling short "
    "of a floor would say that this benchmark cannot carry the test, never that the repression calls "
    "are right. No assay other than the CRISPRi benchmark is read, so nothing here is a statement "
    "about repression evidence in general"
)


# ---- the predicates the counting script applies --------------------------------------------------


def significant(pair: Any) -> bool:
    """Whether a cached pair carries a significant signed effect, read without its sign.

    `repress2.significant` itself, imported rather than reimplemented, so the two ladders cannot apply
    different predicates to the same cached field.
    """
    return rp.significant(pair)


def is_increase(pair: dict[str, Any]) -> bool:
    """Whether a cached pair's outcome is exactly the significant-increase label.

    The eligibility predicate of this lane, applied to the cached row's own `outcome` field. It is read
    only to break a step's denominator down and to build P1 and P2; no position step uses it.
    """
    return pair.get("outcome") == INCREASE


def base_element_id(element: str) -> str:
    """The attributed element a rule sits on, with the `_measured` suffix `not_open_profile.rules`
    gives a measured rule's element removed.

    The same convention `repress2._element_id` applies; `tests/test_increases.py` holds the two
    together on the ids this repository actually produces, so they cannot drift.
    """
    return element[: -len("_measured")] if element.endswith("_measured") else element


def gate(keys: list[cell2.LocusKey], population: str, call: str) -> dict[str, Any]:
    """One population against both imported floors. No direction of any effect is read here.

    `keys` is one `cell2.LocusKey` per link, which is the only path into the grouping and carries
    nothing but the cell, the chromosome, the element's interval and the gene.
    """
    links = len(set(keys))
    loci = cell2.count_loci(sorted(set(keys)), LOCUS_SPAN)
    met = links >= POSITIVE_FLOOR and loci >= LOCUS_FLOOR
    short = []
    if links < POSITIVE_FLOOR:
        short.append(f"links {links}, short of {POSITIVE_FLOOR} by {POSITIVE_FLOOR - links}")
    if loci < LOCUS_FLOOR:
        short.append(f"independent loci {loci}, short of {LOCUS_FLOOR} by {LOCUS_FLOOR - loci}")
    return {
        "population": population,
        "what_it_is": call,
        "links": links,
        "independent_loci": loci,
        "floors": {"links": POSITIVE_FLOOR, "independent_loci": LOCUS_FLOOR},
        "meets_both_floors": met,
        "reading": GO if met else NO_GO,
        "short_by": short,
        "independent_locus_rule": LOCUS_RULE,
        "not_biological_independence": NOT_BIOLOGICAL_INDEPENDENCE,
    }


def breakdown(outcomes: Any) -> dict[str, int]:
    """The exhaustive outcome breakdown of one denominator, every label present even at zero."""
    return {o: int(outcomes.get(o, 0)) for o in ALL_OUTCOMES}


def attribution_of(step: str) -> dict[str, str]:
    """Where a population that falls short at `step` is attributed, as registered before the count."""
    kind, why = ATTRIBUTION_OF_STEPS[step]
    return {"step": step, "data_or_code": kind, "why": why}


def registration() -> dict[str, Any]:
    """Everything this lane fixed before it counted anything, as one payload.

    `scripts/increase_register.py` writes it and `scripts/increase_population.py` carries it into the
    result nested, so the registration and the run cannot drift apart.
    """
    return {
        "result": "increase_registration",
        "lane": "lane-increase",
        "question": (
            "can a repression population be built at all, given that a significant increase raises no "
            "measured rule"
        ),
        "follows_from": {
            "lane": "lane-repress2",
            "commits": ["7a5dd8b", "fd07258", "37a256c", "e9cc033"],
            "registered_no_go_carried_verbatim": INHERITED_NO_GO,
            "their_gate_1_reading": rp.GATE_NO_GO,
            "their_blind_ladder": INHERITED_LADDER,
            "their_measured_layer_count": INHERITED_MEASURED_LAYER,
            "the_structural_cause_they_recorded": STRUCTURAL_CAUSE,
        },
        "the_cause_in_code": CAUSE_IN_CODE,
        "no_extractor_changed": NO_EXTRACTOR_CHANGED,
        "eligibility_predicate": ELIGIBILITY_CALL,
        "outcome_breakdown_at_every_denominator": BREAKDOWN_CALL,
        "independent_locus": {
            "rule": LOCUS_RULE,
            "span": LOCUS_SPAN,
            "imported_from": "genomeos.attribution.cell2.INDEPENDENT_LOCUS_RULE",
            "not_biological_independence": NOT_BIOLOGICAL_INDEPENDENCE,
        },
        "floors": {
            "call": FLOORS_CALL,
            "links": POSITIVE_FLOOR,
            "links_imported_from": "genomeos.attribution.fresh.POSITIVE_FLOOR",
            "independent_loci": LOCUS_FLOOR,
            "independent_loci_imported_from": (
                "genomeos.attribution.fresh.LOCUS_FLOOR, which is "
                "genomeos.attribution.cell2.POOLED_LOCUS_FLOOR"
            ),
            "neither_chosen_here": True,
        },
        "ladder": {
            "call": LADDER_CALL,
            "steps": list(LADDER_STEPS),
            "mirrors_lane_repress2": LADDER_MIRRORS,
            "overlap_rule": OVERLAP_CALL,
            "rules": RULES_CALL,
        },
        "populations": {
            "order": list(POPULATIONS),
            P1: {"call": P1_CALL, "contest_rule": P1_CONTEST_RULE, "split_rule": P1_SPLIT_RULE},
            P2: {"call": P2_CALL, "reconciles_with": P2_RECONCILES_WITH},
            "never_pooled": NEVER_POOLED,
        },
        "readings_before_the_outcome_was_seen": {
            "go": GO,
            "no_go": NO_GO,
            "there_is_no_third": True,
            "what_follows": WHAT_FOLLOWS,
        },
        "data_or_code_attribution": {
            "call": ATTRIBUTION_CALL,
            "per_step": {k: {"data_or_code": v[0], "why": v[1]} for k, v in ATTRIBUTION_OF_STEPS.items()},
        },
        "cannot_establish": CANNOT_ESTABLISH,
        "requests": 0,
        "money": "none: every input is already on disk",
        "strong_effect_for_reference": nop.STRONG_EFFECT,
    }
