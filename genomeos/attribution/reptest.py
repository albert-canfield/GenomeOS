# SPDX-License-Identifier: AGPL-3.0-or-later
"""The answerability count on the second extractor's increase-derived links, taken BEFORE any test.

**What this lane is, and why a count comes first.** lane-inhibits built extractor v2 and emitted 48
increase-derived links, deliberately applying neither floor so that the lane testing them could not
pick its own. The obvious next step is a direction test. It is not the step this lane takes, because
the project already has the counter-example on the record: lane-direction's population held 48
MEASURED increase links and only 21 ANSWERABLE ones, and its gate then refused the comparison. So the
first question about the 48 is not what direction they read but how many of them can be asked at all,
and this module is that question's registration and its verdict rule.

**The one thing this lane does differently from lane-direction, and the reason it is not redundant.**
lane-direction enumerated its own population: a hand-rolled join over the cached CRISPRi pairs, one
link per (attributed element, gene, cell) whose pairs hold a significant effect. That join is
lane-increase's and it is correct, but it is not the extractor. This lane counts THE EXTRACTOR'S OWN
OUTPUT: every record `measured.rule_links_detail(row, measured.EXTRACTOR_V2)` returns, over the rows
`measured.rows` builds. So the answerability figure here is a property of the population v2 actually
emits, and `POPULATION_IS_THE_EXTRACTORS_OUTPUT` says so in the result. If the two populations
coincide that is a finding about the extractor and it is reported as one; it is not assumed.

**The answerability predicate is lane-direction's, imported and not restated**, so that a count here
and a count there mean the same thing: the cached signed deletion value for THAT link's own gene on
THAT link's own cell track, from the finished sweep's per-element response cache, with nothing
substituted for a cell, a gene, an element or a cache the answer is not in.

**What this lane refines.** lane-direction reported 18 of its 39 K562 increase links as `absent`, one
bucket. Three different things hide in one word there: the element may be missing from the cache
archive entirely, the gene may never have been in the scorer's window at that element, or the gene may
have been scored and carry no entry on that cell's own track. `genomeos.attribution.targets` already
names the last two (`NOT_IN_WINDOW`, `NOT_ON_TRACK`) and this lane counts them apart, because a gap
that is a missing run is a different gap from a gene outside a window and only one of them could ever
be closed by reading more of what is already on disk. `CAUSES` is the exhaustive, disjoint list and
the result asserts the causes sum to the population.

**The floors are imported and neither is applied by choice.** `fresh.POSITIVE_FLOOR = 30` links and
`fresh.LOCUS_FLOOR`, which is `cell2.POOLED_LOCUS_FLOOR = 20`, independent loci. They were fixed
before this lane existed, they are compared against the ANSWERABLE count and never against the
measured count, and neither moves after the count is seen. A count short of a floor is a NO-GO and is
reported in those words.

**The statistic, addressed before the count rather than after it.** lane-direction proved, pinned by a
test over 200 cases and written down before any value was read, that the excess of one arm's agreement
minus the other's over a sign-shuffled baseline is identically
`2 * (balanced_accuracy - 0.5) * (n_u - n_d) / (n_d + n_u)`. That algebra is carried here verbatim as
`COLLAPSE` and its premise is CHECKED on this population rather than assumed: `collapse_applies`
reports whether the two arms' answered denominators are unequal and in which direction, from this
lane's own counts. The statistic this lane would register is therefore balanced accuracy against 0.5
and never a between-arm difference, and `IDENTIFIABILITY` states why in advance of knowing whether any
test is taken at all.

**R2, binding and in code.** A measured increase on knockdown is an increase on knockdown. It is never
a silencer, never a repressor and never evidence of a repression mechanism, because an increase under
CRISPRi can be indirect - KRAB spread onto a neighbouring promoter, competition between promoters for
a shared enhancer. `increase_links.check_no_mechanism_claim` is called on every record this lane
handles and is not weakened, reimplemented or bypassed.

**Budget.** No direction is read unless both floors clear. 0 model requests, no money, no downloads,
no network. The cached deletion answers are read; none is bought. An answer that is not cached is a
gap to count and never a request to send.
"""

from __future__ import annotations

from typing import Any

from genomeos.attribution import cell2, fresh, targets
from genomeos.attribution import direction_link as dl
from genomeos.attribution import increase_links as il
from genomeos.attribution import increases as inc
from genomeos.attribution import measured as ms

REGISTERED = "2026-10-02"
LANE = "lane-reptest"

# ---- what is inherited, imported so that two lanes cannot come to mean different things -----------

#: lane-inhibits' registration, which built the population this lane counts.
INHERITED_EXTRACTOR_REGISTRATION = "data/results/increase_links_registration.json"
INHERITED_COMMITS = {
    "increase_links registration code": "fabf43c",
    "increase_links registration result": "9548040",
    "the extractor, v2 added and v1 byte-identical": "ba32829",
    "the record of what v2 adds": "25aac0d",
    "lane-increase's registration, committed before any count": "e1c2671",
    "lane-direction's registration, committed before any direction": "0779c97",
}
INHERITED_FILES = (
    "data/results/increase_links_registration.json",
    "data/results/increase_population.json",
    "data/results/increase_registration.json",
    "data/results/direction_link.json",
    "data/results/direction_link_registration.json",
)

#: Every figure of lane-inhibits' and lane-direction's that bounds this lane. None is recomputed as a
#: finding of this lane, and none may be quoted under a noun its own breakdown does not support.
INHERITED_FIGURES = {
    "v2_links_total": 260,
    "v2_activates_links_from_v1": 212,
    "v2_inhibits_links_added": 48,
    "v1_links_removed_by_v2": 0,
    "conflicts_today": 0,
    "increase_links_by_cell": {"K562": 39, "WTC11": 6, "HCT116": 3},
    "increase_links_by_split": {"training": 25, "heldout": 23},
    "increase_links_independent_loci_as_lane_increase_counted_them": 33,
    "lane_direction_increase_arm_measured": 48,
    "lane_direction_increase_arm_in_the_cell_read": 39,
    "lane_direction_increase_arm_answered": 21,
    "lane_direction_increase_arm_absent": 18,
    "lane_direction_increase_arm_independent_loci_of_the_answered": 13,
}

#: The outcome of lane-direction's gate, carried word for word, because it is the reason this lane
#: exists and it may not be softened into a near miss on the way here.
INHERITED_NO_GO = dl.GATE_NO_GO
INHERITED_INCREASE_GO = inc.GO

# ---- the population: the extractor's own output, not a join rebuilt beside it ---------------------

EXTRACTOR = ms.EXTRACTOR_V2
POPULATION_IS_THE_EXTRACTORS_OUTPUT = (
    "the population counted here is every record `measured.rule_links_detail(row, "
    f"{EXTRACTOR!r})` returns whose `derived_from` is the increase label, over the rows "
    "`measured.rows` builds for each chromosome from `targets.attributed`. It is read out of the "
    "extractor and is not a join rebuilt beside it: lane-direction enumerated an equivalent "
    "population by hand from the cached pairs, and a count taken on a hand-rolled copy is a count "
    "about the copy. Whether the extractor's output and that hand-rolled population agree is "
    "reported as a comparison with both counts named, and it is not assumed in either direction. No "
    "compiled rule is consulted, no predicted value enters the eligibility step, and the extractor is "
    "called exactly as committed with no argument of its own changed"
)

POPULATION_IS_NOT_COVERAGE = il.BY_CONSTRUCTION
DECREASE_ARM_IS_NOT_MODEL_COVERAGE = dl.DECREASE_ARM_IS_NOT_MODEL_COVERAGE
COUNT_IS_NOT_COVERAGE = (
    "neither the 48 nor the answerable subset of them is a coverage figure. Each link exists because a "
    "significant increase was measured on an element the compiled program had already attributed, and "
    "the link then took its gene and its cell FROM THAT MEASUREMENT: the population is the extractor's "
    "construction. An answerable count is smaller still and is a statement about what the cached "
    "deletion sweep holds, not about what the model covers. A reader may not read either number as "
    "coverage of the increase links, of the genes behind them or of anything else"
)

# ---- the answerability predicate, imported from lane-direction -----------------------------------

CACHE = dl.CACHE
CACHED_CELLS = dl.CACHED_CELLS
CELLS_NOT_CACHED = dl.CELLS_NOT_CACHED
ANSWERABLE_CALL = (
    "a link is ANSWERABLE IN ITS OWN CELL when the finished deletion sweep's per-element response "
    f"cache {CACHE} holds, under that link's own element id, a gene entry whose name is exactly that "
    "link's gene, carrying an entry for that link's own cell track, and that entry's value is not "
    "exactly 0.0. This is lane-direction's own predicate, imported from "
    "`direction_link.PREDICTED_CALL` and `direction_link.SIGN_CALL` rather than restated, so a count "
    "here and a count there mean the same thing. NOTHING IS SUBSTITUTED for a missing answer: not "
    "another cell's value, not another gene's value, not a neighbouring element's value, not a pooled "
    "or marginal value, and not the partial second cache "
    "data/knowledge/alphagenome/elements_hct116, which exists over five chromosomes from a different "
    "run and does carry an HCT116 track. Mixing that cache in would widen one arm of a comparison by "
    "substitution under another name, so it is named here and NOT read"
)
PREDICTED_CALL = dl.PREDICTED_CALL
SIGN_CALL = dl.SIGN_CALL

# ---- the causes, exhaustive and disjoint ----------------------------------------------------------

ANSWERABLE = "answerable_in_its_own_cell"
CELL_NOT_CACHED = "the_cache_carries_no_track_of_this_cell"
ELEMENT_ABSENT = "the_element_is_absent_from_the_cache_archive"
GENE_NOT_IN_WINDOW = "the_gene_is_not_in_the_scorers_window_at_this_element"
GENE_NOT_ON_TRACK = "the_gene_was_scored_but_carries_no_entry_on_this_cells_track"
VALUE_EXACTLY_ZERO = "the_cached_value_is_exactly_zero_so_it_carries_no_direction"

#: The exhaustive, disjoint causes. Every link of the population falls in exactly one and the result
#: asserts they sum to the population, because the project has twice quoted a count under a noun its
#: own breakdown did not support (docs/LESSONS.md, "A count is not a measurement of the thing you want
#: to count"); 198 loci where the model says nothing were quoted as 198 places it is wrong when 66
#: held a measured decrease, and 48 measured links cleared a floor while only 21 were answerable.
CAUSES = (
    ANSWERABLE,
    CELL_NOT_CACHED,
    ELEMENT_ABSENT,
    GENE_NOT_IN_WINDOW,
    GENE_NOT_ON_TRACK,
    VALUE_EXACTLY_ZERO,
)
CAUSE_TEXT = {
    ANSWERABLE: (
        "the cache holds a value for this link's own gene on this link's own cell track and the value "
        "is not exactly zero, so the link carries a predicted direction and could be asked"
    ),
    CELL_NOT_CACHED: (
        f"the link's cell is not one of the tracks the cached answers carry ({CACHED_CELLS!r}). The "
        f"sweep holds no track of {CELLS_NOT_CACHED!r} at all, so no amount of reading what is on disk "
        "closes this gap and nothing is substituted for it. These links are absent from any test "
        "rather than negative in it"
    ),
    ELEMENT_ABSENT: (
        "the link's element id is not in the per-chromosome cache archive, so the sweep holds no "
        "record of this element. This is a gap in the run and is reported as one; it is not a "
        "statement about the element"
    ),
    GENE_NOT_IN_WINDOW: (
        f"the element's cached record exists but names no gene matching the link's gene: "
        f"{targets.NOT_IN_WINDOW}. The CRISPRi benchmark tested a pair the scorer's window did not "
        "reach, so there is nothing cached to read and no nearest gene stands in for it"
    ),
    GENE_NOT_ON_TRACK: (
        f"the gene is in the element's cached record but carries no entry for the link's own cell: "
        f"{targets.NOT_ON_TRACK}. The gene was scored and this cell was not scored for it, which is a "
        "different gap from a gene outside the window and is counted apart from it"
    ),
    VALUE_EXACTLY_ZERO: (
        "the cached value is exactly 0.0, which carries no sign under the imported sign rule and is "
        "therefore neither an agreement nor a disagreement. Excluded from the answerable count and "
        "counted on its own, never defaulted to a direction"
    ),
}
CAUSES_ARE_EXHAUSTIVE = (
    "every link of the population falls in exactly one cause and the causes sum to the population "
    "count. The result asserts the sum and fails if it does not hold, so the breakdown cannot drift "
    "from the total it explains. The answerable count is reported BESIDE the population count and "
    "never instead of it"
)

# ---- the floors, imported and never chosen here --------------------------------------------------

POSITIVE_FLOOR = fresh.POSITIVE_FLOOR
LOCUS_FLOOR = fresh.LOCUS_FLOOR
FLOORS = {
    "answerable_links_at_least": POSITIVE_FLOOR,
    "answerable_links_floor_imported_from": "genomeos.attribution.fresh.POSITIVE_FLOOR",
    "independent_loci_at_least": LOCUS_FLOOR,
    "independent_loci_floor_imported_from": (
        "genomeos.attribution.fresh.LOCUS_FLOOR, which is genomeos.attribution.cell2.POOLED_LOCUS_FLOOR"
    ),
    "applied_to": (
        "the ANSWERABLE links and their independent loci, never the measured links. A floor compared "
        "against a count of what exists rather than a count of what can be asked is the error "
        "lane-direction's gate caught: 48 measured links cleared 30 and 21 answerable ones did not"
    ),
    "neither_is_chosen_here": (
        "both floors were fixed before this lane existed. lane-inhibits deliberately applied neither so "
        "that this lane could not pick its own, and neither moves after the count is seen - not down to "
        "admit a population, not up, and not sidestepped by pooling the arms, by widening the "
        "answerable set or by taking a pilot on a short one"
    ),
}
LOCUS_RULE = cell2.INDEPENDENT_LOCUS_RULE
LOCUS_SPAN = cell2.INDEPENDENT_LOCUS_SPAN
NOT_BIOLOGICAL_INDEPENDENCE = inc.NOT_BIOLOGICAL_INDEPENDENCE

GATE_CALL = (
    f"the gate is taken on the ANSWERABLE increase-derived links: their count must be at or above "
    f"{POSITIVE_FLOOR} and their independent loci at or above {LOCUS_FLOOR}, both imported. Both "
    "floors or neither: a population that clears one and misses the other is a no-go. The arms are "
    "never pooled to reach a floor and the measured count is never substituted for the answerable one"
)
GATE_NO_GO = (
    "the answerable increase-derived links are below an imported floor: a NO-GO. No direction is read, "
    "no agreement rate is computed, no interval is taken, no registration of a test is written and no "
    "floor is moved. The answerable set is not widened by substituting a cell, a gene, an element or a "
    "second cache, and the two arms are not pooled to reach a floor. The count goes on the record with "
    "the margin it falls short by, because a reader is owed the margin and NOT because a small margin "
    "is better than a large one. A count short of a floor is reported as a no-go and NEVER as close, "
    "promising, nearly enough, a good start or enough for a pilot"
)
GATE_GO = (
    "the answerable increase-derived links are at or above both imported floors: a direction test may "
    "be registered and then run, in that order, with the registration committed before any direction "
    "is read. This lane takes no direction on the strength of the gate alone"
)
READINGS = (GATE_GO, GATE_NO_GO)
THERE_IS_NO_THIRD = True

#: The words no result, summary or commit message of this lane may use of a count short of a floor.
#: Carried from lane-direction, which carried it from lane-increase.
FORBIDDEN_OF_A_SHORT_COUNT = dl.FORBIDDEN_OF_A_SHORT_OR_UNDETECTED_OUTCOME

# ---- the statistic, addressed before the count ----------------------------------------------------

COLLAPSE = dl.COLLAPSE
TEST_STATISTIC = dl.TEST_STATISTIC
CHANCE = dl.CHANCE
IDENTIFIABILITY = (
    "the statistic a test on this population would be registered on is BALANCED ACCURACY against "
    f"{CHANCE}, and never the difference between the two arms' agreement rates. The reason is not "
    "taste and it is not re-derived here: lane-direction proved it, before reading any value and "
    "pinned by a test over 200 cases, and the algebra is carried verbatim as `COLLAPSE`. The excess of "
    "(one arm's agreement minus the other's) over a sign-shuffled baseline is identically "
    "2 * (balanced_accuracy - 0.5) * (n_u - n_d) / (n_d + n_u), so with equal arms it is exactly zero "
    "and with the decreases arm larger a POSITIVE excess requires balanced accuracy BELOW 0.5 - a "
    "one-sided test of that difference is a test of the model doing worse than chance. That "
    "derivation is not assumed to transfer: `collapse_applies` recomputes its premise from THIS "
    "lane's own answered denominators and reports which arm is larger and by how much, so the claim "
    "rests on this population's counts and not on the last lane's. With one marginal sign rate there "
    "is exactly ONE skill number on a population of this shape and it is balanced accuracy"
)
IDENTIFIABILITY_IS_CHECKED_NOT_ASSUMED = (
    "`collapse_applies` is the check, and it can refuse. It reports the two arms' answered "
    "denominators, whether they are equal, and what the collapse implies for this population. If the "
    "arms were equal the between-arm excess would be identically zero and a difference statistic would "
    "be unidentifiable here, which is a stop condition and not a detail to work around"
)


# ---- the functions the result is built from -------------------------------------------------------


def cause_of(record: dict[str, Any] | None, gene: str, cell: str) -> tuple[str, float | None]:
    """One link's answerability cause and, when it has one, its cached signed value.

    `record` is the link's own element record out of the per-element response cache, or None when the
    archive holds no such element. The gene and the cell are the link's OWN, and no other gene, cell,
    element or cache is consulted: that is the whole predicate, and the branches below are the
    exhaustive causes `CAUSES` names, in the order that makes them disjoint.
    """
    if cell not in CACHED_CELLS:
        return CELL_NOT_CACHED, None
    if record is None:
        return ELEMENT_ABSENT, None
    entry = next((g for g in record.get("genes") or [] if g.get("gene") == gene), None)
    if entry is None:
        return GENE_NOT_IN_WINDOW, None
    value = (entry.get("by_cell") or {}).get(cell)
    if value is None:
        return GENE_NOT_ON_TRACK, None
    if dl.sign_of(value) is None:
        return VALUE_EXACTLY_ZERO, float(value)
    return ANSWERABLE, float(value)


def causes_reconcile(counts: dict[str, int], population: int) -> bool:
    """Whether the per-cause counts are exhaustive over the population. CAUSES_ARE_EXHAUSTIVE."""
    return set(counts) <= set(CAUSES) and sum(counts.values()) == population


def independent_loci(links: list[dict[str, Any]]) -> int:
    """The independent loci of a set of links, by cell2's rule through lane-direction's one path."""
    return len(set(dl.loci_of(links))) if links else 0


def verdict(answerable: list[dict[str, Any]], population: int) -> dict[str, Any]:
    """The gate on the answerable links against both imported floors, with the margin of each miss.

    The floors are compared against the ANSWERABLE count and its loci. `population` is carried into
    the verdict so that the measured count is always reported beside the answerable one and a reader
    can never meet one without the other.
    """
    n = len(answerable)
    loci = independent_loci(answerable)
    short: list[str] = []
    if n < POSITIVE_FLOOR:
        short.append(
            f"answerable links {n}, short of the imported floor {POSITIVE_FLOOR} by {POSITIVE_FLOOR - n}"
        )
    if loci < LOCUS_FLOOR:
        short.append(
            f"independent loci {loci}, short of the imported floor {LOCUS_FLOOR} by {LOCUS_FLOOR - loci}"
        )
    met = not short
    return {
        "population_measured_links": population,
        "answerable_links": n,
        "independent_loci_of_the_answerable_links": loci,
        "floors": dict(FLOORS),
        "gate_call": GATE_CALL,
        "meets_both_floors": met,
        "short_by": short,
        "reading": GATE_GO if met else GATE_NO_GO,
        "independent_locus_rule": LOCUS_RULE,
        "not_biological_independence": NOT_BIOLOGICAL_INDEPENDENCE,
        "what_follows": (
            "register a direction test, commit the registration, then run it once"
            if met
            else "nothing further is registered and no direction is read"
        ),
    }


def collapse_applies(n_down: int, n_up: int) -> dict[str, Any]:
    """Whether lane-direction's collapse identity bites on THIS lane's own answered denominators.

    The identity is `COLLAPSE`, proved there and carried here unchanged. What is computed here is its
    PREMISE on this population: the two answered denominators and whether they are equal. Equal arms
    make the between-arm excess identically zero, which would leave a difference statistic
    unidentifiable, and an unequal pair fixes the sign of the balanced accuracy a positive excess would
    require. Either way the statistic a test would be registered on is balanced accuracy.
    """
    total = n_down + n_up
    factor = None if total == 0 else round((n_up - n_down) / total, 6)
    if total == 0:
        implication = "no answered link in either arm, so no statistic is identifiable at all"
    elif n_up == n_down:
        implication = (
            "the two arms are the same size, so the excess of the between-arm difference over the "
            "sign-shuffled baseline is IDENTICALLY ZERO and a difference statistic is unidentifiable "
            "on this population: the whole difference is the model's own marginal sign rate"
        )
    elif n_up < n_down:
        implication = (
            f"the decreases arm is the larger one ({n_down} against {n_up}), so (n_u - n_d) is "
            "negative and a POSITIVE excess of the between-arm difference requires balanced accuracy "
            "BELOW 0.5. A one-sided test of that difference against the sign shuffle would be a test "
            "of the model doing WORSE than chance, which is not the question"
        )
    else:
        implication = (
            f"the increases arm is the larger one ({n_up} against {n_down}), so (n_u - n_d) is "
            "positive and a positive excess of the between-arm difference requires balanced accuracy "
            "ABOVE 0.5. The excess is still a reparametrisation of balanced accuracy and carries no "
            "second skill number, so the statistic is unchanged"
        )
    return {
        "identity": COLLAPSE,
        "premise_recomputed_on_this_population": True,
        "answered_decrease_links": n_down,
        "answered_increase_links": n_up,
        "arms_are_equal": n_up == n_down,
        "arm_size_factor_n_u_minus_n_d_over_n": factor,
        "implication": implication,
        "statistic_a_test_would_be_registered_on": TEST_STATISTIC,
        "difference_is_identifiable_as_a_skill_statistic": False,
        "why": IDENTIFIABILITY,
        "checked_not_assumed": IDENTIFIABILITY_IS_CHECKED_NOT_ASSUMED,
    }


def registration() -> dict[str, Any]:
    """Everything fixed before the count, as the result states it. Nothing is computed here."""
    return {
        "lane": LANE,
        "registered": REGISTERED,
        "what_is_registered": (
            "the answerability count on extractor v2's increase-derived links and the verdict rule "
            "against both imported floors. No direction test is registered: a test is registered only "
            "if the answerable count clears both floors, and then in its own commit before any "
            "direction is read"
        ),
        "acts_on": {
            "lane": "lane-inhibits",
            "registration": INHERITED_EXTRACTOR_REGISTRATION,
            "files": list(INHERITED_FILES),
            "commits": dict(INHERITED_COMMITS),
            "figures": dict(INHERITED_FIGURES),
            "lane_direction_no_go_carried_verbatim": INHERITED_NO_GO,
            "lane_increase_go_carried_verbatim": INHERITED_INCREASE_GO,
        },
        "population": POPULATION_IS_THE_EXTRACTORS_OUTPUT,
        "extractor": EXTRACTOR,
        "answerable_predicate": ANSWERABLE_CALL,
        "predicted_direction": PREDICTED_CALL,
        "sign_rule_zero_and_absent": SIGN_CALL,
        "causes": {c: CAUSE_TEXT[c] for c in CAUSES},
        "causes_are_exhaustive": CAUSES_ARE_EXHAUSTIVE,
        "floors": dict(FLOORS),
        "gate_call": GATE_CALL,
        "readings": list(READINGS),
        "there_is_no_third": THERE_IS_NO_THIRD,
        "forbidden_of_a_short_count": list(FORBIDDEN_OF_A_SHORT_COUNT),
        "statistic": IDENTIFIABILITY,
        "statistic_is_checked_not_assumed": IDENTIFIABILITY_IS_CHECKED_NOT_ASSUMED,
        "collapse_identity_carried_verbatim": COLLAPSE,
        "binding_wording_r2": il.R2_WORDING,
        "prohibition_is_enforced_by": (
            "genomeos.attribution.increase_links.check_no_mechanism_claim, called on every record this "
            "lane handles. It is not weakened, not reimplemented and not bypassed"
        ),
        "by_construction": POPULATION_IS_NOT_COVERAGE,
        "count_is_not_coverage": COUNT_IS_NOT_COVERAGE,
        "decrease_arm_is_not_model_coverage": DECREASE_ARM_IS_NOT_MODEL_COVERAGE,
        "conflict_is_not_a_finding": il.CONFLICT_IS_NOT_A_FINDING,
        "independent_locus_rule": LOCUS_RULE,
        "not_biological_independence": NOT_BIOLOGICAL_INDEPENDENCE,
        "budget": (
            "0 model requests, no money, no downloads, no network. The cached deletion answers are "
            "read; none is bought. An answer that is not cached is a gap to COUNT and never a request "
            "to send"
        ),
        "order": [
            "1. the population, read out of the extractor under v2 and blind to every cached value",
            "2. the answerability pass: one cause per link, exhaustive and disjoint, summing to the "
            "population",
            "3. the verdict against both imported floors on the ANSWERABLE links and their loci",
            "4. if either floor is missed: stop. Nothing further is registered and no direction is read",
            "5. only if both clear: a direction test registered and committed, then run once",
        ],
    }
