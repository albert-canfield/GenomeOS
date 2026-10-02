# SPDX-License-Identifier: AGPL-3.0-or-later
"""What a re-score would cost to make the increase-derived direction test answerable, at zero requests.

**This is a costing and not a purchase.** lane-reptest closed a no-go: of extractor v2's 48
increase-derived links, 21 are ANSWERABLE IN THEIR OWN CELL over 13 independent loci, against imported
floors of 30 links (`fresh.POSITIVE_FLOOR`) and 20 loci (`cell2.POOLED_LOCUS_FLOOR`) - short by 9 and
by 7. Its causes sum to 48 and three of its six buckets are empty: no link fails for a missing element
record, none for a gene scored on other tracks but not this cell's, and none for an exactly-zero value.
So 27 of the 48 are beyond what the finished sweep contains, 18 because the tested gene was not among
the genes the sweep's record names at that element and 9 because the sweep holds no track of the link's
cell. This module says what it would COST to close those two gaps, in model requests, so that the
decision can be put with a number. Nothing here sends a request, and no figure here is a direction, an
agreement rate, an interval or a test.

**The shortfall is not the population, and that is the first decision this module records.** The no-go
was a BLINDED FEASIBILITY no-go: lane-reptest read no direction, computed no agreement rate and took no
interval, so no scientific outcome exists for a later purchase to have been chosen around, and topping
the data up is admissible. What is not admissible is buying the gap: a set defined as "the links that
were missing" is defined by the count that refused the test. So the costed population is stated by a
rule that can be written down without the no-go - EVERY measured significant link, both directions, in a
cell the sweep kept no track of - and that rule is applied to every link it admits, including links in
the arm that was not short. Which cells those are is computed against `CACHED_CELLS` rather than read
from `direction_link.CELLS_NOT_CACHED`, because that constant names the two cells the increase arm
carried and misses Jurkat, which the decrease arm carries and the sweep covers no better. The window gap
gets the same treatment, and the answer there is harder than an admissibility objection: see
`WINDOW_ARM_IS_CLOSED_ON_A_SERVICE_LIMIT`.

**The request unit is read from the committed chain, not carried over.** `enhancer_target` says "One
request per element" and the chain that produced the sweep says "about 12,000 elements on chr21, one
request each, cached per element"; `score_element` returns a cached answer WITHOUT a request and
otherwise calls the scorer once; and the scorer makes exactly one `client.score_variant` call whose one
response carries every gene in the window on every track. The cells are not a request dimension at all:
`aggregate(effects, cells)` decides which cell tracks are KEPT per gene after the answer is in hand, so
asking for a cell that was not kept costs the element's request again and not a request per cell.
`REQUEST_UNIT` states that and `REQUEST_UNIT_READ_FROM` cites where each clause was read. The HCT116
precedent used the same unit on the same code path - "One deletion request per registry element
overlapping a held-out HCT116 pair", with `KEPT_CELLS = (*CELLS, CELL)` - and the AstroREG registration's
"one request per registry element overlapping a covered pair" is the same unit with a different
denominator, so it is checked rather than assumed.

**The window is already at the service's maximum, and that is what makes most of the first arm
unbuyable.** The committed scorer resizes the variant's reference interval to
`dna_client.SEQUENCE_LENGTH_1MB`, and 1,048,576 is the LARGEST of the four lengths the client accepts;
`score_variant` validates the length against that set and raises for anything else. A link whose tested
gene's transcription start sits further from the element than that window reaches therefore cannot be
bought at any price: there is no wider window to buy. `WINDOW_IS_AT_THE_MAXIMUM` says so and
`window_needed` computes, per link and from the benchmark's own TSS coordinate, the width that would be
required.

**A gene missing from a record is not always a window gap.** `targets.NOT_IN_WINDOW` is the cause
lane-reptest counted, and it is the cause as the predicate defines it: the record names no gene matching
the link's gene. Two different things can produce that. The gene may lie outside the window, which is a
window gap. Or the gene may lie inside the window under a different symbol, which is a naming gap and
costs no request at all. This module separates them by coordinate and reports each count with its
denominator, because the first is a purchase that cannot be made and the second is not a purchase.

**The ceiling is an upper bound and is labelled one.** A request delivers a value; the registered
predicate needs a value that is not exactly zero, on that link's own cell track, in the sweep's own
cache. So the most a purchase can add is one answerable link per purchasable link, and the ceiling is
compared against the imported floors to say whether the money would buy a testable number at all.
`CEILING_IS_NOT_A_COUNT` states what the ceiling is not. Neither floor is moved, recomputed or
reinterpreted here: both are imported through `reptest`, which imports them from where they were fixed.

**Budget.** 0 model requests, no money, no downloads, no network. Every figure is read from disk or
computed from what is on disk, and a figure that cannot be established that way is named as such instead
of being estimated.
"""

from __future__ import annotations

from typing import Any

from genomeos.attribution import increase_links as il
from genomeos.attribution import reptest as rt

REGISTERED = "2026-10-02"
LANE = "lane-rescore"

WHAT_THIS_IS = (
    "a costing of the two arms that lane-reptest's no-go leaves open, in model requests, taken at zero "
    "requests. It is not a purchase, not an authorisation to purchase, not a direction, not an "
    "agreement rate and not a test. No request is sent by this lane under any outcome"
)

# ---- what is inherited, imported so that two lanes cannot come to mean different things -----------

INHERITED_RESULT = "data/results/reptest_answerable.json"
INHERITED_COMMITS = {
    "the answerability module, committed before the count": "365c87e",
    "the answerability runner": "6877f01",
    "the answerability result": "1734d95",
}
#: lane-reptest's own figures. None is recomputed as a finding of this lane: they are the denominators
#: this costing is taken against, and they are read from its result and asserted to match.
INHERITED_FIGURES = {
    "increase_derived_links": 48,
    "answerable_in_its_own_cell": 21,
    "independent_loci_of_the_answerable_links": 13,
    "the_gene_is_not_in_the_scorers_window_at_this_element": 18,
    "the_cache_carries_no_track_of_this_cell": 9,
    "the_element_is_absent_from_the_cache_archive": 0,
    "the_gene_was_scored_but_carries_no_entry_on_this_cells_track": 0,
    "the_cached_value_is_exactly_zero_so_it_carries_no_direction": 0,
}
INHERITED_NO_GO = rt.GATE_NO_GO

#: Both floors, imported through `reptest`, which imports them from where the project fixed them. This
#: module compares a CEILING against them and moves neither.
POSITIVE_FLOOR = rt.POSITIVE_FLOOR
LOCUS_FLOOR = rt.LOCUS_FLOOR
FLOORS_ARE_IMPORTED = (
    "both floors are imported through `reptest`, which imports them from "
    "`genomeos.attribution.fresh.POSITIVE_FLOOR` and from `genomeos.attribution.fresh.LOCUS_FLOOR`, "
    "which is `genomeos.attribution.cell2.POOLED_LOCUS_FLOOR`. This lane compares a ceiling against "
    "them and does not move, recompute, reinterpret or pool its way around either one. A costing that "
    "lowered a floor would be costing a different question"
)

# ---- what may be costed, and why the shortfall may not be --------------------------------------

FOLLOWS_A_NO_GO_WHOSE_OUTCOME_WAS_NEVER_READ = (
    "this costing follows a BLINDED FEASIBILITY NO-GO whose direction outcome was NEVER READ. "
    "lane-reptest stopped at the gate: no direction was taken, no agreement rate computed, no balanced "
    "accuracy formed and no interval drawn, on any arm, before or after the count. That is what makes a "
    "purchase of more data admissible rather than post-hoc - the scientific outcome does not exist to be "
    "chosen around. It is recorded here and not in a covering note, because a reader deciding whether to "
    "spend needs it beside the number"
)
THE_SHORTFALL_MAY_NOT_BE_BOUGHT = (
    "the shortfall is not the population. A set defined as 'the links that were missing' is defined BY "
    "the count that refused the test, and buying it would choose a population by a feasibility outcome "
    "even though no scientific outcome was read. So the population costed here is stated by a rule that "
    "can be written down without reference to the no-go at all, and the rule is applied to every link "
    "it admits - including links that are already answerable and links in the other arm, which cost "
    "requests and do not help the arm that was short. Selecting WHICH links to buy by anything "
    "outcome-related, or buying only the gap, is what is refused"
)
POPULATION_RULE_CELL_ARM = (
    "EVERY measured significant link the committed extractor emits under v2 whose own cell is not one "
    "of the cell tracks the finished sweep kept. The rule names no count, no arm, no floor and no "
    "outcome: it is a property of the cell, and it admits every link in such a cell whichever direction "
    "the measurement had. Both directions are costed and reported APART, because a purchase that mostly "
    "buys decreases is a different proposition from one that buys increases, and the arm that was short "
    "is the increase arm"
)
CELLS_NOT_COVERED_IS_COMPUTED_NOT_INHERITED = (
    "which cells the sweep does not cover is computed from the links themselves against "
    "`direction_link.CACHED_CELLS`, and NOT read from `direction_link.CELLS_NOT_CACHED`. That constant "
    "is the pair of cells the increase arm happened to carry ('WTC11', 'HCT116') and it does not name "
    "every uncovered cell the measured population holds: the decrease arm also holds Jurkat, which the "
    "sweep covers no better. Inheriting the constant would have costed a population defined by the arm "
    "that was short, which is the thing this costing may not do"
)
POPULATION_RULE_WINDOW_ARM = (
    "EVERY measured significant link the extractor emits under v2 whose tested gene's transcription "
    "start lies outside the window the sweep scored at that link's own element. This rule is also "
    "statable without reference to the no-go - it is a geometric property of the element and the "
    "benchmark's own TSS - so the window arm is NOT refused on admissibility. It is refused on a "
    "service limit, which is a harder refusal: see `WINDOW_ARM_IS_CLOSED_ON_A_SERVICE_LIMIT`"
)
WINDOW_ARM_IS_CLOSED_ON_A_SERVICE_LIMIT = (
    "the window arm cannot be bought at any price and for any population. The sweep already asks at "
    "1,048,576 bases, the largest of the four widths the client accepts, and a request outside that set "
    "is refused by the client before it is sent. So there is no wider window to buy for any link, "
    "however the population is defined, and the cost of this arm is not a number of requests: it is NO "
    "AVAILABLE REQUEST. The population is still counted, because the size of what cannot be bought is "
    "the useful figure"
)

# ---- the request unit, read from the committed chain ----------------------------------------------

REQUEST_UNIT = (
    "ONE MODEL REQUEST PER ELEMENT. One request deletes one element from its window and returns every "
    "gene in that window on every track the scorer's track table carries, so a request is per element "
    "and NOT per element per cell and NOT per element per gene. Which cell tracks are kept per gene is "
    "decided after the answer is in hand, by the `cells` argument of `aggregate`, so a cell the sweep "
    "did not keep costs that element's request again - one request, not one per cell. An element whose "
    "answer is already cached costs nothing, because `score_element` returns the cached answer without "
    "a request; an element whose cached answer lacks the cell wanted must be asked again, because the "
    "tracks the sweep did not keep were never written down"
)
REQUEST_UNIT_READ_FROM = {
    "one request per element": (
        "genomeos/predict/enhancer_target.py, module docstring: 'One request per element, about eight "
        "seconds'"
    ),
    "the chain that produced the sweep counts it the same way": (
        "scripts/enhancer_targets_all.py, module docstring: 'about 12,000 elements on chr21, one "
        "request each, cached per element'"
    ),
    "a cached element costs no request": (
        "genomeos/predict/enhancer_target.py, `score_element`: 'Cached per element id; a cached answer "
        "is returned without a request', and the body asks the scorer only when `load_cached` is None"
    ),
    "one request carries every gene and every track": (
        "genomeos/predict/alphagenome_adapter.py, `_live_scorer`: one `client.score_variant` call per "
        "scorer call, and the response is iterated over every gene and every track"
    ),
    "the cells are a post-processing filter and not a request dimension": (
        "genomeos/predict/enhancer_target.py, `aggregate(effects, cells=CELLS)`: the per-cell values "
        "are kept for the names in `cells` out of an answer that already holds all of them, and "
        "`CELLS` is the sweep's four"
    ),
    "the HCT116 precedent used this unit on this code path": (
        "scripts/crispri_hct116.py, module docstring: 'One deletion request per registry element "
        "overlapping a held-out HCT116 pair', with `KEPT_CELLS = (*CELLS, CELL)` - the cells widened, "
        "the unit unchanged"
    ),
}
REQUEST_UNIT_IS_CHECKED_NOT_CARRIED_OVER = (
    "the AstroREG registration counted 'one request per registry element overlapping a covered pair'. "
    "That is the same unit - one request per element - with a different denominator, because it was "
    "counting elements selected by overlap with a screen pair. It is not carried over as a figure and "
    "its waste fraction is not transferred: this lane counts its own elements, reports how many of its "
    "requests serve no link the test could score, and names the two precedents beside its own figure "
    "rather than inheriting either one's ratio"
)

# ---- the window, and why it is the binding constraint --------------------------------------------

#: The four sequence lengths the pinned client accepts, as `alphagenome/models/dna_client.py` defines
#: them in alphagenome 0.9.0 (the version uv.lock pins, wheel sha256
#: a4f35884341ae85b5d2cf088dfe0304961de7ae4a590d4538653069673de32f4). `score_variant` validates the
#: interval's width against this set and raises for a width outside it, so the set is a hard boundary
#: and not a default.
SUPPORTED_SEQUENCE_LENGTHS = (16_384, 131_072, 524_288, 1_048_576)
MAX_SUPPORTED_WINDOW = max(SUPPORTED_SEQUENCE_LENGTHS)
WINDOW_USED_BY_THE_SWEEP = 1_048_576
WINDOW_IS_AT_THE_MAXIMUM = (
    "the committed scorer already asks at the widest window the service accepts. "
    "`AlphaGenomeAdapter._live_scorer` resizes the variant's reference interval to "
    "`dna_client.SEQUENCE_LENGTH_1MB` (1,048,576 bases), and that is the largest of the four lengths "
    "`dna_client.SUPPORTED_SEQUENCE_LENGTHS` holds; `dna_client` checks the width against that set and "
    "raises for any other. So 'score it at a wider window' is not a thing that can be bought: for a "
    "link whose tested gene lies outside the 1 Mb window there is no request to send at any price, and "
    "the cost of that arm is not a large number of requests but NO AVAILABLE REQUEST"
)
WINDOW_CENTRING = (
    "`genome.Interval.resize` centres the new width on the interval's own centre, so the window the "
    "sweep scored at an element is [centre - width/2, centre + width/2) for width 1,048,576, where the "
    "centre is that of the variant's reference interval as `score_element` builds it: the element's "
    "bases with one anchor base before them, [start - 1, end). The window is computed here from the "
    "element's own scored interval as the cache recorded it, not from the row's interval, so the "
    "geometry is the run's and not a reconstruction of it"
)
TSS_IS_THE_BENCHMARK_S = (
    "the gene's position is the benchmark's own `startTSS` for THAT pair, read through "
    "`crispri.Pair.tss`, and the distance is the benchmark's own `distanceToTSS` read through "
    "`crispri.Pair.distance`, both from the table the pair came from. No TSS is looked up elsewhere, "
    "no annotation is joined in, and the benchmark's distance is reported beside the offset this lane "
    "computes so a reader can see they agree"
)
WINDOW_TEST_IS_NECESSARY_NOT_SUFFICIENT = (
    "a TSS inside the window is what the window CAN reach; it is not a guarantee that the service's "
    "gene scorer would emit a row for that gene, and a TSS outside the window is what the window "
    "CANNOT reach. The rule is therefore used in one direction only: outside the widest supported "
    "window means no request can answer it, which is a statement about the service's accepted "
    "interval and not about the gene. Inside the window is reported as 'the window reaches it' and "
    "never as 'a request would answer it'"
)

# ---- the purchasability causes, exhaustive and disjoint over the 27 ------------------------------

PURCHASABLE_CELL_TRACK = "purchasable: one request per element, with this cell kept"
UNPURCHASABLE_BEYOND_MAX_WINDOW = (
    "unpurchasable at any price: the gene lies outside the widest window the service accepts"
)
NOT_A_PURCHASE_NAMING = (
    "not a purchase: the gene is inside the window the sweep already scored, under another symbol"
)

PURCHASABILITY = (PURCHASABLE_CELL_TRACK, UNPURCHASABLE_BEYOND_MAX_WINDOW, NOT_A_PURCHASE_NAMING)
PURCHASABILITY_TEXT = {
    PURCHASABLE_CELL_TRACK: (
        "the sweep's record at this element already names the link's tested gene, and what is missing "
        "is only the link's own cell track, which the sweep did not keep. One request for this element "
        "with that cell among the kept cells would write a value for it. This is the HCT116 arm's own "
        "code path with the cells widened, which is how that arm was built"
    ),
    UNPURCHASABLE_BEYOND_MAX_WINDOW: (
        "the tested gene's transcription start lies further from the element than the widest window "
        "the client accepts reaches, so no request can place that gene in the scorer's window. There "
        "is nothing to buy: this is not an expensive link, it is an unavailable one"
    ),
    NOT_A_PURCHASE_NAMING: (
        "the tested gene's transcription start lies INSIDE the window the sweep already scored at this "
        "element, and the record names no gene by the benchmark's symbol. A request would not change "
        "that, because the window already reached the place. The record is reported with the symbols "
        "it does name, and whether the registered predicate may match a gene under a second symbol is "
        "a question about the predicate and not a cost. This lane does not amend the predicate, does "
        "not match an alias and does not count such a link as answerable"
    ),
}
PURCHASABILITY_IS_EXHAUSTIVE = (
    "every one of the 27 links lane-reptest could not answer falls in exactly one of these three, and "
    "the result asserts they sum to 27 and to the two arms' counts separately. The arms are reported "
    "apart, and a purchasable count is never quoted as an answerable count"
)

# ---- the ceiling, and what it is not -------------------------------------------------------------

CEILING_RULE = (
    "the CEILING of a purchase is the answerable count that would follow if every purchasable link's "
    "request returned a value that is not exactly zero on that link's own cell track: the links "
    "lane-reptest already answers, plus one per purchasable link, and the independent loci of that set "
    "by the same imported rule. It is an upper bound and the gate is taken on it, because a ceiling "
    "that misses a floor settles the question without any purchase: no outcome of the requests can "
    "reach the floor"
)
CEILING_IS_NOT_A_COUNT = (
    "the ceiling is not a count of answerable links and may not be quoted as one. It is what the "
    "answerable count could at most become: every purchasable request is assumed to return a value "
    "that is not exactly zero on the right track, which is not established for any of them and cannot "
    "be established without sending the request. lane-reptest's own count of links excluded for an "
    "exactly-zero value is 0 of 48, which is a fact about the links it could read and NOT a rate, a "
    "probability or a prediction about a link it could not"
)
CEILING_MISSES_A_FLOOR = (
    "the ceiling of the purchase is below an imported floor: buying every request that CAN be bought "
    "would leave the direction test still unanswerable. The purchase would buy a larger number that "
    "still cannot be tested, so the costing's answer is that the money does not buy the test. This is "
    "the headline of the costing and not a caveat on it. The margin is reported because a reader is "
    "owed it and NOT because a small margin is better than a large one; a ceiling short of a floor is "
    "reported as a no-go and NEVER as close, promising, nearly enough, a good start or enough for a "
    "pilot"
)
CEILING_CLEARS_BOTH_FLOORS = (
    "the ceiling of the purchase is at or above both imported floors: the purchase COULD make the test "
    "answerable, and whether it does depends on what the requests return. Nothing here authorises the "
    "purchase, registers a test or reads a direction"
)
CEILING_READINGS = (CEILING_CLEARS_BOTH_FLOORS, CEILING_MISSES_A_FLOOR)
FORBIDDEN_OF_A_SHORT_CEILING = rt.FORBIDDEN_OF_A_SHORT_COUNT

#: A second, deliberately more generous ceiling, so that the obvious follow-up question cannot be put
#: later as though it had not been answered: suppose the purchase is made AND the naming gap is closed
#: for nothing, by a predicate that matched a gene under a current symbol as well as the benchmark's.
#: It assumes something this lane does not do and has no authority to do, and it is reported as an upper
#: bound on an upper bound.
MOST_GENEROUS_CEILING_RULE = (
    "the MOST GENEROUS ceiling: every purchasable request returns a value that is not exactly zero on "
    "the right track, AND the naming gap is closed at no cost by a predicate that matched the tested "
    "gene under a current symbol as well as the benchmark's. Both halves are assumptions. The second "
    "would be an amendment to a registered predicate, which this lane does not make, does not recommend "
    "here and has no authority to make; it is added so that 'but what if you also fixed the symbols' is "
    "answered in the same result rather than later"
)

# ---- what the project has already spent, for comparison ------------------------------------------

#: The two precedents a reader can judge a figure against, each read from the result that records it
#: rather than quoted from memory. The paths are declared inputs of this lane's run.
PRECEDENT_SOURCES = {
    "HCT116 second cell type, delivered": "data/results/crispri_published_v2.json, "
    "second_cell_type_hct116: `requests.delivered_in_the_carried_run.sent` and "
    "`as_carried.covered_pairs`",
    "WTC11 second cell type, costed and not bought": "data/results/crispri_published_v2.json, "
    "second_cell_type/cost/WTC11/requests_frozen_feature",
    "AstroREG request plan, authorised and not sent": "data/results/astroreg_request_plan.json, "
    "`summary.requests`, `summary.requests_serving_at_least_one_scored_pair` and "
    "`summary.requests_serving_no_scored_pair`",
}
PRECEDENTS_ARE_READ_NOT_QUOTED = (
    "each precedent figure is read out of the result that records it, in the run that writes this "
    "costing, and the paths are declared inputs. None is typed in from a brief or from memory, because "
    "a comparison a reader cannot check is not a comparison"
)

R2_WORDING = il.R2_WORDING
PROHIBITION_IS_ENFORCED_BY = (
    "genomeos.attribution.increase_links.check_no_mechanism_claim, called on every link record this "
    "lane handles. It is not weakened, not reimplemented and not bypassed, and this lane adds no word "
    "list of its own"
)
NO_DIRECTION_IS_READ = (
    "no direction, agreement rate, interval, control or test is computed here under any outcome. "
    "Whether a value exists and whether it is exactly zero is the answerability predicate's own "
    "question and is lane-reptest's; this lane reads neither a sign nor a value into any finding, and "
    "where it reports that a delivered value is not exactly zero it reports that fact and not the "
    "value, not the sign and not a comparison"
)


# ---- the functions the result is built from ------------------------------------------------------


def window_centre(scored_start: int, scored_end: int) -> int:
    """The centre every window at this element is taken about. WINDOW_CENTRING."""
    return (scored_start - 1 + scored_end) // 2


def window_bounds(scored_start: int, scored_end: int, width: int) -> tuple[int, int]:
    """The window a request at `width` would cover at an element, as the committed chain centres it.

    `scored_start` and `scored_end` are the element's own interval as the cache recorded it, and the
    variant's reference interval is [scored_start - 1, scored_end): the anchor base and the element's
    bases, which is what `score_element` fetches. `genome.Interval.resize` then centres `width` on that
    interval's centre, start = centre - (width + 1) // 2 and end = centre + width // 2.
    """
    centre = window_centre(scored_start, scored_end)
    return centre - (width + 1) // 2, centre + width // 2


def window_needed(scored_start: int, scored_end: int, tss: int) -> int:
    """The smallest even window width, centred as the chain centres it, that would contain `tss`.

    An even width w covers [centre - w // 2, centre + w // 2), so an offset of o = tss - centre is
    inside when -w // 2 <= o < w // 2: a non-negative offset needs w >= 2 * o + 2 and a negative one
    needs w >= 2 * |o|. The figure is reported per link and is not a request: for most of these links it
    is larger than any width the service accepts, which is the point.
    """
    offset = tss - window_centre(scored_start, scored_end)
    return 2 * offset + 2 if offset >= 0 else 2 * -offset


def smallest_supported_window(width_needed: int) -> int | None:
    """The smallest accepted sequence length at least `width_needed`, or None when none is.

    None is the costing's answer of "unpurchasable at any price": the client refuses any width outside
    `SUPPORTED_SEQUENCE_LENGTHS`, so there is no request to send. WINDOW_IS_AT_THE_MAXIMUM.
    """
    return next((w for w in SUPPORTED_SEQUENCE_LENGTHS if w >= width_needed), None)


def purchasability(gene_named_in_the_record: bool, width_needed: int, cell_is_cached: bool) -> str:
    """One link's purchasability, out of `PURCHASABILITY`, exhaustive and disjoint.

    The order is what makes the three disjoint. A gene the widest accepted window cannot reach is
    unpurchasable whatever else is true of it, because no request can be sent. A gene the sweep's record
    already names, missing only its cell track, is the one purchasable shape: one request per element
    with that cell kept. A gene inside the window that the record does not name is not a purchase at
    all, because the window already reached it.
    """
    if smallest_supported_window(width_needed) is None:
        return UNPURCHASABLE_BEYOND_MAX_WINDOW
    if gene_named_in_the_record and not cell_is_cached:
        return PURCHASABLE_CELL_TRACK
    return NOT_A_PURCHASE_NAMING


def cells_not_covered(cells: list[str]) -> list[str]:
    """The cells among `cells` the sweep kept no track of, computed and not inherited.

    CELLS_NOT_COVERED_IS_COMPUTED_NOT_INHERITED. `direction_link.CELLS_NOT_CACHED` names the two cells
    the increase arm carried and misses Jurkat, which the decrease arm carries and the sweep covers no
    better, so the set is taken against `CACHED_CELLS` from the links in hand.
    """
    return sorted({c for c in cells if c not in rt.CACHED_CELLS})


def requests_for(elements: list[str]) -> int:
    """The request count for a set of elements under `REQUEST_UNIT`: one per DISTINCT element.

    Two links on one element share one request, so the count is over the distinct ids and never over
    the links. This is the AstroREG registration's unit with this lane's own denominator, and the
    result reports both the element count and the link count so neither can be read as the other.
    """
    return len(set(elements))


def ceiling(answerable_now: int, purchasable_links: int, loci_at_the_ceiling: int) -> dict[str, Any]:
    """The gate on the CEILING of a purchase against both imported floors, with each margin.

    `answerable_now` and `loci_at_the_ceiling` are counted elsewhere by the imported rule; this takes
    the verdict and says, in the registered words, whether the money would buy a testable number.
    """
    links = answerable_now + purchasable_links
    short: list[str] = []
    if links < POSITIVE_FLOOR:
        short.append(
            f"links at the ceiling {links}, short of the imported floor {POSITIVE_FLOOR} by "
            f"{POSITIVE_FLOOR - links}"
        )
    if loci_at_the_ceiling < LOCUS_FLOOR:
        short.append(
            f"independent loci at the ceiling {loci_at_the_ceiling}, short of the imported floor "
            f"{LOCUS_FLOOR} by {LOCUS_FLOOR - loci_at_the_ceiling}"
        )
    met = not short
    return {
        "answerable_links_now": answerable_now,
        "purchasable_links": purchasable_links,
        "links_at_the_ceiling": links,
        "independent_loci_at_the_ceiling": loci_at_the_ceiling,
        "floors": {
            "links_at_least": POSITIVE_FLOOR,
            "independent_loci_at_least": LOCUS_FLOOR,
            "imported": FLOORS_ARE_IMPORTED,
        },
        "ceiling_rule": CEILING_RULE,
        "ceiling_is_not_a_count": CEILING_IS_NOT_A_COUNT,
        "meets_both_floors_at_the_ceiling": met,
        "short_by": short,
        "reading": CEILING_CLEARS_BOTH_FLOORS if met else CEILING_MISSES_A_FLOOR,
        "what_follows": (
            "the purchase could make the test answerable; nothing here authorises it"
            if met
            else "the purchase would not make the test answerable, so there is nothing to authorise"
        ),
    }


def per_request_figures(
    requests: int,
    answerable_positives_increase: int,
    answerable_positives_decrease: int,
    loci_added: int,
    requests_serving_no_scorable_link: int,
) -> dict[str, Any]:
    """REQUESTS PER ANSWERABLE POSITIVE, the figure the project costs purchases in, both arms apart.

    The denominator is the number of measured positives the purchase could at most make ANSWERABLE, and
    it is given three ways because the three are different propositions: the increase arm alone, which
    is the arm whose floors gated; the decrease arm alone, which costs requests and moves no floor of
    that arm; and the two together, which is the cheapest-looking ratio and the one a reader is most
    likely to be shown. Every one of them is over a CEILING, so every one is a best case.
    """
    both = answerable_positives_increase + answerable_positives_decrease

    def ratio(n: int) -> float | None:
        return round(requests / n, 4) if n else None

    return {
        "requests": requests,
        "answerable_positives_at_the_ceiling": {
            "increase_arm": answerable_positives_increase,
            "decrease_arm": answerable_positives_decrease,
            "both_arms": both,
        },
        "requests_per_answerable_positive_at_the_ceiling": {
            "increase_arm_only": ratio(answerable_positives_increase),
            "decrease_arm_only": ratio(answerable_positives_decrease),
            "both_arms_pooled": ratio(both),
            "which_one_bears_on_the_decision": (
                "the increase arm's. The floors that refused the test are the increase arm's answerable "
                "links and their loci, so a request that makes a decrease link answerable costs the "
                "same and moves neither of them. The pooled ratio is the flattering one and is given so "
                "that it cannot be quoted without the arm-by-arm figures beside it"
            ),
        },
        "independent_loci_at_most_added_to_the_increase_arm": loci_added,
        "requests_per_locus_at_the_ceiling": round(requests / loci_added, 4) if loci_added else None,
        "requests_serving_no_link_the_test_could_score": requests_serving_no_scorable_link,
        "what_these_are_not": (
            "a yield. Every denominator is a CEILING, so each ratio is the best case per request and "
            "not an observed cost per answerable positive. A request that returns a value of exactly "
            "zero on the link's own cell track buys nothing, and whether any of these would is not "
            "established and cannot be without sending it"
        ),
    }


def registration() -> dict[str, Any]:
    """Everything fixed before the costing, as the result states it. Nothing is computed here."""
    return {
        "lane": LANE,
        "registered": REGISTERED,
        "what_this_is": WHAT_THIS_IS,
        "acts_on": {
            "lane": "lane-reptest",
            "result": INHERITED_RESULT,
            "commits": dict(INHERITED_COMMITS),
            "figures": dict(INHERITED_FIGURES),
            "no_go_carried_verbatim": INHERITED_NO_GO,
        },
        "follows_a_no_go_whose_outcome_was_never_read": FOLLOWS_A_NO_GO_WHOSE_OUTCOME_WAS_NEVER_READ,
        "the_shortfall_may_not_be_bought": THE_SHORTFALL_MAY_NOT_BE_BOUGHT,
        "population_rule_cell_arm": POPULATION_RULE_CELL_ARM,
        "cells_not_covered_is_computed_not_inherited": CELLS_NOT_COVERED_IS_COMPUTED_NOT_INHERITED,
        "population_rule_window_arm": POPULATION_RULE_WINDOW_ARM,
        "window_arm_is_closed_on_a_service_limit": WINDOW_ARM_IS_CLOSED_ON_A_SERVICE_LIMIT,
        "the_registration_that_would_spend_this_is_a_new_one": (
            "nothing here registers a purchase or a test. A purchase would need its own registration, "
            "with its population rule, its floors and its readings fixed BEFORE any request is sent, "
            "and this lane does not write it. This lane costs it"
        ),
        "request_unit": REQUEST_UNIT,
        "request_unit_read_from": dict(REQUEST_UNIT_READ_FROM),
        "request_unit_is_checked_not_carried_over": REQUEST_UNIT_IS_CHECKED_NOT_CARRIED_OVER,
        "window_used_by_the_sweep": WINDOW_USED_BY_THE_SWEEP,
        "supported_sequence_lengths": list(SUPPORTED_SEQUENCE_LENGTHS),
        "max_supported_window": MAX_SUPPORTED_WINDOW,
        "window_is_at_the_maximum": WINDOW_IS_AT_THE_MAXIMUM,
        "window_centring": WINDOW_CENTRING,
        "window_test_is_necessary_not_sufficient": WINDOW_TEST_IS_NECESSARY_NOT_SUFFICIENT,
        "tss_is_the_benchmarks": TSS_IS_THE_BENCHMARK_S,
        "purchasability": {p: PURCHASABILITY_TEXT[p] for p in PURCHASABILITY},
        "purchasability_is_exhaustive": PURCHASABILITY_IS_EXHAUSTIVE,
        "ceiling_rule": CEILING_RULE,
        "most_generous_ceiling_rule": MOST_GENEROUS_CEILING_RULE,
        "ceiling_is_not_a_count": CEILING_IS_NOT_A_COUNT,
        "ceiling_readings": list(CEILING_READINGS),
        "there_is_no_third": True,
        "forbidden_of_a_short_ceiling": list(FORBIDDEN_OF_A_SHORT_CEILING),
        "floors_are_imported": FLOORS_ARE_IMPORTED,
        "independent_locus_rule": rt.LOCUS_RULE,
        "not_biological_independence": rt.NOT_BIOLOGICAL_INDEPENDENCE,
        "precedent_sources": dict(PRECEDENT_SOURCES),
        "precedents_are_read_not_quoted": PRECEDENTS_ARE_READ_NOT_QUOTED,
        "binding_wording_r2": R2_WORDING,
        "prohibition_is_enforced_by": PROHIBITION_IS_ENFORCED_BY,
        "no_direction_is_read": NO_DIRECTION_IS_READ,
        "count_is_not_coverage": rt.COUNT_IS_NOT_COVERAGE,
        "budget": (
            "0 model requests, no money, no downloads, no network. Every figure is read from disk or "
            "computed from what is on disk. A figure that cannot be established that way is NAMED as "
            "one that cannot, and never estimated: this lane counts what a purchase would cost and "
            "makes none"
        ),
        "order": [
            "1. the extractor's own v2 links, both directions, re-enumerated genome-wide and blind to "
            "every cached value",
            "2. the two populations by their own rules: the cells the sweep does not cover, and the "
            "links whose tested gene lies outside the window the sweep scored",
            "3. each link's tested gene's TSS, from the benchmark pairs behind that link",
            "4. the window each would need, against the widest width the service accepts",
            "5. the purchasability of each, exhaustive and disjoint, and the request count by element",
            "6. the ceiling of the purchase against both imported floors, the arms apart",
            "7. the precedents, read from the results that record them",
        ],
    }
