# SPDX-License-Identifier: AGPL-3.0-or-later
"""Sign agreement at the link level: does the predicted deletion direction match the measured sign
more often for decreases than for increases, and above a sign-shuffled control. The registration.

**What this follows from.** lane-increase counted the significant increases of the CRISPRi benchmark
along a ladder registered before the count and gated the population it could build against two
imported floors (`data/results/increase_population.json`, `data/results/increase_registration.json`).
Its reading, carried word for word and never strengthened here, is
`increases.GO`: *"at or above both registered floors: a measured repression population could be built
from the increases the committed extractor discards..."* - **48** links over **33** independent loci,
against floors of 30 links and 20 loci. Behind it stands lane-repress2's registered no-go, also carried
word for word: *"the measured layer holds no repression call to test"*.

**Why this test is at the link level and not at the rule level.** lane-increase's ladder is explicit
about what the compiled program holds for these increases: of the 48, exactly **one** has any compiled
rule naming its own gene, and **zero** have one gated on its own cell
(`increase_population.json`, steps `and_a_rule_on_that_element_names_the_pairs_own_gene` = 1 increase
and `and_a_rule_on_that_element_and_gene_is_gated_on_the_pairs_own_cell` = 0 increases). A rule-level
test of direction on the increases is therefore impossible, and this module does not attempt one. The
model's direction calls are read **directly**, out of the cached deletion answers, for the link's own
gene on the link's own cell track - never through a compiled rule.

**The statistic, and why it is not the difference between the arms.** The question this lane was given
is comparative: does the predicted direction match the measured sign *more often for decreases than for
increases*, above a sign-shuffled control. Writing the control's own unit test, before any value was
read, showed that this difference cannot carry that question. With `n_d` answered decrease links of
which `a` agree and `n_u` answered increase links of which `b` agree, the excess of the difference over
the shuffled baseline is identically `2 * (balanced_accuracy - 0.5) * (n_u - n_d) / (n_d + n_u)`. When
the arms are the same size it is exactly zero, so the whole difference is the shuffled baseline; and
when the decreases arm is the larger one - as it is here, 212 against 48 - a *positive* excess requires
balanced accuracy **below** 0.5, so a one-sided test of the difference against the shuffle is a test of
the model doing worse than chance. The same algebra shows the two arms do not carry two skills: each
arm's agreement against its own shuffled chance level differs from the other's only by that arm-size
factor. With one marginal sign rate there is exactly **one** skill number on this population, and it is
balanced accuracy. So both arms' rates and the difference between them are reported, apart, with their
intervals, the baseline and the excess, as a description of the asymmetry; and the reading is taken on
balanced accuracy against 0.5, which a constant-sign caller of either sign scores exactly. `COLLAPSE`
holds this in the registration rather than in an explanation written afterwards.

**What this module is.** The registration, fixed and committed before any direction is read: the
eligibility predicate; how a predicted value becomes a sign and what happens to a zero and to an
absent value; the sign-shuffled control, its count and what each shuffle holds fixed; the locus
convention, imported from `cell2` with its own wording that it is an operational grouping and not
established biological independence; both floors, imported from `fresh`, which is where the project
fixed them before this lane existed; how the intervals are taken over **independent loci** and not
over links; and the two readings, with no third. The computing code lives in
`scripts/direction_link.py` and runs after this registration is committed.

**The wording that is binding (review item R2).** A measured increase on knockdown is *an increase on
knockdown*. It is never a silencer, never a repressor and never evidence of a repression mechanism.
This test is about **sign agreement**, not mechanism. Nothing here reads a mechanism from a sign, and
no reading of lane-increase's or lane-repress2's is restated in stronger words than they registered.

**What is not done.** No extractor, compiled rule, threshold, floor or registered reading is changed.
No model request is made, nothing is downloaded, no money is spent. No cell is substituted for another:
where the cached answers do not carry the link's own cell the link is unanswerable and is counted as
unanswerable, and the arm it would have joined is reported without it.
"""

from __future__ import annotations

import random
from typing import Any

from genomeos.attribution import cell2, fresh
from genomeos.attribution import increases as inc
from genomeos.attribution import measured as ms
from genomeos.attribution import repress2 as rp

# ---- what this lane inherits, carried word for word with its file --------------------------------

#: lane-repress2's registered no-go, carried through lane-increase's own constant, never restated.
INHERITED_REPRESS2_NO_GO = inc.INHERITED_NO_GO
INHERITED_REPRESS2_FILE = inc.INHERITED_FROM

#: lane-increase's registered GO reading of its P1 gate, imported rather than retyped.
INHERITED_INCREASE_GO = inc.GO
INHERITED_INCREASE_FILE = "data/results/increase_population.json"
INHERITED_INCREASE_REGISTRATION = "data/results/increase_registration.json"
INHERITED_INCREASE_COMMITS = {
    "registration": "e1c2671",
    "population": "f521bdb",
    "module": "6ca8095",
    "attribution_section": "7f2b354",
}

#: Every figure of lane-increase's that constrains this lane, with the file it is read from. None of
#: them is recomputed here and none may be quoted under a noun its own breakdown does not support.
INHERITED_FIGURES = {
    "file": INHERITED_INCREASE_FILE,
    "significant_increases_in_the_benchmark": 159,
    "significant_increases_on_an_attributed_element": 48,
    "significant_decreases_on_an_attributed_element": 212,
    "increase_links": 48,
    "increase_links_independent_loci": 33,
    "increase_links_by_cell": {"K562": 39, "WTC11": 6, "HCT116": 3},
    "increases_behind_the_links_by_split": {"training": 25, "heldout": 23},
    "increases_with_a_rule_naming_their_own_gene": 1,
    "increases_with_a_rule_gated_on_their_own_cell": 0,
    "contested_by_a_regulated_pair": 0,
    "largest_increase_per_link_over_the_40_links_that_file_shows": {
        "min": 0.0102,
        "median": 0.0987,
        "max": 4.2101,
        "at_or_above_0.10": 20,
        "of": 40,
    },
}

#: The population is not balanced across cell types, and this is stated before any rate is read.
NOT_BALANCED = (
    "the increase population is not balanced across cell types: 39 of the 48 links are K562, 6 WTC11 "
    "and 3 HCT116 (data/results/increase_population.json). No rate of this lane's is a genome-wide or "
    "a cell-type-general rate, and none is reported pooled across cells"
)

#: The split rule. The 2026-09-28 split audit permits a held-out row as evaluation only.
SPLIT_RULE = (
    "23 of the 48 increases are held-out rows (data/results/increase_population.json, "
    "increases_behind_the_links_by_split). The 2026-09-28 split audit permits a held-out row as "
    "EVALUATION only, never as a feature and never as a fit. This lane fits nothing, selects nothing "
    "on an outcome and builds no feature: it reads a cached predicted value and a cached measured "
    "outcome and compares their signs, so every eligible link is evaluation and both splits are "
    "admitted. The training/held-out breakdown is reported beside every count"
)

#: The floors are on COUNTS, not on effect sizes, and the distinction is recorded before the gate.
FLOORS_ARE_ON_COUNTS = (
    "both imported floors are on COUNTS - links and independent loci - and neither is a floor on an "
    "effect size. Over the 40 links data/results/increase_population.json shows, the largest increase "
    "per link runs 0.0102 to 4.2101 with a median of 0.0987 and 20 of 40 at or above 0.10. No "
    "magnitude bar is imposed here, no link is dropped for being small, and clearing a count floor is "
    "not a statement that the effects behind it are large"
)

#: The confusion this lane's decrease arm must not inherit, recorded before the arm is built.
DECREASE_ARM_IS_NOT_MODEL_COVERAGE = (
    "the 212-of-212 that the decrease column keeps at lane-increase's last two ladder steps is the "
    "EXTRACTOR'S OWN CONSTRUCTION and not model coverage: a significant decrease on an attributed "
    "element raises the measured rule that then names its gene and its cell, so the decreases survive "
    "those two steps by definition. Reading it as coverage would repeat the error docs/LESSONS.md "
    "records under 'A count is not a measurement of the thing you want to count'. This lane's decrease "
    "arm therefore does NOT sit at those steps: it is taken at the same step as the increase arm, "
    "`and_on_an_attributed_element`, where no compiled rule is consulted at all, so the two arms share "
    "one eligibility rule and the comparison between them is not a comparison of two different gates"
)

#: cell2's convention, imported with its own wording, including the sentence that it is not biological
#: independence. It is carried into every count, interval and summary this lane writes.
LOCUS_RULE = cell2.INDEPENDENT_LOCUS_RULE
LOCUS_SPAN = cell2.INDEPENDENT_LOCUS_SPAN
LOCUS_RULE_IMPORTED_FROM = "genomeos.attribution.cell2.INDEPENDENT_LOCUS_RULE"
NOT_BIOLOGICAL_INDEPENDENCE = inc.NOT_BIOLOGICAL_INDEPENDENCE
INHERITED_33_LOCI_ARE_OPERATIONAL = (
    "lane-increase's 33 independent loci are cell2's operational grouping and not established "
    "biological independence; so is every locus count and every interval this lane reports"
)

#: Both floors, imported from where the project fixed them. Neither is chosen here and neither may
#: move after an outcome is seen.
POSITIVE_FLOOR = fresh.POSITIVE_FLOOR
LOCUS_FLOOR = fresh.LOCUS_FLOOR
FLOORS = {
    "links": POSITIVE_FLOOR,
    "links_imported_from": "genomeos.attribution.fresh.POSITIVE_FLOOR",
    "independent_loci": LOCUS_FLOOR,
    "independent_loci_imported_from": (
        "genomeos.attribution.fresh.LOCUS_FLOOR, which is genomeos.attribution.cell2.POOLED_LOCUS_FLOOR"
    ),
    "neither_chosen_here": True,
    "applied": "to each arm separately, never to the two arms pooled",
    "on_counts_not_effect_sizes": FLOORS_ARE_ON_COUNTS,
}

# ---- eligibility -----------------------------------------------------------------------------

DECREASES = "decreases"
INCREASES = "increases"
ARMS = (DECREASES, INCREASES)

ELIGIBILITY = (
    "a link is one (attributed element, gene, cell) for which that element's reciprocally overlapping "
    "cached CRISPRi pairs hold at least one pair whose outcome is exactly one of "
    f"{(ms.DECREASE, ms.INCREASE)!r}, under the committed join: the element and the pair's tested "
    f"interval reciprocally overlap at measured.RECIPROCAL_OVERLAP = {ms.RECIPROCAL_OVERLAP} or more, "
    f"computed by measured.reciprocal_overlap over the candidates measured.Layer.near returns at "
    f"measured.REACH = {ms.REACH}. That is lane-increase's `and_on_an_attributed_element` step, the "
    "same step for both arms, and no compiled rule enters it. The outcome is read off the cached row's "
    "own outcome field and is never recomputed. A link is in the decreases arm when its pairs hold a "
    "significant decrease and no significant increase, and in the increases arm when they hold a "
    "significant increase and no significant decrease; a link whose pairs hold both is CONTESTED by "
    "the measurement itself, is in neither arm, and is counted and reported on its own. This is "
    "lane-increase's own contest rule applied symmetrically to both arms, and that lane recorded 0 "
    "contested links on the increase side (data/results/increase_population.json)"
)

CONTEST_RULE = (
    "a (element, gene, cell) carrying both a significant decrease and a significant increase is in "
    "NEITHER arm. Its measured sign is contested by the measurement, so an agreement computed on it "
    "would be an agreement with a choice this lane made rather than with a measurement. Contested "
    "links are counted and reported and are in no arm's denominator"
)

#: The measured sign of a link, which follows from the arm it is in and is never read from a magnitude.
MEASURED_SIGN = {DECREASES: -1, INCREASES: +1}
MEASURED_SIGN_CALL = (
    "the measured sign of a link is -1 in the decreases arm and +1 in the increases arm, by the "
    "arm's own definition. No effect size, magnitude or threshold enters it: the sign is the "
    "benchmark's own significance label, and a link is in exactly one arm or is contested"
)

# ---- the predicted direction, and what the cached answers carry -----------------------------------

CACHE = "data/knowledge/alphagenome/elements/<chrom>.json.gz"
PREDICTED_CALL = (
    "the predicted deletion direction of a link is the sign of the signed predicted log2 fold change "
    "the finished genome-wide deletion sweep cached for THAT link's own gene on THAT link's own cell "
    f"track, read from the per-element response cache {CACHE} under the link's own element id, at "
    "`genes[*].by_cell[cell]` for the gene whose name is exactly the link's gene. The link already "
    "names one element, so no second element is consulted and no resolution rule across overlapping "
    "elements is needed or applied. Nothing is recomputed, no model request is made and no value is "
    "derived from any other cell, gene or element"
)

#: Exactly the cell tracks the cached answers carry. Verified from the cache before this registration
#: was written, by reading the first record of the smallest archive only.
CACHED_CELLS = ("HepG2", "IMR-90", "K562", "GM12878")
CELLS_NOT_CACHED = ("WTC11", "HCT116")
CELL_AVAILABILITY = (
    f"the cached deletion answers carry exactly these cell tracks: {CACHED_CELLS!r}. They do NOT carry "
    f"{CELLS_NOT_CACHED!r}. Of lane-increase's 48 increase links, 39 are K562 and are answerable in "
    "their own cell; 6 are WTC11 and 3 are HCT116 and are NOT, because the sweep holds no track of "
    "either cell. Those 9 links are reported as unanswerable with their counts and are in no "
    "denominator. NO SUBSTITUTION IS MADE: no other cell's value, no other gene's value and no pooled "
    "or marginal value stands in for a cell the cache does not carry, because a substituted cell would "
    "make the comparison between the arms meaningless. A separate partial cache "
    "data/knowledge/alphagenome/elements_hct116 exists over five chromosomes from a different run and "
    "does carry an HCT116 track; it is NOT read here, because mixing a second cache into one arm of a "
    "comparison is the same substitution by another route. It is named so that a later lane can "
    "register a test that uses it"
)
PRIMARY_CELL = "K562"
PRIMARY_CELL_CALL = (
    "the test is taken on the K562 links of both arms and on no other cell, because K562 is the one "
    "cell of this population the cached answers carry. Both arms are restricted identically, so the "
    "comparison is within one cell and no cell-type difference can be mistaken for an arm difference. "
    "The counts of every other cell are reported beside it as unanswerable and are never pooled in"
)

# ---- how a value becomes a sign ------------------------------------------------------------------

SIGN_CALL = (
    "a predicted value v becomes a sign: -1 when v < 0, +1 when v > 0, and NO SIGN when v == 0.0 "
    "exactly. A link whose predicted value is exactly zero carries no predicted direction, is "
    "EXCLUDED from its arm's agreement denominator, and is counted and reported as "
    "`predicted_zero_excluded`. A link whose element is absent from the cache, or whose gene carries "
    "no entry for the cell track, is ABSENT: it is excluded from the denominator and counted and "
    "reported as `absent`. Neither a zero nor an absent value is ever scored as an agreement or as a "
    "disagreement, and neither is replaced by a default, a nearest gene, another cell or a marginal "
    "sign. Every denominator is reported with its own exhaustive breakdown of answered, zero and "
    "absent, so no count of this lane's can be quoted under a noun the breakdown does not support "
    "(docs/LESSONS.md, 'A count is not a measurement of the thing you want to count')"
)

AGREEMENT_CALL = (
    "an answered link AGREES when the sign of its predicted value equals its arm's measured sign, and "
    "DISAGREES otherwise. An arm's agreement rate is agreements divided by answered links of that "
    "arm. The two arms are reported APART and are never pooled: a pooled agreement rate over a "
    "lopsided mix is the majority arm's own rate wearing a general noun, and it is not reported as a "
    "result here. The primary statistic is the DIFFERENCE, decreases minus increases"
)

PRIMARY_STATISTIC = "agreement rate in the decreases arm minus agreement rate in the increases arm"

#: The baseline the sign shuffle implies, and the trap it exists to close. Written before any value was
#: read, after the control's own unit test showed that a difference of sensitivities has a permutation
#: null at this baseline and NOT at zero.
MARGINAL_CALL = (
    "the model's own marginal sign rate on the answered links, q = the share of answered links whose "
    "predicted sign is -1, is reported beside every rate, because the difference between the arms has "
    "a baseline built out of it. Under the sign shuffle the decreases arm scores q and the increases "
    "arm scores 1 - q, so the shuffled difference is centred on BASELINE = 2q - 1 and NOT on zero. A "
    "difference compared against zero would therefore manufacture a success out of nothing but a model "
    "that calls one sign more often than the other, which is the trap crispri_direction.py registered "
    "in its own words: 'a sign agreement looks like a coin flip, so 0.5 looks like its chance level. It "
    "is not, and reporting it against 0.5 would manufacture a success'"
)
BASELINE = "2q - 1, where q is the share of answered links whose predicted sign is -1"
EXCESS_CALL = (
    "EXCESS = the difference between the arms minus the baseline the sign shuffle implies, "
    "(decreases rate - increases rate) - (2q - 1). The sign shuffle leaves every predicted sign where "
    "it is, so q is identical in every shuffle and the excess has a permutation expectation of exactly "
    "zero. It is reported with the difference and the baseline, and never one of the three without the "
    "other two"
)

#: The structural fact this lane found from the control's own unit test, BEFORE any value was read, and
#: the reason the comparison between the arms is not the statistic the reading is taken on. It is
#: recorded here as the registration rather than discovered afterwards as an excuse.
COLLAPSE = (
    "the difference between the arms, held against a sign shuffle, cannot answer the question it looks "
    "like it answers, and the algebra is exact: with n_d answered decrease links of which a agree and "
    "n_u answered increase links of which b agree, the excess of the difference over the shuffled "
    "baseline is identically EXCESS = 2 * (balanced_accuracy - 0.5) * (n_u - n_d) / (n_d + n_u). Two "
    "consequences follow and both were written down before the run. First, when the two arms are the "
    "same size the excess is exactly zero, so the whole difference between the arms IS the shuffled "
    "baseline and carries no information beyond the model's own marginal sign rate. Second, when the "
    "decreases arm is the larger one - which it is on this population, 212 against 48 - the factor "
    "(n_u - n_d) is negative, so a positive excess requires balanced accuracy BELOW 0.5: a one-sided "
    "test of the difference against the sign shuffle is a test of the model doing WORSE than chance, "
    "and a model that reads direction well scores below that control rather than above it. The same "
    "algebra says the two arms do not carry two skills: each arm's agreement against its own shuffled "
    "chance level (q for the decreases, 1 - q for the increases) differs from the other's only by that "
    "arm-size factor, and the two sum to 2 * (balanced_accuracy - 0.5). With one marginal sign rate "
    "there is exactly ONE skill number on this population and it is balanced accuracy. So the "
    "difference between the arms is reported, with its interval, its baseline and its excess, as a "
    "DESCRIPTION of the asymmetry and never as evidence that the model reads direction better on one "
    "sign; and the statistic the reading is taken on is balanced accuracy"
)

#: The primary statistic, and the project's own registered reason for it, carried from
#: crispri_direction_both.py word for word.
CHANCE = 0.5
BALANCED_CALL = (
    "the statistic the reading is taken on is BALANCED ACCURACY, the unweighted mean of the two arms' "
    "agreement rates, held against 0.5. crispri_direction_both.py registered the reason in these words "
    "and they are carried here unchanged: 'the registered primary statistic for the combined set is "
    "balanced accuracy, the unweighted mean of the two per-sign sensitivities. Its chance level is 0.5 "
    "for any class mix, and a constant-sign caller - of either sign - scores exactly 0.5 on it by "
    "construction.' And, from the same lane, on why a raw rate is not reported alone: 'the "
    "constant-sign caller's score on this set. Raw agreement is compared to THIS, not to 0.5; balanced "
    "accuracy is the statistic that compares to 0.5 honestly.' Balanced accuracy is a function of BOTH "
    "arms and of neither alone, the two arms' own rates and denominators are always reported apart "
    "beside it, and no pooled raw agreement rate is reported as a result"
)
TEST_STATISTIC = "balanced accuracy, the unweighted mean of the two arms' agreement rates, against 0.5"

#: What the project's two earlier direction lanes registered, carried word for word with their file.
#: Their population is different from this one - held-out in-reach pairs, K562 and GM12878, not links on
#: an attributed element - so nothing of theirs is restated as a finding of this lane, and nothing of
#: theirs is strengthened.
PRIOR_DIRECTION = {
    "file": "data/results/crispri_direction_both.json",
    "population": (
        "held-out signed CRISPRi pairs in reach, K562 and GM12878, joined to the per-element response "
        "cache: 152 pairs, 116 measured down and 36 measured up"
    ),
    "downward_arm": {"k": 96, "n": 116, "rate": 0.8276},
    "upward_arm": {"k": 15, "n": 36, "rate": 0.4167},
    "balanced_accuracy": {"rate": 0.6221, "ci95": [0.5354, None], "band": "0.55 to 0.65"},
    "registered_verdict_carried_verbatim": (
        "balanced accuracy 0.6221 is above chance but under 0.65; direction is faintly readable and not "
        "usable. UNDECIDABLE: the represses half of action stays unsupported, and every action word "
        "GenomeOS writes on a predicted rise is unbacked"
    ),
    "registered_upward_reading_carried_verbatim": (
        "upward sensitivity 0.4167 is at or below chance 0.5. The layer cannot call an increase, and the "
        "passing direction headline is a statement about downward effects only"
    ),
    "not_this_lanes_population": (
        "a different population under a different eligibility rule, so this lane neither confirms nor "
        "contradicts it by construction, and does not quote its figures as its own"
    ),
}

# ---- the intervals, over independent loci ---------------------------------------------------------

DRAWS = 2000
SEED = 20261002
INTERVAL_CALL = (
    f"each arm's rate, BALANCED ACCURACY, the difference between the arms, the baseline and the "
    f"EXCESS of that difference over the baseline all carry a cluster bootstrap interval over "
    f"INDEPENDENT LOCI and not over links: the resampling unit is the independent locus under "
    f"{LOCUS_RULE_IMPORTED_FROM}, loci are drawn with replacement to the observed number of loci, "
    f"every answered link of a drawn locus enters the draw, and the 2.5th and 97.5th percentiles of "
    f"{DRAWS} draws at seed {SEED} are the interval. Loci are grouped over the two arms TOGETHER, so a "
    "locus holding links of both arms is one unit and is drawn or not drawn as one; that is what makes "
    "the difference's interval carry the dependence between the arms rather than hide it. A draw in "
    "which either arm has no answered link yields no difference and is skipped, and the number skipped "
    "is reported. Resampling links instead of loci would treat several links at one locus as several "
    "independent draws, which cell2's convention exists to forbid"
)

# ---- the sign-shuffled control -------------------------------------------------------------------

SHUFFLES = 2000
SHUFFLE_SEED = 20261003
CONTROL_CALL = (
    f"the control is a sign shuffle, {SHUFFLES} of them at seed {SHUFFLE_SEED}. In each shuffle the "
    "measured signs of the answered links are permuted among those links and the primary statistic is "
    "recomputed from the permuted signs. HELD FIXED in every shuffle: the set of answered links; each "
    "link's element, gene, cell, locus and its predicted value and predicted sign; the number of "
    "answered links in each arm, because the permutation is of the existing multiset of measured signs "
    "and so preserves both arm sizes exactly; the locus grouping; and the statistic itself. WHAT "
    "VARIES: only which link carries which measured sign. The control's reported figures are the "
    "median permuted balanced accuracy and permuted difference, their 2.5th and 97.5th percentiles, "
    "and for each the one-sided proportion of shuffles at or above the observed value, which is that "
    "statistic's control p. The reading is taken on the balanced-accuracy p. The permuted difference "
    "is NOT centred on zero: it is centred on 2q - 1, the "
    "baseline a model that calls one sign more often than the other earns for free, so the control is "
    "exactly the comparison that stops that baseline being read as skill. A difference the permutation "
    "reproduces as often as not is a difference this population does not carry"
)
CONTROL_ALPHA = 0.05
CONTROL_LIMIT = (
    "a link-level permutation does not preserve the correlation between links at one locus, so the "
    "control's spread may be narrower than a null that did. The dependence between links is carried by "
    "the locus bootstrap interval and not by the control, which is why BOTH are required for the "
    "established reading and either one alone is not enough"
)

# ---- the gate, taken before any direction is read ------------------------------------------------

GATE_CALL = (
    "the gate is taken first and on each arm separately: an arm's answered links must be at or above "
    f"{POSITIVE_FLOOR} and its answered links' independent loci at or above {LOCUS_FLOOR}, both "
    "imported. The arms are never pooled to reach a floor. A gate reading is not a reading of the "
    "comparison: if either arm is below either floor the comparison is NOT made, no rate of that arm "
    "is read as a result and no floor is moved"
)
GATE_PASS = (
    "both arms are at or above both imported floors on their answered links: the comparison between "
    "the arms may be read"
)
GATE_NO_GO = (
    "an arm is below an imported floor: a no-go. The comparison between the arms is not made, no "
    "agreement rate is read as a result, no floor is moved, the two arms are not pooled to reach a "
    "floor and no answered set is widened by substituting a cell, a gene or another cache. The counts "
    "go on the record with the margin each falls short by, because a reader is owed it and not because "
    "a small margin is better than a large one. A count short of a floor is reported as a no-go and "
    "NEVER as close, promising, nearly enough, a good start or enough for a pilot"
)

# ---- the two readings of the comparison, with no third -------------------------------------------

ESTABLISHED = (
    "sign agreement above the sign-shuffled control is established on this population: the predicted "
    "deletion direction reads the measured sign on the K562 links of the CRISPRi benchmark that lie on "
    "an attributed element, by more than a caller that only prefers one sign. The two arms' own rates "
    "are reported apart beside it and the asymmetry between them is a description, never a second "
    "finding. Both conditions registered here are met - the "
    "cluster bootstrap interval of BALANCED ACCURACY over independent loci lies entirely above 0.5, "
    "and the observed balanced accuracy is above the sign-shuffled control at the one-sided level "
    f"{CONTROL_ALPHA}. This is a statement about SIGN AGREEMENT on this population and it is nothing "
    "else: it is not a statement that the model reads direction in general, not a statement that any "
    "increase is a repression mechanism, not a silencer or repressor call, and not a verdict on any "
    "compiled rule, because no compiled rule was consulted. The locus count it is taken over is "
    "cell2's operational grouping and not established biological independence"
)
NOT_DETECTED = (
    "no sign agreement above the sign-shuffled control is detected on this population. Either the "
    "cluster bootstrap interval of balanced accuracy over independent loci includes 0.5, or the "
    "observed balanced accuracy is not "
    f"above the sign-shuffled control at the registered one-sided level {CONTROL_ALPHA}, or both. The "
    "two arms' rates and denominators go on the record exactly as measured, no floor, predicate, "
    "shuffle count or draw count is moved, no arm is widened, no cell is substituted and the test is "
    "not re-taken under a second statistic to find a difference. A result that does not detect a "
    "difference is reported as not detected and NEVER as close, promising, nearly enough, a good start "
    "or enough for a pilot: it has no encouraging reading here, and the margin by which it misses is "
    "reported because a reader is owed it and not because a small margin is better than a large one. "
    "It is also not a statement that there is no difference: an undetected difference on this "
    "population is an undetected difference, not a measured absence"
)
READINGS = (ESTABLISHED, NOT_DETECTED)
THERE_IS_NO_THIRD = True

#: The words no summary, result or commit message of this lane may use of a short or undetected
#: outcome. lane-increase forbids them by name and this lane carries the prohibition forward.
FORBIDDEN_OF_A_SHORT_OR_UNDETECTED_OUTCOME = (
    "close",
    "promising",
    "nearly enough",
    "a good start",
    "enough for a pilot",
)

WHAT_FOLLOWS = {
    "established": (
        "stop at the reading and report. The difference is reported with both arms' denominators, both "
        "intervals, the control and every unanswerable count beside it; it is not carried into a claim "
        "about any compiled rule, because the rule-level test is impossible on this population and "
        "this lane did not attempt one"
    ),
    "not_detected": (
        "stop at the reading and report. Nothing is re-counted under a different predicate, a different "
        "statistic or a wider answered set to find a difference, and the floors stand"
    ),
}

# ---- what the result cannot establish ------------------------------------------------------------

CANNOT_ESTABLISH = (
    "a sign agreement is not a mechanism. Whatever this test says, nothing here establishes that any "
    "of these increases is a repression mechanism, that an element whose silencing raises a gene "
    "represses it, that any element is a silencer or a repressor, or that any compiled repression call "
    "is right or wrong: no compiled rule is consulted anywhere in this lane, by construction, because "
    "only one of the 48 increases has a rule naming its gene and none has one gated on its cell. "
    "Nothing here is a verdict on the deletion model's direction in general: the test is on one cell "
    "(K562), one assay (the ENCODE CRISPRi benchmark), and only those links that lie on an attributed "
    "element under the committed overlap rule, so it is not a genome-wide, cell-type-general or "
    "assay-general statement. A difference between the arms is a difference in sign agreement and not "
    "a measure of how large either effect is, since no magnitude bar is applied and the floors are on "
    "counts. The locus convention is cell2's operational grouping, never established biological "
    "independence, so an interval taken over loci is an interval over that grouping and not over "
    "established independent observations. The arms are not balanced and are not matched on magnitude, "
    "chromosome, gene, split or dataset, so a difference between them may be carried by any of those "
    "rather than by the sign. An undetected difference is not a measured absence of one. And the 9 "
    "WTC11 and HCT116 links are absent from the test rather than negative in it: this lane says "
    "nothing whatever about direction in a cell the cached answers do not carry"
)

NO_REQUESTS = (
    "0 model requests, no money, no download and no network. Every input is already on disk: the two "
    "benchmark tables, the per-chromosome attribution results and the per-element response cache of "
    "the finished sweep, which was written with threshold=0.0 and so already carries a signed per-cell "
    "value for every gene in the scorer's window"
)

ORDER = (
    "1. the eligibility join, blind to every predicted value: the links of both arms at the shared "
    "`and_on_an_attributed_element` step, with the contested links counted on their own",
    "2. the answerability pass: for each link the cached value for its own gene on its own cell track, "
    "with absent and exactly-zero counted and excluded, and every cell the cache does not carry "
    "reported unanswerable without substitution",
    "3. the gate on each arm separately against both imported floors",
    "4. if and only if the gate passes: the two arms' rates reported apart, their cluster bootstrap "
    "intervals over independent loci, the difference and its interval, and the sign-shuffled control",
    "5. the reading, in the words registered here before the test was run, and the section saying what "
    "the result cannot establish",
)


# ---- the predicates and estimators the computing script applies -----------------------------------


def significant(pair: Any) -> bool:
    """Membership in the two significant outcome labels: repress2's own predicate, imported."""
    return rp.significant({"outcome": getattr(pair, "outcome", None) or pair.get("outcome")})


def sign_of(value: float | None) -> int | None:
    """SIGN_CALL in code: -1, +1, or None for exactly zero. `None` in means absent, `None` out."""
    if value is None:
        return None
    v = float(value)
    if v > 0:
        return 1
    if v < 0:
        return -1
    return None


def arm_of(has_decrease: bool, has_increase: bool) -> str | None:
    """The arm of one link, or None when the measurement contests itself (CONTEST_RULE)."""
    if has_decrease and has_increase:
        return None
    if has_decrease:
        return DECREASES
    if has_increase:
        return INCREASES
    return None


def rate(links: list[dict[str, Any]]) -> float | None:
    """One arm's agreement rate over its ANSWERED links, or None when it has none."""
    if not links:
        return None
    return sum(1 for r in links if r["agrees"]) / len(links)


def difference(by_arm: dict[str, list[dict[str, Any]]]) -> float | None:
    """PRIMARY_STATISTIC: decreases minus increases, or None when either arm is empty."""
    a, b = rate(by_arm.get(DECREASES) or []), rate(by_arm.get(INCREASES) or [])
    if a is None or b is None:
        return None
    return a - b


def balanced_accuracy(links: list[dict[str, Any]]) -> float | None:
    """TEST_STATISTIC: the unweighted mean of the two arms' agreement rates. None if either arm is empty.

    Its chance level is CHANCE for any class mix, and a constant-sign caller of either sign scores
    exactly CHANCE on it by construction - which is the whole reason it is the primary here.
    """
    by_arm: dict[str, list[dict[str, Any]]] = {DECREASES: [], INCREASES: []}
    for r in links:
        if r["arm"] in by_arm:
            by_arm[r["arm"]].append(r)
    a, b = rate(by_arm[DECREASES]), rate(by_arm[INCREASES])
    if a is None or b is None:
        return None
    return (a + b) / 2


def marginal_down(links: list[dict[str, Any]]) -> float | None:
    """q of MARGINAL_CALL: the share of answered links whose PREDICTED sign is -1."""
    if not links:
        return None
    return sum(1 for r in links if r["predicted_sign"] == -1) / len(links)


def baseline(links: list[dict[str, Any]]) -> float | None:
    """BASELINE: 2q - 1, what the sign shuffle earns for free on this answered set."""
    q = marginal_down(links)
    return None if q is None else 2 * q - 1


def excess(links: list[dict[str, Any]]) -> float | None:
    """TEST_STATISTIC: the difference between the arms minus the baseline the shuffle implies."""
    by_arm = {DECREASES: [], INCREASES: []}
    for r in links:
        if r["arm"] in by_arm:
            by_arm[r["arm"]].append(r)
    d, b = difference(by_arm), baseline(links)
    if d is None or b is None:
        return None
    return d - b


def locus_keys(links: list[dict[str, Any]]) -> list[cell2.LocusKey]:
    """One cell2 locus key per link: the only path from a link into the grouping."""
    return [cell2.LocusKey(r["cell"], r["chrom"], int(r["start"]), int(r["end"]), r["gene"]) for r in links]


def loci_of(links: list[dict[str, Any]]) -> list[int]:
    """The locus group id of each link, grouped over the links GIVEN, by cell2's rule."""
    return cell2.group(locus_keys(links), LOCUS_SPAN)


def gate(links: list[dict[str, Any]], arm: str) -> dict[str, Any]:
    """One arm against both imported floors, on its ANSWERED links. GATE_CALL in code."""
    n = len(links)
    loci = len(set(loci_of(links))) if links else 0
    met = n >= POSITIVE_FLOOR and loci >= LOCUS_FLOOR
    short: list[str] = []
    if n < POSITIVE_FLOOR:
        short.append(f"answered links {n}, short of {POSITIVE_FLOOR} by {POSITIVE_FLOOR - n}")
    if loci < LOCUS_FLOOR:
        short.append(f"independent loci {loci}, short of {LOCUS_FLOOR} by {LOCUS_FLOOR - loci}")
    return {
        "arm": arm,
        "answered_links": n,
        "independent_loci": loci,
        "floors": {"links": POSITIVE_FLOOR, "independent_loci": LOCUS_FLOOR},
        "meets_both_floors": met,
        "short_by": short,
        "independent_locus_rule": LOCUS_RULE,
        "not_biological_independence": NOT_BIOLOGICAL_INDEPENDENCE,
    }


def interval(xs: list[float]) -> list[float] | None:
    """The 2.5th and 97.5th percentiles of a draw list, the project's own convention."""
    if not xs:
        return None
    s = sorted(xs)
    return [round(s[int(0.025 * (len(s) - 1))], 4), round(s[int(0.975 * (len(s) - 1))], 4)]


def bootstrap(answered: list[dict[str, Any]], n: int = DRAWS, seed: int = SEED) -> dict[str, Any]:
    """INTERVAL_CALL in code: resample INDEPENDENT LOCI, not links, over both arms together."""
    groups: dict[int, list[dict[str, Any]]] = {}
    for g, r in zip(loci_of(answered), answered, strict=True):
        groups.setdefault(g, []).append(r)
    ids = sorted(groups)
    rng = random.Random(seed)
    draws: dict[str, list[float]] = {
        DECREASES: [],
        INCREASES: [],
        "difference": [],
        "baseline": [],
        "excess": [],
        "balanced_accuracy": [],
    }
    skipped = 0
    for _ in range(n):
        sample = [r for i in (rng.choice(ids) for _ in ids) for r in groups[i]]
        by_arm: dict[str, list[dict[str, Any]]] = {DECREASES: [], INCREASES: []}
        for r in sample:
            by_arm[r["arm"]].append(r)
        d = difference(by_arm)
        if d is None:
            skipped += 1
            continue
        draws[DECREASES].append(rate(by_arm[DECREASES]))  # type: ignore[arg-type]
        draws[INCREASES].append(rate(by_arm[INCREASES]))  # type: ignore[arg-type]
        draws["difference"].append(d)
        draws["baseline"].append(baseline(sample))  # type: ignore[arg-type]
        draws["excess"].append(excess(sample))  # type: ignore[arg-type]
        draws["balanced_accuracy"].append(balanced_accuracy(sample))  # type: ignore[arg-type]
    return {
        "unit": "independent locus",
        "loci": len(ids),
        "draws": n,
        "seed": seed,
        "skipped_draws": skipped,
        "skipped_call": "a draw in which either arm held no answered link yields no difference",
        "ci95": {k: interval(v) for k, v in draws.items()},
        "statistic_the_reading_is_taken_on": TEST_STATISTIC,
        "over_loci_not_links": INTERVAL_CALL,
    }


def shuffle_control(
    answered: list[dict[str, Any]], observed: float, n: int = SHUFFLES, seed: int = SHUFFLE_SEED
) -> dict[str, Any]:
    """CONTROL_CALL in code: permute the measured signs among the answered links, nothing else."""
    signs = [r["measured_sign"] for r in answered]
    preds = [r["predicted_sign"] for r in answered]
    rng = random.Random(seed)
    diffs: list[float] = []
    bals: list[float] = []
    at_or_above = 0
    bal_at_or_above = 0
    observed_balanced = balanced_accuracy(answered)
    for _ in range(n):
        perm = signs[:]
        rng.shuffle(perm)
        agree = {-1: [0, 0], 1: [0, 0]}  # measured sign -> [agreements, links]
        for sg, pr in zip(perm, preds, strict=True):
            agree[sg][1] += 1
            agree[sg][0] += 1 if pr == sg else 0
        if not agree[-1][1] or not agree[1][1]:
            continue
        down = agree[-1][0] / agree[-1][1]
        up = agree[1][0] / agree[1][1]
        diffs.append(down - up)
        bals.append((down + up) / 2)
        at_or_above += 1 if down - up >= observed else 0
        if observed_balanced is not None and (down + up) / 2 >= observed_balanced:
            bal_at_or_above += 1
    s = sorted(diffs)
    mid = len(s) // 2
    return {
        "shuffles": n,
        "seed": seed,
        "held_fixed": CONTROL_CALL,
        "observed_difference": round(observed, 4),
        "centred_on_the_baseline_not_on_zero": MARGINAL_CALL,
        "baseline": None if baseline(answered) is None else round(baseline(answered), 4),
        "marginal_predicted_down_rate": (
            None if marginal_down(answered) is None else round(marginal_down(answered), 4)
        ),
        "median_permuted_difference": (
            None if not s else round(s[mid] if len(s) % 2 else (s[mid - 1] + s[mid]) / 2, 4)
        ),
        "ci95_permuted_difference": interval(diffs),
        "shuffles_at_or_above_observed_difference": at_or_above,
        "p_one_sided_difference": None if not diffs else round(at_or_above / len(diffs), 4),
        "observed_balanced_accuracy": (None if observed_balanced is None else round(observed_balanced, 4)),
        "median_permuted_balanced_accuracy": (
            None
            if not bals
            else round(
                sorted(bals)[len(bals) // 2]
                if len(bals) % 2
                else (sorted(bals)[len(bals) // 2 - 1] + sorted(bals)[len(bals) // 2]) / 2,
                4,
            )
        ),
        "ci95_permuted_balanced_accuracy": interval(bals),
        "shuffles_at_or_above_observed_balanced_accuracy": bal_at_or_above,
        "p_one_sided": None if not bals else round(bal_at_or_above / len(bals), 4),
        "p_one_sided_call": (
            "`p_one_sided` is the balanced-accuracy p, which is the one the reading is taken on. The "
            "difference's own p is reported beside it as a description and COLLAPSE says why it is not "
            "the reading"
        ),
        "alpha": CONTROL_ALPHA,
        "limit": CONTROL_LIMIT,
    }


def reading(ci_balanced: list[float] | None, p_one_sided: float | None) -> dict[str, Any]:
    """The two registered readings, with no third.

    `ci_balanced` is the locus bootstrap interval of TEST_STATISTIC, balanced accuracy, whose chance
    level is CHANCE for any class mix and which a constant-sign caller of either sign scores exactly.
    The difference between the arms is reported with its own interval, its baseline and its excess, but
    the reading is not taken on it: COLLAPSE says why.
    """
    above_zero = ci_balanced is not None and ci_balanced[0] > CHANCE
    beats_control = p_one_sided is not None and p_one_sided < CONTROL_ALPHA
    established = above_zero and beats_control
    return {
        "established": established,
        "reading": ESTABLISHED if established else NOT_DETECTED,
        "balanced_accuracy_interval_above_chance": above_zero,
        "above_the_shuffled_control": beats_control,
        "both_required": (
            "the established reading requires BOTH the locus interval of balanced accuracy entirely "
            f"above {CHANCE} AND its one-sided control p below {CONTROL_ALPHA}; either alone is not "
            "enough"
        ),
        "there_is_no_third": THERE_IS_NO_THIRD,
        "what_follows": WHAT_FOLLOWS["established" if established else "not_detected"],
        "forbidden_wording": list(FORBIDDEN_OF_A_SHORT_OR_UNDETECTED_OUTCOME),
    }


def registration() -> dict[str, Any]:
    """Everything this lane fixes before any direction is read, as the result payload."""
    return {
        "question": (
            "among measured significant effects on attributed elements, does the model's predicted "
            "deletion direction for the link's own gene and cell match the measured sign more often "
            "for decreases than for increases, and above a sign-shuffled control"
        ),
        "level": "link level, not rule level",
        "why_not_rule_level": (
            "of lane-increase's 48 increases exactly 1 has any compiled rule naming its own gene and 0 "
            "have one gated on its own cell (data/results/increase_population.json), so a rule-level "
            "test of direction on the increases is impossible. This lane does not attempt one: the "
            "model's direction calls are read directly out of the cached deletion answers and never "
            "through a compiled rule"
        ),
        "follows_from": {
            "lane": "lane-increase",
            "commits": INHERITED_INCREASE_COMMITS,
            "registered_go_carried_verbatim": INHERITED_INCREASE_GO,
            "figures": INHERITED_FIGURES,
            "and_behind_it": {
                "lane": "lane-repress2",
                "file": INHERITED_REPRESS2_FILE,
                "registered_no_go_carried_verbatim": INHERITED_REPRESS2_NO_GO,
            },
        },
        "binding_wording_r2": (
            "a measured increase on knockdown is 'an increase on knockdown'. It is never a silencer, "
            "never a repressor and never evidence of a repression mechanism. This test is about sign "
            "agreement, not mechanism"
        ),
        "eligibility_predicate": ELIGIBILITY,
        "contest_rule": CONTEST_RULE,
        "arms": list(ARMS),
        "arms_never_pooled": (
            "the two arms are reported apart and are never pooled. The comparison between them is the "
            "result; neither arm alone is. A pooled agreement rate is not reported as a result"
        ),
        "decrease_arm_is_not_model_coverage": DECREASE_ARM_IS_NOT_MODEL_COVERAGE,
        "measured_sign": MEASURED_SIGN_CALL,
        "predicted_direction": PREDICTED_CALL,
        "sign_rule_zero_and_absent": SIGN_CALL,
        "agreement": AGREEMENT_CALL,
        "primary_statistic": PRIMARY_STATISTIC,
        "marginal_sign_rate_and_the_baseline": MARGINAL_CALL,
        "baseline": BASELINE,
        "excess": EXCESS_CALL,
        "why_the_difference_is_not_the_reading": COLLAPSE,
        "test_statistic": TEST_STATISTIC,
        "balanced_accuracy": BALANCED_CALL,
        "chance": CHANCE,
        "prior_direction_lanes": PRIOR_DIRECTION,
        "cells": {
            "cached_cell_tracks": list(CACHED_CELLS),
            "cells_not_cached": list(CELLS_NOT_CACHED),
            "availability": CELL_AVAILABILITY,
            "primary_cell": PRIMARY_CELL,
            "primary_cell_call": PRIMARY_CELL_CALL,
            "not_balanced": NOT_BALANCED,
        },
        "split_rule": SPLIT_RULE,
        "locus_convention": {
            "rule": LOCUS_RULE,
            "span": LOCUS_SPAN,
            "imported_from": LOCUS_RULE_IMPORTED_FROM,
            "not_biological_independence": NOT_BIOLOGICAL_INDEPENDENCE,
            "inherited_33_loci": INHERITED_33_LOCI_ARE_OPERATIONAL,
        },
        "floors": FLOORS,
        "gate": {"call": GATE_CALL, "pass": GATE_PASS, "no_go": GATE_NO_GO},
        "intervals": {
            "call": INTERVAL_CALL,
            "draws": DRAWS,
            "seed": SEED,
            "unit": "independent locus",
        },
        "control": {
            "call": CONTROL_CALL,
            "shuffles": SHUFFLES,
            "seed": SHUFFLE_SEED,
            "alpha": CONTROL_ALPHA,
            "limit": CONTROL_LIMIT,
        },
        "readings_before_the_outcome_was_seen": {
            "established": ESTABLISHED,
            "not_detected": NOT_DETECTED,
            "there_is_no_third": THERE_IS_NO_THIRD,
            "forbidden_wording": list(FORBIDDEN_OF_A_SHORT_OR_UNDETECTED_OUTCOME),
            "what_follows": WHAT_FOLLOWS,
        },
        "cannot_establish": CANNOT_ESTABLISH,
        "order": list(ORDER),
        "requests": 0,
        "money": "none: every input is already on disk",
        "no_requests": NO_REQUESTS,
        "nothing_changed": (
            "no extractor, compiled rule, threshold, floor, predicate, shuffle count, draw count or "
            "registered reading of any lane is changed by this lane, and no registered reading is "
            "restated in stronger words than the lane that registered it used"
        ),
        "what_the_run_will_write": "data/results/direction_link.json",
    }
