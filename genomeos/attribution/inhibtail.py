# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-inhibtail (ROADMAP area I, 2026-10-12): is the compiled layer's ABOVE-BOUND INHIBITOR TAIL a
property of the PREDICTION INSTRUMENT rather than of the biology?

STOP-LIST CHECK (docs/ROADMAP.md §5, "The stop list of 2026-09-27, reinstated"; the six items in
docs/ROADMAP-HISTORY.md, "The stop list, which is as much the plan as the waves"). The six, each named
and each untouched by this lane: (1) nothing more is built on the node model -- no node model is read or
written here; (2) no more biosamples -- no biosample is added, the compiled cell labels already in the
programs are counted and nothing is fetched; (3) nothing is built in area J -- this is area I; (4) no
AlphaGenome requests on the 882 blocks -- 0 requests, no key read, no network; (5) area E's fate ceiling
is not chased -- no fate, no lineage, no contact atlas; (6) BioForge's confidence is not calibrated --
BioForge is not touched. NO RNA VALUE IS READ AND NOTHING ABOUT RNA IS CLAIMED. NOTHING IS FITTED AND NO
FIT CHANGE IS PROPOSED: every number here is arithmetic on compiled rules' own stated removal responses
plus M2h's fixed scale through `scalefit.SCALE`.

THE QUESTION
------------
`lane-inhibitor` (`a4c79ce`, `2816a3e`, `85ca300`, `a1a9a90`, `7d44158`) measured that the 39 genes
fitted by neither form are all refused on `inhibitor_strength_above_1`, that 38 of the 39 carry exactly
ONE offending element out of a median 15 inhibitors, and that ALL 40 OFFENDERS ARE PREDICTED DELETION
NOTES of log2 fold 2.6 to 6.6. It also measured that the refusal is general to the compiled layer: above
the bound, 1,160 genes outside the 39 (1,199 in all, 1,330 elements) unscaled, and 149 outside the 39
(188 in all, 192 elements) at the scale.

So: IS THAT TAIL SOMETHING THE PREDICTOR DOES? If the 1,330 above-bound inhibiting elements of the
17,983 eligible multi-element genes are composed of predicted deletion notes IN A WAY THE REST OF THE
155,055 INHIBITING ELEMENTS OF THE SAME GENES IS NOT, the tail is the instrument's. If they are not, it
is not.

THE BASE-RATE TRAP, WHICH DECIDES WHETHER THIS LANE IS WORTH ANYTHING
--------------------------------------------------------------------
"All 40 offenders are predicted" is UNINFORMATIVE ON ITS OWN if essentially every inhibiting element in
the layer is predicted. `lane-emit7g`'s census (`589b6ca`, `data/results/emit7g_census.json`
`genome_wide.now`) counts 440,589 compiled rule lines = 440,377 predicted + 212 experimental, so
predicted rules are the overwhelming majority BEFORE ANYTHING IS MEASURED. A lane that reports "the tail
is all predicted" without the base rate has measured nothing.

THEREFORE THE BASE RATE IS COMPUTED AND COMMITTED FIRST, as `data/results/inhibtail_registration.json`,
by `scripts/inhibtail_run.py --stage baserate`, BEFORE ANY TAIL NUMBER IS READ. That stage is incapable
of reading the tail: it conditions on nothing but RHO > 1, and `BASE_RATE_FIELDS` is the whole of what it
emits. The tail stage runs only after the registration is committed.

WHAT IS ALREADY DETERMINED BY COMMITTED DATA, before this lane computes anything
--------------------------------------------------------------------------------
Three of the axes a reader would reach for are already settled by committed records or by reading the
code, and saying so in advance is the point of this section -- not computing them and then presenting
the answer as a finding.

* COMPOSITION BY EVIDENCE KIND. `inhibitor_census_result.json`
  `Q5_where_the_rules_came_from.background_all_inhibiting_elements_of_all_eligible_genes.by_evidence_kind`
  is already `{"predicted": 155055}` and `by_how_rho_was_read` already
  `{"evidence note, log2 fold change": 155055}`. THAT FIELD WAS READ BY THIS LANE BEFORE THIS
  REGISTRATION WAS WRITTEN, and the registration says so rather than pretending otherwise (see
  `DISCLOSURE`). It is a committed, inherited field, not a number this lane produced; this lane
  re-measures it in the base-rate stage anyway, because a registration whose base rate is inherited
  cannot check its own arithmetic.
* THE TAUTOLOGY. `bridge.observation_of` returns, for a predicted rule whose evidence note carries an
  effect, `Observation(2.0 ** log2_fold, ...)`; `IMPLIED_STRENGTH(RHO) = (1 - 1/RHO) / H` with
  H = 0.5 is strictly increasing in RHO on RHO > 1. So "strength above 1" IS "RHO above 2" IS "the
  note's log2 fold above 1.0", exactly, and at M2h's scale "RHO above 5.660" IS "the note's log2 fold
  above 2.501". The tail is therefore the upper tail of the predictor's own stated effect sizes BY
  CONSTRUCTION, and the axis "are the tail's elements the ones where the deletion note's effect size is
  largest" is TAUTOLOGICAL. It is registered as tautological, verified as an identity, and it CHOOSES NO
  OUTCOME.
* SHARING. `589b6ca`'s message records `compiled_targets_1` as 440,377 of 440,377 with
  `compiled_targets_2_or_more` 0 -- one element, one target, one rule. So no element is shared between
  genes, in the tail or out of it, and the sharing axis is a property of the compiler's emission rather
  than of the tail. Registered with an EXPECTED VALUE of 0 shared elements, verified, and it CHOOSES NO
  OUTCOME.

WHY THE COMPOSITION AXIS IS DEAD BY DEFINITION AND NOT BY MEASUREMENT -- read off the code, before
the run, and the single most important thing in this file
---------------------------------------------------------------------------------------------------
`combine2.multi_rule_genes`, which builds POPULATION (a) that every downstream lane's "eligible
multi-element genes" means, drops a gene on three tests, and the third is:

    if any(k != "predicted" for _, k in rows):
        why["a_rule_not_predicted"] += 1
        continue

SO A GENE IS ADMITTED TO THIS POPULATION ONLY IF EVERY ONE OF ITS RULES IS PREDICTED. The base rate of
`predicted` over the inhibiting elements of the eligible multi-element genes is therefore EXACTLY 1 BY
THE POPULATION'S DEFINITION, not as an empirical finding, and `inhibitor_census_result.json`'s
`{"predicted": 155055}` is that definition being re-read rather than a measurement. The same holds for
`bridge.observation_of`'s `how`: `classify` drops a gene with any unreadable observation, and the only
`how` a predicted rule with a readable note returns is "evidence note, log2 fold change".

THAT IS THE ANSWER TO THE COMPOSITION AXIS AND IT IS WRITTEN HERE BEFORE THE RUN. The axis is not dead
because the numbers came out one-sided; it is dead because the population was constructed by excluding
every element that could have populated the other arm. The run MEASURES it anyway -- a definition read
from code and never checked against the data is how a lane convinces itself of a thing that is not
there -- but the measurement's expected value is declared here first: fraction 1.0, every field, both
populations, and 0 non-predicted inhibiting elements.

THE CONTROL ARM, named before it is counted
-------------------------------------------
The instrument-versus-biology contrast needs an experimental arm inside the same population. The layer
holds 212 experimental rule lines against 440,377 predicted. `AXES["A9_control_arm"]` counts how many
experimental rules fall on the eligible multi-element genes, how many of those yield a readable
observation, and how many of THOSE are inhibiting (RHO > 1). ITS EXPECTED VALUE IS 0 AND THE REASON IS
THE DEFINITION ABOVE, not a scarcity of experiments: a gene carrying ANY experimental rule is excluded
from the population before an inhibitor of it is ever counted. So there is no experimental inhibitor
anywhere in this population and there CANNOT BE ONE, the contrast has no control arm inside the layer at
all, and that is declared here, in advance, as a COMPLETE AND HONEST ANSWER on the composition axis and
not a failure.

This file is the RECORD. Every question, population, decision rule, outcome wording and expected value
is committed here before `scripts/inhibtail_run.py` computes any of them.
"""

from __future__ import annotations

import bisect
import math
from typing import Any

from genomeos.attribution import bridge, scalefit

LANE = "lane-inhibtail"
AREA = "I"
DATE = "2026-10-12"

#: M2h's fixed scale, taken verbatim; this lane does not refit it and does not propose changing it.
SCALE = scalefit.SCALE
#: The bridge's intact Hill occupancy, the constant the bound is algebra on.
H = bridge.H_INTACT
#: The fit's slack on the bound, taken from the bridge.
STRENGTH_TOLERANCE = bridge.STRENGTH_TOLERANCE

#: RHO above which IMPLIED_STRENGTH exceeds 1 on an UNSCALED observation. Algebra, not measurement:
#: s = (1 - 1/RHO)/H > 1  <=>  RHO > 1/(1 - H) = 2 at H = 0.5.
BOUND_RHO_UNSCALED = 1.0 / (1.0 - H)
#: The same bound after every drop is scaled by SCALE: 1 + SCALE*(RHO-1) > 2  <=>  RHO > 1 + 1/SCALE.
BOUND_RHO_AT_THE_SCALE = 1.0 + 1.0 / SCALE

QUESTION = (
    "Is the compiled layer's above-bound inhibitor tail a property of the PREDICTION INSTRUMENT rather"
    " than of the biology? All 40 offenders of lane-inhibitor's 39 genes were predicted deletion notes."
    " The tail is 1,330 inhibiting elements of the 17,983 eligible multi-element genes. If the tail is"
    " composed of predicted deletion notes in a way the rest of that population's 155,055 inhibiting"
    " elements is not, the tail is something the predictor does; if it is not, it is not. The base rate"
    " is computed and committed FIRST, because 440,377 of the layer's 440,589 rule lines are predicted"
    " before anything is measured and a composition claim read against no base rate measures nothing."
)

# ---------------------------------------------------------------------------------------------------
# POPULATIONS. Every number this lane reports names one of these in the same clause.
# ---------------------------------------------------------------------------------------------------

POPULATIONS = {
    "eligible_genes": (
        "THE ELIGIBLE MULTI-ELEMENT GENES: every gene `combine2.classify` admits across the 24 compiled"
        " non-coding programs, 17,983 of them in lane-inhibitor's committed run. Membership is"
        " re-derived here by the same code path and its count is reported, not assumed."
    ),
    "inhibiting_elements": (
        "THE INHIBITING ELEMENTS OF THE ELIGIBLE MULTI-ELEMENT GENES: every element of such a gene whose"
        " stated removal response RHO exceeds 1, 155,055 of them in lane-inhibitor's committed run. This"
        " is the population the base rate is taken over and the population the tail is a subset of."
    ),
    "non_inhibiting_elements": (
        "THE NON-INHIBITING ELEMENTS OF THE SAME GENES: every element of an eligible multi-element gene"
        " with a readable RHO that does not exceed 1 (an activator, or exactly 1). The base rate is"
        " reported over this population too, so that 'predicted' can be seen to be the layer's answer"
        " rather than the inhibitors'."
    ),
    "the_tail_unscaled": (
        "THE TAIL (UNSCALED): the inhibiting elements of the eligible multi-element genes whose"
        " IMPLIED_STRENGTH on the UNSCALED observation exceeds 1 + STRENGTH_TOLERANCE, equivalently"
        " RHO > 2. 1,330 elements in lane-inhibitor's committed run. THIS IS THE PRIMARY TAIL."
    ),
    "the_body_unscaled": (
        "THE BODY (UNSCALED): the inhibiting elements of the same genes that are not in the tail;"
        " 155,055 - 1,330 = 153,725 in lane-inhibitor's committed run. Every enrichment compares the"
        " tail against THE BODY, never against the whole population, so the two arms are disjoint."
    ),
    "the_tail_at_the_scale": (
        "THE TAIL AT THE SCALE: the same elements judged on the SCALED observation, RHO > 5.660;"
        " 192 elements in lane-inhibitor's committed run. Reported beside the primary tail and used for"
        " no decision rule of its own."
    ),
    "the_39": (
        "THE 39: lane-residuemult's genes fitted by NEITHER form, all on `inhibitor_strength_above_1`"
        " (lane-inhibitor re-derived them and the sets agreed). This lane does not re-derive them; it"
        " PINS their committed figures and reads the tail they sit in."
    ),
}

# ---------------------------------------------------------------------------------------------------
# DISCLOSURE. What makes a registration honest is git, not intention.
# ---------------------------------------------------------------------------------------------------

DISCLOSURE = (
    "WHAT THIS LANE HAD ALREADY READ WHEN THIS REGISTRATION WAS WRITTEN, disclosed because the"
    " registration's credibility rests on git rather than on intention. (1) The whole of"
    " data/results/inhibitor_census_result.json's Q5, Q6, Q7, `implication` and `branch` fields,"
    " including"
    " Q5.background_all_inhibiting_elements_of_all_eligible_genes.by_evidence_kind ="
    " {'predicted': 155055} -- so THE BASE RATE ON THE EVIDENCE-KIND FIELD WAS KNOWN, FROM A COMMITTED"
    " INHERITED FIELD, BEFORE THIS FILE WAS WRITTEN, and the composition axis was already expected to be"
    " dead when its decision rule was fixed. The rule and the threshold DEAD_FRACTION = 0.999 are still"
    " fixed here before this lane's own measurement of it. (2)"
    " data/results/emit7g_census.json genome_wide.now, the 440,589 = 440,377 + 212 split. (3) The"
    " bodies of genomeos/attribution/{bridge,combine2,saturate,scalefit,inhibitor_census}.py and"
    " scripts/inhibitor_census_run.py, from which the tautology and the bound's algebra were read."
    " (4) One timing probe on chr21 alone, in the scratchpad and in no commit, which printed the gate's"
    " kind counts (5,174 predicted, 2 experimental), the element count 5,374, 178 multi-rule genes, 176"
    " eligible and 1,889 inhibiting elements on that one chromosome. (5) AND THIS ONE IS A TAIL NUMBER,"
    " SO IT IS NAMED IN FULL RATHER THAN GLOSSED: a smoke test of `scripts/inhibtail_run.py` on CHR21"
    " ALONE, run in the scratchpad before this file was committed and recorded in no commit, which"
    " printed, over that one program's 176 eligible genes and 1,889 inhibiting elements: a base rate of"
    " fraction 1.0 on all three BASE_RATE_FIELDS (predicted; 'evidence note, log2 fold change';"
    " 'AlphaGenome deletion'), 0 inhibiting elements not predicted, 3,156 non-inhibiting elements also"
    " 1.0 predicted, 0 unreadable -- and THEN 11 ABOVE-BOUND ELEMENTS, a maximum RHO of 8.456, a largest"
    " multiplicative gap of 2.567 at RHO 3.294, and survival slopes of -2.393 below and -1.249 above."
    " THOSE CHR21 FIGURES WERE SEEN BEFORE THIS REGISTRATION WAS COMMITTED. They are one program of 24,"
    " they are not the population any decision rule names, and every window, threshold, tolerance and"
    " outcome wording in this file was written before the smoke test ran -- but the fact that a tail"
    " count was seen at all is what a reader needs in order to judge this registration, and it is"
    " recorded here instead of being left to inference. NOTHING ELSE CONDITIONED ON THE BOUND, ON A"
    " STRENGTH OR ON A FOLD ABOVE 2 HAS BEEN COMPUTED BY THIS LANE AT THE TIME THIS FILE IS COMMITTED,"
    " and no number over the 24 programs together has been computed at all."
)

# ---------------------------------------------------------------------------------------------------
# THE BASE RATE, registered and committed before the tail's composition is read.
# ---------------------------------------------------------------------------------------------------

BASE_RATE_FIELDS = (
    "evidence_kind",  # bridge._kind(rule): predicted / experimental / curated / inferred / none
    "how_rho_was_read",  # the `how` string bridge.observation_of returns
    "evidence_source_prefix",  # the evidence source string up to its first comma
)

BASE_RATE_DECLARATION = (
    "THE BASE RATE, computed and committed BEFORE the tail's composition is read, by"
    " `scripts/inhibtail_run.py --stage baserate` into data/results/inhibtail_registration.json. Over"
    " THE INHIBITING ELEMENTS OF THE ELIGIBLE MULTI-ELEMENT GENES, and separately over THE"
    " NON-INHIBITING ELEMENTS OF THE SAME GENES, the counts of every value present of each of the three"
    " BASE_RATE_FIELDS: the evidence KIND the rule states (`bridge._kind`), the HOW string"
    " `bridge.observation_of` returns for it, and the PREFIX of its evidence source up to the first"
    " comma (which is what separates 'AlphaGenome deletion' from a differently-produced note while"
    " discarding the biosample that follows it). Beside them, the counts of compiled rule lines by kind"
    " over the same programs, and the control-arm counts of AXES['A9_control_arm']. THE BASE-RATE STAGE"
    " CONDITIONS ON NOTHING BUT RHO > 1: it emits no strength, no fold threshold, no above-bound count"
    " and no per-element row, so it CANNOT read the tail even by accident."
)

#: The fraction at or above which a composition field is declared to carry no information.
DEAD_FRACTION = 0.999

ENRICHMENT_STATEMENT = (
    "THE ENRICHMENT STATEMENT, written out here before the tail's numbers are seen. For each of the"
    " three BASE_RATE_FIELDS, and for the field's most common value V in the population: 'Of the"
    " {tail_n} above-bound inhibiting elements of the eligible multi-element genes, {tail_v} are {V}"
    " ({tail_frac}); of the {body_n} inhibiting elements of the same genes below the bound, {body_v} are"
    " ({body_frac}). The risk ratio tail-against-body is {rr}, 95% interval [{lo}, {hi}] by the Katz log"
    " method on the two binomial proportions.' THE ENRICHMENT IS DECLARED PRESENT only when that"
    " interval EXCLUDES 1.0, and absent otherwise. WHERE EITHER ARM IS A PROPORTION OF EXACTLY 1 OR"
    " EXACTLY 0 THE KATZ INTERVAL IS UNDEFINED, THE RATIO IS NOT QUOTED, AND THE FIELD IS REPORTED AS"
    " DEGENERATE -- which is the expected case and is not a failure."
)

DEAD_BY_DEFINITION = (
    "THE COMPOSITION AXIS IS DEAD BY THE POPULATION'S DEFINITION AND NOT BY ITS NUMBERS, declared here"
    " before the run. `combine2.multi_rule_genes`, which builds the POPULATION (a) that every"
    " downstream lane's 'eligible multi-element genes' means, drops a gene when"
    ' `any(k != "predicted" for _, k in rows)`, counting it as `a_rule_not_predicted`. A gene is'
    " admitted only if EVERY one of its rules is predicted, so the fraction of `predicted` over the"
    " inhibiting elements of the eligible multi-element genes is EXACTLY 1 BY CONSTRUCTION, and"
    " inhibitor_census_result.json's {'predicted': 155055} is that construction being re-read rather"
    " than a measurement. The same holds for the `how` field: `combine2.classify` drops a gene with any"
    " unreadable observation, and the one `how` a predicted rule with a readable note returns is"
    " 'evidence note, log2 fold change'. THE AXIS IS NOT DEAD BECAUSE THE NUMBERS CAME OUT ONE-SIDED;"
    " IT IS DEAD BECAUSE THE POPULATION WAS CONSTRUCTED BY EXCLUDING EVERY ELEMENT THAT COULD HAVE"
    " POPULATED THE OTHER ARM. The run measures it regardless -- a definition read from code and never"
    " checked against the data is how a lane convinces itself of a thing that is not there -- and the"
    " measurement's expected values are declared here first: fraction 1.0 on the evidence-kind and the"
    " how fields, in both the inhibiting and the non-inhibiting population, and 0 non-predicted"
    " inhibiting elements. The evidence-source-prefix field has NO such derivation and its value is"
    " genuinely unknown at this point."
)

COMPOSITION_AXIS_DEAD_RULE = (
    "THE RULE THAT KILLS THE COMPOSITION AXIS, fixed before the base rate is measured by this lane. For"
    " each of the three BASE_RATE_FIELDS, let f be the fraction of THE INHIBITING ELEMENTS OF THE"
    " ELIGIBLE MULTI-ELEMENT GENES carrying the field's most common value. The field is DEAD when"
    " f >= DEAD_FRACTION = 0.999, because a two-arm comparison on a field that is 99.9% one value cannot"
    " separate a tail of 1,330 from a body of 153,725 at any ratio worth stating. When ALL THREE fields"
    " are dead the COMPOSITION AXIS IS DEAD, no risk ratio is computed on any of them, and the honest"
    " answer is UNANSWERABLE_ON_THE_COMPOSITION_AXIS, declared complete in advance below."
)

UNANSWERABLE_ANSWER = (
    "UNANSWERABLE ON THE COMPOSITION AXIS -- DECLARED IN ADVANCE AS A COMPLETE ANSWER. If the base rate"
    " is at or near 100% predicted, then 'all 40 offenders are predicted deletion notes' and 'all 1,330"
    " above-bound elements are predicted deletion notes' are both re-statements of the layer's own"
    " composition and say nothing whatever about the tail. The composition question cannot distinguish"
    " anything on this axis, the comparison has no control arm inside the layer, and NO ENRICHMENT IS"
    " CLAIMED IN EITHER DIRECTION -- neither that the tail is the instrument's nor that it is not. This"
    " is a likely outcome, it is written here before the measurement, and IT IS NOT A FAILURE: it is the"
    " measurement of a question that the data cannot be asked. The surviving axes below are what the"
    " lane then reports, each with its own decision."
)

# ---------------------------------------------------------------------------------------------------
# THE AXES. Each one's question, its decision rule and both its wordings, all fixed before any of them
# is computed. `decides_the_outcome` says whether the axis may choose among the four OUTCOMES.
# ---------------------------------------------------------------------------------------------------

AXES: dict[str, dict[str, Any]] = {
    "A1_composition": {
        "question": (
            "Is the TAIL (UNSCALED) composed of predicted deletion notes in a way THE BODY (UNSCALED) is"
            " not, on any of the three BASE_RATE_FIELDS?"
        ),
        "rule": DEAD_BY_DEFINITION + " " + COMPOSITION_AXIS_DEAD_RULE + " " + ENRICHMENT_STATEMENT,
        "readings": {
            "enriched": (
                "the tail is enriched for {value} against the body on {field}: risk ratio {rr}, 95%"
                " interval [{lo}, {hi}], over {tail_n} tail and {body_n} body inhibiting elements of the"
                " eligible multi-element genes"
            ),
            "not_enriched": (
                "no enrichment detected on {field}: the tail's proportion of {value} against the body's"
                " gives risk ratio {rr}, 95% interval [{lo}, {hi}] spanning 1.0, over {tail_n} tail and"
                " {body_n} body inhibiting elements of the eligible multi-element genes"
            ),
            "dead": (
                "{field} carries no information in this population: {frac} of the {n} inhibiting"
                " elements of the eligible multi-element genes carry its most common value {value}, at"
                " or above DEAD_FRACTION {dead}, so no ratio is computed on it"
            ),
        },
        "decides_the_outcome": True,
    },
    "A2_gap_at_the_bound": {
        "question": (
            "Does the bound sit at a GAP in the population's inhibitor folds, or in the body of a"
            " continuous distribution? Over the sorted RHO of every inhibiting element of the eligible"
            " multi-element genes."
        ),
        "rule": (
            "TWO STATISTICS, both fixed here. (i) THE LARGEST CONSECUTIVE MULTIPLICATIVE GAP over all the"
            " sorted RHO values, r[i+1]/r[i], and the RHO at which it occurs: the bound is SEPARATED AT"
            " THE BOUND when that RHO lies within a factor of GAP_NEAR_BAND = 1.2 of"
            " BOUND_RHO_UNSCALED = 2.0, i.e. in [1.667, 2.400], and NOT AT THE BOUND otherwise, in which"
            " case the location of the largest gap is reported instead. (ii) THE WINDOWED COMPARISON: the"
            " largest consecutive multiplicative gap inside the STRADDLING WINDOW RHO in [1.5, 3.0]"
            " against the largest inside the CONTROL WINDOW RHO in [1.2, 1.5], both over the same"
            " population. A GAP AT THE BOUND when the straddling window's largest gap is at least"
            " GAP_RATIO = 2.0 times the control window's, NO GAP AT THE BOUND otherwise. The two"
            " statistics are reported separately and are not combined into one verdict."
        ),
        "readings": {
            "separated": (
                "the bound sits at a gap: the largest consecutive multiplicative gap over the {n}"
                " inhibiting elements of the eligible multi-element genes is {gap} at RHO {at}, within"
                " 1.2 of the bound 2.0"
            ),
            "not_separated": (
                "the bound does not sit at a gap: the largest consecutive multiplicative gap over the"
                " {n} inhibiting elements of the eligible multi-element genes is {gap} at RHO {at}, not"
                " within 1.2 of the bound 2.0"
            ),
            "windowed_gap": (
                "a gap at the bound on the windowed comparison: largest gap {straddling} in RHO [1.5,"
                " 3.0] against {control} in RHO [1.2, 1.5], ratio {ratio}, at or above 2.0"
            ),
            "windowed_no_gap": (
                "no gap at the bound on the windowed comparison: largest gap {straddling} in RHO [1.5,"
                " 3.0] against {control} in RHO [1.2, 1.5], ratio {ratio}, below 2.0"
            ),
        },
        "decides_the_outcome": False,
    },
    "A3_one_regime_or_two": {
        "question": (
            "Is the tail CONTINUOUS with the body, or a separate regime? The slope of the empirical"
            " survival function of RHO - 1 on log-log axes, below the bound and above it, over every"
            " inhibiting element of the eligible multi-element genes."
        ),
        "rule": (
            "WINDOW BELOW: RHO in [1.2, 2.0). WINDOW ABOVE: RHO in [2.0, 20.0]. In each window the slope"
            " is the ordinary least-squares slope of log10 S(r) against log10(r - 1), where S is the"
            " EMPIRICAL SURVIVAL over the WHOLE population (the fraction of all inhibiting elements of"
            " the eligible multi-element genes with RHO >= r) and each window element contributes one"
            " point. ONE CONTINUOUS REGIME when |slope_above / slope_below - 1| <= REGIME_TOLERANCE ="
            " 0.5, TWO REGIMES otherwise. Both windows and the tolerance are fixed here. THIS IS A"
            " DESCRIPTIVE SLOPE OF AN EMPIRICAL SURVIVAL CURVE. It is not a model, it is not fitted to"
            " any outcome, it changes no fit of the project and no fit change is proposed from it."
        ),
        "readings": {
            "one_regime": (
                "one continuous regime: the survival slope of RHO - 1 is {below} on RHO [1.2, 2.0) and"
                " {above} on RHO [2.0, 20.0], a ratio of {ratio}, within 0.5 of 1, over the {n}"
                " inhibiting elements of the eligible multi-element genes"
            ),
            "two_regimes": (
                "two regimes: the survival slope of RHO - 1 is {below} on RHO [1.2, 2.0) and {above} on"
                " RHO [2.0, 20.0], a ratio of {ratio}, further than 0.5 from 1, over the {n} inhibiting"
                " elements of the eligible multi-element genes"
            ),
        },
        "decides_the_outcome": False,
    },
    "A4_chromosome": {
        "question": (
            "Does the tail CONCENTRATE by chromosome? The tail rate of each of the 24 compiled programs"
            " against the genome-wide rate, over the inhibiting elements of the eligible multi-element"
            " genes."
        ),
        "rule": (
            "Per program c, k_c tail elements of n_c inhibiting elements. The genome-wide rate p is the"
            " pooled tail fraction. Only programs with n_c >= MIN_GROUP = 1000 are tested. CONCENTRATED"
            " when some tested program has BOTH an exact two-sided binomial p-value against p below"
            " ALPHA / 24 = 0.001 / 24 = 4.1667e-05 AND a rate ratio (k_c/n_c) / p at or above"
            " RATE_RATIO = 2.0. NOT CONCENTRATED otherwise. The max-over-min rate ratio among tested"
            " programs is reported either way. 24 is the number of programs and is the Bonferroni"
            " divisor whether or not every program is tested."
        ),
        "readings": {
            "concentrated": (
                "the tail concentrates by chromosome: {chrom} carries {k} of its {n} inhibiting elements"
                " above the bound, a rate of {rate} against the genome-wide {p} over the {total}"
                " inhibiting elements of the eligible multi-element genes, rate ratio {ratio}, exact"
                " two-sided binomial p {pval}, below 0.001/24"
            ),
            "not_concentrated": (
                "no concentration by chromosome detected: over the {tested} compiled programs with at"
                " least 1000 inhibiting elements of eligible multi-element genes, the tail rate runs"
                " {lo} to {hi} against the genome-wide {p}, a max-over-min ratio of {spread}, and no"
                " program clears both 0.001/24 and a rate ratio of 2.0"
            ),
        },
        "decides_the_outcome": True,
    },
    "A5_compiled_cell": {
        "question": (
            "Does the tail CONCENTRATE by COMPILED CELL -- the biosample label the predictor's note was"
            " produced in? This is the sharpest instrument axis available: the predictor was run per"
            " biosample, so a tail concentrated in a few labels far beyond their share is something the"
            " predictor does in those labels."
        ),
        "rule": (
            "Exactly A4's rule with the compiled cell label in place of the chromosome, and the"
            " Bonferroni divisor L equal to the NUMBER OF TESTED LABELS (those with n >= MIN_GROUP ="
            " 1000), which is reported. CONCENTRATED BY BIOSAMPLE when some tested label has both an"
            " exact two-sided binomial p below ALPHA / L = 0.001 / L and a rate ratio at or above 2.0;"
            " NOT CONCENTRATED BY BIOSAMPLE otherwise."
        ),
        "readings": {
            "concentrated": (
                "the tail concentrates by compiled cell: {cell} carries {k} of its {n} inhibiting"
                " elements above the bound, a rate of {rate} against the pooled {p} over the {total}"
                " inhibiting elements of the eligible multi-element genes, rate ratio {ratio}, exact"
                " two-sided binomial p {pval}, below 0.001/{labels}"
            ),
            "not_concentrated": (
                "no concentration by compiled cell detected: over the {tested} labels with at least 1000"
                " inhibiting elements of eligible multi-element genes, the tail rate runs {lo} to {hi}"
                " against the pooled {p}, a max-over-min ratio of {spread}, and no label clears both"
                " 0.001/{labels} and a rate ratio of 2.0"
            ),
        },
        "decides_the_outcome": True,
    },
    "A6_gene_inhibitor_count": {
        "question": (
            "Does the tail depend on HOW MANY INHIBITORS the element's gene has? Over the inhibiting"
            " elements of the eligible multi-element genes, each labelled with its own gene's number of"
            " inhibiting elements."
        ),
        "rule": (
            "The quartiles of that per-element label are taken over the whole population. The tail rate"
            " in the TOP quartile is compared with the tail rate in the BOTTOM quartile. DEPENDENT when"
            " the two rates differ by a factor of at least QUARTILE_RATIO = 2.0 in either direction,"
            " INDEPENDENT otherwise. The medians of the label in the tail and in the body are reported"
            " beside it and decide nothing. Ties at a quartile boundary go to the lower quartile, so the"
            " two groups are disjoint and their sizes are reported."
        ),
        "readings": {
            "dependent": (
                "the tail depends on the gene's inhibitor count: the tail rate is {top} in the top"
                " quartile ({top_n} elements) against {bottom} in the bottom quartile ({bottom_n}"
                " elements) of the inhibiting elements of the eligible multi-element genes, a factor of"
                " {ratio}, at or above 2.0"
            ),
            "independent": (
                "no dependence on the gene's inhibitor count detected: the tail rate is {top} in the top"
                " quartile ({top_n} elements) against {bottom} in the bottom quartile ({bottom_n}"
                " elements) of the inhibiting elements of the eligible multi-element genes, a factor of"
                " {ratio}, below 2.0"
            ),
        },
        "decides_the_outcome": False,
    },
    "A7_sharing": {
        "question": (
            "Are tail elements SHARED between genes -- one element inhibiting many -- or one per gene?"
        ),
        "rule": (
            "The number of DISTINCT eligible multi-element genes each element id inhibits is counted over"
            " the whole population. EXPECTED VALUE, DECLARED HERE AND NOT MEASURED FIRST: 0 element ids"
            " inhibit more than one gene, because `589b6ca` records `compiled_targets_1` as 440,377 of"
            " 440,377 and `compiled_targets_2_or_more` as 0 -- one element, one target, one rule. ONE PER"
            " GENE when the measured count of element ids under more than one eligible gene is 0, SHARED"
            " when it is not. EITHER WAY THIS AXIS DECIDES NOTHING ABOUT THE TAIL: at 0 it is a property"
            " of the compiler's emission and holds identically in the tail and in the body, and a"
            " non-zero count would falsify `589b6ca`'s figure rather than say anything about the tail."
        ),
        "readings": {
            "one_per_gene": (
                "one element, one gene: 0 of the {elements} distinct element ids among the inhibiting"
                " elements of the eligible multi-element genes inhibit more than one of those genes, as"
                " 589b6ca's compiled_targets_2_or_more = 0 requires, so sharing is a property of the"
                " compiler's emission and separates nothing in the tail"
            ),
            "shared": (
                "elements are shared: {shared} of the {elements} distinct element ids among the"
                " inhibiting elements of the eligible multi-element genes inhibit more than one of those"
                " genes, which contradicts 589b6ca's compiled_targets_2_or_more = 0 and is a finding"
                " about the compiler, not about the tail"
            ),
        },
        "decides_the_outcome": False,
    },
    "A8_tautology": {
        "question": (
            "Is the tail the set of elements where the predicted deletion note's OWN EFFECT SIZE is"
            " largest -- which would make it the predictor's extreme tail rather than a separate"
            " population?"
        ),
        "rule": (
            "REGISTERED AS TAUTOLOGICAL BEFORE IT IS COMPUTED, from the code and not from the data."
            " `bridge.observation_of` returns RHO = 2 ** (the note's log2 fold change) for every"
            " predicted rule whose note carries an effect, and IMPLIED_STRENGTH(RHO) = (1 - 1/RHO)/H is"
            " strictly increasing in RHO on RHO > 1. So the above-bound set is EXACTLY the set whose"
            " note states a log2 fold above log2(BOUND_RHO_UNSCALED) = 1.0, and at the scale above"
            " log2(BOUND_RHO_AT_THE_SCALE) = 2.501. THE VERIFICATION, with its expected values: the"
            " count of inhibiting elements where (IMPLIED_STRENGTH > 1 + STRENGTH_TOLERANCE) differs"
            " from (RHO > BOUND_RHO_UNSCALED) must be 0; the minimum RHO in the tail must exceed"
            " BOUND_RHO_UNSCALED and the maximum RHO in the body must not. The axis is TAUTOLOGICAL when"
            " those hold and BROKEN when they do not. IT CHOOSES NO OUTCOME: a tautology cannot be"
            " evidence for one."
        ),
        "readings": {
            "tautological": (
                "tautological, and that is the answer to this axis rather than a finding from it: the"
                " above-bound tail IS the set of inhibiting elements of the eligible multi-element genes"
                " whose predicted deletion note states a log2 fold change above 1.0 (above 2.501 at"
                " M2h's scale), because RHO is 2 to that fold and IMPLIED_STRENGTH is strictly"
                " increasing in RHO. {mismatches} of {n} elements disagree between the strength bound"
                " and the fold bound. The tail is the upper tail of the predictor's own stated effect"
                " sizes by construction and is not a separately identified population"
            ),
            "broken": (
                "the identity does not hold: {mismatches} of {n} inhibiting elements of the eligible"
                " multi-element genes disagree between IMPLIED_STRENGTH > 1 and RHO > 2, which is a"
                " defect in this lane's arithmetic or in the bridge and must be resolved before any"
                " reading of the tail is quoted"
            ),
        },
        "decides_the_outcome": False,
    },
    "A9_control_arm": {
        "question": (
            "Is there an EXPERIMENTAL arm inside this population at all -- any inhibiting element of an"
            " eligible multi-element gene whose rule is not predicted?"
        ),
        "rule": (
            "Over the same 24 gated programs: the compiled rule lines by `bridge._kind`; of the"
            " regulatory rules whose target is an eligible multi-element gene, the count by kind; of"
            " those, how many yield a readable `bridge.Observation`; and of THOSE, how many state"
            " RHO > 1. NO CONTROL ARM when the last count is 0 for every kind other than `predicted`,"
            " A CONTROL ARM EXISTS otherwise, with its size. EXPECTED VALUE, DECLARED HERE BEFORE THE"
            " RUN AND WITH ITS REASON: 0, BY THE POPULATION'S DEFINITION, because"
            " `combine2.multi_rule_genes` drops any gene with a rule that is not predicted"
            " (`a_rule_not_predicted`) before an inhibitor of it is ever counted. A non-zero count would"
            " therefore falsify this lane's reading of that function, not reveal a control arm. At 0 the"
            " instrument-versus-biology contrast has nothing to be drawn against inside the layer, and"
            " that is reported as the reason the composition axis cannot be answered rather than as a"
            " result about the tail."
        ),
        "readings": {
            "no_control_arm": (
                "there is no control arm inside the layer: of the {rules} compiled rule lines over the"
                " 24 programs, {experimental} are experimental, and 0 of the {inhibiting} inhibiting"
                " elements of the eligible multi-element genes comes from a rule that is not predicted."
                " The instrument-versus-biology contrast has nothing to be drawn against here, so the"
                " composition axis cannot be answered in either direction on this layer"
            ),
            "control_arm": (
                "a control arm exists: {n} of the {inhibiting} inhibiting elements of the eligible"
                " multi-element genes come from rules that are not predicted, by kind {by_kind}, and the"
                " composition axis can be read against them"
            ),
        },
        "decides_the_outcome": False,
    },
}

#: Axes named in advance as NOT COMPUTED, each with the reason, so that an axis the brief asked for
#: cannot be quietly dropped.
NOT_COMPUTED = {
    "element_structural_class": (
        "THE ELEMENT'S STRUCTURAL CLASS (promoter / enhancer / region / domain) IS NOT COMPUTED. The"
        " compiled rule line carries no class field: the class would have to be read from the element,"
        " region and domain BLOCKS' own properties, which is a separate traversal of the programs this"
        " lane does not declare as an input and does not perform. The compiled CELL label, which the"
        " rule does carry, is taken as the available class and is axis A5. Declared here before the run"
        " so that its absence is a stated limit and not an omission."
    ),
    "a_fit_change": (
        "NO FIT CHANGE IS PROPOSED, EXAMINED OR COSTED. This lane reads a bound that follows from"
        " H = 0.5 by algebra and asks what sits above it. Whether that bound should move is not this"
        " lane's question and no number here should be quoted towards it."
    ),
}

# ---------------------------------------------------------------------------------------------------
# THE CONSTANTS every decision rule above refers to, in one place and fixed before any computation.
# ---------------------------------------------------------------------------------------------------

GAP_NEAR_BAND = 1.2
GAP_STRADDLING_WINDOW = (1.5, 3.0)
GAP_CONTROL_WINDOW = (1.2, 1.5)
GAP_RATIO = 2.0
REGIME_WINDOW_BELOW = (1.2, 2.0)
REGIME_WINDOW_ABOVE = (2.0, 20.0)
REGIME_TOLERANCE = 0.5
MIN_GROUP = 1000
ALPHA = 0.001
RATE_RATIO = 2.0
QUARTILE_RATIO = 2.0
N_PROGRAMS = 24

CONVENTIONS = {
    "scale": SCALE,
    "scale_is_refitted": False,
    "h_intact": H,
    "strength_tolerance": STRENGTH_TOLERANCE,
    "bound_rho_unscaled": BOUND_RHO_UNSCALED,
    "bound_rho_at_the_scale": BOUND_RHO_AT_THE_SCALE,
    "dead_fraction": DEAD_FRACTION,
    "gap_near_band": GAP_NEAR_BAND,
    "gap_straddling_window": list(GAP_STRADDLING_WINDOW),
    "gap_control_window": list(GAP_CONTROL_WINDOW),
    "gap_ratio": GAP_RATIO,
    "regime_window_below": list(REGIME_WINDOW_BELOW),
    "regime_window_above": list(REGIME_WINDOW_ABOVE),
    "regime_tolerance": REGIME_TOLERANCE,
    "min_group": MIN_GROUP,
    "alpha": ALPHA,
    "rate_ratio": RATE_RATIO,
    "quartile_ratio": QUARTILE_RATIO,
    "n_programs": N_PROGRAMS,
    "rna_read": "none",
    "alphagenome_requests": 0,
}

# ---------------------------------------------------------------------------------------------------
# THE FOUR OUTCOMES, all written before any tail number is read, with the rule that selects one.
# ---------------------------------------------------------------------------------------------------

OUTCOMES = {
    "a_the_tail_is_the_instruments": (
        "THE TAIL IS THE INSTRUMENT'S. Licensed by one of two things and nothing else: a composition"
        " field that is NOT dead on which the tail's risk ratio against the body has a 95% Katz interval"
        " EXCLUDING 1.0; or a CONCENTRATED verdict on a concentration axis (A4 chromosome, A5 compiled"
        " cell) with the other concentration axis agreeing. Wording: 'the above-bound inhibitor tail is"
        " a property of the prediction instrument: <the licensing axis's own reading, verbatim>.'"
    ),
    "b_not_distinguishable": (
        "THE TAIL IS NOT DISTINGUISHABLE from the population on these axes. Taken when the composition"
        " axis is ALIVE and its intervals all span 1.0, and both concentration axes return NOT"
        " CONCENTRATED. Wording: 'no property distinguishing the above-bound inhibitor tail from the"
        " body of the same population was detected on the registered axes: <each axis's null reading>."
        " No difference detected is not no difference.'"
    ),
    "c_unanswerable_on_the_composition_axis": (
        "UNANSWERABLE ON THE COMPOSITION AXIS because of the base rate, with the surviving axes"
        " reported. Taken when all three BASE_RATE_FIELDS are DEAD by COMPOSITION_AXIS_DEAD_RULE and"
        " neither concentration axis is CONCENTRATED. Wording: UNANSWERABLE_ANSWER, verbatim, followed"
        " by each surviving axis's own registered reading."
    ),
    "d_mixed_no_verdict_taken": (
        "MIXED / NO VERDICT TAKEN, the lane-regionscale2 shape. Taken when the two concentration axes"
        " DISAGREE -- one CONCENTRATED and the other NOT -- because they are the same kind of evidence"
        " about the same claim and this registration forbids preferring one over the other. Wording:"
        " 'two registered readings of equal standing disagree -- <A4's reading> against <A5's reading> --"
        " and this registration forbids a verdict: none is taken.'"
    ),
}

OUTCOME_RULE = (
    "THE OUTCOME IS SELECTED BY THIS RULE AND NOTHING ELSE, in this order, all of it fixed before any"
    " tail number is read. (1) If A4 and A5 disagree (exactly one of them CONCENTRATED) the outcome is"
    " (d) MIXED / NO VERDICT TAKEN. (2) Else if all three BASE_RATE_FIELDS are DEAD and neither A4 nor"
    " A5 is CONCENTRATED, the outcome is (c) UNANSWERABLE ON THE COMPOSITION AXIS. (3) Else if A4 and A5"
    " are both CONCENTRATED, or some live composition field's 95% Katz interval excludes 1.0, the"
    " outcome is (a) THE TAIL IS THE INSTRUMENT'S. (4) Else the outcome is (b) NOT DISTINGUISHABLE."
    " AXES A2, A3, A6, A7, A8 AND A9 ARE REPORTED IN EVERY OUTCOME AND CHOOSE NONE OF THEM"
    " (`decides_the_outcome` is False on each): in particular the tautology A8 and the empty control arm"
    " A9, which are the two most quotable things this lane will say, are forbidden from manufacturing"
    " outcome (a)."
)

# ---------------------------------------------------------------------------------------------------
# THE INHERITED FIGURES, pinned AS THE SAME OBJECTS they were read from so that this lane's summary
# cannot drift from the committed record. (label, source file, json path, expected value).
# ---------------------------------------------------------------------------------------------------

CENSUS_RESULT = "data/results/inhibitor_census_result.json"
EMIT7G_RESULT = "data/results/emit7g_census.json"

PINNED: tuple[tuple[str, str, tuple[str, ...], Any], ...] = (
    ("the_39", CENSUS_RESULT, ("Q6_outside_the_39", "populations", "the_39"), 39),
    ("all_eligible_genes", CENSUS_RESULT, ("Q6_outside_the_39", "populations", "all_eligible"), 17983),
    ("the_residue", CENSUS_RESULT, ("Q6_outside_the_39", "populations", "the_residue"), 3708),
    (
        "inhibiting_elements",
        CENSUS_RESULT,
        ("Q7_threshold_or_tail", "inhibiting_elements_of_all_eligible_genes"),
        155055,
    ),
    (
        "tail_elements_unscaled",
        CENSUS_RESULT,
        ("Q7_threshold_or_tail", "unscaled", "above_the_bound"),
        1330,
    ),
    (
        "tail_elements_at_the_scale",
        CENSUS_RESULT,
        ("Q7_threshold_or_tail", "at_the_scale", "above_the_bound"),
        192,
    ),
    (
        "just_below_the_bound_at_the_scale",
        CENSUS_RESULT,
        ("Q7_threshold_or_tail", "at_the_scale", "just_below_the_bound"),
        69,
    ),
    ("max_rho", CENSUS_RESULT, ("Q7_threshold_or_tail", "rho", "max"), 145.00913731786218),
    ("median_rho", CENSUS_RESULT, ("Q7_threshold_or_tail", "rho", "quantiles", "0.5"), 1.1227214220736084),
    (
        "median_strength_at_the_scale",
        CENSUS_RESULT,
        ("Q7_threshold_or_tail", "strength_at_the_scale", "quantiles", "0.5"),
        0.05131839300792285,
    ),
    (
        "genes_above_the_bound_unscaled_all",
        CENSUS_RESULT,
        ("Q6_outside_the_39", "above_the_bound_unscaled", "all_eligible", "genes"),
        1199,
    ),
    (
        "genes_above_the_bound_unscaled_outside_the_39",
        CENSUS_RESULT,
        ("Q6_outside_the_39", "above_the_bound_unscaled", "all_eligible", "genes_outside_the_39"),
        1160,
    ),
    (
        "genes_above_the_bound_at_the_scale_all",
        CENSUS_RESULT,
        ("Q6_outside_the_39", "above_the_bound_at_the_scale", "all_eligible", "genes"),
        188,
    ),
    (
        "genes_above_the_bound_at_the_scale_outside_the_39",
        CENSUS_RESULT,
        ("Q6_outside_the_39", "above_the_bound_at_the_scale", "all_eligible", "genes_outside_the_39"),
        149,
    ),
    (
        "residue_at_the_scale_genes_outside_the_39",
        CENSUS_RESULT,
        ("Q6_outside_the_39", "above_the_bound_at_the_scale", "the_residue", "genes_outside_the_39"),
        0,
    ),
    (
        "residue_at_the_scale_elements",
        CENSUS_RESULT,
        ("Q6_outside_the_39", "above_the_bound_at_the_scale", "the_residue", "elements"),
        40,
    ),
    (
        "offending_elements_of_the_39",
        CENSUS_RESULT,
        ("Q5_where_the_rules_came_from", "offending_elements"),
        40,
    ),
    (
        "offenders_by_evidence_kind",
        CENSUS_RESULT,
        ("Q5_where_the_rules_came_from", "by_evidence_kind"),
        {"predicted": 40},
    ),
    (
        "inherited_base_rate_by_evidence_kind",
        CENSUS_RESULT,
        (
            "Q5_where_the_rules_came_from",
            "background_all_inhibiting_elements_of_all_eligible_genes",
            "by_evidence_kind",
        ),
        {"predicted": 155055},
    ),
    ("census_branch", CENSUS_RESULT, ("branch",), "a_separated_tail"),
    ("rules_total", EMIT7G_RESULT, ("genome_wide", "now", "rules_total"), 440589),
    ("rules_predicted", EMIT7G_RESULT, ("genome_wide", "now", "rules_predicted"), 440377),
    ("rules_experimental", EMIT7G_RESULT, ("genome_wide", "now", "rules_experimental"), 212),
)

#: The inherited readings, carried WORD FOR WORD and never strengthened. Each says where it comes from,
#: because one of them is a field of a committed result and one is not.
INHERITED_READINGS = {
    "the_census_implication_committed_field": (
        'data/results/inhibitor_census_result.json `implication`, verbatim: "A SEPARATED TAIL: the'
        " offending folds are not continuous with the bound -- 69 inhibiting elements sit within a"
        " factor of 1.2 just below it against 192 above it -- so the 39 are refused on observations that"
        " stand apart from the population's inhibitor folds rather than on where a bound was drawn. 1160"
        ' eligible genes outside the 39 carry an inhibitor above the bound unscaled."'
    ),
    "the_lanes_reported_implication_not_a_committed_field": (
        "lane-inhibitor's reported implication, carried verbatim from this lane's brief and NOT found as"
        " a string in any committed file of this repository (searched; the committed `implication` field"
        " is the sentence above): \"the compiled layer's inhibitor folds are tiny with a thin tail to"
        " 145; any fold past 2 maps to a strength the fit cannot hold. That tail is the layer's, not the"
        " residue's -- the 39 are its intersection with the residue.\" It is carried as reported and is"
        " not strengthened, and the distinction between a reported reading and a committed field is"
        " stated rather than smoothed over."
    ),
    "the_branch": (
        "lane-inhibitor's registered branch word, carried unchanged and unstrengthened: `a_separated_tail`."
    ),
    "the_row_correction": (
        "THE COMMITTED ROW'S PROSE AND THE RECORD'S FIELDS DISAGREE AND THE FIELDS WIN. `7d44158`'s row"
        ' reads "149 genes OUTSIDE the 39 carry one at the scale and 1,199 unscaled", which pairs an'
        " outside-the-39 figure with an in-all figure. The correct pairing, from"
        " Q6_outside_the_39: 1,160 genes outside the 39 of 1,199 in all, unscaled (1,330 elements); 149"
        " outside the 39 of 188 in all, at the scale (192 elements). This lane quotes the fields, never"
        " the row's prose, and the correction is boarded to lane-inhibitor's row rather than applied by"
        " this lane."
    ),
    "the_population_clause": (
        "Every one of the figures above is a count over THE ELIGIBLE MULTI-ELEMENT GENES of the 24"
        " compiled non-coding programs, or over their inhibiting elements. A figure quoted without that"
        " clause is not ready to quote."
    ),
}

#: The input groups this lane declares, each with the member count DECLARED HERE as this lane's own
#: statement rather than discovered by a glob: a declaration that a filesystem can change is not a
#: declaration. The run REFUSES when a group's member count differs from the number declared here.
INPUT_GROUPS = {
    "the 24 compiled non-coding programs": 24,
}

#: Every reader this lane's run uses, declared. Nothing else under data/ is opened; the run asserts it
#: through `manifest.traced_inputs` and records the reconciliation in the result.
DECLARED_READERS = (
    "data/knowledge/compiled/noncoding_chr{1..22,X,Y}.bio -- read as text, gated by"
    " `compile.gate_program` and parsed by `lang.parser.parse`; the ONLY source of every number this"
    " lane produces",
    f"{CENSUS_RESULT} -- read for the PINNED inherited figures and for nothing else; no count of this"
    " lane is taken from it",
    f"{EMIT7G_RESULT} -- read for the three PINNED rule-line counts and for nothing else",
)

LIMITATIONS = (
    "A PREDICTOR-SIDE CENSUS. Every number is arithmetic on what the compiled rules themselves state."
    " Nothing here is a measurement of biology, and no reading of it is evidence about a real repressor.",
    "NOTHING IS FITTED AND NO FIT CHANGE IS PROPOSED. The bound follows from H = 0.5 by algebra and the"
    " scale is M2h's constant.",
    "THE COMPOSITION AXIS HAS NO CONTROL ARM IF A9 RETURNS 0. A base rate at 100% predicted means the"
    " question cannot be asked of this layer -- not that the tail is biological, and not that it is not.",
    "A NEGATIVE ON A CONCENTRATION AXIS IS NOT AN ABSENCE. MIN_GROUP = 1000 excludes small programs and"
    " small labels from the test, and a concentration confined to them is invisible here.",
    "THE SURVIVAL SLOPE OF A3 IS DESCRIPTIVE. Two windows and one tolerance were chosen in advance;"
    " another pair of windows could give another ratio, and the slope is not a model of anything.",
    "THE ELEMENT'S STRUCTURAL CLASS IS NOT COMPUTED (NOT_COMPUTED['element_structural_class']). The"
    " compiled cell label stands in for it and is not the same thing.",
    "THE TAIL IS DEFINED BY A BOUND THAT IS A PROPERTY OF THE FIT'S ALGEBRA. Any statement about 'the"
    " tail' is a statement about elements above RHO = 2, not about a population nature separated at 2.",
    "NO RNA VALUE IS READ AND NOTHING ABOUT RNA IS CLAIMED.",
    "0 MODEL REQUESTS, NO KEY READ, NO NETWORK, $0.00.",
)


# ---------------------------------------------------------------------------------------------------
# The arithmetic each rule above names. Pure Python: no scipy, no numpy, so this module imports in the
# bare closure and a missing extra cannot silently change a verdict.
# ---------------------------------------------------------------------------------------------------


def implied_strength(rho: float) -> float:
    """The quantity the fit refused: s = (1 - 1/RHO) / H, the multiplicative form's inhibitor strength."""
    return (1.0 - 1.0 / rho) / H


def above_bound(s: float) -> bool:
    """`combine2._strength_ok`'s refusal, with the bridge's own slack."""
    return bool(s > 1.0 + STRENGTH_TOLERANCE)


def katz_interval(a: int, n_a: int, b: int, n_b: int) -> dict[str, Any]:
    """The risk ratio (a/n_a)/(b/n_b) with its 95% Katz log interval, or why it is DEGENERATE.

    Degenerate when either arm is empty or either proportion is exactly 0 or exactly 1: the log variance
    1/a - 1/n_a + 1/b - 1/n_b is then undefined or zero, and ENRICHMENT_STATEMENT forbids quoting a ratio
    there. That is the expected case on a field that is 100% one value and it is reported, not patched.
    """
    if n_a <= 0 or n_b <= 0:
        return {"degenerate": True, "why": "an empty arm", "ratio": None, "lo": None, "hi": None}
    pa, pb = a / n_a, b / n_b
    if a in (0, n_a) or b in (0, n_b):
        return {
            "degenerate": True,
            "why": f"a proportion of exactly {'0' if 0 in (a, b) else '1'}: the Katz log variance is"
            " undefined, so no ratio is quoted",
            "ratio": (pa / pb) if pb > 0 else None,
            "lo": None,
            "hi": None,
            "p_tail": pa,
            "p_body": pb,
        }
    rr = pa / pb
    se = math.sqrt(1.0 / a - 1.0 / n_a + 1.0 / b - 1.0 / n_b)
    return {
        "degenerate": False,
        "why": "",
        "ratio": rr,
        "lo": math.exp(math.log(rr) - 1.959963984540054 * se),
        "hi": math.exp(math.log(rr) + 1.959963984540054 * se),
        "p_tail": pa,
        "p_body": pb,
    }


def _log_binom_pmf(k: int, n: int, p: float) -> float:
    if p <= 0.0:
        return 0.0 if k == 0 else -math.inf
    if p >= 1.0:
        return 0.0 if k == n else -math.inf
    return (
        math.lgamma(n + 1)
        - math.lgamma(k + 1)
        - math.lgamma(n - k + 1)
        + k * math.log(p)
        + (n - k) * math.log1p(-p)
    )


def binom_two_sided(k: int, n: int, p: float, rel_tol: float = 1e-7) -> float:
    """The exact two-sided binomial p-value by the density method: the total probability of every
    outcome no more likely than the observed one. Pure Python, so no extra can change it."""
    if n <= 0:
        return 1.0
    if p <= 0.0:
        return 1.0 if k == 0 else 0.0
    if p >= 1.0:
        return 1.0 if k == n else 0.0
    obs = _log_binom_pmf(k, n, p)
    cut = obs + math.log1p(rel_tol)
    total = 0.0
    for i in range(n + 1):
        lp = _log_binom_pmf(i, n, p)
        if lp <= cut:
            total += math.exp(lp)
    return min(1.0, total)


def largest_multiplicative_gap(
    xs: list[float], lo: float | None = None, hi: float | None = None
) -> dict[str, Any]:
    """The largest r[i+1]/r[i] over the sorted values (optionally restricted to [lo, hi]), and where.

    `at` is the LOWER of the two values the gap sits between, so "the gap at RHO 1.98" means the next
    value above 1.98 is `gap` times it.
    """
    ys = sorted(x for x in xs if (lo is None or x >= lo) and (hi is None or x <= hi))
    if len(ys) < 2:
        return {"n": len(ys), "gap": None, "at": None, "next": None}
    best, at = 1.0, ys[0]
    nxt = ys[1]
    for a, b in zip(ys, ys[1:], strict=False):
        if a > 0 and b / a > best:
            best, at, nxt = b / a, a, b
    return {"n": len(ys), "gap": best, "at": at, "next": nxt}


def survival_slope(xs: list[float], lo: float, hi: float, include_hi: bool) -> dict[str, Any]:
    """OLS slope of log10 S(r) against log10(r - 1) over the window, S the survival on the WHOLE `xs`.

    S(r) is the fraction of ALL of `xs` at or above r, so the two windows read the same curve; a window
    element with r <= 1 or S(r) == 0 cannot be placed on log axes and is counted as `dropped`.
    """
    ys = sorted(xs)
    n = len(ys)
    pts: list[tuple[float, float]] = []
    dropped = 0
    for r in ys:
        if r < lo or (r > hi if include_hi else r >= hi):
            continue
        s = (n - bisect.bisect_left(ys, r)) / n
        if r <= 1.0 or s <= 0.0:
            dropped += 1
            continue
        pts.append((math.log10(r - 1.0), math.log10(s)))
    if len(pts) < 2:
        return {"n": len(pts), "dropped": dropped, "slope": None}
    mx = sum(p[0] for p in pts) / len(pts)
    my = sum(p[1] for p in pts) / len(pts)
    sxx = sum((p[0] - mx) ** 2 for p in pts)
    sxy = sum((p[0] - mx) * (p[1] - my) for p in pts)
    return {
        "n": len(pts),
        "dropped": dropped,
        "slope": (sxy / sxx) if sxx > 0 else None,
        "window": [lo, hi],
        "window_includes_hi": include_hi,
    }


def quantile(xs: list[float], q: float) -> float:
    """The lower-interpolated quantile of a non-empty list, sorted here so callers cannot pass an
    unsorted list by mistake."""
    ys = sorted(xs)
    if not ys:
        raise ValueError("quantile of an empty list")
    if q <= 0:
        return ys[0]
    if q >= 1:
        return ys[-1]
    i = q * (len(ys) - 1)
    lo = int(math.floor(i))
    hi = min(lo + 1, len(ys) - 1)
    return ys[lo] + (i - lo) * (ys[hi] - ys[lo])


def dead_fields(base_rate: dict[str, Any]) -> dict[str, Any]:
    """COMPOSITION_AXIS_DEAD_RULE applied to a base-rate block: per field, its most common value, that
    value's fraction of THE INHIBITING ELEMENTS, and whether the field is dead."""
    out: dict[str, Any] = {}
    for field in BASE_RATE_FIELDS:
        counts = base_rate["inhibiting"][field]
        n = sum(counts.values())
        value, k = max(counts.items(), key=lambda kv: kv[1]) if counts else ("", 0)
        frac = (k / n) if n else 0.0
        out[field] = {
            "n": n,
            "most_common_value": value,
            "count": k,
            "fraction": frac,
            "dead": bool(frac >= DEAD_FRACTION),
        }
    out["composition_axis_is_dead"] = all(out[f]["dead"] for f in BASE_RATE_FIELDS)
    return out


def select_outcome(composition_is_dead: bool, a4: str, a5: str, any_live_interval_excludes_1: bool) -> str:
    """OUTCOME_RULE, as code, so the outcome cannot be chosen by a sentence written afterwards."""
    c4, c5 = a4 == "concentrated", a5 == "concentrated"
    if c4 != c5:
        return "d_mixed_no_verdict_taken"
    if composition_is_dead and not (c4 or c5):
        return "c_unanswerable_on_the_composition_axis"
    if (c4 and c5) or any_live_interval_excludes_1:
        return "a_the_tail_is_the_instruments"
    return "b_not_distinguishable"


def reading(axis: str, which: str, **numbers: Any) -> str:
    """One axis's REGISTERED wording, filled with numbers. The wording cannot be written later: it is
    `AXES[axis]['readings'][which]` and this function only substitutes."""
    return AXES[axis]["readings"][which].format(**numbers)
