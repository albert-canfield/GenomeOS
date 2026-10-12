# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""lane-regionscale (ROADMAP area I, 2026-10-12): does a PER-REGION magnitude correction predict
held-out RNA better than M2h's single global +0.2146, on regions held out of the fit?

M2h (`ac5580b`) fitted ONE global magnitude correction, +0.2145911381525326, and that single number
is what makes the compiled rules read better than no rule at all: held-out Delta MSE
+0.0433912108101358, about 2% of the variance of measured log2(TPM+1) across the thirteen reader
cells. It survives a fresh measurement with nothing fitted (M2i, `f6e71c6`), a chromosome-held-out
split with the scale refitted (lane-scalesplit, `ec9cfd4`) and a 0-difference rebuild.

lane-scalehet (`40b0a8d`) then found the per-chromosome fitted scales NOT CONSISTENT WITH A SINGLE
GLOBAL SCALE: Q = 166.048, p_Q = 0.003996. lane-hetblock2 (`d9e3ddd`) re-tested that against a 1-Mb
locus-block null and IT SURVIVED: p_Q = 0.008991, p_R = 0.025974, the 5-Mb null agreeing at
p_Q = 0.014985, the block unit widening the null's median Q by only 4.8%. THE HONEST OTHER HALF,
CARRIED HERE: at the stricter 0.0167 the primary Q clears and the secondary R does NOT.

Both lanes deliberately PROPOSED NO REGIONAL MODEL, each recording that fitting per-region scales is
"a separate registration and a separate decision". THIS MODULE IS THAT REGISTRATION, and it asks the
one question those lanes did not: a variation that is REAL but UNUSABLE is a different thing from one
that IMPROVES PREDICTION, and only a held-out comparison can tell them apart.

THE WHOLE LANE IS THE HELD-OUT DESIGN, so it is stated here before anything else. A per-region scale
cannot be fitted by leaving the whole region out -- a region's scale needs the region's own genes, so
leave-one-chromosome-out is not a design for this estimator but a refusal of it. The resolution
registered here is a WITHIN-REGION SPLIT on M2h's OWN split rule, carried unchanged: for each
admitted chromosome, that chromosome's scale is fitted on ITS TRAINING GENES ALONE and applied to ITS
HELD-OUT GENES, which contributed nothing to it. The global scale is fitted on ALL training genes --
which is exactly M2h's +0.2146 -- and the SAME held-out genes are scored with it. Under both
predictors, NO HELD-OUT GENE CONTRIBUTES TO THE SCALE THAT SCORES IT, and the two predictors are
scored on the IDENTICAL rows, so the difference is attributable to the regional unit and to nothing
else. The design CAN be made honest at these gene counts, and the counts that say so are gene counts,
which depend on no measured value (`GENE_COUNTS_ARE_NOT_MEASUREMENTS`).

PRIOR EXPOSURE is the weakest position of any lane in this series and is stated TWICE, in
`PRIOR_EXPOSURE` and again in `PRIOR_EXPOSURE_CLAUSE`, which every reading carries.

It re-runs no stage, opens no compiled program, no narrowPeak file and no RNA quantification, and
edits no file of magcal's, scalesplit's, scalehet's or hetblock's. It reads 7h's committed genes
table, pinned by sha256, through magcal's own functions called READ-ONLY.

In the order written and committed: this code, the harness, the registration, then one run. No model
request, no key, no money, 0 AlphaGenome requests.
"""

from __future__ import annotations

from typing import Any

from genomeos.attribution import magcal as mc
from genomeos.attribution import scalehet as sh

REGISTERED = "2026-10-12"
LANE = "lane-regionscale"

#: the one composition form carried, 7h's, the one M2h's global scale was fitted for
FORM = "multiplicative"
FORM_WORDS = mc.FORM_WORDS[FORM]
#: the committed result the P and D are read from, pinned by sha256 by the harness and the run
SOURCE = "data/results/multrna_result.json"

# ---- what is carried verbatim, and from whom -----------------------------------------------------
CELLS = mc.CELLS
MEASURED = mc.MEASURED
NO_RULE_BASELINE = mc.NO_RULE_BASELINE
#: magcal's own gene floor, taken UNCHANGED and applied to EACH HALF of a region (see ADMISSION)
MIN_GENES = mc.MIN_GENES
#: scalehet's own floor on the number of regions, taken UNCHANGED
MIN_REGIONS = sh.MIN_CHROMOSOMES
#: magcal's own significance level, taken UNCHANGED (and deliberately: see ALPHA)
CONTROL_ALPHA = mc.CONTROL_ALPHA
#: magcal's own interval width, taken UNCHANGED, for comparability (see INTERVAL)
FORM_INTERVAL = mc.FORM_INTERVAL
FORM_INTERVAL_WORDS = mc.FORM_INTERVAL_WORDS
#: the nominal 95% width, reported BESIDE the registered one so the WIDTH cannot decide (see INTERVAL)
NOMINAL_INTERVAL = (2.5, 97.5)
NOMINAL_INTERVAL_WORDS = "95% CI"
#: magcal's own bootstrap, resamples and seed taken UNCHANGED
BOOTSTRAP = dict(mc.BOOTSTRAP)
TOLERANCE = mc.TOLERANCE
M2H_TOLERANCE = sh.M2H_TOLERANCE
#: M2h's committed gene-level numbers, as scalehet pinned them, re-pinned here and not re-decided
M2H_COMMITTED = dict(sh.M2H_COMMITTED)
#: 7h's committed table-level numbers, as scalehet pinned them
SEVEN_H_COMMITTED = dict(sh.SEVEN_H_COMMITTED)
#: M2h's global scale, the thing a per-region model must beat, quoted and never refitted as a target
GLOBAL_SCALE_COMMITTED = M2H_COMMITTED["scale"]

CARRIED = {
    "endpoint": MEASURED,
    "cells": list(CELLS),
    "no_rule_baseline": NO_RULE_BASELINE,
    "form": FORM,
    "split_rule": "magcal.SPLIT, sha256 of the gene symbol, odd first 32 bits held out -- UNCHANGED",
    "fit": "magcal.FIT, the one closed-form scale minimising squared error over a group's pairs",
    "metric": "magcal.METRIC, Delta MSE against the no-rule baseline, over gene-cell pairs",
    "bootstrap": "magcal.BOOTSTRAP, its resamples and its seed",
    "min_genes": MIN_GENES,
    "min_regions": "scalehet.MIN_CHROMOSOMES",
    "interval": "magcal.FORM_INTERVAL, 98.33%",
    "alpha": "magcal.CONTROL_ALPHA, 0.05",
    "regional_unit": "scalehet's and hetblock2's, the chromosome",
    "what_is_new_here": (
        "the HELD-OUT COMPARISON of a per-region predictor against the global one. Nothing else of"
        " M2h's, 7h's, scalehet's or hetblock2's design moves"
    ),
}

# ---- the question --------------------------------------------------------------------------------
QUESTION = (
    "Does a PER-REGION magnitude correction predict held-out measured RNA better than M2h's SINGLE"
    " GLOBAL scale +0.2145911381525326, when every region's scale is fitted WITHOUT the genes it is"
    " then scored on? lane-scalehet rejected a single global scale as a description of the"
    " per-chromosome fitted scales (Q = 166.048, p_Q = 0.003996) and lane-hetblock2 showed that"
    " rejection survives a 1-Mb locus-block null (p_Q = 0.008991), and BOTH deliberately proposed no"
    " regional model. A dispersion test says a quantity VARIES; it does not say that modelling the"
    " variation PREDICTS BETTER, because a real per-region difference estimated from a few hundred"
    " genes can be swamped by its own estimation error. This lane asks the prediction question and"
    " only the prediction question: on the SAME held-out genes, scored by two predictors neither of"
    " which saw them, is the per-chromosome magnitude closer to measured log2(TPM+1) across the"
    " thirteen reader cells than the one number M2h fitted for the whole genome?"
)
WHY_NOT_ALREADY_ANSWERED = (
    "scalehet and hetblock2 each fitted per-chromosome scales and each read them IN SAMPLE, as a"
    " dispersion about a global scale, against a permutation null -- and each said so and each refused"
    " to go further, in these words: fitting per-region scales is 'a separate registration and a"
    " separate decision'. Neither computed a held-out Delta MSE under a per-region predictor, and"
    " neither could: an in-sample dispersion statistic and an out-of-sample squared error are"
    " different quantities with different null behaviour. M2h's own chromosome-held-out check"
    " (lane-scalesplit, ec9cfd4) refitted ONE GLOBAL scale on the remaining chromosomes; it never"
    " fitted a scale PER chromosome, so it tested the global scale's portability and not a regional"
    " model. The question is therefore open and is asked here for the first time"
)
WHAT_IS_NOT_ASKED = (
    "NOT whether the per-chromosome scales differ -- scalehet and hetblock2 answered that and their"
    " readings are carried word for word, not re-decided and not strengthened. NOT what CAUSES the"
    " difference. NOT whether some other regional unit (locus block, replication timing domain,"
    " compartment, gene family) would do better: one unit is registered here and the result speaks"
    " about that unit at these gene counts. NOT whether a SHRUNKEN, pooled or hierarchical per-region"
    " estimator would help, which is a different estimator and a different registration. NOT any"
    " change to M2h's global scale, which is quoted and never refitted as a target. AND NOT what to"
    " DO about the answer: this lane PROPOSES NO FURTHER MODEL whatever it finds, and where a"
    " per-region scale helps or hurts is reported as a DESCRIPTIVE, never as a recommendation"
)

# ---- the population ------------------------------------------------------------------------------
POPULATION = (
    "7h's committed genes table for the multiplicative form, UNCHANGED: 16,509 genes, 214,617"
    " gene-cell pairs, thirteen P (the composed shift the compiled rules predict) and thirteen D (the"
    " measured deviation from the no-rule baseline) per gene, one chromosome label per gene as 7h"
    " matched it by gene_id on GRCh38. No gene is added, dropped or re-matched here. The split into"
    " training and held-out genes is M2h's own, reproduced bit-for-bit by the harness: 8,181 training"
    " and 8,328 held out"
)
REGIONAL_UNIT = (
    "THE CHROMOSOME, and the reason is that it is the unit the heterogeneity was MEASURED in. scalehet"
    " grouped by chromosome and rejected a single global scale there; hetblock2 kept the chromosome as"
    " the grouping and changed only the NULL's resampling unit to the 1-Mb locus block, and the"
    " rejection survived. A per-region model tested on any other unit would be answering a question no"
    " lane has yet asked, and a negative result on it could not be read against scalehet's positive"
    " one. So the unit is not chosen here: it is INHERITED, exactly as the two prior lanes defined it,"
    " and the 24 chromosome labels of 7h's committed table are the groups. chrY is in the table with"
    " 13 genes and does not clear the floor; chr21 has 156 and clears scalehet's floor but not this"
    " lane's, which is applied to each HALF of the split (see ADMISSION) -- and that is a consequence"
    " of the held-out design, stated here before any value was read"
)

# ---- THE HELD-OUT DESIGN, which is the whole lane ------------------------------------------------
HELD_OUT_DESIGN = (
    "STATED EXACTLY, AND IT IS THE WHOLE LANE. (1) THE PROBLEM. A per-region scale cannot be tested by"
    " leaving the region out. Leave-one-chromosome-out would require chromosome c's own scale to be"
    " fitted on genes that are not on chromosome c, which is not an estimate of chromosome c's scale"
    " at all -- it is the global scale under another name, and lane-scalesplit already ran exactly that"
    " and reported it (ec9cfd4). So leave-one-region-out is not a design for this estimator; it is a"
    " refusal of it, and saying 'we held the region out' while fitting its scale on the region would be"
    " the error this registration exists to prevent."
    " (2) THE RESOLUTION, REGISTERED: a WITHIN-REGION SPLIT on M2h's OWN split rule, carried unchanged"
    " -- magcal.SPLIT, the sha256 of the gene symbol, odd first 32 bits held out, which partitions"
    " every chromosome's genes into a training half and a held-out half with no reference to any"
    " measured value, to any chromosome and to any outcome. For each admitted chromosome c, the scale"
    " a_c is fitted by magcal's own closed form on THE TRAINING GENES OF CHROMOSOME c AND NOTHING"
    " ELSE, and is then applied to THE HELD-OUT GENES OF CHROMOSOME c, which contributed nothing to"
    " it. THE GENES FITTED ON AND THE GENES SCORED ON ARE DISJOINT, PER REGION, BY CONSTRUCTION."
    " (3) THE BASELINE PREDICTOR, scored on the same rows: the global scale a_g, fitted by the same"
    " closed form on ALL 8,181 training genes -- which IS M2h's +0.2145911381525326, and the harness"
    " proves it is by reproducing M2h's committed number to within 1e-9. It too saw no held-out gene."
    " (4) THE SCORED SET IS IDENTICAL FOR BOTH PREDICTORS: M2h's own 8,328 held-out genes, all of"
    " them, in the same order, with the same D. Nothing is dropped from one arm and kept in the other,"
    " because a comparison of two predictions must be a comparison on the same rows or it is a"
    " comparison of populations. A chromosome that does not clear the floor keeps its held-out genes"
    " in the scored set and is scored with a_g UNDER BOTH PREDICTORS (see ADMISSION), so the"
    " per-region arm differs from the global arm ONLY on the admitted chromosomes and the difference"
    " is attributable to the regional unit and to nothing else."
    " (5) WHAT IS NOT DONE: no scale is fitted on any gene it scores; no held-out value is consulted"
    " to choose a unit, a floor, a width, an alpha or an outcome, all of which are fixed in this"
    " committed module before the run exists; no second split, no re-split and no re-run on a"
    " different seed is permitted (STOP_RULE); and the per-chromosome scales are NOT computed by the"
    " harness, so no value of this design is seen before the registration is committed."
    " (6) THE DESIGN IS HONEST AT THESE GENE COUNTS and the thing that says so is a GENE COUNT, which"
    " depends on no measured value: 22 of the 24 chromosomes carry at least MIN_GENES = 100 training"
    " genes AND at least MIN_GENES = 100 held-out genes, which is above MIN_REGIONS = 10. HAD THEY NOT"
    " -- had fewer than 10 chromosomes carried both halves at the floor -- THE REGISTERED ANSWER WOULD"
    " HAVE BEEN 'a per-region scale cannot be tested out of sample at this gene count' AND THE LANE"
    " WOULD HAVE STOPPED THERE, reporting that as the finding; the outcome INCONCLUSIVE_NOT_TESTABLE"
    " below is that sentence, kept in the registration because the floor is checked by the run and not"
    " only by the harness"
)
ADMISSION = (
    "A chromosome is ADMITTED to the per-region arm when it carries at least MIN_GENES = 100 TRAINING"
    " genes AND at least MIN_GENES = 100 HELD-OUT genes. MIN_GENES is magcal's own floor, taken"
    " unchanged in value; what is new is that it is applied to EACH HALF rather than to the whole"
    " chromosome, and the reason is the held-out design: a scale fitted on fewer than magcal's floor of"
    " genes is below the precision magcal registered for a scale, and a Delta MSE read on fewer than"
    " that floor is below the precision magcal registered for a reading. Applying the floor to the"
    " whole chromosome, as scalehet did for an in-sample dispersion, would admit a chromosome whose"
    " training half is too small to fit on. A NON-ADMITTED CHROMOSOME IS NOT DROPPED: its held-out"
    " genes stay in the scored set and are scored with the GLOBAL scale a_g under BOTH predictors, so"
    " the two arms are scored on identical rows and the per-region arm is a strict modification of the"
    " global one on the admitted chromosomes alone. That fallback is a property of the PREDICTOR,"
    " fixed here, and not a post-hoc exclusion. The admitted set depends only on gene counts and is"
    " therefore INVARIANT under the registered null, which permutes labels while holding every group's"
    " training and held-out counts exactly -- so the null tests the LABELS and not the admission"
)
GENE_COUNTS_ARE_NOT_MEASUREMENTS = (
    "the counts that establish the design's feasibility -- genes per chromosome per half, the admitted"
    " set, the floor cleared or not -- depend on no measured value, which is the reason the harness is"
    " permitted to compute and commit them BEFORE the registration. lane-scalesplit's harness states"
    " the same fact in the same words. No scale, no Delta MSE, no difference of Delta MSEs, no"
    " bootstrap and no permutation is computed before the registration is committed"
)

# ---- the comparison ------------------------------------------------------------------------------
COMPARISON = (
    "ONE primary quantity: DELTA = Delta MSE(held-out genes, per-region predictor) - Delta MSE(held-out"
    " genes, global predictor), in the units of M2h's own metric -- squared error of measured"
    " log2(TPM+1) across the thirteen reader cells, against the no-rule baseline, averaged over"
    " gene-cell pairs. Both terms share the same denominator because both are read on the same"
    " 108,264 pairs, so DELTA is the mean over those pairs of the per-pair difference and the two arms'"
    " Delta MSEs are each reported beside it. POSITIVE DELTA MEANS THE PER-REGION MAGNITUDE IS CLOSER"
    " TO MEASURED RNA than the single global +0.2146; negative means it is further away; and the size"
    " of DELTA is reported against M2h's own +0.0433912108101358, which is about 2% of the variance,"
    " so a reader can see what fraction of the correction that already works a regional refinement"
    " would add or destroy. The fraction of variance is carried beside every 'better' and every"
    " 'worse', as the series requires"
)
STATISTIC = (
    "DELTA, computed through this module's `delta_mse_rowwise`, which takes ONE SCALE PER GENE and"
    " reduces -- when every gene's scale is the same number -- to magcal.delta_mse exactly. The harness"
    " checks that reduction within TOLERANCE on the real table, and tests/test_regionscale.py plants"
    " it, so the per-region arm and the global arm are computed by the same metric magcal registered"
    " and not by a lookalike. Per gene, the metric's numerator is 2 a_i (sum over its cells of P D) -"
    " a_i^2 (sum over its cells of P^2) and its denominator is its pair count; the fit is"
    " (sum over a group's genes of sum P D) / (sum over a group's genes of sum P^2), which is"
    " magcal.FIT's closed form reorganised as group sums and checked against magcal.fit_scale"
)
INTERVAL = (
    "98.33% percentile interval, magcal.FORM_INTERVAL taken UNCHANGED, DECIDED NOW AND NOT AFTER, on"
    " DELTA, by a PAIRED gene bootstrap: magcal's own resamples and magcal's own seed, the genes"
    " resampled with replacement as the unit, THE SAME RESAMPLED GENE INDICES USED FOR BOTH ARMS in"
    " every draw, and both arms' scales HELD FIXED at their fitted values, which is exactly what"
    " magcal.bootstrap_delta does with one scale. Pairing is required: the two arms differ on a"
    " minority of pairs and an unpaired comparison of two intervals would be dominated by the variance"
    " the arms share. WHY 98.33% AND NOT 95%: it is the width M2h's, scalesplit's, scalehet's and"
    " hetblock2's committed readings are stated at, so every number this lane prints beside theirs is"
    " at the same width and a comparison is not also a comparison of widths; and it is what magcal's"
    " own function computes, which this lane calls read-only rather than editing magcal's file. THE"
    " WIDTH IS STRICTER THAN NOMINAL AND THAT DIRECTION FAVOURS THIS LANE'S EXPECTED ANSWER, so the"
    " width is NOT ALLOWED TO DECIDE: the nominal 95% interval is computed from the SAME draws and"
    " reported beside it, and if the 98.33% interval includes 0 while the nominal 95% interval lies"
    " wholly above 0 the outcome is MIXED and NOT 'no improvement' (DECISION step 4). The arms' own"
    " 98.33% intervals are reported as DESCRIPTIVES and decide nothing"
)
ALPHA = (
    "ALPHA = 0.05, one-sided, magcal.CONTROL_ALPHA taken UNCHANGED. CORRECTION, DECIDED NOW: there are"
    " TWO registered readings in this lane -- the primary bootstrap interval of DELTA and the secondary"
    " permutation p -- and NO Bonferroni correction is applied over them, for the reason scalehet gave"
    " and this lane carries forward: THIS LANE'S EXPECTED ANSWER IS THE NULL ONE ('no improvement'),"
    " and splitting 0.05 between the two readings would halve the chance of detecting an improvement"
    " and so INFLATE the chance of the expected answer, making a cheap reading easy to get. The"
    " multiplicity is handled instead by requiring the two readings to AGREE for a positive verdict:"
    " IMPROVES needs the registered interval wholly above 0 AND the permutation p at or below 0.05, and"
    " disagreement reads MIXED, so an uncorrected pair of tests cannot buy a cleaner verdict than the"
    " evidence supports. The p against 0.0167 -- the 0.05/3 divisor M2h's 98.33% width came from, and"
    " the threshold at which hetblock2's secondary statistic did NOT clear -- is reported beside it as"
    " a DESCRIPTIVE so a reader can see both, but THE VERDICT USES 0.05 AND ONLY 0.05"
)
NULL = (
    "the 'no regional structure' control: THE CHROMOSOME LABELS ARE PERMUTED AMONG THE GENES, within"
    " the training half and within the held-out half separately, so every chromosome's training gene"
    " count and held-out gene count are held EXACTLY at their observed values and the admitted set is"
    " unchanged. Then the IDENTICAL pipeline runs on the permuted labels: a scale per permuted group"
    " fitted on that group's training genes, the global scale unchanged (it is fitted on all training"
    " genes and is therefore invariant under any permutation of labels, which the run records), the"
    " held-out genes scored, DELTA recomputed. Each gene carries its own P and D with it, so the null"
    " destroys the ASSOCIATION BETWEEN A GENE AND ITS REGION and nothing else. WHAT IT ANSWERS: how"
    " much held-out improvement does a per-region model produce from labels that carry no regional"
    " information? Under no regional structure a per-region predictor can only add estimation variance"
    " to the global one, so the null DELTA is expected to sit at or below 0, and the observed DELTA is"
    " read against that distribution. This is the lane's one control; M2h's within-gene cell"
    " permutation (magcal.permuted_control) is a control for a BETTER-THAN-NO-RULE reading and is run"
    " here ONCE AND ONLY in the harness, as part of reproducing M2h's committed numbers"
)
PERMUTATIONS = {"count": mc.PERMUTATIONS["count"], "seed": 20261013}
P_VALUE = (
    "one-sided, in the direction a real regional structure moves DELTA, and in magcal's own form:"
    " p = (1 + number of permutations with the null DELTA >= the observed DELTA) / (1 + count). It is"
    " one-sided because only real regional information can make a per-region predictor BEAT the global"
    " one out of sample -- random labels can only add estimation error -- so the low tail of the null"
    " is not a hypothesis this lane entertains. A p at the floor, 0.000999, means no permutation of"
    " 1000 reached the observed value and is reported as that and never as 0. The resolution of a"
    " 1000-draw permutation p is stated with it"
)
CONTROL = (
    "the registered null IS the control: a permutation of the chromosome labels is itself the 'no"
    " regional structure' control and no second control is registered. The harness additionally"
    " reproduces M2h's own within-gene cell permutation control, at M2h's count and seed, as part of"
    " proving M2h's committed numbers come back -- not as a control for this lane's statistic"
)

# ---- the harness ---------------------------------------------------------------------------------
HARNESS = (
    "run and committed BEFORE the registration, and it must pass or nothing is fitted, scored or read."
    " (1) 7h's committed genes table has 16,509 genes and 214,617 pairs, thirteen P and thirteen D per"
    " gene, no empty symbol, no empty chromosome, no repeated symbol. (2) its Delta MSE at scale 1"
    " reproduces 7h's committed -0.5361286982098615 within TOLERANCE through magcal.delta_mse."
    " (3) M2H'S COMMITTED NUMBERS COME BACK FIRST, AND ARE REPORTED BEFORE ANYTHING ELSE, through"
    " M2h's OWN functions at M2h's OWN resamples, seeds and permutation count and under M2h's OWN"
    " sha256-of-symbol split: training 8,181 genes, held out 8,328, fitted scale"
    " +0.2145911381525326 [+0.20074141299120243, +0.2283000359194418], held-out Delta MSE"
    " +0.0433912108101358 [+0.0379045975793137, +0.04894398453032398], control p = 0.000999000999,"
    " outcome BETTER -- every one within M2H_TOLERANCE = 1e-9, because the same code at the same seeds"
    " must agree closely. THIS IS THE CHECK THAT MAKES EVERY LATER NUMBER ATTRIBUTABLE TO THE REGIONAL"
    " UNIT AND TO NOTHING ELSE, and lane-scalesplit and lane-scalehet each passed it bit-identically."
    " (4) THE TWO PATHS ARE ONE PATH: this lane's `delta_mse_rowwise` with every gene given THE SAME"
    " scale equals magcal.delta_mse at that scale within TOLERANCE, and this lane's paired bootstrap"
    " with both arms given M2h's fitted scale returns a DELTA interval of EXACTLY [0, 0] while its"
    " per-arm interval equals magcal.bootstrap_delta's within TOLERANCE -- so the metric and the"
    " interval are magcal's and not a lookalike. WITHIN TOLERANCE AND NOT BIT-FOR-BIT, and the reason"
    " is recorded rather than smoothed over: magcal accumulates a gene's term cell by cell as"
    " 2 a P D - a^2 P^2, while this lane accumulates sum(P D) and sum(P^2) first -- because its"
    " permutation null must re-weight those two group sums at a thousand different scales -- and then"
    " combines them. The two are ALGEBRAICALLY IDENTICAL and differ only in floating-point summation"
    " ORDER, by of order 1e-17 on a bound of a bootstrap interval, far below TOLERANCE = 1e-4. The"
    " harness reports the observed difference as a number so the claim is checkable, and"
    " tests/test_regionscale.py pins both halves: the agreement within TOLERANCE and the fact that the"
    " residual is a summation-order artefact (it vanishes when the per-cell order is reproduced). The"
    " per-chromosome scales are NOT fitted here and the constant scale used for the reduction is M2h's"
    " own committed one, not a value of this design."
    " (5) the GENE COUNTS per chromosome per half, the admitted set under ADMISSION, and whether the"
    " admitted count clears MIN_REGIONS are reported; a count depends on no measured value"
    " (GENE_COUNTS_ARE_NOT_MEASUREMENTS), and these counts are what establish whether the held-out"
    " design is feasible at all. (6) the variance denominator of 'about X% of the variance' is"
    " reported -- the mean of D^2 over the held-out pairs, the no-rule baseline's own held-out MSE --"
    " and M2h's committed +0.04339 is confirmed to be about 2% of it; that denominator is a quantity"
    " M2h's committed reading already states, re-derived and not newly exposed."
    " (7) NO VALUE OF THIS LANE'S DESIGN IS COMPUTED BY THE HARNESS: no per-chromosome scale, no"
    " per-region Delta MSE, no DELTA, no bootstrap of DELTA and no permutation. A miss in any of (1)"
    " to (6) stops the lane with HARNESS_DEFECT and nothing is registered"
)

# ---- the exposure, stated the first of its two times ---------------------------------------------
PRIOR_EXPOSURE = (
    "THIS RNA IS NOT A FRESH MEASUREMENT AND THIS IS THE WEAKEST EXPOSURE POSITION OF ANY LANE IN THE"
    " SERIES. The same ENCODE polyA plus RNA-seq gene quantifications (V29) in the same thirteen reader"
    " biosamples have now been read by: item 7 (the information term), item 7d (lane-combine2, the"
    " runtime mean and promoter competition), item 7h (lane-multrna, the multiplicative composition"
    " whose committed table this lane reads), the capped-sum census, M2h (lane-magcal, the global"
    " magnitude), M2i (the fresh measurement with nothing fitted), lane-scalesplit (the"
    " chromosome-held-out refit), lane-recovered, lane-scalehet (the per-chromosome dispersion) and"
    " lane-hetblock2 (the same dispersion against a locus-block null). AND, SPECIFICALLY FOR THIS"
    " LANE: THE PER-CHROMOSOME SPREAD HAS NOW BEEN LOOKED AT, TESTED, AND TESTED AGAIN. The pattern"
    " this lane models is a pattern two committed lanes have already measured and published a p value"
    " for, and the regional unit is not chosen here but inherited from them. What that costs is"
    " stated plainly: a positive reading here cannot be called a discovery of regional variation,"
    " because the variation was discovered before this lane existed and this lane was opened BECAUSE"
    " of it; what a positive reading would add is only that the variation is USABLE OUT OF SAMPLE,"
    " which is a different and narrower claim. A negative reading is correspondingly cheap to believe"
    " and is the expected one. Nothing here is a fresh test of whether the compiled rules carry"
    " information about RNA at all -- that was M2h's and M2i's question -- and no reading of this lane"
    " may be quoted as one",
    "THE ONE THING THE EXPOSURE DOES NOT TOUCH: the held-out genes of THIS lane's design are M2h's own"
    " held-out genes, and the per-region scales that score them are fitted on the training half alone."
    " Exposure is a reason to discount a POSITIVE reading as a discovery and it is not a reason to"
    " doubt the ARITHMETIC of an out-of-sample comparison, whose two arms are scored on the same rows"
    " under the same metric. The design is stated so a reader can check that separation; it is not"
    " offered as a defence of the exposure",
)

# ---- the decision procedure, written before any value of this lane was computed ------------------
WHICH_GOVERNS = (
    "THE PAIRED BOOTSTRAP INTERVAL OF DELTA DECIDES, at the registered 98.33% width, with the"
    " permutation p required to AGREE for a positive verdict and the nominal 95% interval able to turn"
    " a negative verdict into MIXED. No descriptive, no per-chromosome number, no arm's own interval"
    " and no context number from M2h, 7h, scalesplit, scalehet or hetblock2 may change the verdict."
    " The outcome is the registered function of (the registered interval, the nominal interval, the"
    " permutation p, the admitted region count, the null's degeneracy) and of nothing else, and it is"
    " computed by `decide_region` below, which tests/test_regionscale.py plants in every one of its"
    " directions"
)
NOT_VACUOUS = (
    "the test can come out either way and both outcomes are statements. If the regional variation"
    " scalehet and hetblock2 found is large enough to survive being estimated from a few hundred"
    " training genes per chromosome, the per-region arm will beat the global arm on genes neither saw,"
    " DELTA will be positive, its 98.33% interval will clear 0 and the permutation p will be small --"
    " the floor 0.000999 is attainable at 1000 draws. If the variation is real but swamped by its own"
    " estimation error, DELTA will sit at or below 0 and read NO_IMPROVEMENT or WORSE. A planted test"
    " holds each direction: a table in which one chromosome's genes really do carry a different slope,"
    " at a size and a gene count that can be estimated, must read IMPROVES; a table whose slope is"
    " identical everywhere must read NO_IMPROVEMENT or WORSE and must not read IMPROVES. THE OUTCOME"
    " IS NOT PREDETERMINED BY THE DESIGN, and 'no improvement' is a FULL AND VALUABLE ANSWER, not an"
    " absence of one"
)

#: both outcomes pre-written, with their numbers named, before any held-out value was read
OUTCOMES = {
    "IMPROVES": (
        "A PER-REGION MAGNITUDE PREDICTS HELD-OUT RNA BETTER THAN THE SINGLE GLOBAL SCALE. On M2h's own"
        " {held_genes:,} held-out genes ({held_pairs:,} gene-cell pairs), neither arm having seen any of"
        " them, the per-chromosome magnitude -- {n_admitted} chromosomes each scaled by a number fitted"
        " on its own training genes alone, the remaining {n_fallback} scored by the global scale under"
        " both arms -- lowered the squared error of measured log2(TPM+1) across the thirteen reader"
        " cells below the single global {a_global:+.4f}'s: DELTA = {delta:+.6g}"
        " [{fw} {delta_lo:+.6g}, {delta_hi:+.6g}], permutation p = {p_perm:.6g} against ALPHA 0.05 at"
        " {n_perm:,} draws. In context, the per-region arm reads Delta MSE {dm_region:+.6g} against the"
        " global arm's {dm_global:+.6g}, so the whole correction rises from about {pct_global:.2f}% of"
        " the held-out variance to about {pct_region:.2f}% and THE REGIONAL PART OF IT IS ABOUT"
        " {pct_delta:.2f}% OF THE VARIANCE -- which is the size of the finding and must be carried"
        " beside the word 'better' wherever this reading is quoted. The nominal {nw} interval agrees at"
        " [{nom_lo:+.6g}, {nom_hi:+.6g}]. WHAT THIS DOES AND DOES NOT SAY: the regional variation"
        " lane-scalehet rejected a single global scale over (Q = 166.048, p_Q = 0.003996) and"
        " lane-hetblock2 confirmed against a 1-Mb locus-block null (p_Q = 0.008991, p_R = 0.025974, and"
        " at the stricter 0.0167 the primary Q clears and the secondary R does NOT) is not only REAL but"
        " USABLE OUT OF SAMPLE at the chromosome unit and at these gene counts. It says nothing about"
        " the cause, nothing about any other regional unit, and nothing about what should be built:"
        " WHAT TO DO ABOUT IT IS A SEPARATE DECISION AND THIS LANE PROPOSES NO FURTHER MODEL"
    ),
    "NO_IMPROVEMENT": (
        "NO IMPROVEMENT DETECTED: THE REGIONAL VARIATION IS REAL BUT NOT USABLE, AND THE SINGLE GLOBAL"
        " SCALE REMAINS THE WHOLE CORRECTION IN PRACTICE. On M2h's own {held_genes:,} held-out genes"
        " ({held_pairs:,} gene-cell pairs), neither arm having seen any of them, the per-chromosome"
        " magnitude -- {n_admitted} chromosomes each scaled by a number fitted on its own training genes"
        " alone, the remaining {n_fallback} scored by the global scale under both arms -- did NOT lower"
        " the squared error of measured log2(TPM+1) across the thirteen reader cells below the single"
        " global {a_global:+.4f}'s: DELTA = {delta:+.6g} [{fw} {delta_lo:+.6g}, {delta_hi:+.6g}], an"
        " interval that includes or lies below 0, with permutation p = {p_perm:.6g} against ALPHA 0.05"
        " at {n_perm:,} draws and the nominal {nw} interval agreeing at [{nom_lo:+.6g}, {nom_hi:+.6g}]."
        " In context, the per-region arm reads Delta MSE {dm_region:+.6g} against the global arm's"
        " {dm_global:+.6g}: about {pct_region:.2f}% and about {pct_global:.2f}% of the held-out variance"
        " respectively, so M2h's one number +0.2146 keeps the whole of the roughly 2% of the variance"
        " that the compiled rules explain and a regional refinement adds none of it. THIS IS A COMPLETE"
        " ANSWER AND NOT AN ABSENCE OF ONE, and it is the honest reconciliation of two true statements:"
        " lane-scalehet's rejection of a single global scale (Q = 166.048, p_Q = 0.003996) and"
        " lane-hetblock2's confirmation of it against a 1-Mb locus-block null (p_Q = 0.008991,"
        " p_R = 0.025974, and at the stricter 0.0167 the primary Q clears and the secondary R does NOT)"
        " STAND UNTOUCHED -- the per-chromosome scales really are not all one number -- and the"
        " variation is nonetheless too small, relative to the error of estimating it from a few hundred"
        " training genes per chromosome, to buy any prediction. A quantity can vary REALLY and vary"
        " UNUSABLY, and that is what this series has now measured. THE PRACTICAL CONSEQUENCE, stated"
        " once: the single global magnitude correction stays as it is, and NO PER-REGION MODEL IS"
        " PROPOSED BY THIS LANE"
    ),
    "WORSE": (
        "A PER-REGION MAGNITUDE PREDICTS HELD-OUT RNA WORSE THAN THE SINGLE GLOBAL SCALE, and the"
        " single global scale remains the whole correction in practice. On M2h's own {held_genes:,}"
        " held-out genes ({held_pairs:,} gene-cell pairs), neither arm having seen any of them, the"
        " per-chromosome magnitude -- {n_admitted} chromosomes each scaled by a number fitted on its own"
        " training genes alone, the remaining {n_fallback} scored by the global scale under both arms --"
        " RAISED the squared error of measured log2(TPM+1) across the thirteen reader cells above the"
        " single global {a_global:+.4f}'s, by an interval lying wholly below 0: DELTA = {delta:+.6g}"
        " [{fw} {delta_lo:+.6g}, {delta_hi:+.6g}], nominal {nw} [{nom_lo:+.6g}, {nom_hi:+.6g}],"
        " permutation p = {p_perm:.6g} at {n_perm:,} draws. The per-region arm reads Delta MSE"
        " {dm_region:+.6g} against the global arm's {dm_global:+.6g}: about {pct_region:.2f}% of the"
        " held-out variance against about {pct_global:.2f}%, so a regional refinement DESTROYS about"
        " {pct_delta_abs:.2f}% of the roughly 2% of the variance M2h's one number explains. This is the"
        " same answer as NO_IMPROVEMENT with its mechanism visible: the estimation error of a"
        " per-chromosome scale at these gene counts is larger than the regional signal it estimates."
        " lane-scalehet's and lane-hetblock2's readings STAND UNTOUCHED (Q = 166.048, p_Q = 0.003996;"
        " locus-block p_Q = 0.008991, p_R = 0.025974, and at the stricter 0.0167 the primary Q clears"
        " and the secondary R does NOT): the variation is real and it is not usable. NO PER-REGION MODEL"
        " IS PROPOSED BY THIS LANE"
    ),
    "MIXED": (
        "INCONCLUSIVE, THE TWO REGISTERED READINGS DISAGREE, and no verdict is taken. On M2h's own"
        " {held_genes:,} held-out genes ({held_pairs:,} gene-cell pairs) the per-region arm reads DELTA"
        " = {delta:+.6g} [{fw} {delta_lo:+.6g}, {delta_hi:+.6g}], nominal {nw} [{nom_lo:+.6g},"
        " {nom_hi:+.6g}], permutation p = {p_perm:.6g} against ALPHA 0.05 at {n_perm:,} draws, with the"
        " per-region arm at Delta MSE {dm_region:+.6g} and the global arm at {dm_global:+.6g} -- about"
        " {pct_region:.2f}% and about {pct_global:.2f}% of the held-out variance. {why} The registration"
        " requires the interval and the permutation to AGREE before 'better' may be said, and it"
        " requires the nominal width to agree before 'no improvement' may be said, PRECISELY so that"
        " neither a multiplicity nor a chosen interval width can buy a verdict; the condition is not"
        " met, so none is taken. lane-scalehet's and lane-hetblock2's readings stand untouched and are"
        " not weakened by this: the regional variation is real, and whether it is USABLE is what this"
        " lane leaves open. NO PER-REGION MODEL IS PROPOSED"
    ),
    "INCONCLUSIVE_NOT_TESTABLE": (
        "A PER-REGION SCALE CANNOT BE TESTED OUT OF SAMPLE AT THIS GENE COUNT, and that is the finding."
        " Only {n_admitted} chromosomes carry at least MIN_GENES = 100 training genes AND at least"
        " MIN_GENES = 100 held-out genes, below MIN_REGIONS = {min_regions}, so there is no honest"
        " held-out design for a per-region magnitude on this table: fitting a region's scale on the"
        " genes it is then scored on is not a test, and leaving the whole region out does not estimate"
        " that region's scale at all. NO DELTA IS REPORTED, no per-region scale is fitted and no"
        " comparison is made, because the registration forbids reading a held-out value under a design"
        " it has declared dishonest. lane-scalehet's and lane-hetblock2's readings stand untouched"
        " (Q = 166.048, p_Q = 0.003996; locus-block p_Q = 0.008991, p_R = 0.025974, and at the stricter"
        " 0.0167 the primary Q clears and the secondary R does NOT): the per-chromosome variation is"
        " real, and whether it is USABLE is UNANSWERED and not answerable this way. NO PER-REGION MODEL"
        " IS PROPOSED"
    ),
    "INCONCLUSIVE_DEGENERATE": (
        "INCONCLUSIVE: the fit is degenerate and no comparison is made. {what} A scale is undefined"
        " where the sum of P^2 over a group's training genes is 0, and the registration refuses to"
        " substitute a value for an undefined one. NO DELTA IS REPORTED. lane-scalehet's and"
        " lane-hetblock2's readings stand untouched. NO PER-REGION MODEL IS PROPOSED"
    ),
    "HARNESS_DEFECT": (
        "THE HARNESS DID NOT REPRODUCE AND NOTHING WAS REGISTERED, FITTED OR READ: {what}. No per-region"
        " scale, no DELTA and no comparison exists, and M2h's, 7h's, scalehet's and hetblock2's"
        " committed readings are unaffected"
    ),
}

MIXED_CLAUSES = {
    "interval_without_permutation": (
        "The interval clears 0 but the permutation p does not reach ALPHA, so the improvement is not"
        " distinguishable from what permuted chromosome labels produce."
    ),
    "permutation_without_interval": (
        "The permutation p reaches ALPHA but the registered interval does not clear 0, so the"
        " improvement is distinguishable from permuted labels and is still not distinguishable from 0."
    ),
    "width_decided_it": (
        "The registered 98.33% interval includes 0 while the nominal 95% interval lies wholly above 0,"
        " so the CHOSEN WIDTH and not the evidence would be deciding, and the registration refuses to"
        " let it: 'no improvement' may not be said here."
    ),
}

DECISION = (
    "Applied in this order and in no other, by `decide_region`. (1) If the harness did not reproduce,"
    " HARNESS_DEFECT and the lane stops; nothing below runs. (2) If fewer than MIN_REGIONS = 10"
    " chromosomes clear ADMISSION's two floors, INCONCLUSIVE_NOT_TESTABLE: the finding is that a"
    " per-region scale cannot be tested out of sample at this gene count, no DELTA is computed and the"
    " lane stops. (3) If any admitted group's training fit is undefined, INCONCLUSIVE_DEGENERATE and no"
    " DELTA is reported. (4) Otherwise, with L and H the registered 98.33% bounds of DELTA, l the lower"
    " nominal 95% bound from the same draws, and p the permutation p: IMPROVES if L > 0 AND p <= 0.05;"
    " MIXED if exactly one of (L > 0) and (p <= 0.05) holds; MIXED if neither holds but l > 0, because"
    " then the chosen width and not the evidence would be deciding; WORSE if neither holds and H < 0;"
    " NO_IMPROVEMENT otherwise. (5) The registered sentence for the outcome is emitted with its numbers"
    " filled and nothing else, and the prior-exposure clause is appended ALWAYS, whatever the outcome"
    " (`full_reading`). (6) The fraction of the held-out variance is carried beside every 'better' and"
    " every 'worse' in the sentence itself, so it cannot be dropped in a summary. (7) WHATEVER THE"
    " OUTCOME: no per-region model, no shrunken estimator, no second unit and no next step is proposed"
    " by this lane, and the per-chromosome table is a DESCRIPTIVE"
)
FALSIFIER = (
    "THE CLAIM UNDER TEST is that a per-region (per-chromosome) magnitude correction predicts held-out"
    " measured RNA better than M2h's single global +0.2146. IT IS FALSIFIED IF, on M2h's own held-out"
    " genes with every region's scale fitted on that region's training genes alone, DELTA's registered"
    " 98.33% paired-bootstrap interval includes or lies below 0 AND the nominal 95% interval from the"
    " same draws includes or lies below 0 AND the permutation p exceeds 0.05 -- that is, the outcomes"
    " NO_IMPROVEMENT and WORSE. It is NOT falsified, and nothing is concluded, when the two registered"
    " readings disagree (MIXED) or when the design's gene-count floors are not met"
    " (INCONCLUSIVE_NOT_TESTABLE) or a fit is undefined (INCONCLUSIVE_DEGENERATE)."
    " WHAT 'INCONCLUSIVE' IS, EXACTLY, so it cannot be stretched afterwards: it is one of those three"
    " named outcomes and nothing else. A small DELTA is NOT inconclusive -- a DELTA whose interval"
    " includes 0 is NO_IMPROVEMENT, which is a result. A DELTA that helps on some chromosomes and hurts"
    " on others is NOT inconclusive -- the primary quantity is the pooled DELTA on all held-out genes"
    " and the per-chromosome table is a descriptive carrying no correction over 22 groups and is"
    " therefore NOT 22 verdicts. A result that disagrees with scalehet's or hetblock2's expectation is"
    " NOT inconclusive. AND: NO_IMPROVEMENT DOES NOT FALSIFY scalehet OR hetblock2. Their statistic is"
    " an in-sample dispersion and this lane's is an out-of-sample squared error; a real variation that"
    " is too small to estimate makes both true at once, and that reconciliation is written into the"
    " NO_IMPROVEMENT sentence rather than discovered after the fact"
)

# ---- the exposure, stated the second of its two times, and carried by every reading -------------
PRIOR_EXPOSURE_CLAUSE = (
    "PRIOR EXPOSURE, AND IT IS THE WEAKEST OF ANY LANE IN THIS SERIES: this RNA is not a fresh"
    " measurement. The same ENCODE polyA plus RNA-seq quantifications in the same thirteen reader cells"
    " were read by item 7, item 7d, item 7h, the capped-sum census, M2h, M2i, lane-scalesplit,"
    " lane-recovered, lane-scalehet and lane-hetblock2 -- AND THE PER-CHROMOSOME SPREAD HAS NOW BEEN"
    " LOOKED AT, TESTED, AND TESTED AGAIN, by scalehet and then by hetblock2 against a locus-block"
    " null. The regional unit tested here was not chosen by this lane but inherited from those two, and"
    " this lane was opened because of their result. So a positive reading here is not a discovery of"
    " regional variation and may not be quoted as one; the only thing it could add is that the"
    " variation is USABLE OUT OF SAMPLE. A negative reading is the expected one and is cheap to"
    " believe. The arithmetic of the comparison is unaffected by the exposure -- both arms are scored"
    " on the same held-out rows under the same metric, and neither saw them -- and that separation is"
    " offered as something a reader can check, not as a defence"
)
STOP_RULE = (
    "ONE RUN. The registration is committed with its tests before any held-out value of this design is"
    " read, and the run refuses to start unless it exists, binds every BOUND field word for word, and"
    " finds every pinned input's sha256 unmoved; it refuses to overwrite its own result. No second"
    " seed, no second split, no second interval width, no second regional unit, no shrunken or pooled"
    " estimator and no 'with chr21 included' variant is run, now or later, under this registration."
    " WHATEVER THE OUTCOME, NO FURTHER MODEL IS PROPOSED: if a per-region scale helps, HOW MUCH AND"
    " WHERE is the finding and what to do about it is a separate decision for a separate lane; if it"
    " does not, the single global scale stands and the lane closes. No file of magcal's, scalehet's,"
    " hetblock's, scalesplit's or multrna's is edited by this lane, and magcal's functions are called"
    " read-only"
)
LIMITATIONS = (
    "The exposure above, which is the first and largest limitation and is stated twice.",
    "ONE regional unit, the chromosome, inherited from scalehet and hetblock2 and not chosen here. A"
    " chromosome is a very coarse region, and a negative result says that a per-CHROMOSOME magnitude is"
    " not usable out of sample at these gene counts; it does NOT say that no regional structure of any"
    " scale could be, and the registration forbids this lane from testing a second unit.",
    "ONE estimator, the unshrunken per-group closed form. It is the estimator whose dispersion scalehet"
    " measured, which is why it is the one tested; a shrunken, pooled or hierarchical estimator would"
    " trade bias for variance and might well behave differently, and that is a different registration.",
    "The gene counts are what they are: a few hundred training genes on a typical admitted chromosome"
    " and about a hundred on the smallest admitted one. The design is honest at those counts and it is"
    " not powerful at them, and a NO_IMPROVEMENT reading is therefore a statement about a"
    " per-chromosome scale ESTIMATED FROM THIS MUCH DATA and not about the regional variation's"
    " existence.",
    "chr21 (156 genes) and chrY (13) do not clear the two half-floors, so their held-out genes are"
    " scored by the global scale under both arms and the per-region arm is a modification of the global"
    " arm on the great majority, not the whole, of the scored rows. The exact counts are reported.",
    "One split, M2h's. The sha256-of-symbol partition is fixed and arbitrary with respect to everything"
    " measured, and it is also a single draw: a different partition would give a different DELTA, and"
    " the registration permits no second one, so the bootstrap interval describes gene resampling"
    " within this split and not variation across splits.",
    "Delta MSE across thirteen reader cells is one endpoint among many. It is M2h's and 7h's endpoint,"
    " carried so the comparison is to their committed numbers, and it is not a statement about any"
    " other measure of a rule's quality.",
    "The permutation null destroys gene-to-region association globally. hetblock2 showed that the"
    " per-chromosome dispersion survives a 1-Mb locus-block null, so proximity does not account for it;"
    " this lane's null is the simpler label permutation and it is a control for the PREDICTION"
    " question, not a re-test of the dispersion question hetblock2 settled.",
    "Of 440,558 compiled elements, hetblock2 measured that only 22 drive two or more genes of any kind"
    " and NOT ONE drives two of the 16,509 genes of this table, so gene-level dependence through a"
    " shared element does not occur here and the gene is a sound bootstrap unit. What remains is locus"
    " proximity, which hetblock2 tested and which this lane does not re-test.",
)
DESCRIPTIVES = (
    "Per admitted chromosome: its training gene count, its held-out gene count, the scale fitted on its"
    " training genes, the distance of that scale from the global one, its held-out Delta MSE under its"
    " own scale and under the global scale, and the difference. EACH DECIDING NOTHING, carrying NO"
    " correction over the admitted chromosomes and therefore NOT A TEST of any chromosome and NOT a row"
    " of verdicts. Named a descriptive here, before they are computed, so they cannot later be read as"
    " a finding about a named chromosome.",
    "The null distribution of DELTA: its median, its 98.33% percentile range and its quartiles, so a"
    " reader can see what held-out improvement permuted labels buy at these gene counts.",
    "Each arm's own held-out Delta MSE with its own 98.33% interval, and each as a percentage of the"
    " held-out variance, beside M2h's committed +0.0433912108101358.",
    "The fitted global scale of this run, which must equal M2h's committed +0.2145911381525326 within"
    " M2H_TOLERANCE because it is the same fit on the same training genes, and is reported as a check"
    " rather than as a value.",
)

#: the registration fields the run refuses to differ from: registration field -> constant of this module
BOUND = {
    "question": "QUESTION",
    "why_not_already_answered": "WHY_NOT_ALREADY_ANSWERED",
    "what_is_not_asked": "WHAT_IS_NOT_ASKED",
    "form": "FORM",
    "form_words": "FORM_WORDS",
    "source": "SOURCE",
    "carried_verbatim": "CARRIED",
    "population": "POPULATION",
    "regional_unit": "REGIONAL_UNIT",
    "held_out_design": "HELD_OUT_DESIGN",
    "admission": "ADMISSION",
    "gene_counts_are_not_measurements": "GENE_COUNTS_ARE_NOT_MEASUREMENTS",
    "comparison": "COMPARISON",
    "statistic": "STATISTIC",
    "interval": "INTERVAL",
    "form_interval": "FORM_INTERVAL",
    "form_interval_words": "FORM_INTERVAL_WORDS",
    "nominal_interval": "NOMINAL_INTERVAL",
    "nominal_interval_words": "NOMINAL_INTERVAL_WORDS",
    "bootstrap": "BOOTSTRAP",
    "alpha": "ALPHA",
    "control_alpha": "CONTROL_ALPHA",
    "null": "NULL",
    "permutations": "PERMUTATIONS",
    "p_value": "P_VALUE",
    "control": "CONTROL",
    "harness": "HARNESS",
    "min_genes": "MIN_GENES",
    "min_regions": "MIN_REGIONS",
    "tolerance": "TOLERANCE",
    "m2h_tolerance": "M2H_TOLERANCE",
    "m2h_committed": "M2H_COMMITTED",
    "seven_h_committed": "SEVEN_H_COMMITTED",
    "global_scale_committed": "GLOBAL_SCALE_COMMITTED",
    "prior_exposure": "PRIOR_EXPOSURE",
    "which_governs": "WHICH_GOVERNS",
    "not_vacuous": "NOT_VACUOUS",
    "outcomes": "OUTCOMES",
    "mixed_clauses": "MIXED_CLAUSES",
    "decision": "DECISION",
    "falsifier": "FALSIFIER",
    "prior_exposure_clause": "PRIOR_EXPOSURE_CLAUSE",
    "stop_rule": "STOP_RULE",
    "limitations": "LIMITATIONS",
    "descriptives": "DESCRIPTIVES",
}


def registered_fields() -> dict[str, Any]:
    """The BOUND fields, as the registration writes them and the run re-checks them."""
    out: dict[str, Any] = {}
    for field, name in BOUND.items():
        v = globals()[name]
        out[field] = list(v) if isinstance(v, tuple) else v
    return out


# ---- the within-region split and the admission ---------------------------------------------------
def split_counts(rows: list[dict]) -> dict[str, dict[str, int]]:
    """Per chromosome, the training and held-out gene COUNTS under magcal's own split rule.

    A COUNT, which depends on no measured value (GENE_COUNTS_ARE_NOT_MEASUREMENTS), which is why the
    harness is permitted to compute it before the registration is committed.
    """
    out: dict[str, dict[str, int]] = {}
    for r in rows:
        d = out.setdefault(r["chrom"], {"train": 0, "held": 0})
        d["held" if mc.is_held_out(r["symbol"]) else "train"] += 1
    return out


def admitted_regions(counts: dict[str, dict[str, int]]) -> tuple[list[str], list[str]]:
    """ADMISSION: (admitted, fallback), sorted. Admitted needs MIN_GENES in BOTH halves."""
    adm = sorted(c for c, d in counts.items() if d["train"] >= MIN_GENES and d["held"] >= MIN_GENES)
    fall = sorted(c for c in counts if c not in set(adm))
    return adm, fall


def is_testable(admitted: list[str]) -> bool:
    """HELD_OUT_DESIGN (6): is there an honest held-out design at these gene counts at all?"""
    return len(admitted) >= MIN_REGIONS


# ---- the fit and the metric, one scale per gene ---------------------------------------------------
def gene_sums(rows: list[dict]) -> tuple[list[float], list[float], list[float]]:
    """Per gene, (sum over its cells of P*D, sum of P^2, its pair count): every later sum is of these.

    scalehet.gene_sums is the same first two columns; the pair count is added because this lane's
    metric averages over pairs across groups and must not assume thirteen.
    """
    dps: list[float] = []
    pps: list[float] = []
    ns: list[float] = []
    for r in rows:
        dp = pp = 0.0
        for p, d in zip(r["P"], r["D"], strict=True):
            dp += p * d
            pp += p * p
        dps.append(dp)
        pps.append(pp)
        ns.append(float(len(r["D"])))
    return dps, pps, ns


def scale_from_sums(dp: float, pp: float) -> float | None:
    """magcal.FIT's closed form as group sums: scalehet.scale_from_sums, imported read-only."""
    return sh.scale_from_sums(dp, pp)


def gene_terms_rowwise(dps: list[float], pps: list[float], scales: list[float]) -> list[float]:
    """Per gene, 2 a_i (sum P D) - a_i^2 (sum P^2) with ONE SCALE PER GENE: magcal.METRIC's numerator."""
    return [2.0 * a * dp - a * a * pp for dp, pp, a in zip(dps, pps, scales, strict=True)]


def delta_mse_rowwise(dps: list[float], pps: list[float], ns: list[float], scales: list[float]) -> float:
    """METRIC with one scale per gene. With every scale equal it IS magcal.delta_mse at that scale,
    which HARNESS (4) checks on the real table and tests/test_regionscale.py plants."""
    terms = gene_terms_rowwise(dps, pps, scales)
    total = sum(ns)
    return sum(terms) / total if total else float("nan")


def region_scales(
    dps: list[float], pps: list[float], labels: list[str], admitted: list[str]
) -> tuple[dict[str, float], list[str]]:
    """HELD_OUT_DESIGN (2): one scale per admitted region from THAT region's TRAINING rows alone.

    Returns (scale per admitted region, the admitted regions whose fit is undefined). Call it with the
    TRAINING rows' sums; a region not admitted gets no scale and its genes fall back to the global
    scale under both arms (ADMISSION), which `scales_for` applies.
    """
    acc: dict[str, list[float]] = {c: [0.0, 0.0] for c in admitted}
    for dp, pp, c in zip(dps, pps, labels, strict=True):
        if c in acc:
            acc[c][0] += dp
            acc[c][1] += pp
    out: dict[str, float] = {}
    bad: list[str] = []
    for c in admitted:
        a = scale_from_sums(acc[c][0], acc[c][1])
        if a is None:
            bad.append(c)
        else:
            out[c] = a
    return out, bad


def scales_for(labels: list[str], per_region: dict[str, float], global_scale: float) -> list[float]:
    """The per-region arm's scale for each scored gene: its region's, or the global one (ADMISSION)."""
    return [per_region.get(c, global_scale) for c in labels]


# ---- the paired bootstrap -------------------------------------------------------------------------
def bootstrap_paired(
    dps: list[float],
    pps: list[float],
    ns: list[float],
    scales_region: list[float],
    scales_global: list[float],
    resamples: int,
    seed: int,
) -> dict[str, Any]:
    """INTERVAL: the paired gene bootstrap of DELTA, and of each arm, from the SAME draws.

    The draw sequence is magcal.bootstrap_delta's own -- `rng.integers(0, n, n)` per resample from
    `default_rng(seed)` -- so with a constant scale vector an arm's interval is magcal.bootstrap_delta's
    within TOLERANCE, which HARNESS (4) checks and states the reason for: the residual is a
    floating-point summation-ORDER difference of order 1e-17, not an algebraic one, because magcal
    accumulates a gene's term cell by cell while this lane accumulates sum(P D) and sum(P^2) first.
    Both arms' scales are held FIXED, as magcal holds its one scale fixed, and the SAME resampled
    indices serve both arms in every draw, which is what makes the DELTA interval of two identical arms
    EXACTLY [0, 0].
    """
    import numpy as np

    tr = np.array(gene_terms_rowwise(dps, pps, scales_region))
    tg = np.array(gene_terms_rowwise(dps, pps, scales_global))
    na = np.array(ns)
    n = len(ns)
    rng = np.random.default_rng(seed)
    d = np.empty(resamples)
    ar = np.empty(resamples)
    ag = np.empty(resamples)
    for i in range(resamples):
        idx = rng.integers(0, n, n)
        den = na[idx].sum()
        sr = tr[idx].sum()
        sg = tg[idx].sum()
        d[i] = (sr - sg) / den
        ar[i] = sr / den
        ag[i] = sg / den
    lo, hi = np.percentile(d, list(FORM_INTERVAL))
    nlo, nhi = np.percentile(d, list(NOMINAL_INTERVAL))
    rlo, rhi = np.percentile(ar, list(FORM_INTERVAL))
    glo, ghi = np.percentile(ag, list(FORM_INTERVAL))
    return {
        "delta_ci": [float(lo), float(hi)],
        "delta_ci_nominal": [float(nlo), float(nhi)],
        "region_arm_ci": [float(rlo), float(rhi)],
        "global_arm_ci": [float(glo), float(ghi)],
        "delta_draws_quantiles": quantiles([float(x) for x in d]),
        "resamples": resamples,
        "seed": seed,
    }


# ---- the null -------------------------------------------------------------------------------------
def permuted_delta(
    tr_dp: list[float],
    tr_pp: list[float],
    tr_lab: list[str],
    ho_dp: list[float],
    ho_pp: list[float],
    ho_n: list[float],
    ho_lab: list[str],
    admitted: list[str],
    global_scale: float,
    count: int,
    seed: int,
) -> list[float]:
    """NULL: the chromosome labels permuted within each half, counts held exactly, DELTA recomputed.

    Each gene keeps its own P and D, so only the gene-to-region association is destroyed. The global
    arm is invariant under the permutation -- its scale is fitted on all training genes and its scored
    value depends on no label -- so every permuted DELTA is read against the SAME global arm, which is
    what makes the null a null about the LABELS and about nothing else.
    """
    import numpy as np

    rng = np.random.default_rng(seed)
    idx_of = {c: i for i, c in enumerate(sorted({*tr_lab, *ho_lab}))}
    adm_set = {idx_of[c] for c in admitted}
    tr_dpa, tr_ppa = np.array(tr_dp), np.array(tr_pp)
    ho_dpa, ho_ppa, ho_na = np.array(ho_dp), np.array(ho_pp), np.array(ho_n)
    tr_l = np.array([idx_of[c] for c in tr_lab])
    ho_l = np.array([idx_of[c] for c in ho_lab])
    k = len(idx_of)
    total = float(ho_na.sum())
    g_terms = 2.0 * global_scale * ho_dpa - global_scale * global_scale * ho_ppa
    dm_global = float(g_terms.sum() / total)
    not_admitted = [j for j in range(k) if j not in adm_set]
    out: list[float] = []
    for _ in range(count):
        tl = tr_l[rng.permutation(len(tr_l))]
        hl = ho_l[rng.permutation(len(ho_l))]
        num = np.bincount(tl, weights=tr_dpa, minlength=k)
        den = np.bincount(tl, weights=tr_ppa, minlength=k)
        a = np.full(k, global_scale, dtype=float)
        ok = den != 0.0
        a[ok] = num[ok] / den[ok]
        for j in not_admitted:
            a[j] = global_scale
        av = a[hl]
        r_terms = 2.0 * av * ho_dpa - av * av * ho_ppa
        out.append(float(r_terms.sum() / total) - dm_global)
    return out


def permutation_p(observed: float, null: list[float]) -> float:
    """P_VALUE, magcal's own form, one-sided upward: scalehet.permutation_p, imported read-only."""
    return sh.permutation_p(observed, null)


def percentile_range(xs: list[float]) -> tuple[float, float]:
    """The registered 98.33% percentile range of a null distribution: scalehet's own, read-only."""
    return sh.percentile_range(xs)


# ---- the decision ---------------------------------------------------------------------------------
def decide_region(
    lo: float,
    hi: float,
    nominal_lo: float,
    p_perm: float,
    n_admitted: int,
    degenerate: list[str] | None = None,
) -> str:
    """DECISION steps (2), (3) and (4), applied exactly as written.

    The registered interval decides, the permutation must AGREE for a positive verdict, and the nominal
    width can only turn a negative verdict into MIXED -- never into a positive one.
    """
    if n_admitted < MIN_REGIONS:
        return "INCONCLUSIVE_NOT_TESTABLE"
    if degenerate:
        return "INCONCLUSIVE_DEGENERATE"
    interval_better = lo > 0.0
    perm_rejects = p_perm <= CONTROL_ALPHA
    if interval_better and perm_rejects:
        return "IMPROVES"
    if interval_better != perm_rejects:
        return "MIXED"
    if nominal_lo > 0.0:
        return "MIXED"
    return "WORSE" if hi < 0.0 else "NO_IMPROVEMENT"


def why_mixed(lo: float, p_perm: float) -> str:
    """The registered MIXED clause for the disagreement that occurred, and no other."""
    if lo > 0.0 and p_perm > CONTROL_ALPHA:
        return MIXED_CLAUSES["interval_without_permutation"]
    if lo <= 0.0 and p_perm <= CONTROL_ALPHA:
        return MIXED_CLAUSES["permutation_without_interval"]
    return MIXED_CLAUSES["width_decided_it"]


def reading(kind: str, **numbers: Any) -> str:
    """The registered sentence for an outcome, numbers filled, nothing else."""
    return OUTCOMES[kind].format(**numbers)


def full_reading(kind: str, **numbers: Any) -> str:
    """DECISION step (5): the registered sentence, then the prior-exposure clause, ALWAYS."""
    return reading(kind, **numbers) + ". " + PRIOR_EXPOSURE_CLAUSE


def quantiles(xs: list[float]) -> dict[str, float]:
    """The descriptive quartiles, magcal.quantiles' own convention, imported read-only."""
    return mc.quantiles(xs)
