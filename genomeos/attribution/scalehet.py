# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-scalehet: is M2h's magnitude ONE global scale, or does it vary by region?

M2h (`ac5580b`) fitted ONE scale for the whole genome, +0.2146 [98.33% CI +0.2007, +0.2283], and that
single number is what makes the compiled rules read better than no rule (held-out Delta MSE +0.04339,
about 2% of the variance). lane-scalesplit (`ec9cfd4`) showed the reading survives a
chromosome-held-out split and reported, POST HOC AND DECIDING NOTHING, within-chromosome scales from
+0.1232 (chr13) to +0.2874 (chr14), median +0.2176, q1 +0.2041, q3 +0.2376, no sign change, 10 of 24
inside M2h's interval, chrY flagged and excluded at 13 genes.

Nobody has tested whether that spread is more than sampling noise. A single global error suggests one
mis-specified constant in the compiled layer; regional variation suggests the error depends on
something about the locus, which is a different defect and a different repair. This lane asks the one
falsifiable question -- IS THE PER-CHROMOSOME SCALE VARIATION CONSISTENT WITH A SINGLE GLOBAL SCALE?
-- against a null that holds the per-chromosome GENE COUNTS at their observed values, so the null
carries the same sampling noise the observed spread was computed at.

Nothing of magcal's, scalesplit's, multrna's, combine2's, freshrna's, recovered's or scalefit's files
is edited. Every function that fits a scale, bootstraps, decides an outcome or reads a sentence about
the calibrated form is magcal's, imported by IDENTITY and not by copy. A scale is a MAGNITUDE
CORRECTION and never an explanation of expression.
"""

from __future__ import annotations

from typing import Any

from genomeos.attribution import magcal as mc
from genomeos.attribution import scalesplit as ss

REGISTERED = "2026-10-12"
LANE = "lane-scalehet"

# ---- what is carried, by identity ----------------------------------------------------------------
FORM = "multiplicative"
FORM_WORDS = mc.FORM_WORDS[FORM]
SOURCE = "data/results/multrna_result.json"

CELLS = mc.CELLS
MEASURED = mc.MEASURED
NO_RULE_BASELINE = mc.NO_RULE_BASELINE
#: magcal's own gene floor, taken UNCHANGED as this lane's admission threshold per chromosome
MIN_GENES = mc.MIN_GENES
#: magcal's own significance level, taken UNCHANGED (and deliberately: see ALPHA)
CONTROL_ALPHA = mc.CONTROL_ALPHA
FORM_INTERVAL = mc.FORM_INTERVAL
FORM_INTERVAL_WORDS = mc.FORM_INTERVAL_WORDS
TOLERANCE = mc.TOLERANCE
M2H_TOLERANCE = ss.M2H_TOLERANCE
M2H_COMMITTED = dict(ss.M2H_COMMITTED)
SEVEN_H_COMMITTED = dict(ss.SEVEN_H_COMMITTED)

CARRIED = {
    "endpoint": MEASURED,
    "rna_selection": mc.RNA_SELECTION,
    "pseudocount": mc.PSEUDOCOUNT,
    "no_rule_baseline": NO_RULE_BASELINE,
    "min_genes": MIN_GENES,
    "cells": list(CELLS),
    "selection": mc.SELECTION,
    "metric": mc.METRIC,
    "predictions": dict(mc.PREDICTIONS),
    "fit": mc.FIT,
    "form_interval": list(FORM_INTERVAL),
    "control_alpha": CONTROL_ALPHA,
    "note": (
        "the endpoint, the no-rule baseline, Delta MSE, the closed-form fit, the gene bootstrap, the"
        " within-gene cell permutation control, MIN_GENES, the 98.33% width, the four-way decision table"
        " and every reading word about the CALIBRATED FORM are 7d's, 7h's and M2h's, carried HERE WORD FOR"
        " WORD and imported from genomeos/attribution/magcal.py, which this lane does not edit; the"
        " chromosome-fold design and M2h's committed numbers are imported from"
        " genomeos/attribution/scalesplit.py, which this lane also does not edit. The population is 7h's,"
        " gene for gene, minus the chromosomes this lane's registered gene floor excludes. The ONLY thing"
        " this lane adds is a HETEROGENEITY STATISTIC on the per-chromosome scales and a null for it. No"
        " new endpoint, no new measurement, no new reader, no new form"
    ),
}

# ---- the registered question ---------------------------------------------------------------------
QUESTION = (
    "M2h fitted ONE scale for the whole genome, +0.2146 [98.33% CI +0.2007, +0.2283], and that single"
    " number is what makes the compiled rules read better than no rule (held-out Delta MSE +0.04339,"
    " about 2% of the variance). lane-scalesplit then reported, POST HOC AND DECIDING NOTHING,"
    " within-chromosome scales ranging +0.1232 (chr13) to +0.2874 (chr14), median +0.2176, q1 +0.2041,"
    " q3 +0.2376, with no sign change and 10 of 24 inside M2h's interval. THE FALSIFIABLE QUESTION ASKED"
    " HERE, AND NO OTHER: IS THAT PER-CHROMOSOME SCALE VARIATION CONSISTENT WITH A SINGLE GLOBAL SCALE?"
    " Equivalently: at the gene counts the 24 chromosomes actually carry (chr1 1,700 down to chrY 13), is"
    " the observed dispersion of the per-chromosome fitted scales larger than the dispersion a SINGLE"
    " global scale would itself produce? The question matters for the MODEL and not for the statistics:"
    " one global error is one mis-specified constant in the compiled layer and the repair is one"
    " constant; regional variation means the error depends on something about the locus, which is a"
    " different defect and a different repair. 'CONSISTENT WITH A SINGLE GLOBAL SCALE' IS A FULL ANSWER"
    " AND IS A REACHABLE OUTCOME OF THIS TEST, not a failure to find something"
)
WHY_NOT_ALREADY_ANSWERED = (
    "the spread scalesplit reported is POST HOC and was registered as deciding nothing: no interval was"
    " placed on any per-chromosome scale and no comparison between chromosomes was registered. A spread"
    " that wide is exactly what sampling noise can look like at these gene counts -- chr13 carries 280"
    " genes and chr18 238, and a scale fitted on a few hundred genes is a noisy number -- so the eye"
    " cannot tell a regional magnitude from a small-sample wobble, and nothing in scalesplit's design"
    " could. This lane does not re-report the spread; it tests it against a null that reproduces those"
    " very gene counts"
)
WHAT_IS_NOT_ASKED = (
    "this lane does NOT ask which region has which scale, does NOT fit a per-region scale, does NOT"
    " propose a regional model and does NOT read a per-chromosome scale as a verdict about that"
    " chromosome. If the test rejects a single global scale, THAT VARIATION EXISTS IS THE FINDING and the"
    " lane stops there: fitting per-region scales is a separate registration and a separate decision,"
    " which this lane is forbidden to pre-empt. Whatever the outcome, no held-out claim about a new"
    " endpoint is made: the chromosome-held-out reading is scalesplit's and is quoted, not re-decided"
)

# ---- the population -----------------------------------------------------------------------------
POPULATION = (
    "7h's committed genes table (data/results/multrna_result.json), 16,509 genes x 13 reader cells ="
    " 214,617 pairs, gene for gene the population M2h calibrated and scalesplit re-split -- RESTRICTED"
    " to the chromosomes admitted by MIN_GENES_RULE below. This lane admits no gene 7h did not admit and"
    " drops none it did except by that stated floor. P is the multiplicative composition's prediction and"
    " D the measured difference from the no-rule baseline, both read from 7h's committed rows and neither"
    " recomputed"
)
MIN_GENES_RULE = (
    "A CHROMOSOME IS ADMITTED ONLY IF IT CARRIES AT LEAST MIN_GENES = 100 GENES in 7h's committed table."
    " THE MINIMUM IS FIXED NOW, BEFORE ANY VALUE OF THIS LANE IS READ, and it is magcal's own MIN_GENES"
    " imported unchanged -- the floor M2h registered for a group to be read at all -- so it is not a"
    " number tuned here. WHY A FLOOR AT ALL: the statistic weights each chromosome by the precision of"
    " its own fitted scale, and the null holds the per-chromosome gene counts fixed, so a chromosome with"
    " a handful of genes contributes a wildly noisy scale to both sides; including it would make the test"
    " a statement about that one group's noise rather than about regional structure. WHY 100 AND NOT A"
    " NUMBER CHOSEN TO SUIT THIS TABLE: 100 is magcal's registered floor, it was scalesplit's floor"
    " before this lane existed, and scalesplit's committed harness already reported which chromosomes"
    " fall below it -- chrY alone, at 13 genes. A chromosome's GENE COUNT depends on no measured value"
    " (scalesplit's harness states this in those words), so naming the excluded set here reads no result."
    " EXPECTED EFFECT, stated now: 23 of the 24 chromosomes are admitted and chrY (13 genes) is excluded,"
    " which removes 13 of 16,509 genes. chrY's own fitted scale is reported among the descriptives,"
    " FLAGGED as below the floor and entering no statistic and no verdict. If the admitted count is below"
    " MIN_CHROMOSOMES = 10 the lane reads INCONCLUSIVE and no heterogeneity verdict is given"
)
MIN_CHROMOSOMES = 10

# ---- prior exposure -----------------------------------------------------------------------------
PRIOR_EXPOSURE = (
    "STATED PLAINLY, AND NOT A CAVEAT ADDED AFTERWARDS. THIS IS A RE-ANALYSIS OF A MEASUREMENT ALREADY"
    " READ, NOT A FRESH TEST. The same ENCODE polyA plus RNA-seq of the same thirteen reader biosamples,"
    " under the same log2(TPM+1) endpoint and the same no-rule baseline, has already been read by: item 7"
    " (M2b, lane-readercompile) including its post-hoc information term; item 7d (M2c, lane-combine2) for"
    " the runtime's mean and promoter competition; item 7h (M2f, lane-multrna) for the multiplicative"
    " composition; the capped-sum census (lane-capsumrna), which reached NOT_TESTABLE; M2h (lane-magcal),"
    " which fitted the one scale and read 'better' on a random split; M2i (lane-freshrna), which"
    " replicated it with nothing fitted on RNA fresh to this project; lane-scalesplit, which re-split the"
    " same rows by chromosome; and lane-recovered. The P and D this lane groups by chromosome ARE 7h's"
    " committed rows, which are the rows M2h calibrated and scalesplit re-split, SO EVERY VALUE THIS LANE"
    " READS HAS BEEN READ BEFORE.",
    "AND, SEPARATELY AND MORE SHARPLY: THE PER-CHROMOSOME SPREAD ITSELF WAS ALREADY SEEN. scalesplit"
    " computed and committed all 24 per-chromosome scales as post hoc descriptives, and the range"
    " (+0.1232 chr13 to +0.2874 chr14), the median (+0.2176), the quartiles (+0.2041, +0.2376), the"
    " absence of a sign change and the count inside M2h's interval (10 of 24) are in"
    " data/results/scalesplit_result.json and in this lane's brief. THIS LANE IS THEREFORE TESTING A"
    " PATTERN THAT HAS BEEN LOOKED AT, WHICH IS WEAKER THAN TESTING ONE THAT HAS NOT. Specifically: the"
    " statistic is registered before it is computed, but it is registered by someone who already knows"
    " roughly how wide the observed spread is and which two chromosomes sit at its ends; a statistic"
    " chosen after seeing a pattern can be chosen to flatter it, in either direction. Two things are done"
    " about that and neither removes it. (1) The PRIMARY statistic is a weighted dispersion over ALL"
    " admitted chromosomes, not a statistic about chr13 or chr14, and no chromosome is named in it. (2)"
    " The statistic the eye actually used -- the RANGE, max minus min, which is the +0.1232-to-+0.2874"
    " sentence -- is registered too, as a SECONDARY reading, precisely so that the pattern that was"
    " looked at is the pattern that is tested, and so that a disagreement between the two cannot be"
    " quietly dropped. Any reading from this lane holds as a statement about THIS measurement's"
    " dependence on the chromosome grouping. It is not independent confirmation of M2h, it is not a blind"
    " prediction of an unseen measurement, and a 'consistent' reading is consistency at these gene counts"
    " and this grouping, not a proof that one constant is exactly right",
)

# ---- what is computed ---------------------------------------------------------------------------
GROUPING = (
    "BY CHROMOSOME, as 7h's committed table names it, which is a property of the gene's ADDRESS and of no"
    " measured value. The chromosome is the region unit for one reason stated in advance: it is the unit"
    " scalesplit's post hoc spread was computed on, and this lane's whole purpose is to test that spread"
    " rather than a spread of its own choosing. A gene is wholly in one group with all thirteen of its"
    " pairs. THE CHROMOSOME IS A COARSE AND ARBITRARY REGION UNIT and that is a LIMITATION, written as"
    " one: regional variation at a finer scale -- by compiled element, by locus, by expression level, by"
    " rule count -- would not have to show up between chromosomes, so a 'consistent' reading here is"
    " consistency BETWEEN CHROMOSOMES and is not a statement about every way the error could vary"
)
NO_SPLIT = (
    "NO TRAINING AND HELD-OUT SPLIT IS USED HERE, AND THAT IS DELIBERATE AND IS NOT A WEAKENING. The"
    " question is about the DISPERSION OF FITTED SCALES, not about a held-out prediction, so the quantity"
    " of interest is in-sample by construction: a per-chromosome scale IS an in-sample fit on that"
    " chromosome's genes. A fitted dispersion is biased upward -- 23 scales fitted on 23 groups will"
    " differ even when one scale generated them all -- and THE NULL IS WHAT PRICES THAT BIAS EXACTLY,"
    " because the null statistic is computed by fitting 23 scales on 23 groups of the same sizes in the"
    " same way. Comparing a fitted dispersion to a null of fitted dispersions is the correct comparison;"
    " comparing it to zero would not be. The held-out question about whether the calibrated rules beat no"
    " rule is M2h's and scalesplit's, is already answered, and is QUOTED here, never re-decided"
)
STATISTIC = (
    "PRIMARY, Q: the precision-weighted dispersion of the per-chromosome fitted scales about the global"
    " fitted scale. For an admitted chromosome c, let dp_c = sum over that chromosome's pairs of D*P and"
    " pp_c = sum over that chromosome's pairs of P^2; its fitted scale is a_c = dp_c / pp_c, which is"
    " magcal.FIT applied to that chromosome's rows and nothing else (magcal.fit_scale, imported"
    " read-only). Let a_g = (sum over admitted chromosomes of dp_c) / (sum of pp_c), the ONE global scale"
    " fitted on the admitted population, which is magcal.fit_scale on exactly those rows. Then Q = sum"
    " over admitted chromosomes of pp_c * (a_c - a_g)^2. WHY THAT WEIGHT: under the squared-error fit"
    " this calibration minimises, the variance of a_c is proportional to 1 / pp_c, so pp_c is the"
    " precision of a_c and pp_c * (a_c - a_g)^2 is each chromosome's squared deviation measured in units"
    " of its own noise -- Cochran's homogeneity form. It makes Q large when a PRECISELY fitted chromosome"
    " sits far from the global scale and small when a noisily fitted one does, which is the behaviour the"
    " question needs. Q IS NOT COMPARED TO A CHI-SQUARED DISTRIBUTION and no asymptotic assumption is"
    " made: its null is the permutation null below, so the unknown residual variance cancels and no"
    " normality, independence-of-pairs or large-sample claim is required.",
    "SECONDARY, R: the RANGE of the per-chromosome fitted scales, max a_c - min a_c over the admitted"
    " chromosomes. This is the statistic the post hoc sentence used ('+0.1232 to +0.2874'), registered"
    " here so that the pattern that was looked at is the pattern that is tested. It is REPORTED BESIDE Q"
    " AND IS NOT THE VERDICT, for the reason that a range is driven by its two extreme groups and is the"
    " less powerful and less stable of the two. Its only registered force is the MIXED outcome below: if"
    " Q and R fall on opposite sides of ALPHA the lane says so in the reading and claims neither answer",
)
NULL = (
    "THE NULL IS 'ONE GLOBAL SCALE, AND THE PER-CHROMOSOME DEVIATIONS ARE ONLY WHICH GENES HAPPENED TO"
    " LAND ON WHICH CHROMOSOME.' It is generated by PERMUTING THE CHROMOSOME LABELS AMONG THE ADMITTED"
    " GENES: one permutation draws a uniform random assignment of the admitted genes to the admitted"
    " chromosome labels that holds EACH CHROMOSOME'S GENE COUNT EXACTLY AT ITS OBSERVED VALUE (chr1 keeps"
    " 1,700 slots, chr13 keeps 280, chr18 keeps 238, and so on), recomputes the 23 per-chromosome scales"
    " and the global scale from the permuted grouping by the same formulas, and records Q and R. A GENE"
    " IS PERMUTED AS A BLOCK, carrying all thirteen of its pairs, its own P values, its own D values and"
    " its own pair sums; no pair is broken apart from its gene and no P is detached from its D. Nothing"
    " about any gene's data changes under the null -- only which group it is counted in.",
    "WHY THIS NULL AND NOT ANOTHER, decided now. (1) IT HAS THE SAME PER-CHROMOSOME GENE COUNTS AS THE"
    " OBSERVED GROUPING, BY CONSTRUCTION, so it carries the same sampling noise the observed spread was"
    " computed at -- which is the precise objection the eye could not answer. (2) IT RESAMPLES WITHIN THE"
    " GLOBAL FIT: the whole admitted population, with its one global scale, is the thing being"
    " redistributed, so the null is exactly 'these genes, this one scale, a different accident of"
    " addresses'. (3) IT IS EXACT UNDER GENE EXCHANGEABILITY and assumes no distribution: it needs no"
    " model of the residuals, no variance estimate and no normality, and it prices the upward bias of"
    " fitting 23 scales because it fits 23 scales the same way. (4) It uses the REAL genes, so the"
    " heavy-tailed, heteroscedastic, decidedly non-normal distribution of real P and D values is present"
    " in the null exactly as it is in the observation. WHAT IT ASSUMES, AND THIS IS THE LIMITATION, NOT A"
    " CAVEAT: under the null, genes are EXCHANGEABLE between chromosomes. If the genes of one chromosome"
    " differ from another's in some way that moves the fitted scale but is not 'the compiled rules are"
    " wrong by a different amount there' -- gene density, expression level, rule count per gene, gene"
    " family -- then a rejection says the scale depends on SOMETHING ABOUT THE GROUP and does not say"
    " what. A rejection is therefore a finding of heterogeneity between chromosomes and NOT an"
    " identification of its cause, and this lane will not name one.",
    "PERMUTATION COUNT AND SEED, FIXED NOW: PERMUTATIONS['count'] = 1000 draws at PERMUTATIONS['seed'] ="
    " 20261012, one numpy default_rng seeded once. 1000 is magcal's own permutation count, taken"
    " unchanged. The smallest attainable p is 1/(1+1000) = 0.000999, which is far below ALPHA, so the"
    " test CAN reject; and 1000 draws resolve a p near ALPHA = 0.05 to about +-0.007 (one standard"
    " error), which is stated as the resolution of the p and is not refined after the fact",
)
PERMUTATIONS = {"count": mc.PERMUTATIONS["count"], "seed": 20261012}
P_VALUE = (
    "one-sided, in the direction heterogeneity moves the statistic, and in magcal's own form: p = (1 +"
    " number of permutations with the null statistic >= the observed statistic) / (1 + count). It is"
    " one-sided because both Q and R can only be INFLATED by real regional variation -- Q is a sum of"
    " squared deviations and R a range, and neither has a meaningful low tail to test. A p at the floor,"
    " 0.000999, means no permutation of 1000 reached the observed value and is reported as that and never"
    " as 0. The small-p resolution of a 1000-draw permutation p is stated with it"
)
ALPHA = (
    "ALPHA = 0.05, one-sided, magcal.CONTROL_ALPHA taken UNCHANGED, AND THE CHOICE IS DELIBERATELY THE"
    " ONE THAT MAKES THIS LANE'S EXPECTED ANSWER THE HARDER TO OBTAIN. The expected answer is 'consistent"
    " with a single global scale', which is a FAILURE TO REJECT; a strict threshold would make failing to"
    " reject easy and the reading cheap. 0.05 is the LAXER of the two defensible thresholds available --"
    " the other being 0.0167, the 0.05/3 Bonferroni divisor M2h's 98.33% interval came from -- so"
    " adopting 0.05 gives heterogeneity the EASIER passage and makes a 'consistent' reading the harder"
    " result to reach. The p against 0.0167 is reported beside it as a descriptive so a reader can see"
    " both, but THE VERDICT USES 0.05 AND ONLY 0.05. CORRECTION, DECIDED NOW: there are TWO registered"
    " readings in this lane, Q (primary) and R (secondary), and NO Bonferroni correction is applied to"
    " ALPHA over them -- for the same reason: splitting 0.05 between them would halve the chance of"
    " detecting heterogeneity and so inflate the chance of the expected 'consistent' answer. Instead the"
    " multiplicity is handled by the MIXED outcome, which refuses to pick a side when the two statistics"
    " disagree, so an uncorrected pair of tests cannot buy a cleaner verdict than the evidence supports."
    " This lane registers NO interval that excludes 0 and makes NO claim of an effect's sign, so the"
    " 98.33% width is used only to DESCRIBE the null's spread (see INTERVAL) and never to decide"
)
INTERVAL = (
    "98.33% percentile interval, magcal.FORM_INTERVAL taken UNCHANGED, DECIDED NOW AND NOT AFTER, and"
    " used for DESCRIPTION, with the verdict resting on the permutation p alone. What it is placed on,"
    " fixed here: (1) the NULL DISTRIBUTION OF Q -- its 98.33% percentile range and its median -- so a"
    " reader can see the observed Q against the range a single global scale itself produces; (2) the NULL"
    " DISTRIBUTION OF R, the same way, which is the direct answer in the units of the sentence that"
    " prompted this lane: how wide a min-to-max spread of per-chromosome scales does ONE global scale"
    " generate at these gene counts? (3) per admitted chromosome, the null's 98.33% range of THAT"
    " CHROMOSOME'S OWN fitted scale at its own gene count, beside its observed scale -- 23 bands, each"
    " deciding nothing, carrying NO correction over 23 and therefore NOT A TEST of any chromosome, named"
    " a descriptive here so it cannot later be read as 23 verdicts. WHY 98.33%: it is the width M2h and"
    " scalesplit's committed readings are stated at, so every number this lane prints beside theirs is at"
    " the same width and a comparison is not also a comparison of widths; and it is what magcal's own"
    " functions compute, which this lane calls read-only rather than editing magcal's file"
)
CONTROL = (
    "the registered null IS the control of this lane: a permutation of the chromosome labels is itself the"
    " 'no regional structure' control, and no second control is registered. M2h's within-gene cell"
    " permutation control (magcal.CONTROL, magcal.permuted_control) is a control for a HELD-OUT DELTA MSE"
    " reading and is run here ONCE AND ONLY in the harness, as part of reproducing M2h's committed"
    " numbers; it is not a control for a dispersion statistic and is not applied to one"
)
HARNESS = (
    "run and committed BEFORE the registration, and it must pass or nothing is grouped, fitted or read."
    " (1) 7h's committed genes table has 16,509 genes and 214,617 pairs, thirteen P and thirteen D per"
    " gene, no empty symbol, no empty chromosome, no repeated symbol. (2) its Delta MSE at scale 1"
    " reproduces 7h's committed -0.5361286982098615 within TOLERANCE through magcal.delta_mse. (3) M2H'S"
    " COMMITTED GENE-LEVEL NUMBERS COME BACK, through M2h's OWN functions at M2h's OWN resamples, seeds"
    " and permutation count and under M2h's OWN sha256-of-symbol split: training 8,181 genes, held out"
    " 8,328, fitted scale +0.2145911381525326 [+0.20074141299120243, +0.2283000359194418], held-out Delta"
    " MSE +0.0433912108101358 [+0.0379045975793137, +0.04894398453032398], control p = 0.000999000999,"
    " outcome BETTER -- every one within M2H_TOLERANCE = 1e-9, because the same code at the same seeds"
    " must agree closely. THIS IS THE CHECK THAT MAKES EVERY LATER NUMBER ATTRIBUTABLE TO THE GROUPING"
    " AND TO NOTHING ELSE, and lane-scalesplit passed it bit-identically. M2h's split is reproduced in"
    " THE HARNESS ONLY and is NOT this lane's design. (4) this lane's own vectorised per-chromosome scale"
    " path agrees with magcal.fit_scale applied to each chromosome's rows, within TOLERANCE, for every"
    " chromosome in the table -- so the statistic and its null are computed by the same fit magcal"
    " registered and not by a lookalike. (5) the gene COUNT per chromosome is reported and the admitted"
    " set under MIN_GENES is listed; a count depends on no measured value. (6) THE STATISTIC AND ITS NULL"
    " ARE NOT COMPUTED BY THE HARNESS: no Q, no R, no permutation and no p is produced before the"
    " registration is committed. A miss in any of (1) to (5) stops the lane with HARNESS_DEFECT"
)

# ---- the decision procedure, written before any statistic of this lane was computed --------------
WHICH_GOVERNS = (
    "THE PERMUTATION p OF Q DECIDES, at ALPHA = 0.05, one-sided. R's p can only turn a verdict into"
    " MIXED; it can never, by itself, declare either answer. No interval, no descriptive, no"
    " per-chromosome band and no context number from M2h, 7h or scalesplit may change the verdict. The"
    " outcome is the registered function of (p_Q, p_R, the admitted chromosome count, the null's"
    " degeneracy) and of nothing else, and it is computed by `decide_heterogeneity` below, which is"
    " planted in all of its directions by tests/test_scalehet.py"
)
NOT_VACUOUS = (
    "the test can come out either way and both outcomes are statements. If the per-chromosome scales were"
    " driven by real regional variation, Q would exceed the dispersion that label permutation produces"
    " and p_Q would be small -- the floor 0.000999 is attainable at 1000 draws. If one global scale"
    " generated them, Q would sit inside the null's bulk and p_Q would be large. A planted test holds"
    " each direction: genes on one chromosome given a different slope must read HETEROGENEOUS, and genes"
    " whose slope is identical everywhere must read CONSISTENT_WITH_ONE_SCALE. THE OUTCOME IS NOT"
    " PREDETERMINED BY THE DESIGN, and 'consistent' is reported as a FULL ANSWER and not as an absence"
)
OUTCOMES = {
    "CONSISTENT_WITH_ONE_SCALE": (
        "CONSISTENT WITH A SINGLE GLOBAL SCALE. The dispersion of the {n_admitted} admitted chromosomes'"
        " fitted scales about the global scale {a_global:+.4f} is Q = {q_obs:.6g}, and under the"
        " registered null -- one global scale, the chromosome labels permuted among the same genes with"
        " every chromosome's gene count held at its observed value -- the null Q has median {q_med:.6g}"
        " and 98.33% range [{q_lo:.6g}, {q_hi:.6g}], giving p = {p_q:.6g} against ALPHA 0.05. The"
        " secondary range statistic agrees: the observed max-minus-min spread {r_obs:+.4f} (from"
        " {a_min:+.4f} on {chrom_min} to {a_max:+.4f} on {chrom_max}) lies against a null range with"
        " median {r_med:+.4f} and 98.33% range [{r_lo:+.4f}, {r_hi:+.4f}], p = {p_r:.6g}. SO THE"
        " PER-CHROMOSOME SPREAD scalesplit SAW POST HOC IS WHAT A SINGLE GLOBAL SCALE PRODUCES AT THESE"
        " GENE COUNTS, AND THIS TEST FINDS NO REGIONAL VARIATION IN THE MAGNITUDE BETWEEN CHROMOSOMES:"
        " read as it is registered, THE MAGNITUDE IS ONE CONSTANT the compiled rules are globally wrong"
        " by, and the repair indicated is one constant. The constant itself is M2h's {m2h_scale:+.4f}"
        " [98.33% CI +0.2007, +0.2283], which buys about 2% of the variance and is a MAGNITUDE CORRECTION"
        " AND NOT AN EXPLANATION OF EXPRESSION. Three limits bind this reading and are not caveats added"
        " afterwards: consistency is between CHROMOSOMES, a coarse region unit, and says nothing about"
        " variation by element, locus, expression level or rule count; a failure to reject is not a proof"
        " that one constant is exactly right, and the null's 98.33% Q range above is the honest measure of"
        " what this test could have detected; and the RNA is a measurement already read eight times over,"
        " with this very spread already seen post hoc"
    ),
    "HETEROGENEOUS": (
        "NOT CONSISTENT WITH A SINGLE GLOBAL SCALE. The dispersion of the {n_admitted} admitted"
        " chromosomes' fitted scales about the global scale {a_global:+.4f} is Q = {q_obs:.6g}, against a"
        " null -- one global scale, the chromosome labels permuted among the same genes with every"
        " chromosome's gene count held at its observed value -- whose Q has median {q_med:.6g} and 98.33%"
        " range [{q_lo:.6g}, {q_hi:.6g}]: p = {p_q:.6g} against ALPHA 0.05. The secondary range"
        " statistic agrees: the observed spread {r_obs:+.4f} (from {a_min:+.4f} on {chrom_min} to"
        " {a_max:+.4f} on {chrom_max}) against a null range of median {r_med:+.4f} and 98.33% range"
        " [{r_lo:+.4f}, {r_hi:+.4f}], p = {p_r:.6g}. SO THE MAGNITUDE THE COMPILED RULES ARE WRONG BY"
        " VARIES BETWEEN CHROMOSOMES BY MORE THAN THE SAMPLING NOISE OF THESE GENE COUNTS. THAT IS THE"
        " WHOLE FINDING AND THIS LANE STOPS HERE: no regional model is proposed, no per-region scale is"
        " fitted, no cause is named and no chromosome is read as a verdict about itself. A rejection of"
        " gene exchangeability between chromosomes says the fitted scale depends on SOMETHING about the"
        " group -- which could be gene density, expression level, rules per gene or gene family as easily"
        " as a locus-dependent defect -- and distinguishing those is a separate registration and a"
        " separate decision. M2h's one scale {m2h_scale:+.4f} and its about 2% of the variance are"
        " unaffected as a global average and are not withdrawn by this reading; what is withdrawn is the"
        " assumption that one constant is the whole of the correction. The RNA is a measurement already"
        " read eight times over and this very spread was already seen post hoc, which makes this a weaker"
        " test than one of an unlooked-at pattern"
    ),
    "MIXED": (
        "NEITHER ANSWER IS CLAIMED: the two registered statistics disagree across ALPHA = 0.05. The"
        " primary precision-weighted dispersion gives Q = {q_obs:.6g} with p = {p_q:.6g}, and the"
        " secondary range gives {r_obs:+.4f} with p = {p_r:.6g}, one of the two below 0.05 and the other"
        " above, on the same {n_admitted} admitted chromosomes, the same global scale {a_global:+.4f} and"
        " the same {count} permutations at seed {seed}. The registered rule for this case was fixed before"
        " any statistic was computed and is applied as written: when a dispersion that weights every"
        " chromosome by its own precision and a range driven by its two extreme chromosomes fall on"
        " opposite sides of the threshold, the evidence does not support either 'one global scale' or"
        " 'regional variation', and this lane says so rather than reporting whichever agrees with its"
        " expectation. Both p values, both nulls and all {n_admitted} per-chromosome scales are reported"
        " so a later lane can see exactly which statistic saw what. No regional model is proposed"
    ),
    "INCONCLUSIVE": (
        "INCONCLUSIVE, and nothing about the magnitude is claimed either way: {what}. No verdict on"
        " whether the magnitude is one global scale or regional is given, M2h's and scalesplit's committed"
        " readings stand untouched, and no regional model is proposed"
    ),
    "HARNESS_DEFECT": mc.OUTCOMES["HARNESS_DEFECT"],
}
DECISION = (
    "applied in this order, every step fixed before any statistic of this lane was computed. (1) THE"
    " HARNESS FIRST: if any of HARNESS (1) to (5) misses, the lane reads HARNESS_DEFECT and nothing is"
    " grouped, fitted or read. (2) ADMISSION: chromosomes with fewer than MIN_GENES = 100 genes are"
    " excluded (MIN_GENES_RULE); if fewer than MIN_CHROMOSOMES = 10 remain, the lane reads INCONCLUSIVE."
    " (3) DEGENERACY: if the global fit is undefined (sum of P^2 zero), if any admitted chromosome's fit"
    " is undefined, or if the null Q distribution is degenerate (its 98.33% range has zero width, i.e."
    " the permutation produced no spread at all), the lane reads INCONCLUSIVE with the reason named. (4)"
    " THE PRIMARY TEST: p_Q = (1 + permutations with null Q >= observed Q) / (1 + 1000), one-sided. (5)"
    " THE SECONDARY TEST: p_R, the same form on the range. (6) THE OUTCOME, by this table and no other:"
    " p_Q > 0.05 and p_R > 0.05 -> CONSISTENT_WITH_ONE_SCALE; p_Q <= 0.05 and p_R <= 0.05 ->"
    " HETEROGENEOUS; the two on opposite sides of 0.05 -> MIXED, which claims neither. (7) THE READING IS"
    " THE REGISTERED SENTENCE FOR THE OUTCOME, numbers filled and nothing else, followed by the prior"
    " exposure clause, always. (8) IF THE OUTCOME IS HETEROGENEOUS, NO REGIONAL MODEL IS PROPOSED IN THE"
    " READING, IN THE RESULT OR IN THE ROADMAP ROW: that what varies is a finding is the whole of it"
)
FALSIFIER = (
    "THE CLAIM UNDER TEST IS 'ONE GLOBAL SCALE': that the per-chromosome scales differ only by which"
    " genes landed on which chromosome. IT IS FALSIFIED by p_Q <= 0.05 -- an observed precision-weighted"
    " dispersion larger than 95% of the dispersions a single global scale produces at these very gene"
    " counts -- in which case this lane reports HETEROGENEOUS and 'the magnitude is one constant' may not"
    " be said. IT IS NOT FALSIFIED by p_Q > 0.05, in which case the lane reports"
    " CONSISTENT_WITH_ONE_SCALE, which is a FULL ANSWER to the question asked and not an absence of one."
    " WHAT WOULD HAVE FALSIFIED THE OPPOSITE CLAIM: 'the magnitude varies by region' predicts a Q above"
    " the null's 98.33% range and an observed spread wider than the null's; a Q sitting at the null's"
    " median is that prediction failing. WHAT 'INCONCLUSIVE' IS, exactly and in advance: fewer than"
    " MIN_CHROMOSOMES = 10 chromosomes admitted; an undefined global or per-chromosome fit; a degenerate"
    " null Q distribution; or the harness missing, which is reported as HARNESS_DEFECT. INCONCLUSIVE IS"
    " NOT a p near the threshold, which is reported with its value and decided by the table, and it is"
    " NOT the MIXED case, which is its own named outcome. WHAT NEITHER OUTCOME CAN SAY: neither can"
    " identify what makes a scale differ, neither can speak about regions finer than a chromosome,"
    " neither revises M2h's about 2% of the variance, and neither makes a scale an explanation of"
    " expression rather than a magnitude correction"
)
PRIOR_EXPOSURE_CLAUSE = (
    "carried with every reading, always: this RNA is not a fresh measurement. The same ENCODE polyA plus"
    " RNA-seq of the same thirteen reader biosamples under the same endpoint was read by items 7, 7d and"
    " 7h, by the capped-sum census, by M2h, by M2i, by lane-scalesplit and by lane-recovered, and the P"
    " and D grouped here are 7h's committed rows. AND the per-chromosome spread tested here was ALREADY"
    " SEEN POST HOC by scalesplit, so this is a test of a pattern that has been looked at, which is"
    " weaker than a test of one that has not. A reading here is a statement about THIS measurement's"
    " dependence on the chromosome grouping and a hypothesis for a fresh one, not independent"
    " confirmation of anything"
)
STOP_RULE = (
    "one run of scripts/scalehet_run.py, after this registration is committed. The statistic, the null,"
    " the permutation count, the seed, the gene floor, the admitted-chromosome minimum, the threshold and"
    " the decision table are fixed here and are not revisited after a number is seen. NOTHING IS RE-RUN"
    " AT ANOTHER SEED, another permutation count, another gene floor or another region unit, and the run"
    " refuses to overwrite its own result file. If the outcome is HETEROGENEOUS, no per-region scale is"
    " fitted by this lane under any circumstances"
)
LIMITATIONS = (
    "THE REGION UNIT IS THE CHROMOSOME, which is coarse and biologically arbitrary. It was chosen because"
    " it is the unit the post hoc spread was computed on, so that this lane tests that spread rather than"
    " one of its own choosing. Regional variation by compiled element, by locus, by expression level or by"
    " rules per gene need not show between chromosomes, so a 'consistent' reading is consistency BETWEEN"
    " CHROMOSOMES and nothing wider",
    "THE NULL ASSUMES GENES ARE EXCHANGEABLE BETWEEN CHROMOSOMES. A rejection therefore says the fitted"
    " scale depends on something about the group, not that the compiled rules' error is locus-dependent;"
    " gene density, expression level, rules per gene and gene family would all produce it. This lane names"
    " no cause",
    "GENES THAT SHARE A COMPILED ELEMENT ARE NOT INDEPENDENT. The permutation unit is the GENE, so"
    " correlation between genes within a chromosome is NOT priced into the null, and if it is substantial"
    " the null's Q is too narrow and p_Q is anti-conservative -- the error runs toward DETECTING"
    " heterogeneity, which is against this lane's expected answer and so does not flatter it. An"
    " element-level or locus-block permutation is not registered and is not computed: it would need a"
    " gene-to-element map this lane does not build, and naming it as the fix is a separate registration",
    "A FAILURE TO REJECT IS NOT A PROOF OF IDENTITY. The null's own 98.33% Q range and the null's range of"
    " per-chromosome spreads are reported as the honest measure of what this test could have detected: a"
    " regional variation smaller than that spread would not be seen here",
    "NO HELD-OUT PREDICTION IS MADE BY THIS LANE. The statistic is in-sample by construction and the null"
    " prices that exactly; the held-out reading that the calibrated rules beat no rule is M2h's and"
    " scalesplit's, is quoted verbatim, and is not re-decided here",
    "THE RNA IS A MEASUREMENT ALREADY READ, by items 7, 7d, 7h, the capped-sum census, M2h, M2i,"
    " scalesplit and lane-recovered, AND THE SPREAD TESTED HERE WAS ALREADY SEEN POST HOC. Both are"
    " stated in PRIOR_EXPOSURE and both weaken every reading below what the same test on an unlooked-at"
    " pattern in a fresh measurement would carry",
    "A SCALE IS A MAGNITUDE CORRECTION AND NEVER AN EXPLANATION OF EXPRESSION. Whatever this lane reads,"
    " the calibrated compiled rules buy about 2% of the variance on M2h's held-out side, and that number"
    " is not improved, revised or reinterpreted by anything here",
)
DESCRIPTIVES = (
    "deciding nothing, reported beside the registered reading: every admitted chromosome's fitted scale"
    " with its gene count, pair count, dp_c, pp_c and weighted squared deviation, in the order the"
    " deviation ranks; the excluded chromosomes' scales and counts, FLAGGED as below MIN_GENES and"
    " entering no statistic; the global fitted scale on the admitted population and on all 24"
    " chromosomes, beside M2h's committed +0.2146; the null distributions of Q and R, each with its"
    " median, its 98.33% percentile range and its minimum and maximum over the 1000 draws; per admitted"
    " chromosome, the null's 98.33% range of that chromosome's own scale at its own gene count beside its"
    " observed scale, with NO correction over 23 and therefore NOT a test of any chromosome; the"
    " descriptive p of Q against 0.0167 as well as against 0.05; how many admitted chromosomes lie inside"
    " M2h's committed scale interval and how many have a negative scale; and the median, quartiles,"
    " minimum and maximum of the admitted per-chromosome scales, recomputed here so that scalesplit's"
    " post hoc sentence is reproduced rather than quoted from memory",
    "context, copied verbatim and not re-decided: 7h's registered reading, M2h's registered"
    " multiplicative reading with its committed numbers, and lane-scalesplit's registered reading and its"
    " committed post hoc per-chromosome spread",
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
    "min_genes": "MIN_GENES",
    "min_genes_rule": "MIN_GENES_RULE",
    "min_chromosomes": "MIN_CHROMOSOMES",
    "prior_exposure": "PRIOR_EXPOSURE",
    "grouping": "GROUPING",
    "no_split": "NO_SPLIT",
    "statistic": "STATISTIC",
    "null": "NULL",
    "permutations": "PERMUTATIONS",
    "p_value": "P_VALUE",
    "alpha": "ALPHA",
    "control_alpha": "CONTROL_ALPHA",
    "interval": "INTERVAL",
    "form_interval": "FORM_INTERVAL",
    "form_interval_words": "FORM_INTERVAL_WORDS",
    "control": "CONTROL",
    "harness": "HARNESS",
    "tolerance": "TOLERANCE",
    "m2h_tolerance": "M2H_TOLERANCE",
    "m2h_committed": "M2H_COMMITTED",
    "seven_h_committed": "SEVEN_H_COMMITTED",
    "which_governs": "WHICH_GOVERNS",
    "not_vacuous": "NOT_VACUOUS",
    "outcomes": "OUTCOMES",
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


# ---- the grouping and the admission --------------------------------------------------------------
def gene_counts(rows: list[dict]) -> dict[str, int]:
    """Genes per chromosome. A COUNT, which depends on no measured value."""
    out: dict[str, int] = {}
    for r in rows:
        out[r["chrom"]] = out.get(r["chrom"], 0) + 1
    return out


def admitted_chromosomes(counts: dict[str, int]) -> tuple[list[str], list[str]]:
    """MIN_GENES_RULE: (admitted, excluded) chromosome names, each sorted, by the gene floor alone."""
    adm = sorted(c for c, n in counts.items() if n >= MIN_GENES)
    exc = sorted(c for c, n in counts.items() if n < MIN_GENES)
    return adm, exc


# ---- the per-gene sums the statistic and its null are built from ---------------------------------
def gene_sums(rows: list[dict]) -> tuple[list[float], list[float]]:
    """Per gene, (sum of D*P over its pairs, sum of P^2 over its pairs): the permutation's block.

    A gene's two sums are all the statistic needs from it, and a gene is permuted as one block, so no
    pair is ever separated from its gene and no P from its D. magcal.fit_scale on any set of rows is
    exactly (sum of those rows' dp) / (sum of those rows' pp), which HARNESS (4) checks per chromosome.
    """
    dps, pps = [], []
    for r in rows:
        dp = pp = 0.0
        for p, d in zip(r["P"], r["D"], strict=True):
            dp += p * d
            pp += p * p
        dps.append(dp)
        pps.append(pp)
    return dps, pps


def scale_from_sums(dp: float, pp: float) -> float | None:
    """magcal.FIT in the sums it is written in: dp / pp, or None when the denominator is 0."""
    return None if pp == 0.0 else dp / pp


# ---- the statistics ------------------------------------------------------------------------------
def heterogeneity_q(scales: list[float], weights: list[float], global_scale: float) -> float:
    """STATISTIC, primary: Q = sum over groups of pp_c * (a_c - a_g)^2, Cochran's homogeneity form."""
    return sum(w * (a - global_scale) ** 2 for a, w in zip(scales, weights, strict=True))


def range_stat(scales: list[float]) -> float:
    """STATISTIC, secondary: max a_c - min a_c, the spread the post hoc sentence used."""
    return max(scales) - min(scales) if scales else float("nan")


def permutation_p(observed: float, null: list[float]) -> float:
    """P_VALUE: (1 + null draws >= observed) / (1 + count), one-sided, magcal's own form."""
    return (1 + sum(1 for x in null if x >= observed)) / (1 + len(null))


def observed_statistics(
    dps: list[float], pps: list[float], labels: list[str], groups: list[str]
) -> dict[str, Any]:
    """The observed per-group scales, the global scale, Q and R, on the admitted groups alone.

    `groups` is the admitted chromosome list in its registered order; genes whose label is not in it
    are dropped by MIN_GENES_RULE and enter nothing.
    """
    idx = {c: [] for c in groups}
    for i, lab in enumerate(labels):
        if lab in idx:
            idx[lab].append(i)
    dpc = [sum(dps[i] for i in idx[c]) for c in groups]
    ppc = [sum(pps[i] for i in idx[c]) for c in groups]
    a_g = scale_from_sums(sum(dpc), sum(ppc))
    if a_g is None or any(p == 0.0 for p in ppc):
        return {"degenerate": "a fitted scale is undefined: a group's sum of P^2 is 0"}
    scales = [d / p for d, p in zip(dpc, ppc, strict=True)]
    return {
        "degenerate": None,
        "groups": list(groups),
        "gene_counts": [len(idx[c]) for c in groups],
        "dp": dpc,
        "pp": ppc,
        "scales": scales,
        "global_scale": a_g,
        "q": heterogeneity_q(scales, ppc, a_g),
        "r": range_stat(scales),
    }


def permutation_null(
    dps: list[float],
    pps: list[float],
    sizes: list[int],
    count: int,
    seed: int,
) -> dict[str, Any]:
    """NULL: the chromosome labels permuted among the genes, each group's gene COUNT held exactly.

    `sizes` is the admitted groups' observed gene counts in the registered order and must sum to the
    number of genes passed in. Each draw shuffles the genes, cuts the shuffle into blocks of exactly
    those sizes, refits a scale per block and the global scale, and records Q, R and every block's own
    scale. A GENE IS ONE BLOCK: its two sums travel together, so no pair leaves its gene.
    """
    import numpy as np

    if sum(sizes) != len(dps):
        raise ValueError(f"the group sizes sum to {sum(sizes)}, not the {len(dps)} genes given")
    dpa, ppa = np.array(dps, dtype=float), np.array(pps, dtype=float)
    bounds = list(np.cumsum(sizes)[:-1])
    rng = np.random.default_rng(seed)
    qs: list[float] = []
    rs: list[float] = []
    per_group: list[list[float]] = [[] for _ in sizes]
    globals_: list[float] = []
    for _ in range(count):
        parts = np.split(rng.permutation(len(dps)), bounds)
        dpc = np.array([dpa[p].sum() for p in parts])
        ppc = np.array([ppa[p].sum() for p in parts])
        if (ppc == 0.0).any():
            return {"degenerate": "a permuted group's sum of P^2 is 0, so its scale is undefined"}
        a_g = float(dpc.sum() / ppc.sum())
        a = dpc / ppc
        globals_.append(a_g)
        qs.append(float((ppc * (a - a_g) ** 2).sum()))
        rs.append(float(a.max() - a.min()))
        for i, v in enumerate(a):
            per_group[i].append(float(v))
    return {
        "degenerate": None,
        "q": qs,
        "r": rs,
        "per_group": per_group,
        "global_scales": globals_,
        "count": count,
        "seed": seed,
    }


def percentile_range(xs: list[float]) -> tuple[float, float]:
    """INTERVAL: the 98.33% percentile range of a null distribution, magcal.FORM_INTERVAL unchanged."""
    import numpy as np

    lo, hi = np.percentile(np.array(xs, dtype=float), list(FORM_INTERVAL))
    return float(lo), float(hi)


# ---- the decision --------------------------------------------------------------------------------
def decide_heterogeneity(p_q: float, p_r: float, n_admitted: int, degenerate: str | None = None) -> str:
    """DECISION steps (2), (3) and (6), applied exactly as written. p_Q decides; p_R can only MIX."""
    if degenerate:
        return "INCONCLUSIVE"
    if n_admitted < MIN_CHROMOSOMES:
        return "INCONCLUSIVE"
    q_rejects = p_q <= CONTROL_ALPHA
    r_rejects = p_r <= CONTROL_ALPHA
    if q_rejects != r_rejects:
        return "MIXED"
    return "HETEROGENEOUS" if q_rejects else "CONSISTENT_WITH_ONE_SCALE"


def reading(kind: str, **numbers: Any) -> str:
    """The registered sentence for an outcome, numbers filled, nothing else."""
    return OUTCOMES[kind].format(**numbers)


def full_reading(kind: str, **numbers: Any) -> str:
    """DECISION step (7): the registered sentence, then the prior exposure clause, always."""
    return reading(kind, **numbers) + ". " + PRIOR_EXPOSURE_CLAUSE


def quantiles(xs: list[float]) -> dict[str, float]:
    """The descriptive quartiles, magcal.quantiles' own convention, imported read-only."""
    return mc.quantiles(xs)
