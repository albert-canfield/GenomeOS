# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-hetblock: does scalehet's scale heterogeneity survive a LOCUS-BLOCK null?

lane-scalehet (`b2bddb4` code, `036eec6` registration, `40b0a8d` result) asked whether M2h's single
global magnitude correction (+0.2146) is one number or varies by region, and read NOT CONSISTENT WITH A
SINGLE GLOBAL SCALE: Q = 166.048 against a null of median 78.9561 and 98.33% range [35.4193, 143.929],
p_Q = 0.003996; the secondary range statistic agreeing, +0.1642 observed against a null median +0.0950,
p_R = 0.013986. It reproduced M2h bit-identically first, registered its statistic and its null before
any value, and set ALPHA at the LAXER 0.05 so that "consistent" would be the harder result.

AND IT FLAGGED, AS ITS OWN LIMITATION, THE ONE THING THAT COULD ALONE EXPLAIN ITS REJECTION, in these
words: "THE NULL ASSUMES GENES ARE EXCHANGEABLE BETWEEN CHROMOSOMES... GENES THAT SHARE A COMPILED
ELEMENT ARE NOT INDEPENDENT. The permutation unit is the GENE, so correlation between genes within a
chromosome is NOT priced into the null, and if it is substantial the null's Q is too narrow and p_Q is
anti-conservative". That correlation is established in this project -- lane-freshrna's own POST_HOC
says "genes that share elements or a locus are not independent and the gene-level interval may be too
narrow" and priced it with a 1-Mb LOCUS-CLUSTER bootstrap -- so scalehet's heterogeneity may be an
artefact of a null that is too narrow. THE NAMED, UNREGISTERED FIX IS AN ELEMENT- OR LOCUS-BLOCK
PERMUTATION, and that is this lane and nothing else.

This lane changes ONE thing about scalehet's null: THE EXCHANGEABLE UNIT. The statistic is scalehet's
Q, by identity; the secondary range statistic is scalehet's R, by identity; the p is magcal's one-sided
form, by identity; the gene floor, the 23 admitted chromosomes, the 98.33% description, the draw count
and the seed are scalehet's, unchanged; every chromosome's gene count is held EXACTLY as scalehet held
it. Nothing of scalehet's, magcal's, scalesplit's, recovered's, scalefit's, residue's, freshrna's or
multrna's files is edited; every function that fits a scale, forms a statistic, reads a p or reads a
sentence about the calibrated form is imported BY IDENTITY and not by copy.

NO REGIONAL MODEL IS PROPOSED EITHER WAY. That variation exists, or is not established, is the finding;
fitting per-region scales is a separate registration and a separate decision, which this lane is
forbidden to pre-empt. A scale is a MAGNITUDE CORRECTION and never an explanation of expression, and
M2h's +0.2146 with its about 2% of the variance stands as a GLOBAL AVERAGE whatever this lane reads.
"""

from __future__ import annotations

from typing import Any

from genomeos.attribution import magcal as mc
from genomeos.attribution import scalehet as sh

REGISTERED = "2026-10-12"
LANE = "lane-hetblock"

# ---- what is carried, by identity ----------------------------------------------------------------
FORM = sh.FORM
FORM_WORDS = sh.FORM_WORDS
SOURCE = sh.SOURCE
#: the gene table positions are joined from, repository-relative (d72f959)
GENE_TABLE = "data/cache/gencode_genes.tsv"
#: the directory the compiled non-coding programs are read from for the element-sharing graph
COMPILED_PROGRAMS = "data/knowledge/compiled/noncoding_{chrom}.bio"

CELLS = sh.CELLS
MEASURED = sh.MEASURED
NO_RULE_BASELINE = sh.NO_RULE_BASELINE
#: scalehet's gene floor, which is magcal's, taken UNCHANGED
MIN_GENES = sh.MIN_GENES
#: scalehet's significance level, which is magcal's, taken UNCHANGED at the LAXER 0.05
CONTROL_ALPHA = sh.CONTROL_ALPHA
#: the complement of the 98.33% description, reported beside ALPHA and deciding nothing
STRICTER_ALPHA = 0.0167
FORM_INTERVAL = sh.FORM_INTERVAL
FORM_INTERVAL_WORDS = sh.FORM_INTERVAL_WORDS
TOLERANCE = sh.TOLERANCE
MIN_CHROMOSOMES = sh.MIN_CHROMOSOMES
M2H_COMMITTED = dict(sh.M2H_COMMITTED)
SEVEN_H_COMMITTED = dict(sh.SEVEN_H_COMMITTED)
#: scalehet's committed result, the figures this lane reproduces before anything else
SCALEHET_COMMITTED = {
    "q_observed": 166.04792148316798,
    "r_observed": 0.16420247148502493,
    "global_scale": 0.21501821644770602,
    "p_q": 0.003996003996003996,
    "p_r": 0.013986013986013986,
    "null_q_median": 78.95613453032229,
    "null_q_98_33": [35.419346005111905, 143.9288653653726],
    "null_r_median": 0.09497160409012753,
    "null_r_98_33": [0.052091893649909106, 0.17056618399692858],
    "admitted_count": 23,
    "genes_admitted": 16496,
    "pairs_admitted": 214448,
    "outcome": "HETEROGENEOUS",
}

CARRIED = {
    "endpoint": MEASURED,
    "no_rule_baseline": NO_RULE_BASELINE,
    "min_genes": MIN_GENES,
    "min_chromosomes": MIN_CHROMOSOMES,
    "cells": list(CELLS),
    "metric": mc.METRIC,
    "fit": mc.FIT,
    "form_interval": list(FORM_INTERVAL),
    "control_alpha": CONTROL_ALPHA,
    "statistic_q": sh.STATISTIC,
    "note": (
        "the endpoint, the no-rule baseline, the closed-form fit, MIN_GENES, the 98.33% width, the"
        " one-sided p and every reading word about the CALIBRATED FORM are 7d's, 7h's and M2h's; the"
        " primary statistic Q, the secondary range statistic R, the gene floor, MIN_CHROMOSOMES, the"
        " admitted chromosomes, the draw count and the seed are lane-scalehet's. All of them are IMPORTED"
        " from genomeos/attribution/magcal.py and genomeos/attribution/scalehet.py, neither of which this"
        " lane edits. THE ONLY THING THIS LANE CHANGES IS THE EXCHANGEABLE UNIT OF THE NULL: a locus"
        " block instead of a gene. No new endpoint, no new measurement, no new reader, no new form, no"
        " new statistic and no new threshold"
    ),
}

# ---- the registered question ---------------------------------------------------------------------
QUESTION = (
    "lane-scalehet read NOT CONSISTENT WITH A SINGLE GLOBAL SCALE -- Q = 166.048 against a null of"
    " median 78.9561 and 98.33% range [35.4193, 143.929], p_Q = 0.003996, with the range statistic"
    " agreeing at p_R = 0.013986 -- from a null whose EXCHANGEABLE UNIT IS THE GENE, and it flagged as"
    " its own limitation that genes sharing a compiled element are not independent, so that its null's Q"
    " may be too narrow and its p_Q anti-conservative IN EXACTLY THE REJECTING DIRECTION. THE FALSIFIABLE"
    " QUESTION ASKED HERE, AND NO OTHER: DOES THAT HETEROGENEITY SURVIVE A NULL WHOSE EXCHANGEABLE UNIT"
    " IS A LOCUS BLOCK RATHER THAN A GENE? Equivalently: when the null moves whole blocks of neighbouring"
    " genes -- so that the correlation between genes of one neighbourhood is carried into the null rather"
    " than broken by it -- is the observed dispersion Q = 166.048 still larger than the null produces at"
    " the same per-chromosome gene counts? 'THE HETEROGENEITY IS NOT ESTABLISHED UNDER A BLOCK NULL' IS A"
    " FULL ANSWER AND A REACHABLE OUTCOME OF THIS TEST, AND IT IS THE MORE VALUABLE ONE: it would WITHDRAW"
    " scalehet's rejection to NOT ESTABLISHED, which is a correction of this project's own committed"
    " reading and not a failure to find something"
)
WHY_NOT_ALREADY_ANSWERED = (
    "scalehet named the fix and did not run it, in its own words: 'An element-level or locus-block"
    " permutation would price it and is not done here'. Nobody has run it since. The dependence is not"
    " hypothetical in this project: lane-freshrna's committed POST_HOC (1) says 'genes that share elements"
    " or a locus are not independent and the gene-level interval may be too narrow' and reported a 1-Mb"
    " LOCUS-CLUSTER interval beside its gene-level one for that reason, and lane-scalesplit split by"
    " chromosome in the first place because a gene-level resample does not price between-region"
    " structure. So the one quantity that could alone explain scalehet's rejection has never been priced,"
    " and the direction of the error is known in advance: a gene-level null is TOO NARROW, which pushes"
    " p_Q DOWN, which is the direction that rejects"
)
WHAT_IS_NOT_ASKED = (
    "this lane does NOT ask which region has which scale, does NOT fit a per-region scale, does NOT"
    " propose a regional model, does NOT name a cause and does NOT read a per-block or per-chromosome"
    " scale as a verdict about that block or that chromosome. WHATEVER THE OUTCOME, NO REGIONAL MODEL IS"
    " PROPOSED: that variation exists, or is not established, is the finding and the lane stops there."
    " This lane also does not re-decide M2h's held-out reading, does not revise +0.2146 and does not"
    " revise the about 2% of the variance: those stand as a GLOBAL AVERAGE either way, and what is at"
    " stake here is only whether one constant is the whole of the correction. And this lane does not"
    " re-test scalehet's question against a GENE null: that null is scalehet's, is reproduced here as a"
    " harness check, and is not re-decided"
)

# ---- the population -----------------------------------------------------------------------------
POPULATION = (
    "scalehet's population, gene for gene: 7h's committed genes table"
    " (data/results/multrna_result.json), 16,509 genes x 13 reader cells = 214,617 pairs, restricted by"
    " MIN_GENES to the 23 chromosomes scalehet admitted (16,496 genes, 214,448 pairs; chrY's 13 genes"
    " excluded). P is the multiplicative composition's prediction and D the measured difference from the"
    " no-rule baseline, both read from 7h's committed rows and neither recomputed. NO GENE IS ADDED AND"
    " NONE IS DROPPED that scalehet did not drop: if this lane's population differed from scalehet's by"
    " one gene the two nulls would not be comparable and the comparison is the whole point"
)
MIN_GENES_RULE = (
    "SCALEHET'S FLOOR, UNCHANGED AND NOT RE-CHOSEN: a chromosome is admitted only if it carries at least"
    " MIN_GENES = 100 genes in 7h's committed table, which is magcal's registered floor imported by"
    " identity. 23 of the 24 chromosomes are admitted and chrY (13 genes) is excluded, exactly as"
    " scalehet committed. If the admitted count is below MIN_CHROMOSOMES = 10 this lane reads"
    " INCONCLUSIVE. A chromosome's GENE COUNT depends on no measured value, so naming the admitted set"
    " here reads no result"
)

# ---- the exchangeable unit: the whole of this lane's change ---------------------------------------
BLOCK_UNIT = (
    "THE PERMUTATION UNIT IS A LOCUS BLOCK, NOT A GENE, and that is the whole of this lane's change to"
    " scalehet's null. THE PRIMARY UNIT IS THE 1-Mb LOCUS BLOCK: a gene belongs to the block"
    " (chromosome, gene start // 1,000,000), the gene's start read from data/cache/gencode_genes.tsv and"
    " joined on gene_id with the version suffix stripped, _PAR_Y rows skipped and the first occurrence"
    " kept. THAT IS lane-freshrna'S OWN COMMITTED UNIT, word for word from its POST_HOC (1) -- 'genes"
    " assigned to clusters by (chromosome, gene start // 1,000,000) from data/cache/gencode_genes.tsv'"
    " -- and this lane groups with freshrna.locus_clusters, THE SAME FUNCTION, imported by identity and"
    " not reimplemented. The unit is therefore NOT ONE CHOSEN HERE: it is the unit this project already"
    " committed for precisely this dependence, and choosing it costs this lane the freedom to pick a span"
    " that suits its answer. A gene travels in exactly one block with all thirteen of its pairs, its own"
    " P and its own D, exactly as a gene travelled as one block in scalehet's null"
)
BLOCK_UNIT_WHY_NOT_ELEMENTS = (
    "THE BRIEF PREFERRED THE ELEMENT-SHARING COMPONENTS -- connected components of the graph on genes with"
    " an edge where two genes share a compiled element -- as the stricter and more faithful unit, BECAUSE"
    " THAT IS THE DEPENDENCE SCALEHET NAMED. THEY WERE COMPUTED FIRST, and the graph HAS NO EDGES. Read"
    " from the 24 compiled non-coding programs' rule lines through combine2.rule_lines -- an address-level"
    " parse of the committed programs, running nothing and reading no measured value -- the programs carry"
    " 440,558 elements, of which only 22 drive two or more genes OF ANY KIND, AND NOT ONE DRIVES TWO OF"
    " THE 16,509 GENES OF THIS TABLE. Taking the components over ALL targeted genes, so that a chain"
    " through a gene outside the table would still join two of the table's genes, EVERY ONE OF THE 16,509"
    " GENES IS ITS OWN COMPONENT. So the element-sharing unit IS COMPUTABLE AND IS EXACTLY THE GENE UNIT:"
    " its null is scalehet's null, bit-identically at the same seed, and it can price nothing. THAT IS"
    " ITSELF A FINDING AND IS REPORTED AS ONE: the specific dependence scalehet feared -- two of these"
    " genes driven by one compiled element -- DOES NOT OCCUR IN THIS POPULATION, because the compiler"
    " attributes each element to one gene. WHAT REMAINS, and is not removed by that, is LOCUS PROXIMITY:"
    " neighbouring genes share chromatin, share the reader's open/held calls, share local sequence and are"
    " measured in the same neighbourhood of the same libraries, which is the dependence freshrna priced"
    " with a 1-Mb cluster and the one still unpriced in scalehet's null. HENCE THE 1-Mb LOCUS BLOCK IS THE"
    " PRIMARY UNIT. The element-sharing component null is run and reported anyway, as the degenerate case"
    " it is, so that the claim 'it prices nothing' is a committed number and not an assertion"
)
BLOCK_SELECTION_RULE = (
    "FIXED HERE, BEFORE ANY STATISTIC OF THIS LANE IS COMPUTED, AS A RULE AND NOT AS A PREFERENCE. The"
    " governing unit is the ELEMENT-SHARING COMPONENT if the compiled programs contain at least"
    " MIN_SHARED_ELEMENTS = 1 element driving two or more of the table's genes -- that is, if the"
    " element-sharing graph has ANY edge at all; otherwise it is the 1-Mb LOCUS BLOCK, because a component"
    " graph with no edges IS the gene unit and a null built on it cannot price any dependence. The count"
    " that decides it is a count of rule lines in committed programs and depends on NO MEASURED VALUE,"
    " exactly as scalehet took its per-chromosome gene counts and named its excluded chromosome before"
    " registering. STATED PLAINLY AND NOT HIDDEN: that count was taken before this module was written and"
    " is 0, so the rule selects the 1-Mb locus block; the rule is written as a rule so that the choice is"
    " a stated function of an address-level count rather than this lane's taste, and the harness"
    " recomputes and commits the count before the registration"
)
MIN_SHARED_ELEMENTS = 1
BLOCK_SPAN = 1_000_000
ROBUST_SPAN = 5_000_000
BLOCK_ROBUSTNESS = (
    "A SECOND, STRICTLY COARSER BLOCK NULL IS REGISTERED AND REPORTED: the 5-Mb locus block,"
    " (chromosome, gene start // 5,000,000), the same genes in fewer and larger blocks. It is registered"
    " NOW so that the answer's dependence on the span is visible rather than discovered afterwards, and"
    " it CANNOT DECIDE: the outcome name is fixed by the governing unit's p_Q alone, and the 5-Mb null"
    " enters the reading through one of the two REGISTERED clauses of ROBUST_CLAUSES and nowhere else. A"
    " coarser unit can only widen the null, so a 5-Mb null that does not reject while the 1-Mb null does"
    " means the rejection holds at the finer unit alone, and the registered clause says exactly that"
)
BLOCK_COUNTS_EXPECTED = (
    "STATED NOW, from the same address-level pass: the 16,509 genes of the table fall into 2,452 1-Mb"
    " blocks, the largest holding 52 genes (0.31% of the genes) and 385 holding one; and into 596 5-Mb"
    " blocks, the largest holding 158 (0.96%). ON THE 16,496 GENES OF THE 23 ADMITTED CHROMOSOMES, which"
    " is what the null runs on, that is 2,443 1-Mb blocks and 591 5-Mb blocks -- chrY's 13 genes and its"
    " 9 and 5 blocks being the whole of the difference. Every block lies wholly within one chromosome,"
    " because the chromosome is part of the block key. These are COUNTS OF ADDRESSES and no measured"
    " value enters them"
)
MIN_BLOCKS = 100
MAX_BLOCK_SHARE = 0.10
DEGENERACY_RULE = (
    "FIXED NOW. The block null is DEGENERATE, and the lane reads INCONCLUSIVE with no heterogeneity"
    " verdict, if any of: fewer than MIN_CHROMOSOMES = 10 chromosomes are admitted; a group's sum of P^2"
    " is 0 so a fitted scale is undefined; the governing unit yields fewer than MIN_BLOCKS = 100 blocks;"
    " or one block holds more than MAX_BLOCK_SHARE = 0.10 of the admitted genes. The last two are what"
    " stop a block unit so coarse that the null is a statement about one block's noise: MIN_BLOCKS is"
    " magcal's own floor of 100, the number this project already requires for a group to be read at all,"
    " reused rather than chosen, and 0.10 is a tenth of the population. Both are satisfied by the counts"
    " of BLOCK_COUNTS_EXPECTED with room to spare, and both are written so that a FUTURE coarser unit"
    " cannot be slipped in under this registration"
)

# ---- prior exposure, stated twice ---------------------------------------------------------------
PRIOR_EXPOSURE = (
    "STATED PLAINLY, AND NOT A CAVEAT ADDED AFTERWARDS. THIS IS A RE-ANALYSIS OF A MEASUREMENT ALREADY"
    " READ, NOT A FRESH TEST. The same ENCODE polyA plus RNA-seq of the same thirteen reader biosamples,"
    " under the same log2(TPM+1) endpoint and the same no-rule baseline, has already been read by: item 7"
    " (M2b, lane-readercompile); item 7d (M2c, lane-combine2); item 7h (M2f, lane-multrna); the capped-sum"
    " census (lane-capsumrna), which reached NOT_TESTABLE; M2h (lane-magcal), which fitted the one scale;"
    " M2i (lane-freshrna), which replicated it on RNA fresh to this project; lane-scalesplit, which"
    " re-split the same rows by chromosome; lane-recovered; AND lane-scalehet, which is the reading this"
    " lane re-tests. The P and D grouped here ARE 7h's committed rows, so EVERY VALUE THIS LANE READS HAS"
    " BEEN READ BEFORE, NINE TIMES OVER.",
    "AND, SEPARATELY AND MORE SHARPLY: THE SPREAD WAS ALREADY SEEN POST HOC, AND HAS NOW BEEN TESTED ONCE"
    " AND REJECTED. scalesplit computed and committed all 24 per-chromosome scales as post hoc"
    " descriptives -- the range +0.1232 (chr13) to +0.2874 (chr14), the median +0.2176, the quartiles"
    " +0.2041 and +0.2376, no sign change, 10 of 24 inside M2h's interval -- and then scalehet registered"
    " a statistic for them and read NOT CONSISTENT WITH A SINGLE GLOBAL SCALE at p_Q = 0.003996. SO THIS"
    " LANE IS RE-TESTING A PATTERN THAT HAS BEEN LOOKED AT AND ACTED ON, WHICH IS WEAKER STILL THAN"
    " SCALEHET'S OWN TEST, WHICH WAS ALREADY WEAKER THAN A TEST OF AN UNLOOKED-AT PATTERN. Three things"
    " are done about that and none of them removes it. (1) NOT ONE NUMBER OF THE TEST IS RE-CHOSEN: the"
    " statistic, the secondary statistic, the p form, the gene floor, the admitted set, the draw count,"
    " the seed and ALPHA are scalehet's, imported by identity, and the only change is the exchangeable"
    " unit -- so this lane cannot tune its way to either answer by moving a threshold. (2) THE DIRECTION"
    " OF THE EXPECTED EFFECT WAS WRITTEN DOWN FIRST, by scalehet itself: a gene-level null is too narrow"
    " and its p_Q anti-conservative in the REJECTING direction, so a block null is expected to move p_Q"
    " UP, which is the direction that WITHDRAWS this project's committed reading and costs it a finding."
    " (3) BOTH OUTCOMES ARE WRITTEN IN ADVANCE, in OUTCOMES below, in the words they will be reported in."
    " A reading here holds as a statement about THIS measurement's dependence on the chromosome grouping"
    " at these gene counts and this block unit. It is not independent confirmation of anything, it is not"
    " a blind prediction of an unseen measurement, and a withdrawal is a withdrawal of a REJECTION to NOT"
    " ESTABLISHED and never a demonstration that one constant is exactly right",
)
PRIOR_EXPOSURE_CLAUSE = (
    "carried with every reading, always: this RNA is not a fresh measurement and this pattern is not an"
    " unlooked-at one. The same ENCODE polyA plus RNA-seq of the same thirteen reader biosamples under the"
    " same endpoint was read by items 7, 7d and 7h, by the capped-sum census, by M2h, by M2i, by"
    " lane-scalesplit, by lane-recovered and by lane-scalehet, and the P and D grouped here are 7h's"
    " committed rows. AND THE SPREAD TESTED HERE WAS ALREADY SEEN POST HOC BY SCALESPLIT AND HAS ALREADY"
    " BEEN TESTED ONCE AND REJECTED BY SCALEHET, so this is a re-test of a pattern that has been looked at"
    " AND ACTED ON, which is weaker than a test of one that has not. Nothing here is a fresh measurement,"
    " a new reader or an independent confirmation"
)

# ---- what is computed ---------------------------------------------------------------------------
GROUPING = (
    "BY CHROMOSOME, scalehet's grouping unchanged, as 7h's committed table names it and as a property of"
    " the gene's ADDRESS and of no measured value. The GROUPS being compared are the chromosomes; the"
    " BLOCKS are the units the null exchanges. The two are different things and both are addresses: the"
    " question is still whether the per-chromosome fitted scales disperse more than chance, and the only"
    " change is what chance is taken to be. THE CHROMOSOME IS A COARSE AND ARBITRARY REGION UNIT and that"
    " remains a LIMITATION exactly as scalehet wrote it: variation at a finer scale need not show between"
    " chromosomes"
)
STATISTIC = (
    "SCALEHET'S, BY IDENTITY, BOTH OF THEM, so the two nulls are comparable. PRIMARY, Q = sum over"
    " admitted chromosomes of pp_c * (a_c - a_g)^2, the precision-weighted dispersion of the"
    " per-chromosome fitted scales about the global one, Cochran's homogeneity form, where a_c = dp_c /"
    " pp_c is magcal.FIT on that chromosome's rows, dp_c that chromosome's sum of D*P and pp_c its sum of"
    " P^2. The OBSERVED Q is computed by scalehet.observed_statistics, which calls"
    " scalehet.heterogeneity_q -- the same function objects -- and each NULL DRAW's Q is computed by the"
    " same numpy expression scalehet.permutation_null uses for its own draws, character for character,"
    " which is what makes the singleton-block case bit-identical rather than merely equal to within"
    " floating-point summation order; the two forms are the same sum in a different order and agree"
    " within TOLERANCE, which a test holds. Q is NOT compared to"
    " a chi-squared distribution, so no asymptotic, normality or independence-of-pairs claim is made."
    " SECONDARY, R = max a_c - min a_c, the range, by scalehet.range_stat: the statistic the post hoc"
    " sentence used, carried so that the pattern that was looked at is the pattern re-tested. R can never"
    " declare either answer and its only force is MIXED. NO NEW STATISTIC IS INTRODUCED BY THIS LANE"
)
NULL = (
    "THE CHROMOSOME LABELS PERMUTED AMONG THE GENES BY WHOLE LOCUS BLOCKS, HOLDING EVERY CHROMOSOME'S"
    " GENE COUNT EXACTLY AT ITS OBSERVED VALUE. One draw: the blocks are shuffled; their genes are"
    " concatenated in that shuffled order into one sequence, each block's genes staying adjacent; the"
    " sequence is cut at the OBSERVED cumulative per-chromosome gene counts; each piece is one"
    " chromosome's genes for that draw; a scale is refitted per piece and globally; and Q and R are"
    " recorded. Nothing about any gene's data changes -- a gene moves as one block with all thirteen of"
    " its pairs, its own P and its own D -- and only which group it is counted in changes. FIVE THINGS"
    " ABOUT THIS NULL, DECIDED NOW. (1) IT HOLDS THE PER-CHROMOSOME GENE COUNTS EXACTLY, by construction"
    " of the cut, which is the one property that makes it comparable with scalehet's: the only difference"
    " between the two nulls is the exchangeable unit. (2) IT IS SCALEHET'S OWN NULL WITH THE UNIT"
    " REPLACED, to the letter: scalehet shuffles the gene indices and splits at the same cumulative"
    " bounds, and with one gene per block this lane's draw is scalehet's draw BIT-IDENTICALLY at the same"
    " seed, which the harness proves rather than asserts. (3) A CUT CAN FALL INSIDE A BLOCK, so a block"
    " can straddle two chromosomes in a draw; with 23 groups there are 22 cuts, so AT MOST 22 of the"
    " blocks are split in any draw and the rest travel whole. The mean number split per draw is reported."
    " This is the price of holding the gene counts exact, and it is paid in the CONSERVATIVE direction"
    " only in the sense that it slightly weakens the blocking; it cannot manufacture a wider null. (4)"
    " THE GLOBAL FITTED SCALE IS INVARIANT ACROSS DRAWS, by construction: the same genes are redistributed"
    " so the total sums, and therefore the global scale, are the same in every draw -- exactly, but for"
    " the ORDER OF THE FLOATING-POINT SUMMATION, which moves the last bits and is held within TOLERANCE"
    " rather than claimed bit-equal. (5) IT PRICES THE"
    " CORRELATION BETWEEN NEIGHBOURING GENES AND NOTHING ELSE. It does not price a correlation between"
    " distant genes of one chromosome, it does not price gene density, expression level or rules per gene,"
    " and it cannot turn a true regional magnitude into chance or chance into a regional magnitude"
)
PERMUTATIONS = {"count": sh.PERMUTATIONS["count"], "seed": sh.PERMUTATIONS["seed"]}
PERMUTATIONS_WHY = (
    "SCALEHET'S COUNT AND SCALEHET'S SEED, both unchanged: 1,000 draws at seed 20261012. The count is"
    " magcal's PERMUTATIONS count, which fixes the finest p this test can report at 1/1001 = 0.000999 and"
    " is the same granularity scalehet's p was read at. THE SEED IS DELIBERATELY THE SAME ONE: at one gene"
    " per block this lane's draws ARE scalehet's draws, so the element-sharing null reproduces scalehet's"
    " null bit-identically and the harness can show that the machinery differs from scalehet's in the UNIT"
    " and in nothing else. A seed is not re-drawn to suit an answer and the run is run ONCE"
)
P_VALUE = (
    "scalehet.permutation_p, the same function object, which is magcal's one-sided form: p = (1 + draws"
    " with the statistic >= the observed) / (1 + count). It is one-sided because the question is whether"
    " the observed dispersion is LARGER than chance, and it can never be 0, so a p of 0.000999 means 'no"
    " draw of 1,000 reached it' and not 'impossible'"
)
ALPHA = (
    "ALPHA = CONTROL_ALPHA = 0.05, SCALEHET'S OWN LEVEL, imported by identity and NOT re-chosen here, and"
    " its laxness is deliberate in the SAME direction scalehet chose it: at 0.05 the test rejects a single"
    " global scale more readily, which makes 'the heterogeneity is not established' THE HARDER RESULT TO"
    " REACH -- and 'not established' is the result that withdraws this project's committed reading. So the"
    " threshold is set against the withdrawal, not for it. Both p values are ALSO reported against the"
    " STRICTER 0.0167, the complement of the 98.33% description this project reports intervals at"
    " (0.05/3), BECAUSE SCALEHET REPORTED BOTH and its rejection cleared both; the stricter number is a"
    " report and DECIDES NOTHING, and no outcome of this lane turns on it. Nothing is Bonferroni-adjusted"
    " for the two block spans either: the governing unit is fixed in advance by BLOCK_SELECTION_RULE and"
    " the 5-Mb null cannot change the outcome name"
)
INTERVAL = (
    "every null distribution is summarised at its MEDIAN and its 98.33% PERCENTILE RANGE, magcal's"
    " FORM_INTERVAL unchanged and the width this project reports every interval at, through"
    " scalehet.percentile_range -- the same function object. The range is a DESCRIPTION OF THE NULL and"
    " not a confidence interval on anything observed: it says what dispersion this null produces, which"
    " is the honest measure of what this test could and could not have detected. The minimum and maximum"
    " draw are reported beside it so that a p at the granularity floor can be read for what it is"
)
CONTROL = (
    "THE HARNESS IS THE CONTROL, in three parts, and each is a planted check rather than a claim. (1)"
    " SCALEHET'S COMMITTED FIGURES COME BACK: its Q, its R, its global scale, its null median and 98.33%"
    " range for both statistics, and both p values, recomputed through SCALEHET'S OWN functions at"
    " scalehet's own count and seed and compared with data/results/scalehet_result.json. This is run and"
    " committed BEFORE anything else, and if they do not come back the lane stops. (2) THE SINGLETON-BLOCK"
    " EQUIVALENCE: this lane's block null, given one gene per block, reproduces scalehet's null"
    " BIT-IDENTICALLY -- every one of the 1,000 draws' Q and R equal, not merely close -- which is what"
    " shows the unit is the only thing that changed. (3) THE EXACT-COUNT CHECK: every draw of the block"
    " null is verified to give each chromosome exactly its observed gene count, and the check is run over"
    " the draws rather than asserted from the construction"
)
HARNESS = (
    "(1) the table is 7h's: 16,509 genes, 214,617 pairs, and 7h's committed uncalibrated Delta MSE comes"
    " back at scale 1 within TOLERANCE. (2) the global fitted scale from magcal.fit_scale on the whole"
    " table equals scalehet's committed global scale within TOLERANCE. (3) M2h's committed gene-level"
    " numbers are NOT re-bootstrapped here: lane-scalehet's own committed harness"
    " (data/results/scalehet_harness.json) already reproduced them bit-identically and its outcome"
    " HARNESS_REPRODUCES is read and required, not re-run. (4) SCALEHET'S COMMITTED Q, R, GLOBAL SCALE,"
    " NULL SUMMARIES AND BOTH p VALUES ARE REPRODUCED through scalehet's own observed_statistics and"
    " permutation_null at its committed count and seed, and reported BEFORE ANYTHING ELSE OF THIS LANE."
    " (5) the address-level block structure: the element-sharing graph from the compiled programs (how"
    " many elements, how many drive two or more genes of any kind, how many drive two or more of the"
    " table's genes, and the component sizes over the table's genes), the 1-Mb and 5-Mb block counts, the"
    " largest block's share, and the gene-start join reporting zero genes without a start and zero"
    " chromosome disagreements. All of it depends on no measured value. (6) the singleton-block"
    " equivalence and the exact-count check of CONTROL (2) and (3). (7) NO Q, NO R AND NO p OF ANY BLOCK"
    " NULL OF THIS LANE IS COMPUTED IN THE HARNESS: those are the run's, after the registration is"
    " committed. The harness forms the block MEMBERSHIP and checks the draws' SIZES, and deliberately"
    " computes no statistic from a block draw"
)

# ---- the decision procedure, written before any statistic of this lane was computed --------------
WHICH_GOVERNS = (
    "p_Q UNDER THE GOVERNING BLOCK UNIT DECIDES, and nothing else does. The governing unit is fixed by"
    " BLOCK_SELECTION_RULE before any statistic is computed. p_R under that same unit can only produce"
    " MIXED, never a verdict of its own, exactly as in scalehet. The 5-Mb robustness null and the"
    " element-sharing null are REPORTED ALWAYS and DECIDE NOTHING: they enter the reading through the"
    " registered clauses and cannot change the outcome name. scalehet's own gene-level p values are"
    " reproduced as a harness check and are not re-decided"
)
NOT_VACUOUS = (
    "THIS TEST CAN COME OUT EITHER WAY AND BOTH WAYS ARE WRITTEN DOWN. A block null is WIDER than a gene"
    " null exactly to the extent that neighbouring genes are correlated, so if that correlation is small"
    " the two nulls nearly coincide and scalehet's rejection survives almost unchanged; if it is large the"
    " null widens and the rejection does not survive. Neither is guaranteed by the construction. The"
    " planted tests hold both directions: a planted regional slope survives the block null and reads that"
    " the heterogeneity SURVIVES, and a planted single slope with strong within-block correlation reads"
    " NOT ESTABLISHED, so a test that could only say one thing would fail. AND THE DIRECTION THAT COSTS"
    " THIS PROJECT A FINDING IS THE ONE THE ERROR IS KNOWN TO RUN IN: scalehet wrote down in advance that"
    " its p_Q was anti-conservative toward REJECTING, so a withdrawal is the expected correction and is"
    " the outcome this lane was opened to make reachable"
)
OUTCOMES = {
    "SURVIVES_BLOCK_NULL": (
        "THE HETEROGENEITY SURVIVES A LOCUS-BLOCK NULL. With the exchangeable unit changed from the gene"
        " to the {unit} -- {n_blocks} blocks over the {genes} genes of the {chroms} admitted"
        " chromosomes, every chromosome's gene count held exactly at its observed value and nothing else"
        " of lane-scalehet's test altered -- the observed dispersion of the per-chromosome fitted scales"
        " about the global scale {global_scale} is Q = {q} against a block null whose Q has median"
        " {null_q_median} and 98.33% range [{null_q_lo}, {null_q_hi}]: p_Q = {p_q} against ALPHA"
        " {alpha} (and {p_q} against the stricter {stricter}, reported and deciding nothing). The"
        " secondary range statistic agrees: the observed spread {r} against a block null of median"
        " {null_r_median} and 98.33% range [{null_r_lo}, {null_r_hi}], p_R = {p_r}. {robust} SO"
        " SCALEHET'S REJECTION IS NOT EXPLAINED BY THE DEPENDENCE BETWEEN NEIGHBOURING GENES: its"
        " reading NOT CONSISTENT WITH A SINGLE GLOBAL SCALE STANDS, and the limitation it flagged against"
        " itself -- that its null was too narrow because its unit was the gene -- is now priced and does"
        " not account for the result. THAT IS THE WHOLE FINDING AND THIS LANE STOPS HERE: no regional"
        " model is proposed, no per-region scale is fitted, no cause is named and no chromosome or block"
        " is read as a verdict about itself. M2h's one scale +0.2146 and its about 2% of the variance are"
        " unaffected AS A GLOBAL AVERAGE and are not withdrawn; what remains withdrawn, by scalehet and"
        " not by this lane, is the assumption that one constant is the whole of the correction"
    ),
    "WITHDRAWN_NOT_ESTABLISHED": (
        "THE HETEROGENEITY IS NOT ESTABLISHED, AND SCALEHET'S REJECTION IS WITHDRAWN. With the"
        " exchangeable unit changed from the gene to the {unit} -- {n_blocks} blocks over the {genes}"
        " genes of the {chroms} admitted chromosomes, every chromosome's gene count held exactly at its"
        " observed value and nothing else of lane-scalehet's test altered -- the same observed dispersion"
        " Q = {q} about the same global scale {global_scale} sits inside what the block null itself"
        " produces: median {null_q_median}, 98.33% range [{null_q_lo}, {null_q_hi}], p_Q = {p_q} against"
        " ALPHA {alpha} (and against the stricter {stricter}). The secondary range statistic: observed"
        " {r}, block null median {null_r_median}, 98.33% range [{null_r_lo}, {null_r_hi}], p_R = {p_r}."
        " {robust} SO THE SPREAD OF THE PER-CHROMOSOME SCALES IS WHAT THE CORRELATION BETWEEN"
        " NEIGHBOURING GENES PRODUCES AT THESE GENE COUNTS, AND LANE-SCALEHET'S READING -- NOT CONSISTENT"
        " WITH A SINGLE GLOBAL SCALE, p_Q = 0.003996 -- IS WITHDRAWN TO NOT ESTABLISHED. It was"
        " anti-conservative in exactly the direction it flagged against itself: its null's unit was the"
        " gene, the gene null is too narrow, and its p_Q was too small for that reason. SAID PLAINLY AND"
        " NOT SOFTENED: this project should no longer carry 'the magnitude varies between chromosomes' as"
        " a finding. IT IS NOT SHOWN THAT ONE CONSTANT IS THE WHOLE CORRECTION EITHER -- a failure to"
        " reject is not a proof of identity, and the block null's own 98.33% range is the honest measure"
        " of the variation this test could not have seen -- so the standing of the question returns to"
        " UNTESTED-AND-NOT-ESTABLISHED, not to settled. M2h's one scale +0.2146 and its about 2% of the"
        " variance stand as a GLOBAL AVERAGE, exactly as they did before, neither strengthened nor"
        " weakened by this reading. NO REGIONAL MODEL IS PROPOSED AND NONE IS RULED OUT"
    ),
    "MIXED": (
        "MIXED, AND THE HETEROGENEITY IS THEREFORE NOT ESTABLISHED. Under the {unit} null -- {n_blocks}"
        " blocks, {genes} genes, {chroms} chromosomes, every chromosome's gene count held exactly -- the"
        " two registered statistics fall on OPPOSITE sides of ALPHA {alpha}: Q = {q} against a null of"
        " median {null_q_median} and 98.33% range [{null_q_lo}, {null_q_hi}] gives p_Q = {p_q}, while the"
        " range {r} against a null of median {null_r_median} and 98.33% range [{null_r_lo}, {null_r_hi}]"
        " gives p_R = {p_r}. {robust} THE REGISTERED PROCEDURE CLAIMS NEITHER ANSWER WHEN THEY DISAGREE,"
        " and the consequence is stated plainly: scalehet's rejection is NOT CONFIRMED by this re-test,"
        " so it does not stand as established, and this lane does not assert the opposite either. The"
        " disagreement is reported and not resolved, because resolving it by picking the statistic that"
        " agrees with a prior reading is exactly what registering both in advance was meant to prevent."
        " M2h's +0.2146 and its about 2% of the variance stand as a GLOBAL AVERAGE. No regional model is"
        " proposed"
    ),
    "INCONCLUSIVE": (
        "INCONCLUSIVE: {why}. No heterogeneity verdict is given, scalehet's reading is neither confirmed"
        " nor withdrawn by this lane, and no regional model is proposed. M2h's one scale +0.2146 and its"
        " about 2% of the variance stand as a GLOBAL AVERAGE, untouched"
    ),
}
ROBUST_CLAUSES = {
    "AGREES": (
        "The registered coarser 5-Mb locus-block null AGREES with the governing one at ALPHA: there Q ="
        " {q5} against a null of median {null_q5_median} and 98.33% range [{null_q5_lo}, {null_q5_hi}],"
        " p_Q = {p_q5}, the same side of ALPHA, over {n_blocks5} blocks -- so the reading does not turn"
        " on the span of the block."
    ),
    "DISAGREES": (
        "THE REGISTERED COARSER 5-Mb LOCUS-BLOCK NULL DOES NOT AGREE, and that is part of this reading"
        " and not a caveat added to it: over {n_blocks5} blocks, Q = {q5} against a null of median"
        " {null_q5_median} and 98.33% range [{null_q5_lo}, {null_q5_hi}] gives p_Q = {p_q5}, the OTHER"
        " side of ALPHA. A coarser block can only widen the null, so the governing unit's reading holds"
        " AT THE 1-Mb SPAN ALONE and does not hold once a 5-Mb locus block is priced; which span is the"
        " right one for this dependence is not decided here and was not decided in advance, so the"
        " answer must be read as span-dependent."
    ),
}
DECISION = (
    "APPLIED IN THIS ORDER AND IN NO OTHER, and fixed before any statistic of this lane is computed. (1)"
    " Reproduce lane-scalehet's committed Q, R, global scale, null summaries and both p values; if they"
    " do not come back within TOLERANCE the lane stops and reads nothing. (2) Fix the governing block"
    " unit by BLOCK_SELECTION_RULE from the address-level element-sharing count. (3) If the null is"
    " DEGENERATE by DEGENERACY_RULE, read INCONCLUSIVE and stop. (4) Compute Q and R on the admitted"
    " chromosomes -- which are scalehet's observed values and must equal them -- and the governing block"
    " null's p_Q and p_R. (5) If p_Q and p_R fall on opposite sides of ALPHA, read MIXED. (6) Otherwise"
    " if p_Q <= ALPHA read SURVIVES_BLOCK_NULL; if p_Q > ALPHA read WITHDRAWN_NOT_ESTABLISHED. (7) Fill"
    " the registered sentence for that outcome, selecting the ROBUST_CLAUSES clause by whether the 5-Mb"
    " null's p_Q falls on the same side of ALPHA, and append PRIOR_EXPOSURE_CLAUSE. Always. (8) Report"
    " the element-sharing null as a descriptive, whatever it shows. NO OTHER QUANTITY MAY MOVE THE"
    " OUTCOME, and the outcome name is a function of p_Q, p_R, ALPHA and the degeneracy checks alone"
)
FALSIFIER = (
    "WHAT WOULD HAVE TO HAPPEN FOR EACH READING TO BE WRONG, WRITTEN BEFORE EITHER IS SEEN. A"
    " SURVIVES_BLOCK_NULL reading is wrong if the dependence between the genes grouped here is not"
    " captured by a 1-Mb locus block: if the correlation reaches far across a chromosome, or runs between"
    " genes of the same family or the same expression level wherever they sit, a block null of this span"
    " is still too narrow and a rejection could still be an artefact. It is also wrong if the"
    " correlation runs through the PREDICTOR side rather than the locus -- the same compiled strengths"
    " and the same reader calls drive every P -- which no permutation of addresses can price. A"
    " WITHDRAWN_NOT_ESTABLISHED reading is wrong if the block null is too WIDE for the wrong reason: if"
    " the straddling of blocks by the cuts, or one or two very large blocks, inflates the null's Q"
    " independently of any real correlation. BLOCK_COUNTS_EXPECTED and the reported mean number of split"
    " blocks per draw are what let a reader check that, and MAX_BLOCK_SHARE and MIN_BLOCKS are what refuse"
    " the degenerate case in advance. AND EITHER READING IS WRONG, as a claim about the genome, if the"
    " chromosome is the wrong region unit: both nulls compare CHROMOSOMES, and a magnitude that varies"
    " by element, by locus or by expression level need not vary between chromosomes at all"
)
STOP_RULE = (
    "THE RUN IS RUN ONCE. The registration is committed before the run; the run refuses to start unless"
    " every BOUND field of this module matches the committed registration and every pinned input's sha256"
    " still matches the file on disk; it refuses to overwrite its own result; and the reading is the"
    " registered sentence for the outcome the registered procedure returns, with its numbers filled and"
    " no word changed. NO SECOND SPAN, NO SECOND SEED, NO SECOND UNIT AND NO SECOND STATISTIC IS TRIED"
    " AFTER A VALUE IS SEEN. If a defect is found after the reading is committed it is reported BESIDE"
    " the registered reading and corrected only under a new registration. NO REGIONAL MODEL IS FITTED"
    " WHATEVER THE OUTCOME"
)
LIMITATIONS = (
    "THE BLOCK IS 1 Mb OF SEQUENCE, NOT THE DEPENDENCE ITSELF. A 1-Mb locus block prices the correlation"
    " between genes of one neighbourhood. It does not price a correlation that runs further than 1 Mb,"
    " nor one that runs between genes of the same family, the same expression level or the same rule"
    " count wherever they sit on the chromosome. The 5-Mb null is registered as the one sensitivity to"
    " the span and two spans are not a proof of span-independence",
    "THE ELEMENT-SHARING UNIT, WHICH THE BRIEF PREFERRED, PRICES NOTHING HERE, AND THAT IS A FACT ABOUT"
    " THIS POPULATION RATHER THAN A CHOICE: not one of the programs' 440,558 elements drives two of the"
    " table's 16,509 genes, so every element-sharing component is a single gene and that null IS"
    " scalehet's. The dependence that remains, and that this lane prices, is LOCUS PROXIMITY, which is a"
    " DIFFERENT dependence from the one scalehet named -- so a reading here answers the question"
    " scalehet's limitation should have asked rather than the one it did ask",
    "THE DEPENDENCE ON THE PREDICTOR SIDE IS NOT PRICED BY ANY PERMUTATION OF ADDRESSES. Every P comes"
    " from the same compiled programs, the same attributed elements, the same fitted strengths and the"
    " same ENCODE DNase-seq reader calls, and the scale is M2h's. A permutation null moves genes between"
    " labels; it cannot price a systematic error shared by every prediction",
    "A FAILURE TO REJECT IS NOT A PROOF OF IDENTITY. The block null's own 98.33% Q range is reported as"
    " the honest measure of what this test could have detected: a regional variation smaller than that"
    " spread would not be seen here, so WITHDRAWN_NOT_ESTABLISHED returns the question to not established"
    " and never to settled",
    "THE CUTS SPLIT BLOCKS. Holding every chromosome's gene count EXACTLY -- which is what makes this"
    " null comparable with scalehet's -- forces the cut to fall inside a block sometimes, at most 22"
    " times in a draw of 23 groups. A split block is partly exchanged, so the blocking is very slightly"
    " weaker than a design that let the counts float; the mean number split per draw is reported so the"
    " size of that compromise is on the record",
    "THE REGION UNIT BEING COMPARED IS STILL THE CHROMOSOME, which is coarse and biologically arbitrary,"
    " and is scalehet's unit carried unchanged so that this lane re-tests scalehet's reading rather than a"
    " reading of its own. Regional variation by compiled element, by locus, by expression level or by rule"
    " count need not show BETWEEN CHROMOSOMES at all",
    "THE RNA IS A MEASUREMENT ALREADY READ, by items 7, 7d, 7h, the capped-sum census, M2h, M2i,"
    " scalesplit, recovered and scalehet, AND THE SPREAD RE-TESTED HERE WAS ALREADY SEEN POST HOC AND HAS"
    " ALREADY BEEN TESTED ONCE AND REJECTED. Both are in PRIOR_EXPOSURE and both weaken every reading"
    " below what the same test on an unlooked-at pattern in a fresh measurement would carry",
    "A SCALE IS A MAGNITUDE CORRECTION AND NEVER AN EXPLANATION OF EXPRESSION. Whatever this lane reads,"
    " the calibrated compiled rules buy about 2% of the variance on M2h's held-out side as a GLOBAL"
    " AVERAGE, and that number is not improved, revised or reinterpreted by anything here",
)
DESCRIPTIVES = (
    "reported and DECIDING NOTHING: the element-sharing graph's counts and component sizes; the 1-Mb and"
    " 5-Mb block counts, the largest block's gene share and the number of single-gene blocks; the mean"
    " and maximum number of blocks a draw's cuts split; the per-chromosome fitted scales with their"
    " deviation from the global scale, which are scalehet's committed values unchanged and are NOT read"
    " as verdicts about any chromosome; each admitted chromosome's 98.33% band under the block null,"
    " which carries NO correction over the 23 chromosomes and is NOT a test of any one of them; and the"
    " gene-level p values scalehet committed, quoted beside the block-level ones so the SIZE of the"
    " change in p is visible",
    "AND, AS A DESCRIPTIVE AND NOT A VERDICT: the ratio of the block null's median Q to the gene null's"
    " median Q, which is how much the dependence between neighbouring genes widens the null. It is the"
    " quantity scalehet could not price, it is reported whatever the outcome, and it decides nothing",
)

#: the registration fields the run refuses to differ from: registration field -> constant of this module
BOUND = {
    "lane": "LANE",
    "registered": "REGISTERED",
    "form": "FORM",
    "source": "SOURCE",
    "gene_table": "GENE_TABLE",
    "compiled_programs": "COMPILED_PROGRAMS",
    "cells": "CELLS",
    "measured": "MEASURED",
    "no_rule_baseline": "NO_RULE_BASELINE",
    "min_genes": "MIN_GENES",
    "min_chromosomes": "MIN_CHROMOSOMES",
    "control_alpha": "CONTROL_ALPHA",
    "stricter_alpha": "STRICTER_ALPHA",
    "form_interval": "FORM_INTERVAL",
    "form_interval_words": "FORM_INTERVAL_WORDS",
    "tolerance": "TOLERANCE",
    "scalehet_committed": "SCALEHET_COMMITTED",
    "m2h_committed": "M2H_COMMITTED",
    "seven_h_committed": "SEVEN_H_COMMITTED",
    "carried": "CARRIED",
    "question": "QUESTION",
    "why_not_already_answered": "WHY_NOT_ALREADY_ANSWERED",
    "what_is_not_asked": "WHAT_IS_NOT_ASKED",
    "population": "POPULATION",
    "min_genes_rule": "MIN_GENES_RULE",
    "block_unit": "BLOCK_UNIT",
    "block_unit_why_not_elements": "BLOCK_UNIT_WHY_NOT_ELEMENTS",
    "block_selection_rule": "BLOCK_SELECTION_RULE",
    "min_shared_elements": "MIN_SHARED_ELEMENTS",
    "block_span": "BLOCK_SPAN",
    "robust_span": "ROBUST_SPAN",
    "block_robustness": "BLOCK_ROBUSTNESS",
    "block_counts_expected": "BLOCK_COUNTS_EXPECTED",
    "min_blocks": "MIN_BLOCKS",
    "max_block_share": "MAX_BLOCK_SHARE",
    "degeneracy_rule": "DEGENERACY_RULE",
    "prior_exposure": "PRIOR_EXPOSURE",
    "prior_exposure_clause": "PRIOR_EXPOSURE_CLAUSE",
    "grouping": "GROUPING",
    "statistic": "STATISTIC",
    "null": "NULL",
    "permutations": "PERMUTATIONS",
    "permutations_why": "PERMUTATIONS_WHY",
    "p_value": "P_VALUE",
    "alpha": "ALPHA",
    "interval": "INTERVAL",
    "control": "CONTROL",
    "harness": "HARNESS",
    "which_governs": "WHICH_GOVERNS",
    "not_vacuous": "NOT_VACUOUS",
    "outcomes": "OUTCOMES",
    "robust_clauses": "ROBUST_CLAUSES",
    "decision": "DECISION",
    "falsifier": "FALSIFIER",
    "stop_rule": "STOP_RULE",
    "limitations": "LIMITATIONS",
    "descriptives": "DESCRIPTIVES",
}


def as_written(value: Any) -> Any:
    """`value` as json writes it, NESTED AND NOT ONLY AT THE TOP LEVEL.

    json.dump turns every tuple into a list at EVERY depth, so a top-level-only normalisation makes
    the run's word-for-word check fail on a field that round-trips perfectly -- CARRIED, whose values
    are tuples inside a dict. This LOOSENS NOTHING: no number, no string and no key is touched and the
    comparison stays ==; it only stops a tuple that json cannot represent from reading as drift.
    """
    if isinstance(value, (list, tuple)):
        return [as_written(v) for v in value]
    if isinstance(value, dict):
        return {k: as_written(v) for k, v in value.items()}
    return value


def registered_fields() -> dict[str, Any]:
    """The BOUND fields, as the registration writes them and the run re-checks them."""
    out: dict[str, Any] = {}
    for field, name in BOUND.items():
        v = globals()[name]
        out[field] = list(v) if isinstance(v, tuple) else v
    out = {field: as_written(v) for field, v in out.items()}
    return out


# ---- the blocks: this lane's only addition to scalehet's machinery -------------------------------
def gene_starts(path: str = GENE_TABLE) -> dict[str, int]:
    """gene_id without its version -> start, lane-freshrna's committed convention unchanged.

    The version suffix is stripped, `_PAR_Y` rows are skipped and the first occurrence is kept -- the
    three rules `scripts/freshrna_run.py:gene_starts` applies, which the harness checks by importing
    that function and comparing every gene rather than trusting this sentence.
    """
    from pathlib import Path

    out: dict[str, int] = {}
    with Path(path).open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        gi, si = header.index("gene_id"), header.index("start")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[gi].endswith("_PAR_Y"):
                continue
            out.setdefault(f[gi].split(".")[0], int(f[si]))
    return out


def with_starts(rows: list[dict], starts: dict[str, int]) -> list[str]:
    """The gene_ids of `rows` that the gene table has no start for. REFUSES nothing; it reports."""
    return [r["gene_id"] for r in rows if str(r["gene_id"]).split(".")[0] not in starts]


def locus_blocks(rows: list[dict], span: int = BLOCK_SPAN) -> list[list[int]]:
    """BLOCK_UNIT: the row indices of each (chromosome, start // span) block, blocks in sorted order.

    The grouping is `freshrna.locus_clusters`, THE SAME FUNCTION OBJECT, imported by identity; this
    wrapper only refuses a row without a start, which that function would silently read as 0, and fixes
    the block order so a draw is reproducible.
    """
    from genomeos.attribution import freshrna as fr

    missing = [r.get("gene_id", r.get("gene")) for r in rows if r.get("start") is None]
    if missing:
        raise ValueError(f"{len(missing)} rows have no start, the first being {missing[0]!r}")
    clusters = fr.locus_clusters(rows, span)
    return [clusters[k] for k in sorted(clusters)]


def element_blocks(rows: list[dict], shared: dict[str, list[str]]) -> list[list[int]]:
    """The connected components of the genes under element sharing, as row-index blocks.

    `shared` maps an element to the genes it drives; a component joins every gene an element drives,
    including through a gene that is not in `rows`, so a chain through an unmeasured gene still blocks
    the two measured ones together. With no element driving two of `rows`' genes -- which is this
    population's case -- every component is a single gene and this unit IS the gene unit.
    """
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for genes in shared.values():
        for g in genes[1:]:
            union(genes[0], g)
    comp: dict[str, list[int]] = {}
    for i, r in enumerate(rows):
        comp.setdefault(find(str(r["gene"])), []).append(i)
    return [comp[k] for k in sorted(comp)]


def block_sizes(blocks: list[list[int]]) -> dict[str, Any]:
    """DESCRIPTIVES: the block count, the largest block, its share and the single-gene block count."""
    n = sum(len(b) for b in blocks)
    sizes = sorted((len(b) for b in blocks), reverse=True)
    return {
        "blocks": len(blocks),
        "genes": n,
        "largest": sizes[0] if sizes else 0,
        "largest_share": (sizes[0] / n) if n else float("nan"),
        "single_gene_blocks": sum(1 for s in sizes if s == 1),
        "median_size": sizes[len(sizes) // 2] if sizes else 0,
    }


def block_permutation_null(
    dps: list[float],
    pps: list[float],
    blocks: list[list[int]],
    sizes: list[int],
    count: int,
    seed: int,
) -> dict[str, Any]:
    """NULL: the chromosome labels permuted among the genes BY WHOLE BLOCKS, gene counts held exactly.

    One draw shuffles the BLOCKS, concatenates their genes in that order, cuts the sequence at the
    observed cumulative per-chromosome gene counts, refits a scale per piece and globally, and records Q
    and R. With one gene per block in gene order this is `scalehet.permutation_null` BIT-IDENTICALLY at
    the same seed, which is what the harness checks -- so a draw's Q is formed by THAT function's own
    numpy expression, not by `scalehet.heterogeneity_q`, which is the same sum in a different order and
    would differ in the last bits; R is `scalehet.range_stat`, the function object itself.
    """
    import numpy as np

    if sum(sizes) != len(dps):
        raise ValueError(f"the group sizes sum to {sum(sizes)}, not the {len(dps)} genes given")
    flat = [i for b in blocks for i in b]
    if sorted(flat) != list(range(len(dps))):
        raise ValueError("the blocks are not a partition of the genes given")
    dpa, ppa = np.array(dps, dtype=float), np.array(pps, dtype=float)
    parts_of = [np.array(b, dtype=np.int64) for b in blocks]
    bounds = list(np.cumsum(sizes)[:-1])
    rng = np.random.default_rng(seed)
    qs: list[float] = []
    rs: list[float] = []
    per_group: list[list[float]] = [[] for _ in sizes]
    globals_: list[float] = []
    split_blocks: list[int] = []
    for _ in range(count):
        order = rng.permutation(len(parts_of))
        seq = np.concatenate([parts_of[k] for k in order]) if len(parts_of) else np.empty(0, np.int64)
        parts = np.split(seq, bounds)
        ends = set(np.cumsum([len(parts_of[k]) for k in order]).tolist())
        split_blocks.append(sum(1 for b in bounds if b not in ends))
        dpc = np.array([dpa[p].sum() for p in parts])
        ppc = np.array([ppa[p].sum() for p in parts])
        if (ppc == 0.0).any():
            return {"degenerate": "a permuted group's sum of P^2 is 0, so its scale is undefined"}
        a_g = float(dpc.sum() / ppc.sum())
        a = dpc / ppc
        globals_.append(a_g)
        qs.append(float((ppc * (a - a_g) ** 2).sum()))
        rs.append(sh.range_stat(list(a)))
        for i, v in enumerate(a):
            per_group[i].append(float(v))
    return {
        "degenerate": None,
        "q": qs,
        "r": rs,
        "per_group": per_group,
        "global_scales": globals_,
        "split_blocks_per_draw": {
            "mean": sum(split_blocks) / len(split_blocks) if split_blocks else float("nan"),
            "max": max(split_blocks) if split_blocks else 0,
            "cuts": len(bounds),
        },
        "count": count,
        "seed": seed,
        "blocks": len(blocks),
    }


def draw_group_sizes(blocks: list[list[int]], sizes: list[int], count: int, seed: int) -> list[list[int]]:
    """CONTROL (3): the per-group gene counts of `count` draws, to be checked against `sizes` exactly.

    It reproduces the null's draws with the same rng calls and records ONLY the sizes -- no dp, no pp, no
    scale, no Q and no R -- so the harness can run it before the registration without reading a statistic.
    """
    import numpy as np

    parts_of = [np.array(b, dtype=np.int64) for b in blocks]
    bounds = list(np.cumsum(sizes)[:-1])
    rng = np.random.default_rng(seed)
    out: list[list[int]] = []
    for _ in range(count):
        order = rng.permutation(len(parts_of))
        seq = np.concatenate([parts_of[k] for k in order])
        out.append([len(p) for p in np.split(seq, bounds)])
    return out


# ---- the decision --------------------------------------------------------------------------------
def degeneracy(
    n_admitted: int, n_blocks: int, largest_share: float, undefined: str | None = None
) -> str | None:
    """DEGENERACY_RULE: the reason this lane reads INCONCLUSIVE, or None."""
    if undefined:
        return undefined
    if n_admitted < MIN_CHROMOSOMES:
        return (
            f"only {n_admitted} chromosomes carry at least {MIN_GENES} genes, below MIN_CHROMOSOMES"
            f" {MIN_CHROMOSOMES}"
        )
    if n_blocks < MIN_BLOCKS:
        return f"the governing unit yields {n_blocks} blocks, below MIN_BLOCKS {MIN_BLOCKS}"
    if largest_share > MAX_BLOCK_SHARE:
        return (
            f"one block holds {largest_share:.4f} of the admitted genes, above MAX_BLOCK_SHARE"
            f" {MAX_BLOCK_SHARE}"
        )
    return None


def decide_block(p_q: float, p_r: float, degenerate: str | None = None) -> str:
    """DECISION steps (3), (5) and (6), applied exactly as written. p_Q decides; p_R can only MIX."""
    if degenerate:
        return "INCONCLUSIVE"
    q_rejects = p_q <= CONTROL_ALPHA
    r_rejects = p_r <= CONTROL_ALPHA
    if q_rejects != r_rejects:
        return "MIXED"
    return "SURVIVES_BLOCK_NULL" if q_rejects else "WITHDRAWN_NOT_ESTABLISHED"


def robust_clause(p_q: float, p_q_robust: float, **numbers: Any) -> str:
    """DECISION step (7): the registered 5-Mb clause, chosen by which side of ALPHA it falls on."""
    same = (p_q <= CONTROL_ALPHA) == (p_q_robust <= CONTROL_ALPHA)
    return ROBUST_CLAUSES["AGREES" if same else "DISAGREES"].format(p_q5=p_q_robust, **numbers)


def reading(kind: str, **numbers: Any) -> str:
    """The registered sentence for an outcome, numbers filled, nothing else."""
    return OUTCOMES[kind].format(**numbers)


def full_reading(kind: str, **numbers: Any) -> str:
    """DECISION step (7): the registered sentence, then the prior exposure clause, always."""
    return reading(kind, **numbers) + ". " + PRIOR_EXPOSURE_CLAUSE


def quantiles(xs: list[float]) -> dict[str, float]:
    """The descriptive quartiles, magcal's own convention, imported read-only."""
    return mc.quantiles(xs)
