# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-recovered (ROADMAP §5, area I, 2026-10-12): do the genes scaling makes fittable read better?

`lane-scalefit` (`8469df9` rule and floor committed before counting, `7f65544` census) established, with
NO RNA OPENED, that 7d excluded 12,500 of the 17,983 multi-element genes as `basal_not_positive` -- their
summed single-removal drops reach or exceed the gene's whole intact level -- and that 8,792 of those 12,500
(70.34%) get a positive basal rate back when every drop is multiplied by M2h's fixed constant 0.2146
(`ac5580b`), of which 8,535 pass 7d's full runnable test. Its registered reading claimed NOTHING about
expression and said the RNA test is a separate registration.

THIS MODULE IS THAT REGISTRATION. It asks one question about measured RNA on that recovered population,
and it is built in two committed steps: the population first (POPULATION, `scripts/recovered_population.py`,
`data/results/recovered_population.json`), then the registration (`scripts/recovered_register.py`,
`data/results/recovered_registration.json`), then one run.

The stop list of 2026-09-27 (docs/ROADMAP.md §5) is carried forward: no AlphaGenome request and no
deletion request, no new reader biosample (the 13 covered cells only), nothing in area J, no node-model
work, no BioForge calibration.
"""

from __future__ import annotations

import hashlib
from typing import Any

from genomeos.attribution import combine2 as c2
from genomeos.attribution import magcal as mc
from genomeos.attribution import multrna as mr
from genomeos.attribution import scalefit as sf

LANE = "lane-recovered"
REGISTERED = "2026-10-12"

# ==== step 1: the population, established and committed before the registration ==========================

#: the form whose exclusion defines the population: the MEAN, 7d's runtime default and scalefit's primary
POPULATION_FORM = "mean"

#: the form whose prediction is read against RNA: the MULTIPLICATIVE composition, M2h's arm that read
#: "better". The population form and the prediction form are DIFFERENT and that is deliberate: the
#: scaling recovers genes under the mean's fit, and the question is whether M2h's calibrated
#: multiplicative run -- the one arm of M2h that read better -- is also closer to measured RNA than the
#: no-rule baseline on exactly those genes.
FORM = mr.FORM

POPULATION = (
    "THE GENES SCALING MAKES FITTABLE, re-derived here and not copied from any count. Starting from the"
    " 17,983 multi-element genes 7d/7e found eligible for a fit (`combine2`'s own eligibility, recomputed"
    " from the same 24 compiled non-coding programs), a gene is in this population when all four hold:"
    " (a) under the MEAN form, 7d's fit on the unscaled single-removal observations fails with the reason"
    " `basal_not_positive` (combine2.fit's `if not q < 1.0`), which is the 12,500; (b) under the MEAN"
    " form, 7d's SAME code on observations whose every single-removal effect is multiplied by"
    " scalefit.SCALE (`RHO_i -> 1 - SCALE * (1 - RHO_i)`) returns a positive basal rate, which is"
    " scalefit's PRIMARY `recovered`; (c) that scaled fit also passes 7d's FULL runnable test, i.e. every"
    " fitted strength lies in (0, 1], which is scalefit's SECONDARY `runnable_after`; (d) the gene has a"
    " row in 7h's committed multiplicative genes table (data/results/multrna_result.json, genes_table),"
    " which is what supplies this lane's P and D. (a) to (c) are predictor-side and are recomputed by"
    " scalefit's own functions; (d) is a lookup of 7h's population by SYMBOL ONLY and reads no value."
    " The ENDPOINT is read on the HELD-OUT half of this population under magcal.is_held_out, so that no"
    " gene whose RNA entered M2h's fit of the scale is in the number that decides"
)

#: the genes of (a)+(b)+(c), 7h's table not yet intersected: scalefit's committed `runnable_after`
POPULATION_SOURCE_COUNTS = {
    "source": "data/results/scalefit_census.json, per_form.mean, commit 7f65544 (lane-scalefit)",
    "eligible_for_a_fit": 17_983,
    "excluded_basal_not_positive": 12_500,
    "recovered": 8_792,
    "runnable_after": 8_535,
    "note": (
        "these are scalefit's committed counts, quoted so that the re-derivation has something to"
        " disagree with. THE POPULATION THIS LANE USES IS THE ONE IT RE-DERIVES, and a disagreement with"
        " any of these four numbers stops the lane (POPULATION_HARNESS)"
    ),
}

POPULATION_HARNESS = (
    "before any registration is written and before any RNA value is read: the re-derivation must"
    " reproduce all four of POPULATION_SOURCE_COUNTS' counts exactly, from the same compiled programs and"
    " the same code; the recovered genes must be a subset of the excluded genes and the runnable a subset"
    " of the recovered; every symbol must be non-empty and occur once. A miss in any of them is"
    " POPULATION_DEFECT and nothing else is computed"
)

#: what the overlap with M2h's own population is measured against, and why it must be stated
OVERLAP_QUESTION = (
    "M2h read the multiplicative composition on 16,509 genes, 8,181 training and 8,328 held out, and its"
    " held-out reading of +0.04339 is over those 8,328. The brief that opened this lane stated that the"
    " recovered genes 'were in no previous reading'. THAT IS A CLAIM ABOUT POPULATIONS AND IT IS MEASURED"
    " HERE, NOT ASSUMED: the overlap of the re-derived recovered population with M2h's 16,509, with its"
    " 8,181 training genes and with its 8,328 held-out genes is counted and reported, and whichever way"
    " it comes out is what the registration's exposure sentence says. 7d's mean population (4,725 genes)"
    " is disjoint from the recovered population BY CONSTRUCTION -- a gene 7d excluded is not a gene 7d"
    " fitted -- but 7h's multiplicative population is NOT selected by the mean's basal rate and so may"
    " contain the recovered genes in full"
)


# ==== the population, re-derived ========================================================================


def is_recovered(rho: dict[str, float]) -> tuple[bool, bool, bool]:
    """(excluded, recovered, runnable) for one gene under POPULATION_FORM, by scalefit's own code.

    `excluded` is POPULATION (a), `recovered` is (b) and `runnable` is (c). Nothing here is fitted:
    `scalefit.SCALE` is M2h's constant and `scalefit.recovery` calls `combine2.fit` twice.
    """
    rec = sf.recovery(POPULATION_FORM, rho)
    excluded = rec.before_reason == "basal_not_positive"
    return excluded, rec.recovered, rec.runnable


def held_out_split(symbols: list[str]) -> tuple[list[str], list[str]]:
    """(training, held out) under magcal's SPLIT, carried unchanged: sha256 parity of the symbol."""
    held = sorted(s for s in symbols if mc.is_held_out(s))
    train = sorted(s for s in symbols if not mc.is_held_out(s))
    return train, held


def sha256_of_symbols(symbols: list[str]) -> str:
    """A pin on the population itself, so the registration and the run cannot drift from step 1."""
    h = hashlib.sha256()
    for s in sorted(symbols):
        h.update(s.encode("utf-8"))
        h.update(b"\n")
    return h.hexdigest()


def overlap(mine: set[str], theirs: set[str]) -> dict[str, Any]:
    """The measured relation between two symbol sets, with the word for it."""
    both = mine & theirs
    return {
        "mine": len(mine),
        "theirs": len(theirs),
        "in_both": len(both),
        "mine_only": len(mine - theirs),
        "theirs_only": len(theirs - mine),
        "share_of_mine_in_theirs": (len(both) / len(mine)) if mine else None,
        "disjoint": not both,
        "mine_is_a_subset_of_theirs": mine <= theirs,
    }


def overlap_sentence(ov: dict[str, Any], theirs_name: str) -> str:
    """The plain words for an overlap, used in the registration's exposure sentence."""
    if ov["disjoint"]:
        return (
            f"the recovered population is DISJOINT from {theirs_name} ({ov['mine']} against {ov['theirs']})"
        )
    if ov["mine_is_a_subset_of_theirs"]:
        return (
            f"the recovered population is WHOLLY CONTAINED in {theirs_name}: all {ov['in_both']} of its"
            f" genes are among those {ov['theirs']}"
        )
    return (
        f"the recovered population OVERLAPS {theirs_name}: {ov['in_both']} of its {ov['mine']} genes are"
        f" among those {ov['theirs']} ({100 * ov['share_of_mine_in_theirs']:.2f}%), and"
        f" {ov['mine_only']} are not"
    )


#: the engine constants the population depends on, recorded so a change is visible
CONVENTIONS = {
    "scale": sf.SCALE,
    "scale_rounded": sf.SCALE_ROUNDED,
    "scale_is_refitted": False,
    "population_form": POPULATION_FORM,
    "prediction_form": FORM,
    "theta": c2.THETA,
    "selection": c2.SELECTION,
    "cells": list(c2.CELLS),
    "split_rule": mc.SPLIT,
}


# ==== step 2: the registration, written AFTER the population and BEFORE any RNA value ===================
# Everything below this line was committed before a single P, D or Delta MSE of this population was read.
# The one run that follows re-checks every BOUND field against this module and refuses to differ.

QUESTION = (
    "lane-scalefit showed, with no RNA opened, that 8,535 of the 12,500 genes 7d excluded as"
    " basal_not_positive become fittable when every single-removal drop is multiplied by M2h's fixed"
    " 0.2146. ON THOSE GENES, does the scaled run come closer to measured log2(TPM+1) across the thirteen"
    " reader cell types than the no-rule baseline? One form (the multiplicative composition, M2h's arm"
    " that read better), one scale (M2h's constant +0.2146, NOT refitted), one endpoint, read on the"
    " genes of this population that M2h held out of its own fit"
)

#: 7h's words for the form, carried unchanged
FORM_WORDS = mr.FORM_WORDS[FORM]

# ---- what is carried verbatim from 7d, 7h, M2h and M2i -------------------------------------------------
CELLS = mc.CELLS
MEASURED = mc.MEASURED
NO_RULE_BASELINE = mc.NO_RULE_BASELINE
RNA_SELECTION = mc.RNA_SELECTION
PSEUDOCOUNT = mc.PSEUDOCOUNT
MIN_GENES = mc.MIN_GENES
CONTROL_ALPHA = mc.CONTROL_ALPHA
PERMUTATIONS = mc.PERMUTATIONS
TOLERANCE = mc.TOLERANCE
CARRIED = {
    "endpoint": MEASURED,
    "no_rule_baseline": NO_RULE_BASELINE,
    "rna_selection": RNA_SELECTION,
    "pseudocount": PSEUDOCOUNT,
    "min_genes": MIN_GENES,
    "control_alpha": CONTROL_ALPHA,
    "permutations": PERMUTATIONS,
    "cells": list(CELLS),
    "selection": c2.SELECTION,
    "split": mc.SPLIT,
    "form_words": FORM_WORDS,
    "note": (
        "the endpoint, the baseline, the RNA selection, the pseudocount, the cells, the minimum size, the"
        " control, the bootstrap, the split and the form's words are 7d's, 7h's and M2h's, carried WORD"
        " FOR WORD. This lane defines no endpoint of its own and reads no cell, no experiment and no"
        " quantification those lanes did not. What is new here is ONLY the population: the genes 7d's"
        " mean form could not fit and the scaling makes fittable"
    ),
}

# ---- the form, the prediction and the scale ------------------------------------------------------------
PREDICTIONS = {
    "form": mc.PREDICTIONS["form"],
    "calibrated": (
        "P^a(G, X) = SCALE * P(G, X) with SCALE the CONSTANT +0.2146, M2h's fitted scale for the"
        " multiplicative composition taken from data/results/magcal_result.json"
        " (scales.multiplicative, commit ac5580b) at full precision. NOTHING IS FITTED, REFITTED,"
        " RE-ESTIMATED, TUNED OR SELECTED IN THIS LANE. There is no per-gene parameter, no per-cell"
        " parameter, no intercept and no transformation of P: the whole calibration is one number that"
        " was fixed before this lane existed"
    ),
    "no_rule_baseline": NO_RULE_BASELINE,
}
SCALE = sf.SCALE
SCALE_ROUNDED = sf.SCALE_ROUNDED
SCALE_SOURCE = sf.SCALE_SOURCE
#: M2h's committed TRAINING-gene interval for that scale, carried because the decision table needs one
SCALE_CI = sf.SCALE_CI
NO_REFIT = (
    "THE SCALE IS NOT REFITTED ON THIS POPULATION, AND THAT IS A CHOICE MADE HERE AND STATED, NOT AN"
    " OMISSION. The brief allowed either the fixed constant or a refit with its own held-out split, and"
    " not both. The fixed constant is taken, for three reasons. (1) The question is whether M2h's"
    " magnitude -- one number for the whole form, already fitted and already read -- also brings these"
    " genes closer to measured RNA; a scale refitted on them would answer a different question, namely"
    " whether SOME magnitude helps on them. (2) Refitting would spend the only genes of this population"
    " that M2h did not fit (its 4,009 held-out members) on a fit, leaving the endpoint either on fitted"
    " genes or on a second-order split of 2,000. (3) A fixed scale makes the number directly comparable"
    " with M2h's +0.04339 and M2i's +0.03727, which are the two numbers this reading has to be read"
    " against. The consequence is declared in advance: this lane CANNOT distinguish 'these genes need a"
    " different magnitude' from 'no magnitude helps on these genes', and a WORSE or INCONCLUSIVE reading"
    " here leaves that open (LIMITATIONS)"
)

# ---- the population's RNA rows -------------------------------------------------------------------------
ROWS = (
    "per gene G of the population, the thirteen P and the thirteen D of G's row in 7h's committed"
    " multiplicative genes table (data/results/multrna_result.json, genes_table, pinned by sha256), taken"
    " unchanged and in the committed order. No P and no D is recomputed, recentred, reweighted or"
    " filtered here: they are the same rows 7h's registered reading and M2h's calibration were computed"
    " from, restricted to this population's symbols. A gene of the population with no row there is not"
    " in the population (POPULATION (d))"
)
ENDPOINT_SIDE = (
    "the HELD-OUT half of the population under magcal.is_held_out, carried unchanged: the genes whose"
    " measured RNA did NOT enter M2h's fit of the scale. The training half is reported beside it as a"
    " descriptive and decides nothing, exactly as M2h reports its training Delta MSE beside its held-out"
    " one"
)

# ---- the metric, the interval, the control -------------------------------------------------------------
METRIC = (
    "Delta MSE = mean over the pairs of D^2 - mean of (D - SCALE*P)^2 = mean of (2*SCALE*D*P -"
    " SCALE^2*P^2) (combine2.METRIC, magcal.METRIC and freshrna.METRIC with the fixed scale), positive"
    " when the scaled run's squared error is below the no-rule baseline's. THE REGISTERED ENDPOINT is"
    " this quantity on the HELD-OUT genes of the population at SCALE, with its INTERVAL_WORDS gene"
    " bootstrap interval. Reported beside it and deciding nothing: the same quantity at SCALE on the"
    " training half and on the whole population, and the same quantity at scale 1 (the uncalibrated form"
    " 7h read) on the held-out half"
)
#: M2h's and M2i's width, kept for comparability and not for a multiplicity this lane has
INTERVAL = mc.FORM_INTERVAL
INTERVAL_WORDS = mc.FORM_INTERVAL_WORDS
INTERVAL_REASON = (
    "98.33%, DECIDED NOW AND FOR A STATED REASON. It is M2h's width -- its Bonferroni adjustment of 95%"
    " over its three forms -- and M2i kept it for the same reason this lane does: ONE form and ONE"
    " endpoint are read here, so 95% is the width this lane's own multiplicity would license, and 98.33%"
    " is WIDER than that. Keeping it is therefore the conservative choice, and it is made so that this"
    " number, M2h's +0.04339 [98.33% CI +0.0379, +0.04894] and M2i's +0.03727 [+0.03202, +0.04299] are"
    " read at one width and can be put side by side without an interval conversion. No second interval is"
    " computed and no interval is chosen after a value is read"
)
BOOTSTRAP = {
    "resamples": mc.BOOTSTRAP["resamples"],
    "seed": mc.BOOTSTRAP["seed"],
    "unit": mc.BOOTSTRAP["unit"],
    "interval": (
        "the 98.33% percentile interval over GENE resamples with replacement, each gene carrying its"
        " thirteen pairs, with SCALE held fixed -- M2h's unit, resample count, seed and width"
        " (magcal.bootstrap_delta, called unchanged). No scale interval is bootstrapped here because no"
        " scale is fitted here; the scale interval the decision table needs is M2h's committed training"
        " interval, carried (SCALE_CI)"
    ),
}
CONTROL = (
    "M2h's within-gene cell permutation, called through magcal.permuted_control unchanged"
    " (magcal.CONTROL): for each permutation, each gene's thirteen P values are permuted among its own"
    " thirteen cells by an independent draw of the seeded generator, every gene keeping its own"
    " prediction values, its own range and its own side of the split, and every pair keeping its own D;"
    " only the tie between a prediction and the cell it was made for is broken. PERMUTATIONS['count']"
    " permutations at PERMUTATIONS['seed']; control p = (1 + permutations with held-out Delta MSE >="
    " observed) / (1 + count). It is reported in every reading and changes no verdict"
)
CONTROL_DIFFERENCE = (
    "ONE DIFFERENCE BETWEEN THE CONTROL AND THE ENDPOINT, STATED IN ADVANCE. magcal.permuted_control"
    " REFITS the scale on each permutation's training pairs, while the endpoint takes SCALE fixed. The"
    " function is used as it stands rather than reimplemented, so the control is allowed a fit the"
    " observed reading is not. That can only make the control HARDER to beat, never easier, so it is the"
    " conservative direction; it is named here because an unexplained asymmetry between a reading and its"
    " control is exactly the kind of thing that should not be discovered afterwards"
)

# ---- the harness, run before any value of this population is read --------------------------------------
HARNESS = (
    "before any P, D or Delta MSE of this population is read: (1) data/results/recovered_population.json"
    " is present, its outcome is POPULATION_ESTABLISHED and its population_harness agrees; (2) its"
    " population_sha256 and held_out_sha256 recompute from its own symbol lists, so the population this"
    " run uses is byte-identical to the committed one; (3) 7h's genes table has its committed 16,509"
    " genes and 214,617 pairs, every gene with thirteen P and thirteen D, and no symbol twice; (4) Delta"
    " MSE at scale 1 recomputed over ALL of 7h's pairs equals 7h's committed delta_mse within TOLERANCE,"
    " which is the check that these are the rows 7h's registered reading was computed on; (5) SCALE"
    " equals magcal_result.json's scales.multiplicative exactly and SCALE_CI its scale_cis.multiplicative"
    " exactly; (6) M2h's committed held-out gene count (8,328) and pair count (108,264) are reproduced"
    " from that table and magcal.is_held_out, which is the check that the split is M2h's; (7) every"
    " symbol of this population's held-out half satisfies magcal.is_held_out and none is in M2h's"
    " training set. A miss in any of them stops the run (HARNESS_DEFECT) and nothing is read"
)

# ---- THE EXPOSURE, STATED PLAINLY ----------------------------------------------------------------------
EXPOSURE = (
    "STATED PLAINLY, AS THE FIRST THING THIS REGISTRATION SAYS ABOUT ITS EVIDENCE, AND MEASURED RATHER"
    " THAN ASSERTED. The RNA is the same measurement read by item 7 (M2b), item 7d (M2c), item 7h (M2f),"
    " the capped-sum census (lane-capsumrna), M2h (lane-magcal) and -- for the fresh replication -- M2i"
    " (lane-freshrna): the same ENCODE polyA plus RNA-seq of the same thirteen reader biosamples, the"
    " same log2(TPM+1) endpoint, the same no-rule baseline, and literally the same committed P and D"
    " rows. THE BRIEF THAT OPENED THIS LANE SAID THESE GENES ENTERED NONE OF THOSE READINGS. THAT IS"
    " FALSE, AND data/results/recovered_population.json MEASURED IT BEFORE THIS REGISTRATION WAS"
    " WRITTEN: 7,999 of the 8,535 (93.72%) have a row in 7h's committed multiplicative table, so they"
    " were inside 7h's 16,509-gene reading and inside M2h's 16,509-gene calibration; 3,990 are in M2h's"
    " 8,181 TRAINING genes, so their measured RNA helped fit the constant +0.2146 this lane applies; and"
    " 4,009 are in M2h's 8,328 HELD-OUT genes, which is the population M2h's +0.04339 and M2i's +0.03727"
    " were both read over. The reason is structural: 7h's multiplicative form is not selected by the mean"
    " form's basal rate, so a gene the mean form could not fit may still have been fitted and read under"
    " the multiplicative one. WHAT IS THEREFORE NEW HERE IS THE POPULATION AND NOT THE MEASUREMENT, AND"
    " NOT EVEN THE POPULATION'S MEMBERSHIP OF AN EARLIER READING: this is a SUBPOPULATION READING of"
    " M2h's held-out reading, on the 4,009 of its 8,328 genes that 7d's mean form could not fit and the"
    " scaling makes fittable. It is NOT a fresh test of the endpoint, it is NOT independent of M2h's"
    " reading, and no reading from it may be worded as one. The only genes of this population that no"
    " previous reading of this endpoint touched are the 536 with no row in 7h's table, and they are"
    " excluded here for want of a P and a D"
)
EXPOSURE_CLAUSE = (
    "carried with every reading: the RNA is not a fresh measurement and this population is not a fresh"
    " population. It is the same ENCODE polyA plus RNA-seq of the same thirteen reader biosamples read by"
    " items 7, 7d, 7h, the capped-sum census, M2h and M2i, and the P and D are 7h's committed rows; 7,999"
    " of the 8,535 recovered genes, including all 4,009 the endpoint is read on, were inside M2h's"
    " 16,509-gene calibration, and those 4,009 were inside the 8,328 M2h's +0.04339 and M2i's +0.03727"
    " were read over. This is a SUBPOPULATION reading of M2h's held-out reading, not independent"
    " confirmation of it; what it adds is that the genes 7d could not fit at all are in it"
)

# ---- the decision, written before any value of this population was read --------------------------------
#: the supervisor's four-way table, carried verbatim from M2h and applied as written
WHICH_GOVERNS = mc.WHICH_GOVERNS
WHICH_SCALE_INTERVAL = (
    "THE SCALE INTERVAL THE TABLE'S STRADDLE BRANCH USES IS M2h's COMMITTED TRAINING INTERVAL, SCALE_CI ="
    " [+0.20074, +0.22830], because no scale is fitted in this lane and there is nothing here to"
    " bootstrap. IT EXCLUDES 0, SO THE 'calibration removes the rules' BRANCH IS UNREACHABLE HERE AND"
    " THAT IS DECLARED NOW, BEFORE ANY VALUE IS READ, not discovered afterwards. A straddling held-out"
    " interval therefore reads INCONCLUSIVE. The branch is kept in the table rather than deleted so that"
    " the table carried is the one the supervisor wrote, and so that the reason this lane cannot reach"
    " that branch is on the record: a dead-scale straddle is a thing a lane that FITS a scale can"
    " discover, and a lane that takes the scale as a constant cannot. M2i reached the same conclusion for"
    " the same reason (freshrna.DECISION: 'There is no fourth branch ... nothing is fitted here')"
)
DECISION = (
    "1. If the HARNESS misses, the run stops before any value of this population is read"
    " (HARNESS_DEFECT), and no Delta MSE is computed. 2. Fewer than MIN_GENES held-out genes in the"
    " population is NOT_TESTABLE and no comparison is read. 3. Otherwise the verdict is WHICH_GOVERNS'"
    " four-way table applied as written, through magcal.decide unchanged, on the HELD-OUT Delta MSE at"
    " SCALE with its 98.33% gene bootstrap interval and with SCALE_CI as the scale interval: lower bound"
    " > 0 -> BETTER; upper bound < 0 -> WORSE (the falsifier); a straddle -> INCONCLUSIVE, because"
    " SCALE_CI excludes 0 (WHICH_SCALE_INTERVAL). 4. The reading sentence is magcal.reading's for that"
    " outcome, the registered vocabulary 7d, 7h and M2h all use, with this lane's numbers in it and"
    " nothing else; it is composed by magcal.combined_reading, which appends M2h's own prior-exposure"
    " clause, and then POPULATION_CLAUSE, EXPOSURE_CLAUSE and VARIANCE_CLAUSE are appended in that order"
    " and always. 5. The control p is reported in every verdict sentence and changes no verdict. 6. A"
    " WORSE or an INCONCLUSIVE reading IS REPORTED FIRST AND UNSOFTENED. It would mean the scaling"
    " recovers genes the FIT can accept but not genes the MODEL predicts, which is a sharper finding than"
    " a confirmation and the more useful one for what to build next"
)
VARIANCE_CLAUSE = (
    "about {pct}% of the variance of measured log2(TPM+1) on these genes (Delta MSE divided by the"
    " held-out mean(D^2) = {d2}), carried beside the reading as M2h's and M2i's readings carry theirs"
    " (M2h 2.0%, M2i 2.3%)"
)
POPULATION_CLAUSE = (
    "the population is the {n} genes of the {total} that 7d's mean form excluded as basal_not_positive and"
    " M2h's fixed 0.2146 makes fittable, restricted to those M2h held out of its own fit; a gene is in it"
    " because the FIT can accept it at that scale, which lane-scalefit established and which says"
    " nothing by itself about expression"
)
FALSIFIER = (
    "the reading under test is that the genes scaling makes FITTABLE are also genes the scaled run"
    " PREDICTS: that on the 4,009 held-out recovered genes the multiplicative composition at the fixed"
    " +0.2146 is closer to measured log2(TPM+1) than running no rule, as M2h found on all 8,328 of its"
    " held-out genes. It is FALSIFIED when the held-out Delta MSE interval EXCLUDES 0 ON THE BASELINE'S"
    " SIDE (upper bound < 0, outcome WORSE): the scaling then recovers genes the fit can accept but not"
    " genes the model predicts, M2h's gain does not reach into the population 7d could not fit, and the"
    " recovery is an arithmetic property of the drops and not a gain in prediction. A straddle"
    " (INCONCLUSIVE) falsifies nothing and confirms nothing and is reported as a failure to show the"
    " gain on this population, not as 'no difference'. An interval excluding 0 on the rules' side"
    " (BETTER) is the registered positive, in EXPOSURE's narrow terms and with VARIANCE_CLAUSE attached"
)
STOP_RULE = (
    "ONE run, from this registration's committed code, on the files it pins by sha256. No form, scale,"
    " population, split, endpoint, metric, interval, control, decision, falsifier or wording is changed"
    " after a value of this population is read; a defect found afterwards is reported beside the"
    " registered reading and corrected only under a new registration. There is no post-hoc block in this"
    " registration: anything computed after the reading is committed is labelled post hoc, decides"
    " nothing, and is not required"
)
LIMITATIONS = (
    "PRIOR EXPOSURE, FIRST AND LARGEST: see EXPOSURE. This is a subpopulation reading of M2h's held-out"
    " reading on the same committed P and D rows, not a fresh test. 4,009 of the 4,009 genes the endpoint"
    " is read on were inside the 8,328 M2h's +0.04339 and M2i's +0.03727 were read over. No number here"
    " is a blind prediction",
    "THE SCALE IS NOT REFITTED (NO_REFIT), so this lane cannot tell 'these genes need a different"
    " magnitude' from 'no magnitude helps on these genes'. A WORSE or INCONCLUSIVE reading here leaves"
    " that question open and naming it would need a second registration with its own split",
    "THE TWO FORMS DIFFER. The population is defined by the MEAN form's exclusion and the recovery is"
    " computed at M2h's MULTIPLICATIVE scale applied to the mean's drops, which lane-scalefit already"
    " recorded as a TRANSFER ACROSS FORMS that M2h does not license (scalefit.LIMITATIONS[0]: M2h's own"
    " mean scale read +0.05032). The prediction read against RNA here is the MULTIPLICATIVE one. So"
    " 'the genes scaling makes fittable' is a mean-form statement and 'the scaled run' is a"
    " multiplicative-form prediction, and a reading here does not say that the mean form's scaled fit"
    " predicts anything",
    "RECOVERY IS NOT A PROPERTY OF THE PREDICTION READ. A gene is in this population because the mean"
    " form's basal rate is positive again at 0.2146. Its multiplicative P was computed by 7h without"
    " reference to that, so this lane is asking whether a population selected on one form's fittability"
    " behaves differently under another form's prediction, and a difference from M2h's 8,328 could come"
    " from anything that distinguishes these genes -- their rule counts (median 15 for the recovered"
    " against 58 for the 3,708 that stay unfittable), their expression range, their loci -- and not from"
    " the recovery as such",
    "536 of the 8,535 recovered-and-runnable genes have no row in 7h's table and are excluded for want of"
    " a P and a D. They are the only genes of this population no previous reading of this endpoint"
    " touched, so the one part of it that would have been a fresh population is exactly the part that"
    " cannot be read",
    "the P and D are 7h's committed tables' five-decimal roundings, which is why the harness admits"
    " TOLERANCE rather than exact equality",
    "the gene-level bootstrap treats genes as independent. Genes that share compiled elements or a locus"
    " are not, so the interval may be too narrow; M2i's post-hoc 1-Mb locus-cluster bootstrap moved its"
    " interval very little ([+0.03216, +0.04302] against [+0.03202, +0.04299]) but that is M2i's"
    " population and not this one",
    "the control refits where the endpoint does not (CONTROL_DIFFERENCE), in the conservative direction",
)
DESCRIPTIVES = (
    "deciding nothing: the held-out Delta MSE at scale 1 (the uncalibrated form 7h read) with its"
    " interval; the Delta MSE at SCALE on the training half and on the whole population; the gene and"
    " pair counts of each side; mean(D^2), mean(P^2), mean(SCALE^2 P^2) and mean(D*P) on each side; the"
    " held-out Delta MSE per reader cell; and M2h's and M2i's committed readings quoted verbatim beside"
    " this one, because the only honest way to read a subpopulation number is against the population's"
    " number"
)

#: the registration fields the run refuses to differ from: registration field -> constant of this module
BOUND = {
    "question": "QUESTION",
    "form": "FORM",
    "form_words": "FORM_WORDS",
    "population_form": "POPULATION_FORM",
    "population": "POPULATION",
    "population_harness": "POPULATION_HARNESS",
    "overlap_question": "OVERLAP_QUESTION",
    "carried_verbatim": "CARRIED",
    "predictions": "PREDICTIONS",
    "scale": "SCALE",
    "scale_source": "SCALE_SOURCE",
    "scale_ci": "SCALE_CI",
    "no_refit": "NO_REFIT",
    "rows": "ROWS",
    "endpoint_side": "ENDPOINT_SIDE",
    "metric": "METRIC",
    "interval": "INTERVAL",
    "interval_words": "INTERVAL_WORDS",
    "interval_reason": "INTERVAL_REASON",
    "bootstrap": "BOOTSTRAP",
    "control": "CONTROL",
    "control_difference": "CONTROL_DIFFERENCE",
    "permutations": "PERMUTATIONS",
    "control_alpha": "CONTROL_ALPHA",
    "min_genes": "MIN_GENES",
    "harness": "HARNESS",
    "tolerance": "TOLERANCE",
    "exposure": "EXPOSURE",
    "exposure_clause": "EXPOSURE_CLAUSE",
    "which_governs": "WHICH_GOVERNS",
    "which_scale_interval": "WHICH_SCALE_INTERVAL",
    "decision": "DECISION",
    "variance_clause": "VARIANCE_CLAUSE",
    "population_clause": "POPULATION_CLAUSE",
    "falsifier": "FALSIFIER",
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


def variance_clause(delta: float, mean_d2: float) -> str:
    """VARIANCE_CLAUSE, filled: the share of the measured variance the gain is."""
    return VARIANCE_CLAUSE.format(pct=round(100.0 * delta / mean_d2, 1), d2=round(mean_d2, 4))


def population_clause(held: int, total: int) -> str:
    """POPULATION_CLAUSE, filled."""
    return POPULATION_CLAUSE.format(n=held, total=total)


def full_reading(kind: str, delta: float, mean_d2: float, held: int, total: int, **numbers: Any) -> str:
    """DECISION step 4: M2h's sentence, M2h's prior-exposure clause, then this lane's three clauses.

    `magcal.reading` and `magcal.combined_reading` are called unchanged, so the verdict vocabulary is
    the one 7d, 7h and M2h all use. Nothing is added before them and nothing is removed from them.
    """
    sentence = mc.reading(kind, **numbers)
    composed = mc.combined_reading({FORM: {"reading": sentence}})
    return (
        f"{composed} {population_clause(held, total)}. {EXPOSURE_CLAUSE}. {variance_clause(delta, mean_d2)}."
    )
