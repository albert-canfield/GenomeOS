# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""M2h (ROADMAP §5 item, lane-magcal, 2026-10-10): the magnitude calibration of the compiled rules
against measured RNA.

Three composition forms have now been read against measured log2(TPM+1) in the thirteen reader cells
on one design: 7d's runtime mean and promoter competition (lane-combine2; registration 917ce25,
result 948cfc9) and 7h's multiplicative composition (lane-multrna; registration 5c620d3, result
90e5c86). All three read "worse" than running no rule at all. Beside them, 7's post-hoc information
term (+0.0089, not registered) beat 992 of 1,000 shuffles while the cost of the shifts was about
eight times as large. Those two facts together say the compiled rules may carry information at the
wrong SIZE: a prediction that points the right way but moves too far pays mean(P^2) and reads worse.

This module tests that reading the only way one number can test it: ONE SCALE PER FORM, fitted on
training genes alone and applied to held-out genes. Nothing else of 7d's or 7h's design moves. The
endpoint, the baseline, the metric, the bootstrap, the interval convention, MIN_GENES and the reading
vocabulary are theirs, carried verbatim. No per-gene and no per-cell free parameter exists: a form's
whole calibration is one real number, and the registration says so.

It re-runs no stage. 7d's and 7h's committed results carry, per form, a genes table of the very P and
D their registered readings were computed from; this lane reads those tables, pinned by sha256,
reproduces each form's committed uncalibrated Delta MSE as its harness, and then splits, fits and
reads. It therefore opens no compiled program, no narrowPeak file and no RNA quantification of its
own, and edits none of 7d's or 7h's code or results.

PRIOR EXPOSURE is stated plainly, in the registration and beside every reading: this RNA is not a
fresh measurement for this project.

In the order written and committed: this code, the registration, then one run. No model request, no
key, no money, 0 AlphaGenome requests.
"""

from __future__ import annotations

import hashlib
from typing import Any

from genomeos.attribution import combine2 as c2
from genomeos.attribution import multrna as mr

REGISTERED = "2026-10-10"
LANE = "lane-magcal"

#: the three forms calibrated, with the committed result each one's P and D are read from
FORMS = ("mean", "competition", "multiplicative")
SOURCES = {
    "mean": "data/results/combine2_result.json",
    "competition": "data/results/combine2_result.json",
    "multiplicative": "data/results/multrna_result.json",
}
#: each form's words, 7d's and 7h's, unchanged
FORM_WORDS = {**c2.FORM_WORDS, **mr.FORM_WORDS}

# ---- what is carried verbatim from 7d and 7h ---------------------------------------------------
CELLS = c2.CELLS
SELECTION = c2.SELECTION
MEASURED = c2.MEASURED
PSEUDOCOUNT = c2.PSEUDOCOUNT
RNA_SELECTION = c2.RNA_SELECTION
NO_RULE_BASELINE = c2.PREDICTIONS["no_rule_baseline"]
MIN_GENES = c2.MIN_GENES
CONTROL_ALPHA = c2.CONTROL_ALPHA
PERMUTATIONS = c2.PERMUTATIONS
CARRIED = {
    "endpoint": MEASURED,
    "rna_selection": RNA_SELECTION,
    "pseudocount": PSEUDOCOUNT,
    "no_rule_baseline": NO_RULE_BASELINE,
    "min_genes": MIN_GENES,
    "cells": list(CELLS),
    "selection": SELECTION,
    "seven_d_control": c2.CONTROL,
    "seven_h_control": mr.CONTROL,
    "note": (
        "endpoint, baseline, population, pseudocount, cells and MIN_GENES are 7d's and 7h's, carried word"
        " for word; this lane defines no population of its own and admits no gene those two did not admit."
        " 7d's and 7h's shuffled-reader control is quoted here because it is part of the design being"
        " calibrated, and its registered p values are carried as context (7d: mean 0.2997, competition"
        " 0.7343; 7h: 1.0); it is NOT this lane's control, because the element-level open/held rows it"
        " permutes are not in the committed genes tables and re-deriving them would be a re-run of 7d and"
        " 7h, not a calibration. This lane's own control is CONTROL below, named differently for that"
        " reason"
    ),
}

# ---- the registered question ---------------------------------------------------------------------
QUESTION = (
    "Each of the three composition forms already read against measured RNA in the thirteen reader cells"
    " (7d's runtime mean, 7d's promoter competition, 7h's multiplicative composition) read WORSE than"
    " running no rule at all, while 7's post-hoc information term beat 992 of 1,000 shuffles. Do the"
    " compiled rules therefore carry information at the wrong SIZE: with ONE scale per form, fitted on"
    " training genes alone and applied to held-out genes, does the scaled run come closer to measured"
    " log2(TPM+1) across the thirteen reader cell types than the no-rule baseline; and what is that"
    " scale?"
)
POPULATION = (
    "per form F, exactly the population of F's committed result, gene for gene and in the committed"
    " order: 7d's mean (4,725 genes), 7d's promoter competition (2,501 genes) and 7h's multiplicative"
    " composition (16,509 genes), each gene with its thirteen (gene, cell) pairs. 7d's POPULATION and 7h's"
    " POPULATION are carried unchanged; this lane applies no filter of its own, admits no gene they did"
    " not admit and drops none they did. Each form's rows are read from the committed genes table of"
    " SOURCES[F] (pinned by sha256), which carries, per gene, the thirteen P and the thirteen D the"
    " committed registered reading was computed from. The three populations are nested only in part and"
    " are not combined: each form is calibrated and read on its own genes, as 7d and 7h read each form on"
    " its own genes"
)
PRIOR_EXPOSURE = (
    "STATED PLAINLY, AND NOT A CAVEAT ADDED AFTERWARDS. This measurement is not fresh for this project."
    " The same ENCODE polyA plus RNA-seq of the same thirteen reader biosamples, under the same"
    " log2(TPM+1) endpoint and the same no-rule baseline, has already been read by: item 7 (M2b,"
    " lane-readercompile, data/results/readergate_result.json) including its post-hoc information term;"
    " item 7d (M2c, lane-combine2) for the runtime's mean and promoter competition; item 7h (M2f,"
    " lane-multrna) for the multiplicative composition; and today's capped-sum census (lane-capsumrna),"
    " which reached NOT_TESTABLE. The P and D this lane calibrates ARE 7d's and 7h's committed rows, so"
    " every value it reads has been read before, and the choice to look for a magnitude at all was made"
    " BECAUSE those readings came out worse than no rule. The held-out split limits that exposure and"
    " does not eliminate it: the split is fixed before any held-out number is read and no held-out gene"
    " enters any fit, but the held-out genes' RNA values were in the earlier readings' populations, the"
    " endpoint and baseline were chosen with those readings in hand, and nothing here is a blind"
    " prediction of an unseen measurement. Any reading from this lane therefore holds as a statement"
    " about THIS measurement and is not independent confirmation of it; a magnitude read here is a"
    " hypothesis for a fresh measurement, not a result on one"
)
SPLIT = (
    "held out when int(sha256(symbol as UTF-8).hexdigest()[:8], 16) is odd, training when even"
    " (scripts/saturable_gate.is_held_out, the convention lane-saturable registered today, taken"
    " unchanged); the split is on the gene SYMBOL, so a gene falls the same side in every form, is fixed"
    " by this registration before any held-out number is read, and depends on no measured value. A gene"
    " is wholly on one side with all thirteen of its pairs: no cell of a held-out gene enters any fit"
)
PREDICTIONS = {
    "form": (
        "per form F and gene G, P_F(G, X) for the thirteen cells X: the centred prediction of F's"
        " committed genes table, exactly as 7d's and 7h's PREDICTIONS['form'] define it (the gated run's"
        " log2 level in X minus its mean over the thirteen cells), rounded to five decimals by the"
        " committed table"
    ),
    "calibrated": (
        "P_F^a(G, X) = a_F * P_F(G, X), one real number a_F per form and nothing else. No per-gene"
        " parameter, no per-cell parameter, no intercept, no transformation of P and no reweighting of"
        " pairs or genes: a form's whole calibration is one number, fitted on the training genes alone"
    ),
    "no_rule_baseline": NO_RULE_BASELINE,
}
FIT = (
    "a_F is the one scale that minimises the squared error over the TRAINING pairs of form F, in closed"
    " form: a_F = (sum over training pairs of D*P_F) / (sum over training pairs of P_F^2); it is"
    " undefined, and the form is NOT_FITTABLE, when the denominator is 0. Equivalently it is the scale"
    " that maximises Delta MSE on the training pairs. It is fitted on training pairs only; no held-out"
    " value enters it. a_F = 1 is the uncalibrated form, which is what 7d and 7h read"
)
METRIC = (
    "Delta MSE(F, a) = mean over the pairs of D^2 - mean of (D - a*P_F)^2 = mean of (2a*D*P_F -"
    " a^2*P_F^2) (combine2.METRIC with the calibrated prediction in place of P), positive when the"
    " scaled form's squared error is below the no-rule baseline's. The registered endpoint is this"
    " quantity on the HELD-OUT pairs at a = a_F, with a_F fitted on the training pairs. Reported beside"
    " it, deciding nothing: the same quantity on the held-out pairs at a = 1 (the uncalibrated form on"
    " the held-out genes) and on the training pairs at a = a_F"
)
#: three forms are judged where 7d judged two, so the per-form interval is Bonferroni-adjusted from 95%
FORM_INTERVAL = (0.8333333333333334, 99.16666666666667)
FORM_INTERVAL_WORDS = "98.33% CI"
BOOTSTRAP = {
    "resamples": c2.BOOTSTRAP["resamples"],
    "seed": c2.BOOTSTRAP["seed"],
    "unit": c2.BOOTSTRAP["unit"],
    "interval": (
        "per form, the 98.33% percentile interval (Bonferroni over the three forms, from 95%; 7d"
        " Bonferroni-adjusted two forms to 97.5% on the same convention). The held-out Delta MSE interval"
        " resamples the HELD-OUT genes with replacement, each carrying its thirteen pairs, with a_F held"
        " fixed at the value fitted on the training genes; the SCALE interval resamples the TRAINING"
        " genes with replacement, each carrying its thirteen pairs, and refits a_F on each resample"
    ),
}
CONTROL = (
    "a within-gene cell permutation, this lane's own and NOT 7d's or 7h's shuffled-reader control (see"
    " CARRIED['note']): for each permutation, each gene's thirteen P values are permuted among its own"
    " thirteen cells by an independent draw of the seeded generator, every gene keeping its own"
    " prediction values, its own range and its own side of the split, and every pair keeping its own D."
    " Only the tie between a prediction and the cell it was made for is broken. The scale is REFITTED on"
    " the permuted training pairs and the permuted held-out Delta MSE is computed, so the control asks"
    " exactly what the registered endpoint asks of a form whose magnitudes survive and whose cell"
    " alignment does not. PERMUTATIONS['count'] permutations at PERMUTATIONS['seed']; control p = (1 +"
    " permutations with held-out Delta MSE >= observed) / (1 + count)"
)
HARNESS = (
    "before any split, scale, held-out value or control is computed, for each of the three forms: (1) the"
    " genes table of SOURCES[F] has F's committed gene count and its committed pair count; (2) Delta"
    " MSE(F, 1) recomputed over all of F's pairs from that table equals the committed delta_mse of F's"
    " result within TOLERANCE; (3) every gene's P and D have thirteen entries and every symbol is"
    " non-empty; (4) no symbol occurs twice in one form. A miss in any of them stops the run"
    " (HARNESS_DEFECT) and nothing is split, fitted or read. (2) is the check that this lane is reading"
    " the same population, the same endpoint and the same baseline that 7d's and 7h's registered"
    " readings were computed on"
)
#: the committed tables round P and D to five decimals, so the recomputed Delta MSE moves in that digit
TOLERANCE = 1e-4

# ---- the decision procedure, written before any held-out number was read -------------------------
#: the supervisor's science ruling of 2026-10-10, carried verbatim and applied as written
WHICH_GOVERNS = (
    "THE HELD-OUT INTERVAL DECIDES. The scale only splits the straddle case:\n"
    '- held-out Delta MSE excludes 0 on the rules\' side -> "better" (shrinking the rules to 0 cannot'
    " produce this);\n"
    '- excludes 0 on the baseline\'s side -> "worse" (the falsifier);\n'
    '- straddles 0, AND the training-set scale interval includes 0 -> "calibration removes the rules: no'
    ' magnitude of the compiled rules improves on no rule" (a negative);\n'
    '- straddles 0, AND the scale interval excludes 0 -> "inconclusive".\n'
    'A negative scale is reported beside whichever reading applies, as "the rules point the wrong way on'
    ' average". Every reading carries the prior-exposure sentence (this RNA was read by items 7, 7d, 7h'
    " and today's capped-sum census)."
)
#: why this four-way table is not vacuous, stated in advance
NOT_VACUOUS = (
    "a scale fitted to near 0 makes the calibrated prediction near-identical to the no-rule baseline, so"
    " its held-out Delta MSE comes out near 0 and the interval straddles 0. A dead scale can therefore"
    " MANUFACTURE A STRADDLE but it cannot manufacture a win: a held-out interval wholly above 0 requires"
    " the scaled prediction to lower the squared error on genes that entered no fit, which a prediction"
    " shrunk to nothing cannot do. 'better' is the one reading the vacuity trap cannot fake, which is why"
    " the held-out interval decides and the scale interval is what tells a dead-scale straddle"
    " (CALIBRATION_REMOVES, a negative) from a live-scale straddle (INCONCLUSIVE)"
)
#: the scale, always reported with its training-set interval, and never a verdict of its own
SCALE_CLAUSES = {
    "fitted": (
        "fitted scale {a} [{ci} {alo}, {ahi}], one number for the whole form, fitted on {tg} training"
        " genes ({tn} gene-cell pairs) and applied unchanged to the held-out genes; the uncalibrated size"
        " 7d and 7h read is 1"
    ),
    "negative": (
        "the rules point the wrong way on average: the scale that best fits the training genes is negative"
        " (scale {a} [{ci} {alo}, {ahi}], fitted on {tg} training genes), so the measured level moves"
        " against the compiled prediction"
    ),
    "not_fittable": (
        "no scale is defined: the training pairs carry no prediction variance (sum of P^2 is 0)"
    ),
}
OUTCOMES = {
    "NOT_TESTABLE": (
        "{form}: not testable at the registered size; {g} held-out genes, fewer than the registered"
        " {min}; no comparison was read ({scale})"
    ),
    "NOT_FITTABLE": "{form}: not fitted and not read; {scale}",
    "BETTER": (
        "{form}: better; on genes held out of the fit, the reader-gated run scaled by the one factor"
        " fitted on the training genes lowered the squared error of measured log2(TPM+1) across the"
        " thirteen reader cell types below the no-rule baseline's, so it is closer to measured RNA than"
        " the no-rule baseline (held-out Delta MSE {d} [{ci} {lo}, {hi}] over {n} gene-cell pairs, {g}"
        " genes; control p = {p}; {scale})"
    ),
    "WORSE": (
        "{form}: worse; on genes held out of the fit, the reader-gated run scaled by the one factor"
        " fitted on the training genes raised the squared error of measured log2(TPM+1) across the"
        " thirteen reader cell types above the no-rule baseline's, so it is further from measured RNA"
        " than the no-rule baseline even at its best single magnitude (held-out Delta MSE {d} [{ci} {lo},"
        " {hi}] over {n} gene-cell pairs, {g} genes; control p = {p}; {scale})"
    ),
    "CALIBRATION_REMOVES": (
        "{form}: calibration removes the rules: no magnitude of the compiled rules improves on no rule;"
        " the held-out interval straddles 0 and the training-set scale interval includes 0, so the"
        " calibration is not distinguishable from switching the compiled rules off (held-out Delta MSE"
        " {d} [{ci} {lo}, {hi}] over {n} gene-cell pairs, {g} genes; control p = {p}; {scale}). This is a"
        " negative result and is not 'no difference detected'"
    ),
    "INCONCLUSIVE": (
        "{form}: inconclusive; the held-out interval straddles 0 while the training-set scale interval"
        " excludes 0, so a magnitude is estimated and this measurement does not resolve whether the"
        " scaled run is closer to or further from measured RNA than the no-rule baseline (held-out Delta"
        " MSE {d} [{ci} {lo}, {hi}] over {n} gene-cell pairs, {g} genes; control p = {p}; {scale})"
    ),
    "HARNESS_DEFECT": (
        "stopped before anything was split, fitted or read: the committed genes tables did not reproduce"
        " their committed readings' numbers ({what})"
    ),
}
DECISION = (
    "1. If any form misses the HARNESS, the run stops before anything is split, fitted or read"
    " (HARNESS_DEFECT), and no scale and no held-out value is computed. 2. Per form, the scale a_F is"
    " fitted on the training genes alone (FIT) and reported with its 98.33% training-gene bootstrap"
    " interval (a_F refitted on each resample) in every reading, whatever the reading is; a training"
    " denominator of 0 is NOT_FITTABLE and the form is read no further. 3. Per form, fewer than MIN_GENES"
    " held-out genes is NOT_TESTABLE and no comparison is read. 4. Otherwise the verdict is WHICH_GOVERNS'"
    " four-way table, applied as written, on the held-out Delta MSE at a = a_F with its 98.33%"
    " held-out-gene bootstrap interval (a_F fixed): lower bound > 0 -> BETTER; upper bound < 0 -> WORSE"
    " (the falsifier); the interval straddles 0 and the scale interval includes 0 ->"
    " CALIBRATION_REMOVES, a negative; the interval straddles 0 and the scale interval excludes 0 ->"
    " INCONCLUSIVE. 5. A negative scale (scale interval wholly below 0) is reported BESIDE whichever of"
    " those applies, in SCALE_CLAUSES['negative']'s words, never instead of it. 6. The control p is"
    " reported in every verdict sentence and changes no verdict; it is this lane's own within-gene cell"
    " permutation (CONTROL), not 7d's or 7h's shuffled-reader control. 7. The combined reading is the"
    " three forms' sentences, in FORMS order, joined by '; ', then PRIOR_EXPOSURE_CLAUSE, always"
)
FALSIFIER = (
    "the reading under test is that the compiled rules carry the right information at the wrong size"
    " ('right information, wrong size'), which is what 7d's, 7h's and the capped sum's readings beside"
    " 7's post-hoc information term suggest. It is FALSIFIED for a form when that form's held-out Delta"
    " MSE interval excludes 0 on the BASELINE's side (upper bound < 0, outcome WORSE): with one scale"
    " fitted freely on the training genes, the form is still further from measured RNA than running no"
    " rule at all, so no single magnitude of it carries the information on held-out genes. An interval"
    " straddling 0 falsifies nothing and is read by the scale: CALIBRATION_REMOVES when the scale"
    " interval includes 0 (a negative, not a null), INCONCLUSIVE when it excludes 0. An interval"
    " excluding 0 on the rules' side (BETTER) is a POSITIVE for that form and is stated as such in the"
    " registered words, with PRIOR_EXPOSURE_CLAUSE attached. NOT_VACUOUS states why 'better' is the"
    " reading a dead scale cannot fake"
)
PRIOR_EXPOSURE_CLAUSE = (
    "carried with every reading: this RNA is not a fresh measurement. The same ENCODE polyA plus RNA-seq"
    " of the same thirteen reader biosamples under the same endpoint was read by items 7, 7d and 7h and"
    " by today's capped-sum census, and the P and D calibrated here are 7d's and 7h's committed rows; the"
    " held-out split limits that exposure and does not eliminate it, so a reading here is a statement"
    " about this measurement and a hypothesis for a fresh one, not independent confirmation"
)
STOP_RULE = (
    "one run, from this registration's committed code, on the files the registration pins by sha256 (7d's"
    " and 7h's registrations and results). No form, population, split rule, scale definition, endpoint,"
    " metric, interval, control, decision, falsifier or wording is changed after a held-out value is"
    " read; a defect found afterwards is reported beside the registered reading and corrected only under"
    " a new registration"
)
LIMITATIONS = (
    mr.LIMITATIONS[0],
    "PRIOR EXPOSURE, repeated as a limitation because it is the largest one: see PRIOR_EXPOSURE. The"
    " held-out genes' RNA values were inside the populations of items 7, 7d, 7h and the capped-sum"
    " census, and the decision to look for a magnitude was made because those readings came out worse"
    " than no rule. No number here is a blind prediction",
    "one scale is the smallest calibration there is and also the weakest: it can only shrink, grow or"
    " flip a form's predictions as a block. It cannot fix a form whose error is in the SHAPE of its"
    " prediction across cells or across genes, and a form that reads FURTHER even at its best scale is"
    " not thereby shown to carry no information — only that no single magnitude of it helps on these"
    " held-out genes",
    "the calibration is a statistical rescaling and not a mechanism: a scale of, say, 0.1 does not say"
    " that an element's measured removal fold is ten times too large, because the same number would"
    " arise from a prediction that is right on a tenth of the genes and empty on the rest. Nothing here"
    " attributes the size to any stage of the compiled pipeline",
    "the three populations differ (4,725, 2,501 and 16,509 genes) and are read separately, as 7d and 7h"
    " read them; no scale is compared between forms and no head-to-head is registered here",
    "the P and D are the committed tables' five-decimal roundings, which is why the harness admits"
    " TOLERANCE rather than exact equality",
    "this lane's control is a within-gene cell permutation and not 7d's or 7h's shuffled-reader control"
    " (CARRIED['note']); their registered control p values are carried as context and are not re-decided"
    " here",
    c2.LIMITATIONS[4],
    c2.LIMITATIONS[5],
)
DESCRIPTIVES = (
    "per form, deciding nothing: the held-out Delta MSE at a = 1 (the uncalibrated form on the held-out"
    " genes) with its 98.33% interval, and the training Delta MSE at a = a_F; the training and held-out"
    " gene and pair counts; mean(D^2), mean(P^2), mean(a_F^2 P^2) and mean(D*P) on each side, so that the"
    " size the calibration keeps is visible as a number; the held-out Delta MSE per reader cell at a ="
    " a_F; the predicted range quantiles before and after scaling",
    "context, copied verbatim and not re-decided: 7d's registered reading, 7h's registered reading and"
    " their control p values",
)
#: the registration fields the run refuses to differ from: registration field -> constant of this module
BOUND = {
    "question": "QUESTION",
    "forms": "FORMS",
    "sources": "SOURCES",
    "form_words": "FORM_WORDS",
    "carried_verbatim": "CARRIED",
    "population": "POPULATION",
    "prior_exposure": "PRIOR_EXPOSURE",
    "split": "SPLIT",
    "predictions": "PREDICTIONS",
    "fit": "FIT",
    "metric": "METRIC",
    "form_interval": "FORM_INTERVAL",
    "form_interval_words": "FORM_INTERVAL_WORDS",
    "bootstrap": "BOOTSTRAP",
    "control": "CONTROL",
    "permutations": "PERMUTATIONS",
    "control_alpha": "CONTROL_ALPHA",
    "min_genes": "MIN_GENES",
    "harness": "HARNESS",
    "tolerance": "TOLERANCE",
    "which_governs": "WHICH_GOVERNS",
    "not_vacuous": "NOT_VACUOUS",
    "scale_clauses": "SCALE_CLAUSES",
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


# ---- the split -----------------------------------------------------------------------------------
def is_held_out(symbol: str) -> bool:
    """SPLIT, lane-saturable's convention: odd first 32 bits of sha256(symbol) are held out."""
    return int(hashlib.sha256(symbol.encode("utf-8")).hexdigest()[:8], 16) % 2 == 1


def split_rows(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    """(training, held out) by SPLIT, on the gene symbol."""
    held = [r for r in rows if is_held_out(r["symbol"])]
    train = [r for r in rows if not is_held_out(r["symbol"])]
    return train, held


# ---- the scale and the metric --------------------------------------------------------------------
def _sums(rows: list[dict], key_p: str = "P") -> tuple[float, float, float, int]:
    """(sum D*P, sum P^2, sum D^2, pairs) over every pair of the rows."""
    dp = pp = dd = 0.0
    n = 0
    for r in rows:
        for p, d in zip(r[key_p], r["D"], strict=True):
            dp += p * d
            pp += p * p
            dd += d * d
            n += 1
    return dp, pp, dd, n


def fit_scale(rows: list[dict], key_p: str = "P") -> float | None:
    """FIT: the one scale minimising squared error over the rows' pairs, or None when undefined."""
    dp, pp, _dd, _n = _sums(rows, key_p)
    return None if pp == 0.0 else dp / pp


def delta_mse(rows: list[dict], scale: float, key_p: str = "P") -> float:
    """METRIC: Delta MSE at one scale over every pair of the rows."""
    dp, pp, _dd, n = _sums(rows, key_p)
    return (2.0 * scale * dp - scale * scale * pp) / n if n else float("nan")


def _gene_terms(rows: list[dict], scale: float, key_p: str = "P") -> tuple[list[float], list[float]]:
    """Per gene, (sum of 2 a P D - a^2 P^2 over its pairs, its pair count): the bootstrap's unit."""
    s, n = [], []
    for r in rows:
        t = 0.0
        for p, d in zip(r[key_p], r["D"], strict=True):
            t += 2.0 * scale * p * d - scale * scale * p * p
        s.append(t)
        n.append(float(len(r["D"])))
    return s, n


def bootstrap_delta(rows: list[dict], scale: float, resamples: int, seed: int) -> tuple[float, float]:
    """BOOTSTRAP: the percentile interval of Delta MSE over gene resamples, the scale held fixed."""
    import numpy as np

    s, n = _gene_terms(rows, scale)
    sa, na = np.array(s), np.array(n)
    rng = np.random.default_rng(seed)
    stats = np.empty(resamples)
    for i in range(resamples):
        idx = rng.integers(0, len(rows), len(rows))
        stats[i] = sa[idx].sum() / na[idx].sum()
    lo, hi = np.percentile(stats, list(FORM_INTERVAL))
    return float(lo), float(hi)


def bootstrap_scale(rows: list[dict], resamples: int, seed: int) -> tuple[float, float]:
    """BOOTSTRAP: the percentile interval of the fitted scale over TRAINING gene resamples, refitted."""
    import numpy as np

    dps, pps = [], []
    for r in rows:
        dp = pp = 0.0
        for p, d in zip(r["P"], r["D"], strict=True):
            dp += p * d
            pp += p * p
        dps.append(dp)
        pps.append(pp)
    dpa, ppa = np.array(dps), np.array(pps)
    rng = np.random.default_rng(seed)
    stats = np.empty(resamples)
    for i in range(resamples):
        idx = rng.integers(0, len(rows), len(rows))
        den = ppa[idx].sum()
        stats[i] = dpa[idx].sum() / den if den else np.nan
    ok = stats[~np.isnan(stats)]
    lo, hi = np.percentile(ok, list(FORM_INTERVAL))
    return float(lo), float(hi)


def permuted_control(train: list[dict], held: list[dict], count: int, seed: int) -> list[float]:
    """CONTROL: per permutation, each gene's P is permuted among its own cells, the scale refitted on
    the permuted training pairs and the permuted held-out Delta MSE computed."""
    import numpy as np

    rng = np.random.default_rng(seed)
    tr_p = np.array([r["P"] for r in train], dtype=float)
    tr_d = np.array([r["D"] for r in train], dtype=float)
    ho_p = np.array([r["P"] for r in held], dtype=float)
    ho_d = np.array([r["D"] for r in held], dtype=float)
    k = tr_p.shape[1]
    out: list[float] = []
    for _ in range(count):
        ti = np.argsort(rng.random(tr_p.shape), axis=1)
        hi_ = np.argsort(rng.random(ho_p.shape), axis=1)
        tp = np.take_along_axis(tr_p, ti, axis=1)
        hp = np.take_along_axis(ho_p, hi_, axis=1)
        den = float((tp * tp).sum())
        a = float((tp * tr_d).sum()) / den if den else 0.0
        n = ho_p.shape[0] * k
        out.append(float((2 * a * hp * ho_d - a * a * hp * hp).sum() / n) if n else float("nan"))
    return out


# ---- the decisions -------------------------------------------------------------------------------
def decide(lo: float, hi: float, scale_lo: float, scale_hi: float) -> str:
    """DECISION, step 4: WHICH_GOVERNS' four-way table, applied as written. The held-out interval
    decides; the scale interval only splits the straddle case."""
    if lo > 0.0:
        return "BETTER"
    if hi < 0.0:
        return "WORSE"
    return "CALIBRATION_REMOVES" if scale_lo <= 0.0 <= scale_hi else "INCONCLUSIVE"


def scale_is_negative(scale_lo: float, scale_hi: float) -> bool:
    """DECISION, step 5: the scale interval lies wholly below 0."""
    return scale_hi < 0.0


def scale_clause(kind: str, **numbers: Any) -> str:
    """The registered scale clause, numbers filled, nothing else."""
    return SCALE_CLAUSES[kind].format(**numbers)


def reading(kind: str, **numbers: Any) -> str:
    """The registered sentence for an outcome, numbers filled, nothing else."""
    return OUTCOMES[kind].format(**numbers)


def combined_reading(per_form: dict[str, dict[str, Any]]) -> str:
    """DECISION, step 7: the three forms' sentences, then the prior-exposure clause, always."""
    clauses = [per_form[f]["reading"] for f in FORMS if f in per_form]
    return "; ".join(clauses) + ". " + PRIOR_EXPOSURE_CLAUSE


def quantiles(xs: list[float]) -> dict[str, float]:
    s = sorted(xs)
    return {str(q): s[min(len(s) - 1, int(q * len(s)))] for q in (0.1, 0.25, 0.5, 0.75, 0.9)} if s else {}
