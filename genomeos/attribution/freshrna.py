# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""M2i (ROADMAP §5 item, lane-freshrna, 2026-10-10): the fitted magnitude of M2h tested, unchanged,
on measured RNA of the same thirteen reader cell types that this project has never read.

M2h (`ac5580b`, registration `fbce775`, code `2b0b604`) read the multiplicative composition BETTER
than the no-rule baseline once one scale was fitted on training genes and applied to held-out genes:
held-out Delta MSE +0.04339 [98.33% CI +0.0379, +0.04894] at a fitted scale of +0.2146. It carried
one caveat, in its own words: "this RNA is not a fresh measurement ... the held-out split limits that
exposure and does not eliminate it, so the reading is a statement about this measurement and a
hypothesis for a fresh one, not independent confirmation".

This lane runs that fresh test, and only that. Everything that could be chosen after seeing a value
is already fixed in a committed file: the form is 7h's multiplicative composition unchanged, the
prediction rows P are 7h's committed table unchanged, and THE SCALE IS THE CONSTANT +0.2145911381525326
TAKEN FROM `ac5580b`. NOTHING IS FITTED IN THIS LANE: `fit_scale` exists only to refuse, and the run
calls it nowhere. The one thing that is new is the measured side D, recomputed from ENCODE RNA-seq
gene quantifications of the same thirteen biosample ontology terms whose values no file, script,
result or document of this project has opened, proved at the accession level (experiment, file and
library) against the record, not against memory.

NO READER BIOSAMPLE IS ADDED. `context_evidence.READER_TERMS` is imported and untouched; the reader
stays at its thirteen cells. What is fresh is the MEASUREMENT, not the cell type, and the registration
says in advance exactly how narrow that makes a confirmation.

In the order written and committed: this code, the sources (metadata and bytes only, no value read),
the registration, then one run. No model request, no key, no money, 0 AlphaGenome requests.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

from genomeos.attribution import combine2 as c2
from genomeos.attribution import context_evidence as ce
from genomeos.attribution import magcal as mc
from genomeos.attribution import multrna as mr
from genomeos.attribution import readergate as rg

REGISTERED = "2026-10-10"
LANE = "lane-freshrna"
FORM = "multiplicative"
FORM_WORDS = mr.FORM_WORDS["multiplicative"]

#: where the fresh quantifications are cached, and the record of them (written before the registration)
RNA_DIR = Path("data/cache/freshrna/rna")
SOURCES = Path("data/results/freshrna_sources.json")
#: the committed result the predictions P are read from, and the committed result the scale is taken from
P_SOURCE = "data/results/multrna_result.json"
SCALE_SOURCE = "data/results/magcal_result.json"
#: ac5580b's fitted scale for the multiplicative composition, as a FIXED CONSTANT
SCALE = 0.2145911381525326
SCALE_ORIGIN = (
    "+0.2145911381525326 (+0.2146 to four places), the scale M2h fitted on its 8,181 TRAINING genes and"
    " committed in data/results/magcal_result.json at `ac5580b` (per_form.multiplicative.scale). It"
    " enters this lane as a literal constant of this module and is not recomputed, re-estimated,"
    " rounded, refined or chosen from anything measured here"
)

# ---- NOTHING IS FITTED HERE ----------------------------------------------------------------------
NOTHING_IS_FITTED = (
    "NOTHING IS FITTED IN THIS LANE. There is no free parameter of any kind: no scale, no intercept, no"
    " per-gene and no per-cell parameter, no threshold, no weighting and no transformation is estimated"
    " from the fresh RNA or from anything else. The single number the comparison needs, the scale, is the"
    " constant SCALE taken from `ac5580b` (SCALE_ORIGIN). The code REFUSES TO FIT: `fit_scale` raises"
    " NoFitError when called, and the run calls it nowhere, so a fit cannot be introduced by accident"
    " and cannot be introduced quietly. Because nothing is fitted, no training/held-out split is needed"
    " to protect the fresh measurement: no fresh value enters any estimate, so every gene of the"
    " population is out of sample by construction. The split is nevertheless used to DEFINE the primary"
    " population (POPULATION), for the one reason that M2h's +0.04339 was read on its held-out genes and"
    " the primary reading here must be the same genes, so the two numbers are comparable"
)


class NoFitError(RuntimeError):
    """Raised whenever anything in this lane is asked to fit a parameter. See NOTHING_IS_FITTED."""


def fit_scale(*_args: Any, **_kwargs: Any) -> float:
    """Refuse. This lane fits nothing; the scale is the constant SCALE from `ac5580b`."""
    raise NoFitError(NOTHING_IS_FITTED)


# ---- what is carried verbatim from 7d, 7h and M2h ------------------------------------------------
CELLS = rg.CELLS
READER_TERMS = ce.READER_TERMS
SELECTION = c2.SELECTION
MEASURED = c2.MEASURED
PSEUDOCOUNT = c2.PSEUDOCOUNT
NO_RULE_BASELINE = c2.PREDICTIONS["no_rule_baseline"]
MIN_GENES = c2.MIN_GENES
PERMUTATIONS = mc.PERMUTATIONS
CONTROL_ALPHA = mc.CONTROL_ALPHA
#: M2h's interval, kept so that the number here and ac5580b's +0.04339 are read at the same width
FORM_INTERVAL = mc.FORM_INTERVAL
FORM_INTERVAL_WORDS = mc.FORM_INTERVAL_WORDS
CARRIED = {
    "endpoint": MEASURED,
    "pseudocount": PSEUDOCOUNT,
    "no_rule_baseline": NO_RULE_BASELINE,
    "min_genes": MIN_GENES,
    "cells": list(CELLS),
    "reader_terms": list(READER_TERMS),
    "selection": SELECTION,
    "form_words": FORM_WORDS,
    "seven_h_prediction": mr.PREDICTIONS["form"],
    "m2h_split": mc.SPLIT,
    "m2h_metric": mc.METRIC,
    "note": (
        "the endpoint y(G, X) = log2(TPM(G, X) + 1) and the within-gene centring, the no-rule baseline,"
        " the thirteen reader cells and their ontology terms, MIN_GENES, the bootstrap unit and seed, the"
        " interval width, the permutation count and seed, the form and the reading vocabulary are 7d's,"
        " 7h's and M2h's, carried word for word. READER_TERMS is imported from context_evidence and is"
        " untouched: NO READER BIOSAMPLE IS ADDED, and the reader's thirteen cells and its open/held"
        " calls are exactly those of items 7, 7d, 7h and M2h. Only the MEASUREMENT behind D changes"
    ),
}

# ---- the registered question ---------------------------------------------------------------------
QUESTION = (
    "M2h read the multiplicative composition closer to measured RNA than the no-rule baseline once the"
    " compiled rules were scaled by one fitted factor (+0.2146): held-out Delta MSE +0.04339 [98.33% CI"
    " +0.0379, +0.04894]. That reading was taken on RNA items 7, 7d and 7h had already read. Applying the"
    " SAME form and the SAME scale as a fixed constant, with nothing fitted, does the scaled run still"
    " come closer to measured log2(TPM+1) across the thirteen reader cell types than the no-rule baseline"
    " when the measured side is a DIFFERENT ENCODE RNA-seq measurement of those same thirteen cell"
    " types, one whose values this project has never opened?"
)

# ---- the claim's scope, pre-stated and narrow -----------------------------------------------------
SCOPE = (
    "PRE-STATED AND DELIBERATELY NARROW. This is a fresh MEASUREMENT of the SAME thirteen cell types, so"
    " the only thing a 'better' reading here adds to M2h is that the gain SURVIVES FRESH MEASUREMENT"
    " NOISE: it REPLICATES ON A FRESH MEASUREMENT OF THE SAME CELLS, and nothing wider. Written in"
    " advance, a 'better' reading here would NOT license any of the following: (1) that the gain or the"
    " scale TRANSPORTS TO NEW CELL TYPES -- none is added here, the stop list forbids adding one, and no"
    " cell type outside the reader's thirteen is measured; (2) that 'the model predicts expression' --"
    " the registered words stay 'closer to measured RNA than the no-rule baseline on this measurement',"
    " and the effect is about 2% of the variance (SIZE); (3) independence from the ELEMENTS, the"
    " ATTRIBUTION or the READER -- the predictions P are 7h's committed rows, built from the same"
    " compiled programs, the same attributed elements and the same ENCODE DNase-seq reader calls as"
    " before, and none of that is re-measured here; (4) independence from the compiled programs' own"
    " provenance, including AlphaGenome's training on ENCODE and GTEx RNA-seq of these same biosamples;"
    " (5) that the scale is a mechanism, a dosage or a correction of any stage of the pipeline. What it"
    " would license is exactly one sentence: the magnitude M2h fitted on one measurement of these cells"
    " still beats no rule on a different measurement of the same cells"
)
SIZE = (
    "carried with every 'better' reading, in the supervisor's arithmetic of 2026-10-10 and not softened:"
    " M2h's gain is about 2% OF THE VARIANCE. mean(D*P) = +0.2013 and mean(P^2) = +0.9337 on M2h's"
    " held-out side give the best scale 0.2013/0.9337 = 0.216 (against the fitted 0.2146) and a gain of"
    " mean(D*P)^2/mean(P^2) = 0.0434, which against mean(D^2) = 2.152 is 2.0% of the measured variance,"
    " a correlation of about r = 0.14 between the rules' predicted and the measured deviation. M2i is a"
    " NOTABLE result if it confirms, not a breakthrough: the same 2% figure is computed here on the fresh"
    " measurement and reported beside the reading"
)

# ---- the measurement: RNA this project has never read ---------------------------------------------
RNA_ASSAYS = ("polyA plus RNA-seq", "total RNA-seq")
RNA_SELECTION = (
    "per reader biosample ontology term (READER_TERMS, unchanged, thirteen terms, none added), from the"
    " ENCODE portal: every released GRCh38 'gene quantifications' TSV with genome_annotation V29 (the"
    " ENCODE4 uniform RNA-seq pipeline, RSEM), on the biosample's ontology term name, from a human"
    " experiment the portal lists as not perturbed (perturbed=false) and that carries no ERROR audit, and"
    " whose file biosample summary names no subcellular fraction, treatment, disease or genetic"
    " modification (none of 'fraction', 'treated', 'with ', 'genetically', 'transfected', 'knockdown') --"
    " item 7's RNA_SELECTION, clause for clause -- with exactly THREE changes, each metadata-only and"
    " each stated here before any value was read. (1) EVERY EXPERIMENT THIS PROJECT HAS READ IS EXCLUDED:"
    " a candidate is dropped if its experiment accession, any of its gene-quantification file accessions,"
    " any of its replicate library accessions or any of its replicate biosample accessions appears in the"
    " READ set (READ_SET_RULE). (2) THE ASSAY: 'polyA plus RNA-seq' where an unread qualifying experiment"
    " of that term exists, and 'total RNA-seq' where none does; item 7 took polyA plus for all thirteen,"
    " and for four terms (IMR-90, astrocyte, hepatocyte, cardiac muscle cell) its single qualifying polyA"
    " plus experiment is the one it read, so there is no unread polyA plus one. The assay is recorded per"
    " cell. (3) THE TIE-BREAK IS PREFIXED by one preference: of the unread qualifying experiments of a"
    " term, those sharing NO replicate biosample accession with the READ set come first; within that,"
    " item 7's tie-break is unchanged (the most qualifying files, then the latest date_released, then the"
    " lowest accession). The gene's TPM in the biosample is the mean of its TPM over the chosen"
    " experiment's qualifying files (its replicates), as in item 7. Chosen from portal metadata alone; no"
    " value of any file was read to choose, and the whole choice is committed in"
    " data/results/freshrna_sources.json, pinned by sha256, BEFORE this registration"
)
READ_SET_RULE = (
    "the READ set is built FROM THE RECORD, mechanically, not from memory: (a) every ENCODE accession"
    " (ENCSR, ENCFF, ENCLB, ENCBS) occurring in any git-tracked file of this repository, plus every one"
    " occurring in any of the last 3,000 commit messages, plus every one in the untracked"
    " data/knowledge/rna/tracks_*.json measured-RNA track records; MINUS (b) the accessions that occur"
    " ONLY as a REJECTED CANDIDATE of item 7's own selection, that is as a key of"
    " readergate_rna_sources.json's per-biosample `candidates` map other than the experiment it chose, or"
    " as a key of its `files_excluded` map. (b) is subtracted because item 7 enumerated those accessions"
    " from portal metadata and downloaded and opened none of them (that result's own words:"
    " values_read = 'none: each file was hashed and its header line read; no TPM value was opened', and"
    " only the chosen experiment's files were fetched at all), so their presence in the record is a"
    " metadata mention and not exposure to a value. Everything in the READ set is treated as read. The"
    " two lists -- the READ set and this lane's own accessions -- are committed in the sources file and"
    " in the result, and their intersection is shown"
)
NON_OVERLAP = (
    "PROVED AT ACCESSION LEVEL, not by name, and required before the registration: the chosen thirteen"
    " experiments share NO experiment accession, NO gene-quantification file accession and NO REPLICATE"
    " LIBRARY ACCESSION with the twenty-seven quantifications of the thirteen experiments items 7, 7d, 7h"
    " and M2h read (nor with any other accession of the READ set). Every library is a separate library"
    " preparation and a separate sequencing run, so no replicate is shared. BIOSAMPLE accessions are a"
    " separate and weaker matter and are reported honestly: eleven of the thirteen terms are served by an"
    " experiment that shares no ENCODE biosample accession with the read set either; for TWO --"
    " hepatocyte and cardiac muscle cell -- the chosen experiment draws on the same ENCODE biosample as"
    " the experiment item 7 read, because the portal carries no other unread qualifying experiment of"
    " those terms under the clauses above, and adding a biosample is forbidden. For those two cells the"
    " fresh measurement is therefore a fresh LIBRARY PREPARATION AND SEQUENCING RUN of the same"
    " biological material, not a fresh biological sample; its VALUES are still ones this project has"
    " never read, which is the exposure this lane exists to remove. The per-cell table is in the result"
)
WHAT_IS_FRESH = (
    "the MEASURED side only. D is recomputed from the fresh quantifications; everything on the predictor"
    " side is unchanged and is NOT fresh: the compiled programs, the attributed elements, the fitted"
    " strengths, the reader's ENCODE DNase-seq open/held calls, the thirteen cell types and the"
    " prediction rows P are 7h's committed ones, and the scale is M2h's fitted constant. This test is"
    " therefore independent of the EARLIER RNA and of nothing else"
)

# ---- the population, fixed before any fresh value was read ----------------------------------------
POPULATION = (
    "PRIMARY: the 8,328 genes of M2h's HELD-OUT side of the multiplicative composition -- the genes of"
    " 7h's committed population (data/results/multrna_result.json's genes table, 16,509 genes, pinned by"
    " sha256) whose symbol is held out under M2h's committed split, int(sha256(symbol as"
    " UTF-8).hexdigest()[:8], 16) odd -- intersected with the genes the fresh measurement resolves, which"
    " is item 7's POPULATION clause (d) applied to the fresh tables unchanged: the gene's gene_id without"
    " its version must be on exactly one row of every one of the thirteen chosen experiments' tables"
    " (rows suffixed _PAR_Y not counted), and present in all of the chosen experiment's files. Each gene"
    " contributes its thirteen (gene, cell) pairs. No filter on expression, no filter on any measured"
    " value, and no gene is admitted that 7h did not admit. M2h's held-out side is the primary population"
    " for ONE reason: it is the population ac5580b's +0.04339 was read on, so the two numbers are"
    " comparable, and the scale was fitted without any of these genes. The split, the gene list and this"
    " rule are fixed by this registration and by the two committed results it reads, before any fresh"
    " value is opened. REPORTED BESIDE, deciding nothing: the same quantity on M2h's TRAINING genes"
    " (whose OLD RNA the scale was fitted on, and whose FRESH RNA is still unread) and on all of 7h's"
    " genes. Genes that drop out for want of a fresh row are counted and reported as attrition"
)

# ---- the prediction, the metric, the endpoint -----------------------------------------------------
PREDICTIONS = {
    "form": mr.PREDICTIONS["form"],
    "scaled": (
        "P*(G, X) = SCALE * P(G, X) with SCALE the constant +0.2145911381525326 of SCALE_ORIGIN and"
        " P(G, X) the centred prediction of 7h's committed genes table, unchanged and rounded to five"
        " decimals by that table. One multiplication by a committed constant; nothing else"
    ),
    "no_rule_baseline": NO_RULE_BASELINE,
}
MEASURED_FRESH = (
    "y'(G, X) = log2(TPM'(G, X) + PSEUDOCOUNT), TPM'(G, X) the mean over the chosen fresh experiment's"
    " qualifying files of G's TPM in cell X (RNA_SELECTION). D'(G, X) = y'(G, X) - mean over the thirteen"
    " cells X' of y'(G, X'): exactly MEASURED, with the fresh quantification in place of item 7's"
)
METRIC = (
    "Delta MSE = mean over the pairs of D'^2 - mean of (D' - SCALE*P)^2 = mean of (2*SCALE*D'*P -"
    " SCALE^2*P^2) (combine2.METRIC and magcal.METRIC with the fresh measured deviation and the fixed"
    " scale), positive when the scaled run's squared error is below the no-rule baseline's. The"
    " REGISTERED ENDPOINT is this quantity on the PRIMARY population at SCALE, with its"
    " FORM_INTERVAL_WORDS gene bootstrap interval"
)
BOOTSTRAP = {
    "resamples": c2.BOOTSTRAP["resamples"],
    "seed": c2.BOOTSTRAP["seed"],
    "unit": c2.BOOTSTRAP["unit"],
    "interval": (
        "the 98.33% percentile interval over gene resamples with replacement, each gene carrying its"
        " thirteen pairs, with SCALE fixed -- M2h's unit, resample count, seed and width, kept unchanged"
        " so that this number and ac5580b's +0.04339 [98.33% CI +0.0379, +0.04894] are read at the same"
        " width. The width is M2h's Bonferroni adjustment of 95% over its three forms; one form is read"
        " here, so keeping it is the conservative choice and is made for comparability, not for a"
        " multiplicity this lane has"
    ),
}
CONTROL = (
    "M2h's within-gene cell permutation, with the one change that nothing is refitted (there is nothing"
    " to refit): for each permutation, each gene's thirteen P values are permuted among its own thirteen"
    " cells by an independent draw of the seeded generator, every gene keeping its own prediction values"
    " and its own range and every pair keeping its own fresh D', and the Delta MSE is recomputed at the"
    " same constant SCALE. Only the tie between a prediction and the cell it was made for is broken."
    " PERMUTATIONS['count'] permutations at PERMUTATIONS['seed']; control p = (1 + permutations with"
    " Delta MSE >= observed) / (1 + count). It is reported in every reading and changes no verdict"
)
HARNESS = (
    "before any fresh value is opened: (1) every file of data/results/freshrna_sources.json is on disk"
    " with the sha256 the sources file pins, and the sources file's own accession lists reproduce the"
    " non-overlap (empty intersection of experiments, files and libraries with the READ set); (2) the"
    " genes table of P_SOURCE has 7h's committed 16,509 genes and 214,617 pairs, every gene with"
    " thirteen P and thirteen D and no symbol twice; (3) Delta MSE at scale 1 recomputed over all of 7h's"
    " pairs from that table equals 7h's committed delta_mse within TOLERANCE, which is the check that the"
    " P rows are the ones 7h's registered reading was computed on; (4) the constant SCALE equals"
    " SCALE_SOURCE's committed per_form.multiplicative.scale exactly; (5) M2h's committed held-out gene"
    " count and pair count are reproduced from the table and the split. THEN, with the fresh tables"
    " opened: (6) each of the thirteen tables has the row count and column names its sources entry"
    " records. A miss in any of them stops the run (HARNESS_DEFECT) and nothing is read"
)
TOLERANCE = mc.TOLERANCE

# ---- the decision, written before any fresh value was read ----------------------------------------
DECISION = (
    "1. If the HARNESS misses, the run stops before any fresh value is read (HARNESS_DEFECT). 2. Fewer"
    " than MIN_GENES genes in the primary population is NOT_TESTABLE and no comparison is read. 3."
    " Otherwise, on the primary population's Delta MSE at the constant SCALE with its 98.33% gene"
    " bootstrap interval: lower bound > 0 -> BETTER; upper bound < 0 -> WORSE (the falsifier); the"
    " interval straddles 0 -> INCONCLUSIVE. There is no fourth branch: M2h's CALIBRATION_REMOVES branch"
    " existed to tell a dead fitted scale from a live one, and nothing is fitted here, so a straddle is"
    " simply INCONCLUSIVE. 4. The control p is reported in every verdict sentence and changes no verdict."
    " 5. Every reading carries SCOPE's narrow statement of what it does and does not license, SIZE's"
    " 'about 2% of the variance' and NON_OVERLAP's two-cell biosample qualification. 6. A WORSE or an"
    " INCONCLUSIVE reading is reported FIRST and UNSOFTENED: it says that M2h's positive does not"
    " transport to a fresh measurement of the same cells, which is as valuable as a confirmation"
)
FALSIFIER = (
    "the reading under test is M2h's: that at the one magnitude +0.2146 the multiplicative composition of"
    " the compiled rules is closer to measured RNA than running no rule. It is FALSIFIED here when the"
    " Delta MSE interval on the primary population EXCLUDES 0 ON THE BASELINE'S SIDE (upper bound < 0,"
    " outcome WORSE): on a measurement of the same cells that this project had never read, the same form"
    " at the same magnitude is further from measured RNA than no rule at all, and M2h's positive is a"
    " property of the measurement it was read on. An interval straddling 0 (INCONCLUSIVE) falsifies"
    " nothing and confirms nothing, and is also reported as a failure to transport. An interval excluding"
    " 0 on the rules' side (BETTER) is the registered confirmation, in SCOPE's narrow terms"
)
OUTCOMES = {
    "NOT_TESTABLE": (
        "{form}: not testable at the registered size; {g} genes in the primary population, fewer than the"
        " registered {min}; no comparison was read"
    ),
    "BETTER": (
        "{form}: better; on a measurement of the same thirteen reader cell types that this project had"
        " never read, the reader-gated run scaled by M2h's fitted constant +0.2146 lowered the squared"
        " error of measured log2(TPM+1) across those cell types below the no-rule baseline's, so it is"
        " closer to measured RNA than the no-rule baseline on this measurement, and M2h's gain REPLICATES"
        " ON A FRESH MEASUREMENT OF THE SAME CELLS (Delta MSE {d} [{ci} {lo}, {hi}] over {n} gene-cell"
        " pairs, {g} genes; control p = {p}; the scale is the constant +0.2146 taken from ac5580b and"
        " NOTHING IS FITTED IN THIS LANE)"
    ),
    "WORSE": (
        "{form}: worse; on a measurement of the same thirteen reader cell types that this project had"
        " never read, the reader-gated run scaled by M2h's fitted constant +0.2146 RAISED the squared"
        " error of measured log2(TPM+1) across those cell types above the no-rule baseline's, so it is"
        " further from measured RNA than the no-rule baseline on this measurement and M2h's positive does"
        " NOT transport to a fresh measurement of the same cells (Delta MSE {d} [{ci} {lo}, {hi}] over"
        " {n} gene-cell pairs, {g} genes; control p = {p}; the scale is the constant +0.2146 taken from"
        " ac5580b and NOTHING IS FITTED IN THIS LANE)"
    ),
    "INCONCLUSIVE": (
        "{form}: inconclusive; on a measurement of the same thirteen reader cell types that this project"
        " had never read, the interval straddles 0, so this measurement does not resolve whether the"
        " reader-gated run scaled by M2h's fitted constant +0.2146 is closer to or further from measured"
        " log2(TPM+1) than the no-rule baseline, and M2h's gain is NOT shown to survive a fresh"
        " measurement of the same cells (Delta MSE {d} [{ci} {lo}, {hi}] over {n} gene-cell pairs, {g}"
        " genes; control p = {p}; the scale is the constant +0.2146 taken from ac5580b and NOTHING IS"
        " FITTED IN THIS LANE)"
    ),
    "HARNESS_DEFECT": (
        "stopped before any fresh value was read: the committed inputs did not reproduce their committed"
        " numbers ({what})"
    ),
}
STOP_RULE = (
    "ONE run, from this registration's committed code, on the files the sources file pins by sha256. No"
    " form, scale, population, split, endpoint, metric, interval, control, decision, falsifier, scope or"
    " wording is changed after a fresh value is read; a defect found afterwards is reported beside the"
    " registered reading and corrected only under a new registration. The post-hoc quantities of"
    " POST_HOC are computed only after the registered reading is written and committed, are labelled post"
    " hoc, and decide nothing"
)
POST_HOC = (
    "computed only AFTER the registered reading is written and committed, labelled post hoc, deciding"
    " nothing and changing no verdict: (1) THE DEPENDENCE-ROBUST INTERVAL the supervisor asked for -- the"
    " same Delta MSE with a 1-Mb LOCUS-CLUSTER bootstrap, genes assigned to clusters by (chromosome,"
    " gene start // 1,000,000) from data/cache/gencode_genes.tsv and CLUSTERS resampled with replacement,"
    " each carrying all of its genes and all of their pairs, at the same resample count and seed, because"
    " genes that share elements or a locus are not independent and the gene-level interval may be too"
    " narrow; the registered interval remains the gene-level one. (2) THE SCALE THE FRESH RNA WOULD IMPLY,"
    " sum(D'*P)/sum(P^2) over the primary population, beside the constant +0.2146 -- a DESCRIPTIVE"
    " comparison of two numbers, computed after the fact, which is NOT fitted into the registered reading"
    " and does not change it. (3) the gain as a fraction of the fresh mean(D'^2) (SIZE's 2% figure"
    " recomputed on this measurement) and the implied correlation. (4) the same reading restricted to the"
    " nine cells whose fresh experiment is polyA plus, with D' recentred over those nine, as a sensitivity"
    " to the mixed assay. (5) the per-cell Delta MSE. (6) the reading on M2h's training genes and on all"
    " of 7h's genes. (7) the correlation between the fresh y' and item 7's y per cell, as a descriptive"
    " check that the two measurements measure the same thing"
)
LIMITATIONS = (
    "NOT A FRESH CELL TYPE. This is a fresh measurement of the SAME thirteen cell types, because"
    " docs/ROADMAP.md §5's reinstated stop list says 'No more biosamples' and no reader biosample is"
    " added. Whether the gain transports to a cell type the reader has never held is untested here and"
    " cannot be tested here (SCOPE)",
    "NOT INDEPENDENT OF THE PREDICTOR SIDE. The predictions P are 7h's committed rows: the same compiled"
    " programs, the same attributed elements, the same fitted strengths and the same ENCODE DNase-seq"
    " reader calls, and the scale is M2h's. Only the measured side is new, so this test removes exposure"
    " to the earlier RNA and no other dependence (WHAT_IS_FRESH)",
    "A MIXED ASSAY. Nine cells' fresh measurement is polyA plus RNA-seq, as item 7's was; four (IMR-90,"
    " astrocyte, hepatocyte, cardiac muscle cell) are total RNA-seq, because the portal carries no unread"
    " qualifying polyA plus experiment of those terms. Total RNA-seq TPM includes non-polyadenylated"
    " transcripts, so for those four cells y' carries an assay difference item 7's y did not. Stated in"
    " advance: because the endpoint is 2*SCALE*mean(D'*P) - SCALE^2*mean(P^2) with SCALE and P fixed, a"
    " difference in D' that is unrelated to P leaves the expected Delta MSE unchanged and only widens the"
    " interval, and a difference that dilutes the biological signal can only shrink mean(D'*P). The mixed"
    " assay can therefore COST this test power or attenuate a real gain; it cannot manufacture one. A"
    " WORSE or INCONCLUSIVE reading must be read with that in mind, and a BETTER reading is not weakened"
    " by it. POST_HOC (4) reports the nine polyA plus cells on their own",
    "DIFFERENT DONORS FOR THE TWO TISSUES, AND ONE EMBRYONIC. Of the thirteen terms, ovary and testis are"
    " ENCODE tissues, not cell lines: item 7 read an adult ovary (47 years) and an adult testis (44"
    " years), and the unread experiments of those terms are a different adult ovary (30 years) and an"
    " EMBRYONIC testis. For those two cells the fresh measurement is a different donor, and for testis a"
    " different developmental stage, which is a larger difference than measurement noise and dilutes"
    " rather than inflates (same argument as the mixed assay)",
    "TWO CELLS SHARE A BIOSAMPLE. For hepatocyte and cardiac muscle cell the fresh experiment draws on"
    " the same ENCODE biosample as the read one, with different libraries; it is a fresh library and"
    " sequencing run of the same biological material (NON_OVERLAP)",
    "FEWER REPLICATES ON AVERAGE. Item 7's rule took, per cell, the experiment with the MOST qualifying"
    " files, so the unread experiments have on average fewer replicate files (22 files here against item"
    " 7's 27, and three cells are served by a single file). Fewer replicates means a noisier y', which"
    " again costs power rather than creating a gain",
    "THE SIZE. Even a confirmation is a gain of about 2% of the measured variance, r about 0.14 (SIZE). It"
    " is a notable result, not a breakthrough, and the wording never becomes 'the model predicts"
    " expression'",
    "THE P ROWS ARE THE COMMITTED TABLE'S FIVE-DECIMAL ROUNDINGS, which is why the harness admits"
    " TOLERANCE rather than exact equality",
    mr.LIMITATIONS[0],
)
#: the registration fields the run refuses to differ from: registration field -> constant of this module
BOUND = {
    "question": "QUESTION",
    "form": "FORM",
    "form_words": "FORM_WORDS",
    "scale": "SCALE",
    "scale_origin": "SCALE_ORIGIN",
    "nothing_is_fitted": "NOTHING_IS_FITTED",
    "scope": "SCOPE",
    "size": "SIZE",
    "p_source": "P_SOURCE",
    "scale_source": "SCALE_SOURCE",
    "rna_assays": "RNA_ASSAYS",
    "rna_selection": "RNA_SELECTION",
    "read_set_rule": "READ_SET_RULE",
    "non_overlap": "NON_OVERLAP",
    "what_is_fresh": "WHAT_IS_FRESH",
    "carried_verbatim": "CARRIED",
    "population": "POPULATION",
    "predictions": "PREDICTIONS",
    "measured_fresh": "MEASURED_FRESH",
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
    "decision": "DECISION",
    "falsifier": "FALSIFIER",
    "outcomes": "OUTCOMES",
    "stop_rule": "STOP_RULE",
    "post_hoc": "POST_HOC",
    "limitations": "LIMITATIONS",
}


def registered_fields() -> dict[str, Any]:
    """The BOUND fields, as the registration writes them and the run re-checks them."""
    out: dict[str, Any] = {}
    for field, name in BOUND.items():
        v = globals()[name]
        out[field] = list(v) if isinstance(v, tuple) else v
    return out


# ---- the split, only to name M2h's two sides -----------------------------------------------------
def is_held_out(symbol: str) -> bool:
    """M2h's committed split, carried unchanged: odd first 32 bits of sha256(symbol) are held out."""
    return int(hashlib.sha256(symbol.encode("utf-8")).hexdigest()[:8], 16) % 2 == 1


# ---- the metric, at the fixed scale --------------------------------------------------------------
def sums(rows: list[dict], key_p: str = "P", key_d: str = "D") -> tuple[float, float, float, int]:
    """(sum D*P, sum P^2, sum D^2, pairs) over every pair of the rows, P unscaled."""
    dp = pp = dd = 0.0
    n = 0
    for r in rows:
        for p, d in zip(r[key_p], r[key_d], strict=True):
            dp += p * d
            pp += p * p
            dd += d * d
            n += 1
    return dp, pp, dd, n


def delta_mse(rows: list[dict], scale: float = SCALE, key_p: str = "P", key_d: str = "D") -> float:
    """METRIC: mean of 2*scale*D*P - scale^2*P^2 over every pair, at the FIXED scale."""
    dp, pp, _dd, n = sums(rows, key_p, key_d)
    return (2.0 * scale * dp - scale * scale * pp) / n if n else float("nan")


def gene_terms(
    rows: list[dict], scale: float = SCALE, key_p: str = "P", key_d: str = "D"
) -> tuple[list[float], list[float]]:
    """Per gene, (its summed 2*scale*D*P - scale^2*P^2, its pair count): the bootstrap's unit."""
    s: list[float] = []
    n: list[float] = []
    for r in rows:
        tot = 0.0
        for p, d in zip(r[key_p], r[key_d], strict=True):
            tot += 2.0 * scale * p * d - scale * scale * p * p
        s.append(tot)
        n.append(float(len(r[key_d])))
    return s, n


def bootstrap_delta(
    rows: list[dict],
    resamples: int,
    seed: int,
    pct: tuple[float, float] = FORM_INTERVAL,
    scale: float = SCALE,
) -> tuple[float, float]:
    """BOOTSTRAP: percentile interval of Delta MSE over gene resamples, the scale held fixed."""
    import numpy as np

    s, n = gene_terms(rows, scale)
    sa, na = np.array(s), np.array(n)
    rng = np.random.default_rng(seed)
    stats = np.empty(resamples)
    for i in range(resamples):
        idx = rng.integers(0, len(rows), len(rows))
        stats[i] = sa[idx].sum() / na[idx].sum()
    lo, hi = np.percentile(stats, list(pct))
    return float(lo), float(hi)


def locus_clusters(rows: list[dict], span: int = 1_000_000) -> dict[tuple[str, int], list[int]]:
    """POST_HOC (1): (chromosome, start // span) -> the indices of the rows in that locus cluster."""
    out: dict[tuple[str, int], list[int]] = {}
    for i, r in enumerate(rows):
        key = (str(r.get("chrom")), int(r.get("start", 0)) // span)
        out.setdefault(key, []).append(i)
    return out


def bootstrap_delta_clustered(
    rows: list[dict],
    resamples: int,
    seed: int,
    pct: tuple[float, float] = FORM_INTERVAL,
    scale: float = SCALE,
    span: int = 1_000_000,
) -> tuple[float, float, int]:
    """POST_HOC (1): the same interval with 1-Mb locus CLUSTERS as the resampling unit."""
    import numpy as np

    s, n = gene_terms(rows, scale)
    sa, na = np.array(s), np.array(n)
    groups = list(locus_clusters(rows, span).values())
    cs = np.array([sa[g].sum() for g in groups])
    cn = np.array([na[g].sum() for g in groups])
    rng = np.random.default_rng(seed)
    stats = np.empty(resamples)
    for i in range(resamples):
        idx = rng.integers(0, len(groups), len(groups))
        stats[i] = cs[idx].sum() / cn[idx].sum()
    lo, hi = np.percentile(stats, list(pct))
    return float(lo), float(hi), len(groups)


def permuted_control(rows: list[dict], count: int, seed: int, scale: float = SCALE) -> list[float]:
    """CONTROL: each gene's P permuted among its own cells, Delta MSE recomputed at the FIXED scale."""
    import numpy as np

    rng = np.random.default_rng(seed)
    p = np.array([r["P"] for r in rows], dtype=float)
    d = np.array([r["D"] for r in rows], dtype=float)
    n = float(p.size)
    out: list[float] = []
    for _ in range(count):
        order = np.argsort(rng.random(p.shape), axis=1)
        pp = np.take_along_axis(p, order, axis=1)
        out.append(float((2 * scale * pp * d - scale * scale * pp * pp).sum() / n) if n else float("nan"))
    return out


def implied_scale(rows: list[dict]) -> float | None:
    """POST_HOC (2), DESCRIPTIVE AND POST HOC: sum(D*P)/sum(P^2) on the fresh rows.

    This is NOT a fit of this lane's endpoint: the registered reading is already written and committed
    at the constant SCALE before this is computed, and nothing here feeds back into it. It is called
    only by the post-hoc stage of the run, never by the registered one.
    """
    dp, pp, _dd, _n = sums(rows)
    return None if pp == 0.0 else dp / pp


def decide(lo: float, hi: float) -> str:
    """DECISION 3, the three-way table: no fourth branch, because nothing is fitted."""
    if lo > 0.0:
        return "BETTER"
    if hi < 0.0:
        return "WORSE"
    return "INCONCLUSIVE"


def reading(kind: str, **numbers: Any) -> str:
    """One of OUTCOMES, with the registered numbers in it."""
    return OUTCOMES[kind].format(**numbers)


def full_reading(sentence: str) -> str:
    """The verdict sentence with SCOPE, SIZE and NON_OVERLAP attached, as DECISION 5 requires."""
    return (
        f"{sentence}. WHAT THIS DOES AND DOES NOT LICENSE: {SCOPE}. SIZE: {SIZE}. NON-OVERLAP: {NON_OVERLAP}"
    )
