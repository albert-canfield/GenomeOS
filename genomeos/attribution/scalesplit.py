# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""lane-scalesplit (ROADMAP §5, area I, 2026-10-12): does M2h's one scale survive a split that cannot
leak through shared compiled elements -- held out by CHROMOSOME?

M2h (lane-magcal; code 2b0b604, registration fbce775, result ac5580b) fitted ONE scale per form on
training genes and read it on held-out genes. Its multiplicative arm read "better": held-out Delta MSE
+0.04339 [98.33% CI +0.0379, +0.04894], fitted scale +0.2146, about 2% of the variance, control
p = 0.000999. M2i (lane-freshrna, f6e71c6) replicated it with nothing fitted, on RNA fresh to this
project, +0.03727 [+0.03202, +0.04299].

THE WEAKNESS THIS LANE TESTS. M2h's split is by sha256 of the gene SYMBOL -- a RANDOM split. Genes
that share compiled elements are not independent, which M2i's own post hoc acknowledged by reporting a
1-Mb locus-cluster interval beside its gene-level one. A random split can therefore LEAK: a training
gene and a held-out gene may be fitted through the same element, so "held out" is weaker than it
sounds. Chromosomes are disjoint in element membership -- the compiler emits one program per
chromosome and an element's targets are chosen inside that program, so no element is shared between
two chromosomes -- so a chromosome-held-out fit is a STRICTLY STRONGER generalisation test of the same
claim, on data already read.

What this lane does and does not do. It imports M2h's machinery read-only (magcal.fit_scale,
bootstrap_scale, bootstrap_delta, permuted_control, delta_mse, decide, scale_is_negative, reading,
combined_reading, quantiles and every registered word it carries) and edits no file of magcal's,
multrna's or combine2's. It reimplements nothing of them. It opens no compiled program, no narrowPeak
file and no RNA quantification: it reads 7h's committed genes table, the same rows M2h calibrated,
pinned by sha256. The scale is REFITTED on each fold's training chromosomes, so what is tested is the
scale's TRANSPORTABILITY and not one number's luck.

PRIOR EXPOSURE is stated plainly and is the largest limitation: this RNA was read by items 7, 7d, 7h,
the capped-sum census, M2h and M2i, so this is a RE-ANALYSIS of a measurement already read, not a
fresh test.

In the order written and committed: this code, the harness reproduction of M2h's committed numbers,
the registration, then one run. No model request, no key, no money, 0 AlphaGenome requests.
"""

from __future__ import annotations

from typing import Any

from genomeos.attribution import magcal as mc

REGISTERED = "2026-10-12"
LANE = "lane-scalesplit"

#: the one form read here: M2h's arm that read "better", and nothing else of M2h's design moves
FORM = "multiplicative"
FORM_WORDS = mc.FORM_WORDS[FORM]
SOURCE = "data/results/multrna_result.json"

# ---- what is carried verbatim from 7d, 7h and M2h ----------------------------------------------
CELLS = mc.CELLS
MEASURED = mc.MEASURED
NO_RULE_BASELINE = mc.NO_RULE_BASELINE
MIN_GENES = mc.MIN_GENES
CONTROL_ALPHA = mc.CONTROL_ALPHA
PERMUTATIONS = mc.PERMUTATIONS
BOOTSTRAP_RESAMPLES = mc.BOOTSTRAP["resamples"]
BOOTSTRAP_SEED = mc.BOOTSTRAP["seed"]
BOOTSTRAP_UNIT = mc.BOOTSTRAP["unit"]
#: magcal.FORM_INTERVAL, taken unchanged; see INTERVAL below for why it is not re-derived
FORM_INTERVAL = mc.FORM_INTERVAL
FORM_INTERVAL_WORDS = mc.FORM_INTERVAL_WORDS
TOLERANCE = mc.TOLERANCE

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
    "control": mc.CONTROL,
    "which_governs": mc.WHICH_GOVERNS,
    "not_vacuous": mc.NOT_VACUOUS,
    "outcomes": dict(mc.OUTCOMES),
    "scale_clauses": dict(mc.SCALE_CLAUSES),
    "bootstrap_unit": BOOTSTRAP_UNIT,
    "note": (
        "the endpoint, the no-rule baseline, Delta MSE, the closed-form fit, the gene bootstrap, the"
        " within-gene cell permutation control, MIN_GENES, the four-way decision table and EVERY reading"
        " word are 7d's, 7h's and M2h's, carried HERE WORD FOR WORD and imported from"
        " genomeos/attribution/magcal.py, which this lane does not edit. The population is 7h's, gene for"
        " gene: this lane defines no population of its own, admits no gene 7h did not admit and drops none"
        " it did. The ONLY thing this lane changes is the SPLIT: by chromosome instead of by sha256 of the"
        " gene symbol. Nothing else moves"
    ),
}

# ---- the registered question ---------------------------------------------------------------------
QUESTION = (
    "M2h fitted ONE scale on training genes and read it on held-out genes, splitting by sha256 of the gene"
    " SYMBOL -- a RANDOM split -- and its multiplicative arm read 'better' (held-out Delta MSE +0.04339"
    " [98.33% CI +0.0379, +0.04894], fitted scale +0.2146, about 2% of the variance, control p ="
    " 0.000999). Genes that share a compiled element are not independent, so a random split can leak: a"
    " training gene and a held-out gene may be fitted through the same element. Does the scale survive a"
    " split that CANNOT leak through shared elements -- held out by CHROMOSOME, which is disjoint in"
    " element membership -- with the scale REFITTED on the training chromosomes; and is the magnitude of"
    " the scale the same?"
)
WHY_STRONGER = (
    "the compiler emits ONE PROGRAM PER CHROMOSOME and an element's target is chosen inside that program,"
    " so an element on chr7 can enter the rules of chr7 genes only. The 24 chromosomes are therefore"
    " DISJOINT IN ELEMENT MEMBERSHIP, and in a chromosome-held-out fold no training gene and no held-out"
    " gene can share a compiled element. A chromosome-held-out fit is a strictly stronger generalisation"
    " test of M2h's claim than its random split, on the same rows. The converse is the thing to keep in"
    " view: a chromosome split ALSO breaks regional structure that a random split preserves (replication"
    " timing, gene density, compartment, chromosome-specific gene families), so a drop between the two"
    " splits is evidence that M2h's number is not free of element-level leakage WITHOUT being proof that"
    " leakage is the whole of the difference. That distinction is registered here and is not softened"
    " afterwards"
)
POPULATION = (
    "exactly the population of 7h's committed result, gene for gene and in the committed order: the"
    " multiplicative composition's 16,509 genes, each with its thirteen (gene, cell) pairs (214,617"
    " pairs). 7h's POPULATION is carried unchanged and so is M2h's choice of it; this lane applies no"
    " filter of its own. The rows are read from the committed genes table of"
    " data/results/multrna_result.json (pinned by sha256), which carries, per gene, its chromosome and the"
    " thirteen P and thirteen D the committed registered reading was computed from. Every one of the"
    " 16,509 genes is held out in exactly one of the two folds and is in the training side of exactly one"
)
PRIOR_EXPOSURE = (
    "STATED PLAINLY, AND NOT A CAVEAT ADDED AFTERWARDS. This is a RE-ANALYSIS of a measurement already"
    " read, not a fresh test. The same ENCODE polyA plus RNA-seq of the same thirteen reader biosamples,"
    " under the same log2(TPM+1) endpoint and the same no-rule baseline, has already been read by: item 7"
    " (M2b, lane-readercompile) including its post-hoc information term; item 7d (M2c, lane-combine2) for"
    " the runtime's mean and promoter competition; item 7h (M2f, lane-multrna) for the multiplicative"
    " composition; the capped-sum census (lane-capsumrna), which reached NOT_TESTABLE; M2h (lane-magcal),"
    " which fitted the one scale and read 'better' on a random split; and M2i (lane-freshrna), which"
    " replicated it with nothing fitted on RNA fresh to this project. The P and D re-split here ARE 7h's"
    " committed rows, which are the rows M2h calibrated, so EVERY value this lane reads has been read"
    " before -- and the chromosome split was chosen BECAUSE M2h's random split read 'better' and a random"
    " split can leak. The chromosome split limits that exposure and does not eliminate it: the partition"
    " is fixed by this registration before any held-out number is read and no held-out gene enters its"
    " fold's fit, but every gene's RNA value was in M2h's population, the endpoint and baseline were"
    " chosen with those readings in hand, and nothing here is a blind prediction of an unseen measurement."
    " Any reading from this lane therefore holds as a statement about THIS measurement's dependence on the"
    " split, and is not independent confirmation of M2h; it is the leakage audit of M2h and nothing more"
)
SPLIT = (
    "by CHROMOSOME, which is the one thing this lane changes. The 24 chromosomes carrying genes in 7h's"
    " committed table are partitioned by the PARITY OF THE CHROMOSOME NUMBER -- a property of the"
    " chromosome's NAME and of no measured value, fixed by this registration before any held-out number is"
    " read. ODD (eleven chromosomes): chr1, chr3, chr5, chr7, chr9, chr11, chr13, chr15, chr17, chr19,"
    " chr21. EVEN AND SEX (thirteen chromosomes): chr2, chr4, chr6, chr8, chr10, chr12, chr14, chr16,"
    " chr18, chr20, chr22, chrX, chrY. TWO COMPLEMENTARY FOLDS are registered and both are read, so that"
    " neither side is a chosen one: FOLD odd_held trains on EVEN AND SEX and holds out ODD; FOLD"
    " even_and_sex_held trains on ODD and holds out EVEN AND SEX. A gene is wholly on one side of a fold,"
    " with all thirteen of its pairs: no cell of a held-out gene enters that fold's fit. WHY PARITY: it is"
    " nameable in advance, needs no gene count, no RNA value and no fitted number to state, splits the"
    " genome into two pieces of comparable size, and puts no chromosome in both. The sex chromosomes are"
    " placed with the even autosomes by the one rule that they are not odd-numbered; chrY carries 13 genes"
    " and chrX 633, so that placement moves 646 of 16,509 genes and is declared here rather than tuned"
)
WHICH_SPLIT_IS_WHOSE = (
    "M2h's split is NOT re-run as this lane's endpoint and is NOT this lane's split: it is reproduced in"
    " the HARNESS only, to prove that this lane reads the same rows through the same functions and gets"
    " M2h's committed numbers back. This lane's registered endpoint is read on the chromosome folds alone"
)
REFITTED = (
    "THE SCALE IS REFITTED ON THE TRAINING CHROMOSOMES OF EACH FOLD. It is NOT taken as M2h's constant"
    " +0.2146. This is deliberate and is the point of the lane: refitting tests whether the SCALE ITSELF"
    " TRANSPORTS from one set of chromosomes to another, which is a property of the compiled rules, rather"
    " than testing whether one number fitted on a random half happens to hold -- which is M2h's question"
    " and is already answered. Each fold's scale is the closed-form minimiser over that fold's TRAINING"
    " pairs alone (magcal.FIT, imported unchanged), bootstrapped on that fold's training genes with the"
    " scale refitted on each resample (magcal.bootstrap_scale). No held-out pair enters any fit"
)
PREDICTIONS = {
    "form": mc.PREDICTIONS["form"],
    "calibrated": mc.PREDICTIONS["calibrated"],
    "no_rule_baseline": NO_RULE_BASELINE,
}
FIT = mc.FIT
METRIC = (
    mc.METRIC + ". Carried verbatim from M2h, with the HELD-OUT side now meaning the held-out CHROMOSOMES"
    " of the fold instead of the held-out gene symbols"
)
INTERVAL = (
    "98.33% percentile interval, magcal.FORM_INTERVAL taken UNCHANGED, DECIDED NOW AND NOT AFTER, for"
    " three reasons stated in advance. (1) COMPARABILITY: M2h's committed multiplicative reading is"
    " +0.04339 [98.33% CI +0.0379, +0.04894], and the whole point of this lane is to put its own number"
    " beside that one, which a different width would make a comparison of widths as well as of effects."
    " (2) MACHINERY: this lane calls magcal.bootstrap_delta and magcal.bootstrap_scale read-only and they"
    " read magcal.FORM_INTERVAL; changing the width would mean editing magcal's file or reimplementing its"
    " bootstrap, both of which this lane is forbidden to do and both of which would be worse than a"
    " conservative width. (3) CONSERVATISM: two folds of one form are judged here, so a Bonferroni"
    " correction from 95% over two readings would be 97.5% -- 7d's own two-arm convention. 98.33% is WIDER"
    " than that, so using it cannot manufacture an exclusion of 0 that the correct correction would deny."
    " M2h's 98.33% was Bonferroni over its three forms; here the same width is a deliberate"
    " over-correction and is labelled as one, never as a Bonferroni over two"
)
BOOTSTRAP = {
    "resamples": BOOTSTRAP_RESAMPLES,
    "seed": BOOTSTRAP_SEED,
    "unit": BOOTSTRAP_UNIT,
    "interval": (
        "per fold, the 98.33% percentile interval (INTERVAL). The held-out Delta MSE interval resamples the"
        " HELD-OUT CHROMOSOMES' genes with replacement, each carrying its thirteen pairs, with the scale"
        " held fixed at the value fitted on the training chromosomes (magcal.bootstrap_delta); the SCALE"
        " interval resamples the TRAINING CHROMOSOMES' genes with replacement, each carrying its thirteen"
        " pairs, and refits the scale on each resample (magcal.bootstrap_scale). The unit is the GENE in"
        " both, as M2h registered it. A chromosome-level bootstrap is NOT registered as the endpoint's"
        " interval and is not computed: with eleven and thirteen chromosomes it would be a resample of"
        " eleven or thirteen units, and M2h's gene unit is what is carried verbatim. That the gene unit"
        " does not price in between-chromosome correlation is written as a LIMITATION, not read as a result"
    ),
}
CONTROL = mc.CONTROL + (
    ". Applied per fold, with that fold's training chromosomes as the training side and its held-out"
    " chromosomes as the held-out side; magcal.permuted_control, imported unchanged"
)
HARNESS = (
    "before any chromosome split, scale, held-out value or control of this lane is computed, all of the"
    " following must hold. (1) the genes table of data/results/multrna_result.json has 7h's committed gene"
    " count (16,509) and pair count (214,617); (2) Delta MSE at scale 1 recomputed over all of its pairs"
    " equals 7h's committed delta_mse within TOLERANCE; (3) every gene's P and D have thirteen entries,"
    " every symbol is non-empty, no symbol occurs twice, and every row's chromosome is non-empty; (4)"
    " M2h'S COMMITTED GENE-LEVEL NUMBERS REPRODUCE: applying magcal.split_rows (sha256 of the symbol) and"
    " magcal's fit_scale, bootstrap_scale, delta_mse, bootstrap_delta and permuted_control at M2h's"
    " registered resamples, seeds and permutation count returns M2h's committed training gene count"
    " (8,181), held-out gene count (8,328), fitted scale, scale interval, held-out Delta MSE, held-out"
    " interval, control p and outcome (BETTER) -- counts and outcome exactly, every real number within"
    " M2H_TOLERANCE; (5) every chromosome in the table falls in exactly one of this lane's two chromosome"
    " sets, the two sets are disjoint, their union is the table's chromosomes, and the two folds' held-out"
    " gene counts sum to the table's gene count. A miss in ANY of them stops the run (HARNESS_DEFECT) and"
    " nothing is split by chromosome, fitted or read. (4) is the check that makes every later difference"
    " attributable to the SPLIT and to nothing else: same rows, same functions, same seeds, M2h's numbers"
    " back"
)
#: M2h's committed numbers are recomputed by the same code at the same seeds, so they must match closely
M2H_TOLERANCE = 1e-9
#: M2h's committed multiplicative numbers, quoted from data/results/magcal_result.json
M2H_COMMITTED = {
    "result": "data/results/magcal_result.json",
    "registration": "data/results/magcal_registration.json",
    "code": "genomeos/attribution/magcal.py",
    "commits": "2b0b604 code, fbce775 registration, ac5580b result",
    "training_genes": 8181,
    "held_out_genes": 8328,
    "scale": 0.2145911381525326,
    "scale_ci": [0.20074141299120243, 0.2283000359194418],
    "held_out_delta_mse": 0.0433912108101358,
    "held_out_ci": [0.0379045975793137, 0.04894398453032398],
    "control_p": 0.000999000999000999,
    "outcome": "BETTER",
    "variance_share": "about 2% of the variance (+0.04339 against mean(D^2) 2.1524 on the held-out side)",
}
#: 7h's committed numbers for the uncalibrated form, quoted from data/results/multrna_result.json
SEVEN_H_COMMITTED = {"genes": 16509, "pairs": 214617, "delta_mse": -0.5361286982098615}

# ---- the decision procedure, written before any held-out number of this lane was read ------------
#: M2h's four-way table, the supervisor's science ruling of 2026-10-10, carried verbatim and applied
#: as written, per fold, through magcal.decide
WHICH_GOVERNS = mc.WHICH_GOVERNS
NOT_VACUOUS = mc.NOT_VACUOUS
SCALE_CLAUSES = mc.SCALE_CLAUSES
OUTCOMES = mc.OUTCOMES
FOLDS = ("odd_held", "even_and_sex_held")
ODD_CHROMOSOMES = (
    "chr1",
    "chr3",
    "chr5",
    "chr7",
    "chr9",
    "chr11",
    "chr13",
    "chr15",
    "chr17",
    "chr19",
    "chr21",
)
EVEN_AND_SEX_CHROMOSOMES = (
    "chr2",
    "chr4",
    "chr6",
    "chr8",
    "chr10",
    "chr12",
    "chr14",
    "chr16",
    "chr18",
    "chr20",
    "chr22",
    "chrX",
    "chrY",
)
FOLD_WORDS = {
    "odd_held": (
        "the multiplicative composition (Zhou et al. 2024), scale fitted on the even-numbered autosomes"
        " and the sex chromosomes and read on the odd-numbered autosomes held out"
    ),
    "even_and_sex_held": (
        "the multiplicative composition (Zhou et al. 2024), scale fitted on the odd-numbered autosomes and"
        " read on the even-numbered autosomes and the sex chromosomes held out"
    ),
}
FOLD_TRAINING = {"odd_held": EVEN_AND_SEX_CHROMOSOMES, "even_and_sex_held": ODD_CHROMOSOMES}
FOLD_HELD_OUT = {"odd_held": ODD_CHROMOSOMES, "even_and_sex_held": EVEN_AND_SEX_CHROMOSOMES}

#: this lane's verdict ACROSS the two folds, the three cases named before anything was read
TRANSPORT = {
    "SURVIVES_BOTH_FOLDS": (
        "M2h's scale SURVIVES a chromosome-held-out split: on both complementary chromosome folds, with"
        " the scale refitted on the training chromosomes, the scaled multiplicative composition is closer"
        " to measured RNA than the no-rule baseline on chromosomes that entered no fit, so the gain is not"
        " an artefact of a random split leaking through shared compiled elements. It remains about 2% of"
        " the variance"
    ),
    "SURVIVES_ONE_FOLD": (
        "M2h's scale DOES NOT TRANSPORT CLEANLY ACROSS CHROMOSOMES: it reads better on one of the two"
        " complementary chromosome folds and not on the other, so the gain M2h measured on a random split"
        " of gene symbols is REGIONAL and is NOT shown to be free of leakage through shared compiled"
        " elements. This is a weakness in M2h's one positive result and is stated as one"
    ),
    "DOES_NOT_SURVIVE": (
        "M2h's scale DOES NOT SURVIVE a chromosome-held-out split: on NEITHER complementary chromosome"
        " fold, with the scale refitted on the training chromosomes, is the scaled multiplicative"
        " composition shown closer to measured RNA than the no-rule baseline on chromosomes that entered"
        " no fit. M2h's gain is therefore at least partly a LEAK through compiled elements shared between"
        " its training and held-out genes, which its random split did not exclude. This is the"
        " consequential outcome and it is not softened"
    ),
}
TRANSPORT_RULE = (
    "the transport verdict is a function of the two folds' outcomes and of nothing else, fixed here: both"
    " folds BETTER -> SURVIVES_BOTH_FOLDS; exactly one fold BETTER -> SURVIVES_ONE_FOLD; neither fold"
    " BETTER -> DOES_NOT_SURVIVE. A fold that is NOT_TESTABLE or NOT_FITTABLE counts as not BETTER and is"
    " named in the reading with its own registered sentence. The transport verdict NEVER overrides a"
    " fold's own sentence: both folds' registered sentences are reported in FOLDS order whatever the"
    " transport verdict is"
)
MAGNITUDE_COMPARISON = (
    "reported beside each fold's reading and DECIDING NOTHING: the fold's held-out Delta MSE and fitted"
    " scale against M2h's committed gene-level +0.04339 [98.33% CI +0.0379, +0.04894] and +0.2146"
    " [+0.2007, +0.2283]. THE RULE, FIXED NOW AND NOT AFTER: if two intervals OVERLAP, the honest reading"
    " is CONSISTENT WITH THE SAME EFFECT, and neither 'stronger', 'weaker', 'strengthened' nor 'weakened'"
    " may be said; only a fold interval lying WHOLLY BELOW M2h's lower bound may be read as the chromosome"
    " split giving a SMALLER number than the random split, and wholly above as a LARGER one. The two"
    " populations are NOT independent -- they are the same 16,509 genes repartitioned -- so no test of the"
    " difference is registered and none will be computed"
)
DECISION = (
    "1. If the HARNESS misses on any of its five points, the run stops before anything of this lane is"
    " split, fitted or read (HARNESS_DEFECT), and no chromosome-fold scale and no held-out value is"
    " computed. 2. Per fold, the scale is REFITTED on the fold's training chromosomes alone (REFITTED,"
    " FIT) and reported with its 98.33% training-gene bootstrap interval (refitted on each resample) in"
    " every reading, whatever the reading is; a training denominator of 0 is NOT_FITTABLE and the fold is"
    " read no further. 3. Per fold, fewer than MIN_GENES held-out genes is NOT_TESTABLE and no comparison"
    " is read. 4. Otherwise the fold's verdict is M2h's WHICH_GOVERNS four-way table, applied as written"
    " by magcal.decide, on the held-out Delta MSE at the refitted scale with its 98.33% held-out-gene"
    " bootstrap interval (the scale fixed): lower bound > 0 -> BETTER; upper bound < 0 -> WORSE (the"
    " falsifier); straddles 0 and the scale interval includes 0 -> CALIBRATION_REMOVES, a negative;"
    " straddles 0 and the scale interval excludes 0 -> INCONCLUSIVE. 5. A negative scale is reported"
    " BESIDE whichever of those applies, in SCALE_CLAUSES['negative']'s words, never instead of it. 6. The"
    " control p is reported in every fold sentence and changes no verdict. 7. The transport verdict is"
    " TRANSPORT_RULE applied to the two folds' outcomes. 8. The registered reading is: the transport"
    " sentence FIRST, then the two folds' registered sentences in FOLDS order joined by '; ', then"
    " MAGNITUDE_COMPARISON's finding, then PRIOR_EXPOSURE_CLAUSE, always"
)
FALSIFIER = (
    "the reading under test is M2h's, that one scale on the compiled rules' multiplicative composition"
    " brings the run closer to measured RNA on data that entered no fit, and that this is a property of"
    " the rules rather than of a random split. It is FALSIFIED for a fold when that fold's held-out Delta"
    " MSE interval excludes 0 on the BASELINE's side (upper bound < 0, outcome WORSE): with the scale"
    " refitted freely on the training chromosomes, the form is still further from measured RNA than"
    " running no rule at all on chromosomes that entered no fit. An interval STRADDLING 0 falsifies"
    " nothing and is read by the scale, exactly as M2h registered: CALIBRATION_REMOVES when the scale"
    " interval includes 0 -- 'calibration removes the rules', A NEGATIVE AND NOT 'no difference detected'"
    " -- and INCONCLUSIVE when the scale interval excludes 0, meaning a magnitude is estimated and this"
    " measurement does not resolve the comparison. 'INCONCLUSIVE' is therefore exactly one case here: a"
    " live scale with a straddling held-out interval, and it is the one case from which nothing is"
    " concluded in either direction. An interval excluding 0 on the rules' side (BETTER) is a positive for"
    " that fold. ACROSS the folds, the thing being falsified is TRANSPORT: anything other than BETTER on"
    " BOTH folds means M2h's gain is not shown to survive a split that cannot leak through shared"
    " elements, and SURVIVES_ONE_FOLD and DOES_NOT_SURVIVE both say so in their own registered words."
    " NOT_VACUOUS, carried from M2h, states why 'better' is the one reading a dead scale cannot fake"
)
PRIOR_EXPOSURE_CLAUSE = mc.PRIOR_EXPOSURE_CLAUSE + (
    " Added by this lane and binding on its own reading: the same rows were then read by M2h"
    " (lane-magcal), which fitted the one scale and read 'better' on a random split, and by M2i"
    " (lane-freshrna), which replicated it with nothing fitted on RNA fresh to this project. This lane is"
    " a RE-ANALYSIS of a measurement already read -- a leakage audit of M2h's split -- and not a fresh"
    " test of anything"
)
STOP_RULE = (
    "one run, from this registration's committed code, on the files the registration pins by sha256 (7h's"
    " and M2h's registrations and results). No chromosome set, fold, population, scale definition,"
    " endpoint, metric, interval, control, decision, transport rule, falsifier or wording is changed after"
    " a held-out value of this lane is read; a defect found afterwards is reported beside the registered"
    " reading and corrected only under a new registration. The post hoc of POST_HOC decides nothing and"
    " cannot change a verdict"
)
LIMITATIONS = (
    "PRIOR EXPOSURE, first because it is the largest: see PRIOR_EXPOSURE. These are the rows M2h"
    " calibrated and 7h, 7d and item 7 read before it. No number here is a blind prediction, and a fold"
    " that reads 'better' is not independent confirmation of M2h -- it is M2h's own measurement surviving"
    " a stronger split of itself",
    "the bootstrap unit is the GENE, carried verbatim from M2h, so the held-out interval prices in"
    " variation between genes and NOT between chromosomes. With eleven and thirteen chromosomes a"
    " chromosome-level resample was considered and deliberately NOT registered, because changing M2h's"
    " unit would change more than the split. A fold interval is therefore narrower than one that priced in"
    " between-chromosome correlation would be, and the post-hoc per-chromosome scales are reported so that"
    " a reader can see that spread directly",
    "a chromosome split removes element sharing AND regional structure at once. A drop from M2h's number"
    " is evidence that M2h's random split is not clean, but it does not apportion the drop between element"
    " leakage and genuine regional variation in the scale, and nothing here attempts to",
    "two folds of one form are read and the interval is 98.33%, wider than a Bonferroni over two would"
    " require. That is conservative for exclusion of 0 and it is NOT conservative for a straddle: an"
    " over-wide interval makes CALIBRATION_REMOVES and INCONCLUSIVE easier to reach than the correct"
    " correction would. Both directions of that are stated here in advance",
    "one scale is the smallest calibration there is and also the weakest: it can only shrink, grow or flip"
    " the predictions as a block, and cannot fix a form whose error is in the SHAPE of its prediction"
    " across cells or across genes. Carried from M2h",
    "the calibration is a statistical rescaling and not a mechanism; nothing here attributes the size to"
    " any stage of the compiled pipeline. Carried from M2h",
    "the P and D are the committed table's five-decimal roundings, which is why the harness admits"
    " TOLERANCE rather than exact equality. Carried from M2h",
    "chromosomes are disjoint in element membership BY CONSTRUCTION of the compiler (one program per"
    " chromosome, cis targets). That is a property of the compiled programs and not of biology: it says"
    " nothing about whether a real regulatory element can act between chromosomes, and if a future"
    " compiler emitted a trans rule WHY_STRONGER would no longer hold",
    "accessibility is not function, and an emitted target is a model attribution and never a measured"
    " regulatory relationship. Carried forward unchanged",
)
POST_HOC = (
    "reported beside the registered reading, DECIDING NOTHING and able to change no verdict: (1) THE"
    " PER-CHROMOSOME SCALES -- for each of the 24 chromosomes, the one scale fitted on THAT CHROMOSOME'S"
    " GENES ALONE with its gene and pair count, so a reader can see whether the magnitude is uniform or"
    " regional; (2) for each of the 24 chromosomes, the scale fitted on the OTHER 23 (a"
    " leave-one-chromosome-out training scale) beside its own, which is what a 24-fold version of this"
    " lane would have fitted; (3) the spread of (1): its minimum, maximum, median and interquartile range,"
    " and how many of the 24 lie inside M2h's committed scale interval. A chromosome with fewer than"
    " MIN_GENES genes is reported with its count and FLAGGED as below the registered floor -- chrY carries"
    " 13 genes in 7h's table -- and no chromosome-level scale is read as a verdict about that chromosome."
    " None of (1), (2) or (3) is a test: no interval is placed on a per-chromosome scale and no comparison"
    " between chromosomes is registered"
)
DESCRIPTIVES = (
    "per fold, deciding nothing: the held-out Delta MSE at scale 1 (the uncalibrated form on the held-out"
    " chromosomes) with its 98.33% interval; the training Delta MSE at the refitted scale and at 1; the"
    " training and held-out gene and pair counts; mean(D^2), mean(P^2), mean(scale^2 P^2) and mean(D*P) on"
    " each side, so the size the calibration keeps is visible as a number; the held-out Delta MSE per"
    " reader cell; the predicted range quantiles before and after scaling; the control's median and 95th"
    " percentile",
    "context, copied verbatim and not re-decided: 7h's registered reading, M2h's registered reading and"
    " its committed numbers, and M2i's replication numbers",
)
#: the registration fields the run refuses to differ from: registration field -> constant of this module
BOUND = {
    "question": "QUESTION",
    "why_stronger": "WHY_STRONGER",
    "form": "FORM",
    "form_words": "FORM_WORDS",
    "source": "SOURCE",
    "carried_verbatim": "CARRIED",
    "population": "POPULATION",
    "prior_exposure": "PRIOR_EXPOSURE",
    "split": "SPLIT",
    "which_split_is_whose": "WHICH_SPLIT_IS_WHOSE",
    "refitted": "REFITTED",
    "folds": "FOLDS",
    "fold_words": "FOLD_WORDS",
    "odd_chromosomes": "ODD_CHROMOSOMES",
    "even_and_sex_chromosomes": "EVEN_AND_SEX_CHROMOSOMES",
    "predictions": "PREDICTIONS",
    "fit": "FIT",
    "metric": "METRIC",
    "interval": "INTERVAL",
    "form_interval": "FORM_INTERVAL",
    "form_interval_words": "FORM_INTERVAL_WORDS",
    "bootstrap": "BOOTSTRAP",
    "control": "CONTROL",
    "permutations": "PERMUTATIONS",
    "control_alpha": "CONTROL_ALPHA",
    "min_genes": "MIN_GENES",
    "harness": "HARNESS",
    "tolerance": "TOLERANCE",
    "m2h_tolerance": "M2H_TOLERANCE",
    "m2h_committed": "M2H_COMMITTED",
    "seven_h_committed": "SEVEN_H_COMMITTED",
    "which_governs": "WHICH_GOVERNS",
    "not_vacuous": "NOT_VACUOUS",
    "scale_clauses": "SCALE_CLAUSES",
    "outcomes": "OUTCOMES",
    "transport": "TRANSPORT",
    "transport_rule": "TRANSPORT_RULE",
    "magnitude_comparison": "MAGNITUDE_COMPARISON",
    "decision": "DECISION",
    "falsifier": "FALSIFIER",
    "prior_exposure_clause": "PRIOR_EXPOSURE_CLAUSE",
    "stop_rule": "STOP_RULE",
    "limitations": "LIMITATIONS",
    "post_hoc": "POST_HOC",
    "descriptives": "DESCRIPTIVES",
}


def registered_fields() -> dict[str, Any]:
    """The BOUND fields, as the registration writes them and the run re-checks them."""
    out: dict[str, Any] = {}
    for field, name in BOUND.items():
        v = globals()[name]
        out[field] = list(v) if isinstance(v, tuple) else v
    return out


# ---- the split: by chromosome, SPLIT's parity rule -----------------------------------------------
def chromosome_set(chrom: str) -> str:
    """SPLIT: which of the two named chromosome sets a chromosome is in, by the parity of its name."""
    if chrom in ODD_CHROMOSOMES:
        return "odd"
    if chrom in EVEN_AND_SEX_CHROMOSOMES:
        return "even_and_sex"
    return "unnamed"


def split_rows_by_chromosome(rows: list[dict], fold: str) -> tuple[list[dict], list[dict]]:
    """(training, held out) for one registered fold, by SPLIT. A gene is wholly on one side."""
    train_set = set(FOLD_TRAINING[fold])
    held_set = set(FOLD_HELD_OUT[fold])
    train = [r for r in rows if r["chrom"] in train_set]
    held = [r for r in rows if r["chrom"] in held_set]
    return train, held


def transport_verdict(outcomes: dict[str, str]) -> str:
    """TRANSPORT_RULE, applied as written: a function of the two folds' outcomes and nothing else."""
    better = sum(1 for f in FOLDS if outcomes.get(f) == "BETTER")
    if better == len(FOLDS):
        return "SURVIVES_BOTH_FOLDS"
    return "SURVIVES_ONE_FOLD" if better else "DOES_NOT_SURVIVE"


def overlaps(a: tuple[float, float], b: tuple[float, float]) -> bool:
    """MAGNITUDE_COMPARISON: do two intervals overlap? Overlap means 'consistent with the same effect'."""
    return not (a[1] < b[0] or b[1] < a[0])


def magnitude_finding(lo: float, hi: float, m2h_lo: float, m2h_hi: float, what: str) -> str:
    """MAGNITUDE_COMPARISON's rule, applied as written and deciding nothing."""
    if overlaps((lo, hi), (m2h_lo, m2h_hi)):
        return (
            f"{what}: CONSISTENT WITH THE SAME EFFECT as M2h's -- the intervals overlap"
            f" ([{lo:+.4g}, {hi:+.4g}] against M2h's [{m2h_lo:+.4g}, {m2h_hi:+.4g}]), so neither"
            " 'stronger' nor 'weaker' is said"
        )
    if hi < m2h_lo:
        return (
            f"{what}: SMALLER than M2h's -- the interval [{lo:+.4g}, {hi:+.4g}] lies wholly below M2h's"
            f" lower bound {m2h_lo:+.4g}, on the same 16,509 genes repartitioned, so no test of the"
            " difference is offered"
        )
    return (
        f"{what}: LARGER than M2h's -- the interval [{lo:+.4g}, {hi:+.4g}] lies wholly above M2h's upper"
        f" bound {m2h_hi:+.4g}, on the same 16,509 genes repartitioned, so no test of the difference is"
        " offered"
    )


def registered_reading(transport: str, fold_readings: list[str], magnitude: list[str]) -> str:
    """DECISION step 8: the transport sentence first, then the folds, the magnitude, the exposure."""
    return (
        TRANSPORT[transport]
        + ". "
        + "; ".join(fold_readings)
        + ". "
        + " ".join(magnitude)
        + ". "
        + PRIOR_EXPOSURE_CLAUSE
    )
