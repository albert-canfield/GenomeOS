# SPDX-License-Identifier: Apache-2.0
# Part of the BioLang engine (language, IR, VM, standard library); see LICENSING.md.
"""Pre-registration for stage 2's gate (b), written before any dataset is fetched.

§5.3 states it: *"Allocation must improve the prediction of absolute protein copy numbers from
transcript abundance, against PaxDb (Wang et al. 2015, Proteomics 15:3163) or another open absolute
set, relative to the current one-to-one mapping."* It also states the trap in advance, and that
sentence is why this gate is not a rank test:

    a shared pool with proportional allocation rescales every gene by the same factor, so it cannot
    change a rank correlation.

So a rank improvement credited to the pool layer would be credited to something else. The gate reads
**absolute error in log copy numbers** and **total protein per cell** (Milo 2013, BioEssays 35:1050),
and a rank improvement is creditable only to terms that differ per gene — length cost, half-life,
translation efficiency.

**The eligibility, established before any number exists.** The gate needs two arms and this project
holds neither:

- *the protein side.* PaxDb is named in the spec. It is reachable and versioned — v6.1, with a
  `paxdb-uniprot-crossreferences.txt` that maps its ids onto UniProt accessions, which is how this
  project's packaged proteome is keyed. Confirmed 2026-09-21; the exact per-species file was not
  resolved and is a fetch detail, not a design one.
- *the transcript side.* There is no proteome-wide transcript abundance set here. What exists is 102
  per-gene GTEx files and a 2.8 KB summary of TCGA BRCA z-scores against normal, which is a
  difference rather than an abundance. §5.3 does not name a transcript source, so one has to be
  chosen — and **it must be declared before it is used**, because choosing it after seeing which one
  makes allocation look better is the garden of forking paths with two doors.

**The coverage trap, registered because this project has been caught by it four times.** PaxDb does
not cover the proteome evenly: a protein is in it because somebody could quantify it, and abundant
proteins are easier to quantify. So the gate is a statement about the covered subset, the denominator
is printed beside every figure, and the covered and uncovered sets have their covariates compared —
length, and whatever the transcript source offers — so a reader can see which population the number
is about. `compare.input_presence` says which arm is *bought*: the protein measurement is, the
sequence length is not.

**What failure means, and §5.3 already agreed it.** *"If (b) fails, the pool layer stays optional and
this document says so."* So this gate can genuinely fail, the consequence is written down in advance,
and passing it is not the only acceptable outcome.
"""

from __future__ import annotations

from typing import Any

#: the gate reads error, not rank. Rank is still reported, because reporting it and refusing to use it
#: is what shows the trap was avoided rather than unnoticed.
PRIMARY = "mean absolute error in log10 protein copies per cell"
SECONDARY = "total protein per cell against Milo 2013 (order of magnitude)"

#: the improvement in log10 MAE that counts as the pool layer helping. Fixed before any data: 0.05 of
#: a log10 decade is about 12%, small enough to be reachable and large enough not to be noise on a
#: set of thousands, and the interval has to exclude zero as well.
MIN_IMPROVEMENT_LOG10 = 0.05

#: genes below this are dropped from both arms before anything is computed, because a copy number
#: under a handful per cell is at the quantification floor and its log is mostly noise.
MIN_COPIES = 10.0

#: DECLARED 2026-09-21, before either file was fetched and before any figure existed, because the
#: registration above forbids choosing a transcript source after seeing which one flatters the pool
#: layer. Two candidates were probed; only one answered.
TRANSCRIPT_SOURCE: dict[str, Any] = {
    "chosen": "GTEx v8 gene median TPM",
    "url": (
        "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/"
        "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
    ),
    "size_mb": 7.0,
    "why": (
        "it answered (HTTP 200 at 7 MB) where the Human Protein Atlas download paths did not (404 on"
        " both rna_celline and normal_tissue, so they have moved); it is the RNA source this project"
        " already streams and cites for eQTLs, so one citation covers both uses; and 7 MB distils"
        " without touching a disk that is at 91%"
    ),
    "tissue": "Cells - Cultured fibroblasts",
    "why_that_tissue": (
        "one named population rather than a cross-tissue average, fixed before the fetch. PaxDb's"
        " human data is weighted toward cultured cells, so a cultured-cell transcript arm is the"
        " closest single population on offer. The whole-body median is carried as a declared"
        " sensitivity, not as an alternative to be chosen afterwards"
    ),
    "population_mismatch": (
        "STATED, not resolved: GTEx is post-mortem bulk tissue and PaxDb is an average over"
        " experiments, so even the closest pairing is two populations. This is why the gate asks"
        " whether a shared-capacity model beats a proportionality constant at predicting a"
        " steady-state distribution, and not whether either is right about one cell"
    ),
    "rejected": {
        "Human Protein Atlas": (
            "download paths 404 as of 2026-09-21; it would have paired RNA and protein in one resource"
        ),
        "a cross-tissue GTEx average": "an average over populations is a population nobody sampled",
    },
}

PRE_REGISTRATION: dict[str, Any] = {
    "written": (
        "2026-09-21, before either dataset was fetched and before any figure existed; committed in "
        "this file before the script that runs it"
    ),
    "claim": (
        "the pool layer improves the prediction of absolute protein copy numbers from transcript "
        "abundance, against the one-to-one mapping it replaces"
    ),
    "primary": PRIMARY,
    "secondary": SECONDARY,
    "baseline": (
        "the current one-to-one mapping: protein copies proportional to transcript abundance with a "
        "single constant, fitted on the same genes the allocation arm uses so neither arm gets a "
        "denominator the other does not"
    ),
    "outcomes": {
        "pass": (
            f"log10 MAE falls by at least {MIN_IMPROVEMENT_LOG10} against the one-to-one baseline, "
            "with a bootstrap interval excluding zero, AND total protein per cell stays within an "
            "order of magnitude of Milo 2013"
        ),
        "fail": (
            "no improvement, or an improvement that disappears once the covered subset's covariates "
            "are matched. §5.3's consequence is already written: the pool layer stays optional and "
            "the document says so"
        ),
        "not_creditable": (
            "any improvement in RANK. Proportional allocation rescales every gene equally and cannot "
            "move a rank correlation, so a rank gain is evidence about some other term that differs "
            "per gene, not about allocation. It is reported and excluded from the verdict"
        ),
    },
    "denominators": (
        "genes with a value in both arms and at least "
        f"{MIN_COPIES} copies per cell. The count, the excluded count and the reason are printed "
        "beside every figure, and the covered and uncovered sets have their covariates compared"
    ),
    "cannot_do": (
        "PaxDb is an average over experiments and cell types, and the transcript source will be its "
        "own population; the two are not the same cell. So this tests whether a shared-capacity "
        "model predicts a steady-state abundance distribution better than a proportionality "
        "constant, not whether it is right about any one cell. And it cannot separate allocation "
        "from the cost terms that come with it, because a program states both at once"
    ),
}

# ---------------------------------------------------------------------------------------------
# Amendment 1
# ---------------------------------------------------------------------------------------------

#: the mammalian cell's total protein count, used to turn PaxDb's ppm into copies per cell. Milo 2013
#: (BioEssays 35:1050) puts it at 2-4 million proteins per cubic micron; a cultured human fibroblast
#: is of order 2,000-3,000 cubic microns, which is 4e9 to 1.2e10. The midpoint is taken and the band
#: is carried beside it. THE PRIMARY IS INVARIANT TO THIS CHOICE — in log10 it is an additive
#: constant, and both arms fit a constant — so it is load-bearing only for the secondary, which is
#: exactly where a chosen number should sit rather than hide.
MEASURED_TOTAL_PROTEINS = 5.0e9
MILO_BAND = (4.0e9, 1.2e10)

#: declared before the run, for the secondary's model side. None is fitted: the ribosome capacity is
#: §5.1's, and the two rates are cited measurements.
RIBOSOME_CAPACITY = 5.0e6  # §5.1, HeLa, BioNumbers
ELONGATION_AA_PER_S = 5.6  # mammalian ribosome elongation, Ingolia et al. 2011 / BioNumbers
PROTEIN_HALF_LIFE_H = 46.0  # median, Schwanhausser et al. 2011, Nature 473:337 (NIH3T3)

#: if the derivation in AMENDMENT below is right, the two fitted arms are the same arm and their MAE
#: difference is zero to numerical precision. A difference larger than this falsifies the derivation,
#: not the pool layer, and the amendment is withdrawn rather than defended.
DERIVATION_TOLERANCE_LOG10 = 0.001

AMENDMENT: dict[str, Any] = {
    "written": (
        "2026-09-21, before the transcript arm was fetched, before either arm was joined and before "
        "any figure existed. Reasoned from the runtime source — Economy._share in "
        "genomeos/runtime/economy.py — and not from data. What had been touched when this was "
        "written: PaxDb's header and its first eight rows, to confirm the file exists and carries a "
        "gene symbol and a ppm. No transcript abundance had been fetched, nothing had been joined, "
        "and no error, rank or slope had been computed"
    ),
    "finding": (
        "THE REGISTERED TRAP IS UNDERSTATED, and understated in a way that disqualifies the primary. "
        "§5.3 says a shared pool under proportional cannot move a RANK correlation. Reading the "
        "runtime shows the stronger statement: `_share` returns `capacity / wanted` to EVERY "
        "demander of an oversubscribed pool, one scalar with no index on the gene, and an entity's "
        "factor is the min over the pools it draws on. In a proteome-wide program every "
        "protein-coding gene draws the same two pools — a ribosome per chain and ATP per residue — "
        "so the min is over the same two scalars for everybody and the factor is ONE GLOBAL "
        "CONSTANT F. The allocation arm's prediction is therefore the baseline's prediction "
        "multiplied by a constant: log10 P_alloc = log10 P_base + log10 F"
    ),
    "consequence_for_the_primary": (
        "the baseline fits a constant, so the two arms span the same one-parameter family and the "
        "allocation arm cannot improve mean absolute error in log10 copies. Fitted against fitted "
        "they are the SAME ARM and the difference is exactly zero; unfitted, the allocation arm can "
        "only be worse, by however far its derived constant sits from the fitted optimum. So the "
        "registered pass condition — MAE falls by at least 0.05 with a bootstrap interval excluding "
        "zero — CANNOT BE MET by the pool layer as implemented, whatever the data say. The "
        "registered primary does not measure allocation. What it measures, once the constant is "
        "left underived, is whether the declared capacities put the ABSOLUTE SCALE in the right "
        "place, which is what the secondary already asks"
    ),
    "why_this_is_recorded_before_the_run": (
        "a pass condition that algebra has already closed is not a measurement, and reporting its "
        "'no improvement' afterwards as though the data had spoken would dress a derived identity "
        "up as an empirical negative. The distinction is the same one gate (a) drew when it "
        "registered its magnitude clause as not askable"
    ),
    "still_run_and_why": (
        "a derivation asserted is weaker than a derivation shown, so every registered reading is "
        "still computed on the real join: the two fitted arms are expected to agree to numerical "
        f"precision, and a gap above {DERIVATION_TOLERANCE_LOG10} log10 falsifies the derivation "
        "above rather than the pool layer, at which point this amendment is withdrawn in the open"
    ),
    "what_can_still_fail": {
        "secondary": (
            "the total protein per cell implied by the DECLARED ribosome capacity, nothing fitted: "
            f"{RIBOSOME_CAPACITY:.0e} ribosomes each elongating at {ELONGATION_AA_PER_S} aa/s over "
            "the expression-weighted mean chain length gives a synthesis rate, and a median protein "
            f"half-life of {PROTEIN_HALF_LIFE_H} h turns it into a steady-state count. It is "
            f"compared with Milo 2013's {MILO_BAND[0]:.0e}-{MILO_BAND[1]:.0e} and passes only "
            "within one order of magnitude of that band. It is not fitted, so it can be wrong by "
            "decades, and if it is, §5.3's consequence is the one already written down"
        ),
        "unfitted_allocation_arm": (
            "the allocation arm with its constant DERIVED from capacity rather than fitted, "
            "reported as the absolute-scale reading it actually is rather than as a rival to the "
            "baseline"
        ),
    },
    "instrument_that_would_make_the_primary_able_to_move": (
        "a term that differs PER GENE, and none of the three available ones is allocation. (1) "
        "`priority` is the only implemented policy whose share carries a gene index, and it needs an "
        "ordered list over ~19,000 genes that no measurement supplies, so it is not runnable at "
        "proteome scale. (2) A per-demander saturable share, C / (C + demand_i), instead of the "
        "per-pool C / (C + demand_total) that `competitive` computes TODAY — as implemented, "
        "`competitive` is as gene-blind as `proportional` and differs from it only in the value of "
        "the global constant, which is a finding about this runtime and is reported as one. A "
        "per-demander form would compress the dynamic range, which is the documented direction. (3) "
        "A per-gene degradation term, which is stage 4's, not stage 2's"
    ),
    "new_registered_quantity": (
        "the ordinary-least-squares slope of log10 protein copies on log10 transcript TPM over the "
        "covered set, with a bootstrap interval. It sizes the instrument above in decades per "
        "decade: a slope of 1 would say no per-gene compression is needed and the pool layer has "
        "nothing left to add even in principle; a slope below 1 states how much compression a "
        "future per-demander policy must supply to be worth writing. This number is reported "
        "whatever the verdict, and it is the only thing this gate can newly learn"
    ),
    "verdict_rule": (
        "the gate FAILS by its own registered rule, and the report must say that the failure is "
        "algebraic rather than empirical. §5.3's consequence stands as written: the pool layer stays "
        "optional and the document says so"
    ),
}
