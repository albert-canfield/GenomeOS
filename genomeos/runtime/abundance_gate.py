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
