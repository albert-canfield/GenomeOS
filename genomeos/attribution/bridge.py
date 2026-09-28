# SPDX-License-Identifier: AGPL-3.0-or-later
"""Review R3 (2026-09-28): from executable annotation to a parameterised simulation, and no further.

A compiled program (`attribution/compile.py`) is **executable annotation**: it parses, its rules are
gated on their cell, `bio test` checks its counts, and a runtime will integrate it. It is not a
simulation of anything. Its rule `strength` is an observation's magnitude in that observation's own
unit (|log2 fold change| clipped at 1 for an AlphaGenome deletion, |fractional change| for a CRISPRi
screen), its target genes are stubs with no transcription parameters, and the element it names as a
source is not a species the runtime holds. `runtime/grn.py` read all three silently: a missing source
state as 0.0 and a missing `max_rate` as 0.0, so a compiled rule "ran" and moved nothing.

This module is the bridge, registered here before any of it was built. It says what one observation
may be turned into, what the model must be given before it can be turned into anything, and what is
reported when it cannot.

What it will NOT claim. A deletion log2 fold change is an observed effect of removing a piece of DNA
on one gene's steady-state expression in one cell or track. It is not a rate constant, not a binding
affinity, not a dose response and not a statement about any perturbation other than that removal.
The bridge fits one dimensionless rule strength so that the model reproduces that one observed
response, given transcription parameters and Hill assumptions the caller declares; everything else
the simulation then says is the model's, not the observation's.
"""

from __future__ import annotations

import copy
import re
from dataclasses import dataclass, field

from genomeos.ir import Action, Module, Rule
from genomeos.runtime.grn import Unresolved

REGISTERED = "2026-09-28"

# ---- the observation ------------------------------------------------------------------------------
#: what one compiled rule observed: removing element E changes gene G's steady-state expression in
#: cell C by the fold RHO = expression(E removed) / expression(E intact)
OBSERVABLE = "steady-state expression fold on removing the element: expression(removed) / expression(intact)"
#: how each evidence source's number becomes RHO, and the unit it arrives in
OBSERVATION_UNITS = {
    "predicted": (
        "log2 fold change on deletion (AlphaGenome), read from the evidence note 'effect <x> log2 fold"
        " change'; RHO = 2**x. The rule strength is |x| clipped at 1 and is used only when no note"
        " exists and it is below 1; a strength of exactly 1 is censored and unresolved"
    ),
    "experimental": (
        "fractional change in expression on CRISPRi silencing (ENCODE benchmark EffectSize), carried as"
        " the rule strength |f| with the sign from the action (activates: f < 0); RHO = 1 + f"
        ", used only when the link rests on one regulated pair (EXPERIMENTAL_SUMMARY)"
    ),
}
#: amendment of 2026-09-28, before the build (lane-assay's census A8 of measured.rule_links): a compiled
#: experimental rule's strength is the LARGEST |EffectSize| among the regulated pairs of one (element,
#: gene, cell) in the partition used, a maximum and not an observation whenever there is more than one
#: pair. The compiler writes each link's pair count on the measured element's evidence note
#: ("<gene> in <cell> from <n> pairs"); the bridge maps a link that rests on one pair and reports one
#: that rests on several as `observation_summarised`, never fitting to the maximum. No summary rule
#: (median, mean, meta-analysis) is chosen here; choosing one is its own registration.
EXPERIMENTAL_SUMMARY = "one pair maps; several pairs are unresolved (the compiled strength is their maximum)"

# ---- the model it maps into ----------------------------------------------------------------------
#: the element is a regulator whose state is its presence: dimensionless, 1 intact, 0 removed
SOURCE_STATE = "element presence, dimensionless: 1.0 intact, 0.0 removed; held by clamp, never defaulted"
INTACT, REMOVED = 1.0, 0.0
#: the target gene's transcription parameters, in the runtime's arbitrary units per hour
REQUIRED_GENE_PARAMETERS = {
    "basal_rate": "a.u./h, > 0: transcription with the element's regulation absent",
    "max_rate": "a.u./h, > 0: the regulated term's ceiling",
}
#: Hill assumptions, declared, not measured: threshold K and coefficient n of H(x) = x^n / (K^n + x^n)
HILL_THRESHOLD, HILL_COEFFICIENT = 1.0, 2.0
#: H at the intact state under those assumptions
H_INTACT = INTACT**HILL_COEFFICIENT / (HILL_THRESHOLD**HILL_COEFFICIENT + INTACT**HILL_COEFFICIENT)
#: the runtime's combination rules, named as model assumptions (lane-sign, b1f3405 and 5cbce26): a
#: gene's activators combine as a MEAN of s_i * H(x_i), so a second activator that is low where the
#: first is high HALVES the drive there; adding an activator can lower expression under this rule.
#: Inhibitors multiply. The bridge does not change either rule (one rule change per registration).
ACTIVATOR_COMBINATION = "mean"
INHIBITOR_COMBINATION = "product"
#: the mapping, one mechanism alone on its gene in its cell (b = basal_rate, V = max_rate, h = H_INTACT):
#:   activates (RHO < 1):  RHO = b / (b + V*s*h)            =>  s = b * (1/RHO - 1) / (V*h)
#:   inhibits  (RHO > 1):  RHO = (b + V) / (b + V*(1 - s*h)) =>  s = (b + V) * (1 - 1/RHO) / (V*h)
#: a fitted s outside (0, 1] means the declared parameters cannot produce the observed response
MAPPING = {
    "activates": "s = b * (1/RHO - 1) / (V * h)",
    "inhibits": "s = (b + V) * (1 - 1/RHO) / (V * h)",
}

# ---- measured rates: what a published number becomes (lane-rates, registered 2026-09-28) ----------
#: `basal_rate` and `max_rate` above are a.u./h that nobody has measured, so every compiled mechanism
#: that passes the other checks fails on `gene_parameter_missing`: 161,144 of them in
#: data/results/bridge_audit.json, beside 76,469 gene-cell pairs already unidentifiable. This block
#: says which published quantity becomes which model parameter, in which unit, for which species and
#: cell, and -- the part that matters -- what each one does NOT license. It was written and committed
#: before the rate table was built and before any count was read.
MEASURED_RATE_SOURCES = {
    "schwanhausser2011": (
        "Schwanhausser et al. 2011, Nature 473:337-342, doi 10.1038/nature10098, Supplementary Table 3"
        " (the file replaced 2013-02-13 for the corrigendum, Nature 495:126-127, doi 10.1038/nature11848),"
        " sha256 5343f73b4014a21c2b03630a579e95dabbf8b3c44f5ff0515db0df292abc12cd, 3,394,048 bytes."
        " MOUSE NIH 3T3 fibroblasts, one unperturbed steady state. 5,028 rows; 4,338 transcription"
        " rates, 4,658 mRNA half-lives, 5,028 protein half-lives, 4,309 mRNA copy numbers. Springer"
        " Nature supplementary information, no open licence stated: cached in the git-ignored"
        " data/cache/rates/, read here, never redistributed by this repository"
    ),
    "schofield2018": (
        "Schofield et al. 2018, Nature Methods 15:221-225, doi 10.1038/nmeth.4582, Supplementary Table 2,"
        " sha256 d5c010c1c2b1b22319310644e3194674ac58af059542862639cd43e720e60f19, 726,994 bytes."
        " TimeLapse-seq transcript half-lives in HUMAN K562 (5,419 transcripts, 4 h s4U) and MOUSE"
        " embryonic fibroblasts (2,992, 1 h s4U) by one method in one paper: the only matched"
        " cross-species pair here, and therefore the species falsifier below. Same licence position"
    ),
    "mgi_homology": (
        "Mouse Genome Informatics, HOM_MouseHumanSequence.rpt, retrieved 2026-09-28, 46,522 rows,"
        " sha256 b8220b7689f11066c3d82e0f24c23ba3c7a30ad26a0f24ab639357aebc09582d, 15,112,222 bytes."
        " Free for research use with attribution (MGI conditions of use). Used ONLY to carry a mouse"
        " symbol to a human symbol, and only within a homology class holding exactly one of each"
    ),
    "sender_milo_2021": (
        "already in the repository as bio.std.human_turnover (Sender & Milo 2021, Nat Med 27:45): HUMAN"
        " cell-type lifespans in DAYS. A cell's replacement clock. It is NOT an mRNA or protein"
        " half-life and never becomes delta_m or delta_p; it is listed here so that it is not reached for"
    ),
}
#: each measured quantity: the unit it arrives in, the parameter it becomes, and what it does not license
MEASURED_RATE_UNITS = {
    "transcription rate (vsr) [molecules/(cell*h)]": (
        "becomes T, the gene's TOTAL transcription in the unperturbed cell: T = basal_rate +"
        " max_rate * A * R at the observed state. It does NOT become max_rate. max_rate is the"
        " ceiling of the regulated term, reached only when every activator saturates, and no"
        " measurement in this source observes that state"
    ),
    "mRNA half-life [h]": (
        "becomes the degradation constant delta_m = ln2 / t_half [1/h] (runtime/grn.py's"
        " mrna_half_life). It is NOT a transcription rate and NOT a max_rate: it says only how fast"
        " mRNA disappears, and on its own it sets no level"
    ),
    "protein half-life [h]": (
        "becomes delta_p = ln2 / t_half [1/h] (Protein.half_life_h). It is NOT a translation rate"
    ),
    "translation rate constant (ksp) [molecules/(mRNA*h)]": (
        "becomes k_tl (runtime/grn.py's translation_rate). It is NOT a transcription rate"
    ),
    "mRNA copy number [molecules/cell]": (
        "becomes nothing. It is held as the steady state a run should land on, for the acceptance"
        " check only, and is never supplied as a parameter"
    ),
    "cell lifespan [d] (Sender & Milo)": "becomes nothing here; see MEASURED_RATE_SOURCES",
}
#: THE SPLIT. One measured T and one removal fold RHO determine the basal rate and the element's own
#: contribution exactly, with no invented fraction. For the only case that matters (a gene carrying
#: ONE mechanism in this cell -- a gene with more is already `not_identifiable`), at the intact state
#: h = H_INTACT:
#:   activates: intact = b + V*s*h = T and removed = b, and RHO = removed/intact
#:              =>  b = RHO * T           and  V*s*h = (1 - RHO) * T
#:   inhibits:  intact = b + V*(1 - s*h) = T and removed = b + V = RHO * T
#:              =>  V*s*h = (RHO - 1) * T  and  b = RHO * T - V
#: What is IDENTIFIED is C = |1 - RHO| * T, the element's own contribution to transcription in
#: molecules/(cell*h). V and s are NOT separately identified: only their product is. DECLARED_STRENGTH
#: fixes the split by convention and by nothing else.
DECLARED_STRENGTH = 1.0
#: what DECLARED_STRENGTH does NOT license: it is not a measurement that the element saturates its
#: gene. Every s in (0, 1] with V = C / (s * h) reproduces this observation identically; the
#: simulation's answer to any OTHER perturbation (a dose, a second regulator, a saturating input)
#: depends on V and s separately and is therefore not determined by any measurement used here.
DECLARED_STRENGTH_LICENSES = (
    "the absolute scale of one gene's transcription and its basal fraction, in molecules/(cell*h);"
    " not the ceiling, not the strength, and not the response to any perturbation but this removal"
)
#: a fold is dimensionless and the split is linear in T, so the FITTED STRENGTH does not depend on T at
#: all. A measured rate buys absolute units and the basal fraction; it buys nothing about the
#: regulation. This is why a borrowed rate is a separate, lower tier below and never a headline count.
RATE_TIERS = {
    "measured": "the gene's own measured T, carried to its human symbol through a 1:1 homology class",
    "borrowed_median": (
        "the genome median of the measured T, for a gene with no measurement. It fixes no absolute"
        " number: such a pair is simulable in RELATIVE units only, is counted separately and is never"
        " added to the measured count"
    ),
    "unmeasured": "no rate and no borrowing: `gene_rate_unmeasured`, an explicit marker, never a number",
}
#: EVERY rate here is mouse. A human gene simulated with them is borrowing a mouse fibroblast constant;
#: the result records `species: mouse` on every parameterised gene and the audit repeats it.
RATE_SPECIES = "Mus musculus, NIH 3T3 fibroblasts; every human simulation using them borrows a mouse constant"
#: ALREADY COMPUTED while surveying the source, so disclosed and not presented as a passed test:
#: Schwanhausser's vsr is not the identity N * ln2 / t_half over its own columns -- median relative
#: difference 0.236 across the 4,309 genes carrying all three, because vsr comes from their ODE fit to
#: time courses, not from the ratio. The acceptance check below therefore tests the runtime against
#: T / delta_m, which the model does determine, and reports the distance to the measured copy number
#: as an observation about the source.
#: the declared split puts the refitted strength exactly at DECLARED_STRENGTH, so `fit` returns 1 plus
#: a rounding error and the acceptance window (0, 1] needs a numerical tolerance at its upper end. It
#: is a floating-point allowance and nothing else: a strength above it is still `response_out_of_range`.
STRENGTH_TOLERANCE = 1e-9
RATE_SOURCE_INTERNAL_CONSISTENCY = (
    "median |vsr - N*ln2/t_half| / vsr = 0.236 over 4,309 genes (disclosed, not a test)"
)


def split_measured_rate(rho: float, total_rate: float, strength: float = DECLARED_STRENGTH) -> tuple:
    """(basal_rate, max_rate) in molecules/(cell*h) from one measured T and one removal fold RHO.

    The algebra above. A negative basal rate is a real outcome, not an error: an inhibitory
    mechanism whose removal raises the gene by more than 1 / (1 - strength * H_INTACT) leaves no
    room for a positive basal rate under the declared Hill assumptions, and the caller reports it as
    `rate_split_infeasible` rather than clipping it.
    """
    contribution = abs(1.0 - rho) * total_rate
    vmax = contribution / (strength * H_INTACT)
    basal = rho * total_rate if rho < 1.0 else rho * total_rate - vmax
    return basal, vmax


#: the acceptance test, fixed before the build, as lane-bridge's fixture was
RATE_ACCEPTANCE = (
    "end to end: a one-gene program whose gene carries a measured T and a measured mRNA half-life and"
    " whose one measured regulator is observed at RHO, parameterised through `parameterize(rates=...)`"
    " and integrated to steady state with the element clamped at 1 and at 0, reproduces RHO within 1%",
    "units: the same run's intact steady-state mRNA equals T / delta_m = T * t_half / ln2 within 1%,"
    " so the level the simulation reports is in molecules per cell and not in a.u.",
    "markers: a gene absent from the rate table is reported `gene_rate_unmeasured` and no number is"
    " put in its place; a gene whose split leaves basal <= 0 is reported `rate_split_infeasible`",
)
#: the falsifiers, none of them computed when this was committed
RATE_FALSIFIER = {
    "species_transfer": (
        "Schofield 2018 measured half-lives in mouse fibroblasts and human K562 by one method. Over MGI"
        " 1:1 homology classes the Spearman correlation of the two will be computed ONCE. If rho < 0.5"
        " the mouse constants are declared non-transferable: the rate table stays, every count keeps its"
        " `species: mouse` flag, and the registration states that the absolute human numbers are not"
        " supported by anything measured in a human cell"
    ),
    "inhibitor_bound": (
        "under the declared split an inhibitory mechanism needs RHO < 1 / (1 - DECLARED_STRENGTH *"
        " H_INTACT) = 2 to leave a positive basal rate. If more than half of the inhibitory mechanisms"
        " that reach the split fail it, the binding limit is the declared Hill assumption (K = 1,"
        " n = 2, so H = 0.5 at an intact element) and not the data, and that is the finding reported"
    ),
    "invariant": (
        "supplying rates must not change `not_identifiable` (76,469) or any other reason but"
        " `gene_parameter_missing`: rates are not an observation. A change there is a defect"
    ),
}
#: the expected count, with a direction, before the run
RATE_EXPECTED = (
    "of the 161,144 mechanisms now `gene_parameter_missing`, the number that becomes simulable on"
    " MEASURED rates will be strictly greater than 0 and strictly less than 161,144, and is predicted"
    " to fall between 15,000 and 60,000: the source carries 4,338 measured transcription rates against"
    " roughly 20,000 protein-coding genes, and its genes are the abundant ones, which may be over- or"
    " under-represented among the targets of compiled non-coding elements. Gene-cell pairs simulable"
    " rise from 0 by the same order. If the count lands outside 15,000-60,000 the miss is reported"
)

# ---- the second rate source: a survey for a HUMAN absolute rate (lane-rates2, 2026-09-28) ---------
#: The R3 rate table leaves 138,543 gene-cell pairs on a gene nobody measured, and every rate it does
#: carry is mouse NIH 3T3. Both limits have the same fix: a second source of ABSOLUTE transcription
#: rates, ideally human. This block is the survey, written before anything was applied. Each record is
#: what this lane read in the primary source itself, and -- the part that matters -- why the source
#: cannot supply T. A source that reports only half-lives cannot supply a synthesis rate; saying so is
#: the result, and no rate here is manufactured by pairing a half-life with a copy number from another
#: study. Reading a paywalled article was not bought: where that is the case the record says exactly
#: what was read (the deposit, the abstract) and claims nothing beyond it.
HUMAN_RATE_SOURCE_SURVEY = {
    "schwalb2016": {
        "citation": "Schwalb et al. 2016, Science 352:1225-1228, doi 10.1126/science.aad9841 (TT-seq)",
        "species": "Homo sapiens",
        "cell": "K562",
        "quantity": "RNA synthesis rates and half-lives from a 5 min 4sU pulse with RNA fragmentation",
        "units": "not established here: the article was not readable by this lane (science.org 403)",
        "absolute": "unknown from the article; NOT distributed as an absolute per-gene table",
        "n_genes": "over 10,000 transient RNAs mapped (from the GEO summary), rates per gene not deposited",
        "licence": "subscription article; GEO deposit GSE75792 is public",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75792",
        "why_not": (
            "what is actually distributed, and what this lane read, is the GEO deposit: a binned count"
            " table, a binning annotation, a transcript annotation and exon/intron count tables. Read"
            " counts, no per-gene rate and no cell-number or RNA-per-cell calibration. Deriving"
            " molecules/(cell*h) from it would be a new analysis of raw reads, not the use of a source"
        ),
    },
    "michel2017": {
        "citation": "Michel et al. 2017, Mol Syst Biol 13:920, PMC5371733 (TT-seq)",
        "species": "Homo sapiens",
        "cell": "Jurkat T cells",
        "quantity": "synthesis rate mu_i and degradation rate lambda_i per transcription unit",
        "units": "relative: labelled/total ratios under first-order kinetics",
        "absolute": "no",
        "n_genes": "22,141 transcription units (8,878 mRNA, 590 lincRNA, 12,673 ncRNA)",
        "licence": "CC BY 4.0",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5371733/",
        "why_not": (
            "the spike-ins set sequencing depth sigma_j and cross-contamination epsilon_j per sample,"
            " not a number of cells or an RNA mass per cell. mu_i is therefore a rate in library units"
            " and carries no molecules/cell. Nothing in it fixes the absolute scale T needs"
        ),
    },
    "wachutka2019": {
        "citation": "Wachutka et al. 2019, eLife 8:e45056, PMC6548502 (TT-seq)",
        "species": "Homo sapiens",
        "cell": "K562",
        "quantity": "synthesis and cleavage rates of individual phosphodiester bonds; bond half-lives",
        "units": "minutes for half-lives; synthesis spike-in normalised, no per-cell unit stated",
        "absolute": "no",
        "n_genes": "per splice site, not per gene",
        "licence": "CC BY",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6548502/",
        "why_not": "a bond-level rate is not a gene's transcription rate, and no per-cell unit is given",
    },
    "ietswaart2024": {
        "citation": "Ietswaart et al. 2024, Mol Cell 84:2765-2784, PMC11315470 (subcellular TimeLapse-seq)",
        "species": "Homo sapiens and Mus musculus",
        "cell": "K562 and NIH 3T3",
        "quantity": "chromatin-release, nuclear-export, polysome-loading and degradation rate constants",
        "units": "first-order rate constants (per hour) and half-lives, from the record read here",
        "absolute": "not established: the full text was not readable by this lane (cell.com 403, no OA copy)",
        "n_genes": "all expressed genes, per the abstract",
        "licence": "subscription article",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11315470/",
        "why_not": (
            "only the record and abstract were read, and they describe rate constants and half-lives,"
            " not a synthesis rate in molecules/(cell*h). It is named here as the FIRST PLACE THE NEXT"
            " LANE SHOULD LOOK, for a reason worth more than a rate: it measures human K562 and mouse"
            " NIH 3T3 -- the very cell type the R3 rates come from -- by one method, which would replace"
            " the present species falsifier (mouse MEF against human K562) with a matched one"
        ),
    },
    "shao2022": {
        "citation": "Shao et al. 2022, Mol Syst Biol 18:e10407, PMC8754154 (TT-seq, spike-in scaled)",
        "species": "Mus musculus",
        "cell": "embryonic stem cells (serum/LIF, 2i, mTORi)",
        "quantity": "RNA synthesis rate = labelled rate * transcript copy number per cell",
        "units": "cell^-1 min^-1 (copy/min per cell) -- genuinely absolute, in the source's own words",
        "absolute": "yes",
        "n_genes": "about 10,674 genes carry the TT-seq/Pol II quantities",
        "licence": "CC BY 4.0",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8754154/",
        "why_not": (
            "it is the one source found that reports an absolute synthesis rate per cell, and it is"
            " MOUSE, so it cannot answer the human question at all; it would only swap one borrowed"
            " species for the same one in another cell type. And its per-gene rates are not"
            " distributed: Table EV1 is elongation velocities, and the rates live in GSE168378 as raw"
            " and processed TT-seq. Using it means re-deriving rates from reads, a lane of its own"
        ),
    },
    "hausser2019": {
        "citation": "Hausser, Mayo, Keren & Alon 2019, Nat Commun 10:68, PMC6325141",
        "species": "Homo sapiens (also mouse, yeast, E. coli)",
        "cell": "HeLa",
        "quantity": "transcription rate beta_m for thousands of human genes",
        "units": "mRNA per cell per hour -- absolute in form",
        "absolute": "in form yes, in provenance no: CONSTRUCTED, not measured",
        "n_genes": "thousands (from Eichhorn et al. 2014 mRNA-seq and ribosome profiling)",
        "licence": "CC BY 4.0",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6325141/",
        "why_not": (
            "this is the closest thing to a human Schwanhausser, and reading it settles the question."
            " Their per-gene beta_m is an abundance rescaled by two constants: a total of about 225,000"
            " mRNAs per HeLa cell, which they obtained by taking the 180,000 mRNAs per MOUSE 3T3 cell"
            " and scaling it by a cell-volume ratio 2500/2000, and ONE global decay rate"
            " alpha_m = 0.06 /h from a median HEK293 half-life of 11.4 h. In their own words: 'We could"
            " not find direct measurements of the number of mRNAs per HeLa cell N_m.' So the published"
            " human transcription rates are anchored on the same mouse constant this lane is trying to"
            " escape, and per gene they are proportional to expression and carry no kinetics of their"
            " own. Adopting them would launder a borrowed mouse number into a human-labelled one"
        ),
    },
    "slam_drop_seq2023": {
        "citation": "Liu et al. 2023, Mol Syst Biol 19:e11427, PMC10568207 (SLAM-Drop-seq)",
        "species": "Homo sapiens",
        "cell": "HEK293 (A549 re-analysis; K562 used only as a half-life comparison)",
        "quantity": "transcription, splicing and degradation rates along the cell cycle",
        "units": "counts per hour after CPM normalisation ('from absolute counts to counts per million')",
        "absolute": "no: single-cell capture-limited counts, then rescaled to CPM",
        "n_genes": "a core set of 399 cycling genes",
        "licence": "CC BY 4.0",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC10568207/",
        "why_not": "CPM per hour is not molecules per cell per hour, and 399 genes would move nothing",
    },
    "nascseq2_2024": {
        "citation": "Ramskold et al. 2024, Nat Cell Biol 26:1725-1733, PMC11469958 (NASC-seq2)",
        "species": "Homo sapiens and Mus musculus",
        "cell": "K562 single cells; primary fibroblasts",
        "quantity": "transcriptional on/off rates, burst size and synthesis rate from new-RNA profiles",
        "units": "burst parameters; the per-cell counts reported are DETECTED molecules, capture-limited",
        "absolute": "no per-gene molecules/(cell*h) found in what was read",
        "n_genes": "single-cell, transcriptome-wide",
        "licence": "CC BY 4.0",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11469958/",
        "why_not": "detected molecules per cell depend on capture efficiency; burst rates are not T",
    },
    "schofield2018": {
        "citation": "Schofield et al. 2018, Nature Methods 15:221-225 (already cached; see above)",
        "species": "Homo sapiens and Mus musculus",
        "cell": "K562 and mouse embryonic fibroblasts",
        "quantity": "transcript half-lives",
        "units": "hours",
        "absolute": "n/a: a half-life is not a synthesis rate",
        "n_genes": "5,419 human transcripts, 2,992 mouse",
        "licence": "subscription supplementary, cached and never redistributed",
        "url": "https://www.nature.com/articles/nmeth.4582",
        "why_not": (
            "A HALF-LIFE CANNOT SUPPLY A SYNTHESIS RATE. It sets delta_m and nothing else. The only way"
            " it becomes a T is T = N * ln2 / t_half with a measured copy number N, and that is a"
            " construction, not a measurement; see `HUMAN_RATE_CONSTRUCTION_NOT_TAKEN`"
        ),
    },
}
#: the construction deliberately NOT made, and why it could not be made even if it were wanted
HUMAN_RATE_CONSTRUCTION_NOT_TAKEN = (
    "T = N * ln2 / t_half would turn a human half-life into a human absolute rate. It needs a measured"
    " per-gene human mRNA copy number N. This lane found none genome-wide, and the strongest evidence"
    " that none exists is that Hausser et al. 2019, who needed exactly that number, say they could not"
    " find even the per-cell TOTAL for HeLa and scaled the mouse 3T3 figure instead. There is also a"
    " measured reason to distrust the identity: over the 4,309 Schwanhausser genes carrying vsr, copy"
    " number and half-life together, median |vsr - N*ln2/t_half| / vsr = 0.236. So the construction has"
    " no human input AND a 24% disagreement with the one measurement it can be checked against. It is"
    " not made, and no number in the rate table comes from it"
)
#: the verdict, and the one sentence that is the lane's result
HUMAN_RATE_SURVEY_VERDICT = (
    "NO usable human absolute transcription-rate source exists. Human metabolic labelling (TT-seq,"
    " TimeLapse-seq, SLAM-seq, single-cell new-RNA) reports half-lives and rates in library or CPM"
    " units; the one study found that reports a genuinely absolute rate per cell (copy/min per cell) is"
    " MOUSE embryonic stem cells and does not distribute it per gene; and the one published set of"
    " human genome-wide transcription rates is constructed from abundance times a mouse-derived count"
    " of mRNAs per cell times one global decay constant. The R3 rates therefore stay mouse NIH 3T3, the"
    " 138,543 pairs on an unmeasured gene stay unmeasured, and every rate keeps its borrowed flag"
)
#: what would overturn the verdict, named in advance so the next lane has a target and not an opinion
HUMAN_RATE_SURVEY_FALSIFIER = (
    "a dataset distributing PER GENE, for a human cell type, a transcription or synthesis rate in"
    " molecules per cell per unit time (or a stated conversion to it), with the per-cell calibration"
    " named -- a cell count, an RNA mass per cell, or a spike-in weight per cell. Producing one"
    " overturns this record. Ietswaart et al. 2024 (human K562 and mouse NIH 3T3 by one method) is the"
    " first place to look, and Shao et al. 2022's GSE168378 shows the calibration that would be needed"
)
#: how a human rate would be kept apart from a mouse one if one ever arrives. Registered NOW, and
#: implemented now with a single source, so that a second source cannot be silently pooled with the
#: first: every row of the rate table carries its own `source` and `species_cell`, the audit counts
#: simulable pairs by the species of the rate they used, and a pair on a mouse rate is never reported
#: as a human result. Pooling two sources into one median, or averaging a human and a mouse rate for
#: one gene, is forbidden: with two sources the row keeps both and the audit reports them separately.
RATE_ROW_PROVENANCE = (
    "every rate row carries `source` (a key of MEASURED_RATE_SOURCES) and `species_cell`; the audit"
    " counts `pairs_on_a_human_measured_rate` and `pairs_on_a_borrowed_species_rate` per tier, and"
    " today the first is 0 by construction. Two sources are never averaged and never pooled"
)
#: the expected change, with a direction, before the re-run: NOTHING MOVES. A survey is not a rate.
HUMAN_RATE_SURVEY_EXPECTED = (
    "no count changes. measured stays 22,576 simulable gene-cell pairs, borrowed 160,392,"
    " gene_rate_unmeasured 138,543, not_identifiable 76,469, and the number of pairs on a human rate is"
    " 0. If any of those moves, the labelling changed a number it was not allowed to touch and the"
    " change is a defect, not a finding"
)


# ---- one mechanism, one parameter -----------------------------------------------------------------
#: a mechanism is (element, target gene, cell context); `<id>` and `<id>_measured` are one element
MECHANISM_KEY = ("element id without the _measured suffix", "target gene", "when clauses")
#: which representation a mechanism is simulated from when several exist; the rest are recorded as
#: superseded, never added. Another citation of the same kind with the same observation is the same
#: parameter; with a different observation it is a conflict, reported and not averaged
PRECEDENCE = ("experimental", "curated", "predicted", "inferred", "none")
CITATION_TOLERANCE = 1e-6  # |RHO_a - RHO_b| within this is the same observation

# ---- the diagnostic -------------------------------------------------------------------------------
#: every reason a compiled mechanism does not become a parameter; each is reported by name
UNRESOLVED = {
    "regulator_state_missing": "the source is not a species, not clamped and not a declared zero",
    "gene_parameter_missing": "the target declares no max_rate or no basal_rate > 0",
    "gene_rate_unmeasured": "no measured transcription rate for the target in the rate table supplied",
    "rate_split_infeasible": "the measured rate cannot carry this response and a positive basal rate",
    "observation_null_effect": "the observed removal fold is exactly 1: the mechanism moves nothing",
    "observation_missing": "no effect with a unit could be read from the rule",
    "observation_censored": "the only number is a strength clipped at 1",
    "observation_summarised": "an experimental strength is the maximum of several pairs, not one reading",
    "not_identifiable": "more than one mechanism regulates the gene in this context",
    "conflicting_observations": "two citations of one mechanism report different responses",
    "response_out_of_range": "the fitted strength falls outside (0, 1] for the declared parameters",
}

# ---- the double-counting audit, asked before it was answered ---------------------------------------
AUDIT_QUESTION = (
    "In the 24 compiled programs, how many (element, gene, cell) mechanisms carry both a predicted rule"
    " and a measured rule that one run in that cell would integrate side by side, and how many"
    " (gene, cell) pairs have more than one active regulatory rule, so that no single-deletion"
    " observation identifies a strength under the mean rule?"
)

# ---- acceptance, fixed before the build -----------------------------------------------------------
ACCEPTANCE = (
    "end-to-end: a program compiled from one element whose deletion is observed at log2 fold change -1"
    " in K562, parameterised with declared basal_rate and max_rate, reproduces expression(removed) /"
    " expression(intact) = 0.5 within 1% when the runtime is run with the element clamped at 1 and at 0",
    "diagnostic: running a program whose regulator has no state raises UnresolvedModel naming"
    " regulator_state_missing in strict mode and records it on the trajectory otherwise; a source"
    " renamed '<id>@zero' (network_experiment's edge knockout) is a declared zero and raises nothing",
    "citations: adding a second citation of the same mechanism (a duplicate rule citing another source"
    " with the same observation, or a predicted rule beside a measured one) leaves the simulated"
    " strength unchanged; two representations are never summed or averaged",
)


# ==== the build (after the registration in 06b7557 and its amendment) ==============================

_EFFECT = re.compile(r"effect ([+-]?\d+(?:\.\d*)?(?:e[+-]?\d+)?) log2 fold change")
_LINK = re.compile(r"(\S+) in (\S+) from (\d+) pairs")
MEASURED_SUFFIX = "_measured"


@dataclass(frozen=True, slots=True)
class Observation:
    """One observed removal response: RHO = expression(removed) / expression(intact)."""

    rho: float
    kind: str  # the evidence kind it came from
    how: str  # how RHO was read


@dataclass(slots=True)
class Bridged:
    """A simulation built from annotation, and everything that did not make it in, by name."""

    module: Module  # one fitted rule per resolved mechanism; use it only when `unresolved` is empty
    context: dict[str, str]
    strengths: dict[tuple, float] = field(default_factory=dict)  # mechanism -> fitted s
    observations: dict[tuple, Observation] = field(default_factory=dict)
    required_state: dict[str, float] = field(default_factory=dict)  # clamp these, intact = 1
    rates_used: dict[str, tuple] = field(default_factory=dict)  # gene -> (T, basal, max_rate), tier
    tier: dict[str, str] = field(default_factory=dict)  # gene -> RATE_TIERS key
    superseded: list[Rule] = field(default_factory=list)  # lower-precedence citations, never added
    unresolved: list[Unresolved] = field(default_factory=list)


def element_of(source: str) -> str:
    """`<id>` and `<id>_measured` are one element."""
    return source[: -len(MEASURED_SUFFIX)] if source.endswith(MEASURED_SUFFIX) else source


def mechanism_key(rule: Rule) -> tuple:
    return (element_of(rule.source), rule.target, tuple(sorted(rule.when.items())))


def _kind(rule: Rule) -> str:
    return str(getattr(rule.evidence.kind, "value", rule.evidence.kind))


def observation_of(rule: Rule, module: Module) -> Observation | Unresolved:
    """The removal response a compiled rule states, or why it states none (OBSERVATION_UNITS)."""
    kind = _kind(rule)
    if kind == "predicted":
        m = _EFFECT.search(rule.evidence.note)
        if m:
            return Observation(2.0 ** float(m.group(1)), kind, "evidence note, log2 fold change")
        if 0.0 < rule.strength < 1.0:
            x = -rule.strength if rule.action is Action.ACTIVATE else rule.strength
            return Observation(2.0**x, kind, "rule strength |log2 fold change|, below the clip")
        return Unresolved("observation_censored", rule.source, f"strength {rule.strength} and no effect note")
    if kind == "experimental":
        el = module.entities.get(rule.source)
        cell = rule.when.get("cell_type", "")
        n = None
        for g, c, k in _LINK.findall(el.evidence.note if el is not None else ""):
            if g == rule.target and c == cell:
                n = int(k)
        if n != 1:
            said = "no pair count stated" if n is None else f"the largest of {n} pairs"
            return Unresolved("observation_summarised", rule.source, f"{rule.target} in {cell}: {said}")
        f = -rule.strength if rule.action is Action.ACTIVATE else rule.strength
        return Observation(1.0 + f, kind, "rule strength |EffectSize| of one CRISPRi pair")
    return Unresolved("observation_missing", rule.source, f"{kind} evidence states no removal response")


def fit(rho: float, basal: float, max_rate: float) -> float:
    """The strength that makes one mechanism alone reproduce RHO (MAPPING); may fall outside (0, 1]."""
    if rho < 1.0:
        return basal * (1.0 / rho - 1.0) / (max_rate * H_INTACT)
    return (basal + max_rate) * (1.0 - 1.0 / rho) / (max_rate * H_INTACT)


def parameterize(
    module: Module,
    context: dict[str, str],
    genes: dict[str, dict[str, float]] | None = None,
    build: bool = True,
    rates: dict[str, float] | None = None,
    borrowed_rate: float | None = None,
) -> Bridged:
    """One fitted rule per resolved mechanism active in `context`; every other one reported by name.

    `build=False` classifies without copying the module or writing any rule (the audit's mode).
    `rates` maps a gene symbol to its MEASURED total transcription rate T in molecules/(cell*h);
    `split_measured_rate` turns T and the observed fold into basal_rate and max_rate. `borrowed_rate`
    is the single number a gene with no measurement of its own is given instead, which fixes no
    absolute scale (RATE_TIERS): passing it puts that gene in the `borrowed_median` tier, and leaving
    it None reports the gene as `gene_rate_unmeasured`. An explicit `genes` entry still wins over both.
    """
    genes = genes or {}
    out = copy.deepcopy(module) if build else module
    regulatory = (Action.ACTIVATE, Action.INHIBIT)
    active = [r for r in module.active_rules(context) if r.action in regulatory]
    if build:
        out.rules = [r for r in out.rules if r.action not in regulatory]
    b = Bridged(module=out, context=dict(context))

    by_key: dict[tuple, list[Rule]] = {}
    for r in active:
        by_key.setdefault(mechanism_key(r), []).append(r)
    by_gene: dict[str, list[tuple]] = {}
    for key in by_key:
        by_gene.setdefault(key[1], []).append(key)

    for gene, keys in sorted(by_gene.items()):
        if len(keys) > 1:
            b.unresolved.append(
                Unresolved(
                    "not_identifiable",
                    gene,
                    f"{len(keys)} mechanisms ({', '.join(k[0] for k in keys)}) under the"
                    f" {ACTIVATOR_COMBINATION}/{INHIBITOR_COMBINATION} rules",
                )
            )
            continue
        (key,) = keys
        reps = sorted(by_key[key], key=lambda r: PRECEDENCE.index(_kind(r)) if _kind(r) in PRECEDENCE else 99)
        top = _kind(reps[0])
        chosen = [r for r in reps if _kind(r) == top]
        b.superseded += [r for r in reps if _kind(r) != top]
        obs = [observation_of(r, module) for r in chosen]
        missing = [o for o in obs if isinstance(o, Unresolved)]
        if missing:
            b.unresolved.append(missing[0])
            continue
        rhos = [o.rho for o in obs if isinstance(o, Observation)]
        if max(rhos) - min(rhos) > CITATION_TOLERANCE:
            b.unresolved.append(
                Unresolved(
                    "conflicting_observations",
                    key[0],
                    f"{gene}: {', '.join(f'{x:.4g}' for x in rhos)}",
                )
            )
            continue
        g = out.entities.get(gene)
        rho = obs[0].rho if isinstance(obs[0], Observation) else 1.0
        basal = genes.get(gene, {}).get("basal_rate", getattr(g, "basal_rate", 0.0))
        vmax = genes.get(gene, {}).get("max_rate", (g.attrs.get("max_rate") if g is not None else None))
        tier = None
        if g is not None and not (basal and basal > 0 and vmax and vmax > 0) and rates is not None:
            total = rates.get(gene)
            tier = "measured" if total is not None else ("borrowed_median" if borrowed_rate else None)
            if tier is None:
                b.unresolved.append(
                    Unresolved(
                        "gene_rate_unmeasured", gene, f"no measured transcription rate ({RATE_SPECIES})"
                    )
                )
                continue
            total = float(total if total is not None else borrowed_rate)
            if rho == 1.0:
                b.unresolved.append(
                    Unresolved("observation_null_effect", key[0], f"{gene}: RHO is exactly 1")
                )
                continue
            basal, vmax = split_measured_rate(rho, total)
            if basal <= 0 or vmax <= 0:
                b.unresolved.append(
                    Unresolved(
                        "rate_split_infeasible",
                        key[0],
                        f"{gene}: RHO {rho:.4g} on T = {total:.4g} molecules/(cell*h) leaves"
                        f" basal {basal:.4g} at strength {DECLARED_STRENGTH}",
                    )
                )
                continue
            b.rates_used[gene] = (total, basal, vmax)
            b.tier[gene] = tier
        if g is None or not basal or basal <= 0 or not vmax or vmax <= 0:
            b.unresolved.append(
                Unresolved("gene_parameter_missing", gene, "needs basal_rate > 0 and max_rate > 0 (a.u./h)")
            )
            continue
        ob = obs[0]
        assert isinstance(ob, Observation)
        s = fit(ob.rho, basal, vmax) if ob.rho != 1.0 else 0.0
        if 1.0 < s <= 1.0 + STRENGTH_TOLERANCE:
            s = 1.0
        if not 0.0 < s <= 1.0:
            b.unresolved.append(
                Unresolved("response_out_of_range", key[0], f"{gene}: RHO {ob.rho:.4g} needs s = {s:.4g}")
            )
            continue
        b.strengths[key] = s
        b.observations[key] = ob
        b.required_state[key[0]] = INTACT
        if not build:
            continue
        g.basal_rate = float(basal)
        g.attrs["max_rate"] = float(vmax)
        src = chosen[0]
        out.rules.append(
            Rule(
                id=f"{key[0]}__{gene}",
                source=key[0],
                action=Action.ACTIVATE if ob.rho < 1.0 else Action.INHIBIT,
                target=gene,
                strength=s,
                threshold=HILL_THRESHOLD,
                hill=HILL_COEFFICIENT,
                when=dict(src.when),
                evidence=src.evidence,
            )
        )
    return b
