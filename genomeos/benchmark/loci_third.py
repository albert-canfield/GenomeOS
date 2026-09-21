# SPDX-License-Identifier: AGPL-3.0-or-later
"""A third set of published enhancer-gene loci, chosen for what the first two could not test.

docs/LOCI-BENCHMARK.md now holds two sampling frames. The seventeen of `loci.PANEL` were
hand-curated as famous loci and score 15/17 on `target_derived`; the nine of
`loci_candidates.CANDIDATES` were assembled from perturbation experiments and score 8/9. Section
19's own "left undone" names three gaps, and this module is aimed at two of them:

  **the direction axis is only ever tested in one direction.** Every one of the nine activates, so
  "the deletion moves the target the published way" has never been asked of an element whose
  published job is to hold a gene *down*. Two of the nine got the sign wrong; nobody can say whether
  that is a bias of the model or of the question, because the question only came one way. Two of the
  nine loci here are published **repressors** with a perturbation behind them.

  **a third set, length-matched, is the next thing worth a day.** Section 20 found the panel's
  value-on-syntax reading is a length gradient. This set spans 300 bp to 52 kb on purpose, and the
  52 kb one is a human deletion, not a drawn interval.

The third gap - a published target that is **non-coding**, which would exercise the
`predicted_coding`-first defect H19 exposed - was looked for and **not filled**. It is registered
below as an aim that failed, with the reason, because that is a statement about the literature and
not about this code.

Nothing here is a new scorer. Every reading, hit rule, denominator and matched control comes from
`genomeos.benchmark.loci`, called with a different panel, exactly as `loci_candidates` does. The
reach filter is `loci.read_reach` - a gene-**body** test against the scorer's 1,048,576 bp input,
never a TSS-distance test, for the reason `loci_candidates.BODY_NOT_TSS_CASE` pins.

`PREREGISTRATION` was written before any locus was read through any layer and before any request was
spent. A negative or inconclusive result is the outcome, not a reason to re-cut.
"""

from __future__ import annotations

import contextlib
from pathlib import Path
from typing import Any

from genomeos.benchmark import loci, loci_candidates
from genomeos.benchmark.loci import Expect
from genomeos.results import RESULTS_DIR, load_result, save_result

#: the result this set's stated-interval deletions are written to, labelled `stated_interval` by
#: `loci._deletion_rows` so a hit on an interval this module drew is never pooled with a hit on one
#: ENCODE drew. `loci.STATED_INTERVAL_RESULTS` does not list it yet - see `register_stated_intervals`.
INTERVALS = "loci_third_intervals"
NAME = "loci_third"
#: this set's own distilled GTEx hits. `loci.stream_gtex` clears a directory whose window manifest
#: changed, so three panels sharing one directory would delete each other's cache in a checkout
#: several sessions share. One directory per panel, and none of them touches another.
GTEX_DIR = Path("data/knowledge/loci_third")

# ------------------------------------------------------------------ the registration, before any run
PREREGISTRATION: dict[str, Any] = {
    "written": "2026-09-21, before any locus here was read through any layer or any request was spent",
    "loci": 9,
    "frame": (
        "the THIRD sampling frame in docs/LOCI-BENCHMARK.md. Not merged into loci.PANEL and not"
        " merged into loci_candidates.CANDIDATES; its rate is reported beside 15/17 and 8/9 with each"
        " frame named and is never pooled with either. Where it disagrees with either, the"
        " disagreement is a statement about how the three sets were sampled"
    ),
    "how_they_were_chosen": (
        "published enhancer-gene or element-gene links established by PERTURBATION - CRISPR or mouse"
        " deletion of the element, a human deletion that leaves the gene intact, a transgenic"
        " reporter, or a disease-causing point mutation inside the element - never by association"
        " alone. No element already in loci.PANEL or loci_candidates.CANDIDATES. Chosen to carry the"
        " two properties section 19 said its own set could not: a published REPRESSIVE direction, and"
        " element lengths spanning the panel's range rather than a uniform 500 bp"
    ),
    "deliberate_target_overlaps": {
        "loci": ["HBG1_BCL11A_site", "MYC_PVT1promoter", "SHH_SBE2"],
        "why": (
            "three targets are already panel targets - HBG1/HBG2 through HBB_LCR, MYC through"
            " MYC_8q24, SHH through SHH_ZRS - and the ELEMENTS are different published elements in"
            " every case. This is the point rather than an oversight. At HBG and MYC the new element"
            " acts in the OPPOSITE direction to the panel's, which is the cleanest form the direction"
            " question can take: same target, same cell, opposite published sign. At SHH the panel's"
            " element is 979 kb away and unaskable (section 18), and this one is 456 kb away and"
            " askable, so the gene the benchmark could never ask about becomes a question"
        ),
        "risk": (
            "a shared target means a shared annotation neighbourhood, so `looked_up` and `eqtl` hits"
            " at these three are less independent of the panel's than the rest. Registered here; the"
            " provenance split already reports looked_up on its own line"
        ),
    },
    "coordinates": {
        "liftover": "SHH_SBE2 is an Ensembl GRCh37->GRCh38 assembly map of the published GRCh37"
        " interval chr7:156,061,051-156,061,848, the mechanism loci_candidates used for HBA_HS40",
        "rsid": "CDKN2A_9p21 and KITLG_blond are placed on an rsID resolved from Ensembl GRCh38 at"
        " authoring time and re-checked against Ensembl on every run",
        "anchored": "GATA2_plus95, TAL1_MuTE, SOST_VanBuchem, HBG1_BCL11A_site, MYC_PVT1promoter and"
        " TERT_promoter are anchored on a published offset from a GENCODE gene boundary, the"
        " mechanism docs/LOCI-BENCHMARK.md section 1 uses for SOX9 and H19 and section 19 for PTF1A."
        " Every one carries its own element_citation naming the gene boundary and the offset",
        "exhibits_only": "none. Every locus here has a coordinate this module can cite or derive, so"
        " unlike loci_candidates there is no ungraded exhibit and the graded set is all nine",
    },
    "reach_filter": {
        "how": "loci.read_reach, called not reimplemented: the published target's GENCODE gene BODY"
        " against the scorer's 1,048,576 bp input centred on the element. Free, no request",
        "rule": "a locus with no target in reach is unaskable, is never sent a request, and is"
        " reported as a reach fatality rather than as a miss",
        "expected": "0 of 9. Every element was placed before the filter was run and every published"
        " distance here is under 524 kb, so the filter is expected to kill nothing. That is the"
        " prediction and it is registered so a surprise counts",
    },
    "hit_rules": {
        "source": "loci.score_target, loci.score_cell, loci.score_direction, unchanged and unwrapped",
        "target": "strict: a derived layer's FIRST-ranked gene is one of the published targets. The"
        " lenient reading is reported as `among` and is never counted in any rate",
        "direction": "the deletion layer's action for a published target equals the published"
        " direction, and is judged only where a deletion names a published target",
        "cell": "the reader has the element open in a declared cell type, or GTEx ties the target to"
        " a declared tissue",
        "provenance": "derived, heuristic and looked_up on separate lines, as in both other sets",
    },
    "denominators": [
        "all nine loci",
        "loci where the model could answer (loci.read_reach askable), the same second denominator"
        " the panel reports as target_derived_where_the_model_could_answer",
        "both quoted beside the panel's 15/17 and 13/15 and the candidates' 8/9. No third denominator"
        " is invented, and no rate pools the three frames",
    ],
    "negative_controls": (
        "loci.pick_negatives, five matched windows per locus, the same four covariates - length"
        " exactly, GC within 0.04, distance to the nearest coding TSS within 35%, Zoonomia constrained"
        " fraction closest - drawn away from VISTA elements and GWAS hits and away from every locus in"
        " ALL THREE sets with 200 kb of flank. loci.candidate_windows reads loci.PANEL for its"
        " keep-outs, so `keep_out_all_three_sets` extends that global for the draw and restores it"
        " after; nothing is standardised here and every two-group comparison goes through"
        " genomeos.compare.standardised with input_presence before any coverage stratum"
    ),
    "requests": (
        "zero where the finished all-chromosome sweep has already deleted an element inside the"
        " published interval; one stated-interval deletion for each element it has not. Counted"
        " before anything is spent, reported per locus, and the plan refuses to spend on the rest"
    ),
    "predictions": {
        "target_derived": "6 to 8 of 9. Six of the nine elements sit inside, or within 10 kb of, their"
        " own published target, so the set is closer to the panel's flattered shape than to the"
        " candidates'. This is registered as a WEAKNESS of the set, not as a result",
        "heuristic": "7 of 9, and it is registered because it is the same weakness seen from the"
        " other side. Only two loci have a nearest-TSS trap - RNF32 at SHH_SBE2 (371 kb against"
        " SHH's 456 kb) and MEOX1 at SOST_VanBuchem (4.8 kb against SOST's 40 kb) - and the traps"
        " are geometry over GENCODE computed before scoring, exactly as LMBR1 and HBS1L were. So this"
        " set is expected to make the nearest-gene rule look BETTER than the candidates' 3/9 did and"
        " about as good as the panel's 11/17, and a benchmark's verdict on that rule is therefore a"
        " statement about which loci somebody chose",
        "direction_activators": "right at most of the seven activators, as at both other sets",
        "direction_repressors": (
            "THE REGISTERED QUESTION OF THIS SET, and n = 2 cannot decide it. Two published"
            " repressors is not a rate; it is two readings. Registered in advance: if both come back"
            " with the published sign, that is weak support that the model carries direction and the"
            " candidates' two sign errors were about those loci; if both come back inverted, that is"
            " weak support that the deletion layer's sign is not a reading of direction at all; one"
            " of each decides nothing. No claim about the sign will be made from this set alone, and"
            " no denominator pooling the repressors with the activators will be quoted"
        ),
        "cell": "judged at 2 of 9 only. Seven published tissues here - rostral diencephalon, T-ALL"
        " thymocyte, osteocyte, hair follicle, vascular wall, melanoma - are not reader cell types and"
        " mostly not GTEx tissues either, so the cell axis is close to unmeasurable on this set and"
        " is registered as such in advance",
        "naming_some_target": "near 90% at the matched controls, as everywhere else in this"
        " benchmark, so the claim will not separate and is not expected to",
        "value_on_syntax": "section 20 showed the panel's 9/17 is carried by its long elements. This"
        " set has one long element (SOST_VanBuchem, 52 kb) and eight short ones, so the registered"
        " expectation is that it looks like the candidates' 1/9 and not like the panel's 9/17, and"
        " that a set of eight short elements plus one long one cannot settle a length gradient. n is"
        " stated as too small to decide it BEFORE the number exists",
        "reach_fatalities": "0 of 9. See reach_filter above",
    },
    "what_would_falsify_the_run": (
        "a derived target rate at or below the chance floor (a coding gene drawn at random from the"
        " locus window); or TERT_promoter failing to name TERT, which is this set's built-in positive"
        " control - deleting a gene's own core promoter must move that gene, and if it does not, the"
        " reading is broken and nothing else in the run is interpretable"
    ),
    "aims_that_failed_before_the_run": {
        "a_non_coding_published_target": (
            "section 19 left the `predicted_coding`-first defect untested because every candidate"
            " target was protein coding. A published, perturbation-backed enhancer whose stated"
            " target is a lncRNA or miRNA was looked for and none was found clean enough to state:"
            " the well-known non-coding cases either act through parent-of-origin methylation"
            " (KCNQ1OT1, and H19 is already in the panel and already pending for that reason) or name"
            " the coding neighbour as the causal target when perturbed. So the defect is still"
            " untested after three sets, and that is a fact about the literature"
        ),
        "more_than_two_repressors": (
            "section 19 reported that published repressive elements with a perturbation behind them"
            " were looked for and none was clean enough to state. A second search found two. Two is"
            " the measurement: the direction axis cannot be balanced against seven activators, and"
            " what would balance it is a CRISPRi screen reporting de-repression, not more curation"
        ),
    },
}

# ------------------------------------------------------------------------------- the third set
THIRD: tuple[Expect, ...] = (
    Expect(
        locus="SHH_SBE2",
        chrom="chr7",
        element=(156_268_356, 156_269_154),
        window=(155_770_000, 156_770_000),
        classes=("program",),
        targets=("SHH",),
        direction="activates",
        tissues=("rostral diencephalon", "ventral forebrain", "zona limitans intrathalamica"),
        gtex_tissues=(),
        cells=(),
        distance=455_893,
        variants=(
            {
                "rsid": None,
                "pos": None,
                "effect": "a 444C>T transition inside SBE2 in a holoprosencephaly patient lowers"
                " enhancer activity; SIX3 binds the element and drives forebrain SHH",
            },
        ),
        nearest_gene_trap="RNF32",
        element_citation=(
            "Jeong et al. 2008, Nat Genet 40:1348: SBE2, 460 kb upstream of SHH. The interval is an"
            " Ensembl GRCh37->GRCh38 assembly map of the published GRCh37 element"
            " chr7:156,061,051-156,061,848, which lands 455,893 bp from GENCODE's SHH gene end"
            " (chr7:155,812,463) and so reproduces the published 460 kb offset to within 5 kb"
        ),
        citations=(
            "Jeong et al. 2008, Nat Genet 40:1348 (SIX3 activates SHH through SBE2; a point mutation"
            " in SBE2 from a holoprosencephaly patient abolishes the activation)",
            "Jeong et al. 2006, Development 133:761 (the SBE enhancers drive Shh in the ventral"
            " forebrain in transgenic mouse)",
            "Yao et al. 2016, Open Biology 6:160197 (the SHH regulatory domain extends 1 Mb upstream"
            " and carries several independent brain enhancers)",
        ),
        answer_from=("transgenic mouse reporters", "human pedigree", "reporter assays"),
        note=(
            "the locus section 18 said could not exist: an SHH enhancer the model can actually see."
            " SHH_ZRS is 979 kb away and UNASKABLE, and the panel graded the model wrong there for"
            " months. SBE2 is 456 kb away, SHH's body is inside the 1 Mb input, and the nearest"
            " coding TSS is RNF32 at 371 kb - a genuine trap at a genuine distance, which the panel's"
            " megabase loci could never provide because the model cannot reach them"
        ),
    ),
    Expect(
        locus="GATA2_plus95",
        chrom="chr3",
        element=(128_483_458, 128_483_958),
        window=(127_985_000, 128_985_000),
        classes=("program",),
        targets=("GATA2",),
        direction="activates",
        tissues=("haematopoietic stem cell", "haemogenic endothelium", "bone marrow"),
        gtex_tissues=("Whole_Blood",),
        cells=("K562",),
        distance=9_500,
        variants=(
            {
                "rsid": None,
                "pos": None,
                "effect": "germline point mutations and small deletions in the +9.5 element cause"
                " MonoMAC and Emberger syndrome with the GATA2 coding sequence intact",
            },
        ),
        nearest_gene_trap=None,
        element_source="stated",
        element_citation=(
            "Johnson et al. 2012, J Clin Invest 122:3692 and Hsu et al. 2013, Blood 121:3830: the"
            " +9.5 intronic enhancer, 9.5 kb into GATA2 from its transcription start. Anchored at"
            " that published offset from GENCODE's GATA2 gene end (chr3:128,493,208, minus strand),"
            " which places it in intron 5 as published, and widened to 500 bp"
        ),
        citations=(
            "Johnson et al. 2012, J Clin Invest 122:3692 (germline +9.5 enhancer mutations cause"
            " GATA2 deficiency with an intact coding sequence)",
            "Gao et al. 2013, J Exp Med 210:2833 (deleting the +9.5 element in mouse abolishes"
            " haematopoietic stem cell emergence)",
            "Hsu et al. 2013, Blood 121:3830 (a mutational hotspot in the +9.5 enhancer in MonoMAC"
            " and Emberger syndrome)",
        ),
        answer_from=("mouse enhancer deletion", "human pedigrees", "reporter assays"),
        note=(
            "an enhancer inside its own target's fifth intron, so the nearest coding TSS IS the"
            " published target and the locus tests nothing about reach. It is here for the"
            " perturbation - a targeted mouse deletion of 500 bp that removes the haematopoietic"
            " system - and it is one of the six loci that make this set's heuristic prediction high"
        ),
    ),
    Expect(
        locus="TAL1_MuTE",
        chrom="chr1",
        element=(47_239_475, 47_239_975),
        window=(46_740_000, 47_740_000),
        classes=("program",),
        targets=("TAL1",),
        direction="activates",
        tissues=("T-cell acute lymphoblastic leukaemia", "thymocyte"),
        gtex_tissues=(),
        cells=(),
        distance=7_500,
        variants=(
            {
                "rsid": None,
                "pos": None,
                "effect": "recurrent somatic 2-18 bp indels introduce a MYB binding motif and build a"
                " super-enhancer where none exists in normal cells; deleting it abolishes TAL1",
            },
        ),
        nearest_gene_trap=None,
        element_source="stated",
        element_citation=(
            "Mansour et al. 2014, Science 346:1373: the MuTE element, 7.5 kb upstream of the TAL1"
            " transcription start. Anchored at that published offset from GENCODE's TAL1 gene end"
            " (chr1:47,232,225, minus strand) and widened to 500 bp"
        ),
        citations=(
            "Mansour et al. 2014, Science 346:1373 (somatic mutations 7.5 kb upstream of TAL1 create"
            " a MYB site and an oncogenic super-enhancer; CRISPR deletion of the mutant element"
            " abolishes TAL1 expression in Jurkat)",
            "Navarro et al. 2015, Blood 126:2870 (the same element in further T-ALL primary samples)",
        ),
        answer_from=("CRISPR deletion in a human cell line", "somatic mutation in patients"),
        note=(
            "the only locus in any of the three sets where the published element does not exist in"
            " the reference at all: the enhancer is CREATED by a somatic indel. The model reads the"
            " reference, so what it is asked here is whether the sequence is poised to become an"
            " enhancer, which is a different question from the published one. Registered as a locus"
            " whose miss would be uninformative"
        ),
    ),
    Expect(
        locus="CDKN2A_9p21",
        chrom="chr9",
        element=(22_124_228, 22_124_728),
        window=(21_625_000, 22_625_000),
        classes=("program",),
        targets=("CDKN2A", "CDKN2B"),
        direction="activates",
        tissues=("vascular smooth muscle", "endothelium", "aorta"),
        gtex_tissues=("Artery_Aorta", "Artery_Tibial", "Artery_Coronary"),
        cells=(),
        distance=114_903,
        variants=(
            {
                "rsid": "rs10757278",
                "pos": 22_124_478,
                "ref": "A",
                "alt": "G",
                "effect": "the G risk allele disrupts a STAT1 binding site in the ECAD9 enhancer and"
                " abolishes the interferon-gamma-dependent contact with the CDKN2A/B locus",
            },
        ),
        nearest_gene_trap=None,
        citations=(
            "Visel et al. 2010, Nature 464:409 (deleting the 70 kb mouse interval orthologous to the"
            " human 9p21 risk interval lowers Cdkn2a and Cdkn2b)",
            "Harismendy et al. 2011, Nature 470:264 (rs10757278 sits in a STAT1-binding enhancer;"
            " the risk allele abolishes STAT1 binding and the enhancer's contact with CDKN2B)",
            "Helgadottir et al. 2007, Science 316:1491 (the 9p21 coronary artery disease interval)",
        ),
        answer_from=("mouse interval deletion", "chromatin conformation", "human association"),
        note=(
            "the best-replicated non-coding risk interval in the genome and the one place in this set"
            " where the perturbation is a large mouse deletion rather than the element itself. The"
            " nearest coding TSS is CDKN2B at 115 kb, which IS a published target, so the heuristic"
            " hits without the reach question ever being asked"
        ),
    ),
    Expect(
        locus="KITLG_blond",
        chrom="chr12",
        element=(88_934_308, 88_934_808),
        window=(88_435_000, 89_435_000),
        classes=("program",),
        targets=("KITLG",),
        direction="activates",
        tissues=("hair follicle", "melanocyte", "skin"),
        gtex_tissues=("Skin_Sun_Exposed_Lower_leg", "Skin_Not_Sun_Exposed_Suprapubic"),
        cells=(),
        distance=353_457,
        variants=(
            {
                "rsid": "rs12821256",
                "pos": 88_934_558,
                "ref": "T",
                "alt": "C",
                "effect": "the derived C allele weakens a LEF1 site and lowers enhancer output about"
                " 20%; mice carrying the human element with C have lighter coats",
            },
        ),
        nearest_gene_trap=None,
        citations=(
            "Guenther et al. 2014, Nat Genet 46:748 (a knock-in mouse carrying the human KITLG"
            " enhancer shows the derived allele lowers activity and lightens coat colour)",
            "Sulem et al. 2007, Nat Genet 39:1443 (rs12821256 and blond hair in northern Europeans)",
        ),
        answer_from=("mouse enhancer knock-in", "transgenic reporters", "human association"),
        note=(
            "a 354 kb link where the nearest coding TSS is STILL the published target, because the"
            " neighbourhood is empty. So it is a long-range locus that the nearest-gene rule gets"
            " right for free - the mirror image of MYB_HMIP2 in the second set, and a reminder that"
            " reach and difficulty are different axes"
        ),
    ),
    Expect(
        locus="SOST_VanBuchem",
        chrom="chr17",
        element=(43_666_737, 43_718_737),
        window=(43_193_000, 44_193_000),
        classes=("program",),
        targets=("SOST",),
        direction="activates",
        tissues=("osteocyte", "bone", "cartilage"),
        gtex_tissues=(),
        cells=(),
        distance=35_000,
        variants=(
            {
                "rsid": None,
                "pos": None,
                "effect": "a homozygous 52 kb deletion downstream of an intact SOST silences it in"
                " bone and causes van Buchem disease",
            },
        ),
        nearest_gene_trap="MEOX1",
        element_source="stated",
        element_citation=(
            "Balemans et al. 2002, J Med Genet 39:91: a 52 kb deletion beginning 35 kb downstream of"
            " SOST. Anchored at that published offset and length from GENCODE's SOST 3' end"
            " (chr17:43,753,737, minus strand), which places the interval in the SOST-MEOX1"
            " intergenic region where the publication puts it. The 250 bp ECR5 enhancer that carries"
            " the activity sits inside it (Loots 2005), and a 52 kb interval absorbs the placement"
            " error a 250 bp one could not"
        ),
        citations=(
            "Balemans et al. 2002, J Med Genet 39:91 (the 52 kb van Buchem deletion, with SOST intact)",
            "Loots et al. 2005, Genome Res 15:928 (ECR5 inside the deletion is a bone-specific SOST"
            " enhancer in transgenic mouse)",
            "Collette et al. 2012, PNAS 109:14092 (targeted deletion of the Sost distal enhancer in"
            " mouse raises bone mass)",
        ),
        answer_from=("human deletion", "mouse enhancer deletion", "transgenic reporters"),
        note=(
            "the length-matched locus section 20 asked for: 52 kb, which is the panel's long-element"
            " regime and not the candidates' uniform 500 bp. It is also the set's second"
            " nearest-gene trap - MEOX1's TSS is 4.8 kb from the interval and SOST's is 40 kb - and"
            " the interval is a real human deletion rather than one anybody drew"
        ),
    ),
    Expect(
        locus="HBG1_BCL11A_site",
        chrom="chr11",
        element=(5_249_724, 5_250_224),
        window=(4_750_000, 5_750_000),
        classes=("program",),
        targets=("HBG1", "HBG2"),
        direction="represses",
        tissues=("erythroid", "bone marrow", "fetal liver"),
        gtex_tissues=("Whole_Blood",),
        cells=("K562",),
        distance=115,
        variants=(
            {
                "rsid": None,
                "pos": None,
                "effect": "the -115 and -114 hereditary persistence of fetal haemoglobin point"
                " mutations destroy the BCL11A TGACCA motif, and erythroid cells edited at the motif"
                " raise HbF - so the element's published job is to hold HBG down",
            },
        ),
        nearest_gene_trap=None,
        element_source="stated",
        element_citation=(
            "Martyn et al. 2018, Nat Genet 50:498: BCL11A binds TGACCA at -115 of the HBG1 and HBG2"
            " promoters. Anchored 115 bp upstream of GENCODE's HBG1 gene end (chr11:5,249,859, minus"
            " strand) and widened to 500 bp, which covers -365 to +135 and so holds the motif even if"
            " the annotated cap site is tens of bases off the published one"
        ),
        citations=(
            "Martyn et al. 2018, Nat Genet 50:498 (BCL11A represses HBG through the -115 TGACCA"
            " motif; the HPFH mutations there destroy the motif)",
            "Liu et al. 2018, Cell 173:430 (BCL11A binds the HBG promoters directly at -115)",
            "Traxler et al. 2016, Nat Med 22:987 (editing the HBG promoters at the site raises fetal"
            " haemoglobin in human erythroid cells)",
        ),
        answer_from=("genome editing in human erythroid cells", "human pedigrees", "protein binding"),
        note=(
            "one of the two published REPRESSORS this set exists for, and the cleanest form the"
            " direction question can take: HBG1 and HBG2 are also HBB_LCR's published targets in the"
            " panel, in the same cell, with the opposite published sign. The element is"
            " promoter-proximal, so it tests direction and nothing about reach, and its nearest"
            " coding TSS is HBG1 itself"
        ),
    ),
    Expect(
        locus="MYC_PVT1promoter",
        chrom="chr8",
        element=(127_794_262, 127_794_762),
        window=(127_295_000, 128_295_000),
        classes=("program",),
        targets=("MYC",),
        direction="represses",
        tissues=("breast carcinoma", "colorectal carcinoma", "epithelium"),
        gtex_tissues=(),
        cells=(),
        distance=58_829,
        variants=(
            {
                "rsid": None,
                "pos": None,
                "effect": "no point variant: the perturbation is CRISPR deletion of the promoter,"
                " which RAISES MYC because the promoter competes for the shared enhancers",
            },
        ),
        nearest_gene_trap=None,
        element_source="stated",
        element_citation=(
            "Cho et al. 2018, Cell 173:1398: the PVT1 promoter is a tumour-suppressor boundary"
            " element that competes with the MYC promoter for the same enhancers. Anchored on"
            " GENCODE's PVT1 gene start (chr8:127,794,512, plus strand), widened to 500 bp"
        ),
        citations=(
            "Cho et al. 2018, Cell 173:1398 (CRISPR deletion of the PVT1 promoter raises MYC"
            " expression and tumour growth; the lncRNA transcript is not what does it)",
            "Tseng et al. 2014, Nature 512:82 (PVT1 and MYC are co-amplified and functionally"
            " coupled at 8q24)",
        ),
        answer_from=("CRISPR deletion in human cell lines", "enhancer contact frequency"),
        note=(
            "the second published REPRESSOR, and the second deliberate target overlap: MYC is also"
            " MYC_8q24's published target in the panel, through an element 335 kb the other side of"
            " it that ACTIVATES. Here deleting the element raises MYC. A deletion scorer that reports"
            " a sign has to disagree with itself across these two elements or it is not reading"
            " direction at all"
        ),
    ),
    Expect(
        locus="TERT_promoter",
        chrom="chr5",
        element=(1_294_950, 1_295_250),
        window=(795_000, 1_795_000),
        classes=("program",),
        targets=("TERT",),
        direction="activates",
        tissues=("melanoma", "glioblastoma", "bladder urothelium"),
        gtex_tissues=(),
        cells=(),
        distance=124,
        variants=(
            {
                "rsid": None,
                "pos": 1_295_113,
                "ref": "G",
                "alt": "A",
                "effect": "C228T, the commonest non-coding somatic mutation in human cancer: it"
                " creates a GABP ETS site and raises TERT; reverting it in a mutant cell line lowers"
                " TERT back",
            },
            {
                "rsid": None,
                "pos": 1_295_135,
                "effect": "C250T, the second hotspot, 22 bp away and the same mechanism",
            },
        ),
        nearest_gene_trap=None,
        element_source="stated",
        element_citation=(
            "Huang et al. 2013, Science 339:957 and Horn et al. 2013, Science 339:959: the C228T and"
            " C250T hotspots at -124 and -146 of TERT, GRCh38 chr5:1,295,113 and chr5:1,295,135. The"
            " interval is the 300 bp core promoter spanning both hotspots and GENCODE's TERT gene end"
            " (chr5:1,295,086, minus strand)"
        ),
        citations=(
            "Huang et al. 2013, Science 339:957 and Horn et al. 2013, Science 339:959 (the recurrent"
            " TERT promoter mutations in melanoma)",
            "Chiba et al. 2015, eLife 4:e07918 (CRISPR reversal of the promoter mutation in a mutant"
            " cell line lowers TERT and shortens telomeres)",
            "Bell et al. 2015, Science 348:1036 (GABP binds the mutant promoter and drives TERT)",
        ),
        answer_from=("CRISPR reversal in a human cell line", "recurrent somatic mutation", "protein binding"),
        note=(
            "this set's POSITIVE CONTROL, and it is registered as one rather than presented as a"
            " result. It is a core promoter and not an enhancer: deleting a gene's own promoter must"
            " move that gene, so a miss here means the reading is broken and the rest of the run is"
            " uninterpretable. Its hit is counted in the rate - excluding it after the fact would be"
            " exactly the tuning this benchmark forbids - and the rate is quoted with the note that"
            " one of the nine is a control"
        ),
    ),
)

#: this set's nearest-coding-TSS traps, computed from GENCODE geometry BEFORE anything was scored, so
#: a heuristic hit or miss at these two is arithmetic and not a result. The other seven loci have a
#: published target as their own nearest coding TSS, which is the registered weakness of the set.
TRAPS: dict[str, dict[str, Any]] = {
    "SHH_SBE2": {"trap": "RNF32", "trap_distance": 371_126, "target_distance": 455_893},
    "SOST_VanBuchem": {"trap": "MEOX1", "trap_distance": 4_808, "target_distance": 40_054},
}

#: the two published repressors. Reported on their own line, never pooled into a direction rate with
#: the seven activators, and n = 2 is declared too small to decide anything in PREREGISTRATION.
REPRESSORS: frozenset[str] = frozenset({"HBG1_BCL11A_site", "MYC_PVT1promoter"})

#: the locus whose failure invalidates the run rather than scoring as a miss.
POSITIVE_CONTROL = "TERT_promoter"

#: the hunk `genomeos/benchmark/loci.py` needs and that this lane is not permitted to make. Until it
#: lands, `register_stated_intervals` applies it at runtime; the two are the same one-line change.
NEEDED_LOCI_HUNK = {
    "file": "genomeos/benchmark/loci.py",
    "constant": "STATED_INTERVAL_RESULTS",
    "from": ("loci_stated_intervals", "loci_candidate_intervals"),
    "to": ("loci_stated_intervals", "loci_candidate_intervals", "loci_third_intervals"),
    "why": "loci._deletion_rows reads the panels' self-stated intervals from that tuple and labels"
    " them `stated_interval`. Without the third name, a deletion this set paid for is written, is"
    " never read back, and the locus grades as though no request had been made",
}


def register_stated_intervals() -> None:
    """Make `loci._deletion_rows` read this set's stated intervals too.

    This is `NEEDED_LOCI_HUNK` applied at runtime rather than in the file, because the lane that
    built this set is not permitted to edit `genomeos/benchmark/loci.py`. It is idempotent, it only
    ever appends, and it is the whole of the change: when the hunk lands, this becomes a no-op.
    """
    if INTERVALS not in loci.STATED_INTERVAL_RESULTS:
        loci.STATED_INTERVAL_RESULTS = (*loci.STATED_INTERVAL_RESULTS, INTERVALS)


@contextlib.contextmanager
def keep_out_all_three_sets():
    """Draw negative controls away from every locus in all three sets, not only the seventeen.

    `loci.candidate_windows` builds its keep-out list from the module-level `loci.PANEL`, so the only
    way to extend it without touching that file is to extend the global for the draw. The extension
    is restored on the way out, including on an exception, and nothing else in `loci` reads `PANEL`
    except `panel_by_name` and `build`'s own default argument, neither of which runs inside this.
    """
    original = loci.PANEL
    extra = tuple(e for e in (*loci_candidates.CANDIDATES, *THIRD) if e not in original)
    loci.PANEL = (*original, *extra)
    try:
        yield loci.PANEL
    finally:
        loci.PANEL = original


def by_name() -> dict[str, Expect]:
    return {e.locus: e for e in THIRD}


def plan(results_dir: Path = RESULTS_DIR) -> list[dict[str, Any]]:
    """Which loci the model could answer and what each would cost. No request, no network.

    Two free questions decide the budget, in this order: can the model see the target at all
    (`loci.read_reach`), and has anybody already deleted an element inside the published interval - if
    so the finished all-chromosome sweep holds the answer and it costs nothing.
    """
    register_stated_intervals()
    rows: list[dict[str, Any]] = []
    for chrom in sorted({e.chrom for e in THIRD}):
        ch = loci.Chromosome(chrom, results_dir)
        try:
            for e in (x for x in THIRD if x.chrom == chrom):
                reach = loci.read_reach(ch, e)
                s, en = e.element
                scored = loci._deletion_rows(chrom, s, en, results_dir)  # noqa: SLF001
                rows.append(
                    {
                        "locus": e.locus,
                        "chrom": chrom,
                        "element": [s, en],
                        "length": en - s,
                        "targets": list(e.targets),
                        "direction": e.direction,
                        "askable": reach["askable"],
                        "targets_in_reach": reach["targets_in_reach"],
                        "targets_out_of_reach": reach["targets_out_of_reach"],
                        "graded": reach["askable"],
                        "ccres_over_element": len(ch.ccres_in(s, en)),
                        "already_scored_over_element": len(scored),
                        "element_source": e.element_source,
                        "requests": 0 if (not reach["askable"] or scored) else 1,
                    }
                )
        finally:
            ch.close()
    return sorted(rows, key=lambda r: r["locus"])


def graded(rows: list[dict[str, Any]]) -> tuple[Expect, ...]:
    """The loci the run grades: every one the reach filter keeps. There are no exhibits here."""
    keep = {r["locus"] for r in rows if r["graded"]}
    return tuple(e for e in THIRD if e.locus in keep)


def reach_fatalities(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """How many loci the free filter killed - a number in its own right, whatever it is."""
    dead = [r for r in rows if not r["askable"]]
    return {
        "loci": len(rows),
        "died_at_the_reach_filter": len(dead),
        "registered_expectation": PREREGISTRATION["reach_filter"]["expected"],
        "which": [
            {
                "locus": r["locus"],
                "target": r["targets_out_of_reach"][0]["target"] if r["targets_out_of_reach"] else None,
                "distance_kb": (r["targets_out_of_reach"][0].get("distance_bp", 0) // 1000)
                if r["targets_out_of_reach"]
                else None,
            }
            for r in dead
        ],
        "requests_saved": len(dead),
        "beside_the_other_sets": {
            "the_seventeen": "3 of 17 unaskable (SHH_ZRS, SOX9_PierreRobin, FTO_IRX3's IRX5)",
            "the_nine": "2 of 11 died, both exhibits placed on a bare offset",
            "the_third_set": f"{len(dead)} of {len(rows)}",
        },
        "reading": (
            "a locus outside the scorer's 1,048,576 bp input is unaskable, not a miss. Across the"
            " three frames the filter kills the famous megabase loci and almost nothing else, which"
            " is the same measurement each time: a published enhancer-gene link is usually well"
            " inside the model's reach, and the panel's megabase cases were chosen for being"
            " remarkable rather than for being typical"
        ),
    }


def direction_readings(result: dict[str, Any]) -> dict[str, Any]:
    """The direction axis split by published sign - the question this set was assembled to ask.

    Repressors and activators are reported apart and never pooled. `PREREGISTRATION` says in advance
    that n = 2 on the repressor side decides nothing, and that sentence is carried here beside the
    numbers so it cannot be read off without it.
    """
    rows = []
    for r in result["loci"]:
        d = ((r.get("score") or {}).get("scored") or {}).get("direction") or {}
        rows.append(
            {
                "locus": r["locus"],
                "published": r["expected"]["direction"],
                "repressor": r["locus"] in REPRESSORS,
                "judged": bool(d.get("judged")),
                "hit": bool(d.get("hit_derived")),
                "model_action": d.get("model_action"),
            }
        )
    judged = [r for r in rows if r["judged"]]

    def tally(which: list[dict[str, Any]]) -> dict[str, Any]:
        return {"k": sum(1 for r in which if r["hit"]), "n": len(which)}

    return {
        "per_locus": rows,
        "activators": tally([r for r in judged if not r["repressor"]]),
        "repressors": tally([r for r in judged if r["repressor"]]),
        "registered": PREREGISTRATION["predictions"]["direction_repressors"],
        "not_pooled": (
            "no rate combining the two lines is quoted. The published sign is the thing under test,"
            " so pooling the arms would hide exactly the asymmetry the set was built to expose"
        ),
    }


def positive_control(result: dict[str, Any]) -> dict[str, Any]:
    """Did deleting TERT's own core promoter name TERT. A no here invalidates the rest of the run."""
    row = next((r for r in result["loci"] if r["locus"] == POSITIVE_CONTROL), None)
    if row is None:
        return {"locus": POSITIVE_CONTROL, "ran": False, "passed": None}
    target = ((row.get("score") or {}).get("scored") or {}).get("target") or {}
    deletion = (row.get("readings") or {}).get("deletion") or {}
    return {
        "locus": POSITIVE_CONTROL,
        "ran": True,
        "passed": bool(target.get("hit_derived")),
        "hit_by": target.get("hit_derived_by") or [],
        "deletion_named_first": deletion.get("target"),
        "deletion_log2_fold_change": (deletion.get("targets") or [{}])[0].get("best"),
        "registered": PREREGISTRATION["what_would_falsify_the_run"],
    }


def compare_with_controls(result: dict[str, Any]) -> dict[str, Any]:
    """The four claims at this set against its matched windows, standardised.

    `genomeos.compare` does the standardising; this only shapes the rows. `input_presence` runs
    before any coverage stratum, because the trap section 19 walked into is that GC, distance and
    constraint say nothing about whether anybody ever spent a measurement on the window.
    """
    from genomeos.compare import Strata, imbalance, input_presence, standardised

    def row(gc, dist, con, claims) -> dict[str, Any] | None:
        if gc is None or dist is None:
            return None
        return {"gc": gc, "log_distance": (dist or 1) ** 0.5, "constrained": con or 0.0, **claims}

    def presence(con, claims) -> dict[str, Any]:
        # present-or-None, which is what input_presence counts. A boolean False would read as an
        # input that is there and says no, and that is the substitution the function exists to catch.
        return {
            "deletion_scored": True if claims.get("deletion_scored") else None,
            "eqtl_present": True if claims.get("eqtl_present") else None,
            "constrained": con,
        }

    claimed = [(r, loci.claims(r["readings"])) for r in result["loci"]]
    targets = [
        row(r.get("gc"), r.get("distance_to_coding_tss"), r.get("constrained_fraction"), c)
        for r, c in claimed
    ]
    controls = [
        row(w.get("gc"), w.get("distance_to_coding_tss"), w.get("constrained_fraction"), w["claims"])
        for w in result["negatives"]
    ]
    t_presence = [presence(r.get("constrained_fraction"), c) for r, c in claimed]
    c_presence = [presence(w.get("constrained_fraction"), w["claims"]) for w in result["negatives"]]
    targets = [r for r in targets if r]
    controls = [r for r in controls if r]
    strata = Strata(gc=(0.40, 0.48, 0.56), log_distance=(100.0, 250.0, 500.0), constrained=(0.02, 0.10))
    covered = Strata(
        gc=(0.40, 0.48, 0.56),
        log_distance=(100.0, 250.0, 500.0),
        constrained=(0.02, 0.10),
        deletion_scored=(0.5,),
    )
    for r in targets + controls:
        r["deletion_scored"] = float(bool(r["deletion_scored"]))
    out: dict[str, Any] = {
        "reading": "direct standardisation on the covariates the windows were matched on; the"
        " `*_given_coverage` lines add whether an element inside the window has actually been"
        " deleted, which is the stratum sections 8 and 19 both needed and neither covariate carries",
        "input_presence": input_presence(
            t_presence,
            c_presence,
            {
                "a deletion spent inside the window": "deletion_scored",
                "a GTEx eQTL distilled for the window": "eqtl_present",
                "the mammalian constraint track over the window": "constrained",
            },
        ),
        "why_input_presence_first": (
            "docs/LOCI-BENCHMARK.md section 19: standardising on GC, distance and constraint does NOT"
            " remove a coverage artefact, because none of them says whether anybody spent a"
            " measurement on that window. input_presence is asked before any coverage stratum, so a"
            " null under `*_given_coverage` is known to be a test the claim passed rather than a test"
            " nobody could administer because the input was everywhere already"
        ),
        "imbalance": imbalance(targets, controls, strata),
    }
    for claim in ("target", "cell", "direction", "storage"):
        out[claim] = standardised(targets, controls, strata, hit=claim)
        out[f"{claim}_given_coverage"] = standardised(targets, controls, covered, hit=claim)
    return out


def beside_the_other_two(result: dict[str, Any], results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """This set's headline rates beside the panel's 15/17 and the candidates' 8/9, all three named."""
    panel = (load_result("loci_benchmark", results_dir) or {}).get("aggregate") or {}
    cands = (load_result(loci_candidates.NAME, results_dir) or {}).get("aggregate") or {}
    mine = result["aggregate"]

    def three(key: str) -> dict[str, Any]:
        return {
            "the_third_set": {k: (mine.get(key) or {}).get(k) for k in ("k", "n", "rate")},
            "the_nine_candidates": {k: (cands.get(key) or {}).get(k) for k in ("k", "n", "rate")},
            "the_seventeen_panel": {k: (panel.get(key) or {}).get(k) for k in ("k", "n", "rate")},
        }

    return {
        "frames": {
            "the_seventeen_panel": "loci.PANEL, hand-curated famous loci, 2026-09-13",
            "the_nine_candidates": "loci_candidates.CANDIDATES, perturbation-sourced, 2026-09-17",
            "the_third_set": "loci_third.THIRD, perturbation-sourced and chosen for a published"
            " repressive direction and for element lengths spanning the panel's range, 2026-09-21",
        },
        "target_derived": three("target_derived"),
        "target_derived_where_the_model_could_answer": three("target_derived_where_the_model_could_answer"),
        "target_heuristic": three("target_heuristic"),
        "target_looked_up": three("target_looked_up"),
        "chance_floor": {
            "the_third_set": mine.get("target_by_chance"),
            "the_nine_candidates": cands.get("target_by_chance"),
            "the_seventeen_panel": panel.get("target_by_chance"),
        },
        "reading": (
            "three frames, three rates, never one. The three sets were chosen by different people on"
            " different days for different reasons - fame, perturbation, and direction - and a pooled"
            " rate would hide that. A disagreement between them is a measurement of how much the"
            " benchmark's headline depends on who picked the loci, which is the only thing three"
            " small sets can measure that one large one cannot"
        ),
    }


def run(
    results_dir: Path = RESULTS_DIR,
    network: bool = True,
    negatives: bool = True,
    gtex_dir: Path = GTEX_DIR,
    progress=None,
) -> dict[str, Any]:
    """Plan, filter by reach, then read and score the graded loci through loci.build."""
    register_stated_intervals()
    rows = plan(results_dir)
    panel = graded(rows)
    with keep_out_all_three_sets():
        out = loci.build(
            panel=panel,
            results_dir=results_dir,
            gtex_dir=gtex_dir,
            network=network,
            negatives=negatives,
            progress=progress,
        )
    out["result"] = NAME
    out["preregistration"] = PREREGISTRATION
    out["plan"] = rows
    out["reach"] = reach_fatalities(rows)
    out["traps"] = TRAPS
    out["repressors"] = sorted(REPRESSORS)
    out["directions"] = direction_readings(out)
    out["positive_control"] = positive_control(out)
    out["needed_loci_hunk"] = NEEDED_LOCI_HUNK
    out["against_controls"] = compare_with_controls(out)
    out["beside_the_other_two"] = beside_the_other_two(out, results_dir)
    out["note"] = (
        "A THIRD, separately registered set of published enhancer-gene loci, graded by"
        " genomeos.benchmark.loci with no change to any scorer or hit rule. Reported beside the"
        " seventeen and the nine with each frame named, never pooled with either. " + loci.NOTE
    )
    save_result(NAME, out, results_dir)
    return out
