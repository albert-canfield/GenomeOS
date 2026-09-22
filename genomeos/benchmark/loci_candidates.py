# SPDX-License-Identifier: AGPL-3.0-or-later
"""Eleven more published enhancer-gene loci, pre-registered, reach-filtered and graded.

docs/LOCI-BENCHMARK.md's own conclusion, reached twice (section 8 at twelve loci and section 16 at
one rearrangement), is that what a reading of the same locus buys is nothing and what more loci buy
is everything. This module is that: candidates assembled from published *perturbation* experiments -
mouse or human enhancer deletions, CRISPR interference, transgenic reporters, and disease-causing
enhancer point mutations - never from association alone.

It adds no scorer. Every reading, every hit rule, every matched control and every rate comes from
`genomeos.benchmark.loci`, called with a different panel; the only thing here is the expectations,
the filter that decides which of them are worth a request, and the registration below.

**The reach filter comes first and it is free.** `loci.read_reach` asks whether the published
target's GENCODE *body* overlaps the scorer's 1,048,576 bp input centred on the element. It is a
body test, not a TSS-distance test, and the difference is load-bearing: on the finished chr21 sweep,
79 elements name a coding gene whose TSS is 524 kb or more away, and every one of them is a long
gene (mostly RUNX1, 1.22 Mb) whose body runs through the window. A filter written in TSS distance
would throw those away and would fail silently in the direction that looks safe. This module calls
`read_reach` rather than writing a second one, and `tests/test_loci_candidates.py` pins the RUNX1
case so a rewrite in TSS distance cannot pass.

The candidate list, the hit rules, the denominators, the controls and the predictions are fixed in
`PREREGISTRATION` below and were written before anything was scored. Nothing here is adjusted
afterwards: a negative or an inconclusive result is the outcome, not a reason to re-cut.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from genomeos.benchmark import loci
from genomeos.benchmark.loci import Expect
from genomeos.results import RESULTS_DIR, load_result, save_result

#: the result these candidates' stated-interval deletions are written to. `loci._deletion_rows`
#: reads it beside the panel's own, and labels every row `stated_interval`, so a hit on an interval
#: this module drew can never be pooled with a hit on one ENCODE drew.
INTERVALS = "loci_candidate_intervals"
NAME = "loci_candidates"
#: this panel's own distilled GTEx hits. `loci.stream_gtex` clears a directory whose window manifest
#: has changed, so sharing the seventeen-locus panel's cache would delete it every time either set
#: ran, in a checkout several sessions share. A directory per panel, and neither touches the other.
GTEX_DIR = Path("data/knowledge/loci_candidates")

# ------------------------------------------------------------------- the registration, before any run
PREREGISTRATION: dict[str, Any] = {
    "written": "2026-09-17, before any candidate was read through any layer or any request was spent",
    "candidates": 11,
    "how_they_were_chosen": (
        "published enhancer-gene links established by perturbation - CRISPR or mouse deletion of the"
        " element, CRISPR interference, transgenic reporter, or a disease-causing point mutation in"
        " the element - never by association alone. No locus already in loci.PANEL, and no locus whose"
        " published element the panel already covers"
    ),
    "coordinates": (
        "six elements are placed on an rsID resolved from Ensembl GRCh38 at authoring time and"
        " re-checked against Ensembl on every run; two are Ensembl assembly-map liftovers of a"
        " published GRCh37 interval (HBA_HS40, PAX6_SIMO); one is anchored on a published offset from"
        " a GENCODE gene end (PTF1A_enhancer), the mechanism docs/LOCI-BENCHMARK.md section 1 already"
        " uses for SOX9 and H19; two are anchored approximately and are exhibits only, see below"
    ),
    "reach_filter": {
        "how": "loci.read_reach: the published target's GENCODE gene BODY against the scorer's"
        " 1,048,576 bp input centred on the element. Free, no request, no network",
        "why_body_not_tss": "of 5,174 chr21 elements with a coding target in the finished sweep, 79"
        " name a gene whose TSS is 524 kb or more away and all 79 are genes whose body overlaps the"
        " window. A TSS-distance filter would discard answerable loci and look safe doing it",
        "rule": "a candidate with no target in reach is unaskable, is never sent a request, and is"
        " reported as a reach fatality rather than as a miss",
    },
    "exhibits_only": {
        "loci": ["POU3F4_DFN3", "MYC_BENC"],
        "why": "their elements are placed on a published offset alone (900 kb and 1.7 Mb) with no"
        " citable coordinate, and docs/LOCI-BENCHMARK.md section 13 withdrew a case for exactly that."
        " They are kept because the reach verdict is the one question an approximate placement can"
        " still answer - both targets are short genes and both offsets exceed 524 kb by a factor -"
        " and they are excluded from every scored rate whatever the filter says",
    },
    "graded_set": "the candidates that pass the reach filter and are not exhibits-only",
    "hit_rules": {
        "source": "loci.score_target, loci.score_cell, loci.score_direction, unchanged and unwrapped",
        "target": "strict: a derived layer's FIRST-ranked gene is one of the published targets. The"
        " lenient reading (a published target anywhere in the layer's list) is reported as `among`"
        " and is never counted in any rate",
        "direction": "the deletion layer's action for a published target equals the published"
        " direction, and is judged only where a deletion names a published target",
        "cell": "the reader has the element open in a declared cell type, or GTEx ties the target to"
        " a declared tissue",
        "provenance": "derived, heuristic (nearest coding TSS in the node) and looked_up are counted"
        " on separate lines, as in the panel; a looked-up hit alone is never a derived hit",
    },
    "denominators": [
        "all graded candidates",
        "graded candidates where the model could answer (loci.read_reach askable), reported as the"
        " panel reports target_derived_where_the_model_could_answer",
        "both quoted beside the panel's own 15/17 and 13/15. No third denominator is invented",
    ],
    "negative_controls": (
        "loci.pick_negatives, five matched windows per graded locus, the same four covariates the"
        " panel matches on - length exactly, GC within 0.04, distance to the nearest coding TSS"
        " within 35%, Zoonomia constrained fraction closest of the candidates - drawn away from every"
        " panel locus AND every candidate locus with 200 kb of flank, away from VISTA elements and"
        " GWAS hits, and read through exactly the same layers. Two-group comparisons go through"
        " genomeos.compare.standardised; nothing is standardised here"
    ),
    "requests": (
        "zero where an annotation already drew an interval over the published element, because the"
        " all-element sweep has finished every chromosome and its answer is already cached; one"
        " stated-interval deletion per element no annotation covers. The count is checked before"
        " anything is spent and reported per locus"
    ),
    "predictions": {
        "target_derived": "5 to 7 of the graded loci. The panel's 15/17 is flattered by four loci"
        " where the target is the gene the element sits inside; this set has fewer of those",
        "heuristic": "4 of 9, and this is registered because it is a weakness of the set, not a"
        " result: RET, IL2RA, PTF1A and PAX6 have the published target as their nearest coding TSS,"
        " so a short-range candidate set makes the nearest-TSS baseline look better than the panel"
        " does. The traps (HBS1L, PSRC1, NPRL3, FAM83C, ENSG00000289700) are geometry, computed"
        " before scoring, exactly as the panel computed LMBR1 and KCNJ2",
        "naming_some_target": "near 90% at the matched controls, as everywhere else in this"
        " benchmark, so the claim will not separate and is not expected to",
        "direction": "right wherever it is judged: every candidate is an activator, and the panel"
        " gets 5 of 5 on direction where a deletion names the target",
        "value_on_syntax": "the only claim that has ever separated the panel from its controls, so"
        " it is the one to watch: about half the candidates against about a sixth of the controls",
        "reach_fatalities": "2 of 11, both exhibits. The panel has three megabase-reach loci out of"
        " seventeen and this set has almost none, which says the panel's long-range cases were"
        " chosen for being remarkable rather than for being typical",
    },
    "what_would_falsify_the_run": (
        "a derived target rate at or below the chance floor (a coding gene drawn at random from the"
        " locus window), or a control rate on the value-on-syntax claim that matches the candidates'"
    ),
}

# ------------------------------------------------------------------------------- the candidates
CANDIDATES: tuple[Expect, ...] = (
    Expect(
        locus="RET_MCS97",
        chrom="chr10",
        element=(43_086_358, 43_086_858),
        window=(42_590_000, 43_590_000),
        classes=("program",),
        targets=("RET",),
        direction="activates",
        tissues=("enteric neural crest", "developing gut"),
        gtex_tissues=(),
        cells=(),
        distance=9_540,
        variants=(
            {
                "rsid": "rs2435357",
                "pos": 43_086_608,
                "ref": "C",
                "alt": "T",
                "effect": "T (the risk allele) disrupts a SOX10 site and lowers RET enhancer activity",
            },
        ),
        nearest_gene_trap=None,
        citations=(
            "Emison et al. 2005, Nature 434:857 (a common non-coding RET variant in a conserved"
            " intron 1 enhancer, MCS+9.7, carries most of the Hirschsprung risk)",
            "Emison et al. 2010, Am J Hum Genet 87:60 (rs2435357 is the causal allele and it reduces"
            " enhancer activity)",
            "Chatterjee et al. 2016, Cell 167:355 (deleting the enhancer in mouse lowers Ret; the"
            " variant abolishes SOX10 binding)",
        ),
        answer_from=("mouse enhancer deletion", "reporter assays", "human pedigrees"),
        note=(
            "an enhancer inside its own target's first intron: the nearest coding TSS IS the published"
            " target, so the heuristic is expected to hit and the locus tests nothing about reach"
        ),
    ),
    Expect(
        locus="IRF6_MCS97",
        chrom="chr1",
        element=(209_815_675, 209_816_175),
        window=(209_315_000, 210_315_000),
        classes=("program",),
        targets=("IRF6",),
        direction="activates",
        tissues=("oral periderm", "embryonic ectoderm", "keratinocyte"),
        gtex_tissues=("Skin_Sun_Exposed_Lower_leg", "Skin_Not_Sun_Exposed_Suprapubic"),
        cells=("keratinocyte",),
        distance=9_784,
        variants=(
            {
                "rsid": "rs642961",
                "pos": 209_815_925,
                "ref": "G",
                "alt": "A",
                "effect": "A abolishes an AP-2alpha binding site and lowers MCS-9.7 enhancer activity",
            },
        ),
        nearest_gene_trap="ENSG00000289700",
        citations=(
            "Rahimov et al. 2008, Nat Genet 40:1341 (rs642961 in the IRF6 enhancer MCS-9.7 disrupts"
            " AP-2alpha binding and carries the cleft lip risk)",
            "Fakhouri et al. 2014, Hum Mol Genet 23:2711 (MCS-9.7 drives oral periderm expression in"
            " transgenic mouse and its loss gives the Irf6 phenotype)",
        ),
        answer_from=("transgenic mouse reporters", "electrophoretic mobility shift", "human pedigrees"),
        note=(
            "the nearest coding TSS is not IRF6 but a novel GENCODE gene 37 bp closer, which is the"
            " trap: a heuristic that answers `the nearest gene` is at the mercy of the annotation's"
            " novel entries as much as of the biology"
        ),
    ),
    Expect(
        locus="MYB_HMIP2",
        chrom="chr6",
        element=(135_097_000, 135_098_000),
        window=(134_597_000, 135_597_000),
        classes=("program",),
        targets=("MYB",),
        direction="activates",
        tissues=("erythroid", "bone marrow"),
        gtex_tissues=("Whole_Blood",),
        cells=("K562",),
        distance=83_807,
        variants=(
            {
                "rsid": "rs66650371",
                "pos": 135_097_495,
                "ref": "TACTA",
                "alt": "TA",
                "effect": "the 3 bp deletion removes a TAL1/GATA1 site, lowers MYB and raises fetal"
                " haemoglobin",
            },
        ),
        nearest_gene_trap="HBS1L",
        citations=(
            "Farrell et al. 2011, Mol Cell Biol 31:3298 (the HBS1L-MYB intergenic elements are"
            " erythroid enhancers of MYB)",
            "Stadhouders et al. 2014, J Clin Invest 124:1699 (deleting the -84 kb element in erythroid"
            " cells lowers MYB; rs66650371 disrupts a TAL1 site and is the causal allele)",
        ),
        answer_from=("enhancer deletion in erythroid cells", "chromatin conformation", "human pedigrees"),
        note=(
            "the nearest coding TSS is HBS1L at 43 kb and the published target is MYB at 84 kb: a real"
            " nearest-gene trap at a distance the model can reach, which the panel's two megabase"
            " cases could not offer"
        ),
    ),
    Expect(
        locus="SORT1_1p13",
        chrom="chr1",
        element=(109_274_718, 109_275_218),
        window=(108_774_000, 109_774_000),
        classes=("program",),
        targets=("SORT1",),
        direction="activates",
        tissues=("liver", "hepatocyte"),
        gtex_tissues=("Liver",),
        cells=("HepG2", "hepatocyte"),
        distance=122_949,
        variants=(
            {
                "rsid": "rs12740374",
                "pos": 109_274_968,
                "ref": "G",
                "alt": "T",
                "effect": "T creates a C/EBP binding site, raising hepatic SORT1 and lowering LDL",
            },
        ),
        nearest_gene_trap="PSRC1",
        citations=(
            "Musunuru et al. 2010, Nature 466:714 (rs12740374 creates a C/EBP site in a hepatic"
            " enhancer; altering Sort1 in mouse liver by AAV changes LDL, so SORT1 is the causal gene)",
            "Kathiresan et al. 2008, Nat Genet 40:189 (the 1p13 LDL locus)",
        ),
        answer_from=("AAV expression in mouse liver", "reporter assays", "association"),
        note=(
            "the published paper reports coordinated changes at SORT1, CELSR2 and PSRC1 and concludes"
            " SORT1 is causal for the phenotype, so the expectation names SORT1 alone and PSRC1 -"
            " nearest at 8 kb, and the gene an AlphaGenome variant scan in this repository already"
            " names first - is written down as the trap before scoring"
        ),
    ),
    Expect(
        locus="IL2RA_CaRE4",
        chrom="chr10",
        element=(6_052_484, 6_052_984),
        window=(5_552_000, 6_552_000),
        classes=("program",),
        targets=("IL2RA",),
        direction="activates",
        tissues=("stimulated T cell", "T lymphocyte"),
        gtex_tissues=("Cells_EBV-transformed_lymphocytes", "Whole_Blood"),
        cells=(),
        distance=9_632,
        variants=(
            {
                "rsid": "rs61839660",
                "pos": 6_052_734,
                "ref": "C",
                "alt": "T",
                "effect": "the autoimmune-protective allele in CaRE4 lowers IL2RA induction in"
                " stimulated T cells",
            },
        ),
        nearest_gene_trap=None,
        citations=(
            "Simeonov et al. 2017, Nature 549:111 (CRISPR activation tiling finds CaRE4 in IL2RA"
            " intron 1; deleting it lowers IL2RA, and rs61839660 sits in it)",
        ),
        answer_from=("CRISPR activation screen", "CRISPR deletion of the element", "association"),
        note=(
            "the cleanest CRISPR-perturbation link in the set, and the published cell type is a"
            " stimulated primary T cell, which none of the eleven reader cells is: the cell axis is"
            " expected to be judged by GTEx alone here"
        ),
    ),
    Expect(
        locus="GDF5_GROW1",
        chrom="chr20",
        element=(35_319_108, 35_319_608),
        window=(34_819_000, 35_819_000),
        classes=("program",),
        targets=("GDF5",),
        direction="activates",
        tissues=("growth plate", "chondrocyte", "developing joint"),
        gtex_tissues=(),
        cells=(),
        distance=118_869,
        variants=(
            {
                "rsid": "rs6060369",
                "pos": 35_319_358,
                "ref": "T",
                "alt": "C",
                "effect": "the risk allele lowers GROW1 enhancer activity in the growth plate and"
                " changes GDF5 expression, knee shape and height",
            },
        ),
        nearest_gene_trap="FAM83C",
        citations=(
            "Capellini et al. 2017, Nat Genet 49:1202 (GROW1 is a growth-plate enhancer of GDF5"
            " carrying rs6060369; a humanised mouse allele changes Gdf5 expression and joint shape)",
        ),
        answer_from=("humanised mouse allele", "transgenic reporter", "association"),
        note=(
            "a 119 kb reach with two coding TSSs closer than the target: the same shape as MYB, in a"
            " tissue neither the reader nor GTEx carries"
        ),
    ),
    Expect(
        locus="PTF1A_enhancer",
        chrom="chr10",
        element=(23_219_045, 23_219_445),
        window=(22_719_000, 23_719_000),
        classes=("program",),
        targets=("PTF1A",),
        direction="activates",
        tissues=("pancreatic progenitor", "embryonic pancreas"),
        gtex_tissues=("Pancreas",),
        cells=(),
        distance=26_934,
        variants=(
            {
                "rsid": None,
                "pos": None,
                "effect": "six recessive point mutations in the 400 bp element abolish enhancer"
                " activity and cause isolated pancreatic agenesis",
            },
        ),
        nearest_gene_trap=None,
        element_source="stated",
        element_citation=(
            "Weedon et al. 2014, Nat Genet 46:61: the element is an uncharacterised ~400 bp sequence"
            " 25 kb downstream of PTF1A. No ENCODE cCRE, VISTA element or lentiMPRA tile covers it,"
            " so the interval is stated at that published offset from GENCODE's PTF1A gene end"
            " (chr10:23,194,245), the mechanism section 1 uses for SOX9 and H19"
        ),
        citations=(
            "Weedon et al. 2014, Nat Genet 46:61 (recessive mutations in a distal PTF1A enhancer cause"
            " isolated pancreatic agenesis and abolish enhancer activity)",
        ),
        answer_from=("human pedigrees", "reporter assays", "pancreatic progenitor epigenomes"),
        note=(
            "the only candidate no annotation covers, so the only one that costs a request. Its"
            " coordinate is ANCHORED on a published offset, not quoted: the anchor lands within 500 bp"
            " of an Ensembl assembly-map liftover of the published mutation positions, and the element"
            " is widened to absorb the difference"
        ),
    ),
    Expect(
        locus="PAX6_SIMO",
        chrom="chr11",
        element=(31_664_197, 31_664_597),
        window=(31_164_000, 32_164_000),
        classes=("program",),
        targets=("PAX6",),
        direction="activates",
        tissues=("eye", "lens", "developing forebrain"),
        gtex_tissues=(),
        cells=(),
        distance=146_924,
        variants=(
            {
                "rsid": None,
                "pos": None,
                "effect": "a de novo point mutation in SIMO destroys a PAX6 autoregulatory binding"
                " site and causes aniridia",
            },
        ),
        nearest_gene_trap=None,
        citations=(
            "Bhatia et al. 2013, Am J Hum Genet 93:1126 (a de novo variant in the ultraconserved SIMO"
            " element, 150 kb from PAX6 inside ELP4 intron 9, breaks PAX6 autoregulation and causes"
            " aniridia)",
            "Kleinjan et al. 2001, Hum Mol Genet 10:2049 (downstream deletions remove PAX6 regulatory"
            " elements while leaving the gene intact)",
        ),
        answer_from=("transgenic mouse and zebrafish reporters", "human deletions", "human pedigrees"),
        note=(
            "the element sits inside ELP4 and the nearest coding TSS is nonetheless PAX6, 147 kb away,"
            " because ELP4's start is further still: a long reach with no trap, which is unusual and"
            " is the reason it is kept. Coordinate from an Ensembl assembly-map liftover of the"
            " published GRCh37 position chr11:31,685,945"
        ),
    ),
    Expect(
        locus="HBA_HS40",
        chrom="chr16",
        element=(113_381, 113_908),
        window=(0, 613_000),
        classes=("program",),
        targets=("HBA1", "HBA2"),
        direction="activates",
        tissues=("erythroid", "bone marrow"),
        gtex_tissues=("Whole_Blood",),
        cells=("K562",),
        distance=63_035,
        variants=(
            {
                "rsid": None,
                "pos": None,
                "effect": "deletions removing MCS-R2 leave the alpha-globin genes intact and silence"
                " them, causing alpha-thalassaemia",
            },
        ),
        nearest_gene_trap="NPRL3",
        citations=(
            "Higgs et al. 1990, Genes Dev 4:1588 (a deletion removing HS-40 abolishes alpha-globin"
            " expression with the genes intact)",
            "Hay et al. 2016, Nat Genet 48:895 (deleting the individual MCS-R elements in mouse shows"
            " MCS-R2 carries nearly all the enhancer activity)",
            "Coelho et al. 2010, Blood 116:e46 and Sollaino et al. 2010 (human MCS-R2 deletions)",
        ),
        answer_from=("human deletions", "mouse enhancer deletion", "DNase hypersensitivity"),
        note=(
            "the alpha-globin answer to the panel's beta-globin locus control region, and the same"
            " shape as the ZRS at a twentieth of the distance: the enhancer sits inside NPRL3, which"
            " is the nearest coding TSS and the trap, and the published targets are two globin genes"
            " 60 kb away. Coordinate from an Ensembl assembly-map liftover of the published GRCh37"
            " interval chr16:163,380-163,907"
        ),
    ),
    # ---------------------------------------------------------------- exhibits, never graded
    Expect(
        locus="POU3F4_DFN3",
        chrom="chrX",
        element=(82_607_789, 82_608_789),
        window=(82_100_000, 84_100_000),
        classes=("program",),
        targets=("POU3F4",),
        direction="activates",
        tissues=("otic vesicle", "inner ear mesenchyme"),
        gtex_tissues=(),
        cells=(),
        distance=900_000,
        variants=(
            {
                "rsid": None,
                "pos": None,
                "effect": "microdeletions in a hotspot 900 kb proximal to POU3F4 cause X-linked"
                " deafness with the gene itself intact",
            },
        ),
        nearest_gene_trap=None,
        element_source="stated",
        element_citation=(
            "de Kok et al. 1996, Hum Mol Genet 5:1229: a microdeletion hotspot 900 kb proximal to"
            " POU3F4. The publication gives an offset and no citable coordinate, so this interval is"
            " anchored at that offset from GENCODE's POU3F4 start and is an EXHIBIT ONLY: it is"
            " never graded, and the only question asked of it is the reach one, whose answer a 900 kb"
            " offset settles whatever the placement"
        ),
        citations=(
            "de Kok et al. 1996, Hum Mol Genet 5:1229 (the DFN3 microdeletion hotspot 900 kb proximal"
            " to POU3F4)",
            "Naranjo et al. 2010, Development 137:1721 (otic enhancers in the 1 Mb upstream region)",
        ),
        answer_from=("human microdeletions", "transgenic mouse reporters"),
        note="exhibit: a reach fatality by construction, kept so the filter is counted against"
        " something rather than asserted",
    ),
    Expect(
        locus="MYC_BENC",
        chrom="chr8",
        element=(129_435_230, 129_437_230),
        window=(127_500_000, 130_000_000),
        classes=("program",),
        targets=("MYC",),
        direction="activates",
        tissues=("haematopoietic stem cell", "erythroid"),
        gtex_tissues=("Whole_Blood",),
        cells=("K562",),
        distance=1_700_000,
        variants=(
            {
                "rsid": None,
                "pos": None,
                "effect": "deleting individual modules of the blood enhancer cluster lowers MYC and"
                " collapses the haematopoietic hierarchy",
            },
        ),
        nearest_gene_trap="GSDMC",
        element_source="stated",
        element_citation=(
            "Bahr et al. 2018, Nature 553:515 and Fulco et al. 2016, Science 354:769: a blood enhancer"
            " cluster about 1.7 Mb downstream of MYC, dissected by CRISPR interference and by module"
            " deletion in mouse. The publications give an offset rather than a coordinate this module"
            " can cite, so the interval is anchored at that offset from GENCODE's MYC start and is an"
            " EXHIBIT ONLY, never graded"
        ),
        citations=(
            "Fulco et al. 2016, Science 354:769 (CRISPR interference tiling finds seven MYC enhancers"
            " in K562, up to 1.9 Mb downstream)",
            "Bahr et al. 2018, Nature 553:515 (the blood enhancer cluster at +1.7 Mb; deleting modules"
            " lowers Myc)",
        ),
        answer_from=("CRISPR interference tiling", "mouse enhancer deletion"),
        note="exhibit: the second reach fatality, and the one that shows the filter is not only about"
        " rare developmental loci - a CRISPRi screen in the reader's own cell type produces them too",
    ),
)

#: candidates whose element is placed on a published offset with no citable coordinate. Section 13 of
#: docs/LOCI-BENCHMARK.md withdrew a case for exactly this, so these are never graded; the reach
#: question is the only one an approximate placement can still answer, and it is answered.
EXHIBITS_ONLY: frozenset[str] = frozenset({"POU3F4_DFN3", "MYC_BENC"})

#: the chr21 case the reach filter must keep. A TSS-distance filter would drop RUNX1, whose body is
#: 1.22 Mb long, and would look right doing it. Pinned in the tests.
BODY_NOT_TSS_CASE = {
    "gene": "RUNX1",
    "chrom": "chr21",
    "why": "1.22 Mb of gene body, so elements 524 kb or more from its TSS still sit inside it",
    "measured": "79 of 5,174 chr21 elements with a coding target name a gene whose TSS is 524 kb or"
    " more away, up to 986 kb, and every one is a gene whose body overlaps the window",
}


def by_name() -> dict[str, Expect]:
    return {e.locus: e for e in CANDIDATES}


def plan(results_dir: Path = RESULTS_DIR) -> list[dict[str, Any]]:
    """Which candidates the model could answer, and what each would cost. No request, no network.

    Two free questions decide the budget, and in this order: can the model see the target at all
    (`loci.read_reach`), and did anybody already draw an interval over the published element (if so
    the finished sweep has already deleted it and the answer costs nothing).
    """
    rows: list[dict[str, Any]] = []
    for chrom in sorted({e.chrom for e in CANDIDATES}):
        ch = loci.Chromosome(chrom, results_dir)
        try:
            for e in (x for x in CANDIDATES if x.chrom == chrom):
                reach = loci.read_reach(ch, e)
                s, en = e.element
                scored = loci._deletion_rows(chrom, s, en, results_dir)  # noqa: SLF001
                rows.append(
                    {
                        "locus": e.locus,
                        "chrom": chrom,
                        "element": [s, en],
                        "targets": list(e.targets),
                        "askable": reach["askable"],
                        "targets_in_reach": reach["targets_in_reach"],
                        "targets_out_of_reach": reach["targets_out_of_reach"],
                        "exhibit_only": e.locus in EXHIBITS_ONLY,
                        "graded": reach["askable"] and e.locus not in EXHIBITS_ONLY,
                        "ccres_over_element": len(ch.ccres_in(s, en)),
                        "already_scored_over_element": len(scored),
                        "element_source": e.element_source,
                        "requests": 0 if (not reach["askable"] or e.locus in EXHIBITS_ONLY or scored) else 1,
                    }
                )
        finally:
            ch.close()
    return sorted(rows, key=lambda r: r["locus"])


def graded(rows: list[dict[str, Any]]) -> tuple[Expect, ...]:
    """The candidates the run grades: in reach, and placed precisely enough to grade."""
    keep = {r["locus"] for r in rows if r["graded"]}
    return tuple(e for e in CANDIDATES if e.locus in keep)


def reach_fatalities(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """How many candidates the free filter killed, and why - a number in its own right."""
    dead = [r for r in rows if not r["askable"]]
    return {
        "candidates": len(rows),
        "died_at_the_reach_filter": len(dead),
        "loci": [
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
        "reading": (
            "a locus outside the scorer's 1 Mb input is unaskable, not a miss. Both fatalities here"
            " are exhibits placed on an offset, so the graded set loses nothing; the number worth"
            " carrying is the other one - 9 of 11 published enhancer-gene links sit well inside the"
            " model's reach, so the panel's three megabase cases are not what a published enhancer"
            " usually looks like"
        ),
    }


def compare_with_controls(result: dict[str, Any]) -> dict[str, Any]:
    """The four claims at the candidates against their matched windows, standardised.

    `genomeos.compare` does the standardising; this only shapes the rows. The covariates are the ones
    the windows were matched on, so the standardisation is a check that the matching worked rather
    than a repair of it.
    """
    from genomeos.compare import Strata, imbalance, standardised

    def row(gc, dist, con, claims) -> dict[str, Any] | None:
        if gc is None or dist is None:
            return None
        return {"gc": gc, "log_distance": (dist or 1) ** 0.5, "constrained": con or 0.0, **claims}

    targets = [
        row(
            r.get("gc"),
            r.get("distance_to_coding_tss"),
            r.get("constrained_fraction"),
            loci.claims(r["readings"]),
        )
        for r in result["loci"]
    ]
    controls = [
        row(w.get("gc"), w.get("distance_to_coding_tss"), w.get("constrained_fraction"), w["claims"])
        for w in result["negatives"]
    ]
    targets = [r for r in targets if r]
    controls = [r for r in controls if r]
    strata = Strata(gc=(0.40, 0.48, 0.56), log_distance=(100.0, 250.0, 500.0), constrained=(0.02, 0.10))
    # the same strata plus one: whether an element inside the window has actually been deleted.
    # docs/LOCI-BENCHMARK.md section 8 withdrew the direction claim because it was a statement about
    # where the requests had been spent, and standardising on GC, distance and constraint does not
    # touch that - a window nobody scored cannot produce a direction at all, whatever its GC is.
    covered = Strata(
        gc=(0.40, 0.48, 0.56),
        log_distance=(100.0, 250.0, 500.0),
        constrained=(0.02, 0.10),
        deletion_scored=(0.5,),
    )
    for r in targets + controls:
        r["deletion_scored"] = float(bool(r["deletion_scored"]))
    out: dict[str, Any] = {
        "reading": "direct standardisation on the covariates the windows were matched on; `raw` is"
        " what an unmatched claim would have said, `matched` holds the covariates fixed",
        "imbalance": imbalance(targets, controls, strata),
    }
    for claim in ("target", "cell", "direction", "storage"):
        out[claim] = standardised(targets, controls, strata, hit=claim)
        out[f"{claim}_given_coverage"] = standardised(targets, controls, covered, hit=claim)
    out["why_given_coverage"] = (
        "the four covariates the windows were matched on say nothing about whether anybody spent a"
        " request on them, so standardising on those alone reproduces the coverage artefact section 8"
        " of docs/LOCI-BENCHMARK.md retracted. `*_given_coverage` adds `has an element inside it been"
        " deleted` as a fifth stratum, and the direction claim is the one to read there"
    )
    return out


def beside_the_panel(result: dict[str, Any], results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """The candidates' headline rates beside the panel's, with both denominators and no third one."""
    panel = (load_result("loci_benchmark", results_dir) or {}).get("aggregate") or {}
    mine = result["aggregate"]

    def pair(key: str) -> dict[str, Any]:
        return {
            "candidates": {k: mine[key][k] for k in ("k", "n", "rate")},
            "panel": {k: (panel.get(key) or {}).get(k) for k in ("k", "n", "rate")},
        }

    return {
        "target_derived": pair("target_derived"),
        "target_derived_where_the_model_could_answer": pair("target_derived_where_the_model_could_answer"),
        "target_heuristic": pair("target_heuristic"),
        "target_looked_up": pair("target_looked_up"),
        "chance_floor": {
            "candidates": mine["target_by_chance"],
            "panel": panel.get("target_by_chance"),
        },
        "reading": (
            "the two denominators the panel already reports, and nothing else. The candidates are a"
            " separate panel and their rate is not pooled with the panel's: the two sets were chosen"
            " by different people on different days for different reasons, and a pooled rate would"
            " hide that this set is short-range and the panel's is not"
        ),
    }


def run(
    results_dir: Path = RESULTS_DIR,
    network: bool = True,
    negatives: bool = True,
    gtex_dir: Path = GTEX_DIR,
    progress=None,
) -> dict[str, Any]:
    """Plan, filter by reach, then read and score the graded candidates through loci.build."""
    rows = plan(results_dir)
    panel = graded(rows)
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
    out["body_not_tss"] = BODY_NOT_TSS_CASE
    out["exhibits_only"] = sorted(EXHIBITS_ONLY)
    out["against_controls"] = compare_with_controls(out)
    out["beside_the_panel"] = beside_the_panel(out, results_dir)
    out["note"] = (
        "A second, separately registered panel of published enhancer-gene loci, graded by"
        " genomeos.benchmark.loci with no change to any scorer or hit rule. Reported beside"
        " docs/LOCI-BENCHMARK.md's seventeen, never pooled with them. " + loci.NOTE
    )
    save_result(NAME, out, results_dir)
    return out
