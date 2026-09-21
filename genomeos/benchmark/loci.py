# SPDX-License-Identifier: AGPL-3.0-or-later
"""The known-locus benchmark: the machinery against loci whose answer is published.

Albert's brief (2026-09-13): take the few places where the biology is already known -
their variants, their syntax, their enhancers, their context - and use them as a
coherence check on everything else. If the pipeline cannot reproduce what is settled at
a dozen well-studied loci, its genome-wide numbers are not to be believed.

Each locus carries a machine-readable *expectation* with citations: the target gene, the
tissue or cell type, the causal variant where there is one, the direction of effect, the
distance from element to target, and which of Albert's three functional classes it
belongs to (fixed syntax; a stored value, a position holding one of a few recurring
alleles like a database column; a program that executes according to those values). A
fourth property is asked of two loci only, *order*: the beta-globin cluster switches
fetal to adult under one locus control region, and the HOXD cluster is read in the order
its genes sit on the chromosome, so between them the panel covers ordered-in-time and
ordered-in-space.

Then the project's own layers are read over the same coordinates and compared. The main
threat is circularity, and handling it is the difference between a benchmark and a
mirror: several of these answers are already inside the annotations we read, so every
hit is labelled by where it came from.

  derived    a genome-wide measurement or model read blind at this locus: the AlphaGenome
             deletion effect, GTEx eQTLs, ENCODE DNase (the reader), lentiMPRA activity,
             Zoonomia phyloP, gnomAD Gnocchi, the imported people's own variants, 1000
             Genomes frequencies. None of these names a function; the reading does.
  targeted   a measurement made *at this locus because it was already interesting*:
             Kircher et al.'s saturation mutagenesis. It can confirm a causal base, never
             discover one, so it is counted apart.
  heuristic  the node model's nearest coding transcription start: geometry over GENCODE,
             the baseline the rest has to beat.
  looked_up  an annotation that already contains the answer: VISTA's tissue for the
             element, the GWAS Catalog's mapped gene, ClinVar's gene, the registry's
             element class. A hit here proves nothing about the approach and is reported
             on its own line.

Negative controls are required, not optional: for every locus, windows of the same
length, GC, distance to the nearest coding transcription start and mammalian constraint,
with no known function, scored through exactly the same readers. If the pipeline names a
target, a cell and a direction at the negatives as readily as at the positives, the
panel's hit rate means nothing, and that is the finding.

Nothing here calls the AlphaGenome API: a scoring chain owns that quota. Model readings
come from what is already computed (`data/knowledge/alphagenome/all_elements/*.json`, the
per-element cache, the committed sampled runs); where a locus needs a reading nobody has
computed, it is recorded as `pending` with the reason and the cost, never skipped in
silence. Design and results: docs/LOCI-BENCHMARK.md.
"""

from __future__ import annotations

import bisect
import gzip
import json
import random
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from genomeos.results import RESULTS_DIR, load_result, save_result

# ---------------------------------------------------------------------------- vocabulary
CLASSES = ("syntax", "storage", "program", "order")
PROVENANCE = ("derived", "targeted", "heuristic", "looked_up")
#: the eleven ENCODE cell types the reader has peaks for, genome-wide (genome/reader.py)
READER_CELLS = (
    "K562",
    "HepG2",
    "GM12878",
    "H1",
    "IMR-90",
    "SK-N-SH",
    "cardiac muscle cell",
    "keratinocyte",
    "hepatocyte",
    "astrocyte",
    "CD14-positive monocyte",
)
#: the deletion scorer's whole input, resized around the element (dna_client.SEQUENCE_LENGTH_1MB).
#: A published target whose gene body falls outside it cannot be named however the element is scored.
MODEL_WINDOW = 1_048_576
GTEX_DIR = Path("data/knowledge/loci_benchmark")
#: results holding deletions of intervals a panel STATED from a publication because no registry drew
#: one. Read first by `_deletion_rows` and labelled `stated_interval`, so a hit on one is never
#: pooled with a hit on somebody else's annotation. One per panel: the seventeen, the candidates, the
#: third set, then the fourth frame. A panel missing from this tuple pays for a deletion that is
#: written, never read back, and graded as though no request had been made.
STATED_INTERVAL_RESULTS = (
    "loci_stated_intervals",
    "loci_candidate_intervals",
    "loci_third_intervals",
    "loci_fourth_intervals",
)
NEGATIVES_PER_LOCUS = 5
CANDIDATES_PER_LOCUS = 40  # sampled, then narrowed to the best matches on four covariates
MATCH_GC = 0.04  # matched window: GC within this of the positive's
MATCH_DISTANCE = 0.35  # and distance to the nearest coding TSS within this fraction
NEGATIVE_KEEP_OUT = 200_000  # a negative window stays this far from any panel locus
K562_STRONG_PEAK = 100.0  # DNase signalValue of an unambiguous erythroid hypersensitive site
COMMON_AF = 0.05  # a common value: gnomAD allele frequency at or above this
#: value-domain size from common variable positions per kb, on a decade scale. Set AFTER reading
#: two background windows (1 to 2 per kb) and the storage loci, so the value-domain agreement is
#: descriptive, not a test; the matched negatives are what test it.
DOMAIN_FEW_PER_KB = 10.0
DOMAIN_MANY_PER_KB = 50.0
EVIDENCE = {
    "expectation": "curated: the published answer per locus, with citations, written down before the run",
    "derived": (
        "the project's own layers read at the locus: AlphaGenome deletion effects (already computed),"
        " GTEx v8 eQTLs, ENCODE DNase (the reader), ENCODE4 lentiMPRA, Zoonomia phyloP, gnomAD Gnocchi,"
        " the imported people's variants, 1000 Genomes frequencies"
    ),
    "targeted": "experimental: Kircher et al. 2019 saturation mutagenesis, measured at these loci",
    "heuristic": "inferred: nearest coding TSS inside the CTCF node (confidence 0.4)",
    "looked_up": "curated: VISTA tissues, GWAS Catalog mapped gene, ClinVar gene, the ENCODE cCRE class",
}


@dataclass(frozen=True)
class Expect:
    """What the literature says about one locus, written down before anything is read."""

    locus: str
    chrom: str
    element: tuple[int, int]  # the element or variant-carrying window, 0-based half-open
    window: tuple[int, int]  # the locus: the space a target may sit in
    classes: tuple[str, ...]
    targets: tuple[str, ...]  # the published target gene(s)
    direction: str  # activates | represses | coding | none
    tissues: tuple[str, ...]  # the published tissue or cell type, in words
    gtex_tissues: tuple[str, ...] = ()  # GTEx v8 tissue names that count as the right tissue
    cells: tuple[str, ...] = ()  # of READER_CELLS, the ones that count; () = the panel has none
    distance: int | None = None  # published element-to-target distance, bp
    variants: tuple[dict[str, Any], ...] = ()  # rsid, pos (1-based), effect
    nearest_gene_trap: str | None = None  # the gene the nearest-TSS heuristic is expected to name wrongly
    value_domain: str | None = None  # two | few | many, for the storage loci
    frequency: str | None = None  # the published population structure, in words
    #: published mechanisms that are not in the reference sequence at all, so no layer here can carry
    #: them: parent of origin, copy number, developmental stage. Read as pending, never as a miss.
    beyond_sequence: tuple[str, ...] = ()
    #: "annotated" when some registry drew an interval over the published element, "stated" when the
    #: panel draws it itself because nobody else did. A stated interval must carry its own citation,
    #: and every rate computed from one is reported apart from the registry-derived ones.
    element_source: str = "annotated"
    element_citation: str | None = None
    citations: tuple[str, ...] = ()
    answer_from: tuple[str, ...] = ()  # the kinds of data the published answer itself came from
    note: str = ""

    def as_dict(self) -> dict[str, Any]:
        return {
            "locus": self.locus,
            "chrom": self.chrom,
            "element": list(self.element),
            "window": list(self.window),
            "length": self.element[1] - self.element[0],
            "classes": list(self.classes),
            "targets": list(self.targets),
            "direction": self.direction,
            "tissues": list(self.tissues),
            "gtex_tissues": list(self.gtex_tissues),
            "cells": list(self.cells),
            "distance": self.distance,
            "variants": [dict(v) for v in self.variants],
            "nearest_gene_trap": self.nearest_gene_trap,
            "value_domain": self.value_domain,
            "frequency": self.frequency,
            "beyond_sequence": list(self.beyond_sequence),
            "element_source": self.element_source,
            "element_citation": self.element_citation,
            "citations": list(self.citations),
            "answer_from": list(self.answer_from),
            "note": self.note,
        }


# ------------------------------------------------------------------------------- the panel
PANEL: tuple[Expect, ...] = (
    Expect(
        locus="SHH_ZRS",
        chrom="chr7",
        element=(156_791_086, 156_791_875),
        window=(155_700_000, 157_000_000),
        classes=("program",),
        targets=("SHH",),
        direction="activates",
        tissues=("limb bud", "zone of polarising activity"),
        gtex_tissues=(),
        cells=(),
        distance=979_000,
        variants=(
            {
                "rsid": None,
                "pos": None,
                "effect": "point mutations across the ZRS cause preaxial polydactyly (many, mostly private)",
            },
        ),
        nearest_gene_trap="LMBR1",
        citations=(
            "Lettice et al. 2003, Hum Mol Genet 12:1725 (the ZRS, 1 Mb from SHH, point mutations"
            " cause preaxial polydactyly)",
            "Sagai et al. 2005, Development 132:797 (deleting the ZRS removes limb Shh entirely)",
            "VISTA hs2496, positive in limb",
        ),
        answer_from=("transgenic mouse reporters", "targeted deletion in mouse", "human pedigrees"),
        note=(
            "the long-range case: the element sits inside LMBR1 intron 5, so the nearest-TSS heuristic"
            " must fail and the node plus deletion evidence must carry it"
        ),
    ),
    Expect(
        locus="HERC2_OCA2",
        chrom="chr15",
        element=(28_120_106, 28_121_063),
        window=(27_700_000, 28_400_000),
        classes=("storage", "program"),
        targets=("OCA2",),
        direction="activates",
        tissues=("melanocyte", "skin"),
        gtex_tissues=("Skin_Sun_Exposed_Lower_leg", "Skin_Not_Sun_Exposed_Suprapubic"),
        cells=("keratinocyte",),
        distance=21_000,
        variants=(
            {
                "rsid": "rs12913832",
                "pos": 28_120_472,
                "ref": "A",
                "alt": "G",
                "effect": "G lowers enhancer activity and OCA2 expression: blue eyes; A: brown",
            },
        ),
        nearest_gene_trap="HERC2",
        value_domain="two",
        frequency="G near 0.8 in Europeans, rare in Africa and East Asia",
        citations=(
            "Sturm et al. 2008, Am J Hum Genet 82:424 (rs12913832 in HERC2 intron 86 sets blue-brown)",
            "Eiberg et al. 2008, Hum Genet 123:177",
            "Visser et al. 2012, Genome Res 22:446 (rs12913832 attenuates the enhancer-OCA2 promoter loop)",
        ),
        answer_from=("reporter assays", "chromatin conformation", "association in pedigrees"),
        note=(
            "the storage case: one base holding one of two values inside a program. Gnocchi does not score"
            " this kilobase at all, so the human axis is blind here"
        ),
    ),
    Expect(
        locus="MCM6_LCT",
        chrom="chr2",
        element=(135_850_876, 135_851_276),
        window=(135_700_000, 136_000_000),
        classes=("storage", "program"),
        targets=("LCT",),
        direction="activates",
        tissues=("small intestine", "intestinal epithelium"),
        gtex_tissues=("Small_Intestine_Terminal_Ileum", "Colon_Transverse", "Colon_Sigmoid"),
        cells=(),
        distance=13_900,
        variants=(
            {
                "rsid": "rs4988235",
                "pos": 135_851_076,
                "ref": "G",
                "alt": "A",
                "effect": "A (T-13910 on the LCT strand) keeps LCT transcribed into adulthood",
            },
        ),
        nearest_gene_trap="MCM6",
        value_domain="few",
        frequency="A near 0.75 in northern Europe, near 0 in East Asia; other alleles in Africa and Arabia",
        citations=(
            "Enattah et al. 2002, Nat Genet 30:233 (C/T-13910 in MCM6 intron 13 with lactase persistence)",
            "Olds and Sibley 2003, Hum Mol Genet 12:2333; Lewinsky et al. 2005, Hum Mol Genet 14:3945"
            " (the T allele raises LCT promoter activity in intestinal cells)",
            "Tishkoff et al. 2007, Nat Genet 39:31 (convergent adaptation, different African alleles)",
        ),
        answer_from=("reporter assays", "expression in biopsies", "association and selection scans"),
        note="allele-dependent activity with the sharpest population structure in the panel",
    ),
    Expect(
        locus="HBB_LCR",
        chrom="chr11",
        element=(5_269_000, 5_295_000),
        window=(5_200_000, 5_320_000),
        classes=("syntax", "order", "program"),
        targets=("HBB", "HBD", "HBG1", "HBG2", "HBE1"),
        direction="activates",
        tissues=("erythroid", "bone marrow", "fetal liver"),
        gtex_tissues=("Whole_Blood",),
        cells=("K562",),
        distance=None,
        variants=(
            {
                "rsid": None,
                "pos": None,
                "effect": "deleting the LCR (Hispanic thalassaemia) silences every beta-like globin gene",
            },
        ),
        beyond_sequence=("developmental stage",),
        citations=(
            "Tuan et al. 1985, PNAS 82:6384 (the erythroid-specific hypersensitive sites 5' of epsilon)",
            "Grosveld et al. 1987, Cell 51:975 (position-independent high-level expression from the LCR)",
            "Forrester et al. 1990, Genes Dev 4:1637 (the Hispanic thalassaemia deletion removes the LCR)",
            "Sankaran and Orkin 2013, Cold Spring Harb Perspect Med 3:a011643 (the fetal-to-adult switch)",
        ),
        answer_from=("DNase hypersensitivity", "transgenic mice", "human deletions"),
        note=(
            "five genes in one node under one control region, read in order through development:"
            " HBE1 and HBG1/2 in the fetus, HBB and HBD after birth. The expectation is five"
            " erythroid-specific hypersensitive sites 3 to 5 kb apart and the globin genes in one node"
        ),
    ),
    Expect(
        locus="FTO_IRX3",
        chrom="chr16",
        element=(53_766_842, 53_767_242),
        window=(53_600_000, 54_400_000),
        classes=("program",),
        targets=("IRX3", "IRX5"),
        direction="represses",
        tissues=("adipocyte precursor", "adipose", "brain"),
        gtex_tissues=("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Brain_Cortex"),
        cells=(),
        distance=520_000,
        variants=(
            {
                "rsid": "rs1421085",
                "pos": 53_767_042,
                "ref": "T",
                "alt": "C",
                "effect": "C disrupts an ARID5B repressor site and doubles IRX3 and IRX5 in adipocyte"
                " precursors, so the element represses through the ancestral allele",
            },
            {
                "rsid": "rs1558902",
                "pos": 53_769_662,
                "effect": "the lead obesity association in FTO intron 1",
            },
        ),
        nearest_gene_trap="FTO",
        citations=(
            "Smemo et al. 2014, Nature 507:371 (the FTO intron-1 variants connect to IRX3, not FTO)",
            "Claussnitzer et al. 2015, N Engl J Med 373:895 (rs1421085 disrupts ARID5B; IRX3/IRX5 double)",
        ),
        answer_from=("chromatin conformation", "eQTL in brain and adipose", "editing in primary cells"),
        note=(
            "the textbook nearest-gene fallacy, and our own data says the nearest gene is wrong about a"
            " third of the time, so this locus has to be right for the node model to be worth anything"
        ),
    ),
    Expect(
        locus="BCL11A_enhancer",
        chrom="chr2",
        element=(60_490_000, 60_496_000),
        window=(60_200_000, 60_800_000),
        classes=("program", "storage"),
        targets=("BCL11A",),
        direction="activates",
        tissues=("erythroid",),
        gtex_tissues=("Whole_Blood",),
        cells=("K562",),
        distance=59_000,
        variants=(
            {
                "rsid": "rs1427407",
                "pos": 60_490_908,
                "ref": "T",
                "alt": "G",
                "effect": "T disrupts a GATA1 site in the +62 DHS, lowers BCL11A and raises fetal"
                " haemoglobin",
            },
        ),
        citations=(
            "Bauer et al. 2013, Science 342:253 (the erythroid enhancer of BCL11A in intron 2;"
            " rs1427407 is the causal variant of the fetal-haemoglobin association)",
            "Canver et al. 2015, Nature 527:192 (Cas9 tiling finds the +58 core and its GATA1 motif)",
        ),
        answer_from=("GWAS fine-mapping", "enhancer assays", "Cas9 tiling in erythroid cells"),
        note="a program with a known causal variant, in the erythroid lineage the reader holds (K562)",
    ),
    Expect(
        locus="MYC_8q24",
        chrom="chr8",
        element=(127_400_860, 127_401_260),
        window=(126_800_000, 128_200_000),
        classes=("program",),
        targets=("MYC",),
        direction="activates",
        tissues=("colon", "colorectal epithelium"),
        gtex_tissues=("Colon_Transverse", "Colon_Sigmoid"),
        cells=(),
        distance=334_000,
        variants=(
            {
                "rsid": "rs6983267",
                "pos": 127_401_060,
                "ref": "G",
                "alt": "T",
                "effect": "G strengthens a TCF7L2 site; the enhancer contacts MYC 335 kb away",
            },
            {
                "rsid": "rs11986220",
                "pos": 127_519_444,
                "effect": "a second, prostate-specific 8q24 enhancer variant (FOXA1)",
            },
        ),
        nearest_gene_trap="CASC8",
        citations=(
            "Pomerantz et al. 2009, Nat Genet 41:882 (the rs6983267 enhancer physically interacts with MYC)",
            "Tuupanen et al. 2009, Nat Genet 41:885 (rs6983267 enhances TCF7L2 binding and Wnt signalling)",
            "Ahmadiyeh et al. 2010, PNAS 107:9742 (several tissue-specific 8q24 enhancers, one target)",
        ),
        answer_from=("chromatin conformation", "GWAS fine-mapping", "reporter assays"),
        note="many tissue-specific enhancers for one gene: does a gene's whole input assemble",
    ),
    Expect(
        locus="ABO",
        chrom="chr9",
        element=(133_255_000, 133_258_000),
        window=(133_200_000, 133_300_000),
        classes=("storage",),
        targets=("ABO",),
        direction="coding",
        tissues=("red blood cell", "epithelium"),
        gtex_tissues=("Whole_Blood",),
        cells=(),
        distance=0,
        variants=(
            {
                "rsid": "rs8176719",
                "pos": 133_257_521,
                "effect": "a single-base deletion frameshifts the transferase: the O allele",
            },
            {
                "rsid": "rs8176746",
                "pos": 133_255_935,
                "protein_change": "p.Leu266Met",
                "effect": "one of the substitutions separating the B transferase from A",
            },
        ),
        value_domain="few",
        frequency="A, B and O in every population, O commonest; frequencies differ by continent",
        citations=(
            "Yamamoto et al. 1990, Nature 345:229 (the molecular basis of ABO: a frameshift for O and"
            " four substitutions for B)",
        ),
        answer_from=("cDNA sequencing", "serology"),
        note="storage with a small value domain: the control for the many-valued case below",
    ),
    Expect(
        locus="HLA_DRB1",
        chrom="chr6",
        element=(32_577_901, 32_589_848),
        window=(32_400_000, 32_800_000),
        classes=("storage",),
        targets=("HLA-DRB1",),
        direction="coding",
        tissues=("antigen-presenting cell", "B cell", "monocyte"),
        gtex_tissues=("Whole_Blood", "Cells_EBV-transformed_lymphocytes", "Spleen"),
        cells=("GM12878", "CD14-positive monocyte"),
        distance=0,
        value_domain="many",
        frequency="thousands of alleles; the polymorphism is concentrated in exon 2 (the peptide groove)",
        citations=(
            "Barker et al. 2023, Nucleic Acids Res 51:D1053 (IPD-IMGT/HLA: >3,600 DRB1 alleles)",
            "Brown et al. 1993, Nature 364:33 (the class II groove the exon-2 polymorphism lines)",
        ),
        answer_from=("allele sequencing", "structures"),
        note=(
            "storage with a very large value domain: does the machinery tell a two-value slot from a"
            " many-valued one"
        ),
    ),
    Expect(
        locus="APP",
        chrom="chr21",
        element=(25_897_420, 25_897_820),
        window=(25_800_000, 26_200_000),
        classes=("storage", "syntax"),
        targets=("APP",),
        direction="coding",
        tissues=("brain", "neuron"),
        gtex_tissues=("Brain_Cortex", "Brain_Cerebellum"),
        cells=("SK-N-SH", "astrocyte"),
        distance=0,
        variants=(
            {
                "rsid": "rs63750847",
                "pos": 25_897_620,
                "ref": "C",
                "alt": "T",
                "protein_change": "p.Ala673Thr",
                "effect": "A673T, protective against Alzheimer's disease; the same codon region as the"
                " pathogenic A673V",
            },
        ),
        value_domain="two",
        frequency="T at 0.5% in Iceland, essentially absent elsewhere",
        citations=(
            "Jonsson et al. 2012, Nature 488:96 (A673T protects against Alzheimer's disease)",
            "Di Fede et al. 2009, Science 323:1473 (A673V, the recessive pathogenic change at that codon)",
        ),
        answer_from=("whole-genome association in Iceland", "biochemistry of the peptide"),
        note="the coding-side calibration on the chromosome every layer has been run over",
    ),
    Expect(
        locus="TP53",
        chrom="chr17",
        element=(7_675_954, 7_676_354),
        window=(7_600_000, 7_760_000),
        classes=("syntax", "storage"),
        targets=("TP53",),
        direction="coding",
        tissues=("every tissue",),
        gtex_tissues=("Whole_Blood", "Lung", "Liver"),
        cells=("K562", "HepG2", "GM12878", "IMR-90"),
        distance=0,
        variants=(
            {
                "rsid": "rs1042522",
                "pos": 7_676_154,
                "ref": "G",
                "alt": "C",
                "protein_change": "p.Pro72Arg",
                "effect": "P72R, a common coding value in an otherwise invariant frame",
            },
        ),
        value_domain="two",
        frequency="the R72 allele runs from 0.2 to 0.7 depending on the population",
        citations=(
            "Olivier et al. 2010, Cold Spring Harb Perspect Biol 2:a001008 (TP53 mutation across cancers)",
            "Dumont et al. 2003, Nat Genet 33:357 (codon 72 changes the apoptotic potential)",
        ),
        answer_from=("tumour sequencing", "functional assays"),
        note="coding syntax carrying one common value: the class test for a constrained frame",
    ),
    Expect(
        locus="HOXD",
        chrom="chr2",
        element=(176_092_689, 176_190_926),
        window=(175_900_000, 176_400_000),
        classes=("order", "program", "syntax"),
        targets=("HOXD13", "HOXD12", "HOXD11", "HOXD10", "HOXD9", "HOXD8", "HOXD4", "HOXD3", "HOXD1"),
        direction="activates",
        tissues=("limb bud", "embryonic trunk"),
        gtex_tissues=(),
        cells=(),
        distance=None,
        citations=(
            "Duboule 1994, Development Suppl:135 and Deschamps and Duboule 2017, Genes Dev 31:1406"
            " (temporal and spatial collinearity: 3' genes first and anterior, 5' genes later and posterior)",
            "Andrey et al. 2013, Science 340:1234167 (a switch between two topological domains underlies"
            " HoxD collinearity in limbs)",
            "Rodriguez-Carballo et al. 2017, Genes Dev 31:2264 (a CTCF-site boundary between Hoxd11 and"
            " Hoxd12 splits the cluster's two regulatory landscapes)",
        ),
        answer_from=(
            "in-situ hybridisation over development",
            "4C and Hi-C in mouse limbs",
            "mouse deletions",
        ),
        note=(
            "the ordered case: position in the cluster is the order of execution. Our questions are whether"
            " the node model puts a boundary between HOXD11 and HOXD13, whether the cluster's non-coding"
            " content is denser and more conserved than matched windows, and whether anything we read"
            " carries the order at all"
        ),
    ),
    # ---- added 2026-09-15 (genomeos-k2), widening the panel where the first twelve were thin:
    # a second long-range case, an imprinted locus, a dosage locus, and two more value domains.
    Expect(
        locus="SOX9_PierreRobin",
        chrom="chr17",
        element=(70_669_000, 70_672_000),
        window=(70_400_000, 72_400_000),
        classes=("program",),
        targets=("SOX9",),
        direction="activates",
        tissues=("mandibular arch", "cranial neural crest"),
        gtex_tissues=(),
        cells=(),
        distance=1_450_000,
        variants=(
            {
                "rsid": None,
                "pos": None,
                "effect": "translocation breakpoints and point mutations in the element cause"
                " Pierre Robin sequence (individually rare)",
            },
        ),
        nearest_gene_trap="KCNJ2",
        element_source="stated",
        element_citation=(
            "Benko et al. 2009, Nat Genet 41:359: the Pierre Robin breakpoint cluster lies 1.0 to"
            " 1.5 Mb upstream of SOX9 and the conserved element driving mandibular expression sits at"
            " about -1.45 Mb. No ENCODE cCRE, VISTA element or lentiMPRA tile covers it, so the panel"
            " states the interval at that offset from GENCODE's canonical SOX9 TSS rather than"
            " borrowing one somebody else drew"
        ),
        citations=(
            "Benko et al. 2009, Nat Genet 41:359 (translocation breakpoints 1.0 to 1.5 Mb upstream of"
            " SOX9 cause Pierre Robin sequence; a conserved element at about -1.45 Mb drives mandibular"
            " expression)",
            "Gordon et al. 2014, J Med Genet 51:264 (point mutations in that element)",
            "Long et al. 2020, Cell Stem Cell 27:765 (the SOX9 regulatory landscape stretches over 2 Mb)",
        ),
        answer_from=("human translocation breakpoint mapping", "transgenic mouse reporters"),
        note=(
            "the second long-range case, added because one example of a megabase reach is an anecdote."
            " The element's coordinate is ANCHORED, not cited: the published offset is about -1.45 Mb"
            " and the anchor is GENCODE's canonical SOX9 TSS (chr17:72,121,019), so the element is"
            " placed rather than quoted and the window is wide enough to absorb the uncertainty."
            " The nearest coding TSS is KCNJ2, 500 kb away: the same trap the ZRS sets, in a gene desert"
        ),
    ),
    Expect(
        locus="H19_ICR1",
        chrom="chr11",
        element=(1_999_800, 2_001_900),
        window=(1_900_000, 2_250_000),
        classes=("program", "syntax"),
        targets=("IGF2", "H19"),
        direction="represses",
        tissues=("most fetal tissues", "placenta"),
        gtex_tissues=(),
        cells=(),
        distance=138_000,
        variants=(
            {
                "rsid": None,
                "pos": None,
                "effect": "microdeletions of the CTCF sites and loss of methylation cause"
                " Beckwith-Wiedemann and Silver-Russell syndromes",
            },
        ),
        nearest_gene_trap="MRPL23",
        beyond_sequence=("parent of origin",),
        element_source="stated",
        element_citation=(
            "Bell and Felsenfeld 2000, Nature 405:482 and Hark et al. 2000, Nature 405:486: the"
            " imprinting control region is the CTCF-site cluster 2 to 4 kb upstream of the H19"
            " promoter. No ENCODE cCRE, VISTA element or lentiMPRA tile covers it, so the panel"
            " states the interval at that offset from GENCODE's canonical H19 TSS"
        ),
        citations=(
            "Bell and Felsenfeld 2000, Nature 405:482 and Hark et al. 2000, Nature 405:486 (the ICR is a"
            " methylation-sensitive CTCF boundary: unmethylated on the maternal allele it insulates IGF2"
            " from the shared enhancers)",
            "Sparago et al. 2004, Nat Genet 36:958 (microdeletions of the ICR cause Beckwith-Wiedemann)",
            "Gicquel et al. 2005, Nat Genet 37:1003 (loss of ICR methylation causes Silver-Russell)",
        ),
        answer_from=("human imprinting disorders", "mouse deletions", "allele-specific methylation"),
        note=(
            "the imprinted case, added because parent of origin is the clearest information that is not"
            " in the sequence at all and no layer here carries it: the same two alleles behave"
            " oppositely depending on which parent they came from. The element's coordinate is ANCHORED:"
            " the published ICR is 2 to 4 kb upstream of the H19 promoter and the anchor is GENCODE's"
            " canonical H19 TSS (chr11:1,997,841, minus strand). The nearest coding TSS is MRPL23 and"
            " the published target is IGF2, 138 kb away"
        ),
    ),
    Expect(
        locus="PMP22_CMT1A",
        chrom="chr17",
        element=(15_176_315, 15_272_292),
        window=(14_100_000, 15_600_000),
        classes=("program",),
        targets=("PMP22",),
        direction="activates",
        tissues=("Schwann cell", "peripheral nerve"),
        gtex_tissues=("Nerve_Tibial",),
        cells=(),
        distance=0,
        variants=(
            {
                "rsid": None,
                "pos": None,
                "effect": "a 1.4 Mb duplication between the CMT1A-REP repeats causes CMT1A and the"
                " reciprocal deletion causes HNPP: the variable is copy number, not a base",
            },
        ),
        beyond_sequence=("copy number",),
        citations=(
            "Lupski et al. 1991, Cell 66:219 (a 1.4 Mb duplication containing PMP22 causes CMT1A)",
            "Chance et al. 1993, Cell 72:143 (the reciprocal deletion causes hereditary neuropathy with"
            " liability to pressure palsies)",
            "Li et al. 2013, J Neurol Sci 335:29 (dosage, not sequence, sets the severity)",
        ),
        answer_from=("human pedigrees", "dosage mapping", "transgenic rodents"),
        note=(
            "the dosage case, added because the four classes have no term for copy number and the"
            " project has no layer that carries it: three copies of an intact gene cause the disease and"
            " one copy causes a different one. The element is PMP22's GENCODE span, the dosage-sensitive"
            " unit; the window is the CMT1A-REP interval. Everything except the target's identity is"
            " expected to be pending, which is the point of including it"
        ),
    ),
    Expect(
        locus="CYP2D6",
        chrom="chr22",
        element=(42_120_961, 42_131_727),
        window=(41_900_000, 42_400_000),
        classes=("storage", "syntax"),
        targets=("CYP2D6",),
        direction="coding",
        tissues=("liver",),
        gtex_tissues=("Liver",),
        cells=("HepG2",),
        distance=0,
        variants=(
            {
                "rsid": "rs3892097",
                "pos": None,
                "effect": "CYP2D6*4, the commonest European loss-of-function allele",
            },
            {"rsid": "rs1065852", "pos": None, "effect": "CYP2D6*10, common in East Asia, reduced activity"},
        ),
        value_domain="many",
        frequency="*4 is about 0.19 in Europeans and rare in East Asia; *10 is about 0.5 in East Asia",
        citations=(
            "Gaedigk et al. 2018, Clin Pharmacol Ther 103:399 (PharmVar: over a hundred named star"
            " alleles, four metaboliser phenotypes)",
            "Ingelman-Sundberg 2005, Pharmacogenomics J 5:6 (the clinical value domain)",
        ),
        answer_from=("pharmacogenomic allele curation", "enzyme activity assays"),
        note=(
            "a many-valued slot outside the MHC, added to ask whether the value-domain reading"
            " generalises. The prediction written down before the run is that it will NOT read many:"
            " CYP2D6's diversity is star alleles, gene conversion with CYP2D7 and copy number, not a"
            " dense field of common positions, so counting positions per kilobase should underread it"
            " the same way it could not see ABO's frameshift"
        ),
    ),
    Expect(
        locus="MC1R",
        chrom="chr16",
        element=(89_912_118, 89_920_977),
        window=(89_700_000, 90_100_000),
        classes=("storage", "syntax"),
        targets=("MC1R",),
        direction="coding",
        tissues=("melanocyte",),
        gtex_tissues=(),
        cells=(),
        distance=0,
        variants=(
            {
                "rsid": "rs1805007",
                "pos": None,
                "ref": "C",
                "alt": "T",
                "protein_change": "p.Arg151Cys",
                "effect": "R151C, a red-hair R allele",
            },
            {
                "rsid": "rs1805008",
                "pos": None,
                "ref": "C",
                "alt": "T",
                "protein_change": "p.Arg160Trp",
                "effect": "R160W, a red-hair R allele",
            },
            {"rsid": "rs1805009", "pos": None, "effect": "D294H, a red-hair R allele"},
        ),
        value_domain="few",
        frequency="the R alleles are a European specialty: R151C about 0.08 there and near absent elsewhere",
        citations=(
            "Valverde et al. 1995, Nat Genet 11:328 (MC1R variants and red hair)",
            "Box et al. 1997, Hum Mol Genet 6:1891 (R151C, R160W and D294H are the strong R alleles)",
            "Beaumont et al. 2007, Hum Mol Genet 16:2249 (R against r alleles by receptor function)",
        ),
        answer_from=("association with hair and skin colour", "receptor signalling assays"),
        note=(
            "a few-valued slot with a sharp population structure, and a third coding calibration beside"
            " APP and TP53: the translation engine should derive p.Arg151Cys and p.Arg160Trp from the"
            " reference sequence. The element is MC1R's GENCODE span"
        ),
    ),
)
#: the twelve loci whose verdicts were settled before the panel was widened (2026-09-14). The CI gate
#: pins these; the five added on 2026-09-15 are reported but not yet pinned, because a verdict has to
#: be seen before it can be defended.
GATE_PANEL: tuple[str, ...] = (
    "SHH_ZRS",
    "HERC2_OCA2",
    "MCM6_LCT",
    "HBB_LCR",
    "FTO_IRX3",
    "BCL11A_enhancer",
    "MYC_8q24",
    "ABO",
    "HLA_DRB1",
    "APP",
    "TP53",
    "HOXD",
)
#: the published HOXD boundary, between HOXD11 and HOXD13 (Rodriguez-Carballo et al. 2017)
HOXD_BOUNDARY = (176_096_240, 176_109_754)
#: candidates considered and left out, with the reason (docs/LOCI-BENCHMARK.md)
NOT_IN_PANEL = {
    "IRX5": "kept inside FTO_IRX3 as a second published target rather than as a locus of its own",
    "HBB_HS2_alone": (
        "the single HS2 core has no hg38 coordinate we can cite without lifting one over, so the whole"
        " locus control region is the element and the hypersensitive ladder is derived from DNase"
    ),
    "HLA_B": "one HLA gene is enough for the value-domain question; DRB1 is the classical exon-2 case",
    "LCT_African_alleles": (
        "G-13915 and C-14010 are published but their enhancer window overlaps rs4988235's, so the"
        " expectation would not be independent"
    ),
    # considered for the 2026-09-15 widening and dropped, with the reason
    "SMN1_SMN2": (
        "the cleanest dosage phenotype there is, but SMN1 and SMN2 are a segmental duplication"
        " differing at a handful of bases, so no element here has an unambiguous hg38 coordinate and"
        " every layer's answer would be a statement about read mapping. PMP22 carries dosage instead"
    ),
    "RHD_deletion": (
        "Rh negative is a whole-gene deletion, a two-valued slot whose value is copy number: it asks"
        " the same question as PMP22_CMT1A and chr1 is not swept"
    ),
    "FUT2_secretor": (
        "rs601338 is a clean two-valued slot with sharp population structure, but HERC2/OCA2 already"
        " tests a two-valued slot and this adds no failure mode the panel is missing"
    ),
    "AMY1_copy_number": (
        "copy number again, and the published association with starch in the diet is contested, so the"
        " expectation would not be unambiguous"
    ),
    "EPHA4_IHH_TAD_boundaries": (
        "the Lupianez 2015 rearrangements are the best evidence that node boundaries matter, but the"
        " causal unit is a structural variant spanning a boundary and the panel has no way to write"
        " one down as an element"
    ),
}


def panel_by_name() -> dict[str, Expect]:
    return {e.locus: e for e in PANEL}


# ------------------------------------------------------------------------ chromosome context
class Chromosome:
    """One chromosome's local layers, loaded once: elements, nodes, genes, sequence, peaks."""

    def __init__(self, chrom: str, results_dir: Path = RESULTS_DIR) -> None:
        from genomeos.genome import Annotation, IndexedGenome, default_gencode, reference_fasta
        from genomeos.genome.domains import infer_domains
        from genomeos.genome.regulatory import load_ccres

        self.chrom = chrom
        self.results_dir = results_dir
        self.ccres = load_ccres(chrom)
        gff = default_gencode({chrom})
        fasta = reference_fasta(chrom)
        if not self.ccres or gff is None or not fasta.exists():
            raise FileNotFoundError(f"{chrom} is not fetched; run `genomeos data fetch --chrom {chrom}`")
        self.annotation = Annotation.from_gff3(gff, {chrom})
        self.genome = IndexedGenome(str(fasta))
        self.length = self.genome.lengths[chrom]
        self.domains = infer_domains(chrom, self.length, self.ccres, self.annotation)
        self._dom_starts = [d.start for d in self.domains]
        self.coding = [
            (tss_of(g), g.symbol)
            for g in self.annotation.genes.values()
            if g.locus.chrom == chrom and g.type == "protein_coding"
        ]
        self.coding.sort()
        self._coding_tss = [t for t, _ in self.coding]
        # `load_ccres` returns the registry in FILE order, which is several sorted runs concatenated
        # (one per class group), not one sorted list: every chromosome checked has three descending
        # steps in it. `ccres_in` bisects, and a bisect over an unsorted list silently returns a
        # plausible index, so this keeps its own sorted copy. `self.ccres` is left exactly as it came
        # because the domain model is handed it and must see the registry's own order.
        self._ccre_sorted = sorted(self.ccres, key=lambda c: c.start)
        self._ccre_starts = [c.start for c in self._ccre_sorted]
        self._peaks: dict[str, Any] = {}

    def close(self) -> None:
        self.genome.close()

    # -- geometry
    def domain_at(self, pos: int):
        i = bisect.bisect_right(self._dom_starts, pos) - 1
        d = self.domains[i] if 0 <= i < len(self.domains) else None
        return d if d is not None and d.start <= pos < d.end else None

    def nearest_coding(self, pos: int) -> tuple[str, int] | None:
        """(symbol, distance) of the nearest coding TSS anywhere on the chromosome."""
        if not self.coding:
            return None
        i = bisect.bisect_left(self._coding_tss, pos)
        best = None
        for j in (i - 1, i):
            if 0 <= j < len(self.coding):
                t, sym = self.coding[j]
                if best is None or abs(t - pos) < best[1]:
                    best = (sym, abs(t - pos))
        return best

    def nearest_coding_in_node(self, pos: int) -> tuple[str, int] | None:
        d = self.domain_at(pos)
        if d is None:
            return None
        lo = bisect.bisect_left(self._coding_tss, d.start)
        best = None
        for t, sym in self.coding[lo:]:
            if t >= d.end:
                break
            if best is None or abs(t - pos) < best[1]:
                best = (sym, abs(t - pos))
        return best

    def tss(self, symbol: str) -> int | None:
        for t, sym in self.coding:
            if sym == symbol:
                return t
        g = self.annotation.genes.get(symbol) or next(
            (x for x in self.annotation.genes.values() if x.symbol == symbol), None
        )
        return tss_of(g) if g is not None else None

    def gc(self, start: int, end: int) -> float | None:
        from genomeos.coords import Locus

        seq = str(self.genome.fetch(Locus(self.chrom, start, end))).upper()
        acgt = sum(seq.count(b) for b in "ACGT")
        return round((seq.count("G") + seq.count("C")) / acgt, 4) if acgt >= (end - start) * 0.9 else None

    # -- layers
    def ccres_in(self, start: int, end: int) -> list:
        """Every registry element overlapping [start, end), over the sorted copy.

        The 10 kb rewind covers an element that starts before the window and reaches into it. The
        longest cCRE measured across chr1, 2, 7, 11, 16, 19 and X is 350 bp, so the slack is about
        28 times what it needs to be; it is kept wide because being wrong here is silent.
        """
        i = bisect.bisect_left(self._ccre_starts, start - 10_000)
        out = []
        for c in self._ccre_sorted[i:]:
            if c.start >= end:
                break
            if c.end > start:
                out.append(c)
        return out

    def peaks(self, cell: str):
        from genomeos.genome.reader import PeakIndex, load_peaks

        if cell not in self._peaks:
            rows = load_peaks(cell, self.chrom)
            self._peaks[cell] = PeakIndex(rows) if rows else None
        return self._peaks[cell]


def tss_of(g) -> int:
    """The canonical transcript's start, falling back to the gene span.

    GENCODE gene spans include read-through transcripts: HBG2 and HBE1 span to 5.5 Mb on chr11, so a
    gene-span TSS puts two globin genes a quarter of a megabase from where they are. The canonical
    transcript is the only honest start.
    """
    from genomeos.coords import Strand

    minus = g.locus.strand is Strand.MINUS
    ts = list(g.transcripts.values())
    canon = [t for t in ts if "Ensembl_canonical" in t.tags]
    loc = canon[0].locus if canon else g.locus
    return loc.end - 1 if minus else loc.start


# ------------------------------------------------------------------------------- the readers
def read_registry(ch: Chromosome, start: int, end: int) -> dict[str, Any]:
    """The ENCODE registry over the window: an annotation lookup, and nothing more."""
    els = ch.ccres_in(start, end)
    return {
        "layer": "registry",
        "provenance": "looked_up",
        "elements": len(els),
        "classes": sorted({c.cls for c in els}),
        "ids": [c.id for c in els[:8]],
        "evidence": "curated: ENCODE cCRE registry v3 (a chromatin class, never a target)",
    }


def read_node(ch: Chromosome, start: int, end: int) -> dict[str, Any]:
    """The node model's own answer: nearest coding TSS inside the CTCF domain."""
    mid = (start + end) // 2
    d = ch.domain_at(mid)
    in_node = ch.nearest_coding_in_node(mid)
    anywhere = ch.nearest_coding(mid)
    return {
        "layer": "node",
        "provenance": "heuristic",
        "domain": d.id if d else None,
        "domain_span": [d.start, d.end] if d else None,
        "coding_genes_in_node": d.coding_genes if d else None,
        "genes_in_node": sorted(d.genes)[:24] if d else [],
        "target": in_node[0] if in_node else None,
        "distance": in_node[1] if in_node else None,
        "nearest_coding_anywhere": anywhere[0] if anywhere else None,
        "nearest_distance": anywhere[1] if anywhere else None,
        "evidence": EVIDENCE["heuristic"],
    }


def _deletion_rows(chrom: str, start: int, end: int, results_dir: Path) -> list[dict[str, Any]]:
    """Every element inside the window that some AlphaGenome deletion run has already scored."""
    from genomeos.attribution.targets import RUNS, run_elements

    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    # intervals the panel stated itself, because no registry drew one over the published element
    # (docs/LOCI-BENCHMARK.md 11a). Kept first and labelled, so a rate built on one is never pooled
    # with the registry-derived rates without the split showing.
    for name in STATED_INTERVAL_RESULTS:
        for e in (load_result(name, results_dir) or {}).get("elements", []):
            if e.get("chrom") == chrom and e["start"] < end and e["end"] > start:
                seen.add(e["id"])
                rows.append({**e, "run": "stated_interval"})
    # the panel's own deletions first (scripts/loci_score.py: the windows, not whole chromosomes)
    for e in (load_result("loci_deletions", results_dir) or {}).get("elements", []):
        if e.get("chrom") == chrom and e["start"] < end and e["end"] > start:
            seen.add(e["id"])
            rows.append({**e, "run": "loci_windows"})
    for name in RUNS:
        try:
            els = run_elements(name, chrom, results_dir)
        except FileNotFoundError:
            els = []
        for e in els:
            if e.get("start", 0) < end and e.get("end", 0) > start and e["id"] not in seen:
                seen.add(e["id"])
                rows.append({**e, "run": name})
    for stem in (f"vista_{chrom}", f"mpra_{chrom}"):
        r = load_result(stem, results_dir) or {}
        for e in r.get("rows", []):
            inside = e.get("start", 0) < end and e.get("end", 0) > start and e.get("id") not in seen
            if inside and (e.get("predicted") or e.get("predicted_coding")):
                seen.add(e["id"])
                rows.append({**e, "run": stem.split("_")[0]})
    rows.sort(key=lambda e: e["start"])
    return rows


def read_reach(ch: Chromosome, expect: Expect) -> dict[str, Any]:
    """Whether the deletion layer could name this locus's published target *at all*.

    The scorer resizes its input to `MODEL_WINDOW` centred on the element, so a gene whose body
    lies wholly outside that window is not in the model's input and cannot be returned however the
    element is scored. Asking anyway buys a guaranteed negative, and reading one as a miss is the
    benchmark grading a question it never put.

    This is why the ZRS's own deletion named LMBR1, the gene it sits inside, rather than SHH: SHH is
    979 kb from the element and the reach is 524 kb, so SHH was never a candidate. The benchmark
    recorded that as a miss for months. It is not a miss; it is out of range, and the two have to be
    counted apart or the model is blamed for a limit the panel chose.

    Free: no request, no network, GENCODE bodies against arithmetic.
    """
    mid = (expect.element[0] + expect.element[1]) // 2
    half = MODEL_WINDOW // 2
    w0, w1 = mid - half, mid + half
    out: list[dict[str, Any]] = []
    inside: list[str] = []
    for symbol in expect.targets:
        g = ch.annotation.genes.get(symbol) or next(
            (x for x in ch.annotation.genes.values() if x.symbol == symbol), None
        )
        if g is None:
            out.append({"target": symbol, "reason": "no GENCODE record for this symbol"})
            continue
        if g.locus.start < w1 and g.locus.end > w0:
            inside.append(symbol)
            continue
        edge = min(abs(g.locus.start - mid), abs(g.locus.end - mid))
        out.append({"target": symbol, "distance_bp": edge, "reason": "gene body outside the model's input"})
    return {
        "layer": "reach",
        "provenance": "derived",
        "model_window": MODEL_WINDOW,
        "reach_each_way": half,
        "element_midpoint": mid,
        "targets_in_reach": inside,
        "targets_out_of_reach": out,
        "askable": bool(inside),
        "evidence": (
            "derived: GENCODE gene bodies against the scorer's 1 Mb input window, centred on the"
            " element. No request"
        ),
    }


def read_deletion(ch: Chromosome, start: int, end: int, results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """The model layer: which gene moves when an element here is deleted, from runs already made.

    No request is made. When no run covers the window the reading is `pending` with the cost of
    asking, so the locus is never silently dropped.
    """
    rows = _deletion_rows(ch.chrom, start, end, results_dir)
    genes: dict[str, dict[str, Any]] = {}
    for e in rows:
        for key in ("predicted_coding", "predicted"):
            p = e.get(key) or {}
            if not p.get("gene"):
                continue
            g = genes.setdefault(
                p["gene"],
                {"gene": p["gene"], "elements": 0, "best": 0.0, "action": p["action"], "tissue": p["tissue"]},
            )
            g["elements"] += 1
            if abs(p["log2_fold_change"]) > abs(g["best"]):
                g["best"] = p["log2_fold_change"]
                g["action"] = p["action"]
                g["tissue"] = p["tissue"]
            break
    by_cell: dict[str, float] = {}
    for e in rows:
        for cell, v in (e.get("predicted_coding_by_cell") or {}).items():
            by_cell[cell] = round(by_cell.get(cell, 0.0) + v, 4)
    ranked = sorted(genes.values(), key=lambda g: -abs(g["best"]))
    return {
        "layer": "deletion",
        "provenance": "derived",
        "elements_scored": len(rows),
        "runs": sorted({e["run"] for e in rows}),
        "from_stated_interval": any(e["run"] == "stated_interval" for e in rows),
        "targets": ranked[:6],
        "target": ranked[0]["gene"] if ranked else None,
        "action": ranked[0]["action"] if ranked else None,
        "tissue": ranked[0]["tissue"] if ranked else None,
        "summed_by_cell": by_cell or None,
        "pending": None
        if rows
        else (
            "no AlphaGenome deletion has been scored inside this window; the all-element sweep has"
            f" finished {', '.join(swept_chromosomes(results_dir)) or 'no chromosome yet'}."
            " One element is about 0.76 requests"
        ),
        "evidence": "predicted: expression change on deleting an element (AlphaGenome), already computed",
    }


def read_eqtl(
    ch: Chromosome, start: int, end: int, pad: int = 500, knowledge: Path = GTEX_DIR
) -> dict[str, Any]:
    """GTEx v8 eQTLs inside the window: measured evidence for which gene, and in which tissue."""
    hits = _gtex_hits(ch.chrom, knowledge)
    sym = symbols()
    rows = [h for h in hits if start - pad <= h["pos"] <= end + pad]
    genes: dict[str, dict[str, Any]] = {}
    for h in rows:
        g = genes.setdefault(
            sym.get(h["gene_id"], h["gene_id"]), {"tissues": set(), "best_p": 1.0, "slope": 0.0}
        )
        g["tissues"].add(h["tissue"])
        if h["pval"] < g["best_p"]:
            g["best_p"] = h["pval"]
            g["slope"] = h["slope"]
    ranked = sorted(genes.items(), key=lambda kv: kv[1]["best_p"])
    return {
        "layer": "eqtl",
        "provenance": "derived",
        "eqtls": len(rows),
        "egenes": [
            {
                "gene": g,
                "tissues": sorted(v["tissues"])[:8],
                "n_tissues": len(v["tissues"]),
                "best_p": v["best_p"],
                "slope": v["slope"],
            }
            for g, v in ranked[:6]
        ],
        "target": ranked[0][0] if ranked else None,
        "pending": None if hits else f"no GTEx hits distilled for {ch.chrom} under {knowledge}",
        "evidence": "experimental: GTEx v8 significant single-tissue cis-eQTL pairs",
    }


def read_reader(ch: Chromosome, start: int, end: int) -> dict[str, Any]:
    """ENCODE DNase over the window in eleven cell types: which cells hold it open, and how strongly."""
    open_in: dict[str, float] = {}
    missing = []
    peaks: dict[str, list[tuple[int, int, float]]] = {}
    for cell in READER_CELLS:
        idx = ch.peaks(cell)
        if idx is None:
            missing.append(cell)
            continue
        got = idx.overlapping(start, end)
        if got:
            open_in[cell] = round(max(v for _s, _e, v in got), 2)
            peaks[cell] = got
    best = max(open_in, key=lambda c: open_in[c]) if open_in else None
    return {
        "layer": "reader",
        "provenance": "derived",
        "cells_read": len(READER_CELLS) - len(missing),
        "open_in": open_in,
        "cells_open": sorted(open_in),
        "strongest_cell": best,
        "peaks_K562": peaks.get("K562", [])[:12],
        "pending": f"no DNase peaks on {ch.chrom} for {', '.join(missing)}" if missing else None,
        "evidence": "experimental: ENCODE DNase-seq narrowPeak, eleven cell types (the reader)",
    }


def read_mpra(ch: Chromosome, start: int, end: int) -> dict[str, Any]:
    """ENCODE4 lentiMPRA: was a 200 bp piece of this window active in K562, HepG2 or WTC11."""
    p = Path("data/knowledge/mpra") / f"rows_{ch.chrom}.json"
    if not p.exists():
        return {"layer": "mpra", "provenance": "derived", "pending": f"lentiMPRA rows for {ch.chrom} absent"}
    rows = [r for r in json.loads(p.read_text()) if r["start"] < end and r["end"] > start]
    active = {c: sum(1 for r in rows if (r.get("active") or {}).get(c)) for c in ("K562", "HepG2", "WTC11")}
    return {
        "layer": "mpra",
        "provenance": "derived",
        "tested": len(rows),
        "active": active if rows else None,
        "best": max(
            ((r["activity"], r["start"]) for r in rows if r.get("activity")),
            key=lambda x: max(x[0].values()),
            default=(None, None),
        )[0],
        "pending": "no lentiMPRA element was tested inside this window" if not rows else None,
        "evidence": "experimental: ENCODE4 joint lentiMPRA, log2 RNA/DNA in K562, HepG2 and WTC11",
    }


def read_satmut(ch: Chromosome, start: int, end: int, expect: Expect | None = None) -> dict[str, Any]:
    """Saturation mutagenesis over the window: which bases matter, measured.

    Targeted, not derived: Kircher et al. chose these twenty-one elements because they were already
    known, so this layer can confirm a causal base and never discover one.
    """
    from genomeos.knowledge.satmut import DATA_PATH, base_table, by_element, load

    if not DATA_PATH.exists():
        return {"layer": "satmut", "provenance": "targeted", "pending": "satmut table not cached"}
    rows = [r for r in load() if r["chrom"] == ch.chrom and start <= r["pos"] <= end]
    if not rows:
        return {
            "layer": "satmut",
            "provenance": "targeted",
            "pending": "no saturation-mutagenesis element covers this window",
        }
    out: dict[str, Any] = {"layer": "satmut", "provenance": "targeted", "elements": {}}
    for name, rs in by_element(rows).items():
        table = base_table(rs)
        functional = [p for p, b in table.items() if b["functional"]]
        strong = [p for p, b in table.items() if b["strong"]]
        entry = {
            "bases_measured": len(table),
            "functional": len(functional),
            "strong": len(strong),
            "top": sorted(table.items(), key=lambda kv: -abs(kv[1]["effect"]))[:5],
        }
        for v in expect.variants if expect else ():
            if v.get("pos") and v["pos"] in table:
                b = table[v["pos"]]
                entry.setdefault("published_variant", {})[v.get("rsid") or str(v["pos"])] = {
                    "functional": b["functional"],
                    "strong": b["strong"],
                    "effect": b["effect"],
                    "rank_by_effect": 1
                    + sorted((abs(x["effect"]) for x in table.values()), reverse=True).index(
                        abs(b["effect"])
                    ),
                }
        out["elements"][name] = entry
    out["evidence"] = "experimental: Kircher et al. 2019 saturation mutagenesis MPRA (GRCh38)"
    return out


def read_lookups(ch: Chromosome, start: int, end: int, results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """What a pure annotation lookup gives: VISTA's tissues, the GWAS Catalog's gene, ClinVar's gene."""
    vista = []
    r = load_result(f"vista_{ch.chrom}", results_dir) or {}
    for x in r.get("rows", []):
        if x["start"] < end and x["end"] > start:
            vista.append({"id": x["id"], "status": x["status"], "tissues": x.get("tissues") or []})
    gwas = []
    p = Path("data/knowledge/gwas/hits.tsv")
    if p.exists():
        with p.open() as fh:
            next(fh, None)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0] == ch.chrom and start <= int(f[1]) <= end:
                    gwas.append({"rs": f[2], "trait": f[3][:60], "mapped_gene": f[4]})
    clinvar: list[dict[str, Any]] = []
    p = Path("data/knowledge/clinvar/pathogenic.tsv.gz")
    if p.exists():
        with gzip.open(p, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if f[0] == ch.chrom and start <= int(f[1]) <= end:
                    clinvar.append({"gene": f[4], "significance": f[5]})
                    if len(clinvar) >= 50:
                        break
    return {
        "layer": "lookups",
        "provenance": "looked_up",
        "vista": vista,
        "vista_tissues": sorted({t for v in vista for t in v["tissues"]}),
        "gwas": gwas[:10],
        "gwas_genes": sorted({g["mapped_gene"] for g in gwas}),
        "clinvar_genes": sorted({c["gene"] for c in clinvar}),
        "evidence": EVIDENCE["looked_up"],
    }


def read_values(ch: Chromosome, expect: Expect, frequencies: dict[str, Any] | None = None) -> dict[str, Any]:
    """The value slot: who carries which allele among the people we hold, and how common it is.

    The imported genomes are three (the GIAB trio), which is enough to say a position is variable and
    nothing about frequency; the frequency comes from 1000 Genomes through Ensembl, cached by the
    script (`read_frequencies`, gnomAD v4.1.1 genomes).
    """
    from genomeos.genome.individuals import rows_in, sources

    start, end = expect.element
    carriers: dict[int, dict[str, str]] = {}
    people = []
    for name, path, _ev in sources(ch.chrom):
        if name == "demo":
            continue
        people.append(name)
        for f in rows_in(Path(path), start + 1, end):
            gt = f[9].split(":")[0] if len(f) > 9 else "./."
            alleles = [a for a in gt.replace("|", "/").split("/") if a not in (".", "")]
            if not any(a != "0" for a in alleles):
                continue
            zyg = "hom" if alleles and all(a == alleles[0] for a in alleles) else "het"
            carriers.setdefault(int(f[1]), {})[name] = f"{f[3]}>{f[4].split(',')[0]} {zyg}"
    named = {}
    for v in expect.variants:
        if not v.get("pos"):
            continue
        key = v.get("rsid") or str(v["pos"])
        named[key] = {
            "pos": v["pos"],
            "carried_by": carriers.get(v["pos"], {}),
            "variable_in_the_trio": v["pos"] in carriers,
            "frequency": ((frequencies or {}).get("named") or {}).get(key),
        }
    return {
        "layer": "values",
        "provenance": "derived",
        "people": people,
        "variable_positions_in_element": len(carriers),
        "named_variants": named,
        "evidence": (
            "measured: each person's own variant calls (GIAB HG002/HG003/HG004);"
            " curated: 1000 Genomes allele frequencies through Ensembl"
        ),
    }


def read_constraint(chrom: str, windows: list[tuple[int, int]], progress=None) -> list[dict[str, Any]]:
    """Both constraint axes over a chromosome's windows in one pass each (Zoonomia, Gnocchi)."""
    from genomeos.attribution.constraint import PHYLOP_THRESHOLD, phylop_over_blocks
    from genomeos.attribution.variation import case_of, gnocchi_over, stats_dict

    ivs = sorted(windows)
    phy, _c1 = phylop_over_blocks(chrom, ivs, PHYLOP_THRESHOLD, progress=progress)
    gno, _c2 = gnocchi_over(chrom, ivs, progress=progress)
    out = {}
    for iv, p, g in zip(ivs, phy, gno, strict=True):
        gd = stats_dict(g)
        out[iv] = {
            "layer": "constraint",
            "provenance": "derived",
            "mammals": {
                "bases": p.bases,
                "mean": round(p.mean, 3) if p.mean is not None else None,
                "maximum": round(p.maximum, 3) if p.bases else None,
                "fraction_above": round(p.fraction_above, 4) if p.fraction_above is not None else None,
            },
            "humans": gd,
            "case": case_of(p.fraction_above, gd["fraction_above"] if gd else None),
            "evidence": f"{EVIDENCE['derived']}; thresholds phyloP {PHYLOP_THRESHOLD}, Gnocchi 2.18",
        }
    return [out[iv] for iv in windows]


def variable_positions(ch: Chromosome, start: int, end: int) -> dict[int, set[str]]:
    """Positions in the window where one of the people we hold carries a non-reference SNV."""
    from genomeos.genome.individuals import rows_in, sources

    positions: dict[int, set[str]] = {}
    for name, path, _ev in sources(ch.chrom):
        if name == "demo":
            continue
        for f in rows_in(Path(path), start + 1, end):
            ref, alt = f[3], f[4].split(",")[0]
            gt = f[9].split(":")[0] if len(f) > 9 else "./."
            called = any(a not in ("0", ".", "") for a in gt.replace("|", "/").split("/"))
            if len(ref) == 1 and len(alt) == 1 and called:
                positions.setdefault(int(f[1]), set()).add(name)
    return positions


def read_syntax_values_many(
    ch: Chromosome, windows: list[tuple[int, int]], expects: list[Expect | None] | None = None
) -> list[dict[str, Any]]:
    """Syntax against values over many windows of one chromosome, in one bigWig session.

    Every position where one of the people we hold differs from the reference gets its phyloP; a
    value on constrained sequence is a "value in syntax" (attribution/syntax.py's reading). The
    benchmark asks whether the published causal base comes out near the top of that ranking without
    being named, which is how `genomeos syntax --gene HERC2` found rs12913832. One session for the
    whole chromosome, because opening one per window cost six hours on the first run.
    """
    from genomeos.attribution.bigwig import BigWig
    from genomeos.attribution.constraint import PHYLOP_241_URL, PHYLOP_THRESHOLD

    per_window = [variable_positions(ch, s, e) for s, e in windows]
    every = sorted({p for d in per_window for p in d})
    phylop: dict[int, float] = {}
    if every:
        bw = BigWig(PHYLOP_241_URL)
        try:
            stats = bw.summarise(ch.chrom, [(p - 1, p) for p in every], PHYLOP_THRESHOLD)
        finally:
            bw.close()
        phylop = {p: round(s.maximum, 2) for p, s in zip(every, stats, strict=True) if s.bases}
    out = []
    for i, positions in enumerate(per_window):
        expect = (expects or [None] * len(windows))[i]
        rows = sorted(((p, phylop[p]) for p in positions if p in phylop), key=lambda x: -x[1])
        ranks = {p: n + 1 for n, (p, _v) in enumerate(rows)}
        published = {}
        for v in expect.variants if expect else ():
            if v.get("pos"):
                published[v.get("rsid") or str(v["pos"])] = {
                    "variable_in_the_trio": v["pos"] in positions,
                    "phylop": phylop.get(v["pos"]),
                    "rank": ranks.get(v["pos"]),
                    "of": len(rows),
                }
        out.append(
            {
                "layer": "syntax_values",
                "provenance": "derived",
                "values": len(rows),
                "values_in_syntax": sum(1 for _p, x in rows if x >= PHYLOP_THRESHOLD),
                "top": [{"pos": p, "phylop": x, "carriers": sorted(positions[p])} for p, x in rows[:5]],
                "published": published,
                "evidence": "curated: Zoonomia phyloP per base; measured: the GIAB trio's own SNVs",
            }
        )
    return out


def read_syntax_values(ch: Chromosome, start: int, end: int, expect: Expect | None = None) -> dict[str, Any]:
    return read_syntax_values_many(ch, [(start, end)], [expect])[0]


def _number(text: Any) -> float:
    """A bigBed field as a number; gnomAD writes N/A where a group has no call."""
    try:
        return float(text)
    except (TypeError, ValueError):
        return 0.0


def read_frequencies(ch: Chromosome, start: int, end: int, expect: Expect | None = None) -> dict[str, Any]:
    """How common the values are, from gnomAD v4.1.1 genomes read by range as a bigBed.

    Two readings: the density of common variable positions in the window, which is the derived
    size of the value domain, and each named variant's own frequency with the population carrying
    it most often. Read through the bigBed reader of `attribution/human_panel.py` (genomeos-h1's,
    imported and not copied); Ensembl's overlap endpoint carries no frequencies and its variation
    endpoint is used only to check the hard-coded positions.
    """
    from genomeos.attribution.human_panel import open_bigbed

    try:
        bb = open_bigbed("gnomad_snv")
        rows = bb.query(ch.chrom, [(start, end)])
    except (OSError, ValueError, KeyError) as ex:
        return {"layer": "frequencies", "provenance": "derived", "pending": str(ex)[:160]}
    common = {r["chromStart"] for r in rows if _number(r.get("AF")) >= COMMON_AF}
    kb = (end - start) / 1000
    named: dict[str, Any] = {}
    for v in expect.variants if expect else ():
        if not v.get("pos"):
            continue
        hit = next((r for r in rows if r["chromStart"] == v["pos"] - 1), None)
        named[v.get("rsid") or str(v["pos"])] = (
            {
                "af": _number(hit.get("AF")),
                "rsid_in_gnomad": hit.get("rsId"),
                "commonest_in": hit.get("grpmax"),
                "af_there": _number(hit.get("AF_grpmax")),
                "allele": f"{hit.get('ref')}>{hit.get('alt')}",
            }
            if hit
            else {"pending": "not in the gnomAD genomes track at this position"}
        )
    return {
        "layer": "frequencies",
        "provenance": "derived",
        "variants_in_window": len(rows),
        "common": len(common),
        "per_kb": round(len(common) / kb, 2) if kb else None,
        "named": named,
        "threshold": COMMON_AF,
        "evidence": "curated: gnomAD v4.1.1 genomes, frequencies from up to 76,215 genomes (UCSC bigBed)",
    }


def read_coding(ch: Chromosome, expect: Expect) -> dict[str, Any]:
    """The coding calibration: the gene holding each named variant, and the protein change the engine
    derives by translating the canonical transcript over the local sequence."""
    from genomeos.flow import trace_gene

    out: dict[str, Any] = {"layer": "coding", "provenance": "derived", "variants": {}}
    module = None
    for v in expect.variants:
        if not v.get("pos"):
            continue
        key = v.get("rsid") or str(v["pos"])
        holders = [
            g
            for g in ch.annotation.genes.values()
            if g.type == "protein_coding" and g.locus.start < v["pos"] <= g.locus.end
        ]
        entry: dict[str, Any] = {
            "genes": sorted(g.symbol for g in holders),
            "expected": v.get("protein_change"),
        }
        if v.get("ref") and v.get("alt") and len(v["ref"]) == 1 and len(v["alt"]) == 1:
            for g in holders:
                if g.symbol not in expect.targets:
                    continue
                module = module or ch.annotation.to_module("loci")
                tr = trace_gene(ch.genome, g, module.entities[g.id].transcripts)
                sub = tr.substitute(v["pos"] - 1, v["ref"], v["alt"]) if tr is not None else {}
                entry.update(
                    {
                        "transcript": tr.transcript if tr is not None else None,
                        "consequence": sub.get("consequence"),
                        "hgvs_p": sub.get("hgvs_p"),
                        "matches": bool(v.get("protein_change"))
                        and sub.get("hgvs_p") == v.get("protein_change"),
                    }
                )
        else:
            entry["pending"] = (
                "the local trace substitutes single bases; an indel's consequence is not derived"
            )
        out["variants"][key] = entry
    genes = sorted({g for e in out["variants"].values() for g in e["genes"]})
    out["target"] = next((g for g in genes if g in expect.targets), genes[0] if genes else None)
    out["evidence"] = "derived: GenomeOS trace of the canonical transcript on the local sequence"
    return out


def read_gene_input(ch: Chromosome, expect: Expect, results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """The whole input: every already-scored element in the locus window, grouped by the gene it moves."""
    rows = _deletion_rows(ch.chrom, expect.window[0], expect.window[1], results_dir)
    by_gene: dict[str, dict[str, Any]] = {}
    for e in rows:
        p = e.get("predicted_coding") or {}
        if not p.get("gene"):
            continue
        g = by_gene.setdefault(p["gene"], {"gene": p["gene"], "elements": 0, "activating": 0, "summed": 0.0})
        g["elements"] += 1
        g["activating"] += int(p["action"] == "activates")
        g["summed"] = round(g["summed"] + abs(p["log2_fold_change"]), 3)
    ranked = sorted(by_gene.values(), key=lambda g: -g["summed"])
    mine = [g for g in ranked if g["gene"] in expect.targets]
    return {
        "layer": "gene_input",
        "provenance": "derived",
        "window": list(expect.window),
        "elements_scored": len(rows),
        "naming_a_coding_gene": sum(g["elements"] for g in ranked),
        "genes": ranked[:8],
        "published_targets": mine,
        "target": mine[0]["gene"] if mine else (ranked[0]["gene"] if ranked else None),
        "rank_of_first_published_target": next(
            (i + 1 for i, g in enumerate(ranked) if g["gene"] in expect.targets), None
        ),
        "pending": None if rows else "no already-scored element anywhere in the locus window",
        "evidence": "predicted: AlphaGenome deletions already computed, summed |log2| per coding gene",
    }


# ------------------------------------------------------------------- the two ordered readings
def read_hypersensitive_ladder(ch: Chromosome, start: int, end: int) -> dict[str, Any]:
    """The beta-globin question: does the reader find an erythroid-specific ladder of open sites."""
    k = ch.peaks("K562")
    if k is None:
        return {"layer": "ladder", "provenance": "derived", "pending": "no K562 peaks on this chromosome"}
    sites = []
    for s, e, v in k.overlapping(start, end):
        if v < K562_STRONG_PEAK:
            continue
        others = [c for c in READER_CELLS if c != "K562" and (ch.peaks(c) or _NoPeaks()).overlapping(s, e)]
        sites.append({"start": s, "end": e, "signal": v, "also_open_in": others})
    specific = [s for s in sites if len(s["also_open_in"]) <= 1]
    gaps = [b["start"] - a["start"] for a, b in zip(sites, sites[1:], strict=False)]
    return {
        "layer": "ladder",
        "provenance": "derived",
        "strong_K562_sites": len(sites),
        "erythroid_specific": len(specific),
        "sites": sites,
        "median_spacing": sorted(gaps)[len(gaps) // 2] if gaps else None,
        "evidence": "experimental: ENCODE DNase-seq narrowPeak, K562 against the other ten cell types",
    }


class _NoPeaks:
    """An empty peak index, so a missing cell type reads as closed rather than raising."""

    @staticmethod
    def overlapping(_s: int, _e: int) -> list:
        return []


def read_cluster_order(ch: Chromosome, expect: Expect) -> dict[str, Any]:
    """The HOXD question: is the ordered structure of the cluster visible to the machinery.

    Three parts, each separated from what GENCODE simply tells us. That the genes are in the cluster
    and in which order is annotation. What is asked here is whether the *node* model puts a boundary
    where the published one is, whether the cluster's non-coding content is denser and more conserved
    than matched windows, and whether any measurement we hold carries the order of execution.
    """
    genes = [
        (ch.tss(g), g) for g in expect.targets if ch.tss(g) is not None
    ]  # looked up: GENCODE names and places them
    genes.sort()
    start, end = expect.element
    boundaries = [
        (c.start + c.end) // 2 for c in ch.ccres_in(start, end) if c.cls == "CTCF-only"
    ]  # derived from the registry's CTCF-only class, the node model's own boundary proxy
    nodes = sorted(
        {(d.id, d.start, d.end) for d in (ch.domain_at(t) for t, _ in genes) if d}, key=lambda x: x[1]
    )
    lo, hi = HOXD_BOUNDARY
    inside = [b for b in boundaries if lo <= b <= hi]
    opens = {}
    for t, g in genes:
        cells = [c for c in READER_CELLS if (ch.peaks(c) or _NoPeaks()).overlapping(t - 1000, t + 1000)]
        opens[g] = len(cells)
    order = [g for _t, g in genes]
    return {
        "layer": "order",
        "provenance": "derived",
        "genes_in_order": order,
        "published_order_3_to_5": list(expect.targets)[::-1],
        "nodes_over_the_cluster": [{"id": i, "start": s, "end": e} for i, s, e in nodes],
        "ctcf_only_sites": len(boundaries),
        "boundary_in_published_interval": len(inside),
        "published_boundary": list(HOXD_BOUNDARY),
        "node_split_between": _split_between(nodes, genes),
        "promoters_open_in_n_cells": opens,
        "pending": (
            "no measurement in the project carries developmental time or position along the body axis:"
            " the ordered activation of the cluster cannot be tested with the layers we hold, and the"
            " language has no construct for order"
        ),
        "evidence": "inferred: CTCF-only elements as node boundaries; experimental: DNase over the promoters",
    }


def _split_between(nodes: list[tuple[str, int, int]], genes: list[tuple[int, str]]) -> str | None:
    """Which two genes the node boundary falls between, when the cluster's genes sit in two nodes."""
    if len(nodes) < 2:
        return None
    edges = [n[1] for n in nodes[1:]]
    out = []
    for e in edges:
        before = [g for t, g in genes if t < e]
        after = [g for t, g in genes if t >= e]
        if before and after:
            out.append(f"{before[-1]}|{after[0]}")
    return ", ".join(out) or None


# ------------------------------------------------------------------------- GTEx and frequencies
_GTEX: dict[str, list[dict[str, Any]]] = {}
_SYMBOLS: dict[str, str] = {}


def symbols(path: Path = Path("data/cache/gencode_genes.tsv")) -> dict[str, str]:
    """Ensembl gene id (unversioned) to symbol, from the cached GENCODE gene table."""
    global _SYMBOLS
    if not _SYMBOLS and path.exists():
        with path.open() as fh:
            next(fh, None)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                _SYMBOLS[f[0].split(".")[0]] = f[1]
    return _SYMBOLS


def _gtex_hits(chrom: str, knowledge: Path) -> list[dict[str, Any]]:
    key = f"{knowledge}:{chrom}"
    if key not in _GTEX:
        rows: list[dict[str, Any]] = []
        for p in sorted(knowledge.glob("hits_*.tsv")) if knowledge.is_dir() else []:
            with p.open() as fh:
                next(fh, None)
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if f[1] != chrom:
                        continue
                    rows.append(
                        {
                            "tissue": f[0],
                            "pos": int(f[2]),
                            "gene_id": f[5].split(".")[0],
                            "slope": float(f[6]),
                            "pval": float(f[7]),
                        }
                    )
        _GTEX[key] = rows
    return _GTEX[key]


def gtex_intervals(windows: list[tuple[str, int, int]], chunk: int = 8_000):
    """The eQTL reader searches back one element length, so long windows go in as chunks."""
    from genomeos.attribution.eqtl import Intervals

    ivs = Intervals()
    for chrom, s, e in windows:
        for a in range(s, e, chunk):
            b = min(a + chunk, e)
            ivs.add(chrom, a, b, f"{chrom}:{a}-{b}")
    return ivs.freeze()


def ensembl_variant(rsid: str, timeout: int = 15) -> dict[str, Any] | None:
    """One variant through Ensembl: its GRCh38 position and 1000 Genomes population frequencies."""
    import urllib.request

    url = f"https://rest.ensembl.org/variation/human/{rsid}?pops=1&content-type=application/json"
    req = urllib.request.Request(url, headers={"User-Agent": "GenomeOS/0.9 (loci benchmark)"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:  # noqa: S310
            d = json.load(r)
    except OSError:
        return None
    pos = None
    for m in d.get("mappings", []):
        if m.get("assembly_name") == "GRCh38" and m["location"].count(":") == 1:
            loc = m["location"].split(":")[0]
            if loc.isdigit() or loc in ("X", "Y"):
                pos = int(m["location"].split(":")[1].split("-")[0])
    pops = {}
    for p in d.get("populations", []):
        name = p.get("population", "")
        if name.startswith("1000GENOMES:phase_3:") and name.count(":") == 2:
            code = name.rsplit(":", 1)[1]
            if len(code) == 3 and code.isupper():
                pops.setdefault(code, {})[p["allele"]] = round(p["frequency"], 4)
    return {"position": pos, "populations": pops, "maf": d.get("MAF"), "minor_allele": d.get("minor_allele")}


# --------------------------------------------------------------------------- negative controls
@dataclass
class Window:
    chrom: str
    start: int
    end: int
    gc: float | None = None
    distance: int | None = None
    constrained: float | None = None
    of_locus: str = ""
    reason: str = ""
    readings: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        return {
            "chrom": self.chrom,
            "start": self.start,
            "end": self.end,
            "gc": self.gc,
            "distance_to_coding_tss": self.distance,
            "constrained_fraction": self.constrained,
            "of_locus": self.of_locus,
            **({"readings": self.readings} if self.readings else {}),
        }


def candidate_windows(
    ch: Chromosome,
    expect: Expect,
    positive_gc: float,
    positive_distance: int,
    n: int = CANDIDATES_PER_LOCUS,
    seed: int = 11,
) -> list[Window]:
    """Windows of the positive's length, matched on GC and distance to the nearest coding TSS.

    The keep-out list is every panel locus (and its 200 kb of flank) plus VISTA elements and GWAS
    Catalog hits: a "no known function" control must not be a known element in disguise. Constraint
    is matched afterwards, in one range read per chromosome.
    """
    rng = random.Random(f"{seed}:{expect.locus}")
    length = expect.element[1] - expect.element[0]
    keep_out = [
        (e.window[0] - NEGATIVE_KEEP_OUT, e.window[1] + NEGATIVE_KEEP_OUT)
        for e in PANEL
        if e.chrom == ch.chrom
    ]
    vista = load_result(f"vista_{ch.chrom}") or {}
    keep_out += [(x["start"] - 5_000, x["end"] + 5_000) for x in vista.get("rows", [])]
    gwas = Path("data/knowledge/gwas/hits.tsv")
    if gwas.exists():
        with gwas.open() as fh:
            next(fh, None)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0] == ch.chrom:
                    keep_out.append((int(f[1]) - 5_000, int(f[1]) + 5_000))
    keep_out.sort()
    starts = [k[0] for k in keep_out]

    def blocked(s: int, e: int) -> bool:
        i = bisect.bisect_right(starts, e)
        return any(a < e and s < b for a, b in keep_out[max(0, i - 400) : i])

    out: list[Window] = []
    for _ in range(n * 120):
        if len(out) >= n:
            break
        s = rng.randrange(1_000_000, max(1_000_001, ch.length - length - 1_000_000))
        e = s + length
        if blocked(s, e):
            continue
        gc = ch.gc(s, e)
        if gc is None or abs(gc - positive_gc) > MATCH_GC:
            continue
        near = ch.nearest_coding((s + e) // 2)
        if near is None:
            continue
        d = near[1]
        if positive_distance > 0 and abs(d - positive_distance) > MATCH_DISTANCE * positive_distance:
            continue
        if positive_distance == 0 and d > 5_000:
            continue
        out.append(Window(ch.chrom, s, e, gc=gc, distance=d, of_locus=expect.locus))
    return out


def pick_negatives(
    candidates: list[Window], positive_constrained: float | None, keep: int = NEGATIVES_PER_LOCUS
) -> list[Window]:
    """The `keep` candidates whose constrained fraction is closest to the positive's."""
    scored = [c for c in candidates if c.constrained is not None]
    if positive_constrained is None or not scored:
        return candidates[:keep]
    scored.sort(key=lambda c: abs((c.constrained or 0) - positive_constrained))
    return scored[:keep]


# --------------------------------------------------------------------------------- the scoring
def _named(reading: dict[str, Any], key: str = "target") -> str | None:
    return reading.get(key) if reading else None


def score_target(expect: Expect, readings: dict[str, Any]) -> dict[str, Any]:
    """Which layers name a published target gene, and by what provenance.

    A layer's `hit` is strict: the gene it ranks first is a published target. `among` is the lenient
    reading, a published target anywhere in what the layer names; it is reported and never counted in
    the headline. The coding layer hits when the engine's own translation derives the published
    protein change inside the published gene. A looked-up layer hits if it carries the answer at all,
    which is the point of keeping it on a separate line.
    """
    want = set(expect.targets)
    by_layer: dict[str, Any] = {}
    for layer in ("deletion", "eqtl", "gene_input", "coding", "node", "lookups"):
        r = readings.get(layer)
        if r is None:
            continue
        if layer == "lookups":
            named = list(r.get("gwas_genes") or []) + list(r.get("clinvar_genes") or [])
        elif layer == "deletion":
            named = [g["gene"] for g in r.get("targets") or []]
        elif layer == "eqtl":
            named = [g["gene"] for g in r.get("egenes") or []]
        elif layer == "gene_input":
            named = [g["gene"] for g in r.get("genes") or []]
        else:
            named = [r.get("target")] if r.get("target") else []
        named = [n for n in named if n]
        hit = bool(named) and named[0] in want
        if layer == "coding":
            hit = any(v.get("matches") for v in (r.get("variants") or {}).values())
        elif layer == "lookups":
            hit = bool(want & set(named))
        by_layer[layer] = {
            "provenance": r.get("provenance") or "?",
            "named": named[:6],
            "hit": hit,
            "among": bool(want & set(named)),
            "first_hit": next((n for n in named if n in want), None),
            "pending": r.get("pending"),
            # a layer that could not have named the target is reported apart from one that failed to.
            # The headline rates are left alone deliberately: changing a shared benchmark's
            # denominator is the owner's call, and this makes the split visible so it can be made.
            **({"unaskable": r["unaskable"]} if r.get("unaskable") else {}),
        }
    derived = {k: v for k, v in by_layer.items() if v["provenance"] == "derived"}
    return {
        "field": "target",
        "expected": sorted(want),
        "by_layer": by_layer,
        "derived_unaskable": [k for k, v in derived.items() if v.get("unaskable")],
        "hit_derived": any(v["hit"] for v in derived.values()),
        "hit_derived_by": [k for k, v in derived.items() if v["hit"]],
        "hit_derived_among": any(v["among"] for v in derived.values()),
        "hit_heuristic": bool((by_layer.get("node") or {}).get("hit")),
        "hit_looked_up": bool((by_layer.get("lookups") or {}).get("hit")),
    }


def score_cell(expect: Expect, readings: dict[str, Any]) -> dict[str, Any]:
    """The tissue or cell type: the reader's open cells (at the element and at the promoter), the
    tissues GTEx ties to the target, and the track the model's deletion moved most."""
    reader = readings.get("reader") or {}
    promoter = readings.get("promoter") or {}
    eqtl = readings.get("eqtl") or {}
    deletion = readings.get("deletion") or {}
    open_cells = set(reader.get("cells_open") or [])
    promoter_cells = set(promoter.get("cells_open") or [])
    want_cells = set(expect.cells)
    tissues = {t for g in eqtl.get("egenes") or [] if g["gene"] in expect.targets for t in g["tissues"]}
    want_tissues = set(expect.gtex_tissues)
    model_tissue = (deletion.get("tissue") or "").lower()
    model_hit = any(w.split("_")[0].lower() in model_tissue for w in expect.tissues) if model_tissue else None
    element_hit = bool(want_cells & open_cells)
    promoter_hit = bool(want_cells & promoter_cells)
    eqtl_hit = bool(want_tissues & tissues)
    return {
        "field": "cell",
        "expected_cells": sorted(want_cells),
        "expected_gtex": sorted(want_tissues),
        "reader_open_in_expected": sorted(want_cells & open_cells) if want_cells else None,
        "reader_hit": element_hit if want_cells else None,
        "reader_cells_open": len(open_cells),
        "promoter_open_in_expected": sorted(want_cells & promoter_cells) if want_cells else None,
        "promoter_hit": promoter_hit if want_cells else None,
        "promoter_cells_open": len(promoter_cells) if promoter else None,
        "eqtl_tissues_for_target": sorted(tissues)[:8],
        "eqtl_hit": eqtl_hit if want_tissues else None,
        "model_tissue": deletion.get("tissue"),
        "model_hit": model_hit,
        "hit_derived": bool(element_hit or promoter_hit or eqtl_hit),
        "specific": (
            round(1 - (len(open_cells) - 1) / (len(READER_CELLS) - 1), 2)
            if element_hit and open_cells
            else None
        ),
        "unreachable": None
        if (want_cells or want_tissues)
        else "the published tissue has no counterpart in the eleven reader cells or the GTEx panel",
    }


def score_direction(expect: Expect, readings: dict[str, Any]) -> dict[str, Any]:
    """The sign: does the model's action, or the eQTL slope, agree with the published direction."""
    deletion = readings.get("deletion") or {}
    action = None
    for g in deletion.get("targets") or []:
        if g["gene"] in expect.targets:
            action = g["action"]
            break
    return {
        "field": "direction",
        "expected": expect.direction,
        "model_action": action,
        "hit_derived": bool(
            action and expect.direction in ("activates", "represses") and action == expect.direction
        ),
        "judged": bool(action) and expect.direction in ("activates", "represses"),
    }


def score_distance(expect: Expect, ch: Chromosome, readings: dict[str, Any]) -> dict[str, Any]:
    """The distance from element to the target the machinery named, against the published one."""
    mid = sum(expect.element) // 2
    got = None
    for layer in ("deletion", "eqtl"):
        r = readings.get(layer) or {}
        genes = [g["gene"] for g in (r.get("targets") or r.get("egenes") or [])]
        for g in genes:
            if g in expect.targets:
                t = ch.tss(g)
                if t is not None:
                    got = {"layer": layer, "gene": g, "distance": abs(t - mid)}
                    break
        if got:
            break
    ok = None
    if got and expect.distance:
        ok = 0.5 <= got["distance"] / expect.distance <= 2.0
    return {
        "field": "distance",
        "expected": expect.distance,
        "measured": got,
        "within_twofold": ok,
        "node_distance": (readings.get("node") or {}).get("distance"),
    }


def score_variant(expect: Expect, readings: dict[str, Any]) -> dict[str, Any]:
    """The causal base: is it variable in the people we hold, how common is it where it should be,
    does constraint rank it near the top of the window's variable positions, and does a measurement
    (saturation mutagenesis) call it functional."""
    values = readings.get("values") or {}
    sv = readings.get("syntax_values") or {}
    satmut = readings.get("satmut") or {}
    coding = readings.get("coding") or {}
    named = values.get("named_variants") or {}
    pub = {}
    for el in (satmut.get("elements") or {}).values():
        pub.update(el.get("published_variant") or {})
    ranks = sv.get("published") or {}
    top_decile = {
        k: (v["rank"] is not None and v["of"] and v["rank"] <= max(1, round(0.1 * v["of"])))
        for k, v in ranks.items()
    }
    return {
        "field": "variant",
        "expected": [v.get("rsid") for v in expect.variants if v.get("rsid")],
        "variable_in_the_trio": {k: v["variable_in_the_trio"] for k, v in named.items()},
        "frequencies": {k: v.get("frequency") for k, v in named.items()},
        "constraint_rank": ranks or None,
        "in_the_top_decile_by_constraint": top_decile or None,
        "satmut_functional": pub or None,
        "protein_change": {k: v.get("hgvs_p") for k, v in (coding.get("variants") or {}).items()} or None,
        "hit_derived": bool(
            any(top_decile.values()) or any(v.get("matches") for v in (coding.get("variants") or {}).values())
        ),
        "hit_targeted": any(x.get("functional") for x in pub.values()) if pub else None,
        "reading": (
            "derived means the machinery pointed at the base by itself: constraint ranking it in the"
            " top tenth of the window's variable positions, or the translation engine deriving the"
            " published protein change. Carrying the frequency of a base someone else named is not a hit."
        ),
    }


def score_class(expect: Expect, readings: dict[str, Any]) -> dict[str, Any]:
    """Albert's three classes (plus order) against what the layers say, by a mapping fixed in advance.

    syntax   the window is constrained across mammals (a fifth of its bases or more)
    storage  a position where the people we hold differ sits on constrained sequence: a value in
             syntax, the reading of attribution/syntax.py
    program  a derived layer names a target gene
    order    the published targets sit in one node holding three or more of them (the beta-globin
             case), or the node boundary falls between two of them (the HOXD case)

    Nothing may choose the mapping after seeing the answer, and the value-domain reading is reported
    with the caveat it earns: common variation is everywhere, so a density of common variable
    positions separates a hypervariable locus from an ordinary one and cannot separate a two-value
    slot from a few-value one.
    """
    c = readings.get("constraint") or {}
    mammals = ((c.get("mammals") or {}).get("fraction_above")) or 0.0
    humans = (c.get("humans") or {}).get("fraction_above")
    sv = readings.get("syntax_values") or {}
    per_kb = (readings.get("frequencies") or {}).get("per_kb")
    node = readings.get("node") or {}
    order = readings.get("order") or {}
    in_node = [g for g in expect.targets if g in (node.get("genes_in_node") or [])]
    got: list[str] = []
    if mammals >= 0.20:
        got.append("syntax")
    if sv.get("values_in_syntax"):
        got.append("storage")
    if (readings.get("deletion") or {}).get("target") or (readings.get("eqtl") or {}).get("target"):
        got.append("program")
    if len(in_node) >= 3 or order.get("node_split_between"):
        got.append("order")
    domain = None
    if per_kb is not None:
        domain = "many" if per_kb >= DOMAIN_MANY_PER_KB else ("few" if per_kb >= DOMAIN_FEW_PER_KB else "two")
    return {
        "field": "class",
        "expected": list(expect.classes),
        "read": got,
        "agrees_on": sorted(set(expect.classes) & set(got)),
        "missed": sorted(set(expect.classes) - set(got)),
        "mammal_fraction": round(mammals, 4),
        "human_fraction": humans,
        "case": c.get("case"),
        "values_in_syntax": sv.get("values_in_syntax"),
        "published_targets_in_one_node": in_node,
        "value_domain_expected": expect.value_domain,
        "value_domain_read": domain,
        "common_variable_positions_per_kb": per_kb,
        "value_domain_hit": (domain == expect.value_domain) if expect.value_domain and domain else None,
        "hit_derived": bool(set(expect.classes) & set(got)),
        "hit_derived_all": bool(set(expect.classes) <= set(got)),
    }


# ----------------------------------------------------------------------------------- the run
#: what the project would need before each mechanism outside the reference sequence could be read
BEYOND_SEQUENCE_NEEDS = {
    "parent of origin": (
        "an allele-resolved methylation layer, or phased reads with a parental assignment: the two"
        " alleles are the same sequence and behave oppositely, so nothing read off the reference can"
        " separate them"
    ),
    "copy number": (
        "a dosage axis: read depth or an assembly-based copy count. Every layer here asks what a"
        " sequence says, never how many times it is present"
    ),
    "developmental stage": (
        "a time axis. No layer carries when it is read, so a fetal-to-adult switch is unmeasurable"
    ),
}


def read_beyond_sequence(expect: Expect) -> dict[str, Any]:
    """The published mechanisms this locus turns on that are not in the reference sequence at all.

    Recorded as pending with what it would take, never scored as a miss: a locus the project cannot
    reach in principle is a different thing from a locus it reaches and gets wrong, and the panel is
    worth less if the two are mixed.
    """
    mechanisms = list(expect.beyond_sequence)
    return {
        "layer": "beyond_sequence",
        "provenance": "derived",
        "mechanisms": mechanisms,
        "needs": {m: BEYOND_SEQUENCE_NEEDS.get(m, "unrecorded") for m in mechanisms},
        "pending": None
        if not mechanisms
        else (
            f"{', '.join(mechanisms)}: published here and carried by no layer in the project, so this"
            " locus is scored on its target and its class only"
        ),
        "evidence": "curated: the published mechanism, written down before the run",
    }


LOCAL_LAYERS = ("registry", "node", "deletion", "reader", "mpra", "lookups")


def local_readings(ch: Chromosome, start: int, end: int, results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """Every reading that needs nothing but the machine: the same call for a locus and a negative."""
    return {
        "registry": read_registry(ch, start, end),
        "node": read_node(ch, start, end),
        "deletion": read_deletion(ch, start, end, results_dir),
        "reader": read_reader(ch, start, end),
        "mpra": read_mpra(ch, start, end),
        "lookups": read_lookups(ch, start, end, results_dir),
    }


def coding_target_named(readings: dict[str, Any]) -> bool:
    """Whether a derived layer named any gene at all: the claim a negative must not make as readily."""
    return bool((readings.get("deletion") or {}).get("target") or (readings.get("eqtl") or {}).get("target"))


def direction_named(readings: dict[str, Any]) -> bool:
    return bool((readings.get("deletion") or {}).get("action"))


def cell_named(readings: dict[str, Any]) -> bool:
    return bool((readings.get("reader") or {}).get("cells_open"))


def claims(readings: dict[str, Any]) -> dict[str, bool]:
    """The three claims scored identically at positives and negatives, plus the two separate layers."""
    return {
        "target": coding_target_named(readings),
        "cell": cell_named(readings),
        "direction": direction_named(readings),
        "storage": bool((readings.get("syntax_values") or {}).get("values_in_syntax")),
        "deletion_scored": bool((readings.get("deletion") or {}).get("elements_scored")),
        "eqtl_present": bool((readings.get("eqtl") or {}).get("eqtls")),
    }


def rate(rows: list[dict[str, bool]], key: str) -> dict[str, Any]:
    n = len(rows)
    k = sum(1 for r in rows if r.get(key))
    return {"k": k, "n": n, "rate": round(k / n, 3) if n else None}


def score_locus(expect: Expect, ch: Chromosome, readings: dict[str, Any]) -> dict[str, Any]:
    scored = {
        "target": score_target(expect, readings),
        "cell": score_cell(expect, readings),
        "direction": score_direction(expect, readings),
        "distance": score_distance(expect, ch, readings),
        "variant": score_variant(expect, readings),
        "class": score_class(expect, readings),
    }
    pending = sorted(
        f"{layer}: {r['pending']}"
        for layer, r in readings.items()
        if isinstance(r, dict) and r.get("pending")
    )
    t = scored["target"]
    trap = expect.nearest_gene_trap
    node_target = (readings.get("node") or {}).get("target")
    return {
        "scored": scored,
        "target_hit_derived": t["hit_derived"],
        "target_hit_heuristic": t["hit_heuristic"],
        "target_hit_looked_up": t["hit_looked_up"],
        "nearest_gene_trap": trap,
        "heuristic_fell_in_trap": bool(trap and node_target == trap),
        "cell_hit_derived": scored["cell"]["hit_derived"],
        "direction_hit_derived": scored["direction"]["hit_derived"],
        "class_hit_derived": scored["class"]["hit_derived"],
        "reachable_by_a_derived_target_layer": bool(
            (readings.get("deletion") or {}).get("elements_scored")
            or (readings.get("eqtl") or {}).get("eqtls")
        ),
        # a different sense of reach from the line above, and the names are kept apart on purpose:
        # that one asks whether any data exists at this window, this one whether the published target
        # is inside the model's 1 Mb input at all. A locus can have plenty of data and still be a
        # question the deletion layer could not have answered.
        "deletion_unaskable": (readings.get("deletion") or {}).get("unaskable") or None,
        "pending": pending,
    }


def aggregate(loci: list[dict[str, Any]]) -> dict[str, Any]:
    """Hit rates by provenance, with the misses named; the reachable denominators kept beside them."""

    def hits(key: str, subset: list[dict[str, Any]] | None = None) -> dict[str, Any]:
        rows = subset if subset is not None else loci
        got = [r["locus"] for r in rows if r["score"][key]]
        return {
            "k": len(got),
            "n": len(rows),
            "rate": round(len(got) / len(rows), 3) if rows else None,
            "hits": got,
            "misses": [r["locus"] for r in rows if not r["score"][key]],
        }

    enhancers = [r for r in loci if r["expected"]["direction"] in ("activates", "represses")]
    reachable = [r for r in loci if r["score"]["reachable_by_a_derived_target_layer"]]
    judged_cells = [r for r in loci if r["expected"]["cells"] or r["expected"]["gtex_tissues"]]
    judged_dir = [r for r in enhancers if r["score"]["scored"]["direction"]["judged"]]
    traps = [r for r in loci if r["expected"]["nearest_gene_trap"]]
    chance = [
        min(1.0, len(r["expected"]["targets"]) / r["coding_genes_in_window"])
        for r in loci
        if r.get("coding_genes_in_window")
    ]
    return {
        "loci": len(loci),
        "target_derived": hits("target_hit_derived"),
        "target_by_chance": {
            "expected": round(sum(chance), 2),
            "of": len(chance),
            "reading": (
                "naming one of the published targets by drawing a coding gene at random from the locus"
                " window: the floor any target rate has to clear"
            ),
        },
        "target_derived_where_reachable": hits("target_hit_derived", reachable),
        "target_derived_where_the_model_could_answer": {
            **hits("target_hit_derived", [r for r in loci if not r["score"]["deletion_unaskable"]]),
            "held_out": [r["locus"] for r in loci if r["score"]["deletion_unaskable"]],
            "reading": (
                "the same rate with the loci whose published target lies outside the scorer's 1 Mb"
                " input held out. It is reported beside the headline and does not replace it: the"
                " headline counts every locus the panel chose, and choosing a target the model"
                " cannot see is the panel's decision to own, not the model's failure"
            ),
        },
        "target_heuristic": hits("target_hit_heuristic"),
        "target_looked_up": hits("target_hit_looked_up"),
        "target_only_looked_up": [
            r["locus"]
            for r in loci
            if r["score"]["target_hit_looked_up"] and not r["score"]["target_hit_derived"]
        ],
        "nearest_gene_traps": {
            "loci": [r["locus"] for r in traps],
            "heuristic_fell_in": [r["locus"] for r in traps if r["score"]["heuristic_fell_in_trap"]],
        },
        "cell_derived_where_judged": hits("cell_hit_derived", judged_cells),
        "cell_unreachable": [r["locus"] for r in loci if r not in judged_cells],
        "direction_derived_where_judged": hits("direction_hit_derived", judged_dir),
        "direction_not_judged": [r["locus"] for r in enhancers if r not in judged_dir],
        "class_derived": hits("class_hit_derived"),
        "unreachable_target": [r["locus"] for r in loci if r not in reachable],
        "published_element_coverage": published_element_coverage(loci),
        "target_derived_by_element_source": by_element_source(loci),
    }


def by_element_source(loci: list[dict[str, Any]]) -> dict[str, Any]:
    """The derived target rate split by whether anybody else drew the element we deleted.

    A hit on an interval the panel stated itself is a different kind of claim from a hit on an
    ENCODE cCRE: the first says the machinery found the gene, the second says it found the gene at a
    place the annotation had already picked out. Pooling them hides the second half of that, so the
    split is reported and the pooled rate is only ever quoted beside it.
    """
    out: dict[str, Any] = {
        "reading": "never quote a pooled derived-target rate without this split: a stated interval is"
        " the panel's own claim about where the element is, an annotated one is somebody else's",
    }
    for source in ("annotated", "stated"):
        mine = [r for r in loci if r["expected"].get("element_source", "annotated") == source]
        hits = [r["locus"] for r in mine if r["score"]["target_hit_derived"]]
        out[source] = {
            "k": len(hits),
            "n": len(mine),
            "rate": round(len(hits) / len(mine), 3) if mine else None,
            "loci": [r["locus"] for r in mine],
            "hits": hits,
            "scored_at_the_element": [
                r["locus"] for r in mine if (r["readings"].get("deletion") or {}).get("elements_scored")
            ],
        }
    return out


def published_element_coverage(loci: list[dict[str, Any]]) -> dict[str, Any]:
    """Whether any annotation calls an interval over the published element, and who supplied it.

    The element-level reading can only ask about intervals somebody else drew: the ENCODE cCRE
    registry, VISTA, lentiMPRA. Where none of them drew one over the published element there is
    nothing to delete, and no amount of sweeping will change that - the bound is the annotation's
    universe, not the model's coverage. Recorded per locus so the two are never confused.
    """
    rows = {}
    for r in loci:
        dl = r["readings"].get("deletion") or {}
        rows[r["locus"]] = {
            "ccres_over_element": r.get("ccres_over_element"),
            "elements_scored": dl.get("elements_scored", 0),
            "runs": dl.get("runs") or [],
        }
    uncalled = [k for k, v in rows.items() if not v["ccres_over_element"]]
    unreadable = [k for k, v in rows.items() if not v["elements_scored"]]
    return {
        "per_locus": rows,
        "no_ccre_over_the_published_element": uncalled,
        "no_element_level_reading_at_all": unreadable,
        "reading": (
            "a locus in the second list cannot get an element-level verdict by sweeping harder: no"
            " annotation drew an interval over its published element, so there is nothing to delete."
            " The fix is a project-defined interval, not more requests"
        ),
    }


def swept_chromosomes(results_dir: Path = RESULTS_DIR) -> list[str]:
    """The chromosomes the all-element sweep has finished, so their windows are fully observed."""
    done = []
    for p in sorted(results_dir.glob("enhancer_targets_all_chr*.json")):
        r = load_result(p.stem, results_dir) or {}
        if r.get("complete") and r.get("chrom"):
            done.append(r["chrom"])
    return sorted(done)


def negative_summary(
    positives: list[dict[str, Any]], negatives: list[Window], results_dir: Path = RESULTS_DIR
) -> dict[str, Any]:
    pos = [claims(r["readings"]) for r in positives]
    neg = [claims(w.readings) for w in negatives]
    keys = ("target", "cell", "direction", "storage", "deletion_scored", "eqtl_present")
    by_locus = {}
    for r in positives:
        mine = [claims(w.readings) for w in negatives if w.of_locus == r["locus"]]
        by_locus[r["locus"]] = {
            "positive": claims(r["readings"]),
            "negatives": len(mine),
            **{k: rate(mine, k)["k"] for k in keys},
        }
    return {
        "windows": len(negatives),
        "per_locus": NEGATIVES_PER_LOCUS,
        "matched_on": [
            "length (exact)",
            f"GC within {MATCH_GC}",
            f"distance to the nearest coding TSS within {int(MATCH_DISTANCE * 100)}%",
            "Zoonomia constrained fraction (closest of the candidates)",
        ],
        "excluded": "every panel locus with 200 kb of flank, VISTA elements and GWAS Catalog hits (5 kb)",
        "positives": {k: rate(pos, k) for k in keys},
        "negatives": {k: rate(neg, k) for k in keys},
        "given_deletion_data": conditioned_on_deletion_data(pos, neg),
        "on_swept_chromosomes": on_swept_chromosomes(negatives, results_dir),
        "by_locus": by_locus,
    }


def on_swept_chromosomes(negatives: list[Window], results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """The controls whose chromosome the all-element sweep has finished: coverage is not the variable.

    A control window on an unfinished chromosome can only be under-read, so the rates over the whole
    sixty move as the sweep advances and cannot be compared with the panel. On a finished chromosome
    the window has been seen in full, and about half of them turn out to hold no registry element at
    all - which is the ceiling on every claim only a deletion can make.
    """
    done = set(swept_chromosomes(results_dir))
    rows = [claims(w.readings) for w in negatives if w.chrom in done]
    keys = ("target", "cell", "direction", "storage", "deletion_scored", "eqtl_present")
    return {
        "chromosomes": sorted(done & {w.chrom for w in negatives}),
        "windows": len(rows),
        "rates": {k: rate(rows, k) for k in keys},
        "reading": "matched negatives whose chromosome the sweep has finished: the only controls that"
        " can be compared with the panel on any claim a deletion makes",
    }


def conditioned_on_deletion_data(pos: list[dict[str, bool]], neg: list[dict[str, bool]]) -> dict[str, Any]:
    """The same claims counted only where an element inside the window has actually been deleted.

    The panel's windows were scored on purpose (`scripts/loci_score.py`) and the matched negatives
    were not, so an unconditional rate for any claim a deletion produces - a direction above all -
    measures which windows we chose to spend requests on rather than anything about the sequence.
    Reported beside the unconditional rates so the denominator travels with the number.
    """
    keys = ("target", "cell", "direction", "storage")
    out: dict[str, Any] = {
        "reading": "the same claims at the windows where an element-level deletion exists, and at"
        " the windows where none does: a claim that only a deletion can make is worth nothing until"
        " both sides have been asked",
    }
    for name, rows in (("positives", pos), ("negatives", neg)):
        have = [r for r in rows if r["deletion_scored"]]
        lack = [r for r in rows if not r["deletion_scored"]]
        out[name] = {
            "with_deletion_data": {k: rate(have, k) for k in keys},
            "without_deletion_data": {k: rate(lack, k) for k in keys},
        }
    return out


def build(
    panel: tuple[Expect, ...] = PANEL,
    results_dir: Path = RESULTS_DIR,
    gtex_dir: Path = GTEX_DIR,
    network: bool = True,
    negatives: bool = True,
    progress=None,
) -> dict[str, Any]:
    """Read every locus and its matched negatives through the same layers, then score them."""
    t0 = time.time()
    say = progress or (lambda _m: None)
    chroms: dict[str, Chromosome] = {}
    loci: list[dict[str, Any]] = []
    windows: list[Window] = []
    try:
        for e in panel:
            ch = chroms.get(e.chrom) or chroms.setdefault(e.chrom, Chromosome(e.chrom, results_dir))
            s, en = e.element
            readings = local_readings(ch, s, en, results_dir)
            readings["reach"] = read_reach(ch, e)
            # an out-of-reach target is not a pending request: spending one would buy a certain
            # negative, so the deletion layer says unaskable and the cost line is withdrawn. This
            # applies whether or not the layer has a reading: the ZRS HAS one (LMBR1) and it is
            # still an answer to a question the model was never in a position to be asked.
            if not readings["reach"]["askable"]:
                far = readings["reach"]["targets_out_of_reach"]
                why = (
                    "unaskable rather than unasked: "
                    + "; ".join(
                        f"{x['target']} is {x.get('distance_bp', 0) // 1000} kb from the element" for x in far
                    )
                    + f", and the model's input reaches {MODEL_WINDOW // 2 // 1000} kb each way."
                    " No number of requests can name it"
                )
                readings["deletion"]["unaskable"] = why
                if readings["deletion"].get("pending"):
                    readings["deletion"]["pending"] = why
            readings["satmut"] = read_satmut(ch, s, en, e)
            readings["gene_input"] = read_gene_input(ch, e, results_dir)
            readings["beyond_sequence"] = read_beyond_sequence(e)
            if e.direction == "coding":
                readings["coding"] = read_coding(ch, e)
                tss = ch.tss(e.targets[0])
                if tss is not None:
                    readings["promoter"] = read_reader(ch, tss - 1_000, tss + 1_000)
            if e.locus == "HBB_LCR":
                readings["ladder"] = read_hypersensitive_ladder(ch, s, en)
            if e.locus == "HOXD":
                readings["order"] = read_cluster_order(ch, e)
            gc = ch.gc(s, en)
            near = ch.nearest_coding((s + en) // 2)
            row = {"locus": e.locus, "expected": e.as_dict(), "readings": readings, "gc": gc}
            row["distance_to_coding_tss"] = near[1] if near else None
            row["ccres_over_element"] = len(ch.ccres_in(s, en))
            row["coding_genes_in_window"] = sum(1 for t, _g in ch.coding if e.window[0] <= t < e.window[1])
            loci.append(row)
            if negatives and gc is not None and near is not None:
                n = CANDIDATES_PER_LOCUS if en - s <= 20_000 else CANDIDATES_PER_LOCUS // 2
                cands = candidate_windows(ch, e, gc, near[1], n=n)
                windows.extend(cands)
                say(f"{e.locus}: local layers read, {len(cands)} candidate negatives")

        if network:  # both constraint axes, one range pass per chromosome over positives and candidates
            for chrom in chroms:
                pos = [r for r in loci if r["expected"]["chrom"] == chrom]
                cand = [w for w in windows if w.chrom == chrom]
                ivs = non_overlapping(
                    [tuple(r["expected"]["element"]) for r in pos] + [(w.start, w.end) for w in cand]
                )
                read = _retry(read_constraint, chrom, ivs, say=say)
                got = dict(zip(ivs, read, strict=True)) if read is not None else {}
                if read is None:
                    say(f"{chrom}: constraint unavailable after {RETRIES} tries")
                for r in pos:
                    r["readings"]["constraint"] = got.get(tuple(r["expected"]["element"])) or {
                        "layer": "constraint",
                        "provenance": "derived",
                        "pending": "the range read did not return this window",
                    }
                for w in cand:
                    c = got.get((w.start, w.end))
                    if c:
                        w.constrained = c["mammals"]["fraction_above"]
                        w.readings["constraint"] = c
                say(f"{chrom}: constraint read over {len(ivs)} windows")

        chosen: list[Window] = []
        for r in loci:
            mine = [w for w in windows if w.of_locus == r["locus"]]
            pc = ((r["readings"].get("constraint") or {}).get("mammals") or {}).get("fraction_above")
            r["constrained_fraction"] = pc
            chosen.extend(pick_negatives(mine, pc))
        for w in chosen:
            w.readings.update(local_readings(chroms[w.chrom], w.start, w.end, results_dir))

        gtex_summary = None
        if network:
            spans = [(r["expected"]["chrom"], *r["expected"]["element"]) for r in loci]
            spans += [(w.chrom, w.start, w.end) for w in chosen]
            gtex_summary = stream_gtex(spans, gtex_dir, progress=say)
        for r in loci:
            s, en = r["expected"]["element"]
            r["readings"]["eqtl"] = read_eqtl(chroms[r["expected"]["chrom"]], s, en, knowledge=gtex_dir)
        for w in chosen:
            w.readings["eqtl"] = read_eqtl(chroms[w.chrom], w.start, w.end, knowledge=gtex_dir)

        checked: dict[str, Any] = {}
        resolved_panel: list[Expect] = []
        for r, e in zip(loci, panel, strict=True):
            freqs: dict[str, Any] = {}
            if network:
                for v in e.variants:
                    if v.get("rsid"):
                        checked[v["rsid"]] = position_check(v, ensembl_variant(v["rsid"]))
                was = tuple(v.get("pos") for v in e.variants)
                e = resolve_positions(e, checked)
                r["expected"] = e.as_dict()
                if e.direction == "coding" and was != tuple(v.get("pos") for v in e.variants):
                    # the coding layer ran with the local readings, before Ensembl had been asked:
                    # a locus whose positions the panel does not hard-code needs it read again
                    r["readings"]["coding"] = read_coding(chroms[e.chrom], e)
                freqs = _guarded(read_frequencies, chroms[e.chrom], *e.element, e)
                r["readings"]["frequencies"] = freqs
            resolved_panel.append(e)
            r["readings"]["values"] = read_values(chroms[e.chrom], e, freqs)
        panel = tuple(resolved_panel)
        if network:  # syntax against values, one bigWig session per chromosome, positives and negatives
            for chrom, ch in chroms.items():
                mine = [r for r in loci if r["expected"]["chrom"] == chrom]
                negs = [w for w in chosen if w.chrom == chrom]
                wins = [tuple(r["expected"]["element"]) for r in mine] + [(w.start, w.end) for w in negs]
                by_name = {e.locus: e for e in panel}
                exps = [by_name[r["locus"]] for r in mine] + [None] * len(negs)
                got_sv = _retry(read_syntax_values_many, ch, wins, exps, say=say)
                if got_sv is None:
                    say(f"{chrom}: syntax against values unavailable after {RETRIES} tries")
                    continue
                for r, sv in zip(mine, got_sv[: len(mine)], strict=True):
                    r["readings"]["syntax_values"] = sv
                for w, sv in zip(negs, got_sv[len(mine) :], strict=True):
                    w.readings["syntax_values"] = sv
                say(f"{chrom}: syntax against values over {len(wins)} windows")
                negs_f = [(w, _guarded(read_frequencies, ch, w.start, w.end)) for w in negs]
                for w, fr in negs_f:
                    w.readings["frequencies"] = fr
        for r, e in zip(loci, panel, strict=True):
            r["score"] = score_locus(e, chroms[e.chrom], r["readings"])
            say(f"{e.locus}: scored")
    finally:
        for ch in chroms.values():
            ch.close()

    return {
        "panel": len(panel),
        "classes": list(CLASSES),
        "provenance": list(PROVENANCE),
        "aggregate": aggregate(loci),
        "negative_controls": negative_summary(loci, chosen, results_dir),
        "loci": [trim_locus(r) for r in loci],
        "negatives": [negative_row(w) for w in chosen],
        "variant_positions": checked,
        "not_in_panel": NOT_IN_PANEL,
        "gtex": gtex_summary,
        "evidence": EVIDENCE,
        "network": network,
        "note": NOTE,
        "cost": {"seconds": round(time.time() - t0, 1)},
    }


NOTE = (
    "A coherence check, not a discovery claim. Hits are split by provenance: derived (a genome-wide"
    " measurement or model read blind), targeted (saturation mutagenesis made at these loci because they"
    " were known), heuristic (nearest coding TSS in the node) and looked_up (an annotation that already"
    " carries the answer). Negative controls are matched windows scored through the same readers. No"
    " AlphaGenome request was made: loci without an already computed deletion are pending, not dropped."
)


RETRIES = 3  # a range read over a public track times out now and then; the gate catches a silent gap


def _retry(fn, *args, say=None):
    """Try a range read up to RETRIES times; None when every try failed."""
    for n in range(RETRIES):
        try:
            return fn(*args)
        except OSError as ex:
            if say:
                say(f"try {n + 1} of {RETRIES} failed ({str(ex)[:60]})")
            time.sleep(5 * (n + 1))
    return None


def _guarded(fn, *args) -> dict[str, Any]:
    """A network reading that fails is recorded as pending, never allowed to stop the benchmark."""
    try:
        return fn(*args)
    except (OSError, ValueError, KeyError) as ex:
        return {"layer": fn.__name__.removeprefix("read_"), "provenance": "derived", "pending": str(ex)[:160]}


def non_overlapping(ivs: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """The windows in order with any that overlap an earlier one removed (the range reader needs it)."""
    out: list[tuple[int, int]] = []
    for s, e in sorted(set(ivs)):
        if out and s < out[-1][1]:
            continue
        out.append((s, e))
    return out


def resolve_positions(e: Expect, checked: dict[str, Any]) -> Expect:
    """Fill in a variant the panel deliberately did not hard-code, from Ensembl's answer.

    A locus added later should not need a coordinate typed from memory: `pos: None` says the panel is
    not claiming one, and the position comes from Ensembl at run time. The cost is that the
    position check is not independent for that variant - `agrees` stays null and the result records
    `position_from: ensembl` - so a hard-coded position is still better where one can be cited.
    """
    from dataclasses import replace

    out = []
    changed = False
    for v in e.variants:
        got = checked.get(v.get("rsid") or "")
        if v.get("pos") is None and got and got.get("ensembl"):
            out.append({**v, "pos": got["ensembl"], "position_from": "ensembl"})
            changed = True
        else:
            out.append(dict(v))
    return replace(e, variants=tuple(out)) if changed else e


def position_check(v: dict[str, Any], got: dict[str, Any] | None) -> dict[str, Any]:
    """The panel's hard-coded position against Ensembl's; an indel may differ by its anchor base."""
    if got is None or got.get("position") is None:
        return {"panel": v.get("pos"), "ensembl": None, "agrees": None}
    if v.get("pos") is None:
        # the panel deliberately does not claim a position for this variant: Ensembl supplies it and
        # there is nothing to check. Recording False here would read as a disagreement that never was.
        return {"panel": None, "ensembl": got["position"], "agrees": None, "position_from": "ensembl"}
    return {
        "panel": v["pos"],
        "ensembl": got["position"],
        "agrees": abs(got["position"] - v["pos"]) <= 1,
    }


def stream_gtex(spans: list[tuple[str, int, int]], gtex_dir: Path, pad: int = 500, progress=None):
    """One pass over GTEx's archive keeping the pairs inside the benchmark's windows.

    The hits are git-ignored under data/knowledge/loci_benchmark. A manifest of the windows is kept
    beside them and the directory is cleared when the windows change, because the distiller skips a
    tissue whose file exists and would otherwise answer for a different set of windows.
    """
    from genomeos.attribution.eqtl import distil

    padded = sorted((c, max(0, s - pad), e + pad) for c, s, e in spans)
    manifest = gtex_dir / "windows.json"
    if manifest.exists() and json.loads(manifest.read_text()) != [list(x) for x in padded]:
        for p in gtex_dir.glob("hits_*.tsv"):
            p.unlink()
    gtex_dir.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps([list(x) for x in padded]))
    summary = distil(gtex_intervals(padded), knowledge=gtex_dir, progress=progress)
    _GTEX.clear()
    return {k: summary[k] for k in ("tissues", "pairs_scanned_this_run", "hits_this_run", "seconds")}


def trim_locus(r: dict[str, Any]) -> dict[str, Any]:
    """The committed row: every reading kept, the long lists shortened."""
    rd = dict(r["readings"])
    if rd.get("lookups"):
        rd["lookups"] = {**rd["lookups"], "gwas": rd["lookups"]["gwas"][:5]}
    return {**r, "readings": rd}


def negative_row(w: Window) -> dict[str, Any]:
    """A negative as committed: where it is, what it was matched on, and the claims made at it."""
    rd = w.readings
    return {
        **{k: v for k, v in w.as_dict().items() if k != "readings"},
        "claims": claims(rd),
        "node_target": (rd.get("node") or {}).get("target"),
        "deletion_target": (rd.get("deletion") or {}).get("target"),
        "eqtl_target": (rd.get("eqtl") or {}).get("target"),
        "cells_open": (rd.get("reader") or {}).get("cells_open"),
        "registry_classes": (rd.get("registry") or {}).get("classes"),
        "common_variants_per_kb": (rd.get("frequencies") or {}).get("per_kb"),
        "values_in_syntax": (rd.get("syntax_values") or {}).get("values_in_syntax"),
        "mammal_fraction": ((rd.get("constraint") or {}).get("mammals") or {}).get("fraction_above"),
    }


def run_and_save(
    results_dir: Path = RESULTS_DIR, network: bool = True, progress=None, gate_only: bool = False
) -> dict[str, Any]:
    """The full panel by default; `gate_only` runs the twelve pinned loci for a quick check.

    The gate set is not a different benchmark, only a cheaper pass over the loci whose verdicts are
    already defended, so the two write different files and neither overwrites the other.
    """
    panel = tuple(e for e in PANEL if e.locus in GATE_PANEL) if gate_only else PANEL
    out = build(panel=panel, results_dir=results_dir, network=network, progress=progress)
    save_result("loci_benchmark_gate" if gate_only else "loci_benchmark", out, results_dir)
    return out
