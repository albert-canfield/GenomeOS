# The known-locus benchmark

Albert's brief (2026-09-13): take the few blocks where biology is already
known, their variants, their syntax, their enhancers and their context, and use
them as a coherence check on everything the project builds. If the machinery
cannot reproduce what is settled at a dozen well-studied loci, its genome-wide
numbers are not to be believed.

Code: `genomeos/benchmark/loci.py`, `scripts/loci_benchmark.py`,
`tests/test_loci_benchmark.py`. Result: `data/results/loci_benchmark.json`.
The pattern is the therapeutic benchmark (`tests/test_therapeutic_benchmark.py`).

## 1. The panel

Twelve loci, each chosen because the published answer is unambiguous and
because each tests a different failure mode. The expectations are written down
in the code (`PANEL`) before anything is read, with citations and with the kinds
of data the published answer itself came from, so a reader can see where a
layer might only be restating the discovery.

| locus | class | published target, tissue | causal variant | failure mode it tests |
|---|---|---|---|---|
| SHH_ZRS | program | SHH, limb bud, 979 kb away | many rare point mutations | long-range target naming: the nearest TSS (LMBR1) must fail |
| HERC2_OCA2 | storage, program | OCA2, melanocyte, 21 kb | rs12913832 (A brown, G blue) | a two-value slot Gnocchi does not score |
| MCM6_LCT | storage, program | LCT, small intestine, 14 kb | rs4988235 (T-13910) | allele-dependent activity, population structure |
| HBB_LCR | syntax, order, program | five beta-like globins, erythroid | LCR deletion | several genes in one node, a stage switch, hypersensitive ladder |
| FTO_IRX3 | program | IRX3, IRX5, adipocyte precursors, 520 kb | rs1421085 | the textbook nearest-gene fallacy (FTO) |
| BCL11A_enhancer | program, storage | BCL11A, erythroid, 59 kb | rs1427407 (GATA1 site) | a known causal variant in a cell the reader has (K562) |
| MYC_8q24 | program | MYC, colon, 334 kb | rs6983267 | whether a gene's whole input assembles |
| ABO | storage | ABO | rs8176719 (O), rs8176746 (B) | a few-valued slot |
| HLA_DRB1 | storage | HLA-DRB1, antigen-presenting cells | thousands of alleles | a many-valued slot |
| APP | storage, syntax | APP, brain | rs63750847 (A673T) | coding calibration on a fully scored chromosome |
| TP53 | syntax, storage | TP53, every tissue | rs1042522 (P72R) | a constrained frame with one common value |
| HOXD | order, program, syntax | HOXD1 to HOXD13, limb and trunk | none | position as the order of execution; the boundary between HOXD11 and HOXD13 |

HOXD was added at the coordinator's request (Albert asking in what order DNA is
read): it is the clearest case where position on the chromosome is the order of
execution. With the beta-globin cluster the panel covers ordered in time and
ordered in space.

Considered and left out (`NOT_IN_PANEL`): IRX5 on its own (kept as a second
target of FTO_IRX3); the HS2 core alone (no citable hg38 coordinate without a
liftover, so the whole LCR is the element and its hypersensitive sites are
derived from DNase); HLA-B (one HLA gene answers the value-domain question);
the African lactase alleles (their window overlaps rs4988235's, so the
expectation would not be independent). Nothing on the brief was dropped.

Every hard-coded variant position is checked against Ensembl on each run, and
the two coding positions are checked a second way: the engine's own translation
must derive the published protein change from them.

## 2. Circularity: four provenances

Several answers are already inside the annotations we read, so every hit is
labelled by where it came from, and the rates are reported apart.

| provenance | layers | what a hit proves |
|---|---|---|
| derived | AlphaGenome deletion (already computed), GTEx v8 eQTLs, ENCODE DNase in 11 cell types (the reader), ENCODE4 lentiMPRA, Zoonomia phyloP, gnomAD Gnocchi, the GIAB trio's variants, 1000 Genomes frequencies, the translation engine | the approach works |
| targeted | Kircher et al. 2019 saturation mutagenesis | a measurement made at the locus because it was known: confirms, cannot discover |
| heuristic | nearest coding TSS in the CTCF node | the baseline the rest has to beat |
| looked_up | VISTA tissues, GWAS Catalog mapped gene, ClinVar gene, cCRE class | nothing about the approach |

A layer's target hit is strict: the gene it ranks first must be a published
target. The lenient reading (a published target anywhere in its list) is kept
as `among` and never counted. Two caveats are recorded rather than hidden: the
coding layer derives the protein change, but the gene's name is GENCODE's; and
at FTO and MYC the published answer was itself partly found by eQTLs, so the
eQTL layer there is closer to restating the discovery than elsewhere.

## 3. Negative controls

For every locus, five windows of the same length with no known function,
matched on GC (within 0.04), distance to the nearest coding TSS (within 35%)
and the Zoonomia constrained fraction (closest of 40 candidates, 20 for windows
over 20 kb), drawn from the same chromosome away from every panel locus (200 kb
of flank), VISTA elements and GWAS Catalog hits. They go through the same
readers, and the same four claims are counted at both: a target named by a
derived layer, a cell open in the reader, a direction, a value on syntax.

## 4. The class mapping, fixed in advance

- **syntax**: a fifth or more of the window's bases constrained across mammals;
- **storage**: a position where the people we hold differ sits on constrained
  sequence (a value in syntax, `attribution/syntax.py`'s reading);
- **program**: a derived layer names a target;
- **order**: three or more published targets in one node, or the node boundary
  between two of them.

The value-domain size is read as common variable positions per kilobase
(1000 Genomes MAF of 5% or more through Ensembl).

## 5. The first run (2026-09-14)

12 loci, 60 matched negatives, no AlphaGenome request. GTEx's archive streamed
once (71.5 million pairs scanned, 123,520 kept inside the windows, 132 s),
Zoonomia and Gnocchi read by range over 463 windows, Ensembl asked for the
eleven named variants. Ten of the eleven hard-coded positions agreed with
Ensembl exactly; the eleventh (rs8176746) did not, which is the check working.

| field | derived | heuristic (nearest TSS in node) | looked up |
|---|---|---|---|
| the right target gene | **8/12** | 8/12 | 6/12 |
| the right cell or tissue | 4/10 judged | - | - |
| the direction | 0/1 judged | - | - |
| at least one published class | 10/12 | - | - |

No locus was a looked-up hit without also being a derived one, so the panel is
not measuring annotation lookup. But the derived layers do **not** beat the
nearest-TSS heuristic: both score 8 of 12, on different loci. The derived
misses are HERC2_OCA2, MCM6_LCT, FTO_IRX3 and MYC_8q24, all four programs acting
at a distance or through an allele, all four on chromosomes where no deletion
has been scored, so the only derived target layer there was GTEx, which names
the neighbour (HERC2, MCM6, FTO, CASC8). The heuristic gets OCA2 and LCT right
because they are the nearest genes, and falls into the one trap it was set at
FTO. Of the derived hits, four come from summing every already-scored element in
the locus window (SHH, BCL11A, HLA-DRB1, APP), two from eQTLs (HBG2, ABO), two
from the translation engine (APP, TP53) and one from an element's own deletion
(HOXD10); the only other element-level deletion in the panel, the ZRS itself,
names LMBR1, the textbook wrong answer.

### What the negative controls say, and it is the main finding

| claim | positives | matched negatives |
|---|---|---|
| a derived layer names *some* target | 11/12 (92%) | 52/60 (87%) |
| the reader has it open in some cell | 10/12 (83%) | 32/60 (53%) |
| a direction is produced | 2/12 (17%) | 2/60 (3%) |
| a value sits on constrained sequence | 6/12 (50%) | 9/60 (15%) |

**Naming a target is not evidence of anything.** At windows matched for length,
GC, distance to the nearest coding TSS and mammalian constraint, the pipeline
names a gene almost as often as at the panel, because GTEx ties an eQTL to
something nearly everywhere (51 of 60 negatives carry one). What separates the
panel from its controls is *which* gene is named, whether a direction comes with
it, and whether a variable position sits on syntax. Any genome-wide count of
"elements with a named target" should be read against 87%, not against zero.

### Per locus

| locus | target, derived (layer) | heuristic | cell | causal variant by constraint | notes |
|---|---|---|---|---|---|
| SHH_ZRS | SHH (window input; the element's own deletion says **LMBR1**) | none in node | unreachable (no limb cell) | - | the flagship long-range case fails at the element and passes only at window level |
| HERC2_OCA2 | miss (eQTLs name HERC2) | OCA2 | miss | rs12913832 rank 1 of 4, phyloP 3.41 | no deletion scored on chr15; the 957 bp window reads 100% human-constrained because it touches scored neighbouring kilobases, while rs12913832's own kilobase is unscored |
| MCM6_LCT | miss (eQTLs name MCM6, UBXN4) | LCT | unreachable in the reader; GTEx intestine silent for LCT | rs4988235 rank 1 of 2, phyloP 0.22 (not constrained) | no deletion scored on chr2 here |
| HBB_LCR | HBG2 (eQTL) | HBE1 | K562 | - | the reader finds 8 strong K562 sites, 5 of them erythroid-specific, median spacing 2.9 kb: the hypersensitive ladder, derived; the LCR is only 2.5% mammal-constrained, so "syntax" is not read |
| FTO_IRX3 | strict miss: the eQTLs rank FTO first and IRX3 second, and IRX3 sits at 519.7 kb against the published 520 kb | **FTO, the trap** | unreachable in the reader; the only GTEx tissue for IRX3 here is pancreas | - | the lenient reading is right and the strict one is wrong: the answer is in the layer, one rank down |
| BCL11A_enhancer | BCL11A (window input) | BCL11A | K562 not open over the +62 DHS | rs1427407 rank 1 of 12, phyloP 7.01 | the sharpest causal-variant result of the panel |
| MYC_8q24 | miss: no element on chr8 is scored, and the eQTLs name CASC8 and POU5F1B | none in node | miss | rs6983267 rank 1 of 2, phyloP 4.97 | saturation mutagenesis calls rs6983267 **not** functional in its element (effect 0.0, 41st by effect), a measured disagreement with the published colorectal result |
| ABO | ABO (eQTL, 43 tissues) | ABO | miss: the eQTL tissues are 43 of 49, not blood | rs8176746 16th of 20 by constraint (phyloP -2.27) | the O frameshift is not derivable: the local trace substitutes single bases only |
| HLA_DRB1 | HLA-DRB1 (window input) | HLA-DRB1 | GM12878 and monocyte open | - | no value-domain reading came back, so the many-valued slot is unmeasured |
| APP | APP (window input, 100 of 197 scored elements; and p.Ala673Thr from the translation engine) | none (the variant is 273 kb inside a minus-strand gene) | promoter open in SK-N-SH and astrocyte | - | the coding calibration passes exactly |
| TP53 | TP53 (p.Pro72Arg derived from the reference sequence) | TP53 | K562 | rs1042522 rank 1 of 1, phyloP 3.40 | 23% mammal-constrained, read as syntax |
| HOXD | HOXD10 (element deletion, 2 sampled elements) | HOXD4 | unreachable (no limb cell) | - | see below |

### The two ordered loci

**Beta-globin, ordered in time.** The reader recovers the locus control region
without being told it is there: 8 K562 DNase sites above signal 100 inside
chr11:5,269,000-5,295,000, five of them open in at most one of the other ten
cell types, median spacing 2.9 kb. That is the published HS ladder, derived from
measurement alone. The node holds HBB, HBD and HBG1 together (HBG2 and HBE1 fall in the next node), so
"several genes under one control region" is recovered as structure. The switch
itself is not: nothing in the project carries developmental stage, so fetal
against adult is unmeasurable here and is recorded as pending.

**HOXD, ordered in space.** The node model does not recover the cluster's
architecture. There is no CTCF-only element anywhere inside the cluster in the
registry, so no inferred boundary falls in the published HOXD11-to-HOXD13
interval, and all nine genes sit in one 300 kb node (chr2:D1195) instead of the
two regulatory landscapes the 4C and Hi-C work describes. What is recovered is
that the cluster is unusual sequence: 31% mammal-constrained and 29%
human-constrained (the "syntax" case on both axes, the only panel locus that
reads so), 26 variable positions of the trio sitting on constrained bases, and
a deletion-scored element naming HOXD10. The order of activation is not
testable with anything we hold - there is no time axis and no position along the
body axis in any layer - and the language has no construct for order. That is
the honest answer to the question the locus was added for.

### What could not be reached, and why

| locus | what is missing | cost to fix |
|---|---|---|
| HERC2_OCA2, MCM6_LCT, FTO_IRX3, MYC_8q24, BCL11A, ABO, HLA_DRB1, HBB_LCR | no AlphaGenome deletion has been scored on chr15, chr2, chr16, chr8, chr9, chr6 or chr11 | the sweep has run on chr21, chr22, chrY, chr19, chr20, chr18 and is running on chr17 and chr1; about 0.76 requests per element |
| SHH_ZRS, MCM6_LCT, FTO_IRX3, HOXD | the published cell type (limb bud, intestinal epithelium, adipocyte precursor) is not among the eleven reader cell types | ENCODE DNase for the tissue, or an embryonic panel |
| HBB_LCR | developmental stage: the fetal-to-adult switch | no layer carries time |
| ABO | the O allele is an indel; the local trace substitutes single bases | an indel path in the trace |
| HLA_DRB1, all storage loci | the value-domain size: the Ensembl query returned no allele frequencies, so "a handful against thousands" is unmeasured | the 1000 Genomes common-variant set, or the HPRC panel's value domains (genomeos-h1) |
| every locus | the direction of effect is only judged where a deletion names the target: one locus of twelve | the same sweep |

### The one result that is both derived and sharp

Constraint across mammals, read per base over the positions where the three
people we hold differ, puts the published causal variant **first** at four of
the five loci that have one: rs12913832 (phyloP 3.41, first of 4 variable
positions), rs1427407 (7.01, first of 12), rs6983267 (4.97, first of 2),
rs1042522 (3.40, the only one), and rs4988235 first of 2 at phyloP 0.22, which
is not constrained at all. Read honestly, only BCL11A's rs1427407 is a result
(first of twelve on strongly constrained sequence); where a window holds one or
two variable positions, being first means nothing. The claim that does survive
the controls is the class: a value sitting on constrained sequence at 6 of 12
panel loci against 9 of 60 matched windows.

## 6. Lessons the run produced

- **Name a target and you have said nothing**: 87% of matched negative windows
  get one. Only the identity of the gene, a direction, or a value on syntax
  separates the panel from its controls.
- **GENCODE gene spans are not transcription starts.** HBG2 and HBE1 span to
  chr11:5.5 Mb through read-through transcripts, which puts two globin genes a
  quarter of a megabase from their promoters. The benchmark takes the canonical
  transcript's start (`tss_of`); a gene-span TSS would have made the globin
  cluster unrecognisable.
- **The strict and the lenient reading disagree exactly where it matters.** At
  FTO the right gene is second in the eQTL list at the published distance. A
  pipeline that reports "the target" hides this; reporting a ranked list with
  distances does not.
- **A hard-coded coordinate is a bug waiting.** Ten of eleven positions checked
  out against Ensembl and one did not; the translation engine then confirmed two
  of them a second way by deriving p.Ala673Thr and p.Pro72Arg.
- **Cost.** The run is 20 minutes of layers and constraint plus 2 minutes of
  GTEx; the first version spent six hours opening one bigWig session per
  negative window for the per-base constraint reading. Batch the windows.