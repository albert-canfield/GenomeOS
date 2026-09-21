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

Seventeen loci (twelve until 2026-09-15; see section 9), each chosen because the
published answer is unambiguous and
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
| SOX9_PierreRobin | program | SOX9, mandibular arch, 1.45 Mb away | translocations and rare point mutations | a **second** megabase reach, so the ZRS is not an anecdote; KCNJ2 is the trap |
| H19_ICR1 | program, syntax | IGF2 and H19, most fetal tissues, 138 kb | ICR microdeletions, loss of methylation | **parent of origin**: the same two alleles behave oppositely |
| PMP22_CMT1A | program | PMP22, Schwann cell | the 1.4 Mb CMT1A duplication | **copy number**, not a base, is the variable |
| CYP2D6 | storage, syntax | CYP2D6, liver | rs3892097 (*4), rs1065852 (*10) | a many-valued slot outside the MHC |
| MC1R | storage, syntax | MC1R, melanocyte | rs1805007 (R151C), rs1805008 (R160W) | a few-valued slot with a sharp population structure |

The last five were added on 2026-09-15, after the saturation finding in section 8 showed the sweep
had nothing left to give twelve loci. They were chosen for what the first twelve underweighted:
one more long-range case, and three mechanisms the project cannot represent at all (parent of
origin, copy number, and — already present at the globin locus — developmental stage). Those
mechanisms are written down per locus in `Expect.beyond_sequence` and read as **pending with the
cost of reaching them**, never as a miss: a locus the machinery cannot reach in principle is a
different thing from one it reaches and gets wrong, and mixing the two makes the panel worth less.

Two coordinates in the new five are **anchored rather than quoted**, and the code says so: the
SOX9 element is placed at the published −1.45 Mb offset from GENCODE's canonical SOX9 TSS, and the
H19 ICR at the published 2-to-4 kb upstream of GENCODE's canonical H19 TSS. Both windows are wide
enough to absorb the uncertainty. Variants whose hg38 position the panel does not hard-code carry
`pos: None` and are resolved from Ensembl at run time, labelled `position_from: ensembl`; the cost
is that the position check is not independent for those, so a citable hard-coded position is still
better where one exists.

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

Dropped from the 2026-09-15 widening, and the reasons are the discipline: **SMN1/SMN2**, the
cleanest dosage phenotype there is, because the two genes are a segmental duplication differing at
a handful of bases, so every layer's answer would be a statement about read mapping rather than
about biology — PMP22 carries dosage instead; **RHD**, because Rh negative is a whole-gene deletion
and asks the same question as PMP22 on a chromosome the sweep has not reached; **FUT2** secretor
status, clean but a two-valued slot the panel already has at HERC2/OCA2; **AMY1** copy number,
because the published association with dietary starch is contested; and the **EPHA4/IHH** TAD
rearrangements, the best evidence anywhere that node boundaries matter, because the causal unit is
a structural variant spanning a boundary and the panel has no way to write one down as an element.
That last one is a gap in the benchmark's own vocabulary, not in the literature.

Every hard-coded variant position is checked against Ensembl on each run, and
four coding positions are checked a second way: the engine's own translation
must derive the published protein change from them (APP's p.Ala673Thr, TP53's
p.Pro72Arg, MC1R's p.Arg151Cys and p.Arg160Trp).

## 2. Circularity: four provenances

Several answers are already inside the annotations we read, so every hit is
labelled by where it came from, and the rates are reported apart.

| provenance | layers | what a hit proves |
|---|---|---|
| derived | AlphaGenome deletion (already computed), GTEx v8 eQTLs, ENCODE DNase in 11 cell types (the reader), ENCODE4 lentiMPRA, Zoonomia phyloP, gnomAD Gnocchi, the GIAB trio's variants, gnomAD frequencies, the translation engine | the approach works |
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
(gnomAD v4.1.1 genomes, allele frequency 5% or more, read by range as a bigBed
through `attribution/human_panel.py`'s reader): **many** at 50 per kb or more,
**few** at 10 or more, **two** below. The two boundaries were set after reading
the storage loci and two background windows,
so the agreement with the published domain sizes is descriptive; the matched
negatives are what test the reading.

## 5. The first run, before the quota was reallocated (2026-09-14)

*Superseded by section 6 for the target and direction verdicts, and kept because the
comparison is the argument for the reallocation: this is what the panel said when only
four of twelve loci sat on a chromosome the sweep had reached. The negative-control rates,
the storage results and the two ordered loci are unchanged by it.*

12 loci, 60 matched negatives, no AlphaGenome request, 21 to 48 minutes
depending on the network. GTEx's
archive streamed once (71.5 million pairs scanned, 123,520 kept inside the
windows, 132 s), Zoonomia phyloP and gnomAD Gnocchi read by range over 463
windows plus one base per variable position, gnomAD frequencies read as a
bigBed, Ensembl asked for the named variants' positions. The position check
earned its place on the first pass: ten of eleven agreed and rs8176746 was
seven bases out in the panel, now corrected; on the second pass Ensembl timed
out for most of them, which is recorded as unverified rather than agreed.

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
| a direction is produced | 2/12 (17%) | 6/60 (10%) |
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
| ABO | ABO (eQTL, 43 tissues) | ABO | miss: the eQTL tissues are 43 of 49, not blood | rs8176746 19th of 20 by constraint (phyloP -2.6), at the corrected position | the O frameshift is not derivable: the local trace substitutes single bases only |
| HLA_DRB1 | HLA-DRB1 (window input) | HLA-DRB1 | GM12878 and monocyte open | - | 112.8 common variable positions per kb against a median 2.5 at the controls: the many-valued slot, derived |
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
| ~~HERC2_OCA2, MCM6_LCT, FTO_IRX3, MYC_8q24, BCL11A, ABO, HLA_DRB1, HBB_LCR~~ | no AlphaGenome deletion had been scored on their chromosomes | **fixed** in section 6: 2,283 requests over the windows alone |
| SHH_ZRS, MCM6_LCT, FTO_IRX3, HOXD | the published cell type (limb bud, intestinal epithelium, adipocyte precursor) is not among the eleven reader cell types | ENCODE DNase for the tissue, or an embryonic panel |
| HBB_LCR | developmental stage: the fetal-to-adult switch | no layer carries time |
| every storage locus | two values against a few: the density of common positions cannot separate them | haplotype value domains (the HPRC panel, genomeos-h1) |
| ABO | the O allele is a deletion: neither the frequency reader nor the trace handles an indel | an indel path in both |
| every locus | the direction of effect is only judged where a deletion names the target: one locus of twelve | the same sweep |

### The chance floor, and why 8 of 12 is not a result on its own

A locus window holds few coding genes: 1 at APP, 2 at MYC, BCL11A and ABO, 3 at FTO,
4 at the ZRS, HERC2 and MCM6, 8 at the globin cluster, HLA and TP53, 12 at HOXD.
Drawing a coding gene at random from each window and asking whether it is one of
the published targets gives an expectation of **5.5 hits of 12**. The derived
layers score 8 and the nearest-TSS heuristic 8, so neither is far above the
floor, and the panel is too small for that difference to mean anything. The
result to carry away is not the rate; it is which loci fail and how.

### The storage side: values, frequencies and domains

For every named common value, gnomAD v4.1.1 gives back the published population
structure without being told it: rs12913832 G at 0.76 in non-Finnish Europeans,
rs4988235 A at 0.64 there, rs1421085 C at 0.42, rs1042522 C (R72) at 0.75,
rs6983267 T at 0.60 in East Asians, rs1427407 as the minor allele at 0.19, and
APP's rs63750847 at 0.0003 - rare everywhere, as published for a variant found
in Iceland. That is 6 of 6 common values with the right frequency and the right
commonest population, derived. The seventh, ABO's O allele, is a single-base
deletion, and neither the frequency reader nor the local trace handles an indel:
recorded as pending.

The value-domain size separates the extremes and nothing finer. On the decade
scale, HLA-DRB1 reads 112.8 common positions per kb (many), ABO 17.0 (few),
HERC2 4.2 and APP 0.0 (two) - four of six as published - while TP53 reads 10.0
(few, expected two) and the lactase enhancer 5.0 (two, where the published
domain is a handful of alleles across continents, of which only rs4988235 is
inside this window). The 60 matched negative windows read a median of 2.5 per kb,
90% of them below 8.2 and the highest 25. So the reading
tells a hypervariable slot from an ordinary one, and cannot tell two values from
a few; separating those needs haplotypes, not positions, which is what the HPRC
panel (`attribution/human_panel.py`) is for.

One class result is worth more than the domain sizes. The "value on syntax"
reading - a position where our three people differ that sits on mammal-constrained
sequence - fires at 6 of 12 panel loci and at 9 of 60 matched windows, the only
claim in the benchmark that the controls do not erase. It puts rs1427407 first of
the twelve variable positions in the BCL11A enhancer at phyloP 7.01 and
rs12913832 first of four at 3.41. And it fails in an informative way at the
lactase enhancer: rs4988235 sits at phyloP 0.22, on sequence mammals never held
still, because it is a recent human adaptation. Its kilobase is 50%
human-constrained on the Gnocchi axis. A value slot can be old syntax with a
variable base or a recent change on free sequence, and only reading both axes
tells them apart.

## 6. The quota spent on the windows instead of the chromosomes (2026-09-14)

The chromosome sweep costs about 770 hours of model time genome-wide and buys
+0.07 AUC on direction; eight of the twelve loci were blocked only because their
chromosomes held no deletion data. genomeos-9c reallocated the quota, and
`scripts/loci_score.py` spent it on the panel's windows alone: the same scorer,
the same per-element cache, no whole-chromosome job, so the sweep finds these
answers when it resumes. Elements over the published element go first, then the
rest of the window, smallest window first.

| locus | chromosome | elements | requests | minutes |
|---|---|---|---|---|
| ABO | chr9 | 57 | 57 | 0.1 |
| HBB_LCR | chr11 | 51 | 51 | 0.1 |
| MCM6_LCT | chr2 | 106 | 106 | 0.2 |
| HLA_DRB1 | chr6 | 113 | 113 | 0.2 |
| HOXD | chr2 | 273 | 284 | 1.6 |
| BCL11A_enhancer | chr2 | 318 | 318 | 0.4 |
| SHH_ZRS | chr7 | 405 | 421 | 1.8 |
| MYC_8q24 | chr8 | 917 | 933 | 2.4 |
| HERC2_OCA2, FTO_IRX3, TP53, APP | chr15, chr16, chr17, chr21 | 739 | 0 | 0 | 
| **total** | | **2,979** | **2,283** | **8.4** |

2,283 requests, 43 quota waits, eight and a half minutes, at 16 requests in
flight. Four loci needed nothing because the sweep had finished their
chromosomes. The requests exceed the elements by the retries a quota answer
costs. For comparison, one chromosome of the sweep is 12,000 to 31,000 elements.

### What moved

| field | before (4 loci with deletion data) | after (12) |
|---|---|---|
| right target, derived | 8/12 | **11/12** |
| right target, nearest TSS in node | 8/12 | 8/12 |
| right target, annotation lookup | 6/12 | 6/12 |
| chance floor (random gene in the window) | 5.5/12 | 5.5/12 |
| direction, where a deletion names the target | 0/1 judged | **4/4 judged** |
| right cell or tissue | 4/10 judged | 4/10 judged |
| a published class read | 10/12 | 10/12 |

Every locus that moved moved the right way, and the derived rate now clears both
the heuristic and the chance floor. The direction of effect went from untestable
to right everywhere it could be judged: the lactase enhancer, the globin locus
control region, the BCL11A enhancer and HOXD are all called activating, as
published. Distances land where they should: LCT 13.9 kb against 13.9 published,
BCL11A 60.7 against 59, IRX3 519.7 kb against 520.

**The negative controls did not move, and that is still the headline.** With
element-level deletions now at 10 of the 60 matched windows as well, a derived
layer names *some* target at 52 of 60 (87%) against 11 of 12 loci (92%). What
separates the panel from its controls is which gene, and the direction: 8 of 12
loci produce one against 8 of 60 windows (13%).

### The two element-level failures, with the data in hand

- **The ZRS still names LMBR1.** Its own deletion, scored fresh, moves LMBR1 -
  the gene it sits inside - and calls it repressive in heart tissue. SHH, 979 kb
  away, is named only when every element of the window is summed, where it comes
  first. A 1 Mb reach is not something the element-level reading recovers, and
  this is the clearest statement of that the project has.
- **MYC's 8q24 enhancer names POU5F1B**, with MYC second in the window's summed
  input and the eQTLs naming CASC8 first. It is the one target miss left.

Two element-level misses of the ten loci where an element was scored; the count
is pinned in the CI gate, so a third is a regression.

At HERC2/OCA2 and FTO/IRX3 the new data changed a strict miss into a strict hit
the same way: the summed window input ranks OCA2 above HERC2, and IRX3 above
FTO, where the eQTL layer had ranked the neighbour first. The nearest-gene trap
at FTO is escaped by the model and not by the heuristic, which still answers FTO.

## 7. Lessons the runs produced

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
- **Cost.** The run is 21 to 48 minutes: local layers and two constraint axes over
  523 windows, one GTEx pass of 2 minutes, one bigBed query per window. The
  first version took six hours because it opened a bigWig session per negative
  window for the per-base reading; batching a chromosome's windows into one
  session is the whole difference. Public range reads time out now and then;
  one chr2 timeout silently emptied three loci's syntax reading, the gate test
  caught it, and the reads now retry three times.
- **Score the windows, not the chromosomes.** Twelve loci needed 2,283 requests and eight
  minutes; the same answers by sweeping their eight chromosomes would have been about
  160,000 elements. When a question names its loci, the shape of the job is the loci.
- **A named target needs a floor.** Against a chance expectation of 5.5 of 12,
  8 of 12 is not a result. Any rate quoted for a target-naming layer should come
  with the number of genes it could have chosen from.
- **A claim needs a coverage denominator too.** "A direction is produced" looked
  like the panel's second real discriminator at 8 of 12 against 8 of 60. It was
  a statement about where the requests had been spent; see section 8.
- **A panel that cannot fail is not a panel.** Widening from twelve to seventeen
  added the project's first two loci whose mechanism nothing here can represent
  (parent of origin, copy number) and a second megabase-reach case that behaves
  exactly like the first. A benchmark grows by adding what it is expected to
  fail at, not by adding what it will pass; see section 9.
- **Position-counting cannot see a haplotype domain.** CYP2D6 reads two common
  positions per kilobase where the published domain is a hundred star alleles,
  for the same reason the reader could not see ABO's frameshift. Both new storage
  loci underread, which takes the value-domain reading from 4 of 6 to 4 of 8.

## 8. The reruns as the chromosome sweep lands (2026-09-14, overnight)

The benchmark is rerun whenever the all-element sweep puts a scored deletion
inside one of the 523 windows. The run costs 6 to 11 minutes when GTEx's
distilled hits are already cached (the 21 to 48 minutes of section 5 included
the archive pass) and makes no AlphaGenome request: a locus that still needs
scoring is recorded as pending.

| run | commit | derived | heuristic | looked up | direction (panel) | control target | control direction |
|---|---|---|---|---|---|---|---|
| section 6 | ef80075 | 11/12 | 8/12 | 6/12 | 4/4 judged | 52/60 (87%) | 8/60 (13%) |
| rerun 1 | 02a1b98 | 11/12 | 8/12 | 6/12 | 4/4 judged | 52/60 (87%) | 8/60 (13%) |
| rerun 2 | b7092da | 11/12 | 8/12 | 6/12 | 4/4 judged | 52/60 (87%) | **9/60 (15%)** |
| rerun 3 | 88ae1d2 | 11/12 | 8/12 | 6/12 | 4/4 judged | 52/60 (87%) | **11/60 (18%)** |
| rerun 4 | this section | 11/12 | 8/12 | 6/12 | 4/4 judged | 52/60 (87%) | 11/60 (18%) |

Control windows holding an element-level deletion: 10, then 11, then 13, then 13
of 60. Every one of the three the sweep reached produced a direction on its first
reading. Nothing at the panel moved across any of the four runs.

Rerun 4 was taken when chr11 finished - the first panel chromosome to complete
since the quota was reallocated to the windows. Its last 16,000 elements put
eight more scored deletions inside the benchmark's 523 windows and changed no
verdict anywhere. The chromosome is saturated.

Rerun 1 reproduced every verdict of the committed run exactly: the same sixty
negative windows were redrawn from the same candidates, so the panel is
deterministic and any later movement is new data rather than run-to-run noise.
Its one change was the position check, which had timed out at Ensembl on the
previous pass and this time answered for six of the eleven named variants: all
six agree (rs8176719 within the anchor base its deletion is written on), none
disagree, five remain unverified.

Between the reruns the sweep reached one of the five matched negatives of the
beta-globin locus control region (chr11:100,686,471-100,712,471) and scored six
elements in it. That window immediately produced a direction, and the control
direction rate rose from 8 of 60 to 9 of 60.

### The direction claim is a coverage artefact, and this is the finding

Counting the same claims only where an element inside the window has actually
been deleted - the result now carries this as
`negative_controls.given_deletion_data`, so the number cannot be quoted without
its denominator:

| claim | panel, with deletion data | matched negatives, with deletion data | negatives, without |
|---|---|---|---|
| a derived layer names *some* target | 10/10 | 12/13 (92%) | 40/47 (85%) |
| the reader has it open in some cell | 9/10 (90%) | 9/13 (69%) | 23/47 (49%) |
| **a direction is produced** | **8/10 (80%)** | **11/13 (85%)** | **0/47** |
| **a value sits on constrained sequence** | **5/10 (50%)** | **2/13 (15%)** | **7/47 (15%)** |

The direction separation is gone. No window without a scored element can produce
a direction at all, and among windows that have one the matched negatives produce
a direction slightly *more* often than the panel. The panel looked better only
because `scripts/loci_score.py` deliberately spent 2,283 requests on the twelve
locus windows and nothing on the sixty controls. **Producing a direction is
evidence of having been scored, not of the element doing anything**, and it joins
"names a target" as a claim that must never be quoted on its own.

The `positives.direction >= 3 x negatives.direction` assertion has been removed
from the gate for that reason: it was pinning an artefact, and it would have
broken by arithmetic once about 14 of the 60 controls had deletion data. In its
place the gate requires the run to carry the coverage-conditioned counts, and
requires that no window without a scored element ever produces a direction.

**What survives conditioning is the one claim that survived the controls: a value
on syntax, 50% at the panel against 15% at the matched windows that hold the same
kind of data.** Separating it from the rest is what the benchmark is for. Note
that the value-on-syntax rate at the controls is the same whether or not they
carry a deletion (15% against 15%), which is what a claim independent of model
coverage looks like, and the opposite of what the direction row shows.

### What the sweep can still change

The twelve locus windows were scored end to end in section 6, so the sweep adds
little to them; nearly everything it can still move is on the control side,
where 47 of the 60 windows of that run had no element-level deletion. The
controls sit on the panel's own chromosomes, so the sweep's remaining order
matters to the benchmark only where it reaches chr11 (5 controls), chr9, chr8,
chr7 and chr6 (5 each) and chr2 (15). Each control window the sweep reaches has
about a 5-in-6 chance of producing a direction and about a 1-in-7 chance of
producing a value on syntax. (Section 9 widens the panel and restates these
counts over 85 controls.)

### How far the controls can actually go, measured rather than extrapolated

The first extrapolation here said the control direction rate would converge on
the panel's 80%. Finishing chr11 shows that it will not, and the reason is worth
more than the guess was. **Half of the matched windows contain no registry
element at all**, so no amount of sweeping can give them a deletion: two of the
five beta-globin controls hold nothing even on a finished chromosome.

Restricting to the controls whose chromosome the sweep has completed - chr11,
chr15, chr16, chr17, chr21, twenty-five of the sixty windows, where coverage is
no longer the variable:

| claim | panel | controls on finished chromosomes |
|---|---|---|
| an element inside the window has been deleted | 10/12 (83%) | 12/25 (48%) |
| a derived layer names *some* target | 11/12 (92%) | 20/25 (80%) |
| the reader has it open in some cell | 10/12 (83%) | 9/25 (36%) |
| a direction is produced | 8/12 (67%) | 10/25 (40%) |
| a value sits on constrained sequence | 6/12 (50%) | 2/25 (8%) |

So the honest projected end state, once the sweep reaches chr2, chr6, chr7, chr8
and chr9, is a control direction rate near 40% against the panel's 67% - not the
13% the first run reported, and not the 80% the last extrapolation predicted.
What is left of that gap is mostly the first row: a published functional element
coincides with an ENCODE cCRE more often than a window matched only on length,
GC, TSS distance and constraint does. That is a property of the registry, not of
the model, and it should be read as one.

The value-on-syntax claim is the exception again, and by a wider margin on this
subset than on the whole: 50% against 8%.

## 9. The panel widened to seventeen (2026-09-15)

Section 8 said the sweep had nothing left to give twelve loci, so the panel was widened instead:
a second long-range case (SOX9/Pierre Robin), an imprinted locus (H19/IGF2 ICR1), a dosage locus
(PMP22/CMT1A) and two more value domains (CYP2D6, MC1R). Seventeen loci, 85 matched negatives,
about 10 minutes and no AlphaGenome request. `scripts/loci_benchmark.py --gate` runs only the
twelve pinned loci into `loci_benchmark_gate.json` when the question is whether anything regressed.

| field | twelve loci | seventeen |
|---|---|---|
| right target, derived | 11/12 (92%) | **15/17 (88%)** |
| right target, nearest TSS in node | 8/12 (67%) | 11/17 (65%) |
| right target, annotation lookup | 6/12 (50%) | 8/17 (47%) |
| chance floor (random gene in the window) | 5.5/12 (46%) | 7.1/17 (42%) |
| direction, where a deletion names the target | 4/4 | 5/5 |
| a published class read | 10/12 | 13/17 |

Nothing that passed stopped passing. The derived layers still clear the heuristic and the floor,
by about the same margin, which is the first evidence that the twelve were not a lucky draw.

### What the five new loci say

- **SOX9/Pierre Robin behaves exactly like the ZRS, and that is the result.** The summed window
  input names SOX9 first, ahead of KCNJ2; the nearest-TSS heuristic answers KCNJ2, 500 kb away in a
  gene desert; and no element inside the 3 kb element window has been scored, so there is no
  element-level reading at all. Two megabase-reach loci now behave the same way: **summing a whole
  locus window reaches a megabase and reading one element does not.** One example was an anecdote.
- **H19/IGF2 ICR1 misses everything.** The summed window input names TNNT3, the heuristic names
  MRPL23, the lookups name nothing, and the element holds no common variable position at all
  (0.0 per kb), so no class is read either. It is the panel's second full miss beside MYC. That is
  the honest answer for a locus whose whole mechanism is parent of origin: the two alleles are the
  same sequence, and nothing read off a reference can tell them apart. Recorded as pending with
  what it would take - an allele-resolved methylation layer, or phased reads with a parental
  assignment - rather than as a defeat.
- **PMP22/CMT1A finds its target easily and misses the point.** The element's own deletion, the
  node and the lookups all name PMP22, and the direction is right. But the published variable is
  three copies of an intact gene against one, and no layer here counts copies: every reader asks
  what a sequence says, never how many times it is present. Pending, with the cost written down.
- **CYP2D6 confirms the prediction written into the panel before the run.** The value-domain
  reading gives **two** (5.4 common positions per kb) where the published domain is over a hundred
  star alleles. Counting positions per kilobase cannot see a system whose diversity is haplotypes,
  gene conversion with CYP2D7 and copy number - the same blindness that could not see ABO's
  frameshift. The causal base is found, though: rs1065852 ranks **1 of 57** variable positions by
  constraint at phyloP 8.26, as sharp a result as BCL11A's rs1427407.
- **MC1R also underreads its domain** (two, at 6.0 per kb, against a published few) but gives the
  panel its **third and fourth coding calibrations**: the translation engine derives p.Arg151Cys for
  rs1805007 and p.Arg160Trp for rs1805008 from the reference sequence and both match, beside APP's
  p.Ala673Thr and TP53's p.Pro72Arg. Neither position was typed into the panel; both came from
  Ensembl at run time and were then confirmed a second way by the engine's own translation.

On this run Ensembl answered for every hard-coded variant and **all eleven agree, none disagree** -
the first run where the position check is complete rather than partly timed out.

So the value-domain reading is now right at 4 of 8 storage loci, not 4 of 6. Both new storage loci
underread, in the same direction and for the same reason. **The reading tells a hypervariable slot
from an ordinary one and nothing finer, and it cannot see a domain built from haplotypes at all.**

### The controls moved, and they moved the way section 8 predicted

The five new loci sit on chr11, chr16, chr17 and chr22, all of which the sweep has finished, so
their 25 new control windows arrived with deletion data already in them:

| claim | panel 17 | all 85 controls | the 50 controls on finished chromosomes |
|---|---|---|---|
| a derived layer names *some* target | 14/17 (82%) | 77/85 (91%) | **45/50 (90%)** |
| the reader has it open in some cell | 14/17 (82%) | 54/85 (64%) | 31/50 (62%) |
| a direction is produced | 11/17 (65%) | 31/85 (36%) | **30/50 (60%)** |
| a value sits on constrained sequence | 9/17 (53%) | 15/85 (18%) | **8/50 (16%)** |
| an element inside was deleted | 13/17 (76%) | 34/85 (40%) | 33/50 (66%) |

The control direction rate went from 18% to 36% in one step, purely by adding loci on chromosomes
the sweep had finished, and on the fair subset it is **60% against the panel's 65%**. The
projection in section 8 said 40% from a sample of 25 windows where two thirds held no element;
with 50 windows, two thirds *do* hold one, and the observed 60% is what that revision predicts.
**Naming a target is now worse than worthless: the matched windows do it more often than the panel
does (90% against 82%).**

**The value-on-syntax row has since failed its first independent test (2026-09-17, section 19).** On
nine published loci the project had never seen, it reads **1 of 9 against 3 of 45 matched controls**,
and standardised it sits *below* its controls. It is not element length: on elements of a kilobase or
less the panel reads 4 of 7 and the candidates 1 of 9. Nine loci cannot refute 9 of 17, so the row
above is not struck — but the one claim that ever separated this panel from its controls can no
longer be quoted as established, and the honest form of it is that it failed the first test it was
given. Section 19 also reproduces section 8's retraction on these
independent loci and sharpens it: the apparent separation **survives standardising on GC, distance
and constraint** — none of which says whether anybody ever spent a request on the window — and dies
only when coverage enters as a stratum, taking the direction claim from +0.41 (p 0.026) to −0.11,
with the controls producing a direction more often than the loci. Section 8 was right, and
**standardising on covariates is not a substitute for conditioning on coverage**, which is the part
that generalises past this benchmark.

**A value on constrained sequence is the only claim left standing: 53% at the panel against 16% at
matched windows with the same data behind them.** Every other claim the benchmark scores has now
been shown, on its own controls, to carry no information on its own.

## 10. The element-level reading is bounded by the annotation, not by the sweep (2026-09-15)

> **CORRECTED 2026-09-17, and one sentence below is struck rather than edited.** The counts in this
> section came from `Chromosome.ccres_in`, which bisected an index built straight from `load_ccres`.
> The registry does not arrive in start order: it is several sorted runs concatenated, one per class
> group, with three descending steps on every chromosome checked. A bisect over that returns a
> plausible index and silently drops everything before it. Measured against a brute-force scan, four
> of the panel's seventeen loci were wrong — **MC1R read 1 cCRE where it has 12, and the ZRS, MYC and
> ABO read 0 where they have 1, 2 and 4** — as were 14 of the 85 negative windows, one of them 21
> where the truth is 108.
>
> **"No cCRE over the published element" is 4 of 17, not 7.** SHH_ZRS, MYC_8q24 and ABO leave that
> list; APP, TP53, SOX9_PierreRobin and H19_ICR1 remain, and they are the loci the section's argument
> actually rests on. **The sentence "The flagship long-range element is not an ENCODE cCRE" is
> false** — there is one over the ZRS and the lookup was failing to see it. It is left in place below,
> struck, because a claim this section leaned on should be readable next to its correction.
>
> **What survives, and it is the thesis.** SOX9's element and the H19 ICR still hold zero elements in
> every universe the scorer can ask about, so the section's central point stands: no quota can produce
> an element-level verdict where no annotation drew an interval. The panel-wide count was wrong; the
> reason the two loci could not be scored was not. §18 later found a second and independent reason
> SOX9 cannot be scored — its target is outside the model's 1 Mb input — so that locus is now doubly
> unanswerable.
>
> Found by genomeos-79's more-loci lane and verified here against a brute-force scan before anything
> was changed. Fix and regression test in `benchmark/loci.py` and `tests/test_loci_benchmark.py`;
> `loci_benchmark.json` re-run, and the headline `target_derived` is unmoved at 15/17.


Asked to spend a handful of AlphaGenome requests on the two loci with no element-level reading
(SOX9's element and the H19 ICR), the first thing to check was what there was to spend them on.
The answer is nothing, and it is worth more than the requests would have been.

**Both windows hold zero elements in every universe the scorer can ask about**: 0 ENCODE cCREs,
0 VISTA, 0 lentiMPRA. The sweep did not skip them. Everything nearby is already scored - 11 cCREs
within 20 kb of the SOX9 element, 2 within 20 kb of the H19 ICR - and none of those moves a coding
gene. No quota can produce an element-level verdict where no annotation drew an interval.

Measured across the panel (`aggregate.published_element_coverage` in the result):

| | loci |
|---|---|
| no cCRE over the published element | ~~7 of 17: SHH_ZRS, MYC_8q24, ABO, APP, TP53, SOX9_PierreRobin, H19_ICR1~~ → **4 of 17**: APP, TP53, SOX9_PierreRobin, H19_ICR1 (corrected 2026-09-17) |
| no element-level reading of any kind | **4 of 17**: APP, TP53 (coding loci, scored by translation instead), SOX9_PierreRobin, H19_ICR1 |

And the ZRS's element-level answer - the LMBR1 miss pinned in the CI gate, the single sharpest
result the project has - **does not come from the registry at all**. It comes from the VISTA run.
~~The flagship long-range element is not an ENCODE cCRE.~~ **Struck 2026-09-17: it is one. See the
correction at the head of this section.**

So: **the element-level reading is bounded by which intervals ENCODE and VISTA happened to call,
not by how much of the genome the model has covered.** Sweeping harder cannot fix it. Every rate
this benchmark reports for an element-level layer is conditioned on somebody else having drawn the
element first, which is a looked-up step sitting underneath a derived claim.

## 11. Design: an interval the project defines itself (2026-09-15, before any locus is added)

Two gaps turn out to be the same missing idea. `EPHA4` and `IHH` were dropped from the widening
because the panel cannot write a structural variant down as an element. SOX9 and H19 have no
element-level reading because no annotation called their published element. Both need the same
thing: **a causal unit the project states, rather than one it looks up.** This section is the
design, written before a locus is added, because that is the discipline the rest of the panel was
built with.

### 11a. The easy half: a published interval nobody annotated

`Expect.element` is already a stated interval; what is missing is permission to *delete* it. The
change is to let the scorer take an interval as well as a registry id, so the flagship question -
does deleting the published element name the published gene - can be asked at SOX9 and H19 the way
it is asked at the ZRS. Cost: a handful of requests per locus, and a small change to the scorer's
input, which currently takes registry elements only. **This is a decision for the coordinator, not
something to do while holding a key**, because it changes what the scorer is asked, not just how
much of it runs.

### 11b. The real half: the variable is a rearrangement

At EPHA4 the published answer is not that an element was lost. It is that a deletion, inversion or
duplication **moved a boundary**, so limb enhancers that belonged to EPHA4 now reach IHH, PAX3 or
WNT6, and the limb gets the wrong instruction. The causal unit is a pair of breakpoints and a type,
and the published effect is a **new adjacency**. Nothing in `Expect` can say that.

**What the expectation record holds.** A `rearrangement` field beside `element`:

| field | what it says |
|---|---|
| `kind` | deletion, inversion, duplication, translocation |
| `breakpoints` | two intervals, carrying their published uncertainty rather than pretending to a base |
| `boundary` | the interval of the boundary the rearrangement removes, as published |
| `gains` | the adjacencies created: (donor element or region, recipient gene) |
| `loses` | the adjacency destroyed |
| `phenotype` | what the person or the mouse shows, in words |
| `citations`, `answer_from` | as everywhere else |

**What counts as derived, and it is not what it first looks like.** The tempting reading - our node
model, recomputed on the rearranged coordinates, puts the donor element and the recipient gene in
one node when it did not before - is **inferred, not derived**: the node model is built from
CTCF-only cCREs at confidence 0.4, and it is the very thing under test. So it takes the heuristic
slot that nearest-coding-TSS holds today: the baseline the rest has to beat. That is the honest
place for it, and it makes a rearrangement locus the sharpest test of the node model the project
could have.

| provenance | the reading | cost |
|---|---|---|
| derived | AlphaGenome asked on the rearranged sequence: does the donor element's predicted target change when the boundary is cut | **needs work in the scorer**, which today deletes one element and cannot construct a rearranged input. Requests are the small part |
| heuristic | the node model recomputed on rearranged coordinates: does the boundary disappear and the new adjacency appear | free and local |
| looked_up | DECIPHER and ClinVar on the rearrangement, VISTA on the donor element | nothing about the approach |
| free, and worth asking first | is there a CTCF-only element at the published boundary **at all**? At HOXD the answer was no, and that was one of the panel's sharpest negative results | free |

**What the matched control is.** This is the part that decides whether the locus is worth adding,
because everything the benchmark learned tonight says a claim is worth nothing until a matched
random version of it has been tried. A window is the wrong control for a rearrangement; the unit is
a pair brought together by a cut. The control is therefore **a random rearrangement of the same
kind and size elsewhere on the same chromosome**, five per locus, matched on:

- kind (exact), and span (exact, as window length is matched today);
- coding TSSs inside the span, within 35% - the analogue of the TSS-distance match;
- CTCF-only elements crossed, closest of the candidates - the analogue of the constraint match,
  since boundary density is what decides whether a cut can create an adjacency at all;
- away from every panel locus with 200 kb of flank, and away from known pathogenic CNVs.

Three claims are then counted identically at the published rearrangements and at the random ones:
the node model **loses a boundary**; a **new element-to-gene adjacency** appears; some derived
layer **names the recipient gene**. Only the fourth - whether the gene named is the published one -
belongs to the positives alone.

**The prediction, written down before the run, as CYP2D6's was.** Cutting the genome at random will
destroy boundaries and create new adjacencies constantly, so claims one and two will fire at most
control rearrangements, and claim three will land near the 90% that naming a target already reaches
at matched windows. **What should separate the published rearrangements is which gene, and nothing
else.** If it does not, the node model has failed its sharpest test, and that is the result the
project most needs to know.

### 11c. Order of work

The measurement in section 10 has to come first and is already done: it says how much of the panel
the element-level reading can reach at all. Then 11a, which is cheap and needs a decision on the
scorer's input. Then 11b, whose free half - the node model recomputed on rearranged coordinates,
against matched random rearrangements - can be built and run with no model request whatsoever, and
answers the question on its own. The derived half needs scorer work and should wait until the free
half says whether there is anything there.

## 12. 11a built, and 11b's free half run (2026-09-15)

Both approved by the coordinator. Neither spent a model request; the key went straight to the
executor lane.

### 11a: the scorer already takes a stated interval

**Cost: nothing.** `Context.score_region` in `genomeos/predict/enhancer_target.py` already scores a
region as an ad-hoc element when no ENCODE element overlaps it, caching it under
`{chrom}_{start}_{end}`. It returns `predicted` and `predicted_coding` - exactly what the
benchmark's deletion layer reads. The only thing a stated interval does not get is `verdict`, the
comparison against the domain's nearest-TSS inference, which the benchmark computes itself anyway.
No rewrite was needed and none was attempted.

The three guardrails are in the code, not in this paragraph:

- **a stated interval carries its own citation.** `Expect.element_source` is `annotated` or
  `stated`, and `Expect.element_citation` is required for a stated one - the gate fails a stated
  interval with no citation, because that is the panel inventing an element;
- **every result from one is labelled**, in the row (`readings.deletion.from_stated_interval`) and
  in the aggregate (`target_derived_by_element_source`), and a stated-interval run is read from its
  own result file so its rows can never be mistaken for registry ones;
- **the rates are never pooled without the split showing.** `target_derived_by_element_source`
  reports annotated and stated apart, with the loci in each, and the gate asserts they add up to the
  panel and that the stated set is exactly the loci declaring it.

SOX9_PierreRobin and H19_ICR1 are the two stated intervals. Scoring them is a handful of requests
whenever the key is next free; until then their deletion reading stays pending, as it was.

### 11b free half: the node model against a published rearrangement

`genomeos/benchmark/rearrangements.py`, `scripts/rearrangements.py`,
`data/results/loci_rearrangements.json`. **No model request, no network**: the cCRE registry,
GENCODE and the project's own `infer_domains`, recomputed over a transformed cCRE list and a
transformed chromosome length - the same function, not a simplification of it.

One published case is runnable: the **EPHA4 to PAX3 deletion** (Lupiáñez et al. 2015), where
deleting across the boundary gives PAX3 the EPHA4 limb enhancers and causes polydactyly. Its span
is **stated, not cited**: the paper's breakpoints are hg19 and we have not lifted them over, so the
span is the intergenic interval between GENCODE's EPHA4 and PAX3, which is the interval the
published deletions remove the boundary from. Three more cases are in `NOT_YET_RUNNABLE` with what
each waits for - the IHH duplication and the WNT6 inversion both depend on breakpoints falling
*inside* the domains, so an intergenic stand-in is not faithful the way it is for a deletion.

**The free reading, asked first: there are 7 CTCF-only elements between EPHA4 and PAX3.** Unlike
HOXD, where no CTCF-only element sat anywhere in the cluster and the boundary was a fiction, this
boundary exists in our registry. The node model does separate the pair (`chr2:D1512` and
`chr2:D1517`), which is the precondition the experiment needs.

| claim | published | 5 matched random rearrangements |
|---|---|---|
| the pair is in different nodes to begin with | 1/1 | 5/5 |
| the rearrangement loses a boundary between them | 1/1 | **5/5** |
| a new adjacency appears - they end up in one node | 1/1 | **5/5** |
| a derived layer names the recipient gene | 1/1 | **not judgeable** |

**The registered prediction was right, and the node model failed its sharpest test.** Deleting the
interval between any matched pair of consecutive genes destroys every boundary between them and puts
them in one node - 5 times out of 5. The published rearrangement does exactly what a random one of
the same shape does. Nothing in the node model distinguishes the deletion that causes polydactyly
from a deletion of the same length crossing the same number of boundaries somewhere else on chr2.

The fourth claim - which gene - is the only one that could have separated them, and it **cannot be
judged**: all five control recipients have no deletion scored anywhere near them, because chr2 is
not swept. An ungated 1/1 against 0/5 would have been the coverage artefact of section 8 all over
again, so the result records it as not judgeable and the gate refuses to quote it.

### Two mistakes this run made, both caught by its own controls

Worth recording because both looked like findings first.

1. **The control has to share the construction, not only the size.** The first version drew a random
   span of the published length and took the nearest coding TSS on each side. The published span is
   the whole interval *between* two genes, so deleting it removes every boundary between them by
   construction, while a random span leaves the boundaries between itself and the flanking genes.
   That produced a new adjacency at 1/1 against 0/5 - a beautiful result, and entirely an artefact
   of how the two were built. A control is now the same thing done elsewhere: a consecutive coding
   pair and the whole interval between them.
2. **An endpoint must not sit inside its own deletion.** With the span running from one TSS to the
   next, the lower TSS was inside the interval being removed, so it transformed to nowhere and every
   downstream count was nonsense (boundaries "after" in the hundreds). The span now lies strictly
   between the two starts.

Both were visible only because the claims are scored at controls as well as at the case. A run that
reported the positives alone would have published the first one.

### What this says about the node, beside tonight's other three qualifications

The node claim was already qualified three ways: 2.6 points above random on the full archive, an
excess that changes sign across chromosomes, and a sign that tracks boundary density rather than
biology. This is the fourth and the most direct: **at the one published experiment designed to test
exactly this, the node model's answer is indistinguishable from the answer it gives to a random cut
of the same shape.** One case is an anecdote by the panel's own standard, and the honest next step
is the other three cases, which need the paper's breakpoints lifted to hg38 - a liftover, not a
quota.

## 13. The liftover that was not the problem (2026-09-15)

Asked to lift the other three Lupiáñez breakpoints from hg19 to hg38 and run them, the first step
was to find the breakpoints. **There are none to find.** The accessible paper gives sizes and gene
content and no human coordinates at all: deletions of 1.75 to 1.9 Mb at 2q35-36, an inversion of
about 1.1 Mb whose telomeric breakpoint is about 1.4 Mb from EPHA4, a duplication of about 1.4 Mb
with a breakpoint about 1.2 Mb away, and a polydactyly duplication of about 900 kb. The breakpoints
themselves live in aCGH supplementary data. Two independent routes were tried - the PMC full text
and POSTRE's re-curation of the same variants - and neither carries a coordinate.

**So the blocker was never the liftover, and saying it was sent the next session after the wrong
thing.** UCSC's chain and this project's own chain reader (`attribution/human_panel.py`, already
used for Repli-seq) would do the conversion in seconds. There is no number to convert. The gate now
refuses any `NOT_YET_RUNNABLE` reason containing the word "liftover", so the mistake cannot be
written back in.

### The rule applied to my own case

The instruction was to drop a case rather than place it approximately. Applied honestly, that takes
the case section 12 *ran*: **EPHA4 to PAX3 is now dropped too.** Looking for the breakpoints turned
up two errors in the expectation record I had already committed:

- **the phenotype was wrong.** The EPHA4-to-PAX3 deletion causes **brachydactyly**, not polydactyly;
  polydactyly comes from a duplication of about 900 kb. The record said polydactyly.
- **the stated span was wrong in kind.** The published deletions are 1.75 to 1.9 Mb and **include
  EPHA4 itself**, extending into the non-coding part of the PAX3 domain. The span I stated was the
  626 kb intergenic interval between the two genes, which excludes EPHA4 - three times too short
  and a different construction. It was not the faithful stand-in I claimed it was.

Neither error changed the verdict, for a reason worth more than the case was.

### Two of the three claims could never have discriminated, and that is the finding

The node model puts its boundaries at CTCF-only elements. A deletion removes every element inside
its span. So **any** deletion spanning a boundary loses that boundary, and any deletion of the whole
interval between two genes puts them in one node. That is arithmetic, not biology. Deleting an
interval of the published size (1.75 to 1.9 Mb) between a consecutive coding pair anywhere on chr2
loses a boundary and creates a new adjacency at **4 of 4** attempts, exactly as it did for the
approximately-placed case.

| claim | a deletion of the published size, anywhere on chr2 |
|---|---|
| the pair is in different nodes to begin with | 4/4 |
| loses a boundary between them | **4/4** |
| a new adjacency appears | **4/4** |
| a derived layer names the recipient | 0/4, and **0 of 4 recipients have any deletion scored near them** |

**Two of this module's three claims were incapable of discriminating by construction.** Whatever the
real breakpoints turn out to be, they cannot separate a published rearrangement from a random one on
either claim, because both are guaranteed by how the node model is built. Only the third - which
gene the machinery names - could ever have worked, which is what the registered prediction said, and
it still cannot be scored because no control recipient has a deletion near it until chr2 is swept.

That sharpens what the coordinator is spending chr2 on. It is not that claim three is the last of
three; it is that claim three is **the only one there ever was**.

## 14. chr2 landed, claim 3 scored, and the verdict (2026-09-15)

chr2 completed at 79,639 elements, the largest chromosome in the sweep, promoted deliberately so
that the only claim capable of discriminating could be scored at the controls as well as at the
case. It now can be.

**The case was reinstated, on a stated span, and that word is load-bearing.** The published
breakpoints still do not exist in any accessible form; what is citable is the size. A span anchored
at PAX3's GENCODE gene start and sized across the published 1.75 to 1.9 Mb range removes EPHA4
entire at every size and leaves PAX3 alive to be misexpressed, which is the published construction -
unlike the 626 kb intergenic span this file used before, which excluded EPHA4 and was withdrawn.
The case is run at every size across the range rather than at one, so the published uncertainty is
reported instead of collapsed. Phenotype corrected to brachydactyly. The other three cases stay
withdrawn.

### What a faithful span turns out not to be able to ask

At every size in the published range the deletion **removes the donor gene itself**. So
"donor and recipient end up in one node" is not a question this case can put: EPHA4 does not survive
to share a node with anything. The published mechanism is that EPHA4's surviving *enhancers* reach
PAX3, and this module has no published coordinate for that enhancer cluster. Claims one and two are
therefore recorded as **not applicable** for the case - a question the construction cannot ask is
not a test the model failed - and they remain, as section 13 showed, arithmetic anyway: 5 of 5 at
matched random deletions, 4 of 4 at size-matched ones anywhere on chr2.

### Claim 3, the only claim there ever was

| claim | published | 5 matched random rearrangements |
|---|---|---|
| the pair is in different nodes to begin with | 1/1 | 5/5 |
| loses a boundary between them | **not applicable** | 5/5 |
| a new adjacency appears | **not applicable** | 5/5 |
| **a derived layer names the recipient gene** | **1/1** | **5/5** |

Gated on coverage, as it must be: chr2 is swept, so 4 of 4 control recipients now have deletions
scored near them and the claim is judgeable for the first time. It separates nothing. PAX3 is named
by an already-computed deletion in its neighbourhood, and so is the recipient of every matched
random deletion, 5 times out of 5.

**So: the node model gives a random cut's answer on every claim the experiment can ask.** That is
the verdict. Two claims were arithmetic and could never have discriminated; the third, the one this
whole module was built around and the one chr2 was promoted for, gives the same answer at the
published rearrangement and at matched random cuts of the same size on the same chromosome.

It is the same shape as the panel's oldest finding, one level up: naming a target fires at 90% of
matched windows, and naming the gene at the far end of a boundary-crossing deletion fires at 100% of
matched deletions. What a gene is near, and what a cut puts it next to, are not evidence.

### What this does and does not say

It does not say TAD boundaries are not real. It says **our** boundary model - CTCF-only cCREs, no
Hi-C, confidence 0.4 - carries no information that distinguishes a pathogenic rearrangement from a
random one, on any claim this experiment can put to it. That is the mechanistic counterpart to the
statistical result the fold lane reached the same night: an excess distinguishable from random on
six chromosomes, all cut between 6.55 and 7.00 boundaries per Mb, with the finest-cut chromosome of
all indistinguishable, and ten of sixteen carrying no measurable excess either way.

The honest next step is not another rearrangement. It is a boundary model that is not a proxy:
Hi-C-called domains, which the project already holds for five cell types under
`data/knowledge/hic`, against which the same three claims could be asked without changing a line of
this module.

## 15. The Hi-C test: a positive that did not survive its own scrutiny (2026-09-15)

The negative in section 14 was about **our** boundary model - CTCF-only cCREs, no Hi-C, confidence
0.4 - and that left the interesting question open: does the negative belong to the proxy or to
boundaries themselves? The project holds measured 4D Nucleome boundary calls for five cell types, so
the same three claims can be asked of them. `_domains_from` is what `infer_domains` calls once it has
decided where the boundaries are, so handing it measured calls swaps the boundary source and changes
nothing else: same 50 kb minimum, same merging, same intervals-between-boundaries model. Same matched
random rearrangements, so only one thing changes at a time.

**The reading was registered before the measured boundaries were read**, and it is kept verbatim in
the result. It said: claim three cannot move, because it does not depend on the boundary model at
all; claims one and two are askable at the controls but not at the published case, whose faithful
span deletes the donor; what is left is the precondition, and the yes/no form of it will saturate,
so **the count of boundaries between the pair is the reading to trust**, per cell type, never pooled.

### Coverage first, as it must be

All five cell types have a measured domain over both EPHA4 and PAX3, and all five separate them -
as do all five matched control pairs, in all five cell types. The yes/no form saturates exactly as
registered, so it decides nothing. This is not a coverage failure: the locus is covered everywhere.

### The count, raw, looked like the first positive result about a node in this project

| cell type | boundaries between EPHA4 and PAX3 | at the 5 matched pairs | below every control? |
|---|---|---|---|
| GM12878 | 11 | 24, 26, 32, 37, 37 | **yes** |
| H1-hESC | 3 | 5, 7, 8, 9, 12 | **yes** |
| HepG2 | 2 | 3, 3, 4, 4, 5 | **yes** |
| IMR-90 | 2 | 6, 7, 7, 7, 11 | **yes** |
| K562 | 2 | 3, 3, 4, 4, 4 | **yes** |

Five cell types out of five, the published pair strictly below every matched control. If that had
held it would have been the first positive thing anyone found about a node here.

### It is interval length, and nothing else

**The published pair is 727 kb apart. The control pairs are about 1.8 Mb apart** - because the
controls were matched on the *span of the deletion*, not on the distance between the two genes. A
longer interval holds more boundaries for nothing. Per megabase:

| cell type | published, per Mb | at the 5 matched pairs, per Mb | below every control? |
|---|---|---|---|
| GM12878 | 15.14 | 13.60, 15.73, 16.70, 17.69, 20.43 | no - second lowest of six |
| H1-hESC | 4.13 | 3.02, 3.97, 4.42, 4.70, 5.74 | no - middle |
| HepG2 | 2.75 | 1.43, 1.57, 2.27, 2.42, 2.76 | no - second **highest** |
| IMR-90 | 2.75 | 3.40, 3.65, 3.86, 4.23, 5.26 | yes |
| K562 | 2.75 | 1.57, 1.81, 1.91, 2.21, 2.27 | no - **above every control** |

Below every control in **1 of 5** cell types, above every control in one, inside the range in three.
The 5-of-5 became 1-of-5 on dividing by a length. Both readings are now carried together in the
result and the gate requires it, so the raw count can never be quoted alone again.

### So the negative is about boundaries, not only about our proxy

That is the larger statement and it is worth making carefully. Measured Hi-C boundaries, in five
cell types, separate the published rearrangement's gene pair no more distinctively than they
separate matched random pairs on the same chromosome. The node model's failure in section 14 is not
an artefact of using CTCF-only cCREs as a stand-in: swapping in the measured calls the project holds
does not rescue it.

What this does **not** say is that boundary strength carries nothing. Insulation is a continuous
quantity and this test counted calls; a boundary's insulation score, which 4DN publishes and this
project does not yet read, is the obvious next thing to ask, and it is the reading that could still
come out positive. What it does say is that **the count of boundaries between two genes, measured or
inferred, is not evidence about whether a rearrangement between them matters.**

### The lesson, which is the same one as last time in a new costume

The control shared the deletion's span and not the pair's separation. Section 12 recorded that a
control has to share the construction and not only the size; this is the same mistake one level
further in, and it produced a five-out-of-five positive that survived until it was divided by a
length. **Every count needs its denominator, and the denominator has to be the thing that varies.**

## 16. Boundary strength: the last reading, and it is a null (2026-09-15)

Counting calls is not measuring insulation, so the last thing left to ask was the quantity the calls
are a threshold on. It turns out the project already had it and was throwing it away:
**4DN's boundary files carry five columns - chrom, start, end, a Strong/Weak label and a numeric
strength - and `genome/hic.py` writes only the first three.** The informative column was being
discarded at the point of download. This module re-streams the same accessions named in each cell
type's manifest, keeping the score, so the comparison is against exactly the files the call counts
came from.

### What was registered, before anything was fetched

The denominator problem was already known, so the statistic was chosen to avoid it rather than to be
corrected afterwards. The interval *between the genes* is 727 kb for the published pair and about
1.8 Mb for the controls, and any "deepest point in the interval" statistic is biased towards the
longer one. **So the primary reading uses the deleted span instead**, which the controls were matched
on for length within 35%: length-matched by construction, needing no normalising. It is also the
biologically right question - a deletion removes whatever boundary lies in its span, and the
published claim is that the boundary it removes is the one that mattered. The mean strength inside
the span is reported beside it as a length-independent second reading, and the span lengths are
reported so the match can be checked rather than assumed.

One honest adjustment: the registration was written expecting an insulation profile, where a
boundary is a *minimum*. What 4DN's file actually carries is a per-boundary strength where higher is
stronger. Same question, opposite sign, and it is recorded here rather than quietly rewritten.

The registered prediction was: **no separation.**

### The result

| cell type | boundaries in the span (strong) | strongest, published | strongest, 5 matched control spans | stronger than all? |
|---|---|---|---|---|
| GM12878 | 29 (11) | 1.90 | 1.70, 1.75, 1.99, 2.05, 2.25 | no |
| H1-hESC | 8 (4) | 1.83 | 0.88, 1.31, 1.36, 1.39, 1.73 | **yes** |
| HepG2 | 3 (2) | 0.68 | 0.73, 1.10, 1.28, 1.37, 2.05 | no - weakest of six |
| IMR-90 | 6 (2) | 1.39 | 0.48, 0.74, 0.76, 1.17, 1.46 | no |
| K562 | 3 (1) | 0.59 | 0.57, 0.83, 1.03, 1.08, 1.62 | no - second weakest |

Stronger than every matched control in **1 of 5** cell types. In two of the five the published span's
strongest boundary is the weakest or second weakest of the six. On the mean, published is above every
control in 2 of 5 (H1-hESC and IMR-90). **There is no consistent separation**, and the prediction
registered before the fetch was right.

### The tempting exception, named so it cannot be quietly promoted

The one cell type where the published span's boundary is stronger than every control is **H1-hESC,
the embryonic stem line - and limb malformation is an embryonic phenotype.** That is a good story and
it is not a result. Picking the one cell type of five that fits, after looking, is precisely the
artefact this benchmark has caught three times in a day. It is recorded as a hypothesis with an
obvious test - the other three Lupiáñez cases, and other published limb rearrangements, scored the
same way - and it stays a hypothesis until something like that is run. One locus, one cell type, no
correction for having looked at five.

### The honest summary

**Measured chromatin boundaries, whether read as calls or as strengths, do not distinguish this
published rearrangement from a random cut of the same kind.** Together with sections 14 and 15 that
is the whole of what this experiment could ask: the CTCF-only proxy fails, the measured calls fail,
and the strength behind the calls fails. The node concept as this project uses it carries no
information at this locus.

**And it is one locus**, on a stated span, with the other three published cases still withdrawn for
want of coordinates. That is the load-bearing caveat and it belongs in every quotation of this
result. What would change it is not another reading of the same locus but more loci - which is the
same conclusion the panel reached at twelve, and the reason it is now seventeen.

### What to fix regardless of the verdict

`genome/hic.py` discards columns 4 and 5 of every 4DN boundary file it downloads. Whatever the
answer here, a layer that reads a threshold crossing and throws away the quantity it is a threshold
on is worth repairing, and the fix is three lines in `fetch_boundaries` plus a loader beside
`load_boundaries`. That file is not mine; the finding is passed to its owner with the cache this
module built (`data/knowledge/hic_scored`) as evidence that the column is there and useful.

## 17. The sweep finished, and the coverage artefact of section 8 is confirmed at full coverage (2026-09-16)

The all-elements sweep completed on every chromosome, so this benchmark was rerun with no model
request: the loci it had been queued behind a quota handover for (SOX9, H19) were already scored, and
the handover is not needed.

**This is not a new finding.** Section 8 called the direction claim a coverage artefact on 2026-09-15,
on 13 of 60 control windows that had deletion data, and predicted that the separation would vanish as
the rest were scored. It has, at full coverage, which is the prediction being met rather than a
discovery:

| | 2026-09-15 | now |
|---|---|---|
| control windows with a scored deletion | 34 of 85 | **60 of 85** |
| a direction produced, published loci | 11/17, 0.647 | 11/17, 0.647 |
| a direction produced, matched controls | 31/85, 0.365 | **54/85, 0.635** |

Unconditioned, the direction rate at the controls has risen from 0.365 to 0.635 while the panel's
stayed at 0.647: **0.647 against 0.635 is nothing**, exactly as section 8 said it would be, and by the
arithmetic it named (the assertion it removed would have broken once about 14 of the controls had been
scored; 60 have).

**Everything else is unchanged.** 15 of 17 targets derived, 13 of 14 where the target is reachable, a
cell in 14 of 17, a value on syntax 9 of 17 against 15 of 85 at the controls — the one claim that
separates the panel, and the one whose control rate has not moved with coverage at any point.

**What the rerun is worth, stated plainly.** It cost no requests and confirmed a prediction rather
than producing a result. The reason to record it is that the prediction was made when its own evidence
was 13 windows, and it is the third time this benchmark has had to withdraw something for the same
reason — so a rate whose denominator is "was it measured" is now the first thing to check here, not the
last.

---

## 18. The question the model was never asked: reach, and 11a spent one request (2026-09-17)

11a has been waiting for a free key since section 12. The key came free, and the first thing the
driver did was refuse to spend most of it.

### Two of the panel's targets are outside the model's input, and one of them has been graded as a miss for months

The deletion scorer resizes its input to **1 Mb centred on the element** (`SEQUENCE_LENGTH_1MB`), so
its reach is 524 kb each way. Against GENCODE gene bodies, three published targets fall outside it:

| locus | target | TSS distance from the element | gene body in the model's input |
|---|---|---|---|
| SHH_ZRS | SHH | 979 kb | **no** |
| SOX9_PierreRobin | SOX9 | 1,450 kb | **no** |
| FTO_IRX3 | IRX5 | 1,164 kb | **no** (IRX3 at 520 kb is in, by 4 kb) |

**The test is gene-body overlap, not TSS distance, and the difference is not academic.** The
tempting shorthand — "the reach is 524 kb, so a target further than that cannot be named" — is
wrong, and the finished sweep says how wrong. On chr21, of 5,174 elements with a coding target,
**79 name a gene whose TSS is 524 kb or more away, up to 986 kb**, and **all 79** are genes whose
body overlaps the window. They are almost all **RUNX1**, which is 1.22 Mb long: its promoter is far
outside the window of elements that nonetheless name it, because its body runs right through them.
`read_reach` uses bodies for this reason. The three rows above fail the body test too, which is why
they are listed; the distances are given as TSS distances only because that is what the panel's
expectations record.

**This retires a result the benchmark has been carrying since section 1.** The ZRS's own deletion
names LMBR1, the gene the ZRS sits inside, and that was read as the flagship long-range miss — the
element-level question failing where the summed-window one succeeded. It was not a miss. SHH was
never a candidate: it is not in the model's input, and no number of requests puts it there. The
panel chose a target the model cannot see, and then marked the model wrong for not seeing it.

`read_reach` now computes this per locus, free, from gene bodies and arithmetic. A locus whose
targets are all out of reach has its deletion layer marked **unaskable** rather than **pending**, and
the pending line — which used to quote a price in requests — now says no number of requests will do.

**The headline rates are deliberately unchanged.** `target_derived` stays 15/17. A new line beside
it, `target_derived_where_the_model_could_answer`, holds the two unaskable loci out: **13/15,
0.867** against the headline's 0.882. Note the direction — both held-out loci were counted as
*hits* through other derived layers, so holding them out makes the number slightly **worse**. A
correction that flattered the model would deserve more suspicion than this one does. Which
denominator is the real one is the benchmark owner's call, and the split is now visible so it can be
made.

### 11a: one request, at the only locus where the question survives

SOX9 was the other stated interval, and it is unaskable by the table above, so the run asked
nothing there. That leaves **H19_ICR1: one request.** The prediction was written into
`scripts/loci_stated_intervals.py` before it, with three outcomes — `hit` (IGF2 or H19 first),
`trap` (MRPL23, the nearest coding TSS), `neither` (a third gene) — and the registered expectation
was **not a hit**, because the ICR acts through parent-of-origin methylation, which no layer in this
project carries and which a sequence model reading one allele cannot express.

**The registered outcome was `neither`, and the reason was the registered one.** Deleting the ICR
named INS-IGF2 first among coding genes. But what it named is the whole point:

| rank | gene | effect | what it is |
|---|---|---|---|
| 1 | MIR675 | −0.243 | the miRNA hosted **inside H19** |
| 2 | ENSG00000274866 | −0.241 | |
| 3 | **H19** | −0.178 | a published target |
| — | IGF2-AS | −0.403 (strongest overall) | the antisense of the other published target |
| — | INS-IGF2 | −0.149 | the read-through covering IGF2; first among **coding** genes |

Every gene it moved is in the IGF2/H19 imprinted cluster, and **MRPL23, the nearest-TSS trap, was
not named at all.** So the element-level deletion put its effect in exactly the right place and
still scores as a miss, because the published targets are named as *IGF2* and *H19* and what came
back was their read-through and their antisense.

**A limitation this locus exposed, reported rather than fixed.** The deletion layer's named list
takes `predicted_coding` first and stops, so a published target that is **non-coding** can be hidden
behind a coding read-through. H19 is a lncRNA and is the first locus in the panel where this bites.
It did not change the verdict — the any-gene path names IGF2-AS, which is not IGF2 either — but the
rule is now known to have this failure mode, and changing it is the owner's call.

**The direction is the substantive finding.** The deletion *represses*. Published biology is the
opposite on the maternal allele: the unmethylated ICR binds CTCF and blocks the enhancers from
reaching IGF2, so removing it should *raise* IGF2. The model has no allele and no methylation state,
so it cannot express what the ICR does, and it did not. The registration said this in advance. A
miss here is evidence against the question, not against the model.

`genomeos/benchmark/loci.py` (`MODEL_WINDOW`, `read_reach`), `scripts/loci_stated_intervals.py`,
results `loci_stated_intervals` and `loci_benchmark`. One model request was spent.

---

## 19. Nine more loci, one request, and the last claim standing falls over (2026-09-17)

The roadmap item was "the known-locus benchmark with more loci, not more readings of the one it
has" - the conclusion sections 8 and 16 both reached. Eleven candidates were assembled from
published **perturbation** experiments (mouse or human enhancer deletion, CRISPR interference,
transgenic reporter, or a disease-causing point mutation in the element), never from association
alone; none repeats a panel locus. They are a **separate, separately registered panel**, scored by
`genomeos/benchmark/loci.py` with no change to any scorer or hit rule, and their rates are reported
beside the seventeen and never pooled with them.

`genomeos/benchmark/loci_candidates.py` (the expectations and `PREREGISTRATION`),
`scripts/loci_candidates.py`, `tests/test_loci_candidates.py`, results `loci_candidates` and
`loci_candidate_intervals`. **One model request was spent.**

### The candidates, and what the free filter did to them first

| locus | element | published target | tissue | distance | perturbation |
|---|---|---|---|---|---|
| RET_MCS97 | chr10:43,086,358-43,086,858 (rs2435357) | RET | enteric neural crest | 9.5 kb | mouse enhancer deletion (Chatterjee 2016) |
| IRF6_MCS97 | chr1:209,815,675-209,816,175 (rs642961) | IRF6 | oral periderm | 9.8 kb | transgenic reporter (Fakhouri 2014) |
| MYB_HMIP2 | chr6:135,097,000-135,098,000 (rs66650371) | MYB | erythroid | 83.8 kb | element deletion in erythroid cells (Stadhouders 2014) |
| SORT1_1p13 | chr1:109,274,718-109,275,218 (rs12740374) | SORT1 | liver | 123 kb | AAV in mouse liver (Musunuru 2010) |
| IL2RA_CaRE4 | chr10:6,052,484-6,052,984 (rs61839660) | IL2RA | stimulated T cell | 9.6 kb | CRISPRa tiling and deletion (Simeonov 2017) |
| GDF5_GROW1 | chr20:35,319,108-35,319,608 (rs6060369) | GDF5 | growth plate | 118.9 kb | humanised mouse allele (Capellini 2017) |
| PTF1A_enhancer | chr10:23,219,045-23,219,445 (**stated**) | PTF1A | embryonic pancreas | 26.9 kb | human recessive mutations (Weedon 2014) |
| PAX6_SIMO | chr11:31,664,197-31,664,597 | PAX6 | eye | 146.9 kb | de novo mutation, reporters (Bhatia 2013) |
| HBA_HS40 | chr16:113,381-113,908 | HBA1, HBA2 | erythroid | 63.0 kb | human and mouse deletions (Higgs 1990, Hay 2016) |
| POU3F4_DFN3 | anchored at -900 kb | POU3F4 | inner ear | 900 kb | human microdeletions (de Kok 1996) |
| MYC_BENC | anchored at +1.7 Mb | MYC | haematopoietic | 1,700 kb | CRISPRi tiling, module deletion (Fulco 2016, Bahr 2018) |

Every published direction in the table is *activates*. Six elements are placed on an rsID resolved
from Ensembl GRCh38 and re-checked against Ensembl on every run; two (HBA_HS40, PAX6_SIMO) are
Ensembl assembly-map liftovers of a published GRCh37 interval; PTF1A's is anchored on the published
25 kb offset from GENCODE's gene end, and that anchor lands within 500 bp of a liftover of the
published mutation positions.

**Two candidates died at the reach filter, and the number that matters is the other nine.**
`loci.read_reach` was called - not reimplemented - and it is a gene-**body** test: POU3F4 (900 kb)
and MYC's blood enhancer cluster (1,693 kb measured) have no body anywhere in the model's input, so
they are unaskable and no request was spent on them. Both are also **exhibits only**: their
elements are placed on a published offset with no citable coordinate, which section 13 withdrew a
case for, so they are excluded from every scored rate whatever the filter says.

**Nine of eleven published enhancer-gene links sit well inside the model's 524 kb reach.** The
panel's three megabase cases (SHH, SOX9, IRX5) are not what a published enhancer usually looks like;
they were chosen for being remarkable. A filter that is worth having at the panel costs almost
nothing here, which is itself the measurement.

### What the requests were, and why there was only one

Every chromosome is swept, so an element some annotation already drew has **already been deleted and
the answer is cached**. The plan checks that before anything is spent:

| loci | cCREs over the published element | already deleted | requests |
|---|---|---|---|
| RET, IRF6, SORT1, IL2RA, GDF5, HBA (2 each), MYB, PAX6 (1 each) | 1 to 2 | 1 to 3 | **0** |
| PTF1A_enhancer | **0** | 0 | **1** |
| POU3F4_DFN3, MYC_BENC | - | - | **0** (unaskable) |

**One request, for the whole widening.** Section 6 spent 2,283 on twelve loci; section 18 spent one.
Once the sweep is finished, the cost of a new locus is the cost of the intervals nobody annotated,
and that is 1 of 9 here against 4 of 17 at the panel - because the panel's uncovered elements are
its megabase and coding loci, and ordinary published enhancers are cCREs.

### The scores, beside the panel's

| field | candidates | the panel |
|---|---|---|
| right target, derived | **8/9 (0.889)** | 15/17 (0.882) |
| right target, where the model could answer | **8/9 (0.889)** | 13/15 (0.867) |
| right target, nearest TSS in node | **3/9 (0.333)** | 11/17 (0.647) |
| right target, annotation lookup | **1/9 (0.111)** | 8/17 (0.471) |
| chance floor (random coding gene in the window) | 1.06/9 (0.118) | 7.12/17 (0.419) |
| direction, where a deletion names the target | **3/5 (0.600)** | 5/5 |
| right cell or tissue, where judged | 4/6 (0.667) | - |
| a published class read | 9/9 | 13/17 |

The floors are not comparable: the candidate windows are a uniform 1 Mb and hold more coding genes,
so drawing one at random hits less often. Within each set the derived rate clears its own floor.

**Which layer does the work is not the same here.** The element's own deletion names the published
target first at **5 of 9** (RET, IRF6, IL2RA, GDF5, HBA2); the summed window input at 5 of 9; and
**GTEx eQTLs at 1 of 9**, against the panel where eQTLs carried several loci. The reason is the
tissue list: enteric neural crest, oral periderm, growth plate, embryonic pancreas and a stimulated
T cell are not GTEx tissues. A benchmark built on loci GTEx can see will overrate eQTLs.

**The nearest-TSS heuristic falls from 0.647 to 0.333, and it falls into every trap it was set: 5 of
5** (IRF6 to a novel GENCODE gene 37 bp closer than IRF6, MYB to HBS1L, SORT1 to PSRC1, GDF5 to
FAM83C, HBA to NPRL3). The traps were written down before scoring, and they are geometry, not
results. So a set of ordinary published enhancers is *harder* for the nearest-gene rule than the
panel is, not easier - which was the opposite of what was registered.

### Per locus

| locus | derived hit, by | the element's own deletion | heuristic | direction | notes |
|---|---|---|---|---|---|
| RET_MCS97 | deletion, gene_input | **RET**, -0.45 | RET | right | the enhancer is in RET's own first intron: nothing is at stake |
| IRF6_MCS97 | deletion, eqtl | **IRF6**, -0.57 | trap | right | the only eQTL hit in the set |
| MYB_HMIP2 | gene_input | names nothing | trap (HBS1L) | not judged | one element scored and it moved no gene; the summed window names MYB over HBS1L |
| SORT1_1p13 | **miss** | CELSR2, -2.53 | trap (PSRC1) | not judged | see below |
| IL2RA_CaRE4 | deletion, gene_input | **IL2RA**, +0.24 | IL2RA | **wrong** | deleting the enhancer *raises* IL2RA in the model |
| GDF5_GROW1 | deletion | **GDF5**, +0.15 | trap (FAM83C) | **wrong** | the same sign error at a 119 kb reach |
| PTF1A_enhancer | gene_input | C10orf67, -0.13 | PTF1A | not judged | the one stated interval and the one request; see below |
| PAX6_SIMO | gene_input | names nothing | node names nothing | not judged | a 147 kb reach recovered by summing, not at the element |
| HBA_HS40 | deletion | **HBA2**, -0.77 | trap (NPRL3) | right | the sharpest result of the set |

**HBA_HS40 is the result this widening was for.** Deleting the alpha-globin enhancer names HBA2
first and HBA1 second, 60 kb away, across NPRL3 - the gene the enhancer sits inside and the
nearest-TSS trap. It is the same shape as the ZRS at a twentieth of the distance, and it is the
first time in this benchmark that an element-level deletion has reached a distant published target
past the gene it lives in. The ZRS could never have done it: SHH is not in the model's input.

**The direction is right at 3 of 5, not 5 of 5.** IL2RA and GDF5 both move the *wrong way*: deleting
a published activator raises its target in the model, weakly (+0.24 and +0.15). The panel's 5 of 5
was five loci; this is five more, and pooled the direction reading is 8 of 10. The sign is not
settled.

**SORT1 is the one target miss, and it is the FTO pattern exactly.** SORT1 is second in every
derived layer - second in the eQTLs behind CELSR2, second in the summed window behind ELAPOR1 - so
the lenient reading is right and the strict one is wrong, as at FTO. The element's own deletion
names **CELSR2 at -2.53**, by far the largest effect anywhere in this set, and CELSR2 is the gene
the element sits inside. The published paper reports coordinated changes at SORT1, CELSR2 and PSRC1
and concludes SORT1 is causal for LDL; the expectation named SORT1 alone, before scoring, and the
miss stands.

**The element-level deletion names what it sits inside. This is now three for three.** The ZRS named
LMBR1, the H19 ICR named MIR675 and the IGF2 read-through, and the SORT1 enhancer names CELSR2. Each
time the containing or overlapping transcript takes the top rank and the published target sits below
it. That is not a coincidence of three loci; it is what deleting a piece of a transcript does to that
transcript's own read counts, and it is the same mechanism as the `predicted_coding`-first defect
one level up.

### PTF1A: the one request, and it went to the gene next door

The prediction was registered before the request. The deletion named **C10orf67** first among coding
genes at -0.13, with **PTF1A second** at a mean of +0.003 and its largest drop in cardiac muscle
cell. Every effect is at the noise floor (confidence 0.126). The published tissue is an embryonic
pancreatic progenitor; the model's four cell tracks are K562, HepG2, GM12878 and IMR-90, and its
strongest tissue here was right cardiac atrium. **A miss, and the honest reading is that the question
was asked in the wrong cell**, not that the element does nothing. The summed window names PTF1A
first, so the locus scores as a derived hit through `gene_input` and as an element-level miss.

The `predicted_coding`-first defect that H19 exposed (a non-coding published target hidden behind a
coding read-through) **is not exercised by this set at all**: every candidate target is protein
coding. It is untested here, not absent.

### The controls, standardised - and every claim falls over once coverage is a stratum

45 matched windows, five per graded locus, the same four covariates and the same keep-outs (extended
to keep out the candidate loci too). `genomeos/compare.py` does the standardising; nothing here
reimplements it.

| claim | candidates | controls, covariates only | p | controls, **coverage held fixed** | p |
|---|---|---|---|---|---|
| a derived layer names *some* target | 1.000 | 0.984 | 0.35 | **1.000** | 0.50 |
| the reader has it open in some cell | 0.778 | 0.426 | **0.051** | **0.722** | 0.39 |
| a direction is produced | 0.778 | 0.366 | **0.026** | **0.889** | 0.74 |
| a value sits on constrained sequence | 0.111 | 0.132 | 0.55 | **0.167** | 0.63 |

**Standardising on GC, distance and constraint does not remove the coverage artefact**, because none
of those covariates says whether anybody ever spent a request on the window. Add "has an element
inside it been deleted" as a fifth stratum and the direction claim goes from +0.41 (p=0.026) to
**-0.11**: the matched windows produce a direction *more* often than the candidates do. Section 8
retracted that claim at the panel on 2026-09-15; it is reproduced here on an independent set of
loci, and this time it survived a standardisation, which is the new part. No target was dropped for
want of a control in either reading, and `imbalance` is reported beside both: the constrained
fraction is the one covariate that differs materially between the arms (1.69x).

**The cell claim is a coverage artefact too, and that is new.** At p=0.051 on the matched covariates
it looked like the second real discriminator; holding coverage fixed takes it to +0.06 (p=0.39).
Nobody had conditioned the cell claim before.

### The one claim that survived the panel does not survive this panel

**A value on constrained sequence: 1 of 9 candidates (0.111) against 3 of 45 controls (0.067), and
below the controls once standardised (-0.02, and -0.06 holding coverage fixed).** The registration
predicted about half against about a sixth, on the strength of the panel's 9/17 against 15/85. It
did not happen.

The obvious explanation is element length, and it is wrong. Restricted to elements of 1 kb or less,
the panel reads a value on syntax at **4 of 7** and the candidates at **1 of 9**, and six of the nine
candidate elements are centred on a famous common variant exactly as the panel's short ones are. What
differs is only which variants the three people we hold happen to carry on mammal-constrained bases.

**So the benchmark's last surviving discriminator does not reproduce on a second set of loci.** With
n=9 this is not proof that the panel's 9/17 was noise, and it is not offered as one. It is enough
that the claim can no longer be quoted as established. What would settle it is a third set with
element lengths matched to the panel's, and more than three genomes behind the "value" half of the
reading - the HPRC panel, which was already the recorded cost of the value-domain gap.

### A defect found in this benchmark's own lookup, and what it changes

`Chromosome.ccres_in` bisects the cCRE list, and `load_ccres` **does not return it in start order**,
so the bisect silently missed elements that begin before the window. Measured over the seventeen: it
reported 1 of the 12 cCREs over the MC1R window, 0 of the 1 over the ZRS, 0 of 2 at MYC and 0 of 4 at
ABO. The lookup now sorts its own copy; the domain model is still handed the list exactly as it came,
so nothing about the nodes moves.

**This changes a published count in section 10.** "No cCRE over the published element" is **4 of 17**
(APP, TP53, SOX9_PierreRobin, H19_ICR1), not 7 of 17: ABO, MYC_8q24 and **SHH_ZRS** all have one. So
the sentence "the flagship long-range element is not an ENCODE cCRE" is wrong - the ZRS is covered by
one, and the VISTA run simply reached it first. The conclusion of section 10 is unchanged in
direction (the element-level reading is still bounded by which intervals somebody else drew) and
smaller in size. The panel's committed result was not rerun for this: which of its numbers to refresh
is the benchmark owner's call, and `data/results/loci_benchmark.json` still holds the old count.

### What is left undone

- The candidates are **not** in `loci.PANEL` and the headline stays 15/17. Whether they should be
  merged, and on what denominator, is the owner's call; everything needed to merge them is in
  `loci_candidates.CANDIDATES`.
- The panel has not been rerun since the `ccres_in` fix, so section 10's 7 of 17 is still what the
  committed result says.
- Nine loci is the same size as the evidence that produced the claim it failed to reproduce. A third
  set, length-matched, is the next thing worth a day.
- Every candidate here is an **activator**. Published repressive elements with a perturbation behind
  them were looked for and none was clean enough to state, so the direction axis is only ever tested
  in one direction, at both panels.

---

## 20. The length artefact that the design had already ruled out, and the arm that answers the other question (2026-09-17)

`genomeos/benchmark/loci_length.py`, `scripts/loci_length_controls.py`,
`tests/test_loci_length_controls.py`, result `data/results/loci_length_controls.json`.
**Zero model requests.** 183 seconds, all of it range reads over public tracks.

**This lane is not blinded, and a reader meets that before the numbers rather than after them.** It
was set up as a blinded registration against the value-on-syntax claim, and the blind failed on
timing, not on scope: a paragraph of section 9 above quotes section 19's results, so the lane had a
second lane's figures for this claim, and that lane's own remark about element length, before the
blind was widened to the whole file. No widening can un-read that. What survives is worth stating
plainly and is not nothing: the registration was committed before any number in it existed
(`dbef7b0`), its decision rule is mechanical and symmetric between the outcomes, and its control
selection is seeded and executes rather than being described. That is a **pre-registered unblinded
test** — a weaker instrument than a blinded one and a much stronger one than an unregistered read.

### The result, ahead of every number: the question was already answered by the construction

The suspicion was that "a value sits on constrained sequence" separates the panel from its controls
only because published elements are longer than matched windows, and a longer window is likelier to
overlap constrained sequence for reasons that have nothing to do with being a regulatory element.

`loci.candidate_windows` computes `length = expect.element[1] - expect.element[0]` and draws its
window as `s .. s + length`. **The existing 85 controls are matched to their locus's element length
exactly, window by window: 85 of 85, worst difference 0 bp.** Verified here from the committed
coordinates rather than from a reading of the code (`length_match_check`), and independently by the
coordinating session. Every arm's length distribution is therefore the same five numbers:

| arm | n | min | q1 | median | q3 | max |
|---|---|---|---|---|---|---|
| A panel | 17 | 400 | 400 | 3,000 | 10,766 | 98,237 |
| B existing controls | 85 | 400 | 400 | 3,000 | 10,766 | 98,237 |
| C length-only controls | 85 | 400 | 400 | 3,000 | 10,766 | 98,237 |

So 9 of 17 against 15 of 85 **cannot be a naive length artefact** — not because a test looked and
failed to find one, but because the design makes one impossible. The check was two commands. A
premise nobody had checked survived a brief, a peer's suggestion and a registration before anyone
ran it, which is the transferable part of this section.

### What arm C is for, then: the complement question

Not "does length explain the excess" but **what does matching on length alone buy?** Arm C is five
windows per locus of exactly the element's length, on the same chromosome, outside every panel locus
with 200 kb of flank and outside VISTA elements and GWAS hits with 5 kb — and matched on **nothing
else**: no GC, no distance, no constraint. Seed 4703, 5 of 5 drawn at every one of the seventeen
loci. Nobody should read it as the original test.

### The three arms, with their covariates and their coverage

| | A panel | B existing | C length-only |
|---|---|---|---|
| n | 17 | 85 | 85 |
| **a value on constrained sequence** | **9/17 = 0.529** | **15/85 = 0.176** | **13/85 = 0.153** |
| the same, where the window holds a value at all | 9/13 = 0.692 | 15/62 = 0.242 | 13/60 = 0.217 |
| length median | 3,000 | 3,000 | 3,000 |
| GC median | 0.455 | 0.460 | **0.393** |
| distance to a coding TSS, median | 18,718 | 19,204 | **135,634** |
| constrained fraction, median | 0.0737 | 0.0653 | **0.0159** |
| an element inside was deleted | 14 | 60 | 38 |
| **never looked** (no element scored) | 3 | 25 | 47 |

Arm A reproduces the published 9/17 and arm B the published 15/85 exactly, which is the check that
this lane is reading the same thing the benchmark reads. Arm C is badly unbalanced on distance — its
windows sit 7.2 times further from a coding TSS than the panel's — which is what dropping the
matching costs and why its verdict is taken on the standardised reading, not the raw one.

**The two zeros are counted apart, because one of them is not a measurement.**

| category | A panel | B existing | C length-only |
|---|---|---|---|
| `never_looked_syntax` (no reading came back) | 0 | 0 | 0 |
| `no_variable_position` (looked; nobody in the trio differs here) | 4 | 23 | 25 |
| `values_none_constrained` (looked; values present, none constrained) | 4 | 47 | 47 |
| `values_in_syntax` (the hit) | 9 | 15 | 13 |

A quarter of arm B's zeros and nearly a third of arm C's are windows where the claim had nothing to
be true **of**. That is why the row "where the window holds a value at all" is printed beside the
headline: on that denominator the panel reads 9/13 against 15/62 and 13/60.

### The comparison both ways, covariates only and coverage held fixed

`genomeos.compare.standardised`, direct standardisation, the control arm carrying the matched n.
Arm B decides on `P1_length`; arm C decides on `S1_covariates`, because arm C is matched on nothing
but length and a low rate there could be its distance to a TSS rather than its length. Both
deciding rules were fixed before scoring; the change of deciding stratification for arm C was a
registration revision made in the open and committed (`88d6033`) ahead of the run.

| stratification | panel against B | panel against C |
|---|---|---|
| raw, unstandardised | 0.3529 | 0.3765 |
| `P1_length` (length) | **0.3529**, p 0.010, n 17, dropped 0 | 0.3765, p 0.006, n 17, dropped 0 |
| `P2_length_coverage` (**+ coverage**) | **0.3578**, p 0.009, n 17, dropped 0 | 0.3381, p 0.014, n 17, dropped 0 |
| `P3_length_has_values` (+ the claim's own denominator) | 0.3657, p 0.008, n 17, dropped 0 | 0.3873, p 0.004, n 17, dropped 0 |
| `S1_covariates` (length, GC, distance) | 0.3557, p 0.010, n 17, dropped 0 | **0.2556**, p 0.062, n 14, dropped 3 |
| `S2_covariates_coverage` (**+ coverage**) | 0.3669, p 0.007, n 17, dropped 0 | **0.1421**, p 0.244, n 11, dropped 6 |

**Against the existing controls, coverage does not touch this claim.** 0.3529 without the coverage
stratum, 0.3578 with it; 0.3557 on the full covariate set, 0.3669 with coverage added. That is the
opposite of what coverage did to the direction claim in sections 8, 17 and 19, and it is worth
saying why: a deletion has to be *spent* on a window before that window can produce a direction,
while a value on constrained sequence is read from the trio's own variants and a public phyloP
track, which cover every window whether or not anybody chose it. The coverage denominator that has
withdrawn three readings here is a denominator this one does not have.

**Against the length-only controls it is a different story, and that is the finding of arm C.** The
gap is 0.3765 when only length is held fixed, 0.2556 once GC and distance are standardised back in,
and 0.1421 with coverage beside them — two thirds of it gone, and the matched n falls from 17 to 14
to 11 as strata empty. Matching on length alone buys a control arm that *looks* well separated and
is not: most of that separation is its distance to a coding TSS.

### The verdict, in its registered words

> **undecided at this n**

which is the outcome the registration named in advance as the expected one. It fails "the claim
survives" at arm C on two counts — p 0.062 against the registered 0.05, and 0.1421 with coverage
held fixed against the registered 0.15 — and it does not reach "the claim is an artefact of length"
at either arm, which would have needed a raw excess of 0.25 or more collapsing to 0.10 or less.

The two branches are **not equally hard, on purpose**. Survives needs both arms, two
stratifications, a p-value and a matched-n floor; artefact needs one arm and two point estimates.
The strong form of the negative — bounding the difference below 0.10 with 95% confidence — needs a
standard error near 0.05, which is about 100 published loci. The panel has 17, and at 17 the
standard error of the standardised difference is about 0.13, so nothing below about 0.21 separates
from zero. A durable "survives" at an attenuated effect of 0.20 would need about 35 loci. Those
numbers were in the registration before the run, not fitted to it afterwards.

### The combination rule, fixed before the number existed

There are now three figures about this claim: the panel's own 9/17, section 19's on nine independent
loci, and this section's arms. **They are never pooled, and no "two of the three agree" reading is
to be taken.** The binding reason is arithmetic rather than judgement: **this section's arm A *is*
the panel's seventeen positives**, so pooling it with 9/17 counts the same seventeen observations
twice. Section 19's loci come from a different sampling frame chosen by a different rule, and arms B
and C are control arms, which are not estimates of the claim's rate at published loci at all. They
are reported side by side with their frames named, and a disagreement between them stays a
disagreement.

### What this section does and does not change

It does **not** rehabilitate the value-on-syntax row and it does not retire it. Section 19's reading
stands where it is, on its own loci, and this one stands here, on the panel's. What is settled is
narrower and firmer than either: the **length** explanation is off the table by construction, the
**coverage** explanation — which has taken down every other claim this benchmark scores — does not
touch this one, and what is left at n = 17 is undecided and will stay undecided until the panel is
about twice its present size.

### Read after scoring: where this agrees with section 19 and where it does not

Section 19 was read only after the verdict above was committed. Three things came out of the
comparison and they are kept apart rather than reconciled.

**1. Agreement, by a stronger route.** Section 19 says "the obvious explanation is element length,
and it is wrong", and argues it from a sub-reading: restricted to elements of a kilobase or less the
panel reads 4 of 7 and the candidates 1 of 9. That is two positive arms of seven and nine compared
with each other. This section reaches the same conclusion from the construction — the panel's
controls are length-matched exactly, 85 of 85 at 0 bp — which does not depend on a sample size at
all. Same answer, much better warrant, and the warrant is a property of the code rather than of a
number.

**2. A disagreement about the remedy, and it matters for what gets built next.** Section 19's closing
recommendation is "a third set with element lengths matched to the panel's". Read as *control*
windows, that set already exists and always did, and this section built a second one; neither
settles anything. Read as *loci*, the recommendation is sound but it is a statement about
comparability **between the two panels**, not about the panel-versus-control comparison that produced
9/17 against 15/85. The two must not be run together: nothing about matching will settle this, only n
will. About 35 loci for a durable "survives", about 100 for the strong negative.

**3. A genuine disagreement about coverage, and this one is a finding.** Section 19's heading is
"every claim falls over once coverage is a stratum", and at the candidate panel that holds for all
four claims it scores. **At the seventeen-locus panel it does not hold for this claim.** Holding
`deletion_scored` fixed moves the value-on-syntax difference from 0.3529 to 0.3578, and from 0.3557
to 0.3669 on the full covariate set — no movement at all, where the direction claim went from +0.41
to −0.11 on the same manoeuvre. The mechanism is the reason, and it generalises: **a coverage
stratum only bites a claim that has to be bought.** A direction needs a deletion spent on that
window; a value on constrained sequence is read from the trio's own variants and a public phyloP
track, which cover every window whether or not anyone chose it. Section 19's generalisation —
standardising on covariates is not a substitute for conditioning on coverage — stands, and the
corollary this section adds is that the coverage stratum is not a universal solvent: it is the right
test for a bought claim and a null manoeuvre for a free one, and reporting it both ways is how you
tell which kind you have.

**A post-hoc reading, labelled as post-hoc because it was not registered, that reconciles the two
panels' numbers.** Length does not confound the panel-versus-control comparison, but it is
overwhelmingly the strongest thing in this claim. Over the 170 control windows of arms B and C
pooled:

| window length | a value on constrained sequence | holds a value at all | mean values per window |
|---|---|---|---|
| under 1 kb | **1/70 = 0.014** | 32/70 | 1.0 |
| 1 to 5 kb | 1/30 = 0.033 | 24/30 | 3.6 |
| 5 kb and over | **26/70 = 0.371** | 66/70 | 60.3 |

A 26-fold gradient, and it is arithmetic: a 60-fold longer window holds 60 times as many variable
positions and needs only one of them to land on a constrained base. That is why **the two panels'
figures are not on the same scale.** Section 19's candidate elements are almost all 400 to 1,000 bp,
so its 1/9 and its controls' 3/45 = 0.067 both sit in the regime where the control rate here is
0.014; the panel's 15/85 = 0.176 is lifted by its long elements. Comparing 9/17 with 1/9 directly
compares two different length regimes, which is section 19's own point about matching, arrived at
from the control side. Within the short regime the panel reads 4 of 7 against its own length-matched
controls at 1 of 70 — a wide separation at n = 7, which is offered as a reading and not as a verdict,
because it was computed after the number was known and the registration says what that is worth.

### Left undone

- The registered question needs **loci**, not arms: about 35 for a durable positive, about 100 for a
  negative that bounds the difference below 0.10. No amount of control construction substitutes.
- The "value" half of the reading still rests on three genomes. The length gradient above is partly a
  statement about how few variants the GIAB trio contributes to a short window, and the HPRC panel is
  the recorded cost of fixing that.
- The post-hoc length gradient is unregistered and should be re-run as a registered reading if anyone
  wants to quote it.

## 21. A third set, the direction axis asked both ways, and an answer that was registered as undecidable (2026-09-21)

`genomeos/benchmark/loci_third.py`, `scripts/loci_third.py`, `tests/test_loci_third.py`, results `loci_third` and `loci_third_intervals`. **3 model requests spent.** 237 seconds of range reads over public tracks for the rest.

Section 19's own "left undone" named three gaps. This set is aimed at two of them and records the third as unfillable. **Every one of the nine candidates activates**, so the benchmark has never asked whether the deletion layer's sign means anything when the published element's job is to hold a gene *down*; and section 20 asked for a set whose element lengths span the panel's rather than sitting at a uniform 500 bp. Two of the nine loci here are published **repressors** with a perturbation behind them, and the element lengths run 300 bp to 52,000 bp.

It is a **third sampling frame**, separately registered, scored by `genomeos/benchmark/loci.py` with no change to any scorer or hit rule. Its rates are reported beside the seventeen and the nine with each frame named, and never pooled with either.

### The nine, and what the free filter did to them first

| locus | element | len | published target | direction | distance | perturbation |
|---|---|---|---|---|---|---|
| SHH_SBE2 | chr7:156,268,356-156,269,154 | 798 | SHH | **activates** | 455,893 bp | point mutation in a holoprosencephaly patient; transgenic reporter (Jeong 2008) |
| GATA2_plus95 | chr3:128,483,458-128,483,958 (**stated**) | 500 | GATA2 | **activates** | 9,500 bp | targeted mouse enhancer deletion; germline human mutations (Johnson 2012) |
| TAL1_MuTE | chr1:47,239,475-47,239,975 (**stated**) | 500 | TAL1 | **activates** | 7,500 bp | CRISPR deletion of a somatically created enhancer (Mansour 2014) |
| CDKN2A_9p21 | chr9:22,124,228-22,124,728 | 500 | CDKN2A, CDKN2B | **activates** | 114,903 bp | 70 kb mouse interval deletion (Visel 2010) |
| KITLG_blond | chr12:88,934,308-88,934,808 | 500 | KITLG | **activates** | 353,457 bp | human-enhancer knock-in mouse (Guenther 2014) |
| SOST_VanBuchem | chr17:43,666,737-43,718,737 (**stated**) | 52,000 | SOST | **activates** | 35,000 bp | homozygous 52 kb human deletion; mouse enhancer deletion (Balemans 2002) |
| HBG1_BCL11A_site | chr11:5,249,724-5,250,224 (**stated**) | 500 | HBG1, HBG2 | **represses** | 115 bp | HPFH point mutations; editing in human erythroid cells (Martyn 2018) |
| MYC_PVT1promoter | chr8:127,794,262-127,794,762 (**stated**) | 500 | MYC | **represses** | 58,829 bp | CRISPR deletion of the promoter raises MYC (Cho 2018) |
| TERT_promoter | chr5:1,294,950-1,295,250 (**stated**) | 300 | TERT | **activates** | 124 bp | recurrent somatic mutation; CRISPR reversal (Chiba 2015) |

**0 of 9 died at the reach filter, and that was the registered prediction.** `loci.read_reach` was called, not reimplemented: a gene-**body** test against the scorer's 1,048,576 bp input, never a TSS-distance test. Beside the other two frames - 3 of 17 unaskable (SHH_ZRS, SOX9_PierreRobin, FTO_IRX3's IRX5), 2 of 11 died, both exhibits placed on a bare offset - the filter keeps killing the famous megabase loci and almost nothing else. Three frames have now made the same measurement: a published enhancer-gene link is usually well inside the model's reach.

**The sharpest case is SHH.** Section 18 retired the ZRS as an element-level miss because SHH is 979 kb away and **not in the model's input at all** - no number of requests puts it there. SBE2 is a different published SHH enhancer 455,893 bp away, SHH's body is inside the window, and the question the benchmark could not ask about its flagship gene becomes askable by changing the element rather than the model.

### What the requests were

| loci | cCREs over the element | already deleted | requests |
|---|---|---|---|
| CDKN2A_9p21, GATA2_plus95, HBG1_BCL11A_site, KITLG_blond, SHH_SBE2, SOST_VanBuchem | 1 to 35 | 1 to 30 | **0** |
| MYC_PVT1promoter, TAL1_MuTE, TERT_promoter | 1 to 2 | 0 | **1 each** |

**3 requests for the whole set**, because the finished all-chromosome sweep had already deleted an element inside 6 of the 9 published intervals. Section 6 spent 2,283 on twelve loci; the cost of a new locus now is the cost of the intervals nobody annotated.

**A caveat on those three, from the run's own record rather than from reading it charitably.** `Context.score_region` substitutes an overlapping ENCODE element when there is one, and it did so at MYC_PVT1promoter, TAL1_MuTE, TERT_promoter - every interval this set paid for. So what was deleted is an annotated element *overlapping* the published one, not the published interval itself, and the stated coordinate acted as a pointer rather than as the thing scored. It is recorded per locus in `loci_third_intervals` and is the honest reading of what those three requests bought.

### The scores, beside the other two frames

| field | the third set | the nine candidates | the seventeen panel |
|---|---|---|---|
| right target, derived | **8/9 (0.889)** | 8/9 (0.889) | 15/17 (0.882) |
| right target, where the model could answer | **8/9 (0.889)** | 8/9 (0.889) | 13/15 (0.867) |
| right target, nearest TSS in node | **7/9 (0.778)** | 3/9 (0.333) | 11/17 (0.647) |
| right target, annotation lookup | **3/9 (0.333)** | 1/9 (0.111) | 8/17 (0.471) |
| chance floor (random coding gene in the window) | 2.4/9 | 1.06/9 | 7.12/17 |
| direction, where a deletion names the target | **6/7 (0.857)** | 3/5 (0.600) | 5/5 |
| right cell or tissue, where judged | 0/4 (0.000) | 4/6 (0.667) | - |
| a published class read | 7/9 | 9/9 | 13/17 |

The floors are not comparable across the frames and are printed so that each rate can be read against its own: this set's 2.4 of 9 sits between the candidates' and the panel's. Within each frame the derived rate clears its own floor.

**The three derived rates agree to within 0.007**, which is the least interesting thing here and worth saying first so it cannot be mistaken for the finding. Three sets chosen on different days by different criteria land on 8/9 (0.889), 8/9 (0.889) and 15/17 (0.882). The headline is stable across sampling frames. Everything that moves is underneath it.

### The nearest-gene rule moves further than anything else, and in the opposite direction to section 19

The heuristic reads 7/9 (0.778) here against 3/9 (0.333) at the candidates and 11/17 (0.647) at the panel. **That spread is wider than any spread in the derived rate, and it is a property of the loci, not of the rule.** It was registered in advance as a weakness of this set: 7 of 9 loci here have a published target as their own nearest coding TSS, because an element that sits in its target's intron or promoter is what a clean perturbation experiment usually looks like.

**And the two traps caught nothing, which is not the same as the traps failing.** SHH_SBE2, SOST_VanBuchem were written down before scoring as the two loci where the nearest coding TSS is not the published target - RNF32 at 371,126 bp against SHH's 455,893, MEOX1 at 4,808 bp against SOST's 40,054. The heuristic fell into **0 of them**: at both loci the node model named **no gene at all** rather than naming the wrong one. Its two misses (SHH_SBE2, KITLG_blond) are both silences. A baseline that abstains at the long-range loci and answers at the short-range ones is not the baseline its rate describes, and this is the first set in which that distinction shows.

### The direction axis, asked both ways - and the registered answer is *undecidable*

This is what the set was assembled for, and the registration wrote down in advance what each outcome would be worth: both repressors right is weak support that the layer reads direction, both inverted is weak support that it does not, **one of each decides nothing**.

| arm | right | n |
|---|---|---|
| activators | 5 | 5 |
| **repressors** | **1** | **2** |

**It came back one of each, so by its own registration this set decides nothing about the sign.** No pooled direction rate is quoted here and none should be quoted from it.

| locus | published | the model | |
|---|---|---|---|
| GATA2_plus95 | activates | activates | right |
| TAL1_MuTE | activates | activates | right |
| CDKN2A_9p21 | activates | activates | right |
| SOST_VanBuchem | activates | activates | right |
| HBG1_BCL11A_site (repressor) | represses | activates | **inverted** |
| MYC_PVT1promoter (repressor) | represses | represses | right |
| TERT_promoter | activates | activates | right |

**The two repressors are the two readings worth having, and they disagree with each other.**

At **MYC_PVT1promoter** the model gets it right: deleting the PVT1 promoter *raises* MYC (+0.1392), which is what Cho 2018 reports and is the opposite of what the panel's own MYC element does 335 kb the other side of the gene. A layer that reports a sign has to disagree with itself across those two elements, and it did.

At **HBG1_BCL11A_site** it gets it wrong, and this is the sharper of the two because nothing else about the locus is ambiguous. Deleting the BCL11A motif at -115 *lowers* HBG1 in the model (-1.0195); in human erythroid cells, destroying that motif is one of the best-replicated ways to *raise* fetal haemoglobin. The element is in the reader's own cell type, the target is named first and correctly, the distance is 115 bp - the model has every input it could want and reports the sign backwards.

**What this adds to the ledger.** The candidates got 2 of 5 wrong among activators (IL2RA, GDF5); this set gets 0 of 5 wrong among activators and 1 of 2 wrong among repressors. Pooled across all three frames the sign is right more often than not, and the one clean case where a published *repressor* could be checked with the cell, the target and the distance all in the model's favour came back inverted. **n = 2 was declared too small before the run and is still too small after it.** What would settle it is a CRISPRi screen reporting de-repression, not more curation - published repressive elements with a perturbation behind them were searched for twice now, and the two here are what a second search produced.

### The positive control, and the one miss

**The control passed.** Deleting TERT's own core promoter named TERT first at -5.0034 - by a wide margin the largest effect anywhere in this set. That was registered as the thing whose failure would invalidate the run rather than score as a miss, so it is reported before the rates that depend on it and not after.

**KITLG_blond is the one target miss, and it is the FTO and SORT1 pattern for the third time.** KITLG is *second* in the summed-window layer behind DUSP6 and is named by no other derived layer; the element's own deletion names nothing at all, and neither does the node. So the lenient reading hits and the strict one misses, as at FTO and at SORT1. The expectation named KITLG alone, before scoring, and the miss stands.

**Which layer does the work has moved again.** The element's own deletion names the published target first at **7 of 9** here, against 5 of 9 at the candidates. The difference is not the model: it is that six of these nine elements sit inside or beside their own target, which is the registered weakness of the set and the same fact that lifts the nearest-gene rule.

### The controls, with presence counted before coverage is made a stratum

31 matched windows, the same four covariates and the same hit rules, with the keep-outs extended to every locus in **all three** sets. `genomeos/compare.py` does the standardising; nothing here reimplements it.

**`input_presence` runs first, and it is the reason the table below can be read at all.**

| input | kind | targets | controls |
|---|---|---|---|
| a deletion spent inside the window | **bought** | 9/9 | 17/31 |
| a GTEx eQTL distilled for the window | **bought** | 5/9 | 24/31 |
| the mammalian constraint track over the window | **free** | 9/9 | 31/31 |

A coverage stratum can only bite a claim whose input somebody had to **buy**. The deletion input is bought here - 17 of 31 controls have one - so the direction claim is genuinely under test. The constraint track is **free**, present on every window in both arms, so a null on the storage claim under a coverage stratum is arithmetic and not a test the claim passed. Section 19 did not separate those two cases; this run does, before the strata are run rather than after.

| claim | third set | controls, covariates only | p | controls, **coverage held fixed** | p |
|---|---|---|---|---|---|
| target | 0.778 | 0.875 (-0.125) | 0.742 | **1.000** (-0.250) | 0.949 |
| cell | 0.889 | 0.719 (+0.156) | 0.214 | **0.750** (+0.125) | 0.258 |
| direction | 0.778 | 0.708 (+0.042) | 0.426 | **0.875** (-0.125) | 0.742 |
| storage | 0.111 | 0.412 (-0.287) | 0.915 | **0.412** (-0.287) | 0.915 |

**The coverage artefact reproduces on a third independent set of loci.** Raw, the direction claim separates: 0.778 against 0.484, p = 0.038. Standardising on GC, distance and constraint takes it to +0.042; adding *has an element inside this window actually been deleted* as a fifth stratum takes it to **-0.125** at p = 0.742 - the matched windows produce a direction slightly more often than the loci do. Section 8 retracted this claim at the panel, section 19 reproduced the retraction at the candidates, and it now holds at three frames. It is not a property of which loci anyone picked.

**The value-on-syntax claim does not separate here either:** 0.111 at the loci against 0.194 at the controls, and -0.287 once standardised. That is the candidates' 1/9 shape and not the panel's 9/17, and it was registered in advance as the expected outcome because this set has one long element and eight short ones. **With one long element it cannot settle a length gradient and was declared unable to before the number existed.** What it adds is only that the panel's 9/17 has now failed to reproduce twice.

**The controls are thinner than they should be and that is stated rather than buried.** Five windows per locus were asked for and 31 of 45 were found: the length-and-distance matching could place five for six loci, one for HBG1_BCL11A_site and **none** for MYC_PVT1promoter or TERT_promoter, because a window whose nearest coding TSS is essentially zero away and whose length is 300 to 500 bp is hard to match outside a promoter. One locus is dropped for want of a control in every comparison above (1 of 9), and the constrained fraction is 1.93x between the arms - the one covariate that differs materially, as it did at the candidates.

### Where this set disagrees with the other two, and what the disagreement is

| | the seventeen | the nine | the third set |
|---|---|---|---|
| chosen for | fame | perturbation | direction and length |
| derived target | 15/17 (0.882) | 8/9 (0.889) | 8/9 (0.889) |
| nearest-gene rule | 11/17 (0.647) | 3/9 (0.333) | 7/9 (0.778) |
| annotation lookup | 8/17 (0.471) | 1/9 (0.111) | 3/9 (0.333) |
| unaskable for reach | 3/17 | 2/11 | 0/9 |

**The headline is the stable number and the baselines are the unstable ones.** Across three frames the derived target rate moves by 0.007 and the nearest-gene baseline moves by 0.445. A benchmark reporting a derived rate alone would look reproducible; the thing it is supposed to be measured *against* is what depends on who picked the loci. That is a finding about benchmarks, and it is the main one this set produces.

### Left undone

- **The repressor arm is n = 2 and came back one of each, which was registered as undecidable.** Nothing here licenses a claim about the deletion layer's sign. The HBG inversion is a single reading at a locus where every input favoured the model, and it is worth following up on its own rather than aggregating.
- **A published non-coding target is still missing after three sets**, so the `predicted_coding`-first defect H19 exposed remains untested. The search is recorded in `PREREGISTRATION['aims_that_failed_before_the_run']`: the clean non-coding cases act through parent-of-origin methylation, which no layer here carries.
- **The three paid requests were answered about an overlapping annotated element**, not the stated interval, because `Context.score_region` substitutes one when it can. Whether that substitution should be refused for a stated interval is the benchmark owner's call.
- **`loci.STATED_INTERVAL_RESULTS` does not list `loci_third_intervals`.** Until it does, `loci_third.register_stated_intervals` applies the same one-line change at runtime. The hunk is in `loci_third.NEEDED_LOCI_HUNK`; this lane was not permitted to make it.
- **The cell axis is close to unmeasurable on this set**, as registered: 5 of 9 published tissues have no counterpart in the eleven reader cells or the GTEx panel, and where it could be judged it reads 0/4 (0.000). A benchmark built on loci GTEx can see will keep overrating the tissue layers.
- **The controls are 31 windows, not 45**, and two loci have none at all. A promoter-proximal element cannot be length-and-distance matched against anything that is not itself a promoter, which is a structural limit of `loci.candidate_windows` rather than a shortage of genome.
- **The three sets are still three sets.** Whether any of them should merge, and on what denominator, is the benchmark owner's call; nothing here pools them.

---

## 22. A fourth frame, drawn by a rule instead of chosen, and the headline does not survive it, at n = 50 and in both halves (2026-09-21)

`genomeos/benchmark/loci_fourth.py` (the registration and the rule), `scripts/loci_fourth.py`, `scripts/loci_fourth_section.py`, `tests/test_loci_fourth.py`, results `loci_fourth` and `loci_fourth_intervals`. **4 model requests spent.** 1009 seconds of range reads over public tracks for the rest.

Section 21 ended on a number that reads well and one that does not. Across three frames the derived target rate moves 0.007 and the nearest-gene rule over the same three moves 0.445. Section 21 also said why the baseline moves: at the third set seven of nine published targets are the element's own nearest coding TSS. So the three frames differ mostly in **geometry**, and the question that raises cannot be answered from any of them - is the 0.88 a property of the model, or of the fact that curators pick elements sitting next to their targets?

**A frame chosen by hand cannot answer it, because whoever picks the loci picks the geometry.** So this one was not picked. The rule was written, committed and its tests run before the file it draws from was read (commit 2103c5a), and it contains no loci at all:

> every element in the ENCODE CRISPR benchmark's **held-out** arm whose regulated target is **not** its own nearest coding TSS.

That is the geometry on which the cheap answer is wrong by construction. It reads no model output, no effect size and no chromatin score - only the file's own positive flag and GENCODE arithmetic - and it is adversarial to this benchmark's headline rather than flattering to it. The K562 **training** arm is excluded because `data/results/crispri_benchmark.json` fitted logistic weights on it on 2026-09-16.

### What the rule drew, and what it threw away

| step | elements |
|---|---|
| perturbed elements with a regulated gene in the held-out arm | 175 |
| rejected: the nearest coding TSS IS a published target - the shortcut is right here | 85 |
| rejected: within 200000 bp of a locus in an earlier frame | 16 |
| rejected: a regulated target is not protein coding in GENCODE | 14 |
| **passed the rule** | **60** |
| drawn: every element that passes, the cap raised to 60 | 60 |
| of those, drawn first under the registered cap of 24 | 24 |
| added by the cap raise, same rule | 36 |

**85 of the 175 held-out CRISPR positives have their own nearest coding TSS as the published target** - 49% of them, and 59% of the 145 where the geometry question is even well posed. That number is worth reading before any rate below it: at half the field's own measured enhancer-gene pairs, naming the nearest gene is simply correct. A benchmark drawn without a geometry rule inherits that proportion, whoever draws it.

**14 elements were dropped for having a non-coding regulated target**, and that answers an item section 21 left open. Three frames failed to find a published non-coding target clean enough to state by curation; the held-out arm has that many, and they were dropped here on purpose because the deletion layer ranks `predicted_coding` first, so they test the H19 defect rather than the geometry question. They are a frame of their own, already drawn.

**The cap no longer binds.** It was 24 when the frame was first run (commit 2103c5a) and the first run scored that many, leaving 36 loci that pass the same rule undrawn. Commit f407915 raised it to 60, which is every element the rule passes, so there is now nothing left over for a later hand to choose from. The raise was committed before the second run, with its expected cost and with what each outcome would mean written down first: `loci_fourth.CAP_RAISE`. `assess` and `select` - the rule itself - are untouched, and the two halves are reported separately below because genome order is not random with respect to gene density.

### The 60, and what the free filter did to them first

| locus | element | len | published target | nearest coding TSS | distance | cell | half |
|---|---|---|---|---|---|---|---|
| AGTRAP_K562_chr1_11798k | chr1:11,798,254-11,798,873 | 619 | AGTRAP | MTHFR | 62,478 bp | K562 | first 24 |
| TNFRSF8_K562_chr1_12038k | chr1:12,038,207-12,039,433 | 1,226 | TNFRSF8 | MIIP | 24,555 bp | K562 | first 24 |
| AUNIP_K562_chr1_26377k | chr1:26,377,452-26,379,163 | 1,711 | AUNIP | ZNF683 | 518,848 bp | K562 | first 24 |
| TMEM54_K562_chr1_31928k **(unaskable)** | chr1:31,928,713-31,929,213 | 500 | TMEM54 | ENSG00000288678 | 972,476 bp | K562 | first 24 |
| CCT3_K562_chr1_156119k | chr1:156,119,466-156,120,399 | 933 | CCT3 | LMNA | 218,482 bp | K562 | first 24 |
| DUSP23_K562_chr1_159921k | chr1:159,921,814-159,922,623 | 809 | DUSP23 | TAGLN2 | 141,272 bp | K562 | first 24 |
| TSEN15_K562_chr1_184031k | chr1:184,031,227-184,031,862 | 635 | TSEN15 | COLGALT2 | 20,105 bp | K562 | first 24 |
| PCBP1_WTC11_chr2_70086k | chr2:70,086,003-70,086,502 | 499 | PCBP1 | ENSG00000293615 | 1,199 bp | WTC11 induced pluripotent stem cell | first 24 |
| XPC_K562_chr3_14232k | chr3:14,232,088-14,233,010 | 922 | XPC | LSM3 | 53,876 bp | K562 | first 24 |
| STX18_K562_chr4_5003k | chr4:5,003,097-5,003,889 | 792 | STX18 | CYTL1 | 461,426 bp | K562 | first 24 |
| IL6ST_Jurkat_chr5_56148k | chr5:56,148,470-56,148,970 | 500 | IL6ST | ANKRD55 | 153,726 bp | Jurkat T-lymphoblast | first 24 |
| RGS14_K562_chr5_177390k | chr5:177,390,358-177,390,920 | 562 | RGS14 | SLC34A1 | 32,796 bp | K562 | first 24 |
| GMPR_K562_chr6_15299k **(unaskable)** | chr6:15,299,497-15,299,997 | 500 | GMPR | JARID2 | 938,831 bp | K562 | first 24 |
| MAN1A1_K562_chr6_119137k | chr6:119,137,438-119,137,937 | 499 | MAN1A1 | FAM184A | 212,078 bp | K562 | first 24 |
| STXBP5_K562_chr6_146875k | chr6:146,875,908-146,877,348 | 1,440 | STXBP5 | ADGB | 327,728 bp | K562 | first 24 |
| AP1S1_K562_chr7_101138k | chr7:101,138,234-101,138,832 | 598 | AP1S1 | SERPINE1 | 15,870 bp | K562 | first 24 |
| FAM3C_HCT116_chr7_121301k | chr7:121,301,407-121,302,673 | 1,266 | FAM3C | WNT16 | 94,327 bp | HCT116 colorectal carcinoma | first 24 |
| FAM3C_HCT116_chr7_121303k | chr7:121,303,504-121,304,004 | 500 | FAM3C | WNT16 | 92,613 bp | HCT116 colorectal carcinoma | first 24 |
| MYC_HCT116_chr8_129710k **(unaskable)** | chr8:129,710,072-129,711,049 | 977 | MYC | GSDMC | 1,974,491 bp | HCT116 colorectal carcinoma | first 24 |
| JAK2_K562_chr9_4852k | chr9:4,852,343-4,852,926 | 583 | JAK2 | RCL1 | 132,594 bp | K562 | first 24 |
| DOLPP1_K562_chr9_129742k **(unaskable)** | chr9:129,742,911-129,743,411 | 500 | DOLPP1 | PTGES | 662,057 bp | K562 | first 24 |
| VIM_K562_chr10_17029k | chr10:17,029,282-17,030,231 | 949 | VIM | CUBN | 198,501 bp | K562 | first 24 |
| VIM_K562_chr10_17425k | chr10:17,425,572-17,426,534 | 962 | VIM | ST8SIA6 | 197,794 bp | K562 | first 24 |
| ARID5B_K562_chr10_61765k | chr10:61,765,101-61,766,128 | 1,027 | ARID5B | CABCOCO1 | 135,638 bp | K562 | first 24 |
| ARID5B_K562_chr10_61781k | chr10:61,781,930-61,782,430 | 500 | ARID5B | CABCOCO1 | 119,072 bp | K562 | next 36 |
| EXOC6_K562_chr10_92755k | chr10:92,755,668-92,757,563 | 1,895 | EXOC6 | HHEX | 91,838 bp | K562 | next 36 |
| FADS3_K562_chr11_61833k | chr11:61,833,892-61,835,375 | 1,483 | FADS3, FEN1 | FADS2 | 41,996 bp | K562 | next 36 |
| FADS3_K562_chr11_61841k | chr11:61,841,500-61,842,900 | 1,400 | FADS3 | FADS2 | 49,344 bp | K562 | next 36 |
| CCND1_HCT116_chr11_69298k | chr11:69,298,418-69,298,918 | 500 | CCND1 | MYEOV | 342,435 bp | HCT116 colorectal carcinoma | next 36 |
| CCND1_HCT116_chr11_69352k | chr11:69,352,659-69,353,159 | 500 | CCND1 | MYEOV | 288,194 bp | HCT116 colorectal carcinoma | next 36 |
| NDUFC2_K562_chr11_78290k | chr11:78,290,169-78,290,887 | 718 | NDUFC2 | USP35 | 210,308 bp | K562 | next 36 |
| NDUFC2_K562_chr11_78291k | chr11:78,291,733-78,292,875 | 1,142 | NDUFC2 | USP35 | 212,084 bp | K562 | next 36 |
| SRSF8_K562_chr11_95152k | chr11:95,152,782-95,153,802 | 1,020 | SRSF8 | ENDOD1 | 86,415 bp | K562 | next 36 |
| RBM7_K562_chr11_114216k | chr11:114,216,234-114,217,571 | 1,337 | RBM7 | NNMT | 183,625 bp | K562 | next 36 |
| ENO2_WTC11_chr12_6925k | chr12:6,925,492-6,926,215 | 723 | ENO2 | ATN1 | 11,403 bp | WTC11 induced pluripotent stem cell | next 36 |
| SLC2A3_WTC11_chr12_7995k | chr12:7,995,135-7,995,634 | 499 | SLC2A3 | FOXJ2 | 59,087 bp | WTC11 induced pluripotent stem cell | next 36 |
| PCBP2_K562_chr12_54304k **(unaskable)** | chr12:54,304,142-54,304,642 | 500 | PCBP2 | NFE2 | 852,290 bp | K562 | next 36 |
| ERP29_K562_chr12_111994k | chr12:111,994,126-111,994,726 | 600 | ERP29 | TMEM116 | 18,920 bp | K562 | next 36 |
| ERP29_K562_chr12_111995k | chr12:111,995,026-111,995,626 | 600 | ERP29 | TMEM116 | 18,020 bp | K562 | next 36 |
| DENR_K562_chr12_121811k **(unaskable)** | chr12:121,811,825-121,813,015 | 1,190 | DENR | SETD1B | 940,399 bp | K562 | next 36 |
| HSP90AA1_K562_chr14_102551k | chr14:102,551,749-102,552,488 | 739 | HSP90AA1 | RCOR1 | 464,942 bp | K562 | next 36 |
| SMAD6_K562_chr15_66632k | chr15:66,632,242-66,633,128 | 886 | SMAD6 | LCTL | 69,649 bp | K562 | next 36 |
| GFER_K562_chr16_1828k | chr16:1,828,755-1,829,254 | 499 | GFER | FAHD1 | 155,143 bp | K562 | next 36 |
| ECI1_K562_chr16_2140k | chr16:2,140,522-2,141,021 | 499 | ECI1 | PKD1 | 110,829 bp | K562 | next 36 |
| SEPHS2_Jurkat_chr16_30472k | chr16:30,472,349-30,472,911 | 562 | SEPHS2 | ITGAL | 26,654 bp | Jurkat T-lymphoblast | next 36 |
| TRAPPC2L_K562_chr16_88496k | chr16:88,496,084-88,496,583 | 499 | TRAPPC2L | ZFPM1 | 360,751 bp | K562 | next 36 |
| ZNF276_K562_chr16_88769k **(unaskable)** | chr16:88,769,534-88,771,459 | 1,925 | ZNF276 | PIEZO1 | 950,487 bp | K562 | next 36 |
| TMEM97_K562_chr17_27336k **(unaskable)** | chr17:27,336,191-27,336,690 | 499 | TMEM97 | WSB1 | 982,653 bp | K562 | next 36 |
| TMEM97_K562_chr17_27570k **(unaskable)** | chr17:27,570,114-27,570,613 | 499 | TMEM97 | ENSG00000266728 | 748,730 bp | K562 | next 36 |
| DHRS13_K562_chr17_28864k | chr17:28,864,868-28,866,626 | 1,758 | DHRS13 | ERAL1 | 37,323 bp | K562 | next 36 |
| SOCS3_K562_chr17_78257k | chr17:78,257,711-78,258,455 | 744 | SOCS3 | TMEM235 | 101,995 bp | K562 | next 36 |
| VAPA_K562_chr18_9802k | chr18:9,802,941-9,805,250 | 2,309 | VAPA | TXNDC2 | 109,861 bp | K562 | next 36 |
| VAPA_K562_chr18_9889k | chr18:9,889,157-9,889,957 | 800 | VAPA | TXNDC2 | 24,399 bp | K562 | next 36 |
| VAPA_K562_chr18_9898k | chr18:9,898,666-9,899,734 | 1,068 | VAPA | TXNDC2 | 14,756 bp | K562 | next 36 |
| MRPL4_K562_chr19_10934k **(unaskable)** | chr19:10,934,951-10,936,375 | 1,424 | MRPL4 | TIMM29 | 683,699 bp | K562 | next 36 |
| ADGRE2_K562_chr19_14355k | chr19:14,355,708-14,356,714 | 1,006 | ADGRE2, DDX39A | ADGRE5 | 63,171 bp | K562 | next 36 |
| PDCD5_K562_chr19_32686k | chr19:32,686,989-32,688,162 | 1,173 | PDCD5 | NUDT19 | 106,387 bp | K562 | next 36 |
| EIF3K_K562_chr19_38636k | chr19:38,636,860-38,637,359 | 499 | EIF3K | ACTN4 | 18,035 bp | K562 | next 36 |
| CEBPB_WTC11_chr20_49888k | chr20:49,888,429-49,888,928 | 499 | CEBPB | SPATA2 | 301,903 bp | WTC11 induced pluripotent stem cell | next 36 |
| MSN_K562_chrX_66005k | chrX:66,005,342-66,006,525 | 1,183 | MSN | VSIG4 | 338,289 bp | K562 | next 36 |

**10 of 60 died at the reach filter** - 1 to 4 of the first 24 was the registered prediction, deliberately higher than the third set's registered 0 of 9, because a rule that selects elements skipping over a nearer gene selects for distance. It came back 4 of the first 24 and the same rate holds over all 60. `loci.read_reach` was called, not reimplemented: a gene-**body** test against the scorer's 1,048,576 bp input. They are TMEM54 (965 kb), GMPR (938 kb), MYC (1967 kb), DOLPP1 (652 kb), PCBP2 (819 kb), DENR (940 kb), ZNF276 (949 kb), TMEM97 (982 kb), TMEM97 (748 kb), MRPL4 (675 kb). Beside the other frames: 3 of 17 unaskable (SHH_ZRS, SOX9_PierreRobin, FTO_IRX3's IRX5), 2 of 11 died, both exhibits placed on a bare offset, 0 of 9, which was the registered prediction there.

The element-to-target distances of the 50 graded loci run 1,199 bp to 518,848 bp, median 106,387 bp.

### What it cost

**4 requests for the whole frame, 50 graded loci.** The finished all-chromosome sweep had already deleted an element inside 46 of the 50 graded intervals, and the registered budget of 30 was never approached. Section 6 spent 2,283 requests on twelve loci; once the sweep is finished, a new locus costs only what the intervals nobody annotated cost. The cap raise took n from 20 to 50 for 3 further requests.

**The caveat section 21 recorded applies to these requests too.** `Context.score_region` substitutes an overlapping ENCODE element when there is one, and it did so at PCBP1_WTC11_chr2_70086k, SEPHS2_Jurkat_chr16_30472k. What was deleted is an annotated element *overlapping* the perturbed interval, not the screen's interval itself. It is recorded per locus in `loci_fourth_intervals`.

### The scores, beside the other three frames

The positive control is excluded from every rate below, as registered; it is reported on its own further down.

| field | **the fourth frame** | the third set | the nine candidates | the seventeen panel |
|---|---|---|---|---|
| right target, derived, every locus the frame contains | **10/60 (0.167)** | 8/9 (0.889) | 8/9 (0.889) | 15/17 (0.882) |
| right target, where the model could answer | **10/50 (0.200)** | 8/9 (0.889) | 8/9 (0.889) | 13/15 (0.867) |
| right target, nearest TSS in node | **3/50 (0.060)** | 7/9 (0.778) | 3/9 (0.333) | 11/17 (0.647) |
| right target, annotation lookup | **0/50 (0.000)** | 3/9 (0.333) | 1/9 (0.111) | 8/17 (0.471) |
| chance floor (random coding gene in the window) | **5.48/50 (0.110)** | 2.4/9 | 1.06/9 | 7.12/17 |

**The derived rate is 10/60 (0.167) over every locus the rule returned and 10/50 (0.200) where the model could answer.** The first row is the one comparable with the panel's 15/17, because that rate counts the panel's three unaskable loci as misses; the second is comparable with its 13/15. The three hand-curated frames read 0.882, 0.889 and 0.889 on the first and 0.867, 0.889, 0.889 on the second. This frame is the same benchmark, the same scorers, the same hit rules and the same machine; what changed is that the published target is no longer the nearest gene.

It is **above its own chance floor** - 0.200 against 0.110, so the layers are not guessing - and it is **nowhere near the number this benchmark has been reporting since 2026-09-12**. The registration wrote down in advance what each outcome would mean and this one is the middle case: the number is reported with the per-layer split and without a verdict bolted onto it. The split is where the finding is.

### The two halves, because the cap was applied in genome order

The first 24 were scored on 2026-09-21 under the registered cap; the next 36 were added by raising it, under the same rule, and commit f407915 wrote down before they were scored what they would have to read for the frame to be called stable and what a divergence would mean. The cap ran in genome order, so the first half is chr1 to chr10 and the second is chr11 to chrX - and the later chromosomes carry the gene-dense stretches, where the nearest coding TSS is closer and more coding genes sit in the window. Each half is therefore measured against **its own** chance floor.

| | the first 24 (chr1-chr10) | the next 36 (chr11-chrX) | all 60 |
|---|---|---|---|
| drawn by the rule | 24 | 36 | 60 |
| died at the reach filter | 4 | 6 | 10 |
| graded | 20 | 30 | 50 |
| right target, derived, every drawn locus | 6/24 (0.250) | 4/36 (0.111) | **10/60 (0.167)** |
| right target, where the model could answer | 6/20 (0.300) | 4/30 (0.133) | **10/50 (0.200)** |
| right target, nearest TSS in node | 0/20 (0.000) | 3/30 (0.100) | 3/50 (0.060) |
| chance floor, computed per half | 2.73/20 (0.137) | 2.75/30 (0.092) | 5.48/50 (0.110) |
| deletion layer, the model itself | 1/20 | 2/30 | 3/50 |
| deletion layer named the nearest coding TSS instead | 14/20 | 12/30 | 26/50 |

**The second half reads 4/30 (0.133), BELOW the registered interval and below the first half's 6/20 (0.300) (-0.167).** That is the direction the gene-density gradient predicts, so it is reported against the second half's own floor rather than as a new finding; it does not soften the combined rate and it is not smoothed into it. Against their own floors, which is how the registration said to read them: the first half is 0.300 on a floor of 0.137 and the second is 0.133 on a floor of 0.092, so the margin over chance is 2.2x and 1.4x. The floor itself fell from 0.137 to 0.092 between the halves, which is the gene-density gradient showing up in the arithmetic rather than in the readings: more coding genes in the window makes a random draw likelier to be wrong, not likelier to be right. So the second half is weaker than the first even after its easier floor is taken into account, and the density gradient does not account for the whole of the drop. The honest statement is that the frame's rate is the combined one and that the half that was scored first was the better of the two. Commit f407915 registered, before the second half was read, that a second-half rate between about 0.15 and 0.45 would mean the first half's 0.300 was not an n = 20 accident, that a second half BELOW the first is the expected direction of a gene-density effect and not a new finding, and that a second half materially ABOVE the first would be the interesting divergence because it runs against that gradient. What is not claimed is that the halves are independent replications: they are one rule applied to one file and split by an arbitrary cap, so a difference between them is a statement about the genome's arrangement and about n.

### Which layer, and what it named instead - this is the result

| layer | provenance | names the published target first |
|---|---|---|
| deletion | derived | 3/50 |
| eqtl | derived | 3/50 |
| gene_input | derived | 4/50 |
| node | heuristic | 3/50 |
| lookups | looked_up | 0/50 |

**The deletion layer - the model - names the published target at 3/50.** The derived rate above is a union over three derived layers, and the two that carry it are GTEx eQTL (3/50) and the summed-window reading (4/50), neither of which is the model's own answer about the element.

And what it names instead is the whole point of drawing the frame this way:

| what it named | the deletion layer | the nearest-TSS-in-node rule |
|---|---|---|
| the published target | 3/50 | 3/50 |
| the nearest coding TSS | 26/50 | 41/50 |
| another gene | 19/50 | 2/50 |
| nothing | 2/50 | 4/50 |

**The model names the nearest coding TSS at 26 of 50 loci and the published target at 3.** At these loci the two are different genes by construction, and when they differ the model goes with the nearer one 26 times and with a third gene 19 times - so it is not only proximity, but the published answer is the rarest of the three. The node heuristic names the nearer gene at 41 of 50, which is what a pure proximity rule looks like, and the deletion layer sits between that and the publication rather than at the publication.

**That is the sharpest reading this benchmark has produced, and it is a negative.** It does not say the deletion layer is the nearest-gene rule in general - at the three curated frames it also got directions and cell types right, which proximity alone cannot do. It says that on the axis this frame isolates, where proximity and the published answer disagree, the model names the published target 3 times in 50 and the nearer gene 8 times as often; and that the three frames' agreement at 0.88 was measured almost entirely where the two agree, so it could not have shown this.

### The two nearest-gene rules, kept apart, because one of them is arithmetic

| rule | reads | status |
|---|---|---|
| nearest coding TSS anywhere | 0/50 | **0 by construction** - it is the rule that drew the frame, and it is printed so the construction is visible |
| nearest coding TSS in the CTCF node | 3/50 (0.060) | a genuine reading: a different rule, right wherever the node excludes the nearer gene |

The node rule was registered in advance as low but **not** zero by construction, and it came back 3/50: at 41 loci it named the nearer gene and at 4 it named nothing. So on this frame the node is a proximity rule with abstentions, and what it adds over raw proximity is the handful of loci where the node boundary excludes the nearer gene. The comparison the doc has been making - 0.647, 0.333, 0.778 - is with that rule, and it is 3/50 (0.060) at the one frame where the rule was not allowed to be right by default. **On the where-the-model-could-answer line its spread across four frames is now 0.718 and the derived rate's is 0.689, against 0.022 across the three curated frames alone.** Both numbers move once the geometry is fixed; it is no longer only the baseline that depends on who picked the loci.

### The positive control, and why the run is interpretable at all

**The control passed.** CONTROL_HMGA1_WTC11_chr6_34235k was drawn by the same rule from the same file - among the elements this frame REJECTED for having the shortcut right, the one with the smallest published element-to-TSS distance (1,111 bp), which is the easiest case the data can offer. Deleting it named HMGA1 first at -0.2271, by the deletion layer (deletion). Registered in advance: a failure here would mean the reading is broken and nothing else in the run is interpretable. It is not broken - the same layer that answers correctly at 1.1 kb names the wrong gene at 47 of the 50 loci where the right answer is not the nearest one.

### The matched windows, with presence counted before coverage is made a stratum

239 matched windows, the same four covariates and the same hit rules, with the keep-outs extended to every locus in **all four** frames. `genomeos/compare.py` does the standardising; nothing here reimplements it.

| input | kind | targets | controls |
|---|---|---|---|
| a deletion spent inside the window | **bought** | 50/50 | 97/239 |
| a GTEx eQTL distilled for the window | **bought** | 48/50 | 198/239 |
| the mammalian constraint track over the window | **free** | 50/50 | 239/239 |

| claim | the fourth frame | controls, covariates only | p | controls, **coverage held fixed** | p |
|---|---|---|---|---|---|
| target | 1.000 | 0.886 (+0.114) | 0.006 | **0.983** (+0.018) | 0.188 |
| cell | 0.980 | 0.473 (+0.507) | 0.000 | **0.606** (+0.372) | 0.000 |
| direction | 0.959 | 0.363 (+0.597) | 0.000 | **0.802** (+0.175) | 0.003 |
| storage | 0.061 | 0.049 (+0.012) | 0.395 | **0.082** (-0.037) | 0.760 |

**The coverage artefact reproduces on a fourth independent frame.** The target claim - *some* coding gene is named - separates raw and does not survive the coverage stratum, which is what sections 8, 19 and 21 each found. These claims are about whether a layer says anything at all, not about whether it is right, and the rates above are what say whether it is right.

### What this frame cannot say, stated rather than buried

- **A CRISPRi positive is not always a direct target.** The benchmark file carries a `direct_vs_indirect_negative` column estimating exactly that, and this frame did not use it: the rule takes the file's own `Regulated` flag. Some of these 50 published targets may be indirect consequences of silencing the element, in which case naming the nearer gene is not as wrong as the rate makes it look. Using that column would be a different frame and it would have to be registered as one.
- **3 of the drawn loci have an unnamed novel gene as their nearest coding TSS** (ENSG00000288678, ENSG00000293615, ENSG00000266728), so the shortcut is wrong there because of an annotation rather than because of a genuine gene skip. 1 of those is in the graded 50 (PCBP1_WTC11_chr2_70086k).
- **The frame is one cell line more than it looks.** The drawn cells are HCT116 colorectal carcinoma, Jurkat T-lymphoblast, K562, WTC11 induced pluripotent stem cell, and K562 dominates because the held-out arm does. A cell-type imbalance is a property of the field's screens, not of this rule.
- **The frame is now every element the rule passes, so it is not a sample of anything.** The genome-order cap that made the first run chromosome 1 to 10 has been raised past the last locus; what remains is that the held-out arm is itself whatever the field happened to screen.
- **97 of 239 control windows have a deletion spent inside them against 50/50 of the loci**, so the coverage stratum is a real test here rather than arithmetic.

### Left undone

- **The frame is exhausted.** Every element the rule passes is drawn, so there are no more loci to add without a different rule, a different file, or the K562 training arm this registration excludes. Taking n past 60 means registering a new frame, not raising a cap.
- **The non-coding frame is already drawn and unscored.** The 14 elements dropped at step 4 are the first set this project has ever had of published, perturbation-backed, non-coding targets, which is what section 21 said curation could not produce. They would test the `predicted_coding`-first defect H19 exposed.
- **The direction axis was registered unaskable here and is.** Every CRISPR positive is an element whose silencing lowers its target, so the repressor arm is empty by construction and section 21's sign question stays open at n = 2.
- **`loci.STATED_INTERVAL_RESULTS` now lists `loci_fourth_intervals`** (f407915), so this frame's paid deletions are read back from the file rather than from a global patched at import time. `loci_fourth.register_stated_intervals` is kept because it is idempotent and because a test asserts the two say the same thing.
- **Whether the earlier frames' rates should now be reported with their element-to-target distances beside them** is the benchmark owner's call. This frame's reading is that a derived rate quoted without that distance is not comparable between frames, and sections 9, 19 and 21 do not carry it.
- **The four frames are still four frames.** Nothing here pools them, and the fourth one least of all: it was drawn to disagree with the other three and it does.

---
