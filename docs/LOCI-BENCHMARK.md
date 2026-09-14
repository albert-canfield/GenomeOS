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
where 47 of 60 windows have no element-level deletion. The controls sit on the
panel's own chromosomes, so the sweep's remaining order matters to the benchmark
only where it reaches chr11 (5 controls), chr9, chr8, chr7 and chr6 (5 each) and
chr2 (15). Each control window the sweep reaches has about a 5-in-6 chance of
producing a direction and about a 1-in-7 chance of producing a value on syntax.

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