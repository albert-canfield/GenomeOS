# Simpler genomes as the comparison: C. elegans first

Area I, step 6 of the 2026-09-11 plan ("function for the non-coding 98% via attribution, closure
tests and simpler genomes"). Lane `lane-simplergenome`, started 2026-10-10. No money was spent and
no model request was sent; AlphaGenome is not used anywhere in this document.

Code: `genomeos/attribution/simpler.py`, `scripts/simpler_worm_run.py`, `tests/test_simpler_worm.py`.
Results: `data/results/simpler_worm_inventory.json`, `simpler_worm_registration.json`,
`simpler_worm_result.json`.

## Why a smaller genome

Every registered test of the human compiled attribution so far found nothing beyond its baselines,
and its one positive signal is lost at the bridge (ROADMAP §5). In the human genome a measured
element is the exception: 4,554,153 of the unknown space's 1,008,756,271 bases (0.45%) overlap a
lentiMPRA, VISTA or CRISPRi-benchmark element (docs/ATTRIBUTION.md). The worm has about 100 Mb, a gene every 5 kb, a complete lineage the project
already runs, and one public table in which the same five tissues were measured for both chromatin
accessibility and RNA. There the attribution's logic can be held against measured elements in the
thousands rather than the hundreds.

## 1. What exists (inventory, no outcome read)

Checked on 2026-10-10. Sizes are bytes as served where they were checked; "not checked" where
they were not. Full list with URLs: `data/results/simpler_worm_inventory.json`.

| source | what it measures | assembly | size | licence | fetched |
|---|---|---|---|---|---|
| Serizay et al. 2020, Genome Res 30:1752, Supplemental Table S2 | ATAC accessibility (RPM) of 47,514 elements in sorted nuclei of five young-adult tissues (germline, neurons, muscle, hypodermis, intestine) with a tissue call each; element classes (20,271 putative enhancers, 15,277 coding promoters, ...); nuclear RNA TPM of 20,222 protein-coding genes in the same five tissues | ce11 | 28,664,943 B bundle streamed, 14,777,075 B kept | CC BY 4.0 | yes |
| GEO GSE141213 (the same study's tracks) | per-tissue ATAC and stranded RNA bigWigs | ce11 | 2,386,268,160 B tar; 7 to 237 MB per bigWig | public GEO record | no: Table S2 holds the distilled values |
| Jänes et al. 2018, eLife 7:e37344 | ATAC and nuclear RNA across six stages and ageing, 42,245 typed elements | ce11 | not checked | CC BY 4.0 | no: Table S2 extends the same pipeline |
| UCSC ce11 phyloP135way; phyloP26way, phastCons26way | per-base conservation across 135 or 26 nematodes | ce11 | 356,259,120 B; 180,170,767 B; 159,951,981 B (bigWig) | UCSC download, terms for these tracks not checked | no: readable by range request (`attribution/bigwig.py`) |
| Durham et al. 2021, Genome Res 31:1952 | single-cell accessibility at L2, by cell type | ce11 | not checked | not checked | no |
| Daugherty et al. 2017, Genome Res 27:2096 | ATAC in embryo, L3, young adult, with reporter tests of selected elements | ce11 | not checked | not checked | no |
| WormBase regulatory features (GFF3) | curated enhancer, promoter, regulatory_region, TF_binding_site records | WBcel235 | not checked | not checked | no |
| Genome-wide reporter assays (STARR-seq, MPRA) | none found in one web search | | | | |
| Kaletsky 2018, CeNGEN 2021, modERN 2018 | further tissue and neuron-class RNA; TF binding | ce11 | not checked | not checked | no |

What the project already held for the worm: the WBcel235 genome and Ensembl gene models
(`data/reference/celegans/`, 30.3 MB and 8.2 MB, git-ignored; 19,985 protein-coding genes), the
Ma et al. 2021 TF reporter atlas (88 MB, CC BY 4.0), sixteen `data/results/celegans_*.json`
results (lineage, Packer 2019 cell types, contacts, fates, TF atlas) and the `.bio` organism files
held by lane-wormcontact. None of it is regulatory-element, accessibility, conservation or tissue
RNA data: before this lane the project had no measurement of the worm's non-coding genome.

Downloads for this lane: 28,664,943 bytes, against a budget of 1.5 GB. The disk had 23 GiB free.

## 2. What transfers without AlphaGenome

The human pipeline, piece by piece, and what a worm version needs.

| human piece | transfers? | worm form | measured worm data that can test it |
|---|---|---|---|
| Element-to-gene by deletion in AlphaGenome (`enhancer_targets_all`) | no | none: the project sends no model requests, and AlphaGenome is human and mouse only | — |
| Node model (CTCF domains; an element acts inside its node) | not as it is: C. elegans has no CTCF, and domain structure is clear mainly on the X (Crane et al. 2015) | positional gene domain: an element belongs to the gene it sits in, else to the gene it sits upstream of (the registered rule below) | Serizay 2020 tissue accessibility against tissue RNA: **registered below** |
| Nearest gene / 1/distance (the human baselines) | yes, unchanged | nearest annotated 5' end | the same table; it is the baseline |
| Tiers from constraint (`budget.py`) | yes in form; the threshold does not | nematode phyloP (135-way or 26-way), read by range request; the cut (2.27 at 5% FDR in Zoonomia) must be re-derived for nematodes, and the 135-way tree's branch lengths are k-mer derived (UCSC README) | whether a constrained tier is enriched for elements accessible in some tissue (accessibility is not function); Daugherty 2017's reporter-tested elements for a handful |
| The reader (which cells open an element) | yes, directly | the five-tissue ATAC calls; Durham 2021 by L2 cell type | the same table's RNA |
| Closure at the gene level (element, gene, RNA) | yes, in the form registered below | element tissue against attributed gene's tissue RNA | Serizay 2020 |
| Closure at the cell and organism level | later | the worm's lineage is complete and already in the project: element, gene, cell fate | needs per-cell element activity (Durham 2021) and per-cell RNA (Packer 2019) joined on the lineage |
| Compiled rules that run (`runrules`) | not now | would need worm TF binding (modERN) and rates | — |

Conclusion of this section: without a model, the human attribution reduces to its node model and its
distance baselines. The node does not exist in the worm in its human form, so the first worm test is
whether a gene-structure definition of the node does better than distance alone, on an endpoint the
attribution never reads.

## 3. The registered test

Registered in `genomeos/attribution/simpler.py` (`REGISTRATION`) and
`data/results/simpler_worm_registration.json`, committed before the RNA sheet was read.
`scripts/simpler_worm_run.py score` refuses to run unless both are committed and unchanged.

- **Attribution (the domain rule).** Host protein-coding gene if the element's midpoint lies in one
  (nearest 5' end among overlapping hosts); otherwise the flanking gene whose 5' end faces the
  element, i.e. the gene it sits upstream of; a divergent pair takes the nearer 5' end; a
  convergent pair (between two 3' ends) takes the nearer 5' end of the two and is counted apart.
  It reads WBcel235 positions and strands only: no accessibility, no RNA, and none of the source's
  own element-to-gene columns.
- **Baseline.** The gene with the nearest annotated 5' end (nearest TSS).
- **Population.** The source's putative enhancers accessible in exactly one of the five tissues by
  the source's own call: 5,814 elements (germline 470, neurons 1,434, muscle 1,432, hypodermis
  1,487, intestine 991). The two rules name different genes on 1,248 of them (21.5%): 942 inside
  gene bodies, 300 between genes where another gene's 5' end is nearer than that of the gene the
  element sits upstream of, 6 convergent. 103 1-Mb clusters.
- **Endpoint.** Concordant when the attributed gene's TPM in the element's tissue is at least 5 (the
  source's detection threshold) and higher than in each of the other four tissues.
- **Statistic.** Concordance of the domain rule minus that of nearest TSS on the same elements, in
  points, with a 95% cluster-bootstrap interval over 1-Mb windows.
- **Matched control and falsifier.** Per element, a gene drawn from the same chromosome at the same
  5'-end distance bin as the domain-rule gene; the rule is falsified as a carrier of tissue
  information if its concordance is not above that control (one-sided p ≥ 0.05).
- **Positive control.** Single-tissue coding promoters, nearest-TSS gene, against the same kind of
  control; if it does not separate, the readout is declared unresponsive and the reading withheld.
- **Power in advance.** With 1,248 differing elements of 5,814, the smallest difference with 80%
  power at two-sided 5% is about 1.7 points before any inflation for clustering.
- **Stop rule.** One run; fewer than 200 differing elements would have stopped it unscored. Nothing
  is changed after the RNA sheet is read.
- **What it cannot establish.** Accessibility is not function: concordance says the attributed gene
  is expressed where the element is open, not that the element acts on it. An element inside a gene
  body can be open where its host is transcribed for that reason alone, which favours the domain
  rule by construction of the assay; the intergenic elements are reported apart, and the reading
  carries a registered caveat if the difference is not carried by them. The annotated 5' end is not
  always the start of transcription: most worm transcripts are trans-spliced.

## 4. Result (2026-10-10, registration `595d7c2`)

Run once by `scripts/simpler_worm_run.py score` on the committed registration, 11 seconds, no
request to any model. Outcome **DOMAIN_RULE_AHEAD**; falsifier **NOT_FALSIFIED**. The registered
reading, verbatim:

> On 5,814 C. elegans putative enhancers accessible in exactly one of five tissues (Serizay et al. 2020, young adults), the gene the domain rule attributes is expressed at 5 TPM or more and highest in the element's own tissue for 3,174 of 5,814 (54.6%), against 3,006 of 5,814 (51.7%) for the gene with the nearest 5' end; the rules differ on 1,248 elements; difference +2.89 points, 95% +1.83 to +3.91 over 1-Mb clusters. On this endpoint the domain rule is ahead of nearest TSS. The domain rule's genes match the element's tissue more often than genes at a matched distance (24.1%, one-sided p 0.001), so position carries tissue information here. Accessibility is not function: this says the attributed gene is expressed where the element is open, not that the element acts on it.

The positive control separates, so the readout is responsive: of 5,347 single-tissue coding
promoters, the nearest-TSS gene is concordant for 3,749 (70.1%) against
31.6% for genes at a matched distance (one-sided p 0.001, the floor of
1,000 draws). 103 of the 11,628 attributed gene slots (two per element) name a gene absent
from the RNA sheet; the registration counts them as not concordant.

Beside the primary, deciding nothing (registered strata; percentile intervals over 1-Mb clusters):

| stratum | elements | domain rule | nearest TSS | difference, points [95%] |
|---|---|---|---|---|
| inside a gene body | 3,479 | 55.0% | 52.7% | +2.33 [+0.71, +3.88] |
| between genes | 2,335 | 53.9% | 50.2% | +3.73 [+2.79, +4.74] |
| the 1,248 where the rules differ | 1,248 | 43.9% | 30.4% | +13.46 [+8.74, +18.02] |
| germline | 470 | 60.2% | 59.8% | +0.43 [-3.07, +3.79] |
| neurons | 1,434 | 50.8% | 46.2% | +4.60 [+2.54, +6.85] |
| muscle | 1,432 | 45.5% | 41.3% | +4.12 [+2.11, +6.02] |
| hypodermis | 1,487 | 56.6% | 56.5% | +0.13 [-2.21, +2.17] |
| intestine | 991 | 67.6% | 63.7% | +3.94 [+0.89, +7.18] |

The nearest-TSS gene is also above its own matched-distance control (26.8% expected,
one-sided p 0.001). The difference is not confined to gene bodies: between genes it is +3.73
points, so the registered host caveat was not triggered. That does not remove the confound the
caveat names. Most worm transcripts are trans-spliced, and the outron upstream of the annotated 5'
end is transcribed with the gene, so an element between genes that the rule attributes to the gene
it sits upstream of can also lie in that gene's transcribed sequence. This run does not separate
regulation from transcription in either stratum. Two of the five tissues (germline, hypodermis)
show no difference on their own.

What this does and does not say. On this endpoint the domain rule is ahead of nearest TSS: where
the human node cannot be carried over, a node defined by gene structure names a gene expressed
highest in the element's tissue for 54.6% of 5,814 elements against 51.7% for nearest TSS, on an
endpoint the rule never read. It does not say that any element acts on the gene it is attributed
to: accessibility is not function, and concordance of tissue is not a measured effect of the
element. It is one stage (young adults), five tissues, one source.

## 5. What follows

1. Replication at a second stage or resolution, registered before it is run: Jänes 2018 (six
   stages) or Durham 2021 (L2 cell types), same rule, same baseline, same control.
2. The worm constraint tiers: nematode phyloP by range request (`attribution/bigwig.py`, no
   download), with a threshold derived for nematodes rather than borrowed from Zoonomia; tested
   against which elements are accessible in any tissue, which is not function either.
3. Closure at the cell level on the worm's complete lineage: element, gene, cell fate, once
   per-cell element activity and per-cell RNA can be joined on the lineage.
4. Then Drosophila (genome-wide STARR-seq and embryo enhancer tiles exist; sizes and licences not
   yet checked), fugu as the compact vertebrate, and mouse for transfer, each without a model request.

## 6. Stage 2: the worm's constraint tiers, registered (2026-10-10, lane-simpler2)

Item 2 of section 5. No money, no key, no model request. The phyloP tracks were read by HTTP range
request and nothing of them was stored.

Code: `genomeos/attribution/simpler_tiers.py` (imports `simpler.py`, unchanged),
`scripts/simpler_tiers_run.py` (`tiers`, `register`, `score`), `tests/test_simpler_tiers.py`.
Results: `data/results/simpler_worm_tiers.json` (cuts, histograms, per-element counts; `1e61925`),
`simpler_worm_tiers_registration.json` (`752f7fc`), `simpler_worm_tiers_result.json` (`25e25fd`).
Each was written from a clean worktree of the committed tree.

### How the human tiers were derived, and what transfers

A human base is constrained when its Zoonomia phyloP score is positive and its q value is at most
0.05. Christmas et al. 2023 (Science 380:eabn3943), Methods: "We computed a mammalian phyloP
threshold by converting the p values corresponding to the phyloP scores into q values using a FDR
correction. We considered any column with a resulting q ≤ 0.05 to be significantly evolutionarily
constrained or accelerated, as determined by the sign of the score." On the 241-mammal track that
procedure gives 2.27 (`constraint.PHYLOP_THRESHOLD`). A block of the regulatory class is "under
constraint" when at least 5% of its measured bases pass (`budget.CONSTRAINED_MIN`, `budget._rule`).

The nematode tracks were checked on 2026-10-10. The UCSC ce11 READMEs and track descriptions give
no significance cut-off. They say "the absolute values of the scores represent -log p-values under
a null hypothesis of neutral evolution", computed per column by a likelihood-ratio test
(`--method LRT`, `--mode CONACC`), with the neutral model fitted by phyloFit on 4-fold degenerate
sites. Neither track has a publication of its own. phyloP135way aligns 134 genomes to C. elegans:
111 nematodes, 22 flatworms and Ciona. Its guide tree comes from shared k-mers ("not necessarily
biologically correct", README). phyloP26way aligns 26 nematodes. Neither the number 2.27 nor any
other number was borrowed. Zoonomia's procedure was applied to each track's own scores:
p = 10^-|score| for every scored base, Benjamini-Hochberg at q = 0.05 over both signs, and a
constrained base is one with a positive score at or above the cut.

The distribution was read from a uniform sample of data blocks. The blocks of chrI to V and X, in
file order, were cut into consecutive groups of 16, and 10% of the groups were drawn with seed
20261010. Every scored base sits in exactly one block, so every scored base had the same chance of
being sampled. The cut's interval is a bootstrap over the sample's 1-Mb windows.

| track | 5% FDR cut [95%] | scored bases sampled | pass: all scored | non-exon | CDS | other exon | intron | intergenic | bytes moved |
|---|---|---|---|---|---|---|---|---|---|
| phyloP135way (primary) | 2.08 [2.05, 2.12] | 9,496,101 | 15.7% | 7.2% | 39.6% | 13.7% | 7.7% | 6.6% | 56,378,275 |
| phyloP26way (beside) | 2.15 [2.12, 2.17] | 5,848,517 | 14.4% | 3.1% | 31.3% | 8.4% | 3.9% | 1.7% | 32,712,304 |

The registered instrument check passed on both tracks: the cut places more CDS than intron bases
above it. For comparison, Zoonomia's 5% FDR set is 3.26% of the human genome. Coding sequence
dominates what aligns across nematodes, and that accounts for most of the worm's larger share.

**Element tiers.** Stage 1's 5,814 single-tissue putative enhancers are all 150 bp. Exon bases of
any transcript were removed first, because the human tiers are measured on unknown space, which
excludes exons. The human regulatory-class rule then applies: an element is constrained when at
least 5% of its remaining scored bases pass. On phyloP135way, 3,894 elements are constrained, 1,117
are not, and 803 have no scored base left (not measured). On phyloP26way the counts are 1,290,
2,860 and 1,664. The human 5% rule puts 78% of the measured elements (3,894 of 5,011) in the
constrained tier on the primary track. Under that rule 8 of 150 bases are enough, while 7.2% of
the worm's non-exon scored bases pass the cut. The human rule was set in a genome where 3.26% of
all bases pass. This was seen in the tier counts before anything was registered. The rule was kept
as the primary and a dose view was named beside it.

### The registered test

Registered in `simpler_tiers.REGISTRATION` and committed in `752f7fc` before the RNA sheet was read
against any tier. `score` refuses to run unless the registration, the tiers file and the module are
committed and unchanged, and unless the census recomputes to the registered one.

- **Question.** Does constraint carry tissue information beyond the domain rule? Among the
  single-tissue elements, is the domain-rule gene concordant (stage 1's endpoint, unchanged) more
  often for constrained elements than for unconstrained ones?
- **Matching.** Strata are tissue × inside a gene body or between genes × the distance bin of the
  domain-rule gene's 5' end (stage 1's 8 bins) × measured-length tertile. The length tertile edges
  coincide at 150 bp, so the length strata are in effect "fully measured" and "less". 100 of 113
  strata hold both tiers: 3,865 constrained and 1,109 unconstrained elements.
- **Statistic.** The difference in concordance, standardised to the constrained elements' strata.
  The 95% interval is a percentile cluster bootstrap over 1-Mb windows. The MDE was 4.8 points
  before clustering.
- **Falsifier.** A within-stratum permutation of the tier label (1,000 permutations). The claim is
  FALSIFIED if p ≥ 0.05.
- **Decision.** CONSTRAINT_ADDS, NO_DIFFERENCE_DETECTED or CONSTRAINT_BELOW, read from the interval.
- **Stop rule.** Fewer than 200 elements per tier in support stops the run unscored. The test is run
  once.
- **Confound, stated in advance.** Stage 1 named transcription. Transcribed sequence can be
  constrained for splicing, trans-splicing or RNA-level reasons. Exons are removed and position is
  a stratum, but introns and outrons are not removed, and outrons are not annotated. A difference
  would therefore not separate regulatory constraint from transcription-linked constraint.

### Result (2026-10-10, registration `752f7fc`)

Outcome **NO_DIFFERENCE_DETECTED**. Falsifier **FALSIFIED** (permutation p 0.105). The registered
reading, verbatim:

> On 5,011 C. elegans putative enhancers accessible in exactly one of five tissues (Serizay et al. 2020) with constraint measured on their exon-free bases (UCSC ce11 phyloP135way, cut re-derived at 5% FDR by Zoonomia's procedure), the domain-rule gene is expressed at 5 TPM or more and highest in the element's own tissue for 2,093 of 3,865 constrained elements (54.2%), against 51.5% for the 1,109 not constrained at the same tissue, position class, distance and measured length; difference +2.63 points, 95% -1.35 to +7.03 over 1-Mb clusters (within-stratum permutation p 0.105; FALSIFIED). No difference was detected: on this endpoint constraint adds no detectable tissue information to the domain rule's position. This is not evidence that constrained elements lack function. Conservation and accessibility are not function, and this run does not separate regulatory constraint from transcription-linked constraint.

The views below are reported beside the registered test and decide nothing. Each is matched the
same way, with a percentile interval over 1-Mb clusters and a within-stratum permutation p.

| view | constrained / not, in support | constrained | matched not | difference, points [95%] | p |
|---|---|---|---|---|---|
| inside a gene body | 2,263 / 505 | 54.2% | 51.9% | +2.32 [-3.52, +8.13] | 0.18 |
| between genes | 1,602 / 604 | 54.1% | 51.0% | +3.08 [-2.08, +8.92] | 0.13 |
| germline | 163 / 171 | 56.4% | 60.3% | -3.89 [-13.34, +6.34] | 0.75 |
| neurons | 1,096 / 225 | 52.4% | 46.3% | +6.06 [-1.47, +14.44] | 0.061 |
| muscle | 972 / 221 | 46.4% | 34.1% | +12.31 [+3.49, +21.02] | 0.004 |
| hypodermis | 1,028 / 276 | 56.1% | 59.7% | -3.54 [-11.86, +5.16] | 0.81 |
| intestine | 606 / 216 | 65.8% | 72.7% | -6.84 [-14.76, +2.53] | 0.97 |
| phyloP26way tiers (cut 2.15) | 1,287 / 2,819 | 57.1% | 52.4% | +4.71 [+0.95, +8.36] | 0.007 |
| dose view: constrained fraction ≥ 0.19 against < 0.08 (phyloP135way) | 1,664 / 1,563 | 54.6% | 47.9% | +6.69 [+2.11, +11.68] | 0.001 |

Without matching, 2,104 of 3,894 constrained elements (54.0%) are concordant, against 599 of 1,117
unconstrained (53.6%). 28 domain-rule genes are absent from the RNA sheet and count as not
concordant, as in stage 1.

**What this says and what it does not.** On the registered test, the worm's constraint tier carries
no detectable tissue information beyond the domain rule's position. The difference is +2.63 points
and its interval spans zero, and the permutation does not separate it from chance. This is a
negative result for the registered tier: the 5% rule taken from human, on the 135-genome track.

It is not evidence that constrained elements lack function. Two views beside the test point the
other way, and neither decides anything:

- **The dose view.** The top third of elements by constrained fraction against the bottom third
  gives +6.69 [+2.11, +11.68].
- **The 26-genome track.** Its tiers on the same elements give +4.71 [+0.95, +8.36]. That track
  has a stricter tier, with 31% of measured elements constrained against 78% on the primary.

Both are consistent with a reading that a coarse tier dilutes a graded signal. That reading was not
the registered question. Testing it needs its own registration, on a second source, before it is
run. The tissues disagree: muscle +12.31 [+3.49, +21.02], and germline, hypodermis and intestine
negative, each interval spanning zero. The confound named in advance stands in every view.
Conservation and accessibility are not function, and none of these views separates regulatory
constraint from transcription-linked constraint.

### What follows from stage 2

1. A registered dose test: constrained fraction as a continuous or top-tertile variable, with the
   cut-points fixed before the run, on a second source (Jänes et al. 2018 or Durham et al. 2021).
   This one cannot be re-tested on Serizay 2020 without the views above informing it.
2. An element rule derived for the worm rather than borrowed. The 5% fraction would be re-derived
   against the worm's own non-exon base rate in the way the base cut was, on a population that
   does not read this endpoint.

## 7. Drosophila: what measured data exists (inventory, no test; 2026-10-10, lane-simpler2)

Checked on 2026-10-10. Sizes are bytes as served (HTTP Content-Length or Content-Range). Only files
under 1 MB were opened. Nothing was kept in the project. The Vienna tiles, phyloP124way, phyloP27way
and Calderon peak URLs were re-checked by HEAD request. "Not stated" means the source states no
licence that was found. No licence is inferred.

| source | measures | assembly | where | bytes | licence |
|---|---|---|---|---|---|
| Kvon et al. 2014, Nature 512:91 (Vienna Tiles) | Gal4 reporter tiles in embryos, stages 4–16: 7,793 tiles, 3,604 active, annotated with anatomy terms per stage | dm3 | enhancers.starklab.org/download/tiles/all | 957,907 | not stated (cite the paper) |
| BDGP in situ (gene side for the tiles) | embryo in situ terms by stage, 99,319 rows; presence/absence, not quantitative | gene IDs | insitu.fruitfly.org/insitu-mysql-dump/insitu_annot.csv.gz | 486,763 | "All rights reserved" |
| Arnold et al. 2013, Science 339:1074 (STARR-seq) | STARR-seq peaks (S2 3,966; OSC 3,835), with DHS and RNA in the same two cell lines | dm3 | GEO GSE40739 | 1,730,560 (archive) | GEO: no restrictions |
| Zabidi et al. 2015, Nature 518:556 | STARR-seq with developmental and housekeeping core promoters, S2 and OSC, 5,408 to 6,217 peaks per set | dm3 | GEO GSE57876 | 64,877 to 88,054 per peak file | GEO: no restrictions |
| REDfly v10 (Rivera et al. 2019, NAR 47:D828) | curated reporter-tested CRMs: 58,241 (17,272 in vivo) | Release 6 annotation; coordinate build not checked | redfly.ccr.buffalo.edu; bulk export not verified | not checked | CC BY-NC-ND 4.0 |
| Cusanovich et al. 2018, Nature 555:538 | sci-ATAC on embryos at 2–4, 6–8 and 10–12 h; no RNA | dm3 | GEO GSE101581 | 894,832,640 (archive) | not stated |
| Calderon et al. 2022, Science 377:eabn5800 | sci-ATAC (about 1M nuclei) and sci-RNA (about 0.5M) over 0–20 h embryos, 110,185 peaks, tissue labels in both | dm6 | shendure-web.gs.washington.edu (DEAP) peaks.bed; GEO GSE190147 | 2,532,609 (peaks); 697,773,156 (RNA matrix) | not stated |
| Reddington et al. 2020, Dev Cell 55:648 | DNase on sorted neural and mesoderm nuclei, 4–12 h | dm3 | ArrayExpress E-MTAB-8881 | 328,669,420 per bigWig | not stated |
| Bozek et al. 2019, Genome Res 29:771 | ATAC in anterior–posterior domains of the blastoderm | dm3 | GEO GSE118240 | 744,281 (one peak file) | not stated |
| FlyAtlas2 (Leader et al. 2018, NAR 46:D809) | tissue RNA-seq, adult and larva; no embryo | gene IDs | FlyBase precomputed high-throughput expression table | 25,841,525 | CC BY 4.0 (FlyBase) |
| UCSC dm6 phyloP27way, phyloP124way | per-base conservation | dm6 | hgdownload.soe.ucsc.edu/goldenPath/dm6/ | 415,668,365; 470,242,283 | "freely available for public use"; no cut-off given |
| modERN (Kudron et al. 2018, Genetics 208:937) | 640 TF ChIP-seq experiments | dm6 | ENCODE portal | about 200 KB per peak file | ENCODE: no restrictions |

**What the fly offers and what it lacks.** Calderon et al. 2022 is the one dm6 source that labels
accessibility and RNA in the same embryo tissues over the same windows, which is the worm test's
design. It measures accessibility rather than tested activity, the two assays come from different
nuclei, and no licence is stated. Measured enhancer activity in vivo (the Vienna tiles) pairs only
with presence/absence in situ calls, so "expressed highest" cannot be measured on it. Only 79 of
228 tile terms match BDGP term names exactly, and BDGP reserves all rights. No quantitative embryo
tissue RNA was found, and no adult or larval tissue accessibility was found to pair with
FlyAtlas2. The tiles, both STARR-seq sets, Cusanovich, Reddington and Bozek are on dm3 and would
need lifting to dm6 for the conservation tracks. UCSC gives no phyloP cut-off for dm6 either, so
the cut would be re-derived here in the same way.

## 8. Drosophila: the worm's domain-rule test on Calderon 2022 (2026-10-10, lane-fly)

The worm test of sections 3 and 4, repeated on a second species and a second source. No money, no
key, no model request.

Code: `genomeos/attribution/simpler_fly.py` (imports `simpler.py`, unchanged),
`scripts/simpler_fly_fetch.py`, `scripts/simpler_fly_run.py` (`inventory`, `register`, `score`),
`tests/test_simpler_fly.py`. Results: `data/results/simpler_fly_inventory.json`,
`simpler_fly_registration.json`, `simpler_fly_result.json`, each written from a clean worktree of the
committed tree.

### The data (fetched 2026-10-10)

Calderon et al. 2022 (Science 377:eabn5800) profiled about 1M nuclei by sci-ATAC and about 0.5M cells
by sci-RNA over 0–20 h embryos and annotated clusters in each modality. Disk before the fetch: 24 GiB
free.

| file | what is read | bytes | licence |
|---|---|---|---|
| `SupplementaryTables.xlsx` (authors' site, DEAP) | Table S2: ATAC cluster labels per time window; Table S4: ATAC marker accessible regions per cluster | 17,049,285 | not stated |
| `peaks.bed` (DEAP) | the 110,185-peak set, dm6 | 2,532,609 | not stated |
| `rna_meta.rds` (DEAP) | RNA cell to `manual_annot` label, 547,805 cells, 51 labels | 24,110,965 | not stated |
| GEO GSE190147 `scirnaseq_gene_matrix` mtx, rows, columns | UMI counts, 23,932 genes × 547,805 cells, 212,703,167 non-zero entries | 697,773,156 + 88,224 + 19,958,916 | GEO record; no licence stated by the submitters |
| Ensembl 104 `BDGP6.32` GFF3 | protein-coding genes, transcript 5' ends | 5,713,526 | Ensembl states no restrictions; gene models from FlyBase (CC BY 4.0) |

The matrix was summed over cells by `manual_annot`, all windows pooled, the way the source's own
pseudobulk script does it (`scripts/RNA/misc/make_germ_pseudobulk.R` on DEAP), and then deleted: the
store keeps a 1,001,331-byte gene-by-label table. Every file's sha256 is in the inventory. Fetched
767,226,681 bytes; scouting moved 17,768,866 more (the DEAP page, a range read of the matrix header,
and the GFF3 of Ensembl 99, 104 and 113, compared on names only: release 104's 13,968 protein-coding
gene names all appear among the matrix's row names, against 13,837 for 99 and 13,594 for 113). In
all 784,995,547 bytes of the 1.5 GB budget; 70,454,856 kept. The ATAC peak-by-cell matrices (about
15 GB) and the per-cluster bigWigs were not fetched: Table S4 is the source's own per-cluster call.

### The node: decided before any RNA value was read

Drosophila, unlike the worm, has insulator proteins (dCTCF, BEAF-32, CP190, Su(Hw)) and TADs, so a
node in the human form (an element acts inside its boundary-delimited domain) could be built from a
measured boundary set. **Decision: the worm's gene-structure rule, unchanged** (`simpler.domain_rule`,
imported, not copied). No boundary set is used, and none was fetched. The reasons:

1. This run is a replication. If the node changed together with the species, a disagreement with the
   worm could not be attributed to either change.
2. A boundary set is a choice among several (embryo stage or cell line, Hi-C resolution, caller), and
   each choice is a free parameter the worm test did not have. Fly Hi-C boundaries are also reported
   to sit mostly at active promoters (Ramírez et al. 2018, Nat Commun 9:189), so a boundary node and
   a gene-structure node are not independent in the fly. How far they differ was not measured here.
3. The human form of the node has its own registered tests. What is open is whether the worm's
   positive result holds anywhere else.

A boundary-node test would be its own registration, on a boundary set chosen and cited before it is
run. This decision is committed with the frozen code, before the registration and before any RNA
value was read.

### One tissue vocabulary for both modalities

Each modality was annotated separately by the source. The two label sets were grouped by name into
the 14 tissues both modalities name: neural, ventral midline, epidermis, foregut and hindgut, tracheal
system, salivary gland, Malpighian tubule, muscle, fat body, plasmatocytes, midgut, amnioserosa, germ
cell, yolk (`ATAC_TISSUE`, `RNA_TISSUE`; primordia go with their tissue, the source's one
Proventriculus ATAC cluster sits under its Midgut label and the RNA proventriculus follows it).
Labels that name no tissue both modalities resolve are excluded on both sides: germ-layer and region
anlagen or primordia (ectoderm, mesoderm, endoderm, mesectoderm, trunk mesoderm), blastoderm,
maternal, ubiquitous, unknown, head ectoderm (no RNA label), crystal and procrystal cells (no ATAC
label). The code refuses any label that is neither grouped nor excluded. RNA cells per tissue run
from 1,145 (Malpighian tubule) to 76,961 (neural).

### The registered test

Registered in `simpler_fly.REGISTRATION` and `data/results/simpler_fly_registration.json`, committed
before any RNA value was read. `score` refuses to run unless the registration and all five code files
are committed and unchanged, every input matches its registered sha256, and the census recomputes.
The grouping, the population and the class rule were fixed with the ATAC census in view, as the
worm's were; cells per RNA label were seen, no gene's count was.

- **Attribution and baseline.** The worm's domain rule and nearest annotated 5' end, on Ensembl 104
  protein-coding genes of 2L, 2R, 3L, 3R, 4, X and Y (13,950 genes).
- **Population.** Table S4 marker calls with adjusted p < 0.05 and a positive log fold change, joined
  to Table S2's labels on (time window, cluster). A peak is single-tissue when every such call falls
  in clusters of one tissue and none in an excluded label. Distal: midpoint more than 1 kb from every
  annotated transcript 5' end of any biotype. Funnel: 17,422 marker rows (all at adjusted p < 0.05;
  2,361 with a negative fold change), 6,442 peaks with a positive call, 83 with a call in an excluded
  label, 2,566 in two or more tissues, 3,793 in exactly one; of these 2,518 distal, 460 promoter, 815
  neither. Distal elements by tissue: neural 1,330, epidermis 579, muscle 301, midgut 224, yolk 47,
  foregut and hindgut 21, tracheal system 11, plasmatocytes 2, amnioserosa 2, fat body 1; none for
  ventral midline, salivary gland, Malpighian tubule or germ cell. Table S4 lists markers for 76 of
  the source's 181 ATAC clusters, mostly between 6 and 16 h.
- **Where the rules differ.** On 822 of 2,518 (32.6%): 644 inside gene bodies, 171 between genes
  where another gene's 5' end is nearer, 7 convergent. 1,601 elements lie in a host gene, 917 between
  genes. 102 1-Mb clusters.
- **Endpoint.** Per tissue, UMI counts summed over that tissue's cells, all windows pooled, as counts
  per million. Concordant: the attributed gene is at 5 CPM or more in the element's tissue and higher
  there than in each of the other 13. The 5 CPM threshold copies the worm's 5 TPM; it was set before
  any value was read. A protein-coding gene absent from the RNA rows is not concordant.
- **Statistic, matched control, falsifier, decision.** The worm's, unchanged: the difference in
  points with a 95% percentile cluster bootstrap over 1-Mb windows; a matched-distance control gene
  per element; falsified if the domain rule's concordance is not above that control (one-sided
  p ≥ 0.05); DOMAIN_RULE_AHEAD, NEAREST_AHEAD or NO_DIFFERENCE_DETECTED from the interval. The
  replication is read from the same label: ahead replicates the worm's direction, no difference does
  not replicate it, nearest ahead reverses it.
- **Positive control.** The 460 single-tissue promoter elements (midpoint within 500 bp of a
  protein-coding gene's 5' end), nearest-TSS gene against its matched control; if it does not
  separate, the readout is declared unresponsive and the reading withheld.
- **Power in advance.** With 822 differing elements of 2,518, the smallest difference with 80% power
  at two-sided 5% is about 3.2 points before any inflation for clustering. The worm's difference was
  +2.89. A worm-sized difference may therefore go undetected here, and a NO_DIFFERENCE_DETECTED would
  not by itself say the fly differs from the worm.
- **Stop rule.** One run; fewer than 200 differing elements would have stopped it unscored. Nothing
  is changed after an RNA value is read.
- **What it cannot establish.** Accessibility is not function. The two assays come from different
  nuclei; the RNA is pooled over 0–20 h while a marker call may come from one window; the tissue
  labels were assigned separately in each modality and grouped here by name. Inside gene bodies
  accessibility can follow the host's own transcription. Trans-splicing, the worm's outron confound,
  is rare in the fly, but the host confound stands.

### Result (2026-10-10, registration `fd81e37`)

Run once by `scripts/simpler_fly_run.py score` from a clean worktree of `fd81e37` (stores linked
read-only), 10 seconds, no request to any model. Outcome **DOMAIN_RULE_AHEAD**; falsifier
**NOT_FALSIFIED**. The registered reading, verbatim:

> On 2,518 D. melanogaster embryo distal elements accessible in clusters of exactly one of 14 tissues (Calderon et al. 2022, 0-20 h, the source's own marker call), the gene the domain rule attributes is expressed at 5 counts per million or more and highest in the element's own tissue, in RNA pooled over that tissue's cells, for 299 of 2,518 (11.9%), against 231 of 2,518 (9.2%) for the gene with the nearest 5' end; the rules differ on 822 elements; difference +2.70 points, 95% +1.32 to +4.10 over 1-Mb clusters. On this endpoint the domain rule is ahead of nearest TSS. The worm's result (+2.89 points on Serizay et al. 2020) is replicated in direction on a second species and source. The domain rule's genes match the element's tissue more often than genes at a matched distance (2.8%, one-sided p 0.001), so position carries tissue information here. The difference is not carried by the intergenic elements; inside gene bodies accessibility can follow the host gene's own transcription, so this run does not separate regulation from transcription. Accessibility is not function: this says the attributed gene is expressed where the element is open, not that the element acts on it, and the two assays come from different nuclei.

The positive control separates, so the readout is responsive: of 460 single-tissue promoter
elements, the nearest-TSS gene is concordant for 128 (27.8%) against 5.5% for genes at a matched
distance (one-sided p 0.001, the floor of 1,000 draws). The nearest-TSS gene on the distal elements
is also above its own matched control (9.2% against 3.6%, p 0.001). All 13,950 protein-coding genes
are in the RNA rows, so no attributed gene is missing.

Beside the primary, deciding nothing (registered strata; percentile intervals over 1-Mb clusters;
strata of fewer than 50 elements are not shown, and their intervals mean little):

| stratum | elements | domain rule | nearest TSS | difference, points [95%] |
|---|---|---|---|---|
| inside a gene body | 1,601 | 14.2% | 10.5% | +3.69 [+1.82, +5.72] |
| between genes | 917 | 7.9% | 6.9% | +0.98 [-0.42, +2.31] |
| the 822 where the rules differ | 822 | 11.4% | 3.2% | +8.27 [+3.93, +12.62] |
| neural | 1,330 | 11.9% | 9.2% | +2.63 [+0.48, +4.57] |
| epidermis | 579 | 9.8% | 6.7% | +3.11 [+0.20, +6.83] |
| muscle | 301 | 9.6% | 8.3% | +1.33 [-0.64, +3.97] |
| midgut | 224 | 10.7% | 9.4% | +1.34 [-0.53, +3.90] |

**What this says and what it does not.** The difference has the worm's direction (+2.70 here,
+2.89 in the worm, and the two intervals overlap), on a second species, a second source, a different
assay (single-cell rather than sorted nuclei) and an endpoint that asks for the highest of 14 tissues
rather than of 5.
The concordance rates are low because that endpoint is strict: 2.8% for a gene at a matched
distance, 11.9% for the domain rule.

The fly differs from the worm in where the difference sits. In the worm it held between genes
(+3.73); in the fly it does not (+0.98, interval spanning zero), and the registered host caveat is
part of the reading. The difference is carried by elements inside gene bodies, where accessibility
can follow the host gene's own transcription. In the fly this run therefore does not separate
regulation from transcription. Trans-splicing, which confounded the worm's between-gene stratum, is
rare in the fly, so the two species' caveats are not the same caveat. Two of the four tissues with
more than 200 elements (muscle, midgut) show no difference on their own.

Accessibility is not function. Nothing here says that any element acts on the gene it is attributed
to. It is one source, pooled over 0–20 h, with accessibility and RNA from different nuclei.

### What follows from the fly

1. A boundary node, registered on its own: a measured boundary set (stage, cell source and caller
   named and cited before the run) against the gene-structure rule, on this population.
2. The between-gene question on a second fly source: the worm and the fly disagree on whether the
   difference survives outside gene bodies. A registration that makes the between-gene stratum its
   primary, with the window of the marker call matched to the window of the RNA, would ask it
   directly.

## 9. Drosophila: the boundary node, registered on Stadler 2017 boundaries (2026-10-10, lane-flynode)

Section 8 left the fly's result inside gene bodies: the worm's gene-structure rule beat nearest TSS on
all 2,518 distal elements (+2.70) but not between genes (+0.98 [-0.42, +2.31], 917 elements), so that
run could not separate regulation from the host gene's transcription. The human node is a
boundary-delimited domain, and the fly, unlike the worm, has measured boundaries. This registration
asks the human form of the question on the same elements and endpoint. No money, no key, no model
request.

Code: `genomeos/attribution/simpler_flynode.py` (imports `simpler.py` and `simpler_fly.py`, both
unchanged), `scripts/simpler_flynode_run.py` (`fetch`, `inventory`, `register`, `score`),
`tests/test_simpler_flynode.py`. Results: `data/results/simpler_flynode_inventory.json`,
`simpler_flynode_registration.json`, `simpler_flynode_result.json`, each written from a clean
worktree of the committed tree with the stores linked read-only. Commits: code frozen `6bee65f` and
`90a48c7` (the registration's manifest declares the RNA file it hashes), registration `f5985f7`.

### The boundary set (fetched 2026-10-10)

| file | what is read | bytes | sha256 | licence |
|---|---|---|---|---|
| `elife-29550-supp1-v2.xlsx` (eLife CDN) | Tables S4 (manual), S5 (computational), S7 (representative) boundary lists, dm6 columns | 386,317 | `e944d828...49db0e4` | CC BY 4.0 |

Stadler, Haines and Eisen 2017 (eLife 6:e29550) made Hi-C on hand-sorted nuclear-cycle-14 embryos
(2-3 h) at 500-bp resolution and called boundaries by eye over the whole genome (Table S4, 3,121
distinct boundaries on the used arms), by a directionality score (Table S5, 2,666) and as a merge
(Table S7, 954). Disk before the fetch: 27,064,115,200 bytes free. Scouting moved 11,639,184 bytes more (the GEO
copy of the same tables, a scouting copy of this file, Cavalheiro 2023's supplementary archive),
12,025,501 in all of the 200 MB budget.

**Why this set.** Embryos; dm6 as the authors give it, so nothing is lifted; 500-bp resolution,
finer than the spacing of fly genes; boundaries called by the authors, so no caller or threshold is
chosen here; a stated licence; under 0.4 MB. Not used: Cavalheiro et al. 2023 (Sci Adv, CC BY;
10-kb bins, and its table lists 878 boundaries of a differential analysis against the 772 its text
reports), the GEO copy of Stadler's tables (dm3 only), Ogiyama et al. 2018 (contact matrices, no
boundary list: calling one would choose a caller here), Hug et al. 2017 (not openly licensed),
Ramirez et al. 2018 (Kc167 cells). **Stage:** the boundaries are from 2-3 h embryos, the marker calls
mostly from 6-16 h, the RNA from 0-20 h. That boundaries formed at zygotic genome activation persist
is the sources' report (Hug 2017, Cavalheiro 2023), not checked here.

### The node, decided before any RNA value was read

An element's domain runs from the boundary at or left of its midpoint to the next boundary right of
it, or to the arm end; a gene is in the domain when its 5' end lies there or inside either
boundary's 500-bp bin (a boundary is the bin start + 250, and the bins are shared, since fly
boundaries sit mostly at promoters). The rule names the in-domain gene with the nearest 5' end. When
the arm-wide nearest 5' end is in the domain it is nearest TSS itself; when the domain holds no 5' end
it falls back to nearest TSS. So it differs from nearest TSS exactly where the nearest 5' end lies
across a measured boundary and the domain holds another.

### The census, and the population it fixed

The lane's brief named the 917 between-gene elements as the primary population. The census, from
positions and the element-side tissue call only, showed that between genes the boundary rule changes
the nearest-TSS gene for **90 of 917**, under the stop rule's 200 (median domain around a primary
element 84 kb; 17 of 2,518 elements sit inside a boundary bin). Between genes this boundary set cannot
decide. The primary population was therefore fixed, with that census in view and no RNA value read,
as the **1,634 distal elements outside the body of the gene the boundary rule names**: all 917
between-gene elements and 717 elements inside a gene body for which the rule names another gene. On
it the named gene never contains the element, which is the brief's reason, and the rules differ on
234 (for 115 of these the nearest-TSS gene contains the element, which can favour nearest TSS). Between
genes and all 2,518 distal elements are reported beside. On all distal elements the rules differ on
273 of 2,518; for the computational list 301 of 1,697 primary elements, for the representative 86 of
1,585.

### The registered test

Registered in `simpler_flynode.REGISTRATION` and `data/results/simpler_flynode_registration.json`
(`f5985f7`). `score` refuses to run unless the registration, simpler_fly's registration and the five
code files are committed and unchanged, every input matches its registered sha256, simpler_fly's own
census recomputes, and this census recomputes.

- **Outcome.** Boundary rule against nearest TSS on the 1,634: simpler_fly's statistic, percentile
  cluster bootstrap over 1-Mb windows; BOUNDARY_RULE_AHEAD, NEAREST_AHEAD or NO_DIFFERENCE_DETECTED.
- **Second registered comparison.** Boundary rule against the gene-structure rule on the same
  elements, with its own label, not the outcome; the gene-structure rule names the host for 717 of
  them, where its transcription can favour it. Two comparisons on one population, each at 95%, no
  correction for the pair.
- **Endpoint, matched control, falsifier, positive control.** simpler_fly's, unchanged; the matched
  control and falsifier are applied to the boundary rule's gene on the primary population.
- **Beside, deciding nothing.** The same boundaries rotated along each arm by one uniform offset
  (number and spacing kept, 1,000 draws); the computational and representative lists; the elements
  where the rules differ alone; tissues with 50 or more primary elements; the gene-structure rule
  against nearest TSS.
- **Power in advance.** 234 differing elements of 1,634: about 2.6 points at 80% power, two-sided
  5%, before clustering.
- **Stop rule.** One run; fewer than 200 differing primary elements would have stopped it unscored.

### Result (2026-10-10, registration `f5985f7`)

Run once by `scripts/simpler_flynode_run.py score` from a clean worktree of `f5985f7` (stores linked
read-only), 12 seconds, no request to any model. Outcome **NO_DIFFERENCE_DETECTED**; second
comparison **GENE_STRUCTURE_AHEAD**; falsifier **NOT_FALSIFIED**. The registered reading, verbatim:

> On 1,634 D. melanogaster embryo distal elements accessible in clusters of exactly one of 14 tissues (Calderon et al. 2022, 0-20 h, the source's own marker call) and lying outside the body of the gene the boundary rule names, the gene with the nearest 5' end inside the element's domain, delimited by 3,121 Hi-C boundaries of nuclear-cycle-14 embryos (Stadler et al. 2017, the authors' manual call), is expressed at 5 counts per million or more and highest in the element's own tissue, in RNA pooled over that tissue's cells, for 88 of 1,634 (5.4%), against 84 of 1,634 (5.1%) for the gene with the nearest 5' end on the arm; the rules differ on 234 of 1,634 elements; difference +0.24 points, 95% -0.34 to +0.95 over 1-Mb clusters. No difference detected between the boundary rule and nearest TSS. With 234 elements differing, a difference smaller than about 2.6 points could go undetected. For 115 of the 234 elements where the rules differ, the nearest-TSS gene contains the element, which can favour nearest TSS. Against the gene-structure rule on the same elements (148 of 1,634, 9.1%; the two differ on 943, and the gene-structure rule names the host for 717), the difference is -3.67 points, 95% -5.96 to -1.37: the gene-structure rule is ahead of the boundary rule. The boundary rule's genes match the element's tissue more often than genes at a matched distance (4.4%, one-sided p 0.001), so position carries tissue information here. With the same boundaries rotated along each arm the rule averages 5.5% on these elements (one-sided p 0.612), so the measured positions are not shown to carry more than the boundaries' number and spacing. Between genes (beside, deciding nothing) the boundary rule names a different gene from nearest TSS for only 90 of 917 elements, fewer than the stop rule's 200, so this boundary set cannot decide the between-gene question: there it is concordant for 64 (7.0%) against 63 (6.9%) for nearest TSS, +0.11 points [-0.46, +0.70]. On all 2,518 distal elements (beside, deciding nothing; the rules differ on 273) the boundary rule is concordant for 238 (9.5%) against 231 (9.2%) for nearest TSS, +0.28 points [-0.22, +0.82], and against 299 (11.9%) for the gene-structure rule, -2.42 points [-3.85, -0.94]. Accessibility is not function, and a boundary is not a regulatory unit by itself: this says the attributed gene is expressed where the element is open, not that the element acts on it. The boundaries come from 2-3 h embryos, the marker calls mostly from 6-16 h, and accessibility and RNA from different nuclei.

The positive control separates, so the readout is responsive: of 460 single-tissue promoter
elements the nearest-TSS gene is concordant for 128 (27.8%) against 5.5% at a matched distance
(one-sided p 0.001). No attributed gene is missing from the RNA rows.

Beside the outcome, deciding nothing (percentile intervals over 1-Mb clusters):

| comparison | elements (differ) | first rule | second rule | difference, points [95%] |
|---|---|---|---|---|
| boundary vs nearest, where they differ | 234 (234) | 6.8% | 5.1% | +1.71 [-2.50, +7.01] |
| boundary vs nearest, computational list | 1,697 (301) | 4.8% | 5.3% | -0.47 [-1.23, +0.22] |
| boundary vs nearest, representative list | 1,585 (86) | 5.4% | 5.4% | +0.00 [-0.36, +0.34] |
| boundary vs gene structure, between genes | 917 (226) | 7.0% | 7.9% | -0.87 [-2.36, +0.63] |
| gene structure vs nearest, primary | 1,634 (781) | 9.1% | 5.1% | +3.92 [+1.76, +6.02] |
| boundary vs nearest, neural | 929 (138) | 6.6% | 6.1% | +0.43 [-0.49, +1.64] |
| boundary vs nearest, epidermis | 354 (47) | 1.1% | 1.1% | +0.00 [-0.83, +0.87] |
| boundary vs nearest, muscle | 195 (29) | 6.2% | 5.6% | +0.51 [+0.00, +1.66] |
| boundary vs nearest, midgut | 110 (13) | 7.3% | 6.4% | +0.91 [+0.00, +3.00] |

**What this says and what it does not.** No difference was detected between the boundary rule and
nearest TSS on this endpoint: where the rules differ (234 elements) the in-domain gene is concordant
for 16 and the across-boundary nearest gene for 12, and the difference over all 1,634 is +0.24
[-0.34, +0.95]. A difference under about 2.6
points could go undetected. The rotated boundaries do as well as the measured ones (5.5% against
5.4%), so the measured positions are not shown to carry information here beyond how many boundaries
there are and how they are spaced. Neither of the other two lists changes that.

The gene-structure rule is ahead of the boundary rule on the same elements (-3.67), but on 717 of them
it names the host, whose transcription can favour it; between genes, where no rule names a host, the
two nodes are not told apart (-0.87 [-2.36, +0.63]). So this run does not rank the two nodes on
regulation.

The between-gene question the brief asked stays open: on this boundary set the boundary rule names
the nearest-TSS gene for 827 of 917 between-gene elements, so it cannot separate regulation from
transcription there. The boundaries come from 2-3 h embryos and the marker calls mostly from 6-16 h;
a later-stage boundary set could place boundaries differently, and was not available as an openly
licensed, author-called dm6 list.

Accessibility is not function, and a boundary is not a regulatory unit by itself. Nothing here says
that any element acts on the gene it is attributed to.

### What follows from the boundary node

1. The between-gene question needs a population where the rules differ often enough: a second fly
   source with more between-gene elements (or the marker window matched to the RNA window, section 8
   item 2), not another node on this one.
2. A boundary set from the marker calls' own stages (6-16 h embryos), author-called and dm6, if one is
   published under a stated licence; none was found on 2026-10-10.

## 10. Stage 3: the worm's constraint-dose view, registered on a second source (2026-10-10, lane-wormdose)

Item 1 of "What follows from stage 2". Stage 2's registered tier was FALSIFIED on Serizay 2020, and a
view beside it, deciding nothing, was positive: the top third of elements by constrained fraction
against the bottom third gave +6.69 [+2.11, +11.68]. That view was seen after stage 2's run, so it
cannot be re-tested on that source. Here it is registered as a hypothesis on a second worm source
before any of that source's RNA was read, with stage 1's comparison as a secondary. No money, no key,
no model request.

Code: `genomeos/attribution/simpler_dose.py` (imports `simpler.py` and `simpler_tiers.py`, both
unchanged), `scripts/simpler_dose_run.py` (`fetch`, `inventory`, `tiers`, `register`, `score`),
`tests/test_simpler_dose.py`. Results: `data/results/simpler_dose_inventory.json`,
`simpler_dose_tiers.json`, `simpler_dose_registration.json` and `simpler_dose_result.json` (registration
1 and its run), `simpler_dose_v2_registration.json` and `simpler_dose_v2_result.json` (registration 2
and its run), `simpler_dose_v2_edges_posthoc.json` (post hoc, deciding nothing), each written from a
clean worktree of the committed tree with the stores linked read-only. Commits: code `91c2e75`,
inventory and per-element counts `e4b794a`, registration 1 `ff74444`, run 1 `b8d777a`, code for
registration 2 `7a1b836`, registration 2 `5b930cb`, run 2 `c2adfa3`, post-hoc check `cc86cec` and
`5807027`.

### The source (fetched 2026-10-10)

Jänes et al. 2018 (eLife 7:e37344) made ATAC-seq of whole wild-type animals at six stages (embryo, L1,
L2, L3, L4, young adult) and, at the same six stages, strand-specific long nuclear RNA
(chromatin-associated, over 200 nt, "long cap"), two biological replicates each. It is the same lab as
Serizay 2020 but different samples: whole animals by stage, not sorted nuclei by tissue. The source does
not say whether the ATAC and the RNA came from the same nuclei preparations.

| file | what is read | bytes | sha256 | licence |
|---|---|---|---|---|
| `elife-37344-fig1-data1-v2.txt` | coordinates and the six wild-type ATAC heights (SPMR peak signal) | 3,835,944 | `18a4e546...77a05fc` | CC BY 4.0 |
| `elife-37344-fig2-data1-v2.txt` | ce10 and ce11 coordinates and the element class only | 18,827,378 | `09bede4b...80e4656` | CC BY 4.0 |
| `elife-37344-fig4-data1-v2.txt` | coordinates and the developmental-dynamics call only | 4,739,367 | `91cd4e53...7ac2bd2` | CC BY 4.0 |
| `GSE114483_lcap_wt_<stage>_linear_<fwd,rev>.bw`, 12 files | long-cap coverage, DESeq2 size-factor normalised, averaged over two replicates | 349,972,704 | per file in the inventory | GEO; NCBI places no restriction, the submitters state none |

377,375,393 bytes fetched against the 1.5 GB budget, plus 272,291 bytes of the eLife API record while
choosing; 26.9 GB free before the fetch. The phyloP135way track was read by range request
(61,558,764 bytes moved, nothing stored). Not used: Durham et al. 2021 (single-cell ATAC at L2, but
its cell-type RNA would come from another study, so the endpoint would not be measured with the
accessibility), Daugherty et al. 2017 (three stages; RNA in the same samples not established), GEO
GSE141213 (2,386,268,160 bytes, over budget, and the same source as stages 1 and 2).

### Adaptations, decided before any RNA value was read

- **The group.** Stage 1's group was a tissue the source called ("accessible in exactly one of five
  tissues"). This source gives no per-stage call. Among the putative enhancers the source calls
  developmentally dynamic (`devel_is_dynamic`), an element's stage is the one with its highest ATAC
  height of the six; a tie goes to the earlier stage (46 ties among 24,051 dynamic elements of both
  classes).
- **The RNA.** The source publishes no gene table, so a gene's value at a stage is coverage-TPM: the
  mean long-cap coverage over the union of the exons of all its WBcel235 transcripts, on its own
  strand, scaled per stage so the protein-coding genes of the six nuclear chromosomes sum to one
  million. That is TPM's construction (reads per base, normalised to a million) on coverage instead
  of counts.
- **The endpoint.** Concordant: at least 5 coverage-TPM at the element's stage and strictly higher
  than at each of the other five (stage 1's rule on six stages instead of five tissues).
- **A strand check that can fail the run.** Own-strand coverage summed over the protein-coding genes
  must exceed opposite-strand at all six stages (read as labelled) or fall below it at all six (fwd
  and rev swapped); anything else is INSTRUMENT_FAILED.

The source's own element-to-gene columns and its RNA-derived flag columns are never read.

### The census, from positions, the ATAC calls and the conservation track only

13,852 developmentally dynamic putative enhancers with a domain-rule gene (embryo 4,735, L1 974, L2
1,955, L3 1,160, L4 1,958, young adult 3,070); 11,842 have a scored phyloP135way base left after exon
removal. Tertile edges of the constrained fraction: 0.072848 and 0.19403 (stage 2's were 0.08 and
0.190476). Top third 3,949, bottom third 3,922; in strata holding both arms (127 of 136), 3,941 and
3,905; the smallest difference with 80% power at two-sided 5% is 3.16 points before clustering.
**Overlap:** 3,849 of the 13,852 elements (3,407 of the 11,842 measured) overlap an element of stage
2's population, because Serizay 2020 extends this lab's annotation; without them 2,760 top and 2,902
bottom remain in support. The arms are unequal across stages (young adult 314 top against 1,266
bottom; embryo 1,743 against 1,145), which the stage stratum absorbs. The secondary's rules differ on
2,919 of 13,852 elements (21.1%). Positive control: 10,199 dynamic coding promoters.

### The registered test

Registered in `simpler_dose.REGISTRATION` and `data/results/simpler_dose_registration.json`
(`ff74444`), which hashes all fifteen inputs; registration 2 (below) changes one line of it. `score` refuses unless the registration, the per-element
counts, stage 2's tiers file and the five code files are committed and unchanged, every input matches
its registered sha256 and the census recomputes.

- **Claim.** Among the putative enhancers, the domain-rule gene is concordant with the element's own
  stage more often for the top third by constrained fraction than for the bottom third, at the same
  stage, position class, distance to that gene and measured length.
- **Constraint.** Stage 2's instrument unchanged: phyloP135way by range request, cut 2.08, exon bases
  removed; fraction = bases at or above the cut / scored bases left.
- **Matched control.** The bottom third within the same strata: stage (6) × inside a gene body or
  between genes × distance bin of the domain-rule gene's 5' end (stage 1's 8 bins) × measured-length
  tertile.
- **Statistic and decision.** `simpler_tiers.matched_difference`, unchanged: standardised to the top
  third's strata, 95% percentile interval over 1-Mb windows (2,000 resamples); DOSE_ADDS,
  NO_DIFFERENCE_DETECTED or DOSE_BELOW from the interval.
- **Falsifier.** The arm label shuffled within strata 1,000 times; FALSIFIED if one-sided p ≥ 0.05.
- **Positive control.** Dynamic coding promoters, nearest-TSS gene, against a gene at the same
  distance bin (stage 1's control); if it does not separate, the readout is unresponsive and neither
  reading is given.
- **Secondary, its own label.** Stage 1's comparison on this source: domain rule against nearest TSS
  on all 13,852 elements, stage 1's statistic, matched control and falsifier.
- **Beside, deciding nothing.** Stage 2's numeric edges; the elements overlapping none of stage 2's
  population (a DOSE_ADDS reading carries a registered caveat if that subset's interval does not exclude
  zero); inside gene bodies and between genes; each stage; the middle third against the bottom; stage
  2's 5% tier on this source.
- **Stop rule.** One run; fewer than 200 per arm in support stops it unscored.
- **Confound, stated in advance.** Stages 1 and 2 named transcription; whole-animal RNA adds that a
  stage concordance can come from any tissue. A difference does not separate regulatory constraint from
  transcription-linked constraint.

### Run 1 (registration `ff74444`): stopped at the strand check

Run once by `score` from a clean worktree of `ff74444`, 81 seconds of RNA reading. Outcome
**INSTRUMENT_FAILED**; nothing was scored. The registered reading, verbatim:

> INSTRUMENT_FAILED: the long-cap RNA tracks show neither own-strand coverage above opposite-strand at all six stages nor the reverse at all six, so the RNA's strand cannot be read; nothing is scored.

The per-stage sums in the result were negative wherever a rev track entered. Counted per track over
the six chromosomes, every non-zero value of the six fwd tracks is positive and every non-zero value of
the six rev tracks is negative (43,188,260 to 61,500,850 bases per rev track): the GEO rev tracks store
coverage as negative numbers, a display convention the registration had not anticipated.

### Registration 2 (`5b930cb`): the same test, coverage read as its magnitude

`simpler_dose.REGISTRATION_V2` is registration 1 with one change: each track value is read as its
absolute value, and the strand check applies to those. A test pins that nothing else differs
(population, census, constraint, edges, strata, endpoint, controls, falsifiers, stop rule, wording).
What was seen of the RNA before it was written, recorded in the registration file: run 1's twelve
per-stage sums and the sign counts above. No gene value was saved or read and no concordance was
computed before run 2. `score2` refuses unless run 1's result and registration 2 are committed, and
runs once.

### Result (2026-10-10, registration `5b930cb`)

Run once by `score2` from a clean worktree of `5b930cb`, 86 seconds of RNA reading, no request to any
model. Strand check passed as labelled at all six stages (own-strand sums 1.2 to 3.5 million against
25,000 to 54,000 opposite-strand). Outcome **NO_DIFFERENCE_DETECTED**; falsifier **NOT_FALSIFIED**;
secondary **DOMAIN_RULE_AHEAD**. The registered reading, verbatim:

> On C. elegans putative enhancers whose accessibility changes across development (Jänes et al. 2018, whole animals, six stages from embryo to young adult), each assigned to its stage of highest accessibility, with constraint measured on their exon-free bases (UCSC ce11 phyloP135way, stage 2's cut 2.08), the domain-rule gene is expressed at 5 coverage-TPM or more in long nuclear RNA and highest at the element's own stage for 1,631 of 3,938 elements in the top third by constrained fraction (41.4%), against 38.6% for the 4,094 in the bottom third at the same stage, position class, distance and measured length; difference +2.81 points, 95% -0.17 to +5.71 over 1-Mb clusters (within-stratum permutation p 0.00999; NOT_FALSIFIED). No difference was detected: the dose view seen beside stage 2 does not replicate on this source. This is not evidence that constrained elements lack function. The positive control separates: of 10,199 developmentally dynamic coding promoters, the nearest-TSS gene is stage-concordant for 3,936 (38.6%) against 22.3% at a matched distance (one-sided p 0.001). Secondary, stage 1's comparison on this source: on 13,852 elements the domain-rule gene is stage-concordant for 4,686 (33.8%) against 4,400 (31.8%) for the gene with the nearest 5' end; the rules differ on 2,919; difference +2.06 points, 95% +1.33 to +2.82 over 1-Mb clusters. On this endpoint the domain rule is ahead of nearest TSS, as in stage 1. The domain rule's genes match the element's stage more often than genes at a matched distance (19.3%, one-sided p 0.001). Conservation and accessibility are not function: this says the attributed gene is expressed when the element is open, not that the element acts on it, and it does not separate regulatory constraint from transcription-linked constraint. The RNA is from whole animals, so a stage concordance can come from any tissue.

Every domain-rule gene has an RNA value (19,973 protein-coding genes). Beside the registered test,
deciding nothing (matched the same way; percentile intervals over 1-Mb clusters; within-stratum
permutation p):

| view | top / bottom, in support | top | matched bottom | difference, points [95%] | p |
|---|---|---|---|---|---|
| stage 2's numeric edges (≥ 0.190476 against < 0.08) | 4,091 / 4,403 | 41.2% | 38.2% | +2.97 [+0.08, +5.83] | 0.005 |
| elements overlapping none of stage 2's population | 2,757 / 3,033 | 44.0% | 38.4% | +5.59 [+2.20, +9.00] | 0.001 |
| inside a gene body | 2,214 / 1,815 | 42.5% | 37.8% | +4.78 [+0.31, +9.10] | 0.004 |
| between genes | 1,724 / 2,279 | 40.0% | 39.7% | +0.28 [-3.19, +3.84] | 0.44 |
| middle third against bottom third | 3,776 / 4,097 | 31.7% | 31.5% | +0.24 [-1.69, +2.38] | 0.40 |
| stage 2's 5% tier on this source | 8,888 / 2,930 | 35.4% | 33.3% | +2.13 [+0.14, +4.40] | 0.019 |
| embryo | 1,739 / 1,192 | 59.6% | 49.3% | +10.24 [+4.42, +15.96] | 0.001 |
| L1 | 250 / 317 | 32.0% | 27.5% | +4.53 [-4.96, +12.12] | 0.18 |
| L2 | 819 / 376 | 48.6% | 54.5% | -5.94 [-13.92, +2.95] | 0.97 |
| L3 | 318 / 326 | 11.3% | 7.7% | +3.62 [-2.26, +8.80] | 0.087 |
| L4 | 498 / 573 | 8.4% | 12.4% | -3.95 [-9.26, +1.17] | 0.97 |
| young adult | 314 / 1,310 | 12.4% | 19.4% | -6.96 [-12.12, -1.77] | 0.995 |

The secondary inside gene bodies: +2.26 [+1.08, +3.53] on 7,701 elements; between genes: +1.82
[+1.10, +2.54] on 6,151 (the rules differ on 794 of them).

**A defect of this run, disclosed.** The run read the tertile edges as the census stored them,
rounded to six decimals (0.072848 and 0.19403). The census had counted its arms on the exact edges,
11/151 and 13/67. The stored lower edge lies above 11/151, so the 190 elements whose fraction is
exactly 11/151 (189 of them in strata holding both arms) entered the bottom third in the run, while
the census had them in the middle; 3
elements at 13/67 left the top third. That is why the run's arms in support (3,938 and 4,094) differ
from the census's (3,941 and 3,905). Stage 2's dose view was computed the same way, from its stored
edges. Recomputed post hoc on the exact edges (`simpler_dose_v2_edges_posthoc.json`, deciding nothing):
+3.09 [+0.12, +6.13], p 0.004. The registered outcome stays run 2's. The decision turns on those 190
elements, so neither reading is firm. Future censuses should store the edges unrounded.

**What this says and what it does not.** On a second source the dose view is not shown to replicate
under the registered decision: the difference is +2.81 points and its interval reaches zero
(-0.17). On the whole population the three edge choices give +2.81, +2.97 and +3.09 points, under half
of stage 2's +6.69, with lower bounds of -0.17, +0.08 and +0.12. The permutation, which does not
resample the 1-Mb clusters, gives p 0.00999. The views beside it are consistent with the confound
named in advance: the difference sits inside gene bodies (+4.78) and not between genes (+0.28), and
the stages disagree (embryo +10.24, young adult -6.96). The middle third does no better than the
bottom third (+0.24), so the difference is not graded across the lower two thirds. The elements stage
2 did not hold give +5.59 [+2.20, +9.00], so the estimate is not carried by the elements stage 2 held.

On this endpoint the domain rule is ahead of nearest TSS on this source, as in stage 1: +2.06 [+1.33,
+2.82], between genes (+1.82) as well as inside them (+2.26), and its genes match the element's stage
more often than genes at a matched distance. The endpoint, stage 1's adapted to six stages, responds: the promoters'
nearest genes separate from matched genes (38.6% against 22.3%).

Conservation and accessibility are not function. Nothing here says that any element acts on the gene
it is attributed to, and no view separates regulatory constraint from transcription-linked
constraint.

### What follows from stage 3

1. Whether constraint carries information between genes, where transcription of a host cannot: the
   between-gene stratum here is +0.28 [-3.19, +3.84]. A test of it needs its own registration and a
   population that does not read this endpoint.
2. Censuses that store tertile edges unrounded, so the run and the census count the same arms.
