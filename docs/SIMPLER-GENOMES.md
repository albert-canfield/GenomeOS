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
