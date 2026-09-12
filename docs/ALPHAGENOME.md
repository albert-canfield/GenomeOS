# AlphaGenome in GenomeOS: what it gives us, and what we do with it

Researched 2026-09-11 from the AlphaGenome documentation, the client source,
the Nature paper (Avsec et al., *Nature* 649, 1206-1218, 28 January 2026) and
the model card. Every claim below is from those sources. The feature is
optional and loads disabled without a key (docs/DECISIONS.md D34); how to
enable it is in [DATA.md](DATA.md).

## Why it matters here

GenomeOS decodes the genome into blocks and says what each one is. For coding
blocks the answer is measured: a gene, its transcripts, its protein, verified
against the curators. For the 96.6% of UNKNOWN space we classify by sequence
and by chromatin evidence, the honest answer usually stops at *what kind of
thing this is* (regulatory, repeat, gap) and never reaches *what it does*.
Enhancer targets, in particular, are inferred from the CTCF domain and
labelled inferred at confidence 0.4 (D21), because nothing we have measures
which gene an element reaches.

AlphaGenome is the first model that answers that question at the right scale.
It reads up to 1 Mb of sequence and predicts, at base resolution, what the
cell does with it: expression, splicing, accessibility, transcription-factor
occupancy, and the contact map. That is exactly the layer between our blocks
and our runtime, and it enters as `predicted` evidence, never as measurement.

## What it actually provides

| | |
|---|---|
| Input | one sequence of 16 kb, 100 kb, 500 kb or 1 Mb; human hg38 and mouse mm10; a reference interval, a custom sequence, or a variant |
| Outputs | 11 modalities: RNA-seq, CAGE, PRO-cap, ATAC, DNase, ChIP histone, ChIP transcription factor, splice sites, splice-site usage, splice junctions, contact maps |
| Resolution | 1 bp for expression, accessibility and splicing; 128 bp for ChIP; 2,048 bp for contact maps |
| Tracks | 5,930 human (and 1,128 mouse), each annotated with assay, biosample name, ontology CURIE, GTEx tissue, life stage and data source, so a prediction can be asked of a named tissue or cell type |
| Variant scoring | 19 recommended scorers; results come back as genes × tracks with a raw effect size and a quantile score ranked against about 300,000 common variants |
| Interpretation | in-silico mutagenesis: every alternative base across a window, aggregated into a contribution matrix |

Practical limits that shape the design: the service takes one variant per
request (the client fans out over five workers), the documentation describes
it as suited to thousands of predictions and not millions, and quotas are
unpublished and vary with demand. Results may be saved locally; they may not
be used to train another model; use is non-commercial and never clinical.

Stated limitations, which we repeat wherever we use it: elements more than
100 kb from the gene are weak, cell-type specificity is harder than average
effect, personal genomes are a known weakness, and the model reads one
unphased sequence, so it does not see a diploid individual.

## Seven uses, in the order they are worth building

1. **Variant effect per tissue** (built, feature a). `genomeos predict
   chr:pos REF>ALT` returns the predicted expression change of every gene in
   the window, per tissue track, as `predicted` rules capped at confidence
   0.7. Verified against a published liver eQTL: rs12740374 raises PSRC1 by
   0.81 log2 and CELSR2 by 0.51 in liver tracks, while the coding variant APP
   A673T shows no expression effect at all (`data/results/alphagenome_rs12740374.json`).
2. **Enhancer to gene** (built, feature b). For an enhancer-like element the
   registry says where and what kind; nothing says which gene. `genomeos predict
   --element chr21:START-END` deletes the element as a variant in its 1 Mb
   window and reads the predicted expression of every gene in the window
   across the 371 RNA-seq tracks; the gene that moves most is the predicted
   target, the track where it moves most the tissue, the change the magnitude,
   and a rise on deletion marks the element as silencer-like for that gene.
   Everything is `predicted`, confidence capped at 0.7; answers are cached per
   element under `data/knowledge/alphagenome/elements` so a layer that only
   reads (`genomeos regulation`, the Flow tab) never calls the API. The
   chromosome job `enhancer_targets_chr21` (Progress tab) sampled 200 of the
   6,618 distal enhancers of chr21 evenly along the chromosome, about eight
   seconds each, and holds each answer against the domain inference
   (`data/results/enhancer_targets_chr21.json`):

   | | |
   |---|---|
   | elements deleted | 200 (dELS, 2 kb to 400 kb from the inferred target) |
   | some gene moves by ≥ 0.1 log2 | 127 (63.5%); 40 strong (≥ 0.3), 87 weak |
   | strongest coding gene = nearest TSS in the CTCF node | 61 of 90 named coding targets (67.8%) |
   | strongest coding gene inside the node | 78 of 90 (86.7%); 12 name a gene beyond the boundary |
   | no coding gene moves | 110 (55%); a non-coding gene moves in 37 of these |
   | silencer-like (expression rises on deletion) | 43 of 127 named |
   | median distance when prediction and inference agree | 28 kb |
   | some gene named, element < 20 kb from its inferred target | 52% (n = 63) |
   | some gene named, element ≥ 100 kb from its inferred target | 32% (n = 31) |

   Read carefully: the 87% inside-the-node figure is a prediction agreeing with
   an inference, not a measurement, and the model itself states that elements
   beyond 100 kb are weak, which the last two rows reproduce. What it changes
   for a geneticist is the label on an enhancer: 61 elements move from
   "nearest coding gene in the node, inferred 0.4" to a named gene with a
   tissue and a magnitude, and 43 elements that the registry calls
   enhancer-like behave as silencers for the gene they move. The examples
   worth looking at are PKNOX1 (element EH38E3461530, 56 kb away, −0.80 in a
   neuroectodermal line, agrees) and the elements the inference gets wrong:
   EH38E3454846 is 41 kb from SCAF4 but moves HUNK, beyond the boundary;
   EH38E3457653 sits 234 kb from HLCS and moves SIM2 instead. Other
   chromosomes run through the same job by name (`enhancer_targets_<chrom>`,
   registered for every chromosome); the summary card on the Progress tab
   shows one row per chromosome scored. All 24 chromosomes ran the same day
   (`enhancer_targets_genome_wide`, `data/results/enhancer_targets_genome_wide.json`):

   | | genome-wide (24 chromosomes) |
   |---|---|
   | distal enhancers deleted one by one | 4,800 (200 per chromosome) |
   | some gene moves by ≥ 0.1 log2 | 2,994 (62.4%); 735 strong |
   | a coding gene moves | 2,291 |
   | strongest coding gene = nearest TSS in the CTCF node | 1,631 (71.2%; per chromosome 50% to 83%) |
   | strongest coding gene inside the node | 2,066 (90.2%; per chromosome 79.8% to 91.8%, median 86.4%) |
   | strongest coding gene beyond the boundary | 225 (9.8%) |
   | silencer-like (expression rises on deletion) | 1,157 of 2,994 named |
   | tissues where the strongest effect falls most often | K562, placenta, CD14 monocyte, HepG2, small intestine |

   The chr21 numbers above are one row of this; the shape does not change
   across chromosomes, which is the point: nine in ten predicted coding
   targets fall inside the node the element sits in, and the nearest-TSS
   heuristic names the wrong gene for three in ten.

   **Per cell line (2026-09-12).** The tissue where the deletion moves the
   gene most is rarely the cell a closure test runs in, so the scorer now
   keeps each gene's fold change on the K562, HepG2, GM12878 and IMR-90
   tracks (`predicted_by_cell` per element; `scripts/
   enhancer_targets_by_cell.py` back-fills a chromosome); the gene-level
   closure in ATTRIBUTION.md takes the cell's own magnitude from there.

 The composition budget
   (docs/ATTRIBUTION.md) found that the regulatory tier holds 8.1 Mb of
   constrained sequence genome-wide with no target named. So the deletion
   tool was pointed at constraint: `scripts/constrained_targets.py` ranks a
   chromosome's distal enhancers by their phyloP-constrained fraction
   (Zoonomia, 37 MB of bigWig ranges for chr21's 6,618 elements) and scores
   the 100 most constrained (all ≥ 20% constrained bases; 412 qualify)
   against the uniform sample of 200 (`constrained_targets_chr21.json`):

   | chr21 distal enhancers | constrained top 100 | uniform 200 |
   |---|---|---|
   | some gene moves by ≥ 0.1 log2 | 73% | 63.5% |
   | strong effect (≥ 0.3) among those named | 25 of 73 (34%) | 40 of 127 (31%) |
   | silencer-like | 28 of 73 | 43 of 127 |
   | coding target = nearest TSS in the node | 60.9% | 67.8% |
   | coding target inside the node | 78.1% | 87.4% |

   Constraint predicts function: a constrained enhancer names a gene more
   often and the strongest effects are there (a 46%-constrained element
   moving the testis lncRNA LINC00945 by −1.2 log2; KCNE1, SIM2). It also
   reaches further: the constrained elements land beyond the CTCF-only node
   twice as often as the uniform ones (22% against 13%), which reads as the
   long-range, conserved enhancers being exactly the ones the proxy
   boundary cuts off. That is the attribution the budget asked for, one
   chromosome in: the constrained regulatory tier is where the targets are
   worth naming, and the tool names them for three in four.

   **The genome (2026-09-12).** The same script on the other 23 chromosomes,
   nothing tuned (the `constrained_targets_genome_wide` job; 100 elements per
   chromosome, chrY has only 5 that qualify; 59,415 of 513,972 distal
   enhancers pass the 20% threshold), folded by
   `scripts/constrained_targets_genome_wide.py` into
   `constrained_targets_genome_wide.json`:

   | distal enhancers, all chromosomes | constrained (2,305) | uniform (4,800) |
   |---|---|---|
   | some gene moves by ≥ 0.1 log2 | 75.7% | 62.4% |
   | strong effect (≥ 0.3) among those named | 556 of 1,744 (31.9%) | 735 of 2,995 (24.5%) |
   | silencer-like | 641 of 1,744 (36.8%) | 1,157 of 2,995 (38.6%) |
   | coding target = nearest TSS in the node | 71.9% | 70.7% |
   | coding target inside the node | 89.8% | 86.0% |

   The first finding holds everywhere: constraint predicts function on 22 of
   24 chromosomes (chr9 and chr15 are the exceptions, at 62% and 61% against
   66% and 68%), from chr1's 65% to chr19's 92%, and the strong effects are a
   third of the named ones instead of a quarter (C1QTNF7-AS1 −4.0 log2 from
   a 76%-constrained element; MEF2C-AS2 −3.2; the homeobox gene ARX −3.0 in
   Purkinje cells from a 92%-constrained element on chrX). The second one
   does not: chr21's constrained elements reached beyond the CTCF-only node
   twice as often as uniform ones, but over the genome the constrained
   targets sit inside the node as often as anyone's (89.8% against 86.0%),
   and name the nearest TSS at the same rate. The proxy boundary is not what
   the conserved enhancers are escaping; chr21 was one chromosome's noise.
   What stands is the attribution itself: of the 8.1 Mb of constrained
   regulatory sequence the budget could not name a target for, the most
   constrained elements name one three times in four, and the node they sit
   in predicts which gene nine times in ten.

   **Against measured enhancers (2026-09-12).** VISTA (LBNL) has tested about
   2,400 human sequences for enhancer activity in e11.5 mouse embryos, each
   positive in named tissues or negative; it is the ground truth the parser
   never had. `attribution/vista.py` streams the loci table (120 kB, hg38)
   and `scripts/vista_score.py --chrom C` reads every element blind: is an
   ENCODE enhancer-like element there, which node and nearest coding TSS,
   how constrained, and what the deletion of the whole element does in the
   model (`vista_<chrom>.json`; the `vista_genome_wide` job ran the 23
   chromosomes in five hours, 2,223 non-overlapping elements of 2,442):

   | VISTA elements, genome | positive (1,133) | negative (1,090) | chromosomes leaning this way |
   |---|---|---|---|
   | an ENCODE enhancer-like element sits there | 84.5% | 66.1% | 22 of 23 |
   | constrained (≥ 20% phyloP bases) | 65.0% | 57.7% | 18 of 23 |
   | deletion moves some gene by ≥ 0.1 log2 | 72.7% | 64.9% | 17 of 23 |
   | strong effect (≥ 0.3) among those named | 422 of 806 (52%) | 271 of 695 (39%) | |
   | coding target = nearest TSS in the node | 40.0% | 34.2% | |
   | predicted tissue falls in the measured one | 195 of 286 judged (68%; chance 36%) | no tissue | |

   Every layer leans the right way over the genome, and by different
   amounts. The registry separates best: a VISTA positive sits on an ENCODE
   enhancer-like element five times in six, a negative two times in three.
   Constraint separates least, as it must: VISTA chose every element for
   conservation, so both sides are constrained twice as often as the
   genome's regulatory tier. The deletion model names a gene for positives
   more often than for negatives (the two-chromosome reading of 75% against
   91% was fifty elements' noise; over 2,200 it is 73% against 65%), and the
   difference is in the strong effects: half the positives that name a gene
   move it by 0.3 log2 or more, against two fifths of the negatives. What
   the assay confirms most is the tissue: where the model's strongest track
   is a brain, heart, muscle or liver tissue, it is the tissue the embryo
   stained in two cases of three, against one in three if the labels were
   shuffled (hs1304 on chr21, forebrain to neural tube: −2.7 log2 in neural
   cells; hs2543 on chr22, brain: +2.6 in cerebellum; ARX's enhancer on
   chrX). The caveat stands: a VISTA negative is conserved sequence next to
   a gene that did not drive expression on one embryonic day, and the model
   reads adult and cell-line tracks, so "some gene moves" is a target
   finder with a bias towards real enhancers, not an activity assay; the
   activity assay is lentiMPRA, below.

   **Against measured targets (2026-09-12).** GTEx's cis-eQTLs are the
   measured answer to "which gene": a variant whose alleles go with a gene's
   expression across hundreds of donors, per tissue. `attribution/eqtl.py`
   streams GTEx v8's single-tissue archive (1.56 GB, 49 tissues, 71.5
   million significant pairs) member by member over HTTP ranges in 135 s,
   keeps the 295,139 pairs that fall inside or within 500 bp of an element
   GenomeOS has a prediction for, and discards the rest (34 MB under
   data/knowledge/gtex, local). `scripts/eqtl_targets.py` then asks, for
   every element with an eQTL, whether the gene each attribution names is
   one of the element's eGenes (`eqtl_targets.json`, job `eqtl_targets`):

   | elements | with an eQTL | eGenes each | deletion target is an eGene | coding genes only | nearest coding TSS is an eGene | coding model vs node when they disagree |
   |---|---|---|---|---|---|---|
   | uniform samples (4,800) | 3,725 (78%) | 4.8 | 52.9% (2,372) | 71.7% (1,851) | 66.2% (3,725) | 364 vs 309 of 554 |
   | most constrained (2,305) | 1,394 (60%) | 4.2 | 44.8% (1,054) | 56.0% (921) | 54.9% (1,394) | 148 vs 154 of 292 |

   Read with the model's universe in mind. Asked for any gene, the deletion
   model names an eGene for half the elements; the node's nearest coding TSS
   for two thirds. Asked for a coding gene, the same comparison the enhancer
   deletions were held to, the model is right for 71.7% against 66.2% on the
   uniform sample, and where the two disagree it wins 364 to 309; on the
   constrained tier the two are level. The gap between "any gene" and
   "coding" is the model naming lncRNAs and pseudogenes that GTEx has little
   power to test, not evidence that those calls are wrong, and the 4.8
   eGenes per element mean "is an eGene" is a loose bar that both callers
   clear more often than a random neighbour would. The tissue is harder: for
   the 101 elements whose predicted track is a GTEx tissue, that tissue is
   one where the eQTL for the same gene is significant in 40%. Constrained
   elements carry eQTLs less often (60% against 78%): where selection keeps
   the sequence, the common variants an eQTL needs are fewer, which is why
   the deletion model was pointed there in the first place.
   **Against consequence (2026-09-12).** The last ground truth is the
   weakest and the one a person cares about: does anything that matters to
   a human land on the element? `attribution/gwas.py` streams the GWAS
   Catalog's association table once (69 MB zipped, 982k rows with a GRCh38
   position) and keeps the lead variants within 1 kb of an element, and of
   the same element shifted 100 kb along the chromosome as the chance
   level; `scripts/consequence_targets.py` (job `consequence_targets`)
   reads them over the three element sets, with ClinVar's pathogenic
   variants whose molecular consequence is not a coding change as the
   second, sparser, witness (`consequence_targets.json`):

   | elements | hold a lead variant | shifted control | enrichment | named target / none | strong effect |
   |---|---|---|---|---|---|
   | uniform samples (4,800) | 33.2% | 29.4% | 1.13× | 33.3% / 32.9% | 35.0% |
   | most constrained (2,305) | 37.4% | 30.8% | 1.21× | 38.0% / 35.8% | 39.9% |
   | VISTA elements (2,223) | 47.5% | 42.9% | 1.11× | 50.6% / 41.0% | 55.8% |

   Lead variants land on the attributed elements barely more often than on
   the same stretches 100 kb away: a fifth more on the constrained tier, a
   tenth more elsewhere. Within the sets the attribution moves the rate a
   little in the right direction: VISTA elements with a named target hold
   a lead variant for 51% against 41% without, the strong
   effects 56%, VISTA positives 49% against negatives
   46%; and constrained VISTA elements hold fewer (42% against
   56%), the eQTL lesson again: where selection keeps the sequence the
   common variants a GWAS needs are missing. The traits are the polygenic
   ones (height, body mass index, blood pressure, type 2 diabetes,
   insomnia), which is what a random kilobase of the regulatory genome
   carries. The catalog's mapped gene, its own nearest-gene call, agrees
   with the node's nearest TSS for 63% of the uniform elements and with
   the model's target for 39%, which measures the catalog's heuristic, not
   the truth. ClinVar is the sparse witness: 20, 30 and 24 elements of the
   three sets hold a pathogenic variant with a non-coding consequence, and
   where they do, ClinVar's gene is the model's target for 11, 20 and 19
   of them and the node's for 13, 17 and 12. Consequence is the level
   where the regulatory attribution is least testable today; the four
   measured layers above are where it is tested.

3. **Splicing, at the resolution our grammar lacks** (built, feature c). Our
   learned donor and acceptor matrices reach about 90% recall at seven false
   hits per kilobase, which is why segments are parsed by grammar and never
   asserted (D19). `genomeos segments --chrom C --predicted-sites` feeds the
   parser AlphaGenome's splice-site tracks instead (four probabilities per
   base, donor and acceptor per strand; 1 Mb per request, 45 requests and
   33 s for chr21, cached under `data/knowledge/alphagenome/splice`). The
   tracks were calibrated first: they mark the last exonic base (donor) and
   the first exonic base (acceptor) on both strands, mean probability 0.91 to
   1.00 at canonical boundaries. With nothing else changed, exact CDS
   segments on chr21 go from 4.7% / 5.2% to 40.1% / 46.2% and exact splice
   sites from 7.5% / 6.9% to 51.6% / 49.2% (sensitivity / precision;
   docs/SEQUENCE-GRAMMAR.md has the table and what did not move). Everything
   stays `predicted`; without the key the grammar runs as before.
4. **Domains from a predicted contact map** (built, 2026-09-12; a weak
   agreement, recorded). `genomeos domains --chrom chr21 --contact-map`
   asks for the contact map of every 1 Mb window (45 requests, 28 4DN
   Micro-C and Hi-C cell types at 2 kb), computes the insulation score of
   each bin (mean contact across the diagonal over 20 kb, averaged over the
   cell types) and takes its local minima as predicted boundaries
   (`genomeos/predict/contact_maps.py`, `data/results/domains_chr21_contact.json`).
   Against the 227 CTCF-only boundaries: 83 sit within 20 kb of one of the
   353 predicted minima (37%), and boundaries placed at random get 28%; 24%
   of predicted minima have an inferred boundary; the median distance from
   an inferred boundary to the nearest minimum is 34 kb. The depth threshold
   was swept (0.02 to 0.3; stricter thresholds give fewer minima and less
   support, never more above chance). So the two kinds of evidence agree
   only a little above chance: the CTCF-only proxy draws boundaries where
   the predicted map sees no insulation dip, and the map's dips mostly lack
   a CTCF-only element. The node's confidence stays at 0.4, and the honest
   conclusion is that node *content* is what the project has evidence for
   (the enhancer and mouse results) while node *edges* remain a proxy that
   neither this prediction nor, without Hi-C, a measurement has confirmed.
5. **The reader, predicted where it is not measured**. Reader v1 reads ENCODE
   DNase peaks per biosample, so it only knows the cell types ENCODE assayed.
   Predicted DNase and ATAC extend "which nodes are open in this cell type" to
   tissues with no experiment, clearly marked predicted.
6. **Which bases in a block matter**. In-silico mutagenesis over a few hundred
   base pairs of an unclassified block gives a contribution profile: a direct,
   if expensive, answer to "what is this sequence for". Roughly three requests
   per base, so it is a per-block tool, never a scan.
7. **The twin, with its limitation stated** (built, feature d). The twin
   question is "what do this person's variants do to this gene", not "what
   does this variant do". `genomeos individual predict --name ME --gene APP
   --chrom chr21` takes the elements the regulation layer says reach the gene
   (promoter elements first, then the enhancers AlphaGenome itself tied to the
   gene, then the near ones, then the rest of the node), collects the person's
   variants inside them (SNVs and indels up to 50 bp), scores each with the
   gene-level RNA-seq scorer and sums the effects per haplotype where the
   genotypes are phased (`genomeos/predict/individual_effects.py`). Per-variant
   answers are cached under `data/knowledge/alphagenome/variants`; per-person
   results are never written under `data/results`, so nothing about a named
   person's genome can be committed by accident. Checked on HG002's chr21 rows
   imported as `demo`, gene APP: 39 variants inside 99 elements, scored in 41 s;
   four move APP by ≥ 0.1 log2, the strongest a 6 bp insertion at
   chr21:26,137,053 in enhancer EH38E3453088 (−0.30 in naive CD4 T cells, 23
   tracks), the same element the deletion job (feature b) had already tied to
   APP. The model reads one unphased sequence, so a haplotype sum is a sum of
   single-variant predictions and not a prediction on the whole haplotype; that
   caveat travels with every result. All four features are now built; the
   Twin tab runs this from the "Your genome" card.

## How it fits the architecture

- **Evidence.** Every output becomes a BioIR `Rule` with evidence kind
  `predicted`, the model name and version as source, and confidence derived
  from effect magnitude, capped at 0.7 (`genomeos/predict/alphagenome_adapter.py`).
  A prediction never becomes a measurement, and `genomeos verify` style checks
  stay the arbiter where a measurement exists.
- **Never per tick.** As with every external source (D23), knowledge is
  imported and cached, never called from inside a simulation. Responses cache
  under `data/knowledge/alphagenome/` keyed by variant, scorer and organism;
  what gets committed is the distilled summary in `data/results/`, per the
  stream-distil-discard rule (D7).
- **Quota-shaped jobs.** One element or one variant per request means the unit
  of work is a chromosome's regulatory elements, run as a background job that
  resumes, not a genome-wide sweep. The Progress tab already shows this kind
  of job.
- **Optional, always.** The core never depends on it; without a key the
  features are listed and disabled (D34).

## Open weights

Since 28 January 2026 the model weights are published (all-folds and
individual folds) with JAX research code, under model terms that remain
non-commercial; commercial use goes through self-deployment on Google Cloud.
That matters for GenomeOS later: a self-hosted model removes the quota
ceiling and makes uses 2, 5 and 6 genome-wide instead of per-chromosome. It
does not change the evidence rules, and it does not change the licence
position in [ACKNOWLEDGEMENTS.md](../ACKNOWLEDGEMENTS.md): enabling
AlphaGenome binds the user to non-commercial terms whatever GenomeOS's own
licence says.
