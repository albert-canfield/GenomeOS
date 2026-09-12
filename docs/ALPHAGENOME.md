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
