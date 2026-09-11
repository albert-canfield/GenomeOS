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
   chromosomes run through the same job by name; the summary card on the
   Progress tab shows one row per chromosome scored.
3. **Splicing, at the resolution our grammar lacks** (feature c). Our learned
   donor and acceptor matrices reach about 90% recall at seven false hits per
   kilobase, which is why segments are parsed by grammar and never asserted
   (D19). AlphaGenome predicts splice sites, their usage and the junctions at
   1 bp, so the segment parser can carry a second, independent opinion, and a
   variant can be reported as splice-disrupting with a predicted junction
   change rather than a motif score.
4. **Domains from a predicted contact map**. `genomeos domains` infers a node
   between CTCF-only elements and labels it inferred. A predicted contact map
   at 2 kb over a 1 Mb window is a different kind of evidence for the same
   boundary, and where the two agree the domain earns a higher confidence.
5. **The reader, predicted where it is not measured**. Reader v1 reads ENCODE
   DNase peaks per biosample, so it only knows the cell types ENCODE assayed.
   Predicted DNase and ATAC extend "which nodes are open in this cell type" to
   tissues with no experiment, clearly marked predicted.
6. **Which bases in a block matter**. In-silico mutagenesis over a few hundred
   base pairs of an unclassified block gives a contribution profile: a direct,
   if expensive, answer to "what is this sequence for". Roughly three requests
   per base, so it is a per-block tool, never a scan.
7. **The twin, with its limitation stated**. An individual's non-coding
   variants can be scored for predicted expression change. The model reads one
   unphased sequence, so a diploid twin gets two runs and no interaction
   between haplotypes; that caveat travels with the result.

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
