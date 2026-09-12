# The grammar by comparison: reading the genome as code across species and people

Two conversations on 2026-09-12 proposed a method rather than a metaphor:
treat the genome as a 3 GB undocumented executable of which millions of
slightly different versions exist (people, chimpanzees, mice, fish), and
recover its grammar by comparison:

```
compare many genomes → find what is conserved and what varies → infer the grammar
```

with a five-level grammar (alphabet, lexical, structural, regulatory, module,
organism), a "decompiled view" per locus, and a compiler pipeline (align, lex,
parse, link, infer semantics, mine libraries) whose output is BioIR. This
document holds that proposal against what GenomeOS has already built and
measured, says where the proposal is right, where the data says it is wrong,
and what to build next in order. It is the design document of ROADMAP.md
area J. The reverse-engineering map by consideration (timers, blueprint,
parts, systems) stays in GENOME-AS-CODE.md; the delimiters in
SEQUENCE-GRAMMAR.md; the attribution of the non-coding space in
ATTRIBUTION.md.

## 1. The proposal against the code

Each idea of the conversation, what GenomeOS has for it, and the measured
state. "Done" means a result file exists and a document records the number.

| Idea in the conversation | GenomeOS today | State |
|---|---|---|
| Alphabet A C G T; a compiler does not read characters | the four-letter alphabet as a hypothesis-space constraint (ATTRIBUTION.md, "A note on the alphabet") | done |
| Tokens: start and stop codons, splice sites, polyA, promoter motifs | `genomeos signals`: donor, acceptor, Kozak matrices learned from chr21; stops, polyA, CpG islands measured (SEQUENCE-GRAMMAR.md §1) | done, with one correction: they are probabilities, not tokens (§3) |
| Hierarchy base → motif → element → gene → network → pathway → cell program → module → organism | BioIR types: `signal`/segment, `element`, `gene`/`transcript`/`protein`, `rule`, library, `cell_type`, `organism`; the Blocks tab shows nodes above genes | done as types; the levels are not yet linked bottom to top for one locus |
| Structural parse: genes, exons, introns from signals | `genomeos segments`: twelve scored runs on chr21 and chr22; 92.7% gene precision once RNA over the exons is asked for | done and closed |
| Regulatory DNA as operators (`IF TF_A AND TF_B AND open THEN express X`) | `rule` with `strength`, `threshold`, `hill`, `when:`; enhancer rules compiled from AlphaGenome deletions (`genomeos budget --bio`) | partial: rules exist, but their conditions come from a predictor, not from motifs. The motif scan (JASPAR) is planned in SEQUENCE-GRAMMAR.md §5 and not built |
| The promoter is the dependency list (`requires:`) | designed in SEQUENCE-GRAMMAR.md §2; `cell_type ... expresses:` is the hand-written stand-in | designed, not built |
| Cross-species conservation marks the keywords | Zoonomia phyloP over 241 mammals per base for every UNKNOWN block (1,009 Mb, 24 chromosomes, 2,252 MB of ranges, never stored) plus the 100-vertebrate elements; five tiers per block | done genome-wide |
| Within-species variation marks the tolerant positions | `genomeos syntax GENE` (in progress, `attribution/syntax.py`): per base of one gene, constrained across mammals against variable among the three imported people (HG002 trio), rs12913832 as the named value; gnomAD allele frequency per variant through VEP in the cancer module | partial: three genomes at one gene; the population axis (76,156 genomes, every block) is §4 step 1 |
| Cases A to D: conserved everywhere / variable everywhere / human-specific / clade module | needs both axes above plus a human-lineage layer | **missing**, computable from what exists plus step 1 |
| Orthologues as inherited library functions, `Mammal::GeneX` | MGI curated mouse–human orthology (24,584 symbol pairs) in the node comparison; C. elegans orthologues not used | partial: one species pair; no clade of origin per gene |
| Gene order carries information (synteny) | mouse chr19 and chr11 nodes against human nodes: 379 nodes, 93 to 94% in the same human neighbourhood (NODES-READER-WRITER.md) | done for two chromosomes; this *is* synteny at the node level |
| Gene duplication as copy and paste (paralogues) | `similar_to` between the sixty largest UNKNOWN blocks (shared 20-mers) | **missing** as curated data: no paralogues, no segmental duplications |
| Protein domains as the closest thing to reusable functions | InterPro domains on 99.0% of 19,478 compiled proteins; `has domain` edges in the knowledge graph | done as annotation; no analysis of *reuse* (which domains combine, in which libraries) |
| Pathways as higher-order libraries (`WNT.lib`) | BioLib: 45 libraries with membership computed from GO and Reactome (114 KB shipped) | done; libraries are curated, not mined, and carry no evolutionary depth |
| Context as function arguments | `when:` on every rule; the reader per cell type (eleven ENCODE DNase readers, genome-wide); GTEx per tissue and per isoform | done |
| Variation as natural experiments | HG002 and the trio, ClinVar screen, eQTL targets, GWAS enrichment, the 31 reference truncations | done at gene and element level |
| Covariance across species reveals cross-references (coevolution) | nothing | **missing** |
| Perturbation reveals operators | AlphaGenome deletions (predicted, 4,800 uniform and 2,305 constrained elements); VISTA (2,223 measured), lentiMPRA (K562, HepG2, WTC11) as ground truth; closure at the gene level in four cell lines | done as predicted plus measured; CRISPR screens (DepMap, IMPC) planned in ATTRIBUTION.md |
| 3D contact as a pointer | CTCF nodes; 90.2% of predicted enhancer targets inside the node; predicted contact maps; measured Hi-C built and waiting for a 4DN key | done as inference and prediction |
| Build a graph, not only a browser | knowledge graph, 41,982 nodes and 330,018 edges (proteins, pathways, tissues, domains, writers) | done for molecules; **no non-coding elements and no evolutionary edges** (orthologue-of, paralogue-of) |
| Separate known, inferred and unknown, with provenance | `Evidence.kind` in {experimental, curated, predicted, inferred, none}, a confidence on every entity and rule, `UNKNOWN` legal, predicted capped at 0.7 | done; D4 |
| Mine packages from the graph (community detection) | nothing | **missing** |
| The compiler pipeline: aligner, lexer, parser, linker, semantics, library miner, IR, VM | signals and segments (lexer, parser); regulation and enhancer targets (linker); rules (semantics); lib membership (libraries, from curation); BioIR; BioVM | no **aligner** layer of our own (we read scores others computed from alignments), no library **mining** |
| The decompiled view per locus | `genomeos budget --bio` compiles chr21's non-coding space to BioLang (747 entities, 155 rules); `genomeos protein X --bio`; the Blocks tab | partial: two outputs, not one view of one locus across every layer with orthologues |
| The `BioCodeSignature` per block | the budget record per block: class, repeat coverage, phyloP mean, maximum and constrained fraction, 100-way elements, tier, confidence | partial: half the fields |
| A species ladder rather than human alone | C. elegans grown from one cell; mouse chr19 and chr11; human | three rungs of ten |
| A first developmental locus decompiled end to end | not chosen | **missing** (§4 step 7) |

Reading the table: the conversation independently arrives at the design
GenomeOS has been executing for three days. Of its twenty-six ideas, fifteen
are built and measured, five are partial, six are missing. The six missing
ones share a shape: they are all **comparisons GenomeOS has not yet made**,
across people (variation per base), across the tree (clade of origin,
covariance), and within the genome (duplication). That is the enrichment.

## 2. What the conversation adds that GenomeOS should adopt

1. **The second axis of constraint.** GenomeOS reads constraint across 241
   mammals and calls a block a candidate for function when it passes. The
   conversation's cases A to D need the second axis: constraint *within
   humans*. Albert's framing of 2026-09-12 is the same request in other words:
   the bases that never vary across people are syntax, the bases that vary
   without consequence are value slots, and the two axes together separate
   the two. gnomAD's Gnocchi score (a Z score of observed against expected
   variation per kilobase from 76,156 genomes; Z ≥ 2.18 is the top decile of
   constrained non-coding sequence, Z ≥ 4.0 the top percentile) is published
   as a bigWig at UCSC, and our range reader reads it: three chr21 windows
   cost nine requests and 40 KB in 4.5 s on 2026-09-12. One of the three (chr21:
   40.00–40.02 Mb) is constrained in humans (maximum Z 4.3, 2 kb above the
   threshold); the other two are not.

   | mammals (phyloP) | humans (Gnocchi) | reading | the conversation's case |
   |---|---|---|---|
   | constrained | constrained | syntax: fundamental, deep | A |
   | unconstrained | unconstrained | tolerant, or a value slot | B |
   | unconstrained | constrained | recent: primate or human function, or recent selection | C |
   | constrained | unconstrained | lost or relaxed in the human line; or a slot whose *value* is what varies inside a constrained frame | D, and the eye-colour test |

   The last row is the interesting one for Albert's eye-colour example: the
   OCA2 enhancer in HERC2 intron 86 that carries rs12913832 (hg38
   chr15:28,120,472) is conserved across mammals *and* polymorphic in humans.
   A constrained frame with a variable position inside it is exactly a
   value slot in a syntactic construct, and the two-axis table finds such
   positions genome-wide instead of one at a time.

2. **Phylogenetic depth as the library boundary.** A library's *members* we
   compute from GO and Reactome; its *age* we do not. The conversation's
   `Life.Core / Eukaryote.Core / Animal.Core / Vertebrate.Core / Mammal.Core /
   Primate / Human-specific` is a phylostratigraphy: for each gene, the
   deepest clade that has an orthologue. Ensembl's homology REST endpoint
   answers it per gene with a `target_taxon` (OCA2 has 92 mammalian
   orthologues in one call, checked 2026-09-12); eight taxa (Mammalia,
   Amniota, Vertebrata, Chordata, Metazoa, Opisthokonta, Eukaryota, all) give
   each gene an origin. Aggregated per library, this is a test the catalogue
   has never had: `core.translation` should be universal, `blueprint.axes`
   metazoan, `systems.immune` (adaptive) vertebrate, and a library whose
   members span every stratum is not one library.

3. **Duplication as copy and paste, read from curated data.** Paralogues
   (Ensembl homology, `type=paralogues`) and segmental duplications (UCSC
   `genomicSuperDups`, and the HPRC duplications relative to hg38) replace
   the 20-mer `similar_to` heuristic and give the graph its `paralogue_of`
   edges. This is also the organising evidence ATTRIBUTION.md lists for the
   compiled blocks and has not yet built.

4. **Covariance as cross-reference.** Two libraries whose members are gained
   and lost together across the tree depend on each other. With the origin
   per gene from item 2 (or OrthoDB's hierarchical groups), phylogenetic
   profiling across the 45 libraries is a small computation and the first
   evidence of *dependency between libraries* that does not come from a
   curated pathway.

5. **The motif scan as the source of operators.** JASPAR 2026's CORE
   vertebrate collection has 1,019 profiles and an API that returns them
   (checked 2026-09-12). Scanning promoters and the constrained elements
   gives `requires:` lists with evidence, and, across the genes of one
   library, the recurring *combinations* (GATA plus T-box in the heart
   field, ETS plus RUNX in blood) are the boolean operators the conversation
   describes. SEQUENCE-GRAMMAR.md already plans the scan; the new point is
   to look for recurring operator patterns across genes, not motifs one gene
   at a time, which is what Albert asked for on 2026-09-12.

6. **The decompiled locus as the product.** One command and one view that,
   for a locus, prints the BioLang block with every layer present and the
   evidence table under it: the parse, the elements with targets, the node,
   the reader per cell, constraint on both axes, the orthologues with their
   depth, the paralogues, the domains, the rules, and what is `UNKNOWN`.
   Everything in that list exists as a separate command today.

7. **The species ladder, extended where ground truth exists.** The
   conversation's ladder (E. coli to human) is right; the order should
   follow where GenomeOS can *test*: the worm (grown, lineage as ground
   truth), mouse (nodes compared), then Drosophila (genome-wide STARR-seq
   for enhancer grammar), fugu (compact vertebrate) and zebrafish, as
   ATTRIBUTION.md already lists.

## 3. What the data says the conversation gets wrong

- **Tokens are scores, not symbols.** A learned donor matrix that keeps 90%
  of true donors fires seven times per kilobase of random sequence; the
  median true donor scores 0.57 of the matrix maximum (LESSONS.md). `ATG`
  is a start codon in 220 of 221 chr21 genes, but the chromosome has a
  start-codon-shaped triplet every 64 bases. The lexer is therefore a
  probabilistic gene finder, and a "token" only becomes one inside a
  structure with context (D19). The conversation's own caveat ("a token is
  only a candidate until context confirms its function") is the design;
  the measured cost of ignoring it is 556 candidates for 221 genes.
- **Conserved is "matters", not "reserved code".** Constraint predicts
  function (the most constrained enhancers name a gene in 75.7% of
  deletions against 62.4% for a uniform sample), but it separates measured
  VISTA positives from negatives only weakly (65.0% against 57.7%
  constrained), because VISTA chose its elements for conservation and a
  negative is conserved sequence that did not drive expression on one
  embryonic day. Constraint says a block is under selection; what it does,
  and when, needs the reader and a perturbation.
- **Variable does not mean tolerant.** A position that varies across people
  may be a value slot (eye colour), a recessive or late-onset effect, or
  under recent selection; and Gnocchi's resolution is one kilobase. The
  within-human axis must enter as evidence with a confidence, like every
  other layer, not as a verdict of neutrality (D4, D42).
- **The pointer is mostly local.** The conversation stresses that an
  enhancer 500 kb away can touch its promoter. It can; but for 90.2% of the
  predicted targets on 24 chromosomes the gene sits inside the element's own
  CTCF node, and 71.2% are exactly the nearest TSS. Long-range references
  exist and the model is agnostic about them; the node is the address
  space, and the pointer resolves inside it two times in three.
- **"Function = gene" is the wrong level.** GenomeOS's model is Gene →
  Transcript(s) → Protein isoform(s) → modified state (D23), and one gene
  makes different proteins in different tissues (APP-202 in the cerebellum,
  APP-201 elsewhere). The reusable unit the conversation wants is the
  protein domain, the rule, or the library, not the gene.
- **The knowledge graph is not enough by itself.** Community detection on a
  graph of curated pathways finds the curated pathways. Mined "packages"
  are credible only when several independent layers agree (co-expression,
  co-constraint, co-gain across species, co-perturbation), which is the
  conversation's own triangulation principle and D4's reason to exist.

## 4. What to build, in order

Each step is a chromosome-21 run first, a genome-wide job second, a number
in a table, and evidence in BioIR. Owner: area J in ROADMAP.md.

1. **The within-human axis at population scale.** `attribution/syntax.py`
   (a peer lane, 2026-09-12) reads the axis per base of one gene from the
   three people imported locally; three genomes see common variants only.
   `attribution/variation.py` adds the population: Gnocchi Z per
   UNKNOWN block and per attributed element through the existing bigWig
   reader (1 kb resolution, so a chromosome costs a few megabytes); the
   two-axis table of §2 per block, per tier and per element class; the
   compiled `region` and `element` blocks gain a `variation` evidence line.
   Tests that must hold: mammal-constrained blocks are depleted of human
   variation more than neutral blocks; VISTA positives more than negatives;
   the HERC2 element around rs12913832 lands in row D (constrained frame,
   variable slot), agreeing with what `genomeos syntax OCA2` reads from the
   trio at base resolution. Then genome-wide, as the budget job ran.
2. **Origin per gene, age per library.** `knowledge/homology.py`: Ensembl
   homology per gene at eight target taxa, cached and distilled to one
   stratum per gene (19,478 requests once, or the Compara species-set
   download); `genomeos libs --origin` prints the stratum distribution per
   library and flags libraries that span the tree. Test: ribosomal proteins
   universal, HOX metazoan, MHC and RAG vertebrate. The gene gains an
   `origin:` attribute in BioIR and the graph gains `orthologue_of` edges.
3. **Duplication.** Paralogues from the same endpoint; segmental
   duplications from UCSC `genomicSuperDups` per chromosome through the
   track API (as RepeatMasker is read); `paralogue_of` edges; the block's
   `similar_to` replaced by curated duplication with the heuristic kept as
   fallback (D20's order).
4. **The motif scan and the operators.** JASPAR CORE vertebrates, 1,019
   profiles, cached once; scan promoters (TSS ± 1 kb) and constrained
   elements; `requires:` derived with `curated` evidence for the matrix and
   `predicted` for the hit; across each library's genes, the recurring motif
   pairs and triples with their enrichment against shuffled promoters.
   Test: the heart-field and blood combinations above.
5. **Phylogenetic profiling of libraries.** From step 2's strata, the
   co-gain and co-loss of libraries across the tree; a dependency edge
   between libraries where the profiles agree beyond chance; reported next to
   Reactome's curated hierarchy as `inferred`.
6. **The decompiled locus.** `genomeos decompile chrN:a-b` (or a gene
   symbol): the BioLang block for the locus across every layer, the evidence
   table per layer with its kind and confidence, and the `UNKNOWN` residue
   listed last; a Blocks-tab card. No new inference, only assembly.
7. **One developmental locus end to end.** The conversation's proposed
   experiment. Two candidates, both with every layer available: the OCA2 and
   HERC2 pair (Albert's syntax-versus-value test; VISTA, eQTL, Gnocchi and
   Zoonomia all cover it) and the SHH limb enhancer ZRS (a 1 Mb pointer,
   conserved to fish, with measured mouse phenotypes). Run the full pipeline
   on human, mouse and, where an assembly is on UCSC, chicken and zebrafish,
   and report whether the module reconstructs.

## 5. Sources verified reachable on 2026-09-12

| Source | What | Access | Cost seen |
|---|---|---|---|
| gnomAD Gnocchi | within-human constraint, Z per kb, 76,156 genomes | `https://hgdownload.soe.ucsc.edu/gbdb/hg38/gnomAD/mutConstraint/mutConstraint.bw`, range requests through `attribution/bigwig.py` | 9 requests, 40 KB, 4.5 s for three windows |
| UCSC 470 mammals | phyloP and phastCons over 470 mammals (Hiller lab alignment) | `goldenPath/hg38/phyloP470way/hg38.phyloP470way.bw`, `phastCons470way/` | not yet read; same reader |
| Zoonomia 241 mammals | phyloP per base, the alignment as bigMaf | `goldenPath/hg38/cactus241way/` | in use (D42) |
| Ensembl Compara REST | orthologues and paralogues per gene, `target_taxon` filter, gene trees | `rest.ensembl.org/homology/symbol/human/OCA2?type=orthologues;target_taxon=40674` | one call, 92 mammal orthologues |
| JASPAR 2026 | 1,019 CORE vertebrate profiles as PFM | `jaspar.elixir.no/api/v1/matrix/?collection=CORE&tax_group=vertebrates` | one call per matrix or one bulk download |
| UCSC `genomicSuperDups` | segmental duplications > 1 kb | track API, as RepeatMasker is read | not yet read |
| HPRC 90-way and v2.1 SVs | 90 human assemblies aligned; duplications, inversions, insertions against hg38 from 233 assemblies | UCSC `hprc90way`, `hprcArr*`, `hprc2v21Sv` tracks | not yet read |
| OrthoDB v12 | hierarchical orthologous groups per clade | REST and SPARQL | alternative to Ensembl for step 5 |

Not found: a phyloP bigWig for the 447-way (primate-expanded) alignment
under UCSC's `cactus447way/`, which holds only the alignment and the trees.

## 6. Step 1 built: the human axis beside the mammalian one (2026-09-12)

`genomeos variation --chrom C` (`attribution/variation.py`; results
`variation_chr21`, `variation_chr15`, `variation_vista_genome_wide`) reads
gnomAD's Gnocchi score over every UNKNOWN block of the budget and every
element the target runs scored, next to the Zoonomia constraint the budget
already holds, and files each in one of the four cases of §2 with a
confidence capped at 0.6. Two controls travel with every chromosome: the
canonical coding segments and introns of its genes (which must read as
constrained among people if the axis is real) and VISTA's measured elements
(positives should lean further than negatives). Gnocchi is one value per
kilobase, so an interval shorter than 1 kb counts its whole kilobase and
every share below is a share of kilobases touched. Chromosome 21 cost 4.2 MB
of ranges and 156 s (most of it phyloP for the uniform elements the
constrained run had not scored), chromosome 15 4.5 MB and 222 s.

**The controls hold, and the tiers separate less than the exons do.** Share
of kilobases with Z ≥ 2.18:

| set | chr21 | chr15 |
|---|---|---|
| canonical coding segments | **40.3%** | **30.9%** |
| canonical introns | 25.9% | 17.6% |
| UNKNOWN, regulatory tier | 16.9% | 14.5% |
| UNKNOWN, neutral tier | 2.1% | 17.3% |
| UNKNOWN, fossil tier | 1.7% | 7.6% |
| UNKNOWN, constrained_unknown tier | **1.7%** | **3.2%** |

Coding sequence is the most constrained among people on both chromosomes,
introns next, which is the published shape and says the reader and the
threshold are right. Below that the human axis is noisy at block resolution:
chr21's neutral tier reads 2% and chr15's 17%, with 31 measured blocks. And
the most conserved intergenic blocks, the constrained_unknown tier that
area I named as the real unknown, are *not* depleted of human variation at
this resolution: 1.7% and 3.2%, no more than the neutral tier. Of chr15's 31
measured constrained_unknown blocks, 18 read `relaxed` (held across mammals,
variable among people) and 3 `syntax`. Either their function does not show
at one kilobase in 76,156 genomes, or what varies inside them is a value the
frame permits, which is exactly the question step 7 asks of one locus.

**On the elements the two axes add.** Share of elements whose deletion names
a coding gene in AlphaGenome, by case:

| case (mammals, people) | chr21 | chr15 |
|---|---|---|
| syntax (constrained, constrained) | **75.0%** (44) | **72.2%** (36) |
| relaxed (constrained, free) | 59.6% (47) | 27.3% (66) |
| recent (free, constrained) | 49.1% (55) | 54.5% (33) |
| tolerant (free, free) | 38.8% (85) | 46.1% (102) |
| unmeasured on one axis | 63 | 60 |

The elements constrained on both axes name a gene three times in four on
both chromosomes, above every other case; that is the first result in the
project where the human axis adds to the species axis rather than repeating
it. The order of the middle two rows swaps between chromosomes, so neither
axis alone is a stable second predictor at this sample size; the pair is.
Elements strong in humans (Z ≥ 4.0 somewhere) occur only in the two
human-constrained cases (14 and 15 on chr21, 12 and 6 on chr15), as the
definition requires.

**VISTA.** Share of kilobases constrained among people, positives against
negatives: chr21 32.4% against 25.0% (13 and 6 elements), chr15 18.8% against
19.5% (28 and 34), the genome 20.2% against 16.9% over 1,133 positives and
1,090 negatives on 23 chromosomes (32.9% against 26.9% touch a constrained
kilobase, 6.4% against 5.4% a strong one; 15 MB, 124 s). The lean is the
right way and about the size of the mammalian one (65% against 58% in
ATTRIBUTION.md), for the same reason: VISTA chose conserved sequence on both
sides.

**The known value.** rs12913832, the eye-colour position in HERC2's intron
86 (hg38 chr15:28,120,472), read by the module on both axes: the base is
constrained across mammals (phyloP 3.41, above the 2.27 bar); the kilobase
around it holds 7.7% constrained bases with a maximum of 8.59 inside the
registry element 366 bp upstream (EH38E1749549, 18% constrained), and that
kilobase is one Gnocchi does not score at all; the element 242 bp downstream
sits in a kilobase at Z 3.46. So the textbook value slot is a constrained
base in a constrained frame, and the population axis is silent exactly
there. The per-base form of the same question, the trio's variants against
phyloP at one gene, is `genomeos syntax --gene OCA2 --chrom chr15` in the
peer lane.

**Reading.** The axis is real (the exons say so), agrees with the species
axis where both are strong (the `syntax` elements), and on its own, at one
kilobase, is a weak and chromosome-dependent signal for blocks. The two
axes are reported side by side and never merged (D43); the case is
`inferred` evidence at 0.3 to 0.6. What it changes for area I: the
constrained_unknown blocks now carry a second fact, that people vary in
them, which narrows what they can be.

**The genome (24 chromosomes, 136.9 MB of ranges, 2026-09-12).** The ordering
two chromosomes suggested holds at scale. Share of kilobases constrained
among people: canonical coding segments 38.1%, introns 18.9%, the
regulatory tier 11.9%, fossil 3.0%, neutral 2.7%, constrained_unknown 2.4%,
structural 0.4%; per chromosome the coding share runs from 22% (chr4) to
58% (chr19), and chrY is unmeasured throughout because Gnocchi has no
kilobase there, which the result records as unmeasured rather than as
zero. Elements by case, share whose deletion names a coding gene:

| case (mammals, people) | elements | name a gene | strong in people |
|---|---|---|---|
| syntax (constrained, constrained) | 825 | **73.9%** | 315 |
| relaxed (constrained, free) | 1,337 | 56.0% | 0 |
| recent (free, constrained) | 782 | 50.0% | 123 |
| tolerant (free, free) | 2,235 | 45.3% | 0 |

Both axes beat either alone on every count; per chromosome the syntax
elements name a gene between 51% (chr1) and 87% (chr6). The
constrained_unknown tier splits into 77 `syntax`, 492 `relaxed`, 30
`recent` and 346 `tolerant` blocks of 948 measured: half of area I's real
unknown is held across mammals and variable among people, which is the
signature of the copies §7 finds there and, where it is not a copy, the
value-slot reading. VISTA positives 20.2% against negatives 16.9% over
2,223 elements, as reported above.

## 7. Step 3 built: duplication from curated data (2026-09-12)

`genomeos duplications --chrom C` (`genome/duplications.py`; results
`duplication_chr*` for all 24 chromosomes and `duplication_genome_wide`)
reads UCSC's `genomicSuperDups`, the curated list of segmental duplications
(pairs over 1 kb at 90% identity or more, with the partner locus and the
fraction of matching bases), through the track API the way RepeatMasker is
read, two seconds per chromosome, and lays it over the budget's UNKNOWN
blocks, the classifier's shared-k-mer pairs and the scored elements.

**The genome.** 64,216 pairs, 21,834 within their own chromosome, 166.9 Mb
duplicated (5.4% of the genome); by identity 5,281 young pairs (≥ 98%,
primate-era), 17,825 middle, 41,110 old. Duplicated fraction of the UNKNOWN
tiers, bp-weighted, with the share of blocks that are mostly (≥ 50%) a copy:

| tier | Mb | duplicated | blocks mostly duplicated |
|---|---|---|---|
| neutral | 73.3 | **10.3%** | 693 of 2,632 (26.3%) |
| fossil | 328.5 | 8.5% | 1,240 of 7,302 (17.0%) |
| constrained_unknown | 32.0 | 4.7% | **216 of 1,098 (19.7%)** |
| structural | 229.7 | 3.6% | 56 of 238 |
| regulatory | 345.2 | 2.8% | 689 of 15,536 (4.4%) |

**What it changes for the real unknown.** One block in five of area I's
constrained_unknown tier is mostly a copy of sequence elsewhere, and the
share is a property of the chromosome: chrY 72%, chr21 61% (15 of 20
blocks; the most constrained block of the chromosome, 51 kb at 5.50 Mb, is
100% duplicated with eight partners, some on unplaced scaffolds), chr22
42%, every other chromosome 11% or below. A duplicated block reads as
constrained across mammals because the alignment lands on the paralogue,
and is unscored by gnomAD because variant mapping fails in copies; the two
axes disagreeing (§6) had a cause, and the block organiser now has it as
`duplicated_fraction` per block. Copy first, attribute second.

**The heuristic.** The UNKNOWN classifier's `similar_to` (shared 20-mers
among the sixty largest blocks) fired once genome-wide, and that pair is a
curated duplication. It was never wrong; it was silent. The curated pairs
replace it as the duplication evidence and the heuristic stays as the
fallback for a chromosome without the track (D20's order).

**The elements.** 142 of 7,063 scored regulatory elements (2.0%) sit inside
a duplication: regulatory code pasted twice, of which the deletion model
scored one copy. They are listed per chromosome for the organiser.

**Paralogues** (the gene-level copy-and-paste) come from the Compara stream
of step 2 (§8), not from this track.

## 8. Step 4 built: motifs as the dependency list, and the operator search (2026-09-12)

`genomeos motifs --chrom C [--gene G]` and `--library L` (`genome/motifs.py`;
results `motifs_chr*`, `motifs_genome_wide`) derive what SEQUENCE-GRAMMAR.md
§2 designed and never built: a gene's `requires:` list from its own promoter
sequence. JASPAR 2026's CORE vertebrate collection (1,019 non-redundant
profiles, cached once) is scanned over every canonical promoter (TSS ± 1 kb)
and every scored element of a chromosome.

**Method, and why it is fast.** Each profile becomes a log-odds matrix; a hit
is a window at or above 85% of the way from the matrix minimum to its
maximum, on either strand. The scan is indexed rather than sliding: every
sequence is indexed by its 8-mers (or by its whole width for shorter
motifs), each motif enumerates by branch and bound the k-mers of its most
informative window that can still reach the threshold, and only those
positions are scored in full. Chromosome 21's 221 promoters against every
profile take 16 s; chromosome 19's 1,479 take 117 s.

**Control, and the lesson it taught.** Every promoter is also shuffled and
scanned the same way, and a factor's enrichment is the share of real
promoters with a hit over the share of shuffled ones. With a single-base
shuffle the enriched list was long zinc-finger matrices and nothing else,
because that shuffle destroys the CpG and GC runs real promoters have and
long GC-rich matrices then read as enriched against nothing. The control is
now a dinucleotide-preserving shuffle (Altschul and Erickson), which keeps
every dinucleotide count; the `requires:` lists then carry MEF2, FOX, PBX,
PKNOX1, MLXIP and ZNF143 on real promoters, and the most enriched class is
still, honestly, the long zinc-finger profiles (ZNF143, a promoter-binding
factor, is real; several others are the resolution of an 85% threshold on a
20-column matrix). A `requires:` entry is a factor with a hit whose
enrichment on the chromosome is 1.5 or more, eight per promoter on average.
Nothing here says a factor binds in a cell: that is the reader's claim, and
the decompiled locus (§9) prints both side by side.

**Operators, the genome (20,067 promoters, 24 chromosomes, 42 libraries).**
The genome summary pools every gene's list and, per BioLib library with
eight or more members, reports the factors that hit its members' promoters
at 1.5 times the genome share and the factor *pairs* present together in
members' promoters at twice what the two shares predict and in a quarter of
the members or more. Two readings, one per level.

At the factor level the libraries recover textbook biology without being
told it: REST, the neuronal-gene repressor, is enriched 7.7 times in
`systems.nervous` (28 members) and again in `blueprint.organ_nervous`,
`organ_eye`, `endocrine` and `circadian`; the interferon-response factors
IRF2, IRF3, IRF7 and IRF9 in `systems.immune` (2,821 members) and
`blueprint.organ_blood_immune`; the housekeeping promoter factors ZNF143 and
THAP11 in `core.translation` and `core.replication`, where they belong. The
heart library shows nuclear receptors (NR3C1, NR3C2) and PRDM9 but not the
GATA plus T-box pair the textbook names, which sits in enhancers more than in
promoters, so the promoter scan was the wrong place to look for it.

At the pair level the result is a lesson, not an operator table: the pairs
that recur beyond the two shares are near-identical matrices co-hitting
(MEF2A with MEF2D, MEF2B with MEF2D) and the CGG-repeat binder CGGBP1 with
GC-rich zinc-finger profiles (ZNF93, ZNF131), in every library alike. Both
are the promoter's GC content read twice, and a factor-pair statistic that
assumes independence between two shares cannot see that. An operator test
that means something needs the profiles clustered into families first (one
MEF2, one GC-rich zinc-finger class) and an expectation matched on GC, and
then enhancers, not promoters, for the developmental combinations. That is
the next step of step 4, and the table stands as `inferred` at low
confidence until it is done.

## 9. Step 6 built: the decompiled locus (2026-09-12)

`genomeos decompile GENE --chrom C` (`genomeos/decompile.py`) is the view
the conversation imagined, built by assembly alone: it reads the results on
disk and infers nothing, fetches nothing. For one gene it prints, as a
BioLang-flavoured block with an evidence note per line: the gene and its
canonical structure (GENCODE), the protein with its domains and pathways
(UniProt, InterPro, Reactome), origin and paralogues (Compara, when read),
the promoter's `requires:` (the motif scan), the node it sits in and where
the reader finds its promoter open (eleven ENCODE cell types), every element
whose deletion the model says moves it, with the predicted magnitude and
tissue, both constraint axes and the case, whether the element lies in a
segmental duplication and which motifs it carries, GTEx expression, and
last an `unknown { }` block naming every layer not yet read for that locus.
APP on chr21 decompiles with three elements (one `relaxed`, 46% of its bases
constrained across mammals and none of its kilobases among people) and, at
the time of writing, origin and paralogues as the only residue. The view is
the product form of area J: what the project knows about a locus, layer by
layer, with the gap stated.

## 10. Step 2 built: origin per gene, age per library, paralogues (2026-09-12)

`genomeos origin [--gene G | --library L]` (`knowledge/homology.py`,
`scripts/origin_genome_wide.py`, result `origin_genome_wide`) gives the
ladder the conversation asked for. Every Ensembl species is placed once on a
ladder of 25 clades shared with human, from Life through Eukaryota,
Metazoa, Bilateria, Chordata, the vertebrate and mammal grades down to
Homo, through Ensembl's taxonomy classification with proxy names for the
nodes it omits (Aves stands for Amniota, Amphibia for Tetrapoda,
Laurasiatheria for Boreoeutheria, Bacteria for Life). Then two Compara
homology dumps for human are streamed once over HTTP and never stored: the
vertebrate one (109 MB, 3.9 million rows, 9 s) and the pan-taxonomic one
(46 MB, 2.6 million rows, 5 s; plants, fungi, protists, invertebrates).
Per protein-coding gene the deepest clade with an orthologue is its
origin, the deeper of the two dumps winning; its paralogues come from the
same rows. Species the dumps name that the species list omitted (yeast, fly
and worm are outgroups of the vertebrate trees but not "vertebrate species")
are placed on the way and the stream is run again, which is what turned a
first result with nothing below Chordata into the table below.

**The genome's ladder (19,392 coding genes, 508 species placed).** The
first pass stopped at Chordata for every gene, because Ensembl's vertebrate
species list omits the outgroups its trees use; the second reached the
animals and eukaryotes through the pan-taxonomic dump; the third placed
that dump's 117 bacteria and archaea, whose strain-suffixed names the
taxonomy endpoint resolves only once the suffix is stripped, and gave the
ladder its bottom rung.

| ladder | genes | strata |
|---|---|---|
| Life core | 3,774 | an orthologue in bacteria or archaea |
| Eukaryote core | 8,250 | Eukaryota 7,403, Opisthokonta 847 |
| Animal core | 3,935 | Metazoa 3,267, Bilateria 668 |
| Chordate | 121 | |
| Vertebrate core | 1,815 | Vertebrata 412, Gnathostomata 540, Euteleostomi 863 |
| Tetrapod | 200 | |
| Amniote | 281 | |
| Mammal core | 905 | Mammalia 180, Theria 188, Eutheria 348, Boreoeutheria 176, Euarchontoglires 13 |
| Primate | 92 | |
| Ape | 18 | |
| Human-specific in this species set | 1 | |

One coding gene in five has a relative in bacteria (RPL3 reads Life core,
ACTB and HOXA1 eukaryote core, TP53 opisthokont, APP, RAG1 and HBB animal
core, INS vertebrate, KRTAP10-1 mammal); three in five are older than the
animals; one in ten is vertebrate, one in twenty mammalian; 111 genes are
primate or younger and one has no orthologue anywhere Ensembl looks.
Paralogues: 6,803 genes have none, 3,809 have ten or more; 64,050
`paralogue_of` edges joined the knowledge graph (394,068 edges now), one
per pair of compiled paralogues, with the origin as a node attribute. Mouse
is not a column: Ensembl's human dump lists Mus caroli, spretus, spicilegus
and pahari and the rat but not Mus musculus (its pairs sit in the mouse-side
dump, which `genome/mouse.py` reads for the node comparison), so the grade
is placed by its relatives and the decompiled locus claims no mouse count.

**The libraries, aged.** Every one of the 42 BioLib libraries with eight or
more placed members now carries its origin distribution. The order is the
one biology would give without being told: the oldest are `core.splicing`
(98.5% of members animal-wide or older), `core.replication` (96.0%),
`core.translation` (95.7%) and `timer.telomere` (95.5%), then the
segmentation and circadian timers and the degradation, repair and
chromatin cores (93 to 94%); the youngest are `systems.ligand_receptor_protocol`
(54.4% animal or older, 38.8% vertebrate to tetrapod: the signalling
protocol is largely a vertebrate invention), `blueprint.organ_skin` (15.2%
amniote or younger, the keratin-associated proteins), `systems.immune`
and `blueprint.organ_blood_immune` (6% and 5% amniote or younger). A
library whose members span the whole ladder is a candidate for splitting;
the ligand-receptor protocol is that library.

**Two caveats the data forces.** An origin from the vertebrate dump is a
lower bound on age (a gene lost in every sampled outgroup reads younger
than it is); an origin from the pan-taxonomic dump can overshoot, because
at the deepest levels Compara's "orthologue" is a family-level call: HOXA1
reads eukaryote-wide through plant homeobox proteins, and OCA2 through
transporter relatives. The ladder label is therefore right at the grade
(eukaryote core, animal core, vertebrate core), and a single gene's exact
stratum below Metazoa should be read with the species count beside it.
Ten species the taxonomy endpoint cannot resolve under any name are placed
by an explicit table (mink and sperm whale renamed in NCBI, two birds, six
pan-taxonomic invertebrates and a liverwort), and eight bacteria and
protists with irregular strain names stay unplaced; each is redundant with
others at its stratum, so no origin depends on them. A dump that yields
orthologues in fewer than fifty species now fails loudly rather than read
as a genome of human-specific genes.

The decompiled locus (§9) now prints origin and paralogues for every
coding gene: APP is animal core with one paralogue, APLP1, and has no
unknown residue left among the layers GenomeOS reads.


## 11. Step 5 built: phylogenetic profiling between libraries (2026-09-12)

`genomeos profiling [--library L]` (`knowledge/profiling.py`, result
`profiling_genome_wide`) reads the presence matrix the origin job writes
beside its summary (`origin_presence_genome_wide`: one bit per gene per
species, 355 species with an orthologue of some human gene, ordered by
clade) and gives every library with eight or more placed members its
profile, the share of members with an orthologue in each species; its gain
curve, the share of members whose origin is at or before each stratum; and
its residual profile, what it does beyond the genome's own presence
pattern. Between libraries, the Pearson correlation of residual profiles is
the "gained and lost together" candidate the conversation asked for.

**What it finds.** The developmental blueprint is one history: the axis,
segmentation, germ-layer and organ libraries correlate at 0.97 to 0.99 with
one another (segmentation and its clock at 1.0, which is the same gene list
under two names, and a reminder that shared members inflate the pair). The
least alike are the cores against the vertebrate-born signalling protocol:
`core.replication` against `systems.ligand_receptor_protocol` at −0.87,
`core.splicing` and `core.chromatin_epigenome` against the same at −0.84,
`timer.telomere` against `extracellular_matrix_adhesion` at −0.82. The
immune system's nearest histories are the senescence checkpoint, puberty
and the telomere machinery (0.50 to 0.57), which is the vertebrate-and-later
signature shared with the timers that tune adult life.

**How far to trust it.** A high correlation says two libraries appeared and
persisted in the same species, not that one calls the other; libraries
that share members share history by construction; and 102 of the 355
species are bacteria and archaea, so the profiles weigh the Life boundary
heavily, which is exactly the boundary the cores sit on. The Reactome
hierarchy the libraries were built from is the check, and the residual
correlation is `inferred` evidence. The step the conversation called
covariance is built; making it a dependency claim needs a gene-level
profile and a null over shuffled libraries, which is recorded as its next
step.


## 12. Step 7 built: one locus across species (2026-09-12)

`genomeos across --locus ZRS | HERC2_OCA2` (`knowledge/across.py`, results
`across_ZRS`, `across_HERC2_OCA2`) is the experiment the conversation ended
on: take a developmental locus with every layer available and see whether
its grammar is still there in other vertebrates. For each species Ensembl
Compara's pairwise LASTZ alignment of the human element is fetched (under a
second, cached), giving the other locus, the share of the element that
aligns and the identity over aligned columns; both constraint axes are read
over the human element; and the JASPAR scan runs over both sides of each
alignment.

**The method had to be corrected before it said anything.** The first
version counted a factor as conserved when its motif hit anywhere in the
other species' aligned sequence, and called 170 factors conserved across the
ZRS to fish: at 88% identity over 800 bases nearly every motif is "present".
A factor now counts only when its hit sits at the *same aligned site*: the
human hit's start is mapped through the alignment columns, gaps kept, and the
other species' hit must land within 15 bases of it. Two further defects were
caught in review of the code, not by the tests: the alignment summary had
stripped the gaps the mapping needs, and the stored factor list was capped at
60, so a factor ranked 61st read as "no hit". Both are fixed and tested (an
insertion case, and every factor kept).

**The ZRS** (VISTA hs2496, 789 bp, limb enhancer of SHH inside LMBR1, a
megabase from its target):

| species | element aligned | identity |
|---|---|---|
| mouse | 100% | 88.2% |
| opossum | 100% | 88.6% |
| chicken | 100% | 87.5% |
| zebrafish | 9.8% | 80.5% |
| frog | no pairwise alignment in Compara | |

54.6% of its bases are constrained across mammals (mean phyloP 3.24); among
people its kilobase reads Z 1.97, just under the 2.18 bar. Of 606 factors
whose motif hits the human element, 151 hold their site in none of the four
species and 23 hold it in all four. Those 23 are one class read many times:
homeodomain matrices sharing the TAAT core (CRX, EMX2, GSC2, NKX6-1, NKX6-2,
PDX1, the posterior HOX matrices HOXA9, HOXA10, HOXB5, HOXC8, HOXC10, and
CDX1, CDX4), the HOX cofactors MEIS2, MEIS3 and PBX2, and one ETS composite
(ETV5::FOXO1). Among the factors the literature names for the ZRS, HAND2
holds its site in the three amniotes and not in fish, HOXD13 in mouse and
chicken, ETS1 and ETV5 in chicken and opossum. The site-conserved core is
therefore the HOX, PBX and MEIS class with an ETS site, which is the known
logic of this enhancer, and the tetrapod-only HAND2 site fits a limb program
that fish fins do not share in full.

**The HERC2 enhancer of OCA2** (957 bp around rs12913832) aligns to mouse
(74.4% of the element at 70.6%) and opossum (28.9% at 79.4%) and to no
chicken, frog or fish genome: it is a mammalian element. Only 7.7% of its
bases are constrained across mammals, while its scored kilobase is
constrained among people (Z 3.46), the `recent` case. Of 551 factors hitting
it, 75 hold their site in both mammals; with two close species that is a low
bar. What holds is the E-box class of the MiT family (TFEB, TFEC) with USF1
and MLXIPL, a SOX10 site and RUNX1. MITF and SOX10 drive the melanocyte
program that includes OCA2, and the MiT E-box and SOX10 site holding in
mammals is consistent with it; the MITF matrix itself does not hold its site,
and LEF1 holds in opossum only.

**What the experiment says, and does not.** The abstraction survives its
first test: on the element conserved to fish, the grammar that holds its site
is the one biologists describe for it, found without being told. It does not
yet give a clean operator list, because matrix families (one TAAT core, one
E-box) are counted per matrix, which is the redundancy step 4 found at
promoters; clustering JASPAR into families before counting is the shared fix.
An alignment says a site is kept, not that it is used: VISTA's mouse assay
says the ZRS is active in the limb, and for the rest that is the reader's
claim. Frog is missing because Compara holds no human-to-Xenopus pairwise
alignment for this pair, recorded as unavailable rather than as absent.
