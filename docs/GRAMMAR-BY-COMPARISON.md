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
| HPRC 90-way and v2.1 SVs | 90 human assemblies aligned (hg38, CHM13, 88 haplotypes of 44 people); arrangements against hg38; SVs from 233 assemblies | UCSC `hprc90way` MAF at `gbdb/hg38/hprc/cactus90way/<chrom>.maf` by the block offsets the REST API returns; `hprcArrV1` and `hprc2v21Sv` as bigBed ranges (`attribution/human_panel.py`) | read 2026-09-14 on chr21: 3.28 GB in 236 s, 1.05 M variable columns per assembly; per-base presence, alleles and inserted lengths give block classes against GC- and timing-matched windows, fixed/storage/cannot-place catalogues with value domains, and a kilobase comparison with Gnocchi (weak agreement, Spearman -0.08); genome 266.6 GB (`human_panel_chr21`, ATTRIBUTION.md "Many human genomes") |
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


## 13. Refinement: families, composition and sites (2026-09-12)

Sections 8 and 12 both ended on the same defect: matrices were counted one
by one, so one binding mode read by many JASPAR profiles looked like many
factors. Three corrections, each tested, now sit under both the promoter
operators and the cross-species grammar.

**Families.** JASPAR's TRANSFAC release carries the curated TFClass family
and class of every matrix (`load_families`). A factor is counted as its
family, so MEF2A and MEF2D are one unit and cannot form a pair; a
heterodimer counts as its two families joined; and C2H2 zinc-finger
factors stay individual, because those families are structural ("more
than three adjacent zinc fingers" holds 195 matrices binding unrelated
sequences). 158 matrices in the promoter lists become 110 units.

**Composition.** A pair statistic that assumes independence between two
shares finds the promoter's composition twice. GC was the known case
(CGGBP1 with ZNF93, in 27 libraries). The check for what remained found a
second one: promoters carrying both ZNF135 and ZNF460 sites are 19.4% Alu
on chr1, chr17 and chr19, against 8.8% with one and 1.1% with neither, so
the last surviving "operator" was a transposon matched by two KRAB zinc
fingers. `promoter_composition_genome_wide` now records GC and the
interspersed-repeat share (RepeatMasker SINE, LINE, LTR, DNA) of all 20,067
promoters (median GC 0.56, median repeat share 12%, a quarter repeat-free),
and every expectation is taken within GC-decile by repeat-share strata (39
strata of 30 promoters or more; thinner ones fall back to their GC decile).

**The promoter result, corrected.** No factor-family pair recurs across a
library's promoters beyond its composition-matched expectation, in any of
the 42 libraries; 84 pairs that passed the naive test are explained by GC
or repeats. The single-family enrichments survive the same correction:
REST in `systems.nervous` at 5.7 times its matched expectation (28
members), the THAP and ZNF143 promoter families in the translation and
replication cores, AP-1 in the circadian timer. Promoters say which
families a library's genes share; they do not carry pair logic. That
matches the textbook placement of combinatorial logic in enhancers, and it
is now a measured negative rather than an impression.

**Sites.** Families do not finish the job across species either, because
the homeodomain families (HOX, NK, paired-related, HD-LIM) are separate in
TFClass but share the TAAT core. The honest unit is the site: every held
factor's hit is mapped back to the genome through its alignment block's
start (positions had been recorded in whichever species' concatenated
alignment scored best, frames that differ when coverage is partial), and
hits overlapping on the human genome are one site.

| locus | factors hitting human | held in every species | as families | as sites |
|---|---|---|---|---|
| ZRS | 606 | 23 matrices | 10 | **10** |
| OCA2 enhancer | 551 | 75 matrices | 25 | **11** |

**Correction, the same day: the site table above is superseded.** The
site call behind it had a flaw found while designing the panel below. A
factor counted as held "in every species" when each species held *some*
hit for it, but the site's genome coordinate came from the first species
holding it, so a factor held at one place in mouse and another in zebrafish
(which aligns to 77 of the ZRS's 789 bases) was reported as one site. The
call now places every species' held hit on the genome through that
species' own alignment blocks, clusters the hits into sites, and requires
every core species (those aligning over half the element) to hold the same
site, with one factor holding it in all of them. Under that definition a
single locus says little: the ZRS holds 15 such sites covering 79% of its
bases and 151 families among them, because at 88% identity across three
amniotes nearly every window of conserved sequence keeps some matrix; the
OCA2 enhancer has one core species (mouse) and cannot be judged. Four limb
enhancers and four constraint-matched VISTA negatives held 10.6 to 15.9 and
7.6 to 19.8 sites per kilobase: the count does not separate enhancers from
inactive conserved sequence. The grammar question is therefore a
comparison, positives against negatives family by family, which is the
panel in §15. The paragraphs below describe the first site call and are
kept as the record of what it showed.

For the ZRS the family count already was the site count: the element keeps
two HOX sites, two TALE (MEIS and PBX) sites, two NR2 sites, the ETS
composite and three single zinc-finger sites across 789 bases from human to
zebrafish, a multi-site grammar and the known HOX, PBX, MEIS and ETS logic.
For the OCA2 enhancer the families still overcounted: one 11-base site at
element position 417 is read by 38 matrices from five homeodomain families,
and a 36-base stretch at 246 by 13 matrices from 11 families. Its grammar is
11 sites in two mammals, the E-box (bHLH-ZIP, RUNX and PAS matrices at one
site) and a SOX site among them.


## 14. Refinement: library histories against an age-matched null (2026-09-12)

The library-history table of §11 correlated presence profiles and found
the developmental libraries sharing one history (0.97 to 0.99) and the
cores opposed to the ligand-receptor protocol (−0.87). Three things make
two libraries correlate without a shared history: identical members, the
same age, and the 102 bacterial species dominating the profile. `pair_null`
(`knowledge/profiling.py`) takes the correlation on the members exclusive
to each library, over all species and over the 245 eukaryotic species
alone, and compares it with 200 pairs of random gene sets matching each
side's origin-stratum make-up. Each species is held as a bitmask over genes,
so a profile is 355 popcounts and the 860 testable pairs take 195 s.

**Both headline findings of §11 are age.** Random age-matched sets already
correlate at 0.55 to 0.73 where the developmental libraries sit, so axes
with segmentation reads z 1.9, heart with musculoskeletal 2.2, germ layers
with the nervous system 1.6, none significant. The cores against the
ligand-receptor protocol match their null (replication −0.87 against −0.79,
z −0.9). Segmentation and its clock are the same genes under two names and
are flagged instead of tested.

**What is more alike than age predicts.** 35 of 860 pairs at a false
discovery rate of 5%, led by the immune and signalling libraries: immune
against the ligand-receptor protocol (observed −0.76, age-matched −0.93,
z 9.3), degradation and cytoskeleton against the same protocol (z 8.5),
the cell cycle against immunity (z 6.7), blood with the growth-longevity
axis (observed +0.72 against −0.27). Immune and signalling gene families
expand in lineages, so the obvious confound is gene-family structure; with
every gene of more than four paralogues removed, and again of more than
one, seven of the top eight pairs stay significant (z 3.7 to 12.7) and one
disappears (the protocol with the telomere machinery, z 5.0 to 0.1). The
test is one-sided, "more alike than age"; pairs far less alike than their
age (z down to −35) are recorded and not called.


## 15. The panel: VISTA enhancers against matched negatives (2026-09-12)

The question §13 left: do enhancers keep particular factor families in
place across species more often than conserved sequence that is not an
enhancer? `scripts/across_panel.py` (`knowledge/across.py`, result
`across_panel_vista`) drew 40 limb-only and 40 neural-only VISTA positives
and 80 negatives sampled to the positives' constraint bins, all 500 bp to
2 kb with 30% to 90% of bases constrained across mammals; aligned each to
mouse, opossum, chicken, frog and zebrafish; and recorded the families held
at factor-strict sites (one genome position, the same factor, every species
aligning over half the element). Loci with fewer than two such species were
set aside, leaving 39 limb, 33 neural and 76 negatives, balanced in
constraint (medians 0.53, 0.62, 0.57), length and the share where chicken
or zebrafish align (limb 79% and 15%, neural 91% and 30%, negatives 83% and
18%). 1,258 s for 160 loci.

| comparison | held sites per kb, median | density p (one-sided) | families tested | families at q ≤ 0.05 |
|---|---|---|---|---|
| limb 39 vs negatives 76 | 10.2 vs 10.6 | 0.88 | 274 | 0 |
| neural 33 vs negatives 76 | 10.7 vs 10.6 | 0.37 | 276 | 0 |

The nearest calls do not survive the correction for testing 275 families:
ATF4-related factors held in 13 of 33 neural enhancers against 8 of 76
negatives (p 0.0008, q 0.23), ZBTB11 in 8 of 39 limb enhancers against 2 of
76 (p 0.003, q 0.71). The HOX family, the reading the ZRS invited, is not
among the leading limb families.

**What the negative says.** It does not say enhancers lack a grammar. It
says this readout cannot see one: JASPAR best hits at 85% of the matrix
range, held at aligned positions across amniotes, occur at the same density
and with the same families in active enhancers and in conserved elements
that drove no expression. Four reasons are recorded rather than tried one
after another until something passes. The threshold is permissive (a human
element of 1 kb draws some 600 matrix hits), so conservation, not binding,
decides what is held. Only each factor's best hit is scanned. VISTA
negatives are conserved elements tested at one embryonic day, many of them
enhancers of another tissue or stage, so the control is regulatory-like
sequence rather than inert sequence. And 275 families at 35 to 76 loci per
group leave power only for large differences. The ZRS reading of §12
remains a single-locus observation consistent with the literature, not a
validated signal; the promoter operators of §13 were a measured absence,
and across species the grammar is, for now, not measured.

A useful next step has to change the readout, not the threshold: a
sequence model's in-silico mutagenesis of each element (which bases change
predicted activity, AlphaGenome being the project's model), held against
the same positives and negatives, asks whether the bases that matter are
the ones that are kept.


## 16. Measured: which bases matter, and what marks them (2026-09-13)

Section 15 ended with a negative it could not explain, because nothing in it
measured which bases of an enhancer matter. Kircher et al. 2019 measured
exactly that: saturation mutagenesis read out by a reporter assay over 21
regulatory elements, nearly every single-base substitution with its effect on
activity and a p-value (44,658 measurements on GRCh38; GEO GSE126550, the
table published with the lab's data portal). `genomeos/knowledge/satmut.py`
and `scripts/satmut_grammar.py` read it, call the bases whose substitution
changes activity, and ask what marks them. No AlphaGenome quota, so it runs
beside the deletion-scoring chain.

A base counts as **functional** when some substitution at it is significant
(p < 1e-5, at least 10 barcodes, the portal's defaults) and **strong** when
that substitution also moves activity by at least 0.25 log2, about 19%: at
deep coverage a p-value alone flags very small effects, and the share of
functional bases swings from 1% of FOXE1 to 77% of IRF4. Both definitions are
reported throughout. Loci measured twice or more (SORT1, PKLR, LDLR, TERT,
the ZRS) contribute one primary experiment each; the repeats agree at r 0.63
to 0.97 on shared substitutions, so the measurement is sound.

**The panel's negative, explained.** The share of an element's bases covered
by some JASPAR site, and what that buys:

| match threshold | bases in a site | functional bases recovered | enrichment |
|---|---|---|---|
| 0.80 | 100% | 100% | 1.00 |
| 0.85 | 99.9% | 100% | 1.00 |
| 0.90 | 96.9% | 97.6% | 1.01 |
| 0.95 | 65.0% | 72.0% | 1.11 |

At 0.85, the threshold §8 and §15 used, some profile of the 1,019 covers
essentially every base of a regulatory element. "The site is held across
species" was therefore a statement about conserved sequence, not about
binding, which is why enhancers and conserved negatives looked alike. Only at
0.95 is a site selective, and then it enriches functional bases 1.11 times
(1.16 for strong ones).

**What does predict a functional base.** Area J's two candidate marks, ranked
by the area under the ROC curve per locus (0.5 is chance):

| mark | median AUC, functional | median AUC, strong | above chance |
|---|---|---|---|
| best motif score covering the base | 0.58 | 0.59 | 19 and 18 of 21 loci |
| phyloP across 241 mammals | 0.54 | 0.54 | 15 and 13 of 21 |

Both are weak, and the motif score is the better of the two, which is the
reverse of what the project has leaned on. Per locus the spread is wide:
phyloP reaches 0.74 at the LDLR promoter and sits at 0.45 for SORT1 and 0.40
for FOXE1; the motif score reaches 0.71 at the TERT promoter.

**Where they combine, they do better than either alone.** Among bases
constrained across mammals (phyloP ≥ 2.27), those inside a strict (0.95)
motif site are functional 36.5% of the time against 21.1% outside one
(p 1e-15); for strong effects 24.1% against 13.3% (p 5e-11). That is the
first quantified support in area J for the idea the conversation started
from: conservation plus a recognisable site is a better statement about
function than conservation alone. It is an enrichment of 1.7, not a rule.

**No family stands out.** Of 89 families with at least five sites, none has a
higher share of functional sites after correction; the leaders (KLF1, KLF5,
KLF2, KLF14, SP) are all the same GC box, and the strong definition adds
Ets-related and TEF-1 at q above 0.05. At 21 loci this asks too much of the
data.

**The ZRS.** The element §12 read as textbook HOX, PBX and MEIS logic has 6
to 8 strongly functional bases of 485 measured, and 87 to 102 by the
significance-only definition. Measurement gives that reading no support; it
remains a motif list consistent with the literature.

**What this changes.** The threshold calibration is now a comment where the
scanner is defined: 0.85 stays for per-gene `requires` lists, which take each
factor's best hit, and site-coverage questions must use 0.95. The area's
honest summary is unchanged in shape but sharper: comparison across species
and people places and dates function well, motif sites at a usual threshold
say nothing about it, and the two marks together are modestly informative
about individual bases. A readout that would do better is a model trained to
predict activity, scored the same way against these 9,834 measured bases.

## 17. Grammar as a test that can fail: arrangement against counts (2026-09-16)

Every section before this one used the word "grammar" for something it did not
test: §8 read a promoter as a list of factors, §13 looked for pairs, §15 asked
whether sites are held across species, §16 asked which bases matter. None of
them asked the question the DNA-as-code reading actually rests on — does *word
order* carry information? If an enhancer is a sentence, then how far apart two
factor sites sit, on which strand, and in which turn of the helix should predict
activity beyond which factors are present and how often. That is a claim that
can fail, and this section fails it.

**The instrument.** ENCODE4's joint lentiMPRA library (`attribution/mpra.py`):
51,376 elements of 200 bp, each with measured log2(RNA/DNA) in K562, HepG2 and
WTC11, so the same sequence is read three times in three cells.
`genomeos/attribution/motif_grammar.py` and `scripts/motif_grammar.py` build
three nested feature sets per element and compare them out of sample
(`motif_grammar`, and `motif_grammar_preregistration` for the training-only
stage):

- **(a) composition**: GC, GC squared, CpG observed over expected.
- **(b) plus strict counts**: log1p sites per TFClass family unit — JASPAR 2026
  CORE vertebrates at relative score **0.95**, the threshold §16's measurement
  calibrated, collapsed to one site per family per position (a GC box read by
  KLF1, KLF5 and SP2 is one site, not three) — the 40 most frequent units on the
  training chromosomes, plus total sites and distinct families. Median 34 sites
  per 200 bp element.
- **(c) plus grammar**: for each of the 55 unordered pairs of the ten most
  frequent family units (a family with itself included), six features — the
  number of non-overlapping site pairs with an edge-to-edge gap under 10, 10–20,
  20–50 and 50+ bp, the number on the same strand, and the number whose centres
  sit a whole helical turn apart (10.5 bp, within a quarter turn). 330 grammar
  columns, 376 in all.

Ridge, standard library only; the penalty is chosen by fold-by-chromosome
cross-validation *inside* the 45,228 training elements (chr1–7, 10–20, X, Y).
The falsifier is a **shuffled-grammar control**: the same grammar features after
the site labels and strands are permuted within each element, which keeps every
count and every position and destroys only which family sits where and how it is
turned.

**Pre-registered, in code, and committed before the held-out chromosomes were
read** (`PREREGISTERED` in the module): *(c) beats (b) in Spearman with measured
activity on held-out chromosomes in each cell line — the 95% bootstrap interval
of ρ(c) − ρ(b) over held-out elements lies above 0 in K562, HepG2 and WTC11. The
arrangement reading additionally requires ρ(c) − ρ(c shuffled) above 0 by the
same interval.* Held out: chr8, chr9, chr21, chr22 — 6,146 elements, scored once.

| cell line | (a) GC, CpG | (b) + strict counts | (c) + grammar | (c) shuffled | (c) − (b), 95% CI | (c) − shuffled, 95% CI |
|---|---|---|---|---|---|---|
| K562 | 0.200 | 0.434 | 0.436 | 0.437 | +0.0021 [−0.0049, +0.0085] | −0.0013 [−0.0081, +0.0051] |
| HepG2 | 0.167 | 0.402 | 0.397 | 0.400 | −0.0056 [−0.0120, +0.0014] | −0.0036 [−0.0106, +0.0038] |
| WTC11 | 0.348 | 0.424 | 0.426 | 0.422 | +0.0022 [−0.0050, +0.0097] | +0.0047 [−0.0023, +0.0124] |

Spearman with measured activity on the held-out chromosomes; 1,000 bootstrap
resamples of elements. The same comparison inside the training folds, where
nothing was held out, gives the same picture: (b) 0.420, 0.408, 0.424 against
(c) 0.421, 0.408, 0.425.

**Verdict: the pre-registered claim fails in all three cell lines.** Motif
arrangement — pairwise spacing, relative orientation, helical phase, as read
here — adds nothing detectable over strict motif counts. The point estimates
straddle zero (+0.002, −0.006, +0.002 on a Spearman correlation of 0.4), the
intervals contain zero in every line, and the shuffled control does as well as
the real arrangement, twice out of three times better. Where §16 found a
measured absence at the level of bases, this is a measured absence at the level
of arrangement.

**What did work, and is worth keeping.** Strict counts are the one large effect
in the table: they lift Spearman by +0.234 in K562 and +0.235 in HepG2 over
composition alone (intervals [+0.208, +0.262] and [+0.205, +0.263]), which is
the first time in area J that JASPAR sites have predicted a measured quantity at
scale. *Which* factors have sites matters; *how they are arranged* does not, to
this readout. WTC11 is the exception that says something about the assay: GC and
CpG alone reach 0.348 there against 0.167 in HepG2, so in the iPSC line much of
what the reporter measures is promoter-like composition, and counts add only
+0.077.

**What this negative does and does not license.** It does not say enhancer
grammar is a myth: the literature's clearest cases (the ZRS, the interferon-beta
enhanceosome) are single loci with resolved structure, and the arrangement they
depend on may be too rare to lift a genome-wide correlation. Four limits are
recorded rather than tuned away. A 200 bp element holds few well-separated pairs,
so the long-range arrangement of a full enhancer is out of the assay's reach. A
site call is still a matrix score, not a footprint: at 0.95 the median element
carries 34 of them, more than a 200 bp enhancer plausibly binds, so the pair
features count mostly unbound sites. Ten families give 55 pairs, and the pair
that matters for one cell type may sit outside them. And a reporter measures
episomal activity out of chromatin, where phasing against the nucleosome — the
one mechanism that would make helical phase matter most — is absent by
construction.

The honest reading of area J after this section: comparison across species and
people places and dates function; the presence of strict, family-collapsed motif
sites predicts measured activity at genome scale; their arrangement, at this
resolution and in this assay, does not. The next readout that could overturn it
is not another feature set but a model that reads the sequence itself, scored
against these same held-out elements.

## 18. The count positive, tested for transfer: an assay and a question (2026-09-17)

Section 17 failed the arrangement claim and left one thing standing: strict,
family-collapsed JASPAR site counts lift Spearman with lentiMPRA activity by
+0.234 (K562) and +0.235 (HepG2) over GC and CpG composition. That was measured
on 200 bp fragments carried on a reporter integrated out of their own
chromosome, and §17 said so. This section asks whether the positive is about
enhancers or about the construct, by moving it twice: to a **different assay**
(VISTA, transgenic mouse embryos) and to a **different question** (the 882
constrained-unknown blocks area I calls the real unknown).
`genomeos/attribution/motif_transfer.py`, `scripts/motif_transfer.py`; results
`motif_transfer_preregistration` (written and committed first) and
`motif_transfer`.

### Reading 1: VISTA, in vivo

2,223 human sequences with a curated in-vivo verdict (1,133 positive in named
tissues, 1,090 negative), each with Zoonomia phyloP already measured over it
(`attribution/vista.py`, `attribution/constraint.py`). Six nested ridge models,
counts read **per kb** because a VISTA element runs from 44 bp to 8 kb:

length | + composition | + conservation | + counts | conservation + counts | + shuffled counts

**Pre-registered in `PREREGISTERED_VISTA` and committed before the held-out
chromosomes were read**: counts separate the classes, counts add over
composition, **conservation + counts beats conservation alone** — the transfer
claim, and the central one because VISTA elements were *chosen* for conservation
— and counts beat their own dinucleotide-shuffled control. Held out before
anything was scored: chr4, chr5, chr8, chr9, chr13, chr18, chr21, chr22, chrX —
643 elements, 351 positive, a superset of §17's four chromosomes because VISTA
is 23 times smaller than the lentiMPRA library. Coverage: 2,223 of 2,223
elements scanned, 0 without sequence or sites, 2 without conservation (the
training mean stands in for those two, recorded), 0 with a status other than
positive or negative.

**What is unmatched.** VISTA negatives are sequences that *failed to drive
expression*, not a background drawn to match anything, so the positive-negative
gap below is not an effect size. On the held-out 643, positives sit above
negatives at AUROC 0.566 on length, 0.556 on GC, 0.509 on CpG o/e, 0.550 on
phyloP mean, 0.553 on constrained fraction — and **0.517 on total sites per kb**
(median 182.7 against 182.7). Whatever counts contribute, it is not that
positives carry more sites.

| model | AUROC | 95% CI | AUPRC |
|---|---|---|---|
| length | 0.566 | [0.522, 0.608] | 0.570 |
| + composition | 0.561 | [0.516, 0.606] | 0.589 |
| + conservation | 0.596 | [0.554, 0.638] | 0.626 |
| + counts | 0.630 | [0.584, 0.671] | 0.662 |
| conservation + counts | **0.655** | [0.612, 0.699] | 0.689 |
| + shuffled counts (null) | 0.562 | [0.516, 0.606] | 0.590 |

| difference | point | 95% CI |
|---|---|---|
| counts − composition | +0.069 | [+0.033, +0.104] |
| counts − shuffled counts | +0.067 | [+0.029, +0.111] |
| **(conservation + counts) − conservation** | **+0.060** | **[+0.022, +0.099]** |
| (conservation + counts) − counts | +0.026 | [+0.008, +0.044] |

Label-permutation null for the counts model: mean AUROC 0.4985, 95% interval
[0.454, 0.545], and 0 of 1,000 permutations reached the observed value.

**Verdict: the pre-registered claim holds in all four parts.** Strict
family-collapsed motif counts separate VISTA positives from negatives in
transgenic mouse embryos; they add over composition, they add over conservation,
and they beat a dinucleotide shuffle of the same element that preserves GC, CpG
and every dinucleotide count. §17's positive is not a property of the episomal
reporter. Conservation and counts also add to *each other* in both directions,
which is §16's base-level finding (conservation plus a strict site beats either
alone, enrichment 1.7) reappearing at the level of whole elements.

The size is the thing to keep honest about: AUROC 0.655 on an unmatched
comparison is a weak discriminator, and the whole effect of counts over
conservation is six points of AUROC.

### Reading 1b: the tissue group, where counts add nothing

Among positives only, one-vs-rest on the same nested sets. The group is
predictable — and then the control says by what:

| tissue group | held-out positives | length | + composition | + counts | counts − composition |
|---|---|---|---|---|---|
| neural | 258 | 0.619 | 0.744 | 0.763 | +0.019 [−0.017, +0.056] |
| heart | 43 | 0.716 | 0.839 | 0.846 | +0.006 [−0.032, +0.046] |
| limb and mesenchyme | 98 | 0.575 | 0.545 | 0.584 | +0.039 [−0.010, +0.089] |
| liver | 1 | — | — | — | too few positives |
| pancreas | 5 | — | — | — | too few positives |

60 of the 1,133 positives drive only tissues outside the five groups and carry
no group label. Label-permutation nulls sit at 0.498, 0.498 and 0.499.

**Which tissue a VISTA positive drives is predictable from its sequence — heart
at 0.846, neural at 0.763 — and motif counts add nothing detectable to it.**
Length and GC-CpG composition get there alone. This is the reading the
count-based grammar story would most want, and it fails: the factors whose sites
are present do not say which tissue the element works in, at this resolution.
(The nested comparison for tissue was added as a control after the first
scoring; the pre-registered tissue statement was only that *some* group would be
predicted above a permutation null, and it is. The control is what makes the
number readable, and it is reported whichever way it came out.)

### Reading 2: the real unknown, decided on measurement

The 882 constrained-unknown blocks that are not copies were measured on
2026-09-16 by an AlphaGenome deletion sweep and came out **below** the neutral
tier per element, 0.248 against 0.293, z −5.7 (`docs/ATTRIBUTION.md`,
`constrained_unknown_targets`). The obvious test — do motif counts rank the 331
blocks that carry a lead above the 551 that do not — **was not run as the
primary reading, because that lead is itself an AlphaGenome output**: a motif
model agreeing with it would be two model readings agreeing with each other, not
evidence. So the primary reading is measured, and coverage comes first:

| | real unknown | neutral tier |
|---|---|---|
| blocks | 882 | 2,632 |
| blocks holding at least one lentiMPRA element | **161** | 298 |
| blocks holding none | **721** | 2,334 |
| measured elements | **227** | 387 |
| blocks with an element the sweep tested | 531 | 1,181 |
| blocks with a sweep lead | 331 | 768 |
| median block length | 11.7 kb | 7.9 kb |

The 531 and 331 reproduce the sweep's own arithmetic exactly, which is the check
that both readings are looking at the same blocks. Every one of the 614
lentiMPRA elements inside a real-unknown or neutral block was held out of the
count model's fit **by location**, leaving 50,760 training elements.

**(1) The measured level.** lentiMPRA log2(RNA/DNA) of the 227 elements inside
real-unknown blocks against the 387 inside neutral-tier blocks:

| cell line | AUROC | z | median, real unknown | median, neutral |
|---|---|---|---|---|
| K562 | 0.462 | −1.59 | −0.504 | −0.421 |
| HepG2 | 0.489 | −0.45 | −0.424 | −0.411 |
| WTC11 | 0.472 | −1.18 | −0.857 | −0.758 |

**The measured reporter puts the real unknown at the same level as the neutral
tier, leaning below it in all three cell lines and reaching significance in
none.** It does not contradict the sweep; it points the same way with a tenth of
the evidence. There is no contradiction to report, and 227 measured elements
over 161 of 882 blocks is why: the measured instrument that exists for this
sequence covers 18% of the blocks.

**(2) The measured transfer.** Spearman of the count model (b) with measured
activity on those same 227 held-out-by-location elements, against composition
(a):

| cell line | (a) composition | (b) + counts | (b) − (a) | 95% CI |
|---|---|---|---|---|
| K562 | 0.131 | 0.369 | **+0.238** | [+0.097, +0.378] |
| HepG2 | 0.021 | 0.410 | **+0.389** | [+0.234, +0.528] |
| WTC11 | 0.290 | 0.336 | +0.046 | [−0.056, +0.142] |

+0.238 in K562 against §17's +0.234 over the whole 6,146-element held-out set:
the count positive reproduces inside the real unknown, on sequence chosen for
being constrained and unexplained, to two decimal places. WTC11 behaves as it
did in §17 (composition alone already carries the iPSC line). The same
comparison on the 387 neutral-tier elements gives +0.081 [−0.014, +0.165],
+0.015 and +0.033 — no cell line significant, and the intervals overlap the real
unknown's, so "counts work better here than in neutral sequence" is a
possibility this cannot resolve, not a finding.

**(3) The predicted level** — a prediction, not a measurement. The same model
over every non-overlapping 200 bp window of all 882 blocks (152,556 windows, 68
dropped for N bases, 0 blocks left without a window) against 882 length-matched
neutral blocks (151,810 windows), mean per window per block:

| cell line | real unknown vs matched neutral | vs its own dinucleotide shuffle |
|---|---|---|
| K562 | 0.497 (z −0.25) | 0.514 (z +1.03) |
| HepG2 | 0.499 (z −0.10) | 0.491 (z −0.64) |
| WTC11 | 0.498 (z −0.12) | 0.519 (z +1.35) |

Two flat readings, and the second is the more interesting: **applied blind
across 30.6 Mb of genome, the count model cannot tell real sequence from a
dinucleotide shuffle of the same sequence.** The model that reproduces §17's
+0.238 on 227 assayed elements carries no signal at all when swept over whole
blocks. A lentiMPRA element is a *selected* 200 bp — a candidate cCRE — and the
model inherits that selection; a uniform tiling of a 12 kb block is mostly
sequence the library would never have contained.

**(4) Agreement with the sweep**, reported last and **not validation**: both
sides are model outputs. On the 531 blocks with a tested element the predicted
mean scores AUROC 0.542 [0.494, 0.593] against carrying a lead, permutation null
0.500, and the interval contains 0.5. Over all 882 the mean gives 0.484 while
*block length alone* gives 0.718 and *the number of tested elements alone* gives
0.837 — the label is mostly a statement about how much of a block was tested,
which is exactly why the maximum-over-windows score (0.679) is not evidence
either. Two model readings of these blocks do not agree, and nothing follows
from that in either direction.

One observation, not a claim: the 69 blocks of the syntax case — constrained
across mammals and among people, the slice the sweep also read the other way on
40 elements — carry the highest predicted mean of the five cases (−0.42 against
−0.55 relaxed and −0.60 tolerant), and so do the 29 recent-constraint blocks
(−0.42). Those two cases are also the shortest (2.7 and 3.2 kb median against 11
kb), 5 and 3 of them hold a measured element, and nothing here was
pre-registered on cases. The number to grow, again.

### What section 18 changes

- §17's count positive **transfers to a different assay**. It is not a property
  of the episomal reporter: the same feature set separates in-vivo VISTA
  verdicts, adds over composition, adds over conservation on elements *chosen*
  for conservation, and beats a dinucleotide shuffle. Area J now has one
  sequence-to-function statement that has survived being moved.
- It transfers to **this sequence** too: +0.238 on the 227 measured elements
  inside the real unknown, against +0.234 over the whole library.
- It does **not** transfer to two things it was asked about. It says nothing
  about *which tissue* a positive drives beyond what length and GC say. And
  swept blind over whole blocks it cannot separate real sequence from its own
  dinucleotide shuffle, so it cannot be used as a genome-wide activity track,
  which is what "score the unknown with it" would have meant.
- On the real unknown the honest statement is a **measured tie**: at 227
  elements the reporter puts these blocks level with the neutral tier, leaning
  below it in all three cell lines. The peer's per-element sweep result stands
  unchallenged, and the reason it cannot be challenged here is coverage — 161 of
  882 blocks have any measured element at all. The next thing this reading needs
  is not a better model but an assay over those 721 blocks.

### Correction received after this section was written: the neutral tier is not a control

Reading 2 above compares the real unknown with the neutral tier, and on 2026-09-17 genomeos-9c
withdrew the headline that comparison was built to test (4a2c199). The arms do not overlap in the
covariate that decides the outcome: the median distance to a coding TSS is 81 kb for the 69 syntax
blocks, 240 kb for the relaxed case and 418 kb for the neutral tier, against 42 kb for the genome's
scored elements. Matched on length, GC and TSS distance, the syntax elements are ordinary (−0.033,
p 0.67) and the relaxed case reads *above* neutral rather than below (+0.025 matched, p 0.014,
against −0.038 raw). What survives is only that every unknown tier sits far below the genome's own
elements; the ordering inside the unknown space is not established.

What that does to the three readings here, stated rather than rerun:

- **Reading 2(2), the transfer (counts over composition on the same 227 elements), is unaffected.**
  It compares two feature sets on one set of elements, so a difference between groups cannot
  produce it. The +0.238 (K562) and +0.389 (HepG2) stand.
- **Reading 2(1), the measured level against the neutral tier, inherits the confound.** Its numbers
  (0.462, 0.489, 0.472) were already inside their nulls, so nothing was claimed from them, and they
  should not now be read as agreement with a withdrawn headline either. The comparison to make is
  against the genome's own scored elements with length, GC and TSS distance held fixed.
- **Reading 2(3), the predicted level over whole blocks, was matched on length only.** It is a null
  (0.497, and 0.514 against its own shuffle), and a null measured against a control that does not
  hold is weak evidence, not strong: it says the count model finds nothing when swept blind, which
  is the claim made, and it cannot say the blocks are like or unlike the neutral tier.

The lesson is the peer's and it is general: a tier average is not a control. The neutral tier was
built as a background for constraint and is a background for nothing else. The practice this section owes the next
reading is to print the median TSS distance of every group it compares: that one number, absent
here and absent in the withdrawn headline, would have caught the confound a day earlier in both.

## The kilobase that was called a base: both readings of a binned track (2026-09-17)

The bigWig reader counted every base of every bin an interval touches. On Zoonomia phyloP, one value
per base, that is exact, and the budget and every constraint figure taken from it are unaffected. On
gnomAD Gnocchi, one value per kilobase, it is the share of touched kilobases — a defensible reading
of a measurement with no finer resolution, and the one `aggregate()` documented — but it was also
being summed into absolute fields, and there it is wrong. A 350 bp element reported `"bases": 1000`,
because a cCRE is shorter than a Gnocchi bin: the field was measuring bins and calling them bases.
genomeos-9c found it while reading this lane's results and quantified it from the committed files
without spending a range read; the fix, the test and the measurement below are this lane's.

`IntervalStats` now carries both readings from one pass. `bases`, `total` and `above` count touched
bins, unchanged. `overlap_bases`, `overlap_total` and `overlap_above` count only the bases inside the
interval, crediting each bin the part of it the interval covers, so they can never exceed the
interval's own length. On a per-base track the two are identical, which is asserted rather than
assumed (`tests/test_attribution.py`), as is the binned case where they diverge.

**The size of it, on the chromosome named before the numbers were read** (chr21,
`scripts/bigwig_bin_shift.py`, `bigwig_bin_shift`, one read of the track for both columns):

| intervals | n | own length | touched bins | own bases | overstated | constrained, touched | constrained, own | fraction, touched → own |
|---|---|---|---|---|---|---|---|---|
| UNKNOWN blocks | 310 | 9,514,145 | 7,444,000 | 7,205,452 | ×1.033 | 557,000 | 508,189 | 0.0748 → 0.0705 |
| registry elements | 9,758 | 2,678,600 | 11,645,000 | 2,520,003 | **×4.62** | 3,202,000 | 683,878 | 0.275 → 0.2714 |

Every element is shorter than a bin, so every element-level base count was inflated about fourfold.
Ratios survive, as the arithmetic says they must: both numerator and denominator were scaled by the
same bin width, so `fraction_above` moves only where an edge bin is weighted differently, 0.0748 to
0.0705 and 0.275 to 0.2714.

**What moved in `variation_chr21` when the absolutes were repointed** (`measured_bp` and
`human_constrained_bp` now count the interval's own bases; the old quantity is kept as
`touched_bin_bp` under its own name):

| tier | measured_bp before → after | human_constrained_bp before → after | fraction before → after |
|---|---|---|---|
| regulatory | 2,654,000 → 2,531,725 | 448,989 → 409,940 | 0.1692 → 0.1619 |
| fossil | 2,661,000 → 2,610,966 | 43,999 → 39,671 | 0.0165 → 0.0152 |
| constrained_unknown | 121,000 → 116,760 | 1,999 → 2,000 | 0.0165 → 0.0171 |
| neutral | 1,980,000 → 1,919,278 | 40,990 → 36,270 | 0.0207 → 0.0189 |
| structural | 4,000 → 3,401 | 2,000 → 1,986 | 0.5 → 0.5839 |

No case, no mammal fraction, no mean and no `fraction_above` changed by more than the edge-bin
weighting, so nothing this lane concluded rests on the corrected fields; what was wrong is exactly
the set of figures nobody had drawn a conclusion from. **The other 23 chromosomes' `variation_chr*`
results still carry the old semantics in those two fields** until the lane is re-run; the fields are
named differently now, so a reader can tell which is which, and the re-run is a range-read cost to
schedule rather than a correction to hide.

The lesson, which is the fourth this week of the same family: a quantity measured at one resolution
and reported at another needs two names, not one. The reader had documented the behaviour in a
docstring and summed it into a field anyway, so the documentation was true and the number was still
wrong — a name is enforceable where a docstring is not.

## 19. Phylogenetic profiling at the gene level: the null decides everything (2026-09-17)

The library-level profiling of 2026-09-12 left the gene level open with a condition attached — a
shuffled null before any pair is called a dependency. `scripts/profiling_genes.py` runs it over the
computed membership of all 42 libraries, 19,391 genes with a 355-species presence bitmask, no requests
and no network. It reports **two** nulls, and they disagree so completely that the result is about the
nulls rather than about the libraries.

**Null one, prevalence-preserving shuffle.** Each gene keeps the *number* of species it is present in
and loses *which* ones. Every library beats it, median excess **+0.3067** mean pair Jaccard, not one of
42 with a single draw above the observed value.

**Null two, shuffled library membership.** Draw a random gene set of the same size, each drawn gene
matched on prevalence to one of the members — **real genes, with real phylogenetic profiles**, so the
tree structure that makes species non-independent is present in the null as well as in the library.

| against matched gene sets (20 draws) | libraries |
|---|---|
| clear of the null (0 draws at or above) | **21** |
| marginal (1–9 draws above) | 7 |
| **at or below the null (≥10 draws above)** | **14** |

**The fourteen that fail are the ancient core**: `core.translation`, `core.energy`,
`core.transcription`, `core.splicing`, `core.dna_repair`, `core.chromatin_epigenome`,
`core.folding_proteostasis`, `core.membrane_transport`, `systems.immune`, `systems.endocrine`,
`systems.metabolic_homeostasis`, `timer.circadian`, `timer.telomere`,
`timer.senescence_checkpoint`. Their genes are in nearly every species, and a gene in nearly every
species co-occurs with any other such gene. Nothing about the library produces that.

**The ones that survive are patchy and lineage-restricted**: `timer.developmental_timing` (+0.0413),
`systems.extracellular_matrix_adhesion` (+0.0393), `blueprint.segmentation` and
`timer.segmentation_clock` (+0.0381 each). Co-occurrence carries information exactly where presence
varies, which is the only place it can.

**So the prevalence shuffle measures phylogenetic non-independence and calls it co-evolution.** It is
kept in the result beside the other because the contrast is the lesson: a null that holds the obvious
covariate fixed can still be beaten by an artefact, and the way to tell is to build the null out of
the same kind of thing you are testing. Real genes, not shuffled bits.

**What this still cannot do.** Even the matched-gene null draws genes independently of the tree, so a
clade-wide loss counts once in the library and once in each draw but the *correlation structure*
between the drawn genes is whatever the genome supplies rather than what a model of the tree would
say. Felsenstein contrasts or a birth-death model along a species tree is what would settle a pair,
and this project holds no tree. **Nothing here calls a pair a dependency**; it says which libraries
would survive the weakest honest test, and half of them do not.
