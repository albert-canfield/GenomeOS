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
