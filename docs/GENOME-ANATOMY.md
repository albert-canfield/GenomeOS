# Genome anatomy: how DNA is organised, and what building one from zero means

`genomeos anatomy` counts the blocks and elements of a chromosome from its
sequence and annotation. Run on three very different genomes it shows the
organisation principles directly. Numbers below are from
`data/results/anatomy_comparison.json` (GENCODE 50, Ensembl WBcel235.63).

## Three genomes side by side

| measure | human mtDNA | human chr21 | C. elegans III |
|---|---|---|---|
| length | 16.6 kb | 46.7 Mb (14% gaps) | 13.8 Mb |
| coding genes | 13 (+22 tRNA, 2 rRNA) | 221 (+590 lncRNA, 209 pseudogenes…) | 2,640 (+1,219 non-coding) |
| coding genes per Mb | 785 | 5.5 | 192 |
| bases in protein code | 68.5% | 0.8% | 25.8% |
| bases in introns | 0% | 49.8% | 39.0% |
| bases intergenic | 7.3% | 30.6% | 27.8% |
| other exonic (UTR, ncRNA) | 24.2% | 4.6% | 7.4% |
| exons per transcript (median) | 1 | 6 | 5 |
| intron median | none | 1,953 bp | 95 bp |
| single-exon coding genes | 100% | 25% (the KRTAP keratin-associated cluster is single-exon) | 2.5% |
| coding transcript models per gene (GENCODE/Ensembl) | 1 | 13.1 | 1.6 |
| neighbour orientation | 10 tandem, 1 conv, 1 div (one strand mostly) | 116 / 52 / 52 | 1,394 / 622 / 623 (balanced) |
| overlapping coding pairs | 3 (ATP8/ATP6, ND4L/ND4, …) | 39 | 440 |
| gene median length | 784 bp | 32.0 kb | 2.3 kb |
| intergenic spacing median | 96 bp | 23.0 kb | 718 bp |
| protein median length | 261 aa | 317 aa | 334 aa |
| GC | 44% | 41% (23–60% in 100 kb bins) | 36% (32–39%, flat) |
| CpG islands | 0 | 573 | 220 (the worm does not methylate CpG, so islands are not promoters there) |
| homopolymer runs ≥12 | 0 | 7,919 | 1,145 |

## The whole human genome (GRCh38, all 25 chromosomes, streamed and discarded)

Result: `data/results/anatomy_hg38_by_chromosome.json`. Each chromosome was
downloaded, counted and deleted; total disk use never exceeded one chromosome.

| measure | value |
|---|---|
| bases | 3.088 Gb (4.9% assembly gaps) |
| coding genes | 20,107 (78,733 genes of all types) |
| protein-coding bases | 37.4 Mb = **1.21%** |
| intron and non-coding gene bodies | 1,892 Mb = **61.3%** |
| intergenic | 861 Mb = 27.9% |
| other exonic (UTR, ncRNA exons) | 147 Mb = 4.8% |
| exons per canonical transcript | median 8 |
| intron median by chromosome | 842 bp (chr19) to 2,683 bp |
| coding-gene density | chr19 25.3 per Mb and chr17 14.3 the densest; chr13, chr18 and chrY under 3.5 |
| CpG islands (≥400 bp) | 39,410 |

So the human design spends roughly one base in eighty on protein, six in ten
on the inside of genes that is spliced away, and three in ten between genes.
The gene-dense chromosomes (19, 17) have the shortest introns; the sparse ones
the longest: density and intron size move together, which says that the
budget is one number (bases per gene) partitioned between introns and
spacing.

## What the comparison says

1. **There are three ways to organise a genome, and they are budgets.** The
   mitochondrion is a compact design: almost no intergenic space, no introns,
   genes butted together and even overlapping, one strand carrying most of
   them, every gene a single exon. The worm is a dense eukaryotic design:
   a quarter of the bases code, introns are short (95 bp) and numerous, genes
   sit 192 to the megabase in both orientations. The human chromosome is a
   sparse design: under 1% codes, half the bases are intron, genes are far
   apart (23 kb median spacing) and each carries a dozen annotated transcript models.
2. **Density trades against regulation and variety.** The human genome spends
   its bases on control (regulatory regions in the intergenic and intronic
   space) and on alternative transcripts; the mitochondrion spends nothing on
   either and is read as two long polycistronic transcripts cut into genes.
3. **Exon count is conserved, intron length is not.** Median exons per
   transcript are 5–6 in both worm and human; the difference is intron size
   (95 bp versus 2 kb). Exons are the units; introns are the spacing.
4. **Orientation is balanced in complex genomes.** Worm and human place
   neighbours tandem, convergent and divergent in the ratio a random coin
   would give; the mitochondrion does not. Divergent pairs share a promoter
   region, a design pattern worth reusing.
5. **Sequence-level elements track the organisation.** CpG islands mark
   promoters in mammals (573 on chr21, 58% of protein-coding promoters), not
   in the worm; homopolymer and microsatellite runs scale with intergenic
   and intronic space.

## If you build a genome from zero

A design is a budget and a layout. The inventories give the numbers:

| design choice | compact (mtDNA-like) | dense (worm-like) | sparse (human-like) |
|---|---|---|---|
| bases per coding gene (chromosome ÷ genes) | ~1.3 kb | ~5 kb | ~180 kb (non-gap) |
| introns per transcript | 0 | 4, ~100 bp | 5, ~2 kb, some >100 kb |
| intergenic spacing (median) | ~100 bp | ~700 bp | ~23 kb |
| transcript models per gene | 1 | 1–2 | ~13 |
| where regulation lives | almost nowhere (one control region) | short promoters | promoters + distant enhancers + CpG islands |
| what it can do | run one organelle | build a 959-cell animal | build a 37-trillion-cell animal with 200+ cell types |

So the first design decision is not which genes, but which organisation: the
compact plan when the parts list is small and control is external (the
mitochondrion is controlled by the nucleus), the dense plan for a complete
organism with modest cell-type variety, the sparse plan when many cell types
must read the same genes differently. Everything else (which libraries, in
docs/GENOME-AS-CODE.md) is filled into that budget.

What the anatomy tool cannot yet count, and the next steps: regulatory
elements (needs motif scanning and predicted-evidence models), transposable
elements and repeat families (needs a repeat library), replication origins,
and the 3-D folding that sets which enhancer reaches which promoter.
