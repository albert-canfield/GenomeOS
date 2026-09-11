# The genome's own grammar: how BioLang reads DNA

You asked the right question: raw DNA has its own way of saying where a
block starts, what it needs, and where it ends. BioLang should be built on
that, not on C++ or JavaScript habits. This document records what the data
says the delimiters are, how crisp they are, and the design that follows.

## 1. What the data says (chromosome 21, 221 canonical protein-coding transcripts, 1,822 introns)

| Biological "tag" | Sequence | How crisp |
|---|---|---|
| Block start (translation) | `ATG` in Kozak context `(A/G)nnATGG` | 220 of 221 starts are ATG; A at −3 in 55%, G at +4 in most |
| Sub-block close (exon end → intron) | `GT` plus context `(C/A)AG|GTAAGT` | GT in 99.1%, GC in 0.8%, AT (minor spliceosome) 0.2% |
| Sub-block open (intron → exon) | polypyrimidine tract then `AG` : `(T)n…C AG|` | AG in 99.8%; T-rich for ~14 nt before |
| Block end (translation) | `TGA` 60%, `TAA` 24%, `TAG` 15% | 220 of 221 |
| Message end | polyA signal `AATAAA`/`ATTAAA` in the last 40 nt | 66% (the rest use variants or lie further upstream) |
| Header (transcription start) | CpG island around the TSS | 58% of promoters; the rest are TATA/other |

How separable is a delimiter on its own? Scoring chromosome 21 with the
learned splice-donor matrix: a threshold that keeps 90% of true donors also
fires about 7 times per kilobase of random sequence, and the median true
donor scores only 0.57 of the matrix maximum. The signal is real, learnable
and found again in raw DNA (all three first donors of APP are recovered), but
alone it cannot say where a block is. That is the whole design problem.

Every one of these is **a probability, not a token**. That is the central
difference from HTML: `<div>` is either there or not; a splice donor is a
score. The cell's machinery reads DNA statistically, and so must we.

Three more properties no programming language has:

- **Overlap.** Genes sit on both strands, inside other genes' introns, and one
  block can be read in several frames or spliced several ways.
- **Distance.** The "requirements" of a block (which transcription factors
  must be present) are written in the promoter and in enhancers that may be
  a million bases away, on either strand.
- **Context.** The same block reads differently in different cell types
  because the readers (transcription factors, chromatin) differ, not the text.

## 2. The header is a list of requirements

A promoter is literally the dependency list of a gene: each transcription
factor binding site names a protein that must be present for the block to
run. So the natural BioLang reading of a promoter is

```
gene NKX2-5 {
  requires: GATA4, TBX5, MEF2      # from motif hits in the promoter/enhancers
  context: heart_field             # where those factors are present
  ...
}
```

and `requires` is *derived from the sequence* (motif scan + evidence), not
typed by hand. The cell type supplies the environment in which the
requirements are satisfied. This is the biological form of "imports".

## 3. Design: a grammar of signals, learned from data

BioLang v0.3 gets a sequence layer with two constructs and one rule:

```
signal splice_donor    { learned: chr21 GENCODE 50; consensus: MAG|GTRAGT;   examples: 1822 }
signal splice_acceptor { learned: chr21 GENCODE 50; consensus: (Y)n…CAG|;   examples: 1822 }
signal start_kozak     { learned: chr21 GENCODE 50; consensus: RnnATGG;      examples: 221 }
signal stop            { tokens: TGA, TAA, TAG }
signal polyA           { tokens: AATAAA, ATTAAA; window: last 40 nt }

segment coding {
  open:  start_kozak
  body:  (exon splice_donor intron splice_acceptor)* exon
  close: stop
}
segment transcript {
  open:  promoter            # CpG island or TATA; carries `requires`
  body:  coding | noncoding
  close: polyA
}
```

The rule: **a segment is never asserted, only scored.** The parser is a
probabilistic gene finder (the same idea as HMM gene finders and, now, deep
models). Its output is BioIR `Region`/`Gene` entities with `predicted`
evidence naming the signals that fired and a confidence from their scores.
Where no signal fires, the region stays `UNKNOWN`. Where annotation exists
(GENCODE), it wins with `curated` evidence, and the learned grammar is
tested against it.

## 4. What exists now

- `genomeos/genome/signals.py` learns position weight matrices for the
  donor, acceptor and Kozak signals from an annotation, records the
  dinucleotide, stop and polyA statistics above, saves them as a result, and
  scans raw sequence with them on both strands.
- `genomeos signals learn` / `genomeos signals scan` expose it.
- The test checks that the learned consensus is the known biology (GT…AG,
  ATG) and that scanning a raw window of APP finds its annotated donors.

## 5. What comes next

1. ✅ (v1, 2026-09-11) A segment parser that chains signal hits into
   candidate genes and scores them: `genomeos segments --chrom chr21`
   (`genomeos/genome/segments.py`). Coding potential is a codon log-odds
   table learned from the chromosome's canonical CDS; the parser is a
   Viterbi-style dynamic programme over candidate starts, donors, acceptors
   and stops (prefix sums for exon scores, a next-in-frame-stop table, a
   sliding-window maximum for introns; 46 Mb in about two minutes). Scored
   against every coding transcript of GENCODE on chr21:

   | level | sensitivity | precision |
   |---|---|---|
   | gene (any overlap, same strand) | 89.1% | 21.9% |
   | splice site (exact position) | 7.5% | 6.9% |
   | CDS segment (exact both ends) | 4.7% | 5.2% |

   Tried and dropped the same day: banded exon and intron length models
   learned from the annotation (a sliding-window maximum per intron band).
   They moved nothing (exact exons 4.7% / 5.2% before and after). An in-frame
   hexamer coding model instead of codon usage gained a little on a 2 Mb
   slice (site precision 8% → 11%) and nothing over the whole chromosome
   (4.9% / 5.2%, 186 more candidates), so it was dropped too. The bottleneck
   is what three matrices and a codon table cannot see: the scores would need
   a proper probabilistic model of the whole structure (emission and duration
   models of a gene-finding HMM) or conservation, not tuning.

   That is the honest state of "signals alone": the parser lands on nine in
   ten coding genes but draws their exon boundaries right only rarely, and
   it invents 556 candidates for 221 genes. The scoring is the limit,
   not the search: three PWMs and codon usage carry no exon or intron
   length model, no hexamer frame model, no promoter or polyA signal, no
   conservation. Each of those is a measurable step up on this table. Every
   prediction is saved as `predicted` evidence next to this measurement
   (`data/results/segments_chr21.json`); none enters the block map.
   **Second opinion on the splice sites (feature c, 2026-09-11).** The
   bottleneck named above was tested directly. `genomeos segments --chrom
   chr21 --predicted-sites` keeps the parser exactly as it is (Kozak matrix
   for starts, codon log-odds, Viterbi, the same priors) and replaces only
   the donor and acceptor candidates: instead of the two learned matrices it
   takes AlphaGenome's splice-site tracks, four probabilities per base (donor
   and acceptor on each strand), 45 requests of 1 Mb for the chromosome in
   33 s, cached locally (`genomeos/predict/splice_sites.py`). Calibrated
   against GENCODE first: the donor track peaks on the last exonic base and
   the acceptor track on the first exonic base, on both strands, with mean
   probability 0.91 to 1.00 at canonical exon boundaries. A probability
   enters the Viterbi as log2(p/(1−p)) + 10 bits, so a confident site sits
   where a perfect matrix hit would.

   | level | matrices alone | AlphaGenome sites | + starts at ENCODE promoters |
   |---|---|---|---|
   | CDS segment (exact both ends) | 4.7% / 5.2% | **40.1% / 46.2%** | 30.2% / **61.9%** |
   | canonical CDS segment found (1,907 segments) | 8.7% | **74.8%** | 56.3% |
   | splice site (exact position) | 7.5% / 6.9% | **51.6% / 49.2%** | 38.9% / **66.0%** |
   | gene (any overlap, same strand) | 89.1% / 21.9% | 84.2% / 20.0% | 59.3% / **53.8%** |
   | candidates | 556 | 1,003 | 158 |

   (`data/results/segments_chr21_predicted_sites.json` and
   `segments_chr21_predicted_sites_anchored.json`; sensitivity / precision.)
   The third column adds the one local, curated signal for where a gene may
   begin: `--promoter-anchored` confines start codons to 5 kb downstream of
   an ENCODE promoter-like element (4.2% of the chromosome; both strands,
   since an element does not say which strand it serves). It is a trade, not
   a free gain: candidates fall from 1,003 to 158 and precision doubles or
   triples at every level, while sensitivity drops to the genes whose
   promoter the registry has marked (59% of coding genes). Swept: 2 kb
   confines harder (gene 53.8% / 57.1%), 20 kb looser (68.3% / 37.7%); with
   the matrices alone the windows help precision just as much (gene 67.9% /
   65.5% at 5 kb) but exons stay at 3.5% / 9.1%, because the windows say
   where a gene starts and nothing about where its exons end. Nine times the exact exons from the same parser: the
   delimiters were the limit, as section 1 said. What did not move is as
   telling: gene precision stays at one in five because the start codons and
   the coding model are unchanged, and the 1,003 candidates are mostly
   single-exon open reading frames the Kozak matrix lets through. Thresholds
   were swept on the cached tracks (p ≥ 0.05, 0.2, 0.5; 10 or 13.3 bits):
   the lowest threshold with the 10-bit scale is best on every axis, so
   discarding weak sites loses true exons faster than it removes false ones.
   The remaining gap to the annotation is first and last exons (start and
   stop, no splice signal) and alternative isoforms, since every transcript's
   CDS segments count as truth: scored against the canonical transcript's
   segments alone (second row), three in four are found exactly, so most of
   the "missing" exons are isoform-specific ones the parser was never going
   to produce from one structure per locus. The next measurable steps are therefore a
   predicted start signal of the same quality, and scoring against canonical
   transcripts separately. All of it stays `predicted` evidence; the grammar
   runs unchanged without the key.
2. Promoter motif scanning with JASPAR matrices to derive `requires` lists
   with evidence, and to give `cell_type ... expresses:` a sequence-level
   justification.
3. Deep-model signals (AlphaGenome, Evo 2) for the delimiters no PWM can
   capture (enhancers, chromatin accessibility), entering as `predicted`.
   Splice sites: done above (feature c). Enhancer targets: docs/ALPHAGENOME.md
   feature b.
4. Grammar written in BioLang itself, so a bio-engineer can add a signal
   the way one adds a rule: with its examples and evidence.
