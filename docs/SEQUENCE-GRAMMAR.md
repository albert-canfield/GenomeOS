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

1. A segment parser that chains signal hits into candidate transcripts and
   scores them (Viterbi over the grammar), evaluated against GENCODE as
   sensitivity/precision per signal, per chromosome.
2. Promoter motif scanning with JASPAR matrices to derive `requires` lists
   with evidence, and to give `cell_type ... expresses:` a sequence-level
   justification.
3. Deep-model signals (AlphaGenome, Evo 2) for the delimiters no PWM can
   capture (enhancers, chromatin accessibility), entering as `predicted`.
4. Grammar written in BioLang itself, so a bio-engineer can add a signal
   the way one adds a rule: with its examples and evidence.
