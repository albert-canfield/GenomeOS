# BioLang v0.1

BioLang is the source language of GenomeOS: an organised way to write biology
as code. It is declarative. It states what exists, what acts on what, under
which conditions, and how sure we are. The runtime decides what happens.

## Why a language rather than a database

A database stores facts. A language lets you compose, version, review, diff,
test and reuse them. The goal is that a bio-engineering team works the way a
software team does: modules, imports, code review, tests, continuous
integration against evidence.

## Grammar

```
module <dotted.name>
import <dotted.name>

gene <Id> { <properties> }
protein <Id> { <properties> }
region <Id> { <properties> }
rule <Source> <action> <Target> { <properties> }
param <name> = <number> [unit] { <properties> }
```

`<action>` is one of `activates inhibits modifies produces binds degrades`.

Properties are `key: value`, one per line or separated by `;`. A whole block
may sit on one line. `#` starts a comment.

| Key | Applies to | Meaning |
|---|---|---|
| `evidence` | all | `<kind> "<source>" [note]`, kind in experimental / curated / predicted / inferred / none |
| `confidence` | all | 0..1 |
| `symbol` | gene | HGNC-style symbol |
| `locus` | gene, region | `chr7:1000000-1000500(-)` |
| `produces` | gene | protein id; emits a `produces` rule |
| `basal` | gene | leaky transcription rate |
| `max` | gene | maximal transcription rate |
| `half_life` | protein | hours |
| `sequence` | protein | amino-acid string |
| `role` | region | free text or `unknown` |
| `strength`, `threshold`, `hill` | rule | Hill-function parameters |
| `when` | rule | `key = value, key = value`; `any` matches all |
| `id` | rule | override the generated id |

## Example

`data/demo/repressilator.bio` implements the Elowitz & Leibler 2000
oscillator: three genes each repressing the next. Run it:

```
genomeos check data/demo/repressilator.bio
genomeos run   data/demo/repressilator.bio --hours 60 --init TetR=10
```

## Compile-time guarantees

- undeclared references are errors, so a rule cannot point at a phantom gene
- unknown evidence kinds are errors, so provenance cannot be faked by typo
- confidence outside 0..1 is an error

## Roadmap for the language

v0.2: nested blocks (`transcript` inside `gene`), `cell_type` and `tissue`
declarations, `event` blocks (divide, differentiate, die) with guards,
imports resolved across files, a standard library prelude
(`bio.std.central_dogma`, `bio.std.cell_cycle`).

v0.3: `organism` and `experiment` blocks (see README examples), units checked
at compile time, and a formatter/linter so BioLang files review cleanly.
