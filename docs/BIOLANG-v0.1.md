# BioLang v0.1

**Status on 2026-10-10.** The grammar, property table, example and compile-time guarantees below
still describe the code: the demo compiles and passes its four `bio test` checks, and `genomeos check`
and `genomeos run` take the flags shown. The current grammar is [BIOLANG-GRAMMAR.md](BIOLANG-GRAMMAR.md),
generated from the parser: 27 block kinds plus the `module` and `import` directives, of which v0.1
introduced five (`gene`, `protein`, `region`, `rule`, `param`). The newest, v0.5's `reader`
([BIOLANG-v0.5-READER.md](BIOLANG-v0.5-READER.md), 2026-10-10), lets a cell type's chromatin state
hold a `rule` that acts `through` an element; the grammar header and `bioir_version` still read 0.4.
The roadmap at the end is history; what was built and what was not is marked beside it.

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

*(2026-10-10: built in v0.2 — nested `transcript`, `cell_type`, `event` blocks with guards and
effects, and imports across files ([BIOLANG-v0.2.md](BIOLANG-v0.2.md)); differentiation and death
became v0.3 `decision` actions. Not built — a `tissue` declaration, and the two prelude modules named
here: the prelude is seven other modules, `bio.std.ageing`, `cell_types`, `development`,
`human_stages`, `human_turnover`, `methylation` and `signalling`.)*

v0.3: `organism` and `experiment` blocks (see README examples), units checked
at compile time, and a formatter/linter so BioLang files review cleanly.

*(2026-10-10: `organism` and `experiment` were built ([BIOLANG-v0.3.md](BIOLANG-v0.3.md)). Units are
checked where the language has a closed set — times on timers, stages and decisions (`min`, `h`,
`d`, `wk`, `yr`), volumes, concentrations and cost units — and, since a later 2026-10-10 change
(ROADMAP §5 item 13, A), an event's rate unit (per time) and a `param`'s unit as well; the sets are in
[BIOLANG-GRAMMAR.md](BIOLANG-GRAMMAR.md). No formatter or linter was built: `bio` has `check`,
`compile`, `run`, `test` and `repl`.)*
