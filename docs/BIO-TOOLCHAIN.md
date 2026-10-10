# `bio`: the BioLang toolchain on its own

**Status on 2026-10-10.** Still describes the code: `bio` (`genomeos/bio.py`, a console script in
`pyproject.toml`) has the five verbs below, the REPL commands are as written, and `bio --version`
prints the application's 1.0.0 (the packaged engine prints its own, 0.1.0). Since this was written:
`bio run` also runs located (v0.4) programs and takes `--knockout`; a `# test:` line can name
`P@Compartment` or a methylation context, can be judged `in cell_type=X` and count `active_rules` or
`held_rules` there (v0.5 reader), and its confidence mean counts only rules that state one; and the
separate package exists (`scripts/package_engine.py`, ARCHITECTURE.md §10.1 and §10.2). Over the
committed `data/demo` and `genomeos/std` programs, `bio test` ran 20 files and 63 checks on
2026-10-10, all passing.

The roadmap's 1.x goal is BioLang as an engine like Node: write, check,
compile, run and test biology without the rest of GenomeOS. `bio` is that
entry point, installed next to `genomeos`:

```
bio check   FILE [--context k=v]     compile, resolve references, report evidence and confidence
bio compile FILE [-o out.json]       BioLang → BioIR JSON (the contract between language and VM)
bio run     FILE [--hours H] ...     the network runtime; .xml runs SBML, .bnet a Boolean network
bio test    PATH...                  every .bio file under PATH, judged by its `# test:` lines
bio repl                             type BioLang, run it, inspect it
```

## Tests live in the program

A `# test:` comment is a claim the module makes about itself, checked by
`bio test` after a 48 h run (or without a run, for facts about the module):

```
# test: rules >= 3             module facts: rules, entities, events, unknowns
# test: confidence >= 0.6      mean confidence over the rules
# test: TetR peaks >= 1        a species: final, peaks, min, max
```

An organism program (BioLang v0.3, stages declared) is grown by the Body
runtime to the start of its last stage; the program's own `assert:` lines
become claims, and `# test: alive == 961`, `cells`, `deaths` read the
summary at that point.

Every run also checks that no level went NaN or infinite. The demo
repressilator carries three such lines; `bio test data/demo genomeos/std`
runs eight files and 11 checks. *(2026-10-10: 20 committed files and 63 checks, all passing.)*

## The REPL

Blocks are accepted as they compile. A block that refers to something not
declared yet (a gene that `produces: A` before `protein A` exists) waits,
and folds in when the declaration arrives, so a program can be typed in
any order. `:run 12 TetR=10` runs the module for 12 hours with an initial
level; `:check` prints the confidence report and every rule with its
evidence; `:show` prints the BioIR; `:load FILE` reads a file in.

## What stays separable

`genomeos/bio.py` imports the language (`genomeos.lang`), the IR
(`genomeos.ir`) and the runtime (`genomeos.runtime`); none of those import
the genome, knowledge or web layers. Packaging them as `biolang` with this
entry point is a move, not a rewrite (docs/ARCHITECTURE.md). *(2026-10-10: done since 2026-09-15:
`scripts/package_engine.py` generates the package, ARCHITECTURE.md §10.1.)*
