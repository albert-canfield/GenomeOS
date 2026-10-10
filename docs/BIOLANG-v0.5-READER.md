# BioLang v0.5: the reader

Adds to v0.4 ([BIOLANG-v0.4-ECONOMY.md](BIOLANG-v0.4-ECONOMY.md)) one block kind, `reader`, and one
rule key, `through`. With them, chromatin state per cell type can change what a program does in a
cell. Every v0.4 program compiles and runs unchanged: a rule without `through` is selected exactly as
before, and a program without `reader` blocks has none. Engine only (Apache-2.0): the parser, BioIR and
rule selection in `genomeos/lang/` and `genomeos/ir/`, nothing from the application.
`docs/BIOLANG-GRAMMAR.md` is the generated reference; this file is the prose.

## Why

Area B, Next 1 (ROADMAP §3 B; §5, M2). Every cell carries the same genome and differs in which
parts it can read (NODES-READER-WRITER.md §2; v0.4 §7.2, "specialisation is state"). The project
decodes that per cell type: open, poised or silent, on every chromosome. Until v0.5 a program could
record that reading only as a rule's `context_evidence`, which gates nothing and which the engine does
not read. That is still true of `context_evidence`. The reader is the construct that does gate.

## Syntax

```
reader K562 {                       # header: one declared cell_type id
  open: HS2, HS3                    # elements this cell type can read
  poised: E7                        # stated poised by the measurement
  closed: E9                        # elements it cannot read
  unknown: E12                      # not read: a value, neither open nor closed
  evidence: experimental "DNase-seq, the source file"; confidence: 0.8
}

rule GATA1 activates HBG2 {
  through: HS2                      # the one element this rule acts through
  strength: 0.9; threshold: 5; hill: 2
  evidence: ...
}
```

- A `reader` names one declared `cell_type`, and each list names declared `element` or `region`
  blocks. Several readers may describe one cell type (one per assay, each with its own evidence and
  confidence, like every block). They must not give one element two states.
- `through` names exactly one element, on an `activates` or `inhibits` rule. A rule that acts through
  two elements is written as two rules.

## Semantics: selection, before the network is built

`Module.active_rules(context)` is where a cell's rules are chosen, and every runtime that builds a
network from it inherits the gate. The order is `when`, then `expresses` (both unchanged), then the
reader:

| the cell's reader states the element | the rule |
|---|---|
| `open` | runs as written |
| `poised` | held |
| `closed` | held |
| `unknown`, stated | held |
| nothing: no reader for the cell type, or a reader that does not mention the element | `unknown`, held |
| the context names no `cell_type` at all | `unknown`, held |

A context names its cell type by id or by `name`, as `expresses` already accepts. Readers are not
inherited through `parent`: v0.4 §7.2 has no `extends`, and a parent's reading is not a child's.

**Held means strength 0, not removed.** A held rule stays in the cell's network as a copy at strength
0, so it contributes nothing. It is not dropped, because the network runtime transcribes a gene that
has no activator rule at all at its full `max_rate`. Dropping a closed enhancer's rule would therefore
drive its target to maximum, the opposite of closed. Measured on `data/demo/reader_gating.bio`: with
ENH1 closed, the held rule leaves TGT.mRNA at its basal 0.05, and the same rule dropped would leave it
at 60.05, the reading for ENH1 open. Under the runtime's default `mean` combination a held activator
still counts in the gene's denominator with a zero term. That is how the runtime already reads an
activator whose level is zero. For an `inhibits` rule, holding and dropping are the same.

**UNKNOWN is held, and is said.** An element whose state is unknown in a cell is neither open nor
closed. The engine never treats it as open. It holds the rule, as v0.4 does for a `when` clause on an
unrecorded context, which runs nowhere rather than everywhere (R1). `Module.reader_gate(context)`
returns every gated rule with its state and the reason, and `bio check FILE --context cell_type=X`
prints the held rules split by state (`poised n, closed n, unknown n`), so a rule held for lack of a
reading is never counted as closed.

## What the parser refuses

- a `reader` whose header is not one id, a reader that states no element, or an element stated
  twice in one reader (in one list or in two);
- a `reader` of an undeclared `cell_type`, or one that names something that is not a declared
  `element` or `region`;
- two readers of one cell type that give one element different states (the engine does not choose
  between measurements);
- `through` naming more than one element, naming an undeclared one, or sitting on a `produces`,
  `binds`, `degrades` or `modifies` rule, whose effect no strength carries;
- `through` in a located program (one that declares a compartment). The located runtime chooses its
  rules from `when` alone and does not read the reader, so a held rule would run there.

## What it does not do

- It does not say an element works. A reader states accessibility and marks. An open element is one a
  rule may act through in that cell, not evidence that it does, and a held rule is not shown to be
  inactive.
- It does not read `context_evidence`, which stays a recorded measurement that gates nothing.
- It gates nothing in the Body's per-cell network. `Body` builds its network with no context, so every
  rule with `through` reads `unknown` there and is held. Per-cell contexts in the Body are later work.
- It does not change `when`, `expresses`, or any rule without `through`.

## BioIR

`Reader` (`cell_type`, `states`: element to state, `evidence`, `confidence`), `Rule.through` (default
`""`), and `Module.readers`. BioIR JSON writes them as `"readers"` and `"through"`, and a file written
before v0.5 reads back with none. `bioir_version` is still written as `0.4`. Whether this construct
makes the language v0.5 in the version number is one of the open decisions below.

## Tooling

`bio test` lines can be judged in a context, `# test: TGT.mRNA final > 50 in cell_type=OpenCell`,
with one run per distinct context. Two facts count the gate: `active_rules` (the rules in that cell's
network, held ones included) and `held_rules` (the ones its reader holds). `data/demo/reader_gating.bio`
runs the same two rules in five cell types: open, open with the silencer open too, closed, poised,
and a cell no reader describes.

## What the application half needs (not built here)

`attribution/compile.py` would emit, per compiled program:

1. a `cell_type` block for every reader cell type it uses (`K562`, `HepG2`, ...);
2. one `reader <CellType>` block per cell type, listing the program's own `element` ids under the
   state the decoded chromatin gives them, with the reader's source files as evidence. An element
   the reader did not cover goes under `unknown` or is left out; both hold;
3. `through: <element id>` on each enhancer-to-gene rule that acts through an element.

It must keep its own vocabulary apart. `context_evidence`'s `open_in_reader` / `not_open_in_reader` /
`not_assessable` are a rule's record of one reading. A reader's `open` / `poised` / `closed` /
`unknown` are what the engine gates on. Mapping one onto the other is an application decision, and
`not_assessable` maps to `unknown`, never to `closed`.

## Open decisions for Albert

1. **UNKNOWN: hold, or refuse.** Implemented: the most conservative choice that still runs. A rule
   whose element is unknown in the cell is held and reported, never run. The alternative is to refuse
   the run, as an unbound UNKNOWN parameter does. That would make any program unrunnable in a cell its
   reader does not fully cover.
2. **Poised.** Implemented: poised holds a rule, following the project's reader, which counts a
   poised promoter as not read. A program cannot yet say that a rule runs on poised.
3. **The version number.** The grammar and BioIR still say 0.4. §10 decision 1 of v0.4 asks whether
   v0.5 makes locations mandatory, so the number was not moved without that decision.
