# BioIR v0.1 — Biological Intermediate Representation

BioIR is what BioLang compiles to and what BioVM executes. It is deliberately
independent of any biological database so the compiler can improve while the
runtime stays stable. The reference implementation is `genomeos/ir/model.py`;
the JSON produced by `genomeos compile` is the interchange form.

## Types

### Evidence
```
kind:       experimental | curated | predicted | inferred | none
source:     free text (citation, database id, model name)
organism:   default "Homo sapiens"
note:       free text
```

### Entity (base)
```
id, kind, attrs{}, evidence, confidence
```
Subtypes:

| Type | Extra fields | Notes |
|---|---|---|
| `Gene` | symbol, locus, transcripts[], basal_rate, attrs.max_rate | locus optional in v0.1 |
| `Transcript` | gene_id, exons[Locus], cds | `splice()` builds mRNA from exons |
| `Protein` | sequence, half_life_h (float or UNKNOWN) | |
| `Region` | locus, role (string or UNKNOWN) | the honest placeholder |

### Rule
```
id
source      entity id
action      activates | inhibits | modifies | produces | binds | degrades
target      entity id
strength    0..1
threshold   Hill half-max (K)
hill        Hill coefficient (n)
when        {key: value}; value "any" matches everything
evidence, confidence
```
A rule applies in a context only if every `when` key matches.

### Parameter
```
name, value, unit, evidence, confidence
```
Parameters override runtime defaults (`translation_rate`, `mrna_half_life`,
`protein_half_life`, `noise`) and let a module carry its own calibration.

### Module
```
name, imports[], entities{id: Entity}, rules[], parameters{name: Parameter}
```
Invariants enforced at compile time:
- entity ids are unique within a module
- every rule's source and target are declared entities
- confidence is within 0..1
- evidence kinds are from the closed set above

## JSON form

```json
{
  "bioir_version": "0.1",
  "name": "synthetic.repressilator",
  "imports": [],
  "entities": [{"__type__": "Gene", "id": "tetR", ...}],
  "rules": [{"id": "LacI inhibits tetR", "action": "inhibits", ...}],
  "parameters": [{"name": "translation_rate", "value": 5.0, ...}]
}
```
`UNKNOWN` serialises as the string `"UNKNOWN"`; `Locus` as
`chrom:start-end(strand)`.

## What v0.2 adds (planned)

- `Complex`, `Reaction`, `Signal`, `Receptor` entities for pathway import from
  Reactome and SBML
- `CellType`, `Tissue`, `Event` entities for the cell and developmental layers
- provenance versioning (`source_version`, `compiled_at`)
- a `Model` evidence subtype naming the AI model and its score for
  `predicted` rules (AlphaGenome, Evo 2, STATE)
- cross-module references (`import` resolution) and namespacing

## Unstated confidence (2026-09-28, item 12 S1)

A block that leaves out `confidence:` has not stated one. The parser has read that as `UNSTATED`
since 262ea4e: a float subclass equal to 0.0 in every calculation, told apart from a stated 0.0 by
`confidence_stated()`. A stated `confidence: 0.0` exists on purpose (six parameters in
`bio.std.methylation`, two compiled regions), so the value alone cannot carry the difference. The
second external review found that BioIR JSON lost it: `to_dict` wrote both as 0.0 and `from_dict`
read both back as stated. This section is the census of every path a confidence travels, taken
before any change, then the representation that keeps the difference on each of them.

### Census, before the change

Taken on `e704175` with the reader sweep repeated by hand; line numbers are that revision's.

**Parsing.** `lang/parser.py` `_common` gives a missing key `UNSTATED`: kept. `parser.py:457` gave a
nested transcript `cconf or conf`, so a transcript that states 0.0 took its gene's value instead
(the reverse error: a stated 0.0 read as absent); the tree has one transcript block and it states no
0.0, so no figure depends on it. `parser.py:442` gives a gene's `produces` rule the gene's value
object: kept. `parser.py:532` copies an element's value into each `targets[]` entry: kept in memory,
lost through JSON like every other copy.

**IR.** Every one of the 14 `confidence` fields defaulted to a plain 0.0, so an object built without
a confidence read as a stated 0.0. Ten construction sites build one that way: `attribution/bridge.py:702`
(Rule), `organism/digital_development.py:93` (Experiment), `runtime/division.py:96` and
`runtime/located.py:96` (Regime), `runtime/economy.py:151, 200` (Allocation), and
`scripts/abundance_gate.py:228, 229, 239, 240` (Pool, Allocation). `Module.merge`, `copy.deepcopy`,
pickle and the in-memory copy `Module.from_dict(m.to_dict())` (`runtime/located.py:132`) keep the
object.

**BioIR JSON.** `Module.to_dict` (`ir/model.py:947`) passed the object through, and `json.dumps`
writes a float subclass as its value: 0.0. `Module.from_dict` (`:986`) built a plain float from it,
stated. Every JSON module writer goes through these two: `genomeos compile` (`cli.py:98`), `bio compile`
(`lang/tools.py:84`), the web editor's `/api/compile` module (`web/server.py:503`), BioForge's
`--output` (`cli.py:680`), the REPL's `:show` (`bio.py:297`); and both loaders read a `.json` module
through `from_dict` (`cli.py:89`, `lang/tools.py:17`).

**Other JSON writers of an IR confidence.** `/api/compile`'s `report.confidence` and
`weak_rules[].confidence` (`server.py:508, 515`); the Cell view's `rules[]`, graph edges and
`cell_type` (`server.py:1484, 1523, 1578`); the network debugger's `explain[]` (`server.py:1110`);
an organism experiment's result, `ExperimentResult.to_dict` (`organism/experiment.py:118`), which
`genomeos experiment --json` writes (`cli.py:740`); the uncertainty report,
`UncertaintyReport.to_dict` (`runtime/uncertainty.py:76`), which `genomeos grow` and the web grow,
segmentation and gastrulation endpoints return. All wrote an unstated confidence as 0.0.

**BioLang writers.** No writer turns a `Module` into BioLang text. The one that starts from an IR
object is BioForge's `DesignResult.to_bio` (`organism/forge.py:219`); it writes no `confidence:` key,
so its text parses as unstated, but the Experiment it returns carried `UNSTATED_CONFIDENCE`, a plain
0.0 (`forge.py:65, 209`), so the object and its own text disagreed. Every other `confidence:` writer
writes a value from code or from a result file, never an IR object's: `organism/reference.py:283`,
`organism/fate_rules.py:413`, `organism/human.py:285, 297`, `molecules/compiler.py:578`,
`attribution/compile.py:436, 527, 539`, `organism/tf_atlas.py:209`, `lib/proteome.py:145`,
`knowledge/ligand_receptor.py:141`, `organism/haematopoiesis.py`, `scripts/located_mitochondrion.py:163`.
The compiled chromosome programs leave the key out of every predicted fact (R4), which is the
BioLang form of unstated and already round-trips.

**CSV and the Evidence explorer.** `evidence.py` rows carry `stated`, and `--csv` writes `unstated`
in the confidence column (262ea4e): kept. `scripts/compile_genome_programs.py:70` and
`scripts/measured_layer.py:186` count through `stated` and keep a labelled legacy tally.

**Summaries that averaged an unstated confidence as 0.0.** `Module.confidence_report`
(`ir/model.py:935`), printed by `genomeos check` (`cli.py:125`), `bio check` (`lang/tools.py:53`), the
REPL's `:check` (`bio.py:337`) and the web editor; the REPL's rule list (`bio.py:340`, printed 0.0);
the `# test: confidence` line (`bio.py:193`, mean over all rules); the runtime uncertainty report
(`runtime/uncertainty.py:47-70`, fed by `report_for_network`, `report_for_ageing` and
`runtime/body.py:1238-1245`); the debugger's trace line (`runtime/debugger.py:65`, printed
`conf=0.00`); the located runtime's chain score (`runtime/located.py:505`), which made an unstated
link the weakest at 0 and is printed by `bio run` (`lang/tools.py:236`); and the methylation
program's assumption note (`runtime/methylation.py:266`, "confidence 0").

**Thresholds.** The weak-rule lists at `< 0.5` in `genomeos check` (`cli.py:135`), `bio check`
(`lang/tools.py:65`) and `/api/compile` (`server.py:520`) listed every unstated rule as weak. The
methylation flag at `< 0.4` (`methylation.py:263`) calls an unstated parameter an assumption, which
it is; only its note was wrong. No sort orders by an IR confidence.

**Copies that keep the object.** `min(x.confidence, 0.3)` in `forge/design.py:59-66` and
`organism/forge.py:131, 135` returns `UNSTATED` itself. Derived scores (`located.py:509`,
`uncertainty.py:49`) are plain floats and are fixed where they are computed.

**The test that accepted the loss.** `tests/test_evidence.py::
test_the_distinction_lives_in_the_parsed_ir_and_not_in_bioir_json` (262ea4e) asserted that a
reloaded unstated confidence reads as stated.

**Stored files.** No BioIR JSON file is committed or on disk under the checkout or `data/` (searched
for `bioir_version`, depth 6), and no file in `data/results` stores an uncertainty report. No stored
result can move.

**Measured on the tree at the census.** The 66 programs the Evidence explorer reads (demos,
organisms, `bio.std`, the 24 compiled chromosomes) carry 919,533 unstated confidences on entities,
rules and parameters, in 26 of the programs; after `to_dict`, `json` and `from_dict`, 0 of them are
still unstated. Their `check` weak-rule list held 445,072 rules below 0.5, of which 4,695 state a
confidence and 440,377 state none.

The negatives are `tests/test_unstated_confidence.py` (engine) and
`tests/test_unstated_confidence_writers.py` (application), committed with this census as strict
expected failures: all 18 fail on the revision they were committed to, each on the loss it names.

### The representation, registered before the build

| Where | Unstated | Stated 0.0 |
|---|---|---|
| BioLang text | key left out (unchanged) | `confidence: 0.0` |
| In memory | `UNSTATED`, and the default of every IR `confidence` field | `0.0` |
| BioIR JSON | `"confidence": null` (`UNSTATED_JSON`), file marked `"records_unstated_confidence": true` (`RECORDS_UNSTATED`) | `0.0` |
| Other JSON exports | `null`, through `confidence_to_json` | `0.0` |
| CSV (Evidence explorer) | `unstated` (unchanged) | `0.0` |
| Means | left out, counted beside; `None` where nothing is stated | in the mean |
| Weak lists | never weak; counted apart | weak below the line |
| Weakest-link scores | left out and named; the score is then an upper bound | in the minimum |

`confidence_from_json` is the inverse: null reads back as `UNSTATED`, a number as the value it
states. `from_dict` applies it to every block and to an element's `targets[]`.

**Why null, and not a second key or an absent key.** One key cannot contradict itself; a pair such
as `"confidence": 0.7, "confidence_stated": false` can, and every reader would have to consult both.
A reader that computes with null fails, where a reader that ignores a flag, or reads an absent key
with `.get("confidence", 0.0)`, pools the unstated value with a stated 0.0 without a word, which is
the defect itself. Null is also how the project already writes a number nobody has: a `Certainty`'s
probability, the explorer's mean over nothing stated, an uncertainty level with nothing behind it.
And `to_dict` goes on writing every dataclass field, the list the generated grammar documents. The
cost is JavaScript, where `null < 0.5` is true: the census found no script comparing a BioIR JSON
confidence, and each web view that shows one now tests for null.

**Files written before the mark.** A BioIR file without `records_unstated_confidence` was written
when a stated 0.0 and an unstated confidence were both written as 0.0. Nothing in such a file can say
which a given 0.0 was, and `from_dict` does not guess: it reads the 0.0 the file carries, as stated.
The program the file came from still has the difference; compile it again. No such file exists in
the checkout today.

**No version change.** `bioir_version` stays 0.4: it is the language's version (`lang/grammar.py`),
and the language did not change. The mark carries the change in the JSON form.

**Engine.** Everything above lives in `genomeos/ir`, `genomeos/lang`, `genomeos/runtime` and
`genomeos/bio.py`; none of it imports the application.

### Built (same day), and what moved

Every path in the census now keeps the difference: `to_dict` and `from_dict` as registered; the 14
IR defaults are `UNSTATED`; the transcript takes its gene's confidence only when it states none;
`confidence_report` averages stated values and `confidence_counts` gives both counts;
`lang.tools.confidence_lines` is the one place `genomeos check` and `bio check` print them from; the
web compile report sends `confidence_counts` and `unstated_rules`; the Cell view, the debugger,
`ExperimentResult.to_dict` and BioForge's `UNSTATED_CONFIDENCE` write or carry `UNSTATED` as null;
the uncertainty report keeps `items` as the count of everything behind a level and adds `unstated`,
the mean running over the rest; the located runtime names unstated links in `unstated`; the trace
line, the REPL and the methylation note print "unstated". The 18 negatives pass, and
`tests/test_evidence.py` holds the preservation test that replaces the one that accepted the loss.

Measured before and after on the same tree (66 programs, 42 network programs run through
`report_for_network`, 8 organism programs grown):

| Figure | Before | After |
|---|---|---|
| Evidence explorer, whole tree: facts / stated weak / stated strong / unstated | 986,947 / 19,299 / 48,114 / 919,534 | unchanged |
| Evidence explorer CSV, bytes (with and without the compiled programs) | 243,226,498 and 4,451,303 | unchanged |
| Unstated confidences still unstated after BioIR JSON | 0 of 919,533 | 919,533 of 919,533 |
| Rules in the check reports' weak list (below 0.5) | 445,072 | 4,695 (the 440,377 that state none are counted apart) |
| Check means, 24 compiled chromosomes: `regulatory_element` | 0.026 to 0.041 | 0.750 to 0.784 |
| the same: `rule` | 0.000 to 0.001 | 0.900 in 22; none stated in 2 |
| the same: `gene`, `domain` | 0.00 | none stated |
| `data/organisms/human/noncoding_chr21.bio`: `regulatory_element`; `gene`, `domain` | 0.028; 0.00 | 0.764; none stated |
| `data/demo/repressilator.bio`: `region` | 0.00 | none stated |
| Uncertainty reports, 42 network and 8 organism programs: confidence, label and items per level | | unchanged |

Every moved figure is a live report; none is stored in a result file.
