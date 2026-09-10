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
