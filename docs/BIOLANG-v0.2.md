# BioLang v0.2

Adds to v0.1 (see [BIOLANG-v0.1.md](BIOLANG-v0.1.md)): nested blocks, cell
types, events, imports and a standard prelude. All v0.1 files compile
unchanged.

## New constructs

```
import bio.std.cell_types          # genomeos/std/cell_types.bio
import bio.std.ageing              # genomeos/std/ageing.bio
import organs/heart.bio            # a file relative to this one

gene NKX2-5 {
  max: 50; basal: 0.1; produces: Nkx2-5
  transcript NKX2-5-201 {           # nested
    exons: chr5:173232109-173233366, chr5:173234555-173235310
    cds:   chr5:173232700-173233366, chr5:173234555-173235000
  }
  evidence: curated "GENCODE 50"; confidence: 0.9
}

cell_type Cardiomyocyte2 {
  parent: Cardiomyocyte            # must be declared or imported
  expresses: NKX2-5, MYH6          # genes not listed are silenced in this context
  ontology: CL:0000746
  evidence: curated "Cell Ontology"; confidence: 0.8
}

event divide {
  rate: 0.5 /yr
  when: cell_type = any
  effect: telomere_bp -= 70 bp     # several effect lines allowed
  effect: divisions += 1
  evidence: experimental "Harley 1990"; confidence: 0.8
}
```

### Added 2026-09-11: compiled protein properties and regulatory elements

A `protein` block may carry what the federated compiler knows
(`genomeos protein TP53 --bio` writes one):

```
protein TP53 {
  accession: P04637;
  sequence: MEEPQ…;
  isoforms: P04637-1, P04637-2;
  domains: p53_DNA-bd, p53_TAD2;
  pathways: R-HSA-69541, R-HSA-6804756;
  interactions: MDM2, CREBBP;          # partners with STRING experimental support
  structures: 1TUP, AF-P04637;         # PDB ids; AF- marks the AlphaFold model
  evidence: curated "UniProtKB/Swiss-Prot; InterPro; Reactome; PDB; STRING physical channel";
  confidence: 0.5
}
```

An `element` block is a regulatory element with the genes it reaches and
the node (CTCF domain) that bounds it:

```
element APP_enh1 {
  class: enhancer;                     # promoter | enhancer | insulator | open_chromatin
  locus: chr21:25900000-25900400(+);
  targets: APP, CYYR1;
  domain: chr21:D89;
  evidence: inferred "same CTCF domain"; confidence: 0.4
}
```

Both compile to BioIR entities (`Protein` with the new fields,
`RegulatoryElement`); `bio check` reports them with the rest.

## Semantics

- **Imports** merge the imported module's entities, rules, events and
  parameters into the importing module. Identifier collisions are errors.
  `bio.std.*` resolves inside the package; other dotted names resolve
  relative to the importing file; literal `.bio` paths are allowed.
- **Cell-type context.** `Module.active_rules({"cell_type": X})` returns rules
  whose `when` matches and whose source and target genes are expressed in `X`.
  `Module.silenced_genes(context)` lists what the context turns off. A cell type
  with no `expresses` list silences nothing.
- **Events** are declarative: rate with unit, guards (`when`), and effects
  (`var op value unit`, op in `+= -= *= =`). The cell runtime consumes the
  same facts; `bio.std.ageing` states them with their evidence.
- **Transcripts** nest in genes and carry exons and CDS segments as loci, so
  `genomeos annotate` output and hand-written modules share one shape.

## Standard prelude

| Module | Contents |
|---|---|
| `bio.std.cell_types` | Cell, StemCell, Fibroblast, Neuron, Cardiomyocyte, HematopoieticStemCell with CL ids |
| `bio.std.ageing` | divide, mutate, senesce events and telomere-at-birth parameter, all cited |

## Check it

```
genomeos check data/demo/cell_context.bio --context cell_type=Cardiomyocyte2
genomeos check data/demo/cell_context.bio --context cell_type=Neuron2
```
