# Cancer: healthy cell versus tumour cell

A tumour is the same genome with somatic changes. GenomeOS therefore treats
cancer analysis as a comparison: the person's normal sample against their
tumour sample, the differences annotated with the variant engine and graded
against distilled population knowledge.

## Knowledge source and the API decision

[cBioPortal](https://www.cbioportal.org) hosts hundreds of curated tumour
sequencing studies. Its public REST API (`/api`) needs no key, answers in
milliseconds, and is what GenomeOS uses (`genomeos/cancer/cbioportal.py`,
standard library only). The MCP server at mcp.cbioportal.org requires
authentication and an MCP client; it is meant for conversational agents, not
for a reproducible pipeline, so it is not used.

Following stream-distil-discard, `genomeos cancer distil` pulls the mutations
of a 110-gene driver panel across MSK-IMPACT 2017 (10,945 tumours, 58 cancer
types) in about 30 seconds and keeps only `data/results/cancer_msk_impact_2017.json`:
per gene, the fraction of tumours mutated, the fraction per cancer type, and
the recurrent protein changes (hotspots). The top of that table:

| gene | tumours mutated | top change |
|---|---|---|
| TP53 | 41.5% | R175H, R248Q, R273H … |
| KRAS | 15.0% | G12D 444, G12V 371, G12C 262 |
| TERT | 13.3% | promoter |
| PIK3CA | 12.4% | H1047R, E545K |
| APC | 10.2% | truncations |

KRAS is mutated in most pancreatic tumours in that study; the per-type table
is what lets a sample's driver profile point at a cancer type.

## The comparison

```
genomeos cancer compare --normal blood.vcf --tumour tumour.vcf --chrom chr21 --packet packet.json
```

1. **Somatic set**: variants in the tumour and not in the normal.
2. **Consequence**: each is classified against the gene models (nonsense,
   missense, frameshift, splice, …) by the same engine validated on ClinVar.
3. **Grading**: driver frequency and hotspot match from the distilled
   knowledge, combined with consequence severity into a score. Evidence lines
   name the study and the numbers.
4. **Cancer type**: enrichment of the sample's mutated drivers per cancer
   type, reported as a suggestion (`inferred`), never a diagnosis.
5. **Targets**: altered genes whose products sit on the plasma membrane or
   cell surface (GO cellular component), the only place a designed binder can
   reach from outside the cell.
6. **Agent packet**: a JSON with the ranked variants, suggested types, surface
   targets, hard constraints (use only the inputs, grade every claim, not
   clinical advice), a minimal prompt, and the expected output schema, so a
   focused AI agent can propose a binder target and payload from evidence.

The test builds a RUNX1 truncation on chromosome 21 in a synthetic tumour
VCF, removes a shared germline variant, and checks that it ranks first with
its cBioPortal evidence, that EGFR (and not TP53 or RUNX1) is a surface target,
and that the packet carries the constraints.

## About the "custom protein connector"

The story you describe is the logic of targeted therapy: find a molecule on
the surface of the tumour cell that healthy cells lack or show far less of,
design a protein that binds it (antibodies, nanobodies, and now de-novo
binders designed with tools such as RFdiffusion), and attach a payload: a
toxin (antibody-drug conjugates), a radionuclide, an immune recruiter
(bispecifics), an engineered immune cell (CAR-T, where the "connector" is
on the T cell), or a nucleic-acid cargo. The payload is what kills; the
binder is what makes it selective. GenomeOS's contribution is the front of
that pipeline: from two genomes to a graded list of what is different and
what is reachable. Binder design and anything involving a patient is outside
this project and belongs to laboratories and clinicians.

## Tumour alone: no matched normal

`genomeos cancer tumour --vcf T.vcf --packet packet.json` takes the tumour's
DNA by itself. What changes without the normal sample, and how each gap is
handled:

| step | with a normal | tumour alone | evidence |
|---|---|---|---|
| somatic vs germline | subtraction | gnomAD population frequency ≥ 1% → likely inherited, set aside; rare germline variants cannot be told apart | inferred |
| consequence | local gene models (chr21, chrM) | Ensembl VEP REST on any chromosome: consequence, gene, HGVS c./p., SIFT, PolyPhen, COSMIC ids, gnomAD; cached per variant in `data/knowledge/vep/` | curated / predicted |
| drivers | cBioPortal frequencies and hotspots | same | curated |
| burden | – | coding somatic variants per Mb, assuming a whole exome; the assumption is written into the result | inferred |
| protein | – | mutant peptide window around each top missense change, from the UniProt canonical sequence (neoantigen candidates; HLA binding not predicted) | derived |
| pathways | – | a truncated driver is treated as absent and run through its Reactome pathways: reactions lost, most affected pathway | curated + inferred |

The score adds severity, driver frequency, hotspot, COSMIC presence and
predicted damage; likely-germline variants are kept in the packet but
scaled down and listed separately so the agent does not build on them. The
packet asks for driver events, cancer type, a surface target with payload,
neoantigen candidates, risks and next experiments, each with confidence.
The Cancer tab has a "Tumour alone" form for the same run.

## From alterations to targets

`surface_targets()` here answers one narrow question: which altered genes carry
a Gene Ontology plasma-membrane or cell-surface annotation. That is a first
filter, not a target assessment: it cannot tell a receptor's ectodomain from a
kinase held against the inner leaflet, it says nothing about healthy tissue, and
a mutated gene is not a surface target.

The therapeutic pipeline answers the rest. `genomeos therapeutic --tumour T.vcf`
takes the same ranked variants and works out where each protein sits and whether
a binder can physically reach it, how the tumour differs from healthy tissue,
what happens after binding, whether an intracellular mutation could still be
seen through HLA, and which therapeutic mechanism the biology supports, with the
evidence level and the missing data for each. See
[docs/THERAPEUTICS.md](THERAPEUTICS.md). `genomeos cancer tumour --therapeutic`
runs it from this command.

## The other two ways a gene breaks

A VCF of point mutations carries one of the three ways a tumour breaks a gene.
The other two are read from the same cBioPortal client
(`genomeos cancer alterations`, 110 genes in about a minute,
`data/results/cancer_alterations_msk_impact_2017.json`, 139 KB kept):

| gene | mutated | amplified | deep-deleted | where the frequency lives |
|---|---|---|---|---|
| CDKN2A | 4.3% | 0.06% | **7.6%** | Glioma 32.6%, GIST 21.2% |
| ERBB2 | 3.0% | **4.0%** | 0.01% | Esophagogastric 22.6%, Breast 14.1% |
| CCND1 | 0.5% | **4.3%** | 0.03% | |
| MYC | 0.7% | **4.0%** | 0.02% | |
| PTEN | 6.1% | 0.02% | **2.5%** | |

CCND1 and MYC are passengers by mutation frequency and drivers by copy
number. Structural variants come back with their partners: EML4-ALK in 42
tumours, KIF5B-RET in 15, CCDC6-RET in 10.

The discrete-copy-number endpoint is used rather than the molecular-data one,
because the latter returns a row per gene per sample and 95% of those rows say
"diploid". A panel study calls only what is on its panel and only the deep
events: MSK-IMPACT makes no shallow gain or loss call at all, and the table
records that absence as an absence of a call rather than as diploid.

**These events are alterations, not annotations.** A gene reaches the tumour
comparison and the target ranking through a copy-number or structural call
exactly as it does through a variant (`genomeos/cancer/alterations.py`):
a typed record, the cohort frequency, a score on the same scale as the variant
score, and an evidence line per claim. Supply them with `--cnv` and `--sv` to
`genomeos cancer tumour` or `genomeos therapeutic`.

Two readings are refused. A copy-number table of absolute copies and one of
GISTIC discrete calls disagree about the number 2 — an amplification in one
convention, an untouched gene in the other — so the format is decided by a
stated rule (a negative value can only be a call, a value above 2 can only be a
count) and an ambiguous table is read the way that invents no amplification;
`--cna-format` overrides. And a deep deletion is never read as a reason to aim
at the gene: a homozygously deleted gene makes no product, so it is dropped
from the surface targets, classed `unsuitable`, and every therapeutic mechanism
that has to recognise a product fails on it with that sentence. What such a
finding is worth is the dependency the loss creates, which GenomeOS does not
yet model.

## The `cancer.*` library layer

The BioLib catalogue names the toolkits a healthy genome uses. This layer
names the ones a tumour breaks, and it is the only one whose membership is a
frequency rather than a function: a gene is in `cancer.amplified` because
tumours amplify it, not because of what its product does. Membership is
computed from the two committed tables and never written down
(`genomeos/cancer/libraries.py`, `genomeos cancer libraries`), so a library
cannot drift from its evidence and every member carries the number that put
it there.

| library | genes | at or above | largest members |
|---|---|---|---|
| `cancer.mutated` | 87 | 1% | TP53, KRAS, TERT, PIK3CA, APC |
| `cancer.hotspots` | 8 | 1% | TERT, KRAS, PIK3CA, BRAF, TP53 |
| `cancer.amplified` | 12 | 1% | CCND1, MYC, ERBB2, EGFR, CDK4 |
| `cancer.deleted` | 3 | 1% | CDKN2A, PTEN, RB1 |
| `cancer.rearranged` | 25 | 0.1% | EGFR, ALK, BRAF, ROS1, FGFR2 |

`genomeos cancer libraries --gene CDKN2A` gives the question the layer exists
for — every way the study records one gene breaking, with the cancer types
that carry it:

```
deep_deletion  7.62%  member  [Glioma 32.6%, GIST 21.2%, Melanoma 18.4%]
mutation       4.31%  member  [Skin, non-melanoma 20.3%, Melanoma 14.2%]
hotspot        0.38%  below threshold; changes R80* x42, R58* x34, H83Y x32
fusion         0.23%  member  partners CDKN2B-AS1 x8, MIR548H2 x2, MTAP x2
```

Three things the layer refuses to say. A hotspot gets no per-cancer-type
breakdown, because the distillation counts recurrent changes study-wide and
printing the gene's *mutation* distribution beside a hotspot frequency would
read as the hotspot's. A gene off the panel is reported as never looked at
rather than as unbroken: `breaks_in("CD19")` answers `on_panel: false` and
says that silence is not evidence of a healthy gene. And membership means a
gene is selected, not that it drives, and certainly not that it can be
treated.

The layer does not join the shared `LIBRARIES` catalogue on import —
`genomeos/lib` is another area's file and a layer appearing there as a side
effect of importing `genomeos.cancer` would be a surprise. `register(LIBRARIES,
LAYERS)` wires it in one call, whenever that area wants it.

## Next

- Expression: tumour-versus-normal RNA to find surface proteins that are
  over-expressed rather than mutated, the more common target class. The
  therapeutic pipeline already accepts patient RNA (`--rna`); what is missing is
  a cohort reference to compare it against.
- What a deep deletion is actually worth: the dependency the loss creates
  (MTAP with CDKN2A is in the fusion partners above), which nothing here
  models.
