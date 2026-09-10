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

## Next

- Copy-number and structural variants from cBioPortal (`_cna`,
  `_structural_variants` profiles) alongside mutations.
- Expression: tumour-versus-normal RNA to find surface proteins that are
  over-expressed rather than mutated, the more common target class.
- A Cancer tab in the web UI: choose two VCFs, see the ranked table, the
  suggested types and the targets, download the agent packet.
- Add the distilled driver table as a library layer (`cancer.*`) with
  evidence, so the libraries know which genes break in which cancers.
