# Acknowledgements

GenomeOS stands on public science and open-source software. Nothing here is
our data: every biological fact comes from a curated source that is named in
the output, and every optional engine is somebody else's work. This page
records what we use, under which licence, so a user knows exactly what they
are running and what they are allowed to do with it.

## The core

The GenomeOS core has **no runtime dependencies**: the parser, BioIR, the
engines, the genome reader and the compilers are written against the Python
standard library alone. That is a deliberate constraint (docs/DECISIONS.md
D1), and it is why a `.bio` program runs anywhere Python does.

Development uses **pytest** (MIT) and **ruff** (MIT). Packaging uses
**hatchling** (MIT) and **uv** (Apache 2.0 / MIT).

## Optional engines and models

Each is an extra, never a dependency; without it the feature is listed and
disabled (docs/DECISIONS.md D34).

| Package | Used for | Licence |
|---|---|---|
| [AlphaGenome](https://github.com/google-deepmind/alphagenome) (Google DeepMind) | predicted regulatory effect of a variant: expression, splicing, chromatin | client Apache 2.0; **the API and model are free for non-commercial use only** |
| [process-bigraph](https://github.com/vivarium-collective/process-bigraph) (Vivarium) | composition: several engines on one clock | Apache 2.0 |
| [MHCflurry 2](https://github.com/openvax/mhcflurry) (OpenVax) | peptide/HLA class I binding prediction | Apache 2.0 |
| NetMHCpan 4.1 ([DTU Health Tech](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/)) | peptide/HLA class I binding prediction | **academic licence**; not shipped, not downloaded, called only if the user installs it |
| [biolearn](https://bio-learn.github.io/) | epigenetic clocks (our in-house implementation is the default) | MIT |

AlphaGenome's own model release acknowledges the libraries it is built on,
and they carry through to anyone who enables that feature: Abseil, anndata,
Chex, Einshape, Etils, Haiku, huggingface_hub, JAX, jaxtyping, kagglehub,
NumPy, Optax, Orbax, pandas, pyBigWig, pyarrow, pyfaidx, PyRanges,
TensorFlow, typeguard. MHCflurry brings NumPy, pandas, scikit-learn,
mhcgnomes and PyTorch. We thank all their contributors and maintainers.

## Data sources

Every number GenomeOS reports carries its source and evidence level. The
distilled summaries committed under `data/results/` are derived works of
these; raw files are never redistributed (docs/DECISIONS.md D7).

| Source | Used for | Licence / terms |
|---|---|---|
| GENCODE / Ensembl | gene models, transcripts, CDS, sequence, VEP consequences, paralogues | Ensembl: Apache 2.0 service, EMBL-EBI terms; GENCODE data free for any use |
| UCSC Genome Browser | hg38 sequence, RepeatMasker tracks | free for all uses; data from the UCSC API |
| UniProtKB / Swiss-Prot | protein sequences, isoforms, topology, features | CC BY 4.0 |
| InterPro, PDB, AlphaFold DB | domains, experimental and predicted structures | per-source open terms; AlphaFold DB CC BY 4.0 |
| Reactome | pathways, reachability graphs | CC BY 4.0 |
| STRING | protein associations (association, not interaction) | CC BY 4.0 |
| Gene Ontology, Cell Ontology, Uberon | library membership, cell types, anatomy | CC BY 4.0 |
| ENCODE | candidate cis-regulatory elements, DNase per biosample | free to use, ENCODE data-use policy |
| GTEx | tissue expression | open access summary statistics |
| Human Protein Atlas | healthy-tissue RNA levels, localisation | CC BY-SA 4.0 |
| Open Targets Platform | tractability, approved drugs, safety liabilities | CC0 1.0 |
| cBioPortal | tumour alterations, cohort expression, driver frequencies | ODbL portal, per-study terms |
| ClinVar, dbSNP, gnomAD, COSMIC, PubMed | variant interpretation and citations | public; COSMIC non-commercial terms apply to COSMIC ids |
| Genome in a Bottle (GIAB) / HG002 | the test human, open consent | public domain, PGP open consent |
| CellPhoneDB | ligand-receptor pairs | MIT |
| WormWeb / Sulston 1983, Sulston & Horvitz 1977 | the complete C. elegans lineage | CC BY 1.0 (WormWeb), publications cited |
| Packer et al. 2019 (GEO GSE126954) | worm single-cell lineage annotation | public GEO record |
| Sender & Milo 2021 | human cell counts, lifespans, turnover | publication cited |
| Horvath 2013, Hannum 2013 | epigenetic clock coefficients | publications cited |
| Hart et al. 2017 (CEGv2) | core-essential gene set | publication cited |
| BioModels | SBML reference models | CC0 |

## Use and licence

The GenomeOS source code is released under the **MIT licence** (see
[LICENCE](LICENSE)) and its author is Albert Canfield.

**GenomeOS is intended for research and education, not for commercial use.**
Two things sit behind that statement and they are different:

1. **The code.** MIT permits commercial use. If GenomeOS is to be
   non-commercial as a whole, the licence has to say so (a non-commercial
   licence such as PolyForm Noncommercial or CC BY-NC-SA for the
   documentation). This is an open decision recorded in
   docs/DECISIONS.md; until it is settled the code stays MIT and this
   paragraph states the author's intent.
2. **The data and the optional models.** Several sources GenomeOS reads are
   free only for non-commercial or academic use, whatever the code's licence
   says. Enabling AlphaGenome binds the user to Google DeepMind's
   non-commercial terms; NetMHCpan needs a DTU licence; COSMIC identifiers
   and some cBioPortal studies carry their own restrictions. A commercial
   user must check every source they enable, and GenomeOS names the source
   of every fact precisely so that check is possible.

GenomeOS produces research hypotheses with their evidence and their gaps. It
is not a medical device, it gives no clinical advice, and it never claims a
therapy will work.
