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
| ENCODE and Roadmap Epigenomics, reference epigenomes via the ENCODE portal | histone ChIP-seq (H3K4me3, H3K27ac, H3K4me1, H3K27me3, H3K9me3; replicated peaks and fold change over control) and whole-genome bisulfite sequencing for the reader's eleven biosamples, the epigenome layer (`genome/epigenome.py`); experiments in order H3K4me3, H3K27ac, H3K4me1, H3K27me3, H3K9me3: K562 ENCSR000AKU, ENCSR000AKP, ENCSR000EWC, ENCSR000EWB, ENCSR000APE; WGBS ENCSR765JPC; HepG2 ENCSR000AMP, ENCSR000AMO, ENCSR000APV, ENCSR000AOL, ENCSR000ATD; WGBS ENCSR786DCL; GM12878 ENCSR000AKA, ENCSR000AKC, ENCSR000AKF, ENCSR000DRX, ENCSR000AOX; WGBS ENCSR890UQO; H1 ENCSR814XPE, ENCSR880SUY, ENCSR000ANA, ENCSR186OBR, ENCSR000APZ; WGBS ENCSR617FKV; IMR-90 ENCSR087PFU, ENCSR002YRE, ENCSR831JSP, ENCSR431UUY, ENCSR055ZZY; WGBS ENCSR888FON; SK-N-SH ENCSR975GZA, ENCSR000FCU, ENCSR661BMA, ENCSR914QOK, ENCSR657OGA; WGBS ENCSR145HNT; cardiac muscle cell ENCSR652QNW, ENCSR000NPF, ENCSR276OLB, ENCSR864LRY, ENCSR269ULZ; keratinocyte ENCSR000ALO, ENCSR000ALK, ENCSR000ALI, ENCSR000ALL, ENCSR000ARN; hepatocyte ENCSR442ZOI, ENCSR507UDH, ENCSR689QUB, ENCSR637RLN, ENCSR758GMX; WGBS ENCSR351IPU; astrocyte ENCSR000AOU, ENCSR000AOQ, ENCSR000AOT, ENCSR000AOR, ENCSR000AQR; CD14-positive monocyte ENCSR000ASN, ENCSR000ASJ, ENCSR000ASM, ENCSR000ASK, ENCSR000ASP; WGBS ENCSR017BUL | free to use, ENCODE data-use policy; cite the ENCODE Project Consortium 2020 (Nature 583:699) and the Roadmap Epigenomics Consortium 2015 (Nature 518:317) |
| GTEx | tissue expression | open access summary statistics |
| 1000 Genomes Project phase 3 (through the Ensembl REST API) | population allele frequencies of the known-locus benchmark's named variants | open access; Ensembl EMBL-EBI terms |
| NHGRI-EBI GWAS Catalog | lead variants and mapped genes, read as looked-up evidence in the known-locus benchmark | EMBL-EBI terms, free to use |
| Human Protein Atlas | healthy-tissue RNA levels, localisation | CC BY-SA 4.0 |
| Open Targets Platform | tractability, approved drugs, safety liabilities | CC0 1.0 |
| cBioPortal | tumour alterations, cohort expression, driver frequencies | ODbL portal, per-study terms |
| ClinVar, dbSNP, gnomAD, COSMIC, PubMed | variant interpretation and citations | public; COSMIC non-commercial terms apply to COSMIC ids |
| Genome in a Bottle (GIAB) / HG002 | the test human, open consent | public domain, PGP open consent |
| CellPhoneDB | ligand-receptor pairs | MIT |
| WormWeb / Sulston 1983, Sulston & Horvitz 1977 | the complete C. elegans lineage | CC BY 1.0 (WormWeb), publications cited |
| Packer et al. 2019 (GEO GSE126954) | worm single-cell lineage annotation | public GEO record |
| Ma et al. 2021 (Zenodo 4737593) | protein levels of 266 transcription factors per lineaged C. elegans embryonic cell over time: the worm's reader, terminal-fate rules and the commitment tests | CC BY 4.0 |
| Sender & Milo 2021 | human cell counts, lifespans, turnover | publication cited |
| MitoCarta3.0 (Rath et al. 2021, Nucleic Acids Res 49:D1541), with King & Attardi 1989, Wiedemann & Pfanner 2017 and Alberts et al., *Molecular Biology of the Cell* Table 12-1 | the mitochondrial inventory, sub-compartments, OXPHOS subunits and TargetP signals behind the located OXPHOS program and the BioLang v0.4 stage 1 gate; the rho0 phenotype, import routes and compartment volumes it is checked against | freely available from the Broad Institute; publications cited |
| Horvath 2013, Hannum 2013 | epigenetic clock coefficients | publications cited |
| Hart et al. 2017 (CEGv2) | core-essential gene set | publication cited |
| BioModels | SBML reference models | CC0 |
| Zoonomia (Christmas et al. 2023), UCSC phyloP and phastCons tracks | constraint across 241 placental mammals per base, conserved elements | publication cited; UCSC tracks free for all uses |
| gnomAD genomic constraint (Chen et al. 2024, "Gnocchi") | constraint within 76,156 human genomes per kilobase | gnomAD open terms, publication cited |
| VISTA Enhancer Browser | transgenic mouse enhancer outcomes, the measured ground truth for elements | free for academic and commercial use, publication cited |
| Ensembl Compara | orthologues, paralogues, gene trees, pairwise genome alignments | Apache 2.0 service, EMBL-EBI terms |
| UCSC genomicSuperDups (Bailey and Eichler) | segmental duplications | free for all uses; publication cited |
| JASPAR 2026 | transcription-factor binding profiles, TFClass families | CC BY 4.0 |
| Kircher et al. 2019 (GEO GSE126550) | saturation-mutagenesis MPRA: which bases of 21 regulatory elements change activity | CC BY 4.0 publication; data portal code GPL-3 |
| Human Pangenome Reference Consortium, release 1 (Liao et al. 2023), UCSC `hprc90way`, `hprcArrV1`; release 2 v2.1 structural variants (UCSC `hprc2v21Sv`) | 89 human haplotype assemblies aligned to hg38 by Cactus (presence, alleles, insertions per base), arrangements against hg38, structural variants from 233 assemblies | HPRC open access, publication cited; UCSC tracks free for all uses |
| gnomAD v4.1.1 genomes and v4.1 structural variants (UCSC bigBed tracks) | allele frequencies beside the panel's value domains, each with its allele number | gnomAD open terms |
| TRExplorer v2 tandem-repeat catalogue (UCSC `trexplorer` track) | tandem-repeat loci with HPRC and TenK10K allele-size histograms | as distributed through UCSC |
| ENCODE UW Repli-seq (wavelet-smoothed, 11 cell lines, hg19) and UCSC `hg38ToHg19` liftOver chain | replication timing per kilobase, the mutation-rate control of the human panel | ENCODE data-use policy; UCSC chains free for all uses |
| GTEx v8 fine-mapped eQTLs (DAP-G, UCSC `gtexEqtlDapg`) and MPRAVarDB (UCSC `mpraVarDb`) | read-out evidence attached to storage units: eQTL variants and reporter-assay allele pairs | GTEx open access; MPRAVarDB as distributed through UCSC |

## Use and licence

GenomeOS is open source, written by Albert Canfield, under two licences split
by path: the **BioLang engine** is Apache 2.0 and the **GenomeOS application**
is AGPL-3.0-or-later. [LICENSING.md](LICENSING.md) has the map and the
reasoning; [LICENSE](LICENSE) and [LICENSE-APACHE](LICENSE-APACHE) hold the
texts.

The licence of the code is a separate question from the terms of the data.
Several sources GenomeOS reads are free only for non-commercial or academic
use whatever the code says: **AlphaGenome** (API and weights, non-commercial),
**NetMHCpan** (DTU academic licence, never shipped or downloaded by GenomeOS),
**COSMIC** identifiers and some **cBioPortal** studies. A commercial user must
check every source they enable, and GenomeOS names the source of every fact
precisely so that check is possible.

GenomeOS produces research hypotheses with their evidence and their gaps. It
is not a medical device, it gives no clinical advice, and it never claims a
therapy will work.
