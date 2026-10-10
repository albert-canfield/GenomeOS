# Acknowledgements

**Status on 2026-10-10.** The core still has no runtime dependencies, and the optional
engines and the data sources listed below hold with their licences as stated. The tables are
behind the code, so "full list with licences" means the full list of what this file names,
not of what the project reads: sources in use whose provenance and terms are recorded
elsewhere include DepMap 24Q4 Public, CC BY 4.0 (docs/THERAPEUTICS.md, "Sources"); the CPTAC
phosphoproteomic profiles reached through cBioPortal and the Ochoa et al. 2020 reference
phosphoproteome (docs/PROTEIN.md); and the 4D Nucleome data portal, AlphaMissense, the PGS
Catalog and ENCODE4 lentiMPRA, each with its own section in docs/DATA.md. Folding those into
the tables below is the next edit to this file.
*(Later on 2026-10-10: done; they are the rows marked "row added 2026-10-10" at the end of the
data-source table, with the 2026-10-10 lanes' sources beside them.)*

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
| Ensembl REST variation endpoint | the independent check of every variant position written into the known-locus benchmark | Apache 2.0 service, EMBL-EBI terms |
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
| DepMap 24Q4 Public (Broad Institute; Figshare+, doi:10.25452/figshare.plus.27993248.v1) *(row added 2026-10-10)* | CRISPR gene-effect scores, streamed with only the requested gene columns cached: synthetic lethality and marker selectivity (`therapeutics/depmap.py`; `synthetic_lethality.json`, `marker_selectivity.json`) and K562's dependency row in the cell-level closure test M3 (`attribution/closure_cell.py`; `closure_cell_result.json`) | CC BY 4.0, the licence the Figshare+ record states |
| CPTAC phosphoproteomes through cBioPortal *(row added 2026-10-10)* | relative abundance of curated phosphosites in tumours, per site and per study, no per-patient value (`molecules/ptm.py`, `scripts/cptac_phospho.py`; `phosphosite_tumour_abundance.json`, which carries `luad_cptac_2020`, `gbm_cptac_2021` and `brain_cptac_2020`) | cBioPortal's FAQ: "Unless otherwise noted, data in cBioPortal are available under the ODC Open Database License with no restrictions on the use of the data, as long as you properly give attribution to the original studies", and "There are some studies that restrict the commercial use of the data, but that will be explicitly mentioned in the study information"; the derived summary is offered under ODbL 1.0 (`CPTAC_LICENCE` in `molecules/ptm.py`) |
| Ochoa et al. 2020 (Nat Biotechnol 38:365) reference phosphoproteome, through the authors' funscoR package (`feature_spectral_counts.rda`; raw data PRIDE PXD012174) *(row added 2026-10-10)* | how many cell types or tissues identified each curated phosphosite at 1% site-level FDR, not occupancy (`molecules/ptm.py`, `scripts/phosphosite_observation.py`; `phosphosite_observation.json`) | funscoR's DESCRIPTION: "License: LGPL"; PRIDE under EMBL-EBI's terms, which place "no additional restrictions on the use or redistribution of the data … other than those provided by the original data owners"; publication cited |
| 4D Nucleome data portal: insulation-based boundary calls, and in-situ Hi-C `.hic` matrices read by HTTP range (K562 `4DNFITUOMFUQ`, GM12878 `4DNFI1UEG1HD`) *(row added 2026-10-10)* | measured node edges (`genome/hic.py`; `domains_chr21_hic_*.json`) and measured contact in place of 1/distance (`genome/hic_contact.py`; `crispri_contact.json`); downloads need a free account key (docs/DATA.md) | the portal: "The data at the 4DN Data Portal is generated by the 4DN Network and made freely available to the scientific community"; cite the 4DN White Paper (doi:10.1038/nature23884) and the Data Portal paper (doi:10.1038/s41467-022-29697-4), acknowledge the lab that generated the data, and contact it before a publication |
| AlphaMissense hg38 table (Google DeepMind, public bucket `dm_alphamissense`) *(row added 2026-10-10)* | predicted missense effect for one person's variants (`genome/missense.py`, `genomeos individual coding --predict`), kept under that person's git-ignored directory and never in a committed result (docs/DATA.md) | the two sources disagree, both read 2026-10-10: the table's own header (2023) reads "Licensed under CC BY-NC-SA 4.0 license"; the GitHub README reads "AlphaMissense predictions are licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0)". The code applies the stricter, CC BY-NC-SA 4.0. Code Apache 2.0 |
| PGS Catalog (EMBL-EBI): harmonised GRCh38 scoring files and each score's metadata *(row added 2026-10-10)* | per-person polygenic sums with their coverage (`genome/polygenic.py`, `genomeos individual pgs`), written under the person's directory and never under `data/results` (docs/DATA.md) | the catalog: "Individual PGS obtained from the database should also be cited appropriately, and used in accordance with any licensing restrictions set by the authors", under EMBL-EBI's terms of use; each score's licence is read from its metadata; cite Lambert et al. 2021 and 2024 (Nature Genetics) |
| ENCODE4 lentiMPRA (Nadav Ahituv lab, UCSF): joint library ENCSR106SZM in K562 ENCSR203UFY, HepG2 ENCSR405QCT and WTC11 ENCSR336MKI; element quantifications ENCFF802FUV, ENCFF475FKV, ENCFF769REH *(row added 2026-10-10)* | measured enhancer activity per 200-bp element and cell line (`attribution/mpra.py`; `mpra_chr*.json`, `mpra_genome_wide.json`) | ENCODE: "External data users may freely download, analyze and publish results based on any ENCODE data without restrictions"; acknowledge the Consortium and the production lab and reference the ENCSR and ENCFF accessions |
| ENCODE4 proliferation CRISPR screen (CRISPRi) in K562, ENCSR162HMH (files ENCFF945GIJ, ENCFF995GEZ, ENCFF141MGA, ENCFF100ABB), and its control screen ENCSR625TAF (ENCFF902VHH, ENCFF743IPS, ENCFF482NVE, ENCFF209ITV); Tim Reddy lab, Duke *(row added 2026-10-10)* | the cell-level closure test M3: guide abundance at elements whose attributed gene K562 depends on (`attribution/closure_cell.py`; `closure_cell_result.json`) | ENCODE data-use policy, as for lentiMPRA above |
| ENCODE polyA plus RNA-seq gene quantifications (ENCODE4 uniform pipeline, RSEM, GENCODE V29): one experiment per reader biosample, 13 biosamples, 27 files, accessions in the result *(row added 2026-10-10)* | measured RNA for the reader-gated run M2b (`attribution/readergate.py`; `readergate_rna_sources.json`, `readergate_result.json`) and the multi-element combination forms (`attribution/combine2.py`) | ENCODE data-use policy, as for lentiMPRA above |
| ENCODE knockdown RNA-seq of transcription factors (CRISPRi, siRNA, CRISPR and shRNA RNA-seq; 53 experiments, 50 in K562 and 3 in HepG2, accessions in the result) *(row added 2026-10-10)* | whether a factor's knockdown moves the genes the compiled rules route through its occupied elements, M1b (`attribution/regulators.py`; `regulators_result.json`) | ENCODE data-use policy, as for lentiMPRA above |
| UCSC `encRegTfbsClustered`: ENCODE 3 TF ChIP-seq peak clusters with cell types, hg38 (`encRegTfbsClusteredWithCells.hg38.bed.gz`) *(row added 2026-10-10)* | measured occupancy behind the compiled regulator assignment, M1b (`attribution/regulators.py`; `regulators_census.json`) | the track page: "Users may freely download, analyze and publish results based on any ENCODE data without restrictions", citing the ENCODE Consortium and the production laboratories; UCSC data "freely available for any use" |
| EngreitzLab `CRISPR_comparison` benchmark data (Gschwind et al. 2025), which gathers the CRISPRi enhancer–gene screens of Nasser 2021, Gasperini 2019, Schraivogel 2020 and others *(row added 2026-10-10)* | the measured enhancer–gene pairs of the CRISPRi result and its held-out K562 pairs (`attribution/crispri.py`; `crispri_published_v2.json`, docs/CRISPRI-RESULT.md) | the repository's LICENSE: MIT License, Copyright (c) 2025 Andreas R. Gschwind; the screens' publications cited |
| Gasperini et al. 2019 (Cell), GEO GSE120861: the at-scale gene–gRNA-group pair table *(row added 2026-10-10)* | which public measurement reaches the 882 real-unknown blocks, M5 (`scripts/measured_unknown_inventory.py`; `measured_unknown_inventory.json`) | NCBI: "NCBI itself places no restrictions on the use or distribution of the data contained therein", while "Some submitters of the original data … may claim patent, copyright, or other intellectual property rights"; the record states no licence of its own as the M5 inventory read it; cite the publication |

*(2026-10-10, lane-docfix.)* The rows marked "row added 2026-10-10" fold in the sources the status
block above names and those the 2026-10-10 lanes read. Every quoted term was read that day from
the source's own page or file: the Figshare+ record's licence field, cBioPortal's FAQ, funscoR's
DESCRIPTION, EMBL-EBI's and NCBI's policy pages, the 4DN portal's data-use guidelines, the
AlphaMissense README and the first lines of its hg38 table, the PGS Catalog's about page, ENCODE's
citing page and the `encRegTfbsClustered` track page, and the `CRISPR_comparison` LICENSE. Two
sources the lanes used were already tabled: JASPAR 2026, whose about page reads "licensed under
the Creative Commons Attribution 4.0 International License" and which M1b reads for its motif arm
(`attribution/regulators.py`), and GTEx v8 DAP-G, which M5b's eQTL arm reads
(`attribution/eqtl_arm.py`). AlphaMissense is now listed with two terms because its own two
statements differ; until that is settled the stricter one governs, as the code already applies it.

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
*(2026-10-10.)* **AlphaMissense** belongs on that list while its two statements differ: the hg38
table's own header says CC BY-NC-SA 4.0 and its README says CC BY 4.0 (see its row above).

GenomeOS produces research hypotheses with their evidence and their gaps. It
is not a medical device, it gives no clinical advice, and it never claims a
therapy will work.
