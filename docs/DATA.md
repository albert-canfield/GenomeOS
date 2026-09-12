# Open genome data for testing

All of the following are open access and were verified reachable on
2026-09-10. Fetch with `scripts/fetch_reference.sh <target>`; files land in
`data/reference/` (gitignored).

## A real individual: HG002 (recommended test subject)

HG002 (NA24385) is the son of the Genome in a Bottle Ashkenazi trio and a
Personal Genome Project participant with open consent, so the data are usable
in research and commercially. He is the best-characterised human genome in
existence.

| Target | What | Size | URL |
|---|---|---|---|
| `hg002` | T2T diploid assembly v1.1 (both haplotypes, telomere to telomere) | ~1 GB | https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz |
| `hg002-vcf` | GIAB v4.2.1 benchmark variants against GRCh38 | ~100 MB | https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz |
| `hg001-vcf` | NA12878 / HG001 benchmark variants (the classic CEPH genome) | ~100 MB | https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz |

Raw reads (Illumina, PacBio HiFi, ONT) for the same individuals are under
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ when read-level
work (telomere length estimation, methylation from ONT) begins.

## Reference assemblies

| Target | What | Size | URL |
|---|---|---|---|
| `chrM` | hg38 mitochondrial genome (16,569 bp) | 6 KB | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrM.fa.gz |
| `chr21` | hg38 chromosome 21 (46.7 Mb) | 13 MB | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz |
| `chrN` | any hg38 chromosome | 10–80 MB | same pattern |
| `hg38` | GRCh38 full assembly | ~940 MB | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz |
| `hs1` | T2T-CHM13 v2.0, gapless | ~900 MB | https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz |

## Annotation and knowledge

| Target | What | URL |
|---|---|---|
| `gencode` | GENCODE Release 50 GFF3 (genes, transcripts, exons, CDS) | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gff3.gz |
| — | Gene Ontology basic | https://current.geneontology.org/ontology/go-basic.obo |
| — | Reactome V96 (SBML, BioPAX) | https://reactome.org/download-data/ |
| — | Uberon anatomy, Cell Ontology | http://purl.obolibrary.org/obo/uberon/uberon-ext.obo , http://purl.obolibrary.org/obo/cl.obo |
| — | ClinVar VCF (weekly) | https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz |
| — | dbSNP b157 | https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz |

## Populations

| Target | What | URL |
|---|---|---|
| `1kg-chr21` | 1000 Genomes, 3,202 samples, chr21 genotypes | https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/ |
| — | HPRC Release 2 pangenome (GFA / GBZ) | https://github.com/human-pangenomics/hpp_pangenome_resources |
| — | CELLxGENE Census (single-cell expression, ~33M cells) | `pip install cellxgene-census` |

## Verified today

```
genomeos info data/reference/chrM.fa.gz        # 16,569 bp
genomeos orfs data/reference/chrM.fa.gz --table mito --min-aa 200
    -> MT-CO1 513 aa, MT-CO2 227 aa, MT-ATP6 226 aa (match published lengths)
genomeos info data/reference/chr21.fa.gz       # 46,709,983 bp in 0.5 s
```
Known limitation: the ORF finder only starts at AUG. Human mitochondrial
genes ND1, ND2 and ND5 use AUA/AUU starts, so their ORFs are found at the
first internal AUG. Alternative start codons are on the Phase 1 list.

## The second mammal: mouse mm10 (2026-09-11)

`genomeos mouse --chrom chr19` fetches, once, the mouse chromosome through
the same three doors as a human one: sequence from UCSC (mm10, the assembly
ENCODE's mouse registry is on), gene models from GENCODE vM25 (the
chromosome's rows streamed out of the genome-wide GFF3), and the ENCODE
SCREEN mouse cCREs (rows of the chromosome). Files carry the `mm10_` prefix
under `data/reference` (local) and `data/results/ccres_mm10_<chrom>.bed.gz`
(committed, small); the summary is `data/results/mouse_mm10_<chrom>.json`.
The node comparison maps mouse genes to human through MGI's curated
mouse–human homology report, streamed once (15 MB) and kept as symbol pairs
(`data/results/mgi_mouse_human_orthology.tsv.gz`, 24,584 pairs).

Orthology comes from two curated sources, both streamed once and kept as symbol pairs under
`data/results`: MGI's mouse–human homology classes (`mgi_mouse_human_orthology.tsv.gz`) and Ensembl
Compara 116's gene trees (`compara_mouse_human_orthology.tsv.gz`, 2026-09-12), read from the *mouse*
homology dump (`homologies/mus_musculus/`), since the human dump omits *Mus musculus*. Mouse ids resolve
through the mouse GENCODE files fetched so far, so the Compara pairs cover the mouse chromosomes GenomeOS
has looked at; the file's header counts the unresolved ids. `genomeos mouse` reports both comparisons
and their agreement.

## Optional: AlphaGenome (predicted regulatory effects)

Google DeepMind's AlphaGenome predicts what a sequence change does to
expression, splicing and chromatin in each tissue. It is optional: the core
never depends on it, and its output enters GenomeOS only as `predicted`
evidence capped at confidence 0.7. Without it the features below are loaded
but disabled, and `genomeos predict --status` and the Progress tab say why
and how to enable them.

| | Feature | What the key enables | State |
|---|---|---|---|
| a | variant effect | predicted expression change per tissue for any variant: `genomeos predict chr21:25897620 C>T`, `/api/predict` | built |
| b | regulatory blocks | predicted target gene and tissue for UNKNOWN regulatory elements, tightening the CTCF-domain inference | planned |
| c | sequence grammar | splice-site and chromatin signals the learned matrices cannot capture, entering as predicted evidence | planned |
| d | twin | predicted expression differences between an individual's haplotypes and the reference | planned |

What each feature is worth and how it fits the architecture:
[ALPHAGENOME.md](ALPHAGENOME.md).

To enable:

1. Request a key (free for non-commercial use) at
   https://deepmind.google.com/science/alphagenome/account/terms.
2. Install the client: `uv sync --extra predict`.
3. Put the key in a file named `.env` at the project root (git-ignored,
   never committed) or export it in your shell:

```
ALPHAGENOME_API_KEY=your-key-here
```

4. Check: `uv run genomeos predict --status` prints `enabled`.

The key is read from the environment first and from `.env` second; it is
never written to a result, a log or a job file.

## HG002 on every autosome (2026-09-11)

`scripts/twin_genome_wide.py` applied the GIAB v4.2.1 benchmark calls to
all 22 local autosomes and kept only the statistics
(`hg002_twin_by_chromosome.json`): 4,048,342 PASS variants
(3,460,251 SNVs), 1,619,891 applied on haplotype 1 and
4,032,208 on haplotype 2 (unphased heterozygotes go to
haplotype 2 by policy), 0 reference mismatches and
840 skipped overlaps. Zero mismatches over four million
calls is the check that the fetched hg38 sequences and the benchmark agree
base for base. chrX and chrY are not in the benchmark set. No haplotype
FASTA is kept; `genomeos twin build` makes one on demand.

## Your own genome (`genomeos individual import`, 2026-09-11)

A geneticist's own data enters the way HG002 does. One VCF against GRCh38,
plain or gzipped, from any caller:

    genomeos individual import /path/genome.vcf.gz --name ME --note "GATK 4.5, 2026-08, own consent"
    genomeos individual list
    genomeos individual check --name ME               # apply to GRCh38, count reference mismatches: the assembly check
    genomeos individual screen --name ME              # ClinVar carrier screen: pathogenic alleles carried, offline
    genomeos individual knockouts --name ME           # truncating SNVs genome-wide, homozygous first
    genomeos individual coding --name ME              # coding SNVs by consequence and by gene, missense included
    genomeos individual protein --name ME --gene APP --chrom chr21   # which protein each tissue makes in ME
    genomeos individual report --name ME --out ME.md  # one page from everything computed for the person
    genomeos individual genes --name ME --chrom chr21
    genomeos individual predict --name ME --gene APP --chrom chr21   # AlphaGenome feature d, needs the key
    genomeos lookup chr17:7675088 C>T          # now says whether ME carries it
    genomeos report TP53 --chrom chr17 ...      # lists ME's variants inside the gene
    genomeos twin build --vcf data/individuals/ME/ME_chr21.vcf --sample ME --chrom chr21 --genome data/reference/chr21.fa

The file is streamed once and split into one small PASS file per chromosome
under `data/individuals/<name>/` with a manifest (source file, sample column,
counts per chromosome, phased genotypes, rows skipped as filtered or on
other contigs). Chromosome names `21`, `chr21`, `MT` are normalised to hg38's;
contigs and decoys are dropped and counted. The directory is git-ignored and
no layer sends a byte of it anywhere: every lookup, dossier and gene walk is
local. Genotypes are measured evidence (the caller's); consequences are
derived by the local trace on the canonical transcript. The Twin tab has the
same import, list and gene-by-gene table, and the predicted effect of the
person's regulatory variants on one gene (docs/ALPHAGENOME.md, feature d).
Per-person results never land under `data/results`: the imported files, the
gene walks and the predictions stay under git-ignored directories, so a
named person's genome cannot be committed by accident.

`individual check` is the first thing to run after an import: it applies the
person's variants to the local GRCh38 sequence, haplotype by haplotype, keeps
the statistics only (no FASTA) and counts the calls whose reference allele
disagrees with the reference base. Up to 0.5% mismatches is a match; a file
against GRCh37/hg19 disagrees at most positions and is called out as such,
before anyone reads a consequence off the wrong coordinates. HG002's own
benchmark gives 0 mismatches over 4 million calls, which is what the check
is calibrated against.

`individual screen` is the carrier screen. ClinVar's GRCh38 VCF (190 MB,
4.47 million rows) is streamed once and distilled to its 346,571 pathogenic
and likely pathogenic rows (226,268 and 120,303; conflicting calls left
out), 7.6 MB under `data/knowledge/clinvar`, local; only the counts are
committed (`clinvar_pathogenic_summary.json`). The screen intersects the
person's files with that table, chromosome by chromosome, and reports each
allele carried with genotype, zygosity, ClinVar's review stars and the
conditions. HG002, 4,048,342 variants on 22 autosomes, 7 seconds: two hits,
the F11 c.1716+1G>T-type factor XI deficiency allele heterozygous (two
stars; F11 deficiency is common in the Ashkenazi population HG002 comes
from) and a 9p21 CDKN2B row ClinVar itself labels "likely pathogenic |
protective" with no review stars. Research annotation against a public
database, never a clinical report; the output says so, and the per-person
result stays under the person's directory.

`individual knockouts` walks every coding gene of every chromosome the
person has and lists the SNVs that end the canonical protein early
(nonsense), remove its start or its stop, homozygous first, with the
fraction of the protein lost. HG002: 88 such SNVs in 84 genes (61 nonsense,
11 start lost, 16 stop lost; 30 homozygous), the expected picture of a
healthy genome: olfactory receptors, FUT2 p.Trp154Ter homozygous (the common
non-secretor allele), FCGR2A, SERPINB11. Calls that fall past a stop the
reference itself carries (FCGR2C, a pseudogene in most people) are left out
and counted. SNVs only, so frameshift indels are not in this list; a
truncating variant is a list to look at, not a verdict.

`individual coding` is the whole coding inventory: every coding SNV on a
canonical transcript, by consequence and by gene, protein-changing and
homozygous first. HG002: 21,019 coding SNVs on 22 autosomes, 11,101
synonymous, 9,830 missense, 61 nonsense, 16 stop lost, 11 start lost; 5,438
genes carry a protein-changing variant, 2,529 a homozygous one, and the top
of the list is MUC16 (57 changes in 14,569 residues) and the HLA genes (30
to 37 each), which is what a count of changes per gene measures: length and
polymorphism, not harm. The dossier carries it as a section next to the
truncating list and the ClinVar screen.

The missense variants are ranked by what UniProt says about the residue
they hit, since a count judges nothing: a variant on an annotated site
(active or binding site, modified residue, glycosylation, the two cysteines
of a disulfide bridge, a motif) first, then one inside a domain, then the
rest, homozygous before heterozygous. HG002: 46 of the 9,830 missense
variants sit on an annotated site (5,374 more inside a domain), among them
CD52 p.Asn40Ser homozygous, which removes an N-glycosylation site, and two
cysteines of disulfide bridges lost in GALNTL5 and an olfactory receptor.
A rank from annotation, not a prediction of effect; the residue-level view
of any of them is one click away on the Flow tab.

`individual import` also takes an http(s) URL and streams it, which is how
the GIAB parents arrive: HG003 and HG004 (v4.2.1 benchmark calls, open
consent) as local individuals. `individual regions --name X BED` keeps the
regions a person's calls are trusted in (GIAB's benchmark regions, or a
caller's callable track), and `individual trio --name HG002 --father HG003
--mother HG004` reads inheritance chromosome by chromosome: 4,096,122 child
variants on 22 autosomes, 96.4% inherited (2,336,470 from both parents).
Without trusted regions, 100,114 child calls were absent from both parents
and 48,848 were "Mendelian errors" (homozygous in the child, absent from a
parent); with all three region sets applied, 12,296 and 2,930, and 133,736
child calls sit outside a parent's trusted region and are counted as
untrusted rather than as events. Real de novo variation is about 60 to 100
per genome, so the 12,296 were mostly representation differences between
three separately produced call sets: the same insertion written TGG>TGGG in
one file, T>TG in another, a deletion of two repeat units placed at a
different unit. Since 2026-09-12 the comparison normalises every allele
first (common suffix and prefix trimmed to one anchor base, indels
left-aligned along the local reference), and the excess falls to 1,430
candidates (1,173 SNVs, 257 indels) and 7 Mendelian errors, with 124,881
child calls outside a parent's region. What remains is child PASS calls
inside both parents' high-confidence regions with no parent row within a
base: the true de novos plus the benchmark files' own disagreements. The
imported files carry genotypes only, unphased, so phasing by parent waits
for a phased import. Stored under the child's directory as
`trio_<father>_<mother>.json`, never under `data/results`.

The inventory also surfaces what the reference hides. hg38 itself carries a
frameshift or nonsense allele against the curated protein in 31 genes (the
verified translation disagreements triaged by mechanism); a person who
matches the reference across such a gene carries that truncation, and no
caller will ever list it as a variant. HG002 matches hg38 in 5 of the 30 on
file (OR2T7, OR4C45, OR4K3, OR1P1, SCYGR10, olfactory receptors and a
keratin-associated gene) and differs somewhere inside the other 25, though
none of those differences restores the curated protein. The dossier says so.

`individual protein` is the twin's isoform: for one gene, GTEx's dominant
transcript in each tissue, traced, with the person's coding variants read on
that transcript rather than on the canonical one. HG002 × APP: APP-201 (770
aa) in 33 tissues, APP-204 (751 aa) in 14, APP-202 (695 aa, the neuronal
isoform) in 6 brain regions, APP-203 in the cortex, no coding variant on
any of them; HG002 × FUT2: FUT2-201 in 49 tissues carrying p.Trp154Ter
homozygous, so every tissue that makes FUT2 makes the truncated one. The
isoform per tissue is GTEx's population median, not this person's own
expression, and the output says so.

Checked on HG002's chr21 rows re-imported under another name: 55,210 PASS
variants; 194 of 221 coding genes carry a variant, 269 coding SNVs (135
missense, 133 synonymous, 1 nonsense), the KRTAP10 cluster on top as
expected for a highly polymorphic family.

## Telomere from a BAM read by ranges (2026-09-12)

GIAB's 300x Illumina BAM of HG002 on GRCh38 (601 GB, NCBI FTP over HTTPS,
`NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam`) answers
range requests and comes with its index (`.bai`, 12 MB). `genomeos telomere
<bam-url> --bai <index> --save NAME` never downloads the BAM: the index's
pseudo-bins give the mapped read count per chromosome and the unmapped count
(`genome/bam_range.py`, standard library only); the unmapped tail, where the
pure TTAGGG reads sit, is sampled by range (100 MB by default) and scaled by
that count; the last and first 10 kb of assembled sequence of each chromosome
(found from the local reference, past the terminal N runs) are read whole for
the boundary reads and the coverage; TelSeq's arithmetic then gives a length
per chromosome end. The index is cached under `data/knowledge/telomere/`
(git-ignored); the result for a GIAB person goes under `data/results`.
Evidence: the reads are `experimental`, the estimate `inferred` (0.3), since no
GC normalisation is applied and the sampled tail is a sample. CRAM is not
readable without htslib, so the range route is BAM only.

## Optional: AlphaMissense (predicted missense effect, 2026-09-12)

`genomeos individual coding --name N --predict` streams DeepMind's AlphaMissense
hg38 table (643 MB compressed, 71 million rows, public bucket
`dm_alphamissense`) once and keeps only the rows for that person's missense
variants, under the person's own git-ignored directory
(`data/individuals/N/alphamissense.json`, with the variants asked so a second
call streams nothing). Each variant gets the score on the person's transcript
when the table has it, else the highest across transcripts, with the model's
class (likely pathogenic ≥ 0.564, likely benign ≤ 0.34, ambiguous between) and
confidence capped at 0.7; indels and variants the table lacks stay unscored.
Evidence class `predicted`. Licence: CC BY-NC-SA 4.0, research and personal use
only, so the scores never enter a committed result.

## Optional: 4D Nucleome boundary calls (measured node edges)

The 4DN portal publishes boundary calls (insulation-based BED files) for
hundreds of released GRCh38 Hi-C and Micro-C experiments: GM12878 (11
files), H1-hESC, HFFc6, K562, HCT116, HepG2, IMR-90, HeLa-S3 among them.
The search is open; the downloads return 403 without an account key, and
a free account is all it takes: create one, generate an access key at
https://data.4dnucleome.org and put

```
FOURDN_KEY=...
FOURDN_SECRET=...
```

in the git-ignored `.env`. `genomeos domains --chrom chr21 --hic GM12878`
then streams the biosource's deepest boundary file once into
`data/knowledge/hic/<biosource>/<chrom>.bed` (local) and holds the CTCF-only
boundaries against the measured ones with the same comparison and random
control the predicted contact map got (`genomeos/genome/hic.py`,
`data/results/domains_<chrom>_hic_<biosource>.json`). Without the key the
feature loads disabled and says so; checked 2026-09-12, blocked on the key
only. ENCODE has no open Hi-C domain BEDs and the 3D Genome Browser's hg38
TAD archive is no longer served, so 4DN is the source.

## Measured enhancers: VISTA (2026-09-12)

The VISTA Enhancer Browser (LBNL) publishes its data on GitLab
(`egsb-mfgl/vista-data`): `locus.tsv.gz` (120 kB) lists every tested element
with hg38 coordinates, curation status (positive or negative at e11.5 in the
transgenic mouse) and the tissues where it drove expression. `genomeos`
streams that file once into `data/knowledge/vista/` (git-ignored) the first
time `scripts/vista_score.py` or the `vista_genome_wide` job runs; the
per-element rows it derives stay there too, and only the per-chromosome
summaries with the rows the docs quote are committed (`vista_chr*.json`,
`vista_genome_wide.json`). About 2,440 human elements, 1,267 positive and
1,175 negative, none on chrY. Evidence class `experimental`.

## Measured targets: GTEx v8 cis-eQTLs (2026-09-12)

GTEx's single-tissue eQTL archive (`GTEx_Analysis_v8_eQTL.tar`, 1.56 GB, 49
tissues) sits in a public Google Cloud bucket that answers range requests.
`attribution/eqtl.py` reads the tar's member headers over ranges, streams
each tissue's significant variant-gene pairs (71.5 million rows in all)
through gzip without touching the disk, and keeps only the pairs whose
variant falls inside (± 500 bp) an element a committed result carries a
prediction for: one small TSV per tissue under `data/knowledge/gtex/`
(git-ignored, 34 MB), resumable tissue by tissue. Nothing else of the
archive is stored. `scripts/eqtl_targets.py` (job `eqtl_targets`) distils
when the hits are missing and scores; `eqtl_targets.json` is committed.
Variant ids are hg38 (`chr1_13550_G_A_b38`). Evidence class `experimental`.

## Measured activity: ENCODE4 lentiMPRA (2026-09-12)

The Ahituv lab's joint lentiMPRA library (ENCODE reference ENCSR106SZM:
53,990 200-bp elements assayed in K562 ENCSR203UFY, HepG2 ENCSR405QCT and
WTC11 ENCSR336MKI) is read from the three element-quantification BED files
(ENCFF802FUV, ENCFF475FKV, ENCFF769REH; 1.1 MB each), streamed once into
`data/knowledge/mpra/` (git-ignored) by `attribution/mpra.py`; forward and
reverse copies of an element are averaged. Column 7 is log2(RNA/DNA). The
per-element rows stay local (`data/knowledge/mpra/rows_<chrom>.json`); the
per-chromosome summaries (`mpra_chr*.json`) and the fold
(`mpra_genome_wide.json`) are committed. Predicted DNase per cell line comes
from AlphaGenome (K562 EFO:0002067, HepG2 EFO:0001187) through
`predict/chromatin_tracks.py`, cached as element means under
`data/knowledge/alphagenome/dnase/`. Evidence class `experimental` for the
assay and the reader's peaks, `predicted` for the model's track.

## Consequence: GWAS Catalog and ClinVar non-coding (2026-09-12)

The NHGRI-EBI GWAS Catalog's full association table
(`gwas-catalog-associations-full.zip`, 69 MB, on EBI's FTP under
`releases/latest`) is streamed once by `attribution/gwas.py`; the 1.2 million
rows are read from the zip in memory and only the lead variants within 1 kb of
an element with a prediction, plus the same elements shifted 100 kb as the
chance control, are kept in `data/knowledge/gwas/hits.tsv` (git-ignored, a few
MB). ClinVar's distilled pathogenic table (see the individual's screen) now
carries the molecular consequence as its last column, so
`scripts/consequence_targets.py` can keep the non-coding ones; re-run
`genomeos individual screen --distil` (or `clinvar.distil()`) once to add the
column to an older local copy. `consequence_targets.json` is committed.

## Optional: peptide/HLA binding predictors

The therapeutic pipeline enumerates the peptides a mutation creates; whether
one of them binds the patient's HLA needs a predictor, and GenomeOS ships
none. Two can be enabled, and the analysis runs without either (binding is
reported as unavailable, never invented).

| Predictor | Licence | How to enable |
|---|---|---|
| MHCflurry 2 | Apache 2.0, usable by anyone | `uv sync --extra hla` then `uv run mhcflurry-downloads fetch models_class1_presentation` |
| NetMHCpan 4.1 | academic licence from DTU Health Tech | request the package at https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/, install it, then put `netMHCpan` on `PATH` or set `NETMHCPAN=/path/to/netMHCpan` |

`genomeos therapeutic --hla-predictor mhcflurry|netmhcpan|none` chooses one;
the default is the first installed. `/api/features` and the Progress tab list
both with their licence and their state. GenomeOS never downloads or
redistributes NetMHCpan; a commercial user needs a licence from DTU.
