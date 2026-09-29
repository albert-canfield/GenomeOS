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

`individual diff --name A --against B` (or `--against reference`) is the same
comparison read as a code review: the two people's protein-changing variants
gene by gene, alleles normalised, the lines only A carries marked `+`, only B
`-`, with each variant's consequence on the canonical transcript, its
zygosity, AlphaMissense's class where the person has been scored and ClinVar's
row where the person has been screened; genes are ordered by what matters
(ClinVar first, then truncations, then predicted class). The inventory keeps
every protein-changing variant under the person (`coding_variants.json`) for
it; a person not yet inventoried is traced once. Markdown to the terminal or
`--out`, never under `data/results`.

`--regulatory --chrom C` reads the same two people at the regulatory level:
the variants inside the chromosome's ENCODE elements, each with the element's
class, the target the CTCF node infers for it (with basis and distance), the
AlphaGenome deletion target where the element has been scored, and, with
`--constraint`, the base's own phyloP read by range; grouped by target gene,
constrained and predicted first. A variant in an element is a candidate, not
an effect, and the page says so.

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

## Phasing the trio candidates on the Q100 assembly (2026-09-27)

The trio above left 1,430 de novo candidates and the note that phasing them
by parent "waits for a phased import". The phased resource already exists.
It is the one the resolution below uses.

**The resource, from primary sources (read 2026-09-27).** NIST/GIAB's v5.0q
release,
<https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/>,
is an assembly-based benchmark built on v1.1 of the T2T HG002 Q100 diploid
assembly (<https://github.com/marbl/hg002>,
<https://doi.org/10.1101/2025.09.21.677443>). Its `dipcall_output/` directory
holds the assembly's own variant calls against GRCh38,
`GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.dip.vcf.gz` (40 MB, index 1.5 MB), and
the regions both haplotypes cover, `...dip.bed` (9.5 KB, 2.84 Gb). Beside
them is the curated benchmark: `HG002_GRCh38_v5.0q_smvar.vcf.gz` (40 MB) and
its `.benchmark.bed` (2.74 Gb).

- **Phased and labelled by parent.** Every genotype is written `a|b`. The
  release's own configuration (`defrabb_files/resources.yml`) names the two
  inputs `hg002v1.1.mat.fasta.gz` and `hg002v1.1.pat.fasta.gz`. The pipeline's
  dipcall rule (usnistgov/defrabb, `rules/asm-varcall.smk`) passes
  `paternal.fa` first and `maternal.fa` second, and dipcall's README says
  that for a male sample parent 1 is taken as the father. So haplotype 1 is
  paternal and haplotype 2 is maternal. I checked this on the file itself:
  every PASS call in the chrX non-PAR is `.|1` (136,036 calls), and every
  chrY call is `1|.` (31,822 calls). A son's X comes from his mother and
  his Y from his father.
- **Licence.** The README's NIST Data Use Policy says the data were created by
  NIST employees and are not subject to copyright in the US (17 USC 105).
  They are provided "AS IS". Outside the US, NIST grants a royalty-free right
  to copy, modify, prepare derivative works and distribute them. The one
  condition is acknowledgement of NIST as the source, with a note of any
  change. The files may be downloaded and distilled here. They are stored
  under the git-ignored `data/cache/q100/` and md5-checked against the
  release's `checksum.md5`. Only counts are committed.
- **Sizes.** 42.4 MB for the four files used, far under the 5 GB limit. The
  per-haplotype BAMs (1.0 GB each) are not needed.
- **What the parent label can and cannot do.** Trio binning assigns each of
  the child's haplotypes to the parent it came from. An allele on the
  paternal haplotype sat on the chromosome the father transmitted. That is
  true whether the father carries the allele or it arose in his germline. So
  the assembly says *which parent's chromosome*, not *inherited or de novo*.
  For one candidate, that distinction still rests on the parents' calls,
  which lack the allele by definition. For the set, the distinction shows in
  the paternal share. Germline de novos are about 75 to 80% paternal (Kong et
  al. 2012; Jónsson et al. 2017). Inherited alleles that a parent's call
  missed would split about evenly.

**Pre-registration, written and committed before the phasing ran.** The
code is `phase_trio_by_assembly` in `genomeos/genome/individuals.py`, run by
`scripts/trio_q100_phase.py`. It rebuilds the candidates exactly as `trio`
defines them, and the run is refused unless the count is 1,430. Then it
reads each candidate against the assembly, with alleles normalised the same
way. It also counts overlap with the v5.0q benchmark regions and with NIST's
HG002 de novo and mosaic exclusion regions. That BED is named in the same
`resources.yml` (1,913 intervals, 7.2 Mb including the repeats it was
widened over).

Each candidate falls into one class:

| Class | The assembly | Read as |
|---|---|---|
| outside | the site is outside the diploid regions | not assessable |
| father / mother | the same allele, on one haplotype only | the child's allele is confirmed and sits on that parent's chromosome: a de novo from that parent's germline, or an inherited allele the parent's call missed |
| both | the same allele, on both haplotypes | not a de novo: two independent events at one base are not credible, so at least one parent's call missed it (or both sources share an error) |
| filtered | the same allele, but dipcall's own filter failed | the assembly is uncertain there |
| nearby | no such allele, but a variant within 10 bases | a representation disagreement the normalisation did not reconcile |
| absent | reference within 10 bases | the assembly does not support the child's call: a representation artefact or a false positive in the v4.2.1 child call, or an assembly error |

The control uses the same trio, the same regions and the same code. It takes
the child's heterozygous calls that exactly one parent carries, so the parent
of origin is known from the trio alone.

- **C1.** At least 97% of the control calls inside the diploid regions are
  found as the same allele on exactly one haplotype.
- **C2.** At least 99% of those are on the haplotype of the parent who
  carries the allele.

**Predictions.** The hypothesis is the one the trio note stated: most of the
1,430 are representation artefacts, and the true de novos are the 60 to 100
per genome the literature gives.

- **P1.** father + mother is between 40 and 150. The range is wider than 60
  to 100 because the three trusted-region sets cover part of the genome and
  a v4.2.1 child call can miss a de novo.
- **P2.** Among those, the paternal share is between 0.70 and 0.90.
- **P3.** nearby + absent is at least half of the assessable candidates.
- **P4.** both is under 5%.
- **P5.** Of the father + mother candidates, at least 30% fall in NIST's
  de novo and mosaic exclusion regions, if they are true de novos. This
  would be an independent confirmation.

**The rival reading.** If father + mother is well above 150 and the paternal
share is near 0.5, the excess is inherited variation that the parents' v4.2.1
calls missed, not artefacts on the child's side. P1 to P3 would then fail,
and the note above would be wrong about where the excess comes from.

**What would show the method wrong.** Any of these voids the reading of the
candidates, and it is recorded as void rather than reinterpreted:

- the candidate total is not 1,430;
- C2 is below 95%, which means the haplotype labels are misread;
- C1 is below 90%, which means absence from the assembly does not
  discriminate at this normalisation.

**The assembly is a benchmark too.** The v5.0q README says excluded regions
remain for assembly errors, mosaic variants and places where the assembly
could be aligned another way. A candidate that v4.2.1 calls and the assembly
does not is a disagreement between two measurements. It is not proof that
the child call is wrong. Counts are also reported for the subset inside the
v5.0q benchmark regions, where NIST vouches for the assembly's calls. One
effect is expected there. If the father + mother candidates are true de
novos, they fall outside the benchmark regions more often than the control
does, because NIST excluded known de novo and mosaic sites from them.

**Result.** Filled in below after the run. It goes to
`data/results/trio_q100_phase.json`, which holds counts per class and never
a position, an allele or a genotype.

**Result, run 2026-09-27 after the registration was committed (457702c).**
`scripts/trio_q100_phase.py` took 28 s. It reproduced 1,430 candidates, so
the run stands.

- **Control.** 1,484,060 one-parent heterozygous calls; 1,313 fall outside
  the diploid regions.
  - C1 passes: 1,481,231 of the remaining 1,482,747 (0.99898) are found on
    exactly one haplotype.
  - C2 passes: 1,481,214 of those are on the carrying parent's haplotype,
    with 17 on the other (0.99999).
  - Of the rest, 1,187 are nearby, 287 absent and 42 on both haplotypes.

  The labels are read correctly, and the assembly carries ordinary inherited
  alleles almost without exception.
- **Candidates.** 19 of the 1,430 fall outside the diploid regions. Of the
  1,411 assessable, the classes are:

  | Class | Count |
  |---|---|
  | father | 703 |
  | mother | 679 |
  | absent | 20 |
  | nearby | 9 |
  | both | 0 |
  | filtered | 0 |

  Of the 1,155 assessable SNVs, 584 are father, 564 mother, 6 absent and 1
  nearby. Of the 256 assessable indels, 119 are father, 115 mother, 14 absent
  and 8 nearby. Every candidate is heterozygous in the child.
- **Verdicts.**
  - P1 fails: 1,382 are supported on one haplotype, against 40 to 150.
  - P2 fails: the paternal share is 0.509 (703 / 1,382; 95% interval about
    0.48 to 0.54), against 0.70 to 0.90.
  - P3 fails: nearby + absent is 29 of 1,411 (2.1%), against at least half.
  - P4 holds: both is 0.
  - P5 holds as a number: 953 of the 1,382 (69.0%) lie in NIST's de novo and
    mosaic exclusion regions, which cover about 0.25% of the genome. It is
    not the independent confirmation the registration hoped for, because how
    NIST built that list is not documented in the release.
  - The registered rival reading's signature appeared: many more than 150,
    at a share near 0.5.
- **Benchmark regions.** 420 of the 1,411 (29.8%) lie inside the v5.0q
  benchmark regions: 194 father, 210 mother, 9 nearby, 7 absent. The
  benchmark regions cover 96.4% of the diploid regions' length. The
  candidates avoid them, as expected if NIST excluded these sites as de
  novo or mosaic. The control was not counted against this BED, so this
  comparison is by length, not against the control.

**Exploratory, not registered: the population check.** The rival reading
named a cause, inherited variation that the parents' calls missed.
`scripts/trio_q100_population.py` tests that cause. It asks how often an
allele is carried by any of the 88 HPRC year-1 haplotypes in the local human
panel (UCSC `hprc/cactus90way`). None of those assemblies is HG002: its HPRC
accessions are absent from the panel's list.

| Group | Seen in the panel |
|---|---|
| control (every 100th one-parent inherited SNV) | 11,950 / 12,871 (92.8%) |
| supported candidate SNVs, all | 35 / 1,148 (3.0%) |
| inside NIST's de novo and mosaic regions | 11 / 900 (1.2%) |
| outside those regions | 24 / 248 (9.7%) |

A missed inherited polymorphism would be seen about as often as the control.
The candidates are new alleles. So the rival's cause is contradicted, even
though its signature appeared.

**Reading.** The 1,430 are neither mostly representation artefacts nor mostly
parent-side misses. The assembly independently carries 1,382 of them, and
the child's two measurements agree:

- 97.9% sit on exactly one parental chromosome;
- they split evenly between the father's and the mother's chromosomes;
- they are absent from the population;
- 69% sit where NIST itself set sites aside as de novo or mosaic.

A germline de novo set would be about 80% paternal. An even split is what a
mutation arising after fertilisation gives: in the embryo, or in the
lymphoblastoid cell line from which GIAB's HG002 DNA comes, the assembly's
DNA included. Both measurements would then carry it. This is the
interpretation that fits all four observations, not a measurement of it.
Separating early-embryonic variants from culture-acquired ones needs read
allele fractions or DNA that was never in culture.

The share cannot count the germline de novos mixed in. A 0.80-paternal
component of X in a 0.50 background gives X ≈ (703 − 691) / 0.3 ≈ 40, with
a standard error near 60. That is consistent with anything from none to the
literature's 60 to 100. The trio note above, "a larger count is
representation differences", was wrong for the count after normalisation:
the normalisation had already removed nearly all of those.

**Status.** Area D's phasing item is not blocked: the phased, parent-labelled
resource exists and is used here. Each candidate now has a parental
chromosome, and the set has a reading. What remains open is a germline versus
post-zygotic call for each candidate, which needs reads. Result files hold
counts only: `data/results/trio_q100_phase.json` and
`data/results/trio_q100_population.json`.

## Allele fractions of the trio candidates (2026-09-28)

The phasing above left one question per candidate: germline, or arisen after
fertilisation? Reads answer part of it. A germline heterozygous allele is in
every cell, so about half the reads carry it. An allele present in only some
of the sequenced cells is carried by fewer.

**The resource, from primary sources (read 2026-09-28).** GIAB's v4.2.1
HG002 GRCh38 benchmark VCF,
<https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz>
(156,252,944 bytes, md5 dc750b3807d4af1f7ffec852e9c2f771; the release's
`md5.in` does not list it), is the file the trio's HG002 calls were cut
from. Its FORMAT carries `ADALL`, "net allele depths across all datasets",
and `AD`, the same over the unfiltered datasets with a called genotype. So
every candidate already has read counts, pooled over Illumina, PacBio HiFi,
10x and the other datasets. It is NIST data, not subject to US copyright
(17 USC 105), provided AS IS, used here with acknowledgement. It is stored
under the git-ignored `data/cache/v421/`. No BAM or CRAM is read. The 300x
BAM that `genome/bam_range.py` streams for the telomere would cost several
GB at 1,430 positions, above this lane's 2 GB cap, and is not needed.

**A truncation to state before the run.** The v4.2.1 README (the v3.3
changes, kept since) says heterozygous calls with a net allele fraction
below 0.2 or above 0.8 are excluded from the benchmark regions. The 1,430
candidates were drawn inside those regions. A mosaic allele below 0.2 could
never have become a candidate. This measurement sees only the 0.2 to 0.8
window.

**What allele fraction can and cannot say here.** HG002's DNA comes from a
lymphoblastoid cell line. Such lines grow from one or a few B cells. Any
mutation that the founding cell carried, whether it arose in a parent's
germline, in the early embryo, or in that B cell during the donor's life, is
in every cell of a clonal line and sits at 0.5. Somatic mutations
accumulate in B lymphocytes with age, from hundreds per cell at birth to
over a thousand in adults (Zhang et al. 2019, PNAS 116:9014). So:

- a fraction clearly below 0.5 means the allele is in part of the sequenced
  cells: a mosaic from the embryo, or a subclone of the culture;
- a fraction near 0.5 means clonal in the sequenced cells. It does not mean
  germline.

The phasing already gives an expectation. The Q100 assembly carries 1,382
candidates on one haplotype. A consensus assembly carries an allele on a
haplotype only when most of that haplotype's reads show it. So most of
those candidates should be clonal or close to it.

**Pre-registration, written and committed before the run.** The code is
`trio_allele_fractions` in `genomeos/genome/individuals.py`, run by
`scripts/trio_vaf.py`. It rebuilds the candidates exactly as `trio` defines
them, and the run is refused unless the count is 1,430. It reads `ADALL`
with alleles normalised the same way. The control is the phasing's: the
child's heterozygous calls that exactly one parent carries, from the same
file.

Each call is classed from its reference and alternate read counts:

| Class | Rule | Read as |
|---|---|---|
| shallow | fewer than 30 reads for the two alleles | not assessed |
| low | fraction below 0.40, and a one-sided binomial test against 0.5 gives p < 0.001 | in only part of the sequenced cells |
| high | fraction above 0.60, and the mirror test gives p < 0.001 | more than half the reads (a copy-number or mapping effect) |
| half | anything else | clonal in the sequenced cells |

Pooled depth is high, so the binomial test alone would flag ordinary hets
that sit at 0.47 through reference bias. The 0.40 bound stops that. The
control measures how often an ordinary het is classed low anyway.

**Checks.** The reading is void, and recorded as void, if any fails:

- **C1.** The control SNVs' median fraction is between 0.45 and 0.55.
- **C2.** At most 5% of the control SNVs are classed low.
- **C3.** At least 80% of the control SNVs are assessed (depth 30 or more).
- The candidate total is not 1,430.

**Predictions**, for the father + mother candidates (on one haplotype of
the assembly):

- **P1.** At most 25% are classed low. This follows from the assembly and
  the 0.2 truncation.
- **P2.** More are classed low than the control's rate predicts
  (the control's SNV and indel rates weighted by the candidates' own mix;
  one-sided binomial p < 0.01). The lane-q100 reading, that most candidates
  arose after fertilisation, allows some subclonal ones even inside the
  window.

**Outcomes and readings.** Let s be the share of assessed father + mother
candidates classed low, and s0 the control's rate for the same mix.

| Outcome | Condition | Reading |
|---|---|---|
| clonal | s ≤ s0 + 0.05 | the candidates are clonal in the sequenced DNA. Allele fraction cannot tell germline from the founding cell's own mutations. With the even parental split (germline would be about 80% paternal), the reading becomes: mostly mutations of the cell that founded the line, or of the very early embryo, with germline de novos a minority the split cannot count (about 40 ± 60). The share arisen after fertilisation is not measurable from these reads. |
| clonal majority with a subclonal minority | s0 + 0.05 < s < 0.50 | s is a lower bound on the subclonal share: fractions below 0.2 never entered the set, and clonal candidates can also be post-zygotic |
| subclonal majority | s ≥ 0.50 | most candidates are in only part of the cells: embryonic mosaics or culture subclones, measured rather than inferred. P1 fails. |

**Falsifier of the expectation.** s above 0.25 fails P1: the assembly's
consensus would then carry many alleles that most of the reads lack, and the
argument from the assembly would be wrong.

Also reported, not predicted: the same classes for all 1,430, for SNVs and
indels, for the 48 candidates outside the assembly or absent or nearby, and
for the father + mother candidates inside and outside NIST's de novo and
mosaic exclusion regions.

**Result.** Filled in below after the run. It goes to
`data/results/trio_vaf.json`: counts, medians and 0.05-wide histograms per
group, never a position, an allele or a genotype.

**Result, run 2026-09-28 after the registration was committed (bf789a3).**
`scripts/trio_vaf.py` took 36 s. It reproduced 1,430 candidates, so the run
stands. Negatives first.

- **Neither pure outcome.** The candidates are not clonal like the control,
  and they are not mostly subclonal. The registered outcome is the middle
  one, "clonal majority with a subclonal minority". The next point shows
  that label misdescribes the shape.
- **The registered classes miss the main feature.** The rule expected two
  groups, one at 0.5 and one clearly lower. The candidates form one group,
  shifted: their median fraction is 0.442 against the control's 0.495. The
  0.40 bound puts most of that group in "half". So the 11% classed low are
  mostly the lower tail of the shifted group, not a separate subclonal
  minority.
- **Checks pass.** C1: the control SNVs' median is 0.495. C2: 0.95% of
  control SNVs are classed low. C3: 99.5% are assessed (median depth 322).
  The reading is not void.

| Group | Assessed | Low | Share low | Median fraction | Median depth |
|---|---|---|---|---|---|
| control, one-parent inherited hets | 1,476,629 | 14,653 | 0.99% | 0.495 | 321 |
| control SNVs | 1,282,858 | 12,179 | 0.95% | 0.495 | 322 |
| control indels | 193,771 | 2,474 | 1.28% | 0.499 | 318 |
| candidates, all | 1,414 | 167 | 11.8% | 0.442 | 328 |
| father + mother | 1,375 | 156 | 11.3% (95% 9.8 to 13.1) | 0.442 | 329 |
| father | 699 | 82 | 11.7% | 0.447 | 330 |
| mother | 676 | 74 | 10.9% | 0.439 | 329 |
| absent from the assembly | 15 | 8 | 53% | 0.389 | 328 |
| nearby | 8 | 0 | 0% | 0.478 | 331 |
| outside the assembly | 16 | 3 | 19% | 0.451 | 77 |
| father + mother, in NIST's de novo and mosaic regions | 953 | 112 | 11.8% | 0.440 | 323 |
| father + mother, outside them | 422 | 44 | 10.4% | 0.447 | 350 |

No candidate is classed high. 16 of the 1,430 lack `ADALL` or have under 30
reads. Share of each group by fraction (0.05-wide bins):

| Fraction | 0.35 to 0.40 | 0.40 to 0.45 | 0.45 to 0.50 | 0.50 to 0.55 | 0.55 to 0.60 |
|---|---|---|---|---|---|
| control, all | 0.6% | 7.7% | 46.2% | 38.7% | 3.9% |
| father + mother | 10.2% | 46.5% | 35.1% | 5.8% | 0.7% |

- **Verdicts.**
  - P1 holds: 11.3% are classed low, against at most 25%.
  - P2 holds: the control's rate for the same SNV and indel mix predicts
    1.0%; 156 of 1,375 gives one-sided p ≈ 7 × 10⁻¹⁰⁸.
  - Reported, not predicted: the absent class, the one the assembly does not
    carry, has the highest low share (8 of 15), as the assembly argument
    implies.

**Exploratory, not registered.**

- **A clone at about 88% of the cells.** A heterozygous mutation in a
  fraction c of the cells sits at a read fraction c/2. A median of 0.442
  gives c ≈ 0.88. The depth is the control's, and the shift is the same for
  SNVs (0.441) and indels (0.456), for the father's and the mother's
  chromosomes, and inside and outside NIST's regions. The simplest reading is
  that most candidates are mutations of one clone that makes up about 88% of
  the sequenced cells, with the rest of the cells lacking them. The mutations
  predate the clone's expansion: in the B cell that founded the line, or in
  its lineage in the donor. They are not germline: a germline allele sits at
  0.495 like the control.
- **The germline de novos are in the upper tail.** 89 of the 1,375 sit at
  0.50 or above, where 42.7% of the control sits. Even if all 89 were
  germline, the germline count is at most 89 / 0.427 ≈ 208. Among those 89,
  58 are on the father's chromosome (0.65, two-sided p ≈ 0.006 against 0.5).
  Below 0.50 the split is 641 to 645 (0.498). The paternal excess expected of
  germline de novos appears where germline alleles would sit, and nowhere
  else.

**Reading.** Read fractions replace the phasing's interpretation with a
measurement, but a different one from what was asked. Most of the 1,382
are not germline. They are at about 0.44, not 0.5. They are also not a
scatter of late subclones: they share one fraction, as a single clone's
mutations would. The measured subclonal share, 11.3%, is a lower bound on
post-zygotic events and mostly reflects that shared fraction. A per-candidate
germline call is still not possible: a germline de novo and a clonal somatic
mutation overlap within the binomial spread at depth 330. What the set now
supports is a count: germline de novos are at most about 200, a paternal
excess sits among the candidates at 0.5, and the bulk belong to one clone at
about 88% of the cells. Separating early-embryonic from line-founder
mutations needs DNA that never went through culture, which this lane did not find
among GIAB's public HG002 data.

The result file is `data/results/trio_vaf.json`. It holds counts, medians
and histograms only.

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
Evidence: the reads are `experimental`, the estimate `inferred` (0.3), since
the sampled tail is a sample. Two figures come out. `telomere_bp` is the
telomeric read fraction times the genome over its 92 ends: a derivation, and
it reads low against Southern blots. TelSeq instead divides by the reads whose
GC matches a telomeric read's (48 to 52%, since TTAGGG is half G or C) and
multiplies by a constant fitted to blot lengths; that constant cannot be
derived, so the GC correction is reported as
`telomeric_per_gc_comparable_read`, a dimensionless index, with the GC shares
measured separately on the unmapped tail and on 200,000 mapped reads from the
middle of chromosomes (HG002: 4.3% and 4.1% of reads in the band, index
0.0019 against a plain fraction of 0.00008). Both figures compare people read
the same way; neither is a blot length. CRAM is not
readable without htslib, so the range route is BAM only.

## The deletion cache, packed (2026-09-13)

Every element AlphaGenome has been asked about is cached as one small JSON
under `data/knowledge/alphagenome/elements/<chrom>/<element>.json` (about 14
kB each; a chromosome scored whole is 12,000 to 30,000 of them, and the genome
would be near a million files and 13 GB). A chromosome whose all-elements job
is complete never needs those files individually again, so
`scripts/pack_element_cache.py` folds each into one gzipped archive
(`<chrom>.json.gz`, about eight times smaller) and removes the files only
after every answer is inside; `load_cached` reads the per-element file when it
exists and the archive otherwise, so nothing else changes. chr21, chr22 and
chrY packed from 351 MB in 32,804 files to 34.7 MB in three.

## Polygenic scores: the PGS Catalog (2026-09-12)

`genomeos individual pgs --name N [--score PGS000018 ...]` streams a score's
harmonised GRCh38 weight table from the PGS Catalog's FTP (EBI; a few hundred
variants to a few million per score) once into `data/knowledge/pgs/`
(git-ignored), with the score's metadata from the catalog's REST API (name,
trait, publication, the ancestries it was built and tested in, licence). The
sum runs over the person's genotypes chromosome by chromosome: a called
variant gives its effect-allele dosage; a position with no call inside the
person's trusted regions is homozygous reference (dosage 2 when the effect
allele is the reference base, else 0); a position outside the trusted regions
is missing and counted, never guessed; an allele that matches neither the
person's reference nor alternative is a mismatch and counted. The result is a
raw weighted sum with its coverage, not a percentile: the catalog publishes
weights, not any population's distribution, so the number places the person
only against others scored the same way (`compare` gives a z within the local
people), and the ancestry caveat is carried with every score. Defaults: breast
cancer (PGS000004, 313 variants), LDL cholesterol (PGS000115), height
(PGS000297), coronary artery disease (PGS000018, 1.7 million). Evidence:
weights `curated`, the sum `derived`, confidence 0.3. Results under the
person's directory (`pgs_<id>.json`), never under `data/results`.

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
feature loads disabled and says so. ENCODE has no open Hi-C domain BEDs and
the 3D Genome Browser's hg38 TAD archive is no longer served, so 4DN is the
source. The key went in on 2026-09-12 (the portal answers an authenticated
download with a redirect to S3, and S3 refuses the portal's Authorization
header, so the reader drops it on the way); chr21's 227 CTCF-only boundaries
against five biosources, within 20 kb:

| biosource (assay) | measured boundaries | inferred supported | at random | measured covered | median distance |
|---|---|---|---|---|---|
| GM12878 (in situ Hi-C, 4DNFINFT8NR9) | 557 | 58% | 42% | 27% | 15 kb |
| H1-hESC (Micro-C, 4DNFIM86H1MI) | 149 | 24% | 13% | 36% | 49 kb |
| K562 (in situ Hi-C, 4DNFI4EFYN3Q) | 93 | 14% | 8% | 34% | 99 kb |
| HepG2 (in situ Hi-C, 4DNFIFH6XF5T) | 93 | 15% | 8% | 37% | 99 kb |
| IMR-90 (dilution Hi-C, 4DNFIOSQFOPV) | 126 | 12% | 11% | 22% | 70 kb |

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

## Natural variation as experiments: the C5 probe (2026-09-28)

Item 13 C5 (docs/ROADMAP.md section 5) asks whether human haplotypes can
serve as experiments. Two haplotypes are two versions of the same regulatory
sequence. Where a genotype and a molecular readout are measured in the same
samples, an explanation from the attribution layer ("this element regulates
that gene") predicts an allelic difference that can be checked. VCFs alone
carry no outcome. This lane is the probe: what paired data is open, what it
costs, and how much of the attribution layer it reaches. It measures nothing
and registers nothing that measures. `scripts/c5_paired_variation_probe.py`
writes `data/results/c5_paired_variation_probe.json`, counts only.

**What is paired and reachable (read on the portals, 2026-09-28).**

| Resource | What is paired | Samples | Access | Licence | Size | Allele counts |
| --- | --- | --- | --- | --- | --- | --- |
| GTEx v8 allelic expression | RNA haplotype counts per gene per sample (phASER, with and without WASP) with WGS genotypes | 15,253 samples, 54 tissues, 838 donors | haplotype matrices open on the portal; SNP-level ASE, genotypes and read-back phasing in dbGaP phs000424.v8 (not requested); no haplotype expression in v10 or v11 | portal open access; article CC BY 4.0 | 461.5 to 555.0 MB per matrix | provided, per gene |
| Geuvadis E-GEUV-1 | LCL RNA-seq ASE per individual per site, with 1000 Genomes phased genotypes | 462 (EUR 373, YRI 89), LCL only | open, no registration | EMBL-EBI terms; IGSR: available without embargo | ASE table 918 MB (GRCh37); genotypes 621 MB to 3.5 GB per chromosome | provided (the file is named COV8; the README does not define the floor) |
| EN-TEx | haplotype read counts for ATAC, histone ChIP, CTCF, POLR2A and EP300 ChIP, and RNA-seq, on each donor's own phased diploid genome | 4 donors, 18 to 23 tissues each with both chromatin and RNA | open, "fully open-consented and accessible without registration"; the donors' personalized genomes on the ENCODE portal, 1.7 to 3.6 GB each | ENCODE data use policy (no restrictions); article CC BY 4.0 | `cCREs_default_AS.tsv` 711 MB, `genes_default_AS.tsv` 94 MB, per-SNV table 2.5 GB | provided (AlleleSeq2) |
| ADASTRA Mabel v6.1 | TF ChIP-seq allelic reads at SNVs called from the reads | 1,073 TFs, 649 cell types, 15,970 alignments | open, Zenodo 14174114 (the site's downloads page answered HTTP 502) | CC BY 4.0 | zip 942.7 MB; the four context cells 26 MB | aggregated per SNV per cell type |
| UDACHA IceKing v1.0.3 | DNase, ATAC and FAIRE allelic reads | 5,858 datasets | open, Yandex Disk | none stated for the data; article CC BY-NC-ND 4.0 | zip 826.5 MB; the four context cells 12 MB | aggregated per SNV per cell type |
| AlleleDB (2016) | ASB and ASE on personal genomes | 382 individuals | open | not stated | 0.2 to 96 MB | provided; GRCh37 |
| ENCODE4 Hi-C genophasing | phased variant calls for 43 biosamples, HepG2 and IMR-90 among them; no molecular readout | 43 | open | ENCODE, no restrictions | 122 to 198 MB per VCF | none: from BAMs |
| GIAB HG002 RNA-seq | Illumina, PacBio and ONT RNA on three HG002 cell stocks, with the Q100 assembly phased by parent | 1 individual | open | NIST (17 USC 105) | 7.3 to 9.1 GB per mRNA BAM | none: computable by range reads |

GTEx is the largest and the one that cannot be used. Its open matrix gives
two haplotype counts per gene and nothing about which haplotype carries an
element's allele; that is in the controlled genotypes.

**The element universe.** The all-element archive holds 961,227 elements.
440,377 of them name a predicted coding target (108,599 "strong"); those
element-gene pairs are the explanations a paired readout could check. Every
target symbol is in GENCODE v50. NA12878, HG002 and Geuvadis were all read
from lymphoblastoid lines, so the probe also counts the LCL context:
121,151 elements overlap a GM12878 H3K27ac peak (ENCFF361XMX), 14,174 genes
have a GTEx v8 median TPM of 1 or more in EBV-transformed lymphocytes, and
58,737 pairs meet both. Neither filter is an allelic quantity.

**Overlap, label-blind.** A pair is readable in a genome when its element
carries a heterozygous SNV (otherwise both haplotypes hold the same sequence
and the explanation predicts no difference) and its target's exons carry one
(otherwise expression cannot be split by haplotype). From the allelic tables
only identifier, position, donor, tissue and assay columns were read, never a
count, ratio or significance.

| Resource | Elements reached | Pairs reached |
| --- | --- | --- |
| NA12878, GIAB HG001 v4.2.1 | 142,700 carry a heterozygous SNV (14.85%) | 59,185 heterozygous at both; 7,155 of them in the LCL context (3,143 targets) |
| HG002, Q100 dipcall | 152,427 (15.86%) | 63,645 (63,515 autosomal); 7,681 in the LCL context (3,330 targets) |
| Geuvadis, 449 samples, chr21 only | 11,123 of 12,139 heterozygous in at least one sample, 8,428 in 10 or more; median 38 samples | 3,608 of 5,174 contrastable (at least 10 samples heterozygous at both, and at least 10 with a readable target and a homozygous element); 395 of those in the LCL context (83 targets) |
| EN-TEx, 4 donors | 207,479 accessible in at least one chromatin assay (21.6%); 75,241 to 88,820 per donor | 85,571 with element and target both accessible in the same donor and tissue (19.4%; 23,508 strong); 31,384 to 36,744 per donor |
| ADASTRA, eligible SNVs | GM12878 47,673 (4.96%), K562 51,787, HepG2 48,323, IMR-90 4,058 | none: binding reads the element, never the gene |
| UDACHA, eligible SNVs | DNase: GM12878 16,101, K562 15,171, HepG2 7,931, IMR-90 12,727; ATAC 2,591 to 10,960 | none: accessibility only |

Reached is not powered: whether a reached pair has enough reads is in the
count columns this probe did not read. EN-TEx matches its V2 cCRE regions to
the archive's V3 elements by overlap (198,063 of its 250,722 chromatin
regions overlap an element). The panel drops singletons (MAC 2 or more), and
the GIAB benchmark VCFs cover their confident regions only (854,568 and
940,738 elements lie inside them), so the genome counts are floors. The
UDACHA release carries no readme; that its files list every coverage-passing
SNV and not only the significant ones is not stated by the release.

**Bias controls a build would need, for EN-TEx and Geuvadis.**

- *Reference mapping bias.* EN-TEx maps every read to both of the donor's
  haplotypes (AlleleSeq2) and filters ambiguous mappings, so the reference
  allele is not favoured by construction. The Geuvadis table was mapped to
  GRCh37; its REF_RATIO column is a per-sample, per-allele-pair null, which
  absorbs the average bias but not a site's own. A build drops sites in the
  project's low-mappability track and treats any site that decides a result
  with a WASP-style remap or the panel's personal haplotypes.
- *Imprinted and monoallelic genes.* Among the readable pairs, 684
  (NA12878) and 716 (HG002) have a target geneimprint lists as imprinted,
  407 and 389 one it lists as predicted, and 118 and 84 an immunoglobulin or
  HLA gene, which a clonal B-cell line expresses from one rearranged allele.
  On chr21, 91 of Geuvadis's 3,608 contrastable pairs have an imprinted
  target, none of them in the LCL context. X-linked targets are left out (X
  inactivation in the female line, one X in HG002). HG002's Q100 haplotypes
  are labelled by parent, and EN-TEx publishes the parental origin of each
  phased block (`phased_block.tar.gz`), so parent-of-origin expression can
  be separated in both.
- *Linked variants.* This is the control that decides what a result can
  say. Within 1 Mb either side of a readable target's TSS, one genome
  carries a median of 1,580 (NA12878) and 1,809 (HG002) heterozygous SNVs,
  and a Geuvadis sample 2,012 on chr21. An allelic imbalance belongs to the
  whole haplotype; in one person it cannot be credited to the element. So
  EN-TEx's four donors support a coordination test (does the target's
  expression lean to the haplotype whose copy of the element is the more
  active, more often than a non-target gene at matched distance?), not an
  attribution. The Geuvadis contrast, samples heterozygous at the element
  against samples homozygous at it with the target readable in both, can
  condition on other variants; that is why the probe asks for 10 samples on
  each side.
- *Low counts.* One binomial test at α 0.05 with 80% power needs 47 allelic
  reads to tell 0.70 from 0.50, 85 for 0.65, 194 for 0.60 and 783 for 0.55
  (normal approximation; overdispersion raises each). The Geuvadis file is
  named COV8, which reads as a floor of 8 reads. Counts are summed over a gene's heterozygous exonic sites per
  haplotype, as phASER does, and tested beta-binomially, as EN-TEx does.
- *Phase.* The element's allele and the exon's allele must sit on known
  haplotypes. HG002 is phased by parent genome-wide and EN-TEx's genomes
  with long reads. The 1000 Genomes panel is statistical phasing with a
  pedigree correction (SHAPEIT2-duohmm, its README), whose switch errors grow
  with distance. GIAB HG001 v4.2.1 writes no phase (0 of 2,028,130
  heterozygous SNVs), so NA12878's would come from the panel or from the
  NA12878 personal genome on the EN-TEx portal.
- *Shared provenance.* GTEx eQTLs already enter the attribution layer
  (`eqtl_targets`), so an allelic readout from GTEx donors and reads would
  not be independent of them. EN-TEx's four donors are GTEx donors, and the
  cCRE registry is built from ENCODE assays; whether registry V3 drew on
  EN-TEx datasets was not checked. C4's provenance field (R5) carries both.

**Recommendation.** Open, paired data does reach the element set: EN-TEx
reaches about a fifth of the elements and of the element-gene pairs, with
counts provided, on personal genomes, under no restriction. If the phase B
pilot helps, build C5 on EN-TEx first, as a within-donor coordination test
against matched non-target genes, and on Geuvadis with the GRCh38 panel
second, as the between-individual test in the LCL context, the only one of
the two that can address linked variants. GTEx stays out while its genotypes
are controlled. ADASTRA and UDACHA are the only open allelic data in the
archive's own four cells and can check the element side there, never the
gene. HG002's RNA can check a direction in one person.

**What was read.** Cached in the git-ignored `data/cache/c5/` (151 MB): the
GIAB HG001 v4.2.1 VCF and bed, the Geuvadis sdrf (sample ids only),
geneimprint's human table (symbols only) and GTEx v8 median TPM (the LCL
column only). Streamed and discarded, with the sha256 of the bytes that
passed in the manifest: 1000 Genomes chr21 (427 MB) and EN-TEx's two tables
(745 and 98 MB). Pulled from the release zips by HTTP range: ADASTRA 27 MB,
UDACHA 17 MB. No single download reached 1 GB; the script ran four times
while it settled, streaming each time. 0 AlphaGenome requests: the targets
come from the cached archive.

## EN-TEx allele coordination: the feasibility audit (audit C, checkpoint 1, 2026-09-29)

**The question.** Could EN-TEx support a test of whether a frozen nominated
target's expression leans to the same haplotype as its element's
accessibility more often than distance-, expression- and depth-matched
alternative genes do? That would test coordination, not causality. This is a
new question with its own rationale, not a continuation of the coherence
pilot (stopped at `197c560`). `scripts/entex_feasibility.py` (code
`4f40089`) writes `data/results/entex_feasibility.json`, stamped clean at
`612eb87`. It used depth and counts only, and no signed allelic outcome was
opened.

**Outcome discipline.** The two AS tables were streamed once through one
parsing function, `depth_only`. It keeps the identifier columns and sums
`hap1_count` and `hap2_count` on the line that parses them, keeping only
the sum. It reduces `imbalance_significance` (EN-TEx's 0/1 call) to an
unsigned flag. It never indexes `hap1_allele_ratio` or `p_betabinom`. Only
its output reaches the git-ignored cache (`data/cache/entex/depth_only.tsv.gz`,
sha256 `ce3a9364…6082`). Tests show that swapping the two haplotype columns
changes neither the parser's output nor the cache bytes. The unsigned flag
is reported as per-side marginals only. The joint count (both sides
imbalanced) was never formed, and no flag was computed for an alternative
gene. The result's `outcome_exposure_record` lists every file and column read.

**Sources (sha256).** `cCREs_default_AS.tsv` `679d916c…ba28` (745,298,142 bytes,
5,330,335 rows) and `genes_default_AS.tsv` `0e3f290f…6daeb` (97,973,665 bytes,
793,882 rows) are the same bytes the C5 probe hashed on 2026-09-28. The
other files are `phased_block.tar.gz` `270fc041…588b`,
`individual1_SV.vcf` `134c23fa…9eda`, `individual4_SV.vcf` `f9f4d811…f31f`
and the ENCODE exclusion list v2 `c92e763a…8cf2`.

| Checkpoint item | Finding |
| --- | --- |
| 1. Same donor and tissue, ATAC and RNA | 39 donor-tissues (ENC-001 5, ENC-002 9, ENC-003 11, ENC-004 14) over 22 tissues. Only 3 tissues are in all four donors (gastrocnemius medialis, sigmoid and transverse colon). 45 donor-tissues have ATAC and 89 have RNA. C5's "any chromatin assay" reach is not this number. |
| 2. Phase | All 181,025 ATAC element-target rows have the element and the target's exon span inside one published phased block. The blocks are arm-scale: 24 blocks over 1 Mb cover 2.806 to 2.808 Gb of each donor's autosomes. Switch errors inside a block are not reported and were not assessed. ENC-00k is matched to individual k by number only. |
| 3. Total depth per row | ATAC at the element: median 16 (IQR 12 to 29; the catalogue lists from 9). RNA at the gene: median 92 (IQR 32 to 282; from 8). The element side binds: 22,228 rows reach 47 ATAC reads, against 120,225 that reach 47 RNA reads. |
| 4. Balanced rows (EN-TEx's call, per side) | After exclusions, at any depth: element side 169,526 of 172,947 (98.0%) and gene side 162,342 (93.9%). At 47 reads: 13,780 and 13,240 of 14,538 (94.8% and 91.1%). |
| 5. Exclusions (floor 0; a row can have several) | Deletion under a target exon 2,811 and under the element 153, from EN-TEx's SV calls. These exist for ENC-001 and ENC-004 only, and hold DEL, INS and INV with no DUP. MHC 1,528; imprinted 1,473, predicted 914, other non-null imprinting statuses 310; exclusion list 1,129 (exon) and 387 (element). No IG or TR target remained. Rows fall from 181,025 to 172,947. |
| 6. Matched alternatives | Rows with at least 1, 3 and 5 alternatives matched within 2-fold on distance, GTEx v8 expression and RNA depth, on the element's block: 8,149, 60 and 5 at floor 0; 509, 2 and 0 at 47 reads. |

**Power sensitivity, total depth only (ATAC element side).** The
independent-locus count is taken on the parent universe. Target TSSs are
clustered by single linkage at 1 Mb over every element-target row formed,
before any filter. That gives 566 clusters over 10,243 target genes. A
filtered set counts the parent clusters it touches, and that count is only
descriptive. The first run recomputed the clusters after each filter, which
split connected components (566 became 781 at 47 reads) without adding
independent evidence. The review caught this, and the clusters are now
frozen. A second development run tried looser matchings. Those were
withdrawn unread, and nothing was relaxed after the counts were seen.

| Floor (reads, both sides) | Rows with depth | After phase and exclusions: rows / targets / clusters | With 3 matched alternatives: rows / targets / clusters / donors |
| --- | --- | --- | --- |
| listed (8 to 9) | 181,025 | 172,947 / 9,904 / 558 | 60 / 25 / 22 / 4 |
| 10 | 169,927 | 162,325 / 9,666 / 555 | 52 / 25 / 22 / 4 |
| 20 | 63,765 | 60,943 / 7,038 / 527 | 10 / 5 / 5 / 3 |
| 47 (tells 0.70 from 0.50) | 15,230 | 14,538 / 3,543 / 465 | 2 / 1 / 1 / 1 |
| 85 (0.65) | 4,294 | 4,093 / 1,628 / 371 | 0 |
| 194 (0.60) | 341 | 322 / 207 / 145 | 0 |
| 783 (0.55) | 0 | 0 | 0 |

With H3K27ac read at the element instead (secondary), 3 clusters remain at 47
reads.

**Checkpoint 1: no-go.** The rule was fixed in the script before the first
run. At 47 reads on both sides, after phase, exclusions and three matched
alternatives, it asked for 194 independent loci (the units a one-sample
binomial test needs to tell an agreement rate of 0.60 from 0.50) and for 3
donors with 30 loci each. One locus remains, in one donor (ENC-003). This is
the infeasibility of this matched-control design on EN-TEx, not a finding
that EN-TEx lacks usable allelic evidence. The step that binds is per-row
matching, and it is recorded here, not relaxed. A different control or a
measurement-error design is possible in principle. It would need its own
rationale and registration, and none is started.

**Exposure.** No EN-TEx allelic table or derivative is named in any tracked
file outside the C5 probe, its result and the docs. Of EN-TEx's 691
experiment accessions, two are in the project: ENCSR136ZQZ (testis H3K27ac)
and ENCSR611DJQ (testis H3K4me3), both from ENC-001. They are pinned in
`data/results/epigenome_manifest.json` as total-signal peaks and fold
change. The nominations (AlphaGenome deletion targets) do not read that
layer. The testis and ovary DNase files there are named only by file
accession, so whether they are EN-TEx's was not established.

AlphaGenome is the other route. The nominations come from AlphaGenome
deletion scores. A peer session saved a copy of AlphaGenome's track metadata
(2,841 tracks, `5157dcca…abf7`), read here without a model request. It holds
246 ENCODE tissue tracks named as an EN-TEx tissue (128 histone ChIP, 75
RNA-seq, 25 DNase, 18 ATAC) and 29 GTEx RNA tracks for EN-TEx tissues. The
metadata names no experiment or donor, so whether any of these tracks was
built from an EN-TEx donor cannot be established here. Whether the cCRE
registry V3 drew on EN-TEx experiments was not checked either: every query
to the ENCODE portal on 2026-09-29 timed out or answered 504.

**Independence: unresolved.** EN-TEx's four donors are GTEx donors. GTEx v8
eQTLs enter the attribution layer (`eqtl_targets`), and GTEx RNA-seq tracks
are among AlphaGenome's outputs. So a result on EN-TEx cannot be called an
independent validation of anything fitted on GTEx or by AlphaGenome.

**The later endpoint (proposal only, not registered).** The endpoint would
measure excess haplotype agreement for frozen nominated targets over
distance-, expression- and depth-matched alternative genes, clustered by
locus, with donor-specific sensitivity. It tests coordination, not
causality. A heterozygous SNV alone predicts neither an effect nor its
direction, and balanced RNA does not refute a target. It would need its own
rationale and a registration before any outcome is opened. At this
checkpoint the matched design has one locus, so it is not proposed for
registration as it stands.

**Cost.** The cache build took 67.7 s wall, 14.9 s CPU and 119 MB peak, and
streamed 843,271,807 bytes in 2 requests. The committed run took 75.3 s
wall, 74.6 s CPU and 2.15 GB peak, with 0 bytes downloaded. All six runs
together took 497 s wall and 429 s CPU. The script made 6 requests for
965,991,460 bytes. By hand while designing, 10 requests fetched 2,145,198
bytes. The script's own design-time record says 8 requests and gives the
phased-block archive as 5,388 bytes. Those figures are corrected here:
correcting them in the script would have meant rewriting committed lines.
The total is about 968 MB, under the 3 GB bound. 0 AlphaGenome requests;
nothing paid or controlled-access.

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

## The result registry: what enters it and what is quarantined (item 12 S6, 2026-09-28)

Registered before the build. The second external review (ROADMAP section 5, item 12, S6) confirmed
two defects in the contract built for review R9 (`fe0880a`): `save_result` wrote a new result
before refusing it, and decided "historical" by `not p.exists()`, so a retry found the file it had
just written and only warned; and the revision stamp ran `git status --untracked-files=no`, so a
result could record a clean revision while an untracked script wrote it.

**Census before the change** (`0e50036`, `scripts/manifest_census.py --enforcement`): 221
`save_result` call sites in 148 files; 1,015 results on disk, 955 of them historical at `fe0880a`
(656 tracked in its parent, 299 git-ignored `reader_*_chr*.json` written before it); all 60 names
added since carry a complete manifest, in all 69 of their committed versions, so nothing committed
shows the retry path used. 13 writers put or edit results in `data/results/` without
`save_result`; the contract does not see them and this item does not change them.

**What enters `data/results/`.** A write reaches the registry only when its manifest is complete,
or when its name is on the legacy allowlist (it then warns, as before). Anything else is a failed
new result: it is written to the quarantine and `ManifestError` names that path. `strict=None` in
the registry means exactly "not on the allowlist"; whether the file already exists no longer
matters. `strict=True` enforces in any directory; `strict=False` cannot admit a name that is not
on the allowlist into the registry. Outside the registry (tests, scratch) an incomplete manifest
still only warns unless `strict=True`. A later complete write of the same name removes its
quarantined copy.

**The quarantine** is `data/quarantine/results/` (the rule is `<results_dir's parent>/quarantine/
<results_dir's name>`, so a test's directory gets its own). It is **git-ignored**, by a
`.gitignore` holding `*` that the writer creates in `data/quarantine/`: a quarantined file is a
computation that failed its contract, and tracking it would let the next `git add` publish
unvalidated numbers; its reason is per-machine run state, like `data/jobs/`. The self-ignoring
directory needs no edit to the shared root `.gitignore` and holds in worktrees and test
directories. Each file keeps the whole result plus a `quarantine` block: reason, the problems, the
path it was meant for, and when.

**The legacy allowlist** is `data/results_legacy.txt`: the 955 names historical at `fe0880a`, each
with its source (`tracked` in the tree of `fe0880a^`, or `ignored-local`: a git-ignored
`data/results/*.json` on this disk last written before `fe0880a`). `scripts/manifest_legacy.py`
generated it once from git history and refuses to overwrite it. No code ever adds to it; a name is
added only by a reviewed commit that says why. If it cannot be read, it counts as empty, so every
incomplete write in the registry is quarantined (fails closed).

**The revision stamp** reads `git status --porcelain -z --untracked-files=all`. An untracked file
counts as code when it sits under `genomeos/`, `scripts/` or `tests/` (a writer,
`scripts/grn_clamp_census.py`, imports `tests/test_grn_clamp.py`), or is a BioLang program
(`*.bio`) anywhere, since writers execute `data/organisms/*.bio` and `data/demo/*.bio`. Such files
are listed under `code.untracked_code_paths`, included in `code.dirty_code_paths`, and make
`code.dirty` true. **Untracked data files do not count and are not listed**: what a result read is
pinned by its manifest's `inputs` (path and sha256), which is the contract's mechanism for data;
counting them would mark dirty any result written while a peer has an unrelated file open in this
shared checkout, and the git-ignored stores (`data/reference`, `data/knowledge`, `data/cache`) are
invisible to `git status` anyway. An untracked file under `data/results/` is an earlier writer's
output and is recorded apart under `dirty_result_paths`, as a modified one already is.

**Predicted census after the build**: the same 1,015 results, 71 carrying a manifest, 67
complete (nothing is rewritten); 955 names on the allowlist, 60 on disk not on it, all complete;
0 quarantined. **Falsifiers**: the list generated from history is not 955 names or differs from
the census's historical set; or an existing test fails other than the two that assert the defect
(`test_a_new_result_without_the_contract_is_written_then_refused` and the retry half of
`test_strict_is_the_default_for_a_new_name_in_the_registry_only`), which would mean a legacy
writer was broken.

**Result, as predicted.** `scripts/manifest_legacy.py` generated 955 names from history (656
tracked, 299 ignored-local); `--check` finds 0 differences from the census's historical set. The
census after the build: the same 1,015 results, 71 carrying a manifest, 67 complete (nothing
rewritten); 955 on the allowlist, 60 not, 0 of those incomplete; 0 quarantined. The 76 existing
manifest tests pass, two of them rewritten from asserting the defect to asserting the fix, and
`tests/test_manifest_enforced.py` adds 19. In this shared checkout at build time the stamp named
three untracked files peers had not yet committed (one script, two tests) and marked the revision
dirty; before, it said nothing about them.

### Every writer through `save_result` (item 12 S6 follow-up, lane-contract, 2026-09-29)

**Census** (`d6a77dd`, `scripts/results_writer_census.py`, a static taint analysis of every Python file
under `genomeos/` and `scripts/`, replacing lane-s6's pattern): 36 sites in 28 files wrote into
`data/results/` without `save_result`. lane-s6's 13 writers were a lower bound by seven: `cmd_signals`,
`tf_atlas.save_cells`, `storage._distil_celegans_lineage`, `phosphosite_observation`,
`reader_depth_family` (three results), `trio_q100_population` and `variation_rerun_chain`'s backup. Three
registry files have no code writer at all (`alphagenome_rs12740374`, `manifest_headlines`,
`clinvar_chr21_coding.vcf.gz`); no static check can see a hand edit.

**Converted.** The 17 direct writers and the backup now go through `save_result` with a complete
manifest (sources, inputs with sha256, assembly, coordinates, parameters, exclusions, partitions).
Shared pieces: `genomeos/organism/provenance.py` (the worm's two sources and common inputs),
`manifest.files_entry` (one digest over a set of files read together, named relative to the checkout),
and `save_result(..., compact=True)` for the three per-cell and per-site tables that were written compact.
The two in-place editors became results of their own, because an edit inside a result sits under a
manifest that describes another run and rewriting the file through `save_result` would replace its
committed code stamp: `crispri_hct116` (its block in `crispri_published.json` is kept as committed; that
block's own manifest lacked `sources` and `exclusions`, which no check read while it lived inside another
result, and `save_result` refused it until they were stated) and `clause2_statistical_corrections` (the
four blocks lane-design added under `corrections_after_the_statistical_review_2026_09_28` stay as
committed). `variation_rerun_chain` keeps its backup in `data/cache/variation_rerun/`. No committed result
was rewritten: the legacy names regenerate with their manifest the next time their writer runs.

**No figure moved.** Each converted writer ran at the base commit and with the change, in two clean
worktrees, and every value except `date`, timing keys and the manifest was compared: 0 differences in
every output (22; the table is in the commit message). The only additions are the `result` and `date`
keys `save_result` writes. 0 AlphaGenome requests (no key in the worktrees; the HCT116 ledger unchanged).
Two side findings, before and after alike and not changed here: `celegans_commitment`'s competence nulls
depend on Python's string hash seed (12 values move in the fourth decimal between `PYTHONHASHSEED` 0 and
1), and `celegans_fate_reads` no longer reproduces its committed file (82 values) at this base.

**Out of scope, decided.** Ten sites write fetched reference data into `data/results/`
(`superdups_`, `rmsk_`, `ccres_`, `ccres_mm10_`, `dnase_*_<chrom>.bed.gz`, `gencode_v50_*.gff3.gz`, the two
orthology tables, `HG002_chr21.vcf`) and one deletes `rmsk_` files. They are not results: none is JSON,
the registry lists `*.json` only, and a manifest is a key inside a result's JSON. They do not belong under
`data/results/`; moving them is a separate change (their readers and the tracked chr21 subsets move with
them) and is not done here.

**The guard** is `tests/test_results_writers_guard.py`: it fails on any site under `genomeos/` or
`scripts/` that writes into `data/results/` other than through `save_result`, unless it is on one of two
lists with a reason each: two exceptions (`manifest_rebuild` copies inside its disposable second
checkout; `activator_combination --out` is a false positive of the flow-insensitive analysis) and the
eleven fetched-data sites. An entry that no longer matches a site fails too. The analysis is tested on
each pattern the census found and on reads that must not be flagged; a planted direct write fails it.
