Three tumours whose driver is not a point mutation, used as controls for the
copy-number and structural-variant path. Each VCF carries a variant in a
*different* gene, so the gene under test can reach the candidate list only
through its copy-number or structural call. Run the same inputs against a
build without that path and the gene under test does not appear at all: that
is the control these fixtures exist for.

| input | the driver | the VCF's variant | what it tests |
|---|---|---|---|
| `erbb2_amplification_only.vcf` + `.cnv` | ERBB2 at 12 copies | PIK3CA | a surface target with an approved antibody, reachable only through the amplification |
| `cdkn2a_deleted.vcf` + `.cnv` (GISTIC) | CDKN2A homozygous loss, EGFR amplified | BRAF V600E | a deleted gene must not be offered as a target: no product exists to bind. The pairing is the glioma pattern |
| `alk_eml4_fusion.vcf` + `.sv` | EML4-ALK | PIK3CA | a rearrangement, with no copy-number call and no coding change |
| `fgfr3_tacc3_urothelial.vcf` + `.sv` | FGFR3-TACC3, one patient's record | ARID1A R1287T and TSC1 R190P, the same tumour's own | a caller's own export, read by its column names, whose 5' partner keeps its outward face: the one fusion case where the product's extent decides nothing and has to stay unknown |

Three of the four are written by hand. The fourth, `fgfr3_tacc3_urothelial`, is one
public tumour's record: sample TCGA-CF-A3MG-01 of `blca_tcga_pan_can_atlas_2018`
(TCGA Bladder Urothelial Carcinoma, PanCancer Atlas) as the cBioPortal public REST API
returns it, in that export's own column names and its own GRCh37 coordinates. Its VCF
carries two of that sample's 125 somatic mutations, with the allele fractions its read
counts give, lifted to GRCh38 through the Ensembl assembly map because the other
fixtures here are GRCh38; the lifted positions were checked against the record's
reference alleles and annotate to the protein changes the record names. FGFR3 is not
mutated in this sample, which is what makes the structural call the only route to it.
