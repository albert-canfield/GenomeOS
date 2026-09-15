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
