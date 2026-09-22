A synthetic B-cell lymphoma, for the expression-driven route.

The VCF carries one driver mutation (TP53). The RNA table is where the targets
are: CD19 and MS4A1 (CD20) high, TNFRSF17 (BCMA) moderate, a few housekeeping
genes and several membrane proteins at healthy-tissue levels so the scan has
something to reject. None of these numbers are from a real patient; they are
in the range a B-cell malignancy would give, and every conclusion drawn from
them is a demonstration of the route, not a finding.

    genomeos therapeutic --tumour data/demo/expression/bcell_lymphoma.vcf \
      --rna data/demo/expression/bcell_lymphoma_rna.tsv --scan 6
