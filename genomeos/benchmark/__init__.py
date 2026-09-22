# SPDX-License-Identifier: AGPL-3.0-or-later
"""Benchmarks: the project's machinery against answers that are already known.

A genome-wide number is only worth what the same machinery produces at the few
places biology has already settled. This package holds those coherence checks.
`loci.py` is the known-locus panel (docs/LOCI-BENCHMARK.md): a dozen loci whose
target gene, tissue, causal variant, direction and functional class are published,
read blind through the project's own layers, with every hit labelled by where it
came from (a measurement the pipeline read, or an annotation that already carried
the answer) and with matched windows of no known function scored identically.

The therapeutic benchmark (`scripts/therapeutic_benchmark.py`,
`tests/test_therapeutic_benchmark.py`) is the same idea one lane over.
"""
