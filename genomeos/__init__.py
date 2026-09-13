"""GenomeOS: an executable model of biology.

Layers:
    genomeos.genome   - lossless genome engine (FASTA, sequences, loci)
    genomeos.ir       - BioIR: the intermediate representation
    genomeos.lang     - BioLang: source language compiled to BioIR
    genomeos.runtime  - BioVM: engines that run BioIR through time
"""
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.

from genomeos.version import __version__  # the engine owns it; re-exported here for compatibility

__all__ = ["__version__"]
