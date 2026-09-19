"""The DNA → RNA → protein relationship, as one traceable object.

This is the home of the central-dogma relationship in GenomeOS. `trace()`
takes a transcript and the genome and returns a CentralDogmaTrace that maps
every genomic base to its position in the mature mRNA, its codon and its
residue, and back. `genomeos flow` and the Flow tab of the web UI use it;
so does the variant classifier when it explains a coding change.
"""

from genomeos.flow.trace import CentralDogmaTrace, trace, trace_gene

__all__ = ["CentralDogmaTrace", "trace", "trace_gene"]
