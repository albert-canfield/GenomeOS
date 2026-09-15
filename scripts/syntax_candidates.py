# SPDX-License-Identifier: AGPL-3.0-or-later
"""Read the organiser's syntax candidates one by one from local evidence (no model calls).

Writes data/results/syntax_candidates_genome_wide.json. About ten minutes: GENCODE, DNase peaks,
the registry, conserved elements and repeats per chromosome, the segment parser over every block
and its control windows, and one pass over the local UniProt proteome.

    uv run python scripts/syntax_candidates.py
"""

from genomeos.attribution.candidates import main

if __name__ == "__main__":
    main()
