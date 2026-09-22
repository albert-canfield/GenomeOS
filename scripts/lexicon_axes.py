# SPDX-License-Identifier: AGPL-3.0-or-later
"""The lexicon re-asked at base resolution (phyloP per base, the HPRC panel per base, JASPAR at 0.95
with column-permuted decoys). No model calls; phyloP and Gnocchi are streamed once per chromosome
and kept as byte tracks under data/knowledge/lexicon.

    uv run python scripts/lexicon_axes.py chr21 chr22
"""

from genomeos.attribution.lexicon_axes import main

if __name__ == "__main__":
    main()
