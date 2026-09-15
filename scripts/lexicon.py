# SPDX-License-Identifier: AGPL-3.0-or-later
"""Build the lexicon index for one or more chromosomes (no model calls, nothing streamed).

The full index goes to data/knowledge/lexicon/lexicon_<chrom>.json.gz (git-ignored); the summary,
with the two questions answered and every test's null beside it, is committed as
data/results/lexicon_<chrom>.json.

    uv run python scripts/lexicon.py chr21 chr22
"""

from genomeos.attribution.lexicon import main

if __name__ == "__main__":
    main()
