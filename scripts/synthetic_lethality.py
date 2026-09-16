# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Synthetic lethality across DepMap: does a loss make a dependency?

Streams the DepMap 24Q4 Public CRISPR gene-effect matrix and the matching
mutation, copy-number and MSI calls, tests every pair in the pre-declared
candidate space with a one-sided rank test, corrects with Benjamini-Hochberg,
and reports the seven pre-registered controls and a permuted-label null.

Nothing is stored but the distilled per-gene cache under
data/knowledge/therapeutics/depmap (gzipped, requested genes only); the raw
matrices are streamed and discarded.

Run:     uv run python scripts/synthetic_lethality.py
Offline: uv run python scripts/synthetic_lethality.py --offline   (cache only)
Writes:  data/results/synthetic_lethality.json
"""

from __future__ import annotations

import sys
import time

from genomeos.results import save_result
from genomeos.therapeutics import synthetic_lethality as sl
from genomeos.therapeutics.depmap import DepMapUnavailableError


def main() -> int:
    net = "--offline" not in sys.argv
    log = sys.stdout
    t0 = time.time()
    print(f"synthetic lethality over {sl.depmap.RELEASE}, network={'on' if net else 'off'}", file=log)
    try:
        result = sl.run(net=net, log=log)
    except DepMapUnavailableError as e:
        print(f"DepMap unavailable: {e}", file=sys.stderr)
        return 2
    result["seconds"] = round(time.time() - t0, 1)
    save_result("synthetic_lethality", result)
    print(
        f"\n{result['tests']} pairs tested over {result['data']['cell_lines_with_gene_effect']} "
        f"cell lines; {result['hits']} at q <= {sl.FDR}; "
        f"{result['controls_recovered']}/{len(result['controls'])} controls as pre-registered; "
        f"permuted null mean {result['shuffled_null']['mean_hits']} hits",
        file=log,
    )
    for c in result["controls"]:
        res = c["result"]
        print(
            f"  {c['lost']} {c['loss']} -> {c['dependency']}: {c['verdict']} "
            f"(expected {c['expect']}; n={res.get('lines_lost')} lost vs {res.get('lines_intact')}, "
            f"delta={res.get('cliffs_delta')}, p={res.get('p')}, q={res.get('q_within_controls')})",
            file=log,
        )
    for hit in result["top_hits"][:15]:
        print(
            f"  hit: {hit['lost']} loss -> {hit['dependency']} q={hit['q']} "
            f"delta={hit['cliffs_delta']} ({hit['lines_lost']} lost lines)",
            file=log,
        )
    return 0 if result["controls_recovered"] == len(result["controls"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
