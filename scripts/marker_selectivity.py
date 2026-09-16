# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Marker logic against healthy tissue: A, A and B, A and not C.

Scores Boolean marker gates over the packaged Human Protein Atlas surface
universe: tumour coverage from DepMap cell-line RNA per lineage, healthy burden
from the atlas's essential organs, and the pre-registered benchmark that asks
whether gates containing an approved or clinical target outrank random ones.

Run:     uv run python scripts/marker_selectivity.py
Offline: uv run python scripts/marker_selectivity.py --offline   (cache only)
Writes:  data/results/marker_selectivity.json
"""

from __future__ import annotations

import sys
import time

from genomeos.results import save_result
from genomeos.therapeutics import selectivity
from genomeos.therapeutics.depmap import DepMapUnavailableError


def main() -> int:
    net = "--offline" not in sys.argv
    log = sys.stdout
    t0 = time.time()
    print(f"marker selectivity, network={'on' if net else 'off'}", file=log)
    try:
        result = selectivity.run(net=net, log=log)
    except DepMapUnavailableError as e:
        print(f"DepMap unavailable: {e}", file=sys.stderr)
        return 2
    result["seconds"] = round(time.time() - t0, 1)
    save_result("marker_selectivity", result)
    bench = result["benchmark"]
    print(
        f"\n{result['combinations_scored']} gates scored over {result['data']['cell_lines']} cell lines "
        f"in {len(result['data']['lineages'])} lineages",
        file=log,
    )
    null = bench["placebo_null"]
    print(
        f"pre-registered claim {'supported' if bench['supported'] else 'NOT supported'}: benchmark gates "
        f"median percentile {bench['median_percentile_benchmark']} against "
        f"{bench['median_percentile_random']} for random gates, p={bench['p_one_sided']:.3g}, "
        f"delta={bench['cliffs_delta']}",
        file=log,
    )
    print(
        f"placebo null over {null['gene_sets_drawn']} random gene sets of {null['set_size']}: "
        f"mean delta {null['mean_delta']}, 95th percentile {null['p95_delta']}, max {null['max_delta']} "
        f"({'flat' if null['flat'] else 'NOT FLAT - read the result against this, not against zero'})",
        file=log,
    )
    for row in result["top_combinations"][:15]:
        gate = " and ".join(row["required"]) + (
            "".join(f" and not {g}" for g in row["excluded"]) if row["excluded"] else ""
        )
        print(
            f"  {gate}: {row['lineage']} coverage {row['coverage']:.0%}, "
            f"healthy {row['healthy_burden_ntpm']:g} nTPM, score {row['score']}",
            file=log,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
