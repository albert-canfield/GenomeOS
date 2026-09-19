# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Distil the Human Protein Atlas's healthy-tissue RNA for the surface universe.

Two requests, about ten seconds, and the raw download is discarded. What is
left is a table shipped inside the package, so every candidate the therapeutic
pipeline scores has a healthy-tissue answer with no network at all, and a scan
over the whole membrane class becomes possible instead of five thousand
requests.

Run:    uv run python scripts/normal_tissue.py
Writes: genomeos/therapeutics/data/normal_tissue.json.gz (the table)
        data/results/normal_tissue_atlas.json           (the record of the run)
"""

from __future__ import annotations

import sys
import time

from genomeos.results import save_result
from genomeos.therapeutics.atlas import PACKAGED, distil, summary, write
from genomeos.therapeutics.providers import CRITICAL_TISSUES


def main() -> int:
    log = sys.stdout
    t0 = time.time()
    print(f"Human Protein Atlas: {len(CRITICAL_TISSUES)} tissues", file=log, flush=True)
    table = distil(CRITICAL_TISSUES, log=log)
    path = write(table)
    s = summary(table)
    save_result("normal_tissue_atlas", s)
    print(
        f"\n{s['genes']} genes, {s['genes_with_tissue_values']} with tissue values, "
        f"{path.stat().st_size / 1024:.0f} KB gzipped, {time.time() - t0:.1f}s",
        file=log,
    )
    for gene, ex in s["examples"].items():
        high = ex["highest_queried_tissue"]
        print(
            f"  {gene:<10} {ex['specificity'] or 'unclassed':<24}"
            + (f" highest {high['tissue']} at {high['nTPM']:g} nTPM" if high else " no tissue values"),
            file=log,
        )
    return 0 if PACKAGED.exists() else 1


if __name__ == "__main__":
    raise SystemExit(main())
