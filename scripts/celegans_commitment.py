# SPDX-License-Identifier: AGPL-3.0-or-later
"""Commitment, competence and lateral inhibition tested in the worm's measured factors (area E).

Reads the Ma 2021 presence table (data/results/celegans_tf_atlas_cells.json), the local atlas levels
(data/knowledge/celegans/atlas_levels.json, built by scripts/celegans_fates.py) and the reference lineage,
runs the four tests of genomeos.organism.commitment and the lateral-inhibition reference model, and writes
data/results/celegans_commitment.json. No model calls; a few minutes on one core.

    uv run python scripts/celegans_commitment.py [--permutations 200]
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import date
from pathlib import Path

from genomeos.organism import commitment as cm
from genomeos.organism import lateral
from genomeos.organism.atlas_levels import load_levels, load_strain_peaks
from genomeos.organism.reference import ReferenceLineage
from genomeos.organism.tf_atlas import CELLS_FILE, load_cells

RESULT = Path("data/results/celegans_commitment.json")


def _clean(x):
    """NaN (a bin with nothing to measure) is written as null."""
    if isinstance(x, float) and math.isnan(x):
        return None
    if isinstance(x, dict):
        return {k: _clean(v) for k, v in x.items()}
    if isinstance(x, list | tuple):
        return [_clean(v) for v in x]
    return x


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--permutations", type=int, default=200)
    args = ap.parse_args()
    n = args.permutations
    ref = ReferenceLineage.load()
    atlas = load_cells()
    lifetimes = json.loads(CELLS_FILE.read_text())["lifetimes"]
    levels = load_levels()
    cells = cm.atlas_cells(ref, atlas, lifetimes)
    names = sorted({f for fs in atlas.values() for f in fs})
    programme_factors = {f for fs in cm.PROGRAMMES.values() for f in fs}
    out = {
        "result": "celegans_commitment",
        "date": date.today().isoformat(),
        "source": "Ma et al. 2021 atlas (Zenodo 4737593, CC BY 4.0) x WormWeb lineage (Sulston 1983)",
        "cells_with_complete_lifetimes": len(cells),
        "states": cm.states_test(cells, permutations=n),
        "programmes": cm.programme_test(cells, names, permutations=n),
        "ratchet": cm.ratchet_test(cells, ref, permutations=5 * n),
        "competence": {
            "all_factors": cm.competence_test(cells, ref, permutations=n),
            "all_factors_same_generation": cm.competence_test(cells, ref, permutations=n, by_generation=True),
            "programme_factors": cm.competence_test(
                cells, ref, factors=programme_factors, permutations=5 * n
            ),
        },
        "sisters": cm.sister_test(cells, ref, load_strain_peaks(), levels, permutations=10 * n),
        "lateral_inhibition_model": lateral.selection_statistics(runs=200),
    }
    out = _clean(out)
    RESULT.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    print(
        json.dumps({k: v for k, v in out.items() if k in ("states", "programmes", "ratchet")}, default=str)[
            :3000
        ]
    )


if __name__ == "__main__":
    main()
