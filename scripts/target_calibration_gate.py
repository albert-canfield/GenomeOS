# SPDX-License-Identifier: AGPL-3.0-or-later
"""The calibration's top-target gate: counted first, then the calibration re-fitted without it.

    uv run python scripts/target_calibration_gate.py --census   # the counts only, nothing re-fitted
    uv run python scripts/target_calibration_gate.py            # the census and both arms

`target_calibration` reads the compact all_elements table, which keeps one gene per element, so its
`deletion_drop` is a structural zero for every pair whose measured gene is not that gene and its
predicted-target fit is 245 pairs of 8,796. `attribution/targets.ElementResponses` carries every gene
in the scorer's 1 Mb window with a signed change on the cell's own track, so the same pairs can be
asked about their own gene. `--census` says how many pairs the gate admits, how many it excludes and
which of the excluded ones the sweep did answer; the full run adds the re-fit, with the shipped arm
beside it as the control.

Reads the cached benchmark tables (data/knowledge/crispri), the all-enhancer deletion table
(data/knowledge/alphagenome/all_elements), the per-element response cache
(data/knowledge/alphagenome/elements, 775 MB, one chromosome held at a time) and GENCODE v50. No
model request, no network. Saves data/results/target_calibration_gate.json.
"""

from __future__ import annotations

import argparse
import json
import time
from typing import Any

from genomeos.attribution import crispri
from genomeos.attribution import target_calibration as tc
from genomeos.attribution import targets as tg
from genomeos.results import save_result


def census(
    training: list[crispri.Pair], heldout: list[crispri.Pair], table: crispri.DeletionTable
) -> dict[str, Any]:
    """The gate, counted on the populations the published calibration is fitted and judged on."""
    crispri.annotate(training, table)
    crispri.annotate(heldout, table)
    tc.add_features(training, table)
    tc.add_features(heldout, table)
    train = tc.scored(training)
    held = tc.scored([p for p in heldout if p.cell == tc.CELL])
    responses = tg.ElementResponses()
    return {
        "preregistered": tc.PREREGISTERED_GATE,
        "what": (
            "the pairs the top-target gate admits and excludes, before anything is re-fitted, with "
            "the sweep's own per-element cache asked about every excluded pair's own gene"
        ),
        "training, every feature present": tc.gate_population(train),
        "held-out K562, every feature present": tc.gate_population(held),
        "training, on a deleted element": tc.gate_population([p for p in training if p.covered]),
        "held-out K562, on a deleted element": tc.gate_population(
            [p for p in heldout if p.cell == tc.CELL and p.covered]
        ),
        "what the window reader says, training": tc.gate_silences(train, table, responses),
        "what the window reader says, held-out K562": tc.gate_silences(held, table, responses),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", action="store_true", help="count the gate; re-fit nothing")
    args = ap.parse_args()
    t0 = time.time()

    training, heldout = crispri.load(crispri.TRAINING), crispri.load(crispri.HELDOUT)
    table = crispri.DeletionTable()
    out = census(training, heldout, table)
    for key in ("training, every feature present", "held-out K562, every feature present"):
        print(f"{key}: {json.dumps(out[key])}")
    for key in ("what the window reader says, training", "what the window reader says, held-out K562"):
        print(f"{key}: {json.dumps(out[key]['by_side'])}")
        print(f"    rank: {json.dumps(out[key]['measured_gene_rank_in_the_window'])}")
    if args.census:
        print(f"census only, nothing re-fitted ({time.time() - t0:.0f} s)")
        return 0

    path = save_result("target_calibration_gate", out)
    print(f"saved {path} ({time.time() - t0:.0f} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
