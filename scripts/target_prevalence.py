# SPDX-License-Identifier: AGPL-3.0-or-later
"""The prevalence term: a per-screen intercept, fitted only after PREREGISTERED_PREVALENCE was fixed.

    uv run python scripts/target_prevalence.py --census   # the screen census alone, nothing fitted
    uv run python scripts/target_prevalence.py            # the census, both reads and the genome

Reads the cached benchmark tables (data/knowledge/crispri), the all-enhancer deletion table
(data/knowledge/alphagenome/all_elements), the registry classes (data/results/ccres_*.bed.gz) and
GENCODE v50 (data/reference), and saves data/results/target_prevalence.json. No model request, no
network, and data/results/target_calibration.json is not touched.
"""

from __future__ import annotations

import argparse
import json
import time

from genomeos.attribution import crispri
from genomeos.attribution import target_calibration as tc
from genomeos.results import save_result


def show_census(name: str, census: dict[str, dict[str, object]]) -> None:
    print(f"  {name}")
    for s, c in census.items():
        print(
            f"    {s:24s} pairs {c['pairs']:>5}  regulated {c['regulated']:>4}  rate {c['rate']:<7} "
            f"on gate {c['on_gate_pairs']:>4}/{c['on_gate_rate']}  own term {c['own_term']}"
        )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", action="store_true", help="count the screens and stop")
    ap.add_argument("--chroms", default="", help="comma-separated chromosomes for the genome pass")
    args = ap.parse_args()
    t0 = time.time()

    training, heldout = crispri.load(crispri.TRAINING), crispri.load(crispri.HELDOUT)
    table = crispri.DeletionTable()
    if args.census:
        crispri.annotate(training, table)
        crispri.annotate(heldout, table)
        tc.add_features(training, table)
        tc.add_features(heldout, table)
        train = tc.scored(training)
        held = tc.scored([p for p in heldout if p.cell == tc.CELL])
        show_census("fitted (training K562)", tc.screen_census(train))
        show_census("read (held-out K562)", tc.screen_census(held))
        print(f"  in both: {sorted(set(tc.screen_census(train)) & set(tc.screen_census(held)))}")
        print(f"census only, nothing fitted ({time.time() - t0:.0f} s)")
        return 0

    chroms = tuple(c for c in args.chroms.split(",") if c) or tc.CHROMS
    result = tc.prevalence_terms(training, heldout, table, chroms=chroms)
    show_census("fitted (training K562)", result["census"]["fitted (training K562)"]["screens"])
    show_census("read (held-out K562)", result["census"]["read (held-out K562)"]["screens"])
    print(f"  in both: {result['census']['screens_in_both_populations']}")
    print(f"  on the gate: {json.dumps(result['census']['on_the_gate'], indent=1)}")
    b = result["held_out_before_and_after"]
    print(
        f"  held-out K562 shipped: {b['shipped curve']['bins_consistent']}/"
        f"{b['shipped curve']['bins']} bins, ECE {b['shipped curve']['ece']}"
    )
    print(
        f"  one pooled shift {b['one pooled shift (the published description)']['log_odds_shift']}: "
        f"{b['one pooled shift (the published description)']['bins_consistent']}/10 bins"
    )
    print(
        f"  per-screen offsets: {b['per-screen offsets']['bins_consistent']}/"
        f"{b['per-screen offsets']['bins']} bins, ECE {b['per-screen offsets']['ece']}"
    )
    loso = result["leave_one_screen_out"]
    for s, v in loso["per_screen"].items():
        print(
            f"    {s:24s} pooled {v['pooled_model']['bins_consistent']}/{v['pooled_model']['bins']} "
            f"ECE {v['pooled_model']['ece']:<7} -> offset "
            f"{v['offset_model']['bins_consistent']}/{v['offset_model']['bins']} "
            f"ECE {v['offset_model']['ece']:<7} residual {v['offset_model']['residual_log_odds_shift']}"
        )
    print(
        f"  screens improving {loso['screens_improving']}, weighted ECE {loso['weighted_ece']} "
        f"(fell {loso['weighted_ece_fell_by']}), median residual "
        f"{loso['median_absolute_residual_shift']}, verdict {loso['verdict']}"
    )
    for label, g in result["genome_under_each_screen"].items():
        print(f"  genome, {label}: {g['targets']} targets, top band 0.9-1 {g['top_band_0.9_1']}")
    path = save_result("target_prevalence", result)
    print(f"saved {path} ({time.time() - t0:.0f} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
