# SPDX-License-Identifier: AGPL-3.0-or-later
"""The upward arm of the direction question: the half `crispri_direction.py` could not test.

    uv run python scripts/crispri_direction_both.py

Reads the per-element response cache of the finished genome-wide sweep, which was written with
threshold=0.0 and therefore carries every gene in the scorer's window rather than the one target the
compact table keeps. That makes the upward measured pairs answerable and costs nothing: the 30 to 35
AlphaGenome requests this was costed at would re-buy an answer that is already on disk.

Saves data/results/crispri_direction_both.json. No model request is made.
"""

from __future__ import annotations

import time

from genomeos.attribution import crispri, crispri_direction_both
from genomeos.results import save_result


def main() -> int:
    t0 = time.time()
    for name in (crispri.TRAINING, crispri.HELDOUT):
        crispri.fetch(name)
    result = crispri_direction_both.direction_both()
    assert result["requests"] == 0, "this measurement spends no request"
    path = save_result("crispri_direction_both", result)

    cost = result["costing_of_the_route_not_taken"]
    print("costing of the route not taken:")
    for label in ("K562", "K562+GM12878"):
        c = cost[label]
        print(f"  {label:14s} {c['pairs']} upward pairs on {c['requests']} distinct elements")
    print(f"  spent: {cost['spent']}")

    red = result["rederivation_of_the_44"]
    print(
        f"re-derivation of the 44: passes={red.get('passes')} "
        f"missing={red.get('missing_here')} sign_changed={red.get('sign_changed')}"
    )

    up, down = result["upward_arm"], result["downward_arm"]
    print(f"UPWARD arm   {up['k']}/{up['n']} = {up['rate']}  95% {up['ci95']}")
    print(f"  marginal up-rate (the constant-sign caller): {up['marginal_up_rate']}")
    print(f"  verdict: {up['verdict']}")
    print(f"    {up['why']}")
    print(f"DOWNWARD arm {down['k']}/{down['n']} = {down['rate']}  95% {down['ci95']}")

    comb = result["combined"]
    ba, raw = comb["balanced_accuracy"], comb["raw_agreement"]
    print(
        f"COMBINED balanced accuracy {ba['rate']}  95% {ba['ci95']}  ({ba['n_down']} down, {ba['n_up']} up)"
    )
    print(
        f"  raw agreement {raw['k']}/{raw['n']} = {raw['rate']} "
        f"against a majority-class {raw['majority_class_rate']}"
    )
    print(f"  verdict: {comb['verdict']}  (band {comb['band']})")
    print(f"    {comb['why']}")

    mm = result["magnitude_matched_sub_test"]
    print(f"magnitude-matched, |predicted| > {mm['threshold']}:")
    print(
        f"  down {mm['down']['k']}/{mm['down']['n']} = {mm['down']['rate']}   "
        f"up {mm['up']['k']}/{mm['up']['n']} = {mm['up']['rate']}   "
        f"balanced {mm['balanced']['rate']}"
    )
    for sign in ("down", "up"):
        lad = result["magnitude_ladder"][sign]
        bins = "  ".join(f"{lo}:{s['agree']}/{s['pairs']}" for lo, s in lad["bins"].items())
        weak, strong = lad["at_or_below_0.05"], lad["above_0.1"]
        print(f"  ladder {sign}: {bins}")
        print(
            f"    <= 0.05 {weak['k']}/{weak['n']}   > 0.1 {strong['k']}/{strong['n']}   "
            f"largest error at {lad['largest_error_magnitude']}"
        )

    sel = result["top_target_selection"]
    a, b = sel["was_top_target"], sel["newly_answerable"]
    print(f"top-target selected: {a['k']}/{a['n']} = {a['rate']}")
    print(f"newly answerable:    {b['k']}/{b['n']} = {b['rate']}")

    for cell, m in result["per_cell"].items():
        print(
            f"  {cell:8s} up {m['up']['k']}/{m['up']['n']}  "
            f"down {m['down']['k']}/{m['down']['n']}  balanced {m['balanced']['rate']}"
        )

    tr = result["training_arm_fitted_on_not_a_test"]
    print(
        f"training arm (fitted on, not a test): up {tr['up']['k']}/{tr['up']['n']}  "
        f"down {tr['down']['k']}/{tr['down']['n']}  balanced {tr['balanced']['rate']}"
    )
    print(f"coverage strata: {result['coverage']['heldout']}")
    print(f"unanswered pairs: {result['coverage']['unanswered_pairs']}")
    print(
        "does the pass of 5c842ca survive restatement over both signs: "
        f"{result['pass_of_5c842ca_survives_restatement_over_both_signs']}"
    )
    print(f"  {result['why_survives']}")
    print(f"requests: {result['requests']}  ({time.time() - t0:.0f} s) -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
