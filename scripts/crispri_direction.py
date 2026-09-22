# SPDX-License-Identifier: AGPL-3.0-or-later
"""Does the deletion layer get the *direction* right, held against measured CRISPRi perturbations?

    uv run python scripts/crispri_direction.py

Reads the signed EffectSize of the ENCODE enhancer-gene benchmark against the signed predicted log2
fold change of the finished genome-wide deletion sweep, and saves data/results/crispri_direction.json.
No model request is made: both sides are already on disk.

The headline is judged against the model's own marginal down-rate over the whole sweep, not against
0.5. See genomeos/attribution/crispri_direction.py for why 0.5 would manufacture a success.
"""

from __future__ import annotations

import time

from genomeos.attribution import crispri, crispri_direction
from genomeos.results import save_result


def main() -> int:
    t0 = time.time()
    for name in (crispri.TRAINING, crispri.HELDOUT):
        crispri.fetch(name)
    result = crispri_direction.direction()
    path = save_result("crispri_direction", result)

    head = result["headline"]
    print(
        f"answerable held-out pairs: {head['n']} (predicted zero excluded: {head['predicted_zero_excluded']})"
    )
    print(f"  measured  down {head['measured_down']}  up {head['measured_up']}")
    print(f"  predicted down {head['predicted_down']}  up {head['predicted_up']}")
    print(f"  sign agreement {head['rate']}  95% {head['ci95']}")
    print(f"  both measured signs present: {head['both_measured_signs_present']}")
    for cell, m in result["per_cell"].items():
        print(f"    {cell:8s} {m['k']}/{m['n']} = {m['rate']}  95% {m['ci95']}")
    print(f"marginal down-rate (the constant-sign caller): {result['pooled_marginal_down_rate']}")
    for cell in result["cells"]:
        m = result["marginal_down_rate"].get(cell, {})
        print(f"    {cell:8s} {m.get('k')}/{m.get('n')} = {m.get('rate')}  95% {m.get('ci95')}")
    ph = result["posthoc_magnitude_matched"]
    print(f"post-hoc magnitude-matched marginal: {ph['marginal']['pooled']}")
    for cell in result["cells"]:
        print(f"    {cell:8s} {ph['marginal'][cell]['rate']}  over {ph['marginal'][cell]['pairs']} pairs")
    print(f"  verdict against it: {ph['verdict_against_it']['verdict']}")
    print("  agreement by |predicted log2fc|:")
    for lo, s in sorted(ph["agreement_by_magnitude"].items(), key=lambda t: float(t[0])):
        print(f"    >= {lo:<6} {s['agree']}/{s['pairs']}")
    train = result["training_arm_fitted_on_not_a_test"]["agreement"]
    print(f"training arm (fitted on, not a test): {train['k']}/{train['n']} = {train['rate']}")
    print(f"covered but the gene is not the top target: {result['what_it_named_instead']['pairs']} pairs")
    print(f"coverage strata: {result['coverage']['heldout']}")
    print(f"pre-registered: {result['preregistered']}")
    print(f"verdict: {result['verdict']}")
    print(f"  {result['why']}")
    print(f"requests: {result['requests']}  ({time.time() - t0:.0f} s) -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
