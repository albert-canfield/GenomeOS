# SPDX-License-Identifier: AGPL-3.0-or-later
"""The CRISPRi result on the published pair sets, beside the published figures, and the coverage arms.

    uv run python scripts/crispri_published.py

Scores the frozen model (crispri.FEATURES, weights fitted on the covered K562 training pairs) by the
ENCODE benchmark's own estimator on all 10,356 training pairs (hold-one-chromosome-out) and all
4,378 held-out pairs (pooled, weighted by direct-effect probability), sets the result beside
Gschwind et al. 2026's Supplementary Table 3, and runs the three coverage arms on held-out K562.
The registration is crispri.PREREGISTERED_PUBLISHED, committed before this script existed.
No AlphaGenome request: the deletion values are the sweep's own cache.
"""

from __future__ import annotations

import time

from genomeos.attribution import crispri
from genomeos.results import save_result


def main() -> int:
    t0 = time.time()
    for name in (crispri.TRAINING, crispri.HELDOUT):
        crispri.fetch(name)
    training, heldout = crispri.load(crispri.TRAINING), crispri.load(crispri.HELDOUT)
    result = crispri.score_published(training, heldout, crispri.DeletionTable(), crispri.ElementCache())
    path = save_result("crispri_published", result)
    print(f"estimator check: {result['estimator_check_raw_distance']}")
    tr = result["training_published_split"]
    for name, m in tr["models"].items():
        print(f"  training LOCO, all pairs  {name:32s} AUPRC {m['auprc']}")
    print(f"  training gain {tr['deletion_gain']}; vs ENCODE-rE2G: {tr['against_encode_re2g']}")
    ho = result["heldout_published_pairs"]
    for name, m in ho["models"].items():
        print(f"  held out pooled weighted  {name:32s} AUPRC {m['auprc']}")
    print(f"  held-out gain {ho['deletion_gain']}; vs ENCODE-rE2G: {ho['against_encode_re2g']}")
    print(f"  baseline vs ABC: {ho['baseline_against_abc']}")
    for cell, c in ho["per_cell_type_weighted"].items():
        m = c["models"]
        print(
            f"  {cell:8s} weighted: a+d {m['activity + distance']['auprc']}  "
            f"+del {m['activity + distance + deletion']['auprc']}  gain {c['deletion_gain']}"
        )
    arms = result["coverage_arms_k562_heldout"]
    for k in ("arm1_all_pairs", "arm2_coverage_matched", "arm3_coverage_indicator_control"):
        print(f"  {k}: {arms[k]}")
    print(f"  invariant to coverage: {arms['invariant_to_coverage']}")
    print(f"second cell type: {result['second_cell_type']['cost']}")
    print(f"({time.time() - t0:.0f} s) -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
