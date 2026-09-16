# SPDX-License-Identifier: AGPL-3.0-or-later
"""Hold the enhancer-to-gene attribution against measured CRISPRi perturbations.

    uv run python scripts/crispri_score.py

Streams the ENCODE enhancer-gene benchmark tables into data/knowledge/crispri once, joins them to
the all-enhancer deletion table and saves data/results/crispri_benchmark.json. No model request.
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
    result = crispri.score(training, heldout, crispri.DeletionTable())
    path = save_result("crispri_benchmark", result)
    loco = result["training_leave_chromosome_out"]
    print(f"coverage: {result['coverage']}")
    for name, m in loco["models"].items():
        print(f"  K562 leave-chromosome-out  {name:32s} AUPRC {m['auprc']}  AUROC {m['auroc']}")
    print(f"  deletion gain {loco['deletion_gain']}")
    for cell, h in result["heldout"].items():
        if "refused" in h:
            print(f"  held out {cell}: refused ({h['refused']})")
            continue
        for name, m in h["models"].items():
            counts = f"{m['positives']}/{m['pairs']}"
            print(f"  held out {cell:8s} {name:32s} AUPRC {m['auprc']}  AUROC {m['auroc']}  ({counts})")
        print(f"  held out {cell}: deletion gain {h['deletion_gain']}, passes {h['passes']}")
    for cell, h in result["heldout_elements_not_in_training"].items():
        m = h["models"]["activity + distance + deletion"]
        print(f"  held out {cell}, elements not in training: AUPRC {m['auprc']}, gain {h['deletion_gain']}")
    print(f"pre-registered: {result['preregistered']}")
    print(f"verdict: {result['verdict']}  ({time.time() - t0:.0f} s) -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
