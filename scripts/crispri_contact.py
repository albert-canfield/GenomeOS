# SPDX-License-Identifier: AGPL-3.0-or-later
"""Measured 3D contact instead of 1/distance in the CRISPRi enhancer-to-gene comparison.

    uv run python scripts/crispri_contact.py

Reads the contact between each CRISPRi element and its gene's TSS from a released 4D Nucleome
in-situ Hi-C matrix (K562 and GM12878, GRCh38, 5 kb bins) by HTTP range request: the matrices are
6.7 GB and 22.6 GB and neither is downloaded. The answers are cached under data/knowledge/hic_contact
so a second run reads nothing from the network. Saves data/results/crispri_contact.json.

Needs FOURDN_KEY and FOURDN_SECRET in .env. No model request.
"""

from __future__ import annotations

import time

from genomeos.attribution import crispri
from genomeos.genome import hic_contact
from genomeos.results import save_result


def main() -> int:
    t0 = time.time()
    for name in (crispri.TRAINING, crispri.HELDOUT):
        crispri.fetch(name)
    training, heldout = crispri.load(crispri.TRAINING), crispri.load(crispri.HELDOUT)
    source = hic_contact.ContactSource()
    try:
        result = crispri.score_contact(
            training, heldout, crispri.DeletionTable(), source, progress=lambda m: print(f"  {m}")
        )
    finally:
        for p in source.save():
            print(f"  cache: {p}")
    path = save_result("crispri_contact", result)
    print(f"coverage: {result['coverage']}")
    loco = result["training_leave_chromosome_out"]
    for name, m in loco["models"].items():
        print(f"  K562 leave-chromosome-out  {name:32s} AUPRC {m['auprc']}  AUROC {m['auroc']}")
    for key in ("contact_gain", "contact_gain_with_deletion", "contact_added_to_distance"):
        print(f"  training {key}: {loco[key]}")
    for cell, h in result["heldout"].items():
        if "refused" in h:
            print(f"  held out {cell}: refused ({h['refused']})")
            continue
        for name, m in h["models"].items():
            counts = f"{m['positives']}/{m['pairs']}"
            print(f"  held out {cell:8s} {name:32s} AUPRC {m['auprc']}  AUROC {m['auroc']}  ({counts})")
        for key in ("contact_gain", "contact_gain_with_deletion", "contact_added_to_distance"):
            print(f"  held out {cell}: {key} {h[key]}")
        print(f"  held out {cell}: passes {h['passes']}")
    print(f"pre-registered: {result['preregistered']}")
    print(f"verdict: {result['verdict']}  ({time.time() - t0:.0f} s) -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
