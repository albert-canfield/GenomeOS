# SPDX-License-Identifier: AGPL-3.0-or-later
"""The node model against published rearrangements, free half (data/results/loci_rearrangements.json).

    uv run python scripts/rearrangements.py

Makes no model request and needs no network: the ENCODE cCRE registry, GENCODE and the project's own
node model, read locally. Design, provenance and the registered prediction: docs/LOCI-BENCHMARK.md
section 11b.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from genomeos.benchmark.rearrangements import run_and_save  # noqa: E402


def main() -> int:
    t0 = time.time()

    def say(msg: str) -> None:
        print(f"[{time.time() - t0:6.1f}s] {msg}", flush=True)

    out = run_and_save(progress=say)
    print()
    print("registered before the run:")
    print(f"  {out['prediction_registered_before_the_run']}")
    print()
    print("no published case is runnable:")
    print(f"  {out['no_case_is_runnable']}")
    print()
    b = out["boundary_behaviour"]
    lo, hi = b["published_deletion_size"]
    print(f"a deletion of the published size ({lo / 1e6:.2f}-{hi / 1e6:.2f} Mb), anywhere on {b['chrom']}:")
    for claim, r in b["rates"].items():
        print(f"  {claim:18s} {r['k']}/{r['n']}")
    near, total = b["scored_near_the_recipient"], len(b["deletions"])
    print(f"  recipients with any deletion scored near them: {near}/{total}")
    print()
    print(out["analytic_note"])
    print()
    print(f"saved data/results/loci_rearrangements.json in {out['cost']['seconds']} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
