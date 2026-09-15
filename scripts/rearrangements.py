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
    a = out["aggregate"]
    print()
    print("registered before the run:")
    print(f"  {out['prediction_registered_before_the_run']}")
    print()
    for claim in ("separated_before", "boundary_lost", "new_adjacency"):
        p, q = a["published"][claim], a["controls"][claim]
        print(f"{claim:18s} published {p['k']}/{p['n']}   matched random {q['k']}/{q['n']}")
    n = a["recipient_named_where_a_deletion_was_scored"]
    if n["judgeable"]:
        print(
            f"{'recipient named':18s} published {n['published']['k']}/{n['published']['n']}"
            f"   matched random {n['controls']['k']}/{n['controls']['n']}"
        )
    else:
        print(
            f"{'recipient named':18s} NOT JUDGEABLE: {n['controls']['not_scored']} of"
            f" {n['controls']['of']} control recipients have no deletion scored near them,"
            " so the claim would measure coverage and not biology"
        )
    for c in out["cases"]:
        b = c.get("boundary_census") or {}
        walls = b.get("ctcf_only_between_donor_and_recipient")
        print(f"  {c['locus']}: CTCF-only between donor and recipient = {walls}")
    print(f"saved data/results/loci_rearrangements.json in {out['cost']['seconds']} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
