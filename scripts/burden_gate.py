# SPDX-License-Identifier: AGPL-3.0-or-later
"""Run stage 2's gate (a) against its registration (runtime/burden_gate.py). No network, no quota.

    uv run python scripts/burden_gate.py

The design is `PRE_REGISTRATION`, committed in `genomeos/runtime/burden_gate.py` before this file
existed. Nothing here may change it: this expresses the costly protein at the registered levels,
holds the unrelated gene's own demand fixed, reads the factor each one is scaled by, and applies the
registered bars without interpretation.

The falsifier arm re-runs the same demands against pools multiplied by `UNCONSTRAINED_MULTIPLIER`.
It is the arm a bug passes and the reason the gate can fail: a model that cuts an unrelated gene
whatever the capacity is not reproducing a burden, and only this arm can tell the two apart.
"""

from __future__ import annotations

import argparse
import time
from dataclasses import replace
from itertools import pairwise
from pathlib import Path
from typing import Any

from genomeos.lang import parse_file
from genomeos.results import save_result
from genomeos.runtime.burden_gate import (
    LEVELS,
    PRE_REGISTRATION,
    UNCONSTRAINED_MULTIPLIER,
    UNRELATED_DEMAND,
)
from genomeos.runtime.economy import Economy

PROGRAM = Path("data/demo/stage2_pools.bio")
COSTLY, UNRELATED = "BURDENp", "HOUSEKEEPERp"
#: the lengths the per-residue and per-nucleotide costs multiply. Stated here rather than guessed
#: inside the economy, which refuses a per-length cost with no length rather than charging zero.
LENGTHS = {"BURDEN": 900.0, "HOUSEKEEPER": 900.0, "BURDENp": 300.0, "HOUSEKEEPERp": 300.0}


def curve(economy: Economy) -> list[dict[str, Any]]:
    """The unrelated gene's factor at each registered level, its own demand never changing."""
    rows = []
    for level in LEVELS:
        report = economy.report({COSTLY: level, UNRELATED: UNRELATED_DEMAND})
        rows.append(
            {
                "level": level,
                "unrelated_factor": report["factor"][UNRELATED],
                "costly_factor": report["factor"][COSTLY],
                "oversubscribed": sorted(p for p, v in report["pools"].items() if v["oversubscribed"]),
            }
        )
    return rows


def unconstrained(economy: Economy) -> Economy:
    """The same economy with every pool large enough that no registered level can exhaust it."""
    pools = {
        pid: replace(p, size=float(p.size) * UNCONSTRAINED_MULTIPLIER) for pid, p in economy.pools.items()
    }
    return Economy(pools=pools, costs=economy.costs, allocations=economy.allocations, lengths=economy.lengths)


def verdict(constrained: list[dict], free: list[dict]) -> dict[str, Any]:
    """The registered bars, applied without interpretation."""
    factors = [r["unrelated_factor"] for r in constrained]
    direction = factors[-1] < factors[0]
    # deliberately offset by one, so the lists differ in length: strict= would be wrong here, and
    # itertools.pairwise says the intent rather than leaving a reader to check the slice
    monotonic = all(b <= a + 1e-12 for a, b in pairwise(factors))
    falsifier = all(r["unrelated_factor"] == 1.0 for r in free)
    passed = direction and monotonic and falsifier
    return {
        "direction": direction,
        "direction_reads": f"{factors[0]} at the lowest level, {factors[-1]} at the highest",
        "monotonic": monotonic,
        "falsifier_holds": falsifier,
        "falsifier_reads": (
            "every factor stays 1.0 when no pool is the constraint"
            if falsifier
            else "a reduction appeared with the pool switched off, which §5.3 calls falsified"
        ),
        "magnitude": "not askable: no measured burden curve is held by this project (registered)",
        "passed": passed,
        "verdict": (
            "gate (a) passes on the three askable clauses; the magnitude clause is unasked and "
            "contributes nothing to this"
            if passed
            else "gate (a) FAILS: "
            + ", ".join(
                n
                for n, ok in (
                    ("direction", direction),
                    ("monotonicity", monotonic),
                    ("falsifier", falsifier),
                )
                if not ok
            )
        ),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args(argv)

    t0 = time.time()
    module = parse_file(str(PROGRAM))
    economy = Economy.from_module(module, LENGTHS)
    constrained = curve(economy)
    free = curve(unconstrained(economy))
    got = verdict(constrained, free)

    print(f"{'level':>10} {'unrelated':>10} {'costly':>8}  oversubscribed")
    for r in constrained:
        print(
            f"{r['level']:10.0e} {r['unrelated_factor']:10.4f} {r['costly_factor']:8.4f}  "
            f"{','.join(r['oversubscribed']) or '-'}"
        )
    print(f"\nfalsifier arm (pools x {UNCONSTRAINED_MULTIPLIER:.0e}):")
    print(f"  unrelated factors: {sorted({r['unrelated_factor'] for r in free})}")
    print(f"\n{got['verdict']}")

    out = {
        "result": "burden_gate",
        "pre_registration": PRE_REGISTRATION,
        "program": str(PROGRAM),
        "lengths": LENGTHS,
        "constrained_arm": constrained,
        "unconstrained_arm": free,
        **got,
        "seconds": round(time.time() - t0, 2),
    }
    if args.no_save:
        print("\nnot saved (--no-save)")
    else:
        print(f"\nsaved {save_result(out['result'], out)}")
    return 0 if got["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
