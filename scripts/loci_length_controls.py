#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Score the length-matched control arm against the known-locus panel. No model requests.

    python3 scripts/loci_length_controls.py            # the full run, one range read per chromosome
    python3 scripts/loci_length_controls.py --offline   # arms and comparison from what is committed

The registration is `genomeos.benchmark.loci_length.PREREGISTRATION`, printed first so a run cannot
be read without it.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos.benchmark.loci_length import PREREGISTRATION, run_and_save  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--offline", action="store_true", help="no range reads; arm C keeps no syntax reading")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    print(PREREGISTRATION)
    print("-" * 100)
    say = (lambda _m: None) if args.quiet else (lambda m: print(m, flush=True))
    out = run_and_save(network=not args.offline, progress=say)

    for name, arm in out["arms"].items():
        print(f"\n{name}: {arm['what']}")
        print(
            f"  n {arm['n']}  storage {arm['storage']['k']}/{arm['storage']['n']} = {arm['storage']['rate']}"
        )
        print(f"  length {arm['length']}")
        tss = arm["distance_to_coding_tss"].get("median")
        print(f"  gc median {arm['gc'].get('median')}  tss distance median {tss}")
        print(f"  constrained fraction median {arm['constrained_fraction'].get('median')}")
        print(f"  coverage {arm['coverage']}")
        print(f"  deletion scored {arm['deletion_scored']}  never looked {arm['no_element_scored']}")
    print("\nlength matching:", json.dumps(out["length_matching"], indent=1))
    decides = out["verdict"]["decides_on"]
    print(
        f"\nwhich stratification decides: arm B on {decides['B_existing']},"
        f" arm C on {decides['C_length_only']} (arm C is matched on nothing but length, so its"
        " covariates have to be standardised away before the verdict reads it)"
    )
    for arm, comps in out["comparisons"].items():
        print(f"\npanel against {arm}")
        for name, c in comps.items():
            m, raw = c["matched"], c["raw"]
            mark = "*" if name == decides[arm] else " "
            print(
                f" {mark}{name:24s} raw {raw['difference']}  matched {m['difference']}"
                f" (p {m['p_one_sided']}, upper95 {m['upper_95']})"
                f"  matched n {c['targets_matched']}  dropped {c['dropped_for_want_of_a_control']}"
            )
    print("\nverdict:", out["verdict"]["verdict"])
    for r in out["verdict"]["reasons"]:
        print("  ", r)
    print("  the branches are not equally hard:", out["verdict"]["the_branches_are_not_equally_hard"])
    print("model requests spent:", out["model_requests"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
