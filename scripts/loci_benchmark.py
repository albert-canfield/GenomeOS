# SPDX-License-Identifier: AGPL-3.0-or-later
"""Run the known-locus benchmark and commit its result (data/results/loci_benchmark.json).

    uv run python scripts/loci_benchmark.py            # full run: local layers, constraint, GTEx, Ensembl
    uv run python scripts/loci_benchmark.py --offline  # local layers only; constraint and eQTLs pending

Makes no AlphaGenome request. The network half reads Zoonomia phyloP and gnomAD Gnocchi by range,
streams GTEx v8's eQTL archive once (about 1.4 GB, kept only inside the windows, under
data/knowledge/loci_benchmark), and asks Ensembl for the named variants' positions and frequencies.
The tests (tests/test_loci_benchmark.py) read the committed result, so CI stays offline. Design and
numbers: docs/LOCI-BENCHMARK.md.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from genomeos.benchmark.loci import run_and_save  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--offline", action="store_true", help="local layers only")
    ap.add_argument(
        "--gate",
        action="store_true",
        help="only the twelve pinned loci, into loci_benchmark_gate.json: a quicker pass when the"
        " question is whether anything regressed rather than what the widened panel says",
    )
    args = ap.parse_args()
    t0 = time.time()

    def say(msg: str) -> None:
        print(f"[{time.time() - t0:7.1f}s] {msg}", flush=True)

    out = run_and_save(network=not args.offline, progress=say, gate_only=args.gate)
    a = out["aggregate"]
    n = out["negative_controls"]
    print()
    for key in ("target_derived", "target_derived_where_reachable", "target_heuristic", "target_looked_up"):
        v = a[key]
        print(f"{key:34s} {v['k']}/{v['n']}  hits {v['hits']}  misses {v['misses']}")
    print(f"unreachable target: {a['unreachable_target']}")
    for claim in ("target", "cell", "direction"):
        p, q = n["positives"][claim], n["negatives"][claim]
        print(f"claims a {claim:9s} positives {p['k']}/{p['n']}  negatives {q['k']}/{q['n']}")
    stem = "loci_benchmark_gate" if args.gate else "loci_benchmark"
    print(f"saved data/results/{stem}.json in {out['cost']['seconds']} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
