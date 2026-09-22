# SPDX-License-Identifier: AGPL-3.0-or-later
"""Does motif *arrangement* predict measured enhancer activity better than motif counts?

    uv run python scripts/motif_grammar.py scan [CHROM ...]   # strict sites per element, cached locally
    uv run python scripts/motif_grammar.py preregister         # penalty chosen inside the training set only
    uv run python scripts/motif_grammar.py score               # the held-out chromosomes, scored once

`scan` reads each chromosome's lentiMPRA elements (data/knowledge/mpra/rows_chr*.json), fetches their
200 bp of reference sequence and writes every JASPAR site at relative score 0.95, collapsed per TFClass
family unit, to data/knowledge/mpra/grammar_sites_chr*.json (local, tens of MB). `preregister` fixes the
claim and the ridge penalty using the training chromosomes alone; `score` fits on the training
chromosomes, scores the held-out ones and saves `motif_grammar.json`. Nothing here uses the network
beyond the one-off JASPAR download, and no model quota.
"""

from __future__ import annotations

import sys
import time

from genomeos.attribution.motif_grammar import (
    CHROMS,
    HELD_OUT,
    PREREGISTERED,
    run_and_save,
    scan_chromosome,
)


def scan(chroms: list[str]) -> int:
    for c in chroms or list(CHROMS):
        t = time.time()
        n = scan_chromosome(c)
        print(f"{c}: {n:,} elements scanned in {time.time() - t:.1f}s", flush=True)
    return 0


def report(score_held_out: bool) -> int:
    print(f"pre-registered: {PREREGISTERED}", flush=True)
    print(f"held out: {', '.join(HELD_OUT)}", flush=True)
    r = run_and_save(score_held_out=score_held_out, progress=lambda m: print(f"  {m}", flush=True))
    held = r.get("held_out")
    if not held:
        print(f"training only: {r['training_elements']:,} elements, {r['columns']['c']} columns", flush=True)
        return 0
    for cell, rec in held.items():
        print(
            f"{cell} ({rec['held_out_elements']:,} held-out elements): "
            f"(a) {rec['a']['point']:.3f} (b) {rec['b']['point']:.3f} (c) {rec['c']['point']:.3f} "
            f"(c shuffled) {rec['shuffled']['point']:.3f}; "
            f"c-b {rec['c-b']['point']:+.4f} {rec['c-b']['ci95']}; "
            f"c-shuffled {rec['c-shuffled']['point']:+.4f} {rec['c-shuffled']['ci95']}",
            flush=True,
        )
    print(f"verdict: {r['verdict']['statement']}", flush=True)
    return 0


def main(argv: list[str]) -> int:
    stage = argv[0] if argv else "score"
    if stage == "scan":
        return scan(argv[1:])
    if stage == "preregister":
        return report(False)
    if stage == "score":
        return report(True)
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
