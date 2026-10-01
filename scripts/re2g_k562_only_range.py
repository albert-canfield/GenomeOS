# SPDX-License-Identifier: AGPL-3.0-or-later
"""Reproduce `re2g.K562_ONLY_RANGE`: what a K562-only prediction file can attain as a pooled weighted
AUPRC over all 4,378 held-out pairs.

    uv run --frozen python scripts/re2g_k562_only_range.py

The question this answers is whether the one prediction file that fits the lane's download budget could
serve the gate. It cannot, and not for the reason first written down: the range it can attain brackets
the published figure, so a figure near the target would be no evidence that the scores are the published
ones.

Real labels and real weights, synthetic scores. No comparator score is read, and nothing here is a
measurement of either model: the two ends are a no-skill ranking and a perfectly separating ranking
inside K562, with every non-K562 pair forced to 0 by the benchmark's own fill rule. Writes no result.
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos.attribution import crispri, re2g  # noqa: E402

SEED = 0


def main() -> None:
    pairs = crispri.load(crispri.HELDOUT)
    labels = [p.regulated for p in pairs]
    weights = [p.weight for p in pairs]
    rng = random.Random(SEED)

    def pooled(scores: list[float]) -> float:
        value = crispri.benchmark_auprc(scores, labels, weights)
        assert value is not None
        return round(value, 4)

    # Perfect separation inside K562: positives in the top half of the range, negatives in the bottom,
    # continuous so the estimator sees one point per pair rather than a handful of coarse ties.
    perfect = [
        0.0 if p.cell != "K562" else (0.5 + rng.random() * 0.5 if p.regulated else rng.random() * 0.5)
        for p in pairs
    ]
    # No skill inside K562: the floor any K562-only file sits above.
    no_skill = [0.0 if p.cell != "K562" else rng.random() for p in pairs]

    low, high = pooled(no_skill), pooled(perfect)
    k562 = [p for p in pairs if p.cell == "K562"]
    weight_in_k562 = sum(p.weight for p in k562 if p.regulated)
    weight_all = sum(p.weight for p in pairs if p.regulated)

    print("population: all 4,378 held-out pairs, pooled over the five cell types")
    print(
        f"  weighted positives in K562: {weight_in_k562:.2f} of {weight_all:.2f} "
        f"({100 * weight_in_k562 / weight_all:.1f}%)"
    )
    print(f"  pairs forced to 0 by a K562-only file: {len(pairs) - len(k562)}")
    print(f"  attainable pooled weighted AUPRC: {low} (no skill in K562) to {high} (perfect in K562)")
    print(f"  published target: {re2g.GATE_TARGET:.4f}, tolerance {re2g.GATE_TOLERANCE}")
    brackets = low <= re2g.GATE_TARGET <= high
    print(f"  the attainable range brackets the target: {brackets}")
    print(f"  registered K562_ONLY_RANGE: {re2g.K562_ONLY_RANGE}")
    if (low, high) != re2g.K562_ONLY_RANGE:
        raise SystemExit(f"re2g.K562_ONLY_RANGE is {re2g.K562_ONLY_RANGE}, this run gives {(low, high)}")


if __name__ == "__main__":
    main()
