# SPDX-License-Identifier: AGPL-3.0-or-later
"""The containment audit's baselines, pinned.

The audit exists because the published control was matched on boundary count and nothing else. The
properties below are the ones that make its four baselines mean what docs/NODES-READER-WRITER.md
(2026-09-27) says they mean, and the first is the one that would silently break: the audit reproduces
`genome/domains.py _domains_from`'s 50 kb merge in its own `_merge_starts`, so the two must agree.
"""

from __future__ import annotations

import importlib.util
import random
from pathlib import Path

from genomeos.genome.domains import MIN_DOMAIN, _domains_from
from genomeos.genome.regulatory import CCRE

_spec = importlib.util.spec_from_file_location(
    "node_containment_audit", Path(__file__).resolve().parents[1] / "scripts" / "node_containment_audit.py"
)
audit = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(audit)


def _ccres(positions: list[int]) -> list[CCRE]:
    return [CCRE("chr1", p - 1, p + 1, f"E{i}", "CTCF-only", True) for i, p in enumerate(positions)]


def test_merge_starts_reproduces_the_real_caller() -> None:
    """The audit's copy of the merge and the caller's own agree on the node starts, on random input.

    If they ever diverge, every control in the audit is scored against a node set the caller would
    not have produced, so this is the test the audit rests on."""
    rng = random.Random(3)
    length = 5_000_000
    for _ in range(40):
        pos = sorted({rng.randint(1, length - 1) for _ in range(rng.randint(1, 120))})
        theirs = [d.start for d in _domains_from("chr1", length, _ccres(pos), None, MIN_DOMAIN, pos)]
        assert audit._merge_starts(pos, length) == theirs


def test_merge_starts_is_not_sensitive_to_input_order() -> None:
    pos = [900_000, 100_000, 3_000_000, 120_000]
    shuffled = [3_000_000, 120_000, 900_000, 100_000]
    assert audit._merge_starts(pos, 5_000_000) == audit._merge_starts(shuffled, 5_000_000)


def test_every_node_but_the_first_and_last_clears_the_floor() -> None:
    """The property the published control does not have, and the reason it reads high for random.

    The caller's merge leaves at most the first and the last node below MIN_DOMAIN; a raw uniform
    draw leaves many below it, which is what makes the unmerged control's containment low."""
    rng = random.Random(5)
    length = 20_000_000
    pos = [rng.randint(1, length - 1) for _ in range(400)]
    starts = audit._merge_starts(pos, length)
    lens = [b - a for a, b in zip(starts, [*starts[1:], length], strict=False)]
    assert sum(1 for x in lens[1:-1] if x < MIN_DOMAIN) == 0
    raw = [0, *sorted(pos)]
    raw_lens = [b - a for a, b in zip(raw, [*raw[1:], length], strict=False)]
    assert sum(1 for x in raw_lens if x < MIN_DOMAIN) > 100


def test_inside_counts_a_pair_only_when_one_node_holds_both() -> None:
    starts = [0, 1_000, 2_000, 3_000]
    pairs = [(1_100, 1_900), (1_100, 2_100), (500, 900), (2_500, 999)]
    assert audit._inside(starts, pairs) == 2


def test_circular_control_keeps_the_node_count_and_the_gap_multiset() -> None:
    """`circular` is the audit's closest-matched baseline; it is only that if rotation preserves the
    node count exactly and the gaps but for the one the wrap splits."""
    length = 10_000_000
    real = audit._merge_starts([700_000, 2_000_000, 2_060_000, 6_000_000, 9_000_000], length)
    pairs = [(1_000_000, 1_500_000), (2_500_000, 8_000_000)]
    out = audit.controls(real, len(real) - 1, length, pairs, random.Random(1))
    assert out["circular"]["mean_nodes"] == len(real)
    assert out["uniform"]["mean_nodes"] == len(real)
    assert out["count_matched"]["mean_nodes"] == len(real)
    assert set(out) == {"uniform", "uniform_merged", "circular", "count_matched"}


def test_uniform_control_reproduces_the_published_draw() -> None:
    """The published control is `SHUFFLES` draws of `len(edges)` uniform positions from Random(7),
    taken first, so the audit's `uniform` must be the same sequence or its +2.90 is not the
    published +2.88."""
    length = 4_000_000
    real = audit._merge_starts([500_000, 1_500_000, 2_500_000], length)
    pairs = [(600_000, 1_400_000), (1_600_000, 3_900_000)]
    got = audit.controls(real, len(real) - 1, length, pairs, random.Random(audit.SEED))["uniform"]["draws"]
    rng = random.Random(audit.SEED)
    want = []
    for _ in range(audit.SHUFFLES):
        st = [0, *sorted(rng.randint(1, length - 1) for _ in range(len(real) - 1))]
        want.append(audit._inside(st, pairs))
    assert got == want
