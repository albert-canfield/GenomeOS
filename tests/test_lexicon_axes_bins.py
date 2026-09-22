# SPDX-License-Identifier: AGPL-3.0-or-later
"""The byte-quantised phyloP track: its threshold, its anchor and what each bin edge means.

A duplication sweep flagged `PHYLOP_ANCHOR = 23` as a hardcoded threshold that a change to
`constraint.PHYLOP_THRESHOLD` would silently bypass. Checked, that is half right. The anchor is not a
threshold: `quantise_phylop` subtracts `PHYLOP_THRESHOLD` before adding the anchor, so
`q >= PHYLOP_ANCHOR` means `p >= PHYLOP_THRESHOLD` for any threshold. What the sweep was right about
is that this module redeclares the threshold instead of importing it, and that `PHYLOP_BINS` is a
tuple of fixed `q` values whose float meanings move when the threshold does — and three of the float
labels in its comment were already wrong by up to 0.33.

These tests pin the derivation rather than the labels, so the next person to move the threshold gets
a failure instead of a stale comment.
"""

from __future__ import annotations

import numpy as np

from genomeos.attribution import constraint
from genomeos.attribution import lexicon_axes as lx


def test_the_two_declarations_of_the_threshold_agree() -> None:
    """Two homes for one number is the divergence; equal today, and asserted so it stays equal."""
    assert lx.PHYLOP_THRESHOLD == constraint.PHYLOP_THRESHOLD


def test_the_anchor_means_the_threshold_whatever_the_threshold_is() -> None:
    """`q >= 23` is exactly `p >= 2.27` by construction, not by coincidence."""
    values = np.array([2.2699, 2.27, 2.2701, 0.0, 5.0, -1.3], dtype=np.float32)

    q = lx.quantise_phylop(values)

    assert list(q >= lx.PHYLOP_ANCHOR) == list(values >= lx.PHYLOP_THRESHOLD)


def test_missing_stays_missing_through_a_round_trip() -> None:
    q = lx.quantise_phylop(np.array([np.nan, 1.0], dtype=np.float32))

    assert q[0] == lx.PHYLOP_MISSING
    assert np.isnan(lx.dequantise_phylop(q)[0])
    assert not np.isnan(lx.dequantise_phylop(q)[1])


def test_every_bin_edge_means_the_phylop_its_comment_claims() -> None:
    """The edges are q values; their float meanings are derived, and the comment now states those.

    Recomputed here rather than copied: the old comment said -1.1, 1.3 and 4 where the arithmetic
    gives -0.93, +0.97 and +3.97.
    """
    edges = lx.PHYLOP_BINS[1:-1]
    meant = [(e - lx.PHYLOP_ANCHOR) / lx.PHYLOP_SCALE + lx.PHYLOP_THRESHOLD for e in edges]

    assert [round(m, 2) for m in meant] == [-0.93, -0.03, 0.97, 2.27, 3.97]
    assert lx.PHYLOP_ANCHOR in edges, "the constrained threshold must itself be a bin edge"


def test_the_bins_are_sorted_and_span_the_byte_range() -> None:
    assert list(lx.PHYLOP_BINS) == sorted(lx.PHYLOP_BINS)
    assert lx.PHYLOP_BINS[0] == lx.PHYLOP_MISSING
