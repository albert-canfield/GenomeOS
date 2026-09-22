# SPDX-License-Identifier: AGPL-3.0-or-later
"""Gene-level phylogenetic profiling, and the two nulls that disagree about it."""

from __future__ import annotations

import importlib.util
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("profiling_genes", ROOT / "scripts" / "profiling_genes.py")
assert spec and spec.loader
pg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pg)


def test_jaccard_is_the_overlap_over_the_union() -> None:
    assert pg.jaccard(0b1100, 0b0110) == 1 / 3
    assert pg.jaccard(0b1111, 0b1111) == 1.0
    assert pg.jaccard(0b1000, 0b0001) == 0.0
    assert pg.jaccard(0, 0) == 0.0  # no species either side is not a perfect match


def test_the_shuffle_keeps_prevalence_and_loses_the_species() -> None:
    """A null free to change how many species a gene is in would call two ubiquitous genes a finding."""
    rng = random.Random(3)
    mask = 0b1010101010

    out = pg.shuffled(mask, 10, rng)

    assert out.bit_count() == mask.bit_count()
    assert out < (1 << 10)


def test_a_matched_draw_uses_real_genes_of_similar_prevalence() -> None:
    """The null that matters is built out of the same kind of thing being tested: real profiles."""
    rng = random.Random(5)
    members = [0b111111, 0b111111]  # prevalence 6
    pool = [0b111111, 0b111110, 0b000011, 0b000001]  # two near, two far

    drawn = pg.matched_draw(members, pool, rng)

    assert len(drawn) == 2
    assert all(abs(d.bit_count() - 6) <= 1 for d in drawn)


def test_a_matched_draw_falls_back_rather_than_returning_nothing() -> None:
    """A member with no near-prevalence gene in the pool still gets a control, and it is a real gene."""
    rng = random.Random(7)
    drawn = pg.matched_draw([0b1111111111], [0b1], rng)

    assert drawn == [0b1]


def test_a_library_below_its_matched_null_is_reported_as_such() -> None:
    """Fourteen of 42 real libraries are at or below the matched-gene null; the shape must survive."""
    # ten genes sharing one set of 30 species out of 40: the shuffle can move them, real genes with
    # the same set cannot be moved apart, which is exactly the contrast section 19 turns on
    shared = (1 << 30) - 1
    by_gene = {f"g{i}": shared for i in range(10)}
    pool = [shared for _ in range(50)]

    row = pg.profile_library("shared_set", sorted(by_gene), by_gene, 40, pool)

    assert row is not None
    assert row["mean_pair_jaccard"] == 1.0
    # against the prevalence shuffle they look perfect; against real genes of the same prevalence they
    # are ordinary, which is the whole finding of section 19
    assert row["excess_over_null"] > 0
    assert row["matched_draws_at_or_above_observed"] == row["draws"]


def test_a_library_below_the_member_floor_is_skipped() -> None:
    by_gene = {"a": 0b111, "b": 0b101}
    assert pg.profile_library("tiny", ["a", "b"], by_gene, 3, [0b111]) is None
