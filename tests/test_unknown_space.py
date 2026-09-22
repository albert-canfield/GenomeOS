# SPDX-License-Identifier: AGPL-3.0-or-later
"""The arithmetic behind the unknown-space readings: standardisation, coverage and the library."""

from __future__ import annotations

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


matched = _load("syntax_blocks_matched")
coverage = _load("unknown_coverage")
library = _load("unknown_library")


def _row(arm: str, moved: bool, length: int = 300, gc: float = 0.5, tss: int = 10_000) -> dict:
    return {"arm": arm, "moved": moved, "length": length, "gc": gc, "tss": tss}


# --- standardisation: the defect that was published on 2026-09-17 ----------------------------------


def test_a_large_stratum_cannot_outvote_a_small_one() -> None:
    """The bug: pooling every control once per target let one big stratum decide the answer.

    Two strata, one target in each. The first stratum's controls never move, the second's always do,
    and the second holds a hundred times more controls. Standardised, the control rate is the mean of
    the two strata's rates, 0.5. Pooled the old way it would be about 0.99.
    """
    targets = [_row("t", True, tss=1_000), _row("t", True, tss=1_000_000)]
    controls = [_row("c", False, tss=1_000) for _ in range(2)]
    controls += [_row("c", True, tss=1_000_000) for _ in range(200)]

    out = matched.matched_pair(targets + controls, "t", "c")

    assert out["matched"]["b"] == 0.5
    assert out["left_in_a_shared_stratum"] == 2


def test_the_raw_comparison_still_pools_because_that_is_what_raw_means() -> None:
    """Raw is the unstandardised rate; only the matched figure holds the strata fixed."""
    targets = [_row("t", True, tss=1_000), _row("t", True, tss=1_000_000)]
    controls = [_row("c", False, tss=1_000) for _ in range(2)]
    controls += [_row("c", True, tss=1_000_000) for _ in range(200)]

    out = matched.matched_pair(targets + controls, "t", "c")

    assert out["raw"]["b"] > 0.98  # 200 of 202 controls move
    assert out["matched"]["b"] == 0.5  # the strata say otherwise


def test_a_target_with_no_control_in_its_stratum_is_dropped_and_counted() -> None:
    inside = _row("t", True, tss=1_000)
    outside = _row("t", True, tss=50_000_000)
    control = _row("c", True, tss=1_000)

    out = matched.matched_pair([inside, outside, control], "t", "c")

    assert out["left_in_a_shared_stratum"] == 1


# `stratum_rates` moved to genomeos.compare on 2026-09-17 and is tested in tests/test_compare.py,
# with the large-stratum regression beside it; the script keeps only its reporting shape.


# --- coverage: what counts as measured -------------------------------------------------------------


def test_overlapping_measurements_count_a_base_once() -> None:
    """Two assays over the same sequence must not make a block look twice as measured."""
    spans = [(100, 200), (150, 250)]

    assert coverage.overlap_bp(spans, 0, 1000) == 150


def test_a_measurement_outside_the_block_does_not_count() -> None:
    assert coverage.overlap_bp([(0, 50)], 100, 200) == 0
    assert coverage.overlap_bp([], 100, 200) == 0


def test_a_measurement_is_clipped_to_the_block() -> None:
    assert coverage.overlap_bp([(50, 500)], 100, 200) == 100


def test_the_real_unknown_keeps_its_case_and_copies_are_separated() -> None:
    real = {"tier": "constrained_unknown", "copy": False, "case": "syntax"}
    copied = {"tier": "constrained_unknown", "copy": True, "case": "syntax"}

    assert coverage.label_of(real) == "real_unknown_syntax"
    assert coverage.label_of(copied) == "constrained_unknown_copy"
    assert coverage.label_of({"tier": "neutral"}) == "neutral"


# --- the library -----------------------------------------------------------------------------------


def test_the_shuffle_keeps_every_dinucleotide() -> None:
    """Composition is the one thing that transfers, so the null has to hold pairs, not bases."""
    import random

    seq = "ACGTACGGCTTAGCATCGACTGCATCGAACGT" * 4

    out = library.dinucleotide_shuffle(seq, random.Random(3))

    def pairs(s: str) -> dict:
        got: dict[str, int] = {}
        for a, b in zip(s, s[1:], strict=False):
            got[a + b] = got.get(a + b, 0) + 1
        return got

    assert len(out) == len(seq)
    assert sorted(out) == sorted(seq)  # every base kept
    before, after = pairs(seq), pairs(out)
    assert sum(before.values()) == sum(after.values())
    # every dinucleotide that leaves a base appears the same number of times
    assert {k: v for k, v in before.items() if k[0] != seq[-1]} == {
        k: v for k, v in after.items() if k[0] != out[-1]
    }


def test_the_shuffle_cannot_change_gc() -> None:
    """The scrambled arm is a composition control, so its GC must equal its source's exactly."""
    import random

    rng = random.Random(11)
    for _ in range(50):
        seq = "".join(rng.choice("ACGT") for _ in range(300))
        out = library.dinucleotide_shuffle(seq, rng)
        assert seq.count("G") + seq.count("C") == out.count("G") + out.count("C")


def test_both_low_complexity_rules_are_charged_as_a_union() -> None:
    """One rule reads the sequence and the other reads an annotation; an oligo failing either is out.

    genomeos-79's rule fires on what RepeatMasker saw, mine on the sequence itself, so an unannotated
    simple repeat passes theirs and fails mine - which is why the union is taken rather than a choice.
    """
    clean = "ACGTTGCAACGGTTACGGATCCGTAAGCTTAGCGATCGATTCAGGCATCAG" * 6
    annotated_only = library.family_a(clean, 0.5, {"Simple_repeat": 200})
    neither = library.family_a(clean, 0.5, {"LINE": 300})

    assert annotated_only == ["low_complexity_by_annotation"]
    assert neither == []
    assert library.family_a("A" * 300, 0.0, None)[:2] == ["homopolymer", "gc_outside_25_75"]


def test_a_block_is_tiled_end_to_end_without_running_past_it() -> None:
    windows = library.oligos_of_block({"start": 1_000, "end": 1_950}, library.OLIGO)

    assert windows == [(1000, 1300), (1300, 1600), (1600, 1900)]
    assert all(end <= 1_950 for _, end in windows)


def test_a_block_shorter_than_an_oligo_yields_nothing() -> None:
    assert library.oligos_of_block({"start": 0, "end": 200}, library.OLIGO) == []
