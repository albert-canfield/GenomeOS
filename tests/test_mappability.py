# SPDX-License-Identifier: AGPL-3.0-or-later
"""Measured mappability in place of the repeat proxy: the track read by range from a bigWig built in
memory, missing values counted as the zeros Umap means them to be, the rule at three values, and the
two arms cross-tabulated without ever being added together. No network."""

from __future__ import annotations

import importlib.util
import struct
import zlib
from pathlib import Path

import pytest

from genomeos.attribution import bigwig as bw
from genomeos.attribution import mappability as mp

ROOT = Path(__file__).resolve().parents[1]


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


script = _load("mappability")


# --- a bigWig small enough to build by hand ----------------------------------------------------


def _bigwig(sections: list[tuple[str, int, list[float]]], chrom_size: int = 100_000) -> bytes:
    """A little dense fixedStep bigWig: one section per (chromosome, start, values), step and span 1.

    The same layout the reader's own tests use, widened to several chromosomes so a chromosome that
    is absent from the track can be told from one whose values are absent.
    """
    names = sorted({name for name, _, _ in sections})
    ids = {name: i for i, name in enumerate(names)}
    key_size = max(len(n) for n in names)

    chrom_tree = struct.pack("<IIIIQQ", bw.CHROM_TREE_MAGIC, len(names), key_size, 8, len(names), 0)
    chrom_tree += struct.pack("<BBH", 1, 0, len(names))
    for name in names:
        chrom_tree += name.encode().ljust(key_size, b"\0") + struct.pack("<II", ids[name], chrom_size)

    blobs, meta = [], []
    for name, start, values in sections:
        raw = struct.pack(
            "<IIIIIBBH", ids[name], start, start + len(values), 1, 1, bw.FIXED_STEP, 0, len(values)
        ) + struct.pack(f"<{len(values)}f", *values)
        blobs.append(zlib.compress(raw))
        meta.append((ids[name], start, start + len(values)))

    chrom_tree_at = 64
    data_at = chrom_tree_at + len(chrom_tree)
    data = b"".join(blobs)
    index_at = data_at + len(data)

    offsets, at = [], data_at
    for blob in blobs:
        offsets.append((at, len(blob)))
        at += len(blob)

    block_size = max(2, len(sections))
    rtree = struct.pack(
        "<IIQIIIIQII",
        bw.RTREE_MAGIC,
        block_size,
        len(sections),
        0,
        0,
        len(names) - 1,
        chrom_size,
        index_at,
        1,
        0,
    )
    leaf = struct.pack("<BBH", 1, 0, len(sections))
    for (cid, lo, hi), (off, size) in zip(meta, offsets, strict=True):
        leaf += struct.pack("<IIIIQQ", cid, lo, cid, hi, off, size)
    leaf += b"\0" * (4 + block_size * 32)  # the reader reads a whole node's worth; give it room

    header = struct.pack(
        "<IHHQQQHHQQIQ", bw.BIGWIG_MAGIC, 4, 0, chrom_tree_at, data_at, index_at, 0, 0, 0, 0, 1 << 15, 0
    )
    assert len(header) == 64
    return header + chrom_tree + data + rtree + leaf


@pytest.fixture
def track(tmp_path):
    """chr1 1000-1600: the first 300 bases all 1.0, the next 300 half 1.0 and half 0.5.
    chr2 2000-2150: values only over half of one 300 bp window, so the rest is a measured zero."""
    values = [1.0] * 300 + [1.0] * 150 + [0.5] * 150
    p = tmp_path / "k24.bw"
    p.write_bytes(_bigwig([("chr1", 1000, values), ("chr2", 2000, [1.0] * 150)]))
    return p


def _rows(track, chrom, intervals, **kw):
    res = mp.read_chromosome(24, chrom, intervals, source=track, cache_dir=None, **kw)
    return res, [res.rows[iv] for iv in intervals if iv in res.rows]


# --- reading -----------------------------------------------------------------------------------


def test_the_track_is_read_by_range_and_counts_bases_at_the_cut(track):
    res, rows = _rows(track, "chr1", [(1000, 1300), (1300, 1600)])
    assert res.in_track and not res.cap_reached and res.unread == 0
    assert res.bytes_fetched > 0 and res.requests > 0
    assert rows[0]["bases"] == 300 and rows[0]["above"] == 300
    assert rows[1]["bases"] == 300 and rows[1]["above"] == 150
    assert mp.unique_fraction(rows[0]) == 1.0
    assert mp.unique_fraction(rows[1]) == 0.5


def test_a_base_the_track_has_no_value_for_is_a_measured_zero_not_a_missing_one(track):
    """Umap writes nothing where mappability is zero, so the denominator is the oligo, not the
    bases that happen to carry a value. Half a window covered at 1.0 is half mappable."""
    res, rows = _rows(track, "chr2", [(2000, 2300)])
    assert rows[0]["bases"] == 150 and rows[0]["above"] == 150
    assert mp.unique_fraction(rows[0]) == 0.5
    assert mp.mean_mappability(rows[0]) == 0.5
    assert res.in_track


def test_a_window_with_no_value_at_all_is_unmappable_and_still_assessed(track):
    res, rows = _rows(track, "chr1", [(50_000, 50_300)])
    assert res.in_track and res.unread == 0
    assert rows[0]["bases"] == 0 and mp.unique_fraction(rows[0]) == 0.0
    assert not mp.is_mappable(rows[0], 0.5)


def test_a_chromosome_absent_from_the_track_is_named_not_zeroed(track):
    res = mp.read_chromosome(24, "chrZ", [(1000, 1300)], source=track, cache_dir=None)
    assert res.in_track is False
    assert res.rows == {} and res.unread == 1


def test_the_byte_cap_stops_the_read_and_counts_what_it_did_not_reach(track):
    res = mp.read_chromosome(
        24, "chr1", [(1000, 1300), (1300, 1600)], source=track, cache_dir=None, byte_cap=1
    )
    assert res.cap_reached and res.unread == 2 and res.rows == {}


def test_overlapping_and_repeated_intervals_are_de_duplicated_before_the_reader_sees_them(track):
    res = mp.read_chromosome(
        24, "chr1", [(1000, 1300), (1000, 1300), (1300, 1600)], source=track, cache_dir=None
    )
    assert len(res.rows) == 2


# --- the cache ---------------------------------------------------------------------------------


def test_the_cache_round_trips_and_a_second_read_costs_no_bytes(track, tmp_path):
    cache = tmp_path / "cache"
    first = mp.read_chromosome(24, "chr1", [(1000, 1300), (1300, 1600)], source=track, cache_dir=cache)
    assert first.bytes_fetched > 0 and first.from_cache == 0
    second = mp.read_chromosome(24, "chr1", [(1000, 1300), (1300, 1600)], source=track, cache_dir=cache)
    assert second.bytes_fetched == 0 and second.requests == 0
    assert second.from_cache == 2
    assert {k: dict(v) for k, v in second.rows.items()} == {k: dict(v) for k, v in first.rows.items()}
    assert mp.cache_bytes(cache) == mp.cache_path(24, "chr1", cache).stat().st_size


def test_the_cache_refuses_to_pass_its_cap(track, tmp_path):
    cache = tmp_path / "cache"
    rows = {(1000, 1300): {"bases": 300, "above": 300, "total": 300.0}}
    assert mp.save_cache(24, "chr1", rows, cache, cap=10 << 20) is True
    assert mp.save_cache(24, "chr2", rows, cache, cap=10) is False
    assert not mp.cache_path(24, "chr2", cache).exists()


def test_each_k_has_its_own_track_and_no_two_share_a_file(track):
    assert sorted(mp.TRACKS) == [24, 36, 50, 100]
    for k, t in mp.TRACKS.items():
        assert t.k == k and f"k{k}." in t.url and t.assembly.startswith("GRCh38")
        assert t.as_dict()["k"] == k
    assert len({t.url for t in mp.TRACKS.values()}) == 4
    assert len({t.size_bytes for t in mp.TRACKS.values()}) == 4


# --- the rule ----------------------------------------------------------------------------------


def _oligo(above: int, interspersed: float, segdup: float = 0.0, gc: float = 0.40, tss: int = 50_000):
    return {
        "start": 0,
        "end": 300,
        "bases": 300,
        "above": above,
        "total": float(above),
        "interspersed_fraction": interspersed,
        "segdup_fraction": segdup,
        "block_length": 5_000,
        "gc": gc,
        "tss_distance": tss,
    }


def test_the_rule_is_a_fraction_of_the_oligo_and_moves_with_its_named_constant():
    row = _oligo(above=270, interspersed=0.0)  # 90% of the oligo uniquely readable
    assert mp.is_mappable(row, 0.50) and mp.is_mappable(row, 0.90)
    assert not mp.is_mappable(row, 0.99)
    assert mp.MAPPABLE_FRACTION == 0.90 and mp.BASE_CUT == 1.0
    assert mp.SENSITIVITY_FRACTIONS == (0.50, 0.90, 0.99)


def test_the_proxy_is_family_b_and_either_fraction_can_raise_it():
    assert mp.proxy_flagged(_oligo(300, interspersed=0.6))
    assert mp.proxy_flagged(_oligo(300, interspersed=0.0, segdup=0.7))
    assert not mp.proxy_flagged(_oligo(300, interspersed=0.49, segdup=0.49))


def test_the_sensitivity_count_can_only_fall_as_the_rule_tightens():
    rows = [_oligo(a, 0.0) for a in (300, 290, 250, 150, 0)]
    counts = [s["attributable_by_mappability_track"] for s in mp.sensitivity(rows)]
    assert counts == sorted(counts, reverse=True)
    assert [s["min_fraction"] for s in mp.sensitivity(rows)] == list(mp.SENSITIVITY_FRACTIONS)


# --- the cross-tabulation ----------------------------------------------------------------------


def test_every_oligo_lands_in_exactly_one_of_the_four_named_cells():
    rows = [
        _oligo(300, 0.9),  # proxy flags, track mappable
        _oligo(0, 0.9),  # both
        _oligo(0, 0.0),  # the dangerous cell: proxy passes, track unmappable
        _oligo(300, 0.0),  # neither
    ]
    x = mp.crosstab(rows)
    assert [x[c] for c in mp.CELLS] == [1, 1, 1, 1]
    assert sum(x[c] for c in mp.CELLS) == len(rows) == x["oligos"]
    assert x["agreement"] == 0.5
    assert x["proxy_passes_track_unmappable"] == 1


def test_the_two_arms_are_three_named_numbers_and_never_one():
    rows = [_oligo(300, 0.9), _oligo(0, 0.9), _oligo(0, 0.0), _oligo(300, 0.0)]
    a = mp.arms(rows)
    assert a["attributable_by_repeat_proxy"] == 2
    assert a["attributable_by_mappability_track"] == 2
    assert a["attributable_by_both_arms"] == 1
    assert a["attributable_by_both_arms"] <= min(
        a["attributable_by_repeat_proxy"], a["attributable_by_mappability_track"]
    )
    assert "attributable" not in a  # no merged headline to quote by mistake


def test_the_proxy_swing_is_reported_against_a_track_that_does_not_move_with_it():
    rows = [_oligo(300, f) for f in (0.1, 0.3, 0.6, 0.9)]
    out = mp.proxy_sensitivity(rows)
    assert [o["interspersed_max"] for o in out] == list(mp.PROXY_INTERSPERSED)
    assert [o["attributable_by_repeat_proxy"] for o in out] == [1, 2, 3]
    assert {o["attributable_by_mappability_track"] for o in out} == {4}


# --- groups and the matched comparison ---------------------------------------------------------


def test_every_group_prints_length_gc_and_distance_to_a_coding_tss():
    rows = [_oligo(300, 0.9), _oligo(0, 0.9), _oligo(0, 0.0), _oligo(300, 0.0)]
    g = mp.groups(rows)
    for name in (*mp.CELLS, "track_mappable", "track_unmappable", "proxy_flagged", "proxy_passed"):
        cov = g[name]
        assert cov["oligo_bp"] == mp.OLIGO
        assert set(cov) >= {
            "oligos",
            "median_block_length",
            "block_length_quartiles",
            "median_gc",
            "median_tss_distance",
            "tss_distance_quartiles",
            "without_tss_distance",
        }
    assert g["track_mappable"]["oligos"] + g["track_unmappable"]["oligos"] == len(rows)


def test_the_comparison_comes_from_compare_with_its_matched_n_and_its_imbalance():
    rows = [_oligo(300, 0.0, gc=0.40, tss=50_000) for _ in range(20)]
    rows += [_oligo(0, 0.9, gc=0.40, tss=50_000) for _ in range(10)]
    rows += [_oligo(300, 0.9, gc=0.40, tss=50_000) for _ in range(10)]
    out = mp.proxy_against_track(rows)
    assert out["targets"] == 20 and out["controls"] == 20
    assert out["targets_matched"] == 20 and out["dropped_for_want_of_a_control"] == 0
    assert out["matched"]["difference"] == pytest.approx(0.5)
    assert set(out["imbalance"]) == {"length", "gc", "tss_distance"}
    assert out["excluded_without_a_covariate"] == 0


def test_an_oligo_without_a_covariate_is_excluded_from_the_matching_and_counted():
    rows = [_oligo(300, 0.0) for _ in range(5)] + [_oligo(0, 0.9, tss=None) for _ in range(3)]
    out = mp.proxy_against_track(rows)
    assert out["excluded_without_a_covariate"] == 3
    assert out["targets"] == 0


# --- the script's joins ------------------------------------------------------------------------


def test_the_manifest_reader_counts_the_checks_the_zero_reading_rests_on(tmp_path):
    path = tmp_path / "lib.tsv"
    seq = "ACGT" * 75
    path.write_text(
        "id\tarm\tchrom\tstart\tend\tblock\tgc\ttss_distance\tsequence\n"
        f"o1\ttest\tchr1\t1000\t1300\tchr1:900-9900\t0.5\t4000\t{seq}\n"
        f"o2\ttest\tchr1\t1300\t1600\tchr1:900-9900\t0.5\t\t{seq}\n"
        f"o3\tscrambled\tchr1\t1000\t1300\tchr1:900-9900\t0.5\t4000\t{seq}\n"
    )
    by_chrom, checks = script.load_manifest(path, "test")
    assert checks["oligos"] == 2 and checks["n_bases_in_arm"] == 0
    assert checks["oligos_not_one_oligo_long"] == 0
    assert [o["block_length"] for o in by_chrom["chr1"]] == [9000, 9000]
    assert by_chrom["chr1"][1]["tss_distance"] is None


def test_an_n_base_in_the_arm_is_counted_because_it_would_break_the_zero_reading(tmp_path):
    path = tmp_path / "lib.tsv"
    path.write_text(
        "id\tarm\tchrom\tstart\tend\tblock\tgc\ttss_distance\tsequence\n"
        "o1\ttest\tchr1\t1000\t1300\tchr1:900-9900\t0.5\t4000\t" + "N" * 300 + "\n"
    )
    _, checks = script.load_manifest(path, "test")
    assert checks["n_bases_in_arm"] == 300


def test_a_block_string_parses_back_to_its_span():
    assert script.block_span("chr21:15067837-15071837") == ("chr21", 15067837, 15071837)
