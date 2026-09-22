"""Measurability of the untouched real unknown: the sequence rules, the tiling arithmetic and the
base accounting.

Everything here is synthetic. No network, no reference genome, no cached track.
"""

import gzip
import random

from genomeos.attribution.measurability import (
    ARITHMETIC,
    FAMILY_A,
    FAMILY_B,
    KMER,
    OLIGO,
    REASONS,
    STEP,
    Thresholds,
    assess_block,
    block_counts,
    calibration,
    charged_reason,
    class_groups,
    covariates,
    duplicated_kmer_fraction,
    failed_reasons,
    fraction_covered,
    gc_fraction,
    kmer_counts,
    longest_homopolymer,
    n_fraction,
    nearest,
    pooled,
    real_unknown_label,
    repeat_spans_by_class,
    sample,
    sensitivity,
    tiling,
    verdict,
    window_covariates,
    window_metrics,
)

EMPTY = {"tandem": [], "interspersed": [], "other_repeat": []}


def _metrics(**kw):
    """A window's measurements with everything passing unless overridden."""
    base = {
        "start": 0,
        "end": OLIGO,
        "n_fraction": 0.0,
        "gc": 0.45,
        "homopolymer": 4,
        "tandem_fraction": 0.0,
        "interspersed_fraction": 0.0,
        "segdup_fraction": 0.0,
        "young_segdup_fraction": 0.0,
        "dup_kmer_fraction": 0.0,
    }
    return {**base, **kw}


# --- sequence measurements --------------------------------------------------------------------


def test_n_fraction_counts_every_uncalled_base():
    assert n_fraction("ACGT") == 0.0
    assert n_fraction("ACGN") == 0.25
    assert n_fraction("NNNN") == 1.0
    assert n_fraction("") == 1.0  # nothing is not measurable either


def test_gc_fraction_ignores_n_and_is_none_without_called_bases():
    assert gc_fraction("GCGC") == 1.0
    assert gc_fraction("ATAT") == 0.0
    assert gc_fraction("GCNN") == 1.0  # N is not counted in the denominator
    assert gc_fraction("NNNN") is None


def test_longest_homopolymer_finds_the_run_anywhere():
    assert longest_homopolymer("ACGT") == 1
    assert longest_homopolymer("AAACGT") == 3
    assert longest_homopolymer("ACGTTTTTTTTTT") == 10
    assert longest_homopolymer("") == 0


def test_duplicated_kmer_fraction_separates_a_tandem_copy_from_unique_sequence():
    rng = random.Random(1)
    unique = "".join(rng.choice("ACGT") for _ in range(2000))
    assert duplicated_kmer_fraction(unique, 0, 300, KMER) == 0.0
    # the same 300 bp twice: every k-mer of the first copy occurs again in the second
    doubled = unique[:300] + unique[:300] + unique[300:]
    assert duplicated_kmer_fraction(doubled, 0, 300, KMER) > 0.8
    assert duplicated_kmer_fraction(doubled, 300, 300, KMER) > 0.8
    # a window past both copies is unaffected
    assert duplicated_kmer_fraction(doubled, 900, 300, KMER) == 0.0


def test_kmer_counts_counts_every_position():
    counts = kmer_counts("AAAA", 2)
    assert counts == {"AA": 3}


def test_duplicated_kmer_fraction_is_zero_when_the_window_is_shorter_than_k():
    assert duplicated_kmer_fraction("ACGT" * 5, 0, 10, KMER) == 0.0


# --- interval measurements --------------------------------------------------------------------


def test_fraction_covered_is_the_share_of_the_window():
    merged = [(100, 200), (300, 350)]
    assert fraction_covered(merged, 0, 100) == 0.0
    assert fraction_covered(merged, 100, 200) == 1.0
    assert fraction_covered(merged, 150, 250) == 0.5
    assert fraction_covered(merged, 0, 0) == 0.0


def test_class_groups_splits_tandem_from_interspersed_and_keeps_the_rest():
    groups = class_groups(
        {
            "Simple_repeat": [(0, 10)],
            "Satellite": [(10, 20)],
            "LINE": [(100, 200)],
            "SINE": [(150, 250)],
            "tRNA": [(400, 410)],
        }
    )
    assert groups["tandem"] == [(0, 20)]  # adjacent tandem classes merge
    assert groups["interspersed"] == [(100, 250)]  # overlapping LINE and SINE merge
    assert groups["other_repeat"] == [(400, 410)]


def test_repeat_spans_by_class_streams_only_the_regions_asked_for(tmp_path):
    p = tmp_path / "rmsk_chrT.bed.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("0\t100\tLINE\tL1\tL1X\t0.100\n")
        fh.write("5000\t5100\tSINE\tAlu\tAluY\t0.050\n")
    got = repeat_spans_by_class("chrT", [(0, 200)], tmp_path)
    assert got == {"LINE": [(0, 100)]}  # the distant SINE is not loaded
    assert repeat_spans_by_class("chrT", [], tmp_path) == {}
    assert repeat_spans_by_class("chrMissing", [(0, 10)], tmp_path) is None


def test_nearest_returns_none_without_points():
    assert nearest([], 10) is None
    assert nearest([0, 100, 500], 90) == 10


# --- the exclusion rules ----------------------------------------------------------------------


def test_a_clean_window_is_usable_and_lists_no_reason():
    m = _metrics()
    assert failed_reasons(m) == []
    assert charged_reason(m) == "usable"
    assert verdict(m) == "usable"


def test_every_reason_fires_on_its_own_measurement():
    assert failed_reasons(_metrics(n_fraction=0.5)) == ["assembly_gap"]
    assert failed_reasons(_metrics(gc=0.2)) == ["gc_extreme"]
    assert failed_reasons(_metrics(gc=0.9)) == ["gc_extreme"]
    assert failed_reasons(_metrics(gc=None)) == ["gc_extreme"]
    assert failed_reasons(_metrics(homopolymer=10)) == ["homopolymer"]
    assert failed_reasons(_metrics(tandem_fraction=0.6)) == ["tandem_low_complexity"]
    assert failed_reasons(_metrics(dup_kmer_fraction=0.7)) == ["non_unique_in_block"]
    assert failed_reasons(_metrics(segdup_fraction=0.8)) == ["segmental_duplication"]
    assert failed_reasons(_metrics(interspersed_fraction=0.8)) == ["interspersed_repeat"]


def test_a_window_failing_several_is_charged_to_the_first_in_precedence():
    m = _metrics(n_fraction=0.5, gc=0.1, homopolymer=20, interspersed_fraction=1.0)
    assert failed_reasons(m) == ["assembly_gap", "gc_extreme", "homopolymer", "interspersed_repeat"]
    assert charged_reason(m) == "assembly_gap"


def test_family_b_alone_is_unattributable_not_unsynthesisable():
    assert verdict(_metrics(interspersed_fraction=1.0)) == "unattributable"
    assert verdict(_metrics(segdup_fraction=1.0)) == "unattributable"
    assert verdict(_metrics(homopolymer=30, interspersed_fraction=1.0)) == "unsynthesisable"


def test_thresholds_can_be_moved_without_touching_the_measurements():
    m = _metrics(gc=0.28)
    assert charged_reason(m) == "usable"
    assert charged_reason(m, Thresholds(gc_low=0.30)) == "gc_extreme"


def test_the_reason_families_are_disjoint_and_cover_the_named_reasons():
    assert set(FAMILY_A) | set(FAMILY_B) == set(REASONS)
    assert not set(FAMILY_A) & set(FAMILY_B)
    assert not set(ARITHMETIC) & set(REASONS)


# --- tiling arithmetic ------------------------------------------------------------------------


def test_tiling_only_returns_windows_that_fit_whole():
    assert tiling(0, 1000, 300, 300) == [(0, 300), (300, 600), (600, 900)]
    assert tiling(0, 299, 300, 300) == []
    assert tiling(0, 300, 300, 300) == [(0, 300)]
    assert tiling(0, 600, 300, 150) == [(0, 300), (150, 450), (300, 600)]  # a stated overlap


def test_window_metrics_reads_one_window_against_the_tracks():
    seq = "GC" * 150
    m = window_metrics(
        seq, 1000, 1300, {"tandem": [(1000, 1150)], "interspersed": []}, [(1000, 1300)], [], 0.0
    )
    assert m["gc"] == 1.0
    assert m["n_fraction"] == 0.0
    assert m["tandem_fraction"] == 0.5
    assert m["segdup_fraction"] == 1.0
    assert m["interspersed_fraction"] == 0.0


def _block(start=0, end=1000, chrom="chrT", case="syntax"):
    return {"chrom": chrom, "start": start, "end": end, "case": case}


def test_assess_block_tiles_and_measures_every_window():
    rng = random.Random(7)
    seq = "".join(rng.choice("ACGT") for _ in range(1000))
    row = assess_block(_block(), seq, EMPTY, [], [], 12_345)
    assert row["oligos"] == 3
    assert row["tiled_bp"] == 900
    assert row["remainder_bp"] == 100
    assert row["tss_distance"] == 12_345
    assert len(row["windows"]) == 3
    assert row["repeat_fraction"] == 0.0


def test_assess_block_reports_repeat_content_by_group():
    seq = "ACGT" * 250
    groups = class_groups({"Simple_repeat": [(0, 500)], "LINE": [(500, 750)], "tRNA": [(900, 1000)]})
    row = assess_block(_block(), seq, groups, [(0, 100)], [(0, 50)], None)
    assert row["tandem_fraction"] == 0.5
    assert row["interspersed_fraction"] == 0.25
    assert row["other_repeat_fraction"] == 0.1
    assert row["repeat_fraction"] == 0.85
    assert row["segdup_fraction"] == 0.1
    assert row["young_segdup_fraction"] == 0.05


def test_a_block_shorter_than_one_oligo_yields_no_window():
    row = assess_block(_block(end=250), "A" * 250, EMPTY, [], [], None)
    assert row["oligos"] == 0
    assert row["remainder_bp"] == 250


def test_block_counts_keeps_family_b_synthesisable():
    rng = random.Random(3)
    seq = "".join(rng.choice("ACGT") for _ in range(900))
    row = assess_block(_block(end=900), seq, {"tandem": [], "interspersed": [(0, 900)]}, [], [], None)
    c = block_counts(row)
    assert c["oligos"] == 3
    assert c["synthesisable"] == 3  # a repeat-derived oligo can still be ordered
    assert c["attributable"] == 0  # but its answer belongs to a family, not to a locus
    assert c["charged"]["interspersed_repeat"] == 3


# --- pooled accounting ------------------------------------------------------------------------


def _rows():
    rng = random.Random(11)
    clean = "".join(rng.choice("ACGT") for _ in range(1000))
    return [
        assess_block(_block(chrom="chrT"), clean, EMPTY, [], [], 5_000),
        # all-repeat block: three windows charged to Family B
        assess_block(
            _block(start=2000, end=2900, chrom="chrT"),
            clean[:900],
            {"tandem": [], "interspersed": [(2000, 2900)]},
            [],
            [],
            50_000,
        ),
        # below one oligo: nothing tiles
        assess_block(_block(start=5000, end=5250, chrom="chrT"), clean[:250], EMPTY, [], [], 1_000),
    ]


def test_pooled_accounts_for_every_base_exactly_once():
    p = pooled(_rows())
    assert p["bp"] == 1000 + 900 + 250
    assert p["bp_unaccounted"] == 0
    assert p["bp_accounted"] == p["bp"]
    assert p["excluded_bp"]["shorter_than_one_oligo"] == 250
    assert p["excluded_bp"]["tiling_remainder"] == 100
    assert p["oligos_tiled"] == 6
    assert p["oligos_synthesisable"] == 6
    assert p["oligos_attributable"] == 3
    assert p["charged_windows"]["interspersed_repeat"] == 3
    assert p["overlap"] == OLIGO - STEP


def test_pooled_names_every_reason_even_at_zero():
    p = pooled(_rows())
    assert set(p["charged_windows"]) == set(REASONS)
    assert set(p["excluded_bp"]) == set(REASONS) | set(ARITHMETIC)
    assert p["charged_windows"]["assembly_gap"] == 0  # named and counted, not absent


def test_pooled_flags_blocks_that_cannot_carry_a_tiling():
    p = pooled(_rows(), min_usable=3)
    # the repeat block and the short block both fall below three usable oligos
    assert p["blocks_below_min_usable"] == 2
    assert pooled(_rows(), min_usable=0)["blocks_below_min_usable"] == 0


def test_raw_reason_counts_are_at_least_the_charged_ones():
    p = pooled(_rows())
    for name in REASONS:
        assert p["reasons_raw_windows"][name] >= p["charged_windows"][name]


def test_sensitivity_loosening_a_threshold_never_loses_usable_oligos():
    rows = _rows()
    got = {r["threshold"]: r for r in sensitivity(rows)}
    assert got["gc 0.2-0.8"]["oligos_attributable"] >= got["gc 0.3-0.7"]["oligos_attributable"]
    assert got["homopolymer 12"]["oligos_synthesisable"] >= got["homopolymer 8"]["oligos_synthesisable"]
    assert (
        got["interspersed_max 0.75"]["oligos_attributable"]
        >= got["interspersed_max 0.25"]["oligos_attributable"]
    )


def test_covariates_print_length_gc_and_tss_for_any_group():
    c = covariates(_rows())
    assert c["n"] == 3
    assert c["median_length"] == 900
    assert c["median_tss_distance"] == 5_000
    assert c["without_tss_distance"] == 0
    assert c["length_quartiles"] == [250, 900, 900]  # nearest rank, truncating


def test_covariates_count_a_group_with_no_tss_annotation():
    rows = [assess_block(_block(), "ACGT" * 250, EMPTY, [], [], None)]
    c = covariates(rows)
    assert c["median_tss_distance"] is None
    assert c["without_tss_distance"] == 1


def test_window_covariates_name_every_reason_group():
    w = window_covariates(_rows())
    assert set(w) == {"usable", *REASONS}
    assert w["interspersed_repeat"]["windows"] == 3
    assert w["assembly_gap"]["windows"] == 0
    assert w["usable"]["bp"] == 3 * OLIGO


# --- the calibration against sequence already measured ----------------------------------------


def test_calibration_reports_the_rate_the_rules_would_have_thrown_away():
    rows = [{"tss_distance": 1000, "window": _metrics()} for _ in range(8)]
    rows += [{"tss_distance": 1000, "window": _metrics(homopolymer=15)}]
    rows += [{"tss_distance": 1000, "window": _metrics(interspersed_fraction=1.0)}]
    c = calibration(rows)
    assert c["windows"] == 10
    assert c["excluded_family_a"] == 1
    assert c["false_exclusion_rate"] == 0.1
    assert c["flagged_family_b"] == 1
    assert c["family_b_rate"] == 0.1
    assert c["charged"]["usable"] == 8
    assert c["not_calibrated"] == ["non_unique_in_block"]


def test_calibration_of_nothing_is_none_rather_than_zero():
    c = calibration([])
    assert c["windows"] == 0
    assert c["false_exclusion_rate"] is None


# --- the set this measures --------------------------------------------------------------------


def test_real_unknown_label_excludes_copies_and_other_tiers():
    assert real_unknown_label({"tier": "constrained_unknown", "copy": False, "case": "syntax"}) == (
        "real_unknown_syntax"
    )
    assert real_unknown_label({"tier": "constrained_unknown", "copy": True, "case": "syntax"}) is None
    assert real_unknown_label({"tier": "neutral", "copy": False, "case": "syntax"}) is None
    assert real_unknown_label({"tier": "constrained_unknown", "copy": False, "case": None}) == (
        "real_unknown_uncased"
    )


def test_sample_is_deterministic_and_keeps_everything_under_the_cap():
    items = list(range(50))
    assert sample(items, 100) == items
    assert sample(items, 10) == sample(items, 10)
    assert len(sample(items, 10)) == 10
