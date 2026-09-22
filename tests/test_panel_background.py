# SPDX-License-Identifier: AGPL-3.0-or-later
"""What the panel's matched background is made of: the arithmetic, and the two kinds of zero.

Everything here is synthetic: a store written base by base, tracks written line by line, no network
and no committed result read. The chromosome is built so that one region is quiet and one is noisy,
which makes the direction of every decomposition predictable.
"""

import gzip
import json

import pytest

from genomeos.attribution import panel_background as pb

CHROM = "chrT"
LENGTH = 20_000
QUIET_FROM = 15_000  # kilobases 15..19 carry two recurring events each; 0..14 carry ten
ASSEMBLIES = ["hg38", "a1", "a2", "a3"]


def write_store(d):
    """A four-assembly store: one alignment block over the whole chromosome, GC 0.40 everywhere."""
    d.mkdir(parents=True, exist_ok=True)
    with gzip.open(d / "blocks.tsv.gz", "wt") as fh:
        fh.write(f"0\t{LENGTH}\t{'s' * len(ASSEMBLIES)}\t0\n")
    with gzip.open(d / "sites.tsv.gz", "wt") as fh:
        for kb in range(LENGTH // 1000):
            n = 2 if kb * 1000 >= QUIET_FROM else 10
            for i in range(n):  # spread one per 100-base window; two assemblies carry G, so recurring
                fh.write(f"{kb * 1000 + i * (1000 // n) + 3}\tA\t=G=G\n")
    with gzip.open(d / "insertions.tsv.gz", "wt") as fh:
        fh.write("")
    with gzip.open(d / "gc.tsv.gz", "wt") as fh:
        for b in range(LENGTH // 100):
            fh.write(f"{b}\t40\t100\n")
    (d / "meta.json").write_text(json.dumps({"chrom": CHROM, "assemblies": ASSEMBLIES, "source": "test"}))
    (d / "replication_timing.json").write_text(json.dumps({str(b): 1.0 for b in range(LENGTH // 1000)}))


def write_result(results_dir):
    blocks = [
        {"start": 0, "end": 5_000, "tier": "regulatory", "mammal_fraction": 0.01},
        {"start": QUIET_FROM, "end": LENGTH, "tier": "neutral", "mammal_fraction": None},
    ]
    pooled = {"regulatory": {"ratio_gc_rt": 1.25}, "neutral": {"ratio_gc_rt": 0.25}}
    (results_dir / f"human_panel_{CHROM}.json").write_text(json.dumps({"blocks": blocks, "pooled": pooled}))


def write_tracks(d):
    """Conserved elements over the quiet region; segmental duplication read and empty; no rmsk, no cCREs."""
    with gzip.open(d / f"phastConsElements100way_{CHROM}.bed.gz", "wt") as fh:
        fh.write("# synthetic\n")
        for s in range(QUIET_FROM, LENGTH, 1000):
            fh.write(f"{CHROM}\t{s}\t{s + 900}\t100\n")
    with gzip.open(d / f"superdups_{CHROM}.bed.gz", "wt") as fh:
        fh.write("# synthetic: read, and empty\n")


@pytest.fixture()
def out(tmp_path):
    write_store(tmp_path / CHROM)
    write_result(tmp_path)
    write_tracks(tmp_path)
    return pb.analyse(
        CHROM, cache=tmp_path, results_dir=tmp_path, tracks_dir=tmp_path, constraint_dir=tmp_path
    )


# -- the arithmetic ------------------------------------------------------------------------------
def test_cover_cells_counts_each_base_once():
    a = pb.cover_cells([(0, 150), (100, 250), (400, 401)], cells=6)
    assert list(a) == [100, 100, 50, 0, 1, 0]
    assert pb.cells_bp(a, 0, 300) == 250
    assert pb.span_bp(a, 0, 100) == 100
    # a feature never covers more than the span it is asked about
    assert pb.span_bp(a, 120, 130) == 10


def test_nearest_without_points_is_far_not_zero():
    assert pb.nearest([], 500) == pb.FAR
    assert pb.nearest([100, 900], 500) == 400


def test_features_separate_absent_from_never_measured():
    f = pb._features({"exon_any": (0, 1000), "conserved_elements": (250, 1000)}, ["repeats"])
    # measured and there is none, versus never looked at: different sentences, not the same zero
    assert f["exon_any"] == {"share": 0.0, "bases_assessed": 1000, "reads": "measured, and none of it"}
    assert f["conserved_elements"]["share"] == 0.25
    assert f["repeats"] == {
        "share": None,
        "bases_assessed": 0,
        "reads": "not assessed",
        "why": pb.NO_TRACK,
    }


def test_the_decomposition_cell_is_a_sentence_not_a_zero():
    assert pb.reads_explained(None) == "not assessed"
    assert pb.reads_explained(0.0) == "measured, accounts for none"
    assert pb.reads_explained(0.33) == "measured, accounts for 33% of the offset"
    assert pb.reads_explained(-0.34) == "measured, widens the offset by 34%"


def test_explained_says_which_covariates_were_assessed():
    offset = {
        "as_the_panel_builds_it": {"neutral": 1.10},
        "without_exon_any": {"assessed": True, "neutral": 1.05},
        "without_repeats": {"assessed": False, "why": "no track for this chromosome"},
    }
    e = pb.explained(offset, "neutral")
    assert e["without_exon_any"] == {
        "assessed": True,
        "ratio": 1.05,
        "explains": 0.5,
        "reads": "measured, accounts for 50% of the offset",
    }
    assert e["without_repeats"] == {
        "assessed": False,
        "explains": None,
        "reads": "not assessed",
        "why": "no track for this chromosome",
    }


# -- the rebuild ---------------------------------------------------------------------------------
def test_local_tracks_names_what_is_missing(tmp_path):
    write_tracks(tmp_path)
    t = pb.local_tracks(CHROM, tmp_path, tmp_path)
    assert sorted(t["missing"]) == ["registry", "repeats"]
    assert t["segmental_duplication"] == []  # read, and empty: not the same as missing
    assert len(t["conserved_elements"]) == 5


def test_background_is_the_panels_own_kilobases(out):
    bg = out["composition"]["background"]
    assert bg["kilobases"] == 20 and bg["bases"] == 20_000
    assert out["coverage"]["kilobases_in_panel_blocks"] == 20
    assert out["coverage"]["kilobases_not_measured"] == {}
    assert out["coverage"]["kilobases_measured_without_replication_timing"] == 0
    # 15 kilobases of ten events and 5 of two: 160 over 20 kb
    assert bg["recurring_per_kb"] == pytest.approx(8.0, abs=0.01)


def test_the_quiet_region_reads_below_its_own_background(out):
    built = out["offset"]["as_the_panel_builds_it"]
    assert built["regulatory"] == pytest.approx(10 / 8, abs=0.01)
    assert built["neutral"] == pytest.approx(2 / 8, abs=0.01)


def test_the_rebuild_is_checked_against_the_committed_ratios(out):
    """Nothing below means anything if the rebuilt background is not the panel's own."""
    r = out["reproduces_the_panel"]
    assert r["tiers_checked"] == 2 and r["largest_absolute_difference"] < 0.01
    assert r["tiers_not_in_the_committed_result"] == ["cds", "structural", "fossil", "constrained_unknown"]


def test_taking_the_quiet_region_out_raises_the_background_and_moves_the_ratio(out):
    built = out["offset"]["as_the_panel_builds_it"]
    without = out["offset"]["without_conserved_elements"]
    assert without["assessed"] is True
    assert without["background_per_kb"] > built["background_per_kb"]
    # the conserved elements were the quiet part, so removing them moves the noisy tier towards 1
    assert out["explained"]["regulatory"]["without_conserved_elements"]["explains"] > 0


def test_a_covariate_with_no_data_is_never_an_explanation_of_zero(out):
    cov = out["coverage"]
    assert "repeats" in cov["tracks_absent"] and cov["gencode_read"] is False
    assert cov["covariates_not_assessed"]["exon_any"] == "GENCODE not read for this chromosome"
    # read and empty is a different sentence from never looked at
    assert cov["covariates_not_assessed"]["segmental_duplication"] == "read, and empty on this chromosome"
    for name in ("exon_any", "segmental_duplication"):
        assert out["explained"]["regulatory"][f"without_{name}"]["explains"] is None
    feats = out["composition"]["background"]["features"]
    assert feats["repeats"]["share"] is None and feats["repeats"]["bases_assessed"] == 0
    assert feats["segmental_duplication"]["share"] == 0.0
    assert feats["segmental_duplication"]["bases_assessed"] == 20_000
    assert out["feature_splits"]["repeats"]["kilobases_assessed"] == 0
    assert out["feature_splits"]["conserved_elements"]["kilobases_assessed"] == 20


def test_every_group_prints_length_gc_and_distance_to_a_coding_tss(out):
    for name in ("background", "background_outside_every_block", "regulatory", "neutral"):
        row = out["composition"][name]
        assert row["gc_median"] == pytest.approx(0.4, abs=0.001)
        assert row["tss_median"] == pb.FAR  # no GENCODE here, so the distance is not a small number
        assert row["tss_assessed"] == 0


def test_the_two_group_comparison_goes_through_compare(out):
    s = out["standardised"]["neutral"]["gc_and_timing"]
    assert s["targets"] and s["controls"]
    assert s["dropped_for_want_of_a_control"] == 0
    assert set(s["imbalance"]) == {"length", "gc", "rt"}
    assert s["matched"]["difference"] < 0  # the quiet tier carries a recurring event less often
    assert "tss" in out["standardised"]["neutral"]["gc_timing_and_tss"]["imbalance"]


def test_coverage_is_held_fixed_as_well_as_the_covariates(out):
    """Standardising on covariates is not conditioning on effort; both are reported."""
    arms = out["standardised"]["neutral"]
    assert set(arms) == {"gc_and_timing", "gc_timing_and_tss", "gc_timing_tss_and_coverage"}
    held = arms["gc_timing_tss_and_coverage"]
    assert "coverage" in held["imbalance"] and held["imbalance"]["coverage"]["target_median"] == 1.0
    assert out["standardised"]["coverage_measure"] == pb.COVERAGE_MEASURE
    # the store aligns every assembly everywhere, so conditioning on coverage changes nothing here
    assert held["matched"]["difference"] == pytest.approx(arms["gc_timing_and_tss"]["matched"]["difference"])
    assert out["offset"][pb.COVERAGE_KEY]["assessed"] is True
    assert out["offset"][pb.COVERAGE_KEY]["kilobases_left_out"] == 0
    assert out["claim_available"] == pb.CLAIM_WITH_COVERAGE
    assert out["coverage_conditioned"] == {
        "in_the_window_comparison": True,
        "by_rebuilding_the_background": True,
    }


def test_a_background_too_thin_to_condition_says_so_instead_of_returning_a_number(tmp_path):
    """Where nearly every kilobase fails the coverage cut -- a sex chromosome -- the arm reports why."""
    d = tmp_path / CHROM
    write_store(d)
    write_result(tmp_path)
    write_tracks(tmp_path)
    with gzip.open(d / "blocks.tsv.gz", "wt") as fh:  # one assembly of four informs anything
        fh.write(f"0\t{LENGTH}\ts...\t0\n")
    out = pb.analyse(
        CHROM, cache=tmp_path, results_dir=tmp_path, tracks_dir=tmp_path, constraint_dir=tmp_path
    )
    arm = out["offset"][pb.COVERAGE_KEY]
    assert arm["assessed"] is False and arm["why"] == pb.COVERAGE_TOO_THIN
    assert arm["kilobases_left_out"] == arm["kilobases"] == 20
    assert "neutral" not in arm  # no ratio is returned rather than a ratio built from scraps
    assert out["coverage_conditioned"]["by_rebuilding_the_background"] is False
    # the window comparison still holds coverage, so the stronger sentence is still the available one
    assert out["claim_available"] == pb.CLAIM_WITH_COVERAGE


def test_phylop_is_not_sampled_unless_asked(out):
    sampled = out["constrained_bases"]["background_sampled"]
    assert sampled == {"kilobases_sampled": 0, "kilobases_not_sampled": 20, "why": "not requested"}
    tiers = out["constrained_bases"]["tiers_from_the_budget"]
    assert tiers["neutral"] == {"blocks": 1, "blocks_assessed": 0, "constrained_fraction": None}
    assert tiers["regulatory"]["blocks_assessed"] == 1
