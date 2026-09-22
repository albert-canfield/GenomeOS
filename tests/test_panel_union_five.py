# SPDX-License-Identifier: AGPL-3.0-or-later
"""The five-way union: the covariate that runs the other way, and the sign it must keep.

Everything here is synthetic -- a store written base by base, a GFF3, a conserved-element track and a
segmental-duplication track written line by line, no network and no committed result read. The
fixtures are `test_panel_leftover.py`'s, loaded from that file rather than copied, and the
conserved-element track is `test_panel_union.py`'s, so a chromosome whose classes are fixed by
construction in one place is the same chromosome in all three.

Unlike the four-way lane's fixture, this one fixes the SIGN of the fifth covariate by construction:
the duplication track is laid over sequence the store makes noisy (ten recurring events per kilobase
against the gene bodies' two), so the kilobases holding it must read noisier than the kilobases that
do not, and taking it out of the background must widen the offset rather than close it. A run that
comes out the other way is a defect and not weather.

The two things that could make this lane a lie are the two things tested hardest: a widening
covariate averaged into a share as though it narrowed the offset, and a five-way arm that quietly
stops being the four-way arm plus one.
"""

import gzip
import importlib.util
from inspect import signature
from pathlib import Path

import pytest

from genomeos.attribution import panel_background as pb
from genomeos.attribution import panel_leftover as pl
from genomeos.attribution import panel_union as pu
from genomeos.attribution import panel_union_five as p5

_spec = importlib.util.spec_from_file_location("tpl", Path("tests/test_panel_leftover.py"))
tpl = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(tpl)

_union = importlib.util.spec_from_file_location("tpu", Path("tests/test_panel_union.py"))
tpu = importlib.util.module_from_spec(_union)
_union.loader.exec_module(tpu)

_script = importlib.util.spec_from_file_location("p5s", Path("scripts/panel_union_five.py"))
p5s = importlib.util.module_from_spec(_script)
_script.loader.exec_module(p5s)

CHROM = tpl.CHROM
LENGTH = tpl.LENGTH
# the duplication track lies inside the last block, 35,000..40,000, which the store makes noisy: ten
# recurring events per kilobase there against two inside a gene body. The direction of every reading
# about this covariate is fixed by that and by nothing else.
SUPERDUPS = [(36_000, 39_000)]


def write_superdups(tracks: Path):
    tracks.mkdir(parents=True, exist_ok=True)
    with gzip.open(tracks / f"superdups_{CHROM}.bed.gz", "wt") as fh:
        for s, e in SUPERDUPS:
            fh.write(f"{s}\t{e}\t{CHROM}\tdup\n")


@pytest.fixture()
def gff(tmp_path):
    p = tmp_path / "synthetic.gff3.gz"
    tpl.write_gff(p)
    return p


@pytest.fixture()
def tracks(tmp_path, gff):
    """Both local tracks where `panel_background.local_tracks` looks for them."""
    write_superdups(tmp_path)
    tpu.write_conserved(tmp_path / "constraint")
    return pb.local_tracks(CHROM, tmp_path, tmp_path / "constraint")


@pytest.fixture()
def out(tmp_path, gff, tracks):
    tpl.write_store(tmp_path / CHROM)
    model = pl.gene_model(CHROM, gff)
    cut = pl.replay_the_cut(model["spans_in_start_order"], LENGTH)
    tpl.write_budget(tmp_path, cut["blocks"])
    tpl.write_panel(tmp_path, cut["blocks"])
    return p5.analyse(
        CHROM,
        cache=tmp_path,
        results_dir=tmp_path,
        gff=gff,
        tracks_dir=tmp_path,
        constraint_dir=tmp_path / "constraint",
    )


# -- the fifth covariate is the first lane's own, and the other four are the fourth lane's ---------
def test_the_five_are_the_four_and_the_one_the_first_lane_measured(gff, tracks):
    ivs = p5.covariate_intervals(p5.covariate_model(CHROM, gff), tracks)
    four = pu.covariate_intervals(pu.covariate_model(CHROM, gff), tracks)
    for name in pu.COVARIATES:
        assert ivs[name] == four[name], name
    # the fifth comes from the track panel_background reads under that name, not from a new parse
    assert ivs[p5.FIFTH] == tracks[p5.FIFTH] == pb._ivs(tracks, {}, p5.FIFTH)
    assert (*pu.COVARIATES, p5.FIFTH) == p5.COVARIATES
    assert p5.WIDENING_WHEN_MEASURED == (p5.FIFTH,)


def test_a_missing_track_and_an_empty_one_are_different_answers_for_the_fifth(gff):
    model = p5.covariate_model(CHROM, gff)
    assert p5.why_absent(p5.FIFTH, model, {"missing": [p5.FIFTH]}) == p5.NO_TRACK
    assert p5.why_absent(p5.FIFTH, model, {"missing": []}) == p5.MEASURED_EMPTY
    # and the four still answer through the lane that defined them
    assert p5.why_absent("promoter", {}, {"missing": []}) == pu.why_absent("promoter", {}, {"missing": []})


def test_the_floor_under_a_share_has_one_home(gff):
    """`MOVED_AT_ALL` is the floor `panel_background.reads_explained` draws, not a second one."""
    assert signature(pb.reads_explained).parameters["floor"].default == p5.MOVED_AT_ALL


# -- the overlaps: the four-way lane's own numbers, and the four pairs the fifth adds --------------
def test_the_four_way_overlaps_are_the_four_way_lanes_own(gff, tracks):
    ivs = p5.covariate_intervals(p5.covariate_model(CHROM, gff), tracks)
    mine = p5.pairwise_overlap(ivs, LENGTH)
    theirs = pu.pairwise_overlap({n: ivs[n] for n in pu.COVARIATES}, LENGTH)
    for key, row in theirs["pairs"].items():
        assert mine["pairs"][key] == row, key
    assert mine["the_four_way"]["union_bases"] == theirs["union_bases"]
    assert mine["the_four_way"]["sum_of_the_singles_bases"] == theirs["sum_of_the_singles_bases"]
    # the five-way union is the four-way union plus whatever of the fifth lies outside it
    assert mine["union_bases"] >= theirs["union_bases"]
    assert mine["sum_of_the_singles_bases"] == theirs["sum_of_the_singles_bases"] + pl.bases(ivs[p5.FIFTH])


def test_the_pairs_the_fifth_adds_are_counted_and_the_double_count_follows(gff, tracks):
    ivs = p5.covariate_intervals(p5.covariate_model(CHROM, gff), tracks)
    ov = p5.pairwise_overlap(ivs, LENGTH)
    for name in pu.COVARIATES:
        assert f"{name}_and_{p5.FIFTH}" in ov["pairs"]
    # the duplication track lies in a block, so it touches no gene-model covariate by construction
    assert ov["pairs"][f"coding_gene_intron_and_{p5.FIFTH}"]["shared_bases"] == 0
    assert ov["bases_counted_more_than_once_by_the_sum"] == (
        ov["sum_of_the_singles_bases"] - ov["union_bases"]
    )


def test_a_pair_with_no_track_is_null_and_not_zero(gff):
    ivs = p5.covariate_intervals(p5.covariate_model(CHROM, gff), {"missing": [p5.FIFTH]})
    ov = p5.pairwise_overlap(ivs, LENGTH)
    assert ov["pairs"][f"exon_any_and_{p5.FIFTH}"]["shared_bases"] is None
    assert p5.FIFTH not in ov["covariates_with_bases"]


# -- the arms: the five-way is the four-way plus one, and the arrays are the same arrays -----------
def test_the_arms_are_the_four_way_lanes_arms_plus_the_fifth(gff, tracks, tmp_path):
    tpl.write_store(tmp_path / CHROM)
    model = pl.gene_model(CHROM, gff)
    cut = pl.replay_the_cut(model["spans_in_start_order"], LENGTH)
    tpl.write_budget(tmp_path, cut["blocks"])
    tpl.write_panel(tmp_path, cut["blocks"])
    rb = pb.rebuild(CHROM, tmp_path, tmp_path)
    ivs = p5.covariate_intervals(p5.covariate_model(CHROM, gff), tracks)
    arms = p5.arm_intervals(ivs)
    mine = p5.arm_arrays(rb, arms)
    theirs, their_arms = pu.union_arrays(rb, {n: ivs[n] for n in pu.COVARIATES})
    # every array the four-way lane builds is built here, bit for bit, under the same key
    for key, arr in theirs.items():
        assert mine[key] == arr, key
    assert arms[pu.UNION] == their_arms[pu.UNION]
    assert pl.bases(arms[p5.UNION]) >= pl.bases(arms[pu.UNION])


def test_the_five_way_arm_excludes_the_four_way_arm_and_no_more_than_the_sum(out):
    five = out["offset"][f"without_{p5.UNION}"]
    four = out["offset"][f"without_{pu.UNION}"]
    singles = [out["offset"][f"without_{n}"] for n in p5.COVARIATES]
    singles = [s for s in singles if s.get("excluded_bases")]
    assert five["excluded_bases"] >= four["excluded_bases"]
    assert five["excluded_bases"] <= sum(s["excluded_bases"] for s in singles)
    assert five["covariates_included"] == list(p5.COVARIATES)
    assert four["covariates_included"] == list(pu.COVARIATES)
    assert five["widening_covariates_included"] == [p5.FIFTH]
    assert four["widening_covariates_included"] == []


# -- which way each covariate runs, read off the background rather than asserted -------------------
def test_the_fifth_reads_noisier_than_the_rest_by_construction(out):
    row = out["which_way_each_covariate_runs"][p5.FIFTH]
    assert row["recurring_per_kb_where_it_is"] > row["recurring_per_kb_where_it_is_not"]
    assert row["reads"] == p5.NOISIER_THAN_THE_REST
    assert row["declared_as_widening"] is True


def test_the_fifths_arm_widens_the_offset_on_every_tier_it_has(out):
    """Taking noisier-than-average sequence out of the background lifts every tier away from 1."""
    tiers = out["eligible"]["tiers_with_blocks_on_this_chromosome"]
    signs = {
        t: out["explained"][t][f"without_{p5.FIFTH}"]["explains"]
        for t in tiers
        if out["explained"][t].get(f"without_{p5.FIFTH}", {}).get("explains") is not None
    }
    assert signs, "the fixture must give the fifth covariate an assessed arm on some tier"
    assert all(v < 0 for v in signs.values()), signs
    for t in signs:
        assert "widens" in out["explained"][t][f"without_{p5.FIFTH}"]["reads"]


def test_the_rate_and_the_arm_are_two_readings_and_the_result_says_whether_they_agree(out):
    row = out["which_way_each_covariate_runs"][p5.FIFTH]
    assert row["the_rate_and_the_arm"] in (p5.AGREE, p5.DISAGREE, p5.UNASSESSED)
    assert row["the_rate_and_the_arm"] == p5.AGREE
    # the direction of the ratio is what the rate is checked against; the sign of the share is beside it
    assert row["tiers_whose_ratio_rises_when_it_is_removed"]
    assert row["tiers_whose_ratio_falls_when_it_is_removed"] == []
    assert row["tiers_whose_arm_widens_the_offset"]


def test_the_direction_is_taken_from_the_ratio_and_not_from_the_sign_of_the_share():
    """A tier below 1 flips the sign of every share over it, and must not flip the direction.

    The constrained-unknown tier reads below the panel's background, so removing quiet sequence pushes
    it further below and the share reads "widens" while the background moved exactly as it does under
    every other tier. Checking a rate against the share would call that a disagreement; checking it
    against the ratio does not.
    """
    exp = {
        # above 1: removing the quiet covariate lowers the ratio, and the share reads as narrowing
        "fossil": {
            "ratio_as_built": 1.10,
            "without_exon_any": {"assessed": True, "ratio": 1.05, "explains": 0.5},
        },
        # below 1: the same move in the background, and the share reads as widening
        "constrained_unknown": {
            "ratio_as_built": 0.92,
            "without_exon_any": {"assessed": True, "ratio": 0.88, "explains": -0.5},
        },
    }
    splits = {"exon_any": {"with_per_kb": 5.3, "without_per_kb": 7.4, "share_of_background_bases": 0.1}}
    row = p5.which_way_each_covariate_runs(splits, exp)["exon_any"]
    assert row["reads"] == p5.QUIETER_THAN_THE_REST
    assert row["tiers_whose_ratio_falls_when_it_is_removed"] == ["constrained_unknown", "fossil"]
    assert row["tiers_whose_ratio_rises_when_it_is_removed"] == []
    assert row["the_rate_and_the_arm"] == p5.AGREE  # one rate, one direction, two signs
    assert row["tiers_whose_arm_widens_the_offset"] == ["constrained_unknown"]
    assert row["tiers_whose_arm_narrows_the_offset"] == ["fossil"]


def test_a_covariate_with_no_kilobase_on_either_side_is_not_given_a_direction():
    assert p5.which_way_it_runs(None, 7.0) == p5.RATE_NOT_ASSESSED
    assert p5.which_way_it_runs(7.0, None) == p5.RATE_NOT_ASSESSED
    assert p5.which_way_it_runs(7.0, 0) == p5.RATE_NOT_ASSESSED
    assert p5.which_way_it_runs(7.0, 7.05) == p5.SAME_AS_THE_REST
    assert p5.which_way_it_runs(14.3, 7.4) == p5.NOISIER_THAN_THE_REST
    assert p5.which_way_it_runs(5.3, 7.4) == p5.QUIETER_THAN_THE_REST


# -- the sign, which is the thing this lane exists to keep ----------------------------------------
def _shaped(singles, four, five, degenerate=False):
    out = {
        f"without_{n}": {"assessed": True, "degenerate": False, "explains": v}
        for n, v in zip(p5.COVARIATES, singles, strict=True)
    }
    out[f"without_{pu.UNION}"] = {"assessed": True, "degenerate": False, "explains": four}
    out[f"without_{p5.UNION}"] = {
        "assessed": not degenerate,
        "degenerate": degenerate,
        "explains": five,
    }
    return out


def test_the_widening_single_is_kept_apart_from_the_narrowing_ones():
    """The two sums are separate, the signed sum is both, and no mean mixes them."""
    sit = p5.where_the_five_sit(_shaped([0.56, 0.13, 0.31, 0.01, -0.29], 0.99, 0.70))
    assert sit["narrowing_singles"] == ["coding_gene_intron", "exon_any", "conserved_elements", "promoter"]
    assert sit["widening_singles"] == [p5.FIFTH]
    assert sit["sum_of_the_narrowing_singles"] == 1.01
    assert sit["sum_of_the_widening_singles"] == -0.29
    assert sit["signed_sum_of_all_five"] == 0.72
    assert sit["the_widening_singles_are_not_averaged_in"] is True
    assert sit["per_single_with_its_sign"][p5.FIFTH] == -0.29
    # the largest single is the largest NARROWING one: a negative arm is not a small positive one
    assert sit["largest_narrowing_single"] == 0.56
    assert sit["sign_rule"] == p5.SIGN_RULE


def test_a_five_way_union_smaller_than_the_four_way_is_the_result_and_is_named():
    sit = p5.where_the_five_sit(_shaped([0.56, 0.13, 0.31, 0.01, -0.29], 0.99, 0.70))
    assert sit["four_way_union_explains"] == 0.99
    assert sit["union_explains"] == 0.70
    assert sit["the_fifth_covariate_moves_the_union_by"] == -0.29
    assert sit["the_fifth_covariate_reads"] == p5.FIFTH_LOWERS
    assert "LESS than the same union without it" in p5.WITH_A_WIDENING_TERM
    assert sit["what_a_union_with_a_widening_term_means"] == p5.WITH_A_WIDENING_TERM


def test_a_five_way_union_larger_than_the_four_way_is_named_as_the_arms_interacting():
    sit = p5.where_the_five_sit(_shaped([0.56, 0.13, 0.31, 0.01, -0.29], 0.99, 1.06))
    assert sit["the_fifth_covariate_moves_the_union_by"] == 0.07
    assert sit["the_fifth_covariate_reads"] == p5.FIFTH_RAISES


def test_a_fifth_that_moves_nothing_is_a_sentence_and_not_a_zero():
    sit = p5.where_the_five_sit(_shaped([0.56, 0.13, 0.31, 0.01, -0.29], 0.99, 0.9901))
    assert sit["the_fifth_covariate_reads"] == p5.FIFTH_MOVES_NOTHING
    assert sit["the_fifth_covariate_moves_the_union_by"] == 0.0001


def test_a_share_above_one_is_an_overshoot_in_this_lane_too():
    sit = p5.where_the_five_sit(_shaped([0.66, 0.33, 0.25, 0.06, -0.34], 1.41, 1.12))
    assert sit["union_overshoots_the_offset"] is True
    assert sit["overshoot_reads"] == pu.OVERSHOOT  # the four-way lane's own words, not a new sentence
    below = p5.where_the_five_sit(_shaped([0.3, 0.11, 0.18, 0.03, -0.33], 0.56, 0.30))
    assert below["union_overshoots_the_offset"] is False
    assert below["overshoot_reads"] == pu.NO_OVERSHOOT


def test_the_union_is_placed_against_the_narrowing_interval_and_the_signed_sum():
    sit = p5.where_the_five_sit(_shaped([0.60, 0.10, 0.10, 0.10, -0.30], 0.88, 0.62))
    assert sit["largest_narrowing_single"] == 0.60
    assert sit["sum_of_the_narrowing_singles"] == 0.90
    assert sit["reads"] == pu.NEAR_MAX
    assert sit["position_among_the_narrowing_singles"] == round((0.62 - 0.60) / 0.30, 4)
    assert sit["signed_sum_of_all_five"] == 0.60
    assert sit["union_against_the_signed_sum_of_all_five"] == 0.02
    assert p5.where_the_five_sit(_shaped([0.6, 0.1, 0.1, 0.1, -0.3], 0.88, 0.89))["reads"] == pu.NEAR_SUM
    assert p5.where_the_five_sit(_shaped([0.6, 0.1, 0.1, 0.1, -0.3], 0.88, 0.75))["reads"] == pu.BETWEEN
    assert p5.where_the_five_sit(_shaped([0.6, 0.1, 0.1, 0.1, -0.3], 0.88, 0.40))["reads"] == pu.BELOW_MAX
    assert p5.where_the_five_sit(_shaped([0.6, 0.1, 0.1, 0.1, -0.3], 0.88, 0.99))["reads"] == pu.ABOVE_SUM


def test_a_degenerate_five_way_union_reads_as_unassessed_rather_than_as_a_number():
    sit = p5.where_the_five_sit(_shaped([0.6, 0.1, 0.1, 0.1, -0.3], 0.88, None, degenerate=True))
    assert sit["union_explains"] is None
    assert sit["union_degenerate"] is True
    assert sit["reads"] == pu.UNASSESSED
    assert sit["the_fifth_covariate_reads"] == p5.NOT_COMPARABLE
    assert sit["the_fifth_covariate_moves_the_union_by"] is None


def test_a_single_the_chromosome_does_not_carry_is_named_and_not_summed():
    shaped = _shaped([0.6, 0.1, 0.1, 0.1, -0.3], 0.88, 0.62)
    shaped[f"without_{p5.FIFTH}"] = {"assessed": False, "degenerate": False, "why": p5.NO_TRACK}
    sit = p5.where_the_five_sit(shaped)
    assert sit["singles_not_in_any_sum"] == [p5.FIFTH]
    assert sit["widening_singles"] == []
    assert sit["sum_of_the_widening_singles"] is None
    assert sit["signed_sum_of_all_five"] == 0.9


# -- the four-way reading beside it is the four-way lane's, not a rewrite of it --------------------
def test_the_four_way_row_is_the_four_way_lanes_own_reading(out):
    for tier in out["where_the_union_sits"]:
        four = out["where_the_union_sits"][tier]["the_four_way_as_it_stands"]
        assert four == pu.where_it_sits(out["explained"][tier])
        # the fifth covariate is nowhere in it: its singles are the four and its union is all_four
        assert p5.FIFTH not in four["per_single"]
        assert p5.FIFTH not in four["singles_in_the_sum"]


# -- degeneracy, carried with the imported threshold ----------------------------------------------
def test_every_arm_carries_the_block_share_and_the_imported_threshold(out):
    assert p5.DEGENERATE_MIN == pl.DEGENERATE_MIN
    for key in p5.ARMS:
        arm = out["offset"][key]
        if not arm.get("assessed") and not arm.get("degenerate"):
            continue  # the covariate has no bases: named, not zeroed
        assert arm["degenerate_threshold"] == pl.DEGENERATE_MIN
        assert "block_share_of_the_remaining_background" in arm
        share = arm["block_share_of_the_remaining_background"]
        assert arm["degenerate"] == (share is not None and share >= pl.DEGENERATE_MIN)


def test_a_degenerate_arm_is_given_no_share_by_explained():
    offset = {
        "as_the_panel_builds_it": {"fossil": 1.10},
        f"without_{p5.FIFTH}": {"assessed": True, "degenerate": False, "fossil": 1.15},
        f"without_{p5.UNION}": {
            "assessed": False,
            "degenerate": True,
            "block_share_of_the_remaining_background": 0.99,
            "fossil": 1.01,
            "why": pu.DEGENERATE_WHY,
        },
    }
    exp = pb.explained(offset, "fossil")
    assert exp[f"without_{p5.FIFTH}"]["explains"] == -0.5  # the sign survives the decomposition
    assert "widens" in exp[f"without_{p5.FIFTH}"]["reads"]
    assert exp[f"without_{p5.UNION}"]["explains"] is None


# -- the windows, and the fifth counted into them -------------------------------------------------
def test_the_fifth_is_counted_into_the_four_way_lanes_own_windows(out):
    win = out["coverage"]["windows"]
    assert "windows_free_of_all_four" not in win
    assert win["windows_free_of_all_five"] + win["windows_holding_at_least_one"] == win["windows"]
    assert win["windows_holding_a_widening_covariate"] > 0
    assert win["coverage_measure"] == pb.COVERAGE_MEASURE


def test_add_the_fifth_counts_a_window_that_holds_only_the_fifth():
    arr = pb.cover_cells([(0, 100)], 4)
    rows = [
        {"start": 0, "covariates_held": 0, "covariate_bases": 0},
        {"start": 100, "covariates_held": 2, "covariate_bases": 50},
    ]
    p5.add_the_fifth(rows, {p5.FIFTH: arr})
    assert rows[0]["covariates_held"] == 1 and rows[0]["covariate_bases"] == 100
    assert rows[0]["holds_a_widening_covariate"] is True
    assert rows[1]["covariates_held"] == 2 and rows[1]["holds_a_widening_covariate"] is False


def test_presence_is_counted_before_the_stratum(out):
    row = out["comparisons"]["the_covariates_against_the_covariate_free_background"]
    presence = row["input_presence"]
    assert set(presence) == set(pu.INPUTS)
    for entry in presence.values():
        assert entry["kind"] in ("free", "bought")
        assert "targets_with_the_input" in entry
    for key in p5.STRATIFICATIONS:
        assert key in row


def test_the_control_is_the_sequence_free_of_all_five(out):
    row = out["comparisons"]["the_covariates_against_the_covariate_free_background"]
    assert row["controls"] == out["comparisons"]["control_windows_free_of_all_five"]
    assert "control_windows_free_of_all_four" not in out["comparisons"]
    for tier in out["eligible"]["tiers_with_blocks_on_this_chromosome"]:
        other = out["comparisons"][f"{tier}_against_the_covariate_free_background"]
        assert other["controls"] == row["controls"]


# -- the rest of the run's obligations ------------------------------------------------------------
def test_every_group_prints_length_gc_distance_and_coverage(out):
    for name, g in out["groups"].items():
        if not g.get("bases"):
            continue
        assert g["gc_median"] is not None, name
        assert g["tss_median"] is not None, name
        assert "coverage_median" in g, name
        assert g.get("length_median", 1000), name
    assert "background_free_of_all_five" in out["groups"]
    assert f"kilobases_holding_{p5.FIFTH}" in out["groups"]


def test_the_claim_is_one_of_the_two_and_never_a_third(out):
    assert out["claim_available"] in pu.CLAIMS
    assert "nothing explains it" not in out["claim_available"]
    assert out["sign_rule"] == p5.SIGN_RULE


def test_eligibility_is_reported_before_what_was_measured(out):
    e = out["eligible"]
    assert e["covariates_declared"] == list(p5.COVARIATES)
    assert e["covariates_declared_widening"] == [p5.FIFTH]
    assert set(e["covariates_with_bases_on_this_chromosome"]) <= set(p5.COVARIATES)
    assert e["kilobases_in_the_panel_alignment"] >= out["coverage"]["kilobases_measured"]


def test_the_rebuild_is_checked_against_the_committed_panel(out):
    assert "tiers_checked" in out["reproduces_the_panel"]


def test_the_script_distils_one_chromosome(out):
    runs = [out]
    offsets = p5s.offset_across(runs)
    sits = p5s.union_across(offsets)
    for tier in sits:
        assert "reads" in sits[tier]["the_five_way"]
        assert "reads" in sits[tier]["the_four_way_as_it_stands"]
        assert set(sits[tier]["chromosomes_assessed"]) == set(p5.ARMS)
    ov = p5s.overlap_across(runs)
    assert ov["union_bases"] == out["overlap"]["union_bases"]
    assert ov["bases_the_fifth_adds_to_the_union"] == (
        out["overlap"]["union_bases"] - out["overlap"]["the_four_way"]["union_bases"]
    )
    groups = p5s.groups_across(runs)
    assert groups["background"]["bases"] == out["groups"]["background"]["bases"]
    ways = p5s.which_way_across(runs)
    assert ways[p5.FIFTH]["chromosomes_reading"][p5.NOISIER_THAN_THE_REST] == 1
    assert ways[p5.FIFTH]["chromosomes_where_the_rate_and_the_arm_agree"] == 1
