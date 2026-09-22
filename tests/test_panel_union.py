# SPDX-License-Identifier: AGPL-3.0-or-later
"""The union arm: the four covariates taken out at once, and the two things that could make it a lie.

Everything here is synthetic -- a store written base by base, a GFF3 and a conserved-element track
written line by line, no network and no committed result read. The fixtures are the ones
`test_panel_leftover.py` already lays out, loaded from that file rather than copied, so that a
chromosome whose classes are fixed by construction there is the same chromosome here.

The two things that could make the union arm a lie are the two things tested hardest: a covariate
model that has quietly drifted from the lanes whose definitions it claims to use, and a degeneracy
flag that does not fire when the exclusion leaves a background made of blocks.
"""

import gzip
import importlib.util
from pathlib import Path

import pytest

from genomeos.attribution import panel_background as pb
from genomeos.attribution import panel_leftover as pl
from genomeos.attribution import panel_union as pu

_spec = importlib.util.spec_from_file_location("tpl", Path("tests/test_panel_leftover.py"))
tpl = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(tpl)

_script = importlib.util.spec_from_file_location("pua", Path("scripts/panel_union_arms.py"))
pua = importlib.util.module_from_spec(_script)
_script.loader.exec_module(pua)

CHROM = tpl.CHROM
LENGTH = tpl.LENGTH
# conserved elements: one inside GENE1's intron, one in the block at 25,000..30,000, so the track
# overlaps the introns and the blocks both, which is the arrangement the overlap table has to show
CONSERVED = [(7_000, 7_500), (26_000, 26_400)]


def write_conserved(cache: Path):
    cache.mkdir(parents=True, exist_ok=True)
    with gzip.open(cache / f"phastConsElements100way_{CHROM}.bed.gz", "wt") as fh:
        for s, e in CONSERVED:
            fh.write(f"{CHROM}\t{s}\t{e}\tel\t500\n")


@pytest.fixture()
def gff(tmp_path):
    p = tmp_path / "synthetic.gff3.gz"
    tpl.write_gff(p)
    return p


@pytest.fixture()
def out(tmp_path, gff):
    tpl.write_store(tmp_path / CHROM)
    model = pl.gene_model(CHROM, gff)
    cut = pl.replay_the_cut(model["spans_in_start_order"], LENGTH)
    tpl.write_budget(tmp_path, cut["blocks"])
    tpl.write_panel(tmp_path, cut["blocks"])
    write_conserved(tmp_path / "constraint")
    return pu.analyse(
        CHROM,
        cache=tmp_path,
        results_dir=tmp_path,
        gff=gff,
        tracks_dir=tmp_path,
        constraint_dir=tmp_path / "constraint",
    )


# -- the covariate model is the two lanes' own definitions, not a third one ----------------------
def test_the_gene_model_equals_what_both_lanes_compute(gff, monkeypatch):
    """One walk of the annotation, and it has to give back exactly what the two lanes' walks give.

    `panel_background.gene_features` takes no `gff` and asks `default_gencode` for one, so the
    synthetic annotation is put where it looks. That is the whole point of the test: if either lane
    changes what it means by `exon_any` or by a coding TSS, this fails rather than drifting.
    """
    monkeypatch.setattr("genomeos.genome.annotation.default_gencode", lambda chroms=None: gff)
    mine = pu.covariate_model(CHROM, gff)
    theirs_bg = pb.gene_features(CHROM)
    theirs_pl = pl.gene_model(CHROM, gff)
    assert mine["exon_any"] == theirs_bg["exon_any"]
    assert mine["tss"] == theirs_bg["tss"] == theirs_pl["tss"]
    assert mine["coding_spans"] == theirs_pl["coding_spans"]
    assert mine["coding_exons"] == theirs_pl["coding_exons"]
    assert mine["genes"] == theirs_pl["genes"]


def test_the_intron_covariate_is_panel_leftovers_class(gff):
    """`coding_gene_intron` here must be the same bases `panel_leftover.classify` calls by that name."""
    model = pl.gene_model(CHROM, gff)
    cut = pl.replay_the_cut(model["spans_in_start_order"], LENGTH)
    theirs = pl.classify(model, cut, LENGTH)["intervals"]["coding_gene_intron"]
    mine = pu.covariate_intervals(pu.covariate_model(CHROM, gff), {})["coding_gene_intron"]
    assert mine == theirs


def test_promoter_and_conserved_come_from_the_named_places(gff, tmp_path):
    write_conserved(tmp_path / "constraint")
    tracks = pb.local_tracks(CHROM, tmp_path, tmp_path / "constraint")
    ivs = pu.covariate_intervals(pu.covariate_model(CHROM, gff), tracks)
    assert ivs["conserved_elements"] == CONSERVED
    # the same intervals panel_background._ivs builds for "promoter", from the same TSS
    genes = {"tss": pl.gene_model(CHROM, gff)["tss"]}
    assert ivs["promoter"] == pb._ivs({}, genes, "promoter")


def test_a_missing_track_and_an_empty_one_are_different_answers(gff):
    model = pu.covariate_model(CHROM, gff)
    assert pu.why_absent("conserved_elements", model, {"missing": ["conserved_elements"]}) == pu.NO_TRACK
    assert pu.why_absent("conserved_elements", model, {"missing": []}) == pu.MEASURED_EMPTY
    assert pu.why_absent("promoter", {}, {"missing": []}) == pu.NO_GENCODE


# -- the overlaps, which decide where the union can sit ------------------------------------------
def test_overlap_counts_shared_bases_and_the_double_count(gff, tmp_path):
    write_conserved(tmp_path / "constraint")
    tracks = pb.local_tracks(CHROM, tmp_path, tmp_path / "constraint")
    ivs = pu.covariate_intervals(pu.covariate_model(CHROM, gff), tracks)
    ov = pu.pairwise_overlap(ivs, LENGTH)
    # the element at 7,000..7,500 is inside GENE1's intron; the one at 26,000 is in a block
    assert ov["pairs"]["coding_gene_intron_and_conserved_elements"]["shared_bases"] == 500
    # an exon and an intron of the same gene share no base, by construction of the two sets
    assert ov["pairs"]["coding_gene_intron_and_exon_any"]["shared_bases"] == 0
    assert ov["union_bases"] < ov["sum_of_the_singles_bases"]
    assert ov["bases_counted_more_than_once_by_the_sum"] == ov["sum_of_the_singles_bases"] - ov["union_bases"]


def test_a_pair_with_no_bases_is_null_and_not_zero(gff):
    ivs = pu.covariate_intervals(pu.covariate_model(CHROM, gff), {})
    ov = pu.pairwise_overlap(ivs, LENGTH)
    row = ov["pairs"]["coding_gene_intron_and_conserved_elements"]
    assert row["shared_bases"] is None  # never looked: the track is not there
    assert "conserved_elements" not in ov["covariates_with_bases"]


# -- degeneracy, which is the live risk ----------------------------------------------------------
def test_block_share_after_reads_what_the_arm_leaves():
    """Two kilobases, one all block and one all background; an arm that takes the background leaves blocks."""
    rows = [
        {pu.BLOCKS: 1000, f"{pu.BLOCKS}_and_a": 0, "a": 0},
        {pu.BLOCKS: 0, f"{pu.BLOCKS}_and_a": 0, "a": 1000},
    ]
    assert pu.block_share_after(rows, "a") == 1.0  # the remainder is entirely block sequence
    rows[1]["a"] = 0
    assert pu.block_share_after(rows, "a") == 0.5  # nothing taken: half of it is block
    # an arm that takes the blocks as well leaves neither side over-represented
    rows = [
        {pu.BLOCKS: 1000, f"{pu.BLOCKS}_and_a": 1000, "a": 1000},
        {pu.BLOCKS: 0, f"{pu.BLOCKS}_and_a": 0, "a": 0},
    ]
    assert pu.block_share_after(rows, "a") == 0.0


def test_a_degenerate_arm_is_flagged_and_given_no_share():
    """The flag fires at the line the earlier lanes drew, and `explained` then refuses the arm."""
    offset = {
        "as_the_panel_builds_it": {"fossil": 1.10},
        "without_a": {"assessed": True, "degenerate": False, "fossil": 1.05},
        "without_b": {
            "assessed": False,
            "degenerate": True,
            "block_share_of_the_remaining_background": 0.998,
            "fossil": 1.01,
            "why": pu.DEGENERATE_WHY,
        },
    }
    exp = pb.explained(offset, "fossil")
    assert exp["without_a"]["explains"] == 0.5
    assert exp["without_b"]["explains"] is None  # degenerate: reported, never counted
    assert exp["without_b"]["assessed"] is False
    assert pu.DEGENERATE_MIN == pl.DEGENERATE_MIN  # one line, drawn once, in the earlier lane


# -- where the union sits ------------------------------------------------------------------------
def test_where_it_sits_names_the_two_ends_and_the_middle():
    def exp(singles, union, degenerate=False):
        out = {
            f"without_{n}": {"assessed": True, "degenerate": False, "explains": v}
            for n, v in zip(pu.COVARIATES, singles, strict=True)
        }
        out[f"without_{pu.UNION}"] = {
            "assessed": not degenerate,
            "degenerate": degenerate,
            "explains": union,
        }
        return out

    nested = pu.where_it_sits(exp([0.60, 0.10, 0.10, 0.10], 0.62))
    assert nested["reads"] == pu.NEAR_MAX
    assert nested["largest_single"] == 0.60 and nested["sum_of_the_singles"] == 0.90
    assert pu.where_it_sits(exp([0.60, 0.10, 0.10, 0.10], 0.89))["reads"] == pu.NEAR_SUM
    assert pu.where_it_sits(exp([0.60, 0.10, 0.10, 0.10], 0.75))["reads"] == pu.BETWEEN
    assert pu.where_it_sits(exp([0.60, 0.10, 0.10, 0.10], 0.50))["reads"] == pu.BELOW_MAX
    assert pu.where_it_sits(exp([0.60, 0.10, 0.10, 0.10], 0.95))["reads"] == pu.ABOVE_SUM
    # a share at or under 1 does not overshoot; the flag is on the share, not on the reading
    assert nested["union_overshoots_the_offset"] is False
    assert nested["overshoot_reads"] == pu.NO_OVERSHOOT


def test_a_share_above_one_is_named_an_overshoot_and_not_an_explanation():
    """141% of an offset is a ratio carried past 1, and the result has to say so in its own words."""
    over = pu.where_it_sits(
        {
            "without_coding_gene_intron": {"assessed": True, "degenerate": False, "explains": 0.66},
            f"without_{pu.UNION}": {"assessed": True, "degenerate": False, "explains": 1.41},
        }
    )
    assert over["union_overshoots_the_offset"] is True
    assert over["overshoot_reads"] == pu.OVERSHOOT
    assert "not explanation" in pu.OVERSHOOT
    # and the flag never claims an overshoot where there is no number to claim it from
    none = pu.where_it_sits({f"without_{pu.UNION}": {"assessed": False, "degenerate": True}})
    assert none["union_overshoots_the_offset"] is False
    assert none["overshoot_reads"] == pu.UNASSESSED


def test_a_degenerate_union_reads_as_unassessed_rather_than_as_a_number():
    out = pu.where_it_sits(
        {
            "without_coding_gene_intron": {"assessed": True, "degenerate": False, "explains": 0.5},
            f"without_{pu.UNION}": {"assessed": False, "degenerate": True, "explains": None},
        }
    )
    assert out["reads"] == pu.UNASSESSED
    assert out["union_explains"] is None
    assert out["union_degenerate"] is True
    assert out["position_between_the_largest_and_the_sum"] is None


def test_a_single_left_out_of_the_sum_is_named(gff):
    out = pu.where_it_sits(
        {
            "without_coding_gene_intron": {"assessed": True, "degenerate": False, "explains": 0.5},
            "without_exon_any": {"assessed": False, "degenerate": False, "why": pu.NO_TRACK},
            f"without_{pu.UNION}": {"assessed": True, "degenerate": False, "explains": 0.6},
        }
    )
    assert out["singles_in_the_sum"] == ["coding_gene_intron"]
    assert "exon_any" in out["singles_not_in_the_sum"]


# -- the whole run -------------------------------------------------------------------------------
def test_the_union_arm_takes_the_union_of_the_four(out):
    union = out["offset"][f"without_{pu.UNION}"]
    singles = [
        out["offset"][f"without_{n}"]
        for n in pu.COVARIATES
        if out["offset"][f"without_{n}"].get("excluded_bases")
    ]
    # the union arm excludes at least as much as any single and no more than their sum
    assert union["excluded_bases"] >= max(s["excluded_bases"] for s in singles)
    assert union["excluded_bases"] <= sum(s["excluded_bases"] for s in singles)
    assert union["covariates_included"] == list(pu.COVARIATES)
    assert union["covariates_left_out"] == {}


def test_every_group_prints_the_four_covariates_it_is_compared_on(out):
    for name, g in out["groups"].items():
        if not g.get("bases"):
            continue
        assert g["gc_median"] is not None, name
        assert g["tss_median"] is not None, name
        assert "coverage_median" in g, name
        assert g.get("length_median", 1000), name


def test_presence_is_counted_before_the_stratum(out):
    row = out["comparisons"]["the_covariates_against_the_covariate_free_background"]
    presence = row["input_presence"]
    assert set(presence) == set(pu.INPUTS)
    for entry in presence.values():
        assert entry["kind"] in ("free", "bought")
        assert "targets_with_the_input" in entry
    # the stratified numbers come after, and the ladder is the earlier lanes' ladder
    for key in pu.STRATIFICATIONS:
        assert key in row


def test_the_two_arms_partition_the_non_block_sequence_and_share_no_window(out):
    """The control is what the union arm leaves, and the target is what it takes: no window is both.

    The fixture does not fix the SIGN of this contrast -- its covariate-free remainder is the lncRNA
    body, which it also makes quiet -- so what is tested is the invariant it does fix: the two arms
    are complementary, no tier appears in either, and no window is counted twice.
    """
    row = out["comparisons"]["the_covariates_against_the_covariate_free_background"]
    assert row["assessed"] is True
    assert row["targets"] == out["comparisons"]["background_windows_holding_at_least_one"]
    assert row["controls"] == out["comparisons"]["control_windows_free_of_all_four"]
    win = out["coverage"]["windows"]
    assert row["targets"] and row["controls"]
    # every window in either arm is outside every block and outside the canonical CDS, so the two
    # together can never exceed the windows the run measured minus the ones it put in a block
    assert row["targets"] + row["controls"] <= win["windows"] - win["windows_in_a_block"]
    # and no tier is in either arm: the tier comparisons below are the only place a tier appears
    for tier in out["eligible"]["tiers_with_blocks_on_this_chromosome"]:
        assert (
            out["comparisons"][f"{tier}_against_the_covariate_free_background"]["controls"]
            == (row["controls"])
        )


def test_the_claim_is_one_of_the_two_and_never_a_third(out):
    assert out["claim_available"] in pu.CLAIMS
    assert "nothing explains it" not in out["claim_available"]


def test_the_rebuild_is_checked_against_the_committed_panel(out):
    assert "reproduces_the_panel" in out
    assert "tiers_checked" in out["reproduces_the_panel"]


def test_eligibility_is_reported_before_what_was_measured(out):
    e = out["eligible"]
    assert e["covariates_declared"] == list(pu.COVARIATES)
    assert set(e["covariates_with_bases_on_this_chromosome"]) <= set(pu.COVARIATES)
    assert e["kilobases_in_the_panel_alignment"] >= out["coverage"]["kilobases_measured"]


def test_the_script_distils_one_chromosome(out, tmp_path):
    runs = [out]
    offsets = pua.offset_across(runs)
    sits = pua.union_across(offsets)
    for tier in sits:
        assert "reads" in sits[tier]
        assert "chromosomes_assessed" in sits[tier]
    ov = pua.overlap_across(runs)
    assert ov["union_bases"] == out["overlap"]["union_bases"]
    groups = pua.groups_across(runs)
    assert groups["background"]["bases"] == out["groups"]["background"]["bases"]
