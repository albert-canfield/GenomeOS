# SPDX-License-Identifier: AGPL-3.0-or-later
"""The never-cut background: the cut replayed, the classes, and the two kinds of nothing.

Everything here is synthetic: a store written base by base, a GFF3 written line by line, no network
and no committed result read. The chromosome is laid out so that every class exists and the direction
of every comparison is fixed by construction -- gene bodies are quiet, blocks and the sub-kilobase
inter-gene gap are noisy -- so a decomposition that comes out the other way is a defect, not weather.
"""

import gzip
import importlib.util
import json
from pathlib import Path

import pytest

from genomeos.attribution import panel_leftover as pl

_spec = importlib.util.spec_from_file_location("plc", Path("scripts/panel_leftover_classes.py"))
plc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(plc)

CHROM = "chrT"
LENGTH = 40_000
ASSEMBLIES = ["hg38", "a1", "a2", "a3"]

# Two coding genes and one lncRNA, laid out so that the cut leaves every class behind:
#   0.. 5,000   intergenic, cut into a block (>= min_size)
#   5,000..15,000   GENE1, protein coding: CDS, UTR and intron
#   15,600..25,000  GENE2, protein coding (the 600-base gap before it is under min_size)
#   25,000..30,000  intergenic, cut into a block
#   30,000..35,000  LNC1, a lncRNA: a non-coding gene body
#   35,000..40,000  intergenic, cut into a block
GENES = [
    ("GENE1", "protein_coding", 5_000, 15_000),
    ("GENE2", "protein_coding", 15_600, 25_000),
    ("LNC1", "lncRNA", 30_000, 35_000),
]
QUIET = [(5_000, 15_000), (15_600, 25_000), (30_000, 35_000)]  # the gene bodies


def write_store(d):
    """A four-assembly store: one alignment block over the whole chromosome, GC 0.40 everywhere.

    Gene bodies carry two recurring events per kilobase, everything else ten, so the never-cut
    sequence is quiet and the blocks are noisy by construction.
    """
    d.mkdir(parents=True, exist_ok=True)
    with gzip.open(d / "blocks.tsv.gz", "wt") as fh:
        fh.write(f"0\t{LENGTH}\t{'s' * len(ASSEMBLIES)}\t0\n")
    with gzip.open(d / "sites.tsv.gz", "wt") as fh:
        for kb in range(LENGTH // 1000):
            start = kb * 1000
            n = 2 if any(s <= start < e for s, e in QUIET) else 10
            for i in range(n):
                fh.write(f"{start + i * (1000 // n) + 3}\tA\t=G=G\n")
    with gzip.open(d / "insertions.tsv.gz", "wt") as fh:
        fh.write("")
    with gzip.open(d / "gc.tsv.gz", "wt") as fh:
        for b in range(LENGTH // 100):
            fh.write(f"{b}\t40\t100\n")
    (d / "meta.json").write_text(json.dumps({"chrom": CHROM, "assemblies": ASSEMBLIES, "source": "test"}))
    (d / "replication_timing.json").write_text(json.dumps({str(b): 1.0 for b in range(LENGTH // 1000)}))


def write_gff(path):
    """GENCODE-dialect GFF3, 1-based inclusive, with exons and CDS on the two coding genes."""
    rows = []
    for name, gtype, s, e in GENES:
        rows.append(
            (
                "gene",
                s + 1,
                e,
                f"ID={name};gene_name={name};gene_type={gtype}",
            )
        )
        rows.append(
            (
                "transcript",
                s + 1,
                e,
                f"ID={name}.1;Parent={name};transcript_type={gtype};tag=Ensembl_canonical",
            )
        )
        if gtype != "protein_coding":
            continue
        # two exons of 500 bases at each end; the CDS is the inner 300 of each, so both UTR and
        # intron exist and are not the same bases
        for a, b in ((s, s + 500), (e - 500, e)):
            rows.append(("exon", a + 1, b, f"ID={name}.e{a};Parent={name}.1"))
            rows.append(("CDS", a + 101, b - 100, f"ID={name}.c{a};Parent={name}.1;phase=0"))
    with gzip.open(path, "wt") as fh:
        fh.write("##gff-version 3\n")
        for kind, s, e, attrs in rows:
            fh.write(f"{CHROM}\tTEST\t{kind}\t{s}\t{e}\t.\t+\t.\t{attrs}\n")


def write_budget(results_dir, blocks):
    (results_dir / f"budget_{CHROM}.json").write_text(
        json.dumps({"chromosome_length": LENGTH, "blocks": [{"start": s, "end": e} for s, e in blocks]})
    )


def write_panel(results_dir, blocks):
    tiers = ["regulatory", "neutral", "fossil"]
    rows = [
        {"start": s, "end": e, "tier": tiers[i % len(tiers)], "mammal_fraction": 0.01}
        for i, (s, e) in enumerate(blocks)
    ]
    (results_dir / f"human_panel_{CHROM}.json").write_text(json.dumps({"blocks": rows, "pooled": {}}))


@pytest.fixture()
def gff(tmp_path):
    p = tmp_path / "synthetic.gff3.gz"
    write_gff(p)
    return p


@pytest.fixture()
def out(tmp_path, gff):
    write_store(tmp_path / CHROM)
    model = pl.gene_model(CHROM, gff)
    cut = pl.replay_the_cut(model["spans_in_start_order"], LENGTH)
    write_budget(tmp_path, cut["blocks"])
    write_panel(tmp_path, cut["blocks"])
    return pl.analyse(CHROM, cache=tmp_path, results_dir=tmp_path, gff=gff)


# -- interval arithmetic -------------------------------------------------------------------------
def test_subtract_removes_only_the_overlap():
    assert pl.subtract([(0, 100)], [(20, 40)]) == [(0, 20), (40, 100)]
    assert pl.subtract([(0, 100)], [(0, 100)]) == []
    assert pl.subtract([(0, 100)], [(200, 300)]) == [(0, 100)]
    assert pl.subtract([(0, 100), (200, 300)], [(50, 250)]) == [(0, 50), (250, 300)]
    # overlapping inputs are merged first, so no base is removed twice or counted twice
    assert pl.bases([(0, 100), (50, 150)]) == 150
    assert pl.subtract([(0, 100), (50, 150)], [(10, 20), (15, 30)]) == [(0, 10), (30, 150)]


# -- the cut, replayed ---------------------------------------------------------------------------
def test_replay_partitions_the_chromosome(gff):
    model = pl.gene_model(CHROM, gff)
    cut = pl.replay_the_cut(model["spans_in_start_order"], LENGTH)
    total = (
        pl.bases(cut["blocks"])
        + pl.bases(cut["gene_spans"])
        + pl.bases(cut["intergenic_gap_under_min_size"])
        + pl.bases(cut["intergenic_tail_under_min_size"])
    )
    assert total == LENGTH  # exhaustive and disjoint: every base has exactly one reason
    assert cut["blocks"] == [(0, 5_000), (25_000, 30_000), (35_000, 40_000)]
    # the 600-base gap between GENE1 and GENE2 is declined for its length, not for being in a gene
    assert cut["intergenic_gap_under_min_size"] == [(15_000, 15_600)]
    assert cut["intergenic_tail_under_min_size"] == []


def test_replay_is_the_same_loop_as_unknown_blocks(gff):
    from genomeos.genome.annotation import Annotation
    from genomeos.genome.unknown import unknown_blocks

    ann = Annotation.from_gff3(gff, {CHROM})
    theirs = sorted((b.start, b.end) for b in unknown_blocks(ann, CHROM, LENGTH, pl.MIN_SIZE))
    model = pl.gene_model(CHROM, gff)
    assert pl.replay_the_cut(model["spans_in_start_order"], LENGTH)["blocks"] == theirs


def test_reproduces_the_budget_reports_a_disagreement():
    cut = {"blocks": [(0, 100), (200, 300)]}
    same = [{"start": 0, "end": 100}, {"start": 200, "end": 300}]
    assert pl.reproduces_the_budget(cut, same)["identical"] is True
    other = [{"start": 0, "end": 100}, {"start": 200, "end": 400}]
    r = pl.reproduces_the_budget(cut, other)
    assert r["identical"] is False
    assert r["only_in_the_replay"] == 1 and r["only_in_the_committed_budget"] == 1


# -- the classes ---------------------------------------------------------------------------------
def test_the_classes_partition_the_never_cut_sequence(gff):
    model = pl.gene_model(CHROM, gff)
    cut = pl.replay_the_cut(model["spans_in_start_order"], LENGTH)
    cls = pl.classify(model, cut, LENGTH)
    ivs = cls["intervals"]
    assert sum(pl.bases(v) for v in ivs.values()) == LENGTH - pl.bases(cut["blocks"])
    assert cls["bases_in_a_block_and_in_a_class"] == 0
    assert ivs[pl.UNACCOUNTED] == []
    # CDS, UTR and intron are three different sets of bases, and they add up to the gene bodies
    assert pl.bases(ivs["coding_gene_cds"]) == 2 * 2 * 300
    assert pl.bases(ivs["coding_gene_utr"]) == 2 * 2 * 200
    assert pl.bases(ivs["coding_gene_intron"]) == (10_000 - 1_000) + (9_400 - 1_000)
    assert pl.bases(ivs["noncoding_gene_body"]) == 5_000
    assert pl.bases(ivs["pseudogene_body"]) == 0  # measured, and this chromosome has none


def test_every_class_carries_the_reason_the_budget_declined_it():
    assert set(pl.WHY_NOT_CUT) == {*pl.CLASSES, pl.UNACCOUNTED}
    assert "min_size" in pl.WHY_NOT_CUT[pl.LENGTH_RULE_ONLY]
    assert "gene span" in pl.WHY_NOT_CUT["coding_gene_intron"]


# -- separability: a property of the rule, not of the sample --------------------------------------
def test_separability_counts_only_the_length_rule():
    split = {
        "by_class": {
            "coding_gene_intron": {"bases": 990},
            "intergenic_gap_under_min_size": {"bases": 10},
            "intergenic_tail_under_min_size": {"bases": 0},
        }
    }
    s = pl.separability(split)
    assert s["bases_declined_for_the_length_rule_alone"] == 10
    assert s["share_that_can_separate_the_two_readings"] == 0.01
    assert s["separable"] is False
    split["by_class"]["intergenic_gap_under_min_size"]["bases"] = 500
    assert pl.separability(split)["separable"] is True


# -- the verdicts, which are data and not prose ---------------------------------------------------
def test_verdict_needs_a_gap_before_it_says_anything():
    assert pl.verdict(True, None, 0.0, 10_000)["verdict"] == pl.VERDICT_UNASSESSED
    assert pl.verdict(True, -0.01, 0.0, 10_000)["verdict"] == pl.VERDICT_NEITHER


def test_verdict_is_blocked_when_the_two_readings_are_the_same_predicate():
    out = pl.verdict(False, 0.02, 0.06, 10_000)
    assert out["verdict"] == pl.VERDICT_SAME_PARTITION
    # the mass question is blocked and the mechanism question is not: they are separate sentences
    assert out["mechanism"] == pl.MECHANISM_GENE_BODY


def test_verdict_separates_the_cut_from_the_class_when_it_can():
    assert pl.verdict(True, 0.02, 0.001, 10_000)["verdict"] == pl.VERDICT_INSTRUMENT
    assert pl.verdict(True, 0.02, 0.05, 10_000)["verdict"] == pl.VERDICT_CLASS


def test_mechanism_is_unassessed_below_its_own_sample_size():
    out = pl.verdict(False, 0.02, 0.06, pl.MECHANISM_MIN_WINDOWS - 1)
    assert out["mechanism"] == pl.MECHANISM_UNASSESSED
    assert out["length_rule_windows"] == pl.MECHANISM_MIN_WINDOWS - 1
    assert out["mechanism_on_the_line"] is False


def test_a_mechanism_decided_by_a_hair_says_so_in_the_result():
    # a hand-picked threshold can decide a sentence by a hair, and the margin is data, not prose
    clear = pl.verdict(False, 0.02, 0.06, 10_000)
    assert clear["mechanism_margin"] == 0.05 and clear["mechanism_on_the_line"] is False
    hair = pl.verdict(False, 0.02, 0.0101, 10_000)
    assert hair["mechanism"] == pl.MECHANISM_GENE_BODY
    assert hair["mechanism_on_the_line"] is True
    assert "would flip if the threshold moved" in hair["mechanism_caveat"]
    # and the same number as a position on the line from the never-cut sequence to the blocks
    assert hair["length_rule_share_of_the_block_gap"] == 0.505


# -- end to end ----------------------------------------------------------------------------------
def test_the_cut_is_checked_before_anything_is_read_from_it(out):
    assert out["assessed"] is True
    assert out["reproduces_the_budget"]["identical"] is True
    assert out["reproduces_the_budget"]["blocks_replayed"] == 3


def test_what_was_eligible_is_reported_before_what_was_measured(out):
    e = out["eligible"]
    assert e["classes_declared"] == [*pl.CLASSES, pl.UNACCOUNTED]
    # declared, and the ones this chromosome actually carries: different lists, both printed
    assert "pseudogene_body" not in e["classes_with_bases_on_this_chromosome"]
    assert "coding_gene_intron" in e["classes_with_bases_on_this_chromosome"]
    assert e["kilobases_in_the_panel_alignment"] == LENGTH // 1000


def test_the_split_accounts_for_every_never_cut_base(out):
    s = out["split"]
    assert s["bases_in_no_class"] == 0
    assert s["by_class"][pl.UNACCOUNTED]["bases"] == 0
    assert s["separability"]["separable"] is False
    assert s["separability"]["bases_declined_for_the_length_rule_alone"] == 600


def test_a_class_with_no_bases_is_not_the_same_as_a_class_that_is_quiet(out):
    # measured, and this chromosome has none: a zero with a denominator beside it
    assert out["split"]["by_class"]["pseudogene_body"]["bases"] == 0
    assert out["coverage"]["classes_not_assessed"]["pseudogene_body"] == (
        "read, and empty on this chromosome"
    )
    # never looked: the arm reports itself unassessed with the reason, rather than a ratio
    assert out["offset"]["without_pseudogene_body"]["assessed"] is False
    assert out["explained"]["neutral"]["without_pseudogene_body"]["reads"] == "not assessed"


def test_every_group_prints_its_covariates_and_its_coverage(out):
    for g in (out["background"], out["leftover"]):
        for key in ("bases", "gc_median", "tss_median", "coverage_median", "recurring_per_kb"):
            assert g[key] is not None
    # the never-cut sequence is the quiet half by construction, and the measurement says so
    assert out["leftover"]["recurring_per_kb"] < out["background"]["recurring_per_kb"]


def test_the_degenerate_arm_names_itself_rather_than_being_counted(out):
    arm = out["offset"]["without_every_never_cut_base"]
    assert arm["degenerate"] is True
    assert arm["assessed"] is False
    assert arm["block_share_of_the_remaining_background"] >= pl.DEGENERATE_MIN
    assert "tier average" in arm["why"]


def test_the_block_comparison_runs_on_two_groups_of_sequence(out):
    c = out["comparisons"]["every_block_against_the_never_cut"]
    assert c["assessed"] is True
    # presence of every input is counted before any coverage stratum is applied to it
    for row in c["input_presence"].values():
        assert row["kind"] in ("free", "bought")
    d = c["gc_timing_tss_and_coverage"]["matched"]["difference"]
    assert d > 0  # blocks are noisy and the never-cut sequence is quiet, by construction
    # targets with no control in their stratum are counted rather than quietly kept
    dropped = c["gc_timing_tss_and_coverage"]["dropped_for_want_of_a_control"]
    assert dropped + c["gc_timing_tss_and_coverage"]["targets_matched"] == c["targets"]
    # this chromosome is 40 kb and its strata are coarse, so some targets have no control at all;
    # the point of the assertion is that the number is reported, not that it is nil
    assert dropped < 0.5 * c["targets"]


def test_the_length_rule_class_is_compared_against_the_rest_of_the_never_cut(out):
    key = f"{pl.LENGTH_RULE_ONLY}_against_the_rest_of_the_never_cut"
    c = out["comparisons"][key]
    assert c["assessed"] is True
    # the gap is noisy like a block although the budget declined it: the class, not the cut
    assert c["gc_timing_tss_and_coverage"]["matched"]["difference"] > 0


def test_the_claim_is_one_of_two_strings_and_never_a_third(out):
    assert out["claim_available"] in (pl.CLAIM_WITH_COVERAGE, pl.CLAIM_WITHOUT_COVERAGE)
    assert pl.CLAIM_WITHOUT_COVERAGE == "the covariates I measured do not explain it"


# -- the distillation, which must not need the stores again ---------------------------------------
def test_the_top_level_can_be_redistilled_without_reading_a_store(tmp_path, out):
    """Reading the stores is the expensive half; correcting the distillation must not repeat it."""
    from genomeos.results import save_result

    result = {"per_chromosome": {out["chrom"]: {k: out[k] for k in plc.PER_CHROMOSOME}}}
    save_result("probe", result, tmp_path)
    runs = plc.from_result("probe", tmp_path)
    assert len(runs) == 1
    assert runs[0]["chrom"] == CHROM and runs[0]["assessed"] is True
    # every field the across-chromosome aggregators read is carried, so nothing is re-measured
    assert plc.split_across(runs)["never_cut_bases"] == out["split"]["bases_accounted_by_a_class"]
    assert (
        plc.comparisons_across(runs)[plc.BLOCK_ARM]["targets"]
        == (out["comparisons"][plc.BLOCK_ARM]["targets"])
    )
    assert plc.offset_across(runs)["neutral"]["as_the_panel_builds_it"]["chromosomes_assessed"] == 1


def test_redistilling_a_result_that_has_no_runs_says_so(tmp_path):
    from genomeos.results import save_result

    save_result("probe", {"chromosomes": []}, tmp_path)
    with pytest.raises(FileNotFoundError, match="no per-chromosome runs"):
        plc.from_result("probe", tmp_path)
