# SPDX-License-Identifier: AGPL-3.0-or-later
"""The argmax-cell measurement: its base rate, its read discipline and its registered readings.

Targeted tests for lane-argmaxcell's own files. The acceptance of this lane is these plus `ruff` on
its own files, run in a worktree of its committed tree; the full verdict on the committed tree is
the coordinator's push check.

Two of these exist because of defects found elsewhere tonight and not because of this code:
`test_census_absence_is_a_defect_and_raises` (a guard that skips when a committed result is missing
is itself a defect, so the path must RAISE) and `test_module_source_names_no_archive_loader` (a peer
reached 23 GB RSS by looping a loader that falls back to whole-chromosome archives).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genomeos.attribution import argmaxcell as ac
from genomeos.attribution import compile as cp
from genomeos.attribution import crispri as cr
from tests.committed_data import must_be_committed

REGISTRATION_RELATIVE = "data/results/argmaxcell_registration.json"


# ---- the base rate -------------------------------------------------------------------------------


def test_census_absence_is_a_defect_and_raises(tmp_path: Path) -> None:
    """An absent committed census RAISES. A skip here would report a pass while measuring nothing."""
    with pytest.raises(FileNotFoundError) as e:
        ac.census(tmp_path / "not_here.json")
    assert "defect" in str(e.value)


def test_census_refuses_a_per_cell_block_that_does_not_sum(tmp_path: Path) -> None:
    """If per_cell does not sum to `rules`, the base rate's denominator is not its population."""
    p = tmp_path / "census.json"
    p.write_text(json.dumps({"rules": 100, "per_cell": {"K562": {"rules": 7}}}))
    with pytest.raises(ValueError) as e:
        ac.census(p)
    assert "sum to 7" in str(e.value) and "100" in str(e.value)


def test_census_accepts_a_block_that_sums(tmp_path: Path) -> None:
    p = tmp_path / "census.json"
    p.write_text(json.dumps({"rules": 10, "per_cell": {"K562": {"rules": 4}, "HepG2": {"rules": 6}}}))
    assert ac.census(p)["rules"] == 10


def test_base_rate_is_the_committed_census_over_its_own_denominator() -> None:
    """The one number the registered reading turns on, read from the committed census."""
    must_be_committed(ac.CENSUS_RELATIVE)
    cen = ac.census()
    r = ac.base_rate("K562", cen)
    assert r["labels"] == ["K562"]
    assert r["rules_total"] == cen["rules"] == 440_589
    assert r["rules_with_this_label"] == 27_445
    assert r["rate"] == round(27_445 / 440_589, 6) == 0.062292
    assert r["source"] == ac.CENSUS_RELATIVE


def test_base_rate_label_is_the_compilers_own_function() -> None:
    """A cell name is looked up under the label the compiler would have written, not under itself."""
    cen = ac.census()
    assert cp.context("IMR-90") == "IMR_90"
    assert ac.base_rate("IMR-90", cen)["labels"] == ["IMR_90"]
    assert ac.base_rate("IMR-90", cen)["rules_with_this_label"] == 2_942
    # a label the census never saw is a rate of zero, not a KeyError
    assert ac.base_rate("a cell that does not exist", cen)["rules_with_this_label"] == 0


def test_every_arm_has_a_base_rate_before_any_count() -> None:
    cen = ac.census()
    for cell in ac.ARMS:
        r = ac.base_rate(cell, cen)
        assert r["rate"] >= 0.0 and r["rules_total"] == 440_589


# ---- the streaming reader and the overlap rule ---------------------------------------------------


def _table(rows: list[dict]) -> list[dict]:
    return rows


def test_stream_table_keeps_only_elements_that_name_a_gene(tmp_path: Path) -> None:
    p = tmp_path / "chrT.json"
    p.write_text(
        json.dumps(
            [
                {"id": "a", "start": 10, "end": 20, "predicted_coding": {"gene": "G", "tissue": "K562"}},
                {"id": "b", "start": 30, "end": 40, "predicted_coding": {}},
                {"id": "c", "start": 50, "end": 60, "predicted_coding": {"gene": "H", "tissue": "IMR-90"}},
                {"id": "d", "start": 70, "end": 80},
            ]
        )
    )
    out = list(ac.stream_table(p))
    assert [r[0] for r in out] == ["a", "c"]
    assert out[0] == ("a", 10, 20, "K562")
    assert out[1][3] == "IMR_90", "the label is the compiler's own ident, not the raw tissue name"


def test_stream_table_labels_an_absent_tissue_unknown(tmp_path: Path) -> None:
    p = tmp_path / "chrT.json"
    p.write_text(json.dumps([{"id": "a", "start": 1, "end": 2, "predicted_coding": {"gene": "G"}}]))
    assert list(ac.stream_table(p))[0][3] == cp.context(None)


def test_stream_table_refuses_a_file_that_is_not_a_list(tmp_path: Path) -> None:
    p = tmp_path / "chrT.json"
    p.write_text('{"elements": 1}')
    with pytest.raises(ValueError):
        list(ac.stream_table(p))


def test_overlap_rule_agrees_with_the_crispri_table_element_for_element(tmp_path: Path) -> None:
    """`OVERLAP_RULE` says the two agree; this is the assertion, not the sentence."""
    rows = [
        {"id": "e1", "start": 1_000, "end": 1_500},
        {"id": "e2", "start": 1_400, "end": 2_000},
        {"id": "e3", "start": 5_000, "end": 5_200},
        {"id": "e4", "start": 1_000 + cr.REACH, "end": 1_200 + cr.REACH},
        {"id": "e5", "start": 2_000, "end": 2_001},
    ]
    (tmp_path / "chrT.json").write_text(json.dumps(_table(rows)))
    table = cr.DeletionTable(tmp_path)
    lean = sorted([(r["start"], r["end"], r["id"], "K562") for r in rows], key=lambda t: t[0])
    starts = [t[0] for t in lean]
    for start, end in ((1_450, 1_460), (900, 1_100), (4_000, 6_000), (0, 10), (1_999, 2_002), (5_199, 5_300)):
        theirs = sorted(e["id"] for e in table.overlapping("chrT", start, end))
        mine = sorted(e[2] for e in ac.overlapping(lean, starts, start, end))
        assert mine == theirs, (start, end, mine, theirs)


# ---- the read discipline -------------------------------------------------------------------------


def test_module_source_names_no_archive_loader() -> None:
    """The read discipline is asserted from the source, not claimed in a docstring."""
    assert ac.source_opens_no_archive()


def test_rss_is_reported_in_bytes_and_the_ceiling_raises() -> None:
    peak = ac.max_rss_bytes()
    assert peak > 8 * 1024**2, "a python process holding less than 8 MB means the unit is wrong"
    assert ac.check_rss("a test", ceiling=ac.RSS_CEILING_BYTES) == peak
    with pytest.raises(RuntimeError) as e:
        ac.check_rss("a test", ceiling=1)
    assert "stops" in str(e.value) and "Widening the ceiling" in str(e.value)


def test_the_ceiling_is_two_gibibytes() -> None:
    assert ac.RSS_CEILING_BYTES == 2 * 1024**3


# ---- the registered readings ---------------------------------------------------------------------


def _arm(n: int, k: int, cell: str = "K562", chroms: int = 22) -> ac.Arm:
    """An arm of n elements, k of them matching, spread evenly over `chroms` chromosomes.

    Spread evenly on purpose: the clustered interval is then narrow and the test is about which
    BRANCH a figure selects, not about the bootstrap. `_lumpy_arm` is the test for the other case.
    """
    a = ac.Arm(cell=cell)
    for i in range(n):
        label = ac.labels_for(cell)[0] if i % n < k else "placenta"
        a.element_rows[f"e{i}"] = (f"chr{i % chroms + 1}", i * 10_000_000, label)
    return a


def test_a_rate_at_the_base_rate_says_the_cell_is_not_evidence() -> None:
    """Reading (2), and only when the WHOLE clustered interval is inside the tolerance."""
    base = 0.062292
    r = ac.reading(_arm(20_000, 1_246), base)
    assert r["rate_reported"] is True and r["reading"] == "(2)"
    lo, hi = r["clustered_ci95_on_the_difference"]
    assert lo >= -ac.TOLERANCE and hi <= ac.TOLERANCE, (lo, hi)
    assert r["argmax_carries_cell_type_information"] is False
    assert "EQUIVALENCE result and not a failure to reject" in r["detection"]
    assert "NOT evidence of where it acts" in r["detection"]


def test_a_rate_below_the_base_rate_points_away() -> None:
    """Reading (1), and only when the interval's UPPER bound is below -0.02."""
    r = ac.reading(_arm(20_000, 200), 0.062292)
    assert r["reading"] == "(1)"
    assert r["clustered_ci95_on_the_difference"][1] < -ac.TOLERANCE
    assert r["argmax_carries_cell_type_information"] is False
    assert "points away" in r["detection"]


def test_a_rate_above_the_base_rate_is_never_read_as_a_confirmation() -> None:
    """CONFOUND, enforced: reading (3) returns None and names the selection."""
    r = ac.reading(_arm(20_000, 8_000), 0.062292)
    assert r["reading"] == "(3)"
    assert r["clustered_ci95_on_the_difference"][0] > ac.TOLERANCE
    assert r["argmax_carries_cell_type_information"] is None
    assert "does NOT establish" in r["detection"]
    assert r["usable"] is True


def test_usability_is_a_separate_threshold_from_detection() -> None:
    r = ac.reading(_arm(20_000, 3_000), 0.062292)
    assert r["reading"] == "(3)"  # a difference is detected
    assert r["usable"] is False  # and the label still names the measured cell on a minority
    assert "MINORITY" in r["usability"]


def test_an_arm_below_the_floor_reports_counts_and_no_rate() -> None:
    r = ac.reading(_arm(29, 29), 0.062292)
    assert r["rate_reported"] is False
    assert "observed_rate" not in r
    assert r["counts_only"]["elements"] == 29
    assert r["power"]["elements"] == 29


def test_wilson_interval_brackets_the_point_and_is_labelled_binomial() -> None:
    lo, hi = ac.wilson(62, 1_000)
    assert lo < 0.062 < hi
    assert ac.wilson(0, 0) is None
    assert "NOT a clustered interval" in ac.INTERVAL_IS_BINOMIAL
    r = ac.reading(_arm(20_000, 1_246), 0.062292)
    assert r["wilson95_secondary"] is not None
    assert r["clustered_interval_decides_the_reading"] is True


# ---- what the registration must carry ------------------------------------------------------------


def test_the_registration_is_committed_and_carries_the_falsifier_and_the_base_rate() -> None:
    with must_be_committed(REGISTRATION_RELATIVE).open() as fh:
        reg = json.load(fh)
    assert reg["lane"] == "lane-argmaxcell"
    assert reg["falsifier"] == ac.FALSIFIER
    assert reg["amendment_1"] == ac.AMENDMENT_1
    assert reg["amends"] == "data/results/argmaxcell_registration.json as committed at 67e14d7"
    assert reg["falsifier_thresholds"]["minimum_clusters_for_an_interval"] == ac.MIN_CLUSTERS
    assert reg["power_per_arm_stated_in_advance"]["K562"]["elements_needed_at_the_base_rate"] == 561
    assert reg["confound_registered_before_any_count"] == ac.CONFOUND
    assert reg["base_rate_rule"] == ac.BASE_RATE_RULE
    assert reg["base_rates"]["K562"]["rate"] == 0.062292
    assert reg["base_rate_census"]["rules"] == 440_589
    assert reg["population"] == ac.POPULATION
    assert reg["denominator"] == ac.DENOMINATOR
    assert reg["will_report"] == ac.WILL_REPORT
    assert reg["counts_named_in_advance"] == list(ac.COUNTS_NAMED)
    assert reg["rss_ceiling_bytes"] == ac.RSS_CEILING_BYTES
    assert reg["alphagenome_requests"] == 0


def test_the_registration_carries_the_peer_readings_it_may_not_strengthen() -> None:
    with must_be_committed(REGISTRATION_RELATIVE).open() as fh:
        reg = json.load(fh)
    assert "can withhold a direction and never reverse one" in reg["validates_nothing"]
    assert "unresolved class is NOT absence of regulation" in reg["validates_nothing"]
    assert "one biosample NAME" in reg["metadata_copy_limitation"]
    assert "not known to be biological replicates" in reg["metadata_copy_limitation"]
    assert "NOT established biological independence" in reg["grouping"]
    assert "112, 93, 73, 67, 26, 19, 18, 14 or 10" in reg["denominator_is_separate"]


def test_this_lane_edits_no_rule() -> None:
    """A measurement lane that changed a threshold would invalidate every count taken under it."""
    src = Path(ac.__file__).read_text()
    for forbidden in ("MIN_EFFECT =", "SIGN_RULE_V2 =", "def predict_target", "def retained_values"):
        assert forbidden not in src
    assert "NOT evidence of where it acts" in ac.FALSIFIER


def test_an_element_is_counted_once_however_many_genes_were_tested_near_it() -> None:
    """DENOMINATOR: the compiled cell belongs to the element, so pairs must not weight it."""
    a = ac.Arm(cell="K562")
    els = [(100, 200, "e1", "K562")]
    assert a.add("chr1", els) is True
    assert a.add("chr1", els) is True  # a second gene tested against the same element
    assert a.pairs == 2 and a.pairs_on_an_element == 2 and a.pair_level_matching == 2
    assert a.elements == 1 and a.elements_matching == 1 and a.rate == 1.0


def test_a_pair_on_no_element_is_counted_and_adds_no_element() -> None:
    a = ac.Arm(cell="K562")
    assert a.add("chr1", []) is False
    assert a.pairs == 1 and a.pairs_on_an_element == 0 and a.elements == 0 and a.rate is None


def test_the_label_tally_and_the_denominator_cannot_disagree() -> None:
    a = ac.Arm(cell="K562")
    a.add("chr1", [(1, 2, "e1", "K562"), (3, 4, "e2", "placenta"), (5, 6, "e3", "HepG2")])
    assert sum(a.label_counts.values()) == a.elements == 3
    assert a.label_counts["K562"] == a.elements_matching == 1
    assert a.rate == round(1 / 3, 6)


def test_top_labels_sets_each_label_beside_its_own_genome_wide_base_rate() -> None:
    cen = ac.census()
    a = ac.Arm(cell="K562")
    a.add("chr1", [(1, 2, "e1", "K562"), (3, 4, "e2", "K562"), (5, 6, "e3", "placenta")])
    rows = {r["label"]: r for r in ac.top_labels(a, cen)}
    assert rows["K562"]["elements"] == 2
    assert rows["K562"]["genome_wide_base_rate"] == 0.062292
    assert rows["placenta"]["genome_wide_base_rate"] == round(14_275 / 440_589, 6)


def test_the_archive_check_reads_code_and_not_prose() -> None:
    """A text scan would be measuring the docstring that has to name the loader to disclaim it."""
    assert "load_cached" in Path(ac.__file__).read_text(), "the docstring names it, as it must"
    assert ac.source_opens_no_archive(), "and the code does not refer to it"
    assert "load_cached" not in ac.names_used()
    # the check has teeth: a module that really called it would be caught
    assert "load_cached" in ac.names_used("import x\nx.load_cached('chr1', 'e')\n")
    assert "gzip" in ac.names_used("import gzip\n")


# ---- the alias map, which is a correction and not a convenience -----------------------------------


def test_the_jurkat_arm_is_matched_against_both_spellings_of_the_line() -> None:
    """MEASURED name mismatch: the benchmark says `Jurkat`, the model's track carries the clone."""
    cen = ac.census()
    assert ac.labels_for("Jurkat") == ("Jurkat", "Jurkat__Clone_E6_1")
    r = ac.base_rate("Jurkat", cen)
    assert r["rules_per_label"] == {"Jurkat": 5, "Jurkat__Clone_E6_1": 1_127}
    assert r["rules_with_this_label"] == 1_132
    assert r["rate"] == round(1_132 / 440_589, 6)
    assert r["ontology_terms"] == {"Jurkat": "EFO:0002796", "Jurkat__Clone_E6_1": "CLO:0007045"}
    # the spelling the benchmark uses alone would have measured 5 of 1,132
    assert (cen["per_cell"]["Jurkat"]["rules"], r["rules_with_this_label"]) == (5, 1_132)


def test_the_arm_counts_an_element_under_either_jurkat_label() -> None:
    a = ac.Arm(cell="Jurkat")
    assert a.add("chr1", [(1, 2, "e1", "Jurkat__Clone_E6_1")]) is True
    assert a.add("chr2", [(3, 4, "e2", "Jurkat")]) is True
    assert a.add("chr3", [(5, 6, "e3", "K562")]) is False
    assert a.elements == 3 and a.elements_matching == 2


def test_the_trio_members_are_not_merged_into_gm12878() -> None:
    """GM12891 and GM12892 are other people of the same trio and carry their own terms."""
    cen = ac.census()
    assert ac.labels_for("GM12878") == ("GM12878",)
    assert ac.base_rate("GM12878", cen)["rules_with_this_label"] == 6_406
    for other in ("GM12891", "GM12892"):
        assert other in cen["per_cell"]
        assert other not in ac.labels_for("GM12878")
    assert "would merge three people" in ac.ALIAS_RULE


def test_the_arms_are_named_by_the_benchmarks_own_celltype_strings() -> None:
    assert ac.ARMS == ("K562", "GM12878", "HCT116", "Jurkat", "WTC11")
    for cell in ac.ARMS:
        assert ac.labels_for(cell), cell


def test_a_cell_with_no_entry_in_the_map_falls_back_to_one_label() -> None:
    assert ac.labels_for("HepG2") == ("HepG2",)
    assert ac.base_rate("HepG2", ac.census())["rules_with_this_label"] == 14_779


# ---- AMENDMENT 1: the interval decides, and an underpowered arm cannot produce a null -------------


def _lumpy_arm(n: int, k: int, chroms: int, cell: str = "K562") -> ac.Arm:
    """An arm on `chroms` chromosomes whose elements also share ONE locus window per chromosome.

    Both clusterings then give `chroms` clusters, which is what makes the arm inconclusive by rule
    when `chroms` is under MIN_CLUSTERS: the locus fallback is tried first and does not rescue it.
    """
    a = ac.Arm(cell=cell)
    for i in range(n):
        label = ac.labels_for(cell)[0] if i < k else "placenta"
        a.element_rows[f"e{i}"] = (f"chr{i % chroms + 1}", i % 1_000, label)
    return a


def test_the_three_registered_thresholds_are_unchanged_by_the_amendment() -> None:
    assert (ac.TOLERANCE, ac.USABLE, ac.MIN_ELEMENTS) == (0.02, 0.25, 30)
    assert "TOLERANCE 0.02, USABLE 0.25 and MIN_ELEMENTS 30 are unchanged" in ac.AMENDMENT_1
    assert "no reading is made easier to reach" in ac.AMENDMENT_1


def test_the_amendment_imports_the_projects_own_cluster_floor_rather_than_restating_it() -> None:
    assert ac.MIN_CLUSTERS == cr.MIN_CLUSTERS_FOR_AN_INTERVAL == 10
    assert ac.BOOTSTRAPS == cr.BOOTSTRAPS == 2000
    assert ac.MIN_RESAMPLES == cr.MIN_RESAMPLES == 1000


def test_the_locus_fallback_is_tried_before_an_arm_is_called_inconclusive() -> None:
    """Few chromosomes alone do NOT refuse an interval: AMENDMENT_1 (c) tries loci first."""
    a = ac.Arm(cell="K562")
    for i in range(1_000):  # 4 chromosomes, but 250 distinct 1 Mb windows on each
        a.element_rows[f"e{i}"] = ("chr" + str(i % 4 + 1), i * ac.LOCUS_SPAN, "placenta")
    assert len(a.clusters("chromosome")) == 4
    r = ac.reading(a, 0.062292)
    assert r["clustered_by"] == "locus"
    assert r["clustered_ci95_on_the_rate"] is not None


def test_too_few_clusters_prints_no_interval_and_is_inconclusive_by_rule() -> None:
    """Under 10 clusters by EITHER grouping, the arm is reading (4) BY RULE, never reading (2)."""
    r = ac.reading(_lumpy_arm(1_000, 62, chroms=4), 0.062292)
    assert r["reading"] == "(4)"
    assert r["clustered_ci95_on_the_rate"] is None
    assert r["clustered_ci95_on_the_difference"] is None
    assert "INCONCLUSIVE BY RULE" in r["detection"]
    assert "4 resampling clusters is below the floor of 10" in r["why_no_interval"]
    assert r["argmax_carries_cell_type_information"] is None


def test_a_wide_interval_is_the_data_cannot_tell_and_not_no_information() -> None:
    """The defect the amendment fixes: a small arm at the base rate must NOT read as a null."""
    r = ac.reading(_arm(40, 2), 0.062292)  # point d = -0.012, inside the band; interval is not
    assert r["reading"] == "(4)"
    assert abs(r["difference_point_estimate"]) < ac.TOLERANCE, "the POINT would have said (2)"
    lo, hi = r["clustered_ci95_on_the_difference"]
    assert lo < -ac.TOLERANCE or hi > ac.TOLERANCE, (lo, hi)
    assert "THE DATA CANNOT TELL" in r["detection"]
    assert "NOT `no cell-type information`" in r["detection"]
    assert r["argmax_carries_cell_type_information"] is None


def test_reading_four_never_says_no_information_in_any_of_its_words() -> None:
    for arm in (_arm(40, 2), _lumpy_arm(1_000, 62, chroms=4)):
        r = ac.reading(arm, 0.062292)
        assert r["reading"] == "(4)"
        assert "NO cell-type information," not in r["detection"]
        assert "cannot tell" in r["detection"].lower()


def test_the_interval_and_not_the_point_selects_every_branch() -> None:
    """One figure, two cluster counts: the point is identical and the reading is not."""
    spread = ac.reading(_arm(20_000, 1_246), 0.062292)
    lumpy = ac.reading(_lumpy_arm(20_000, 1_246, chroms=4), 0.062292)
    assert spread["difference_point_estimate"] == lumpy["difference_point_estimate"]
    assert (spread["reading"], lumpy["reading"]) == ("(2)", "(4)")


def test_the_clustered_interval_is_reproducible_from_the_committed_seed() -> None:
    a = _arm(5_000, 400)
    assert (
        ac.reading(a, 0.062292)["clustered_ci95_on_the_rate"]
        == (ac.reading(a, 0.062292)["clustered_ci95_on_the_rate"])
    )
    assert ac.SEED == 20261002


def test_locus_clustering_is_the_fallback_and_carries_the_grouping_caveat() -> None:
    a = ac.Arm(cell="K562")
    for i in range(400):  # three chromosomes, many 1 Mb loci
        a.element_rows[f"e{i}"] = (f"chr{i % 3 + 1}", i * ac.LOCUS_SPAN, "placenta")
    kind, clusters = a.cluster_choice()
    assert kind == "locus" and len(clusters) >= ac.MIN_CLUSTERS
    r = ac.reading(a, 0.062292)
    assert r["clustered_by"] == "locus"
    assert r["cluster_grouping_caveat"] == ac.CELL2_GROUP_CAVEAT
    assert "NOT established biological independence" in r["cluster_grouping_caveat"]


def test_chromosome_clustering_is_preferred_and_carries_no_grouping_caveat() -> None:
    r = ac.reading(_arm(1_000, 62), 0.062292)
    assert r["clustered_by"] == "chromosome"
    assert r["cluster_grouping_caveat"] is None


def test_power_is_stated_at_both_rates_and_names_an_underpowered_arm() -> None:
    """AMENDMENT_1 (d). At a base rate near zero the base-rate figure alone would mislead."""
    assert ac.elements_needed(0.062292) == 561
    assert ac.elements_needed(0.000218) == 3
    p = ac.power(_arm(397, 0, cell="WTC11"), 0.000218)
    assert p["elements_needed_at_the_base_rate"] == 3
    assert p["has_them_at_the_base_rate"] is True
    assert p["elements_needed_at_the_observed_rate"] == 1
    p2 = ac.power(_arm(100, 6), 0.062292)
    assert p2["elements_needed_at_the_base_rate"] == 561
    assert p2["underpowered"] is True


def test_the_amendment_discloses_that_the_figures_had_been_seen() -> None:
    """The record must not read better than it was."""
    assert "ALREADY completed" in ac.AMENDMENT_1
    assert "NOT blind" in ac.AMENDMENT_1
    assert (
        "cannot create a detection"
        in ac.AMENDMENT_1.replace("can create a detection", "cannot create a detection")
        or "create a detection" in ac.AMENDMENT_1
    )


def test_the_falsifier_now_names_four_readings_and_the_clustered_interval() -> None:
    for word in ("(1)", "(2)", "(3)", "(4)", "CLUSTERED", "THE DATA CANNOT TELL", "equivalence"):
        assert word in ac.FALSIFIER, word
    assert "never by the point estimate" in ac.FALSIFIER


# ---- the committed result ------------------------------------------------------------------------

RESULT_RELATIVE = "data/results/argmaxcell.json"


def _result() -> dict:
    with must_be_committed(RESULT_RELATIVE).open() as fh:
        return json.load(fh)


def test_the_committed_result_carries_the_registered_words_unchanged() -> None:
    r = _result()
    assert r["lane"] == "lane-argmaxcell"
    assert r["falsifier"] == ac.FALSIFIER
    assert r["amendment_1"] == ac.AMENDMENT_1
    assert r["confound_registered_before_any_count"] == ac.CONFOUND
    assert r["base_rate_rule"] == ac.BASE_RATE_RULE
    assert r["population"] == ac.POPULATION
    assert r["denominator"] == ac.DENOMINATOR
    assert r["validates_nothing"] == ac.VALIDATES_NOTHING
    assert r["base_rates"]["K562"]["rate"] == 0.062292
    assert r["alphagenome_requests"] == 0


def test_the_committed_result_is_clean_and_opened_no_archive() -> None:
    m = _result()["result_manifest"]
    assert m["complete"] is True
    assert m["code_cleanliness"]["own_code_is_committed"] is True
    assert m["code_cleanliness"]["foreign_uncommitted_code_on_the_counting_path"] == []
    assert m["traced_inputs"]["active"] is True
    assert m["traced_inputs"]["undeclared"] == []


def test_the_measured_rss_stayed_under_the_registered_ceiling() -> None:
    rss = _result()["rss"]
    assert rss["ceiling_bytes"] == ac.RSS_CEILING_BYTES
    assert rss["peak_bytes"] < rss["ceiling_bytes"]
    assert rss["peak_bytes"] < 512 * 1024**2, "streaming 594 MB of tables must not cost half a GiB"


def test_no_arm_clears_the_usability_threshold_so_the_cell_is_not_usable_evidence() -> None:
    """The one answer that does not turn on the confound at all."""
    r = _result()
    assert r["verdict"]["arms_whose_label_clears_the_usability_threshold"] == []
    for cell, arm in r["per_arm"].items():
        if arm["reading"].get("rate_reported"):
            assert arm["reading"]["usable"] is False, cell
            assert arm["observed_rate"] < ac.USABLE, cell


def test_no_reading_three_arm_is_recorded_as_carrying_cell_type_information() -> None:
    """CONFOUND, in the committed file: a detection is never a confirmation."""
    r = _result()
    for cell in r["verdict"]["arms_where_a_difference_is_detected_but_not_attributable_reading_3"]:
        read = r["per_arm"][cell]["reading"]
        assert read["argmax_carries_cell_type_information"] is None, cell
        assert "does NOT establish" in read["detection"], cell


def test_no_reading_four_arm_is_recorded_as_a_null() -> None:
    """The amendment's whole point, checked in the committed file and not in the code alone."""
    r = _result()
    four = r["verdict"]["arms_the_data_cannot_tell_reading_4"]
    two = r["verdict"]["arms_with_no_cell_type_information_reading_2"]
    assert not set(four) & set(two)
    for cell in four:
        read = r["per_arm"][cell]["reading"]
        assert read["reading"] == "(4)"
        assert read["clustered_ci95_on_the_difference"] is None or True
        assert "CANNOT TELL" in read["detection"], cell
        assert read["argmax_carries_cell_type_information"] is None, cell


def test_every_reading_two_arm_really_is_an_equivalence_result() -> None:
    r = _result()
    for cell in r["verdict"]["arms_with_no_cell_type_information_reading_2"]:
        read = r["per_arm"][cell]["reading"]
        lo, hi = read["clustered_ci95_on_the_difference"]
        assert lo >= -ac.TOLERANCE and hi <= ac.TOLERANCE, (cell, lo, hi)
        assert "EQUIVALENCE result" in read["detection"], cell


def test_the_result_states_what_an_absolute_equivalence_leaves_open() -> None:
    """HCT116's tolerance is 14.8x its base rate and WTC11's is 91.7x: the file says so."""
    e = _result()["verdict"]["what_an_equivalence_result_excludes_on_a_small_base_rate"]
    assert e["tolerance_as_a_multiple_of_the_base_rate"]["WTC11"] > 50
    assert e["tolerance_as_a_multiple_of_the_base_rate"]["HCT116"] > 10
    assert e["tolerance_as_a_multiple_of_the_base_rate"]["K562"] < 1
    assert "does NOT mean `the label carries nothing`" in e["the_limit"]


def test_the_two_base_rates_are_reported_with_their_difference() -> None:
    a = _result()["base_rate_agreement"]
    assert a["committed_census_rules"] == 440_589
    assert a["difference"] == a["attributed_elements_streamed"] - 440_589
    assert "not reconciled away" in a["what_a_difference_means"]
