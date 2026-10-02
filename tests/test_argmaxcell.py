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
    assert r["label"] == "K562"
    assert r["rules_total"] == cen["rules"] == 440_589
    assert r["rules_with_this_label"] == 27_445
    assert r["rate"] == round(27_445 / 440_589, 6) == 0.062292
    assert r["source"] == ac.CENSUS_RELATIVE


def test_base_rate_label_is_the_compilers_own_function() -> None:
    """A cell name is looked up under the label the compiler would have written, not under itself."""
    cen = ac.census()
    assert cp.context("IMR-90") == "IMR_90"
    assert ac.base_rate("IMR-90", cen)["label"] == "IMR_90"
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


def _arm(n: int, k: int, cell: str = "K562") -> ac.Arm:
    a = ac.Arm(cell=cell, label=cp.context(cell))
    for i in range(n):
        a.element_label[f"e{i}"] = cp.context(cell) if i < k else "placenta"
    return a


def test_a_rate_at_the_base_rate_says_the_cell_is_not_evidence() -> None:
    """The branch the question turns on: at the base rate, the cell carries no information."""
    base = 0.062292
    r = ac.reading(_arm(1_000, 62), base)
    assert r["rate_reported"] is True
    assert abs(r["difference"]) < ac.TOLERANCE
    assert r["argmax_carries_cell_type_information"] is False
    assert "NO cell-type information is detected" in r["detection"]
    assert "NOT evidence of where it acts" in r["detection"]


def test_a_rate_below_the_base_rate_points_away() -> None:
    r = ac.reading(_arm(1_000, 10), 0.062292)
    assert r["difference"] <= -ac.TOLERANCE
    assert r["argmax_carries_cell_type_information"] is False
    assert "points away" in r["detection"]


def test_a_rate_above_the_base_rate_is_never_read_as_a_confirmation() -> None:
    """CONFOUND, enforced: the positive branch returns None and names the selection."""
    r = ac.reading(_arm(1_000, 400), 0.062292)
    assert r["difference"] >= ac.TOLERANCE
    assert r["argmax_carries_cell_type_information"] is None
    assert "does NOT establish" in r["detection"]
    assert r["usable"] is True


def test_usability_is_a_separate_threshold_from_detection() -> None:
    r = ac.reading(_arm(1_000, 150), 0.062292)
    assert r["difference"] >= ac.TOLERANCE  # a difference is detected
    assert r["usable"] is False  # and the label still names the measured cell on a minority
    assert "MINORITY" in r["usability"]


def test_an_arm_below_the_floor_reports_counts_and_no_rate() -> None:
    r = ac.reading(_arm(29, 29), 0.062292)
    assert r["rate_reported"] is False
    assert "observed_rate" not in r
    assert r["counts_only"]["elements"] == 29


def test_wilson_interval_brackets_the_point_and_is_labelled_binomial() -> None:
    lo, hi = ac.wilson(62, 1_000)
    assert lo < 0.062 < hi
    assert ac.wilson(0, 0) is None
    assert "NOT a clustered interval" in ac.INTERVAL_IS_BINOMIAL


# ---- what the registration must carry ------------------------------------------------------------


def test_the_registration_is_committed_and_carries_the_falsifier_and_the_base_rate() -> None:
    with must_be_committed(REGISTRATION_RELATIVE).open() as fh:
        reg = json.load(fh)
    assert reg["lane"] == "lane-argmaxcell"
    assert reg["falsifier"] == ac.FALSIFIER
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
    a = ac.Arm(cell="K562", label="K562")
    els = [(100, 200, "e1", "K562")]
    assert a.add(els) is True
    assert a.add(els) is True  # a second gene tested against the same element
    assert a.pairs == 2 and a.pairs_on_an_element == 2 and a.pair_level_matching == 2
    assert a.elements == 1 and a.elements_matching == 1 and a.rate == 1.0


def test_a_pair_on_no_element_is_counted_and_adds_no_element() -> None:
    a = ac.Arm(cell="K562", label="K562")
    assert a.add([]) is False
    assert a.pairs == 1 and a.pairs_on_an_element == 0 and a.elements == 0 and a.rate is None


def test_the_label_tally_and_the_denominator_cannot_disagree() -> None:
    a = ac.Arm(cell="K562", label="K562")
    a.add([(1, 2, "e1", "K562"), (3, 4, "e2", "placenta"), (5, 6, "e3", "HepG2")])
    assert sum(a.labels.values()) == a.elements == 3
    assert a.labels["K562"] == a.elements_matching == 1
    assert a.rate == round(1 / 3, 6)


def test_top_labels_sets_each_label_beside_its_own_genome_wide_base_rate() -> None:
    cen = ac.census()
    a = ac.Arm(cell="K562", label="K562")
    a.add([(1, 2, "e1", "K562"), (3, 4, "e2", "K562"), (5, 6, "e3", "placenta")])
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
