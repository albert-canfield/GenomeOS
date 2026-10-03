# SPDX-License-Identifier: AGPL-3.0-or-later
"""The rebuilt CRISPRi benchmark: its transcription of the committed result, and the rules it applies.

Two kinds of check live here. The first reads `data/results/crispri_benchmark.json` as committed and
holds every figure `scripts/crispri_benchmark_v2.py` quotes from it against it, so a mistyped number
fails here rather than being published as history. The second checks the rules on synthetic strata,
where the number of chromosomes and the cell types are set by the test: an interval below
`crispri.MIN_CLUSTERS_FOR_AN_INTERVAL` is withdrawn and not replaced, and a stratum carrying no
deletion value reports no number at all.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
import tracked_paths as tp

from genomeos import manifest as mf
from genomeos.attribution import crispri
from scripts import crispri_benchmark_v2 as v2

OLD = Path(v2.OLD)


def committed() -> dict[str, Any]:
    # The `pragma: no cover` here said it outright -- "the file is committed" -- and then skipped on
    # its absence anyway. `data/results/crispri_benchmark.json` is TRACKED, so a checkout without it
    # is not a stripped machine but a broken tree, and that is a failure and not a skip.
    tp.must_be_present(OLD, was=f"{OLD} is not in this checkout")
    return json.loads(OLD.read_text())


# --- the transcription ----------------------------------------------------------------------------


def test_the_committed_file_is_the_one_this_script_quotes() -> None:
    """Named so a reader knows which file the quotations below are checked against."""
    old = committed()
    assert old["result"] == "crispri_benchmark"
    assert old["date"] == v2.OLD_DATE
    assert old["result"] != v2.NAME, "the rebuild must not be written over the committed result"


@pytest.mark.parametrize("spec", v2.BEFORE, ids=[s["key"] for s in v2.BEFORE])
def test_every_quoted_interval_matches_the_committed_file(spec: dict[str, Any]) -> None:
    old = committed()
    found = v2.at(old, spec["path"])
    assert found is not v2.MISSING, f"{spec['path']} is not in {OLD}"
    for field, quoted in spec["before"].items():
        assert found[field] == quoted, f"{spec['path']}/{field}: quoted {quoted}, file has {found[field]}"
    assert v2.at(old, spec["pairs_path"]) == spec["pairs"]
    assert v2.at(old, spec["positives_path"]) == spec["positives"]


def test_the_quoted_population_counts_appear_in_the_population_sentence() -> None:
    """A figure without its population is unreadable, so the counts are in the prose as well as the field."""
    for spec in v2.BEFORE:
        assert f"{spec['pairs']:,}" in spec["population"] or str(spec["pairs"]) in spec["population"]
        assert str(spec["positives"]) in spec["population"]


def test_the_defect_being_rebuilt_is_the_one_the_committed_file_actually_has() -> None:
    """Five intervals, all at 200 draws, none recording a cluster count."""
    old = committed()
    found: list[tuple[str, dict[str, Any]]] = []

    def walk(node: Any, where: str = "") -> None:
        if isinstance(node, dict):
            for k, v in node.items():
                if k == "deletion_gain":
                    found.append((f"{where}/{k}", v))
                else:
                    walk(v, f"{where}/{k}")

    walk(old)
    assert len(found) == len(v2.BEFORE) == 5
    assert {p.lstrip("/") for p, _ in found} == {s["path"] for s in v2.BEFORE}
    for path, gain in found:
        assert gain["resamples"] == v2.OLD_DRAWS, path
        assert "clusters" not in gain, path
        assert "enough_clusters" not in gain, path


@pytest.mark.parametrize("path", sorted(v2.UNCHANGED_BEFORE))
def test_every_field_expected_unchanged_matches_the_committed_file(path: str) -> None:
    old = committed()
    found = v2.at(old, path)
    assert found is not v2.MISSING, f"{path} is not in {OLD}"
    assert found == v2.UNCHANGED_BEFORE[path]


def test_the_unchanged_set_is_as_wide_as_the_standard_asks() -> None:
    assert len(v2.UNCHANGED_BEFORE) >= 46


# --- the rules the rebuild applies ---------------------------------------------------------------


def _pairs(chroms: int, cell: str = "K562", per_chrom: int = 6) -> list[crispri.Pair]:
    """One synthetic stratum: `chroms` chromosomes, each with a regulated pair and some that are not."""
    out = []
    for c in range(chroms):
        for i in range(per_chrom):
            out.append(
                crispri.Pair(
                    chrom=f"chr{c + 1}",
                    start=1_000 * i,
                    end=1_000 * i + 500,
                    gene=f"g{c}_{i}",
                    cell=cell,
                    dataset="synthetic",
                    distance=10_000.0,
                    dhs=1.0,
                    h3k27ac=1.0,
                    regulated=i == 0,
                )
            )
    return out


def test_an_interval_on_too_few_clusters_is_withdrawn_and_not_replaced() -> None:
    rows = _pairs(crispri.MIN_CLUSTERS_FOR_AN_INTERVAL - 3)
    a = [1.0 if p.regulated else 0.0 for p in rows]
    b = [0.5 for _ in rows]
    g = crispri.gain_interval(a, b, rows, n=200)
    assert g["clusters"] == crispri.MIN_CLUSTERS_FOR_AN_INTERVAL - 3
    assert g["enough_clusters"] is False
    assert g["ci95"] is None, "no interval may be reported below the cluster minimum"
    assert "interval_unreliable" in g
    assert isinstance(g["gain"], float), "the point estimate stands; only the interval is withheld"


def test_an_interval_on_enough_clusters_is_reported_with_its_provenance() -> None:
    rows = _pairs(crispri.MIN_CLUSTERS_FOR_AN_INTERVAL + 2)
    a = [1.0 if p.regulated else 0.0 for p in rows]
    b = [0.5 for _ in rows]
    g = crispri.gain_interval(a, b, rows, n=crispri.BOOTSTRAPS)
    assert g["clusters"] == crispri.MIN_CLUSTERS_FOR_AN_INTERVAL + 2
    assert g["enough_clusters"] is True
    assert g["draws_requested"] == crispri.BOOTSTRAPS == 2000
    assert g["draws_dropped"] == g["draws_requested"] - g["resamples"]
    assert isinstance(g["ci95"], list) and len(g["ci95"]) == 2
    assert "interval_unreliable" not in g


def test_a_stratum_without_a_deletion_value_reports_no_number() -> None:
    rows = _pairs(12, cell="WTC11")
    assert crispri.deletion_available(rows) is False
    g = crispri.stratum_gain(rows, lambda: pytest.fail("the gain must not be computed here"))
    assert g["gain"] is None and g["ci95"] is None
    assert g["unavailable"] == crispri.UNAVAILABLE_GAIN


def test_an_empty_stratum_is_refused_for_its_own_reason_not_the_other_one() -> None:
    g = crispri.stratum_gain([], lambda: pytest.fail("the gain must not be computed here"))
    assert g["gain"] is None
    assert g["unavailable"] == crispri.NO_PAIRS_GAIN != crispri.UNAVAILABLE_GAIN


def test_a_stratum_with_a_deletion_value_is_measured() -> None:
    rows = _pairs(12)
    assert crispri.deletion_available(rows) is True
    assert crispri.stratum_gain(rows, lambda: {"gain": 0.1407})["gain"] == 0.1407


def test_a_stratum_mixing_a_scored_cell_with_an_unscored_one_is_refused() -> None:
    """A gain is one number over the whole stratum, so part of it without the feature refuses the whole."""
    assert crispri.deletion_available(_pairs(6) + _pairs(6, cell="Jurkat")) is False


# --- the two comparison blocks -------------------------------------------------------------------


def _fake_run(**over: Any) -> dict[str, Any]:
    """A result shaped like `crispri.score`'s, carrying only what the comparison blocks read."""
    out: dict[str, Any] = {}

    def put(path: str, value: Any, sep: str = "/") -> None:
        node = out
        keys = path.split(sep)
        for k in keys[:-1]:
            node = node.setdefault(k, {})
        node[keys[-1]] = value

    for spec in v2.BEFORE:
        put(spec["path"], {"gain": spec["before"]["gain"], "ci95": [0.0, 0.3], "clusters": 22})
        put(spec["pairs_path"], spec["pairs"])
        put(spec["positives_path"], spec["positives"])
    for path, value in over.items():
        put(path, value, sep="|")
    return out


def test_the_comparison_names_one_population_per_row_and_says_when_it_changed() -> None:
    block = v2.beside_the_committed_result(_fake_run())
    assert set(block["rows"]) == {s["key"] for s in v2.BEFORE}
    assert block["intervals_in_the_committed_result"] == 5
    for key, row in block["rows"].items():
        assert row["population"], key
        assert row["population_unchanged"] is True
        assert row["gain_moved"] == 0.0
        assert row[f"result_of_{v2.OLD_DATE}"]["clusters_recorded"] is False


def test_no_difference_is_taken_when_the_two_rows_are_about_different_pairs() -> None:
    run = _fake_run(**{"heldout|K562|models|distance|pairs": 1700})
    row = v2.beside_the_committed_result(run)["rows"]["heldout_K562"]
    assert row["population_unchanged"] is False
    assert "gain_moved" not in row, "a gain on other pairs is not a difference in this gain"
    assert "not about the same pairs" in row["gain_not_compared"]


def test_a_withdrawn_interval_and_a_refused_gain_are_reported_apart() -> None:
    run = _fake_run(
        **{
            "heldout|GM12878|deletion_gain": {
                "gain": 0.0352,
                "ci95": None,
                "clusters": 7,
                "enough_clusters": False,
                "interval_unreliable": "too few clusters",
            },
            "heldout|K562|deletion_gain": {"gain": None, "ci95": None, "unavailable": "no value"},
        }
    )
    block = v2.beside_the_committed_result(run)
    assert block["intervals_now_withdrawn"] == ["heldout_GM12878"]
    assert block["gains_now_refused"] == ["heldout_K562"]
    assert "heldout_GM12878" not in block["intervals_that_survive"]


def test_the_unchanged_block_names_every_field_that_moved() -> None:
    block = v2.what_did_not_change({"verdict": "failed"})
    assert block["fields_checked"] == len(v2.UNCHANGED_BEFORE)
    assert block["all_unchanged"] is False
    assert "verdict" in block["fields_that_moved"]
    assert block["fields"]["verdict"]["this_run"] == "failed"
    assert block["fields"]["coverage/training_pairs"]["present"] is False


# --- provenance ----------------------------------------------------------------------------------


def test_the_cleanliness_block_sets_the_uncommitted_files_against_the_counting_path() -> None:
    block = mf.code_cleanliness(v2.ENTRY, v2.OWN_CODE, v2.ROOT)
    assert block["counting_path_count"] == len(block["counting_path"]) > 1
    assert v2.ENTRY in block["counting_path"]
    assert "genomeos/attribution/crispri.py" in block["counting_path"]
    # Until 2026-10-02 this test also asserted `own_code_is_committed is True` and
    # `foreign_uncommitted_code_on_the_counting_path == []`. Those assert the state of this shared
    # working tree, not the behaviour of any code: four lanes commit here, so the test went red for a
    # peer's uncommitted file and green in a quiet moment, while the mechanism behind the fields was
    # exercised neither way. Measured: in two consecutive full runs that assertion was the only red in
    # about 4,160 tests, naming a different lane's file each time. What they stand on is now asserted
    # as invariants that hold in any tree state, in the style of tests/test_context_evidence.py, and
    # the three states they claimed to check are made deterministically in a planted repository in
    # tests/test_code_cleanliness_hermetic.py, with a counterfactual that strips the foreign-on-path
    # computation and must then fail.
    #
    # NO ENFORCEMENT MOVED, and the refusal that matters never lived here. scripts/check_staged.py
    # refuses the commit of any staged data/results/*.json whose block has `own_code_is_committed`
    # false or `foreign_uncommitted_code_on_the_counting_path` non-empty; save_result quarantines a
    # block that did not come from the shared function at all, which it checks by key set and not by
    # these two values. This test was a readiness check for writing a result, filed as a unit test.
    for p in block["foreign_uncommitted_code_on_the_counting_path"]:
        assert p in block["counting_path"] and p in block["foreign_uncommitted_code"]
    # 2026-10-02: the ENTRY SCRIPT is counted in this block now (genomeos.manifest.entry_script), and
    # when the commit does not hold it the script is NAMED here -- the lane ran it, so it is the lane's
    # own uncommitted code whatever OWN_CODE declared, which is precisely how 24 results asserted
    # `own_code_is_committed: True` while the bytes came from /private/tmp
    # (working tree at d4d0449, 19:41; superseded by e1def2e, 20:41).
    # Under pytest `sys.argv[0]`
    # is a test runner, so the one member beyond OWN_CODE is that and nothing else; the claim is
    # otherwise unchanged and `own_code_is_committed == (own_uncommitted_code == [])` still holds below.
    # The widened set is keyed to a value the SUBJECT supplies, so the containment alone would admit
    # a block naming ANY path at all -- a guard keyed to what the caller hands it is not a guard. The
    # form is asserted beside it, so the one admitted member must really be a test runner.
    assert block["entry_script"]["form"] == "test_runner"
    assert set(block["own_uncommitted_code"]) <= v2.OWN_CODE | {block["entry_script"]["argv0"]}
    assert block["own_code_is_committed"] == (block["own_uncommitted_code"] == [])
    assert not set(block["foreign_uncommitted_code"]) & v2.OWN_CODE


def test_the_two_fields_the_rebuild_requires_are_present_under_their_exact_names() -> None:
    block = mf.code_cleanliness(v2.ENTRY, v2.OWN_CODE, v2.ROOT)
    for field in ("own_code_is_committed", "foreign_uncommitted_code_on_the_counting_path"):
        assert field in block


# `v2.manifest([], [])` sha256-pins the EPCrispr benchmark table, so it cannot be composed without it,
# and it failed in a worktree of the committed tree because the table is git-ignored. The same store is
# already guarded elsewhere in the suite (tests/test_increase_links.py) and at line 32 of this file on a
# different input; it was missing here. Nothing below is weakened.
#
# ONLY the training_K562 table is named, and the narrowing is measured rather than assumed: a first
# draft also named a training_HCT116 twin, which is NOT in this checkout, and the test PASSES here --
# so naming it would have turned a running test into a silent skip. A marker that over-names is a
# weakening by the back door, and the only way to tell is to run it where the store is.
@pytest.mark.needs_local_data(
    "data/knowledge/crispri/EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz",
    how="the EPCrispr benchmark table is distilled into data/knowledge/crispri",
)
def test_the_manifest_states_the_interval_rules_and_no_model_request() -> None:
    m = v2.manifest([], [])
    assert mf.validate({**m, "code": {"git_sha": "x", "dirty": False}}) == []
    p = m["parameters"]
    assert p["bootstraps"] == crispri.BOOTSTRAPS == 2000
    assert p["min_clusters_for_an_interval"] == crispri.MIN_CLUSTERS_FOR_AN_INTERVAL
    assert p["alphagenome_requests"] == 0
    assert p["gain_where_unavailable"] == crispri.UNAVAILABLE_GAIN
    assert m["supersedes"]["file"] == v2.OLD and m["supersedes"]["date"] == v2.OLD_DATE
    assert "kept" in m["supersedes"]
