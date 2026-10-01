# SPDX-License-Identifier: AGPL-3.0-or-later
"""The rebuild of the CRISPRi published-benchmark result under a new name (scripts/crispri_published_v2).

Three things are checked without running the measurement: that the 2026-09-27 figures this script
quotes are the ones the committed result actually holds, that the side-by-side block names a
population for every row and refuses a number exactly where the feature is absent, and that the
import closure is found by following imports rather than by listing files.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from genomeos.attribution import crispri  # noqa: E402
from scripts import crispri_published_v2 as v2  # noqa: E402

OLD = ROOT / "data" / "results" / "crispri_published.json"


@pytest.fixture(scope="module")
def old() -> dict:
    if not OLD.is_file():
        pytest.skip(f"{OLD} is not in this checkout")
    return json.loads(OLD.read_text())


# --- the quoted figures are the committed ones ----------------------------------------------------


def test_every_quoted_interval_is_the_one_the_2026_09_27_file_holds(old: dict) -> None:
    for key, path, _population, before in v2.BEFORE:
        was = v2.at(old, path)
        assert isinstance(was, dict), f"{key}: {path} is not in {OLD}"
        assert was["gain"] == before["gain"], key
        assert was["ci95"] == before["ci95"], key
        assert was["resamples"] == before["resamples"] == v2.OLD_DRAWS, key


def test_every_unchanged_field_is_transcribed_from_the_same_file(old: dict) -> None:
    for path, before in v2.UNCHANGED_BEFORE.items():
        assert v2.at(old, path) == before, path


def test_the_old_file_is_the_one_with_two_hundred_draws(old: dict) -> None:
    assert old["date"] == v2.OLD_DATE
    assert old["result"] == "crispri_published"


def test_the_three_strata_without_a_value_are_the_ones_that_reported_a_number(old: dict) -> None:
    per_cell = old["heldout_published_pairs"]["per_cell_type_weighted"]
    reported = sorted(c for c, v in per_cell.items() if not v["deletion_available"])
    assert reported == sorted(v2.WITHOUT_A_DELETION_VALUE)
    for cell in reported:
        assert per_cell[cell]["deletion_gain"]["gain"] is not None, cell


def test_the_rebuild_writes_a_different_name_from_the_result_it_sits_beside() -> None:
    assert v2.NAME != "crispri_published"
    assert v2.OLD.endswith("crispri_published.json")


def test_the_draw_count_this_run_asks_for_is_two_thousand() -> None:
    assert crispri.BOOTSTRAPS == 2000
    assert crispri.BOOTSTRAPS > v2.OLD_DRAWS


# --- the side-by-side block -----------------------------------------------------------------------


def _interval(gain: float, clusters: int = 23, requested: int = 2000, kept: int = 1990) -> dict:
    return {
        "gain": gain,
        "ci95": [gain - 0.01, gain + 0.01],
        "resamples": kept,
        "clusters": clusters,
        "draws_requested": requested,
        "draws_dropped": requested - kept,
        "met_minimum": kept >= crispri.MIN_RESAMPLES,
    }


def _result() -> dict:
    """A result shaped like the real one, with one interval per BEFORE row."""
    out: dict = {}
    for i, (_key, path, _pop, before) in enumerate(v2.BEFORE):
        node = out
        keys = path.split("/")
        for k in keys[:-1]:
            node = node.setdefault(k, {})
        gain = round(before["gain"] + 0.001 * i, 4)
        node[keys[-1]] = (
            crispri.gain_unavailable()
            if any(c in path for c in v2.WITHOUT_A_DELETION_VALUE)
            else _interval(gain)
        )
    return out


def test_every_row_names_a_population_and_quotes_the_old_file() -> None:
    block = v2.beside_the_old_result(_result())
    assert block["rows_not_found_in_this_run"] == []
    assert len(block["rows"]) == len(v2.BEFORE)
    for key, row in block["rows"].items():
        assert row["population"].strip(), key
        old = row[f"result_of_{v2.OLD_DATE}"]
        assert old["draws"] == v2.OLD_DRAWS, key
        assert old["file"] == v2.OLD, key


def test_the_two_held_out_populations_are_named_apart() -> None:
    rows = v2.beside_the_old_result(_result())["rows"]
    assert "all 4,378 held-out pairs" in rows["heldout_pooled"]["population"]
    assert "the 1,918 K562 held-out pairs" in rows["heldout_K562"]["population"]
    assert rows["heldout_pooled"]["population"] != rows["heldout_K562"]["population"]


def test_a_refused_gain_is_reported_as_refused_and_nothing_else() -> None:
    block = v2.beside_the_old_result(_result())
    assert block["gains_now_refused"] == ["heldout_HCT116", "heldout_Jurkat", "heldout_WTC11"]
    for key in block["gains_now_refused"]:
        row = block["rows"][key]
        assert row["this_run"]["gain"] is None
        assert row["this_run"]["ci95"] is None
        assert row["this_run"]["unavailable"] == crispri.UNAVAILABLE_GAIN
        assert row["gain_moved"] is None
        assert row["ci95_width"]["this_run"] is None
        assert row[f"result_of_{v2.OLD_DATE}"]["gain"] is not None


def test_an_available_row_carries_the_provenance_of_its_interval() -> None:
    row = v2.beside_the_old_result(_result())["rows"]["heldout_K562"]
    assert row["this_run"]["draws_requested"] == 2000
    assert row["this_run"]["clusters"] == 23
    assert row["this_run"]["draws_dropped"] == 10
    assert row["gain_moved"] is not None


def test_a_row_whose_field_is_missing_is_named_rather_than_passed_over() -> None:
    result = _result()
    del result["heldout_published_pairs"]["per_cell_type_weighted"]["K562"]
    block = v2.beside_the_old_result(result)
    assert block["rows_not_found_in_this_run"] == ["heldout_K562"]
    assert block["rows"]["heldout_K562"]["this_run"] == {}


def test_the_reading_sentence_states_the_draws_of_both_runs() -> None:
    reading = v2.beside_the_old_result(_result())["reading"]
    assert str(v2.OLD_DRAWS) in reading
    assert str(crispri.BOOTSTRAPS) in reading


# --- what did not change --------------------------------------------------------------------------


def test_an_unchanged_field_that_moved_is_listed(old: dict) -> None:
    block = v2.what_did_not_change(old)
    assert block["fields_checked"] == len(v2.UNCHANGED_BEFORE)
    assert block["fields_that_moved"] == []
    assert block["all_unchanged"] is True

    moved = json.loads(json.dumps(old))
    moved["heldout_published_pairs"]["models"]["activity + distance + deletion"]["auprc"] = 0.9
    after = v2.what_did_not_change(moved)
    assert after["fields_that_moved"] == [
        "heldout_published_pairs/models/activity + distance + deletion/auprc"
    ]
    assert after["all_unchanged"] is False


def test_a_missing_field_counts_as_not_the_same(old: dict) -> None:
    gone = json.loads(json.dumps(old))
    del gone["alphagenome_requests"]
    block = v2.what_did_not_change(gone)
    assert block["fields"]["alphagenome_requests"]["present"] is False
    assert block["fields"]["alphagenome_requests"]["same"] is False


def test_the_matched_arm_draws_are_not_what_this_run_changed() -> None:
    assert crispri.MATCH_DRAWS == 1000
    assert v2.UNCHANGED_BEFORE["coverage_arms_k562_heldout/arm2_coverage_matched/draws"] == 1000


def test_the_request_count_stays_zero() -> None:
    assert v2.UNCHANGED_BEFORE["alphagenome_requests"] == 0


# --- the import closure ---------------------------------------------------------------------------


def test_the_closure_holds_the_entry_and_what_it_imports() -> None:
    files = v2.import_closure()
    assert v2.ENTRY in files
    for expected in (
        "scripts/crispri_published.py",
        "genomeos/attribution/crispri.py",
        "genomeos/manifest.py",
        "genomeos/results.py",
        "genomeos/__init__.py",
        "genomeos/attribution/__init__.py",
    ):
        assert expected in files, expected


def test_the_closure_reaches_an_import_two_steps_away() -> None:
    """crispri.py imports genomeos.genome.epigenome, which this script never names itself."""
    source = (ROOT / v2.ENTRY).read_text()
    assert "epigenome" not in source
    assert "genomeos/genome/epigenome.py" in v2.import_closure()


def test_the_closure_is_sorted_and_free_of_duplicates() -> None:
    files = v2.import_closure()
    assert files == sorted(files)
    assert len(files) == len(set(files))


def test_the_closure_names_nothing_outside_the_repository() -> None:
    for f in v2.import_closure():
        assert not f.startswith("/")
        assert (ROOT / f).is_file(), f
        assert f.endswith(".py")


def test_the_closure_is_computed_not_listed(tmp_path: Path) -> None:
    """A tree of three modules it has never seen: the closure must follow them, and stop at the leaf."""
    (tmp_path / "pkg").mkdir()
    (tmp_path / "pkg" / "__init__.py").write_text("")
    (tmp_path / "entry.py").write_text("import json\nfrom pkg import a\n")
    (tmp_path / "pkg" / "a.py").write_text("from .b import thing\n")
    (tmp_path / "pkg" / "b.py").write_text("import os\n\nthing = 1\n")
    (tmp_path / "pkg" / "unused.py").write_text("raise AssertionError\n")
    assert v2.import_closure("entry.py", tmp_path) == [
        "entry.py",
        "pkg/__init__.py",
        "pkg/a.py",
        "pkg/b.py",
    ]


def test_a_relative_import_above_the_package_is_resolved(tmp_path: Path) -> None:
    (tmp_path / "pkg" / "sub").mkdir(parents=True)
    for p in ("pkg/__init__.py", "pkg/sub/__init__.py"):
        (tmp_path / p).write_text("")
    (tmp_path / "entry.py").write_text("from pkg.sub import leaf\n")
    (tmp_path / "pkg" / "sub" / "leaf.py").write_text("from ..top import x\n")
    (tmp_path / "pkg" / "top.py").write_text("x = 1\n")
    assert v2.import_closure("entry.py", tmp_path) == [
        "entry.py",
        "pkg/__init__.py",
        "pkg/sub/__init__.py",
        "pkg/sub/leaf.py",
        "pkg/top.py",
    ]


def test_a_cycle_does_not_stop_the_closure(tmp_path: Path) -> None:
    (tmp_path / "entry.py").write_text("import one\n")
    (tmp_path / "one.py").write_text("import two\n")
    (tmp_path / "two.py").write_text("import one\n")
    assert v2.import_closure("entry.py", tmp_path) == ["entry.py", "one.py", "two.py"]


def test_the_block_sets_the_uncommitted_files_against_the_closure() -> None:
    block = v2.closure_block()
    assert block["entry"] == v2.ENTRY
    assert block["count"] == len(block["files"]) > 1
    assert block["uncommitted_code_on_the_closure"] == [], (
        "an uncommitted file in this shared checkout is on this script's import closure: the result "
        f"cannot claim a clean rebuild until it is committed: {block['uncommitted_code_on_the_closure']}"
    )


# --- the manifest ---------------------------------------------------------------------------------


def test_the_manifest_says_why_the_rebuild_exists_and_keeps_the_old_file() -> None:
    assert "2,000" in v2.WHY_V2
    assert "kept unchanged" in v2.WHY_V2
    assert v2.OLD in v2.WHY_V2
