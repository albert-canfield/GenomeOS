# SPDX-License-Identifier: AGPL-3.0-or-later
"""Review item R9, third pass: the two contract gaps the second pass named, and the CRISPRi scorers'
own record carrying the benchmark's power columns.

1. `code.dirty` said true whenever any tracked file differed from the commit, so a chain of writers
   in one checkout marked every result after the first dirty with its predecessor's output. A
   modified file under data/results/ is output, not code: it is recorded apart and does not make
   the code dirty.
2. `manifest_rebuild` ignored only the top-level date, so a rebuild of a result with wall-clock
   "seconds" keys reported timings as differences. Timing keys are ignored at any depth, and the
   paths ignored are listed in the report, so nothing is hidden.
3. `crispri.Pair` dropped the PowerAtEffectSize columns that `measured.CrispriPair` carries.
"""

import io
import subprocess
import sys
from pathlib import Path

import pytest

from genomeos import manifest as mf
from genomeos.attribution import crispri

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import manifest_rebuild as mr  # noqa: E402


def _repo(tmp_path: Path) -> Path:
    def git(*a):
        subprocess.run(["git", "-C", str(tmp_path), *a], check=True, capture_output=True)

    git("init", "-q")
    git("config", "user.email", "t@example.org")
    git("config", "user.name", "t")
    (tmp_path / "data" / "results").mkdir(parents=True)
    (tmp_path / "genomeos").mkdir()
    (tmp_path / "data" / "results" / "a.json").write_text("{}")
    (tmp_path / "genomeos" / "x.py").write_text("X = 1\n")
    git("add", ".")
    git("commit", "-q", "-m", "base")
    return tmp_path


# 1. dirty ------------------------------------------------------------------------------------------


def test_a_modified_earlier_result_does_not_make_the_code_dirty(tmp_path):
    root = _repo(tmp_path)
    (root / "data" / "results" / "a.json").write_text('{"an": "earlier writer"}')
    rev = mf.code_revision(root)
    assert rev["dirty"] is False
    assert rev["dirty_result_paths"] == ["data/results/a.json"]
    assert rev["dirty_code_paths"] == []


def test_a_modified_code_file_still_makes_it_dirty(tmp_path):
    root = _repo(tmp_path)
    (root / "genomeos" / "x.py").write_text("X = 2\n")
    (root / "data" / "results" / "a.json").write_text('{"an": "earlier writer"}')
    rev = mf.code_revision(root)
    assert rev["dirty"] is True
    assert rev["dirty_code_paths"] == ["genomeos/x.py"]
    assert rev["dirty_result_paths"] == ["data/results/a.json"]


def test_a_clean_tree_is_clean(tmp_path):
    rev = mf.code_revision(_repo(tmp_path))
    assert rev["dirty"] is False and rev["dirty_result_paths"] == [] and rev["git_sha"]


def test_outside_git_dirty_stays_unknown(tmp_path):
    assert mf.code_revision(tmp_path)["dirty"] is None


# 2. timing fields ------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "key", ["seconds", "maf_seconds", "seconds_this_run", "seconds_per_megabase", "chains_per_second"]
)
def test_timing_keys_are_ignored_at_any_depth(key):
    a = {"date": "2026-09-27", "x": 1, "deep": [{"y": 2, key: 10.5}]}
    b = {"date": "2026-09-28", "x": 1, "deep": [{"y": 2, key: 99.1}]}
    assert mr.diff(mr.comparable(a), mr.comparable(b)) == []
    assert mr.timing_paths(a) == [f"/deep[0]/{key}"]


@pytest.mark.parametrize(
    "key", ["second_endpoint_meets_success", "second", "duration", "case_for_a_second_rate"]
)
def test_a_key_that_only_looks_like_timing_is_still_compared(key):
    a, b = {key: 1}, {key: 2}
    assert mr.diff(mr.comparable(a), mr.comparable(b)) == [f"/{key}: 1 vs 2"]


def test_a_changed_figure_beside_a_timing_key_is_still_reported():
    a = {"auprc": 0.724, "seconds": 3.0}
    b = {"auprc": 0.725, "seconds": 4.0}
    assert mr.diff(mr.comparable(a), mr.comparable(b)) == ["/auprc: 0.724 vs 0.725"]


def test_the_run_record_is_ignored_but_the_rest_of_the_manifest_is_compared():
    a = {mf.KEY: {"code": {"git_sha": "a"}, "inputs": [{"sha256": "1"}]}}
    b = {mf.KEY: {"code": {"git_sha": "b"}, "inputs": [{"sha256": "2"}]}}
    assert mr.diff(mr.comparable(a), mr.comparable(b)) == [f"/{mf.KEY}/inputs[0]/sha256: '1' vs '2'"]


# 3. crispri.Pair power columns ------------------------------------------------------------------------

POWER = [f"PowerAtEffectSize{n}" for n in (10, 15, 20, 25, 50)]
HEADER = (
    "chrom\tchromStart\tchromEnd\tmeasuredGeneSymbol\tValidConnection\tCellType\tRegulated\tDataset\t"
    "distanceToTSS\tDHS.RPM\tH3K27ac.RPM\t" + "\t".join(POWER) + "\n"
)


def test_pair_carries_the_power_columns_and_none_for_an_empty_cell():
    text = HEADER + "chr1\t100\t600\tA\tTRUE\tK562\tFALSE\tD\t5000\t1\t1\t0.1\t0.3\t0.85\tNA\t\n"
    (p,) = crispri.parse(io.StringIO(text))
    assert (
        p.power_at_effect_size_10,
        p.power_at_effect_size_15,
        p.power_at_effect_size_20,
        p.power_at_effect_size_25,
        p.power_at_effect_size_50,
    ) == (0.1, 0.3, 0.85, None, None)


def test_a_table_without_power_columns_parses_as_before():
    text = (
        "chrom\tchromStart\tchromEnd\tmeasuredGeneSymbol\tValidConnection\tCellType\tRegulated\tDataset\t"
        "distanceToTSS\tDHS.RPM\tH3K27ac.RPM\nchr1\t100\t600\tA\tTRUE\tK562\tTRUE\tD\t5000\t1\t1\n"
    )
    (p,) = crispri.parse(io.StringIO(text))
    assert p.power_at_effect_size_20 is None and p.regulated and p.features == {}


def test_no_scorer_feature_reads_a_power_column():
    names = {n for cols in crispri.FEATURES.values() for n in cols}
    names |= {n for cols in crispri.CONTACT_FEATURES.values() for n in cols}
    names |= {n for cols in crispri.DNASE_FEATURES.values() for n in cols}
    assert not {n for n in names if "power" in n.lower()}


def test_the_two_records_name_the_same_power_fields():
    from dataclasses import fields

    from genomeos.attribution.measured import CrispriPair

    ours = {f.name for f in fields(crispri.Pair) if f.name.startswith("power_")}
    theirs = {f.name for f in fields(CrispriPair) if f.name.startswith("power_")}
    assert ours == theirs and len(ours) == 5
