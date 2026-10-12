# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-execscale: the calibrated magnitude as an opt-in runtime parameter, census first.

The census (`scripts/execscale_census.py`) is registered BEFORE the parameter exists, because the
whole safety property of the change -- that at the default of 1.0 nothing committed moves -- cannot be
checked after the fact by anyone who did not write the committed state down first. These tests hold
the census instrument itself: that its claim names byte digests, that its comparison notices a moved
byte and a moved count, and that the magnitude it guards against is the fixed constant of `ac5580b`
described as the global AVERAGE it is.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

from genomeos import manifest as mf
from tests.tracked_paths import must_be_present

ROOT = Path(__file__).resolve().parents[1]


def _load(name: str):
    path = ROOT / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"_execscale_{name}", path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ec = _load("execscale_census")

TINY = """module tiny.one
gene g { max: 10; basal: 0.0 }
protein S {}
rule S activates g { strength: 1.0; threshold: 1.0; hill: 2 }
# test: rules == 1
"""


def test_the_claim_is_about_bytes_and_counts_and_says_the_default_moves_nothing():
    assert "byte-identical" in ec.CLAIM
    assert "1.0" in ec.CLAIM and "opt-in" in ec.CLAIM
    # the claim names the count it is held to, so no later reader has to guess which counts
    assert "# test: rules == 5176" in ec.CLAIM
    assert "sha256" in ec.VERIFICATION
    # the method is a digest comparison, stated as such, and not a reading of a diff
    assert "never a reading of a diff" in ec.VERIFICATION


def test_the_magnitude_is_the_fixed_constant_of_ac5580b_and_is_called_a_global_average():
    assert ec.MAGNITUDE["value"] == 0.2146
    assert ec.MAGNITUDE["source_commit"] == "ac5580b"
    assert "+0.04339" in ec.MAGNITUDE["reading"]
    assert "AVERAGE" in ec.MAGNITUDE["it_is_a_global_average"]
    assert "d9e3ddd" in ec.MAGNITUDE["it_is_a_global_average"]
    assert "fits nothing" in ec.MAGNITUDE["nothing_is_fitted_here"]
    # and the census claims nothing about RNA
    assert "says nothing about RNA" in ec.DOES_NOT_CLAIM


def test_the_census_reads_a_programmes_own_test_lines_verbatim(tmp_path: Path):
    p = tmp_path / "tiny.bio"
    p.write_text(TINY)
    assert ec.declared_tests(p) == ["# test: rules == 1"]


def test_the_committed_chr21_programme_still_declares_the_rule_count_the_claim_names():
    chr21 = ROOT / "data" / "organisms" / "human" / "noncoding_chr21.bio"
    assert "# test: rules == 5176" in ec.declared_tests(chr21)


def test_PLANT_one_moved_byte_moves_the_group_digest(tmp_path: Path):
    a, b = tmp_path / "a.bio", tmp_path / "b.bio"
    a.write_text(TINY)
    b.write_text(TINY)
    before = ec.group_digest([a, b])
    assert ec.group_digest([a, b]) == before  # the digest is stable over the same bytes
    b.write_text(TINY.replace("strength: 1.0", "strength: 0.9"))
    assert ec.group_digest([a, b]) != before


def test_PLANT_the_group_digest_notices_a_reordering_and_a_renaming(tmp_path: Path):
    a, b = tmp_path / "a.bio", tmp_path / "b.bio"
    a.write_text(TINY)
    b.write_text(TINY.replace("tiny.one", "tiny.two"))
    assert ec.group_digest([a, b]) != ec.group_digest([b, a])


def test_PLANT_the_field_comparison_catches_a_moved_count_and_a_dropped_field():
    before = {"bio_test": {"files": 58, "checks_passed": 410}, "programs": [{"rules": 5176}]}
    after_moved = {"bio_test": {"files": 58, "checks_passed": 409}, "programs": [{"rules": 5176}]}
    after_dropped = {"bio_test": {"files": 58}, "programs": [{"rules": 5176}]}
    flat = ec._flat(before)
    assert flat["bio_test.checks_passed"] == 410
    assert flat["programs[0].rules"] == 5176
    assert ec._flat(after_moved) != flat
    assert ec._flat(after_dropped) != flat
    # a comparison that only looked at the keys both records share would call the drop identical
    assert set(ec._flat(after_dropped)) < set(flat)


def test_the_census_keys_cover_every_counted_field_and_no_prose_field():
    assert ec.CENSUS_KEYS == (
        "programs",
        "program_count",
        "group_sha256",
        "compiled_rule_totals",
        "bio_test",
        "chr21",
    )
    for prose in ("claim", "verification", "magnitude", "does_not_claim", "question"):
        assert prose not in ec.CENSUS_KEYS


@pytest.mark.parametrize("tree", ec.PROGRAM_TREES)
def test_every_tree_the_claim_covers_exists_and_holds_programmes(tree: str):
    assert sorted((ROOT / tree).rglob("*.bio"))


# ---- the committed record ------------------------------------------------------------------
# These read data/results/execscale_census.json, which is committed in the same commit as they are.

RECORD = "data/results/execscale_census.json"


def _record() -> dict:
    must_be_present(RECORD)
    return json.loads((ROOT / RECORD).read_text())


def test_the_record_registers_the_claim_and_the_counts_it_is_held_to():
    d = _record()
    assert d["claim"] == ec.CLAIM
    assert d["verification"] == ec.VERIFICATION
    assert d["magnitude"]["value"] == 0.2146
    assert d["result_manifest"]["complete"] is True
    assert d["bio_test"]["totals"] == {
        "exit_code": 0,
        "files": 58,
        "checks_passed": 410,
        "checks_total": 410,
        "paths": list(ec.TEST_PATHS),
    }
    assert d["program_count"] == 58
    assert len(d["programs"]) == 58


def test_the_recorded_chr21_rule_total_is_the_one_the_programme_declares():
    d = _record()
    for path, lines in d["chr21"].items():
        assert "# test: rules == 5176" in lines
        assert d["compiled_rule_totals"][path] == 5176


def test_the_committed_programmes_are_byte_identical_to_the_record():
    """The claim's own check: every digest in the record against the bytes on disk now."""
    d = _record()
    for row in d["programs"]:
        must_be_present(row["path"])
        digest, size, _ = mf.sha256_of(ROOT / row["path"])
        assert (digest, size) == (row["sha256"], row["bytes"]), row["path"]
    assert ec.group_digest([ROOT / r["path"] for r in d["programs"]]) == d["group_sha256"]
