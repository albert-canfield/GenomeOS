# SPDX-License-Identifier: AGPL-3.0-or-later
"""Item 12 S6 (second external review, 2026-09-28): the evidence contract enforced.

Two defects in the contract built for R9 (fe0880a), negatives first:

1. `save_result` wrote a new result before refusing it, and decided "historical" by whether the
   file existed, so a retry found its own file and only warned. Now a failed new result goes to a
   quarantine outside the registry, and "historical" is an explicit committed list.
2. The revision stamp ran `git status --untracked-files=no`, so an untracked script could write a
   result stamped clean. Now untracked code makes the stamp dirty and is named.
"""

import json
import subprocess
import sys
from pathlib import Path

import pytest

from genomeos import manifest as mf
from genomeos import results
from genomeos.results import ManifestWarning, save_result

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import manifest_legacy as ml  # noqa: E402


def full(tmp_path: Path) -> dict:
    src = tmp_path / "in.tsv"
    src.write_text("chrom\tstart\tend\nchr21\t0\t10\n")
    return {
        "sources": [{"accession": "ENCSR000XXX", "version": "v1"}],
        "inputs": [mf.input_entry(src, partition="heldout")],
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {"threshold": 0.5},
        "exclusions": [],
        "partitions": {"heldout": "the test"},
    }


@pytest.fixture
def registry(tmp_path, monkeypatch):
    """A registry at tmp/data/results with a two-name legacy allowlist beside it."""
    reg = tmp_path / "data" / "results"
    reg.mkdir(parents=True)
    allow = tmp_path / "data" / "results_legacy.txt"
    allow.write_text("# a header line\nold_result\ttracked\nold_reader_chr1\tignored-local\n")
    monkeypatch.setattr(results, "RESULTS_DIR", reg)
    monkeypatch.setattr(results, "LEGACY_ALLOWLIST", allow)
    return reg


# 1. the retry path -------------------------------------------------------------------------------------


def test_a_retried_failed_result_never_reaches_the_registry(registry):
    for attempt in (1, 2, 3):
        with pytest.raises(mf.ManifestError, match="quarantined at"):
            save_result("brand_new", {"value": attempt}, registry)
        assert not (registry / "brand_new.json").exists()
    q = results.quarantine_dir(registry) / "brand_new.json"
    kept = json.loads(q.read_text())
    assert kept["value"] == 3 and kept[mf.KEY]["complete"] is False
    assert kept["quarantine"]["meant_for"] == str(registry / "brand_new.json")
    assert "not on the legacy allowlist" in kept["quarantine"]["reason"]
    assert "missing sources" in kept["quarantine"]["problems"]
    assert q.parent == registry.parent / "quarantine" / "results"


def test_an_existing_file_does_not_make_a_name_historical(registry):
    """The retry path's cause: a file left in the registry by any means is not an allowlist entry."""
    (registry / "brand_new.json").write_text('{"left": "by an earlier failed run"}')
    with pytest.raises(mf.ManifestError):
        save_result("brand_new", {"value": 1}, registry)
    assert json.loads((registry / "brand_new.json").read_text()) == {"left": "by an earlier failed run"}


def test_strict_false_cannot_admit_a_name_off_the_allowlist(registry):
    with pytest.raises(mf.ManifestError):
        save_result("brand_new", {"value": 1}, registry, strict=False)
    assert not (registry / "brand_new.json").exists()


def test_an_unreadable_allowlist_fails_closed(registry, monkeypatch):
    monkeypatch.setattr(results, "LEGACY_ALLOWLIST", registry.parent / "absent.txt")
    with pytest.raises(mf.ManifestError):
        save_result("old_result", {"value": 1}, registry)
    assert not (registry / "old_result.json").exists()


def test_the_quarantine_is_git_ignored(tmp_path):
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    reg = tmp_path / "data" / "results"
    with pytest.raises(mf.ManifestError):
        save_result("failed", {"value": 1}, reg, strict=True)
    q = results.quarantine_dir(reg) / "failed.json"
    assert q.exists()
    shown = subprocess.run(
        ["git", "-C", str(tmp_path), "status", "--porcelain", "--untracked-files=all"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    assert "quarantine" not in shown
    ignored = subprocess.run(["git", "-C", str(tmp_path), "check-ignore", "-q", str(q)])
    assert ignored.returncode == 0


# 2. what still works ------------------------------------------------------------------------------------


def test_an_allowlisted_historical_result_still_regenerates(registry):
    (registry / "old_result.json").write_text('{"value": 1}')
    with pytest.warns(ManifestWarning, match="manifest incomplete"):
        p = save_result("old_result", {"value": 2}, registry)
    assert p == registry / "old_result.json" and json.loads(p.read_text())["value"] == 2
    # on a machine where the historical file was never written (a git-ignored one) it regenerates too
    with pytest.warns(ManifestWarning):
        save_result("old_reader_chr1", {"value": 3}, registry)
    assert (registry / "old_reader_chr1.json").exists()
    assert not results.quarantine_dir(registry).exists()


def test_a_complete_new_result_enters_and_supersedes_its_quarantined_attempt(registry, tmp_path):
    with pytest.raises(mf.ManifestError):
        save_result("brand_new", {"value": 1}, registry)
    q = results.quarantine_dir(registry) / "brand_new.json"
    assert q.exists()
    p = save_result("brand_new", {"value": 2}, registry, manifest=full(tmp_path))
    assert json.loads(p.read_text())[mf.KEY]["complete"] and not q.exists()


def test_the_committed_allowlist_is_the_registered_list():
    listed = ml.read(ROOT / results.LEGACY_ALLOWLIST)
    assert len(listed) == results.LEGACY_COUNT == 955
    assert set(listed.values()) <= set(ml.SOURCES)
    assert sum(v == "tracked" for v in listed.values()) == 656
    names = [
        line.split("\t")[0]
        for line in (ROOT / results.LEGACY_ALLOWLIST).read_text().splitlines()
        if not line.startswith("#")
    ]
    assert names == sorted(set(names))
    assert results.legacy_names(ROOT / results.LEGACY_ALLOWLIST) == frozenset(listed)


def test_the_tracked_part_of_the_allowlist_matches_git_history():
    have = subprocess.run(
        ["git", "-C", str(ROOT), "cat-file", "-e", f"{mf.CONTRACT_COMMIT}^"], capture_output=True
    )
    if have.returncode:
        pytest.skip("the contract commit is not in this clone's history")
    tree = subprocess.run(
        ["git", "-C", str(ROOT), "ls-tree", "--name-only", f"{mf.CONTRACT_COMMIT}^", "data/results/"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.split()
    tracked = {Path(x).stem for x in tree if x.endswith(".json")}
    listed = ml.read(ROOT / results.LEGACY_ALLOWLIST)
    assert {n for n, s in listed.items() if s == "tracked"} == tracked


def test_the_generator_never_rewrites_the_list(tmp_path):
    out = tmp_path / "legacy.txt"
    out.write_text("# an existing list\nx\ttracked\n")
    assert ml.main(["--out", str(out)]) == 2
    assert out.read_text() == "# an existing list\nx\ttracked\n"


def test_no_result_off_the_allowlist_in_the_registry_is_incomplete():
    """The invariant the enforcement keeps: every result in data/results is complete or allowlisted."""
    legacy = results.legacy_names(ROOT / results.LEGACY_ALLOWLIST)
    bad = []
    for p in sorted((ROOT / "data" / "results").glob("*.json")):
        if p.stem not in legacy and not mf.read(json.loads(p.read_text()))["complete"]:
            bad.append(p.stem)
    assert bad == []


# 3. the revision stamp -----------------------------------------------------------------------------------


def _repo(tmp_path: Path) -> Path:
    def git(*a):
        subprocess.run(["git", "-C", str(tmp_path), *a], check=True, capture_output=True)

    git("init", "-q")
    git("config", "user.email", "t@example.org")
    git("config", "user.name", "t")
    for d in ("data/results", "genomeos", "scripts", "tests", "data/organisms"):
        (tmp_path / d).mkdir(parents=True)
    (tmp_path / "data" / "results" / "a.json").write_text("{}")
    (tmp_path / "genomeos" / "x.py").write_text("X = 1\n")
    git("add", ".")
    git("commit", "-q", "-m", "base")
    return tmp_path


@pytest.mark.parametrize(
    "path",
    ["scripts/new_writer.py", "genomeos/sub/pkg/new.py", "tests/test_helper.py", "data/organisms/new.bio"],
)
def test_untracked_code_makes_the_stamp_dirty_and_is_named(tmp_path, path):
    root = _repo(tmp_path)
    assert mf.code_revision(root)["dirty"] is False
    (root / path).parent.mkdir(parents=True, exist_ok=True)
    (root / path).write_text("print('not committed')\n")
    rev = mf.code_revision(root)
    assert rev["dirty"] is True
    assert rev["untracked_code_paths"] == [path]
    assert rev["dirty_code_paths"] == [path]


def test_a_path_with_a_space_is_named_exactly(tmp_path):
    root = _repo(tmp_path)
    (root / "scripts" / "my writer.py").write_text("")
    assert mf.code_revision(root)["untracked_code_paths"] == ["scripts/my writer.py"]


def test_untracked_data_does_not_make_the_stamp_dirty(tmp_path):
    root = _repo(tmp_path)
    (root / "data" / "table.tsv").write_text("a\tb\n")
    (root / "notes.md").write_text("draft\n")
    (root / "data" / "results" / "b.json").write_text("{}")
    rev = mf.code_revision(root)
    assert rev["dirty"] is False and rev["untracked_code_paths"] == [] and rev["dirty_code_paths"] == []
    assert rev["dirty_result_paths"] == ["data/results/b.json"]  # an earlier writer's output, apart


def test_a_renamed_code_file_is_named_once_by_its_new_path(tmp_path):
    root = _repo(tmp_path)
    subprocess.run(["git", "-C", str(root), "mv", "genomeos/x.py", "genomeos/y.py"], check=True)
    rev = mf.code_revision(root)
    assert rev["dirty"] is True and rev["dirty_code_paths"] == ["genomeos/y.py"]
    assert rev["untracked_code_paths"] == []


def test_a_result_written_by_an_untracked_script_records_it(tmp_path, monkeypatch):
    root = _repo(tmp_path)
    (root / "scripts" / "writer.py").write_text("")
    monkeypatch.chdir(root)
    p = save_result("demo", {}, root / "out", manifest=full(tmp_path), strict=True)
    code = json.loads(p.read_text())[mf.KEY]["code"]
    assert code["dirty"] is True and "scripts/writer.py" in code["untracked_code_paths"]
