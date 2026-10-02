# SPDX-License-Identifier: AGPL-3.0-or-later
"""The one framework read `manifest_rebuild` exempts, and the condition that makes the exemption sound.

`data/results_legacy.txt` (`genomeos/results.py: LEGACY_ALLOWLIST`) is read by `save_result` itself on
every registry write, so once a rebuild reconciles a run's recorded reads against its declared inputs,
every post-`fe0880a` result reads a file it does not declare -- some 1,100 of them at once. It is
therefore exempt by EXACT PATH (`scripts/manifest_rebuild.py: FRAMEWORK_READS`).

AN EXEMPTION IS ONLY AS GOOD AS ITS CONDITION, and this one's condition is that the allowlist can never
change a result's VALUES -- only whether `save_result` ADMITS it. That is a claim about behaviour, so it
is PROVEN here rather than argued: the allowlist is changed under a temporary registry and a writer is
re-run, and the result's bytes come back IDENTICAL, or the write is QUARANTINED and nothing reaches the
registry. Never written with different content.

If any value could move with the allowlist, the file would be an INPUT and would have to be declared
with a sha256 like any other, and the exemption would be wrong. These tests are what would say so.

What the exemption gives up, which is not nothing and is written down where it is taken: a rebuild is
blind to a result that was admitted because a name was ADDED to that list. That is a change in what the
registry accepts, not in what a result says.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from typing import Any

import planted_cleanliness
import pytest

from genomeos import manifest as mf
from genomeos import results
from genomeos.results import save_result

_spec = importlib.util.spec_from_file_location("mr3", Path("scripts/manifest_rebuild.py"))
mr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mr)

NAME = "planted_exemption_probe"


@pytest.fixture
def registry(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    """A registry at tmp/data/results with its own allowlist beside it, as the real pair sit."""
    reg = tmp_path / "data" / "results"
    reg.mkdir(parents=True)
    allow = tmp_path / "data" / "results_legacy.txt"
    allow.write_text("# a header line\n")
    monkeypatch.setattr(results, "RESULTS_DIR", reg)
    monkeypatch.setattr(results, "LEGACY_ALLOWLIST", allow)
    results.legacy_names.cache_clear() if hasattr(results.legacy_names, "cache_clear") else None
    return {"reg": reg, "allow": allow}


def _complete(tmp_path: Path) -> dict[str, Any]:
    """A manifest that meets the whole contract, cleanliness block included, so the only thing that can
    change the outcome below is the allowlist."""
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
        # 2026-10-02: the entry script is counted now (genomeos.manifest.entry_script), so under pytest
        # `sys.argv[0]` is a test runner and no in-process call can return a passing block. This
        # scaffolding block is built in a planted repository whose entry script really is committed and
        # clean, which is what it always meant to assert, instead of naming a file of this shared
        # checkout and reddening when a peer edits it (tests/planted_cleanliness.py).
        "code_cleanliness": planted_cleanliness.planted_clean_block(tmp_path, "allowlist_clean"),
    }


def _write(registry: dict[str, Path], tmp_path: Path) -> dict[str, Any]:
    p = save_result(NAME, {"value": 7, "derived": [1, 2, 3]}, registry["reg"], manifest=_complete(tmp_path))
    return json.loads(p.read_text())


# --- the condition: the allowlist moves admission, never a value -----------------------------------


def test_changing_the_allowlist_moves_no_value_of_an_admitted_result(
    registry: dict[str, Path], tmp_path: Path
) -> None:
    """The exemption's whole justification, RUN rather than asserted, and it was run before it was
    written: the allowlist was changed under a temporary registry, the writer re-run, and the two
    results compared leaf by leaf. Exactly ONE leaf moved and it is not a value.

    `/result_manifest/traced_inputs/opens` went 13 to 16, which is the tracer's own total count of
    `open` calls -- already classified a RUN FIELD by this project, by exact path, with what exempting
    it costs written down in `scripts/manifest_rebuild.py: RUN_FIELDS`. It moved here because the two
    writes happen in one process in sequence, so the second write's tracer window held the first
    write's reads; it did not move because of the allowlist.

    The assertion is anchored to the project's OWN exact-path classification rather than to a carve-out
    made here, so a leaf moving anywhere else -- including a value -- fails this test.
    """
    before = _write(registry, tmp_path)
    registry["allow"].write_text(f"# a header line\n{NAME}\ttracked\nsome_other_name\ttracked\n")
    after = _write(registry, tmp_path)
    found = mr.diff(before, after)
    rest, run = mr.run_differences(found)
    assert rest == [], (
        "a value moved with the legacy allowlist: it is an INPUT, not a framework read, and the "
        f"exemption in manifest_rebuild.FRAMEWORK_READS is wrong. Moved: {rest}"
    )
    assert before["value"] == after["value"] == 7
    assert before["derived"] == after["derived"] == [1, 2, 3]
    assert len(run) <= 1, run
    for d in run:
        assert mr._base_path(d) in mr.RUN_FIELDS


def test_removing_a_name_from_the_allowlist_quarantines_rather_than_writing_different_content(
    registry: dict[str, Path], tmp_path: Path
) -> None:
    """The other half of "identical or quarantined". An incomplete manifest under a listed name warns and
    writes; the same manifest under an unlisted name is REFUSED, and nothing reaches the registry -- so
    the allowlist can stop a write, and cannot change what a write says."""
    thin = {"sources": [{"accession": "ENCSR000XXX", "version": "v1"}]}
    registry["allow"].write_text(f"# a header line\n{NAME}\ttracked\n")
    with pytest.warns(results.ManifestWarning):
        listed = save_result(NAME, {"value": 7}, registry["reg"], manifest=dict(thin)).read_bytes()
    assert b'"value": 7' in listed

    registry["allow"].write_text("# a header line\n")
    with pytest.raises(mf.ManifestError, match="quarantined at"):
        save_result(NAME, {"value": 7}, registry["reg"], manifest=dict(thin))
    # the earlier admitted bytes are still there and unchanged: the refusal wrote nothing over them
    assert (registry["reg"] / f"{NAME}.json").read_bytes() == listed
    q = results.quarantine_dir(registry["reg"]) / f"{NAME}.json"
    assert q.exists(), "the refused payload is kept, so hours of compute are not lost"
    assert json.loads(q.read_text())["value"] == 7


def test_the_allowlist_is_read_by_save_result_itself(registry: dict[str, Path], tmp_path: Path) -> None:
    """Why it is a FRAMEWORK read: no writer asks for it, `save_result` does, on every registry write."""
    assert registry["allow"] == results.LEGACY_ALLOWLIST
    assert results.legacy_names() == frozenset()
    registry["allow"].write_text(f"# a header line\n{NAME}\ttracked\n")
    assert NAME in results.legacy_names()


# --- the exemption is an exact path, and the loss is named where it is taken -------------------------


def test_the_exemption_is_an_exact_path_and_not_a_pattern() -> None:
    assert mr.FRAMEWORK_READS == ("data/results_legacy.txt",)
    for p in mr.FRAMEWORK_READS:
        assert not p.endswith("/") and "*" not in p


def test_a_result_under_the_same_directory_is_not_exempted() -> None:
    """A prefix or a directory match would have swept in every committed result. It must not."""
    got = mr.reads_against_declared(
        Path(),
        {"inputs": []},
        ["data/results_legacy.txt", "data/results/unknown_chr21.json"],
        [],
    )
    assert got["framework_reads_ignored"] == ["data/results_legacy.txt"]
    assert got["read_but_not_declared"] == ["data/results/unknown_chr21.json"]
    assert got["every_read_is_declared"] is False


def test_the_measured_blindness_is_what_the_reconciliation_now_catches() -> None:
    """The live case: 193 of 193 declared inputs hashed and a 194th file read, pinned by no sha256."""
    declared = [
        {"path": f"data/results/declared_{i}.json", "sha256": "x", "partition": None} for i in range(193)
    ]
    reads = [i["path"] for i in declared] + ["data/results/unknown_chr21.json"]
    got = mr.reads_against_declared(Path(), {"inputs": declared}, reads, [])
    assert got["reads_under_data_count"] == 194
    assert got["read_but_not_declared"] == ["data/results/unknown_chr21.json"]


def test_a_read_that_is_declared_is_not_reported() -> None:
    declared = [{"path": "data/cache/x.tsv", "sha256": "x", "partition": None}]
    got = mr.reads_against_declared(Path(), {"inputs": declared}, ["data/cache/x.tsv"], [])
    assert got["read_but_not_declared"] == []
    assert got["every_read_is_declared"] is True


def test_a_read_outside_data_is_not_a_result_s_input() -> None:
    """The interpreter, the virtual environment and the standard library are not inputs."""
    got = mr.reads_against_declared(Path(), {"inputs": []}, ["/usr/lib/python3.12/json/__init__.py"], [])
    assert got["reads_under_data"] == []
