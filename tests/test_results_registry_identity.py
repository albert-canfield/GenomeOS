# SPDX-License-Identifier: AGPL-3.0-or-later
"""`genomeos/results.py: _is_registry` decides identity by FILE, not by spelling (2026-10-02).

That one flag gates FOUR decisions in `save_result` -- `legacy`, `enforce` (where a registry write
stops being enforced and becomes a function of whatever the caller passed as `strict`), the
requirement to carry a cleanliness block at all, and the entry-script refusal -- so a wrong answer
turns four refusals off together. It used to read

    try:
        return results_dir.resolve() == RESULTS_DIR.resolve()
    except OSError:
        return True

and that compares SPELLINGS. Three holes, each measured on this machine rather than reasoned about:

1. CASE. From the repository root, `Path('data/RESULTS').resolve() != Path('data/results').resolve()`
   while `os.path.samefile('data/RESULTS', 'data/results')` is True: this volume is
   case-INSENSITIVE, so `data/RESULTS` IS the registry and the comparison denied it. No OSError
   needed. On a case-SENSITIVE volume the variant simply does not exist and is correctly not the
   registry -- so the tests below branch on what the filesystem actually is, measured by writing one
   spelling and probing for another, never by platform name.

2. A PATH THAT DOES NOT EXIST YET can still be INSIDE the registry: `data/results/newsub`. Closed
   here as a door rather than a leak -- `git ls-files data/results` and
   `find data/results -mindepth 1 -type d` both report no subdirectory today.

3. The handler itself. `except OSError` was measured UNREACHABLE through the filesystem on
   CPython 3.12/APFS: a symlink loop makes non-strict `Path.resolve()` raise RuntimeError
   ("Symlink loop from ..."), which that clause does not catch, while EACCES on an unsearchable
   parent and ENAMETOOLONG make it raise NOTHING at all, so it returned the path unchanged and the
   spelling comparison answered False. `os.path.samefile` raises OSError 62/13/63 for the same
   three, so the handler now fires where it was always meant to.

Three outcomes, not two. `registry = _is_registry(results_dir)` runs 104 lines BEFORE
`results_dir.mkdir(...)`, and every writer in the suite passes `tmp_path / "results"`, which does not
exist yet. Reading "missing" as the registry would refuse all of them. So: identity by file;
FileNotFoundError walks UP to the nearest existing ancestor and asks whether THAT is, or lies under,
the registry; any other OSError leaves the question UNANSWERED and unanswered counts as the registry,
because a false-loud refusal gets read and a false-silent pass does not.

Each plant has a twin that runs the same scenario through a copy of the module carrying the OLD
guarded return, because a test that passes before and after has tested nothing.
"""

import errno
import importlib.util
import json
import os
import tempfile
from pathlib import Path

import pytest

from genomeos import manifest as mf
from genomeos import results

# --- the two bodies, read out of the module so a future edit fails here rather than vacuously ------

#: The guarded return as `genomeos/results.py` now writes it.
NEW_GUARDED = """    try:
        return os.path.samefile(results_dir, RESULTS_DIR)
    except FileNotFoundError:
        return _nearest_existing_is_registry(results_dir)
    except OSError:
        return True
"""

#: The body this lane replaced: a spelling comparison with one handler.
OLD_GUARDED = """    try:
        return results_dir.resolve() == RESULTS_DIR.resolve()
    except OSError:
        return True
"""

_OLD: list = []


def old_body():
    """A copy of `genomeos/results.py` whose `_is_registry` is the spelling comparison again.

    Patched from the file's own source, so if that body stops being spelled as `NEW_GUARDED` the
    assertion fires instead of the twins passing for the wrong reason."""
    if _OLD:
        return _OLD[0]
    whole = Path(results.__file__).read_text()
    assert whole.count(NEW_GUARDED) == 1, (
        "the guarded return these twins patch is no longer in genomeos/results.py as written here: "
        "update NEW_GUARDED"
    )
    out = Path(tempfile.mkdtemp(prefix="registryid-")) / "old_results.py"
    out.write_text(whole.replace(NEW_GUARDED, OLD_GUARDED, 1))
    spec = importlib.util.spec_from_file_location("old_body_results", out)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    _OLD.append(module)
    return module


# --- the scaffolding -------------------------------------------------------------------------------


def contract(tmp_path):
    """A manifest that meets every part of the contract except the cleanliness block."""
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


def allowlist(tmp_path) -> Path:
    """A one-name legacy allowlist beside the planted registry."""
    allow = tmp_path / "data" / "results_legacy.txt"
    allow.parent.mkdir(parents=True, exist_ok=True)
    allow.write_text("# a header line\nold_result\ttracked\n")
    return allow


def point(module, tmp_path, registry_at) -> None:
    """Point the OLD-body copy's registry at `registry_at`. Plain assignment, because `monkeypatch`
    cannot be asked to restore a module it does not know about; the `old` fixture restores it."""
    module.RESULTS_DIR = registry_at
    module.LEGACY_ALLOWLIST = allowlist(tmp_path)


@pytest.fixture
def tree(tmp_path, monkeypatch):
    """A checkout-shaped registry at tmp/data/results, existing, with its allowlist beside it."""
    reg = tmp_path / "data" / "results"
    reg.mkdir(parents=True)
    monkeypatch.setattr(results, "RESULTS_DIR", reg)
    monkeypatch.setattr(results, "LEGACY_ALLOWLIST", allowlist(tmp_path))
    return reg


@pytest.fixture
def old(tmp_path):
    """The same module with the old guarded return, its globals restored after the test."""
    module = old_body()
    was = (module.RESULTS_DIR, module.LEGACY_ALLOWLIST)
    yield module
    module.RESULTS_DIR, module.LEGACY_ALLOWLIST = was


def case_insensitive(where: Path) -> bool:
    """Whether `where`'s volume is case-insensitive, MEASURED: write one spelling, probe another.

    Not `sys.platform` -- a case-sensitive APFS volume and a case-insensitive one are both darwin,
    and the right assertion for plant (a) is a different one on each."""
    probe = where / "CaseSensitivityProbe"
    probe.write_text("")
    try:
        return (where / "casesensitivityprobe").exists()
    finally:
        probe.unlink()


def raises_errno(path: Path, wanted: int) -> bool:
    """Whether `os.stat(path)` really fails with `wanted` here, so a plant that depends on it skips
    by name instead of asserting something this machine cannot produce (root bypasses EACCES)."""
    try:
        os.stat(path)
    except OSError as e:
        return e.errno == wanted
    return False


# --- (a) a case variant of the registry ------------------------------------------------------------


def test_a_case_variant_is_decided_by_the_file_and_not_by_the_spelling(tree, tmp_path):
    """On a case-INSENSITIVE volume `data/RESULTS` IS `data/results`, so the write must refuse. On a
    case-SENSITIVE one it is a different, non-existent directory whose nearest existing ancestor is
    `data/`, which is not under the registry, so the write is correctly ALLOWED. The test asserts the
    property the filesystem actually has rather than one it has on half the machines."""
    variant = tmp_path / "data" / "RESULTS"

    if case_insensitive(tmp_path / "data"):
        assert os.path.samefile(variant, tree), "the premise: the variant is the same directory"
        assert results._is_registry(variant) is True
        with pytest.raises(mf.ManifestError) as e:
            results.save_result("brand_new", {"value": 7}, variant, manifest=contract(tmp_path))
        assert str(variant) in str(e.value), "the refusal names the path"
        assert not (tree / "brand_new.json").exists(), "nothing entered the registry"
    else:
        assert not variant.exists(), "the premise: on a case-sensitive volume the variant is absent"
        assert results._is_registry(variant) is False
        p = results.save_result("brand_new", {"value": 7}, variant, manifest=contract(tmp_path))
        assert p.exists(), "a different directory is not the registry and the write happens"


def test_the_old_body_got_the_case_variant_wrong_on_a_case_insensitive_volume(old, tree, tmp_path):
    """The twin. A spelling comparison: `data/RESULTS` != `data/results` whatever the volume, so on
    this one the old body let an unenforced write into the registry itself. On a case-sensitive
    volume both bodies answer False and the fix changed nothing there -- stated, not hidden."""
    point(old, tmp_path, tree)
    variant = tmp_path / "data" / "RESULTS"
    insensitive = case_insensitive(tmp_path / "data")

    assert old._is_registry(variant) is False, (
        "the old body compared resolved spellings, which differ whatever the volume"
    )
    p = old.save_result("brand_new", {"value": 7}, variant, manifest=contract(tmp_path))
    if insensitive:
        assert p.exists() and os.path.samefile(p.parent, tree), (
            "and so an unenforced write landed in the registry itself"
        )
    else:
        assert not os.path.samefile(p.parent, tree), (
            "on a case-sensitive volume the variant is genuinely elsewhere and both bodies agree"
        )


# --- (b) a path the question cannot be answered for ------------------------------------------------


def test_a_symlink_loop_refuses_with_the_path_named_and_the_payload_kept(tree, tmp_path):
    """ELOOP is not FileNotFoundError: the question is unanswered, unanswered counts as the registry,
    and the refusal takes the path that was already there -- quarantined with its reason,
    ManifestError naming where it went, nothing written."""
    loop = tmp_path / "data" / "loop"
    loop.symlink_to("loop")
    assert raises_errno(loop, errno.ELOOP), "the premise: this path really raises ELOOP here"

    assert results._is_registry(loop) is True
    with pytest.raises(mf.ManifestError) as e:
        results.save_result("brand_new", {"value": 7}, loop, manifest=contract(tmp_path))

    assert str(loop) in str(e.value), "the refusal names the path it could not classify"
    kept = json.loads((results.quarantine_dir(loop) / "brand_new.json").read_text())
    assert kept["value"] == 7, "the payload is kept, so hours of compute are not lost"


def test_the_old_body_did_not_even_raise_oserror_on_a_symlink_loop(old, tree, tmp_path):
    """The twin, and the sharpest of the three measurements. The old handler was written for this
    case -- its docstring named ELOOP -- but CPython 3.12's non-strict `Path.resolve` converts ELOOP
    into RuntimeError, which `except OSError` does not catch. So the write was neither refused nor
    recorded: it raised an error of a kind no caller is told about, and nothing was quarantined."""
    point(old, tmp_path, tree)
    loop = tmp_path / "data" / "loop"
    loop.symlink_to("loop")

    with pytest.raises(RuntimeError, match="Symlink loop"):
        old._is_registry(loop)
    assert not (old.quarantine_dir(loop) / "brand_new.json").exists()


@pytest.mark.parametrize("err", [errno.EACCES, errno.ENAMETOOLONG])
def test_the_other_unanswerable_paths_count_as_the_registry_where_the_old_body_read_them_as_not(
    old, tree, tmp_path, err
):
    """EACCES on an unsearchable parent, and ENAMETOOLONG. `os.path.samefile` raises OSError for
    both, so the write refuses; non-strict `Path.resolve` raises NOTHING for either and handed back
    the path unchanged, so the old body fell through to the spelling comparison and answered False --
    not the registry, four refusals off. Asserted on `_is_registry` rather than through
    `save_result`, because the quarantine for such a path sits under the same unusable parent."""
    point(old, tmp_path, tree)
    if err == errno.EACCES:
        locked = tmp_path / "locked"
        (locked / "results").mkdir(parents=True)
        path = locked / "results"
        locked.chmod(0o000)
    else:
        path = tmp_path / "data" / ("x" * 300)

    try:
        if not raises_errno(path, err):
            pytest.skip(f"this machine does not produce {errno.errorcode[err]} for {path}")
        assert results._is_registry(path) is True, "unanswered is treated as the registry"
        assert old._is_registry(path) is False, "the old body answered 'not the registry'"
    finally:
        if err == errno.EACCES:
            (tmp_path / "locked").chmod(0o755)


# --- (c) the near miss, which must keep passing ----------------------------------------------------


def test_the_near_miss_a_scratch_directory_not_yet_created_records_the_false_verdict_and_writes():
    """How the whole suite writes: `tmp_path / "results"`, which does not exist when `_is_registry`
    runs, 104 lines before the `mkdir`. Its nearest existing ancestor is `tmp_path`, which is not
    under the registry, so it is NOT a registry write: the false verdict is recorded and the bytes
    are written. This passed on the old body too -- it is the regression guard that makes the fix a
    fix rather than a refusal of everything.

    The real `RESULTS_DIR` is left in place on purpose: this is the suite's own case, not a planted
    one."""
    with tempfile.TemporaryDirectory() as d:
        tmp_path = Path(d)
        out = tmp_path / "results"
        assert not out.exists(), "the premise: it does not exist when the flag is decided"
        assert results._is_registry(out) is False
        assert old_body()._is_registry(out) is False, "and the old body agreed, as it had to"

        m = contract(tmp_path)
        m["code_cleanliness"] = mf.code_cleanliness(__file__, ())
        p = results.save_result("anything", {"value": 1}, out, manifest=m)

        assert p.exists() and json.loads(p.read_text())["value"] == 1
        block = json.loads(p.read_text())[mf.KEY]["code_cleanliness"]
        assert block["own_code_is_committed"] is False, (
            "the false verdict is RECORDED, not hidden: under pytest the entry script is the runner"
        )


# --- (d) a subdirectory of the registry that does not exist yet ------------------------------------


def test_a_subdirectory_of_the_registry_that_does_not_exist_yet_is_a_registry_write(tree, tmp_path):
    """Missing is not the same as outside. The nearest existing ancestor of `data/results/newsub` is
    `data/results` itself, so the write is a registry write and refuses -- and the `mkdir` 104 lines
    down never runs, so no result directory is created under the tracked registry.

    What the refusal DOES create, and it is a finding rather than this lane's doing: `quarantine_dir`
    is `results_dir.parent / "quarantine" / results_dir.name`, so the quarantine for a registry
    SUBdirectory lands at `data/results/quarantine/newsub/` -- inside the registry, under a
    `.gitignore` of its own. Named here instead of asserted away."""
    sub = tree / "newsub"
    assert not sub.exists()

    assert results._is_registry(sub) is True
    with pytest.raises(mf.ManifestError) as e:
        results.save_result("brand_new", {"value": 7}, sub, manifest=contract(tmp_path))

    assert str(sub) in str(e.value), "the refusal names the path"
    assert not sub.exists(), "no result directory is created inside the registry"
    assert list(tree.glob("*.json")) == [], "and no result entered the registry"


def test_the_old_body_wrote_into_a_new_registry_subdirectory_unenforced(old, tree, tmp_path):
    """The twin: a different resolved spelling, so not the registry, so `enforce` fell back to the
    caller's `strict` and the bytes landed under the tracked registry with no cleanliness block and
    no entry-script refusal."""
    point(old, tmp_path, tree)
    sub = tree / "newsub"

    assert old._is_registry(sub) is False
    p = old.save_result("brand_new", {"value": 7}, sub, manifest=contract(tmp_path))
    assert p.exists() and p.parent.parent == tree, "written inside the registry, unenforced"
    stamped = json.loads(p.read_text())[mf.KEY]
    assert "code_cleanliness" not in stamped, (
        "and carrying no cleanliness block at all, while reading `complete: True` -- the four "
        f"refusals were off, so nothing asked for one (complete: {stamped['complete']})"
    )


# --- (e) the registry itself missing ---------------------------------------------------------------


def test_a_registry_that_is_not_there_cannot_be_answered_for_and_refuses(old, tmp_path, monkeypatch):
    """If RESULTS_DIR itself is the missing one, no path can be placed relative to it and every
    answer would be a guess. Unanswerable refuses.

    The spelling is deliberately a scratch directory that plainly exists: handed the registry's own
    spelling the old body answered True by accident, the two strings being equal, and the plant would
    have passed whatever the body did."""
    missing = tmp_path / "data" / "results"  # never created
    monkeypatch.setattr(results, "RESULTS_DIR", missing)
    monkeypatch.setattr(results, "LEGACY_ALLOWLIST", allowlist(tmp_path))
    point(old, tmp_path, missing)
    scratch = tmp_path / "scratch"
    scratch.mkdir()

    assert not missing.exists(), "the premise"
    assert results._is_registry(scratch) is True
    assert results._is_registry(missing) is True, "and the registry's own spelling too"
    assert old._is_registry(scratch) is False, "the twin: the old body answered for it anyway"

    with pytest.raises(mf.ManifestError) as e:
        results.save_result("brand_new", {"value": 7}, scratch, manifest=contract(tmp_path))
    assert str(scratch) in str(e.value), "the refusal names the path"
    assert not (scratch / "brand_new.json").exists()


def test_a_registry_whose_own_identity_cannot_be_read_refuses(tmp_path, monkeypatch):
    """The other half of the walk's first step: the path is merely missing, but the registry it would
    have to be compared against cannot be stat-ed at all. Unanswered, so it refuses -- rather than
    reading the error as 'the registry is not there, therefore this is not a registry write'."""
    locked = tmp_path / "locked"
    (locked / "results").mkdir(parents=True)
    monkeypatch.setattr(results, "RESULTS_DIR", locked / "results")
    monkeypatch.setattr(results, "LEGACY_ALLOWLIST", allowlist(tmp_path))
    locked.chmod(0o000)
    try:
        if not raises_errno(locked / "results", errno.EACCES):
            pytest.skip("this machine does not produce EACCES here")
        assert results._is_registry(tmp_path / "scratch" / "results") is True
    finally:
        locked.chmod(0o755)


# --- the ancestor walk itself ----------------------------------------------------------------------


def test_the_walk_stops_at_the_nearest_existing_ancestor_and_reads_it(tree, tmp_path):
    """Several levels missing at once: `data/results/a/b/c` is still inside the registry, and
    `tmp/x/y/z` is still outside it, because the walk climbs until something exists and asks its
    question about that."""
    assert results._is_registry(tree / "a" / "b" / "c") is True
    assert results._is_registry(tmp_path / "x" / "y" / "z") is False


def test_the_walk_reads_the_ancestor_by_file_too_so_a_case_variant_of_it_counts(tree, tmp_path):
    """The two holes together: a missing subdirectory of a CASE VARIANT of the registry. Neither the
    child nor the variant can be compared by spelling, and on a case-insensitive volume both are the
    registry."""
    target = tmp_path / "data" / "RESULTS" / "newsub"
    expected = case_insensitive(tmp_path / "data")
    assert results._is_registry(target) is expected
