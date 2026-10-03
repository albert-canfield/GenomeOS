# SPDX-License-Identifier: AGPL-3.0-or-later
"""A refused registry write is quarantined OUTSIDE the registry, for every input (2026-10-03).

`genomeos/results.py: quarantine_dir` used to read

    return results_dir.parent / QUARANTINE / results_dir.name

which is right for the registry itself -- `data/results` -> `data/quarantine/results`, what the
docstring and docs/DATA.md describe -- and wrong for anything BELOW it. For a SUBdirectory of the
registry the parent is the registry, so the quarantine landed at

    data/results/newsub   ->   data/results/quarantine/newsub

INSIDE the tracked registry it exists to keep failed numbers out of, under a `.gitignore` of its own
so that `git status` would never mention it. Measured, not reasoned about: the mapping above is this
file's twin assertion against a copy of the module carrying the old one-liner.

WHY IT BECAME REACHABLE. `d33146a` (lane-registryid) made `_is_registry` answer by FILE and walk up
to the nearest existing ancestor, so `data/results/newsub` is a registry write although it does not
exist, and the refusal path -- which is what builds the quarantine path -- runs for subdirectories
for the first time. `git ls-files data/results` and `find data/results -mindepth 1 -type d` both
report no subdirectory today, so this is a closed door rather than a live leak, which is when to fix
it. The finding is named in tests/test_results_registry_identity.py (plant (d)) by the lane that
found it and routed it here.

THE FIX IS (a) OF THE TWO OFFERED: put the quarantine outside the registry for every input, by
splitting `results_dir` at the registry and MIRRORING the part below it under `data/quarantine/`,
rather than (b) refusing to quarantine at all below the registry. The quarantine exists so that
"hours of compute are kept and the run still fails" (the `save_result` docstring); (b) would have
thrown the payload away for exactly the writes that are new enough to be aiming at a subdirectory,
and would have given `quarantine_dir` -- which `scripts/manifest_census.py` calls for a path it does
not control -- a way to raise.

Each plant has its verdict on BOTH bodies in one run: `old_body()` is a copy of the module with the
one-liner restored, because a test that passes before and after has tested nothing.
"""

import importlib.util
import json
import re
import tempfile
from pathlib import Path

import pytest
import tracked_paths as tp

from genomeos import manifest as mf
from genomeos import results

# --- the two bodies, read out of the module so a future edit fails here rather than vacuously ------

#: The split-at-the-registry mapping as `genomeos/results.py` now writes it.
NEW_MAPPING = """    anchor, below = _registry_anchor(results_dir)
    out = anchor.parent / QUARANTINE / anchor.name
    for part in below:
        out = out / part
    return out
"""

#: The body this lane replaced: the parent of whatever it was handed.
OLD_MAPPING = """    return results_dir.parent / QUARANTINE / results_dir.name
"""

#: Where the `.gitignore` goes. The quarantine is one level deep no longer, so `q.parent.parent` is
#: no longer the quarantine root; both spellings are pinned for the same reason as the mapping.
NEW_IGNORE_ROOT = "    root = _quarantine_root(q)\n"
OLD_IGNORE_ROOT = "    root = q.parent.parent\n"

_OLD: list = []


def old_body():
    """A copy of `genomeos/results.py` whose quarantine mapping is the one-liner again.

    Patched from the file's own source, so if either body stops being spelled as pinned above the
    assertion fires instead of the twins passing for the wrong reason."""
    if _OLD:
        return _OLD[0]
    whole = Path(results.__file__).read_text()
    for new in (NEW_MAPPING, NEW_IGNORE_ROOT):
        assert whole.count(new) == 1, (
            f"no longer in genomeos/results.py as written here, so the twins cannot be built: {new!r}"
        )
    whole = whole.replace(NEW_MAPPING, OLD_MAPPING, 1).replace(NEW_IGNORE_ROOT, OLD_IGNORE_ROOT, 1)
    out = Path(tempfile.mkdtemp(prefix="quarantine-")) / "old_results.py"
    out.write_text(whole)
    spec = importlib.util.spec_from_file_location("old_body_quarantine_results", out)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    _OLD.append(module)
    return module


# --- the scaffolding, the same shapes tests/test_results_registry_identity.py plants --------------


def contract(tmp_path):
    """A manifest that meets every part of the contract except the cleanliness block, so a registry
    write of a name that is not on the allowlist is refused and therefore quarantined."""
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
    allow = tmp_path / "data" / "results_legacy.txt"
    allow.parent.mkdir(parents=True, exist_ok=True)
    allow.write_text("# a header line\nold_result\ttracked\n")
    return allow


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
    """The same module with the one-liner, its globals restored after the test. Plain assignment:
    `monkeypatch` cannot be asked to restore a module it does not know about."""
    module = old_body()
    was = (module.RESULTS_DIR, module.LEGACY_ALLOWLIST)
    yield module
    module.RESULTS_DIR, module.LEGACY_ALLOWLIST = was


def point(module, tmp_path, registry_at) -> None:
    module.RESULTS_DIR = registry_at
    module.LEGACY_ALLOWLIST = allowlist(tmp_path)


# --- plant 1: a refused write below the registry lands nowhere under it ---------------------------


def test_a_refused_write_to_a_registry_subdirectory_lands_outside_the_registry(tree, tmp_path):
    """The defect's own case. Nothing the refusal creates may be anywhere under the registry, and
    the payload must still be kept: the quarantine mirrors the part below the registry, so
    `<reg>/newsub` is quarantined at `<reg>/../quarantine/results/newsub`."""
    sub = tree / "newsub"
    assert results._is_registry(sub) is True, "the premise: d33146a reads this as a registry write"

    q = results.quarantine_dir(sub)
    assert tree not in q.parents and q != tree, f"the quarantine is inside the registry: {q}"
    assert q == tmp_path / "data" / "quarantine" / "results" / "newsub"

    with pytest.raises(mf.ManifestError) as e:
        results.save_result("brand_new", {"value": 7}, sub, manifest=contract(tmp_path))
    assert str(q / "brand_new.json") in str(e.value), "the refusal names where the bytes went"

    assert list(tree.rglob("*")) == [], f"the registry gained files: {list(tree.rglob('*'))}"
    kept = json.loads((q / "brand_new.json").read_text())
    assert kept["value"] == 7, "hours of compute are kept, which is what the quarantine is for"
    assert kept["quarantine"]["meant_for"] == str(sub / "brand_new.json")
    root = tmp_path / "data" / "quarantine"
    assert (root / ".gitignore").read_text() == results._QUARANTINE_IGNORE, (
        "and the ignore file is at the root of the quarantine, which is the only place it covers the "
        "whole of it"
    )
    assert {p.name for p in tree.parent.iterdir()} == {"results", "quarantine", "results_legacy.txt"}


def test_the_old_body_quarantined_it_inside_the_registry(old, tree, tmp_path):
    """The twin, and the measurement the top of this file reports: the one-liner put the bytes at
    `<reg>/quarantine/newsub/`, inside the tracked registry, under a `.gitignore` of its own."""
    point(old, tmp_path, tree)
    sub = tree / "newsub"
    assert old._is_registry(sub) is True

    q = old.quarantine_dir(sub)
    assert q == tree / "quarantine" / "newsub", "the old mapping, for the record"
    assert tree in q.parents, "inside the registry it was protecting"

    with pytest.raises(mf.ManifestError):
        old.save_result("brand_new", {"value": 7}, sub, manifest=contract(tmp_path))

    landed = sorted(p.relative_to(tree).as_posix() for p in tree.rglob("*"))
    assert landed == [
        "quarantine",
        "quarantine/.gitignore",
        "quarantine/newsub",
        "quarantine/newsub/brand_new.json",
    ], f"the registry gained a quarantine of its own: {landed}"
    assert json.loads((q / "brand_new.json").read_text())["value"] == 7


# --- plant 2: the registry's own path is unchanged, byte for byte ---------------------------------


def test_the_registrys_own_quarantine_path_is_unchanged(tree, tmp_path):
    """`data/quarantine/results` is depended upon (docs/DATA.md, `scripts/manifest_census.py`,
    `data/quarantine/results/` on this disk), so the registry's own answer may not move by a
    character. Asserted against the literal, on both bodies, and for the real `RESULTS_DIR`
    spelling as well as the planted one."""
    assert results.quarantine_dir(tree) == tmp_path / "data" / "quarantine" / "results"
    assert str(results.quarantine_dir(tree)) == str(old_body().quarantine_dir(tree)), (
        "the registry's own answer is the one-liner's answer, character for character"
    )


def test_the_real_registry_spelling_still_answers_data_quarantine_results(monkeypatch):
    """Unpatched, from the repository root: the documented answer, as a string. `chdir` rather than
    the ambient directory so the relative spelling is what it is in the docstring."""
    root = Path(__file__).resolve().parent.parent
    # `data/results` holds hundreds of committed files, so git TRACKS the directory: "no
    # data/results" is not a state a checkout of this commit can be in and the skip never fired
    tp.must_be_present(root / "data" / "results", was="no data/results")
    monkeypatch.chdir(root)
    assert str(results.quarantine_dir(Path("data/results"))) == "data/quarantine/results"
    assert str(results.quarantine_dir(results.RESULTS_DIR)) == "data/quarantine/results"


def test_a_refused_registry_write_is_quarantined_byte_for_byte_as_before(old, tree, tmp_path):
    """Not just the path: the FILE. The same refusal under both bodies, compared character for
    character with only the three run-dependent fields normalised -- `written_at`, and the two
    open-counters of `traced_inputs`, which count the enforcement machinery's own reads and so differ
    between the first run of a process and the second whatever the body is (measured: 8 against 12)."""
    m = contract(tmp_path)
    q = results.quarantine_dir(tree) / "brand_new.json"
    with pytest.raises(mf.ManifestError):
        results.save_result("brand_new", {"value": 7}, tree, manifest=m)
    now = q.read_text()
    q.unlink()

    point(old, tmp_path, tree)
    q_old = old.quarantine_dir(tree) / "brand_new.json"
    assert q_old == q, "the same path"
    with pytest.raises(mf.ManifestError):
        old.save_result("brand_new", {"value": 7}, tree, manifest=m)
    before = q_old.read_text()

    def normalised(text: str) -> str:
        text = re.sub(r'"written_at": "[^"]*"', '"written_at": "T"', text)
        return re.sub(r'"(opens|opens_under_data)": \d+', r'"\1": N', text)

    assert normalised(now) == normalised(before), "the quarantined file is what it was"
    assert '"value": 7' in now and '"meant_for"' in now, "and it is not vacuously empty"


# --- plant 3: the near-miss, a scratch directory, quarantines where it does now -------------------


def test_a_scratch_directory_still_quarantines_beside_itself(tmp_path):
    """`tmp_path/"results"` is how the whole suite writes and is NOT the registry (plant (c) of
    tests/test_results_registry_identity.py). Its quarantine must stay beside it, on both bodies:
    this is the near-miss that a fix anchored on `RESULTS_DIR` alone would have moved into
    `data/quarantine/`."""
    out = tmp_path / "results"
    assert results._is_registry(out) is False, "the premise"

    assert results.quarantine_dir(out) == tmp_path / "quarantine" / "results"
    assert old_body().quarantine_dir(out) == tmp_path / "quarantine" / "results"

    m = contract(tmp_path)
    del m["partitions"]  # incomplete, so strict=True refuses outside the registry
    with pytest.raises(mf.ManifestError):
        results.save_result("scratch_thing", {"value": 1}, out, manifest=m, strict=True)
    q = tmp_path / "quarantine" / "results" / "scratch_thing.json"
    assert json.loads(q.read_text())["value"] == 1
    assert (tmp_path / "quarantine" / ".gitignore").exists(), "the root ignore, where it always was"


def test_a_scratch_subdirectory_outside_the_registry_is_untouched(tmp_path):
    """Nothing below a non-registry directory is mirrored anywhere: the mapping only splits at the
    registry, so a scratch tree keeps its own shape."""
    deep = tmp_path / "results" / "newsub"
    assert results.quarantine_dir(deep) == tmp_path / "results" / "quarantine" / "newsub"
    assert results.quarantine_dir(deep) == old_body().quarantine_dir(deep)


# --- the mapping itself, for the shapes a writer can hand it --------------------------------------


def test_the_mapping_splits_at_the_registry_however_the_registry_is_spelled(tree, tmp_path):
    """A case variant and a `./` detour name the same directory, so `samefile` answers for them and
    the split happens at the ANCESTOR AS SPELLED -- which is what keeps the registry's own answer
    byte-identical for every spelling of it."""
    quarantine = tmp_path / "data" / "quarantine"
    assert results.quarantine_dir(tree / "a" / "b") == quarantine / "results" / "a" / "b"

    detour = tree.parent / "." / "results"
    assert results.quarantine_dir(detour) == quarantine / "results", "a './' detour is the registry"
    assert results.quarantine_dir(detour / "newsub") == quarantine / "results" / "newsub"


def test_an_unanswerable_ancestor_is_climbed_past_rather_than_stopping_the_split(tree, tmp_path):
    """A directory that cannot be stat'ed between the registry and the write does not hide the
    registry: the walk keeps climbing, so the split still happens at the registry and the quarantine
    is still outside it."""
    blocked = tree / "noperm"
    blocked.mkdir()
    (blocked / "sub").mkdir()
    blocked.chmod(0o000)
    try:
        q = results.quarantine_dir(blocked / "sub")
        assert q == tmp_path / "data" / "quarantine" / "results" / "noperm" / "sub"
        assert tree not in q.parents
    finally:
        blocked.chmod(0o755)
