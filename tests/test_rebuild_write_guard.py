# SPDX-License-Identifier: AGPL-3.0-or-later
"""The rebuild may not write to the stores it verifies a result against (2026-10-02).

The defect, measured on this machine before the guard was written: `scripts/manifest_rebuild.py` SYMLINKS
`data/reference`, `data/knowledge` and `data/cache` into its worktree, and a symlink carries the target's
mode. `open(<worktree>/data/reference/HG002_chr1.vcf.gz, "r+b")` through the link SUCCEEDED, the real
store's directories are mode 0755, and a file created through the link is created in the one copy on this
machine. A verification tool that can modify the inputs it compares against could produce a clean verdict
by changing the thing it compared to.

The two probes that found it are the two cases planted here -- an open-for-update of an existing store
file, and the creation of a new file in the store -- and each is paired with its counterfactual: with the
hook disarmed the same call SUCCEEDS. So removing the hook does not make these tests pass; it makes the
counterfactual assertions the real ones and the refusals fail.

The guard is deliberately NARROWER than "no write outside the worktree": a rebuild legitimately writes
temporary files and package caches outside its worktree, and refusing those would break rebuilds rather
than protect stores. Both halves of that narrowing are tested.
"""

from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from genomeos import manifest as mf

_spec = importlib.util.spec_from_file_location("rwg", Path("scripts/rebuild_write_guard.py"))
guard = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(guard)


@pytest.fixture
def planted(tmp_path: Path) -> dict[str, Path]:
    """A real store with a file in it, and a worktree that reaches it only through a symlink."""
    store = tmp_path / "real" / "data" / "reference"
    store.mkdir(parents=True)
    (store / "HG002_chr1.vcf.gz").write_bytes(b"original bytes\n")
    wt = tmp_path / "worktree"
    (wt / "data").mkdir(parents=True)
    (wt / "data" / "reference").symlink_to(store)
    (wt / "data" / "results").mkdir()
    return {"store": store, "wt": wt, "link": wt / "data" / "reference"}


@pytest.fixture
def armed(planted: dict[str, Path]):
    """The guard armed over the planted store, and emptied again however the test ends."""
    guard.arm([str(planted["store"])])
    try:
        yield planted
    finally:
        guard.disarm()


# --- the write test is pinned to the tracer's, so the two copies cannot drift --------------------


WRITE_CASES = [
    ("r", None),
    ("rb", None),
    ("w", None),
    ("wb", None),
    ("a", None),
    ("x", None),
    ("r+", None),
    ("r+b", None),
    ("rt", None),
    (None, os.O_RDONLY),
    (None, os.O_WRONLY),
    (None, os.O_RDWR),
    (None, os.O_RDONLY | os.O_CREAT),
    (None, os.O_WRONLY | os.O_APPEND),
    (None, os.O_RDONLY | os.O_CLOEXEC),
    (None, None),
]


@pytest.mark.parametrize(("mode", "flags"), WRITE_CASES)
def test_the_guard_and_the_tracer_agree_on_what_counts_as_a_write(mode, flags) -> None:
    """Two copies of a predicate that nothing compared is how ten copies of the cleanliness block drifted
    apart. The guard cannot import `genomeos` (the worktree's copy is at an older revision), so the two
    are pinned here instead."""
    assert guard.is_write(mode, flags) is mf._is_write(mode, flags), (mode, flags)


# --- probe one: an open-for-update of an existing store file, through the symlink ------------------


def test_the_probe_that_found_the_defect_is_refused(armed: dict[str, Path]) -> None:
    target = armed["link"] / "HG002_chr1.vcf.gz"
    with pytest.raises(PermissionError, match="may not write to the store"), open(target, "r+b"):
        pass


def test_the_same_probe_succeeds_with_the_hook_disarmed(planted: dict[str, Path]) -> None:
    """The counterfactual: this is what every rebuild before 2026-10-02 could do."""
    guard.disarm()
    with open(planted["link"] / "HG002_chr1.vcf.gz", "r+b") as f:
        assert f.readable() and f.writable()


# --- probe two: a file created in the real store, through the symlink -----------------------------


def test_creating_a_file_in_the_store_through_the_link_is_refused_and_nothing_appears(
    armed: dict[str, Path],
) -> None:
    probe = armed["link"] / "promogate_write_probe.tmp"
    with pytest.raises(PermissionError, match="may not write to the store"), open(probe, "w") as f:
        f.write("x")
    assert not (armed["store"] / "promogate_write_probe.tmp").exists(), (
        "the refusal must come before a byte moves, so nothing may be left in the real store"
    )


def test_creating_the_same_file_succeeds_with_the_hook_disarmed(planted: dict[str, Path]) -> None:
    guard.disarm()
    probe = planted["link"] / "promogate_write_probe.tmp"
    with open(probe, "w") as f:
        f.write("x")
    assert (planted["store"] / "promogate_write_probe.tmp").exists()
    probe.unlink()


# --- the reads a rebuild depends on are untouched --------------------------------------------------


def test_reading_through_the_symlink_still_works(armed: dict[str, Path]) -> None:
    assert (armed["link"] / "HG002_chr1.vcf.gz").read_bytes() == b"original bytes\n"
    with open(armed["link"] / "HG002_chr1.vcf.gz", "rb") as f:
        assert f.read() == b"original bytes\n"
    assert sorted(p.name for p in armed["link"].iterdir()) == ["HG002_chr1.vcf.gz"]


def test_the_sha256_of_a_store_file_can_still_be_taken(armed: dict[str, Path]) -> None:
    """`check_inputs` hashes every declared input; a guard that broke that would break every rebuild."""
    assert mf.sha256_of(armed["link"] / "HG002_chr1.vcf.gz")[1] == len(b"original bytes\n")


# --- the narrowing, in both directions -------------------------------------------------------------


def test_a_write_inside_the_worktree_is_allowed(armed: dict[str, Path]) -> None:
    out = armed["wt"] / "data" / "results" / "rebuilt.json"
    out.write_text("{}")
    assert out.read_text() == "{}"


def test_a_write_outside_the_worktree_but_outside_every_store_is_allowed(
    armed: dict[str, Path], tmp_path: Path
) -> None:
    """A rebuild writes temporary files and package caches outside its worktree. Refusing every outside
    write would break the rebuilds instead of protecting the stores."""
    elsewhere = tmp_path / "not-a-store.txt"
    elsewhere.write_text("fine")
    assert elsewhere.read_text() == "fine"


def test_an_absolute_path_straight_to_the_real_store_is_refused_too(armed: dict[str, Path]) -> None:
    """The symlink is one route in, not the only one: the check is on the resolved real path."""
    with (
        pytest.raises(PermissionError, match="may not write to the store"),
        open(armed["store"] / "straight_in.tmp", "w"),
    ):
        pass


# --- the other calls that change a path without opening it ----------------------------------------


def test_removing_a_store_file_is_refused(armed: dict[str, Path]) -> None:
    with pytest.raises(PermissionError, match="may not write to the store"):
        os.remove(armed["link"] / "HG002_chr1.vcf.gz")
    assert (armed["store"] / "HG002_chr1.vcf.gz").exists()


def test_renaming_a_store_file_out_of_the_store_is_refused(armed: dict[str, Path], tmp_path: Path) -> None:
    with pytest.raises(PermissionError, match="may not write to the store"):
        os.rename(armed["link"] / "HG002_chr1.vcf.gz", tmp_path / "stolen")
    assert (armed["store"] / "HG002_chr1.vcf.gz").exists()


def test_making_a_directory_in_the_store_is_refused(armed: dict[str, Path]) -> None:
    with pytest.raises(PermissionError, match="may not write to the store"):
        os.mkdir(armed["link"] / "newdir")
    assert not (armed["store"] / "newdir").exists()


def test_the_watched_events_are_a_named_set_and_not_a_pattern() -> None:
    assert "os.remove" in guard.WATCHED and "os.rename" in guard.WATCHED
    assert "open" not in guard.WATCHED, "open is decided by its mode and flags, not by its name"
    # a read-only call must not be on the list, or the guard would refuse the rebuild's own hashing
    for benign in ("os.listdir", "os.scandir", "os.stat", "os.walk"):
        assert benign not in guard.WATCHED


# --- the classification, asked of the pure function without arming anything -------------------------


def test_offence_is_none_when_no_root_is_armed(planted: dict[str, Path]) -> None:
    """Empty roots is the resting state the hook is left in; it must pass everything."""
    args = (str(planted["link"] / "x"), "w", None)
    assert guard.offence("open", args, ()) is None


def test_offence_names_the_path_and_the_store(planted: dict[str, Path]) -> None:
    roots = (str(planted["store"].resolve()),)
    found = guard.offence("open", (str(planted["link"] / "x"), "w", None), roots)
    assert found is not None
    path, root = found
    assert path.endswith("x") and root == str(planted["store"].resolve())


def test_an_integer_file_descriptor_is_not_mistaken_for_a_path(planted: dict[str, Path]) -> None:
    roots = (str(planted["store"].resolve()),)
    assert guard.offence("open", (3, None, os.O_WRONLY), roots) is None
    assert guard.offence("os.truncate", (3, 0), roots) is None


# --- the child process, which is where the rebuilt command actually runs ----------------------------


def _child(planted: dict[str, Path], tmp_path: Path, body: str) -> subprocess.CompletedProcess[str]:
    guard_dir = tmp_path / "guarddir"
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    env.update(guard.write_sitecustomize(guard_dir, Path("scripts/rebuild_write_guard.py").resolve()))
    env[guard.ENV_ROOTS] = str(planted["store"].resolve())
    env["PYTHONPATH"] = str(guard_dir)
    env["LINK"] = str(planted["link"])
    r = subprocess.run([sys.executable, "-c", body], capture_output=True, text=True, env=env)
    r.guard_dir = guard_dir  # type: ignore[attr-defined]
    return r


CHILD_WRITE = """
import os
with open(os.path.join(os.environ["LINK"], "promogate_write_probe.tmp"), "w") as f:
    f.write("x")
print("WROTE")
"""

CHILD_READ = """
import os
print(open(os.path.join(os.environ["LINK"], "HG002_chr1.vcf.gz"), "rb").read().decode().strip())
"""


def test_the_child_process_arms_the_guard_and_leaves_its_evidence(
    planted: dict[str, Path], tmp_path: Path
) -> None:
    """An audit hook in the parent does not reach the subprocess, and the subprocess is where the
    rebuilt command runs. The parent does not assume the child armed: the child records its pid."""
    r = _child(planted, tmp_path, CHILD_READ)
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == "original bytes"
    assert guard.marker_pids(r.guard_dir) == [int(x) for x in guard.marker_pids(r.guard_dir)]
    assert guard.marker_pids(r.guard_dir), "the child must leave evidence that the hook armed"


def test_a_write_in_the_child_process_is_refused(planted: dict[str, Path], tmp_path: Path) -> None:
    r = _child(planted, tmp_path, CHILD_WRITE)
    assert r.returncode != 0, f"the child wrote to the store: {r.stdout}"
    assert "may not write to the store" in r.stderr
    assert "WROTE" not in r.stdout
    assert not (planted["store"] / "promogate_write_probe.tmp").exists()


def test_the_same_child_write_succeeds_without_the_sitecustomize(
    planted: dict[str, Path], tmp_path: Path
) -> None:
    """The counterfactual for the child half: remove the injection and the write goes through."""
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    env.pop(guard.ENV_GUARD, None)
    env.pop(guard.ENV_ROOTS, None)
    env["LINK"] = str(planted["link"])
    r = subprocess.run([sys.executable, "-c", CHILD_WRITE], capture_output=True, text=True, env=env)
    assert r.returncode == 0 and "WROTE" in r.stdout, r.stderr
    (planted["store"] / "promogate_write_probe.tmp").unlink()


# --- the resting state the rest of the process is left in -------------------------------------------


def test_disarming_leaves_the_hook_installed_and_inert(planted: dict[str, Path], tmp_path: Path) -> None:
    guard.arm([str(planted["store"])])
    guard.disarm()
    assert guard.armed() == ()
    # an audit hook cannot be uninstalled, so this is the only state that can be restored
    (planted["store"] / "after_disarm.tmp").write_text("x")
    (planted["store"] / "after_disarm.tmp").unlink()


# --- the tool's own cleanliness: no verdict from verification code nobody can reproduce -------------
#
# Measured by a peer: the guard above sat uncommitted, 232 insertions against HEAD, and a rebuild run in
# that state is judged by code in no commit. MUST_HOLD cannot catch it, because
# `foreign_uncommitted_code_on_the_counting_path` is read off the REBUILT manifest, composed inside the
# worktree, which knows nothing about the checkout the tool itself ran from.
#
# The cases are planted in a real repository holding a real copy of the tool, so the closure is computed
# by the shared function from real import statements and the dirty state is real git state.

_mr_spec = importlib.util.spec_from_file_location("mr2", Path("scripts/manifest_rebuild.py"))
mr = importlib.util.module_from_spec(_mr_spec)
_mr_spec.loader.exec_module(mr)


def _tool_repo(tmp_path: Path, dirty: str | None) -> Path:
    """A repository holding the tool's own files, committed, with `dirty` then edited and left so."""
    repo = tmp_path / "toolrepo"
    for rel in ("scripts/manifest_rebuild.py", "scripts/rebuild_write_guard.py", "genomeos/manifest.py"):
        (repo / rel).parent.mkdir(parents=True, exist_ok=True)
        (repo / rel).write_text(Path(rel).read_text())
    (repo / "genomeos" / "__init__.py").write_text(Path("genomeos/__init__.py").read_text())
    for args in (
        ("init", "-b", "main"),
        ("config", "user.email", "lane@example.invalid"),
        ("config", "user.name", "planted"),
        ("config", "commit.gpgsign", "false"),
        ("add", "-A"),
        ("commit", "-m", "the tool as committed"),
    ):
        subprocess.run(["git", "-C", str(repo), *args], check=True, capture_output=True, env=dict(os.environ))
    if dirty:
        (repo / dirty).write_text((repo / dirty).read_text() + "\n# an uncommitted edit\n")
    return repo


@pytest.fixture(autouse=True)
def _no_inherited_git_environment(monkeypatch: pytest.MonkeyPatch) -> None:
    for name in ("GIT_INDEX_FILE", "GIT_DIR", "GIT_WORK_TREE", "GIT_OBJECT_DIRECTORY"):
        monkeypatch.delenv(name, raising=False)


def test_a_committed_tool_passes_its_own_cleanliness(tmp_path: Path) -> None:
    t = mr.tool_cleanliness(_tool_repo(tmp_path, None))
    assert t["tool_is_committed"] is True
    assert t["uncommitted_code_on_the_tool_s_own_closure"] == []
    assert t["git_sha"]


def test_an_uncommitted_edit_to_the_tool_itself_is_named(tmp_path: Path) -> None:
    t = mr.tool_cleanliness(_tool_repo(tmp_path, "scripts/manifest_rebuild.py"))
    assert t["tool_is_committed"] is False
    assert t["uncommitted_code_on_the_tool_s_own_closure"] == ["scripts/manifest_rebuild.py"]


def test_an_uncommitted_edit_to_the_write_guard_is_named_although_no_import_points_at_it(
    tmp_path: Path,
) -> None:
    """The recursion the ruling warned about. `counting_path` follows import STATEMENTS and the guard is
    loaded by path, so the computed closure does NOT contain it: it is named in TOOL_ENTRIES instead.
    Drop it from TOOL_ENTRIES and this test fails, which is the only way the choice stays visible."""
    repo = _tool_repo(tmp_path, "scripts/rebuild_write_guard.py")
    t = mr.tool_cleanliness(repo)
    assert t["tool_is_committed"] is False
    assert t["uncommitted_code_on_the_tool_s_own_closure"] == ["scripts/rebuild_write_guard.py"]
    computed = mf.counting_path("scripts/manifest_rebuild.py", repo)
    assert "scripts/rebuild_write_guard.py" not in computed, (
        "if the import walk ever does find the guard, this test's reason for existing has changed and "
        "the naming in TOOL_ENTRIES should be reconsidered rather than left as a duplicate"
    )


def test_an_uncommitted_edit_to_a_module_the_tool_imports_is_named(tmp_path: Path) -> None:
    t = mr.tool_cleanliness(_tool_repo(tmp_path, "genomeos/manifest.py"))
    assert t["uncommitted_code_on_the_tool_s_own_closure"] == ["genomeos/manifest.py"]
    assert t["tool_is_committed"] is False


def test_a_dirty_file_off_the_tool_s_closure_does_not_make_the_tool_dirty(tmp_path: Path) -> None:
    """The counterfactual that keeps the refusal from being "any dirty file anywhere"."""
    repo = _tool_repo(tmp_path, None)
    (repo / "scripts" / "unrelated_lane_script.py").write_text("VALUE = 1\n")
    subprocess.run(
        ["git", "-C", str(repo), "add", "--", "scripts/unrelated_lane_script.py"],
        check=True,
        capture_output=True,
        env=dict(os.environ),
    )
    t = mr.tool_cleanliness(repo)
    assert t["tool_is_committed"] is True, t["uncommitted_code_on_the_tool_s_own_closure"]


def _planted_result(tmp_path: Path) -> Path:
    """A well formed result, so the only thing that can stop a rebuild of it is the tool."""
    result = tmp_path / "planted.json"
    result.write_text(
        json.dumps(
            {
                "result": "planted",
                mf.KEY: {"code": {"git_sha": "0" * 40, "dirty": False, "argv": ["scripts/x.py"]}},
            }
        )
    )
    return result


def test_the_refusal_is_a_refusal_and_not_a_field(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A dirty tool gives NO VERDICT: nothing is built, no command runs, `rebuilt` is False.

    The dirty state is real git state in a real repository holding a real copy of the tool, so nothing
    here is patched and a refusal that stopped working would fail this.
    """
    repo = _tool_repo(tmp_path, "scripts/rebuild_write_guard.py")
    result = _planted_result(tmp_path)
    monkeypatch.chdir(repo)
    r = mr.rebuild(result, tmp_path / "where")
    assert r["rebuilt"] is False
    assert any("no verdict: the verification tool is uncommitted" in u for u in r["unavailable"])
    assert any("scripts/rebuild_write_guard.py" in u for u in r["unavailable"]), "the refusal names it"
    assert "differences" not in r, "a refusal must not carry a comparison"
    assert not (tmp_path / "where").exists(), "nothing may be built when the tool is uncommitted"


def test_with_the_tool_committed_the_rebuild_gets_past_the_tool_check(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The counterfactual. With the same result and a CLEAN tool the refusal does not fire; the rebuild
    goes on to fail at `git worktree add` on the planted result's unknown sha, which is a different
    failure in a different place. So the refusal above is the tool check and not the result."""
    repo = _tool_repo(tmp_path, None)
    result = _planted_result(tmp_path)
    monkeypatch.chdir(repo)
    with pytest.raises(subprocess.CalledProcessError):
        mr.rebuild(result, tmp_path / "where")


# --- read recording, and the two ways a recorder gives a pass it has not earned ---------------------
#
# Measured on a live result: `constrained_unknown_targets` rebuilt at 193 of 193 inputs declared and
# checked, 0 differences after one leaf set aside -- while a tracer on the same run saw a 194th read,
# `data/results/unknown_chr21.json`, pinned by no declared sha256. Two instruments disagreed about one
# run and the flattering one was the clean bill of health.
#
# Both failure cases below really happened in a peer's first two drafts, and each is planted with the
# OLD WAY beside it, failing. A counterfactual that is not run is an argument.

READ_FOUR = """
import io, os, sys
d = sys.argv[1]
open(os.path.join(d, "a.txt")).read()                 # builtins.open
io.open(os.path.join(d, "b.txt")).read()              # io.open
os.close(os.open(os.path.join(d, "c.txt"), os.O_RDONLY))
__import__("pathlib").Path(os.path.join(d, "d.txt")).read_text()
"""


@pytest.fixture
def four_files(tmp_path: Path) -> Path:
    d = tmp_path / "reads"
    d.mkdir()
    for name in ("a.txt", "b.txt", "c.txt", "d.txt"):
        (d / name).write_text(name)
    return d


@pytest.fixture
def recording(tmp_path: Path):
    """The recorder armed into its own directory, and closed again however the test ends."""
    where = tmp_path / "record"
    guard.record_to(where)
    try:
        yield where
    finally:
        guard.record_to(None)


def test_the_recorder_sees_every_way_python_opens_a_file(recording: Path, four_files: Path) -> None:
    """The reason the instrument is a hook: it sees io.open, os.open and pathlib, not only builtins."""
    for name in ("a.txt", "b.txt", "c.txt", "d.txt"):
        (four_files / name).read_text()
    got = guard.reads_recorded(recording)
    assert got["recorded"] is True
    names = {Path(p).name for p in got["paths"]}
    assert {"a.txt", "b.txt", "c.txt", "d.txt"} <= names, names


def test_a_builtins_open_patch_records_one_of_the_four_the_hook_records(
    tmp_path: Path, four_files: Path
) -> None:
    """Case (a), with the old way beside it. A patch of `builtins.open` leaves `io.open` -- and so all of
    pathlib, and pandas -- invisible: 25 reads recorded where there were 3,259. Measured here at the
    same shape rather than quoted: the patch sees 1 of 4, the hook sees 4 of 4."""
    patched = """
import builtins, io, os, sys
seen = []
real = builtins.open
builtins.open = lambda *a, **k: (seen.append(a[0]), real(*a, **k))[1]
d = sys.argv[1]
open(os.path.join(d, "a.txt")).read()
io.open(os.path.join(d, "b.txt")).read()
os.close(os.open(os.path.join(d, "c.txt"), os.O_RDONLY))
__import__("pathlib").Path(os.path.join(d, "d.txt")).read_text()
print(len({os.path.basename(str(s)) for s in seen}))
"""
    r = subprocess.run([sys.executable, "-c", patched, str(four_files)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert int(r.stdout.strip()) == 1, (
        f"the builtins patch saw {r.stdout.strip()} of 4; if it ever sees all four, this test's reason "
        f"for existing has changed"
    )

    guard.record_to(tmp_path / "hookrecord")
    try:
        for name in ("a.txt", "b.txt", "c.txt", "d.txt"):
            (four_files / name).read_text()
        got = guard.reads_recorded(tmp_path / "hookrecord")
    finally:
        guard.record_to(None)
    assert len({Path(p).name for p in got["paths"]} & {"a.txt", "b.txt", "c.txt", "d.txt"}) == 4


def _recording_child(record_dir: Path, body: str, *args: str) -> subprocess.CompletedProcess[str]:
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    env[guard.ENV_RECORD_DIR] = str(record_dir)
    env["GUARD"] = str(Path("scripts/rebuild_write_guard.py").resolve())
    pre = (
        "import importlib.util, os\n"
        "s = importlib.util.spec_from_file_location('g', os.environ['GUARD'])\n"
        "g = importlib.util.module_from_spec(s); s.loader.exec_module(g)\n"
        "g.arm([])\n"
        "g.record_to(os.environ['GENOMEOS_REBUILD_READ_RECORD_DIR'])\n"
    )
    return subprocess.run([sys.executable, "-c", pre + body, *args], capture_output=True, text=True, env=env)


def test_a_process_that_read_nothing_cannot_overwrite_a_reader_s_record(
    tmp_path: Path, four_files: Path
) -> None:
    """Case (b). Under `--workers 2` a pool child that read nothing reached exit last and OVERWROTE the
    single record, and the result came back with 0 paths opened -- "every read is declared", vacuously.

    Here the silent child runs SECOND, which is the losing order for the old design.
    """
    record = tmp_path / "record"
    reader = _recording_child(record, "import sys\nopen(sys.argv[1]).read()\n", str(four_files / "a.txt"))
    assert reader.returncode == 0, reader.stderr
    silent = _recording_child(record, "pass\n")
    assert silent.returncode == 0, silent.stderr

    got = guard.reads_recorded(record)
    assert "a.txt" in {Path(p).name for p in got["paths"]}, (
        "the reader's record was lost after a process that read nothing ran: that is case (b)"
    )
    assert len(got["record_files"]) == 2, got["record_files"]
    assert len(got["processes"]) >= 1


def test_the_old_single_file_design_really_does_lose_it(tmp_path: Path, four_files: Path) -> None:
    """The counterfactual for case (b), run rather than asserted: one shared file written at exit, the
    silent child last, and the reader's path is gone. This is what the per-process append-only record
    is for, and if this test ever passes the way the one above does, the design no longer matters."""
    shared = tmp_path / "one-record.tsv"
    old = (
        "import atexit, os\n"
        "atexit.register(lambda: open(os.environ['R'], 'w').write(os.environ.get('P', '')))\n"
    )
    env = dict(os.environ)
    env["R"] = str(shared)
    env["P"] = str(four_files / "a.txt")
    subprocess.run([sys.executable, "-c", old], capture_output=True, text=True, env=env, check=True)
    assert "a.txt" in shared.read_text()
    env["P"] = ""  # the child that read nothing
    subprocess.run([sys.executable, "-c", old], capture_output=True, text=True, env=env, check=True)
    assert "a.txt" not in shared.read_text(), (
        "the single-file design must lose the reader's record here, or it is not the failure being "
        "guarded against"
    )


def test_an_absent_record_is_not_the_same_claim_as_no_reads(tmp_path: Path) -> None:
    got = guard.reads_recorded(tmp_path / "never-written")
    assert got["recorded"] is False
    assert got["paths"] == [] and got["record_files"] == []


def test_a_write_is_not_recorded_as_a_read(recording: Path, tmp_path: Path) -> None:
    (tmp_path / "written.txt").write_text("x")
    got = guard.reads_recorded(recording)
    assert "written.txt" not in {Path(p).name for p in got["paths"]}


def test_the_same_path_read_twice_is_recorded_once(recording: Path, four_files: Path) -> None:
    for _ in range(5):
        (four_files / "a.txt").read_text()
    lines = [p for p in guard.reads_recorded(recording)["paths"] if Path(p).name == "a.txt"]
    assert len(lines) == 1


def test_the_recorder_states_its_limit_where_it_is_read(tmp_path: Path) -> None:
    """The limit travels with the record, not only with a report: a C library opening a file itself is
    invisible to an audit hook, and this project streams BAMs through pysam."""
    limit = guard.reads_recorded(tmp_path)["limit"]
    assert "pysam" in limit and "htslib" in limit
    assert "Python-level reads" in limit
    assert "counts_not_recorded_because" in guard.reads_recorded(tmp_path)


def test_disarming_stops_the_recording(tmp_path: Path, four_files: Path) -> None:
    where = tmp_path / "record"
    guard.record_to(where)
    guard.disarm()
    (four_files / "a.txt").read_text()
    assert guard.reads_recorded(where)["recorded"] is False
