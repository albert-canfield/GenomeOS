# SPDX-License-Identifier: AGPL-3.0-or-later
"""The entry script is counted: a result cannot certify code that was never committed.

THE MEASURED DEFECT, not a worry. Until 2026-10-02 `genomeos.manifest.code_cleanliness` compared
git's dirty list against the paths the CALLER DECLARED and never looked at the script that was
actually running, so a script outside the repository could write a result asserting that its code was
committed. It did. On `data/results/organised_chr21.json`
(working tree at d4d0449, 19:41; superseded by e1def2e, 20:41,
so these values are NOT reproducible from the tree; the cause is recorded in `2ebe9df` and the
scratchpad script is the live artefact):

    code.argv[0]             /private/tmp/claude-502/.../scratchpad/write.py
    code_cleanliness.git_sha d4d0449bb0d7aaa114b3c8b3d553974e863d86b9
    code_cleanliness.dirty   True
    own_code_is_committed    True        <- false
    argv[0] in counting_path  no

Twenty-four `organised_chr*.json` files were written that way at 19:41-19:42, each with a 39-entry
counting path holding only `genomeos/` modules and no script at all. Nothing in those manifests
describes the code that produced the bytes, and `scripts/manifest_rebuild.py` has nothing to run.
That is worse than a stamp that has gone out of date, which was the morning's other finding: this
claim was affirmatively false.

WHAT IS PLANTED HERE, and the two cases differ in ONE thing. A git repository is planted with a
writer script committed at `scripts/planted_writer.py`; a byte-identical copy of that script is put
in a directory outside the repository. Both runs declare the SAME entry (`scripts/planted_writer.py`,
the committed path), carry the same manifest and write the same payload with `strict=True`; the only
difference is which copy the interpreter is pointed at. The outside copy must be REFUSED and the
inside copy must PASS. A check that refuses everything looks safe and is useless, so the near-miss is
asserted as hard as the refusal, and both run as real subprocesses: `sys.argv[0]` is read from the
process, so a test that called the function in-process would be testing an argument instead.

The odd forms are decided in `genomeos.manifest.ENTRY_FORMS` rather than when a census needs an
answer, and each decision is exercised below: `-c`, stdin, a pytest-driven call, an empty argv, a
path that is not a file, and a `-m` run -- that last one as a real subprocess, because the claim that
`-m` needs no case of its own rests on `runpy` setting `sys.argv[0]` to the module's file, which is a
premise about the interpreter and not about this repository.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from genomeos import manifest as mf

REPO = Path(__file__).resolve().parent.parent

#: The committed location of the planted writer, inside the planted repository. `mf.is_code` counts a
#: path as code only under CODE_ROOTS, so a planted writer at the repository root would be invisible
#: to the revision stamp and both cases below could pass for the wrong reason.
WRITER = "scripts/planted_writer.py"

#: The writer. It declares `WRITER` as its entry in BOTH cases -- that declaration is exactly what the
#: defect trusted -- and lets `code_cleanliness` read the truth from `sys.argv[0]` itself.
WRITER_SOURCE = """\
import json
import sys

from genomeos import manifest as mf
from genomeos import results as rs

ENTRY = "scripts/planted_writer.py"

manifest = {
    "sources": [{"accession": "planted", "version": "1"}],
    "inputs": [mf.input_entry("data/planted_input.txt")],
    "assembly": "n/a: the planted input is one line of text",
    "coordinates": "n/a: the planted input is one line of text",
    "parameters": {"k": 1},
    "exclusions": [],
    "partitions": "n/a: one arm",
    "code_cleanliness": mf.code_cleanliness(ENTRY, (ENTRY,)),
}
# GENOMEOS_PLANTED_OUT sends the write OUTSIDE the registry, which is what the exemption's
# counterfactual needs: the same bytes, the same block, a different destination.
import os

out = os.environ.get("GENOMEOS_PLANTED_OUT")
where = {"results_dir": __import__("pathlib").Path(out)} if out else {}
try:
    p = rs.save_result("planted", {"value": 1}, manifest=manifest, strict=True, **where)
except mf.ManifestError as e:
    print(json.dumps({"refused": True, "why": str(e)}))
    sys.exit(0)
print(json.dumps({"refused": False, "path": str(p)}))
"""


def _git(repo: Path, *args: str) -> None:
    subprocess.run(["git", "-C", str(repo), *args], check=True, capture_output=True, text=True)


def _planted_repo(tmp_path: Path) -> Path:
    """A git repository holding the writer, committed. `-c` keeps this independent of any user config."""
    repo = tmp_path / "planted"
    (repo / "scripts").mkdir(parents=True)
    (repo / WRITER).write_text(WRITER_SOURCE)
    (repo / "data" / "results").mkdir(parents=True)
    (repo / "data" / "planted_input.txt").write_text("1\n")
    _git(repo.parent, "init", "-q", "planted")
    _git(repo, "-c", "user.email=p@p", "-c", "user.name=p", "add", WRITER)
    _git(repo, "-c", "user.email=p@p", "-c", "user.name=p", "commit", "-q", "-m", "planted writer")
    return repo


def _run(repo: Path, script: Path, out_dir: Path | None = None) -> dict:
    """Run `script` with the planted repository as the working directory and read its one JSON line.

    `out_dir` sends the write outside the planted registry, for the exemption's counterfactual."""
    env = {**os.environ, "PYTHONPATH": str(REPO), "GIT_CONFIG_GLOBAL": "/dev/null"}
    if out_dir is not None:
        env["GENOMEOS_PLANTED_OUT"] = str(out_dir)
    done = subprocess.run(
        [sys.executable, str(script)], cwd=repo, env=env, capture_output=True, text=True, timeout=180
    )
    assert done.returncode == 0, f"the planted writer failed: {done.stdout}\n{done.stderr}"
    return json.loads(done.stdout.strip().splitlines()[-1])


# --- the refusal: a script in a temp directory ----------------------------------------------------


def test_a_writer_outside_the_repository_is_refused_and_nothing_enters_the_registry(
    tmp_path: Path,
) -> None:
    repo = _planted_repo(tmp_path)
    outside = tmp_path / "scratchpad"
    outside.mkdir()
    (outside / "write.py").write_text(WRITER_SOURCE)  # byte-identical to the committed copy

    out = _run(repo, outside / "write.py")

    assert out["refused"] is True, "a script outside the repository wrote a result into the registry"
    assert "entry script is not committed" in out["why"]
    assert "outside_repository" in out["why"]
    assert not (repo / "data" / "results" / "planted.json").exists()
    quarantined = json.loads((repo / "data" / "quarantine" / "results" / "planted.json").read_text())
    block = quarantined["result_manifest"]["code_cleanliness"]
    assert block["own_code_is_committed"] is False
    assert block["entry_script"]["form"] == "outside_repository"
    assert block["entry_script"]["is_committed"] is False
    # The declaration the defect trusted is still there, and now it is visibly contradicted.
    assert block["entry_script"]["declared_entry"] == WRITER
    assert str(outside) in block["entry_script"]["argv0"]
    assert block["entry_script"]["argv0_is_on_the_counting_path"] is False
    # The script the lane ran is named as the lane's own uncommitted code, which is what
    # scripts/check_staged.py reads, so the commit of such a result is refused as well.
    assert block["entry_script"]["argv0"] in block["own_uncommitted_code"]


# --- the near-miss: the same script, committed, inside the repository -----------------------------


def test_the_same_writer_committed_inside_the_repository_passes(tmp_path: Path) -> None:
    """The near-miss. A check that refuses every form would pass the test above and be useless."""
    repo = _planted_repo(tmp_path)

    out = _run(repo, repo / WRITER)

    assert out["refused"] is False, f"a committed script inside the repository was refused: {out}"
    written = json.loads((repo / "data" / "results" / "planted.json").read_text())
    block = written["result_manifest"]["code_cleanliness"]
    assert block["own_code_is_committed"] is True
    assert block["own_uncommitted_code"] == []
    assert block["entry_script"]["form"] == "committed_repository_file"
    assert block["entry_script"]["is_committed"] is True
    assert block["entry_script"]["path"] == WRITER
    assert block["entry_script"]["declared_entry_matches_argv0"] is True
    assert written["result_manifest"]["complete"] is True


def test_the_two_cases_differ_only_in_where_the_script_lives(tmp_path: Path) -> None:
    """The bytes are the same, so the verdict is about the location and not about the code."""
    repo = _planted_repo(tmp_path)
    outside = tmp_path / "scratchpad"
    outside.mkdir()
    (outside / "write.py").write_text(WRITER_SOURCE)
    assert (outside / "write.py").read_bytes() == (repo / WRITER).read_bytes()


def test_a_committed_script_edited_in_place_is_refused(tmp_path: Path) -> None:
    """The third condition. Tracked and modified is not committed, which is how a lane actually gets
    there: the file is in the commit, so a path-only check would call it clean."""
    repo = _planted_repo(tmp_path)
    (repo / WRITER).write_text(WRITER_SOURCE + "\n# edited after the commit\n")

    out = _run(repo, repo / WRITER)

    assert out["refused"] is True
    assert "uncommitted_repository_file" in out["why"]


def test_an_untracked_script_inside_the_repository_is_refused(tmp_path: Path) -> None:
    repo = _planted_repo(tmp_path)
    (repo / "scripts" / "untracked_writer.py").write_text(WRITER_SOURCE)

    out = _run(repo, repo / "scripts" / "untracked_writer.py")

    assert out["refused"] is True
    assert "uncommitted_repository_file" in out["why"]


# --- the odd forms, each decided in advance -------------------------------------------------------


@pytest.mark.parametrize(
    ("argv0", "form"),
    [
        ("-c", "inline_source"),
        ("-", "stdin"),
        ("", "no_entry_script"),
        ("/usr/local/bin/pytest", "test_runner"),
        (".venv/bin/pytest", "test_runner"),
        (".venv/lib/python3.12/site-packages/pytest/__main__.py", "test_runner"),
        ("scripts/this_file_does_not_exist_anywhere.py", "not_a_file"),
    ],
)
def test_each_odd_form_is_classified_as_registered_and_none_of_them_is_clean(
    tmp_path: Path, argv0: str, form: str
) -> None:
    repo = _planted_repo(tmp_path)
    block = mf.entry_script(repo, argv0)
    assert block["form"] == form
    assert block["is_committed"] is False
    assert block["why"] == mf.ENTRY_FORMS[form]


def test_a_console_script_shim_is_refused_as_an_uncommitted_repository_file(tmp_path: Path) -> None:
    """The decision for `uv run genomeos ...`, planted rather than reasoned about. The shim the
    installer generates is git-ignored and in no commit, so it refuses; the committed module it
    dispatches to is reached by `python -m`, which the test above covers."""
    repo = _planted_repo(tmp_path)
    shim = repo / ".venv" / "bin" / "genomeos"
    shim.parent.mkdir(parents=True)
    shim.write_text('#!/bin/sh\nexec python -m genomeos "$@"\n')
    block = mf.entry_script(repo, ".venv/bin/genomeos")
    assert block["form"] == "uncommitted_repository_file"
    assert block["is_committed"] is False


def test_the_registered_forms_and_the_classifier_name_the_same_set(tmp_path: Path) -> None:
    """ENTRY_FORMS is the decision, so a form the classifier can return and the table does not
    describe would be a decision nobody made. `entry_script` reads `why` out of the table, so such a
    form raises rather than returning an undescribed verdict -- asserted here as set equality."""
    import ast
    import inspect

    returned = {
        node.value.elts[0].value
        for node in ast.walk(ast.parse(textwrap.dedent(inspect.getsource(mf._classify_entry))))
        if isinstance(node, ast.Return) and isinstance(node.value, ast.Tuple)
        for element in [node.value.elts[0]]
        if isinstance(element, ast.Constant)
    }
    named = {mf.ENTRY_FORM_CLEAN, "uncommitted_repository_file"}
    assert returned | named == set(mf.ENTRY_FORMS), (
        f"classifier returns {sorted(returned | named)}, ENTRY_FORMS describes {sorted(mf.ENTRY_FORMS)}"
    )


def test_exactly_one_form_is_clean() -> None:
    """There is no exemption to argue about: every form but one refuses."""
    assert mf.ENTRY_FORM_CLEAN in mf.ENTRY_FORMS
    assert sum(mf.entry_script(REPO, f)["is_committed"] for f in ("-c", "-", "", "pytest")) == 0


def test_a_module_run_with_dash_m_is_classified_by_the_module_s_own_file(tmp_path: Path) -> None:
    """The premise behind having no `-m` case: `runpy` sets `sys.argv[0]` to the module's file. Read
    from a real subprocess, so if that ever stopped being true this fails instead of the census
    silently reclassifying every `-m` run as something else."""
    repo = _planted_repo(tmp_path)
    (repo / "scripts" / "__init__.py").write_text("")
    (repo / "scripts" / "planted_form.py").write_text(
        "import json\nfrom genomeos import manifest as mf\nprint(json.dumps(mf.entry_script(None, None)))\n"
    )
    _git(
        repo,
        "-c",
        "user.email=p@p",
        "-c",
        "user.name=p",
        "add",
        "scripts/__init__.py",
        "scripts/planted_form.py",
    )
    _git(repo, "-c", "user.email=p@p", "-c", "user.name=p", "commit", "-q", "-m", "form probe")
    env = {**os.environ, "PYTHONPATH": f"{REPO}{os.pathsep}{repo}", "GIT_CONFIG_GLOBAL": "/dev/null"}
    done = subprocess.run(
        [sys.executable, "-m", "scripts.planted_form"],
        cwd=repo,
        env=env,
        capture_output=True,
        text=True,
        timeout=180,
    )
    assert done.returncode == 0, done.stderr
    block = json.loads(done.stdout.strip().splitlines()[-1])
    assert block["argv0"].endswith("scripts/planted_form.py"), block["argv0"]
    assert block["form"] == "committed_repository_file"
    assert block["path"] == "scripts/planted_form.py"


# --- the counterfactual ---------------------------------------------------------------------------

#: The line of `genomeos/results.py` that refuses. Stripped from a copy of the function's own source
#: below rather than paraphrased, so the counterfactual cannot drift from the implementation: if this
#: text stops matching, the test fails instead of passing vacuously.
REFUSAL = "unentered = mf.entry_script_problems(stamped)"


def test_without_the_refusal_the_outside_script_would_be_published(tmp_path: Path) -> None:
    """The counterfactual for the refusal itself: with that one line gone, the planted outside script
    writes `planted.json` into the registry with `own_code_is_committed` false and nothing stops it.
    That is the state the 24 committed results are in."""
    import inspect

    from genomeos import results as rs

    source = inspect.getsource(rs.save_result)
    assert source.count(REFUSAL) == 1, (
        "the line this counterfactual strips is no longer in genomeos/results.py as written here, so "
        "the counterfactual would pass without testing anything: update REFUSAL"
    )
    repo = _planted_repo(tmp_path)
    outside = tmp_path / "scratchpad"
    outside.mkdir()
    patched = tmp_path / "patched_results.py"
    whole = Path(rs.__file__).read_text()
    patched.write_text(whole.replace(REFUSAL, "unentered = []  # counterfactual: refusal removed"))
    (outside / "write.py").write_text(
        WRITER_SOURCE.replace(
            "from genomeos import results as rs",
            "import importlib.util as _u\n"
            f"_s = _u.spec_from_file_location('patched_results', {str(patched)!r})\n"
            "rs = _u.module_from_spec(_s)\n_s.loader.exec_module(rs)",
        )
    )

    out = _run(repo, outside / "write.py")

    assert out["refused"] is False, "the counterfactual did not reproduce the defect"
    written = json.loads((repo / "data" / "results" / "planted.json").read_text())
    block = written["result_manifest"]["code_cleanliness"]
    assert block["own_code_is_committed"] is False
    assert block["entry_script"]["form"] == "outside_repository"


# --- the census reads and never writes ------------------------------------------------------------


def _census_module():
    import importlib.util

    spec = importlib.util.spec_from_file_location("esc", "scripts/entry_script_census.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize(
    ("argv0", "verdict"),
    [
        (None, "unrecorded"),
        ("-c", "inline_source"),
        ("-", "stdin"),
        ("scripts/organise.py", "inside_repository"),
        ("/private/tmp/claude-502/x/scratchpad/write.py", "outside_repository"),
        (f"{REPO}/scripts/organise.py", "inside_recorded_absolute"),
        ("/usr/local/bin/pytest", "test_runner"),
    ],
)
def test_the_census_reads_the_recorded_argv0_as_registered(argv0: str | None, verdict: str) -> None:
    """The census's verdict comes from the recorded string, not from the present tree, because
    `manifest._argv` relativises argv[0] exactly when the script was under the root. An absolute path
    that nonetheless lies under this checkout is counted apart rather than read as outside."""
    census = _census_module()
    assert census.verdict_for(argv0) == verdict
    assert verdict in census.VERDICTS


def test_the_census_does_not_write_any_result(tmp_path: Path) -> None:
    """It reports a finding; it does not repair one. Several of these results' sha256 are
    registrations in data/results/manifest_headlines.json and a rewrite would break a pin."""
    census = _census_module()
    before = {p: p.stat().st_mtime_ns for p in Path("data/results").glob("*.json")}
    out = census.census(sorted(census.readme_cited()[0]))
    assert out["examined"] > 0
    assert {p: p.stat().st_mtime_ns for p in Path("data/results").glob("*.json")} == before
    assert "writes none of them" in out["nothing_was_written"]


def test_the_census_population_is_the_readme_definition_already_in_use() -> None:
    """Not a new definition: the names `tests/test_headline_registry_matches_readme.py` uses."""
    census = _census_module()
    record = json.loads(Path("data/results/manifest_headlines.json").read_text())
    registered = {e["result"] for e in record["rebuilt"] + record.get("pending", [])}
    names, _ = census.readme_cited()
    assert registered <= set(names)


def test_the_census_separates_a_flagged_result_s_own_bytes_from_its_claim(tmp_path: Path) -> None:
    """A working-tree file about to be discarded is not a published result that certified itself
    falsely. `bytes_state` is the column that keeps a reader from counting the second where only the
    first happened: the 24 files that prompted this lane are modifications whose committed versions
    are older, legacy-shaped results carrying no `result_manifest` at all."""
    census = _census_module()
    untracked = Path("data/results/this_result_is_not_tracked_anywhere.json")
    assert census.bytes_state(untracked) == "untracked"
    assert census.bytes_state(Path("README.md")) in ("committed", "modified")
    rows = census.census(sorted(census.readme_cited()[0]))["rows"]
    assert all(r["bytes_state"] in ("committed", "modified", "untracked", "unknown") for r in rows)


def test_a_form_git_cannot_answer_for_is_unknown_and_refuses(tmp_path: Path) -> None:
    """The last row of the table, planted. Where git cannot say whether the file is committed -- here
    a directory that is not a repository at all -- the verdict is `unknown` and it refuses. A check
    that fell back to clean when it could not tell would be worse than no check."""
    outside_any_repo = tmp_path / "not_a_repository"
    (outside_any_repo / "scripts").mkdir(parents=True)
    (outside_any_repo / WRITER).write_text(WRITER_SOURCE)
    block = mf.entry_script(outside_any_repo, WRITER)
    assert block["form"] == "unknown"
    assert block["is_committed"] is False
    assert mf.entry_script_problems({"code_cleanliness": {"entry_script": block}}) != []


# --- the enforcement point, and the counterfactual for its one exemption --------------------------
#
# The refusal fires on a write into the REGISTRY, legacy name or not, and not on `strict` alone. That
# is stronger than the brief's wording where it matters -- `organised_chr*`, the 24 names this check
# exists for, are LEGACY names, for which `enforce` is `bool(strict)` -- and weaker only outside the
# registry, where no published result ever lives. What it buys, measured: with `enforce` as the
# condition, nine positive controls across tests/test_traced_inputs.py, tests/test_organise_inputs.py
# and tests/test_attribution.py went red, every one of them a real writer driven into a temporary
# directory under pytest. The exemption is not argued harmless; the pair below is its counterfactual.


def test_outside_the_registry_the_same_write_is_recorded_and_not_refused(tmp_path: Path) -> None:
    """The exemption. The same script, the same block, a scratch destination: it writes.

    And it is NOT certified -- that is the half that makes the exemption safe to have, so it is
    asserted rather than assumed. `own_code_is_committed` is false and the form is recorded, which is
    exactly what `scripts/check_staged.py` reads to refuse the commit of such a file.
    """
    repo = _planted_repo(tmp_path)
    outside = tmp_path / "scratchpad"
    outside.mkdir()
    (outside / "write.py").write_text(WRITER_SOURCE)
    scratch = tmp_path / "scratch_results"
    scratch.mkdir()

    out = _run(repo, outside / "write.py", out_dir=scratch)

    assert out["refused"] is False, "a scratch write must not be refused for the entry script"
    block = json.loads((scratch / "planted.json").read_text())["result_manifest"]["code_cleanliness"]
    assert block["own_code_is_committed"] is False, "a scratch write must still not be certified"
    assert block["entry_script"]["form"] == "outside_repository"


def test_into_the_registry_that_same_write_refuses(tmp_path: Path) -> None:
    """The other side of the pair, so the exemption cannot be the whole rule: destination is the only
    difference between this and the test above."""
    repo = _planted_repo(tmp_path)
    outside = tmp_path / "scratchpad"
    outside.mkdir()
    (outside / "write.py").write_text(WRITER_SOURCE)

    out = _run(repo, outside / "write.py")  # no out_dir: the planted repository's data/results

    assert out["refused"] is True
    assert "entry script must be committed for any write into the registry" in out["why"]


def test_a_legacy_registry_name_is_refused_too(tmp_path: Path) -> None:
    """The case the brief's wording would have missed, and the case that actually happened: a LEGACY
    registry name, for which `enforce` is `bool(strict)` and the rest of the contract only warns.
    `organised_chr21` is such a name. Planted by putting `planted` on the allowlist."""
    repo = _planted_repo(tmp_path)
    (repo / "data" / "results_legacy.txt").write_text("# a header\nplanted\ttracked\n")
    outside = tmp_path / "scratchpad"
    outside.mkdir()
    (outside / "write.py").write_text(WRITER_SOURCE.replace("strict=True", "strict=False"))

    out = _run(repo, outside / "write.py")

    assert out["refused"] is True, "a legacy name written from outside the repository still refuses"
    assert "entry script" in out["why"]
    assert not (repo / "data" / "results" / "planted.json").exists()


# --- the gate the refusal hangs on must not fail open ---------------------------------------------
#
# `genomeos/results.py: _is_registry` returned False on an OSError from `Path.resolve`. That flag
# gates FOUR decisions in `save_result` and they all go off together on that one exception: `legacy`,
# `enforce` (where a registry write stops being enforced at all and becomes a function of whatever the
# caller passed as `strict`), the requirement to carry a cleanliness block, and the entry-script
# refusal added here. Pre-existing, and harmless while the flag only chose between severities of
# warning; load-bearing the moment a refusal hangs on it. An unclassifiable path now counts AS the
# registry. Both halves are planted, because a fix that refused everything would have broken every
# test writer instead -- which is the failure mode the non-registry exemption exists to avoid.

#: The exact body of `_is_registry`'s guarded return, taken from genomeos/results.py as written. The
#: patched copies below replace this whole block rather than a line that two handlers share, so a
#: future edit makes the assertion fail instead of letting the patch land in the wrong handler.
_IS_REGISTRY_BODY = """    try:
        return results_dir.resolve() == RESULTS_DIR.resolve()
    except OSError:
        return True
"""


def _results_module(handler_returns: str):
    """A copy of `genomeos/results.py` whose `_is_registry` raises OSError and whose handler returns
    `handler_returns`. Patched from the file's own source, so if that body stops being spelled this way
    the test fails instead of passing vacuously."""
    import importlib.util

    from genomeos import results as rs

    whole = Path(rs.__file__).read_text()
    assert whole.count(_IS_REGISTRY_BODY) == 1, (
        "the guarded return this test patches is no longer in genomeos/results.py as written here: "
        "update _IS_REGISTRY_BODY"
    )
    body = whole.replace(
        _IS_REGISTRY_BODY,
        "    try:\n"
        '        raise OSError(62, "ELOOP planted by tests/test_entry_script_counted.py")\n'
        "    except OSError:\n"
        f"        return {handler_returns}\n",
        1,
    )
    out = Path(_scratch()) / f"results_{handler_returns}.py"
    out.write_text(body)
    spec = importlib.util.spec_from_file_location(f"patched_results_{handler_returns}", out)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_SCRATCH: list[Path] = []


def _scratch() -> Path:
    if not _SCRATCH:
        import tempfile

        _SCRATCH.append(Path(tempfile.mkdtemp(prefix="entryscript-")))
    return _SCRATCH[0]


def _write_with(module, repo: Path, script: Path, out_dir: Path | None = None) -> dict:
    """Run the planted outside writer against `module`'s `save_result` instead of the installed one."""
    outside = script.parent
    patched = outside / "write_patched.py"
    patched.write_text(
        WRITER_SOURCE.replace(
            "from genomeos import results as rs",
            "import importlib.util as _u\n"
            f"_s = _u.spec_from_file_location('patched_results', {str(module.__file__)!r})\n"
            "rs = _u.module_from_spec(_s)\n_s.loader.exec_module(rs)",
        )
    )
    return _run(repo, patched, out_dir=out_dir)


def test_a_results_dir_that_cannot_be_resolved_counts_as_the_registry_and_refuses(
    tmp_path: Path,
) -> None:
    """The plant. `resolve` raises, so nothing can say whether this is the registry, and the write is
    REFUSED with the path named -- not recorded and let through."""
    repo = _planted_repo(tmp_path)
    outside = tmp_path / "scratchpad"
    outside.mkdir()
    (outside / "write.py").write_text(WRITER_SOURCE)
    module = _results_module("True")

    out = _write_with(module, repo, outside / "write.py")

    assert out["refused"] is True, "an unresolvable results_dir must be treated as the registry"
    assert "entry script must be committed for any write into the registry" in out["why"]
    assert "data/results" in out["why"], "the refusal must name the path it could not classify"


def test_that_same_test_fails_against_heads_body_which_returned_false(tmp_path: Path) -> None:
    """The fix is a one-line change, so the test must be shown to FAIL on the old body rather than
    reasoned about: a test that passes both ways has tested nothing. HEAD's handler returned False, so
    the same unresolvable write was NOT the registry and went through."""
    repo = _planted_repo(tmp_path)
    outside = tmp_path / "scratchpad"
    outside.mkdir()
    (outside / "write.py").write_text(WRITER_SOURCE)
    module = _results_module("False")

    out = _write_with(module, repo, outside / "write.py")

    assert out["refused"] is False, (
        "with the old handler the unresolvable write must get through, or the test above would pass "
        "whatever the handler returns"
    )


def test_the_near_miss_a_scratch_directory_that_resolves_normally_still_writes(tmp_path: Path) -> None:
    """The other half: the fix must not refuse everything. A scratch directory resolves fine, so it is
    not the registry, and the write happens -- with the false verdict recorded, not hidden."""
    repo = _planted_repo(tmp_path)
    outside = tmp_path / "scratchpad"
    outside.mkdir()
    (outside / "write.py").write_text(WRITER_SOURCE)
    scratch = tmp_path / "scratch_results"
    scratch.mkdir()

    out = _run(repo, outside / "write.py", out_dir=scratch)

    assert out["refused"] is False
    block = json.loads((scratch / "planted.json").read_text())["result_manifest"]["code_cleanliness"]
    assert block["own_code_is_committed"] is False
    assert block["entry_script"]["form"] == "outside_repository"
