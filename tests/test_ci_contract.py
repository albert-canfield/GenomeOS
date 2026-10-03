# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The contract between this checkout and the CI run, so the four defects of 2026-10-02 cannot return.

CI went red on dev on 2026-10-02 at 12:08 UTC and stayed red. Four separate defects, and what the
four have in common is that NONE of them was a wrong result: each was a way for the verification to
say less than its reader thought. They are guarded here, in one file, because they are one subject.

1. THE SHALLOW CHECKOUT. `actions/checkout` defaults to `fetch-depth: 1`, and none of ci.yml's three
   checkouts set it, so no job had any history. A growing family of tests reads the record BY
   REVISION -- `git show <sha>:<path>`, `<sha>^`, `git log -L` -- and every one of them fails with
   exit 128 on a depth-1 clone. In the 12:08 run five appeared as ERRORS at SETUP, all five from the
   one command `git show a39073d:data/results/crispri_published.json`. The exposure GREW the night
   before the run: tests/test_respmap_v2.py pins `ad07bf5^` through `carriers_in_tree`, the
   registration-lineage guards do the same, and scripts/entry_script_census.py walks `git log -L`.
   The fix committed at 4f44dbf uses exactly the pattern that was reddening CI, and had not been run
   by CI when it landed. Guarded two ways below: statically, that every checkout step names
   `fetch-depth: 0`, and dynamically, that this checkout is not shallow.

2. SHELLCHECK WAS A DIFFERENT PROGRAM IN THE TWO PLACES. scripts/check.sh fetches a pinned
   shellcheck-py wheel; ci.yml ran the ubuntu-24.04 runner image's apt shellcheck, 0.9.0-1. Neither
   set a `--severity`, so this was never a severity difference: both default to `style`, the lowest,
   and both reported everything they had. 0.9.0 simply HAS findings 0.11.0 does not -- on the same
   bytes it emitted SC2015 and SC2317 and exited 1 where 0.11.0 emitted nothing and exited 0. The
   consequence was worse than a red lint: `test` died at that step and never ran pytest at all, so
   the run's only finding was the lint and the whole suite went unrun in that job. Guarded below by
   requiring the two files to name the SAME pin, because an unpinned linter in CI can redden the gate
   with no commit at all.

3. uv.lock DRIFTED FROM pyproject.toml. c9e330b added `anndata>=0.13` to the dev group and committed
   pyproject.toml without the matching lock, so `uv sync --frozen` in `test-bare` installed from a
   lock that did not know about it. Guarded below by `uv lock --check`, which is the same question
   CI's `--frozen` asks, asked here where it is cheap -- 3 ms, offline, reading the lock only.

4. THE SKIP SUMMARY WAS FOLDED, so CI's log could not say which tests skipped either. Found by
   lane-skipids and credited to it. `-rs` prints one line per (file, line, reason) --
   `SKIPPED [3] tests/test_x.py:12: reason` -- a count and a location and never an id, so parametrised
   skips collapse and no id appears anywhere. `-rs` is in ci.yml precisely so that what CI did not run
   is on the record; folded, it is half on the record. Guarded below statically, and PLANTED: a
   parametrised skip is run under both spellings in a subprocess, and its id must be absent from one
   and present from the other.

WHAT THIS FILE IS NOT. It holds no allowlist and exempts nothing. Every assertion is about the bytes
of .github/workflows/ci.yml, scripts/check.sh and uv.lock, or about this checkout's own git state.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
CI_YML = ROOT / ".github" / "workflows" / "ci.yml"
CHECK_SH = ROOT / "scripts" / "check.sh"

#: The one shellcheck the project runs. Named here so a reader sees the version in one place; the
#: tests below do not trust this constant, they require ci.yml and check.sh to AGREE with each other.
SHELLCHECK_PIN_PATTERN = re.compile(r"shellcheck-py==(\d+\.\d+\.\d+(?:\.\d+)?)")

#: What the runner image shipped, and the reason the pin exists. From actions/runner-images,
#: images/ubuntu/Ubuntu2404-Readme.md, read 2026-10-03: `| shellcheck | 0.9.0-1 |`.
RUNNER_IMAGE_APT_SHELLCHECK = "0.9.0-1"


# ---------------------------------------------------------------------------
# 1. The shallow checkout
# ---------------------------------------------------------------------------


def test_this_checkout_is_not_shallow() -> None:
    """THE GUARD for defect 1: a shallow clone fails here instead of five tests failing obscurely.

    `git rev-parse --is-shallow-repository` answers `true` or `false` and is the only reliable
    question: a depth-1 clone can be partly unshallowed, and `.git/shallow` existing is not the same
    as the repository still being shallow. On a depth-1 clone this says `true` and this test fails,
    which is the point -- the five SETUP errors of the 12:08 run named a missing object and said
    nothing about why it was missing, and a reader had to infer the checkout depth from the workflow.

    It is a SKIP, with its reason, only where there is no git work tree at all: an unpacked sdist has
    no history to be shallow, and nothing that reads a revision can run there either. CI always has
    one, so the guard is live in the place it exists for.
    """
    git = shutil.which("git")
    if git is None:
        pytest.skip(
            "git is not on PATH, so the depth of this checkout cannot be read; what goes unproven "
            "is only this checkout's own depth -- the static guard on ci.yml's fetch-depth still runs"
        )
    inside = subprocess.run(
        [git, "-C", str(ROOT), "rev-parse", "--is-inside-work-tree"],
        capture_output=True,
        text=True,
    )
    if inside.returncode != 0 or inside.stdout.strip() != "true":
        pytest.skip(
            "this tree is not a git work tree (an unpacked sdist has no history), so there is no "
            "depth to check and no test that reads a revision can run here either"
        )
    shallow = subprocess.run(
        [git, "-C", str(ROOT), "rev-parse", "--is-shallow-repository"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert shallow.stdout.strip() == "false", (
        "this checkout is SHALLOW. Every test that pins a claim to a revision -- `git show "
        "<sha>:<path>` in tests/test_crispri_published.py, `ad07bf5^` in tests/test_respmap_v2.py, "
        "the registration-lineage guards, `git log -L` in scripts/entry_script_census.py -- fails "
        "with exit 128 here, and in CI that reads as five mysterious SETUP errors. In a workflow the "
        "cure is `fetch-depth: 0` on the actions/checkout step; locally, `git fetch --unshallow`."
    )


def test_every_ci_checkout_fetches_the_whole_history() -> None:
    """THE STATIC HALF of defect 1's guard, and the half that can never be skipped.

    actions/checkout's default is depth 1, so a checkout step that says nothing is a shallow one.
    Counting the steps matters as much as checking them: on 2026-10-02 there were three and all three
    were silent, and a guard that checked only the first would have passed over two.
    """
    text = CI_YML.read_text()
    steps = text.split("uses: actions/checkout@")
    occurrences = len(steps) - 1
    assert occurrences >= 3, (
        f"ci.yml names actions/checkout {occurrences} times; the three pytest-running jobs (test, "
        "test-bare, live-readers) each need one, so a smaller number means a job lost its checkout"
    )
    for n, following in enumerate(steps[1:], start=1):
        window = following[:400]
        assert "fetch-depth: 0" in window, (
            f"checkout step {n} of {occurrences} in ci.yml does not set `fetch-depth: 0` within its "
            "own step. actions/checkout defaults to depth 1 and every test that reads a revision "
            "fails with exit 128 on a depth-1 clone; silence here is the defect, not a default."
        )


# ---------------------------------------------------------------------------
# 2. One shellcheck, not two
# ---------------------------------------------------------------------------


def test_ci_and_check_sh_run_the_same_pinned_shellcheck() -> None:
    """THE GUARD for defect 2, and it compares the two files rather than trusting either.

    The defect was not a flag: it was that `scripts/check.sh` ran shellcheck 0.11.0 from a pinned
    wheel and ci.yml ran the runner image's apt shellcheck 0.9.0-1, and the two programs disagree
    about this tree. So the requirement is that both name a pin and that it is the SAME pin. An
    unpinned linter in CI is the thing being forbidden: a runner-image bump would otherwise redden
    the gate with no commit at all, which is the most expensive kind of red this project gets.
    """
    ci_pins = set(SHELLCHECK_PIN_PATTERN.findall(CI_YML.read_text()))
    check_pins = set(SHELLCHECK_PIN_PATTERN.findall(CHECK_SH.read_text()))
    assert len(check_pins) == 1, (
        f"scripts/check.sh names {sorted(check_pins)} shellcheck pins; exactly one is the contract"
    )
    assert len(ci_pins) == 1, (
        f".github/workflows/ci.yml names {sorted(ci_pins)} shellcheck pins. Exactly one is the "
        "contract, and none at all means CI is back to the runner image's apt shellcheck "
        f"({RUNNER_IMAGE_APT_SHELLCHECK}), which is a different program from the one check.sh runs."
    )
    assert ci_pins == check_pins, (
        f"ci.yml pins shellcheck {sorted(ci_pins)} and scripts/check.sh pins {sorted(check_pins)}. "
        "Two versions are two checks, and the one that is not run locally is the one that reddens "
        "CI. They have to be the same string."
    )


def test_no_shell_script_carries_a_bare_shellcheck_disable() -> None:
    """A disable with no reason beside it is the thing the fix for defect 2 was told not to do.

    SC2317 in scripts/pre-push.sh is a false positive -- shellcheck follows no `trap`, so it calls
    the commands inside a trap-invoked function unreachable -- and the resolution was to name it in a
    directive whose stated reason is written beside it. This test is what keeps that an annotation
    with a reason rather than a silencer.

    The reason may sit in any of three places, because this tree already uses two of them: on the
    directive's own line (scripts/pre-push.sh, before this lane), on the comment lines immediately
    ABOVE it (scripts/suite_lock.sh:41-42 explains SC2034 that way), or on the comment lines
    immediately below (how the widened SC2329,SC2317 directive reads now). What is refused is a
    directive with prose in none of the three.
    """
    offenders = []
    for script in sorted((ROOT / "scripts").glob("*.sh")):
        lines = script.read_text().splitlines()
        for i, line in enumerate(lines):
            match = re.search(r"#\s*shellcheck\s+disable=([A-Za-z0-9,=]+)(.*)$", line)
            if not match:
                continue
            reasons = [match.group(2)]
            for neighbour in (i - 1, i + 1):
                if 0 <= neighbour < len(lines):
                    text = lines[neighbour].strip()
                    # A neighbouring DIRECTIVE is not a reason; a neighbouring sentence is, even
                    # when the sentence is about shellcheck, which both of this tree's are.
                    is_directive = re.match(r"#\s*shellcheck\s+\w+=", text)
                    if text.startswith("#") and not is_directive:
                        reasons.append(text)
            if not any(len(r.strip().lstrip("#:").strip()) >= 10 for r in reasons):
                offenders.append(f"{script.name}:{i + 1}: {line.strip()}")
    assert offenders == [], (
        "a bare `# shellcheck disable=` says what to ignore and never why, so the next reader cannot "
        "tell a false positive from a silenced defect. Put the reason beside it:\n" + "\n".join(offenders)
    )


# ---------------------------------------------------------------------------
# 3. The lock and pyproject.toml say the same thing
# ---------------------------------------------------------------------------


def test_uv_lock_is_in_sync_with_pyproject() -> None:
    """THE GUARD for defect 3, asked here because CI asks it too late to be cheap.

    `test-bare` installs with `uv sync --frozen`, which is this question asked at the point where the
    answer costs a red run. `uv lock --check` is the same question, offline and in milliseconds: it
    reads uv.lock and pyproject.toml and nothing else, and it is what c9e330b needed. That commit
    added `anndata>=0.13` to the dev group and committed pyproject.toml alone; `uv lock --check`
    against that pair reports `the lockfile at uv.lock needs to be updated` and exits non-zero --
    measured against the committed pair on 2026-10-03, which is how this guard was shown to work
    rather than assumed to.

    A skip only where `uv` is absent, naming what goes unproven. Both CI jobs invoke pytest through
    uv, so it is present where this matters.
    """
    uv = shutil.which("uv")
    if uv is None:
        pytest.skip(
            "uv is not on PATH, so pyproject.toml and uv.lock cannot be compared here; what goes "
            "unproven is the lock's freshness, which CI's `uv sync --frozen` then asks at full price"
        )
    done = subprocess.run(
        [uv, "lock", "--check", "--offline"],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert done.returncode == 0, (
        "uv.lock and pyproject.toml disagree. `uv sync --frozen` -- which is how ci.yml's `test-bare` "
        "job installs -- fails on exactly this, and the usual cause is a dependency added to "
        "pyproject.toml and committed without the lock beside it (c9e330b, `anndata>=0.13`). Run "
        "`uv lock` and COMMIT uv.lock in the same commit as the pyproject change.\n"
        f"stdout: {done.stdout.strip()}\nstderr: {done.stderr.strip()}"
    )


# ---------------------------------------------------------------------------
# 4. The skip summary names its tests
# ---------------------------------------------------------------------------


def test_every_ci_pytest_invocation_reports_skips_unfolded() -> None:
    """THE STATIC HALF of defect 4's guard: `-rs` without `--no-fold-skipped` is half a report."""
    text = CI_YML.read_text()
    invocations = [
        line.strip()
        for line in text.splitlines()
        if re.search(r"\buv run\b.*\bpytest\b", line) and "--help" not in line
    ]
    assert len(invocations) >= 3, (
        f"ci.yml holds {len(invocations)} pytest invocations; the three jobs each run one, so a "
        f"smaller number means one stopped being matched or stopped existing: {invocations}"
    )
    for invocation in invocations:
        assert "-rs" in invocation, f"no skip report at all in: {invocation}"
        assert "--no-fold-skipped" in invocation, (
            "`-rs` alone prints pytest's FOLDED summary -- one line per (file, line, reason), a "
            "count and a location and never an id -- so parametrised skips collapse into one line "
            f"and CI's log cannot say which tests it did not run: {invocation}"
        )


@pytest.mark.parametrize("case", [1, 2, 3])
def test_planted_skip_folded(case: int) -> None:
    """Three skips at one (file, line, reason). Folded they are one line; unfolded they are three ids.

    This is not a check of anything. It is the PLANT the next test reads, left in the suite on
    purpose so that the folding behaviour is demonstrated against real pytest output rather than
    asserted from the changelog. It also means this project's own CI log carries three ids that a
    folded run would have hidden.
    """
    pytest.skip("planted: three parametrised skips share one file, line and reason")


def test_the_unfolded_summary_names_the_planted_skip_by_id() -> None:
    """THE PLANTED HALF of defect 4's guard: the id appears under one spelling and not the other.

    A copy of the plant above is run in a subprocess twice, with `-rs` and with
    `-rs --no-fold-skipped`, and the test is that `test_planted_skip_folded[1]` is ABSENT from the
    first and PRESENT in the second. Asserting only the second would pass even if pytest never
    folded anything, and the claim being made is about the difference.

    The copy is written into a tmp_path outside the repository so that pytest's rootdir is that
    directory and this project's conftest.py is not collected; `-p no:cacheprovider` keeps it from
    writing anything.
    """
    import tempfile

    source = (
        "import pytest\n"
        "\n"
        "\n"
        '@pytest.mark.parametrize("case", [1, 2, 3])\n'
        "def test_planted_skip_folded(case):\n"
        '    pytest.skip("planted")\n'
    )
    with tempfile.TemporaryDirectory() as box:
        planted = Path(box) / "test_planted.py"
        planted.write_text(source)

        def run(*extra: str) -> str:
            done = subprocess.run(
                [sys.executable, "-m", "pytest", "-q", "-rs", "-p", "no:cacheprovider", *extra, "."],
                cwd=box,
                capture_output=True,
                text=True,
            )
            return done.stdout + done.stderr

        folded = run()
        unfolded = run("--no-fold-skipped")

    assert "3 skipped" in folded, f"the plant did not skip three times:\n{folded}"
    assert "test_planted_skip_folded[1]" not in folded, (
        "pytest did not fold the skip summary here, so the difference this guard rests on does not "
        f"exist in this pytest and the guard has to be rewritten rather than trusted:\n{folded}"
    )
    for case in (1, 2, 3):
        assert f"test_planted_skip_folded[{case}]" in unfolded, (
            f"--no-fold-skipped did not name test_planted_skip_folded[{case}]. ci.yml passes this "
            "flag so that the ids of what CI skipped are on the record; if the flag no longer does "
            f"that, ci.yml is carrying a promise it cannot keep:\n{unfolded}"
        )


def test_this_file_exempts_nothing() -> None:
    """No allowlist, no threshold, no name-based escape. Every assertion above reads bytes or git."""
    import ast

    tree = ast.parse(Path(__file__).read_text())
    offenders = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = [node.target] if isinstance(node, ast.AnnAssign) else node.targets
            offenders += [
                ast.unparse(t)
                for t in targets
                # Whole words: `shallow` contains "allow" and this test once reported itself.
                if re.search(r"\b(allow\w*|exempt\w*|ignore\w*|skip\w*)\b", ast.unparse(t).lower())
            ]
        if isinstance(node, ast.FunctionDef):
            offenders += [
                ast.unparse(d)
                for d in node.decorator_list
                if "skip" in ast.unparse(d) and "parametrize" not in ast.unparse(d)
            ]
    assert offenders == [], offenders
