# SPDX-License-Identifier: AGPL-3.0-or-later
"""The installed pre-push hook must BE `scripts/pre-push-hook.sh`, byte for byte.

WHY. `scripts/install-hooks.sh` COPIES a file into `.git/hooks/pre-push`, so the hook is not the
script and an edit to the script changes nothing the push runs. On 2026-10-02 the installed copy was
taken at 04:34 and still symlinked the git-ignored data stores into the push's verification worktree,
while `scripts/pre-push.sh` had moved at 17:52 to read-only APFS clones. `git check-ignore` cannot
answer past a symbolic link, 20 tests ERRORED at setup, and two RED status files were written for
trees that were nothing of the kind -- one of them refusing a push and leaving origin/dev thirty
commits behind. The clone switch was correct and could not have fixed this: it never reached what the
push ran.

The answer is a mechanism and not a note. What is installed is now a FIXED WRAPPER that extracts
`HEAD:scripts/pre-push.sh` and runs that, and refuses when the working copy differs from HEAD's -- so
an edit takes effect when it is committed, and the wrapper itself never needs re-copying. This test is
the part that notices when a checkout is nevertheless carrying something else.

IT SKIPS BY NAME WHERE NO HOOK IS INSTALLED, which is CI and a fresh clone: there is nothing to push
from and nothing to check, and a skip that names the file and the installer tells a reader "not run
here" rather than "passed". It FAILS where a hook IS installed and is not the wrapper, because that
checkout can push and its pushes are running something nobody can name.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

WRAPPER = Path("scripts/pre-push-hook.sh")
INSTALLER = "scripts/install-hooks.sh"


def installed_hook() -> Path:
    """`.git/hooks/pre-push` for this checkout, the COMMON git dir so a worktree sees the real one."""
    common = subprocess.run(
        ["git", "rev-parse", "--path-format=absolute", "--git-common-dir"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    return Path(common) / "hooks" / "pre-push"


def test_the_installed_hook_is_the_wrapper_byte_for_byte() -> None:
    hook = installed_hook()
    if not hook.exists():
        pytest.skip(
            f"NOT RUN HERE, not passed: no pre-push hook is installed at {hook}, so this checkout "
            f"does not push and there is nothing to compare; `{INSTALLER}` installs it"
        )
    assert hook.read_bytes() == WRAPPER.read_bytes(), (
        f"{hook} is not {WRAPPER}. A hook is a COPY, so it goes stale silently and the push then runs "
        f"a checker nobody can name -- which on 2026-10-02 symlinked the data stores and produced two "
        f"red verdicts about nothing. Run `{INSTALLER}`"
    )


def test_the_wrapper_is_executable_and_runs_the_committed_checker_not_the_working_copy() -> None:
    """The two properties the wrapper exists for, read off its own text so neither can be dropped."""
    assert WRAPPER.stat().st_mode & 0o111, f"{WRAPPER} is not executable"
    text = WRAPPER.read_text()
    assert 'git show "HEAD:$script"' in text, "the wrapper must run the COMMITTED checker"
    assert "git hash-object" in text, "the wrapper must compare the working copy against HEAD"
    assert "REFUSED" in text, "a difference must be refused by name, not reported and continued"


def test_the_installer_installs_the_wrapper_and_not_the_checker() -> None:
    """The regression itself: installing pre-push.sh directly is what went stale."""
    text = Path(INSTALLER).read_text()
    assert "cp scripts/pre-push-hook.sh .git/hooks/pre-push" in text, text
    assert "cp scripts/pre-push.sh .git/hooks/pre-push" not in text, (
        "the installer is copying the CHECKER into .git/hooks again; that copy is what went stale"
    )


def test_the_wrapper_records_what_it_does_not_close() -> None:
    """A guard read as closing more than it does is worse than none. The trade-off must be written in.

    The hook runs the pre-push.sh of the commit being pushed, so a commit CAN weaken its own check.
    That case is the commit review's -- it is a visible diff in a tracked file. The wrapper closes the
    SILENT case, where nobody changed anything and the push ran last week's checker.
    """
    text = WRAPPER.read_text()
    assert "WHAT IT DOES NOT CLOSE" in text
    assert "SILENT" in text


def test_the_wrapper_refuses_an_uncommitted_checker(tmp_path: Path) -> None:
    """Measured, not asserted from the text: a modified working copy must stop the push by name."""
    repo = tmp_path / "repo"
    (repo / "scripts").mkdir(parents=True)
    checker = repo / "scripts" / "pre-push.sh"
    checker.write_text("#!/usr/bin/env bash\necho THE-COMMITTED-CHECKER-RAN\nexit 0\n")
    (repo / "scripts" / "pre-push-hook.sh").write_bytes(WRAPPER.read_bytes())

    def git(*args: str) -> subprocess.CompletedProcess:
        return subprocess.run(["git", *args], cwd=str(repo), capture_output=True, text=True)

    git("init", "-q", "-b", "main", ".")
    git("add", "scripts/pre-push.sh")
    tree = git("write-tree").stdout.strip()
    sha = git(
        "-c", "user.email=t@example.invalid", "-c", "user.name=t", "commit-tree", tree, "-m", "init"
    ).stdout.strip()
    git("update-ref", "refs/heads/main", sha, "")

    clean = subprocess.run(
        ["bash", "scripts/pre-push-hook.sh", "origin", "https://example.invalid/x.git"],
        cwd=str(repo),
        input="",
        capture_output=True,
        text=True,
    )
    assert clean.returncode == 0, clean.stderr
    assert "THE-COMMITTED-CHECKER-RAN" in clean.stdout, clean

    checker.write_text("#!/usr/bin/env bash\nexit 0   # a draft nobody committed\n")
    dirty = subprocess.run(
        ["bash", "scripts/pre-push-hook.sh", "origin", "https://example.invalid/x.git"],
        cwd=str(repo),
        input="",
        capture_output=True,
        text=True,
    )
    assert dirty.returncode == 1, dirty
    assert "REFUSED" in dirty.stderr and "scripts/pre-push.sh" in dirty.stderr, dirty.stderr
    assert "no verdict" in dirty.stderr, dirty.stderr
    assert "THE-COMMITTED-CHECKER-RAN" not in dirty.stdout


def test_the_wrapper_passes_stdin_through_to_the_committed_checker(tmp_path: Path) -> None:
    """git sends four fields per ref ON STDIN. A wrapper that ate them would check nothing, quietly."""
    repo = tmp_path / "repo"
    (repo / "scripts").mkdir(parents=True)
    (repo / "scripts" / "pre-push.sh").write_text(
        '#!/usr/bin/env bash\nwhile read -r a b c d; do echo "GOT $a $b $c $d"; done\n'
    )
    (repo / "scripts" / "pre-push-hook.sh").write_bytes(WRAPPER.read_bytes())

    def git(*args: str) -> subprocess.CompletedProcess:
        return subprocess.run(["git", *args], cwd=str(repo), capture_output=True, text=True)

    git("init", "-q", "-b", "main", ".")
    git("add", "scripts/pre-push.sh")
    tree = git("write-tree").stdout.strip()
    sha = git(
        "-c", "user.email=t@example.invalid", "-c", "user.name=t", "commit-tree", tree, "-m", "init"
    ).stdout.strip()
    git("update-ref", "refs/heads/main", sha, "")

    r = subprocess.run(
        ["bash", "scripts/pre-push-hook.sh", "origin", "https://example.invalid/x.git"],
        cwd=str(repo),
        input="refs/heads/dev aaa refs/heads/dev bbb\n",
        capture_output=True,
        text=True,
    )
    assert r.returncode == 0, r.stderr
    assert "GOT refs/heads/dev aaa refs/heads/dev bbb" in r.stdout, r.stdout
