# SPDX-License-Identifier: AGPL-3.0-or-later
"""A cleanliness block that meets the whole contract, built in a planted repository.

Why this exists (2026-10-02, lane-entrypoint). Since the entry script is counted
(`genomeos.manifest.entry_script`), `sys.argv[0]` under pytest is a test runner and no form but a
tracked, unmodified file inside the repository is clean, so NO in-process call to
`code_cleanliness` can produce a passing block while the suite is running. Four fixtures needed one
only as scaffolding -- they were testing the allowlist, the revision race and the key-set refusal, not
cleanliness -- and each built it by naming a file of THIS shared checkout, which is the pattern
`tests/test_code_cleanliness_hermetic.py` was written to remove: a verdict that depends on what some
lane happens to be editing.

So the block is built in a repository planted for the purpose, with one committed script, and the
script is both the declared entry and the `argv0`. Nothing here weakens a check: the block is a real
`code_cleanliness` return value over a real git repository in which the entry script really is
committed and clean, which is the state the contract asks for. A fixture that merely asserted the
keys would have passed whatever the function did.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any

from genomeos import manifest as mf

#: Under `scripts/`, because `mf.is_code` counts a path as code only under CODE_ROOTS: a planted
#: script at a repository root would be invisible to the revision stamp and the block would be clean
#: for the wrong reason.
SCRIPT = "scripts/planted_clean_writer.py"

SOURCE = "VALUE = 1\n"


def planted_clean_block(tmp_path: Path, name: str = "planted_clean") -> dict[str, Any]:
    """A `code_cleanliness` block over a planted repository whose entry script is committed and clean."""
    repo = tmp_path / name
    (repo / "scripts").mkdir(parents=True, exist_ok=True)
    (repo / SCRIPT).write_text(SOURCE)
    if not (repo / ".git").exists():
        subprocess.run(["git", "init", "-q", str(repo)], check=True, capture_output=True)
    git = ["git", "-C", str(repo), "-c", "user.email=p@p", "-c", "user.name=p"]
    subprocess.run([*git, "add", "--", SCRIPT], check=True, capture_output=True)
    # Idempotent: a caller may ask for the same planted repository twice in one test (the allowlist
    # exemption writes the result before and after changing the list), and `git commit` with nothing
    # staged exits 1. Committing only when something is staged keeps the second call a no-op rather
    # than an error, and the state the block is read from is the same either way.
    staged = subprocess.run(
        [*git, "diff", "--cached", "--name-only"], check=True, capture_output=True, text=True
    )
    if staged.stdout.strip():
        subprocess.run([*git, "commit", "-q", "-m", "planted"], check=True, capture_output=True)
    return mf.code_cleanliness(SCRIPT, (SCRIPT,), repo, argv0=SCRIPT)
