# SPDX-License-Identifier: AGPL-3.0-or-later
"""Stage one section of a shared markdown file: the HEAD copy, with only your section spliced in.

    GIT_INDEX_FILE=/tmp/idx-$$ git read-tree HEAD
    uv run python scripts/stage_section.py docs/ATTRIBUTION.md "## My heading" "## The next heading"

Several sessions write this checkout at once, so a shared file's working copy holds other lanes'
half-finished text as well as yours. `git add docs/ROADMAP.md` stages all of it, and twice now a lane
has found its unfinished wording published in another lane's commit under another lane's message.

This removes the hazard structurally rather than by care. The base is the file as the index already
holds it (`git read-tree HEAD` puts HEAD's copy there), the only edit applied to it is the one section
whose heading you name, lifted from the working copy, and the result is written straight into your
private index with `hash-object` and `update-index`. Nothing outside your section can reach the index,
because everything outside it is the base verbatim. Run it once per section you own; each run builds on
the last, so two sections of the same file compose.

A section runs from its heading to the next heading of the same level or shallower, so naming a `###`
subsection takes that subsection and not the rest of its parent. Headings inside fenced code blocks are
not headings.

It refuses rather than guesses: if your heading is not in the working copy there is nothing to stage,
and if it is not in the base it has no home, so you must name the heading it goes before. A refusal
stages nothing and exits 2.

The third argument is optional and only used when the section is new to the base.
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path


class RefusedError(Exception):
    """Staging nothing is always safer than staging the wrong thing."""


def _level(line: str) -> int:
    """Heading depth, or 0 for a line that is not a heading. `### A. BioLang` is 3."""
    n = len(line) - len(line.lstrip("#"))
    return n if 1 <= n <= 6 and line[n : n + 1] in (" ", "\t") else 0


def _headings(lines: list[str]) -> list[int]:
    """Heading depth per line, counting nothing inside a fenced code block: a shell comment at the
    start of a line in a ``` block is not a heading, and treating one as a heading would end a
    section early and drop the rest of it."""
    out, fenced = [], False
    for line in lines:
        if line.lstrip().startswith(("```", "~~~")):
            fenced = not fenced
            out.append(0)
            continue
        out.append(0 if fenced else _level(line))
    return out


def _bounds(lines: list[str], heading: str) -> tuple[int, int] | None:
    """The half-open line range of the `heading` section: from its heading to the next heading of the
    same level or shallower.

    Ending at the next `## ` regardless of depth was wrong, and wrong in the dangerous direction. Every
    roadmap area is a `### ` heading, so a lane staging `### E. From one cell to an organism` got
    everything from E to the next `## ` — areas E through J in one splice, nine other lanes' text
    carried out under its commit, which is the exact hazard this tool exists to prevent.
    """
    depth = _level(heading)
    if not depth:
        raise RefusedError(f"{heading!r} is not a markdown heading: it must start with # and a space")
    levels = _headings(lines)
    start = next((i for i, x in enumerate(lines) if levels[i] and x.startswith(heading)), None)
    if start is None:
        return None
    end = next(
        (j for j in range(start + 1, len(lines)) if 0 < levels[j] <= depth),
        len(lines),
    )
    return start, end


def splice(base: list[str], work: list[str], heading: str, anchor: str | None = None) -> list[str]:
    """`base` with the `heading` section replaced by the working copy's, and nothing else changed."""
    mine = _bounds(work, heading)
    if mine is None:
        raise RefusedError(f"no section starting {heading!r} in the working copy: nothing of yours to stage")
    body = work[mine[0] : mine[1]]

    there = _bounds(base, heading)
    if there is not None:
        return base[: there[0]] + body + base[there[1] :]
    if not anchor:
        raise RefusedError(
            f"{heading!r} is new to the committed file; name the heading it goes before, so it is not "
            "appended to the end of a file another lane is also writing"
        )
    levels = _headings(base)
    at = next((i for i, x in enumerate(base) if levels[i] and x.startswith(anchor)), None)
    if at is None:
        raise RefusedError(f"the anchor {anchor!r} is not in the committed file")
    return base[:at] + body + base[at:]


def _base(path: str) -> list[str]:
    """The file as the index holds it, falling back to HEAD when it is not staged yet."""
    for args in (["git", "cat-file", "-p", f":{path}"], ["git", "show", f"HEAD:{path}"]):
        done = subprocess.run(args, capture_output=True, text=True)
        if done.returncode == 0:
            return done.stdout.splitlines(keepends=True)
    raise RefusedError(f"{path} is in neither the index nor HEAD")


def stage(path: str, heading: str, anchor: str | None = None) -> str:
    out = splice(_base(path), Path(path).read_text().splitlines(keepends=True), heading, anchor)
    with tempfile.NamedTemporaryFile("w", suffix=".md", delete=False) as f:
        f.writelines(out)
        tmp = f.name
    try:
        sha = subprocess.run(
            ["git", "hash-object", "-w", tmp], capture_output=True, text=True, check=True
        ).stdout.strip()
        subprocess.run(["git", "update-index", "--add", "--cacheinfo", f"100644,{sha},{path}"], check=True)
    finally:
        Path(tmp).unlink(missing_ok=True)
    return sha


def main(argv: list[str]) -> int:
    if not 3 <= len(argv) <= 4:
        print(__doc__, file=sys.stderr)
        return 2
    path, heading, anchor = argv[1], argv[2], argv[3] if len(argv) > 3 else None
    try:
        stage(path, heading, anchor)
    except RefusedError as e:
        print(f"refused: {e}", file=sys.stderr)
        return 2
    shown = subprocess.run(
        ["git", "diff", "--cached", "--stat", "HEAD", "--", path], capture_output=True, text=True
    ).stdout.strip()
    print(f"staged {path}: HEAD plus {heading!r} alone\n{shown or '  (no change)'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
