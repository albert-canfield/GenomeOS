"""Check docs/ROADMAP.md's current-state claims against machine-checkable values.

A narrative sentence about a number is not the number. Prose in an area section is
evidence about its own date; the value lives in a test pin, a committed result under
``data/results/``, a constant in a module, or git itself. This script checks the three
classes where the claim's value is mechanically locatable, so the pass is re-runnable
instead of a reading someone has to repeat:

1. **Cited commits.** Every backticked 7-to-10 hex token in the document that names a
   commit in this repository must be an ancestor of HEAD. A commit that exists but is
   not an ancestor is either on an unmerged branch or was relabelled before a push, and
   the row citing it is making a claim the history does not carry. Tokens that are not
   commits here at all (external repositories, tree and digest prefixes) are reported
   separately rather than counted as failures, because they are not claims about this
   history.

2. **Result names in the ``- **Data.**`` bullets.** Each area's Data bullet is an
   inventory of committed results. Every plain backticked name in it must exist under
   ``data/results/``; every ``name_chr*`` glob is counted on disk and the count printed
   beside the count the prose states, so a row that understates its own area is visible.

3. **Milestones.** ``roadmap.parse_milestones`` is the parser the web status card uses,
   so the ✅/◑ markers have one reading. The script prints the states and the reached
   count for comparison with any prose that quotes it.

One limitation to read the output with: the stated count taken is the first one after
the name, and area C's Data row deliberately quotes the figure it superseded before
giving the live one, so that row reports a difference it has already corrected in prose.

Counts are reported, never asserted: the prose figure was correct on its own date and a
drifted figure is stale rather than wrong, so the additive correction belongs to whoever
owns the row. What *is* asserted is existence and ancestry, which cannot be stale and
only ever indicate a defect.

Run: ``uv run python scripts/roadmap_truth_pass.py`` (``--json`` for the machine form).
Exit status is 1 when an asserted invariant fails.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ROADMAP = ROOT / "docs" / "ROADMAP.md"
RESULTS = ROOT / "data" / "results"

SHA = re.compile(r"`([0-9a-f]{7,10})`")
BACKTICKED = re.compile(r"`([^`]+)`")
DATA_BULLET = re.compile(r"^- \*\*Data\.\*\*")
# The count a Data row states for a family, in either form the document uses:
# "`budget_chr*` (24 of 24)" and "`proteome_chr*` **25 of 25**".
STATED_COUNT = re.compile(r"^[\s,;]*(?:\(|\*\*)\s*(\d+)(?:\s+of\s+(\d+))?")

# Commits that are deliberately not ancestors of dev's HEAD, with the reason. Each is a
# claim the document makes accurately; the checker would otherwise report them forever.
NOT_ANCESTOR_EXPECTED = {
    "7ff4e37": "main's merge commit for PR #11 (promotion), not on dev",
    "f7016ff": "main's merge commit for PR #12 (promotion), not on dev",
    "a1ce607": "main's merge commit for PR #13 (promotion), not on dev",
    "9c64371": "merge of main into promote-bdee236, not on dev",
    "054eb46": "merge of main into promote-747b7ff, not on dev",
    "cc677e8": "relabelled to a3f1c04 before it was pushed; the citing row says so",
}


def _git(*args: str) -> tuple[int, str]:
    proc = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=False)
    return proc.returncode, proc.stdout.strip()


def check_commits(text: str) -> dict:
    """Every cited commit of this repository must be an ancestor of HEAD."""
    cited = sorted(set(SHA.findall(text)))
    ancestors: list[str] = []
    not_commits: list[str] = []
    stray: list[dict] = []
    for sha in cited:
        if _git("cat-file", "-e", f"{sha}^{{commit}}")[0] != 0:
            not_commits.append(sha)
            continue
        if _git("merge-base", "--is-ancestor", sha, "HEAD")[0] == 0:
            ancestors.append(sha)
            continue
        subject = _git("log", "-1", "--format=%ad %s", "--date=short", sha)[1]
        stray.append(
            {
                "commit": sha,
                "commit line": subject,
                "expected": NOT_ANCESTOR_EXPECTED.get(sha),
            }
        )
    return {
        "cited": len(cited),
        "ancestors of HEAD": len(ancestors),
        "not commits here": not_commits,
        "not ancestors": stray,
        "unexplained": [s for s in stray if s["expected"] is None],
    }


def _data_bullets(lines: list[str]) -> list[dict]:
    """Each ``- **Data.**`` bullet with its line number and its continuation lines."""
    out = []
    for i, line in enumerate(lines):
        if not DATA_BULLET.match(line):
            continue
        body = [line]
        for nxt in lines[i + 1 :]:
            if nxt.startswith("- **") or not nxt.startswith("  "):
                break
            body.append(nxt)
        out.append({"line": i + 1, "text": "\n".join(body)})
    return out


CHROM = r"chr(?:[0-9]+|X|Y|M)"


def _on_disk(stem: str) -> tuple[int, int]:
    """How many chromosomes ``stem`` covers on disk, and how many files carry it.

    The prose glob ``unknown_chr*`` means one result per chromosome, so the first
    number counts files named exactly ``<stem>chrN.json`` — not ``rmsk_chr1.bed.gz``
    beside its JSON, and not ``domains_chr21_hic_K562.json``, which is a different
    family under the same prefix. The second number is every file with the prefix, so
    a family that has grown variants is visible rather than silently equal.
    """
    base = re.sub(r"_?chr$", "", stem)
    exact = re.compile(rf"^{re.escape(base)}_?({CHROM})\.[a-z0-9.]+$")
    prefix = re.compile(rf"^{re.escape(base)}_?{CHROM}")
    names = [p.name for p in RESULTS.iterdir()]
    # Distinct chromosomes, not files: rmsk ships a JSON and a BED per chromosome, and
    # counting files would read 50 where the prose means 25.
    covered = {m.group(1) for n in names if (m := exact.match(n))}
    return len(covered), sum(1 for n in names if prefix.match(n))


def check_data_bullets(lines: list[str]) -> list[dict]:
    """Named results must exist; globbed families are counted beside the stated count."""
    rows = []
    for bullet in _data_bullets(lines):
        text = bullet["text"]
        for name in BACKTICKED.findall(text):
            if "/" in name or " " in name or name.endswith(".bio"):
                continue  # a path or a prose aside, not a result name
            after = text.split(f"`{name}`", 1)[1][:24]
            stated = STATED_COUNT.search(after)
            if name.endswith("*"):
                chromosomes, files = _on_disk(name[:-1])
                rows.append(
                    {
                        "line": bullet["line"],
                        "claim": name,
                        "chromosomes on disk": chromosomes,
                        "files with the prefix": files,
                        "stated": int(stated.group(1)) if stated else None,
                        "kind": "family",
                    }
                )
            else:
                rows.append(
                    {
                        "line": bullet["line"],
                        "claim": name,
                        "exists": (RESULTS / f"{name}.json").exists(),
                        "kind": "file",
                    }
                )
    return rows


def check_milestones(text: str) -> dict:
    sys.path.insert(0, str(ROOT))
    from genomeos import roadmap  # noqa: PLC0415  (import after the path is set)

    parsed = roadmap.parse_milestones(text)
    return {
        "milestones": [{"milestone": m.get("milestone"), "state": m.get("state")} for m in parsed],
        "done": sum(1 for m in parsed if m.get("state") == "done"),
        "total": len(parsed),
    }


def run() -> dict:
    text = ROADMAP.read_text()
    lines = text.splitlines()
    return {
        "commits": check_commits(text),
        "data bullets": check_data_bullets(lines),
        "milestones": check_milestones(text),
    }


def report(found: dict) -> int:
    commits = found["commits"]
    print(f"commits cited {commits['cited']}, ancestors of HEAD {commits['ancestors of HEAD']}")
    for stray in commits["not ancestors"]:
        why = stray["expected"] or "UNEXPLAINED — the history does not carry this claim"
        print(f"  not an ancestor: {stray['commit']}  {stray['commit line']}  [{why}]")
    if commits["not commits here"]:
        print(
            "  not commits in this repository (external shas, tree and digest "
            f"prefixes): {', '.join(commits['not commits here'])}"
        )

    print("\nData bullets")
    missing = []
    for row in found["data bullets"]:
        if row["kind"] == "file":
            mark = "ok" if row["exists"] else "MISSING"
            if not row["exists"]:
                missing.append(row)
            print(f"  L{row['line']:>5}  {row['claim']:<44} {mark}")
        else:
            stated = "-" if row["stated"] is None else str(row["stated"])
            chroms = row["chromosomes on disk"]
            drift = "" if row["stated"] in (None, chroms) else "   <- differs"
            print(
                f"  L{row['line']:>5}  {row['claim']:<44} "
                f"stated {stated:>4}  chromosomes {chroms:>4}  "
                f"files {row['files with the prefix']:>4}{drift}"
            )

    ms = found["milestones"]
    print(f"\nMilestones: {ms['done']} done of {ms['total']}")
    for m in ms["milestones"]:
        print(f"  {m['state']:<8} {m['milestone']}")

    failed = bool(commits["unexplained"]) or bool(missing)
    print(
        "\nasserted invariants: "
        + ("FAILED" if failed else "hold")
        + " (counts are reported, not asserted: a drifted count is stale rather than "
        "wrong, and the correction belongs to the row's owner)"
    )
    return 1 if failed else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="the machine form")
    args = parser.parse_args()
    found = run()
    if args.json:
        print(json.dumps(found, indent=2))
        return 1 if found["commits"]["unexplained"] else 0
    return report(found)


if __name__ == "__main__":
    raise SystemExit(main())
