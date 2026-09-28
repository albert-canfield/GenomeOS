# SPDX-License-Identifier: Apache-2.0
"""Which `when` values on rules and events use a construct the old equality test could not match.

Decisions and timers read `when` with `genomeos.ir.model.matches`, which knows `any`, `unknown`,
`absent`, alternatives `a|b` and the comparisons `>=n` `<=n` `>n` `<n`. Until 2026-09-28 `Rule.applies`
and `Event.applies` compared by plain equality (with `any` and, for rules, `unknown` special-cased),
so a rule written `when: cell_type = K562|HepG2` compiled and never fired. This census lists every
rule and event `when` clause in the repo's BioLang programs, with the program sources of
scripts/biolang_key_census.py (tracked .bio, fenced Markdown, string literals in code, the
untracked compiled chromosomes, and with `--tests DIR` every program the test suite parses), and
says for each clause whether the old test failed it silently.

Writes data/results/when_census.json.
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import biolang_key_census as bkc  # noqa: E402

from genomeos import manifest as mf  # noqa: E402
from genomeos.lang import parser as lang_parser  # noqa: E402
from genomeos.results import save_result  # noqa: E402

NAME = "when_census"
KINDS = ("rule", "event")  # the two blocks whose `applies` compared by equality


def construct(value: str) -> str:
    """The grammar construct a `when` value uses (docs/BIOLANG-GRAMMAR.md, lexical rules)."""
    if value == "any":
        return "any"
    if value == "unknown":
        return "unknown"
    if value == "absent":
        return "absent"
    if value[:2] == "!=":
        return "not_equal"  # the parser accepts `k != v` and stores `!=v`; the matcher has no case for it
    if value[:2] in (">=", "<=") or value[:1] in (">", "<"):
        return "comparison"
    if "|" in value:
        return "alternatives"
    return "plain"


# what the equality test did with each construct, against what `matches` does
SILENT = {
    "plain": False,
    "any": False,
    "unknown": False,  # rules: special-cased by R1; events: equality, but Event.applies has no caller
    "absent": True,  # equality matched only a context whose value is the string "absent"
    "alternatives": True,  # matched only a context whose value is the string "a|b"
    "comparison": True,  # matched only a context whose value is the string ">=n"
    "not_equal": True,  # neither test implements it: both match only the string "!=v"
}


def when_hits(text: str) -> list[dict]:
    _, top = lang_parser._parse_blocks(text.splitlines())
    hits = []

    def walk(b: lang_parser.Block) -> None:
        if b.kind in KINDS and "when" in b.props:
            try:
                when = lang_parser._parse_when(b.props["when"])
            except lang_parser.BioLangError as e:
                hits.append({"kind": b.kind, "block": b.header, "line": b.line, "parse_error": str(e)})
                when = {}
            for key, value in when.items():
                c = construct(value)
                hits.append(
                    {
                        "kind": b.kind,
                        "block": b.header,
                        "line": b.line,
                        "key": key,
                        "value": value,
                        "construct": c,
                        "failed_silently_before": SILENT[c],
                    }
                )
        for c in b.children:
            walk(c)

    for b in top:
        walk(b)
    return hits


def summarise(hits: list[dict], programs: int) -> dict:
    """Counts for every clause; only the clauses that are not plain equality are listed."""
    by = Counter((h["kind"], h["construct"]) for h in hits if "construct" in h)
    special = [h for h in hits if h.get("construct") != "plain"]
    return {
        "programs": programs,
        "clauses": sum(by.values()),
        "by_kind_and_construct": {f"{k}/{c}": n for (k, c), n in sorted(by.items())},
        "non_plain_hits": special,
    }


def main(argv: list[str]) -> int:
    bkc.unknown_keys = when_hits  # the key census's source walkers, with this census's check
    tracked = [
        ROOT / p
        for p in bkc.subprocess.run(
            ["git", "ls-files", "*.bio"], cwd=ROOT, capture_output=True, text=True
        ).stdout.split()
    ]
    generated = sorted((ROOT / "data/knowledge/compiled").glob("*.bio"))
    parts: dict = {}
    parts["tracked_bio"] = summarise(*bkc.census_files(tracked))
    parts["markdown_fenced"] = summarise(*bkc.census_markdown())
    parts["code_literals"] = summarise(*bkc.census_code())
    parts["generated_bio_untracked"] = summarise(*bkc.census_files(generated))
    if "--tests" in argv:
        scratch = Path(argv[argv.index("--tests") + 1])
        h, n, pytest_summary = bkc.census_tests(scratch)
        parts["test_programs"] = {**summarise(h, n), "pytest": pytest_summary}
    silent = [
        {"source_part": name, **h}
        for name, part in parts.items()
        for h in part["non_plain_hits"]
        if h.get("failed_silently_before")
    ]
    payload = {
        "question": "which rule and event `when` clauses used a construct Rule/Event.applies could not match",
        "matcher": "genomeos.ir.model.matches (decisions, timers, competence windows, commitments, signals)",
        "grammar_constructs": {
            "plain": "k = v, equality",
            "any": "matches anything, present or absent",
            "unknown": "the context was not recorded: matches no context (R1, 2026-09-28)",
            "absent": "the key is not in the context",
            "alternatives": "a|b, any listed value",
            "comparison": ">=n <=n >n <n, numeric; `k >= n` is parsed as `k = >=n`",
            "not_equal": "`k != v` is parsed to `!=v` but neither the grammar nor the matcher defines it",
        },
        "old_behaviour": {
            "Rule.applies": "equality; `any` passes; `unknown` fails (R1)",
            "Event.applies": "equality; `any` passes; no caller in genomeos, scripts or tests",
        },
        **parts,
        "silent_failures_before": silent,
        "silent_failure_count": len(silent),
    }
    manifest = {
        "sources": [
            {
                "accession": "git ls-files *.bio, *.md fences, *.py/*.js/*.html literals",
                "version": "HEAD of dev",
            },
            {"accession": "data/knowledge/compiled/*.bio", "version": "untracked, as built in this checkout"},
        ],
        "inputs": [mf.input_entry("genomeos/lang/parser.py"), mf.input_entry("genomeos/ir/model.py")],
        "assembly": "n/a: program text, no genome coordinates",
        "coordinates": "n/a: no genomic intervals",
        "parameters": {"kinds": list(KINDS), "tests": "--tests" in argv},
        "exclusions": [],
        "partitions": "n/a: no evaluation split",
    }
    save_result(NAME, payload, manifest=manifest)
    for name, part in parts.items():
        print(name, part["programs"], "programs", part["clauses"], "clauses", part["by_kind_and_construct"])
    print("silent failures before:", len(silent))
    for h in silent:
        print("  ", h["source_part"], h.get("source"), h["kind"], h["block"], h["key"], "=", h["value"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
