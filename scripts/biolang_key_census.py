# SPDX-License-Identifier: Apache-2.0
"""Which `key: value` lines in the repo's BioLang programs name a key their block does not read.

Before this census the parser stored every key it was given and read only the ones it knew, so a
misspelt key (`basal_rate:` for `basal:`) compiled to a program that silently lacked it. The census
lexes every program with the parser's own block reader and checks each key against the grammar
table (`genomeos.lang.grammar.BLOCKS` plus the keys every block takes). Sources:

- every tracked `.bio` file;
- every fenced code block in the repo's Markdown that lexes into BioLang blocks;
- every string literal in tracked Python, and every backtick template literal in tracked JavaScript
  and HTML, that lexes into BioLang blocks (smoke programs, REPL inputs, web examples);
- with `--tests`, every program the test suite hands the parser (recorded by wrapping the block
  reader for one pytest run), because most test programs are inline strings;
- untracked generated `.bio` files (data/knowledge/compiled), reported separately.

Writes data/results/biolang_key_census.json.
"""

from __future__ import annotations

import difflib
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos.lang import grammar  # noqa: E402
from genomeos.lang import parser as lang_parser  # noqa: E402

OUT = ROOT / "data/results/biolang_key_census.json"  # written through save_result (item 12 S6)
FENCE = re.compile(r"^```[^\n]*\n(.*?)^```", re.S | re.M)


def accepted(kind: str) -> set[str]:
    return set(grammar.BLOCKS[kind].get("props", {})) | set(grammar.COMMON)


def unknown_keys(text: str) -> list[dict]:
    _, top = lang_parser._parse_blocks(text.splitlines())
    hits = []

    def walk(b: lang_parser.Block) -> None:
        for key in b.props:
            if key not in accepted(b.kind):
                close = difflib.get_close_matches(key, sorted(accepted(b.kind)), n=1)
                hits.append(
                    {
                        "kind": b.kind,
                        "block": b.header,
                        "line": b.line,
                        "key": key,
                        "closest": close[0] if close else "",
                    }
                )
        for c in b.children:
            walk(c)

    for b in top:
        walk(b)
    return hits


def census_files(paths: list[Path]) -> tuple[list[dict], int]:
    hits, n = [], 0
    for path in paths:
        try:
            found = unknown_keys(path.read_text())
        except lang_parser.BioLangError as e:
            hits.append({"source": str(path.relative_to(ROOT)), "lex_error": str(e)})
            continue
        n += 1
        hits += [{"source": str(path.relative_to(ROOT)), **h} for h in found]
    return hits, n


def census_markdown() -> tuple[list[dict], int]:
    hits, n = [], 0
    files = subprocess.run(
        ["git", "ls-files", "*.md"], cwd=ROOT, capture_output=True, text=True
    ).stdout.split()
    for name in files:
        for m in FENCE.finditer((ROOT / name).read_text(errors="replace")):
            try:
                _, top = lang_parser._parse_blocks(m.group(1).splitlines())
            except lang_parser.BioLangError:
                continue
            if not top:
                continue
            n += 1
            hits += [{"source": f"{name} (fenced block)", **h} for h in unknown_keys(m.group(1))]
    return hits, n


TEMPLATE = re.compile(r"`([^`]*\{[^`]*)`", re.S)


def _lexes(text: str) -> bool:
    try:
        _, top = lang_parser._parse_blocks(text.splitlines())
    except lang_parser.BioLangError:
        return False
    return bool(top)


def census_code() -> tuple[list[dict], int]:
    import ast

    hits, n = [], 0
    names = subprocess.run(
        ["git", "ls-files", "*.py", "*.js", "*.html"], cwd=ROOT, capture_output=True, text=True
    ).stdout.split()
    for name in names:
        text = (ROOT / name).read_text(errors="replace")
        if name.endswith(".py"):
            try:
                tree = ast.parse(text)
            except SyntaxError:
                continue
            found = [
                (node.lineno, node.value)
                for node in ast.walk(tree)
                if isinstance(node, ast.Constant) and isinstance(node.value, str) and "{" in node.value
            ]
        else:
            found = [(text.count("\n", 0, m.start()) + 1, m.group(1)) for m in TEMPLATE.finditer(text)]
        for lineno, program in found:
            if "{{" in program:  # a str.format template: undouble the braces, fill each field with 1
                program = re.sub(r"(?<!\{)\{\w*\}(?!\})", "1", program).replace("{{", "{").replace("}}", "}")
            # a literal that does not lex whole (a REPL script with a deliberately broken line) is
            # read line by line, so its one-line blocks are still checked
            pieces = [program] if _lexes(program) else [x for x in program.splitlines() if _lexes(x)]
            for piece in pieces:
                n += 1
                hits += [{"source": f"{name}:{lineno}", **h} for h in unknown_keys(piece)]
    return hits, n


PLUGIN = """
import json, os
from genomeos.lang import parser as lang_parser
_orig = lang_parser._parse_blocks
_seen = []
def _wrapped(lines):
    directives, top = _orig(lines)
    _seen.append("\\n".join(lines))
    return directives, top
lang_parser._parse_blocks = _wrapped
def pytest_runtest_setup(item):
    _seen.append("@@" + item.nodeid)
def pytest_sessionfinish(session):
    with open(os.environ["BIOLANG_CENSUS_OUT"], "w") as f:
        json.dump(_seen, f)
"""


def census_tests(scratch: Path) -> tuple[list[dict], int, str]:
    scratch.mkdir(parents=True, exist_ok=True)
    (scratch / "biolang_census_plugin.py").write_text(PLUGIN)
    out = scratch / "seen.json"
    env = {**__import__("os").environ, "BIOLANG_CENSUS_OUT": str(out), "PYTHONPATH": str(scratch)}
    r = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "-p",
            "biolang_census_plugin",
            "--continue-on-collection-errors",
            "--no-header",
            "-p",
            "no:cacheprovider",
            "tests",
        ],
        cwd=ROOT,
        env=env,
        capture_output=True,
        text=True,
    )
    summary = (r.stdout.strip().splitlines() or [""])[-1]
    seen = json.loads(out.read_text())
    hits, programs, node, done = [], 0, "", set()
    for s in seen:
        if s.startswith("@@"):
            node = s[2:]
            continue
        if s in done:
            continue
        done.add(s)
        programs += 1
        hits += [{"source": f"test program first parsed in {node}", **h} for h in unknown_keys(s)]
    return hits, programs, summary


def main(argv: list[str]) -> int:
    tracked = [
        ROOT / p
        for p in subprocess.run(
            ["git", "ls-files", "*.bio"], cwd=ROOT, capture_output=True, text=True
        ).stdout.split()
    ]
    generated = sorted((ROOT / "data/knowledge/compiled").glob("*.bio"))
    result: dict = {"accepted_from": "genomeos.lang.grammar.BLOCKS + COMMON"}
    h, n = census_files(tracked)
    result["tracked_bio"] = {"programs": n, "hits": h}
    h, n = census_markdown()
    result["markdown_fenced"] = {"programs": n, "hits": h}
    h, n = census_code()
    result["code_literals"] = {"programs": n, "hits": h}
    h, n = census_files(generated)
    result["generated_bio_untracked"] = {"programs": n, "hits": h}
    if "--tests" in argv:
        h, n, summary = census_tests(
            Path(argv[argv.index("--tests") + 1])
            if len(argv) > argv.index("--tests") + 1
            else ROOT / ".census"
        )
        result["test_programs"] = {"programs": n, "pytest": summary, "hits": h}
    keys = sorted(
        {
            (x["kind"], x["key"], x["closest"])
            for part in result.values()
            if isinstance(part, dict)
            for x in part.get("hits", [])
            if "key" in x
        }
    )
    result["distinct_unknown_keys"] = [{"kind": k, "key": key, "closest": c} for k, key, c in keys]
    result["total_hits"] = sum(
        len(part["hits"]) for part in result.values() if isinstance(part, dict) and "hits" in part
    )
    # scripts/when_census.py, this census's sibling, writes the same way
    from genomeos import manifest as mf
    from genomeos.results import save_result

    def ls(*patterns: str) -> list[Path]:
        out = subprocess.run(["git", "ls-files", *patterns], cwd=ROOT, capture_output=True, text=True)
        return [ROOT / p for p in out.stdout.split()]

    manifest = {
        "sources": [
            {
                "accession": "git ls-files *.bio, *.md fences, *.py/*.js/*.html literals",
                "version": "the checkout's tracked files at the run (code.git_sha; digests in inputs)",
            },
            {"accession": "data/knowledge/compiled/*.bio", "version": "untracked, as built in this checkout"},
        ],
        "inputs": [
            mf.input_entry("genomeos/lang/parser.py"),
            mf.input_entry("genomeos/lang/grammar.py"),
            mf.files_entry("tracked *.bio", tracked),
            mf.files_entry("tracked *.md", ls("*.md")),
            mf.files_entry("tracked *.py, *.js, *.html", ls("*.py", "*.js", "*.html")),
            mf.files_entry("data/knowledge/compiled/*.bio (untracked)", generated),
        ],
        "assembly": "n/a: program text, no genome coordinates",
        "coordinates": "n/a: no genomic intervals",
        "parameters": {"tests": "--tests" in argv},
        "exclusions": [],
        "partitions": "n/a: no evaluation split",
    }
    save_result(OUT.stem, result, manifest=manifest)
    for part, v in result.items():
        if isinstance(v, dict):
            print(part, v.get("programs"), "programs,", len(v["hits"]), "hits", v.get("pytest", ""))
    for k in result["distinct_unknown_keys"]:
        print("  ", k)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
