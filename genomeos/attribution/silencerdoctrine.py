# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The guard on the PREDICTED `silencer_like` key, and on prose about it.

`silencer_like` is a DIRECTION and never a function. `genomeos/lang/grammar.py` defines
`molecular_role: silencer` as "a repressive element; never inferred from a direction of effect
alone", and `silencer_like` is inferred from exactly a direction of effect alone: per element, on
deleting it in one deterministic run, the single largest-moving gene in its 1 Mb window moved UP.

Every guard the project had on that doctrine before this module sat on the MEASURED path
(`attribution/increase_links.py`, `attribution/hct116_count.py`, `attribution/direction_link.py`,
`attribution/measured.py`) or on the COMPILED label (the programme headers, the grammar axis). The
layer that generates the 230,839 calls had none, and neither did prose: before 7072f30 eight
published sentences gave elements the repressive FUNCTION, one of them in the predicted module's
own docstring. This module is the two guards that were missing.

PART 1, over the committed results that carry the key. The key is a count of elements, taken from a
sign, and these are the properties that make it one:

  * every occurrence is a non-negative `int` -- never a rate, never a share, never a label;
  * it never exceeds the population counted beside it in the same object;
  * where that object also states `strong`, `weak` and `with_predicted_target`, those agree, so the
    population the count is a subset of is the named-target population and not the scored one;
  * where the result also carries its per-element rows, the count is RECOMPUTED from them and must
    equal both the number of rows the model calls `represses` AND the number whose recorded
    log2 fold change is POSITIVE. The second equality is the one that makes it a direction: it says
    the count counts elements whose deletion RAISES the target, in the committed bytes;
  * no object carrying the key carries a sibling that names a role or a function.

PART 2, over tracked text. A sentence may not assert that an element HAS the repressive function.
Two tiers, because the two carry different weight:

  * `FUNCTION_ASSERTION` -- the X-as-the-thing form, in which one of `behave`, `act`, `function`,
    `serve`, `work` or `operate` is followed by `as` and then the bare noun. This is the shape all
    eight overstatements took and the shape the grammar forbids outright, so it is refused everywhere
    in the corpus with no subject test at all. Exactly one line in the tree is excused, by path and
    line number, in `ALLOWED`.
  * `ROLE_ASSERTION` -- the copula form, in which a form of `be` or `become` is followed by the bare
    noun, or by `repressive`. Refused only where the enclosing block is about the predicted layer
    (`PREDICTED_SUBJECT`) and the sentence is not a REFUSAL of the inference. The project's existing
    sentences of this shape are all refusals -- the doctrine being stated, which must stay.

A sentence is assembled across line breaks inside a block, because one of the eight ran over a line
break and a line-at-a-time scan would have missed it. For a `.py` file the corpus is its comments
and all of its string literals, which is a superset of its docstrings.
"""

from __future__ import annotations

import io
import json
import re
import subprocess
import tokenize
from collections.abc import Iterator
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent.parent

#: The name of the predicted key. One spelling, in one place, so nothing here can drift from it.
KEY = "silencer_like"

#: The three result families the predicted summary writes the key into. `genomeos/predict/
#: enhancer_target.py` emits it from `summarise()`; these are the committed names that reach it.
PREDICTED_FAMILIES = ("enhancer_targets_all", "enhancer_targets", "constrained_targets")

#: Sibling keys that name a population the count must be a subset of.
POPULATION_SIBLINGS = (
    "with_predicted_target",
    "with_a_named_predicted_target",
    "named",
    "elements_scored",
    "scored",
    "elements",
    "n",
)

#: A sibling key or value that would make the key a role rather than a count. `molecular_role` is
#: the grammar axis the word `silencer` belongs to and may never appear beside this count.
ROLE_SIBLING = re.compile(r"^(silencer|silencers|is_silencer|repressor|repressors|molecular_role|role)$")

# ---- part 2: the prose patterns ------------------------------------------------------------------

#: Tier A. The X-as-the-thing form: an element is said to HAVE the function. Refused everywhere.
FUNCTION_ASSERTION = re.compile(
    r"\b(?:behav|act|function|serv|work|operat|double|perform)\w*\s+as\s+"
    r"(?:a\s+|an\s+|the\s+)?(?:silencer|repressor)s?\b",
    re.IGNORECASE,
)

#: Tier B. The copula form. Refused where the block is about the predicted layer and the sentence is
#: not a refusal. `silencer-like` and `silencer_like` are excluded by the trailing class: the hedged
#: word is the correct word and must not be flagged.
ROLE_ASSERTION = re.compile(
    r"\b(?:is|are|was|were|be|being|becomes?|became|remains?)\s+"
    r"(?:(?:a|an)\s+)?(?:(?:silencer|repressor)s?\b(?![-_\w])|repressive\b)",
    re.IGNORECASE,
)

#: A sentence that REFUSES the inference rather than making it. The doctrine is written in sentences
#: of exactly tier B's shape, denied -- the project's own "never a silencer, never a repressor" --
#: and those are the guard's own subject matter, not its target.
REFUSAL = re.compile(
    r"\b(?:never|not|n't|cannot|can't|no\s|nor\s|neither|without|forbid\w*|refus\w*|withheld|"
    r"withhold\w*|unsupported|denies|deny|rather\s+than|instead\s+of)\b",
    re.IGNORECASE,
)

#: The subject test for tier B: the enclosing block is about the predicted layer, the layer that
#: emits the key. A measured or published repressor is a different subject and stays allowed.
PREDICTED_SUBJECT = re.compile(
    r"silencer_like|silencer-like|enhancer_targets|constrained_targets|enhancer-like|"
    r"AlphaGenome|cCRE|EH38E|\bdelet\w+|\bpredict\w+|1\s*Mb\s+window",
    re.IGNORECASE,
)

#: Lines excused from tier A, by path and line number, never by a pattern.
#:
#: genomeos/predict/enhancer_target.py line 11 is the SOURCE sentence that the five corrected
#: documents inherited, and it sits inside the frozen paid-study closure: it is on the computed
#: import closure of the AstroREG study (`counting_path` in data/results/astroreg2_registration.json,
#: asserted by this module's tests), so editing it -- even uncommitted -- moves a digest a paid run
#: is waiting on. docs/ROADMAP.md names it as one of the three one-target consumers inside that
#: closure. A ROADMAP row queues the correction for after the closure thaws; until then the sentence
#: is wrong and excused, not right.
ALLOWED: tuple[tuple[str, int], ...] = (
    ("genomeos/predict/enhancer_target.py", 11),  # frozen closure; correction queued, see above
)

#: The result that holds the closure this allowlist's reason appeals to.
CLOSURE_REGISTRATION = "data/results/astroreg2_registration.json"

#: The corpus. Tracked paths only, so the guard reads the tree git holds and not a lane's scratch.
CORPUS = (
    "docs/*.md",
    "README.md",
    "CONTRIBUTING.md",
    "genomeos/**/*.py",
    "genomeos/web/static/*",
    "data/demo/*.bio",
)

# ---- part 1 --------------------------------------------------------------------------------------


def _tracked(root: Path, patterns: tuple[str, ...]) -> list[str]:
    out = subprocess.run(
        ["git", "-C", str(root), "ls-files", "--", *patterns],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.split("\n")
    return [p for p in out if p]


def results_carrying_the_key(root: Path = ROOT) -> dict[str, Any]:
    """Every tracked result under data/results whose JSON contains the key, and where."""
    found: dict[str, Any] = {}
    for rel in _tracked(root, ("data/results/*.json",)):
        text = (root / rel).read_text()
        if f'"{KEY}"' not in text:
            continue
        found[rel] = json.loads(text)
    return found


def occurrences(doc: Any, pointer: str = "") -> Iterator[tuple[str, Any, dict[str, Any]]]:
    """Each (pointer, value, enclosing object) at which the key appears."""
    if isinstance(doc, dict):
        for k, v in doc.items():
            if k == KEY:
                yield f"{pointer}/{k}", v, doc
            yield from occurrences(v, f"{pointer}/{k}")
    elif isinstance(doc, list):
        for i, v in enumerate(doc):
            yield from occurrences(v, f"{pointer}[{i}]")


def family(rel: str) -> str | None:
    """The predicted-summary family a result belongs to, or None if it is not one of them."""
    stem = Path(rel).stem
    for fam in PREDICTED_FAMILIES:
        if stem.startswith(fam):
            return fam
    return None


def key_violations(rel: str, doc: Any) -> list[str]:
    """Every way one committed result breaks the count-not-a-function invariants."""
    bad: list[str] = []
    for pointer, value, parent in occurrences(doc):
        where = f"{rel}{pointer}"
        if isinstance(value, bool) or not isinstance(value, int):
            bad.append(f"{where}: {KEY} is {type(value).__name__} {value!r}, not an int count")
            continue
        if value < 0:
            bad.append(f"{where}: {KEY} is negative ({value})")
        for sib in POPULATION_SIBLINGS:
            cap = parent.get(sib)
            if isinstance(cap, int) and not isinstance(cap, bool) and value > cap:
                bad.append(f"{where}: {KEY} {value} exceeds its population {sib}={cap}")
        if {"strong", "weak", "with_predicted_target"} <= set(parent) and (
            parent["strong"] + parent["weak"] != parent["with_predicted_target"]
        ):
            bad.append(
                f"{where}: strong+weak={parent['strong'] + parent['weak']} but "
                f"with_predicted_target={parent['with_predicted_target']}, so the population "
                f"{KEY} is a subset of is not the named-target population"
            )
        for sib in parent:
            if sib != KEY and ROLE_SIBLING.match(sib):
                bad.append(f"{where}: {KEY} sits beside {sib!r}, which names a role, not a count")
    bad.extend(recount_violations(rel, doc))
    return bad


def recount_violations(rel: str, doc: Any) -> list[str]:
    """Where the result carries its element rows, recompute the count from them.

    Both equalities are asserted. Against `action == "represses"` it says the summary counts what the
    model's own direction field says. Against `log2_fold_change > 0` it says that field means the
    target ROSE, which is the whole content of calling the key a direction.
    """
    if not isinstance(doc, dict):
        return []
    rows = doc.get("elements")
    summary = doc.get("summary")
    if not isinstance(rows, list) or not isinstance(summary, dict) or KEY not in summary:
        return []
    preds = [r["predicted"] for r in rows if isinstance(r, dict) and r.get("predicted")]
    if not preds:
        return []
    by_action = sum(1 for p in preds if p.get("action") == "represses")
    by_sign = sum(
        1 for p in preds if isinstance(p.get("log2_fold_change"), (int, float)) and p["log2_fold_change"] > 0
    )
    bad = []
    if summary[KEY] != by_action:
        bad.append(f"{rel}: summary {KEY}={summary[KEY]} but {by_action} rows say represses")
    if summary[KEY] != by_sign:
        bad.append(
            f"{rel}: summary {KEY}={summary[KEY]} but {by_sign} rows record a POSITIVE "
            f"log2_fold_change, so the count is not a count of elements whose deletion raises "
            f"the target"
        )
    return bad


def key_report(root: Path = ROOT) -> dict[str, Any]:
    """What the key guard sees across the committed tree."""
    docs = results_carrying_the_key(root)
    fams: dict[str, list[str]] = {}
    for rel in sorted(docs):
        fams.setdefault(family(rel) or "other", []).append(rel)
    recounted = [
        rel
        for rel, doc in docs.items()
        if isinstance(doc, dict)
        and isinstance(doc.get("elements"), list)
        and isinstance(doc.get("summary"), dict)
        and KEY in doc["summary"]
        and any(isinstance(r, dict) and r.get("predicted") for r in doc["elements"])
    ]
    violations: list[str] = []
    for rel, doc in sorted(docs.items()):
        violations.extend(key_violations(rel, doc))
    return {
        "results_carrying_the_key": len(docs),
        "by_family": {k: len(v) for k, v in sorted(fams.items())},
        "occurrences": sum(1 for doc in docs.values() for _ in occurrences(doc)),
        "recounted_from_their_own_element_rows": len(recounted),
        "violations": violations,
    }


# ---- part 2 --------------------------------------------------------------------------------------


def blocks(path: Path) -> Iterator[tuple[int, str]]:
    """(first line number, text) for each prose block of a corpus file.

    For a `.py` file the blocks are its comment runs and its string literals, read with `tokenize`,
    so code is never scanned. For anything else a block is a run of non-blank lines, which keeps a
    markdown paragraph and a table row whole and lets a sentence cross a line break.
    """
    text = path.read_text(errors="replace")
    if path.suffix == ".py":
        yield from _python_blocks(text)
        return
    start = None
    buf: list[str] = []
    for n, line in enumerate(text.splitlines(), 1):
        if line.strip():
            if start is None:
                start = n
            buf.append(line)
        elif buf:
            yield start or 1, "\n".join(buf)
            start, buf = None, []
    if buf:
        yield start or 1, "\n".join(buf)


def _python_blocks(text: str) -> Iterator[tuple[int, str]]:
    try:
        toks = list(tokenize.generate_tokens(io.StringIO(text).readline))
    except (tokenize.TokenError, IndentationError, SyntaxError):
        return
    run: list[str] = []
    run_start = 0
    for tok in toks:
        if tok.type == tokenize.COMMENT:
            if not run:
                run_start = tok.start[0]
            run.append(tok.string.lstrip("# "))
            continue
        if run:
            yield run_start, "\n".join(run)
            run = []
        if tok.type == tokenize.STRING:
            yield tok.start[0], tok.string
    if run:
        yield run_start, "\n".join(run)


_SENTENCE_END = re.compile(r"(?<=[.!?;])\s+(?=[^\s])")


def sentences(block: str) -> Iterator[tuple[int, str]]:
    """(line offset inside the block, sentence). A decimal point does not end a sentence."""
    flat = re.sub(r"(\d)\.(\d)", r"\1․\2", block)
    pos = 0
    for part in _SENTENCE_END.split(flat):
        offset = flat.count("\n", 0, flat.index(part, pos))
        pos = flat.index(part, pos) + len(part)
        yield offset, part.replace("․", ".")


def scan_text(rel: str, path: Path, allow: tuple[tuple[str, int], ...] = ALLOWED) -> list[dict[str, Any]]:
    """Every offending sentence in one corpus file, with the tier that caught it.

    The line reported is the line the forbidden phrase itself sits on, not the line its sentence
    starts on, so an allowlist entry names the offending line and not a line above it.
    """
    out: list[dict[str, Any]] = []
    for first, block in blocks(path):
        low = block.lower()
        if "silencer" not in low and "repress" not in low:
            continue
        block_is_predicted = bool(PREDICTED_SUBJECT.search(block))
        for offset, sentence in sentences(block):
            for tier, pattern in (
                ("FUNCTION_ASSERTION", FUNCTION_ASSERTION),
                ("ROLE_ASSERTION", ROLE_ASSERTION),
            ):
                if tier == "ROLE_ASSERTION" and (not block_is_predicted or REFUSAL.search(sentence)):
                    continue
                m = pattern.search(sentence)
                if not m:
                    continue
                line = first + offset + sentence[: m.start()].count("\n")
                if (rel, line) in allow:
                    break
                out.append(
                    {
                        "path": rel,
                        "line": line,
                        "tier": tier,
                        "phrase": m.group(0),
                        "sentence": " ".join(sentence.split())[:300],
                    }
                )
                break
    return out


def scan_tree(root: Path = ROOT, allow: tuple[tuple[str, int], ...] = ALLOWED) -> list[dict[str, Any]]:
    """Every offending sentence in the tracked corpus."""
    out: list[dict[str, Any]] = []
    for rel in _tracked(root, CORPUS):
        p = root / rel
        if p.is_file():
            out.extend(scan_text(rel, p, allow))
    return out


def allowlist_still_earns_its_place(root: Path = ROOT) -> list[str]:
    """An allowlist entry that no longer excuses an offence is a hole, and is reported as one.

    The check is the scan run with the allowlist OFF: each entry must appear in it, as a
    FUNCTION_ASSERTION, at exactly the path and line the entry names.
    """
    unfiltered = {(h["path"], h["line"]): h for h in scan_tree(root, allow=())}
    bad = []
    for rel, line in ALLOWED:
        hit = unfiltered.get((rel, line))
        if hit is None:
            bad.append(f"{rel}:{line}: allowlisted but the scan finds no offence there to excuse")
        elif hit["tier"] != "FUNCTION_ASSERTION":
            bad.append(f"{rel}:{line}: allowlisted but the offence is {hit['tier']}")
    return bad
