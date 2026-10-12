#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""globscope: the scan that applies the committed rule. WRITTEN AND COMMITTED BEFORE IT RUNS.

    uv run python scripts/globscope_classify.py

The rule is `genomeos/attribution/globscope.py`, committed at `d31dd3d` before any call site was read.
This script gathers the FACTS the rule's tests need and hands them to `globscope.classify_site` and
`globscope.classify_variation`. It decides nothing itself: every class in the result comes out of those
two functions, so the classification cannot be bent here.

HOW THE FACTS ARE GATHERED, declared before any site is opened:

1. AST ONLY, never a grep. The grep undercount is the recorded reason the inherited number is 18 and
   not 17 (`scripts/trio_q100_phase.py:72` spans four lines with its filter on line 74).

2. `refusal_in_enclosing_scope`: the nearest enclosing `FunctionDef` is walked for an `Assert`, a
   `Raise` under an `If`, a `sys.exit`/`exit` call, or a `must_be_present` call, and the names in its
   governing test are intersected with the FEEDERS of the glob -- the names appearing in the
   comprehension's source iterable, plus the variable the filtered list is bound to.

3. `glob_root_resolved`, `all_matched_files_tracked`: literal path strings are resolved out of the
   paths argument by following module-level constants and same-module function definitions to a depth
   of 3, including f-string literal parts. A directory is then taken from each literal. TWO NARROWINGS
   OF THE RULE, both stated here before the scan and both STRICTER than the rule requires, never looser:
   (a) tracked-ness is measured over the files ON DISK UNDER THE RESOLVED DIRECTORY rather than over the
   glob pattern's own matches, because the pattern is often not statically resolvable; a wholly tracked
   directory proves CANNOT_SHRINK for any subset of it, and a wholly untracked one proves the opposite,
   (b) a MIXED directory -- some files tracked, some not -- is recorded as `None` and the site becomes
   UNDECIDABLE_FROM_THE_CODE, because the pattern might match only the tracked files and this scan
   cannot tell. Neither narrowing can move a site INTO CAN_SHRINK_SILENTLY.

4. `independent_scope_recorded`: searched over the WHOLE MODULE, not only the function holding the
   glob, because the manifest is conventionally written in a separate `manifest()` function in the same
   file and the rule's subject is "the result's own artefact states the expected scope". This is a
   WIDENING of where the scan looks, and it can only move a site OUT of CAN_SHRINK_SILENTLY and into
   SCOPE_IS_DECLARED_ELSEWHERE -- that is, it works against this lane's own finding, which is the
   direction a widening has to work in to be allowed.

5. FACTS THE MACHINE CANNOT RESOLVE may be supplied by READING the site, from the optional file
   `data/results/globscope_read_facts.json`, and only under three conditions enforced by
   `_apply_read_facts`: the mechanical scan must have marked that fact unresolved, every supplied fact
   must carry `evidence` naming a `path:line`, and no supplied fact may name a class. The class is
   still computed by `globscope.classify_site` from the facts. This procedure is declared HERE, before
   any site is read, precisely so that reading a site can add evidence and cannot add a verdict.

Input paths go to `input_entry` REPO-RELATIVE, for the reason `globcensus_count.py` records.
"""

from __future__ import annotations

import ast
import json
import subprocess
from pathlib import Path
from typing import Any

from genomeos import manifest as mf
from genomeos.attribution import globcensus as gc
from genomeos.attribution import globscope as gs
from genomeos.results import save_result

ROOT = Path(__file__).resolve().parent.parent
RESULT = "globscope"
OWN_CODE = (
    "genomeos/attribution/globscope.py",
    "scripts/globscope_classify.py",
    "tests/test_globscope.py",
)
READ_FACTS = ROOT / "data/results/globscope_read_facts.json"

#: Top-level directories a literal path may start with for this scan to treat it as a repo path.
REPO_TOPS = ("data", "genomeos", "scripts", "tests", "docs", "web")

#: Dict keys and manifest helpers that count as an independent record of expected scope.
SCOPE_KEYS = (
    "chromosomes",
    "chroms",
    "chromosome",
    "chromosome_list",
    "n_chromosomes",
    "chromosomes_scanned",
    "samples",
    "n_samples",
    "sample_list",
    "expected_files",
    "expected",
    "scope",
    "panel",
)

#: The six varying groups, inherited from globcensus and not re-derived.
VARYING_SOURCE = "data/results/globcensus.json : groups_whose_count_varies"

#: lane-shrinkimpact's settled fourth condition for the two groups it recomputed (8d620bc, 9a68bf2).
#: Carried, not re-measured: "realised in 0 of 6".
SHRINKIMPACT_SETTLED = {
    "per_element_response_cache": {
        "missing_members_carry_eligible_data": False,
        "evidence": "globcensus.STAGE_2_REFUTED: in the four results of 2026-10-02, chr13 and chrY"
        " carry ZERO ELIGIBLE LINKS, so 22 IS the benchmark's own scope rather than a shortfall of it",
    },
    "GENCODE 50 per-chromosome annotation, for the nearest protein-coding TSS": {
        "missing_members_carry_eligible_data": False,
        "evidence": "globcensus.STAGE_2_REFUTED: closure_cell_inventory's chr21 WAS READ and was filed"
        " under another group's label -- a prefix-grouping artefact, not an absence; chrY is not in the"
        " CRISPRi library at all",
    },
}


# ------------------------------------------------------------------------------- AST fact gathering


def _parents(tree: ast.AST) -> dict[int, ast.AST]:
    out: dict[int, ast.AST] = {}
    for node in ast.walk(tree):
        for child in ast.iter_child_nodes(node):
            out[id(child)] = node
    return out


def _enclosing_function(node: ast.AST, parents: dict[int, ast.AST]) -> ast.AST | None:
    cur = parents.get(id(node))
    while cur is not None:
        if isinstance(cur, (ast.FunctionDef, ast.AsyncFunctionDef)):
            return cur
        cur = parents.get(id(cur))
    return None


def _names(node: ast.AST) -> set[str]:
    return {n.id for n in ast.walk(node) if isinstance(n, ast.Name)} | {
        n.attr for n in ast.walk(node) if isinstance(n, ast.Attribute)
    }


def _paths_argument(call: ast.Call) -> ast.AST | None:
    if len(call.args) >= 2:
        return call.args[1]
    for kw in call.keywords:
        if kw.arg == "paths":
            return kw.value
    return None


def _source_iterable(expr: ast.AST) -> ast.AST:
    """The pre-filter iterable of a comprehension, or the expression itself."""
    if isinstance(expr, (ast.ListComp, ast.GeneratorExp, ast.SetComp)) and expr.generators:
        return expr.generators[0].iter
    return expr


def _module_level(tree: ast.Module) -> dict[str, list[ast.AST]]:
    out: dict[str, list[ast.AST]] = {}
    for stmt in tree.body:
        if isinstance(stmt, ast.Assign):
            for t in stmt.targets:
                if isinstance(t, ast.Name):
                    out.setdefault(t.id, []).append(stmt.value)
        elif isinstance(stmt, ast.AnnAssign) and isinstance(stmt.target, ast.Name) and stmt.value:
            out.setdefault(stmt.target.id, []).append(stmt.value)
        elif isinstance(stmt, (ast.FunctionDef, ast.AsyncFunctionDef)):
            out.setdefault(stmt.name, []).append(stmt)
    return out


def _literal_strings(node: ast.AST, defs: dict[str, list[ast.AST]], depth: int = 3) -> set[str]:
    """Every string literal reachable from `node` through module-level names, to `depth` hops."""
    found: set[str] = set()
    if depth < 0:
        return found
    for sub in ast.walk(node):
        if isinstance(sub, ast.Constant) and isinstance(sub.value, str):
            found.add(sub.value)
        elif isinstance(sub, ast.Name) and sub.id in defs:
            for d in defs[sub.id]:
                found |= _literal_strings(d, defs, depth - 1)
        elif isinstance(sub, ast.Attribute) and sub.attr in defs:
            for d in defs[sub.attr]:
                found |= _literal_strings(d, defs, depth - 1)
    return found


def _directory_of(literal: str) -> str | None:
    """The repo-relative directory a path literal names, or None if it is not a repo path."""
    s = literal.strip().lstrip("./")
    if "/" not in s:
        return None
    if s.split("/")[0] not in REPO_TOPS:
        return None
    if "{" in s:
        s = s.split("{")[0]
        return s.rsplit("/", 1)[0] if "/" in s else None
    last = s.rsplit("/", 1)[-1]
    if "." in last or "*" in last:
        return s.rsplit("/", 1)[0]
    return s


# ------------------------------------------------------------------------------------- git facts


def _git(args: list[str]) -> str:
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=False).stdout


def _directory_tracked(directory: str) -> dict[str, Any]:
    """On-disk files under `directory`, and how many of them git tracks."""
    d = ROOT / directory
    if not d.is_dir():
        return {"directory": directory, "exists": False, "on_disk": 0, "tracked": 0, "verdict": None}
    on_disk = {str(p.relative_to(ROOT)) for p in d.rglob("*") if p.is_file() and ".git" not in p.parts}
    tracked = set(_git(["ls-files", "--", directory]).split("\n")) - {""}
    n_tracked = len(on_disk & tracked)
    if not on_disk:
        verdict = None
    elif n_tracked == len(on_disk):
        verdict = True
    elif n_tracked == 0:
        verdict = False
    else:
        verdict = None  # mixed: narrowing (b) of the docstring
    return {
        "directory": directory,
        "exists": True,
        "on_disk": len(on_disk),
        "tracked": n_tracked,
        "verdict": verdict,
    }


# ------------------------------------------------------------------------------------- per site


def site_facts(path: str, line: int) -> dict[str, Any]:
    """Every fact the rule's tests need, for the `files_entry` call at `path:line`."""
    tree = ast.parse((ROOT / path).read_text())
    parents = _parents(tree)
    defs = _module_level(tree)
    call = None
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and node.lineno == line:
            name = getattr(node.func, "attr", None) or getattr(node.func, "id", None)
            if name == "files_entry":
                call = node
                break
    if call is None:
        return {"path": path, "line": line, "error": "no files_entry call at that line"}

    label = call.args[0].value if call.args and isinstance(call.args[0], ast.Constant) else None
    arg = _paths_argument(call)
    src_iter = _source_iterable(arg) if arg is not None else None
    feeders = _names(src_iter) if src_iter is not None else set()

    # the variable the filtered list is bound to, if the call takes a bare Name
    if isinstance(arg, ast.Name):
        feeders |= {arg.id}
        for d in defs.get(arg.id, []):
            feeders |= _names(_source_iterable(d))

    fn = _enclosing_function(call, parents)
    refusal_at: list[str] = []
    scope_body = fn if fn is not None else tree
    fn_parents = _parents(scope_body)
    for node in ast.walk(scope_body):
        test = None
        if isinstance(node, ast.Assert):
            test = node.test
        elif isinstance(node, ast.Raise):
            cur = fn_parents.get(id(node))
            while cur is not None:
                if isinstance(cur, ast.If):
                    test = cur.test
                    break
                cur = fn_parents.get(id(cur))
        elif isinstance(node, ast.Call):
            fname = getattr(node.func, "attr", None) or getattr(node.func, "id", None)
            if fname in ("exit", "must_be_present"):
                test = node
        if test is not None and _names(test) & feeders:
            refusal_at.append(f"{path}:{getattr(node, 'lineno', line)}")

    # independent scope, over the whole module (widening 4 of the docstring)
    scope_at: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Dict):
            for k, v in zip(node.keys, node.values, strict=False):
                if isinstance(k, ast.Constant) and k.value in SCOPE_KEYS and _names(v) & feeders:
                    scope_at.append(f"{path}:{k.lineno} key={k.value!r}")
        elif isinstance(node, ast.keyword) and node.arg in SCOPE_KEYS and _names(node.value) & feeders:
            scope_at.append(f"{path}:{getattr(node.value, 'lineno', line)} kw={node.arg}")

    literals = _literal_strings(arg, defs) if arg is not None else set()
    dirs = sorted({d for lit in literals if (d := _directory_of(lit))})
    dir_facts = [_directory_tracked(d) for d in dirs]
    verdicts = {f["verdict"] for f in dir_facts if f["exists"]}
    if not dir_facts or not any(f["exists"] for f in dir_facts):
        tracked: bool | None = None
        on_disk = 0
    elif verdicts == {True}:
        tracked, on_disk = True, sum(f["on_disk"] for f in dir_facts)
    elif verdicts == {False}:
        tracked, on_disk = False, sum(f["on_disk"] for f in dir_facts)
    else:
        tracked, on_disk = None, sum(f["on_disk"] for f in dir_facts)

    return {
        "path": path,
        "line": line,
        "group_label": label,
        "enclosing_function": getattr(fn, "name", "<module>"),
        "feeders": sorted(feeders),
        "refusal_in_enclosing_scope": bool(refusal_at),
        "refusal_at": refusal_at,
        "glob_root_resolved": bool(dirs),
        "resolved_directories": dir_facts,
        "matched_files_on_disk": on_disk,
        "all_matched_files_tracked": tracked,
        "independent_scope_recorded": bool(scope_at),
        "independent_scope_at": scope_at,
        "unresolved": sorted(
            k
            for k, v in {
                "glob_root_resolved": bool(dirs),
                "all_matched_files_tracked": tracked is not None,
            }.items()
            if not v
        ),
    }


def _apply_read_facts(facts: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Facts supplied by READING a site, under the three conditions declared in the docstring."""
    if not READ_FACTS.exists():
        return []
    supplied = json.loads(READ_FACTS.read_text())
    applied: list[dict[str, Any]] = []
    index = {(f["path"], f["line"]): f for f in facts}
    for entry in supplied["facts"]:
        key = (entry["path"], entry["line"])
        if key not in index:
            raise SystemExit(f"read fact names no scanned site: {key}")
        if "class" in entry or "classification" in entry:
            raise SystemExit(f"a read fact may not name a class: {key}")
        if ":" not in str(entry.get("evidence", "")):
            raise SystemExit(f"read fact carries no path:line evidence: {key}")
        site = index[key]
        for field, value in entry["sets"].items():
            if field not in site:
                raise SystemExit(f"read fact names no such field: {field}")
            mechanical_unresolved = field in site["unresolved"] or (
                field == "independent_scope_recorded" and not site["independent_scope_recorded"]
            )
            if not mechanical_unresolved:
                raise SystemExit(
                    f"read fact overrides a fact the scan RESOLVED, which is not allowed: {key} {field}"
                )
            site[field] = value
        site.setdefault("read_evidence", []).append(entry["evidence"])
        applied.append({"site": f"{key[0]}:{key[1]}", "sets": entry["sets"], "evidence": entry["evidence"]})
    return applied


# ------------------------------------------------------------------------- the six varying groups


def _scope_from_result(result_stem: str) -> dict[str, Any]:
    f = ROOT / "data/results" / f"{result_stem}.json"
    try:
        d = json.loads(f.read_text())
    except Exception:
        return {"result": result_stem, "readable": False}
    params = (d.get("result_manifest") or {}).get("parameters") or {}
    scope_count, scope_field = None, None
    for key in SCOPE_KEYS:
        if key in params:
            v = params[key]
            if isinstance(v, (list, tuple)):
                scope_count, scope_field = len(v), f"parameters.{key} (list of {len(v)})"
            elif isinstance(v, int) and not isinstance(v, bool):
                scope_count, scope_field = v, f"parameters.{key} = {v}"
            else:
                scope_field = f"parameters.{key} is not countable: {str(v)[:60]!r}"
            if scope_count is not None:
                break
    excl = (d.get("result_manifest") or {}).get("exclusions") or []
    return {
        "readable": True,
        "scope_count": scope_count,
        "scope_field": scope_field,
        "exclusions_mention_a_member": any(
            "chr" in str(e).lower() or "sample" in str(e).lower() for e in excl
        ),
    }


def varying_group_rows(groups: dict[str, list[tuple[str, Any]]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for group, decls in gc.varying(groups).items():
        settled = SHRINKIMPACT_SETTLED.get(group)
        per_result = []
        for result_stem, count in sorted(groups[group]):
            scope = _scope_from_result(result_stem)
            fourth = settled["missing_members_carry_eligible_data"] if settled else None
            cls = gs.classify_variation(
                artefact_states_a_scope=scope.get("scope_count") is not None,
                declared_count=count,
                scope_count=scope.get("scope_count"),
                missing_members_present_in_store_today=None,
                missing_members_named_in_exclusions=scope.get("exclusions_mention_a_member"),
                missing_members_carry_eligible_data=fourth,
            )
            per_result.append({"result": result_stem, "declared": count, "class": cls, **scope})
        classes = {r["class"] for r in per_result}
        rows.append(
            {
                "group": group,
                "counts_seen": decls,
                "declarations": len(per_result),
                "per_result": per_result,
                "group_class": (
                    classes.pop()
                    if len(classes) == 1
                    else "MIXED_ACROSS_RESULTS: " + ", ".join(sorted(classes))
                ),
                "fourth_condition": (
                    settled["evidence"]
                    if settled
                    else "NOT ESTABLISHED for this group: no recomputation on the complete store"
                    " exists, so a shrink verdict is unavailable"
                ),
            }
        )
    return rows


# ------------------------------------------------------------------------------------- manifest


def manifest(sources: list[Path]) -> dict:
    inputs = [mf.input_entry(p, partition=None) for p in sources]
    return {
        "sources": [
            {"accession": "data/results/globcensus.json, the inherited 18", "version": "as committed"},
            {"accession": "the genomeos/ and scripts/ trees, read by AST", "version": "as committed"},
            {"accession": "git ls-files, for tracked-ness of each resolved directory", "version": "HEAD"},
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "rule_committed_at": "d31dd3d, before any call site was read",
            "classes": list(gs.CLASSES),
            "decision_order": list(gs.DECISION_ORDER),
            "method": "ast.parse only, never a grep; the class comes from globscope.classify_site",
            "narrowings": "tracked-ness over the resolved directory not the pattern; a MIXED directory"
            " is UNDECIDABLE; independent scope searched over the whole module",
        },
        "exclusions": [
            "static only: nothing was run on any store, no result's bytes were touched, and nothing"
            " was repaired",
            "globcensus.py, tests/test_globcensus.py and data/results/globcensus.json are not touched",
        ],
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
        "partitions": {
            "none": "one population, the 18 inherited sites; no partition is taken and the absence of"
            " one is not a claim of independence"
        },
    }


def main() -> int:
    inherited = gs.inherited_site_count()
    calls = gc.files_entry_calls()
    filtered = [c for c in calls if c["exists_filtered"]]
    if len(filtered) != inherited:
        raise SystemExit(
            f"the scan sees {len(filtered)} exists()-filtered sites, globcensus.json pins {inherited};"
            " the inherited object and this scan must agree before anything is classified"
        )

    facts = [site_facts(c["path"], c["line"]) for c in filtered]
    applied = _apply_read_facts(facts)
    for f in facts:
        f["class"] = gs.classify_site(
            refusal_in_enclosing_scope=f["refusal_in_enclosing_scope"],
            glob_root_resolved=f["glob_root_resolved"],
            matched_files_on_disk=f["matched_files_on_disk"],
            all_matched_files_tracked=f["all_matched_files_tracked"],
            independent_scope_recorded=f["independent_scope_recorded"],
        )
    counts = {c: sum(1 for f in facts if f["class"] == c) for c in gs.DECISION_ORDER}
    if sum(counts.values()) != inherited:
        raise SystemExit("the classes do not sum to the inherited count")

    groups = gc.declared_groups()
    rows = varying_group_rows(groups)

    payload = {
        "question": (
            "which of the 18 exists()-filtered input-group call sites CAN silently shrink and which"
            " cannot, and for the 6 varying groups whether each variation is legitimate scope or a"
            " silent shrink"
        ),
        "rule": {
            "committed_at": "d31dd3d, before any call site was read",
            "classes": gs.CLASSES,
            "decision_order": list(gs.DECISION_ORDER),
            "why_undecidable_is_second": gs.WHY_UNDECIDABLE_IS_SECOND,
            "variation_classes": gs.VARIATION_CLASSES,
            "prewritten_outcomes": gs.PREWRITTEN_OUTCOMES,
        },
        "inherited": {
            "site_count": inherited,
            "source": gs.INHERITED_COUNT_SOURCE,
            "reading_carried_verbatim": gs.INHERITED_READING,
            "limit_carried_verbatim": gs.INHERITED_LIMIT,
            "stage_2_refuted": gc.STAGE_2_REFUTED,
            "stage_2_summary": gc.STAGE_2_SUMMARY,
        },
        "class_counts": counts,
        "sites": facts,
        "read_facts_applied": applied,
        "varying_groups": rows,
        "varying_source": VARYING_SOURCE,
        "does_not_claim": gs.DOES_NOT_CLAIM,
        "negative_boundary": gs.NEGATIVE_BOUNDARY,
        "stop_list": gs.STOP_LIST_READ,
        "model_requests": 0,
        "money_spent": "$0.00",
    }
    sources = [Path("data/results/globcensus.json")]
    path = save_result(RESULT, payload, manifest(sources))
    print(f"wrote {path}")
    print(json.dumps(counts, indent=1))
    for r in rows:
        print(f"  {r['group_class']:<28} {r['group'][:70]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
