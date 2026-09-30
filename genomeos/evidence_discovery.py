"""The discovery view: set-valued CRISPRi evidence from A's discovery-only review, read-only.

A's discovery-only review (registration a9ab0db, result 66051ee) asked which screened CRISPRi
observations audit A's policy would make discoverable as a set-valued attachment beside C4's current
rule. Its answer is a committed result file, `data/results/discovery_review.json`. This module reads
that file once per version, checks it against itself, and arranges it for the Evidence tab:

- `records`: one record per observation, keyed by the observation id `split|chrom:start-end|gene|cell`,
  each holding the file's fields unchanged plus the tested interval parsed from the id, a context per
  attached claim, one mutually exclusive row context, and the display strings the page shows;
- section (a), discovery: the Gasperini2019/K562 observations whose attached links differ between the
  two rules, counted against the file's own stratum count, with the file's denominator beside it;
- section (b), resolution change: the unique assignments that become unresolved under the new rule;
- section (c), other strata: labelled counts only, the registered screen readings unchanged.

Sections refer to records by id, so an observation in two sections is one record. Nothing is scored,
matched or regenerated here: every number is the file's, and the view refuses to serve when the file
disagrees with itself. The wording is fixed by the approved specification and kept verbatim.

    from genomeos.evidence_discovery import view
    view(Path(".")) -> {"source": {...}, "records": {...}, "sections": {...}, ...}
"""

from __future__ import annotations

import copy
import hashlib
import json
import re
from pathlib import Path
from typing import Any

SOURCE = Path("data/results/discovery_review.json")
# Section (a)'s study and cell and the change classes it reads, named by the approved specification.
STRATUM = "Gasperini2019|K562"
DISCOVERY_CHANGES = ("newly_attached", "attached_links_added")
SPARSE_MAX = 16  # the fixed label below states it; the view refuses a passing stratum above it

# Fixed wording, verbatim from the approved specification.
HEADER = "The discovery view adds access to evidence, not new validation."
CONTEXT_NOTE = (
    "The stated-cell subset is selected using model predictions and is not representative of all "
    "screened pairs."
)
VERDICT_NOTE = (
    "Historical C4 reading (any-cell attachment, reciprocal-0.5 rule, b7e4bf0); "
    "not a v3 stated-context verdict."
)
VERDICT_RULE = "C4 any-cell attachment, reciprocal-0.5 rule (b7e4bf0)"
SCREEN_NOTE = "passing the screen is not sufficient evidence for implementation"
COVERAGE = (
    "{n} admitted registry candidates; {pct}% of the tested interval lies inside them; "
    "the rest is unresolved."
)
SPARSE_LABEL = f"screen passed, sparse (n ≤ {SPARSE_MAX})"

NEW_RULE = "min_side_half"  # audit A's policy: overlap at least half of the smaller interval
CURRENT_RULE = "reciprocal_half"  # C4's reciprocal 0.5 rule

CLAIM_CONTEXTS = {
    "match": "a stated cell equals the assay cell",
    "other_cell": "the claim states cells, none of them the assay cell",
    "no_stated_cell": "the claim states no cell",
}
ROW_CONTEXTS = {
    "all_attached_claims_match": "all attached claims match the assay cell",
    "some_match": "some attached claims match the assay cell",
    "none_match": "no attached claim matches the assay cell",
    "no_assessable_claim_context": "attached claims exist, but none states a cell",
    "no_attached_claims": "no attached claims",
}
# The file's own per-observation reading under the new rule, which each row context must refine.
_FILE_AGREEMENT = {
    "all_attached_claims_match": "stated_cell_match",
    "some_match": "stated_cell_match",
    "none_match": "other_cell_only",
    "no_assessable_claim_context": "no_stated_cell",
    "no_attached_claims": "not_attached",
}
# The screen's registered reading for each list of decision_1_reading.
_SCREEN_OF_LIST = {
    "strata_passing": "passes",
    "strata_failing": "fails",
    "strata_undefined_shifted_zero": "undefined_shifted_count_zero",
    "strata_nothing_discovered": "nothing_discovered",
}
_ID = re.compile(
    r"^(?P<split>[^|]+)\|(?P<chrom>[^:|]+):(?P<start>\d+)-(?P<end>\d+)\|(?P<gene>[^|]+)\|(?P<cell>[^|]+)$"
)


class MissingError(FileNotFoundError):
    """The result file is not in this checkout."""


class RefusedError(ValueError):
    """The result file disagrees with itself, so the view serves nothing rather than a wrong row."""


def _words(token: str) -> str:
    return str(token).replace("_", " ")


def _pct(fraction: float) -> str:
    """A fraction as a percentage at the file's precision (four decimals of a fraction), no trailing 0."""
    return f"{fraction * 100:.2f}".rstrip("0").rstrip(".")


def parse_id(obs_id: str) -> dict[str, Any]:
    """`split|chrom:start-end|gene|cell` into its parts; the interval is the tested one."""
    m = _ID.match(obs_id)
    if not m:
        raise RefusedError(f"observation id {obs_id!r} is not split|chrom:start-end|gene|cell")
    return {
        "split": m["split"],
        "chrom": m["chrom"],
        "start": int(m["start"]),
        "end": int(m["end"]),
        "gene": m["gene"],
        "cell": m["cell"],
    }


def claim_context(claim: dict[str, Any], assay_cell: str) -> str:
    """`match`, `other_cell` or `no_stated_cell`: the claim's stated cells against the assay cell.

    Cells are compared the way the result compared them (correctness.norm_cell), and a cell counts as
    stated the way it did (correctness.stated), so the row context refines the file's own reading.
    """
    from genomeos.attribution.correctness import norm_cell, stated

    cells = [c for c in claim.get("stated_cells") or [] if stated(c)]
    if not cells:
        return "no_stated_cell"
    if any(norm_cell(c) == norm_cell(assay_cell) for c in cells):
        return "match"
    return "other_cell"


def row_context(record: dict[str, Any]) -> str:
    """One mutually exclusive category per observation, from its attached claims' contexts."""
    contexts = [claim_context(a, record["cell"]) for a in record.get("attached") or []]
    if not contexts:
        return "no_attached_claims"
    if all(c == "match" for c in contexts):
        return "all_attached_claims_match"
    if "match" in contexts:
        return "some_match"
    if "other_cell" in contexts:
        return "none_match"
    return "no_assessable_claim_context"


def _display(rec: dict[str, Any], iv: dict[str, Any]) -> dict[str, Any]:
    """The strings the Evidence tab shows for one observation, verbatim."""
    n = rec["multiplicity"][NEW_RULE]
    links = len(rec["attached"])
    candidates = "; ".join(
        f"{c['id']} (covered fraction {c['fraction_covered']}"
        + (", admitted under the current rule too)" if c["current"] else ")")
        for c in rec["candidates"]
    )
    claims = [
        f"{a['element']} {a['action']} {a['gene']} · stated cell{'' if len(a['stated_cells']) == 1 else 's'}"
        f": {', '.join(a['stated_cells']) or 'none'} · context: {a['context']}"
        + (" · attached under the current rule too" if a["current"] else "")
        for a in rec["attached"]
    ]
    return {
        "observation": (
            f"tested interval {iv['chrom']}:{iv['start']}-{iv['end']} ({rec['tested_width_bp']} bp) · "
            f"measured gene {rec['gene']} · assay cell {rec['cell']} · {rec['study']} · "
            f"{rec['split']} split · outcome: {_words(rec['outcome'])}"
        ),
        "coverage": COVERAGE.format(n=n, pct=_pct(rec["union_coverage"][NEW_RULE])),
        "candidates": f"admitted registry candidates: {candidates}"
        if candidates
        else "no admitted registry candidate",
        "links": f"{links} attached link{'s' if links != 1 else ''} under the new rule",
        "claims": claims,
        "context": f"row context: {ROW_CONTEXTS[rec['row_context']]} ({rec['row_context']})",
        "verdict": f"historical verdict: {_words(rec['verdict_current_rule'])} [{VERDICT_RULE}]",
        "verdict_note": VERDICT_NOTE,
        "status": f"{_words(rec['status_under_new_rule'])} under the new rule",
        "direction": (
            "direction relative to the attached claims (a description only): "
            f"{_words(rec['direction_description_policy'])}"
        ),
    }


def _record(src: dict[str, Any]) -> dict[str, Any]:
    """The served record: the file's fields unchanged, plus the parsed interval, contexts and display."""
    iv = parse_id(src["id"])
    if iv["split"] != src["split"] or iv["cell"] != src["cell"]:
        raise RefusedError(f"{src['id']}: the id's split or cell differs from the record's")
    if iv["end"] - iv["start"] != src["tested_width_bp"]:
        raise RefusedError(f"{src['id']}: the id's interval is not {src['tested_width_bp']} bp wide")
    if len(src["candidates"]) != src["multiplicity"][NEW_RULE]:
        raise RefusedError(
            f"{src['id']}: {len(src['candidates'])} candidates listed, multiplicity says otherwise"
        )
    rec = {
        "id": src["id"],
        "interval": {"chrom": iv["chrom"], "start": iv["start"], "end": iv["end"]},
        "gene": iv["gene"],
        "cell": src["cell"],
        "study": src["study"],
        "split": src["split"],
        "outcome": src["outcome"],
        "tested_width_bp": src["tested_width_bp"],
        "change": src["change"],
        "candidates": copy.deepcopy(src["candidates"]),
        "multiplicity": dict(src["multiplicity"]),
        "union_coverage": dict(src["union_coverage"]),
        "attached": [{**a, "context": claim_context(a, src["cell"])} for a in src["attached"]],
        "resolution": dict(src["resolution"]),
        "stated_cell_agreement": dict(src["stated_cell_agreement"]),
        "verdict_current_rule": src["verdict_current_rule"],
        "verdict_rule": VERDICT_RULE,
        "verdict_note": VERDICT_NOTE,
        "status_under_new_rule": src["status_under_new_rule"],
        "direction_description_policy": src["direction_description_policy"],
    }
    rec["row_context"] = row_context(rec)
    if _FILE_AGREEMENT[rec["row_context"]] != src["stated_cell_agreement"][NEW_RULE]:
        raise RefusedError(
            f"{src['id']}: row context {rec['row_context']} contradicts the file's "
            f"{src['stated_cell_agreement'][NEW_RULE]}"
        )
    rec["display"] = _display(rec, iv)
    return rec


def _counts(records: dict[str, dict[str, Any]], ids: list[str]) -> dict[str, int]:
    out = dict.fromkeys(ROW_CONTEXTS, 0)
    for i in ids:
        out[records[i]["row_context"]] += 1
    return out


def _stratum_line(name: str, s: dict[str, Any]) -> str:
    enrichment = "undefined" if s["enrichment"] is None else s["enrichment"]
    plural = "s" if s["real"] != 1 else ""
    return (
        f"{name}: n = {s['real']:,} observation{plural} whose attached links differ, of "
        f"{s['observations']:,} screened observations; enrichment {enrichment} "
        f"(registered reading: {s['screen']})"
    )


def build(data: dict[str, Any], path: str, sha256: str) -> dict[str, Any]:
    """The view of one parsed result file; raises RefusedError where the file disagrees with itself."""
    rows = data["changed_observations"]
    by_id: dict[str, dict[str, Any]] = {}
    for r in rows:
        if r["id"] in by_id:
            raise RefusedError(f"observation {r['id']} appears twice in changed_observations")
        by_id[r["id"]] = r
    screen = data["screen"]["attached_links_differ"]
    decision = data["decision_1_reading"]
    universe = data["universe"]["observations"]
    study, cell = STRATUM.split("|")

    # (a) discovery: counted against the file's own stratum count, or nothing is served
    a_ids = [
        r["id"]
        for r in rows
        if r["study"] == study and r["cell"] == cell and r["change"] in DISCOVERY_CHANGES
    ]
    stratum = screen[STRATUM]
    if len(a_ids) != stratum["real"]:
        raise RefusedError(
            f"section (a) has {len(a_ids)} rows but the file's {STRATUM} stratum counts {stratum['real']}"
        )
    if STRATUM not in decision["strata_passing"] or stratum["screen"] != "passes":
        raise RefusedError(f"{STRATUM} is not a passing stratum in the file's decision reading")

    # (b) resolution change: every id of newly_ambiguous.each, with its verdict kept
    each = data["newly_ambiguous"]["each"]
    b_ids = [e["id"] for e in each]
    if len(b_ids) != data["newly_ambiguous"]["n"] or len(set(b_ids)) != len(b_ids):
        raise RefusedError("newly_ambiguous.each does not list newly_ambiguous.n distinct observations")
    for e in each:
        src = by_id.get(e["id"])
        if src is None:
            raise RefusedError(f"newly ambiguous {e['id']} is not among the changed observations")
        if src["verdict_current_rule"] != e["historical_verdict_current_rule"]:
            raise RefusedError(f"{e['id']}: its historical verdict differs between the two lists")
        if src["status_under_new_rule"] != "unresolved" or e["status_under_new_rule"] != "unresolved":
            raise RefusedError(f"{e['id']}: not unresolved under the new rule")

    records: dict[str, dict[str, Any]] = {}
    for i in a_ids + b_ids:
        if i not in records:
            records[i] = _record(by_id[i])

    # (c) other strata: labelled counts only, the registered readings unchanged
    groups = {
        "strata_passing": SPARSE_LABEL,
        "strata_failing": "screen failed",
        "strata_undefined_shifted_zero": "screen undefined",
        "strata_nothing_discovered": "nothing discovered",
    }
    seen: set[str] = set()
    c_groups = []
    for key, label in groups.items():
        strata = []
        for name in decision[key]:
            if name in seen or name not in screen:
                raise RefusedError(f"stratum {name} is listed twice or has no screen reading")
            seen.add(name)
            s = screen[name]
            if s["screen"] != _SCREEN_OF_LIST[key]:
                raise RefusedError(f"{name} is listed under {key} but its screen reads {s['screen']}")
            if name == STRATUM:
                continue  # served as section (a)
            if key == "strata_passing" and s["real"] > SPARSE_MAX:
                raise RefusedError(
                    f"{name} passes with n = {s['real']}, above the sparse label's {SPARSE_MAX}"
                )
            strata.append(
                {
                    "stratum": name,
                    "n": s["real"],
                    "n_unit": "observations whose attached links differ",
                    "observations": s["observations"],
                    "observations_unit": "screened observations",
                    "enrichment": s["enrichment"],
                    "screen": s["screen"],
                    "display": _stratum_line(name, s),
                }
            )
        c_groups.append({"label": label, "list": key, "strata": strata})
    unlisted = sorted(set(screen) - seen - {"all_observations"})
    if unlisted:
        raise RefusedError(f"strata with a screen reading but no place in the decision reading: {unlisted}")
    pooled = screen["all_observations"]

    reg = data["registration"]
    a_links = sum(len(records[i]["attached"]) for i in a_ids)
    also = [i for i in b_ids if i in set(a_ids)]
    return {
        "view": "discovery",
        "header": HEADER,
        "source": {
            "path": path,
            "sha256": sha256,
            "result": data["result"],
            "date": data["date"],
            "status": data["status"],
            "question": data["question"],
        },
        "notes": {
            "context": CONTEXT_NOTE,
            "verdict": VERDICT_NOTE,
            "screen": SCREEN_NOTE,
            "unique": reg["wording"],
        },
        "rules": {
            "current": f"{CURRENT_RULE}: C4's reciprocal 0.5 rule",
            "new": f"{NEW_RULE}: audit A's policy, overlap at least half of the smaller interval",
        },
        "claim_contexts": CLAIM_CONTEXTS,
        "row_contexts": ROW_CONTEXTS,
        "records": records,
        "sections": {
            "discovery": {
                "title": "(a) Discovery: Gasperini2019 in K562",
                "stratum": STRATUM,
                "changes": list(DISCOVERY_CHANGES),
                "ids": a_ids,
                "rows": {
                    "count": len(a_ids),
                    "unit": "observations",
                    "checked_against_file": stratum["real"],
                },
                "denominator": {
                    "count": stratum["observations"],
                    "unit": "observations",
                    "definition": (
                        f"every screened {study} observation in {cell} in the {universe:,}-observation "
                        "universe, attached or not, whatever its outcome"
                    ),
                },
                "attached_links": {"count": a_links, "unit": "links"},
                "context_counts": {"counts": _counts(records, a_ids), "unit": "observations"},
                "screen": {
                    "observations": stratum["observations"],
                    "real": stratum["real"],
                    "shifted": stratum["shifted"],
                    "enrichment": stratum["enrichment"],
                    "screen": stratum["screen"],
                },
                "display": {
                    "rows": (
                        f"{len(a_ids):,} observations whose attached links differ between the two rules "
                        f"(change {' or '.join(DISCOVERY_CHANGES)})"
                    ),
                    "denominator": (
                        f"denominator: {stratum['observations']:,} observations, every screened {study} "
                        f"observation in {cell} in the {universe:,}-observation universe, attached or not"
                    ),
                    "links": f"{a_links:,} attached links under the new rule across these observations",
                    "screen": _stratum_line(STRATUM, stratum) + f"; {SCREEN_NOTE}",
                },
            },
            "resolution_change": {
                "title": (
                    "(b) Resolution change: unique under the current rule, unresolved under the new rule"
                ),
                "ids": b_ids,
                "rows": {"count": len(b_ids), "unit": "observations"},
                "historical_verdicts": {
                    "counts": dict(data["newly_ambiguous"]["historical_verdicts"]),
                    "unit": "observations",
                },
                "also_in_discovery": also,
                "context_counts": {"counts": _counts(records, b_ids), "unit": "observations"},
                "display": {
                    "rows": (
                        f"{len(b_ids):,} observations, across studies, each a unique-element assignment "
                        "under the current rule and unresolved under the new rule; each keeps its "
                        "historical verdict"
                    ),
                },
            },
            "other_strata": {
                "title": "(c) Other strata: counts only",
                "groups": c_groups,
                "pooled": {
                    "observations": pooled["observations"],
                    "real": pooled["real"],
                    "enrichment": pooled["enrichment"],
                    "screen": pooled["screen"],
                    "display": _stratum_line("all observations pooled", pooled),
                },
                "reading": decision["reading"],
                "meaning": decision["meaning"],
                "note": SCREEN_NOTE,
            },
        },
    }


_CACHE: dict[str, tuple[tuple[int, int], dict[str, Any]]] = {}


def view(root: Path) -> dict[str, Any]:
    """The discovery view of `root`'s committed result, read once per version of the file."""
    p = Path(root) / SOURCE
    try:
        st = p.stat()
    except FileNotFoundError:
        raise MissingError(f"{SOURCE} is not in this checkout") from None
    stamp = (st.st_mtime_ns, st.st_size)
    hit = _CACHE.get(str(p))
    if hit is None or hit[0] != stamp:
        raw = p.read_bytes()
        out = build(json.loads(raw), str(SOURCE), hashlib.sha256(raw).hexdigest())
        _CACHE[str(p)] = hit = (stamp, out)
    return copy.deepcopy(hit[1])
