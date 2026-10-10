# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Assay-validated silencers joined beside the derived direction call, cell line by cell line.

WHAT THIS ADDS. `molecular_role: silencer` is defined in `genomeos/lang/grammar.py` as "a repressive
element; never inferred from a direction of effect alone", so no scored element could carry it: the
only call the project held was the derived one, a direction read off a model deletion. This module
attaches the one curated source that can satisfy that definition -- the ReSE screen of Pang and
Snyder (Nat Genet 2020, doi 10.1038/s41588-020-0578-5), as curated into RefSeq Functional Elements
and read from NCBI Gene -- to the project's scored elements, and gives the role ONLY where the assay
measured repression in the same cell line.

WHAT IT NEVER TOUCHES. The derived call is read, never written. Nothing here opens a scored element
file for writing, and the join itself is handed `(start, end, id)` tuples with no direction in them,
so the derived call cannot steer which elements receive the role. The derived direction enters only
afterwards, as a description of the elements the assay already chose.

THE RULE, decided and committed before the join was run (`decision()`):

  * Overlap. A fragment attaches to an element under the measured layer's production rule,
    `reciprocal_half`: the shared bases cover at least half of the element AND at least half of the
    fragment (`attribution/measured.py`, RECIPROCAL_OVERLAP). One shared base is reported and sets
    nothing; containment either way is reported and sets nothing.
  * One fragment, one element. A fragment that attaches to two or more elements under the rule sets
    the role on none of them: the assay measured one piece of DNA and cannot say which element
    carried the repression. It is reported.
  * Cell line. ReSE names K562 and HepG2; those map to the scored lines of the same names. A fragment
    validated in K562 gives the role in K562 only. The same element in HepG2 reads UNKNOWN with the
    reason that the source does not say whether the fragment was tested there.
  * Lines the source does not cover. GM12878 and IMR-90 have no ReSE record. Every element in them
    reads UNKNOWN with that reason, never a negative.
  * Absence. The source lists validated fragments only, not the fragments the screen tested without
    effect, so an element with no attached fragment reads UNKNOWN in a covered line too. No reading
    anywhere in this module means "no repressive role".

WHAT THE ROLE LICENSES. That a fragment overlapping the element under the rule repressed a reporter
in that cell line, integrated by lentivirus outside its native locus. It names no target gene, says
nothing about the element's effect in its own chromosome, and attaches no probability (review R4).
"""

from __future__ import annotations

import bisect
import json
import re
from collections import Counter, defaultdict
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import Any

from genomeos.attribution import silenceragree as sa
from genomeos.attribution.measured import min_side_overlap, reciprocal_overlap

ROOT = sa.ROOT
RESE_CACHE = sa.RESE_CACHE
ELEMENTS_DIR = sa.ELEMENTS_DIR
DECISION = "data/results/silencerjoin_decision.json"

# ---- the decision ----------------------------------------------------------------------------------

#: The role-setting overlap: the measured layer's production rule, at its own value. Pinned here so a
#: change to that constant is a visible break in this module's tests and not a silent re-decision.
ROLE_RULE = "reciprocal_half"
ROLE_OVERLAP = 0.5
#: Reported beside the rule, never used to set the role.
RECIPROCAL_SENSITIVITY = (0.25, 0.75)
MIN_SIDE_SENSITIVITY = 0.5

#: ReSE's cell-line name -> the project's scored line of the same cell line. Nothing else maps.
CELL_LINE_MAP = {"HepG2": "HepG2", "K562": "K562"}
SCORED_CELL_LINES = sa.PROJECT_CELL_LINES
COVERED = tuple(sorted(set(CELL_LINE_MAP.values())))

ROLE = "silencer"
UNKNOWN = "unknown"
#: R7 `evidence_status` with its qualifier: `measured` is the grammar's "an assay measured it", and the
#: qualifier names the assay. The source is named in full on every reading under `evidence`.
EVIDENCE_STATUS = "measured/ReSE"
EVIDENCE = (
    "ReSE lentiviral screen of repressive ability (cell survival under transcriptional repression of "
    "an apoptosis-inducing fusion), Pang and Snyder, Nat Genet 2020, doi 10.1038/s41588-020-0578-5, "
    "curated into RefSeq Functional Elements and read from NCBI Gene"
)

#: Why a reading is UNKNOWN. Each is a reason the role was not READ, and none is a negative.
LINE_NOT_COVERED = "line_not_covered"
ANOTHER_LINE_ONLY = "validated_in_another_line_only"
SEVERAL_ELEMENTS = "fragment_attaches_to_several_elements"
BELOW_RULE = "overlap_below_the_rule"
NOT_LISTED = "no_validated_fragment_attached"
REASONS = {
    LINE_NOT_COVERED: "ReSE validated silencers in K562 and HepG2 only and no record names this cell "
    "line, so the source cannot read a role here; this is UNKNOWN, not a negative",
    ANOTHER_LINE_ONLY: "a fragment attached under the rule was validated in another cell line only; "
    "the source does not say whether it was tested in this one, so the role here is not read",
    SEVERAL_ELEMENTS: "a fragment validated in this line overlaps this element but attaches to more "
    "than one element under the rule, so the assay cannot say which element carried the repression",
    BELOW_RULE: "a fragment validated in this line shares bases with this element but below the rule "
    "(under half of the element or under half of the fragment)",
    NOT_LISTED: "no fragment validated in this line attaches to this element; the source lists "
    "validated fragments only and not the fragments the screen tested without effect, so absence "
    "here is not a negative",
}
#: The precedence among the reasons of a covered line, most informative first.
REASON_ORDER = (ANOTHER_LINE_ONLY, SEVERAL_ELEMENTS, BELOW_RULE, NOT_LISTED)

#: The derived call, read only to describe the elements the assay chose.
REPRESSES, ACTIVATES = "represses", "activates"

NOT_A_TEST = (
    "A description and not a test. The registered check of the derived call against these same "
    "assay-validated silencers already ran (registration c3ebf70, result a5a43b3: 420 of 1,181, "
    "35.563% against a matched control's 37.707%, NO_SIGNAL). It is not re-run here, no control is "
    "drawn and no interval is computed, so these counts carry no reading of agreement."
)
WHAT_THE_ROLE_LICENSES = (
    "That a ReSE fragment attached to the element under the rule repressed a reporter in that cell "
    "line, in a lentiviral screen outside the fragment's native locus. That is what molecular_role "
    "silencer with evidence_status measured/ReSE states, in that line and no other."
)
WHAT_THE_ROLE_DOES_NOT_LICENSE = (
    "It names no target gene and no effect of the element in its own chromosome. It confirms nothing "
    "about the derived direction call, which is unchanged and is reported beside it. It attaches no "
    "probability (review R4). An element without the role is not shown to lack one: every reading "
    "without the role is UNKNOWN with its reason."
)


def decision() -> dict[str, Any]:
    """The join rule, the cell-line mapping and what is reported, as committed BEFORE the join.

    Holds no join outcome: it is built from constants and from the ReSE records alone, and never
    opens a scored element.
    """
    kept, accounting = sa.rese_intervals(sa.rese_records())
    lengths = sorted(f["end"] - f["start"] + 1 for f in kept)
    return {
        "question": "Which scored elements may carry molecular_role silencer from an assay, in which "
        "cell line, beside (never in place of) the derived direction call?",
        "source": sa.SOURCE,
        "source_reuse": "the probe's own cache, fetched 2026-10-03 by scripts/silenceragree_fetch.py "
        "and pinned by sha256 in data/results/silenceragree.json; read in place, not re-fetched, so "
        "this join and the registered check read the same bytes",
        "silencerdb_not_used": "SilencerDB's validated tier is coordinate-identical to ReSE and states "
        "no licence anywhere; ReSE alone is used",
        "fragments_before_the_join": {
            "joinable": len(kept),
            "by_validation_line": dict(Counter("+".join(f["cells"]) or "none" for f in kept).most_common()),
            "length_bp": {
                "min": lengths[0],
                "median": lengths[len(lengths) // 2],
                "max": lengths[-1],
            },
            "accounting": accounting,
        },
        "rule": rule(),
        "cell_line_mapping": {
            "rese_to_scored_line": CELL_LINE_MAP,
            "scored_lines": list(SCORED_CELL_LINES),
            "covered": list(COVERED),
            "not_covered": [c for c in SCORED_CELL_LINES if c not in COVERED],
            "parsed_from": "the record summary's clause 'identified as a functional silencer in <cell "
            "line>' (silenceragree.CELL_CLAUSE); a fragment naming both lines gives the role in both, "
            "a fragment naming neither gives it in none and is counted",
        },
        "readings": {
            "role": {"molecular_role": ROLE, "evidence_status": EVIDENCE_STATUS, "evidence": EVIDENCE},
            "unknown": {"molecular_role": UNKNOWN, "evidence_status": UNKNOWN, "reasons": REASONS},
            "reason_precedence": list(REASON_ORDER),
        },
        "reported": REPORTED,
        "derived_call": NOT_A_TEST,
        "what_the_role_licenses": WHAT_THE_ROLE_LICENSES,
        "what_the_role_does_not_license": WHAT_THE_ROLE_DOES_NOT_LICENSE,
    }


def rule() -> dict[str, Any]:
    return {
        "role_rule": ROLE_RULE,
        "role_overlap": ROLE_OVERLAP,
        "statement": "a ReSE fragment [fs, fe) attaches to a scored element [es, ee) when the shared "
        "bases are at least ROLE_OVERLAP of the element's width AND at least ROLE_OVERLAP of the "
        "fragment's width: the measured layer's production rule (attribution/measured.py "
        "RECIPROCAL_OVERLAP, reciprocal_half), adopted rather than invented",
        "one_fragment_one_element": "a fragment attaching to two or more elements under the rule sets "
        "the role on none of them and is counted",
        "fragment_coordinates": "half-open [chrstart, chrstop + 1) from NCBI genomicinfo, which gives "
        "the 0-based first and last base: 4,760 of 4,992 GRCh38 spans have chrstop - chrstart + 1 "
        "equal to the 1-based inclusive length of the record's GRCh37 description (silenceragree "
        "accounting). The GRCh37 description is never joined against",
        "element_coordinates": "half-open [start, end) as stored in the scored element files",
        "does_not_set_the_role": [
            "a shared base below the rule",
            "a fragment containing the element, or contained in it, below the rule",
            "a fragment validated in another cell line",
            "a fragment attaching to several elements",
            "the sensitivity thresholds below",
        ],
        "sensitivity_reported_only": {
            "reciprocal": list(RECIPROCAL_SENSITIVITY),
            "min_side": MIN_SIDE_SENSITIVITY,
            "any_shared_base": True,
        },
        "joined_population": "every scored element in data/knowledge/alphagenome/all_elements, whether "
        "or not a target is named: the role is about the element's sequence, not about the prediction",
    }


#: What the result reports, as (field, what) pairs, so the field names a join fills appear here only as
#: values and the decision cannot be mistaken for, or refused as, a result.
REPORTED = [
    {
        "field": "by_cell_line",
        "what": "per scored line: elements given the role; elements read UNKNOWN, split by reason; "
        "fragments validated in the line, placed on one element, on several, below the rule only, and "
        "on no scored element",
    },
    {
        "field": "derived_call_on_the_role_elements",
        "what": "per covered line, of the elements given the role: how many have a named predicted "
        "target, how many of those carry the derived repressive direction (the derived silencer_like "
        "call) and how many the activating one, and how many have no named target",
    },
    {
        "field": "coverage_gap",
        "what": "per scored line: the scored elements whose role the source cannot read, and why, and "
        "the fragments validated there that no scored element carries; for a line the source does not "
        "cover, all of them",
    },
    {
        "field": "sensitivity",
        "what": "elements that would receive the role at reciprocal 0.25 and 0.75, at min-side 0.5, and "
        "at any shared base, per covered line; reported, never applied",
    },
    {
        "field": "elements",
        "what": "one row per scored element sharing at least one base with any ReSE fragment: its "
        "overlaps, its derived direction as committed, and its reading in each of the four scored lines",
    },
]

# ---- the fragments ---------------------------------------------------------------------------------


def fragments(
    records: dict[str, dict[str, Any]] | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """The joinable ReSE fragments, half-open, with their validation lines mapped to scored lines."""
    kept, accounting = sa.rese_intervals(records if records is not None else sa.rese_records())
    out = []
    for f in kept:
        out.append(
            {
                "fragment": f["uid"],
                "chrom": f["chrom"],
                "start": f["start"],
                "end": f["end"] + 1,
                "lines": sorted({CELL_LINE_MAP[c] for c in f["cells"] if c in CELL_LINE_MAP}),
            }
        )
    return out, accounting


# ---- the scored elements -----------------------------------------------------------------------------


def element_spans(chrom: str, elements_dir: Path = ELEMENTS_DIR) -> list[tuple[int, int, str]]:
    """`(start, end, id)` for one chromosome, sorted. NO DIRECTION: the join cannot see the derived call."""
    path = elements_dir / f"{chrom}.json"
    if not path.exists():
        return []
    return sorted((r["start"], r["end"], r["id"]) for r in json.loads(path.read_text()))


def derived_directions(
    chrom: str, ids: Iterable[str], elements_dir: Path = ELEMENTS_DIR
) -> dict[str, str | None]:
    """The committed derived direction of each named element: `represses`, `activates`, or None when no
    target is named. Read only AFTER the join, to describe what the assay chose."""
    want = set(ids)
    if not want:
        return {}
    path = elements_dir / f"{chrom}.json"
    out: dict[str, str | None] = {}
    for r in json.loads(path.read_text()):
        if r["id"] in want:
            pred = r.get("predicted")
            out[r["id"]] = pred.get("action") if pred else None
    return out


# ---- the join ----------------------------------------------------------------------------------------


def pairs_on(frags: list[dict[str, Any]], spans: list[tuple[int, int, str]]) -> list[dict[str, Any]]:
    """Every (fragment, element) pair sharing at least one base on one chromosome, with its fractions."""
    if not spans:
        return []
    starts = [s for s, _, _ in spans]
    widest = max(e - s for s, e, _ in spans)
    out = []
    for f in frags:
        fs, fe = f["start"], f["end"]
        j = bisect.bisect_left(starts, fs - widest)
        while j < len(spans) and spans[j][0] < fe:
            es, ee, eid = spans[j]
            j += 1
            shared = min(ee, fe) - max(es, fs)
            if shared <= 0:
                continue
            out.append(
                {
                    "fragment": f["fragment"],
                    "lines": f["lines"],
                    "element": eid,
                    "shared_bp": shared,
                    "reciprocal": round(reciprocal_overlap(es, ee, fs, fe), 4),
                    "min_side": round(min_side_overlap(es, ee, fs, fe), 4),
                    "fragment_within_element": es <= fs and fe <= ee,
                    "element_within_fragment": fs <= es and ee <= fe,
                }
            )
    return out


def attaches_at(fraction: float) -> Callable[[dict[str, Any]], bool]:
    return lambda p: p["reciprocal"] >= fraction


def attaches_min_side(fraction: float) -> Callable[[dict[str, Any]], bool]:
    return lambda p: p["min_side"] >= fraction


def any_shared_base(p: dict[str, Any]) -> bool:
    return p["shared_bp"] > 0


ROLE_TEST = attaches_at(ROLE_OVERLAP)


def label_pairs(pairs: list[dict[str, Any]], test: Callable[[dict[str, Any]], bool] = ROLE_TEST) -> None:
    """Mark each pair `unique`, `several` or `below_rule` under `test`, fragment by fragment, in place."""
    by_fragment: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for p in pairs:
        by_fragment[p["fragment"]].append(p)
    for ps in by_fragment.values():
        hits = [p for p in ps if test(p)]
        for p in ps:
            if not test(p):
                p["attaches"] = "below_rule"
            else:
                p["attaches"] = "unique" if len(hits) == 1 else "several"


def role_elements(
    pairs: list[dict[str, Any]], line: str, test: Callable[[dict[str, Any]], bool] = ROLE_TEST
) -> set[str]:
    """Elements a fragment validated in `line` attaches to, and to no other element, under `test`."""
    by_fragment: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for p in pairs:
        by_fragment[p["fragment"]].append(p)
    out: set[str] = set()
    for ps in by_fragment.values():
        if line not in ps[0]["lines"]:
            continue
        hits = [p for p in ps if test(p)]
        if len(hits) == 1:
            out.add(hits[0]["element"])
    return out


def unknown(reason: str) -> dict[str, Any]:
    return {"molecular_role": UNKNOWN, "evidence_status": UNKNOWN, "reason": reason}


def reading(overlaps: list[dict[str, Any]], line: str) -> dict[str, Any]:
    """One element's reading in one scored line, from the labelled pairs it is part of.

    The role is given only by a fragment validated in THIS line that attaches to this element alone.
    Every other case is UNKNOWN with its reason; there is no negative value.
    """
    if line not in COVERED:
        return unknown(LINE_NOT_COVERED)
    same = sorted(p["fragment"] for p in overlaps if line in p["lines"] and p["attaches"] == "unique")
    if same:
        return {
            "molecular_role": ROLE,
            "evidence_status": EVIDENCE_STATUS,
            "evidence": EVIDENCE,
            "cell_line": line,
            "fragments": same,
        }
    reasons = set()
    for p in overlaps:
        if line not in p["lines"]:
            if p["attaches"] == "unique" and p["lines"]:
                reasons.add(ANOTHER_LINE_ONLY)
        elif p["attaches"] == "several":
            reasons.add(SEVERAL_ELEMENTS)
        elif p["attaches"] == "below_rule":
            reasons.add(BELOW_RULE)
    return unknown(next((r for r in REASON_ORDER if r in reasons), NOT_LISTED))


def element_rows(
    pairs: list[dict[str, Any]],
    spans: dict[str, tuple[str, int, int]],
    directions: dict[str, str | None],
) -> list[dict[str, Any]]:
    """One row per element sharing a base with any fragment: overlaps, derived direction, four readings."""
    by_element: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for p in pairs:
        by_element[p["element"]].append(p)
    rows = []
    for eid in sorted(by_element, key=lambda e: (_chrom_key(spans[e][0]), spans[e][1], e)):
        chrom, start, end = spans[eid]
        ovs = sorted(by_element[eid], key=lambda p: p["fragment"])
        rows.append(
            {
                "id": eid,
                "chrom": chrom,
                "start": start,
                "end": end,
                "derived_direction": directions.get(eid),
                "overlaps": [
                    {
                        k: p[k]
                        for k in ("fragment", "lines", "shared_bp", "reciprocal", "min_side", "attaches")
                    }
                    for p in ovs
                ],
                "by_cell_line": {line: reading(ovs, line) for line in SCORED_CELL_LINES},
            }
        )
    return rows


def _chrom_key(chrom: str) -> tuple[int, str]:
    tail = chrom.removeprefix("chr")
    return (int(tail), "") if tail.isdigit() else (100, tail)


# ---- the counts --------------------------------------------------------------------------------------


def by_cell_line(
    rows: list[dict[str, Any]],
    pairs: list[dict[str, Any]],
    frags: list[dict[str, Any]],
    elements_scored: int,
) -> dict[str, Any]:
    """Per scored line: the role given, the UNKNOWN readings by reason, and where the fragments went."""
    out: dict[str, Any] = {}
    paired = defaultdict(set)
    for p in pairs:
        paired[p["fragment"]].add(p["attaches"])
    for line in SCORED_CELL_LINES:
        given = sum(1 for r in rows if r["by_cell_line"][line]["molecular_role"] == ROLE)
        if line not in COVERED:
            out[line] = {
                "source_covers_the_line": False,
                "elements_scored": elements_scored,
                "elements_given_the_role": 0,
                "elements_read_unknown": elements_scored,
                "unknown_by_reason": {LINE_NOT_COVERED: elements_scored},
                "fragments_validated_in_the_line": 0,
            }
            continue
        reasons = Counter(
            r["by_cell_line"][line]["reason"]
            for r in rows
            if r["by_cell_line"][line]["molecular_role"] != ROLE
        )
        reasons[NOT_LISTED] += elements_scored - len(rows)
        mine = [f for f in frags if line in f["lines"]]
        placed = Counter()
        for f in mine:
            kinds = paired.get(f["fragment"], set())
            if "unique" in kinds:
                placed["placed_on_one_element"] += 1
            elif "several" in kinds:
                placed["attached_to_several_elements"] += 1
            elif "below_rule" in kinds:
                placed["sharing_bases_below_the_rule_only"] += 1
            else:
                placed["on_no_scored_element"] += 1
        out[line] = {
            "source_covers_the_line": True,
            "elements_scored": elements_scored,
            "elements_given_the_role": given,
            "elements_read_unknown": elements_scored - given,
            "unknown_by_reason": dict(sorted(reasons.items())),
            "fragments_validated_in_the_line": len(mine),
            "fragments": {
                k: placed.get(k, 0)
                for k in (
                    "placed_on_one_element",
                    "attached_to_several_elements",
                    "sharing_bases_below_the_rule_only",
                    "on_no_scored_element",
                )
            },
        }
    return out


def coverage_gap(lines: dict[str, Any]) -> dict[str, Any]:
    """Per scored line, the scored elements whose role the source cannot read, and why."""
    out: dict[str, Any] = {}
    for line, c in lines.items():
        gap = {
            "source_covers_the_line": c["source_covers_the_line"],
            "elements_scored": c["elements_scored"],
            "elements_whose_role_is_not_read": c["elements_read_unknown"],
            "share_not_read": round(c["elements_read_unknown"] / c["elements_scored"], 6)
            if c["elements_scored"]
            else None,
            "by_reason": c["unknown_by_reason"],
        }
        if c["source_covers_the_line"]:
            gap["validated_fragments_no_scored_element_carries"] = c["fragments"]["on_no_scored_element"]
        out[line] = gap
    return out


def derived_call_on_role_elements(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Per covered line, the committed derived direction of the elements given the role. Not a test."""
    out: dict[str, Any] = {"reading": NOT_A_TEST}
    for line in COVERED:
        got = [r for r in rows if r["by_cell_line"][line]["molecular_role"] == ROLE]
        named = [r for r in got if r["derived_direction"] is not None]
        out[line] = {
            "elements_given_the_role": len(got),
            "with_a_named_predicted_target": len(named),
            "carrying_the_derived_repressive_direction": sum(
                r["derived_direction"] == REPRESSES for r in named
            ),
            "carrying_the_derived_activating_direction": sum(
                r["derived_direction"] == ACTIVATES for r in named
            ),
            "with_no_named_predicted_target": len(got) - len(named),
        }
    both = [r for r in rows if all(r["by_cell_line"][line]["molecular_role"] == ROLE for line in COVERED)]
    out["given_the_role_in_both_covered_lines"] = len(both)
    return out


def sensitivity(pairs: list[dict[str, Any]]) -> dict[str, Any]:
    """Elements that WOULD receive the role at other thresholds, per covered line. Reported, never applied."""
    tests = {f"reciprocal_{x}": attaches_at(x) for x in sorted((*RECIPROCAL_SENSITIVITY, ROLE_OVERLAP))}
    tests[f"min_side_{MIN_SIDE_SENSITIVITY}"] = attaches_min_side(MIN_SIDE_SENSITIVITY)
    tests["any_shared_base"] = any_shared_base
    return {
        "applied": f"reciprocal_{ROLE_OVERLAP}",
        "by_threshold": {
            name: {line: len(role_elements(pairs, line, t)) for line in COVERED} for name, t in tests.items()
        },
    }


def join(
    frags: list[dict[str, Any]], elements_dir: Path = ELEMENTS_DIR
) -> tuple[list[dict[str, Any]], dict[str, tuple[str, int, int]], int, list[str]]:
    """The pairs over every scored chromosome, the spans of the paired elements, the count of elements
    scored, and the chromosomes read. Each chromosome is read once and released."""
    by_chrom: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for f in frags:
        by_chrom[f["chrom"]].append(f)
    pairs: list[dict[str, Any]] = []
    spans_of: dict[str, tuple[str, int, int]] = {}
    scored = 0
    chroms = sorted((p.stem for p in elements_dir.glob("chr*.json")), key=_chrom_key)
    for chrom in chroms:
        spans = element_spans(chrom, elements_dir)
        scored += len(spans)
        got = pairs_on(by_chrom.get(chrom, []), spans)
        hit = {p["element"] for p in got}
        for s, e, eid in spans:
            if eid in hit:
                spans_of[eid] = (chrom, s, e)
        pairs.extend(got)
    label_pairs(pairs)
    return pairs, spans_of, scored, chroms


def build(elements_dir: Path = ELEMENTS_DIR) -> dict[str, Any]:
    """The whole layer: the join, the rows, the counts. Reads the derived direction only after the join."""
    frags, accounting = fragments()
    pairs, spans_of, scored, chroms = join(frags, elements_dir)
    directions: dict[str, str | None] = {}
    by_chrom = defaultdict(set)
    for eid, (chrom, _, _) in spans_of.items():
        by_chrom[chrom].add(eid)
    for chrom, ids in sorted(by_chrom.items()):
        directions.update(derived_directions(chrom, ids, elements_dir))
    rows = element_rows(pairs, spans_of, directions)
    lines = by_cell_line(rows, pairs, frags, scored)
    on_scored = {f["chrom"] for f in frags} & set(chroms)
    return {
        "rule": rule(),
        "population": {
            "elements_scored": scored,
            "chromosomes_read": chroms,
            "elements_sharing_a_base_with_a_fragment": len(rows),
        },
        "fragments": {
            "accounting": accounting,
            "joinable": len(frags),
            "on_a_scored_chromosome": sum(1 for f in frags if f["chrom"] in on_scored),
            "with_no_validation_line_parsed": sum(1 for f in frags if not f["lines"]),
            "pairs_sharing_a_base": len(pairs),
            "pairs_by_label": dict(sorted(Counter(p["attaches"] for p in pairs).items())),
        },
        "by_cell_line": lines,
        "coverage_gap": coverage_gap(lines),
        "derived_call_on_the_role_elements": derived_call_on_role_elements(rows),
        "sensitivity": sensitivity(pairs),
        "elements": rows,
    }


# ---- the layer as BioLang ----------------------------------------------------------------------------


def bio_block(row: dict[str, Any], line: str) -> str | None:
    """The reading as a BioLang `element` block named `<id>_rese_<line>`, beside the predicted `<id>`.

    Written only where the role was given. Two names, never one, as the measured layer does it: the
    predicted block and its derived alternatives are not touched.
    """
    r = row["by_cell_line"][line]
    if r["molecular_role"] != ROLE:
        return None
    name = re.sub(r"[^A-Za-z0-9_]", "_", f"{row['id']}_rese_{line}")
    note = f"ReSE {line}, Pang and Snyder 2020, NCBI Gene {' '.join(r['fragments'])}"
    return (
        f"element {name} {{\n"
        f"  locus: {row['chrom']}:{row['start']}-{row['end']}\n"
        f"  molecular_role: {ROLE}\n"
        f"  evidence_status: {EVIDENCE_STATUS}\n"
        f'  evidence: experimental "{note}"\n'
        f"}}"
    )
