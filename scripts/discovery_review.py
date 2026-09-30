# SPDX-License-Identifier: AGPL-3.0-or-later
"""Audit A's adoption review, discovery only: which screened CRISPRi observations audit A's policy would
make discoverable as a set-valued attachment, with union coverage, multiplicity and a shifted screen.

    uv run python scripts/discovery_review.py

Lane-discover, 2026-09-29, for the coordinator (genomeos-9f), under the external reviewer's amended rules
for A's adoption review. It is bounded, read-only with respect to production, recorded apart, and it
stops after its registered report. It informs one decision only, Decision 1 (discovery): whether
evidence may be discoverable as a set-valued attachment. It proposes no implementation, changes no
matcher and proposes no broader verdict eligibility: verdicts stay under the current rule.

The universe is every screened observation of the placement census (`1253d69`, `32af2cb`), under both
of its rules: C4's reciprocal 0.5 rule and audit A's registered policy (an overlap of at least half of
the smaller interval). The census and audit A are imported, never edited; so are ablation.py,
holdout.py, measured.py and correctness.py. Every outcome it reads was read before (PRIOR_EXPOSURE), so
it is internal development evidence. No model request is made and nothing is downloaded.

The review's definitions (REGISTRATION) were committed with this file before any count of it was taken.
Writes data/results/discovery_review.json through the result contract.
"""

from __future__ import annotations

import importlib.util
import json
import resource
import sys
import time
from collections import Counter, defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np

from genomeos import manifest as mf
from genomeos.attribution import correctness as cr
from genomeos.attribution import holdout as ho
from genomeos.attribution import measured as ms
from genomeos.results import RESULTS_DIR, save_result


def _census() -> Any:
    """The placement census (scripts/placement_census.py), imported and never edited; it imports audit A
    (scripts/placement_audit.py) the same way. Its rules, observations and verdicts are the ones used."""
    if "placement_census" in sys.modules:
        return sys.modules["placement_census"]
    spec = importlib.util.spec_from_file_location(
        "placement_census", Path(__file__).resolve().parent / "placement_census.py"
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


pc = _census()
pa = pc.pa

NAME = "discovery_review"
REGISTERED = "2026-09-29"

# ==================================================================================================
# The registration (2026-09-29, lane-discover), committed before any count of this review was taken.
# ==================================================================================================
CURRENT, POLICY = pc.CURRENT, pc.POLICY
SHIFTS = (pa.SHIFT_BP, -pa.SHIFT_BP)  # audit A's positional control: +20 kb and -20 kb, the mean of the two
SCREEN_THRESHOLD = 1.5

# the change classes of the changed set, mutually exclusive, first match wins
LOST = "attachment_or_candidate_lost"
NEWLY_ATTACHED = "newly_attached"
LINKS_ADDED = "attached_links_added"
CANDIDATES_ADDED_ATTACHED = "candidates_added_attached_links_unchanged"
CANDIDATES_ADDED_UNATTACHED = "candidates_added_unattached_under_both"
CHANGES = (LOST, NEWLY_ATTACHED, LINKS_ADDED, CANDIDATES_ADDED_ATTACHED, CANDIDATES_ADDED_UNATTACHED)
#: the classes whose attached links differ between the rules: the evidence discovery would add
DISCOVERED = (NEWLY_ATTACHED, LINKS_ADDED)
CHANGE_MEANING = {
    LOST: "a candidate or an attached link of the current rule is missing under the policy (must be empty: "
    "discovery condition 1)",
    NEWLY_ATTACHED: "no attached link under the current rule, at least one under the policy",
    LINKS_ADDED: "attached under both rules, with links the current rule does not attach",
    CANDIDATES_ADDED_ATTACHED: "the same attached links under both rules, more admitted registry candidates "
    "under the policy (the resolution may move from unique to ambiguous)",
    CANDIDATES_ADDED_UNATTACHED: "no attached link under either rule, more admitted registry candidates "
    "under the policy",
}

# the stated-cell agreement of one observation under one rule, over its attached links
MATCH = "stated_cell_match"
OTHER_ONLY = "other_cell_only"
NO_STATED = "no_stated_cell"
NOT_ATTACHED = "not_attached"
AGREEMENTS = (MATCH, OTHER_ONLY, NO_STATED, NOT_ATTACHED)
AGREEMENT_MEANING = {
    MATCH: "a compiled predicted rule of an attached link (its element and the measured gene) states the "
    "cell the screen was run in (correctness.norm_cell equal)",
    OTHER_ONLY: "the attached links' rules state cells, none of them the screened cell",
    NO_STATED: "no rule of an attached link states a cell (correctness.stated false)",
    NOT_ATTACHED: "no attached link",
}

# the direction description of one observation against its attached links under the policy: the
# census's verdict renamed, because under this review it describes and authorises nothing
DESCRIPTION = {
    pc.SUPPORTED: "agrees",
    pc.REFUTED_OPPOSITE: "opposite",
    pc.SPLIT: "split_direction",
    pc.REFUTED_NULL: "well_powered_null",
    pc.UNASSESSED_UNDERPOWERED: "underpowered_null",
    pc.UNASSESSED_MISSING: "missing",
    pc.UNATTACHED: "not_attached",
}
UNRESOLVED = "unresolved"
UNIQUE_AMONG_ADMITTED = "unique_among_all_admitted_registry_candidates"

# the screen's per-stratum reading
PASS, FAIL = "passes", "fails"
UNDEFINED = "undefined_shifted_count_zero"
NOTHING = "nothing_discovered"
SCREEN_QUANTITY = "attached_links_differ"
ALSO_SCREENED = ("attached_under_the_policy", "in_the_changed_set")

# the registered readings of Decision 1
MERITS = "merits_a_separately_reviewed_implementation_proposal"
RESTRICTED = "merits_a_separately_reviewed_proposal_only_in_the_passing_strata"
DOES_NOT = "does_not_merit_an_implementation_proposal"

FRACTION_BINS = tuple(round(k / 10, 1) for k in range(11))

#: the coordinator's in-memory reading, never committed, recorded here to be reproduced (or not)
PROVISIONAL = {
    "newly_supported_unique": {"n": 37, MATCH: 24, OTHER_ONLY: 13, "matches_all_in": "K562"},
    "newly_supported_ambiguous": {"n": 79, "with_a_stated_cell_match": 39},
    "current_rule_supported": {"n": 93, MATCH: 39, "not_a_match": 54},
}

REGISTRATION = {
    "registered": REGISTERED,
    "approved_as": "A's adoption review, discovery only: bounded, read-only with respect to production, "
    "recorded apart, stopping after its registered report (the external reviewer's amended rules, relayed "
    "by the coordinator genomeos-9f)",
    "status": "internal development evidence: every outcome read here was read before (prior_exposure)",
    "decision_informed": "Decision 1 (discovery) only: whether evidence may be discoverable as a set-valued "
    "attachment. No other decision; no implementation is proposed and no matcher is changed",
    "universe": "all screened CRISPRi observations of the placement census (every valid row of both "
    "benchmark files, 14,734 in the census), each under both of the census's rules and each keeping its "
    "experimental observation id (placement_census.observation_id). The census's observe() is imported "
    "and applied unchanged; the rules are C4's reciprocal 0.5 rule and audit A's policy (overlap at least "
    "half of the smaller interval)",
    "changed_set": "every observation whose candidate set, attached links or resolution differs between the "
    "two rules. Each is placed in one change class, first match wins",
    "change_classes": CHANGE_MEANING,
    "coverage": {
        "union_coverage": "for one observation under one rule: the admitted registry candidates' intervals, "
        "each clipped to the tested interval, merged where they overlap or abut, and the merged length "
        "divided by the tested width. Overlapping candidates are counted once, never summed. 1 minus it is "
        "the fraction of the tested interval no admitted registry candidate covers, which stays unresolved",
        "candidate_fraction": "for each admitted registry candidate: its overlap with the tested interval "
        "divided by the candidate's own width",
        "multiplicity": "the number of admitted registry candidates (the candidate set's size) under "
        "each rule",
        "reported_for": "every changed observation, under both rules (the current rule's candidates are a "
        "subset of the policy's, checked); per-observation rows carry the full admitted-candidate set",
        "summaries": "distributions and summaries only: n, mean, minimum, 10th, 25th, 50th, 75th and 90th "
        "percentiles (numpy linear interpolation), maximum, and counts in tenths of the fraction with "
        "exactly 1 apart; multiplicity as counts per value. By change class, and for the discovered "
        "evidence by resolution. Tested width is summarised in the same continuous way. No width "
        "boundary in base pairs is used anywhere",
    },
    "descriptions_for_the_changed_set": {
        "stated_cell_agreement": AGREEMENT_MEANING,
        "stated_cell_source": "correctness.compiled_claims, target-axis claims of the compiled predicted "
        "rules, keyed by (element, gene), their cells read with correctness.stated and matched with "
        "correctness.norm_cell. A predicted link stores no cell; the rule beside it does",
        "direction_relative_to_the_claim": "the census's verdict under the policy, renamed: agrees, "
        "opposite, split direction; with well-powered nulls, underpowered nulls, missing effect sizes and "
        "significant increases each counted apart (outcome by description)",
        "breakdowns": "by study (the Dataset column), by cell, and by study and cell",
        "as": "descriptions only: supports against contradictions are reported descriptively and authorise "
        "nothing",
    },
    "verdict_eligibility": "unchanged: the only eligible verdict of an observation is the one under the "
    "current rule. No broader eligibility is proposed",
    "newly_ambiguous": "an observation that is a unique-element assignment under the current rule and "
    "ambiguous under the policy keeps its historical verdict under the current rule and is marked "
    f"'{UNRESOLVED}' under the new rule; both are shown, one row per observation",
    "discovery_conditions": [
        "no current attachment is lost (every current candidate and attached link is kept under the "
        "policy) and the denominators are unchanged (the same observations, the same predicted links as "
        "C4's committed count)",
        "each observation appears once, with its full admitted-candidate set and its coverage fractions "
        "(the per-observation rows: ids unique, one row per changed observation, one coverage fraction per "
        "admitted candidate, a union coverage per row)",
    ],
    "screen": {
        "what": "a diagnostic screen for discovery, not an estimate of false attribution",
        "quantity": f"'{SCREEN_QUANTITY}': the observation's attached links under the policy differ from "
        "those under the current rule (the change classes newly_attached and attached_links_added), the "
        "evidence discovery would add",
        "real": "the fraction of all observations of the stratum, regardless of outcome, with the quantity "
        "true at the tested interval",
        "shifted": f"the same fraction with every tested interval moved by +{pa.SHIFT_BP} bp and by "
        f"-{pa.SHIFT_BP} bp (the same gene, cell and width; both rules re-applied), the mean of the two",
        "enrichment": "real divided by shifted",
        "threshold": f"the stratum passes when the enrichment is at least {SCREEN_THRESHOLD}; fails below "
        f"it; reads '{UNDEFINED}' when the shifted count is zero and the real one is not (the counts are "
        f"reported); reads '{NOTHING}' when the real count is zero",
        "strata": "each study (Dataset column) and cell combination present, and all observations pooled",
        "also_reported_without_a_threshold": [
            "the same enrichment for attachment under the policy at all",
            "the same enrichment for membership of the changed set",
        ],
        "caveat": "moving an interval by 20 kb can change its distance to the measured gene and the local "
        "density of registry elements and predicted links, so the shifted rate is a positional comparison, "
        "not a null for the attachment's truth",
    },
    "reading": {
        MERITS: "both discovery conditions hold, the pooled screen passes, and every stratum where the "
        "screen is defined passes",
        RESTRICTED: "both discovery conditions hold and the pooled screen passes, but a stratum where it is "
        "defined fails; the failing strata are named",
        DOES_NOT: "a discovery condition fails, or the pooled screen fails or is undefined",
        "note": "a proposal, if merited, is separate and separately reviewed; this review is not it",
    },
    "stated_cell_record": {
        "what": "the coordinator's in-memory reading (never committed), recomputed and recorded: the "
        "census's newly supported observations (supported under the policy, not under the current rule) "
        "split by resolution and stated-cell agreement, and the current rule's supported observations "
        "split by stated-cell agreement, checked against the v2/v3 direction split "
        "(attribution_correctness_v3.json, in_context_side_by_side, v2 activity)",
        "provisional_values": PROVISIONAL,
        "scope": "C4's any-cell attachment and direction agreement; not v3 stated-context support",
    },
    "wording": "'unique' means one among all admitted registry candidates (ENCODE cCRE v3 elements that "
    "meet the rule), not an experimentally resolved single cause; unregistered sequence inside a "
    "perturbation remains unresolved",
    "gates": [
        "the census reproduces: both rules' observation tables and the verdict transitions equal "
        "placement_census.json as committed (32af2cb)",
        "C4's predicted-link count reproduces (440,377) and every compiled element is a registry element "
        "with the registry's coordinates",
    ],
    "cost": "wall time, CPU time and peak memory of the run; model requests and downloads (zero)",
    "out_of_scope": [
        "an implementation proposal, any matcher change, any regenerated result",
        "any broader verdict eligibility",
        "any other endpoint (lentiMPRA, VISTA, saturation mutagenesis, GTEx)",
    ],
}
PRIOR_EXPOSURE = [
    *pc.PRIOR_EXPOSURE,
    {"commit": "1253d69", "read": "the placement census's registration (2026-09-29)"},
    {
        "commit": "32af2cb",
        "read": "the placement census's result: both rules' observation tables and the difference (180 "
        "newly attached, 116 newly supported, 13 unique assignments becoming ambiguous) (2026-09-29)",
    },
    {
        "commit": "a55f540",
        "read": "the census counts with their transitions: 116 newly supported, 37 unique and 79 "
        "ambiguous, and their any-cell scope (2026-09-29)",
    },
    {
        "commit": "4567219",
        "read": "the v2 direction split: 39 activity claims judged in the stated cell, 54 agreeing only in "
        "another cell (2026-09-29)",
    },
    {"commit": "d1e09ae", "read": "the v3 target split: 39 in context, 59 elsewhere only (2026-09-29)"},
    {
        "commit": None,
        "read": "the coordinator's in-memory join of the census to the stated cells, never committed: of "
        "the 37 newly supported unique observations 24 match the stated cell (all K562) and 13 are other "
        "cell only; of the 79 newly supported ambiguous 39 have a stated-cell match; under the current rule "
        "39 of 93 supported match and 54 do not (2026-09-29, genomeos-9f)",
    },
    {
        "commit": None,
        "read": "this lane read placement_census.json's committed difference and rule summaries, and the "
        "v2 and v3 side-by-side counts, before this registration (2026-09-29, lane-discover)",
    },
]


def registration() -> dict[str, Any]:
    return REGISTRATION


# --- coverage -------------------------------------------------------------------------------------
def union_coverage(start: int, end: int, intervals: Iterable[tuple[int, int]]) -> float | None:
    """The fraction of [start, end) covered by the union of `intervals`: overlapping or abutting
    intervals are merged, so a base two candidates share is counted once. None for a zero-width interval."""
    width = end - start
    if width <= 0:
        return None
    clipped = sorted((max(a, start), min(b, end)) for a, b in intervals if min(b, end) > max(a, start))
    covered = 0
    run_s = run_e = None
    for a, b in clipped:
        if run_e is None or a > run_e:
            if run_e is not None:
                covered += run_e - run_s
            run_s, run_e = a, b
        else:
            run_e = max(run_e, b)
    if run_e is not None:
        covered += run_e - run_s
    return covered / width


def candidate_fraction(start: int, end: int, a: int, b: int) -> float | None:
    """The fraction of the candidate [a, b) that the tested interval [start, end) covers."""
    if b <= a:
        return None
    return pa.overlap(a, b, start, end) / (b - a)


def candidate_intervals(registry: Any, o: Any) -> list[tuple[str, int, int]]:
    """The observation's admitted registry candidates with their coordinates, in id order."""
    want = set(o.candidates)
    got = sorted(
        (i, a, b) for a, b, i in registry.touching(o.pair.chrom, o.pair.start, o.pair.end) if i in want
    )
    if [g[0] for g in got] != sorted(want):
        raise ValueError(f"{o.id}: a candidate is not in the registry index")
    return got


# --- the changed set ------------------------------------------------------------------------------
def change_class(c: Any, p: Any) -> str | None:
    """The change class of one observation (current `c`, policy `p`), or None when it is unchanged."""
    if c.id != p.id:
        raise ValueError(f"{c.id} against {p.id}: not the same observation")
    c_att, p_att = set(c.attached), set(p.attached)
    c_cand, p_cand = set(c.candidates), set(p.candidates)
    if not c_att <= p_att or not c_cand <= p_cand:
        return LOST
    if c_att != p_att:
        return NEWLY_ATTACHED if not c_att else LINKS_ADDED
    if c_cand != p_cand:
        return CANDIDATES_ADDED_ATTACHED if p_att else CANDIDATES_ADDED_UNATTACHED
    if c.resolution != p.resolution:
        raise ValueError(f"{c.id}: the resolution changed with neither set changing")
    return None


def is_changed(c: Any, p: Any) -> bool:
    return change_class(c, p) is not None


# --- stated cells ---------------------------------------------------------------------------------
def stated_cells(claims: Iterable[Any]) -> tuple[dict[tuple[str, str], frozenset[str]], set[tuple[str, str]]]:
    """(element, gene) -> the stated cells of its compiled predicted rules, and every (element, gene) that
    has a rule at all. Target-axis claims only: one per compiled predicted rule."""
    cells: dict[tuple[str, str], set[str]] = defaultdict(set)
    ruled: set[tuple[str, str]] = set()
    for c in claims:
        if c.axis != cr.TARGET:
            continue
        key = (c.element, c.gene)
        ruled.add(key)
        if cr.stated(c.cell):
            cells[key].add(c.cell)
    return {k: frozenset(v) for k, v in cells.items()}, ruled


def agreement(o: Any, links: Sequence[Any], stated: dict[tuple[str, str], frozenset[str]]) -> str:
    """The stated-cell agreement of one observation over its attached links (AGREEMENT_MEANING)."""
    if not o.attached:
        return NOT_ATTACHED
    cells = {cr.norm_cell(x) for i in o.attached for x in stated.get((links[i].element, links[i].gene), ())}
    if not cells:
        return NO_STATED
    return MATCH if cr.norm_cell(o.pair.cell) in cells else OTHER_ONLY


# --- summaries ------------------------------------------------------------------------------------
def summary(values: Sequence[float], fraction: bool = True) -> dict[str, Any]:
    """A distribution's summary: quantiles and, for a fraction, counts in tenths with exactly 1 apart."""
    v = np.asarray([x for x in values if x is not None], dtype=float)
    if not len(v):
        return {"n": 0}
    qs = np.quantile(v, [0.1, 0.25, 0.5, 0.75, 0.9])
    nd = 4 if fraction else 1
    out: dict[str, Any] = {
        "n": len(v),
        "mean": round(float(v.mean()), nd),
        "min": round(float(v.min()), nd),
        "p10": round(float(qs[0]), nd),
        "p25": round(float(qs[1]), nd),
        "median": round(float(qs[2]), nd),
        "p75": round(float(qs[3]), nd),
        "p90": round(float(qs[4]), nd),
        "max": round(float(v.max()), nd),
    }
    if fraction:
        if (v < 0).any() or (v > 1).any():
            raise ValueError("a fraction outside [0, 1]")
        tenths: Counter = Counter()
        for x in v:
            if x == 1.0:
                tenths["1.0"] += 1
            else:
                k = min(int(x * 10), 9)
                tenths[f"{FRACTION_BINS[k]:.1f}-{FRACTION_BINS[k + 1]:.1f}"] += 1
        out["tenths"] = dict(sorted(tenths.items()))
        if sum(tenths.values()) != len(v):
            raise ValueError("the tenths do not sum to n")
    return out


def value_counts(values: Iterable[int]) -> dict[str, int]:
    return {str(k): n for k, n in sorted(Counter(values).items())}


# --- one changed observation, once ----------------------------------------------------------------
def row(
    c: Any, p: Any, registry: Any, links: Sequence[Any], stated: dict[tuple[str, str], frozenset[str]]
) -> dict[str, Any]:
    """One changed observation with its full admitted-candidate set under the policy and its coverage."""
    cls = change_class(c, p)
    if cls is None:
        raise ValueError(f"{p.id} is not changed")
    pair = p.pair
    cand = candidate_intervals(registry, p)
    cur = set(c.candidates)
    u_pol = union_coverage(pair.start, pair.end, ((a, b) for _, a, b in cand))
    u_cur = union_coverage(pair.start, pair.end, ((a, b) for i, a, b in cand if i in cur)) if cur else 0.0
    new_status = (
        UNRESOLVED
        if p.resolution in (pc.AMBIGUOUS_ONE, pc.AMBIGUOUS_SEVERAL)
        else UNIQUE_AMONG_ADMITTED
        if p.resolution == pc.UNIQUE
        else pc.UNATTACHED
    )
    return {
        "id": p.id,
        "study": pair.dataset,
        "cell": pair.cell,
        "split": pair.split,
        "outcome": pair.outcome,
        "tested_width_bp": pair.end - pair.start,
        "change": cls,
        "multiplicity": {CURRENT: len(c.candidates), POLICY: len(p.candidates)},
        "union_coverage": {CURRENT: _r(u_cur), POLICY: _r(u_pol)},
        "candidates": [
            {
                "id": i,
                "fraction_covered": _r(candidate_fraction(pair.start, pair.end, a, b)),
                "current": i in cur,
            }
            for i, a, b in cand
        ],
        "attached": [
            {
                "element": links[i].element,
                "gene": links[i].gene,
                "action": links[i].action,
                "stated_cells": sorted(stated.get((links[i].element, links[i].gene), ())),
                "current": i in set(c.attached),
            }
            for i in p.attached
        ],
        "resolution": {CURRENT: c.resolution, POLICY: p.resolution},
        "verdict_current_rule": c.verdict,
        "status_under_new_rule": new_status,
        "direction_description_policy": DESCRIPTION[p.verdict],
        "stated_cell_agreement": {CURRENT: agreement(c, links, stated), POLICY: agreement(p, links, stated)},
    }


def _r(x: float | None) -> float | None:
    return None if x is None else round(x, 4)


def changed_rows(
    cur: Sequence[Any],
    pol: Sequence[Any],
    registry: Any,
    links: Sequence[Any],
    stated: dict[tuple[str, str], frozenset[str]],
) -> list[dict[str, Any]]:
    """One row per changed observation, in the universe's order; an id appearing twice is refused."""
    if [o.id for o in cur] != [o.id for o in pol]:
        raise ValueError("the two rules were not applied to the same observations")
    rows = [row(c, p, registry, links, stated) for c, p in zip(cur, pol, strict=True) if is_changed(c, p)]
    ids = [r["id"] for r in rows]
    if len(set(ids)) != len(ids):
        raise ValueError("an observation appears twice in the changed set")
    return rows


# --- the discovery conditions ---------------------------------------------------------------------
def conditions(
    cur: Sequence[Any], pol: Sequence[Any], rows: Sequence[dict[str, Any]], n_links: int, c4_links: int
) -> dict[str, Any]:
    lost_c = sum(1 for c, p in zip(cur, pol, strict=True) if not set(c.candidates) <= set(p.candidates))
    lost_a = sum(1 for c, p in zip(cur, pol, strict=True) if not set(c.attached) <= set(p.attached))
    same_obs = [o.id for o in cur] == [o.id for o in pol]
    changed = {c.id for c, p in zip(cur, pol, strict=True) if is_changed(c, p)}
    ids = [r["id"] for r in rows]
    once = len(ids) == len(set(ids)) and set(ids) == changed
    full = all(
        len(r["candidates"]) == r["multiplicity"][POLICY]
        and r["multiplicity"][POLICY] > 0
        and all(x["fraction_covered"] is not None for x in r["candidates"])
        and r["union_coverage"][POLICY] is not None
        for r in rows
    )
    one = {
        "current_attachments_lost": {"candidate_sets": lost_c, "attached_links": lost_a},
        "denominators": {
            "observations": {CURRENT: len(cur), POLICY: len(pol), "same_observation_ids": same_obs},
            "predicted_links": n_links,
            "c4_committed_predicted_links": c4_links,
        },
    }
    one["holds"] = lost_c == 0 and lost_a == 0 and same_obs and len(cur) == len(pol) and n_links == c4_links
    two = {
        "rows": len(rows),
        "changed_observations": len(changed),
        "each_once": once,
        "full_candidate_set_and_coverage_on_every_row": full,
    }
    two["holds"] = once and full
    return {
        "1_nothing_lost_denominators_unchanged": one,
        "2_each_observation_once_with_its_candidates_and_coverage": two,
        "hold": one["holds"] and two["holds"],
    }


# --- the shifted screen ---------------------------------------------------------------------------
def stratum_key(p: ms.CrispriPair) -> str:
    return f"{p.dataset}|{p.cell}"


def flags(cur: Sequence[Any], pol: Sequence[Any]) -> list[dict[str, bool]]:
    """Per observation: the screen's quantity and the two quantities reported beside it."""
    return [
        {
            SCREEN_QUANTITY: set(c.attached) != set(p.attached),
            "attached_under_the_policy": bool(p.attached),
            "in_the_changed_set": is_changed(c, p),
        }
        for c, p in zip(cur, pol, strict=True)
    ]


def shifted(pairs: Sequence[ms.CrispriPair], shift: int) -> list[ms.CrispriPair]:
    return [replace(p, start=p.start + shift, end=p.end + shift) for p in pairs]


def enrichment(real: int, shifted_counts: Sequence[int], n: int) -> dict[str, Any]:
    """One stratum's screen: the real and shifted counts over the same n, their rates and ratio."""
    sh = sum(shifted_counts) / len(shifted_counts)
    ratio = real / sh if sh > 0 else None
    if real == 0:
        verdict = NOTHING
    elif sh == 0:
        verdict = UNDEFINED
    else:
        verdict = PASS if ratio >= SCREEN_THRESHOLD else FAIL
    return {
        "observations": n,
        "real": real,
        "shifted": list(shifted_counts),
        "shifted_mean": sh,
        "real_rate": round(real / n, 6) if n else None,
        "shifted_rate": round(sh / n, 6) if n else None,
        "enrichment": None if ratio is None else round(ratio, 3),
        "screen": verdict,
    }


def screen(
    pairs: Sequence[ms.CrispriPair],
    real: Sequence[dict[str, bool]],
    moved: Sequence[Sequence[dict[str, bool]]],
) -> dict[str, Any]:
    """The screen per study and cell stratum and pooled, for the screen's quantity and the two beside it."""
    groups: dict[str, list[int]] = defaultdict(list)
    for i, p in enumerate(pairs):
        groups[stratum_key(p)].append(i)
    groups_all = {"all_observations": list(range(len(pairs))), **dict(sorted(groups.items()))}
    out: dict[str, Any] = {}
    for q in (SCREEN_QUANTITY, *ALSO_SCREENED):
        rows = {}
        for k, idx in groups_all.items():
            rows[k] = enrichment(
                sum(real[i][q] for i in idx), [sum(m[i][q] for i in idx) for m in moved], len(idx)
            )
        if q != SCREEN_QUANTITY:  # reported beside the screen, without a threshold
            for v in rows.values():
                del v["screen"]
        out[q] = rows
    if sum(len(v) for k, v in groups_all.items() if k != "all_observations") != len(pairs):
        raise ValueError("the strata do not sum to the observations")
    return out


def reading(conds: dict[str, Any], scr: dict[str, Any]) -> dict[str, Any]:
    rows = scr[SCREEN_QUANTITY]
    pooled = rows["all_observations"]["screen"]
    strata = {k: v["screen"] for k, v in rows.items() if k != "all_observations"}
    failing = sorted(k for k, v in strata.items() if v == FAIL)
    if not conds["hold"] or pooled != PASS:
        r = DOES_NOT
    elif failing:
        r = RESTRICTED
    else:
        r = MERITS
    return {
        "reading": r,
        "meaning": REGISTRATION["reading"][r],
        "conditions_hold": conds["hold"],
        "pooled_screen": pooled,
        "strata_passing": sorted(k for k, v in strata.items() if v == PASS),
        "strata_failing": failing,
        "strata_undefined_shifted_zero": sorted(k for k, v in strata.items() if v == UNDEFINED),
        "strata_nothing_discovered": sorted(k for k, v in strata.items() if v == NOTHING),
    }


# --- the descriptions -----------------------------------------------------------------------------
def crosstab(rows: Iterable[dict[str, Any]], a: Any, b: Any) -> dict[str, dict[str, int]]:
    t: dict[str, Counter] = defaultdict(Counter)
    for r in rows:
        t[a(r)][b(r)] += 1
    return {k: dict(sorted(v.items())) for k, v in sorted(t.items())}


def descriptions(rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """The changed set described: change class, direction, outcome, stated-cell agreement, study, cell."""
    disc = [r for r in rows if r["change"] in DISCOVERED]

    def ch(r):
        return r["change"]

    def dr(r):
        return r["direction_description_policy"]

    def oc(r):
        return r["outcome"]

    def ag(r):
        return r["stated_cell_agreement"][POLICY]

    def res(r):
        return "unique" if r["resolution"][POLICY] == pc.UNIQUE else "ambiguous"

    return {
        "note": "descriptions only; supports against contradictions authorise nothing, and verdicts stay "
        "under the current rule",
        "changed_set_by_class": {k: sum(r["change"] == k for r in rows) for k in CHANGES},
        "by_class_and_direction": crosstab(rows, ch, dr),
        "by_class_and_outcome": crosstab(rows, ch, oc),
        "by_class_and_stated_cell_agreement": crosstab(rows, ch, ag),
        "discovered_evidence": {
            "n": len(disc),
            "by_resolution": dict(sorted(Counter(res(r) for r in disc).items())),
            "by_direction": dict(sorted(Counter(dr(r) for r in disc).items())),
            "by_outcome_and_direction": crosstab(disc, oc, dr),
            "by_stated_cell_agreement_and_direction": crosstab(disc, ag, dr),
            "by_resolution_and_direction": crosstab(disc, res, dr),
            "by_study_and_direction": crosstab(disc, lambda r: r["study"], dr),
            "by_cell_and_direction": crosstab(disc, lambda r: r["cell"], dr),
            "by_cell_and_stated_cell_agreement": crosstab(disc, lambda r: r["cell"], ag),
            "increases": sum(r["outcome"] == ms.INCREASE for r in disc),
            "well_powered_nulls": sum(r["outcome"] == ms.NULL_INFORMATIVE for r in disc),
            "underpowered_nulls": sum(r["outcome"] == ms.NULL_INCONCLUSIVE for r in disc),
        },
        "by_study_and_class": crosstab(rows, lambda r: r["study"], ch),
        "by_cell_and_class": crosstab(rows, lambda r: r["cell"], ch),
        "by_study_cell_and_class": crosstab(rows, lambda r: f"{r['study']}|{r['cell']}", ch),
    }


def coverage_summaries(rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """Union coverage, candidate fractions, multiplicity and tested width, as distributions."""
    groups: dict[str, list[dict[str, Any]]] = {"all_changed": list(rows)}
    for k in CHANGES:
        groups[k] = [r for r in rows if r["change"] == k]
    disc = [r for r in rows if r["change"] in DISCOVERED]
    groups["discovered_unique"] = [r for r in disc if r["resolution"][POLICY] == pc.UNIQUE]
    groups["discovered_ambiguous"] = [r for r in disc if r["resolution"][POLICY] != pc.UNIQUE]
    out = {}
    for k, g in groups.items():
        out[k] = {
            "observations": len(g),
            "union_coverage_policy": summary([r["union_coverage"][POLICY] for r in g]),
            "union_coverage_current": summary([r["union_coverage"][CURRENT] for r in g]),
            "fraction_of_each_candidate_covered": summary(
                [x["fraction_covered"] for r in g for x in r["candidates"]]
            ),
            "fraction_of_each_candidate_covered_added_by_the_policy": summary(
                [x["fraction_covered"] for r in g for x in r["candidates"] if not x["current"]]
            ),
            "multiplicity_policy": {
                "summary": summary([r["multiplicity"][POLICY] for r in g], fraction=False),
                "counts": value_counts(r["multiplicity"][POLICY] for r in g),
            },
            "multiplicity_current": value_counts(r["multiplicity"][CURRENT] for r in g),
            "tested_width_bp": summary([r["tested_width_bp"] for r in g], fraction=False),
        }
    return out


def newly_ambiguous(rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """Previously unique assignments that become ambiguous: the historical verdict beside 'unresolved'."""
    each = [
        {
            "id": r["id"],
            "historical_verdict_current_rule": r["verdict_current_rule"],
            "status_under_new_rule": r["status_under_new_rule"],
            "resolution": r["resolution"],
            "outcome": r["outcome"],
            "multiplicity": r["multiplicity"],
            "union_coverage": r["union_coverage"],
            "attached_links": {CURRENT: sum(x["current"] for x in r["attached"]), POLICY: len(r["attached"])},
        }
        for r in rows
        if r["resolution"][CURRENT] == pc.UNIQUE
        and r["resolution"][POLICY] in (pc.AMBIGUOUS_ONE, pc.AMBIGUOUS_SEVERAL)
    ]
    if any(e["status_under_new_rule"] != UNRESOLVED for e in each):
        raise ValueError("a newly ambiguous assignment is not marked unresolved")
    return {
        "n": len(each),
        "historical_verdicts": dict(
            sorted(Counter(e["historical_verdict_current_rule"] for e in each).items())
        ),
        "each": each,
    }


def stated_cell_record(
    cur: Sequence[Any],
    pol: Sequence[Any],
    links: Sequence[Any],
    stated: dict[tuple[str, str], frozenset[str]],
    v2_split: dict[str, int] | None,
) -> dict[str, Any]:
    """The 24/13 split, the 39 of 79 and the 39/54 check, recomputed from the census and the stated cells."""
    newly = [
        p for c, p in zip(cur, pol, strict=True) if p.verdict == pc.SUPPORTED and c.verdict != pc.SUPPORTED
    ]
    uniq = [p for p in newly if p.resolution == pc.UNIQUE]
    amb = [p for p in newly if p.resolution != pc.UNIQUE]
    cur_sup = [c for c in cur if c.verdict == pc.SUPPORTED]

    def split(obs):
        t = Counter(agreement(o, links, stated) for o in obs)
        return {a: t.get(a, 0) for a in AGREEMENTS if a != NOT_ATTACHED}

    u, a, s = split(uniq), split(amb), split(cur_sup)
    got = {
        "newly_supported_unique": {
            "n": len(uniq),
            **u,
            "cells_of_the_matches": dict(
                sorted(Counter(o.pair.cell for o in uniq if agreement(o, links, stated) == MATCH).items())
            ),
        },
        "newly_supported_ambiguous": {"n": len(amb), "with_a_stated_cell_match": a[MATCH], **a},
        "current_rule_supported": {"n": len(cur_sup), **s, "not_a_match": s[OTHER_ONLY] + s[NO_STATED]},
    }
    cells = got["newly_supported_unique"]["cells_of_the_matches"]
    reproduces = (
        got["newly_supported_unique"]["n"] == PROVISIONAL["newly_supported_unique"]["n"]
        and u[MATCH] == PROVISIONAL["newly_supported_unique"][MATCH]
        and u[OTHER_ONLY] == PROVISIONAL["newly_supported_unique"][OTHER_ONLY]
        and set(cells) <= {PROVISIONAL["newly_supported_unique"]["matches_all_in"]}
        and len(amb) == PROVISIONAL["newly_supported_ambiguous"]["n"]
        and a[MATCH] == PROVISIONAL["newly_supported_ambiguous"]["with_a_stated_cell_match"]
        and len(cur_sup) == PROVISIONAL["current_rule_supported"]["n"]
        and s[MATCH] == PROVISIONAL["current_rule_supported"][MATCH]
        and got["current_rule_supported"]["not_a_match"]
        == PROVISIONAL["current_rule_supported"]["not_a_match"]
    )
    check = None
    if v2_split is not None:
        check = {
            "v2_activity_supported_in_context": v2_split["in_context"],
            "v2_activity_supported_elsewhere_only": v2_split["elsewhere_only"],
            "equal": s[MATCH] == v2_split["in_context"]
            and got["current_rule_supported"]["not_a_match"] == v2_split["elsewhere_only"],
        }
    return {
        "scope": REGISTRATION["stated_cell_record"]["scope"],
        "recomputed": got,
        "provisional_in_memory_reading": PROVISIONAL,
        "reproduces_the_provisional_reading": reproduces,
        "check_against_the_v2_v3_direction_split": check,
    }


# --- the manifest ---------------------------------------------------------------------------------
@mf.depends_on_models("alphagenome")  # it reads the compiled predicted layer
def manifest(programs: list[Path]) -> dict[str, Any]:
    inputs = [mf.input_entry(p, partition=None) for p in programs]
    for name in ms.CRISPRI_FILES:
        p = ms.CRISPRI_KNOWLEDGE / name
        if p.exists():
            inputs.append(mf.input_entry(p, partition=ms.CRISPRI_SPLIT_OF[name]))
    inputs += [
        mf.input_entry(p, partition=None)
        for c in ho.CHROMS
        if (p := RESULTS_DIR / f"ccres_{c}.bed.gz").exists()
    ]
    for name in (
        "c4_alphagenome_ablation",
        "placement_audit",
        "placement_census",
        "attribution_correctness_v3",
    ):
        inputs.append(mf.input_entry(RESULTS_DIR / f"{name}.json", partition=None))
    return {
        "sources": [
            {
                "accession": "EngreitzLab/CRISPR_comparison EPCrisprBenchmark (Gschwind et al. 2025)",
                "version": "main, as fetched 2026-09-16; pinned by sha256",
            },
            {
                "accession": "ENCODE SCREEN cCREs v3, https://downloads.wenglab.org/V3/GRCh38-cCREs.bed",
                "version": "per-chromosome subsets as distilled; pinned by sha256",
            },
            {
                "accession": "this repository: the 24 compiled programs "
                "data/knowledge/compiled/noncoding_<chrom>.bio (generated, untracked)",
                "version": "pinned by sha256",
            },
            {"accession": "this repository: data/results/c4_alphagenome_ablation.json", "version": "b7e4bf0"},
            {"accession": "this repository: data/results/placement_audit.json", "version": "6c0d39f"},
            {"accession": "this repository: data/results/placement_census.json", "version": "32af2cb"},
            {
                "accession": "this repository: data/results/attribution_correctness_v3.json",
                "version": "d1e09ae",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "current_rule_reciprocal_overlap": ms.RECIPROCAL_OVERLAP,
            "policy": POLICY,
            "policy_threshold": pa.POLICY_THRESHOLD,
            "shift_bp": list(SHIFTS),
            "screen_threshold": SCREEN_THRESHOLD,
            "power_for_negatives": ms.POWER_FOR_NEGATIVES,
            "well_powered": ms.WELL_POWERED,
            "registered": REGISTERED,
        },
        "exclusions": [
            "CRISPRi pairs the benchmark marks as not a valid connection are not read (none in these files)",
            "the per-element response cache is not opened; the model is not called; nothing is downloaded",
            "no other endpoint is read",
        ],
        "partitions": {
            ms.TRAINING: "the CRISPRi benchmark's training file (K562); read as measurement, nothing fitted",
            ms.HELDOUT: "the CRISPRi benchmark's held-out file; read as measurement, nothing fitted",
        },
    }


# --- the run --------------------------------------------------------------------------------------
def _gate_census(obs: dict[str, list[Any]], committed: dict[str, Any]) -> dict[str, Any]:
    for name in (CURRENT, POLICY):
        mine = json.loads(json.dumps(pc.observation_table(obs[name])))
        if mine != committed["rules"][name]["observations"]:
            raise SystemExit(f"the census's {name} observation table does not reproduce")
    trans = dict(
        sorted(
            Counter(
                f"{c.verdict}->{p.verdict}" for c, p in zip(obs[CURRENT], obs[POLICY], strict=True)
            ).items()
        )
    )
    if trans != committed["difference_policy_minus_current"]["verdict_transitions"]:
        raise SystemExit("the census's verdict transitions do not reproduce")
    return {
        "observation_tables_equal_to_committed": [CURRENT, POLICY],
        "verdict_transitions_equal_to_committed": True,
    }


def main() -> None:
    t0, r0 = time.perf_counter(), resource.getrusage(resource.RUSAGE_SELF)
    programs = ho.compiled_programs()
    pairs, invalid = ms.load_crispri()
    links = pc.predicted_links(programs)
    if len({x.element for x in links}) != len(links):
        raise SystemExit("a compiled element id repeats")
    reg_items, _ = pa.registry_elements()
    where = {i: (c, s, e) for c, s, e, i in reg_items}
    off = sum(1 for x in links if where.get(x.element) != (x.chrom, x.start, x.end))
    if off:
        raise SystemExit(f"{off} compiled elements are not registry elements with the registry's coordinates")
    del where
    registry = pa.ElementSet(reg_items)
    link_set = pa.ElementSet((x.chrom, x.start, x.end, x.index) for x in links)
    c4 = json.loads((RESULTS_DIR / "c4_alphagenome_ablation.json").read_text())["elements"]
    census = json.loads((RESULTS_DIR / "placement_census.json").read_text())
    v3 = json.loads((RESULTS_DIR / "attribution_correctness_v3.json").read_text())
    print("observations", len(pairs), "links", len(links), "registry", len(registry), flush=True)

    obs = {name: pc.observe(pairs, links, link_set, registry, rule) for name, rule in pc.RULES.items()}
    gate = _gate_census(obs, census)
    cur, pol = obs[CURRENT], obs[POLICY]

    claims = cr.compiled_claims()
    stated, ruled = stated_cells(claims)
    attached_keys = {(links[i].element, links[i].gene) for o in pol for i in o.attached}
    rows = changed_rows(cur, pol, registry, links, stated)
    conds = conditions(cur, pol, rows, len(links), c4["predicted_elements"])
    print("changed", len(rows), "conditions hold", conds["hold"], flush=True)

    real = flags(cur, pol)
    moved = []
    for shift in SHIFTS:
        sp = shifted(pairs, shift)
        moved.append(
            flags(*(pc.observe(sp, links, link_set, registry, pc.RULES[n]) for n in (CURRENT, POLICY)))
        )
    scr = screen(pairs, real, moved)
    read = reading(conds, scr)
    print("screen", scr[SCREEN_QUANTITY]["all_observations"], read["reading"], flush=True)

    ic = v3["in_context_side_by_side"]["v2"]["activity"]
    v2_split = {"in_context": ic["supported_in_context"], "elsewhere_only": ic["supported_elsewhere_only"]}
    record = stated_cell_record(cur, pol, links, stated, v2_split)

    r1 = resource.getrusage(resource.RUSAGE_SELF)
    peak = r1.ru_maxrss if sys.platform == "darwin" else r1.ru_maxrss * 1024
    payload = {
        "question": "which screened CRISPRi observations audit A's policy would make discoverable as a "
        "set-valued attachment beside C4's current rule, with union coverage, multiplicity, a shifted "
        "screen per study and cell, and verdicts kept under the current rule (Decision 1, discovery only)",
        "status": "internal development evidence: A's adoption review, discovery only, bounded, read-only "
        "with respect to production, recorded apart, one pass. Every outcome read here was read before "
        "(prior_exposure). It proposes no implementation, changes no matcher and regenerates nothing",
        "registration": registration(),
        "alphagenome_requests": 0,
        "downloads": 0,
        "per_element_response_cache_opened": False,
        "census_reproduced": gate,
        "universe": {
            "observations": len(pairs),
            "invalid_rows_not_read": invalid,
            "predicted_links": len(links),
            "changed_set": len(rows),
            "unchanged": len(pairs) - len(rows),
        },
        "discovery_conditions": conds,
        "screen": scr,
        "decision_1_reading": read,
        "verdict_eligibility": REGISTRATION["verdict_eligibility"],
        "coverage": coverage_summaries(rows),
        "descriptions": descriptions(rows),
        "newly_ambiguous": newly_ambiguous(rows),
        "stated_cell_record": record,
        "stated_cell_join": {
            "attached_links_under_the_policy": len(attached_keys),
            "with_a_compiled_rule": len(attached_keys & ruled),
            "with_a_stated_cell": sum(1 for k in attached_keys if stated.get(k)),
        },
        "changed_observations": rows,
        "prior_exposure": PRIOR_EXPOSURE,
        "cost": {
            "wall_seconds": round(time.perf_counter() - t0, 1),
            "cpu_seconds": round((r1.ru_utime - r0.ru_utime) + (r1.ru_stime - r0.ru_stime), 1),
            "peak_memory_mb": round(peak / 2**20, 1),
            "downloads": 0,
            "requests": 0,
        },
    }
    p = save_result(NAME, payload, manifest=manifest(programs))
    print("saved", p, flush=True)


if __name__ == "__main__":
    main()
