# SPDX-License-Identifier: AGPL-3.0-or-later
"""The placement census: C4's CRISPRi evidence attachment under two interval policies, in one pass.

    uv run python scripts/placement_census.py

Lane-census, 2026-09-29, for the coordinator (genomeos-9f), approved by the external reviewer as a bounded,
read-only census recorded apart, that stops after one pass. Item 13 C4 (`b7e4bf0`) attaches CRISPRi
evidence to the 440,377 predicted links of the 24 compiled programs under the measured layer's reciprocal
0.5 overlap rule (`ablation.CrispriIndex.of`, `ablation.link_class`). Audit A (`dd8c49c`, `6c0d39f`)
found that rule fails most decrease links on width and registered one other policy before scoring it:
an overlap of at least half of the smaller interval. This census counts C4's attachment under both, with
each screened CRISPRi observation kept as one set-valued observation.

It replaces no headline result and changes no production matching: the second policy lives only in this
script (through placement_audit.policy_rule), and ablation.py, holdout.py and measured.py are imported,
never edited. It adopts nothing and regenerates nothing. Every outcome it reads was read before
(PRIOR_EXPOSURE), so it is internal development evidence. No model request is made and nothing is
downloaded.

The census definition (REGISTRATION) was committed with this file before any count of it was taken.
Writes data/results/placement_census.json through the result contract.
"""

from __future__ import annotations

import importlib.util
import json
import resource
import sys
import time
from collections import Counter, defaultdict
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from genomeos import manifest as mf
from genomeos.attribution import ablation as ab
from genomeos.attribution import holdout as ho
from genomeos.attribution import measured as ms
from genomeos.results import RESULTS_DIR, save_result


def _audit_a() -> Any:
    """Audit A's module (scripts/placement_audit.py), imported and never edited: its two rules, its
    element index and its registry reader are the ones this census uses."""
    if "placement_audit" in sys.modules:
        return sys.modules["placement_audit"]
    spec = importlib.util.spec_from_file_location(
        "placement_audit", Path(__file__).resolve().parent / "placement_audit.py"
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


pa = _audit_a()

NAME = "placement_census"
REGISTERED = "2026-09-29"

# ==================================================================================================
# The registration (2026-09-29, lane-census), committed before any count of this census was taken.
# ==================================================================================================
CURRENT = "reciprocal_half"
POLICY = "min_side_half"
RULES: dict[str, Callable[[int, int, int, int], bool]] = {
    CURRENT: pa.current_rule,
    POLICY: pa.policy_rule,
}
RULE_MEANING = {
    CURRENT: f"C4's rule, unchanged: measured.measures, the overlap at least {ms.RECIPROCAL_OVERLAP} of "
    "BOTH widths (ablation.CrispriIndex.of)",
    POLICY: "audit A's registered policy (dd8c49c), unchanged: the overlap at least "
    f"{pa.POLICY_THRESHOLD} of the SMALLER width (placement_audit.policy_rule). Applied here only",
}

# the observation and its candidate set
UNATTACHED = "unattached"
UNIQUE = "unique_element"
AMBIGUOUS_ONE = "ambiguous_one_link"
AMBIGUOUS_SEVERAL = "ambiguous_several_links"
RESOLUTIONS = (UNIQUE, AMBIGUOUS_ONE, AMBIGUOUS_SEVERAL, UNATTACHED)
RESOLUTION_MEANING = {
    UNIQUE: "exactly one registry element meets the rule with the tested interval, and it carries a "
    "predicted link to the measured gene: a unique-element assignment",
    AMBIGUOUS_ONE: "one predicted link to the measured gene is attached, but more than one registry element "
    "meets the rule with the tested interval, so the response could be explained by an element the "
    "link does not name: an ambiguous candidate set",
    AMBIGUOUS_SEVERAL: "more than one predicted link to the measured gene is attached: an ambiguous "
    "candidate set, one observation over several links",
    UNATTACHED: "no predicted link to the measured gene meets the rule with the tested interval: the "
    "observation assesses no predicted link",
}

# the observation's verdict against its attached links
SUPPORTED = "supported"
REFUTED_NULL = "refuted_well_powered_null"
REFUTED_OPPOSITE = "refuted_opposite_sign"
SPLIT = "split_direction"
UNASSESSED_UNDERPOWERED = "unassessed_underpowered"
UNASSESSED_MISSING = "unassessed_missing"
VERDICTS = (
    SUPPORTED,
    REFUTED_NULL,
    REFUTED_OPPOSITE,
    SPLIT,
    UNASSESSED_UNDERPOWERED,
    UNASSESSED_MISSING,
    UNATTACHED,
)
WANT = {"activates": ms.DECREASE, "inhibits": ms.INCREASE}
VERDICT_MEANING = {
    SUPPORTED: "attached, and a significant change in the direction every attached link predicts "
    "(activates: a decrease; inhibits: an increase)",
    REFUTED_NULL: "attached, and a well-powered null (PowerAtEffectSize20 >= 0.8, measured.outcome)",
    REFUTED_OPPOSITE: "attached, and a significant change opposite to the direction every attached link "
    "predicts",
    SPLIT: "attached to links of both directions, and a significant change: it supports some and "
    "contradicts others, counted as neither",
    UNASSESSED_UNDERPOWERED: "attached, and an underpowered null",
    UNASSESSED_MISSING: "attached, and no effect size in the table",
    UNATTACHED: "no attached link (RESOLUTION_MEANING)",
}

# the strata
WIDTH_EDGE = 700  # the widest tested interval the reciprocal rule can pair with a 350 bp registry element
WIDTH_FINE = (
    ("zero_width", None, 0),
    ("1_to_74", 1, 74),
    ("75_to_350", 75, 350),
    ("351_to_500", 351, 500),
    ("501_to_700", 501, 700),
    ("701_to_1000", 701, 1000),
    ("1001_to_2000", 1001, 2000),
    ("over_2000", 2001, None),
)
STRATA = ("study", "cell", "split", "width", "width_fine")

REGISTRATION = {
    "registered": REGISTERED,
    "approved_as": "a bounded, read-only census, recorded apart, that stops after one pass (the external "
    "reviewer, relayed by the coordinator genomeos-9f)",
    "rules": RULE_MEANING,
    "rules_note": "the current rule implies the policy (half of both widths is at least half of the "
    "smaller), so every attachment under the current rule is one under the policy; this is checked on "
    "every observation, not assumed",
    "observation": "one row of the CRISPRi benchmark (both files, ValidConnection TRUE, "
    "measured.parse_crispri): its tested interval, measured gene, cell and file. Its id is "
    "split|chrom:start-end|gene|cell, checked unique; its outcome is measured's (decrease, increase, "
    "well-powered null, underpowered null, missing)",
    "predicted_link": "one predicted element block of the 24 compiled programs (not `_measured`), with its "
    "first target and its action, exactly as ablation.program_census reads it",
    "candidate_set": "the ENCODE cCRE v3 registry elements that meet the rule with the tested interval. "
    "Every compiled predicted element is a registry element with the registry's id and coordinates "
    "(checked), so the attached links are a subset of the candidate set",
    "attached_links": "the predicted links whose element meets the rule with the tested interval and whose "
    "target is the measured gene, in any cell: C4's element-label attachment (ablation.link_class with "
    "cell None)",
    "set_valued_rule": "one observation is one count. An observation attached to several links, or whose "
    "candidate set holds several elements, is one ambiguous observation, never several validations. Counts "
    "of observations and counts of element-level attachments (observation-link pairs, and C4's per-link "
    "census) are reported apart, and every supported count of the census is a count of observations",
    "resolution": RESOLUTION_MEANING,
    "verdict": VERDICT_MEANING,
    "per_rule": [
        "observations by verdict, and by verdict and resolution: supported, refuted (well-powered null and "
        "opposite sign apart), split direction, unassessed (underpowered and missing apart), unattached",
        "observations by outcome and verdict, so decreases, increases, well-powered nulls and underpowered "
        "nulls are each read apart",
        "element-level attachments apart: observation-link pairs, links with an attached observation, and "
        "C4's per-link census (ablation.link_class) with the distinct observations behind its supported "
        "links and how many of those links rest only on ambiguous observations",
        "the same observation table per stratum: study (the Dataset column), cell, split, tested width "
        f"(<= {WIDTH_EDGE} bp, > {WIDTH_EDGE} bp, and a zero-width class) and finer width bins "
        f"({', '.join(n for n, _, _ in WIDTH_FINE)})",
    ],
    "difference": [
        "observations newly attached under the policy, by resolution, verdict, width, study and cell",
        f"new observation-link attachments, and how many come from tested intervals over {WIDTH_EDGE} bp",
        "observations attached under both rules whose candidate set or attached links grew, and resolution "
        "transitions",
        "the verdict transition table (current to policy) and C4's per-link class transition table",
    ],
    "denominator": "unchanged under both rules and checked: the same 440,377 predicted links (C4's "
    "committed count) and the same screened observations (every valid benchmark row); every table sums "
    "to its denominator",
    "gates": [
        "C4's committed per-link census (link_by_class and link_by_detail) reproduces exactly under the "
        "current rule, both through ablation.CrispriIndex.of and through this census's own index",
        "audit A's committed placements reproduce: observations with a compiled element meeting the rule, "
        "and with one only in the registry, under both rules; and the policy's multi-element placements",
    ],
    "cost": "wall time, CPU time and peak memory of the run",
    "out_of_scope": [
        "C4's rule family (R1, the rule's own cell) and its headline family: not censused",
        "any other endpoint (lentiMPRA, VISTA, saturation mutagenesis, GTEx): not read",
        "any recommendation to adopt the policy: the census stops after one pass, changes no matcher and "
        "regenerates nothing; adopting a matching change waits until the evidence-resolution implications "
        "are reviewed",
    ],
}
PRIOR_EXPOSURE = [
    *pa.PRIOR_EXPOSURE,
    {"commit": "dd8c49c", "read": "audit A's registration of the policy (2026-09-29)"},
    {
        "commit": "6c0d39f",
        "read": "audit A's result: every screened pair placed under both rules, positives and the "
        "non-decrease comparator apart, with the policy's multi-element placements (2026-09-29)",
    },
    {"commit": "0183ef5", "read": "audit A's per-width table of links, placements and rescues (2026-09-29)"},
    {"commit": "44d95d3", "read": "audit A's checkpoint (2026-09-29)"},
    {"commit": "8f87bf6", "read": "the reviewer's correction beside audit A (2026-09-29)"},
]


def registration() -> dict[str, Any]:
    return REGISTRATION


# --- the inputs -----------------------------------------------------------------------------------
@dataclass(frozen=True)
class PredictedLink:
    """One predicted element block and its link, as ablation.program_census reads it."""

    index: int
    chrom: str
    start: int
    end: int
    element: str
    gene: str
    action: str


def predicted_links(programs: Iterable[Path]) -> list[PredictedLink]:
    out: list[PredictedLink] = []
    for path in programs:
        for kind, name, f in ho.program_blocks(path):
            if kind != "element" or name.endswith("_measured"):
                continue
            m = ho.LOCUS_RE.match(f.get("locus", ""))
            if not m:
                continue
            gene = f.get("targets", "").split(",")[0].strip()
            action = "activates" if "activates_target" in f.get("activity", "") else "inhibits"
            out.append(
                PredictedLink(len(out), m.group(1), int(m.group(2)), int(m.group(3)), name, gene, action)
            )
    return out


def observation_id(p: ms.CrispriPair) -> str:
    return f"{p.split}|{p.chrom}:{p.start}-{p.end}|{p.gene}|{p.cell}"


def width_class(width: int) -> str:
    if width <= 0:
        return "zero_width"
    return f"le_{WIDTH_EDGE}" if width <= WIDTH_EDGE else f"gt_{WIDTH_EDGE}"


def width_fine(width: int) -> str:
    for name, lo, hi in WIDTH_FINE:
        if (lo is None or width >= lo) and (hi is None or width <= hi):
            return name
    raise ValueError(width)


def stratum(p: ms.CrispriPair, kind: str) -> str:
    if kind == "study":
        return p.dataset
    if kind == "cell":
        return p.cell
    if kind == "split":
        return p.split
    if kind == "width":
        return width_class(p.end - p.start)
    if kind == "width_fine":
        return width_fine(p.end - p.start)
    raise ValueError(kind)


# --- one observation under one rule ---------------------------------------------------------------
@dataclass(frozen=True)
class Observation:
    """One screened observation under one rule: its candidate set, its attached links, its reading."""

    id: str
    pair: ms.CrispriPair
    candidates: tuple[str, ...]  # registry ids meeting the rule
    compiled: tuple[int, ...]  # predicted link indices meeting the rule, any gene
    attached: tuple[int, ...]  # the subset whose target is the measured gene
    resolution: str
    verdict: str


def resolve(n_candidates: int, n_attached: int) -> str:
    if n_attached == 0:
        return UNATTACHED
    if n_attached > 1:
        return AMBIGUOUS_SEVERAL
    return UNIQUE if n_candidates == 1 else AMBIGUOUS_ONE


def verdict(outcome: str, actions: Iterable[str]) -> str:
    """The reading of one attached observation against the actions of its attached links."""
    acts = set(actions)
    if not acts:
        return UNATTACHED
    if outcome in (ms.DECREASE, ms.INCREASE):
        agree = {a for a in acts if WANT[a] == outcome}
        if agree == acts:
            return SUPPORTED
        if not agree:
            return REFUTED_OPPOSITE
        return SPLIT
    if outcome == ms.NULL_INFORMATIVE:
        return REFUTED_NULL
    if outcome == ms.NULL_INCONCLUSIVE:
        return UNASSESSED_UNDERPOWERED
    if outcome == ms.MISSING:
        return UNASSESSED_MISSING
    raise ValueError(outcome)


def observe(
    pairs: list[ms.CrispriPair],
    links: list[PredictedLink],
    link_set: Any,
    registry: Any,
    rule: Callable[[int, int, int, int], bool],
) -> list[Observation]:
    """Every observation under one rule. `link_set` and `registry` are placement_audit.ElementSet indexes
    (the links keyed by their index, the registry by its id)."""
    ids = [observation_id(p) for p in pairs]
    if len(set(ids)) != len(ids):
        raise ValueError("an observation id repeats: one observation would be counted twice")
    out = []
    for oid, p in zip(ids, pairs, strict=True):
        cands = tuple(
            sorted(i for a, b, i in registry.touching(p.chrom, p.start, p.end) if rule(a, b, p.start, p.end))
        )
        comp = tuple(
            sorted(i for a, b, i in link_set.touching(p.chrom, p.start, p.end) if rule(a, b, p.start, p.end))
        )
        att = tuple(i for i in comp if links[i].gene == p.gene)
        missing = {links[i].element for i in att} - set(cands)
        if missing:
            raise ValueError(f"{oid}: attached elements outside the candidate set: {sorted(missing)}")
        out.append(
            Observation(
                id=oid,
                pair=p,
                candidates=cands,
                compiled=comp,
                attached=att,
                resolution=resolve(len(cands), len(att)),
                verdict=verdict(p.outcome, (links[i].action for i in att)),
            )
        )
    return out


# --- the tables -----------------------------------------------------------------------------------
def observation_table(obs: list[Observation]) -> dict[str, Any]:
    """Counts of observations: each observation id enters every table exactly once."""
    if len({o.id for o in obs}) != len(obs):
        raise ValueError("an observation id repeats")
    by_verdict = Counter(o.verdict for o in obs)
    by_res = Counter(o.resolution for o in obs)
    vr: dict[str, Counter] = defaultdict(Counter)
    ov: dict[str, Counter] = defaultdict(Counter)
    for o in obs:
        vr[o.verdict][o.resolution] += 1
        ov[o.pair.outcome][o.verdict] += 1
    attached = sum(1 for o in obs if o.attached)
    out = {
        "observations": len(obs),
        "attached": attached,
        "unattached": len(obs) - attached,
        "by_resolution": {r: by_res.get(r, 0) for r in RESOLUTIONS},
        "by_verdict": {v: by_verdict.get(v, 0) for v in VERDICTS},
        "summary": {
            "supported": by_verdict.get(SUPPORTED, 0),
            "refuted": by_verdict.get(REFUTED_NULL, 0) + by_verdict.get(REFUTED_OPPOSITE, 0),
            "split_direction": by_verdict.get(SPLIT, 0),
            "unassessed_attached": by_verdict.get(UNASSESSED_UNDERPOWERED, 0)
            + by_verdict.get(UNASSESSED_MISSING, 0),
            "unattached": by_verdict.get(UNATTACHED, 0),
        },
        "by_verdict_and_resolution": {
            v: {r: vr[v].get(r, 0) for r in RESOLUTIONS if vr[v].get(r, 0)} for v in VERDICTS if v in vr
        },
        "by_outcome_and_verdict": {k: dict(sorted(ov[k].items())) for k in ms.OUTCOMES if k in ov},
    }
    if sum(out["by_verdict"].values()) != len(obs) or sum(out["by_resolution"].values()) != len(obs):
        raise ValueError("the observation table does not sum to its denominator")
    return out


def strata_tables(obs: list[Observation]) -> dict[str, dict[str, dict[str, Any]]]:
    out: dict[str, dict[str, dict[str, Any]]] = {}
    for kind in STRATA:
        groups: dict[str, list[Observation]] = defaultdict(list)
        for o in obs:
            groups[stratum(o.pair, kind)].append(o)
        out[kind] = {k: _compact(observation_table(v)) for k, v in sorted(groups.items())}
        if sum(t["observations"] for t in out[kind].values()) != len(obs):
            raise ValueError(f"stratum {kind} does not sum to the observations")
    return out


def _compact(t: dict[str, Any]) -> dict[str, Any]:
    """A stratum's row: the counts only, without the outcome breakdown."""
    return {
        "observations": t["observations"],
        "attached": t["attached"],
        **{f"resolution_{r}": n for r, n in t["by_resolution"].items() if r != UNATTACHED},
        **t["by_verdict"],
        "increases_attached": sum(
            n for v, n in t["by_outcome_and_verdict"].get(ms.INCREASE, {}).items() if v != UNATTACHED
        ),
    }


def link_census(
    links: list[PredictedLink], obs: list[Observation]
) -> tuple[Counter, Counter, dict[int, list[Observation]], list[tuple[str, str | None]]]:
    """C4's per-link census under one rule, through ablation.link_class on the observations whose tested
    interval meets the rule with the link's element (any gene), i.e. CrispriIndex.of's pairs. Returns the
    counts by class and by detail, the observations on each link, and each link's (class, reason)."""
    on: dict[int, list[Observation]] = defaultdict(list)
    for o in obs:
        for i in o.compiled:
            on[i].append(o)
    by_class: Counter = Counter()
    by_detail: Counter = Counter()
    per_link: list[tuple[str, str | None]] = []
    for x in links:
        cls, reason, detail = ab.link_class([o.pair for o in on.get(x.index, ())], x.gene, x.action)
        by_class[(cls, reason)] += 1
        by_detail[detail if cls != ab.SURVIVES else "survives"] += 1
        per_link.append((cls, reason))
    return by_class, by_detail, on, per_link


def element_attachments(
    links: list[PredictedLink], obs: list[Observation], on: dict[int, list[Observation]]
) -> dict[str, Any]:
    """Element-level attachments, apart from the observation counts: pairs, links, and the observations
    behind C4's supported links."""
    pairs_n = sum(len(o.attached) for o in obs)
    linked = {i for o in obs for i in o.attached}
    supported_links = 0
    behind: set[str] = set()
    only_ambiguous = 0
    with_unique = 0
    per_link: Counter = Counter()
    for x in links:
        sup = [o for o in on.get(x.index, ()) if o.pair.gene == x.gene and o.pair.outcome == WANT[x.action]]
        if not sup:
            continue
        supported_links += 1
        behind.update(o.id for o in sup)
        per_link[min(len(sup), 3)] += 1
        if any(o.resolution == UNIQUE for o in sup):
            with_unique += 1
        else:
            only_ambiguous += 1
    return {
        "observation_link_pairs": pairs_n,
        "links_with_an_attached_observation": len(linked),
        "supported_links": supported_links,
        "distinct_observations_behind_supported_links": len(behind),
        "supported_links_with_a_unique_element_observation": with_unique,
        "supported_links_resting_only_on_ambiguous_observations": only_ambiguous,
        "supported_links_by_supporting_observations": {
            ("3_or_more" if k == 3 else str(k)): n for k, n in sorted(per_link.items())
        },
    }


def class_rows(c: Counter) -> list[dict[str, Any]]:
    return [{"class": k[0], "reason": k[1], "n": n} for k, n in c.most_common()]


# --- the difference between the rules -------------------------------------------------------------
def difference(
    cur: list[Observation],
    pol: list[Observation],
    link_cur: list[tuple[str, str | None]],
    link_pol: list[tuple[str, str | None]],
) -> dict[str, Any]:
    if [o.id for o in cur] != [o.id for o in pol]:
        raise ValueError("the two rules were not applied to the same observations")
    lost_candidates = sum(
        1 for c, p in zip(cur, pol, strict=True) if not set(c.candidates) <= set(p.candidates)
    )
    lost_attached = sum(1 for c, p in zip(cur, pol, strict=True) if not set(c.attached) <= set(p.attached))
    new = [p for c, p in zip(cur, pol, strict=True) if not c.attached and p.attached]
    both = [(c, p) for c, p in zip(cur, pol, strict=True) if c.attached and p.attached]
    grown = [
        (c, p) for c, p in both if len(p.attached) > len(c.attached) or len(p.candidates) > len(c.candidates)
    ]
    new_pairs = [(p, i) for c, p in zip(cur, pol, strict=True) for i in set(p.attached) - set(c.attached)]
    new_ids = {o.id for o in new}
    wide = f"gt_{WIDTH_EDGE}"

    def tally(rows: list[Observation], key: Callable[[Observation], str]) -> dict[str, int]:
        return dict(sorted(Counter(key(o) for o in rows).items()))

    newly_supported = [
        p for c, p in zip(cur, pol, strict=True) if p.verdict == SUPPORTED and c.verdict != SUPPORTED
    ]
    return {
        "attachments_lost_under_the_policy": {
            "candidate_sets": lost_candidates,
            "attached_links": lost_attached,
        },
        "newly_attached_observations": {
            "n": len(new),
            "by_resolution": {r: sum(o.resolution == r for o in new) for r in RESOLUTIONS if r != UNATTACHED},
            "unique": sum(o.resolution == UNIQUE for o in new),
            "ambiguous": sum(o.resolution in (AMBIGUOUS_ONE, AMBIGUOUS_SEVERAL) for o in new),
            "by_width": tally(new, lambda o: width_class(o.pair.end - o.pair.start)),
            "by_width_fine": tally(new, lambda o: width_fine(o.pair.end - o.pair.start)),
            "by_verdict": tally(new, lambda o: o.verdict),
            "by_outcome": tally(new, lambda o: o.pair.outcome),
            "by_study": tally(new, lambda o: o.pair.dataset),
            "by_cell": tally(new, lambda o: o.pair.cell),
            "by_split": tally(new, lambda o: o.pair.split),
        },
        "newly_supported_observations": {
            "n": len(newly_supported),
            "by_resolution": tally(newly_supported, lambda o: o.resolution),
            "by_width": tally(newly_supported, lambda o: width_class(o.pair.end - o.pair.start)),
            "by_study": tally(newly_supported, lambda o: o.pair.dataset),
        },
        "new_observation_link_attachments": {
            "n": len(new_pairs),
            f"from_tested_intervals_over_{WIDTH_EDGE}_bp": sum(
                width_class(o.pair.end - o.pair.start) == wide for o, _ in new_pairs
            ),
            "on_newly_attached_observations": sum(1 for o, _ in new_pairs if o.id in new_ids),
            "by_width_fine": dict(
                sorted(Counter(width_fine(o.pair.end - o.pair.start) for o, _ in new_pairs).items())
            ),
        },
        "attached_under_both_rules": {
            "n": len(both),
            "candidate_set_or_attached_links_grew": len(grown),
            "resolution_transitions": dict(
                sorted(Counter(f"{c.resolution}->{p.resolution}" for c, p in both).items())
            ),
            "verdict_changed": sum(c.verdict != p.verdict for c, p in both),
        },
        "verdict_transitions": dict(
            sorted(Counter(f"{c.verdict}->{p.verdict}" for c, p in zip(cur, pol, strict=True)).items())
        ),
        "link_class_transitions": dict(
            sorted(
                Counter(
                    f"{a[0]}/{a[1]}->{b[0]}/{b[1]}" for a, b in zip(link_cur, link_pol, strict=True) if a != b
                ).items()
            )
        ),
        "links_whose_class_changed": sum(a != b for a, b in zip(link_cur, link_pol, strict=True)),
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
    inputs.append(mf.input_entry(RESULTS_DIR / "c4_alphagenome_ablation.json", partition=None))
    inputs.append(mf.input_entry(RESULTS_DIR / "placement_audit.json", partition=None))
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
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "current_rule_reciprocal_overlap": ms.RECIPROCAL_OVERLAP,
            "policy": POLICY,
            "policy_threshold": pa.POLICY_THRESHOLD,
            "width_edge_bp": WIDTH_EDGE,
            "power_for_negatives": ms.POWER_FOR_NEGATIVES,
            "well_powered": ms.WELL_POWERED,
            "registered": REGISTERED,
        },
        "exclusions": [
            "CRISPRi pairs the benchmark marks as not a valid connection are not read (none in these files)",
            "the per-element response cache is not opened; the model is not called; nothing is downloaded",
            "no other endpoint is read; C4's rule and headline families are not censused",
        ],
        "partitions": {
            ms.TRAINING: "the CRISPRi benchmark's training file (K562); read as measurement, nothing fitted",
            ms.HELDOUT: "the CRISPRi benchmark's held-out file; read as measurement, reported apart as a "
            "stratum",
        },
    }


# --- the run --------------------------------------------------------------------------------------
def _gate_c4(by_class: Counter, by_detail: Counter, committed: dict[str, Any]) -> None:
    want_class = {(r["class"], r["reason"]): r["n"] for r in committed["link_by_class"]}
    if dict(by_class) != want_class or dict(by_detail) != committed["link_by_detail"]:
        raise SystemExit(f"C4's per-link census does not reproduce: {dict(by_class)} / {dict(by_detail)}")


def _gate_audit_a(obs: list[Observation], committed: dict[str, Any], multi: int | None) -> dict[str, int]:
    placed = sum(1 for o in obs if o.compiled)
    registry_only = sum(1 for o in obs if not o.compiled and o.candidates)
    want_placed = sum(committed["coverage"][g]["placed_compiled"] for g in ("positive", "null"))
    want_reg = sum(committed["coverage"][g]["registry_would"] for g in ("positive", "null"))
    got = {"placed_compiled": placed, "registry_only": registry_only}
    if (placed, registry_only) != (want_placed, want_reg):
        raise SystemExit(f"audit A's placements do not reproduce: {got} against {want_placed}, {want_reg}")
    if multi is not None:
        many = sum(1 for o in obs if len(o.compiled) > 1)
        if many != multi:
            raise SystemExit(f"audit A's multi-element placements do not reproduce: {many} against {multi}")
        got["carried_by_more_than_one_compiled_element"] = many
    return got


def main() -> None:
    t0, r0 = time.perf_counter(), resource.getrusage(resource.RUSAGE_SELF)
    programs = ho.compiled_programs()
    pairs, invalid = ms.load_crispri()
    links = predicted_links(programs)
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
    if len(links) != c4["predicted_elements"]:
        raise SystemExit(f"{len(links)} predicted links against C4's {c4['predicted_elements']}")
    print("observations", len(pairs), "links", len(links), "registry", len(registry), flush=True)

    # gate 1: C4's census through its own index, before anything of this census is read
    crispri = ab.CrispriIndex(pairs)
    direct_class: Counter = Counter()
    direct_detail: Counter = Counter()
    direct_on: list[set[str]] = []
    for x in links:
        got = crispri.of(x.chrom, x.start, x.end)
        direct_on.append({observation_id(p) for p in got})
        cls, reason, detail = ab.link_class(got, x.gene, x.action)
        direct_class[(cls, reason)] += 1
        direct_detail[detail if cls != ab.SURVIVES else "survives"] += 1
    _gate_c4(direct_class, direct_detail, c4)

    audit_a = json.loads((RESULTS_DIR / "placement_audit.json").read_text())
    per_rule: dict[str, Any] = {}
    observed: dict[str, list[Observation]] = {}
    classes: dict[str, list[tuple[str, str | None]]] = {}
    for name, rule in RULES.items():
        obs = observe(pairs, links, link_set, registry, rule)
        by_class, by_detail, on, per_link = link_census(links, obs)
        if name == CURRENT:
            _gate_c4(by_class, by_detail, c4)
            mine = [{o.id for o in on.get(x.index, ())} for x in links]
            if mine != direct_on:
                raise SystemExit("this census's index does not attach what ablation.CrispriIndex.of does")
            gate_a = _gate_audit_a(obs, audit_a["current_rule"], None)
        else:
            amb = audit_a["policy"]["ambiguity"]
            multi = sum(amb[g]["carried_by_more_than_one_element"] for g in ("positive", "null"))
            gate_a = _gate_audit_a(obs, audit_a["policy"], multi)
        if sum(by_class.values()) != len(links):
            raise SystemExit(f"{name}: the per-link census does not sum to {len(links)}")
        observed[name], classes[name] = obs, per_link
        attach = element_attachments(links, obs, on)
        if attach["supported_links"] != by_class[(ab.SURVIVES, None)]:
            raise SystemExit(f"{name}: supported links disagree with ablation.link_class")
        per_rule[name] = {
            "rule": RULE_MEANING[name],
            "observations": observation_table(obs),
            "element_level": {
                "attachments": attach,
                "c4_per_link_census": {
                    "predicted_links": sum(by_class.values()),
                    "by_class": class_rows(by_class),
                    "by_detail": dict(by_detail.most_common()),
                },
            },
            "strata": strata_tables(obs),
            "audit_a_reproduced": gate_a,
        }
        print(name, per_rule[name]["observations"]["summary"], flush=True)

    diff = difference(observed[CURRENT], observed[POLICY], classes[CURRENT], classes[POLICY])
    if any(diff["attachments_lost_under_the_policy"].values()):
        raise SystemExit(f"the policy loses attachments the current rule makes: {diff}")
    denominators = {
        "predicted_links": len(links),
        "c4_committed_predicted_links": c4["predicted_elements"],
        "screened_observations": len(pairs),
        "invalid_rows_not_read": invalid,
        "same_under_both_rules": all(
            per_rule[r]["observations"]["observations"] == len(pairs)
            and per_rule[r]["element_level"]["c4_per_link_census"]["predicted_links"] == len(links)
            for r in RULES
        ),
    }
    if not denominators["same_under_both_rules"] or len(links) != c4["predicted_elements"]:
        raise SystemExit(f"the denominators moved: {denominators}")

    r1 = resource.getrusage(resource.RUSAGE_SELF)
    peak = r1.ru_maxrss if sys.platform == "darwin" else r1.ru_maxrss * 1024
    payload = {
        "question": "C4's CRISPRi evidence attachment under the current reciprocal rule and audit A's "
        "registered policy, with each screened observation counted once as a set-valued observation",
        "status": "internal development evidence: a bounded, read-only census recorded apart, one pass. "
        "Every outcome read here was read before (prior_exposure). It replaces no headline result, "
        "changes no matcher, regenerates nothing and recommends nothing",
        "registration": registration(),
        "alphagenome_requests": 0,
        "downloads": 0,
        "per_element_response_cache_opened": False,
        "c4_reproduced": {"equal_to_committed": True, "through": ["ablation.CrispriIndex.of", "this census"]},
        "denominators": denominators,
        "rules": per_rule,
        "difference_policy_minus_current": diff,
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
