# SPDX-License-Identifier: AGPL-3.0-or-later
"""Audit A, candidate coverage: why measured CRISPRi links fail placement, cause by cause.

    uv run python scripts/placement_audit.py

Lane-place, 2026-09-29, for the coordinator (genomeos-9f). Item 13 C4 (`b7e4bf0`) sorted the 661 measured
CRISPRi decrease links with `ablation.placement`: 212 meet the overlap rule in a compiled predicted element,
62 more would meet it in an ENCODE cCRE, 387 meet it in neither. The 449 not placed are placement
failures, not 449 proven absent enhancers. This audit asks why each link fails, one cause at a time,
over every screened pair (decreases and the rest, reported apart at every step), with the ranking and
the intervals frozen.

It is a bounded audit and internal development evidence: every outcome it reads was read by earlier
lanes (PRIOR_EXPOSURE), and no resplit restores independence. No model request is made and nothing is
downloaded. It reads the compiled programs and the whole-chromosome deletion run's local element tables
(membership only: whether a registry element was scored and whether a target gene was named), never the
per-element response cache.

The one interval policy it scores (POLICY) and the rule that stops it (STOP_RULE) were registered in the
commit that added this file, before any per-outcome count of this audit was computed. Writes
data/results/placement_audit.json through the result contract.
"""

from __future__ import annotations

import bisect
import csv
import gzip
import hashlib
import json
import re
import resource
import sys
import time
from collections import Counter, defaultdict
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from genomeos import manifest as mf
from genomeos.attribution import ablation as ab
from genomeos.attribution import holdout as ho
from genomeos.attribution import measured as ms
from genomeos.results import RESULTS_DIR, save_result

NAME = "placement_audit"
REGISTERED = "2026-09-29"
ALL_ELEMENTS = Path("data/knowledge/alphagenome/all_elements")

# ==================================================================================================
# The registration (2026-09-29, lane-place), fixed before any per-outcome count of this audit.
# The only outcome figures seen before it were C4's 212 / 62 / 387 of the 661 decrease links.
# ==================================================================================================
POLICY_THRESHOLD = 0.5  # the measured layer's own value (ms.RECIPROCAL_OVERLAP); no other value is scored
POLICY = {
    "name": "min_side_half",
    "rule": "an element carries a tested interval when their overlap is at least half of the smaller of "
    "the two widths: overlap >= 0.5 * min(tested width, element width). Everything else in C4's cascade "
    "is unchanged: compiled predicted elements first, the ENCODE cCRE v3 registry second, the same "
    "intervals, the same links, no padding, no ranking, no second threshold",
    "rationale_from_conventions": [
        "measured.RECIPROCAL_OVERLAP = 0.5 asks the overlap to be half of BOTH widths, so it can pair two "
        "intervals only when their widths are within a factor of two",
        "the registry's intervals are 150 to 350 bp by its construction (ENCODE cCREs v3), and the compiled "
        "predicted elements are registry elements with the registry's coordinates",
        "the benchmark's tested intervals are set by each screen's design: most are 500 bp, the widest "
        "4,181 bp and the narrowest 52 bp (read from the files without outcomes)",
        "so under the reciprocal rule a tested interval wider than 700 bp or narrower than 75 bp cannot be "
        "carried by any registry element wherever it lies: a property of the two width conventions, not "
        "of the DNA",
        "the fraction of the smaller interval is the one of the two fractions that does not depend on "
        "which source's width convention is the wider; at equal widths it is the reciprocal rule, and it "
        "never admits a graze of a few bases",
        "the benchmark's own comparison and this project's eligibility rule (measured.ELIGIBILITY_RULE) "
        "accept any shared base; this policy is stricter than both",
    ],
    "what_it_is_not": "not a wider interval: no interval is padded, moved or merged; not tuned: one "
    "threshold, the measured layer's, chosen before any per-outcome count",
}
SHIFT_BP = 20_000  # the matched positional control: every tested interval moved by +20 kb and by -20 kb
BOOTSTRAP_DRAWS = 2000
BOOTSTRAP_SEED = 20260929
STOP_RULE = {
    "unit": "a link as C4 keys it: (chrom, start, end, gene, cell, split); positive when any record of the "
    "key is a significant decrease, null otherwise",
    "rescue": "among links the current rule leaves unplaced in the compiled programs, the fraction the "
    "policy places there; computed for positives and for nulls",
    "excess": "rescue of positives minus rescue of nulls, with a gene-clustered bootstrap 95% interval "
    f"({BOOTSTRAP_DRAWS} draws, seed {BOOTSTRAP_SEED})",
    "stop_a_nulls_as_fast": "stop when the excess interval's lower bound is <= 0: the policy gains nulls as "
    "fast as positives",
    "stop_b_not_positional": f"stop when the interval of (excess - the same excess on the same links with "
    f"every tested interval shifted by +{SHIFT_BP} and -{SHIFT_BP} bp, the mean of the two) has lower bound "
    "<= 0: the gain is explained by interval width and local element density, not by where the tested "
    "DNA lies",
    "stop_c_wide_intervals": "stop when the excess interval within tested intervals of 75 to 700 bp (where "
    "the reciprocal rule can be met at all) has lower bound <= 0: the gain rests on the intervals the "
    "reciprocal rule cannot reach, i.e. on wide intervals",
    "on_stop": "the policy is reported as denominator or overlap inflation, beside the current rule, and "
    "nothing is built on it",
    "errors": "a placement is ambiguous when more than one element of the set carries the tested interval "
    "(the rule cannot say which element the measurement is of); chance placements are the placements of "
    "the shifted intervals. Both are reported for positives and nulls under both rules",
}
BASELINE = (
    "distance to TSS: the benchmark's distanceToTSS column, score = -distance, AUROC and average "
    "precision of decreases against every other pair, on the same links under each rule"
)

# the cause cascade, applied to every link, first match wins
PLACED = "placed_compiled"
REGISTRY_WOULD = "registry_would"
ZERO_WIDTH = "coordinate_zero_width"
UNREACHABLE = "rule_width_unreachable"
WIDTH_MISMATCH = "rule_width_mismatch"
PARTIAL = "boundary_partial_overlap"
ABSENT = "absent_from_registry"
CAUSES = (PLACED, REGISTRY_WOULD, ZERO_WIDTH, UNREACHABLE, WIDTH_MISMATCH, PARTIAL, ABSENT)
CAUSE_MEANING = {
    PLACED: "a compiled predicted element meets the rule",
    REGISTRY_WOULD: "no compiled element meets the rule, a registry element does: candidate exclusion",
    ZERO_WIDTH: "the tested interval holds no base (end <= start): a coordinate defect in the source row",
    UNREACHABLE: "a registry element holds at least half of the smaller interval, but the tested interval "
    "is wider than twice the widest registry element or narrower than half the narrowest, so no registry "
    "element can meet the reciprocal rule wherever it lies: the rule and the width conventions",
    WIDTH_MISMATCH: "a registry element holds at least half of the smaller interval and the widths could "
    "meet the reciprocal rule, but these do not: the rule and the width conventions",
    PARTIAL: "a registry element shares bases with the tested interval, less than half of the smaller: a "
    "boundary offset",
    ABSENT: "no registry element shares a base with the tested interval: absent from both element sets",
}
#: why a registry element that meets the rule is not a compiled element, furthest along first. The first
#: two would be defects of this repository (a compiled element under other coordinates; a named coding
#: target the compiler dropped); the last three are the model's element selection and the sweep's scope.
EXCLUSION = (
    "compiled_under_other_coordinates",
    "coding_target_not_compiled",
    "scored_noncoding_target_only",
    "scored_no_target_named",
    "not_scored_by_the_sweep",
)
PRIOR_EXPOSURE = [
    {
        "commit": "138824f",
        "read": "the CRISPRi benchmark files against the all-element deletion (2026-09-16)",
    },
    {
        "commit": "ecc3452",
        "read": "every screened pair against the compiled elements under the reciprocal rule: the measured "
        "layer's eligibility (2026-09-17)",
    },
    {"commit": "f48b909", "read": "the Regulated pairs, for node_containment_measured (2026-09-28)"},
    {"commit": "a39073d", "read": "every pair's outcome, for crispri_published (2026-09-28)"},
    {"commit": "b7e4bf0", "read": "C4's placement of the 661 decrease links: 212 / 62 / 387 (2026-09-28)"},
    {"commit": "cc824f4", "read": "the pairs as a held-out source of the coherence programme (2026-09-28)"},
    {
        "commit": "b7cb814",
        "read": "the pairs as a held-out source of the correctness judge v1 to v3 (2026-09-29)",
    },
]


def registration() -> dict[str, Any]:
    return {
        "registered": REGISTERED,
        "policy": POLICY,
        "policy_threshold": POLICY_THRESHOLD,
        "stop_rule": STOP_RULE,
        "shift_bp": SHIFT_BP,
        "bootstrap": {"draws": BOOTSTRAP_DRAWS, "seed": BOOTSTRAP_SEED, "cluster": "measured gene"},
        "baseline": BASELINE,
        "causes": CAUSE_MEANING,
        "cause_order": list(CAUSES),
        "candidate_exclusion_order": list(EXCLUSION),
        "current_rule": f"measured.measures: reciprocal overlap >= {ms.RECIPROCAL_OVERLAP} (both fractions)",
    }


# --- the two rules --------------------------------------------------------------------------------
def overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def min_side(a_start: int, a_end: int, b_start: int, b_end: int) -> float:
    """The overlap as a fraction of the smaller interval; 0.0 when either holds no base."""
    wa, wb = a_end - a_start, b_end - b_start
    if wa <= 0 or wb <= 0:
        return 0.0
    return overlap(a_start, a_end, b_start, b_end) / min(wa, wb)


def current_rule(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return ms.measures(a_start, a_end, b_start, b_end)


def policy_rule(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return min_side(a_start, a_end, b_start, b_end) >= POLICY_THRESHOLD


Rule = Callable[[int, int, int, int], bool]


# --- an element set, searchable by a tested interval ----------------------------------------------
class ElementSet:
    """Intervals per chromosome, each with an id, found by any tested interval they share a base with."""

    def __init__(self, items: Iterable[tuple[str, int, int, str]]) -> None:
        by: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
        for chrom, s, e, i in items:
            by[chrom].append((s, e, i))
        self.by = {c: sorted(v) for c, v in by.items()}
        self.starts = {c: [s for s, _, _ in v] for c, v in self.by.items()}
        self.reach = {c: max((e - s for s, e, _ in v), default=0) for c, v in self.by.items()}
        widths = [e - s for v in self.by.values() for s, e, _ in v]
        self.min_width = min(widths, default=0)
        self.max_width = max(widths, default=0)

    def __len__(self) -> int:
        return sum(len(v) for v in self.by.values())

    def touching(self, chrom: str, s: int, e: int) -> list[tuple[int, int, str]]:
        v = self.by.get(chrom)
        if not v or e <= s:
            return []
        st = self.starts[chrom]
        lo = bisect.bisect_left(st, s - self.reach[chrom])
        hi = bisect.bisect_left(st, e)
        return [x for x in v[lo:hi] if overlap(x[0], x[1], s, e) > 0]


@dataclass(frozen=True)
class Hit:
    """How one element set meets one tested interval."""

    touching: int
    min_side_half: int  # elements holding at least half of the smaller interval
    by_rule: tuple[str, ...]  # ids of the elements that meet the rule asked for

    @property
    def status(self) -> str:
        if self.by_rule:
            return "rule"
        if self.min_side_half:
            return "min_side"
        if self.touching:
            return "partial"
        return "none"


def hit(es: ElementSet, chrom: str, s: int, e: int, rule: Rule) -> Hit:
    t = es.touching(chrom, s, e)
    return Hit(
        touching=len(t),
        min_side_half=sum(1 for a, b, _ in t if min_side(a, b, s, e) >= POLICY_THRESHOLD),
        by_rule=tuple(i for a, b, i in t if rule(a, b, s, e)),
    )


def cause(comp: Hit, reg: Hit, width: int, reg_min: int, reg_max: int) -> str:
    """The first cause of the cascade that holds for one link (CAUSES, CAUSE_MEANING)."""
    if comp.status == "rule":
        return PLACED
    if reg.status == "rule":
        return REGISTRY_WOULD
    if width <= 0:
        return ZERO_WIDTH
    if reg.status == "min_side":
        return UNREACHABLE if width > 2 * reg_max or 2 * width < reg_min else WIDTH_MISMATCH
    if reg.status == "partial":
        return PARTIAL
    return ABSENT


def exclusion(ids: Iterable[str], compiled_ids: set[str], swept: dict[str, str]) -> str:
    """Why the registry elements that meet the rule are not compiled elements: the state furthest along."""
    states = [EXCLUSION[0] if i in compiled_ids else swept.get(i, EXCLUSION[4]) for i in ids]
    return min(states, key=EXCLUSION.index) if states else EXCLUSION[4]


# --- the screened pairs ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Link:
    chrom: str
    start: int
    end: int
    gene: str
    cell: str
    split: str
    positive: bool
    outcomes: tuple[str, ...]
    dataset: str
    distance: float  # the benchmark's distanceToTSS; NaN when absent

    @property
    def key(self) -> tuple[str, int, int, str, str, str]:
        return (self.chrom, self.start, self.end, self.gene, self.cell, self.split)

    @property
    def group(self) -> str:
        """positive, or the null's outcome when all its records agree, else mixed_null."""
        if self.positive:
            return "positive"
        return self.outcomes[0] if len(set(self.outcomes)) == 1 else "mixed_null"


def _float(v: str | None) -> float:
    try:
        return float(v) if v not in (None, "", "NA") else float("nan")
    except ValueError:
        return float("nan")


def read_pairs(knowledge: Path = ms.CRISPRI_KNOWLEDGE) -> tuple[list[ms.CrispriPair], list[dict[str, str]]]:
    """Every valid pair through measured.parse_crispri, with its raw row beside it (same order, checked)."""
    pairs: list[ms.CrispriPair] = []
    raw: list[dict[str, str]] = []
    for name in ms.CRISPRI_FILES:
        p = knowledge / name
        if not p.exists():
            continue
        with gzip.open(p, "rt") as fh:
            got, _ = ms.parse_crispri(fh, None, ms.CRISPRI_SPLIT_OF[name], name)
        with gzip.open(p, "rt") as fh:
            rows = [r for r in csv.DictReader(fh, delimiter="\t") if r.get("ValidConnection") == "TRUE"]
        if len(rows) != len(got) or any(
            (r["chrom"], int(r["chromStart"]), int(r["chromEnd"]), r["measuredGeneSymbol"], r["CellType"])
            != (g.chrom, g.start, g.end, g.gene, g.cell)
            for r, g in zip(rows, got, strict=True)
        ):
            raise ValueError(f"{name}: raw rows do not align with measured.parse_crispri")
        pairs += got
        raw += rows
    return pairs, raw


def links_of(pairs: list[ms.CrispriPair], raw: list[dict[str, str]]) -> list[Link]:
    """One link per C4 key; positive when any of its records is a significant decrease (as C4 reads it)."""
    recs: dict[tuple, list[tuple[ms.CrispriPair, dict[str, str]]]] = defaultdict(list)
    for p, r in zip(pairs, raw, strict=True):
        recs[(p.chrom, p.start, p.end, p.gene, p.cell, p.split)].append((p, r))
    out = []
    for k, rs in recs.items():
        outs = tuple(sorted({p.outcome for p, _ in rs}))
        dist = [_float(r.get("distanceToTSS")) for _, r in rs]
        dist = [d for d in dist if d == d]
        out.append(
            Link(
                *k,
                positive=ms.DECREASE in outs,
                outcomes=outs,
                dataset="+".join(sorted({p.dataset for p, _ in rs})),
                distance=float(np.median(dist)) if dist else float("nan"),
            )
        )
    out.sort(key=lambda x: x.key)
    return out


def fold_hash(links: list[Link], split: str) -> str:
    h = hashlib.sha256()
    for x in links:
        if x.split == split:
            h.update(("\t".join(map(str, (*x.key, int(x.positive)))) + "\n").encode())
    return h.hexdigest()


# --- the element sets -----------------------------------------------------------------------------
def compiled_elements(programs: list[Path]) -> list[tuple[str, int, int, str]]:
    """Every predicted element block's locus (not `_measured`), as ablation.program_census reads it."""
    out = []
    for path in programs:
        for kind, name, f in ho.program_blocks(path):
            if kind != "element" or name.endswith("_measured"):
                continue
            m = ho.LOCUS_RE.match(f.get("locus", ""))
            if m:
                out.append((m.group(1), int(m.group(2)), int(m.group(3)), name))
    return out


def registry_elements(results_dir: Path = RESULTS_DIR) -> tuple[list[tuple[str, int, int, str]], Counter]:
    out = []
    classes: Counter = Counter()
    for chrom in ho.CHROMS:
        p = results_dir / f"ccres_{chrom}.bed.gz"
        if not p.exists():
            continue
        with gzip.open(p, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                out.append((f[0], int(f[1]), int(f[2]), f[3]))
                classes[f[4]] += 1
    return out, classes


def swept_states(ids_by_chrom: dict[str, set[str]], where: Path = ALL_ELEMENTS) -> dict[str, str]:
    """For the registry ids asked for: whether the whole-chromosome deletion run scored them and named a
    target. Membership only; no effect value is kept."""
    out: dict[str, str] = {}
    for chrom, ids in ids_by_chrom.items():
        p = where / f"{chrom}.json"
        if not ids or not p.exists():
            continue
        t = json.loads(p.read_text())
        t = t if isinstance(t, list) else t.get("elements", [])
        for e in t:
            if e.get("id") not in ids:
                continue
            if (e.get("predicted_coding") or {}).get("gene"):
                out[e["id"]] = EXCLUSION[1]
            elif (e.get("predicted") or {}).get("gene"):
                out[e["id"]] = EXCLUSION[2]
            else:
                out[e["id"]] = EXCLUSION[3]
        del t
    return out


# --- scoring --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Placed:
    link: Link
    cause: str
    comp: Hit
    reg: Hit


def place(links: list[Link], comp: ElementSet, reg: ElementSet, rule: Rule, shift: int = 0) -> list[Placed]:
    out = []
    for x in links:
        s, e = x.start + shift, x.end + shift
        c, r = hit(comp, x.chrom, s, e, rule), hit(reg, x.chrom, s, e, rule)
        out.append(Placed(x, cause(c, r, e - s, reg.min_width, reg.max_width), c, r))
    return out


GROUPS = (
    "positive",
    "null",
    ms.NULL_INFORMATIVE,
    ms.NULL_INCONCLUSIVE,
    ms.INCREASE,
    ms.MISSING,
    "mixed_null",
)


def attrition(
    rows: list[Placed], extra: dict[int, str] | None = None
) -> dict[str, dict[str, dict[str, int]]]:
    """Counts per group, split and cause; `extra` adds a sub-cause per row index (candidate exclusion)."""
    out: dict[str, dict[str, Counter]] = {g: defaultdict(Counter) for g in GROUPS}
    for i, p in enumerate(rows):
        groups = ["positive"] if p.link.positive else ["null", p.link.group]
        label = p.cause if not (extra and i in extra) else f"{p.cause}/{extra[i]}"
        for g in groups:
            for sp in (p.link.split, ms.ALL):
                out[g][sp][label] += 1
                out[g][sp]["total"] += 1
    return {g: {sp: dict(sorted(c.items())) for sp, c in v.items()} for g, v in out.items() if v}


def auroc(scores: np.ndarray, labels: np.ndarray) -> float | None:
    """The Mann-Whitney AUROC, ties counted half."""
    pos, neg = scores[labels], scores[~labels]
    if len(pos) == 0 or len(neg) == 0:
        return None
    allv = np.concatenate([pos, neg])
    order = np.argsort(allv, kind="mergesort")
    sorted_v = allv[order]
    ranks = np.empty(len(allv))
    i = 0
    while i < len(sorted_v):
        j = i
        while j + 1 < len(sorted_v) and sorted_v[j + 1] == sorted_v[i]:
            j += 1
        ranks[order[i : j + 1]] = (i + j) / 2 + 1
        i = j + 1
    return float((ranks[: len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def average_precision(scores: np.ndarray, labels: np.ndarray) -> float | None:
    if labels.sum() == 0:
        return None
    order = np.argsort(-scores, kind="mergesort")
    y = labels[order]
    prec = np.cumsum(y) / np.arange(1, len(y) + 1)
    return float((prec * y).sum() / y.sum())


def baseline(links: list[Link]) -> dict[str, Any]:
    have = [x for x in links if x.distance == x.distance]
    if not have:
        return {"links": 0, "without_distance": len(links)}
    d = np.array([abs(x.distance) for x in have], dtype=float)
    y = np.array([x.positive for x in have], dtype=bool)
    a, ap = auroc(-d, y), average_precision(-d, y)
    return {
        "links": len(have),
        "without_distance": len(links) - len(have),
        "positives": int(y.sum()),
        "prevalence": round(float(y.mean()), 4),
        "auroc": None if a is None else round(a, 4),
        "average_precision": None if ap is None else round(ap, 4),
        "median_distance_positive": float(np.median(d[y])) if y.any() else None,
        "median_distance_null": float(np.median(d[~y])) if (~y).any() else None,
    }


def cluster_bootstrap(
    genes: list[str], columns: dict[str, np.ndarray], stat: Callable[[dict[str, float]], float]
) -> tuple[float, float]:
    """A 95% interval for `stat` of column sums over links, resampling measured genes with replacement."""
    uniq = sorted(set(genes))
    idx_of = {g: i for i, g in enumerate(uniq)}
    gi = np.array([idx_of[g] for g in genes])
    sums = {k: np.bincount(gi, weights=v.astype(float), minlength=len(uniq)) for k, v in columns.items()}
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    vals = []
    for _ in range(BOOTSTRAP_DRAWS):
        w = np.bincount(rng.integers(0, len(uniq), len(uniq)), minlength=len(uniq))
        v = stat({k: float((s * w).sum()) for k, s in sums.items()})
        if v == v:
            vals.append(v)
    if not vals:
        return (float("nan"), float("nan"))
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return (round(float(lo), 4), round(float(hi), 4))


def _rate(a: float, b: float) -> float:
    return a / b if b else float("nan")


def rescue_columns(
    cur: list[Placed], pol: list[Placed], cur_s: list[Placed] | None = None, pol_s: list[Placed] | None = None
) -> dict[str, np.ndarray]:
    """Per-link 0/1 columns for the rescue statistics; with shifted rows, the same on the shifted links."""
    pos = np.array([p.link.positive for p in cur], dtype=bool)
    un = np.array([p.cause != PLACED for p in cur], dtype=bool)
    res = un & np.array([p.cause == PLACED for p in pol], dtype=bool)
    cols = {"pos_un": pos & un, "null_un": ~pos & un, "pos_res": pos & res, "null_res": ~pos & res}
    if cur_s is not None and pol_s is not None:
        un_s = np.array([p.cause != PLACED for p in cur_s], dtype=bool)
        res_s = un_s & np.array([p.cause == PLACED for p in pol_s], dtype=bool)
        cols |= {
            "pos_un_s": pos & un_s,
            "null_un_s": ~pos & un_s,
            "pos_res_s": pos & res_s,
            "null_res_s": ~pos & res_s,
        }
    return cols


def excess(s: dict[str, float], suffix: str = "") -> float:
    """Rescue of positives minus rescue of nulls, from column sums (`suffix` picks a shifted copy)."""
    return _rate(s["pos_res" + suffix], s["pos_un" + suffix]) - _rate(
        s["null_res" + suffix], s["null_un" + suffix]
    )


def ambiguity(rows: list[Placed]) -> dict[str, dict[str, int]]:
    """Placed links carried by more than one compiled element, positives and nulls."""
    out: dict[str, Counter] = {"positive": Counter(), "null": Counter()}
    for p in rows:
        if p.cause != PLACED:
            continue
        g = "positive" if p.link.positive else "null"
        out[g]["placed"] += 1
        out[g]["carried_by_more_than_one_element"] += len(p.comp.by_rule) > 1
    return {g: dict(c) for g, c in out.items()}


def coverage(rows: list[Placed]) -> dict[str, Any]:
    out = {}
    for g, want in (("positive", True), ("null", False)):
        sel = [p for p in rows if p.link.positive is want]
        n = len(sel)
        placed = sum(p.cause == PLACED for p in sel)
        reg = sum(p.cause == REGISTRY_WOULD for p in sel)
        out[g] = {
            "links": n,
            "placed_compiled": placed,
            "registry_would": reg,
            "coverage_compiled": round(_rate(placed, n), 4),
            "coverage_compiled_or_registry": round(_rate(placed + reg, n), 4),
            "abstention": round(1 - _rate(placed, n), 4),
        }
    return out


# --- the conventions of each input ----------------------------------------------------------------
NAME_RE = re.compile(r"\|(chr[0-9XY]+):(\d+)-(\d+)")


def crispri_conventions(raw: list[dict[str, str]], reg: ElementSet) -> dict[str, Any]:
    """Build and coordinate checks on the benchmark rows, outcome-blind (no outcome column is read)."""
    by_ds: dict[str, Counter] = defaultdict(Counter)
    widths: dict[str, list[int]] = defaultdict(list)
    uniq_now: set[tuple[str, int, int]] = set()
    uniq_name: set[tuple[str, int, int]] = set()
    zero = []
    for r in raw:
        s, e, ds = int(r["chromStart"]), int(r["chromEnd"]), r["Dataset"]
        c = by_ds[ds]
        c["rows"] += 1
        widths[ds].append(e - s)
        c["odd_width"] += (e - s) % 2
        uniq_now.add((r["chrom"], s, e))
        if e <= s:
            zero.append({"chrom": r["chrom"], "start": s, "end": e, "name": r["name"], "dataset": ds})
        m = NAME_RE.search(r["name"])
        if not m:
            c["name_without_coordinates"] += 1
        else:
            ns, ne = int(m.group(2)), int(m.group(3))
            c["name_coords_equal" if (m.group(1), ns) == (r["chrom"], s) else "name_coords_differ"] += 1
            c["width_kept" if ne - ns == e - s else "width_changed"] += 1
            if ns != s:
                uniq_name.add((m.group(1), ns, ne))
        c["tss_one_base" if int(r["endTSS"]) - int(r["startTSS"]) == 1 else "tss_not_one_base"] += 1

    def touching_rate(ivs: set[tuple[str, int, int]], shift: int = 0) -> float:
        n = sum(1 for c, s, e in ivs if reg.touching(c, s + shift, e + shift))
        return round(_rate(n, len(ivs)), 4)

    per = {}
    for ds, c in sorted(by_ds.items()):
        w = sorted(widths[ds])
        per[ds] = {
            **dict(sorted(c.items())),
            "width_min": w[0],
            "width_median": w[len(w) // 2],
            "width_max": w[-1],
            "wider_than_700": sum(x > 700 for x in w),
            "narrower_than_75": sum(x < 75 for x in w),
        }
    return {
        "per_dataset_rows": per,
        "zero_width_rows": zero,
        "distinct_tested_intervals": len(uniq_now),
        "registry_touching_rate": {
            "as_read_grch38_columns": touching_rate(uniq_now),
            f"shifted_plus_{SHIFT_BP}": touching_rate(uniq_now, SHIFT_BP),
            f"shifted_minus_{SHIFT_BP}": touching_rate(uniq_now, -SHIFT_BP),
            "name_field_coordinates_where_they_differ": touching_rate(uniq_name),
            "name_field_intervals": len(uniq_name),
        },
    }


def tss_check(raw: list[dict[str, str]]) -> dict[str, Any]:
    """The benchmark's TSS against GENCODE v50's gene TSS (0-based) by Ensembl id, one row per gene."""
    tss = {g.gene_id: g.tss for g in ho.genes()}
    seen: dict[str, int] = {}
    for r in raw:
        gid = r["measuredGeneEnsemblId"].split(".")[0]
        if gid and int(r["endTSS"]) - int(r["startTSS"]) == 1:
            seen.setdefault(gid, int(r["startTSS"]))
    diffs = np.array([abs(v - tss[g]) for g, v in seen.items() if g in tss])
    return {
        "genes": len(seen),
        "in_gencode_v50": int(len(diffs)),
        "exact": int((diffs == 0).sum()),
        "within_1bp": int((diffs <= 1).sum()),
        "within_1kb": int((diffs <= 1000).sum()),
        "median_abs_difference_bp": float(np.median(diffs)) if len(diffs) else None,
    }


def compiled_conventions(comp: list[tuple[str, int, int, str]], reg: list[tuple[str, int, int, str]]) -> dict:
    """Whether every compiled predicted element is a registry element with the registry's coordinates."""
    where = {i: (c, s, e) for c, s, e, i in reg}
    same = differ = absent = 0
    widths = []
    for c, s, e, i in comp:
        widths.append(e - s)
        if i not in where:
            absent += 1
        elif where[i] == (c, s, e):
            same += 1
        else:
            differ += 1
    return {
        "compiled_predicted_elements": len(comp),
        "same_id_same_coordinates_in_registry": same,
        "same_id_other_coordinates": differ,
        "id_not_in_registry": absent,
        "width_min": min(widths, default=None),
        "width_max": max(widths, default=None),
    }


def off_by_one(links: list[Link], comp: ElementSet, reg: ElementSet, base: list[Placed]) -> dict[str, Any]:
    """How many links change cause under the current rule when one end moves by one base."""
    out = {}
    for label, ds, de in (("start_minus_1", -1, 0), ("start_plus_1", 1, 0), ("end_plus_1", 0, 1)):
        moved: Counter = Counter()
        for x, b in zip(links, base, strict=True):
            s, e = x.start + ds, x.end + de
            c, r = hit(comp, x.chrom, s, e, current_rule), hit(reg, x.chrom, s, e, current_rule)
            if cause(c, r, e - s, reg.min_width, reg.max_width) != b.cause:
                moved["positive" if x.positive else "null"] += 1
        out[label] = {"positive": moved["positive"], "null": moved["null"]}
    return out


# --- the manifest ---------------------------------------------------------------------------------
@mf.depends_on_models("alphagenome")  # it reads the compiled predicted layer and the sweep's membership
def manifest(programs: list[Path], swept: list[Path]) -> dict[str, Any]:
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
    inputs += [mf.input_entry(p, partition=None) for p in ho.gencode_paths()]
    if swept:
        inputs.append(mf.files_entry("data/knowledge/alphagenome/all_elements (membership only)", swept))
    inputs.append(mf.input_entry(RESULTS_DIR / "c4_alphagenome_ablation.json", partition=None))
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
            {
                "accession": "this repository: the whole-chromosome deletion run's element tables "
                "(generated, untracked), read for membership only",
                "version": "pinned by sha256",
            },
            {"accession": "GENCODE v50", "version": "pinned by sha256"},
            {"accession": "this repository: data/results/c4_alphagenome_ablation.json", "version": "b7e4bf0"},
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "current_rule_reciprocal_overlap": ms.RECIPROCAL_OVERLAP,
            "policy": POLICY["name"],
            "policy_threshold": POLICY_THRESHOLD,
            "shift_bp": SHIFT_BP,
            "bootstrap_draws": BOOTSTRAP_DRAWS,
            "bootstrap_seed": BOOTSTRAP_SEED,
            "registered": REGISTERED,
        },
        "exclusions": [
            "CRISPRi pairs the benchmark marks as not a valid connection are not read (none in these files)",
            "the per-element response cache is not opened; the model is not called; nothing is downloaded",
            "no effect value of the deletion run is read, only whether an element was scored and named "
            "a gene",
        ],
        "partitions": {
            ms.TRAINING: "the CRISPRi benchmark's training file (K562); read as measurement, nothing fitted",
            ms.HELDOUT: "the CRISPRi benchmark's held-out file; read as measurement, reported apart",
        },
    }


# --- the run --------------------------------------------------------------------------------------
def _stop(lo: float) -> bool:
    return lo != lo or lo <= 0


def policy_block(
    links: list[Link],
    comp: ElementSet,
    reg: ElementSet,
    cur: list[Placed],
    swept: dict[str, str],
    compiled_ids: set[str],
) -> dict[str, Any]:
    """The registered policy beside the current rule, with the three registered stops."""
    pol = place(links, comp, reg, policy_rule)
    genes = [x.gene for x in links]
    shifted = {
        sh: (place(links, comp, reg, current_rule, sh), place(links, comp, reg, policy_rule, sh))
        for sh in (SHIFT_BP, -SHIFT_BP)
    }
    cols = rescue_columns(cur, pol)
    both = dict(cols)
    for sh, tag in ((SHIFT_BP, "p"), (-SHIFT_BP, "m")):
        c = rescue_columns(cur, pol, *shifted[sh])
        for k in ("pos_un", "null_un", "pos_res", "null_res"):
            both[f"{k}_s{tag}"] = c[f"{k}_s"]

    def diff_shift(s: dict[str, float]) -> float:
        return excess(s) - (excess(s, "_sp") + excess(s, "_sm")) / 2

    point = {k: float(v.sum()) for k, v in both.items()}
    ex, ex_ci = excess(point), cluster_bootstrap(genes, cols, excess)
    ds_point, ds_ci = diff_shift(point), cluster_bootstrap(genes, both, diff_shift)
    reachable = np.array(
        [2 * (x.end - x.start) >= reg.min_width and (x.end - x.start) <= 2 * reg.max_width for x in links]
    )
    narrow = {k: v & reachable for k, v in cols.items()}
    narrow_point = {k: float(v.sum()) for k, v in narrow.items()}
    narrow_ci = cluster_bootstrap(genes, narrow, excess)
    stops = {
        "a_nulls_as_fast": _stop(ex_ci[0]),
        "b_not_positional": _stop(ds_ci[0]),
        "c_wide_intervals": _stop(narrow_ci[0]),
    }
    pol_ex = {
        i: exclusion(p.reg.by_rule, compiled_ids, swept)
        for i, p in enumerate(pol)
        if p.cause == REGISTRY_WOULD
    }

    def counts(c: dict[str, np.ndarray]) -> dict[str, Any]:
        s = {k: int(v.sum()) for k, v in c.items()}
        return {
            "positives_unplaced_by_current_rule": s["pos_un"],
            "positives_rescued": s["pos_res"],
            "rescue_positive": round(_rate(s["pos_res"], s["pos_un"]), 4),
            "nulls_unplaced_by_current_rule": s["null_un"],
            "nulls_rescued": s["null_res"],
            "rescue_null": round(_rate(s["null_res"], s["null_un"]), 4),
        }

    return {
        "scored": True,
        "attrition": attrition(pol, pol_ex),
        "coverage": coverage(pol),
        "ambiguity": ambiguity(pol),
        "rescue": {
            **counts(cols),
            "excess": round(ex, 4),
            "excess_ci95": list(ex_ci),
            "shifted_excess": {
                f"+{SHIFT_BP}": round(excess(point, "_sp"), 4),
                f"-{SHIFT_BP}": round(excess(point, "_sm"), 4),
            },
            "excess_minus_shifted": round(ds_point, 4),
            "excess_minus_shifted_ci95": list(ds_ci),
            "within_75_to_700_bp": {
                **counts(narrow),
                "excess": round(excess(narrow_point), 4),
                "excess_ci95": list(narrow_ci),
            },
            "by_tested_width_descriptive": by_width(cur, pol),
        },
        "shifted_controls": {
            f"{sh:+d}": {
                "current_rule": coverage(cs),
                "policy": coverage(ps),
                "policy_ambiguity": ambiguity(ps),
            }
            for sh, (cs, ps) in shifted.items()
        },
        "baseline_distance_to_tss": {
            "placed_policy": baseline([p.link for p in pol if p.cause == PLACED]),
            "rescued_by_policy": baseline(
                [p.link for p, q in zip(cur, pol, strict=True) if p.cause != PLACED and q.cause == PLACED]
            ),
        },
        "stop": stops,
        "stopped": any(stops.values()),
    }


FAMILIES = {
    "placed": (PLACED,),
    "candidate_exclusion": (REGISTRY_WOULD,),
    "coordinate_defect": (ZERO_WIDTH,),
    "rule_and_width_conventions": (UNREACHABLE, WIDTH_MISMATCH),
    "boundary_offset": (PARTIAL,),
    "absent_from_both_sets": (ABSENT,),
}


def cause_families(rows: list[Placed]) -> dict[str, dict[str, int]]:
    """The causes grouped into the audit's four questions, positives and nulls apart."""
    fam_of = {c: f for f, cs in FAMILIES.items() for c in cs}
    out: dict[str, Counter] = {"positive": Counter(), "null": Counter()}
    for p in rows:
        out["positive" if p.link.positive else "null"][fam_of[p.cause]] += 1
    return {g: {f: c[f] for f in FAMILIES} for g, c in out.items()}


def width_class(width: int, reg_min: int = 150, reg_max: int = 350) -> str:
    if 2 * width < reg_min:
        return "narrower_than_half_the_narrowest"
    if width > 2 * reg_max:
        return "wider_than_twice_the_widest"
    return "reachable_width"


def by_width(cur: list[Placed], pol: list[Placed]) -> dict[str, dict[str, dict[str, int]]]:
    """Descriptive, added after the first run and not registered: links, placed under the current rule and
    rescued by the second rule, per tested-width class, positives and nulls apart."""
    out: dict[str, dict[str, Counter]] = {"positive": defaultdict(Counter), "null": defaultdict(Counter)}
    for p, q in zip(cur, pol, strict=True):
        c = out["positive" if p.link.positive else "null"][width_class(p.link.end - p.link.start)]
        c["links"] += 1
        c["placed_current_rule"] += p.cause == PLACED
        c["unplaced_current_rule"] += p.cause != PLACED
        c["placed_second_rule"] += q.cause == PLACED
        c["rescued"] += p.cause != PLACED and q.cause == PLACED
        c["ambiguous_second_rule"] += q.cause == PLACED and len(q.comp.by_rule) > 1
    return {g: {w: dict(c) for w, c in sorted(v.items())} for g, v in out.items()}


def group_counts(links: list[Link]) -> dict[str, dict[str, int]]:
    out = {}
    for g in ("positive", "null"):
        sel = [x for x in links if x.positive is (g == "positive")]
        out[g] = {
            "links": len(sel),
            "distinct_genes": len({x.gene for x in sel}),
            "distinct_tested_intervals": len({(x.chrom, x.start, x.end) for x in sel}),
        }
    return out


def main() -> None:
    t0, r0 = time.perf_counter(), resource.getrusage(resource.RUSAGE_SELF)
    programs = ho.compiled_programs()
    pairs, raw = read_pairs()
    links = links_of(pairs, raw)
    comp_items = compiled_elements(programs)
    reg_items, reg_classes = registry_elements()
    comp, reg = ElementSet(comp_items), ElementSet(reg_items)
    print("links", len(links), "compiled", len(comp), "registry", len(reg), flush=True)

    # C4's own numbers, reproduced before anything else is read off this run
    c4 = ab.placement(pairs, {c: [(s, e) for s, e, _ in v] for c, v in comp.by.items()}, ab.ccre_intervals())
    committed = json.loads((RESULTS_DIR / "c4_alphagenome_ablation.json").read_text())["placement"]
    if {k: v for k, v in committed.items() if k != "what"} != c4:
        raise SystemExit(f"C4's placement does not reproduce: {c4} against {committed}")
    cur = place(links, comp, reg, current_rule)
    pos_cur = [p for p in cur if p.link.positive]
    on_decreases = {
        "placed": sum(p.cause == PLACED for p in pos_cur),
        "registry_would": sum(p.cause == REGISTRY_WOULD for p in pos_cur),
        "neither": sum(p.cause not in (PLACED, REGISTRY_WOULD) for p in pos_cur),
    }
    want = {
        "placed": sum(v for k, v in c4.items() if k.startswith("placed_")),
        "registry_would": sum(v for k, v in c4.items() if k.startswith("not_placed_registry_would_")),
        "neither": sum(v for k, v in c4.items() if k.startswith("not_placed_no_registry_element_")),
    }
    if on_decreases != want:
        raise SystemExit(f"the cascade does not reproduce C4 on the decrease links: {on_decreases} vs {want}")

    need: dict[str, set[str]] = defaultdict(set)
    for p in cur:
        for _, _, i in reg.touching(p.link.chrom, p.link.start, p.link.end):
            need[p.link.chrom].add(i)
    swept_paths = [ALL_ELEMENTS / f"{c}.json" for c in ho.CHROMS if (ALL_ELEMENTS / f"{c}.json").exists()]
    swept = swept_states(need)
    compiled_ids = {i for _, _, _, i in comp_items}
    cur_ex = {
        i: exclusion(p.reg.by_rule, compiled_ids, swept)
        for i, p in enumerate(cur)
        if p.cause == REGISTRY_WOULD
    }
    touched_state = {
        g: Counter(
            "compiled" if i in compiled_ids else swept.get(i, EXCLUSION[4])
            for i in {
                t[2]
                for p in cur
                if p.link.positive is (g == "positive")
                for t in reg.touching(p.link.chrom, p.link.start, p.link.end)
            }
        )
        for g in ("positive", "null")
    }
    fails = Counter(
        ("positive" if p.link.positive else "null", p.cause, p.comp.status)
        for p in cur
        if p.cause not in (PLACED, REGISTRY_WOULD)
    )
    reg_paths = [
        RESULTS_DIR / f"ccres_{c}.bed.gz" for c in ho.CHROMS if (RESULTS_DIR / f"ccres_{c}.bed.gz").exists()
    ]
    conventions = {
        "crispri": {
            "sha256": {n: mf.input_entry(ms.CRISPRI_KNOWLEDGE / n)["sha256"] for n in ms.CRISPRI_FILES},
            "declared_build": "GRCh38 (the file names; chrom/chromStart/chromEnd are the perturbed element, "
            "the name column keeps each source study's own coordinates)",
            "declared_convention": "BED: 0-based, half-open (chromStart/chromEnd; the TSS as a one-base "
            "startTSS/endTSS interval)",
            "found": crispri_conventions(raw, reg),
            "tss_against_gencode_v50": tss_check(raw),
        },
        "registry": {
            "sha256_of_24_files": mf.files_entry("data/results/ccres_<chrom>.bed.gz", reg_paths)["sha256"],
            "declared_build": "GRCh38 (each file's header: ENCODE cCREs v3 from "
            "downloads.wenglab.org/V3/GRCh38-cCREs.bed)",
            "declared_convention": "BED: 0-based, half-open, written unchanged by "
            "genome/regulatory.save_ccres",
            "elements": len(reg),
            "by_class": dict(sorted(reg_classes.items())),
            "width_min": reg.min_width,
            "width_max": reg.max_width,
            "chromosomes": len(reg.by),
        },
        "compiled": {
            "sha256_of_24_programs": mf.files_entry(
                "data/knowledge/compiled/noncoding_<chrom>.bio", programs
            )["sha256"],
            "declared_build": "GRCh38 (the registry's coordinates, carried by id)",
            "declared_convention": "locus chrom:start-end, the registry's 0-based half-open numbers "
            "unchanged",
            "found": compiled_conventions(comp_items, reg_items),
        },
    }
    obo = off_by_one(links, comp, reg, cur)
    genes = [x.gene for x in links]
    pos = np.array([x.positive for x in links], dtype=bool)
    placed = np.array([p.cause == PLACED for p in cur], dtype=bool)
    cov_cols = {"pos": pos, "null": ~pos, "pos_pl": pos & placed, "null_pl": ~pos & placed}

    def cov_diff(s: dict[str, float]) -> float:
        return _rate(s["pos_pl"], s["pos"]) - _rate(s["null_pl"], s["null"])

    policy = policy_block(links, comp, reg, cur, swept, compiled_ids)

    r1 = resource.getrusage(resource.RUSAGE_SELF)
    peak = r1.ru_maxrss if sys.platform == "darwin" else r1.ru_maxrss * 1024
    payload = {
        "question": "why measured CRISPRi links fail placement in C4, cause by cause, over every "
        "screened pair",
        "status": "internal development evidence: every outcome read here was read by earlier lanes "
        "(prior_exposure); no resplit restores independence. A bounded audit, not a biological result",
        "registration": registration(),
        "alphagenome_requests": 0,
        "downloads": 0,
        "per_element_response_cache_opened": False,
        "c4_reproduced": {"placement": c4, "equal_to_committed": True, "cascade_on_decreases": on_decreases},
        "folds": {
            ms.TRAINING: fold_hash(links, ms.TRAINING),
            ms.HELDOUT: fold_hash(links, ms.HELDOUT),
            "what": "sha256 over each split's sorted link keys and positive flag",
        },
        "conventions": conventions,
        "off_by_one_sensitivity_current_rule": obo,
        "current_rule": {
            "attrition": attrition(cur, cur_ex),
            "coverage": coverage(cur),
            "coverage_difference_positive_minus_null": round(
                cov_diff({k: float(v.sum()) for k, v in cov_cols.items()}), 4
            ),
            "coverage_difference_ci95": list(cluster_bootstrap(genes, cov_cols, cov_diff)),
            "ambiguity": ambiguity(cur),
            "failures_by_compiled_status": [
                {"group": g, "cause": c, "compiled": s, "n": n} for (g, c, s), n in sorted(fails.items())
            ],
            "registry_elements_touched_by_a_tested_interval": {
                g: dict(sorted(c.items())) for g, c in touched_state.items()
            },
            "cause_families": cause_families(cur),
            "by_tested_width_descriptive": by_width(cur, cur),
            "baseline_distance_to_tss": {
                "all_links": baseline(links),
                "placed_current_rule": baseline([p.link for p in cur if p.cause == PLACED]),
                "not_placed_current_rule": baseline([p.link for p in cur if p.cause != PLACED]),
            },
        },
        "policy": policy,
        "independent_loci": group_counts(links),
        "positive_links_with_another_outcome_too": sum(
            1 for x in links if x.positive and len(x.outcomes) > 1
        ),
        "prior_exposure": PRIOR_EXPOSURE,
        "cost": {
            "wall_seconds": round(time.perf_counter() - t0, 1),
            "cpu_seconds": round((r1.ru_utime - r0.ru_utime) + (r1.ru_stime - r0.ru_stime), 1),
            "peak_memory_mb": round(peak / 2**20, 1),
            "downloads": 0,
            "requests": 0,
        },
    }
    p = save_result(NAME, payload, manifest=manifest(programs, swept_paths))
    print("saved", p, flush=True)


if __name__ == "__main__":
    main()
