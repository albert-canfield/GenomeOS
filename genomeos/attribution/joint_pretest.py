# SPDX-License-Identifier: AGPL-3.0-or-later
"""R8's pretest: does a coupling term beat independent per-element scoring at naming CRISPRi targets?

Review item R8 (docs/ROADMAP.md section 5 item 11) allows a joint inference engine only if a coupling
term carries information that independent per-element scoring does not. This module is that test and
nothing more; no search is built here. Registered 2026-09-28 before any number was produced
(docs/ATTRIBUTION.md, "Pre-registration: R8's coupling pretest").

What is compared, on the CRISPRi benchmark's training file only (`measured.load_crispri(split=
"training")`, passed through `development_only`; the held-out file is never opened):

- **independent**: a pair's predicted drop, the largest predicted fall of the measured gene on the
  cell's own track over the cached elements overlapping the tested interval, floored at zero
  (`crispri.deletion_drop`), read from the per-element response cache (`targets.ElementResponses`).
  An interval no cached element overlaps, or where no overlapping element scored the gene in that
  cell, is missing and is never imputed.
- **element competition**: the pair's share of its perturbation's predicted drops, `d / (sum d + TAU)`
  over every gene the screen tested against the same interval in the same cell (the within-
  perturbation contrast of genomeos-8a's brief).
- **gene budget**: the pair's share of its gene's predicted drops, `d / (sum d + TAU)` over every
  interval the screen tested against that gene in that cell.

Each coupling variant enters a logistic model `[d, share]` fitted leave-one-chromosome-out; the
independent score is `d` itself (a one-feature logistic ranks exactly as `d` does). A self-only share
`d / (d + TAU)` is the control: it has the same nonlinearity and no partner, so a gain it also shows
is not coupling. Outcomes follow R2: a significant decrease is a positive, a well-powered null a
negative, and increases and underpowered nulls are counted apart and never scored.
"""

from __future__ import annotations

import bisect
import json
import random
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from genomeos.attribution import crispri, measured
from genomeos.attribution.targets import ELEMENT_CACHE, ElementResponses

# --- the registration (fixed 2026-09-28, before the run) ------------------------------------------
REGISTERED = "2026-09-28"
TAU = 0.1  # the share's floor, in log2 units: the scorer's own MIN_EFFECT bar, not fitted
VARIANTS = ("element_competition", "gene_budget")
CONTROL = "self_share"  # d / (d + TAU): the share's nonlinearity with no partner
BOOTSTRAPS = 2000
SEED = 0
# two variants are tested, so each gain's interval is two-sided at 97.5% (Bonferroni over two)
INTERVAL = (0.0125, 0.9875)
RIDGE = 1e-3  # crispri.logistic_fit's own default
CELLS = crispri.MODEL_CELLS  # the cells the sweep scored on their own track
ALL_ELEMENTS = crispri.ELEMENTS
PASS_RULE = (
    "a coupling variant passes if the lower end of its 97.5% component-bootstrap interval of AUPRC "
    "gain over the independent score, on out-of-fold scores pooled over the leave-one-chromosome-out "
    "folds, lies above 0, and its point gain exceeds the self-share control's point gain. R8 proceeds "
    "if at least one variant passes"
)
FALSIFIER = (
    "a pass is void if the self-share control, which has the share's nonlinearity but no partner, "
    "gains at least as much as the variant: the gain is then a reshaping of d, not coupling. A "
    "positive that holds on components but whose chromosome-bootstrap interval includes 0 is reported "
    "as a weak pass and names the chromosome that carries it"
)
READINGS = {
    "neither_passes": (
        "a negative: coupling carries no information about CRISPRi targets beyond independent "
        "per-element scoring on the data held, so a joint search over this score cannot beat per-block "
        "scoring. R8 closes; no search engine is built"
    ),
    "one_passes": (
        "the passing coupling term is kept, frozen at the form registered here, and R8 proceeds to its "
        "build (fixed-boundary search on a synthetic two-change case first); the other term is dropped"
    ),
    "both_pass": "both terms are kept, frozen as registered, and R8 proceeds to its build",
    "void": "the falsifier fired: reported as a negative with the control's gain beside it; R8 closes",
}

POSITIVE, NEGATIVE = measured.DECREASE, measured.NULL_INFORMATIVE


def label_of(outcome: str) -> bool | None:
    """R2's reading: True for a significant decrease, False for a well-powered null, else not scored."""
    if outcome == POSITIVE:
        return True
    if outcome == NEGATIVE:
        return False
    return None


@dataclass
class Row:
    """One training pair and what the cache predicted for it; `drop` None means missing."""

    chrom: str
    start: int
    end: int
    gene: str
    cell: str
    outcome: str
    drop: float | None = None
    reason: str = ""
    features: dict[str, float] = field(default_factory=dict)

    @property
    def label(self) -> bool | None:
        return label_of(self.outcome)

    @property
    def element(self) -> tuple[str, int, int, str]:
        return (self.chrom, self.start, self.end, self.cell)

    @property
    def gene_key(self) -> tuple[str, str]:
        return (self.gene, self.cell)


MISSING_CELL = "the cell was not scored on its own track by the sweep"
MISSING_ELEMENT = "no cached element overlaps the tested interval"
MISSING_GENE = "no overlapping element scored this gene in this cell"


def load_training() -> list[measured.CrispriPair]:
    """The training split only; `load_crispri` skips the held-out file before opening it."""
    pairs, _invalid = measured.load_crispri(split=measured.TRAINING)
    return measured.development_only(pairs)


def _elements(root: Path, chrom: str) -> tuple[list[dict[str, Any]], list[int]]:
    p = root / f"{chrom}.json"
    els = sorted(json.loads(p.read_text()), key=lambda e: e["start"]) if p.exists() else []
    return [{"id": e["id"], "start": e["start"], "end": e["end"]} for e in els], [e["start"] for e in els]


def _overlapping(els: list[dict[str, Any]], starts: list[int], start: int, end: int) -> list[dict[str, Any]]:
    out = []
    j = bisect.bisect_left(starts, end) - 1
    while j >= 0 and els[j]["start"] > start - crispri.REACH:
        if els[j]["end"] > start:
            out.append(els[j])
        j -= 1
    return out


def score(
    pairs: list[measured.CrispriPair],
    cache: ElementResponses | None = None,
    elements_root: Path = ALL_ELEMENTS,
) -> list[Row]:
    """Each pair with its predicted drop, one chromosome at a time (one cache reader, one archive)."""
    measured.development_only(pairs)
    cache = cache if cache is not None else ElementResponses(ELEMENT_CACHE)
    by_chrom: dict[str, list[measured.CrispriPair]] = defaultdict(list)
    for p in pairs:
        by_chrom[p.chrom].append(p)
    rows: list[Row] = []
    for chrom in sorted(by_chrom):
        els, starts = _elements(elements_root, chrom)
        for p in by_chrom[chrom]:
            r = Row(p.chrom, p.start, p.end, p.gene, p.cell, p.outcome)
            rows.append(r)
            if p.cell not in CELLS:
                r.reason = MISSING_CELL
                continue
            hits = _overlapping(els, starts, p.start, p.end)
            if not hits:
                r.reason = MISSING_ELEMENT
                continue
            values = [v for e in hits if (v := cache.value(chrom, e["id"], p.gene, p.cell)) is not None]
            if not values:
                r.reason = MISSING_GENE
                continue
            r.drop = crispri.deletion_drop(values)
    return rows


def couple(rows: list[Row], tau: float = TAU) -> None:
    """Fill each scored row's shares. Partners are every scored row of the same perturbation (element
    competition) or of the same gene (gene budget), whatever their outcome: the partners' predictions
    are read, never their measured labels."""
    by_el: dict[tuple, float] = defaultdict(float)
    by_gene: dict[tuple, float] = defaultdict(float)
    n_el: Counter = Counter()
    n_gene: Counter = Counter()
    for r in rows:
        if r.drop is None:
            continue
        by_el[r.element] += r.drop
        by_gene[r.gene_key] += r.drop
        n_el[r.element] += 1
        n_gene[r.gene_key] += 1
    for r in rows:
        if r.drop is None:
            continue
        d = r.drop
        r.features = {
            "drop": d,
            "element_competition": d / (by_el[r.element] + tau),
            "gene_budget": d / (by_gene[r.gene_key] + tau),
            "self_share": d / (d + tau),
            "element_partners": float(n_el[r.element] - 1),
            "gene_partners": float(n_gene[r.gene_key] - 1),
        }


def components(rows: list[Row]) -> list[int]:
    """Connected component of each row in the perturbation-gene graph (an edge per tested pair)."""
    parent: dict[Any, Any] = {}

    def find(x: Any) -> Any:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for r in rows:
        a, b = find(("e", r.element)), find(("g", r.gene_key))
        if a != b:
            parent[a] = b
    ids: dict[Any, int] = {}
    return [ids.setdefault(find(("e", r.element)), len(ids)) for r in rows]


def out_of_fold(rows: list[Row], names: tuple[str, ...], lam: float = RIDGE) -> list[float]:
    """Logistic scores for each row from a model fitted on every other chromosome."""
    out = [0.0] * len(rows)
    for chrom in sorted({r.chrom for r in rows}):
        fit = [r for r in rows if r.chrom != chrom]
        w = crispri.logistic_fit(
            [[r.features[n] for n in names] for r in fit], [bool(r.label) for r in fit], lam
        )
        for i, r in enumerate(rows):
            if r.chrom == chrom:
                out[i] = w[0] + sum(a * r.features[n] for a, n in zip(w[1:], names, strict=True))
    return out


def gain_interval(
    a: list[float],
    b: list[float],
    labels: list[bool],
    groups: list[Any],
    n: int = BOOTSTRAPS,
    seed: int = SEED,
    q: tuple[float, float] = INTERVAL,
) -> dict[str, Any]:
    """AUPRC(a) - AUPRC(b) with an interval from resampling whole groups (components or chromosomes)."""
    ap = crispri.average_precision
    point = (ap(a, labels) or 0.0) - (ap(b, labels) or 0.0)
    members: dict[Any, list[int]] = defaultdict(list)
    for i, g in enumerate(groups):
        members[g].append(i)
    keys = sorted(members, key=str)
    rng = random.Random(seed)
    diffs = []
    for _ in range(n):
        idx = [i for k in (rng.choice(keys) for _ in keys) for i in members[k]]
        lab = [labels[i] for i in idx]
        if not any(lab) or all(lab):
            continue
        diffs.append((ap([a[i] for i in idx], lab) or 0.0) - (ap([b[i] for i in idx], lab) or 0.0))
    diffs.sort()
    if not diffs:
        return {"gain": round(point, 4), "interval": None, "resamples": 0, "groups": len(keys)}
    lo = diffs[int(q[0] * len(diffs))]
    hi = diffs[min(len(diffs) - 1, int(q[1] * len(diffs)))]
    return {
        "gain": round(point, 4),
        "interval": [round(lo, 4), round(hi, 4)],
        "resamples": len(diffs),
        "groups": len(keys),
    }


def evaluate(rows: list[Row], n: int = BOOTSTRAPS, seed: int = SEED) -> dict[str, Any]:
    """The registered comparison on the scored, labelled rows, pooled over the chromosome folds."""
    couple(rows)
    comp_all = components(rows)
    ev = [i for i, r in enumerate(rows) if r.drop is not None and r.label is not None]
    er = [rows[i] for i in ev]
    labels = [bool(r.label) for r in er]
    comp = [comp_all[i] for i in ev]
    chroms = [r.chrom for r in er]
    indep = [r.features["drop"] for r in er]
    ap = crispri.average_precision
    scores = {"independent": indep}
    for name in (*VARIANTS, CONTROL):
        scores[name] = out_of_fold(er, ("drop", name))
    auprc = {k: round(ap(v, labels) or 0.0, 4) for k, v in scores.items()}
    auprc["share_alone_element_competition"] = round(
        ap([r.features["element_competition"] for r in er], labels) or 0.0, 4
    )
    auprc["share_alone_gene_budget"] = round(ap([r.features["gene_budget"] for r in er], labels) or 0.0, 4)
    gains = {}
    for name in (*VARIANTS, CONTROL):
        gains[name] = {
            "by_component": gain_interval(scores[name], indep, labels, comp, n, seed),
            "by_chromosome": gain_interval(scores[name], indep, labels, chroms, n, seed),
        }
    # secondary: only rows whose coupling is non-trivial (at least one scored partner)
    secondary = {}
    for name, part in (("element_competition", "element_partners"), ("gene_budget", "gene_partners")):
        sub = [i for i, r in enumerate(er) if r.features[part] >= 1]
        sl = [labels[i] for i in sub]
        secondary[name] = {
            "pairs": len(sub),
            "positives": sum(sl),
            "auprc_independent": round(ap([indep[i] for i in sub], sl) or 0.0, 4) if any(sl) else None,
            "auprc_coupled": round(ap([scores[name][i] for i in sub], sl) or 0.0, 4) if any(sl) else None,
            "gain_by_component": gain_interval(
                [scores[name][i] for i in sub], [indep[i] for i in sub], sl, [comp[i] for i in sub], n, seed
            )
            if any(sl)
            else None,
        }
    verdict = {}
    ctrl = gains[CONTROL]["by_component"]["gain"]
    for name in VARIANTS:
        g = gains[name]["by_component"]
        iv = g["interval"]
        above = iv is not None and iv[0] > 0
        beats_ctrl = g["gain"] > ctrl
        chrom_iv = gains[name]["by_chromosome"]["interval"]
        verdict[name] = {
            "interval_above_zero": above,
            "beats_self_share_control": beats_ctrl,
            "passes": above and beats_ctrl,
            "weak": above and beats_ctrl and not (chrom_iv is not None and chrom_iv[0] > 0),
        }
    passed = [v for v in VARIANTS if verdict[v]["passes"]]
    any_above = any(verdict[v]["interval_above_zero"] for v in VARIANTS)
    reading = (
        "both_pass"
        if len(passed) == 2
        else "one_passes"
        if passed
        else "void"
        if any_above
        else "neither_passes"
    )
    return {
        "pairs_scored": len(er),
        "positives": sum(labels),
        "negatives": len(labels) - sum(labels),
        "components": len(set(comp)),
        "chromosomes": len(set(chroms)),
        "auprc": auprc,
        "gain_over_independent": gains,
        "secondary_with_partners": secondary,
        "verdict": verdict,
        "reading": reading,
        "reading_text": READINGS[reading],
        "r8": "proceeds" if passed else "closes",
    }


def counts(rows: list[Row]) -> dict[str, Any]:
    """What was scored, what was missing and why, and the outcomes kept out of the metric."""
    out: dict[str, Any] = {"pairs": len(rows), "by_cell": dict(Counter(r.cell for r in rows))}
    out["outcomes"] = dict(Counter(r.outcome for r in rows))
    out["missing"] = dict(Counter(r.reason for r in rows if r.drop is None))
    scored = [r for r in rows if r.drop is not None]
    out["scored"] = len(scored)
    out["scored_outcomes"] = dict(Counter(r.outcome for r in scored))
    out["excluded_from_metric"] = {
        "significant_increase": sum(r.outcome == measured.INCREASE for r in scored),
        "not_significant_underpowered": sum(r.outcome == measured.NULL_INCONCLUSIVE for r in scored),
        "missing_outcome": sum(r.outcome == measured.MISSING for r in scored),
    }
    return out


def per_cell(
    rows: list[Row], pooled: dict[str, Any] | None = None, n: int = BOOTSTRAPS, seed: int = SEED
) -> dict[str, Any]:
    """The comparison per cell where the training file holds that cell with at least
    measured.MIN_FOR_A_COMPARISON positives; smaller cells are counted and not computed. A cell that
    is every row is the pooled result itself, which is passed in rather than computed twice."""
    out: dict[str, Any] = {}
    cells = {r.cell for r in rows}
    if pooled is not None and len(cells) == 1:
        return {next(iter(cells)): {"computed": True, "same_as_pooled": True}}
    for cell in sorted({r.cell for r in rows}):
        sub = [r for r in rows if r.cell == cell]
        pos = sum(1 for r in sub if r.drop is not None and r.label is True)
        if pos < measured.MIN_FOR_A_COMPARISON:
            out[cell] = {"positives": pos, "computed": False}
            continue
        out[cell] = {"computed": True, **evaluate(sub, n, seed)}
    return out
