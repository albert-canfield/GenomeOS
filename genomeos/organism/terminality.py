# SPDX-License-Identifier: AGPL-3.0-or-later
"""The terminality gate: can a rule set say WHEN a cell is done dividing (area E).

Area E's worm program scores 522 of 555 terminal fates, but 522 is the program's score and not the
factors' knowledge. Credited only where a rule decided the type it is 184 of 555 with nothing fitted,
381 in sample and 374 held out with the fitted neuron list, 136 free, and **0 fates decided that the
lineage had not already written, 414 of 414** (`542d613`, `d5215bf`). That last zero is fixed by a
guard, not by a read: every generated fate rule is written `when: cell_type = <one of the 17 terminal
types>`, and only the observed-lineage lookup ever puts a terminal type on a cell. So a factor rule in
this program chooses *which* of 17 types a cell that has already been called terminal becomes. It never
says *that* the cell is terminal. `d5215bf` closed the read question and named this as the only thing
left that a rule change can still move.

This module asks the gate question directly, off the runtime and against the atlas: **given what a cell
can measure about itself, can a rule call terminal against dividing?** If it cannot, that is not a
failure of this lane. It is a measurement of what the atlas contains, and it bounds every fate claim
area E can make: a fate rule set that cannot say when division stops is a classifier that runs only
after something else has told it which cells to classify.

============================================================================================
CIRCULARITY AUDIT, established from the code and the atlas before any rule was written
============================================================================================

A signal derived from the lineage lookup is the answer sheet in a second handwriting. This area has
been caught by that once already -- the 522 headline itself -- so every candidate signal is ruled on
first, and the ruling is part of the registration.

  1. `cell_type`  THE ANSWER RESTATED. EXCLUDED. Only `fate_<cell>` from the lookup ever sets a
     terminal type (`body.py` `_pick_fate`, `_commit`), so a guard on it is a guard on the lookup
     having already decided. This is exactly the guard that pins the count at 0 of 414.

  2. `cell` (the lineage name)  THE LOOKUP'S OWN ADDRESS. EXCLUDED AS A RULE TERM. `ABalaaaalal`
     spells out eleven divisions, and `when: cell = X` is the form of every row of the lookup itself,
     so a rule keyed on a name is a lookup row wearing a rule's syntax. The name is retained only as
     the address at which a measurement is looked up -- which is what `reader.bio` already does, and
     what the shipped fate rules already borrow -- and that borrowing is stated, not credited.

  3. Tracked lifetime (`_lifetimes` in the Ma 2021 table)  CIRCULAR. EXCLUDED. A tracked cell's
     lifetime ends either when it divides or when the movie stops. The imaging stops at frame 320 =
     400.0 min; 454 cells end at that frame, 322 of them terminal against 78 dividing. A tracked
     lifetime shorter than the run is therefore a division that was seen to happen, and reading it
     would be reading the answer. This is the trap this audit exists to catch, and it is the one that
     would have looked most like a measurement.

  4. The cell's own cycle length (`divides - born`)  THE ANSWER ITSELF. EXCLUDED. It is defined only
     for a cell that divides.

  5. Atlas coverage -- whether a cell is named in Ma 2021 at all  LINEAGE-DERIVED. EXCLUDED AS A RULE
     TERM. 455 of 555 terminal cells are named (82.0%) against 717 of 771 dividing cells (93.0%), so
     "absent from the atlas" carries eleven points of signal about the class. That is a fact about
     which cells the imaging resolved, not about the cell, and the imaging stops before the lineage
     does (deepest atlas-named cell generation 10 and born 616.0 min; deepest terminal cell generation
     11 and born 785.0 min). No arm branches on coverage: a cell the atlas does not name inherits the
     factor set of its nearest tracked ancestor, which is what the reader gives it in the Body today
     and which, if anything, works against the factors.

  6. `generation` (lineage depth)  PERMITTED, BUT LOOKUP-DERIVED, SO A BASELINE AND NEVER A RESULT.
     It is the cell's own division history and it is not the answer restated: generation 9 holds 283
     terminal cells and 46 dividing ones. But in this program every division is scheduled by a
     `div_<cell>` decision of the lookup, so a rule's `generation` is a number the lookup produced.
     Depth arms are registered below as baselines to beat.

  7. Birth time, stage, and the cycle timers  SAME STATUS AS DEPTH, and for the same reason: the
     timers are Sulston cycle means per (founder, generation), distilled from the same tree.

  8. `lineage` (the founder)  LOOKUP-DERIVED. The founders are named by the `lineages:` clause of the
     lookup's own division decisions. Used only to define held-out folds and the per-sublineage depth
     baseline, which is the most generous baseline registered here.

  9. Measured transcription factors (Ma et al. 2021, presence at >= 20% of the factor's maximum;
     instantaneous or integrated along the lineage path)  A MEASUREMENT OF THE CELL. This is the only
     candidate in the list that is not a restatement of the tree, and it is what the arms are for.

 10. Factor COUNT, identities erased  A MEASUREMENT, BUT A CLOCK. Terminal cells carry more measured
     factors (74.2 against 46.9 among atlas-named cells) and are born later (mean 468.7 min against
     308.5). A rule set cannot be credited for what a bare count of factors already gives, so the
     count is registered as its own arm.

============================================================================================
PRE-REGISTRATION.  Committed before any factor arm was run and before the result file existed.
============================================================================================

THE UNIVERSE.  The 1,326 embryonic cells of the reference lineage (born < 800 min) that are not
programmed deaths: 555 terminal and 771 dividing. The 113 embryonic deaths are a third class and are
reported separately, in a secondary universe of 1,439 where a death counts as not-dividing.

THE METRIC.  Balanced accuracy -- the mean of sensitivity on the 555 terminal cells and specificity on
the 771 dividing ones -- with **the lookup's terminality deleted**: no arm may read `cell_type`,
children, division time or tracked lifetime. **Every cell is called.** An arm that abstains has its
abstentions called `Dividing`, the majority class, and those calls are scored. This is registered
deliberately: the cited fate rules today claim more non-terminal cells than terminal ones (105 against
79, `542d613`), so a rule set that merely abstains more would look better on any claimed-cell metric.
It cannot look better on this one. Accuracy, Matthews correlation and the claimed count are published
beside the balanced accuracy for every arm, so abstention stays visible without being rewarded.

THE BASELINES, computed before any factor arm was run:

    T1  always-dividing                       balanced 0.5000   accuracy 0.5814
    T2  always-terminal                       balanced 0.5000   accuracy 0.4186
    T3  best single depth threshold, gen >= 9 balanced 0.7815   accuracy 0.8069   (in sample)
    T4  a depth threshold per founder         balanced 0.8998   accuracy 0.8967   (in sample)
        ABal 9, ABar 9, ABpl 9, ABpr 9, MS 6, E 5, C 6, D 4, other 1

T3 and T4 are fitted on the whole universe, which is generous to them on purpose: a depth rule is the
thing a factor rule has to be worth more than, and it should be met at its strongest. Held-out forms
(T3cv, T4cv, one founder sublineage out at a time) are reported too.

WHAT WOULD COUNT AS A RULE SET THAT KNOWS SOMETHING.  All three clauses, held out by founder
sublineage:

    (a) the best factor-ONLY arm reaches max(T1..T4cv) balanced accuracy + 0.05;
    (b) the best factor-PLUS-DEPTH arm beats the best depth-only arm by at least 0.03 balanced
        accuracy -- the factors must add something to the clock rather than restate it;
    (c) pooled WITHIN-GENERATION-BAND balanced accuracy reaches 0.55, with at least two individual
        bands at 0.55 or better. Inside a band depth is constant, so this clause is the one that
        cannot be passed by a clock.

WHAT WOULD COUNT AS FAILURE.  (a) and (b) both failing. Registered in advance as the reportable
result: it would mean nothing measured in this atlas distinguishes a cell that will divide again from
one that will not, and that every fate number area E publishes is conditional on a terminality call it
did not make.

A LOW NUMBER IS NOT A FAILURE OF THIS LANE. BROKEN ACCOUNTING IS.

PREDICTIONS, all open when this was committed:

    P1  every factor-only arm (learned rules on either read, and the logistic ceiling) is BELOW T3cv.
    P2  clause (b) fails: factors plus depth beat depth alone by less than 0.03 balanced accuracy.
    P3  the pooled within-band figure is below 0.60.
    P4  the shuffle null -- the same learner on factor sets permuted among the 1,326 cells -- sits
        within 0.03 of 0.500 balanced accuracy.
    P5  the identity-free factor-COUNT threshold comes within 0.03 of the best learned factor rule
        set, i.e. what the rules read is the clock and not the identity of any factor.
    P6  the logistic ceiling on factors alone does not reach T3cv either, so the shortfall is in the
        atlas and not in the rule language.

VOID CONDITIONS.  The run stops and reports nothing if (V1) the universe is not 555 terminal + 771
dividing = 1,326, or (V2) any arm's called count is not 1,326.

HYPERPARAMETERS, fixed in advance and not tuned: the learner is `fate_rules.learn` unchanged --
the same greedy decision list that produced the shipped fate rules, pointed at terminality instead of
type -- with `min_support=3`, `min_precision=0.8`, `max_rules=40`. The logistic ceiling is L2 with
`lam=1.0`, 2,000 full-batch gradient steps at `1/L` for the loss's own Lipschitz constant, so it
converges monotonically and has no step to damp; its decision threshold is chosen on the TRAINING fold
to maximise training balanced accuracy and is never touched on the held-out fold.

0 model requests. Nothing is fetched.

============================================================================================
ADDENDUM, written after the registered arms were run. The registration above keeps every number
it was committed with (`de7a135`); nothing in it is edited.
============================================================================================

CLAUSE (c) WAS REGISTERED ON THE WRONG STATISTIC, AND THIS IS THE CORRECTION.  "Pooled
within-generation-band balanced accuracy" was meant to be the clause a clock cannot pass, because
depth is constant inside a band. Pooling the four counts across bands does not do that. Strata with
different class priors, pooled, let a predictor that names **each stratum's own majority** -- which
reads nothing but the stratum, that is nothing but depth and founder -- score far above 0.5 while
scoring exactly 0.500 inside every single stratum. Measured, on these strata:

    stratum-majority caller, within band            pooled 0.7575   macro 0.5000
    stratum-majority caller, within band + founder  pooled 0.7886   macro 0.5000

So the registered clause (c) figure was mostly the prior it was meant to remove. The corrected
statistic is the MACRO AVERAGE of the per-stratum balanced accuracies, whose floor is exactly 0.500
by construction and is what the stratum-majority baseline scores. Both are published, the registered
one first, and the registered clause (c) verdict stands as computed as well as being corrected.

THREE POST-HOC CIRCULARITY CONTROLS on the within-band result, each labelled as post-hoc:

    D1  no per-band arm selection: one arm fixed across all bands, so the number is not the best of
        four chosen on the held-out fold.
    D2  atlas-named cells only, so a cell that inherits a distant ancestor's factor set cannot let
        the excluded coverage signal back in through the read.
    D3  depth AND founder sublineage both held constant, with the folds taken by GRANDPARENT so a
        cell is never tested on a model that saw its sister. This is the strictest form available:
        within one stratum both lookup-derived variables are fixed and the only thing left varying
        between the cells is what was measured in them.

0 model requests in the addendum either.
"""

from __future__ import annotations

import random
from collections import Counter, defaultdict
from dataclasses import dataclass

from .fate_rules import Rule, apply_rules, learn, sublineage
from .reference import HATCH_MIN, ReferenceLineage

TERMINAL = "Terminal"
DIVIDING = "Dividing"
DEFAULT_CALL = DIVIDING  # registered: an abstention is called for the majority class and is scored

MIN_SUPPORT = 3
MIN_PRECISION = 0.8
MAX_RULES = 40
LOGISTIC_LAMBDA = 1.0
LOGISTIC_STEPS = 2000


# ---- the universe -----------------------------------------------------------------------------


def universe(ref: ReferenceLineage, deaths: bool = False) -> dict[str, str]:
    """cell -> Terminal or Dividing for every embryonic cell, with the lookup's terminality kept only
    here, as the label. `deaths=False` drops the 113 programmed deaths; `deaths=True` counts them as
    Terminal (they are done dividing, but they are not a fate)."""
    out: dict[str, str] = {}
    for c in ref.cells.values():
        if c.born >= HATCH_MIN:
            continue
        if c.dies is not None:
            if deaths:
                out[c.id] = TERMINAL
            continue
        out[c.id] = TERMINAL if not c.children else DIVIDING
    return out


def depths(ref: ReferenceLineage, cells) -> dict[str, int]:
    return {cid: ref.cells[cid].generation for cid in cells}


# ---- scoring: every cell is called -------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Score:
    arm: str
    balanced: float
    accuracy: float
    mcc: float
    sensitivity: float
    specificity: float
    claimed: int  # cells an arm actually decided; the rest were called by DEFAULT_CALL
    called_terminal: int
    called_dividing: int
    n: int

    def row(self) -> dict:
        return {
            "arm": self.arm,
            "balanced": round(self.balanced, 4),
            "accuracy": round(self.accuracy, 4),
            "mcc": round(self.mcc, 4),
            "sensitivity": round(self.sensitivity, 4),
            "specificity": round(self.specificity, 4),
            "claimed": self.claimed,
            "called_terminal": self.called_terminal,
            "called_dividing": self.called_dividing,
            "n": self.n,
        }


def score(arm: str, pred: dict[str, str], labels: dict[str, str]) -> Score:
    """`pred` may abstain on a cell (absent, or None): the abstention is called DEFAULT_CALL and is
    scored. So an arm cannot raise this metric by claiming fewer cells."""
    tp = tn = fp = fn = 0
    claimed = 0
    for cid, truth in labels.items():
        call = pred.get(cid) or DEFAULT_CALL
        if pred.get(cid):
            claimed += 1
        if call == TERMINAL:
            if truth == TERMINAL:
                tp += 1
            else:
                fp += 1
        elif truth == TERMINAL:
            fn += 1
        else:
            tn += 1
    n = len(labels)
    sens = tp / (tp + fn) if tp + fn else 0.0
    spec = tn / (tn + fp) if tn + fp else 0.0
    den = ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) ** 0.5
    mcc = ((tp * tn - fp * fn) / den) if den else 0.0
    return Score(
        arm, (sens + spec) / 2, (tp + tn) / n if n else 0.0, mcc, sens, spec, claimed, tp + fp, tn + fn, n
    )


# ---- depth baselines ---------------------------------------------------------------------------


def _threshold_balanced(gen: dict[str, int], labels: dict[str, str], k: int) -> float:
    return score("t", {cid: (TERMINAL if g >= k else DIVIDING) for cid, g in gen.items()}, labels).balanced


def best_threshold(gen: dict[str, int], labels: dict[str, str], hi: int = 13) -> int:
    """The depth threshold with the best balanced accuracy on the cells given."""
    return max(range(hi + 1), key=lambda k: (_threshold_balanced(gen, labels, k), -k))


def depth_call(gen: dict[str, int], k: int) -> dict[str, str]:
    return {cid: (TERMINAL if g >= k else DIVIDING) for cid, g in gen.items()}


def per_sublineage_thresholds(gen: dict[str, int], labels: dict[str, str]) -> dict[str, int]:
    groups: dict[str, list[str]] = defaultdict(list)
    for cid in labels:
        groups[sublineage(cid)].append(cid)
    return {
        s: best_threshold({c: gen[c] for c in cs}, {c: labels[c] for c in cs}) for s, cs in groups.items()
    }


def per_sublineage_call(gen: dict[str, int], ks: dict[str, int], fallback: int) -> dict[str, str]:
    return {cid: (TERMINAL if g >= ks.get(sublineage(cid), fallback) else DIVIDING) for cid, g in gen.items()}


def count_call(counts: dict[str, int], k: int) -> dict[str, str]:
    return {cid: (TERMINAL if n >= k else DIVIDING) for cid, n in counts.items()}


def best_count_threshold(counts: dict[str, int], labels: dict[str, str]) -> int:
    hi = max(counts.values(), default=0)
    return max(range(hi + 2), key=lambda k: (score("c", count_call(counts, k), labels).balanced, -k))


# ---- learned rules, held out by founder sublineage ----------------------------------------------


def depth_tokens(gen: int) -> set[str]:
    """Depth offered to the rule learner in the only form its language has: presence tokens. Marked
    `gen>=k` so a rule that uses one is visible as a depth rule in the written rule set."""
    return {f"gen>={k}" for k in range(1, gen + 1)}


def with_depth(feats: dict[str, set[str]], gen: dict[str, int]) -> dict[str, set[str]]:
    return {cid: set(fs) | depth_tokens(gen[cid]) for cid, fs in feats.items()}


def learn_rules(feats: dict[str, set[str]], labels: dict[str, str]) -> list[Rule]:
    return learn(feats, labels, min_support=MIN_SUPPORT, min_precision=MIN_PRECISION, max_rules=MAX_RULES)


def held_out_rules(feats: dict[str, set[str]], labels: dict[str, str]) -> dict[str, str]:
    """One founder sublineage out at a time: rules learned on the other seven call the held-out cells."""
    pred: dict[str, str] = {}
    groups: dict[str, list[str]] = defaultdict(list)
    for cid in labels:
        groups[sublineage(cid)].append(cid)
    for g, held in groups.items():
        train = {c: feats[c] for c in labels if sublineage(c) != g}
        rules = learn_rules(train, {c: labels[c] for c in train})
        pred.update(apply_rules(rules, {c: feats[c] for c in held}))
    return pred


def in_sample_rules(feats: dict[str, set[str]], labels: dict[str, str]) -> tuple[dict[str, str], list[Rule]]:
    rules = learn_rules({c: feats[c] for c in labels}, labels)
    return apply_rules(rules, {c: feats[c] for c in labels}), rules


def shuffle_null(feats: dict[str, set[str]], labels: dict[str, str], n: int = 20, seed: int = 0) -> dict:
    """The same learner on factor sets permuted among the cells: what the machinery scores when the
    measurement has been detached from the cell it was measured in."""
    rng = random.Random(seed)
    cells = sorted(labels)
    out = []
    for _ in range(n):
        order = cells[:]
        rng.shuffle(order)
        permuted = {c: feats[o] for c, o in zip(cells, order, strict=True)}
        out.append(score("null", held_out_rules(permuted, labels), labels).balanced)
    mean = sum(out) / len(out)
    sd = (sum((x - mean) ** 2 for x in out) / len(out)) ** 0.5
    return {"draws": n, "balanced_mean": round(mean, 4), "balanced_sd": round(sd, 4)}


# ---- the ceiling: what ANY function of these factors could know ---------------------------------


def _fit_logistic(
    rows: list[list[float]], y: list[int], lam: float = LOGISTIC_LAMBDA, steps: int = LOGISTIC_STEPS
):
    """L2 logistic regression by full-batch gradient descent at 1/L for the loss's own Lipschitz
    constant, so the step is never larger than the curvature allows and there is nothing to damp."""
    import numpy as np

    x = np.asarray(rows, dtype=float)
    t = np.asarray(y, dtype=float)
    n, d = x.shape if x.ndim == 2 else (0, 0)
    if n == 0 or d == 0:
        return np.zeros(max(d, 1)), 0.0
    w = np.zeros(d)
    b = 0.0
    sigma = float(np.linalg.norm(x, 2))
    lip = 0.25 * (sigma * sigma + n) / n + lam  # +n for the intercept column of ones
    step = 1.0 / lip
    for _ in range(steps):
        p = 1.0 / (1.0 + np.exp(-(x @ w + b)))
        w -= step * (x.T @ (p - t) / n + lam * w)
        b -= step * float((p - t).mean())
    return w, b


def _predict_logistic(w, b, rows: list[list[float]]):
    import numpy as np

    x = np.asarray(rows, dtype=float)
    if x.size == 0:
        return np.zeros(0)
    return 1.0 / (1.0 + np.exp(-(x @ w + b)))


def _best_cut(p, y: list[int]) -> float:
    pos = sum(y)
    neg = len(y) - pos
    if not pos or not neg:
        return 0.5
    best, cut = -1.0, 0.5
    for c in sorted({round(float(x), 4) for x in p} | {0.5}):
        tp = sum(1 for a, b in zip(p, y, strict=True) if a >= c and b == 1)
        tn = sum(1 for a, b in zip(p, y, strict=True) if a < c and b == 0)
        v = (tp / pos + tn / neg) / 2
        if v > best:
            best, cut = v, c
    return cut


def logistic_held_out(
    feats: dict[str, set[str]],
    labels: dict[str, str],
    extra: dict[str, list[float]] | None = None,
    groups: dict[str, str] | None = None,
) -> dict[str, str]:
    """Leave one founder sublineage out. The decision threshold is chosen on the TRAINING fold to
    maximise training balanced accuracy and applied unchanged to the held-out fold."""
    by: dict[str, list[str]] = defaultdict(list)
    for cid in labels:
        by[(groups or {}).get(cid) or sublineage(cid)].append(cid)
    pred: dict[str, str] = {}
    for held in by.values():
        train = [c for c in labels if c not in set(held)]
        if not train or len({labels[c] for c in train}) < 2:
            continue
        cols = sorted({f for c in train for f in feats[c]})
        idx = {f: i for i, f in enumerate(cols)}
        ncols = len(cols)
        width = ncols + (len(extra[train[0]]) if extra else 0)

        def row(c: str, idx=idx, ncols=ncols, width=width) -> list[float]:
            v = [0.0] * width
            for f in feats[c]:
                if f in idx:
                    v[idx[f]] = 1.0
            if extra:
                for j, x in enumerate(extra[c]):
                    v[ncols + j] = x
            return v

        xtr = [row(c) for c in train]
        ytr = [1 if labels[c] == TERMINAL else 0 for c in train]
        w, b = _fit_logistic(xtr, ytr)
        cut = _best_cut(_predict_logistic(w, b, xtr), ytr)
        for c, p in zip(held, _predict_logistic(w, b, [row(c) for c in held]), strict=True):
            pred[c] = TERMINAL if p >= cut else DIVIDING
    return pred


# ---- within a generation band: depth held constant ----------------------------------------------


def bands(gen: dict[str, int], labels: dict[str, str], min_each: int = 20) -> list[int]:
    """Generation bands with at least `min_each` cells of each class: the only places where a
    terminality call can be made without the clock making it."""
    per: dict[int, Counter] = defaultdict(Counter)
    for cid, lab in labels.items():
        per[gen[cid]][lab] += 1
    return sorted(g for g, c in per.items() if c[TERMINAL] >= min_each and c[DIVIDING] >= min_each)


def strata(
    gen: dict[str, int],
    labels: dict[str, str],
    founder: bool = False,
    min_cells: int = 40,
    min_each: int = 10,
) -> list[tuple[str, list[str]]]:
    """Cells grouped by generation, and optionally by founder sublineage too. Inside one stratum both
    lookup-derived variables are constant, so nothing a clock or a founder knows can call a cell."""
    by: dict[tuple, list[str]] = defaultdict(list)
    for cid in labels:
        by[(gen[cid], sublineage(cid) if founder else "")].append(cid)
    out = []
    for key, cs in sorted(by.items()):
        n_t = sum(1 for c in cs if labels[c] == TERMINAL)
        if len(cs) >= min_cells and n_t >= min_each and len(cs) - n_t >= min_each:
            out.append((f"{key[0]}/{key[1]}" if founder else str(key[0]), cs))
    return out


def stratum_majority(cells: list[str], labels: dict[str, str]) -> dict[str, str]:
    """The baseline the macro average exists to expose: name the stratum's own majority. It reads
    nothing but the stratum -- that is, nothing but depth and founder -- and it scores exactly 0.500
    inside every stratum while scoring far above 0.5 when the strata are pooled."""
    n_t = sum(1 for c in cells if labels[c] == TERMINAL)
    return dict.fromkeys(cells, TERMINAL if n_t * 2 > len(cells) else DIVIDING)


def stratified(
    groups: list[tuple[str, list[str]]], labels: dict[str, str], predict
) -> tuple[dict[str, float], float, float]:
    """Per-stratum balanced accuracy (its floor is exactly 0.500), the MACRO average of those, and the
    naive pooled figure that the macro average replaces."""
    per: dict[str, float] = {}
    tot = [0, 0, 0, 0]
    for key, cells in groups:
        lab = {c: labels[c] for c in cells}
        if len(set(lab.values())) < 2:
            continue
        pred = predict(cells, lab)
        tp = fn = tn = fp = 0
        for c, truth in lab.items():
            call = pred.get(c) or DEFAULT_CALL
            if truth == TERMINAL:
                tp, fn = tp + (call == TERMINAL), fn + (call != TERMINAL)
            else:
                fp, tn = fp + (call == TERMINAL), tn + (call != TERMINAL)
        per[key] = (tp / (tp + fn) + tn / (tn + fp)) / 2
        tot = [tot[0] + tp, tot[1] + fn, tot[2] + tn, tot[3] + fp]
    macro = sum(per.values()) / len(per) if per else 0.0
    tp, fn, tn, fp = tot
    pooled = (tp / (tp + fn) + tn / (tn + fp)) / 2 if (tp + fn) and (tn + fp) else 0.0
    return per, macro, pooled


def grandparent_folds(parents: dict[str, str], cells: list[str], k: int = 5) -> dict[str, int]:
    """Fold assignment by grandparent, so a cell is never tested on a model that saw its sister. A
    parent link says nothing about whether this cell divides, so it is not the answer restated."""
    gp = {}
    for c in cells:
        p = parents.get(c) or c
        gp[c] = parents.get(p) or p
    order = sorted(set(gp.values()))
    at = {g: i % k for i, g in enumerate(order)}
    return {c: at[gp[c]] for c in cells}


def within_stratum(
    cells: list[str],
    labels: dict[str, str],
    feats: dict[str, set[str]],
    parents: dict[str, str],
    kind: str = "rules",
    k: int = 5,
) -> dict[str, str]:
    """One stratum, held out by grandparent fold. Depth and founder are already constant here, so what
    is left for a rule to read is what was measured in the cell."""
    fold = grandparent_folds(parents, cells, k)
    sub = {c: feats[c] for c in cells}
    lab = {c: labels[c] for c in cells}
    pred: dict[str, str] = {}
    for i in range(k):
        held = [c for c in cells if fold[c] == i]
        train = [c for c in cells if fold[c] != i]
        if not held or len({lab[c] for c in train}) < 2:
            continue
        if kind == "logistic":
            pred.update(
                logistic_held_out(sub, lab, groups={c: ("H" if fold[c] == i else "T") for c in cells})
            )
        else:
            rules = learn_rules({c: sub[c] for c in train}, {c: lab[c] for c in train})
            pred.update(apply_rules(rules, {c: sub[c] for c in held}))
    return pred
