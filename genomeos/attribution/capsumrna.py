# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""lane-capsumrna (ROADMAP area I, the follow-on 97ff71f names, 2026-10-10): the capped sum against RNA.

WHAT WAS ASKED. lane-strengthform (7d23bc3, 10ca354, fc328c5) ran the strength-dependent combination form
on held-out enhancer pairs (Martinez-Ara et al. 2024) and found one comparator the form did not beat: the
additive envelope, the interval the capped sum and `or` give on a pair. The follow-on is the capped sum
against measured RNA in the thirteen reader cells, as its own registration, on 7d's design (combine2.py,
917ce25 / 948cfc9) and 7h's (multrna.py, 5c620d3 / 90e5c86): 7d's population, endpoint, baseline,
control, interval and wording, with head-to-heads against the mean, competition and the multiplicative
form on the genes all four run.

THE FORM, AS THE ENGINE AND THE TWO LANES DEFINE IT. `runtime/grn.py` already carries it: `sum_capped`,
A = min(1, W), W the sum of the terms s_i*H(x_i) (`combine_activators`), one of the four rules lane-combine
compared on 2026-09-28 (`grn.ACTIVATOR_RULES`, `grn.COMBINATION_TEST`: "sum_capped max <= P <= S1+S2"),
and selectable for one run with `NetworkRuntime(..., activator_combination='sum_capped')`. Nothing in the
engine is changed here. lane-strengthform's `additive_envelope` is the interval it and `or` give on a pair
in boost units, max(b1, b2) <= boost <= log2(2**b1 + 2**b2 - 1); its upper end is the capped sum below
its cap. On RNA the engine's own units hold, and combine2's conventions (`CONVENTIONS`) fix them; no
number is added.

WHAT THE ALGEBRA SAYS BEFORE ANY COUNT (`STRUCTURE`). Under 7d's design an element the reader holds is a
strength-0 copy of its rule, kept, and the bridge's removal is a clamp at 0, kept. Under both the mean is
LINEAR in the element terms (combine2.STRUCTURE['mean_is_linear_under_the_hold']), and so is the capped
sum wherever its cap does not bind. At combine2's conventions the cap binds for no gene that has a fit
with a unique root: a gene the capped sum runs has, in every reader cell, the level the mean gives it.
Its RNA reading would be 7d's mean on a subset of 7d's genes, whose RNA 7d's run has read and committed
(data/results/combine2_result.json, genes_tables). Where the cap could bind (a second, capped root on
genes with four or more activating elements) the definition does not choose the root, and the root sets
the cap through the strength scale `bridge.DECLARED_STRENGTH` fixes by convention, which the bridge says
licenses no prediction but the single removal it was fitted to.

CORRECTED BEFORE ANY CENSUS WAS SAVED (2026-10-10; written as an addition because the commit guard
refuses removing 43dd4ac's lines). The paragraph above is wrong in one clause: "the cap binds for no
gene that has a fit with a unique root". The cap binds on no gene's UNCAPPED root, which is the mean on
every gene. But the capped root's larger element-driven share (q = Z > D) fits an inhibitor the mean's
cannot, so a gene the mean refuses with inhibitor_strength_above_1 can have the capped root ALONE: its
one fit, with a prediction the mean does not give (`STRUCTURE['a_capped_root_alone']`). The first full
census run, from 43dd4ac in a clean worktree, stopped at such a gene (SPMIP4) because
`difference_from_mean` raised there; nothing was written. DECISION's text is unchanged and reads such a
gene as a unique fit whose prediction is not the mean's, against its condition (1).

So this module holds the fit, the closed form, the run and the census, and the rule (`DECISION`) that
the census, not a registration, is the last step unless it finds a gene whose fit is unique and whose
prediction differs from the mean's. No RNA value is read here. No model request, no key, no money.
"""

from __future__ import annotations

import copy
import math
from dataclasses import dataclass, field
from typing import Any

from genomeos.attribution import bridge
from genomeos.attribution import combine2 as c2
from genomeos.ir import Module
from genomeos.runtime import grn

DATE = "2026-10-10"
LANE = "lane-capsumrna"
ROADMAP_ITEM = (
    "area I, the follow-on 97ff71f names: the capped-sum (additive envelope) combination against measured"
    " RNA in the thirteen reader cells, as its own registration, on 7d's and 7h's design"
)
#: the form under test, as runtime/grn.py names it (in ACTIVATOR_RULES since 2026-09-28; engine unchanged)
FORM = "sum_capped"
SELECTION = c2.SELECTION
CELLS = c2.CELLS
H = bridge.H_INTACT
INFORMATIVE_LOG2 = c2.INFORMATIVE_LOG2
CLOSED_FORM_TOLERANCE = c2.CLOSED_FORM_TOLERANCE
RUNTIME = c2.RUNTIME

FORM_SOURCE = (
    "runtime/grn.py combine_activators 'sum_capped': A = min(1, W), W = sum of s_i*H(x_i) over the gene's"
    " activating rules, a held rule counted at strength 0; inhibitors multiply, R = prod (1 - s_j*H(x_j))."
    " One of the four rules lane-combine compared (grn.ACTIVATOR_RULES; grn.COMBINATION_TEST, registered"
    " 2026-09-28: 'sum_capped max <= P <= S1+S2'), selectable for one run with"
    " NetworkRuntime(activator_combination='sum_capped'). lane-strengthform's 'additive_envelope'"
    " (strengthform.COMPARATORS) is the interval sum_capped and or give on an enhancer pair in boost units,"
    " max(b1, b2) <= boost <= log2(2**b1 + 2**b2 - 1), the upper end being the capped sum below its cap;"
    " this lane runs the capped sum itself, in the engine's units, and not or"
)
CONVENTIONS = (
    "combine2's, unchanged (FIT['mean']): every element at its intact state x = bridge.INTACT, so H(x) ="
    " bridge.H_INTACT = 0.5; the strongest activating element at bridge.DECLARED_STRENGTH (s = 1), the"
    " others in proportion to their drops; an inhibitor's s_j*h = e_j/(q + e_j). The cap is the engine's"
    " 1 in those units. No number is added"
)
FIT = {
    "principle": c2.FIT["principle"],
    FORM: (
        "d_i = 1 - RHO_i (activating), e_j = RHO_j - 1 (inhibiting), D = sum d_i, top = max d_i. The"
        " strongest at s = 1 fixes Z = V*R0/T = top/H and each term t_i = d_i/Z, so W = H*D/top at the"
        " intact state."
        " UNCAPPED ROOT (W <= 1): the cap binds in no cell, since a hold only lowers W; q = D, b = T(1 - q),"
        " the inhibitors and every s as combine2.fit('mean') gives them, and V = T*top/(H*R0) (the mean's V"
        " without its n). CAPPED ROOT (the cap binding at the intact state, W = 1 + c, c > 0): removing i"
        " lowers A from 1 to 1 - d_i/Z, so t_i = d_i/Z + c, and the strongest at s = 1 gives"
        " Z = ((k-1)*top - D)/(H*(k-1) - 1), c = (1 - D/Z)/(k-1), q = Z, k the number of activating"
        " elements. Level relative to T with the open set O: (1 - q) + Z*min(1, sum of t_i over open"
        " activators) * prod over closed inhibitors (1 + e_j/q)"
    ),
    "runnable": (
        "a gene is runnable under the capped sum when combine2.fit('mean') is (q < 1, every s in (0, 1]),"
        " its W <= 1 (+ bridge.STRENGTH_TOLERANCE), and no capped root exists beside the uncapped one;"
        " otherwise it is reported by reason: the mean's reason, sum_past_the_cap, or fit_not_unique"
        "; a gene with a capped root and no uncapped one is capped_root_only, its one fit the capped root"
    ),
}
STRUCTURE = {
    "identity_with_the_mean": (
        "the reader's hold is a strength-0 copy of the rule, kept, and the bridge's removal a clamp at 0,"
        " kept; under both the mean is linear in the element terms (closing element i lowers transcription"
        " by d_i of T whatever the others do), and so is the capped sum while W <= 1. On its uncapped root"
        " the capped sum therefore gives, in every reader cell, the mean's level, and its centred prediction"
        " is 7d's P_mean. Its population is a subset of the mean's (the mean's runnable genes with W <= 1)"
    ),
    "no_fit_past_the_cap": (
        "a capped root needs c > 0, i.e. Z > D, and with H*(k-1) > 1 that is (k-1)*top - D > D*(H*(k-1) - 1),"
        " i.e. D < top/H: the uncapped W = H*D/top below 1. With H*(k-1) <= 1 (k <= 3 at H = 0.5) no capped"
        " root has Z > 0 and c > 0 at all. So a gene whose uncapped W exceeds 1 has no fit at these"
        " conventions (sum_past_the_cap): the strongest term is H, and with the cap binding no removal can"
        " lower A by as much as its observation needs"
    ),
    "a_second_root": (
        "with four or more activating elements (H*(k-1) > 1) and D < top/H, the capped root exists beside"
        " the uncapped one when its q = Z < 1 and its inhibitors fit: both reproduce every removal and T at"
        " the same conventions, so the fit is not unique (fit_not_unique). The capped root's prediction"
        " differs from the mean's; it is counted and run, never chosen, because the definition does not"
        " choose between the roots"
    ),
    "a_capped_root_alone": (
        "the capped root's element-driven share q = Z exceeds the uncapped root's q = D, and an inhibitor's"
        " fitted s_j*h = e_j/(q + e_j) falls as q rises, so a gene the mean refuses with"
        " inhibitor_strength_above_1 can have a capped root and no uncapped one (capped_root_only). Its one"
        " fit is the capped root, and the mean has no prediction for it. A gene the mean refuses with"
        " basal_not_positive (D >= 1) has none, since Z > D"
    ),
    "what_the_bridge_licenses": (
        "where the cap binds, how far it binds is set by the strength scale, which bridge.DECLARED_STRENGTH"
        " fixes by convention: " + bridge.DECLARED_STRENGTH_LICENSES + ". At a smaller declared scale every"
        " gene the mean runs is uncapped, and the capped sum is the mean on every one of them"
    ),
    "what_is_assumed": c2.STRUCTURE["what_is_assumed"],
    "inhibitor_only": c2.STRUCTURE["inhibitor_only"],
    "network": c2.STRUCTURE["network"],
}
CENSUS = (
    "a count on the predictor side, before any RNA value is read (data/results/capsumrna_census.json), over"
    " the 24 compiled programs gated as 7d gated them: of the multi-element genes combine2.classify admits,"
    " how many the capped sum fits, by root, and why the rest are not fitted; 7d's mean and competition"
    " counts recomputed beside them and checked against data/results/combine2_census.json; on every gene"
    " with an uncapped root the largest difference, over the thirteen cells, between its closed-form"
    " prediction and the mean's; every informative gene run by the runtime under"
    " activator_combination='sum_capped' with combine2's HARNESS checks, on each root it has; the genes"
    " informative under all four forms (mean, competition, multiplicative, the capped sum's unique root)."
    " It reads no RNA value"
)
CENSUS_SENTENCE = (
    "of the {n} multi-element genes combine2.classify admits, the capped sum at combine2's conventions fits"
    " {fitted} on its uncapped root ({x} are not fitted: {reasons}); on every one of them its prediction in"
    " the thirteen reader cells is the runtime's mean's (largest difference {dmax} in log2), and the runtime"
    " under activator_combination='sum_capped' gives that closed form on the {ran} informative genes it ran"
    " ({miss} harness failures). {two} of the {fitted} also have a second, capped root, on which the"
    " capped sum differs from the mean ({two_inf} informative under it); {unique} have a unique fit, all"
    " with the mean's prediction, {unique_inf} of them informative and {common} informative under all four"
    " forms, before the RNA tables' row check"
)
#: CENSUS_SENTENCE above says every gene with a unique fit has the mean's prediction, which
#: STRUCTURE['a_capped_root_alone'] refutes; the census fills this one whenever a gene has the capped root
#: alone (scripts/capsumrna_census.py, `sentence`)
CENSUS_SENTENCE_CAPPED_ROOT_ALONE = (
    "of the {n} multi-element genes combine2.classify admits, the capped sum at combine2's conventions fits"
    " {fitted} on its uncapped root ({x} are not fitted: {reasons}); on every one of them its prediction in"
    " the thirteen reader cells is the runtime's mean's (largest difference {dmax} in log2), and the runtime"
    " under activator_combination='sum_capped' gives the closed form on every informative gene it ran, on"
    " either root ({ran} runs of a gene on a root, {miss} of them failing the harness). {two} of the"
    " {fitted} also have a second, capped root, on which the capped sum differs from the mean ({two_inf}"
    " informative under it). {only} genes have only a capped root ({only_inf} informative under it,"
    " {only_sym} of those with a symbol on one gene_id), the mean having no fit for them ({only_reasons})."
    " {unique} genes have a unique fit, {differ} of them with a prediction other than the mean's;"
    " {unique_inf} are informative and {common} informative under all four forms, before the RNA tables'"
    " row check"
)
DECISION = (
    "fixed before the census was counted, from STRUCTURE: no registration and no RNA run if (1) every gene"
    " with a unique fit has the mean's prediction in every reader cell within CLOSED_FORM_TOLERANCE, and"
    " (2) the runtime under activator_combination='sum_capped' gives the closed form on every informative"
    " gene run. Then 7d's run has already read the RNA of every gene the capped sum runs with a unique fit,"
    " under the same prediction (data/results/combine2_result.json, genes_tables), so a registration on"
    " 7d's design would re-read a subset of 7d's committed mean with its outcome in the record, and its"
    " head-to-head against the mean would be zero by construction. The capped root, where it exists, is a"
    " choice the definition does not make, set by a strength scale the bridge says licenses nothing beyond"
    " the single removal; it is counted for the coordinator's decision and not registered here. If either"
    " condition fails, the registration goes ahead on the genes where the two differ"
)
DECISIONS = ("no_registration_identity_with_the_mean", "register_on_the_genes_that_differ")
MIN_GENES = c2.MIN_GENES
#: what "on 7d's design" carries for the genes DECISION would register on: 7d's DECISION 2, copied
SIZE_CHECK = (
    "7d's DECISION 2, part of the design this lane was asked to register on: a population of fewer than"
    " MIN_GENES genes is NOT_TESTABLE, no comparison read. Applied to the genes whose unique fit differs from"
    " the mean's, informative, with a symbol on one gene_id (POPULATION (f)'s row check, which needs the RNA"
    " tables' rows, can only lower the count): below MIN_GENES no registration is written and no RNA table"
    " is opened. It was written after an unsaved census run showed the count (EXPLORATION_PRECEDED_THIS)"
)
#: what this lane read before the census was saved, stated rather than claimed away
EXPLORATION_PRECEDED_THIS = (
    "combine2.py, multrna.py, saturate.py, strengthform.py and their census scripts, and the per-form counts"
    " of data/results/combine2_census.json and multrna_census.json (predictor side, no RNA value); one"
    " unsaved smoke run of scripts/capsumrna_census.py on chr21, after DECISION and STRUCTURE were written"
    " and before this code was committed, which changed neither. No RNA table and no genes table of 7d or"
    " 7h was opened; their registered readings were known from the lane brief and docs/ROADMAP.md"
    ". After 43dd4ac: one full census run from it in a clean worktree, which stopped with an exception at"
    " SPMIP4 and wrote nothing (the correction in the module docstring); one unsaved full census run of the"
    " corrected code before it was committed, whose counts led to SIZE_CHECK being written (7d's DECISION"
    " 2, not a new rule)"
)


# ==== the fit (closed form; no RNA) =======================================================================


@dataclass
class CFit:
    """One gene's capped-sum fit on one root, relative to its intact transcription T (fractions of T).

    `weights` holds each activator's d_i so that combine2.fitted_module reads it as combine2's fits."""

    ok: bool
    reason: str = ""
    root: str = ""  # uncapped | capped
    q: float = math.nan  # element-driven share of intact transcription; basal fraction 1 - q
    z: float = math.nan  # V*R0/T
    w: float = math.nan  # W at the intact state (uncapped reading for a gene not fitted)
    terms: dict[str, float] = field(default_factory=dict)  # activator -> t_i = s_i*H at the intact state
    weights: dict[str, float] = field(default_factory=dict)  # activator -> d_i
    e: dict[str, float] = field(default_factory=dict)  # inhibitor -> RHO_j - 1
    s: dict[str, float] = field(default_factory=dict)  # element -> fitted strength
    vmax_frac: float = math.nan  # max_rate / T
    form: str = FORM
    half_saturation: float = math.nan  # unused: combine2.fitted_module reads it only under competition


def _drops(rho: dict[str, float]) -> tuple[dict[str, float], dict[str, float]]:
    return {el: 1.0 - r for el, r in rho.items() if r < 1.0}, {
        el: r - 1.0 for el, r in rho.items() if r > 1.0
    }


def _inhibitors(e: dict[str, float], q: float) -> tuple[dict[str, float], float] | None:
    """combine2's inhibitor strengths at element-driven share q, and R0; None when one exceeds 1."""
    s, r0 = {}, 1.0
    for el, x in e.items():
        sj = c2._strength_ok(x / (H * (q + x)))
        if sj is None:
            return None
        s[el] = sj
        r0 *= 1.0 - sj * H
    return s, r0


def fit(rho: dict[str, float]) -> CFit:
    """The uncapped root: combine2.fit('mean')'s numbers with Z = top/H, refused past the cap."""
    mean = c2.fit("mean", rho)
    if not mean.ok:
        return CFit(False, mean.reason, q=mean.q)
    d, e = dict(mean.weights), dict(mean.e)
    top, total = max(d.values()), sum(d.values())
    w = H * total / top
    if w > 1.0 + bridge.STRENGTH_TOLERANCE:
        return CFit(False, "sum_past_the_cap", q=mean.q, w=w)
    z = top / H
    r0 = 1.0
    for el in e:
        r0 *= 1.0 - mean.s[el] * H
    return CFit(
        True,
        "",
        "uncapped",
        q=mean.q,
        z=z,
        w=w,
        terms={el: x / z for el, x in d.items()},
        weights=d,
        e=e,
        s=dict(mean.s),
        vmax_frac=z / r0,
    )


def capped_root(rho: dict[str, float]) -> CFit:
    """The root with the cap binding at the intact state (FIT[FORM]), or why there is none."""
    d, e = _drops(rho)
    if not d or len(d) + len(e) != len(rho):
        return CFit(False, "no_capped_root")
    k, top, total = len(d), max(d.values()), sum(d.values())
    den = H * (k - 1) - 1.0
    if not den > 0:
        return CFit(False, "no_capped_root")
    z = ((k - 1) * top - total) / den
    if not z > total:
        return CFit(False, "no_capped_root")
    c = (1.0 - total / z) / (k - 1)
    if not z < 1.0:
        return CFit(False, "capped_root_basal_not_positive", q=z, z=z)
    inh = _inhibitors(e, z)
    if inh is None:
        return CFit(False, "capped_root_inhibitor_strength_above_1", q=z, z=z)
    s, r0 = inh
    terms = {el: x / z + c for el, x in d.items()}
    for el, t in terms.items():
        si = c2._strength_ok(t / H)
        if si is None:
            return CFit(False, "capped_root_strength_out_of_range", q=z, z=z)
        s[el] = si
    return CFit(True, "", "capped", q=z, z=z, w=1.0 + c, terms=terms, weights=d, e=e, s=s, vmax_frac=z / r0)


def level(f: CFit, open_: dict[str, bool]) -> float:
    """Transcription relative to T with the elements in `open_` open and the rest held (the closed form)."""
    w = sum(t for el, t in f.terms.items() if open_.get(el, False))
    boost = 1.0
    for el, x in f.e.items():
        if not open_.get(el, False):
            boost *= 1.0 + x / f.q
    return (1.0 - f.q) + f.z * min(1.0, w) * boost


def predicted_log2(f: CFit, gr: c2.GeneRules) -> list[float]:
    """log2 of the closed-form level relative to T in each of CELLS."""
    return [math.log2(level(f, gr.open_in(j))) for j in range(len(CELLS))]


def difference_from_mean(f: CFit, gr: c2.GeneRules) -> float:
    """The largest |log2 level under the capped sum - log2 level under the mean| over the thirteen cells."""
    # a_capped_root_alone (after 43dd4ac): no mean's prediction to compare with, which DECISION (1)
    # reads as a difference; None, so that the raise below is reached by no gene
    if not c2.fit("mean", gr.rho).ok:
        return None
    mean = c2.fit("mean", gr.rho)
    if not mean.ok:
        raise ValueError(f"{gr.gene}: the mean has no fit ({mean.reason})")
    return max(abs(a - b) for a, b in zip(predicted_log2(f, gr), c2.predicted_log2(mean, gr), strict=True))


def classify(gr: c2.GeneRules) -> dict[str, Any]:
    """Both roots of one eligible gene, and its capped-sum status (FIT['runnable'])."""
    un, cap = fit(gr.rho), capped_root(gr.rho)
    if un.ok and cap.ok:
        status = "fit_not_unique"
    elif un.ok:
        status = "runnable"
    elif cap.ok:
        status = "capped_root_only"
    else:
        status = un.reason
    return {"uncapped": un, "capped": cap, "status": status}


# ==== the run (no RNA) =====================================================================================


def run_gene(small: Module, gr: c2.GeneRules, f: CFit, rate: float, label: dict[str, str]) -> dict[str, Any]:
    """The gated run of one gene under the capped sum in every reader cell, and combine2's HARNESS checks
    (combine2.run_gene_form's, with this module's fit and closed form)."""
    m = c2.fitted_module(small, gr.gene, f, rate, label)  # type: ignore[arg-type]
    present = dict.fromkeys(gr.elements, bridge.INTACT)
    delta = math.log(2) / grn.DEFAULTS["mrna_half_life"]
    levels, closed, miss = [], [], []
    for j, cid in enumerate(CELLS):
        lv = c2._run_level(m, gr.gene, {"cell_type": cid}, present, FORM)
        want = rate * level(f, gr.open_in(j)) / delta
        levels.append(lv)
        closed.append(math.log2(want))
        if not abs(math.log2(lv / want)) <= CLOSED_FORM_TOLERANCE:
            miss.append(f"cell {cid}: run {lv:.12g}, closed form {want:.12g}")
    ungated = copy.deepcopy(m)
    for r in ungated.rules:
        r.through = ""
    ctx = {"cell_type": CELLS[0]}
    full = c2._run_level(ungated, gr.gene, ctx, present, FORM)
    for el in gr.elements:
        lv = c2._run_level(ungated, gr.gene, ctx, {**present, el: bridge.REMOVED}, FORM)
        if not abs(math.log2(lv / full) - math.log2(gr.rho[el])) <= CLOSED_FORM_TOLERANCE:
            miss.append(f"removal {el}: run {lv / full:.12g}, observed {gr.rho[el]:.12g}")
    return {"levels": levels, "s": [math.log2(x) for x in levels], "closed_s": closed, "harness_miss": miss}


# ==== the decision (counts only) ===========================================================================


def decide(summary: dict[str, Any]) -> str:
    """DECISION, from the census counts: whether a registration on 7d's design would test anything."""
    if summary["unique_fit"].get("differs_from_the_mean", {}).get("genes", 0):
        return DECISIONS[1]  # a unique fit with a prediction the mean does not give (a_capped_root_alone)
    unique = summary["unique_fit"]
    same = unique["max_abs_difference_from_mean_log2"] is None or (
        unique["max_abs_difference_from_mean_log2"] <= CLOSED_FORM_TOLERANCE
    )
    runs = [summary["roots"][r] for r in ("uncapped", "capped")]
    harness = all(r["harness_failing_genes"] == 0 for r in runs) and all(
        r["run"] == r["informative"] for r in runs
    )
    return DECISIONS[0] if same and harness else DECISIONS[1]


def size_check(summary: dict[str, Any]) -> dict[str, Any]:
    """SIZE_CHECK on the census counts: whether the genes DECISION would register on reach MIN_GENES."""
    n = summary["unique_fit"].get("differs_from_the_mean", {}).get("informative_symbol_on_one_gene_id", 0)
    return {
        "genes_whose_unique_fit_differs_from_the_mean": n,
        "min_genes": MIN_GENES,
        "under_7d_decision_2": "NOT_TESTABLE" if n < MIN_GENES else "at_or_above_the_registered_size",
        "registration_written": n >= MIN_GENES and decide(summary) == DECISIONS[1],
        "rna_opened": False,
    }
