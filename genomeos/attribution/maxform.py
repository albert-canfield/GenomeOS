# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""lane-maxform (ROADMAP area I, 2026-10-10): the `max` combination form, and `or`, against measured RNA.

STOP LIST (docs/ROADMAP.md section 5, "The stop list of 2026-09-27, reinstated"): no AlphaGenome request and
no deletion request, no new biosample (the thirteen reader cells only; READER_TERMS unchanged), nothing in
area J, no node-model work, no BioForge calibration. Nothing here fetches, sends or spends.

WHY THIS LANE EXISTS. Three lanes established today, independently, that this RNA endpoint cannot test an
ADDITIVE rule on 7d's population:

* lane-capsumrna (43dd4ac, 8994fca, f74dca5; data/results/capsumrna_census.json): the capped sum's
  prediction equals the runtime mean's on every one of the 2,288 genes it fits (largest difference
  6.66e-16 log2); only 3 genes have a fit whose prediction differs, under 7d's MIN_GENES of 100, so
  NOT_TESTABLE;
* lane-envelope (b0eefdc, df02887; data/results/envelope_census.json): `additive_envelope` is an INTERVAL,
  not a point form -- its lower end IS grn's `max`, its upper end with its cap IS grn's `sum_capped` (four
  identity checks returning 0.0 exactly over 20,000 draws). As its own form it is lane-capsumrna's form;
* lane-rootpick (55a2cff, eea998b, c0f0285): under the registered root convention the capped sum is not a
  distinct form on this population at all.

THE PREMISE, AND WHAT THIS LANE FOUND OF IT (`PREMISE`). The structural reason offered for those three
results is that the runtime mean is LINEAR in the element terms under the reader's hold, and so is any sum
below its cap, so no additive rule can separate from the mean at 7d's size; and therefore that only a rule
NON-LINEAR under the hold -- `max`, `or` -- or a capped root could be tested on this endpoint. This module
checks that premise from the algebra rather than taking it on trust, and the check splits it in two
(`PREMISE`): the linearity half HOLDS, and the conclusion drawn from it is right about `or` and WRONG about
`max`, for a reason that has nothing to do with linearity. `max` IS non-linear under the hold, and that is
exactly why it has NO FIT on any gene with two or more activating elements: under `max` removing a
non-maximal activator changes nothing, so it cannot reproduce two measured drops at once. grn's own
registered case list already says so (`grn.COMBINATION_TEST`, "crispri_training_two_decreases": "impossible
under max (removing a non-maximal activator changes nothing)"). Non-linearity is necessary for a form to
separate from the mean here; it is not sufficient, because a form must first FIT the bridge's removal
observations, and `max` cannot.

`or` is the form this premise was reaching for: A = 1 - prod(1 - t_i) is non-linear under the hold, every
activator's removal lowers it, its fit at combine2's conventions exists, is unique in closed form and gives
a prediction that is not the mean's. So this lane's registered arm is `or`, the one form of
`grn.ACTIVATOR_RULES` never held against measured RNA, and `max` is reported NOT_TESTABLE from the census
with the structural reason, not from a size.

Nothing in `runtime/grn.py` is changed: both forms are already there (`grn.ACTIVATOR_RULES`, since
2026-09-28), selectable for one run with `NetworkRuntime(activator_combination=...)`.

No RNA value is read in the census. No model request, no key, no money.
"""

from __future__ import annotations

import copy
import math
from dataclasses import dataclass, field
from typing import Any

from genomeos.attribution import bridge
from genomeos.attribution import combine2 as c2
from genomeos.attribution import magcal as mc
from genomeos.attribution import multrna as mr
from genomeos.ir import Module
from genomeos.runtime import grn

DATE = "2026-10-10"
LANE = "lane-maxform"
ROADMAP_ITEM = (
    "area I: the `max` combination form (and `or`, the other form of grn.ACTIVATOR_RULES never held"
    " against measured RNA) on 7d's population, registered"
)
#: the two forms this lane examines, as runtime/grn.py names them (both in ACTIVATOR_RULES since
#: 2026-09-28; the engine is unchanged). `max` is censused and reported; `or` is the arm registered.
FORMS = ("max", "or")
FORM = "or"  # the arm the registration judges, if the census reaches 7d's MIN_GENES
SELECTION = c2.SELECTION
CELLS = c2.CELLS
H = bridge.H_INTACT
INFORMATIVE_LOG2 = c2.INFORMATIVE_LOG2
CLOSED_FORM_TOLERANCE = c2.CLOSED_FORM_TOLERANCE
RUNTIME = c2.RUNTIME
MIN_GENES = c2.MIN_GENES

FORM_WORDS = {
    "max": "the strongest-element form (grn `max`)",
    "or": "the independent-binding form (grn `or`)",
}
FORM_SOURCES = {
    "max": (
        "runtime/grn.py combine_activators 'max': A = max over the gene's activating rules of s_i*H(x_i),"
        " a held rule counted at strength 0 (so A = 0 with every activator held). One of the four rules"
        " lane-combine compared on 2026-09-28 (grn.ACTIVATOR_RULES), selectable for one run with"
        " NetworkRuntime(activator_combination='max'). Inhibitors multiply, R = prod (1 - s_j*H(x_j)), as"
        " under the mean. It is the lower end of lane-strengthform's additive_envelope (lane-envelope's"
        " identity checks, df02887)"
    ),
    "or": (
        "runtime/grn.py combine_activators 'or': A = 1 - prod over the gene's activating rules of"
        " (1 - s_i*H(x_i)), independent binding of any one activator sufficing; a held rule is counted at"
        " strength 0 and so drops out of the product. One of the four rules lane-combine compared on"
        " 2026-09-28 (grn.ACTIVATOR_RULES; grn.COMBINATION_TEST, 'or max <= P <= S1+S2'), selectable for"
        " one run with NetworkRuntime(activator_combination='or'). Inhibitors multiply, as under the mean"
    ),
}
CONVENTIONS = (
    "combine2's, unchanged (FIT['principle'] and FIT['mean']): every element at its intact state"
    " x = bridge.INTACT, so H(x) = bridge.H_INTACT = 0.5; the strongest activating element at"
    " bridge.DECLARED_STRENGTH (s = 1); an inhibitor's s_j*h = e_j/(q + e_j) at the form's own"
    " element-driven share q; the bridge's Hill assumptions (K = 1, n = 2, h = 0.5). No number is added"
    " and no parameter is tuned to any measured RNA value"
)

# ---- the algebra, written before any count -------------------------------------------------------------
PREMISE = {
    "as_handed_to_this_lane": (
        "lane-envelope's conclusion, carried in this lane's brief: on 7d's design the runtime mean is linear"
        " under the reader's hold, and so is any sum below its cap, so no additive rule can separate from"
        " the mean at 7d's size; only a rule non-linear under the hold -- `max`, `or` -- or the capped root"
        " could be tested on this endpoint"
    ),
    "the_linear_half_holds": (
        "VERIFIED here from the same algebra, not taken on trust. Module.active_rules returns a held rule as"
        " a strength-0 copy and keeps it, and the bridge's removal is a clamp of the element at 0 with the"
        " rule kept (combine2.STRUCTURE['mean_is_linear_under_the_hold']). Under both, the mean's drive"
        " A = (1/n) sum t_i is linear in the element terms, so closing element i lowers transcription by"
        " d_i of T whatever the other elements do, and a capped sum below its cap is the same function of"
        " the open set up to the scale the fit absorbs. Hence an additive rule's centred prediction IS the"
        " mean's on every gene both fit, which is what lane-capsumrna measured (largest difference"
        " 6.66e-16 log2 over 2,288 genes) rather than argued"
    ),
    "the_conclusion_is_wrong_about_max": (
        "non-linearity under the hold is NECESSARY for a form to separate from the mean here and NOT"
        " SUFFICIENT, because a form must first reproduce every element's own removal observation"
        " (combine2.FIT['principle']). Under `max`, A = max t_i, so removing a non-maximal activator leaves"
        " A unchanged: d_i = Z*(A_full - A_without_i) = 0 for every activator but the strongest. 7d's"
        " POPULATION (b) requires every rule's observation to be readable and not exactly 1, so every"
        " activating element has d_i > 0, and `max` therefore has NO FIT on any gene with two or more"
        " activating elements -- not a worse fit, no fit. grn.COMBINATION_TEST already registered the same"
        " impossibility from a different measurement ('crispri_training_two_decreases': 'impossible under"
        " max (removing a non-maximal activator changes nothing)'). On a gene with exactly ONE activating"
        " element (and one or more inhibitors, which is how such a gene is in 7d's multi-element"
        " population) `max` reduces to that single term and is the mean element for element: the same q,"
        " the same fitted strengths, the same V, the same prediction in every reader cell. grn says as much"
        " in ACTIVATOR_RULES' own comment, 'all four agree when a gene has one activator'. So on 7d's"
        " population `max` either has no fit or gives the runtime mean's prediction exactly, and the count"
        " of genes where its prediction differs from the mean's is 0 BEFORE any RNA value is read"
    ),
    "or_is_the_form_the_premise_was_reaching_for": (
        "`or` is non-linear under the hold and every activator's removal lowers it, so it can reproduce k"
        " drops at once: A = 1 - P with P = prod(1 - t_i), and removing i gives A = 1 - P/(1 - t_i), so"
        " d_i = Z*P*t_i/(1 - t_i) > 0 for every activator. FIT['or'] solves it in closed form, the solution"
        " is unique at combine2's conventions, and its centred prediction is not the mean's wherever the"
        " gene has two or more activating elements. It is the arm this registration judges"
    ),
    "what_this_lane_does_not_claim": (
        "that no form could test this endpoint. It establishes that an additive form cannot (three lanes"
        " today) and that `max` cannot be fitted at all on this population, and it tests `or`. A form"
        " outside grn.ACTIVATOR_RULES, a different population, or a fresh measurement are each untouched by"
        " any of it"
    ),
}
FIT = {
    "principle": c2.FIT["principle"],
    "max": (
        "d_i = 1 - RHO_i (activating), e_j = RHO_j - 1 (inhibiting). With two or more activating elements"
        " there is no fit (max_needs_one_activating_element): A = max t_i, so removing any activator but the"
        " strongest leaves A, and hence the level, unchanged, while every activating element of 7d's"
        " population has an observed drop. With exactly one activating element, A = t_1 with it open and 0"
        " with it held, so d_1 = Z*H at s_1 = 1: Z = d_1/H, q = Z*H = d_1, b = T(1 - q),"
        " s_j*h = e_j/(q + e_j), V = T*Z/R0 -- combine2.fit('mean')'s numbers for n = 1, element for"
        " element. Level relative to T with the open set O: (1 - q) + Z*max(t_i over open activators, else"
        " 0) * prod over closed inhibitors (1 + e_j/q), which for one activator is the mean's"
        " (1 - q) + [d_1 if open] * that product"
    ),
    "or": (
        "d_i = 1 - RHO_i (activating), e_j = RHO_j - 1 (inhibiting), top = max d_i. Writing u_i ="
        " t_i/(1 - t_i) and P = prod(1 - t_i), the removal of activator i gives d_i = Z*P*u_i, so u_i is"
        " proportional to d_i; the DECLARED_STRENGTH convention puts the strongest activator at s = 1, i.e."
        " t_top = H = 0.5 and u_top = H/(1 - H) = 1, which fixes the constant: u_i = d_i/top,"
        " t_i = d_i/(top + d_i), P = prod top/(top + d_i), Z = top/P and the element-driven share"
        " q = Z*(1 - P) = Z - top. Every activator's s_i = t_i/H = 2*d_i/(top + d_i) is in (0, 1] by"
        " construction, so only q < 1 and the inhibitors can refuse a gene. b = T(1 - q),"
        " s_j*h = e_j/(q + e_j), V = T*Z/R0 with R0 = prod(1 - s_j*h). Level relative to T with the open"
        " set O: (1 - q) + Z*(1 - prod over open activators of (1 - t_i)) * prod over closed inhibitors"
        " (1 + e_j/q)"
    ),
    "the_or_root_is_unique": (
        "u_i proportional to d_i is forced by the k removal equations, and the one free constant is fixed by"
        " the DECLARED_STRENGTH convention, whose 'strongest' is unambiguous because t_i = d_i/(top + d_i)"
        " is increasing in d_i, so argmax s = argmax d. There is no second root, and lane-rootpick's"
        " registered root convention (both roots exist -> the uncapped; the capped root alone -> the gene is"
        " UNKNOWN and leaves the population) is therefore not triggered by either form here: it is cited"
        " and applied vacuously. Under `max` the question does not arise either -- a gene with two or more"
        " activating elements has no root at all"
    ),
    "where_the_forms_can_differ": (
        "MEASURED, not assumed, and it is what the census counts: every form here is fitted to the same"
        " single-removal observations, so on any state in which AT MOST ONE of a gene's activating elements"
        " is closed, every form gives the same level by construction. A form can only separate from the mean"
        " in a cell where TWO OR MORE of the gene's activating elements are closed at once. So a gene whose"
        " thirteen reader cells never close two of its activating elements together has `or`'s prediction"
        " equal to the mean's in all thirteen, however non-linear the form is, and the census counts it as"
        " not differing. The census question is therefore a question about the reader's open/closed"
        " pattern as much as about the form, and it is answered on the thirteen actual reader cells and not"
        " on hypothetical states"
    ),
    "q_rises_with_the_product": (
        "q_or = top*(prod(1 + d_i/top) - 1) >= sum d_i = q_mean, with equality only at one activating"
        " element: `or` saturates, so a single removal moves a nearly-saturated drive very little and the"
        " dynamic range Z needed to reproduce the measured drop is larger. Genes with several strong"
        " activating elements therefore leave no positive basal rate under `or` (basal_not_positive) and"
        " are not run, exactly as 7d's own populations exclude genes whose summed drops exceed T. The"
        " opposite also follows: q_or above q_mean makes an inhibitor's fitted s_j = 2*e_j/(q + e_j)"
        " SMALLER, so `or` admits inhibitors the mean refuses, and its population is not a subset of 7d's"
    ),
    "runnable": (
        "a gene is runnable under a form when the fit gives b > 0 (q < 1) and every fitted s in (0, 1]"
        " (bridge.STRENGTH_TOLERANCE above 1 allowed), 7d's FIT['runnable'] unchanged; otherwise it is"
        " reported by reason and not run"
    ),
}
STRUCTURE = {
    "what_is_assumed": c2.STRUCTURE["what_is_assumed"],
    "inhibitor_only": c2.STRUCTURE["inhibitor_only"],
    "network": c2.STRUCTURE["network"],
}
CENSUS = (
    "a count on the predictor side, before any RNA value is read (data/results/maxform_census.json), over"
    " the 24 compiled programs gated as 7d gated them: of the multi-element genes combine2.classify admits,"
    " how many each of `max` and `or` fits and why the rest are not fitted; on every gene a form fits, the"
    " largest difference over the thirteen reader cells between its closed-form prediction and the runtime"
    " mean's, and the count of genes where that difference exceeds CLOSED_FORM_TOLERANCE -- the census"
    " question this lane was given; 7d's mean and competition counts and 7h's multiplicative counts"
    " recomputed beside them and checked against data/results/combine2_census.json and"
    " multrna_census.json; every informative gene RUN by the runtime under the form with combine2's"
    " HARNESS checks; and the genes each form shares with 7d's and 7h's populations. It reads no RNA value"
)
CENSUS_SENTENCE = (
    "of the {n} multi-element genes combine2.classify admits, `max` fits {max_fitted} ({max_differ} with a"
    " prediction other than the runtime mean's; not fitted: {max_reasons}) and `or` fits {or_fitted}"
    " ({or_differ} with a prediction other than the runtime mean's, largest difference {or_dmax} in log2;"
    " not fitted: {or_reasons}). Informative and with a symbol on exactly one gene_id: {max_pop} under"
    " `max`, {or_pop} under `or`, against 7d's MIN_GENES of {min}. The runtime under"
    " activator_combination gives each form's closed form on every informative gene it ran ({ran} runs,"
    " {miss} harness failures), before the RNA tables' row check"
)
DECISION_BEFORE_THE_CENSUS = (
    "fixed before the census was counted, from PREMISE: per form, the genes whose fit gives a prediction"
    " differing from the runtime mean's by more than CLOSED_FORM_TOLERANCE in any reader cell, informative"
    " and with a symbol on exactly one gene_id, are the population a registration could test. Fewer than"
    " MIN_GENES of them is NOT_TESTABLE for that form (7d's DECISION 2, carried unchanged): no registration"
    " is written for it, no RNA table is opened for it, and the census is the result. At or above MIN_GENES"
    " the form is registered and run once. If BOTH forms are NOT_TESTABLE, the finding is that this"
    " endpoint separates no form tried so far -- additive (three lanes today), strongest-element, or"
    " independent-binding -- from the runtime mean on 7d's population, which is itself the result"
)
DECISIONS = ("NOT_TESTABLE", "register_and_run")
EXPLORATION_PRECEDED_THIS = (
    "before this code was committed this lane read, predictor-side only: runtime/grn.py's"
    " combine_activators, ACTIVATOR_RULES and COMBINATION_TEST; combine2.py, capsumrna.py, magcal.py and"
    " multrna.py with their census, registration and run scripts; and the lane notes of lane-capsumrna,"
    " lane-envelope, lane-rootpick and lane-magcal from `genomeos work list`. It also read lane-magcal's"
    " committed M2h reading (ac5580b), which is a measured result on this RNA and is named in"
    " PRIOR_EXPOSURE for that reason. No RNA quantification table and no genes table of 7d or 7h was"
    " opened before the census was saved; PREMISE and FIT were written before any count was run"
)


# ==== the fits (closed form; no RNA) =====================================================================


@dataclass
class MFit:
    """One gene's fit under one of this lane's forms, relative to its intact transcription T.

    `weights` holds each activator's d_i and `terms` its t_i, so combine2.fitted_module reads it as it
    reads combine2's own fits (it needs `weights` only to tell an activator from an inhibitor)."""

    form: str
    ok: bool
    reason: str = ""
    q: float = math.nan  # element-driven share of intact transcription; basal fraction 1 - q
    z: float = math.nan  # V*R0/T
    terms: dict[str, float] = field(default_factory=dict)  # activator -> t_i = s_i*H at the intact state
    weights: dict[str, float] = field(default_factory=dict)  # activator -> d_i
    e: dict[str, float] = field(default_factory=dict)  # inhibitor -> RHO_j - 1
    s: dict[str, float] = field(default_factory=dict)  # element -> fitted strength
    vmax_frac: float = math.nan  # max_rate / T
    half_saturation: float = math.nan  # unused: combine2.fitted_module reads it only under competition


def _split(rho: dict[str, float]) -> tuple[dict[str, float], dict[str, float]] | None:
    """(activating d_i, inhibiting e_j), or None when an observation is exactly 1."""
    d = {el: 1.0 - r for el, r in rho.items() if r < 1.0}
    e = {el: r - 1.0 for el, r in rho.items() if r > 1.0}
    return None if len(d) + len(e) != len(rho) else (d, e)


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


def fit(form: str, rho: dict[str, float]) -> MFit:
    """The joint fit of one gene's elements under `form` (FIT), from their removal folds RHO."""
    if form not in FORMS:
        raise ValueError(f"form {form!r} is not one of {FORMS}")
    split = _split(rho)
    if split is None:
        return MFit(form, False, "observation_null_effect")
    d, e = split
    if not d:
        return MFit(form, False, "no_activating_rule")
    top = max(d.values())
    if form == "max":
        if len(d) > 1:
            return MFit(form, False, "max_needs_one_activating_element")
        z, q, terms = top / H, top, {el: H for el in d}
    else:
        p = 1.0
        for x in d.values():
            p *= top / (top + x)
        z = top / p
        q = z - top
        terms = {el: x / (top + x) for el, x in d.items()}
    if not q < 1.0:
        return MFit(form, False, "basal_not_positive", q=q, z=z)
    inh = _inhibitors(e, q)
    if inh is None:
        return MFit(form, False, "inhibitor_strength_above_1", q=q, z=z)
    s, r0 = inh
    for el, t in terms.items():
        si = c2._strength_ok(t / H)
        if si is None:
            return MFit(form, False, "activator_strength_out_of_range", q=q, z=z)
        s[el] = si
    return MFit(form, True, "", q, z, terms, dict(d), e, s, z / r0)


def drive(f: MFit, open_: dict[str, bool]) -> float:
    """The activator drive A under the form, with the elements in `open_` open and the rest held."""
    ts = [t for el, t in f.terms.items() if open_.get(el, False)]
    if f.form == "max":
        return max(ts, default=0.0)
    miss = 1.0
    for t in ts:
        miss *= 1.0 - t
    return 1.0 - miss


def level(f: MFit, open_: dict[str, bool]) -> float:
    """Transcription relative to T with the elements in `open_` open and the rest held (the closed form)."""
    boost = 1.0
    for el, x in f.e.items():
        if not open_.get(el, False):
            boost *= 1.0 + x / f.q
    return (1.0 - f.q) + f.z * drive(f, open_) * boost


def predicted_log2(f: MFit, gr: c2.GeneRules) -> list[float]:
    """log2 of the closed-form level relative to T in each of CELLS."""
    return [math.log2(level(f, gr.open_in(j))) for j in range(len(CELLS))]


def difference_from_mean(f: MFit, gr: c2.GeneRules) -> float | None:
    """The largest |log2 level under the form - log2 level under the runtime mean| over the thirteen cells,
    or None when the mean has no fit for the gene (`or` admits inhibitors the mean refuses, FIT)."""
    mean = c2.fit("mean", gr.rho)
    if not mean.ok:
        return None
    return max(abs(a - b) for a, b in zip(predicted_log2(f, gr), c2.predicted_log2(mean, gr), strict=True))


def differs_from_mean(f: MFit, gr: c2.GeneRules) -> bool:
    """The census question: does this fit's prediction differ from the runtime mean's at all?

    A gene the mean does not fit counts as differing: the mean has no prediction to coincide with."""
    diff = difference_from_mean(f, gr)
    return True if diff is None else diff > CLOSED_FORM_TOLERANCE


# ==== the run (no RNA) ===================================================================================


def run_gene(small: Module, gr: c2.GeneRules, f: MFit, rate: float, label: dict[str, str]) -> dict[str, Any]:
    """The gated run of one gene under the form in every reader cell, and combine2's HARNESS checks
    (combine2.run_gene_form's, with this module's fit and closed form)."""
    m = c2.fitted_module(small, gr.gene, f, rate, label)  # type: ignore[arg-type]
    present = dict.fromkeys(gr.elements, bridge.INTACT)
    delta = math.log(2) / grn.DEFAULTS["mrna_half_life"]
    levels, closed, miss = [], [], []
    for j, cid in enumerate(CELLS):
        lv = c2._run_level(m, gr.gene, {"cell_type": cid}, present, f.form)
        want = rate * level(f, gr.open_in(j)) / delta
        levels.append(lv)
        closed.append(math.log2(want))
        if not abs(math.log2(lv / want)) <= CLOSED_FORM_TOLERANCE:
            miss.append(f"cell {cid}: run {lv:.12g}, closed form {want:.12g}")
    ungated = copy.deepcopy(m)
    for r in ungated.rules:
        r.through = ""
    ctx = {"cell_type": CELLS[0]}
    full = c2._run_level(ungated, gr.gene, ctx, present, f.form)
    for el in gr.elements:
        lv = c2._run_level(ungated, gr.gene, ctx, {**present, el: bridge.REMOVED}, f.form)
        if not abs(math.log2(lv / full) - math.log2(gr.rho[el])) <= CLOSED_FORM_TOLERANCE:
            miss.append(f"removal {el}: run {lv / full:.12g}, observed {gr.rho[el]:.12g}")
    return {"levels": levels, "s": [math.log2(x) for x in levels], "closed_s": closed, "harness_miss": miss}


# ==== the census decision (counts only) ==================================================================


def decide_from_census(per_form: dict[str, dict[str, Any]], form: str) -> str:
    """DECISION_BEFORE_THE_CENSUS, from the census counts for one form."""
    n = per_form[form].get("differs_from_the_mean", {}).get("informative_symbol_on_one_gene_id", 0)
    return DECISIONS[1] if n >= MIN_GENES else DECISIONS[0]


# ==== the identity checks (predictor side, no RNA) =======================================================


def max_identity_checks(draws: int = 20000, seed: int = 20261010) -> dict[str, float]:
    """PREMISE['the_conclusion_is_wrong_about_max'], measured rather than argued, on random observations.

    Returns the largest violation of each claim over `draws` random RHO vectors; each should be 0.0:
    (1) `max` has no fit whenever two or more observations are activating; (2) with exactly one activating
    observation its prediction equals the runtime mean's in every reader-cell open set; (3) `or`'s fit,
    where it exists, reproduces every single removal observation exactly."""
    import random

    rng = random.Random(seed)
    worst = {"max_fits_two_activators": 0.0, "max_differs_from_mean_at_one": 0.0, "or_removal_error": 0.0}
    for _ in range(draws):
        k = rng.randint(2, 5)
        rho = {f"E{i}": rng.choice([rng.uniform(0.3, 0.999), rng.uniform(1.001, 2.0)]) for i in range(k)}
        n_act = sum(1 for r in rho.values() if r < 1.0)
        fm, fo = fit("max", rho), fit("or", rho)
        if n_act >= 2 and fm.ok:
            worst["max_fits_two_activators"] = 1.0
        if n_act == 1 and fm.ok:
            mean = c2.fit("mean", rho)
            if mean.ok:
                for _ in range(4):
                    op = {el: rng.random() < 0.5 for el in rho}
                    err = abs(math.log2(level(fm, op)) - math.log2(c2.level(mean, op)))
                    worst["max_differs_from_mean_at_one"] = max(worst["max_differs_from_mean_at_one"], err)
        if fo.ok:
            allopen = dict.fromkeys(rho, True)
            base = level(fo, allopen)
            for el in rho:
                op = dict(allopen)
                op[el] = False
                worst["or_removal_error"] = max(
                    worst["or_removal_error"], abs(level(fo, op) / base - rho[el])
                )
    return worst


# ==== the registered test (written before any RNA value of this lane's population was read) ==============

QUESTION = (
    "For genes the compiled programs drive through two or more elements, does a reader-gated run under the"
    " independent-binding combination form (grn `or`: A = 1 - prod(1 - t_i)), at one magnitude scale fitted"
    " on training genes and applied unchanged to held-out genes, track measured RNA across the thirteen"
    " reader cell types better than a no-rule baseline on the held-out genes? And is `max` testable on this"
    " endpoint at all?"
)
#: 7d's POPULATION, word for word, with the form bound to this lane's arm
POPULATION = c2.POPULATION
POPULATION_NOTE = (
    "7d's POPULATION above is carried verbatim and its clause (d), 'the joint fit under F is runnable', is"
    " read with F = `or` and FIT['or']; clause (e)'s informativeness and clause (f)'s symbol and row checks"
    " are unchanged. The population is therefore NOT 7d's gene set: `or`'s q exceeds the mean's, so genes"
    " whose summed drops saturate it leave no positive basal rate and drop out, while `or` admits"
    " inhibitors the mean refuses and keeps genes 7d's mean does not run (FIT['q_rises_with_the_product'])."
    " The census counts the overlap, and the head-to-head is on the genes the compared forms both run"
)
PRIOR_EXPOSURE = (
    "STATED PLAINLY, AND NOT A CAVEAT ADDED AFTERWARDS. This measurement is not fresh for this project. The"
    " same ENCODE polyA plus RNA-seq of the same thirteen reader biosamples, under the same log2(TPM+1)"
    " endpoint and the same no-rule baseline, has already been read by: item 7 (M2b, lane-readercompile,"
    " data/results/readergate_result.json) including its post-hoc information term; item 7d (M2c,"
    " lane-combine2) for the runtime's mean and promoter competition; item 7h (M2f, lane-multrna) for the"
    " multiplicative composition; item M2h (lane-magcal, ac5580b) for all three of those forms at one"
    " fitted magnitude scale on a held-out split; and today's capped-sum census (lane-capsumrna) and"
    " additive-envelope census (lane-envelope), both of which reached NOT_TESTABLE without opening an RNA"
    " table. The genes this lane reads are largely the same genes, and the choice to try another"
    " combination form at all was made BECAUSE the earlier readings came out worse than no rule"
    " uncalibrated and because M2h's calibrated multiplicative arm came out better. The held-out split"
    " limits that exposure and does not eliminate it: the split is fixed by this registration before any"
    " held-out number is read and no held-out gene enters any fit, but those genes' RNA values were in the"
    " earlier readings' populations, the endpoint and baseline were chosen with those readings in hand, and"
    " nothing here is a blind prediction of an unseen measurement. Any reading from this lane therefore"
    " holds as a statement about THIS measurement and is not independent confirmation of it; a magnitude"
    " read here is a hypothesis for a fresh measurement, not a result on one"
)
PRIOR_EXPOSURE_CLAUSE = mc.PRIOR_EXPOSURE_CLAUSE
#: lane-magcal's split, taken unchanged so that a gene falls the same side as it did in M2h
SPLIT = mc.SPLIT
PREDICTIONS = {
    "form": (
        "s_F(G, X) = log2(L_F(G, X)), L_F(G, X) the level the runtime gives G in reader cell X under"
        " F = `or` (NetworkRuntime(..., activator_combination='or'), the fitted rules acting through their"
        " elements, the reader blocks in force, every element clamped present so that only the reader holds"
        " a rule). P_F(G, X) = s_F(G, X) - mean over the thirteen cells of s_F(G, X'), 7d's centring"
        " unchanged"
    ),
    "scaled_form": (
        "a*P_F(G, X), one scale a per form for the whole population, fitted on the TRAINING genes' pairs"
        " alone and applied unchanged to the held-out genes (magcal.FIT, taken unchanged)"
    ),
    "no_rule_baseline": c2.PREDICTIONS["no_rule_baseline"],
}
FIT_SCALE = mc.FIT
METRIC = mc.METRIC
#: 95%/4 -> 98.75%. RULED by the supervisor through the coordinator on 2026-10-10 and reasoned here: this
#: is the same question, the same RNA and the same split as M2h's three arms, so this lane's arm joins them
#: as ONE FAMILY OF FOUR -- {7d's runtime mean, 7d's promoter competition, 7h's multiplicative composition,
#: this lane's arm} -- and the correction is Bonferroni over four. M2h's three readings stay as registered
#: at 98.33% and are NOT recomputed; nothing of lane-magcal's is touched. Fixed here before any held-out
#: value exists and not changed afterwards.
FORM_INTERVAL = (0.625, 99.375)
FORM_INTERVAL_WORDS = "98.75% CI"
FAMILY = (
    "the family corrected over is four: 7d's runtime mean, 7d's promoter-competition form, 7h's"
    " multiplicative composition and this lane's arm. The same RNA, the same endpoint, the same no-rule"
    " baseline and the same held-out split, so a fourth arm changes the correction; 95%/4 gives 98.75%"
    " two-sided, percentiles 0.625 and 99.375. M2h's three arms keep their registered 98.33% intervals and"
    " are not recomputed here"
)
BOOTSTRAP = {
    "resamples": mc.BOOTSTRAP["resamples"],
    "seed": mc.BOOTSTRAP["seed"],
    "unit": mc.BOOTSTRAP["unit"],
    "interval": (
        "the 98.75% percentile interval (Bonferroni over the family of four, FAMILY) of the held-out Delta"
        " MSE at the fitted scale, the scale held fixed; and the 98.75% percentile interval of the scale"
        " itself over TRAINING gene resamples, refitted in each. magcal's `_gene_terms` and `_sums` build"
        " the resampling statistic unchanged; only the percentile pair is this registration's, because"
        " magcal's bootstrap_delta and bootstrap_scale read magcal's own FORM_INTERVAL and lane-magcal's"
        " files are not edited by this lane. tests/test_maxform.py checks that this lane's bootstrap, run at"
        " magcal.FORM_INTERVAL, reproduces magcal.bootstrap_delta and magcal.bootstrap_scale exactly"
    ),
}
CONTROL = mc.CONTROL
CONTROL_NOTE = (
    "magcal's control, taken unchanged: a WITHIN-GENE cell permutation, which asks whether the rules'"
    " cell-to-cell pattern carries anything beyond a gene-level offset, the question the supervisor's audit"
    " of M2h called the right one. control p = (1 + permutations with held-out Delta MSE >= observed) /"
    " (1 + PERMUTATIONS['count']), so with 1,000 permutations the smallest p it can report is 0.000999: a p"
    " of 0.000999 is the FLOOR of the control and not a smaller number"
)
PERMUTATIONS = mc.PERMUTATIONS
CONTROL_ALPHA = mc.CONTROL_ALPHA
HARNESS = (
    "for every gene run under the form: (1) in each of the thirteen reader cells the gated run's level"
    " equals T times the form's closed form for that cell's open set; (2) the same fitted program with its"
    " rules ungated, every element present, then each element clamped absent in turn, reproduces that"
    " element's RHO. Each within CLOSED_FORM_TOLERANCE in log2 (7d's HARNESS, unchanged). Any miss stops"
    " the run before the RNA is opened (HARNESS_DEFECT). Beside it, the census's own checks:"
    " max_identity_checks must return 0.0 on all three claims, and the recomputed 7d and 7h census counts"
    " must equal the committed ones"
)
TOLERANCE = mc.TOLERANCE
WHICH_GOVERNS = mc.WHICH_GOVERNS
NOT_VACUOUS = mc.NOT_VACUOUS
SCALE_CLAUSES = mc.SCALE_CLAUSES
#: magcal's four-way table, applied as written (magcal.decide); its OUTCOMES sentences are reused word for
#: word, with this lane's form words and interval words filled in
OUTCOMES = mc.OUTCOMES
DECISION = (
    "1. The census is counted first and reads no RNA. Per form, fewer than MIN_GENES genes whose prediction"
    " differs from the runtime mean's (informative, symbol on exactly one gene_id) is NOT_TESTABLE and the"
    " census is that form's result: no registration, no RNA table opened"
    " (DECISION_BEFORE_THE_CENSUS). 2. If any gene misses the HARNESS, the run stops before the RNA is"
    " opened (HARNESS_DEFECT). 3. The split is applied (SPLIT). 4. One scale is fitted on the TRAINING"
    " pairs alone and bootstrapped on training genes. 5. ONLY THEN is the held-out Delta MSE computed at"
    " that scale, with its 98.75% interval and the within-gene cell permutation control, and magcal's"
    " four-way table (WHICH_GOVERNS, magcal.decide) applied as written: held-out interval wholly above 0 ->"
    " BETTER, wholly below -> WORSE, straddling 0 with the scale interval also straddling 0 ->"
    " CALIBRATION_REMOVES, straddling 0 with the scale interval excluding 0 -> INCONCLUSIVE. 6. The"
    " uncalibrated head-to-heads against the runtime mean and the multiplicative form, over the genes both"
    " forms run, are reported beside it as DESCRIPTIVE and decide nothing, because this lane's registered"
    " endpoint is the calibrated held-out one. 7. The reading is magcal's sentence for the outcome with the"
    " prior-exposure clause attached, always"
)
INCONCLUSIVE = (
    "'inconclusive' is magcal's INCONCLUSIVE and means exactly this: the held-out Delta MSE interval"
    " includes 0 while the fitted scale's interval excludes it, so the form carries a magnitude this"
    " measurement can see but the held-out gain at that magnitude is not separated from zero at 98.75%."
    " It is not evidence that the form tracks the RNA and not evidence that it does not. It is distinct"
    " from CALIBRATION_REMOVES, where the scale interval includes 0 as well and the best scale this"
    " measurement supports may be no rule at all, and from NOT_TESTABLE, which is a size verdict reached"
    " before any RNA value is read"
)
FALSIFIER = (
    "the form's held-out reading does not beat the no-rule baseline at the registered interval (its 98.75%"
    " held-out Delta MSE interval includes 0 or lies below it): run under this form at one fitted"
    " magnitude, the compiled rules of multi-element genes do not make a run track measured RNA across the"
    " reader's cell types in this population better than no rule. The narrower claim that the rules'"
    " CELL-TO-CELL pattern carries it, rather than a gene-level offset, is falsified when the within-gene"
    " permutation control is not beaten (control p > CONTROL_ALPHA). The claim that the form carries a"
    " magnitude at all is falsified when the fitted scale's interval includes 0. A negative fitted scale"
    " whose interval lies wholly below 0 falsifies the form's DIRECTION on this population"
    " (magcal.scale_is_negative), and is reported as such and never as a gain"
)
HONEST_SIZE = (
    "required of any 'better' reading, as the supervisor requires of M2h's: the gain is stated as a"
    " fraction of mean(D^2) over the held-out pairs beside the Delta MSE, so that a positive Delta MSE is"
    " read against the variance it is explaining and not alone"
)
STOP_RULE = (
    "one run, from this registration's committed code, on the files the registration pins by sha256 (the 24"
    " compiled programs, the rate table, the gene table, the narrowPeak files, this lane's census and the"
    " thirteen RNA quantifications). No form, parameter, population rule, split, scale definition,"
    " endpoint, pseudocount, metric, interval, control, decision, falsifier or wording is changed after a"
    " held-out value is read; a defect found afterwards is reported beside the registered reading and"
    " corrected only under a new registration"
)
LIMITATIONS = (
    c2.LIMITATIONS[0],
    "the fit treats a gene's k observations as taken at one intact state, but each is the argmax over 371"
    " tracks in its own track, and each reader cell's run transports them there. A maximum per element"
    " inflates each element's fold; under `or` that inflation raises q and so removes genes from the"
    " population rather than widening their predicted range",
    POPULATION_NOTE,
    "`or` is a form the engine carries and lane-combine compared against published cases on 2026-09-28; it"
    " is not a published model of these genes. Independent binding of any one of up to dozens of compiled"
    " elements, in thirteen cell types, is a declared form, not a measured mechanism",
    *c2.LIMITATIONS[3:],
)
DESCRIPTIVES = (
    "the uncalibrated held-out and whole-population Delta MSE at scale 1, beside the calibrated one",
    "the head-to-heads against 7d's runtime mean and 7h's multiplicative composition over the genes both"
    " forms run, uncalibrated and at each form's own scale: descriptive, deciding nothing",
    "a POST HOC 1-Mb locus-cluster bootstrap of the held-out Delta MSE beside the registered gene-level"
    " one, labelled post hoc, deciding nothing and never replacing the registered interval: genes sharing"
    " elements are not independent, which the supervisor's audit of M2h raised, and a cluster bootstrap"
    " resamples 1-Mb bins of the genome instead of genes",
    "Delta MSE by number of rules (2-3, 4-9, 10 or more) and by whether the gene has an inhibiting rule;"
    " Delta MSE per reader cell; the predicted range quantiles",
    "the census counts: each form's fitted, informative and differing-from-the-mean genes, and the overlap"
    " with 7d's and 7h's populations",
)
BOUND = {
    "question": "QUESTION",
    "forms": "FORMS",
    "form": "FORM",
    "form_words": "FORM_WORDS",
    "form_sources": "FORM_SOURCES",
    "conventions": "CONVENTIONS",
    "premise": "PREMISE",
    "fit": "FIT",
    "structure": "STRUCTURE",
    "population": "POPULATION",
    "population_note": "POPULATION_NOTE",
    "prior_exposure": "PRIOR_EXPOSURE",
    "prior_exposure_clause": "PRIOR_EXPOSURE_CLAUSE",
    "split": "SPLIT",
    "predictions": "PREDICTIONS",
    "fit_scale": "FIT_SCALE",
    "metric": "METRIC",
    "form_interval": "FORM_INTERVAL",
    "form_interval_words": "FORM_INTERVAL_WORDS",
    "family": "FAMILY",
    "bootstrap": "BOOTSTRAP",
    "control": "CONTROL",
    "control_note": "CONTROL_NOTE",
    "permutations": "PERMUTATIONS",
    "control_alpha": "CONTROL_ALPHA",
    "min_genes": "MIN_GENES",
    "harness": "HARNESS",
    "tolerance": "TOLERANCE",
    "which_governs": "WHICH_GOVERNS",
    "not_vacuous": "NOT_VACUOUS",
    "scale_clauses": "SCALE_CLAUSES",
    "outcomes": "OUTCOMES",
    "decision": "DECISION",
    "decision_before_the_census": "DECISION_BEFORE_THE_CENSUS",
    "inconclusive": "INCONCLUSIVE",
    "falsifier": "FALSIFIER",
    "honest_size": "HONEST_SIZE",
    "stop_rule": "STOP_RULE",
    "limitations": "LIMITATIONS",
    "descriptives": "DESCRIPTIVES",
    "exploration_preceded_this": "EXPLORATION_PRECEDED_THIS",
}


def registered_fields() -> dict[str, Any]:
    """The BOUND fields, as the registration writes them and the run re-checks them."""
    out: dict[str, Any] = {}
    for field_name, name in BOUND.items():
        v = globals()[name]
        out[field_name] = list(v) if isinstance(v, tuple) else v
    return out


# ==== the score (magcal's machinery; only the percentile pair is this lane's) =============================


def bootstrap_delta(rows: list[dict], scale: float, resamples: int, seed: int) -> tuple[float, float]:
    """magcal.bootstrap_delta's body with this registration's FORM_INTERVAL (BOOTSTRAP['interval']).

    The resampling statistic is magcal's `_gene_terms` unchanged; at magcal.FORM_INTERVAL this function
    reproduces magcal.bootstrap_delta exactly (tests/test_maxform.py)."""
    import numpy as np

    s, n = mc._gene_terms(rows, scale)
    sa, na = np.array(s), np.array(n)
    rng = np.random.default_rng(seed)
    stats = np.empty(resamples)
    for i in range(resamples):
        idx = rng.integers(0, len(rows), len(rows))
        stats[i] = sa[idx].sum() / na[idx].sum()
    lo, hi = np.percentile(stats, list(FORM_INTERVAL))
    return float(lo), float(hi)


def bootstrap_scale(rows: list[dict], resamples: int, seed: int) -> tuple[float, float]:
    """magcal.bootstrap_scale's body with this registration's FORM_INTERVAL (BOOTSTRAP['interval'])."""
    import numpy as np

    dps, pps = [], []
    for r in rows:
        dp = pp = 0.0
        for p, d in zip(r["P"], r["D"], strict=True):
            dp += p * d
            pp += p * p
        dps.append(dp)
        pps.append(pp)
    dpa, ppa = np.array(dps), np.array(pps)
    rng = np.random.default_rng(seed)
    stats = np.empty(resamples)
    for i in range(resamples):
        idx = rng.integers(0, len(rows), len(rows))
        den = ppa[idx].sum()
        stats[i] = dpa[idx].sum() / den if den else np.nan
    ok = stats[~np.isnan(stats)]
    lo, hi = np.percentile(ok, list(FORM_INTERVAL))
    return float(lo), float(hi)


def cluster_bootstrap_delta(
    rows: list[dict], scale: float, resamples: int, seed: int, span: int = 1_000_000
) -> dict[str, Any]:
    """POST HOC (DESCRIPTIVES), deciding nothing: the held-out Delta MSE resampled over 1-Mb locus bins
    instead of genes, because genes sharing elements are not independent. Never replaces the registered
    gene-level interval."""
    import numpy as np

    s, n = mc._gene_terms(rows, scale)
    keys = [f"{r.get('chrom', '?')}:{int(r.get('start', 0)) // span}" for r in rows]
    order: dict[str, int] = {}
    for k in keys:
        order.setdefault(k, len(order))
    groups: list[list[int]] = [[] for _ in order]
    for i, k in enumerate(keys):
        groups[order[k]].append(i)
    gs = np.array([sum(s[i] for i in g) for g in groups])
    gn = np.array([sum(n[i] for i in g) for g in groups])
    rng = np.random.default_rng(seed)
    stats = np.empty(resamples)
    for i in range(resamples):
        idx = rng.integers(0, len(groups), len(groups))
        stats[i] = gs[idx].sum() / gn[idx].sum()
    lo, hi = np.percentile(stats, list(FORM_INTERVAL))
    return {
        "post_hoc": True,
        "decides": "nothing",
        "span_bp": span,
        "clusters": len(groups),
        "ci": [float(lo), float(hi)],
    }


def head_to_head(
    own: dict[str, dict], other: dict[str, dict], resamples: int, seed: int, scales: tuple[float, float]
) -> dict[str, Any]:
    """DESCRIPTIVE (DECISION 6), deciding nothing: MSE(other) - MSE(own) over the genes both forms run, at
    the given scales, with a percentile interval at FORM_INTERVAL."""
    import numpy as np

    common = sorted(set(own) & set(other))
    if len(common) < MIN_GENES:
        return {"genes": len(common), "compared": False}
    a_own, a_other = scales

    def terms(rs: dict[str, dict], scale: float) -> Any:
        return np.array(
            [
                sum(
                    2.0 * scale * p * d - scale * scale * p * p
                    for p, d in zip(rs[g]["P"], rs[g]["D"], strict=True)
                )
                for g in common
            ]
        )

    diff = terms(own, a_own) - terms(other, a_other)  # per gene: SSE(other) - SSE(own)
    n = np.array([len(own[g]["D"]) for g in common], dtype=float)
    rng = np.random.default_rng(seed)
    stats = np.empty(resamples)
    for i in range(resamples):
        idx = rng.integers(0, len(common), len(common))
        stats[i] = diff[idx].sum() / n[idx].sum()
    lo, hi = np.percentile(stats, list(FORM_INTERVAL))
    return {
        "genes": len(common),
        "pairs": int(n.sum()),
        "compared": True,
        "descriptive": True,
        "decides": "nothing",
        "scales": [a_own, a_other],
        "mse_other_minus_mse_own": float(diff.sum() / n.sum()),
        "ci": [float(lo), float(hi)],
        "clause": mr.HEAD_TO_HEAD["closer" if lo > 0 else "further" if hi < 0 else "same"],
    }
