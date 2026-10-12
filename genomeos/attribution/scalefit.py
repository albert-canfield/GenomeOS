# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-scalefit (ROADMAP §5, area I, 2026-10-12): do 7d's unfittable genes become fittable at M2h's scale?

7d/7e (`lane-combine2`, `948cfc9`; census `data/results/combine2_census.json`) excluded 12,500 of the
17,983 multi-element genes eligible for a fit with the reason `basal_not_positive`: under the mean form
the summed single-removal drops reach or exceed the gene's whole intact level (q >= 1), so no basal rate
is left and the gene has no fit at all. M2h (`lane-magcal`, `ac5580b`) separately measured that the
compiled rules run about five times too large: one scale per form, fitted on training genes and held out
by sha256 of the symbol, read +0.2146 [98.33% CI +0.2007, +0.2283] for the multiplicative composition.

The question this module answers is arithmetic and predictor-side only: if every single-removal drop is
multiplied by M2h's fixed scale, how many of the excluded genes have a positive basal rate again? NO RNA
VALUE IS OPENED, NOTHING IS FITTED, AND M2h's SCALE IS TAKEN VERBATIM FROM `ac5580b` AND NEVER
RE-ESTIMATED. A recovered gene is a gene 7d's fit can now accept, nothing more; this module makes no
claim about expression, about any form being right, or about any gene's level.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from genomeos.attribution import bridge, combine2

LANE = "lane-scalefit"

# ==== the scale, taken fixed =============================================================================

#: M2h's fitted scale for the multiplicative composition, copied from `data/results/magcal_result.json`
#: (`scales.multiplicative`, commit `ac5580b`) at full precision. IT IS A CONSTANT HERE. Nothing in this
#: module fits, refits, re-estimates, tunes or selects it, and no measured value is read to check it.
SCALE = 0.2145911381525326
SCALE_ROUNDED = 0.2146
SCALE_CI = (0.20074141299120243, 0.2283000359194418)
SCALE_SOURCE = (
    "data/results/magcal_result.json, field scales.multiplicative, commit ac5580b (lane-magcal, M2h,"
    " ROADMAP §5 item 7j): the one scale per form fitted on 8,181 training genes and held out by"
    " int(sha256(symbol)[:8], 16) parity, reading +0.2146 [98.33% CI +0.2007, +0.2283] with a held-out"
    " Delta MSE of +0.04339 [98.33% CI +0.0379, +0.04894]. TAKEN FIXED, NOT REFITTED, NOT RE-ESTIMATED"
)

#: The exclusion threshold on the UNSCALED ratio. A gene recovers exactly when its unscaled q is below it.
RECOVERY_THRESHOLD = 1.0 / SCALE


# ==== the rule, committed before the count ===============================================================

RULE = {
    "population": (
        "the genes 7d excluded with the reason `basal_not_positive`, read first-hand from 7d's committed"
        " census `data/results/combine2_census.json`, field reader_selection.forms.<form>.not_runnable."
        " The exclusion is computed in genomeos/attribution/combine2.py:fit as `if not q < 1.0: return"
        " Fit(form, False, 'basal_not_positive', ...)`, where for the mean q = sum over activating"
        " elements of d_i = 1 - RHO_i and for promoter competition q = Z * THETA with Z the one positive"
        " root of sum d_i (1+K)^2 / (Z + d_i (1+K)) = K. The primary population is the MEAN form's, the"
        " runtime default and 7d's 'roughly 13,000'; competition is reported beside it, not pooled"
    ),
    "the_scaling": (
        "every single-removal observation of every element of the gene is scaled by the fixed constant"
        " SCALE: RHO_i -> 1 - SCALE * (1 - RHO_i). This multiplies an activating drop d_i = 1 - RHO_i by"
        " SCALE and an inhibiting lift e_j = RHO_j - 1 by SCALE, in one expression and with no sign"
        " special case. The scaled observations are passed to combine2.fit UNCHANGED, so the recovery"
        " test is the same code 7d applied, at THETA = 0.5 and bridge.H_INTACT = 0.5"
    ),
    "recovered": (
        "PRIMARY. A gene counts as recovered when, under the same test 7d applied (combine2.fit's"
        " `q < 1.0`), its basal rate is positive on the scaled observations where it was not on the"
        " unscaled ones. Because q is exactly homogeneous of degree 1 in the drops under both forms"
        " (proved and tested: for the mean q = sum d_i; for competition substituting Z -> SCALE * Z"
        " leaves the root equation unchanged), this is equivalent to the unscaled q being below"
        " RECOVERY_THRESHOLD = 1 / SCALE = 4.6600, and the count is computed by calling combine2.fit, not"
        " by that identity"
    ),
    "runnable": (
        "SECONDARY. Of the recovered genes, how many pass 7d's FULL runnable test, i.e. also have every"
        " fitted strength in (0, 1] and so are not excluded at `inhibitor_strength_above_1`. A gene that"
        " recovers its basal rate but fails this is NOT fittable, and is reported separately"
    ),
    "informative": (
        "DESCRIPTIVE, deciding nothing. Of the recovered-and-runnable genes, how many have a predicted"
        " range above combine2.INFORMATIVE_LOG2 across the thirteen reader cells, which is 7d's further"
        " filter from runnable to run"
    ),
    "the_ratio": (
        "the ratio the brief asks for, sum of drops / gene level, IS q under the mean form, because the"
        " gene's intact level is its rate T and q is the element-driven share of it (combine2's FIT:"
        " basal b = T(1 - q)). It is reported as quantiles over the excluded population before scaling"
        " and after, plus the share at or above 1 in each, for both forms"
    ),
    "nothing_is_fitted": (
        "SCALE is a module constant copied from ac5580b. No measured RNA value, no expression table and"
        " no model request enters this census; no number in it is estimated from data this lane reads"
    ),
}

#: The floor, declared before the count. Below it the answer is "the magnitude error does not explain the
#: exclusions".
FLOOR = 6_251
FLOOR_SHARE = 0.5
FLOOR_REASON = (
    "MORE THAN HALF of the mean form's 12,500 `basal_not_positive` genes, i.e. 6,251. The hypothesis is"
    " about the POPULATION: that the same magnitude error M2h measured is what makes those genes"
    " unfittable. A scale that leaves the majority of them still unfittable has not explained the"
    " population — the remainder would still need a different cause, and the recovered minority would be"
    " the genes whose drops happened to be small, which is a statement about the tail and not about the"
    " magnitude error. Half is also the point at which the answer stops depending on where the threshold"
    " falls in the tail: the excluded genes all have q >= 1 and recovery needs q < 4.66, so a floor at"
    " half asks that the bulk of the excluded mass sit inside one scale-width of the boundary"
)

#: A second, stricter bar, also declared before the count. It changes the WORDING, never the verdict.
STRONG_BAR = 10_000
STRONG_BAR_REASON = (
    "80% of 12,500. At or above it the magnitude error is THE explanation of the exclusions; between the"
    " floor and it, a PARTIAL explanation with a named residue that needs another cause; below the floor,"
    " not an explanation. The verdict is decided by FLOOR alone"
)

#: Read before the floor was set, and therefore disclosed: 7d's census already publishes
#: reader_selection.summed_single_removal_drop_quantiles over all 17,983 eligible genes (0.1 0.4353,
#: 0.25 0.8345, 0.5 1.8518, 0.75 4.0577, 0.9 6.7135). Those are the ratio's quantiles over the eligible
#: population, not the excluded one, but they bound the answer loosely, and the 0.75 quantile sitting
#: below 4.66 while the 0.9 sits above it indicated BEFORE the floor was chosen that a floor at half
#: would very likely be cleared. The floor is therefore NOT blind, and its value as a blind test is
#: weakened accordingly. It was not raised in response, because the reasons above are the reasons for
#: half, and raising it after the exposure would be the same choice made invisible.
PRIOR_EXPOSURE = (
    "7d's committed census was read to FIND the population and to record its exact count, and the same"
    " census publishes the summed-drop quantiles over all 17,983 eligible genes (median 1.8518, 0.75"
    " 4.0577, 0.9 6.7135). Those quantiles indicated, before the floor was chosen, that a floor at half"
    " of the excluded genes would very likely be cleared, so this floor is not blind. It was kept at half"
    " for the reasons in FLOOR_REASON rather than raised, and STRONG_BAR was declared in the same breath"
    " so that the stronger reading also has a bar it had to clear in advance. No outcome of THIS census"
    " was computed before this rule was committed"
)

LIMITATIONS = (
    "M2h's scale was fitted for the MULTIPLICATIVE composition's predictor, not for the mean's or"
    " competition's. M2h's own mean scale read +0.05032 [+0.02144, +0.08224] and its competition scale"
    " +0.007315 [-0.02249, +0.03855]. Applying +0.2146 to the mean's and competition's drops is a"
    " TRANSFER ACROSS FORMS, and it is the brief's instruction, not something M2h licenses. At the mean's"
    " own scale the threshold would be 1/0.05032 = 19.9 and more genes would recover; at competition's,"
    " 136.7. The constant is the one the brief names and the only one used",
    "recovery is arithmetic. It says the fit's basal rate is positive again, so 7d's code would accept the"
    " gene. It does not say the gene's compiled rules are right, that its drops are right at this scale,"
    " or anything whatever about its measured expression. NO EXPRESSION CLAIM IS MADE HERE",
    "scaling every drop by one constant is M2h's model of the magnitude error: one number for the whole"
    " form. If the error is instead per gene, per element or size-dependent, the recovered count is not"
    " the number of genes a correct magnitude would recover",
    "the scaled drops are NOT the drops the bridge observed. The bridge's principle is that a run"
    " reproduces each element's own removal observation exactly (combine2's FIT['principle']); a scaled"
    " run reproduces 0.2146 of it. A gene that recovers here is fittable only under a bridge that admits"
    " the scale, which is M2h's finding and not 7d's convention",
    "the inhibitor-strength and activator-strength tests are scale-INVARIANT (s_j = e_j / (h (q + e_j))"
    " and s_i = w_i / max w are both ratios of quantities that scale alike), so the secondary count can"
    " only lose genes that 7d never reached, never gain any. This is a property of the algebra, tested",
)


# ==== the computation ====================================================================================


def scaled_rho(rho: dict[str, float], scale: float = SCALE) -> dict[str, float]:
    """Every single-removal observation with its effect multiplied by `scale` (RULE['the_scaling'])."""
    return {el: 1.0 - scale * (1.0 - r) for el, r in rho.items()}


def ratio(rho: dict[str, float]) -> float:
    """Sum of single-removal drops divided by the gene's intact level: the mean form's q exactly."""
    return sum(1.0 - r for r in rho.values() if r < 1.0)


@dataclass
class Recovery:
    """One excluded gene under one form, before and after the scaling."""

    form: str
    before_reason: str
    after_reason: str
    q_before: float
    q_after: float
    recovered: bool  # the basal rate is positive again (PRIMARY)
    runnable: bool  # and every fitted strength is in (0, 1] (SECONDARY)

    def as_row(self) -> dict[str, Any]:
        return {
            "form": self.form,
            "before_reason": self.before_reason,
            "after_reason": self.after_reason,
            "q_before": self.q_before,
            "q_after": self.q_after,
            "recovered": self.recovered,
            "runnable": self.runnable,
        }


def recovery(form: str, rho: dict[str, float], scale: float = SCALE) -> Recovery:
    """7d's fit on the observations, then on the scaled ones, with the same code and conventions."""
    before = combine2.fit(form, rho)
    after = combine2.fit(form, scaled_rho(rho, scale))
    was_excluded = (not before.ok) and before.reason == "basal_not_positive"
    return Recovery(
        form=form,
        before_reason="" if before.ok else before.reason,
        after_reason="" if after.ok else after.reason,
        q_before=before.q,
        q_after=after.q,
        recovered=was_excluded and after.q < 1.0,
        runnable=was_excluded and after.ok,
    )


def q_is_homogeneous(form: str, rho: dict[str, float], scale: float = SCALE) -> float:
    """|q(scaled) - scale*q(unscaled)|: 0 up to floating point, for both forms (RULE['recovered'])."""
    before, after = combine2.fit(form, rho), combine2.fit(form, scaled_rho(rho, scale))
    if math.isnan(before.q) or math.isnan(after.q):
        return math.nan
    return abs(after.q - scale * before.q)


def quantiles(xs: list[float], qs: tuple[float, ...] = (0.1, 0.25, 0.5, 0.75, 0.9, 0.99)) -> dict[str, float]:
    s = sorted(xs)
    return {str(q): s[min(len(s) - 1, int(q * len(s)))] for q in qs} if s else {}


def verdict(recovered: int, excluded: int) -> str:
    """The registered wording, decided by FLOOR alone (RULE, FLOOR_REASON, STRONG_BAR_REASON)."""
    if recovered < FLOOR:
        return (
            f"the magnitude error does not explain the exclusions: scaling every single-removal drop by"
            f" M2h's fixed {SCALE_ROUNDED} leaves {excluded - recovered} of the {excluded} genes 7d"
            f" excluded as basal_not_positive still without a positive basal rate, and recovers"
            f" {recovered}, below the floor of {FLOOR} declared before the count"
        )
    if recovered >= STRONG_BAR:
        return (
            f"the magnitude error explains the exclusions: scaling every single-removal drop by M2h's"
            f" fixed {SCALE_ROUNDED} gives {recovered} of the {excluded} genes 7d excluded as"
            f" basal_not_positive a positive basal rate again, clearing both the floor of {FLOOR} and the"
            f" strong bar of {STRONG_BAR} declared before the count. A recovered gene is a gene the fit"
            f" can now accept and nothing more: this says nothing about expression"
        )
    return (
        f"the magnitude error partly explains the exclusions: scaling every single-removal drop by M2h's"
        f" fixed {SCALE_ROUNDED} gives {recovered} of the {excluded} genes 7d excluded as"
        f" basal_not_positive a positive basal rate again, clearing the floor of {FLOOR} but not the"
        f" strong bar of {STRONG_BAR}; the remaining {excluded - recovered} need another cause. A"
        f" recovered gene is a gene the fit can now accept and nothing more: this says nothing about"
        f" expression"
    )


# the two engine constants this census depends on, recorded so a change is visible
CONVENTIONS = {
    "theta": combine2.THETA,
    "h_intact": bridge.H_INTACT,
    "strength_tolerance": bridge.STRENGTH_TOLERANCE,
    "informative_log2": combine2.INFORMATIVE_LOG2,
    "forms": list(combine2.FORMS),
}
