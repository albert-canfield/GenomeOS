# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-residue (ROADMAP §5, area I, 2026-10-12): why the residue stays unfittable at M2h's scale.

`lane-scalefit` (rule and floor `8469df9`, census `7f65544`) measured that the mean form excluded 12,500
of 17,983 multi-element genes as `basal_not_positive` — the sum of their single-removal drops reaches or
exceeds the gene's whole intact level, so no basal rate is left — and that 8,792 of those get a positive
basal rate back when every drop is multiplied by M2h's fixed scale 0.2146. The remainder, which that
census registered as needing another cause, is the population of THIS lane: the genes whose summed drops
still exceed the gene's level even at the scale.

The question is exactly one: **is the excess carried by one or a few elements, or is it spread across
many?** A cap on how much any single element may claim of a gene's level can only help in the first case.
This module declares the statistics that answer it, the per-element cap levels at which the recovery is
counted, and the thresholds that map the measurement onto the implication — ALL BEFORE ANY COUNT IS
COMPUTED, so that no summary statistic and no cap level is chosen after the distribution has been seen.

PREDICTOR-SIDE ONLY. NO RNA VALUE IS OPENED, NOTHING IS FITTED, TUNED, SELECTED OR RE-ESTIMATED, no model
request is sent and no key is read. M2h's scale is taken verbatim from `genomeos.attribution.scalefit`,
which copied it from `ac5580b`. This lane PROPOSES NO FORM, registers no RNA test, and makes no claim
whatever about expression: a gene a fit can accept is not a gene the model predicts, which `lane-recovered`
(`f775e9b`) had to establish as a separate registration on separate data.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from genomeos.attribution import bridge, combine2, scalefit

LANE = "lane-residue"

#: M2h's scale, taken from `scalefit` and NOT re-declared here, so the two lanes cannot drift apart.
SCALE = scalefit.SCALE
SCALE_ROUNDED = scalefit.SCALE_ROUNDED

#: A gene's unscaled q at or above this still exceeds its level after the scaling: 1 / SCALE = 4.660025.
RESIDUE_THRESHOLD = scalefit.RECOVERY_THRESHOLD


# ==== the questions, committed BEFORE they are computed ==================================================

QUESTION = (
    "Why can the mean form not fit the residue of `lane-scalefit`'s census -- the genes 7d excluded as"
    " basal_not_positive whose summed single-removal drops STILL reach or exceed the gene's whole intact"
    " level after every drop is multiplied by M2h's fixed 0.2146 -- and is the excess carried by one or a"
    " few elements, or spread across many? The answer decides what the next combination form should be: a"
    " PER-ELEMENT CAP if one or a few elements dominate, something else if the sum is diffuse."
)

POPULATION = (
    "RE-DERIVED, not copied. The 24 compiled non-coding programs are gated and parsed again with 7d's own"
    " code (`combine2.multi_rule_genes`, `combine2.gene_rules`, `combine2.classify`, `combine2.fit` at"
    " THETA = 0.5 and bridge.H_INTACT = 0.5, selection `combine2.SELECTION`), the mean form's"
    " `basal_not_positive` genes are taken from that run, and the RESIDUE is those of them whose fit is"
    " still `basal_not_positive` on the scaled observations (`scalefit.scaled_rho`). The count is computed"
    " by calling `combine2.fit` twice per gene, exactly as scalefit did, never by the identity q_scaled ="
    " SCALE * q. scalefit's committed counts (12,500 excluded and its residue) are read AFTERWARDS and"
    " reconciled against this run; a disagreement is a finding of this lane and not something to repair."
    " The RECOVERED population, the excluded genes that are NOT in the residue, is measured beside it in"
    " every statistic, because a concentration number means nothing without the population it is"
    " contrasted with. The two are never pooled. The promoter-competition form is NOT this lane's"
    " population: the mean form is the runtime default and scalefit's primary, and the drop shares below"
    " are the mean form's q by construction."
)

#: The drops the statistics are over, declared so the choice is not made later.
DROPS = (
    "every statistic below is over a gene's ACTIVATING single-removal drops d_i = 1 - RHO_i for RHO_i < 1,"
    " and over no other quantity. That is the set `combine2.fit` sums into q under the mean form, so the"
    " shares below are shares of exactly the quantity that exceeds the gene's level. An inhibiting"
    " element's lift (RHO_j - 1) does not enter q and is not counted in any share; the number of"
    " inhibiting elements is recorded per gene as a descriptive column and decides nothing."
)

QUESTIONS: dict[str, str] = {
    "Q1_re_derivation": (
        "How many genes are in the residue, re-derived from the compiled programs by this lane's own run?"
        " Reported with the re-derived count of the excluded population it is a subset of, and reconciled"
        " against scalefit's committed numbers after the fact."
    ),
    "Q2_how_many_elements": (
        "For each residue gene, the number of elements (rules) and the number of activating elements."
        " Reported as the quantiles 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99 and the min, max and mean, over"
        " the residue and over the recovered beside it -- the DISTRIBUTIONS, not only the medians."
    ),
    "Q3_share_of_the_largest_element": (
        "For each residue gene, top1 = (the largest activating drop) / (the sum of the activating drops):"
        " the share of the summed drop carried by its largest single element. Quantiles as in Q2, over the"
        " residue and the recovered. This is the statistic that decides the question: a residue whose"
        " median top1 is near 1 is a residue of genes with one dominant element; a residue whose median"
        " top1 is near 1/n is a residue of genes whose elements all contribute alike."
    ),
    "Q4_share_of_the_largest_three": (
        "The same for top3 = (the sum of the three largest activating drops) / (the sum of all of them),"
        " which is 1 for a gene with three or fewer activating elements. Quantiles as in Q2, over both"
        " populations, and separately over the residue genes with four or more activating elements, for"
        " which the number is not forced to 1."
    ),
    "Q5_ratio_after_scaling": (
        "For each residue gene, the ratio after scaling, q_after = SCALE * (sum of the activating drops),"
        " which is the mean form's q on the scaled observations. All residue genes are at or above 1 by"
        " construction, so what the quantiles say is HOW FAR above: a residue sitting just above 1 is a"
        " different problem from one sitting at 5. Quantiles as in Q2, plus the count at or above 2 and at"
        " or above 5. The unscaled ratio's quantiles are reported beside it."
    ),
    "Q6_how_few_elements_carry_the_excess": (
        "For each residue gene, k_min = the smallest number of its LARGEST activating drops whose removal"
        " from the sum would bring the scaled sum of the rest below 1, i.e. the smallest k with SCALE *"
        " (sum of all drops - the k largest) < 1; and NaN if no k < n achieves it (a gene whose drops"
        " cannot be brought under the level by dropping all but one). Quantiles as in Q2, and the counts"
        " at k_min = 1, = 2, = 3, <= 3 and > 3. REMOVING an element is not a form this lane proposes and"
        " not a fit: it is a concentration statistic, the direct reading of 'carried by one or a few'."
    ),
    "Q7_recovery_under_a_per_element_cap": (
        "Under a PER-ELEMENT CAP at level L, how many residue genes get a positive basal rate back? The"
        " cap rule, the ladder of L, and the arithmetic bound on what a cap can do are in CAP_RULE,"
        " CAP_LEVELS and CAP_BOUND below, all declared before the count. Reported as the count and share"
        " recovered at EVERY L of the ladder -- the whole curve, not a chosen point -- with the primary L"
        " named in advance."
    ),
    "Q8_diffuseness": (
        "For each residue gene, n_eff = (sum of the drops)^2 / (sum of the drops squared), the inverse"
        " Herfindahl index of the drop shares: the effective number of elements the summed drop is spread"
        " over, which is 1 when one element carries all of it and n when all n contribute equally. Also"
        " n_eff / n_activating. Quantiles as in Q2, over both populations. This is Q3 and Q4's statistic"
        " read on the whole vector rather than its head, and it is declared so that a flat top1 cannot be"
        " reported as concentration or a peaked one as diffusion."
    ),
}


# ==== the per-element cap, its ladder and its bound, all declared before the count ========================

CAP_RULE = (
    "A per-element cap at level L replaces every scaled activating drop by min(SCALE * d_i, L) and leaves"
    " the inhibiting lifts untouched. The gene has a positive basal rate again -- it becomes fittable at"
    " the step 7d excluded it on -- exactly when the capped sum is below 1. L is in units of the gene's"
    " INTACT LEVEL, which is the unit q is already in, so L = 0.1 reads 'no single element may claim more"
    " than a tenth of the gene'. Capping is applied on top of M2h's scale and never instead of it, and the"
    " count is computed by summing the capped drops, with `combine2.fit` called on the capped observations"
    " as a cross-check on every gene."
)

#: The cap levels, declared BEFORE the count. A ladder, not a choice: the whole curve is reported.
CAP_LEVELS: tuple[float, ...] = (1.0, 0.5, 0.25, 0.1, 0.05, 0.02, 0.01)

#: The one level named in advance as primary, for a reason that is not about this distribution.
CAP_PRIMARY = 0.5

CAP_LEVELS_REASON = (
    "Seven round levels -- halves and decades -- spanning two orders of magnitude, from the whole of the"
    " gene's level down to a hundredth of it, fixed in this file before a single gene was counted. THE"
    " WHOLE CURVE IS REPORTED AND NO LEVEL IS SELECTED AFTERWARDS, because a cap chosen after seeing the"
    " distribution is a fitted parameter wearing a census's clothes: it would be the one number that makes"
    " the count come out well, reported as though it had been a prediction. The ladder is round decades"
    " and halves precisely so that it cannot have been tuned -- there is no free digit in it. The PRIMARY"
    " level is 0.5, for a reason about the model and not about the data: an element may claim at most half"
    " of the gene's intact level, which is the same half the engine already uses as bridge.H_INTACT = 0.5"
    " and combine2.THETA = 0.5, and it is the weakest cap that forces a multi-element gene's level to be"
    " shared. The primary decides the WORDING of Q7, never the verdict on the question, which is decided"
    " by Q3, Q4, Q6 and Q8 through DECISION_RULE."
)

CAP_BOUND = (
    "Stated in advance, from the algebra: a cap at L leaves the capped sum at most n_activating * L, so a"
    " cap at L CANNOT fit any gene with n_activating >= 1/L, and it fits EVERY gene with n_activating <"
    " 1/L whatever its drops are. L = 1.0 and L = 0.5 therefore recover at most the residue genes with one"
    " activating element -- of which there are none, the population being multi-element and the drops"
    " summing above 4.66 -- unless a gene's capped sum falls below 1 some other way, so the top of the"
    " ladder is expected to be near zero and its being zero is not a finding. The bound is the reason the"
    " ladder runs down two decades: the count must be able to go from nothing to nearly everything inside"
    " it for the curve to say anything, and 1/0.01 = 100 sits above the rule counts scalefit published for"
    " the residue (see PRIOR_EXPOSURE). A cap also destroys the bridge's exact-reproduction principle"
    " further than the scale already does: a capped element's run reproduces neither its own removal"
    " observation nor 0.2146 of it. THIS LANE PROPOSES NO FORM; it counts what a cap could and could not"
    " do, and stops."
)


# ==== the decision rule, declared before the count =======================================================

#: The thresholds that map the measurement onto the implication, so the reading is not chosen afterwards.
CONCENTRATED_TOP1 = 0.5
FEW_TOP3 = 0.75
FEW_KMIN = 3
DIFFUSE_TOP1 = 0.25

DECISION_RULE = (
    "Decided on the RESIDUE's medians, in this order, with the thresholds fixed here before the count."
    " CONCENTRATED, if the median top1 >= 0.5: the largest single element carries most of the summed drop,"
    " and a per-element cap is the indicated next form. FEW, if the median top3 >= 0.75 or the median"
    " k_min <= 3: no one element dominates but three do, and a correction on a few elements per gene --"
    " of which a cap is one -- can reach the residue. DIFFUSE, if the median top1 < 0.25 AND the median"
    " k_min > 3: what exceeds the gene is the sum over many comparable elements, no element is the"
    " culprit, and a per-element cap is NOT the indicated form. MIXED otherwise, reported as mixed and not"
    " rounded to either side. The thresholds are the plain readings of the words: a half is 'most of', a"
    " largest-of-three at three quarters is 'a few carry it', a quarter is 'the largest carries little',"
    " and three is the brief's own 'one or a few'. Q7's curve is reported in every case and changes no"
    " branch of this rule."
)

PRIOR_EXPOSURE = (
    "Disclosed, because it was read before these questions were written. scalefit's committed census"
    " (`7f65544`, data/results/scalefit_census.json) already publishes, for the residue, the RULE-COUNT"
    " quantiles 0.25 41, median 58, 0.75 84, 0.99 194, against 0.25 10, median 15, 0.75 23 for the"
    " recovered, and its commit message argues from them that 'summing per-element maxima over many"
    " elements -- not the magnitude of any one element -- is what exceeds the gene'. So Q2's answer is"
    " ALREADY KNOWN in outline and Q2 is a re-derivation and not a test, and the DIFFUSE branch of"
    " DECISION_RULE is the branch that prior reading points at. What is NOT known from it, and is what"
    " this lane measures, is how the summed drop is distributed WITHIN a gene: a gene with 58 elements can"
    " still have one element carrying 90% of the sum, and the published rule counts say nothing either"
    " way. Q3, Q4, Q6, Q7 and Q8 have no published precedent that this lane is aware of. The thresholds"
    " and the cap ladder were nevertheless chosen to be round and symmetric rather than tuned to make the"
    " expected branch come out, and the ladder's lower end was set at 0.01 because 1/0.01 = 100 clears the"
    " published 0.99 rule count of 194 only partly -- a level chosen to be round, not to be flattering."
)

LIMITATIONS = (
    "the residue is defined by M2h's scale for the MULTIPLICATIVE composition, transferred to the mean"
    " form. M2h's own mean scale read +0.05032, at which the threshold would be 19.9 instead of 4.66 and"
    " the residue would be much smaller. The transfer is scalefit's registered choice and this lane"
    " inherits it unchanged rather than re-opening it; every count here is conditional on it",
    "fittable is not predicted. Every count in this census is about whether 7d's fit has a positive basal"
    " rate left, which is a property of the arithmetic. Whether any of these genes would then read closer"
    " to measured RNA is a SEPARATE registration on separate data, as `lane-recovered` had to be for the"
    " genes that recovered at the scale alone. NO RNA VALUE IS OPENED IN THIS LANE AND NO EXPRESSION CLAIM"
    " IS MADE",
    "a cap is counted, not proposed. Q7 measures what a per-element cap could do to the residue's"
    " arithmetic. It does not say a cap is right, it does not say the capped drops are the drops the bridge"
    " observed -- they are not -- and this lane registers no form and no test of one",
    "the shares are of the compiled drops, which are themselves predicted rules carrying whatever error"
    " the compilation carries. A concentration measured on wrong drops is a concentration of the wrong"
    " drops. This census measures the predictor as it stands and nothing behind it",
    "k_min and n_eff are concentration statistics computed on one gene at a time. Neither is a model, and"
    " neither licenses deleting or reweighting any element in any run",
)


# ==== the computation ====================================================================================


def activating_drops(rho: dict[str, float]) -> list[float]:
    """The activating single-removal drops d_i = 1 - RHO_i, largest first (DROPS)."""
    return sorted((1.0 - r for r in rho.values() if r < 1.0), reverse=True)


def top_share(drops: list[float], k: int = 1) -> float:
    """The share of the summed drop carried by the `k` largest elements (Q3, Q4)."""
    total = sum(drops)
    if not drops or total <= 0.0:
        return math.nan
    return sum(drops[:k]) / total


def n_effective(drops: list[float]) -> float:
    """(sum d)^2 / sum d^2: the effective number of elements the summed drop is spread over (Q8)."""
    total, sq = sum(drops), sum(x * x for x in drops)
    if not drops or sq <= 0.0:
        return math.nan
    return total * total / sq


def k_min_to_fit(drops: list[float], scale: float = SCALE) -> float:
    """Smallest k with scale * (sum of all drops but the k largest) < 1, NaN if no k < n does (Q6)."""
    rest = sum(drops)
    for k in range(len(drops)):
        if scale * rest < 1.0:
            return float(k)
        rest -= drops[k]
    return math.nan


def capped_sum(drops: list[float], cap: float, scale: float = SCALE) -> float:
    """The summed scaled drop with every element capped at `cap` of the gene's level (CAP_RULE)."""
    return sum(min(scale * d, cap) for d in drops)


def fits_under_cap(drops: list[float], cap: float, scale: float = SCALE) -> bool:
    """Whether the capped sum leaves a positive basal rate, the test 7d excluded the gene on."""
    return capped_sum(drops, cap, scale) < 1.0


def capped_rho(rho: dict[str, float], cap: float, scale: float = SCALE) -> dict[str, float]:
    """The scaled observations with each activating drop capped at `cap`; lifts untouched (CAP_RULE).

    Only used to cross-check `fits_under_cap` through `combine2.fit` itself, never to propose a form.
    """
    out: dict[str, float] = {}
    for el, r in rho.items():
        if r < 1.0:
            out[el] = 1.0 - min(scale * (1.0 - r), cap)
        else:
            out[el] = 1.0 + scale * (r - 1.0)
    return out


QUANTILES: tuple[float, ...] = (0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)


def distribution(xs: list[float]) -> dict[str, Any]:
    """The quantiles Q2 declares, plus min, max, mean and n -- the distribution, never only a median."""
    finite = sorted(x for x in xs if not math.isnan(x))
    out: dict[str, Any] = {"n": len(finite), "n_nan": len(xs) - len(finite)}
    if not finite:
        return out
    out["quantiles"] = {str(q): finite[min(len(finite) - 1, int(q * len(finite)))] for q in QUANTILES}
    out["min"], out["max"] = finite[0], finite[-1]
    out["mean"] = sum(finite) / len(finite)
    return out


def median(xs: list[float]) -> float:
    """The 0.5 quantile exactly as `distribution` computes it, so the decision reads the reported number."""
    finite = sorted(x for x in xs if not math.isnan(x))
    return finite[min(len(finite) - 1, int(0.5 * len(finite)))] if finite else math.nan


def decide(median_top1: float, median_top3: float, median_kmin: float) -> str:
    """The branch of DECISION_RULE the residue's medians fall in. Declared before the count."""
    if not math.isnan(median_top1) and median_top1 >= CONCENTRATED_TOP1:
        return "concentrated"
    few = (not math.isnan(median_top3) and median_top3 >= FEW_TOP3) or (
        not math.isnan(median_kmin) and median_kmin <= FEW_KMIN
    )
    if few:
        return "few"
    diffuse = (not math.isnan(median_top1) and median_top1 < DIFFUSE_TOP1) and (
        math.isnan(median_kmin) or median_kmin > FEW_KMIN
    )
    return "diffuse" if diffuse else "mixed"


READINGS: dict[str, str] = {
    "concentrated": (
        "the excess is CARRIED BY ONE ELEMENT: over the residue the median gene's largest activating"
        " single-removal drop is at least half of its summed drop, so what exceeds the gene's level is one"
        " element's size and not the number of elements. A PER-ELEMENT CAP is the indicated next"
        " combination form."
    ),
    "few": (
        "the excess is CARRIED BY A FEW ELEMENTS: no single element carries most of the residue gene's"
        " summed drop, but its three largest do, or removing at most three of them would bring the sum"
        " under the gene's level. A correction that acts on a few elements per gene, of which a"
        " per-element cap is one, can reach the residue."
    ),
    "diffuse": (
        "the excess is SPREAD ACROSS MANY ELEMENTS: over the residue the median gene's largest activating"
        " drop is under a quarter of its summed drop and more than three elements would have to be removed"
        " to bring the sum under the gene's level. What exceeds the gene is the SUM over many comparable"
        " elements, no element is the culprit, and a per-element cap is NOT the indicated next form."
    ),
    "mixed": (
        "the excess is NEITHER clearly concentrated nor clearly diffuse by the thresholds declared before"
        " the count: the residue's medians fall between them. Reported as mixed and not rounded to either"
        " side; what the next form should be is not decided by this census."
    ),
}


def reading(branch: str, **numbers: Any) -> str:
    """The registered wording for a branch, with the numbers appended and the wording never softened."""
    tail = ", ".join(f"{k} {v}" for k, v in numbers.items())
    return f"{READINGS[branch]} ({tail})" if tail else READINGS[branch]


@dataclass
class GeneResidue:
    """One excluded gene's concentration statistics, residue or recovered (POPULATION)."""

    symbol: str
    chrom: str
    n_rules: int
    n_activating: int
    n_inhibiting: int
    sum_d: float
    q_after: float
    in_residue: bool
    top1: float
    top3: float
    n_eff: float
    k_min: float
    fits_under: dict[str, bool]

    def as_row(self) -> dict[str, Any]:
        return {
            "symbol": self.symbol,
            "chrom": self.chrom,
            "n_rules": self.n_rules,
            "n_activating": self.n_activating,
            "n_inhibiting": self.n_inhibiting,
            "sum_d": self.sum_d,
            "q_after": self.q_after,
            "in_residue": self.in_residue,
            "top1": self.top1,
            "top3": self.top3,
            "n_eff": self.n_eff,
            "k_min": self.k_min,
            "fits_under": self.fits_under,
        }


def gene_residue(symbol: str, chrom: str, rho: dict[str, float], in_residue: bool) -> GeneResidue:
    """Every declared per-gene statistic for one excluded gene, from its removal folds alone."""
    drops = activating_drops(rho)
    return GeneResidue(
        symbol=symbol,
        chrom=chrom,
        n_rules=len(rho),
        n_activating=len(drops),
        n_inhibiting=sum(1 for r in rho.values() if r > 1.0),
        sum_d=sum(drops),
        q_after=SCALE * sum(drops),
        in_residue=in_residue,
        top1=top_share(drops, 1),
        top3=top_share(drops, 3),
        n_eff=n_effective(drops),
        k_min=k_min_to_fit(drops),
        fits_under={str(c): fits_under_cap(drops, c) for c in CAP_LEVELS},
    )


# the engine constants this census depends on, recorded so a change is visible
CONVENTIONS = {
    "theta": combine2.THETA,
    "h_intact": bridge.H_INTACT,
    "strength_tolerance": bridge.STRENGTH_TOLERANCE,
    "forms": list(combine2.FORMS),
    "scale": SCALE,
    "scale_is_refitted": False,
    "residue_threshold_on_the_unscaled_ratio": RESIDUE_THRESHOLD,
}
