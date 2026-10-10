# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""lane-rootpick (ROADMAP area A, from lane-capsumrna's census of 2026-10-10): which root the capped
sum's fit takes, decided from the definition and the bridge and NOT from which root fits better.

THE AMBIGUITY, as lane-capsumrna measured and registered it. `runtime/grn.py combine_activators`
carries `sum_capped`: A = min(1, W), W = sum of s_i*H(x_i). Below the cap that is linear in the element
terms, so the capped sum gives the mean's level there. Where four or more activating elements make
H*(k-1) > 1 and D < top/H, a SECOND solution exists with the cap binding -- an element-driven share
q = Z greater than the uncapped root's q = D. In that lane's own words, kept verbatim because it is the
registered reading and must not be strengthened:

    "both reproduce every removal and T at the same conventions, so the fit is not unique
    (fit_not_unique). The capped root's prediction differs from the mean's; it is counted and run,
    never chosen, because the definition does not choose between the roots"

So the observables the fit is built from -- every removal, and T -- do not distinguish the roots. That
is genuine underdetermination and it was measured, not assumed. The previous lane stopped there and
counted, which was right: it had no ground to choose, and it said so.

THIS LANE SUPPLIES THE GROUND, AND IT IS NOT IN min(1, W). The definition alone really does not choose;
anyone looking only at the form will stay stuck where that lane stopped. The ground is in the BRIDGE,
which states in its own words what the strength convention does and does not license
(`attribution/bridge.py`, DECLARED_STRENGTH and DECLARED_STRENGTH_LICENSES, quoted in GROUND below):

  - "V and s are NOT separately identified: only their product is. DECLARED_STRENGTH fixes the split by
    convention and by nothing else."
  - "what DECLARED_STRENGTH does NOT license: it is not a measurement that the element saturates its
    gene."
  - it licenses "the absolute scale of one gene's transcription and its basal fraction ... not the
    ceiling, not the strength, and not the response to any perturbation but this removal"

THE CAPPED ROOT IS A SATURATION CLAIM. It holds only when the summed drive reaches the ceiling, W > 1.
The bridge says the measurements used here license "not the ceiling" and are "not a measurement that the
element saturates its gene". So the capped root asserts exactly the quantity the bridge names as
unlicensed, while the uncapped root's q = D is the share the declared split does fix. The roots are
therefore NOT symmetric, even though the fit cannot tell them apart.

WHAT THE RULE CLAIMS, AND IT IS NARROW. It does not say the cap never binds in nature, nor that the
capped root is numerically wrong. It says this evidence cannot place a gene on the capped branch,
because placing it there requires a ceiling the bridge does not identify.
"""

from __future__ import annotations

from typing import Any

#: The decision, fixed before any RNA or outcome value was read by this lane.
ROOT_RULE = (
    "where both roots exist, the fit takes the UNCAPPED root, whose element-driven share q = D is the"
    " one the declared split fixes; where only the CAPPED root exists, the gene is UNKNOWN and leaves"
    " the population, because its single fit rests on a ceiling the bridge does not identify. The cap"
    " is a ceiling and never a target: a fit equally consistent with an unbound cap does not assert"
    " saturation"
)

#: Why, in the bridge's own words. Quoted and not paraphrased, because the whole rule rests on them.
GROUND = (
    "attribution/bridge.py: 'V and s are NOT separately identified: only their product is."
    " DECLARED_STRENGTH fixes the split by convention and by nothing else'; 'what DECLARED_STRENGTH"
    " does NOT license: it is not a measurement that the element saturates its gene'; and"
    " DECLARED_STRENGTH_LICENSES = 'the absolute scale of one gene's transcription and its basal"
    " fraction, in molecules/(cell*h); not the ceiling, not the strength, and not the response to any"
    " perturbation but this removal'. The capped root exists only when W > 1, which is a statement"
    " about the ceiling; the uncapped root's q = D is the declared split itself"
)

#: What would show the rule wrong. It is a licensing claim, so a licensing fact refutes it.
FALSIFIER = (
    "the rule is WRONG if any measurement in the bridge's own tier table identifies the ceiling or"
    " identifies s separately from V -- then the capped branch is licensed and the rule's ground"
    " dissolves. It is also wrong if a gene's two roots can be distinguished by an observable the fit"
    " already uses (a removal or T), because then the fit was not underdetermined and no rule was"
    " needed. Neither is a question about which root fits the RNA better, and a better fit on the"
    " capped branch would NOT refute it"
)

#: The branch where this lane decides nothing.
INCONCLUSIVE = (
    "if the bridge's tiers neither identify nor exclude the ceiling, the rule is not established and"
    " the genes stay as lane-capsumrna left them: counted, run, never chosen. This lane then reports"
    " 'the ground is not in the record' and changes no classification"
)

#: Read before the rule was fixed, so it is exposure and is disclosed rather than claimed as blind.
PRIOR_EXPOSURE = (
    "this lane read, before fixing the rule: runtime/grn.py's combine_activators; "
    "attribution/capsumrna.py's docstring, FORM, fit and STRUCTURE strings; attribution/bridge.py's"
    " DECLARED_STRENGTH block; and the PREDICTOR-SIDE counts of data/results/capsumrna_census.json,"
    " which that file states are 'a count on the predictor side, before any RNA value is read'."
    " NO RNA value, no measured level, no fit residual and no outcome of any kind was read. The"
    " coordinator's brief supplied the figures 126 and 3 before any of this"
)

#: Stated here because it is the strongest evidence the rule was not chosen for its effect.
FORECLOSES = (
    "taking the CAPPED root on the ambiguous genes would have given 125 informative capped roots plus"
    " the 3 capped_root_only genes = 128, ABOVE lane-capsumrna's MIN_GENES gate of 100, and would have"
    " made the capped sum a testable form against RNA. This rule gives ZERO genes whose prediction"
    " differs from the mean, so it FORECLOSES that test entirely. The rule chosen is the one that"
    " loses the result; that asymmetry is the point and is recorded before any count"
)

#: What the rule does not say, so a later reader does not widen it.
DOES_NOT_CLAIM = (
    "not that the cap never binds in nature; not that the capped root is numerically wrong; not that"
    " sum_capped is a bad form; and nothing at all about the mean, competition or multiplicative,"
    " which are single-valued by construction (mean and sum-below-cap are linear, competition"
    " W/(h+W) and multiplicative exp(gain*W) are strictly monotone in W and so invertible)"
)

#: The statuses lane-capsumrna's census assigns, and what this rule does with each.
VERDICT_BY_STATUS: dict[str, str] = {
    "runnable": "admit, uncapped root (unchanged)",
    "fit_not_unique": "admit, uncapped root (was: counted, run, never chosen)",
    "capped_root_only": "UNKNOWN, leaves the population (was: admitted as its one fit)",
    "sum_past_the_cap": "no fit at these conventions (unchanged)",
    "inhibitor_strength_above_1": "the mean's refusal (unchanged)",
    "basal_not_positive": "the mean's refusal (unchanged)",
}

#: The forms this ambiguity can touch at all, with the reason each is or is not single-valued.
FORMS = {
    "sum_capped": "AMBIGUOUS: min(1, W) is flat above W = 1, so the inverse is not a function there",
    "mean": "single-valued: sum(t)/n is linear in the terms",
    "competition": "single-valued: W/(h + W) is strictly monotone in W, so invertible",
    "multiplicative": "single-valued: exp(gain*W) is strictly monotone in W, so invertible",
    "max": "not fitted as a share in this family; outside this lane's scope and said so",
    "or": "not fitted as a share in this family; outside this lane's scope and said so",
}


def licensed_root(status: str) -> str:
    """`uncapped`, `unknown`, or `no_fit` for one census status, under ROOT_RULE."""
    if status in ("runnable", "fit_not_unique"):
        return "uncapped"
    if status == "capped_root_only":
        return "unknown"
    if status in ("sum_past_the_cap", "inhibitor_strength_above_1", "basal_not_positive"):
        return "no_fit"
    raise ValueError(f"unknown capped-sum status {status!r}; one of {sorted(VERDICT_BY_STATUS)}")


def admits(status: str) -> bool:
    """Does a gene with this status carry a capped-sum prediction under ROOT_RULE?"""
    return licensed_root(status) == "uncapped"


def differs_from_the_mean(status: str) -> bool:
    """Does its prediction differ from the mean's under ROOT_RULE?

    No status does. On the uncapped root the capped sum IS the mean (both linear below the cap), and
    the only status whose prediction differed -- capped_root_only -- is UNKNOWN here. That is the
    finding, and it is arithmetic from the rule rather than a measurement.
    """
    return False


def reclassify(status_counts: dict[str, int]) -> dict[str, Any]:
    """Apply ROOT_RULE to lane-capsumrna's committed status counts. Structural; reads no RNA."""
    unknown = {s: n for s, n in status_counts.items() if licensed_root(s) == "unknown"}
    admitted = {s: n for s, n in status_counts.items() if admits(s)}
    return {
        "admitted": sum(admitted.values()),
        "admitted_by_status": admitted,
        "unknown_under_this_rule": sum(unknown.values()),
        "unknown_by_status": unknown,
        "differs_from_the_mean": 0,
        "why_zero": (
            "the uncapped root gives the mean's level on every admitted gene, and the only status whose"
            " prediction differed from the mean's is UNKNOWN under this rule"
        ),
    }
