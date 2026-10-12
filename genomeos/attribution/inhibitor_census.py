# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-inhibitor (ROADMAP §5, area I, 2026-10-12): why do 39 genes carry an inhibitor of strength above
1, and what does that say about the COMPILED LAYER rather than about those genes?

STOP-LIST CHECK (docs/ROADMAP.md §5, "The stop list of 2026-09-27, reinstated"): no AlphaGenome request,
no deletion request, no biosample added, nothing in area J, no node-model work, no BioForge calibration.
THIS LANE READS NO RNA VALUE, FITS NOTHING, registers no test and makes no claim about expression. It is a
PREDICTOR-SIDE CENSUS: every number in it is arithmetic on the compiled rules' own stated removal
responses, and the only outside constant is M2h's scale, taken verbatim through `scalefit.SCALE`.

WHAT IS ALREADY IN THE RECORD
-----------------------------
`lane-scalefit` (`7f65544`) measured that the runtime's MEAN form excludes 12,500 of the 17,983
multi-element genes `combine2.classify` admits, as `basal_not_positive`; 8,792 recover when every drop is
scaled by M2h's fixed 0.2146 and 8,535 pass the full runnable test. `lane-residue` (`5d2c8f5`) measured
that the 3,708 residue is DIFFUSE -- no element carries over 21.0% of a gene's summed drop and a
per-element cap recovers 0 genes at L >= 0.25. `lane-residuemult` (`d2065a1`) then looked the residue up
in the MULTIPLICATIVE form and found 3,669 of the 3,708 (98.9%) already fitted at the same fixed scale,
so the residue is the ADDITIVE form's structural failure. It also recorded, in
`informative_under_neither_by_reason`, that the 39 genes fitted by NEITHER form fail on one thing and the
same thing: `inhibitor_strength_above_1`. Not one is `basal_not_positive`. Nobody has looked at why.

WHAT THE REFUSAL IS, read off the code before any count
-------------------------------------------------------
An element whose removal RAISES the gene (RHO > 1) is an inhibitor. Both fits turn that observation into a
strength through the bridge's intact Hill occupancy H = `bridge.H_INTACT` = 0.5 and then require the
strength to lie in (0, 1] (`combine2._strength_ok`, with `bridge.STRENGTH_TOLERANCE` = 1e-9 of slack):

* the MULTIPLICATIVE form, `saturate.fit`, uses s = (1 - 1/RHO) / H -- the strength that makes the
  inhibitor alone, at intact occupancy, account for the whole fold RHO;
* the MEAN form, `combine2.fit`, uses s = (RHO - 1) / (H * (q + RHO - 1)) with q the summed drop, which
  depends on the gene's activators as well and is therefore NOT the quantity the 39 were refused on.

The 39 were refused by `saturate.fit` on the SCALED observations, so the quantity that refused them is
IMPLIED_STRENGTH on `scalefit.scaled_rho(RHO)`, and this file's census is about that quantity alone. Two
thresholds follow from the formula by algebra, not from any measurement, and are committed here before any
number is read (`derived_thresholds`): with H = 0.5 the refusal is RHO > 2 on an unscaled observation --
removal more than DOUBLING the gene -- and, because scaling an inhibitor moves its fold towards 1 as
1 + SCALE * (RHO - 1), the refusal at M2h's scale is RHO > 1 + 1/SCALE = 5.66, removal multiplying the
gene by more than 5.66. Scaling can therefore only ever RELIEVE this refusal, which is why a gene refused
at the scale is necessarily refused unscaled too, and why the 39 are a subset of a larger unscaled set.

WHAT IS AT STAKE, and it is not the genes
-----------------------------------------
A fitted inhibitor strength above 1 is one of three things and the census is meant to tell them apart:
(i) a rule whose stated removal response really is that large, i.e. the compiled layer asserting a
repressor whose removal multiplies its gene several-fold; (ii) a SIGN or UNIT convention read the wrong
way round somewhere between the predictor's output and `bridge.observation_of`; or (iii) a threshold
artefact -- the bound s <= 1 cutting through the body of a continuous distribution of inhibitor folds, in
which case the 39 are not special, many genes outside them carry the same thing, and the refusal is a
property of the BOUND rather than of those 39 genes. QUESTIONS Q6 and Q7 are what decide (iii), and
BRANCHES is the wording each outcome gets, committed before the count.

This file is the RECORD: the questions, the definitions, the derived thresholds and the branch rule, all
committed before `scripts/inhibitor_census_run.py` computes any of them.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

from genomeos.attribution import bridge, combine2, residuemult, saturate, scalefit

LANE = "lane-inhibitor"

#: M2h's scale, taken through `scalefit` (which took it from `ac5580b`) and NOT re-declared, so no drift.
SCALE = scalefit.SCALE
SCALE_ROUNDED = scalefit.SCALE_ROUNDED

#: The intact Hill occupancy the bridge maps a removal fold through. Taken from `bridge`, not re-declared.
H = bridge.H_INTACT

#: The slack `combine2._strength_ok` allows above 1 before it refuses. Taken from `bridge`.
STRENGTH_TOLERANCE = bridge.STRENGTH_TOLERANCE

#: The refusal this lane is about, as `saturate.fit` and `combine2.fit` both spell it.
REASON = "inhibitor_strength_above_1"

#: The committed result the 39 are re-derived against (never copied from).
RESIDUEMULT_LOOKUP = "data/results/residuemult_lookup.json"

#: The one identity the tables share, the same convention `lane-residuemult` fixed in its IDENTITY.
MATCH_FIELD = "symbol"

#: A fold is "near" the threshold when it is within this factor of it, on the strength axis (Q7).
NEAR_BAND = 1.2

#: The quantiles every distribution in this census is reported at, chosen before any of them is seen.
QUANTILES = (0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0)

QUESTION = (
    "Why do 39 of the 3,708 residue genes carry an inhibitor of fitted strength above 1 -- the one reason"
    " `saturate.fit` refuses them at M2h's scale -- and what does that imply about the COMPILED LAYER"
    " rather than about those 39 genes? Is it a rule stating a removal response that large, a sign or"
    " unit convention read the wrong way round, or the bound s <= 1 cutting through the body of a"
    " continuous distribution of inhibitor folds?"
)

#: THE CENSUS QUESTIONS, committed before any of them is computed, so that no summary statistic in this
#: lane is chosen after its value is known. Each is answered once, as a distribution and not as a median.
QUESTIONS = {
    "Q1_re_derivation": (
        "Q1. RE-DERIVATION. Running `lane-residuemult`'s own code path again from the 24 pinned compiled"
        " programs -- `compile.gate_program` at `combine2.SELECTION`, `combine2.multi_rule_genes`,"
        " `combine2.gene_rules`, `combine2.classify`, the mean fit before and after"
        " `scalefit.scaled_rho`, then `saturate.fit` on both -- how many residue genes are informative"
        " under NEITHER form, and how many of those are refused for"
        " `inhibitor_strength_above_1`? The count is taken from THIS run; `residuemult_lookup.json`'s"
        " committed per-gene table is then read on the `symbol` column alone and RECONCILED against it,"
        " set against set. A disagreement is a finding of this lane and is reported, never repaired."
    ),
    "Q2_strength_distribution": (
        "Q2. THE DISTRIBUTION OF INHIBITOR STRENGTH among the 39. For each gene, the MAXIMUM over its"
        " inhibiting elements of IMPLIED_STRENGTH at the scale (the quantity that refused it), the same"
        " maximum unscaled, and the maximum raw removal fold RHO. Reported at QUANTILES with min, max and"
        " mean for each, over the 39 and never as a median alone. Also the SAME three distributions over"
        " every inhibiting element of the 39, not only each gene's worst, so that one gene with many"
        " inhibitors cannot stand for the set."
    ),
    "Q3_how_many_inhibitors": (
        "Q3. HOW MANY INHIBITORS EACH OF THE 39 CARRIES: the distribution of `n_inhibiting` (elements with"
        " RHO > 1) and of `n_rules` over the 39, at QUANTILES with min and max, and the same two"
        " distributions over the whole residue and over every eligible multi-element gene, so that the"
        " 39's counts can be read against the population's instead of in isolation."
    ),
    "Q4_one_inhibitor_or_several": (
        "Q4. ONE INHIBITOR OR SEVERAL. For each of the 39, how many of its inhibiting elements are"
        " themselves above the bound at the scale (IMPLIED_STRENGTH > 1 + STRENGTH_TOLERANCE). The"
        " distribution of that count, how many genes have exactly one offender, and -- for the genes with"
        " exactly one -- the share of the gene's inhibiting elements that the offender is. A single"
        " offending element on a gene with many inhibitors is a different object from a gene whose"
        " inhibitors are all that large, and the census must not merge them."
    ),
    "Q5_where_the_rules_came_from": (
        "Q5. WHERE THE OFFENDING RULES CAME FROM. For every offending element of the 39, the compiled"
        " rule's evidence KIND (`bridge._kind`: predicted / experimental / curated / inferred / none),"
        " the evidence SOURCE string as the program states it, the `how` string"
        " `bridge.observation_of` returns (which says whether RHO was read from an evidence note's log2"
        " fold change, from the rule strength below the clip, or from a CRISPRi EffectSize), the compiled"
        " cell label, and the chromosome. Reported as counts per category, plus the same counts over all"
        " inhibiting elements of all eligible genes as the background to read them against."
    ),
    "Q6_outside_the_39": (
        "Q6. HOW MANY GENES OUTSIDE THE 39 CARRY AN INHIBITOR ABOVE 1. The question that decides whether"
        " the 39 are special at all, asked at both thresholds and over three nested populations: (a) the"
        " 3,708 residue; (b) the 12,500 the mean form excludes; (c) every eligible multi-element gene."
        " In each, how many genes carry at least one inhibiting element with IMPLIED_STRENGTH above the"
        " bound AT THE SCALE, and how many unscaled, with the 39 removed from each count. If many genes"
        " outside the 39 carry one, then the 39 are not a distinguished set and this refusal is a general"
        " property of the compiled layer -- which would be the more important answer and is why this"
        " question is in the record before the count."
    ),
    "Q7_threshold_or_tail": (
        "Q7. IS THE THRESHOLD AT A GAP OR IN THE BODY? Over every inhibiting element of every eligible"
        " gene, the distribution of IMPLIED_STRENGTH at the scale and unscaled at QUANTILES, and the"
        " NEAR-BAND counts: how many elements sit just below the bound (strength in [1/NEAR_BAND, 1]) and"
        " how many above it. If the just-below count is at least the above count, the bound cuts through"
        " the body of a continuous distribution and the threshold is doing the work; if the offenders are"
        " separated from it, they are a tail. BRANCHES fixes that wording before the counts are known."
    ),
}

IMPLIED_STRENGTH_DEFINITION = (
    "IMPLIED_STRENGTH(RHO) = (1 - 1/RHO) / H, `saturate.fit`'s own expression for an inhibitor's strength,"
    " with H = bridge.H_INTACT = 0.5 taken from the bridge and not re-declared. It is defined only for an"
    " inhibitor (RHO > 1), where it is the strength that makes that element alone, at intact occupancy,"
    " account for the whole observed fold. `combine2.fit`'s mean-form inhibitor strength is a DIFFERENT"
    " expression -- (RHO - 1) / (H * (q + RHO - 1)), which also depends on the gene's activators -- and it"
    " is reported beside this one for the 39 but it is not the quantity they were refused on."
)

BRANCHES = {
    "threshold_is_doing_the_work": (
        "THE THRESHOLD IS DOING THE WORK, NOT THE BIOLOGY: the bound s <= 1 cuts through the body of a"
        " continuous distribution of compiled inhibitor folds -- {near_below} inhibiting elements sit"
        " within a factor of {near_band} just below it against {above} above it, and {outside_genes}"
        " eligible genes outside the 39 carry an inhibitor above the bound unscaled. The 39 are not a"
        " distinguished set of genes: they are the {residue} residue genes whose worst inhibitor happens"
        " to fall on the far side of a bound the compiled layer's folds cross continuously."
    ),
    "a_separated_tail": (
        "A SEPARATED TAIL: the offending folds are not continuous with the bound -- {near_below}"
        " inhibiting elements sit within a factor of {near_band} just below it against {above} above it --"
        " so the 39 are refused on observations that stand apart from the population's inhibitor folds"
        " rather than on where a bound was drawn. {outside_genes} eligible genes outside the 39 carry an"
        " inhibitor above the bound unscaled."
    ),
}

IMPLICATION_RULE = (
    "THE BRANCH IS DECIDED BY THE NEAR-BAND COUNTS ALONE, over every inhibiting element of every eligible"
    " multi-element gene: `threshold_is_doing_the_work` when the number of elements with IMPLIED_STRENGTH"
    " in [1/NEAR_BAND, 1] at the scale is at least the number above 1 + STRENGTH_TOLERANCE, and"
    " `a_separated_tail` otherwise. NEAR_BAND = 1.2 and the comparison are both fixed here, before the"
    " count. The number of genes outside the 39 is reported in either wording and does not choose it."
)

BOUND_CHECKS = (
    "Checks whose violation is reported and never repaired. (a) MONOTONICITY: scaling moves an inhibitor's"
    " fold towards 1, so IMPLIED_STRENGTH at the scale is at most IMPLIED_STRENGTH unscaled for every"
    " inhibiting element -- a violation would mean `scalefit.scaled_rho` does not do what its docstring"
    " says. (b) THE REFUSAL IS THE OFFENDER: every gene this lane finds refused by `saturate.fit` at the"
    " scale for `inhibitor_strength_above_1` has at least one inhibiting element with IMPLIED_STRENGTH"
    " above 1 + STRENGTH_TOLERANCE, and every gene with such an element is refused -- the census quantity"
    " and the code's own refusal must agree gene for gene, or the census is measuring something else."
    " (c) THE DERIVED THRESHOLDS: every offending element unscaled has RHO > 2 and every offending element"
    " at the scale has RHO > 1 + 1/SCALE, the two bounds `derived_thresholds` states by algebra."
    " (d) SUBSET: the genes refused at the scale are a subset of those refused unscaled."
    " (e) SIGN: every inhibiting element's compiled action is `inhibits`, which `combine2.classify`"
    " already enforces, re-asserted here so that a sign convention cannot pass unnoticed."
)

WHAT_THIS_IS_NOT = (
    "Not a test, not a registration, not a fit and not a claim about expression. No RNA value is opened."
    " THIS LANE PROPOSES NO CHANGE TO ANY FIT, no new form, no new threshold and no repair: it states what"
    " the compiled layer's inhibitor folds are and where the existing bound falls among them, and stops"
    " there. `lane-multrna`'s outcome on measured RNA was WORSE and is carried forward unchanged: the"
    " multiplicative form being the one that CAN RUN these genes is not the multiplicative form reading"
    " closer to measured RNA on them."
)

LIMITATIONS = (
    "A CENSUS of the compiled layer's own numbers. Every value is arithmetic on stated removal responses;"
    " nothing here is evidence that any of those responses is correct.",
    "The compiled rules are overwhelmingly PREDICTED (`lane-residuemult` recorded 43,681 predicted against"
    " 19 experimental rules gated on chr1 alone). A census of predicted folds says what the predictor"
    " asserts, not what a cell does.",
    "IMPLIED_STRENGTH assumes the bridge's intact occupancy H = 0.5 and the single-mechanism mapping. The"
    " strength of an inhibitor is not observable; it is this mapping's image of an observed fold.",
    "The 39 are a predictor-side refusal count. Whether those genes' RNA is well or badly predicted is not"
    " asked and cannot be read off anything here.",
    "Q6's populations are nested and the counts are of GENES, not elements; a gene with several offending"
    " elements counts once. Both are reported so neither can be mistaken for the other.",
    "The NEAR_BAND of 1.2 is a choice, fixed before the count. A different band would move the just-below"
    " count; the quantiles of the full distribution are reported so the band cannot hide the shape.",
    "Nothing in this lane distinguishes a correct large fold from a sign or unit error upstream. It can"
    " only say which rules the large folds come from, by evidence kind and `how` string.",
    "The scale is M2h's multiplicative-form scale, taken as a constant and not re-estimated; M2h's own"
    " heterogeneity finding (`40b0a8d`: not consistent with a single global scale) applies to it here too.",
)

CONVENTIONS = {
    "match_field": MATCH_FIELD,
    "scale": SCALE,
    "scale_is_refitted": False,
    "h_intact": H,
    "strength_tolerance": STRENGTH_TOLERANCE,
    "near_band": NEAR_BAND,
    "quantiles": list(QUANTILES),
    "form_refusing_the_39": "`saturate.fit` on `scalefit.scaled_rho(rho)` (the multiplicative form)",
    "informative_log2": residuemult.INFORMATIVE_LOG2,
    "rna_read": "none",
}


def derived_thresholds() -> dict[str, float]:
    """The two folds the bound s <= 1 becomes, by algebra on IMPLIED_STRENGTH. Committed before any count.

    Unscaled: (1 - 1/RHO)/H > 1 <=> RHO > 1/(1 - H) = 2 at H = 0.5. At the scale, an inhibitor's fold
    becomes 1 + SCALE*(RHO - 1), so the same bound is RHO > 1 + 1/(SCALE * (1/H - 1)) = 1 + 1/SCALE.
    """
    unscaled = 1.0 / (1.0 - H)
    return {
        "fold_above_which_an_inhibitor_is_refused_unscaled": unscaled,
        "fold_above_which_an_inhibitor_is_refused_at_the_scale": 1.0 + (unscaled - 1.0) / SCALE,
        "h_intact": H,
        "scale": SCALE,
    }


def implied_strength(rho: float) -> float:
    """`saturate.fit`'s own inhibitor strength for one removal fold (IMPLIED_STRENGTH_DEFINITION)."""
    if not rho > 1.0:
        return math.nan
    return (1.0 - 1.0 / rho) / H


def mean_form_strength(rho: float, q: float) -> float:
    """`combine2.fit`'s mean-form inhibitor strength, reported beside the other and never confused with it."""
    if not rho > 1.0:
        return math.nan
    x = rho - 1.0
    return x / (H * (q + x))


def above_bound(s: float) -> bool:
    """What `combine2._strength_ok` refuses, as a predicate: strength past 1 and past its tolerance."""
    return bool(s > 1.0 + STRENGTH_TOLERANCE)


@dataclass
class ElementRow:
    """One inhibiting element of one gene, as the compiled program states it."""

    symbol: str
    chrom: str
    element: str
    rho: float
    rho_scaled: float
    strength: float
    strength_scaled: float
    mean_strength: float
    offends_unscaled: bool
    offends_scaled: bool
    kind: str
    how: str
    source: str
    compiled_cell: str

    def as_row(self) -> list[Any]:
        return [
            self.symbol,
            self.chrom,
            self.element,
            self.rho,
            self.rho_scaled,
            self.strength,
            self.strength_scaled,
            self.mean_strength,
            int(self.offends_unscaled),
            int(self.offends_scaled),
            self.kind,
            self.how,
            self.source,
            self.compiled_cell,
        ]


ELEMENT_COLUMNS = (
    "symbol",
    "chrom",
    "element",
    "rho",
    "rho_scaled",
    "strength",
    "strength_scaled",
    "mean_form_strength",
    "offends_unscaled",
    "offends_scaled",
    "evidence_kind",
    "how_rho_was_read",
    "evidence_source",
    "compiled_cell",
)


@dataclass
class GeneRow:
    """One gene's inhibitor census: the counts and maxima Q2, Q3, Q4 and Q6 are built from."""

    symbol: str
    chrom: str
    n_rules: int
    n_inhibiting: int
    n_offending_unscaled: int
    n_offending_scaled: int
    max_rho: float
    max_strength: float
    max_strength_scaled: float
    refused_unscaled: bool
    refused_scaled: bool
    reason_unscaled: str
    reason_scaled: str
    in_residue: bool
    excluded_by_the_mean: bool
    elements: list[ElementRow] = field(default_factory=list)

    def as_row(self) -> list[Any]:
        return [
            self.symbol,
            self.chrom,
            self.n_rules,
            self.n_inhibiting,
            self.n_offending_unscaled,
            self.n_offending_scaled,
            self.max_rho,
            self.max_strength,
            self.max_strength_scaled,
            int(self.refused_unscaled),
            int(self.refused_scaled),
            self.reason_unscaled,
            self.reason_scaled,
            int(self.in_residue),
            int(self.excluded_by_the_mean),
        ]


GENE_COLUMNS = (
    "symbol",
    "chrom",
    "n_rules",
    "n_inhibiting",
    "n_offending_unscaled",
    "n_offending_scaled",
    "max_rho",
    "max_strength",
    "max_strength_scaled",
    "refused_unscaled",
    "refused_scaled",
    "reason_unscaled",
    "reason_scaled",
    "in_residue",
    "excluded_by_the_mean",
)


def element_rows(
    gr: combine2.GeneRules, kinds: dict[str, tuple[str, str, str]], q: float
) -> list[ElementRow]:
    """Every inhibiting element of one gene, with both strengths and the rule's provenance.

    `kinds` is element -> (evidence kind, the `how` string of `bridge.observation_of`, evidence source),
    read by the caller from the gated module's own rules. `q` is the mean form's summed drop, used only
    for the mean-form strength reported beside the one that refused.
    """
    out = []
    scaled = scalefit.scaled_rho(gr.rho, SCALE)
    for el in gr.elements:
        r = gr.rho[el]
        if not r > 1.0:
            continue
        kind, how, source = kinds.get(el, ("", "", ""))
        s, s_scaled = implied_strength(r), implied_strength(scaled[el])
        out.append(
            ElementRow(
                symbol=gr.symbol,
                chrom=gr.chrom,
                element=el,
                rho=r,
                rho_scaled=scaled[el],
                strength=s,
                strength_scaled=s_scaled,
                mean_strength=mean_form_strength(r, q),
                offends_unscaled=above_bound(s),
                offends_scaled=above_bound(s_scaled),
                kind=kind,
                how=how,
                source=source,
                compiled_cell=gr.compiled.get(el, ""),
            )
        )
    return out


def gene_row(
    gr: combine2.GeneRules,
    rows: list[ElementRow],
    in_residue: bool,
    excluded_by_the_mean: bool,
) -> GeneRow:
    """One gene's census row. NOTHING IS FITTED: `saturate.fit` has no parameter and SCALE is a constant."""
    unscaled = saturate.fit(gr.rho)
    at_scale = saturate.fit(scalefit.scaled_rho(gr.rho, SCALE))
    return GeneRow(
        symbol=gr.symbol,
        chrom=gr.chrom,
        n_rules=len(gr.elements),
        n_inhibiting=sum(1 for r in gr.rho.values() if r > 1.0),
        n_offending_unscaled=sum(1 for e in rows if e.offends_unscaled),
        n_offending_scaled=sum(1 for e in rows if e.offends_scaled),
        max_rho=max((e.rho for e in rows), default=float("nan")),
        max_strength=max((e.strength for e in rows), default=float("nan")),
        max_strength_scaled=max((e.strength_scaled for e in rows), default=float("nan")),
        refused_unscaled=not unscaled.ok,
        refused_scaled=not at_scale.ok,
        reason_unscaled=unscaled.reason,
        reason_scaled=at_scale.reason,
        in_residue=in_residue,
        excluded_by_the_mean=excluded_by_the_mean,
    )


def gene_bound_violations(row: GeneRow) -> list[str]:
    """BOUND_CHECKS (a) to (e) on one gene's rows. Reported, never repaired."""
    bad = []
    t = derived_thresholds()
    for e in row.elements:
        if e.strength_scaled > e.strength + 1e-12:
            bad.append(f"{row.symbol}/{e.element}: strength at the scale above the unscaled one (a)")
        if e.rho_scaled > e.rho + 1e-12:
            bad.append(f"{row.symbol}/{e.element}: scaled fold above the raw fold (a)")
        if e.offends_unscaled and not e.rho > t["fold_above_which_an_inhibitor_is_refused_unscaled"]:
            bad.append(f"{row.symbol}/{e.element}: offends unscaled with RHO {e.rho} (c)")
        if e.offends_scaled and not e.rho > t["fold_above_which_an_inhibitor_is_refused_at_the_scale"]:
            bad.append(f"{row.symbol}/{e.element}: offends at the scale with RHO {e.rho} (c)")
        if e.offends_scaled and not e.offends_unscaled:
            bad.append(f"{row.symbol}/{e.element}: offends at the scale but not unscaled (d)")
    refused = row.reason_scaled == REASON
    if refused != (row.n_offending_scaled > 0):
        bad.append(
            f"{row.symbol}: `saturate.fit` at the scale reads {row.reason_scaled!r} with"
            f" {row.n_offending_scaled} offending element(s) -- the census quantity and the code's own"
            f" refusal disagree (b)"
        )
    if row.refused_scaled and not row.refused_unscaled:
        bad.append(f"{row.symbol}: refused at the scale but runnable unscaled, against the bound (d)")
    return bad


def distribution(xs: list[float]) -> dict[str, Any]:
    """One distribution as this census reports it: n, min, max, mean and QUANTILES. Never a median alone."""
    if not xs:
        return {"n": 0}
    s = sorted(xs)
    out: dict[str, Any] = {"n": len(s), "min": s[0], "max": s[-1], "mean": sum(s) / len(s)}
    out["quantiles"] = {f"{q:g}": s[min(len(s) - 1, int(q * (len(s) - 1) + 0.5))] for q in QUANTILES}
    return out


def counts_of(xs: list[int]) -> dict[str, int]:
    """The exact value counts of a small integer distribution, as well as its quantiles (Q3, Q4)."""
    return {str(v): sum(1 for x in xs if x == v) for v in sorted(set(xs))}


def branch(near_below: int, above: int) -> str:
    """IMPLICATION_RULE, applied. The near-band counts alone decide; nothing else may."""
    return "threshold_is_doing_the_work" if near_below >= above else "a_separated_tail"


def reading(which: str, **numbers: Any) -> str:
    """BRANCHES, filled. The wording is the committed one; only the numbers are substituted."""
    return BRANCHES[which].format(**numbers)


#: The finished ROADMAP row for this lane, written AFTER the run and carrying the registered wording of
#: `BRANCHES` word for word. docs/ROADMAP.md has one editor (genomeos-f7); this is the row text to insert.
ROADMAP_ROW = (
    "**Area I, why 39 genes carry an inhibitor above 1 (lane-inhibitor, 2026-10-12; a4c79ce, 2816a3e,"
    " 85ca300, a1a9a90).** Re-derived from the 24 compiled programs, matched on `symbol`: the 39 reconcile"
    " exactly with lane-residuemult's set, all on `inhibitor_strength_above_1`. 38 carry ONE offending"
    " element of a median 15 inhibitors; all 40 offenders are PREDICTED deletion notes, log2 fold 2.6"
    " to 6.6. Worst strength at the scale 1.037 to 1.907 (median 1.178), ceiling 2. 149 genes OUTSIDE"
    " the 39"
    " carry one at the scale and 1,199 unscaled: the refusal is general to the compiled layer, not to"
    " these genes. Branch `a_separated_tail`. 0 violations, no RNA, nothing fitted, no money."
)
