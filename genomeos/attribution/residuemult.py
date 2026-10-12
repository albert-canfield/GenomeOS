# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-residuemult (ROADMAP §5, area I, 2026-10-12): are the 3,708 the mean cannot fit already fitted
multiplicatively? A LOOKUP in two committed tables, not a new form.

STOP-LIST CHECK (docs/ROADMAP.md §5, "The stop list of 2026-09-27, reinstated"): no AlphaGenome request,
no deletion request, no biosample added, nothing in area J, no node-model work, no BioForge calibration.
This lane reads NO RNA VALUE, fits NOTHING, registers NO TEST, and makes no claim about expression.

WHAT IS ALREADY IN THE RECORD, and why this is the right question
-----------------------------------------------------------------
`lane-scalefit` (rule and floor `8469df9`, census `7f65544`) and `lane-residue` (`5d2c8f5`) measured, with
no RNA opened, that the runtime's MEAN form excludes 12,500 of the 17,983 multi-element genes
`combine2.classify` admits as `basal_not_positive` -- their single-removal drops sum to or past the gene's
whole intact level, so no basal rate is left -- that 8,792 of them get a positive basal rate back at M2h's
fixed scale 0.2146, and that 3,708 do NOT. `lane-residue` then measured that the residue's excess is
DIFFUSE: the median residue gene's largest activating drop is 8.88% of its summed drop, no gene in the
population has one above 21.0%, n_eff is 25.9 of about 37 activating elements, and a PER-ELEMENT CAP
recovers 0 genes at every level L >= 0.25 of the declared ladder.

A saturating whole-set form was then proposed and refused as ALREADY DISPROVED IN THE RECORD, by
`genomeos/attribution/saturate.py`'s `WHY_NOT_SATURATING`, which `lane-saturate` (7f) wrote before any fit:
under ANY concave, nondecreasing composition of the summed element input with f(0) = 0 and a nonnegative
basal rate, the single-removal drops sum to AT MOST the gene's whole level -- "saturation makes them sum to
LESS than addition does (redundancy), not more". The residue's drops sum to 1.42 of the level even after
the scale. So NO saturating and no sub-additive whole-set form can fit a single one of the 3,708, and a
census of one would count to zero by theorem, not by measurement.

Drops that sum PAST the whole need a composition that is CONVEX where the gene sits -- multiplicative, i.e.
log-additive, which is synergy -- and the engine already has one: `saturate.fit`'s multiplicative form, the
form `lane-multrna` (7h) ran across the thirteen reader cells. The same file records 7f fitting 199 of 210
genes multiplicatively where concave forms could fit at most 12.

THE QUESTION OF THIS LANE, and it is only a lookup
--------------------------------------------------
How many of the 3,708 are FITTED (runnable and informative) in 7h's COMMITTED multiplicative table, and how
many under that same form with M2h's fixed scale applied? Both counts are matched identities and a
runnability test that already exists; nothing here is fitted, tuned, selected or re-estimated, M2h's scale
is taken verbatim from `scalefit.SCALE` (which copied it from `ac5580b`), and the multiplicative fit and
its runnability test are taken verbatim from `saturate.fit` and `multrna.predicted_log2`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from genomeos.attribution import combine2, multrna, residue, saturate, scalefit

LANE = "lane-residuemult"

#: M2h's scale, taken from `scalefit` (which took it from `ac5580b`) and NOT re-declared, so no drift.
SCALE = scalefit.SCALE
SCALE_ROUNDED = scalefit.SCALE_ROUNDED

#: The informativeness threshold, taken from 7d/7h and not re-declared.
INFORMATIVE_LOG2 = multrna.INFORMATIVE_LOG2

#: The committed tables this lane looks up. Nothing in them is recomputed except where RE_DERIVED says.
COMMITTED = {
    "multrna_census": "data/results/multrna_census.json",
    "multrna_result": "data/results/multrna_result.json",
    "residue_census": "data/results/residue_census.json",
    "residue_genes": "data/results/residue_genes.json",
    "magcal_result": "data/results/magcal_result.json",
}

QUESTION = (
    "How many of the 3,708 multi-element genes the runtime's MEAN form cannot fit even at M2h's fixed scale"
    " 0.2146 -- `lane-residue`'s residue -- are FITTED (runnable and informative) under the MULTIPLICATIVE"
    " form of `lane-multrna`'s committed table, and how many under that same form with M2h's scale applied?"
    " A lookup: two identity matches and a runnability test that already exist. No RNA value is opened,"
    " nothing is fitted, and no form is proposed, because the form in question is already validated."
)


# ==== the identity convention, stated before the match ===================================================

IDENTITY = (
    "The MATCH FIELD IS THE GENE SYMBOL, `combine2.GeneRules.symbol`, which is the `symbol` column of"
    " `lane-multrna`'s committed per-gene table (`data/results/multrna_result.json`, field"
    " `genes_table.columns`) and the same string `lane-residue` passed to `residue.gene_residue` as the"
    " gene's name. It is the only identity the two tables share: 7h's table also carries `gene` and"
    " `gene_id` (ENSG), but `lane-residue`'s committed per-gene table"
    " (`data/results/residue_genes.json`) carries NEITHER -- its rows hold only `chrom` and the"
    " distribution columns, with no gene identity at all. THAT IS WHY THE RESIDUE'S IDENTITIES ARE"
    " RE-DERIVED here rather than read (see RE_DERIVED), and why the re-derivation is reconciled against"
    " `lane-residue`'s committed counts as a check of this lane."
)

RE_DERIVED = (
    "The residue's 3,708 IDENTITIES are re-derived exactly as `lane-residue` derived its counts, by the"
    " same code on the same pinned programs: the 24 compiled non-coding programs are gated and parsed with"
    " `compile.gate_program` at `combine2.SELECTION`, the multi-element genes come from"
    " `combine2.multi_rule_genes` and `combine2.gene_rules`, a gene is ELIGIBLE when `combine2.classify`"
    " gives no reason, EXCLUDED when `combine2.fit('mean', rho)` fails with `basal_not_positive`, and in"
    " the RESIDUE when `combine2.fit('mean', scalefit.scaled_rho(rho))` fails the same way -- two calls to"
    " 7d's own fit per gene, never the identity q_scaled = SCALE * q. `lane-residue`'s committed"
    " Q1_counts are read AFTERWARDS and reconciled; a disagreement is a finding of this lane and not"
    " something to repair. The only thing this lane adds to that run is the gene's SYMBOL, which"
    " `lane-residue` computed and did not serialise."
)


# ==== the two counts, defined before they are computed ===================================================

COUNT_COMMITTED = (
    "COUNT 1, 'fitted in the committed table': the residue symbols that appear in"
    " `multrna_result.json`'s `genes_table` rows. That table IS 7h's committed multiplicative population:"
    " the genes the multiplicative form made runnable and informative AND that passed 7h's RNA-table row"
    " checks, which is 16,509 of the 16,678 its census counts informative (169 fell to"
    " `population_exclusions.multiplicative`: 139 `gene_id_not_on_exactly_one_row_of_every_table` and 30"
    " `symbol_not_on_exactly_one_gene_id`). ONLY THE `symbol` COLUMN OF THAT TABLE IS READ. Its `P` and"
    " `D` columns are predictions and RNA and are NOT opened, loaded or touched by this lane. Reported"
    " beside COUNT 1b, the same question asked PREDICTOR-SIDE ONLY -- `saturate.fit(rho).ok` and the"
    " closed form's range above INFORMATIVE_LOG2 across the thirteen reader cells, which is 7h's census"
    " definition of informative before any RNA table -- so that the row checks' 169 cannot be mistaken"
    " for a failure of the form."
)

COUNT_SCALED = (
    "COUNT 2, 'fitted at 0.2146': the residue genes for which `saturate.fit(scalefit.scaled_rho(rho))` is"
    " ok and the multiplicative closed form's predicted log2 level is not the same in all thirteen reader"
    " cells (its range exceeds INFORMATIVE_LOG2) -- THE SAME runnability and informativeness test 7h's"
    " committed table used, applied to the scaled observations. `scaled_rho` is `scalefit`'s own:"
    " RHO_i -> 1 - SCALE * (1 - RHO_i), every element's effect multiplied by M2h's fixed scale, taken as a"
    " CONSTANT from `ac5580b` via `scalefit.SCALE`. NOTHING IS FITTED: the scale is not re-estimated, no"
    " level is chosen after seeing a count, and the multiplicative form has no parameter of its own"
    " (`saturate.FIT['nothing_fitted']`). The scale is M2h's MULTIPLICATIVE-form scale"
    " (`magcal_result.json`, `scales.multiplicative`), which is the form it is applied to here."
)

COMPLEMENT = (
    "THE COMPLEMENT: how many of the 3,708 are informative under NEITHER -- neither in 7h's committed"
    " table nor under the scaled multiplicative test -- reported with each gene's reason from"
    " `saturate.fit`, so that a count of genes the multiplicative form cannot run is never reported as a"
    " count of genes it fits. A gene is in the complement when it is in neither count."
)

WHY_THE_SCALE_CAN_ONLY_HELP_THE_MULTIPLICATIVE_FORM = (
    "Stated in advance, from the algebra, so that neither count can be read as a surprise. The"
    " multiplicative form has no basal share to exhaust -- with every element removed the level is"
    " T*exp(-sum beta)/R0 > 0 -- so `basal_not_positive`, the reason the mean excludes all 12,500, CANNOT"
    " arise under it at all. Its only predictor-side refusals are `inhibitor_strength_above_1` (an"
    " inhibiting element with RHO_j > 1/(1 - H) = 2) and `exponent_beyond_float_range` (sum of -ln RHO_i"
    " above `saturate.EXPONENT_LIMIT` = 700). Multiplying every effect by SCALE = 0.2146 shrinks an"
    " inhibitor's lift, so the inhibitor bound becomes RHO_j <= 1 + 1/SCALE = 5.660, and shrinks every"
    " beta_i by more than the log does, so the exponent bound gets no harder. THEREFORE COUNT 2 CANNOT BE"
    " BELOW THE PREDICTOR-SIDE COUNT 1b, and the only genes either count can lose are the ones with a"
    " strong inhibitor. That is an arithmetic prediction, declared here before the count, and the run"
    " CHECKS it: `BOUND_CHECKS` below. ONE CAVEAT, stated rather than discovered: RUNNABILITY cannot be"
    " lost to the scale, but INFORMATIVENESS could be, because scaling shrinks every |log2 RHO_i| and so"
    " shrinks the closed form's range across the thirteen reader cells; a gene whose unscaled range was"
    " already within a factor 1/SCALE of INFORMATIVE_LOG2 = 1e-9 could fall below it. So the claim is"
    " exact for runnability and is a near-certainty, not a theorem, for informativeness. Every gene on"
    " which it fails is REPORTED, by name, and nothing is repaired."
)

BOUND_CHECKS = (
    "Checked on every one of the 3,708, and a violation is reported, not repaired: (a) every residue gene"
    " refused by `saturate.fit` is refused for `inhibitor_strength_above_1` or"
    " `exponent_beyond_float_range` and NEVER for `basal_not_positive`, which the multiplicative form has"
    " no state for; (b) every gene runnable unscaled is runnable scaled, which is exact; (c) every gene"
    " informative unscaled is informative scaled, which is the near-certainty the caveat above names, so"
    " a failure is a named finding and not an error; (d) `combine2.fit('mean', rho)` fails"
    " `basal_not_positive` and `combine2.fit('mean', scaled_rho(rho))` fails it too, which is the residue's"
    " own definition, re-asserted per gene; (e) the mean form's q on the scaled observations equals SCALE"
    " times its q on the unscaled, the identity `lane-scalefit` proved and `lane-residue` re-checked."
)


# ==== the implication, with its threshold, declared BEFORE the count =====================================

#: "Most" is a majority of the residue. Fixed here before a single gene was matched.
MOST = 0.5

IMPLICATION_RULE = (
    "Decided on the UNION of COUNT 1 and COUNT 2 as a share of the 3,708, with the threshold MOST = 0.5"
    " fixed here before the match. ALREADY_ANSWERED, if the share exceeds 0.5: most of the residue is"
    " fitted multiplicatively, so THE RESIDUE IS THE ADDITIVE FORM'S STRUCTURAL FAILURE AND IS ALREADY"
    " ANSWERED BY THE FORM M2h VALIDATED, AND NO NEW FORM IS NEEDED. NOT_ANSWERED, if the share is 0.5 or"
    " below: most of the residue is not fitted multiplicatively either, which is reported and NOTHING IS"
    " PROPOSED. The threshold is the plain reading of the word: a majority is 'most'. Neither branch is a"
    " claim about expression: `multrna`'s own outcome on measured RNA was WORSE, which LIMITATIONS carries"
    " forward word for word, and a gene a fit can accept is not a gene the model predicts -- which"
    " `lane-recovered` (`f775e9b`) had to establish as a separate registration on separate data."
)

READINGS = {
    "ALREADY_ANSWERED": (
        "ALREADY_ANSWERED: {union} of the {residue} genes the runtime's mean form cannot fit even at M2h's"
        " fixed 0.2146 ({share}) are fitted under the multiplicative form -- {committed} of them are in"
        " `lane-multrna`'s committed table and {scaled} are runnable and informative under the same test"
        " with the scale applied -- so THE RESIDUE IS THE ADDITIVE FORM'S STRUCTURAL FAILURE AND IS"
        " ALREADY ANSWERED BY THE FORM M2h VALIDATED, AND NO NEW FORM IS NEEDED. {complement} are"
        " informative under neither. No RNA value was opened and no test is registered here"
    ),
    "NOT_ANSWERED": (
        "NOT_ANSWERED: only {union} of the {residue} genes the runtime's mean form cannot fit even at"
        " M2h's fixed 0.2146 ({share}) are fitted under the multiplicative form -- {committed} in"
        " `lane-multrna`'s committed table and {scaled} under the same test with the scale applied --"
        " {complement} are informative under neither, and THIS LANE PROPOSES NOTHING. No RNA value was"
        " opened and no test is registered here"
    ),
}


def reading(branch: str, **numbers: Any) -> str:
    return READINGS[branch].format(**numbers)


def implication(union: int, residue_genes: int) -> str:
    """ALREADY_ANSWERED or NOT_ANSWERED by IMPLICATION_RULE, on the union share alone."""
    if residue_genes <= 0:
        raise ValueError("the residue is empty: there is nothing to decide")
    return "ALREADY_ANSWERED" if union / residue_genes > MOST else "NOT_ANSWERED"


# ==== the per-gene lookup ================================================================================


@dataclass
class ResidueGene:
    """One residue gene under the multiplicative form, by lookup and by the committed runnability test."""

    symbol: str
    chrom: str
    n_rules: int
    n_inhibiting: int
    in_committed_table: bool  # COUNT 1: its symbol is a row of multrna_result.json's genes_table
    runnable_unscaled: bool  # saturate.fit(rho).ok
    informative_unscaled: bool  # COUNT 1b: and the closed form varies across the thirteen reader cells
    runnable_scaled: bool  # saturate.fit(scaled_rho(rho)).ok
    informative_scaled: bool  # COUNT 2
    reason_unscaled: str
    reason_scaled: str
    range_log2_scaled: float

    @property
    def fitted(self) -> bool:
        """In COUNT 1 or COUNT 2: the union IMPLICATION_RULE decides on."""
        return self.in_committed_table or self.informative_scaled


def informative_under(fit: saturate.MFit, gr: combine2.GeneRules) -> tuple[bool, float]:
    """7h's own test: runnable, and the multiplicative closed form's log2 level is not the same in all
    thirteen reader cells. Returns (informative, the range in log2), the range being nan when not run."""
    if not fit.ok:
        return False, float("nan")
    s = multrna.predicted_log2(fit, gr)
    rng = max(s) - min(s)
    return rng > INFORMATIVE_LOG2, rng


def lookup_gene(gr: combine2.GeneRules, in_committed_table: bool) -> ResidueGene:
    """One residue gene's row. NOTHING IS FITTED: `saturate.fit` has no parameter and SCALE is a
    constant taken from `ac5580b` through `scalefit`."""
    unscaled = saturate.fit(gr.rho)
    scaled = saturate.fit(scalefit.scaled_rho(gr.rho, SCALE))
    inf_u, _ = informative_under(unscaled, gr)
    inf_s, rng_s = informative_under(scaled, gr)
    return ResidueGene(
        symbol=gr.symbol,
        chrom=gr.chrom,
        n_rules=len(gr.elements),
        n_inhibiting=sum(1 for r in gr.rho.values() if r > 1.0),
        in_committed_table=in_committed_table,
        runnable_unscaled=unscaled.ok,
        informative_unscaled=inf_u,
        runnable_scaled=scaled.ok,
        informative_scaled=inf_s,
        reason_unscaled=unscaled.reason,
        reason_scaled=scaled.reason,
        range_log2_scaled=rng_s,
    )


def bound_violations(row: ResidueGene) -> list[str]:
    """BOUND_CHECKS (a), (b), (c) on one gene. Reported, never repaired."""
    bad = []
    if not row.runnable_unscaled and row.reason_unscaled not in saturate_refusals():
        bad.append(f"{row.symbol}: multiplicative refusal {row.reason_unscaled!r} is not a declared one")
    if not row.runnable_scaled and row.reason_scaled not in saturate_refusals():
        bad.append(f"{row.symbol}: scaled refusal {row.reason_scaled!r} is not a declared one")
    if row.runnable_unscaled and not row.runnable_scaled:
        bad.append(f"{row.symbol}: runnable unscaled but not at the scale, against the stated bound")
    if row.informative_unscaled and not row.informative_scaled:
        bad.append(f"{row.symbol}: informative unscaled but not at the scale, against the stated bound")
    return bad


def saturate_refusals() -> frozenset[str]:
    """The only predictor-side reasons the multiplicative form can refuse a gene the mean EXCLUDED as
    `basal_not_positive` for: it has no basal share to exhaust, so that reason cannot arise under it."""
    return frozenset({"inhibitor_strength_above_1", "exponent_beyond_float_range"})


LIMITATIONS = (
    "A LOOKUP, not a test. Every number here is whether a fit accepts a gene, which is arithmetic on the"
    " elements' own removal observations. No RNA value was opened in this lane and no claim about"
    " expression is made or implied.",
    "`lane-multrna`'s outcome on measured RNA was WORSE, carried forward word for word: the multiplicative"
    " form being the one that CAN RUN these genes is not the multiplicative form reading closer to"
    " measured RNA on them. Whether any of the 3,708 would read closer under it is a SEPARATE"
    " registration on separate data, as `lane-recovered` (`f775e9b`) had to be for the genes that"
    " recovered at the scale alone.",
    "COUNT 1 is a membership test against a committed table, so it inherits that table's population"
    " exclusions: 169 of 7h's 16,678 informative genes are absent from it for RNA-table row reasons"
    " (139 `gene_id_not_on_exactly_one_row_of_every_table`, 30 `symbol_not_on_exactly_one_gene_id`) and"
    " not for anything about the form. COUNT 1b is the same question without those exclusions and is"
    " reported beside it.",
    "The match is on the gene SYMBOL, which 7h itself found is not always one gene_id: its census counts"
    " 16,648 of the 16,678 informative genes whose symbol names exactly one. A residue symbol that names"
    " two gene_ids could match a 7h row that is not the same locus. The count of residue symbols that are"
    " not unique within this run is reported.",
    "The residue identities are RE-DERIVED, not read, because `lane-residue`'s committed per-gene table"
    " carries no gene identity. The re-derivation is reconciled against its committed Q1_counts and a"
    " disagreement is reported as a finding.",
    "M2h's scale is applied to the multiplicative form's observations as a fixed constant. That is what"
    " `magcal` fitted it for (`scales.multiplicative`), but it was fitted on 8,181 TRAINING genes, and"
    " nothing here establishes that it is the right magnitude for THIS population -- `lane-scalehet` holds"
    " the question of whether one global scale is even the right shape.",
    "The multiplicative form keeps the bridge's exact-reproduction principle at the intact state"
    " (`saturate.FIT['principle']`); the SCALED form does not, exactly as `scalefit` recorded for the"
    " mean. A scaled element's run reproduces 0.2146 of its own removal observation, not the"
    " observation.",
    "This lane registers no RNA test, proposes no form and changes no default. The runtime's"
    " `activator_combination` default is untouched.",
)

PRIOR_EXPOSURE = (
    "What this lane had seen before it was written, stated so that no number here can be mistaken for a"
    " blind one. `lane-scalefit`'s and `lane-residue`'s counts (17,983 eligible, 12,500 excluded, 8,792"
    " recovered, 3,708 in the residue; median top1 0.0888, max top1 0.210, n_eff 25.9, 0 recovered by a"
    " cap at L >= 0.25) are in the brief and in the commit messages `7f65544` and `5d2c8f5`."
    " `lane-multrna`'s committed census sentence is in `multrna_census.json` and was read before this"
    " file was written: 16,784 of the 17,983 runnable multiplicatively, 16,678 informative, 11,802 of"
    " them informative under neither 7d form, of which 11,396 are `basal_not_positive` under the mean."
    " THAT LAST NUMBER ALREADY BOUNDS THIS LOOKUP from one side and was known: the 3,708 are a subset of"
    " the 12,500 the mean excludes as `basal_not_positive`, and 11,396 of those 12,500 are informative"
    " multiplicatively, so at most 1,104 of the 3,708 can fail to be. The counts below are therefore NOT"
    " a surprise and are not reported as one; what they add is the exact figure per population and the"
    " scaled form's, neither of which is in the record. `saturate.WHY_NOT_SATURATING` and its 199-of-210"
    " were read too. The IMPLICATION_RULE's threshold was fixed knowing all of this, which is why it is"
    " a plain majority and not a number chosen to be cleared."
)

CONVENTIONS = {
    "form_mean": "the runtime default, `combine2.fit('mean', ...)`, at THETA = 0.5",
    "form_multiplicative": "`saturate.fit`, the form `lane-multrna` ran, at bridge.H_INTACT = 0.5",
    "scale": SCALE,
    "scale_is_refitted": False,
    "informative_log2": INFORMATIVE_LOG2,
    "match_field": "symbol",
    "most": MOST,
    "residue_threshold_on_the_unscaled_ratio": residue.RESIDUE_THRESHOLD,
}

BOUND = {
    "question": "QUESTION",
    "identity": "IDENTITY",
    "re_derived": "RE_DERIVED",
    "count_committed": "COUNT_COMMITTED",
    "count_scaled": "COUNT_SCALED",
    "complement": "COMPLEMENT",
    "why_the_scale_can_only_help": "WHY_THE_SCALE_CAN_ONLY_HELP_THE_MULTIPLICATIVE_FORM",
    "bound_checks": "BOUND_CHECKS",
    "implication_rule": "IMPLICATION_RULE",
    "most": "MOST",
    "readings": "READINGS",
    "limitations": "LIMITATIONS",
    "prior_exposure": "PRIOR_EXPOSURE",
    "conventions": "CONVENTIONS",
    "committed": "COMMITTED",
}
