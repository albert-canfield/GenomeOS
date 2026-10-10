# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""M2e (ROADMAP §5 item 7f, lane-saturate, 2026-10-10): the bridge stage under a composition that fits.

M1c (`signaltrace.py`, 9637584 / 1103719) traced the 1,918 held-out ENCODE CRISPRi pairs in K562
through the compiled layer: target choice kept 39 pairs and a gain over activity and distance; the bridge
kept none (every reached gene carries two or more rules, `not_identifiable`). M2d (`s4joint.py`, a9055e0 /
86ddf73) replaced the bridge with combine2's joint fit and still lost the signal there, because the fit
could run few of the genes: of the 210 genes the traced pairs reach, 10 under the runtime's mean and 2
under promoter competition, the rest mostly `basal_not_positive` - a gene's single-removal drops add up
to more than its whole expression, which leaves no basal rate.

The brief for this item proposed a SATURATING composition as the explanation. It is not one, and
`WHY_NOT_SATURATING` says why in writing, before any fit: under a concave composition the single-removal
drops add up to LESS than the whole, never more. The composition registered here is the one that can
fit such genes and has a citation in K562 CRISPRi: the MULTIPLICATIVE (log-additive) composition of Zhou
et al. 2024, added to the runtime as the opt-in rule `multiplicative` (1f8f12b). It has no parameter of
its own; nothing is fitted on any CRISPRi pair.

This module re-implements no stage: the target and the lookups are signaltrace's, the gene reading and
the exclusions are combine2's (`gene_rules`, `classify`), the comparators and the gene records are
s4joint's. Its own code is the multiplicative fit, its runs, the join and the readings.

In the order written and committed: the runtime rule (1f8f12b), this code, the registration
(data/results/saturate_registration.json), then one run. No model request, no key, no money.
"""

from __future__ import annotations

import copy
import json
import math
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any

from genomeos.attribution import bridge, crispri
from genomeos.attribution import combine2 as c2
from genomeos.attribution import compile as cp
from genomeos.attribution import readergate as rgate
from genomeos.attribution import s4joint as sj
from genomeos.attribution import signaltrace as st
from genomeos.ir import Action, Module
from genomeos.runtime import grn

REGISTERED = "2026-10-10"
LANE = "lane-saturate"
CELL = st.CELL
K562_INDEX = sj.K562_INDEX
H = bridge.H_INTACT
FORM = "multiplicative"
#: exp(G * W) must stay a finite double: a gene whose summed log-folds exceed this is not run
EXPONENT_LIMIT = 700.0

# ---- the form, chosen in writing before any fit -----------------------------------------------------------
WHY_NOT_SATURATING = (
    "The brief proposed a saturating composition: level = b + f(W), W the summed element input and f"
    " saturating in it (Michaelis-Menten or Hill in W, or a logistic). A concave f cannot give what M2d"
    " found. Let f be concave and nondecreasing with f(0) = 0, b >= 0 and the inhibitors' product fixed."
    " Removing element i lowers the level by f(W) - f(W - w_i). Ordering the elements 1..k, f(W) - f(W -"
    " w_i) <= f(w_1 + .. + w_i) - f(w_1 + .. + w_(i-1)), because a concave function's rise over an interval"
    " of fixed length falls as the interval moves right; summed over i, the single-removal drops add up to"
    " at most f(W) - f(0), which is at most the whole level. So under any concave composition with a"
    " nonnegative basal rate the drops sum to at most the whole: saturation makes them sum to LESS than"
    " addition does (redundancy), not more. combine2's competition form is such a composition"
    " (Michaelis-Menten in W) and ran fewer of M2d's genes than the mean (2 against 10 of 210). Drops that"
    " sum past the whole need a composition that is convex where the gene sits: multiplicative"
    " (log-additive), or a sigmoid below its inflection, i.e. away from saturation. The form registered"
    " here is the multiplicative one; no saturating form is registered, and the census reports how many of"
    " the reached genes any concave composition could fit at most"
)
FORM_SOURCE = (
    "Zhou JL, Guruvayurappan K, Toneyan S, Chen HV, Chen AR, Koo P, McVicker G (2024) Analysis of"
    " single-cell CRISPR perturbations indicates that enhancers predominantly act multiplicatively. Cell"
    " Genomics 100672, doi 10.1016/j.xgen.2024.100672 (bioRxiv 2023.04.26.538501; cited in"
    " docs/DESIGN-MINIMAL-CELL.md since 2026-09-28). A generalized linear model with a log link, log"
    " E[expression] = b0 + sum over enhancers of b_i x_i (x_i = 1 when enhancer i is perturbed), fitted to"
    " enhancer pairs in Gasperini et al. 2019's K562 single-cell CRISPRi screen, in which enhancer effects"
    " combined multiplicatively. Read here as a composition rule: a gene's level with the set O of its"
    " elements present is T times the product, over the elements not in O, of RHO_i, each element's"
    " single-removal fold at the state with every element present. In the runtime (grn opt-in"
    " 'multiplicative', 1f8f12b): transcription = max_rate * exp(G * W) * R, W = sum of s_i*H(x_i) over the"
    " activators, R = product of (1 - s_j*H(x_j)) over the inhibitors (the runtime's inhibitor rule,"
    " unchanged), basal_rate 0"
)
PARAMETER_SOURCE = (
    "none: the form has no parameter of its own, and nothing is fitted on any CRISPRi pair, training or"
    " held-out. A gene's numbers are the bridge's, as combine2.FIT['principle'] states them: each"
    " element's removal observation RHO_i (bridge.observation_of on the predicted program) reproduced"
    " exactly at the state with every element present, and the gene's rate T (data/results/gene_rates.json,"
    " else its genome median), on which no fold depends. The scale conventions are combine2's: the"
    " strongest activator has s = 1 (DECLARED_STRENGTH) and the bridge's Hill assumptions hold (K = 1,"
    " n = 2, H(intact) = 0.5). The CRISPRi training pairs enter only as the fit rows of each stage's"
    " logistic model, as in M1c and M2d"
)
FORM_WORDS = {FORM: "the multiplicative composition (Zhou et al. 2024)"}
FIT = {
    "principle": c2.FIT["principle"],
    FORM: (
        "beta_i = -ln RHO_i for an activating element (RHO_i < 1); s_i = beta_i / max beta (the strongest"
        " has s = 1) and the gene's gain G = max beta / H, so that G*s_i*H = beta_i. An inhibiting element"
        " (RHO_j > 1) keeps the runtime's product with s_j = (1 - 1/RHO_j)/H. basal_rate = 0 and max_rate"
        " V = T*exp(-sum beta)/R0, R0 = product of (1 - s_j*H). Level relative to T with the open set O:"
        " the product of RHO_i over the elements not in O. With every element removed the level is"
        " T*exp(-sum beta)/R0 > 0: no basal share is subtracted, so the single-removal drops 1 - RHO_i may"
        " sum past 1"
    ),
    "runnable": (
        "a gene combine2.classify admits is runnable when every inhibitor's s_j is in (0, 1]"
        " (bridge.STRENGTH_TOLERANCE above 1 allowed), i.e. RHO_j <= 1/(1 - H) = 2, and sum beta <="
        " EXPONENT_LIMIT, so that exp(G*W) is a finite double; otherwise it is reported by reason and not"
        " run. Every activator's s_i is in (0, 1] by construction"
    ),
}
CENSUS = (
    "a count on the predictor side, read before any stage: of the genes the traced pairs (fit rows and"
    " population) reach, how many the multiplicative fit makes runnable (FIT['runnable']) and why the rest"
    " are not; beside it M2d's committed counts under the runtime's mean and promoter competition"
    " (data/results/s4joint_result.json), and how many of the genes combine2.classify admits any concave"
    " composition with a positive basal rate could fit at most: those whose activating single-removal"
    " drops d_i = 1 - RHO_i sum to less than 1 (WHY_NOT_SATURATING). It reads no label and decides"
    " nothing"
)
CENSUS_SENTENCE = (
    "under the multiplicative composition {r} of the {n} genes the traced pairs reach are identifiable"
    " ({x} are not: {reasons}); M2d ran {m} of them under the runtime's mean and {c} under promoter"
    " competition, and a concave composition with a positive basal rate, whatever its parameters, could"
    " fit at most the {b} whose activating single-removal drops sum to less than 1"
)
JOINT_FIT = (
    "per chromosome, the program with its experimental rule lines removed (PREDICTED_ONLY) is gated as"
    " signaltrace's s4_bridge gated it (compile.gate_program, selection reader_cells: every predicted rule"
    " is offered in K562 whatever its compiled cell) and parsed. Every gene that a rule reached by the"
    " traced pairs (fit rows and population, signaltrace.chromosome_lookup's reached rules) targets is read"
    " by combine2.gene_rules on that module: all of the gene's rules in that program, each with its removal"
    " observation and its reader states. combine2.classify's exclusions apply unchanged (every observation"
    " readable, none exactly 1, each rule's action agreeing with its observation's sign, at least one"
    " activating); its two forms' fits are not used. The multiplicative fit (FIT) gives the gene's"
    " strengths, gain and max rate; T as M2d reads it. A runnable gene is run by NetworkRuntime with"
    " activator_combination='multiplicative' in each of combine2's thirteen reader cells and ungated with"
    " each element removed in turn (HARNESS), and, in K562, once more for each of its elements open there"
    " with that element clamped absent (bridge.REMOVED), which under the form equals holding the rule"
)
HARNESS = (
    "for every gene run: (1) in each of the thirteen reader cells the gated run's level equals T times the"
    " closed form (FIT) for that cell's open set; (2) the same fitted program with its rules ungated, every"
    " element present, then each element clamped absent in turn, reproduces that element's RHO; (3) in"
    " K562, each open element's hold run equals the closed form with that element closed as well. Each"
    " within CLOSED_FORM_TOLERANCE in log2"
)
S4, S5, S6, S7 = (f"{t}_{FORM}" for t in ("s4", "s5_reader", "s6_regulator", "s7_run"))
STAGES: dict[str, dict[str, str]] = {
    **{s: st.STAGES[s] for s in ("s1_anchor", "s2_target", "s3_compiled_value")},
    S4: {
        "name": "the strength the multiplicative joint fit gives in K562",
        "follows": "s3_compiled_value",
        "kept": (
            "a rule kept by s3_compiled_value (an overlapping element's predicted rule naming the measured"
            " gene) whose gene the multiplicative fit makes runnable (JOINT_FIT), the rule among the gene's"
            " fitted rules"
        ),
        "value": (
            "max(0, -log2 RHO) of the rule's removal observation, which the fitted run reproduces at the"
            " state where every element of the gene is present (HARNESS (2)); the largest over the kept"
            " rules. It is signaltrace's s4_bridge value and M2d's s4_joint value"
        ),
    },
    S5: {
        "name": "the reader gate in K562 after the multiplicative fit",
        "follows": S4,
        "kept": (
            f"a rule kept by {S4} whose element the gated program's K562 reader block states open (an ENCODE"
            " DNase-seq narrowPeak overlaps it); closed and unknown hold the rule"
        ),
        "value": f"as {S4}, over the open rules",
    },
    S6: {
        "name": "the regulator assignment in K562 after the multiplicative fit",
        "follows": S5,
        "kept": (
            f"a rule kept by {S5} whose element signaltrace's lookup gives a source factor in K562"
            " (regulators.candidates with the rule's cell set to K562: an ENCODE TF ChIP-seq peak cluster"
            " naming K562 overlaps the element)"
        ),
        "value": f"as {S5}, over those rules",
    },
    S7: {
        "name": "the multiplicatively fitted run's predicted change in K562",
        "follows": S5,
        "kept": (
            f"a rule kept by {S5}, run: the gene's fitted program (fitted_module on"
            " readergate.gated_gene_module) run by runtime.grn.NetworkRuntime with"
            " activator_combination='multiplicative' in K562, every element clamped present so that the"
            " reader holds the closed ones"
        ),
        "value": (
            "max(0, log2(L_K562 / L_held)): the gene's level in K562 with the reader applied, over its level"
            " with this rule held as well (its element clamped absent); the largest over the run rules"
        ),
    },
}
MULT_STAGES = (S4, S5, S6, S7)
CHAIN = ("s1_anchor", "s2_target", "s3_compiled_value", S4, S5, S7)
BRANCHES = (S6,)
COMPARATORS = sj.COMPARATORS
ANCHOR = st.ANCHOR
COMPARATOR = (
    "the anchor and signaltrace's s2_target and s3_compiled_value, recomputed by signaltrace's own code on"
    " the same pairs and fit rows (s4joint.comparators); the run stops (HARNESS_DEFECT) unless the anchor"
    " equals its committed arm-1 figures and s2_target and s3_compiled_value equal signaltrace's committed"
    " readings exactly (data/results/signaltrace_result.json: kept counts, gain over activity and distance,"
    " interval, and against the anchor with its interval). Every multiplicative stage is read against the"
    " anchor and against the stage it follows by the same paired resampling; signaltrace's s4_bridge and"
    " chain readings and M2d's per-form bridge and chain readings are quoted beside it, verbatim"
)
POPULATION = sj.POPULATION
PREDICTED_ONLY = st.PREDICTED_ONLY
BASE_COLUMNS = st.BASE_COLUMNS
METRIC = st.METRIC
BOOTSTRAP = st.BOOTSTRAP
DECISION = st.DECISION
LOSS = st.LOSS
BRIDGE_OF = sj.BRIDGE_OF
BRIDGE_DECISION = (
    f"from {S4}'s outcome under DECISION: SURVIVES if CARRIES or WEAKENED (the 95% interval of its gain over"
    " activity and distance lies above zero, so the chain passes the bridge), LOST if NONE_KEPT, LOST or"
    " WORSE (signaltrace's LOSS), UNDECIDED if UNDECIDED. 'Survives' means exactly that; whether the gain is"
    " below the anchor's is the stage reading's own clause, kept beside it. The chain (CHAIN) is walked as"
    f" signaltrace walks it, with {S4} in s4_bridge's place, {S5} in s5_reader's and {S7} in s7_run's, and"
    f" {S6} read beside it on {S5}'s pairs"
)
FALSIFIER = (
    "the claim under test: the bridge's loss in M1c (NONE_KEPT) and M2d (LOST under both forms) came from"
    " a composition that could not fit the reached genes, not from the signal the compiled rules carry. It"
    f" is falsified if {S4}'s bridge answer is LOST (outcome NONE_KEPT, LOST or WORSE): with the genes"
    " fittable, the measured signal still does not pass the bridge. It is not falsified if the answer is"
    " SURVIVES and not decided if UNDECIDED. No answer tests the multiplicative form against a measurement"
    " (BY_CONSTRUCTION (4))"
)
BY_CONSTRUCTION = (
    "(1) every gene combine2.classify admits whose inhibitors each have RHO_j <= 2 and whose summed"
    " log-folds stay within EXPONENT_LIMIT is runnable, whatever its observations: the census counts, it"
    " does not test",
    f"(2) {S4}'s value is the rule's own observation, max(0, -log2 RHO), as s4_bridge's and s4_joint's are;"
    " it differs from s3_compiled_value's (the compiled |log2 fold|, clipped at 1) only above the clip,"
    f" and {S4} keeps s3's pairs less those whose gene is not runnable. A reading at {S4} close to"
    " s3_compiled_value's is therefore expected and is not evidence for the form",
    "(3) under the form, holding an element changes the level by that element's own factor whatever the"
    f" other elements' states, so {S7}'s two columns equal {S5}'s to within CLOSED_FORM_TOLERANCE: the run"
    " checks the runtime, it adds no information to the reader stage",
    "(4) the CRISPRi labels are single removals and every fit reproduces every single removal, so no"
    " reading here can tell the multiplicative composition from another composition that fits the same"
    " genes; the evidence for the form is Zhou et al. 2024's, not this run's",
)
CLOSED_FORM_TOLERANCE = c2.CLOSED_FORM_TOLERANCE
STOP_RULE = (
    "one run, from this registration's committed code, on the cached files the registration pins by"
    " sha256, signaltrace's registration and result and M2d's result among them. No stage, column,"
    " population, fit row, form, metric or wording is changed after any multiplicative stage is read; a"
    " defect found afterwards is reported beside the registered readings and corrected only under a new"
    " registration. The run stops with HARNESS_DEFECT, before any multiplicative stage is read, if: the"
    " population differs from signaltrace's pinned list; signaltrace's registration differs from"
    " signaltrace.py; signaltrace's own stage columns, recomputed, differ on any population pair from its"
    " committed pairs table; signaltrace's own harness fails; a reached gene is not read from the gated"
    " program, or its K562 reader state differs from signaltrace's lookup; any gene's run misses HARNESS; a"
    " multiplicative stage keeps a pair the stage it follows dropped; the anchor does not reproduce the"
    " committed arm-1 figures; or s2_target or s3_compiled_value does not reproduce signaltrace's committed"
    " reading"
)
OUTCOMES = st.OUTCOMES
CHAIN_OUTCOMES = st.CHAIN_OUTCOMES
BEFORE_CLAUSE = st.BEFORE_CLAUSE
BRIDGE_OUTCOMES = sj.BRIDGE_OUTCOMES
FORM_READING = sj.FORM_READING
LIMITATIONS = (
    *st.LIMITATIONS,
    c2.LIMITATIONS[1],
    "the form is chosen in writing on the grounds of WHY_NOT_SATURATING and of its citation, not on a"
    " measurement here: Zhou et al. 2024 is evidence from enhancer pairs in another K562 screen (Gasperini"
    " 2019), and BY_CONSTRUCTION (4) says why this run cannot add to it",
    "the multiplicative fit sets the basal rate to 0: a gene's level with every element removed is"
    " T*exp(-sum beta)/R0, for a gene with many elements a small fraction of T. No reading depends on that"
    " level and no level in molecules per cell is quoted",
    "a pair whose gene the form cannot fit is dropped at the multiplicative bridge with 0 in the stage's two"
    " columns, as the anchor treats an uncovered pair",
    f"under a composition that admits nearly every gene, {S4} differs from s3_compiled_value only by the"
    f" pairs it drops and the values above the compiled clip (BY_CONSTRUCTION (2)); what no earlier"
    f" registered stage has read on the target-choice pairs is {S5} and {S6}, the reader and the regulator"
    " in K562",
    "this lane read M1c's and M2d's modules, scripts and committed results before registering"
    " (EXPLORATION_PRECEDED_THIS in the registration): the comparators and M2d's census are known; the"
    " multiplicative stages on the real pairs are not",
)
DESCRIPTIVES = (
    "the census (CENSUS) with its reasons, per chromosome; over the reached genes classify admits, the"
    " quantiles of the summed activating drops (sum of 1 - RHO_i, the mean's q) and of the summed"
    " log-folds (sum beta), and the rule counts of the runnable genes",
    "per stage: the pairs kept of the 1,918 and of the fit rows, how many of them regulated, how many with"
    " a value above 0, the reasons the pairs it drops were dropped, the fitted weights of its model, and its"
    " paired interval against the stage it follows",
    "the run checks: gene runs, K562 hold runs, the largest difference from the closed form in log2",
    "signaltrace's committed s4_bridge and chain readings and M2d's per-form readings, quoted verbatim for"
    " context",
)
#: the registration fields the run refuses to differ from: registration field -> constant of this module
BOUND = {
    "question": "QUESTION",
    "why_not_saturating": "WHY_NOT_SATURATING",
    "form": "FORM",
    "form_source": "FORM_SOURCE",
    "parameter_source": "PARAMETER_SOURCE",
    "form_words": "FORM_WORDS",
    "fit": "FIT",
    "exponent_limit": "EXPONENT_LIMIT",
    "census": "CENSUS",
    "census_sentence": "CENSUS_SENTENCE",
    "population": "POPULATION",
    "predicted_only": "PREDICTED_ONLY",
    "cell": "CELL",
    "joint_fit": "JOINT_FIT",
    "harness": "HARNESS",
    "stages": "STAGES",
    "chain": "CHAIN",
    "branches": "BRANCHES",
    "anchor": "ANCHOR",
    "comparator": "COMPARATOR",
    "base_columns": "BASE_COLUMNS",
    "metric": "METRIC",
    "bootstrap": "BOOTSTRAP",
    "decision": "DECISION",
    "loss": "LOSS",
    "bridge_of": "BRIDGE_OF",
    "bridge_decision": "BRIDGE_DECISION",
    "falsifier": "FALSIFIER",
    "by_construction": "BY_CONSTRUCTION",
    "closed_form_tolerance": "CLOSED_FORM_TOLERANCE",
    "stop_rule": "STOP_RULE",
    "outcomes": "OUTCOMES",
    "chain_outcomes": "CHAIN_OUTCOMES",
    "before_clause": "BEFORE_CLAUSE",
    "bridge_outcomes": "BRIDGE_OUTCOMES",
    "form_reading": "FORM_READING",
    "limitations": "LIMITATIONS",
    "descriptives_not_verdicts": "DESCRIPTIVES",
}
QUESTION = (
    "M1c lost the measured CRISPRi signal at the bridge in K562 (no reached gene identifiable) and M2d lost"
    " it there again with combine2's joint fit, which could run 10 (mean) and 2 (promoter competition) of"
    " the 210 reached genes, the rest mostly because a gene's single-removal drops sum past its whole"
    " expression. Under a composition that can fit such genes - the multiplicative composition of Zhou et"
    " al. 2024, with no parameter of its own - (a) how many of the reached genes become identifiable, and"
    " (b) on the same 1,918 pinned held-out K562 pairs, with the same metric, bootstrap and decision as M1c"
    " and M2d, does the agreement with the measured CRISPRi outcome survive the bridge, and where along"
    " the reader, the regulator and the run does it go?"
)


def registered_fields() -> dict[str, Any]:
    """Every BOUND constant as JSON stores it, so a registration and this module compare equal or not."""
    return {f: json.loads(json.dumps(globals()[c])) for f, c in BOUND.items()}


# ==== the multiplicative fit (closed form; no label) =======================================================


@dataclass
class MFit:
    """One gene's multiplicative fit, relative to its intact transcription T."""

    ok: bool
    reason: str = ""
    rho: dict[str, float] = field(default_factory=dict)  # element -> RHO
    beta: dict[str, float] = field(default_factory=dict)  # activator -> -ln RHO
    s: dict[str, float] = field(default_factory=dict)  # element -> fitted strength
    gain: float = math.nan  # G
    vmax_frac: float = math.nan  # max_rate / T
    summed_drop: float = math.nan  # sum of 1 - RHO_i over the activators (the mean's q)
    summed_beta: float = math.nan


def summed_drop(rho: dict[str, float]) -> float:
    return sum(1.0 - r for r in rho.values() if r < 1.0)


def fit(rho: dict[str, float]) -> MFit:
    """The multiplicative fit of one gene's elements from their removal folds (FIT)."""
    beta = {el: -math.log(r) for el, r in rho.items() if r < 1.0}
    inh = {el: r for el, r in rho.items() if r > 1.0}
    base = {"rho": dict(rho), "summed_drop": summed_drop(rho), "summed_beta": sum(beta.values())}
    if len(beta) + len(inh) != len(rho):
        return MFit(False, "observation_null_effect", **base)
    if not beta:
        return MFit(False, "no_activating_rule", **base)
    s: dict[str, float] = {}
    r0 = 1.0
    for el, r in inh.items():
        s_j = c2._strength_ok((1.0 - 1.0 / r) / H)
        if s_j is None:
            return MFit(False, "inhibitor_strength_above_1", **base)
        s[el] = s_j
        r0 *= 1.0 - s_j * H
    if not base["summed_beta"] <= EXPONENT_LIMIT:
        return MFit(False, "exponent_beyond_float_range", **base)
    top = max(beta.values())
    for el, b in beta.items():
        s[el] = b / top
    return MFit(True, "", beta=beta, s=s, gain=top / H, vmax_frac=math.exp(-base["summed_beta"]) / r0, **base)


def level(f: MFit, open_: dict[str, bool]) -> float:
    """Transcription relative to T with the elements in `open_` open and the rest held (the closed form)."""
    return math.prod(r for el, r in f.rho.items() if not open_.get(el, False))


def fitted_module(small: Module, gene: str, f: MFit, rate: float, label: dict[str, str]) -> Module:
    """The gene's gated module with the fit written in: basal 0, max rate, gain, every rule's strength and
    action, the bridge's Hill assumptions and the rate label. `through` is kept."""
    m = copy.deepcopy(small)
    g = m.entities[gene]
    g.basal_rate = 0.0
    g.attrs["max_rate"] = float(rate * f.vmax_frac)
    g.attrs[grn.MULTIPLICATIVE_GAIN] = float(f.gain)
    for attr, k in zip(grn.RATE_PROVENANCE_ATTRS, ("source", "species_cell", "tier"), strict=True):
        g.attrs[attr] = label[k]
    for r in m.rules:
        if r.target != gene or r.action not in c2.REGULATORY:
            continue
        r.strength = f.s[r.through]
        r.threshold = bridge.HILL_THRESHOLD
        r.hill = bridge.HILL_COEFFICIENT
        r.action = Action.ACTIVATE if r.through in f.beta else Action.INHIBIT
    return m


def _off(run: float, want: float) -> float:
    return abs(math.log2(run / want))


def run_gene(small: Module, gr: c2.GeneRules, f: MFit, rate: float, label: dict[str, str]) -> dict[str, Any]:
    """The gene's runs under the form (HARNESS (1)-(3)) and its K562 hold changes."""
    m = fitted_module(small, gr.gene, f, rate, label)
    present = dict.fromkeys(gr.elements, bridge.INTACT)
    delta = math.log(2) / grn.DEFAULTS["mrna_half_life"]
    miss: list[str] = []
    worst = 0.0
    levels = []
    for j, cid in enumerate(c2.CELLS):
        lv = c2._run_level(m, gr.gene, {"cell_type": cid}, present, FORM)
        want = rate * level(f, gr.open_in(j)) / delta
        levels.append(lv)
        worst = max(worst, _off(lv, want))
        if not _off(lv, want) <= CLOSED_FORM_TOLERANCE:
            miss.append(f"cell {cid}: run {lv:.12g}, closed form {want:.12g}")
    ungated = copy.deepcopy(m)
    for r in ungated.rules:
        r.through = ""
    ctx = {"cell_type": c2.CELLS[0]}
    full = c2._run_level(ungated, gr.gene, ctx, present, FORM)
    for el in gr.elements:
        lv = c2._run_level(ungated, gr.gene, ctx, {**present, el: bridge.REMOVED}, FORM)
        off = abs(math.log2(lv / full) - math.log2(gr.rho[el]))
        worst = max(worst, off)
        if not off <= CLOSED_FORM_TOLERANCE:
            miss.append(f"removal {el}: run {lv / full:.12g}, observed {gr.rho[el]:.12g}")
    open_k = gr.open_in(K562_INDEX)
    level_k = levels[K562_INDEX]
    change: dict[str, float] = {}
    for el in gr.elements:
        if not open_k[el]:
            continue
        held = c2._run_level(m, gr.gene, {"cell_type": CELL}, {**present, el: bridge.REMOVED}, FORM)
        want = rate * level(f, {**open_k, el: False}) / delta
        worst = max(worst, _off(held, want))
        if not _off(held, want) <= CLOSED_FORM_TOLERANCE:
            miss.append(f"K562 hold {el}: run {held:.12g}, closed form {want:.12g}")
        change[el] = math.log2(level_k / held)
    return {
        "ok": True,
        "reason": "",
        "tier": label["tier"],
        "summed_drop": f.summed_drop,
        "summed_beta": f.summed_beta,
        "change": change,
        "harness_miss": miss,
        "max_abs_run_minus_closed_form_log2": worst,
        "gene_runs": len(c2.CELLS) + 1 + len(gr.elements),
        "hold_runs": len(change),
    }


def joint_lookup(
    chrom: str,
    text: str,
    look: st.ChromLookup,
    rates: dict[str, float],
    borrowed: float,
    prov: dict | None,
    lent: dict | None,
    readers: Any,
    results_dir: Any,
) -> sj.JointLookup:
    """Every reached gene of one chromosome fitted and run under the form (no label is read)."""
    from genomeos.lang.parser import parse

    out = sj.JointLookup(chrom)
    reached = set(look.bridge_reader) | set(look.why_reader)
    out.reached_genes = {gene for _el, gene in reached}
    out.counts["reached_rules"] = len(reached)
    out.counts["reached_genes"] = len(out.reached_genes)
    if not reached:
        return out
    pred, _removed = st.predicted_text(text)
    gated, _info = cp.gate_program(
        pred, chrom, results_dir, selection=cp.SELECTION_READER_CELLS, readers=readers
    )
    module = parse(gated, f"noncoding_{chrom}_reader")
    for gr in c2.gene_rules(module, out.reached_genes, chrom, rgate.compiled_cells(pred)):
        _fits, why = c2.classify(gr)
        g = sj.GeneFit(
            gr.gene,
            gr.symbol,
            list(gr.elements),
            dict(gr.rho),
            {el: gr.states[el][K562_INDEX] for el in gr.elements},
            why,
        )
        if why:
            g.forms[FORM] = {"ok": False, "reason": why}
        else:
            f = fit(gr.rho)
            if not f.ok:
                g.forms[FORM] = {
                    "ok": False,
                    "reason": f.reason,
                    "summed_drop": f.summed_drop,
                    "summed_beta": f.summed_beta,
                }
            else:
                rate, label = sj.rate_of(gr.gene, rates, borrowed, prov, lent)
                g.forms[FORM] = run_gene(rgate.gated_gene_module(module, gr.gene), gr, f, rate, label)
                out.counts["runs"] += 1
        out.genes[gr.gene] = g
    del module
    return out


# ==== one pair through the multiplicative stages ===========================================================


def trace(
    p: crispri.Pair, look: st.ChromLookup | None, joint: sj.JointLookup | None
) -> tuple[dict[str, float], dict[str, str]]:
    """The two columns of every multiplicative stage for one pair, and why each stage that drops it did.

    Reads the chromosome's two lookups and never the pair's label."""
    f: dict[str, float] = {}
    why: dict[str, str] = {}
    named = look.named(p.start, p.end, p.gene) if look else []
    fitted = [(lk, joint.fitted(FORM, lk)) for lk in named] if joint else []
    b = [(lk, g) for lk, g in fitted if g is not None]
    f[f"{S4}_kept"] = float(bool(b))
    f[f"{S4}_drop"] = max((st._drop(g.rho[lk.element]) for lk, g in b), default=0.0)
    if named and not b:
        why[S4] = sorted(joint.reason(FORM, lk) for lk in named)[0] if joint else "no_joint_lookup"
    o = [(lk, g) for lk, g in b if g.k562[lk.element] == "open"]
    f[f"{S5}_kept"] = float(bool(o))
    f[f"{S5}_drop"] = max((st._drop(g.rho[lk.element]) for lk, g in o), default=0.0)
    if b and not o:
        why[S5] = "reader_" + sorted(g.k562[lk.element] for lk, g in b)[0]
    h = [(lk, g) for lk, g in o if look.factors.get(lk.element)]
    f[f"{S6}_kept"] = float(bool(h))
    f[f"{S6}_drop"] = max((st._drop(g.rho[lk.element]) for lk, g in h), default=0.0)
    if o and not h:
        why[S6] = "no factor cluster names K562 at the element"
    runs = [g.forms[FORM]["change"][lk.element] for lk, g in o if lk.element in g.forms[FORM]["change"]]
    f[f"{S7}_kept"] = float(bool(runs))
    f[f"{S7}_drop"] = max((max(0.0, c) for c in runs), default=0.0)
    if b and not runs:
        why[S7] = "held in K562"
    return f, why


# ==== the harness and the census, before any multiplicative stage is read ===================================


def harness(
    lookups: dict[str, st.ChromLookup], joints: dict[str, sj.JointLookup], pairs: Sequence[crispri.Pair]
) -> list[str]:
    """Every predictor-side check of STOP_RULE on the multiplicative stages that reads no label."""
    problems = []
    for chrom, joint in sorted(joints.items()):
        look = lookups[chrom]
        missing = sorted(joint.reached_genes - set(joint.genes))
        if missing:
            problems.append(f"{chrom}: reached genes not read from the gated program: {missing[:5]}")
        for el, gene in sorted(set(look.bridge_reader) | set(look.why_reader)):
            g = joint.genes.get(gene)
            if g is not None and el in g.k562 and g.k562[el] != look.reader.get(el):
                problems.append(f"{chrom} {el}: K562 state {g.k562[el]} against {look.reader.get(el)}")
        for gene, g in sorted(joint.genes.items()):
            fr = g.forms.get(FORM, {})
            if fr.get("harness_miss"):
                problems.append(f"{chrom} {gene}: {'; '.join(fr['harness_miss'][:3])}")
    for stage in MULT_STAGES:
        parent = STAGES[stage]["follows"]
        stray = sum(1 for p in pairs if p.features[f"{stage}_kept"] and not p.features[f"{parent}_kept"])
        if stray:
            problems.append(f"{stage} keeps {stray} pairs that {parent} dropped")
    return problems


def _ranked(c: Counter) -> list[tuple[str, int]]:
    """(reason, count), the most common first and ties in alphabetical order, so the words are stable."""
    return sorted(c.items(), key=lambda kv: (-kv[1], kv[0]))


def _reasons(c: Counter) -> str:
    return ", ".join(f"{n} {why}" for why, n in _ranked(c)) or "none"


def census(joints: dict[str, sj.JointLookup], m2d: dict[str, Any]) -> dict[str, Any]:
    """CENSUS: identifiable genes under the form, M2d's counts beside it, and the concave bound."""
    genes = [g for j in joints.values() for g in j.genes.values()]
    ok = [g for g in genes if g.forms[FORM]["ok"]]
    why = Counter(g.forms[FORM]["reason"] for g in genes if not g.forms[FORM]["ok"])
    admitted = [g for g in genes if not g.excluded]
    bound = sum(1 for g in admitted if summed_drop(g.rho) < 1.0)
    per = m2d["descriptives_not_verdicts"]["per_form"]
    out = {
        "reached_genes": len(genes),
        "identifiable": len(ok),
        "not_identifiable_by_reason": dict(_ranked(why)),
        "admitted_by_combine2_classify": len(admitted),
        "concave_bound_summed_drop_below_1": bound,
        "m2d_runnable": {form: per[form]["runnable"] for form in ("mean", "competition")},
    }
    out["sentence"] = CENSUS_SENTENCE.format(
        r=f"{len(ok):,}",
        n=f"{len(genes):,}",
        x=f"{len(genes) - len(ok):,}",
        reasons=_reasons(why),
        m=f"{out['m2d_runnable']['mean']:,}",
        c=f"{out['m2d_runnable']['competition']:,}",
        b=f"{bound:,}",
    )
    return out


# ==== the readings ========================================================================================


def stage_row(
    stage: str,
    fit_rows: Sequence[crispri.Pair],
    pop: Sequence[crispri.Pair],
    anchor: dict[str, Any],
    scores: dict[str, list[float]],
    weights: dict[str, list[float]],
) -> dict[str, Any]:
    """One stage's reading, built as signaltrace.score builds it (METRIC, DECISION, OUTCOMES)."""
    spec = STAGES[stage]
    row: dict[str, Any] = {
        "name": spec["name"],
        "follows": spec["follows"],
        "columns": list(st.columns(stage)),
        "weights": [round(x, 6) for x in weights[stage]],
        "population": st.kept_counts(pop, stage),
        "fit_rows": st.kept_counts(fit_rows, stage),
        "gain_over_activity_distance": crispri.gain_interval(scores[stage], anchor["base"], list(pop)),
        "against_anchor": crispri.gain_interval(scores[stage], scores["s1_anchor"], list(pop)),
        "against_the_stage_it_follows": crispri.gain_interval(
            scores[stage], scores[spec["follows"]], list(pop)
        ),
    }
    row["outcome"] = st.decide(
        row["population"]["kept"], row["gain_over_activity_distance"]["ci95"], row["against_anchor"]["ci95"]
    )
    row["reading"] = st.reading(row["outcome"], **st.stage_numbers(row))
    return row


def chain_verdict(stages: dict[str, Any]) -> tuple[str, str, str | None]:
    """(outcome, its sentence, the stage it names) for CHAIN: signaltrace.chain_verdict's walk."""
    first_weak = None
    for stage in CHAIN[1:]:
        row = stages[stage]
        if row["outcome"] == "WEAKENED" and first_weak is None:
            first_weak = stage
        if row["outcome"] in ("CARRIES", "WEAKENED"):
            continue
        before = BEFORE_CLAUSE.format(stage=first_weak, name=STAGES[first_weak]["name"]) if first_weak else ""
        kind = "LOST_AT" if row["outcome"] in LOSS else "UNDECIDED_AT"
        text = CHAIN_OUTCOMES[kind].format(
            stage=stage, name=STAGES[stage]["name"], stage_reading=row["reading"], before=before
        )
        return kind, text, stage
    run = stages[CHAIN[-1]]["reading"]
    if first_weak:
        text = CHAIN_OUTCOMES["WEAKENED_NOT_LOST"].format(
            stage=first_weak, name=STAGES[first_weak]["name"], run_reading=run
        )
        return "WEAKENED_NOT_LOST", text, first_weak
    return "CARRIED_TO_THE_RUN", CHAIN_OUTCOMES["CARRIED_TO_THE_RUN"].format(run_reading=run), None


def bridge_verdict(stages: dict[str, Any]) -> tuple[str, str]:
    row = stages[S4]
    kind = BRIDGE_OF[row["outcome"]]
    text = BRIDGE_OUTCOMES[kind].format(
        form=FORM_WORDS[FORM], stage=S4, name=STAGES[S4]["name"], stage_reading=row["reading"]
    )
    return kind, text


def score(
    fit_rows: Sequence[crispri.Pair],
    pop: Sequence[crispri.Pair],
    anchor: dict[str, Any],
    comp_rows: dict[str, Any],
    comp_scores: dict[str, list[float]],
    comp_weights: dict[str, list[float]],
) -> dict[str, Any]:
    """Every multiplicative stage's reading, the bridge answer, the chain verdict and the reading."""
    weights = {"s1_anchor": anchor["weights"], **comp_weights}
    scores = {"s1_anchor": anchor["scores"], **comp_scores}
    for s in MULT_STAGES:
        weights[s], scores[s] = st.fit_scores(fit_rows, list(pop), st.columns(s))
    stages = dict(comp_rows)
    for s in MULT_STAGES:
        stages[s] = stage_row(s, fit_rows, pop, anchor, scores, weights)
    bkind, btext = bridge_verdict(stages)
    ckind, ctext, at = chain_verdict(stages)
    return {
        "stages": stages,
        "bridge": bkind,
        "bridge_reading": btext,
        "chain": ckind,
        "chain_reading": ctext,
        "stage_where_it_stops": at,
        "reading": FORM_READING.format(bridge=btext, form=FORM_WORDS[FORM], chain=ctext),
    }
