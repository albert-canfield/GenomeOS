# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-combine2 (ROADMAP §5 items 1 and 7, 2026-10-10): genes driven by several elements, run.

M1 (`runrules.py`, 747b207) and M2b (`readergate.py`, ff35b09) ran only the 738 of 18,973 genes that
carry ONE compiled rule: `bridge.parameterize` refuses a gene with two or more mechanisms in a cell
(`not_identifiable`), and under the reader selection every predicted rule of a gene is offered in every
reader cell, so every multi-element gene was outside both tests. R3 (lane-combine, 93caf61) left the
activator combination "mean retained, undetermined" and the strength-dependent (promoter-competition)
form open as its own registration. This module is that registration and the code that runs it.

What makes a multi-element gene runnable (`STRUCTURE`). The runtime's removal state is a clamp (the
bridge's REMOVED, the rule kept with H(0) = 0) and the reader's hold is a strength-0 copy of the rule,
kept. Under both, a held or removed activator stays in the mean's denominator, so the mean is LINEAR in
the element terms: removing element i lowers transcription by V*h*s_i*R0/n whatever the other terms are.
k single-removal observations taken at one intact state therefore identify the k contributions jointly
(`fit`), with the bridge's own conventions for what they leave free. The 2026-09-28 non-identifiability
(dA_i = (A - t_i)/(n - 1)) is the change when a rule is DROPPED from the mean, which neither the clamp nor
the reader does. What is not identified, and is an assumption here, is that a gene's k AlphaGenome
observations were taken at one state: each is the argmax over 371 tracks, in its own track.

This module holds, in the order they were written and committed:

1. the joint fit under each candidate form, the census and the code that runs and scores (no RNA);
2. the registered test (`QUESTION` .. `OUTCOMES`), committed before any RNA value was read.

No model request, no key, no money. The measured RNA is M2b's: ENCODE polyA plus RNA-seq of the thirteen
reader biosamples, fetched by scripts/readergate_fetch.py and pinned by sha256.
"""

from __future__ import annotations

import copy
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Any

from genomeos.attribution import bridge
from genomeos.attribution import readergate as rg
from genomeos.attribution import runrules as rr
from genomeos.ir import Action, Module
from genomeos.runtime import grn
from genomeos.runtime.grn import NetworkRuntime

REGISTERED = "2026-10-10"
LANE = "lane-combine2"
SELECTION = rg.SELECTION
CELLS = rg.CELLS
REGULATORY = rg.REGULATORY
H = bridge.H_INTACT

# ---- the candidate combination forms, fixed before any RNA value was read ---------------------------
FORMS = ("mean", "competition")
#: competition's one declared number: the promoter's occupancy with every element of the gene present,
#: W/(K + W). Not measured for any human gene. Set to H_INTACT, the half-saturation the bridge already
#: declares for every element at its intact state, which puts it midway between the two regimes Bothma
#: et al. 2015 describe (theta -> 0: weak engagement, the elements add, which is the mean under the
#: hold; theta -> 1: saturated, the strongest element dominates)
THETA = 0.5
FORM_SOURCES = {
    "mean": (
        "runtime/grn.py's default since the first commit (ACTIVATOR_COMBINATION = 'mean'), as the runtime"
        " runs it with the reader: A = (1/n) * sum of s_i*H(x_i) over the gene's n activating rules, a held"
        " rule counted at strength 0; inhibitors multiply, R = prod (1 - s_j*H(x_j)). Under the hold the"
        " mean is linear in the element terms: closing elements adds their contributions"
    ),
    "competition": (
        "grn OPT_IN_ACTIVATOR_RULES 'competition': A = W/(K + W), W = sum of s_i*H(x_i). Source: Bothma JP,"
        " Garcia HG, Ng S, Perry MW, Gregor T, Levine M (2015) Enhancer additivity and non-additivity are"
        " determined by enhancer strength in the Drosophila embryo. eLife 4:e07956, doi"
        " 10.7554/eLife.07956, PMC4532966, Methods 'Mathematical modeling' eqs. 7-9 (a promoter engaged"
        " by at most one enhancer at a time; read here for n enhancers). Declared, not from the paper:"
        " one transcriptional efficiency r for every element of a gene (the paper's eq. 9 lets each"
        " have its own; with one r the form is monotone and cannot give snail's pair below its stronger"
        " single); the promoter occupancy THETA = 0.5 with every element present; and the runtime's"
        " basal rate, which the paper's model does not have. Inhibitors multiply, as under the mean"
    ),
}
FIT = {
    "principle": (
        "the bridge's: the run reproduces each element's own removal observation RHO_i = 2**(log2 fold"
        " change on deletion) exactly at the state where every element of the gene is present (S5,"
        " 'direction and fold ... are the input'), and the gene's intact transcription is its rate T"
        " (measured, else the borrowed median; no fold prediction depends on T)"
    ),
    "mean": (
        "d_i = 1 - RHO_i (activating), e_j = RHO_j - 1 (inhibiting). q = sum d_i is the element-driven"
        " share of intact transcription, basal b = T(1 - q); an inhibitor's s_j*h = e_j/(q + e_j); an"
        " activator's s_i = d_i/max d (DECLARED_STRENGTH convention: the strongest has s = 1) with"
        " max_rate V = n * max d * T / (h * R0). Level relative to T with the open set O:"
        " (1 - q) + [sum of d_i over open activators] * prod over closed inhibitors (1 + e_j/q)"
    ),
    "competition": (
        "K = THETA/(1 - THETA); weights kappa_i (sum K) and Z = V*R0/T solve"
        " d_i = Z*kappa_i/((1 + K)(1 + K - kappa_i)), i.e. sum over i of d_i(1 + K)^2/(Z + d_i(1 + K)) = K,"
        " one root Z > 0 (brentq). q = Z*THETA, b = T(1 - q), s_j*h = e_j/(q + e_j); s_i = kappa_i/max"
        " kappa and the gene's promoter_half_saturation = h/max kappa (DECLARED_STRENGTH convention)."
        " Level relative to T: (1 - q) + Z * W/(1 + W) * prod over closed inhibitors (1 + e_j/q), W the"
        " sum of kappa_i over open activators"
    ),
    "runnable": (
        "a gene is runnable under a form when the fit gives b > 0 (q < 1) and every fitted s in (0, 1]"
        " (bridge.STRENGTH_TOLERANCE above 1 allowed); otherwise it is reported by reason and not run"
    ),
}
STRUCTURE = {
    "mean_is_linear_under_the_hold": (
        "Module.active_rules returns a held rule as a strength-0 copy and keeps it, and the bridge's"
        " removal is a clamp of the element at 0 with the rule kept; under both the activator stays in"
        " the mean's n. So removing or holding element i lowers transcription by V*h*s_i*R0/n, independent"
        " of the other terms, and the mean's joint fit (FIT['mean']) is exact. The 2026-09-28"
        " non-identifiability (dA_i = (A - t_i)/(n - 1)) is the DROP reading, which neither does"
    ),
    "what_is_assumed": (
        "that a gene's k observations were taken at one state. Under the reader selection a gene's rules"
        " come from k different argmax tracks; the joint fit treats them as one intact state, and each"
        " reader cell's run transports them there, as M2b transported one rule's m"
    ),
    "inhibitor_only": (
        "a gene with no activating rule has A = 1 under every activator form, so the forms cannot differ"
        " on it and its q is fixed by no observation; it is not run"
    ),
    "network": "no compiled rule's source is a gene (runrules census): one independent system per gene",
}

# ---- the registered test (written before any RNA value was read) --------------------------------------
QUESTION = (
    "For genes the compiled programs drive through two or more elements, does a reader-gated run under a"
    " registered combination form (the runtime's mean; the promoter-competition form) track measured RNA"
    " across the thirteen reader cell types better than a no-rule baseline, and does the strength-"
    "dependent form track it better than the mean?"
)
POPULATION = (
    "for each form F, every gene G such that, over all 24 compiled programs (data/knowledge/compiled,"
    " pinned by sha256) gated with compile.gate_program(selection='reader_cells'): (a) two or more"
    " activates/inhibits rule lines target G, all in one program, all of predicted evidence, each acting"
    " through its own element; (b) every rule's removal observation is readable (bridge.observation_of)"
    " and none is exactly 1; (c) at least one rule is activating; (d) the joint fit under F is runnable"
    " (FIT['runnable']); (e) the gated prediction under F is not the same in all thirteen reader cells"
    " (its range exceeds INFORMATIVE_LOG2 in log2); (f) G's symbol names exactly one gene_id in"
    " data/cache/gencode_genes.tsv and that gene_id is on exactly one row of every one of the thirteen"
    " quantification tables. (a)-(e) are predictor-side: no measured value enters them. Each G"
    " contributes its thirteen (G, cell) pairs. No filter on expression"
)
INFORMATIVE_LOG2 = 1e-9
MIN_GENES = 100
PSEUDOCOUNT = rg.PSEUDOCOUNT
MEASURED = rg.MEASURED
RNA_SELECTION = rg.RNA_SELECTION
PREDICTIONS = {
    "form": (
        "s_F(G, X) = log2(L_F(G, X)), L_F(G, X) the level the runtime gives G in reader cell X under form F"
        " (NetworkRuntime(..., activator_combination=F), the fitted rules acting through their elements,"
        " the reader blocks in force, every element clamped present so that only the reader holds a"
        " rule). P_F(G, X) = s_F(G, X) - mean over the thirteen cells of s_F(G, X')"
    ),
    "no_rule_baseline": "P(G, X) = 0: G's level in X equals its mean over the thirteen reader cells",
}
METRIC = (
    "Delta MSE(F) = mean over the population's pairs of D^2 - mean of (D - P_F)^2 (rg.METRIC), positive"
    " when the form's squared error is below the baseline's; an uninformative prediction still costs"
    " mean(P^2), so it reads WORSE once the population resolves that cost. The head-to-head is"
    " MSE(mean) - MSE(competition) over the genes in both forms' populations, positive when the"
    " competition form is closer"
)
#: two forms are judged, so each form's interval is Bonferroni-adjusted: 97.5%, percentiles 1.25-98.75
FORM_INTERVAL = (1.25, 98.75)
HEAD_TO_HEAD_INTERVAL = (2.5, 97.5)
BOOTSTRAP = {
    "resamples": 2000,
    "seed": 20261010,
    "unit": "gene: genes are resampled with replacement and carry all thirteen of their pairs",
    "interval": "per form, the 97.5% percentile interval (Bonferroni over the two forms); head-to-head 95%",
}
CONTROL = (
    "shuffled-reader control, per form: the thirteen-cell open/held rows of the population's elements are"
    " permuted among all its elements (PERMUTATIONS times); every element keeps its own observation and"
    " role and every gene its fit, so the control keeps each cell's open count and the distribution of"
    " rows and breaks only the tie between a gene's elements and their own chromatin. The permuted"
    " prediction is the form's closed form (FIT), which the run must equal (HARNESS). control p = (1 +"
    " permutations with Delta MSE(F) >= observed) / (1 + PERMUTATIONS)"
)
PERMUTATIONS = {"count": 1000, "seed": 20261011}
CONTROL_ALPHA = 0.05
RUNTIME = rr.RUNTIME
#: the run must equal its closed form, and reproduce every single-removal observation, in log2
CLOSED_FORM_TOLERANCE = rr.CLOSED_FORM_TOLERANCE
HARNESS = (
    "for every gene run under every form: (1) in each of the thirteen reader cells the gated run's level"
    " equals T times the form's closed form for that cell's open set; (2) the same fitted program with"
    " its rules ungated, every element present, then each element clamped absent in turn, reproduces"
    " that element's RHO. Each within CLOSED_FORM_TOLERANCE in log2"
)
DECISION = (
    "1. If any gene misses the HARNESS under any form, the run stops before the RNA is opened"
    " (HARNESS_DEFECT). 2. Per form: a population of fewer than MIN_GENES genes is NOT_TESTABLE (no"
    " comparison read for it). 3. Per form, on Delta MSE(F) and its 97.5% bootstrap interval: lower bound"
    " > 0 -> BEATS (control p <= CONTROL_ALPHA -> BEATS_AND_CONTROL, else BEATS_NOT_CONTROL); upper bound"
    " < 0 -> WORSE; otherwise NO_DIFFERENCE. 4. Head-to-head, over the genes in both populations (if at"
    " least MIN_GENES), from the 95% interval of MSE(mean) - MSE(competition): lower bound > 0 -> 'closer"
    " to measured RNA than the mean', upper bound < 0 -> 'further from measured RNA than the mean',"
    " otherwise 'not detectably different from the mean'"
)
FALSIFIER = (
    "per form: the form's gated run does not beat the no-rule baseline (its 97.5% interval of Delta MSE"
    " includes 0 or lies below it): run under that form, the compiled rules of multi-element genes do"
    " not make a run track measured RNA across the reader's cell types in this population. The narrower"
    " claim that a gene's OWN elements' chromatin carries what it tracks is falsified when the"
    " shuffled-reader control is not beaten (control p > CONTROL_ALPHA). The claim that the strength-"
    "dependent form improves on the mean is falsified when the head-to-head is not 'closer than'"
)
STOP_RULE = (
    "one run, from this registration's committed code, on the files the registration pins by sha256 (the"
    " 24 compiled programs, the rate table, the gene table, the narrowPeak files, the census and the RNA"
    " quantifications). No form, parameter, population rule, endpoint, pseudocount, metric, control or"
    " wording is changed after an RNA value is read; a defect found afterwards is reported beside the"
    " registered reading and corrected only under a new registration"
)
FORM_WORDS = {
    "mean": "the runtime's mean",
    "competition": "the promoter-competition form (Bothma 2015, declared half-occupied promoter)",
}
#: the exact sentence of each per-form outcome; `reading` fills the numbers and the form's words only
OUTCOMES = {
    "HARNESS_DEFECT": (
        "stopped before the measured RNA was opened: the run differed from its closed form or missed an"
        " observation it was fitted to at {bad} of {n} gene-form runs"
    ),
    "NOT_TESTABLE": (
        "{form}: not testable at the registered size; {n} genes, fewer than the registered {min}; no"
        " comparison was read"
    ),
    "BEATS_AND_CONTROL": (
        "{form}: beats baseline; the reader-gated run of genes with two or more compiled elements lowered"
        " the squared error of measured log2(TPM+1) across the thirteen reader cell types below the no-rule"
        " baseline's (Delta MSE {d} [97.5% CI {lo}, {hi}] over {n} gene-cell pairs, {g} genes), by more"
        " than the shuffled-reader control (control p = {p})"
    ),
    "BEATS_NOT_CONTROL": (
        "{form}: beats baseline, but not the shuffled-reader control; the reader-gated run of genes with"
        " two or more compiled elements lowered the squared error of measured log2(TPM+1) across the"
        " thirteen reader cell types below the no-rule baseline's (Delta MSE {d} [97.5% CI {lo}, {hi}] over"
        " {n} gene-cell pairs, {g} genes), and other elements' chromatin did as well (control p = {p}); the"
        " gain is not carried by the genes' own elements"
    ),
    "NO_DIFFERENCE": (
        "{form}: no difference detected; the reader-gated run of genes with two or more compiled elements"
        " did not change the squared error of measured log2(TPM+1) across the thirteen reader cell types"
        " against the no-rule baseline by a detectable amount (Delta MSE {d} [97.5% CI {lo}, {hi}] over {n}"
        " gene-cell pairs, {g} genes; control p = {p})"
    ),
    "WORSE": (
        "{form}: worse; the reader-gated run of genes with two or more compiled elements raised the squared"
        " error of measured log2(TPM+1) across the thirteen reader cell types above the no-rule baseline's"
        " (Delta MSE {d} [97.5% CI {lo}, {hi}] over {n} gene-cell pairs, {g} genes; control p = {p})"
    ),
}
HEAD_TO_HEAD = {
    "closer": "closer to measured RNA than the mean",
    "further": "further from measured RNA than the mean",
    "same": "not detectably different from the mean",
    "sentence": (
        "Head-to-head over the {g} genes in both populations, the promoter-competition form was {versus}"
        " (MSE(mean) - MSE(competition) {u} [95% CI {ulo}, {uhi}])"
    ),
    "not_compared": (
        "Head-to-head not compared: {g} genes are in both populations, fewer than the registered {min}"
    ),
}
#: carried with every outcome, word for word. The first is the independence caveat, stated up front.
LIMITATIONS = (
    "independence: each element's removal fold comes from AlphaGenome's RNA-seq gene scorer (an in-silico"
    " deletion, the argmax over 371 RNA-seq tracks), a model trained on ENCODE and GTEx RNA-seq, with"
    " RNA-seq tracks named for the reader biosamples among them, and on ENCODE DNase-seq; the reader's"
    " open/closed calls are ENCODE DNase-seq narrowPeak files of the same biosamples; the measured"
    " outcome is ENCODE polyA plus RNA-seq of those biosamples, data of the kind, and possibly the very"
    " experiments, the model was trained on. The cross-cell pattern of a gated prediction comes from the"
    " DNase calls alone and the size of each element's part from the model",
    "the joint fit treats a gene's k observations as taken at one intact state, but each is the argmax"
    " over 371 tracks in its own track, and each reader cell's run transports them there. A maximum per"
    " element inflates the summed single-removal effects, which is one reason most multi-element genes"
    " are not runnable under either form (the census)",
    "the population is the genes whose summed single-removal effects leave a positive basal rate under the"
    " form: a selected subset of the multi-element genes, different for each form; a gene that is not"
    " runnable is not run, and nothing here says how those genes would read",
    "declared, not measured: the competition form's promoter occupancy THETA = 0.5 and one transcriptional"
    " efficiency per gene; the DECLARED_STRENGTH convention fixing the scale of s; the bridge's Hill"
    " assumptions (K = 1, n = 2, h = 0.5), which decide which inhibitors are feasible (s <= 1)",
    "accessibility near a gene tracks its expression across cell types whatever any compiled rule says;"
    " an open element is one a rule may act through, not evidence that it does, and no reading here is"
    " evidence of a mechanism",
    "one polyA plus experiment per biosample, from several labs; two of the thirteen are tissues (ovary,"
    " testis). Every rate under the run is mouse NIH 3T3 or their genome median; no fold prediction"
    " depends on it and no level in molecules per cell is quoted",
)
DESCRIPTIVES = (
    "per form: Delta MSE by number of rules (2-3, 4-9, 10 or more) and by whether the gene has an"
    " inhibiting rule; Delta MSE per reader cell",
    "per form: Spearman correlation of P_F with D over pairs with P_F != 0; per gene, the Spearman"
    " correlation of P_F with D across the thirteen cells, its median and the share of genes above 0",
    "accessibility only: Spearman correlation of D with the centred share of the gene's elements open in"
    " the cell (no rule sign, no size), per form population",
    "per form: the predicted range max - min of s_F over the thirteen cells, quantiles (predictor-side)",
    "context, not computed here: M2b's registered reading for the single-rule genes, copied verbatim from"
    " data/results/readergate_result.json",
)


# ==== the joint fit (closed form; no RNA) ===============================================================


@dataclass
class Fit:
    """One gene's joint fit under one form, relative to its intact transcription T (all fractions of T)."""

    form: str
    ok: bool
    reason: str = ""
    q: float = math.nan  # element-driven share of intact transcription; basal fraction is 1 - q
    z: float = math.nan  # competition: V*R0/T; mean: 1
    weights: dict[str, float] = field(default_factory=dict)  # activator -> d_i (mean) or kappa_i
    e: dict[str, float] = field(default_factory=dict)  # inhibitor -> RHO_j - 1
    s: dict[str, float] = field(default_factory=dict)  # element -> fitted strength
    vmax_frac: float = math.nan  # max_rate / T
    half_saturation: float = math.nan  # competition: the gene's promoter_half_saturation


def _strength_ok(s: float) -> float | None:
    if 1.0 < s <= 1.0 + bridge.STRENGTH_TOLERANCE:
        return 1.0
    return s if 0.0 < s <= 1.0 else None


def solve_z(d: list[float], k: float) -> float:
    """The one Z > 0 with sum of d_i (1 + K)^2 / (Z + d_i (1 + K)) = K (FIT['competition'])."""
    from scipy.optimize import brentq

    def g(z: float) -> float:
        return sum(x * (1 + k) ** 2 / (z + x * (1 + k)) for x in d) - k

    hi = 2.0 * len(d) * max(d) * (1 + k) ** 2 / k + 1.0
    return float(brentq(g, 1e-300, hi, xtol=1e-300, rtol=1e-15, maxiter=500))


def fit(form: str, rho: dict[str, float], theta: float = THETA) -> Fit:
    """The joint fit of one gene's elements under `form`, from their removal folds RHO (element -> RHO)."""
    if form not in FORMS:
        raise ValueError(f"form {form!r} is not one of {FORMS}")
    d = {el: 1.0 - r for el, r in rho.items() if r < 1.0}
    e = {el: r - 1.0 for el, r in rho.items() if r > 1.0}
    if len(d) + len(e) != len(rho):
        return Fit(form, False, "observation_null_effect")
    if not d:
        return Fit(form, False, "no_activating_rule")
    if form == "mean":
        q, z = sum(d.values()), 1.0
        weights = dict(d)
    else:
        k = theta / (1.0 - theta)
        z = solve_z(list(d.values()), k)
        weights = {el: x * (1 + k) ** 2 / (z + x * (1 + k)) for el, x in d.items()}
        q = z * theta
    if not q < 1.0:
        return Fit(form, False, "basal_not_positive", q=q, z=z)
    s: dict[str, float] = {}
    r0 = 1.0
    for el, x in e.items():
        sj = _strength_ok(x / (H * (q + x)))
        if sj is None:
            return Fit(form, False, "inhibitor_strength_above_1", q=q, z=z)
        s[el] = sj
        r0 *= 1.0 - sj * H
    top = max(weights.values())
    for el, w in weights.items():
        s[el] = w / top
    if form == "mean":
        vmax, half = len(d) * top / (H * r0), math.nan
    else:
        vmax, half = z / r0, H / top
    return Fit(form, True, "", q, z, weights, e, s, vmax, half)


def level(f: Fit, open_: dict[str, bool]) -> float:
    """Transcription relative to T with the elements in `open_` open and the rest held (the closed form)."""
    w = sum(x for el, x in f.weights.items() if open_.get(el, False))
    a = w if f.form == "mean" else f.z * w / (1.0 + w)
    boost = 1.0
    for el, x in f.e.items():
        if not open_.get(el, False):
            boost *= 1.0 + x / f.q
    return (1.0 - f.q) + a * boost


# ==== reading the programs (no RNA) =========================================================================

_RULE = re.compile(
    r"^rule (\S+) (activates|inhibits) (\S+) \{ strength: ([^;]*); when: cell_type = ([^;]*); .*?"
    r"evidence: (predicted|experimental) \"(.*)$"
)


def rule_lines(text: str) -> list[tuple[str, str, str, str, str, str]]:
    """(element, action, gene, compiled cell, evidence kind, rest of line) for every rule line of an
    UNGATED program; refuses a rule line it cannot read."""
    out = []
    for ln in text.splitlines():
        if ln.startswith("rule "):
            m = _RULE.match(ln)
            if m is None:
                raise ValueError(f"a rule line this census cannot read: {ln[:120]}")
            el, act, gene, _st, cell, kind, rest = m.groups()
            out.append((el, act, gene, cell, kind, rest))
    return out


def multi_rule_genes(texts: dict[str, str]) -> tuple[dict[str, str], Counter]:
    """gene -> its chromosome, for genes targeted by two or more rule lines, all in one program and all
    predicted (POPULATION (a)); exclusions counted. Read from the ungated text."""
    seen: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for chrom, text in texts.items():
        for _el, _act, gene, _cell, kind, _rest in rule_lines(text):
            seen[gene].append((chrom, kind))
    out: dict[str, str] = {}
    why: Counter = Counter()
    for g, rows in seen.items():
        if len(rows) < 2:
            why["one_rule"] += 1
            continue
        why["two_or_more_rules"] += 1
        if len({c for c, _ in rows}) > 1:
            why["rules_in_several_programs"] += 1
            continue
        if any(k != "predicted" for _, k in rows):
            why["a_rule_not_predicted"] += 1
            continue
        out[g] = rows[0][0]
    why["genes_with_any_rule"] = len(seen)
    why["multi_rule_all_predicted_one_program"] = len(out)
    return out, why


@dataclass
class GeneRules:
    """One multi-element gene as the gated program states it."""

    gene: str
    symbol: str
    chrom: str
    elements: list[str]
    actions: list[str]
    rho: dict[str, float]
    unreadable: list[str]
    compiled: dict[str, str]  # element -> compiled cell label
    states: dict[str, list[str]]  # element -> engine reader state in each of CELLS

    def open_in(self, j: int) -> dict[str, bool]:
        return {el: self.states[el][j] == "open" for el in self.elements}


def gene_rules(
    module: Module, genes: set[str], chrom: str, compiled: dict[tuple[str, str], str]
) -> list[GeneRules]:
    """The multi-element genes of one gated, parsed program, with each rule's RHO and reader states."""
    by_gene: dict[str, list] = defaultdict(list)
    for r in module.rules:
        if r.action in REGULATORY and r.target in genes:
            by_gene[r.target].append(r)
    states = {cid: module.reader_states(cid) for cid in CELLS}
    out = []
    for g in sorted(by_gene):
        rules = by_gene[g]
        elements = [r.through for r in rules]
        own = all(el and el == r.source for el, r in zip(elements, rules, strict=True))
        if not own or len(set(elements)) != len(elements):
            raise ValueError(f"{chrom}: {g} has a rule that does not act through its own, distinct element")
        rho: dict[str, float] = {}
        bad = []
        for r in rules:
            ob = bridge.observation_of(r, module)
            if isinstance(ob, bridge.Observation):
                rho[r.through] = ob.rho
            else:
                bad.append(f"{r.through}: {ob.reason}")
        out.append(
            GeneRules(
                gene=g,
                symbol=str(getattr(module.entities.get(g), "symbol", "") or g),
                chrom=chrom,
                elements=elements,
                actions=["activates" if r.action is Action.ACTIVATE else "inhibits" for r in rules],
                rho=rho,
                unreadable=bad,
                compiled={r.through: compiled.get((r.source, g), "") for r in rules},
                states={el: [states[cid].get(el, ("unknown", None))[0] for cid in CELLS] for el in elements},
            )
        )
    return out


def predicted_log2(f: Fit, gr: GeneRules) -> list[float]:
    """log2 of the closed-form level relative to T in each of CELLS."""
    return [math.log2(level(f, gr.open_in(j))) for j in range(len(CELLS))]


def classify(gr: GeneRules) -> tuple[dict[str, Fit], str]:
    """Each form's fit for one gene, or the predictor-side reason it has none (POPULATION (b), (c))."""
    if gr.unreadable:
        return {}, "observation_unreadable"
    if any(r == 1.0 for r in gr.rho.values()):
        return {}, "observation_null_effect"
    sign = {el: ("activates" if gr.rho[el] < 1.0 else "inhibits") for el in gr.elements}
    if any(sign[el] != a for el, a in zip(gr.elements, gr.actions, strict=True)):
        return {}, "action_disagrees_with_observation"
    if all(r > 1.0 for r in gr.rho.values()):
        return {}, "no_activating_rule"
    return {form: fit(form, gr.rho) for form in FORMS}, ""


# ==== the census: the compiled contexts, no reader (same-track joint fit) =================================


def compiled_context_census(texts: dict[str, str]) -> dict[str, Any]:
    """Per compiled `when` context (the argmax track, no reader): gene-cell pairs with two or more
    mechanisms (bridge.mechanism_key: an element and its `_measured` twin are one), and how many each
    form's joint fit makes runnable there, where all of a pair's observations come from ONE track."""
    pairs: dict[tuple[str, str], dict[str, list[tuple[str, str, str]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for text in texts.values():
        for el, act, gene, cell, kind, rest in rule_lines(text):
            pairs[(gene, cell)][bridge.element_of(el)].append((act, kind, rest))
    tot: Counter = Counter()
    per: dict[str, Counter] = defaultdict(Counter)
    for (_gene, cell), mech in pairs.items():
        if len(mech) < 2:
            continue
        c = per[cell]
        tot["pairs_two_or_more_mechanisms"] += 1
        c["pairs_two_or_more_mechanisms"] += 1
        rho: dict[str, float] = {}
        why = ""
        for el, reps in mech.items():
            if any(k != "predicted" for _a, k, _r in reps):
                why = "a_measured_rule"
                break
            m = bridge._EFFECT.search(reps[0][2])
            if m is None:
                why = "observation_unreadable"
                break
            rho[el] = 2.0 ** float(m.group(1))
        if not why and any(r == 1.0 for r in rho.values()):
            why = "observation_null_effect"
        if not why and all(r > 1.0 for r in rho.values()):
            why = "no_activating_rule"
        if why:
            tot[f"excluded_{why}"] += 1
            continue
        for form in FORMS:
            f = fit(form, rho)
            tot[f"runnable_{form}" if f.ok else f"not_runnable_{form}_{f.reason}"] += 1
            if f.ok:
                c[f"runnable_{form}"] += 1
    return {
        "totals": dict(sorted(tot.items())),
        "contexts": len(per),
        "per_context": {
            cell: dict(c)
            for cell, c in sorted(per.items(), key=lambda kv: -kv[1]["pairs_two_or_more_mechanisms"])
        },
    }


# ==== running (the runtime; no RNA) =======================================================================


def fitted_module(small: Module, gene: str, f: Fit, rate: float, label: dict[str, str]) -> Module:
    """The gene's gated module with the fit written in: basal and max rate, every rule's strength and
    action, the bridge's Hill assumptions, the rate label, and K under competition. `through` is kept."""
    m = copy.deepcopy(small)
    g = m.entities[gene]
    g.basal_rate = float(rate * (1.0 - f.q))
    g.attrs["max_rate"] = float(rate * f.vmax_frac)
    if f.form == "competition":
        g.attrs[grn.PROMOTER_HALF_SATURATION] = float(f.half_saturation)
    for attr, k in zip(grn.RATE_PROVENANCE_ATTRS, ("source", "species_cell", "tier"), strict=True):
        g.attrs[attr] = label[k]
    for r in m.rules:
        if r.target != gene or r.action not in REGULATORY:
            continue
        r.strength = f.s[r.through]
        r.threshold = bridge.HILL_THRESHOLD
        r.hill = bridge.HILL_COEFFICIENT
        r.action = Action.ACTIVATE if r.through in f.weights else Action.INHIBIT
    return m


def _run_level(module: Module, gene: str, ctx: dict[str, str], clamp: dict[str, float], form: str) -> float:
    rt = NetworkRuntime(module, context=ctx, strict=RUNTIME["strict"], activator_combination=form)
    steps = int(round(RUNTIME["hours"] / RUNTIME["dt"]))
    traj = rt.run(hours=RUNTIME["hours"], dt=RUNTIME["dt"], clamp=clamp, record_every=steps)
    return traj.final()[f"{gene}.mRNA"]


def run_gene_form(small: Module, gr: GeneRules, f: Fit, rate: float, label: dict[str, str]) -> dict[str, Any]:
    """The gated run of one gene under one form in every reader cell, and the HARNESS checks."""
    m = fitted_module(small, gr.gene, f, rate, label)
    present = dict.fromkeys(gr.elements, bridge.INTACT)
    delta = math.log(2) / grn.DEFAULTS["mrna_half_life"]
    levels, closed, miss = [], [], []
    for j, cid in enumerate(CELLS):
        lv = _run_level(m, gr.gene, {"cell_type": cid}, present, f.form)
        want = rate * level(f, gr.open_in(j)) / delta
        levels.append(lv)
        closed.append(math.log2(want))
        if not abs(math.log2(lv / want)) <= CLOSED_FORM_TOLERANCE:
            miss.append(f"cell {cid}: run {lv:.12g}, closed form {want:.12g}")
    ungated = copy.deepcopy(m)
    for r in ungated.rules:
        r.through = ""
    ctx = {"cell_type": CELLS[0]}
    full = _run_level(ungated, gr.gene, ctx, present, f.form)
    for el in gr.elements:
        lv = _run_level(ungated, gr.gene, ctx, {**present, el: bridge.REMOVED}, f.form)
        if not abs(math.log2(lv / full) - math.log2(gr.rho[el])) <= CLOSED_FORM_TOLERANCE:
            miss.append(f"removal {el}: run {lv / full:.12g}, observed {gr.rho[el]:.12g}")
    return {"levels": levels, "s": [math.log2(x) for x in levels], "closed_s": closed, "harness_miss": miss}


# ==== the score (RNA) =====================================================================================


def assemble(
    runs: list[dict[str, Any]], tpm: dict[str, dict[str, float]], ids: dict[str, str]
) -> tuple[list[dict], Counter]:
    """Rows for one form's population: D and the centred prediction over CELLS; exclusions by reason.
    `runs` are that form's informative gene runs; applies POPULATION (f)."""
    out, why = [], Counter()
    for g in runs:
        gid = ids.get(g["symbol"])
        if gid is None:
            why["symbol_not_on_exactly_one_gene_id"] += 1
            continue
        if any(gid not in tpm[c] for c in CELLS):
            why["gene_id_not_on_exactly_one_row_of_every_table"] += 1
            continue
        y = [rg._y(tpm[c][gid]) for c in CELLS]
        out.append(
            {
                **{
                    k: g[k]
                    for k in ("gene", "symbol", "chrom", "n_rules", "n_inhibiting", "elements", "states")
                },
                **{k: g[k] for k in ("weights", "e", "q", "z")},
                "gene_id": gid,
                "D": rg._centre(y),
                "P": rg._centre(g["s"]),
                "open_share": rg._centre(g["open_share"]),
            }
        )
    return out, why


def _gene_sums(rows: list[dict], key: str = "P") -> list[float]:
    return [sum(2.0 * p * d - p * p for p, d in zip(r[key], r["D"], strict=True)) for r in rows]


def dmse(rows: list[dict], key: str = "P") -> float:
    """Delta MSE of one prediction over every pair of the rows (mean of 2 P D - P^2)."""
    n = sum(len(r["D"]) for r in rows)
    return sum(_gene_sums(rows, key)) / n if n else float("nan")


def bootstrap(rows: list[dict], resamples: int, seed: int, pct: tuple[float, float]) -> tuple[float, float]:
    """Percentile interval of Delta MSE over gene resamples (each gene carries its thirteen pairs)."""
    import numpy as np

    s = np.array(_gene_sums(rows))
    n = np.array([len(r["D"]) for r in rows], dtype=float)
    rng = np.random.default_rng(seed)
    stats = np.empty(resamples)
    for i in range(resamples):
        idx = rng.integers(0, len(rows), len(rows))
        stats[i] = s[idx].sum() / n[idx].sum()
    lo, hi = np.percentile(stats, list(pct))
    return float(lo), float(hi)


def head_to_head(mean_rows: dict[str, dict], comp_rows: dict[str, dict], resamples: int, seed: int) -> dict:
    """MSE(mean) - MSE(competition) over the genes in both populations, with its 95% interval."""
    import numpy as np

    common = sorted(set(mean_rows) & set(comp_rows))
    if len(common) < MIN_GENES:
        return {"genes": len(common), "compared": False}
    a = np.array(_gene_sums([mean_rows[g] for g in common]))
    b = np.array(_gene_sums([comp_rows[g] for g in common]))
    n = np.array([len(mean_rows[g]["D"]) for g in common], dtype=float)
    diff = b - a  # per gene: SSE(mean) - SSE(competition)
    rng = np.random.default_rng(seed)
    stats = np.empty(resamples)
    for i in range(resamples):
        idx = rng.integers(0, len(common), len(common))
        stats[i] = diff[idx].sum() / n[idx].sum()
    lo, hi = np.percentile(stats, list(HEAD_TO_HEAD_INTERVAL))
    return {
        "genes": len(common),
        "pairs": int(n.sum()),
        "compared": True,
        "mse_mean_minus_mse_competition": float(diff.sum() / n.sum()),
        "ci95": [float(lo), float(hi)],
        "clause": HEAD_TO_HEAD["closer" if lo > 0 else "further" if hi < 0 else "same"],
    }


def _closed_arrays(rows: list[dict], form: str) -> dict[str, Any]:
    """Per-element arrays, grouped by gene, for the vectorised closed form under permuted open rows."""
    import numpy as np

    w, le, open_, starts, q, z, d = [], [], [], [], [], [], []
    for r in rows:
        starts.append(len(w))
        for el in r["elements"]:
            w.append(r["weights"].get(el, 0.0))
            le.append(math.log1p(r["e"][el] / r["q"]) if el in r["e"] else 0.0)
            open_.append([s == "open" for s in r["states"][el]])
        q.append(r["q"])
        z.append(r["z"])
        d.append(r["D"])
    return {
        "w": np.array(w)[:, None],
        "le": np.array(le)[:, None],
        "open": np.array(open_, dtype=float),
        "starts": np.array(starts),
        "q": np.array(q)[:, None],
        "z": np.array(z)[:, None],
        "D": np.array(d),
        "competition": form == "competition",
    }


def closed_dmse(arr: dict[str, Any], open_rows: Any) -> float:
    """Delta MSE of the closed form with the given element open rows (elements x cells)."""
    import numpy as np

    wsum = np.add.reduceat(arr["w"] * open_rows, arr["starts"], axis=0)
    lsum = np.add.reduceat(arr["le"] * (1.0 - open_rows), arr["starts"], axis=0)
    a = arr["z"] * wsum / (1.0 + wsum) if arr["competition"] else wsum
    s = np.log2((1.0 - arr["q"]) + a * np.exp(lsum))
    p = s - s.mean(axis=1, keepdims=True)
    d = arr["D"]
    return float((2.0 * p * d - p * p).sum() / d.size)


def shuffled_reader(rows: list[dict], form: str, count: int, seed: int) -> tuple[list[float], float]:
    """Delta MSE(F) with the element open rows permuted among all the population's elements, `count`
    times; and the unpermuted closed form's Delta MSE (which the run's must equal within the HARNESS)."""
    import numpy as np

    arr = _closed_arrays(rows, form)
    rng = np.random.default_rng(seed)
    base = closed_dmse(arr, arr["open"])
    n_el = arr["open"].shape[0]
    return [closed_dmse(arr, arr["open"][rng.permutation(n_el)]) for _ in range(count)], base


def decide(lo: float, hi: float, p: float) -> str:
    """DECISION, step 3."""
    if lo > 0:
        return "BEATS_AND_CONTROL" if p <= CONTROL_ALPHA else "BEATS_NOT_CONTROL"
    if hi < 0:
        return "WORSE"
    return "NO_DIFFERENCE"


def reading(kind: str, **numbers: Any) -> str:
    """The registered sentence for an outcome, numbers and the form's words filled, nothing else."""
    return OUTCOMES[kind].format(**numbers)


def score_form(rows: list[dict], form: str) -> dict[str, Any]:
    """The registered comparison and descriptives for one form's assembled population."""
    genes = len(rows)
    words = FORM_WORDS[form]
    if genes < MIN_GENES:
        return {
            "outcome": "NOT_TESTABLE",
            "reading": reading("NOT_TESTABLE", form=words, n=genes, min=MIN_GENES),
            "genes": genes,
        }
    n = sum(len(r["D"]) for r in rows)
    dm = dmse(rows)
    lo, hi = bootstrap(rows, BOOTSTRAP["resamples"], BOOTSTRAP["seed"], FORM_INTERVAL)
    perm, closed = shuffled_reader(rows, form, PERMUTATIONS["count"], PERMUTATIONS["seed"])
    p_ctrl = (1 + sum(1 for x in perm if x >= dm)) / (1 + len(perm))
    kind = decide(lo, hi, p_ctrl)
    nums = {
        "form": words,
        "d": f"{dm:+.4g}",
        "lo": f"{lo:+.4g}",
        "hi": f"{hi:+.4g}",
        "n": f"{n:,}",
        "g": f"{genes:,}",
        "p": f"{p_ctrl:.3g}",
    }
    d_all = [x for r in rows for x in r["D"]]
    p_all = [x for r in rows for x in r["P"]]
    nz = [(p, d) for p, d in zip(p_all, d_all, strict=True) if p != 0]
    per_gene = [x for x in (rr._spearman(r["P"], r["D"]) for r in rows) if x is not None]

    def sub(rs: list[dict]) -> dict[str, Any]:
        return {"genes": len(rs), "delta_mse": dmse(rs) if rs else None}

    by_cell = {}
    for j, c in enumerate(CELLS):
        pairs = [(r["P"][j], r["D"][j]) for r in rows]
        by_cell[c] = sum(2 * p * d - p * p for p, d in pairs) / len(pairs)
    ranges = sorted(max(r["P"]) - min(r["P"]) for r in rows)
    return {
        "outcome": kind,
        "reading": reading(kind, **nums),
        "genes": genes,
        "pairs": n,
        "delta_mse": dm,
        "ci97_5": [lo, hi],
        "closed_form_delta_mse": closed,
        "mse_baseline": sum(x * x for x in d_all) / n,
        "mse_form": sum((d - p) ** 2 for p, d in zip(p_all, d_all, strict=True)) / n,
        "control_p": p_ctrl,
        "control_delta_mse_median": rr._median(perm),
        "control_delta_mse_95th_percentile": sorted(perm)[int(0.95 * len(perm)) - 1],
        "descriptives": {
            "by_rule_count": {
                "2-3": sub([r for r in rows if r["n_rules"] <= 3]),
                "4-9": sub([r for r in rows if 4 <= r["n_rules"] <= 9]),
                "10+": sub([r for r in rows if r["n_rules"] >= 10]),
            },
            "by_inhibitor": {
                "has_inhibiting_rule": sub([r for r in rows if r["n_inhibiting"] > 0]),
                "activating_only": sub([r for r in rows if r["n_inhibiting"] == 0]),
            },
            "by_cell": by_cell,
            "spearman_P_D_nonzero": rr._spearman([p for p, _ in nz], [d for _, d in nz]),
            "per_gene_spearman": {
                "genes_defined": len(per_gene),
                "median": rr._median(per_gene) if per_gene else None,
                "share_above_0": sum(1 for x in per_gene if x > 0) / len(per_gene) if per_gene else None,
            },
            "accessibility_only_spearman_open_share_D": rr._spearman(
                [x for r in rows for x in r["open_share"]], d_all
            ),
            "predicted_range_log2_quantiles": {
                str(qq): ranges[min(len(ranges) - 1, int(qq * len(ranges)))]
                for qq in (0.1, 0.25, 0.5, 0.75, 0.9)
            },
        },
    }


def combined_reading(per_form: dict[str, dict], h2h: dict[str, Any]) -> str:
    """The registered reading: each form's sentence, then the head-to-head sentence."""
    parts = [per_form[f]["reading"] for f in FORMS]
    if h2h.get("compared"):
        tail = HEAD_TO_HEAD["sentence"].format(
            g=f"{h2h['genes']:,}",
            versus=h2h["clause"],
            u=f"{h2h['mse_mean_minus_mse_competition']:+.4g}",
            ulo=f"{h2h['ci95'][0]:+.4g}",
            uhi=f"{h2h['ci95'][1]:+.4g}",
        )
    else:
        tail = HEAD_TO_HEAD["not_compared"].format(g=f"{h2h.get('genes', 0):,}", min=MIN_GENES)
    return "; ".join(parts) + ". " + tail + "."
