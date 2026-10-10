# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""M2f (ROADMAP §5 item 7h, lane-multrna, 2026-10-10): the multiplicative composition against measured RNA.

7d (lane-combine2, `combine2.py`; registration 917ce25, result 948cfc9) ran the genes the compiled
programs drive through two or more elements under the runtime's mean and the promoter-competition form,
reader-gated in the thirteen reader cells, against ENCODE polyA plus RNA-seq of those cells: both
"worse" than a no-rule baseline, and competition further from the RNA than the mean. Most multi-element
genes were outside both fits (`basal_not_positive`: a gene's single-removal drops add up past its
whole level). 7f (lane-saturate, `saturate.py`; runtime rule 1f8f12b, registration fcdd4f2, result
416aaa5) registered the multiplicative composition of Zhou et al. 2024 (A = exp(G*W), basal 0, no
parameter of its own), which fits such genes, and read it only against held-out CRISPRi pairs in K562,
where every fit reproduces every single removal and so cannot tell one composition from another
(saturate.BY_CONSTRUCTION (4)). It has not been read against an expression measurement.

This module is that reading, on 7d's design exactly: 7d's population rules with the multiplicative
fit in place of 7d's, 7d's endpoint, baseline, metric, bootstrap, control, decision and sentences, and
a head-to-head against 7d's two forms on the genes all three run. It re-implements no stage: the gene
reading and exclusions are combine2's (`gene_rules`, `classify`, `assemble`), the fit and the fitted
module are saturate's (`fit`, `level`, `fitted_module`), the runtime is runtime/grn.py's opt-in rule.
Its own code is the run under the form, the multiplicative closed form of the control, the
head-to-heads and the readings. Nothing is fitted: the form has no parameter.

In the order written and committed: this code, the census (no RNA read), the registration, then one
run. No model request, no key, no money.
"""

from __future__ import annotations

import copy
import json
import math
from typing import Any

from genomeos.attribution import bridge
from genomeos.attribution import combine2 as c2
from genomeos.attribution import runrules as rr
from genomeos.attribution import saturate as sa
from genomeos.ir import Module
from genomeos.runtime import grn

REGISTERED = "2026-10-10"
LANE = "lane-multrna"
FORM = sa.FORM
COMPARATOR_FORMS = c2.FORMS
SELECTION = c2.SELECTION
CELLS = c2.CELLS
EXPONENT_LIMIT = sa.EXPONENT_LIMIT

# ---- the form, as 7f registered it ---------------------------------------------------------------------
FORM_SOURCE = sa.FORM_SOURCE
FORM_WORDS = {FORM: sa.FORM_WORDS[FORM]}
FIT = {
    "principle": c2.FIT["principle"],
    FORM: sa.FIT[FORM],
    "runnable": sa.FIT["runnable"],
    "nothing_fitted": (
        "the form has no parameter of its own and nothing is fitted on any RNA value: a gene's numbers"
        " are its elements' removal observations RHO_i and its rate T, on which no fold depends"
    ),
}
STRUCTURE = {
    "closed_form_under_the_hold": (
        "under the multiplicative rule a held activator (the reader's strength-0 copy) contributes exp(0)"
        " = 1 where it contributed exp(G*s_i*H) = 1/RHO_i, and a held inhibitor contributes 1 where it"
        " contributed 1 - s_j*H = 1/RHO_j; the bridge's removal (a clamp at 0) does the same. So the level"
        " in a reader cell with the open set O is T times the product of RHO_i over the elements not in O,"
        " and the centred prediction is the centred sum, over the elements held in that cell, of log2"
        " RHO_i. No basal share, no T and no fitted scale enter it"
    ),
    "what_differs_from_7d": (
        "under the runtime's mean the drops of activators held together add on the level (each lowers it"
        " by its own d_i = 1 - RHO_i of T, and q = sum of d_i must stay below 1); under the multiplicative"
        " form they add on its logarithm (each multiplies it by RHO_i), so no sum of drops is bounded. With"
        " every activator held the multiplicative level is T times the product of their RHO_i, at or above"
        " the mean's T(1 - q) wherever both run (the product of 1 - d_i is at least 1 - sum of d_i): the"
        " multiplicative form runs genes whose drops sum past 1, and on a gene both run its prediction falls"
        " less where several activators are held together. One activator held with every other element"
        " open gives its own RHO_i under both"
    ),
    "what_is_assumed": c2.STRUCTURE["what_is_assumed"],
    "inhibitor_only": c2.STRUCTURE["inhibitor_only"],
    "network": c2.STRUCTURE["network"],
}

# ---- the registered test (written before any RNA value was read) --------------------------------------
QUESTION = (
    "For genes the compiled programs drive through two or more elements, does a reader-gated run under the"
    " multiplicative composition (Zhou et al. 2024: A = exp(G*W), basal 0, no parameter of its own) track"
    " measured RNA across the thirteen reader cell types better than a no-rule baseline, as 7d asked of the"
    " runtime's mean and the promoter-competition form; and on the genes all three forms run, is it closer"
    " to measured RNA than either of them?"
)
POPULATION = (
    "every gene G such that, over all 24 compiled programs (data/knowledge/compiled, pinned by sha256)"
    " gated with compile.gate_program(selection='reader_cells'): (a) two or more activates/inhibits rule"
    " lines target G, all in one program, all of predicted evidence, each acting through its own element;"
    " (b) every rule's removal observation is readable (bridge.observation_of) and none is exactly 1; (c)"
    " at least one rule is activating, and each rule's action agrees with its observation's sign"
    " (combine2.classify); (d) the multiplicative fit is runnable (FIT['runnable']: every inhibitor's"
    " RHO_j <= 2 and the summed log-folds within EXPONENT_LIMIT); (e) the gated prediction under the"
    " multiplicative form is not the same in all thirteen reader cells (its range exceeds INFORMATIVE_LOG2"
    " in log2); (f) G's symbol names exactly one gene_id in data/cache/gencode_genes.tsv and that gene_id"
    " is on exactly one row of every one of the thirteen quantification tables. (a)-(c) and (f) are"
    " combine2.POPULATION's, unchanged; (d) and (e) are the multiplicative form's in place of 7d's."
    " (a)-(e) are predictor-side: no measured value enters them. Each G contributes its thirteen (G, cell)"
    " pairs. No filter on expression. 7d's two populations are combine2.POPULATION's, recomputed by 7d's"
    " committed code (COMPARATOR)"
)
COMMON = (
    "the genes in all three populations (7d's mean, 7d's competition, the multiplicative form's), each with"
    " its thirteen pairs; a gene enters every form with the one D combine2.assemble gives it"
)
INFORMATIVE_LOG2 = c2.INFORMATIVE_LOG2
MIN_GENES = c2.MIN_GENES
PSEUDOCOUNT = c2.PSEUDOCOUNT
MEASURED = c2.MEASURED
RNA_SELECTION = c2.RNA_SELECTION
PREDICTIONS = {
    "form": (
        "s_M(G, X) = log2(L_M(G, X)), L_M(G, X) the level the runtime gives G in reader cell X under the"
        " multiplicative form (NetworkRuntime(..., activator_combination='multiplicative'), the fit written"
        " in by saturate.fitted_module, the fitted rules acting through their elements, the reader blocks in"
        " force, every element clamped present so that only the reader holds a rule). P_M(G, X) = s_M(G, X)"
        " - mean over the thirteen cells of s_M(G, X')"
    ),
    "comparators": (
        "7d's P_mean and P_competition, recomputed in this run by 7d's committed code"
        " (scripts/combine2_census.collect, scripts/combine2_run.form_runs, combine2.assemble), as 7d's"
        " PREDICTIONS['form'] defines them"
    ),
    "no_rule_baseline": c2.PREDICTIONS["no_rule_baseline"],
}
METRIC = (
    "Delta MSE(M) = mean over the population's pairs of D^2 - mean of (D - P_M)^2 (combine2.METRIC),"
    " positive when the form's squared error is below the baseline's; an uninformative prediction still"
    " costs mean(P^2), so it reads WORSE once the population resolves that cost. The head-to-heads, over"
    " the genes in all three populations (COMMON): MSE(mean) - MSE(multiplicative) and MSE(competition) -"
    " MSE(multiplicative), each positive when the multiplicative form is closer"
)
#: 7d's per-form interval (Bonferroni over its two forms), kept so the three forms stand on one interval
FORM_INTERVAL = c2.FORM_INTERVAL
#: two head-to-heads are judged where 7d judged one, so each is Bonferroni-adjusted: 97.5%
HEAD_TO_HEAD_INTERVAL = (1.25, 98.75)
BOOTSTRAP = {
    "resamples": c2.BOOTSTRAP["resamples"],
    "seed": c2.BOOTSTRAP["seed"],
    "unit": c2.BOOTSTRAP["unit"],
    "interval": (
        "the multiplicative form's 97.5% percentile interval, 7d's per-form interval (Bonferroni over 7d's"
        " two forms), kept so that the three forms' readings stand on one interval; each of the two"
        " head-to-heads 97.5% (Bonferroni over the two; 7d's one head-to-head was at 95%)"
    ),
}
CONTROL = (
    "shuffled-reader control, 7d's: the thirteen-cell open/held rows of the population's elements are"
    " permuted among all its elements (PERMUTATIONS times); every element keeps its own observation and"
    " role and every gene its fit, so the control keeps each cell's open count and the distribution of"
    " rows and breaks only the tie between a gene's elements and their own chromatin. The permuted"
    " prediction is the form's closed form (STRUCTURE['closed_form_under_the_hold']), which the run must"
    " equal (HARNESS). control p = (1 + permutations with Delta MSE(M) >= observed) / (1 + PERMUTATIONS)"
)
PERMUTATIONS = c2.PERMUTATIONS
CONTROL_ALPHA = c2.CONTROL_ALPHA
RUNTIME = c2.RUNTIME
CLOSED_FORM_TOLERANCE = c2.CLOSED_FORM_TOLERANCE
HARNESS = (
    "for every gene run under the multiplicative form: (1) in each of the thirteen reader cells the gated"
    " run's level equals T times the closed form for that cell's open set; (2) the same fitted program with"
    " its rules ungated, every element present, then each element clamped absent in turn, reproduces that"
    " element's RHO. Each within CLOSED_FORM_TOLERANCE in log2. 7d's two forms are run in the same run by"
    " 7d's code under 7d's HARNESS"
)
#: the committed genes tables round P and D to 5 decimals
COMPARATOR_TOLERANCE = 1e-5
COMPARATOR = (
    "7d's two forms enter only through the head-to-heads, recomputed in this run by 7d's committed code"
    " (scripts/combine2_census.collect, scripts/combine2_run.form_runs, combine2.assemble,"
    " combine2.score_form, combine2.head_to_head, combine2.combined_reading) on the same pinned inputs."
    " The run stops (COMPARATOR_DEFECT), before the multiplicative form's comparison is computed, unless:"
    " each of 7d's two populations is, gene for gene and in order, the committed one"
    " (data/results/combine2_result.json, genes_tables); every committed P and D is reproduced within"
    " COMPARATOR_TOLERANCE; and 7d's registered reading, recomputed, equals the committed reading word for"
    " word. A gene in all three populations must carry the same D in each"
)
CENSUS = (
    "a count on the predictor side, before any RNA value is read (data/results/multrna_census.json): of"
    " the multi-element genes combine2.classify admits, how many the multiplicative fit makes runnable"
    " (FIT['runnable']) and informative, and why the rest are not; beside them 7d's counts under the"
    " runtime's mean and promoter competition, recomputed by combine2.fit and combine2.predicted_log2 and"
    " checked against data/results/combine2_census.json; which of the multiplicative form's informative"
    " genes neither of 7d's forms could run (runnable and informative), by 7d's reasons; and how many are"
    " informative under all three forms. It reads no RNA value and decides nothing"
)
CENSUS_SENTENCE = (
    "of the {n} multi-element genes combine2.classify admits, the multiplicative composition makes {r}"
    " runnable ({x} are not: {reasons}) and {i} informative across the thirteen reader cells; 7d's forms"
    " made {m} (the runtime's mean) and {c} (promoter competition) informative, and {new} of the"
    " multiplicative form's informative genes are informative under neither 7d form (under the runtime's"
    " mean: {new_reasons}); {common} genes are informative under all three forms, before the RNA tables'"
    " row check"
)
DECISION = (
    "1. If any gene misses the HARNESS under the multiplicative form, or under either of 7d's forms in this"
    " run, the run stops before the RNA is opened (HARNESS_DEFECT). 2. After the RNA is opened and the"
    " populations assembled, if 7d's comparators do not reproduce (COMPARATOR), the run stops before the"
    " multiplicative form's comparison is computed (COMPARATOR_DEFECT). 3. A multiplicative population of"
    " fewer than MIN_GENES genes is NOT_TESTABLE (no comparison read for it). 4. On Delta MSE(M) and its"
    " 97.5% bootstrap interval: lower bound > 0 -> BEATS (control p <= CONTROL_ALPHA -> BEATS_AND_CONTROL,"
    " else BEATS_NOT_CONTROL); upper bound < 0 -> WORSE; otherwise NO_DIFFERENCE (7d's DECISION, step 3)."
    " 5. Head-to-heads, over the genes in all three populations (if at least MIN_GENES), each from its"
    " 97.5% interval of MSE(F) - MSE(multiplicative): lower bound > 0 -> 'closer to measured RNA than F',"
    " upper bound < 0 -> 'further from measured RNA than F', otherwise 'not detectably different from F'"
)
FALSIFIER = (
    "the claim that, run under the multiplicative composition, the compiled rules of multi-element genes"
    " make a reader-gated run track measured RNA across the reader's cell types in this population is"
    " falsified when the run does not beat the no-rule baseline (its 97.5% interval of Delta MSE includes 0"
    " or lies below it). The narrower claim that a gene's OWN elements' chromatin carries what it tracks is"
    " falsified when the shuffled-reader control is not beaten (control p > CONTROL_ALPHA). The claim that"
    " the multiplicative composition improves on 7d's mean (on 7d's competition form) on the genes all"
    " three run is falsified when that head-to-head is not 'closer than'"
)
STOP_RULE = (
    "one run, from this registration's committed code, on the files the registration pins by sha256 (the"
    " 24 compiled programs, the rate table, the gene table, the narrowPeak files, this lane's census, 7d's"
    " census, registration and result, the RNA source record and the RNA quantifications). No form,"
    " population rule, endpoint, pseudocount, metric, interval, control, comparator or wording is changed"
    " after an RNA value is read; a defect found afterwards is reported beside the registered reading and"
    " corrected only under a new registration"
)
#: 7d's per-form sentences, word for word, and the one stop this design adds
OUTCOMES = {
    **c2.OUTCOMES,
    "COMPARATOR_DEFECT": (
        "stopped before the multiplicative form's comparison was computed: 7d's comparators, recomputed in"
        " this run, did not reproduce the committed result ({what})"
    ),
}
HEAD_TO_HEAD = {
    "closer": "closer to measured RNA than {versus}",
    "further": "further from measured RNA than {versus}",
    "same": "not detectably different from {versus}",
    "sentence": (
        "Head-to-head over the {g} genes in all three populations ({n} gene-cell pairs), the multiplicative"
        " composition was {versus_mean} (MSE(mean) - MSE(multiplicative) {u} [97.5% CI {ulo}, {uhi}]) and"
        " {versus_competition} (MSE(competition) - MSE(multiplicative) {v} [97.5% CI {vlo}, {vhi}])"
    ),
    "not_compared": (
        "Head-to-head not compared: {g} genes are in all three populations, fewer than the registered {min}"
    ),
}
#: carried with every outcome, word for word. The first is 7d's independence caveat, stated up front.
LIMITATIONS = (
    c2.LIMITATIONS[0],
    "the fit treats a gene's k observations as taken at one intact state, but each is the argmax over 371"
    " tracks in its own track, and each reader cell's run transports them there. A maximum per element"
    " inflates each element's fold, and under the multiplicative form every held element's fold enters the"
    " prediction in full, so the inflation widens the predicted range",
    "the population is the genes the multiplicative fit runs, a different and larger set than either of"
    " 7d's (the census): the comparison against the baseline is not on 7d's genes, and only the"
    " head-to-heads are. The head-to-heads are on the genes all three forms run, those whose summed"
    " activating drops leave a positive basal rate under both of 7d's fits; nothing in them says how 7d's"
    " forms would read on the genes only the multiplicative form runs",
    "the multiplicative form is Zhou et al. 2024's reading of enhancers perturbed together in one K562"
    " CRISPRi screen (Gasperini 2019); here it composes every compiled element of a gene, up to dozens, in"
    " thirteen cell types, which that evidence does not cover. It sets the basal rate to 0: a gene's level"
    " with every element held is T*exp(-sum beta)/R0, and no level in molecules per cell is quoted",
    "declared, not measured: the DECLARED_STRENGTH convention and the bridge's Hill assumptions (K = 1,"
    " n = 2, h = 0.5), which under this form decide only which inhibitors are feasible (RHO_j <= 2); no"
    " fold prediction depends on them",
    c2.LIMITATIONS[4],
    c2.LIMITATIONS[5],
    "7d's two forms are read here only through the head-to-heads; their own readings against the baseline"
    " are 7d's registered readings, quoted verbatim and not re-decided",
)
DESCRIPTIVES = (
    "the multiplicative form, as 7d's per form: Delta MSE by number of rules (2-3, 4-9, 10 or more) and by"
    " whether the gene has an inhibiting rule; Delta MSE per reader cell; Spearman correlation of P_M with"
    " D over pairs with P_M != 0; per gene, the Spearman correlation of P_M with D across the thirteen"
    " cells, its median and the share of genes above 0; accessibility only (Spearman of D with the centred"
    " share of the gene's elements open); the predicted range quantiles (predictor-side)",
    "the multiplicative form's Delta MSE on the genes new under it (informative under neither 7d form) and"
    " on the genes of 7d's mean population, with the mean's own Delta MSE on the latter",
    "on the genes in all three populations: each form's Delta MSE and predicted range quantiles",
    "context, not decided here: 7d's registered reading and M2e's (data/results/saturate_result.json)"
    " census sentence, copied verbatim",
)
#: the registration fields the run refuses to differ from: registration field -> constant of this module
BOUND = {
    "question": "QUESTION",
    "form": "FORM",
    "comparator_forms": "COMPARATOR_FORMS",
    "form_source": "FORM_SOURCE",
    "form_words": "FORM_WORDS",
    "fit": "FIT",
    "exponent_limit": "EXPONENT_LIMIT",
    "structure": "STRUCTURE",
    "selection": "SELECTION",
    "cells": "CELLS",
    "rna_selection": "RNA_SELECTION",
    "population": "POPULATION",
    "common": "COMMON",
    "informative_log2": "INFORMATIVE_LOG2",
    "min_genes": "MIN_GENES",
    "pseudocount": "PSEUDOCOUNT",
    "measured": "MEASURED",
    "predictions": "PREDICTIONS",
    "metric": "METRIC",
    "form_interval": "FORM_INTERVAL",
    "head_to_head_interval": "HEAD_TO_HEAD_INTERVAL",
    "bootstrap": "BOOTSTRAP",
    "control": "CONTROL",
    "permutations": "PERMUTATIONS",
    "control_alpha": "CONTROL_ALPHA",
    "runtime": "RUNTIME",
    "closed_form_tolerance": "CLOSED_FORM_TOLERANCE",
    "harness": "HARNESS",
    "comparator": "COMPARATOR",
    "comparator_tolerance": "COMPARATOR_TOLERANCE",
    "census": "CENSUS",
    "census_sentence": "CENSUS_SENTENCE",
    "decision": "DECISION",
    "falsifier": "FALSIFIER",
    "stop_rule": "STOP_RULE",
    "outcomes": "OUTCOMES",
    "head_to_head": "HEAD_TO_HEAD",
    "limitations": "LIMITATIONS",
    "descriptives_not_verdicts": "DESCRIPTIVES",
}


def registered_fields() -> dict[str, Any]:
    """Every BOUND constant as JSON stores it, so a registration and this module compare equal or not."""
    return {f: json.loads(json.dumps(globals()[c])) for f, c in BOUND.items()}


# ==== the form on one gene (closed form and the runtime; no RNA) =========================================


def predicted_log2(f: sa.MFit, gr: c2.GeneRules) -> list[float]:
    """log2 of the closed-form level relative to T in each of CELLS (STRUCTURE)."""
    return [math.log2(sa.level(f, gr.open_in(j))) for j in range(len(CELLS))]


def run_gene(
    small: Module, gr: c2.GeneRules, f: sa.MFit, rate: float, label: dict[str, str]
) -> dict[str, Any]:
    """The gated run of one gene under the form in every reader cell, and the HARNESS checks
    (combine2.run_gene_form's, with saturate's fit and fitted module)."""
    m = sa.fitted_module(small, gr.gene, f, rate, label)
    present = dict.fromkeys(gr.elements, bridge.INTACT)
    delta = math.log(2) / grn.DEFAULTS["mrna_half_life"]
    levels, closed, miss = [], [], []
    for j, cid in enumerate(CELLS):
        lv = c2._run_level(m, gr.gene, {"cell_type": cid}, present, FORM)
        want = rate * sa.level(f, gr.open_in(j)) / delta
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


# ==== the score (RNA) =====================================================================================


def assemble(
    runs: list[dict[str, Any]], tpm: dict[str, dict[str, float]], ids: dict[str, str]
) -> tuple[list[dict], Any]:
    """The multiplicative population's rows: combine2.assemble's D, P and POPULATION (f), with each
    element's log2 RHO attached for the control."""
    recs = [{**r, "weights": {}, "e": {}, "q": math.nan, "z": math.nan} for r in runs]
    rows, why = c2.assemble(recs, tpm, ids)
    lr = {r["gene"]: r["log2_rho"] for r in runs}
    for row in rows:
        for k in ("weights", "e", "q", "z"):
            row.pop(k)
        row["log2_rho"] = lr[row["gene"]]
    return rows, why


def _closed_arrays(rows: list[dict]) -> dict[str, Any]:
    """Per-element log2 RHO and open rows, grouped by gene, for the vectorised closed form."""
    import numpy as np

    lr, open_, starts, d = [], [], [], []
    for r in rows:
        starts.append(len(lr))
        for el in r["elements"]:
            lr.append(r["log2_rho"][el])
            open_.append([s == "open" for s in r["states"][el]])
        d.append(r["D"])
    return {
        "lr": np.array(lr)[:, None],
        "open": np.array(open_, dtype=float),
        "starts": np.array(starts),
        "D": np.array(d),
    }


def closed_dmse(arr: dict[str, Any], open_rows: Any) -> float:
    """Delta MSE of the multiplicative closed form with the given element open rows (elements x cells)."""
    import numpy as np

    s = np.add.reduceat(arr["lr"] * (1.0 - open_rows), arr["starts"], axis=0)
    p = s - s.mean(axis=1, keepdims=True)
    d = arr["D"]
    return float((2.0 * p * d - p * p).sum() / d.size)


def shuffled_reader(rows: list[dict], count: int, seed: int) -> tuple[list[float], float]:
    """Delta MSE(M) with the element open rows permuted among all the population's elements, `count`
    times (combine2.shuffled_reader's draw); and the unpermuted closed form's Delta MSE."""
    import numpy as np

    arr = _closed_arrays(rows)
    rng = np.random.default_rng(seed)
    base = closed_dmse(arr, arr["open"])
    n_el = arr["open"].shape[0]
    return [closed_dmse(arr, arr["open"][rng.permutation(n_el)]) for _ in range(count)], base


def reading(kind: str, **numbers: Any) -> str:
    """The registered sentence for an outcome, numbers and the form's words filled, nothing else."""
    return OUTCOMES[kind].format(**numbers)


def _sub(rows: list[dict]) -> dict[str, Any]:
    return {"genes": len(rows), "delta_mse": c2.dmse(rows) if rows else None}


def quantiles(xs: list[float]) -> dict[str, float]:
    s = sorted(xs)
    return {str(q): s[min(len(s) - 1, int(q * len(s)))] for q in (0.1, 0.25, 0.5, 0.75, 0.9)} if s else {}


def score_form(rows: list[dict]) -> dict[str, Any]:
    """The registered comparison and 7d's descriptives for the multiplicative population
    (combine2.score_form's, with the multiplicative closed form in the control)."""
    genes = len(rows)
    words = FORM_WORDS[FORM]
    if genes < MIN_GENES:
        return {
            "outcome": "NOT_TESTABLE",
            "reading": reading("NOT_TESTABLE", form=words, n=genes, min=MIN_GENES),
            "genes": genes,
        }
    n = sum(len(r["D"]) for r in rows)
    dm = c2.dmse(rows)
    lo, hi = c2.bootstrap(rows, BOOTSTRAP["resamples"], BOOTSTRAP["seed"], FORM_INTERVAL)
    perm, closed = shuffled_reader(rows, PERMUTATIONS["count"], PERMUTATIONS["seed"])
    p_ctrl = (1 + sum(1 for x in perm if x >= dm)) / (1 + len(perm))
    kind = c2.decide(lo, hi, p_ctrl)
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
    by_cell = {}
    for j, c in enumerate(CELLS):
        pairs = [(r["P"][j], r["D"][j]) for r in rows]
        by_cell[c] = sum(2 * p * d - p * p for p, d in pairs) / len(pairs)
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
                "2-3": _sub([r for r in rows if r["n_rules"] <= 3]),
                "4-9": _sub([r for r in rows if 4 <= r["n_rules"] <= 9]),
                "10+": _sub([r for r in rows if r["n_rules"] >= 10]),
            },
            "by_inhibitor": {
                "has_inhibiting_rule": _sub([r for r in rows if r["n_inhibiting"] > 0]),
                "activating_only": _sub([r for r in rows if r["n_inhibiting"] == 0]),
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
            "predicted_range_log2_quantiles": quantiles([max(r["P"]) - min(r["P"]) for r in rows]),
        },
    }


def head_to_heads(rows: dict[str, dict[str, dict]], resamples: int, seed: int) -> dict[str, Any]:
    """MSE(F) - MSE(multiplicative) for each of 7d's forms over the genes in all three populations, each
    with its 97.5% interval (genes resampled with their thirteen pairs). `rows`: form -> gene -> row."""
    import numpy as np

    common = sorted(set.intersection(*(set(rows[f]) for f in (FORM, *COMPARATOR_FORMS))))
    if len(common) < MIN_GENES:
        return {"genes": len(common), "compared": False}
    for g in common:
        if any(rows[f][g]["D"] != rows[FORM][g]["D"] for f in COMPARATOR_FORMS):
            raise ValueError(f"{g} enters the forms with different D")
    gm = np.array(c2._gene_sums([rows[FORM][g] for g in common]))
    n = np.array([len(rows[FORM][g]["D"]) for g in common], dtype=float)
    out: dict[str, Any] = {"genes": len(common), "pairs": int(n.sum()), "compared": True}
    for f in COMPARATOR_FORMS:
        diff = gm - np.array(c2._gene_sums([rows[f][g] for g in common]))  # per gene: SSE(F) - SSE(M)
        rng = np.random.default_rng(seed)
        stats = np.empty(resamples)
        for i in range(resamples):
            idx = rng.integers(0, len(common), len(common))
            stats[i] = diff[idx].sum() / n[idx].sum()
        lo, hi = np.percentile(stats, list(HEAD_TO_HEAD_INTERVAL))
        kind = "closer" if lo > 0 else "further" if hi < 0 else "same"
        out[f] = {
            "mse_form_minus_mse_multiplicative": float(diff.sum() / n.sum()),
            "ci97_5": [float(lo), float(hi)],
            "kind": kind,
            "clause": HEAD_TO_HEAD[kind].format(versus=c2.FORM_WORDS[f]),
        }
    return out


def combined_reading(mult: dict[str, Any], h2h: dict[str, Any]) -> str:
    """The registered reading: the multiplicative form's sentence, then the head-to-head sentence."""
    if h2h.get("compared"):
        a, b = h2h["mean"], h2h["competition"]
        tail = HEAD_TO_HEAD["sentence"].format(
            g=f"{h2h['genes']:,}",
            n=f"{h2h['pairs']:,}",
            versus_mean=a["clause"],
            u=f"{a['mse_form_minus_mse_multiplicative']:+.4g}",
            ulo=f"{a['ci97_5'][0]:+.4g}",
            uhi=f"{a['ci97_5'][1]:+.4g}",
            versus_competition=b["clause"],
            v=f"{b['mse_form_minus_mse_multiplicative']:+.4g}",
            vlo=f"{b['ci97_5'][0]:+.4g}",
            vhi=f"{b['ci97_5'][1]:+.4g}",
        )
    else:
        tail = HEAD_TO_HEAD["not_compared"].format(g=f"{h2h.get('genes', 0):,}", min=MIN_GENES)
    return mult["reading"] + ". " + tail + "."


def comparator_problems(rows: dict[str, list[dict]], committed: dict[str, Any], recomputed: str) -> list[str]:
    """COMPARATOR: 7d's populations, P and D, and its reading, recomputed, against the committed result."""
    out = []
    for f in COMPARATOR_FORMS:
        table = committed["genes_tables"][f]
        cols = table["columns"]
        gi, pi, di = cols.index("gene"), cols.index("P"), cols.index("D")
        want = [r[gi] for r in table["rows"]]
        got = [r["gene"] for r in rows[f]]
        if got != want:
            out.append(
                f"{f}: population of {len(got)} genes against the committed {len(want)}, or another order"
            )
            continue
        worst = max(
            (
                abs(a - b)
                for r, c in zip(rows[f], table["rows"], strict=True)
                for k, j in (("P", pi), ("D", di))
                for a, b in zip(r[k], c[j], strict=True)
            ),
            default=0.0,
        )
        if not worst <= COMPARATOR_TOLERANCE:
            out.append(f"{f}: P or D differs from the committed table by {worst:.3g}")
    if recomputed != committed["reading"]:
        out.append("7d's reading, recomputed, differs from the committed reading")
    return out
