# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Multi-target compiled choice (ROADMAP §5, after item 7c M1c; lane-targetchoice, 2026-10-10).

M1c (`signaltrace.py`, 9637584 registration, 1103719 result) traced the 1,918 held-out ENCODE CRISPRi
pairs in K562 through each compiled stage with the anchor's metric. The anchor, the AlphaGenome
deletion feature, gains +0.1361 [0.0749, 0.2258] in average precision over activity and distance. The
compiled programme keeps ONE target gene per element, and that stage (s2_target) kept 39 of the 1,918
pairs and weakened the gain to +0.0792 [0.0282, 0.1562]: 1,209 pairs sit on compiled elements whose
single rule names another gene, 670 on no compiled element.

The single target is a property of the compact projection, not of the sweep: the response cache the
compiled rules were read from holds every gene in the scorer's window (`targets.ElementResponses`).
This module asks whether a compiled target choice that keeps every candidate target the compiler's own
bar names - every protein-coding gene whose predicted change on deletion reaches
`enhancer_target.MIN_EFFECT`, by the size rule `predict_target` names the compiled target with - ranks
the same held-out pairs better than the one-target choice does. Nothing is tuned: the bar is the
compiler's. A secondary choice, the first k candidate targets per element, takes its k from the
TRAINING split only (`TOP_K`, fixed by the census before the registration).

Only which pairs a choice keeps changes between the two: the value column is the anchor's own
`deletion_drop`, as in signaltrace's s2_target, so a difference is a difference of target choice.

In the order written and committed: the code (`chromosome_choice`, `trace_pair`, `census`, `score`),
then the training census (data/results/targetchoice_census.json), then the registration (`QUESTION` ..
`OUTCOMES`, data/results/targetchoice_registration.json), then one run on the held-out pairs.

No model request, no key, no money: every input is on this machine.
"""

from __future__ import annotations

import json
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any

from genomeos.attribution import compile as cp
from genomeos.attribution import crispri
from genomeos.attribution import signaltrace as st
from genomeos.predict.enhancer_target import MIN_EFFECT

REGISTERED = "2026-10-10"
LANE = "lane-targetchoice"
CELL = st.CELL

#: the choices, by the column prefix each writes
ONE = "one_target"
MULTI = "multi_target"
TOPK = "top_k"
#: the k of the secondary choice, chosen on the training split by CENSUS_RULE before the registration;
#: None until the census has run, and the secondary choice is then not read
TOP_K: int | None = None
#: the k the census tries, smallest first (ties go to the smaller k)
K_TRIED = tuple(range(1, 11))

# ---- the comparator and the anchor, as committed by M1c ----------------------------------------------------
COMPARATOR = {
    "result": "data/results/signaltrace_result.json",
    "registration": "data/results/signaltrace_registration.json",
    "stage": "s2_target",
    "kept": 39,
    "kept_regulated": 35,
    "gain": 0.0792,
    "ci95": [0.0282, 0.1562],
    "against_anchor": -0.0569,
    "against_anchor_ci95": [-0.097, -0.0248],
}
ANCHOR = st.ANCHOR
BASE_COLUMNS = st.BASE_COLUMNS
ANCHOR_COLUMNS = st.ANCHOR_COLUMNS

# ---- the registered comparison -----------------------------------------------------------------------------
QUESTION = (
    "M1c found that the compiled programme's target choice, one target gene per element, keeps 39 of the"
    " 1,918 held-out ENCODE CRISPRi pairs in K562 and weakens the anchor's gain over activity and distance"
    " from +0.1361 [95% CI 0.0749, 0.2258] to +0.0792 [0.0282, 0.1562]; 1,209 of the pairs sit on compiled"
    " elements whose single rule names another gene. Does a compiled target choice that keeps every"
    " protein-coding gene the compiler's own bar names at an element (enhancer_target.MIN_EFFECT, 0.1 log2,"
    " by the size rule predict_target names the compiled target with), nothing tuned, rank the same"
    " held-out pairs better than the one-target choice, by the anchor's metric?"
)
POPULATION = (
    "the 1,918 held-out K562 pairs M1c registered (data/results/signaltrace_registration.json,"
    " population_ids and population_sha256), as crispri.load reads them from the held-out file pinned by"
    " sha256, in file order; the label is the file's Regulated column. Every choice is read on all of"
    " them: a pair a choice does not keep stays in the population with 0 in that choice's two columns. The"
    " fit rows are M1c's and the anchor's: the K562 training pairs (EPCrisprBenchmark training_K562, the"
    " split the project's CRISPRi benchmark fits on) that crispri.annotate marks covered"
)
TRAINING_SPLIT = (
    "EPCrisprBenchmark_combined_data.training_K562.GRCh38 (crispri.TRAINING): every K562 pair of the"
    " training file, of which the covered pairs are the fit rows. The census, the choice of TOP_K and"
    " every fit read only these; no held-out pair enters a fit or a choice"
)
CANDIDATES = (
    "the candidate targets of a compiled element, in order: the target its predicted rule names (the"
    " compiled programme with its experimental rule lines removed, signaltrace.PREDICTED_ONLY), then every"
    " gene of the response cache's window at that element (targets.ElementResponses, threshold 0.0 sweep)"
    " that is a protein-coding gene of the chromosome in GENCODE v50 (genome.Annotation.protein_coding,"
    " the kind of set the sweep named its compiled target from) and whose predicted size on deletion"
    " reaches the compiler's bar, MIN_EFFECT = 0.1 log2 fold change; the size is targets.genes_at_bar's,"
    " the larger of the predicted fall and the predicted rise over every track, ties to the fall,"
    " strongest first, which is the rule predict_target names the compiled target with. A gene already"
    " listed is not listed twice. An element the cache does not hold keeps its compiled target alone"
)
CHOICES: dict[str, dict[str, str]] = {
    ONE: {
        "name": "the one-target choice (M1c's s2_target)",
        "phrase": "keeping the one compiled target per element",
        "kept": (
            "an element of the compiled programme that overlaps the CRISPRi element (half-open intervals,"
            " signaltrace.ChromLookup.overlapping) carries a predicted rule whose target gene's symbol is"
            " the measured gene: signaltrace.trace_pair's s2_target, reproduced pair for pair"
        ),
        "value": "the anchor's deletion_drop on kept pairs, 0 otherwise",
    },
    MULTI: {
        "name": "the multi-target choice at the compiler's bar (primary)",
        "phrase": "keeping every coding gene at the compiler's bar",
        "kept": (
            "an element of the compiled programme that overlaps the CRISPRi element (the same elements and"
            " the same overlap search as the one-target choice) lists the measured gene among its candidate"
            " targets (CANDIDATES)"
        ),
        "value": "the anchor's deletion_drop on kept pairs, 0 otherwise",
    },
    TOPK: {
        "name": "the top-k choice, k from the training split (secondary)",
        "phrase": "keeping the first k candidate targets per element",
        "kept": (
            "an element of the compiled programme that overlaps the CRISPRi element lists the measured gene"
            " among its first TOP_K candidate targets (CANDIDATES, compiled target first)"
        ),
        "value": "the anchor's deletion_drop on kept pairs, 0 otherwise",
    },
}
CENSUS_RULE = (
    "on the fit rows only (the covered K562 training pairs): for each k of K_TRIED, the top-k choice's"
    " model (METRIC) is fitted on the fit rows and its average precision is read on the same fit rows;"
    " TOP_K is the smallest k with the largest average precision. If TOP_K is 1 the top-k choice is the"
    " one-target choice and is not read; if no compiled element reached by a fit row lists more than TOP_K"
    " candidates it is the multi-target choice and is not read. Written before the census was computed"
)
METRIC = (
    "the anchor's metric, unchanged (signaltrace.METRIC). A choice's model is activity + distance"
    " (log_distance, log_activity, activity_over_distance) plus the choice's two columns, kept (0 or 1)"
    " and value, both 0 for a pair the choice does not keep. Each model is fitted by crispri.logistic_fit"
    " (lam 1e-3, 25 rounds) on the fit rows and scored on the 1,918 population pairs; the score is average"
    " precision (crispri.average_precision). The primary reading is the multi-target choice's average"
    " precision minus the one-target choice's, with the 95% interval of crispri.gain_interval; beside it,"
    " by the same paired resampling, each choice's gain over activity + distance and its difference from"
    " the anchor"
)
BOOTSTRAP = st.BOOTSTRAP
DECISION = (
    "for the multi-target choice (primary) and, when it is read, the top-k choice (secondary), with [dlo,"
    " dhi] the 95% interval of its average precision minus the one-target choice's and [alo, ahi] that of"
    " its average precision minus the anchor's: NONE_ADDED if it keeps no population pair the one-target"
    " choice does not keep; else BETTER if dlo > 0, worded by the anchor: BETTER_BELOW_ANCHOR if ahi < 0,"
    " BETTER_ABOVE_ANCHOR if alo > 0, otherwise BETTER_NO_DIFFERENCE_FROM_ANCHOR; WORSE if dhi < 0;"
    " otherwise NO_BETTER. The falsifier is the primary's NO_BETTER, WORSE or NONE_ADDED: then the"
    " multi-target choice is no better than the one-target choice on these pairs. The question is answered"
    " by the primary alone; the secondary is read beside it and decides nothing. No margin beyond the"
    " intervals is applied"
)
STOP_RULE = (
    "one run, from this registration's committed code, on the cached files the registration pins by"
    " sha256. No choice, column, population, fit row, metric, parameter or rule is changed after the"
    " multi-target choice is read on any held-out pair; a defect found afterwards is reported beside the"
    " registered readings and corrected only under a new registration. The run stops with HARNESS_DEFECT,"
    " before the multi-target or top-k choice is read on the population, if: a pinned input moved; the"
    " population's ids or sha256 differ from M1c's registration; the removed experimental lines of a"
    " program differ from readergate.rule_targets' count; a pair kept by the one-target choice is not kept"
    " by the top-k choice, or one kept by the top-k choice is not kept by the multi-target choice; the"
    " one-target columns differ from M1c's committed s2_target columns for any population pair; the anchor"
    " does not reproduce the committed arm-1 figures exactly (1,918 pairs, 118 positives, gain 0.1361,"
    " interval [0.0749, 0.2258]); or the one-target choice does not reproduce M1c's committed s2_target"
    " reading exactly (39 kept, 35 regulated, gain 0.0792 [0.0282, 0.1562], against the anchor -0.0569"
    " [-0.097, -0.0248])"
)
#: the exact sentence of each outcome; `reading` fills the numbers and nothing else
OUTCOMES = {
    "HARNESS_DEFECT": "stopped before the multi-target choice was read on the held-out pairs: {what}",
    "NONE_ADDED": (
        "no pair added: {phrase} keeps no held-out K562 pair beyond the one-target choice's ({k} of the {n}"
        " pairs kept, {kp} of them regulated), so on these pairs target choice cannot rank the measured"
        " pairs differently; against the one-target choice {d} [95% CI {dlo}, {dhi}] in average precision"
    ),
    "BETTER_BELOW_ANCHOR": (
        "better than the one-target choice and below the anchor: {phrase} ranks the {n} held-out K562 pairs"
        " better than keeping one target per element, {d} [95% CI {dlo}, {dhi}] in average precision ({k}"
        " pairs kept, {kp} of them regulated, against {k1} and {kp1}); its gain over activity and distance"
        " is {g} [95% CI {lo}, {hi}], the one-target choice's {g1} [95% CI {lo1}, {hi1}]; it stays below"
        " the anchor, {a} [95% CI {alo}, {ahi}]"
    ),
    "BETTER_NO_DIFFERENCE_FROM_ANCHOR": (
        "better than the one-target choice: {phrase} ranks the {n} held-out K562 pairs better than keeping"
        " one target per element, {d} [95% CI {dlo}, {dhi}] in average precision ({k} pairs kept, {kp} of"
        " them regulated, against {k1} and {kp1}); its gain over activity and distance is {g} [95% CI {lo},"
        " {hi}], the one-target choice's {g1} [95% CI {lo1}, {hi1}]; against the anchor {a} [95% CI {alo},"
        " {ahi}], no difference detected"
    ),
    "BETTER_ABOVE_ANCHOR": (
        "better than the one-target choice and above the anchor: {phrase} ranks the {n} held-out K562 pairs"
        " better than keeping one target per element, {d} [95% CI {dlo}, {dhi}] in average precision ({k}"
        " pairs kept, {kp} of them regulated, against {k1} and {kp1}), and better than the anchor, {a} [95%"
        " CI {alo}, {ahi}]; its gain over activity and distance is {g} [95% CI {lo}, {hi}], the one-target"
        " choice's {g1} [95% CI {lo1}, {hi1}]"
    ),
    "NO_BETTER": (
        "no better than the one-target choice: between {phrase} and keeping one target per element no"
        " difference was detected on the {n} held-out K562 pairs, {d} [95% CI {dlo}, {dhi}] in average"
        " precision ({k} pairs kept, {kp} of them regulated, against {k1} and {kp1}); its gain over activity"
        " and distance is {g} [95% CI {lo}, {hi}], the one-target choice's {g1} [95% CI {lo1}, {hi1}];"
        " against the anchor {a} [95% CI {alo}, {ahi}]"
    ),
    "WORSE": (
        "worse than the one-target choice: {phrase} ranks the {n} held-out K562 pairs below keeping one"
        " target per element, {d} [95% CI {dlo}, {dhi}] in average precision ({k} pairs kept, {kp} of them"
        " regulated, against {k1} and {kp1}); its gain over activity and distance is {g} [95% CI {lo},"
        " {hi}], the one-target choice's {g1} [95% CI {lo1}, {hi1}]; against the anchor {a} [95% CI {alo},"
        " {ahi}]"
    ),
}
BETTER = ("BETTER_BELOW_ANCHOR", "BETTER_NO_DIFFERENCE_FROM_ANCHOR", "BETTER_ABOVE_ANCHOR")
FALSIFIED_BY = ("NO_BETTER", "WORSE", "NONE_ADDED")
LIMITATIONS = (
    "both choices carry the anchor's deletion_drop as their value, so the comparison isolates which pairs"
    " each keeps; it says nothing about the value a compiled rule states or about the stages after it,"
    " where M1c lost every pair at the bridge",
    "the protein-coding set is GENCODE v50 as the project's annotation reads it now; the sweep named its"
    " compiled targets from a coding set of its own day, so a compiled target that is not the first"
    " coding gene at the bar is possible and is counted beside the readings; the compiled target is kept"
    " whatever the set says",
    "the candidate targets and the anchor's two columns come from one AlphaGenome sweep, and AlphaGenome"
    " was trained on ENCODE; the multi-target choice is not independent of the anchor's value column",
    "the held-out file's positives were filtered on chromatin at the tested element (crispri.POSITIVE_FILTER)"
    " and the activity column reads H3K27ac from the same file; both choices share that baseline, so it"
    " cancels from the choice-to-choice difference and not from a gain over activity and distance",
    "the metric ranks all 1,918 pairs; a choice that keeps few of them can move the ranking of those few"
    " only. The kept counts stand beside every reading for that reason",
    "no choice is a test of a mechanism: a better ranking says only that the compiled layer, at target"
    " choice, ranks the measured pairs better when it keeps more targets",
)
DESCRIPTIVES = (
    "per choice: the pairs kept of the 1,918 and of the fit rows, how many of them regulated, how many"
    " with a value above 0, and its fitted weights",
    "the population pairs the multi-target choice keeps and the one-target choice does not, and how many"
    " of them regulated; where the measured gene sits among an element's candidate targets",
    "the candidate-target census of the compiled elements on the chromosomes read: how many hold 2 or more,"
    " how many are not cached, how many compiled targets are not the first coding gene at the bar",
)
#: the registration fields the run refuses to differ from: registration field -> constant of this module
BOUND = {
    "question": "QUESTION",
    "population": "POPULATION",
    "training_split": "TRAINING_SPLIT",
    "candidates": "CANDIDATES",
    "choices": "CHOICES",
    "min_effect": "MIN_EFFECT",
    "top_k": "TOP_K",
    "k_tried": "K_TRIED",
    "census_rule": "CENSUS_RULE",
    "comparator": "COMPARATOR",
    "anchor": "ANCHOR",
    "base_columns": "BASE_COLUMNS",
    "anchor_columns": "ANCHOR_COLUMNS",
    "metric": "METRIC",
    "bootstrap": "BOOTSTRAP",
    "decision": "DECISION",
    "stop_rule": "STOP_RULE",
    "outcomes": "OUTCOMES",
    "limitations": "LIMITATIONS",
    "descriptives_not_verdicts": "DESCRIPTIVES",
}


def registered_fields() -> dict[str, Any]:
    """Every BOUND constant as JSON stores it, so a registration and this module compare equal or not."""
    return {f: json.loads(json.dumps(globals()[c])) for f, c in BOUND.items()}


def read_choices(top_k: int | None) -> tuple[str, ...]:
    """The choices a run reads: the one-target comparator, the primary, and the secondary when it is read."""
    return (ONE, MULTI, TOPK) if top_k is not None and top_k > 1 else (ONE, MULTI)


# ==== the compiled programme's candidate targets, one chromosome at a time =================================


@dataclass
class ChoiceLookup:
    """The compiled elements of one chromosome, their one compiled target and their candidate targets."""

    chrom: str
    look: st.ChromLookup  # links and the overlap index, as signaltrace builds them
    candidates: dict[str, tuple[str, ...]] = field(default_factory=dict)  # element -> CANDIDATES
    counts: Counter = field(default_factory=Counter)

    def overlapping(self, start: int, end: int) -> list[str]:
        return self.look.overlapping(start, end)


def candidate_list(
    heads: Sequence[str], at_bar: Sequence[tuple[str, float]] | None, coding: set[str]
) -> tuple[str, ...]:
    """CANDIDATES for one element: its compiled targets, then the coding genes at the bar, strongest first."""
    window = [g for g, _ in (at_bar or ()) if g in coding]
    return tuple(dict.fromkeys([*heads, *window]))


def chromosome_choice(chrom: str, text: str, coding: set[str], responses: Any) -> ChoiceLookup:
    """Every compiled element of one chromosome with its candidate targets (no label is read).

    `responses` is a `targets.ElementResponses`, or anything with its `at_bar(chrom, id, min_effect)`."""
    from genomeos.lang.parser import parse

    elements = {e: (s, t) for e, _c, s, t in cp.program_elements(text)}
    pred, removed = st.predicted_text(text)
    module = parse(pred, f"noncoding_{chrom}")
    look = st.ChromLookup(chrom)
    look.links = st.links_of(module)
    del module
    look.index({e: elements[e] for e in look.links})
    out = ChoiceLookup(chrom, look)
    out.counts["experimental_rule_lines_removed"] = removed
    for el in sorted(look.links):
        heads = [lk.symbol for lk in look.links[el]]
        bar = responses.at_bar(chrom, el, MIN_EFFECT)
        out.candidates[el] = candidate_list(heads, bar, coding)
        out.counts["elements"] += 1
        out.counts["compiled_targets_1" if len(set(heads)) == 1 else "compiled_targets_2_or_more"] += 1
        if bar is None:
            out.counts["not_cached"] += 1
            continue
        window = [g for g, _ in bar if g in coding]
        if not window:
            out.counts["no_coding_gene_at_the_bar"] += 1
        elif window[0] != heads[0]:
            out.counts["compiled_target_not_first_at_the_bar"] += 1
        if heads[0] not in window:
            out.counts["compiled_target_not_at_the_bar"] += 1
    return out


def size_band(s: int) -> str:
    if s <= 3:
        return str(s)
    return "4-5" if s <= 5 else "6-10" if s <= 10 else "11+"


def size_histogram(sizes: Sequence[int]) -> dict[str, int]:
    bands = Counter(size_band(s) for s in sizes)
    return {b: bands.get(b, 0) for b in ("1", "2", "3", "4-5", "6-10", "11+")}


# ==== one pair through the choices ========================================================================


def rank_of(p: crispri.Pair, cl: ChoiceLookup | None) -> int | None:
    """The measured gene's best position (1 = the compiled target) among an overlapping element's
    candidate targets, or None when no overlapping compiled element lists it."""
    if cl is None:
        return None
    ranks = [
        cl.candidates[el].index(p.gene) + 1
        for el in cl.overlapping(p.start, p.end)
        if p.gene in cl.candidates[el]
    ]
    return min(ranks, default=None)


def trace_pair(p: crispri.Pair, cl: ChoiceLookup | None, top_k: int | None = None) -> dict[str, float]:
    """The two columns of every choice for one pair. Reads the anchor's deletion_drop (crispri.annotate
    must have run) and the chromosome's lookup; never the pair's label."""
    els = cl.overlapping(p.start, p.end) if cl else []
    one = cl is not None and any(lk.symbol == p.gene for el in els for lk in cl.look.links.get(el, []))
    rank = rank_of(p, cl)
    kept = {ONE: one, MULTI: rank is not None}
    if top_k is not None:
        kept[TOPK] = rank is not None and rank <= top_k
    f: dict[str, float] = {}
    for choice, k in kept.items():
        f[f"{choice}_kept"] = float(k)
        f[f"{choice}_drop"] = p.features["deletion_drop"] if k else 0.0
    return f


def columns(choice: str) -> tuple[str, ...]:
    return (*BASE_COLUMNS, f"{choice}_kept", f"{choice}_drop")


# ==== the harness, before the multi-target choice is read on the population ================================


def harness(pairs: Sequence[crispri.Pair], choices: Sequence[str]) -> list[str]:
    """The predictor-side subset checks of STOP_RULE, reading no label; an empty list passes."""
    problems = []
    chain = [c for c in (ONE, TOPK, MULTI) if c in choices]
    for a, b in zip(chain, chain[1:], strict=False):
        stray = sum(1 for p in pairs if p.features[f"{a}_kept"] and not p.features[f"{b}_kept"])
        if stray:
            problems.append(f"{stray} pairs kept by {a} are not kept by {b}")
    return problems


def comparator_mismatches(pop: Sequence[crispri.Pair], table: dict[str, Any]) -> list[str]:
    """Population pairs whose one-target columns differ from M1c's committed s2_target columns."""
    cols = table["columns"]
    i_id, i_k, i_d = cols.index("pair"), cols.index("s2_target_kept"), cols.index("s2_target_drop")
    committed = {row[i_id]: (float(row[i_k]), float(row[i_d])) for row in table["rows"]}
    out = []
    for p in pop:
        mine = (p.features[f"{ONE}_kept"], round(p.features[f"{ONE}_drop"], 6))
        if committed.get(st.pair_id(p)) != mine:
            out.append(st.pair_id(p))
    return out


# ==== the training census (fit rows only) =================================================================


def k_reading(fit_rows: Sequence[crispri.Pair], choice: str, ap_base: float) -> dict[str, Any]:
    """One choice's in-sample reading on the fit rows (CENSUS_RULE)."""
    labels = [p.regulated for p in fit_rows]
    w, s = st.fit_scores(fit_rows, fit_rows, columns(choice))
    ap = crispri.average_precision(s, labels) or 0.0
    kept = [p for p in fit_rows if p.features[f"{choice}_kept"]]
    return {
        "kept": len(kept),
        "kept_regulated": sum(p.regulated for p in kept),
        "average_precision_in_sample": round(ap, 4),
        "gain_over_activity_distance_in_sample": round(ap - ap_base, 4),
        "weights": [round(x, 6) for x in w],
    }


def choose_k(by_k: dict[str, dict[str, Any]]) -> int:
    """CENSUS_RULE: the smallest k of K_TRIED with the largest in-sample average precision."""
    best = max(by_k[str(k)]["average_precision_in_sample"] for k in K_TRIED)
    return min(k for k in K_TRIED if by_k[str(k)]["average_precision_in_sample"] == best)


def secondary_status(top_k: int, largest_reached: int) -> str:
    if top_k == 1:
        return "not read: TOP_K is 1, so the top-k choice is the one-target choice"
    if largest_reached <= top_k:
        return (
            "not read: no compiled element reached by a fit row lists more than TOP_K candidates, so the"
            " top-k choice is the multi-target choice"
        )
    return "read beside the primary"


def element_census(sizes: Sequence[int]) -> dict[str, Any]:
    return {
        "elements": len(sizes),
        "candidate_targets": size_histogram(sizes),
        "with_2_or_more": sum(1 for s in sizes if s >= 2),
        "largest": max(sizes, default=0),
    }


def census(lookups: dict[str, ChoiceLookup], fit_rows: Sequence[crispri.Pair]) -> dict[str, Any]:
    """CENSUS_RULE on the fit rows, with the candidate-target census of the compiled elements. Reads the
    training labels and nothing of the held-out file."""
    fit_rows = list(fit_rows)
    all_sizes = [len(c) for cl in lookups.values() for c in cl.candidates.values()]
    reached: set[tuple[str, str]] = set()
    for p in fit_rows:
        cl = lookups.get(p.chrom)
        if cl:
            reached.update((p.chrom, el) for el in cl.overlapping(p.start, p.end))
    reached_sizes = [len(lookups[c].candidates[el]) for c, el in sorted(reached)]
    ranks: Counter = Counter()
    for p in fit_rows:
        r = rank_of(p, lookups.get(p.chrom))
        ranks["not listed" if r is None else size_band(r)] += 1
    counts: Counter = Counter()
    for cl in lookups.values():
        counts.update(cl.counts)
    _, base = st.fit_scores(fit_rows, fit_rows, BASE_COLUMNS)
    ap_base = crispri.average_precision(base, [p.regulated for p in fit_rows]) or 0.0
    by_k: dict[str, Any] = {}
    for k in K_TRIED:
        for p in fit_rows:
            p.features.update(trace_pair(p, lookups.get(p.chrom), top_k=k))
        by_k[str(k)] = k_reading(fit_rows, TOPK, ap_base)
    by_k["all"] = k_reading(fit_rows, MULTI, ap_base)
    one = k_reading(fit_rows, ONE, ap_base)
    top_k = choose_k(by_k)
    reached_census = element_census(reached_sizes)
    keys = (
        "not_cached",
        "no_coding_gene_at_the_bar",
        "compiled_target_not_first_at_the_bar",
        "compiled_target_not_at_the_bar",
        "compiled_targets_1",
        "compiled_targets_2_or_more",
        "experimental_rule_lines_removed",
    )
    return {
        "census_rule": CENSUS_RULE,
        "fit_rows": {"pairs": len(fit_rows), "positives": sum(p.regulated for p in fit_rows)},
        "average_precision_activity_distance_in_sample": round(ap_base, 4),
        "compiled_elements_on_the_chromosomes_read": {
            **element_census(all_sizes),
            **{k: counts.get(k, 0) for k in keys},
        },
        "compiled_elements_reached_by_a_fit_row": reached_census,
        "fit_rows_by_the_measured_genes_position_among_candidates": {
            b: ranks.get(b, 0) for b in ("1", "2", "3", "4-5", "6-10", "11+", "not listed")
        },
        "one_target_in_sample": one,
        "by_k_in_sample": by_k,
        "top_k": top_k,
        "secondary": secondary_status(top_k, reached_census["largest"]),
    }


# ==== the readings ========================================================================================


def kept_counts(rows: Sequence[crispri.Pair], choice: str) -> dict[str, int]:
    kept = [p for p in rows if p.features[f"{choice}_kept"]]
    return {
        "pairs": len(rows),
        "kept": len(kept),
        "kept_regulated": sum(p.regulated for p in kept),
        "kept_with_value_above_0": sum(1 for p in kept if p.features[f"{choice}_drop"] > 0),
        "kept_regulated_with_value_above_0": sum(
            1 for p in kept if p.regulated and p.features[f"{choice}_drop"] > 0
        ),
    }


def decide(added: int, vs_one: Sequence[float] | None, vs_anchor: Sequence[float] | None) -> str:
    if added == 0:
        return "NONE_ADDED"
    if vs_one is None or vs_anchor is None:
        raise ValueError(
            "an interval is missing: the population spans fewer chromosomes than an interval needs"
        )
    dlo, dhi = vs_one
    alo, ahi = vs_anchor
    if dlo > 0:
        if ahi < 0:
            return "BETTER_BELOW_ANCHOR"
        return "BETTER_ABOVE_ANCHOR" if alo > 0 else "BETTER_NO_DIFFERENCE_FROM_ANCHOR"
    if dhi < 0:
        return "WORSE"
    return "NO_BETTER"


def reading(kind: str, **numbers: Any) -> str:
    return OUTCOMES[kind].format(**numbers)


def _signed(x: float) -> str:
    return f"{x:+.4f}"


def _triple(names: tuple[str, str, str], interval: dict[str, Any]) -> dict[str, str]:
    """A gain_interval's point and bounds under the three placeholder names an OUTCOMES sentence uses."""
    point, lo, hi = names
    return {
        point: _signed(interval["gain"]),
        lo: _signed(interval["ci95"][0]),
        hi: _signed(interval["ci95"][1]),
    }


def numbers(row: dict[str, Any], one: dict[str, Any], phrase: str) -> dict[str, Any]:
    return {
        "phrase": phrase,
        **_triple(("g", "lo", "hi"), row["gain_over_activity_distance"]),
        **_triple(("d", "dlo", "dhi"), row["against_one_target"]),
        **_triple(("a", "alo", "ahi"), row["against_anchor"]),
        **_triple(("g1", "lo1", "hi1"), one["gain_over_activity_distance"]),
        "n": f"{row['population']['pairs']:,}",
        "k": f"{row['population']['kept']:,}",
        "kp": f"{row['population']['kept_regulated']:,}",
        "k1": f"{one['population']['kept']:,}",
        "kp1": f"{one['population']['kept_regulated']:,}",
    }


def phrase_of(choice: str, top_k: int | None) -> str:
    if choice == TOPK:
        return f"keeping the first {top_k} candidate targets per element"
    return CHOICES[choice]["phrase"]


def comparator_reproduced(one: dict[str, Any]) -> bool:
    c = COMPARATOR
    return (
        one["population"]["kept"] == c["kept"]
        and one["population"]["kept_regulated"] == c["kept_regulated"]
        and one["gain_over_activity_distance"]["gain"] == c["gain"]
        and one["gain_over_activity_distance"]["ci95"] == c["ci95"]
        and one["against_anchor"]["gain"] == c["against_anchor"]
        and one["against_anchor"]["ci95"] == c["against_anchor_ci95"]
    )


def score_one(fit_rows: Sequence[crispri.Pair], pop: Sequence[crispri.Pair], anchor: dict[str, Any]) -> dict:
    """The one-target choice's reading, computed first: it must reproduce M1c's s2_target exactly."""
    pop = list(pop)
    w, s = st.fit_scores(fit_rows, pop, columns(ONE))
    return {
        "name": CHOICES[ONE]["name"],
        "columns": list(columns(ONE)),
        "weights": [round(x, 6) for x in w],
        "population": kept_counts(pop, ONE),
        "fit_rows": kept_counts(fit_rows, ONE),
        "gain_over_activity_distance": crispri.gain_interval(s, anchor["base"], pop),
        "against_anchor": crispri.gain_interval(s, anchor["scores"], pop),
        "scores": s,
    }


def score(
    fit_rows: Sequence[crispri.Pair],
    pop: Sequence[crispri.Pair],
    anchor: dict[str, Any],
    one: dict[str, Any],
    choices: Sequence[str],
    top_k: int | None,
) -> dict[str, Any]:
    """Every choice beyond the one-target one (METRIC), its outcome (DECISION), and the primary verdict."""
    pop = list(pop)
    rows: dict[str, Any] = {}
    for choice in choices:
        if choice == ONE:
            continue
        w, s = st.fit_scores(fit_rows, pop, columns(choice))
        added = [p for p in pop if p.features[f"{choice}_kept"] and not p.features[f"{ONE}_kept"]]
        row: dict[str, Any] = {
            "name": CHOICES[choice]["name"],
            "columns": list(columns(choice)),
            "weights": [round(x, 6) for x in w],
            "population": kept_counts(pop, choice),
            "fit_rows": kept_counts(fit_rows, choice),
            "population_pairs_added_over_one_target": len(added),
            "population_pairs_added_regulated": sum(p.regulated for p in added),
            "gain_over_activity_distance": crispri.gain_interval(s, anchor["base"], pop),
            "against_one_target": crispri.gain_interval(s, one["scores"], pop),
            "against_anchor": crispri.gain_interval(s, anchor["scores"], pop),
        }
        row["outcome"] = decide(len(added), row["against_one_target"]["ci95"], row["against_anchor"]["ci95"])
        row["reading"] = reading(row["outcome"], **numbers(row, one, phrase_of(choice, top_k)))
        rows[choice] = row
    primary = rows[MULTI]
    return {
        "choices": rows,
        "outcome": primary["outcome"],
        "reading": primary["reading"],
        "falsified": primary["outcome"] in FALSIFIED_BY,
    }
