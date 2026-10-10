# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""M1c (ROADMAP §5 item 7c, lane-signaltrace, 2026-10-10): where the compiled layer loses its signal.

Four registered tests found nothing the compiled attribution layer predicts beyond its baselines: M1
(runrules, 747b207) no difference detected on GTEx RNA; M2b (readergate, ff35b09) worse on ENCODE
polyA RNA; M3 (closure_cell, 80bf3db) no difference detected on a CRISPRi growth screen; M1b
(regulators, a951f52) not factor-specific. The layer's one positive independent result sits where its
input is closest to a measurement: the AlphaGenome deletion feature on the held-out ENCODE CRISPRi
pairs in K562 (+0.1361 [0.0749, 0.2258] on 1,918 pairs, data/results/crispri_published_v2.json).

This module follows that ONE measured population through each stage of the compiled layer, with the
anchor's own metric at every stage, and names the first stage at which the agreement with the
measured outcome is no longer detectable. It re-implements no stage: the target is what the compiled
programs (data/knowledge/compiled) name, the strength is what `bridge.parameterize` fits, the gate is
what `compile.gate_program`'s reader block states, the regulator is what `regulators.candidates`
assigns, the run is `readergate.run_gene`. Its own code is the join of a pair to the compiled elements
it overlaps, the two columns each stage contributes, and the readings.

In the order written and committed: the code (`chromosome_lookup`, `trace_pair`, `score`), then the
registration (`QUESTION` .. `OUTCOMES`, data/results/signaltrace_registration.json), then one run.

No model request, no key, no money: every input is on this machine.
"""

from __future__ import annotations

import bisect
import hashlib
import json
import math
from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from typing import Any

from genomeos.attribution import bridge, crispri
from genomeos.attribution import compile as cp
from genomeos.attribution import readergate as rgate
from genomeos.attribution import regulators as regs
from genomeos.attribution import runrules as rr
from genomeos.ir import Action

REGISTERED = "2026-10-10"
LANE = "lane-signaltrace"
CELL = "K562"
CONTEXT = {"cell_type": CELL}
REGULATORY = (Action.ACTIVATE, Action.INHIBIT)

# ---- the anchor, as committed ----------------------------------------------------------------------------
ANCHOR = {
    "result": "data/results/crispri_published_v2.json",
    "field": "coverage_arms_k562_heldout/arm1_all_pairs",
    "pairs": 1918,
    "positives": 118,
    "uncovered_pairs": 174,
    "gain": 0.1361,
    "ci95": [0.0749, 0.2258],
    "model": "activity + distance + deletion against activity + distance, weights fitted on the covered"
    " K562 training pairs (crispri.score_published), average precision, crispri.gain_interval",
}
BASE_COLUMNS = crispri.FEATURES["activity + distance"]
ANCHOR_COLUMNS = crispri.FEATURES["activity + distance + deletion"]

# ---- the registered trace --------------------------------------------------------------------------------
QUESTION = (
    "The compiled layer's one positive independent result is the AlphaGenome deletion feature on the"
    " held-out ENCODE CRISPRi element-gene pairs in K562 (+0.1361 [95% CI 0.0749, 0.2258] on 1,918 pairs,"
    " crispri_published_v2.json, arm 1). Followed through each stage of the compiled layer - the target"
    " the compiled programme names, the value it states, the strength the bridge fits in K562, the reader"
    " gate in K562, the regulator assignment in K562 and the run's predicted change in K562 - at which"
    " stage does the agreement with the measured CRISPRi outcome stop being detectable over activity and"
    " distance?"
)
POPULATION = (
    "the pairs of the ENCODE held-out file EPCrisprBenchmark_combined_data.heldout_5_cell_types.GRCh38"
    " (pinned by sha256) whose CellType is K562, as crispri.load reads them (ValidConnection TRUE), in file"
    " order: 1,918 pairs, each identified as chrom:chromStart-chromEnd:measuredGeneSymbol (unique). The list"
    " and the sha256 of its newline-joined form are in the registration. The label is the file's Regulated"
    " column (crispri.Pair.regulated), the anchor's own label. Every stage is read on all of them: a pair a"
    " stage drops stays in the population with 0 in that stage's two columns. The fit rows are the anchor's:"
    " the K562 training pairs crispri.annotate marks covered, unchanged for every stage"
)
PREDICTED_ONLY = (
    "every stage after the anchor reads the compiled programs with their experimental rule lines removed."
    " Those rules carry the CRISPRi screens' own EffectSize, the held-out file's pairs included (marked"
    " 'held-out split, evaluation only, never a feature'), and the bridge prefers an experimental rule to a"
    " predicted one for the same mechanism (bridge.PRECEDENCE), so a stage that read them would read the"
    " outcome. The removed lines are counted per program and must equal readergate.rule_targets' count of"
    " experimental rule lines; nothing else in a program is changed"
)
#: the stages, in the order the compiled layer applies them. `follows` is the stage whose kept pairs this
#: one starts from; the chain is CHAIN, and a branch is read beside it and feeds nothing downstream.
STAGES: dict[str, dict[str, str]] = {
    "s1_anchor": {
        "name": "the deletion feature (anchor)",
        "follows": "",
        "kept": (
            "crispri top_target: the measured gene is the single top predicted target of a deleted element"
            " overlapping the CRISPRi element, in the compact all-element table"
        ),
        "value": (
            "crispri deletion_drop: the largest predicted fall of the measured gene on deleting an"
            " overlapping element, on K562's own track, from the per-element cache (any gene the sweep"
            " scored)"
        ),
    },
    "s2_target": {
        "name": "target choice",
        "follows": "s1_anchor",
        "kept": (
            "an element of the compiled programme that overlaps the CRISPRi element (half-open intervals)"
            " carries a predicted rule whose target gene's symbol is the measured gene"
        ),
        "value": (
            "the anchor's deletion_drop, on kept pairs only, so this stage changes only which pairs are kept"
        ),
    },
    "s3_compiled_value": {
        "name": "the value the compiled rule states",
        "follows": "s2_target",
        "kept": "as s2_target",
        "value": (
            "the compiled rule's own strength (|log2 fold change| in its compiled cell, the argmax track over"
            " the scorer's RNA-seq tracks, clipped at 1) for an activating rule and 0 for an inhibiting one,"
            " the largest over the kept rules"
        ),
    },
    "s4_bridge": {
        "name": "the strength the bridge fits in K562",
        "follows": "s3_compiled_value",
        "kept": (
            "a kept rule whose mechanism bridge.parameterize resolves in context cell_type=K562 on the gated"
            " program (compile.gate_program, selection reader_cells: every predicted rule is offered in K562"
            " whatever its compiled cell), on the borrowed tier (data/results/gene_rates.json and its genome"
            " median), one program at a time as runrules.run_all bridges"
        ),
        "value": (
            "max(0, -log2 RHO) of the bridged observation; the fitted strength itself is DECLARED_STRENGTH,"
            " 1.0, for every resolved mechanism, so RHO is what the bridge carries"
        ),
    },
    "s4a_bridge_as_compiled": {
        "name": "the strength the bridge fits in K562, rules as compiled",
        "follows": "s3_compiled_value",
        "kept": (
            "a kept rule whose mechanism bridge.parameterize resolves in context cell_type=K562 on the"
            " predicted program as compiled (each rule active only in its compiled cell, the path M1 and M1b"
            " ran), on the borrowed tier"
        ),
        "value": "max(0, -log2 RHO) of the bridged observation",
    },
    "s5_reader": {
        "name": "the reader gate in K562",
        "follows": "s4_bridge",
        "kept": (
            "a rule kept by s4_bridge whose element the gated program's K562 reader block states open (an"
            " ENCODE DNase-seq narrowPeak overlaps it); closed and unknown hold the rule"
        ),
        "value": "as s4_bridge, over the open rules",
    },
    "s6_regulator": {
        "name": "the regulator assignment in K562",
        "follows": "s5_reader",
        "kept": (
            "a rule kept by s5_reader that gains a source in K562: regulators.candidates, with the rule's"
            " cell set to K562 (the cell it is offered in), gives at least one factor a source kind"
            " (regulators.SOURCE_KINDS), i.e. an ENCODE TF ChIP-seq peak cluster naming K562 overlaps the"
            " element"
        ),
        "value": "as s5_reader, over those rules",
    },
    "s7_run": {
        "name": "the run's predicted change in K562",
        "follows": "s5_reader",
        "kept": (
            "a rule kept by s4_bridge whose gated run in K562 (readergate.run_gene on"
            " readergate.gated_gene_module) has its element open"
        ),
        "value": (
            "max(0, log2(L_K562 / L_held)): the runtime's level of the gene in K562 with the element present"
            " and the reader applied, over its level with the rule held, which is its level with the element"
            " removed; the largest over the run rules"
        ),
    },
}
CHAIN = ("s1_anchor", "s2_target", "s3_compiled_value", "s4_bridge", "s5_reader", "s7_run")
BRANCHES = ("s4a_bridge_as_compiled", "s6_regulator")
WHY_THE_RUN_FOLLOWS_THE_READER = (
    "the run reads no regulator: compile.py emits no regulator rule (regulators.COMPILE_SPEC), so the run"
    " follows s5_reader and s6_regulator is read beside the chain on s5_reader's pairs"
)
METRIC = (
    "the anchor's metric, unchanged. A stage's model is activity + distance (log_distance, log_activity,"
    " activity_over_distance) plus the stage's two columns, kept (0 or 1) and value, both 0 for a pair the"
    " stage drops (the benchmark's minimum-score rule, as the anchor treats an uncovered pair); for the"
    " anchor the two columns are top_target and deletion_drop. Each model is fitted by crispri.logistic_fit"
    " (lam 1e-3, 25 rounds) on the fit rows and scored on the 1,918 population pairs; the score is average"
    " precision (crispri.average_precision). A stage's reading is its average precision minus the activity"
    " + distance model's (its gain), with the 95% interval of crispri.gain_interval (2,000 resamples of"
    " whole chromosomes, seed 0); beside it, by the same paired resampling, the stage minus the anchor and"
    " the stage minus the stage it follows"
)
BOOTSTRAP = {
    "function": "crispri.gain_interval",
    "resamples": crispri.BOOTSTRAPS,
    "seed": 0,
    "unit": "chromosome: whole chromosomes resampled with replacement, the anchor's unit",
    "interval": "95% percentile interval",
}
LOSS = ("NONE_KEPT", "LOST", "WORSE")
DECISION = (
    "per stage after the anchor, with [lo, hi] the 95% interval of its gain over activity + distance and"
    " [dlo, dhi] that of its average precision minus the anchor's: NONE_KEPT if the stage keeps none of the"
    " 1,918 pairs; else WORSE if hi < 0; WEAKENED if lo > 0 and dhi < 0; CARRIES if lo > 0; LOST if dhi < 0;"
    " otherwise UNDECIDED. The chain (s2_target, s3_compiled_value, s4_bridge, s5_reader, s7_run) is walked"
    " in order and stops at the first stage whose outcome is neither CARRIES nor WEAKENED: LOST_AT that"
    " stage if its outcome is NONE_KEPT, LOST or WORSE, UNDECIDED_AT it if UNDECIDED. If no stage stops it:"
    " WEAKENED_NOT_LOST if any chain stage is WEAKENED, else CARRIED_TO_THE_RUN. 'Loses the signal' means"
    " exactly this: the stage's interval of gain over activity and distance no longer lies above zero, and"
    " its interval against the anchor lies below zero. No margin beyond the intervals is applied"
)
CLOSED_FORM_TOLERANCE = rr.CLOSED_FORM_TOLERANCE
STOP_RULE = (
    "one run, from this registration's committed code, on the cached files the registration pins by"
    " sha256. No stage, column, population, fit row, metric or rule is changed after any stage beyond the"
    " anchor is read; a defect found afterwards is reported beside the registered readings and corrected"
    " only under a new registration. The run stops with HARNESS_DEFECT, before any stage beyond the anchor"
    " is read, if: the population's sha256 differs from the registration's; the removed experimental lines"
    " of a program differ from readergate.rule_targets' count; a run differs from its closed form"
    " (-log2 RHO where its element is open in K562, 0 where held) by more than CLOSED_FORM_TOLERANCE log2,"
    " states another K562 state than the reader block, or returns no run for a resolved rule; a stage keeps"
    " a pair the stage it follows dropped; or the anchor does not reproduce the committed arm-1 figures"
    " exactly (1,918 pairs, 118 positives, gain 0.1361, interval [0.0749, 0.2258])"
)
#: the exact sentence of each outcome; `reading` fills the numbers and nothing else
OUTCOMES = {
    "HARNESS_DEFECT": "stopped before any stage beyond the anchor was read: {what}",
    "ANCHOR": (
        "anchor reproduced: the deletion feature's gain over activity and distance is {g} [95% CI {lo}, {hi}]"
        " on the {n} held-out K562 pairs ({k} with the measured gene as an overlapping element's top"
        " predicted target, {kp} of them regulated), as committed"
    ),
    "NONE_KEPT": (
        "none kept: the stage keeps none of the {n} held-out K562 pairs ({kf} of the {nf} fit rows), so no"
        " compiled prediction reaches the measurement at this stage; gain over activity and distance {g}"
        " [95% CI {lo}, {hi}], against the anchor {d} [95% CI {dlo}, {dhi}]"
    ),
    "CARRIES": (
        "carries the anchor's signal: gain over activity and distance {g} [95% CI {lo}, {hi}] on the {n}"
        " held-out K562 pairs ({k} kept by the stage, {kp} of them regulated); against the anchor {d}"
        " [95% CI {dlo}, {dhi}], no difference detected"
    ),
    "WEAKENED": (
        "weakened: a gain over activity and distance remains, {g} [95% CI {lo}, {hi}] on the {n} held-out"
        " K562 pairs ({k} kept by the stage, {kp} of them regulated), and it is below the anchor's, {d}"
        " [95% CI {dlo}, {dhi}]"
    ),
    "LOST": (
        "lost: no gain over activity and distance detected, {g} [95% CI {lo}, {hi}] on the {n} held-out K562"
        " pairs ({k} kept by the stage, {kp} of them regulated), and below the anchor, {d} [95% CI {dlo},"
        " {dhi}]"
    ),
    "UNDECIDED": (
        "undecided: no gain over activity and distance detected, {g} [95% CI {lo}, {hi}] on the {n} held-out"
        " K562 pairs ({k} kept by the stage, {kp} of them regulated), and no difference from the anchor"
        " detected either, {d} [95% CI {dlo}, {dhi}]"
    ),
    "WORSE": (
        "worse: below activity and distance alone, {g} [95% CI {lo}, {hi}] on the {n} held-out K562 pairs"
        " ({k} kept by the stage, {kp} of them regulated); against the anchor {d} [95% CI {dlo}, {dhi}]"
    ),
}
CHAIN_OUTCOMES = {
    "LOST_AT": "the measured signal is lost at {stage} ({name}): {stage_reading}{before}",
    "UNDECIDED_AT": (
        "the chain cannot say where the signal goes from {stage} ({name}) on: {stage_reading}{before}"
    ),
    "WEAKENED_NOT_LOST": (
        "the measured signal is weakened, first at {stage} ({name}), and not lost before the run: at the"
        " run, {run_reading}"
    ),
    "CARRIED_TO_THE_RUN": (
        "the measured signal is carried through every stage to the run: at the run, {run_reading}"
    ),
}
BEFORE_CLAUSE = "; it was first weakened at {stage} ({name})"
LIMITATIONS = (
    "the metric ranks all 1,918 pairs; a stage that keeps few of them can move the ranking of those few"
    " only, so a loss at such a stage says that the stage's predictions reach too few of the measured"
    " pairs, as much as anything about the values it carries. The kept counts stand beside every reading"
    " for that reason",
    "a K562 reader state is ENCODE DNase-seq peak overlap and a regulator is ENCODE TF ChIP-seq occupancy:"
    " accessibility and occupancy, not function. AlphaGenome was trained on ENCODE, so neither is"
    " independent of the deletion feature",
    "the held-out file's positives were filtered on chromatin at the tested element (crispri.POSITIVE_FILTER,"
    " found post hoc on 2026-09-27) and the activity column reads H3K27ac from the same file; every stage"
    " shares that baseline, so it cancels from the stage-to-stage comparisons and not from a gain over"
    " activity and distance",
    "the bridge's rates are mouse NIH 3T3 or their genome median; the predicted change does not depend on"
    " them (runrules.READING), and no level in molecules per cell is quoted",
    "no stage is a test of a mechanism: a stage that carries the signal says only that the compiled layer,"
    " at that stage, still ranks the measured pairs better than activity and distance do",
)
DESCRIPTIVES = (
    "per stage: the pairs kept of the 1,918 and of the fit rows, how many of them regulated, how many with"
    " a value above 0, and the reasons the pairs it drops were dropped",
    "per stage: the fitted weights of its model",
    "the bridge's unresolved reasons for the rules s4_bridge and s4a_bridge_as_compiled drop",
    "how many population pairs the compiled programs' experimental layer already names in K562 (the lines"
    " PREDICTED_ONLY removes)",
)
#: the registration fields the run refuses to differ from: registration field -> constant of this module
BOUND = {
    "question": "QUESTION",
    "population": "POPULATION",
    "predicted_only": "PREDICTED_ONLY",
    "cell": "CELL",
    "stages": "STAGES",
    "chain": "CHAIN",
    "branches": "BRANCHES",
    "why_the_run_follows_the_reader": "WHY_THE_RUN_FOLLOWS_THE_READER",
    "anchor": "ANCHOR",
    "base_columns": "BASE_COLUMNS",
    "anchor_columns": "ANCHOR_COLUMNS",
    "metric": "METRIC",
    "bootstrap": "BOOTSTRAP",
    "decision": "DECISION",
    "loss": "LOSS",
    "closed_form_tolerance": "CLOSED_FORM_TOLERANCE",
    "stop_rule": "STOP_RULE",
    "outcomes": "OUTCOMES",
    "chain_outcomes": "CHAIN_OUTCOMES",
    "before_clause": "BEFORE_CLAUSE",
    "limitations": "LIMITATIONS",
    "descriptives_not_verdicts": "DESCRIPTIVES",
}


def registered_fields() -> dict[str, Any]:
    """Every BOUND constant as JSON stores it, so a registration and this module compare equal or not."""
    return {f: json.loads(json.dumps(globals()[c])) for f, c in BOUND.items()}


# ==== the compiled programme, one chromosome at a time ====================================================


def pair_id(p: crispri.Pair) -> str:
    return f"{p.chrom}:{p.start}-{p.end}:{p.gene}"


def population_ids(pairs: Sequence[crispri.Pair]) -> list[str]:
    ids = [pair_id(p) for p in pairs]
    if len(set(ids)) != len(ids):
        raise ValueError("the population's pair ids are not unique")
    return ids


def ids_sha256(ids: Sequence[str]) -> str:
    return hashlib.sha256("\n".join(ids).encode()).hexdigest()


def is_experimental_rule(line: str) -> bool:
    return line.startswith("rule ") and 'evidence: experimental "' in line


def predicted_text(text: str) -> tuple[str, int]:
    """The program with every experimental rule line removed (PREDICTED_ONLY), and how many were removed.

    Refuses when the count differs from `readergate.rule_targets`, the committed reader of rule lines."""
    kept, removed = [], 0
    for ln in text.splitlines():
        if is_experimental_rule(ln):
            removed += 1
            continue
        kept.append(ln)
    expected = sum(1 for _, _, kind in rgate.rule_targets(text) if kind == "experimental")
    if removed != expected:
        raise ValueError(f"removed {removed} experimental rule lines where rule_targets counts {expected}")
    return "\n".join(kept) + "\n", removed


@dataclass(frozen=True, slots=True)
class Link:
    """One predicted compiled rule: element -> gene, as the parsed program states it."""

    element: str
    gene_id: str
    symbol: str
    action: str  # "activates" | "inhibits"
    strength: float
    cell: str  # its compiled `when: cell_type`

    @property
    def key(self) -> tuple[str, str]:
        return (self.element, self.gene_id)


@dataclass
class ChromLookup:
    """What the compiled layer says on one chromosome, for the rules the traced pairs can reach."""

    chrom: str
    starts: list[int] = field(default_factory=list)
    spans: list[tuple[int, int, str]] = field(default_factory=list)  # (start, end, element), by start
    links: dict[str, list[Link]] = field(default_factory=dict)  # element -> its predicted rules
    bridge_reader: dict[tuple[str, str], float] = field(default_factory=dict)  # resolved -> RHO
    bridge_compiled: dict[tuple[str, str], float] = field(default_factory=dict)
    why_reader: dict[tuple[str, str], str] = field(default_factory=dict)  # unresolved -> reason
    why_compiled: dict[tuple[str, str], str] = field(default_factory=dict)
    reader: dict[str, str] = field(default_factory=dict)  # element -> K562 reader state
    factors: dict[str, list[str]] = field(default_factory=dict)  # element -> source factors in K562
    runs: dict[tuple[str, str], dict[str, Any]] = field(default_factory=dict)
    measured: list[tuple[int, int, str, str]] = field(default_factory=list)  # (start, end, symbol, cell)
    counts: Counter = field(default_factory=Counter)

    def index(self, elements: dict[str, tuple[int, int]]) -> None:
        self.spans = sorted((s, e, el) for el, (s, e) in elements.items())
        self.starts = [s for s, _, _ in self.spans]
        longest = max((e - s for s, e, _ in self.spans), default=0)
        if longest > crispri.REACH:
            raise ValueError(f"{self.chrom}: a compiled element spans {longest} bp, beyond crispri.REACH")

    def overlapping(self, start: int, end: int) -> list[str]:
        """Elements sharing at least one base with [start, end), as crispri.DeletionTable searches."""
        out = []
        j = bisect.bisect_left(self.starts, end) - 1
        while j >= 0 and self.spans[j][0] > start - crispri.REACH:
            if self.spans[j][1] > start:
                out.append(self.spans[j][2])
            j -= 1
        return sorted(out)

    def named(self, start: int, end: int, symbol: str) -> list[Link]:
        return [
            lk for el in self.overlapping(start, end) for lk in self.links.get(el, []) if lk.symbol == symbol
        ]

    def measured_names(self, start: int, end: int, symbol: str, cell: str = CELL) -> bool:
        return any(s < end and e > start and g == symbol and c == cell for s, e, g, c in self.measured)


def _kind(rule: Any) -> str:
    return str(getattr(rule.evidence.kind, "value", rule.evidence.kind))


def _symbol(module: Any, gene_id: str) -> str:
    g = module.entities.get(gene_id)
    return str(getattr(g, "symbol", "") or gene_id)


def links_of(module: Any) -> dict[str, list[Link]]:
    """element -> its predicted activates/inhibits rules, from a parsed (ungated) compiled program."""
    out: dict[str, list[Link]] = {}
    for r in module.rules:
        if r.action not in REGULATORY or _kind(r) != "predicted":
            continue
        out.setdefault(r.source, []).append(
            Link(
                r.source,
                r.target,
                _symbol(module, r.target),
                r.action.value,
                float(r.strength),
                r.when.get("cell_type", "") if r.when else "",
            )  # fmt: skip
        )
    return out


def bridged(b: bridge.Bridged) -> tuple[dict[tuple[str, str], float], dict[str, str]]:
    """(element, gene) -> RHO for every mechanism the bridge resolved; subject -> reason for the rest."""
    rho = {(k[0], k[1]): b.observations[k].rho for k in b.strengths}
    why: dict[str, str] = {}
    for u in b.unresolved:
        why.setdefault(u.subject, u.reason)
    return rho, why


def unresolved_reason(link: Link, why: dict[str, str], active_in_k562: bool) -> str:
    if not active_in_k562:
        return "rule_not_active_in_K562"
    return why.get(link.gene_id) or why.get(link.element) or "not_reported_by_the_bridge"


def k562_change(run: dict[str, Any]) -> float:
    """The run's predicted log2 change of the gene on removing the element in K562 (s7_run's value)."""
    return math.log2(run["levels"][CELL] / run["L_held"])


def chromosome_lookup(
    chrom: str,
    text: str,
    pairs: Iterable[crispri.Pair],
    rates: dict[str, float],
    borrowed: float,
    prov: dict | None,
    lent: dict | None,
    clusters: regs.Clusters | None,
    readers: Any,
    results_dir: Any,
) -> ChromLookup:
    """Every stage's answer on one chromosome, for the rules the given pairs reach (no label is read)."""
    from genomeos.lang.parser import parse

    look = ChromLookup(chrom)
    elements = {e: (s, t) for e, c, s, t in cp.program_elements(text)}
    pred, removed = predicted_text(text)
    look.counts["experimental_rule_lines_removed"] = removed

    # s2, s3 and s4a: the predicted program as compiled
    module = parse(pred, f"noncoding_{chrom}")
    # the experimental layer, recorded for the descriptive and never read by a stage
    cells = rgate.compiled_cells(text)
    for src, gene, kind in rgate.rule_targets(text):
        if kind == "experimental" and src in elements:
            s, t = elements[src]
            look.measured.append((s, t, _symbol(module, gene), cells.get((src, gene), "")))
    look.links = links_of(module)
    look.index({e: elements[e] for e in look.links})
    reached: set[tuple[str, str]] = set()
    for p in pairs:
        for lk in look.named(p.start, p.end, p.gene):
            reached.add(lk.key)
    look.counts["rules_reached"] = len(reached)
    b = bridge.parameterize(
        module, CONTEXT, build=False, rates=rates, borrowed_rate=borrowed, rate_provenance=prov or {},
        borrowed_provenance=lent or {},
    )  # fmt: skip
    rho, why = bridged(b)
    by_key = {lk.key: lk for links in look.links.values() for lk in links}
    for key in reached:
        if key in rho:
            look.bridge_compiled[key] = rho[key]
        else:
            look.why_compiled[key] = unresolved_reason(by_key[key], why, by_key[key].cell == CELL)
    del module

    # s4, s5 and s7: the gated program (reader cells)
    gated, info = cp.gate_program(
        pred, chrom, results_dir, selection=cp.SELECTION_READER_CELLS, readers=readers
    )
    look.counts["gated_rules"] = sum(info["rules_gated"].values())
    module = parse(gated, f"noncoding_{chrom}_reader")
    (rd,) = [r for r in module.readers if r.cell_type == CELL]
    b = bridge.parameterize(
        module, CONTEXT, build=False, rates=rates, borrowed_rate=borrowed, rate_provenance=prov or {},
        borrowed_provenance=lent or {},
    )  # fmt: skip
    rho, why = bridged(b)
    for key in sorted(reached):
        look.reader[key[0]] = rd.states.get(key[0], "unknown")
        if key not in rho:
            look.why_reader[key] = unresolved_reason(by_key[key], why, True)
            continue
        look.bridge_reader[key] = rho[key]
        out = rgate.run_gene(rgate.gated_gene_module(module, key[1]), key[1], rates, borrowed, prov, lent)
        if out is None:
            look.counts["run_returned_none"] += 1
            continue
        look.runs[key] = {
            "state": out["states"][CELL],
            "levels_k562": out["levels"][CELL],
            "L_held": out["L_held"],
            "L_open": out["L_open"],
            "rho": out["rho"],
            "change": k562_change(out),
            "gate_mismatch": list(out["gate_mismatch"]),
        }
    del module

    # s6: the regulator assignment in K562, for the elements the reached rules act through
    for el in sorted({k[0] for k in reached}):
        s, t = elements[el]
        if clusters is None:
            look.factors[el] = []
            continue
        occ = regs.element_occupancy(clusters, chrom, s, t)
        lk = look.links[el][0]
        rule = regs.CompiledRule(
            f"noncoding_{chrom}", el, lk.gene_id, lk.action, CELL, "predicted", chrom, s, t
        )
        cands = regs.candidates(rule, occ, frozenset(), {})
        look.factors[el] = sorted(c["factor"] for c in cands if c["kind"] in regs.SOURCE_KINDS)
    return look


# ==== one pair through the stages =========================================================================


def _drop(rho: float) -> float:
    return max(0.0, -math.log2(rho))


def trace_pair(p: crispri.Pair, look: ChromLookup | None) -> tuple[dict[str, float], dict[str, str]]:
    """The two columns of every stage for one pair, and why each stage that drops it did.

    Reads the pair's anchor columns (crispri.annotate must have run) and the chromosome's lookup; never
    its label."""
    f: dict[str, float] = {
        "s1_anchor_kept": p.features["top_target"],
        "s1_anchor_drop": p.features["deletion_drop"],
    }
    why: dict[str, str] = {}
    overlap = look.overlapping(p.start, p.end) if look else []
    named = look.named(p.start, p.end, p.gene) if look else []
    f["s2_target_kept"] = float(bool(named))
    f["s2_target_drop"] = p.features["deletion_drop"] if named else 0.0
    if not named:
        why["s2_target"] = (
            "no compiled element overlaps" if not overlap else "the overlapping elements name another gene"
        )
    f["s3_compiled_value_kept"] = f["s2_target_kept"]
    f["s3_compiled_value_drop"] = max((lk.strength for lk in named if lk.action == "activates"), default=0.0)

    a = [look.bridge_compiled[lk.key] for lk in named if look and lk.key in look.bridge_compiled]
    f["s4a_bridge_as_compiled_kept"] = float(bool(a))
    f["s4a_bridge_as_compiled_drop"] = max((_drop(r) for r in a), default=0.0)
    if named and not a:
        why["s4a_bridge_as_compiled"] = sorted(look.why_compiled[lk.key] for lk in named)[0]

    b = [lk for lk in named if lk.key in look.bridge_reader] if look else []
    f["s4_bridge_kept"] = float(bool(b))
    f["s4_bridge_drop"] = max((_drop(look.bridge_reader[lk.key]) for lk in b), default=0.0)
    if named and not b:
        why["s4_bridge"] = sorted(look.why_reader[lk.key] for lk in named)[0]

    g = [lk for lk in b if look.reader.get(lk.element) == "open"]
    f["s5_reader_kept"] = float(bool(g))
    f["s5_reader_drop"] = max((_drop(look.bridge_reader[lk.key]) for lk in g), default=0.0)
    if b and not g:
        why["s5_reader"] = "reader_" + sorted(look.reader.get(lk.element, "unknown") for lk in b)[0]

    h = [lk for lk in g if look.factors.get(lk.element)]
    f["s6_regulator_kept"] = float(bool(h))
    f["s6_regulator_drop"] = max((_drop(look.bridge_reader[lk.key]) for lk in h), default=0.0)
    if g and not h:
        why["s6_regulator"] = "no factor cluster names K562 at the element"

    runs = [look.runs[lk.key] for lk in b if lk.key in look.runs]
    f["s7_run_kept"] = float(any(r["state"] == "open" for r in runs))
    f["s7_run_drop"] = max((max(0.0, r["change"]) for r in runs if r["state"] == "open"), default=0.0)
    if b and not f["s7_run_kept"]:
        why["s7_run"] = "held in K562"
    return f, why


def columns(stage: str) -> tuple[str, ...]:
    if stage == "s1_anchor":
        return ANCHOR_COLUMNS
    return (*BASE_COLUMNS, f"{stage}_kept", f"{stage}_drop")


# ==== the harness, before any stage beyond the anchor is read =============================================


def harness(lookups: dict[str, ChromLookup], pairs: Sequence[crispri.Pair]) -> list[str]:
    """Every predictor-side check of STOP_RULE that reads no label; an empty list passes."""
    problems = []
    for chrom, look in sorted(lookups.items()):
        if look.counts.get("run_returned_none"):
            problems.append(f"{chrom}: {look.counts['run_returned_none']} resolved rules returned no run")
        for key, run in sorted(look.runs.items()):
            expected = -math.log2(run["rho"]) if run["state"] == "open" else 0.0
            if abs(run["change"] - expected) > CLOSED_FORM_TOLERANCE:
                problems.append(f"{chrom} {key}: run {run['change']:.9g} against closed form {expected:.9g}")
            if run["state"] != look.reader.get(key[0]):
                problems.append(
                    f"{chrom} {key}: run state {run['state']} against reader {look.reader.get(key[0])}"
                )
            if run["gate_mismatch"]:
                problems.append(
                    f"{chrom} {key}: gated level differs from the clamped one in {run['gate_mismatch']}"
                )
    for stage, spec in STAGES.items():
        parent = spec["follows"]
        if not parent or parent == "s1_anchor":
            continue  # the compiled target is not a subset of the compact table's top target by construction
        stray = sum(1 for p in pairs if p.features[f"{stage}_kept"] and not p.features[f"{parent}_kept"])
        if stray:
            problems.append(f"{stage} keeps {stray} pairs that {parent} dropped")
    return problems


# ==== the readings ========================================================================================


def fit_scores(fit_rows: Sequence[crispri.Pair], pop: Sequence[crispri.Pair], cols: tuple[str, ...]) -> tuple:
    w = crispri.logistic_fit(crispri.matrix(list(fit_rows), cols), [p.regulated for p in fit_rows])
    return w, crispri.logistic_score(w, crispri.matrix(list(pop), cols))


def decide(kept: int, ci: Sequence[float] | None, vs_anchor: Sequence[float] | None) -> str:
    if kept == 0:
        return "NONE_KEPT"
    if ci is None or vs_anchor is None:
        raise ValueError(
            "an interval is missing: the population spans fewer chromosomes than an interval needs"
        )
    lo, hi = ci
    dhi = vs_anchor[1]
    if hi < 0:
        return "WORSE"
    if lo > 0:
        return "WEAKENED" if dhi < 0 else "CARRIES"
    return "LOST" if dhi < 0 else "UNDECIDED"


def reading(kind: str, **numbers: Any) -> str:
    return OUTCOMES[kind].format(**numbers)


def _signed(x: float) -> str:
    return f"{x:+.4f}"


def stage_numbers(row: dict[str, Any]) -> dict[str, Any]:
    g, d = row["gain_over_activity_distance"], row["against_anchor"]
    return {
        "g": _signed(g["gain"]),
        "lo": _signed(g["ci95"][0]),
        "hi": _signed(g["ci95"][1]),
        "d": _signed(d["gain"]),
        "dlo": _signed(d["ci95"][0]),
        "dhi": _signed(d["ci95"][1]),
        "n": f"{row['population']['pairs']:,}",
        "k": f"{row['population']['kept']:,}",
        "kp": f"{row['population']['kept_regulated']:,}",
        "kf": f"{row['fit_rows']['kept']:,}",
        "nf": f"{row['fit_rows']['pairs']:,}",
    }


def kept_counts(rows: Sequence[crispri.Pair], stage: str) -> dict[str, int]:
    kept = [p for p in rows if p.features[f"{stage}_kept"]]
    return {
        "pairs": len(rows),
        "kept": len(kept),
        "kept_regulated": sum(p.regulated for p in kept),
        "kept_with_value_above_0": sum(1 for p in kept if p.features[f"{stage}_drop"] > 0),
        "kept_regulated_with_value_above_0": sum(
            1 for p in kept if p.regulated and p.features[f"{stage}_drop"] > 0
        ),
    }


def anchor_reading(fit_rows: Sequence[crispri.Pair], pop: Sequence[crispri.Pair]) -> dict[str, Any]:
    """The anchor recomputed by its own recipe, and whether it equals the committed figures."""
    _, base = fit_scores(fit_rows, pop, BASE_COLUMNS)
    w, anchor = fit_scores(fit_rows, pop, ANCHOR_COLUMNS)
    g = crispri.gain_interval(anchor, base, list(pop))
    row = {
        "pairs": len(pop),
        "positives": sum(p.regulated for p in pop),
        "uncovered_pairs": sum(not p.covered for p in pop),
        **g,
    }
    same = all(row[k] == ANCHOR[k] for k in ("pairs", "positives", "uncovered_pairs", "gain", "ci95"))
    return {
        "recomputed": row,
        "committed": ANCHOR,
        "reproduced": same,
        "weights": w,
        "scores": anchor,
        "base": base,
    }


def score(
    fit_rows: Sequence[crispri.Pair], pop: Sequence[crispri.Pair], anchor: dict[str, Any]
) -> dict[str, Any]:
    """Every stage's reading (METRIC), its outcome (DECISION) and the chain's verdict."""
    pop = list(pop)
    scores = {"s1_anchor": anchor["scores"]}
    weights = {"s1_anchor": anchor["weights"]}
    for stage in STAGES:
        if stage != "s1_anchor":
            weights[stage], scores[stage] = fit_scores(fit_rows, pop, columns(stage))
    stages: dict[str, Any] = {}
    for stage, spec in STAGES.items():
        row: dict[str, Any] = {
            "name": spec["name"],
            "follows": spec["follows"],
            "columns": list(columns(stage)),
            "weights": [round(x, 6) for x in weights[stage]],
            "population": kept_counts(pop, stage),
            "fit_rows": kept_counts(fit_rows, stage),
            "gain_over_activity_distance": crispri.gain_interval(scores[stage], anchor["base"], pop),
        }
        if stage == "s1_anchor":
            row["outcome"] = "ANCHOR"
            num = stage_numbers({**row, "against_anchor": {"gain": 0.0, "ci95": [0.0, 0.0]}})
            row["reading"] = reading("ANCHOR", **num)
        else:
            row["against_anchor"] = crispri.gain_interval(scores[stage], scores["s1_anchor"], pop)
            row["against_the_stage_it_follows"] = crispri.gain_interval(
                scores[stage], scores[spec["follows"]], pop
            )
            row["outcome"] = decide(
                row["population"]["kept"],
                row["gain_over_activity_distance"]["ci95"],
                row["against_anchor"]["ci95"],
            )
            row["reading"] = reading(row["outcome"], **stage_numbers(row))
        stages[stage] = row
    verdict, sentence, at = chain_verdict(stages)
    return {"stages": stages, "chain": list(CHAIN), "branches": list(BRANCHES), "outcome": verdict,
            "reading": sentence, "stage_where_it_stops": at}  # fmt: skip


def chain_verdict(stages: dict[str, Any]) -> tuple[str, str, str | None]:
    """(outcome, its sentence, the stage it names) for the chain, walked in order (DECISION)."""
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


def dropped_because(whys: Sequence[dict[str, str]], stage: str) -> dict[str, int]:
    return dict(Counter(w[stage] for w in whys if stage in w).most_common())
