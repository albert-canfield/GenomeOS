# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""M2d (ROADMAP §5 item 7e, lane-s4joint, 2026-10-10): signaltrace's bridge stage re-run with the joint fit.

M1c (`signaltrace.py`, registration 9637584, result 1103719) traced the 1,918 held-out ENCODE CRISPRi
pairs in K562 through each compiled stage with the anchor's metric. Target choice kept 39 pairs and a
gain over activity and distance (+0.0792 [0.0282, 0.1562]); the bridge in K562 kept none, because every
gene those pairs reach carries two or more compiled rules and `bridge.parameterize` refuses such a gene
(`not_identifiable`). Registered reading: LOST_AT s4_bridge, first weakened at s2_target.

lane-combine2 (`combine2.py`, bc9fc7a .. 948cfc9) showed that under the reader's hold a held rule stays in
the mean at strength 0, so the mean is linear in the element terms and a JOINT fit of a gene's k
single-removal observations identifies it; it registered two forms (the runtime's mean; the promoter-
competition form, Bothma 2015, THETA = 0.5 declared).

This module asks the one question left between the two: with the bridge replaced by combine2's joint fit,
does the measured signal survive the bridge once the genes are fittable? Same pairs (signaltrace's ids,
reused by sha256), same fit rows, same metric, same bootstrap, same decision and the same outcome
sentences; the forms and their parameters exactly as combine2 committed them; nothing new is fitted. It
re-implements no stage: the target is signaltrace's, the joint fit is `combine2.classify` and
`combine2.fit`, the run is `combine2.run_gene_form` and the runtime's own level with one more element held,
the regulator lookup is signaltrace's. Its own code is the join and the readings.

In the order written and committed: the code (`joint_lookup`, `trace_joint`, `score_forms`), then the
registration (`QUESTION` .. `DESCRIPTIVES`, data/results/s4joint_registration.json), then one run.

No model request, no key, no money: every input is on this machine.
"""

from __future__ import annotations

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
from genomeos.attribution import signaltrace as st
from genomeos.runtime import grn

REGISTERED = "2026-10-10"
LANE = "lane-s4joint"
CELL = st.CELL
K562_INDEX = c2.CELLS.index(CELL)

# ---- the two forms, exactly as combine2 committed them (data/results/combine2_registration.json) --------
FORMS = c2.FORMS
THETA = c2.THETA
FIT = c2.FIT
FORM_SOURCES = c2.FORM_SOURCES
FORM_WORDS = c2.FORM_WORDS
#: combine2_registration.json field -> combine2 constant; the run refuses if any differs (nothing refitted)
COMBINE2_BOUND = {
    "forms": "FORMS",
    "theta": "THETA",
    "fit": "FIT",
    "form_sources": "FORM_SOURCES",
    "structure": "STRUCTURE",
    "selection": "SELECTION",
    "cells": "CELLS",
    "runtime": "RUNTIME",
    "closed_form_tolerance": "CLOSED_FORM_TOLERANCE",
    "harness": "HARNESS",
    "form_words": "FORM_WORDS",
}

# ---- the registered re-run ------------------------------------------------------------------------------
QUESTION = (
    "M1c (signaltrace) followed the held-out ENCODE CRISPRi pairs in K562 through the compiled layer and"
    " lost the measured signal at the bridge: target choice kept 39 of 1,918 pairs with a gain over"
    " activity and distance (+0.0792 [95% CI 0.0282, 0.1562]), and the bridge in K562 kept none, because"
    " every gene those pairs reach carries two or more compiled rules, which bridge.parameterize cannot"
    " identify. With the bridge replaced by combine2's joint fit of a gene's elements - once under the"
    " runtime's mean, once under the promoter-competition form, each as combine2 committed it - and the"
    " reader gate, the regulator assignment and the run following it as in M1c, does the agreement with"
    " the measured CRISPRi outcome survive the bridge once the genes are fittable?"
)
POPULATION = (
    "signaltrace's population, reused and not re-chosen: the 1,918 K562 pairs of the ENCODE held-out file,"
    " each identified as chrom:chromStart-chromEnd:measuredGeneSymbol, in file order, whose newline-joined"
    " list has the sha256 signaltrace's registration (data/results/signaltrace_registration.json) pins;"
    " the run refuses if the pairs it reads differ from that list. The label is the file's Regulated"
    " column, the anchor's own. Every stage is read on all 1,918: a pair a stage drops stays with 0 in that"
    " stage's two columns. The fit rows are the anchor's and signaltrace's: the K562 training pairs"
    " crispri.annotate marks covered, unchanged for every stage"
)
PREDICTED_ONLY = st.PREDICTED_ONLY
JOINT_FIT = (
    "per chromosome, the program with its experimental rule lines removed (PREDICTED_ONLY) is gated as"
    " signaltrace's s4_bridge gated it (compile.gate_program, selection reader_cells: every predicted rule"
    " is offered in K562 whatever its compiled cell) and parsed. Every gene that a rule reached by the"
    " traced pairs (fit rows and population, signaltrace.chromosome_lookup's reached rules) targets is read"
    " by combine2.gene_rules on that module: all of the gene's rules in that program, each with its removal"
    " observation (bridge.observation_of) and its reader states. combine2.classify admits it (every"
    " observation readable, none exactly 1, each rule's action agreeing with its observation's sign, at"
    " least one activating) and combine2.fit gives each form's joint fit; a gene is runnable under a form"
    " as combine2.FIT['runnable'] states. The rate T is the gene's rate in data/results/gene_rates.json,"
    " else that table's genome median (the borrowed tier, as combine2's census and signaltrace's bridge"
    " read it); no fold depends on T. A runnable gene is run by combine2.run_gene_form (its thirteen gated"
    " reader cells and its single-removal reproduction, combine2.HARNESS) and, in K562, once more for each"
    " of its elements open there with that element clamped absent (bridge.REMOVED), which under both forms"
    " equals holding the rule. The joint fit applies to a gene of any rule count; signaltrace found every"
    " one of the 377 reached rules not_identifiable in K562, so every gene read here has two or more"
)
#: the stages this registration adds, per form; s1_anchor, s2_target and s3_compiled_value are signaltrace's
STAGE_TEMPLATES: dict[str, dict[str, str]] = {
    "s4_joint": {
        "name": "the strength combine2's joint fit gives in K562, {words}",
        "follows": "s3_compiled_value",
        "kept": (
            "a rule kept by s3_compiled_value (an overlapping element's predicted rule naming the measured"
            " gene) whose gene the joint fit makes runnable under {form} (JOINT_FIT), the rule among the"
            " gene's fitted rules"
        ),
        "value": (
            "max(0, -log2 RHO) of the rule's removal observation, which the fitted run reproduces at the"
            " state where every element of the gene is present (combine2.HARNESS (2)); the largest over the"
            " kept rules. It is signaltrace's s4_bridge value, now carried for every rule of a fitted gene"
        ),
    },
    "s5_reader": {
        "name": "the reader gate in K562 after the joint fit, {words}",
        "follows": "s4_joint_{form}",
        "kept": (
            "a rule kept by s4_joint_{form} whose element the gated program's K562 reader block states open"
            " (an ENCODE DNase-seq narrowPeak overlaps it); closed and unknown hold the rule"
        ),
        "value": "as s4_joint_{form}, over the open rules",
    },
    "s6_regulator": {
        "name": "the regulator assignment in K562 after the joint fit, {words}",
        "follows": "s5_reader_{form}",
        "kept": (
            "a rule kept by s5_reader_{form} whose element signaltrace's lookup gives a source factor in K562"
            " (regulators.candidates with the rule's cell set to K562: an ENCODE TF ChIP-seq peak cluster"
            " naming K562 overlaps the element)"
        ),
        "value": "as s5_reader_{form}, over those rules",
    },
    "s7_run": {
        "name": "the jointly fitted run's predicted change in K562, {words}",
        "follows": "s5_reader_{form}",
        "kept": (
            "a rule kept by s5_reader_{form}, run: the gene's fitted program (combine2.fitted_module on"
            " readergate.gated_gene_module) run by runtime.grn.NetworkRuntime with"
            " activator_combination={form} in K562, every element clamped present so that the reader holds"
            " the closed ones"
        ),
        "value": (
            "max(0, log2(L_K562 / L_held)): the gene's level in K562 with the reader applied, over its level"
            " with this rule held as well (its element clamped absent); the largest over the run rules. The"
            " gene's other elements are open or held as K562's reader states them, so this is the rule's"
            " own observation only when no other element of the gene changes the level it acts on"
        ),
    },
}


def stage_name(template: str, form: str) -> str:
    return f"{template}_{form}"


def form_stages(form: str) -> dict[str, dict[str, str]]:
    """The four stages of one form, their words filled from the templates."""
    words = FORM_WORDS[form]
    return {
        stage_name(t, form): {k: v.format(form=form, words=words) for k, v in spec.items()}  # fmt: skip
        for t, spec in STAGE_TEMPLATES.items()
    }


STAGES: dict[str, dict[str, str]] = {
    **{s: st.STAGES[s] for s in ("s1_anchor", "s2_target", "s3_compiled_value")},
    **{s: spec for form in FORMS for s, spec in form_stages(form).items()},
}
COMPARATORS = ("s2_target", "s3_compiled_value")
CHAINS = {
    form: ("s1_anchor", *COMPARATORS, *(stage_name(t, form) for t in ("s4_joint", "s5_reader", "s7_run")))
    for form in FORMS
}
BRANCHES = {form: (stage_name("s6_regulator", form),) for form in FORMS}
JOINT_STAGES = tuple(s for form in FORMS for s in form_stages(form))
ANCHOR = st.ANCHOR
COMPARATOR = (
    "the anchor and signaltrace's s2_target and s3_compiled_value, recomputed by signaltrace's own code on"
    " the same pairs and fit rows; the run stops (HARNESS_DEFECT) unless the anchor equals its committed"
    " arm-1 figures and s2_target and s3_compiled_value equal signaltrace's committed readings exactly"
    " (data/results/signaltrace_result.json: kept counts, gain over activity and distance, interval, and"
    " against the anchor with its interval). Every joint stage is read against the anchor and against the"
    " stage it follows by the same paired resampling; signaltrace's s4_bridge reading (NONE_KEPT) and"
    " chain reading are quoted beside it, verbatim"
)
BASE_COLUMNS = st.BASE_COLUMNS
METRIC = st.METRIC
BOOTSTRAP = st.BOOTSTRAP
DECISION = st.DECISION
PER_FORM = (
    "signaltrace's DECISION is applied to each form F on its own chain (s2_target, s3_compiled_value,"
    " s4_joint_F, s5_reader_F, s7_run_F), with s4_joint_F in s4_bridge's place, s5_reader_F in"
    " s5_reader's and s7_run_F in s7_run's, and s6_regulator_F read beside it on s5_reader_F's pairs as"
    " signaltrace reads s6_regulator. The two forms are read separately, each at the anchor's 95%"
    " interval; neither reading is pooled with the other, and no difference between the forms is tested"
)
#: s4_joint_F's stage outcome -> the bridge answer for F
BRIDGE_OF = {
    "CARRIES": "SURVIVES",
    "WEAKENED": "SURVIVES",
    "NONE_KEPT": "LOST",
    "LOST": "LOST",
    "WORSE": "LOST",
    "UNDECIDED": "UNDECIDED",
}
BRIDGE_DECISION = (
    "per form, from s4_joint_F's outcome under DECISION: SURVIVES if CARRIES or WEAKENED (the 95% interval"
    " of its gain over activity and distance lies above zero, so the chain passes the bridge), LOST if"
    " NONE_KEPT, LOST or WORSE (signaltrace's LOSS), UNDECIDED if UNDECIDED. 'Survives' means exactly"
    " that; whether the gain is below the anchor's is the stage reading's own clause, kept beside it"
)
CLOSED_FORM_TOLERANCE = c2.CLOSED_FORM_TOLERANCE
STOP_RULE = (
    "one run, from this registration's committed code, on the cached files the registration pins by"
    " sha256, signaltrace's and combine2's registrations among them. No stage, column, population, fit"
    " row, form, parameter, metric or wording is changed after any joint stage is read; a defect found"
    " afterwards is reported beside the registered readings and corrected only under a new registration."
    " The run stops with HARNESS_DEFECT, before any joint stage is read, if: the population differs from"
    " signaltrace's pinned list; a field of combine2's registration in COMBINE2_BOUND differs from"
    " combine2.py; signaltrace's own stage columns, recomputed, differ on any population pair from its"
    " committed pairs table; signaltrace's own harness fails; a reached gene is not read from the gated"
    " program, or its K562 reader state differs from signaltrace's lookup; any gene-form run misses"
    " combine2.HARNESS or a K562 run with one more element held differs from the form's closed form by more"
    " than CLOSED_FORM_TOLERANCE in log2; a joint stage keeps a pair the stage it follows dropped; the anchor"
    " does not reproduce the committed arm-1 figures; or s2_target or s3_compiled_value does not reproduce"
    " signaltrace's committed reading"
)
OUTCOMES = st.OUTCOMES
CHAIN_OUTCOMES = st.CHAIN_OUTCOMES
BEFORE_CLAUSE = st.BEFORE_CLAUSE
LOSS = st.LOSS
BRIDGE_OUTCOMES = {
    "SURVIVES": (
        "with the genes fitted jointly under {form}, the measured signal survives the bridge: at {stage}"
        " ({name}), {stage_reading}"
    ),
    "LOST": (
        "with the genes fitted jointly under {form}, the measured signal does not survive the bridge: at"
        " {stage} ({name}), {stage_reading}"
    ),
    "UNDECIDED": (
        "with the genes fitted jointly under {form}, whether the measured signal survives the bridge is not"
        " decided: at {stage} ({name}), {stage_reading}"
    ),
}
FORM_READING = "{bridge}. Along the chain under {form}, {chain}."
LIMITATIONS = (
    *st.LIMITATIONS,
    c2.LIMITATIONS[1],
    c2.LIMITATIONS[2],
    c2.LIMITATIONS[3],
    "a pair whose gene no form can fit is dropped at s4_joint with 0 in the stage's two columns, as the"
    " anchor treats an uncovered pair; the runnable genes are those whose summed single-removal effects"
    " leave a positive basal rate, a selected subset that differs between the forms",
    "at s4_joint the value is the rule's own observation under both forms, because each fit reproduces"
    " every single-removal observation; s4_joint_mean and s4_joint_competition differ only in the pairs"
    " they keep. The forms differ in value only at s7_run, where the gene's other elements are open or held"
    " as K562's reader states them",
    "two forms are read, each at the anchor's 95% interval with no adjustment for reading two; each form's"
    " reading stands alone, nothing here pools them, and no difference between the forms is tested",
    "this lane read signaltrace's committed readings and combine2's census and result before registering"
    " (exploration_preceded_this): the comparators are known, the joint stages on the real pairs are not",
)
DESCRIPTIVES = (
    "per form: the reached genes by whether the joint fit runs them and, where not, why (combine2's"
    " reasons); the quantiles of q, the element-driven share of intact transcription, over the reached"
    " runnable genes",
    "per stage: the pairs kept of the 1,918 and of the fit rows, how many of them regulated, how many with"
    " a value above 0, the reasons the pairs it drops were dropped, the fitted weights of its model, and its"
    " paired interval against the stage it follows",
    "the run checks: gene-form runs, K562 hold runs, the largest difference from the closed form in log2",
    "signaltrace's committed s4_bridge and chain readings, quoted verbatim for context",
)
#: the registration fields the run refuses to differ from: registration field -> constant of this module
BOUND = {
    "question": "QUESTION",
    "population": "POPULATION",
    "predicted_only": "PREDICTED_ONLY",
    "cell": "CELL",
    "forms": "FORMS",
    "theta": "THETA",
    "fit": "FIT",
    "form_sources": "FORM_SOURCES",
    "form_words": "FORM_WORDS",
    "combine2_bound": "COMBINE2_BOUND",
    "joint_fit": "JOINT_FIT",
    "stages": "STAGES",
    "chains": "CHAINS",
    "branches": "BRANCHES",
    "anchor": "ANCHOR",
    "comparator": "COMPARATOR",
    "base_columns": "BASE_COLUMNS",
    "metric": "METRIC",
    "bootstrap": "BOOTSTRAP",
    "decision": "DECISION",
    "per_form": "PER_FORM",
    "bridge_of": "BRIDGE_OF",
    "bridge_decision": "BRIDGE_DECISION",
    "closed_form_tolerance": "CLOSED_FORM_TOLERANCE",
    "stop_rule": "STOP_RULE",
    "outcomes": "OUTCOMES",
    "chain_outcomes": "CHAIN_OUTCOMES",
    "before_clause": "BEFORE_CLAUSE",
    "loss": "LOSS",
    "bridge_outcomes": "BRIDGE_OUTCOMES",
    "form_reading": "FORM_READING",
    "limitations": "LIMITATIONS",
    "descriptives_not_verdicts": "DESCRIPTIVES",
}


def registered_fields() -> dict[str, Any]:
    """Every BOUND constant as JSON stores it, so a registration and this module compare equal or not."""
    return {f: json.loads(json.dumps(globals()[c])) for f, c in BOUND.items()}


def combine2_mismatch(reg: dict[str, Any]) -> list[str]:
    """The COMBINE2_BOUND fields where combine2's committed registration and combine2.py differ."""
    return [
        k for k, const in COMBINE2_BOUND.items() if reg.get(k) != json.loads(json.dumps(getattr(c2, const)))
    ]


# ==== the joint fit, one chromosome at a time =============================================================


@dataclass
class GeneFit:
    """One reached gene as combine2 reads, classifies, fits and runs it in the gated program."""

    gene: str
    symbol: str
    elements: list[str]
    rho: dict[str, float]
    k562: dict[str, str]  # element -> K562 reader state
    excluded: str  # combine2.classify's reason, "" when admitted
    forms: dict[str, dict[str, Any]] = field(default_factory=dict)


@dataclass
class JointLookup:
    """The joint fit's answer on one chromosome, for the genes the traced pairs reach."""

    chrom: str
    genes: dict[str, GeneFit] = field(default_factory=dict)
    reached_genes: set[str] = field(default_factory=set)
    counts: Counter = field(default_factory=Counter)

    def fitted(self, form: str, link: st.Link) -> GeneFit | None:
        g = self.genes.get(link.gene_id)
        if g is None or g.excluded or not g.forms[form]["ok"] or link.element not in g.rho:
            return None
        return g

    def reason(self, form: str, link: st.Link) -> str:
        g = self.genes.get(link.gene_id)
        if g is None:
            return "gene_not_read_from_the_gated_program"
        if g.excluded:
            return g.excluded
        if not g.forms[form]["ok"]:
            return g.forms[form]["reason"]
        return "rule_not_among_the_genes_fitted_rules"


def rate_of(gene: str, rates: dict[str, float], borrowed: float, prov: dict | None, lent: dict | None):
    """(T, its label): the gene's measured rate, else the borrowed median (combine2's census, verbatim)."""
    if gene in rates:
        rate, label = rates[gene], {**(prov or {}).get(gene, {}), "tier": "measured"}
    else:
        rate, label = borrowed, {**(lent or {}), "tier": "borrowed_median"}
    return rate, {k: str(label.get(k, "unstated")) for k in ("source", "species_cell", "tier")}


def run_k562(small: Any, gr: c2.GeneRules, f: c2.Fit, rate: float, label: dict[str, str]) -> dict[str, Any]:
    """combine2's run of one gene under one form, and in K562 one more run per open element held."""
    got = c2.run_gene_form(small, gr, f, rate, label)
    m = c2.fitted_module(small, gr.gene, f, rate, label)
    present = dict.fromkeys(gr.elements, bridge.INTACT)
    delta = math.log(2) / grn.DEFAULTS["mrna_half_life"]
    open_k = gr.open_in(K562_INDEX)
    level = got["levels"][K562_INDEX]
    change: dict[str, float] = {}
    miss = list(got["harness_miss"])
    worst = 0.0
    for el in gr.elements:
        if not open_k[el]:
            continue
        held = c2._run_level(m, gr.gene, {"cell_type": CELL}, {**present, el: bridge.REMOVED}, f.form)
        want = rate * c2.level(f, {**open_k, el: False}) / delta
        off = abs(math.log2(held / want))
        worst = max(worst, off)
        if not off <= CLOSED_FORM_TOLERANCE:
            miss.append(f"K562 hold {el}: run {held:.12g}, closed form {want:.12g}")
        change[el] = math.log2(level / held)
    worst = max([worst, *(abs(a - b) for a, b in zip(got["s"], got["closed_s"], strict=True))])
    return {
        "ok": True,
        "reason": "",
        "q": f.q,
        "tier": label["tier"],
        "change": change,
        "harness_miss": miss,
        "max_abs_run_minus_closed_form_log2": worst,
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
) -> JointLookup:
    """Every reached gene of one chromosome fitted and run under each form (no label is read)."""
    from genomeos.lang.parser import parse

    out = JointLookup(chrom)
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
        fits, why = c2.classify(gr)
        g = GeneFit(
            gr.gene,
            gr.symbol,
            list(gr.elements),
            dict(gr.rho),
            {el: gr.states[el][K562_INDEX] for el in gr.elements},
            why,
        )
        small = None
        for form in FORMS:
            if why:
                g.forms[form] = {"ok": False, "reason": why}
                continue
            f = fits[form]
            if not f.ok:
                g.forms[form] = {"ok": False, "reason": f.reason, "q": f.q}
                continue
            if small is None:
                small = rgate.gated_gene_module(module, gr.gene)
            rate, label = rate_of(gr.gene, rates, borrowed, prov, lent)
            g.forms[form] = run_k562(small, gr, f, rate, label)
            out.counts[f"runs_{form}"] += 1
        out.genes[gr.gene] = g
    del module
    return out


# ==== one pair through the joint stages ===================================================================


def trace_joint(
    p: crispri.Pair, look: st.ChromLookup | None, joint: JointLookup | None
) -> tuple[dict[str, float], dict[str, str]]:
    """The two columns of every joint stage for one pair, and why each stage that drops it did.

    Reads the chromosome's two lookups and never the pair's label."""
    f: dict[str, float] = {}
    why: dict[str, str] = {}
    named = look.named(p.start, p.end, p.gene) if look else []
    for form in FORMS:
        s4, s5, s6, s7 = (stage_name(t, form) for t in STAGE_TEMPLATES)
        fitted = [(lk, joint.fitted(form, lk)) for lk in named] if joint else []
        b = [(lk, g) for lk, g in fitted if g is not None]
        f[f"{s4}_kept"] = float(bool(b))
        f[f"{s4}_drop"] = max((st._drop(g.rho[lk.element]) for lk, g in b), default=0.0)
        if named and not b:
            why[s4] = sorted(joint.reason(form, lk) for lk in named)[0] if joint else "no_joint_lookup"
        o = [(lk, g) for lk, g in b if g.k562[lk.element] == "open"]
        f[f"{s5}_kept"] = float(bool(o))
        f[f"{s5}_drop"] = max((st._drop(g.rho[lk.element]) for lk, g in o), default=0.0)
        if b and not o:
            why[s5] = "reader_" + sorted(g.k562[lk.element] for lk, g in b)[0]
        h = [(lk, g) for lk, g in o if look.factors.get(lk.element)]
        f[f"{s6}_kept"] = float(bool(h))
        f[f"{s6}_drop"] = max((st._drop(g.rho[lk.element]) for lk, g in h), default=0.0)
        if o and not h:
            why[s6] = "no factor cluster names K562 at the element"
        runs = [g.forms[form]["change"][lk.element] for lk, g in o if lk.element in g.forms[form]["change"]]
        f[f"{s7}_kept"] = float(bool(runs))
        f[f"{s7}_drop"] = max((max(0.0, c) for c in runs), default=0.0)
        if b and not runs:
            why[s7] = "held in K562"
    return f, why


# ==== the harness, before any joint stage is read =========================================================


def harness(
    lookups: dict[str, st.ChromLookup], joints: dict[str, JointLookup], pairs: Sequence[crispri.Pair]
) -> list[str]:
    """Every predictor-side check of STOP_RULE on the joint stages that reads no label; empty passes."""
    problems = []
    for chrom, joint in sorted(joints.items()):
        look = lookups[chrom]
        missing = sorted(joint.reached_genes - set(joint.genes))
        if missing:
            problems.append(f"{chrom}: reached genes not read from the gated program: {missing[:5]}")
        for el, gene in sorted(set(look.bridge_reader) | set(look.why_reader)):
            g = joint.genes.get(gene)
            if g is not None and el in g.k562 and g.k562[el] != look.reader.get(el):
                problems.append(f"{chrom} {el}: joint K562 state {g.k562[el]} against {look.reader.get(el)}")
        for gene, g in sorted(joint.genes.items()):
            for form, fr in g.forms.items():
                if fr.get("harness_miss"):
                    problems.append(f"{chrom} {gene} {form}: {'; '.join(fr['harness_miss'][:3])}")
    for stage in JOINT_STAGES:
        parent = STAGES[stage]["follows"]
        stray = sum(1 for p in pairs if p.features[f"{stage}_kept"] and not p.features[f"{parent}_kept"])
        if stray:
            problems.append(f"{stage} keeps {stray} pairs that {parent} dropped")
    return problems


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


def fit_stages(
    stages: Sequence[str],
    fit_rows: Sequence[crispri.Pair],
    pop: Sequence[crispri.Pair],
    anchor: dict[str, Any],
) -> tuple[dict[str, list[float]], dict[str, list[float]]]:
    weights, scores = {"s1_anchor": anchor["weights"]}, {"s1_anchor": anchor["scores"]}
    for s in stages:
        weights[s], scores[s] = st.fit_scores(fit_rows, list(pop), st.columns(s))
    return weights, scores


def comparators(
    fit_rows: Sequence[crispri.Pair], pop: Sequence[crispri.Pair], anchor: dict[str, Any]
) -> tuple[dict[str, Any], dict[str, list[float]], dict[str, list[float]]]:
    """s2_target and s3_compiled_value recomputed: rows, scores and weights (anchor included)."""
    weights, scores = fit_stages(COMPARATORS, fit_rows, pop, anchor)
    rows = {s: stage_row(s, fit_rows, pop, anchor, scores, weights) for s in COMPARATORS}
    return rows, scores, weights


COMPARED_FIELDS = ("population", "fit_rows")


def comparator_mismatch(rows: dict[str, Any], committed: dict[str, Any]) -> list[str]:
    """Where a recomputed comparator differs from signaltrace's committed stage (signaltrace_result.json)."""
    bad = []
    for s in COMPARATORS:
        mine, theirs = rows[s], committed[s]
        for k in COMPARED_FIELDS:
            if mine[k] != theirs[k]:
                bad.append(f"{s} {k}: {mine[k]} against {theirs[k]}")
        for k in ("gain_over_activity_distance", "against_anchor"):
            for part in ("gain", "ci95"):
                if mine[k][part] != theirs[k][part]:
                    bad.append(f"{s} {k} {part}: {mine[k][part]} against {theirs[k][part]}")
        if mine["outcome"] != theirs["outcome"] or mine["reading"] != theirs["reading"]:
            bad.append(f"{s}: outcome or reading differs")
    return bad


def chain_verdict(stages: dict[str, Any], chain: Sequence[str]) -> tuple[str, str, str | None]:
    """(outcome, its sentence, the stage it names) for one form's chain: signaltrace.chain_verdict's walk."""
    first_weak = None
    for stage in chain[1:]:
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
    run = stages[chain[-1]]["reading"]
    if first_weak:
        text = CHAIN_OUTCOMES["WEAKENED_NOT_LOST"].format(
            stage=first_weak, name=STAGES[first_weak]["name"], run_reading=run
        )
        return "WEAKENED_NOT_LOST", text, first_weak
    return "CARRIED_TO_THE_RUN", CHAIN_OUTCOMES["CARRIED_TO_THE_RUN"].format(run_reading=run), None


def bridge_verdict(stages: dict[str, Any], form: str) -> tuple[str, str]:
    stage = stage_name("s4_joint", form)
    row = stages[stage]
    kind = BRIDGE_OF[row["outcome"]]
    text = BRIDGE_OUTCOMES[kind].format(
        form=FORM_WORDS[form], stage=stage, name=STAGES[stage]["name"], stage_reading=row["reading"]
    )
    return kind, text


def score_forms(
    fit_rows: Sequence[crispri.Pair],
    pop: Sequence[crispri.Pair],
    anchor: dict[str, Any],
    comp_rows: dict[str, Any],
    comp_scores: dict[str, list[float]],
    comp_weights: dict[str, list[float]],
) -> dict[str, Any]:
    """Every joint stage's reading, each form's bridge answer and chain verdict, and the joint reading."""
    weights, scores = fit_stages(JOINT_STAGES, fit_rows, pop, anchor)
    scores.update(comp_scores)
    weights.update(comp_weights)
    stages = dict(comp_rows)
    for s in JOINT_STAGES:
        stages[s] = stage_row(s, fit_rows, pop, anchor, scores, weights)
    per_form: dict[str, Any] = {}
    for form in FORMS:
        bkind, btext = bridge_verdict(stages, form)
        ckind, ctext, at = chain_verdict(stages, CHAINS[form])
        per_form[form] = {
            "bridge": bkind,
            "bridge_reading": btext,
            "chain": ckind,
            "chain_reading": ctext,
            "stage_where_it_stops": at,
            "reading": FORM_READING.format(bridge=btext, form=FORM_WORDS[form], chain=ctext),
        }
    return {
        "stages": stages,
        "per_form": per_form,
        "reading": " ".join(per_form[f]["reading"] for f in FORMS),
    }
