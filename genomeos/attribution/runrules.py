# SPDX-License-Identifier: AGPL-3.0-or-later
"""M1 (ROADMAP §5 item 14, lane-runrules, 2026-10-10): do the compiled rules, RUN, predict measured RNA?

The project compiles its non-coding attributions into 24 BioLang programs (440,589 regulatory rules)
and `attribution/bridge.py` turns one observation per gene and cell into one fitted strength, so that
`runtime/grn.py` can integrate the result (R3: 22,576 pairs on a measured mouse rate, 160,392 in
relative units on the borrowed median). Nobody had asked whether a run of them predicts anything.

This module holds three things, in the order they were written and committed:

1. The census (no outcome read): per cell context, how many compiled rules, mechanisms and gene-cell
   pairs there are, how many the bridge resolves on each rate tier, and why the rest are not; plus the
   structural facts that decide what a run CAN predict (`STRUCTURE`).
2. The registered test (`QUESTION` .. `OUTCOMES`), fixed and committed before the measured RNA was
   opened: the population, the prediction a run makes, the no-rule baseline, the shuffled-rule
   control, the metric, the decision rule, the stop rule and the exact wording of every outcome.
3. The code that runs it (`run_pair`, `score`, `reading`), committed before the registration.

What the census found, and why the registered test is not the on/off test ROADMAP §5 item 14 M1
named. The bridge fits each mechanism so that the run reproduces its own input observation (RHO) in
its own cell, and sets the gene's intact level to its one rate T in every cell it resolves; every
resolved gene has a basal rate > 0, because the split refuses one that does not
(`rate_split_infeasible`). So a run predicts "on" for every resolved gene in every cell, and the same
level in every cell: neither on/off nor level can differ from a no-rule baseline, by construction.
No compiled rule's source is a gene, so a run integrates independent one-gene systems and nothing
emerges from the network. What a run does compute that is not its own input restated in its own cell
is the gene's level where the rule's element is present against where it is absent, read across
cells (`READING`). That is the one prediction about measured RNA a run of the compiled programs
makes today, and it is what is registered.

No model request, no network, no money: every input is on this machine.
"""

from __future__ import annotations

import copy
import gzip
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from genomeos.attribution import bridge
from genomeos.ir import Action, Module, Rule
from genomeos.runtime.grn import NetworkRuntime

REGISTERED = "2026-10-10"
LANE = "lane-runrules"

# ---- the measured RNA and the panel ------------------------------------------------------------------
#: GTEx v8 gene median TPM, the one per-tissue measured RNA table cached on this machine. ENCODE's
#: cell-line RNA-seq for the reader's cells (data/knowledge/rna/tracks_*.json) is held only as bigWig
#: URLs, not as values, so it is not used; nothing is downloaded.
GTEX = Path("data/cache/c5/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz")
#: the 54 tissue columns of that file, in its order, copied from its header line (the only line of it
#: read before the registration). The run refuses if the header differs.
GTEX_V8_TISSUES = (
    "Adipose - Subcutaneous",
    "Adipose - Visceral (Omentum)",
    "Adrenal Gland",
    "Artery - Aorta",
    "Artery - Coronary",
    "Artery - Tibial",
    "Bladder",
    "Brain - Amygdala",
    "Brain - Anterior cingulate cortex (BA24)",
    "Brain - Caudate (basal ganglia)",
    "Brain - Cerebellar Hemisphere",
    "Brain - Cerebellum",
    "Brain - Cortex",
    "Brain - Frontal Cortex (BA9)",
    "Brain - Hippocampus",
    "Brain - Hypothalamus",
    "Brain - Nucleus accumbens (basal ganglia)",
    "Brain - Putamen (basal ganglia)",
    "Brain - Spinal cord (cervical c-1)",
    "Brain - Substantia nigra",
    "Breast - Mammary Tissue",
    "Cells - Cultured fibroblasts",
    "Cells - EBV-transformed lymphocytes",
    "Cervix - Ectocervix",
    "Cervix - Endocervix",
    "Colon - Sigmoid",
    "Colon - Transverse",
    "Esophagus - Gastroesophageal Junction",
    "Esophagus - Mucosa",
    "Esophagus - Muscularis",
    "Fallopian Tube",
    "Heart - Atrial Appendage",
    "Heart - Left Ventricle",
    "Kidney - Cortex",
    "Kidney - Medulla",
    "Liver",
    "Lung",
    "Minor Salivary Gland",
    "Muscle - Skeletal",
    "Nerve - Tibial",
    "Ovary",
    "Pancreas",
    "Pituitary",
    "Prostate",
    "Skin - Not Sun Exposed (Suprapubic)",
    "Skin - Sun Exposed (Lower leg)",
    "Small Intestine - Terminal Ileum",
    "Spleen",
    "Stomach",
    "Testis",
    "Thyroid",
    "Uterus",
    "Vagina",
    "Whole Blood",
)


def tissue_label(name: str) -> str:
    """A GTEx column name as the compiled programs spell the AlphaGenome track named for it.

    Every run of characters other than letters and digits becomes one underscore, and the ends are
    stripped: "Brain - Spinal cord (cervical c-1)" -> "Brain_Spinal_cord_cervical_c_1". Matching is
    exact and case-sensitive, so the ENCODE biosample `liver` is NOT the GTEx track `Liver`.
    """
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")


#: compiled context label -> GTEx column. Only these 54 contexts are in the panel: they are the
#: compiled labels of AlphaGenome's GTEx-sourced RNA-seq tracks, so the mapping is by name and needs
#: no judgement. An ENCODE or ENTEx biosample (psoas_muscle, left_lobe_of_liver, ...) would need a
#: hand-made mapping and is left out rather than mapped after the fact.
PANEL = {tissue_label(t): t for t in GTEX_V8_TISSUES}

# ---- the census's structural facts (what decides what a run can predict) ------------------------------
STRUCTURE = {
    "network": (
        "a compiled rule's source is an element, never a gene: the census counts rules whose source is a"
        " gene entity, and if that count is 0 a run integrates one independent one-gene system per"
        " regulated gene, so nothing a run reports depends on any other gene"
    ),
    "level": (
        "the bridge fits each resolved gene so that its intact steady state is T / delta_m in every cell"
        " it resolves (split_measured_rate: intact transcription = T), so the predicted level of a gene"
        " does not depend on any compiled rule and is the same in every cell"
    ),
    "on_off": (
        "every resolved gene has basal_rate > 0, because the split reports one that does not as"
        " rate_split_infeasible; a run therefore predicts 'on' for every resolved gene in every cell,"
        " with its element present or absent, and cannot predict 'off' anywhere"
    ),
    "in_its_own_cell": (
        "with the element clamped present and absent in the rule's own cell, a run returns the"
        " observation it was fitted to (bridge.RATE_ACCEPTANCE; S5: 'direction and fold for all (they are"
        " the input)')"
    ),
}

# ---- the registered test ----------------------------------------------------------------------------
QUESTION = (
    "Does running the compiled rules predict measured RNA? Specifically: across the 54 GTEx tissues,"
    " is a gene's measured expression in the tissue its compiled rule is gated on better predicted by"
    " the level a run of that rule computes than by a baseline that ignores the compiled rules?"
)

READING = (
    "element-state reading of the compiled gate. A run of gene G's bridged program in tissue C with"
    " its element clamped PRESENT (1.0) gives L_present; the same program with the element clamped"
    " ABSENT (0.0) gives L_absent. A tissue C' in which G has no compiled rule at all is read as a"
    " tissue in which the element does not act, i.e. the ABSENT state - the bridge's own SOURCE_STATE"
    " ('element presence, dimensionless: 1.0 intact, 0.0 removed'). The run's prediction for G in C"
    " against G elsewhere is therefore m = log2(L_present / L_absent): positive for a rule whose action"
    " is 'activates', negative for 'inhibits'. m does not depend on the rate T, on the Hill threshold"
    " or coefficient, or on the mRNA half-life (the split is linear in T and the fold is"
    " dimensionless), so the borrowed median is used only to make a gene with no measured rate"
    " executable and moves no prediction"
)
READING_NOT_TAKEN = (
    "the runtime's literal reading of an UNGATED rule is not registered: a gene with no active"
    " activator takes A = 1.0 and transcribes at basal + max_rate (runtime/grn.py, `a = 1.0`), so"
    " with the C-fitted parameters a run in a tissue with no rule gives (2 - RHO) * T for an"
    " activating rule and RHO * T for an inhibiting one, and predicts G LOWER in its rule's tissue for"
    " both actions. grn.py's own comment calls that convention 'the opposite of silencing', and"
    " lane-mesoderm (d340e5c) recorded that no paper can source it. Its sign agreement is not a"
    " separate test: it is one minus the gate-only figure reported under DESCRIPTIVES"
)

POPULATION = (
    "every (gene G, tissue C) with C one of the 54 panel contexts such that, pooled over all 24"
    " compiled programs, G has EXACTLY ONE mechanism (bridge.mechanism_key) active in C, and"
    " bridge.parameterize resolves it on the borrowed tier (rates = data/results/gene_rates.json,"
    " borrowed_rate = its genome median), i.e. the pair is runnable at least in relative units; AND G's"
    " symbol names exactly one row of the GTEx table; AND G has at least MIN_REFERENCE_TISSUES panel"
    " tissues in which it has no compiled regulatory rule of any kind (its reference tissues R(G)). No"
    " filter on expression is applied: no measured value selects a pair"
)
MIN_REFERENCE_TISSUES = 20
MIN_POPULATION = 500
PSEUDOCOUNT = 1.0
MEASURED = (
    "y(G, T) = log2(median TPM of G in GTEx tissue T + PSEUDOCOUNT). The measured deviation of G in its"
    " rule's tissue is D(G, C) = y(G, C) - median over C' in R(G) of y(G, C')"
)
BASELINE = (
    "no-rule baseline: predicts D(G, C) = 0, i.e. G's level in C equals its median level over the"
    " tissues where it has no compiled rule. It uses no compiled rule, no gate and no sign, and it is"
    " the same level prediction as the run's except for the run's computed shift m"
)
METRIC = (
    "Delta MSE = mean over the population of D^2 - mean of (D - m)^2 (= mean of 2*m*D - m^2). Positive"
    " means the run's predicted level is closer to the measured level than the no-rule baseline's."
    " Stated before the run because it decides how a null reads: a shift that carries no information"
    " about D still costs mean(m^2), so a run whose shifts are unrelated to the measured deviation"
    " reads WORSE, not NO_DIFFERENCE, once the population resolves mean(m^2); NO_DIFFERENCE means the"
    " information the shifts carry and their cost could not be told apart"
)
BOOTSTRAP = {
    "resamples": 2000,
    "seed": 20261010,
    "unit": "gene: genes are resampled with replacement and carry all their (gene, tissue) pairs",
    "interval": "95% percentile interval of Delta MSE",
}
CONTROL = (
    "shuffled-rule control: within each tissue, the run's shifts m are permuted among that tissue's"
    " pairs (PERMUTATIONS times), which keeps every tissue's gate and its distribution of shifts and"
    " breaks only the tie between a gene and its own rule's sign and size. control p = (1 + number of"
    " permutations with Delta MSE >= observed) / (1 + PERMUTATIONS)"
)
PERMUTATIONS = {"count": 1000, "seed": 20261011}
CONTROL_ALPHA = 0.05
#: the runtime as run: RK4 from zero to steady state with the element held by clamp in every stage
RUNTIME = {"hours": 40.0, "dt": 0.05, "strict": True, "start": "every species at 0"}
#: the run must agree with the closed form m = -log2(RHO) at every pair before the RNA is opened
CLOSED_FORM_TOLERANCE = 1e-6

DECISION = (
    "1. If any pair's run differs from the closed form by more than CLOSED_FORM_TOLERANCE in log2, the"
    " run stops before the measured RNA is opened (HARNESS_DEFECT). 2. If the population holds fewer"
    " than MIN_POPULATION pairs, the run stops before any comparison is read (NOT_TESTABLE). 3. Else:"
    " lower bound of the bootstrap interval > 0 -> BEATS (then control p <= CONTROL_ALPHA ->"
    " BEATS_AND_CONTROL, else BEATS_NOT_CONTROL); upper bound < 0 -> WORSE; otherwise NO_DIFFERENCE"
)
FALSIFIER = (
    "the run does not beat the no-rule baseline: the 95% interval of Delta MSE includes 0 or lies below"
    " it. And the narrower claim that the gene's OWN rule carries information is falsified when the"
    " shuffled-rule control is not beaten (control p > CONTROL_ALPHA)"
)
STOP_RULE = (
    "one run, from this registration's committed code, on the cached files whose sha256 the"
    " registration pins. No parameter, panel, population rule, pseudocount or metric is changed after"
    " the measured RNA is opened; a defect found afterwards is reported beside the registered reading"
    " and corrected only under a new registration"
)
#: the exact sentence of each outcome; `reading` fills the numbers and nothing else
OUTCOMES = {
    "HARNESS_DEFECT": (
        "stopped before the measured RNA was opened: the run differed from its closed form by more"
        " than {tol} log2 at {bad} of {n} pairs"
    ),
    "NOT_TESTABLE": (
        "not testable at the registered size: the population holds {n} (gene, tissue) pairs, fewer"
        " than the registered {min}; no comparison was read"
    ),
    "BEATS_AND_CONTROL": (
        "beats baseline: running the compiled rules lowered the squared error of measured"
        " log2(TPM+1) below the no-rule baseline's (Delta MSE {d} [95% CI {lo}, {hi}] over {n}"
        " gene-tissue pairs, {g} genes, {t} GTEx tissues), and by more than the shuffled-rule control"
        " (control p = {p})"
    ),
    "BEATS_NOT_CONTROL": (
        "beats baseline, but not the shuffled-rule control: running the compiled rules lowered the"
        " squared error of measured log2(TPM+1) below the no-rule baseline's (Delta MSE {d} [95% CI"
        " {lo}, {hi}] over {n} gene-tissue pairs, {g} genes, {t} GTEx tissues), and shifts shuffled"
        " among the genes of the same tissue did as well (control p = {p}); the gain is not carried by"
        " a gene's own rule"
    ),
    "NO_DIFFERENCE": (
        "no difference detected: running the compiled rules did not change the squared error of"
        " measured log2(TPM+1) against the no-rule baseline by a detectable amount (Delta MSE {d}"
        " [95% CI {lo}, {hi}] over {n} gene-tissue pairs, {g} genes, {t} GTEx tissues; control p ="
        " {p})"
    ),
    "WORSE": (
        "worse: running the compiled rules raised the squared error of measured log2(TPM+1) above"
        " the no-rule baseline's (Delta MSE {d} [95% CI {lo}, {hi}] over {n} gene-tissue pairs, {g}"
        " genes, {t} GTEx tissues; control p = {p})"
    ),
}
#: carried with every outcome, word for word
LIMITATIONS = (
    "every rule in the population was compiled and gated from AlphaGenome RNA-seq tracks named for"
    " GTEx tissues; the saved track metadata records data_source 'gtex' for such a track (the"
    " Whole_Blood row was read), so the measured outcome is the kind of data the rules' source model"
    " was trained on and is not independent of it. No reading here is evidence of a mechanism",
    "the run's shift m equals the input observation's -log2(RHO): what is tested is the compiled"
    " observation read across tissues under READING, not dynamics that a run adds; there are none to"
    " add (STRUCTURE)",
    "a GTEx median over donors is not a cell type and not the AlphaGenome track's own sample set",
    "every rate under the run is mouse NIH 3T3 or the genome median of those rates; m does not depend"
    " on it, and no level in molecules per cell is quoted",
)
DESCRIPTIVES = (
    "sign table: rule action (activates, inhibits) against sign of D (+, 0, -), counts",
    "gate-only agreement: the share of pairs with D > 0, which is the agreement of a predictor that"
    " uses only the gate and always says 'higher in the rule's tissue'; the runtime-literal reading"
    " (READING_NOT_TAKEN) agrees on exactly the share with D < 0",
    "Delta MSE by action, by tissue (count, Delta MSE, mean D, mean m), and on the subset of genes"
    " standing on their own measured (mouse) rate",
    "Spearman correlation of m with D",
    "on/off, descriptive only and not a test of the run (which predicts 'on' everywhere): share of"
    " pairs with TPM >= 1 in the rule's tissue, and mean share of reference tissues with TPM >= 1",
)


# ==== census ===========================================================================================

REGULATORY = (Action.ACTIVATE, Action.INHIBIT)
TIERS = ("measured", "borrowed")


def _tier_kwargs(tier: str, rates: dict, borrowed: float, prov: dict, lent: dict) -> dict[str, Any]:
    if tier == "measured":
        return {"rates": rates, "rate_provenance": prov}
    return {"rates": rates, "borrowed_rate": borrowed, "rate_provenance": prov, "borrowed_provenance": lent}


def census_module(
    module: Module, rates: dict[str, float], borrowed: float, prov: dict, lent: dict
) -> dict[str, Any]:
    """Per cell context of one program: its shape, what the bridge resolves on each tier, and why not.

    Classify-only (`build=False`): no rule is written and nothing is run. Reads no measured RNA.
    """
    genes = {g.id for g in module.genes()}
    out: dict[str, Any] = {
        "rules": sum(1 for r in module.rules if r.action in REGULATORY),
        "rules_with_gene_source": sum(
            1 for r in module.rules if r.action in REGULATORY and r.source in genes
        ),
        "contexts": {},
    }
    contexts = sorted({r.when.get("cell_type") for r in module.rules if r.when.get("cell_type")})
    for cell in contexts:
        ctx = {"cell_type": cell}
        by_key: dict[tuple, list[Rule]] = {}
        for r in module.active_rules(ctx):
            if r.action in REGULATORY:
                by_key.setdefault(bridge.mechanism_key(r), []).append(r)
        per_gene = Counter(k[1] for k in by_key)
        c: dict[str, Any] = {
            "rules": sum(len(v) for v in by_key.values()),
            "mechanisms": len(by_key),
            "gene_cell_pairs": len(per_gene),
            "single_mechanism_pairs": sum(1 for n in per_gene.values() if n == 1),
            "several_mechanism_pairs": sum(1 for n in per_gene.values() if n > 1),
        }
        for tier in TIERS:
            b = bridge.parameterize(
                module, ctx, build=False, **_tier_kwargs(tier, rates, borrowed, prov, lent)
            )
            c[f"resolved_{tier}"] = len(b.strengths)
            c[f"resolved_{tier}_activates"] = sum(1 for o in b.observations.values() if o.rho < 1.0)
            c[f"resolved_{tier}_inhibits"] = sum(1 for o in b.observations.values() if o.rho > 1.0)
            c[f"resolved_{tier}_basal_not_positive"] = sum(
                1 for _t, bas, _v in b.rates_used.values() if bas <= 0
            )
            c[f"resolved_{tier}_on_measured_rate"] = sum(1 for t in b.tier.values() if t == "measured")
            c[f"unresolved_{tier}"] = dict(Counter(i.reason for i in b.unresolved))
        out["contexts"][cell] = c
    return out


def pool_census(per_program: dict[str, dict]) -> dict[str, Any]:
    """Sum the per-program census by context; numbers only, reasons merged."""
    by_ctx: dict[str, Counter] = defaultdict(Counter)
    reasons: dict[str, dict[str, Counter]] = defaultdict(lambda: {t: Counter() for t in TIERS})
    totals = Counter()
    for prog in per_program.values():
        totals["rules"] += prog["rules"]
        totals["rules_with_gene_source"] += prog["rules_with_gene_source"]
        for cell, c in prog["contexts"].items():
            for k, v in c.items():
                if k.startswith("unresolved_"):
                    reasons[cell][k[len("unresolved_") :]].update(v)
                else:
                    by_ctx[cell][k] += v
    contexts = {
        cell: {**dict(by_ctx[cell]), **{f"unresolved_{t}": dict(reasons[cell][t]) for t in TIERS}}
        for cell in sorted(by_ctx, key=lambda x: -by_ctx[x]["rules"])
    }
    return {"totals": dict(totals), "contexts": contexts}


def summarise(contexts: dict[str, dict], cells: set[str] | None = None) -> dict[str, Any]:
    """Pooled counts over a set of contexts (all when None)."""
    s = Counter()
    why = {t: Counter() for t in TIERS}
    for cell, c in contexts.items():
        if cells is not None and cell not in cells:
            continue
        s["contexts"] += 1
        for k, v in c.items():
            if k.startswith("unresolved_"):
                why[k[len("unresolved_") :]].update(v)
            else:
                s[k] += v
    return {**dict(s), **{f"unresolved_{t}": dict(why[t]) for t in TIERS}}


# ==== the run ==========================================================================================


def rule_tissues(module: Module, panel: dict[str, str] = PANEL) -> dict[str, set[str]]:
    """gene -> the panel contexts in which it has ANY compiled regulatory rule (its non-reference tissues)."""
    out: dict[str, set[str]] = defaultdict(set)
    for r in module.rules:
        cell = r.when.get("cell_type")
        if r.action in REGULATORY and cell in panel:
            out[r.target].add(cell)
    return out


def mechanisms_in(module: Module, ctx: dict[str, str]) -> Counter:
    """gene -> number of mechanisms active in this context (pooled later over programs)."""
    keys = {bridge.mechanism_key(r) for r in module.active_rules(ctx) if r.action in REGULATORY}
    return Counter(k[1] for k in keys)


def one_gene_module(module: Module, gene: str) -> Module:
    """The gene, every regulatory rule that targets it (all contexts) and those rules' source entities.

    `bridge.parameterize` groups by gene and reads only the gene's own rules and their sources' evidence,
    so on this module it gives the gene exactly what it gives it on the whole program; the run asserts
    the fitted strength is the same on both.
    """
    rules = [copy.deepcopy(r) for r in module.rules if r.action in REGULATORY and r.target == gene]
    ents = {gene: copy.deepcopy(module.entities[gene])}
    for r in rules:
        if r.source in module.entities:
            ents[r.source] = copy.deepcopy(module.entities[r.source])
    return Module(name=f"{module.name}.{gene}", entities=ents, rules=rules)


def run_pair(
    small: Module,
    gene: str,
    ctx: dict[str, str],
    rates: dict[str, float],
    borrowed: float,
    prov: dict | None = None,
    lent: dict | None = None,
    runtime: dict[str, Any] = RUNTIME,
) -> dict[str, Any] | None:
    """Bridge one gene in one context on the borrowed tier and run it with its element present and absent.

    Returns None when the bridge does not resolve the gene there. Otherwise the run's own numbers:
    L_present, L_absent, m = log2(L_present / L_absent), the observation RHO and the closed form
    -log2(RHO) it must equal.
    """
    b = bridge.parameterize(
        small, ctx, build=True, rates=rates, borrowed_rate=borrowed, rate_provenance=prov or {},
        borrowed_provenance=lent or {},
    )  # fmt: skip
    keys = [k for k in b.strengths if k[1] == gene]
    if not keys:
        return None
    (key,) = keys
    element = key[0]
    species = f"{gene}.mRNA"
    levels = {}
    for state in (bridge.INTACT, bridge.REMOVED):
        rt = NetworkRuntime(b.module, context=ctx, strict=runtime["strict"])
        steps = int(round(runtime["hours"] / runtime["dt"]))
        traj = rt.run(hours=runtime["hours"], dt=runtime["dt"], clamp={element: state}, record_every=steps)
        levels[state] = traj.final()[species]
    rho = b.observations[key].rho
    present, absent = levels[bridge.INTACT], levels[bridge.REMOVED]
    m = math.log2(present / absent)
    return {
        "gene": gene,
        "tissue": ctx["cell_type"],
        "element": element,
        "action": "activates" if rho < 1.0 else "inhibits",
        "rho": rho,
        "strength": b.strengths[key],
        "tier": b.tier.get(gene, ""),
        "L_present": present,
        "L_absent": absent,
        "m": m,
        "closed_form": -math.log2(rho),
    }


# ==== the measured RNA and the score ===================================================================


def load_gtex(path: Path = GTEX) -> tuple[list[str], dict[str, list[float]], dict[str, int]]:
    """(tissue columns, symbol -> TPM row for symbols naming exactly one row, counts).

    Refuses a header that is not GTEX_V8_TISSUES. A symbol on several rows is ambiguous and dropped.
    """
    rows: dict[str, list[list[float]]] = defaultdict(list)
    with gzip.open(path, "rt") as fh:
        fh.readline()
        fh.readline()
        header = fh.readline().rstrip("\n").split("\t")
        tissues = header[2:]
        if tuple(tissues) != GTEX_V8_TISSUES:
            raise ValueError(f"{path}: the tissue header is not the registered one")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows[f[1]].append([float(x) for x in f[2:]])
    unique = {s: r[0] for s, r in rows.items() if len(r) == 1}
    counts = {
        "rows": sum(len(r) for r in rows.values()),
        "symbols": len(rows),
        "symbols_on_one_row": len(unique),
        "symbols_on_several_rows": sum(1 for r in rows.values() if len(r) > 1),
    }
    return list(tissues), unique, counts


def _y(tpm: float) -> float:
    return math.log2(tpm + PSEUDOCOUNT)


def _median(xs: list[float]) -> float:
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def assemble(
    pairs: list[dict], gene_rule_tissues: dict[str, set[str]], tpm: dict[str, list[float]], tissues: list[str]
) -> tuple[list[dict], Counter]:
    """The population with each pair's D, and the count of every exclusion by reason."""
    col = {tissue_label(t): i for i, t in enumerate(tissues)}
    panel = list(col)
    out, why = [], Counter()
    for p in pairs:
        row = tpm.get(p["gene"])
        if row is None:
            why["symbol_not_on_exactly_one_gtex_row"] += 1
            continue
        ref = [c for c in panel if c not in gene_rule_tissues.get(p["gene"], set())]
        if len(ref) < MIN_REFERENCE_TISSUES:
            why["fewer_than_min_reference_tissues"] += 1
            continue
        y_c = _y(row[col[p["tissue"]]])
        d = y_c - _median([_y(row[col[c]]) for c in ref])
        on_ref = sum(1 for c in ref if row[col[c]] >= 1.0) / len(ref)
        out.append(
            {**p, "D": d, "on_in_rule_tissue": row[col[p["tissue"]]] >= 1.0, "on_share_reference": on_ref}
        )
    return out, why


def delta_mse(m: list[float], d: list[float]) -> float:
    n = len(m)
    return sum(2.0 * a * b - a * a for a, b in zip(m, d, strict=True)) / n if n else float("nan")


def cluster_bootstrap(pop: list[dict], resamples: int, seed: int) -> tuple[float, float]:
    """95% percentile interval of Delta MSE, resampling genes with all their pairs."""
    import numpy as np

    by_gene: dict[str, list[int]] = defaultdict(list)
    for i, p in enumerate(pop):
        by_gene[p["gene"]].append(i)
    genes = sorted(by_gene)
    s = np.array([sum(2 * pop[i]["m"] * pop[i]["D"] - pop[i]["m"] ** 2 for i in by_gene[g]) for g in genes])
    n = np.array([len(by_gene[g]) for g in genes], dtype=float)
    rng = np.random.default_rng(seed)
    stats = np.empty(resamples)
    for k in range(resamples):
        idx = rng.integers(0, len(genes), len(genes))
        stats[k] = s[idx].sum() / n[idx].sum()
    lo, hi = np.percentile(stats, [2.5, 97.5])
    return float(lo), float(hi)


def within_tissue_permutation(pop: list[dict], count: int, seed: int) -> list[float]:
    """Delta MSE with each tissue's shifts permuted among its own pairs, `count` times."""
    import numpy as np

    rng = np.random.default_rng(seed)
    groups: dict[str, list[int]] = defaultdict(list)
    for i, p in enumerate(pop):
        groups[p["tissue"]].append(i)
    m = np.array([p["m"] for p in pop])
    d = np.array([p["D"] for p in pop])
    m2 = float((m * m).sum())
    out = []
    for _ in range(count):
        mp = m.copy()
        for idx in groups.values():
            ix = np.array(idx)
            mp[ix] = m[rng.permutation(ix)]
        out.append(float((2.0 * (mp * d).sum() - m2) / len(pop)))
    return out


def reading(kind: str, **numbers: Any) -> str:
    """The registered sentence for an outcome, numbers filled and nothing else changed."""
    return OUTCOMES[kind].format(**numbers)


def decide(d: float, lo: float, hi: float, p: float) -> str:
    """DECISION, step 3."""
    if lo > 0:
        return "BEATS_AND_CONTROL" if p <= CONTROL_ALPHA else "BEATS_NOT_CONTROL"
    if hi < 0:
        return "WORSE"
    return "NO_DIFFERENCE"


def _spearman(a: list[float], b: list[float]) -> float | None:
    from scipy.stats import spearmanr

    if len(a) < 3:
        return None
    r = spearmanr(a, b).statistic
    return None if r != r else float(r)


def score(pop: list[dict]) -> dict[str, Any]:
    """The registered comparison and the registered descriptives over an assembled population."""
    n = len(pop)
    genes = len({p["gene"] for p in pop})
    tissues = len({p["tissue"] for p in pop})
    if n < MIN_POPULATION:
        return {
            "outcome": "NOT_TESTABLE",
            "reading": reading("NOT_TESTABLE", n=n, min=MIN_POPULATION),
            "pairs": n,
        }
    m = [p["m"] for p in pop]
    d = [p["D"] for p in pop]
    dm = delta_mse(m, d)
    lo, hi = cluster_bootstrap(pop, BOOTSTRAP["resamples"], BOOTSTRAP["seed"])
    perm = within_tissue_permutation(pop, PERMUTATIONS["count"], PERMUTATIONS["seed"])
    p_ctrl = (1 + sum(1 for x in perm if x >= dm)) / (1 + len(perm))
    kind = decide(dm, lo, hi, p_ctrl)
    nums = {
        "d": f"{dm:+.4g}",
        "lo": f"{lo:+.4g}",
        "hi": f"{hi:+.4g}",
        "n": f"{n:,}",
        "g": f"{genes:,}",
        "t": tissues,
        "p": f"{p_ctrl:.3g}",
    }
    signs = Counter((p["action"], "+" if p["D"] > 0 else "-" if p["D"] < 0 else "0") for p in pop)

    def sub(rows: list[dict]) -> dict[str, Any]:
        if not rows:
            return {"pairs": 0}
        return {
            "pairs": len(rows),
            "genes": len({r["gene"] for r in rows}),
            "delta_mse": delta_mse([r["m"] for r in rows], [r["D"] for r in rows]),
            "mean_D": sum(r["D"] for r in rows) / len(rows),
            "mean_m": sum(r["m"] for r in rows) / len(rows),
        }

    by_tissue = defaultdict(list)
    for p in pop:
        by_tissue[p["tissue"]].append(p)
    measured = [p for p in pop if p["tier"] == "measured"]
    measured_ci = {}
    if len(measured) >= MIN_POPULATION:
        lo_m, hi_m = cluster_bootstrap(measured, BOOTSTRAP["resamples"], BOOTSTRAP["seed"])
        measured_ci = {"ci95": [lo_m, hi_m]}
    return {
        "outcome": kind,
        "reading": reading(kind, **nums),
        "limitations": list(LIMITATIONS),
        "pairs": n,
        "genes": genes,
        "tissues": tissues,
        "delta_mse": dm,
        "ci95": [lo, hi],
        "mse_baseline": sum(x * x for x in d) / n,
        "mse_run": sum((x - y) ** 2 for x, y in zip(d, m, strict=True)) / n,
        "control_p": p_ctrl,
        "control_delta_mse_median": _median(perm),
        "control_delta_mse_95th_percentile": sorted(perm)[int(0.95 * len(perm)) - 1],
        "descriptives": {
            "sign_table": {f"{a} D{s}": c for (a, s), c in sorted(signs.items())},
            "gate_only_agreement_D_positive": sum(1 for x in d if x > 0) / n,
            "runtime_literal_agreement_D_negative": sum(1 for x in d if x < 0) / n,
            "run_sign_agreement": sum(1 for p in pop if (p["m"] > 0 < p["D"]) or (p["m"] < 0 > p["D"])) / n,
            "by_action": {a: sub([p for p in pop if p["action"] == a]) for a in ("activates", "inhibits")},
            "measured_rate_subset": {**sub(measured), **measured_ci},
            "spearman_m_D": _spearman(m, d),
            "on_off": {
                "share_on_in_rule_tissue": sum(1 for p in pop if p["on_in_rule_tissue"]) / n,
                "mean_share_on_in_reference_tissues": sum(p["on_share_reference"] for p in pop) / n,
                "note": "descriptive only: the run predicts 'on' in every tissue (STRUCTURE['on_off'])",
            },
            "by_tissue": {t: sub(rows) for t, rows in sorted(by_tissue.items(), key=lambda kv: -len(kv[1]))},
        },
    }
