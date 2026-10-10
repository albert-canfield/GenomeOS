# SPDX-License-Identifier: AGPL-3.0-or-later
"""M2b (ROADMAP §5 item 14, lane-readercompile, 2026-10-10): does chromatin gating make a run track RNA?

M1 (`attribution/runrules.py`, 747b207) found that a run of the compiled rules makes one cross-tissue
prediction and that it restates its own input: the rule's cell is the argmax AlphaGenome track, so
"higher where the rule is" is the model's observation read back. Registered reading: no difference
detected against a no-rule baseline on GTEx.

The reader construct (BioLang v0.5, 5682e6c) and its application half (`compile.gate_program`, this
lane) change what a run can say. Every predicted rule acts `through` its element, the rule is offered
in each of the thirteen cell types reader v1 read (`SELECTION_READER_CELLS`), and `Module.active_rules`
runs it only where that cell's reader states the element open, holding it at strength 0 elsewhere. A
gene's element-driven contribution is then on in some reader cells and held in others, and the
pattern comes from measured chromatin (ENCODE DNase-seq), not from the model's track argmax. That is
the first cross-cell prediction a compiled run makes that does not restate its input, and this module
registers the test of it against measured RNA in the same thirteen biosamples.

This module holds, in the order they were written and committed:

1. the code that runs the gated programs (no RNA opened) and the code that scores a run;
2. the registered test (`QUESTION` .. `OUTCOMES`), committed before any RNA value was read.

No model request, no key, no money. The RNA quantifications are public ENCODE files fetched once by
scripts/readergate_fetch.py and pinned by checksum before the registration.
"""

from __future__ import annotations

import copy
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from genomeos.attribution import bridge
from genomeos.attribution import compile as cp
from genomeos.attribution import runrules as rr
from genomeos.ir import Action, Module, Reader
from genomeos.runtime.grn import NetworkRuntime

REGISTERED = "2026-10-10"
LANE = "lane-readercompile"
REGULATORY = (Action.ACTIVATE, Action.INHIBIT)
SELECTION = cp.SELECTION_READER_CELLS

#: the thirteen reader biosamples as the gated programs declare them (compile.reader_cells)
CELLS = tuple(cid for _, cid in cp.reader_cells())

# ---- the measured RNA ---------------------------------------------------------------------------------
RNA_DIR = Path("data/cache/readergate/rna")
RNA_SOURCES = Path("data/results/readergate_rna_sources.json")
GENE_TABLE = Path("data/cache/gencode_genes.tsv")
RNA_ASSAY = "polyA plus RNA-seq"
RNA_EXCLUDE_WORDS = ("fraction", "treated", "with ", "genetically", "transfected", "knockdown")
RNA_SELECTION = (
    "per reader biosample, from the ENCODE portal: every released GRCh38 'gene quantifications' TSV with"
    " genome_annotation V29 (the ENCODE4 uniform RNA-seq pipeline, RSEM) of assay 'polyA plus RNA-seq',"
    " on the biosample's ontology term name, from a human experiment the portal lists as not perturbed"
    " (perturbed=false) and that carries no ERROR audit, and whose"
    " file biosample summary names no subcellular fraction, treatment, disease or genetic modification"
    " (none of 'fraction', 'treated', 'with ', 'genetically', 'transfected', 'knockdown'). polyA plus is"
    " the one assay with such an experiment for all thirteen (SK-N-SH has no whole-cell untreated total"
    " RNA-seq in this query). Of the qualifying experiments the one with the most qualifying files is"
    " taken, ties broken by the latest date_released and then the lowest accession; the gene's TPM in"
    " the biosample is the mean of its TPM over that experiment's qualifying files (its replicates)."
    " Chosen from portal metadata alone; no value of any file was read to choose"
)

# ---- what the gated programs can predict (written from the code, before any RNA was read) -----------
STRUCTURE = {
    "identifiable": (
        "the bridge fits one strength per mechanism and refuses a gene regulated by more than one in a"
        " cell (not_identifiable). Under SELECTION_READER_CELLS every predicted rule of a gene is offered"
        " in every reader cell, so a gene is runnable only if exactly ONE regulatory rule targets it over"
        " all 24 programs; a gene with several elements has no fitted parameters to run"
    ),
    "two_levels": (
        "for such a gene the gated run has two levels: L_open, where the cell's reader states the element"
        " open and the rule runs, and L_held, where it is closed or unknown and the rule is held at"
        " strength 0 (for one activating rule the basal rate; for one inhibiting rule basal + max_rate)."
        " m = log2(L_open / L_held) equals the observation's -log2(RHO) and does not depend on the"
        " borrowed rate, the Hill parameters or the mRNA half-life (runrules.READING)"
    ),
    "offered_everywhere": (
        "a run that offers the rule in every reader cell WITHOUT a reader runs it everywhere: one level in"
        " all thirteen cells, which after centring within the gene predicts no difference between cells."
        " Under the registered endpoint it is the no-rule baseline exactly, so it is not a separate test"
    ),
    "network": "no compiled rule's source is a gene (runrules census): one independent system per gene",
    "bridge_drops_through": (
        "bridge.parameterize rebuilds each fitted rule without `through`, so a bridged program runs"
        " ungated unless `through` is put back; run_gene puts it back and checks the gated level against"
        " the clamped level of the cell's state at every gene and cell"
    ),
}

QUESTION = (
    "With chromatin gating, a compiled run predicts a gene's element-driven contribution ON only in the"
    " reader cell types whose reader states the element open. Across the thirteen reader cell types,"
    " does that gated prediction track measured RNA better than a no-rule baseline, and better than the"
    " ungated run as compiled (the rule on in its compiled cell only)?"
)
POPULATION = (
    "every gene G such that, over all 24 compiled programs (data/knowledge/compiled, the files M1 ran,"
    " pinned by sha256) gated with compile.gate_program(selection='reader_cells'): (a) exactly one"
    " activates/inhibits rule line targets G, and its evidence is predicted (one element E, one compiled"
    " cell C); (b) bridge.parameterize resolves G on the borrowed tier (data/results/gene_rates.json"
    " and its genome median) in each of the thirteen reader cells; (c) E's reader state is open in at"
    " least one reader cell and not open (closed or unknown) in at least one, so the gated run predicts"
    " a difference between cells (a predictor-side condition; no measured value enters it); (d) G's"
    " symbol names exactly one gene_id in data/cache/gencode_genes.tsv, and that gene_id without its"
    " version is on exactly one row of every one of the thirteen quantification tables (rows suffixed"
    " _PAR_Y are not counted). Each G contributes its thirteen (G, cell) pairs. No filter on expression"
)
MIN_GENES = 100
PSEUDOCOUNT = 1.0
MEASURED = (
    "y(G, X) = log2(TPM(G, X) + PSEUDOCOUNT), TPM(G, X) the mean over the chosen experiment's files"
    " (RNA_SELECTION). The measured deviation is D(G, X) = y(G, X) - mean over the thirteen cells X' of"
    " y(G, X'): G's level in X relative to its own mean across the reader cells"
)
PREDICTIONS = {
    "gated": (
        "s(G, X) = log2(L(G, X) / L_held(G)), L(G, X) the level the gated run of G gives in cell X (the"
        " runtime with the bridged rule acting through E and the reader blocks in force): m where E is"
        " open in X, 0 where the rule is held. P(G, X) = s(G, X) - mean over X' of s(G, X')"
    ),
    "ungated_as_compiled": (
        "M1's READING applied in the reader cells: s(G, X) = m where X is G's compiled cell C, 0 in every"
        " other reader cell (a cell where G has no compiled rule read as the element-absent state); 0 in"
        " all thirteen when C is not a reader cell. Centred the same way"
    ),
    "no_rule_baseline": "P(G, X) = 0: G's level in X equals its mean over the thirteen reader cells",
}
METRIC = (
    "Delta MSE(pred) = mean over the population's pairs of D^2 - mean of (D - P_pred)^2, positive when"
    " the prediction's squared error is below the baseline's. A shift unrelated to D still costs"
    " mean(P^2), so an uninformative prediction reads WORSE once the population resolves that cost."
    " The ungated comparison is Delta MSE(gated) - Delta MSE(ungated_as_compiled), i.e. MSE(ungated) -"
    " MSE(gated), positive when the gated run is closer"
)
BOOTSTRAP = {
    "resamples": 2000,
    "seed": 20261010,
    "unit": "gene: genes are resampled with replacement and carry all thirteen of their pairs",
    "interval": "95% percentile interval",
}
CONTROL = (
    "shuffled-reader control: the thirteen-cell open/held rows of E's reader states are permuted among"
    " the population's genes (PERMUTATIONS times); every gene keeps its own m and compiled cell, so the"
    " control keeps the reader's per-cell open counts and the distribution of rows and breaks only the"
    " tie between a gene and its own element's chromatin. control p = (1 + permutations with Delta"
    " MSE(gated) >= observed) / (1 + PERMUTATIONS)"
)
PERMUTATIONS = {"count": 1000, "seed": 20261011}
CONTROL_ALPHA = 0.05
RUNTIME = rr.RUNTIME
CLOSED_FORM_TOLERANCE = rr.CLOSED_FORM_TOLERANCE
#: the gated run's level in a cell must equal the same bridged program run with the element clamped to
#: its state's value (open: present; held: absent), relative to the larger of the two levels
GATE_TOLERANCE = 1e-9

DECISION = (
    "1. If any gene's m differs from -log2(RHO) by more than CLOSED_FORM_TOLERANCE in log2, or the gated"
    " run's level in any cell differs from the clamped level of its state (open: element present; held:"
    " element absent) by more than GATE_TOLERANCE relative, the run stops before the RNA is opened"
    " (HARNESS_DEFECT). 2. If the population holds fewer than MIN_GENES genes, the run stops before any"
    " comparison is read (NOT_TESTABLE). 3. Else, on Delta MSE(gated) and its bootstrap interval: lower"
    " bound > 0 -> BEATS (control p <= CONTROL_ALPHA -> BEATS_AND_CONTROL, else BEATS_NOT_CONTROL);"
    " upper bound < 0 -> WORSE; otherwise NO_DIFFERENCE. 4. The ungated clause, from the interval of"
    " MSE(ungated) - MSE(gated): lower bound > 0 -> 'closer than', upper bound < 0 -> 'further than',"
    " otherwise 'not detectably different from'"
)
FALSIFIER = (
    "the gated run does not beat the no-rule baseline (the 95% interval of Delta MSE(gated) includes 0"
    " or lies below it): chromatin gating of the compiled rules does not make a run track measured RNA"
    " across the reader's cell types in this population. The narrower claim that the gene's OWN element's"
    " chromatin carries what it tracks is falsified when the shuffled-reader control is not beaten"
    " (control p > CONTROL_ALPHA); the claim over the ungated run is falsified when the ungated clause"
    " is not 'closer than'"
)
STOP_RULE = (
    "one run, from this registration's committed code, on the files the registration pins by sha256 (the"
    " 24 compiled programs, the rate table, the gene table, the 312 narrowPeak files and the RNA"
    " quantifications). No population rule, endpoint, pseudocount, metric, control or wording is changed"
    " after an RNA value is read; a defect found afterwards is reported beside the registered reading"
    " and corrected only under a new registration"
)
UNGATED_CLAUSES = {
    "closer": "closer to measured RNA than the ungated run as compiled",
    "further": "further from measured RNA than the ungated run as compiled",
    "same": "not detectably different from the ungated run as compiled",
}
#: the exact sentence of each outcome; `reading` fills the numbers and the ungated clause, nothing else
OUTCOMES = {
    "HARNESS_DEFECT": (
        "stopped before the measured RNA was opened: the gated run differed from its closed form or from"
        " its clamped state at {bad} of {n} genes"
    ),
    "NOT_TESTABLE": (
        "not testable at the registered size: the population holds {n} genes, fewer than the registered"
        " {min}; no comparison was read"
    ),
    "BEATS_AND_CONTROL": (
        "beats baseline: the reader-gated run lowered the squared error of measured log2(TPM+1) across the"
        " thirteen reader cell types below the no-rule baseline's (Delta MSE {d} [95% CI {lo}, {hi}] over"
        " {n} gene-cell pairs, {g} genes), by more than the shuffled-reader control (control p = {p}), and"
        " was {versus} (MSE difference {u} [95% CI {ulo}, {uhi}])"
    ),
    "BEATS_NOT_CONTROL": (
        "beats baseline, but not the shuffled-reader control: the reader-gated run lowered the squared"
        " error of measured log2(TPM+1) across the thirteen reader cell types below the no-rule baseline's"
        " (Delta MSE {d} [95% CI {lo}, {hi}] over {n} gene-cell pairs, {g} genes), and other genes'"
        " elements' chromatin did as well (control p = {p}); the gain is not carried by the gene's own"
        " element. It was {versus} (MSE difference {u} [95% CI {ulo}, {uhi}])"
    ),
    "NO_DIFFERENCE": (
        "no difference detected: the reader-gated run did not change the squared error of measured"
        " log2(TPM+1) across the thirteen reader cell types against the no-rule baseline by a detectable"
        " amount (Delta MSE {d} [95% CI {lo}, {hi}] over {n} gene-cell pairs, {g} genes; control p = {p})"
        " and was {versus} (MSE difference {u} [95% CI {ulo}, {uhi}])"
    ),
    "WORSE": (
        "worse: the reader-gated run raised the squared error of measured log2(TPM+1) across the thirteen"
        " reader cell types above the no-rule baseline's (Delta MSE {d} [95% CI {lo}, {hi}] over {n}"
        " gene-cell pairs, {g} genes; control p = {p}) and was {versus} (MSE difference {u} [95% CI {ulo},"
        " {uhi}])"
    ),
}
#: carried with every outcome, word for word. The first is the independence caveat, stated up front.
LIMITATIONS = (
    "independence: each rule's sign and size m come from AlphaGenome's RNA-seq gene scorer (an in-silico"
    " deletion, the argmax over 371 RNA-seq tracks), a model trained on ENCODE and GTEx RNA-seq, with"
    " RNA-seq tracks named for the reader biosamples among them, and on ENCODE DNase-seq; the reader's"
    " open/closed calls are ENCODE DNase-seq narrowPeak files, a measurement of the same biosamples; the"
    " measured outcome is ENCODE polyA plus RNA-seq of those biosamples, data of the kind, and possibly"
    " the very experiments, the model was trained on. The cross-cell pattern of the gated prediction"
    " comes from the DNase calls alone and its sign and size from the model",
    "accessibility near a gene tracks its expression across cell types whatever any compiled rule says,"
    " so a gated run that beats the baseline may carry only that association; the shuffled-reader"
    " control and the unsigned-openness descriptive bound how much of a gain belongs to the gene's own"
    " element and to the rule's sign. An open element is one a rule may act through, not evidence that it"
    " does; no reading here is evidence of a mechanism",
    "the population is the genes with exactly one compiled rule over the genome, a small and selected"
    " subset of the genes with any rule: a gene with several elements is not identifiable under the"
    " bridge and is not run",
    "m transports the AlphaGenome effect of the rule's compiled cell (the track where the deletion moved"
    " the gene most) to every reader cell where the element is open",
    "one polyA plus experiment per biosample, from several labs; two of the thirteen are tissues (ovary,"
    " testis) and the rest cell lines or primary or derived cells. Every rate under the run is mouse NIH"
    " 3T3 or their genome median; m does not depend on it and no level in molecules per cell is quoted",
)
DESCRIPTIVES = (
    "Delta MSE by rule action (activates, inhibits)",
    "unsigned openness: P = |m| where E is open, centred, which ignores the rule's sign (the"
    " accessibility-only association); Delta MSE and its interval",
    "Delta MSE per reader cell; Spearman correlation of P(gated) with D over pairs with P != 0",
    "the gated run as compiled (selection 'compiled_cell': on in C only where E is open in C), Delta MSE",
    "reader state counts over the population per cell (open, closed, unknown); genes whose compiled"
    " cell is a reader cell",
)


# ==== running the gated programs (no RNA) ================================================================

_COMPILED_WHEN = re.compile(
    r"^rule (\S+) (?:activates|inhibits) (\S+) \{ strength: [^;]*; when: cell_type = ([^;]*); "
)


def compiled_cells(text: str) -> dict[tuple[str, str], str]:
    """(source element, target) -> the compiled `when: cell_type` label, from an UNGATED program's text."""
    out = {}
    for ln in text.splitlines():
        m = _COMPILED_WHEN.match(ln)
        if m:
            out[(m.group(1), m.group(2))] = m.group(3)
    return out


_RULE_KIND = re.compile(
    r'^rule (\S+) (?:activates|inhibits) (\S+) \{ .*?evidence: (predicted|experimental) "'
)


def rule_targets(text: str) -> list[tuple[str, str, str]]:
    """(source element, target gene id, evidence kind) for every activates/inhibits rule line of a program."""
    out = []
    for ln in text.splitlines():
        if ln.startswith("rule "):
            m = _RULE_KIND.match(ln)
            if m is None:
                raise ValueError(f"a rule line this census cannot read: {ln[:120]}")
            out.append((m.group(1), m.group(2), m.group(3)))
    return out


def single_rule_genes(texts: dict[str, str]) -> tuple[dict[str, tuple[str, str]], Counter]:
    """gene id -> (chromosome, its one rule's element), over every program, for genes targeted by exactly
    one activates/inhibits rule line and that rule predicted (POPULATION (a)); exclusions counted."""
    seen: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    for chrom, text in texts.items():
        for src, gene, kind in rule_targets(text):
            seen[gene].append((chrom, src, kind))
    out: dict[str, tuple[str, str]] = {}
    why: Counter = Counter()
    for g, rows in seen.items():
        if len(rows) != 1:
            why["several_rules"] += 1
            continue
        chrom, src, kind = rows[0]
        if kind != "predicted":
            why["one_rule_not_predicted"] += 1
            continue
        out[g] = (chrom, src)
    why["genes_with_any_rule"] = len(seen)
    why["single_predicted_rule"] = len(out)
    return out, why


def run_chromosome(
    chrom: str,
    text: str,
    singles: dict[str, tuple[str, str]],
    rates: dict[str, float],
    borrowed: float,
    prov: dict | None = None,
    lent: dict | None = None,
    results_dir: Path | None = None,
) -> tuple[dict[str, Any], list[dict], Counter, Module]:
    """Gate one ungated program (compile.gate_program, SELECTION), parse it, and run every single-rule
    gene it holds (`run_gene`). Returns the gating counts, the gene rows, the counts on the way and
    the parsed gated program."""
    from genomeos.lang.parser import parse
    from genomeos.results import RESULTS_DIR

    cells_of = compiled_cells(text)
    gated, info = cp.gate_program(text, chrom, results_dir or RESULTS_DIR, selection=SELECTION)
    module = parse(gated, f"noncoding_{chrom}_reader")
    counts: Counter = Counter()
    rows = []
    for gene, (c, element) in sorted(singles.items()):
        if c != chrom:
            continue
        rules = [r for r in module.rules if r.action in REGULATORY and r.target == gene]
        if len(rules) != 1 or rules[0].through != element:
            raise ValueError(
                f"{chrom}: {gene} has {len(rules)} parsed rules where the text has one on {element}"
            )
        counts["genes_run"] += 1
        out = run_gene(gated_gene_module(module, gene), gene, rates, borrowed, prov, lent)
        if out is None:
            counts["not_resolved_in_every_reader_cell"] += 1
            continue
        out["chrom"] = chrom
        out["compiled_cell"] = cells_of[(element, gene)]
        rows.append(out)
    counts["genes_resolved"] = len(rows)
    info["gated_text_bytes"] = len(gated.encode())
    return info, rows, counts, module


def gated_gene_module(module: Module, gene: str) -> Module:
    """The gene, its rules and their sources (runrules.one_gene_module) plus the reader cell types and a
    reader per cell stating only the elements those rules act through. On this module the gene's
    selection and gate in any cell are what they are on the whole gated program."""
    small = rr.one_gene_module(module, gene)
    through = {r.through for r in small.rules if r.through}
    for c in module.cell_types():
        small.entities[c.id] = copy.deepcopy(c)
    for rd in module.readers:
        st = {e: s for e, s in rd.states.items() if e in through}
        if st:
            small.readers.append(Reader(cell_type=rd.cell_type, states=st, evidence=rd.evidence))
    return small


def _level(module: Module, gene: str, ctx: dict[str, str], clamp: dict[str, float]) -> float:
    rt = NetworkRuntime(module, context=ctx, strict=RUNTIME["strict"])
    steps = int(round(RUNTIME["hours"] / RUNTIME["dt"]))
    traj = rt.run(hours=RUNTIME["hours"], dt=RUNTIME["dt"], clamp=clamp, record_every=steps)
    return traj.final()[f"{gene}.mRNA"]


def run_gene(
    small: Module,
    gene: str,
    rates: dict[str, float],
    borrowed: float,
    prov: dict | None = None,
    lent: dict | None = None,
    cells: tuple[str, ...] = CELLS,
) -> dict[str, Any] | None:
    """The gated run of one single-rule gene in every reader cell, and the two clamped levels it must equal.

    In each cell: bridge the gene (borrowed tier), put `through` back on the fitted rule (the bridge
    rebuilds rules without it, STRUCTURE['bridge_drops_through']), and run the runtime with the element
    held present, so `Module.active_rules` applies the reader. Once, the same fitted program with its
    rule ungated is run with the element clamped present and absent: L_open and L_held. None when the
    bridge does not resolve the gene in every cell, or fits it differently in two cells.
    """
    (rule,) = [r for r in small.rules if r.action in REGULATORY and r.target == gene]
    element = rule.through
    if not element:
        raise ValueError(f"{gene}: its rule acts through no element; gate the program first")
    out: dict[str, Any] = {
        "gene": gene,
        "symbol": str(getattr(small.entities.get(gene), "symbol", "") or gene),
        "element": element,
        "levels": {},
        "states": {},
    }
    fitted: set[float] = set()
    clamped: tuple[float, float] | None = None
    for cell in cells:
        ctx = {"cell_type": cell}
        b = bridge.parameterize(
            small, ctx, build=True, rates=rates, borrowed_rate=borrowed, rate_provenance=prov or {},
            borrowed_provenance=lent or {},
        )  # fmt: skip
        keys = [k for k in b.strengths if k[1] == gene]
        if len(keys) != 1:
            return None
        (key,) = keys
        fitted.add(round(b.strengths[key], 12))
        for r in b.module.rules:
            if r.source == element and r.target == gene:
                r.through = element
        (state,) = [s for r, s, _w in b.module.reader_gate(ctx) if r.target == gene]
        out["states"][cell] = state
        out["levels"][cell] = _level(b.module, gene, ctx, {element: bridge.INTACT})
        if clamped is None:
            ungated = copy.deepcopy(b.module)
            for r in ungated.rules:
                r.through = ""
            clamped = (
                _level(ungated, gene, ctx, {element: bridge.INTACT}),
                _level(ungated, gene, ctx, {element: bridge.REMOVED}),
            )
            out["rho"] = b.observations[key].rho
            out["strength"] = b.strengths[key]
            out["tier"] = b.tier.get(gene, "")
    if clamped is None or len(fitted) != 1:
        return None
    present, absent = clamped
    out["L_open"], out["L_held"] = present, absent
    out["m"] = math.log2(present / absent)
    out["closed_form"] = -math.log2(out["rho"])
    out["action"] = "activates" if out["rho"] < 1.0 else "inhibits"
    scale = max(abs(present), abs(absent))
    out["gate_mismatch"] = [
        c
        for c in cells
        if abs(out["levels"][c] - (present if out["states"][c] == "open" else absent))
        > GATE_TOLERANCE * scale
    ]
    return out


# ==== the measured RNA and the score =====================================================================


def load_quant(path: Path) -> tuple[dict[str, float], dict[str, int]]:
    """gene_id without version -> TPM, from one ENCODE RSEM gene quantification TSV, with counts.

    Columns are found by name (`gene_id`, `TPM`); a gene_id on several rows (after dropping _PAR_Y rows)
    is ambiguous and dropped.
    """
    rows: dict[str, list[float]] = defaultdict(list)
    par = 0
    with Path(path).open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        gi, ti = header.index("gene_id"), header.index("TPM")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            gid = f[gi]
            if gid.endswith("_PAR_Y"):
                par += 1
                continue
            rows[gid.split(".")[0]].append(float(f[ti]))
    unique = {g: v[0] for g, v in rows.items() if len(v) == 1}
    return unique, {"rows": sum(len(v) for v in rows.values()), "par_y_rows": par, "unique": len(unique)}


def mean_tpm(tables: list[dict[str, float]]) -> dict[str, float]:
    """gene_id -> mean TPM over the replicate tables, for gene_ids on every one of them."""
    common = set.intersection(*(set(t) for t in tables)) if tables else set()
    return {g: sum(t[g] for t in tables) / len(tables) for g in common}


def symbol_gene_ids(path: Path = GENE_TABLE) -> dict[str, str]:
    """symbol -> gene_id without version, for symbols naming exactly one gene_id in the gene table."""
    ids: dict[str, set[str]] = defaultdict(set)
    with Path(path).open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        gi, si = header.index("gene_id"), header.index("symbol")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[gi].endswith("_PAR_Y"):
                continue
            ids[f[si]].add(f[gi].split(".")[0])
    return {s: next(iter(v)) for s, v in ids.items() if len(v) == 1}


def _y(tpm: float) -> float:
    return math.log2(tpm + PSEUDOCOUNT)


def _centre(xs: list[float]) -> list[float]:
    mu = sum(xs) / len(xs)
    return [x - mu for x in xs]


def assemble(
    genes: list[dict], tpm: dict[str, dict[str, float]], ids: dict[str, str], cells: tuple[str, ...] = CELLS
) -> tuple[list[dict], Counter]:
    """One row per population gene: D over the cells and the centred predictions; exclusions by reason.

    `tpm` is cell -> gene_id -> TPM (the mean over the chosen files). Applies population rules (c) and
    (d); (a) and (b) were applied when the genes were run.
    """
    out, why = [], Counter()
    for g in genes:
        opened = [g["states"][c] == "open" for c in cells]
        if all(opened) or not any(opened):
            why["element_open_in_all_or_none_of_the_cells"] += 1
            continue
        gid = ids.get(g["symbol"])
        if gid is None:
            why["symbol_not_on_exactly_one_gene_id"] += 1
            continue
        if any(gid not in tpm[c] for c in cells):
            why["gene_id_not_on_exactly_one_row_of_every_table"] += 1
            continue
        y = [_y(tpm[c][gid]) for c in cells]
        m = g["m"]
        cc = g.get("compiled_cell", "")
        gated = [m if o else 0.0 for o in opened]
        ungated = [m if c == cc else 0.0 for c in cells]
        unsigned = [abs(m) if o else 0.0 for o in opened]
        as_compiled = [m if (c == cc and o) else 0.0 for c, o in zip(cells, opened, strict=True)]
        out.append(
            {
                **{k: g[k] for k in ("gene", "symbol", "element", "action", "m", "rho", "tier")},
                "chrom": g.get("chrom", ""),
                "compiled_cell": cc,
                "gene_id": gid,
                "open": opened,
                "states": [g["states"][c] for c in cells],
                "D": _centre(y),
                "P_gated": _centre(gated),
                "P_ungated": _centre(ungated),
                "P_unsigned": _centre(unsigned),
                "P_compiled_cell_gated": _centre(as_compiled),
            }
        )
    return out, why


def _gene_sums(rows: list[dict], key: str) -> list[float]:
    return [sum(2.0 * p * d - p * p for p, d in zip(r[key], r["D"], strict=True)) for r in rows]


def dmse(rows: list[dict], key: str) -> float:
    """Delta MSE of one prediction over every pair of the rows (mean of 2 P D - P^2)."""
    n = sum(len(r["D"]) for r in rows)
    return sum(_gene_sums(rows, key)) / n if n else float("nan")


def bootstrap(
    rows: list[dict], keys: tuple[str, ...], resamples: int, seed: int
) -> dict[str, tuple[float, float]]:
    """95% percentile intervals of Delta MSE for each key and of gated minus ungated, from ONE set of
    gene resamples (each gene carries its thirteen pairs)."""
    import numpy as np

    sums = {k: np.array(_gene_sums(rows, k)) for k in keys}
    n = np.array([len(r["D"]) for r in rows], dtype=float)
    rng = np.random.default_rng(seed)
    stats = {k: np.empty(resamples) for k in keys}
    for i in range(resamples):
        idx = rng.integers(0, len(rows), len(rows))
        tot = n[idx].sum()
        for k in keys:
            stats[k][i] = sums[k][idx].sum() / tot
    out: dict[str, tuple[float, float]] = {}
    for k, v in stats.items():
        lo, hi = np.percentile(v, [2.5, 97.5])
        out[k] = (float(lo), float(hi))
    lo, hi = np.percentile(stats["P_gated"] - stats["P_ungated"], [2.5, 97.5])
    out["gated_minus_ungated"] = (float(lo), float(hi))
    return out


def shuffled_reader(rows: list[dict], count: int, seed: int) -> list[float]:
    """Delta MSE(gated) with the open/held rows permuted among genes, each gene keeping its own m."""
    import numpy as np

    rng = np.random.default_rng(seed)
    opened = np.array([r["open"] for r in rows], dtype=float)
    m = np.array([r["m"] for r in rows])
    d = np.array([r["D"] for r in rows])
    out = []
    for _ in range(count):
        s = opened[rng.permutation(len(rows))] * m[:, None]
        p = s - s.mean(axis=1, keepdims=True)
        out.append(float((2.0 * p * d - p * p).sum() / d.size))
    return out


def reading(kind: str, **numbers: Any) -> str:
    """The registered sentence for an outcome, numbers and the ungated clause filled, nothing else."""
    return OUTCOMES[kind].format(**numbers)


def decide(lo: float, hi: float, p: float) -> str:
    """DECISION, step 3."""
    if lo > 0:
        return "BEATS_AND_CONTROL" if p <= CONTROL_ALPHA else "BEATS_NOT_CONTROL"
    if hi < 0:
        return "WORSE"
    return "NO_DIFFERENCE"


def versus(lo: float, hi: float) -> str:
    """DECISION, step 4."""
    return UNGATED_CLAUSES["closer" if lo > 0 else "further" if hi < 0 else "same"]


def score(rows: list[dict], cells: tuple[str, ...] = CELLS) -> dict[str, Any]:
    """The registered comparison and the registered descriptives over an assembled population."""
    genes = len(rows)
    n = sum(len(r["D"]) for r in rows)
    if genes < MIN_GENES:
        return {
            "outcome": "NOT_TESTABLE",
            "reading": reading("NOT_TESTABLE", n=genes, min=MIN_GENES),
            "genes": genes,
        }
    keys = ("P_gated", "P_ungated", "P_unsigned", "P_compiled_cell_gated")
    dm = {k: dmse(rows, k) for k in keys}
    ci = bootstrap(rows, keys, BOOTSTRAP["resamples"], BOOTSTRAP["seed"])
    perm = shuffled_reader(rows, PERMUTATIONS["count"], PERMUTATIONS["seed"])
    p_ctrl = (1 + sum(1 for x in perm if x >= dm["P_gated"])) / (1 + len(perm))
    lo, hi = ci["P_gated"]
    ulo, uhi = ci["gated_minus_ungated"]
    kind = decide(lo, hi, p_ctrl)
    nums = {
        "d": f"{dm['P_gated']:+.4g}",
        "lo": f"{lo:+.4g}",
        "hi": f"{hi:+.4g}",
        "n": f"{n:,}",
        "g": f"{genes:,}",
        "p": f"{p_ctrl:.3g}",
        "versus": versus(ulo, uhi),
        "u": f"{dm['P_gated'] - dm['P_ungated']:+.4g}",
        "ulo": f"{ulo:+.4g}",
        "uhi": f"{uhi:+.4g}",
    }
    d_all = [x for r in rows for x in r["D"]]
    p_all = [x for r in rows for x in r["P_gated"]]
    nz = [(p, d) for p, d in zip(p_all, d_all, strict=True) if p != 0]
    by_cell = {}
    for j, c in enumerate(cells):
        pairs = [(r["P_gated"][j], r["D"][j]) for r in rows]
        by_cell[c] = {
            "open": sum(1 for r in rows if r["open"][j]),
            "delta_mse_gated": sum(2 * p * d - p * p for p, d in pairs) / len(pairs),
            "mean_D": sum(d for _, d in pairs) / len(pairs),
        }

    def sub(rs: list[dict]) -> dict[str, Any]:
        return {"genes": len(rs), "delta_mse_gated": dmse(rs, "P_gated") if rs else None}

    return {
        "outcome": kind,
        "reading": reading(kind, **nums),
        "limitations": list(LIMITATIONS),
        "genes": genes,
        "pairs": n,
        "delta_mse": dm["P_gated"],
        "ci95": list(ci["P_gated"]),
        "mse_baseline": sum(x * x for x in d_all) / n,
        "mse_gated": sum((d - p) ** 2 for p, d in zip(p_all, d_all, strict=True)) / n,
        "control_p": p_ctrl,
        "control_delta_mse_median": rr._median(perm),
        "control_delta_mse_95th_percentile": sorted(perm)[int(0.95 * len(perm)) - 1],
        "ungated_as_compiled": {
            "delta_mse": dm["P_ungated"],
            "ci95": list(ci["P_ungated"]),
            "gated_minus_ungated": dm["P_gated"] - dm["P_ungated"],
            "gated_minus_ungated_ci95": list(ci["gated_minus_ungated"]),
            "clause": nums["versus"],
        },
        "descriptives": {
            "by_action": {a: sub([r for r in rows if r["action"] == a]) for a in ("activates", "inhibits")},
            "unsigned_openness": {"delta_mse": dm["P_unsigned"], "ci95": list(ci["P_unsigned"])},
            "compiled_cell_gated": {
                "delta_mse": dm["P_compiled_cell_gated"],
                "ci95": list(ci["P_compiled_cell_gated"]),
            },
            "spearman_P_gated_D_nonzero": rr._spearman([p for p, _ in nz], [d for _, d in nz]),
            "by_cell": by_cell,
            "states_by_cell": {c: dict(Counter(r["states"][j] for r in rows)) for j, c in enumerate(cells)},
            "compiled_cell_is_a_reader_cell": sum(1 for r in rows if r["compiled_cell"] in cells),
        },
    }
