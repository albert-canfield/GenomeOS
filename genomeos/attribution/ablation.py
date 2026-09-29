# SPDX-License-Identifier: AGPL-3.0-or-later
"""The AlphaGenome ablation: which of the project's conclusions survive on experimental evidence alone.

Item 13 C4's diagnostic, registered 2026-09-28 before it was run. Every AlphaGenome-derived input is
removed (`ABLATED`): the all-element deletion sweep and its per-element response cache, the sample
target runs, and the predicted links and rules compiled from them. The model is never called; this
module reads the compiled programs, the measured assay files and the recorded results.

Three families of conclusions are censused (`FAMILIES`) and each is classified by one rule
(`CLASSIFICATION_RULE`) into `SURVIVES`, `DOES_NOT_SURVIVE` or `NEVER_MODEL_DEPENDENT`, with one
reason for every conclusion that does not survive (`REASONS`).

The rule is strict on purpose. A predicted link survives only when an experiment of the same endpoint
- an endogenous CRISPRi perturbation of the same element and the same gene - establishes the same
direction. A reporter tile, a transgenic embryo, a saturation-mutagenesis base or an eQTL can sit on
the same element and still say nothing about which gene it regulates (review S4: activity,
perturbation, contact and association never substitute for one another), so they are counted beside a
link that does not survive and never make it survive.
"""

from __future__ import annotations

import bisect
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from genomeos.attribution import holdout as ho
from genomeos.attribution import measured as ms
from genomeos.results import RESULTS_DIR

# ==================================================================================================
# The registration (2026-09-28, item 13 C4, lane-c4), fixed before the classification was applied.
# ==================================================================================================
REGISTERED = "2026-09-28"
SURVIVES = "survives"
DOES_NOT_SURVIVE = "does_not_survive"
NEVER_MODEL_DEPENDENT = "was_never_model_dependent"
CLASSES = (SURVIVES, DOES_NOT_SURVIVE, NEVER_MODEL_DEPENDENT)
ABOUT_THE_MODEL = "about_the_model"
CONTRADICTED = "contradicted"
INCONCLUSIVE = "inconclusive"
NEVER_MEASURED = "never_measured"
REASONS = {
    ABOUT_THE_MODEL: "the statement's subject is a model output (its value, its gain or its agreement)",
    CONTRADICTED: "the same endpoint measured the same element and gene: a well-powered null or a "
    "significant change in the opposite direction",
    INCONCLUSIVE: "the same endpoint measured it without deciding: an underpowered null, a missing "
    "effect, or an experimental arm below its own floor or with an interval across zero",
    NEVER_MEASURED: "no measurement of the same endpoint: the element and gene were never perturbed "
    "together (for a rule: not in the rule's cell), or the result has no experimental arm",
}
ABLATED = (
    "the all-element deletion sweep (results enhancer_targets_all_<chrom> and their local tables)",
    "its per-element response cache (data/knowledge/alphagenome/elements, all_elements)",
    "the sample target runs of the same model (constrained_targets_<chrom>, enhancer_targets_<chrom>)",
    "the compiled predicted layer: every `element` block without `_measured` and every `rule` with "
    "`evidence: predicted` in data/knowledge/compiled/noncoding_<chrom>.bio",
    "any value computed from these",
)
#: an input path containing one of these is model output
ABLATED_PATH_MARKS = ("alphagenome", "enhancer_targets", "constrained_targets", "knowledge/compiled/")
CLASSIFICATION_RULE = {
    NEVER_MODEL_DEPENDENT: "no ablated input enters the stated value or selects the units it is "
    "computed on. A model file read only for a quantity the conclusion does not state (chromosome names "
    "and lengths; a column the stated figure does not use) does not count, and is named. For a "
    "per-element label, the label's value is classified; that the element carries a block at all "
    "because the model named a target is counted separately as placement",
    SURVIVES: "an ablated input enters it, and with every ablated input removed the same statement is "
    "established on experimental evidence alone. A predicted link (element, gene, direction): a CRISPRi "
    "pair on the same element under measured.RECIPROCAL_OVERLAP, on the same gene, significant in the "
    "stated direction (activates: significant_decrease; inhibits: significant_increase), in any cell "
    "and either benchmark file. A compiled rule: the same, and in the rule's own `when` cell (R1). A "
    "result: an experimental arm of the same statistic is in the record, computed without ablated input, "
    "with the model figure's sign and a 95% interval (the result's own headline interval) excluding "
    "zero; the surviving figure is the arm's, never the model's",
    DOES_NOT_SURVIVE: "otherwise, with one reason: " + "; ".join(f"{k} ({v})" for k, v in REASONS.items()),
    "other_endpoints": "reporter activity, VISTA, saturation mutagenesis and GTEx association are other "
    "endpoints: they never make a link survive (S4); where they touch a link that does not survive they "
    "are counted beside it",
}
FAMILIES = (
    "headline results: the six rebuilt in data/results/manifest_headlines.json",
    "per-element labels: the five axes (and class, targets) of every element block of the 24 compiled "
    "programs, predicted and measured, and the region blocks",
    "compiled rules: every rule of the 24 programs, by evidence kind, after R1",
)
#: the registered reading of each headline result's derivation (read from its writer before the run)
HEADLINES: dict[str, dict[str, Any]] = {
    "node_containment_audit": {
        "conclusion": "the default node caller keeps an element's most-moved coding gene inside the "
        "element's own node 2.90 points above the published uniform control",
        "model_input": "value: the most-moved coding gene of each element in the deletion archive "
        "(scripts/node_containment_audit.py stage 1 reads data/knowledge/alphagenome/all_elements)",
        "arm": {
            "result": "node_containment_measured",
            "value": ["measured", "controls", "uniform", "excess_points"],
            "interval": ["measured", "controls", "uniform", "bootstrap_chromosomes", "ci95"],
            "what": "the same statistic against the same control on CRISPRi Regulated=TRUE pairs",
        },
        "model_value": ["modelled", "controls", "uniform", "excess_points"],
    },
    "node_containment_measured": {
        "conclusion": "on measured CRISPRi Regulated=TRUE pairs the node caller keeps the pair inside "
        "one node 5.89 points above the published uniform control",
        "model_input": "bookkeeping: stage 2 reuses stage 1's per-chromosome rows for the chromosome name "
        "and length only (measured_arm reads r['chrom'] and r['length']); the statistic reads CRISPRi "
        "pairs, ENCODE cCREs and GENCODE",
        "bookkeeping_only": True,
    },
    "constrained_unknown_targets": {
        "conclusion": "real-unknown blocks hold an element with a predicted coding mover less often than "
        "length-matched random windows (62.3% against 86.0%)",
        "model_input": "value: the mover is the deletion model's",
        "arm": {
            "result": "clause2_measured_arm",
            "reading": ["primary_reading", "outcome"],
            "what": "the same block-against-window question with CRISPRi-measured regulation of a "
            "coding gene",
        },
    },
    "therapeutic_benchmark": {
        "conclusion": "9 of 9 targets recovered, verdicts correct, top mechanisms defensible",
        "model_input": None,
    },
    "unknown_coverage": {
        "conclusion": "160,447 of 30,602,182 bp of the real unknown (0.52%) are read by any assay held",
        "model_input": "a column the figure does not state: organise.blocks reads the target runs for "
        "attributed_elements and targets; the blocks come from the budget and the spans from the "
        "assay files",
        "bookkeeping_only": True,
    },
    "crispri_published": {
        "conclusion": "the AlphaGenome deletion feature adds to ENCODE-rE2G on the published pair sets "
        "(training AUPRC 0.724 against 0.507; held-out DNase-only 0.476 to 0.639; K562 gain +0.136 "
        "[+0.081, +0.225])",
        "model_input": "value: the deletion feature",
        "subject_is_model_output": True,
    },
}
#: the axis values of a predicted element block that carry the model's link, and those about the model
REPRESSION_GROUP = "silencer|insulator_like|competing_promoter|unknown"
LINK_VALUES = {
    "activity": {"activates_target", "represses_target"},
    "target_relation": {"predicted_deletion_target", "nearest_tss_in_domain"},
    "molecular_role": {REPRESSION_GROUP},
}
ABOUT_MODEL_VALUES = {"evidence_status": {"predicted_model", "conflicting"}}
AXES = ("class", "targets", "origin", "molecular_role", "activity", "target_relation", "evidence_status")


def _norm_cell(c: str) -> str:
    return re.sub(r"[^a-z0-9]", "", c.lower())


# --- the same-endpoint evidence: CRISPRi pairs by element -----------------------------------------
class CrispriIndex:
    """Every valid CRISPRi pair (both benchmark files), searchable by the element it measured."""

    def __init__(self, pairs: list[ms.CrispriPair]) -> None:
        by: dict[str, list[ms.CrispriPair]] = defaultdict(list)
        for p in pairs:
            by[p.chrom].append(p)
        self.by = {c: sorted(v, key=lambda p: p.start) for c, v in by.items()}
        self.starts = {c: [p.start for p in v] for c, v in self.by.items()}
        self.reach = max((p.end - p.start for p in pairs), default=0)

    def of(self, chrom: str, start: int, end: int) -> list[ms.CrispriPair]:
        """Pairs whose element is this element under the measured layer's overlap rule."""
        items = self.by.get(chrom)
        if not items:
            return []
        st = self.starts[chrom]
        lo = bisect.bisect_left(st, start - self.reach)
        hi = bisect.bisect_left(st, end)
        return [p for p in items[lo:hi] if ms.measures(start, end, p.start, p.end)]


def link_class(
    pairs: list[ms.CrispriPair], gene: str, action: str, cell: str | None = None
) -> tuple[str, str | None, str]:
    """(class, reason, detail) of one predicted link against the pairs measured on its element.

    `cell` None classifies the link as an element label (any cell); a cell classifies a rule (R1)."""
    want = ms.DECREASE if action == "activates" else ms.INCREASE
    opposite = ms.INCREASE if action == "activates" else ms.DECREASE
    on_gene = [p for p in pairs if p.gene == gene]
    if cell is not None:
        on_gene = [p for p in on_gene if _norm_cell(p.cell) == _norm_cell(cell)]
    outcomes = {p.outcome for p in on_gene}
    if want in outcomes:
        return SURVIVES, None, f"CRISPRi {want} on {gene}"
    if ms.NULL_INFORMATIVE in outcomes or opposite in outcomes:
        return (
            DOES_NOT_SURVIVE,
            CONTRADICTED,
            "well-powered null" if ms.NULL_INFORMATIVE in outcomes else opposite,
        )
    if outcomes:
        return DOES_NOT_SURVIVE, INCONCLUSIVE, "/".join(sorted(outcomes))
    if pairs:
        return (
            DOES_NOT_SURVIVE,
            NEVER_MEASURED,
            "element screened, this gene not tested" + (" in this cell" if cell is not None else ""),
        )
    return DOES_NOT_SURVIVE, NEVER_MEASURED, "element never screened"


# --- the other endpoints, counted beside a link that does not survive -----------------------------
class OtherEndpoints:
    """lentiMPRA, VISTA, saturation mutagenesis and GTEx over one element, under the measured layer's
    overlap rule (GTEx: an indexed element sharing a base, same gene)."""

    def __init__(self, units: dict[str, tuple[ho.Unit, ...]]) -> None:
        tiles: dict[tuple[str, int, int], dict[str, float]] = defaultdict(dict)
        for name in [s for s in units if ho.assay_of(s) == "lentimpra"]:
            for u in units[name]:
                if u.value is not None:
                    tiles[(u.chrom, u.start, u.end)][u.cell] = u.value
        self.tiles = ho._Index(_Tile(c, s, e, v) for (c, s, e), v in tiles.items())
        self.vista = ho._Index(units.get("vista", ()))
        self.satmut = ho._Index(units.get("satmut", ()))
        self.gtex = ho._Index(u for u in units.get("gtex", ()) if u.outcome == "associated")
        self.ids_of: dict[str, set[str]] = defaultdict(set)
        for g in ho.genes():
            self.ids_of[g.symbol].add(g.gene_id)

    def at(self, chrom: str, start: int, end: int, gene: str) -> list[str]:
        out = []
        hits = [t for t in self.tiles.near(chrom, start, end) if ms.measures(start, end, t.start, t.end)]
        if hits:
            cells = {c for t in hits for c in t.values}
            if any(
                ms.reporter_label([t.values[c] for t in hits if c in t.values]) == ms.LABEL_ACTIVE
                for c in cells
            ):
                out.append("lentimpra_active")
        if any(
            v.outcome == "positive" and ms.measures(start, end, v.start, v.end)
            for v in self.vista.near(chrom, start, end)
        ):
            out.append("vista_positive")
        if any(b.outcome == "functional" for b in self.satmut.near(chrom, start, end)):
            out.append("satmut_functional_base")
        ids = self.ids_of.get(gene, set())
        if ids and any(h.gene_id in ids for h in self.gtex.near(chrom, start, end)):
            out.append("gtex_association_same_gene")
        return out


class _Tile:
    __slots__ = ("chrom", "start", "end", "values")

    def __init__(self, chrom: str, start: int, end: int, values: dict[str, float]) -> None:
        self.chrom, self.start, self.end, self.values = chrom, start, end, values


# --- family 2 and 3: the compiled programs ---------------------------------------------------------
def _groups(value: str) -> list[str]:
    return [g.strip() for g in value.split(", ") if g.strip()]


def axis_class(axis: str, value: str, measured: bool, link: tuple[str, str | None]) -> tuple[str, str | None]:
    """The class of one axis value of an element block under CLASSIFICATION_RULE."""
    if value in ABOUT_MODEL_VALUES.get(axis, ()):
        return DOES_NOT_SURVIVE, ABOUT_THE_MODEL
    if not measured and (axis == "targets" or value in LINK_VALUES.get(axis, ())):
        return link
    return NEVER_MODEL_DEPENDENT, None


def program_census(
    programs: list[Path],
    crispri: CrispriIndex,
    other: OtherEndpoints | None = None,
) -> dict[str, Any]:
    """Families 2 and 3 over the compiled programs, streamed one program at a time."""
    elements: Counter = Counter()
    element_detail: Counter = Counter()
    beside: dict[str, Counter] = defaultdict(Counter)
    axes: dict[str, dict[str, Counter]] = {
        "predicted": defaultdict(Counter),
        "measured": defaultdict(Counter),
    }
    rules: Counter = Counter()
    rule_detail: Counter = Counter()
    experimental: Counter = Counter()
    regions = Counter()
    predicted_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for path in programs:
        loci: dict[str, tuple[str, int, int]] = {}
        for kind, name, f in ho.program_blocks(path):
            if kind == "region":
                regions["blocks"] += 1
                regions["axis_values"] += sum(len(_groups(f[a])) for a in AXES if a in f)
                continue
            if kind == "element":
                m = ho.LOCUS_RE.match(f.get("locus", ""))
                if not m:
                    continue
                chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
                is_measured = name.endswith("_measured")
                link: tuple[str, str | None] = (NEVER_MODEL_DEPENDENT, None)
                if not is_measured:
                    loci[name] = (chrom, start, end)
                    predicted_intervals[chrom].append((start, end))
                    gene = f.get("targets", "").split(",")[0].strip()
                    action = "activates" if "activates_target" in f.get("activity", "") else "inhibits"
                    cls, reason, detail = link_class(crispri.of(chrom, start, end), gene, action)
                    link = (cls, reason)
                    elements[(cls, reason)] += 1
                    element_detail[detail if cls != SURVIVES else "survives"] += 1
                    if cls != SURVIVES and other is not None:
                        found = other.at(chrom, start, end, gene)
                        for k in found:
                            beside[reason][k] += 1
                        beside[reason]["any_other_endpoint" if found else "no_other_endpoint"] += 1
                side = "measured" if is_measured else "predicted"
                for axis in AXES:
                    if axis not in f:
                        continue
                    vals = [f[axis]] if axis == "targets" else _groups(f[axis])
                    for v in vals:
                        cls, reason = axis_class(axis, v, is_measured, link)
                        key = "gene" if axis == "targets" else v
                        axes[side][axis][(key, cls, reason)] += 1
                continue
            if kind == "rule":
                if f["evidence"] == "experimental":
                    experimental["rules"] += 1
                    experimental["heldout_marked" if f["heldout"] else "training"] += 1
                    rules[(NEVER_MODEL_DEPENDENT, None)] += 1
                    continue
                if name not in loci:
                    rules[("unmatched", None)] += 1
                    continue
                chrom, start, end = loci[name]
                cls, reason, detail = link_class(
                    crispri.of(chrom, start, end), f["gene"], f["action"], f["cell"]
                )
                rules[(cls, reason)] += 1
                rule_detail[detail if cls != SURVIVES else "survives"] += 1
    return {
        "elements": elements,
        "element_detail": element_detail,
        "beside": beside,
        "axes": axes,
        "rules": rules,
        "rule_detail": rule_detail,
        "experimental_rules": experimental,
        "regions": regions,
        "predicted_intervals": predicted_intervals,
    }


def placement(
    pairs: list[ms.CrispriPair],
    predicted_intervals: dict[str, list[tuple[int, int]]],
    ccres: dict[str, list[tuple[int, int]]],
) -> dict[str, Any]:
    """The measured regulatory links (significant decreases, as element-gene-cell) and whether a
    model-compiled element met the overlap rule to carry them; for the rest, whether an ENCODE cCRE
    would, with no model involved."""

    def index(iv: dict[str, list[tuple[int, int]]]):
        out = {}
        for c, v in iv.items():
            v = sorted(v)
            out[c] = (v, [s for s, _ in v], max((e - s for s, e in v), default=0))
        return out

    pi, ci = index(predicted_intervals), index(ccres)

    def meets(idx, chrom: str, s: int, e: int) -> bool:
        if chrom not in idx:
            return False
        v, st, reach = idx[chrom]
        lo = bisect.bisect_left(st, s - reach)
        hi = bisect.bisect_left(st, e)
        return any(ms.measures(a, b, s, e) for a, b in v[lo:hi])

    links = {(p.chrom, p.start, p.end, p.gene, p.cell, p.split) for p in pairs if p.outcome == ms.DECREASE}
    out = Counter()
    for chrom, s, e, _, _, split in links:
        if meets(pi, chrom, s, e):
            out[f"placed_{split}"] += 1
        elif meets(ci, chrom, s, e):
            out[f"not_placed_registry_would_{split}"] += 1
        else:
            out[f"not_placed_no_registry_element_{split}"] += 1
    return {"measured_links": len(links), **dict(sorted(out.items()))}


def ccre_intervals(results_dir: Path = RESULTS_DIR) -> dict[str, list[tuple[int, int]]]:
    out: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for chrom in ho.CHROMS:
        p = results_dir / f"ccres_{chrom}.bed.gz"
        if not p.exists():
            continue
        with gzip.open(p, "rt") as fh:
            for line in fh:
                if not line.startswith("#"):
                    f = line.split("\t", 3)
                    out[f[0]].append((int(f[1]), int(f[2])))
    return out


# --- family 1: the headline results ---------------------------------------------------------------
def _get(d: Any, path: list[str]) -> Any:
    for k in path:
        d = d[k]
    return d


def model_inputs(result: dict[str, Any]) -> list[str]:
    """The result's manifest inputs whose path marks them as model output."""
    rm = result.get("result_manifest") or {}
    paths = [str(i.get("path", "")) for i in rm.get("inputs", []) if isinstance(i, dict)]
    return sorted(p for p in paths if any(m in p for m in ABLATED_PATH_MARKS))


def headline_census(results_dir: Path = RESULTS_DIR) -> list[dict[str, Any]]:
    """Family 1: each headline result classified under the registered reading of its derivation."""
    heads = json.loads((results_dir / "manifest_headlines.json").read_text())
    names = [r["result"] for r in heads["rebuilt"]]
    out = []
    for name in names:
        reg = HEADLINES.get(name)
        res = json.loads((results_dir / f"{name}.json").read_text())
        found = model_inputs(res)
        row: dict[str, Any] = {
            "result": name,
            "conclusion": reg["conclusion"] if reg else None,
            "model_input_registered": reg.get("model_input") if reg else None,
            "model_inputs_in_manifest": len(found),
            "model_input_examples": found[:3],
            "carries_model_dependencies_block": bool(
                (res.get("result_manifest") or {}).get("model_dependencies")
            ),
        }
        if reg is None:
            row.update(cls="unregistered", reason=None, why="no registered reading; not classified")
        elif reg.get("model_input") is None:
            if found:
                row.update(
                    cls="registration_contradicted",
                    reason=None,
                    why="registered as model-free but the manifest names model inputs; not classified",
                )
            else:
                row.update(cls=NEVER_MODEL_DEPENDENT, reason=None, why="no model input in its manifest")
        elif reg.get("bookkeeping_only"):
            row.update(
                cls=NEVER_MODEL_DEPENDENT,
                reason=None,
                why="model file read for bookkeeping only: " + reg["model_input"],
            )
        elif reg.get("subject_is_model_output"):
            row.update(cls=DOES_NOT_SURVIVE, reason=ABOUT_THE_MODEL, why="its subject is the model's output")
        elif "arm" in reg:
            arm = reg["arm"]
            ar = json.loads((results_dir / f"{arm['result']}.json").read_text())
            if "reading" in arm:
                reading = _get(ar, arm["reading"])
                row["arm"] = {"result": arm["result"], "what": arm["what"], "reading": reading}
                if reading == "cannot_decide":
                    row.update(
                        cls=DOES_NOT_SURVIVE,
                        reason=INCONCLUSIVE,
                        why=f"its experimental arm ({arm['result']}) reads cannot_decide",
                    )
                else:
                    row.update(cls="unhandled_reading", reason=None, why=f"arm reading {reading!r}")
            else:
                v = float(_get(ar, arm["value"]))
                lo, hi = (float(x) for x in _get(ar, arm["interval"]))
                mv = float(_get(res, reg["model_value"])) if "model_value" in reg else None
                row["arm"] = {
                    "result": arm["result"],
                    "what": arm["what"],
                    "value": v,
                    "interval": [lo, hi],
                    "model_value": mv,
                }
                same_sign = mv is None or (v > 0) == (mv > 0)
                excludes_zero = lo > 0 or hi < 0
                if same_sign and excludes_zero and (lo > 0) == (v > 0):
                    row.update(
                        cls=SURVIVES,
                        reason=None,
                        why=f"its experimental arm {arm['result']} gives {v} [{lo}, {hi}]; the "
                        "surviving figure is the arm's",
                    )
                elif excludes_zero and not same_sign:
                    row.update(
                        cls=DOES_NOT_SURVIVE,
                        reason=CONTRADICTED,
                        why=f"arm {v} [{lo}, {hi}] has the other sign",
                    )
                else:
                    row.update(
                        cls=DOES_NOT_SURVIVE,
                        reason=INCONCLUSIVE,
                        why=f"arm interval [{lo}, {hi}] includes zero",
                    )
        else:
            row.update(cls=DOES_NOT_SURVIVE, reason=NEVER_MEASURED, why="model value, no experimental arm")
        out.append(row)
    return out


def registration() -> dict[str, Any]:
    return {
        "registered": REGISTERED,
        "classes": list(CLASSES),
        "reasons": REASONS,
        "ablated": list(ABLATED),
        "ablated_path_marks": list(ABLATED_PATH_MARKS),
        "rule": CLASSIFICATION_RULE,
        "families": list(FAMILIES),
        "headlines": HEADLINES,
        "link_values": {k: sorted(v) for k, v in LINK_VALUES.items()},
        "about_model_values": {k: sorted(v) for k, v in ABOUT_MODEL_VALUES.items()},
        "overlap_rule": f"measured.RECIPROCAL_OVERLAP = {ms.RECIPROCAL_OVERLAP}",
    }
