# SPDX-License-Identifier: AGPL-3.0-or-later
"""N1, frozen before any measurement: do the committed motif-to-gene predictions anticipate a factor's
knockdown response better than promoter proximity? (docs/ATTRIBUTION.md, "N1, the Perturb-seq response
test, registered before any measurement", lane-n1, 2026-10-01.)

Data: Replogle et al. 2022 (Cell 185:2559, doi 10.1016/j.cell.2022.05.013), genome-scale CRISPRi
Perturb-seq in K562 at day 8. The response comes from `K562_gwps_normalized_bulk_01.h5ad`, the
gemgroup-Z-normalized pseudobulk; the calibration gate's expression strata come from the non-targeting
rows of `K562_gwps_raw_bulk_01.h5ad`. Both are on Figshare+ (article 20029387, CC BY 4.0).

The claim is narrow. The element arm is the committed sample scan (300 elements per chromosome, 7,063
elements), so the test concerns those sampled predictions and at most 66 candidate factors before any
quality filter, not genome-wide attribution. A knockdown can act indirectly, so a response anticipated by
a motif is not a validated connection.

Every rule is a module constant, and no function takes a threshold as an argument, so a run applies each
rule once and has nothing to tune. `run` reads in a fixed order and logs every read:

1. the identities of both files, checked against the digest frozen at registration;
2. the calibration gate, on the 585 non-targeting rows only (`calibration_gate`); if it fails the run
   stops there, before any factor row is read, and records the failure;
3. the knockdown fields of the candidate rows only (`knockdown_status`); a field-semantics mismatch
   stops the run, and so do fewer than 30 eligible factors;
4. the expression rows of the eligible factors only; then the per-factor AUROCs (`factor_result`) and
   the paired comparison (`paired_analysis`).

The functions above `run` are pure: inputs are arguments and nothing is read from disk. `H5adPseudobulk`
is the only reader of the two files. It checks Figshare's md5 before its first read, reads datasets by
name from a fixed list, and never opens a statistic outside that list (the energy test, the
Anderson-Darling and Mann-Whitney counts, the leverage scores, any `.var` statistic).
"""

from __future__ import annotations

import hashlib
import json
import math
import random
from collections import Counter, defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Protocol

# --- the frozen rules ------------------------------------------------------------------------------------

RESPONDS_T = 3.0  # a gene responds when |T| >= 3, T = X * sqrt(num_cells_filtered)
NOMINAL_RATE = math.erfc(RESPONDS_T / math.sqrt(2.0))  # two-sided N(0, 1) tail at 3: 0.0026998
GATE_MAX_RATIO = 1.5  # observed |T| >= 3 share over the nominal share, overall and in every stratum
GATE_STRATA = 10  # deciles of control-row expression over the universe genes
GATE_MIN_EXPECTED = 20.0  # a stratum with fewer expected exceedances cannot be tested, and fails
KD_MAX_FOLD = 0.40  # eligible only at >= 60% on-target knockdown (fold_expr <= 0.40)
KD_MIN_EXPECTED_UMI = 10.0  # knockdown assessable only when n * control_expr >= 10 target UMIs
MIN_CELLS = 25  # num_cells_filtered of an eligible row
MIN_POSITIVES = 10  # a candidate has >= 10 genes with a nonzero score in each arm of its universe
NEIGHBOUR_BP = 10_000  # universe genes whose canonical TSS lies this close to the factor's are excluded
FLOOR_FACTORS = 30  # defined paired differences needed for an estimate; a floor, not shown power
MIN_GAIN = 0.02  # criterion 1: the equal-weight mean paired AUROC gain is at least +0.02
MIN_CLUSTERS = 10  # criterion 2 is assessable only with at least 10 TFClass clusters
BOOTSTRAP_B = 10_000
BOOTSTRAP_SEED = 20261001
LOWER_INDEX = int(0.025 * BOOTSTRAP_B)  # the 251st of 10,000 sorted means is the 95% lower bound
UPPER_INDEX = int(0.975 * BOOTSTRAP_B) - 1  # the 9,750th is the upper bound
SEMANTICS_REL_TOL = 1e-4

NON_TARGETING = "non-targeting"

SOURCE = {
    "article": "Figshare+ 20029387",
    "doi": "10.25452/figshare.plus.20029387.v1",
    "license": "CC BY 4.0",
    "paper_doi": "10.1016/j.cell.2022.05.013",
    "files": {
        "normalized": {
            "name": "K562_gwps_normalized_bulk_01.h5ad",
            "figshare_file_id": 35773217,
            "bytes": 374_587_922,
            "md5": "a3dfaa94ea8724217f5ecb1e14a5f0c8",
            "url": "https://ndownloader.figshare.com/files/35773217",
        },
        "raw": {
            "name": "K562_gwps_raw_bulk_01.h5ad",
            "figshare_file_id": 35774443,
            "bytes": 374_587_922,
            "md5": "4570b53c9d62ff6df281e622f0350060",
            "url": "https://ndownloader.figshare.com/files/35774443",
        },
    },
}

CONSTANTS = {
    "responds_t": RESPONDS_T,
    "nominal_rate": NOMINAL_RATE,
    "gate_max_ratio": GATE_MAX_RATIO,
    "gate_strata": GATE_STRATA,
    "gate_min_expected": GATE_MIN_EXPECTED,
    "kd_max_fold": KD_MAX_FOLD,
    "kd_min_expected_umi": KD_MIN_EXPECTED_UMI,
    "min_cells": MIN_CELLS,
    "min_positives": MIN_POSITIVES,
    "neighbour_bp": NEIGHBOUR_BP,
    "floor_factors": FLOOR_FACTORS,
    "min_gain": MIN_GAIN,
    "min_clusters": MIN_CLUSTERS,
    "bootstrap_b": BOOTSTRAP_B,
    "bootstrap_seed": BOOTSTRAP_SEED,
    "lower_index": LOWER_INDEX,
    "upper_index": UPPER_INDEX,
    "semantics_rel_tol": SEMANTICS_REL_TOL,
}


def _finite(v: Any) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)


# --- identities --------------------------------------------------------------------------------------------


def parse_row(label: str) -> dict[str, Any]:
    """A pseudobulk row label, `<n>_<symbol>_<tss>_<ENSG>` or `<n>_non-targeting_...`, as fields."""
    parts = label.split("_")
    if NON_TARGETING in parts:
        return {"label": label, "non_targeting": True}
    if len(parts) < 4 or not parts[-1].startswith("ENSG"):
        return {"label": label, "non_targeting": False, "unparsed": True}
    return {
        "label": label,
        "non_targeting": False,
        "symbol": parts[1],
        "tss": "_".join(parts[2:-1]),
        "gene_id": parts[-1],
    }


def identity_digest(obs_index: Sequence[str], var_gene_id: Sequence[str]) -> str:
    """sha256 of the row labels and measured-gene IDs in file order, the file's identity at registration."""
    body = json.dumps({"obs_index": list(obs_index), "var_gene_id": list(var_gene_id)}, separators=(",", ":"))
    return hashlib.sha256(body.encode()).hexdigest()


def monomer_name(jaspar_name: str) -> str | None:
    """A JASPAR profile name as a factor: upper case; a heterodimer (`MAX::MYC`) is None, never a factor."""
    n = jaspar_name.strip().upper()
    return None if "::" in n else n


def resolve_factor(
    name: str,
    symbol_ids: dict[str, set[str]],
    perturbed_ids: set[str],
    obs_symbol_ids: dict[str, set[str]],
) -> dict[str, Any]:
    """The perturbed gene a JASPAR monomer name stands for, by one declared path.

    1. The name is a GENCODE v50 protein-coding gene name with exactly one gene ID: that ID, matched to the
       perturbation index by Ensembl ID, never by symbol. Two IDs: ambiguous, excluded.
    2. The name is no GENCODE v50 protein-coding gene name: the perturbation index's own symbol, when it
       names exactly one gene ID.
    A name GENCODE resolves is never re-matched by symbol, so a renamed gene cannot borrow another's rows.
    """
    ids = symbol_ids.get(name, set())
    if len(ids) > 1:
        return {"status": "ambiguous_gencode_name"}
    if len(ids) == 1:
        gid = next(iter(ids))
        if gid in perturbed_ids:
            return {"status": "resolved", "gene_id": gid, "via": "gencode_v50_name"}
        return {"status": "not_perturbed", "gene_id": gid, "via": "gencode_v50_name"}
    obs = obs_symbol_ids.get(name, set())
    if len(obs) == 1:
        return {"status": "resolved", "gene_id": next(iter(obs)), "via": "perturbation_index_symbol"}
    if len(obs) > 1:
        return {"status": "ambiguous_perturbation_symbol"}
    return {"status": "unresolved_name"}


# --- the frozen predictions ------------------------------------------------------------------------------


@dataclass
class Predictions:
    """Both arms' scores by Ensembl gene ID, and which genes each arm scanned.

    A gene in `promoter_scanned` but absent from `promoter[f]` scores 0 for f: its promoter was scanned and
    f is not in its frozen `requires` list (no hit at 85% of the matrix range, or a hit below the list's
    top 8 enriched factors). A gene not in `promoter_scanned` has no proximity score at all: missing, not
    zero. The same holds for elements.
    """

    promoter_scanned: set[str] = field(default_factory=set)
    element_attributed: set[str] = field(default_factory=set)
    promoter: dict[str, dict[str, float]] = field(default_factory=lambda: defaultdict(dict))
    element: dict[str, dict[str, float]] = field(default_factory=lambda: defaultdict(dict))
    dropped: Counter = field(default_factory=Counter)


def predictions(results: dict[str, dict], symbol_index: dict[str, dict[str, list[str]]]) -> Predictions:
    """Both arms from the committed `motifs_<chrom>` results.

    `symbol_index[chrom][symbol]` lists the protein-coding gene IDs of that name on that chromosome (GENCODE
    v50, the annotation `motifs.promoter_loci` used). The proximity score of factor f for gene g is f's hit
    score in g's promoter `requires` list (TSS +- 1,000 bp, canonical transcript). The attribution score is
    the largest of f's hit scores over the sampled elements whose `target` is g. An element without a
    target, or whose target does not map to exactly one gene ID, is attributed to nothing and counted.
    """
    pred = Predictions()
    for chrom in sorted(results):
        res = results[chrom]
        index = symbol_index.get(chrom, {})
        for sym, rec in res.get("genes", {}).items():
            ids = index.get(sym, [])
            if len(ids) != 1:
                pred.dropped["promoter_symbol_ambiguous" if ids else "promoter_symbol_unmapped"] += 1
                continue
            g = ids[0]
            pred.promoter_scanned.add(g)
            for r in rec.get("requires", []):
                f = r["factor"].upper()
                pred.promoter[f][g] = max(pred.promoter[f].get(g, 0.0), float(r["score"]))
        for e in res.get("elements", []):
            target = e.get("target")
            if not target:
                pred.dropped["element_without_target"] += 1
                continue
            ids = index.get(target, [])
            if len(ids) != 1:
                pred.dropped["element_target_ambiguous" if ids else "element_target_unmapped"] += 1
                continue
            g = ids[0]
            pred.element_attributed.add(g)
            for r in e.get("requires", []):
                f = r["factor"].upper()
                pred.element[f][g] = max(pred.element[f].get(g, 0.0), float(r["score"]))
    return pred


def shared_universe(measured: Iterable[str], pred: Predictions) -> dict[str, Any]:
    """The genes both arms can score and the file measures: measured, promoter scanned, element attributed.

    Every measured gene outside it is counted by the reason it is out; none is scored as a zero.
    """
    m = set(measured)
    ps, ea = pred.promoter_scanned, pred.element_attributed
    u = m & ps & ea
    return {
        "genes": sorted(u),
        "coverage": {
            "measured": len(m),
            "measured_with_scanned_promoter": len(m & ps),
            "measured_with_attributed_element": len(m & ea),
            "universe": len(u),
            "universe_share_of_measured": round(len(u) / len(m), 4) if m else None,
        },
        "exclusions": {
            "measured_without_scanned_promoter": len(m - ps),
            "measured_with_promoter_without_sampled_element": len((m & ps) - ea),
            "measured_with_element_without_promoter": len((m & ea) - ps),
            "scanned_promoter_not_measured": len(ps - m),
            "attributed_element_target_not_measured": len(ea - m),
        },
    }


def excluded_genes(gene_id: str, tss: dict[str, tuple[str, int]], universe: Iterable[str]) -> list[str]:
    """The perturbed gene itself and every universe gene whose canonical TSS lies within 10 kb of its own:
    CRISPRi can silence a neighbour (Replogle 2022, Figure S3), which is not regulation by the factor."""
    out = {gene_id}
    if gene_id in tss:
        chrom, pos = tss[gene_id]
        out |= {
            g for g in universe if g in tss and tss[g][0] == chrom and abs(tss[g][1] - pos) <= NEIGHBOUR_BP
        }
    return sorted(out)


def factor_scores(factor: str, pred: Predictions, genes: Iterable[str]) -> dict[str, dict[str, float]]:
    """The nonzero scores of one factor over the given genes, per arm; every other gene there scores 0."""
    gs = set(genes)
    return {
        "attribution": {g: s for g, s in sorted(pred.element.get(factor, {}).items()) if g in gs and s > 0},
        "proximity": {g: s for g, s in sorted(pred.promoter.get(factor, {}).items()) if g in gs and s > 0},
    }


PENDING = "pending measurement"


def candidate_ledger(
    jaspar_names: Iterable[str],
    symbol_ids: dict[str, set[str]],
    rows: Sequence[dict[str, Any]],
    pred: Predictions,
    universe: Sequence[str],
    tss: dict[str, tuple[str, int]],
    cluster_of: dict[str, str],
) -> list[dict[str, Any]]:
    """One line per distinct JASPAR name, with the outcome of every step decided from predictions and
    identities, and the measurement steps left pending.

    `rows` are `parse_row` records with their position in the file as `index`. `cluster_of` maps a factor
    to its TFClass unit (`motifs.family_unit`). Steps, in order: monomer, resolved, perturbed, predicted
    (>= 10 nonzero genes in each arm of the factor's universe), then knockdown eligibility and a defined
    AUROC, both pending.
    """
    perturbed: dict[str, list[dict]] = defaultdict(list)
    obs_symbol_ids: dict[str, set[str]] = defaultdict(set)
    for r in rows:
        if r.get("non_targeting") or r.get("unparsed"):
            continue
        perturbed[r["gene_id"]].append(r)
        obs_symbol_ids[r["symbol"].upper()].add(r["gene_id"])
    perturbed_ids = set(perturbed)
    out: list[dict[str, Any]] = []
    for jn in sorted({n.strip().upper() for n in jaspar_names}):
        line: dict[str, Any] = {"jaspar_name": jn, "steps": {}}
        out.append(line)
        name = monomer_name(jn)
        if name is None:
            line["steps"]["monomer"] = "excluded: heterodimer profile, not attributable to one knockdown"
            line["status"] = "excluded"
            continue
        line["steps"]["monomer"] = "passed"
        res = resolve_factor(name, symbol_ids, perturbed_ids, obs_symbol_ids)
        if res["status"] in ("ambiguous_gencode_name", "ambiguous_perturbation_symbol", "unresolved_name"):
            line["steps"]["resolved"] = f"excluded: {res['status']}"
            line["status"] = "excluded"
            continue
        line["gene_id"], line["via"] = res["gene_id"], res["via"]
        line["steps"]["resolved"] = f"passed ({res['via']})"
        if res["status"] == "not_perturbed":
            line["steps"]["perturbed"] = "excluded: no row in the genome-wide screen"
            line["status"] = "excluded"
            continue
        line["rows"] = [
            {"index": r["index"], "label": r["label"], "tss": r["tss"]} for r in perturbed[res["gene_id"]]
        ]
        line["steps"]["perturbed"] = f"passed ({len(line['rows'])} row(s))"
        excl = excluded_genes(res["gene_id"], tss, universe)
        skip = set(excl)
        own = [g for g in universe if g not in skip]
        sc = factor_scores(name, pred, own)
        pos = {arm: len(v) for arm, v in sc.items()}
        line["positives"] = pos
        if min(pos.values()) < MIN_POSITIVES:
            line["steps"]["predicted"] = (
                f"excluded: {pos['attribution']} attribution and {pos['proximity']} proximity genes with a "
                f"nonzero score in its universe; {MIN_POSITIVES} needed in each"
            )
            line["status"] = "excluded"
            continue
        line["steps"]["predicted"] = "passed"
        line["steps"]["knockdown_eligible"] = PENDING
        line["steps"]["auroc_defined"] = PENDING
        line.update(
            {
                "status": "candidate",
                "factor": name,
                "cluster": cluster_of.get(name, name),
                "excluded_genes": excl,
                "universe_genes": len(own),
                "scores": sc,
            }
        )
    return out


# --- the calibration gate --------------------------------------------------------------------------------


def expression_strata(expression: Sequence[float]) -> list[int]:
    """Deciles: genes ranked by control-row expression (ties by position), cut into GATE_STRATA groups whose
    sizes differ by at most one. Stratum 0 holds the least expressed genes."""
    n = len(expression)
    order = sorted(range(n), key=lambda j: (expression[j], j))
    strata = [0] * n
    for rank, j in enumerate(order):
        strata[j] = rank * GATE_STRATA // n
    return strata


def calibration_gate(
    control_x: Sequence[Sequence[float]], control_n: Sequence[float], raw_x: Sequence[Sequence[float]]
) -> dict[str, Any]:
    """Whether T behaves as N(0, 1) where nothing was perturbed, overall and at every expression level.

    Inputs are the 585 non-targeting rows only, over the universe genes: normalized X, num_cells_filtered,
    and raw X for the strata. A row with a non-finite or sub-1 cell count is left out and counted. A gene
    with any non-finite T or a non-finite expression cannot be calibrated: it leaves the universe before
    any factor row is read, and is listed. Expression is the cell-weighted mean of raw X over the usable
    rows. Pass: the share of |T| >= 3 is at most 1.5 times the nominal 0.0027 overall and in each of the
    ten deciles, and each decile expects at least 20 exceedances. Anything else fails, and a failed gate
    ends the run.
    """
    usable = [i for i, n in enumerate(control_n) if _finite(n) and n >= 1]
    k = len(control_x[0]) if control_x else 0
    out: dict[str, Any] = {
        "rows_given": len(control_n),
        "rows_used": len(usable),
        "rows_excluded_bad_cell_count": len(control_n) - len(usable),
        "nominal_rate": NOMINAL_RATE,
        "max_ratio": GATE_MAX_RATIO,
        "min_expected": GATE_MIN_EXPECTED,
    }
    if not usable or k == 0:
        return {**out, "passed": False, "reason": "no usable control row", "genes_uncalibrated": []}
    root = {i: math.sqrt(control_n[i]) for i in usable}
    total_n = sum(control_n[i] for i in usable)
    t_cols: list[list[float]] = []
    expr: list[float] = []
    bad: list[int] = []
    for j in range(k):
        col = [control_x[i][j] * root[i] for i in usable]
        e = sum(raw_x[i][j] * control_n[i] for i in usable) / total_n
        if not all(math.isfinite(t) for t in col) or not math.isfinite(e):
            bad.append(j)
            col, e = [], float("nan")
        t_cols.append(col)
        expr.append(e)
    badset = set(bad)
    good = [j for j in range(k) if j not in badset]
    out["genes_uncalibrated"] = bad
    out["genes_calibrated"] = len(good)
    if not good:
        return {**out, "passed": False, "reason": "no calibratable gene"}
    strata = expression_strata([expr[j] for j in good])

    def tally(cols: list[int]) -> dict[str, Any]:
        values = sum(len(t_cols[j]) for j in cols)
        exceed = sum(1 for j in cols for t in t_cols[j] if abs(t) >= RESPONDS_T)
        expected = NOMINAL_RATE * values
        share = exceed / values if values else None
        ratio = share / NOMINAL_RATE if share is not None else None
        ok = expected >= GATE_MIN_EXPECTED and ratio is not None and ratio <= GATE_MAX_RATIO
        return {
            "genes": len(cols),
            "values": values,
            "exceedances": exceed,
            "expected": round(expected, 3),
            "share": share,
            "ratio": ratio,
            "passed": ok,
        }

    overall = tally(good)
    per: list[dict[str, Any]] = []
    for s in range(GATE_STRATA):
        cols = [good[p] for p, st in enumerate(strata) if st == s]
        ex = [expr[j] for j in cols]
        per.append(
            {
                "stratum": s,
                "expression_min": min(ex) if ex else None,
                "expression_max": max(ex) if ex else None,
                **tally(cols),
            }
        )
    allt = [t for j in good for t in t_cols[j]]
    mean = sum(allt) / len(allt)
    var = sum((t - mean) ** 2 for t in allt) / len(allt)
    passed = overall["passed"] and all(p["passed"] for p in per)
    reason = None if passed else "observed |T| >= 3 share exceeds 1.5x nominal, or a stratum cannot be tested"
    return {
        **out,
        "overall": overall,
        "strata": per,
        "t_mean": mean,
        "t_variance": var,
        "passed": passed,
        "reason": reason,
    }


# --- knockdown eligibility and the multi-row rule ----------------------------------------------------------


def semantics_consistent(fold: float, pct: Any) -> bool:
    """The registered reading of the fields: fold_expr is perturbed over non-targeting mean expression of
    the target, and pct_expr the fractional change, fold - 1 (or the same in percent). A row where pct_expr
    is neither contradicts the reading."""
    if not _finite(pct):
        return False
    return math.isclose(pct, fold - 1.0, rel_tol=SEMANTICS_REL_TOL, abs_tol=1e-6) or math.isclose(
        pct, 100.0 * (fold - 1.0), rel_tol=SEMANTICS_REL_TOL, abs_tol=1e-4
    )


def knockdown_status(fields: dict[str, Any]) -> dict[str, Any]:
    """One row's eligibility from num_cells_filtered, control_expr, fold_expr and pct_expr, in this order:

    1. num_cells_filtered finite and >= 25, else `too_few_cells`;
    2. control_expr finite, > 0, and num_cells_filtered * control_expr >= 10 expected target UMIs, else
       `knockdown_unassessable` (a target too weakly expressed to measure its knockdown is ineligible,
       never assumed knocked down);
    3. fold_expr finite and >= 0, else `knockdown_unassessable`;
    4. pct_expr agrees with fold_expr - 1 (`semantics_consistent`), else `semantics_mismatch`, which stops
       the run;
    5. fold_expr <= 0.40, else `insufficient_knockdown`.
    """
    n, c = fields.get("num_cells_filtered"), fields.get("control_expr")
    f, p = fields.get("fold_expr"), fields.get("pct_expr")
    if not _finite(n) or n < MIN_CELLS:
        return {"eligible": False, "reason": "too_few_cells"}
    if not _finite(c) or c <= 0 or n * c < KD_MIN_EXPECTED_UMI:
        return {"eligible": False, "reason": "knockdown_unassessable"}
    if not _finite(f) or f < 0:
        return {"eligible": False, "reason": "knockdown_unassessable"}
    if not semantics_consistent(f, p):
        return {"eligible": False, "reason": "semantics_mismatch"}
    if f > KD_MAX_FOLD:
        return {"eligible": False, "reason": "insufficient_knockdown"}
    return {"eligible": True, "reason": None}


def pooled_t(rows: Sequence[tuple[Sequence[float], float]]) -> list[float]:
    """The multi-row rule: every eligible row of a factor is pooled, weighted by its cells. X is a mean over
    each row's cells, so the pooled mean is sum(n_r * X_r) / sum(n_r) and T is that times sqrt(sum(n_r)).
    One row gives X * sqrt(n). No row is chosen by its response."""
    total = sum(n for _, n in rows)
    k = len(rows[0][0])
    root = math.sqrt(total)
    return [sum(x[j] * n for x, n in rows) / total * root for j in range(k)]


# --- per-factor AUROC and the paired comparison ------------------------------------------------------------


def auroc(scores: Sequence[float], labels: Sequence[bool]) -> tuple[float | None, str | None]:
    """Mann-Whitney AUROC with midranks for ties; None with its reason when undefined."""
    pos = sum(1 for x in labels if x)
    neg = len(labels) - pos
    if pos == 0:
        return None, "no_responders"
    if neg == 0:
        return None, "all_responders"
    if len(set(scores)) == 1:
        return None, "no_score_variation"
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    ranks = [0.0] * len(scores)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        mid = (i + j) / 2 + 1
        for q in range(i, j + 1):
            ranks[order[q]] = mid
        i = j + 1
    rank_sum = sum(r for r, x in zip(ranks, labels, strict=True) if x)
    return (rank_sum - pos * (pos + 1) / 2) / (pos * neg), None


def factor_result(candidate: dict[str, Any], genes: Sequence[str], t: Sequence[float]) -> dict[str, Any]:
    """One factor: responders by |T| >= 3 over its universe (the analysis universe minus its excluded
    genes), and the AUROC of each arm's score, a gene absent from an arm's scores scoring 0."""
    excl = set(candidate["excluded_genes"])
    keep = [j for j, g in enumerate(genes) if g not in excl and math.isfinite(t[j])]
    nonfinite = sum(1 for j, g in enumerate(genes) if g not in excl and not math.isfinite(t[j]))
    labels = [abs(t[j]) >= RESPONDS_T for j in keep]
    att = [candidate["scores"]["attribution"].get(genes[j], 0.0) for j in keep]
    prox = [candidate["scores"]["proximity"].get(genes[j], 0.0) for j in keep]
    a_att, why_att = auroc(att, labels)
    a_prox, why_prox = auroc(prox, labels)
    defined = a_att is not None and a_prox is not None
    return {
        "factor": candidate["factor"],
        "cluster": candidate["cluster"],
        "genes": len(keep),
        "genes_nonfinite_response": nonfinite,
        "responders": sum(labels),
        "auroc_attribution": a_att,
        "auroc_proximity": a_prox,
        "undefined_attribution": why_att,
        "undefined_proximity": why_prox,
        "difference": a_att - a_prox if defined else None,
    }


def undefined_reason(r: dict[str, Any]) -> str:
    """Why a factor's paired difference is undefined: no responders or all responders (both arms alike),
    else which arm's scores do not vary over its universe."""
    a, p = r["undefined_attribution"], r["undefined_proximity"]
    if a in ("no_responders", "all_responders"):
        return a
    if a and p:
        return "no_score_variation_both"
    return "no_score_variation_attribution" if a else "no_score_variation_proximity"


def cluster_bootstrap(clusters: dict[str, Sequence[float]]) -> list[float]:
    """Sorted bootstrap means: each replicate draws as many TFClass clusters as there are, with
    replacement, and averages every factor in the drawn clusters with equal weight."""
    names = sorted(clusters)
    sums = [sum(clusters[c]) for c in names]
    counts = [len(clusters[c]) for c in names]
    rng = random.Random(BOOTSTRAP_SEED)
    k = len(names)
    out: list[float] = []
    for _ in range(BOOTSTRAP_B):
        s, n = 0.0, 0
        for _ in range(k):
            c = rng.randrange(k)
            s += sums[c]
            n += counts[c]
        out.append(s / n)
    out.sort()
    return out


def paired_analysis(results: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """The primary comparison: the equal-weight mean over factors of AUROC(attribution) - AUROC(proximity).

    Every factor is counted in coverage, with the reason when its difference is undefined. Fewer than 30
    defined differences stop the analysis with no estimate. The two criteria are reported apart and never
    merged: the point gain against +0.02, and the 95% cluster-bootstrap lower bound against 0 (not
    assessable with fewer than 10 clusters).
    """
    defined = [r for r in results if r["difference"] is not None]
    why = Counter(undefined_reason(r) for r in results if r["difference"] is None)
    coverage = {
        "factors_analysed": len(results),
        "defined": len(defined),
        "undefined": len(results) - len(defined),
        "undefined_reasons": dict(sorted(why.items())),
        "floor": FLOOR_FACTORS,
    }
    if len(defined) < FLOOR_FACTORS:
        return {
            "status": "insufficient_coverage",
            "coverage": coverage,
            "estimate": None,
            "criteria": None,
            "note": f"{len(defined)} defined paired differences, below the floor of {FLOOR_FACTORS}; no "
            "estimate is reported and no rule is relaxed",
        }
    d = [r["difference"] for r in defined]
    est = sum(d) / len(d)
    clusters: dict[str, list[float]] = defaultdict(list)
    for r in defined:
        clusters[r["cluster"]].append(r["difference"])
    lower = upper = None
    if len(clusters) >= MIN_CLUSTERS:
        boot = cluster_bootstrap(clusters)
        lower, upper = boot[LOWER_INDEX], boot[UPPER_INDEX]
    return {
        "status": "analysed",
        "coverage": coverage,
        "estimate": est,
        "clusters": len(clusters),
        "interval_95": [lower, upper] if lower is not None else None,
        "criteria": {
            "point_gain": {"rule": f"estimate >= +{MIN_GAIN}", "value": est, "met": est >= MIN_GAIN},
            "lower_bound": {
                "rule": "95% cluster-bootstrap lower bound > 0",
                "value": lower,
                "met": None if lower is None else lower > 0,
                "assessable": lower is not None,
            },
        },
    }


# --- the run: gate first, then eligibility, then responses ------------------------------------------------


class Reader(Protocol):
    def identities(self) -> dict[str, dict[str, list[str]]]: ...

    def control_rows(self, rows: list[int], cols: list[int]) -> dict[str, Any]: ...

    def knockdown_fields(self, rows: list[int]) -> dict[int, dict[str, Any]]: ...

    def response_rows(self, rows: list[int], cols: list[int]) -> dict[int, list[float]]: ...


def run(plan: dict[str, Any], reader: Reader) -> dict[str, Any]:
    """Apply the registration once. `plan` is the registration's `plan`; `reader` serves the two files."""
    reads: list[dict[str, Any]] = []

    def stop(status: str, **extra: Any) -> dict[str, Any]:
        return {"status": status, "reads": reads, **extra}

    ident = reader.identities()
    reads.append({"step": "identities", "datasets": ["obs index", "var gene_id"], "files": sorted(ident)})
    norm = ident["normalized"]
    if identity_digest(norm["obs_index"], norm["var_gene_id"]) != plan["identity_digest"]:
        return stop("identity_mismatch", note="the normalized file's rows or genes differ from registration")
    if ident["raw"]["obs_index"] != norm["obs_index"] or ident["raw"]["var_gene_id"] != norm["var_gene_id"]:
        return stop(
            "identity_mismatch", note="the raw file's rows or genes differ from the normalized file's"
        )
    col_of = {g: j for j, g in enumerate(norm["var_gene_id"])}
    universe = list(plan["universe"])
    cols = [col_of[g] for g in universe]
    control = list(plan["control_rows"])
    c = reader.control_rows(control, cols)
    reads.append(
        {"step": "calibration gate", "rows": len(control), "row_kind": "non-targeting", "genes": len(cols)}
    )
    gate = calibration_gate(c["x"], c["n"], c["raw_x"])
    if not gate["passed"]:
        return stop("gate_failed", gate=gate, note="the run stops; no factor row is read")
    drop = set(gate["genes_uncalibrated"])
    genes = [g for j, g in enumerate(universe) if j not in drop]
    gcols = [col_of[g] for g in genes]
    cands = list(plan["candidates"])
    cand_rows = sorted({r["index"] for cand in cands for r in cand["rows"]})
    kd = reader.knockdown_fields(cand_rows)
    reads.append({"step": "knockdown eligibility", "rows": len(cand_rows), "row_kind": "candidate factor"})
    row_status = {i: knockdown_status(kd[i]) for i in cand_rows}
    mismatched = sorted(i for i, s in row_status.items() if s["reason"] == "semantics_mismatch")
    ledger = []
    eligible = []
    for cand in cands:
        rs = [{"index": r["index"], "label": r["label"], **row_status[r["index"]]} for r in cand["rows"]]
        ok = [r for r in rs if r["eligible"]]
        ledger.append({"factor": cand["factor"], "rows": rs, "knockdown_eligible": bool(ok)})
        if ok:
            eligible.append((cand, [r["index"] for r in ok]))
    if mismatched:
        return stop("semantics_check_failed", gate=gate, eligibility=ledger, rows_mismatched=mismatched)
    if len(eligible) < FLOOR_FACTORS:
        return stop(
            "insufficient_coverage",
            gate=gate,
            eligibility=ledger,
            note=f"{len(eligible)} knockdown-eligible factors, below the floor of {FLOOR_FACTORS}; no "
            "response row is read",
        )
    resp_rows = sorted({i for _, idx in eligible for i in idx})
    x = reader.response_rows(resp_rows, gcols)
    reads.append(
        {"step": "responses", "rows": len(resp_rows), "row_kind": "eligible factor", "genes": len(gcols)}
    )
    results = []
    for cand, idx in eligible:
        t = pooled_t([(x[i], kd[i]["num_cells_filtered"]) for i in idx])
        results.append({**factor_result(cand, genes, t), "rows_pooled": len(idx)})
    analysis = paired_analysis(results)
    return {
        "status": analysis["status"],
        "reads": reads,
        "gate": gate,
        "eligibility": ledger,
        "factors": results,
        "analysis": analysis,
        "genes_analysed": len(genes),
    }


# --- the one reader of the two files ----------------------------------------------------------------------


def file_digests(path: Path) -> tuple[str, str]:
    """(md5, sha256) of a file, read once."""
    md5, sha = hashlib.md5(usedforsecurity=False), hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            md5.update(chunk)
            sha.update(chunk)
    return md5.hexdigest(), sha.hexdigest()


class H5adPseudobulk:
    """Serves `run` from the two downloaded pseudobulk files, by dataset name from a fixed list.

    Construction checks each file's md5 against Figshare's before anything is read and keeps its sha256.
    Rows are read whole from `X` (a dense float32 matrix) and only the universe columns are kept.
    """

    OBS_INDEX = "obs/gene_transcript"
    VAR_INDEX = "var/gene_id"
    KNOCKDOWN = ("num_cells_filtered", "control_expr", "fold_expr", "pct_expr")

    def __init__(self, normalized: Path, raw: Path, expected: dict[str, dict] | None = None):
        expected = expected or SOURCE["files"]
        self.paths = {"normalized": Path(normalized), "raw": Path(raw)}
        self.sha256: dict[str, str] = {}
        for kind, p in self.paths.items():
            md5, sha = file_digests(p)
            if md5 != expected[kind]["md5"]:
                raise ValueError(f"{p}: md5 {md5} is not Figshare's {expected[kind]['md5']}")
            self.sha256[kind] = sha
        self.datasets_read: list[str] = []

    def _read(self, kind: str, name: str, rows: list[int] | None = None):
        import h5py

        self.datasets_read.append(f"{kind}:{name}" + ("" if rows is None else f"[{len(rows)} rows]"))
        with h5py.File(self.paths[kind], "r") as h:
            ds = h[name]
            return ds[()] if rows is None else ds[sorted(rows)]

    def identities(self) -> dict[str, dict[str, list[str]]]:
        out = {}
        for kind in self.paths:
            dec = [x.decode() if isinstance(x, bytes) else str(x) for x in self._read(kind, self.OBS_INDEX)]
            var = [x.decode() if isinstance(x, bytes) else str(x) for x in self._read(kind, self.VAR_INDEX)]
            out[kind] = {"obs_index": dec, "var_gene_id": var}
        return out

    def _rows(self, kind: str, rows: list[int], cols: list[int]) -> dict[int, list[float]]:
        block = self._read(kind, "X", rows)
        return {i: [float(v) for v in block[p][cols]] for p, i in enumerate(sorted(rows))}

    def control_rows(self, rows: list[int], cols: list[int]) -> dict[str, Any]:
        xs = self._rows("normalized", rows, cols)
        raw = self._rows("raw", rows, cols)
        n = self._read("normalized", "obs/num_cells_filtered", rows)
        n_of = {i: float(v) for i, v in zip(sorted(rows), n, strict=True)}
        return {"x": [xs[i] for i in rows], "raw_x": [raw[i] for i in rows], "n": [n_of[i] for i in rows]}

    def knockdown_fields(self, rows: list[int]) -> dict[int, dict[str, Any]]:
        cols = {f: self._read("normalized", f"obs/{f}", rows) for f in self.KNOCKDOWN}
        return {i: {f: float(cols[f][p]) for f in self.KNOCKDOWN} for p, i in enumerate(sorted(rows))}

    def response_rows(self, rows: list[int], cols: list[int]) -> dict[int, list[float]]:
        return self._rows("normalized", rows, cols)


# =========================================================================================================
# Amendment 1 (2026-10-01, lane-n1), additive: everything above is the code the registration of f9a9130
# froze (b359867) and is unchanged. data/results/n1_registration_amendment_1.json states each change and why.
#
# N1 is an exploratory comparison against an operational response label (|T| >= 3). T = X * sqrt(n) is a
# proposed statistic, not a calibrated test, and passing the gate does not validate per-gene biological
# responses. The knockdown fields' units are not documented for this file, so the knockdown rule is
# unsupported and `run_v2` stops before applying it (ELIGIBILITY_SUPPORTED). The functions after that stop
# are frozen and tested so that a later, documented amendment changes one constant, not the analysis.
# =========================================================================================================

AMENDMENT = "n1_registration_amendment_1"
ELIGIBILITY_SUPPORTED = (
    False  # item 2: control_expr's units are undocumented, so the knockdown rule is unsupported
)

CONSTANTS_V2 = {
    **CONSTANTS,
    "amendment": 1,
    "auroc": "auroc_v2: ties count half; undefined only when a class is empty",
    "primary_auroc": "stratified by the gate's control-expression deciles, strata weighted by pairs",
    "secondary_auroc": "unstratified auroc_v2, reported beside the primary, no criterion applied",
    "eligibility_supported": ELIGIBILITY_SUPPORTED,
    "reader": "H5adPseudobulkV2: X decoded at the selected columns only",
}


def auroc_v2(scores: Sequence[float], labels: Sequence[bool]) -> tuple[float | None, str | None]:
    """AUROC as the Mann-Whitney probability that a responder outscores a non-responder, ties counting half
    (the standard definition; a constant score gives 0.5). Undefined only when a class is empty."""
    pos = sum(1 for x in labels if x)
    neg = len(labels) - pos
    if pos == 0:
        return None, "no_responders"
    if neg == 0:
        return None, "all_responders"
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    ranks = [0.0] * len(scores)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        mid = (i + j) / 2 + 1
        for q in range(i, j + 1):
            ranks[order[q]] = mid
        i = j + 1
    rank_sum = sum(r for r, x in zip(ranks, labels, strict=True) if x)
    return (rank_sum - pos * (pos + 1) / 2) / (pos * neg), None


def stratified_auroc(
    scores: Sequence[float], labels: Sequence[bool], strata: Sequence[int]
) -> tuple[float | None, str | None, int]:
    """Item 4: AUROC with responders compared only to non-responders of the same control-expression decile.

    Each stratum holding both classes contributes its auroc_v2, weighted by its responder x non-responder
    pairs, so the result is the probability over within-stratum pairs. A stratum with one class contributes
    nothing. Undefined with no responders, with all responders, or when no stratum holds both classes.
    Returns (value, reason, strata used)."""
    pos = sum(1 for x in labels if x)
    if pos == 0:
        return None, "no_responders", 0
    if pos == len(labels):
        return None, "all_responders", 0
    groups: dict[int, list[int]] = defaultdict(list)
    for j, s in enumerate(strata):
        groups[s].append(j)
    num = den = 0.0
    used = 0
    for s in sorted(groups):
        idx = groups[s]
        lab = [labels[j] for j in idx]
        p = sum(lab)
        pairs = p * (len(lab) - p)
        if pairs == 0:
            continue
        a, _ = auroc_v2([scores[j] for j in idx], lab)
        num += a * pairs
        den += pairs
        used += 1
    if den == 0:
        return None, "no_stratum_with_both_classes", 0
    return num / den, None, used


def control_expression(raw_x: Sequence[Sequence[float]], control_n: Sequence[float]) -> list[float]:
    """The gate's per-gene control expression: the cell-weighted mean of raw X over the usable rows."""
    usable = [i for i, n in enumerate(control_n) if _finite(n) and n >= 1]
    total = sum(control_n[i] for i in usable)
    k = len(raw_x[0]) if raw_x else 0
    return [sum(raw_x[i][j] * control_n[i] for i in usable) / total for j in range(k)]


def factor_result_v2(
    candidate: dict[str, Any], genes: Sequence[str], t: Sequence[float], strata: Sequence[int]
) -> dict[str, Any]:
    """One factor under amendment 1: the primary AUROC is stratified by control-expression decile and the
    unstratified auroc_v2 is reported beside it. Both arms share the labels and strata, so an undefined
    factor has one reason for both."""
    excl = set(candidate["excluded_genes"])
    keep = [j for j, g in enumerate(genes) if g not in excl and math.isfinite(t[j])]
    labels = [abs(t[j]) >= RESPONDS_T for j in keep]
    st = [strata[j] for j in keep]
    att = [candidate["scores"]["attribution"].get(genes[j], 0.0) for j in keep]
    prox = [candidate["scores"]["proximity"].get(genes[j], 0.0) for j in keep]
    a_att, why, used = stratified_auroc(att, labels, st)
    a_prox, _, _ = stratified_auroc(prox, labels, st)
    u_att, _ = auroc_v2(att, labels)
    u_prox, _ = auroc_v2(prox, labels)
    defined = a_att is not None and a_prox is not None
    return {
        "factor": candidate["factor"],
        "cluster": candidate["cluster"],
        "genes": len(keep),
        "genes_nonfinite_response": sum(
            1 for j, g in enumerate(genes) if g not in excl and not math.isfinite(t[j])
        ),
        "responders": sum(labels),
        "strata_used": used,
        "auroc_attribution": a_att,
        "auroc_proximity": a_prox,
        "difference": a_att - a_prox if defined else None,
        "undefined_reason": None if defined else why,
        "unstratified": {
            "attribution": u_att,
            "proximity": u_prox,
            "difference": u_att - u_prox if u_att is not None and u_prox is not None else None,
        },
    }


def paired_analysis_v2(results: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """The comparison of amendment 1: the registered estimate, floor, cluster bootstrap and two separate
    criteria (`paired_analysis`, unchanged) over the factors whose stratified difference is defined; the
    coverage counts every factor by its one reason; the unstratified mean is reported, with no criterion."""
    defined = [r for r in results if r["difference"] is not None]
    out = paired_analysis(defined)
    why = Counter(r["undefined_reason"] for r in results if r["difference"] is None)
    out["coverage"] = {
        "factors_analysed": len(results),
        "defined": len(defined),
        "undefined": len(results) - len(defined),
        "undefined_reasons": dict(sorted(why.items())),
        "floor": FLOOR_FACTORS,
    }
    u = [r["unstratified"]["difference"] for r in defined if r["unstratified"]["difference"] is not None]
    out["secondary_unstratified"] = {"estimate": sum(u) / len(u) if u else None, "factors": len(u)}
    return out


def run_v2(plan: dict[str, Any], reader: Reader) -> dict[str, Any]:
    """Amendment 1's run: identities, then the gate, then a stop before the knockdown rule while it is
    unsupported. No candidate row is read. The analysis after the stop is `_analyse_v2`."""
    reads: list[dict[str, Any]] = []

    def stop(status: str, **extra: Any) -> dict[str, Any]:
        return {"status": status, "amendment": 1, "reads": reads, **extra}

    ident = reader.identities()
    reads.append({"step": "identities", "datasets": ["obs index", "var gene_id"], "files": sorted(ident)})
    norm = ident["normalized"]
    if identity_digest(norm["obs_index"], norm["var_gene_id"]) != plan["identity_digest"]:
        return stop("identity_mismatch", note="the normalized file's rows or genes differ from registration")
    if ident["raw"]["obs_index"] != norm["obs_index"] or ident["raw"]["var_gene_id"] != norm["var_gene_id"]:
        return stop(
            "identity_mismatch", note="the raw file's rows or genes differ from the normalized file's"
        )
    col_of = {g: j for j, g in enumerate(norm["var_gene_id"])}
    universe = list(plan["universe"])
    cols = [col_of[g] for g in universe]
    control = list(plan["control_rows"])
    c = reader.control_rows(control, cols)
    reads.append(
        {"step": "calibration gate", "rows": len(control), "row_kind": "non-targeting", "genes": len(cols)}
    )
    gate = calibration_gate(c["x"], c["n"], c["raw_x"])
    if not gate["passed"]:
        return stop("gate_failed", gate=gate, note="the run stops; no factor row is read")
    drop = set(gate["genes_uncalibrated"])
    keep = [j for j in range(len(universe)) if j not in drop]
    genes = [universe[j] for j in keep]
    expr = control_expression(c["raw_x"], c["n"])
    strata = expression_strata([expr[j] for j in keep])
    if not ELIGIBILITY_SUPPORTED:
        return stop(
            "eligibility_unsupported",
            gate=gate,
            note="the knockdown rule rests on undocumented units of control_expr and fold_expr (amendment 1, "
            "item 2); the run stops before reading any candidate row",
        )
    return _analyse_v2(plan, reader, reads, gate, genes, [col_of[g] for g in genes], strata)


def _analyse_v2(
    plan: dict[str, Any],
    reader: Reader,
    reads: list[dict[str, Any]],
    gate: dict[str, Any],
    genes: list[str],
    gcols: list[int],
    strata: list[int],
) -> dict[str, Any]:
    """After the gate, only once a documented amendment sets ELIGIBILITY_SUPPORTED: the registered knockdown
    rules, the multi-row rule, factor_result_v2 and paired_analysis_v2."""
    cands = list(plan["candidates"])
    cand_rows = sorted({r["index"] for cand in cands for r in cand["rows"]})
    kd = reader.knockdown_fields(cand_rows)
    reads.append({"step": "knockdown eligibility", "rows": len(cand_rows), "row_kind": "candidate factor"})
    row_status = {i: knockdown_status(kd[i]) for i in cand_rows}
    mismatched = sorted(i for i, s in row_status.items() if s["reason"] == "semantics_mismatch")
    ledger, eligible = [], []
    for cand in cands:
        rs = [{"index": r["index"], "label": r["label"], **row_status[r["index"]]} for r in cand["rows"]]
        ok = [r for r in rs if r["eligible"]]
        ledger.append({"factor": cand["factor"], "rows": rs, "knockdown_eligible": bool(ok)})
        if ok:
            eligible.append((cand, [r["index"] for r in ok]))
    base = {"amendment": 1, "reads": reads, "gate": gate, "eligibility": ledger}
    if mismatched:
        return {"status": "semantics_check_failed", **base, "rows_mismatched": mismatched}
    if len(eligible) < FLOOR_FACTORS:
        return {
            "status": "insufficient_coverage",
            **base,
            "note": f"{len(eligible)} knockdown-eligible factors",
        }
    resp_rows = sorted({i for _, idx in eligible for i in idx})
    x = reader.response_rows(resp_rows, gcols)
    reads.append(
        {"step": "responses", "rows": len(resp_rows), "row_kind": "eligible factor", "genes": len(gcols)}
    )
    results = []
    for cand, idx in eligible:
        t = pooled_t([(x[i], kd[i]["num_cells_filtered"]) for i in idx])
        results.append({**factor_result_v2(cand, genes, t, strata), "rows_pooled": len(idx)})
    analysis = paired_analysis_v2(results)
    return {
        "status": analysis["status"],
        **base,
        "factors": results,
        "analysis": analysis,
        "genes_analysed": len(genes),
    }


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def freeze_problems(
    amendment: dict[str, Any], original_path: Path, module_path: Path, runner_path: Path
) -> list[str]:
    """Item 5: what differs from the freeze, checked before any measurement is opened; empty when nothing.

    The original registration's bytes, this module's and the v2 runner's sha256, and CONSTANTS_V2, against
    the amendment registration."""
    frozen = amendment.get("frozen_code", {})
    checks = [
        ("original registration", original_path, amendment.get("amends", {}).get("sha256")),
        ("analysis module", module_path, frozen.get("module", {}).get("sha256")),
        ("v2 runner", runner_path, frozen.get("runner", {}).get("sha256")),
    ]
    problems = [
        f"{what}: sha256 {sha256_file(p)} is not the registered {want}"
        for what, p, want in checks
        if sha256_file(p) != want
    ]
    if json.loads(json.dumps(CONSTANTS_V2)) != amendment.get("constants"):
        problems.append("constants: CONSTANTS_V2 differ from the registered constants")
    return problems


class H5adPseudobulkV2(H5adPseudobulk):
    """Item 5: as H5adPseudobulk, but X is decoded at the selected columns only, one row at a time
    (`X[i, sorted columns]`); whole rows are never decoded. The operating system may still read surrounding
    bytes from disk."""

    def _rows(self, kind: str, rows: list[int], cols: list[int]) -> dict[int, list[float]]:
        import h5py

        order = sorted(range(len(cols)), key=lambda q: cols[q])
        scols = [cols[q] for q in order]
        if len(set(scols)) != len(scols):
            raise ValueError("a column is selected twice")
        self.datasets_read.append(f"{kind}:X[{len(rows)} rows x {len(cols)} columns]")
        out: dict[int, list[float]] = {}
        with h5py.File(self.paths[kind], "r") as h:
            ds = h["X"]
            for i in sorted(rows):
                vals = ds[i, scols]
                row = [0.0] * len(cols)
                for p, q in enumerate(order):
                    row[q] = float(vals[p])
                out[i] = row
        return out
