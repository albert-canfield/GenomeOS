# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-inhibitor: the census of compiled inhibitor strengths above 1, answering the committed QUESTIONS.

    uv run python scripts/inhibitor_census_run.py [--chroms chr21] [--no-save] [--workers 4]

A CENSUS. Every question, definition, derived threshold and branch rule was committed first, in
`genomeos/attribution/inhibitor_census.py`, and this script computes them once. It re-derives the residue
and the 39 with `lane-residuemult`'s own code path (`combine2` + `saturate` + `scalefit`), reconciles the
set it finds against `residuemult_lookup.json`'s committed per-gene table on the `symbol` column, and then
counts the inhibiting elements of EVERY eligible multi-element gene so that the 39 can be read against the
population instead of in isolation.

NO RNA VALUE IS OPENED. Of `multrna_result.json`'s per-gene table only the `symbol` column is read; its
`P` and `D` columns are never touched. NOTHING IS FITTED, TUNED, SELECTED OR RE-ESTIMATED: the scale is a
constant from `ac5580b` through `scalefit.SCALE` and the multiplicative form has no parameter of its own.
No model request is sent and no key is read. This lane registers no test, proposes no change to any fit
and makes no claim about expression.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from multiprocessing import get_context
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import bridge  # noqa: E402
from genomeos.attribution import combine2 as c2  # noqa: E402
from genomeos.attribution import compile as cp  # noqa: E402
from genomeos.attribution import inhibitor_census as ic  # noqa: E402
from genomeos.attribution import readergate as rg  # noqa: E402
from genomeos.attribution import residuemult as rm  # noqa: E402
from genomeos.attribution import scalefit as sf  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.results import save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readergate_census import peak_files  # noqa: E402

RESULT = "inhibitor_census_result"
MULTRNA_RESULT = Path("data/results/multrna_result.json")
RESIDUEMULT = Path("data/results/residuemult_lookup.json")
RESIDUE_CENSUS = Path("data/results/residue_census.json")
MAGCAL = Path("data/results/magcal_result.json")
CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]
FORM = "mean"  # the form whose residue this is; the one that refused the 39 is saturate.fit
OWN_CODE = (
    "genomeos/attribution/inhibitor_census.py",
    "scripts/inhibitor_census_run.py",
    "tests/test_inhibitor_census.py",
)


def committed_symbols() -> set[str]:
    """The `symbol` column of `lane-multrna`'s committed per-gene table, AND NOTHING ELSE OF IT.

    The table's `P` and `D` columns are predictions and measured RNA. They are not returned, not kept and
    not looked at: this function reads one column out of the rows and discards the rest.
    """
    table = json.loads(MULTRNA_RESULT.read_text())["genes_table"]
    i = table["columns"].index("symbol")
    return {row[i] for row in table["rows"]}


def residuemult_neither() -> dict[str, Any]:
    """`lane-residuemult`'s committed answer, read for RECONCILIATION only (Q1), on `symbol` alone.

    The set of residue symbols that are in neither count, and the committed reason map. This lane's own
    count comes from its own run; this is what that count is checked against.
    """
    d = json.loads(RESIDUEMULT.read_text())
    t = d["per_gene_table"]
    i = {c: t["columns"].index(c) for c in t["columns"]}
    neither = {
        r[i["symbol"]] for r in t["rows"] if not (r[i["in_committed_table"]] or r[i["informative_scaled"]])
    }
    return {
        "source": str(RESIDUEMULT),
        "match_field": ic.MATCH_FIELD,
        "committed_residue_rows": len(t["rows"]),
        "committed_neither": sorted(neither),
        "committed_counts": d["counts"],
        "committed_reasons": d["informative_under_neither_by_reason"],
    }


def _provenance(module: Any, genes: set[str]) -> dict[tuple[str, str], tuple[str, str, str]]:
    """(gene, element) -> (evidence kind, how RHO was read, evidence source), from the gated module's
    own rules through `bridge.observation_of` -- the same call `combine2.gene_rules` reads RHO with."""
    out: dict[tuple[str, str], tuple[str, str, str]] = {}
    for r in module.rules:
        if r.action in c2.REGULATORY and r.target in genes:
            ob = bridge.observation_of(r, module)
            how = ob.how if isinstance(ob, bridge.Observation) else f"unresolved: {ob.reason}"
            out[(r.target, r.through)] = (bridge._kind(r), how, str(r.evidence.source))
    return out


def _chromosome(args: tuple[str, set[str], set[str]]) -> dict[str, Any]:
    """One program: the residue re-derived as `lane-residuemult` derived it, and every eligible gene's
    inhibiting elements censused. Nothing is fitted; `saturate.fit` has no parameter and SCALE is a
    constant."""
    from genomeos.lang.parser import parse

    chrom, genes, committed = args
    t0 = time.time()
    text = (Path(COMPILED_DIR) / f"noncoding_{chrom}.bio").read_text()
    gated, info = cp.gate_program(text, chrom, selection=c2.SELECTION)
    module = parse(gated, f"noncoding_{chrom}_reader")
    prov = _provenance(module, genes)
    out: dict[str, Any] = {
        "chrom": chrom,
        "program": {k: info[k] for k in ("base_sha256", "gated_sha256", "elements", "rules_gated")},
        "eligible": 0,
        "excluded": 0,
        "recovered": 0,
        "residue": 0,
        "neither": [],  # the 39, as this run finds them
        "neither_gene_rows": [],
        "neither_element_rows": [],
        "offending_gene_rows": [],  # every eligible gene with an offending element unscaled
        "strength_scaled": [],  # every inhibiting element of every eligible gene
        "strength_unscaled": [],
        "rho": [],
        "sign_disagreements": 0,
        "n_inhibiting_all": [],
        "n_rules_all": [],
        "n_inhibiting_residue": [],
        "n_rules_residue": [],
        "kinds_all": {},
        "hows_all": {},
        "violations": [],
    }
    kinds: Counter = Counter()
    hows: Counter = Counter()
    for gr in c2.gene_rules(module, genes, chrom, rg.compiled_cells(text)):
        _fits, why = c2.classify(gr)
        if why:
            continue
        out["eligible"] += 1
        before = c2.fit(FORM, gr.rho)
        excluded = (not before.ok) and before.reason == "basal_not_positive"
        in_residue = False
        if excluded:
            out["excluded"] += 1
            after = c2.fit(FORM, sf.scaled_rho(gr.rho, ic.SCALE))
            if after.q < 1.0:
                out["recovered"] += 1
            else:
                in_residue = True
                out["residue"] += 1
        kmap = {el: prov.get((gr.gene, el), ("", "", "")) for el in gr.elements}
        els = ic.element_rows(gr, kmap, before.q)
        row = ic.gene_row(gr, els, in_residue, excluded)
        row.elements = els
        out["violations"].extend(ic.gene_bound_violations(row))
        out["n_inhibiting_all"].append(row.n_inhibiting)
        out["n_rules_all"].append(row.n_rules)
        if in_residue:
            out["n_inhibiting_residue"].append(row.n_inhibiting)
            out["n_rules_residue"].append(row.n_rules)
        for e in els:
            out["strength_scaled"].append(e.strength_scaled)
            out["strength_unscaled"].append(e.strength)
            out["rho"].append(e.rho)
            kinds[e.kind] += 1
            hows[e.how] += 1
            # BOUND_CHECKS (e): the sign, re-asserted. classify() already enforces it.
            j = gr.elements.index(e.element)
            if gr.actions[j] != "inhibits":
                out["sign_disagreements"] += 1
                out["violations"].append(f"{gr.symbol}/{e.element}: RHO > 1 with action {gr.actions[j]} (e)")
        if row.n_offending_unscaled:
            out["offending_gene_rows"].append(row.as_row())
        if in_residue:
            lk = rm.lookup_gene(gr, gr.symbol in committed)
            if not (lk.in_committed_table or lk.informative_scaled):
                out["neither"].append(gr.symbol)
                out["neither_gene_rows"].append(row.as_row())
                out["neither_element_rows"].extend(e.as_row() for e in els)
                if lk.reason_scaled != row.reason_scaled:
                    out["violations"].append(
                        f"{gr.symbol}: residuemult's reason {lk.reason_scaled!r} against this run's"
                        f" {row.reason_scaled!r}"
                    )
    out["kinds_all"] = dict(kinds)
    out["hows_all"] = dict(hows)
    out["seconds"] = round(time.time() - t0, 1)
    return out


def census(chroms: list[str], workers: int) -> dict[str, Any]:
    """The whole census: re-derive, reconcile, then count (the committed QUESTIONS, in their order)."""
    texts = {c: (Path(COMPILED_DIR) / f"noncoding_{c}.bio").read_text() for c in chroms}
    chrom_of, _why = c2.multi_rule_genes(texts)
    del texts
    by_chrom: dict[str, set[str]] = {c: set() for c in chroms}
    for g, c in chrom_of.items():
        if c in by_chrom:
            by_chrom[c].add(g)
    committed = committed_symbols()
    ctx = get_context("fork")
    with ctx.Pool(workers) as pool:
        parts = pool.map(_chromosome, [(c, by_chrom[c], committed) for c in chroms])

    programs: dict[str, Any] = {}
    neither: list[str] = []
    neither_gene_rows: list[list[Any]] = []
    neither_element_rows: list[list[Any]] = []
    offending_gene_rows: list[list[Any]] = []
    s_scaled: list[float] = []
    s_unscaled: list[float] = []
    rhos: list[float] = []
    n_inh_all: list[int] = []
    n_rules_all: list[int] = []
    n_inh_res: list[int] = []
    n_rules_res: list[int] = []
    violations: list[str] = []
    kinds_all: Counter = Counter()
    hows_all: Counter = Counter()
    eligible = excluded = recovered = residue = signs = 0
    for p in parts:
        programs[p["chrom"]] = p["program"]
        neither.extend(p["neither"])
        neither_gene_rows.extend(p["neither_gene_rows"])
        neither_element_rows.extend(p["neither_element_rows"])
        offending_gene_rows.extend(p["offending_gene_rows"])
        s_scaled.extend(p["strength_scaled"])
        s_unscaled.extend(p["strength_unscaled"])
        rhos.extend(p["rho"])
        n_inh_all.extend(p["n_inhibiting_all"])
        n_rules_all.extend(p["n_rules_all"])
        n_inh_res.extend(p["n_inhibiting_residue"])
        n_rules_res.extend(p["n_rules_residue"])
        violations.extend(p["violations"])
        kinds_all.update(p["kinds_all"])
        hows_all.update(p["hows_all"])
        eligible += p["eligible"]
        excluded += p["excluded"]
        recovered += p["recovered"]
        residue += p["residue"]
        signs += p["sign_disagreements"]

    gi = {c: ic.GENE_COLUMNS.index(c) for c in ic.GENE_COLUMNS}
    ei = {c: ic.ELEMENT_COLUMNS.index(c) for c in ic.ELEMENT_COLUMNS}
    the_39 = set(neither)

    # ---- Q1: the re-derivation and its reconciliation --------------------------------------------------
    rmc = residuemult_neither()
    committed_neither = set(rmc["committed_neither"])
    by_reason = Counter(r[gi["reason_scaled"]] for r in neither_gene_rows)
    q1 = {
        "question": ic.QUESTIONS["Q1_re_derivation"],
        "match_field": ic.MATCH_FIELD,
        "eligible_for_a_fit": eligible,
        "excluded_basal_not_positive": excluded,
        "recovered_at_the_scale": recovered,
        "the_residue": residue,
        "informative_under_neither": len(neither_gene_rows),
        "by_reason_at_the_scale": dict(by_reason),
        "all_on_the_one_reason": set(by_reason) == {ic.REASON},
        "reconciliation": {
            "source": rmc["source"],
            "committed_neither_count": len(committed_neither),
            "re_derived_count": len(the_39),
            "sets_agree": the_39 == committed_neither,
            "in_this_run_not_committed": sorted(the_39 - committed_neither),
            "in_committed_not_this_run": sorted(committed_neither - the_39),
            "committed_counts": rmc["committed_counts"],
            "committed_reasons": rmc["committed_reasons"],
        },
        "symbols": sorted(the_39),
    }

    # ---- Q2: the distribution of inhibitor strength among the 39 ---------------------------------------
    q2 = {
        "question": ic.QUESTIONS["Q2_strength_distribution"],
        "definition": ic.IMPLIED_STRENGTH_DEFINITION,
        "per_gene_worst": {
            "max_strength_scaled": ic.distribution([r[gi["max_strength_scaled"]] for r in neither_gene_rows]),
            "max_strength_unscaled": ic.distribution([r[gi["max_strength"]] for r in neither_gene_rows]),
            "max_rho": ic.distribution([r[gi["max_rho"]] for r in neither_gene_rows]),
        },
        "per_element": {
            "strength_scaled": ic.distribution([r[ei["strength_scaled"]] for r in neither_element_rows]),
            "strength_unscaled": ic.distribution([r[ei["strength"]] for r in neither_element_rows]),
            "rho": ic.distribution([r[ei["rho"]] for r in neither_element_rows]),
            "mean_form_strength": ic.distribution(
                [r[ei["mean_form_strength"]] for r in neither_element_rows]
            ),
        },
        "offending_elements_only": {
            "strength_scaled": ic.distribution(
                [r[ei["strength_scaled"]] for r in neither_element_rows if r[ei["offends_scaled"]]]
            ),
            "rho": ic.distribution([r[ei["rho"]] for r in neither_element_rows if r[ei["offends_scaled"]]]),
        },
    }

    # ---- Q3: how many inhibitors each carries ----------------------------------------------------------
    q3 = {
        "question": ic.QUESTIONS["Q3_how_many_inhibitors"],
        "the_39": {
            "n_inhibiting": ic.distribution([float(r[gi["n_inhibiting"]]) for r in neither_gene_rows]),
            "n_inhibiting_counts": ic.counts_of([r[gi["n_inhibiting"]] for r in neither_gene_rows]),
            "n_rules": ic.distribution([float(r[gi["n_rules"]]) for r in neither_gene_rows]),
        },
        "the_residue": {
            "n_inhibiting": ic.distribution([float(x) for x in n_inh_res]),
            "n_rules": ic.distribution([float(x) for x in n_rules_res]),
        },
        "all_eligible": {
            "n_inhibiting": ic.distribution([float(x) for x in n_inh_all]),
            "n_rules": ic.distribution([float(x) for x in n_rules_all]),
            "genes_with_no_inhibitor": sum(1 for x in n_inh_all if x == 0),
        },
    }

    # ---- Q4: one inhibitor or several ------------------------------------------------------------------
    off = [r[gi["n_offending_scaled"]] for r in neither_gene_rows]
    exactly_one = [r for r in neither_gene_rows if r[gi["n_offending_scaled"]] == 1]
    q4 = {
        "question": ic.QUESTIONS["Q4_one_inhibitor_or_several"],
        "n_offending_at_the_scale": ic.distribution([float(x) for x in off]),
        "n_offending_at_the_scale_counts": ic.counts_of(off),
        "genes_with_exactly_one_offender": len(exactly_one),
        "genes_with_several_offenders": sum(1 for x in off if x > 1),
        "n_offending_unscaled_counts": ic.counts_of(
            [r[gi["n_offending_unscaled"]] for r in neither_gene_rows]
        ),
        "for_the_genes_with_exactly_one": {
            "share_of_the_gene_s_inhibitors_that_offends": ic.distribution(
                [1.0 / r[gi["n_inhibiting"]] for r in exactly_one if r[gi["n_inhibiting"]]]
            ),
            "n_inhibiting": ic.distribution([float(r[gi["n_inhibiting"]]) for r in exactly_one]),
        },
    }

    # ---- Q5: where the offending rules came from -------------------------------------------------------
    offenders = [r for r in neither_element_rows if r[ei["offends_scaled"]]]
    q5 = {
        "question": ic.QUESTIONS["Q5_where_the_rules_came_from"],
        "offending_elements": len(offenders),
        "by_evidence_kind": dict(Counter(r[ei["evidence_kind"]] for r in offenders)),
        "by_how_rho_was_read": dict(Counter(r[ei["how_rho_was_read"]] for r in offenders)),
        "by_evidence_source": dict(Counter(r[ei["evidence_source"]] for r in offenders)),
        "by_compiled_cell": dict(Counter(r[ei["compiled_cell"]] for r in offenders)),
        "by_chromosome": dict(Counter(r[ei["chrom"]] for r in offenders)),
        "background_all_inhibiting_elements_of_all_eligible_genes": {
            "elements": sum(kinds_all.values()),
            "by_evidence_kind": dict(kinds_all),
            "by_how_rho_was_read": dict(hows_all),
        },
    }

    # ---- Q6: how many genes outside the 39 carry an inhibitor above 1 ---------------------------------
    def outside(rows: list[list[Any]], key: str, pop: str | None) -> dict[str, int]:
        sel = [r for r in rows if r[gi[key]] > 0 and (pop is None or r[gi[pop]])]
        return {
            "genes": len(sel),
            "genes_outside_the_39": sum(1 for r in sel if r[gi["symbol"]] not in the_39),
            "elements": sum(r[gi[key]] for r in sel),
        }

    q6 = {
        "question": ic.QUESTIONS["Q6_outside_the_39"],
        "populations": {
            "all_eligible": eligible,
            "excluded_by_the_mean": excluded,
            "the_residue": residue,
            "the_39": len(the_39),
        },
        "above_the_bound_unscaled": {
            "all_eligible": outside(offending_gene_rows, "n_offending_unscaled", None),
            "excluded_by_the_mean": outside(
                offending_gene_rows, "n_offending_unscaled", "excluded_by_the_mean"
            ),
            "the_residue": outside(offending_gene_rows, "n_offending_unscaled", "in_residue"),
        },
        "above_the_bound_at_the_scale": {
            "all_eligible": outside(offending_gene_rows, "n_offending_scaled", None),
            "excluded_by_the_mean": outside(
                offending_gene_rows, "n_offending_scaled", "excluded_by_the_mean"
            ),
            "the_residue": outside(offending_gene_rows, "n_offending_scaled", "in_residue"),
        },
        "refused_by_saturate_unscaled_for_this_reason": sum(
            1 for r in offending_gene_rows if r[gi["reason_unscaled"]] == ic.REASON
        ),
    }

    # ---- Q7: is the threshold at a gap or in the body? ------------------------------------------------
    lo = 1.0 / ic.NEAR_BAND
    near_below = sum(1 for s in s_scaled if lo <= s <= 1.0)
    above = sum(1 for s in s_scaled if ic.above_bound(s))
    near_below_u = sum(1 for s in s_unscaled if lo <= s <= 1.0)
    above_u = sum(1 for s in s_unscaled if ic.above_bound(s))
    q7 = {
        "question": ic.QUESTIONS["Q7_threshold_or_tail"],
        "inhibiting_elements_of_all_eligible_genes": len(s_scaled),
        "strength_at_the_scale": ic.distribution(s_scaled),
        "strength_unscaled": ic.distribution(s_unscaled),
        "rho": ic.distribution(rhos),
        "near_band": ic.NEAR_BAND,
        "at_the_scale": {"just_below_the_bound": near_below, "above_the_bound": above},
        "unscaled": {"just_below_the_bound": near_below_u, "above_the_bound": above_u},
        "derived_thresholds": ic.derived_thresholds(),
    }

    which = ic.branch(near_below, above)
    out_u = q6["above_the_bound_unscaled"]["all_eligible"]["genes_outside_the_39"]
    return {
        "result": RESULT,
        "lane": ic.LANE,
        "question": ic.QUESTION,
        "questions": ic.QUESTIONS,
        "conventions": ic.CONVENTIONS,
        "derived_thresholds": ic.derived_thresholds(),
        "implied_strength_definition": ic.IMPLIED_STRENGTH_DEFINITION,
        "Q1_re_derivation": q1,
        "Q2_strength_distribution": q2,
        "Q3_how_many_inhibitors": q3,
        "Q4_one_inhibitor_or_several": q4,
        "Q5_where_the_rules_came_from": q5,
        "Q6_outside_the_39": q6,
        "Q7_threshold_or_tail": q7,
        "implication_rule": ic.IMPLICATION_RULE,
        "branch": which,
        "implication": ic.reading(
            which,
            near_below=near_below,
            near_band=ic.NEAR_BAND,
            above=above,
            outside_genes=out_u,
            residue=residue,
        ),
        "bound_checks": ic.BOUND_CHECKS,
        "bound_violations": violations,
        "sign_disagreements": signs,
        "what_this_is_not": ic.WHAT_THIS_IS_NOT,
        "limitations": list(ic.LIMITATIONS),
        "programs": programs,
        "the_39_gene_table": {"columns": list(ic.GENE_COLUMNS), "rows": neither_gene_rows},
        "the_39_element_table": {"columns": list(ic.ELEMENT_COLUMNS), "rows": neither_element_rows},
        "genes_with_an_offending_inhibitor_table": {
            "columns": list(ic.GENE_COLUMNS),
            "rows": offending_gene_rows,
        },
        "rna_read": "none",
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
    }


def manifest(chroms: list[str]) -> dict[str, Any]:
    inputs = [
        *(mf.input_entry(Path(COMPILED_DIR) / f"noncoding_{c}.bio", partition=None) for c in chroms),
        mf.input_entry(MULTRNA_RESULT, partition=None),
        mf.input_entry(RESIDUEMULT, partition=None),
        mf.input_entry(RESIDUE_CENSUS, partition=None),
        mf.input_entry(MAGCAL, partition=None),
        mf.input_entry(rg.GENE_TABLE, partition=None),
        mf.files_entry("reader v1 DNase narrowPeak files", [p for p in peak_files(chroms) if p.exists()]),
    ]
    return {
        "sources": [
            {
                "accession": "the 24 compiled non-coding programs 7d, 7h, lane-residue and"
                " lane-residuemult ran",
                "version": "pinned by sha256",
            },
            {"accession": "data/results/multrna_result.json (7h, the committed table)", "version": "dev"},
            {"accession": "data/results/residuemult_lookup.json (lane-residuemult)", "version": "d2065a1"},
            {"accession": "data/results/residue_census.json (lane-residue)", "version": "5d2c8f5"},
            {"accession": "data/results/magcal_result.json (M2h, the scale)", "version": "ac5580b"},
            {"accession": "reader v1 ENCODE DNase-seq narrowPeak files, 13 biosamples", "version": "pinned"},
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "form_of_the_residue": FORM,
            "form_that_refused_the_39": "multiplicative (saturate.fit)",
            "selection": c2.SELECTION,
            **ic.CONVENTIONS,
        },
        "exclusions": [
            "7d's own: one rule, rules in several programs, a measured rule, an unreadable or exactly-1"
            " observation, an action that disagrees with it, no activating rule",
            "an element with RHO < 1 is an activator and is not in any inhibitor distribution here",
            "a gene with no inhibiting element contributes to the gene counts and to no strength"
            " distribution",
        ],
        "partitions": "n/a: a predictor-side census; no outcome read and no RNA opened",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", nargs="*", default=CHROMS)
    ap.add_argument("--no-save", action="store_true")
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    t0 = time.time()
    out = census(args.chroms, args.workers)
    out["seconds"] = round(time.time() - t0, 1)
    q1 = out["Q1_re_derivation"]
    print(
        "eligible / excluded / recovered / residue:",
        q1["eligible_for_a_fit"],
        q1["excluded_basal_not_positive"],
        q1["recovered_at_the_scale"],
        q1["the_residue"],
        flush=True,
    )
    print(
        "informative under neither (re-derived):",
        q1["informative_under_neither"],
        q1["by_reason_at_the_scale"],
        flush=True,
    )
    print(
        "reconciliation:",
        {k: q1["reconciliation"][k] for k in ("committed_neither_count", "re_derived_count", "sets_agree")},
        flush=True,
    )
    print(
        "Q2 per-gene worst strength at the scale:",
        out["Q2_strength_distribution"]["per_gene_worst"]["max_strength_scaled"],
        flush=True,
    )
    print(
        "Q3 the 39's n_inhibiting:",
        out["Q3_how_many_inhibitors"]["the_39"]["n_inhibiting_counts"],
        flush=True,
    )
    print(
        "Q4 offenders per gene:",
        out["Q4_one_inhibitor_or_several"]["n_offending_at_the_scale_counts"],
        flush=True,
    )
    print("Q5 by evidence kind:", out["Q5_where_the_rules_came_from"]["by_evidence_kind"], flush=True)
    print("Q5 by how:", out["Q5_where_the_rules_came_from"]["by_how_rho_was_read"], flush=True)
    print("Q6 unscaled:", out["Q6_outside_the_39"]["above_the_bound_unscaled"], flush=True)
    print("Q6 at the scale:", out["Q6_outside_the_39"]["above_the_bound_at_the_scale"], flush=True)
    print(
        "Q7 near band:",
        out["Q7_threshold_or_tail"]["at_the_scale"],
        out["Q7_threshold_or_tail"]["unscaled"],
        flush=True,
    )
    print("violations:", len(out["bound_violations"]), out["bound_violations"][:5], flush=True)
    print("BRANCH", out["branch"], flush=True)
    print("IMPLICATION", out["implication"], flush=True)
    if not args.no_save:
        print("saved", save_result(RESULT, out, manifest=manifest(args.chroms), compact=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
