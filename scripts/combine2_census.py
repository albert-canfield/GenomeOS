# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-combine2 census: genes driven by several elements, and what each combination form can run.

    uv run python scripts/combine2_census.py [--chroms chr21] [--no-save] [--workers 4]

Over the 24 generated programs M1 and M2b ran (data/knowledge/compiled, pinned by sha256):

1. as compiled, no reader: per compiled `when` context, the gene-cell pairs with two or more mechanisms
   (the bridge's `not_identifiable`), and how many each form's joint fit makes runnable there, where all
   of a pair's observations come from one track (`combine2.compiled_context_census`);
2. under the reader selection (compile.gate_program, 'reader_cells', in memory, as M2b): every gene with
   two or more rules, why each is or is not runnable under each form, and per reader cell, ungated and
   gated, how many rules run;
3. every runnable, informative gene is RUN under each form by `runtime.grn.NetworkRuntime` in all thirteen
   reader cells with the HARNESS checks the registered run applies.

No RNA value is opened, no model request is sent, nothing is fetched.
"""

from __future__ import annotations

import argparse
import sys
import time
from collections import Counter
from multiprocessing import get_context
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import combine2 as c2  # noqa: E402
from genomeos.attribution import compile as cp  # noqa: E402
from genomeos.attribution import readergate as rg  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.results import save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readergate_census import peak_files  # noqa: E402
from runrules_census import RATE_TABLE, load_rates  # noqa: E402

RESULT = "combine2_census"
CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]
OWN_CODE = (
    "genomeos/runtime/grn.py",
    "genomeos/attribution/combine2.py",
    "scripts/combine2_census.py",
    "scripts/combine2_register.py",
    "scripts/combine2_run.py",
    "tests/test_combine2.py",
    "tests/test_grn_competition.py",
)


def texts_of(chroms: list[str]) -> dict[str, str]:
    return {c: (Path(COMPILED_DIR) / f"noncoding_{c}.bio").read_text() for c in chroms}


def _chromosome(args: tuple[str, set[str], bool]) -> dict[str, Any]:
    """Gate, parse and classify one program's multi-element genes; run the runnable, informative ones."""
    from genomeos.lang.parser import parse

    chrom, genes, run = args
    t0 = time.time()
    rates, borrowed, prov, lent = load_rates()
    text = (Path(COMPILED_DIR) / f"noncoding_{chrom}.bio").read_text()
    gated, info = cp.gate_program(text, chrom, selection=c2.SELECTION)
    module = parse(gated, f"noncoding_{chrom}_reader")
    out: dict[str, Any] = {
        "chrom": chrom,
        "program": {k: info[k] for k in ("base_sha256", "gated_sha256", "elements", "rules_gated")},
        "genes": [],
    }
    for gr in c2.gene_rules(module, genes, chrom, rg.compiled_cells(text)):
        fits, why = c2.classify(gr)
        n_open = [sum(1 for el in gr.elements if gr.states[el][j] == "open") for j in range(len(c2.CELLS))]
        row: dict[str, Any] = {
            "gene": gr.gene,
            "symbol": gr.symbol,
            "chrom": chrom,
            "n_rules": len(gr.elements),
            "n_inhibiting": sum(1 for r in gr.rho.values() if r > 1.0),
            "excluded": why,
            "n_open": n_open,
            "sum_d": sum(1.0 - r for r in gr.rho.values() if r < 1.0),
            "forms": {},
        }
        for form, f in fits.items():
            fr: dict[str, Any] = {"ok": f.ok, "reason": f.reason, "q": f.q}
            if f.ok:
                s = c2.predicted_log2(f, gr)
                fr["informative"] = max(s) - min(s) > c2.INFORMATIVE_LOG2
                fr["range_log2"] = max(s) - min(s)
                if fr["informative"] and run:
                    if gr.gene in rates:
                        rate, label = rates[gr.gene], {**prov.get(gr.gene, {}), "tier": "measured"}
                    else:
                        rate, label = borrowed, {**lent, "tier": "borrowed_median"}
                    label = {k: str(label.get(k, "unstated")) for k in ("source", "species_cell", "tier")}
                    got = c2.run_gene_form(rg.gated_gene_module(module, gr.gene), gr, f, rate, label)
                    fr.update(
                        {
                            "s": got["s"],
                            "harness_miss": got["harness_miss"],
                            "closed_s": got["closed_s"],
                            "tier": label["tier"],
                            "weights": f.weights,
                            "e": f.e,
                            "z": f.z,
                        }
                    )
            row["forms"][form] = fr
        if not why:
            row["elements"] = gr.elements
            row["states"] = gr.states
            row["open_share"] = [x / len(gr.elements) for x in n_open]
        out["genes"].append(row)
    out["seconds"] = round(time.time() - t0, 1)
    return out


def collect(chroms: list[str], run: bool = True, workers: int = 4) -> tuple[list[dict], Counter, dict]:
    """Every multi-element gene of the programs, classified, and (when `run`) run; exclusions; programs."""
    multi, why = c2.multi_rule_genes(texts_of(chroms))
    by_chrom: dict[str, set[str]] = {c: set() for c in chroms}
    for g, c in multi.items():
        by_chrom[c].add(g)
    order = sorted(chroms, key=lambda c: -(Path(COMPILED_DIR) / f"noncoding_{c}.bio").stat().st_size)
    jobs = [(c, by_chrom[c], run) for c in order]
    res: dict[str, dict] = {}
    if workers > 1:
        with get_context("spawn").Pool(workers) as pool:
            for r in pool.imap_unordered(_chromosome, jobs):
                res[r["chrom"]] = r
                print(r["chrom"], len(r["genes"]), r["seconds"], flush=True)
    else:
        for j in jobs:
            r = _chromosome(j)
            res[r["chrom"]] = r
            print(r["chrom"], len(r["genes"]), r["seconds"], flush=True)
    genes = [g for c in chroms for g in res[c]["genes"]]
    programs = {c: res[c]["program"] for c in chroms}
    return genes, why, programs


def per_cell(genes: list[dict]) -> dict[str, Any]:
    """Per reader cell: ungated (every offered rule runs) and gated (the reader holds), by form."""
    elig = [g for g in genes if not g["excluded"]]
    out = {}
    for j, cid in enumerate(c2.CELLS):
        cell: dict[str, Any] = {
            "ungated": {
                "genes_two_or_more_rules_offered": len(elig),
                **{f"runnable_{f}": sum(1 for g in elig if g["forms"][f]["ok"]) for f in c2.FORMS},
            },
            "gated": {
                "genes_by_open_rules": {
                    "0": sum(1 for g in elig if g["n_open"][j] == 0),
                    "1": sum(1 for g in elig if g["n_open"][j] == 1),
                    "2+": sum(1 for g in elig if g["n_open"][j] >= 2),
                }
            },
        }
        for f in c2.FORMS:
            ok = [g for g in elig if g["forms"][f]["ok"]]
            cell["gated"][f"runnable_{f}_by_open_rules"] = {
                "0": sum(1 for g in ok if g["n_open"][j] == 0),
                "1": sum(1 for g in ok if g["n_open"][j] == 1),
                "2+": sum(1 for g in ok if g["n_open"][j] >= 2),
            }
        out[cid] = cell
    return out


def _quantiles(xs: list[float]) -> dict[str, float]:
    s = sorted(xs)
    return {str(q): s[min(len(s) - 1, int(q * len(s)))] for q in (0.1, 0.25, 0.5, 0.75, 0.9)} if s else {}


def summarise(genes: list[dict], why: Counter) -> dict[str, Any]:
    ids = rg.symbol_gene_ids()
    elig = [g for g in genes if not g["excluded"]]
    forms: dict[str, Any] = {}
    for f in c2.FORMS:
        ok = [g for g in elig if g["forms"][f]["ok"]]
        inf = [g for g in ok if g["forms"][f]["informative"]]
        ran = [g for g in inf if "s" in g["forms"][f]]
        forms[f] = {
            "runnable": len(ok),
            "not_runnable": dict(Counter(g["forms"][f]["reason"] for g in elig if not g["forms"][f]["ok"])),
            "informative": len(inf),
            "informative_symbol_on_one_gene_id": sum(1 for g in inf if g["symbol"] in ids),
            "informative_with_an_inhibiting_rule": sum(1 for g in inf if g["n_inhibiting"]),
            "informative_rule_count_quantiles": _quantiles([g["n_rules"] for g in inf]),
            "predicted_range_log2_quantiles": _quantiles([g["forms"][f]["range_log2"] for g in inf]),
            "run": len(ran),
            "harness_failing_genes": sum(1 for g in ran if g["forms"][f]["harness_miss"]),
            "max_abs_run_minus_closed_form_log2": max(
                (
                    abs(a - b)
                    for g in ran
                    for a, b in zip(g["forms"][f]["s"], g["forms"][f]["closed_s"], strict=True)
                ),
                default=None,
            ),
        }
    both = [g for g in elig if all(g["forms"][f]["ok"] and g["forms"][f]["informative"] for f in c2.FORMS)]
    return {
        **dict(why),
        "classified": len(genes),
        "excluded_by_reason": dict(Counter(g["excluded"] for g in genes if g["excluded"])),
        "eligible_for_a_fit": len(elig),
        "summed_single_removal_drop_quantiles": _quantiles([g["sum_d"] for g in elig]),
        "forms": forms,
        "informative_under_both_forms": len(both),
        "informative_under_both_symbol_on_one_gene_id": sum(1 for g in both if g["symbol"] in ids),
        "rna_read": "none",
    }


def census(chroms: list[str], workers: int) -> dict:
    t = time.time()
    compiled = c2.compiled_context_census(texts_of(chroms))
    print("as compiled", compiled["totals"], round(time.time() - t, 1), flush=True)
    genes, why, programs = collect(chroms, run=True, workers=workers)
    return {
        "lane": c2.LANE,
        "date": time.strftime("%Y-%m-%d"),
        "question": "genes driven by several elements, and what each combination form can run (no RNA read)",
        "forms": list(c2.FORMS),
        "form_sources": c2.FORM_SOURCES,
        "theta": c2.THETA,
        "fit": c2.FIT,
        "structure": c2.STRUCTURE,
        "as_compiled_no_reader": compiled,
        "reader_selection": summarise(genes, why),
        "per_reader_cell": per_cell(genes),
        "programs": programs,
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
    }


def manifest(chroms: list[str]) -> dict:
    inputs = [
        *(mf.input_entry(Path(COMPILED_DIR) / f"noncoding_{c}.bio", partition=None) for c in chroms),
        mf.input_entry(RATE_TABLE, partition=None),
        mf.input_entry(rg.GENE_TABLE, partition=None),
        mf.files_entry("reader v1 DNase narrowPeak files", [p for p in peak_files(chroms) if p.exists()]),
    ]
    return {
        "sources": [
            {
                "accession": "the 24 compiled non-coding programs M1 and M2b ran",
                "version": "pinned by sha256",
            },
            {"accession": "reader v1 ENCODE DNase-seq narrowPeak files, 13 biosamples", "version": "pinned"},
            {"accession": "data/results/gene_rates.json (mouse NIH 3T3, e921fae)", "version": "pinned"},
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "forms": list(c2.FORMS),
            "theta": c2.THETA,
            "selection": c2.SELECTION,
            "runtime": c2.RUNTIME,
            "closed_form_tolerance": c2.CLOSED_FORM_TOLERANCE,
        },
        "exclusions": [
            "genes with one rule, rules in several programs, or a measured rule",
            "observations unreadable or exactly 1, actions that disagree with the observation",
            "genes with no activating rule; fits that leave no positive basal rate or an inhibitor s > 1",
        ],
        "partitions": "n/a: a census; no outcome read",
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
    print({k: v for k, v in out["reader_selection"].items() if k != "forms"})
    print(out["reader_selection"]["forms"])
    if not args.no_save:
        print("saved", save_result(RESULT, out, manifest=manifest(args.chroms)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
