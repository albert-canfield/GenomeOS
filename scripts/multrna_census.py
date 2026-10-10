# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-multrna census: the multi-element genes the multiplicative composition runs, beside 7d's forms.

    uv run python scripts/multrna_census.py [--chroms chr21] [--no-save] [--no-run] [--workers 4]

Over the 24 generated programs 7d ran (data/knowledge/compiled, pinned by sha256), gated with the reader
selection exactly as 7d gated them (compile.gate_program, 'reader_cells', in memory):

1. every gene with two or more rules, classified by combine2.classify (7d's exclusions);
2. under 7d's two forms, runnable and informative by combine2.fit and combine2.predicted_log2 (the closed
   form 7d's census used), checked against data/results/combine2_census.json;
3. under the multiplicative form, runnable by saturate.fit and informative by its closed form; every
   runnable, informative gene is RUN by `runtime.grn.NetworkRuntime` (activator_combination=
   'multiplicative') in all thirteen reader cells with the HARNESS checks the registered run applies;
4. the genes the multiplicative form makes informative that neither 7d form could, and the genes all
   three forms make informative.

No RNA value is opened, no model request is sent, nothing is fetched.
"""

from __future__ import annotations

import argparse
import json
import math
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
from genomeos.attribution import multrna as mr  # noqa: E402
from genomeos.attribution import readergate as rg  # noqa: E402
from genomeos.attribution import saturate as sa  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from combine2_census import CHROMS, texts_of  # noqa: E402
from readergate_census import peak_files  # noqa: E402
from runrules_census import RATE_TABLE, load_rates  # noqa: E402

RESULT = "multrna_census"
C2_CENSUS = RESULTS_DIR / "combine2_census.json"
FORM = mr.FORM
ALL_FORMS = (*c2.FORMS, FORM)
OWN_CODE = (
    "genomeos/attribution/multrna.py",
    "scripts/multrna_census.py",
    "scripts/multrna_register.py",
    "scripts/multrna_run.py",
    "tests/test_multrna.py",
)
#: the 7d census counts this census must reproduce, per form
SEVEN_D_KEYS = ("runnable", "not_runnable", "informative", "informative_symbol_on_one_gene_id")
SEVEN_D_TOTALS = ("classified", "eligible_for_a_fit", "excluded_by_reason", "informative_under_both_forms")


def rate_of(gene: str, rates: dict, borrowed: float, prov: dict, lent: dict) -> tuple[float, dict[str, str]]:
    """The gene's rate T and its label, as scripts/combine2_census.py reads them."""
    if gene in rates:
        rate, label = rates[gene], {**prov.get(gene, {}), "tier": "measured"}
    else:
        rate, label = borrowed, {**lent, "tier": "borrowed_median"}
    return rate, {k: str(label.get(k, "unstated")) for k in ("source", "species_cell", "tier")}


def _chromosome(args: tuple[str, set[str], bool]) -> dict[str, Any]:
    """Gate, parse and classify one program's multi-element genes under all three forms; run the
    multiplicative form's runnable, informative ones."""
    from genomeos.lang.parser import parse

    chrom, genes, run = args
    t0 = time.time()
    rates, borrowed, prov, lent = load_rates()
    text = (Path(COMPILED_DIR) / f"noncoding_{chrom}.bio").read_text()
    gated, info = cp.gate_program(text, chrom, selection=mr.SELECTION)
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
            "forms": {},
        }
        if not why:
            for form, f in fits.items():
                fr: dict[str, Any] = {"ok": f.ok, "reason": f.reason}
                if f.ok:
                    s = c2.predicted_log2(f, gr)
                    fr["informative"] = max(s) - min(s) > c2.INFORMATIVE_LOG2
                row["forms"][form] = fr
            m = sa.fit(gr.rho)
            fr = {"ok": m.ok, "reason": m.reason}
            if m.ok:
                s = mr.predicted_log2(m, gr)
                fr["informative"] = max(s) - min(s) > mr.INFORMATIVE_LOG2
                fr["range_log2"] = max(s) - min(s)
                if fr["informative"] and run:
                    rate, label = rate_of(gr.gene, rates, borrowed, prov, lent)
                    got = mr.run_gene(rg.gated_gene_module(module, gr.gene), gr, m, rate, label)
                    fr.update(
                        {
                            "s": got["s"],
                            "closed_s": got["closed_s"],
                            "harness_miss": got["harness_miss"],
                            "tier": label["tier"],
                            "log2_rho": {el: math.log2(r) for el, r in gr.rho.items()},
                        }
                    )
            row["forms"][FORM] = fr
            row["sum_d"] = m.summed_drop
            row["sum_beta"] = m.summed_beta
            row["elements"] = gr.elements
            row["states"] = gr.states
            row["open_share"] = [x / len(gr.elements) for x in n_open]
        out["genes"].append(row)
    out["seconds"] = round(time.time() - t0, 1)
    return out


def collect(chroms: list[str], run: bool = True, workers: int = 4) -> tuple[list[dict], Counter, dict]:
    """Every multi-element gene of the programs, classified under all three forms, and (when `run`) run
    under the multiplicative form; 7d's exclusion counts; the gated programs' digests."""
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


def informative(g: dict, form: str) -> bool:
    fr = g["forms"].get(form, {})
    return bool(fr.get("ok") and fr.get("informative"))


def _why(g: dict, form: str) -> str:
    fr = g["forms"][form]
    return fr["reason"] if not fr["ok"] else "runnable_not_informative"


def summarise(genes: list[dict], why: Counter) -> dict[str, Any]:
    """The census counts (CENSUS); reads the gene table for the symbol check and no RNA."""
    ids = rg.symbol_gene_ids()
    elig = [g for g in genes if not g["excluded"]]
    forms: dict[str, Any] = {}
    for f in ALL_FORMS:
        ok = [g for g in elig if g["forms"][f]["ok"]]
        inf = [g for g in ok if g["forms"][f]["informative"]]
        forms[f] = {
            "runnable": len(ok),
            "not_runnable": dict(Counter(g["forms"][f]["reason"] for g in elig if not g["forms"][f]["ok"])),
            "informative": len(inf),
            "informative_symbol_on_one_gene_id": sum(1 for g in inf if g["symbol"] in ids),
        }
    mult = [g for g in elig if informative(g, FORM)]
    ran = [g for g in mult if "s" in g["forms"][FORM]]
    forms[FORM].update(
        {
            "informative_with_an_inhibiting_rule": sum(1 for g in mult if g["n_inhibiting"]),
            "informative_rule_count_quantiles": mr.quantiles([g["n_rules"] for g in mult]),
            "informative_summed_beta_quantiles": mr.quantiles([g["sum_beta"] for g in mult]),
            "predicted_range_log2_quantiles": mr.quantiles([g["forms"][FORM]["range_log2"] for g in mult]),
            "run": len(ran),
            "harness_failing_genes": sum(1 for g in ran if g["forms"][FORM]["harness_miss"]),
            "max_abs_run_minus_closed_form_log2": max(
                (
                    abs(a - b)
                    for g in ran
                    for a, b in zip(g["forms"][FORM]["s"], g["forms"][FORM]["closed_s"], strict=True)
                ),
                default=None,
            ),
        }
    )
    new = [g for g in mult if not any(informative(g, f) for f in c2.FORMS)]
    lost = [g for g in elig if any(informative(g, f) for f in c2.FORMS) and not informative(g, FORM)]
    common = [g for g in elig if all(informative(g, f) for f in ALL_FORMS)]
    return {
        **dict(why),
        "classified": len(genes),
        "excluded_by_reason": dict(Counter(g["excluded"] for g in genes if g["excluded"])),
        "eligible_for_a_fit": len(elig),
        "forms": forms,
        "informative_under_both_forms": sum(1 for g in elig if all(informative(g, f) for f in c2.FORMS)),
        "new_under_the_multiplicative_form": {
            "genes": len(new),
            "symbol_on_one_gene_id": sum(1 for g in new if g["symbol"] in ids),
            "why_not_under_the_mean": dict(sa._ranked(Counter(_why(g, "mean") for g in new))),
            "why_not_under_competition": dict(sa._ranked(Counter(_why(g, "competition") for g in new))),
            "rule_count_quantiles": mr.quantiles([g["n_rules"] for g in new]),
            "summed_drop_quantiles": mr.quantiles([g["sum_d"] for g in new]),
            "with_an_inhibiting_rule": sum(1 for g in new if g["n_inhibiting"]),
        },
        "informative_under_a_7d_form_but_not_the_multiplicative": {
            "genes": len(lost),
            "why": dict(sa._ranked(Counter(_why(g, FORM) for g in lost))),
        },
        "informative_under_the_multiplicative_and_the_mean": sum(1 for g in mult if informative(g, "mean")),
        "informative_under_the_multiplicative_and_competition": sum(
            1 for g in mult if informative(g, "competition")
        ),
        "informative_under_all_three": len(common),
        "informative_under_all_three_symbol_on_one_gene_id": sum(1 for g in common if g["symbol"] in ids),
        "rna_read": "none",
    }


def sentence(s: dict[str, Any]) -> str:
    """CENSUS_SENTENCE, filled from the counts only."""
    m = s["forms"][FORM]
    new = s["new_under_the_multiplicative_form"]
    return mr.CENSUS_SENTENCE.format(
        n=f"{s['eligible_for_a_fit']:,}",
        r=f"{m['runnable']:,}",
        x=f"{s['eligible_for_a_fit'] - m['runnable']:,}",
        reasons=sa._reasons(Counter(m["not_runnable"])),
        i=f"{m['informative']:,}",
        m=f"{s['forms']['mean']['informative']:,}",
        c=f"{s['forms']['competition']['informative']:,}",
        new=f"{new['genes']:,}",
        new_reasons=sa._reasons(Counter(new["why_not_under_the_mean"])),
        common=f"{s['informative_under_all_three']:,}",
    )


def seven_d_problems(s: dict[str, Any], committed: dict[str, Any]) -> list[str]:
    """Where this census's counts for 7d's forms differ from 7d's committed census."""
    rs = committed["reader_selection"]
    out = [f"{k}: {s[k]} against {rs[k]}" for k in SEVEN_D_TOTALS if s[k] != rs[k]]
    for f in c2.FORMS:
        out += [
            f"{f} {k}: {s['forms'][f][k]} against {rs['forms'][f][k]}"
            for k in SEVEN_D_KEYS
            if s["forms"][f][k] != rs["forms"][f][k]
        ]
    return out


def census(chroms: list[str], workers: int, run: bool) -> dict:
    t = time.time()
    genes, why, programs = collect(chroms, run=run, workers=workers)
    s = summarise(genes, why)
    out = {
        "lane": mr.LANE,
        "date": time.strftime("%Y-%m-%d"),
        "question": (
            "the multi-element genes the multiplicative composition runs, beside 7d's forms (no RNA read)"
        ),
        "census": mr.CENSUS,
        "form": FORM,
        "form_source": mr.FORM_SOURCE,
        "fit": mr.FIT,
        "structure": mr.STRUCTURE,
        "reader_selection": s,
        "sentence": sentence(s),
        "programs": programs,
        "ran_the_runtime": run,
        "seconds_collect": round(time.time() - t, 1),
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
    }
    if list(chroms) == CHROMS:
        problems = seven_d_problems(s, json.loads(C2_CENSUS.read_text()))
        out["seven_d_census_reproduced"] = {
            "source": "data/results/combine2_census.json",
            "problems": problems,
        }
    return out


def manifest(chroms: list[str]) -> dict:
    inputs = [
        *(mf.input_entry(Path(COMPILED_DIR) / f"noncoding_{c}.bio", partition=None) for c in chroms),
        mf.input_entry(RATE_TABLE, partition=None),
        mf.input_entry(rg.GENE_TABLE, partition=None),
        mf.files_entry("reader v1 DNase narrowPeak files", [p for p in peak_files(chroms) if p.exists()]),
        mf.input_entry(C2_CENSUS, partition="7d's committed census, read for the reproduction check"),
    ]
    return {
        "sources": [
            {"accession": "the 24 compiled non-coding programs 7d ran", "version": "pinned by sha256"},
            {"accession": "reader v1 ENCODE DNase-seq narrowPeak files, 13 biosamples", "version": "pinned"},
            {"accession": "data/results/gene_rates.json (mouse NIH 3T3, e921fae)", "version": "pinned"},
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "form": FORM,
            "comparator_forms": list(c2.FORMS),
            "theta": c2.THETA,
            "selection": mr.SELECTION,
            "runtime": mr.RUNTIME,
            "closed_form_tolerance": mr.CLOSED_FORM_TOLERANCE,
            "exponent_limit": mr.EXPONENT_LIMIT,
        },
        "exclusions": [
            "genes with one rule, rules in several programs, or a measured rule",
            "observations unreadable or exactly 1, actions that disagree with the observation",
            "genes with no activating rule; under the multiplicative form an inhibitor with RHO > 2 or summed"
            " log-folds beyond EXPONENT_LIMIT",
        ],
        "partitions": "n/a: a census; no outcome read",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", nargs="*", default=CHROMS)
    ap.add_argument("--no-save", action="store_true")
    ap.add_argument("--no-run", action="store_true", help="classify only; do not run the runtime")
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    t0 = time.time()
    out = census(args.chroms, args.workers, run=not args.no_run)
    out["seconds"] = round(time.time() - t0, 1)
    print({k: v for k, v in out["reader_selection"].items() if k != "forms"})
    print(json.dumps(out["reader_selection"]["forms"], indent=1))
    print(out["sentence"])
    if "seven_d_census_reproduced" in out:
        print("7d census problems:", out["seven_d_census_reproduced"]["problems"])
    if not args.no_save:
        print("saved", save_result(RESULT, out, manifest=manifest(args.chroms)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
