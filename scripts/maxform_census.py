# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-maxform census: how many genes `max` and `or` separate from the runtime mean, before any RNA.

    uv run python scripts/maxform_census.py [--chroms chr21] [--no-save] [--no-run] [--workers 4]

Over the 24 generated programs 7d ran (data/knowledge/compiled, pinned by sha256), gated with the reader
selection exactly as 7d gated them (compile.gate_program, 'reader_cells', in memory):

1. every gene with two or more rules, classified by combine2.classify (7d's exclusions);
2. under 7d's two forms and 7h's multiplicative form, runnable and informative by their closed forms,
   checked against data/results/combine2_census.json and data/results/multrna_census.json;
3. under `max` and `or` (grn.ACTIVATOR_RULES, the engine unchanged), maxform.fit, and for every gene a
   form fits the largest difference over the thirteen reader cells between its closed-form prediction and
   the runtime mean's -- THE CENSUS QUESTION: how many genes does the form's prediction differ on;
4. every informative gene RUN by `runtime.grn.NetworkRuntime` (activator_combination=form) in all thirteen
   reader cells with combine2's HARNESS checks;
5. maxform.max_identity_checks, the algebra of maxform.PREMISE measured on random observations;
6. maxform.DECISION_BEFORE_THE_CENSUS, per form, from those counts.

No RNA value is opened, no model request is sent, nothing is fetched.
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
from genomeos.attribution import combine2 as c2  # noqa: E402
from genomeos.attribution import compile as cp  # noqa: E402
from genomeos.attribution import maxform as xf  # noqa: E402
from genomeos.attribution import multrna as mr  # noqa: E402
from genomeos.attribution import readergate as rg  # noqa: E402
from genomeos.attribution import saturate as sa  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from combine2_census import CHROMS, texts_of  # noqa: E402
from readergate_census import peak_files  # noqa: E402
from runrules_census import RATE_TABLE, load_rates  # noqa: E402

RESULT = "maxform_census"
C2_CENSUS = RESULTS_DIR / "combine2_census.json"
MR_CENSUS = RESULTS_DIR / "multrna_census.json"
MULT = mr.FORM
OWN_CODE = (
    "genomeos/attribution/maxform.py",
    "scripts/maxform_census.py",
    "scripts/maxform_register.py",
    "scripts/maxform_run.py",
    "tests/test_maxform.py",
)


def rate_of(gene: str, rates: dict, borrowed: float, prov: dict, lent: dict) -> tuple[float, dict[str, str]]:
    """The gene's rate T and its label, as scripts/combine2_census.py reads them."""
    if gene in rates:
        rate, label = rates[gene], {**prov.get(gene, {}), "tier": "measured"}
    else:
        rate, label = borrowed, {**lent, "tier": "borrowed_median"}
    return rate, {k: str(label.get(k, "unstated")) for k in ("source", "species_cell", "tier")}


def _informative(s: list[float]) -> bool:
    return max(s) - min(s) > xf.INFORMATIVE_LOG2


def _chromosome(args: tuple[str, set[str], bool]) -> dict[str, Any]:
    """Gate, parse and classify one program's multi-element genes under every form; run this lane's forms
    on their informative genes."""
    from genomeos.lang.parser import parse

    chrom, genes, run = args
    t0 = time.time()
    rates, borrowed, prov, lent = load_rates()
    text = (Path(COMPILED_DIR) / f"noncoding_{chrom}.bio").read_text()
    gated, info = cp.gate_program(text, chrom, selection=xf.SELECTION)
    module = parse(gated, f"noncoding_{chrom}_reader")
    out: dict[str, Any] = {
        "chrom": chrom,
        "program": {k: info[k] for k in ("base_sha256", "gated_sha256", "elements", "rules_gated")},
        "genes": [],
    }
    for gr in c2.gene_rules(module, genes, chrom, rg.compiled_cells(text)):
        fits, why = c2.classify(gr)
        d = [1.0 - r for r in gr.rho.values() if r < 1.0]
        row: dict[str, Any] = {
            "gene": gr.gene,
            "symbol": gr.symbol,
            "chrom": chrom,
            "n_rules": len(gr.elements),
            "n_activating": len(d),
            "n_inhibiting": sum(1 for r in gr.rho.values() if r > 1.0),
            "excluded": why,
            "forms": {},
        }
        if not why:
            for form, f in fits.items():
                fr: dict[str, Any] = {"ok": f.ok, "reason": f.reason}
                if f.ok:
                    fr["informative"] = _informative(c2.predicted_log2(f, gr))
                row["forms"][form] = fr
            m = sa.fit(gr.rho)
            row["forms"][MULT] = {"ok": m.ok, "reason": m.reason}
            if m.ok:
                row["forms"][MULT]["informative"] = _informative(mr.predicted_log2(m, gr))
            for form in xf.FORMS:
                mfit = xf.fit(form, gr.rho)
                fr = {"ok": mfit.ok, "reason": mfit.reason, "q": mfit.q, "z": mfit.z}
                if mfit.ok:
                    s = xf.predicted_log2(mfit, gr)
                    fr["informative"] = _informative(s)
                    fr["range_log2"] = max(s) - min(s)
                    fr["difference_from_mean_log2"] = xf.difference_from_mean(mfit, gr)
                    fr["differs_from_mean"] = xf.differs_from_mean(mfit, gr)
                    fr["mean_has_a_fit"] = fr["difference_from_mean_log2"] is not None
                    if fr["informative"] and run:
                        rate, label = rate_of(gr.gene, rates, borrowed, prov, lent)
                        got = xf.run_gene(rg.gated_gene_module(module, gr.gene), gr, mfit, rate, label)
                        fr["harness_miss"] = got["harness_miss"]
                        fr["max_abs_run_minus_closed_form_log2"] = max(
                            abs(a - b) for a, b in zip(got["s"], got["closed_s"], strict=True)
                        )
                        fr["tier"] = label["tier"]
                row["forms"][form] = fr
            row["sum_d"] = sum(d)
            row["top_d"] = max(d) if d else None
        out["genes"].append(row)
    out["seconds"] = round(time.time() - t0, 1)
    return out


def collect(chroms: list[str], run: bool = True, workers: int = 4) -> tuple[list[dict], Counter, dict]:
    """Every multi-element gene of the programs, classified under every form, and (when `run`) this lane's
    forms' informative genes run; 7d's exclusion counts; the gated programs' digests."""
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


def _max(xs: list[Any]) -> float | None:
    vals = [x for x in xs if x is not None]
    return max(vals, default=None)


def summarise(genes: list[dict], why: Counter) -> dict[str, Any]:
    """The census counts (maxform.CENSUS); reads the gene table for the symbol check and no RNA."""
    ids = rg.symbol_gene_ids()
    elig = [g for g in genes if not g["excluded"]]
    forms: dict[str, Any] = {}
    for f in (*c2.FORMS, MULT):
        ok = [g for g in elig if g["forms"][f]["ok"]]
        inf = [g for g in ok if g["forms"][f]["informative"]]
        forms[f] = {
            "runnable": len(ok),
            "not_runnable": dict(Counter(g["forms"][f]["reason"] for g in elig if not g["forms"][f]["ok"])),
            "informative": len(inf),
            "informative_symbol_on_one_gene_id": sum(1 for g in inf if g["symbol"] in ids),
        }
    own: dict[str, Any] = {}
    for f in xf.FORMS:
        ok = [g for g in elig if g["forms"][f]["ok"]]
        inf = [g for g in ok if g["forms"][f]["informative"]]
        ran = [g for g in inf if "harness_miss" in g["forms"][f]]
        diff = [g for g in ok if g["forms"][f]["differs_from_mean"]]
        dinf = [g for g in diff if g["forms"][f]["informative"]]
        own[f] = {
            "fitted": len(ok),
            "not_fitted": dict(Counter(g["forms"][f]["reason"] for g in elig if not g["forms"][f]["ok"])),
            "informative": len(inf),
            "informative_symbol_on_one_gene_id": sum(1 for g in inf if g["symbol"] in ids),
            "run": len(ran),
            "harness_failing_genes": sum(1 for g in ran if g["forms"][f]["harness_miss"]),
            "max_abs_run_minus_closed_form_log2": _max(
                [g["forms"][f]["max_abs_run_minus_closed_form_log2"] for g in ran]
            ),
            # THE CENSUS QUESTION
            "differs_from_the_mean": {
                "genes": len(diff),
                "informative": len(dinf),
                "informative_symbol_on_one_gene_id": sum(1 for g in dinf if g["symbol"] in ids),
                "the_mean_has_no_fit": sum(1 for g in diff if not g["forms"][f]["mean_has_a_fit"]),
                "max_abs_difference_from_mean_log2": _max(
                    [g["forms"][f]["difference_from_mean_log2"] for g in diff]
                ),
                "difference_from_mean_log2_quantiles": mr.quantiles(
                    [x for g in diff if (x := g["forms"][f]["difference_from_mean_log2"]) is not None]
                )
                or None,
            },
            "max_abs_difference_from_mean_log2_over_all_fitted": _max(
                [g["forms"][f]["difference_from_mean_log2"] for g in ok]
            ),
            "activating_count_quantiles": mr.quantiles([g["n_activating"] for g in ok]) if ok else None,
            "q_quantiles": mr.quantiles([g["forms"][f]["q"] for g in ok]) if ok else None,
            "informative_and_mean_informative": sum(1 for g in inf if informative(g, "mean")),
            "informative_and_multiplicative_informative": sum(1 for g in inf if informative(g, MULT)),
            "informative_under_this_form_and_all_three_of_m2h": sum(
                1 for g in inf if all(informative(g, x) for x in (*c2.FORMS, MULT))
            ),
        }
    return {
        **dict(why),
        "classified": len(genes),
        "excluded_by_reason": dict(Counter(g["excluded"] for g in genes if g["excluded"])),
        "eligible_for_a_fit": len(elig),
        "forms": forms,
        "this_lanes_forms": own,
        "rna_read": "none",
    }


def sentence(s: dict[str, Any]) -> str:
    """CENSUS_SENTENCE, filled from the counts only."""
    mx, orr = s["this_lanes_forms"]["max"], s["this_lanes_forms"]["or"]
    dmax = orr["differs_from_the_mean"]["max_abs_difference_from_mean_log2"]
    return xf.CENSUS_SENTENCE.format(
        n=f"{s['eligible_for_a_fit']:,}",
        max_fitted=f"{mx['fitted']:,}",
        max_differ=f"{mx['differs_from_the_mean']['genes']:,}",
        max_reasons=sa._reasons(Counter(mx["not_fitted"])),
        or_fitted=f"{orr['fitted']:,}",
        or_differ=f"{orr['differs_from_the_mean']['genes']:,}",
        or_dmax="none" if dmax is None else f"{dmax:.4g}",
        or_reasons=sa._reasons(Counter(orr["not_fitted"])),
        max_pop=f"{mx['differs_from_the_mean']['informative_symbol_on_one_gene_id']:,}",
        or_pop=f"{orr['differs_from_the_mean']['informative_symbol_on_one_gene_id']:,}",
        min=xf.MIN_GENES,
        ran=f"{mx['run'] + orr['run']:,}",
        miss=f"{mx['harness_failing_genes'] + orr['harness_failing_genes']:,}",
    )


SEVEN_D_KEYS = ("runnable", "informative", "informative_symbol_on_one_gene_id")


def committed_problems(s: dict[str, Any], c2c: dict, mrc: dict) -> list[str]:
    """7d's and 7h's committed per-form counts, recomputed here and compared (HARNESS)."""
    bad = []
    want_c2 = c2c["reader_selection"]["forms"]
    for f in c2.FORMS:
        for k in SEVEN_D_KEYS:
            got, exp = s["forms"][f][k], want_c2[f].get(k)
            if exp is not None and got != exp:
                bad.append(f"7d {f}.{k}: recomputed {got}, committed {exp}")
    if s["eligible_for_a_fit"] != c2c["reader_selection"].get("eligible_for_a_fit"):
        bad.append(
            f"eligible_for_a_fit: recomputed {s['eligible_for_a_fit']}, committed"
            f" {c2c['reader_selection'].get('eligible_for_a_fit')}"
        )
    want_mr = mrc.get("reader_selection", {}).get("forms", {}).get(MULT, {})
    for k in ("runnable", "informative"):
        got, exp = s["forms"][MULT][k], want_mr.get(k)
        if exp is not None and got != exp:
            bad.append(f"7h {MULT}.{k}: recomputed {got}, committed {exp}")
    return bad


def census(chroms: list[str], workers: int, run: bool = True) -> dict[str, Any]:
    t = time.time()
    genes, why, programs = collect(chroms, run=run, workers=workers)
    s = summarise(genes, why)
    out = {
        "lane": xf.LANE,
        "date": time.strftime("%Y-%m-%d"),
        "roadmap_item": xf.ROADMAP_ITEM,
        "question": (
            "on 7d's population, how many genes does the prediction of grn's `max`, and of grn's `or`,"
            " differ from the runtime mean's on, before any RNA value is read"
        ),
        "census": xf.CENSUS,
        "forms": list(xf.FORMS),
        "form_sources": xf.FORM_SOURCES,
        "conventions": xf.CONVENTIONS,
        "premise": xf.PREMISE,
        "fit": xf.FIT,
        "structure": xf.STRUCTURE,
        "decision_rule": xf.DECISION_BEFORE_THE_CENSUS,
        "root_convention_cited": xf.FIT["the_or_root_is_unique"],
        "exploration_preceded_this": xf.EXPLORATION_PRECEDED_THIS,
        "reader_selection": s,
        "sentence": sentence(s),
        "decision": {f: xf.decide_from_census(s["this_lanes_forms"], f) for f in xf.FORMS},
        "identity_checks": xf.max_identity_checks(),
        "identity_checks_note": (
            "the largest violation of each claim of PREMISE over 20,000 random observation vectors;"
            " max_fits_two_activators and max_differs_from_mean_at_one are exactly 0, or_removal_error is"
            " floating-point noise on an exact closed form"
        ),
        "programs": programs,
        "ran_the_runtime": run,
        "seconds_collect": round(time.time() - t, 1),
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
    }
    if list(chroms) == CHROMS:
        out["committed_censuses_reproduced"] = {
            "sources": ["data/results/combine2_census.json", "data/results/multrna_census.json"],
            "problems": committed_problems(
                s, json.loads(C2_CENSUS.read_text()), json.loads(MR_CENSUS.read_text())
            ),
        }
    return out


def manifest(chroms: list[str]) -> dict:
    inputs = [
        *(mf.input_entry(Path(COMPILED_DIR) / f"noncoding_{c}.bio", partition=None) for c in chroms),
        mf.input_entry(RATE_TABLE, partition=None),
        mf.input_entry(rg.GENE_TABLE, partition=None),
        mf.files_entry("reader v1 DNase narrowPeak files", [p for p in peak_files(chroms) if p.exists()]),
        mf.input_entry(C2_CENSUS, partition="7d's committed census, read for the reproduction check"),
        mf.input_entry(MR_CENSUS, partition="7h's committed census, read for the reproduction check"),
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
            "forms": list(xf.FORMS),
            "comparator_forms": [*c2.FORMS, MULT],
            "h_intact": xf.H,
            "selection": xf.SELECTION,
            "runtime": xf.RUNTIME,
            "closed_form_tolerance": xf.CLOSED_FORM_TOLERANCE,
            "informative_log2": xf.INFORMATIVE_LOG2,
        },
        "exclusions": [
            "genes with one rule, rules in several programs, or a measured rule",
            "observations unreadable or exactly 1, actions that disagree with the observation",
            "genes with no activating rule; under `max`, every gene with two or more activating elements"
            " (max_needs_one_activating_element, maxform.PREMISE); under `or`, basal_not_positive and"
            " inhibitor_strength_above_1",
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
    s = out["reader_selection"]
    print(json.dumps({"eligible_for_a_fit": s["eligible_for_a_fit"], **s["this_lanes_forms"]}, indent=1))
    print(json.dumps(s["forms"], indent=1))
    print(out["sentence"])
    print("decision:", out["decision"])
    print("identity checks:", out["identity_checks"])
    if "committed_censuses_reproduced" in out:
        print("committed census problems:", out["committed_censuses_reproduced"]["problems"])
    if not args.no_save:
        print("saved", save_result(RESULT, out, manifest=manifest(args.chroms)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
