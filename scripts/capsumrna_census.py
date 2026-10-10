# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-capsumrna census: the multi-element genes the capped sum runs on 7d's design, beside 7d's and 7h's.

    uv run python scripts/capsumrna_census.py [--chroms chr21] [--no-save] [--no-run] [--workers 4]

Over the 24 generated programs 7d ran (data/knowledge/compiled, pinned by sha256), gated with the reader
selection exactly as 7d gated them (compile.gate_program, 'reader_cells', in memory):

1. every gene with two or more rules, classified by combine2.classify (7d's exclusions);
2. under 7d's two forms and 7h's multiplicative form, runnable and informative by their closed forms,
   checked against data/results/combine2_census.json and data/results/multrna_census.json;
3. under the capped sum (grn 'sum_capped'), both roots of capsumrna.FIT by capsumrna.classify, each
   root's closed-form prediction against the mean's, and every informative gene on each root RUN by
   `runtime.grn.NetworkRuntime` (activator_combination='sum_capped') in all thirteen reader cells with
   combine2's HARNESS checks;
4. capsumrna.DECISION, from those counts.

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
from genomeos.attribution import capsumrna as cs  # noqa: E402
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

RESULT = "capsumrna_census"
C2_CENSUS = RESULTS_DIR / "combine2_census.json"
MR_CENSUS = RESULTS_DIR / "multrna_census.json"
MULT = mr.FORM
ROOTS = ("uncapped", "capped")
OWN_CODE = ("genomeos/attribution/capsumrna.py", "scripts/capsumrna_census.py", "tests/test_capsumrna.py")
SEVEN_D_KEYS = ("runnable", "not_runnable", "informative", "informative_symbol_on_one_gene_id")
SEVEN_D_TOTALS = ("classified", "eligible_for_a_fit", "excluded_by_reason")
SEVEN_H_KEYS = ("runnable", "informative")


def rate_of(gene: str, rates: dict, borrowed: float, prov: dict, lent: dict) -> tuple[float, dict[str, str]]:
    """The gene's rate T and its label, as scripts/combine2_census.py reads them."""
    if gene in rates:
        rate, label = rates[gene], {**prov.get(gene, {}), "tier": "measured"}
    else:
        rate, label = borrowed, {**lent, "tier": "borrowed_median"}
    return rate, {k: str(label.get(k, "unstated")) for k in ("source", "species_cell", "tier")}


def _informative(s: list[float]) -> bool:
    return max(s) - min(s) > cs.INFORMATIVE_LOG2


def _chromosome(args: tuple[str, set[str], bool]) -> dict[str, Any]:
    """Gate, parse and classify one program's multi-element genes under every form; run the capped sum's
    informative genes on each root."""
    from genomeos.lang.parser import parse

    chrom, genes, run = args
    t0 = time.time()
    rates, borrowed, prov, lent = load_rates()
    text = (Path(COMPILED_DIR) / f"noncoding_{chrom}.bio").read_text()
    gated, info = cp.gate_program(text, chrom, selection=cs.SELECTION)
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
            cl = cs.classify(gr)
            row["status"] = cl["status"]
            row["roots"] = {}
            for root in ROOTS:
                f = cl[root]
                fr = {"ok": f.ok, "reason": f.reason, "w": f.w}
                if f.ok:
                    s = cs.predicted_log2(f, gr)
                    fr["informative"] = _informative(s)
                    fr["range_log2"] = max(s) - min(s)
                    fr["difference_from_mean_log2"] = cs.difference_from_mean(f, gr)
                    if fr["informative"] and run:
                        rate, label = rate_of(gr.gene, rates, borrowed, prov, lent)
                        got = cs.run_gene(rg.gated_gene_module(module, gr.gene), gr, f, rate, label)
                        fr["harness_miss"] = got["harness_miss"]
                        fr["max_abs_run_minus_closed_form_log2"] = max(
                            abs(a - b) for a, b in zip(got["s"], got["closed_s"], strict=True)
                        )
                        fr["tier"] = label["tier"]
                row["roots"][root] = fr
            row["sum_d"] = sum(d)
            row["top_d"] = max(d)
        out["genes"].append(row)
    out["seconds"] = round(time.time() - t0, 1)
    return out


def collect(chroms: list[str], run: bool = True, workers: int = 4) -> tuple[list[dict], Counter, dict]:
    """Every multi-element gene of the programs, classified under every form, and (when `run`) the capped
    sum's informative genes run on each root; 7d's exclusion counts; the gated programs' digests."""
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


def _root_inf(g: dict, root: str) -> bool:
    fr = g["roots"][root]
    return bool(fr["ok"] and fr["informative"])


def _max(xs: list[float]) -> float | None:
    return max(xs, default=None)


def summarise(genes: list[dict], why: Counter) -> dict[str, Any]:
    """The census counts (capsumrna.CENSUS); reads the gene table for the symbol check and no RNA."""
    if any(g.get("status") == "capped_root_only" for g in genes):  # after 43dd4ac: a_capped_root_alone
        return summarise_capped_root_alone(genes, why)
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
    roots: dict[str, Any] = {}
    for root in ROOTS:
        ok = [g for g in elig if g["roots"][root]["ok"]]
        inf = [g for g in ok if g["roots"][root]["informative"]]
        ran = [g for g in inf if "harness_miss" in g["roots"][root]]
        roots[root] = {
            "fitted": len(ok),
            "informative": len(inf),
            "informative_and_mean_informative": sum(1 for g in inf if informative(g, "mean")),
            "run": len(ran),
            "harness_failing_genes": sum(1 for g in ran if g["roots"][root]["harness_miss"]),
            "max_abs_run_minus_closed_form_log2": _max(
                [g["roots"][root]["max_abs_run_minus_closed_form_log2"] for g in ran]
            ),
            "max_abs_difference_from_mean_log2": _max(
                [g["roots"][root]["difference_from_mean_log2"] for g in ok]
            ),
            "difference_from_mean_log2_quantiles": mr.quantiles(
                [g["roots"][root]["difference_from_mean_log2"] for g in ok]
            )
            if ok
            else None,
            "activating_count_quantiles": mr.quantiles([g["n_activating"] for g in ok]) if ok else None,
        }
    roots["uncapped"]["not_fitted"] = dict(
        Counter(g["roots"]["uncapped"]["reason"] for g in elig if not g["roots"]["uncapped"]["ok"])
    )
    unique = [g for g in elig if g["status"] == "runnable"]
    uinf = [g for g in unique if _root_inf(g, "uncapped")]
    common = [g for g in uinf if all(informative(g, f) for f in (*c2.FORMS, MULT))]
    two = [g for g in elig if g["status"] == "fit_not_unique"]
    return {
        **dict(why),
        "classified": len(genes),
        "excluded_by_reason": dict(Counter(g["excluded"] for g in genes if g["excluded"])),
        "eligible_for_a_fit": len(elig),
        "forms": forms,
        "capped_sum_status": dict(sa._ranked(Counter(g["status"] for g in elig))),
        "roots": roots,
        "unique_fit": {
            "genes": len(unique),
            "informative": len(uinf),
            "informative_symbol_on_one_gene_id": sum(1 for g in uinf if g["symbol"] in ids),
            "max_abs_difference_from_mean_log2": _max(
                [g["roots"]["uncapped"]["difference_from_mean_log2"] for g in unique]
            ),
            "informative_under_all_four": len(common),
            "informative_under_all_four_symbol_on_one_gene_id": sum(1 for g in common if g["symbol"] in ids),
            "share_of_the_means_informative_genes": round(
                len(uinf) / max(1, forms["mean"]["informative"]), 4
            ),
        },
        "two_roots": {
            "genes": len(two),
            "informative_on_the_capped_root": sum(1 for g in two if _root_inf(g, "capped")),
            "informative_on_the_capped_root_symbol_on_one_gene_id": sum(
                1 for g in two if _root_inf(g, "capped") and g["symbol"] in ids
            ),
            "summed_drop_quantiles": mr.quantiles([g["sum_d"] for g in two]) if two else None,
        },
        "rna_read": "none",
    }


def sentence(s: dict[str, Any]) -> str:
    """CENSUS_SENTENCE, filled from the counts only."""
    if "capped_root_only" in s:  # after 43dd4ac: CENSUS_SENTENCE's "all with the mean's prediction" is wrong
        return sentence_capped_root_alone(s)
    un = s["roots"]["uncapped"]
    cap = s["roots"]["capped"]
    dmax = un["max_abs_difference_from_mean_log2"]
    return cs.CENSUS_SENTENCE.format(
        n=f"{s['eligible_for_a_fit']:,}",
        fitted=f"{un['fitted']:,}",
        x=f"{s['eligible_for_a_fit'] - un['fitted']:,}",
        reasons=sa._reasons(Counter(un["not_fitted"])),
        dmax="none" if dmax is None else f"{dmax:.1e}",
        ran=f"{un['run']:,}",
        miss=f"{un['harness_failing_genes'] + cap['harness_failing_genes']:,}",
        two=f"{s['two_roots']['genes']:,}",
        two_inf=f"{s['two_roots']['informative_on_the_capped_root']:,}",
        unique=f"{s['unique_fit']['genes']:,}",
        unique_inf=f"{s['unique_fit']['informative']:,}",
        common=f"{s['unique_fit']['informative_under_all_four']:,}",
    )


#: a gene's one fit, by its capped-sum status
UNIQUE = {"runnable": "uncapped", "capped_root_only": "capped"}


def _differs(g: dict) -> bool:
    """capsumrna.DECISION (1): the unique fit's prediction is not the mean's (or the mean has none)."""
    x = g["roots"][UNIQUE[g["status"]]]["difference_from_mean_log2"]
    return x is None or x > cs.CLOSED_FORM_TOLERANCE


def summarise_capped_root_alone(genes: list[dict], why: Counter) -> dict[str, Any]:
    """`summarise` for a census in which some gene has the capped root alone (capsumrna.STRUCTURE
    ['a_capped_root_alone']), whose difference from the mean is None; reads no RNA."""
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
    roots: dict[str, Any] = {}
    for root in ROOTS:
        ok = [g for g in elig if g["roots"][root]["ok"]]
        inf = [g for g in ok if g["roots"][root]["informative"]]
        ran = [g for g in inf if "harness_miss" in g["roots"][root]]
        diffs = [x for g in ok if (x := g["roots"][root]["difference_from_mean_log2"]) is not None]
        roots[root] = {
            "fitted": len(ok),
            "informative": len(inf),
            "informative_and_mean_informative": sum(1 for g in inf if informative(g, "mean")),
            "run": len(ran),
            "harness_failing_genes": sum(1 for g in ran if g["roots"][root]["harness_miss"]),
            "max_abs_run_minus_closed_form_log2": _max(
                [g["roots"][root]["max_abs_run_minus_closed_form_log2"] for g in ran]
            ),
            "mean_has_no_fit": len(ok) - len(diffs),
            "max_abs_difference_from_mean_log2": _max(diffs),
            "difference_from_mean_log2_quantiles": mr.quantiles(diffs) if diffs else None,
            "activating_count_quantiles": mr.quantiles([g["n_activating"] for g in ok]) if ok else None,
        }
    roots["uncapped"]["not_fitted"] = dict(
        Counter(g["roots"]["uncapped"]["reason"] for g in elig if not g["roots"]["uncapped"]["ok"])
    )
    unique = [g for g in elig if g["status"] in UNIQUE]
    uinf = [g for g in unique if _root_inf(g, UNIQUE[g["status"]])]
    differ = [g for g in unique if _differs(g)]
    common = [g for g in uinf if all(informative(g, f) for f in (*c2.FORMS, MULT))]
    two = [g for g in elig if g["status"] == "fit_not_unique"]
    alone = [g for g in elig if g["status"] == "capped_root_only"]
    ainf = [g for g in alone if _root_inf(g, "capped")]
    return {
        **dict(why),
        "classified": len(genes),
        "excluded_by_reason": dict(Counter(g["excluded"] for g in genes if g["excluded"])),
        "eligible_for_a_fit": len(elig),
        "forms": forms,
        "capped_sum_status": dict(sa._ranked(Counter(g["status"] for g in elig))),
        "roots": roots,
        "unique_fit": {
            "genes": len(unique),
            "by_status": dict(sa._ranked(Counter(g["status"] for g in unique))),
            "informative": len(uinf),
            "informative_symbol_on_one_gene_id": sum(1 for g in uinf if g["symbol"] in ids),
            "max_abs_difference_from_mean_log2": _max(
                [
                    x
                    for g in unique
                    if (x := g["roots"][UNIQUE[g["status"]]]["difference_from_mean_log2"]) is not None
                ]
            ),
            "differs_from_the_mean": {
                "genes": len(differ),
                "informative": sum(1 for g in differ if g in uinf),
                "informative_symbol_on_one_gene_id": sum(
                    1 for g in differ if g in uinf and g["symbol"] in ids
                ),
            },
            "informative_under_all_four": len(common),
            "informative_under_all_four_symbol_on_one_gene_id": sum(1 for g in common if g["symbol"] in ids),
            "uncapped_informative_share_of_the_means_informative_genes": round(
                sum(1 for g in uinf if g["status"] == "runnable") / max(1, forms["mean"]["informative"]), 4
            ),
        },
        "two_roots": {
            "genes": len(two),
            "informative_on_the_capped_root": sum(1 for g in two if _root_inf(g, "capped")),
            "informative_on_the_capped_root_symbol_on_one_gene_id": sum(
                1 for g in two if _root_inf(g, "capped") and g["symbol"] in ids
            ),
            "summed_drop_quantiles": mr.quantiles([g["sum_d"] for g in two]) if two else None,
        },
        "capped_root_only": {
            "genes": len(alone),
            "informative": len(ainf),
            "informative_symbol_on_one_gene_id": sum(1 for g in ainf if g["symbol"] in ids),
            "why_the_mean_has_no_fit": dict(sa._ranked(Counter(g["forms"]["mean"]["reason"] for g in alone))),
            "informative_under_competition": sum(1 for g in ainf if informative(g, "competition")),
            "informative_under_the_multiplicative_form": sum(1 for g in ainf if informative(g, MULT)),
            "activating_count_quantiles": mr.quantiles([g["n_activating"] for g in alone]) if alone else None,
            "with_an_inhibiting_rule": sum(1 for g in alone if g["n_inhibiting"]),
        },
        "rna_read": "none",
    }


def sentence_capped_root_alone(s: dict[str, Any]) -> str:
    """capsumrna.CENSUS_SENTENCE_CAPPED_ROOT_ALONE, filled from the counts only."""
    un = s["roots"]["uncapped"]
    cap = s["roots"]["capped"]
    dmax = un["max_abs_difference_from_mean_log2"]
    alone = s["capped_root_only"]
    return cs.CENSUS_SENTENCE_CAPPED_ROOT_ALONE.format(
        n=f"{s['eligible_for_a_fit']:,}",
        fitted=f"{un['fitted']:,}",
        x=f"{s['eligible_for_a_fit'] - un['fitted']:,}",
        reasons=sa._reasons(Counter(un["not_fitted"])),
        dmax="none" if dmax is None else f"{dmax:.1e}",
        ran=f"{un['run'] + cap['run']:,}",
        miss=f"{un['harness_failing_genes'] + cap['harness_failing_genes']:,}",
        two=f"{s['two_roots']['genes']:,}",
        two_inf=f"{s['two_roots']['informative_on_the_capped_root']:,}",
        only=f"{alone['genes']:,}",
        only_inf=f"{alone['informative']:,}",
        only_sym=f"{alone['informative_symbol_on_one_gene_id']:,}",
        only_reasons=sa._reasons(Counter(alone["why_the_mean_has_no_fit"])),
        unique=f"{s['unique_fit']['genes']:,}",
        differ=f"{s['unique_fit']['differs_from_the_mean']['genes']:,}",
        unique_inf=f"{s['unique_fit']['informative']:,}",
        common=f"{s['unique_fit']['informative_under_all_four']:,}",
    )


def committed_problems(
    s: dict[str, Any], c2_committed: dict[str, Any], mr_committed: dict[str, Any]
) -> list[str]:
    """Where this census's counts for 7d's and 7h's forms differ from their committed censuses."""
    rs = c2_committed["reader_selection"]
    out = [f"{k}: {s[k]} against {rs[k]}" for k in SEVEN_D_TOTALS if s[k] != rs[k]]
    for f in c2.FORMS:
        out += [
            f"{f} {k}: {s['forms'][f][k]} against {rs['forms'][f][k]}"
            for k in SEVEN_D_KEYS
            if s["forms"][f][k] != rs["forms"][f][k]
        ]
    ms = mr_committed["reader_selection"]["forms"][MULT]
    out += [
        f"{MULT} {k}: {s['forms'][MULT][k]} against {ms[k]}"
        for k in SEVEN_H_KEYS
        if s["forms"][MULT][k] != ms[k]
    ]
    return out


def census(chroms: list[str], workers: int, run: bool) -> dict:
    t = time.time()
    genes, why, programs = collect(chroms, run=run, workers=workers)
    s = summarise(genes, why)
    out = {
        "lane": cs.LANE,
        "date": time.strftime("%Y-%m-%d"),
        "roadmap_item": cs.ROADMAP_ITEM,
        "question": (
            "the multi-element genes the capped sum runs on 7d's design, beside 7d's and 7h's forms"
            " (no RNA read)"
        ),
        "census": cs.CENSUS,
        "form": cs.FORM,
        "form_source": cs.FORM_SOURCE,
        "conventions": cs.CONVENTIONS,
        "fit": cs.FIT,
        "structure": cs.STRUCTURE,
        "decision_rule": cs.DECISION,
        "exploration_preceded_this": cs.EXPLORATION_PRECEDED_THIS,
        "reader_selection": s,
        "sentence": sentence(s),
        "decision": cs.decide(s),
        "size_check_rule": cs.SIZE_CHECK,
        "size_check": cs.size_check(s),
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
            "form": cs.FORM,
            "comparator_forms": [*c2.FORMS, MULT],
            "h_intact": cs.H,
            "selection": cs.SELECTION,
            "runtime": cs.RUNTIME,
            "closed_form_tolerance": cs.CLOSED_FORM_TOLERANCE,
            "informative_log2": cs.INFORMATIVE_LOG2,
        },
        "exclusions": [
            "genes with one rule, rules in several programs, or a measured rule",
            "observations unreadable or exactly 1, actions that disagree with the observation",
            "genes with no activating rule; under the capped sum the mean's reasons, sum_past_the_cap, and"
            " fit_not_unique (counted, both roots run)",
            "capped_root_only genes (after 43dd4ac) are counted as unique fits whose prediction is not the"
            " mean's",
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
    print(
        json.dumps(
            {
                k: s[k]
                for k in ("eligible_for_a_fit", "capped_sum_status", "roots", "unique_fit", "two_roots")
            },
            indent=1,
        )
    )
    print(json.dumps(s["forms"], indent=1))
    print(out["sentence"])
    print("decision:", out["decision"])
    print("capped root alone:", s.get("capped_root_only"), "size check:", out["size_check"])
    if "committed_censuses_reproduced" in out:
        print("committed census problems:", out["committed_censuses_reproduced"]["problems"])
    if not args.no_save:
        print("saved", save_result(RESULT, out, manifest=manifest(args.chroms)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
