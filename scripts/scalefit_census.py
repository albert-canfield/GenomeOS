# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-scalefit census: how many of 7d's `basal_not_positive` genes get a basal rate back at M2h's scale.

    uv run python scripts/scalefit_census.py [--chroms chr21] [--no-save] [--workers 4]

The rule and the floor were committed first, in `genomeos/attribution/scalefit.py` and
`data/results/scalefit_rule.json`. This script computes the one answer that record asks for: it gates and
parses the same 24 compiled programs 7d ran, classifies the same multi-element genes with `combine2`'s own
code, and for every gene 7d excluded as `basal_not_positive` refits it on observations whose single-removal
effects are multiplied by M2h's FIXED constant `scalefit.SCALE`.

NO RNA VALUE IS OPENED. NOTHING IS FITTED, TUNED OR SELECTED. No model request is sent and no key is read.
A recovered gene is a gene 7d's fit can now accept, nothing more: this census makes no claim whatever
about expression.
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
from genomeos.attribution import readergate as rg  # noqa: E402
from genomeos.attribution import scalefit as sf  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.results import save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readergate_census import peak_files  # noqa: E402

RESULT = "scalefit_census"
RULE_FILE = Path("data/results/scalefit_rule.json")
CENSUS_7D = Path("data/results/combine2_census.json")
MAGCAL = Path("data/results/magcal_result.json")
GENE_TABLE = Path("data/results/gene_rates.json")
CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]
OWN_CODE = (
    "genomeos/attribution/scalefit.py",
    "scripts/scalefit_register.py",
    "scripts/scalefit_census.py",
    "tests/test_scalefit.py",
)


def texts_of(chroms: list[str]) -> dict[str, str]:
    return {c: (Path(COMPILED_DIR) / f"noncoding_{c}.bio").read_text() for c in chroms}


def _chromosome(args: tuple[str, set[str]]) -> dict[str, Any]:
    """One program's multi-element genes: 7d's classification, then the scaled refit of the excluded."""
    from genomeos.lang.parser import parse

    chrom, genes = args
    t0 = time.time()
    text = (Path(COMPILED_DIR) / f"noncoding_{chrom}.bio").read_text()
    gated, info = cp.gate_program(text, chrom, selection=c2.SELECTION)
    module = parse(gated, f"noncoding_{chrom}_reader")
    out: dict[str, Any] = {
        "chrom": chrom,
        "program": {k: info[k] for k in ("base_sha256", "gated_sha256", "elements", "rules_gated")},
        "genes": [],
        "max_abs_q_minus_scale_times_q": 0.0,
    }
    for gr in c2.gene_rules(module, genes, chrom, rg.compiled_cells(text)):
        fits, why = c2.classify(gr)
        row: dict[str, Any] = {
            "symbol": gr.symbol,
            "n_rules": len(gr.elements),
            "n_inhibiting": sum(1 for r in gr.rho.values() if r > 1.0),
            "excluded": why,
            "sum_d": sf.ratio(gr.rho) if not why else None,
            "forms": {},
        }
        if not why:
            scaled = sf.scaled_rho(gr.rho)
            for form in c2.FORMS:
                rec = sf.recovery(form, gr.rho)
                fr = rec.as_row()
                if rec.recovered:
                    gap = sf.q_is_homogeneous(form, gr.rho)
                    out["max_abs_q_minus_scale_times_q"] = max(out["max_abs_q_minus_scale_times_q"], gap)
                if rec.runnable:
                    after = c2.fit(form, scaled)
                    s = c2.predicted_log2(after, gr)
                    fr["range_log2_after"] = max(s) - min(s)
                    fr["informative_after"] = max(s) - min(s) > c2.INFORMATIVE_LOG2
                row["forms"][form] = fr
        out["genes"].append(row)
    out["seconds"] = round(time.time() - t0, 1)
    return out


def collect(chroms: list[str], workers: int) -> tuple[list[dict], dict, float]:
    multi, _why = c2.multi_rule_genes(texts_of(chroms))
    by_chrom: dict[str, set[str]] = {c: set() for c in chroms}
    for g, c in multi.items():
        by_chrom[c].add(g)
    order = sorted(chroms, key=lambda c: -(Path(COMPILED_DIR) / f"noncoding_{c}.bio").stat().st_size)
    jobs = [(c, by_chrom[c]) for c in order]
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
    gap = max(res[c]["max_abs_q_minus_scale_times_q"] for c in chroms)
    return genes, {c: res[c]["program"] for c in chroms}, gap


def per_form(genes: list[dict], form: str) -> dict[str, Any]:
    """Everything the rule declares for one form, over the genes 7d excluded as basal_not_positive."""
    elig = [g for g in genes if not g["excluded"]]
    excl = [g for g in elig if g["forms"][form]["before_reason"] == "basal_not_positive"]
    rec = [g for g in excl if g["forms"][form]["recovered"]]
    run = [g for g in rec if g["forms"][form]["runnable"]]
    inf = [g for g in run if g["forms"][form].get("informative_after")]
    return {
        "eligible_for_a_fit": len(elig),
        "excluded_basal_not_positive": len(excl),
        "recovered": len(rec),
        "recovered_share_of_the_excluded": (len(rec) / len(excl)) if excl else None,
        "runnable_after": len(run),
        "lost_to_inhibitor_strength_above_1_after_recovery": len(rec) - len(run),
        "informative_after": len(inf),
        "still_basal_not_positive": len(excl) - len(rec),
        "after_reason_of_the_not_recovered": dict(
            Counter(g["forms"][form]["after_reason"] for g in excl if not g["forms"][form]["recovered"])
        ),
        "after_reason_of_the_recovered": dict(Counter(g["forms"][form]["after_reason"] for g in rec)),
        "ratio_sum_of_drops_over_level": {
            "over": "the genes this form excluded as basal_not_positive",
            "before_quantiles": sf.quantiles([g["sum_d"] for g in excl]),
            "after_quantiles": sf.quantiles([sf.SCALE * g["sum_d"] for g in excl]),
            "mean_before": (sum(g["sum_d"] for g in excl) / len(excl)) if excl else None,
            "mean_after": (sf.SCALE * sum(g["sum_d"] for g in excl) / len(excl)) if excl else None,
            "at_or_above_1_before": sum(1 for g in excl if g["sum_d"] >= 1.0),
            "at_or_above_1_after": sum(1 for g in excl if sf.SCALE * g["sum_d"] >= 1.0),
            "max_before": max((g["sum_d"] for g in excl), default=None),
        },
        "q_of_the_fit": {
            "note": (
                "q is the element-driven share of the intact level: the sum of drops itself under the"
                " mean, Z*THETA under competition. The exclusion test is q >= 1 in both"
            ),
            "before_quantiles": sf.quantiles([g["forms"][form]["q_before"] for g in excl]),
            "after_quantiles": sf.quantiles([g["forms"][form]["q_after"] for g in excl]),
            "at_or_above_1_after": sum(1 for g in excl if g["forms"][form]["q_after"] >= 1.0),
        },
        "rule_count_quantiles_of_the_recovered": sf.quantiles([float(g["n_rules"]) for g in rec]),
        "rule_count_quantiles_of_the_not_recovered": sf.quantiles(
            [float(g["n_rules"]) for g in excl if not g["forms"][form]["recovered"]]
        ),
    }


def reconcile(forms: dict[str, Any]) -> dict[str, Any]:
    """This census's excluded counts against 7d's committed ones: they must be identical."""
    rs = json.loads(CENSUS_7D.read_text())["reader_selection"]
    out: dict[str, Any] = {"source": str(CENSUS_7D), "agrees": True}
    for f in c2.FORMS:
        theirs = rs["forms"][f]["not_runnable"]["basal_not_positive"]
        mine = forms[f]["excluded_basal_not_positive"]
        out[f] = {"7d": theirs, "here": mine, "same": theirs == mine}
        out["agrees"] = out["agrees"] and theirs == mine
    out["eligible_7d"] = rs["eligible_for_a_fit"]
    out["eligible_here"] = forms[c2.FORMS[0]]["eligible_for_a_fit"]
    out["agrees"] = out["agrees"] and out["eligible_7d"] == out["eligible_here"]
    return out


def census(chroms: list[str], workers: int) -> dict:
    genes, programs, gap = collect(chroms, workers)
    forms = {f: per_form(genes, f) for f in c2.FORMS}
    mean = forms["mean"]
    return {
        "lane": sf.LANE,
        "date": time.strftime("%Y-%m-%d"),
        "registration": str(RULE_FILE),
        "question": json.loads(RULE_FILE.read_text())["question"],
        "scale": sf.SCALE,
        "scale_rounded": sf.SCALE_ROUNDED,
        "scale_source": sf.SCALE_SOURCE,
        "scale_is_refitted": False,
        "recovery_threshold_on_the_unscaled_ratio": sf.RECOVERY_THRESHOLD,
        "floor": sf.FLOOR,
        "strong_bar": sf.STRONG_BAR,
        "verdict": sf.verdict(mean["recovered"], mean["excluded_basal_not_positive"]),
        "verdict_is_on": "the mean form, the runtime default and 7d's primary population",
        "per_form": forms,
        "reconciliation_with_7d": reconcile(forms),
        "max_abs_q_after_minus_scale_times_q_before": gap,
        "what_this_is_not": (
            "not a claim about expression. A recovered gene is a gene 7d's fit can now accept. Whether"
            " the scaled run is closer to measured RNA is a SEPARATE registration and a separate lane,"
            " and this lane opened no RNA value at any point"
        ),
        "prior_exposure": sf.PRIOR_EXPOSURE,
        "limitations": list(sf.LIMITATIONS),
        "conventions": sf.CONVENTIONS,
        "programs": programs,
        "rna_read": "none",
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
    }


def manifest(chroms: list[str]) -> dict:
    inputs = [
        *(mf.input_entry(Path(COMPILED_DIR) / f"noncoding_{c}.bio", partition=None) for c in chroms),
        mf.input_entry(RULE_FILE, partition=None),
        mf.input_entry(CENSUS_7D, partition=None),
        mf.input_entry(MAGCAL, partition=None),
        mf.input_entry(rg.GENE_TABLE, partition=None),
        mf.files_entry("reader v1 DNase narrowPeak files", [p for p in peak_files(chroms) if p.exists()]),
    ]
    return {
        "sources": [
            {
                "accession": "the 24 compiled non-coding programs M1, M2b and 7d ran",
                "version": "pinned by sha256",
            },
            {"accession": "data/results/combine2_census.json (7d/7e)", "version": "948cfc9"},
            {"accession": "data/results/magcal_result.json (M2h, the scale)", "version": "ac5580b"},
            {"accession": "reader v1 ENCODE DNase-seq narrowPeak files, 13 biosamples", "version": "pinned"},
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "scale": sf.SCALE,
            "scale_is_refitted": False,
            "floor": sf.FLOOR,
            "strong_bar": sf.STRONG_BAR,
            "selection": c2.SELECTION,
            **sf.CONVENTIONS,
        },
        "exclusions": [
            "7d's own: one rule, rules in several programs, a measured rule, an unreadable or exactly-1"
            " observation, an action that disagrees with it, no activating rule",
            "genes 7d excluded for a reason other than basal_not_positive are not counted as recovered",
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
    for f, d in out["per_form"].items():
        print(f, {k: v for k, v in d.items() if isinstance(v, (int, float, type(None)))}, flush=True)
        print("  ratio before", d["ratio_sum_of_drops_over_level"]["before_quantiles"], flush=True)
        print("  ratio after ", d["ratio_sum_of_drops_over_level"]["after_quantiles"], flush=True)
    print("reconciliation", out["reconciliation_with_7d"], flush=True)
    print("VERDICT", out["verdict"], flush=True)
    if not args.no_save:
        print("saved", save_result(RESULT, out, manifest=manifest(args.chroms)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
