# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-residue census: how the summed drop is distributed WITHIN the genes that stay unfittable.

    uv run python scripts/residue_census.py [--chroms chr21] [--no-save] [--workers 4]

The questions, the per-element cap ladder, the cap's arithmetic bound and the decision rule were committed
first, in `genomeos/attribution/residue.py` and `data/results/residue_questions.json`. This script computes
the answers that record asks for and nothing else: it gates and parses the same 24 compiled programs 7d and
`lane-scalefit` ran, re-derives the mean form's `basal_not_positive` population and the residue of it that
is still `basal_not_positive` at M2h's fixed scale, and reports for every one of those genes how much of
its summed activating drop its largest element carries, its largest three carry, how diffusely the sum is
spread, how few elements would have to go to bring it under the gene's level, and how many genes a
per-element cap at each declared level would make fittable.

NO RNA VALUE IS OPENED. NOTHING IS FITTED, TUNED OR SELECTED. No model request is sent and no key is read.
A gene a fit can accept is not a gene the model predicts: this census makes no claim about expression, and
it proposes no combination form.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from multiprocessing import get_context
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import combine2 as c2  # noqa: E402
from genomeos.attribution import compile as cp  # noqa: E402
from genomeos.attribution import readergate as rg  # noqa: E402
from genomeos.attribution import residue as rs  # noqa: E402
from genomeos.attribution import scalefit as sf  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.results import save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readergate_census import peak_files  # noqa: E402

RESULT = "residue_census"
GENES_RESULT = "residue_genes"
QUESTIONS_FILE = Path("data/results/residue_questions.json")
SCALEFIT_CENSUS = Path("data/results/scalefit_census.json")
CENSUS_7D = Path("data/results/combine2_census.json")
MAGCAL = Path("data/results/magcal_result.json")
CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]
FORM = "mean"  # the runtime default and scalefit's primary population (residue.POPULATION)
OWN_CODE = (
    "genomeos/attribution/residue.py",
    "scripts/residue_register.py",
    "scripts/residue_census.py",
    "tests/test_residue.py",
)
COLUMNS = (
    "chrom",
    "n_rules",
    "n_activating",
    "n_inhibiting",
    "sum_d",
    "q_after",
    "in_residue",
    "top1",
    "top3",
    "n_eff",
    "k_min",
    *(f"fits_under_{c}" for c in rs.CAP_LEVELS),
)


def texts_of(chroms: list[str]) -> dict[str, str]:
    return {c: (Path(COMPILED_DIR) / f"noncoding_{c}.bio").read_text() for c in chroms}


def _chromosome(args: tuple[str, set[str]]) -> dict[str, Any]:
    """One program: 7d's classification, the scaled refit, and the declared statistics of the excluded."""
    from genomeos.lang.parser import parse

    chrom, genes = args
    t0 = time.time()
    text = (Path(COMPILED_DIR) / f"noncoding_{chrom}.bio").read_text()
    gated, info = cp.gate_program(text, chrom, selection=c2.SELECTION)
    module = parse(gated, f"noncoding_{chrom}_reader")
    out: dict[str, Any] = {
        "chrom": chrom,
        "program": {k: info[k] for k in ("base_sha256", "gated_sha256", "elements", "rules_gated")},
        "eligible": 0,
        "excluded": 0,
        "rows": [],
        "max_abs_q_after_minus_scale_times_q_before": 0.0,
        "max_abs_capped_sum_minus_combine2_q": 0.0,
        "cap_disagreements": 0,
    }
    for gr in c2.gene_rules(module, genes, chrom, rg.compiled_cells(text)):
        _fits, why = c2.classify(gr)
        if why:
            continue
        out["eligible"] += 1
        before = c2.fit(FORM, gr.rho)
        if before.ok or before.reason != "basal_not_positive":
            continue
        out["excluded"] += 1
        after = c2.fit(FORM, sf.scaled_rho(gr.rho))
        in_residue = not (after.q < 1.0)
        row = rs.gene_residue(gr.symbol, chrom, gr.rho, in_residue)
        # the identity scalefit proved, checked here on every gene rather than relied on
        out["max_abs_q_after_minus_scale_times_q_before"] = max(
            out["max_abs_q_after_minus_scale_times_q_before"], abs(after.q - rs.SCALE * before.q)
        )
        # CAP_RULE's cross-check: the capped sum IS 7d's q on the capped observations
        drops = rs.activating_drops(gr.rho)
        for cap in rs.CAP_LEVELS:
            capped = c2.fit(FORM, rs.capped_rho(gr.rho, cap))
            out["max_abs_capped_sum_minus_combine2_q"] = max(
                out["max_abs_capped_sum_minus_combine2_q"],
                abs(capped.q - rs.capped_sum(drops, cap)),
            )
            if rs.fits_under_cap(drops, cap) != (capped.q < 1.0):
                out["cap_disagreements"] += 1
        out["rows"].append(
            [
                chrom,
                row.n_rules,
                row.n_activating,
                row.n_inhibiting,
                row.sum_d,
                row.q_after,
                int(row.in_residue),
                row.top1,
                row.top3,
                row.n_eff,
                row.k_min,
                *(int(row.fits_under[str(c)]) for c in rs.CAP_LEVELS),
            ]
        )
    out["seconds"] = round(time.time() - t0, 1)
    return out


def collect(chroms: list[str], workers: int) -> dict[str, Any]:
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
                print(r["chrom"], r["eligible"], r["excluded"], r["seconds"], flush=True)
    else:
        for j in jobs:
            r = _chromosome(j)
            res[r["chrom"]] = r
            print(r["chrom"], r["eligible"], r["excluded"], r["seconds"], flush=True)
    return {
        "rows": [row for c in chroms for row in res[c]["rows"]],
        "eligible": sum(res[c]["eligible"] for c in chroms),
        "excluded": sum(res[c]["excluded"] for c in chroms),
        "programs": {c: res[c]["program"] for c in chroms},
        "max_abs_q_after_minus_scale_times_q_before": max(
            res[c]["max_abs_q_after_minus_scale_times_q_before"] for c in chroms
        ),
        "max_abs_capped_sum_minus_combine2_q": max(
            res[c]["max_abs_capped_sum_minus_combine2_q"] for c in chroms
        ),
        "cap_disagreements": sum(res[c]["cap_disagreements"] for c in chroms),
    }


def col(rows: list[list], name: str) -> list[float]:
    i = COLUMNS.index(name)
    return [float(r[i]) for r in rows]


def statistics(rows: list[list], label: str) -> dict[str, Any]:
    """Every declared distribution for one population (Q2 to Q6, Q8). The shape, never only a median."""
    kmin = col(rows, "k_min")
    n_act = col(rows, "n_activating")
    n_eff = col(rows, "n_eff")
    q_after = col(rows, "q_after")
    top3_four_plus = [float(r[COLUMNS.index("top3")]) for r in rows if r[COLUMNS.index("n_activating")] >= 4]
    return {
        "population": label,
        "genes": len(rows),
        "Q2_n_rules": rs.distribution(col(rows, "n_rules")),
        "Q2_n_activating": rs.distribution(n_act),
        "Q2_n_inhibiting": rs.distribution(col(rows, "n_inhibiting")),
        "Q3_top1_share": rs.distribution(col(rows, "top1")),
        "Q4_top3_share": rs.distribution(col(rows, "top3")),
        "Q4_top3_share_of_the_genes_with_four_or_more_activating_elements": rs.distribution(top3_four_plus),
        "Q5_q_after_scaling": rs.distribution(q_after),
        "Q5_sum_of_drops_before_scaling": rs.distribution(col(rows, "sum_d")),
        "Q5_at_or_above_1_after": sum(1 for x in q_after if x >= 1.0),
        "Q5_at_or_above_2_after": sum(1 for x in q_after if x >= 2.0),
        "Q5_at_or_above_5_after": sum(1 for x in q_after if x >= 5.0),
        "Q6_k_min": rs.distribution(kmin),
        "Q6_k_min_counts": {
            "0": sum(1 for x in kmin if x == 0.0),
            "1": sum(1 for x in kmin if x == 1.0),
            "2": sum(1 for x in kmin if x == 2.0),
            "3": sum(1 for x in kmin if x == 3.0),
            "between_1_and_3": sum(1 for x in kmin if 1.0 <= x <= 3.0),
            "above_3": sum(1 for x in kmin if x > 3.0),
            "no_k_below_n_works": sum(1 for x in kmin if x != x),
        },
        "Q8_n_eff": rs.distribution(n_eff),
        "Q8_n_eff_over_n_activating": rs.distribution(
            [e / n for e, n in zip(n_eff, n_act, strict=True) if n > 0]
        ),
    }


def cap_curve(rows: list[list]) -> dict[str, Any]:
    """Q7: the whole curve over the declared ladder, no level selected after the fact."""
    out: dict[str, Any] = {
        "rule": rs.CAP_RULE,
        "levels_declared_before_the_count": list(rs.CAP_LEVELS),
        "primary_level": rs.CAP_PRIMARY,
        "bound": rs.CAP_BOUND,
        "genes": len(rows),
        "per_level": {},
    }
    n_act = [float(r[COLUMNS.index("n_activating")]) for r in rows]
    for cap in rs.CAP_LEVELS:
        i = COLUMNS.index(f"fits_under_{cap}")
        n_fit = sum(1 for r in rows if r[i])
        out["per_level"][str(cap)] = {
            "recovered": n_fit,
            "share": (n_fit / len(rows)) if rows else None,
            "still_not_fittable": len(rows) - n_fit,
            "genes_below_the_bound_1_over_L": sum(1 for n in n_act if n < 1.0 / cap),
            "bound_note": (
                f"a cap at {cap} cannot fit a gene with {round(1.0 / cap)} or more activating elements"
                f" and fits every gene with fewer, whatever its drops"
            ),
        }
    return out


def reconcile(eligible: int, excluded: int, residue: int, recovered: int) -> dict[str, Any]:
    """This run's re-derived counts against scalefit's and 7d's committed ones."""
    m = json.loads(SCALEFIT_CENSUS.read_text())["per_form"][FORM]
    sel = json.loads(CENSUS_7D.read_text())["reader_selection"]
    out: dict[str, Any] = {
        "note": (
            "every count on the right was re-derived in THIS run from the 24 compiled programs with"
            " combine2's own code; the committed numbers are read only to compare. A disagreement is a"
            " finding of this lane and is reported, not repaired"
        ),
        "eligible_for_a_fit": {"committed_7d": sel["eligible_for_a_fit"], "re_derived": eligible},
        "excluded_basal_not_positive": {
            "committed_scalefit": m["excluded_basal_not_positive"],
            "committed_7d": sel["forms"][FORM]["not_runnable"]["basal_not_positive"],
            "re_derived": excluded,
        },
        "recovered_at_the_scale": {"committed_scalefit": m["recovered"], "re_derived": recovered},
        "the_residue": {"committed_scalefit": m["still_basal_not_positive"], "re_derived": residue},
    }
    out["agrees"] = (
        out["eligible_for_a_fit"]["committed_7d"] == eligible
        and out["excluded_basal_not_positive"]["committed_scalefit"] == excluded
        and out["excluded_basal_not_positive"]["committed_7d"] == excluded
        and out["recovered_at_the_scale"]["committed_scalefit"] == recovered
        and out["the_residue"]["committed_scalefit"] == residue
    )
    return out


def census(chroms: list[str], workers: int) -> tuple[dict[str, Any], dict[str, Any]]:
    got = collect(chroms, workers)
    rows = got["rows"]
    i_res = COLUMNS.index("in_residue")
    residue = [r for r in rows if r[i_res]]
    recovered = [r for r in rows if not r[i_res]]
    med_top1 = rs.median(col(residue, "top1"))
    med_top3 = rs.median(col(residue, "top3"))
    med_kmin = rs.median(col(residue, "k_min"))
    branch = rs.decide(med_top1, med_top3, med_kmin)
    out = {
        "lane": rs.LANE,
        "date": time.strftime("%Y-%m-%d"),
        "registration": str(QUESTIONS_FILE),
        "question": rs.QUESTION,
        "form": FORM,
        "form_is": "the runtime default and scalefit's primary population; competition is not counted here",
        "scale": rs.SCALE,
        "scale_rounded": rs.SCALE_ROUNDED,
        "scale_is_refitted": False,
        "residue_threshold_on_the_unscaled_ratio": rs.RESIDUE_THRESHOLD,
        "Q1_counts": {
            "eligible_for_a_fit": got["eligible"],
            "excluded_basal_not_positive": got["excluded"],
            "recovered_at_the_scale": len(recovered),
            "the_residue": len(residue),
        },
        "decision_rule": rs.DECISION_RULE,
        "decision_thresholds": {
            "concentrated_median_top1_at_or_above": rs.CONCENTRATED_TOP1,
            "few_median_top3_at_or_above": rs.FEW_TOP3,
            "few_median_k_min_at_or_below": rs.FEW_KMIN,
            "diffuse_median_top1_below": rs.DIFFUSE_TOP1,
        },
        "decision_inputs": {
            "median_top1_of_the_residue": med_top1,
            "median_top3_of_the_residue": med_top3,
            "median_k_min_of_the_residue": med_kmin,
        },
        "branch": branch,
        "verdict": rs.reading(
            branch,
            residue=len(residue),
            median_top1=round(med_top1, 4),
            median_top3=round(med_top3, 4),
            median_k_min=med_kmin,
        ),
        "residue": statistics(residue, "the residue: still basal_not_positive at the scale"),
        "recovered": statistics(recovered, "the recovered: a positive basal rate again at the scale"),
        "Q7_per_element_cap": cap_curve(residue),
        "Q7_per_element_cap_on_the_recovered_for_comparison": cap_curve(recovered),
        "reconciliation": reconcile(got["eligible"], got["excluded"], len(residue), len(recovered)),
        "identity_checks": {
            "max_abs_q_after_minus_scale_times_q_before": got["max_abs_q_after_minus_scale_times_q_before"],
            "max_abs_capped_sum_minus_combine2_q": got["max_abs_capped_sum_minus_combine2_q"],
            "cap_verdict_disagreements_with_combine2": got["cap_disagreements"],
        },
        "what_this_is_not": (
            "not a claim about expression and not a proposed form. Every count here is about whether 7d's"
            " fit has a positive basal rate left, which is arithmetic. Whether any of these genes would"
            " read closer to measured RNA under any form is a SEPARATE registration on separate data, as"
            " lane-recovered had to be for the genes that recovered at the scale alone. NO RNA VALUE WAS"
            " OPENED IN THIS LANE"
        ),
        "prior_exposure": rs.PRIOR_EXPOSURE,
        "limitations": list(rs.LIMITATIONS),
        "conventions": rs.CONVENTIONS,
        "per_gene_table": f"data/results/{GENES_RESULT}.json",
        "programs": got["programs"],
        "rna_read": "none",
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
    }
    genes = {
        "lane": rs.LANE,
        "date": out["date"],
        "of": f"data/results/{RESULT}.json",
        "form": FORM,
        "columns": list(COLUMNS),
        "note": (
            "one row per gene the mean form excluded as basal_not_positive, re-derived in this run."
            " in_residue is 1 when the gene is STILL basal_not_positive at M2h's fixed scale. The"
            " fits_under_* columns are the per-element cap verdicts at the levels declared before the"
            " count. No RNA value is in this table and none was opened"
        ),
        "rows": rows,
    }
    return out, genes


def manifest(chroms: list[str]) -> dict[str, Any]:
    inputs = [
        *(mf.input_entry(Path(COMPILED_DIR) / f"noncoding_{c}.bio", partition=None) for c in chroms),
        mf.input_entry(QUESTIONS_FILE, partition=None),
        mf.input_entry(SCALEFIT_CENSUS, partition=None),
        mf.input_entry(CENSUS_7D, partition=None),
        mf.input_entry(MAGCAL, partition=None),
        mf.input_entry(rg.GENE_TABLE, partition=None),
        mf.files_entry("reader v1 DNase narrowPeak files", [p for p in peak_files(chroms) if p.exists()]),
    ]
    return {
        "sources": [
            {
                "accession": "the 24 compiled non-coding programs M1, M2b, 7d and scalefit ran",
                "version": "pinned by sha256",
            },
            {"accession": "data/results/residue_questions.json (this lane's registration)", "version": "dev"},
            {"accession": "data/results/scalefit_census.json (lane-scalefit)", "version": "7f65544"},
            {"accession": "data/results/combine2_census.json (7d/7e)", "version": "948cfc9"},
            {"accession": "data/results/magcal_result.json (M2h, the scale)", "version": "ac5580b"},
            {"accession": "reader v1 ENCODE DNase-seq narrowPeak files, 13 biosamples", "version": "pinned"},
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "form": FORM,
            "cap_levels": list(rs.CAP_LEVELS),
            "cap_primary": rs.CAP_PRIMARY,
            "quantiles": list(rs.QUANTILES),
            "selection": c2.SELECTION,
            **rs.CONVENTIONS,
        },
        "exclusions": [
            "7d's own: one rule, rules in several programs, a measured rule, an unreadable or exactly-1"
            " observation, an action that disagrees with it, no activating rule",
            "genes the mean form excluded for a reason other than basal_not_positive are not in this"
            " census's population at all",
            "the promoter-competition form is not counted here",
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
    out, genes = census(args.chroms, args.workers)
    out["seconds"] = round(time.time() - t0, 1)
    print("Q1", out["Q1_counts"], flush=True)
    keys = ("Q2_n_activating", "Q3_top1_share", "Q4_top3_share", "Q5_q_after_scaling", "Q6_k_min", "Q8_n_eff")
    for name in ("residue", "recovered"):
        d = out[name]
        print(name, d["genes"], flush=True)
        for key in keys:
            print("  ", key, d[key].get("quantiles"), flush=True)
        print("   Q6 counts", d["Q6_k_min_counts"], flush=True)
    print("Q7 cap curve on the residue:", flush=True)
    for cap, v in out["Q7_per_element_cap"]["per_level"].items():
        print(f"   L={cap}: {v['recovered']} of {out['Q7_per_element_cap']['genes']}", flush=True)
    print("reconciliation", out["reconciliation"]["agrees"], flush=True)
    print("identity checks", out["identity_checks"], flush=True)
    print("BRANCH", out["branch"], flush=True)
    print("VERDICT", out["verdict"], flush=True)
    if not args.no_save:
        print("saved", save_result(RESULT, out, manifest=manifest(args.chroms)))
        print("saved", save_result(GENES_RESULT, genes, manifest=manifest(args.chroms), compact=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
