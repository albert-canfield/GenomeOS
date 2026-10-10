# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-multrna, the registered run: the multiplicative composition against measured RNA, beside 7d.

    uv run python scripts/multrna_run.py [--workers 4]

Runs only from the committed registration (data/results/multrna_registration.json); refuses if any
registered field differs from genomeos/attribution/multrna.py, any pinned file's sha256 differs from
the file on disk, or 7d's own registration check (scripts/combine2_run.check_registration) refuses.

Order, as registered (multrna.DECISION):
  1. 7d's two forms are classified, fitted and run by 7d's committed code (scripts/combine2_census.collect);
     the multiplicative form by this lane's census code (scripts/multrna_census.collect), on the same
     gated programs;
  2. the HARNESS is checked on every gene run of all three forms, and the multiplicative side must equal
     the committed census; any miss stops the run here, before the RNA opens (HARNESS_DEFECT);
  3. the RNA quantifications are opened and the three populations assembled; 7d's populations, P, D and
     registered reading must reproduce the committed result (COMPARATOR_DEFECT otherwise);
  4. only then is the multiplicative form scored and the two head-to-heads read.

No model request, no network, no money.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import combine2 as c2  # noqa: E402
from genomeos.attribution import multrna as mr  # noqa: E402
from genomeos.attribution import readergate as rg  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import combine2_census as c2census  # noqa: E402
import combine2_run as c2run  # noqa: E402
import multrna_census as mcensus  # noqa: E402
from readergate_census import peak_files  # noqa: E402
from readergate_run import load_rna  # noqa: E402
from runrules_census import RATE_TABLE  # noqa: E402

RESULT = "multrna_result"
FORM = mr.FORM
REGISTRATION = RESULTS_DIR / "multrna_registration.json"
CENSUS = RESULTS_DIR / "multrna_census.json"
C2_REGISTRATION = c2run.REGISTRATION
C2_CENSUS = c2run.CENSUS
C2_RESULT = RESULTS_DIR / "combine2_result.json"
M2B_RESULT = c2run.M2B_RESULT
M2E_RESULT = RESULTS_DIR / "saturate_result.json"
CHROMS = c2census.CHROMS
#: the census counts the run's predictor side must equal before the RNA is opened
CENSUS_KEYS = (
    "classified",
    "eligible_for_a_fit",
    "forms",
    "informative_under_all_three",
    "new_under_the_multiplicative_form",
)


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def check_registration() -> dict:
    """The committed registration, refused if a bound field differs from the code or a pin from disk."""
    reg = json.loads(REGISTRATION.read_text())
    fields = mr.registered_fields()
    bad = [k for k in mr.BOUND if reg.get(k) != fields[k]]
    if bad:
        raise SystemExit(f"refused: the registration and multrna.py disagree on {', '.join(bad)}")
    group = mf.files_entry("reader v1 DNase narrowPeak files", [p for p in peak_files(CHROMS) if p.exists()])
    moved = [
        p
        for p, digest in reg["pinned_sha256"].items()
        if p != group["path"] and (not Path(p).is_file() or _sha(Path(p)) != digest)
    ]
    if reg["pinned_sha256"].get(group["path"]) != group["sha256"]:
        moved.append(group["path"])
    if moved:
        raise SystemExit(f"refused: pinned inputs changed or missing since the registration: {moved}")
    return reg


def mult_runs(genes: list[dict]) -> list[dict[str, Any]]:
    """The records multrna.assemble reads, for the genes runnable and informative under the form."""
    out = []
    for g in genes:
        fr = g["forms"].get(FORM)
        if g["excluded"] or not fr or not fr["ok"] or not fr["informative"]:
            continue
        out.append(
            {
                **{
                    k: g[k]
                    for k in ("gene", "symbol", "chrom", "n_rules", "n_inhibiting", "elements", "states")
                },
                "open_share": g["open_share"],
                "s": fr["s"],
                "log2_rho": fr["log2_rho"],
            }
        )
    return out


def harness_check(genes7: list[dict], genes_m: list[dict]) -> dict[str, Any]:
    ran7 = [(g, f) for g in genes7 for f, fr in g["forms"].items() if "s" in fr]
    ran_m = [g for g in genes_m if "s" in g["forms"].get(FORM, {})]
    bad = [(g["gene"], f) for g, f in ran7 if g["forms"][f]["harness_miss"]]
    bad += [(g["gene"], FORM) for g in ran_m if g["forms"][FORM]["harness_miss"]]
    worst = max(
        (
            abs(a - b)
            for fr in [g["forms"][f] for g, f in ran7] + [g["forms"][FORM] for g in ran_m]
            for a, b in zip(fr["s"], fr["closed_s"], strict=True)
        ),
        default=None,
    )
    return {
        "gene_form_runs": len(ran7) + len(ran_m),
        "seven_d_gene_form_runs": len(ran7),
        "multiplicative_gene_runs": len(ran_m),
        "failing": len(bad),
        "failing_examples": bad[:20],
        "max_abs_run_minus_closed_form_log2": worst,
    }


def predictor_problems(
    genes7: list[dict], genes_m: list[dict], p7: dict, pm: dict, summary: dict
) -> list[str]:
    """Before the RNA: both collections read the same programs and genes, and the multiplicative side
    equals the committed census."""
    out = []
    if p7 != pm:
        out.append("7d's and this lane's gated programs differ")
    if [(g["gene"], g["excluded"]) for g in genes7] != [(g["gene"], g["excluded"]) for g in genes_m]:
        out.append("7d's and this lane's multi-element genes or exclusions differ")
    committed = json.loads(CENSUS.read_text())["reader_selection"]
    out += [
        f"census {k} differs from the committed census"
        for k in CENSUS_KEYS
        if json.loads(json.dumps(summary[k])) != committed[k]
    ]
    return out


def manifest(rna_read: bool) -> dict:
    reg = json.loads(REGISTRATION.read_text())
    inputs = [
        *(mf.input_entry(Path(COMPILED_DIR) / f"noncoding_{c}.bio", partition=None) for c in CHROMS),
        mf.input_entry(RATE_TABLE, partition=None),
        mf.input_entry(rg.GENE_TABLE, partition=None),
        mf.files_entry("reader v1 DNase narrowPeak files", [p for p in peak_files(CHROMS) if p.exists()]),
        mf.input_entry(REGISTRATION, partition=None),
        mf.input_entry(CENSUS, partition=None),
        mf.input_entry(C2_REGISTRATION, partition="7d's registration, checked by 7d's own check"),
        mf.input_entry(C2_CENSUS, partition="hashed by 7d's registration check"),
        mf.input_entry(
            C2_RESULT, partition="the comparator: 7d's committed result, read after the RNA was opened"
        ),
        mf.input_entry(M2B_RESULT, partition="hashed by 7d's registration check"),
        mf.input_entry(M2E_RESULT, partition="context: M2e's census sentence, copied verbatim"),
        mf.input_entry(rg.RNA_SOURCES, partition=None),
    ]
    opened = (
        "measured outcome, values opened after the harness checks"
        if rna_read
        else "measured outcome, hashed by the pin check; no value opened"
    )
    inputs += [
        mf.input_entry(p, partition=opened)
        for p in reg["pinned_sha256"]
        if p.startswith(rg.RNA_DIR.as_posix())
    ]
    return {
        "sources": [
            {"accession": "the 24 compiled non-coding programs 7d ran", "version": "pinned by sha256"},
            {"accession": "reader v1 ENCODE DNase-seq narrowPeak files, 13 biosamples", "version": "pinned"},
            {"accession": "ENCODE polyA plus RNA-seq gene quantifications, V29", "version": "pinned"},
            {"accession": "data/results/gene_rates.json (mouse NIH 3T3, e921fae)", "version": "pinned"},
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38; RNA genes matched by gene_id through data/cache/gencode_genes.tsv",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "form": FORM,
            "comparator_forms": list(mr.COMPARATOR_FORMS),
            "theta": c2.THETA,
            "selection": mr.SELECTION,
            "runtime": mr.RUNTIME,
            "min_genes": mr.MIN_GENES,
            "pseudocount": mr.PSEUDOCOUNT,
            "bootstrap": mr.BOOTSTRAP,
            "head_to_head_interval": list(mr.HEAD_TO_HEAD_INTERVAL),
            "permutations": mr.PERMUTATIONS,
            "control_alpha": mr.CONTROL_ALPHA,
            "comparator_tolerance": mr.COMPARATOR_TOLERANCE,
        },
        "exclusions": [
            "genes with one rule, rules in several programs, or a measured rule",
            "observations unreadable or exactly 1; no activating rule",
            "multiplicative: an inhibitor with RHO > 2 or summed log-folds beyond EXPONENT_LIMIT; a"
            " prediction the same in all thirteen cells",
            "genes whose symbol is not one gene_id, or whose gene_id is not on one row of every table",
        ],
        "partitions": "n/a: every qualifying gene is scored",
        "code_cleanliness": mf.code_cleanliness(__file__, mcensus.OWN_CODE),
    }


def _table(rows: list[dict]) -> dict[str, Any]:
    cols = ("gene", "symbol", "gene_id", "chrom", "n_rules", "n_inhibiting", "P", "D")
    out = [
        [*(r[c] for c in cols[:6]), [round(x, 5) for x in r["P"]], [round(x, 5) for x in r["D"]]]
        for r in rows
    ]
    return {"columns": list(cols), "rows": out}


def descriptives(rows: dict[str, list[dict]], new_genes: set[str], h2h: dict[str, Any]) -> dict[str, Any]:
    """DESCRIPTIVES 2 and 3: the multiplicative form on the genes new under it and on 7d's mean genes;
    each form on the genes in all three populations. Not verdicts."""
    by = {f: {r["gene"]: r for r in rs} for f, rs in rows.items()}
    m = by[FORM]
    on_mean = sorted(set(m) & set(by["mean"]))
    new = [m[g] for g in sorted(new_genes & set(m))]
    common = sorted(set(m) & set(by["mean"]) & set(by["competition"]))

    def rng_q(rs: list[dict]) -> dict[str, float]:
        return mr.quantiles([max(r["P"]) - min(r["P"]) for r in rs])

    return {
        "new_under_the_multiplicative_form": {
            "genes": len(new),
            "multiplicative_delta_mse": c2.dmse(new) if new else None,
            "predicted_range_log2_quantiles": rng_q(new),
        },
        "on_7d_mean_population": {
            "genes": len(on_mean),
            "multiplicative_delta_mse": c2.dmse([m[g] for g in on_mean]) if on_mean else None,
            "mean_delta_mse": c2.dmse([by["mean"][g] for g in on_mean]) if on_mean else None,
        },
        "on_the_genes_in_all_three": {
            "genes": len(common),
            **{
                f: {
                    "delta_mse": c2.dmse([by[f][g] for g in common]) if common else None,
                    "predicted_range_log2_quantiles": rng_q([by[f][g] for g in common]),
                }
                for f in (FORM, *mr.COMPARATOR_FORMS)
            },
        },
        "head_to_heads_compared": h2h.get("compared", False),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    t0 = time.time()
    reg = check_registration()
    c2run.check_registration()
    genes7, why7, programs7 = c2census.collect(CHROMS, run=True, workers=args.workers)
    print("7d collected", round(time.time() - t0, 1), flush=True)
    genes_m, why_m, programs_m = mcensus.collect(CHROMS, run=True, workers=args.workers)
    print("multiplicative collected", round(time.time() - t0, 1), flush=True)
    summary = mcensus.summarise(genes_m, why_m)
    hc = harness_check(genes7, genes_m)
    pre = predictor_problems(genes7, genes_m, programs7, programs_m, summary)
    base = {
        "lane": mr.LANE,
        "registration": "data/results/multrna_registration.json",
        "registration_date": reg.get("date"),
        "question": mr.QUESTION,
        "run_counts": dict(why_m),
        "programs": programs_m,
        "harness_check": hc,
        "predictor_side_problems": pre,
        "census_sentence": mcensus.sentence(summary),
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
    }
    if hc["failing"] or pre:
        kind = "HARNESS_DEFECT" if hc["failing"] else "PREDICTOR_DEFECT"
        nums = (
            {"bad": hc["failing"], "n": hc["gene_form_runs"]} if hc["failing"] else {"what": "; ".join(pre)}
        )
        out = {
            **base,
            "outcome": kind,
            "reading": mr.reading(kind, **nums),
            "seconds": round(time.time() - t0, 1),
        }
        print("saved", save_result(RESULT, out, manifest=manifest(rna_read=False)))
        print(out["reading"], pre)
        return 1
    tpm, rna_counts = load_rna()
    ids = rg.symbol_gene_ids()
    rows: dict[str, list[dict]] = {}
    excl: dict[str, dict] = {}
    for f in mr.COMPARATOR_FORMS:
        rows[f], why_f = c2.assemble(c2run.form_runs(genes7, f), tpm, ids)
        excl[f] = dict(why_f)
    rows[FORM], why_f = mr.assemble(mult_runs(genes_m), tpm, ids)
    excl[FORM] = dict(why_f)
    per7 = {f: c2.score_form(rows[f], f) for f in mr.COMPARATOR_FORMS}
    h7 = c2.head_to_head(
        {r["gene"]: r for r in rows["mean"]},
        {r["gene"]: r for r in rows["competition"]},
        c2.BOOTSTRAP["resamples"],
        c2.BOOTSTRAP["seed"],
    )
    reading7 = c2.combined_reading(per7, h7)
    committed = json.loads(C2_RESULT.read_text())
    problems = mr.comparator_problems({f: rows[f] for f in mr.COMPARATOR_FORMS}, committed, reading7)
    by = {f: {r["gene"]: r for r in rs} for f, rs in rows.items()}
    common = set(by[FORM]) & set(by["mean"]) & set(by["competition"])
    moved_d = sorted(g for g in common if any(by[f][g]["D"] != by[FORM][g]["D"] for f in mr.COMPARATOR_FORMS))
    if moved_d:
        problems.append(f"{len(moved_d)} genes in all three populations carry different D")
    if problems:
        out = {
            **base,
            "outcome": "COMPARATOR_DEFECT",
            "reading": mr.reading("COMPARATOR_DEFECT", what="; ".join(problems)),
            "population_exclusions": excl,
            "seconds": round(time.time() - t0, 1),
        }
        print("saved", save_result(RESULT, out, manifest=manifest(rna_read=True)))
        print(out["reading"])
        return 1
    print("comparators reproduced", round(time.time() - t0, 1), flush=True)
    mult = mr.score_form(rows[FORM])
    h2h = mr.head_to_heads(by, mr.BOOTSTRAP["resamples"], mr.BOOTSTRAP["seed"])
    reading = mr.combined_reading(mult, h2h)
    new_genes = {
        g["gene"]
        for g in genes_m
        if mcensus.informative(g, FORM) and not any(mcensus.informative(g, f) for f in c2.FORMS)
    }
    m2e = json.loads(M2E_RESULT.read_text())
    out = {
        **base,
        "reading": reading,
        "outcome": mult["outcome"],
        "head_to_head_clauses": {f: h2h[f]["kind"] for f in mr.COMPARATOR_FORMS}
        if h2h.get("compared")
        else None,
        "multiplicative": mult,
        "head_to_heads": h2h,
        "comparators_reproduced": {
            "source": "data/results/combine2_result.json",
            "populations": {f: len(rows[f]) for f in mr.COMPARATOR_FORMS},
            "reading_recomputed_equals_committed": True,
            "comparator_tolerance": mr.COMPARATOR_TOLERANCE,
        },
        "descriptives_not_verdicts": descriptives(rows, new_genes, h2h),
        "limitations": list(mr.LIMITATIONS),
        "context_verbatim": {
            "seven_d_reading": committed["reading"],
            "m2e_census_sentence": m2e["census"]["sentence"],
        },
        "rna_tables": rna_counts,
        "population_exclusions": excl,
        "cells": list(mr.CELLS),
        "genes_table": _table(rows[FORM]),
        "seconds": round(time.time() - t0, 1),
    }
    print("saved", save_result(RESULT, out, manifest=manifest(rna_read=True), compact=True))
    print(reading)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
