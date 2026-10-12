# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-recovered step 1: re-derive the genes scaling makes fittable, and measure the overlap with M2h's.

    uv run python scripts/recovered_population.py [--chroms chr21] [--no-save] [--workers 4]

THIS STEP OPENS NO RNA VALUE. It recomputes `combine2`'s eligibility and `scalefit`'s recovery over the
same 24 compiled non-coding programs 7d and lane-scalefit ran, so the 12,500 / 8,792 / 8,535 are
RE-DERIVED here and not copied; `scalefit.SCALE` is M2h's constant and nothing is fitted. From 7h's
committed multiplicative genes table it reads THE SYMBOL COLUMN ONLY, to find which recovered genes have
a P and a D at all and to count the overlap with M2h's own training and held-out halves. No P and no D
value is read, summed, scaled or compared in this script.

It runs BEFORE the registration, because the registration's exposure sentence has to say what the
overlap actually is (`recovered.OVERLAP_QUESTION`) instead of assuming it.
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
from genomeos.attribution import recovered as rc  # noqa: E402
from genomeos.attribution import scalefit as sf  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.results import save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readergate_census import peak_files  # noqa: E402

RESULT = "recovered_population"
SCALEFIT_CENSUS = Path("data/results/scalefit_census.json")
CENSUS_7D = Path("data/results/combine2_census.json")
MULTRNA = Path("data/results/multrna_result.json")
MAGCAL = Path("data/results/magcal_result.json")
CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]
OWN_CODE = (
    "genomeos/attribution/recovered.py",
    "scripts/recovered_population.py",
    "scripts/recovered_register.py",
    "scripts/recovered_run.py",
    "tests/test_recovered.py",
)


def texts_of(chroms: list[str]) -> dict[str, str]:
    return {c: (Path(COMPILED_DIR) / f"noncoding_{c}.bio").read_text() for c in chroms}


def _chromosome(args: tuple[str, set[str]]) -> dict[str, Any]:
    """One program's multi-element genes, classified by 7d's code then by scalefit's recovery."""
    from genomeos.lang.parser import parse

    chrom, genes = args
    t0 = time.time()
    text = (Path(COMPILED_DIR) / f"noncoding_{chrom}.bio").read_text()
    gated, info = cp.gate_program(text, chrom, selection=c2.SELECTION)
    module = parse(gated, f"noncoding_{chrom}_reader")
    out: dict[str, Any] = {
        "chrom": chrom,
        "program": {k: info[k] for k in ("base_sha256", "gated_sha256", "elements", "rules_gated")},
        "eligible": [],
        "excluded": [],
        "recovered": [],
        "runnable": [],
        "n_rules": {},
    }
    for gr in c2.gene_rules(module, genes, chrom, rg.compiled_cells(text)):
        _fits, why = c2.classify(gr)
        if why:
            continue
        out["eligible"].append(gr.symbol)
        excluded, recovered, runnable = rc.is_recovered(gr.rho)
        if excluded:
            out["excluded"].append(gr.symbol)
        if recovered:
            out["recovered"].append(gr.symbol)
            out["n_rules"][gr.symbol] = len(gr.elements)
        if runnable:
            out["runnable"].append(gr.symbol)
    out["seconds"] = round(time.time() - t0, 1)
    return out


def collect(chroms: list[str], workers: int) -> tuple[dict[str, list[str]], dict, dict[str, int]]:
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
                print(r["chrom"], len(r["eligible"]), len(r["runnable"]), r["seconds"], flush=True)
    else:
        for j in jobs:
            r = _chromosome(j)
            res[r["chrom"]] = r
            print(r["chrom"], len(r["eligible"]), len(r["runnable"]), r["seconds"], flush=True)
    sets = {
        key: [s for c in chroms for s in res[c][key]]
        for key in ("eligible", "excluded", "recovered", "runnable")
    }
    rules = {s: n for c in chroms for s, n in res[c]["n_rules"].items()}
    return sets, {c: res[c]["program"] for c in chroms}, rules


def multrna_symbols() -> set[str]:
    """7h's multiplicative population, SYMBOL COLUMN ONLY. No P and no D value is read."""
    table = json.loads(MULTRNA.read_text())["genes_table"]
    i = table["columns"].index("symbol")
    return {row[i] for row in table["rows"]}


def harness(sets: dict[str, list[str]]) -> dict[str, Any]:
    """POPULATION_HARNESS, run before anything else is computed."""
    want = rc.POPULATION_SOURCE_COUNTS
    excl, rec, run = set(sets["excluded"]), set(sets["recovered"]), set(sets["runnable"])
    checks = {
        "eligible_for_a_fit": (len(sets["eligible"]), want["eligible_for_a_fit"]),
        "excluded_basal_not_positive": (len(sets["excluded"]), want["excluded_basal_not_positive"]),
        "recovered": (len(sets["recovered"]), want["recovered"]),
        "runnable_after": (len(sets["runnable"]), want["runnable_after"]),
    }
    out: dict[str, Any] = {
        "counts": {k: {"here": a, "scalefit": b, "same": a == b} for k, (a, b) in checks.items()},
        "recovered_is_a_subset_of_excluded": rec <= excl,
        "runnable_is_a_subset_of_recovered": run <= rec,
        "every_symbol_non_empty": all(s for s in sets["eligible"]),
        "no_symbol_twice": all(len(sets[k]) == len(set(sets[k])) for k in sets),
    }
    out["agrees"] = all(v["same"] for v in out["counts"].values()) and all(
        out[k]
        for k in (
            "recovered_is_a_subset_of_excluded",
            "runnable_is_a_subset_of_recovered",
            "every_symbol_non_empty",
            "no_symbol_twice",
        )
    )
    return out


def population(chroms: list[str], workers: int) -> dict:
    sets, programs, rules = collect(chroms, workers)
    hc = harness(sets)
    if not hc["agrees"]:
        return {
            "lane": rc.LANE,
            "date": time.strftime("%Y-%m-%d"),
            "outcome": "POPULATION_DEFECT",
            "population_harness": hc,
            "note": "the re-derivation disagrees with scalefit's committed counts; nothing else computed",
            "rna_read": "none",
            "alphagenome_requests": 0,
            "money": "none: no request is sent and no key is read",
        }

    run = set(sets["runnable"])
    theirs = multrna_symbols()
    with_rows = sorted(run & theirs)
    tr_mine, ho_mine = rc.held_out_split(with_rows)
    tr_theirs, ho_theirs = rc.held_out_split(sorted(theirs))

    ov_all = rc.overlap(run, theirs)
    ov_train = rc.overlap(run, set(tr_theirs))
    ov_held = rc.overlap(run, set(ho_theirs))
    mean_run = json.loads(CENSUS_7D.read_text())["reader_selection"]["forms"]["mean"]["run"]

    return {
        "lane": rc.LANE,
        "date": time.strftime("%Y-%m-%d"),
        "step": "1 of 3: the population, established BEFORE the registration and before any RNA value",
        "outcome": "POPULATION_ESTABLISHED",
        "population": rc.POPULATION,
        "population_form": rc.POPULATION_FORM,
        "prediction_form": rc.FORM,
        "scale": sf.SCALE,
        "scale_rounded": sf.SCALE_ROUNDED,
        "scale_is_refitted": False,
        "scale_source": sf.SCALE_SOURCE,
        "population_harness": hc,
        "population_harness_words": rc.POPULATION_HARNESS,
        "re_derived_counts": {
            "eligible_for_a_fit": len(sets["eligible"]),
            "excluded_basal_not_positive": len(sets["excluded"]),
            "recovered": len(sets["recovered"]),
            "runnable_after": len(sets["runnable"]),
            "runnable_with_a_row_in_7h_s_multiplicative_table": len(with_rows),
            "runnable_without_a_row_in_7h_s_multiplicative_table": len(run - theirs),
            "training_of_this_population": len(tr_mine),
            "held_out_of_this_population": len(ho_mine),
        },
        "overlap_question": rc.OVERLAP_QUESTION,
        "overlap_with_m2h": {
            "all_16509": ov_all,
            "m2h_training_8181": ov_train,
            "m2h_held_out_8328": ov_held,
            "sentences": [
                rc.overlap_sentence(ov_all, "M2h's multiplicative population"),
                rc.overlap_sentence(ov_train, "M2h's training genes"),
                rc.overlap_sentence(ov_held, "M2h's held-out genes"),
            ],
        },
        "overlap_with_7d_s_mean_population": {
            "7d_mean_run_genes_recorded_in_its_census": mean_run,
            "note": (
                "disjoint BY CONSTRUCTION: this population is exactly the genes 7d's mean form EXCLUDED,"
                " so no gene 7d's mean reading ran on can be in it. Stated as a construction, not as a"
                " count, because 7d's census records its run population as a count and not as symbols"
            ),
        },
        "population_sha256": rc.sha256_of_symbols(with_rows),
        "held_out_sha256": rc.sha256_of_symbols(ho_mine),
        "symbols": {
            "note": (
                "the population itself, so the registration and the run cannot drift from this step."
                " `all` is the genes of POPULATION (a)-(d); `held_out` is the half the endpoint is read"
                " on under magcal.is_held_out"
            ),
            "all": with_rows,
            "held_out": ho_mine,
            "training": tr_mine,
            "runnable_without_a_row": sorted(run - theirs),
        },
        "rule_count_quantiles_of_the_population": sf.quantiles(
            [float(rules[s]) for s in with_rows if s in rules]
        ),
        "what_this_is_not": (
            "not a reading. No RNA value is opened, no prediction is evaluated and no comparison with the"
            " no-rule baseline is made here. Those belong to the registration and the one run that"
            " follows it"
        ),
        "conventions": rc.CONVENTIONS,
        "programs": programs,
        "rna_read": "none: 7h's genes table is read for its SYMBOL column only",
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
    }


def manifest(chroms: list[str]) -> dict:
    inputs = [
        *(mf.input_entry(Path(COMPILED_DIR) / f"noncoding_{c}.bio", partition=None) for c in chroms),
        mf.input_entry(SCALEFIT_CENSUS, partition=None),
        mf.input_entry(CENSUS_7D, partition=None),
        mf.input_entry(MULTRNA, partition=None),
        mf.input_entry(MAGCAL, partition=None),
        mf.input_entry(rg.GENE_TABLE, partition=None),
        mf.files_entry("reader v1 DNase narrowPeak files", [p for p in peak_files(chroms) if p.exists()]),
    ]
    return {
        "sources": [
            {
                "accession": "the 24 compiled non-coding programs M1, M2b, 7d and lane-scalefit ran",
                "version": "pinned by sha256",
            },
            {"accession": "data/results/scalefit_census.json (lane-scalefit)", "version": "7f65544"},
            {"accession": "data/results/combine2_census.json (7d/7e)", "version": "948cfc9"},
            {
                "accession": "data/results/multrna_result.json (7h, symbol column only)",
                "version": "pinned by sha256",
            },
            {
                "accession": "data/results/magcal_result.json (M2h, the scale and the split)",
                "version": "ac5580b",
            },
            {"accession": "reader v1 ENCODE DNase-seq narrowPeak files, 13 biosamples", "version": "pinned"},
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {k: v for k, v in rc.CONVENTIONS.items() if k != "split_rule"},
        "exclusions": [
            "7d's own: one rule, rules in several programs, a measured rule, an unreadable or exactly-1"
            " observation, an action that disagrees with it, no activating rule",
            "genes 7d excluded for a reason other than basal_not_positive",
            "recovered genes that do not pass 7d's full runnable test on the scaled observations",
            "recovered-and-runnable genes with no row in 7h's committed multiplicative genes table, which"
            " is what would supply their P and D",
        ],
        "partitions": "n/a: a predictor-side population step; no outcome read and no RNA opened",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", nargs="*", default=CHROMS)
    ap.add_argument("--no-save", action="store_true")
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    t0 = time.time()
    out = population(args.chroms, args.workers)
    out["seconds"] = round(time.time() - t0, 1)
    print("OUTCOME", out["outcome"], flush=True)
    print("harness", json.dumps(out["population_harness"]), flush=True)
    if out["outcome"] == "POPULATION_ESTABLISHED":
        print("counts", json.dumps(out["re_derived_counts"], indent=1), flush=True)
        for s in out["overlap_with_m2h"]["sentences"]:
            print("OVERLAP:", s, flush=True)
    if not args.no_save:
        print("saved", save_result(RESULT, out, manifest=manifest(args.chroms)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
