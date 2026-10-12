# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-residuemult: how many of the 3,708 the mean form cannot fit are already fitted multiplicatively.

    uv run python scripts/residuemult_lookup.py [--chroms chr21] [--no-save] [--workers 4]

A LOOKUP. The definitions, the identity convention, the two counts, the bound checks and the implication
rule with its threshold were all committed first, in `genomeos/attribution/residuemult.py`. This script
re-derives the residue's gene SYMBOLS with `lane-residue`'s own code path (its committed per-gene table
carries no gene identity), matches them against `lane-multrna`'s committed multiplicative table, and
applies 7h's own runnability and informativeness test to the same genes with M2h's fixed scale on the
observations.

NO RNA VALUE IS OPENED. Of `multrna_result.json`'s per-gene table only the `symbol` column is read; its
`P` and `D` columns are never touched. NOTHING IS FITTED, TUNED, SELECTED OR RE-ESTIMATED: the scale is a
constant from `ac5580b` and the multiplicative form has no parameter of its own. No model request is sent
and no key is read. This lane registers no test and makes no claim about expression.
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
from genomeos.attribution import residuemult as rm  # noqa: E402
from genomeos.attribution import scalefit as sf  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.results import save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readergate_census import peak_files  # noqa: E402

RESULT = "residuemult_lookup"
MULTRNA_RESULT = Path("data/results/multrna_result.json")
MULTRNA_CENSUS = Path("data/results/multrna_census.json")
RESIDUE_CENSUS = Path("data/results/residue_census.json")
RESIDUE_GENES = Path("data/results/residue_genes.json")
MAGCAL = Path("data/results/magcal_result.json")
CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]
FORM = "mean"  # the form whose residue this is; the multiplicative one is saturate.fit
OWN_CODE = (
    "genomeos/attribution/residuemult.py",
    "scripts/residuemult_lookup.py",
    "tests/test_residuemult.py",
)

COLUMNS = (
    "symbol",
    "chrom",
    "n_rules",
    "n_inhibiting",
    "in_committed_table",
    "runnable_unscaled",
    "informative_unscaled",
    "runnable_scaled",
    "informative_scaled",
    "reason_unscaled",
    "reason_scaled",
)


def committed_symbols() -> set[str]:
    """The `symbol` column of `lane-multrna`'s committed per-gene table, AND NOTHING ELSE OF IT.

    The table's `P` and `D` columns are predictions and measured RNA. They are not returned, not kept and
    not looked at: this function reads one column out of the rows and discards the rest.
    """
    d = json.loads(MULTRNA_RESULT.read_text())
    table = d["genes_table"]
    i = table["columns"].index("symbol")
    return {row[i] for row in table["rows"]}


def _chromosome(args: tuple[str, set[str], set[str]]) -> dict[str, Any]:
    """One program: `lane-residue`'s re-derivation of the residue, with each gene's symbol and its
    multiplicative lookup. Nothing is fitted; `saturate.fit` has no parameter and SCALE is a constant."""
    from genomeos.lang.parser import parse

    chrom, genes, committed = args
    t0 = time.time()
    text = (Path(COMPILED_DIR) / f"noncoding_{chrom}.bio").read_text()
    gated, info = cp.gate_program(text, chrom, selection=c2.SELECTION)
    module = parse(gated, f"noncoding_{chrom}_reader")
    out: dict[str, Any] = {
        "chrom": chrom,
        "program": {k: info[k] for k in ("base_sha256", "gated_sha256", "elements", "rules_gated")},
        "eligible": 0,
        "excluded": 0,
        "recovered": 0,
        "rows": [],
        "symbols": [],
        "violations": [],
        "max_abs_q_after_minus_scale_times_q_before": 0.0,
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
        # BOUND_CHECKS (e): the identity lane-scalefit proved, re-checked on every gene
        out["max_abs_q_after_minus_scale_times_q_before"] = max(
            out["max_abs_q_after_minus_scale_times_q_before"], abs(after.q - rm.SCALE * before.q)
        )
        if after.q < 1.0:
            out["recovered"] += 1
            continue
        # BOUND_CHECKS (d): the residue's own definition, re-asserted per gene
        if after.ok or after.reason != "basal_not_positive":
            out["violations"].append(
                f"{gr.symbol}: in the residue by q but the scaled mean fit reads {after.reason!r}"
            )
        row = rm.lookup_gene(gr, gr.symbol in committed)
        out["violations"].extend(rm.bound_violations(row))
        out["symbols"].append(gr.symbol)
        out["rows"].append(
            [
                row.symbol,
                chrom,
                row.n_rules,
                row.n_inhibiting,
                int(row.in_committed_table),
                int(row.runnable_unscaled),
                int(row.informative_unscaled),
                int(row.runnable_scaled),
                int(row.informative_scaled),
                row.reason_unscaled,
                row.reason_scaled,
            ]
        )
    out["seconds"] = round(time.time() - t0, 1)
    return out


def reconcile(eligible: int, excluded: int, recovered: int, residue: int) -> dict[str, Any]:
    """This run's re-derived counts against `lane-residue`'s committed ones. A disagreement is reported."""
    d = json.loads(RESIDUE_CENSUS.read_text())["Q1_counts"]
    got = {
        "eligible_for_a_fit": eligible,
        "excluded_basal_not_positive": excluded,
        "recovered_at_the_scale": recovered,
        "the_residue": residue,
    }
    return {
        "source": str(RESIDUE_CENSUS),
        "committed": d,
        "re_derived": got,
        "agrees": got == d,
        "differences": {k: [d[k], got[k]] for k in d if d[k] != got[k]},
    }


def lookup(chroms: list[str], workers: int) -> dict[str, Any]:
    """The whole lookup: re-derive the residue's identities, match, test at the scale, decide."""
    texts = {c: (Path(COMPILED_DIR) / f"noncoding_{c}.bio").read_text() for c in chroms}
    chrom_of, _why = c2.multi_rule_genes(texts)
    del texts
    by_chrom: dict[str, set[str]] = {c: set() for c in chroms}
    for g, c in chrom_of.items():
        if c in by_chrom:
            by_chrom[c].add(g)
    committed = committed_symbols()
    jobs = [(c, by_chrom[c], committed) for c in chroms]
    ctx = get_context("fork")
    with ctx.Pool(workers) as pool:
        parts = pool.map(_chromosome, jobs)

    rows: list[list[Any]] = []
    symbols: list[str] = []
    violations: list[str] = []
    programs: dict[str, Any] = {}
    eligible = excluded = recovered = 0
    max_identity = 0.0
    for p in parts:
        rows.extend(p["rows"])
        symbols.extend(p["symbols"])
        violations.extend(p["violations"])
        programs[p["chrom"]] = p["program"]
        eligible += p["eligible"]
        excluded += p["excluded"]
        recovered += p["recovered"]
        max_identity = max(max_identity, p["max_abs_q_after_minus_scale_times_q_before"])

    i = {c: COLUMNS.index(c) for c in COLUMNS}
    n = len(rows)
    in_table = [r for r in rows if r[i["in_committed_table"]]]
    inf_unscaled = [r for r in rows if r[i["informative_unscaled"]]]
    inf_scaled = [r for r in rows if r[i["informative_scaled"]]]
    union = [r for r in rows if r[i["in_committed_table"]] or r[i["informative_scaled"]]]
    neither = [r for r in rows if not (r[i["in_committed_table"]] or r[i["informative_scaled"]])]
    branch = rm.implication(len(union), n)
    dupes = sorted(s for s, k in Counter(symbols).items() if k > 1)

    return {
        "result": RESULT,
        "lane": rm.LANE,
        "question": rm.QUESTION,
        "identity": rm.IDENTITY,
        "match_field": "symbol",
        "re_derived": rm.RE_DERIVED,
        "count_definitions": {
            "count_1_committed_table": rm.COUNT_COMMITTED,
            "count_2_at_the_scale": rm.COUNT_SCALED,
            "complement": rm.COMPLEMENT,
        },
        "residue_genes": n,
        "counts": {
            "count_1_in_multrna_s_committed_table": len(in_table),
            "count_1b_informative_multiplicatively_predictor_side": len(inf_unscaled),
            "count_2_informative_multiplicatively_at_the_scale": len(inf_scaled),
            "union_of_count_1_and_count_2": len(union),
            "informative_under_neither": len(neither),
            "runnable_unscaled": sum(1 for r in rows if r[i["runnable_unscaled"]]),
            "runnable_at_the_scale": sum(1 for r in rows if r[i["runnable_scaled"]]),
        },
        "shares": {
            "count_1": len(in_table) / n if n else float("nan"),
            "count_2": len(inf_scaled) / n if n else float("nan"),
            "union": len(union) / n if n else float("nan"),
            "neither": len(neither) / n if n else float("nan"),
        },
        "reasons_not_runnable": {
            "unscaled": dict(Counter(r[i["reason_unscaled"]] for r in rows if not r[i["runnable_unscaled"]])),
            "at_the_scale": dict(Counter(r[i["reason_scaled"]] for r in rows if not r[i["runnable_scaled"]])),
        },
        "informative_under_neither_by_reason": dict(
            Counter(
                r[i["reason_scaled"]] if not r[i["runnable_scaled"]] else "runnable_not_informative"
                for r in neither
            )
        ),
        "in_the_committed_table_but_not_informative_at_the_scale": sum(
            1 for r in in_table if not r[i["informative_scaled"]]
        ),
        "informative_at_the_scale_but_not_in_the_committed_table": sum(
            1 for r in inf_scaled if not r[i["in_committed_table"]]
        ),
        "implication_rule": rm.IMPLICATION_RULE,
        "most": rm.MOST,
        "branch": branch,
        "implication": rm.reading(
            branch,
            union=len(union),
            residue=n,
            share=f"{100 * len(union) / n:.1f}%" if n else "n/a",
            committed=len(in_table),
            scaled=len(inf_scaled),
            complement=len(neither),
        ),
        "why_the_scale_can_only_help": rm.WHY_THE_SCALE_CAN_ONLY_HELP_THE_MULTIPLICATIVE_FORM,
        "bound_checks": rm.BOUND_CHECKS,
        "bound_violations": violations,
        "identity_checks": {
            "max_abs_q_after_minus_scale_times_q_before": max_identity,
            "residue_symbols_not_unique_in_this_run": dupes,
        },
        "reconciliation": reconcile(eligible, excluded, recovered, n),
        "committed_tables": rm.COMMITTED,
        "multrna_committed_table_rows": len(committed),
        "multrna_census_sentence": json.loads(MULTRNA_CENSUS.read_text())["sentence"],
        "what_this_is_not": (
            "not a test, not a registration and not a claim about expression. Every count here is whether"
            " a fit accepts a gene, which is arithmetic on the elements' own removal observations."
            " lane-multrna's outcome on measured RNA was WORSE and that is carried forward unchanged: the"
            " multiplicative form being the one that CAN RUN these genes is not the multiplicative form"
            " reading closer to measured RNA on them. NO RNA VALUE WAS OPENED IN THIS LANE"
        ),
        "prior_exposure": rm.PRIOR_EXPOSURE,
        "limitations": list(rm.LIMITATIONS),
        "conventions": rm.CONVENTIONS,
        "programs": programs,
        "per_gene_table": {"columns": list(COLUMNS), "rows": rows},
        "rna_read": "none",
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
    }


def manifest(chroms: list[str]) -> dict[str, Any]:
    inputs = [
        *(mf.input_entry(Path(COMPILED_DIR) / f"noncoding_{c}.bio", partition=None) for c in chroms),
        mf.input_entry(MULTRNA_RESULT, partition=None),
        mf.input_entry(MULTRNA_CENSUS, partition=None),
        mf.input_entry(RESIDUE_CENSUS, partition=None),
        mf.input_entry(RESIDUE_GENES, partition=None),
        mf.input_entry(MAGCAL, partition=None),
        mf.input_entry(rg.GENE_TABLE, partition=None),
        mf.files_entry("reader v1 DNase narrowPeak files", [p for p in peak_files(chroms) if p.exists()]),
    ]
    return {
        "sources": [
            {
                "accession": "the 24 compiled non-coding programs 7d, 7h and lane-residue ran",
                "version": "pinned by sha256",
            },
            {"accession": "data/results/multrna_result.json (7h, the committed table)", "version": "dev"},
            {"accession": "data/results/multrna_census.json (7h)", "version": "dev"},
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
            "form_looked_up": "multiplicative",
            "selection": c2.SELECTION,
            **rm.CONVENTIONS,
        },
        "exclusions": [
            "7d's own: one rule, rules in several programs, a measured rule, an unreadable or exactly-1"
            " observation, an action that disagrees with it, no activating rule",
            "genes the mean form excluded for a reason other than basal_not_positive are not in this"
            " lookup's population at all",
            "genes the mean form recovered at M2h's scale are counted and then dropped: the population is"
            " the residue alone",
            "COUNT 1 inherits 7h's table's own RNA-row exclusions (139 + 30); COUNT 1b is the same"
            " question without them",
        ],
        "partitions": "n/a: a predictor-side lookup; no outcome read and no RNA opened",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", nargs="*", default=CHROMS)
    ap.add_argument("--no-save", action="store_true")
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    t0 = time.time()
    out = lookup(args.chroms, args.workers)
    out["seconds"] = round(time.time() - t0, 1)
    print("residue genes (re-derived):", out["residue_genes"], flush=True)
    for k, v in out["counts"].items():
        print(f"  {k}: {v}", flush=True)
    print("shares", {k: round(v, 4) for k, v in out["shares"].items()}, flush=True)
    print("not runnable, by reason:", out["reasons_not_runnable"], flush=True)
    print("neither, by reason:", out["informative_under_neither_by_reason"], flush=True)
    print("reconciliation agrees:", out["reconciliation"]["agrees"], out["reconciliation"]["differences"])
    print("identity checks:", out["identity_checks"]["max_abs_q_after_minus_scale_times_q_before"])
    print("duplicate symbols:", len(out["identity_checks"]["residue_symbols_not_unique_in_this_run"]))
    print("bound violations:", len(out["bound_violations"]), out["bound_violations"][:5], flush=True)
    print("BRANCH", out["branch"], flush=True)
    print("IMPLICATION", out["implication"], flush=True)
    if not args.no_save:
        print("saved", save_result(RESULT, out, manifest=manifest(args.chroms), compact=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
