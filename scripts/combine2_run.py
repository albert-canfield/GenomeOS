# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-combine2, the registered run: multi-element genes under each combination form, against RNA.

    uv run python scripts/combine2_run.py [--workers 4]

Runs only from the committed registration (data/results/combine2_registration.json); refuses if any
registered field differs from the constants in genomeos/attribution/combine2.py, or any pinned file's
sha256 differs from the file on disk.

Order, as registered (combine2.DECISION):
  1. each of the 24 programs is gated (compile.gate_program, 'reader_cells') and parsed; every
     multi-element gene is classified and fitted under each form, and every runnable, informative one
     is RUN by `runtime.grn.NetworkRuntime` in all thirteen reader cells (scripts/combine2_census.py's
     `collect`, the code the census ran);
  2. the HARNESS is checked on every gene-form run; any miss stops the run here, before the RNA opens;
  3. only then are the RNA quantifications opened, each form's population assembled and scored, and
     the head-to-head read.

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
from genomeos.attribution import readergate as rg  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from combine2_census import CHROMS, OWN_CODE, collect  # noqa: E402
from readergate_census import peak_files  # noqa: E402
from readergate_run import load_rna  # noqa: E402
from runrules_census import RATE_TABLE  # noqa: E402

RESULT = "combine2_result"
REGISTRATION = RESULTS_DIR / "combine2_registration.json"
CENSUS = RESULTS_DIR / "combine2_census.json"
M2B_RESULT = RESULTS_DIR / "readergate_result.json"
#: registration field -> the combine2 constant it must equal, or the run refuses
BOUND = {
    "question": "QUESTION",
    "forms": "FORMS",
    "theta": "THETA",
    "form_sources": "FORM_SOURCES",
    "fit": "FIT",
    "structure": "STRUCTURE",
    "selection": "SELECTION",
    "cells": "CELLS",
    "rna_selection": "RNA_SELECTION",
    "population": "POPULATION",
    "informative_log2": "INFORMATIVE_LOG2",
    "min_genes": "MIN_GENES",
    "pseudocount": "PSEUDOCOUNT",
    "measured": "MEASURED",
    "predictions": "PREDICTIONS",
    "metric": "METRIC",
    "form_interval": "FORM_INTERVAL",
    "head_to_head_interval": "HEAD_TO_HEAD_INTERVAL",
    "bootstrap": "BOOTSTRAP",
    "control": "CONTROL",
    "permutations": "PERMUTATIONS",
    "control_alpha": "CONTROL_ALPHA",
    "runtime": "RUNTIME",
    "closed_form_tolerance": "CLOSED_FORM_TOLERANCE",
    "harness": "HARNESS",
    "decision": "DECISION",
    "falsifier": "FALSIFIER",
    "stop_rule": "STOP_RULE",
    "form_words": "FORM_WORDS",
    "outcomes": "OUTCOMES",
    "head_to_head": "HEAD_TO_HEAD",
    "limitations": "LIMITATIONS",
    "descriptives_not_verdicts": "DESCRIPTIVES",
}


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def check_registration() -> dict:
    """The committed registration, refused if a bound field differs from the code or a pin from disk."""
    reg = json.loads(REGISTRATION.read_text())
    bad = [k for k, const in BOUND.items() if reg.get(k) != json.loads(json.dumps(getattr(c2, const)))]
    if bad:
        raise SystemExit(f"refused: the registration and combine2.py disagree on {', '.join(bad)}")
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


def form_runs(genes: list[dict], form: str) -> list[dict[str, Any]]:
    """The records `combine2.assemble` reads, for the genes runnable and informative under `form`."""
    out = []
    for g in genes:
        fr = g["forms"].get(form)
        if g["excluded"] or not fr or not fr["ok"] or not fr["informative"]:
            continue
        out.append(
            {
                **{
                    k: g[k]
                    for k in ("gene", "symbol", "chrom", "n_rules", "n_inhibiting", "elements", "states")
                },
                "open_share": g["open_share"],
                **{k: fr[k] for k in ("weights", "e", "q", "z", "s")},
            }
        )
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
        mf.input_entry(rg.RNA_SOURCES, partition=None),
        mf.input_entry(M2B_RESULT, partition="context: M2b's registered reading, copied verbatim"),
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
            {
                "accession": "the 24 compiled non-coding programs M1 and M2b ran",
                "version": "pinned by sha256",
            },
            {"accession": "reader v1 ENCODE DNase-seq narrowPeak files, 13 biosamples", "version": "pinned"},
            {"accession": "ENCODE polyA plus RNA-seq gene quantifications, V29", "version": "pinned"},
            {"accession": "data/results/gene_rates.json (mouse NIH 3T3, e921fae)", "version": "pinned"},
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38; RNA genes matched by gene_id through data/cache/gencode_genes.tsv",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "forms": list(c2.FORMS),
            "theta": c2.THETA,
            "selection": c2.SELECTION,
            "runtime": c2.RUNTIME,
            "min_genes": c2.MIN_GENES,
            "pseudocount": c2.PSEUDOCOUNT,
            "bootstrap": c2.BOOTSTRAP,
            "permutations": c2.PERMUTATIONS,
            "control_alpha": c2.CONTROL_ALPHA,
        },
        "exclusions": [
            "genes with one rule, rules in several programs, or a measured rule",
            "observations unreadable or exactly 1; no activating rule",
            "per form: a fit leaving no positive basal rate or an inhibitor s > 1; a prediction the same in"
            " all thirteen cells",
            "genes whose symbol is not one gene_id, or whose gene_id is not on one row of every table",
        ],
        "partitions": "n/a: every qualifying gene is scored",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def _table(rows: list[dict]) -> dict[str, Any]:
    cols = ("gene", "symbol", "gene_id", "chrom", "n_rules", "n_inhibiting", "q", "P", "D")
    out = [
        [
            *(r[c] for c in cols[:6]),
            round(r["q"], 6),
            [round(x, 5) for x in r["P"]],
            [round(x, 5) for x in r["D"]],
        ]
        for r in rows
    ]
    return {"columns": list(cols), "rows": out}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    t0 = time.time()
    reg = check_registration()
    genes, why, programs = collect(CHROMS, run=True, workers=args.workers)
    ran = [(g, f) for g in genes for f, fr in g["forms"].items() if "s" in fr]
    bad = [(g["gene"], f) for g, f in ran if g["forms"][f]["harness_miss"]]
    base = {
        "lane": c2.LANE,
        "registration": "data/results/combine2_registration.json",
        "registration_date": reg.get("date"),
        "question": c2.QUESTION,
        "run_counts": dict(why),
        "programs": programs,
        "harness_check": {
            "gene_form_runs": len(ran),
            "failing": len(bad),
            "max_abs_run_minus_closed_form_log2": max(
                (
                    abs(a - b)
                    for g, f in ran
                    for a, b in zip(g["forms"][f]["s"], g["forms"][f]["closed_s"], strict=True)
                ),
                default=None,
            ),
        },
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
    }
    if bad:
        out = {
            **base,
            "outcome": "HARNESS_DEFECT",
            "reading": c2.reading("HARNESS_DEFECT", bad=len(bad), n=len(ran)),
            "failing_examples": bad[:20],
            "seconds": round(time.time() - t0, 1),
        }
        print("saved", save_result(RESULT, out, manifest=manifest(rna_read=False)))
        print(out["reading"])
        return 1
    tpm, rna_counts = load_rna()
    ids = rg.symbol_gene_ids()
    per_form, rows_by_form, excl = {}, {}, {}
    for form in c2.FORMS:
        rows, why_f = c2.assemble(form_runs(genes, form), tpm, ids)
        rows_by_form[form] = rows
        excl[form] = dict(why_f)
        per_form[form] = c2.score_form(rows, form)
        print(form, per_form[form]["outcome"], flush=True)
    h2h = c2.head_to_head(
        {r["gene"]: r for r in rows_by_form["mean"]},
        {r["gene"]: r for r in rows_by_form["competition"]},
        c2.BOOTSTRAP["resamples"],
        c2.BOOTSTRAP["seed"],
    )
    reading = c2.combined_reading(per_form, h2h)
    m2b = json.loads(M2B_RESULT.read_text())
    out = {
        **base,
        "reading": reading,
        "outcomes": {f: per_form[f]["outcome"] for f in c2.FORMS},
        "per_form": per_form,
        "head_to_head": h2h,
        "limitations": list(c2.LIMITATIONS),
        "context_single_rule_m2b": {
            "source": "data/results/readergate_result.json",
            "outcome": m2b.get("outcome"),
            "reading": m2b.get("reading"),
        },
        "rna_tables": rna_counts,
        "population_exclusions": excl,
        "cells": list(c2.CELLS),
        "genes_tables": {f: _table(rows_by_form[f]) for f in c2.FORMS},
        "seconds": round(time.time() - t0, 1),
    }
    print("saved", save_result(RESULT, out, manifest=manifest(rna_read=True)))
    print(reading)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
