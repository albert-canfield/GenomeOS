# SPDX-License-Identifier: AGPL-3.0-or-later
"""M2b, the registered run: reader-gated compiled rules against measured RNA in the reader's 13 cells.

    uv run python scripts/readergate_run.py

Runs only from the committed registration (data/results/readergate_registration.json); refuses if any
registered field differs from the constants in genomeos/attribution/readergate.py, or any pinned file's
sha256 differs from the file on disk.

Order, as registered (readergate.DECISION):
  1. each of the 24 programs is gated (compile.gate_program, selection 'reader_cells') and parsed, and
     every single-rule gene is bridged and RUN by `runtime.grn.NetworkRuntime` in all thirteen reader
     cells, the reader in force;
  2. every gene is checked against its closed form and every cell's level against the clamped level
     of its state; any miss stops the run here, before the RNA is opened;
  3. only then are the RNA quantifications opened, the population assembled and the comparison scored.

No model request, no network, no money.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import compile as cp  # noqa: E402
from genomeos.attribution import readergate as rg  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readergate_census import CHROMS, peak_files  # noqa: E402
from readergate_fetch import OWN_CODE  # noqa: E402
from runrules_census import RATE_TABLE, load_rates  # noqa: E402

RESULT = "readergate_result"
REGISTRATION = RESULTS_DIR / "readergate_registration.json"
#: the registration fields that must equal the module's constants, or the run refuses
BOUND = {
    "question": rg.QUESTION,
    "selection": rg.SELECTION,
    "cells": list(rg.CELLS),
    "reader_state": cp.READER_STATE,
    "rna_selection": rg.RNA_SELECTION,
    "population": rg.POPULATION,
    "measured": rg.MEASURED,
    "predictions": rg.PREDICTIONS,
    "metric": rg.METRIC,
    "bootstrap": rg.BOOTSTRAP,
    "control": rg.CONTROL,
    "permutations": rg.PERMUTATIONS,
    "control_alpha": rg.CONTROL_ALPHA,
    "decision": rg.DECISION,
    "falsifier": rg.FALSIFIER,
    "stop_rule": rg.STOP_RULE,
    "outcomes": rg.OUTCOMES,
    "ungated_clauses": rg.UNGATED_CLAUSES,
    "limitations": list(rg.LIMITATIONS),
    "min_genes": rg.MIN_GENES,
    "pseudocount": rg.PSEUDOCOUNT,
    "runtime": rg.RUNTIME,
    "closed_form_tolerance": rg.CLOSED_FORM_TOLERANCE,
    "gate_tolerance": rg.GATE_TOLERANCE,
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
    bad = [k for k, v in BOUND.items() if reg.get(k) != json.loads(json.dumps(v))]
    if bad:
        raise SystemExit(f"refused: the registration and readergate.py disagree on {', '.join(bad)}")
    moved = [p for p, digest in reg["pinned_sha256"].items() if Path(p).is_file() and _sha(Path(p)) != digest]
    group = mf.files_entry("reader v1 DNase narrowPeak files", [p for p in peak_files(CHROMS) if p.exists()])
    if reg["pinned_sha256"].get(group["path"]) != group["sha256"]:
        moved.append(group["path"])
    if moved:
        raise SystemExit(f"refused: pinned inputs changed since the registration: {moved}")
    return reg


def run_all() -> tuple[list[dict], Counter, dict]:
    rates, borrowed, prov, lent = load_rates()
    texts = {c: (Path(COMPILED_DIR) / f"noncoding_{c}.bio").read_text() for c in CHROMS}
    singles, why = rg.single_rule_genes(texts)
    counts = Counter(why)
    genes, programs = [], {}
    for chrom in CHROMS:
        info, rows, c, _module = rg.run_chromosome(chrom, texts[chrom], singles, rates, borrowed, prov, lent)
        programs[chrom] = {k: info[k] for k in ("base_sha256", "gated_sha256", "elements", "rules_gated")}
        counts.update(c)
        genes += rows
        print(chrom, dict(c), flush=True)
    return genes, counts, programs


def load_rna() -> tuple[dict[str, dict[str, float]], dict]:
    """cell -> gene_id -> mean TPM over the chosen experiment's files; per-file counts."""
    src = json.loads(rg.RNA_SOURCES.read_text())
    tpm, counts = {}, {}
    for cid in rg.CELLS:
        tables = []
        for f in src["biosamples"][cid]["files"]:
            t, c = rg.load_quant(Path(f["path"]))
            tables.append(t)
            counts[f["accession"]] = {**c, "cell": cid}
        tpm[cid] = rg.mean_tpm(tables)
    return tpm, counts


def manifest(rna_read: bool) -> dict:
    reg = json.loads(REGISTRATION.read_text())
    inputs = [
        *(mf.input_entry(Path(COMPILED_DIR) / f"noncoding_{c}.bio", partition=None) for c in CHROMS),
        mf.input_entry(RATE_TABLE, partition=None),
        mf.input_entry(rg.GENE_TABLE, partition=None),
        mf.files_entry("reader v1 DNase narrowPeak files", [p for p in peak_files(CHROMS) if p.exists()]),
        mf.input_entry(REGISTRATION, partition=None),
        mf.input_entry(rg.RNA_SOURCES, partition=None),
    ]
    if rna_read:
        inputs += [
            mf.input_entry(p, partition="measured outcome, opened after the harness checks")
            for p in reg["pinned_sha256"]
            if p.startswith(rg.RNA_DIR.as_posix())
        ]
    return {
        "sources": [
            {"accession": "the 24 compiled non-coding programs M1 ran", "version": "pinned by sha256"},
            {"accession": "reader v1 ENCODE DNase-seq narrowPeak files, 13 biosamples", "version": "pinned"},
            {"accession": "ENCODE polyA plus RNA-seq gene quantifications, V29", "version": "pinned"},
            {"accession": "data/results/gene_rates.json (mouse NIH 3T3, e921fae)", "version": "pinned"},
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38; RNA genes matched by gene_id through data/cache/gencode_genes.tsv",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "selection": rg.SELECTION,
            "runtime": rg.RUNTIME,
            "min_genes": rg.MIN_GENES,
            "pseudocount": rg.PSEUDOCOUNT,
            "bootstrap": rg.BOOTSTRAP,
            "permutations": rg.PERMUTATIONS,
            "control_alpha": rg.CONTROL_ALPHA,
        },
        "exclusions": [
            "genes with several compiled rules or one measured rule",
            "genes the bridge does not resolve in every reader cell",
            "genes whose element is open in all thirteen reader cells or in none",
            "genes whose symbol is not one gene_id, or whose gene_id is not on one row of every table",
        ],
        "partitions": "n/a: every qualifying gene is scored",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main() -> int:
    t0 = time.time()
    reg = check_registration()
    genes, counts, programs = run_all()
    bad = [
        g["gene"]
        for g in genes
        if abs(g["m"] - g["closed_form"]) > rg.CLOSED_FORM_TOLERANCE or g["gate_mismatch"]
    ]
    base = {
        "lane": rg.LANE,
        "registration": "data/results/readergate_registration.json",
        "registration_date": reg.get("date"),
        "question": rg.QUESTION,
        "structure": rg.STRUCTURE,
        "run_counts": dict(counts),
        "programs": programs,
        "harness_check": {
            "genes": len(genes),
            "failing": len(bad),
            "max_abs_m_minus_closed_form": max((abs(g["m"] - g["closed_form"]) for g in genes), default=None),
            "gate_mismatch_cells": sum(len(g["gate_mismatch"]) for g in genes),
        },
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
    }
    if bad:
        out = {
            **base,
            "outcome": "HARNESS_DEFECT",
            "reading": rg.reading("HARNESS_DEFECT", bad=len(bad), n=len(genes)),
            "seconds": round(time.time() - t0, 1),
        }
        print("saved", save_result(RESULT, out, manifest=manifest(rna_read=False)))
        print(out["reading"])
        return 1
    tpm, rna_counts = load_rna()
    rows, why = rg.assemble(genes, tpm, rg.symbol_gene_ids())
    scored = rg.score(rows)
    cols = ("gene", "symbol", "gene_id", "chrom", "element", "action", "m", "compiled_cell", "open", "D")
    out = {
        **base,
        "rna_tables": rna_counts,
        "population_exclusions": dict(why),
        **scored,
        "cells": list(rg.CELLS),
        "genes_table": {"columns": list(cols), "rows": [[r[c] for c in cols] for r in rows]},
        "seconds": round(time.time() - t0, 1),
    }
    print("saved", save_result(RESULT, out, manifest=manifest(rna_read=True)))
    print(json.dumps({k: out.get(k) for k in ("outcome", "reading", "genes", "pairs")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
