# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register M2b before it is run (data/results/readergate_registration.json).

    uv run python scripts/readergate_register.py

Every binding word is imported from `genomeos.attribution.readergate`, the module that applies it, and
scripts/readergate_run.py refuses to run if any registered field differs from that module or any
pinned file's sha256 differs from the file on disk. The census (data/results/readergate_census.json, no
RNA read) supplies the counts the design rests on. The RNA quantifications are pinned by the sha256 the
fetch recorded; this script reads no value of them.

No request is sent, no key is read, no outcome is computed.
"""

from __future__ import annotations

import json
import sys
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import compile as cp  # noqa: E402
from genomeos.attribution import readergate as rg  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from readergate_census import CHROMS, peak_files  # noqa: E402
from readergate_fetch import OWN_CODE  # noqa: E402
from runrules_census import RATE_TABLE  # noqa: E402

RESULT = "readergate_registration"
CENSUS = RESULTS_DIR / "readergate_census.json"
WILL_WRITE = "data/results/readergate_result.json"

EXPLORATION_PRECEDED_THIS = (
    "stated so that nothing here is read as a blind prediction of figures this lane had seen. Before"
    " writing it the lane read: the count of rule lines per target gene over the 24 generated programs"
    " by grep (18,973 genes with any rule, 744 with exactly one rule line); a prototype gate, parse and"
    " run of chr21 and chr22 (24 single-rule genes, 17 with an element open in some reader cells and not"
    " others, 9 activating and 8 inhibiting, closed form and gate checks clean); this lane's census"
    " (data/results/readergate_census.json, no RNA read); the ENCODE portal's metadata for every polyA"
    " plus and total RNA-seq experiment and gene-quantification file on the thirteen biosamples (titles,"
    " dates, labs, biosample summaries, replicate counts, sizes, md5); the header line of each fetched"
    " quantification file; and the header and the K562 RNA-seq rows of the saved AlphaGenome track"
    " metadata. No TPM value of any gene in any biosample was read, so no D, no sign agreement and no"
    " outcome was seen before this file was committed"
)
DEPARTURES = (
    "from the coordinator's brief: the reader layer is applied to the 24 generated programs in memory and"
    " is not written into the generated store, because five scripts and a census glob that store or"
    " data/**/*.bio and would count every rule twice; no gated chr21 is committed, because the tracked"
    " corpus is the population of three census tests (rule_evidence_tier, rule_cell_provenance,"
    " compiled_defaults) and a fourth chr21 program would move their figures for nothing the in-memory"
    " build lacks. The sha256 of every gated program is recorded in the census instead",
    "from the engine lane's list of what compile.py emits: under the registered selection every"
    " predicted rule's `when` names the reader cell types instead of its compiled cell (the argmax track,"
    " which data/results/argmaxcell.json measured as the screened cell for 19.0% of K562 elements, below"
    " its registered 25% usability bar). With `when` kept the reader could only hold a rule in its"
    " compiled cell, and a run could make no cross-cell prediction that does not restate its input."
    " Both selections are implemented (compile.SELECTIONS); the as-compiled one is reported as a"
    " descriptive",
    "not_open_in_reader is emitted as `closed`, reader v1's own word for a window no narrowPeak"
    " overlaps (compile.READER_CLOSED_MEANS); not_assessable as `unknown`. The engine holds both alike",
)


def census_facts() -> dict[str, Any]:
    c = json.loads(CENSUS.read_text())
    e = c["emitted"]
    return {
        "census": "data/results/readergate_census.json",
        "programs": e["programs"],
        "elements": e["elements"],
        "rules_gated": e["rules_gated"],
        "element_states_by_cell": e["element_states_by_cell"],
        "population_census": c["population_census"],
    }


def rna_inputs() -> list[dict[str, Any]]:
    src = json.loads(rg.RNA_SOURCES.read_text())
    out = []
    for cid in rg.CELLS:
        b = src["biosamples"][cid]
        for f in b["files"]:
            out.append(
                mf.input_entry(
                    f["path"],
                    partition="the measured outcome, pinned by sha256 only; no value read here",
                    cell=cid,
                    experiment=b["experiment"],
                    accession=f["accession"],
                    fetched_sha256=f["sha256"],
                )
            )
    bad = [x["path"] for x in out if x["sha256"] != x["fetched_sha256"]]
    if bad:
        raise SystemExit(f"refused: the files on disk are not the ones fetched: {bad}")
    return out


def inputs() -> list[dict[str, Any]]:
    return [
        *(mf.input_entry(Path(COMPILED_DIR) / f"noncoding_{c}.bio", partition=None) for c in CHROMS),
        mf.input_entry(RATE_TABLE, partition=None),
        mf.input_entry(rg.GENE_TABLE, partition=None),
        mf.files_entry("reader v1 DNase narrowPeak files", [p for p in peak_files(CHROMS) if p.exists()]),
        mf.input_entry(CENSUS, partition=None),
        mf.input_entry(rg.RNA_SOURCES, partition=None),
        *rna_inputs(),
    ]


def payload(entries: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "lane": rg.LANE,
        "date": date.today().isoformat(),
        "roadmap_item": "§5 item 14, M2b: the reader construct, application half",
        "question": rg.QUESTION,
        "structure": rg.STRUCTURE,
        "census_facts": census_facts(),
        "exploration_preceded_this": EXPLORATION_PRECEDED_THIS,
        "departures": list(DEPARTURES),
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
        "descriptives_not_verdicts": list(rg.DESCRIPTIVES),
        "min_genes": rg.MIN_GENES,
        "pseudocount": rg.PSEUDOCOUNT,
        "runtime": rg.RUNTIME,
        "closed_form_tolerance": rg.CLOSED_FORM_TOLERANCE,
        "gate_tolerance": rg.GATE_TOLERANCE,
        "pinned_sha256": {x["path"]: x["sha256"] for x in entries},
        "what_the_run_will_write": WILL_WRITE,
        "the_code_that_runs_it": "genomeos/attribution/readergate.py, genomeos/attribution/compile.py and"
        " scripts/readergate_run.py, committed before this registration",
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
    }


def main() -> int:
    entries = inputs()
    out = payload(entries)
    manifest = {
        "sources": [
            {"accession": "the 24 compiled non-coding programs M1 ran", "version": "pinned by sha256"},
            {"accession": "reader v1 ENCODE DNase-seq narrowPeak files, 13 biosamples", "version": "pinned"},
            {"accession": "ENCODE polyA plus RNA-seq gene quantifications, V29", "version": "pinned"},
            {"accession": "data/results/gene_rates.json (mouse NIH 3T3, e921fae)", "version": "pinned"},
        ],
        "inputs": entries,
        "input_count": len(entries),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {"min_genes": rg.MIN_GENES, "pseudocount": rg.PSEUDOCOUNT, "selection": rg.SELECTION},
        "exclusions": ["as POPULATION states"],
        "partitions": "n/a: a registration; no outcome read",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }
    print("saved", save_result(RESULT, out, manifest=manifest))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
