# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register M1 before it is run (data/results/runrules_registration.json).

    uv run python scripts/runrules_register.py

Every binding word is imported from `genomeos.attribution.runrules`, the module that applies it, and
`scripts/runrules_run.py` refuses to run if any registered field differs from that module. The census
(data/results/runrules_census.json, no outcome read) is read for the structural facts the design rests
on. The GTEx table and the compiled programs are pinned by sha256 here, before the run; no value of
the GTEx table is read by this script (hashing reads its bytes, not its numbers).

No request is sent, no key is read, no outcome is computed. 0 model requests, no network, no money.
"""

from __future__ import annotations

import json
import sys
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import runrules as rr  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from runrules_census import CHROMS, OWN_CODE, RATE_TABLE  # noqa: E402

RESULT = "runrules_registration"
CENSUS = RESULTS_DIR / "runrules_census.json"
WILL_WRITE = "data/results/runrules_result.json"

EXPLORATION_PRECEDED_THIS = (
    "stated so that nothing here is read as a blind prediction of figures this lane had seen. Before"
    " writing it the lane read: the rule count per compiled cell context over all 24 programs (318"
    " contexts, by grep of the program text; Whole_Blood 11,325 rules, all 54 GTEx-named contexts"
    " present); data/results/bridge_audit.json's pooled counts (237,613 gene-cell pairs, 76,469 with"
    " several mechanisms, 22,576 resolved on a measured rate, 160,392 on the borrowed median); a"
    " prototype on chr21 only (33 of the 54 panel contexts present, 147 panel pairs resolved on the"
    " borrowed tier, 30 Whole_Blood pairs run with the runtime, largest |m - closed form| 9.4e-16);"
    " the first three lines of the GTEx table (its version line, its dimensions 56,200 x 54, and the"
    " tissue header) and NO value of it; one row of the saved AlphaGenome track metadata (the"
    " Whole_Blood track, data_source 'gtex'); and data/results/argmaxcell.json's verdict headline."
    " No measured expression value of any gene in any tissue was read, so no outcome, no D and no"
    " sign agreement was seen before this file was committed"
)

DEPARTURES = (
    "from ROADMAP §5 item 14 M1 ('the on/off state per cell type of genes the compiled rules drive,"
    " against measured RNA in the reader's cell types'), each with its reason. (1) Not on/off: the"
    " census shows a run predicts 'on' for every resolved gene in every cell (STRUCTURE['on_off']), so"
    " an on/off test cannot separate the run from any baseline; on/off is kept as a descriptive."
    " (2) Not the reader's cells: the measured RNA this machine holds for K562, HepG2, GM12878, IMR-90,"
    " A549 and MCF-7 is six JSON files of bigWig URLs, no values; the one cached per-tissue table is"
    " GTEx v8, and the 54 compiled contexts named for GTEx tissues map to it by name with no judgement."
    " (3) The strength-dependent activator form R3 left open is not attempted here: every pair in this"
    " population has one mechanism, so no combination rule enters it, and that form needs its own"
    " registration",
    "from the coordinator's suggested shape: the baseline is the gene's own level in the tissues where"
    " it has no compiled rule (its 'expression in other cells'), not its basal rate, because the basal"
    " rate is the run's own number and the rate cancels from the prediction (READING)",
)


def census_facts() -> dict[str, Any]:
    """The census numbers the design rests on, read from the committed census result."""
    c = json.loads(CENSUS.read_text())
    keys = (
        "contexts", "rules", "gene_cell_pairs", "single_mechanism_pairs", "several_mechanism_pairs",
        "resolved_measured", "resolved_borrowed", "resolved_borrowed_activates", "resolved_borrowed_inhibits",
        "resolved_borrowed_basal_not_positive",
    )  # fmt: skip
    return {
        "census": "data/results/runrules_census.json",
        "totals": c["totals"],
        "all_contexts": {k: c["all_contexts"].get(k) for k in keys},
        "panel_contexts": {k: c["panel_contexts"].get(k) for k in keys},
        "panel_unresolved_borrowed": c["panel_contexts"].get("unresolved_borrowed"),
        "panel_labels_absent": c["panel_labels_absent"],
    }


def inputs() -> list[dict[str, Any]]:
    out = [mf.input_entry(Path(COMPILED_DIR) / f"noncoding_{c}.bio", partition=None) for c in CHROMS]
    out.append(mf.input_entry(RATE_TABLE, partition=None))
    out.append(mf.input_entry(CENSUS, partition=None))
    out.append(
        mf.input_entry(rr.GTEX, partition="the measured outcome, pinned by sha256 only; no value read here")
    )
    return out


def payload() -> dict[str, Any]:
    entries = inputs()
    return {
        "lane": rr.LANE,
        "date": date.today().isoformat(),
        "roadmap_item": "§5 item 14, M1: compiled rules that run",
        "question": rr.QUESTION,
        "structure": rr.STRUCTURE,
        "census_facts": census_facts(),
        "exploration_preceded_this": EXPLORATION_PRECEDED_THIS,
        "departures": list(DEPARTURES),
        "reading": rr.READING,
        "reading_not_taken": rr.READING_NOT_TAKEN,
        "population": rr.POPULATION,
        "measured": rr.MEASURED,
        "baseline": rr.BASELINE,
        "metric": rr.METRIC,
        "bootstrap": rr.BOOTSTRAP,
        "control": rr.CONTROL,
        "permutations": rr.PERMUTATIONS,
        "control_alpha": rr.CONTROL_ALPHA,
        "decision": rr.DECISION,
        "falsifier": rr.FALSIFIER,
        "stop_rule": rr.STOP_RULE,
        "outcomes": rr.OUTCOMES,
        "limitations": list(rr.LIMITATIONS),
        "descriptives_not_verdicts": list(rr.DESCRIPTIVES),
        "min_reference_tissues": rr.MIN_REFERENCE_TISSUES,
        "min_population": rr.MIN_POPULATION,
        "pseudocount": rr.PSEUDOCOUNT,
        "runtime": rr.RUNTIME,
        "closed_form_tolerance": rr.CLOSED_FORM_TOLERANCE,
        "gtex_tissues": list(rr.GTEX_V8_TISSUES),
        "panel": rr.PANEL,
        "what_the_run_will_write": WILL_WRITE,
        "the_code_that_runs_it": "genomeos/attribution/runrules.py and scripts/runrules_run.py, committed"
        " before this registration",
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
        "network": "none: every input is already on this machine",
        "result_manifest": {
            "sources": [
                {"accession": "the 24 compiled non-coding programs", "version": "pinned by sha256 in inputs"},
                {"accession": "data/results/gene_rates.json (mouse NIH 3T3, e921fae)", "version": "pinned"},
                {"accession": "data/results/runrules_census.json (this lane)", "version": "pinned"},
                {
                    "accession": "GTEx v8 gene median TPM, GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9",
                    "version": "the cached file, pinned by sha256; only its header line was read",
                },
            ],
            "inputs": entries,
            "input_count": len(entries),
            "assembly": "GRCh38 for the compiled programs; GTEx rows are matched by symbol only",
            "coordinates": "n/a: no interval is read, written or compared",
            "parameters": {
                "min_reference_tissues": rr.MIN_REFERENCE_TISSUES,
                "min_population": rr.MIN_POPULATION,
                "pseudocount": rr.PSEUDOCOUNT,
                "bootstrap": rr.BOOTSTRAP,
                "permutations": rr.PERMUTATIONS,
                "control_alpha": rr.CONTROL_ALPHA,
                "runtime": rr.RUNTIME,
                "closed_form_tolerance": rr.CLOSED_FORM_TOLERANCE,
            },
            "exclusions": [
                "every compiled context not named for a GTEx tissue: no hand-made mapping is made",
                "pairs with several mechanisms in their tissue, pooled over programs",
                "genes whose symbol is not on exactly one GTEx row",
                f"genes with fewer than {rr.MIN_REFERENCE_TISSUES} panel tissues free of any compiled rule",
            ],
            "partitions": "n/a: no training partition; every compiled panel pair that qualifies is scored",
            "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
        },
    }


def main() -> int:
    path = save_result(RESULT, payload())
    print(f"registered: {path}")
    print(f"the run will write: {WILL_WRITE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
