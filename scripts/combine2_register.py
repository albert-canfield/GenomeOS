# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register lane-combine2 before it is run (data/results/combine2_registration.json).

    uv run python scripts/combine2_register.py

Every binding word is imported from `genomeos.attribution.combine2`, the module that applies it, and
scripts/combine2_run.py refuses to run if any registered field differs from that module or any pinned
file's sha256 differs from the file on disk. The census (data/results/combine2_census.json, no RNA read)
supplies the counts the design rests on. The RNA quantifications are M2b's, pinned by the sha256 its
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
from genomeos.attribution import combine2 as c2  # noqa: E402
from genomeos.attribution import readergate as rg  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from combine2_census import CHROMS, OWN_CODE  # noqa: E402
from combine2_run import BOUND, M2B_RESULT  # noqa: E402
from readergate_census import peak_files  # noqa: E402
from runrules_census import RATE_TABLE  # noqa: E402

RESULT = "combine2_registration"
CENSUS = RESULTS_DIR / "combine2_census.json"
WILL_WRITE = "data/results/combine2_result.json"

EXPLORATION_PRECEDED_THIS = (
    "stated so that nothing here is read as a blind prediction of figures this lane had seen. Before"
    " writing it the lane read: the rule lines of the 24 generated programs by a scratch regex (18,973"
    " genes with any rule line, 18,229 with two or more; the summed single-removal drop sum(1 - RHO) over"
    " a gene's activating rules, median 1.84; a prototype of both fits, 4,937 and 2,687 of 18,066 genes"
    " runnable under the mean and the competition form; 76,482 compiled gene-cell pairs with two or more"
    " rules in one compiled context, 52,536 and 45,492 runnable); a prototype gate, parse and run of chr21"
    " and chr22 with both forms (harness clean); this lane's census (data/results/combine2_census.json,"
    " no RNA read); the Methods of Bothma et al. 2015 (Europe PMC full text of PMC4532966); M1's and"
    " M2b's registrations and registered readings (M2b: 'worse', the reader-gated single-rule run further"
    " from measured RNA than the no-rule baseline), which the coordinator's brief quoted. This lane read"
    " no TPM value of any gene in any biosample. The RNA files were opened by M2b's run"
    " (lane-readercompile), whose committed result holds values only for single-rule genes, a set"
    " disjoint from this one by construction (one rule line against two or more)"
)
DECISIONS = (
    "endpoint: M2b's ENCODE polyA plus RNA-seq in the thirteen reader biosamples, not M1's GTEx path. A"
    " reader-gated run is the one cross-cell prediction of a compiled run that does not restate its"
    " input (the open/held pattern is ENCODE DNase, not the argmax track); on the GTEx path a gene is"
    " compared in its rule's compiled tissue, the very track whose deletion the fit reproduces; and M2b's"
    " result on this endpoint is the single-rule context the multi-element run is read beside",
    "forms: the runtime's mean, as the runtime runs it under the reader (a held activator counted at"
    " strength 0), and the promoter-competition form from Bothma et al. 2015 (eLife 4:e07956, Methods"
    " eqs. 7-9), its one free number THETA = 0.5 declared (not measured, not fitted to anything). A"
    " super-additive form (each removal a fixed fold whatever else is held) would carry more genes, since"
    " the summed single-removal drops of most multi-element genes exceed the gene's intact level, but"
    " under the engine's hold (a held rule is a strength-0 copy) an activator's held state cannot be a"
    " fold below 1 without a per-rule floor the runtime does not have; it is left open, not chosen by"
    " how many genes it would admit",
    "the competition form is added to runtime/grn.py as an OPT-IN rule (OPT_IN_ACTIVATOR_RULES, selected"
    " per run by NetworkRuntime(activator_combination=...)); the default stays the mean, ACTIVATOR_RULES"
    " and the 2026-09-28 choice test are untouched, and the engine imports nothing from the application."
    " grn.py is in the code closure the AstroREG-2 registration pins, so this edit moves that digest",
    "inhibitor-only genes are not run: their activator drive is 1 under every activator form and no"
    " observation fixes their element-driven share",
    "two forms are judged, so each form's interval is the Bonferroni 97.5% interval; the head-to-head is"
    " one further comparison at 95%",
)


def census_facts() -> dict[str, Any]:
    c = json.loads(CENSUS.read_text())
    rs = c["reader_selection"]
    return {
        "census": "data/results/combine2_census.json",
        "as_compiled_no_reader_totals": c["as_compiled_no_reader"]["totals"],
        "reader_selection": {k: v for k, v in rs.items() if k != "forms"},
        "forms": rs["forms"],
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
        mf.input_entry(M2B_RESULT, partition="context: M2b's registered reading, copied verbatim by the run"),
        *rna_inputs(),
    ]


def payload(entries: list[dict[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {
        "lane": c2.LANE,
        "date": date.today().isoformat(),
        "roadmap_item": "§5 item 7 'Next' and the open part of item 1: a registered combination form",
    }
    for field, const in BOUND.items():
        out[field] = getattr(c2, const)
    out.update(
        {
            "census_facts": census_facts(),
            "exploration_preceded_this": EXPLORATION_PRECEDED_THIS,
            "decisions": list(DECISIONS),
            "pinned_sha256": {x["path"]: x["sha256"] for x in entries},
            "what_the_run_will_write": WILL_WRITE,
            "the_code_that_runs_it": (
                "genomeos/attribution/combine2.py, genomeos/runtime/grn.py, scripts/combine2_census.py and"
                " scripts/combine2_run.py, committed before this registration"
            ),
            "alphagenome_requests": 0,
            "money": "none: no request is sent and no key is read",
        }
    )
    return out


def main() -> int:
    entries = inputs()
    out = payload(entries)
    manifest = {
        "sources": [
            {
                "accession": "the 24 compiled non-coding programs M1 and M2b ran",
                "version": "pinned by sha256",
            },
            {"accession": "reader v1 ENCODE DNase-seq narrowPeak files, 13 biosamples", "version": "pinned"},
            {"accession": "ENCODE polyA plus RNA-seq gene quantifications, V29", "version": "pinned"},
            {"accession": "data/results/gene_rates.json (mouse NIH 3T3, e921fae)", "version": "pinned"},
        ],
        "inputs": entries,
        "input_count": len(entries),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {"forms": list(c2.FORMS), "theta": c2.THETA, "min_genes": c2.MIN_GENES},
        "exclusions": ["as POPULATION states"],
        "partitions": "n/a: a registration; no outcome read",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }
    print("saved", save_result(RESULT, out, manifest=manifest))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
