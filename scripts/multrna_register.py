# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register lane-multrna before it is run (data/results/multrna_registration.json).

    uv run python scripts/multrna_register.py

Every binding word is imported from `genomeos.attribution.multrna` (its BOUND constants), and
scripts/multrna_run.py refuses to run if any registered field differs from that module or any pinned
file's sha256 differs from the file on disk. The census (data/results/multrna_census.json, no RNA read)
supplies the counts the design rests on, and must have reproduced 7d's census with a clean harness. The
RNA quantifications are 7d's (M2b's), pinned by the sha256 their fetch recorded; this script reads no
value of them.

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
from genomeos.attribution import multrna as mr  # noqa: E402
from genomeos.attribution import readergate as rg  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.results import save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from combine2_register import rna_inputs  # noqa: E402
from multrna_census import OWN_CODE  # noqa: E402
from multrna_run import (  # noqa: E402
    C2_CENSUS,
    C2_REGISTRATION,
    C2_RESULT,
    CENSUS,
    CHROMS,
    M2B_RESULT,
    M2E_RESULT,
)
from readergate_census import peak_files  # noqa: E402
from runrules_census import RATE_TABLE  # noqa: E402

RESULT = "multrna_registration"
WILL_WRITE = "data/results/multrna_result.json"

EXPLORATION_PRECEDED_THIS = (
    "stated so that nothing here is read as a blind prediction of figures this lane had seen. Before"
    " writing it the lane read: combine2.py, saturate.py and their census, register and run scripts;"
    " 7d's committed census (data/results/combine2_census.json, no RNA read) and, from 7d's committed"
    " result, its registered reading, its two outcomes, its per-form gene and pair counts, its"
    " head-to-head counts and interval, its population exclusion counts and its run time; 7f's commit"
    " messages and M2e's census sentence; a chr21 dry run of this lane's census (159 multiplicative genes"
    " run, harness clean, not saved); and this lane's census. The registered readings 7d and 7f reached"
    " were in the coordinator's brief. This lane looked at no TPM value of any gene in any biosample and"
    " at no row of 7d's genes tables (the P and D of 7d's populations); this script loads 7d's result"
    " file only to quote its reading and outcomes"
)
DECISIONS = (
    "design: 7d's, exactly. The endpoint (ENCODE polyA plus RNA-seq of the thirteen reader biosamples,"
    " log2(TPM+1)), POPULATION (a)-(c) and (f), the no-rule baseline, Delta MSE, the gene bootstrap, the"
    " shuffled-reader control, MIN_GENES, the decision and the outcome sentences are combine2's,"
    " unchanged; only (d) runnable and (e) informative are the multiplicative form's",
    "form: the multiplicative composition 7f registered (Zhou et al. 2024; saturate.fit,"
    " saturate.fitted_module, runtime/grn.py's opt-in 'multiplicative' of 1f8f12b). It has no parameter"
    " and nothing is fitted on any RNA value. combine2.py, saturate.py and grn.py are imported, not"
    " edited",
    "intervals: 7d judged two forms and Bonferroni-adjusted each to 97.5%; the multiplicative form keeps"
    " 97.5% so that the three readings stand on one interval. Two head-to-heads are judged here where 7d"
    " judged one, so each is Bonferroni-adjusted to 97.5% (7d's single head-to-head was at 95%)",
    "head-to-heads: over the genes in all three populations, the only genes on which the three forms can"
    " be compared on equal terms; the multiplicative form's own comparison against the baseline is on its"
    " larger population, as 7d's were on theirs",
    "comparators: 7d's two forms are recomputed in the run by 7d's committed code and checked against the"
    " committed result (populations, P and D, the registered reading word for word), so that each"
    " head-to-head uses the same D the multiplicative rows use and a defect on either side shows before"
    " the multiplicative comparison is computed. 7d's readings against the baseline are quoted, not"
    " re-decided",
    "no K562 hold runs (7f's stage s7) are run: this test needs only the thirteen gated levels and the"
    " single-removal harness, as 7d's",
)


def census_facts() -> dict[str, Any]:
    c = json.loads(CENSUS.read_text())
    rs = c["reader_selection"]
    m = rs["forms"][mr.FORM]
    refusals = []
    if c.get("seven_d_census_reproduced", {}).get("problems", ["not checked"]):
        refusals.append(f"7d's census not reproduced: {c.get('seven_d_census_reproduced')}")
    if not c["ran_the_runtime"] or m["harness_failing_genes"] or m["run"] != m["informative"]:
        refusals.append("the census did not run every informative gene with a clean harness")
    if refusals:
        raise SystemExit(f"refused: {refusals}")
    return {
        "census": "data/results/multrna_census.json",
        "sentence": c["sentence"],
        "reader_selection": rs,
        "seven_d_census_reproduced": c["seven_d_census_reproduced"],
    }


def inputs() -> list[dict[str, Any]]:
    return [
        *(mf.input_entry(Path(COMPILED_DIR) / f"noncoding_{c}.bio", partition=None) for c in CHROMS),
        mf.input_entry(RATE_TABLE, partition=None),
        mf.input_entry(rg.GENE_TABLE, partition=None),
        mf.files_entry("reader v1 DNase narrowPeak files", [p for p in peak_files(CHROMS) if p.exists()]),
        mf.input_entry(CENSUS, partition=None),
        mf.input_entry(C2_CENSUS, partition=None),
        mf.input_entry(C2_REGISTRATION, partition=None),
        mf.input_entry(
            C2_RESULT,
            partition="the comparator: 7d's committed result; its reading quoted, its genes tables not read",
        ),
        mf.input_entry(M2B_RESULT, partition="pinned by 7d's registration"),
        mf.input_entry(M2E_RESULT, partition="context: M2e's census sentence, copied verbatim"),
        mf.input_entry(rg.RNA_SOURCES, partition=None),
        *rna_inputs(),
    ]


def payload(entries: list[dict[str, Any]]) -> dict[str, Any]:
    c2res = json.loads(C2_RESULT.read_text())
    m2e = json.loads(M2E_RESULT.read_text())
    out: dict[str, Any] = {
        "lane": mr.LANE,
        "date": date.today().isoformat(),
        "roadmap_item": "§5 item 7h, M2f: the multiplicative composition against measured RNA, registered",
        **mr.registered_fields(),
    }
    out.update(
        {
            "census_facts": census_facts(),
            "seven_d_committed_verbatim": {
                "source": "data/results/combine2_result.json",
                "reading": c2res["reading"],
                "outcomes": c2res["outcomes"],
            },
            "m2e_census_sentence_verbatim": m2e["census"]["sentence"],
            "exploration_preceded_this": EXPLORATION_PRECEDED_THIS,
            "decisions": list(DECISIONS),
            "pinned_sha256": {x["path"]: x["sha256"] for x in entries},
            "what_the_run_will_write": WILL_WRITE,
            "the_code_that_runs_it": (
                "genomeos/attribution/multrna.py, scripts/multrna_census.py and scripts/multrna_run.py,"
                " committed before this registration, with the committed combine2.py, saturate.py,"
                " runtime/grn.py, scripts/combine2_census.py and scripts/combine2_run.py they import"
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
            {"accession": "the 24 compiled non-coding programs 7d ran", "version": "pinned by sha256"},
            {"accession": "reader v1 ENCODE DNase-seq narrowPeak files, 13 biosamples", "version": "pinned"},
            {"accession": "ENCODE polyA plus RNA-seq gene quantifications, V29", "version": "pinned"},
            {"accession": "data/results/gene_rates.json (mouse NIH 3T3, e921fae)", "version": "pinned"},
        ],
        "inputs": entries,
        "input_count": len(entries),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "form": mr.FORM,
            "comparator_forms": list(mr.COMPARATOR_FORMS),
            "min_genes": mr.MIN_GENES,
        },
        "exclusions": ["as POPULATION states"],
        "partitions": "n/a: a registration; no outcome read",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }
    print("saved", save_result(RESULT, out, manifest=manifest))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
