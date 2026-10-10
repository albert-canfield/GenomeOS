# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register lane-maxform before it is run (data/results/maxform_registration.json).

    uv run python scripts/maxform_register.py

Every binding word is imported from `genomeos.attribution.maxform` (its BOUND constants), and
scripts/maxform_run.py refuses to run if any registered field differs from that module or any pinned
file's sha256 differs from the file on disk. The census (data/results/maxform_census.json, no RNA read)
supplies the counts the design rests on: `max` reached NOT_TESTABLE there at 0 genes whose prediction
differs from the runtime mean's -- structurally, because it has no fit on a gene with two or more
activating elements -- so the arm registered here is `or`, which differs from the mean on enough genes.
The RNA quantifications are 7d's (M2b's), pinned by the sha256 their fetch recorded; this script reads no
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
from genomeos.attribution import maxform as xf  # noqa: E402
from genomeos.attribution import readergate as rg  # noqa: E402
from genomeos.evidence import COMPILED_DIR  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from combine2_register import rna_inputs  # noqa: E402
from maxform_census import CHROMS, OWN_CODE  # noqa: E402
from readergate_census import peak_files  # noqa: E402
from runrules_census import RATE_TABLE  # noqa: E402

RESULT = "maxform_registration"
WILL_WRITE = "data/results/maxform_result.json"
CENSUS = RESULTS_DIR / "maxform_census.json"
C2_CENSUS = RESULTS_DIR / "combine2_census.json"
C2_REGISTRATION = RESULTS_DIR / "combine2_registration.json"
C2_RESULT = RESULTS_DIR / "combine2_result.json"
MR_RESULT = RESULTS_DIR / "multrna_result.json"
MC_REGISTRATION = RESULTS_DIR / "magcal_registration.json"
MC_RESULT = RESULTS_DIR / "magcal_result.json"

DECISIONS = (
    "the arm: `or`, not `max`. The census (data/results/maxform_census.json, no RNA read) found `max`'s"
    " prediction differs from the runtime mean's on 0 of the 17,983 multi-element genes, because it fits"
    " only the 283 with ONE activating element -- where it IS the mean element for element -- and has no"
    " fit at all on the other 17,432 (maxform.PREMISE). That is NOT_TESTABLE under 7d's MIN_GENES and it"
    " is structural, not a small sample. `or` differs on 1,927 fitted genes, 1,843 of them informative"
    " with a symbol on exactly one gene_id, so it is the arm this registration judges",
    "design: 7d's population and endpoint with lane-magcal's M2h analysis. POPULATION (a)-(c) and (f),"
    " the endpoint (ENCODE polyA plus RNA-seq of the thirteen reader biosamples, log2(TPM+1)), the no-rule"
    " baseline, the centring, MIN_GENES and the harness are combine2's, unchanged; (d) runnable and (e)"
    " informative are `or`'s. The scale fit, the held-out split, the metric, the within-gene cell"
    " permutation control, the four-way decision table and the outcome sentences are magcal's, imported"
    " and not reimplemented; no file of lane-magcal's is edited",
    "the interval: 98.75%, percentiles 0.625 and 99.375, Bonferroni over a family of FOUR -- 7d's runtime"
    " mean, 7d's promoter competition, 7h's multiplicative composition and this arm. It is the same"
    " question, the same RNA and the same split as M2h's three arms, so a fourth arm changes the"
    " correction. Ruled by the supervisor through the coordinator on 2026-10-10 before any value was read."
    " M2h's three readings keep their registered 98.33% intervals and are NOT recomputed",
    "magcal's bootstrap_delta and bootstrap_scale read magcal's own FORM_INTERVAL, so this lane calls"
    " magcal's `_sums` and `_gene_terms` and applies its own percentile pair; tests/test_maxform.py checks"
    " that, run at magcal.FORM_INTERVAL, it reproduces magcal's two functions exactly",
    "head-to-heads against the runtime mean and the multiplicative form are DESCRIPTIVE and decide"
    " nothing: this lane's registered endpoint is the calibrated held-out one, and the comparators are"
    " 7d's and 7h's committed rows, which M2h already calibrated. They are reported because the lane was"
    " asked for them, labelled as deciding nothing",
    "the 1-Mb locus-cluster bootstrap is POST HOC, labelled as such, deciding nothing and never replacing"
    " the registered gene-level interval. It is computed because genes sharing elements are not"
    " independent, which the supervisor's audit of M2h raised",
    "no root is picked anywhere: `or`'s fit is unique in closed form (maxform.FIT"
    " ['the_or_root_is_unique']), so lane-rootpick's registered root convention is cited and applies"
    " vacuously; `max` has no root at all on a gene with two or more activating elements",
    "no AlphaGenome request, no deletion request, no new biosample, nothing in area J, no node-model work,"
    " no BioForge calibration (the stop list of 2026-09-27, reinstated)",
)


def census_facts() -> dict[str, Any]:
    """The census's counts, and the refusals that stop a registration that should not be written."""
    c = json.loads(CENSUS.read_text())
    rs = c["reader_selection"]
    own = rs["this_lanes_forms"][xf.FORM]
    refusals = []
    if c.get("committed_censuses_reproduced", {}).get("problems", ["not checked"]):
        refusals.append(f"7d's and 7h's censuses not reproduced: {c.get('committed_censuses_reproduced')}")
    if not c["ran_the_runtime"] or own["harness_failing_genes"] or own["run"] != own["informative"]:
        refusals.append("the census did not run every informative gene with a clean harness")
    if c["decision"].get(xf.FORM) != xf.DECISIONS[1]:
        refusals.append(f"the census decision for {xf.FORM} is {c['decision'].get(xf.FORM)}")
    checks = c["identity_checks"]
    if checks["max_fits_two_activators"] or checks["max_differs_from_mean_at_one"]:
        refusals.append(f"the identity checks of PREMISE did not return 0: {checks}")
    if refusals:
        raise SystemExit(f"refused: {refusals}")
    return {
        "census": "data/results/maxform_census.json",
        "sentence": c["sentence"],
        "decision": c["decision"],
        "identity_checks": checks,
        "reader_selection": rs,
        "committed_censuses_reproduced": c["committed_censuses_reproduced"],
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
            C2_RESULT, partition="the comparators: 7d's committed genes tables, read in the run only"
        ),
        mf.input_entry(MR_RESULT, partition="the comparator: 7h's committed genes table"),
        mf.input_entry(
            MC_REGISTRATION, partition="M2h's registration, whose split and analysis this lane uses"
        ),
        mf.input_entry(MC_RESULT, partition="M2h's committed reading, quoted and not re-decided"),
        mf.input_entry(rg.RNA_SOURCES, partition=None),
        *rna_inputs(),
    ]


def payload(entries: list[dict[str, Any]]) -> dict[str, Any]:
    mcres = json.loads(MC_RESULT.read_text())
    out: dict[str, Any] = {
        "lane": xf.LANE,
        "date": date.today().isoformat(),
        "roadmap_item": xf.ROADMAP_ITEM,
        **xf.registered_fields(),
    }
    out.update(
        {
            "census_facts": census_facts(),
            "m2h_committed_verbatim": {
                "source": "data/results/magcal_result.json",
                "reading": mcres.get("reading"),
                "per_form": {
                    f: {k: v.get(k) for k in ("outcome", "reading")}
                    for f, v in (mcres.get("per_form") or {}).items()
                },
            },
            "decisions": list(DECISIONS),
            "pinned_sha256": {x["path"]: x["sha256"] for x in entries},
            "what_the_run_will_write": WILL_WRITE,
            "the_code_that_runs_it": (
                "genomeos/attribution/maxform.py, scripts/maxform_census.py, scripts/maxform_register.py"
                " and scripts/maxform_run.py, committed before this registration, with the committed"
                " combine2.py, magcal.py, multrna.py, saturate.py, runtime/grn.py and the census and run"
                " scripts of 7d and 7h they import"
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
            "form": xf.FORM,
            "forms_censused": list(xf.FORMS),
            "form_interval": list(xf.FORM_INTERVAL),
            "min_genes": xf.MIN_GENES,
            "bootstrap": xf.BOOTSTRAP,
            "permutations": xf.PERMUTATIONS,
        },
        "exclusions": ["as POPULATION and POPULATION_NOTE state"],
        "partitions": "n/a: a registration; no outcome read",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }
    print("saved", save_result(RESULT, out, manifest=manifest))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
