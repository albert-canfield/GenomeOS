# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register the rule-number source census before a single source is fetched
(data/results/rule_number_sources_registration.json).

    uv run --frozen python scripts/rule_number_sources_register.py

It reads nothing but its own code and the digest of the one program the census is about. The
question, the population, the three axes, the cascade, the search plan, what will be reported and
the denominator all come from `genomeos.provenance.rule_number_sources`, so the registration and the code
that applies it cannot drift apart.

No class is assigned here, no source is fetched, no program is written to, no number is changed and
no model request is made.
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.provenance import rule_number_sources as rns  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "rule_number_sources_registration"
WILL_WRITE = "data/results/rule_number_sources_census.json"

#: The measurements this census follows from, by the commits that made them.
DESCRIBES = (
    "cd43fce / 00e7028 / 2e1f676 - the predecessor census: data/demo/gastrulation.bio has 42 scoped "
    "slots and ZERO sourced numbers (S_quoted 0, S_derived 0, E_existence_only 30, A_asserted_bare "
    "8, D_language_default 4), and corpus-wide tier 1 is 161 numbers of which 8 are quoted from a "
    "cited source",
    "53d9936 / b1e8839 - three different sourced repressors of SOX17 each gave an answer identical "
    "to the digit, because the module's own `threshold 2.0; hill 3` makes any added inhibitor a "
    "switch and not a brake",
    "8bb9123 / 19423bd - tuning two free numbers to the unsourced expected_* proportions, refused; "
    "this lane inherits that refusal",
    "5cbce26 - the Tbxt-SOX17 sign correction, which stands",
)

#: This lane's own files.
OWN_CODE = (
    "genomeos/provenance/rule_number_sources.py",
    "scripts/rule_number_sources_register.py",
    "scripts/rule_number_sources_census.py",
    "tests/test_rule_number_sources.py",
)


def inputs() -> list[dict[str, Any]]:
    """The one program the census reads, digested before any source is fetched."""
    return [mf.input_entry(Path(rns.PROGRAM), partition="the hand-authored BioLang program under census")]


def main() -> None:
    payload = rns.registration()
    payload["result"] = RESULT
    payload["date"] = date.today().isoformat()
    payload["lane"] = "lane-kinetics"
    payload["describes"] = list(DESCRIBES)
    payload["what_the_run_will_write"] = WILL_WRITE
    entries = inputs()
    payload["frozen_inputs"] = [
        {k: e[k] for k in ("path", "label", "sha256", "bytes") if k in e} for e in entries
    ]
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": rns.PROGRAM,
                "version": "the working tree at the code revision stamped below, digested in inputs",
            },
            {
                "accession": "genomeos/runtime/grn.py, read for what each rule field IS in the runtime",
                "version": "the same code revision; transcribed into "
                "rule_number_sources.QUANTITY_DEFINITIONS",
            },
            {
                "accession": "the published literature, reached by free web fetch only",
                "version": "each finding of the census names the URL it fetched and what it retrieved",
            },
        ],
        "inputs": entries,
        "method": (
            "Pre-registration of a source census. The population is enumerated from the program by "
            "`rule_number_sources.slots`, never asserted. Three axes are fixed here and not "
            "afterwards: the class (what the literature holds for the quantity), the attribution of "
            "the program's own citation, and whether a literature value is commensurable with the "
            "field at all. A parameter obtained by fitting a model is recorded as fitted and is "
            "never counted as measured; Hill coefficients reach that class by construction, which "
            "is declared here before the search."
        ),
        "assembly": "n/a: no genomic sequence or coordinate is read; the population is a "
        "hand-authored program.",
        "coordinates": "n/a: the census addresses declarations by file, line and field name.",
        "parameters": {
            "program": rns.PROGRAM,
            "rule_fields": list(rns.RULE_FIELDS),
            "population_size": payload["population_size"],
            "classes": dict(rns.CLASSES),
            "cascade": list(rns.CASCADE),
            "cited_source_attribution_values": dict(rns.ATTRIBUTION),
            "commensurability_values": dict(rns.COMMENSURABILITY),
            "required_finding_keys": list(rns.REQUIRED_FINDING_KEYS),
            "search_plan": list(rns.SEARCH_PLAN),
            "denominator": payload["denominator"],
            "reports": list(rns.REPORTS),
            "findings_at_registration": len(rns.FINDINGS),
            "model_requests": 0,
            "money_spent": "none; no paid API is called. Web fetches of the cited works are free.",
        },
        "exclusions": list(rns.EXCLUSIONS),
        "partitions": dict(rns.CLASSES),
        "code_cleanliness": mf.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
    }
    path = save_result(RESULT, payload)
    print(f"{RESULT}: {path}")
    print(f"  population: {payload['population_size']} rule numbers in {rns.PROGRAM}")
    print(f"  classes, fixed now: {', '.join(rns.CLASSES)}")
    print(f"  attribution values: {', '.join(rns.ATTRIBUTION)}")
    print(f"  commensurability values: {', '.join(rns.COMMENSURABILITY)}")
    print(f"  findings in the table at registration: {len(rns.FINDINGS)}")
    print(f"  denominator: {payload['denominator']}")
    print(f"  the run will write: {WILL_WRITE}")
    print("  no slot classified, no source fetched, no number changed, 0 model requests")


if __name__ == "__main__":
    main()
