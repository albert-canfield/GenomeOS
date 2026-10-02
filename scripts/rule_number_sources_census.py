# SPDX-License-Identifier: AGPL-3.0-or-later
"""The rule-number source census of data/demo/gastrulation.bio
(data/results/rule_number_sources_census.json).

    uv run --frozen python scripts/rule_number_sources_census.py

It enumerates the program's 21 rule numbers from the program itself and joins each one to the
finding recorded for it in `genomeos.lang.rule_number_sources.FINDINGS`, which holds, per slot, the
URL that was fetched, what was retrieved from it, what the source claims, what it does NOT claim,
and the class. Every count printed is a recount of those findings.

It classifies; it changes nothing. No `.bio` file is written to, no number is adjusted toward any
value found, no run of the program is made, and no model request is sent.
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.lang import rule_number_sources as rns  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "rule_number_sources_census"
REGISTERED_BEFORE = "data/results/rule_number_sources_registration.json"

OWN_CODE = (
    "genomeos/lang/rule_number_sources.py",
    "scripts/rule_number_sources_register.py",
    "scripts/rule_number_sources_census.py",
    "tests/test_rule_number_sources.py",
)


def findings() -> list[dict[str, Any]]:
    """One record per rule number, in the program's own order, with nothing filled in by default."""
    out: list[dict[str, Any]] = []
    for s in rns.slots():
        key = (s["rule"], s["field"])
        if key not in rns.FINDINGS:
            raise SystemExit(f"no finding recorded for {key}: the census may not leave a slot unclassified")
        record = {**s, **rns.FINDINGS[key]}
        missing = [k for k in rns.REQUIRED_FINDING_KEYS if not record.get(k)]
        if missing:
            raise SystemExit(f"finding for {key} is missing {missing}")
        out.append(record)
    return out


def main() -> None:
    recs = findings()
    counts = {c: 0 for c in rns.CLASSES}
    for r in recs:
        counts[r["class"]] += 1
    measured_strict = counts["M_measured_quoted"]
    measured_loose = sum(1 for r in recs if r["class"] == "M_measured_quoted" or r.get("measured_elsewhere"))
    payload: dict[str, Any] = {
        "result": RESULT,
        "date": date.today().isoformat(),
        "lane": "lane-kinetics",
        "registered_before_the_search": REGISTERED_BEFORE,
        "program": rns.PROGRAM,
        "question": rns.registration()["question"],
        "population_size": len(recs),
        "findings": recs,
        "counts": counts,
        "measured_strict": measured_strict,
        "measured_loose": measured_loose,
        "sources_fetched": rns.SOURCES_FETCHED,
        "verdict": rns.VERDICT,
        "attribution_findings": rns.ATTRIBUTION_FINDINGS,
        "contradictions": rns.CONTRADICTIONS,
        "reports": list(rns.REPORTS),
        "whichever_way_it_falls": rns.WHICHEVER_WAY_IT_FALLS,
        "exclusions": list(rns.EXCLUSIONS),
    }
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": rns.PROGRAM,
                "version": "the working tree at the code revision stamped below, digested in inputs",
            },
            {
                "accession": "the published literature, reached by free web fetch only",
                "version": "every URL fetched is named in sources_fetched and on the "
                "findings that rest on it",
            },
        ],
        "inputs": [
            mf.input_entry(Path(rns.PROGRAM), partition="the hand-authored BioLang program under census"),
            mf.input_entry(Path(REGISTERED_BEFORE), partition="this census's own pre-registration"),
        ],
        "assembly": "n/a: no genomic sequence or coordinate is read.",
        "coordinates": "n/a: the census addresses declarations by file, line and field name.",
        "method": (
            "Every one of the program's rule numbers is enumerated from the program and joined to a "
            "hand-recorded finding that names the URL fetched, what was retrieved, what the source "
            "claims and what it does not. The classes, the cascade, the two other axes and the "
            "denominator are the registration's, unchanged. A value obtained by fitting a model is "
            "recorded as fitted and is never counted as measured. Every count is a recount of the "
            "findings; the test suite recounts them independently."
        ),
        "parameters": {
            "classes": dict(rns.CLASSES),
            "cascade": list(rns.CASCADE),
            "cited_source_attribution_values": dict(rns.ATTRIBUTION),
            "commensurability_values": dict(rns.COMMENSURABILITY),
            "denominator": rns.registration()["denominator"],
            "model_requests": 0,
            "money_spent": "none; no paid API is called. Every fetch was a free web fetch.",
        },
        "exclusions": list(rns.EXCLUSIONS),
        "partitions": dict(rns.CLASSES),
        "code_cleanliness": mf.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
    }
    path = save_result(RESULT, payload)
    print(f"{RESULT}: {path}")
    print(f"  population: {len(recs)} rule numbers in {rns.PROGRAM}")
    for c in rns.CASCADE:
        print(f"  {c}: {counts[c]}")
    print(f"  measured with a quoted source, strictly: {measured_strict} of {len(recs)}")
    print(f"  measured with a quoted source, any organism: {measured_loose} of {len(recs)}")
    print(f"  sources fetched: {len(rns.SOURCES_FETCHED)}")
    print("  no number changed, no .bio file written to, 0 model requests")


if __name__ == "__main__":
    main()
