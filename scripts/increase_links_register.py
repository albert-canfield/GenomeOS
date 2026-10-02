# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register the versioned measured-layer extractor before `measured.rule_links` changes
(data/results/increase_links_registration.json).

    uv run --frozen python scripts/increase_links_register.py

It reads nothing but its own code and the digests of the files the change will be proved against. The
eligibility predicate, the version rule and its default, the conflict rule, the R2 wording, both floors
and the by-construction statement all come from `genomeos.attribution.increase_links`, so the
registration and the code that applies it cannot drift apart.

No extractor is changed by this script, no link is emitted, no direction is read and no model request is
made. The inputs are listed by their own paths with their digests, one entry per path.
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri  # noqa: E402
from genomeos.attribution import increase_links as il  # noqa: E402
from genomeos.attribution import measured as ms  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "increase_links_registration"
WILL_CHANGE = "genomeos/attribution/measured.py"

#: This lane's own files. Everything else uncommitted in this shared checkout belongs to another lane.
OWN_CODE = (
    "genomeos/attribution/increase_links.py",
    "scripts/increase_links_register.py",
    "tests/test_increase_links.py",
)


def inputs() -> list[dict[str, Any]]:
    """The named files this registration digests before the extractor changes, one entry per path."""
    out: list[dict[str, Any]] = []
    for name in il.INHERITED_FILES:
        p = Path(name)
        if p.exists():
            out.append(mf.input_entry(p, partition="the reading this lane acts on"))
    for name in ms.CRISPRI_FILES:
        p = crispri.KNOWLEDGE / name
        if p.exists():
            out.append(mf.input_entry(p, partition=ms.CRISPRI_SPLIT_OF[name]))
    return out


def main() -> None:
    payload = il.registration()
    payload["date"] = date.today().isoformat()
    payload["what_the_change_will_touch"] = WILL_CHANGE
    entries = inputs()
    payload["frozen_inputs"] = [{k: e[k] for k in ("path", "sha256", "bytes") if k in e} for e in entries]
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "ENCODE CRISPRi enhancer-gene benchmark "
                "(EngreitzLab/CRISPR_comparison, Gschwind et al.), as "
                "genomeos.attribution.measured already caches it",
                "version": "the two benchmark tables on this machine, digested in inputs",
            },
        ],
        "inputs": entries,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "eligibility_predicate": il.ELIGIBILITY_CALL,
            "extractor_versions": list(il.EXTRACTORS),
            "extractor_default": il.DEFAULT_EXTRACTOR,
            "action_of_an_increase_derived_link": il.ACTION,
            "evidence_status_of_an_increase_derived_link": il.EVIDENCE_STATUS,
            "outcome_of_an_increase_derived_link": il.OUTCOME_TEXT,
            "molecular_role_of_an_increase_derived_link": il.MOLECULAR_ROLE_TEXT,
            "conflict_rule": il.CONFLICT_RULE,
            "link_floor": il.POSITIVE_FLOOR,
            "locus_floor": il.LOCUS_FLOOR,
            "reciprocal_overlap": ms.RECIPROCAL_OVERLAP,
            "reach": ms.REACH,
            "aggregate_function": "count",
            "fill_value": 0,
        },
        "exclusions": [
            il.CONFLICT_RULE,
            "no pair is excluded by this registration beyond what the committed extractor already "
            "excludes: the benchmark's own invalid rows (ValidConnection FALSE) stay out of every "
            "count, as measured.parse_crispri already returns them separately, and are never read as "
            "measured negatives",
        ],
        "partitions": {
            "v1": il.EXTRACTOR_RULE[il.EXTRACTOR_V1],
            "v2": il.EXTRACTOR_RULE[il.EXTRACTOR_V2],
            "conflicts": il.CONFLICT_RULE,
        },
        "code_cleanliness": mf.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
    }
    path = save_result(RESULT, payload)
    print(f"{RESULT}: {path}")
    print(f"  versions: {', '.join(il.EXTRACTORS)}; default {il.DEFAULT_EXTRACTOR} (byte-identical)")
    print(f"  an increase-derived link: action {il.ACTION}, {il.EVIDENCE_STATUS}, {il.OUTCOME_TEXT!r}")
    forbidden = ", ".join(il.FORBIDDEN_OF_AN_INCREASE_LINK)
    print(f"  molecular role: {il.MOLECULAR_ROLE_TEXT}; forbidden: {forbidden}")
    print(f"  floors recorded and not applied: {il.POSITIVE_FLOOR} links, {il.LOCUS_FLOOR} loci")
    print(f"  inputs digested: {len(entries)}")
    print(f"  the change will touch: {WILL_CHANGE}")
    print("  no extractor changed yet, no link emitted, no direction read, 0 model requests")


if __name__ == "__main__":
    main()
