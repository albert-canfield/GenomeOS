# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register the significant-increase feasibility lane before any count is taken
(data/results/increase_registration.json).

    uv run --frozen python scripts/increase_register.py

It reads nothing but its own code and the digests of the files the run will read. The eligibility
predicate, the locus convention, both floors, the ladder's steps, the two populations, the reading of
each outcome and the data-or-code attribution of every step all come from
`genomeos.attribution.increases`, so the registration and the code that applies it cannot drift apart.

No pair is opened, no increase is counted, no direction is read, no extractor is changed and no model
request is made. The inputs the run will read are listed by their own paths with their digests, one
entry per path: `manifest.files_entry` is not used, because resolving an input as a group reports it
absent without a file being opened.
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri  # noqa: E402
from genomeos.attribution import increases as inc  # noqa: E402
from genomeos.attribution import measured as ms  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "increase_registration"
WILL_WRITE = "data/results/increase_population.json"
#: The reading this lane follows from, by its own file and path.
DESCRIBES = ("data/results/repress2_population.json", "data/results/repress2_registration.json")

#: This lane's own files. Everything else uncommitted in this shared checkout belongs to another lane.
OWN_CODE = (
    "genomeos/attribution/increases.py",
    "scripts/increase_register.py",
    "scripts/increase_population.py",
    "tests/test_increases.py",
)


def inputs() -> list[dict[str, Any]]:
    """The named files this registration can digest before any run, one entry per path.

    The compiled rules are re-enumerated from the per-chromosome results on disk by the run itself, so
    here the two readings this lane follows from and the two benchmark tables behind every count are
    frozen by their own paths.
    """
    out: list[dict[str, Any]] = []
    for name in DESCRIBES:
        p = Path(name)
        if p.exists():
            out.append(mf.input_entry(p, partition="the reading this lane follows from"))
    for name in ms.CRISPRI_FILES:
        p = crispri.KNOWLEDGE / name
        if p.exists():
            out.append(mf.input_entry(p, partition=ms.CRISPRI_SPLIT_OF[name]))
    return out


def main() -> None:
    payload = inc.registration()
    payload["date"] = date.today().isoformat()
    payload["describes"] = list(DESCRIBES)
    payload["what_the_run_will_write"] = WILL_WRITE
    payload["order"] = [
        "1. the ladder, blind to direction at every position step, with the exhaustive outcome "
        "breakdown reported at each step's own denominator",
        "2. P1 and P2 against both imported floors, separately and never pooled",
        "3. the verdict, in the words registered here before the count, and the data-or-code "
        "attribution of the step a short population fell at",
        "4. stop: on either outcome this lane reports and tests nothing. A clearing count is handed on "
        "as a proposal and a separate lane, not spent here",
    ]
    entries = inputs()
    payload["frozen_inputs"] = [{k: e[k] for k in ("path", "sha256", "bytes") if k in e} for e in entries]
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "ENCODE CRISPRi enhancer-gene benchmark (EngreitzLab/CRISPR_comparison, "
                "Gschwind et al.), as genomeos.attribution.measured already caches it",
                "version": "the two benchmark tables on this machine, digested in inputs",
            },
            {
                "accession": "the compiled non-coding programs of the 24 chromosomes, as "
                "genomeos.attribution.compile emits them",
                "version": "the attributed element results on this machine; the run digests what it reads",
            },
        ],
        "inputs": entries,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "eligibility_predicate": inc.ELIGIBILITY_CALL,
            "independent_locus_rule": inc.LOCUS_RULE,
            "independent_locus_span": inc.LOCUS_SPAN,
            "locus_floor": inc.LOCUS_FLOOR,
            "link_floor": inc.POSITIVE_FLOOR,
            "reciprocal_overlap": ms.RECIPROCAL_OVERLAP,
            "reach": ms.REACH,
            "ladder_steps": list(inc.LADDER_STEPS),
            "populations": list(inc.POPULATIONS),
            "aggregate_function": "count",
            "fill_value": 0,
        },
        "exclusions": [
            inc.P1_CONTEST_RULE,
            "no pair is excluded from a ladder denominator: the benchmark's own invalid rows "
            "(ValidConnection FALSE) are counted and reported in their own row, as "
            "measured.parse_crispri already returns them, and never read as measured negatives",
            "no rule is excluded from the rule enumeration: every rule the compiler emits is "
            "enumerated, and a rule outside a reported population is counted in its own row",
        ],
        "partitions": {
            "ladder": "cached CRISPRi pairs, one count per step, blind to direction at every step",
            "P1": inc.P1_CALL,
            "P2": inc.P2_CALL,
            "training_and_heldout": inc.P1_SPLIT_RULE,
        },
        "code_cleanliness": mf.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
    }
    path = save_result(RESULT, payload)
    print(f"{RESULT}: {path}")
    print(f"  floors: {inc.LOCUS_FLOOR} independent loci, {inc.POSITIVE_FLOOR} links, both imported")
    print(f"  ladder: {' -> '.join(inc.LADDER_STEPS)}")
    print(f"  populations: {', '.join(inc.POPULATIONS)}")
    print(f"  inputs digested: {len(entries)}")
    print(f"  the run will write: {WILL_WRITE}")
    print("  nothing counted, no direction read, no extractor changed, 0 model requests")


if __name__ == "__main__":
    main()
