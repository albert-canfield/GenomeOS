# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register the link-level direction test before any direction is read
(data/results/direction_link_registration.json).

    uv run --frozen python scripts/direction_link_register.py

It reads nothing but its own code and the digests of the files the run will read. The eligibility
predicate, the sign rule, the arms and their contest rule, the locus convention, both floors, the
interval, the sign-shuffled control and the two readings all come from
`genomeos.attribution.direction_link`, so the registration and the code that applies it cannot drift
apart.

No predicted value is read, no agreement is computed, no extractor is changed and no model request is
made. The inputs the run will read are listed by their own paths with their digests, one entry per
path, except the per-element response cache, whose per-chromosome archives are declared as one group
by `manifest.files_entry` so that every member is named.
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri  # noqa: E402
from genomeos.attribution import direction_link as dl  # noqa: E402
from genomeos.attribution import measured as ms  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "direction_link_registration"
WILL_WRITE = "data/results/direction_link.json"
#: The readings this lane follows from, by their own files.
DESCRIBES = (
    "data/results/increase_population.json",
    "data/results/increase_registration.json",
    "data/results/repress2_population.json",
)
CACHE_DIR = Path("data/knowledge/alphagenome/elements")

#: This lane's own files. Everything else uncommitted in this shared checkout belongs to another lane.
OWN_CODE = (
    "genomeos/attribution/direction_link.py",
    "scripts/direction_link_register.py",
    "scripts/direction_link.py",
    "tests/test_direction_link.py",
)


def inputs() -> list[dict[str, Any]]:
    """The named files this registration can digest before any run."""
    out: list[dict[str, Any]] = []
    for name in DESCRIBES:
        p = Path(name)
        if p.exists():
            out.append(mf.input_entry(p, partition="the reading this lane follows from"))
    for name in ms.CRISPRI_FILES:
        p = crispri.KNOWLEDGE / name
        if p.exists():
            out.append(mf.input_entry(p, partition=ms.CRISPRI_SPLIT_OF[name]))
    archives = sorted(q.as_posix() for q in CACHE_DIR.glob("chr*.json.gz"))
    if archives:
        out.append(
            mf.files_entry(
                "per_element_response_cache",
                archives,
                partition="the finished deletion sweep's own per-element cache, read for a sign only",
            )
        )
    return out


def main() -> None:
    payload = dl.registration()
    payload["result"] = RESULT
    payload["date"] = date.today().isoformat()
    payload["lane"] = "lane-direction"
    payload["describes"] = list(DESCRIBES)
    payload["what_the_run_will_write"] = WILL_WRITE
    entries = inputs()
    payload["frozen_inputs"] = [
        {k: e[k] for k in ("path", "label", "sha256", "bytes", "files") if k in e} for e in entries
    ]
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "ENCODE CRISPRi enhancer-gene benchmark (EngreitzLab/CRISPR_comparison, "
                "Gschwind et al.), as genomeos.attribution.measured already caches it",
                "version": "the two benchmark tables on this machine, digested in inputs",
            },
            {
                "accession": "the finished genome-wide AlphaGenome deletion sweep's per-element "
                "response cache, written with threshold=0.0",
                "version": "the per-chromosome archives on this machine, digested in inputs as a group",
            },
            {
                "accession": "the attributed elements of the 24 chromosomes, as "
                "genomeos.attribution.targets.attributed enumerates them",
                "version": "the per-chromosome results on this machine; the run digests what it reads",
            },
        ],
        "inputs": entries,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "eligibility_predicate": dl.ELIGIBILITY,
            "measured_sign": dl.MEASURED_SIGN_CALL,
            "predicted_direction": dl.PREDICTED_CALL,
            "sign_rule_zero_and_absent": dl.SIGN_CALL,
            "primary_statistic": dl.TEST_STATISTIC,
            "chance": dl.CHANCE,
            "difference_between_the_arms": dl.PRIMARY_STATISTIC,
            "why_the_difference_is_not_the_reading": dl.COLLAPSE,
            "baseline": dl.BASELINE,
            "independent_locus_rule": dl.LOCUS_RULE,
            "independent_locus_span": dl.LOCUS_SPAN,
            "locus_floor": dl.LOCUS_FLOOR,
            "link_floor": dl.POSITIVE_FLOOR,
            "reciprocal_overlap": ms.RECIPROCAL_OVERLAP,
            "reach": ms.REACH,
            "bootstrap_unit": "independent locus",
            "bootstrap_draws": dl.DRAWS,
            "bootstrap_seed": dl.SEED,
            "shuffles": dl.SHUFFLES,
            "shuffle_seed": dl.SHUFFLE_SEED,
            "alpha_one_sided": dl.CONTROL_ALPHA,
            "cells_read": [dl.PRIMARY_CELL],
            "cached_cell_tracks": list(dl.CACHED_CELLS),
            "cells_not_cached": list(dl.CELLS_NOT_CACHED),
            "aggregate_function": "mean of a 0/1 agreement per answered link, per arm",
            "fill_value": "none: an absent or exactly-zero predicted value is excluded and counted",
        },
        "exclusions": [
            dl.CONTEST_RULE,
            dl.SIGN_CALL,
            dl.CELL_AVAILABILITY,
            "no pooled raw agreement rate over the two arms is reported as a result: on a lopsided mix "
            "it is the majority arm's own rate wearing a general noun",
            "no pair is excluded for its magnitude: no effect-size bar is applied anywhere, and the "
            "imported floors are on counts",
        ],
        "partitions": {
            "decreases": "links whose overlapping pairs hold a significant decrease and no significant "
            "increase; measured sign -1",
            "increases": "links whose overlapping pairs hold a significant increase and no significant "
            "decrease; measured sign +1",
            "contested": dl.CONTEST_RULE,
            "training_and_heldout": dl.SPLIT_RULE,
        },
        "code_cleanliness": mf.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
    }
    path = save_result(RESULT, payload)
    print(f"{RESULT}: {path}")
    print(f"  floors: {dl.LOCUS_FLOOR} independent loci, {dl.POSITIVE_FLOOR} links, both imported")
    print(f"  arms: {', '.join(dl.ARMS)}, reported apart and never pooled")
    print(f"  statistic the reading is taken on: {dl.TEST_STATISTIC}")
    print(f"  cached cell tracks: {', '.join(dl.CACHED_CELLS)}; NOT cached: {', '.join(dl.CELLS_NOT_CACHED)}")
    print(f"  interval: {dl.DRAWS} draws over independent loci at seed {dl.SEED}")
    print(f"  control: {dl.SHUFFLES} sign shuffles at seed {dl.SHUFFLE_SEED}, one-sided {dl.CONTROL_ALPHA}")
    print(f"  inputs digested: {len(entries)}")
    print(f"  the run will write: {WILL_WRITE}")
    print("  no direction read, no agreement computed, no extractor changed, 0 model requests")


if __name__ == "__main__":
    main()
