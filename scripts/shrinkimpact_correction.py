#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The one additive correction this lane owes: closure_cell_inventory's three DepMap gene tallies.

    uv run python scripts/shrinkimpact_correction.py

Of everything this lane recomputed, exactly three committed leaves do not reproduce on the complete
stores, and they are all in `closure_cell_inventory.json`'s `dependency` block. The lane's rule is
that a number that moves gets an ADDITIVE correction result naming the old and the new value, as its
own file, with its own manifest, and that no committed result is ever rewritten. This is that file.
It is named for what it corrects so a reader finds it from the result, and
`data/results/closure_cell_inventory.json` is opened read-only and left exactly as committed.

WHAT IS BEING CORRECTED IS A FIELD'S MEANING, NOT AN ARITHMETIC ERROR. The three figures were true
of this machine's distilled DepMap cache when they were written on the morning of 2026-10-10. They
are not counts over the census' population, and the census presents them as though they were. The
correction states both numbers, which quantity each is, and what did not change.

NOTHING IS RECOMPUTED HERE. Every number in this file is read from a committed result:
`shrinkimpact_closure_cell_inventory.json` (the recount, 04df48d's script, which measured them),
`closure_cell_inventory.json` (the committed census) and `closure_cell2_inventory.json` (the sibling
of the same day whose own `cache_denominator_note` named the defect first). That is deliberate: a
correction that re-derived its own numbers would be a second measurement needing its own check.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.results import save_result  # noqa: E402

NAME = "shrinkimpact_correction_closure_cell_inventory_dependency"
ENTRY = "scripts/shrinkimpact_correction.py"
OWN_CODE = (
    "scripts/shrinkimpact_correction.py",
    "tests/test_shrinkimpact.py",
)
CORRECTS = "data/results/closure_cell_inventory.json"
RECOUNT = "data/results/shrinkimpact_closure_cell_inventory.json"
SIBLING = "data/results/closure_cell2_inventory.json"
TALLIES = ("genes_asked", "genes_with_a_value", "genes_absent_from_the_release")


def _read(rel: str) -> dict[str, Any]:
    return json.loads((ROOT / rel).read_text())


def correction() -> dict[str, Any]:
    """The correction, every figure of it read from a committed result."""
    census = _read(CORRECTS)
    recount = _read(RECOUNT)
    sibling = _read(SIBLING)
    block = recount["dependency_block_is_a_property_of_the_cache_and_not_of_the_census"]
    moved = next(c for c in recount["comparisons"] if c["label"] == "dependency gene tallies")
    unchanged = next(
        c for c in recount["comparisons"] if c["label"].startswith("dependency.genes_asked, against")
    )
    population = next(c for c in recount["comparisons"] if c["label"].startswith("the frozen population"))
    k562 = next(c for c in sibling["candidates"] if c["cell"] == "K562" and c.get("dependency"))
    return {
        "corrects": CORRECTS,
        "corrects_fields": [f"dependency.{t}" for t in TALLIES],
        "the_committed_result_is_not_rewritten": "data/results/closure_cell_inventory.json is "
        "opened read-only by this lane and stands exactly as committed on 2026-10-10. This file is "
        "additive and carries both values of every figure it touches.",
        "old_and_new": [
            {
                "field": f"dependency.{entry['path']}",
                "committed_2026_10_10": entry["committed"],
                "on_the_complete_stores_2026_10_12": entry["recomputed"],
            }
            for entry in moved["moved"]
        ],
        "which_quantity_each_is": {
            "what_the_committed_figures_are": "a reading of this machine's distilled DepMap cache, "
            "data/knowledge/closure_cell/depmap_k562_gene_effect.json.gz, at the moment the census "
            "was written. `closure_cell.k562_gene_effect` returns the cache's ACCUMULATED UNION "
            "whenever the genes a run asks for are already covered -- the whole cached record comes "
            "back, `requested` included -- and the census writes these three fields from "
            "`len(held['requested'])`, `len(held['values'])` and `len(held['missing'])`.",
            "what_the_field_names_claim": "a count over the genes THIS census asked DepMap about.",
            "the_census_own_ask": block["the_census_own_ask_now"],
            "and_it_is_unchanged": unchanged["reading"],
        },
        "why_the_cache_grew": {
            "reading": "`closure_cell2_inventory`, written later the same day, asked for more genes "
            "through the same persisted file, and its own record carries the three numbers the cache "
            "now holds.",
            "the_siblings_own_figures": {
                k: k562["dependency"][k]
                for k in (
                    "genes_in_the_cache_for_this_cell",
                    "genes_with_a_value",
                    "genes_absent_from_the_release",
                )
            },
            "the_sibling_named_the_defect_first": k562["dependency"]["cache_denominator_note"],
            "what_was_never_done": "nothing carried that note back to the census of the same "
            "morning, which is why a faithful rebuild of `closure_cell_inventory` from its manifest "
            "prints the larger figures into its `dependency` block.",
        },
        "what_does_not_change_and_is_measured_not_argued": {
            "the_frozen_population": population["reading"],
            "every_dependency_call": "the frozen population is compared row by row in the recount "
            "with `dependency_attributed` and `dependency_nearest` on every one of its 47,912 rows, "
            "and every leaf matches. No element's dependency label moves.",
            "the_arm_sizes": census["arms"]["sizes"],
            "m3s_two_registered_populations": {
                "attributed_essential_neighbourhood_not": census["arms"]["sizes"][
                    "attributed_essential_neighbourhood_not"
                ],
                "attributed_not_neighbourhood_not": census["arms"]["sizes"][
                    "attributed_not_neighbourhood_not"
                ],
                "promoter_of_a_gene_k562_depends_on": census["arms"]["sizes"][
                    "promoter_of_a_gene_k562_depends_on"
                ],
            },
            "reading": "35 and 5,427 are unchanged on the complete stores, and so is the positive "
            "control's 757. A recomputed count is what the complete stores give; it does not by "
            "itself overturn a registered reading, and no registered test was re-run.",
        },
        "this_is_not_the_input_shrink": "the defect this lane was opened on is a declared input "
        "group that shrinks silently on a partial store. This is a different defect found on the "
        "way: three fields that are not a function of the declared inputs' content alone. It is "
        "reported, not repaired -- `genomeos/attribution/closure_cell.py` is not this lane's to "
        "change, and the repair is a decision for whoever owns M3.",
        "alphagenome_requests": 0,
        "model_requests": 0,
        "money": "none: every figure here is read from a committed result",
    }


def result_manifest(read_now: list[str]) -> dict[str, Any]:
    """The contract's block. Three committed results, declared one by one, repository-relative."""
    inputs = [
        mf.input_entry(p, partition="read by the correction")
        for p in sorted(set(read_now))
        if (ROOT / p).exists()
    ]
    return {
        "sources": [
            {"accession": CORRECTS, "version": "as committed 2026-10-10; the result corrected"},
            {
                "accession": RECOUNT,
                "version": "as committed 2026-10-12; the recount that measured the movement",
            },
            {
                "accession": SIBLING,
                "version": "as committed 2026-10-10; the sibling whose cache_denominator_note "
                "named the defect first",
            },
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "n/a: no coordinate is read or written here. Every figure is copied from a "
        "committed result, and the results it reads are GRCh38 throughout.",
        "coordinates": {
            "base": None,
            "interval": "n/a",
            "note": "no interval is read or written: this file states two values of three counts "
            "and which quantity each is",
        },
        "parameters": {
            "corrects": CORRECTS,
            "fields": [f"dependency.{t}" for t in TALLIES],
            "recomputes_nothing": True,
            "rewrites_nothing": True,
        },
        "exclusions": [
            "every recomputation: this file derives no number of its own, by design, because a "
            "correction that re-measured would need its own check",
            "every repair: closure_cell.py and closure_cell_inventory.py are not this lane's to "
            "change, and no committed result is rewritten",
            "every registered test: re-running one on a changed population needs a new "
            "registration and is not this lane's to decide",
        ],
        "partitions": "n/a: an additive correction, with no outcome measured",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    mf.trace_begin()
    started = time.time()
    record = correction()
    record["seconds"] = round(time.time() - started, 1)
    record["result_manifest"] = result_manifest(mf.traced_reads())
    path = save_result(NAME, record)
    print(f"wrote {path}", flush=True)
    for entry in record["old_and_new"]:
        print(
            f"  {entry['field']}: {entry['committed_2026_10_10']} -> "
            f"{entry['on_the_complete_stores_2026_10_12']}",
            flush=True,
        )
    print("  nothing recomputed, nothing rewritten", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
