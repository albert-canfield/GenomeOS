#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-rootpick: apply the registered root rule to lane-capsumrna's committed counts, once.

    uv run python scripts/rootpick_count.py

READS NO RNA. Its only input is data/results/capsumrna_census.json, whose own text states its figures
are "a count on the predictor side, before any RNA value is read". The rule was registered and its code
frozen at 55a2cff before this ran, and `rootpick` imports no outcome source (a static guard in
tests/test_rootpick.py asserts that).

THE INPUT PATH IS RELATIVE ON PURPOSE. `genomeos.manifest.input_entry` records `str(path)` verbatim
while its sibling `files_entry` relativises, so an absolute path here would make this result
unrebuildable from any other checkout -- the defect that leaves 11 committed results unrebuildable
(lane-silencerebuild, 2026-10-03). manifest.py is inside the paid study's frozen closure and cannot be
fixed, so the way not to add a twelfth is to hand it a repo-relative path.
"""

from __future__ import annotations

import json
from pathlib import Path

from genomeos import manifest as mf
from genomeos.attribution import rootpick as rp
from genomeos.results import save_result

RESULT = "rootpick_result"
CENSUS_REL = Path("data/results/capsumrna_census.json")
MIN_GENES = 100

#: The paths this lane is answerable for. Everything else uncommitted in this shared checkout is a
#: peer's and `code_cleanliness` reports it apart, which is the point of the split.
OWN_CODE = (
    "genomeos/attribution/rootpick.py",
    "scripts/rootpick_count.py",
    "tests/test_rootpick.py",
)


def manifest() -> dict:
    inputs = [mf.input_entry(CENSUS_REL, partition=None)]
    return {
        "sources": [
            {
                "accession": "data/results/capsumrna_census.json (lane-capsumrna, 2026-10-10)",
                "version": "pinned by sha256",
            }
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "root_rule": rp.ROOT_RULE,
            "ground": rp.GROUND,
            "min_genes_gate": MIN_GENES,
            "verdict_by_status": rp.VERDICT_BY_STATUS,
        },
        "exclusions": [
            "no RNA, no measured level, no fit residual and no outcome: this is a reclassification of"
            " committed predictor-side counts under a rule registered before it ran"
        ],
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
        "partitions": {
            "none": "a single population, lane-capsumrna's reader selection; no partition is taken, and"
            " the absence of one is not a claim of independence"
        },
    }


def main() -> int:
    census = json.loads(CENSUS_REL.read_text())["reader_selection"]
    counts = census["capped_sum_status"]
    out = rp.reclassify(counts)

    would_be_capped = census["two_roots"]["informative_on_the_capped_root"] + counts["capped_root_only"]
    payload = {
        "question": (
            "which root the capped sum's fit takes where two exist, and what that does to the"
            " population -- decided from the bridge's licensing and not from which root fits better"
        ),
        "root_rule": rp.ROOT_RULE,
        "ground": rp.GROUND,
        "falsifier": rp.FALSIFIER,
        "inconclusive": rp.INCONCLUSIVE,
        "prior_exposure": rp.PRIOR_EXPOSURE,
        "does_not_claim": rp.DOES_NOT_CLAIM,
        "forecloses": rp.FORECLOSES,
        "forms": rp.FORMS,
        "as_committed_by_capsumrna": counts,
        "under_this_rule": out,
        "the_rejected_choice": {
            "if_the_capped_root_were_taken": would_be_capped,
            "min_genes_gate": MIN_GENES,
            "would_have_crossed_the_gate": would_be_capped > MIN_GENES,
            "note": (
                "recorded because the rule chosen is the one that LOSES the testable population;"
                " the rejected alternative is the one that would have bought a result"
            ),
        },
        "committed_readings_that_would_move": {
            "capsumrna_census.reader_selection.capped_sum_status.capped_root_only": {
                "committed": counts["capped_root_only"],
                "under_this_rule": 0,
                "why": "its single fit rests on a ceiling the bridge does not identify",
            },
            "capsumrna_census.reader_selection.unique_fit.differs_from_the_mean.genes": {
                "committed": census["unique_fit"]["differs_from_the_mean"]["genes"],
                "under_this_rule": 0,
                "why": (
                    "the uncapped root gives the mean's level on every admitted gene, and the only"
                    " status whose prediction differed is UNKNOWN under this rule"
                ),
            },
            "the_NOT_TESTABLE_verdict_s_basis": {
                "committed": (
                    "NOT_TESTABLE because only 3 unique genes differ from the mean, below MIN_GENES 100"
                ),
                "under_this_rule": (
                    "stronger and different in kind: ZERO genes differ from the mean, so on this"
                    " population the capped sum is not a DISTINCT FORM at all rather than a distinct"
                    " form with too few genes to test"
                ),
            },
        },
        "nothing_was_rewritten": (
            "no committed result's bytes were touched. The three readings above are reported for the"
            " coordinator to place with the supervisor, as instructed, and this lane stopped here"
        ),
    }
    print("saved", save_result(RESULT, payload, manifest=manifest()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
