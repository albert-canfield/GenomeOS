# SPDX-License-Identifier: AGPL-3.0-or-later
"""Write the registration of the rule-evidence-tier PROPOSAL
(data/results/rule_evidence_tier_proposal_registration.json).

    uv run --frozen python scripts/rule_evidence_tier_register.py

It writes the proposal's design AS REGISTERED and nothing else: the axis and its five values, where
the mark lives at each level, what the runtime does, the default, the storage decision with the
option it rejects, and six falsifiable predictions - all before any of it is implemented and before
any corpus count exists. The adoption cost, the prediction outcomes and the recommendation are
written as empty, and their emptiness is written into the payload as counts so the artefact proves
it rather than claiming it.

It edits no `.bio` file, changes no number, runs no program and sends no model request.
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.provenance import rule_evidence_tier as ret  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "rule_evidence_tier_proposal_registration"
CENSUS = "data/results/rule_number_sources_census.json"
PROGRAM = "data/demo/gastrulation.bio"

OWN_CODE = (
    "genomeos/provenance/rule_evidence_tier.py",
    "scripts/rule_evidence_tier_register.py",
)


def main() -> None:
    reg = ret.registration()
    empty = {
        k: reg[k] for k in ("adoption_cost_entries", "prediction_outcomes_entries", "recommendation_length")
    }
    if any(v != 0 for v in empty.values()):
        raise SystemExit(f"a registration may not carry a finding: {empty}")
    payload: dict[str, Any] = {
        "result": RESULT,
        "date": date.today().isoformat(),
        "kind": (
            "a registered PROPOSAL, not a measurement. It proposes a mark; it adopts nothing, "
            "recommends nothing and decides nothing. Adopting it changes every hand-authored "
            "program, so the decision is Albert's."
        ),
        **reg,
    }
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": CENSUS,
                "version": "result 4d4003d, registered at 4392d24; cited, and its counts are not "
                "restated in this artefact",
            },
            {
                "accession": PROGRAM,
                "version": "the working tree at the code revision stamped below; read only, never written to",
            },
        ],
        "inputs": [
            mf.input_entry(Path(CENSUS), partition="the census this proposal answers"),
            mf.input_entry(Path(PROGRAM), partition="the program whose rule numbers the census covers"),
        ],
        "assembly": "n/a: no genomic sequence or coordinate is read.",
        "coordinates": "n/a: the proposal addresses declarations by rule id and field name.",
        "method": (
            "The design is written down before it is implemented. The four tiers and their "
            "definitions are IMPORTED from genomeos.provenance.rule_number_sources as the same Python "
            "objects (TIERS is CASCADE, TIER_DEFINITIONS is CLASSES), so the language cannot grow a "
            "divergent second copy of the census's vocabulary; the suite asserts the identity with "
            "`is`. A fifth axis value, `not_assessed`, is the default and is not a tier: it marks a "
            "slot nobody has searched, which is a different state of the world from a slot searched "
            "and found empty, and defaulting to either neighbour would assert something no artefact "
            "supports. Six predictions are registered with the method that will test each one, "
            "including two about whether an EXISTING module already answers the question, so the "
            "proposal can be found redundant by its own tests."
        ),
        "parameters": {
            "axis_values": list(ret.AXIS_VALUES),
            "tiers": list(ret.TIERS),
            "default": ret.NOT_ASSESSED,
            "tiered_fields": list(ret.TIERED_FIELDS),
            "predictions_registered": len(ret.PREDICTIONS),
            "bio_files_written": 0,
            "numbers_changed": 0,
            "model_requests": 0,
            "money_spent": "none; no paid API is called and no network request is made.",
        },
        "exclusions": list(ret.EXCLUSIONS),
        "partitions": {
            "tier": "one of the four census classes, per rule number",
            ret.NOT_ASSESSED: "the absence of a tier, per rule number; the default",
            "interaction_standing": "one of the census's second-axis values, per rule",
        },
        "code_cleanliness": mf.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
    }
    path = save_result(RESULT, payload)
    print(f"{RESULT}: {path}")
    print(f"  axis values: {len(ret.AXIS_VALUES)} ({len(ret.TIERS)} tiers + {ret.NOT_ASSESSED})")
    print(f"  default: {ret.NOT_ASSESSED}")
    print(f"  predictions registered: {len(ret.PREDICTIONS)}")
    print("  adoption cost, prediction outcomes and recommendation: all registered EMPTY")
    print("  0 .bio files written, 0 numbers changed, 0 model requests")


if __name__ == "__main__":
    main()
