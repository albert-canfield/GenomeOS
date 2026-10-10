# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register lane-strengthform before it is run (data/results/strengthform_registration.json).

    uv run python scripts/strengthform_register.py

Every binding word comes from `genomeos.attribution.strengthform.registered_fields`, and
scripts/strengthform_run.py refuses to run if any registered field differs from that module or the
table's sha256 differs from the pinned one. The census below is computed with
`read_table(read_outcome=False)`, which never touches the activity column of an enhancer-enhancer
row, so the population is counted and committed before any pair outcome exists in this lane.

No request is sent, no key is read, no outcome is computed.
"""

from __future__ import annotations

import sys
from collections import Counter
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import combine2 as c2  # noqa: E402
from genomeos.attribution import strengthform as sf  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "strengthform_registration"
WILL_WRITE = "data/results/strengthform_result.json"
OWN_CODE = (
    "genomeos/attribution/strengthform.py",
    "scripts/strengthform_fetch.py",
    "scripts/strengthform_register.py",
    "scripts/strengthform_run.py",
    "tests/test_strengthform.py",
)

EXPLORATION_PRECEDED_THIS = (
    "stated so that nothing here is read as a blind prediction. Before fixing the form the lane read:"
    " lane-combine's census, registration and result (`a6fc5c4`, `7f3e597`, `93caf61`) and their"
    " committed reading; runtime/grn.py's four fixed rules and its two opt-in forms; the 7d, 7f and 7h"
    " rows of ROADMAP section 5; the Methods and the DATA AVAILABILITY of Martinez-Ara et al. 2024 and"
    " that paper's own one-sentence summary of its finding (that enhancer pairs combine near-additively"
    " and that weaker promoters can be supra-additive); and from the table itself the design columns,"
    " the eight promoter baselines, the control-control rows and the single-enhancer boost indices,"
    " which are this test's INPUTS. It read no activity_all value of any enhancer-enhancer row, and"
    " none of the authors' own expectation or cooperativity columns (`columns_never_read`). The form's"
    " shape is Bothma et al. 2015's and its one constant is 7d's committed THETA; neither was chosen"
    " after seeing an outcome of this dataset."
)

DATA_DECISION = (
    "lane-combine's own five Bothma pairs cannot test this form independently: that lane read their"
    " outcomes before writing its rule and recorded that it had. A held-out population exists that it"
    " never read, Martinez-Ara et al. 2024 (GSE240586), whose combinations are designed singles and"
    " pairs with one baseline per promoter in one assay and one cell context, so the form is tested"
    " there instead, with its shape and its one constant fixed beforehand from two other publications."
)


def census() -> dict[str, Any]:
    """The population, counted with no pair outcome read."""
    table = sf.read_table(read_outcome=False)
    combos = table["combinations"]
    clusters = {c.cluster for c in combos}
    per_promoter = Counter(c.promoter for c in combos)
    return {
        "counts": table["counts"],
        "informative_combinations": len(combos),
        "unordered_enhancer_pairs": len(clusters),
        "promoters": sorted(per_promoter),
        "combinations_per_promoter": dict(sorted(per_promoter.items())),
        "promoter_baselines": {k: table["baselines"][k] for k in sorted(table["baselines"])},
        "promoter_largest_single_boost": {
            k: table["promoter_maxima"][k] for k in sorted(table["promoter_maxima"])
        },
        "singles_read": len(table["singles"]),
        "outcome_read": False,
    }


def inputs() -> list[dict[str, Any]]:
    return [mf.input_entry(sf.TABLE, partition="held out: never read by any earlier lane")]


def main() -> int:
    entries = inputs()
    out: dict[str, Any] = {
        "lane": sf.LANE,
        "date": date.today().isoformat(),
        "roadmap_item": sf.ROADMAP_ITEM,
        **sf.registered_fields(),
        "theta_source": (
            "attribution/combine2.THETA, declared by 7d on 2026-10-10 and committed in `bc9fc7a`;"
            f" this lane reads it as {c2.THETA} and refuses to run if it differs"
        ),
        "data_decision": DATA_DECISION,
        "exploration_preceded_this": EXPLORATION_PRECEDED_THIS,
        "census": census(),
        "engine_unchanged": (
            "runtime/grn.py is an Apache-2.0 engine file and is not edited by this lane: the form is"
            " its existing opt-in rule `competition`, and tests/test_strengthform.py pins the file's"
            " licence header and every existing form's identity"
        ),
        "stop_list_checked": (
            "ROADMAP section 5, 'The stop list of 2026-09-27, reinstated': no AlphaGenome request, no"
            " deletion request on the 882 blocks, no new biosample, nothing built in area J, nothing"
            " on the node model, no fate ceiling, no BioForge calibration"
        ),
        "what_the_run_will_write": WILL_WRITE,
        "the_code_that_runs_it": (
            "genomeos/attribution/strengthform.py and scripts/strengthform_run.py, committed before"
            " this registration, with the committed runtime/grn.py they import"
        ),
        "alphagenome_requests": 0,
        "money": "none: one public GEO file already fetched; no request is sent and no key is read",
    }
    out["pinned_sha256"] = {x["path"]: x["sha256"] for x in entries}
    manifest = {
        "sources": [
            {"accession": sf.SOURCE["accession"], "version": "pinned by sha256"},
            {"accession": sf.SOURCE["publication"], "version": "2024"},
        ],
        "inputs": entries,
        "input_count": len(entries),
        "assembly": "n/a: a barcoded episomal reporter library, not genomic coordinates",
        "coordinates": {"base": "n/a", "interval": "n/a"},
        "parameters": {
            "form": sf.FORM,
            "theta": sf.THETA,
            "half_saturation": sf.HALF_SATURATION,
            "comparators": list(sf.COMPARATORS),
            "bootstrap_draws": sf.BOOTSTRAP_DRAWS,
            "bootstrap_seed": sf.BOOTSTRAP_SEED,
        },
        "exclusions": list(sf.INFORMATIVE),
        "partitions": "held out: this dataset has never been read by any earlier lane of the project",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }
    print("saved", save_result(RESULT, out, manifest=manifest))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
