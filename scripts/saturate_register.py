# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register M2e before any multiplicative stage is computed on a held-out pair
(data/results/saturate_registration.json).

    uv run python scripts/saturate_register.py

Every binding word is imported from `genomeos.attribution.saturate`, the module that applies it, and
`scripts/saturate_run.py` refuses to run if any registered field differs from that module, if the frozen
code or the stage code it calls moved, or if a pinned input's sha256 moved. This script reads the
population's identifiers (chrom, start, end, gene) and no label, checks them against signaltrace's pinned
list, quotes signaltrace's and M2d's committed readings, and computes no stage and no census.

No request is sent, no key is read, no outcome is computed. 0 model requests, no network, no money.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import saturate as sa  # noqa: E402
from genomeos.attribution import signaltrace as st  # noqa: E402
from genomeos.results import save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import saturate_run as run  # noqa: E402
import signaltrace_run as sr  # noqa: E402

RESULT = "saturate_registration"
WILL_WRITE = "data/results/saturate_result.json"

EXPLORATION_PRECEDED_THIS = (
    "stated so that nothing here is read as blind to figures this lane had seen. Before writing it the"
    " lane read: signaltrace's module, scripts and committed result (s2_target 39 of 1,918 kept, +0.0792"
    " [0.0282, 0.1562], WEAKENED; s3_compiled_value WEAKENED; s4_bridge NONE_KEPT; chain LOST_AT s4_bridge);"
    " combine2's module (the two forms, classify, fit, run_gene_form); s4joint's module, scripts, tests and"
    " committed result: s4_joint LOST under the mean (1 of 1,918 kept) and NONE_KEPT under competition,"
    " and its census of the 210 reached genes (mean: 10 runnable, 198 basal_not_positive, 2"
    " inhibitor_strength_above_1; competition: 2 runnable, 206 basal_not_positive, 2"
    " inhibitor_strength_above_1; median 28 rules per reached gene). It also read"
    " docs/DESIGN-MINIMAL-CELL.md's account of Zhou et al. 2024 and confirmed the paper's title, authors,"
    " journal and DOI by one web search. Before freezing the code the lane ran one predictor-side dry run"
    " of the joint lookup, the trace and the harness on the TRAINING pairs of chr21 and chr22 only (no"
    " held-out pair, no stage scored, no label read): it reached 6 genes (2 and 4), the multiplicative fit"
    " ran all 6 (7 to 111 rules each) and the harness found no problem. Beyond those 6, which of the"
    " reached genes the multiplicative fit makes runnable, and every multiplicative stage on a held-out"
    " pair, were NOT looked at before this file was committed; otherwise the stages were exercised only on"
    " the synthetic programs of tests/test_saturate.py and tests/test_grn_multiplicative.py"
)
DEPARTURES = (
    "the brief named a saturating composition (Hill or Michaelis-Menten in the summed input, or a"
    " logistic) with its parameters cited or fitted on the training split. WHY_NOT_SATURATING shows that"
    " no concave composition with a nonnegative basal rate can fit a gene whose single-removal drops sum"
    " past its whole level, which is what M2d's 198 and 206 genes are; the registered form is therefore"
    " the multiplicative one (Zhou et al. 2024). It has no free parameter, so nothing is fitted on the"
    " training split either",
    "combine2.classify's no_activating_rule exclusion is kept although the multiplicative runtime could"
    " carry an inhibitor-only gene, so that the admitted genes are M2d's and the census compares with its",
    "s4a_bridge_as_compiled is not re-run, as in M2d: signaltrace's s4a reading stands as committed",
    "as in M2d, the program is read with its experimental rule lines removed (signaltrace's"
    " PREDICTED_ONLY), so a gene whose program also carries a CRISPRi line is fitted on its predicted"
    " rules",
    "combine2's thirteen-cell runs are kept as the HARNESS, not as a filter; no RNA is read and the stages"
    " are read in K562 only",
)


def git_head() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True
    ).stdout.strip()


def file_sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def payload() -> dict[str, Any]:
    training, pop = sr.load_pairs()
    ids = st.population_ids(pop)
    st_reg = json.loads(run.ST_REGISTRATION.read_text())
    if st.ids_sha256(ids) != st_reg["population_sha256"] or ids != st_reg["population_ids"]:
        raise SystemExit("refused: the pairs read differ from signaltrace's pinned population")
    if [f for f, v in st.registered_fields().items() if st_reg.get(f) != v]:
        raise SystemExit("refused: signaltrace's registration and signaltrace.py disagree")
    st_result = json.loads(run.ST_RESULT.read_text())
    m2d = json.loads(run.M2D_RESULT.read_text())
    chroms = sr.chromosomes(training + pop)
    entries = run.inputs(chroms, training, pop)
    return {
        "lane": sa.LANE,
        "date": date.today().isoformat(),
        "roadmap_item": (
            "§5 item 7f, M2e: the bridge stage under a composition that can fit the reached genes, on the"
            " pinned CRISPRi pairs"
        ),
        **sa.registered_fields(),
        "population_record": {
            "sha256": st_reg["population_sha256"],
            "pairs": len(ids),
            "ids_from": "data/results/signaltrace_registration.json (population_ids)",
            "ids_from_sha256": file_sha256(run.ST_REGISTRATION),
            "chromosomes": sorted({p.chrom for p in pop}),
        },
        "comparator_committed_verbatim": {
            "source": "data/results/signaltrace_result.json",
            "source_sha256": file_sha256(run.ST_RESULT),
            "anchor": st_result["anchor"]["committed"],
            **{s: st_result["readings_verbatim"][s] for s in ("s2_target", "s3_compiled_value", "s4_bridge")},
            "chain_outcome": st_result["outcome"],
            "chain_reading": st_result["reading"],
        },
        "m2d_committed_verbatim": {
            "source": "data/results/s4joint_result.json",
            "source_sha256": file_sha256(run.M2D_RESULT),
            "outcome": m2d["outcome"],
            "reading": m2d["reading"],
            "census": {
                f: {
                    k: m2d["descriptives_not_verdicts"]["per_form"][f][k]
                    for k in ("reached_genes", "runnable")
                }
                | {
                    "not_runnable_by_reason": m2d["descriptives_not_verdicts"]["per_form"][f][
                        "not_runnable_by_reason"
                    ]
                }
                for f in ("mean", "competition")
            },
        },
        "chromosomes_read": chroms,
        "exploration_preceded_this": EXPLORATION_PRECEDED_THIS,
        "departures": list(DEPARTURES),
        "frozen_code": {
            "commit": git_head(),
            "sha256": sr.code_sha256(run.FROZEN),
            "tests": "tests/test_saturate.py and tests/test_grn_multiplicative.py (synthetic programs; not on"
            " the run's path, not pinned)",
        },
        "stage_code": {"commit": git_head(), "sha256": sr.code_sha256(run.STAGE_CODE)},
        "what_the_run_will_write": WILL_WRITE,
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
        "network": "none: every input is already on this machine",
        "result_manifest": {
            **run.manifest(entries),
            "code_cleanliness": mf.code_cleanliness(__file__, run.OWN_CODE),
        },
    }


def differs_from_head(paths: tuple[str, ...]) -> list[str]:
    """The paths whose working copy differs from HEAD: pinned code must be committed code."""
    out = subprocess.run(["git", "diff", "--name-only", "HEAD", "--", *paths], capture_output=True, text=True)
    return [ln for ln in out.stdout.splitlines() if ln]


def main() -> int:
    differs = differs_from_head((*run.FROZEN, *run.STAGE_CODE))
    if differs:
        raise SystemExit(f"refused: pinned code differs from HEAD: {differs}")
    body = payload()
    clean = body["result_manifest"]["code_cleanliness"]
    if not clean.get("own_code_is_committed"):
        raise SystemExit(f"refused: the frozen code is not committed: {clean.get('own_uncommitted_code')}")
    path = save_result(RESULT, body)
    print(f"registered: {path}")
    rec = body["population_record"]
    print(f"population: {rec['pairs']} pairs, sha256 {rec['sha256']}")
    print(f"frozen code at {body['frozen_code']['commit']}; the run will write {WILL_WRITE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
