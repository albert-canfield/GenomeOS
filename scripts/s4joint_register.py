# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register M2d before any joint stage is computed on a real pair (data/results/s4joint_registration.json).

    uv run python scripts/s4joint_register.py

Every binding word is imported from `genomeos.attribution.s4joint`, the module that applies it, and
`scripts/s4joint_run.py` refuses to run if any registered field differs from that module, if the frozen
code or the stage code it calls moved, or if a pinned input's sha256 moved. This script reads the
population's identifiers (chrom, start, end, gene) and no label, checks them against signaltrace's pinned
list, checks combine2's registration against combine2.py, quotes signaltrace's committed comparator
readings, and computes no stage.

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
from genomeos.attribution import s4joint as sj  # noqa: E402
from genomeos.attribution import signaltrace as st  # noqa: E402
from genomeos.results import save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import s4joint_run as run  # noqa: E402
import signaltrace_run as sr  # noqa: E402

RESULT = "s4joint_registration"
WILL_WRITE = "data/results/s4joint_result.json"

EXPLORATION_PRECEDED_THIS = (
    "stated so that nothing here is read as blind to figures this lane had seen. Before writing it the"
    " lane read: signaltrace's module, scripts and committed result (every stage's kept counts of the 1,918"
    " pairs and the 9,237 fit rows, gains, intervals, outcomes and drop reasons: s2_target 39 kept, 35"
    " regulated, +0.0792 [0.0282, 0.1562], WEAKENED; s3_compiled_value +0.0704 [0.0225, 0.1412], WEAKENED;"
    " s4_bridge NONE_KEPT, all 39 not_identifiable and all 377 reached rules not_identifiable; s4a 6 kept,"
    " +0.0260 [0.0027, 0.0711]; chain LOST_AT s4_bridge); combine2's module, census script and run script,"
    " its census (17,983 multi-element genes eligible for a fit; runnable 4,937 under the mean and 2,687"
    " under competition, the rest mostly basal_not_positive; median summed single-removal drop 1.85; no"
    " gene with rules in several programs) and its registered RNA reading (both forms worse than the"
    " no-rule baseline). Which of the reached genes the joint fit makes runnable was NOT looked at: no joint"
    " stage was computed on any real pair before this file was committed; the joint stages were exercised"
    " only on the synthetic program of tests/test_s4joint.py"
)
DEPARTURES = (
    "s4a_bridge_as_compiled is not re-run: the brief names s4 and the chain after it; signaltrace's s4a"
    " reading stands as committed",
    "combine2's POPULATION (a) asks that every rule of a gene be predicted; here the program is read with"
    " its experimental rule lines removed (signaltrace's PREDICTED_ONLY), so a gene whose program also"
    " carries a CRISPRi line is fitted on its predicted rules. Excluding it instead would select pairs by"
    " the CRISPRi screens' own layer, the held-out file included",
    "combine2's POPULATION (e) (a prediction that varies across the thirteen cells) and (f) (an RNA row)"
    " are not applied: no RNA is read and the stages are read in K562 only; combine2's 13-cell runs are"
    " kept as the HARNESS, not as a filter",
    "the joint fit is applied whatever a gene's rule count; combine2 ran genes with two or more rules, and"
    " signaltrace's committed result shows every reached rule's gene has two or more",
    "the comparators are recomputed by signaltrace's own code and must reproduce its committed readings"
    " exactly, so the joint stages are read on the same footing as M1c's",
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
    c2_bad = sj.combine2_mismatch(json.loads(run.C2_REGISTRATION.read_text()))
    if c2_bad:
        raise SystemExit(f"refused: combine2's registration and combine2.py disagree on {c2_bad}")
    st_result = json.loads(run.ST_RESULT.read_text())
    chroms = sr.chromosomes(training + pop)
    entries = run.inputs(chroms, training, pop)
    return {
        "lane": sj.LANE,
        "date": date.today().isoformat(),
        "roadmap_item": (
            "§5 item 7e, M2d: the bridge stage re-run with the joint fit on the pinned CRISPRi pairs"
        ),
        **sj.registered_fields(),
        "population_record": {
            "sha256": st_reg["population_sha256"],
            "pairs": len(ids),
            "ids_from": "data/results/signaltrace_registration.json (population_ids)",
            "ids_from_sha256": file_sha256(run.ST_REGISTRATION),
            "chromosomes": sorted({p.chrom for p in pop}),
        },
        "population": sj.POPULATION,
        "comparator_committed_verbatim": {
            "source": "data/results/signaltrace_result.json",
            "source_sha256": file_sha256(run.ST_RESULT),
            "anchor": st_result["anchor"]["committed"],
            **{s: st_result["readings_verbatim"][s] for s in ("s2_target", "s3_compiled_value", "s4_bridge")},
            "chain_outcome": st_result["outcome"],
            "chain_reading": st_result["reading"],
        },
        "combine2_registration": {
            "path": "data/results/combine2_registration.json",
            "sha256": file_sha256(run.C2_REGISTRATION),
            "fields_equal_to_combine2_py": list(sj.COMBINE2_BOUND),
        },
        "chromosomes_read": chroms,
        "exploration_preceded_this": EXPLORATION_PRECEDED_THIS,
        "departures": list(DEPARTURES),
        "frozen_code": {
            "commit": git_head(),
            "sha256": sr.code_sha256(run.FROZEN),
            "tests": "tests/test_s4joint.py (synthetic programs; not on the run's path, not pinned)",
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
