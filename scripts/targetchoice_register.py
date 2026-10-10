# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register the multi-target compiled choice before it is read on any held-out pair.

    uv run python scripts/targetchoice_register.py

Writes data/results/targetchoice_registration.json. Every binding word is imported from
`genomeos.attribution.targetchoice`, the module that applies it, and `scripts/targetchoice_run.py`
refuses to run if any registered field differs from that module, if the frozen code or the code it calls
moved, or if a pinned input's sha256 moved. The population is M1c's, taken from M1c's committed
registration by its ids; the held-out file is hashed and not parsed here, and no choice is computed.

No request is sent, no key is read, no outcome is computed. 0 model requests, no network, no money.
"""

from __future__ import annotations

import json
import subprocess
import sys
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri  # noqa: E402
from genomeos.attribution import targetchoice as tc  # noqa: E402
from genomeos.results import save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import targetchoice_run as run  # noqa: E402

RESULT = "targetchoice_registration"
WILL_WRITE = "data/results/targetchoice_result.json"
FROZEN = run.FROZEN

EXPLORATION_PRECEDED_THIS = (
    "stated so that nothing here is read as blind to figures this lane had seen. Before writing it the"
    " lane read: M1c's committed figures (signaltrace_result.json: the anchor's 1,918 pairs, 118"
    " positives, top_target kept 40 with 36 regulated; s2_target kept 39 with 35 regulated, gain +0.0792"
    " [0.0282, 0.1562], against the anchor -0.0569 [-0.0970, -0.0248], dropped because 'the overlapping"
    " elements name another gene' 1,209 and 'no compiled element overlaps' 670; 9,237 fit rows, 245 kept"
    " by s2_target with 188 regulated) and M1c's code; the shape of the response cache (one chr21"
    " record); a parse of chr21's compiled program (5,174 predicted rules, one per element); the training"
    " census (data/results/targetchoice_census.json, 3d85122), and before it a dry run of the same census"
    " on the chr21 and chr22 training pairs only; and, after TOP_K was fixed from the census, a dry run of"
    " the run's code path with a TRAINING stand-in for the population (every fifth training pair scored,"
    " the other four fifths as fit rows, the anchor and comparator checks bypassed, the held-out file not"
    " parsed, written to a scratch directory), whose reading on that training stand-in was NO_BETTER,"
    " +0.0576 [95% CI -0.0177, +0.1163] for the multi-target choice against the one-target choice (72"
    " pairs kept against 43). No parameter or rule was changed after it. This lane computed no choice on"
    " any held-out pair before this file was committed (the one-target choice's held-out figures above are"
    " M1c's committed ones); the choices were exercised only on the training split and on the synthetic"
    " programs of tests/test_targetchoice.py"
)
DEPARTURES = (
    "the primary multi-target choice takes the compiler's own bar (enhancer_target.MIN_EFFECT) and no"
    " parameter, as the brief preferred; the top-k choice is secondary, its k fixed on the training split"
    " by CENSUS_RULE (in-sample average precision on the fit rows), and decides nothing",
    "the training split is the one the project's CRISPRi benchmark fits on, EPCrisprBenchmark"
    " training_K562: 9,237 covered K562 pairs, 451 regulated, the anchor's and M1c's fit rows",
    "the value column of every choice is the anchor's deletion_drop, as in M1c's s2_target, so the"
    " comparison changes only which pairs are kept; the value a compiled rule would state for a further"
    " target is not read",
    "candidate targets are protein-coding genes, the kind of target the compiled programme names"
    " (compile.compile_chromosome reads predicted_coding)",
)


def git_head() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True
    ).stdout.strip()


def payload() -> dict[str, Any]:
    m1c = json.loads(run.M1C_REGISTRATION.read_text())
    ids = list(m1c["population_ids"])
    training = crispri.load(crispri.TRAINING)
    chroms = run.chromosomes_of(training, ids)
    entries = run.inputs(chroms, len(training), len(ids))
    census = json.loads(run.CENSUS.read_text())
    if census["top_k"] != tc.TOP_K:
        raise SystemExit(f"refused: the census names k = {census['top_k']}, targetchoice.TOP_K is {tc.TOP_K}")
    return {
        "lane": tc.LANE,
        "date": date.today().isoformat(),
        "roadmap_item": "§5 after item 7c (M1c): a multi-target compiled choice on the held-out pairs",
        **tc.registered_fields(),
        "population_ids": ids,
        "population_sha256": m1c["population_sha256"],
        "population_pairs": len(ids),
        "population_from": "data/results/signaltrace_registration.json (M1c, 9637584)",
        "chromosomes_read": chroms,
        "choices_read": list(tc.read_choices(tc.TOP_K)),
        "falsifier": (
            "the primary's NO_BETTER, WORSE or NONE_ADDED: the multi-target choice is then no better than"
            " the one-target choice on these pairs"
        ),
        "training_census": {
            "result": tc.CENSUS,
            "top_k": census["top_k"],
            "secondary": census["secondary"],
            "fit_rows": census["fit_rows"],
            "one_target_in_sample": {
                k: census["one_target_in_sample"][k] for k in ("kept", "kept_regulated")
            },
            "multi_target_in_sample": {
                k: census["by_k_in_sample"]["all"][k] for k in ("kept", "kept_regulated")
            },
            "compiled_elements_reached_by_a_fit_row": census["compiled_elements_reached_by_a_fit_row"],
        },
        "exploration_preceded_this": EXPLORATION_PRECEDED_THIS,
        "departures": list(DEPARTURES),
        "frozen_code": {
            "commit": git_head(),
            "sha256": run.code_sha256(FROZEN),
            "tests": "tests/test_targetchoice.py (synthetic programs; not on the run's path, not pinned)",
        },
        "stage_code": {"commit": git_head(), "sha256": run.code_sha256(run.STAGE_CODE)},
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
    differs = differs_from_head((*FROZEN, *run.STAGE_CODE))
    if differs:
        raise SystemExit(f"refused: pinned code differs from HEAD: {differs}")
    body = payload()
    clean = body["result_manifest"]["code_cleanliness"]
    if not clean.get("own_code_is_committed"):
        raise SystemExit(f"refused: the frozen code is not committed: {clean.get('own_uncommitted_code')}")
    path = save_result(RESULT, body)
    print(f"registered: {path}")
    print(f"population: {body['population_pairs']} pairs, sha256 {body['population_sha256']}")
    print(f"choices: {body['choices_read']}, TOP_K {body['top_k']}")
    print(f"frozen code at {body['frozen_code']['commit']}; the run will write {WILL_WRITE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
