# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register M1c before any stage beyond the anchor is computed (data/results/signaltrace_registration.json).

    uv run python scripts/signaltrace_register.py

Every binding word is imported from `genomeos.attribution.signaltrace`, the module that applies it, and
`scripts/signaltrace_run.py` refuses to run if any registered field differs from that module, if the
frozen code or the stage code it calls moved, or if a pinned input's sha256 moved. This script reads
the population's identifiers (chrom, start, end, gene) and no label, and computes no stage.

No request is sent, no key is read, no outcome is computed. 0 model requests, no network, no money.
"""

from __future__ import annotations

import subprocess
import sys
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import signaltrace as st  # noqa: E402
from genomeos.results import save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import signaltrace_run as run  # noqa: E402

RESULT = "signaltrace_registration"
WILL_WRITE = "data/results/signaltrace_result.json"
FROZEN = run.FROZEN

EXPLORATION_PRECEDED_THIS = (
    "stated so that nothing here is read as blind to figures this lane had seen. Before writing it the"
    " lane read: the anchor's committed figures (crispri_published_v2.json arm 1: 1,918 pairs, 118"
    " positives, 174 uncovered, gain 0.1361 [0.0749, 0.2258]; crispri_benchmark_v2.json: K562 1,744"
    " covered pairs, 114 positives, gain 0.1407 [0.0786, 0.2321]) and their code (crispri.score_published,"
    " coverage_arms, gain_interval); the designs and registered readings of M1 (runrules), M2b"
    " (readergate), M1b (regulators) and M3 (closure_cell) as the work board and their module docstrings"
    " state them; the compiled program format on chr21 (a parse: 5,176 rules, 2 of them experimental; one"
    " element block and its all-element table row, whose coordinates are the same numbers); the compiled"
    " header's quotation of data/results/argmaxcell.json (the compiled cell is the measured cell on 612 of"
    " 3,228 K562 elements); that the compiled programs carry CRISPRi measured rules including held-out"
    " ones (two rule lines of chr21, read as text); and the population's identifiers only (1,918 unique"
    " ids on 22 chromosomes, 424 genes, none shared with the 10,356 training ids, 2,113 training genes)."
    " No stage beyond the anchor was computed on any real pair before this file was committed; the"
    " stages were exercised only on the synthetic program of tests/test_signaltrace.py"
)
DEPARTURES = (
    "from the brief's six stages: 'target choice' is split into s2_target (which pairs the compiled"
    " programme keeps, the value held at the anchor's) and s3_compiled_value (the value the compiled rule"
    " states), so that each stage changes one thing and a loss can be placed on one of them",
    "s4a_bridge_as_compiled is read beside the chain: M1 and M1b ran each rule in its compiled cell, which"
    " in K562 keeps only rules whose compiled cell is K562, while the chain follows the reader-cell"
    " selection M2b ran, where every predicted rule is offered in K562",
    "the run (s7_run) follows the reader gate (s5_reader), not the regulator stage: "
    + st.WHY_THE_RUN_FOLLOWS_THE_READER,
    "'loses the signal' is decided by intervals, not by a fixed margin: DECISION",
)


def git_head() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True
    ).stdout.strip()


def payload() -> dict[str, Any]:
    training, pop = run.load_pairs()
    ids = st.population_ids(pop)
    chroms = run.chromosomes(training + pop)
    entries = run.inputs(chroms, training, pop)
    return {
        "lane": st.LANE,
        "date": date.today().isoformat(),
        "roadmap_item": "§5 item 7c, M1c: where the compiled layer loses its signal",
        **st.registered_fields(),
        "population_ids": ids,
        "population_sha256": st.ids_sha256(ids),
        "population_pairs": len(ids),
        "population_chromosomes": sorted({p.chrom for p in pop}),
        "chromosomes_read": chroms,
        "exploration_preceded_this": EXPLORATION_PRECEDED_THIS,
        "departures": list(DEPARTURES),
        "frozen_code": {
            "commit": git_head(),
            "sha256": run.code_sha256(FROZEN),
            "tests": "tests/test_signaltrace.py (synthetic programs; not on the run's path, not pinned)",
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


def unchanged_from_head(paths: tuple[str, ...]) -> list[str]:
    """The paths whose working copy differs from HEAD: pinned code must be committed code."""
    out = subprocess.run(["git", "diff", "--name-only", "HEAD", "--", *paths], capture_output=True, text=True)
    return [ln for ln in out.stdout.splitlines() if ln]


def main() -> int:
    differs = unchanged_from_head((*FROZEN, *run.STAGE_CODE))
    if differs:
        raise SystemExit(f"refused: pinned code differs from HEAD: {differs}")
    body = payload()
    clean = body["result_manifest"]["code_cleanliness"]
    if not clean.get("own_code_is_committed"):
        raise SystemExit(f"refused: the frozen code is not committed: {clean.get('own_uncommitted_code')}")
    path = save_result(RESULT, body)
    print(f"registered: {path}")
    print(f"population: {body['population_pairs']} pairs, sha256 {body['population_sha256']}")
    print(f"frozen code at {body['frozen_code']['commit']}; the run will write {WILL_WRITE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
