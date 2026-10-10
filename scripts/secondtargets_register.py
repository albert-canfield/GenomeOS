# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register the second-targets test before any readout of Gasperini 2019 is read.

    uv run python scripts/secondtargets_register.py

Writes data/results/secondtargets_registration.json. Every binding word is imported from
`genomeos.attribution.secondtargets`, the module that applies it, and `scripts/secondtargets_run.py`
refuses to run if any registered field differs from that module, if the frozen code or the code it calls
moved, or if a pinned input's sha256 moved. The registration carries the label-free census (the screen's
tests, the independence exclusion, the shortlist, the matching) and the ids of the matched shortlist
tests; no readout is opened (the readout table is not on this machine).

No request is sent, no key is read, no outcome is computed. 0 model requests, no network, no money.
"""

from __future__ import annotations

import subprocess
import sys
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos.attribution import secondtargets as s2  # noqa: E402
from genomeos.results import save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import secondtargets_run as run  # noqa: E402

RESULT = "secondtargets_registration"
WILL_WRITE = "data/results/secondtargets_result.json"
FROZEN = run.FROZEN

INDEPENDENCE_CHECK = (
    "the CRISPRi benchmark the project uses DOES include Gasperini 2019: its training file"
    " (EPCrisprBenchmark training_K562) carries 5,299 Gasperini2019 K562 rows of its 10,356, and its"
    " held-out file carries none. Of the at-scale screen's 128,962 candidate-enhancer tests (guide group"
    " top_two), 187 are on elements selected for a pilot hit and 38 do not lift, leaving 128,737; 6,015 of"
    " those name the same gene at an element overlapping a benchmark row (5,466 through Gasperini2019 rows"
    " alone, 549 through other datasets or several), and are excluded. 122,722 are independent pair by"
    " pair: 64,574 on an element the benchmark tested against another gene and 58,148 on an element it"
    " never tested. The screen's readout for them is not on this machine"
)
EXPLORATION_PRECEDED_THIS = (
    "stated so that nothing here is read as blind to what this lane saw. Before writing this file the lane"
    " read: the at-scale pair table's header (12 columns, no effect, p-value or call column), five"
    " non-targeting rows and three candidate rows; counts of its guide-group suffixes and of its selection"
    " labels (general_group), which include labels ending '-hit' (5,620 rows 'test_5k-hit', 374"
    " 'not_model_selected_K1119hit-hit') that may mark pilot-screen hits by element; they are used only"
    " for the registered exclusion of not_model_selected_*. Counts of the benchmark files' rows by"
    " dataset, cell type and ValidConnection. An exploratory overlap count that also included the"
    " 'second_two' groups (7,492 of 138,057 tests the same gene at an overlapping benchmark element). The"
    " first lines of one ENCODE CRISPRi-FlowFISH element-quantification table (ENCFF005IVC.tsv, to learn"
    " its columns), which displayed one PPIF row's effect size; FlowFISH is not used here. Then, after"
    " MIN_MATCHED, DISTANCE_FACTOR, the hit rule and the outcomes were written, one dry run of this"
    " registration's label-free census (the counts registered below); after it one census key was"
    " renamed and nothing else changed. No Gasperini readout was read: the readout table is not on this"
    " machine. This lane did not read the benchmark's Regulated, Significant or EffectSize columns; the"
    " project has read the benchmark's Gasperini2019 labels in training for weeks, which is why those"
    " tests are excluded"
)
DEPARTURES = (
    "the independent population is what the benchmark's filters left out of the screen (power,"
    " expression, distance), not a random sample of the screen; LIMITATIONS says what that does to power",
    "the hit rule's column names are those of the published readout table as this lane recalls them; the"
    " table is not on this machine to check, so a header that differs stops the run before any value is"
    " read, and a column-name-only amendment is the one change STOP_RULE allows",
    "the nearest-gene baseline is the nearest independent protein-coding gene the screen tested at the"
    " element, which may be a compiled target or a shortlist gene itself; it is read beside the primary"
    " and decides nothing",
    "the shortlist is targetchoice's candidate list as committed; the matching is within the same element,"
    " so the element's own activity is shared by the shortlist test and its controls",
)


#: the numbers INDEPENDENCE_CHECK states, checked against the census before anything is written
STATED = {
    ("pair_table", "candidate_enhancer_tests"): 128962,
    ("pair_table", "selected_for_a_pilot_hit"): 187,
    ("pair_table", "not_lifted"): 38,
    ("candidate_enhancer_tests_read",): 128737,
    ("excluded_same_gene_in_benchmark",): 6015,
    ("excluded_by_benchmark_dataset", "Gasperini2019"): 5466,
    ("independent",): 122722,
    ("independence", "element_other_gene"): 64574,
    ("independence", "element_untested"): 58148,
}
STATED_BENCHMARK = {"Gasperini2019": 5299}


def stated_mismatches(census: dict[str, Any], by_dataset: dict[str, int]) -> list[str]:
    """INDEPENDENCE_CHECK's numbers that the census does not reproduce."""
    out = []
    for path, want in STATED.items():
        got: Any = census
        for k in path:
            got = got.get(k) if isinstance(got, dict) else None
        if got != want:
            out.append(f"{'.'.join(path)}: stated {want}, census {got}")
    out += [
        f"{d}: stated {n}, benchmark {by_dataset.get(d)}"
        for d, n in STATED_BENCHMARK.items()
        if by_dataset.get(d) != n
    ]
    return out


def git_head() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True
    ).stdout.strip()


def payload() -> dict[str, Any]:
    built = run.build()
    census = built["census"]
    entries = run.inputs(built["chroms"])
    ids = built["ids"]
    wrong = stated_mismatches(census, built["benchmark_rows_by_dataset"])
    if wrong:
        raise SystemExit(f"refused: INDEPENDENCE_CHECK states numbers the census does not give: {wrong}")
    return {
        "lane": s2.LANE,
        "date": date.today().isoformat(),
        "roadmap_item": "§3 area I, Next 3: the model's second targets as a shortlist, held against an"
        " independent perturbation",
        **s2.registered_fields(),
        "independence_check": INDEPENDENCE_CHECK,
        "census": census,
        "population_ids": ids,
        "population_sha256": s2.ids_sha256(ids),
        "population_tests": len(ids),
        "shortlist_genome_wide": built["shortlist"],
        "benchmark_rows_by_dataset": built["benchmark_rows_by_dataset"],
        "positive_control_groups": len(built["tss"]),
        "chromosomes_read": built["chroms"],
        "readout": {"path": str(s2.READOUT), "url": s2.READOUT_URL, "held": s2.READOUT.exists()},
        "prediction": (
            "a second or further target of a compiled element is more often a screen hit than the genes the"
            " screen tested at the same element at a matched distance that the shortlist does not name"
        ),
        "baselines": [
            "distance-matched non-shortlist genes at the same element (the primary comparison)",
            "the nearest tested protein-coding gene at the same element (read beside, decides nothing)",
        ],
        "falsifier": (
            "NO_DIFFERENCE or LESS_OFTEN: second targets were then not detected to be measured hits more"
            " often than distance-matched genes the shortlist does not name"
        ),
        "exploration_preceded_this": EXPLORATION_PRECEDED_THIS,
        "departures": list(DEPARTURES),
        "frozen_code": {
            "commit": git_head(),
            "sha256": run.code_sha256(FROZEN),
            "tests": "tests/test_secondtargets.py (synthetic tables; not on the run's path, not pinned)",
        },
        "stage_code": {"commit": git_head(), "sha256": run.code_sha256(run.STAGE_CODE)},
        "what_the_run_will_write": WILL_WRITE,
        "alphagenome_requests": 0,
        "readout_values_read": 0,
        "money": "none: no request is sent and no key is read",
        "network": "none: every input is already on this machine",
        "result_manifest": run.manifest(entries, own_entry=__file__),
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
    print(f"matched shortlist tests: {body['population_tests']}, sha256 {body['population_sha256']}")
    print(f"readout held: {body['readout']['held']}; frozen code at {body['frozen_code']['commit']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
