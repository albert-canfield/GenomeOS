# SPDX-License-Identifier: AGPL-3.0-or-later
"""The base rate the 0.3252 has to be read against: the same rules, read in every other biosample.

    uv run --frozen python scripts/context_evidence_baserate.py [--chroms chr21]
    uv run --frozen python scripts/context_evidence_baserate.py --result context_evidence_baserate_v2 \
        --supersedes

`--result` changes the name written and nothing else; `--supersedes` records which committed result the
run is written beside. The committed data/results/context_evidence_baserate.json is kept byte for byte:
its sha256 is a declared input of four other committed results.

`data/results/context_evidence.json` reports that 26,551 of the 81,635 assessable compiled rules are
`open_in_reader`, 0.3252 of that population. That share means nothing on its own. The elements are
ENCODE registry cCREs, built from DNase over many biosamples, so **any** biosample opens some share of
them whatever cell a rule names. Without the base rate, 0.3252 could be a finding or could be
arithmetic.

So this reads **exactly the same assessable rules** in each of the thirteen reader biosamples, not
only in the one the compiler assigned, and reports the assigned-cell share beside the mean over the
twelve others. If the assigned share sits near that mean, the cell the compiler assigns - the most
extreme AlphaGenome track - carries little information about where the element is actually open. If it
sits well above, the assignment carries real information.

**This comparison is descriptive and was not registered.** `context_evidence_registration.json` fixed
the states, the openness call, the mapping and the denominator, and none of them is touched here: the
same `state_for` is used, the same threshold, the same table, the same rule set. The only thing added
is the same reading taken in twelve more biosamples. The tolerance the reading uses was chosen before
the numbers were computed and is stated in the result.

`open_in_reader` remains a consistency check between two readings of the same ENCODE chromatin and is
never independent evidence that a rule is right. This comparison tests whether the cell assignment is
**informative**, not whether it is **correct**.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import context_evidence as ce  # noqa: E402
from genomeos.genome import reader  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

RESULT = "context_evidence_baserate"
CHROMS = tuple(f"chr{c}" for c in [*range(1, 23), "X", "Y"])
ENTRY = Path(__file__).resolve()

#: Chosen before any share was computed: how far the assigned cell's share must sit from the mean over
#: the other twelve biosamples before the difference is read as a difference at all. Not registered,
#: because this whole comparison is descriptive; stated so that the reading cannot be softened after
#: the numbers are in.
TOLERANCE = 0.02

NEAR = (
    "the cell the compiler assigns carries little information about where the element is actually "
    "open: the same rules are open in the other reader biosamples at about the same rate"
)
ABOVE = (
    "the cell the compiler assigns carries information about where the element is actually open: the "
    "same rules are open in that cell more often than in the other reader biosamples"
)
BELOW = (
    "the cell the compiler assigns is open LESS often than the other reader biosamples over the same "
    "rules, which would mean the assignment points away from the open chromatin"
)

OWN_CODE = (
    "genomeos/attribution/context_evidence.py",
    "genomeos/attribution/compile.py",
    "genomeos/lang/grammar.py",
    "genomeos/lang/parser.py",
    "genomeos/ir/model.py",
    "scripts/context_evidence_register.py",
    "scripts/context_evidence_census.py",
    "scripts/context_evidence_baserate.py",
    "tests/test_context_evidence.py",
    "docs/BIOLANG-GRAMMAR.md",
    "docs/ATTRIBUTION.md",
)


#: Why this base rate is re-recorded under a name of its own rather than over the committed file. The
#: committed manifest declares the 312 cached DNase peak sets as one group under a label with no member
#: list, so a rebuild in a second environment has no path to open or hash: `files_entry` did not name the
#: members of a group until 2026-10-02. The committed bytes cannot be replaced to fix it: the sha256 of
#: data/results/context_evidence_baserate.json is a declared input of four other committed results
#: (not_open_profile, not_open_profile_chr21, repress2_registration, repress2_registration_amendment),
#: which new bytes under the same name would make unrebuildable.
WHY_A_NEW_NAME = (
    "re-recorded under a new name on 2026-10-02 so that every group of input files names its members, "
    "each with its own sha256 and byte count, which is what a rebuild in a second environment needs to "
    "open and hash them; no figure of the run differs. The committed result is kept unchanged because "
    "its sha256 is a declared input of four other committed results"
)
#: The committed results that declare this one's sha256 as an input, which is why it is kept byte for byte.
DECLARED_AS_AN_INPUT_BY = (
    "data/results/not_open_profile.json",
    "data/results/not_open_profile_chr21.json",
    "data/results/repress2_registration.json",
    "data/results/repress2_registration_amendment.json",
)


#: What naming every member of every group does NOT fix, recorded in each new file because a rebuild
#: that passes without it would read as reproduction. The 24 `enhancer_targets_all_chr*.json` files this
#: manifest declares are 1.6 KB each: they are pointers. `genomeos.attribution.targets.run_elements`
#: reads the path in their `elements_where` field and opens it, which is
#: data/knowledge/alphagenome/all_elements/<chrom>.json, about 1.4 GB over the 24 chromosomes, and no
#: manifest declares those tables. A rebuild in a clean worktree finds them anyway, because data/knowledge
#: is linked read-only, so it PASSES without having hashed the largest thing the run read. Raised to the
#: supervisor as an item of its own on 2026-10-02; not closed here, because declaring them changes the
#: input list of several results materially.
ELEMENT_TABLES_UNDECLARED = (
    "the 24 compiler_inputs:enhancer_targets_all_chr*.json files declared above are 1.6 KB pointers: "
    "genomeos.attribution.targets.run_elements opens the path in each one's `elements_where` field, "
    "data/knowledge/alphagenome/all_elements/<chrom>.json, about 1.4 GB in all, and no manifest here "
    "declares those tables. Naming the members of every group does not close that gap, and a rebuild in "
    "a clean worktree passes without hashing them because data/knowledge is linked read-only. Raised to "
    "the supervisor as its own item on 2026-10-02"
)


def supersedes() -> dict[str, Any]:
    """The committed result this run is written beside, never over, with the bytes it is kept at."""
    p = RESULTS_DIR / f"{RESULT}.json"
    entry = mf.input_entry(p, partition=None)
    old = json.loads(p.read_text())
    return {
        "file": p.as_posix(),
        "date": old.get("date"),
        "sha256": entry["sha256"],
        "bytes": entry["bytes"],
        "kept": "unchanged; this run is written beside it under a new name, not over it",
        "why": WHY_A_NEW_NAME,
        "declared_as_an_input_by": list(DECLARED_AS_AN_INPUT_BY),
        "element_tables_behind_the_pointer_files": ELEMENT_TABLES_UNDECLARED,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", default=",".join(CHROMS))
    ap.add_argument("--result", default=RESULT)
    ap.add_argument(
        "--supersedes",
        action="store_true",
        help=(
            f"record, beside the manifest, that this run is written beside the committed {RESULT} "
            f"rather than over it, with that file's date and sha256"
        ),
    )
    args = ap.parse_args()
    chroms = [c for c in args.chroms.split(",") if c]

    started = time.time()
    table = ce.mapping()
    biosamples = tuple(ce.READER_TERMS)

    population = 0
    assigned: Counter[str] = Counter()  # the state in the rule's own assigned biosample
    per_biosample: dict[str, Counter[str]] = defaultdict(Counter)  # over the same rules, every cell
    # per assigned biosample: its own rules, its own openness, and the others' over those same rules
    per_assignment: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    others_open = 0
    others_assessable = 0
    others_outside = 0

    for chrom in chroms:
        readers = ce.Readers()
        for label, start, end in ce.rule_loci(chrom):
            mine, _term, _why = ce.reader_for(label, table)
            if mine is None:
                continue
            span = readers.span(mine, chrom)
            if span is None or not (end > span[0] and start < span[1]):
                continue  # not assessable in its own cell, so not in this population
            population += 1
            row = per_assignment[mine]
            row["rules"] += 1
            for b in biosamples:
                sb = readers.span(b, chrom)
                if sb is None or not (end > sb[0] and start < sb[1]):
                    state = ce.REASON_OUTSIDE_SPAN
                else:
                    idx = readers.index(b, chrom)
                    assert idx is not None
                    state = ce.STATE_OPEN if idx.overlapping(start, end) else ce.STATE_NOT_OPEN
                per_biosample[b][state] += 1
                if b == mine:
                    assigned[state] += 1
                    row["own_" + state] += 1
                    continue
                if state == ce.REASON_OUTSIDE_SPAN:
                    others_outside += 1
                    row["others_outside_span"] += 1
                    continue
                others_assessable += 1
                row["others_assessable"] += 1
                if state == ce.STATE_OPEN:
                    others_open += 1
                    row["others_open"] += 1
        print(f"  {chrom}: population now {population:,}", flush=True)

    def share(num: int, den: int) -> float | None:
        return None if not den else round(num / den, 4)

    assigned_share = share(assigned[ce.STATE_OPEN], population)
    others_mean = share(others_open, others_assessable)
    difference = (
        None if assigned_share is None or others_mean is None else round(assigned_share - others_mean, 4)
    )
    if difference is None:
        verdict = "not computable: the population is empty"
    elif abs(difference) <= TOLERANCE:
        verdict = NEAR
    elif difference > TOLERANCE:
        verdict = ABOVE
    else:
        verdict = BELOW

    assignments = {}
    for b, row in sorted(per_assignment.items(), key=lambda kv: -kv[1]["rules"]):
        own = share(row["own_" + ce.STATE_OPEN], row["rules"])
        other = share(row["others_open"], row["others_assessable"])
        assignments[b] = {
            "rules_assigned_to_this_cell": row["rules"],
            "open_in_this_cell": row["own_" + ce.STATE_OPEN],
            "share_open_in_this_cell": own,
            "open_in_the_other_twelve": row["others_open"],
            "rule_cell_pairs_assessable_in_the_other_twelve": row["others_assessable"],
            "mean_share_open_in_the_other_twelve": other,
            "difference": None if own is None or other is None else round(own - other, 4),
        }

    payload: dict[str, Any] = {
        "question": (
            "over exactly the assessable rules of context_evidence.json, how often is the element "
            "open in the cell the compiler assigned, against how often it is open in the other "
            "reader biosamples?"
        ),
        "status": (
            "descriptive and NOT registered: context_evidence_registration.json fixed the states, "
            "the openness call, the mapping table and the denominator, and this changes none of "
            "them. It is the same reading taken in twelve more biosamples"
        ),
        "why": (
            "the elements are ENCODE registry cCREs, built from DNase over many biosamples, so any "
            "biosample opens some share of them whatever cell a rule names. The assigned cell's "
            "share is uninterpretable until the base rate is beside it"
        ),
        "chromosomes": chroms,
        "population": population,
        "population_is": (
            "every compiled rule whose cell maps to a reader biosample by ontology term and whose "
            "locus lies inside that biosample's measured span on its chromosome: the "
            "open_in_reader plus not_open_in_reader rules of context_evidence.json"
        ),
        "assigned_cell": {
            "open": assigned[ce.STATE_OPEN],
            "not_open": assigned[ce.STATE_NOT_OPEN],
            "share_open": assigned_share,
            "population": f"the {population} assessable rules, each read in its own assigned cell",
        },
        "other_twelve_biosamples": {
            "open": others_open,
            "rule_cell_pairs_assessable": others_assessable,
            "rule_cell_pairs_outside_a_measured_span": others_outside,
            "mean_share_open": others_mean,
            "population": (
                f"the {others_assessable} rule-and-biosample pairs made by reading each of the "
                f"{population} assessable rules in the twelve reader biosamples that are not its "
                "own; a pair whose locus lies outside that biosample's measured span is excluded "
                "and counted separately above"
            ),
        },
        "difference": difference,
        "tolerance": TOLERANCE,
        "tolerance_was_chosen": (
            "before any share was computed, so the reading could not be softened once the numbers "
            "were in. It is not registered, because this comparison is not"
        ),
        "reading": verdict,
        "per_biosample_over_the_same_rules": {
            b: {
                "open": per_biosample[b][ce.STATE_OPEN],
                "not_open": per_biosample[b][ce.STATE_NOT_OPEN],
                "outside_its_measured_span": per_biosample[b][ce.REASON_OUTSIDE_SPAN],
                "share_open_of_the_whole_population": share(per_biosample[b][ce.STATE_OPEN], population),
                "share_open_of_the_rules_it_could_read": share(
                    per_biosample[b][ce.STATE_OPEN],
                    per_biosample[b][ce.STATE_OPEN] + per_biosample[b][ce.STATE_NOT_OPEN],
                ),
                "population": (
                    f"the same {population} assessable rules, every one of them read in this "
                    "biosample whether or not any rule names it"
                ),
            }
            for b in biosamples
        },
        "per_assigned_cell": assignments,
        "no_interval_is_reported": (
            "these are exhaustive counts over a fixed rule set, not an estimate from a sample, so no "
            "confidence interval is computed and none would mean anything here"
        ),
        "not_validation": ce.NOT_VALIDATION,
        "what_this_tests": (
            "whether the cell the compiler assigns is informative about where the element is open. "
            "It does not test whether the rule is right, and open_in_reader is still a consistency "
            "check between two readings of the same ENCODE chromatin"
        ),
        "code_cleanliness": ce.code_cleanliness(ENTRY, OWN_CODE),
        "alphagenome_requests": 0,
        "money": "none: every input was already on disk",
        "seconds": round(time.time() - started, 1),
    }

    inputs = [mf.input_entry(ce.TRACK_METADATA, partition=None)]
    peaks = sorted(
        p.as_posix()
        for cell in ce.READER_TERMS
        for chrom in chroms
        if (p := RESULTS_DIR / reader.peaks_path(cell, chrom).name).exists()
    )
    if len(peaks) <= 200:
        inputs += [mf.input_entry(p, partition=None) for p in peaks]
    else:
        inputs.append(mf.files_entry("reader_v1_dnase_peak_sets", peaks, partition=None))
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "ENCODE DNase-seq narrowPeak, GRCh38, released, 13 biosamples",
                "version": "as reader v1 cached them under data/results/dnase_*_chr*.bed.gz",
            },
            {"accession": "AlphaGenome output track metadata", "version": str(ce.TRACK_METADATA)},
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "openness_call": ce.OPENNESS_CALL,
            "measured_span": ce.SPAN_CALL,
            "tolerance": TOLERANCE,
            "chromosomes": chroms,
        },
        "exclusions": [
            "a rule not assessable in its own assigned cell is not in the population",
            "a rule-and-biosample pair outside that biosample's measured span is excluded from that "
            "biosample's share and counted separately",
        ],
        "partitions": {
            "assigned_cell": "each rule read in the cell the compiler assigned it",
            "other_twelve": "the same rules read in the twelve reader biosamples that are not its own",
            "per_biosample": "the same rules read in each of the thirteen, one biosample at a time",
        },
        "code_cleanliness": ce.code_cleanliness(ENTRY, OWN_CODE),
    }
    if args.supersedes:
        payload["result_manifest"]["supersedes"] = supersedes()
    path = save_result(args.result, payload)
    print(f"{args.result}: {path}")
    print(f"  population {population:,} assessable rules")
    print(f"  assigned cell {assigned_share} open")
    print(f"  other twelve  {others_mean} open over {others_assessable:,} rule-cell pairs")
    print(f"  difference {difference} against a tolerance of {TOLERANCE}")
    print(f"  reading: {verdict}")


if __name__ == "__main__":
    main()
