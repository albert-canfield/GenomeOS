# SPDX-License-Identifier: AGPL-3.0-or-later
"""Count the context-evidence state of every compiled rule, per state and per cell.

    uv run --frozen python scripts/context_evidence_census.py
    uv run --frozen python scripts/context_evidence_census.py --chroms chr21 --result context_evidence_chr21

The denominator is **every compiled rule**: each chromosome is compiled and the rule lines are read
back from the compiled text, so the census counts the rules the compiler emits and cannot drift from
them. No model is called and nothing is downloaded; the reading was registered first, in
`data/results/context_evidence_registration.json`.

    uv run python scripts/context_evidence_census.py
    uv run python scripts/context_evidence_census.py --result context_evidence_v2 --supersedes

`--result` changes the name written and nothing else; `--supersedes` records which committed result the
run is written beside. The committed data/results/context_evidence.json is kept byte for byte: its sha256
is a declared input of four other committed results.

What the states mean, what they do not mean, and why `open_in_reader` is never validation is in
`genomeos.attribution.context_evidence` and is copied into this result, so the result cannot be read
without them.
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
from genomeos.attribution import compile as cp  # noqa: E402
from genomeos.attribution import context_evidence as ce  # noqa: E402
from genomeos.genome import reader  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

RESULT = "context_evidence"
CHROMS = tuple(f"chr{c}" for c in [*range(1, 23), "X", "Y"])
RULES_ON_THE_RECORD = 440_589
#: Above this many input files the manifest records them as groups rather than one by one.
INPUTS_PER_FILE_MAX = 200
#: The input a census of many files used to declare to say that its groups could not be checked.
INPUTS_AS_GROUPS = "inputs_recorded_as_groups"

#: This lane's own files. Everything else uncommitted in this shared checkout belongs to another lane.
#: This script, as the entry whose transitive import closure is the counting path of its results
#: (genomeos.manifest.counting_path, the one implementation). Named rather than derived from
#: __file__ so the closure is the same however the script is invoked.
ENTRY = "scripts/context_evidence_census.py"

OWN_CODE = (
    "genomeos/attribution/context_evidence.py",
    "genomeos/attribution/compile.py",
    "genomeos/lang/grammar.py",
    "genomeos/lang/parser.py",
    "genomeos/ir/model.py",
    "scripts/context_evidence_register.py",
    "scripts/context_evidence_census.py",
    "tests/test_context_evidence.py",
    "docs/BIOLANG-GRAMMAR.md",
)
ROOT = Path(__file__).resolve().parents[1]

#: The result files the compiler reads per chromosome, as globs, so the manifest names its inputs
#: without a hand-written list of several hundred paths going stale.
INPUT_GLOBS = (
    "budget_chr*.json",
    "budget_axes_chr*.json",
    "variation_chr*.json",
    "duplication_chr*.json",
    "domains_chr*.json",
    "unknown_chr*.json",
    "reader_*_chr*.json",
    "ccres_chr*.bed.gz",
    "enhancer_targets_chr*.json",
    "enhancer_targets_all_chr*.json",
    "constrained_targets_chr*.json",
)


def rules_of(text: str) -> list[tuple[str, str]]:
    """(cell label, context-evidence value) for every rule line of a compiled program."""
    out = []
    for line in text.splitlines():
        if not line.startswith("rule "):
            continue
        if "when: cell_type = " not in line:
            raise ValueError(f"a rule with no cell: {line[:120]}")
        cell = line.split("when: cell_type = ", 1)[1].split(";", 1)[0].strip()
        if "context_evidence: " not in line:
            raise ValueError(f"a rule with no context evidence: {line[:120]}")
        value = line.split("context_evidence: ", 1)[1].split(";", 1)[0].strip()
        out.append((cell, value))
    return out


#: Why this census is re-recorded under a name of its own rather than over the committed file. The
#: committed manifest declares 13 inputs of which 12 are group labels with no member list, so a rebuild
#: in a second environment has no path to open or hash: `files_entry` did not name the members of a group
#: until 2026-10-02. The committed bytes cannot be replaced to fix it: the sha256 of
#: data/results/context_evidence.json is a declared input of four other committed results
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
#: What the committed file would have declared had it been written by the code of 6a8c794 or later: a
#: fourteenth input `inputs_recorded_as_groups`, sha256 "n/a", whose note said that
#: scripts/manifest_rebuild.py "reports every group absent". That stopped being true at ebbded6, which
#: opens and hashes a group member by member, and an input with no bytes behind it is now a named reason
#: to refuse a verdict (manifest_rebuild.NOT_HASHED), so the entry is not written at all.
THE_FOURTEENTH_INPUT_THAT_IS_NOT_WRITTEN = (
    "the writer of 6a8c794 appended a fourteenth input `inputs_recorded_as_groups` with sha256 'n/a' "
    "and a note saying scripts/manifest_rebuild.py 'reports every group absent'; the committed file at "
    "data/results/context_evidence.json predates that writer and declares 13 inputs without it. The "
    "note was untrue from ebbded6 onward, which opens and hashes a group one member at a time, and an "
    "input whose sha256 is 'n/a' is itself a reason to refuse a verdict, so no such entry is written "
    "here and every group names its members instead"
)


def supersedes() -> dict[str, Any]:
    """The committed result this run is written beside, never over, with the bytes it is kept at."""
    p = ROOT / RESULTS_DIR / f"{RESULT}.json"
    entry = mf.input_entry(p, partition=None)
    old = json.loads(p.read_text())
    return {
        "file": (RESULTS_DIR / f"{RESULT}.json").as_posix(),
        "date": old.get("date"),
        "sha256": entry["sha256"],
        "bytes": entry["bytes"],
        "inputs_it_declared": len(((old.get(mf.KEY) or {}).get("inputs")) or []),
        "kept": "unchanged; this run is written beside it under a new name, not over it",
        "why": WHY_A_NEW_NAME,
        "declared_as_an_input_by": list(DECLARED_AS_AN_INPUT_BY),
        "inputs_recorded_as_groups": THE_FOURTEENTH_INPUT_THAT_IS_NOT_WRITTEN,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", default=",".join(CHROMS))
    ap.add_argument(
        "--result",
        default=RESULT,
        help=(
            "the result name to write. A run over fewer than all 24 chromosomes writes under a name "
            "of its own, so a partial census can never overwrite the genome-wide one"
        ),
    )
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
    name = args.result
    if name == RESULT and tuple(chroms) != CHROMS:
        raise SystemExit(
            f"{RESULT} is the genome-wide census; a run over {len(chroms)} chromosomes needs --result"
        )

    started = time.time()
    per_state: Counter[str] = Counter()
    per_reason: Counter[str] = Counter()
    per_biosample: dict[str, Counter[str]] = defaultdict(Counter)
    per_cell: dict[str, Counter[str]] = defaultdict(Counter)
    per_chrom: dict[str, dict[str, int]] = {}
    total = 0
    for chrom in chroms:
        text = cp.compile_chromosome(chrom)
        rows = rules_of(text)
        here: Counter[str] = Counter()
        for cell, raw in rows:
            state, rest = ce.parse_value(raw)
            per_state[state] += 1
            here[state] += 1
            per_cell[cell][state] += 1
            if state == ce.STATE_NOT_ASSESSABLE:
                per_reason[rest[0]] += 1
                if len(rest) > 1:
                    per_biosample[rest[1]][rest[0]] += 1
            else:
                per_biosample[rest[0]][state] += 1
        per_chrom[chrom] = {"rules": len(rows), **{s: here[s] for s in ce.STATES}}
        total += len(rows)
        print(f"  {chrom}: {len(rows):,} rules, {here[ce.STATE_OPEN]:,} open_in_reader", flush=True)

    table = ce.mapping()
    assessed = per_state[ce.STATE_OPEN] + per_state[ce.STATE_NOT_OPEN]
    # every cell label that carries rules, with the term it resolved to and why it could not be read
    unmapped: dict[str, list[dict[str, Any]]] = {
        "no_registered_ontology_term": [],
        "term_is_not_one_of_the_thirteen_biosamples": [],
        "label_resolves_to_more_than_one_term": [],
    }
    cells: dict[str, dict[str, Any]] = {}
    for cell, counts in sorted(per_cell.items(), key=lambda kv: (-sum(kv[1].values()), kv[0])):
        biosample, term, why = ce.reader_for(cell, table)
        rules_here = sum(counts.values())
        cells[cell] = {
            "rules": rules_here,
            "ontology_term": term,
            "reader_biosample": biosample,
            **{s: counts[s] for s in ce.STATES},
        }
        if biosample is not None:
            continue
        cells[cell]["not_assessable_because"] = why
        row = {"cell": cell, "rules": rules_here, "ontology_term": term, "why": why}
        if cell in table.ambiguous:
            unmapped["label_resolves_to_more_than_one_term"].append(row)
        elif term is None:
            unmapped["no_registered_ontology_term"].append(row)
        else:
            unmapped["term_is_not_one_of_the_thirteen_biosamples"].append(row)

    payload: dict[str, Any] = {
        "question": (
            "for every compiled rule, is the element open in the cell the rule's `when:` names, at "
            "reader v1's own call, or is there a reason no reading can be taken?"
        ),
        "registration": "data/results/context_evidence_registration.json",
        "chromosomes": chroms,
        "scope": (
            "every compiled rule of all 24 chromosomes"
            if tuple(chroms) == CHROMS
            else f"the compiled rules of {len(chroms)} of 24 chromosomes: {', '.join(chroms)}"
        ),
        "rules": total,
        "rules_on_the_record": RULES_ON_THE_RECORD,
        "denominator_agrees_with_the_record": (
            total == RULES_ON_THE_RECORD
            if tuple(chroms) == CHROMS
            else "n/a: the record is the figure over all 24 chromosomes and this run covers fewer"
        ),
        "per_state": {s: per_state[s] for s in ce.STATES},
        "not_assessable_by_reason": {r: per_reason[r] for r in ce.NOT_ASSESSABLE_REASONS},
        "assessable_rules": assessed,
        "open_share_of_assessable_rules": (
            None if not assessed else round(per_state[ce.STATE_OPEN] / assessed, 4)
        ),
        "open_share_population": (
            "the rules whose cell maps to a reader biosample and whose locus lies inside that "
            "biosample's measured span on the chromosome: open_in_reader + not_open_in_reader"
        ),
        "no_interval_is_reported": (
            "these are exhaustive counts over every compiled rule, not an estimate from a sample, so "
            "no confidence interval is computed and none would mean anything here"
        ),
        "per_reader_biosample": {
            b: dict(sorted(counts.items())) for b, counts in sorted(per_biosample.items())
        },
        "per_cell": cells,
        "cells_without_a_reading": {
            "counts": {k: len(v) for k, v in unmapped.items()},
            "rules": {k: sum(r["rules"] for r in v) for k, v in unmapped.items()},
            "cells": {k: v for k, v in unmapped.items()},
        },
        "per_chromosome": per_chrom,
        **ce.registration(),
        "reading": {
            "moves": "honest: an asserted cell context becomes an evidenced or a labelled one",
            "secondary": (
                "correct: a not_open_in_reader rule is a located candidate defect, a place where "
                "the program asserts activity in a cell whose own chromatin shows no open element"
            ),
            "does_not_move": (
                "complete: no rule is added, none is removed, and no part of the genome becomes "
                "attributed that was not attributed before"
            ),
            "never": (
                "open_in_reader is never reported as validation, support or confirmation of a rule, "
                "and not_open_in_reader is never reported as the element being closed or the rule "
                "being contradicted"
            ),
        },
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
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
    groups = {"reader_v1_dnase_peak_sets": peaks}
    for glob in INPUT_GLOBS:
        for chrom in chroms:
            found = sorted(p.as_posix() for p in RESULTS_DIR.glob(glob.replace("chr*", chrom)))
            if found:
                groups.setdefault(f"compiler_inputs:{glob}", []).extend(found)
    every = sorted({p for v in groups.values() for p in v})
    # A `files_entry` records a group of files under a label, and `scripts/manifest_rebuild.py` checks
    # each input by its `path`, so a label is reported absent whatever is on disk. Where the run reads
    # few enough files, each one is recorded by its own path instead and the rebuild can check them.
    #
    # Amended 2026-10-02: the paragraph above describes the world before ebbded6. A group now names
    # every member with its own sha256 and byte count (`genomeos.manifest.files_entry`), and the
    # rebuild opens and hashes them one file at a time, so a label is no longer unresolvable and the
    # entry below - which said it was, under sha256 "n/a" - is written only in the case it actually
    # describes: a group that names no member. That case no longer arises here, and the committed
    # data/results/context_evidence.json declares 13 inputs without any such entry, having been written
    # before 6a8c794 added it; what it said is carried into the `supersedes` note of a result written
    # under another name. An input whose sha256 is "n/a" is itself a reason to refuse a verdict
    # (manifest_rebuild.NOT_HASHED), so it is not declared where it would not be true.
    if len(every) <= INPUTS_PER_FILE_MAX:
        inputs += [mf.input_entry(p, partition=None) for p in every]
    else:
        for label, found in groups.items():
            if found:
                inputs.append(mf.files_entry(label, sorted(set(found)), partition=None))
        inputs.append(
            {
                "path": "inputs_recorded_as_groups",
                "sha256": "n/a",
                "partition": None,
                "note": (
                    f"{len(every)} files over {len(chroms)} chromosomes, too many to record one by "
                    "one, so each group above carries one sha256 over its files in sorted path "
                    "order. scripts/manifest_rebuild.py checks an input by its path and a group's "
                    "label is not one, so it reports every group absent; a run over fewer "
                    "chromosomes records each file by its own path"
                ),
            }
        )
        # ... and taken back out again when every group did name its members, which is the only state
        # `files_entry` can now produce. Written this way round, rather than by not appending it, so
        # that the entry above keeps the shape and the indentation 6a8c794 gave it: it is a statement
        # about a group that cannot be checked, and it stands only while one cannot be.
        if all(e.get("members") for e in inputs if e.get("group")):
            inputs = [e for e in inputs if e["path"] != INPUTS_AS_GROUPS]
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "ENCODE DNase-seq narrowPeak, GRCh38, released, 13 biosamples",
                "version": "as reader v1 cached them under data/results/dnase_*_chr*.bed.gz",
            },
            {"accession": "AlphaGenome output track metadata", "version": str(ce.TRACK_METADATA)},
            {
                "accession": "the compiled non-coding programs of the 24 chromosomes",
                "version": "recompiled in this run from the results on disk, not read from a file",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "openness_call": ce.OPENNESS_CALL,
            "measured_span": ce.SPAN_CALL,
            "states": list(ce.STATES),
            "not_assessable_reasons": list(ce.NOT_ASSESSABLE_REASONS),
            "chromosomes": chroms,
        },
        "exclusions": [
            "no rule is excluded: the denominator is every rule the compiler emits",
            "no new cut-off is introduced; reader v1's registered call is used unchanged",
        ],
        "partitions": {
            "per_state": "the three states, with the two reasons under not_assessable",
            "per_cell": "the rule's own `when: cell_type` label, every label counted",
            "per_biosample": "the reader biosample a label maps to by ontology term",
        },
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }
    if args.supersedes:
        payload["result_manifest"]["supersedes"] = supersedes()
    path = save_result(name, payload)
    print(f"{name}: {path}")
    print(f"  rules {total:,} (on the record {RULES_ON_THE_RECORD:,})")
    for s in ce.STATES:
        print(f"  {s}: {per_state[s]:,}")
    for r in ce.NOT_ASSESSABLE_REASONS:
        print(f"    {r}: {per_reason[r]:,}")


if __name__ == "__main__":
    main()
