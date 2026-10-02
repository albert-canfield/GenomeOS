# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register the context-evidence reading before any genome-wide count is read.

    uv run --frozen python scripts/context_evidence_register.py

What is fixed here, and cannot move afterwards: the three state names, the two reasons a reading may
be `not_assessable`, the openness call (reader v1's own, unchanged), the definition of a biosample's
measured span, the ontology mapping table with everything it cannot map, the denominator (every
compiled rule) and the wording the result is allowed to be read in.

Everything it writes comes from `genomeos.attribution.context_evidence`, so the registration and the
code that applies it cannot drift apart. It reads the mapping table and the peak-set inventory, and
no rule: no state is computed here and no count over the compiled programs is taken.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import context_evidence as ce  # noqa: E402
from genomeos.genome import reader  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "context_evidence_registration"
CHROMS = tuple(f"chr{c}" for c in [*range(1, 23), "X", "Y"])
RULES_ON_THE_RECORD = 440_589


def peak_inventory(results_dir: Path) -> dict[str, Any]:
    """Which of the thirteen biosamples has a cached peak set for which chromosome, and how many peaks.

    Counted before any state is read, so the denominator of every `region_outside_measured_span` is
    on the record rather than inferred from the outcome.
    """
    per_cell: dict[str, dict[str, Any]] = {}
    for cell in ce.READER_TERMS:
        files, peaks, missing = 0, 0, []
        for chrom in CHROMS:
            path = results_dir / reader.peaks_path(cell, chrom).name
            if not path.exists():
                missing.append(chrom)
                continue
            files += 1
            peaks += sum(1 for _ in reader.load_peaks(cell, chrom))
        per_cell[cell] = {
            "chromosomes_with_a_cached_peak_set": files,
            "chromosomes_without_one": missing,
            "peaks": peaks,
        }
    return per_cell


def main() -> None:
    table = ce.mapping()
    results_dir = reader.RESULTS
    inventory = peak_inventory(results_dir)

    payload: dict[str, Any] = {
        "status": (
            "registered before any genome-wide count was read; no state over the compiled programs "
            "is counted here"
        ),
        "lane": "lane-context2",
        "question": (
            "for every compiled rule, is the element open in the cell the rule's `when:` names, at "
            "reader v1's own call, or is there a reason no reading can be taken?"
        ),
        "why": (
            "R1 (2026-09-28) made every rule carry its cell as an executable `when:`. That cell is "
            "asserted: it is the winning AlphaGenome track's tissue, or the cell a CRISPRi screen "
            "silenced in, and nothing checked whether the element is open there. This gives each "
            "rule's context chromatin evidence, or a label saying why it has none"
        ),
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
        "denominator_on_the_record": {
            "rules": RULES_ON_THE_RECORD,
            "source": (
                "data/results/when_census.json, generated_bio_untracked: 24 programs, 440,589 "
                "rule/plain `when` clauses"
            ),
            "the_census_counts": (
                "the rules the compiler emits at the sha of its own run, and reports that total "
                "beside this one rather than assuming they agree"
            ),
        },
        "peak_sets_cached_before_the_run": inventory,
        "mapping_table_as_read": {
            "labels": len(table.terms),
            "ambiguous_labels": len(table.ambiguous),
            "per_ontology": dict(sorted(table.prefixes.items())),
        },
        "built_while_registering": (
            "the compiler was run once on chr21 during development, to check that the field is "
            "written on every rule line and parses back; 5,176 rules. No genome-wide count has been "
            "read at the time of writing, and no threshold here was chosen after seeing an outcome"
        ),
        "code": {
            "module": "genomeos/attribution/context_evidence.py",
            "compiler": "genomeos/attribution/compile.py",
            "grammar": "genomeos/lang/grammar.py, genomeos/lang/parser.py, genomeos/ir/model.py",
            "tests": "tests/test_context_evidence.py",
            "register": "scripts/context_evidence_register.py",
            "census": "scripts/context_evidence_census.py",
        },
        "alphagenome_requests": 0,
        "money": "none: every input was already on disk",
    }

    inputs = [mf.input_entry(ce.TRACK_METADATA, partition=None)]
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "ENCODE DNase-seq narrowPeak, GRCh38, released, 13 biosamples",
                "version": "as reader v1 cached them under data/results/dnase_*_chr*.bed.gz",
            },
            {
                "accession": "AlphaGenome output track metadata",
                "version": str(ce.TRACK_METADATA),
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
    }
    path = save_result(RESULT, payload)
    print(f"{RESULT}: {path}")
    print(f"  states {', '.join(ce.STATES)}")
    print(f"  mapping {len(table.terms)} labels, {len(table.ambiguous)} ambiguous")
    print(f"  denominator on the record {RULES_ON_THE_RECORD:,} compiled rules")


if __name__ == "__main__":
    main()
