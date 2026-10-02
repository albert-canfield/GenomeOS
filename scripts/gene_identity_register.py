#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register the gene-identity check before any count of it is read.

Writes `data/results/gene_identity_registration.json`: what counts as a mis-resolution, what is done
with a symbol whose id cannot be recovered, the precedence the two are decided in, the population,
the annotation by sha256, and the reconciliation and correlation this lane will report. Nothing here
is a count over the compiled programs, and no threshold in it may move afterwards.

    uv run --frozen python scripts/gene_identity_register.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import context_evidence as ce  # noqa: E402
from genomeos.attribution import gene_identity as gi  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "gene_identity_registration"
OWN_CODE = (
    "genomeos/attribution/gene_identity.py",
    "scripts/gene_identity_register.py",
    "scripts/gene_identity.py",
    "tests/test_gene_identity.py",
)


def main() -> None:
    annotation = [mf.input_entry(gi.annotation_path(c)) for c in gi.CHROMS if gi.annotation_path(c).exists()]
    programs = [mf.input_entry(gi.program_path(c)) for c in gi.CHROMS if gi.program_path(c).exists()]
    window = [
        mf.input_entry(gi.WINDOW_CACHE / f"{c}.json.gz")
        for c in gi.CHROMS
        if (gi.WINDOW_CACHE / f"{c}.json.gz").exists()
    ]
    payload = {
        "result": RESULT,
        "date": "2026-10-02",
        "lane": "lane-identity",
        "status": (
            "registered before any count over the compiled programs was read. No outcome of the "
            "comparison this fixes has been computed at the time of writing"
        ),
        "question": (
            "a compiled rule names its target gene by symbol and the AlphaGenome deletion answer it "
            "is built from measured one gene. For every compiled rule: does the target token resolve "
            "to the same Ensembl gene id as the gene the deletion answer measured?"
        ),
        "why": (
            "lane-notopen registered that 1,678 rules carry a TSS further from the element than half "
            "the 1 Mb deletion window, so a repeated gene symbol resolved to another locus - a fault "
            "in the distance, not the rule. That is a distance test over a subset. An id comparison "
            "asks the identity question directly and over the whole compiled set"
        ),
        "population": {
            "rules": (
                "every `rule` line of the 24 compiled programs, data/knowledge/compiled/noncoding_chr*.bio"
            ),
            "why_the_programs_and_not_a_re_enumeration": (
                "the compiled text is what the model asserts, and its rule lines are countable "
                "without any git-ignored intermediate. The count must equal 440,589, the denominator "
                "context_evidence_registration.json fixed, and the run refuses to write anything if "
                "it does not"
            ),
            "expected_total": 440589,
            "per_chromosome_check": "chr21 must hold 5,176 rules, the figure the chr21 census fixed",
        },
        "mis_resolution": gi.MIS_RESOLUTION,
        "compiled_side": gi.COMPILED_RESOLUTION,
        "measured_side": gi.MEASURED_RESOLUTION,
        "what_the_deletion_answer_records": (
            "a gene NAME and no Ensembl gene id: `predicted_coding` in "
            "data/knowledge/alphagenome/all_elements/<chrom>.json is "
            "{gene, action, log2_fold_change, tissue, strength, confidence, basis}, where `gene` is a "
            "symbol or, for a gene with no symbol, a bare versionless Ensembl id. The measured side "
            "is therefore recovered from the record and never read off it"
        ),
        "unresolvable_symbols": {
            "rule": (
                "a rule whose two sides cannot both be resolved is counted under a named cause and "
                "is never resolved on a preference, never dropped, and never swept into a residue. "
                "agree + differ + the causes sum to the population"
            ),
            "causes": gi.CAUSES,
            "precedence": [
                "target_absent_from_the_chromosome_annotation",
                "rule_of_the_experimental_layer",
                "element_absent_from_the_window_cache",
                "target_not_in_the_recorded_window_list",
                "target_named_more_than_once_in_the_recorded_window_list",
                "no_annotated_locus_of_that_name_overlaps_the_scorer_window",
                "two_or_more_annotated_loci_of_that_name_in_the_scorer_window",
                "then agree or differ",
            ],
            "precedence_is_fixed": (
                "this order is fixed here and may not be re-ordered after a count is read. A rule "
                "that would fall under two causes is counted under the first of them"
            ),
        },
        "the_1678_reconciliation": {
            "will_be_reported": (
                "the same distance test lane-notopen ran - the compiled side's element-midpoint-to-TSS "
                "distance above not_open_profile.SCORER_HALF_WINDOW - computed over the whole compiled "
                "set, and cross-tabulated against the outcome classes of this check"
            ),
            "populations_are_named_and_not_equated": (
                "lane-notopen's 1,678 is that distance test restricted to the 55,084 rules in state "
                "not_open_in_reader. This lane's population is all 440,589 rules and it does not "
                "re-derive the openness state. The two numbers are stated with their populations and "
                "the 1,678 is never asserted to be the same set as any count here"
            ),
            "the_containment_is_computed_not_assumed": (
                "a compiled TSS outside the window does not by itself put the compiled and the "
                "measured gene apart: a long gene's body can overlap the window while its TSS does "
                "not. The overlap between the distance test and `differ` is therefore counted, not "
                "argued"
            ),
        },
        "correlations": {
            "axes": [
                "the activity axis, the element block's own `activity:` value as the compiler wrote it",
                "the effect band, " + gi.EFFECT_BAND_CALL,
                "the distance band, executor._band over the compiled side's distance with the "
                "residual row above executor.WIDE_MAX_DISTANCE",
            ],
            "denominator": (
                "within each level, the resolvable rules of that level: agree + differ. The "
                "unresolvable rules of a level are reported beside the share and never inside it"
            ),
            "reported_descriptively": (
                "a share that differs between levels is reported as it came out. No cause is "
                "inferred, no axis is said to explain another, and no statement is made that the "
                "model's predictions are wrong by any of these proportions"
            ),
        },
        "sensitivity_fixed_in_advance": {
            "half_window": gi.HALF_WINDOW,
            "second_half_window": gi.SENSITIVITY_HALF_WINDOW,
            "why": (
                "the project fixed half the scorer window at 500,000; AlphaGenome's 1 Mb interval is "
                "2**20 bp, half of which is 524,288. Every count is reported at both, fixed here "
                "before either was read, so no figure can be read as chosen for the rounder number"
            ),
        },
        "the_experimental_layer": (
            "the rules whose evidence is a CRISPRi screen have no deletion answer and are counted "
            "under `rule_of_the_experimental_layer`. The screen does record its own Ensembl gene id, "
            "so that comparison is reported as a separate figure with its own denominator and is "
            "never added into the deletion-answer counts"
        ),
        "what_this_cannot_establish": gi.limitations(),
        "imported_readings": gi.imported_readings(),
        "alphagenome_requests": 0,
        "network_requests": 0,
        "money": "none: every input was already on disk",
        "code": {
            "module": "genomeos/attribution/gene_identity.py",
            "register": "scripts/gene_identity_register.py",
            "count": "scripts/gene_identity.py",
            "tests": "tests/test_gene_identity.py",
        },
        "code_cleanliness": ce.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
        "inputs_named_by_sha256": {
            "annotation": annotation,
            "compiled_programs": programs,
            "window_cache": window,
            "all_three_are_git_ignored": (
                "data/reference/, data/knowledge/ and data/cache/ are git-ignored, so a manifest "
                "rebuild of any result of this lane stops at the inputs rather than rebuilding them"
            ),
        },
    }
    path = save_result(RESULT, payload)
    print(f"wrote {path}")
    print(f"annotation files {len(annotation)}, programs {len(programs)}, window caches {len(window)}")


if __name__ == "__main__":
    main()
