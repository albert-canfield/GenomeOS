#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Amendment 2: the denominator a count of zero mis-resolutions has to be read against.

Written after the genome-wide count and naming it. It adds one descriptive breakdown and changes no
definition, no threshold, no denominator and no precedence. It exists because the headline figure came
out at zero, and a zero without the number of rules the comparison could have separated at all would
read as a far stronger statement than the evidence carries.

    uv run --frozen python scripts/gene_identity_amend_2.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import gene_identity as gi  # noqa: E402
from genomeos.results import load_result, save_result  # noqa: E402

RESULT = "gene_identity_registration_amendment_2"
AMENDS = "gene_identity_registration"
OWN_CODE = (
    "genomeos/attribution/gene_identity.py",
    "scripts/gene_identity_register.py",
    "scripts/gene_identity_amend_1.py",
    "scripts/gene_identity_amend_2.py",
    "scripts/gene_identity.py",
    "tests/test_gene_identity.py",
)


def main() -> None:
    if load_result(AMENDS) is None:
        raise SystemExit(f"{AMENDS} is not on disk: an amendment amends something")
    result = load_result("gene_identity")
    if result is None:
        raise SystemExit("gene_identity is not on disk: this amendment names the count it follows")
    counts = result["counts"]
    caught = result["what_the_comparison_could_have_caught"]
    payload = {
        "result": RESULT,
        "date": "2026-10-02",
        "lane": "lane-identity",
        "amends": AMENDS,
        "after": ["gene_identity_registration_amendment_1"],
        "why": (
            "the genome-wide count came out at 0 rules whose target resolves to a different Ensembl "
            "gene id from the gene the deletion answer measured. A zero is the one outcome that can "
            "be misread as stronger than it is, because the two sides resolve to the same locus by "
            "construction wherever the target name has only one locus on its chromosome. Without the "
            "count of rules the comparison could have separated at all, the figure would look like a "
            "test the compiled set passed rather than a test that mostly could not fire"
        ),
        "counts_already_read_when_this_was_written": [
            f"the genome-wide run: {counts['differ']} differ, {counts['agree']} agree, "
            f"{counts['unresolvable']} unresolvable over {result['population']['rules']} rules",
            "the unresolvable causes: "
            + ", ".join(f"{k} {v}" for k, v in counts["unresolvable_by_cause"].items()),
            "the distance test over the whole population: "
            f"{result['the_1678']['the_same_distance_test_over_this_whole_population']} rules",
            "the chr21 run, where all 5,176 rules name a target with exactly one locus on the "
            "chromosome, which is what pointed at this denominator",
        ],
        "what_is_added": {
            "breakdown": (
                "every rule is counted by how many loci of its target name its chromosome carries - "
                "0, 1, 2, or 3 and above - crossed with its outcome, and the number of rules a "
                "differing id was possible for at all is reported from it"
            ),
            "what_it_establishes": (
                "where the name has one locus, both sides resolve to that locus whatever the window "
                "says, so `agree` there is a property of the annotation and not a check that could "
                "have failed. Only the rows with two or more loci are rules the comparison could "
                "have separated"
            ),
            "the_figures_it_came_out_at": {
                "the_rules_a_differing_id_was_possible_for": caught[
                    "the_rules_a_differing_id_was_possible_for"
                ],
                "of_those_the_resolvable_ones": caught["of_those_the_resolvable_ones"],
                "by_loci_of_that_name_on_the_chromosome": caught["by_loci_of_that_name_on_the_chromosome"],
            },
        },
        "definitions_that_did_not_move": [
            "what a mis-resolution is: unchanged, word for word",
            "the half-window of 500,000 and the second half-window of 524,288: unchanged",
            "the denominator of 440,589 rules, and chr21's 5,176 as the per-chromosome check: both unchanged",
            "the seven causes of the registration and the one amendment 1 added: unchanged",
            "the precedence between them: unchanged",
            "the correlation axes and their denominator: unchanged",
            "no rule moved from one outcome to another because of this breakdown: it is read off the "
            "same classification and the outcome counts are identical",
        ],
        "the_registration_is_not_edited": (
            f"data/results/{AMENDS}.json is left exactly as it was committed at 87d5f9c, and "
            "amendment 1 as it was committed at d108236. This file sits beside them"
        ),
        "alphagenome_requests": 0,
        "network_requests": 0,
        "money": "none: every input was already on disk",
        "code_cleanliness": mf.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
    }
    payload["result_manifest"] = {
        "sources": [
            {"accession": "GENCODE", "version": "v50, data/reference/gencode_v50_chr*.gff3.gz"},
            {
                "accession": "the compiled non-coding programs of the 24 chromosomes",
                "version": "data/knowledge/compiled/noncoding_chr*.bio, read as they are on disk",
            },
        ],
        "inputs": [
            mf.input_entry(Path(f"data/results/{AMENDS}.json")),
            mf.input_entry(Path("data/results/gene_identity_registration_amendment_1.json")),
            mf.input_entry(Path("data/results/gene_identity.json")),
        ],
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "loci_of_a_name": (
                "how many GENCODE v50 genes of the rule's target name the element's chromosome "
                "carries, counted under the gene_name and the compile.ident index amendment 1 fixed"
            ),
            "half_window": gi.HALF_WINDOW,
            "second_half_window": gi.SENSITIVITY_HALF_WINDOW,
        },
        "exclusions": [
            "nothing is excluded here: this file adds a breakdown and changes no definition",
            "no rule is reclassified: the outcome counts are the same before and after it",
        ],
        "partitions": {
            "by_loci_on_the_chromosome": "how many loci of the target name the chromosome carries"
        },
        "code_cleanliness": mf.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
    }
    path = save_result(RESULT, payload)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
