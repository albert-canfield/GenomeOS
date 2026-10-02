#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Amendment 1 to the gene-identity registration: the compiled target is an identifier, not a symbol.

Written after the chr21 run and before any genome-wide count, and it says exactly what had been read
when it was written. It adds one cause and one lookup; it moves no threshold, no denominator and no
precedence.

    uv run --frozen python scripts/gene_identity_amend_1.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import gene_identity as gi  # noqa: E402
from genomeos.results import load_result, save_result  # noqa: E402

RESULT = "gene_identity_registration_amendment_1"
AMENDS = "gene_identity_registration"
OWN_CODE = (
    "genomeos/attribution/gene_identity.py",
    "scripts/gene_identity_register.py",
    "scripts/gene_identity_amend_1.py",
    "scripts/gene_identity.py",
    "tests/test_gene_identity.py",
)

ADDED = "the_token_matches_two_or_more_symbols_in_the_annotation"


def main() -> None:
    base = load_result(AMENDS)
    if base is None:
        raise SystemExit(f"{AMENDS} is not on disk: an amendment amends something")
    payload = {
        "result": RESULT,
        "date": "2026-10-02",
        "lane": "lane-identity",
        "amends": AMENDS,
        "why": (
            "the registration said the compiled side resolves the rule's target token as a gene_name "
            "or a versionless gene_id. It is neither when the gene's symbol carries a character "
            "BioLang identifiers do not: `compile.ident` replaces every non-word character with an "
            "underscore, so GENCODE's `KRTAP10-1` is written `KRTAP10_1` in the compiled text and an "
            "exact match on the symbol finds nothing. The registration's resolution would have "
            "counted those rules as a token absent from the annotation, which is a statement about "
            "the annotation and would have been the wrong one"
        ),
        "counts_already_read_when_this_was_written": [
            "the chr21 run, data/results/gene_identity_chr21.json at sha 87d5f9c: 5,176 rules, "
            "0 differ, 5,022 agree, 154 unresolvable, of which 152 were "
            "`target_absent_from_the_chromosome_annotation` and 2 the experimental layer",
            "the 152 are tokens such as ERVH48_1, KRTAP12_1 and GET1_SH3BGR, which is what pointed at "
            "compile.ident",
            "the distance test over chr21's 5,176 rules selects 115 rules and all 115 came out agree",
            "no genome-wide count of any kind had been read when this was written",
        ],
        "what_changes": {
            "lookup": (
                "a gene is indexed under its GENCODE gene_name AND under compile.ident of that name, "
                "so a raw symbol still finds its own gene - pilot_bio.gene_tss's choice is unchanged "
                "and the chr21 test that holds the two together still passes - and a compiled token "
                "finds the gene the compiler wrote it from"
            ),
            "the_names_looked_for_in_the_recorded_window_list": (
                "the cache records AlphaGenome's own gene names, which carry the hyphen the compiled "
                "token has lost, so the names looked for are the gene_names the token matches plus "
                "the token itself"
            ),
            "one_cause_added": {ADDED: gi.CAUSES[ADDED]},
            "why_it_is_a_cause_and_not_a_resolution": (
                "compile.ident is not injective, so a token can be matched by two distinct GENCODE "
                "gene_names. The compiled text then does not say which symbol it meant and no id is "
                "recoverable from it. It is refused rather than picked between, and counted on its own"
            ),
        },
        "causes_added": {ADDED: gi.CAUSES[ADDED]},
        "definitions_that_did_not_move": [
            "what a mis-resolution is: unchanged, word for word",
            "the half-window of 500,000 and the second half-window of 524,288: unchanged",
            "the denominator of 440,589 rules and chr21's 5,176: unchanged",
            "the precedence of the causes: unchanged, with the new cause taking the place of "
            "`target_absent_from_the_chromosome_annotation` only where that cause would have been "
            "wrong about why nothing resolved",
            "the gene-body overlap that decides a window candidate: unchanged",
            "the correlation axes and their denominator: unchanged",
        ],
        "also_added_after_the_chr21_run": {
            "what": (
                "within the rules the distance test selects, how many loci of that name the "
                "chromosome carries, and whether the compiled gene's own body reaches into the "
                "scorer window"
            ),
            "why": (
                "the chr21 run selected 115 rules by distance and every one agreed, so why the "
                "distance is large is worth counting rather than assuming. NCAM2 on chr21 is one "
                "gene spanning 544 kb: an element 908 kb from its TSS is still inside the window its "
                "body occupies"
            ),
            "what_it_is_not": (
                "it is a descriptive breakdown of a population the registration already fixed. It "
                "moved no definition, no threshold and no precedence, and it is marked in the result "
                "as added after the chr21 run so no reader can take it for a pre-registered figure"
            ),
        },
        "the_registration_is_not_edited": (
            f"data/results/{AMENDS}.json is left exactly as it was committed at 87d5f9c. This file "
            "sits beside it and the result cites both"
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
            mf.input_entry(Path("data/results/gene_identity_chr21.json")),
        ],
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "cause_added": ADDED,
            "ident": (
                "genomeos.attribution.compile.ident: every non-word character becomes an underscore, "
                "and a leading digit takes a `g_` prefix"
            ),
            "ident_collisions_in_gencode_v50": (
                "counted over all 24 chromosomes before this was written: 0 of 83,343 tokens are "
                "matched by more than one gene_name, so the cause is named and comes out empty "
                "rather than being left unnamed"
            ),
            "half_window": gi.HALF_WINDOW,
            "second_half_window": gi.SENSITIVITY_HALF_WINDOW,
        },
        "exclusions": [
            "nothing is excluded here: this file amends definitions and counts nothing",
            "no threshold, denominator or precedence of the registration is changed",
        ],
        "partitions": {
            "none": "this amendment partitions nothing",
        },
        "code_cleanliness": mf.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
    }
    path = save_result(RESULT, payload)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
