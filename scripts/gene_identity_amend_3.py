#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Amendment 3: the registration's own reason for one lookup was vacuous, and it is corrected here.

The registration called the direct Ensembl-id lookup "the one extension beyond the project's
resolver", because `pilot_bio.gene_tss` skips genes with no symbol. GENCODE v50 has no gene with an
empty `gene_name`: an unnamed gene carries its own versionless id as that name. So `gene_tss` resolves
those tokens too, there is no extension in effect, and no rule's outcome turns on the difference. The
figures are measured here rather than asserted.

No count changes. This amendment corrects a justification, not a definition or a number.

    uv run --frozen python scripts/gene_identity_amend_3.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import gene_identity as gi  # noqa: E402
from genomeos.results import load_result, save_result  # noqa: E402

RESULT = "gene_identity_registration_amendment_3"
AMENDS = "gene_identity_registration"
OWN_CODE = (
    "genomeos/attribution/gene_identity.py",
    "scripts/gene_identity_register.py",
    "scripts/gene_identity_amend_1.py",
    "scripts/gene_identity_amend_2.py",
    "scripts/gene_identity_amend_3.py",
    "scripts/gene_identity.py",
    "tests/test_gene_identity.py",
)


def measure() -> dict[str, int]:
    """The three figures this amendment rests on, counted over every chromosome's annotation."""
    rows_total = empty_name = name_is_own_id = name_is_another_id = 0
    for chrom in gi.CHROMS:
        rows = gi.annotation(chrom)
        ids = {g.gene_id for g in rows}
        rows_total += len(rows)
        empty_name += sum(1 for g in rows if not g.symbol)
        name_is_own_id += sum(1 for g in rows if g.symbol == g.gene_id)
        name_is_another_id += sum(1 for g in rows if g.symbol in ids and g.symbol != g.gene_id)
    return {
        "gene_rows": rows_total,
        "rows_with_an_empty_gene_name": empty_name,
        "rows_whose_gene_name_is_their_own_versionless_id": name_is_own_id,
        "rows_whose_gene_name_is_another_genes_id": name_is_another_id,
    }


def main() -> None:
    if load_result(AMENDS) is None:
        raise SystemExit(f"{AMENDS} is not on disk: an amendment amends something")
    result = load_result("gene_identity")
    if result is None:
        raise SystemExit("gene_identity is not on disk: this amendment names the count it follows")
    figures = measure()
    payload = {
        "result": RESULT,
        "date": "2026-10-02",
        "lane": "lane-identity",
        "amends": AMENDS,
        "after": [
            "gene_identity_registration_amendment_1",
            "gene_identity_registration_amendment_2",
        ],
        "what_the_registration_said": (
            "that a target token which is itself a versionless Ensembl id resolves to that gene "
            "directly, 'which gene_tss cannot do because it skips genes with no symbol', and that this "
            "is 'the one extension beyond the project's resolver'"
        ),
        "why_that_reason_does_not_hold": (
            "GENCODE v50 gives every gene a gene_name, and an unnamed gene's gene_name is its own "
            "versionless Ensembl id. pilot_bio.gene_tss skips a gene whose symbol is empty, and no "
            "gene's is, so it resolves those tokens as readily as any other name. There is no "
            "extension in effect: the direct lookup and the name lookup return the same gene"
        ),
        "measured_rather_than_asserted": figures,
        "so_the_two_lookups_cannot_disagree": (
            f"{figures['rows_whose_gene_name_is_another_genes_id']} of the "
            f"{figures['gene_rows']} gene rows carry another gene's versionless id as their "
            "gene_name, so a token cannot resolve to one gene by id and a different gene by name"
        ),
        "what_this_changes_in_the_counts": (
            "nothing. The direct lookup is kept because it is the exact one, and it returns the gene "
            f"the name lookup would return. The result's counts - {result['counts']['differ']} differ, "
            f"{result['counts']['agree']} agree, {result['counts']['unresolvable']} unresolvable - are "
            "the same either way, and no rule moves between outcomes"
        ),
        "what_it_does_change": (
            "a claim this lane's own registration made about another module's behaviour, and the "
            "reading that followed from it. The first draft of the ATTRIBUTION section said the 4,340 "
            "rules naming their target by a bare Ensembl id 'carry no distance in any reading built "
            "on' gene_tss. That was wrong and is corrected in the section before it was pushed: those "
            "rules do carry a distance, and any reading built on gene_tss includes them"
        ),
        "the_4340": (
            f"{result['target_token_was_itself_an_ensembl_id']} rules name their target by a bare "
            "versionless Ensembl id rather than a symbol, which is how the deletion answer names a "
            "gene that has none. For those rules the answer does record the id, so both sides are "
            "exact rather than recovered. That remains a true descriptive count"
        ),
        "definitions_that_did_not_move": [
            "what a mis-resolution is: unchanged, word for word",
            "the half-window of 500,000 and the second half-window of 524,288: unchanged",
            "the denominator of 440,589 rules: unchanged",
            "the eight causes and their precedence: unchanged",
            "every count in data/results/gene_identity.json: unchanged",
        ],
        "the_registration_is_not_edited": (
            f"data/results/{AMENDS}.json stands as committed at 87d5f9c, amendment 1 at d108236 and "
            "amendment 2 beside them. This file corrects the reason, not the file"
        ),
        "alphagenome_requests": 0,
        "network_requests": 0,
        "money": "none: every input was already on disk",
        "code_cleanliness": mf.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
    }
    payload["result_manifest"] = {
        "sources": [
            {"accession": "GENCODE", "version": "v50, data/reference/gencode_v50_chr*.gff3.gz"},
        ],
        "inputs": [
            mf.input_entry(Path(f"data/results/{AMENDS}.json")),
            mf.input_entry(Path("data/results/gene_identity.json")),
            *[mf.input_entry(gi.annotation_path(c)) for c in gi.CHROMS if gi.annotation_path(c).exists()],
        ],
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "counted_over": "every gene row of every chromosome's GENCODE v50 annotation",
            "gene_name_convention": (
                "holdout.genes and element_types.gene_starts both read gene_name and fall back to the "
                "versionless gene_id; this amendment counts how often that fallback can ever be "
                "reached"
            ),
        },
        "exclusions": [
            "nothing is excluded: the three figures are counted over every gene row",
            "no count of data/results/gene_identity.json is changed by this file",
        ],
        "partitions": {"none": "this amendment partitions nothing"},
        "code_cleanliness": mf.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
    }
    path = save_result(RESULT, payload)
    print(f"wrote {path}")
    print(figures)


if __name__ == "__main__":
    main()
