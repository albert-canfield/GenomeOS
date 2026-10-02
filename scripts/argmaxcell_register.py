# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register the argmax-cell question before a single count is taken.

    uv run --frozen python scripts/argmaxcell_register.py

Writes data/results/argmaxcell_registration.json. This file is committed ON ITS OWN, before
scripts/argmaxcell_count.py is run, because "registered first" has to rest on git and not on a
lane's word: a peer put its registration in the same commit as its result tonight and then rested
on the sentence rather than on the history.

What it binds, all of it before the observed side exists:

1. THE FALSIFIER, with its thresholds and both directions of the reading, including the branch that
   says plainly that a rate at or near the base rate means a rule's cell is NOT evidence of where it
   acts, and the asymmetry (CONFOUND) that stops the opposite branch being read as a confirmation.
2. THE BASE RATE AND HOW IT IS COMPUTED, as a number, from the COMMITTED
   data/results/context_evidence.json: per_cell[label]['rules'] / rules over 440,589 compiled rules.
   The number is in this file, so it cannot have been derived after seeing an observed figure.
3. THE POPULATION, THE DENOMINATOR AND WHAT WILL BE REPORTED.

It also records the read discipline and the RSS ceiling, because a peer paged the machine tonight by
looping a cached-element loader that falls back to whole-chromosome archives.

No request is sent, no key is read, no element table is opened and no observed figure is computed.
0 model requests, no network, no money.
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import argmaxcell as ac  # noqa: E402
from genomeos.attribution import crispri as cr  # noqa: E402
from genomeos.attribution import direction_v2 as dv2  # noqa: E402

RESULT = "argmaxcell_registration"
WILL_WRITE = "data/results/argmaxcell.json"

#: This lane's own files. Everything else uncommitted in this shared checkout belongs to a peer.
OWN_CODE = (
    "genomeos/attribution/argmaxcell.py",
    "scripts/argmaxcell_register.py",
    "scripts/argmaxcell_count.py",
    "tests/test_argmaxcell.py",
)


def inputs() -> list[dict[str, Any]]:
    """The named files the run will read, digested before it runs.

    The 24 compact element tables are NOT digested: they total 594 MB and digesting them would read
    every byte of them here, which is the heavy read this lane is registering a discipline against.
    They are named instead, by the committed summary that points at each one.
    """
    out: list[dict[str, Any]] = []
    if ac.CENSUS.exists():
        out.append(
            mf.input_entry(
                ac.CENSUS,
                partition="the committed genome-wide census of every compiled rule's `when: "
                "cell_type` label. The BASE RATE of every arm is read from its per_cell block and "
                "from nowhere else; this lane computes no part of it",
            )
        )
    for name in (cr.TRAINING, cr.HELDOUT):
        p = cr.KNOWLEDGE / name
        if p.exists():
            out.append(
                mf.input_entry(
                    p,
                    partition="a cached ENCODE CRISPRi benchmark table. Read for each valid pair's "
                    "CellType, chrom, chromStart, chromEnd and measuredGeneSymbol only: no outcome "
                    "column, no DHS or H3K27ac value and no power column is read",
                )
            )
    return out


def payload() -> dict[str, Any]:
    entries = inputs()
    cen = ac.census()
    rates = {cell: ac.base_rate(cell, cen) for cell in ac.ARMS}
    return {
        "result": RESULT,
        "date": date.today().isoformat(),
        "lane": "lane-argmaxcell",
        "registered_in_its_own_commit": (
            "this file is committed alone, before scripts/argmaxcell_count.py is run and before any "
            "count exists. A registration in the same commit as its result is a claim resting on a "
            "sentence; this one rests on the history, and the result commit names this file's sha"
        ),
        "question": ac.QUESTION,
        "argmax_rule": ac.ARGMAX_RULE,
        "argmax_is_a_selection": ac.ARGMAX_IS_A_SELECTION,
        "assigned_cell": getattr(dv2, "ASSIGNED_CELL", ""),
        # --- 1. the falsifier, before any count ---
        "falsifier": ac.FALSIFIER,
        "falsifier_thresholds": {
            "tolerance_on_the_difference": ac.TOLERANCE,
            "usable_rate": ac.USABLE,
            "minimum_elements_for_a_rate": ac.MIN_ELEMENTS,
        },
        "confound_registered_before_any_count": ac.CONFOUND,
        "cross_arm_check": ac.CROSS_ARM,
        "what_a_rate_at_the_base_rate_means": (
            "stated here so it cannot be softened later: if the observed rate is at or near the "
            "base rate, a rule's `cell` is NOT evidence of where it acts. That is the finding and "
            "it will be written in those words. It is not a failure of the measurement and it is "
            "not an absence of regulation anywhere: it is a fact about what the cell label carries"
        ),
        # --- 2. the base rate, and how it is computed, before the observed rate ---
        "base_rate_rule": ac.BASE_RATE_RULE,
        "base_rate_is_not_uniform": ac.BASE_RATE_IS_NOT_UNIFORM,
        "base_rates": rates,
        "base_rate_census": {
            "source": ac.CENSUS_RELATIVE,
            "rules": cen["rules"],
            "distinct_labels": len(cen["per_cell"]),
            "per_cell_sums_to_rules": True,
            "checked_here": "the sum is recomputed by argmaxcell.census() and the registration "
            "refuses to be written if it does not equal the census's own `rules`",
        },
        "second_base_rate": ac.SECOND_BASE_RATE,
        # --- 3. the population, the denominator, what will be reported ---
        "population": ac.POPULATION,
        "denominator": ac.DENOMINATOR,
        "denominator_is_separate": ac.DENOMINATOR_IS_SEPARATE,
        "will_report": ac.WILL_REPORT,
        "arms": list(ac.ARMS),
        "overlap_rule": ac.OVERLAP_RULE,
        "label_rule": ac.LABEL_RULE,
        "alias_rule": ac.ALIAS_RULE,
        "labels_for_each_arm": {c: list(ac.labels_for(c)) for c in ac.ARMS},
        "interval_is_binomial": ac.INTERVAL_IS_BINOMIAL,
        "grouping": ac.GROUPING,
        "counts_named_in_advance": list(ac.COUNTS_NAMED),
        "refusals": list(ac.REFUSALS),
        # --- the read discipline ---
        "read_discipline": ac.READ_DISCIPLINE,
        "rss_ceiling_bytes": ac.RSS_CEILING_BYTES,
        "no_archive_is_opened": (
            "no element answer is read through the cached-element loader and no chromosome archive "
            "is opened. That loader falls back to the archive and caches every archive it opens in "
            "a module-level dict that never evicts; a peer's inline loop over 93 elements reached "
            "23 GB RSS tonight and drove free disk under the floor. The compact per-chromosome "
            "element table is read instead, one element at a time, and argmaxcell's own source is "
            "asserted by test to name no archive loader"
        ),
        # --- what this lane may not claim ---
        "exploration_preceded_this": ac.EXPLORATION_PRECEDED_THIS,
        "edits_no_rule": ac.NO_RECOMMENDATION,
        "validates_nothing": ac.VALIDATES_NOTHING,
        "name_match_is_not_a_confirmation": ac.NAME_MATCH_IS_NOT_A_CONFIRMATION,
        "metadata_copy_limitation": ac.METADATA_COPY_LIMITATION,
        "retention_is_not_this_question": ac.RETENTION_IS_NOT_THIS_QUESTION,
        "retained_cells_today": list(ac.RETAINED_CELLS_TODAY),
        "cannot_establish": list(ac.CANNOT_ESTABLISH),
        "what_the_run_will_write": WILL_WRITE,
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
        "result_manifest": {
            "sources": [
                {
                    "accession": "the committed genome-wide context-evidence census of compiled "
                    "rule cell labels",
                    "version": "the committed data/results/context_evidence.json on this machine, "
                    "digested in inputs",
                },
                {
                    "accession": cr.PUBLISHED_SOURCE,
                    "version": "the two benchmark tables cached under data/knowledge/crispri on "
                    "this machine, digested in inputs",
                },
                {
                    "accession": "the all-enhancer deletion sweep's compact per-chromosome element "
                    "tables, each named by its committed enhancer_targets_all_<chrom> summary",
                    "version": "the tables on this machine. NOT digested here: 594 MB, and "
                    "digesting them would be the heavy read this registration forbids",
                },
            ],
            "inputs": entries,
            "input_count": len(entries),
            "assembly": "GRCh38: the assembly both benchmark tables and the deletion sweep are on. "
            "No coordinate of this registration's own is in any build",
            "coordinates": {
                "base": 0,
                "interval": "half-open",
                "note": "the benchmark tables' own chromStart/chromEnd and the element tables' own "
                "start/end, both 0-based half-open, used only for the overlap test in OVERLAP_RULE. "
                "This registration itself reads no interval",
            },
            "parameters": {
                "tolerance": ac.TOLERANCE,
                "usable_rate": ac.USABLE,
                "minimum_elements_for_a_rate": ac.MIN_ELEMENTS,
                "reach": cr.REACH,
                "arms": list(ac.ARMS),
                "rss_ceiling_bytes": ac.RSS_CEILING_BYTES,
                "no_threshold_of_any_rule_is_set_or_moved": "MIN_EFFECT, the direction clauses and "
                "every magnitude floor stay exactly as committed. The only numbers this lane sets "
                "are its own reading thresholds, and they are set here before any figure exists",
            },
            "exclusions": [
                "every benchmark pair whose ValidConnection is not TRUE: crispri.parse drops them, "
                "which also drops the promoter and exon overlaps",
                "every element with no predicted coding target: it emits no compiled rule, so it is "
                "in neither the base rate nor any arm",
                "every per-element archive and every cached-element loader call: none is made",
                "every outcome, DHS, H3K27ac and power column of the benchmark tables",
            ],
            "partitions": {
                "per_arm_elements": "the DISTINCT attributed elements overlapping that arm's pairs: "
                "this lane's headline denominator, one per arm",
                "per_arm_pairs": "that arm's valid pairs: the secondary denominator, never the headline",
                "the_census_440589": "the base rate's denominator, read from the committed census "
                "and computed by no code of this lane",
            },
            "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
        },
    }


def main() -> int:
    """Every read happens BEFORE save_result closes the trace window, so what is reconciled is this
    run's own reads. Printing from the payload rather than re-reading the census is the reason."""
    from genomeos.results import save_result

    body = payload()
    path = save_result(RESULT, body)
    print(f"registered: {path}")
    for cell, r in body["base_rates"].items():
        print(f"  base rate {cell:18} {r['rules_with_this_label']:>7} / {r['rules_total']} = {r['rate']}")
    print(f"the run will write: {WILL_WRITE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
