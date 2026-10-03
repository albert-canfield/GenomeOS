# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register item (h)'s astrocyte-argmax question before any answer is bought.

    uv run --frozen python scripts/astroargmax_register.py

Writes data/results/astroargmax_registration.json, which is committed ALONE. Committing it
satisfies ONE of the five clauses of Albert's second approval -- `astrorun`'s
`check_astroargmax_registration_is_committed` -- and authorises nothing: the cap, the separate
ledger, today's CI green and the supervisor's "dry run reviewed" are four separate refusals that
this file does not touch.

THE ORDER OF THE PAYLOAD IS THE POINT, and it is the order things were decided in:

1. THE BANDS, relative to the control, as numbers with their reason. Before the control, because a
   band width chosen after seeing a control can be chosen to clear it.
2. THE LABEL-SET RULE in words, with every borderline call decided and reasoned.
3. THE ENUMERATION under that rule: 27 CURIEs with labels and panel row counts, three registered
   alternatives, and every excluded borderline term with its reason. Not a descendant query -- there
   is no UBERON file on this machine and `data/knowledge` has 0 tracked files.
4. THE CONTROL, computed here from the committed census as the label set's share of COMPILED RULES,
   never of track names. The panel share is printed beside it and LABELLED a panel share.
5. THE DISCLOSURE of what was seen before this was written, and when.

0 AlphaGenome requests, no network to the service, no key read, no money. No astrocyte answer
exists, so no observed rate is computed and the fields a run would fill are ABSENT rather than
empty: nothing here can be read as a preview.
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import astroargmax as aa  # noqa: E402

RESULT = "astroargmax_registration"

#: This lane's own files. Everything else uncommitted in this shared checkout belongs to a peer.
OWN_CODE = (
    "genomeos/attribution/astroargmax.py",
    "scripts/astroargmax_register.py",
    "tests/test_astroargmax.py",
)


def inputs() -> list[dict[str, Any]]:
    """The three files this registration reads, digested. One of them is git-ignored by design."""
    out: list[dict[str, Any]] = []
    if aa.CENSUS.exists():
        out.append(
            mf.input_entry(
                aa.CENSUS,
                partition="the committed genome-wide census of every compiled rule's `when: "
                "cell_type` label. The CONTROL is its per_cell block's share of 440,589 compiled "
                "rules for the label set, and nothing else; no part of it is computed here",
            )
        )
    if aa.PLAN.exists():
        out.append(
            mf.input_entry(
                aa.PLAN,
                partition="the committed AstroREG-2 dry run. Read for ONE thing: the chromosome of "
                "each of the 1,232 requests, which gives the cluster sizes the power figure is "
                "stated in. No outcome, no label and no gene of any request is read for a rate",
            )
        )
    if aa.PANEL.exists():
        out.append(
            mf.input_entry(
                aa.PANEL,
                partition="the assay panel the answer's track axis is drawn from, read for its 667 "
                "rna_seq rows' name, ontology_curie, biosample_name and data_source. GIT-IGNORED by "
                "Albert's rule that data/cache stays local, so it is DECLARED BY HASH and every "
                "figure taken from it is written into the registration",
                git_ignored=True,
                declared_by_hash_not_committed=aa.THE_PANEL_IS_DECLARED_BY_HASH,
                declared_sha256=aa.PANEL_SHA256,
                declared_bytes=aa.PANEL_BYTES,
            )
        )
    return out


def payload() -> dict[str, Any]:
    entries = inputs()
    rows = aa.panel_rows()
    cen = aa.census()
    registered = aa.enumerate_set(aa.REGISTERED_SET, rows)
    alternatives = {name: aa.enumerate_set(terms, rows) for name, terms in aa.ALTERNATIVE_SETS.items()}
    ctl = aa.control(aa.REGISTERED_SET, cen)
    sizes = aa.cluster_sizes()
    return {
        "result": RESULT,
        "date": date.today().isoformat(),
        "lane": "lane-astroargmax",
        "status": (
            "A REGISTRATION AND NOTHING ELSE. No astrocyte answer exists, so no observed rate, no "
            "interval, no verdict and no count of any arm exists or can. 0 AlphaGenome requests, no "
            "network to the service, no key read, no money"
        ),
        "registered_in_its_own_commit": (
            "this file is committed ALONE, before any of the 1,232 requests is sent and before any "
            "astrocyte answer exists. A registration in the same commit as its result is a claim "
            "resting on a sentence; this one rests on the history, and astrorun's clause "
            "`astroargmax_registration_committed` reads git and not this file's own word for it"
        ),
        "question": aa.QUESTION,
        "this_file_does_not_authorise_the_run": aa.NOT_AN_AUTHORISATION,
        "ordering_is_load_bearing": aa.ORDERING_IS_LOAD_BEARING,
        # ---------------------------------------------- 1. THE BANDS, decided before the control
        "bands": {
            "decided_before_the_control": aa.BANDS_BEFORE_THE_CONTROL,
            "ratio_enriched_floor": aa.RATIO_ENRICHED_FLOOR,
            "absolute_level_floor": aa.ABSOLUTE_LEVEL_FLOOR,
            "equivalence_band_relative": aa.EQUIVALENCE_BAND,
            "equivalence_interval_on_the_ratio": [
                round(1 / (1 + aa.EQUIVALENCE_BAND), 6),
                1 + aa.EQUIVALENCE_BAND,
            ],
            "ratio_depleted_ceiling": aa.RATIO_DEPLETED_CEILING,
            "why_relative_and_not_absolute": aa.WHY_RELATIVE_AND_NOT_ABSOLUTE,
            "why_a_ratio_floor_is_not_enough": aa.WHY_A_RATIO_FLOOR_IS_NOT_ENOUGH,
            "why_two_point_zero": aa.WHY_TWO_POINT_ZERO,
            "minimum_clusters_for_an_interval": aa.MIN_CLUSTERS,
            "locus_span_when_chromosomes_are_too_few": aa.LOCUS_SPAN,
            "bootstrap_draws": aa.BOOTSTRAP_DRAWS,
            "bootstrap_seed": aa.BOOTSTRAP_SEED,
            "minimum_distinct_resamples": aa.MIN_DISTINCT_RESAMPLES,
            "readings": aa.READINGS,
            "instrument_assumptions_printed_beside_every_interval": aa.INSTRUMENT_ASSUMPTIONS,
            "every_ratio_carries_its_absolute_level": aa.EVERY_RATIO_CARRIES_ITS_LEVEL,
        },
        # ------------------------------------- 2. THE LABEL-SET RULE, in words before any CURIE
        "label_set_rule": aa.LABEL_SET_RULE,
        "borderline_calls": aa.BORDERLINE_CALLS,
        "not_by_ontology_descendant_resolution": aa.NOT_BY_DESCENDANT_RESOLUTION,
        # --------------------------------------------------------------- 3. THE ENUMERATION
        "enumeration_preceded_the_control": aa.ENUMERATION_PRECEDED_THE_CONTROL,
        "which_set_decides": aa.WHICH_SET_DECIDES,
        "label_set_CNS_PRIMARY": registered,
        "label_sets_registered_as_alternatives": alternatives,
        "panel": aa.panel_summary(rows),
        "panel_declared_by_hash": {
            "path": aa.PANEL_RELATIVE,
            "sha256": aa.PANEL_SHA256,
            "bytes": aa.PANEL_BYTES,
            "git_ignored": True,
            "a_fresh_checkout_cannot_recompute_this": aa.THE_PANEL_IS_DECLARED_BY_HASH,
        },
        # ------------------------------------------------ 4. THE CONTROL: compiled rules only
        "control": ctl,
        "control_for_each_alternative_set": {
            name: aa.control(terms, cen) for name, terms in aa.ALTERNATIVE_SETS.items()
        },
        "census": {
            "source": aa.CENSUS_RELATIVE,
            "rules": cen["rules"],
            "distinct_labels": cen["distinct_labels"],
            "per_cell_sums_to_rules": True,
            "checked_here": "the sum is recomputed by astroargmax.census() and this registration "
            "refuses to be written if it does not equal the census's own `rules`",
        },
        "what_a_rate_at_the_control_means": aa.WHAT_A_RATE_AT_THE_CONTROL_MEANS,
        "panel_share_beside_the_control": aa.panel_substring_crosscheck(aa.REGISTERED_SET, rows),
        "power_stated_in_advance": aa.power_in_cluster_units(ctl["rate"], sizes),
        "cluster_sizes": sizes,
        "design": aa.design(sizes),
        "comparison_b_the_K562_arm": aa.K562_ARM,
        "training_exposure": aa.TRAINING_EXPOSURE,
        # ------------------------------------------------------------------- 5. THE DISCLOSURE
        "disclosure": aa.DISCLOSURE,
        "cannot_establish": list(aa.CANNOT_ESTABLISH),
        "refusals": list(aa.REFUSALS),
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
                    "accession": "the committed AstroREG-2 request plan, whose `requests` list IS the 1,232",
                    "version": "the committed data/results/astroreg2_request_plan.json on this "
                    "machine, digested in inputs",
                },
                {
                    "accession": "AlphaGenome's own track metadata, as saved by the model client",
                    "version": f"the local copy at {aa.PANEL_RELATIVE}, sha256 {aa.PANEL_SHA256}, "
                    f"{aa.PANEL_BYTES} bytes. GIT-IGNORED and declared by that hash",
                },
            ],
            "inputs": entries,
            "input_count": len(entries),
            "assembly": "n/a: this registration reads no coordinate and no sequence",
            "coordinates": "n/a: no interval is read. The only genomic quantity read is each of the "
            "1,232 requests' CHROMOSOME, which is the cluster label and not a coordinate",
            "parameters": {
                "ratio_enriched_floor": aa.RATIO_ENRICHED_FLOOR,
                "absolute_level_floor": aa.ABSOLUTE_LEVEL_FLOOR,
                "equivalence_band_relative": aa.EQUIVALENCE_BAND,
                "ratio_depleted_ceiling": aa.RATIO_DEPLETED_CEILING,
                "minimum_clusters_for_an_interval": aa.MIN_CLUSTERS,
                "locus_span_when_chromosomes_are_too_few": aa.LOCUS_SPAN,
                "bootstrap_draws": aa.BOOTSTRAP_DRAWS,
                "bootstrap_seed": aa.BOOTSTRAP_SEED,
                "minimum_distinct_resamples": aa.MIN_DISTINCT_RESAMPLES,
                "label_set_CNS_PRIMARY": [c for c, _ in aa.REGISTERED_SET],
                "label_sets_registered_as_alternatives": {
                    name: [c for c, _ in terms] for name, terms in aa.ALTERNATIVE_SETS.items()
                },
                "panel_share_patterns": list(aa.PANEL_SHARE_PATTERNS),
                "panel_share_fields": list(aa.PANEL_SHARE_FIELDS),
                "no_threshold_of_any_other_lane_is_set_or_moved": "astrorun's 1,232 cap, its five "
                "approval clauses, its CI date floor, its disk floor and every registration this "
                "lane cites stay exactly as committed. The only numbers this lane sets are its own "
                "reading bands, and they are set before any figure of the question exists",
            },
            "exclusions": [
                "every panel row whose `output` is not rna_seq: 2,174 of 2,841. The gene-axis argmax "
                "this question is about is taken over rna_seq columns",
                "every borderline term decided OUT in borderline_calls, each with its reason there",
                "every astrocyte answer: none exists. No request is sent and none is read",
                "every winning-name distribution of the all-element sweep: none is read, so no "
                "K562-arm rate exists in this file",
            ],
            "partitions": {
                "the_1232_answers": "the question's denominator, one per paid astrocyte element. NOT "
                "the 667 track columns and NOT the 371 names: those are the response's AXIS",
                "the_census_440589": "the control's denominator, read from the committed census and "
                "computed by no code of this lane",
                "the_667_rna_seq_rows": "the PANEL, the axis the argmax is taken over. Used for the "
                "enumeration's row counts and the labelled panel share, never as a denominator for "
                "the question",
            },
            "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
        },
    }


def main() -> int:
    from genomeos.results import save_result

    body = payload()
    path = save_result(RESULT, body)
    print(f"registered: {path}")
    print(
        f"  bands (before the control): ratio floor {aa.RATIO_ENRICHED_FLOOR}, level floor "
        f"{aa.ABSOLUTE_LEVEL_FLOOR}, equivalence +/-{aa.EQUIVALENCE_BAND} relative"
    )
    s = body["label_set_CNS_PRIMARY"]
    print(
        f"  label set CNS_PRIMARY: {s['terms']} terms, {s['panel_rows']} panel rows "
        f"({s['panel_row_share']}), prefixes {s['prefixes_represented']}"
    )
    c = body["control"]
    print(
        f"  CONTROL (share of COMPILED RULES): {c['rules_with_a_label_in_the_set']} / "
        f"{c['rules_total']} = {c['rate']} over {c['census_labels_matched_count']} census labels"
    )
    p = body["panel_share_beside_the_control"]
    print(
        f"  PANEL SHARE (a panel share, not the control): {p['matching_rows']} / {p['panel_rows']} "
        f"= {p['panel_share']}"
    )
    print("  requests this registers nothing about: 0 sent, 1232 unbought")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
