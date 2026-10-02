# SPDX-License-Identifier: AGPL-3.0-or-later
"""Count the label concordance the registration at 47835d1, amended at AMENDMENT_1, asked for.

    uv run --frozen python scripts/label_gene_count.py

Reads data/organisms/human/noncoding_chr21.bio line by line - one committed element table, no bulk
chromosome archive, no cached-element loader - and the GTEx v8 median-TPM matrix for the SECONDARY
arm only. Both declared input digests are re-checked and a mismatch RAISES, so a result can never
be read against a different input than the registered one. Peak RSS is checked against the
registered 4 GiB ceiling and the check RAISES rather than warns.

What it reports, in the registration's own terms: the structural census this lane verified rather
than took on trust; control 1 as the sum of p squared over this population's own label
distribution; the primary within-gene pair concordance; the distance-matched different-gene arm;
the three deciding ratio intervals with their bands, their identical-resample shares and the
instrument's assumptions printed beside each; which of (a)/(a')/(b)/(c) the comparison reads;
power in genes; AMENDMENT_1's printed concentration and gene-equal-weight sensitivity with the
binding disagreement sentence; and the GTEx secondary arm under its training-exposure caveat.

0 model requests, no network, no AlphaGenome request of any kind, no money.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import label_gene as lg  # noqa: E402

RESULT = "label_gene"
REGISTRATION = Path("data/results/label_gene_registration.json")
AMENDMENT = Path("data/results/label_gene_registration_amendment_1.json")
BIO = Path("data/organisms/human/noncoding_chr21.bio")
GTEX = Path("data/cache/c5/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz")

OWN_CODE = (
    "genomeos/attribution/label_gene.py",
    "scripts/label_gene_register.py",
    "scripts/label_gene_amend_1.py",
    "scripts/label_gene_count.py",
    "tests/test_label_gene.py",
)

#: What the coordinator measured and asked this lane to verify rather than adopt.
COORDINATOR_FIGURES = {
    "rule_lines": 5176,
    "predicted_rules": 5174,
    "experimental_rules": 2,
    "distinct_elements": 5176,
    "genes": 193,
    "genes_with_two_or_more_elements": 180,
}


def recheck_inputs() -> dict[str, str]:
    """Re-digest both declared inputs and RAISE on a mismatch with the registration."""
    declared = {}
    if REGISTRATION.exists():
        import json

        declared = json.loads(REGISTRATION.read_text(encoding="utf-8"))["declared_inputs"]
    seen = {}
    for path in (BIO, GTEX):
        digest, _, _ = mf.sha256_of(path)
        seen[str(path)] = digest
        want = declared.get(str(path))
        if want and want != digest:
            raise ValueError(
                f"{path} digests {digest} but the registration declared {want}: a result read "
                "against an input other than the registered one is refused here"
            )
    return seen


def structure(rows: list[lg.RuleRow], census: dict[str, int]) -> dict[str, Any]:
    """The structural facts, verified by this lane, beside the coordinator's own figures."""
    identity = lg.verify_one_gene_per_element(rows)
    per_gene: dict[str, int] = {}
    for row in rows:
        per_gene[row.gene] = per_gene.get(row.gene, 0) + 1
    with_two = sum(1 for n in per_gene.values() if n >= 2)
    mine = {
        "rule_lines": census["rule_lines"],
        "predicted_rules": census["predicted_rules"],
        "experimental_rules": census["experimental_rules"],
        "distinct_elements": identity["distinct_elements"],
        "genes": len(per_gene),
        "genes_with_two_or_more_elements": with_two,
    }
    disagreements = {
        key: {"coordinator": COORDINATOR_FIGURES[key], "verified_here": mine[key]}
        for key in COORDINATOR_FIGURES
        if COORDINATOR_FIGURES[key] != mine[key]
    }
    return {
        "verified_here": mine,
        "coordinator_figures": COORDINATOR_FIGURES,
        "disagreements": disagreements,
        "element_blocks_in_the_file": census["element_blocks"],
        "elements_per_gene_mean": round(census["predicted_rules"] / max(1, len(per_gene)), 4),
        "elements_per_gene_max": max(per_gene.values()) if per_gene else 0,
        "every_element_appears_exactly_once": True,
        "symmetric_test_cannot_be_run": lg.registration()["cannot_be_answered"],
    }


def payload() -> dict[str, Any]:
    digests = recheck_inputs()
    rows, census = lg.parse_population(str(BIO))
    lg.check_rss()
    structural = structure(rows, census)
    control1 = lg.control_one(rows)
    labels = lg.label_distribution(rows)
    counts = lg.aggregate_pairs(rows)
    lg.check_rss()
    bins = lg.contributing_bins(counts)
    keep = bins["contributing"]
    point = lg.point_estimates(counts, keep)
    conc = lg.concentration(counts)
    draws = lg.resample(counts, keep)
    lg.check_rss()
    series = draws["series"]

    primary_vs_control = lg.ratio_interval(series["within_primary"], control1, "within_primary / control_1")
    within_vs_control = lg.ratio_interval(series["within_matched"], control1, "within_matched / control_1")
    diff_vs_control = lg.ratio_interval(
        series["different_matched"], control1, "different_gene_matched / control_1"
    )
    within_vs_diff = lg.ratio_interval(
        series["within_matched"], series["different_matched"], "within_matched / different_matched"
    )

    # Every registered route to a reading is evaluated and ALL of them are reported, so the verdict
    # cannot be the one that happened to be tested first.
    band_reading, band_words = lg.read_the_comparison(
        within_vs_diff["band"], within_vs_control["band"], diff_vs_control["band"]
    )
    floor_fired = bins["inconclusive_by_rule"] or draws["clusters"] < lg.MIN_CLUSTERS
    floor_words = (
        None
        if not floor_fired
        else "INCONCLUSIVE BY RULE - THE DATA CANNOT TELL. "
        + (
            f"the dropped within-gene matching weight {bins['dropped_within_gene_weight']} "
            f"exceeds {lg.MAX_DROPPED_WEIGHT}"
            if bins["inconclusive_by_rule"]
            else f"{draws['clusters']} gene clusters is below the floor of {lg.MIN_CLUSTERS}"
        )
    )
    power_fired = any(r["wider_than_the_band"] for r in (within_vs_diff, within_vs_control, diff_vs_control))
    power_words = (
        None
        if not power_fired
        else "INCONCLUSIVE BY CONSTRUCTION - THE DATA CANNOT TELL. A deciding interval's achieved "
        f"relative half-width exceeds the band's own half-width of "
        f"{within_vs_diff['band_half_width']}, so by the registered power rule the comparison "
        "cannot separate 'at control' from 'far above' whatever its point estimate."
    )
    # The weakest reading any registered route selects is the one reported: a lane may not pick
    # whichever route reads better, and (c) is the weaker reading wherever it is reached.
    if floor_fired or power_fired or band_reading == "c":
        reading = "c"
        words = " ALSO: ".join(
            w for w in (band_words if band_reading == "c" else None, floor_words, power_words) if w
        )
    else:
        reading, words = band_reading, band_words
    links = {
        "within_matched_far_above_different_matched": within_vs_diff["band"] == "far_above",
        "within_matched_far_above_control_1": within_vs_control["band"] == "far_above",
        "different_matched_position_against_control_1": diff_vs_control["band"],
        "what_is_established_and_what_is_not": (
            "stated as a fact about the three intervals and NOT as a reading, because the reading "
            "is (c) and may not be read upward from here: the two conditions clause (a) names are "
            "both at `far_above`, while the matched different-gene arm's own position against "
            "control 1 crosses a band edge, so the data cannot tell (a-i) from (a-ii) and the "
            "registered procedure returns (c)."
            if (
                within_vs_diff["band"] == "far_above"
                and within_vs_control["band"] == "far_above"
                and reading == "c"
            )
            else "the three intervals' bands are printed above; no link is claimed beyond them"
        ),
        "a_defect_in_this_lanes_own_registration_disclosed_rather_than_resolved": (
            "TWO registered clauses pull against each other on this outcome and neither is "
            "softened here. (i) READINGS['a'] names two conditions and says its sub-cases "
            "'neither denies (a)'; (ii) READINGS['c'] says ANY deciding interval crossing a band "
            "edge reads (c), and the registered procedure read_the_comparison - committed at "
            "e13992a before any count - tests (ii) first. The procedure therefore returns (c) and "
            "this lane reports (c), taking the WEAKER reading by rule rather than the clause that "
            "reads better. Separately, the registered POWER rule is miscalibrated for a ratio far "
            "from 1: it compares an achieved half-width against the band's half-width of 0.25 "
            "without reference to where the interval sits, so it can fire on an interval nowhere "
            "near a band edge. It is applied as registered and not weakened; it is named here as "
            "a defect for the owner, and it changes nothing, because the band rule reaches (c) on "
            "its own."
        ),
    }

    equal = lg.gene_equal_weight(counts, keep)
    equal_within = lg.ratio_interval(
        equal["series"]["within"], control1, "within_gene_equal_weight / control_1"
    )
    equal_diff = lg.ratio_interval(
        equal["series"]["different"], control1, "different_gene_equal_weight / control_1"
    )
    equal_point = (
        None if equal["within_gene_equal_weight"] is None else equal["within_gene_equal_weight"] / control1
    )
    disagree = lg.disagreement(within_vs_control, equal_within, equal_point)
    equal.pop("series", None)

    gtex = lg.gtex_top_tissue(str(GTEX), {row.gene for row in rows})
    secondary = lg.gtex_arm(rows, gtex, control1)
    peak = lg.check_rss()

    ranked = sorted(labels.items(), key=lambda kv: -kv[1])
    body: dict[str, Any] = {
        "result": RESULT,
        "date": "2026-10-02",
        "lane": "lane-labelgene",
        "question": lg.QUESTION,
        "registered_at": (
            "data/results/label_gene_registration.json, committed ALONE at 47835d1 (design code "
            "e13992a) before any concordance existed, and amended additively and still pre-count "
            "by data/results/label_gene_registration_amendment_1.json"
        ),
        "structure": structural,
        "control_1": {
            "value": round(control1, 6),
            "rule": lg.CONTROL_ONE,
            "is_a_committed_constant": lg.CONTROL_ONE_IS_A_COMMITTED_CONSTANT,
            "distinct_labels": len(labels),
            "one_over_the_label_count_for_comparison_only": round(1 / len(labels), 6),
            "top_labels": {name: n for name, n in ranked[:10]},
        },
        "primary": {"rule": lg.PRIMARY, **point},
        "concentration": conc,
        "control_2": {"rule": lg.CONTROL_TWO, "matching_rule": lg.MATCHING_RULE, **bins},
        "deciding_intervals": {
            "within_matched_over_different_matched": within_vs_diff,
            "within_matched_over_control_1": within_vs_control,
            "different_matched_over_control_1": diff_vs_control,
        },
        "primary_ratio_interval": primary_vs_control,
        "resampling": {k: v for k, v in draws.items() if k != "series"},
        "reading": {
            "which": reading,
            "in_the_registrations_own_words": words,
            "selected_by": (
                "the three deciding intervals above, by the registered band rule and the "
                "registered floor and power rules, never by a point estimate. Every route is "
                "evaluated and all of them are reported, so the verdict is not whichever route "
                "was tested first"
            ),
            "every_registered_route": {
                "band_rule": {"reads": band_reading, "words": band_words},
                "floor_rule": {"fired": bool(floor_fired), "words": floor_words},
                "power_rule": {"fired": bool(power_fired), "words": power_words},
                "the_weakest_reading_any_route_selects_is_reported": (
                    "(c) wherever it is reached: a lane may not pick the route that reads better"
                ),
            },
            "links": links,
            "all_four_readings_as_registered": lg.READINGS,
        },
        "power": {
            "rule": lg.POWER_RULE,
            "gene_clusters": draws["clusters"],
            "gene_clusters_informative_for_the_within_gene_arm": conc["genes_with_at_least_one_pair"],
            "cluster_floor": lg.MIN_CLUSTERS,
            "band_half_width": within_vs_diff["band_half_width"],
            "achieved_relative_half_width": {
                "within_matched_over_different_matched": within_vs_diff["relative_half_width"],
                "within_matched_over_control_1": within_vs_control["relative_half_width"],
                "different_matched_over_control_1": diff_vs_control["relative_half_width"],
            },
            "wider_than_the_band": {
                "within_matched_over_different_matched": within_vs_diff["wider_than_the_band"],
                "within_matched_over_control_1": within_vs_control["wider_than_the_band"],
                "different_matched_over_control_1": diff_vs_control["wider_than_the_band"],
            },
        },
        "amendment_1_sensitivity": {
            **equal,
            "within_over_control_1": equal_within,
            "different_over_control_1": equal_diff,
            "disagreement": disagree,
        },
        "secondary_gtex": secondary,
        "degeneracy": {
            "rule": lg.DEGENERACY_RULE,
            "within_gene_pair_counts_degenerate": lg.bootstrap_degenerate(
                point["within_gene_concordant_pairs"], point["within_gene_pairs"]
            ),
            "different_gene_pair_counts_degenerate": lg.bootstrap_degenerate(
                point["different_gene_concordant_pairs"], point["different_gene_pairs"]
            ),
            "identical_resample_share_per_interval": {
                "within_matched_over_different_matched": within_vs_diff["identical_resample_share"],
                "within_matched_over_control_1": within_vs_control["identical_resample_share"],
                "different_matched_over_control_1": diff_vs_control["identical_resample_share"],
                "within_primary_over_control_1": primary_vs_control["identical_resample_share"],
                "equal_weight_within_over_control_1": equal_within["identical_resample_share"],
                "equal_weight_different_over_control_1": equal_diff["identical_resample_share"],
                "secondary_gtex": secondary["identical_resample_share"],
            },
            "wilson_on_the_cluster_count_if_degenerate": lg.wilson(
                round(
                    (point["within_gene_concordance_primary"] or 0.0) * conc["genes_with_at_least_one_pair"]
                ),
                conc["genes_with_at_least_one_pair"],
            ),
            "the_wilson_bound_is_only_the_reading_when_degenerate": (
                "printed unconditionally so a reader can see it; it is the DECIDING bound only "
                "when the clustered bootstrap is degenerate, which the two flags above state"
            ),
        },
        "model_consequence_proposal_only": lg.MODEL_CONSEQUENCE,
        "cannot_be_answered": lg.registration()["cannot_be_answered"],
        "peak_rss_bytes": peak,
        "rss_ceiling_bytes": lg.RSS_CEILING_BYTES,
        "rss_rule": lg.RSS_RULE,
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
        "never_weakened": lg.NEVER_WEAKENED,
    }
    body["result_manifest"] = {
        "sources": [
            {
                "accession": "the committed compiled program data/organisms/human/"
                "noncoding_chr21.bio, generated by `genomeos budget --chrom chr21 --bio`",
                "version": f"sha256 {digests[str(BIO)]}, re-checked against the registration",
            },
            {
                "accession": "GTEx v8 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm",
                "version": f"sha256 {digests[str(GTEX)]}, SECONDARY arm only",
            },
        ],
        "inputs": [
            mf.input_entry(BIO, partition="the whole committed chr21 non-coding program"),
            mf.input_entry(GTEX, partition="GTEx v8 median TPM, SECONDARY arm only"),
            *(
                [mf.input_entry(REGISTRATION, partition="this lane's registration at 47835d1")]
                if REGISTRATION.exists()
                else []
            ),
            *(
                [mf.input_entry(AMENDMENT, partition="AMENDMENT_1, additive and pre-count")]
                if AMENDMENT.exists()
                else []
            ),
        ],
        "input_count": 2 + REGISTRATION.exists() + AMENDMENT.exists(),
        "assembly": "GRCh38: the assembly the element loci are on",
        "coordinates": {
            "base": 0,
            "interval": "half-open",
            "note": "the element blocks' own `locus: chr21:start-end`, read only for the midpoint "
            "that defines element-to-element distance",
        },
        "parameters": {
            "band_at_control": list(lg.BAND_AT_CONTROL),
            "band_far_above": lg.BAND_FAR_ABOVE,
            "distance_bin_edges_bp": [
                int(e) if e != float("inf") else "inf" for e in lg.DISTANCE_BIN_EDGES_BP
            ],
            "minimum_pairs_per_bin": lg.MIN_PAIRS_PER_BIN,
            "maximum_dropped_matching_weight": lg.MAX_DROPPED_WEIGHT,
            "minimum_clusters_genes": lg.MIN_CLUSTERS,
            "bootstraps": lg.BOOTSTRAPS,
            "seed": lg.SEED,
            "rss_ceiling_bytes": lg.RSS_CEILING_BYTES,
            "every_one_of_them_fixed_before_this_run": "at 47835d1 and AMENDMENT_1; this script "
            "sets none of them",
        },
        "exclusions": [
            "the two rule lines carrying `experimental`",
            "every `_measured` element block",
            "every gene carrying one element, from the within-gene arm only",
            "every distance bin below the registered pair minimum in either arm",
            "every rule whose normalised label matches no GTEx column, from the SECONDARY denominator only",
            "every bulk chromosome archive and every cached-element loader call: none is made",
        ],
        "partitions": {
            "within_gene_pairs": str(point["within_gene_pairs"]),
            "different_gene_pairs": str(point["different_gene_pairs"]),
            "gene_clusters": str(draws["clusters"]),
            "gtex_secondary_rules": str(secondary["rules_in_the_arm"]),
        },
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }
    return body


def main() -> int:
    from genomeos.results import save_result

    body = payload()
    path = save_result(RESULT, body)
    print(f"wrote {path}")
    s = body["structure"]["verified_here"]
    print(
        f"  verified: {s['predicted_rules']} predicted rules, {s['experimental_rules']} "
        f"experimental, {s['distinct_elements']} distinct elements, {s['genes']} genes, "
        f"{s['genes_with_two_or_more_elements']} with two or more"
    )
    if body["structure"]["disagreements"]:
        print(f"  DISAGREEMENTS with the coordinator: {body['structure']['disagreements']}")
    print(f"  control 1 (sum of p squared): {body['control_1']['value']}")
    print(
        f"  within-gene concordance (primary, pair-pooled): "
        f"{body['primary']['within_gene_concordance_primary']} over "
        f"{body['primary']['within_gene_pairs']} pairs of "
        f"{s['genes_with_two_or_more_elements']} genes"
    )
    print(
        f"  within matched {body['primary']['within_gene_concordance_matched']} vs "
        f"different-gene matched {body['primary']['different_gene_concordance_matched']}"
    )
    for name, block in body["deciding_intervals"].items():
        print(
            f"  {name}: ci95 {block['ci95']} band {block['band']} "
            f"identical-resample share {block['identical_resample_share']}"
        )
    print(
        f"  READING ({body['reading']['which']}): {body['reading']['in_the_registrations_own_words'][:200]}"
    )
    print(
        f"  concentration: top {conc_n(body)} genes hold "
        f"{body['concentration']['share_of_pairs_from_the_top_genes']} of all within-gene pairs"
    )
    print(f"  {body['amendment_1_sensitivity']['disagreement']['statement']}")
    print(
        f"  SECONDARY GTEx: {body['secondary_gtex']['agree']}/"
        f"{body['secondary_gtex']['rules_in_the_arm']} = {body['secondary_gtex']['share']} "
        f"(chance {body['secondary_gtex']['expected_share_if_labels_were_drawn_from_this_population']})"
    )
    print(f"  peak RSS {body['peak_rss_bytes']} of ceiling {body['rss_ceiling_bytes']}")
    print("  0 model requests, no network, no money")
    return 0


def conc_n(body: dict[str, Any]) -> int:
    return int(body["concentration"]["top_n"])


if __name__ == "__main__":
    raise SystemExit(main())
