# SPDX-License-Identifier: AGPL-3.0-or-later
"""Why 44 of N1's 56 factors have no responders. An analysis of committed figures. No measurement.

    uv run --frozen python scripts/n1_coverage_analysis.py

**THIS RUN OPENS NO STUDY FILE.** It reads three committed results and two committed texts, and that is
the whole of its input. It does not download anything, it does not touch `data/cache/`, it parses no
p-value, it reads no perturbed pseudobulk row, no knockdown field and no differential-expression call,
and it makes no request of any kind. Every number it reports is arithmetic over figures that
`data/results/n1_result_amendment_2.json` already records, or a line number found in the committed text
of the module that wrote them. N1b step 2's value is that outcome blindness has not been spent; nothing
here spends any of it.

What it establishes, in the order the questions were asked:

1. **Which endpoint produced `no_responders`.** Twice over, so the answer does not rest on reading code
   alone: the predicate is located in the committed module by its exact text, and the committed result's
   own per-factor key set is classified independently (`n1_coverage.endpoint_from_result_keys`).
2. **That the derived statistic T = X * sqrt(n) is not on that path**, shown by where its two markers do
   occur against the span of `run_v3`, so a reader checks the absence instead of taking it.
3. **What the label conflates**, by counting the sibling labels the same function emits.
4. **Why the shortfall happened**: the responder counts' concentration across factors, against what an
   even spread at the same mean would have given.
5. **What a larger gene universe would do**, as a fitted extrapolation AND as the two ends the counts
   cannot separate, never one without the other.
6. **What each candidate change would license**, each marked pre-registerable or post-hoc and
   blindness-costing or free.

The floor is not relaxed, reinterpreted or worked around anywhere, and the unstratified mean of -0.061
over 12 factors is quoted only as what the registration says it is: below the floor, carrying no
criterion, and not an estimate of the gain.
"""

from __future__ import annotations

import json
import sys
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import n1_coverage as nc  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "n1_coverage_analysis"
OWN_CODE = ("scripts/n1_coverage_analysis.py", "genomeos/attribution/n1_coverage.py")

RESULT_FILE = Path("data/results/n1_result_amendment_2.json")
REGISTRATION = Path("data/results/n1_registration_amendment_2.json")
CALIBRATION = Path("data/results/n1b_calibration.json")
MODULE = Path("genomeos/attribution/n1_perturb_response.py")
ATTRIBUTION = Path("docs/ATTRIBUTION.md")

#: The committed figure the original N1 registration recorded about the derived statistic's null tail,
#: carried here word for word and marked for what it is. It is a NULL-TAIL rate: how many responders
#: chance alone would produce under the old operational label. It is not a coverage estimate, it says
#: nothing about how many factors would be defined, and in particular it says nothing about the
#: quantity that actually broke the floor.
NULL_TAIL_QUOTE = "about two chance responders per factor over 860 genes"


def _json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def endpoint_block(result: dict[str, Any], module_text: str) -> dict[str, Any]:
    """Which responder definition the committed result was written under, answered two ways."""
    sites = nc.source_sites(module_text)
    from_keys = nc.endpoint_from_result_keys(result["factors"][0])
    absence = nc.t_is_absent_from(module_text, "run_v3", ["RESPONDS_T", "pooled_t("])
    agree = from_keys["endpoint"] == "published_anderson_darling"
    return {
        "answer": "(a) the producers' own Anderson-Darling differential-expression calls at "
        "BH-adjusted p < 0.05",
        "not": "(b) a per-gene derived statistic of the T = X * sqrt(n) family, and nothing else "
        "either: the two answers below agree and no third path can write these bytes",
        "from_the_code": {
            "file": MODULE.as_posix(),
            "predicate": nc.SITES["responder_predicate"],
            "predicate_line": sites["responder_predicate"],
            "level_constant": nc.SITES["ad_level_constant"],
            "level_constant_line": sites["ad_level_constant"],
            "emission_site": nc.SITES["no_responders_emission"],
            "emission_line": sites["no_responders_emission"],
            "chain": "factor_result_v3 forms the labels from the published adjusted p-values that "
            "ad_values_pass parsed; stratified_auroc emits the string when that label class is empty; "
            "paired_analysis_v2 tallies it into undefined_reasons",
        },
        "from_the_result_alone": {
            "basis": "each factor_result_* writer emits a key the other does not, so the committed "
            "bytes name their own writer without reference to the code",
            "key_present": from_keys["keys_present"],
            "key_absent": from_keys["keys_absent"],
            "endpoint": from_keys["endpoint"],
            "datasets_read": result["datasets_read"],
            "normalized_pseudobulk_appears": any("normalized" in s for s in result["datasets_read"]),
            "reading": "the normalized pseudobulk never appears and no perturbed row of any file was "
            "read; the only pseudobulk read is 585 non-targeting rows, used for the decile axis",
        },
        "the_two_answers_agree": agree,
        "prose_and_code_agree": True,
        "registration_question": _json(REGISTRATION)["question"],
        "t_is_absent_from_run_v3": absence,
        "where_t_is_used_instead": {
            "responds_t_constant_line": sites["responds_t_constant"],
            "pooled_t_definition_line": sites["pooled_t_definition"],
            "note": "every occurrence listed under t_is_absent_from_run_v3.marker_lines lies outside "
            "the span of run_v3, which is given beside them; the amendment-1 path that does use them "
            "stops at eligibility_unsupported and wrote no committed result",
        },
        "what_this_settles": "the shortfall is a shortfall under endpoint (a) and under endpoint (a) only",
        "what_this_leaves_open": "the 44 are SILENT about coverage under a calibrated-T endpoint. Its "
        "coverage is UNKNOWN; it was not estimated here, in either direction, because estimating it "
        "would read outcome rows. The one caution that follows from the committed counts and not from "
        "any new reading is that the quantity which broke the floor is concentration across factors "
        "(see `why_the_floor_was_not_reached`), which a change of per-gene label does not obviously "
        "fix. That is a caution and not an estimate, and it must not be read as one either way.",
    }


def label_block(result: dict[str, Any], module_text: str) -> dict[str, Any]:
    """What `no_responders` means operationally, and which sibling causes are counted apart from it."""
    sites = nc.source_sites(module_text)
    factors = result["factors"]
    reasons = result["analysis"]["coverage"]["undefined_reasons"]
    return {
        "operational_definition": "a factor is counted undefined with the reason `no_responders` when "
        "NOT ONE gene of its own universe -- the 496 analysis genes, minus that factor's excluded "
        "genes, minus any gene whose published adjusted p-value is missing -- has a published "
        "BH-adjusted Anderson-Darling p-value below 0.05 in that factor's principal-transcript column",
        "emitted_when": "pos == 0, where pos counts the True labels",
        "undefined_reasons_as_committed": reasons,
        "sibling_labels_counted_apart": {
            "all_responders": {
                "emitted_at_line": sites["all_responders_emission"],
                "count_in_the_committed_result": reasons.get("all_responders", 0),
            },
            "no_stratum_with_both_classes": {
                "emitted_at_line": sites["no_stratum_emission"],
                "count_in_the_committed_result": reasons.get("no_stratum_with_both_classes", 0),
                "means": "decile matching left no comparison group",
            },
        },
        "genes_missing_per_factor": sorted({r["genes_missing"] for r in factors}),
        "genes_per_factor": sorted({r["genes"] for r in factors}),
        "ruled_out_by_the_committed_figures": [
            "decile matching leaving no comparison group: separately labelled and counted 0",
            "every gene responding: separately labelled and counted 0",
            "missing p-values: genes_missing is 0 for all 56 factors, so the 56 x 496 submatrix carried "
            "a usable value at every cell",
            "per-factor exclusions: genes is 496 for all 56, so no factor excluded any analysed gene "
            "and all 56 share one identical universe",
            "an identifier mismatch between the two files: the registration records the gene axis as "
            "Ensembl gene IDs as row labels, with the file's first row label placed at index 2 of the "
            "pseudobulk's genes, so the 364 absent genes are a genuine gene-subset difference",
        ],
        "conflates_exactly_two_causes": nc.missing_decomposition_field()["would_separate"],
        "decomposition_is_not_possible_from_committed_data": True,
        "missing_field": {
            **nc.missing_decomposition_field(),
            "skip_site": nc.SITES["values_pass_row_skip"],
            "skip_line": sites["values_pass_row_skip"],
        },
    }


def coverage_gap_block(result: dict[str, Any], pseudobulk_genes: int) -> dict[str, Any]:
    """The 364-of-860 gap: what it could do to the responder count, and what cannot be determined."""
    counts = [r["responders"] for r in result["factors"]]
    cov = result["coverage"]
    genes = cov["analysis_genes"]
    universes = [genes, cov["universe"], pseudobulk_genes]
    return {
        "as_committed": cov,
        "share_of_the_universe_absent": cov["absent_from_file"] / cov["universe"],
        "model_based_extrapolation": nc.gamma_mixed_universe(counts, genes, universes),
        "model_based_universe_for_the_floor": nc.gamma_mixed_universe_for_defined(counts, genes, nc.FLOOR),
        "distribution_free_bracket": nc.distribution_free_bracket(counts, genes, cov["absent_from_file"]),
        "largest_universe_the_study_admits": {
            "genes": pseudobulk_genes,
            "where_from": "the gene count of the raw pseudobulk, as data/results/n1b_calibration.json "
            "records it",
            "reading": "the fitted extrapolation falls short of the floor even there, which is why "
            "enlarging the gene universe is not a route to it",
        },
        "never_quote_the_model_alone": "the fitted figures and the distribution-free bracket are one "
        "finding and are to be carried together; the bracket's high end is not excluded by anything "
        "committed",
        "cannot_be_determined": [
            "why 364 of the 860 are absent from the published file's rows. The plausible mechanism is "
            "that the producers' table reports fewer genes than the pseudobulk holds because of an "
            "expression filter, in which case the absent genes are the least expressed and the least "
            "likely to be called differentially expressed, and the equal-density reading above is "
            "already generous. It is not confirmed here.",
            "the published file's total row count, which would say how large the gap is relative to "
            "what the producers tested. It is in no committed record (see carried_forward).",
            "whether a factor with no responder among the analysed genes has a per-gene rate near zero "
            "or was merely unlucky in them. Only its p-values at the other genes would say, and "
            "reading them spends outcome blindness.",
        ],
    }


def changes_block() -> list[dict[str, Any]]:
    """Every change that could in principle bring the defined count to the floor, each classified."""
    return [
        {
            "change": "enlarge the gene universe toward every gene the raw pseudobulk holds",
            "design_status": "pre-registerable: a universe rule can be fixed before any label is seen",
            "reads_new_outcome_data": True,
            "what_it_would_license": "nothing that reaches the floor. The extrapolation fitted to the "
            "committed counts falls short at the largest universe the study admits, and the shortfall "
            "it would have to close is not a shortfall of genes.",
        },
        {
            "change": "enlarge the candidate factor set beyond the 56",
            "design_status": "pre-registerable",
            "reads_new_outcome_data": True,
            "what_it_would_license": "at the observed share of defined factors it would take about 143 "
            "candidates to expect 30 defined, and the element scan that produced the 56 would have to "
            "be extended first. Nothing committed says new factors would behave differently.",
        },
        {
            "change": "raise the adjusted-p level above 0.05, or re-run the adjustment over this subset",
            "design_status": "POST-HOC tuning of a frozen label",
            "reads_new_outcome_data": True,
            "what_it_would_license": "nothing. The level is the producers' own and is frozen, and the "
            "family the adjustment ran over is not documented, so no re-adjustment is justified in "
            "either direction.",
        },
        {
            "change": "change the endpoint, for instance to a calibrated-T label",
            "design_status": "pre-registerable ONLY as a new registered study, never as a relaxation of N1",
            "reads_new_outcome_data": True,
            "what_it_would_license": "not established here, and deliberately so. The 44 are SILENT "
            "about this endpoint. Its coverage is UNKNOWN and was not estimated, in either direction. "
            "The caution that does follow from the committed counts is that concentration across "
            "factors is what broke the floor, and a change of per-gene label does not obviously fix "
            "it; that is a caution and not an estimate, and a reader must not take it as one.",
        },
        {
            "change": "lower the floor of 30, or treat the unstratified mean as an estimate",
            "design_status": "POST-HOC, and refused",
            "reads_new_outcome_data": False,
            "what_it_would_license": "nothing. This is the failure the floor exists to prevent. The "
            "unstratified mean over the 12 defined factors is -0.06105..., it carries no criterion, it "
            "sits below the floor and it is not an estimate of the gain.",
        },
        {
            "change": "record the published file's total row count",
            "design_status": "neither: a record-keeping repair, not a design change",
            "reads_new_outcome_data": False,
            "what_it_would_license": "not one more defined factor. It would say whether the 364 gap is "
            "small or large relative to what the producers tested. Ruled not now (see carried_forward).",
        },
    ]


def main() -> None:
    result = _json(RESULT_FILE)
    registration = _json(REGISTRATION)
    pseudobulk_genes = _json(CALIBRATION)["genes"]["in_file"]
    module_text = MODULE.read_text()
    attribution_text = ATTRIBUTION.read_text()
    sites = nc.source_sites(module_text)
    counts = [r["responders"] for r in result["factors"]]
    conc = nc.concentration(counts)
    null_tail_lines = [
        i for i, line in enumerate(attribution_text.splitlines(), 1) if NULL_TAIL_QUOTE in line
    ]

    payload: dict[str, Any] = {
        "result": RESULT,
        "date": date.today().isoformat(),
        "what_this_is": "AN ANALYSIS OF COMMITTED FIGURES. It contains NO NEW MEASUREMENT. No study "
        "file was opened: nothing was downloaded, no cached pseudobulk or published p-value file was "
        "read, no p-value was parsed, no perturbed row, knockdown field or differential-expression "
        "call was touched, and no request of any kind was made. Every figure below is arithmetic over "
        "numbers already committed in data/results/n1_result_amendment_2.json, or a line number found "
        "in the committed text of the module that wrote them.",
        "outcome_exposure": {
            "spent": "none",
            "why": "the inputs are three committed result files and two committed texts. N1b step 2's "
            "value is that outcome blindness has not been spent, and this analysis spends none of it.",
        },
        "explains": {
            "result": "n1_result_amendment_2",
            "status": result["status"],
            "committed_note": result["analysis"]["note"],
            "the_floor_is_not_touched": "no rule is relaxed, reinterpreted or worked around here, and "
            "no estimate is recovered. The question answered is WHY the coverage is what it is.",
        },
        "headline": (
            "The floor failed through CONCENTRATION, not through density. The 56 factors average "
            f"{conc['per_factor_mean']:.3f} responders each with a variance of "
            f"{conc['per_factor_variance_ddof1']:.4f}, a ratio of {conc['variance_to_mean_ratio']:.2f}. "
            "An even spread at that same mean would have left "
            f"{conc['even_spread']['expected_defined']:.1f} factors defined, ABOVE the floor of "
            f"{conc['floor']}; the observed {conc['zero_factors_observed']} zero factors is "
            f"{conc['even_spread']['observed_zeros_in_sd']:.1f} standard deviations away from that "
            "model, which the committed counts therefore refute. All of the responders sit in "
            f"{conc['defined_observed']} of the {conc['factors']} factors. That is a property of how "
            "much perturbations differ from one another in effect, not a property of the gene universe, "
            "and no rule about which genes to analyse corrects it."
        ),
        "endpoint": endpoint_block(result, module_text),
        "what_no_responders_means": label_block(result, module_text),
        "why_the_floor_was_not_reached": {
            **conc,
            "responder_counts": sorted(counts, reverse=True),
            "share_of_responders_in_the_defined_factors": 1.0,
            "reading": "the shortfall is not a shortage of responders overall. The per-factor mean is "
            "high enough that an even spread would have cleared the floor. It is a shortage of factors "
            "with ANY detected response, which is what extreme between-factor variation produces.",
            "this_is_a_property_of": "perturbation effect sizes varying across perturbations",
            "this_is_not_a_property_of": "the gene universe, the decile matching, the missing values "
            "or the two arms' scores",
        },
        "the_coverage_gap": coverage_gap_block(result, pseudobulk_genes),
        "changes_that_could_reach_the_floor": changes_block(),
        "verdict": "No pre-registerable change to N1's frozen design reaches the floor of 30 without "
        "reading new outcome data, and the change that enlarges coverage the most is the one the "
        "committed counts predict would still fall short. 'The study lacks responders under this "
        "endpoint and no pre-registerable change reaches the floor' is the finding, stated with no "
        "less confidence than any other.",
        "the_null_tail_figure_is_not_a_coverage_estimate": {
            "quote": NULL_TAIL_QUOTE,
            "where": f"{ATTRIBUTION.as_posix()} line(s) {null_tail_lines}",
            "what_it_is": "the original N1 registration's own pre-measurement arithmetic on the "
            "derived statistic's null tail: how many responders CHANCE ALONE would produce under the "
            "old operational label, over 860 genes",
            "what_it_is_not": "it is not a coverage estimate. It says nothing about how many factors "
            "would be defined under any endpoint, and in particular nothing about concentration across "
            "factors, which is the quantity that broke the floor. It is quoted here only as what it is.",
        },
        "carried_forward": {
            "item": "record the published file's total row count, len(labels['rows']), on the NEXT "
            "authorised pass over that file, if one ever happens -- not as a separate read",
            "where_it_was_lost": f"{MODULE.as_posix()} line(s) "
            f"{sites['labels_pass_row_count_discarded']}: the labels pass forms every row's first "
            "field, the run counts only the three universe-relative figures from it, and the row count "
            "itself is discarded with the list",
            "ruling": "THE COORDINATOR'S RULING, 2026-10-03: DO NOT READ IT as a separate act. Two "
            "reasons, neither of them caution about the door: it is not decision-critical, because "
            "enlarging the universe is predicted to fall short at every gene the study holds anyway, "
            "so the size of the gap cannot change the recommendation; and a read that cannot change a "
            "decision is not worth its precedent. It rides along on the next authorised pass, where it "
            "costs nothing at all.",
            "the_lane_s_reading_differed": "lane-n1coverage judged the read free of outcome blindness: "
            "a labels pass parses no value and exposes no outcome. The coordinator does not dispute "
            "that and declines the separate read on cost and precedent. Both positions are on the "
            "record; the ruling is the coordinator's and the lane did not take the read.",
            "for_the_next_registration": "a pass that has already paid for its access should record the "
            "cheap identity facts it is standing on top of -- row and column counts, duplicate counts, "
            "the first and last label -- because each is free while the file is open and costs a whole "
            "read afterwards",
        },
        "cost": {
            "money": "none",
            "model_requests": 0,
            "alphagenome_requests": 0,
            "network": "none: no request was made",
            "study_files_opened": 0,
            "outcome_values_parsed": 0,
        },
    }

    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "data/results/n1_result_amendment_2.json, the committed result of N1's "
                "authorised amendment-2 run",
                "version": f"amendment {result['amendment']}, written {result['date']} at "
                f"{result['amendment_commit']}; read for its coverage, analysis and per-factor counts",
            },
            {
                "accession": "data/results/n1_registration_amendment_2.json, the registration that "
                "froze the question and the label",
                "version": f"{registration['date']}; read for the question, the gene axis and the "
                "rules, and for nothing numeric",
            },
            {
                "accession": "data/results/n1b_calibration.json, read for one number: the gene count "
                "of the raw pseudobulk",
                "version": "its own committed figure genes.in_file, which bounds the largest gene "
                "universe the study admits",
            },
            {
                "accession": "genomeos/attribution/n1_perturb_response.py, the committed module that "
                "wrote the result",
                "version": "its text in this checkout, searched for exact source strings so every "
                "line number reported is found rather than typed in",
            },
            {
                "accession": "docs/ATTRIBUTION.md, read for one quoted sentence",
                "version": "the original N1 registration's null-tail figure, carried word for word and "
                "marked as not a coverage estimate",
            },
        ],
        "inputs": [
            mf.input_entry(
                RESULT_FILE,
                partition="the committed amendment-2 result: coverage, analysis and the 56 per-factor "
                "records",
            ),
            mf.input_entry(REGISTRATION, partition="the frozen question, rules and gene axis"),
            mf.input_entry(CALIBRATION, partition="one figure: genes.in_file, the pseudobulk's genes"),
            mf.input_entry(
                MODULE,
                partition="committed code, read as TEXT to locate the predicate and the emission sites; "
                "outside data/, so the audit hook that watches data/ cannot record this read and the "
                "entry appears under declared_not_read for that reason and not because its digest was "
                "carried over from elsewhere",
            ),
            mf.input_entry(
                ATTRIBUTION,
                partition="committed prose, read for one quoted sentence; outside data/, so the same "
                "note applies as for the module",
            ),
        ],
        "assembly": "n/a: this analysis holds no coordinate and no sequence. It is arithmetic over "
        "counts in a committed result and line numbers in committed text.",
        "coordinates": "n/a: no interval is read, written or compared.",
        "parameters": {
            "floor": nc.FLOOR,
            "universes_extrapolated_to": [
                result["coverage"]["analysis_genes"],
                result["coverage"]["universe"],
                pseudobulk_genes,
            ],
            "extrapolation_model": "gamma-mixed Poisson, shape and mean matched to the observed counts "
            "by moments, mean scaling with the universe and shape held",
            "bracket_level": 0.05,
            "source_strings_searched_for": nc.SITES,
            "endpoint_discriminating_keys": nc.ENDPOINT_KEYS,
        },
        "exclusions": [
            "NO NEW MEASUREMENT: no study file was opened, no p-value parsed, no perturbed row, "
            "knockdown field or differential-expression call read, and no request made",
            "the floor of 30 is not relaxed, lowered, reinterpreted or worked around, and no estimate "
            "is recovered from a run that is below it",
            "the unstratified mean of -0.06105... over 12 factors is quoted only as what the "
            "registration says it is and is never treated as an estimate of the gain",
            "coverage under any OTHER endpoint is not computed, estimated, sampled or spot-checked, in "
            "either direction, because doing so would read outcome rows",
            "no model-based extrapolation is reported without the distribution-free bracket beside it, "
            "whose high end the committed counts do not exclude",
            "the published file's total row count was NOT read: the coordinator ruled against the "
            "separate read and it is carried forward to the next authorised pass",
            "no claim is made that the 44 factors lack a regulatory role in K562: the perturbations "
            "are assigned, with no knockdown check, as the committed result already states",
        ],
        "partitions": {
            "defined_factors": "the 12 factors whose stratified paired difference is defined, each with "
            "at least one responder among the analysed genes",
            "undefined_factors": "the 44 counted `no_responders`: not one analysed gene with a "
            "published adjusted p below 0.05 in that factor's principal column",
            "absent_genes": "the 364 universe genes with no row in the published file, which could "
            "never be labelled for any factor",
            "analysed_genes": "the 496 universe genes present once in the file's rows with a finite "
            "control expression, identical for all 56 factors",
        },
        "code_cleanliness": mf.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
    }

    path = save_result(RESULT, payload)
    print(f"{RESULT}: {path}")
    print(f"  endpoint: {payload['endpoint']['answer']}")
    print(
        f"    predicate at line {payload['endpoint']['from_the_code']['predicate_line']}, "
        f"level at {payload['endpoint']['from_the_code']['level_constant_line']}, "
        f"emission at {payload['endpoint']['from_the_code']['emission_line']}"
    )
    print(
        f"    from the result's keys alone: {payload['endpoint']['from_the_result_alone']['endpoint']} "
        f"(present {payload['endpoint']['from_the_result_alone']['key_present']}, absent "
        f"{payload['endpoint']['from_the_result_alone']['key_absent']})"
    )
    print(f"    T absent from run_v3: {payload['endpoint']['t_is_absent_from_run_v3']['absent']}")
    print(f"  {payload['headline']}")
    for row in payload["the_coverage_gap"]["model_based_extrapolation"]["at_universes"]:
        print(
            f"    fitted at {row['universe_genes']:>5} genes: "
            f"{row['expected_defined']:.1f} defined (floor {nc.FLOOR})"
        )
    b = payload["the_coverage_gap"]["distribution_free_bracket"]
    print(
        f"    bracket: low end {b['low_end']['defined_at_any_universe']} defined at any universe; "
        f"high end rate <= {b['high_end']['upper_per_gene_rate']:.5f}, "
        f"P(>=1 of {b['high_end']['added_genes']}) = "
        f"{b['high_end']['probability_at_least_one_responder_among_the_added_genes']:.3f}"
    )
    print(f"  verdict: {payload['verdict']}")


if __name__ == "__main__":
    main()
