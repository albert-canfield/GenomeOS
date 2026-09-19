# SPDX-License-Identifier: AGPL-3.0-or-later
"""Are the project's own confidences calibrated? The measured layer as a reliability table.

    uv run python scripts/confidence_calibration.py [--chroms chr21,chr22] [--no-save] [--markdown]

GenomeOS states an evidence kind and a confidence on every fact and the Evidence explorer counts
them; nobody has asked whether a fact stated at 0.4 is right 40% of the time. `measured.py` made the
question askable on 2026-09-17 by attaching CRISPRi, lentiMPRA, VISTA and saturation-mutagenesis
outcomes to compiled elements. This bins those elements by the confidence the compiler wrote on them
and reports the observed agreement rate per band with a 95% Wilson interval - per assay, never
pooled, and with both of the layer's agreement denominators.

The failure patterns are named in `genomeos/attribution/confidence_calibration.py` and so are the
bands, which are imported unchanged from `target_calibration` where a different lane fixed them for a
different quantity. Nothing in this script chooses a cut.

`--markdown` prints the documentation section from the saved result, so no figure in prose is typed
beside a computed one.
"""

from __future__ import annotations

import argparse
import time
from typing import Any

from genomeos.attribution import confidence_calibration as cc
from genomeos.attribution import measured
from genomeos.attribution.targets import attributed
from genomeos.results import load_result, save_result

CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]
RESULT = "confidence_calibration_genome"


def read_chromosome(chrom: str) -> tuple[list[cc.Element], dict[str, Any]]:
    """One chromosome's compiled elements, with the measured layer and the covariates attached."""
    from genomeos.attribution.unknown_scoring import annotate

    elements = attributed(chrom)
    layer = measured.Layer.load(chrom)
    rows = measured.rows(chrom, elements, layer)
    cov = {r["id"]: r for r in annotate(chrom, elements)}
    got = cc.elements_of(chrom, elements, layer, rows, cov)
    per = {
        # the three denominators in order, on every chromosome as in the pooled result
        "compiled_elements_with_a_stated_confidence": len(got),
        "elements_in_any_assay_footprint": sum(1 for e in got if e.eligible),
        "elements_measured": sum(1 for e in got if e.verdicts),
        "elements_without_covariates": sum(1 for e in got if not e.covariates_present),
        "stated_confidence_median": (
            round(sorted(e.stated_confidence for e in got)[len(got) // 2], 3) if got else None
        ),
        "measured_by_assay": {a: sum(1 for e in got if a in e.verdicts) for a in measured.ASSAYS},
    }
    return got, per


def markdown(payload: dict[str, Any]) -> str:
    """The documentation section, every figure read from the result rather than typed beside it."""
    slice_ = payload["the_measured_slice"]
    im = slice_["imbalance"]
    ans = payload["the_answer"]
    out: list[str] = []
    a = out.append
    a(
        "## Are the project's own confidences calibrated? The stated level sits outside its own "
        f"interval in {ans['bands_where_it_is_not']} of "
        f"{ans['bands_populated_across_every_table']} bands, and the ordering holds in "
        f"{len(ans['tables_where_the_rate_rises_with_the_confidence'])} of the {ans['tables']} "
        f"reliability tables ({payload['date']})"
    )
    a("")
    a(
        "Every compiled fact carries a confidence and the Evidence explorer counts them. Nobody had "
        "asked whether a fact stated at 0.4 is right about 40% of the time. The measured layer made "
        "the question askable, and the answer is reported per assay because *agreement* is a "
        "different event in each of the four."
    )
    a("")
    a(f"**What the number is.** `{payload['stated_confidence_formula']}`")
    a("")
    b = payload["bands"]
    a(
        "**What a failure would look like, written down before the rates were computed.** The four "
        "patterns are named constants in `genomeos/attribution/confidence_calibration.py` and the "
        "bands are imported unchanged from `target_calibration`, where another lane fixed them on the "
        "same day for a different quantity. A band holds a verdict only at "
        f"{b['min_for_a_band']} measured elements or more, an assay gets a verdict only at "
        f"{b['min_bands_for_a_verdict']} judged bands or more, and the level clears the bar at "
        f"{b['min_share_consistent']:.0%} of judged bands consistent."
    )
    a("")
    for name, text in payload["pre_stated_patterns"].items():
        a(f"- **`{name}`** — {text}.")
    a("")
    a(
        "**The answer, per assay and never pooled.** The rate rises with the stated confidence in "
        + (", ".join(f"`{k}`" for k in ans["tables_where_the_rate_rises_with_the_confidence"]) or "none")
        + "; it does not in "
        + (", ".join(f"`{k}`" for k in ans["tables_where_it_does_not"]) or "none")
        + "; the pre-registered rule refuses a verdict in "
        + (", ".join(f"`{k}`" for k in ans["tables_refused"]) or "none")
        + "."
    )
    a("")
    a("| table | n | mean stated | observed | median offset | ECE | pattern |")
    a("|---|---|---|---|---|---|---|")
    for k, r in ans["per_denominator"].items():
        if "table" in r:
            continue
        a(
            f"| `{k}` | {r['elements']:,} | {r['mean_stated_confidence_overall']} | "
            f"{r['observed_overall']} | **{r['median_gap_over_populated_bands']:+}** | "
            f"{r['expected_calibration_error']} | `{r['pattern']}` |"
        )
    a("")
    a(
        "**The offsets have opposite signs, and that is the result of record.** "
        + ", ".join(f"`{k}` {v:+}" for k, v in ans["offsets_by_table"].items())
        + ". "
        + ans["level_is_a_property_of_the_population"][0].upper()
        + ans["level_is_a_property_of_the_population"][1:]
        + "."
    )
    a("")
    a(
        "**Three denominators, in order.** Of "
        f"**{slice_['compiled_elements_with_a_stated_confidence']:,}** compiled elements that state a "
        f"confidence, **{slice_['measured_by_any_assay']:,}** are measured by any assay "
        f"(**{slice_['coverage'] * 100:.2f}%**) and **{slice_['never_measured']:,}** are not. "
        f"{slice_['elements_without_covariates']:,} have no covariates and enter no comparison."
    )
    a("")
    a("**The measured slice against the rest, on the covariates and on the confidence itself:**")
    a("")
    a("| | measured | never measured | ratio |")
    a("|---|---|---|---|")
    for key, label in (
        ("length", "length (bp)"),
        ("gc", "GC"),
        ("nearest_coding_tss", "distance to a coding TSS (bp)"),
        ("stated_confidence", "**stated confidence**"),
    ):
        r = im[key]
        a(f"| {label} | {r['target_median']} | {r['control_median']} | {r['ratio']} |")
    a("")
    a("**Coverage by band — the selection this lane cannot undo, as a number:**")
    a("")
    a("| band | measured | never measured | coverage |")
    a("|---|---|---|---|")
    for b in cc.BAND_NAMES:
        m, u = slice_["band_counts_measured"][b], slice_["band_counts_never_measured"][b]
        if not (m or u):
            continue
        cov = slice_["coverage_by_band"][b]
        a(f"| {b} | {m:,} | {u:,} | {cov * 100:.2f}% |")
    a("")
    for assay in measured.ASSAYS:
        rep = payload["per_assay"][assay]
        a(f"### {assay}")
        a("")
        a(f"*{rep['what_agreement_means']['question']}* — {rep['what_agreement_means']['why']}.")
        a("")
        if "refused" in rep:
            a(
                f"**No table.** {rep['elements_in_the_footprint']} elements in the footprint, "
                f"{rep['elements_measured']} measured, {rep['agrees']} agreeing. {rep['refused']}."
            )
            a("")
            continue
        a(
            f"{rep['elements_in_the_footprint']:,} compiled elements in this assay's footprint and "
            f"{rep['elements_measured']:,} measured; "
            + (
                f"{rep['elements_where_the_predicted_gene_was_tested']:,} of them are tested by it, "
                "which is every one, so the two denominators coincide - this assay never asks about "
                "the predicted gene at all, and the wide rate is kept as a named key so the narrow "
                "one is never the only rate on the page."
                if rep["denominators_coincide"]
                else f"{rep['elements_where_the_predicted_gene_was_tested']:,} of them are ones where "
                "the screen tested the very gene the deletion named."
            )
        )
        a("")
        for key in ("where_the_predicted_gene_was_tested", "over_all_matched_elements"):
            block = rep[key]
            if "identical_to" in block:
                continue
            t, v = block["table"], block["verdict"]
            a(
                f"**{key}** — n = {t['elements']:,}, observed {t['observed_overall']}, mean stated "
                f"{t['mean_stated_confidence_overall']}, expected calibration error "
                f"{t['expected_calibration_error']}."
            )
            a("")
            a(
                "| band | compiled | in footprint | measured | coverage | mean stated | agrees | "
                "observed | 95% CI | inside | length | GC | TSS |"
            )
            a("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
            for r in t["rows"]:
                if not r["measured_by_this_assay"]:
                    continue
                ci = f"{r['ci95'][0]:.3f}–{r['ci95'][1]:.3f}"
                mark = "yes" if r["stated_inside_the_interval"] else "no"
                if not r["judged"]:
                    mark += " *(not judged)*"
                a(
                    f"| {r['band']} | {r['compiled_elements_in_the_band']:,} | "
                    f"{r['in_this_assays_footprint']:,} | {r['measured_by_this_assay']:,} | "
                    f"{r['coverage_of_the_band'] * 100:.2f}% | {r['mean_stated_confidence']} | "
                    f"{r['agrees']:,} | **{r['observed']}** | {ci} | {mark} | "
                    f"{r['median_length']} | {r['median_gc']} | {r['median_nearest_coding_tss']} |"
                )
            a("")
            if v["pattern"] == "refused":
                short = v.get("bands_populated_but_below_the_bar") or {}
                a(
                    f"**Verdict: refused.** {v['reading']}."
                    + (
                        " Short of the bar: "
                        + "; ".join(
                            f"band {b} holds {d['measured']}, {d['short_of_the_bar_by']} below "
                            f"{cc.MIN_FOR_A_BAND}"
                            for b, d in short.items()
                        )
                        + "."
                        if short
                        else ""
                    )
                    + " The gaps are still described: "
                    + ", ".join(
                        f"{b} {g:+}" for b, g in v["described_not_judged"]["gap_per_populated_band"].items()
                    )
                    + f", median {v['described_not_judged']['median_gap_over_populated_bands']:+}."
                )
            else:
                a(
                    f"**Verdict: `{v['pattern']}`.** Ordering: {v['ordering']} — "
                    f"{v['lowest_judged_band']} at {v['lowest_observed']} to "
                    f"{v['highest_judged_band']} at {v['highest_observed']}, intervals "
                    f"{'disjoint' if v['intervals_disjoint'] else 'overlapping'}. Level: the stated "
                    f"confidence lies inside its band's interval in "
                    f"{v['bands_where_the_stated_confidence_is_inside_the_interval']} of "
                    f"{v['bands_judged']} judged bands "
                    f"({v['share_consistent'] * 100:.0f}%, bar {v['min_share_consistent'] * 100:.0f}%); "
                    f"median signed offset {v['median_signed_offset_observed_minus_stated']:+}, range "
                    f"{v['offset_range'][0]:+} to {v['offset_range'][1]:+}."
                )
            a("")
            ci_ = block["carries_information"]
            pres = ci_["input_presence"]
            a(
                "Is the band reading the confidence or the covariate? Bought or free first: "
                + "; ".join(
                    f"{name[0].upper() + name[1:]} is **{d['kind']}** "
                    f"({d['targets_with_the_input']:,}/{d['targets']:,} of the measured arm, "
                    f"{d['controls_with_the_input']:,}/{d['controls']:,} of the unmeasured)"
                    for name, d in pres.items()
                )
                + "."
            )
            if "comparison_refused" in ci_:
                a("")
                a(f"*Comparison refused:* {ci_['comparison_refused']}")
            else:
                s = ci_["standardised"]
                a("")
                a(
                    f"Standardised on length, GC and distance to a coding TSS: "
                    f"{s['targets']:,} high-confidence elements against {s['controls']:,} low, "
                    f"{s['targets_matched']:,} matched and "
                    f"{s['dropped_for_want_of_a_control']:,} dropped for want of a control; raw "
                    f"{s['raw']['difference']:+}, matched **{s['matched']['difference']:+}** at "
                    f"one-sided p {s['matched']['p_one_sided']}."
                )
            a("")
    pooled = payload["why_pooling_is_refused"]
    a("**Why there is no pooled number.** " + pooled[0].upper() + pooled[1:] + ".")
    a("")
    a("**What this cannot say.** " + payload["scope"] + ".")
    a("")
    a("**What would falsify the transfer.** " + payload["falsifies_transfer"] + ".")
    a("")
    a(
        f"**Cost.** {len(payload['per_chromosome'])} chromosomes read, binned and compared in "
        f"**{payload['seconds']} seconds**, no requests and no model. "
        "`uv run python scripts/confidence_calibration.py`; `--markdown` prints this section back out "
        "of `data/results/confidence_calibration_genome.json`, which is how every figure above got "
        "here rather than being typed beside a computed one."
    )
    a("")
    return "\n".join(out)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--chroms", default="", help="comma separated; default every chromosome")
    ap.add_argument("--no-save", action="store_true")
    ap.add_argument("--markdown", action="store_true", help="print the doc section from the saved result")
    args = ap.parse_args(argv)

    if args.markdown:
        got = load_result(RESULT)
        if got is None:
            print(f"no {RESULT}.json; run this script without --markdown first")
            return 2
        print(markdown(got))
        return 0

    chroms = [c.strip() for c in args.chroms.split(",") if c.strip()] or CHROMS
    t0 = time.time()
    population: list[cc.Element] = []
    per_chromosome: dict[str, Any] = {}
    refused: dict[str, str] = {}
    for chrom in chroms:
        try:
            got, per = read_chromosome(chrom)
        except (FileNotFoundError, KeyError, ValueError) as exc:
            refused[chrom] = f"{type(exc).__name__}: {exc}"[:160]
            print(f"{chrom}: not read ({refused[chrom]})", flush=True)
            continue
        population.extend(got)
        per_chromosome[chrom] = per
        print(
            f"{chrom}: {per['compiled_elements_with_a_stated_confidence']:,} compiled elements state a "
            f"confidence, {per['elements_in_any_assay_footprint']:,} are in an assay's footprint, "
            f"{per['elements_measured']:,} are measured",
            flush=True,
        )

    payload = cc.report(population, per_chromosome)
    payload["chromosomes_not_read"] = refused
    payload["seconds"] = round(time.time() - t0, 1)
    print()
    for assay in measured.ASSAYS:
        print(cc.summary_line(payload, assay), flush=True)

    if not args.no_save:
        print(f"\nwrote {save_result(RESULT, payload)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
