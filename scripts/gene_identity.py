#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Count the compiled rules whose target resolves to a different gene from the one that was measured.

Reads the 24 compiled programs, resolves each rule's target token and the gene its deletion answer
measured to an Ensembl gene id, and counts the three outcomes the registration fixed:
`data/results/gene_identity_registration.json` must be on disk, and every definition is imported from
`genomeos.attribution.gene_identity` rather than restated here.

    uv run --frozen python scripts/gene_identity.py                    # all 24 chromosomes
    uv run --frozen python scripts/gene_identity.py --chrom chr21 --result gene_identity_chr21

The run refuses to write a genome-wide result whose rule count is not the 440,589 the denominator was
fixed at, and a chr21 result whose count is not 5,176.
"""

from __future__ import annotations

import argparse
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import gene_identity as gi  # noqa: E402
from genomeos.results import load_result, save_result  # noqa: E402

REGISTRATION = "gene_identity_registration"
EXPECTED_GENOME_WIDE = 440_589
EXPECTED_CHR21 = 5_176
OWN_CODE = (
    "genomeos/attribution/gene_identity.py",
    "scripts/gene_identity_register.py",
    "scripts/gene_identity.py",
    "tests/test_gene_identity.py",
)
EXAMPLES = 12


def share(num: int, den: int) -> float | None:
    return None if den == 0 else round(num / den, 4)


def run(chroms: list[str]) -> dict[str, Any]:
    """Every rule of every named chromosome, classified at both registered half-windows."""
    outcomes: Counter[str] = Counter()
    second: Counter[str] = Counter()
    per_chrom: dict[str, Counter[str]] = {}
    by_activity: dict[str, Counter[str]] = defaultdict(Counter)
    by_effect: dict[str, Counter[str]] = defaultdict(Counter)
    by_distance: dict[str, Counter[str]] = defaultdict(Counter)
    by_cell: dict[str, Counter[str]] = defaultdict(Counter)
    distance_test: Counter[str] = Counter()
    distance_test_outcomes: Counter[str] = Counter()
    token_was_an_id = 0
    ambiguous_with_one_coding: Counter[str] = Counter()
    window_records: dict[str, int] = {}
    examples: list[dict[str, Any]] = []
    crispri_rules: list[dict[str, Any]] = []
    total = 0

    for chrom in chroms:
        rows = gi.annotation(chrom)
        if not rows:
            raise SystemExit(f"no GENCODE v50 annotation on disk for {chrom}: {gi.annotation_path(chrom)}")
        symbols, ids = gi.by_symbol(rows), gi.by_id(rows)
        rules = list(gi.rules(chrom))
        names: dict[str, frozenset[str]] = {}
        wanted: dict[str, set[str]] = defaultdict(set)
        for r in rules:
            if r.source == gi.SOURCE_PREDICTED:
                if r.target not in names:
                    names[r.target] = gi.recorded_names(r.target, symbols, ids)
                wanted[r.element] |= names[r.target]
        counts, seen, in_cache = gi.window_gene_counts(chrom, {k: frozenset(v) for k, v in wanted.items()})
        window_records[chrom] = seen
        here: Counter[str] = Counter()
        for r in rules:
            total += 1
            if r.source != gi.SOURCE_PREDICTED:
                occ: int | None = None
            elif r.element not in in_cache:
                occ = None
            else:
                got = counts.get(r.element, Counter())
                occ = sum(got.get(n, 0) for n in names[r.target])
            out = gi.classify(r, symbols, ids, occ)
            out2 = gi.classify(r, symbols, ids, occ, half_window=gi.SENSITIVITY_HALF_WINDOW)
            outcomes[out.outcome] += 1
            second[out2.outcome] += 1
            here[out.outcome] += 1
            if out.token_was_an_id:
                token_was_an_id += 1
            by_activity[r.activity or "none"][out.outcome] += 1
            by_effect[r.effect_band or "no band"][out.outcome] += 1
            by_cell[r.cell][out.outcome] += 1
            compiled = gi.compiled_locus(r.target, symbols, ids)
            d = gi.compiled_distance(r, compiled)
            by_distance[gi.distance_band(d)][out.outcome] += 1
            if d is not None and d > gi.HALF_WINDOW:
                distance_test["rules"] += 1
                distance_test_outcomes[out.outcome] += 1
                # Added after the chr21 run, and said so where it is reported: the chr21 run showed
                # this test selecting names with one locus on the chromosome, so how many loci the
                # name has is counted rather than left to be assumed from the distance.
                loci = len(symbols.get(r.target, ()))
                distance_test[
                    "name_has_one_locus_on_the_chromosome"
                    if loci == 1
                    else "name_has_more_than_one_locus_on_the_chromosome"
                    if loci > 1
                    else "name_has_no_locus_on_the_chromosome"
                ] += 1
                if (
                    compiled is not None
                    and compiled.end > r.midpoint - gi.HALF_WINDOW
                    and compiled.start < r.midpoint + gi.HALF_WINDOW
                ):
                    distance_test["compiled_gene_body_reaches_into_the_scorer_window"] += 1
            if out.outcome == "two_or_more_annotated_loci_of_that_name_in_the_scorer_window":
                cand = gi.window_candidates(r.target, r.midpoint, symbols)
                ambiguous_with_one_coding[
                    "exactly_one_is_protein_coding"
                    if sum(1 for g in cand if g.coding) == 1
                    else "not_exactly_one_is_protein_coding"
                ] += 1
            if out.outcome == gi.DIFFER and len(examples) < EXAMPLES:
                m = ids.get(out.measured_id or "")
                examples.append(
                    {
                        "element": r.element,
                        "locus": f"{r.chrom}:{r.start}-{r.end}",
                        "target_token": r.target,
                        "cell": r.cell,
                        "activity": r.activity,
                        "compiled_gene_id": out.compiled_id,
                        "compiled_tss_distance_bp": d,
                        "measured_gene_id": out.measured_id,
                        "measured_tss_distance_bp": gi.compiled_distance(r, m) if m else None,
                    }
                )
            if r.source == gi.SOURCE_MEASURED:
                crispri_rules.append(
                    {
                        "element": r.element,
                        "chrom": r.chrom,
                        "target_token": r.target,
                        "cell": r.cell,
                        "compiled_gene_id": out.compiled_id,
                    }
                )
        per_chrom[chrom] = here
        del rows, symbols, ids, rules, counts, wanted, in_cache

    return {
        "rules": total,
        "outcomes": dict(sorted(outcomes.items())),
        "outcomes_at_the_second_half_window": dict(sorted(second.items())),
        "per_chromosome": {c: dict(sorted(v.items())) for c, v in per_chrom.items()},
        "by_activity": {k: dict(sorted(v.items())) for k, v in sorted(by_activity.items())},
        "by_effect_band": {k: dict(sorted(v.items())) for k, v in sorted(by_effect.items())},
        "by_distance_band": {k: dict(sorted(v.items())) for k, v in sorted(by_distance.items())},
        "by_cell": {k: dict(sorted(v.items())) for k, v in sorted(by_cell.items())},
        "distance_test": {
            "rules": distance_test["rules"],
            "outcomes": dict(sorted(distance_test_outcomes.items())),
            "breakdown": {k: v for k, v in sorted(distance_test.items()) if k != "rules"},
        },
        "target_token_was_itself_an_ensembl_id": token_was_an_id,
        "ambiguous_window_loci_by_type": dict(sorted(ambiguous_with_one_coding.items())),
        "window_cache_records_per_chromosome": window_records,
        "examples_of_differ": examples,
        "crispri_rules": crispri_rules,
    }


def crispri_comparison(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """The experimental layer's own comparison: the screen records its Ensembl gene id, so compare it.

    Reported on its own and never added into the deletion-answer counts. If the benchmark table is
    not on disk this says so by name rather than reporting a smaller denominator.
    """
    out: dict[str, Any] = {
        "what_it_is": (
            "the rules of the experimental layer have no deletion answer. The CRISPRi screen does "
            "record the Ensembl gene id of the gene it measured (`measuredGeneEnsemblId`, carried as "
            "holdout.Unit.gene_id), so the compiled side is compared against that instead. This is a "
            "different test on a different population and is never added into the counts above"
        ),
        "rules": len(rows),
    }
    try:
        from genomeos.attribution import holdout as ho

        units = ho.all_units()
    except Exception as exc:  # the benchmark table is not on disk, or cannot be read
        out["not_computed"] = f"{type(exc).__name__}: {exc}"
        return out
    ids_for: dict[tuple[str, str], set[str]] = defaultdict(set)
    for name, us in units.items():
        if not name.startswith("crispri:"):
            continue
        for u in us:
            if u.gene and u.gene_id:
                ids_for[(u.gene, u.cell)].add(u.gene_id)
                ids_for[(u.gene, "")].add(u.gene_id)
    tally: Counter[str] = Counter()
    for r in rows:
        recorded = ids_for.get((r["target_token"], r["cell"])) or ids_for.get((r["target_token"], ""))
        if not recorded:
            tally["no_recorded_gene_id_for_that_name_in_the_screen"] += 1
        elif r["compiled_gene_id"] is None:
            tally["target_absent_from_the_chromosome_annotation"] += 1
        elif r["compiled_gene_id"] in recorded:
            tally[gi.AGREE] += 1
        else:
            tally[gi.DIFFER] += 1
    out["outcomes"] = dict(sorted(tally.items()))
    out["sums_to_the_rules"] = sum(tally.values()) == len(rows)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--chrom", action="append", help="one chromosome; repeat for several")
    ap.add_argument("--result", default="gene_identity")
    args = ap.parse_args()
    started = time.time()

    if load_result(REGISTRATION) is None:
        raise SystemExit(
            f"{REGISTRATION} is not on disk: the registration is written before the counts, "
            "never after. Run scripts/gene_identity_register.py first"
        )
    chroms = args.chrom or list(gi.CHROMS)
    data = run(chroms)

    expected = (
        EXPECTED_GENOME_WIDE
        if chroms == list(gi.CHROMS)
        else (EXPECTED_CHR21 if chroms == ["chr21"] else None)
    )
    if expected is not None and data["rules"] != expected:
        raise SystemExit(
            f"the population does not reproduce: {data['rules']} rules read, {expected} expected. "
            "Nothing is written"
        )
    classified = sum(data["outcomes"].values())
    if classified != data["rules"]:
        raise SystemExit(f"{classified} classified against {data['rules']} rules: nothing is written")

    outcomes = data["outcomes"]
    agree, differ = outcomes.get(gi.AGREE, 0), outcomes.get(gi.DIFFER, 0)
    unresolvable = {k: v for k, v in outcomes.items() if k in gi.CAUSES}
    resolvable = agree + differ

    def differ_share(table: dict[str, dict[str, int]]) -> dict[str, Any]:
        out = {}
        for level, counts in table.items():
            a, d = counts.get(gi.AGREE, 0), counts.get(gi.DIFFER, 0)
            out[level] = {
                "differ": d,
                "agree": a,
                "resolvable": a + d,
                "differ_share_of_the_resolvable": share(d, a + d),
                "unresolvable_in_this_level": sum(v for k, v in counts.items() if k in gi.CAUSES),
                "population_of_that_share": f"the {a + d} resolvable rules of this level: agree + differ",
            }
        return out

    crispri = crispri_comparison(data.pop("crispri_rules"))

    payload: dict[str, Any] = {
        "result": args.result,
        "date": "2026-10-02",
        "lane": "lane-identity",
        "question": (
            "for every compiled rule, does the target token resolve to the same Ensembl gene id as "
            "the gene the deletion answer measured?"
        ),
        "status": (
            "descriptive. The definitions, the precedence and both half-windows were fixed in "
            f"data/results/{REGISTRATION}.json before any of these counts was read, and none of them "
            "moved afterwards"
        ),
        "registration": f"data/results/{REGISTRATION}.json",
        "chromosomes": chroms,
        "population": {
            "rules": data["rules"],
            "source": "every `rule` line of data/knowledge/compiled/noncoding_chr*.bio",
            "reproduces_the_fixed_denominator": expected is not None and data["rules"] == expected,
            "fixed_denominator": expected,
        },
        "counts": {
            "differ": differ,
            "agree": agree,
            "unresolvable": sum(unresolvable.values()),
            "unresolvable_by_cause": dict(sorted(unresolvable.items())),
            "sums_to_the_population": agree + differ + sum(unresolvable.values()) == data["rules"],
            "differ_share_of_the_resolvable": share(differ, resolvable),
            "population_of_that_share": f"the {resolvable} resolvable rules: agree + differ",
            "what_differ_means": gi.MIS_RESOLUTION,
            "what_differ_does_not_mean": (
                "it is not a share of the model's predictions that are wrong, and it says nothing "
                "about any one rule"
            ),
        },
        "at_the_second_half_window": {
            "half_window": gi.SENSITIVITY_HALF_WINDOW,
            "outcomes": data["outcomes_at_the_second_half_window"],
            "differ": data["outcomes_at_the_second_half_window"].get(gi.DIFFER, 0),
            "agree": data["outcomes_at_the_second_half_window"].get(gi.AGREE, 0),
            "moved_from_the_registered_half_window": {
                k: data["outcomes_at_the_second_half_window"].get(k, 0) - outcomes.get(k, 0)
                for k in sorted(set(outcomes) | set(data["outcomes_at_the_second_half_window"]))
            },
        },
        "the_1678": {
            "imported_reading": gi.imported_readings()["lane_notopen_1678"],
            "imported_from": gi.imported_readings()["lane_notopen_1678_where"],
            "its_population": gi.imported_readings()["lane_notopen_1678_population"],
            "the_same_distance_test_over_this_whole_population": data["distance_test"]["rules"],
            "its_outcomes_here": data["distance_test"]["outcomes"],
            "breakdown_of_those_rules": data["distance_test"]["breakdown"],
            "the_breakdown_was_added_after_the_chr21_run": (
                "`name_has_one_locus_on_the_chromosome` and the three figures beside it were added "
                "after the chr21 run, not fixed in the registration, because the chr21 run showed "
                "the distance test selecting names with a single locus. They are descriptive counts "
                "over a population the registration already fixed and they moved no definition, no "
                "threshold and no precedence; the amendment records them"
            ),
            "the_two_populations": (
                "the 1,678 is the distance test restricted to the 55,084 rules in state "
                "not_open_in_reader. This lane's distance-test figure is the same test over all "
                f"{data['rules']} rules and does not re-derive the openness state, so the 1,678 is a "
                "subset of it and is not equated with any count here"
            ),
        },
        "correlations": {
            "by_activity": differ_share(data["by_activity"]),
            "by_effect_band": differ_share(data["by_effect_band"]),
            "by_distance_band": differ_share(data["by_distance_band"]),
            "by_cell": differ_share(data["by_cell"]),
            "reported_descriptively": (
                "these shares are reported as they came out. No cause is inferred from a difference "
                "between levels"
            ),
            "the_activity_figures_this_is_read_beside": gi.imported_readings()["lane_notopen_activity_axis"],
        },
        "per_chromosome": data["per_chromosome"],
        "target_token_was_itself_an_ensembl_id": data["target_token_was_itself_an_ensembl_id"],
        "ambiguous_window_loci_by_type": {
            "counts": data["ambiguous_window_loci_by_type"],
            "what_it_is": (
                "of the rules left unresolvable because two or more annotated loci of that name lie "
                "in the scorer window, how many have exactly one protein_coding locus among them. "
                "`predicted_coding` was chosen from the symbols that are protein_coding somewhere on "
                "the chromosome, which is a property of the SYMBOL and not of the locus, so this is "
                "reported as a count and is NOT used to resolve any of them"
            ),
        },
        "window_cache_records_per_chromosome": data["window_cache_records_per_chromosome"],
        "examples_of_differ": data["examples_of_differ"],
        "the_experimental_layer": crispri,
        "what_this_cannot_establish": gi.limitations(),
        "no_interval_is_reported": (
            "these are exhaustive counts over a fixed rule set, not an estimate from a sample, so no "
            "confidence interval is computed and none would mean anything here"
        ),
        "nothing_was_changed": (
            "no compiled rule was deleted, rewritten or relabelled by this run. It reads the compiled "
            "programs and writes a count beside them"
        ),
        "alphagenome_requests": 0,
        "network_requests": 0,
        "money": "none: every input was already on disk",
        "seconds": round(time.time() - started, 1),
    }

    inputs = [mf.input_entry(gi.annotation_path(c)) for c in chroms if gi.annotation_path(c).exists()]
    inputs += [mf.input_entry(gi.program_path(c)) for c in chroms if gi.program_path(c).exists()]
    inputs += [
        mf.input_entry(gi.WINDOW_CACHE / f"{c}.json.gz")
        for c in chroms
        if (gi.WINDOW_CACHE / f"{c}.json.gz").exists()
    ]
    inputs.append(mf.input_entry(Path(f"data/results/{REGISTRATION}.json")))
    inputs.sort(key=lambda e: e["path"])

    payload["result_manifest"] = {
        "sources": [
            {"accession": "GENCODE", "version": "v50, data/reference/gencode_v50_chr*.gff3.gz"},
            {
                "accession": "AlphaGenome deletion answers per element, as the project cached them",
                "version": (
                    "data/knowledge/alphagenome/elements/chr*.json.gz, the per-element record of "
                    "every gene the scorer read in the element's window"
                ),
            },
            {
                "accession": "the compiled non-coding programs of the 24 chromosomes",
                "version": "data/knowledge/compiled/noncoding_chr*.bio, read as they are on disk",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "mis_resolution": gi.MIS_RESOLUTION,
            "compiled_side": gi.COMPILED_RESOLUTION,
            "measured_side": gi.MEASURED_RESOLUTION,
            "half_window": gi.HALF_WINDOW,
            "second_half_window": gi.SENSITIVITY_HALF_WINDOW,
            "effect_band_call": gi.EFFECT_BAND_CALL,
            "causes": gi.CAUSES,
            "chromosomes": chroms,
        },
        "exclusions": [
            "no rule is excluded: agree, differ and the named causes sum to the population",
            "no new cut-off is introduced: the half-window, the effect band and the distance bands "
            "are imported from where the project already fixed them",
            "the experimental layer's own comparison is reported separately and is never added into "
            "the deletion-answer counts",
        ],
        "partitions": {
            "by_activity": "the element block's own `activity:` value as the compiler wrote it",
            "by_effect_band": "the rule line's own `strength:`, banded at STRONG_EFFECT",
            "by_distance_band": "executor._band over the compiled side's element-midpoint-to-TSS distance",
            "by_cell": "the rule's own `when: cell_type` label, unchanged",
        },
        "code_cleanliness": mf.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
    }

    path = save_result(args.result, payload)
    print(f"wrote {path}")
    print(f"rules {data['rules']}  differ {differ}  agree {agree}  unresolvable {sum(unresolvable.values())}")
    for k, v in sorted(unresolvable.items(), key=lambda kv: -kv[1]):
        print(f"  {k}: {v}")
    print(f"distance test over this population: {data['distance_test']['rules']}")


if __name__ == "__main__":
    main()
