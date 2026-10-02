# SPDX-License-Identifier: AGPL-3.0-or-later
"""Take the argmax-cell count registered before it (data/results/argmaxcell.json).

    uv run --frozen python scripts/argmaxcell_count.py --gate-sample chrY chr21   # measure RSS only
    uv run --frozen python scripts/argmaxcell_count.py                            # the run

The question: a compiled rule's cell is an argmax over the scorer's tracks, so does it carry
cell-type information at all? For element-gene pairs the ENCODE CRISPRi benchmark MEASURED in a
named cell, at what rate is the element's compiled cell label the measured cell, against the rate at
which that label appears over every compiled rule genome-wide?

Every threshold, population, denominator, base rate and reading is imported from
`genomeos.attribution.argmaxcell` and was committed in data/results/argmaxcell_registration.json
before any count of this lane existed. This script chooses nothing.

READ DISCIPLINE, which is why this is a committed script and not an inline loop. The cached-element
loader falls back to whole-chromosome archives and caches every archive it opens; a peer's loop over
93 elements reached 23 GB RSS tonight. Nothing here opens an archive. One chromosome's COMPACT
element table is decoded one element at a time, four scalars are kept per element, the text is
released before the next chromosome, and `argmaxcell.check_rss` RAISES at the registered ceiling
after every chromosome. `--gate-sample` measures the figure on the smallest and a middling
chromosome and writes no result, so the full run starts against a measured number.

This lane edits no rule, no threshold and no `.bio` file.
0 model requests, no network, no money.
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import argmaxcell as ac  # noqa: E402
from genomeos.attribution import crispri as cr  # noqa: E402

RESULT = "argmaxcell"
REGISTRATION = "data/results/argmaxcell_registration.json"

OWN_CODE = (
    "genomeos/attribution/argmaxcell.py",
    "scripts/argmaxcell_register.py",
    "scripts/argmaxcell_count.py",
    "tests/test_argmaxcell.py",
)

#: The committed summaries this run opens through `load_result`, declared so the traced reads
#: reconcile. `enhancer_targets_all_<chrom>` is the pointer at the compact table; the other two runs
#: carry their elements inline, and `targets.attributed` reads all three.
SUMMARIES = tuple(
    f"data/results/{run}_{chrom}.json"
    for chrom in ac.CHROMS
    for run in ("enhancer_targets_all", "constrained_targets", "enhancer_targets")
)


def pairs_by_chromosome() -> tuple[dict[str, list[Any]], dict[str, int], dict[str, int]]:
    """Every valid benchmark pair, grouped by chromosome, with the per-cell and per-table tallies.

    Both cached tables are read. A pair is valid by the benchmark's own `ValidConnection`, which
    `crispri.parse` applies; no outcome column is read here and none is carried.
    """
    by_chrom: dict[str, list[Any]] = defaultdict(list)
    per_cell: dict[str, int] = defaultdict(int)
    per_table: dict[str, int] = {}
    for name in (cr.TRAINING, cr.HELDOUT):
        rows = cr.load(name)
        per_table[name] = len(rows)
        for p in rows:
            by_chrom[p.chrom].append(p)
            per_cell[p.cell] += 1
    return by_chrom, dict(per_cell), per_table


def run(chroms: tuple[str, ...], arms: tuple[str, ...]) -> dict[str, Any]:
    """One streaming pass: the arms' elements, the second base rate and the peak RSS."""
    by_chrom, pairs_per_cell, per_table = pairs_by_chromosome()
    armed = {c: ac.Arm(cell=c) for c in arms}
    genome_labels: dict[str, int] = defaultdict(int)
    streamed = 0
    per_chrom: dict[str, dict[str, int]] = {}
    peak = ac.max_rss_bytes()
    for chrom in chroms:
        table = ac.attributed_table(chrom)
        starts = [t[0] for t in table]
        for _s, _e, _i, label in table:
            genome_labels[label] += 1
        streamed += len(table)
        here = {"attributed_elements": len(table), "pairs": 0}
        for p in by_chrom.get(chrom, []):
            if p.cell not in armed:
                continue
            here["pairs"] += 1
            armed[p.cell].add(chrom, ac.overlapping(table, starts, p.start, p.end))
        per_chrom[chrom] = here
        del table, starts
        peak = max(peak, ac.check_rss(chrom))
        print(
            f"  {chrom}: {here['attributed_elements']:>6,} attributed elements, "
            f"{here['pairs']:>5,} armed pairs, rss {peak / 1024**2:.0f} MB",
            flush=True,
        )
    return {
        "arms": armed,
        "genome_labels": dict(genome_labels),
        "attributed_elements_streamed": streamed,
        "per_chromosome": per_chrom,
        "pairs_per_cell": pairs_per_cell,
        "pairs_per_table": per_table,
        "max_rss_bytes": peak,
    }


def cross_arm(armed: dict[str, ac.Arm], cen: dict[str, Any]) -> dict[str, Any]:
    """`CROSS_ARM`: every arm's rate on every arm's cell, each column's base rate beside it."""
    labels = {c: ac.labels_for(c) for c in armed}
    return {
        "columns": [
            {
                "cell": c,
                "labels": list(labels[c]),
                "genome_wide_base_rate": ac.base_rate(c, cen)["rate"],
            }
            for c in armed
        ],
        "rows": [
            {
                "arm": a,
                "elements": arm.elements,
                "rate_of_each_column": {
                    c: (
                        round(sum(arm.label_counts.get(x, 0) for x in labels[c]) / arm.elements, 6)
                        if arm.elements
                        else None
                    )
                    for c in armed
                },
            }
            for a, arm in armed.items()
        ],
        "what_it_is": ac.CROSS_ARM,
    }


def second_base_rates(genome_labels: dict[str, int], arms: tuple[str, ...]) -> dict[str, Any]:
    """`SECOND_BASE_RATE`: this run's own share of each arm's label over every element it streamed."""
    total = sum(genome_labels.values())
    return {
        "attributed_elements_streamed": total,
        "distinct_labels_seen": len(genome_labels),
        "per_arm": {
            c: {
                "labels": list(ac.labels_for(c)),
                "elements_with_this_label": sum(genome_labels.get(x, 0) for x in ac.labels_for(c)),
                "rate": (
                    round(sum(genome_labels.get(x, 0) for x in ac.labels_for(c)) / total, 6)
                    if total
                    else None
                ),
            }
            for c in arms
        },
        "what_it_is": ac.SECOND_BASE_RATE,
    }


def verdict(armed: dict[str, ac.Arm], readings: dict[str, Any]) -> dict[str, Any]:
    """The plain answer, decided by the measured figures through the registered branches.

    Written after the counts; the branch each arm takes was fixed before them by `FALSIFIER`, and
    `argmax_carries_cell_type_information` is `reading`'s own field and not a sentence chosen here.
    """
    read = {a: r for a, r in readings.items() if r.get("rate_reported")}
    no_info = [a for a, r in read.items() if r["reading"] == "(2)"]
    anti = [a for a, r in read.items() if r["reading"] == "(1)"]
    detected = [a for a, r in read.items() if r["reading"] == "(3)"]
    cannot_tell = [a for a, r in read.items() if r["reading"] == "(4)"]
    unresolved = detected + cannot_tell
    floor = [a for a, r in readings.items() if not r.get("rate_reported")]
    usable = [a for a, r in read.items() if r.get("usable")]
    underpowered = [a for a, r in readings.items() if (r.get("power") or {}).get("underpowered")]
    if read and not unresolved:
        answer = (
            "NO. Every arm with a rate is an EQUIVALENCE result: its whole clustered interval on d "
            "lies inside the tolerance, so each has excluded a difference larger than the "
            "tolerance and a rule's `cell` is NOT evidence of where it acts. The confound pushes "
            "the other way - the benchmark tested elements active in the cell it was testing - so "
            "this reading is not an artefact of the selection"
        )
    elif unresolved and not no_info:
        answer = (
            "A DIFFERENCE IS DETECTED on every arm with a rate, and it does NOT establish that the "
            "argmax carries cell-type information. The benchmark chose elements active in the cell "
            "it was testing, which produces this sign on its own; this lane cannot separate the two "
            "and does not claim to (CONFOUND). What can be said is the rate itself, with its "
            "denominator, and whether it clears the usability threshold"
        )
    elif cannot_tell and not (no_info or detected or anti):
        answer = (
            "THE DATA CANNOT TELL on any arm with a rate. Every clustered interval crosses a "
            "tolerance bound, so no arm detects a difference larger than the tolerance and no arm "
            "excludes one. This is NOT a negative result and may not be reported as one"
        )
    elif read:
        answer = (
            f"SPLIT, and reported as split, with each arm named by the reading its own interval "
            f"selected. Reading (2), an equivalence result - the whole clustered interval on d "
            f"inside the tolerance, so a difference larger than the tolerance is EXCLUDED and the "
            f"compiled cell is not evidence of where the rule acts: {sorted(no_info)}. Reading "
            f"(3), a difference detected whose sign the benchmark's selection of tested elements "
            f"explains as readily as the model does, so it establishes nothing about cell-type "
            f"information: {sorted(detected)}. Reading (1), anti-correlated: {sorted(anti)}. "
            f"Reading (4), THE DATA CANNOT TELL - which is not a null and must never be read as "
            f"one: {sorted(cannot_tell)}. No set is pooled with another and no arm's reading is "
            "carried to another arm"
        )
    else:
        answer = (
            "NOT ANSWERED on these populations: every arm fell below the registered element floor, "
            "so counts are reported and no rate is. An absent rate is not a negative result"
        )
    return {
        "written_after_the_counts": (
            "this block was written after the counts and is a reading of them. The branch each arm "
            "takes was fixed BEFORE them, by FALSIFIER under AMENDMENT_1, and is applied by "
            "`argmaxcell.reading` from the CLUSTERED interval; nothing here chooses a sentence that "
            "a figure did not select"
        ),
        "does_the_argmax_carry_cell_type_information": answer,
        "arms_with_no_cell_type_information_reading_2": sorted(no_info),
        "reading_2_is_an_equivalence_result": (
            "an arm is in the list above only because its WHOLE clustered interval on d lies inside "
            "the tolerance, so it has EXCLUDED a difference larger than the tolerance. An arm that "
            "merely failed to show one is in `arms_the_data_cannot_tell_reading_4` and the two lists "
            "must never be merged"
        ),
        "arms_anti_correlated_reading_1": sorted(anti),
        "arms_where_a_difference_is_detected_but_not_attributable_reading_3": sorted(detected),
        "arms_the_data_cannot_tell_reading_4": sorted(cannot_tell),
        "arms_below_the_element_floor": sorted(floor),
        "arms_declared_underpowered": sorted(underpowered),
        "arms_whose_label_clears_the_usability_threshold": sorted(usable),
        "what_an_equivalence_result_excludes_on_a_small_base_rate": {
            "the_limit": (
                "stated because reading (2) is an ABSOLUTE equivalence against a tolerance of "
                f"{ac.TOLERANCE}, and an arm whose base rate is far below that can satisfy it "
                "while still carrying a large RELATIVE enrichment. On such an arm reading (2) "
                "means `no difference larger than the tolerance` and does NOT mean `the label "
                "carries nothing`. The multiple each arm's tolerance is of its own base rate is "
                "below, so a reader can see which arms the reading is strong on and which it is "
                "weak on without recomputing anything"
            ),
            "tolerance_as_a_multiple_of_the_base_rate": {
                a: (
                    round(ac.TOLERANCE / r["committed_base_rate"], 1)
                    if r.get("committed_base_rate")
                    else None
                )
                for a, r in read.items()
            },
            "observed_as_a_multiple_of_the_base_rate": {
                a: (
                    round(r["observed_rate"] / r["committed_base_rate"], 2)
                    if r.get("committed_base_rate")
                    else None
                )
                for a, r in read.items()
            },
        },
        "what_this_does_not_say": ac.VALIDATES_NOTHING,
        "cannot_establish": list(ac.CANNOT_ESTABLISH),
        "no_rule_is_changed_by_this": ac.NO_RECOMMENDATION,
    }


def inputs() -> list[dict[str, Any]]:
    """Declared so the traced reads reconcile: a declared DIRECTORY covers the files under it.

    The 24 compact tables are declared as their directory, one digest over all of them, because the
    run reads every one. The 72 committed run summaries are declared individually: they are the
    pointer files, and a pointer declared without the bytes it points at is the defect
    `manifest.undeclared_reads` exists to catch.
    """
    out = [
        mf.input_entry(
            ac.CENSUS,
            partition="the committed genome-wide census of compiled rule cell labels: the BASE RATE "
            "of every arm, read from its per_cell block and computed by no code of this lane",
        )
    ]
    for name in (cr.TRAINING, cr.HELDOUT):
        out.append(
            mf.input_entry(
                cr.KNOWLEDGE / name,
                partition="a cached ENCODE CRISPRi benchmark table, read for each valid pair's "
                "CellType, chrom, chromStart, chromEnd and gene only",
            )
        )
    out.append(
        mf.input_entry(
            ac.COMPACT,
            partition="the all-enhancer sweep's 24 compact per-chromosome element tables: the "
            "attributed elements and their predicted_coding tissue. Declared as the directory, one "
            "digest over every table the run streams. NOT the per-element archives, which are not "
            "opened",
        )
    )
    for rel in SUMMARIES:
        p = Path(rel)
        if p.exists():
            out.append(
                mf.input_entry(
                    p,
                    partition="a committed deletion-run summary read through load_result: the "
                    "pointer at a compact table, or a sampled run's inline elements",
                )
            )
    return out


def payload(state: dict[str, Any]) -> dict[str, Any]:
    cen = ac.census()
    armed: dict[str, ac.Arm] = state["arms"]
    bases = {c: ac.base_rate(c, cen) for c in armed}
    # AMENDMENT_2: the design effect is measured from this result's own estimable arms FIRST, then
    # applied to any arm whose bootstrap is degenerate. Two passes, so a degenerate arm's reading
    # never rests on an assumption about its own clustering.
    deff = ac.measured_design_effect(armed)
    readings = {c: ac.reading(arm, bases[c]["rate"], deff=deff["taken"]) for c, arm in armed.items()}
    entries = inputs()
    return {
        "result": RESULT,
        "date": date.today().isoformat(),
        "lane": "lane-argmaxcell",
        "registration": (
            f"{REGISTRATION}, committed in a commit of its own before any count of this lane existed"
        ),
        "question": ac.QUESTION,
        "argmax_rule": ac.ARGMAX_RULE,
        "argmax_is_a_selection": ac.ARGMAX_IS_A_SELECTION,
        "population": ac.POPULATION,
        "denominator": ac.DENOMINATOR,
        "denominator_is_separate": ac.DENOMINATOR_IS_SEPARATE,
        "falsifier": ac.FALSIFIER,
        "amendment_1": ac.AMENDMENT_1,
        "amendment_2": ac.AMENDMENT_2,
        "measured_design_effect": deff,
        "amendment_timeline": ac.AMENDMENT_TIMELINE,
        "reporting_code_changed_after_a_run": ac.REPORTING_CODE_CHANGED_AFTER_A_RUN,
        "relative_bands_from_now_on": ac.RELATIVE_BANDS_FROM_NOW_ON,
        "confound_registered_before_any_count": ac.CONFOUND,
        "base_rate_rule": ac.BASE_RATE_RULE,
        "base_rate_is_not_uniform": ac.BASE_RATE_IS_NOT_UNIFORM,
        "base_rates": bases,
        "per_arm": {
            c: {
                "cell": c,
                "labels": list(arm.labels),
                "pairs_in_cell": state["pairs_per_cell"].get(c, 0),
                "pairs_of_this_arm_seen": arm.pairs,
                "pairs_on_an_attributed_element": arm.pairs_on_an_element,
                "pairs_whose_element_carries_the_measured_cell": arm.pair_level_matching,
                "pair_level_rate_is_secondary": (
                    "reported because a reader will ask for it; it is NOT the headline, because the "
                    "compiled cell belongs to the element and a pair count weights an element by "
                    "how many genes the screen happened to assay near it (DENOMINATOR)"
                ),
                "pair_level_rate": (
                    round(arm.pair_level_matching / arm.pairs_on_an_element, 6)
                    if arm.pairs_on_an_element
                    else None
                ),
                "elements": arm.elements,
                "elements_whose_label_is_the_measured_cell": arm.elements_matching,
                "observed_rate": arm.rate,
                "reading": readings[c],
                "power": readings[c]["power"],
                "clusters": {
                    "kind": arm.cluster_choice()[0],
                    "count": len(arm.cluster_choice()[1]),
                    "minimum": ac.MIN_CLUSTERS,
                },
                "label_distribution_top10": ac.top_labels(arm, cen),
            }
            for c, arm in armed.items()
        },
        "cross_arm_matrix": cross_arm(armed, cen),
        "second_base_rate": second_base_rates(state["genome_labels"], tuple(armed)),
        "base_rate_agreement": {
            "committed_census_rules": cen["rules"],
            "attributed_elements_streamed": state["attributed_elements_streamed"],
            "difference": state["attributed_elements_streamed"] - cen["rules"],
            "what_a_difference_means": (
                "the committed census counted RULE lines emitted by compile_chromosome and this run "
                "counted the attributed ELEMENTS those lines are written from, one rule per "
                "element. A difference is reported as a difference and is not reconciled away: it "
                "would mean the two populations are not identical, and in that case the registered "
                "reading still uses the COMMITTED base rate, which is the one fixed in advance"
            ),
        },
        "per_chromosome": state["per_chromosome"],
        "pairs_per_table": state["pairs_per_table"],
        "pairs_per_cell_in_the_benchmark_tables": state["pairs_per_cell"],
        "verdict": verdict(armed, readings),
        "label_rule": ac.LABEL_RULE,
        "alias_rule": ac.ALIAS_RULE,
        "overlap_rule": ac.OVERLAP_RULE,
        "interval_is_binomial": ac.INTERVAL_IS_BINOMIAL,
        "grouping": ac.GROUPING,
        "name_match_is_not_a_confirmation": ac.NAME_MATCH_IS_NOT_A_CONFIRMATION,
        "metadata_copy_limitation": ac.METADATA_COPY_LIMITATION,
        "retention_is_not_this_question": ac.RETENTION_IS_NOT_THIS_QUESTION,
        "validates_nothing": ac.VALIDATES_NOTHING,
        "no_recommendation": ac.NO_RECOMMENDATION,
        "cannot_establish": list(ac.CANNOT_ESTABLISH),
        "read_discipline": ac.READ_DISCIPLINE,
        "rss": {
            "ceiling_bytes": ac.RSS_CEILING_BYTES,
            "peak_bytes": state["max_rss_bytes"],
            "peak_mib": round(state["max_rss_bytes"] / 1024**2, 1),
            "no_archive_was_opened": (
                "no per-element archive and no cached-element loader call. The peak above is the "
                "measured figure for streaming all 24 compact tables one at a time"
            ),
        },
        "alphagenome_requests": 0,
        "money": "none: no request was sent and no key was read",
        "result_manifest": {
            "sources": [
                {
                    "accession": "the committed genome-wide context-evidence census of compiled "
                    "rule cell labels",
                    "version": "the committed data/results/context_evidence.json on this machine",
                },
                {
                    "accession": cr.PUBLISHED_SOURCE,
                    "version": "the two benchmark tables cached under data/knowledge/crispri",
                },
                {
                    "accession": "the all-enhancer deletion sweep's compact per-chromosome element "
                    "tables and the 72 committed run summaries that name them",
                    "version": "the tables and summaries on this machine, digested in inputs",
                },
            ],
            "inputs": entries,
            "input_count": len(entries),
            "assembly": "GRCh38: the assembly both benchmark tables and the deletion sweep are on",
            "coordinates": {
                "base": 0,
                "interval": "half-open",
                "note": "the benchmark tables' own chromStart/chromEnd and the element tables' own "
                "start/end, compared only by OVERLAP_RULE",
            },
            "parameters": {
                "tolerance": ac.TOLERANCE,
                "usable_rate": ac.USABLE,
                "minimum_elements_for_a_rate": ac.MIN_ELEMENTS,
                "reach": cr.REACH,
                "arms": list(armed),
                "rss_ceiling_bytes": ac.RSS_CEILING_BYTES,
                "minimum_clusters_for_an_interval": ac.MIN_CLUSTERS,
                "bootstrap_draws": ac.BOOTSTRAPS,
                "bootstrap_seed": ac.SEED,
                "locus_span": ac.LOCUS_SPAN,
                "no_threshold_of_any_rule_is_set_or_moved": "MIN_EFFECT, the direction clauses and "
                "every magnitude floor stay exactly as committed",
            },
            "exclusions": [
                "every benchmark pair whose ValidConnection is not TRUE",
                "every benchmark pair whose CellType is not one of the arms",
                "every element with no predicted coding target: it emits no compiled rule",
                "every per-element archive and every cached-element loader call",
                "every outcome, DHS, H3K27ac and power column of the benchmark tables",
            ],
            "partitions": {
                "per_arm_elements": "the DISTINCT attributed elements overlapping that arm's pairs: "
                "the headline denominator, one per arm, never pooled across arms",
                "per_arm_pairs": "that arm's valid pairs: the secondary denominator",
                "the_census_440589": "the base rate's denominator, from the committed census",
            },
            "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="the argmax-cell count")
    ap.add_argument("--gate-sample", nargs="*", default=None, metavar="CHROM")
    ap.add_argument("--chroms", nargs="*", default=None, metavar="CHROM")
    a = ap.parse_args()
    if a.gate_sample is not None:
        chroms = tuple(a.gate_sample) or ("chrY", "chr21")
        print(f"gate sample over {len(chroms)} chromosomes, no result is written")
        state = run(chroms, ac.ARMS)
        ceiling = ac.RSS_CEILING_BYTES / 1024**2
        print(f"peak rss {state['max_rss_bytes'] / 1024**2:.0f} MiB of a ceiling of {ceiling:.0f} MiB")
        print(f"attributed elements streamed {state['attributed_elements_streamed']:,}")
        return 0
    chroms = tuple(a.chroms) if a.chroms else ac.CHROMS
    if chroms != ac.CHROMS:
        raise SystemExit(
            f"{RESULT} is the genome-wide count and its base rate is a genome-wide base rate; a run "
            f"over {len(chroms)} chromosomes would not be comparable with it"
        )
    state = run(chroms, ac.ARMS)
    from genomeos.results import save_result

    body = payload(state)
    path = save_result(RESULT, body)
    print(f"wrote: {path}")
    for cell, row in body["per_arm"].items():
        r = row["reading"]
        if r.get("rate_reported"):
            print(
                f"  {cell:18} {row['elements_whose_label_is_the_measured_cell']:>5} / "
                f"{row['elements']:<6} = {r['observed_rate']:.4f}  base {r['committed_base_rate']:.6f}"
                f"  d {r['difference_point_estimate']:+.4f}  dCI {r['clustered_ci95_on_the_difference']}"
                f"  {r['clustered_by']}x{r['clustered_interval_provenance']['clusters']}"
                f"  reading {r['reading']}  carries={r['argmax_carries_cell_type_information']}"
            )
        else:
            print(f"  {cell:18} {row['elements']:>5} elements: below the floor, no rate reported")
    print(f"  peak rss {body['rss']['peak_mib']} MiB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
