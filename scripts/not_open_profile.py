# SPDX-License-Identifier: AGPL-3.0-or-later
"""Type and locate the `not_open_in_reader` rules: four breakdowns and the sharp cross-tabulation.

    uv run --frozen python scripts/not_open_profile.py
    uv run --frozen python scripts/not_open_profile.py --chroms chr21 --result not_open_profile_chr21
    uv run --frozen python scripts/not_open_profile.py --dump /tmp/sharp.tsv

55,084 of the 81,635 assessable compiled rules are `not_open_in_reader`. This run takes the same rules,
in the same order, and breaks that number down by the cell the rule asserts in, the class the compiler
already writes on the element, the effect band the AlphaGenome prediction already carries, and the
element-to-gene distance in the bands the project already uses. It then crosses the three axes that
make a candidate sharp: a strong predicted effect, a cell whose assignment the base-rate comparison
found informative, and the element not detected open there.

Nothing is refitted, no threshold is moved, no verdict is moved and nothing is deleted. The state per
rule is `genomeos.attribution.context_evidence.state_for` unchanged, and this run's per-biosample open
and not-open counts are checked against the census result's own, so the two cannot disagree about the
population being described. No model is called and nothing is downloaded.

What `not_open_in_reader` does and does not mean is in `genomeos.attribution.not_open_profile` and is
copied into the result, so the result cannot be read without it.
"""

from __future__ import annotations

import argparse
import gzip
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import context_evidence as ce  # noqa: E402
from genomeos.attribution import holdout as ho  # noqa: E402
from genomeos.attribution import not_open_profile as prof  # noqa: E402
from genomeos.genome import reader  # noqa: E402
from genomeos.results import RESULTS_DIR, load_result, save_result  # noqa: E402

RESULT = "not_open_profile"
CENSUS = "context_evidence"
CHROMS = tuple(f"chr{c}" for c in [*range(1, 23), "X", "Y"])
#: The census figures this run must reproduce, from data/results/context_evidence.json.
RULES_ON_THE_RECORD = 440_589
NOT_OPEN_ON_THE_RECORD = 55_084
OPEN_ON_THE_RECORD = 26_551
ASSESSABLE_ON_THE_RECORD = 81_635
#: How many of the sharpest rules the result names per cell, strongest predicted effect first.
EXAMPLES_PER_CELL = 3
EXAMPLES_OVERALL = 12

#: This lane's own files. Everything else uncommitted in this shared checkout belongs to another lane.
OWN_CODE = (
    "genomeos/attribution/not_open_profile.py",
    "scripts/not_open_profile.py",
    "tests/test_not_open_profile.py",
)

#: The result files the compiler reads per chromosome, as globs: the same list the census used, so the
#: two runs record the same inputs. Each matched file is recorded by its own path below.
INPUT_GLOBS = (
    "budget_chr*.json",
    "budget_axes_chr*.json",
    "variation_chr*.json",
    "duplication_chr*.json",
    "domains_chr*.json",
    "unknown_chr*.json",
    "reader_*_chr*.json",
    "ccres_chr*.bed.gz",
    "enhancer_targets_chr*.json",
    "enhancer_targets_all_chr*.json",
    "constrained_targets_chr*.json",
)

TSV_HEADER = (
    "chromosome",
    "element",
    "locus",
    "gene",
    "cell",
    "reader_biosample",
    "element_class",
    "activity_axis",
    "origin_axis",
    "alphagenome_log2_fold_change",
    "effect_band",
    "distance_to_tss_bp",
    "distance_band",
)


def share(part: int, whole: int) -> float | None:
    return None if not whole else round(part / whole, 4)


def _split(counts: dict[str, Counter[str]], key: str) -> dict[str, Any]:
    """One breakdown: per value, the not-open rules, the assessable rules and the share, with the
    denominator named beside every share because a share without one is not a figure."""
    out: dict[str, Any] = {}
    for value, c in sorted(counts.items(), key=lambda kv: (-kv[1][ce.STATE_NOT_OPEN], kv[0])):
        assessable = c[ce.STATE_OPEN] + c[ce.STATE_NOT_OPEN]
        out[value] = {
            "not_open_in_reader": c[ce.STATE_NOT_OPEN],
            "open_in_reader": c[ce.STATE_OPEN],
            "assessable_rules": assessable,
            "not_open_share_of_this_" + key: share(c[ce.STATE_NOT_OPEN], assessable),
            "population_of_that_share": (
                f"the {assessable} assessable rules of this {key}: open_in_reader + not_open_in_reader"
            ),
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", default=",".join(CHROMS))
    ap.add_argument("--result", default=RESULT, help="a run over fewer chromosomes writes its own name")
    ap.add_argument(
        "--dump",
        default=None,
        help=(
            "write the whole sharp cross-tabulation as a TSV (gzipped when the name ends .gz). The "
            "result itself names a handful; this writes every row for working through by hand"
        ),
    )
    args = ap.parse_args()
    chroms = [c for c in args.chroms.split(",") if c]
    if args.result == RESULT and tuple(chroms) != CHROMS:
        raise SystemExit(f"{RESULT} is the genome-wide profile; a run over {len(chroms)} needs --result")
    whole_genome = tuple(chroms) == CHROMS

    started = time.time()
    info = prof.informativeness()
    table = ce.mapping()
    readers = ce.Readers(RESULTS_DIR)

    per_state: Counter[str] = Counter()
    by_biosample: dict[str, Counter[str]] = defaultdict(Counter)
    by_label: dict[str, Counter[str]] = defaultdict(Counter)
    by_class: dict[str, Counter[str]] = defaultdict(Counter)
    by_activity: dict[str, Counter[str]] = defaultdict(Counter)
    by_origin: dict[str, Counter[str]] = defaultdict(Counter)
    by_effect: dict[str, Counter[str]] = defaultdict(Counter)
    by_distance: dict[str, Counter[str]] = defaultdict(Counter)
    by_source: dict[str, Counter[str]] = defaultdict(Counter)
    by_chrom: dict[str, Counter[str]] = defaultdict(Counter)
    beyond_window: Counter[str] = Counter()
    sharp: list[tuple[prof.Rule, str]] = []
    total = 0

    for chrom in chroms:
        here = Counter()
        for r in prof.rules(chrom, RESULTS_DIR):
            total += 1
            state, _rest = ce.parse_value(ce.state_for(r.cell, chrom, r.start, r.end, readers, table))
            per_state[state] += 1
            here[state] += 1
            by_chrom[chrom][state] += 1
            if state == ce.STATE_NOT_ASSESSABLE:
                continue
            biosample, _term, _why = ce.reader_for(r.cell, table)
            assert biosample is not None  # an assessable rule has a reader biosample by definition
            by_biosample[biosample][state] += 1
            by_label[r.cell][state] += 1
            by_class[r.element_class][state] += 1
            by_activity[r.activity_axis][state] += 1
            by_origin[r.origin_axis][state] += 1
            by_effect[r.effect_band or prof.SOURCE_MEASURED][state] += 1
            by_distance[r.distance_band][state] += 1
            by_source[r.source][state] += 1
            if state != ce.STATE_NOT_OPEN:
                continue
            if r.beyond_the_scorer_window:
                beyond_window[biosample] += 1
            if r.effect_band == prof.STRONG and info.informative(biosample):
                sharp.append((r, biosample))
        print(
            f"  {chrom}: {sum(here.values()):,} rules, {here[ce.STATE_NOT_OPEN]:,} not_open_in_reader",
            flush=True,
        )

    assessable = per_state[ce.STATE_OPEN] + per_state[ce.STATE_NOT_OPEN]
    census = load_result(CENSUS, RESULTS_DIR) or {}
    agrees: Any
    if whole_genome:
        if (total, per_state[ce.STATE_NOT_OPEN], per_state[ce.STATE_OPEN]) != (
            RULES_ON_THE_RECORD,
            NOT_OPEN_ON_THE_RECORD,
            OPEN_ON_THE_RECORD,
        ):
            raise SystemExit(
                f"this run reads {total} rules, {per_state[ce.STATE_NOT_OPEN]} not_open_in_reader and "
                f"{per_state[ce.STATE_OPEN]} open_in_reader, against the census's "
                f"{RULES_ON_THE_RECORD}, {NOT_OPEN_ON_THE_RECORD} and {OPEN_ON_THE_RECORD}. The "
                "population being described is not the one on the record; nothing is written."
            )
        mine = {b: {s: c[s] for s in (ce.STATE_NOT_OPEN, ce.STATE_OPEN)} for b, c in by_biosample.items()}
        theirs = {
            b: {s: v.get(s, 0) for s in (ce.STATE_NOT_OPEN, ce.STATE_OPEN)}
            for b, v in (census.get("per_reader_biosample") or {}).items()
        }
        if theirs and mine != theirs:
            raise SystemExit(f"per-biosample counts differ from the census: {mine} against {theirs}")
        agrees = True
    else:
        agrees = f"n/a: the record is over all 24 chromosomes and this run covers {len(chroms)}"

    sharp.sort(key=lambda rb: (-abs(rb[0].effect or 0.0), rb[0].chrom, rb[0].start, rb[0].gene))
    strict = [rb for rb in sharp if rb[1] not in prof.CELLS_LANE_CONTEXT2_SINGLED_OUT]
    per_cell_sharp: dict[str, list[tuple[prof.Rule, str]]] = defaultdict(list)
    for rb in sharp:
        per_cell_sharp[rb[1]].append(rb)

    def rows(pairs: list[tuple[prof.Rule, str]]) -> list[dict[str, Any]]:
        return [{**r.row(), "reader_biosample": b} for r, b in pairs]

    payload: dict[str, Any] = {
        "question": (
            "of the 55,084 compiled rules whose element is not detected open in the cell they assert "
            "in, which cells, which element classes, which predicted effect sizes and which distances "
            "to the gene carry them - and how many are a strong predicted effect in a cell whose "
            "assignment the base-rate comparison found informative?"
        ),
        "status": (
            "descriptive and diagnostic, and NOT registered: context_evidence_registration.json fixed "
            "the states, the openness call, the mapping table and the denominator, and this changes "
            "none of them. Every state is context_evidence.state_for's, unchanged."
        ),
        "describes": "data/results/context_evidence.json",
        "chromosomes": chroms,
        "rules": total,
        "per_state": {s: per_state[s] for s in ce.STATES},
        "assessable_rules": assessable,
        "assessable_rules_on_the_record": ASSESSABLE_ON_THE_RECORD if whole_genome else None,
        "not_open_in_reader": per_state[ce.STATE_NOT_OPEN],
        "population": (
            f"the {per_state[ce.STATE_NOT_OPEN]} rules in state not_open_in_reader, out of the "
            f"{assessable} assessable rules (open_in_reader + not_open_in_reader), out of the {total} "
            "rules the compiler emits"
        ),
        "reproduces_the_census": agrees,
        "no_interval_is_reported": (
            "these are exhaustive counts over a fixed rule set, not an estimate from a sample, so no "
            "confidence interval is computed and none would mean anything here"
        ),
        "by_cell": {
            "reader_biosample": {
                b: {
                    **v,
                    "assignment_difference_over_the_base_rate": info.of(b),
                    "clears_the_registered_tolerance": info.informative(b),
                    "singled_out_by_lane_context2": b in prof.CELLS_LANE_CONTEXT2_SINGLED_OUT,
                }
                for b, v in _split(by_biosample, "biosample").items()
            },
            "compiled_cell_label": _split(by_label, "label"),
            "tolerance": info.tolerance,
            "informativeness_source": info.source,
            "reading": (
                "a not_open_in_reader rule in a cell whose assignment carries no information about "
                "openness means something different from one in a cell where it does: in the first, "
                "the cell the compiler chose is not telling us where the element is open, so the "
                "element not being open there is close to uninformative. The three cells "
                "lane-context2 singled out are flagged on every row"
            ),
        },
        "by_element_class": {
            **_split(by_class, "class"),
            "call": prof.registration()["element_class_call"],
        },
        "by_activity_axis": _split(by_activity, "activity"),
        "by_origin_axis": _split(by_origin, "origin"),
        "by_effect_band": {
            **_split(by_effect, "band"),
            "call": prof.EFFECT_BAND_CALL,
            "measured_rules_carry_no_alphagenome_effect": prof.NO_ALPHAGENOME_EFFECT,
        },
        "by_distance_band": {
            **_split(by_distance, "band"),
            "call": prof.DISTANCE_CALL,
            "beyond_the_scorer_window": {
                "rules": sum(beyond_window.values()),
                "per_biosample": dict(sorted(beyond_window.items())),
                "what_it_is": prof.BEYOND_THE_WINDOW,
                "population": (
                    f"the {per_state[ce.STATE_NOT_OPEN]} not_open_in_reader rules; these are counted "
                    "in every other breakdown unchanged"
                ),
            },
        },
        "by_source_of_the_rule": _split(by_source, "source"),
        "by_chromosome": {c: {s: by_chrom[c][s] for s in ce.STATES} for c in chroms},
        "the_sharp_cross_tabulation": {
            "what_it_crosses": (
                "not_open_in_reader AND a strong predicted effect AND a cell whose assignment clears "
                "the base-rate tolerance. These are the sharpest candidate defects the three axes "
                "can name together; sharper is not the same as wrong"
            ),
            "rules": len(sharp),
            "share_of_the_not_open_rules": share(len(sharp), per_state[ce.STATE_NOT_OPEN]),
            "population_of_that_share": f"the {per_state[ce.STATE_NOT_OPEN]} not_open_in_reader rules",
            "also_dropping_the_three_cells_lane_context2_singled_out": {
                "rules": len(strict),
                "cells_dropped": list(prof.CELLS_LANE_CONTEXT2_SINGLED_OUT),
                "why_both_are_given": (
                    "docs/ATTRIBUTION.md states SK-N-SH as the cell whose assignment carries no "
                    "information, and astrocyte and ovary as clearing the tolerance only just. The "
                    "first count applies the registered tolerance as arithmetic; the second also "
                    "drops the two that only just clear it, so the figure does not depend on reading "
                    "that sentence one way"
                ),
            },
            "per_cell": {
                b: {
                    "rules": len(v),
                    "assignment_difference_over_the_base_rate": info.of(b),
                    "share_of_this_cell_s_not_open_rules": share(len(v), by_biosample[b][ce.STATE_NOT_OPEN]),
                    "population_of_that_share": (
                        f"the {by_biosample[b][ce.STATE_NOT_OPEN]} not_open_in_reader rules assigned to {b}"
                    ),
                }
                for b, v in sorted(per_cell_sharp.items(), key=lambda kv: (-len(kv[1]), kv[0]))
            },
            "per_element_class": dict(Counter(r.element_class for r, _ in sharp).most_common()),
            "per_activity_axis": dict(Counter(r.activity_axis for r, _ in sharp).most_common()),
            "per_origin_axis": dict(Counter(r.origin_axis for r, _ in sharp).most_common()),
            "per_distance_band": dict(Counter(r.distance_band for r, _ in sharp).most_common()),
            "per_chromosome": dict(Counter(r.chrom for r, _ in sharp).most_common()),
            "effect_size": {
                "largest_absolute_log2_fold_change": (
                    round(abs(sharp[0][0].effect or 0.0), 4) if sharp else None
                ),
                "smallest_absolute_log2_fold_change": (
                    round(abs(sharp[-1][0].effect or 0.0), 4) if sharp else None
                ),
                "note": (
                    f"every rule here is in the strong band, so the smallest is at or just above "
                    f"STRONG_EFFECT = {prof.STRONG_EFFECT}"
                ),
            },
            "strongest_rules": rows(sharp[:EXAMPLES_OVERALL]),
            "strongest_rules_per_cell": {
                b: rows(v[:EXAMPLES_PER_CELL])
                for b, v in sorted(per_cell_sharp.items(), key=lambda kv: (-len(kv[1]), kv[0]))
            },
            "the_whole_list": (
                "the rows above are a handful for looking at by hand. The whole cross-tabulation is "
                "written as a TSV by this script's --dump option, which is reproducible from the "
                "committed code in one run and is not carried in this result"
            ),
        },
        **prof.registration(),
        "reading": {
            "moves": (
                "one number becomes a list: the 55,084 not_open_in_reader rules are located by cell, "
                "element class, predicted effect size and distance to the gene, and the three axes "
                "are crossed to name the sharpest candidates"
            ),
            "does_not_move": (
                "no verdict moved, no rule was deleted, no compiled label changed, no threshold or "
                "registered definition changed, and nothing was refitted. No rule is established to "
                "be wrong by anything here"
            ),
            "never": (
                "not_open_in_reader is never reported as the element being closed or the rule being "
                "contradicted, and open_in_reader is never reported as validation"
            ),
        },
        "code_cleanliness": ce.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
        "money": "none: every input was already on disk",
        "seconds": round(time.time() - started, 1),
    }

    # Every input by its own path: `manifest.files_entry` records a group under a label and
    # `scripts/manifest_rebuild.py` resolves an input by its `path`, so a grouped input is reported
    # absent without any file being checked. Recording each file by itself keeps the rebuild able to
    # compare what this run read.
    paths = {ce.TRACK_METADATA, Path(f"data/results/{CENSUS}.json"), Path(info.source)}
    paths.update(ho.gencode_paths())
    for cell in ce.READER_TERMS:
        for chrom in chroms:
            p = RESULTS_DIR / reader.peaks_path(cell, chrom).name
            if p.exists():
                paths.add(p)
    for glob in INPUT_GLOBS:
        for chrom in chroms:
            paths.update(RESULTS_DIR.glob(glob.replace("chr*", chrom)))
    inputs = [mf.input_entry(p, partition=None) for p in sorted(paths, key=lambda q: q.as_posix())]

    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "ENCODE DNase-seq narrowPeak, GRCh38, released, 13 biosamples",
                "version": "as reader v1 cached them under data/results/dnase_*_chr*.bed.gz",
            },
            {"accession": "AlphaGenome output track metadata", "version": str(ce.TRACK_METADATA)},
            {"accession": "GENCODE", "version": "v50, data/reference/gencode_v50_chr*.gff3.gz, for the TSS"},
            {
                "accession": "the compiled non-coding programs of the 24 chromosomes",
                "version": "re-enumerated in this run from the results on disk, not read from a file",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "openness_call": ce.OPENNESS_CALL,
            "measured_span": ce.SPAN_CALL,
            "effect_band_call": prof.EFFECT_BAND_CALL,
            "distance_call": prof.DISTANCE_CALL,
            "informative_call": prof.INFORMATIVE_CALL,
            "tolerance": info.tolerance,
            "chromosomes": chroms,
        },
        "exclusions": [
            "no rule is excluded from the breakdowns: the denominator is every assessable rule",
            "the sharp cross-tabulation is a selection, not an exclusion: every rule outside it is "
            "still counted in all four breakdowns",
            "no new cut-off is introduced; the openness call, the effect band, the distance bands and "
            "the informativeness tolerance are all imported from where the project already fixed them",
        ],
        "partitions": {
            "by_cell": "the reader biosample the rule's label maps to by ontology term, and the label",
            "by_element_class": (
                "`class:` as the compiler writes it on the element, with the compiler's own "
                "`activity:` and `origin:` axis values beside it"
            ),
            "by_effect_band": "`strength` as the cached AlphaGenome prediction carries it",
            "by_distance_band": "executor._band over the element-midpoint-to-TSS distance",
        },
        "code_cleanliness": ce.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
    }

    path = save_result(args.result, payload)
    if args.dump:
        out = Path(args.dump)
        out.parent.mkdir(parents=True, exist_ok=True)
        opener = gzip.open if out.name.endswith(".gz") else open
        with opener(out, "wt") as fh:
            fh.write("\t".join(TSV_HEADER) + "\n")
            for r, b in sharp:
                row = {**r.row(), "reader_biosample": b}
                fh.write("\t".join("" if row[k] is None else str(row[k]) for k in TSV_HEADER) + "\n")
        print(f"  the whole sharp cross-tabulation: {out} ({len(sharp):,} rows)")

    print(f"{args.result}: {path}")
    print(
        f"  rules {total:,}, assessable {assessable:,}, not_open_in_reader {per_state[ce.STATE_NOT_OPEN]:,}"
    )
    print(f"  sharp (strong x informative cell x not open): {len(sharp):,}; stricter {len(strict):,}")


if __name__ == "__main__":
    main()
