# SPDX-License-Identifier: AGPL-3.0-or-later
"""The count data/results/increase_registration.json fixed, taken after it was committed
(data/results/increase_population.json).

    uv run --frozen python scripts/increase_population.py
    uv run --frozen python scripts/increase_population.py --chroms chr21 --result increase_population_chr21

In the registered order:

1. **The ladder**, blind to direction at every position step: how many cached CRISPRi pairs are in the
   benchmark, how many are significant, how many are on an attributed element under the committed
   reciprocal-overlap rule, how many have a compiled rule on that element naming their own gene, and
   how many have one gated on their own cell. Every step reports the exhaustive outcome breakdown at
   its own denominator, and the last two steps also report the repression-axis column, which is the
   other side of lane-repress2's 11 and 0.
2. **P1 and P2 against both imported floors**, separately and never pooled.
3. **The verdict**, in the words `increases.GO` and `increases.NO_GO` registered before this run, with
   the data-or-code attribution of the step a short population fell at.

Nothing is changed. `measured.rule_links`, `measured.Layer.near`, `measured.reciprocal_overlap` and
`not_open_profile.rules` are called exactly as committed; P1 is counted beside the extractor, never
through a modified one. No threshold is moved, no verdict moves, no rule is deleted, nothing is
refitted, nothing is downloaded and no model request is made.
"""

from __future__ import annotations

import argparse
import sys
import time
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import (
    cell2,  # noqa: E402
    crispri,  # noqa: E402
)
from genomeos.attribution import increases as inc  # noqa: E402
from genomeos.attribution import measured as ms  # noqa: E402
from genomeos.attribution import not_open_profile as nop  # noqa: E402
from genomeos.attribution import repress2 as rp  # noqa: E402
from genomeos.attribution.targets import attributed  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

RESULT = "increase_population"
REGISTRATION = "increase_registration"
CHROMS = tuple(f"chr{c}" for c in [*range(1, 23), "X", "Y"])

#: This lane's own files. Everything else uncommitted in this shared checkout belongs to another lane.
OWN_CODE = (
    "genomeos/attribution/increases.py",
    "scripts/increase_register.py",
    "scripts/increase_population.py",
    "tests/test_increases.py",
)

#: The same per-chromosome results `not_open_profile.rules` reads, each digested by its own path.
INPUT_GLOBS = (
    "budget_chr*.json",
    "budget_axes_chr*.json",
    "variation_chr*.json",
    "duplication_chr*.json",
    "domains_chr*.json",
    "unknown_chr*.json",
    "ccres_chr*.bed.gz",
    "enhancer_targets_chr*.json",
    "enhancer_targets_all_chr*.json",
    "constrained_targets_chr*.json",
)

S1, S2, S3, S4, S5 = inc.LADDER_STEPS
#: The two axis columns reported at the last two steps, so each reads against lane-repress2's own step.
AXIS_COLUMNS = (rp.REPRESSES, rp.ACTIVATES)


class Tally:
    """One ladder step's denominator: the count, and the exhaustive outcome breakdown at it."""

    def __init__(self) -> None:
        self.n = 0
        self.outcomes: Counter[str] = Counter()

    def add(self, outcome: str) -> None:
        self.n += 1
        self.outcomes[outcome] += 1

    def payload(self) -> dict[str, Any]:
        return {"pairs": self.n, "outcome_breakdown": inc.breakdown(self.outcomes)}


class Index:
    """What a pair has to be found in to clear the last two ladder steps.

    `on_gene` holds every (attributed element, gene) a compiled rule targets. `axes_of_gene` maps each
    of those to the activity axes of every rule there, whatever cell it is gated on; `axes_of_cell`
    maps (element, gene, cell) to the axes of the rules gated on that cell. The element is the
    attributed one with the `_measured` suffix removed, so a predicted rule and a measured rule on the
    same element index under one key.
    """

    def __init__(self, rules: list[nop.Rule]) -> None:
        self.on_gene: set[tuple[str, str]] = set()
        self.axes_of_gene: dict[tuple[str, str], set[str]] = defaultdict(set)
        self.axes_of_cell: dict[tuple[str, str, str], set[str]] = defaultdict(set)
        for r in rules:
            e = inc.base_element_id(r.element)
            self.on_gene.add((e, r.gene))
            self.axes_of_gene[(e, r.gene)].add(r.activity_axis)
            self.axes_of_cell[(e, r.gene, r.cell)].add(r.activity_axis)


def overlap_index(
    elements: list[dict[str, Any]], layer: ms.Layer
) -> tuple[dict[int, list[str]], dict[str, tuple[int, int]]]:
    """Which attributed elements each cached pair is a measurement *of*, by the committed rule.

    `measured.rows` keeps a pair against an element when the two reciprocally overlap at
    `measured.RECIPROCAL_OVERLAP` or more, over the candidates `Layer.near` returns. The same two calls
    are made here from the element side, and the answer is indexed by the pair's position in
    `layer.crispri` so that a pair covered by several elements is still one pair.
    """
    where: dict[int, list[str]] = defaultdict(list)
    at = {id(p): i for i, p in enumerate(layer.crispri)}
    spans: dict[str, tuple[int, int]] = {}
    for e in elements:
        start, end = e["start"], e["end"]
        spans[e["id"]] = (start, end)
        for p in layer.near("crispri", start, end):
            if ms.reciprocal_overlap(start, end, p.start, p.end) >= ms.RECIPROCAL_OVERLAP:
                where[at[id(p)]].append(e["id"])
    return where, spans


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", default=",".join(CHROMS))
    ap.add_argument("--result", default=RESULT, help="a run over fewer chromosomes writes its own name")
    args = ap.parse_args()
    chroms = [c for c in args.chroms.split(",") if c]
    if args.result == RESULT and tuple(chroms) != CHROMS:
        raise SystemExit(f"{RESULT} is the genome-wide reading; a run over {len(chroms)} needs --result")

    started = time.time()
    steps: dict[str, Tally] = {s: Tally() for s in inc.LADDER_STEPS}
    axis_column: dict[str, dict[str, Tally]] = {s: {axis: Tally() for axis in AXIS_COLUMNS} for s in (S4, S5)}
    invalid = 0
    rules_seen = 0
    elements_seen = 0
    where: dict[str, Counter[str]] = {k: Counter() for k in ("cell", "chrom", "split", "dataset")}
    on_element: dict[str, Counter[str]] = {k: Counter() for k in ("cell", "chrom", "split")}
    p1_keys: list[cell2.LocusKey] = []
    p2_keys: list[cell2.LocusKey] = []
    p2_contested = 0
    contested: list[dict[str, Any]] = []
    p1_splits: Counter[str] = Counter()
    p1_cells: Counter[str] = Counter()
    p1_links: list[dict[str, Any]] = []

    for chrom in chroms:
        layer = ms.Layer.load(chrom)
        invalid += layer.crispri_invalid
        pairs = layer.crispri
        if not pairs:
            continue
        elements = attributed(chrom, RESULTS_DIR)
        elements_seen += len(elements)
        rules = nop.rules(chrom, RESULTS_DIR, layer)
        rules_seen += len(rules)
        index = Index(rules)
        at_element, spans = overlap_index(elements, layer)

        for i, p in enumerate(pairs):
            outcome = p.outcome
            up = outcome == inc.INCREASE
            steps[S1].add(outcome)
            if up:
                where["cell"][p.cell] += 1
                where["chrom"][chrom] += 1
                where["split"][p.split] += 1
                where["dataset"][p.dataset] += 1
            if not rp.significant({"outcome": outcome}):
                continue
            steps[S2].add(outcome)
            ids = at_element.get(i) or []
            if not ids:
                continue
            steps[S3].add(outcome)
            if up:
                on_element["cell"][p.cell] += 1
                on_element["chrom"][chrom] += 1
                on_element["split"][p.split] += 1
            named = [e for e in ids if (e, p.gene) in index.on_gene]
            if not named:
                continue
            steps[S4].add(outcome)
            seen = {ax for e in named for ax in index.axes_of_gene[(e, p.gene)]}
            for axis, tally in axis_column[S4].items():
                if axis in seen:
                    tally.add(outcome)
            gated = [e for e in named if (e, p.gene, p.cell) in index.axes_of_cell]
            if not gated:
                continue
            steps[S5].add(outcome)
            seen = {ax for e in gated for ax in index.axes_of_cell[(e, p.gene, p.cell)]}
            for axis, tally in axis_column[S5].items():
                if axis in seen:
                    tally.add(outcome)

        # ---- P1 and P2: one link per (element, gene, cell), counted beside the extractor ----------
        grouped: dict[tuple[str, str, str], list[ms.CrispriPair]] = defaultdict(list)
        for i, p in enumerate(pairs):
            for e in at_element.get(i) or []:
                grouped[(e, p.gene, p.cell)].append(p)
        for (e, gene, cell), gp in sorted(grouped.items()):
            ups = [p for p in gp if p.outcome == inc.INCREASE]
            if not ups:
                continue
            start, end = spans[e]
            key = cell2.LocusKey(cell, chrom, start, end, gene)
            down = sum(1 for p in gp if p.regulated)
            if rp.REPRESSES in index.axes_of_cell.get((e, gene, cell), ()):
                p2_keys.append(key)
                p2_contested += 1 if down else 0
            if down:
                contested.append(
                    {
                        "element": e,
                        "chrom": chrom,
                        "gene": gene,
                        "cell": cell,
                        "significant_increases": len(ups),
                        "regulated_pairs": down,
                    }
                )
                continue
            p1_keys.append(key)
            p1_cells[cell] += 1
            for p in ups:
                p1_splits[p.split] += 1
            if len(p1_links) < 40:
                p1_links.append(
                    {
                        "element": e,
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "gene": gene,
                        "cell": cell,
                        "increases": len(ups),
                        "largest_effect_size": round(max(p.effect_size for p in ups), 4),
                        "splits": sorted({p.split for p in ups}),
                    }
                )

    ladder = {
        "call": inc.LADDER_CALL,
        "overlap_rule": inc.OVERLAP_CALL,
        "rules": inc.RULES_CALL,
        "breakdown_call": inc.BREAKDOWN_CALL,
        "mirrors_lane_repress2": inc.LADDER_MIRRORS,
        "steps": {s: steps[s].payload() for s in inc.LADDER_STEPS},
        "reconciles": [
            {
                "step": s,
                "pairs": steps[s].n,
                "step_above": prev,
                "pairs_above": steps[prev].n,
                "at_or_below_the_step_above": steps[s].n <= steps[prev].n,
                "lost_here": steps[prev].n - steps[s].n,
                "increases_here": steps[s].outcomes.get(inc.INCREASE, 0),
                "increases_above": steps[prev].outcomes.get(inc.INCREASE, 0),
                "increases_at_or_below_the_step_above": (
                    steps[s].outcomes.get(inc.INCREASE, 0) <= steps[prev].outcomes.get(inc.INCREASE, 0)
                ),
            }
            for prev, s in zip(inc.LADDER_STEPS[:-1], inc.LADDER_STEPS[1:], strict=True)
        ],
        "activity_axis_column": {
            step: {axis: tally.payload() for axis, tally in cols.items()}
            for step, cols in axis_column.items()
        },
        "activity_axis_column_call": (
            "of the pairs that cleared this step, how many had at least one matching compiled rule "
            f"whose activity axis is exactly {rp.REPRESSES!r}, and how many exactly {rp.ACTIVATES!r}. "
            "A pair can appear in both columns when two rules match it, so the columns do not sum to "
            f"the step. The {rp.REPRESSES!r} column of these two steps is the other side of "
            "lane-repress2's `and_a_pair_on_the_rules_own_gene` (11) and `and_in_the_rules_own_cell` (0)"
        ),
        "benchmark_rows_the_benchmark_itself_marks_invalid": invalid,
        "invalid_rows_call": (
            "ValidConnection FALSE: promoter or exon overlaps the benchmark says are not a test of an "
            "enhancer-to-gene link. measured.parse_crispri returns them as a count and they are in no "
            "denominator above; they are not measured negatives"
        ),
    }

    gates = {
        inc.P1: inc.gate(p1_keys, inc.P1, inc.P1_CALL),
        inc.P2: inc.gate(p2_keys, inc.P2, inc.P2_CALL),
    }
    cleared = [name for name, g in gates.items() if g["meets_both_floors"]]
    killed = _killed_at(steps, gates)
    verdict = {
        "go_or_no_go": "go" if cleared else "no-go",
        "reading": inc.GO if cleared else inc.NO_GO,
        "populations_at_or_above_both_floors": cleared,
        "floors_quoted_from_the_code_that_defines_them": {
            "links": inc.POSITIVE_FLOOR,
            "links_defined_in": "genomeos/attribution/fresh.py: POSITIVE_FLOOR = 30",
            "independent_loci": inc.LOCUS_FLOOR,
            "independent_loci_defined_in": (
                "genomeos/attribution/cell2.py: POOLED_LOCUS_FLOOR = 20, imported into "
                "genomeos/attribution/fresh.py as LOCUS_FLOOR"
            ),
        },
        "what_follows": inc.WHAT_FOLLOWS["go" if cleared else "no_go"],
        "the_step_that_killed_the_population": killed,
        "data_or_code": inc.ATTRIBUTION_CALL,
    }

    payload: dict[str, Any] = {
        "result": args.result,
        "date": date.today().isoformat(),
        "lane": "lane-increase",
        "registered_before_this_run": f"data/results/{REGISTRATION}.json",
        "follows_from": inc.INHERITED_FROM,
        "chromosomes": chroms,
        "attributed_elements_enumerated": elements_seen,
        "rules_enumerated": rules_seen,
        "ladder": ladder,
        "where_the_increases_are": {
            "call": (
                "every significant increase of the benchmark, by the field named. "
                "`on_an_attributed_element` is the same count restricted to the third ladder step"
            ),
            "by_cell": dict(where["cell"].most_common()),
            "by_chromosome": dict(sorted(where["chrom"].items())),
            "by_split": dict(where["split"].most_common()),
            "by_dataset": dict(where["dataset"].most_common()),
            "on_an_attributed_element": {
                "by_cell": dict(on_element["cell"].most_common()),
                "by_chromosome": dict(sorted(on_element["chrom"].items())),
                "by_split": dict(on_element["split"].most_common()),
            },
        },
        "populations": {
            "order": list(inc.POPULATIONS),
            "never_pooled": inc.NEVER_POOLED,
            "per_population": gates,
            inc.P1: {
                "contest_rule": inc.P1_CONTEST_RULE,
                "split_rule": inc.P1_SPLIT_RULE,
                "increases_behind_the_links_by_split": dict(p1_splits.most_common()),
                "links_by_cell": dict(p1_cells.most_common()),
                "contested_by_a_regulated_pair": len(contested),
                "contested_rows": contested[:40],
                "contested_rows_shown": min(len(contested), 40),
                "first_links": p1_links,
                "first_links_shown": len(p1_links),
            },
            inc.P2: {
                "reconciles_with": inc.P2_RECONCILES_WITH,
                "lane_repress2_and_in_the_rules_own_cell": inc.INHERITED_LADDER["and_in_the_rules_own_cell"],
                "at_or_below_it": len(set(p2_keys)) <= inc.INHERITED_LADDER["and_in_the_rules_own_cell"],
                "of_those_links_contested_by_a_regulated_pair": p2_contested,
            },
        },
        "verdict": verdict,
        "cannot_establish": inc.CANNOT_ESTABLISH,
        "no_extractor_changed": inc.NO_EXTRACTOR_CHANGED,
        "the_cause_in_code": inc.CAUSE_IN_CODE,
        # nested, never spread: the registration carries its own `ladder` and `populations` keys, and
        # spreading it over this payload would replace the counts above with their descriptions
        "registration_as_committed": inc.registration(),
        "reading": {
            "moves": (
                "the significant increases of the CRISPRi benchmark are counted along the ladder "
                "registered before the count, and the two populations a changed extractor could reach "
                "are gated against both imported floors"
            ),
            "does_not_move": (
                "no extractor changed, no threshold or registered definition changed, no verdict "
                "moved, no compiled label changed, no rule was deleted and nothing was refitted. No "
                "repression call is established to be right or wrong by anything here"
            ),
        },
        "requests": 0,
        "money": "none: every input was already on disk",
        "seconds": round(time.time() - started, 1),
    }

    paths: set[Path] = {Path(f"data/results/{REGISTRATION}.json")}
    for name in ms.CRISPRI_FILES:
        p = crispri.KNOWLEDGE / name
        if p.exists():
            paths.add(p)
    for glob in INPUT_GLOBS:
        for chrom in chroms:
            paths.update(RESULTS_DIR.glob(glob.replace("chr*", chrom)))
    inputs = [mf.input_entry(p, partition=None) for p in sorted(paths, key=lambda q: q.as_posix())]

    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "ENCODE CRISPRi enhancer-gene benchmark "
                "(EngreitzLab/CRISPR_comparison, Gschwind et al.)",
                "version": "the two benchmark tables cached under data/knowledge, digested in inputs",
            },
            {
                "accession": "the compiled non-coding programs of the chromosomes this run covers",
                "version": "re-enumerated in this run from the results on disk, not read from a file",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "eligibility_predicate": inc.ELIGIBILITY_CALL,
            "independent_locus_rule": inc.LOCUS_RULE,
            "independent_locus_span": inc.LOCUS_SPAN,
            "locus_floor": inc.LOCUS_FLOOR,
            "link_floor": inc.POSITIVE_FLOOR,
            "reciprocal_overlap": ms.RECIPROCAL_OVERLAP,
            "reach": ms.REACH,
            "ladder_steps": list(inc.LADDER_STEPS),
            "populations": list(inc.POPULATIONS),
            "chromosomes": chroms,
            "aggregate_function": "count",
            "fill_value": 0,
        },
        "exclusions": [
            inc.P1_CONTEST_RULE,
            "the benchmark's own invalid rows (ValidConnection FALSE) are in no ladder denominator and "
            "are reported in their own row, as measured.parse_crispri returns them",
            "no rule is excluded from the enumeration: every rule the compiler emits is indexed",
        ],
        "partitions": {
            "ladder": "cached CRISPRi pairs, blind to direction at every step",
            "P1": inc.P1_CALL,
            "P2": inc.P2_CALL,
            "training_and_heldout": inc.P1_SPLIT_RULE,
        },
        "code_cleanliness": mf.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
    }

    path = save_result(args.result, payload)
    print(f"{args.result}: {path}")
    for s in inc.LADDER_STEPS:
        b = steps[s].payload()["outcome_breakdown"]
        print(f"  {s}: {steps[s].n} pairs (increase {b[inc.INCREASE]}, decrease {b[ms.DECREASE]})")
    for name, g in gates.items():
        print(
            f"  {name}: {g['links']} links, {g['independent_loci']} independent loci -> "
            f"{'at or above both floors' if g['meets_both_floors'] else 'no-go'}"
            + (f" ({'; '.join(g['short_by'])})" if g["short_by"] else "")
        )
    print(f"  contested by a regulated pair: {len(contested)}")
    print(f"  verdict: {verdict['go_or_no_go']}")
    print(f"  killed at: {killed['step']} -> {killed['data_or_code']}")


def _killed_at(steps: dict[str, Tally], gates: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """The step a short population is attributed to, by the map registered before the count.

    The earliest ladder step whose significant-increase count is already below the link floor, or the
    contest rule when every step clears the floor but P1 does not. `increases.ATTRIBUTION_OF_STEPS`
    then says whether that step is a property of the data or of this repository's code.
    """
    if all(g["meets_both_floors"] for g in gates.values()):
        return {"step": None, "data_or_code": None, "why": "no population fell short"}
    for s in inc.LADDER_STEPS[1:]:
        if steps[s].outcomes.get(inc.INCREASE, 0) < inc.POSITIVE_FLOOR:
            out = inc.attribution_of(s)
            out["increases_at_this_step"] = steps[s].outcomes.get(inc.INCREASE, 0)
            out["link_floor"] = inc.POSITIVE_FLOOR
            return out
    out = inc.attribution_of("the_contest_rule")
    out["increases_at_the_last_step"] = steps[inc.LADDER_STEPS[-1]].outcomes.get(inc.INCREASE, 0)
    out["link_floor"] = inc.POSITIVE_FLOOR
    return out


if __name__ == "__main__":
    main()
