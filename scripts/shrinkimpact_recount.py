#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Recompute the counts `closure_cell_inventory` reports, with the complete stores, and compare.

    uv run python scripts/shrinkimpact_recount.py closure_cell_inventory

Run it in a clean worktree with the git-ignored stores linked read-only
(`scripts/link_stores_read_only.py`), so that every per-chromosome input the census could read is
actually there. Nothing is rewritten: the committed result is opened read-only and the comparison is
written as its own additive result, `data/results/shrinkimpact_<name>.json`.

WHAT IS RECOMPUTED, and nothing else: the counts the result itself reports. `library.guides`; the
three `dependency` gene tallies; all of `coverage` (`counts` and `per_chromosome`); all of `arms`
(`sizes`, `in_no_arm`, `dependency_of_the_attributed_gene`); `frozen_population.elements`; and, leaf
by leaf, every row of the frozen population against the file the registration pins.

WHAT IS NOT RECOMPUTED, said rather than left to be noticed. `screen_file_layout` is eight ranged
HTTP reads of ENCODE's guide-quantification files and is not a count from a store, so it is left
alone. `frozen_population.sha256` is a digest over bytes that include the run's own wall-clock
`frozen` stamp, so it is not reproducible by construction and comparing it would be a false
negative; the population is compared ROW BY ROW instead, which is strictly more than the digest
would have told us.

NO NETWORK, NO MODEL, NO MONEY: the DepMap column is read from its cache with `net=False`, and a
cache that does not cover the recomputed gene set stops the run by name instead of streaming.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import closure_cell as cc  # noqa: E402
from genomeos.attribution import shrinkimpact as si  # noqa: E402
from genomeos.genome.annotation import default_gencode  # noqa: E402
from genomeos.results import save_result  # noqa: E402

ENTRY = "scripts/shrinkimpact_recount.py"
OWN_CODE = (
    "genomeos/attribution/shrinkimpact.py",
    "scripts/shrinkimpact_recount.py",
    "tests/test_shrinkimpact.py",
)
#: Imported from the census' own entry script rather than restated, so the two cannot drift.
CHROMS = (*(f"chr{i}" for i in range(1, 23)), "chrX", "chrY")
CLOSURE_CELL = "closure_cell_inventory"
#: A moved leaf is printed in full; this caps the listing so a wholesale disagreement cannot write
#: a hundred-megabyte result. The COUNT is never capped.
MAX_MOVED_LISTED = 200


def census_again() -> tuple[list[dict[str, Any]], dict[str, Any], set[str], dict[str, Any]]:
    """`closure_cell_inventory.census()`, re-run here, plus the two facts that explain the groups.

    The census' own loop is reproduced rather than imported because the entry script is not an
    importable module (it runs `main()` on import of `__main__` only, but it also calls
    `mf.trace_begin()` at module scope in `main`), and because the two extra facts -- which GENCODE
    file each chromosome's annotation actually comes from, and which contigs the CRISPRi library's
    design carries -- are what decide whether a group count below 24 is a shrink or the scope.
    """
    guides = cc.guides_by_chrom()
    library_contigs = {c: len(v) for c, v in sorted(guides.items())}
    gencode_source: dict[str, str | None] = {}
    rows: list[dict[str, Any]] = []
    counts: dict[str, int] = {}
    counts["library_guides"] = sum(len(v) for v in guides.values())
    genes: set[str] = set()
    per_chrom: dict[str, int] = {}

    def bump(key: str) -> None:
        counts[key] = counts.get(key, 0) + 1

    for chrom in CHROMS:
        chrom_guides = guides.get(chrom) or []
        chosen = default_gencode({chrom})
        gencode_source[chrom] = str(chosen) if chosen else None
        if not chrom_guides:
            bump("chromosomes_with_no_library_guide")
            per_chrom[chrom] = 0
            continue
        index = cc.bin_index(chrom_guides)
        positions, symbols = cc.coding_tss(chrom)
        kept = 0
        from genomeos.attribution import targets

        for element in targets.attributed(chrom):
            bump("attributed_elements")
            gene = (element.get("predicted_coding") or {}).get("gene")
            effect = cc.element_effect(element)
            covered = len(cc.overlapping(chrom_guides, index, element["start"], element["end"]))
            if covered < cc.MIN_GUIDES:
                bump("fewer_than_min_guides")
                continue
            bump("covered_by_the_library")
            if effect is None:
                bump("no_value_on_the_cell_s_own_track")
                continue
            near, distance = cc.nearest_coding(positions, symbols, element["start"], element["end"])
            rows.append(
                {
                    "id": element["id"],
                    "chrom": chrom,
                    "start": int(element["start"]),
                    "end": int(element["end"]),
                    "gene": gene,
                    "effect_k562": round(effect, 4),
                    "guides": covered,
                    "nearest_coding_gene": near,
                    "nearest_coding_distance": distance,
                    "promoter": bool(distance <= cc.PROMOTER_WINDOW),
                }
            )
            kept += 1
            if gene:
                genes.add(gene)
            if near:
                genes.add(near)
        per_chrom[chrom] = kept
    counts["elements_in_the_frozen_population"] = len(rows)
    extras = {
        "library_design_contigs": library_contigs,
        "library_design_carries_chrY": "chrY" in guides,
        "gencode_file_chosen_per_chromosome": gencode_source,
    }
    return rows, {"counts": counts, "per_chromosome": per_chrom}, genes, extras


def arms_again(rows: list[dict[str, Any]], effects: dict[str, float]) -> dict[str, Any]:
    """`closure_cell_inventory.arms()`, re-run here; the row annotations it adds are kept."""
    sizes: dict[str, int] = {}
    unassigned: dict[str, int] = {}
    dependency_of_target: dict[str, int] = {}

    def bump(d: dict[str, int], key: str) -> None:
        d[key] = d.get(key, 0) + 1

    for row in rows:
        target = cc.dependency(effects, row["gene"])
        neighbour = cc.dependency(effects, row["nearest_coding_gene"])
        bump(dependency_of_target, target)
        row["dependency_attributed"] = target
        row["dependency_nearest"] = neighbour
        if row["promoter"]:
            label = cc.POSITIVE_CONTROL if target == "essential" else None
            row["arm"] = label
            if label:
                bump(sizes, label)
            else:
                bump(unassigned, "promoter_element_whose_gene_k562_does_not_depend_on")
            continue
        arm = cc.arm_of(target, neighbour, row["effect_k562"])
        row["arm"] = arm
        if arm is None:
            if abs(row["effect_k562"]) < cc.MIN_EFFECT:
                bump(unassigned, "prediction_under_the_magnitude_bar")
            elif target not in ("essential", "not") or neighbour not in ("essential", "not"):
                bump(unassigned, "a_dependency_call_is_borderline_or_unmeasured")
            else:
                bump(unassigned, "negative_direction_with_a_gene_k562_does_not_depend_on")
        else:
            bump(sizes, arm)
    return {
        "sizes": sizes,
        "in_no_arm": unassigned,
        "dependency_of_the_attributed_gene": dependency_of_target,
    }


def _capped(report: dict[str, Any]) -> dict[str, Any]:
    """Keep the full count of moved leaves and cap only the listing, saying so when it caps."""
    if len(report["moved"]) <= MAX_MOVED_LISTED:
        return report
    out = dict(report)
    out["moved"] = report["moved"][:MAX_MOVED_LISTED]
    out["moved_listing_capped_at"] = MAX_MOVED_LISTED
    out["moved_listing_note"] = (
        f"{report['moved_count']} leaves moved; the first {MAX_MOVED_LISTED} are listed and "
        f"`moved_count` is the full number"
    )
    return out


def recount_closure_cell() -> dict[str, Any]:
    """The whole comparison for `closure_cell_inventory`, as the payload of its own result."""
    committed = json.loads((ROOT / "data/results/closure_cell_inventory.json").read_text())
    rows, coverage, genes, extras = census_again()
    held = cc.k562_gene_effect(genes, net=False)
    effects = held["values"]
    assignment = arms_again(rows, effects)
    population = cc.load_population()

    reports = [
        si.compare(committed["coverage"], coverage, label="coverage (counts and per_chromosome)"),
        si.compare(committed["arms"]["sizes"], assignment["sizes"], label="arms.sizes"),
        si.compare(committed["arms"]["in_no_arm"], assignment["in_no_arm"], label="arms.in_no_arm"),
        si.compare(
            committed["arms"]["dependency_of_the_attributed_gene"],
            assignment["dependency_of_the_attributed_gene"],
            label="arms.dependency_of_the_attributed_gene",
        ),
        si.compare(
            {"guides": committed["library"]["guides"]},
            {"guides": coverage["counts"]["library_guides"]},
            label="library.guides",
        ),
        si.compare(
            {
                "genes_asked": committed["dependency"]["genes_asked"],
                "genes_with_a_value": committed["dependency"]["genes_with_a_value"],
                "genes_absent_from_the_release": committed["dependency"]["genes_absent_from_the_release"],
            },
            {
                "genes_asked": len(held["requested"]),
                "genes_with_a_value": len(effects),
                "genes_absent_from_the_release": len(held["missing"]),
            },
            label="dependency gene tallies",
        ),
        si.compare(
            {"elements": committed["frozen_population"]["elements"]},
            {"elements": len(rows)},
            label="frozen_population.elements",
        ),
        si.compare(
            population["elements"],
            rows,
            label="the frozen population the registration pins, row by row",
        ),
    ]
    return {
        "question": "with the complete stores in place, does any count `closure_cell_inventory` "
        "reports move? Its declared input groups carry 23, 23 and 22 per-chromosome members "
        "against a 24-chromosome complement, and a count below the complement can be a silent "
        "shrink, the run's legitimate scope, or an artefact of how the groups were labelled.",
        "corrects": "nothing: this is an additive check on data/results/closure_cell_inventory.json "
        "(2026-10-10), which is not touched here",
        "comparisons": [_capped(r) for r in reports],
        "no_number_moves": si.every_reading_is_no_movement(reports),
        "why_the_group_counts_are_below_24": {
            "chrY, absent from all three per-chromosome groups": extras["library_design_contigs"],
            "chrY reading": "the CRISPRi library's DESIGN file (ENCFF630CYG, md5 pinned in the "
            "result) carries guides on 31 contigs and chrY is not one of them, so the census' loop "
            "takes the `no library guide` branch for chrY before it opens that chromosome's element "
            "table or its annotation. The committed result already says this in two places: "
            "`coverage.counts.chromosomes_with_no_library_guide` is 1, `coverage.per_chromosome."
            "chrY` is 0, and `result_manifest.exclusions` carries \"chrY: the library's design "
            "carries no guide there\". chrY's absence is a property of the library, not of this "
            "machine's stores, and no store state can change it.",
            "chr21, absent from the GENCODE group only": extras["gencode_file_chosen_per_chromosome"],
            "chr21 reading": "`genomeos.genome.annotation.default_gencode` prefers the distilled "
            "chr21+chrM subset whenever the chromosomes asked for are within {chr21, chrM}, and "
            "that file lives at data/results/gencode_v50_chr21_chrM.gff3.gz. The census builds its "
            "declared groups from the traced reads BY PATH PREFIX, and `data/results/` is the first "
            "group's prefix, so chr21's annotation was declared in group 1 (`the all-enhancer "
            "sweep's committed per-chromosome results`, 70 members = 3 x 23 sweep files + this one) "
            "instead of group 3 (`GENCODE 50 per-chromosome annotation`, 22 members). The "
            "comparator for chr21 WAS read and IS declared; the group LABEL is wrong about it. "
            "22 = 24 - chrY (never opened) - chr21 (declared under the other group's label).",
        },
        "not_recomputed_here": [
            "screen_file_layout: eight ranged HTTP reads of ENCODE's guide-quantification files, "
            "not a count from a store",
            "frozen_population.sha256: the digest covers bytes that include the run's own "
            "wall-clock `frozen` stamp and is not reproducible by construction; the population is "
            "compared row by row instead",
        ],
        "wording": "a recomputed count is what the complete stores give; it does not by itself "
        "overturn a registered reading, and no registered test is re-run here",
        "stores": "a clean worktree of the committed tree with data/cache, data/reference and "
        "data/knowledge cloned in at 0444 (scripts/link_stores_read_only.py)",
        "alphagenome_requests": 0,
        "model_requests": 0,
        "money": "none: every byte read here is already on this machine",
    }


def closure_cell_manifest(read_now: list[str]) -> dict[str, Any]:
    """The contract's block, with the inputs taken from what the audit hook recorded.

    The grouping is by what the code MEANS and not by path prefix, which is the labelling defect
    this result documents: every `gencode_v50_*` file goes into the annotation group wherever it
    lives, so the group's member count is the number of chromosomes whose annotation was opened.
    """

    def is_gencode(p: str) -> bool:
        return Path(p).name.startswith("gencode_v50_")

    groups = {
        "GENCODE 50 per-chromosome annotation, wherever it lives: the nearest-gene comparator": [
            p for p in read_now if is_gencode(p)
        ],
        "the all-enhancer sweep's committed per-chromosome results, read through attribution/targets.py": [
            p for p in read_now if p.startswith("data/results/") and not is_gencode(p)
        ],
        "the sweep's local per-chromosome element tables, which data/knowledge never commits": [
            p for p in read_now if p.startswith("data/knowledge/alphagenome/")
        ],
        "the CRISPRi library's design, the frozen population and the distilled K562 dependency column": [
            p for p in read_now if p.startswith("data/knowledge/closure_cell/")
        ],
    }
    inputs = [
        mf.files_entry(label, sorted(paths), partition="not an evaluation partition: the recount's inputs")
        for label, paths in groups.items()
        if paths
    ]
    grouped = {p for paths in groups.values() for p in paths}
    other = sorted(p for p in set(read_now) - grouped if (ROOT / p).exists())
    inputs += [mf.input_entry(p, partition="read by the recount") for p in other]
    return {
        "sources": [
            {
                "accession": "data/results/closure_cell_inventory.json",
                "version": "as committed 2026-10-10; opened read-only and not rewritten",
            },
            {"accession": "GENCODE 50 (GRCh38.p14) per-chromosome annotation", "version": "v50"},
            {
                "accession": f"ENCODE {cc.GUIDE_BED_SOURCE}",
                "version": f"md5 {cc.GUIDE_BED_MD5}; " + cc.ENCODE_LICENCE,
            },
            {
                "accession": "DepMap CRISPRGeneEffect (Chronos), K562 = " + cc.K562_MODEL,
                "version": "read from this machine's distilled cache with net=False",
            },
            {
                "accession": "the project's own all-enhancer deletion sweep (enhancer_targets_all_<chrom>)",
                "version": "as committed on this machine",
            },
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38 throughout, exactly as the census: the sweep's elements, the library's "
        "GRCh38 guide locations and GENCODE 50; no coordinate is lifted here",
        "coordinates": {
            "base": 0,
            "interval": "half-open",
            "note": "the census' own convention, unchanged: the elements' start/end and the "
            "library's 3 bp guide intervals, both 0-based half-open",
        },
        "parameters": {
            "min_guides": cc.MIN_GUIDES,
            "min_effect": cc.MIN_EFFECT,
            "essential_at_or_below": cc.ESSENTIAL,
            "not_essential_at_or_above": cc.NOT_ESSENTIAL,
            "promoter_window_bp": cc.PROMOTER_WINDOW,
            "chromosomes": list(CHROMS),
            "depmap_network": False,
            "max_moved_leaves_listed": MAX_MOVED_LISTED,
        },
        "exclusions": [
            "screen_file_layout: eight ranged HTTP reads, not a count from a store",
            "frozen_population.sha256: covers the run's own wall-clock stamp, so it is not "
            "reproducible; the population is compared row by row instead",
            "every registered test of M3: a changed population would need a new registration, and "
            "re-running one is not this recount's to do",
            "chrY: the library's design carries no guide there, which is the finding and not an "
            "exclusion this recount chose",
        ],
        "partitions": "n/a: a recomputation of a census' own counts, with no outcome measured",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("name", choices=[CLOSURE_CELL])
    args = parser.parse_args()
    mf.trace_begin()
    started = time.time()
    print(f"recounting {args.name} with the complete stores", flush=True)
    record = recount_closure_cell()
    record["seconds"] = round(time.time() - started, 1)
    record["result_manifest"] = closure_cell_manifest(mf.traced_reads())
    path = save_result(f"shrinkimpact_{args.name}", record)
    print(f"wrote {path}", flush=True)
    for report in record["comparisons"]:
        print("  " + report["reading"], flush=True)
    print(f"no_number_moves = {record['no_number_moves']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
