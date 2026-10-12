#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Recompute the counts `closure_cell2_inventory` reports, with the complete stores, and compare.

    uv run python scripts/shrinkimpact_cc2.py

Run it in a clean worktree with the git-ignored stores linked read-only
(`scripts/link_stores_read_only.py`). Nothing is rewritten: the committed result is opened
read-only and the comparison is written as its own additive result,
`data/results/shrinkimpact_closure_cell2_inventory.json`.

WHY A SECOND SCRIPT and not a mode of `scripts/shrinkimpact_recount.py`: adding a mode there means
changing its `choices` line, and `commit_own.sh`'s stale-base guard refuses a commit that removes a
line an earlier commit of this lane added. This lane does not use --force, so each result family
gets its own entry script. That also keeps `code.argv` in each result naming exactly the code that
wrote it.

WHAT IS RECOMPUTED, and nothing else: the counts the result itself reports. `elements`
(`attributed_elements`, `touched_by_some_library`, `genes_asked_of_depmap`); for every candidate,
its `library` row census, median interval width and chromosome count, and its whole `coverage`
block -- `counts`, `sizes`, `in_no_arm`, `arm_a`, `floor_arm`, `smallest_detectable_separation` and
`qualifies`; and for every reachable candidate the three `dependency` tallies that are over THIS
SEARCH'S OWN gene set (`genes_this_search_asked_about`, `of_those_with_a_value_in_this_cell`,
`of_those_absent_from_the_release`).

THE TWO CACHE TALLIES ARE RECORDED, NOT COMPARED. `dependency.genes_with_a_value`,
`genes_absent_from_the_release` and `genes_in_the_cache_for_this_cell` are read off the distilled
DepMap cache, which is reused between runs and accumulates -- the committed result says so itself in
its own `cache_denominator_note`. They are not counts over this search's population, and
`cache_tallies_recorded_not_compared` carries the committed and the current values side by side.

NO NETWORK, NO MODEL, NO MONEY: every DepMap row is read from its cache with `net=False`, and a
cache that does not cover the recomputed gene set stops the run by name instead of streaming.
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import closure_cell as cc  # noqa: E402
from genomeos.attribution import closure_cell2 as cc2  # noqa: E402
from genomeos.attribution import shrinkimpact as si  # noqa: E402
from genomeos.attribution import targets  # noqa: E402
from genomeos.genome.annotation import default_gencode  # noqa: E402
from genomeos.results import save_result  # noqa: E402

NAME = "shrinkimpact_closure_cell2_inventory"
COMMITTED = "data/results/closure_cell2_inventory.json"
ENTRY = "scripts/shrinkimpact_cc2.py"
OWN_CODE = (
    "genomeos/attribution/shrinkimpact.py",
    "scripts/shrinkimpact_cc2.py",
    "tests/test_shrinkimpact.py",
)
CHROMS = (*(f"chr{i}" for i in range(1, 23)), "chrX", "chrY")
MAX_MOVED_LISTED = 200


def _inventory_module() -> Any:
    """`scripts/closure_cell2_inventory.py` as a module, so its own `covered_elements` is reused.

    Loaded by path rather than imported, because `scripts/` is not a package. Its `main()` runs
    only under `__main__`, so importing it is side-effect free apart from its module-level
    constants.
    """
    path = ROOT / "scripts/closure_cell2_inventory.py"
    spec = importlib.util.spec_from_file_location("closure_cell2_inventory", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def dependency_columns_offline(genes: set[str]) -> dict[str, dict[str, Any]]:
    """`closure_cell2_inventory.dependency_columns`, with the network taken away.

    The committed run let both readers stream if the cache fell short. Here they may not: a lane
    recomputing a count must not be able to buy a different answer, and a cache that does not cover
    the ask is a finding to report rather than a gap to fill.
    """
    out: dict[str, dict[str, Any]] = {}
    for cell, model in cc2.CELLS.items():
        if cell == cc.CELL:
            out[cell] = cc.k562_gene_effect(genes, net=False)
        else:
            out[cell] = cc2.gene_effect_row(model, genes, net=False)
    return out


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


#: The three `dependency` tallies that are over this search's own gene set, and so are counts of
#: the search. The other three are readings of the shared cache and are recorded, not compared.
OWN_SET_TALLIES = (
    "genes_this_search_asked_about",
    "of_those_with_a_value_in_this_cell",
    "of_those_absent_from_the_release",
)
CACHE_TALLIES = (
    "genes_with_a_value",
    "genes_absent_from_the_release",
    "genes_in_the_cache_for_this_cell",
)


def recount() -> dict[str, Any]:
    """The whole comparison, as the payload of its own result."""
    committed = json.loads((ROOT / COMMITTED).read_text())
    inventory = _inventory_module()

    designs: dict[str, dict[str, list[tuple[int, int, None]]]] = {}
    design_records: list[dict[str, Any]] = []
    for accession in cc2.LIBRARIES:
        # net=False: the design files are already on this machine, and a recount must not be able
        # to fetch a different answer than the one it is checking.
        record = cc2.fetch_design(accession, net=False)
        intervals, census = cc2.read_design(accession)
        design_records.append(
            {
                **record,
                "row_census": census,
                "median_interval_width_bp": cc2.median_width(intervals),
                "chromosomes": len(intervals),
            }
        )
        designs[accession] = intervals

    rows, counts, genes = inventory.covered_elements(designs)
    columns = dependency_columns_offline(genes)

    reports = [
        si.compare(
            committed["elements"],
            {
                "attributed_elements": counts["attributed_elements"],
                "touched_by_some_library": counts["touched_by_some_library"],
                "genes_asked_of_depmap": len(genes),
            },
            label="elements (the three census counts)",
        )
    ]
    cache_tallies: list[dict[str, Any]] = []
    by_accession = {d["accession"]: d for d in design_records}
    for candidate in cc2.CANDIDATES:
        key = candidate["key"]
        was = next(c for c in committed["candidates"] if c["key"] == key)
        design = by_accession[candidate["library"]]
        reports.append(
            si.compare(
                {
                    "row_census": was["library"]["row_census"],
                    "median_interval_width_bp": was["library"]["median_interval_width_bp"],
                    "chromosomes": was["library"]["chromosomes"],
                },
                {
                    "row_census": design["row_census"],
                    "median_interval_width_bp": design["median_interval_width_bp"],
                    "chromosomes": design["chromosomes"],
                },
                label=f"{key}: the library design's own census",
            )
        )
        if not candidate["reachable"]:
            continue
        cell = candidate["cell"]
        counted = cc2.arms_over(rows, designs[candidate["library"]], columns[cell]["values"], cell)
        reports.append(si.compare(was["coverage"], counted, label=f"{key}: coverage, sizes and the gate"))
        held = columns[cell]["values"]
        asked_with_a_value = sum(1 for g in genes if g.upper() in held)
        reports.append(
            si.compare(
                {k: was["dependency"][k] for k in OWN_SET_TALLIES},
                {
                    "genes_this_search_asked_about": len(genes),
                    "of_those_with_a_value_in_this_cell": asked_with_a_value,
                    "of_those_absent_from_the_release": len(genes) - asked_with_a_value,
                },
                label=f"{key}: the dependency tallies over this search's OWN gene set",
            )
        )
        cache_tallies.append(
            {
                "candidate": key,
                "cell": cell,
                "committed": {k: was["dependency"][k] for k in CACHE_TALLIES},
                "the_cache_now": {
                    "genes_with_a_value": len(columns[cell]["values"]),
                    "genes_absent_from_the_release": len(columns[cell]["missing"]),
                    "genes_in_the_cache_for_this_cell": len(columns[cell]["requested"]),
                },
            }
        )

    chr_y_elements = sum(1 for _ in targets.attributed("chrY"))
    gencode_source = {c: str(default_gencode({c}) or "") for c in CHROMS}
    return {
        "question": "with the complete stores in place, does any count `closure_cell2_inventory` "
        "reports move? Its GENCODE input group carries 23 per-chromosome members against a "
        "24-chromosome complement, and chr21 -- the project's reference chromosome -- is the one "
        "the group does not name.",
        "corrects": f"nothing: this is an additive check on {COMMITTED} (2026-10-10), which is not "
        "touched here",
        "comparisons": [_capped(r) for r in reports],
        "no_number_moves": si.every_reading_is_no_movement(reports),
        "cache_tallies_recorded_not_compared": {
            "why": "`dependency.genes_with_a_value`, `genes_absent_from_the_release` and "
            "`genes_in_the_cache_for_this_cell` are read off the distilled DepMap cache, which is "
            "reused between runs and accumulates. They are not counts over this search's "
            "population, and the committed result says exactly that in its own "
            "`cache_denominator_note`. Recorded here so the figures are on the record either way.",
            "per_candidate": cache_tallies,
        },
        "why_the_gencode_group_holds_23_and_not_24": {
            "gencode_file_chosen_per_chromosome": gencode_source,
            "reading": "`genomeos.genome.annotation.default_gencode` prefers the distilled "
            "chr21+chrM subset whenever the chromosomes asked for are within {chr21, chrM}, and "
            "that file lives at data/results/gencode_v50_chr21_chrM.gff3.gz. "
            "`scripts/closure_cell2_inventory.py: result_manifest` sorts its traced reads into "
            "groups BY PATH PREFIX with `data/results/` first, so chr21's annotation is declared "
            "in the sweep's group -- which holds 73 members, that is 3 x 24 committed sweep files "
            "plus this one -- and the GENCODE group holds the other 23. 23 = 24 - 1, and the one "
            "is chr21 declared under the other group's label. The comparator for the project's "
            "reference chromosome WAS read; the group label is wrong about where it came from. "
            "This is the same labelling artefact as in closure_cell_inventory, one chromosome "
            "further along because this search covers chrY as well.",
            "chrY_attributed_elements": chr_y_elements,
            "chrY_reading": "this search's libraries DO reach chrY, so all 24 chromosomes are in "
            "every per-chromosome group here except GENCODE's. That is also the arithmetic behind "
            "the two censuses differing: this result's 440,377 attributed elements against "
            "closure_cell_inventory's 440,252 is exactly chrY's own count, measured above.",
        },
        "not_recomputed_here": [
            "every guide count and phenotype value in every candidate screen: the committed run "
            "fetched not one byte of them and neither does this",
            "the ENCODE portal metadata in `readout_files`: sizes and accessions, not counts from a store",
        ],
        "wording": "a recomputed count is what the complete stores give; it does not by itself "
        "overturn a registered reading, and no registered test is re-run here",
        "stores": "a clean worktree of the committed tree with data/cache, data/reference and "
        "data/knowledge cloned in at 0444 (scripts/link_stores_read_only.py)",
        "alphagenome_requests": 0,
        "model_requests": 0,
        "money": "none: every byte read here is already on this machine",
    }


def result_manifest(read_now: list[str]) -> dict[str, Any]:
    """The contract's block, grouped by what the code MEANS and not by path prefix."""

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
        "the candidate libraries' design files and the distilled dependency columns": [
            p
            for p in read_now
            if p.startswith(("data/knowledge/closure_cell2/", "data/knowledge/closure_cell/"))
        ],
        "the DepMap model table, for which cells have a dependency column at all": [
            p for p in read_now if p.startswith("data/knowledge/therapeutics/depmap/")
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
                "accession": COMMITTED,
                "version": "as committed 2026-10-10; opened read-only and not rewritten",
            },
            {"accession": "GENCODE 50 (GRCh38.p14) per-chromosome annotation", "version": "v50"},
            {
                "accession": "the candidate CRISPR libraries' ENCODE design files: "
                + ", ".join(cc2.LIBRARIES),
                "version": "md5 as published by ENCODE, digested in inputs; " + cc2.ENCODE_LICENCE,
            },
            {
                "accession": "DepMap CRISPRGeneEffect (Chronos), one row per reachable cell: "
                + ", ".join(f"{k}={v}" for k, v in cc2.CELLS.items()),
                "version": "read from this machine's distilled caches with net=False",
            },
            {
                "accession": "the project's own all-enhancer deletion sweep (enhancer_targets_all_<chrom>)",
                "version": "as committed on this machine",
            },
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38 throughout, exactly as the search: the sweep's elements and every "
        "candidate library's published intervals; no coordinate is lifted here",
        "coordinates": {
            "base": 0,
            "interval": "half-open",
            "note": "the search's own convention, unchanged: the elements' start/end and the "
            "libraries' intervals, both 0-based half-open",
        },
        "parameters": {
            "min_guides": cc.MIN_GUIDES,
            "min_effect": cc.MIN_EFFECT,
            "essential_at_or_below": cc.ESSENTIAL,
            "not_essential_at_or_above": cc.NOT_ESSENTIAL,
            "promoter_window_bp": cc.PROMOTER_WINDOW,
            "gate_arm_a_at_least": cc2.M3_ARM_A,
            "gate_detectable_at_most": cc2.M3_DETECTABLE,
            "alpha": cc2.ALPHA,
            "power": cc2.POWER,
            "chromosomes": list(CHROMS),
            "depmap_network": False,
            "max_moved_leaves_listed": MAX_MOVED_LISTED,
        },
        "exclusions": [
            "every guide count and phenotype value in every candidate screen, exactly as the "
            "committed run: not one byte of a quantification file is fetched",
            "the three `dependency` tallies read off the shared DepMap cache: recorded under "
            "`cache_tallies_recorded_not_compared`, not compared, because they are not counts "
            "over this search's population and the committed result says so itself",
            "every registered test: a changed population would need a new registration, and "
            "re-running one is not this recount's to do",
        ],
        "partitions": "n/a: a recomputation of a search's own counts, with no outcome measured",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    mf.trace_begin()
    started = time.time()
    print("recounting closure_cell2_inventory with the complete stores", flush=True)
    record = recount()
    record["seconds"] = round(time.time() - started, 1)
    record["result_manifest"] = result_manifest(mf.traced_reads())
    path = save_result(NAME, record)
    print(f"wrote {path}", flush=True)
    for report in record["comparisons"]:
        print("  " + report["reading"], flush=True)
    print(f"no_number_moves = {record['no_number_moves']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
