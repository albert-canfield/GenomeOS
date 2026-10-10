#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The search for what follows M3: every public element-level fitness screen, counted and gated.

What this opens: the project's own attributed elements, seven CRISPR library DESIGN files (guide and
deletion positions, which carry no measurement), GENCODE for the nearest coding gene, and one
DepMap 24Q4 row per reachable cell -- the chain's premise ("the gene matters to this cell"), not its
outcome.

What this does NOT open: any guide count or phenotype value in any of the screens. Not one byte of a
guide-quantification file is fetched here. The accessions and sizes in the record come from the
ENCODE portal's metadata.

What it writes: `data/results/closure_cell2_inventory.json` -- every candidate with its licence,
size, cell, assay, the compiled elements its library covers and the arm sizes M3's rule would give
in that cell, each read against the gate fixed in `attribution/closure_cell2.py` before the counting.
If no candidate meets the gate the inventory is the lane's result and nothing is registered.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import closure_cell as cc  # noqa: E402
from genomeos.attribution import closure_cell2 as cc2  # noqa: E402
from genomeos.attribution import targets  # noqa: E402
from genomeos.results import save_result  # noqa: E402
from genomeos.therapeutics import depmap  # noqa: E402

NAME = "closure_cell2_inventory"
ENTRY = "scripts/closure_cell2_inventory.py"
OWN_CODE = (
    "genomeos/attribution/closure_cell2.py",
    "scripts/closure_cell2_inventory.py",
    "tests/test_closure_cell2.py",
)
CHROMS = (*(f"chr{i}" for i in range(1, 23)), "chrX", "chrY")


def covered_elements(
    designs: dict[str, dict[str, list[tuple[int, int, None]]]], log=None
) -> tuple[list[dict[str, Any]], dict[str, int], set[str]]:
    """Every attributed element at least one candidate library touches, with everything but outcome.

    One pass over the chromosomes, because the alternative -- one pass per library -- reads the
    sweep's per-chromosome tables seven times. An element is kept when any library has a row
    overlapping it; the per-library guide counts are recomputed per candidate from these rows, which
    is cheap once the 440,000 are down to the covered few.
    """
    rows: list[dict[str, Any]] = []
    counts = {"attributed_elements": 0, "touched_by_some_library": 0}
    genes: set[str] = set()
    for chrom in CHROMS:
        held = {
            key: (design[chrom], cc.bin_index(design[chrom]))
            for key, design in designs.items()
            if design.get(chrom)
        }
        if not held:
            continue
        positions, symbols = cc.coding_tss(chrom)
        for element in targets.attributed(chrom):
            counts["attributed_elements"] += 1
            if not any(
                cc.overlapping(intervals, index, element["start"], element["end"])
                for intervals, index in held.values()
            ):
                continue
            counts["touched_by_some_library"] += 1
            near, distance = cc.nearest_coding(positions, symbols, element["start"], element["end"])
            gene = (element.get("predicted_coding") or {}).get("gene")
            rows.append(
                {
                    "id": element["id"],
                    "chrom": chrom,
                    "start": int(element["start"]),
                    "end": int(element["end"]),
                    "gene": gene,
                    "effects": {c: cc2.effect_in(element, c) for c in cc2.CELLS},
                    "nearest_coding_gene": near,
                    "nearest_coding_distance": distance,
                }
            )
            for name in (gene, near):
                if name:
                    genes.add(name)
        if log:
            print(f"  {chrom}: {counts['touched_by_some_library']} covered so far", file=log, flush=True)
    return rows, counts, genes


def dependency_columns(genes: set[str], log=None) -> dict[str, dict[str, Any]]:
    """One DepMap Chronos row per reachable cell, K562's through M3's own reader."""
    out: dict[str, dict[str, Any]] = {}
    for cell, model in cc2.CELLS.items():
        held = cc.k562_gene_effect(genes) if cell == cc.CELL else cc2.gene_effect_row(model, genes)
        out[cell] = held
        if log:
            print(
                f"  {cell} ({model}): {len(held['values'])} of {len(genes)} genes, "
                f"from_cache={held.get('from_cache')}",
                file=log,
                flush=True,
            )
    return out


def result_manifest(read_now: list[str], designs: list[dict[str, Any]]) -> dict[str, Any]:
    """The contract's block, with the inputs taken from what the audit hook recorded."""
    groups = {
        "the all-enhancer sweep's committed per-chromosome results, read through attribution/targets.py": [
            p for p in read_now if p.startswith("data/results/")
        ],
        "the sweep's local per-chromosome element tables, which data/knowledge never commits": [
            p for p in read_now if p.startswith("data/knowledge/alphagenome/")
        ],
        "GENCODE 50 per-chromosome annotation, for the nearest protein-coding TSS": [
            p for p in read_now if p.startswith("data/reference/gencode")
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
        mf.files_entry(label, sorted(paths), partition="not an evaluation partition: the search's inputs")
        for label, paths in groups.items()
        if paths
    ]
    grouped = {p for paths in groups.values() for p in paths}
    other = sorted(p for p in set(read_now) - grouped if (ROOT / p).exists())
    inputs += [mf.input_entry(p, partition="read by the search") for p in other]
    sources = [
        {
            "accession": f"{d['reference']}, file {d['accession']}",
            "version": f"md5 {d['md5']}; {cc2.ENCODE_LICENCE}",
        }
        for d in designs
    ]
    sources += [
        {
            "accession": "ENCODE4 functional characterisation: the candidate screens' "
            "guide-quantification files, named by accession and NOT fetched here ("
            + ", ".join(f"{c['key']}={c['screen']}" for c in cc2.CANDIDATES)
            + ")",
            "version": "ENCODE portal metadata, released; " + cc2.ENCODE_LICENCE,
        },
        {
            "accession": "DepMap CRISPRGeneEffect (Chronos), one row per reachable cell: "
            + ", ".join(f"{k}={v}" for k, v in cc2.CELLS.items()),
            "version": f"{depmap.RELEASE}, doi:{depmap.DOI}; {depmap.LICENCE}",
        },
        {
            "accession": "the project's own all-enhancer deletion sweep "
            "(enhancer_targets_all_<chrom>), per-cell values on each cell's own RNA track",
            "version": "as committed on this machine",
        },
        {"accession": "GENCODE 50 (GRCh38.p14) per-chromosome annotation", "version": "v50"},
        {
            "accession": "M3's own result, for the numbers the gate is set against",
            "version": cc2.M3_RESULT + " (committed 80bf3db)",
        },
    ]
    return {
        "sources": sources,
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38: the sweep's elements and every candidate library's published guide or "
        "deletion positions are on one build, and no coordinate is lifted here",
        "coordinates": {
            "base": 0,
            "interval": "half-open",
            "note": "the elements' own start/end and the libraries' intervals, both 0-based "
            "half-open; a library row counts for an element when the two overlap",
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
        },
        "exclusions": [
            "every guide count and every phenotype value in every candidate screen: this search "
            "reads positions and the dependency column only, and fetches no readout file",
            "every element with fewer than MIN_GUIDES library rows inside it, per candidate",
            "every screen in a cell the sweep never scored: the project's own claim cannot be made "
            "there, and those candidates are recorded with that reason instead of a count",
            "the sweep's two cells without a DepMap 24Q4 column (GM12878, IMR-90)",
        ],
        "partitions": "n/a: a search and a gate, with no outcome measured",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    mf.trace_begin()
    started = time.time()
    log = sys.stdout
    designs: dict[str, dict[str, list[tuple[int, int, None]]]] = {}
    design_records: list[dict[str, Any]] = []
    for accession in cc2.LIBRARIES:
        record = cc2.fetch_design(accession)
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
        print(
            f"{accession}: {census['counted']} rows counted of {census['rows']} "
            f"({census['not_targeting']} not targeting)",
            file=log,
            flush=True,
        )
    rows, counts, genes = covered_elements(designs, log=log)
    print(
        f"{len(rows)} attributed elements touched by some library; {len(genes)} genes",
        file=log,
        flush=True,
    )
    columns = dependency_columns(genes, log=log)

    results: list[dict[str, Any]] = []
    for candidate in cc2.CANDIDATES:
        entry: dict[str, Any] = {k: v for k, v in candidate.items() if k != "library"}
        entry["library"] = next(d for d in design_records if d["accession"] == candidate["library"])
        if not candidate["reachable"]:
            entry["coverage"] = "not counted: the chain cannot be run in this cell"
            entry["qualifies"] = False
            results.append(entry)
            continue
        cell = candidate["cell"]
        counted = cc2.arms_over(rows, designs[candidate["library"]], columns[cell]["values"], cell)
        held = columns[cell]["values"]
        asked_with_a_value = sum(1 for g in genes if g.upper() in held)
        entry["dependency"] = {
            "cell": cell,
            "model": cc2.CELLS[cell],
            "release": columns[cell]["release"],
            "genes_with_a_value": len(columns[cell]["values"]),
            "genes_absent_from_the_release": len(columns[cell]["missing"]),
            "genes_this_search_asked_about": len(genes),
            "of_those_with_a_value_in_this_cell": asked_with_a_value,
            "of_those_absent_from_the_release": len(genes) - asked_with_a_value,
            "genes_in_the_cache_for_this_cell": len(columns[cell]["requested"]),
            "cache_denominator_note": "the two counts named genes_with_a_value and "
            "genes_absent_from_the_release are over the CACHE, which is reused between runs and "
            "holds more genes for K562 (M3's census asked about them) than for HepG2; the three "
            "counts over this search's own gene set are the ones to compare between cells",
            "streamed_sha256": columns[cell].get("streamed_sha256"),
        }
        entry["coverage"] = counted
        entry["qualifies"] = counted["qualifies"]
        results.append(entry)
        print(
            f"{candidate['key']}: {counted['counts']['elements_with_min_guides']} elements, "
            f"arm A {counted['arm_a']}, floor {counted['floor_arm']}, "
            f"detectable {counted['smallest_detectable_separation']}, "
            f"qualifies {counted['qualifies']}",
            file=log,
            flush=True,
        )

    qualifying = [r["key"] for r in results if r.get("qualifies")]
    record: dict[str, Any] = {
        "question": "is there a public element-level fitness screen, other than the one M3 read, on "
        "which the cell-level chain could be run over a LARGER or BETTER-DECIDED population than "
        "M3's 35 elements? Counted from positions and the dependency column alone, before any guide "
        "count or phenotype value is read.",
        "what_the_chain_needs_in_one_cell": [
            "the project's own claim on that cell's RNA track: the sweep scored K562, HepG2, "
            "GM12878 and IMR-90 and no other cell",
            "a DepMap 24Q4 Chronos column for that cell: K562 (ACH-000551) and HepG2 (ACH-000739) "
            "have one, GM12878 and IMR-90 do not",
            "a public element-level fitness or growth readout in that cell that is not M3's screen",
        ],
        "cells_the_sweep_scored_without_a_dependency_column": cc2.CELLS_WITHOUT_A_DEPENDENCY_COLUMN,
        "other_cells_of_the_myc_myb_screen_the_chain_cannot_reach": (cc2.OTHER_CELLS_OF_THE_MYC_MYB_SCREEN),
        "gate_fixed_before_the_counting": cc2.GATE,
        "arm_keys_are_m3s": cc2.POSITIVE_CONTROL_READS,
        "m3_for_comparison": {
            "result": cc2.M3_RESULT,
            "arm_a": cc2.M3_ARM_A,
            "floor_arm": 5427,
            "smallest_detectable_separation": cc2.M3_DETECTABLE,
            "reading": "No difference was detected at the cell level",
        },
        "elements": {**counts, "genes_asked_of_depmap": len(genes)},
        "candidates": results,
        "candidates_considered": len(results),
        "candidates_the_chain_can_reach": sum(1 for c in cc2.CANDIDATES if c["reachable"]),
        "qualifying": qualifying,
        "outcome": "A_CANDIDATE_QUALIFIES" if qualifying else "NO_CANDIDATE_QUALIFIES",
        "what_would_change_the_answer": [
            "a public element-level fitness screen in HepG2 beyond the MYC and MYB loci, which "
            "would carry the chain in a second cell the sweep scored",
            "a sweep run over a cell that already has both a fitness screen and a DepMap column "
            "(OCI-AML2, HCT116 and the eight other cells of the MYC/MYB screen are each blocked "
            "only by the missing track)",
            "a dependency call that admits the borderline band, which is the largest single reason "
            "an element covered by a library sits in no arm",
        ],
        "not_read_here": [
            "every guide count and phenotype value in every candidate screen",
            "M3's own eight guide-quantification files",
            "the benchmark's expression effect sizes under data/knowledge/crispri",
        ],
        "alphagenome_requests": 0,
        "money": "none: every byte here is public and already paid for",
        "seconds": round(time.time() - started, 1),
    }
    record["result_manifest"] = result_manifest(mf.traced_reads(), design_records)
    path = save_result(NAME, record)
    print(f"wrote {path}", file=log)
    print(f"outcome: {record['outcome']}; qualifying: {qualifying or 'none'}", file=log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
