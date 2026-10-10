#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The census for M3, taken before any cell-level outcome is read, and the population it freezes.

What this opens: the project's own attributed elements, the CRISPRi library's DESIGN (guide target
positions, which carry no measurement), GENCODE for the nearest coding gene, and one row of DepMap
24Q4 -- K562's Chronos gene effect, which is the premise of the chain ("the gene matters to this
cell") and not its outcome.

What this does NOT open: any count in the screen's guide-quantification files. Those eight files
are touched only for their first `PROBE_BYTES` bytes, and only to fix the column layout -- the
number of fields per row, that field 5 parses as a number and that field 16 holds a guide type. No
count value, no element total and no arm figure is formed here, and the registration that comes
next is written and committed before the first count is read.

What it writes: `data/results/closure_cell_inventory.json`, the census with the arm sizes known in
advance, and `data/knowledge/closure_cell/closure_cell_population.json.gz`, the frozen population,
whose sha256 the registration pins so the measured set cannot move afterwards.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
import urllib.request
from collections import Counter
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import closure_cell as cc  # noqa: E402
from genomeos.attribution import targets  # noqa: E402
from genomeos.results import save_result  # noqa: E402
from genomeos.therapeutics import depmap  # noqa: E402

NAME = "closure_cell_inventory"
ENTRY = "scripts/closure_cell_inventory.py"
OWN_CODE = (
    "genomeos/attribution/closure_cell.py",
    "scripts/closure_cell_inventory.py",
    "scripts/closure_cell_register.py",
    "scripts/closure_cell_run.py",
    "tests/test_closure_cell.py",
)
#: Enough of each file to see the layout repeat on thousands of rows, and far too little to form
#: any aggregate: the first 2 MB of a 173 MB file.
PROBE_BYTES = 2 << 20
CHROMS = (*(f"chr{i}" for i in range(1, 23)), "chrX", "chrY")


def probe_layout(accession: str) -> dict[str, Any]:
    """The column layout of one guide-quantification file, from its first `PROBE_BYTES` bytes.

    Deliberately narrow. What comes back is the distribution of field counts, whether field 5
    parses as a number on every row, whether every target interval is 3 bp and which guide types
    appear -- facts about the SHAPE of the file. No value from field 5 is kept or returned.
    """
    url = cc.ENCODE_FILE.format(a=accession)
    request = urllib.request.Request(
        url, headers={"User-Agent": "GenomeOS/0.1 (closure-cell)", "Range": f"bytes=0-{PROBE_BYTES - 1}"}
    )
    with urllib.request.urlopen(request, timeout=300) as response:  # noqa: S310
        raw = response.read(PROBE_BYTES)
    lines = raw.decode("utf-8", "replace").split("\n")[:-1]  # the last line may be cut in half
    fields = Counter(len(line.split("\t")) for line in lines if line.strip())
    types = Counter()
    numeric = 0
    three_bp = 0
    parsed = 0
    for line in lines:
        row = cc.parse_guide_row(line)
        if row is None:
            continue
        parsed += 1
        chrom, start, end, _count, guide_type = row
        types[guide_type] += 1
        numeric += 1
        three_bp += int(end - start == 3)
    return {
        "accession": accession,
        "url": url,
        "bytes_read": len(raw),
        "rows_examined": len([line for line in lines if line.strip()]),
        "field_counts": dict(sorted(fields.items())),
        "rows_matching_the_layout": parsed,
        "rows_whose_count_field_is_numeric": numeric,
        "rows_whose_target_interval_is_3bp": three_bp,
        "guide_types": dict(types),
        "note": "the shape of the file only: no value of the count field is kept here",
    }


def census() -> tuple[list[dict[str, Any]], dict[str, Any], set[str]]:
    """Every attributed element covered by the library, with everything but its outcome."""
    guides = cc.guides_by_chrom()
    rows: list[dict[str, Any]] = []
    counts: Counter[str] = Counter()
    counts["library_guides"] = sum(len(v) for v in guides.values())
    genes: set[str] = set()
    per_chrom: dict[str, int] = {}
    for chrom in CHROMS:
        chrom_guides = guides.get(chrom) or []
        if not chrom_guides:
            counts["chromosomes_with_no_library_guide"] += 1
            per_chrom[chrom] = 0
            continue
        index = cc.bin_index(chrom_guides)
        positions, symbols = cc.coding_tss(chrom)
        kept = 0
        for element in targets.attributed(chrom):
            counts["attributed_elements"] += 1
            gene = (element.get("predicted_coding") or {}).get("gene")
            effect = cc.element_effect(element)
            covered = len(cc.overlapping(chrom_guides, index, element["start"], element["end"]))
            if covered < cc.MIN_GUIDES:
                counts["fewer_than_min_guides"] += 1
                continue
            counts["covered_by_the_library"] += 1
            if effect is None:
                counts["no_value_on_the_cell_s_own_track"] += 1
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
    return rows, {"counts": dict(counts), "per_chromosome": per_chrom}, genes


def arms(rows: list[dict[str, Any]], effects: dict[str, float]) -> dict[str, Any]:
    """The arm sizes, stated in advance, and every reason an element is in none of them."""
    sizes: Counter[str] = Counter()
    unassigned: Counter[str] = Counter()
    dependency_of_target: Counter[str] = Counter()
    for row in rows:
        target = cc.dependency(effects, row["gene"])
        neighbour = cc.dependency(effects, row["nearest_coding_gene"])
        dependency_of_target[target] += 1
        row["dependency_attributed"] = target
        row["dependency_nearest"] = neighbour
        if row["promoter"]:
            label = cc.POSITIVE_CONTROL if target == "essential" else None
            row["arm"] = label
            if label:
                sizes[label] += 1
            else:
                unassigned["promoter_element_whose_gene_k562_does_not_depend_on"] += 1
            continue
        arm = cc.arm_of(target, neighbour, row["effect_k562"])
        row["arm"] = arm
        if arm is None:
            if abs(row["effect_k562"]) < cc.MIN_EFFECT:
                unassigned["prediction_under_the_magnitude_bar"] += 1
            elif target not in ("essential", "not") or neighbour not in ("essential", "not"):
                unassigned["a_dependency_call_is_borderline_or_unmeasured"] += 1
            else:
                unassigned["negative_direction_with_a_gene_k562_does_not_depend_on"] += 1
        else:
            sizes[arm] += 1
    return {
        "sizes": dict(sizes),
        "in_no_arm": dict(unassigned),
        "dependency_of_the_attributed_gene": dict(dependency_of_target),
        "meaning": cc.ARM_MEANING,
    }


def result_manifest(read_now: list[str]) -> dict[str, Any]:
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
        "the CRISPRi library's design and the distilled K562 dependency column": [
            p for p in read_now if p.startswith("data/knowledge/closure_cell/")
        ],
    }
    inputs = [
        mf.files_entry(label, sorted(paths), partition="not an evaluation partition: the census' inputs")
        for label, paths in groups.items()
        if paths
    ]
    grouped = {p for paths in groups.values() for p in paths}
    other = sorted(p for p in set(read_now) - grouped if (ROOT / p).exists())
    inputs += [mf.input_entry(p, partition="read by the census") for p in other]
    return {
        "sources": [
            {"accession": cc.SCREEN_SOURCE, "version": "ENCODE portal, released; " + cc.ENCODE_LICENCE},
            {"accession": cc.GUIDE_BED_SOURCE, "version": f"md5 {cc.GUIDE_BED_MD5}; " + cc.ENCODE_LICENCE},
            {
                "accession": "DepMap CRISPRGeneEffect (Chronos), K562 = " + cc.K562_MODEL,
                "version": f"{depmap.RELEASE}, doi:{depmap.DOI}; {depmap.LICENCE}",
            },
            {
                "accession": "the project's own all-enhancer deletion sweep "
                "(enhancer_targets_all_<chrom>), per-cell values on K562's own RNA track",
                "version": "as committed on this machine",
            },
            {"accession": "GENCODE 50 (GRCh38.p14) per-chromosome annotation", "version": "v50"},
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38: the sweep's elements, the library's GRCh38 guide locations and the "
        "GRCh38 guide quantifications are all on one build, and no coordinate is lifted here",
        "coordinates": {
            "base": 0,
            "interval": "half-open",
            "note": "the elements' own start/end and the library's 3 bp guide intervals, both "
            "0-based half-open; a guide counts for an element when the two overlap",
        },
        "parameters": {
            "min_guides": cc.MIN_GUIDES,
            "min_effect": cc.MIN_EFFECT,
            "essential_at_or_below": cc.ESSENTIAL,
            "not_essential_at_or_above": cc.NOT_ESSENTIAL,
            "promoter_window_bp": cc.PROMOTER_WINDOW,
            "probe_bytes_per_screen_file": PROBE_BYTES,
            "chromosomes": list(CHROMS),
        },
        "exclusions": [
            "every element with fewer than MIN_GUIDES of the library's guides inside it: the screen "
            "cannot speak about it",
            "every element with no value on K562's own track in the sweep",
            "every count in the screen's guide-quantification files: this census reads the first "
            f"{PROBE_BYTES} bytes of each for the column layout and forms no aggregate from them",
            "chrY: the library's design carries no guide there",
        ],
        "partitions": "n/a: a census and a frozen population, with no outcome measured",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    mf.trace_begin()
    started = time.time()
    log = sys.stdout
    print("census over the project's attributed elements and the library's design", file=log, flush=True)
    rows, counts, genes = census()
    print(f"  {len(rows)} elements covered, {len(genes)} genes to ask DepMap about", file=log, flush=True)
    held = cc.k562_gene_effect(genes)
    effects = held["values"]
    print(
        f"  K562 dependency for {len(effects)} of {len(genes)} genes "
        f"({len(held['missing'])} not in the release)",
        file=log,
        flush=True,
    )
    assignment = arms(rows, effects)
    layouts = [probe_layout(a) for a in (*cc.SCREEN["files"].values(), *cc.CONTROL["files"].values())]
    digest = cc.save_population(
        rows,
        {
            "frozen": time.strftime("%Y-%m-%d %H:%M:%S"),
            "cell": cc.CELL,
            "thresholds": {
                "min_guides": cc.MIN_GUIDES,
                "min_effect": cc.MIN_EFFECT,
                "essential": cc.ESSENTIAL,
                "not_essential": cc.NOT_ESSENTIAL,
                "promoter_window": cc.PROMOTER_WINDOW,
            },
            "depmap_release": held["release"],
            "depmap_streamed_sha256": held.get("streamed_sha256"),
        },
    )
    record: dict[str, Any] = {
        "question": "when the project says element e drives gene G on K562's own track, and K562 "
        "depends on G, does silencing e in K562 cost the cell its place in a pooled culture? This "
        "file is the census and the frozen population, taken and written before any count is read.",
        "cell": cc.CELL,
        "screen": {**cc.SCREEN, "control": cc.CONTROL, "licence": cc.ENCODE_LICENCE},
        "library": {
            "design_file": str(cc.GUIDE_BED),
            "md5": hashlib.md5(cc.GUIDE_BED.read_bytes()).hexdigest(),  # noqa: S324
            "md5_published_by_encode": cc.GUIDE_BED_MD5,
            "guides": counts["counts"].get("library_guides"),
            "source": cc.GUIDE_BED_SOURCE,
        },
        "dependency": {
            "release": held["release"],
            "doi": held["doi"],
            "licence": held["licence"],
            "model": cc.K562_MODEL,
            "genes_asked": len(held["requested"]),
            "genes_with_a_value": len(effects),
            "genes_absent_from_the_release": len(held["missing"]),
            "streamed_bytes": held.get("streamed_bytes"),
            "streamed_sha256": held.get("streamed_sha256"),
            "from_cache": held.get("from_cache"),
        },
        "coverage": counts,
        "arms": assignment,
        "screen_file_layout": layouts,
        "frozen_population": {
            "path": str(cc.POPULATION),
            "sha256": digest,
            "elements": len(rows),
            "why": "the registration pins this digest, so the set that is measured cannot be the "
            "set the measurement suggested",
        },
        "not_read_here": [
            "every count in the eight guide-quantification files beyond the layout probe",
            "every other ENCODE CRISPR screen",
            "the benchmark's expression effect sizes under data/knowledge/crispri",
        ],
        "alphagenome_requests": 0,
        "money": "none: every byte here is public and already paid for",
        "seconds": round(time.time() - started, 1),
    }
    record["result_manifest"] = result_manifest(mf.traced_reads())
    path = save_result(NAME, record)
    print(f"wrote {path}", file=log)
    print(f"arm sizes: {json.dumps(assignment['sizes'], sort_keys=True)}", file=log)
    print(f"population {digest} over {len(rows)} elements", file=log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
