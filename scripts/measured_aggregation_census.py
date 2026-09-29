# SPDX-License-Identifier: AGPL-3.0-or-later
"""Review item R6, step 1: every aggregation in the measured layer, what reads it, and how often it bites.

    uv run python scripts/measured_aggregation_census.py [--no-save]

`attribution/measured.py` attaches assay results to compiled elements. Several of its fields are
summaries over more than one observation (a maximum, a strongest outcome, a first match), and the
external review of 2026-09-28 asks that none of them silently decide an element's label. This script
does not change any of them. It lists each one with the consumers `git grep` finds, and it counts,
genome-wide and against the local assay files, how many elements each summary is actually taken over
more than one observation for, and how many of those observations disagree. No requests, no model.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from collections import Counter
from typing import Any

from genomeos import manifest as mf
from genomeos.attribution import measured, mpra
from genomeos.attribution.targets import RUNS, attributed
from genomeos.knowledge import satmut as satmut_knowledge
from genomeos.results import RESULTS_DIR, save_result

CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]

# Read by hand from measured.py at 0924a5c and checked against `git grep` for every consumer.
AGGREGATIONS: list[dict[str, Any]] = [
    {
        "id": "A1",
        "field": "measured.lentimpra.activity (and cells_active, cells_silent from it)",
        "where": "measured.py Layer.for_element, lentiMPRA block",
        "operation": "maximum over every tile that meets the overlap rule, per cell",
        "consumers": [
            "measured.agreement: lentimpra verdict is agrees when any cell is in cells_active",
            "agreement verdict and assays_agreeing/disagreeing -> census agrees/disagrees -> "
            "data/results/measured_layer_genome.json",
            "confidence_calibration.elements_of reads agreement['lentimpra'] -> "
            "data/results/confidence_calibration_genome.json",
            "compile.py programme header counts agreeing/disagreeing elements (verdict)",
            "measured.basis_text prints the per-cell value -> the basis: line of every <id>_measured "
            "block in data/knowledge/compiled/*.bio and data/organisms/human/noncoding_chr21.bio",
            "tests/test_measured_layer.py (cells_active on a one-tile fixture)",
        ],
        "status": "R6 changes it: the label comes from a declared rule over the tiles, the maximum is "
        "kept only under a descriptive name",
    },
    {
        "id": "A2",
        "field": "measured.lentimpra.max_activity",
        "where": "measured.py Layer.for_element",
        "operation": "maximum over cells of the per-cell maximum over tiles",
        "consumers": ["none in code; carried in the row only"],
        "status": "R6 keeps it, renamed to say it is descriptive",
    },
    {
        "id": "A3",
        "field": "measured.<assay>.overlap for crispri, lentimpra, vista, satmut",
        "where": "measured.py Layer.for_element and _satmut_of",
        "operation": "maximum reciprocal overlap over the matched intervals",
        "consumers": ["none in code; carried in the row only (CRISPRi pairs keep their own overlap)"],
        "status": "descriptive; R6 adds each lentiMPRA tile's own overlap beside it",
    },
    {
        "id": "A4",
        "field": "measured.lentimpra.elements",
        "where": "measured.py Layer.for_element",
        "operation": "names of the matched tiles, without their values, strand or overlap",
        "consumers": ["none in code"],
        "status": "R6 adds a tiles list: name, interval, strand, overlap, per-cell value",
    },
    {
        "id": "A5",
        "field": "mpra.Element.activity",
        "where": "mpra.parse",
        "operation": "mean of the rows that share an interval in one cell's file (forward and reverse)",
        "consumers": [
            "measured.load_mpra, and every mpra.load caller (measurability, unknown_coverage, "
            "mpra_score, motif_grammar, motif_transfer, syntax_tiling, candidates, benchmark.loci)",
        ],
        "status": "counted below: the mean never fires in the three ENCODE files (one row per interval "
        "per file), and the strand column is dropped; R6 keeps the strand beside the value",
    },
    {
        "id": "A6",
        "field": "measured.crispri.genes_regulated / genes_increased / genes_no_effect_*",
        "where": "measured.py Layer.for_element, CRISPRi block",
        "operation": "strongest outcome any pair of the gene has, across cells and datasets",
        "consumers": [
            "measured.agreement (crispri verdict)",
            "compile.py targets: (genes_regulated_training)",
            "basis_text, census",
        ],
        "status": "out of R6's scope: R1 (6a4a8a9) keeps outcomes_by_cell and every pair, and "
        "basis_text names where cells differ",
    },
    {
        "id": "A7",
        "field": "measured.crispri.outcomes_by_cell",
        "where": "measured.outcomes_by_cell",
        "operation": "strongest outcome within one cell across its datasets",
        "consumers": ["basis_text via context_differences"],
        "status": "out of R6's scope; the pairs stay in the row",
    },
    {
        "id": "A8",
        "field": "rule strength of a compiled experimental rule",
        "where": "measured.rule_links",
        "operation": "largest |effect size| among the training pairs of one (element, gene, cell)",
        "consumers": ["compile.py _measured_blocks: strength: of every experimental rule"],
        "status": "out of R6's scope (the compiler is lane-bridge's); reported to the coordinator",
    },
    {
        "id": "A9",
        "field": "measured.agreement vista",
        "where": "measured.agreement",
        "operation": "agrees when any matched VISTA element is positive",
        "consumers": ["agreement verdict, census, confidence_calibration"],
        "status": "counted below: positive and negative lists are both kept; no element has both",
    },
    {
        "id": "A10",
        "field": "measured.satmut (every field)",
        "where": "measured._satmut_primaries and Layer._satmut_of",
        "operation": "the primary experiment of a locus only (repeats counted, never read); the "
        "first experiment wins a base two loci share; strongest_effect is the largest |effect|",
        "consumers": ["agreement satmut verdict, basis_text, census satmut block, base_level_rows"],
        "status": "R6 reads the repeat experiments per base (SORT1's group holds SORT1-flip, the "
        "same element in the other orientation) and lists where they disagree",
    },
    {
        "id": "A11",
        "field": "row confidence",
        "where": "measured.confidence_of",
        "operation": "the strongest assay present decides (0.9 perturbation or embryo, else 0.75)",
        "consumers": ["compile.py confidence: of every <id>_measured block and rule"],
        "status": "out of R6's scope (R4's confidence lanes); the docstring names the rule",
    },
    {
        "id": "A12",
        "field": "measured.agreement verdict",
        "where": "measured.agreement",
        "operation": "agrees / disagrees / mixed over the assays' own verdicts",
        "consumers": ["compile.py header, census"],
        "status": "keeps conflict as 'mixed'; unchanged",
    },
]

EPISOMAL = (
    "the word 'episomal' describes the ENCODE4 lentiMPRA in measured.py (module docstring), "
    "confidence_calibration.py (ASSAY_QUESTION lentimpra), scripts/unknown_coverage.py, "
    "docs/ATTRIBUTION.md (four places), docs/GRAMMAR-BY-COMPARISON.md, docs/ROADMAP.md and the stored "
    "confidence_calibration_genome.json; the cited assay integrates its reporters by lentivirus "
    "(Agarwal et al. 2025, Nature, https://www.nature.com/articles/s41586-024-08430-9)"
)


def file_census() -> dict[str, Any]:
    """What the three ENCODE element files hold per row, read directly."""
    out: dict[str, Any] = {}
    for cell, acc in mpra.FILES.items():
        p = mpra.KNOWLEDGE / f"{acc}.bed.gz"
        if not p.exists():
            out[cell] = "file not present"
            continue
        keys: Counter = Counter()
        strands: Counter = Counter()
        pq = 0
        with gzip.open(p, "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 7:
                    continue
                keys[(f[0], f[1], f[2])] += 1
                strands[f[5]] += 1
                pq += len(f) > 10 and (f[9] != "-1" or f[10] != "-1")
        out[cell] = {
            "accession": acc,
            "rows": sum(keys.values()),
            "intervals": len(keys),
            "intervals_with_more_than_one_row": sum(1 for v in keys.values() if v > 1),
            "rows_by_strand": dict(strands),
            "rows_with_a_p_or_q_value": pq,
        }
    out["columns_read"] = (
        "7: log2(RNA/DNA). Columns 8 and 9 hold the normalised DNA and RNA counts the ratio is made "
        "from; 10 and 11 (p and q) are -1 throughout; the files carry no replicate-level value, so no "
        "per-tile uncertainty can be read from them"
    )
    return out


def element_census() -> dict[str, Any]:
    """Per chromosome, how many matched elements each summary is taken over more than one observation."""
    lenti: Counter = Counter()
    vista_c: Counter = Counter()
    sat: Counter = Counter()
    examples: list[dict[str, Any]] = []
    t = measured.MPRA_ACTIVE
    for ch in CHROMS:
        els = attributed(ch)
        layer = measured.Layer(
            chrom=ch,
            lentimpra=measured.load_mpra(ch),
            vista=measured.load_vista(ch),
            satmut=measured.load_satmut(ch),
        )
        for e in els:
            s, en = e["start"], e["end"]
            tiles = [
                x for x in layer.near("lentimpra", s, en) if measured.measures(s, en, x.start, x.end, 0.5)
            ]
            if tiles:
                lenti["elements"] += 1
                lenti[f"elements_with_{min(len(tiles), 3)}{'+' if len(tiles) >= 3 else ''}_tiles"] += 1
                if len(tiles) > 1:
                    for cell in sorted({c for x in tiles for c in x.activity}):
                        vals = [x.activity[cell] for x in tiles if cell in x.activity]
                        a = sum(v >= t for v in vals)
                        lenti["cells_over_more_than_one_tile"] += 1
                        if 0 < a < len(vals):
                            lenti["cells_where_tiles_disagree"] += 1
                            kind = (
                                "tie"
                                if 2 * a == len(vals)
                                else "minority"
                                if 2 * a < len(vals)
                                else "majority"
                            )
                            lenti[f"cells_where_tiles_disagree_{kind}_active"] += 1
                            if len(examples) < 5:
                                examples.append({"chrom": ch, "id": e["id"], "cell": cell, "tiles": vals})
            vs = [x for x in layer.near("vista", s, en) if measured.measures(s, en, x.start, x.end, 0.5)]
            if vs:
                vista_c["elements"] += 1
                vista_c["elements_over_more_than_one_vista_element"] += len(vs) > 1
                pos = [x for x in vs if x.status == "positive"]
                vista_c["elements_with_positive_and_negative"] += bool(pos) and len(pos) < len(vs)
            sm = [x for x in layer.near("satmut", s, en) if measured.measures(s, en, x.start, x.end, 0.5)]
            if sm:
                sat["elements"] += 1
                sat["repeat_experiments_not_read"] += sum(x.repeats - 1 for x in sm)
                sat["elements_over_more_than_one_locus"] += len(sm) > 1
    groups = {}
    if satmut_knowledge.DATA_PATH.exists():
        loci = satmut_knowledge.loci(
            satmut_knowledge.by_element(satmut_knowledge.load(satmut_knowledge.DATA_PATH))
        )
        groups = {k: v for k, v in loci.items() if len(v) > 1}
    return {
        "lentimpra": {**dict(lenti), "active_threshold_log2": t, "examples": examples},
        "vista": dict(vista_c),
        "satmut": {**dict(sat), "loci_with_repeat_experiments": groups},
    }


def inputs() -> list[dict[str, Any]]:
    out = [
        mf.input_entry(mpra.KNOWLEDGE / f"{acc}.bed.gz", partition=None)
        for acc in mpra.FILES.values()
        if (mpra.KNOWLEDGE / f"{acc}.bed.gz").exists()
    ]
    for ch in CHROMS:
        for name in RUNS:
            p = RESULTS_DIR / f"{name}_{ch}.json"
            if p.exists():
                out.append(mf.input_entry(p, partition=None))
    if satmut_knowledge.DATA_PATH.exists():
        out.append(mf.input_entry(satmut_knowledge.DATA_PATH, partition=None))
    return out


def manifest() -> dict[str, Any]:
    return {
        "sources": [
            {
                "accession": f"ENCODE {mpra.LIBRARY} ({', '.join(mpra.FILES.values())})",
                "version": "as fetched 2026-09-12",
            },
            {
                "accession": "VISTA Enhancer Browser loci (LBNL)",
                "version": "as cached in data/knowledge/vista",
            },
            {"accession": "GEO GSE126550 (Kircher et al. 2019)", "version": "as cached in data/knowledge"},
            {
                "accession": "this repository, the compiled element runs " + ", ".join(RUNS),
                "version": "pinned by sha256",
            },
        ],
        "inputs": inputs(),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "reciprocal_overlap": measured.RECIPROCAL_OVERLAP,
            "active_threshold_log2": measured.MPRA_ACTIVE,
            "chromosomes": CHROMS,
        },
        "exclusions": [],
        "partitions": "n/a: a census of summaries, no evaluation",
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args()
    payload = {
        "review_item": "R6",
        "aggregations": AGGREGATIONS,
        "files": file_census(),
        "elements": element_census(),
        "episomal_wording": EPISOMAL,
        "reading": (
            "one summary decides a label: the lentiMPRA maximum (A1). It is taken over more than one tile "
            "for the elements counted under elements.lentimpra, and where the tiles disagree the maximum "
            "has so far always read active. The other maxima are descriptive or already sit beside the "
            "observations they summarise"
        ),
    }
    lt = payload["elements"]["lentimpra"]
    print(
        f"lentiMPRA: {lt.get('elements', 0)} elements, "
        f"{lt.get('cells_over_more_than_one_tile', 0)} cell readings over more than one tile, "
        f"{lt.get('cells_where_tiles_disagree', 0)} where the tiles disagree"
    )
    print(
        "vista",
        payload["elements"]["vista"],
        "satmut",
        {k: v for k, v in payload["elements"]["satmut"].items() if k != "loci_with_repeat_experiments"},
    )
    if not args.no_save:
        print(f"saved {save_result('measured_aggregation_census', payload, manifest=manifest())}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
