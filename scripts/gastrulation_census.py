# SPDX-License-Identifier: AGPL-3.0-or-later
"""Gastrulation against a measured census (pre-registered 2026-09-28).

    uv run python scripts/gastrulation_census.py distill   # cache files -> gastrulation_census_counts
    uv run python scripts/gastrulation_census.py compare   # counts + one default model run -> comparison

`distill` reads two public per-cell annotation tables fetched into the git-ignored data/cache/gastrulation/:

    E-MTAB-9388.sdrf.txt          Tyser et al. 2021 human CS7 gastrula, 1,195 cells (1.3 MB)
        https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/388/E-MTAB-9388/Files/E-MTAB-9388.sdrf.txt
    pijuan_sala_2019_meta.tab.gz  Pijuan-Sala et al. 2019 mouse atlas cell metadata (4.9 MB)
        https://content.cruk.cam.ac.uk/jmlab/atlas_data/meta.tab.gz

and keeps only cell counts per label (and per dissection site for the human embryo). `compare` applies
the registered mappings and tolerance in genomeos/runtime/gastrulation.py, runs the model once with
CENSUS_MODEL_RUN (the command-line defaults), and writes the verdict. Nothing is tuned.
"""

from __future__ import annotations

import csv
import gzip
import json
import sys
from collections import Counter
from pathlib import Path

from genomeos import manifest as mf
from genomeos.results import RESULTS_DIR, save_result
from genomeos.runtime import gastrulation as g

CACHE = Path("data/cache/gastrulation")
HUMAN = CACHE / "E-MTAB-9388.sdrf.txt"
MOUSE = CACHE / "pijuan_sala_2019_meta.tab.gz"
COUNTS = "gastrulation_census_counts"
COMPARISON = "gastrulation_census_comparison"
LAYERS = ("ectoderm", "mesoderm", "endoderm")

SOURCES = [
    {
        "accession": "E-MTAB-9388",
        "version": "ArrayExpress SDRF, public release 2021-08-23",
        "paper": "Tyser RCV et al. 2021, Nature 600:285-289, doi:10.1038/s41586-021-04158-y (PMC7615353)",
        "counts": "1,195 Smart-seq2 cells after QC from one CS7 human embryo (16-19 days), FACS-isolated "
        "from three dissected regions (665 caudal, 340 rostral, 190 yolk sac); author cluster labels",
        "licence": "EMBL-EBI ArrayExpress terms of use (public, no restriction stated); no personal data "
        "beyond stage and sex is in the SDRF and none is kept",
    },
    {
        "accession": "E-MTAB-6967",
        "version": "meta.tab.gz, content.cruk.cam.ac.uk/jmlab/atlas_data (MouseGastrulationData 1.x source)",
        "paper": "Pijuan-Sala B et al. 2019, Nature 566:490-495, doi:10.1038/s41586-019-0933-9 (PMC6522369)",
        "counts": "116,312 10x cells from pooled whole C57BL/6 embryos at nine time points E6.5-E8.5; "
        "author celltype labels; doublets and stripped nuclei carry no label and are dropped",
        "licence": "ArrayExpress terms of use; the Bioconductor wrapper MouseGastrulationData is GPL-3",
    },
]


def _human() -> dict:
    rows = list(csv.reader(HUMAN.open(), delimiter="\t"))
    h = rows[0]
    i_id = h.index("Source Name")
    i_lab = h.index("Characteristics[inferred cell type - authors labels]")
    i_site = h.index("Characteristics[sampling site]")
    cells = {r[i_id]: (r[i_lab], r[i_site]) for r in rows[1:]}  # two rows per cell (paired reads)
    by: dict[str, Counter] = {}
    for lab, site in cells.values():
        by.setdefault(lab, Counter())[site] += 1
    return {
        "cells": len(cells),
        "by_label_and_site": {k: dict(v) for k, v in sorted(by.items())},
    }


def _mouse() -> dict:
    by: dict[str, Counter] = {}
    dropped = Counter()
    with gzip.open(MOUSE, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["celltype"] == "NA":
                dropped["doublet" if row["doublet"] == "TRUE" else "stripped"] += 1
                continue
            by.setdefault(row["stage"], Counter())[row["celltype"]] += 1
    return {
        "cells": sum(sum(c.values()) for c in by.values()),
        "unlabelled_dropped": dict(dropped),
        "by_stage_and_label": {s: dict(c.most_common()) for s, c in sorted(by.items())},
    }


def distill() -> None:
    payload = {
        "question": "measured cell-type counts in a human and a mouse gastrula, the census "
        "genomeos/runtime/gastrulation.py is compared with",
        "human_cs7": _human(),
        "mouse_atlas": _mouse(),
        "sampling_notes": [
            "Both count sampled cells, not cells in the embryo. Human: one embryo, cells FACS-sorted into "
            "plates per region, so the region shares (665/340/190) reflect how many were sorted, not tissue "
            "size. Mouse: pooled whole embryos per sample, 10x capture; per-stage totals vary 2,075-16,909.",
            "Human 'endodermal cell' holds hypoblast and yolk-sac endoderm (extraembryonic) with definitive "
            "endoderm; the subclusters are not in the SDRF. Human 'ectodermal cell' holds amnion with "
            "non-neural ectoderm. Primordial germ cells sit inside the human 'primitive streak' cluster.",
        ],
    }
    manifest = {
        "sources": SOURCES,
        "inputs": [mf.input_entry(HUMAN), mf.input_entry(MOUSE)],
        "assembly": "n/a: cell annotations, no genome coordinates",
        "coordinates": "n/a: no genomic intervals",
        "parameters": {
            "human_label_column": "inferred cell type - authors labels",
            "mouse_label_column": "celltype",
        },
        "exclusions": [
            "mouse cells with celltype NA (doublets, stripped nuclei), counted in unlabelled_dropped"
        ],
        "partitions": "n/a: no evaluation split",
    }
    print(save_result(COUNTS, payload, manifest=manifest))


def _shares(counts: dict[str, int], mapping: dict[str, tuple[str, ...]]) -> dict:
    n = {layer: sum(counts.get(lab, 0) for lab in mapping[layer]) for layer in LAYERS}
    tot = sum(n.values())
    return {"cells": n, "total": tot, "shares": {k: round(v / tot, 4) for k, v in n.items()}}


def compare() -> None:
    counts = json.loads((RESULTS_DIR / f"{COUNTS}.json").read_text())
    hum = counts["human_cs7"]["by_label_and_site"]
    human = {}
    for mname, mapping in g.CENSUS_HUMAN_MAPPINGS.items():
        for sname, sites in g.CENSUS_HUMAN_SITES.items():
            c = {lab: sum(v for s, v in by.items() if s in sites) for lab, by in hum.items()}
            human[f"{mname}/{sname}"] = _shares(c, mapping)
    mouse = {}
    for stage in g.CENSUS_MOUSE_STAGES:
        c = counts["mouse_atlas"]["by_stage_and_label"][stage]
        for mname, mapping in g.CENSUS_MOUSE_MAPPINGS.items():
            mouse[f"{stage}/{mname}"] = _shares(c, mapping)

    r = g.run_gastrulation(**g.CENSUS_MODEL_RUN)
    model = {k: round(v, 4) for k, v in r.proportions().items()}

    def judge(variants: dict) -> dict:
        out = {}
        for layer in LAYERS:
            vals = [v["shares"][layer] for v in variants.values()]
            lo, hi = min(vals), max(vals)
            m = model[layer]
            gap = round(max(lo - m, m - hi, 0.0), 4)
            out[layer] = {
                "interval": [lo, hi],
                "model": m,
                "outside_by": gap,
                "pass": gap <= g.CENSUS_TOLERANCE,
            }
        return out

    human_j = judge(human)
    verdict = "not contradicted" if all(v["pass"] for v in human_j.values()) else "falsified"
    payload = {
        "question": "do the gastrulation model's default germ-layer proportions fall within the registered "
        "tolerance of a measured human CS7 gastrula census",
        "registration": "docs/DESIGN-MINIMAL-CELL.md, 'Gastrulation against a measured census (registered "
        "2026-09-28)'; constants in genomeos/runtime/gastrulation.py",
        "tolerance": g.CENSUS_TOLERANCE,
        "model_run": g.CENSUS_MODEL_RUN,
        "model": model,
        "human_cs7_variants": human,
        "human_cs7_judgement": human_j,
        "verdict": verdict,
        "mouse_variants_no_verdict": mouse,
        "mouse_judgement_no_verdict": judge(mouse),
    }
    manifest = {
        "sources": SOURCES,
        "inputs": [
            mf.input_entry(RESULTS_DIR / f"{COUNTS}.json"),
            mf.input_entry("data/demo/gastrulation.bio"),
            mf.input_entry("genomeos/runtime/gastrulation.py"),
        ],
        "assembly": "n/a: cell annotations and a simulated network",
        "coordinates": "n/a: no genomic intervals",
        "parameters": {"tolerance": g.CENSUS_TOLERANCE, **g.CENSUS_MODEL_RUN},
        "exclusions": [
            f"human labels never mapped: {', '.join(g.CENSUS_HUMAN_EXTRAEMBRYONIC)}; "
            "'differentiated' also leaves out epiblast cell and primitive streak",
            "mouse extraembryonic, blood, PGC, caudal epiblast and NMP labels are never mapped",
        ],
        "partitions": "n/a: one pre-registered comparison, nothing fitted",
    }
    print(save_result(COMPARISON, payload, manifest=manifest))
    print(json.dumps({"model": model, "human": human_j, "verdict": verdict}, indent=1))


if __name__ == "__main__":
    {"distill": distill, "compare": compare}[sys.argv[1]]()
