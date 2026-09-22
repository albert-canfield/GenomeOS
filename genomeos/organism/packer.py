# SPDX-License-Identifier: AGPL-3.0-or-later
"""Packer et al. 2019 (Science, GEO GSE126954): 86,024 embryonic single cells with their lineage
and terminal cell type. Streamed once and distilled to one table: lineage → cell type, plus an
agreement score against the reference lineage's tissue per cell. Two independent sources
(Sulston's microscope, Packer's transcriptomes) checking each other."""

from __future__ import annotations

import csv
import gzip
import io
import urllib.request
from collections import Counter

from .reference import ReferenceLineage

GEO_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126954/suppl/GSE126954_cell_annotation.csv.gz"

# Packer cell.type -> WormWeb tissue label
PACKER_TISSUE = {
    "Body_wall_muscle": "muscle",
    "Pharyngeal_muscle": "muscle",
    "Intestinal_and_rectal_muscle": "muscle",
    "Ciliated_amphid_neuron": "neuron",
    "Ciliated_non_amphid_neuron": "neuron",
    "Pharyngeal_neuron": "neuron",
    "Hypodermis": "hypoderm",
    "Seam_cell": "hypoderm",
    "Glia": "glia",
    "Intestine": "intestine",
    "Germline": "repro",
    "Z1_Z4": "repro",
    "Coelomocyte": "coelomocyte",
    "Pharyngeal_marginal_cell": "marginal",
    "Pharyngeal_gland": "gland",
    "Excretory_gland": "gland",
    "Rectal_gland": "gland",
    "Pharyngeal_intestinal_valve": "valve",
    "Rectal_cell": "rectal",
    "Excretory_cell": "excretory",
    "Excretory_duct_and_pore": "excretory",
    "Arcade_cell": "epithelium",
    "M_cell": "mesoderm",
    "GLR": "other",
    "hmc": "other",
    "hmc_homolog": "other",
    "hmc_and_homolog": "other",
}
GLIA = {"socket", "sheath"}


def expand_lineage(lin: str) -> list[str]:
    """`ABpxaapapap` -> ABplaapapap, ABpraapapap; `A/B` alternatives are split."""
    out: list[str] = []
    for alt in lin.split("/"):
        names = [""]
        for ch in alt:
            names = [n + c for n in names for c in (("l", "r") if ch == "x" else (ch,))]
        out.extend(names)
    return out


def stream_annotation(url: str = GEO_URL, timeout: float = 120.0) -> list[dict]:
    with urllib.request.urlopen(url, timeout=timeout) as r:  # noqa: S310 (fixed public URL)
        raw = r.read()
    text = gzip.decompress(raw).decode()
    rows = list(csv.DictReader(io.StringIO(text)))
    return [
        r for r in rows if r.get("lineage") not in (None, "", "NA") and r.get("cell.type") not in ("", "NA")
    ]


def _reference_tissue(ref: ReferenceLineage, cell_id: str) -> str:
    """Majority tissue among the terminal, surviving descendants of a reference cell."""
    if cell_id not in ref.cells:
        return ""
    stack, counts = [cell_id], Counter()
    while stack:
        c = ref.cells[stack.pop()]
        if c.children:
            stack.extend(c.children)
        elif c.dies is None:
            counts[c.tissue] += 1
    if not counts:
        return ""
    t = counts.most_common(1)[0][0]
    return "glia" if t in GLIA else t


def _single_tissue(ref: ReferenceLineage, cell_id: str) -> bool:
    """True when every surviving terminal descendant of the reference cell has the same tissue."""
    stack, seen = [cell_id], set()
    while stack:
        c = ref.cells[stack.pop()]
        if c.children:
            stack.extend(c.children)
        elif c.dies is None:
            seen.add("glia" if c.tissue in GLIA else c.tissue)
            if len(seen) > 1:
                return False
    return len(seen) == 1


def distil(rows: list[dict], ref: ReferenceLineage) -> dict:
    by_lineage: dict[str, Counter] = {}
    for r in rows:
        by_lineage.setdefault(r["lineage"], Counter())[r["cell.type"]] += 1
    table = {}
    agree = disagree = unchecked = 0
    pure_agree = pure_total = (
        0  # ids whose reference descendants are one tissue: labelling depth cannot confound
    )
    disagreements = []
    for lin, ctr in sorted(by_lineage.items()):
        cell_type, n = ctr.most_common(1)[0]
        packer_tissue = PACKER_TISSUE.get(cell_type, "")
        ref_tissues = {t for t in (_reference_tissue(ref, x) for x in expand_lineage(lin)) if t}
        row = {"cell_type": cell_type, "cells": sum(ctr.values()), "majority": n}
        if packer_tissue and ref_tissues:
            ok = packer_tissue in ref_tissues
            row["reference_tissue"] = sorted(ref_tissues)
            row["agree"] = ok
            agree += ok
            disagree += not ok
            pure = all(_single_tissue(ref, x) for x in expand_lineage(lin) if x in ref.cells)
            row["single_tissue_descendants"] = pure
            if pure:
                pure_total += 1
                pure_agree += ok
            if not ok:
                disagreements.append({"lineage": lin, "packer": cell_type, "reference": sorted(ref_tissues)})
        else:
            unchecked += 1
        table[lin] = row
    checked = agree + disagree
    return {
        "source": "Packer et al. 2019, Science 365:eaax1971; GEO GSE126954 cell annotation",
        "url": GEO_URL,
        "cells_with_lineage_and_type": len(rows),
        "lineage_ids": len(by_lineage),
        "cell_types": len({r["cell.type"] for r in rows}),
        "agreement_with_reference": {
            "checked": checked,
            "agree": agree,
            "rate": round(agree / checked, 4) if checked else None,
            "unchecked": unchecked,
            "single_tissue_ids": pure_total,
            "single_tissue_agree": pure_agree,
            "single_tissue_rate": round(pure_agree / pure_total, 4) if pure_total else None,
            "note": "internal ids whose descendants span several tissues are a labelling-depth question, "
            "not a disagreement between the atlas and the lineage; the single-tissue rate excludes them",
        },
        "disagreements": disagreements,
        "by_lineage": table,
    }
