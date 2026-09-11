# SPDX-License-Identifier: AGPL-3.0-or-later
"""Ma et al. 2021 (Nat Methods 18:893; Zenodo 10.5281/zenodo.4737593, CC BY 4.0): a 4D protein atlas of
266 transcription factors in ~1,200 lineaged C. elegans embryonic cells, 1.25 min per frame to the bean
stage. This is the measured reader configuration of the embryo: which factors each cell carries.

Streamed once (88 MB) and distilled to a per-cell presence table. Expression is a reporter intensity with
its own scale per factor, so presence is called per factor: a cell expresses a factor when its maximum
adjusted expression reaches PRESENCE_FRACTION of that factor's maximum over all cells. The distilled table
generates data/organisms/celegans/reader.bio, one `express` decision per cell, so every cell of the
organism program carries its measured factors in its context; and a textbook check reports how well
known factors predict the tissues their cells go on to form.
"""

from __future__ import annotations

import csv
import io
import json
import urllib.request
import zipfile
from collections import defaultdict
from pathlib import Path

from .reference import ReferenceLineage

ZENODO_URL = "https://zenodo.org/api/records/4737593/files/Time_series_expression_dataset.zip/content"
SOURCE = "Ma et al. 2021, Nat Methods 18:893 (Zenodo 4737593)"
PRESENCE_FRACTION = 0.2  # of the factor's maximum adjusted expression over all cells
FRAME_MIN = 1.25  # minutes per frame
CELLS_FILE = Path("data/results/celegans_tf_atlas_cells.json")

# Textbook factor -> tissue relations to check the atlas and the lineage against each other.
# WormWeb tissue labels; "pharynx" is checked against Packer 2019 cell types instead (WormWeb has no class).
TEXTBOOK = [
    ("ELT-2", "intestine", "Fukushige et al. 1998, Dev Biol 198:286"),
    ("END-1", "intestine", "Zhu et al. 1997, Genes Dev 11:2883"),
    ("ELT-7", "intestine", "Sommermann et al. 2010, Dev Biol 347:154"),
    ("HLH-1", "muscle", "Krause et al. 1990, Cell 63:907"),
    ("UNC-120", "muscle", "Fukushige et al. 2006, Genes Dev 20:3395"),
    ("ELT-1", "hypoderm", "Page et al. 1997, Genes Dev 11:1651"),
    ("LIN-26", "hypoderm", "Labouesse et al. 1994, Development 120:2359"),
    ("NHR-25", "hypoderm", "Gissendanner & Sluder 2000, Dev Biol 221:259"),
    ("PHA-4", "pharynx", "Mango et al. 1994, Development 120:3019"),
    ("CEH-22", "pharynx", "Okkema & Fire 1994, Development 120:2175"),
]
PACKER_PHARYNX = (
    "Pharyngeal_muscle",
    "Pharyngeal_neuron",
    "Pharyngeal_marginal_cell",
    "Pharyngeal_gland",
    "Pharyngeal_intestinal_valve",
    "Arcade_cell",
)


def fetch(url: str = ZENODO_URL, timeout: float = 900.0) -> bytes:
    with urllib.request.urlopen(url, timeout=timeout) as r:  # noqa: S310 (fixed public URL)
        return r.read()


def distil(archive: bytes, fraction: float = PRESENCE_FRACTION) -> dict:
    """Per factor: the cells that express it (max adjusted expression, first frame at or above threshold)."""
    z = zipfile.ZipFile(io.BytesIO(archive))
    peak: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))  # tf -> cell -> max
    first: dict[str, dict[str, int]] = defaultdict(dict)  # tf -> cell -> first frame with expression
    strains: dict[str, list[str]] = defaultdict(list)
    seen_first: dict[str, int] = {}  # cell -> first frame it is tracked in any strain (its birth)
    seen_last: dict[str, int] = {}  # cell -> last frame (its division or the end of imaging)
    for name in z.namelist():
        if not name.endswith(".csv"):
            continue
        tf = Path(name).stem.split("_")[0]
        strains[tf].append(Path(name).stem)
        for row in csv.DictReader(io.StringIO(z.read(name).decode("utf-8", "replace"))):
            cell = row["Cell-name"]
            try:
                v = float(row["Adjustment-expression"])
                t = int(row["Time"])
            except ValueError:
                continue
            if v > peak[tf][cell]:
                peak[tf][cell] = v
            if v > 0 and (cell not in first[tf] or t < first[tf][cell]):
                first[tf][cell] = t
            if t < seen_first.get(cell, 10**9):
                seen_first[cell] = t
            if t > seen_last.get(cell, -1):
                seen_last[cell] = t
    table: dict[str, dict] = {}
    for tf, cells in peak.items():
        top = max(cells.values(), default=0.0)
        thr = top * fraction
        expressing = {c: [round(v, 1), first[tf].get(c, 0)] for c, v in cells.items() if top > 0 and v >= thr}
        table[tf] = {
            "max": round(top, 1),
            "threshold": round(thr, 1),
            "strains": strains[tf],
            "cells": expressing,
        }
    table["_lifetimes"] = {c: [seen_first[c], seen_last[c]] for c in sorted(seen_first)}
    return table


def cells_expressing(table: dict) -> dict[str, list[str]]:
    """cell -> factors present, from the per-factor table."""
    out: dict[str, set[str]] = defaultdict(set)
    for tf, rec in table.items():
        if tf.startswith("_"):
            continue
        for cell in rec["cells"]:
            out[cell].add(tf)
    return {c: sorted(v) for c, v in sorted(out.items())}


def _descendant_tissues(ref: ReferenceLineage, cell: str) -> list[str]:
    """Tissues of the surviving terminal descendants of a reference cell (the cell itself if terminal)."""
    if cell not in ref.cells:
        return []
    out, stack = [], [cell]
    while stack:
        c = ref.cells[stack.pop()]
        if c.children:
            stack.extend(c.children)
        elif c.dies is None and c.tissue:
            out.append(c.tissue)
    return out


def textbook_check(table: dict, ref: ReferenceLineage, packer: dict | None) -> list[dict]:
    """For each textbook factor: of the terminal descendants of expressing cells, the share in the expected
    tissue (precision), and of all cells of that tissue, the share descending from an expressing cell
    (recall). Pharynx uses Packer 2019 lineage ids."""
    checks = []
    packer_lineages = (packer or {}).get("by_lineage", {})
    for tf, tissue, source in TEXTBOOK:
        rec = table.get(tf)
        if rec is None:
            checks.append({"factor": tf, "tissue": tissue, "source": source, "status": "factor not in atlas"})
            continue
        expressing = list(rec["cells"])
        if tissue == "pharynx":
            hits = total = 0
            for cell in expressing:
                for lin, row in packer_lineages.items():
                    if any(
                        alt.startswith(cell) or cell.startswith(alt.replace("x", "l"))
                        for alt in lin.split("/")
                    ):
                        total += 1
                        hits += row["cell_type"] in PACKER_PHARYNX
            checks.append(
                {
                    "factor": tf,
                    "tissue": tissue,
                    "source": source,
                    "expressing_cells": len(expressing),
                    "checked_against": "Packer 2019 lineage ids",
                    "descendants_checked": total,
                    "precision": round(hits / total, 3) if total else None,
                }
            )
            continue
        tissues = [t for cell in expressing for t in _descendant_tissues(ref, cell)]
        hits = sum(1 for t in tissues if t == tissue)
        tissue_total = sum(1 for c in ref.terminal() if c.tissue == tissue)
        covered = {
            c.id for cell in expressing for c in _terminal_descendants(ref, cell) if c.tissue == tissue
        }
        checks.append(
            {
                "factor": tf,
                "tissue": tissue,
                "source": source,
                "expressing_cells": len(expressing),
                "descendants_checked": len(tissues),
                "precision": round(hits / len(tissues), 3) if tissues else None,
                "recall": round(len(covered) / tissue_total, 3) if tissue_total else None,
            }
        )
    return checks


def _terminal_descendants(ref: ReferenceLineage, cell: str):
    if cell not in ref.cells:
        return []
    out, stack = [], [cell]
    while stack:
        c = ref.cells[stack.pop()]
        if c.children:
            stack.extend(c.children)
        elif c.dies is None:
            out.append(c)
    return out


def to_bio_reader(table: dict, module_name: str = "organism.celegans.reader") -> str:
    """One `express` decision per cell: the measured factors it carries."""
    by_cell = cells_expressing(table)
    lines = [
        "# Generated by `genomeos data distil --only celegans_tf_atlas` from the Ma et al. 2021 atlas.",
        f"# Presence: max adjusted expression at or above {PRESENCE_FRACTION:.0%} of the factor's maximum.",
        "# Every cell named here carries these measured factors in its context; nothing is inferred.",
        f"module {module_name}",
        "",
    ]
    for cell, tfs in by_cell.items():
        lines.append(
            f"decision read_{cell} {{ action: express; when: cell = {cell}; sets: {', '.join(tfs)}; "
            f'evidence: experimental "{SOURCE}"; confidence: 0.8 }}'
        )
    return "\n".join(lines) + "\n"


def summary(table: dict, checks: list[dict]) -> dict:
    real = {tf: rec for tf, rec in table.items() if not tf.startswith("_")}
    by_cell = cells_expressing(table)
    per_tf = {tf: len(rec["cells"]) for tf, rec in real.items()}
    firsts = [FRAME_MIN * min(v[1] for v in rec["cells"].values()) for rec in real.values() if rec["cells"]]
    return {
        "source": SOURCE,
        "url": ZENODO_URL,
        "licence": "CC BY 4.0",
        "presence_fraction_of_max": PRESENCE_FRACTION,
        "factors": len(real),
        "strains": sum(len(rec["strains"]) for rec in real.values()),
        "cells_with_a_factor": len(by_cell),
        "factors_per_cell_median": sorted(len(v) for v in by_cell.values())[len(by_cell) // 2]
        if by_cell
        else 0,
        "cells_per_factor_median": sorted(per_tf.values())[len(per_tf) // 2] if per_tf else 0,
        "earliest_expression_min": round(min(firsts), 1) if firsts else None,
        "textbook": checks,
        "cells_file": str(CELLS_FILE),
    }


def save_cells(table: dict, path: Path = CELLS_FILE) -> Path:
    """Compact per-cell table: factor names once, then for each cell the indices of the factors it carries
    and the frame each was first seen (about a fifth of the size of a per-factor listing)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    factors = sorted(tf for tf in table if not tf.startswith("_"))
    index = {tf: i for i, tf in enumerate(factors)}
    cells: dict[str, list[list[int]]] = defaultdict(list)
    for tf, rec in table.items():
        if tf.startswith("_"):
            continue
        for cell, (_, first_frame) in rec["cells"].items():
            cells[cell].append([index[tf], int(first_frame)])
    compact = {
        "source": SOURCE,
        "presence_fraction_of_max": PRESENCE_FRACTION,
        "frame_min": FRAME_MIN,
        "factors": factors,
        "thresholds": [table[tf]["threshold"] for tf in factors],
        "cells": {c: sorted(v) for c, v in sorted(cells.items())},
        "lifetimes": table.get("_lifetimes", {}),
    }
    path.write_text(json.dumps(compact, separators=(",", ":")))
    return path


def load_cells(path: Path = CELLS_FILE) -> dict[str, list[str]]:
    """cell -> factor names, from the compact table."""
    data = json.loads(Path(path).read_text())
    names = data["factors"]
    return {cell: [names[i] for i, _ in pairs] for cell, pairs in data["cells"].items()}
