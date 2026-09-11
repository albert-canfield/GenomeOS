# SPDX-License-Identifier: AGPL-3.0-or-later
"""The reference lineage of C. elegans: the ground truth a grown organism is scored against.

Source: the complete timed cell lineage compiled by Nikhil Bhatla (WormWeb,
CC BY 1.0) from Sulston et al. 1983 (embryo) and Sulston & Horvitz 1977
(post-embryonic): 2,183 named cells from P0 to the adult hermaphrodite,
with birth and division times, 131 programmed deaths, and the tissue of
every terminal cell. It is streamed once, distilled to a compact JSON under
data/results/ (committed) and the raw file discarded (stream, distil, discard).

The time axis is the one given by the source, in minutes from first
cleavage. Its absolute scale has not been re-verified against the 20 °C
chart of Sulston 1983, so timing comparisons are made on this axis and said
to be so; topology (who divides into whom), fates and deaths are exact.
"""

from __future__ import annotations

import json
import re
import statistics
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

WORMWEB_URL = "http://wormweb.org/js/json-celllineage.js"
KNOWLEDGE = Path("data/results/celegans_lineage_cells.json")
EMBRYONIC_FOUNDERS = ("AB", "MS", "E", "C", "D", "P4")
HATCH_MIN = 800.0  # cells born after this on the source axis are post-embryonic

# WormWeb tissue labels -> GenomeOS cell types (data/organisms/celegans/cell_types.bio)
TISSUE_CELL_TYPE = {
    "neuron": "Neuron",
    "hypoderm": "Hypodermis",
    "muscle": "Muscle",
    "repro": "Reproductive",
    "intestine": "Intestine",
    "epithelium": "Epithelium",
    "socket": "GliaSocket",
    "sheath": "GliaSheath",
    "marginal": "PharyngealMarginal",
    "valve": "Valve",
    "coelomocyte": "Coelomocyte",
    "excretory": "Excretory",
    "gland": "Gland",
    "rectal": "Rectal",
    "tail": "TailSpike",
    "mesoderm": "Mesoderm",
    "other": "Other",
}


@dataclass(slots=True)
class RefCell:
    id: str  # lineage name (ABalaaaalal); unique
    name: str  # display name: same, or the terminal name (AINL) or blast name (P1)
    parent: str | None
    born: float
    divides: float | None = None
    dies: float | None = None
    tissue: str = ""  # WormWeb tissue label for terminal cells
    founder: str = ""  # AB, MS, E, C, D, P (germline P1..P4, Z2/Z3) or a blast class (V, Pn, H, ...)
    generation: int = 0  # divisions since the founder
    children: list[str] = field(default_factory=list)

    @property
    def terminal(self) -> bool:
        return not self.children

    @property
    def end(self) -> float:
        if self.dies is not None:
            return self.dies
        return self.divides if self.divides is not None else float("inf")

    @property
    def cell_type(self) -> str:
        return TISSUE_CELL_TYPE.get(self.tissue, "UNKNOWN")


def fetch_wormweb(url: str = WORMWEB_URL, timeout: float = 60.0) -> str:
    with urllib.request.urlopen(url, timeout=timeout) as r:  # noqa: S310 (fixed public URL)
        return r.read().decode("utf-8", "replace")


def parse_wormweb(js: str) -> dict:
    """The file is `function getJson() { var json; json = {...}; return json; }`: extract the literal."""
    m = re.search(r"=\s*\{", js)
    if not m:
        raise ValueError("no object literal found in the WormWeb file")
    start = m.end() - 1
    depth = 0
    end = start
    for i, ch in enumerate(js[start:], start):
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                end = i + 1
                break
    obj = js[start:end]
    obj = re.sub(r"(\{|,)\s*([A-Za-z_][A-Za-z0-9_]*)\s*:", r'\1"\2":', obj)  # quote keys
    obj = re.sub(r",\s*([\]}])", r"\1", obj)  # trailing commas
    return json.loads(obj)


def _blast_class(name: str) -> str:
    """V1L -> V, P1 -> Pn, TL -> T, Md -> M, QR -> Q, Z1 -> Z."""
    m = re.match(r"^([A-Z]+)", name)
    base = m.group(1) if m else name
    if len(base) > 1 and base[-1] in "LR":
        base = base[:-1]
    return "Pn" if base == "P" else base


def distil(tree: dict) -> dict:
    """WormWeb tree -> compact cell table with founders and generations."""
    cells: dict[str, RefCell] = {}
    order: list[str] = []

    def walk(node: dict, parent: RefCell | None) -> None:
        d = node["data"]
        cid = node.get("did") or node["name"]
        if cid == "P0a":
            cid = "AB"  # the source keys AB by its P0 daughter position
        born = float(d["totalDistance"] - d["levelDistance"])
        c = RefCell(cid, node["name"], parent.id if parent else None, born)
        kids = node.get("children", [])
        if kids:
            c.divides = float(d["totalDistance"])
        elif d.get("deathDistance", 0):
            c.dies = born + float(d["deathDistance"])
        else:
            c.tissue = d.get("type", "")
        if parent:
            parent.children.append(cid)
        cells[cid] = c
        order.append(cid)
        for k in kids:
            walk(k, c)

    walk(tree, None)
    # founders and generations: embryonic founders by name, post-embryonic blast lineages by class
    for cid in order:
        c = cells[cid]
        p = cells[c.parent] if c.parent else None
        if c.name in EMBRYONIC_FOUNDERS and (p is None or p.founder in ("P", "EMS")):
            c.founder, c.generation = c.name, 0
        elif c.name in ("P0", "P1", "P2", "P3") and (p is None or p.founder == "P"):
            c.founder, c.generation = "P", 0
        elif c.name == "EMS":
            c.founder, c.generation = "EMS", 0
        elif p is not None and c.children and c.born < HATCH_MIN <= (c.divides or 0):
            c.founder, c.generation = _blast_class(c.name), 0  # a larval blast cell
        elif p is not None:
            c.founder, c.generation = p.founder, p.generation + 1
    return {
        "source": "WormWeb C. elegans interactive cell lineage (Nikhil Bhatla), from Sulston et al. 1983 "
        "and Sulston & Horvitz 1977",
        "url": WORMWEB_URL,
        "licence": "CC BY 1.0",
        "time_axis": "minutes from first cleavage as given by the source; absolute scale not re-verified",
        "cells": [
            {
                "id": c.id,
                "name": c.name,
                "parent": c.parent,
                "born": c.born,
                "divides": c.divides,
                "dies": c.dies,
                "tissue": c.tissue,
                "founder": c.founder,
                "generation": c.generation,
            }
            for c in (cells[i] for i in order)
        ],
    }


@dataclass(slots=True)
class ReferenceLineage:
    cells: dict[str, RefCell]
    meta: dict

    @classmethod
    def from_dict(cls, data: dict) -> ReferenceLineage:
        cells: dict[str, RefCell] = {}
        for r in data["cells"]:
            cells[r["id"]] = RefCell(
                r["id"],
                r["name"],
                r["parent"],
                r["born"],
                r["divides"],
                r["dies"],
                r["tissue"],
                r["founder"],
                r["generation"],
            )
        for c in cells.values():
            if c.parent:
                cells[c.parent].children.append(c.id)
        return cls(cells, {k: v for k, v in data.items() if k != "cells"})

    @classmethod
    def load(cls, path: str | Path = KNOWLEDGE) -> ReferenceLineage:
        return cls.from_dict(json.loads(Path(path).read_text()))

    # ---- queries ---------------------------------------------------------

    def alive_at(self, t: float) -> list[RefCell]:
        return [c for c in self.cells.values() if c.born <= t < c.end]

    def count_at(self, t: float) -> int:
        return len(self.alive_at(t))

    def deaths(self) -> list[RefCell]:
        return [c for c in self.cells.values() if c.dies is not None]

    def terminal(self) -> list[RefCell]:
        return [c for c in self.cells.values() if c.terminal and c.dies is None]

    def embryonic(self) -> list[RefCell]:
        return [c for c in self.cells.values() if c.born < HATCH_MIN]

    def fates(self) -> dict[str, int]:
        out: dict[str, int] = {}
        for c in self.terminal():
            out[c.cell_type] = out.get(c.cell_type, 0) + 1
        return dict(sorted(out.items(), key=lambda kv: -kv[1]))

    def cycle_stats(self) -> dict[tuple[str, int], dict]:
        """Cell-cycle length (birth to division) per (founder, generation): mean, sd, n."""
        groups: dict[tuple[str, int], list[float]] = {}
        for c in self.cells.values():
            if c.divides is not None and c.parent is not None:
                groups.setdefault((c.founder, c.generation), []).append(c.divides - c.born)
        out = {}
        for key, xs in sorted(groups.items()):
            out[key] = {
                "mean": round(statistics.fmean(xs), 1),
                "sd": round(statistics.pstdev(xs), 1) if len(xs) > 1 else 0.0,
                "n": len(xs),
            }
        return out

    def summary(self) -> dict:
        term = self.terminal()
        adult_end = max(c.end for c in self.cells.values() if c.end != float("inf"))
        return {
            "cells": len(self.cells),
            "terminal": len(term),
            "deaths": len(self.deaths()),
            "embryonic_deaths": sum(1 for c in self.deaths() if c.dies is not None and c.dies < HATCH_MIN),
            "last_event_min": adult_end,
            "alive_at": {str(t): self.count_at(t) for t in (0.5, 50, 100, 200, 350, 500, 800, 2000, 4000)},
            "alive_adult": self.count_at(adult_end + 1),
            "fates": self.fates(),
            "founders": sorted({c.founder for c in self.cells.values()}),
        }

    # ---- BioLang generation ----------------------------------------------

    def to_bio_timers(self) -> str:
        """Cell-cycle timers per founder and generation, distilled from the tree with their spread."""
        lines = [
            "# Generated by `genomeos data distil --only celegans_lineage` from the reference lineage.",
            "# One timer per founder lineage and generation: mean cycle length and its spread (sd).",
            "module organism.celegans.timers",
            "",
        ]
        for (founder, gen), s in self.cycle_stats().items():
            conf = 0.9 if s["n"] >= 8 else 0.7 if s["n"] >= 3 else 0.5
            lines.append(
                f"timer cycle_{founder}_{gen} {{ duration: {s['mean']} min; sd: {s['sd']}; "
                f"when: lineage = {founder}, generation = {gen}; "
                f'evidence: experimental "Sulston et al. 1983 via WormWeb lineage (n={s["n"]})"; '
                f"confidence: {conf} }}"
            )
        return "\n".join(lines) + "\n"

    def to_bio_program(self, since: float = 0.0, until: float = float("inf")) -> str:
        """The observed lineage as decisions: every division with its daughters' names, every death,
        every terminal fate. This is the lineage-lookup layer; mechanistic modules imported before it
        take precedence for the cells they cover."""
        ev = 'evidence: experimental "Sulston 1983 via WormWeb"'
        lines = [
            "# Generated by `genomeos data distil --only celegans_lineage` from the reference lineage.",
            "# Divisions (with daughter names), programmed deaths and terminal fates as observed.",
            "module organism.celegans.lineage",
            "import cell_types.bio",
            "",
        ]
        for c in self.cells.values():
            if not since <= c.born < until:
                continue
            if c.children:
                kids = [self.cells[k] for k in c.children[:2]]
                a, b = (k.id for k in kids)
                founders = ", ".join(f"{k.id} = {k.founder}" for k in kids if k.generation == 0)
                lineages = f" lineages: {founders};" if founders else ""
                head = f"decision div_{c.id} {{ action: divide; when: cell = {c.id}; daughters: {a}, {b};"
                lines.append(f"{head}{lineages} {ev}; confidence: 0.95 }}")
            elif c.dies is not None:
                after = c.dies - c.born
                head = f"decision die_{c.id} {{ action: die; when: cell = {c.id}; after: {after:.0f} min"
                lines.append(f"{head}; {ev}; confidence: 0.9 }}")
            else:
                name = f" name: {c.name};" if c.name != c.id else ""
                lines.append(
                    f"decision fate_{c.id} {{ action: differentiate; when: cell = {c.id}; to: {c.cell_type};"
                    f"{name} {ev}; confidence: 0.9 }}"
                )
        return "\n".join(lines) + "\n"
