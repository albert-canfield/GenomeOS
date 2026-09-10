"""Minimal organism (Phase 5): the invariant early lineage of C. elegans.

C. elegans is the right first organism: 959 somatic cells, an invariant,
fully mapped lineage (Sulston 1983), and every cell named. This engine runs
the early embryo as discrete division events from a zygote bootstrap state,
with cycle lengths taken from a BioLang module carrying their evidence.

Nomenclature follows Sulston: AB → ABa/ABp → ABal/ABar/ABpl/ABpr …;
P0 → AB + P1; P1 → EMS + P2; EMS → MS + E; P2 → C + P3; P3 → D + P4.
Daughters of AB/MS/E/C/D alternate anterior/posterior then left/right
naming (a/p, l/r) as in the real lineage's first rounds.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from genomeos.ir import Module
from genomeos.lang import parse_file

FOUNDERS = ("AB", "MS", "E", "C", "D")


@dataclass(slots=True)
class LineageCell:
    name: str
    lineage: str  # founder lineage: AB, MS, E, C, D, or P
    generation: int  # divisions since the founder
    born: float  # minutes after first cleavage
    divides_at: float | None = None
    parent: str = ""
    children: list[str] = field(default_factory=list)

    @property
    def alive_until(self) -> float:
        return self.divides_at if self.divides_at is not None else float("inf")


@dataclass(slots=True)
class Lineage:
    cells: dict[str, LineageCell] = field(default_factory=dict)
    module: Module | None = None

    def alive_at(self, t: float) -> list[LineageCell]:
        return [c for c in self.cells.values() if c.born <= t < c.alive_until]

    def count_at(self, t: float) -> int:
        return len(self.alive_at(t))

    def lineage_count_at(self, lineage: str, t: float) -> int:
        return sum(1 for c in self.alive_at(t) if c.lineage == lineage)

    def tree(self, root: str = "P0", depth: int = 4) -> list[str]:
        out: list[str] = []

        def walk(name: str, d: int) -> None:
            c = self.cells[name]
            div = f"  divides {c.divides_at:.0f} min" if c.divides_at is not None else ""
            out.append(f"{'  ' * d}{name} (born {c.born:.0f} min){div}")
            if d < depth:
                for ch in c.children:
                    walk(ch, d + 1)

        walk(root, 0)
        return out


def _cycle(module: Module, lineage: str, generation: int) -> float:
    p = module.parameters
    if lineage == "AB":
        return p["ab_cycle_min"].value * (p["ab_lengthening"].value ** max(0, generation - 1))
    if lineage == "MS":
        return p["ms_cycle_min"].value * (1.1 ** max(0, generation - 1))
    if lineage == "E":
        return p["e_cycle_min"].value * (1.15 ** max(0, generation - 1))
    if lineage == "C":
        return p["c_cycle_min"].value * (1.1 ** max(0, generation - 1))
    if lineage == "D":
        return p["d_cycle_min"].value * (1.1 ** max(0, generation - 1))
    return p["p_cycle_min"].value


def _daughter_names(cell: LineageCell) -> tuple[str, str]:
    """Sulston-style names: alternate a/p and l/r suffixes for somatic founders."""
    if cell.name == "EMS":
        return ("MS", "E")
    if cell.lineage == "P":
        return {"P0": ("AB", "P1"), "P1": ("EMS", "P2"), "P2": ("C", "P3"), "P3": ("D", "P4")}.get(
            cell.name, (cell.name + "a", cell.name + "p")
        )
    suffix = ("a", "p") if cell.generation % 2 == 0 else ("l", "r")
    return (cell.name + suffix[0], cell.name + suffix[1])


def run_lineage(
    module_path: str | Path = "data/demo/celegans_lineage.bio", until_min: float = 150.0, max_cells: int = 600
) -> Lineage:
    module = parse_file(module_path)
    lin = Lineage(module=module)
    p0 = LineageCell("P0", "P", 0, -1.0, divides_at=0.0)
    lin.cells["P0"] = p0
    queue = [p0]
    while queue:
        queue.sort(key=lambda c: c.divides_at if c.divides_at is not None else float("inf"))
        cell = queue.pop(0)
        if cell.divides_at is None or cell.divides_at > until_min or len(lin.cells) >= max_cells:
            continue
        t = cell.divides_at
        a_name, b_name = _daughter_names(cell)
        for name in (a_name, b_name):
            if cell.name == "P0":
                lineage = "AB" if name == "AB" else "P"
            elif cell.name == "P1":
                lineage = "P" if name == "P2" else "P"  # EMS is still a P-derived precursor
            elif cell.name == "EMS":
                lineage = name  # MS or E
            elif cell.name == "P2":
                lineage = "C" if name == "C" else "P"
            elif cell.name == "P3":
                lineage = "D" if name == "D" else "P"
            else:
                lineage = cell.lineage
            generation = 0 if name in FOUNDERS or lineage == "P" else cell.generation + 1
            child = LineageCell(name, lineage, generation, t, parent=cell.name)
            if name == "P4":
                child.divides_at = None  # germline precursor: no further divisions in the embryo window
            elif name == "EMS":
                child.divides_at = t + module.parameters["ems_cycle_min"].value
            elif name == "P1":
                child.divides_at = t + module.parameters["p1_cycle_min"].value
            elif lineage == "P":
                child.divides_at = t + module.parameters["p_cycle_min"].value
            elif name == "AB":
                child.divides_at = t + module.parameters["ab_cycle_min"].value
            else:
                child.divides_at = t + _cycle(module, lineage, generation)
            lin.cells[name] = child
            cell.children.append(name)
            queue.append(child)
    return lin
