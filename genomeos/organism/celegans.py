"""Minimal organism (Phase 5): the invariant early lineage of C. elegans.

Since v0.3 the lineage is a BioLang program (data/demo/celegans_lineage.bio:
one zygote, cited cycle timers, founder divisions by name) run by the Body
runtime; this module keeps the earlier `Lineage` view over the result.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from genomeos.ir import Module
from genomeos.lang import parse_file
from genomeos.runtime.body import Body

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
    body: Body | None = None

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


def run_lineage(
    module_path: str | Path = "data/demo/celegans_lineage.bio", until_min: float = 150.0, max_cells: int = 600
) -> Lineage:
    module = parse_file(module_path)
    body = Body(module, max_cells=max_cells).run(until=until_min)
    lin = Lineage(module=module, body=body)
    for c in body.cells.values():
        lineage = "P" if c.lineage in ("P", "P4", "EMS", "") else c.lineage
        lc = LineageCell(c.name, lineage, c.generation, c.born, parent=c.parent)
        if c.divides_at is not None and c.divides_at <= until_min:
            lc.divides_at = c.divides_at
        lc.children = list(c.children)
        lin.cells[c.name] = lc
    return lin
