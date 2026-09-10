"""Germ-layer specification (task 4.3): a population of cells, each running the
gastrulation regulatory module, reading a NODAL gradient from one pole."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

from genomeos.ir import Module
from genomeos.lang import parse_file

from .grn import NetworkRuntime
from .uncertainty import report_for_network


@dataclass(slots=True)
class GastrulationResult:
    fates: list[str]
    positions: list[float]
    module: Module

    def proportions(self) -> dict[str, float]:
        n = len(self.fates)
        return {f: self.fates.count(f) / n for f in ("ectoderm", "mesoderm", "endoderm")}

    def uncertainty(self):
        return report_for_network(self.module, self.module.rules)


def run_gastrulation(
    module_path: str | Path = "data/demo/gastrulation.bio",
    cells: int = 120,
    hours: float = 40.0,
    dt: float = 0.05,
    nodal_max: float = 6.0,
    decay_length: float = 0.35,
) -> GastrulationResult:
    """Cells are spread along a normalised axis x in [0, 1]; NODAL falls off
    exponentially from the x = 0 pole (the node / primitive streak)."""
    module = parse_file(module_path)
    fates: list[str] = []
    positions: list[float] = []
    for i in range(cells):
        x = (i + 0.5) / cells
        nodal = nodal_max * math.exp(-x / decay_length)
        vm = NetworkRuntime(module)
        # NODAL is an external signal here (its gene has max 0): hold it at the local level
        traj = vm.run(hours=hours, dt=dt, initial={"Sox2": 1.0}, record_every=10**9, clamp={"Nodal": nodal})
        final = traj.final()
        levels = {"ectoderm": final["Sox2"], "mesoderm": final["Tbxt"], "endoderm": final["Sox17"]}
        fates.append(max(levels, key=levels.get))
        positions.append(x)
    return GastrulationResult(fates, positions, module)
