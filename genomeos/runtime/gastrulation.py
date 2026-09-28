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

# Census comparison, pre-registered 2026-09-28 (docs/DESIGN-MINIMAL-CELL.md, "Gastrulation against a
# measured census"). Constants only: the model below is not changed by them. The primary census is the
# one measured human gastrula, Carnegie stage 7 (Tyser et al. 2021, Nature 600:285, ArrayExpress
# E-MTAB-9388, 1,195 Smart-seq2 cells, author cell-type labels); the module's genes and proteins are
# human (GENCODE, UniProt P48431, O15178, Q9H6I2, Q96S42). The mouse atlas (Pijuan-Sala et al. 2019,
# Nature 566:490, E-MTAB-6967) is reported beside it at E7.0, E7.25 and E7.5 without a verdict.
# Each census proportion is a share of sampled cells renormalised over the cells a mapping assigns to
# a layer; the interval is the lowest and highest share over every registered variant (two mappings,
# and for human two site filters), and the model passes a layer when its share lies within
# CENSUS_TOLERANCE of that interval. Any layer outside it on the human census falsifies the model's
# default proportions against a measured gastrula.
CENSUS_TOLERANCE = 0.05
CENSUS_MODEL_RUN = {"cells": 120, "hours": 40.0, "dt": 0.05, "nodal_max": 6.0, "decay_length": 0.35}
# Tyser 2021 author labels ("Characteristics[inferred cell type - authors labels]" in the SDRF).
CENSUS_HUMAN_EXTRAEMBRYONIC = ("yolk sac mesoderm", "hemogenic endothelial progenitor", "erythrocyte")
CENSUS_HUMAN_MAPPINGS = {
    # differentiated: only cells that have left the epiblast and the streak
    "differentiated": {
        "ectoderm": ("ectodermal cell",),
        "mesoderm": ("nascent mesoderm", "emergent mesoderm", "advanced mesoderm", "axial mesoderm"),
        "endoderm": ("endodermal cell",),
    },
    # marker: the model's own reading, the SOX2 default holds the epiblast and TBXT holds the streak
    "marker": {
        "ectoderm": ("epiblast cell", "ectodermal cell"),
        "mesoderm": (
            "primitive streak",
            "nascent mesoderm",
            "emergent mesoderm",
            "advanced mesoderm",
            "axial mesoderm",
        ),
        "endoderm": ("endodermal cell",),
    },
}
CENSUS_HUMAN_SITES = {"all": ("rostral", "caudal", "yolk sac"), "disc": ("rostral", "caudal")}
# Pijuan-Sala 2019 `celltype` labels (meta.tab), secondary, no verdict.
CENSUS_MOUSE_STAGES = ("E7.0", "E7.25", "E7.5")
_MOUSE_ECTO = (
    "Rostral neurectoderm",
    "Caudal neurectoderm",
    "Surface ectoderm",
    "Forebrain/Midbrain/Hindbrain",
    "Spinal cord",
    "Neural crest",
)
_MOUSE_MESO = (
    "Nascent mesoderm", "Mixed mesoderm", "Intermediate mesoderm", "Paraxial mesoderm", "Somitic mesoderm",
    "Pharyngeal mesoderm", "Caudal Mesoderm", "Notochord", "Mesenchyme", "Cardiomyocytes",
)  # fmt: skip
_MOUSE_ENDO = ("Def. endoderm", "Gut")
CENSUS_MOUSE_MAPPINGS = {
    "differentiated": {"ectoderm": _MOUSE_ECTO, "mesoderm": _MOUSE_MESO, "endoderm": _MOUSE_ENDO},
    "marker": {
        "ectoderm": ("Epiblast", *_MOUSE_ECTO),
        "mesoderm": ("Primitive Streak", *_MOUSE_MESO),
        "endoderm": ("Anterior Primitive Streak", *_MOUSE_ENDO),
    },
}


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
    exponentially from the x = 0 pole (the node / primitive streak).

    The module's expected_* proportions are unsourced guesses, and the defaults
    here miss the mesoderm one: 0.10 against 0.35 (endoderm 0.43 against 0.20).
    Each fate holds a fixed band of NODAL level, so decay_length and nodal_max
    (neither sourced, both in model units) only move the band edges along the
    axis, and two free numbers can fit any three proportions. Meeting the
    expectation by setting them would therefore show nothing."""
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
