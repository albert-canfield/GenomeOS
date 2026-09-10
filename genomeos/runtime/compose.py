"""Composition spike (task 2.5): GenomeOS engines as process-bigraph processes.

process-bigraph (Vivarium 2.0) is an optional extra (`uv sync --extra compose`).
This module wraps the network runtime and the cell-ageing runtime as
Processes so that one Composite can advance both on a shared clock.
"""

from __future__ import annotations

from typing import Any

try:  # optional dependency
    from process_bigraph import Composite, Process, allocate_core

    AVAILABLE = True
except ImportError:  # pragma: no cover
    AVAILABLE = False
    Process = object  # type: ignore[assignment,misc]

from genomeos.ir import Module
from genomeos.lang import parse

from .cell import CellRuntime, Environment
from .grn import NetworkRuntime


class NetworkProcess(Process):  # type: ignore[misc]
    """A BioLang regulatory module advanced by RK4 inside a Composite. State: species levels."""

    config_schema = {"source": "string", "dt": {"_type": "float", "_default": "0.01"}}

    def initialize(self, config: dict[str, Any] | None = None) -> None:
        self.module: Module = parse(self.config["source"], name_hint="composite")
        self.vm = NetworkRuntime(self.module)

    def inputs(self) -> dict:
        return {"levels": {s: "float" for s in self.vm.species}}

    def outputs(self) -> dict:
        return {"levels": {s: "float" for s in self.vm.species}}

    def initial_state(self) -> dict:
        return {"levels": dict.fromkeys(self.vm.species, 0.0)}

    def update(self, state: dict, interval: float) -> dict:
        levels = dict(state["levels"])
        traj = self.vm.run(
            hours=interval, dt=min(self.config["dt"], interval), initial=levels, record_every=10**9
        )
        final = traj.final()
        return {"levels": {s: final[s] - levels[s] for s in self.vm.species}}


class AgeingProcess(Process):  # type: ignore[misc]
    """A cell population's ageing clocks advanced per interval (hours -> years)."""

    config_schema = {
        "cell_type": {"_type": "string", "_default": "fibroblast"},
        "cells": {"_type": "integer", "_default": "200"},
        "seed": {"_type": "integer", "_default": "1"},
    }

    def initialize(self, config: dict[str, Any] | None = None) -> None:
        self.rt = CellRuntime(seed=int(self.config["seed"]), env=Environment())
        self.cells = [self.rt.new_cell(self.config["cell_type"]) for _ in range(int(self.config["cells"]))]

    def inputs(self) -> dict:
        return {"years": "float"}

    def outputs(self) -> dict:
        return {"years": "float", "mean_telomere_bp": "float", "senescent_fraction": "float"}

    def initial_state(self) -> dict:
        r = self.rt._report(self.cells, 0.0)
        return {
            "years": 0.0,
            "mean_telomere_bp": r.mean_telomere_bp,
            "senescent_fraction": r.senescent_fraction,
        }

    def update(self, state: dict, interval: float) -> dict:
        dt_years = interval / (24 * 365.25)
        for c in self.cells:
            self.rt.tick(c, dt_years)
        r = self.rt._report(self.cells, state["years"] + dt_years)
        # process-bigraph applies float updates additively
        return {
            "years": dt_years,
            "mean_telomere_bp": r.mean_telomere_bp - state.get("mean_telomere_bp", r.mean_telomere_bp),
            "senescent_fraction": r.senescent_fraction - state.get("senescent_fraction", 0.0),
        }


def build_composite(source: str, cell_type: str = "fibroblast", cells: int = 100, seed: int = 1):
    """One organism-like composite: a regulatory network and an ageing population on one clock."""
    if not AVAILABLE:
        raise ImportError("process-bigraph is not installed; run `uv sync --extra compose`")
    core = allocate_core()
    core.register_link("genomeos.NetworkProcess", NetworkProcess)
    core.register_link("genomeos.AgeingProcess", AgeingProcess)
    state = {
        "network": {
            "_type": "process",
            "address": "local:genomeos.NetworkProcess",
            "config": {"source": source, "dt": 0.01},
            "inputs": {"levels": ["levels"]},
            "outputs": {"levels": ["levels"]},
        },
        "ageing": {
            "_type": "process",
            "address": "local:genomeos.AgeingProcess",
            "config": {"cell_type": cell_type, "cells": cells, "seed": seed},
            "inputs": {"years": ["years"]},
            "outputs": {
                "years": ["years"],
                "mean_telomere_bp": ["mean_telomere_bp"],
                "senescent_fraction": ["senescent_fraction"],
            },
        },
    }
    return Composite({"state": state}, core=core)
