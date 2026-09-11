"""Composition: GenomeOS engines as process-bigraph processes.

process-bigraph (Vivarium 2.0) is an optional extra (`uv sync --extra compose`).
This module wraps the network runtime, the cell-ageing runtime and the Body
runtime (an organism grown from one cell) as Processes so that one Composite
can advance them on a shared clock. A network can gate an organism: when a
species crosses a threshold, a named factor appears in the organism's context
and its decisions read it.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

try:  # optional dependency
    from process_bigraph import Composite, Process, allocate_core

    AVAILABLE = True
except ImportError:  # pragma: no cover
    AVAILABLE = False
    Process = object  # type: ignore[assignment,misc]

from genomeos.ir import Module
from genomeos.lang import parse, parse_file

from .body import Body
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


class BodyProcess(Process):  # type: ignore[misc]
    """An organism program (BioLang v0.3, source text or a .bio path) advanced by the Body runtime per
    interval (hours).

    Outputs are cumulative counts reported as additive deltas; `fates` is the census by cell type.
    `gates` couple a network to the organism: {factor: {"species": name, "threshold": level}} sets the
    factor in the organism's context while the species is at or above the threshold."""

    config_schema = {
        "source": "string",
        "seed": {"_type": "integer", "_default": "-1"},
        "means": {"_type": "boolean", "_default": "false"},
        "knockouts": {"_type": "list[string]", "_default": "[]"},
        "gates": {"_type": "map[map[string]]", "_default": "{}"},
    }

    def initialize(self, config: dict[str, Any] | None = None) -> None:
        source = str(self.config["source"])
        if source.strip().endswith(".bio") and Path(source.strip()).is_file():
            self.module: Module = parse_file(source.strip())  # a file keeps its relative imports
        else:
            self.module = parse(source, name_hint="composite")
        seed = None if int(self.config["seed"]) < 0 else int(self.config["seed"])
        self.body = Body(
            self.module, seed=seed, means=bool(self.config["means"]), knockouts=set(self.config["knockouts"])
        )
        self.gates: dict[str, dict[str, str]] = dict(self.config["gates"])
        self.minutes = 0.0
        # a map update only reaches keys that already exist, and outputs are not visible in `state`, so the
        # process keeps its own census with every cell type of the program present from the start
        self.last_fates = {c.id: 0.0 for c in self.module.cell_types()} | {"UNKNOWN": 0.0}

    def inputs(self) -> dict:
        return {"levels": "map[float]", "minutes": "float"}

    def outputs(self) -> dict:
        return {
            "minutes": "float",
            "cells": "float",
            "deaths": "float",
            "born": "float",
            "turnover_per_day": "float",
            "fates": "map[float]",
        }

    def initial_state(self) -> dict:
        s = self.body.summary()  # at time 0: a root that divides at once already counts as not alive
        return {
            "minutes": 0.0,
            "cells": float(s["alive"]),
            "deaths": float(s["deaths"]),
            "born": float(s["cells_born"]),
            "turnover_per_day": float(s["turnover_per_day"]),
            "fates": dict(self.last_fates),
        }

    def _apply_gates(self, levels: dict[str, float]) -> None:
        for factor, gate in self.gates.items():
            level = levels.get(str(gate.get("species", "")))
            if level is None:
                continue
            on = float(level) >= float(gate.get("threshold", 0.0))
            self.body.set_factor(factor, "present" if on else None)

    def update(self, state: dict, interval: float) -> dict:
        self._apply_gates(dict(state.get("levels") or {}))
        before = self.body.summary()
        self.minutes += interval * 60.0
        self.body.run(until=self.minutes)
        after = self.body.summary()
        census = after["fates"]
        delta = {k: census.get(k, 0.0) - v for k, v in self.last_fates.items()}
        self.last_fates = {k: census.get(k, 0.0) for k in self.last_fates}
        return {
            "minutes": interval * 60.0,
            "cells": after["alive"] - before["alive"],
            "deaths": after["deaths"] - before["deaths"],
            "born": after["cells_born"] - before["cells_born"],
            "turnover_per_day": after["turnover_per_day"] - before["turnover_per_day"],
            "fates": delta,
        }


def build_body_composite(
    organism_source: str,
    network_source: str | None = None,
    gates: dict[str, dict[str, Any]] | None = None,
    seed: int | None = None,
    means: bool = False,
    knockouts: list[str] | None = None,
    dt: float = 0.01,
):
    """An organism grown by the Body runtime, optionally beside a regulatory network whose species gate
    its factors, on one composite clock (hours)."""
    if not AVAILABLE:
        raise ImportError("process-bigraph is not installed; run `uv sync --extra compose`")
    core = allocate_core()
    core.register_link("genomeos.NetworkProcess", NetworkProcess)
    core.register_link("genomeos.BodyProcess", BodyProcess)
    state: dict[str, Any] = {
        "body": {
            "_type": "process",
            "address": "local:genomeos.BodyProcess",
            "config": {
                "source": organism_source,
                "seed": -1 if seed is None else seed,
                "means": means,
                "knockouts": list(knockouts or []),
                "gates": {k: {kk: str(vv) for kk, vv in v.items()} for k, v in (gates or {}).items()},
            },
            "inputs": {"levels": ["levels"], "minutes": ["minutes"]},
            "outputs": {
                "minutes": ["minutes"],
                "cells": ["cells"],
                "deaths": ["deaths"],
                "born": ["born"],
                "turnover_per_day": ["turnover_per_day"],
                "fates": ["fates"],
            },
        }
    }
    if network_source:
        state["network"] = {
            "_type": "process",
            "address": "local:genomeos.NetworkProcess",
            "config": {"source": network_source, "dt": dt},
            "inputs": {"levels": ["levels"]},
            "outputs": {"levels": ["levels"]},
        }
    return Composite({"state": state}, core=core)


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
