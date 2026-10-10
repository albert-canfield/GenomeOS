# SPDX-License-Identifier: Apache-2.0
# Part of the BioLang engine (language, IR, VM, standard library); see LICENSING.md.
"""libRoadRunner adapter for SBML: the standard simulator behind the stand-in's interface.

`runtime/sbml.py` integrates SBML with its own MathML evaluator and a fixed-step RK4. That is
enough for small, non-stiff models and fails on stiff ones. libRoadRunner (Apache-2.0, an optional
extra: `uv sync --extra roadrunner`) compiles the same SBML and integrates it with CVODE.

`RoadRunnerRuntime` has the interface of `SbmlRuntime`: built from an `SbmlModel`, it exposes
`model` and `dynamic` and returns a `Trajectory` from `run(duration, dt, record_every, initial)` on the
same time grid (0, then every `record_every * dt`). What a caller changed on the `SbmlModel` after
parsing -- an initial value, a species held constant (a knockout) -- is carried over to libRoadRunner;
everything else comes from the SBML file itself, so initial assignments, events and the parts of
MathML the stand-in does not read are libRoadRunner's to interpret.

It is chosen by name (`sbml.runtime(model, engine="roadrunner")`) and never replaces the stand-in
silently: without the extra, building one raises `EngineUnavailableError` naming the extra.
"""

from __future__ import annotations

from pathlib import Path

from .grn import Trajectory
from .sbml import SbmlModel

ENGINE = "roadrunner"
EXTRA = "roadrunner"
HOW_TO_ENABLE = f"install libRoadRunner with `uv sync --extra {EXTRA}`"


class EngineUnavailableError(ImportError):
    """The named engine's optional extra is not installed."""


def _import():
    try:
        import roadrunner  # type: ignore[import-not-found]
    except ImportError as e:  # the extra is absent: say which one, never fall back
        raise EngineUnavailableError(f"libRoadRunner is not installed; {HOW_TO_ENABLE}") from e
    return roadrunner


def available() -> bool:
    """True when libRoadRunner can be imported in this environment."""
    try:
        _import()
    except EngineUnavailableError:
        return False
    return True


def version() -> str:
    return str(_import().__version__)


class RoadRunnerRuntime:
    """An `SbmlModel` run on libRoadRunner (CVODE), with `SbmlRuntime`'s interface."""

    engine = ENGINE

    def __init__(self, model: SbmlModel, integrator: str = "cvode") -> None:
        self._rr_module = _import()
        if not model.source:
            raise ValueError("the SbmlModel has no source file; build it with SbmlModel.from_file")
        self.model = model
        self.integrator = integrator
        self.sbml = Path(model.source).read_text()
        pristine = SbmlModel.from_file(model.source)
        ruled = {v for v, _ in model.rate_rules}
        # the same set, in the same order, as SbmlRuntime.dynamic
        self.dynamic = [
            s
            for s in model.species
            if (s not in model.boundary and s not in model.constant_species) or s in ruled
        ] + [v for v in ruled if v not in model.species]
        # what the caller changed after parsing: initial values, and species now held constant
        self._changed = {
            s: v for s, v in model.species.items() if s not in pristine.species or pristine.species[s] != v
        }
        self._held = sorted(
            (model.constant_species | model.boundary) - (pristine.constant_species | pristine.boundary)
        )

    def _symbol(self, sid: str) -> str:
        """libRoadRunner's name for the value the stand-in reports: a species' concentration, or its
        amount when the SBML says hasOnlySubstanceUnits; any other variable by its id."""
        if sid in self.model.species and sid not in self.model.amount_units:
            return f"[{sid}]"
        return sid

    def _load(self, initial: dict[str, float] | None):
        """A fresh RoadRunner for one run, so that nothing set for one run leaks into the next."""
        rr = self._rr_module.RoadRunner(self.sbml)
        rr.setIntegrator(self.integrator)
        for sid in self._held:
            rr.setBoundary(sid, True, False)
        if self._held:
            rr.regenerateModel()
        for sid, value in {**self._changed, **(initial or {})}.items():
            rr.setValue(f"init({self._symbol(sid)})", float(value))
        rr.reset()
        return rr

    def settings(self) -> dict[str, object]:
        """The integrator and its settings as libRoadRunner reports them, for a result to record."""
        rr = self._load(None)  # held while the integrator is read: it belongs to this instance
        integ = rr.getIntegrator()
        names = ("relative_tolerance", "absolute_tolerance", "stiff", "maximum_num_steps")
        out = {
            "engine": ENGINE,
            "libroadrunner": version(),
            "integrator": integ.getName(),
            **{n: integ.getValue(n) for n in names if n in integ.getSettings()},
        }
        del integ, rr
        return out

    def run(
        self,
        duration: float,
        dt: float = 0.01,
        record_every: int = 10,
        initial: dict[str, float] | None = None,
    ) -> Trajectory:
        rr = self._load(initial)
        records = int(round(duration / dt)) // record_every
        symbols = [self._symbol(s) for s in self.dynamic]
        traj = Trajectory(times=[], species=list(self.dynamic), levels={s: [] for s in self.dynamic})
        if records == 0:  # nothing to integrate: the initial state is the whole trajectory
            rows = [[0.0, *(rr.getValue(sym) for sym in symbols)]]
        else:
            rows = rr.simulate(0.0, records * record_every * dt, records + 1, ["time", *symbols])
        for row in rows:
            traj.times.append(float(row[0]))
            for i, s in enumerate(self.dynamic, start=1):
                traj.levels[s].append(float(row[i]))
        return traj
