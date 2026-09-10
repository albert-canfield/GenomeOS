"""Biological debugger (task 4.4): step, breakpoints, and an evidence trace.

Works on the network runtime (species levels) and the cell-ageing runtime
(per-cell clocks). A breakpoint is a comparison `name op value` on a state
variable; when it fires, `explain()` follows the change back to the rules or
parameters that produced it, with their evidence.
"""

from __future__ import annotations

import operator
import re
from dataclasses import dataclass, field

from genomeos.ir import Action, Module, Parameter, Rule

from .cell import CELL_TYPES, GLOBAL_PARAMS, CellRuntime, CellState
from .grn import NetworkRuntime, _hill

_OPS = {
    ">": operator.gt,
    ">=": operator.ge,
    "<": operator.lt,
    "<=": operator.le,
    "==": operator.eq,
    "!=": operator.ne,
}
_BP = re.compile(r"^\s*([\w.\-]+)\s*(>=|<=|==|!=|>|<)\s*([-+0-9.eE]+|true|false)\s*$")


@dataclass(frozen=True, slots=True)
class Breakpoint:
    variable: str
    op: str
    value: float

    @classmethod
    def parse(cls, text: str) -> Breakpoint:
        m = _BP.match(text)
        if not m:
            raise ValueError(f"bad breakpoint {text!r}; expected 'name op number' e.g. 'TetR > 50'")
        raw = m.group(3).lower()
        val = 1.0 if raw == "true" else 0.0 if raw == "false" else float(raw)
        return cls(m.group(1), m.group(2), val)

    def hit(self, state: dict[str, float]) -> bool:
        if self.variable not in state:
            return False
        return _OPS[self.op](float(state[self.variable]), self.value)

    def __str__(self) -> str:
        return f"{self.variable} {self.op} {self.value:g}"


@dataclass(slots=True)
class TraceLine:
    time: float
    subject: str
    message: str
    evidence: str = ""
    confidence: float | None = None

    def format(self) -> str:
        ev = f"  [{self.evidence}]" if self.evidence else ""
        conf = f" conf={self.confidence:.2f}" if self.confidence is not None else ""
        return f"t={self.time:8.3f}  {self.subject:<18} {self.message}{ev}{conf}"


@dataclass(slots=True)
class Session:
    kind: str
    time: float = 0.0
    state: dict[str, float] = field(default_factory=dict)
    breakpoints: list[Breakpoint] = field(default_factory=list)
    trace: list[TraceLine] = field(default_factory=list)
    stopped_at: Breakpoint | None = None


class NetworkDebugger:
    """Step a BioLang regulatory module and explain what drives each species."""

    def __init__(self, module: Module, context: dict[str, str] | None = None, dt: float = 0.01) -> None:
        self.module = module
        self.vm = NetworkRuntime(module, context=context)
        self.dt = dt
        self.session = Session("network", state=dict.fromkeys(self.vm.species, 0.0))

    def set_initial(self, initial: dict[str, float]) -> None:
        self.session.state.update(initial)

    def add_breakpoint(self, text: str) -> Breakpoint:
        bp = Breakpoint.parse(text)
        if bp.variable not in self.session.state:
            raise ValueError(f"unknown species {bp.variable!r}; known: {self.vm.species}")
        self.session.breakpoints.append(bp)
        return bp

    def step(self, hours: float | None = None) -> Breakpoint | None:
        """Advance by `hours` (default one dt); stop early at the first breakpoint hit."""
        n = max(1, int(round((hours or self.dt) / self.dt)))
        s = self.session
        s.stopped_at = None
        for _ in range(n):
            traj = self.vm.run(hours=self.dt, dt=self.dt, initial=s.state, record_every=1)
            s.state = traj.final()
            s.time += self.dt
            for bp in s.breakpoints:
                if bp.hit(s.state):
                    s.stopped_at = bp
                    s.trace.append(
                        TraceLine(
                            s.time, bp.variable, f"breakpoint {bp} hit (value {s.state[bp.variable]:.3f})"
                        )
                    )
                    return bp
        return None

    def run_until_break(self, max_hours: float = 100.0) -> Breakpoint | None:
        return self.step(max_hours)

    def contributions(self, species: str) -> list[tuple[Rule | Parameter | str, float, str]]:
        """Terms driving d[species]/dt right now: (rule-or-parameter, value, description)."""
        st = self.session.state
        out: list[tuple[Rule | Parameter | str, float, str]] = []
        if species.endswith(".mRNA"):
            gid = species[: -len(".mRNA")]
            gene = self.module.entities[gid]
            max_rate = float(gene.attrs.get("max_rate", 0.0))
            out.append(("basal", gene.basal_rate, f"basal transcription of {gid}"))
            for r in self.vm.activators.get(gid, []):
                h = _hill(st.get(r.source, 0.0), r.threshold, r.hill)
                out.append((r, r.strength * h * max_rate, f"{r.source} activates {gid}: Hill={h:.3f}"))
            for r in self.vm.inhibitors.get(gid, []):
                h = _hill(st.get(r.source, 0.0), r.threshold, r.hill)
                out.append(
                    (
                        r,
                        -r.strength * h * max_rate,
                        f"{r.source} inhibits {gid}: Hill={h:.3f} ({r.source}={st.get(r.source, 0):.2f})",
                    )
                )
            p = self.module.parameters.get("mrna_half_life")
            out.append(
                (
                    p or "mrna_half_life",
                    -st[species] * 0.693 / (p.value if p else self.vm.params["mrna_half_life"]),
                    "mRNA decay",
                )
            )
        else:
            for gid in self.vm.produces.get(species, []):
                r = next(
                    (
                        x
                        for x in self.vm.active_rules
                        if x.action is Action.PRODUCE and x.source == gid and x.target == species
                    ),
                    None,
                )
                out.append(
                    (
                        r or f"{gid} produces {species}",
                        self.vm.params["translation_rate"] * st[f"{gid}.mRNA"],
                        f"translation from {gid}.mRNA",
                    )
                )
            prot = self.module.entities.get(species)
            hl = (
                prot.half_life_h
                if prot is not None
                and not isinstance(prot.half_life_h, type(None))
                and prot.half_life_h.__class__.__name__ != "_Unknown"
                else None
            )
            out.append(
                (
                    prot or species,
                    -st[species] * 0.693 / (hl or self.vm.params["protein_half_life"]),
                    "protein decay",
                )
            )
        out.sort(key=lambda x: -abs(x[1]))
        return out

    def explain(self, species: str) -> list[TraceLine]:
        lines = []
        for term, value, desc in self.contributions(species):
            if isinstance(term, Rule | Parameter) or hasattr(term, "evidence"):
                lines.append(
                    TraceLine(
                        self.session.time,
                        species,
                        f"{value:+.3f}  {desc}",
                        f"{term.evidence.kind.value}: {term.evidence.source or 'no source'}",
                        term.confidence,
                    )
                )
            else:
                lines.append(
                    TraceLine(self.session.time, species, f"{value:+.3f}  {desc}", "runtime default", None)
                )
        self.session.trace.extend(lines)
        return lines


class AgeingDebugger:
    """Step one cell's clocks and explain which parameters moved them."""

    def __init__(self, cell_type: str = "fibroblast", seed: int = 1, dt_years: float = 0.5) -> None:
        self.rt = CellRuntime(seed=seed)
        self.cell: CellState = self.rt.new_cell(cell_type)
        self.dt = dt_years
        self.session = Session("ageing", state=self._state())

    def _state(self) -> dict[str, float]:
        c = self.cell
        return {
            "age": c.age_years,
            "telomere_bp": c.telomere_bp,
            "divisions": float(c.divisions),
            "somatic_mutations": float(c.somatic_mutations),
            "epigenetic_age": c.epigenetic_age,
            "senescent": 1.0 if c.senescent else 0.0,
        }

    def add_breakpoint(self, text: str) -> Breakpoint:
        bp = Breakpoint.parse(text)
        if bp.variable not in self.session.state:
            raise ValueError(f"unknown variable {bp.variable!r}; known: {list(self.session.state)}")
        self.session.breakpoints.append(bp)
        return bp

    def step(self, years: float | None = None) -> Breakpoint | None:
        n = max(1, int(round((years or self.dt) / self.dt)))
        s = self.session
        s.stopped_at = None
        for _ in range(n):
            before = self._state()
            self.rt.tick(self.cell, self.dt)
            s.state = self._state()
            s.time = self.cell.age_years
            for key in ("divisions", "senescent"):
                if s.state[key] != before[key]:
                    s.trace.append(self._event_line(key, before, s.state))
            for bp in s.breakpoints:
                if bp.hit(s.state):
                    s.stopped_at = bp
                    s.trace.append(
                        TraceLine(
                            s.time, bp.variable, f"breakpoint {bp} hit (value {s.state[bp.variable]:.3f})"
                        )
                    )
                    return bp
        return None

    def _event_line(self, key: str, before: dict[str, float], after: dict[str, float]) -> TraceLine:
        p = CELL_TYPES[self.cell.cell_type]
        if key == "divisions":
            par = p.telomere_loss_per_division_bp
            return TraceLine(
                self.cell.age_years,
                "divide",
                f"{int(after['divisions'] - before['divisions'])} division(s); "
                f"telomere {before['telomere_bp']:.0f} -> {after['telomere_bp']:.0f} bp",
                f"{par.evidence.kind.value}: {par.evidence.source}",
                par.confidence,
            )
        par = GLOBAL_PARAMS["telomere_senescence_bp"]
        return TraceLine(
            self.cell.age_years,
            "senesce",
            f"cell became senescent at telomere {after['telomere_bp']:.0f} bp, "
            f"{int(after['somatic_mutations'])} mutations",
            f"{par.evidence.kind.value}: {par.evidence.source}",
            par.confidence,
        )

    def explain(self) -> list[TraceLine]:
        """The parameters behind each clock right now, with evidence."""
        p = CELL_TYPES[self.cell.cell_type]
        params = [
            p.divisions_per_year,
            p.telomere_loss_per_division_bp,
            p.telomerase_compensation,
            p.mutations_per_year,
            GLOBAL_PARAMS["telomere_senescence_bp"],
            GLOBAL_PARAMS["senescence_hazard_base"],
            GLOBAL_PARAMS["epigenetic_pace_sd"],
        ]
        lines = [
            TraceLine(
                self.cell.age_years,
                par.name,
                f"= {par.value:g} {par.unit}",
                f"{par.evidence.kind.value}: {par.evidence.source}",
                par.confidence,
            )
            for par in params
        ]
        self.session.trace.extend(lines)
        return lines
