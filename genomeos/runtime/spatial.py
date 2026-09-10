"""Spatial engine (task 4.1): a 2-D field of cells with diffusing morphogens.

The minimal machinery development needs: cells with positions and state,
morphogen fields that diffuse and decay on a grid (explicit finite
differences), sources that emit into the field, and fate rules that read local
concentration. This is the in-house stand-in for CompuCell3D / Tissue Forge
(neither has Python 3.14 wheels); the API is kept small so an adapter can
replace it.

Units are arbitrary: grid spacing 1, time in hours.
"""

from __future__ import annotations

import math
import random
from collections.abc import Callable
from dataclasses import dataclass, field

from genomeos.ir import Evidence, EvidenceKind, Parameter

E_WOLPERT = Evidence(
    EvidenceKind.CURATED, "Wolpert 1969, J Theor Biol 25:1 (positional information / French flag)"
)


@dataclass(slots=True)
class Field2D:
    """A scalar concentration field on a W×H grid with diffusion and decay."""

    name: str
    width: int
    height: int
    diffusion: float = 0.2  # grid^2 / h
    decay: float = 0.01  # 1 / h
    values: list[list[float]] = field(default_factory=list)

    def __post_init__(self) -> None:
        if not self.values:
            self.values = [[0.0] * self.width for _ in range(self.height)]

    def at(self, x: int, y: int) -> float:
        return self.values[y][x]

    def add(self, x: int, y: int, amount: float) -> None:
        self.values[y][x] += amount

    def step(self, dt: float) -> None:
        """One explicit diffusion/decay step (stable for diffusion*dt <= 0.25)."""
        d = self.diffusion * dt
        if d > 0.25:
            raise ValueError("diffusion * dt must be <= 0.25 for a stable explicit step")
        v = self.values
        w, h = self.width, self.height
        new = [[0.0] * w for _ in range(h)]
        for y in range(h):
            row = v[y]
            up = v[y - 1] if y > 0 else row
            down = v[y + 1] if y < h - 1 else row
            for x in range(w):
                c = row[x]
                left = row[x - 1] if x > 0 else c
                right = row[x + 1] if x < w - 1 else c
                lap = left + right + up[x] + down[x] - 4 * c
                new[y][x] = max(0.0, c + d * lap - self.decay * dt * c)
        self.values = new

    def profile_x(self, y: int | None = None) -> list[float]:
        """Mean concentration along x (or a single row)."""
        if y is not None:
            return list(self.values[y])
        return [sum(self.values[yy][x] for yy in range(self.height)) / self.height for x in range(self.width)]


@dataclass(slots=True)
class SpatialCell:
    id: int
    x: int
    y: int
    cell_type: str = "undifferentiated"
    state: dict[str, float] = field(default_factory=dict)


FateRule = Callable[[SpatialCell, dict[str, float]], str | None]


@dataclass(slots=True)
class Source:
    field: str
    x: int
    y: int
    rate: float  # amount per hour


class SpatialRuntime:
    def __init__(self, width: int, height: int, seed: int | None = None) -> None:
        self.width, self.height = width, height
        self.fields: dict[str, Field2D] = {}
        self.cells: list[SpatialCell] = []
        self.sources: list[Source] = []
        self.fate_rules: list[tuple[str, FateRule]] = []
        self.rng = random.Random(seed)
        self.time = 0.0
        self.parameters: list[Parameter] = []

    def add_field(self, name: str, diffusion: float = 0.2, decay: float = 0.01) -> Field2D:
        f = Field2D(name, self.width, self.height, diffusion, decay)
        self.fields[name] = f
        return f

    def add_source(self, field_name: str, x: int, y: int, rate: float) -> None:
        self.sources.append(Source(field_name, x, y, rate))

    def fill_cells(self, cell_type: str = "undifferentiated") -> None:
        self.cells = [
            SpatialCell(y * self.width + x, x, y, cell_type)
            for y in range(self.height)
            for x in range(self.width)
        ]

    def add_fate_rule(self, name: str, rule: FateRule) -> None:
        self.fate_rules.append((name, rule))

    def local(self, cell: SpatialCell) -> dict[str, float]:
        return {name: f.at(cell.x, cell.y) for name, f in self.fields.items()}

    def step(self, dt: float) -> None:
        for s in self.sources:
            self.fields[s.field].add(s.x, s.y, s.rate * dt)
        for f in self.fields.values():
            f.step(dt)
        for cell in self.cells:
            env = self.local(cell)
            for _, rule in self.fate_rules:
                fate = rule(cell, env)
                if fate:
                    cell.cell_type = fate
                    break
        self.time += dt

    def run(self, hours: float, dt: float = 0.5) -> None:
        for _ in range(int(round(hours / dt))):
            self.step(dt)

    def census(self) -> dict[str, int]:
        out: dict[str, int] = {}
        for c in self.cells:
            out[c.cell_type] = out.get(c.cell_type, 0) + 1
        return out

    def type_map(self) -> list[str]:
        """One character per cell, row by row: first letter of the cell type."""
        rows = []
        for y in range(self.height):
            rows.append("".join(self.cells[y * self.width + x].cell_type[0] for x in range(self.width)))
        return rows

    def bands_along_x(self) -> list[tuple[str, int, int]]:
        """Contiguous bands of cell type along x in the middle row: (type, start, end)."""
        y = self.height // 2
        types = [self.cells[y * self.width + x].cell_type for x in range(self.width)]
        bands: list[tuple[str, int, int]] = []
        start = 0
        for x in range(1, self.width + 1):
            if x == self.width or types[x] != types[start]:
                bands.append((types[start], start, x))
                start = x
        return bands


def french_flag(
    width: int = 60,
    height: int = 10,
    hours: float = 200.0,
    seed: int = 0,
    high: float = 2.0,
    low: float = 0.5,
) -> SpatialRuntime:
    """Wolpert's French flag: a morphogen source at the left edge, cells read
    their local concentration against two thresholds and adopt blue/white/red."""
    rt = SpatialRuntime(width, height, seed)
    rt.add_field("morphogen", diffusion=0.4, decay=0.02)
    for y in range(height):
        rt.add_source("morphogen", 0, y, rate=1.0)
    rt.fill_cells()
    rt.parameters = [
        Parameter("threshold_high", high, "a.u.", E_WOLPERT, 0.6),
        Parameter("threshold_low", low, "a.u.", E_WOLPERT, 0.6),
    ]

    def fate(cell: SpatialCell, env: dict[str, float]) -> str | None:
        m = env["morphogen"]
        return "blue" if m >= high else "white" if m >= low else "red"

    rt.add_fate_rule("french_flag", fate)
    rt.run(hours, dt=0.5)
    return rt


def gradient_decay_length(diffusion: float, decay: float) -> float:
    """Characteristic length of a steady diffusion/decay gradient, sqrt(D/k)."""
    return math.sqrt(diffusion / decay)
