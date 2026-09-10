"""Segmentation clock and wavefront (task 4.2).

Clock-and-wavefront (Cooke & Zeeman 1976): every presomitic cell carries an
oscillator (HES7/Notch) while it sits in the FGF-rich posterior zone; as the
growth zone moves away and FGF falls below threshold, the cell freezes at its
current phase. Cells frozen in the same half-cycle form one segment, so the
segment count is elongation time / period.

The oscillator here is a phase oscillator (period from the BioLang module);
the FGF wavefront is a moving source on the 1-D spatial field. Everything
tunable comes from `data/demo/segmentation_clock.bio` with its evidence.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

from genomeos.ir import Module
from genomeos.lang import parse_file

from .spatial import SpatialRuntime
from .uncertainty import report_for_network


@dataclass(slots=True)
class SegmentationResult:
    hours: float
    period_h: float
    segments: list[tuple[int, int]]  # (start_x, end_x) of each segment
    frozen_phase: list[float]
    module: Module
    fields: dict[str, list[float]] = field(default_factory=dict)

    @property
    def count(self) -> int:
        return len(self.segments)

    @property
    def frozen_cells(self) -> int:
        return sum(1 for x in self.frozen_phase if x == x)  # not NaN

    @property
    def expected(self) -> float:
        """Segments expected in the frozen region: cells / (wavefront speed × period)."""
        speed = self.module.parameters["wavefront_speed"].value
        return self.frozen_cells / (speed * self.period_h)

    def uncertainty(self):
        return report_for_network(self.module, self.module.rules)


def run_segmentation(
    module_path: str | Path = "data/demo/segmentation_clock.bio",
    length: int = 120,
    hours: float = 60.0,
    dt: float = 0.1,
) -> SegmentationResult:
    module = parse_file(module_path)
    period = module.parameters["period_h"].value
    speed = module.parameters["wavefront_speed"].value  # cells per hour
    threshold = module.parameters["fgf_threshold"].value

    rt = SpatialRuntime(length, 1)
    fgf = rt.add_field("fgf", diffusion=0.5, decay=0.15)
    rt.fill_cells("presomitic")
    phase = [0.0] * length
    frozen = [False] * length
    frozen_phase = [math.nan] * length
    omega = 2 * math.pi / period

    # the wavefront (growth zone) starts at x = 0 and moves toward +x; cells are
    # laid down behind it. FGF is emitted at the current front position.
    front = 0.0
    t = 0.0
    steps = int(round(hours / dt))
    for _ in range(steps):
        fx = min(length - 1, int(front))
        fgf.add(fx, 0, 4.0 * dt)
        fgf.step(dt)
        for x in range(length):
            if frozen[x]:
                continue
            if x > front:  # not yet formed tissue: sits in the zone, oscillating in sync
                phase[x] = (phase[x] + omega * dt) % (2 * math.pi)
                continue
            if fgf.at(x, 0) >= threshold:
                phase[x] = (phase[x] + omega * dt) % (2 * math.pi)
            else:
                frozen[x] = True
                frozen_phase[x] = phase[x]
                rt.cells[x].cell_type = "anterior" if phase[x] < math.pi else "posterior"
        front += speed * dt
        t += dt

    # segments: a new segment starts at each anterior→... boundary where the
    # frozen phase wraps (posterior half followed by anterior half)
    segments: list[tuple[int, int]] = []
    start = None
    prev = None
    for x in range(length):
        if not frozen[x]:
            break
        half = "A" if frozen_phase[x] < math.pi else "P"
        if start is None:
            start = x
        elif prev == "P" and half == "A":
            segments.append((start, x))
            start = x
        prev = half
    if start is not None and prev is not None:
        segments.append((start, x))
    return SegmentationResult(hours, period, segments, frozen_phase, module, {"fgf": fgf.profile_x(0)})
