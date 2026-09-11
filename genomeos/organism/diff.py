"""LineageDiff: how far a grown body is from the reference lineage, as numbers.

born      reference cells (within the window) that the body also produced, by name
topology  matched cells whose parent is the same
timing    division-time error on matched cells (source time axis)
fates     terminal cells whose cell type matches the reference tissue
deaths    programmed deaths reproduced
counts    alive-cell curve at checkpoints
"""

from __future__ import annotations

import statistics
from dataclasses import dataclass, field

from genomeos.runtime.body import Body

from .reference import ReferenceLineage

CHECKPOINTS = (0.5, 50, 100, 200, 350, 500, 800, 2000, 4000, 5700)


@dataclass(slots=True)
class LineageDiff:
    until: float
    reference_cells: int = 0
    matched: int = 0
    missing: list[str] = field(default_factory=list)
    extra: list[str] = field(default_factory=list)
    parent_mismatch: list[str] = field(default_factory=list)
    timing_errors: list[float] = field(default_factory=list)
    fates_checked: int = 0
    fates_correct: int = 0
    fates_unknown: int = 0
    fate_confusion: dict[str, int] = field(default_factory=dict)
    deaths_expected: int = 0
    deaths_matched: int = 0
    counts: list[dict] = field(default_factory=list)

    @property
    def born_fraction(self) -> float:
        return self.matched / self.reference_cells if self.reference_cells else 0.0

    @property
    def fate_accuracy(self) -> float | None:
        return self.fates_correct / self.fates_checked if self.fates_checked else None

    @property
    def timing_median(self) -> float | None:
        return statistics.median(self.timing_errors) if self.timing_errors else None

    @property
    def timing_p90(self) -> float | None:
        if not self.timing_errors:
            return None
        xs = sorted(self.timing_errors)
        return xs[min(len(xs) - 1, int(0.9 * len(xs)))]

    def to_dict(self) -> dict:
        return {
            "until_min": self.until,
            "reference_cells": self.reference_cells,
            "matched": self.matched,
            "born_fraction": round(self.born_fraction, 4),
            "missing": len(self.missing),
            "missing_examples": self.missing[:10],
            "extra": len(self.extra),
            "extra_examples": self.extra[:10],
            "parent_mismatch": len(self.parent_mismatch),
            "timing_median_min": self.timing_median,
            "timing_p90_min": self.timing_p90,
            "fates_checked": self.fates_checked,
            "fates_correct": self.fates_correct,
            "fates_unknown": self.fates_unknown,
            "fate_accuracy": None if self.fate_accuracy is None else round(self.fate_accuracy, 4),
            "fate_confusion": dict(sorted(self.fate_confusion.items(), key=lambda kv: -kv[1])[:10]),
            "deaths_expected": self.deaths_expected,
            "deaths_matched": self.deaths_matched,
            "counts": self.counts,
        }

    def format(self) -> str:
        acc = "n/a" if self.fate_accuracy is None else f"{self.fate_accuracy:.1%}"
        tm = (
            "n/a"
            if self.timing_median is None
            else f"median {self.timing_median:.0f} min, p90 {self.timing_p90:.0f} min"
        )
        out = [
            f"against the reference lineage to {self.until:.0f} min:",
            f"  cells born      {self.matched}/{self.reference_cells} ({self.born_fraction:.1%}); "
            f"missing {len(self.missing)}, extra {len(self.extra)}, "
            f"parent mismatches {len(self.parent_mismatch)}",
            f"  division timing {tm}",
            f"  terminal fates  {self.fates_correct}/{self.fates_checked} correct ({acc}); "
            f"{self.fates_unknown} UNKNOWN",
            f"  deaths          {self.deaths_matched}/{self.deaths_expected}",
            "  alive cells     "
            + "  ".join(f"t={c['t']:g}: {c['body']:g}/{c['reference']}" for c in self.counts),
        ]
        if self.fate_confusion:
            worst = ", ".join(f"{k} ×{v}" for k, v in list(self.to_dict()["fate_confusion"].items())[:5])
            out.append(f"  confusions      {worst}")
        return "\n".join(out)


def compare(body: Body, ref: ReferenceLineage, until: float | None = None) -> LineageDiff:
    until = body.time if until is None else until
    d = LineageDiff(until)
    ref_cells = {c.id: c for c in ref.cells.values() if c.born <= until}
    d.reference_cells = len(ref_cells)
    body_cells = {c.name: c for c in body.cells.values() if c.born <= until}
    for cid, rc in ref_cells.items():
        bc = body_cells.get(cid)
        if bc is None:
            d.missing.append(cid)
            continue
        d.matched += 1
        if (rc.parent or "") != (bc.parent or ""):
            d.parent_mismatch.append(cid)
        if rc.divides is not None and bc.divides_at is not None and rc.divides <= until:
            d.timing_errors.append(abs(rc.divides - bc.divides_at))
        if rc.terminal and rc.dies is None:
            d.fates_checked += 1
            want = rc.cell_type
            got = bc.cell_type if (bc.terminal or bc.quiescent) else "UNKNOWN"
            if got == want:
                d.fates_correct += 1
            elif got in ("", "UNKNOWN", "Blastomere", "Progenitor", "Zygote"):
                d.fates_unknown += 1
            else:
                key = f"{want}->{got}"
                d.fate_confusion[key] = d.fate_confusion.get(key, 0) + 1
        if rc.dies is not None and rc.dies <= until:
            d.deaths_expected += 1
            if bc.dies_at is not None:
                d.deaths_matched += 1
    d.extra = sorted(set(body_cells) - set(ref_cells))
    for t in CHECKPOINTS:
        if t <= until:
            d.counts.append({"t": t, "body": body.count_at(t), "reference": ref.count_at(t)})
    return d
