# SPDX-License-Identifier: AGPL-3.0-or-later
"""Experiments: a perturbed organism run against the wild type.

An `experiment` block names what is removed (maternal factors, signals), added or
changed; the runner grows the wild type and the mutant with the same program, seed
and horizon, and reports every cell that decided differently: which decisions fired
in one and not the other, which cell type it took, whether it died. The earliest
differing cells are the founders whose mechanism the knockout touched; the rest are
their descendants. `expect:` carries the published phenotype so the reader can judge
whether the mechanism written in the program reproduces it.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from genomeos.ir import Experiment, Module, to_minutes
from genomeos.runtime.body import Body, Cell


@dataclass(slots=True)
class CellChange:
    cell: str
    wild_type: str
    mutant: str
    gained: list[str] = field(default_factory=list)  # decisions fired only in the mutant
    lost: list[str] = field(default_factory=list)  # decisions fired only in the wild type
    death: str = ""  # "gained" | "lost" | ""


@dataclass(slots=True)
class ExperimentResult:
    experiment: Experiment
    until: float
    wild_type: Body
    mutant: Body
    changes: list[CellChange] = field(default_factory=list)
    only_wild_type: list[str] = field(default_factory=list)
    only_mutant: list[str] = field(default_factory=list)

    @property
    def founders(self) -> list[CellChange]:
        """The earliest changed cells: none of their ancestors changed."""
        changed = {c.cell for c in self.changes}
        out = []
        for ch in self.changes:
            p = self.mutant.cells.get(ch.cell) or self.wild_type.cells.get(ch.cell)
            anc, root = p.parent if p else "", False
            while anc:
                if anc in changed:
                    root = True
                    break
                anc = (self.mutant.cells.get(anc) or self.wild_type.cells.get(anc)).parent
            if not root:
                out.append(ch)
        return out

    def early(self, before: float = 120.0) -> list[CellChange]:
        """Changed cells born in the founder window (the first ~100 minutes for the worm)."""
        out = []
        for ch in self.changes:
            c = self.mutant.cells.get(ch.cell) or self.wild_type.cells.get(ch.cell)
            if c is not None and c.born <= before:
                out.append(ch)
        return out

    def descendants_changed(self, cell: str) -> int:
        stack, n = [cell], 0
        while stack:
            c = (self.mutant.cells.get(stack.pop()) or Cell("", "", 0, 0, "")).children
            n += len(c)
            stack.extend(c)
        return n

    def asserts(self) -> list[dict]:
        from genomeos.runtime.body import evaluate_assert

        return [evaluate_assert(self.mutant, a) for a in self.experiment.asserts]

    def to_dict(self) -> dict:
        return {
            "experiment": self.experiment.name,
            "knockouts": self.experiment.knockouts,
            "adds": self.experiment.adds,
            "expect": self.experiment.expect,
            "until_min": self.until,
            "changed_cells": len(self.changes),
            "early": [
                {
                    "cell": f.cell,
                    "wild_type": f.wild_type,
                    "mutant": f.mutant,
                    "gained": f.gained,
                    "lost": f.lost,
                }
                for f in self.early()
            ],
            "founders": [
                {
                    "cell": f.cell,
                    "wild_type": f.wild_type,
                    "mutant": f.mutant,
                    "gained": f.gained,
                    "lost": f.lost,
                    "death": f.death,
                    "descendants": self.descendants_changed(f.cell),
                }
                for f in self.founders
            ],
            "only_wild_type": len(self.only_wild_type),
            "only_mutant": len(self.only_mutant),
            "mutant_summary": self.mutant.summary(),
            "wild_type_summary": self.wild_type.summary(),
            "asserts": self.asserts(),
            "evidence": {
                "kind": self.experiment.evidence.kind.value,
                "source": self.experiment.evidence.source,
                "confidence": self.experiment.confidence,
            },
        }

    def format(self) -> str:
        e = self.experiment
        what = ", ".join(["-" + k for k in e.knockouts] + ["+" + a for a in e.adds]) or "no change"
        head = f"experiment {e.name} ({what}) to {self.until:g} min"
        out = [f"{head}: {len(self.changes)} cells decided differently"]
        if e.expect:
            out.append(f"  published: {e.expect}")
        rows = self.early() or self.founders
        for f in rows[:12]:
            arrow = (
                f"{f.wild_type or 'Blastomere'} -> {f.mutant or 'Blastomere'}"
                if f.wild_type != f.mutant
                else ""
            )
            parts = [
                x
                for x in (
                    arrow,
                    f"lost {', '.join(f.lost)}" if f.lost else "",
                    f"gained {', '.join(f.gained)}" if f.gained else "",
                    f"death {f.death}" if f.death else "",
                )
                if x
            ]
            out.append(f"  {f.cell:<12} {'; '.join(parts)}  ({self.descendants_changed(f.cell)} descendants)")
        if len(rows) > 12:
            out.append(f"  ... {len(rows) - 12} more")
        if self.only_wild_type or self.only_mutant:
            a, b = len(self.only_wild_type), len(self.only_mutant)
            out.append(f"  cells only in the wild type: {a}; only in the mutant: {b}")
        wt, mu = self.wild_type.summary(), self.mutant.summary()
        out.append(f"  alive {wt['alive']:g} -> {mu['alive']:g}; deaths {wt['deaths']} -> {mu['deaths']}")
        for chk in self.asserts():
            mark = "ok  " if chk["ok"] else "FAIL"
            out.append(f"  assert {mark} {chk['assert']}  (got {chk.get('value', chk.get('error'))})")
        if e.evidence.source:
            out.append(f"  evidence: [{e.evidence.kind.value}] {e.evidence.source}")
        return "\n".join(out)


def _horizon(module: Module, ex: Experiment) -> float:
    if ex.until is not None:
        return ex.until
    if module.stages:
        return max(to_minutes(s.start, s.unit) for s in module.stages)
    return 800.0


def run_experiment(
    module: Module, ex: Experiment, seed: int | None = None, means: bool = True
) -> ExperimentResult:
    until = _horizon(module, ex)
    wt = Body(module, seed=seed, means=means).run(until=until)
    mu = Body(
        module,
        seed=seed,
        means=means,
        knockouts=set(ex.knockouts),
        adds=set(ex.adds),
        environment=ex.environment,
    ).run(until=until)
    res = ExperimentResult(ex, until, wt, mu)
    for name, w in wt.cells.items():
        m = mu.cells.get(name)
        if m is None:
            res.only_wild_type.append(name)
            continue
        gained = sorted(set(m.fired) - set(w.fired))
        lost = sorted(set(w.fired) - set(m.fired))
        death = (
            "gained"
            if m.dies_at is not None and w.dies_at is None
            else "lost"
            if w.dies_at is not None and m.dies_at is None
            else ""
        )
        if gained or lost or w.cell_type != m.cell_type or death:
            res.changes.append(CellChange(name, w.cell_type, m.cell_type, gained, lost, death))
    res.only_mutant = sorted(set(mu.cells) - set(wt.cells))
    return res


def run_experiments(
    module: Module, names: list[str] | None = None, seed: int | None = None
) -> list[ExperimentResult]:
    return [
        run_experiment(module, ex, seed=seed) for ex in module.experiments if not names or ex.name in names
    ]
