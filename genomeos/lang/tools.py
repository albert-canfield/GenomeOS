# SPDX-License-Identifier: Apache-2.0
"""The BioLang toolchain's own check, compile and run, written against the engine only
(language, IR, runtime), so `bio` stands without the GenomeOS application."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from genomeos.ir import Module
from genomeos.lang import parse_file


def load_module(path: str | Path) -> Module:
    p = Path(path)
    if p.suffix == ".json":
        return Module.from_dict(json.loads(p.read_text()))
    return parse_file(p)


def spark(values: list[float], width: int = 60) -> str:
    bars = "▁▂▃▄▅▆▇█"
    if not values:
        return ""
    step = max(1, len(values) // width)
    sample = values[::step][:width]
    lo, hi = min(sample), max(sample)
    if hi - lo < 1e-12:
        return bars[0] * len(sample)
    return "".join(bars[min(7, int((v - lo) / (hi - lo) * 7.999))] for v in sample)


def check_module(module: Module, context: dict[str, str] | None = None) -> str:
    """The compile report as text: counts, active rules in a context, confidence by kind, weak rules."""
    lines = [f"module {module.name}"]
    lines.append(
        f"  entities: {len(module.entities)}   rules: {len(module.rules)}   "
        f"parameters: {len(module.parameters)}   events: {len(module.events)}   "
        f"cell types: {len(module.cell_types())}"
    )
    if context:
        active = module.active_rules(context)
        silenced = sorted(module.silenced_genes(context))
        lines.append(
            f"  context {context}: {len(active)}/{len(module.rules)} rules active; "
            f"silenced: {silenced or 'none'}"
        )
        lines += [f"    {r.id}" for r in active]
    lines.append("  mean confidence by kind:")
    for kind, val in sorted(module.confidence_report().items()):
        lines.append(f"    {kind:<10} {'█' * int(val * 20):<20} {val:.2f}")
    unknowns = module.unknowns()
    if unknowns:
        lines.append(f"  UNKNOWN regions: {len(unknowns)}")
        lines += [f"    {u.id}  {u.locus}" for u in unknowns]
    weak = [r for r in module.rules if r.confidence < 0.5]
    if weak:
        lines.append(f"  rules with confidence < 0.5: {len(weak)}")
        lines += [f"    {r.id}  ({r.evidence.kind.value}: {r.evidence.source or 'no source'})" for r in weak]
    return "\n".join(lines)


def compile_module(module: Module) -> str:
    return json.dumps(module.to_dict(), indent=2)


def run_module(
    module: Module,
    hours: float = 48.0,
    dt: float = 0.05,
    context: dict[str, str] | None = None,
    initial: dict[str, float] | None = None,
    seed: int = 0,
    csv: str | None = None,
) -> tuple[str, Any]:
    """Run on the network runtime; returns (report text, trajectory)."""
    from genomeos.runtime import NetworkRuntime
    from genomeos.runtime.uncertainty import report_for_network

    vm = NetworkRuntime(module, context=context or {}, seed=seed)
    traj = vm.run(hours=hours, dt=dt, initial=initial or {})
    lines = [
        f"module {module.name}: {len(vm.active_rules)}/{len(module.rules)} rules active in context "
        f"{context or '{}'}; {hours} h simulated"
    ]
    for s in traj.species:
        xs = traj.levels[s]
        lines.append(f"  {s:<16} {spark(xs)}  final={xs[-1]:8.2f}  peaks={traj.peaks(s)}")
    lines.append(report_for_network(module, vm.active_rules).format())
    if csv:
        with open(csv, "w") as fh:
            fh.write("time," + ",".join(traj.species) + "\n")
            for i, t in enumerate(traj.times):
                fh.write(f"{t:.4f}," + ",".join(f"{traj.levels[s][i]:.6g}" for s in traj.species) + "\n")
        lines.append(f"  wrote {csv}")
    return "\n".join(lines), traj


def run_sbml(path: str, hours: float = 48.0, dt: float = 0.05) -> str:
    """An SBML model on the in-house engine (an engine module, so no application import)."""
    from genomeos.runtime.sbml import SbmlModel, SbmlRuntime

    model = SbmlModel.from_file(path)
    rt = SbmlRuntime(model)
    traj = rt.run(duration=hours, dt=dt, record_every=max(1, int(hours / dt / 600)))
    lines = [
        f"SBML {model.id} {model.name!r}: {len(model.species)} species, {len(model.reactions)} reactions, "
        f"{len(model.assignments)} assignment rules; {hours} time units simulated"
    ]
    for sp in traj.species:
        xs = traj.levels[sp]
        lines.append(f"  {sp:<16} {spark(xs)}  final={xs[-1]:10.3f}  peaks={traj.peaks(sp)}")
    return "\n".join(lines)


def run_boolean(path: str, init: dict[str, bool] | None = None, seed: int = 0) -> str:
    """A Boolean network (.bnet): one attractor from an initial state, or a sample of attractors."""
    from genomeos.runtime.boolean import BooleanNetwork

    net = BooleanNetwork.from_file(path)
    lines = [f"Boolean network {path}: {len(net.nodes)} nodes"]
    if init:
        att = net.attractor(init)
        kind = "fixed point" if len(att) == 1 else f"cycle of length {len(att)}"
        lines.append(f"  from {init}: {kind}")
        lines += ["    " + " ".join(f"{n}={int(v)}" for n, v in st.items()) for st in att]
    else:
        atts = net.attractors(samples=300, seed=seed)
        lines.append(f"  {len(atts)} distinct synchronous attractors from 300 random starts")
        for att in atts:
            on = sorted(n for n, v in att[0].items() if v)
            lines.append(f"    length {len(att)}: first state on={on}")
    return "\n".join(lines)
