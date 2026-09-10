"""Gene regulatory network runtime.

Executes a BioIR Module as a continuous-time system: mRNA and protein levels
change according to PRODUCE / ACTIVATE / INHIBIT rules. Regulation uses Hill
functions, the standard phenomenological form for transcription-factor binding.

    transcription(gene) = basal + max_rate * A(gene) * R(gene)

    A = mean over activators of  s_i * H(x_i)         (1 if no activators)
    R = product over inhibitors of (1 - s_j * H(x_j))
    H(x) = x^n / (K^n + x^n)

    d[mRNA]/dt    = transcription - delta_m * [mRNA]
    d[protein]/dt = k_tl * [mRNA] - delta_p * [protein]

Integration: RK4 when noise == 0, Euler-Maruyama otherwise. Everything is in
arbitrary concentration units and hours; parameters in the module override the
defaults below and carry their own evidence.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, field

from genomeos.ir import UNKNOWN, Action, Gene, Module, Protein, Rule

LN2 = math.log(2)

DEFAULTS = {
    "translation_rate": 5.0,  # protein per mRNA per hour
    "mrna_half_life": LN2,  # hours  (delta_m = 1 / hour)
    "protein_half_life": LN2 / 5,  # hours (delta_p = 5 / hour) unless protein declares one
    "noise": 0.0,  # multiplicative noise amplitude
}


@dataclass(slots=True)
class Trajectory:
    times: list[float]
    species: list[str]
    levels: dict[str, list[float]] = field(default_factory=dict)

    def final(self) -> dict[str, float]:
        return {s: self.levels[s][-1] for s in self.species}

    def peaks(self, species: str) -> int:
        """Count local maxima; a crude oscillation detector."""
        xs = self.levels[species]
        return sum(1 for i in range(1, len(xs) - 1) if xs[i - 1] < xs[i] >= xs[i + 1])


def _hill(x: float, k: float, n: float) -> float:
    if x <= 0:
        return 0.0
    xn = x**n
    return xn / (k**n + xn)


class NetworkRuntime:
    def __init__(
        self, module: Module, context: dict[str, str] | None = None, seed: int | None = None
    ) -> None:
        self.module = module
        self.context = context or {}
        self.rng = random.Random(seed)
        self.params = dict(DEFAULTS)
        for name, p in module.parameters.items():
            if name in self.params:
                self.params[name] = p.value

        self.genes: list[Gene] = module.genes()
        self.proteins: list[Protein] = module.proteins()
        self.active_rules: list[Rule] = [r for r in module.rules if r.applies(self.context)]

        # index regulators per gene, translation per protein
        self.activators: dict[str, list[Rule]] = {g.id: [] for g in self.genes}
        self.inhibitors: dict[str, list[Rule]] = {g.id: [] for g in self.genes}
        self.produces: dict[str, list[str]] = {p.id: [] for p in self.proteins}  # protein <- genes
        for r in self.active_rules:
            if r.action is Action.ACTIVATE and r.target in self.activators:
                self.activators[r.target].append(r)
            elif r.action is Action.INHIBIT and r.target in self.inhibitors:
                self.inhibitors[r.target].append(r)
            elif r.action is Action.PRODUCE and r.target in self.produces:
                self.produces[r.target].append(r.source)

        self.species = [f"{g.id}.mRNA" for g in self.genes] + [p.id for p in self.proteins]

    # ---- dynamics --------------------------------------------------------

    def _derivatives(self, state: dict[str, float]) -> dict[str, float]:
        d: dict[str, float] = {}
        delta_m = LN2 / self.params["mrna_half_life"]
        k_tl = self.params["translation_rate"]

        for g in self.genes:
            acts = self.activators[g.id]
            inhs = self.inhibitors[g.id]
            a = 1.0
            if acts:
                a = sum(
                    r.strength * _hill(state.get(r.source, 0.0), r.threshold, r.hill) for r in acts
                ) / len(acts)
            rep = 1.0
            for r in inhs:
                rep *= 1.0 - r.strength * _hill(state.get(r.source, 0.0), r.threshold, r.hill)
            max_rate = float(g.attrs.get("max_rate", 0.0))
            key = f"{g.id}.mRNA"
            d[key] = g.basal_rate + max_rate * a * rep - delta_m * state[key]

        for p in self.proteins:
            hl = p.half_life_h if p.half_life_h is not UNKNOWN else self.params["protein_half_life"]
            delta_p = LN2 / float(hl)
            production = sum(k_tl * state[f"{gid}.mRNA"] for gid in self.produces[p.id])
            d[p.id] = production - delta_p * state[p.id]
        return d

    def run(
        self,
        hours: float,
        dt: float = 0.01,
        initial: dict[str, float] | None = None,
        record_every: int = 10,
        clamp: dict[str, float] | None = None,
    ) -> Trajectory:
        """Integrate for `hours`. `clamp` holds species at fixed levels (external
        signals such as a morphogen the module does not itself produce)."""
        state = dict.fromkeys(self.species, 0.0)
        if initial:
            state.update(initial)
        if clamp:
            state.update(clamp)
        noise = self.params["noise"]
        steps = int(round(hours / dt))
        traj = Trajectory(times=[], species=list(self.species), levels={s: [] for s in self.species})

        def record(t: float) -> None:
            traj.times.append(t)
            for s in self.species:
                traj.levels[s].append(state[s])

        record(0.0)
        for step in range(1, steps + 1):
            if noise > 0:
                d = self._derivatives(state)
                for s in self.species:
                    state[s] = max(
                        0.0,
                        state[s]
                        + d[s] * dt
                        + noise * math.sqrt(dt) * self.rng.gauss(0, 1) * math.sqrt(max(state[s], 1e-9)),
                    )
            else:
                k1 = self._derivatives(state)
                s2 = {s: state[s] + 0.5 * dt * k1[s] for s in self.species}
                k2 = self._derivatives(s2)
                s3 = {s: state[s] + 0.5 * dt * k2[s] for s in self.species}
                k3 = self._derivatives(s3)
                s4 = {s: state[s] + dt * k3[s] for s in self.species}
                k4 = self._derivatives(s4)
                for s in self.species:
                    state[s] = max(0.0, state[s] + dt / 6 * (k1[s] + 2 * k2[s] + 2 * k3[s] + k4[s]))
            if clamp:
                state.update(clamp)
            if step % record_every == 0:
                record(step * dt)
        if steps % record_every != 0:  # always record the final state
            record(steps * dt)
        return traj
