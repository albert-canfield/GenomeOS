# SPDX-License-Identifier: Apache-2.0
# Part of the BioLang engine (language, IR, VM, standard library); see LICENSING.md.
"""The located runtime: the central dogma with places (BioLang v0.4 §4, stage 1).

Every species is a pair of a molecule and a compartment, so `SDHB.mRNA@Nucleus`,
`SDHB.mRNA@Cytosol`, `SDHBp@Cytosol` and `SDHBp@Mitochondrion` are four state variables. A gene
is transcribed where it is read, its mRNA is translated where it meets ribosomes, and a protein
only arrives somewhere if a transport carries it: transports are saturable and shared, so cargos
compete for one finite capacity.

    d[mRNA@c]/dt    = transcription(c) + imports - exports - delta_m [mRNA@c]
    d[P@c]/dt       = translation(c) + imports - exports - delta_p [P@c]
    d[complex@c]/dt = assembly_rate * min(subunit levels at c) - delta_p [complex@c]
                      (and one of each subunit is used up per complex: nothing is made from nothing)

    flux of cargo i through transport t = capacity * (x_i/K) / (1 + sum_j x_j/K)

What the runtime refuses to do quietly: a species with no place, a rule across compartments and a
protein that reaches a compartment its declaration does not name are reported, never ignored. The
regime (continuous or stochastic per species, the update scheme, the allocation policy, the seed)
is recorded in the result, because a result without its regime is not a result.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, field

from genomeos.ir import UNKNOWN, Action, Module, Regime, molecules_in, to_molar
from genomeos.lang.located import layout, rule_sites

from .uncertainty import KIND_WEIGHT

LN2 = math.log(2)

DEFAULTS = {
    "translation_rate": 5.0,  # protein per mRNA per hour
    "mrna_half_life": LN2,  # hours
    "protein_half_life": LN2 / 5,  # hours, unless the protein declares one
    "assembly_rate": 1.0,  # complex per hour per unit of the limiting subunit
}


@dataclass(slots=True)
class LocatedResult:
    """Levels per molecule and compartment, what went where, how sure it is, and how it was run."""

    hours: float
    times: list[float] = field(default_factory=list)
    levels: dict[str, list[float]] = field(default_factory=dict)  # "P@Compartment" -> trajectory
    stranded: list[dict] = field(default_factory=list)  # declared somewhere it never arrived
    transit: list[dict] = field(default_factory=list)  # present on the route to a declared place
    ectopic: list[dict] = field(default_factory=list)  # present where nothing declared or leads on
    transports_used: dict[str, float] = field(default_factory=dict)  # total cargo moved per transport
    confidence: dict[str, dict] = field(default_factory=dict)  # species -> {score, weakest}
    regime: dict = field(default_factory=dict)

    def final(self, species: str) -> float:
        xs = self.levels.get(species)
        return xs[-1] if xs else 0.0

    def where(self, molecule: str, threshold: float = 1e-6) -> list[str]:
        """Compartments where a molecule ends the run above `threshold`."""
        pre = molecule + "@"
        return sorted(
            s[len(pre) :] for s, xs in self.levels.items() if s.startswith(pre) and xs[-1] > threshold
        )


def _hill(x: float, k: float, n: float) -> float:
    if x <= 0:
        return 0.0
    xn = x**n
    return xn / (k**n + xn)


class LocatedRuntime:
    """Run a located module. `knockouts` may name genes (no transcription), proteins (no translation),
    transports (closed), chromosomes such as `chrM` (every gene read from it: a rho0 cell) or a
    signal as `PROTEIN:signal` (the protein loses that targeting signal)."""

    def __init__(
        self,
        module: Module,
        context: dict[str, str] | None = None,
        knockouts: set[str] | frozenset[str] = frozenset(),
        regime: Regime | None = None,
    ) -> None:
        if not module.located:
            raise ValueError(f"module {module.name} declares no compartment; use the network runtime")
        self.module = module
        self.layout = layout(module)
        if self.layout.errors:
            raise ValueError("located program is invalid: " + "; ".join(self.layout.errors[:3]))
        self.regime = regime or module.regime or Regime()
        self._check_regime()
        self.knockouts = set(knockouts)
        self.params = dict(DEFAULTS)
        for name, p in module.parameters.items():
            if name in self.params:
                self.params[name] = p.value
        self.proteins = {p.id: p for p in module.proteins()}
        self.genes = {g.id: g for g in module.genes()}
        self.rules = [r for r in module.rules if r.applies(context or {})]
        self._knock_out()
        self._build()

    def _check_regime(self) -> None:
        rg = self.regime
        if rg.update not in ("continuous",):
            raise ValueError(
                f"regime update {rg.update!r} is not implemented by the located runtime (stage 1)"
            )
        if rg.allocation != "competitive":
            raise ValueError(
                f"allocation {rg.allocation!r} arrives with pools (stage 2); transports are competitive"
            )
        if rg.treatment != "continuous" and rg.units != "copies":
            raise ValueError(
                "stochastic treatment needs `units: copies`; arbitrary units have no copy number"
            )

    def _knock_out(self) -> None:
        chroms = {c for c in self.knockouts if c.startswith("chr") or c == "MT"}
        signals = [k for k in self.knockouts if ":" in k]
        known = set(self.genes) | set(self.proteins) | {t.id for t in self.module.transports()}
        unknown = sorted(k for k in self.knockouts if k not in known and k not in chroms and k not in signals)
        if unknown:
            raise ValueError(f"knockouts name nothing in the program: {unknown}")
        if signals:  # a protein without its signal: re-resolve the routes on a copy of the program
            module = Module.from_dict(self.module.to_dict())
            for k in signals:
                pid, _, sig = k.partition(":")
                if pid not in module.entities or sig not in getattr(module.entities[pid], "signals", []):
                    raise ValueError(f"knockout {k!r}: {pid} carries no signal {sig!r}")
                module.entities[pid].signals.remove(sig)
            self.module = module
            self.layout = layout(module)
            self.proteins = {p.id: p for p in module.proteins()}
            self.genes = {g.id: g for g in module.genes()}
        pl = self.layout.places
        self.genes_off = {g for g in self.genes if g in self.knockouts}
        for gid, g in self.genes.items():
            site = pl.by_id[self.layout.gene_site[gid]]
            chrom = g.locus.chrom if g.locus is not None else ""
            if chrom in chroms or (not chrom and any(site.genome == [c] for c in chroms)):
                self.genes_off.add(gid)
        self.proteins_off = {p for p in self.proteins if p in self.knockouts}
        self.transports_off = {t.id for t in self.module.transports() if t.id in self.knockouts}

    def _build(self) -> None:
        lay, pl = self.layout, self.layout.places
        self.mrna = [(g, c) for g in self.genes for c in lay.mrna_sites.get(g, ())]
        self.protein_sites = {
            p.id: sorted(lay.reach.get(p.id, set()) | {c for c in p.location if c in pl.by_id})
            for p in self.proteins.values()
        }
        self.species = [f"{g}.mRNA@{c}" for g, c in self.mrna]
        self.species += [f"{p}@{c}" for p, cs in self.protein_sites.items() for c in cs]
        self.transcription = []
        for g, c in self.mrna:
            if c != lay.gene_site[g]:
                continue  # mRNA only appears where the gene is read; elsewhere it arrives by transport
            acts = [r for r in self.rules if r.action is Action.ACTIVATE and r.target == g]
            inhs = [r for r in self.rules if r.action is Action.INHIBIT and r.target == g]
            keep = lambda r, site=c: site in rule_sites(self.module, lay, r, self.proteins)  # noqa: E731
            self.transcription.append((g, c, [r for r in acts if keep(r)], [r for r in inhs if keep(r)]))
        self.thresholds = self._thresholds()
        self.translation = []
        for r in self.rules:
            if r.action is not Action.PRODUCE or r.source not in self.genes or r.target not in self.proteins:
                continue
            sites = [c for c in lay.mrna_sites.get(r.source, ()) if pl.by_id[c].translation]
            if sites and r.target not in self.proteins_off:
                self.translation.append((r.source, sites[0], r.target))
        have = set(self.species)
        self.moves = []  # (transport, [(from species, to species)])
        for t in self.module.transports():
            pairs = [
                (f"{g}.mRNA@{t.from_compartment}", f"{g}.mRNA@{t.to_compartment}")
                for g in self.genes
                if t.carries(g, mrna=True)
            ]
            pairs += [
                (f"{p.id}@{t.from_compartment}", f"{p.id}@{t.to_compartment}")
                for p in self.proteins.values()
                if t.carries(p.id, signals=tuple(p.signals))
            ]
            pairs = [(a, b) for a, b in pairs if a in have and b in have]
            if pairs:
                self.moves.append((t, pairs))
        self.assembly = []  # (complex species, [(subunit, site)])
        by_complex: dict[str, list] = {}
        for r in self.rules:
            if r.action is Action.BIND and r.target in self.proteins and r.source in self.proteins:
                for c in rule_sites(self.module, lay, r, self.proteins):
                    by_complex.setdefault(f"{r.target}@{c}", []).append((r.source, c))
        self.assembly = list(by_complex.items())
        self.delta_m = LN2 / self.params["mrna_half_life"]
        self.delta_p = {}
        for p in self.proteins.values():
            hl = p.half_life_h if p.half_life_h is not UNKNOWN else self.params["protein_half_life"]
            self.delta_p[p.id] = LN2 / float(hl)

    def _thresholds(self) -> dict[tuple[str, str], float]:
        """Each regulatory rule's threshold as an amount at each compartment it acts in.

        A threshold with no unit is already an amount and passes through unchanged, which is every
        program committed before 2026-09-19. A threshold stated as a concentration is divided by the
        compartment's absolute volume here and nowhere else, so the same rule in a big cell and a
        small one is a different number of molecules - which is the whole point of the field (§4.1)."""
        out: dict[tuple[str, str], float] = {}
        for _g, c, acts, inhs in self.transcription:
            for r in (*acts, *inhs):
                if not r.threshold_unit:
                    out[(r.id, c)] = r.threshold
                    continue
                if self.regime.units != "copies":
                    raise ValueError(
                        f"rule {r.id!r} states a threshold in {r.threshold_unit}, a concentration, but "
                        f"the regime declares units {self.regime.units!r}: molecules over a volume are "
                        "only comparable with a counted amount, so state `regime { units: copies }`"
                    )
                volume = self.layout.places.by_id[c].absolute_volume_fl
                if volume is UNKNOWN:  # also a compile error; here for a module built without the parser
                    raise ValueError(
                        f"rule {r.id!r} acts in {c!r}, which declares no absolute_volume, so its "
                        f"threshold of {r.threshold} {r.threshold_unit} has nothing to divide by"
                    )
                out[(r.id, c)] = molecules_in(to_molar(r.threshold, r.threshold_unit), float(volume))
        return out

    # ---- dynamics --------------------------------------------------------------------

    def _read(self, state: dict[str, float], molecule: str, site: str) -> float:
        """A molecule's level as read at a site: there, or from a membrane facing it."""
        here = state.get(f"{molecule}@{site}")
        if here is not None:
            return here
        p = self.proteins.get(molecule)
        if p is None:
            return 0.0
        pl = self.layout.places
        return sum(
            state.get(f"{molecule}@{c}", 0.0) for c in p.location if c in pl.by_id and site in pl.faces(c)
        )

    def _stock(self, state: dict[str, float], molecule: str, site: str) -> str:
        """The species a reaction at `site` draws a molecule from: there, or the membrane facing it."""
        here = f"{molecule}@{site}"
        if here in state:
            return here
        pl = self.layout.places
        for c in self.proteins[molecule].location:
            if c in pl.by_id and site in pl.faces(c) and f"{molecule}@{c}" in state:
                return f"{molecule}@{c}"
        return here

    def fluxes(self, state: dict[str, float]) -> list[tuple[float, tuple, str]]:
        """Every reaction as (rate, ((species, change), ...), label); both treatments use these."""
        out: list[tuple[float, tuple, str]] = []
        for g, c, acts, inhs in self.transcription:
            if g in self.genes_off:
                continue
            gene = self.genes[g]
            a = 1.0
            k = self.thresholds
            if acts:
                a = sum(
                    r.strength * _hill(self._read(state, r.source, c), k[(r.id, c)], r.hill) for r in acts
                )
                a /= len(acts)
            rep = 1.0
            for r in inhs:
                rep *= 1.0 - r.strength * _hill(self._read(state, r.source, c), k[(r.id, c)], r.hill)
            rate = gene.basal_rate + float(gene.attrs.get("max_rate", 0.0)) * a * rep
            out.append((rate, ((f"{g}.mRNA@{c}", 1),), "transcription"))
        for g, c, p in self.translation:
            rate = self.params["translation_rate"] * state[f"{g}.mRNA@{c}"]
            out.append((rate, ((f"{p}@{c}", 1),), "translation"))
        for t, pairs in self.moves:
            if t.id in self.transports_off:
                continue
            cap = t.capacity
            if t.via:
                v = max(
                    self._read(state, t.via, t.from_compartment), self._read(state, t.via, t.to_compartment)
                )
                cap *= v / (t.via_threshold + v) if v > 0 else 0.0
            load = sum(state[a] for a, _ in pairs) / t.affinity
            if cap <= 0 or load <= 0:
                continue
            for a, b in pairs:
                rate = cap * (state[a] / t.affinity) / (1.0 + load)
                if rate > 0:
                    out.append((rate, ((a, -1), (b, 1)), t.id))
        for sp, subunits in self.assembly:
            limiting = min(self._read(state, sub, site) for sub, site in subunits)
            if limiting > 0:  # one of each subunit is used up per complex made (stoichiometry: stage 2)
                used = tuple((self._stock(state, sub, site), -1) for sub, site in subunits)
                out.append((self.params["assembly_rate"] * limiting, ((sp, 1), *used), "assembly"))
        for sp in self.species:
            x = state[sp]
            if x <= 0:
                continue
            molecule = sp.split("@", 1)[0]
            delta = self.delta_m if molecule.endswith(".mRNA") else self.delta_p[molecule]
            out.append((delta * x, ((sp, -1),), "decay"))
        return out

    def _derivatives(self, state: dict[str, float]) -> dict[str, float]:
        d = dict.fromkeys(self.species, 0.0)
        for rate, changes, _ in self.fluxes(state):
            for sp, change in changes:
                d[sp] += change * rate
        return d

    def _treatment_of(self, state: dict[str, float]) -> dict[str, str]:
        """Which species run stochastically: as declared, or under `auto` those below the threshold."""
        mode = self.regime.treatment
        if mode == "continuous":
            return dict.fromkeys(self.species, "continuous")
        if mode == "stochastic":
            return dict.fromkeys(self.species, "stochastic")
        return {s: ("stochastic" if state[s] < self.regime.threshold else "continuous") for s in self.species}

    def _poisson(self, lam: float) -> float:
        if lam <= 0:
            return 0.0
        if lam > 50:  # normal approximation, clipped at zero
            return max(0.0, round(self.rng.gauss(lam, math.sqrt(lam))))
        limit, k, prod = math.exp(-lam), 0, self.rng.random()
        while prod > limit:
            k += 1
            prod *= self.rng.random()
        return float(k)

    def _step_leap(self, state: dict[str, float], dt: float, moved: dict[str, float]) -> dict[str, int]:
        """One tau-leap step: reactions touching a stochastic species fire a Poisson number of times."""
        treat = self._treatment_of(state)
        counts = {"stochastic": 0, "continuous": 0}
        delta = dict.fromkeys(self.species, 0.0)
        for rate, changes, label in self.fluxes(state):
            noisy = any(treat[sp] == "stochastic" for sp, _ in changes)
            n = self._poisson(rate * dt) if noisy else rate * dt
            counts["stochastic" if noisy else "continuous"] += 1
            for sp, change in changes:
                delta[sp] += change * n
            if label in moved:
                moved[label] += n
        for s in self.species:
            state[s] = max(0.0, state[s] + delta[s])
        return counts

    def run(
        self, hours: float, dt: float = 0.05, initial: dict[str, float] | None = None, record_every: int = 10
    ) -> LocatedResult:
        """Integrate for `hours`. `initial` sets species levels by `molecule@compartment`."""
        state = dict.fromkeys(self.species, 0.0)
        for p in self.proteins.values():  # the bootstrap state the program declares
            for c in p.location:
                if p.initial and f"{p.id}@{c}" in state:
                    state[f"{p.id}@{c}"] = p.initial
        for k, v in (initial or {}).items():
            if k not in state:
                raise ValueError(f"initial level for {k!r}, which is not a species of this program")
            state[k] = float(v)
        self.rng = random.Random(self.regime.seed)
        res = LocatedResult(hours=hours, levels={s: [] for s in self.species})
        moved = {t.id: 0.0 for t, _ in self.moves}
        reactions = {"stochastic": 0, "continuous": 0}
        steps = int(round(hours / dt))

        def record(t: float) -> None:
            res.times.append(t)
            for s in self.species:
                res.levels[s].append(state[s])

        record(0.0)
        for step in range(1, steps + 1):
            if self.regime.treatment == "continuous":
                self._step_rk4(state, dt, moved)
            else:
                for k, n in self._step_leap(state, dt, moved).items():
                    reactions[k] += n
            if step % record_every == 0 or step == steps:
                record(step * dt)
        res.transports_used = {k: round(v, 6) for k, v in moved.items()}
        res.regime = self._regime_record(state, dt, reactions)
        self._diagnose(res)
        return res

    def _stage(self, state: dict[str, float], moved_rate: dict[str, float]) -> dict[str, float]:
        d = dict.fromkeys(self.species, 0.0)
        for rate, changes, label in self.fluxes(state):
            for sp, change in changes:
                d[sp] += change * rate
            if label in moved_rate:
                moved_rate[label] += rate
        return d

    def _step_rk4(self, state: dict[str, float], dt: float, moved: dict[str, float]) -> None:
        rates = [dict.fromkeys(moved, 0.0) for _ in range(4)]
        k1 = self._stage(state, rates[0])
        k2 = self._stage({s: state[s] + 0.5 * dt * k1[s] for s in self.species}, rates[1])
        k3 = self._stage({s: state[s] + 0.5 * dt * k2[s] for s in self.species}, rates[2])
        k4 = self._stage({s: state[s] + dt * k3[s] for s in self.species}, rates[3])
        for s in self.species:
            state[s] = max(0.0, state[s] + dt / 6 * (k1[s] + 2 * k2[s] + 2 * k3[s] + k4[s]))
        for label in moved:
            moved[label] += (
                dt / 6 * (rates[0][label] + 2 * rates[1][label] + 2 * rates[2][label] + rates[3][label])
            )

    def _regime_record(self, state: dict[str, float], dt: float, reactions: dict[str, int]) -> dict:
        rg = self.regime
        final = self._treatment_of(state)
        return {
            "name": rg.name,
            "declared": rg.name != "default" or self.module.regime is not None,
            "treatment": rg.treatment,
            "threshold": rg.threshold,
            "units": rg.units,
            "update": rg.update,
            "integrator": "RK4" if rg.treatment == "continuous" else "tau-leap (Poisson per reaction)",
            "dt_hours": dt,
            "allocation": rg.allocation,
            "seed": rg.seed,
            "species_stochastic_at_end": sorted(s for s, v in final.items() if v == "stochastic"),
            "reaction_firings": reactions if rg.treatment != "continuous" else {},
        }

    # ---- what went where, and how sure it is ----------------------------------------------

    def _route(self, starts: list[str], goal: str, allowed) -> list | None:
        """The transports on the shortest open route from any start to `goal` (None if there is none)."""
        prev: dict[str, tuple] = {s: () for s in starts}
        queue = list(starts)
        while queue:
            here = queue.pop(0)
            if here == goal:
                path, node = [], here
                while prev[node]:
                    t, node = prev[node]
                    path.append(t)
                return path[::-1]
            for t in self.module.transports():
                if t.id in self.transports_off or t.from_compartment != here or not allowed(t):
                    continue
                if t.to_compartment not in prev:
                    prev[t.to_compartment] = (t, here)
                    queue.append(t.to_compartment)
        return None

    def _diagnose(self, res: LocatedResult, eps: float = 1e-6) -> None:
        for p in self.proteins.values():
            sites = self.protein_sites[p.id]
            levels = {c: res.final(f"{p.id}@{c}") for c in sites}
            if not any(v > eps for v in levels.values()):
                continue  # not made at all: absent, not mislocalised
            carries = lambda t, pid=p.id, sig=tuple(p.signals): t.carries(pid, signals=sig)  # noqa: E731
            for c in p.location:
                if (
                    max(res.levels.get(f"{p.id}@{c}", [0.0])) <= eps
                ):  # never arrived (depletion is not mislocalisation)
                    held = sorted(k for k, v in levels.items() if v > eps)
                    res.stranded.append(
                        {"protein": p.id, "declared": c, "found_in": held, "signals": list(p.signals)}
                    )
            for c, v in levels.items():
                if v <= eps or c in p.location:
                    continue
                onward = any(self._route([c], d, carries) for d in p.location)
                entry = {
                    "protein": p.id,
                    "compartment": c,
                    "level": round(v, 6),
                    "declared": list(p.location),
                }
                (res.transit if onward else res.ectopic).append(entry)
        self._confidence(res)

    def _confidence(self, res: LocatedResult) -> None:
        """Each protein level is as grounded as the weakest fact on the chain that put it there."""
        by_id = self.module.entities
        rules = {(r.source, r.target): r for r in self.rules if r.action is Action.PRODUCE}
        for p in self.proteins.values():
            gid, produces = "", None
            for (src, tgt), r in rules.items():
                if tgt == p.id:
                    gid, produces = src, r
                    break
            start = self.layout.synthesis.get(p.id) or list(p.location)
            carries = lambda t, pid=p.id, sig=tuple(p.signals): t.carries(pid, signals=sig)  # noqa: E731
            for c in self.protein_sites[p.id]:
                chain = [(f"compartment {c}", by_id[c]), (f"protein {p.id}", p)]
                if produces is not None:
                    chain += [(f"gene {gid}", self.genes[gid]), (f"rule {produces.id}", produces)]
                    mrna = self._route(
                        [self.layout.gene_site[gid]], start[0], lambda t, g=gid: t.carries(g, mrna=True)
                    )
                    chain += [(f"transport {t.id}", t) for t in mrna or []]
                chain += [(f"transport {t.id}", t) for t in self._route(start, c, carries) or []]
                scored = [(n, o.confidence * KIND_WEIGHT[o.evidence.kind]) for n, o in chain]
                if self._route(start, c, carries) is None:
                    scored.append((f"no open route to {c}", 0.0))
                weakest = min(scored, key=lambda kv: kv[1])
                res.confidence[f"{p.id}@{c}"] = {
                    "score": round(weakest[1], 3),
                    "weakest": weakest[0],
                    "chain": [n for n, _ in scored],
                }
