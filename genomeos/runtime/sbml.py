"""SBML engine (task 2.3): run reaction-network models from BioModels.

Supports the SBML core used by most curated ODE models: compartments,
species (amounts or concentrations), global and local parameters,
assignment rules, reactions with kinetic laws in MathML, and a subset of
MathML (arithmetic, power, exp/ln/log, piecewise, comparisons, constants).
Integration is RK4. This is the zero-dependency stand-in for libRoadRunner;
a roadrunner adapter can replace `SbmlRuntime.run` without changing callers.
"""

from __future__ import annotations

import math
import xml.etree.ElementTree as ET
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

from genomeos.ir import Action, Entity, Evidence, EvidenceKind, Module, Rule

from .grn import Trajectory

MATHML = "http://www.w3.org/1998/Math/MathML"
Env = dict[str, float]
Fn = Callable[[Env], float]

_BINARY = {
    "plus": lambda a, b: a + b,
    "minus": lambda a, b: a - b,
    "times": lambda a, b: a * b,
    "divide": lambda a, b: a / b if b != 0 else math.inf,
    "power": lambda a, b: a**b,
    "lt": lambda a, b: float(a < b),
    "leq": lambda a, b: float(a <= b),
    "gt": lambda a, b: float(a > b),
    "geq": lambda a, b: float(a >= b),
    "eq": lambda a, b: float(a == b),
    "neq": lambda a, b: float(a != b),
}
_UNARY = {
    "exp": math.exp,
    "ln": lambda x: math.log(x) if x > 0 else -math.inf,
    "abs": abs,
    "floor": math.floor,
    "ceiling": math.ceil,
    "sin": math.sin,
    "cos": math.cos,
    "tan": math.tan,
    "not": lambda x: float(not x),
}
_CONST = {"pi": math.pi, "exponentiale": math.e, "true": 1.0, "false": 0.0}


def _tag(el: ET.Element) -> str:
    return el.tag.split("}", 1)[1] if "}" in el.tag else el.tag


def compile_math(el: ET.Element) -> Fn:
    """Compile a MathML element to a closure over an environment."""
    tag = _tag(el)
    if tag == "math":
        children = list(el)
        return compile_math(children[0]) if children else (lambda env: 0.0)
    if tag == "cn":
        text = (el.text or "").strip()
        sep = el.find(f"{{{MATHML}}}sep")
        if sep is not None:  # e-notation / rational
            mant = float(text)
            exp = float((sep.tail or "0").strip())
            val = mant * 10**exp if el.get("type") == "e-notation" else mant / exp
        else:
            val = float(text)
        return lambda env, v=val: v
    if tag == "ci":
        name = (el.text or "").strip()
        return lambda env, n=name: env[n]
    if tag in _CONST:
        return lambda env, v=_CONST[tag]: v
    if tag == "piecewise":
        pieces = []

        def otherwise(env: Env) -> float:
            return 0.0

        for child in el:
            if _tag(child) == "piece":
                val, cond = (compile_math(c) for c in list(child)[:2])
                pieces.append((cond, val))
            elif _tag(child) == "otherwise":
                otherwise = compile_math(list(child)[0])

        def pw(env: Env, pieces=pieces, otherwise=otherwise) -> float:
            for cond, val in pieces:
                if cond(env):
                    return val(env)
            return otherwise(env)

        return pw
    if tag == "apply":
        children = list(el)
        op = _tag(children[0])
        args = [compile_math(c) for c in children[1:]]
        if op in ("plus", "times") and len(args) != 2:
            if op == "plus":
                return lambda env, args=args: sum(a(env) for a in args)
            return lambda env, args=args: math.prod(a(env) for a in args)
        if op == "minus" and len(args) == 1:
            return lambda env, a=args[0]: -a(env)
        if op in _BINARY:
            f = _BINARY[op]
            return lambda env, f=f, a=args[0], b=args[1]: f(a(env), b(env))
        if op in _UNARY:
            f = _UNARY[op]
            return lambda env, f=f, a=args[0]: f(a(env))
        if op == "log":
            if _tag(children[1]) == "logbase":
                base = compile_math(list(children[1])[0])
                a = compile_math(children[2])
                return lambda env, base=base, a=a: math.log(a(env), base(env))
            return lambda env, a=args[0]: math.log10(a(env))
        if op == "root":
            if _tag(children[1]) == "degree":
                deg = compile_math(list(children[1])[0])
                a = compile_math(children[2])
                return lambda env, deg=deg, a=a: a(env) ** (1.0 / deg(env))
            return lambda env, a=args[0]: math.sqrt(a(env))
        if op == "and":
            return lambda env, args=args: float(all(a(env) for a in args))
        if op == "or":
            return lambda env, args=args: float(any(a(env) for a in args))
        raise ValueError(f"unsupported MathML operator {op!r}")
    raise ValueError(f"unsupported MathML element {tag!r}")


@dataclass(slots=True)
class Reaction:
    id: str
    reactants: list[tuple[str, float]]
    products: list[tuple[str, float]]
    modifiers: list[str]
    rate: Fn
    local_params: dict[str, float] = field(default_factory=dict)


@dataclass(slots=True)
class SbmlModel:
    id: str
    name: str
    compartments: dict[str, float]
    species: dict[str, float]  # initial values
    species_compartment: dict[str, str]
    boundary: set[str]
    constant_species: set[str]
    parameters: dict[str, float]
    assignments: list[tuple[str, Fn]]  # in dependency order
    reactions: list[Reaction]
    notes: str = ""

    @classmethod
    def from_file(cls, path: str | Path) -> SbmlModel:
        root = ET.parse(path).getroot()
        ns = root.tag.split("}")[0] + "}"
        model = root.find(f"{ns}model")
        if model is None:
            raise ValueError("no <model> element")

        def items(list_name: str) -> list[ET.Element]:
            lst = model.find(f"{ns}{list_name}")
            return list(lst) if lst is not None else []

        comps = {c.get("id"): float(c.get("size") or 1.0) for c in items("listOfCompartments")}
        species, species_comp, boundary, const = {}, {}, set(), set()
        for sp in items("listOfSpecies"):
            sid = sp.get("id")
            val = sp.get("initialAmount") or sp.get("initialConcentration") or "0"
            species[sid] = float(val)
            species_comp[sid] = sp.get("compartment", "")
            if sp.get("boundaryCondition") == "true":
                boundary.add(sid)
            if sp.get("constant") == "true":
                const.add(sid)
        params = {
            p.get("id"): float(p.get("value"))
            for p in items("listOfParameters")
            if p.get("value") not in (None, "")
        }
        for p in items("listOfParameters"):
            params.setdefault(p.get("id"), 0.0)

        raw_assign: dict[str, ET.Element] = {}
        for r in items("listOfRules"):
            if _tag(r) == "assignmentRule":
                raw_assign[r.get("variable")] = r.find(f"{{{MATHML}}}math")
        # order assignment rules by dependency
        deps = {
            v: {ci.text.strip() for ci in m.iter(f"{{{MATHML}}}ci")} & set(raw_assign)
            for v, m in raw_assign.items()
        }
        ordered: list[str] = []
        while len(ordered) < len(raw_assign):
            progress = False
            for v in raw_assign:
                if v not in ordered and deps[v] <= set(ordered):
                    ordered.append(v)
                    progress = True
            if not progress:
                raise ValueError("cyclic assignment rules")
        assignments = [(v, compile_math(raw_assign[v])) for v in ordered]

        reactions = []
        for r in items("listOfReactions"):

            def refs(name: str, r=r) -> list[tuple[str, float]]:
                lst = r.find(f"{ns}{name}")
                return (
                    []
                    if lst is None
                    else [(x.get("species"), float(x.get("stoichiometry") or 1.0)) for x in lst]
                )

            kl = r.find(f"{ns}kineticLaw")
            local = {}
            if kl is not None:
                lp = kl.find(f"{ns}listOfParameters") or kl.find(f"{ns}listOfLocalParameters")
                if lp is not None:
                    local = {p.get("id"): float(p.get("value")) for p in lp}
                rate = compile_math(kl.find(f"{{{MATHML}}}math"))
            else:
                rate = lambda env: 0.0  # noqa: E731
            reactions.append(
                Reaction(
                    r.get("id"),
                    refs("listOfReactants"),
                    refs("listOfProducts"),
                    [x[0] for x in refs("listOfModifiers")],
                    rate,
                    local,
                )
            )
        notes = model.find(f"{ns}notes")
        return cls(
            model.get("id", ""),
            model.get("name", ""),
            comps,
            species,
            species_comp,
            boundary,
            const,
            params,
            assignments,
            reactions,
            "".join(notes.itertext()).strip()[:500] if notes is not None else "",
        )

    # ---- BioIR ------------------------------------------------------------

    def to_module(self, source: str = "BioModels") -> Module:
        ev = Evidence(EvidenceKind.CURATED, f"{source}: {self.id} {self.name}".strip())
        m = Module(name=f"sbml.{self.id or 'model'}")
        for sid in self.species:
            m.add(
                Entity(
                    id=sid, kind="species", attrs={"initial": self.species[sid]}, evidence=ev, confidence=0.8
                )
            )
        for rx in self.reactions:
            for s, _ in rx.reactants:
                for p, _ in rx.products:
                    m.rules.append(
                        Rule(
                            id=f"{rx.id}:{s}->{p}",
                            source=s,
                            action=Action.MODIFY,
                            target=p,
                            evidence=ev,
                            confidence=0.8,
                        )
                    )
            for mod in rx.modifiers:
                for p, _ in rx.products:
                    m.rules.append(
                        Rule(
                            id=f"{rx.id}:{mod}~{p}",
                            source=mod,
                            action=Action.MODIFY,
                            target=p,
                            evidence=ev,
                            confidence=0.8,
                        )
                    )
        return m


class SbmlRuntime:
    def __init__(self, model: SbmlModel) -> None:
        self.model = model
        self.dynamic = [
            s for s in model.species if s not in model.boundary and s not in model.constant_species
        ]

    def _env(self, state: dict[str, float]) -> Env:
        env: Env = dict(self.model.parameters)
        env.update(self.model.compartments)
        env.update(state)
        for var, fn in self.model.assignments:
            env[var] = fn(env)
        return env

    def _derivatives(self, state: dict[str, float]) -> dict[str, float]:
        env = self._env(state)
        d = dict.fromkeys(self.dynamic, 0.0)
        for rx in self.model.reactions:
            local_env = {**env, **rx.local_params} if rx.local_params else env
            v = rx.rate(local_env)
            for s, st in rx.reactants:
                if s in d:
                    d[s] -= st * v
            for s, st in rx.products:
                if s in d:
                    d[s] += st * v
        return d

    def run(
        self,
        duration: float,
        dt: float = 0.01,
        record_every: int = 10,
        initial: dict[str, float] | None = None,
    ) -> Trajectory:
        state = dict(self.model.species)
        if initial:
            state.update(initial)
        steps = int(round(duration / dt))
        traj = Trajectory(times=[], species=list(self.dynamic), levels={s: [] for s in self.dynamic})

        def record(t: float) -> None:
            traj.times.append(t)
            for s in self.dynamic:
                traj.levels[s].append(state[s])

        record(0.0)
        for step in range(1, steps + 1):
            k1 = self._derivatives(state)
            s2 = {**state, **{s: state[s] + 0.5 * dt * k1[s] for s in self.dynamic}}
            k2 = self._derivatives(s2)
            s3 = {**state, **{s: state[s] + 0.5 * dt * k2[s] for s in self.dynamic}}
            k3 = self._derivatives(s3)
            s4 = {**state, **{s: state[s] + dt * k3[s] for s in self.dynamic}}
            k4 = self._derivatives(s4)
            for s in self.dynamic:
                state[s] = max(0.0, state[s] + dt / 6 * (k1[s] + 2 * k2[s] + 2 * k3[s] + k4[s]))
            if step % record_every == 0:
                record(step * dt)
        return traj
