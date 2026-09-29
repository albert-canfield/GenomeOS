# SPDX-License-Identifier: AGPL-3.0-or-later
"""Network knockouts as an experiment: hold nodes or edges of a network at zero and read the consequence.

The `experiment` block perturbs organisms; a network (a Boolean model or a rule module run by the GRN
runtime) had knockouts only inside hand-written tests. Here an experiment names nodes (a species held
at zero for the whole run) and edges (`A>B`: B reads A as zero, everything else reads A as it is), the
perturbed run is compared against the unperturbed one from the same start, and the size of the change
is ranked against knockouts of the same size drawn at random, so a consequence is read against chance.
"""

from __future__ import annotations

import ast
import copy
import itertools
import math
import random
from dataclasses import dataclass, field

from genomeos.ir import Action, Module
from genomeos.runtime.boolean import BooleanNetwork, _py
from genomeos.runtime.grn import NetworkRuntime

# Registered 2026-09-27, before the code below was written.
EDGE_SEPARATOR = ">"  # "A>B": the edge from A into B's rule
CONTROL_DRAWS = 200  # random matched knockouts; all of them when fewer matched sets exist
CONTROL_SEED = 0
# consequence: mean over the nodes knocked out in NEITHER run of |fraction of the attractor (Boolean) or
# of the second half of the run (continuous) spent ON / the mean level, knockout minus unperturbed|, on
# a 0..1 scale for Boolean and as |log2((ko + 1e-3) / (wt + 1e-3))| for continuous levels; plus whether
# the dynamics changed kind (fixed point <-> cycle; oscillating <-> not).
# chance: p = (1 + #matched random knockouts with a consequence >= the observed) / (1 + #draws). Matched
# = the same number of nodes and the same number of edges, drawn from nodes and edges not chosen.
NETWORK_KNOCKOUT_SPEC = (
    "run(network, knockouts) returns the unperturbed and the perturbed dynamics from the same start, "
    "the per-node change, whether the dynamics changed kind, and the consequence ranked against "
    "matched random knockouts that never include the chosen nodes or edges."
)
NETWORK_KNOCKOUT_FALSIFIER = (
    "A knockout of a node no rule reads must have consequence 0 and chance p = 1.0; if the control "
    "calls it exceptional, or ever draws a chosen node, the control is broken. And the known case "
    "below must hold."
)
# the known case, from an existing fixture: data/models/mammalian_cell_cycle.bnet (Faure, Naldi,
# Chaouiya & Thieffry 2006, Bioinformatics 22:e124) and tests/test_boolean.py, which pins that with
# CycD (growth factor) on the network cycles and with CycD off it rests in a G1 fixed point with Rb,
# p27 and Cdh1 on and CycA, CycB off.
NETWORK_KNOCKOUT_KNOWN_CASE = {
    "network": "data/models/mammalian_cell_cycle.bnet",
    "start": "CycD on, every other node off (the fixture's start)",
    "knockout": "CycD",
    "unperturbed": "a cyclic attractor (length > 1)",
    "perturbed": "a fixed point with Rb, p27, Cdh1 on and CycA, CycB off; kind changed",
    "source": "Faure et al. 2006, as pinned by tests/test_boolean.py::test_faure_mammalian_cell_cycle",
}

# ---------------------------------------------------------------------------
# The code that honours the registration above.

Edge = tuple[str, str]
_LEVEL_FLOOR = 1e-3  # the registered floor inside the continuous log2 ratio


def parse_knockouts(spec: list[str] | tuple[str, ...]) -> tuple[list[str], list[Edge]]:
    """Split a knockout list into nodes and `A>B` edges, in the order given, without duplicates."""
    nodes: list[str] = []
    edges: list[Edge] = []
    for raw in spec:
        item = raw.strip()
        if not item:
            continue
        if EDGE_SEPARATOR in item:
            src, _, dst = item.partition(EDGE_SEPARATOR)
            edge = (src.strip(), dst.strip())
            if edge not in edges:
                edges.append(edge)
        elif item not in nodes:
            nodes.append(item)
    return nodes, edges


class _Zero(ast.NodeTransformer):
    """Rewrite the named inputs of a Boolean expression to False: B reads A as zero."""

    def __init__(self, names: set[str]) -> None:
        self.names = names

    def visit_Name(self, node: ast.Name) -> ast.AST:  # noqa: N802  (ast visitor name)
        return ast.copy_location(ast.Constant(False), node) if node.id in self.names else node


class _BooleanEngine:
    """A Boolean network: a node held at zero has the rule `False`; an edge reads its source as False."""

    kind_words = ("fixed point", "cycle")

    def __init__(self, net: BooleanNetwork, start: dict[str, bool] | None = None) -> None:
        self.net = net
        self.start = {n: bool((start or {}).get(n, False)) for n in net.nodes}

    @property
    def nodes(self) -> list[str]:
        return self.net.nodes

    @property
    def edges(self) -> list[Edge]:
        return [(src, node) for node in self.net.nodes for src in sorted(self.net.inputs_of(node))]

    def _perturbed(self, nodes_off: list[str], edges_off: list[Edge]) -> BooleanNetwork:
        off = set(nodes_off)
        reads_zero: dict[str, set[str]] = {}
        for src, dst in edges_off:
            reads_zero.setdefault(dst, set()).add(src)
        out = BooleanNetwork()
        for node, expr in self.net.rules.items():
            if node in off:
                out.add(node, "False")
                continue
            if node in reads_zero:
                tree = _Zero(reads_zero[node]).visit(ast.parse(_py(expr), mode="eval"))
                expr = ast.unparse(ast.fix_missing_locations(tree))
            out.add(node, expr)
        return out

    def measure(self, nodes_off: list[str], edges_off: list[Edge]) -> tuple[dict[str, float], str]:
        """Per node, the fraction of the attractor it spends on; the kind is fixed point or cycle."""
        net = self._perturbed(nodes_off, edges_off)
        start = dict(self.start)
        for n in nodes_off:
            start[n] = False
        att = net.attractor(start)
        levels = {n: sum(1 for s in att if s[n]) / len(att) for n in net.nodes}
        return levels, "fixed point" if len(att) == 1 else f"cycle of {len(att)}"

    def change(self, wt: float, ko: float) -> float:
        return abs(ko - wt)  # both are fractions of the attractor: already on a 0..1 scale


class _ContinuousEngine:
    """A rule module under the GRN runtime: a node is held at zero, an edge's target reads it as zero."""

    def __init__(
        self,
        module: Module,
        start: dict[str, float] | None = None,
        hours: float = 20.0,
        dt: float = 0.02,
        record_every: int = 5,
        seed: int | None = 0,
    ) -> None:
        self.module = module
        self.start = dict(start or {})
        self.hours, self.dt, self.record_every, self.seed = hours, dt, record_every, seed
        self._species = {}
        for g in module.genes():
            self._species[g.id] = f"{g.id}.mRNA"
        for p in module.proteins():
            self._species.setdefault(p.id, p.id)

    @property
    def nodes(self) -> list[str]:
        return list(self._species)

    @property
    def edges(self) -> list[Edge]:
        return [(r.source, r.target) for r in self.module.active_rules({})]

    def _perturbed(self, nodes_off: list[str], edges_off: list[Edge]) -> Module:
        m = copy.deepcopy(self.module)
        genes = {g.id: g for g in m.genes()}
        # a node held at zero is held through the integrator's own sub-steps, not clamped between them:
        # `NetworkRuntime.run`'s clamp is re-applied only after a whole step, so a clamped mRNA still
        # takes its basal rate inside the Runge-Kutta stages and its protein reads that. Measured on the
        # three-gene ring: the protein of a clamped gene settled at 1.87, a tenth of its 19.2, not zero.
        # Zeroing what makes the species instead holds it at exactly zero at every stage.
        # (The runtime's clamp holds every stage since 2026-09-28, grn_clamp_census.json; zeroing stays,
        # because it also cuts the protein's `produces` rules, which a clamp on the mRNA alone would not.)
        off = set(nodes_off)
        for name in nodes_off:
            if name in genes:
                genes[name].basal_rate = 0.0
                genes[name].attrs["max_rate"] = 0.0
        wanted = set(edges_off)
        kept = []
        for r in m.rules:
            if r.action is Action.PRODUCE and (r.target in off or (r.source, r.target) in wanted):
                continue  # the protein is not made: from its own knockout, or from this edge's
            if (r.source, r.target) not in wanted:
                kept.append(r)
            else:
                # the target keeps its rule and reads a source that is not in the state, which the
                # runtime reads as 0: dropping the rule instead would hand a gene its full max_rate
                r.source = f"{r.source}@zero"
                kept.append(r)
        m.rules = kept
        return m

    def measure(self, nodes_off: list[str], edges_off: list[Edge]) -> tuple[dict[str, float], str]:
        """Per node, the mean level over the second half of the run; the kind is steady or oscillating."""
        m = self._perturbed(nodes_off, edges_off) if (nodes_off or edges_off) else self.module
        off = {self._species[n] for n in nodes_off if n in self._species}
        traj = NetworkRuntime(m, seed=self.seed).run(
            hours=self.hours,
            dt=self.dt,
            initial={k: v for k, v in self.start.items() if k not in off},
            record_every=self.record_every,
        )
        half = len(traj.times) // 2
        levels, oscillating = {}, False
        for node, sp in self._species.items():
            xs = traj.levels[sp][half:]
            mean = sum(xs) / len(xs)
            levels[node] = mean
            peaks = sum(1 for i in range(1, len(xs) - 1) if xs[i - 1] < xs[i] > xs[i + 1])
            # a swing of at least 1% of the mean, or a flat line's floating-point wobble counts as a cycle
            oscillating = oscillating or (peaks >= 2 and max(xs) - min(xs) > 0.01 * mean)
        return levels, "oscillating" if oscillating else "steady"

    def change(self, wt: float, ko: float) -> float:
        return abs(math.log2((ko + _LEVEL_FLOOR) / (wt + _LEVEL_FLOOR)))


@dataclass(slots=True)
class NetworkKnockout:
    """One knockout of a network, its consequence, and where that consequence sits against chance."""

    nodes: list[str]
    edges: list[Edge]
    consequence: float
    per_node: dict[str, float]
    baseline: dict[str, float]  # per node, unperturbed: the fraction of the attractor on, or the level
    levels: dict[str, float]  # the same reading of the perturbed run
    unperturbed: str
    perturbed: str
    changed_kind: bool
    p_value: float
    draws: int
    exhaustive: bool  # every matched knockout was run, not a sample
    control: list[float] = field(default_factory=list)

    @property
    def label(self) -> str:
        parts = [f"-{n}" for n in self.nodes] + [f"-{a}{EDGE_SEPARATOR}{b}" for a, b in self.edges]
        return ", ".join(parts) or "unperturbed"

    @property
    def measured(self) -> list[str]:
        return sorted(self.per_node)

    def to_dict(self) -> dict:
        return {
            "knockout": self.label,
            "nodes": list(self.nodes),
            "edges": [f"{a}{EDGE_SEPARATOR}{b}" for a, b in self.edges],
            "consequence": self.consequence,
            "per_node": dict(self.per_node),
            "baseline": dict(self.baseline),
            "levels": dict(self.levels),
            "measured_nodes": self.measured,
            "unperturbed": self.unperturbed,
            "perturbed": self.perturbed,
            "changed_kind": self.changed_kind,
            "p": self.p_value,
            "control_draws": self.draws,
            "control_exhaustive": self.exhaustive,
            "control_max": max(self.control) if self.control else None,
        }

    def format(self) -> str:
        top = sorted(self.per_node.items(), key=lambda kv: -kv[1])[:5]
        moved = ", ".join(f"{n} {v:+.3g}" for n, v in top if v > 0) or "nothing moved"
        kind = f"{self.unperturbed} -> {self.perturbed}" + (" (kind changed)" if self.changed_kind else "")
        how = "all matched knockouts" if self.exhaustive else f"{self.draws} matched draws"
        return "\n".join(
            [
                f"knockout {self.label}: consequence {self.consequence:.4g} over "
                f"{len(self.per_node)} node(s) read in both runs",
                f"  dynamics  {kind}",
                f"  changed   {moved}",
                f"  chance    p = {self.p_value:.4g} against {how}",
            ]
        )


def _matched_sets(
    nodes: list[str], edges: list[Edge], n_nodes: int, n_edges: int, draws: int, seed: int
) -> tuple[list[tuple[list[str], list[Edge]]], bool]:
    """Knockouts of the same size, drawn only from nodes and edges the observed knockout left alone."""
    if n_nodes > len(nodes) or n_edges > len(edges):
        return [], False
    total = math.comb(len(nodes), n_nodes) * math.comb(len(edges), n_edges)
    if total <= draws:  # few enough matched knockouts exist: run all of them instead of sampling
        out = [
            (list(ns), list(es))
            for ns in itertools.combinations(nodes, n_nodes)
            for es in itertools.combinations(edges, n_edges)
        ]
        return out, True
    rng = random.Random(seed)
    return [
        (rng.sample(nodes, n_nodes), [tuple(e) for e in rng.sample(edges, n_edges)]) for _ in range(draws)
    ], False


def run(
    network: BooleanNetwork | Module,
    knockouts: list[str] | tuple[str, ...],
    start: dict | None = None,
    draws: int = CONTROL_DRAWS,
    seed: int = CONTROL_SEED,
    **engine: object,
) -> NetworkKnockout:
    """Hold `knockouts` (node names and `A>B` edges) at zero and read the consequence against chance.

    `network` is a Boolean model or a rule module for the GRN runtime; `start` is the initial state both
    runs share. The consequence is the mean change over the nodes knocked out in NEITHER run, and the
    p-value is its rank among matched random knockouts: the same number of nodes and edges, drawn only
    from the nodes and edges this knockout left alone (`NETWORK_KNOCKOUT_SPEC`, `..._FALSIFIER`).
    """
    nodes_off, edges_off = parse_knockouts(knockouts)
    eng: _BooleanEngine | _ContinuousEngine = (
        _BooleanEngine(network, start)  # type: ignore[arg-type]
        if isinstance(network, BooleanNetwork)
        else _ContinuousEngine(network, start, **engine)  # type: ignore[arg-type]
    )
    unknown = [n for n in nodes_off if n not in eng.nodes]
    bad_edges = [e for e in edges_off if e not in eng.edges]
    if unknown or bad_edges:
        named = unknown + [f"{a}{EDGE_SEPARATOR}{b}" for a, b in bad_edges]
        raise ValueError(f"the knockout names nothing in the network: {', '.join(named)}")

    wt, wt_kind = eng.measure([], [])  # the one unperturbed run both the knockout and the controls use

    def consequence(
        nodes: list[str], edges: list[Edge]
    ) -> tuple[float, dict[str, float], str, dict[str, float]]:
        ko, ko_kind = eng.measure(nodes, edges)
        read = [n for n in eng.nodes if n not in set(nodes)]  # knocked out in neither run
        per = {n: eng.change(wt[n], ko[n]) for n in read}
        return (sum(per.values()) / len(per) if per else 0.0), per, ko_kind, ko

    observed, per_node, ko_kind, ko_levels = consequence(nodes_off, edges_off)
    free_nodes = [n for n in eng.nodes if n not in set(nodes_off)]
    free_edges = [e for e in eng.edges if e not in set(edges_off)]
    matched, exhaustive = _matched_sets(free_nodes, free_edges, len(nodes_off), len(edges_off), draws, seed)
    control = [consequence(ns, es)[0] for ns, es in matched]
    ge = sum(1 for c in control if c >= observed - 1e-12)
    return NetworkKnockout(
        nodes=nodes_off,
        edges=edges_off,
        consequence=observed,
        per_node=per_node,
        baseline=wt,
        levels=ko_levels,
        unperturbed=wt_kind,
        perturbed=ko_kind,
        changed_kind=wt_kind.split()[0] != ko_kind.split()[0],
        p_value=(1 + ge) / (1 + len(control)),
        draws=len(control),
        exhaustive=exhaustive,
        control=control,
    )


def run_experiment(network: BooleanNetwork | Module, experiment: object, **kwargs: object) -> NetworkKnockout:
    """The same run, taking the knockout list off a parsed `experiment` block's `knockout:` line."""
    return run(network, list(getattr(experiment, "knockouts", [])), **kwargs)  # type: ignore[arg-type]
