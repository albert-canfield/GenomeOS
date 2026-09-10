"""Boolean (logical) network engine (task 2.4).

Reads BoolNet-style `.bnet` rules ("target, factors") and simulates
synchronous or asynchronous updates, finds attractors, and exports the
network as BioIR rules so logical models sit in the same representation as
everything else. This is the in-house stand-in for MaBoSS; a MaBoSS adapter
can replace the simulator without changing the file format or the IR.
"""

from __future__ import annotations

import ast
import random
from dataclasses import dataclass, field
from pathlib import Path

from genomeos.ir import Action, Entity, Evidence, EvidenceKind, Module, Rule

_ALLOWED = (
    ast.BoolOp,
    ast.UnaryOp,
    ast.Name,
    ast.Constant,
    ast.And,
    ast.Or,
    ast.Not,
    ast.Expression,
    ast.Load,
)


def _py(expression: str) -> str:
    """BoolNet operators to Python; stripped because eval-mode parsing rejects leading whitespace."""
    return expression.replace("&", " and ").replace("|", " or ").replace("!", " not ").strip()


def _edge_signs(expression: str) -> dict[str, set[bool]]:
    """For each input name, the set of negation states it appears under."""
    signs: dict[str, set[bool]] = {}

    def walk(node: ast.AST, negated: bool) -> None:
        if isinstance(node, ast.Name):
            signs.setdefault(node.id, set()).add(negated)
        elif isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.Not):
            walk(node.operand, not negated)
        elif isinstance(node, ast.BoolOp):
            for v in node.values:
                walk(v, negated)
        elif isinstance(node, ast.Expression):
            walk(node.body, negated)

    walk(ast.parse(_py(expression), mode="eval"), False)
    return signs


@dataclass(slots=True)
class BooleanNetwork:
    rules: dict[str, str] = field(default_factory=dict)  # node -> BoolNet expression
    _code: dict[str, object] = field(default_factory=dict, repr=False)

    @classmethod
    def from_bnet(cls, text: str) -> BooleanNetwork:
        net = cls()
        for raw in text.splitlines():
            line = raw.split("#", 1)[0].strip()
            if not line or line.lower().startswith("targets"):
                continue
            target, _, expr = line.partition(",")
            net.add(target.strip(), expr.strip())
        return net

    @classmethod
    def from_file(cls, path: str | Path) -> BooleanNetwork:
        return cls.from_bnet(Path(path).read_text())

    def add(self, node: str, expression: str) -> None:
        tree = ast.parse(_py(expression), mode="eval")
        for n in ast.walk(tree):
            if not isinstance(n, _ALLOWED):
                raise ValueError(f"unsupported syntax in rule for {node}: {expression!r}")
        self.rules[node] = expression
        self._code[node] = compile(tree, f"<rule {node}>", "eval")

    @property
    def nodes(self) -> list[str]:
        return list(self.rules)

    def inputs_of(self, node: str) -> set[str]:
        return set(_edge_signs(self.rules[node]))

    def _eval(self, node: str, env: dict[str, bool]) -> bool:
        return bool(eval(self._code[node], {"__builtins__": {}}, env))  # noqa: S307  (AST-validated)

    def step(self, state: dict[str, bool]) -> dict[str, bool]:
        env = {k: bool(v) for k, v in state.items()}
        return {n: self._eval(n, env) for n in self.rules}

    def step_async(self, state: dict[str, bool], rng: random.Random) -> dict[str, bool]:
        node = rng.choice(self.nodes)
        new = dict(state)
        new[node] = self._eval(node, {k: bool(v) for k, v in state.items()})
        return new

    def trajectory(self, initial: dict[str, bool], steps: int = 50) -> list[dict[str, bool]]:
        state = {n: bool(initial.get(n, False)) for n in self.rules}
        out = [state]
        for _ in range(steps):
            state = self.step(state)
            out.append(state)
        return out

    def attractor(self, initial: dict[str, bool], max_steps: int = 10_000) -> list[dict[str, bool]]:
        """Synchronous attractor reached from `initial`: a list of states (length 1 = fixed point)."""
        state = {n: bool(initial.get(n, False)) for n in self.rules}
        seen: dict[tuple[bool, ...], int] = {}
        history: list[dict[str, bool]] = []
        for i in range(max_steps):
            key = tuple(state[n] for n in self.rules)
            if key in seen:
                return history[seen[key] :]
            seen[key] = i
            history.append(state)
            state = self.step(state)
        raise RuntimeError("no attractor found within max_steps")

    def attractors(self, samples: int = 200, seed: int = 0) -> list[list[dict[str, bool]]]:
        """Distinct synchronous attractors found from random initial states."""
        rng = random.Random(seed)
        found: dict[frozenset[tuple[bool, ...]], list[dict[str, bool]]] = {}
        for _ in range(samples):
            init = {n: rng.random() < 0.5 for n in self.rules}
            att = self.attractor(init)
            key = frozenset(tuple(s[n] for n in self.rules) for s in att)
            found.setdefault(key, att)
        return sorted(found.values(), key=len)

    def to_module(self, name: str, evidence: Evidence | None = None, confidence: float = 0.5) -> Module:
        """Each node becomes an entity; each input edge a rule signed by how it appears in the expression."""
        ev = evidence or Evidence(EvidenceKind.CURATED, "Boolean model")
        m = Module(name=name)
        for n in self.rules:
            m.add(
                Entity(id=n, kind="node", attrs={"rule": self.rules[n]}, evidence=ev, confidence=confidence)
            )
        for target, expr in self.rules.items():
            for src, negs in _edge_signs(expr).items():
                action = (
                    Action.INHIBIT
                    if negs == {True}
                    else Action.ACTIVATE
                    if negs == {False}
                    else Action.MODIFY
                )
                m.rules.append(
                    Rule(
                        id=f"{src} {action.value} {target}",
                        source=src,
                        action=action,
                        target=target,
                        evidence=ev,
                        confidence=confidence,
                    )
                )
        return m
