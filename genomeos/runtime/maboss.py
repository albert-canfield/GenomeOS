# SPDX-License-Identifier: Apache-2.0
# Part of the BioLang engine (language, IR, VM, standard library); see LICENSING.md.
"""MaBoSS adapter for Boolean networks: continuous-time stochastic simulation behind an extra.

`runtime/boolean.py` updates a BoolNet-style network in discrete steps, synchronously or one random
node at a time, and finds attractors. MaBoSS (Stoll et al. 2012, BMC Syst Biol 6:116; BSD-3-Clause,
an optional extra: `uv sync --extra maboss`, which installs its Python engine `cmaboss`) treats the
same logic as a continuous-time Markov process: each node flips at a rate while its rule disagrees
with its value, and the output is the probability of each network state over time, estimated over
many trajectories.

Two ways in, both refused where the translation would not be faithful:

  * `MaBoSSModel.from_network(net, ...)` writes a `BooleanNetwork` (the subset `boolean.py` runs:
    names, `&`, `|`, `!`, parentheses, 0/1) as a MaBoSS `.bnd`/`.cfg` pair. Every node gets MaBoSS's
    default rates (up 1 while the rule is on, down 1 while it is off). Anything else in a rule raises
    `TranslationError`.
  * `MaBoSSModel.from_files(bnd, cfg)` takes a published pair as it was published.

`MaBoSSModel.held({...})` holds nodes at a value from time 0 (a knockout or a forced gain): the
node's rates are both set to 0 and its initial state to the value, so it never moves.

`to_network(model)` goes the other way, for the stand-in: it returns a `BooleanNetwork` only when the
two semantics coincide -- every node's rates in MaBoSS's default form and one common positive value,
held nodes aside. Then MaBoSS's process and the stand-in's `step_async` under an exponential clock of
rate (number of nodes) x (that value) are the same Markov chain. Unequal rates, rates that read other
nodes, joint or fractional initial states and logic outside the subset raise `TranslationError`.

It is chosen by name (`boolean.runtime(net, engine="maboss")`) and never replaces the stand-in
silently: without the extra, building a `MaBoSSRuntime` raises `EngineUnavailableError` naming it.
"""

from __future__ import annotations

import ast
import re
from dataclasses import dataclass, field
from pathlib import Path

from .boolean import BooleanNetwork, _py
from .roadrunner import EngineUnavailableError

ENGINE = "maboss"
EXTRA = "maboss"
HOW_TO_ENABLE = f"install MaBoSS's Python engine (cmaboss) with `uv sync --extra {EXTRA}`"

#: MaBoSS's name for a state with no non-internal node on
NIL = "<nil>"

_NAME = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
#: words MaBoSS reads as operators or keywords, so a node may not be called by them
_RESERVED = frozenset({"node", "logic", "rate_up", "rate_down", "and", "or", "not", "xor", "true", "false"})

__all__ = [
    "ENGINE",
    "EXTRA",
    "NIL",
    "EngineUnavailableError",
    "MaBoSSModel",
    "MaBoSSResult",
    "MaBoSSRuntime",
    "TranslationError",
    "available",
    "to_network",
    "version",
]


class TranslationError(ValueError):
    """A construct the adapter does not translate faithfully; it refuses rather than guess."""


def _import():
    try:
        import cmaboss  # type: ignore[import-not-found]
    except ImportError as e:  # the extra is absent: say which one, never fall back
        raise EngineUnavailableError(f"MaBoSS is not installed; {HOW_TO_ENABLE}") from e
    return cmaboss


def available() -> bool:
    """True when MaBoSS's Python engine can be imported in this environment."""
    try:
        _import()
    except EngineUnavailableError:
        return False
    return True


def version() -> str:
    from importlib.metadata import version as dist_version

    _import()
    return dist_version("cmaboss")


# --------------------------------------------------------------------------------------------------
# BioLang (the BoolNet subset boolean.py runs) -> MaBoSS


def _check_name(name: str) -> str:
    if not _NAME.match(name) or name.lower() in _RESERVED:
        raise TranslationError(f"node name {name!r} is not one MaBoSS can read")
    return name


def _render(node: ast.AST) -> str:
    """One validated rule as MaBoSS logic; anything outside names, and, or, not and 0/1 refuses."""
    if isinstance(node, ast.Expression):
        return _render(node.body)
    if isinstance(node, ast.Name):
        return _check_name(node.id)
    if isinstance(node, ast.Constant) and type(node.value) in (bool, int) and node.value in (0, 1):
        return "1" if node.value else "0"
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.Not):
        return f"!{_render(node.operand)}"
    if isinstance(node, ast.BoolOp) and isinstance(node.op, (ast.And, ast.Or)):
        joiner = " & " if isinstance(node.op, ast.And) else " | "
        return "(" + joiner.join(_render(v) for v in node.values) + ")"
    raise TranslationError(f"not in the translated subset: {ast.dump(node)[:80]}")


def logic_of(expression: str) -> str:
    """A BoolNet rule as `boolean.py` reads it, written as MaBoSS logic."""
    try:
        tree = ast.parse(_py(expression), mode="eval")
    except SyntaxError as e:
        raise TranslationError(f"rule does not parse: {expression!r}") from e
    return _render(tree)


# --------------------------------------------------------------------------------------------------
# the .bnd / .cfg subset this module reads back


def _strip_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
    return re.sub(r"//[^\n]*", " ", text)


_NODE_BLOCK = re.compile(r"\b[Nn]ode\s+([A-Za-z_][A-Za-z0-9_]*)\s*\{(.*?)\}", re.S)


def parse_bnd(text: str) -> dict[str, dict[str, str]]:
    """Node -> {field: expression} for every `Node X { field = expr; ... }` block."""
    nodes: dict[str, dict[str, str]] = {}
    for name, body in _NODE_BLOCK.findall(_strip_comments(text)):
        fields: dict[str, str] = {}
        for stmt in body.split(";"):
            if not stmt.strip():
                continue
            key, eq, value = stmt.partition("=")
            if not eq:
                raise TranslationError(f"node {name}: statement without '=': {stmt.strip()!r}")
            fields[key.strip()] = value.strip()
        nodes[name] = fields
    return nodes


_ONE_NODE_ISTATE = re.compile(r"^\[\s*([A-Za-z_][A-Za-z0-9_]*)\s*\]\s*\.\s*istate$")
_DISTRIBUTION = re.compile(
    r"^\s*([0-9.eE+-]+)\s*\[\s*([01])\s*\]\s*(?:,\s*([0-9.eE+-]+)\s*\[\s*([01])\s*\]\s*)?$"
)


def _one_node_istate(value: str) -> str:
    """`p [0], q [1]` for one node as 0, 1, -1 (random, 1/2 each) or, for anything else, the text."""
    m = _DISTRIBUTION.match(value)
    if not m:
        return value
    weights = {m.group(2): float(m.group(1))}
    if m.group(3) is not None:
        weights[m.group(4)] = weights.get(m.group(4), 0.0) + float(m.group(3))
    on, off = weights.get("1", 0.0), weights.get("0", 0.0)
    if on > 0 and off == 0:
        return "1"
    if off > 0 and on == 0:
        return "0"
    return "-1" if on == off else value


def parse_cfg(text: str) -> dict[str, object]:
    """Parameters, per-node istate and is_internal, joint istates, and the run settings. A one-node
    bracketed istate (`[A].istate = 0.5 [0], 0.5 [1];`) is read as that node's istate."""
    out: dict[str, object] = {"params": {}, "istate": {}, "internal": {}, "joint_istate": [], "settings": {}}
    for stmt in _strip_comments(text).split(";"):
        stmt = stmt.strip()
        if not stmt:
            continue
        key, eq, value = stmt.partition("=")
        key, value = key.strip(), value.strip()
        if not eq:
            raise TranslationError(f"cfg statement without '=': {stmt!r}")
        if key.startswith("$"):
            out["params"][key] = value  # type: ignore[index]
        elif _ONE_NODE_ISTATE.match(key):
            out["istate"][_ONE_NODE_ISTATE.match(key).group(1)] = _one_node_istate(value)  # type: ignore[index,union-attr]
        elif key.startswith("["):
            out["joint_istate"].append(stmt)  # type: ignore[union-attr]
        elif key.endswith(".istate"):
            out["istate"][key[: -len(".istate")]] = value  # type: ignore[index]
        elif key.endswith(".is_internal"):
            out["internal"][key[: -len(".is_internal")]] = value  # type: ignore[index]
        else:
            out["settings"][key] = value  # type: ignore[index]
    return out


# A per-node istate written as a distribution without brackets (`A.istate = 0.5 [0], 0.5 [1];`) is
# not MaBoSS syntax, and cmaboss 1.0.0b32 does not reject it: the interpreter dies with SIGSEGV
# (measured 2026-10-10). The bracketed form `[A].istate = ...` is the valid one.
_BAD_ISTATE = re.compile(r"(^|;)\s*[A-Za-z_][A-Za-z0-9_]*\s*\.\s*istate\s*=[^;]*\[", re.M)


@dataclass(slots=True)
class MaBoSSModel:
    """A MaBoSS `.bnd`/`.cfg` pair, as text, and where it came from."""

    bnd: str
    cfg: str
    source: str = ""

    def __post_init__(self) -> None:
        if _BAD_ISTATE.search(_strip_comments(self.cfg)):
            raise TranslationError(
                "a per-node istate written as a distribution; MaBoSS reads that form only as "
                "`[NODE].istate = p [0], q [1];` (and cmaboss crashes on the other form)"
            )

    @classmethod
    def from_files(cls, bnd: str | Path, cfg: str | Path) -> MaBoSSModel:
        return cls(Path(bnd).read_text(), Path(cfg).read_text(), source=f"{bnd} + {cfg}")

    @classmethod
    def from_network(
        cls,
        net: BooleanNetwork,
        initial: dict[str, bool] | None = None,
        random_nodes: list[str] | tuple[str, ...] = (),
        outputs: list[str] | tuple[str, ...] | None = None,
        max_time: float = 50.0,
        time_tick: float = 0.5,
        sample_count: int = 10_000,
        seed: int = 0,
        thread_count: int = 1,
    ) -> MaBoSSModel:
        """`net` with MaBoSS's default rates. A node starts at `initial`'s value, or random (1/2) when
        listed in `random_nodes`, or off -- the stand-in's own default. `outputs` are the nodes the
        state probabilities are reported over (all of them when None)."""
        initial = dict(initial or {})
        unknown = (set(initial) | set(random_nodes) | set(outputs or ())) - set(net.rules)
        if unknown:
            raise TranslationError(f"not nodes of the network: {sorted(unknown)}")
        both = set(initial) & set(random_nodes)
        if both:
            raise TranslationError(f"both fixed and random at the start: {sorted(both)}")
        blocks = []
        for node, expression in net.rules.items():
            blocks.append(
                f"Node {_check_name(node)} {{\n  logic = {logic_of(expression)};\n"
                "  rate_up = @logic ? 1 : 0;\n  rate_down = @logic ? 0 : 1;\n}\n"
            )
        shown = set(net.rules) if outputs is None else set(outputs)
        cfg = [
            f"{n}.istate = {int(bool(initial.get(n, False)))};" for n in net.rules if n not in random_nodes
        ]
        cfg += [f"[{n}].istate = 0.5 [0], 0.5 [1];" for n in net.rules if n in random_nodes]
        cfg += [f"{n}.is_internal = {0 if n in shown else 1};" for n in net.rules]
        cfg += [
            "discrete_time = 0;",
            "use_physrandgen = 0;",
            f"seed_pseudorandom = {int(seed)};",
            f"sample_count = {int(sample_count)};",
            f"max_time = {float(max_time)};",
            f"time_tick = {float(time_tick)};",
            f"thread_count = {int(thread_count)};",
        ]
        return cls("\n".join(blocks), "\n".join(cfg) + "\n", source="BooleanNetwork")

    @property
    def nodes(self) -> list[str]:
        return list(parse_bnd(self.bnd))

    def held(self, values: dict[str, bool]) -> MaBoSSModel:
        """This model with each named node held at its value from time 0 (rates 0, istate the value)."""
        blocks = parse_bnd(self.bnd)
        missing = set(values) - set(blocks)
        if missing:
            raise TranslationError(f"not nodes of the model: {sorted(missing)}")
        joint = " ".join(parse_cfg(self.cfg)["joint_istate"])  # type: ignore[arg-type]
        for node in values:
            if re.search(rf"[\[,\s]{re.escape(node)}[\],\s]", joint):
                raise TranslationError(f"{node} has a joint initial state; holding it would change others")

        def rewrite(match: re.Match[str]) -> str:
            name, body = match.group(1), match.group(2)
            if name not in values:
                return match.group(0)
            kept = [
                s
                for s in body.split(";")
                if s.strip() and s.split("=")[0].strip() not in ("rate_up", "rate_down")
            ]
            stmts = "".join(f"{s.rstrip()};" for s in kept)
            return f"Node {name} {{{stmts}\n  rate_up = 0;\n  rate_down = 0;\n}}"

        bnd = _NODE_BLOCK.sub(rewrite, _strip_comments(self.bnd))
        cfg = _strip_comments(self.cfg)
        for node, value in values.items():
            own = rf"(^|;)\s*(?:{re.escape(node)}|\[\s*{re.escape(node)}\s*\])\s*\.\s*istate\s*=[^;]*;"
            cfg = re.sub(own, r"\1", cfg, flags=re.M)
            cfg += f"\n{node}.istate = {int(bool(value))};"
        label = ", ".join(f"{n}={int(bool(v))}" for n, v in values.items())
        return MaBoSSModel(bnd, cfg + "\n", source=f"{self.source} held {label}")


# --------------------------------------------------------------------------------------------------
# running it


@dataclass(slots=True)
class MaBoSSResult:
    """The last time window's state probabilities over the non-internal nodes, as MaBoSS reports them.

    `states` maps each state, written as the sorted tuple of its active non-internal nodes (`()` is
    MaBoSS's `<nil>`), to its probability over the window [`window`[0], `window`[1]); `errors` is
    MaBoSS's own error estimate for each; `nodes` is each non-internal node's probability of being on.
    """

    states: dict[tuple[str, ...], float]
    errors: dict[tuple[str, ...], float]
    nodes: dict[str, float]
    window: tuple[float, float]
    settings: dict[str, object] = field(default_factory=dict)


def _state(label: str) -> tuple[str, ...]:
    return () if label.strip() == NIL else tuple(sorted(part.strip() for part in label.split("--")))


class MaBoSSRuntime:
    """A `MaBoSSModel` run on MaBoSS (cmaboss). Each run builds a fresh simulation."""

    engine = ENGINE

    def __init__(self, model: MaBoSSModel) -> None:
        self._cm = _import()
        self.model = model

    @classmethod
    def from_network(cls, net: BooleanNetwork, **settings: object) -> MaBoSSRuntime:
        return cls(MaBoSSModel.from_network(net, **settings))  # type: ignore[arg-type]

    def _sim(self):
        return self._cm.MaBoSSSim(network_str=self.model.bnd, config_str=self.model.cfg)

    def settings(self) -> dict[str, object]:
        """The run settings as MaBoSS itself reads them back from the cfg, for a result to record."""
        read = parse_cfg(self._sim().str_cfg())["settings"]
        keys = (
            "discrete_time",
            "use_physrandgen",
            "seed_pseudorandom",
            "sample_count",
            "max_time",
            "time_tick",
            "thread_count",
        )
        return {"engine": ENGINE, "cmaboss": version(), **{k: read.get(k) for k in keys}}  # type: ignore[union-attr]

    def run(self) -> MaBoSSResult:
        res = self._sim().run()
        probs, times, labels, errors = res.get_last_probtraj()
        node_probs, _, node_labels, _ = res.get_last_nodes_probtraj()
        tick = float(parse_cfg(self.model.cfg)["settings"].get("time_tick", 1))  # type: ignore[union-attr]
        states: dict[tuple[str, ...], float] = {}
        errs: dict[tuple[str, ...], float] = {}
        for i, label in enumerate(labels):
            key = _state(label)
            states[key] = states.get(key, 0.0) + float(probs[0][i])
            errs[key] = float(errors[0][i])
        nodes = {str(n): float(node_probs[0][i]) for i, n in enumerate(node_labels)}
        t0 = float(times[0])
        return MaBoSSResult(states, errs, nodes, (t0, t0 + tick), self.settings())


# --------------------------------------------------------------------------------------------------
# MaBoSS -> the stand-in, only where the two are the same process

_RATE_UP = re.compile(r"^\(?\s*@logic\s*\?\s*([^:]+?)\s*:\s*0(?:\.0*)?\s*\)?$")
_RATE_DOWN = re.compile(r"^\(?\s*@logic\s*\?\s*0(?:\.0*)?\s*:\s*([^)]+?)\s*\)?$")
_LOGIC_TOKENS = re.compile(r"\s*(?:(AND|OR|NOT)\b|([A-Za-z_][A-Za-z0-9_]*)|(&&|\|\||[&|!()])|([01])\b)")


def _number(text: str, params: dict[str, str], where: str) -> float:
    text = text.strip()
    if text in params:
        text = params[text].strip()
    try:
        return float(text)
    except ValueError:
        raise TranslationError(f"{where}: rate {text!r} is not a number or a parameter set to one") from None


def _bnet_logic(expression: str, node: str) -> str:
    """MaBoSS logic in the subset (names, AND/OR/NOT or &,|,!, parentheses, 0/1) as a BoolNet rule."""
    out, pos = [], 0
    text = expression.strip()
    while pos < len(text):
        m = _LOGIC_TOKENS.match(text, pos)
        if not m or m.end() == pos:
            raise TranslationError(f"node {node}: logic outside the subset: {expression!r}")
        word, name, op, const = m.groups()
        if word:
            out.append({"AND": " & ", "OR": " | ", "NOT": "!"}[word])
        elif name:
            out.append(name)
        elif op:
            out.append({"&&": " & ", "||": " | "}.get(op, op))
        else:
            out.append(const)
        pos = m.end()
        while pos < len(text) and text[pos].isspace():
            pos += 1
    return "".join(out)


def to_network(model: MaBoSSModel) -> tuple[BooleanNetwork, dict[str, bool], list[str], float]:
    """The stand-in's network for a model whose MaBoSS process is the stand-in's asynchronous update.

    Returns (network, fixed initial values, nodes that start random at 1/2, the common rate). Held
    nodes (both rates 0, or the one rate that could move them 0) become constant rules at their
    initial value. Everything else must be MaBoSS's default form with one common positive rate.
    """
    nodes = parse_bnd(model.bnd)
    cfg = parse_cfg(model.cfg)
    if cfg["joint_istate"]:
        raise TranslationError(f"joint initial states are not translated: {cfg['joint_istate']}")
    params: dict[str, str] = cfg["params"]  # type: ignore[assignment]
    istate: dict[str, str] = cfg["istate"]  # type: ignore[assignment]
    initial: dict[str, bool] = {}
    random_nodes: list[str] = []
    for node in nodes:
        value = istate.get(node, "-1").strip()
        if value in ("0", "1"):
            initial[node] = value == "1"
        elif value == "-1":
            random_nodes.append(node)
        else:
            raise TranslationError(f"{node}: initial state {value!r} is neither 0, 1 nor random")
    rates: set[float] = set()
    net = BooleanNetwork()
    for node, fields in nodes.items():
        extra = set(fields) - {"logic", "rate_up", "rate_down"}
        if extra:
            raise TranslationError(f"node {node}: fields not translated: {sorted(extra)}")
        if "logic" not in fields:
            raise TranslationError(f"node {node} has no logic")
        rule = _bnet_logic(fields["logic"], node)
        up = _RATE_UP.match(fields.get("rate_up", "@logic ? 1 : 0").strip())
        down = _RATE_DOWN.match(fields.get("rate_down", "@logic ? 0 : 1").strip())
        up_rate = 0.0 if fields.get("rate_up", "").strip() in ("0", "0.0") else None
        down_rate = 0.0 if fields.get("rate_down", "").strip() in ("0", "0.0") else None
        if up_rate is None:
            if not up:
                raise TranslationError(
                    f"node {node}: rate_up {fields.get('rate_up')!r} is not @logic ? r : 0"
                )
            up_rate = _number(up.group(1), params, node)
        if down_rate is None:
            if not down:
                raise TranslationError(
                    f"node {node}: rate_down {fields.get('rate_down')!r} is not @logic ? 0 : r"
                )
            down_rate = _number(down.group(1), params, node)
        if up_rate == 0 and down_rate == 0:
            held = node in initial
            value = initial.get(node, False)
        elif up_rate == 0:
            held, value = initial.get(node) is False, False
        elif down_rate == 0:
            held, value = initial.get(node) is True, True
        else:
            held, value = False, False
        if held:
            net.add(node, "1" if value else "0")
            continue
        if up_rate == 0 or down_rate == 0:
            raise TranslationError(f"node {node}: one rate 0 and a start that lets it move one way only")
        rates |= {up_rate, down_rate}
        net.add(node, rule)
    if len(rates) > 1:
        raise TranslationError(
            f"unequal rates {sorted(rates)}: the stand-in picks nodes uniformly, which is MaBoSS's process "
            "only when every rate is the same"
        )
    if rates and min(rates) <= 0:
        raise TranslationError(f"a negative rate: {sorted(rates)}")
    return net, initial, random_nodes, (rates.pop() if rates else 1.0)
