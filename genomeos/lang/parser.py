"""BioLang v0.2 parser: source text -> BioIR Module.

Grammar (see docs/BIOLANG-v0.2.md):

    module <dotted.name>
    import <dotted.name> | <relative/path.bio>

    gene <Id> { <props>  transcript <Id> { exons: ...; cds: ... } ... }
    protein <Id> { <props> }
    region <Id> { locus: chr7:1000-2000; role: unknown }
    cell_type <Id> { expresses: a, b; parent: X; ontology: CL:0000540 }
    event <Id> { rate: 0.5 /yr; when: cell_type = X; effect: telomere_bp -= 70 bp; ... }
    rule <Source> activates|inhibits|produces|binds|degrades|modifies <Target> { <props> }
    param <name> = <number> [unit] { evidence: ...; confidence: ... }

Properties are `key: value`, one per line or `;`-separated; a block may sit on
one line. Blocks may nest (transcript inside gene). `#` starts a comment.
Common keys on any block: evidence, confidence.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from genomeos.coords import Locus
from genomeos.ir import (
    UNKNOWN,
    Action,
    CellType,
    Effect,
    Event,
    Evidence,
    EvidenceKind,
    Gene,
    Module,
    Parameter,
    Protein,
    Region,
    Rule,
    Transcript,
)

_ACTIONS = {a.value: a for a in Action}
_KINDS = ("gene", "protein", "region", "rule", "param", "cell_type", "event", "transcript")
_HEADER = re.compile(r"^(" + "|".join(_KINDS) + r")\s+([^{]*?)\s*\{(.*)$")
_PARAM = re.compile(r"^(\w[\w.]*)\s*=\s*([-+0-9.eE]+)\s*([^\s{]*)\s*$")
_EFFECT = re.compile(r"^(\w+)\s*(\+=|-=|\*=|=)\s*([-+0-9.eE]+)\s*(.*)$")
STD_DIR = Path(__file__).resolve().parent.parent / "std"


class BioLangError(SyntaxError):
    pass


@dataclass
class Block:
    kind: str
    header: str
    props: dict[str, str] = field(default_factory=dict)
    children: list[Block] = field(default_factory=list)
    line: int = 0


# ---------------------------------------------------------------- lexing ---


def _split_statements(text: str) -> list[str]:
    """Split on ';' but not inside double-quoted strings."""
    out: list[str] = []
    buf: list[str] = []
    quoted = False
    for ch in text:
        if ch == '"':
            quoted = not quoted
        if ch == ";" and not quoted:
            out.append("".join(buf))
            buf = []
        else:
            buf.append(ch)
    out.append("".join(buf))
    return out


def _strip(line: str) -> str:
    return line.split("#", 1)[0].strip()


def _parse_blocks(lines: list[str]) -> tuple[list[tuple[str, str]], list[Block]]:
    """Return (directives, top-level blocks). Directives are (module|import, value)."""
    directives: list[tuple[str, str]] = []
    top: list[Block] = []
    stack: list[Block] = []

    def add_props(block: Block, text: str, line_no: int) -> None:
        for stmt in _split_statements(text):
            stmt = stmt.strip()
            if not stmt:
                continue
            k, sep, v = stmt.partition(":")
            if not sep:
                raise BioLangError(f"line {line_no}: expected 'key: value', got {stmt!r}")
            key = k.strip()
            if key == "effect" and key in block.props:
                block.props[key] += " ; " + v.strip()  # several effects per event
            else:
                block.props[key] = v.strip()

    for i, raw in enumerate(lines, 1):
        line = _strip(raw)
        if not line:
            continue
        if not stack and line.startswith("module "):
            directives.append(("module", line[7:].strip()))
            continue
        if not stack and line.startswith("import "):
            directives.append(("import", line[7:].strip()))
            continue
        # a line may contain: header { ... [}] , or props, or }
        while line:
            m = _HEADER.match(line)
            if m:
                blk = Block(m.group(1), m.group(2).strip(), line=i)
                (stack[-1].children if stack else top).append(blk)
                stack.append(blk)
                line = m.group(3).strip()
                continue
            if "}" in line:
                before, _, after = line.partition("}")
                if not stack:
                    raise BioLangError(f"line {i}: unexpected '}}'")
                add_props(stack[-1], before, i)
                stack.pop()
                line = after.strip()
                continue
            if not stack:
                raise BioLangError(f"line {i}: cannot parse {line!r}")
            add_props(stack[-1], line, i)
            line = ""
    if stack:
        raise BioLangError(f"line {stack[-1].line}: unterminated block {stack[-1].header!r}")
    return directives, top


# ----------------------------------------------------------- value parsing ---


def _parse_evidence(value: str) -> Evidence:
    m = re.match(r'^(\w+)\s*(?:"([^"]*)")?\s*(.*)$', value)
    if not m:
        raise BioLangError(f"bad evidence: {value!r}")
    kind, source, note = m.group(1).lower(), m.group(2) or "", m.group(3).strip()
    try:
        return Evidence(EvidenceKind(kind), source, note=note)
    except ValueError:
        raise BioLangError(
            f"unknown evidence kind {kind!r}; use one of {[k.value for k in EvidenceKind]}"
        ) from None


def _parse_when(value: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for clause in value.split(","):
        if not clause.strip():
            continue
        k, _, v = clause.partition("=")
        if not v:
            raise BioLangError(f"bad when clause: {clause!r}")
        out[k.strip()] = v.strip()
    return out


def _float(value: str, key: str, line_no: int) -> float:
    try:
        return float(value.split()[0])
    except (ValueError, IndexError):
        raise BioLangError(f"line {line_no}: {key} expects a number, got {value!r}") from None


def _list(value: str) -> list[str]:
    return [x.strip() for x in value.split(",") if x.strip()]


def _common(b: Block) -> tuple[Evidence, float]:
    ev = _parse_evidence(b.props["evidence"]) if "evidence" in b.props else Evidence()
    conf = _float(b.props["confidence"], "confidence", b.line) if "confidence" in b.props else 0.0
    if not 0.0 <= conf <= 1.0:
        raise BioLangError(f"line {b.line}: confidence must be within 0..1")
    return ev, conf


# ------------------------------------------------------------- compiling ---


def _compile_block(b: Block, module: Module) -> None:
    ev, conf = _common(b)
    p = b.props
    if b.kind == "gene":
        g = Gene(id=b.header, kind="gene", symbol=p.get("symbol", b.header), evidence=ev, confidence=conf)
        if "locus" in p:
            g.locus = Locus.parse(p["locus"])
        if "basal" in p:
            g.basal_rate = _float(p["basal"], "basal", b.line)
        if "max" in p:
            g.attrs["max_rate"] = _float(p["max"], "max", b.line)
        module.add(g)
        if "produces" in p:
            for target in _list(p["produces"]):
                module.rules.append(
                    Rule(
                        id=f"{b.header}->{target}",
                        source=b.header,
                        action=Action.PRODUCE,
                        target=target,
                        evidence=ev,
                        confidence=conf,
                    )
                )
        for child in b.children:
            if child.kind != "transcript":
                raise BioLangError(f"line {child.line}: only transcript blocks may nest inside a gene")
            cev, cconf = _common(child)
            tx = Transcript(
                id=child.header, kind="transcript", gene_id=g.id, evidence=cev or ev, confidence=cconf or conf
            )
            if "exons" in child.props:
                tx.exons = [Locus.parse(x) for x in _list(child.props["exons"])]
            if "cds" in child.props:
                tx.cds_segments = [Locus.parse(x) for x in _list(child.props["cds"])]
                tx.cds = Locus(
                    tx.cds_segments[0].chrom,
                    min(c.start for c in tx.cds_segments),
                    max(c.end for c in tx.cds_segments),
                    tx.cds_segments[0].strand,
                )
            module.add(tx)
            g.transcripts.append(tx)
        return
    if b.children:
        raise BioLangError(f"line {b.line}: {b.kind} blocks cannot contain nested blocks")
    if b.kind == "protein":
        pr = Protein(id=b.header, kind="protein", evidence=ev, confidence=conf)
        if "half_life" in p:
            pr.half_life_h = _float(p["half_life"], "half_life", b.line)
        if "sequence" in p:
            pr.sequence = p["sequence"]
        module.add(pr)
    elif b.kind == "region":
        role = p.get("role", "unknown")
        r = Region(
            id=b.header,
            kind="region",
            evidence=ev,
            confidence=conf,
            role=UNKNOWN if role.lower() == "unknown" else role,
        )
        if "locus" in p:
            r.locus = Locus.parse(p["locus"])
        module.add(r)
    elif b.kind == "cell_type":
        module.add(
            CellType(
                id=b.header,
                kind="cell_type",
                name=p.get("name", b.header),
                parent=p.get("parent", ""),
                expresses=_list(p.get("expresses", "")),
                ontology_id=p.get("ontology", ""),
                evidence=ev,
                confidence=conf,
            )
        )
    elif b.kind == "event":
        e = Event(id=b.header, evidence=ev, confidence=conf)
        if "rate" in p:
            parts = p["rate"].split()
            e.rate = _float(parts[0], "rate", b.line)
            if len(parts) > 1:
                e.rate_unit = parts[1]
        if "when" in p:
            e.when = _parse_when(p["when"])
        for eff in p.get("effect", "").split(" ; "):
            eff = eff.strip()
            if not eff:
                continue
            m = _EFFECT.match(eff)
            if not m:
                raise BioLangError(f"line {b.line}: bad effect {eff!r} (expected 'var op number [unit]')")
            e.effects.append(Effect(m.group(1), m.group(2), float(m.group(3)), m.group(4).strip()))
        module.events.append(e)
    elif b.kind == "rule":
        parts = b.header.split()
        if len(parts) != 3 or parts[1] not in _ACTIONS:
            raise BioLangError(
                f"line {b.line}: rule header must be '<source> <action> <target>' "
                f"with action in {sorted(_ACTIONS)}; got {b.header!r}"
            )
        src, act, tgt = parts
        rule = Rule(
            id=p.get("id", f"{src} {act} {tgt}"),
            source=src,
            action=_ACTIONS[act],
            target=tgt,
            evidence=ev,
            confidence=conf,
        )
        for key in ("strength", "threshold", "hill"):
            if key in p:
                setattr(rule, key, _float(p[key], key, b.line))
        if "when" in p:
            rule.when = _parse_when(p["when"])
        module.rules.append(rule)
    elif b.kind == "param":
        pm = _PARAM.match(b.header)
        if not pm:
            raise BioLangError(f"line {b.line}: param must be 'name = value [unit]'")
        module.parameters[pm.group(1)] = Parameter(
            name=pm.group(1), value=float(pm.group(2)), unit=pm.group(3), evidence=ev, confidence=conf
        )
    elif b.kind == "transcript":
        raise BioLangError(f"line {b.line}: transcript blocks must be nested inside a gene")


def _check_references(module: Module) -> None:
    for r in module.rules:
        for ref in (r.source, r.target):
            if ref not in module.entities:
                raise BioLangError(f"rule {r.id!r} refers to undeclared entity {ref!r}")
    for c in module.cell_types():
        if c.parent and c.parent not in module.entities:
            raise BioLangError(f"cell_type {c.id!r} has undeclared parent {c.parent!r}")


def resolve_import(name: str, base_dir: Path | None) -> Path:
    """`bio.std.ageing` -> genomeos/std/ageing.bio; `a.b` -> <base>/a/b.bio; or a literal path."""
    if name.endswith(".bio"):
        candidates = [Path(name)] + ([base_dir / name] if base_dir else [])
    else:
        rel = Path(*name.split(".")[2:]) if name.startswith("bio.std.") else Path(*name.split("."))
        candidates = [STD_DIR / rel.with_suffix(".bio")] + (
            [base_dir / rel.with_suffix(".bio")] if base_dir else []
        )
    for c in candidates:
        if c.is_file():
            return c
    raise BioLangError(f"cannot resolve import {name!r} (looked in {[str(c) for c in candidates]})")


def parse(
    text: str, name_hint: str = "unnamed", base_dir: Path | None = None, _seen: frozenset[str] = frozenset()
) -> Module:
    directives, blocks = _parse_blocks(text.splitlines())
    name = next((v for k, v in directives if k == "module"), name_hint)
    imports = [v for k, v in directives if k == "import"]
    module = Module(name=name, imports=imports)
    for imp in imports:
        path = resolve_import(imp, base_dir)
        key = str(path.resolve())
        if key in _seen:
            raise BioLangError(f"circular import of {imp!r}")
        module.merge(parse(path.read_text(), name_hint=path.stem, base_dir=path.parent, _seen=_seen | {key}))
    for b in blocks:
        _compile_block(b, module)
    if not blocks and not imports:
        raise BioLangError("empty source: no module, entities or rules found")
    _check_references(module)
    return module


def parse_file(path: str | Path) -> Module:
    path = Path(path)
    return parse(path.read_text(), name_hint=path.stem, base_dir=path.parent)
