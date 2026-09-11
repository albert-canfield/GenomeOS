"""BioLang v0.3 parser: source text -> BioIR Module.

Grammar (see docs/BIOLANG-v0.2.md and docs/BIOLANG-v0.3.md):

    module <dotted.name>
    import <dotted.name> | <relative/path.bio>

    gene <Id> { <props>  transcript <Id> { exons: ...; cds: ... } ... }
    protein <Id> { <props> }
    region <Id> { locus: chr7:1000-2000; role: unknown }
    cell_type <Id> { expresses: a, b; parent: X; ontology: CL:0000540 }
    event <Id> { rate: 0.5 /yr; when: cell_type = X; effect: telomere_bp -= 70 bp; ... }
    rule <Source> activates|inhibits|produces|binds|degrades|modifies <Target> { <props> }
    param <name> = <number> [unit] { evidence: ...; confidence: ... }
    domain <Id> { locus: chr21:a-b; genes: A, B; boundaries: E1, E2 }
    organism <Id> { species: ...; genome: ...; tempo: 1.0; cell_type: Zygote; factors: A, B; ... }
    stage <Id> { from: 0 min; to: 100 min }
    timer <Id> { duration: 20 min; sd: 2; lengthening: 1.1; when: lineage = AB }
    signal <Id> { mode: contact; ligand: APX-1; receptor: GLP-1; from: cell = P2; to: cell = ABp; sets: N }
    decision <Id> { action: divide|differentiate|migrate|quiesce|die; when: ...; daughters: A, B; to: T; ... }
    experiment <Id> { knockout: POP-1; until: 800 min; expect: "..."; assert: ... }

Properties are `key: value`, one per line or `;`-separated; a block may sit on
one line. Blocks may nest (transcript inside gene). `#` starts a comment.
Common keys on any block: evidence, confidence.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

from genomeos.coords import Locus
from genomeos.ir import (
    DECISION_ACTIONS,
    UNKNOWN,
    Action,
    CellType,
    Decision,
    Domain,
    Effect,
    Event,
    Evidence,
    EvidenceKind,
    Experiment,
    Gene,
    Module,
    Organism,
    Parameter,
    Protein,
    Region,
    RegulatoryElement,
    Rule,
    Signal,
    Stage,
    Timer,
    Transcript,
    to_minutes,
)

_ACTIONS = {a.value: a for a in Action}
_KINDS = (
    "gene",
    "protein",
    "region",
    "element",
    "rule",
    "param",
    "cell_type",
    "event",
    "transcript",
    "domain",
    "organism",
    "stage",
    "timer",
    "signal",
    "decision",
    "experiment",
)
_REPEATABLE = ("effect", "assert", "observe")
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
            if key in _REPEATABLE and key in block.props:
                block.props[key] += " ; " + v.strip()  # several effect / assert / observe lines per block
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


def _quantity(value: str, key: str, line_no: int, default_unit: str = "min") -> tuple[float, str]:
    parts = value.split()
    return _float(parts[0], key, line_no), (parts[1] if len(parts) > 1 else default_unit)


def _arrows(value: str, line_no: int) -> dict[str, str]:
    """`POP-1 -> MS, PIE-1 -> P2` -> {"POP-1": "MS", "PIE-1": "P2"}."""
    out: dict[str, str] = {}
    for clause in _list(value):
        a, sep, b = clause.partition("->")
        if not sep:
            raise BioLangError(f"line {line_no}: expected 'factor -> daughter', got {clause!r}")
        out[a.strip()] = b.strip()
    return out


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
        if "accession" in p:
            pr.accession = p["accession"]
        if "isoforms" in p:
            pr.isoforms = _list(p["isoforms"])
        if "domains" in p:
            pr.domains = [{"name": d} for d in _list(p["domains"])]
        if "pathways" in p:
            pr.pathways = _list(p["pathways"])
        if "interactions" in p:
            pr.interactions = _list(p["interactions"])
        if "structures" in p:
            pr.structures = [
                {"source": "PDB" if not x.startswith("AF-") else "AlphaFold", "id": x}
                for x in _list(p["structures"])
            ]
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
    elif b.kind == "element":
        el = RegulatoryElement(
            id=b.header,
            kind="regulatory_element",
            evidence=ev,
            confidence=conf,
            cls=p.get("class", "unknown"),
        )
        if "locus" in p:
            el.locus = Locus.parse(p["locus"])
        if "domain" in p:
            el.domain = p["domain"]
        if "targets" in p:
            el.targets = [
                {"gene": t, "basis": p.get("basis", "stated"), "confidence": conf}
                for t in _list(p["targets"])
            ]
        module.add(el)
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
    elif b.kind == "domain":
        dm = Domain(id=b.header, kind="domain", evidence=ev, confidence=conf)
        if "locus" in p:
            dm.locus = Locus.parse(p["locus"])
        dm.genes = _list(p.get("genes", ""))
        dm.boundaries = _list(p.get("boundaries", ""))
        module.add(dm)
    elif b.kind == "signal":
        sg = Signal(id=b.header, kind="signal", evidence=ev, confidence=conf, mode=p.get("mode", "contact"))
        if sg.mode not in ("contact", "gradient", "systemic"):
            raise BioLangError(f"line {b.line}: signal mode must be contact, gradient or systemic")
        sg.ligand, sg.receptor = p.get("ligand", ""), p.get("receptor", "")
        sg.sender, sg.receiver = _parse_when(p.get("from", "")), _parse_when(p.get("to", ""))
        factor, _, value = p.get("sets", "").partition("=")
        sg.sets, sg.value = factor.strip(), value.strip() or "active"
        module.add(sg)
    elif b.kind == "organism":
        if module.organism is not None:
            raise BioLangError(
                f"line {b.line}: a program may declare one organism ({module.organism.name!r} exists)"
            )
        org = Organism(name=b.header, evidence=ev, confidence=conf)
        org.species, org.genome = p.get("species", ""), p.get("genome", "")
        if "tempo" in p:
            org.tempo = _float(p["tempo"], "tempo", b.line)
        org.root = p.get("root", org.root)
        org.resolution = p.get("resolution", org.resolution)
        if "seed" in p:
            org.seed = int(_float(p["seed"], "seed", b.line))
        if org.resolution not in ("cells", "populations"):
            raise BioLangError(f"line {b.line}: resolution must be cells or populations")
        org.cell_type = p.get("cell_type", "")
        org.factors = _list(p.get("factors", ""))
        org.environment = _parse_when(p.get("environment", ""))
        org.observe = [x for part in p.get("observe", "").split(" ; ") for x in _list(part)]
        org.asserts = [x.strip() for x in p.get("assert", "").split(" ; ") if x.strip()]
        org.reference = p.get("reference", "")
        module.organism = org
    elif b.kind == "experiment":
        ex = Experiment(name=b.header, evidence=ev, confidence=conf)
        ex.knockouts = _list(p.get("knockout", ""))
        ex.adds = _list(p.get("add", ""))
        ex.environment = _parse_when(p.get("environment", ""))
        if "until" in p:
            val, unit = _quantity(p["until"], "until", b.line)
            ex.until = to_minutes(val, unit)
        ex.asserts = [x.strip() for x in p.get("assert", "").split(" ; ") if x.strip()]
        ex.expect = p.get("expect", "")
        module.experiments.append(ex)
    elif b.kind == "stage":
        st = Stage(name=b.header, evidence=ev, confidence=conf)
        if "from" in p:
            st.start, st.unit = _quantity(p["from"], "from", b.line)
        if "to" in p:
            end, unit = _quantity(p["to"], "to", b.line)
            st.end = to_minutes(end, unit) / to_minutes(1.0, st.unit)
        module.stages.append(st)
    elif b.kind == "timer":
        if "duration" not in p:
            raise BioLangError(f"line {b.line}: timer {b.header!r} needs a duration")
        dur, unit = _quantity(p["duration"], "duration", b.line)
        tm = Timer(name=b.header, duration=dur, unit=unit, evidence=ev, confidence=conf)
        to_minutes(1.0, unit)  # validates the unit
        if "sd" in p:
            tm.sd = _float(p["sd"], "sd", b.line)
        if "lengthening" in p:
            tm.lengthening = _float(p["lengthening"], "lengthening", b.line)
        if "when" in p:
            tm.when = _parse_when(p["when"])
        module.timers.append(tm)
    elif b.kind == "decision":
        action = p.get("action", "")
        if action not in DECISION_ACTIONS:
            raise BioLangError(
                f"line {b.line}: decision action must be one of {DECISION_ACTIONS}; got {action!r}"
            )
        dc = Decision(id=b.header, action=action, evidence=ev, confidence=conf)
        if "when" in p:
            dc.when = _parse_when(p["when"])
        dc.daughters = _list(p.get("daughters", ""))
        dc.to, dc.name, dc.timer = p.get("to", ""), p.get("name", ""), p.get("timer", "")
        if "asymmetric" in p:
            dc.asymmetric = _arrows(p["asymmetric"], b.line)
        if "lineages" in p:
            dc.lineages = _parse_when(p["lineages"])
        if "after" in p:
            val, unit = _quantity(p["after"], "after", b.line)
            dc.after = to_minutes(val, unit)
        if "fraction" in p:
            dc.fraction = _float(p["fraction"], "fraction", b.line)
            top = math.inf if action == "divide" else 1.0  # a population may more than double per step
            if not 0.0 <= dc.fraction <= top:
                raise BioLangError(
                    f"line {b.line}: fraction must be within 0..{'1' if top == 1.0 else 'any'}"
                )
        if action == "divide" and len(dc.daughters) not in (0, 2):
            raise BioLangError(f"line {b.line}: a division names two daughters or none")
        if action == "differentiate" and not dc.to:
            raise BioLangError(f"line {b.line}: differentiate needs 'to: <cell_type>'")
        module.decisions.append(dc)
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
    cell_types = {c.id for c in module.cell_types()}
    timers = {t.name for t in module.timers}
    for d in module.decisions:
        if d.to and d.to not in cell_types:
            raise BioLangError(f"decision {d.id!r} differentiates to undeclared cell_type {d.to!r}")
        if d.timer and d.timer not in timers:
            raise BioLangError(f"decision {d.id!r} waits on undeclared timer {d.timer!r}")
    org = module.organism
    if org and org.cell_type and org.cell_type not in cell_types:
        raise BioLangError(f"organism {org.name!r} starts as undeclared cell_type {org.cell_type!r}")


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
    text: str,
    name_hint: str = "unnamed",
    base_dir: Path | None = None,
    _seen: frozenset[str] = frozenset(),
    _done: set[str] | None = None,
) -> Module:
    directives, blocks = _parse_blocks(text.splitlines())
    name = next((v for k, v in directives if k == "module"), name_hint)
    imports = [v for k, v in directives if k == "import"]
    module = Module(name=name, imports=imports)
    done = set() if _done is None else _done  # files already merged into this program (diamond imports)
    for imp in imports:
        path = resolve_import(imp, base_dir)
        key = str(path.resolve())
        if key in _seen:
            raise BioLangError(f"circular import of {imp!r}")
        if key in done:
            continue
        done.add(key)
        module.merge(
            parse(
                path.read_text(), name_hint=path.stem, base_dir=path.parent, _seen=_seen | {key}, _done=done
            )
        )
    for b in blocks:
        _compile_block(b, module)
    if not blocks and not imports:
        raise BioLangError("empty source: no module, entities or rules found")
    if _done is None:  # references are checked once, over the whole program
        _check_references(module)
    return module


def parse_file(path: str | Path) -> Module:
    path = Path(path)
    return parse(path.read_text(), name_hint=path.stem, base_dir=path.parent)
