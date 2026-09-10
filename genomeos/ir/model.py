"""BioIR v0.1 - the Biological Intermediate Representation.

Everything the compiler emits and the VM executes is one of these types.
Three principles are enforced by the types themselves:

1. Every rule carries Evidence and a Confidence.
2. UNKNOWN is a legal value for a role; the compiler must not invent function.
3. Rules are contextual: they apply only when their `when` conditions hold.

Serialisation is plain JSON (to_dict / from_dict) so that BioIR can be produced
by other tools and diffed by humans.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, StrEnum
from typing import Any

from genomeos.coords import Locus


class EvidenceKind(StrEnum):
    EXPERIMENTAL = "experimental"  # measured in the lab
    CURATED = "curated"  # from a curated database (GO, Reactome, ...)
    PREDICTED = "predicted"  # from a computational model (AlphaGenome, ...)
    INFERRED = "inferred"  # by analogy / homology / expert estimate
    NONE = "none"  # no evidence; the value is a placeholder


@dataclass(frozen=True, slots=True)
class Evidence:
    kind: EvidenceKind = EvidenceKind.NONE
    source: str = ""  # citation, database id, model name
    organism: str = "Homo sapiens"
    note: str = ""


Confidence = float  # 0.0 (pure guess) .. 1.0 (settled fact)


class _Unknown:
    """Singleton marking a role, function or value we do not know."""

    _instance: _Unknown | None = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __repr__(self) -> str:
        return "UNKNOWN"


UNKNOWN = _Unknown()


class Action(StrEnum):
    ACTIVATE = "activates"
    INHIBIT = "inhibits"
    MODIFY = "modifies"
    PRODUCE = "produces"
    BIND = "binds"
    DEGRADE = "degrades"


@dataclass(slots=True)
class Entity:
    """Base for anything the VM can hold state for."""

    id: str
    kind: str  # "gene", "protein", "region", ...
    attrs: dict[str, Any] = field(default_factory=dict)
    evidence: Evidence = field(default_factory=Evidence)
    confidence: Confidence = 0.0


@dataclass(slots=True)
class Region(Entity):
    """A genomic interval whose role may be UNKNOWN."""

    locus: Locus | None = None
    role: str | _Unknown = UNKNOWN

    def __post_init__(self) -> None:
        self.kind = "region"


@dataclass(slots=True)
class Gene(Entity):
    symbol: str = ""
    locus: Locus | None = None
    transcripts: list[Transcript] = field(default_factory=list)
    basal_rate: float = 0.0  # transcription with no regulators (a.u./h)

    def __post_init__(self) -> None:
        self.kind = "gene"


@dataclass(slots=True)
class Transcript(Entity):
    gene_id: str = ""
    exons: list[Locus] = field(default_factory=list)
    cds: Locus | None = None  # span of the coding sequence
    cds_segments: list[Locus] = field(default_factory=list)  # CDS pieces in genomic order
    cds_phase: int = 0  # phase of the first CDS segment
    tags: list[str] = field(default_factory=list)  # e.g. GENCODE "cds_start_NF"

    def __post_init__(self) -> None:
        self.kind = "transcript"


@dataclass(slots=True)
class Protein(Entity):
    sequence: str = ""
    half_life_h: float | _Unknown = UNKNOWN

    def __post_init__(self) -> None:
        self.kind = "protein"


@dataclass(slots=True)
class CellType(Entity):
    """A cell type: which genes it expresses (context for rules) and its parent type."""

    name: str = ""
    parent: str = ""
    expresses: list[str] = field(default_factory=list)
    ontology_id: str = ""  # e.g. CL:0000540

    def __post_init__(self) -> None:
        self.kind = "cell_type"


@dataclass(slots=True)
class Effect:
    target: str  # a state variable, e.g. "telomere_bp"
    op: str  # "+=", "-=", "=", "*="
    value: float
    unit: str = ""


@dataclass(slots=True)
class Event:
    """A discrete event (divide, differentiate, die, ...) with a rate, guards and effects."""

    id: str
    rate: float = 0.0
    rate_unit: str = "1/yr"
    when: dict[str, str] = field(default_factory=dict)
    effects: list[Effect] = field(default_factory=list)
    evidence: Evidence = field(default_factory=Evidence)
    confidence: Confidence = 0.0

    def applies(self, context: dict[str, str]) -> bool:
        return all(v == "any" or context.get(k) == v for k, v in self.when.items())


@dataclass(slots=True)
class Parameter:
    """A numeric fact used by the runtime, with provenance attached."""

    name: str
    value: float
    unit: str = ""
    evidence: Evidence = field(default_factory=Evidence)
    confidence: Confidence = 0.0


@dataclass(slots=True)
class Rule:
    """`source ACTION target` under conditions `when`.

    strength: 0..1 for regulatory rules. For ACTIVATE/INHIBIT the runtime uses a
    Hill function with `threshold` (the source level giving half-maximal effect)
    and `hill` (cooperativity).
    """

    id: str
    source: str
    action: Action
    target: str
    strength: float = 1.0
    threshold: float = 1.0
    hill: float = 2.0
    when: dict[str, str] = field(default_factory=dict)
    evidence: Evidence = field(default_factory=Evidence)
    confidence: Confidence = 0.0

    def applies(self, context: dict[str, str]) -> bool:
        for key, wanted in self.when.items():
            if wanted == "any":
                continue
            if context.get(key) != wanted:
                return False
        return True


@dataclass(slots=True)
class Module:
    """A compiled unit: entities plus rules plus parameters. This is what
    BioLang files compile to and what the VM loads."""

    name: str
    entities: dict[str, Entity] = field(default_factory=dict)
    rules: list[Rule] = field(default_factory=list)
    parameters: dict[str, Parameter] = field(default_factory=dict)
    imports: list[str] = field(default_factory=list)
    events: list[Event] = field(default_factory=list)

    def add(self, entity: Entity) -> None:
        if entity.id in self.entities:
            raise ValueError(f"duplicate entity id {entity.id!r} in module {self.name}")
        self.entities[entity.id] = entity

    def genes(self) -> list[Gene]:
        return [e for e in self.entities.values() if isinstance(e, Gene)]

    def proteins(self) -> list[Protein]:
        return [e for e in self.entities.values() if isinstance(e, Protein)]

    def cell_types(self) -> list[CellType]:
        return [e for e in self.entities.values() if isinstance(e, CellType)]

    def merge(self, other: Module) -> None:
        """Import another module's content (ids must not collide)."""
        for e in other.entities.values():
            self.add(e)
        self.rules.extend(other.rules)
        self.events.extend(other.events)
        for k, v in other.parameters.items():
            self.parameters.setdefault(k, v)

    def active_rules(self, context: dict[str, str]) -> list[Rule]:
        """Rules whose `when` matches and whose genes are expressed in the context's cell type."""
        silenced = self.silenced_genes(context)
        return [
            r
            for r in self.rules
            if r.applies(context) and r.source not in silenced and r.target not in silenced
        ]

    def silenced_genes(self, context: dict[str, str]) -> set[str]:
        ct = context.get("cell_type")
        if not ct:
            return set()
        cell = next((c for c in self.cell_types() if c.id == ct or c.name == ct), None)
        if cell is None or not cell.expresses:
            return set()
        expressed = set(cell.expresses)
        return {g.id for g in self.genes() if g.id not in expressed and g.symbol not in expressed}

    def unknowns(self) -> list[Entity]:
        return [e for e in self.entities.values() if isinstance(e, Region) and e.role is UNKNOWN]

    def confidence_report(self) -> dict[str, float]:
        """Mean confidence by entity kind and for rules. Zero if nothing of that kind."""
        by_kind: dict[str, list[float]] = {}
        for e in self.entities.values():
            by_kind.setdefault(e.kind, []).append(e.confidence)
        report = {k: sum(v) / len(v) for k, v in by_kind.items()}
        if self.rules:
            report["rule"] = sum(r.confidence for r in self.rules) / len(self.rules)
        return report

    # ---- serialisation -------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        def conv(obj: Any) -> Any:
            if obj is UNKNOWN:
                return "UNKNOWN"
            if isinstance(obj, Locus):
                return str(obj)
            if isinstance(obj, Enum):
                return obj.value
            if isinstance(obj, dict):
                return {k: conv(v) for k, v in obj.items()}
            if isinstance(obj, list):
                return [conv(v) for v in obj]
            if hasattr(obj, "__dataclass_fields__"):
                d = {f: conv(getattr(obj, f)) for f in obj.__dataclass_fields__}
                d["__type__"] = type(obj).__name__
                return d
            return obj

        return {
            "bioir_version": "0.1",
            "name": self.name,
            "imports": list(self.imports),
            "entities": [conv(e) for e in self.entities.values()],
            "rules": [conv(r) for r in self.rules],
            "events": [conv(e) for e in self.events],
            "parameters": [conv(p) for p in self.parameters.values()],
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Module:
        types = {
            "Gene": Gene,
            "Protein": Protein,
            "Region": Region,
            "Transcript": Transcript,
            "Entity": Entity,
            "CellType": CellType,
        }

        def evid(d: dict) -> Evidence:
            return Evidence(
                EvidenceKind(d.get("kind", "none")),
                d.get("source", ""),
                d.get("organism", "Homo sapiens"),
                d.get("note", ""),
            )

        def locus(v: Any) -> Locus | None:
            return Locus.parse(v) if isinstance(v, str) else None

        def unk(v: Any) -> Any:
            return UNKNOWN if v == "UNKNOWN" else v

        m = cls(name=data["name"], imports=list(data.get("imports", [])))
        for ed in data.get("entities", []):
            t = types[ed.pop("__type__", "Entity")]
            ed["evidence"] = evid(ed.get("evidence", {}))
            for key in ("locus", "cds"):
                if key in ed:
                    ed[key] = locus(ed[key])
            if "exons" in ed:
                ed["exons"] = [Locus.parse(x) for x in ed["exons"]]
            if "cds_segments" in ed:
                ed["cds_segments"] = [Locus.parse(x) for x in ed["cds_segments"]]
            if "transcripts" in ed:
                ed["transcripts"] = []  # transcripts are top-level entities too
            for key in ("role", "half_life_h"):
                if key in ed:
                    ed[key] = unk(ed[key])
            m.add(t(**ed))
        for rd in data.get("rules", []):
            rd.pop("__type__", None)
            rd["evidence"] = evid(rd.get("evidence", {}))
            rd["action"] = Action(rd["action"])
            m.rules.append(Rule(**rd))
        for evd in data.get("events", []):
            evd.pop("__type__", None)
            evd["evidence"] = evid(evd.get("evidence", {}))
            evd["effects"] = [
                Effect(**{k: v for k, v in ef.items() if k != "__type__"}) for ef in evd.get("effects", [])
            ]
            m.events.append(Event(**evd))
        for pd in data.get("parameters", []):
            pd.pop("__type__", None)
            pd["evidence"] = evid(pd.get("evidence", {}))
            m.parameters[pd["name"]] = Parameter(**pd)
        return m
