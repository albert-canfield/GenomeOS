"""Cell types from the Cell Ontology (task 2.2)."""

from __future__ import annotations

from pathlib import Path

from genomeos.ir import Entity, Evidence, EvidenceKind, Module

from .obo import Ontology

CL_EVIDENCE = Evidence(EvidenceKind.CURATED, "Cell Ontology (cl-basic.obo)")


class CellTypes:
    def __init__(self, ontology: Ontology) -> None:
        self.ontology = ontology

    @classmethod
    def from_obo(cls, path: str | Path) -> CellTypes:
        return cls(Ontology.from_obo(path))

    def __len__(self) -> int:
        return sum(1 for t in self.ontology.terms.values() if t.id.startswith("CL:") and not t.obsolete)

    def search(self, text: str, limit: int = 20):
        return [t for t in self.ontology.search(text, limit * 3) if t.id.startswith("CL:")][:limit]

    def lineage(self, cl_id: str) -> list[str]:
        """Names from the term up to the root, following is_a only (a readable lineage)."""
        out = []
        cur = cl_id
        seen = set()
        while cur and cur not in seen:
            seen.add(cur)
            term = self.ontology.terms.get(cur)
            if term is None:
                break
            out.append(f"{term.id} {term.name}")
            cur = term.is_a[0] if term.is_a else None
        return out

    def to_module(self, name: str = "bio.parts.cell_types", limit: int | None = None) -> Module:
        m = Module(name=name)
        n = 0
        for t in self.ontology.terms.values():
            if not t.id.startswith("CL:") or t.obsolete:
                continue
            m.add(
                Entity(
                    id=t.id,
                    kind="cell_type",
                    attrs={"name": t.name, "is_a": list(t.is_a), "part_of": list(t.part_of)},
                    evidence=CL_EVIDENCE,
                    confidence=0.9,
                )
            )
            n += 1
            if limit and n >= limit:
                break
        return m
