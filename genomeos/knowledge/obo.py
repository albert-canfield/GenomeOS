"""Minimal OBO parser for GO, Cell Ontology and Uberon (is_a and part_of)."""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cache
from pathlib import Path


@dataclass(slots=True)
class Term:
    id: str
    name: str = ""
    namespace: str = ""
    definition: str = ""
    is_a: list[str] = field(default_factory=list)
    part_of: list[str] = field(default_factory=list)
    synonyms: list[str] = field(default_factory=list)
    obsolete: bool = False


class Ontology:
    def __init__(self, terms: dict[str, Term], title: str = "") -> None:
        self.terms = terms
        self.title = title
        self._children: dict[str, set[str]] | None = None

    @classmethod
    def from_obo(cls, path: str | Path) -> Ontology:
        terms: dict[str, Term] = {}
        cur: Term | None = None
        in_term = False
        with open(path, encoding="utf-8") as fh:
            for raw in fh:
                line = raw.strip()
                if line == "[Term]":
                    in_term = True
                    cur = None
                    continue
                if line.startswith("["):
                    in_term = False
                    cur = None
                    continue
                if not in_term or not line:
                    continue
                key, _, val = line.partition(": ")
                if key == "id":
                    cur = Term(id=val)
                    terms[val] = cur
                elif cur is None:
                    continue
                elif key == "name":
                    cur.name = val
                elif key == "namespace":
                    cur.namespace = val
                elif key == "def":
                    cur.definition = val.split('"')[1] if '"' in val else val
                elif key == "is_a":
                    cur.is_a.append(val.split(" ")[0])
                elif key == "relationship" and val.startswith("part_of "):
                    cur.part_of.append(val.split(" ")[1])
                elif key == "synonym":
                    cur.synonyms.append(val.split('"')[1] if '"' in val else val)
                elif key == "is_obsolete" and val == "true":
                    cur.obsolete = True
        return cls(terms, Path(path).name)

    def __len__(self) -> int:
        return len(self.terms)

    def __contains__(self, term_id: str) -> bool:
        return term_id in self.terms

    def name(self, term_id: str) -> str:
        t = self.terms.get(term_id)
        return t.name if t else term_id

    @cache  # noqa: B019  (ontology objects are few and long-lived)
    def ancestors(self, term_id: str) -> frozenset[str]:
        """All is_a / part_of ancestors, including the term itself."""
        seen: set[str] = set()
        stack = [term_id]
        while stack:
            t = stack.pop()
            if t in seen:
                continue
            seen.add(t)
            term = self.terms.get(t)
            if term:
                stack.extend(term.is_a)
                stack.extend(term.part_of)
        return frozenset(seen)

    def descendants(self, term_id: str) -> set[str]:
        if self._children is None:
            self._children = {}
            for t in self.terms.values():
                for p in t.is_a + t.part_of:
                    self._children.setdefault(p, set()).add(t.id)
        out: set[str] = set()
        stack = [term_id]
        while stack:
            t = stack.pop()
            if t in out:
                continue
            out.add(t)
            stack.extend(self._children.get(t, ()))
        return out

    def search(self, text: str, limit: int = 50) -> list[Term]:
        q = text.lower()
        hits = [
            t
            for t in self.terms.values()
            if not t.obsolete and (q in t.name.lower() or any(q in s.lower() for s in t.synonyms))
        ]
        hits.sort(key=lambda t: (not t.name.lower().startswith(q), len(t.name)))
        return hits[:limit]
