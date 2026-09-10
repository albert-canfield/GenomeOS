"""Gene Ontology annotations (GAF 2.x) for human: gene symbol <-> GO terms."""

from __future__ import annotations

import gzip
from pathlib import Path

from .obo import Ontology


class GeneOntologyAnnotations:
    def __init__(self, direct: dict[str, set[str]], ontology: Ontology) -> None:
        self.direct = direct  # symbol -> directly annotated GO ids (NOT-qualified excluded)
        self.ontology = ontology
        self._closure: dict[str, frozenset[str]] = {}
        self._members: dict[str, frozenset[str]] = {}

    @classmethod
    def from_gaf(
        cls, gaf: str | Path, ontology: Ontology, taxon: str = "taxon:9606"
    ) -> GeneOntologyAnnotations:
        direct: dict[str, set[str]] = {}
        opener = gzip.open if str(gaf).endswith(".gz") else open
        with opener(gaf, "rt") as fh:
            for line in fh:
                if line.startswith("!"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 13 or "NOT" in f[3].split("|"):
                    continue
                if taxon and not f[12].startswith(taxon):
                    continue
                direct.setdefault(f[2], set()).add(f[4])
        return cls(direct, ontology)

    def terms(self, symbol: str) -> frozenset[str]:
        """GO terms of a gene including all ancestors."""
        if symbol not in self._closure:
            out: set[str] = set()
            for t in self.direct.get(symbol, ()):
                out |= self.ontology.ancestors(t)
            self._closure[symbol] = frozenset(out)
        return self._closure[symbol]

    def members(self, term_id: str) -> frozenset[str]:
        """Genes annotated to a term or any of its descendants."""
        if term_id not in self._members:
            wanted = self.ontology.descendants(term_id)
            self._members[term_id] = frozenset(s for s, ts in self.direct.items() if ts & wanted)
        return self._members[term_id]

    def has(self, symbol: str, term_id: str) -> bool:
        return term_id in self.terms(symbol)

    def __len__(self) -> int:
        return len(self.direct)
