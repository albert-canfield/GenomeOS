"""Data-driven library membership (task 1.5).

The hand-written catalogue names GO terms and Reactome pathways for each
library; this module computes the actual member genes from the downloaded
knowledge files and verifies the catalogue's representative genes against
them. The catalogue therefore becomes a test, not a source of truth.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path

from genomeos.knowledge import GeneOntologyAnnotations, Ontology, Reactome

from .catalog import LIBRARIES, Library

DEFAULT_KNOWLEDGE = Path("data/knowledge")
DEFAULT_GENES = Path("data/cache/gencode_genes.tsv")


@dataclass(slots=True)
class Verification:
    library: str
    checked: tuple[str, ...]
    missing: tuple[str, ...]
    unknown_terms: tuple[str, ...]
    members: int

    @property
    def ok(self) -> bool:
        return not self.missing and not self.unknown_terms


@dataclass
class KnowledgeBase:
    knowledge_dir: Path = DEFAULT_KNOWLEDGE
    genes_tsv: Path = DEFAULT_GENES
    _symbol_of: dict[str, str] = field(default_factory=dict, init=False)

    @classmethod
    def available(cls, knowledge_dir: Path = DEFAULT_KNOWLEDGE, genes_tsv: Path = DEFAULT_GENES) -> bool:
        """GO files present, plus either the Reactome mapping or the distilled membership summary."""
        go_ok = all((knowledge_dir / f).exists() for f in ("go-basic.obo", "goa_human.gaf.gz"))
        return go_ok and ((knowledge_dir / "Ensembl2Reactome.txt").exists() or cls.distilled() is not None)

    @staticmethod
    def distilled() -> dict | None:
        """The distilled membership: data/results if present, else the copy shipped in the package."""
        import gzip
        import json

        from genomeos.results import load_result

        r = load_result("library_members")
        if r is not None:
            return r
        packaged = Path(__file__).resolve().parent / "data" / "members.json.gz"
        if packaged.exists():
            with gzip.open(packaged, "rt") as fh:
                return json.load(fh)
        return None

    @cached_property
    def go(self) -> Ontology:
        return Ontology.from_obo(self.knowledge_dir / "go-basic.obo")

    @cached_property
    def annotations(self) -> GeneOntologyAnnotations:
        return GeneOntologyAnnotations.from_gaf(self.knowledge_dir / "goa_human.gaf.gz", self.go)

    @cached_property
    def reactome(self) -> Reactome:
        return Reactome.from_file(self.knowledge_dir / "Ensembl2Reactome.txt")

    @cached_property
    def symbol_of(self) -> dict[str, str]:
        """Unversioned Ensembl gene id -> HGNC symbol (from the GENCODE gene table)."""
        out: dict[str, str] = {}
        if self.genes_tsv.exists():
            with open(self.genes_tsv) as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    out[row["gene_id"].split(".")[0]] = row["symbol"]
        return out

    # ---- membership ------------------------------------------------------

    def go_members(self, lib: Library) -> set[str]:
        out: set[str] = set()
        for t in lib.go_terms:
            out |= self.annotations.members(t)
        return out

    def reactome_members(self, lib: Library) -> set[str]:
        out: set[str] = set()
        if not (self.knowledge_dir / "Ensembl2Reactome.txt").exists():
            # the raw mapping was discarded after distillation: use the summary (GO ∪ Reactome members)
            d = self.distilled()
            return set(d["members"].get(lib.id, [])) - self.go_members(lib) if d else set()
        for needle in lib.reactome:
            pid = self.reactome.id_of(needle)
            if pid:
                out |= {self.symbol_of.get(g, g) for g in self.reactome.pathway_genes.get(pid, ())}
        return out

    def members(self, lib: Library) -> set[str]:
        return self.go_members(lib) | self.reactome_members(lib)

    def verify(self, lib: Library) -> Verification:
        unknown = tuple(t for t in lib.go_terms if t not in self.go or self.go.terms[t].obsolete)
        if (self.knowledge_dir / "Ensembl2Reactome.txt").exists():
            unknown += tuple(n for n in lib.reactome if self.reactome.id_of(n) is None)
        checked = tuple(g for g in lib.genes if g not in lib.noncoding)
        members = self.members(lib)
        missing = tuple(g for g in checked if g not in members)
        return Verification(lib.id, checked, missing, unknown, len(members))

    def verify_all(self) -> list[Verification]:
        return [self.verify(lib) for lib in LIBRARIES.values() if lib.go_terms or lib.reactome]

    def libraries_of(self, symbol: str) -> list[str]:
        """Which libraries a gene belongs to, by data."""
        return [lib.id for lib in LIBRARIES.values() if symbol in self.members(lib)]
