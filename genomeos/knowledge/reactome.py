"""Reactome pathway membership.

Ensembl2Reactome.txt maps genes to lowest-level pathways only, so membership
is propagated up the hierarchy (ReactomePathwaysRelation.txt) to give every
pathway, including top-level ones such as "Signaling by WNT", its full gene set.
"""

from __future__ import annotations

from pathlib import Path


class Reactome:
    def __init__(self) -> None:
        self.gene_pathways: dict[str, set[str]] = {}  # ENSG (unversioned) -> pathway ids (all levels)
        self.pathway_name: dict[str, str] = {}
        self.pathway_genes: dict[str, set[str]] = {}
        self.parents: dict[str, set[str]] = {}

    @classmethod
    def from_files(
        cls,
        ensembl2reactome: str | Path,
        relations: str | Path | None = None,
        names: str | Path | None = None,
        species: str = "Homo sapiens",
        prefix: str = "R-HSA-",
    ) -> Reactome:
        r = cls()
        if names and Path(names).exists():
            with open(names, encoding="utf-8") as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) >= 3 and f[2] == species:
                        r.pathway_name[f[0]] = f[1]
        if relations and Path(relations).exists():
            with open(relations, encoding="utf-8") as fh:
                for line in fh:
                    parent, _, child = line.rstrip("\n").partition("\t")
                    if parent.startswith(prefix):
                        r.parents.setdefault(child, set()).add(parent)
        with open(ensembl2reactome, encoding="utf-8") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 6 or f[5] != species or not f[0].startswith("ENSG"):
                    continue
                gene = f[0].split(".")[0]
                r.pathway_name.setdefault(f[1], f[3])
                for pid in r._with_ancestors(f[1]):
                    r.gene_pathways.setdefault(gene, set()).add(pid)
                    r.pathway_genes.setdefault(pid, set()).add(gene)
        return r

    @classmethod
    def from_file(cls, path: str | Path, species: str = "Homo sapiens") -> Reactome:
        d = Path(path).parent
        return cls.from_files(path, d / "ReactomePathwaysRelation.txt", d / "ReactomePathways.txt", species)

    def _with_ancestors(self, pid: str) -> set[str]:
        out: set[str] = set()
        stack = [pid]
        while stack:
            x = stack.pop()
            if x in out:
                continue
            out.add(x)
            stack.extend(self.parents.get(x, ()))
        return out

    def id_of(self, name: str) -> str | None:
        for pid, n in self.pathway_name.items():
            if n.lower() == name.lower():
                return pid
        return None

    def search(self, text: str, limit: int = 30) -> list[tuple[str, str, int]]:
        q = text.lower()
        hits = [
            (pid, name, len(self.pathway_genes.get(pid, ())))
            for pid, name in self.pathway_name.items()
            if q in name.lower()
        ]
        hits.sort(key=lambda x: -x[2])
        return hits[:limit]

    def __len__(self) -> int:
        return len(self.pathway_name)
