"""Ligand-receptor interactions from CellPhoneDB (task 2.6).

Produces BioIR `binds` rules between gene symbols so that the communication
protocol between cells is an executable rule set, not a spreadsheet.
Complexes (e.g. integrin_a2b1_complex) are expanded to their subunits using
the `interactors` column when possible.
"""

from __future__ import annotations

import csv
from pathlib import Path

from genomeos.ir import Action, Entity, Evidence, EvidenceKind, Module, Rule

CPDB_EVIDENCE = Evidence(EvidenceKind.CURATED, "CellPhoneDB v5 interaction_input.csv")


class LigandReceptorTable:
    def __init__(self) -> None:
        self.pairs: list[
            tuple[str, str, str, str]
        ] = []  # (partner_a_symbol, partner_b_symbol, classification, source)
        self.symbol_of_uniprot: dict[str, str] = {}

    @classmethod
    def from_cellphonedb(cls, interactions: str | Path, genes: str | Path) -> LigandReceptorTable:
        t = cls()
        with open(genes, newline="") as fh:
            for row in csv.DictReader(fh):
                t.symbol_of_uniprot[row["uniprot"]] = row["hgnc_symbol"]
        with open(interactions, newline="") as fh:
            for row in csv.DictReader(fh):
                a, b = row["partner_a"], row["partner_b"]
                inter = row.get("interactors", "")
                if "-" in inter and "+" not in inter.split("-")[0]:
                    left, _, right = inter.partition("-")
                    a_syms = left.split("+")
                    b_syms = right.split("+")
                else:
                    a_syms = [t.symbol_of_uniprot.get(a, a)]
                    b_syms = [t.symbol_of_uniprot.get(b, b)]
                for sa in a_syms:
                    for sb in b_syms:
                        if sa and sb and sa != sb:
                            t.pairs.append((sa, sb, row.get("classification", ""), row.get("source", "")))
        return t

    def __len__(self) -> int:
        return len(self.pairs)

    def partners(self, symbol: str) -> list[tuple[str, str]]:
        out = []
        for a, b, cls_, _ in self.pairs:
            if a == symbol:
                out.append((b, cls_))
            elif b == symbol:
                out.append((a, cls_))
        return out

    def to_module(self, name: str = "bio.systems.ligand_receptor") -> Module:
        m = Module(name=name)
        seen: set[tuple[str, str]] = set()
        for a, b, cls_, source in self.pairs:
            for sym in (a, b):
                if sym not in m.entities:
                    m.add(Entity(id=sym, kind="protein", evidence=CPDB_EVIDENCE, confidence=0.8))
            if (a, b) in seen:
                continue
            seen.add((a, b))
            ev = Evidence(EvidenceKind.CURATED, f"CellPhoneDB: {source}" if source else "CellPhoneDB")
            m.rules.append(
                Rule(
                    id=f"{a} binds {b}",
                    source=a,
                    action=Action.BIND,
                    target=b,
                    when={"classification": cls_ or "any"},
                    evidence=ev,
                    confidence=0.8,
                )
            )
        return m
