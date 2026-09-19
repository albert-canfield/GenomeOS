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


# The developmental signalling pathways worth a `bio.std` module: ligand -> receptors, with the reviews
# that make each pair textbook; CellPhoneDB supplies the pair and its curation, the review the role.
DEVELOPMENTAL = {
    "DLL1": ("Notch", "Bray 2016, Nat Rev Mol Cell Biol 17:722"),
    "DLL4": ("Notch", "Bray 2016"),
    "JAG1": ("Notch", "Bray 2016"),
    "JAG2": ("Notch", "Bray 2016"),
    "WNT3A": ("Wnt", "Clevers & Nusse 2012, Cell 149:1192"),
    "WNT5A": ("Wnt", "Clevers & Nusse 2012"),
    "WNT1": ("Wnt", "Clevers & Nusse 2012"),
    "BMP4": ("BMP", "Wang et al. 2014, Genes Dis 1:87"),
    "BMP2": ("BMP", "Wang et al. 2014"),
    "BMP7": ("BMP", "Wang et al. 2014"),
    "NODAL": ("Nodal", "Schier 2009, Cold Spring Harb Perspect Biol 1:a003459"),
    "TGFB1": ("TGF-beta", "Massagué 2012, Nat Rev Mol Cell Biol 13:616"),
    "FGF8": ("FGF", "Ornitz & Itoh 2015, WIREs Dev Biol 4:215"),
    "FGF2": ("FGF", "Ornitz & Itoh 2015"),
    "FGF4": ("FGF", "Ornitz & Itoh 2015"),
    "FGF10": ("FGF", "Ornitz & Itoh 2015"),
    "SHH": ("Hedgehog", "Briscoe & Thérond 2013, Nat Rev Mol Cell Biol 14:416"),
    "IHH": ("Hedgehog", "Briscoe & Thérond 2013"),
    "EGF": ("EGF", "Schneider & Wolf 2009, J Cell Physiol 218:460"),
    "CXCL12": ("Chemokine", "Nagasawa 2014, Int Immunol 26:589"),
    "KITLG": ("Kit", "Lennartsson & Rönnstrand 2012, Physiol Rev 92:1619"),
    "EPO": ("Erythropoietin", "Bunn 2013, Cold Spring Harb Perspect Med 3:a011619"),
    "THPO": ("Thrombopoietin", "Kaushansky 2005, J Clin Invest 115:3339"),
    "CSF3": ("G-CSF", "Panopoulos & Watowich 2008, Cytokine 42:277"),
    "CSF1": ("M-CSF", "Stanley & Chitu 2014, Cold Spring Harb Perspect Biol 6:a021857"),
    "IL7": ("IL-7", "Mazzucchelli & Durum 2007, Nat Rev Immunol 7:144"),
    "GDF11": ("GDF", "Walker et al. 2016, Circ Res 118:1125"),
    "RSPO1": ("Wnt", "de Lau et al. 2014, Genes Dev 28:305"),
}


def to_std_signalling(table: LigandReceptorTable, module_name: str = "bio.std.signalling") -> str:
    """The developmental ligand-receptor pairs as BioLang `signal` blocks (contact mode; the receiver gains
    `<Pathway>_signal = received`), generated from the CellPhoneDB table so the pairs are curated."""
    lines = [
        "# bio.std.signalling — developmental signalling pairs as BioLang signals, generated from",
        "# CellPhoneDB (curated ligand-receptor pairs) restricted to the pathways of development and",
        "# haematopoiesis; the receiver of each gains `<Pathway>_signal = received`. A program binds",
        "# the pair to cells with its own `from:` and `to:` conditions, e.g. by importing and restating.",
        f"module {module_name}",
        "",
    ]
    seen: set[tuple[str, str]] = set()
    for a, b, _cls, source in table.pairs:
        for ligand, receptor in ((a, b), (b, a)):
            if ligand in DEVELOPMENTAL and (ligand, receptor) not in seen and receptor not in DEVELOPMENTAL:
                seen.add((ligand, receptor))
                pathway, review = DEVELOPMENTAL[ligand]
                factor = pathway.replace("-", "").replace(" ", "") + "_signal"
                cite = f"CellPhoneDB v5 ({source}); {review}" if source else f"CellPhoneDB v5; {review}"
                lines.append(
                    f"signal {ligand}_{receptor} {{ mode: contact; ligand: {ligand}; receptor: {receptor}; "
                    f'sets: {factor} = received; evidence: curated "{cite}"; confidence: 0.8 }}'
                )
    return "\n".join(lines) + "\n"
