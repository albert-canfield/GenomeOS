"""Design budgets from the libraries: what a cell type needs, and what it costs.

Answers questions like "how much of the genome does a neuron need?" by
composing library memberships (data-driven, from GO/Reactome) with the gene
table (spans) and the anatomy budgets (compact / dense / sparse). The
answer is a parts list with counts and base-pair budgets, each line with its
evidence. It is a design estimate, never a claim that the list is sufficient.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

from genomeos.lib import LIBRARIES, KnowledgeBase

GENES_TSV = Path("data/cache/gencode_genes.tsv")
ESSENTIAL_FILE = Path(__file__).resolve().parent / "lib" / "data" / "essential_cegv2.txt"
ESSENTIAL_EVIDENCE = "Hart et al. 2017, G3 7:2719 (CEGv2: 684 core-essential genes from CRISPR screens)"


def essential_genes() -> set[str]:
    """Genes every human cell line needs to survive (CEGv2), shipped with the package."""
    return {line.strip() for line in ESSENTIAL_FILE.read_text().splitlines() if line.strip()}


# what every cell needs, by library; plus the identity libraries chosen per cell type
CORE = [lib for lib in LIBRARIES if lib.startswith("core.")]
IDENTITY = {
    "neuron": ["blueprint.organ_nervous", "systems.nervous", "blueprint.signalling_toolkit"],
    "cardiomyocyte": [
        "blueprint.organ_heart",
        "systems.nervous",
        "core.energy",
        "blueprint.signalling_toolkit",
    ],
    "hepatocyte": [
        "blueprint.organ_endoderm_gut",
        "systems.metabolic_homeostasis",
        "blueprint.signalling_toolkit",
    ],
    "fibroblast": ["systems.extracellular_matrix_adhesion", "blueprint.signalling_toolkit"],
    "neural_progenitor": [
        "blueprint.pluripotency",
        "timer.stem_cell_niches",
        "blueprint.organ_nervous",
        "blueprint.signalling_toolkit",
    ],
}
# the master regulators that specify the identity (literature; the switch, not the machinery)
MASTER = {
    "neuron": ("NEUROG2", "ASCL1", "NEUROD1", "MYT1L", "POU3F2"),
    "neural_progenitor": ("SOX2", "PAX6", "NES", "HES1"),
    "cardiomyocyte": ("NKX2-5", "GATA4", "TBX5", "MEF2C", "HAND2"),
    "hepatocyte": ("HNF4A", "FOXA2", "HNF1A"),
    "fibroblast": ("PRRX1", "TWIST1"),
}
MASTER_EVIDENCE = {
    "neuron": (
        "Vierbuchen et al. 2010 (Ascl1/Brn2/Myt1l convert fibroblasts to neurons); "
        "Zhang et al. 2013 (NGN2 alone converts stem cells to neurons in days)"
    ),
    "neural_progenitor": "SOX2/PAX6 neural induction; NES marks progenitors",
    "cardiomyocyte": "Ieda et al. 2010 (Gata4/Mef2c/Tbx5 reprogramming)",
    "hepatocyte": "Huang et al. 2011 (Hnf4a/Foxa reprogramming)",
    "fibroblast": "mesenchymal identity factors",
}


@dataclass(slots=True)
class DesignBudget:
    cell_type: str
    essential_genes: int
    core_genes: int
    identity_genes: int
    total_genes: int
    minimal_genes: int  # essential + master regulators + identity libraries
    master_regulators: tuple[str, ...]
    span_bp: int  # sum of genomic gene spans (human layout)
    budgets_bp: dict[str, int]  # compact / dense / sparse estimates
    libraries: dict[str, int]
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "cell_type": self.cell_type,
            "core_genes": self.core_genes,
            "identity_genes": self.identity_genes,
            "total_genes": self.total_genes,
            "master_regulators": list(self.master_regulators),
            "human_layout_span_bp": self.span_bp,
            "budgets_bp": self.budgets_bp,
            "libraries": self.libraries,
            "notes": self.notes,
        }


def _spans() -> dict[str, int]:
    out: dict[str, int] = {}
    if GENES_TSV.exists():
        with open(GENES_TSV) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row["gene_type"] == "protein_coding":
                    out[row["symbol"]] = int(row["end"]) - int(row["start"])
    return out


def design(cell_type: str, kb: KnowledgeBase | None = None) -> DesignBudget:
    if cell_type not in IDENTITY:
        raise KeyError(f"unknown cell type {cell_type!r}; known: {sorted(IDENTITY)}")
    kb = kb or KnowledgeBase()
    core: set[str] = set()
    per_lib: dict[str, int] = {}
    for lib in CORE:
        m = kb.members(LIBRARIES[lib])
        per_lib[lib] = len(m)
        core |= m
    identity: set[str] = set()
    for lib in IDENTITY[cell_type]:
        m = kb.members(LIBRARIES[lib])
        per_lib[lib] = len(m)
        identity |= m
    identity -= core
    essential = essential_genes()
    masters = set(MASTER[cell_type])
    total = core | identity | masters
    minimal = essential | identity | masters
    spans = _spans()
    span = sum(spans.get(g, 0) for g in total)
    n = len(minimal)
    # budgets for the *minimal* set, from the measured anatomies (bases per coding gene incl. spacing)
    budgets = {
        "compact (mtDNA-like, 1.3 kb/gene)": n * 1_300,
        "dense (worm-like, 5 kb/gene)": n * 5_000,
        "sparse (human-like, 150 kb/gene)": n * 150_000,
    }
    notes = [
        f"master regulators: {', '.join(MASTER[cell_type])} [{MASTER_EVIDENCE[cell_type]}]",
        f"essential = {ESSENTIAL_EVIDENCE}",
        "core = union of the 13 core libraries (data-driven GO/Reactome membership): what all cells use, "
        "a superset of what they strictly need",
        "identity = chosen blueprint/systems libraries minus core; minimal = essential + identity + masters",
        "a eukaryotic cell also needs parts the genome does not encode as genes: rRNA/tRNA arrays, "
        "centromeres, telomeres, replication origins, and the maternal machinery of the first cell",
    ]
    return DesignBudget(
        cell_type,
        len(essential),
        len(core),
        len(identity),
        len(total),
        n,
        MASTER[cell_type],
        span,
        budgets,
        per_lib,
        notes,
    )
