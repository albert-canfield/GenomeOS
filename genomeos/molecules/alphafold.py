"""AlphaFold DB (public): predicted 3-D structure with per-residue confidence.

The model's pLDDT (0-100) sits in the B-factor column of the PDB file, so
every residue arrives with its own confidence: high (>90), confident (70-90),
low (50-70), very low (<50). Structures are `predicted` evidence by nature.
Files are ~100-500 KB and cached under data/cache/af (ignored by git).
"""

from __future__ import annotations

import json
import math
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
CACHE = Path("data/cache/af")
EVIDENCE = "AlphaFold DB (EMBL-EBI / DeepMind), predicted structure; pLDDT per residue"


@dataclass(slots=True)
class Structure:
    entry: str
    accession: str
    ca: list[tuple[float, float, float]] = field(default_factory=list)
    plddt: list[float] = field(default_factory=list)
    residues: list[str] = field(default_factory=list)  # one-letter codes
    version: str = ""

    @property
    def length(self) -> int:
        return len(self.ca)

    def mean_plddt(self) -> float:
        return sum(self.plddt) / len(self.plddt) if self.plddt else 0.0

    def confident_fraction(self, threshold: float = 70.0) -> float:
        return sum(1 for p in self.plddt if p >= threshold) / len(self.plddt) if self.plddt else 0.0

    def radius_of_gyration(self) -> float:
        if not self.ca:
            return 0.0
        cx = sum(p[0] for p in self.ca) / len(self.ca)
        cy = sum(p[1] for p in self.ca) / len(self.ca)
        cz = sum(p[2] for p in self.ca) / len(self.ca)
        return math.sqrt(
            sum((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2 for x, y, z in self.ca) / len(self.ca)
        )

    def to_dict(self, coords: bool = True) -> dict:
        d = {
            "entry": self.entry,
            "accession": self.accession,
            "version": self.version,
            "length": self.length,
            "mean_plddt": round(self.mean_plddt(), 1),
            "confident_fraction": round(self.confident_fraction(), 3),
            "radius_of_gyration_A": round(self.radius_of_gyration(), 1),
            "evidence": EVIDENCE,
            "confidence": round(self.mean_plddt() / 100, 2),
        }
        if coords:
            d["ca"] = [
                [round(x, 2), round(y, 2), round(z, 2), round(p, 1)]
                for (x, y, z), p in zip(self.ca, self.plddt, strict=True)
            ]
            d["residues"] = "".join(self.residues)
        return d


AA3 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def parse_pdb(text: str, entry: str = "", accession: str = "") -> Structure:
    s = Structure(entry, accession)
    for line in text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            s.ca.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            s.plddt.append(float(line[60:66]))
            s.residues.append(AA3.get(line[17:20], "X"))
    return s


def fetch_structure(accession: str, timeout: int = 120) -> Structure | None:
    CACHE.mkdir(parents=True, exist_ok=True)
    cached = next(CACHE.glob(f"AF-{accession}-F1-*.pdb"), None)
    if cached is None:
        req = urllib.request.Request(API.format(acc=accession), headers={"User-Agent": "GenomeOS/0.1"})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:  # noqa: S310
                meta = json.load(r)
        except Exception:  # noqa: BLE001
            return None
        if not meta:
            return None
        url = meta[0]["pdbUrl"]
        cached = CACHE / url.rsplit("/", 1)[-1]
        with urllib.request.urlopen(
            urllib.request.Request(url, headers={"User-Agent": "GenomeOS/0.1"}), timeout=timeout
        ) as r:  # noqa: S310
            cached.write_bytes(r.read())
    s = parse_pdb(cached.read_text(), cached.stem, accession)
    s.version = cached.stem.rsplit("-", 1)[-1]
    return s
