"""Reactome pathways as executable graphs.

Reactome exports every pathway as SBML: species (proteins, complexes, small
molecules, each annotated with the UniProt parts it contains) and reactions
(inputs, outputs, catalysts, regulators). The export has no rate laws, so
it cannot be integrated like a BioModels file; what it *can* do exactly is
reachability: given which molecules are present, which reactions can fire
and which products can be made. That is enough for the first executable
question of the protein layer: what stops working when a protein is lost.

Evidence: curated (Reactome, version recorded from the file); the knockout
logic itself is inferred (a reaction needs all its inputs and catalysts).
"""

from __future__ import annotations

import re
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

CACHE = Path("data/knowledge/pathways")
EVIDENCE = "curated: Reactome SBML export (ContentService/exporter/event)"
SBO_REACTANT, SBO_PRODUCT, SBO_CATALYST, SBO_INHIBITOR = (
    "SBO:0000010",
    "SBO:0000011",
    "SBO:0000013",
    "SBO:0000020",
)


def fetch_pathway(pathway_id: str, cache_dir: Path = CACHE, timeout: int = 120) -> Path:
    """Download a pathway's SBML once; later calls read the cached file."""
    if not re.fullmatch(r"R-HSA-\d+", pathway_id):
        raise ValueError(f"not a human Reactome pathway id: {pathway_id}")
    cache_dir.mkdir(parents=True, exist_ok=True)
    out = cache_dir / f"{pathway_id}.sbml"
    if out.exists():
        return out
    url = f"https://reactome.org/ContentService/exporter/event/{pathway_id}.sbml"
    req = urllib.request.Request(url, headers={"User-Agent": "GenomeOS/0.1"})
    with urllib.request.urlopen(req, timeout=timeout) as r:  # noqa: S310
        out.write_bytes(r.read())
    return out


@dataclass(slots=True)
class Species:
    id: str
    name: str
    compartment: str
    parts: set[str] = field(default_factory=set)  # UniProt accessions this entity contains
    kind: str = "entity"  # protein | complex | chemical | set | entity


@dataclass(slots=True)
class PathwayReaction:
    id: str
    name: str
    inputs: list[str]
    outputs: list[str]
    catalysts: list[str]
    inhibitors: list[str]
    reactome_id: str = ""


@dataclass(slots=True)
class PathwayModel:
    id: str
    name: str
    species: dict[str, Species]
    reactions: list[PathwayReaction]
    version: str = ""

    @classmethod
    def from_sbml(cls, path: str | Path) -> PathwayModel:
        root = ET.parse(path).getroot()
        ns = root.tag.split("}")[0] + "}"
        model = root.find(f"{ns}model")
        if model is None:
            raise ValueError("no <model> element")
        notes = root.find(f"{ns}notes")
        version = ""
        if notes is not None:
            m = re.search(r"Reactome version (\d+)", "".join(notes.itertext()))
            version = m.group(1) if m else ""

        def items(list_name: str) -> list[ET.Element]:
            lst = model.find(f"{ns}{list_name}")
            return list(lst) if lst is not None else []

        comps = {c.get("id"): c.get("name", c.get("id")) for c in items("listOfCompartments")}
        species: dict[str, Species] = {}
        for sp in items("listOfSpecies"):
            sid = sp.get("id")
            text = ET.tostring(sp, encoding="unicode")
            parts = set(re.findall(r"identifiers\.org/uniprot[:/]([A-Z0-9]+)", text))
            sbo = sp.get("sboTerm", "")
            kind = {"SBO:0000297": "protein", "SBO:0000253": "complex", "SBO:0000247": "chemical"}.get(
                sbo, "entity"
            )
            if "Derived from a Reactome DefinedSet" in text or "CandidateSet" in text:
                kind = "set"
            species[sid] = Species(
                sid, sp.get("name", sid), comps.get(sp.get("compartment"), ""), parts, kind
            )

        def refs(rx: ET.Element, list_name: str, sbo: str | None = None) -> list[str]:
            lst = rx.find(f"{ns}{list_name}")
            if lst is None:
                return []
            return [r.get("species") for r in lst if sbo is None or r.get("sboTerm") == sbo]

        reactions = []
        for rx in items("listOfReactions"):
            text = ET.tostring(rx, encoding="unicode")
            rid = re.search(r"identifiers\.org/reactome[:/](R-HSA-\d+)", text)
            reactions.append(
                PathwayReaction(
                    rx.get("id"),
                    rx.get("name", rx.get("id")),
                    refs(rx, "listOfReactants"),
                    refs(rx, "listOfProducts"),
                    refs(rx, "listOfModifiers", SBO_CATALYST),
                    refs(rx, "listOfModifiers", SBO_INHIBITOR),
                    rid.group(1) if rid else "",
                )
            )
        pid = re.search(
            r"identifiers\.org/reactome[:/](R-HSA-\d+)", ET.tostring(model, encoding="unicode")[:5000]
        )
        return cls(
            pid.group(1) if pid else model.get("id", ""), model.get("name", ""), species, reactions, version
        )

    # ---- execution: reachability ----------------------------------------------------
    def sources(self) -> set[str]:
        produced = {s for r in self.reactions for s in r.outputs}
        return set(self.species) - produced

    def species_with(self, accession: str) -> set[str]:
        return {sid for sid, sp in self.species.items() if accession in sp.parts}

    def reach(self, present: set[str] | None = None, absent: set[str] | None = None) -> dict[str, Any]:
        """Fixpoint: a reaction fires when every input and catalyst is present; its outputs
        become present. Species in `absent` can never be present (a knockout)."""
        absent = absent or set()
        have = set(present if present is not None else self.sources()) - absent
        fired: list[str] = []
        changed = True
        while changed:
            changed = False
            for r in self.reactions:
                if r.id in fired:
                    continue
                needs = set(r.inputs) | set(r.catalysts)
                if needs and needs <= have:
                    fired.append(r.id)
                    new = set(r.outputs) - absent - have
                    if new:
                        have |= new
                        changed = True
                    else:
                        changed = changed or False
        return {
            "present": have,
            "fired": fired,
            "blocked": [r.id for r in self.reactions if r.id not in fired],
        }

    def knockout(self, accession: str) -> dict[str, Any]:
        """What is lost when every entity containing this protein is absent."""
        gone = self.species_with(accession)
        base = self.reach()
        ko = self.reach(absent=gone)
        lost_rx = [r for r in self.reactions if r.id in base["fired"] and r.id not in ko["fired"]]
        lost_sp = sorted((base["present"] - ko["present"]) - self.sources())
        return {
            "pathway": self.id,
            "name": self.name,
            "knockout": accession,
            "entities_containing": sorted(self.species[s].name for s in gone),
            "reactions_total": len(self.reactions),
            "reactions_reachable_baseline": len(base["fired"]),
            "reactions_lost": [{"id": r.reactome_id or r.id, "name": r.name} for r in lost_rx],
            "products_unreachable": [self.species[s].name for s in lost_sp],
            "fraction_lost": round(len(lost_rx) / len(base["fired"]), 3) if base["fired"] else 0.0,
            "evidence": EVIDENCE + (f" v{self.version}" if self.version else ""),
            "logic": "inferred: a reaction needs all inputs and catalysts; inhibitors reported, not applied",
            "confidence": 0.6,
        }

    def summary(self) -> dict[str, Any]:
        kinds: dict[str, int] = {}
        for s in self.species.values():
            kinds[s.kind] = kinds.get(s.kind, 0) + 1
        proteins = {p for s in self.species.values() for p in s.parts}
        base = self.reach()
        return {
            "pathway": self.id,
            "name": self.name,
            "version": self.version,
            "species": len(self.species),
            "species_kinds": kinds,
            "proteins": len(proteins),
            "reactions": len(self.reactions),
            "reactions_reachable_from_sources": len(base["fired"]),
            "catalysed": sum(1 for r in self.reactions if r.catalysts),
            "inhibited": sum(1 for r in self.reactions if r.inhibitors),
            "evidence": EVIDENCE,
        }
