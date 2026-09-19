from pathlib import Path

from genomeos.lang.parser import parse
from genomeos.molecules.reactome import PathwayModel, PathwayReaction, Species

FIXTURE = Path("data/models/reactome_R-HSA-69541.sbml")


def _toy() -> PathwayModel:
    sp = {
        "A": Species("A", "A", "cytosol", {"P1"}, "protein"),
        "B": Species("B", "B", "cytosol", {"P2"}, "protein"),
        "E": Species("E", "enzyme", "cytosol", {"P3"}, "protein"),
        "AB": Species("AB", "A:B complex", "cytosol", {"P1", "P2"}, "complex"),
        "C": Species("C", "C", "nucleus", {"P4"}, "protein"),
        "D": Species("D", "D", "nucleus", set(), "chemical"),
        "D2": Species("D2", "D2", "nucleus", set(), "chemical"),
    }
    rx = [
        PathwayReaction("r1", "A binds B", ["A", "B"], ["AB"], [], []),
        PathwayReaction("r2", "AB makes C", ["AB"], ["C"], ["E"], []),
        PathwayReaction("r3", "D on its own", ["D"], ["D2"], [], []),
    ]
    return PathwayModel("R-HSA-0", "toy", sp, rx, "97")


def test_reachability_and_knockout_on_a_toy_pathway():
    m = _toy()
    assert m.sources() == {"A", "B", "E", "D"}
    base = m.reach()
    assert set(base["fired"]) == {"r1", "r2", "r3"} and "C" in base["present"]
    k = m.knockout("P2")
    assert k["entities_containing"] == ["A:B complex", "B"]
    assert [r["name"] for r in k["reactions_lost"]] == ["A binds B", "AB makes C"]
    assert k["products_unreachable"] == ["A:B complex", "C"]
    assert k["fraction_lost"] == round(2 / 3, 3)
    # losing the catalyst blocks only the catalysed step
    k2 = m.knockout("P3")
    assert [r["name"] for r in k2["reactions_lost"]] == ["AB makes C"]
    assert m.knockout("P9")["reactions_lost"] == []


def test_reactome_export_parses_with_uniprot_parts():
    m = PathwayModel.from_sbml(FIXTURE)
    assert m.id == "R-HSA-69541" and m.name.startswith("Stabilization of p53") and m.version == "97"
    assert len(m.reactions) == 15 and len(m.species) == 28
    assert m.species_with("P04637")  # TP53 appears in several entities
    binds = next(r for r in m.reactions if r.name == "MDM2 binds TP53")
    assert len(binds.inputs) == 2 and len(binds.outputs) == 1 and binds.inhibitors
    k = m.knockout("P04637")
    assert "MDM2 binds TP53" in {r["name"] for r in k["reactions_lost"]}
    assert 0 < k["fraction_lost"] < 1


def test_biolang_protein_block_carries_compiled_properties():
    src = """module t
protein TP53 {
  accession: P04637; sequence: MEEPQ; isoforms: P04637-1, P04637-2; domains: p53_DNA-bd, p53_TAD2;
  pathways: R-HSA-69541; interactions: MDM2, CREBBP; structures: 1TUP, AF-P04637;
  evidence: curated "UniProt"; confidence: 0.9
}
"""
    m = parse(src)
    p = m.entities["TP53"]
    assert p.accession == "P04637" and p.isoforms == ["P04637-1", "P04637-2"]
    assert [d["name"] for d in p.domains] == ["p53_DNA-bd", "p53_TAD2"]
    assert p.pathways == ["R-HSA-69541"] and p.interactions == ["MDM2", "CREBBP"]
    assert p.structures == [{"source": "PDB", "id": "1TUP"}, {"source": "AlphaFold", "id": "AF-P04637"}]
