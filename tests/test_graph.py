import json

from genomeos.molecules.graph import build, summarise


def _defn(gene, acc, pathways=(), partners=(), domains=(), tissues=None):
    return {
        "id": f"UniProt:{acc}",
        "gene": gene,
        "sections": {
            "identity": {"items": {"accession": acc, "name": gene + " protein", "length": 100}},
            "pathways": {"items": [{"id": p, "name": p + " pathway"} for p in pathways]},
            "domains": {"items": {"interpro": [{"id": d, "name": d} for d in domains], "features": []}},
            "interactions": {
                "items": [{"partner": p, "score": s, "physical_evidence": phys} for p, s, phys in partners]
            },
            "expression": {"items": {"tissue_ntpm": tissues or {}}},
        },
    }


def test_graph_from_cached_definitions(tmp_path):
    (tmp_path / "A.json").write_text(
        json.dumps(
            _defn(
                "A",
                "P1",
                ["R-1"],
                [("B", 0.9, True), ("C", 0.75, False), ("D", 0.5, True)],
                ["IPR1"],
                {"liver": 12.0},
            )
        )
    )
    (tmp_path / "B.json").write_text(json.dumps(_defn("B", "P2", ["R-1", "R-2"], [("A", 0.9, True)], [])))
    g = build(tmp_path)
    kinds = {n["kind"] for n in g.nodes.values()}
    assert kinds == {"protein", "pathway", "domain", "tissue"}
    assert g.nodes["C"]["compiled"] is False and g.nodes["B"]["compiled"] is True
    # the A–B pair appears once even though both definitions list it; D is below the score cut
    assoc = [e for e in g.edges if e["rel"] == "associates"]
    assert sorted((e["a"], e["b"]) for e in assoc) == [("A", "B"), ("A", "C")]
    assert g.degree("A", "member_of") == 1 and g.degree("B", "member_of") == 2
    n = g.neighbourhood("A")
    assert {x["id"] for x in n["nodes"]} == {"A", "B", "C", "R-1", "IPR1", "tissue:liver"}
    s = summarise(g)
    assert s["compiled_proteins"] == 2 and s["physical_associations"] == 1
    assert s["hubs"][0]["gene"] == "A" and s["largest_component"] == 2
    assert s["biggest_pathways"][0] == ("R-1 pathway", 2)
