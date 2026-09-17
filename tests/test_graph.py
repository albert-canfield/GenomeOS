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


def test_orthologue_counts_are_per_clade_and_absent_where_unread(tmp_path):
    """One total over 355 species cannot tell a vertebrate-wide gene from a primate expansion.

    TP53 and APP sit in all 65 Euteleostomi; OR5H1, an olfactory receptor, is Boreoeutheria-restricted
    and absent from the deeper clades. A clade a gene has no reading in is left out rather than set to
    0, because "not read there" and "missing there" are different statements.
    """
    import json

    from genomeos.molecules import graph as gmod

    (tmp_path / "origin_presence_genome_wide.json").write_text(
        json.dumps(
            {
                "species": ["a", "b", "c", "d"],
                "strata": {"a": "Deep", "b": "Deep", "c": "Shallow", "d": "Shallow"},
                # bit i from the LEFT: 1100 means present in a and b only
                "genes": {"DEEPGENE": "c", "SHALLOWGENE": "3", "NOWHERE": "0"},
            }
        )
    )
    out = gmod._orthologues_by_clade(tmp_path)

    assert out["DEEPGENE"] == {"Deep": 2}
    assert out["SHALLOWGENE"] == {"Shallow": 2}
    assert "NOWHERE" not in out  # no species at all is no counts, not a row of zeros


def test_a_missing_presence_result_is_a_no_op(tmp_path):
    from genomeos.molecules import graph as gmod

    assert gmod._orthologues_by_clade(tmp_path) == {}
