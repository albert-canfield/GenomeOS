"""Compiler normalisation on a fixture entry (no network)."""

from genomeos.molecules.compiler import ProteinState, coverage, normalise_uniprot, states_from_definition

FIXTURE = {
    "primaryAccession": "P00001",
    "proteinDescription": {"recommendedName": {"fullName": {"value": "Test protein"}}},
    "genes": [{"geneName": {"value": "TST1"}}],
    "sequence": {"length": 6, "value": "MKDEFG"},
    "proteinExistence": "1: Evidence at protein level",
    "keywords": [{"name": "Nucleus"}],
    "comments": [
        {"commentType": "FUNCTION", "texts": [{"value": "Does a thing."}]},
        {"commentType": "SUBCELLULAR LOCATION", "subcellularLocations": [{"location": {"value": "Nucleus"}}]},
        {
            "commentType": "ALTERNATIVE PRODUCTS",
            "isoforms": [
                {"isoformIds": ["P00001-1"], "name": {"value": "1"}, "isoformSequenceStatus": "Displayed"},
                {"isoformIds": ["P00001-2"], "name": {"value": "2"}, "isoformSequenceStatus": "Described"},
            ],
        },
        {
            "commentType": "DISEASE",
            "disease": {
                "diseaseId": "Test syndrome",
                "acronym": "TS",
                "diseaseCrossReference": {"id": "123456"},
                "description": "d",
            },
        },
    ],
    "features": [
        {"type": "Domain", "description": "Kinase", "location": {"start": {"value": 1}, "end": {"value": 4}}},
        {
            "type": "Modified residue",
            "description": "Phosphoserine",
            "location": {"start": {"value": 3}, "end": {"value": 3}},
        },
        {"type": "Signal", "description": "", "location": {"start": {"value": 1}, "end": {"value": 2}}},
    ],
    "uniProtKBCrossReferences": [
        {
            "database": "PDB",
            "id": "1ABC",
            "properties": [
                {"key": "Method", "value": "X-ray"},
                {"key": "Resolution", "value": "2.10 A"},
                {"key": "Chains", "value": "A=1-6"},
            ],
        },
        {
            "database": "PDB",
            "id": "2DEF",
            "properties": [
                {"key": "Method", "value": "NMR"},
                {"key": "Resolution", "value": "-"},
                {"key": "Chains", "value": "A=1-6"},
            ],
        },
        {
            "database": "InterPro",
            "id": "IPR000001",
            "properties": [{"key": "EntryName", "value": "Kinase_dom"}],
        },
        {
            "database": "Reactome",
            "id": "R-HSA-1",
            "properties": [{"key": "PathwayName", "value": "Signalling"}],
        },
        {"database": "AlphaFoldDB", "id": "P00001", "properties": []},
        {
            "database": "Ensembl",
            "id": "ENST1",
            "properties": [{"key": "ProteinId", "value": "ENSP1"}],
            "isoformId": "P00001-1",
        },
    ],
}


def test_normalise_uniprot_splits_into_evidence_tagged_sections():
    s = normalise_uniprot(FIXTURE)
    assert s["identity"]["items"]["accession"] == "P00001" and s["identity"]["evidence"] == "curated"
    assert s["function"]["items"]["summary"] == ["Does a thing."]
    assert [i["id"] for i in s["isoforms"]["items"]] == ["P00001-1", "P00001-2"]
    assert s["isoforms"]["ensembl_products"][0]["protein"] == "ENSP1"
    exp = s["structures_experimental"]
    assert exp["evidence"] == "experimental" and exp["count"] == 2
    assert exp["items"][0]["method"] == "X_RAY" and exp["items"][0]["resolution_A"] == 2.1
    assert exp["items"][1]["method"] == "NMR" and exp["items"][1]["resolution_A"] is None
    assert s["domains"]["items"]["interpro"] == [{"id": "IPR000001", "name": "Kinase_dom"}]
    assert [f["type"] for f in s["domains"]["items"]["features"]] == ["Domain"]
    assert [f["type"] for f in s["modifications"]["items"]] == ["Modified residue"]
    assert [f["type"] for f in s["processing"]["items"]] == ["Signal"]
    assert s["pathways"]["items"] == [{"id": "R-HSA-1", "name": "Signalling"}]
    assert s["diseases"]["items"][0]["mim"] == "123456"
    assert s["alphafold_ids"] == ["P00001"]
    # predicted and experimental never share a section
    assert "structures_predicted" not in s


def test_coverage_and_states():
    s = normalise_uniprot(FIXTURE)
    s.pop("alphafold_ids")
    defn = {"id": "UniProt:P00001", "gene": "TST1", "sections": s}
    cov = coverage(defn)
    assert cov["sequence"] and cov["structure_experimental"] and cov["pathways"] and cov["disease"]
    assert not cov["interactions"] and not cov["expression"] and not cov["structure_predicted"]
    defn["sections"]["expression"] = {
        "items": {"subcellular_main": ["Nucleoplasm"], "tissue_ntpm": {"liver": 40.0, "brain": 5.5}}
    }
    st = states_from_definition(defn)
    assert len(st) == 2 and isinstance(st[0], ProteinState)
    assert st[0].protein == "UniProt:P00001" and st[0].tissue == "liver" and st[0].level == 40.0
    assert st[0].localisation == "Nucleoplasm" and st[0].evidence.startswith("experimental")


def test_cached_failures_are_refetched(tmp_path, monkeypatch):
    import json

    from genomeos.molecules import compiler

    bad = {
        "id": "UniProt:P00001",
        "gene": "TST1",
        "sections": {
            "identity": {
                "items": {"accession": "P00001", "sequence": "MK", "name": "t"},
                "evidence": "curated",
                "source": "x",
                "confidence": 0.9,
            },
            "genomic_origin": {
                "items": None,
                "evidence": "none",
                "confidence": 0,
                "source": "x",
                "error": "HTTP 500",
            },
        },
    }
    (tmp_path / "TST1.json").write_text(json.dumps(bad))
    calls = []
    monkeypatch.setattr(
        compiler,
        "ensembl_gene",
        lambda sym: (
            calls.append(sym)
            or {
                "gene_id": "ENSG1",
                "biotype": "protein_coding",
                "locus": "chr1:1-2(+)",
                "description": "",
                "transcripts": [],
                "protein_products": 0,
            }
        ),
    )  # noqa: E501
    monkeypatch.setattr(
        compiler, "uniprot_raw", lambda sym: (_ for _ in ()).throw(AssertionError("must not refetch"))
    )
    monkeypatch.setattr(compiler, "hpa_entry", lambda gid: {"tissue_specificity": "x"})
    d = compiler.compile_protein("TST1", cache_dir=tmp_path)
    assert calls == ["TST1"]
    assert d["sections"]["genomic_origin"]["evidence"] == "curated"
    assert d["sections"]["identity"]["items"]["accession"] == "P00001"
    assert d["coverage"]["genomic_origin"]
    # second load: nothing failed, nothing refetched
    compiler.compile_protein("TST1", cache_dir=tmp_path)
    assert calls == ["TST1"]


def test_origin_from_local_models_skips_ensembl(tmp_path, monkeypatch):
    from genomeos.molecules import compiler

    monkeypatch.setattr(compiler, "uniprot_raw", lambda sym: FIXTURE)
    monkeypatch.setattr(
        compiler, "ensembl_gene", lambda sym: (_ for _ in ()).throw(AssertionError("no Ensembl"))
    )
    monkeypatch.setattr(compiler, "string_network", lambda sym: [])
    monkeypatch.setattr(compiler, "hpa_entry", lambda gid: {"tissue_specificity": "x", "gene_id_seen": gid})
    origin = {
        "gene_id": "ENSG1",
        "biotype": "protein_coding",
        "locus": "chr1:1-2(+)",
        "description": None,
        "transcripts": [],
        "protein_products": 1,
    }
    d = compiler.compile_protein("TST1", cache_dir=tmp_path, origin=origin)
    go = d["sections"]["genomic_origin"]
    assert go["evidence"] == "curated" and "GENCODE" in go["source"] and go["items"]["gene_id"] == "ENSG1"
    assert d["sections"]["expression"]["items"]["gene_id_seen"] == "ENSG1"
    assert d["coverage"]["genomic_origin"] and d["coverage"]["expression"]
