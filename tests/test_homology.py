"""Origin per gene and age per library from a streamed Compara dump (area J step 2)."""

import gzip

from genomeos.knowledge.homology import (
    LADDER,
    PROXIES,
    RANK,
    STRATA,
    deepest,
    distil,
    merge_genes,
    place,
    presence_matrix,
    stream,
    strip_species_sets,
    taxonomy_names,
)

HEADER = (
    "gene_stable_id\tprotein_stable_id\tspecies\tidentity\thomology_type\thomology_gene_stable_id\t"
    "homology_protein_stable_id\thomology_species\thomology_identity\tdn\tds\tgoc_score\twga_coverage\t"
    "is_high_confidence\thomology_id\n"
)


def row(gene, typ, partner, species):
    return (
        f"{gene}\tP\thomo_sapiens\t50\t{typ}\t{partner}\tPP\t{species}\t50\tNULL\tNULL\tNULL\tNULL\tNULL\t1\n"
    )


def test_ladder_is_total_and_deepest_picks_the_oldest():
    assert len(STRATA) == len(RANK) == len(LADDER)
    assert deepest(["Mammalia", "Eukaryota", None, "Homo"]) == "Eukaryota"
    assert deepest([None]) is None and deepest([]) is None


def test_stream_and_distil(tmp_path):
    dump = tmp_path / "compara.tsv.gz"
    with gzip.open(dump, "wt") as fh:
        fh.write(HEADER)
        # RPL3-like: orthologues to yeast, fly, fish and mouse, one paralogue
        fh.write(row("ENSG1", "ortholog_one2one", "Y1", "saccharomyces_cerevisiae"))
        fh.write(row("ENSG1", "ortholog_one2one", "F1", "drosophila_melanogaster"))
        fh.write(row("ENSG1", "ortholog_one2many", "Z1", "danio_rerio"))
        fh.write(row("ENSG1", "ortholog_one2one", "M1", "mus_musculus"))
        fh.write(row("ENSG1", "within_species_paralog", "ENSG2", "homo_sapiens"))
        # a mammal-only gene, with a species the ladder has not placed
        fh.write(row("ENSG2", "ortholog_one2one", "M2", "mus_musculus"))
        fh.write(row("ENSG2", "ortholog_one2one", "X2", "unplaced_species"))
        fh.write(row("ENSG2", "other_paralog", "ENSG1", "homo_sapiens"))
        # a gene with paralogues only: human-specific as far as this set sees
        fh.write(row("ENSG3", "within_species_paralog", "ENSG4", "homo_sapiens"))
        # a non-coding gene that must drop out
        fh.write(row("ENSG5", "ortholog_one2one", "M5", "mus_musculus"))
    strata = {
        "saccharomyces_cerevisiae": "Eukaryota",
        "drosophila_melanogaster": "Metazoa",
        "danio_rerio": "Euteleostomi",
        "mus_musculus": "Euarchontoglires",
    }
    genes, cost = stream(str(dump), strata)
    assert cost["rows"] == 10 and cost["unplaced_species"] == {"unplaced_species": 1}
    assert genes["ENSG1"]["one2one"] == 3 and len(genes["ENSG1"]["species"]) == 4
    symbols = {
        "ENSG1": ("RPL3", "protein_coding"),
        "ENSG2": ("MAMM1", "protein_coding"),
        "ENSG3": ("NEW1", "protein_coding"),
        "ENSG4": ("NEW2", "protein_coding"),
        "ENSG5": ("LNC1", "lncRNA"),
    }
    members = {"core.translation": ["RPL3", "MAMM1", "MISSING"], "systems.immune": ["NEW1"]}
    d = distil(genes, symbols, members)
    g = d["genes"]
    assert g["RPL3"]["origin"] == "Eukaryota" and g["RPL3"]["paralogues"] == ["MAMM1"]
    assert g["RPL3"]["ladder"] == "Eukaryote core"
    assert g["MAMM1"]["origin"] == "Euarchontoglires" and g["MAMM1"]["paralogue_types"] == {
        "other_paralog": 1
    }
    assert g["NEW1"]["origin"] == "Homo" and g["NEW1"]["species"] == 0
    assert "LNC1" not in g and d["coding_genes"] == 3
    assert d["by_stratum"] == {"Eukaryota": 1, "Euarchontoglires": 1, "Homo": 1}
    lib = d["libraries"]["core.translation"]
    assert lib["members_placed"] == 2 and lib["members"] == 3 and lib["deepest"] == "Eukaryota"
    assert lib["median_origin"] == "Eukaryota" and lib["share_animal_or_older"] == 0.5
    assert lib["share_amniote_or_younger"] == 0.5 and lib["distribution"] == {
        "Eukaryota": 1,
        "Euarchontoglires": 1,
    }
    assert d["libraries"]["systems.immune"]["median_origin"] == "Homo"


def test_place_uses_proxies_for_the_nodes_ensembl_omits():
    assert set(PROXIES.values()) <= set(STRATA)
    chicken = {
        "Gallus",
        "Aves",
        "Archosauria",
        "Euteleostomi",
        "Vertebrata",
        "Chordata",
        "Metazoa",
        "Eukaryota",
    }
    assert place(chicken) == "Amniota"
    assert place({"Xenopus", "Amphibia", "Euteleostomi", "Vertebrata", "Eukaryota"}) == "Tetrapoda"
    assert place({"Bos", "Laurasiatheria", "Eutheria", "Mammalia", "Eukaryota"}) == "Boreoeutheria"
    assert place({"Monodelphis", "Marsupialia", "Mammalia", "Eukaryota"}) == "Theria"
    assert place({"Pan", "Homininae", "Hominidae", "Primates", "Eukaryota"}) == "Homininae"
    assert place({"Saccharomyces", "Fungi", "Eukaryota"}) == "Opisthokonta"
    assert place({"Eptatretus", "Cyclostomata", "Vertebrata", "Chordata", "Eukaryota"}) == "Vertebrata"
    # human's own list must not place a species at Homo; an unknown list places nothing
    assert place({"Homo", "Hominidae"}) == "Hominidae" and place({"Nothing"}) is None


def test_taxonomy_names_strip_assembly_and_strain_suffixes():
    assert taxonomy_names("bos_taurus") == ["bos_taurus"]
    assert taxonomy_names("bos_taurus_gca963921495v1") == ["bos_taurus_gca963921495v1", "bos_taurus"]
    assert taxonomy_names("mus_musculus_129s1svimj") == ["mus_musculus_129s1svimj", "mus_musculus"]
    assert taxonomy_names("canis_lupus_familiaris") == ["canis_lupus_familiaris", "canis_lupus"]
    bug = "escherichia_coli_str_k_12_substr_mg1655_gca_000005845"
    assert taxonomy_names(bug) == [bug, "escherichia_coli_str_k_12_substr_mg1655", "escherichia_coli"]


def test_life_stratum_and_merge_of_two_streams():
    from collections import Counter

    assert STRATA[0] == "Life" and LADDER["Life"] == "Life core"
    assert (
        place({"Escherichia", "Bacteria"}) == "Life"
        and place({"Arabidopsis", "Viridiplantae", "Eukaryota"}) == "Eukaryota"
    )
    a = {
        "G1": {
            "rank": RANK["Chordata"],
            "species": {"ciona"},
            "one2one": 1,
            "paralogues": {"G2"},
            "paralogue_types": Counter({"within_species_paralog": 1}),
        }
    }
    b = {
        "G1": {
            "rank": RANK["Life"],
            "species": {"e_coli"},
            "one2one": 1,
            "paralogues": set(),
            "paralogue_types": Counter(),
        },
        "G3": {
            "rank": None,
            "species": set(),
            "one2one": 0,
            "paralogues": set(),
            "paralogue_types": Counter(),
        },
    }
    m = merge_genes(a, b)
    assert (
        m["G1"]["rank"] == RANK["Life"]
        and m["G1"]["species"] == {"ciona", "e_coli"}
        and m["G1"]["one2one"] == 2
    )
    assert m["G1"]["paralogues"] == {"G2"} and "G3" in m


def test_species_fallback_places_renamed_species():
    from genomeos.knowledge.homology import SPECIES_FALLBACK

    assert set(SPECIES_FALLBACK.values()) <= set(STRATA)
    assert SPECIES_FALLBACK["physeter_catodon"] == "Boreoeutheria"


def test_presence_matrix_orders_species_by_stratum_and_packs_bits():
    per_gene = {
        "A": {"_species_set": {"mus_musculus", "danio_rerio", "saccharomyces_cerevisiae"}},
        "B": {"_species_set": {"mus_musculus"}},
        "C": {"_species_set": set()},
    }
    strata = {
        "mus_musculus": "Euarchontoglires",
        "danio_rerio": "Euteleostomi",
        "saccharomyces_cerevisiae": "Opisthokonta",
    }
    m = presence_matrix(per_gene, strata)
    assert m["species"] == ["saccharomyces_cerevisiae", "danio_rerio", "mus_musculus"]  # deep to shallow
    assert m["genes"]["A"] == "7" and m["genes"]["B"] == "1" and m["genes"]["C"] == "0"
    strip_species_sets(per_gene)
    assert "_species_set" not in per_gene["A"]
