# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The cancer.* library layer: which genes break, how, and in which cancers.

Membership is computed from the two committed distillations rather than
written down, so these tests are about the computation and about the three
things the layer refuses to say.
"""

from __future__ import annotations

import pytest

from genomeos.cancer.libraries import (
    CANCER_LIBRARIES,
    Knowledge,
    breaks_in,
    cancer_libraries_of,
    library_members,
    register,
    summary,
)

K = Knowledge.load()
pytestmark = pytest.mark.skipif(
    not K.available, reason="run genomeos cancer distil and genomeos cancer alterations"
)


def test_every_library_is_computed_and_names_its_source():
    for lib in CANCER_LIBRARIES.values():
        assert lib.id.startswith("cancer.")
        assert lib.purpose and lib.note and lib.source.startswith("cBioPortal ")
        assert library_members(lib, K), lib.id


def test_the_layer_separates_the_ways_a_gene_breaks():
    """Each library holds the genes that break that way, and only those."""
    amplified = library_members(CANCER_LIBRARIES["cancer.amplified"], K)
    deleted = library_members(CANCER_LIBRARIES["cancer.deleted"], K)
    rearranged = library_members(CANCER_LIBRARIES["cancer.rearranged"], K)
    assert {"ERBB2", "MYC", "CCND1", "EGFR", "CDK4"} <= set(amplified)
    assert {"CDKN2A", "PTEN", "RB1"} <= set(deleted)
    assert {"ALK", "RET", "ROS1"} <= set(rearranged)
    assert "CDKN2A" not in amplified and "ERBB2" not in deleted


def test_a_gene_the_mutation_table_calls_a_passenger_is_found_by_the_layer():
    """MYC and CCND1 are in cancer.amplified and not in cancer.mutated."""
    for gene in ("MYC", "CCND1"):
        libs = cancer_libraries_of(gene, K)
        assert "cancer.amplified" in libs
        assert "cancer.mutated" not in libs


def test_membership_carries_the_number_that_put_the_gene_there():
    members = library_members(CANCER_LIBRARIES["cancer.deleted"], K)
    assert members["CDKN2A"] > 0.05
    assert list(members) == sorted(members, key=lambda g: -members[g]), "ordered by frequency"
    assert all(f >= CANCER_LIBRARIES["cancer.deleted"].threshold for f in members.values())


def test_breaks_in_gives_the_cancer_type_where_the_break_actually_lives():
    b = breaks_in("ERBB2", K)
    amp = next(w for w in b["breaks"] if w["kind"] == "amplification")
    assert amp["frequency"] < 0.05
    assert amp["by_cancer_type"]["Breast Cancer"] > 0.10
    assert "Esophagogastric" in amp["evidence"]["claim"]
    assert amp["evidence"]["kind"] == "experimental"


def test_a_hotspot_never_borrows_the_mutation_distribution():
    """The distillation counts recurrent changes study-wide, not per cancer type.

    Reporting the gene's mutation frequencies beside a hotspot frequency would
    read as the hotspot's distribution, which nothing measured.
    """
    b = breaks_in("KRAS", K)
    hotspot = next(w for w in b["breaks"] if w["kind"] == "hotspot")
    assert hotspot["by_cancer_type"] == {}
    assert "no per-type breakdown of a hotspot exists" in hotspot["meaning"]
    assert hotspot["recurrent_changes"][0][0] in ("G12D", "G12V", "G12C")
    mutation = next(w for w in b["breaks"] if w["kind"] == "mutation")
    assert mutation["by_cancer_type"]["Pancreatic Cancer"] > 0.5


def test_a_gene_off_the_panel_is_not_called_healthy():
    """CD19 is a real target and this study never looked at it."""
    b = breaks_in("CD19", K)
    assert b["on_panel"] is False
    assert b["libraries"] == []
    assert b["breaks"] == []
    assert "silence is not evidence" in b["note"]
    assert set(b["not_called"]) == {lib.kind for lib in CANCER_LIBRARIES.values()}


def test_the_layer_states_its_own_limits():
    s = summary(K)
    assert s["layer"] == "cancer" and s["samples"] > 10_000
    assert len(s["limits"]) == 3
    assert any("silence is not evidence" in x for x in s["limits"])


def test_registering_into_the_biolib_catalogue_is_explicit_and_additive():
    """It is another area's catalogue, so nothing joins it as an import side effect."""
    from genomeos.lib import LAYERS, LIBRARIES

    assert not any(k.startswith("cancer.") for k in LIBRARIES), (
        "importing genomeos.cancer must not mutate the shared catalogue"
    )
    sandbox: dict = {}
    layers = register(sandbox, LAYERS)
    assert "cancer" in layers and "cancer" not in LAYERS
    assert set(sandbox) == set(CANCER_LIBRARIES)
    lib = sandbox["cancer.deleted"]
    assert lib.layer == "cancer" and "CDKN2A" in lib.genes and lib.source and lib.purpose
