# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The packaged healthy-tissue atlas, and the two tissues it recovered.

The distillation needs the network (`scripts/normal_tissue.py`); these tests
read the table it ships, so the surface universe, the tissue coverage and the
offline safety answer are all checked in CI without a request.
"""

from __future__ import annotations

import pytest

from genomeos.results import load_result
from genomeos.therapeutics.atlas import PACKAGED, load, row, universe
from genomeos.therapeutics.expression import normal_profile, normal_tissue_safety
from genomeos.therapeutics.providers import CRITICAL_TISSUES, DiskCache, HpaExpressionProvider, hpa_column

pytestmark = pytest.mark.skipif(not PACKAGED.exists(), reason="run scripts/normal_tissue.py")

#: no per-gene cache and no network: only the packaged table can answer
OFFLINE = HpaExpressionProvider(net=False, cache=DiskCache("hpa_absent_in_tests"))


def test_the_surface_universe_is_shipped_and_large():
    genes = universe()
    assert len(genes) > 5_000
    assert {"ERBB2", "EGFR", "MSLN", "FOLH1"} <= set(genes)


def test_the_antigens_of_the_approved_cell_therapies_are_present():
    """CD19 and BCMA are the point of the CD-marker query: neither is mutated."""
    assert row("CD19") and row("TNFRSF17")
    assert row("CD19").get("cd_marker") is True


def test_every_weighed_tissue_is_actually_read():
    """The defect this table found: two of the twenty were silently never read.

    The Atlas returns the `t_RNA_skin_1` and `t_RNA_stomach_1` fields under the
    titles "skin 1" and "stomach 1"; the lookup asked for "skin" and "stomach"
    and got nothing, so skin — where an EGFR antibody's classic toxicity shows
    — never reached the safety score.
    """
    assert hpa_column("t_RNA_skin_1") == "Tissue RNA - skin 1 [nTPM]"
    assert hpa_column("t_RNA_heart_muscle") == "Tissue RNA - heart muscle [nTPM]"
    egfr = row("EGFR")
    for field in CRITICAL_TISSUES:
        assert hpa_column(field) in egfr, field
    profile = normal_profile("EGFR", OFFLINE.normal_expression("EGFR"))
    measured = [t for t in profile.tissues if t.value is not None]
    assert len(measured) == len(CRITICAL_TISSUES)
    skin = next(t for t in profile.tissues if t.tissue == "skin")
    assert skin.value and skin.value > 10


def test_a_surface_gene_gets_its_safety_answered_with_no_network():
    answer = OFFLINE.normal_expression("ERBB2")
    assert answer.available
    profile = normal_profile("ERBB2", answer)
    score, basis = normal_tissue_safety(profile)
    assert score is not None and "Human Protein Atlas" in basis
    assert profile.on_target_off_tumour_risk in ("low", "moderate", "high")


def test_an_intracellular_gene_is_absent_rather_than_answered_as_safe():
    """The table covers the membrane classes; outside them it says nothing."""
    assert row("BRAF") is None
    answer = OFFLINE.normal_expression("BRAF")
    assert not answer.available
    assert "no cached HPA record" in answer.reason


def test_the_table_states_its_own_limits():
    meta = load()["meta"]
    assert meta["genes"] > 5_000
    assert len(meta["tissues"]) == len(CRITICAL_TISSUES)
    assert any("not the body" in x for x in meta["limits"])
    assert any("RNA is not protein" in x for x in meta["limits"])
    result = load_result("normal_tissue_atlas")
    if result:
        assert result["genes"] == meta["genes"]
        assert result["protein_class_flags"]["cd_marker"] > 300
