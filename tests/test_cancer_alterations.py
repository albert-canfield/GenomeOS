# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Copy-number and structural-variant knowledge, and the alterations it grades.

The distillation itself needs the network (`genomeos cancer alterations`);
these tests read its committed result, so a change that loses the events or
mis-grades them fails in CI within a second.
"""

from __future__ import annotations

import pytest

from genomeos.results import load_result

K = load_result("cancer_alterations_msk_impact_2017")
MUT = load_result("cancer_msk_impact_2017")
pytestmark = pytest.mark.skipif(not K, reason="run genomeos cancer alterations")


def freq(gene: str, kind: str) -> float:
    return ((K["genes"].get(gene) or {}).get(kind) or {}).get("frequency", 0.0)


def test_the_study_has_both_profiles():
    assert K["cna_profile"] == "msk_impact_2017_cna"
    assert K["sv_profile"] == "msk_impact_2017_structural_variants"
    assert K["samples"] > 10_000


def test_the_canonical_deep_deletion_and_amplification_are_found():
    """CDKN2A is deleted and ERBB2 amplified, or the table is not describing cancer."""
    assert freq("CDKN2A", "deep_deletion") > 0.05
    assert freq("CDKN2A", "deep_deletion") > freq("CDKN2A", "amplification") * 20
    assert freq("ERBB2", "amplification") > 0.03
    assert freq("PTEN", "deep_deletion") > 0.02
    assert freq("RB1", "deep_deletion") > 0.01


def test_amplification_finds_genes_the_mutation_table_calls_passengers():
    """The reason this table exists: frequency by mutation alone under-counts.

    CCND1 and MYC are mutated in well under 1% of these tumours and amplified
    in several percent. A pipeline reading only the mutation table sees two
    passengers.
    """
    if not MUT:
        pytest.skip("run genomeos cancer distil")
    for gene in ("CCND1", "MYC"):
        mutated = MUT["genes"][gene]["frequency"]
        assert mutated < 0.01
        assert freq(gene, "amplification") > 4 * mutated


def test_the_cancer_type_is_where_the_frequency_lives():
    """ERBB2 amplification is 4% of all tumours and 14% of breast tumours."""
    breast = K["genes"]["ERBB2"]["amplification"]["by_cancer_type"]["Breast Cancer"]
    assert breast > 0.10 > freq("ERBB2", "amplification")
    glioma = K["genes"]["CDKN2A"]["deep_deletion"]["by_cancer_type"]["Glioma"]
    assert glioma > 0.25


def test_the_recurrent_fusions_have_their_known_partners():
    partners = dict(K["genes"]["ALK"]["fusion_partners"])
    assert partners.get("EML4", 0) > 20, "EML4-ALK is the fusion of lung adenocarcinoma"
    ret = dict(K["genes"]["RET"]["fusion_partners"])
    assert {"KIF5B", "CCDC6"} <= set(ret)


def test_a_panel_studys_silence_is_recorded_as_silence():
    """MSK-IMPACT calls only deep deletions and amplifications, not shallow ones."""
    kinds = {k for entry in K["genes"].values() for k in entry}
    assert "gain" not in kinds and "shallow_deletion" not in kinds
    assert "absence of a call" in K["note"]
