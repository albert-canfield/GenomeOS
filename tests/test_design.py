from pathlib import Path

import pytest

from genomeos.design import IDENTITY, design, essential_genes


def test_essential_list_ships_with_the_package():
    e = essential_genes()
    assert len(e) == 684 and {"POLR2A", "RPL3", "PCNA"} <= e


@pytest.mark.skipif(not Path("data/knowledge/go-basic.obo").exists(), reason="GO ontology not downloaded")
def test_neuron_budget_is_a_small_fraction_of_the_genome():
    d = design("neuron")
    assert d.essential_genes == 684
    assert 1000 < d.minimal_genes < 5000 < d.total_genes < 20000
    assert d.budgets_bp["dense (worm-like, 5 kb/gene)"] < 30_000_000  # < 1% of the 3.1 Gb genome
    assert "NEUROG2" in d.master_regulators and "Zhang" in d.notes[0]
    assert set(IDENTITY) >= {"neuron", "neural_progenitor", "cardiomyocyte"}
    assert (
        design("cardiomyocyte").identity_genes < d.identity_genes
    )  # heart identity is a smaller library set
