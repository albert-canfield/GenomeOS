# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
# ruff: noqa: E501  (BioLang source kept one block per line)
"""BioLang v0.4 stage 1: compartments, locations, targeting signals, transports, the regime record.

The semantics under test are the three facts a located program must obey (docs/BIOLANG-v0.4-ECONOMY.md
§4): everything is somewhere, a protein only arrives where a transport carries it, and a rule acts
where its target is. What must never happen is a silent answer: a missing location, a rule across
compartments and a route that does not exist are all reported.
"""

from __future__ import annotations

import pytest

from genomeos.ir import Module, Regime
from genomeos.lang import BioLangError, parse
from genomeos.lang.located import layout
from genomeos.runtime.located import LocatedRuntime

CELL = """
module test.located
compartment Extracellular { }
compartment PlasmaMembrane { parent: Extracellular; membrane: yes }
compartment Cytosol { parent: PlasmaMembrane; volume: 0.54; translation: yes }
compartment Nucleus { parent: Cytosol; volume: 0.06; genome: nuclear }
compartment Mitochondrion { parent: Cytosol; volume: 0.22; genome: chrM; translation: yes; copies: 1000 }
transport Pore { from: Nucleus; to: Cytosol; cargo: mRNA; capacity: 1e4; affinity: 1e3 }
transport TOM { from: Cytosol; to: Mitochondrion; cargo: signal = presequence; capacity: 1e4; affinity: 1e3 }
"""

TWO_GENOMES = (
    CELL
    + """
gene MT-ND1 { locus: chrM:3306-4262(+); basal: 1; produces: ND1p }
gene SDHB { location: Nucleus; basal: 1; produces: SDHBp }
protein ND1p { location: Mitochondrion }
protein SDHBp { location: Mitochondrion; signals: presequence }
protein CII { location: Mitochondrion }
rule SDHBp binds CII { }
"""
)


def test_compartments_form_a_tree_and_carry_their_facts():
    m = parse(TWO_GENOMES)
    assert m.located
    mito = m.entities["Mitochondrion"]
    assert mito.parent == "Cytosol" and mito.genome == ["chrM"] and mito.translation and mito.copies == 1000
    assert mito.volume == 0.22 and m.entities["PlasmaMembrane"].membrane
    pl = layout(m).places
    assert pl.roots == ["Extracellular"] and pl.adjacent("Cytosol", "Mitochondrion")
    assert not pl.adjacent("Nucleus", "Mitochondrion")
    assert pl.faces("PlasmaMembrane") == ["PlasmaMembrane", "Extracellular", "Cytosol"]
    m2 = Module.from_dict(m.to_dict())  # the IR round-trips
    assert m2.entities["Mitochondrion"].genome == ["chrM"] and m2.entities["SDHBp"].signals == ["presequence"]


def test_a_gene_is_placed_by_its_genome_and_its_mrna_finds_a_ribosome():
    lay = layout(parse(TWO_GENOMES))
    assert lay.errors == []
    assert lay.gene_site["MT-ND1"] == "Mitochondrion"  # from the locus, through the compartment's genome
    assert lay.gene_site["SDHB"] == "Nucleus"
    assert lay.mrna_sites["SDHB"] == ["Nucleus", "Cytosol"]  # the nuclear pore
    assert lay.synthesis == {"ND1p": ["Mitochondrion"], "SDHBp": ["Cytosol"]}
    assert lay.reach["SDHBp"] == {"Cytosol", "Mitochondrion"}


@pytest.mark.parametrize(
    ("body", "message"),
    [
        ("gene X { produces: Xp }\nprotein Xp { location: Cytosol }", "has no location"),
        ("gene X { location: Nucleus; produces: Xp }\nprotein Xp { }", "protein 'Xp' has no location"),
        ("gene MT-ND1 { locus: chrM:3306-4262(+); location: Nucleus }", "does not read chrM"),
        (
            "gene X { location: Nucleus; produces: Xp }\nprotein Xp { location: Nucleus }\nprotein T { location: Mitochondrion }\nrule T activates X { }",
            "acts across compartments",
        ),
        ("transport Bad { from: Nucleus; to: Mitochondrion; cargo: mRNA; capacity: 1 }", "not adjacent"),
        (
            "transport Bad { from: Nucleus; to: Cytosol; cargo: mRNA; capacity: 1; via: Nope }",
            "undeclared protein",
        ),
    ],
)
def test_the_compiler_refuses_a_program_that_cannot_mean_anything(body, message):
    with pytest.raises(BioLangError, match=message):
        parse(CELL + body)


def test_a_nuclear_gene_needs_its_mrna_exported():
    body = "gene X { location: Nucleus; produces: Xp }\nprotein Xp { location: Cytosol }"
    without_pore = CELL.replace(
        "transport Pore { from: Nucleus; to: Cytosol; cargo: mRNA; capacity: 1e4; affinity: 1e3 }\n", ""
    )
    with pytest.raises(BioLangError, match="reaches no compartment with translation"):
        parse(without_pore + body)
    assert layout(parse(CELL + body)).errors == []


def test_locations_and_transports_need_compartments():
    with pytest.raises(BioLangError, match="need declared compartments"):
        parse("module t\ngene X { location: Nucleus }")


def test_two_genomes_one_cell_and_what_breaks_when_a_route_is_cut():
    m = parse(TWO_GENOMES)
    wt = LocatedRuntime(m).run(hours=24)
    assert wt.where("ND1p") == ["Mitochondrion"]  # made inside, never anywhere else
    assert wt.final("SDHBp@Mitochondrion") > 0 and wt.final("CII@Mitochondrion") > 0
    assert wt.stranded == [] and wt.ectopic == []
    assert [t["compartment"] for t in wt.transit] == ["Cytosol"]  # the precursor on its way in
    rho0 = LocatedRuntime(m, knockouts={"chrM"}).run(hours=24)
    assert rho0.where("ND1p") == [] and rho0.final("CII@Mitochondrion") > 0
    closed = LocatedRuntime(m, knockouts={"TOM"}).run(hours=24)
    assert closed.final("CII@Mitochondrion") == 0 and closed.transports_used["TOM"] == 0
    assert closed.where("ND1p") == ["Mitochondrion"]  # the mitochondrial genome does not need import here
    assert [s["protein"] for s in closed.stranded] == ["SDHBp"]


def test_a_mislocalisation_changes_the_prediction_and_is_named():
    m = parse(TWO_GENOMES)
    lost = LocatedRuntime(m, knockouts={"SDHBp:presequence"}).run(hours=24)
    assert lost.final("CII@Mitochondrion") == 0  # the prediction moves
    assert lost.stranded[0]["protein"] == "SDHBp" and lost.stranded[0]["declared"] == "Mitochondrion"
    assert lost.ectopic[0]["compartment"] == "Cytosol"  # where it piled up instead
    chain = lost.confidence["SDHBp@Mitochondrion"]["chain"]
    assert "no open route to Mitochondrion" in chain  # the chain says why, at confidence 0


def test_a_knockout_that_names_nothing_is_an_error_not_a_no_op():
    with pytest.raises(ValueError, match="name nothing in the program"):
        LocatedRuntime(parse(TWO_GENOMES), knockouts={"NDUFS1"})


def test_one_transport_is_shared_by_its_cargos():
    """A second cargo through the same pore slows the first: capacity is finite, not per molecule."""
    one = (
        CELL
        + "gene A { location: Nucleus; basal: 10; produces: Ap }\nprotein Ap { location: Mitochondrion; signals: presequence }\n"
    )
    two = (
        one
        + "gene B { location: Nucleus; basal: 1000; produces: Bp }\nprotein Bp { location: Mitochondrion; signals: presequence }\n"
    )
    alone = LocatedRuntime(parse(one)).run(hours=6)
    crowded = LocatedRuntime(parse(two)).run(hours=6)
    assert crowded.final("Ap@Mitochondrion") < alone.final("Ap@Mitochondrion")
    assert crowded.final("Ap@Cytosol") > alone.final("Ap@Cytosol")  # it queues in the cytosol instead


def test_the_regime_is_declared_and_recorded():
    src = TWO_GENOMES + "regime single_cell { treatment: auto; threshold: 100; units: copies; seed: 3 }\n"
    m = parse(src)
    assert m.regime.treatment == "auto" and m.regime.units == "copies" and m.regime.seed == 3
    res = LocatedRuntime(m).run(hours=6, dt=0.05)
    rec = res.regime
    assert rec["declared"] and rec["treatment"] == "auto" and rec["seed"] == 3
    assert rec["integrator"].startswith("tau-leap") and rec["reaction_firings"]["stochastic"] > 0
    assert rec["species_stochastic_at_end"]  # low-copy species, named in the output
    m2 = Module.from_dict(m.to_dict())
    assert m2.regime.threshold == 100


def test_stochastic_treatment_is_refused_in_arbitrary_units():
    m = parse(TWO_GENOMES)
    with pytest.raises(ValueError, match="units: copies"):
        LocatedRuntime(m, regime=Regime(treatment="stochastic"))
    with pytest.raises(ValueError, match="not implemented"):
        LocatedRuntime(m, regime=Regime(update="asynchronous"))
    with pytest.raises(ValueError, match="stage 2"):
        LocatedRuntime(m, regime=Regime(allocation="priority"))


def test_a_cell_with_no_genome_still_runs():
    """The red blood cell: no nucleus, no DNA, no gene; the runtime must not assume a stored program."""
    src = """
module test.red_cell
compartment Extracellular { }
compartment PlasmaMembrane { parent: Extracellular; membrane: yes }
compartment Cytosol { parent: PlasmaMembrane }
protein HBA { location: Cytosol; initial: 100; half_life: 2856 h }
protein HBB { location: Cytosol; initial: 100; half_life: 2856 h }
protein HbA { location: Cytosol; half_life: 2856 h }
rule HBA binds HbA { }
rule HBB binds HbA { }
"""
    m = parse(src)
    assert m.genes() == [] and m.located
    res = LocatedRuntime(m).run(hours=48)
    assert res.final("HbA@Cytosol") > 90 and res.stranded == []
    assert res.final("HBA@Cytosol") + res.final("HbA@Cytosol") <= 100.001  # nothing made from nothing


def test_confidence_follows_the_whole_chain_to_the_weakest_link():
    src = CELL.replace(
        "transport TOM { from: Cytosol; to: Mitochondrion; cargo: signal = presequence; capacity: 1e4; affinity: 1e3 }",
        "transport TOM { from: Cytosol; to: Mitochondrion; cargo: signal = presequence; capacity: 1e4; affinity: 1e3\n"
        '  evidence: inferred "capacity not measured"; confidence: 0.3 }',
    )
    src = src.replace(
        "compartment Mitochondrion {",
        'compartment Mitochondrion { evidence: curated "Alberts"; confidence: 0.9;',
    )
    src = src.replace(
        "compartment Cytosol {", 'compartment Cytosol { evidence: curated "Alberts"; confidence: 0.9;'
    )
    src = src.replace(
        "compartment Nucleus {", 'compartment Nucleus { evidence: curated "Alberts"; confidence: 0.9;'
    )
    src += (
        'gene SDHB { location: Nucleus; basal: 1; produces: SDHBp; evidence: curated "MitoCarta3.0"; confidence: 0.9 }\n'
        'protein SDHBp { location: Mitochondrion; signals: presequence; evidence: curated "UniProt"; confidence: 0.9 }\n'
    )
    src = src.replace("transport Pore {", 'transport Pore { evidence: curated "Alberts"; confidence: 0.8;')
    res = LocatedRuntime(parse(src)).run(hours=12)
    inside = res.confidence["SDHBp@Mitochondrion"]
    assert inside["weakest"] == "transport TOM" and inside["score"] == pytest.approx(0.3 * 0.4, abs=1e-6)
    assert res.confidence["SDHBp@Cytosol"]["score"] > inside["score"]  # the precursor is better grounded
