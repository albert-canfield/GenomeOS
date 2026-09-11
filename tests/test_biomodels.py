# SPDX-License-Identifier: AGPL-3.0-or-later
"""Kinetic pathways from BioModels on the in-house engine: run, knock out, compare (no network)."""

from genomeos.molecules.biomodels import compare, run, species_matching

XML = """<?xml version="1.0"?>
<sbml xmlns="http://www.sbml.org/sbml/level2/version4" level="2" version="4">
<model id="toy" name="toy cascade">
 <listOfCompartments><compartment id="c" size="1"/></listOfCompartments>
 <listOfSpecies>
  <species id="Raf" name="Raf kinase" compartment="c" initialConcentration="1">
   <annotation><rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
    <rdf:Description rdf:about="#Raf"><bqbiol:is><rdf:Bag><rdf:li rdf:resource="http://identifiers.org/uniprot/P04049"/></rdf:Bag></bqbiol:is></rdf:Description>
   </rdf:RDF></annotation>
  </species>
  <species id="MEK" name="MEK" compartment="c" initialConcentration="1"/>
  <species id="MEKp" name="MEK-P" compartment="c" initialConcentration="0"/>
 </listOfSpecies>
 <listOfParameters><parameter id="k" value="1"/><parameter id="kd" value="0.1"/></listOfParameters>
 <listOfReactions>
  <reaction id="phos" reversible="false">
   <listOfReactants><speciesReference species="MEK"/></listOfReactants>
   <listOfProducts><speciesReference species="MEKp"/></listOfProducts>
   <listOfModifiers><modifierSpeciesReference species="Raf"/></listOfModifiers>
   <kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML">
    <apply><times/><ci>k</ci><ci>Raf</ci><ci>MEK</ci></apply></math></kineticLaw>
  </reaction>
  <reaction id="dephos" reversible="false">
   <listOfReactants><speciesReference species="MEKp"/></listOfReactants>
   <listOfProducts><speciesReference species="MEK"/></listOfProducts>
   <kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML">
    <apply><times/><ci>kd</ci><ci>MEKp</ci></apply></math></kineticLaw>
  </reaction>
 </listOfReactions>
</model></sbml>"""


def test_run_knockout_and_compare(tmp_path):
    f = tmp_path / "toy.xml"
    f.write_text(XML)
    base = run(f, duration=50, dt=0.01)
    assert base["model"] == "toy" and base["levels"]["MEKp"]["final"] > 0.8  # Raf drives MEK to MEK-P
    from genomeos.runtime.sbml import SbmlModel

    m = SbmlModel.from_file(f)
    assert species_matching(m, "raf") == ["Raf"]
    assert m.annotations["Raf"] == ["http://identifiers.org/uniprot/P04049"]
    assert species_matching(m, "RAF1", accession="P04049") == ["Raf"] and species_matching(m, "RAF1") == []
    ko = run(f, duration=50, dt=0.01, knockout=["RAF"])
    assert ko["knocked_out"] == ["Raf"] and ko["levels"]["MEKp"]["final"] < 1e-6
    changed = compare(base, ko)
    names = {r["species"]: r for r in changed}
    assert "MEKp" in names and names["MEKp"]["relative_change"] < -0.8
    assert names["MEK"]["relative_change"] > 0.8
