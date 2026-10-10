# SPDX-License-Identifier: AGPL-3.0-or-later
"""The libRoadRunner adapter (runtime/roadrunner.py): chosen by name, the stand-in's interface, and an
honest refusal when the `roadrunner` extra is absent. ROADMAP §5 item 13, area C."""

import sys

import pytest

from genomeos.molecules import biomodels
from genomeos.runtime import roadrunner as rr_adapter
from genomeos.runtime.sbml import ENGINES, SbmlModel, SbmlRuntime, runtime

needs_roadrunner = pytest.mark.requires_extra("roadrunner")

CASCADE = """<?xml version="1.0"?>
<sbml xmlns="http://www.sbml.org/sbml/level2/version4" level="2" version="4">
<model id="toy" name="toy cascade">
 <listOfCompartments><compartment id="c" size="2"/></listOfCompartments>
 <listOfSpecies>
  <species id="Raf" name="Raf kinase" compartment="c" initialConcentration="1"/>
  <species id="MEK" name="MEK" compartment="c" initialAmount="2"/>
  <species id="MEKp" name="MEK-P" compartment="c" initialConcentration="0"/>
  <species id="Q" name="counted" compartment="c" initialAmount="0" hasOnlySubstanceUnits="true"/>
 </listOfSpecies>
 <listOfParameters><parameter id="k" value="1"/><parameter id="kd" value="0.1"/></listOfParameters>
 <listOfReactions>
  <reaction id="phos" reversible="false">
   <listOfReactants><speciesReference species="MEK"/></listOfReactants>
   <listOfProducts><speciesReference species="MEKp"/><speciesReference species="Q"/></listOfProducts>
   <listOfModifiers><modifierSpeciesReference species="Raf"/></listOfModifiers>
   <kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML">
    <apply><times/><ci>c</ci><ci>k</ci><ci>Raf</ci><ci>MEK</ci></apply></math></kineticLaw>
  </reaction>
  <reaction id="dephos" reversible="false">
   <listOfReactants><speciesReference species="MEKp"/></listOfReactants>
   <listOfProducts><speciesReference species="MEK"/></listOfProducts>
   <kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML">
    <apply><times/><ci>c</ci><ci>kd</ci><ci>MEKp</ci></apply></math></kineticLaw>
  </reaction>
 </listOfReactions>
</model></sbml>"""

FUNCTIONS = """<?xml version="1.0"?>
<sbml xmlns="http://www.sbml.org/sbml/level2/version4" level="2" version="4">
<model id="fd">
 <listOfFunctionDefinitions>
  <functionDefinition id="mm"><math xmlns="http://www.w3.org/1998/Math/MathML">
   <lambda><bvar><ci>v</ci></bvar><bvar><ci>k</ci></bvar><bvar><ci>s</ci></bvar>
    <apply><divide/><apply><times/><ci>v</ci><ci>s</ci></apply><apply><plus/><ci>k</ci><ci>s</ci></apply></apply>
   </lambda></math></functionDefinition>
 </listOfFunctionDefinitions>
 <listOfCompartments><compartment id="c" size="1"/></listOfCompartments>
 <listOfSpecies>
  <species id="S" compartment="c" initialConcentration="10"/>
  <species id="P" compartment="c" initialConcentration="0"/>
  <species id="T" compartment="c" initialConcentration="0"/>
 </listOfSpecies>
 <listOfParameters><parameter id="Vmax" value="1"/><parameter id="Km" value="2"/></listOfParameters>
 <listOfRules>
  <rateRule variable="T"><math xmlns="http://www.w3.org/1998/Math/MathML"><cn>1</cn></math></rateRule>
 </listOfRules>
 <listOfReactions>
  <reaction id="r1" reversible="false">
   <listOfReactants><speciesReference species="S"/></listOfReactants>
   <listOfProducts><speciesReference species="P"/></listOfProducts>
   <kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML">
    <apply><ci>mm</ci><ci>Vmax</ci><ci>Km</ci><ci>S</ci></apply></math></kineticLaw>
  </reaction>
 </listOfReactions>
</model></sbml>"""


@pytest.fixture
def cascade(tmp_path):
    path = tmp_path / "toy.xml"
    path.write_text(CASCADE)
    return path


def test_the_stand_in_stays_the_default_and_engines_are_named(cascade):
    model = SbmlModel.from_file(cascade)
    assert model.source == str(cascade)
    assert ENGINES[0] == "standin" and "roadrunner" in ENGINES
    assert type(runtime(model)) is SbmlRuntime
    with pytest.raises(ValueError, match="unknown SBML engine"):
        runtime(model, "copasi")


def test_a_missing_extra_is_reported_never_substituted(cascade, monkeypatch):
    """Without libRoadRunner the adapter says which extra to install; the stand-in is not run in its place."""
    monkeypatch.setitem(sys.modules, "roadrunner", None)  # importing it now raises ImportError
    assert rr_adapter.available() is False
    model = SbmlModel.from_file(cascade)
    with pytest.raises(rr_adapter.EngineUnavailableError, match=r"uv sync --extra roadrunner"):
        runtime(model, "roadrunner")
    with pytest.raises(ImportError, match=r"--extra roadrunner"):
        biomodels.run(cascade, duration=1, dt=0.01, engine="roadrunner")
    assert biomodels.run(cascade, duration=1, dt=0.01)["engine"] == "standin"


@needs_roadrunner
def test_same_interface_same_grid_same_answer_on_a_non_stiff_model(cascade):
    """Concentrations, an amount-only species and the compartment volume read the same on both engines."""
    model = SbmlModel.from_file(cascade)
    standin = SbmlRuntime(model).run(duration=20, dt=0.001, record_every=500)
    adapter = runtime(model, "roadrunner")
    assert adapter.engine == "roadrunner" and adapter.dynamic == SbmlRuntime(model).dynamic
    traj = adapter.run(duration=20, dt=0.001, record_every=500)
    assert traj.species == standin.species and len(traj.times) == len(standin.times) == 41
    assert max(abs(a - b) for a, b in zip(traj.times, standin.times, strict=True)) < 1e-9
    for s in traj.species:
        for a, b in zip(traj.levels[s], standin.levels[s], strict=True):
            assert abs(a - b) <= 1e-4 * max(1.0, abs(b)), s
    assert abs(traj.levels["MEK"][0] - 1.0) < 1e-12  # initialAmount 2 in a compartment of 2
    assert traj.levels["Q"][-1] > 1.0  # an amount, not divided by the volume


@needs_roadrunner
def test_function_definitions_and_rate_rules(tmp_path):
    path = tmp_path / "fd.xml"
    path.write_text(FUNCTIONS)
    traj = runtime(SbmlModel.from_file(path), "roadrunner").run(duration=5, dt=0.01, record_every=100)
    s, p, t = traj.levels["S"][-1], traj.levels["P"][-1], traj.levels["T"][-1]
    assert abs(s + p - 10) < 1e-6 and p > 3 and abs(t - 5) < 1e-6


@needs_roadrunner
def test_a_knockout_and_an_initial_value_reach_libroadrunner(cascade):
    base = biomodels.run(cascade, duration=50, dt=0.01, engine="roadrunner")
    assert base["engine"] == "roadrunner" and "libRoadRunner" in base["evidence"]
    assert base["levels"]["MEKp"]["final"] > 0.8
    ko = biomodels.run(cascade, duration=50, dt=0.01, knockout=["RAF"], engine="roadrunner")
    assert ko["knocked_out"] == ["Raf"] and "Raf" not in ko["levels"]
    assert ko["levels"]["MEKp"]["final"] < 1e-9 and abs(ko["levels"]["MEK"]["final"] - 1.0) < 1e-9
    adapter = runtime(SbmlModel.from_file(cascade), "roadrunner")
    held = adapter.run(duration=10, dt=0.01, initial={"Raf": 0.0})
    assert held.levels["MEKp"][-1] < 1e-9
    again = adapter.run(duration=10, dt=0.01)  # nothing set for one run leaks into the next
    assert again.levels["MEKp"][-1] > 0.5
    settings = adapter.settings()
    assert settings["integrator"] == "cvode" and settings["libroadrunner"] == rr_adapter.version()
