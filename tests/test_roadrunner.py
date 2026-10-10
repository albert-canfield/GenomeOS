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


# ---- the score against BioModels' curation figures ------------------------------------------------
# Fixed on 2026-10-10 BEFORE either engine was run on the scored models (genomeos/molecules/
# sbml_reference.py has the reading). A change to any number below, or to any reference feature, fails
# here: the tolerance is part of the question, not a knob for the answer.
FROZEN_TOLERANCE = {"relative": 0.10, "pixels": 3, "period_relative": 0.03}
FROZEN_AGREEMENT = 0.02
FROZEN_REFERENCES_SHA256 = "66322616dcfb533c195a0313a6560a47f4c06d55b309b443a9c27c99cefe4884"
SCHOEBERL = "data/knowledge/biomodels/BIOMD0000000019.xml"
REPRESSILATOR = "data/models/BIOMD0000000012.xml"


def test_the_tolerance_and_the_references_are_the_ones_fixed_before_the_runs():
    import hashlib
    import json

    from genomeos.molecules import sbml_reference as ref

    assert ref.TOLERANCE == FROZEN_TOLERANCE and ref.AGREEMENT == FROZEN_AGREEMENT
    digest = hashlib.sha256(json.dumps(ref.REFERENCES, sort_keys=True).encode()).hexdigest()
    assert digest == FROZEN_REFERENCES_SHA256
    features = sum(len(d["features"]) for m in ref.REFERENCES.values() for d in m["doses"].values())
    assert features == 52 and set(ref.REFERENCES) == {
        "BIOMD0000000019",
        "BIOMD0000000010",
        "BIOMD0000000012",
        "BIOMD0000000084",
        "BIOMD0000000116",
    }


def test_each_feature_kind_reads_a_known_series_and_applies_the_rule():
    """A cosine of period 10 and amplitude 1 around 2, on a grid of 0.01: every kind, and both bounds."""
    import math

    from genomeos.molecules import sbml_reference as ref

    times = [i * 0.01 for i in range(4001)]
    xs = [2 + math.cos(2 * math.pi * t / 10) for t in times]
    axes = {"t": 0.1, "X": 0.01}

    def one(feature):
        return ref.score_feature(feature, times, xs, axes)

    peak = one(("peak", "X", (1.0, 15.0), 3.0, 10.0))
    assert [c["check"] for c in peak] == ["level", "time"] and all(c["pass"] for c in peak)
    assert abs(peak[1]["simulated"] - 10.0) < 1e-9
    assert one(("at", "X", 5.0, 1.0))[0]["pass"]
    assert one(("min", "X", (0.0, 40.0), 1.0))[0]["pass"] and one(("max", "X", (0.0, 40.0), 3.0))[0]["pass"]
    period = one(("period", "X", (0.0, 40.0), 10.0))[0]
    assert period["pass"] and abs(period["simulated"] - 10.0) < 1e-9 and period["bound"] == 0.3
    assert not one(("period", "X", (0.0, 40.0), 10.4))[0]["pass"]  # 4% off: outside 3%
    assert one(("first_peak", "X", 12.0, 20.0))[0]["simulated"] == 20.0
    # the bound is 10% of the REFERENCE: a simulated 3.0 passes against 3.33 (bound 0.333), fails against 3.34
    assert (
        one(("max", "X", (0.0, 40.0), 3.33))[0]["pass"]
        and not one(("max", "X", (0.0, 40.0), 3.34))[0]["pass"]
    )
    # near zero the pixel bound rules: 3 px x 0.01 = 0.03 around a reference of 1.02
    assert one(("min", "X", (0.0, 40.0), 1.02))[0]["bound"] == 0.102
    assert one(("at", "X", 5.0, 0.02))[0]["bound"] == 0.03 and not one(("at", "X", 5.0, 0.02))[0]["pass"]


def _score(path, model_id):
    from genomeos.molecules import sbml_reference as ref

    spec = ref.REFERENCES[model_id]
    out = {}
    for dose_name, dose in spec["doses"].items():
        model = SbmlModel.from_file(path)
        dt = spec["standin_dt"]
        traj = runtime(model, "roadrunner").run(
            duration=spec["duration"],
            dt=dt,
            record_every=max(1, round(spec["record_dt"] / dt)),
            initial=dose["initial"],
        )
        checks = []
        for feature in dose["features"]:
            checks += ref.score_feature(
                feature, traj.times, ref.series(model, traj, feature[1], dose["initial"]), spec["axes"]
            )
        out[dose_name] = [c for c in checks if not c["pass"]]
    return out


@needs_roadrunner
@pytest.mark.needs_local_data(SCHOEBERL, how="genomeos kinetic BIOMD0000000019 caches it from BioModels")
def test_schoeberl2002_reproduces_on_libroadrunner_at_all_three_egf_doses():
    """The registered expectation: every feature of the curation figure, at 50, 0.5 and 0.125 ng/ml."""
    failed = _score(SCHOEBERL, "BIOMD0000000019")
    assert len(failed) == 3 and failed == {d: [] for d in failed}, failed


@needs_roadrunner
def test_the_repressilator_reproduces_on_libroadrunner():
    failed = _score(REPRESSILATOR, "BIOMD0000000012")
    assert failed == {"as curated": []}, failed
