"""Task 2.3: SBML models from BioModels run unchanged."""

from pathlib import Path

import pytest

from genomeos.runtime.sbml import SbmlModel, SbmlRuntime, compile_math

MODEL = Path("data/models/BIOMD0000000012.xml")


def test_mathml_evaluator():
    import xml.etree.ElementTree as ET

    m = ET.fromstring(
        '<math xmlns="http://www.w3.org/1998/Math/MathML"><apply><divide/><apply><times/><ci>a</ci><cn>3</cn></apply>'
        "<apply><power/><ci>b</ci><cn>2</cn></apply></apply></math>"
    )
    f = compile_math(m)
    assert f({"a": 4.0, "b": 2.0}) == 3.0
    e = ET.fromstring(
        '<math xmlns="http://www.w3.org/1998/Math/MathML"><cn type="e-notation">1<sep/>-3</cn></math>'
    )
    assert abs(compile_math(e)({}) - 1e-3) < 1e-12


@pytest.mark.skipif(not MODEL.exists(), reason="download BIOMD0000000012 first")
def test_biomodels_repressilator_oscillates():
    model = SbmlModel.from_file(MODEL)
    assert model.name.startswith("Elowitz2000")
    assert len(model.reactions) == 12 and set(model.species) == {"PX", "PY", "PZ", "X", "Y", "Z"}
    rt = SbmlRuntime(model)
    env = rt._env(dict(model.species))
    assert abs(env["beta"] - 0.2) < 1e-9  # assignment rule tau_mRNA / tau_prot = 2 / 10
    assert abs(env["alpha"] - 216.404) < 1.0  # the file stores a rounded snapshot of the assignment rule
    traj = rt.run(duration=1000.0, dt=0.05, record_every=20)  # minutes
    for s in ("PX", "PY", "PZ"):
        assert traj.peaks(s) >= 3, s
        xs = traj.levels[s][len(traj.levels[s]) // 2 :]
        assert max(xs) - min(xs) > 100  # sustained oscillation in protein copies
    m = model.to_module()
    assert len(m.entities) == 6 and m.rules and abs(m.confidence_report()["species"] - 0.8) < 1e-9
