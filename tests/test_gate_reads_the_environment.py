# SPDX-License-Identifier: AGPL-3.0-or-later
"""An optional-feature gate must REPORT, and it must report about the environment.

Two defects of one shape, found on 2026-10-03. `importlib.util.find_spec` consults `sys.modules`
before it asks the filesystem, so the answer to "is this optional client installed" moved with
whatever had been left under that name -- and when what was left there has no spec, `find_spec` does
not answer at all, it RAISES `ValueError`. A gate that throws is not a gate that says no.

    genomeos/predict/alphagenome_adapter.py:181   `status()`, the gate in front of a paid client
    genomeos/therapeutics/binding.py:125          `MhcflurryPredictor.installed()`

The reading is now `genomeos.importable.importable`, which asks `sys.meta_path`'s finders -- what an
actual `import` consults AFTER the `sys.modules` shortcut -- and is injectable at each call site so
that these tests can prove a gate answers about what it was handed. Every production call passes
nothing and gets the live reading, which is pinned by source below.

EVERY PLANT HERE IS RUN AGAINST THE UNREPAIRED CODE AND SHOWN TO FAIL. The counterfactual is the
module's own current source with the old reading put back in place of the new one, loaded under a
package-qualified name so its relative imports still resolve. It therefore keeps the repaired
signature while ignoring the argument, which is the one thing a passing injection test must be able
to see: a function still consulting `sys.modules` must not be able to pass itself off as one that
does not.

Nothing here installs mhcflurry, imports alphagenome, makes a network call or spends anything.
"""

from __future__ import annotations

import ast
import importlib.util
import inspect
import sys
import types
from pathlib import Path

import pytest

from genomeos.importable import importable
from genomeos.therapeutics import binding

REPO = Path(__file__).resolve().parents[1]

#: The reading each gate took until 2026-10-03, as the counterfactual puts it back. Written out so a
#: rename cannot leave these plants measuring something that is no longer the defect.
OLD_READING = b'return importlib.util.find_spec("mhcflurry") is not None'
NEW_READING = b'return (importable if present is None else present)("mhcflurry")'


def _load(source: bytes, name: str, tmp_path: Path):
    """Import `source` as its own module beside the real one, under a package-qualified name.

    The name matters: `binding.py` uses relative imports (`from .evidence import prediction`), which
    resolve through the module's package. Loaded as a bare top-level name the copy would not import
    at all, and a counterfactual that cannot be loaded proves nothing.
    """
    path = tmp_path / f"{name.rsplit('.', 1)[-1]}.py"
    path.write_bytes(source)
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    try:
        spec.loader.exec_module(mod)
    except BaseException:
        sys.modules.pop(name, None)
        raise
    return mod


@pytest.fixture
def unrepaired(tmp_path):
    """`binding.py` with the old reading put back: the repaired signature, ignoring the argument."""
    source = (REPO / "genomeos/therapeutics/binding.py").read_bytes()
    assert NEW_READING in source, "the line the counterfactual replaces is no longer there"
    planted = source.replace(NEW_READING, b"import importlib.util\n\n        " + OLD_READING)
    assert planted != source
    mod = _load(planted, "genomeos.therapeutics.binding_unrepaired", tmp_path)
    yield mod
    sys.modules.pop("genomeos.therapeutics.binding_unrepaired", None)


# --- the reading itself ------------------------------------------------------------------------------


def test_a_spec_less_stub_does_not_move_the_reading(tmp_path, monkeypatch):
    """The defect, exercised and not described: `find_spec` must really raise, or this is no plant.

    The module is REAL and its own -- a file on `sys.path` -- so the environment genuinely provides it
    and the only thing the stub changes is `sys.modules`. Stubbing an installed library instead would
    answer the same question while displacing something this suite stands on.
    """
    name = "a_module_the_gate_tests_really_provide"
    (tmp_path / f"{name}.py").write_text("VALUE = 1\n")
    monkeypatch.syspath_prepend(str(tmp_path))
    assert importable(name) is True, "the finders see a file on sys.path"

    monkeypatch.setitem(sys.modules, name, types.ModuleType(name))
    with pytest.raises(ValueError, match="__spec__"):
        importlib.util.find_spec(name)
    assert importable(name) is True, "and still answer about the environment"


def test_a_stub_carrying_a_spec_does_not_make_an_absent_module_read_as_present(monkeypatch):
    """The other direction: a module nothing provides must read absent however `sys.modules` looks."""
    name = "no_distribution_provides_this_name_either"
    spec = importlib.util.spec_from_loader(name, loader=None)
    monkeypatch.setitem(sys.modules, name, importlib.util.module_from_spec(spec))

    assert importlib.util.find_spec(name) is not None, "the old reading is satisfied by the stub"
    assert importable(name) is False


def test_a_finder_that_raises_is_answered_and_not_propagated(monkeypatch):
    """An unaskable question is answered False, because the whole defect was a gate that raised."""

    class Hostile:
        def find_spec(self, name, path, target=None):
            raise ValueError("this finder refuses to answer")

    monkeypatch.setattr(sys, "meta_path", [Hostile()])
    assert importable("anything_at_all") is False


def test_a_later_finder_is_still_asked_after_an_earlier_one_raises(monkeypatch):
    """`continue`, not `return`: one broken finder must not hide a module another one provides."""

    class Hostile:
        def find_spec(self, name, path, target=None):
            raise TypeError("nor this one")

    class Helpful:
        def find_spec(self, name, path, target=None):
            return importlib.util.spec_from_loader(name, loader=None)

    monkeypatch.setattr(sys, "meta_path", [Hostile(), Helpful()])
    assert importable("provided_by_the_second_finder") is True


# --- the MHCflurry gate ------------------------------------------------------------------------------


def test_the_mhcflurry_gate_reports_under_a_spec_less_stub(monkeypatch, unrepaired):
    """The repaired gate answers; the unrepaired one raises `ValueError` instead of answering."""
    monkeypatch.setitem(sys.modules, "mhcflurry", types.ModuleType("mhcflurry"))

    assert binding.MhcflurryPredictor.installed() in (True, False)
    with pytest.raises(ValueError, match="__spec__"):
        unrepaired.MhcflurryPredictor.installed()


def test_the_whole_predictor_still_answers_under_a_spec_less_stub(monkeypatch, unrepaired):
    """Through the public surface, because `installed` is read by `available` and by `predict`."""
    monkeypatch.setitem(sys.modules, "mhcflurry", types.ModuleType("mhcflurry"))
    p = binding.MhcflurryPredictor()

    assert p.available in (True, False)
    assert p.predict(["KLVFFAEDV"], ["HLA-A*02:01"]).available is False

    with pytest.raises(ValueError, match="__spec__"):
        _ = unrepaired.MhcflurryPredictor().available


def test_the_injected_reading_decides_and_sys_modules_does_not(monkeypatch, unrepaired):
    """The plant the repair has to survive: a reading that is IGNORED must be visible as ignored.

    A spec-bearing stub is put under the name, so `sys.modules` says present while the injected
    reading says absent. The repaired gate must follow what it was handed; the unrepaired one answers
    from `sys.modules` and so contradicts its own argument, which is how it is caught.
    """
    spec = importlib.util.spec_from_loader("mhcflurry", loader=None)
    monkeypatch.setitem(sys.modules, "mhcflurry", importlib.util.module_from_spec(spec))
    assert importlib.util.find_spec("mhcflurry") is not None, "sys.modules says present"

    assert binding.MhcflurryPredictor.installed(lambda _m: False) is False
    assert unrepaired.MhcflurryPredictor.installed(lambda _m: False) is True, (
        "the counterfactual ignores the argument, which is what makes the line above a plant"
    )


def test_the_injected_reading_is_asked_about_mhcflurry_and_nothing_else():
    """What the gate asks for, so an injection cannot pass by being asked the wrong question."""
    asked: list[str] = []

    def watched(module: str) -> bool:
        asked.append(module)
        return True

    assert binding.MhcflurryPredictor.installed(watched) is True
    assert asked == ["mhcflurry"]


def test_every_production_call_takes_the_live_reading_and_is_told_nothing(unrepaired):
    """A guard keyed to what its caller supplies is not a guard, so the call sites are pinned.

    The injection exists for the tests above. If `available` or `predict` ever passed a reading of its
    own, every MHCflurry decision in the application would be whatever that argument said, and the
    plants above would be measuring a parameter nothing production reaches.
    """
    source = (REPO / "genomeos/therapeutics/binding.py").read_text()
    assert source.count("self.installed()") == 2, source
    assert "self.installed(" not in source.replace("self.installed()", ""), (
        "some production call site now supplies its own reading"
    )
    # asked of the AST and not of the text, because the docstring NAMES the old reading in order to
    # record what the defect was, and a text search cannot tell a history note from a live import
    imported = {
        alias.name
        for node in ast.walk(ast.parse(source))
        if isinstance(node, ast.Import)
        for alias in node.names
    }
    assert "importlib.util" not in imported, "the old reading's import is back in the production file"

    signature = inspect.signature(binding.MhcflurryPredictor.installed)
    assert signature.parameters["present"].default is None, signature
    assert binding.MhcflurryPredictor.installed() is importable("mhcflurry"), (
        "the default is the live reading"
    )
