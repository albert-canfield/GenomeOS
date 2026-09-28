# SPDX-License-Identifier: AGPL-3.0-or-later
"""R9 follow-up: every AlphaGenome client asks for one named model version, and answers from before the
pin (none requested) and after stay side by side, labelled, with a result that reads both saying so.
Stub clients only: no request, no network."""

import ast
import enum
import inspect
import json
import sys
import types
from pathlib import Path

import pytest

from genomeos import manifest as mf
from genomeos.predict import alphagenome_adapter as ag
from genomeos.predict import splice_sites
from genomeos.predict.enhancer_target import model_version_of, model_versions, score_element, summarise


class _Version(enum.Enum):
    ALL_FOLDS = enum.auto()
    FOLD_0 = enum.auto()


class _Client:
    def __init__(self, model_version):
        self._model_version = model_version.name if model_version is not None else None


@pytest.fixture
def fake_alphagenome(monkeypatch):
    """A stand-in `alphagenome` package whose create records its arguments and returns a stub client."""
    calls: list[tuple[tuple, dict]] = []
    dna_client = types.ModuleType("alphagenome.models.dna_client")
    dna_client.ModelVersion = _Version
    dna_client.SEQUENCE_LENGTH_1MB = 1_048_576

    def create(*a, **k):
        calls.append((a, k))
        return _Client(k.get("model_version"))

    dna_client.create = create
    models = types.ModuleType("alphagenome.models")
    models.dna_client = dna_client
    models.variant_scorers = types.ModuleType("alphagenome.models.variant_scorers")
    data = types.ModuleType("alphagenome.data")
    data.genome = types.ModuleType("alphagenome.data.genome")
    pkg = types.ModuleType("alphagenome")
    pkg.models, pkg.data = models, data
    for name, mod in {
        "alphagenome": pkg,
        "alphagenome.models": models,
        "alphagenome.models.dna_client": dna_client,
        "alphagenome.models.variant_scorers": models.variant_scorers,
        "alphagenome.data": data,
        "alphagenome.data.genome": data.genome,
    }.items():
        monkeypatch.setitem(sys.modules, name, mod)
    return calls


def test_create_client_requests_the_project_model_version(fake_alphagenome):
    client = ag.create_client("key", timeout=300)
    ((args, kwargs),) = fake_alphagenome
    assert args == ("key",) and kwargs == {"model_version": _Version.ALL_FOLDS, "timeout": 300}
    assert client._model_version == ag.ALPHAGENOME_MODEL_VERSION == "ALL_FOLDS"


def test_the_live_scorer_and_splice_client_go_through_the_pin(fake_alphagenome, monkeypatch):
    score = ag.AlphaGenomeAdapter(api_key="key")._live_scorer(threshold=0.0)
    assert score.model["model_version"] == "ALL_FOLDS" and score.model["model_version_note"] is None
    monkeypatch.setenv(ag.KEY_VAR, "key")
    assert splice_sites.client_factory()._model_version == "ALL_FOLDS"
    assert all(k["model_version"] is _Version.ALL_FOLDS for _, k in fake_alphagenome)
    assert len(fake_alphagenome) == 2


def test_run_metadata_without_a_client_records_the_constant():
    m = ag.run_metadata(date="2026-09-28")
    assert m["model_version"] == ag.ALPHAGENOME_MODEL_VERSION and m["model_version_note"] is None


def test_no_client_is_created_outside_the_adapter():
    """Every dna_client.create call goes through create_client, so none can skip the version."""
    offenders = []
    for root in ("genomeos", "scripts"):
        for f in Path(root).rglob("*.py"):
            for node in ast.walk(ast.parse(f.read_text(encoding="utf-8"))):
                fn = getattr(node, "func", None) if isinstance(node, ast.Call) else None
                if (
                    isinstance(fn, ast.Attribute)
                    and fn.attr == "create"
                    and isinstance(fn.value, ast.Name)
                    and fn.value.id == "dna_client"
                ):
                    offenders.append(f"{f}:{node.lineno}")
    assert offenders == ["genomeos/predict/alphagenome_adapter.py:" + str(_create_line())]


def _create_line() -> int:
    src = inspect.getsource(ag)
    return next(i for i, ln in enumerate(src.splitlines(), 1) if "return dna_client.create(" in ln)


def test_the_installed_client_takes_model_version_and_names_the_constant():
    dna_client = pytest.importorskip("alphagenome.models.dna_client")
    if not hasattr(dna_client, "ModelVersion"):
        pytest.skip("a stand-in alphagenome is loaded")
    assert "model_version" in inspect.signature(dna_client.create).parameters
    assert ag.ALPHAGENOME_MODEL_VERSION in dna_client.ModelVersion.__members__


# --- the cache keeps both, labelled ----------------------------------------------------------------


def _write(cache: Path, element: str, model: dict | None) -> None:
    hit = {"id": element, "chrom": "chr21", "start": 10, "end": 14, "length": 4, "genes": []}
    if model is not None:
        hit["model"] = model
    (cache / "chr21").mkdir(parents=True, exist_ok=True)
    (cache / "chr21" / f"{element}.json").write_text(json.dumps(hit))


def test_cached_answers_before_and_after_the_pin_are_kept_and_labelled(tmp_path):
    _write(tmp_path, "OLD", None)  # the sweep: no run record at all
    _write(tmp_path, "MID", {"model_version": None})  # 2026-09-28 before the pin: record, none requested

    def never(*a):
        raise AssertionError("a cached answer must not be asked again")

    def pinned(chrom, pos, ref, alt):
        return [("G", "liver", -0.3)]

    pinned.model = {"model_version": "ALL_FOLDS"}
    rows = [
        score_element(never, lambda locus: "ACGTAC", "chr21", "OLD", 10, 14, cache=tmp_path),
        score_element(never, lambda locus: "ACGTAC", "chr21", "MID", 10, 14, cache=tmp_path),
        score_element(pinned, lambda locus: "ACGTAC", "chr21", "NEW", 10, 14, cache=tmp_path),
    ]
    assert [r["model_version"] for r in rows] == ["unrequested", "unrequested", "ALL_FOLDS"]
    assert model_version_of(json.loads((tmp_path / "chr21" / "NEW.json").read_text())) == "ALL_FOLDS"
    assert json.loads((tmp_path / "chr21" / "OLD.json").read_text()).get("model") is None  # untouched
    assert summarise(rows)["answers_by_model_version"] == {"unrequested": 2, "ALL_FOLDS": 1}


def test_a_result_mixing_both_says_so_in_its_manifest():
    block = mf.answers_model_dependency(
        "alphagenome", model_versions([{}, {}, {"model_version": "ALL_FOLDS"}])
    )
    assert block["answers_by_model_version"] == {"unrequested": 2, "ALL_FOLDS": 1}
    assert block["mixed"].startswith("mixed: ") and "not guaranteed" in block["mixed"]
    assert block["model_version"] is None and block["unpinned"].startswith("unpinned: ")
    assert mf.validate({"model_dependencies": [block]}) == [m for m in mf.validate({}) if m]


def test_a_result_of_pinned_answers_only_pins_and_one_of_old_answers_stays_unpinned():
    pinned = mf.answers_model_dependency("alphagenome", {"ALL_FOLDS": 5})
    assert pinned["model_version"] == "ALL_FOLDS" and pinned["mixed"] is None and pinned["unpinned"] is None
    old = mf.answers_model_dependency("alphagenome", {"unrequested": 5})
    assert old["model_version"] is None and old["mixed"] is None and old["unpinned"].startswith("unpinned: ")
    for block in (pinned, old):
        assert not [p for p in mf.validate({"model_dependencies": [block]}) if "model dependency" in p]
