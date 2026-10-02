# SPDX-License-Identifier: AGPL-3.0-or-later
"""The gene-axis record is additive: the frozen deletion path computes what it computed before.

The question these tests answer is not "does the new code look additive" but "does the pre-change module,
given the same response objects at the same moment in the same process, return the same effects". So the
module as it stood before this lane is loaded from its own git blob -- pinned twice, by blob id and by the
sha256 of the bytes that come back -- and run beside the current module against one shared set of stub
responses. Two runs at different times could not tell this lane's change from anything else that moved.

The frozen quantities compared are the scorer's effects list, `enhancer_target.aggregate` over it, and
`predict_target` over that: the deletion feature's whole chain from the response to the named target. A
planted counterfactual (`test_the_comparison_catches_a_real_change`) perturbs the emitting threshold by one
character and asserts the same comparison fails, so a passing comparison cannot be passing vacuously.

No request is made and no key is read: the `alphagenome` package is not installed here, so its modules are
stubbed and the client is a fake that hands back fixed pandas-backed responses.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import subprocess
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from genomeos.predict import alphagenome_adapter as new
from genomeos.predict import enhancer_target as et

REPO = Path(__file__).resolve().parents[1]

#: The adapter exactly as it stood before this lane, pinned as a git blob so that it stays the
#: pre-change module after this lane is committed (pinning HEAD would start pinning this lane's own
#: change the moment it lands) and so that no branch, worktree or working-tree state can move it.
FROZEN_BLOB = "c1ed397eb165c5b9a3ae9ce277abaf2ce6fb0eb7"
#: The sha256 of those bytes, so the content is pinned independently of the blob id as well.
FROZEN_SHA256 = "7f52071ae218964312bfc2bf2aebf1fb42a034df6a1d7e669c556bd842de5aca"

THRESHOLD = 0.05


# --- the stub response -------------------------------------------------------------------------------


def _frame(cols: dict) -> pd.DataFrame:
    """An axis as the client builds it: a pandas frame whose index is the row's position as a string."""
    df = pd.DataFrame(cols)
    df.index = df.index.map(str)
    return df


def responses() -> list:
    """Two gene-scorer responses, shaped like the ones the pinned client constructs.

    The first carries the case the whole lane is about: two gene-axis rows, positions 0 and 2, under one
    symbol AAA with two different gene ids, which the symbol-keyed accumulator pools into one row. It also
    carries an exact zero on a K562 track (kept by the record, dropped by the effects), a value exactly at
    the emitting threshold (the counterfactual's lever), a track whose gtex_tissue is the string "nan" and
    one whose gtex_tissue is empty, and a row every one of whose values is below the threshold.

    The second carries no gene_name column at all, which is possible because the client sets that column
    only when the gene proto has the field. The effects loop iterates the symbol column, so such a
    response yields no effects whatever its values say.
    """
    var = _frame(
        {
            "biosample_name": ["K562", "K562", "HepG2", "liver", "GM12878"],
            "gtex_tissue": [None, "nan", None, "Liver", ""],
        }
    )
    first = types.SimpleNamespace(
        obs=_frame(
            {
                "gene_id": ["ENSG1.1", "ENSG2.3", "ENSG3.2", "ENSG4.1"],
                "strand": ["+", "-", "+", "+"],
                "gene_name": ["AAA", "BBB", "AAA", "CCC"],
                "gene_type": ["protein_coding", "lncRNA", "protein_coding", "protein_coding"],
            }
        ),
        var=var,
        X=np.array(
            [
                [0.60, 0.00, -0.30, 0.02, 0.05],
                [-0.80, -0.70, 0.10, 0.00, 0.40],
                [0.20, 0.15, -0.90, 0.30, -0.10],
                [0.01, 0.00, 0.02, 0.03, 0.04],
            ]
        ),
        layers={"quantiles": np.zeros((4, 5))},
    )
    second = types.SimpleNamespace(
        obs=_frame({"gene_id": ["ENSG9.1", "ENSG8.1"], "strand": ["+", "-"]}),
        var=var,
        X=np.array([[0.9, -0.9, 0.5, -0.5, 0.1], [0.3, 0.3, 0.3, 0.3, 0.3]]),
        layers=None,
    )
    return [first, second]


class FakeClient:
    """Hands back the same response objects on every call, so both modules read one shared input."""

    _model_version = "ALL_FOLDS-stub"

    def __init__(self, shared: list) -> None:
        self.shared = shared
        self.calls = 0

    def score_variant(self, interval, variant, variant_scorers):  # noqa: ARG002
        self.calls += 1
        return self.shared


@pytest.fixture
def stub_alphagenome(monkeypatch):
    """The parts of the uninstalled client that `_live_scorer` imports, and nothing more."""

    class Interval:
        def resize(self, n):
            return self

    class Variant:
        def __init__(self, **kw):
            self.kw = kw

        @property
        def reference_interval(self):
            return Interval()

    root = types.ModuleType("alphagenome")
    data = types.ModuleType("alphagenome.data")
    genome = types.ModuleType("alphagenome.data.genome")
    models = types.ModuleType("alphagenome.models")
    dna_client = types.ModuleType("alphagenome.models.dna_client")
    scorers = types.ModuleType("alphagenome.models.variant_scorers")
    genome.Variant = Variant
    dna_client.SEQUENCE_LENGTH_1MB = 1_048_576
    scorers.RECOMMENDED_VARIANT_SCORERS = {"RNA_SEQ": object()}
    data.genome = genome
    models.dna_client = dna_client
    models.variant_scorers = scorers
    root.data = data
    root.models = models
    for name, mod in [
        ("alphagenome", root),
        ("alphagenome.data", data),
        ("alphagenome.data.genome", genome),
        ("alphagenome.models", models),
        ("alphagenome.models.dna_client", dna_client),
        ("alphagenome.models.variant_scorers", scorers),
    ]:
        monkeypatch.setitem(sys.modules, name, mod)
    return root


# --- loading the pre-change module -------------------------------------------------------------------


def _frozen_source() -> bytes:
    out = subprocess.run(
        ["git", "cat-file", "blob", FROZEN_BLOB],
        cwd=REPO,
        capture_output=True,
        check=True,
    ).stdout
    got = hashlib.sha256(out).hexdigest()
    assert got == FROZEN_SHA256, f"the pinned pre-change adapter is not the expected bytes: {got}"
    return out


def _load(source: bytes, name: str, tmp_path: Path):
    """Import `source` as its own module in this process, beside the current one."""
    path = tmp_path / f"{name}.py"
    path.write_bytes(source)
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    # registered before execution because `dataclass(slots=True)` resolves its string annotations
    # through sys.modules; the name is unique per load, so the current module is never displaced
    sys.modules[name] = mod
    try:
        spec.loader.exec_module(mod)
    except BaseException:
        sys.modules.pop(name, None)
        raise
    return mod


def _frozen_quantities(module, shared: list) -> dict:
    """Everything AstroREG-2's registered quantity is computed from, for one module and one input."""
    adapter = module.AlphaGenomeAdapter(api_key="not-used-no-request-is-made")
    adapter._client = FakeClient(shared)  # set directly: create_client is never reached
    score = adapter._live_scorer(THRESHOLD)
    effects = score("chr1", 1_000_000, "ACGT", "A")
    rows = et.aggregate(effects)
    return {
        "effects": [[g, t, repr(v)] for g, t, v in effects],
        "aggregate": rows,
        "predict_target": et.predict_target(rows),
        "last_scan": adapter.last_scan,
        "model_tracks": score.model.get("tracks"),
        "model_tracks_sha256": score.model.get("tracks_sha256"),
        "_model": score.model,
    }


def _digest(q: dict) -> str:
    frozen = {k: v for k, v in q.items() if k != "_model"}
    return hashlib.sha256(json.dumps(frozen, sort_keys=True, default=repr).encode()).hexdigest()


# --- (b) the byte-identity proof ---------------------------------------------------------------------


def test_the_frozen_deletion_path_is_byte_identical_with_the_record_present(stub_alphagenome, tmp_path):
    """One shared input, both modules, one process: the effects, their aggregate and the named target are
    the same bytes, and the current module carries the additive record while the pre-change one does not.
    """
    shared = responses()
    old = _load(_frozen_source(), "alphagenome_adapter_frozen", tmp_path)
    before = _frozen_quantities(old, shared)
    after = _frozen_quantities(new, shared)
    assert _digest(after) == _digest(before)
    # and the comparison was not run on an empty chain
    assert before["effects"], "the stub emitted no effects, so nothing was compared"
    assert before["predict_target"] is not None
    # the record is present on the new module only, which is what makes the identity above load-bearing
    assert "gene_axis_outputs" not in before["_model"]
    assert after["_model"]["gene_axis_schema"] == new.GENE_AXIS_SCHEMA
    assert after["_model"]["gene_merge"]["merged_keys"] == 1


def test_the_pre_existing_run_metadata_is_unchanged(stub_alphagenome, tmp_path):
    """Every key the run record carried before still carries the same value; only new keys are added."""
    shared = responses()
    old = _load(_frozen_source(), "alphagenome_adapter_frozen_meta", tmp_path)
    before = _frozen_quantities(old, shared)["_model"]
    after = _frozen_quantities(new, shared)["_model"]
    # `date` is the day the run record was built and is not a product of this change
    for k, v in before.items():
        if k == "date":
            continue
        assert after[k] == v, k
    assert set(after) - set(before) == {"gene_axis_schema", "gene_axis_outputs", "gene_merge"}


def test_the_comparison_catches_a_real_change(stub_alphagenome, tmp_path):
    """The planted counterfactual. The current module with the emitting threshold loosened by one
    character -- `>` to `>=`, which admits the value sitting exactly at the threshold -- must be caught by
    the same comparison that passes above. Without this, a passing identity test would say nothing about
    whether the comparison can see anything at all."""
    shared = responses()
    source = (REPO / "genomeos/predict/alphagenome_adapter.py").read_bytes()
    planted = source.replace(b"if abs(val) > threshold:", b"if abs(val) >= threshold:")
    assert planted != source, "the line the counterfactual perturbs is no longer there"
    old = _load(_frozen_source(), "alphagenome_adapter_frozen_cf", tmp_path)
    bad = _load(planted, "alphagenome_adapter_planted", tmp_path)
    assert _digest(_frozen_quantities(bad, shared)) != _digest(_frozen_quantities(old, shared))


def test_a_record_only_difference_does_not_move_the_frozen_digest(stub_alphagenome, tmp_path):
    """The other half of the counterfactual: a change confined to the additive record must leave the
    frozen digest alone, or the digest would be reporting the record rather than the feature."""
    shared = responses()
    source = (REPO / "genomeos/predict/alphagenome_adapter.py").read_bytes()
    gutted = source.replace(b'rec["cell_values"] = {', b'rec["cell_values_removed"] = {')
    assert gutted != source
    mod = _load(gutted, "alphagenome_adapter_record_only", tmp_path)
    old = _load(_frozen_source(), "alphagenome_adapter_frozen_ro", tmp_path)
    assert _digest(_frozen_quantities(mod, shared)) == _digest(_frozen_quantities(old, shared))


# --- (a) what the record holds -----------------------------------------------------------------------


def test_the_record_keeps_every_gene_axis_column_the_response_carries(stub_alphagenome, tmp_path):
    q = _frozen_quantities(new, responses())
    axis = q["_model"]["gene_axis_outputs"][0]
    assert axis["fields"] == ["gene_id", "strand", "gene_name", "gene_type"]
    row = axis["rows"][0]
    assert row["gene_id"] == "ENSG1.1"
    assert row["strand"] == "+"
    assert row["gene_type"] == "protein_coding"
    assert row["gene_name"] == "AAA"
    # the row's position in the axis, which the effects path does not keep
    assert row["row"] == 0 and row["obs_index"] == "0"
    assert axis["rows"][2]["row"] == 2 and axis["rows"][2]["gene_id"] == "ENSG3.2"


def test_every_recorded_field_is_a_column_of_the_response_and_none_is_invented(stub_alphagenome, tmp_path):
    """The field list is read off the response, not declared: a column the code does not name is still
    kept, and the code names no column the response lacks."""
    q = _frozen_quantities(new, responses())
    axis = q["_model"]["gene_axis_outputs"][0]
    obs_columns = list(responses()[0].obs.columns)
    assert axis["fields"] == [str(c) for c in obs_columns]
    for row in axis["rows"]:
        extra = set(row) - set(axis["fields"])
        assert extra == {"row", "obs_index", "tracks_total", "values_emitted", "cell_values"}


def test_a_cells_whole_multiset_is_kept_including_the_exact_zero(stub_alphagenome, tmp_path):
    """K562 has two tracks in this response. The effects keep one value for it and drop the exact zero;
    the record keeps both, in track order, unthresholded."""
    q = _frozen_quantities(new, responses())
    row = q["_model"]["gene_axis_outputs"][0]["rows"][0]
    assert row["cell_values"]["K562"] == [0.60, 0.00]
    assert row["cell_values"]["HepG2"] == [-0.30]
    assert row["cell_values"]["GM12878"] == [0.05]
    assert "IMR-90" not in row["cell_values"]  # no track for it, so no entry is invented
    # The zero is absent from the effects, which is the loss the record repairs. The three values that
    # are present under one name and one cell are themselves the merge: 0.60 is row 0's K562 track and
    # 0.20 and 0.15 are row 2's, and the effects do not say which row any of them came from.
    aaa_k562 = [e[2] for e in q["effects"] if e[0] == "AAA" and e[1] == "K562"]
    assert aaa_k562 == [repr(0.60), repr(0.20), repr(0.15)]
    assert repr(0.0) not in aaa_k562


def test_a_row_whose_values_are_all_below_the_threshold_is_still_recorded(stub_alphagenome, tmp_path):
    """CCC emits nothing, so the effects do not mention it at all; the record still holds its identity
    and its values, which is the difference between an absent gene and a gene that did not move."""
    q = _frozen_quantities(new, responses())
    rows = q["_model"]["gene_axis_outputs"][0]["rows"]
    ccc = next(r for r in rows if r["gene_name"] == "CCC")
    assert ccc["values_emitted"] == 0
    assert ccc["cell_values"]["K562"] == [0.01, 0.00]
    assert "CCC" not in {e[0] for e in q["effects"]}


def test_values_emitted_reconciles_with_what_the_effects_received(stub_alphagenome, tmp_path):
    """Per row, the number of values clearing the threshold equals the number of effects that row
    contributed, so a pooled track count can be split back across the rows that made it."""
    q = _frozen_quantities(new, responses())
    for axis in q["_model"]["gene_axis_outputs"]:
        if not axis["gene_name_present"]:
            continue
        for row in axis["rows"]:
            expected = sum(1 for v in responses()[axis["output"]].X[row["row"]] if abs(v) > THRESHOLD)
            assert row["values_emitted"] == expected


def test_the_quantiles_layer_is_recorded_as_present_and_not_read(stub_alphagenome, tmp_path):
    q = _frozen_quantities(new, responses())
    axes = q["_model"]["gene_axis_outputs"]
    assert axes[0]["quantiles_present"] is True
    assert axes[1]["quantiles_present"] is False
    assert "quantiles" not in json.dumps(axes[0]["rows"])
    assert "not read" in axes[0]["quantiles_note"]


# --- (c) the merge audit -----------------------------------------------------------------------------


def test_the_audit_names_the_pooled_rows_and_their_gene_ids(stub_alphagenome, tmp_path):
    """What a reader of a new answer can recover: that AAA is two gene-axis rows, which positions they
    held, which gene ids they were, and how many values each contributed to the pooled track count."""
    q = _frozen_quantities(new, responses())
    audit = q["_model"]["gene_merge"]
    assert audit["key"] == "gene_name"
    assert audit["merged_keys"] == 1
    assert audit["merged_rows"] == 2
    merged = audit["merged"][0]
    assert merged["gene_name"] == "AAA"
    assert merged["positions"] == [[0, 0], [0, 2]]
    assert merged["gene_ids"] == ["ENSG1.1", "ENSG3.2"]
    assert merged["values_emitted"] == [2, 5]


def test_the_pooled_value_itself_is_not_repaired(stub_alphagenome, tmp_path):
    """The audit makes the pooling visible and leaves it in place: AAA is still one aggregated row whose
    track count is the sum of both rows and whose largest drop came from the second one unlabelled."""
    q = _frozen_quantities(new, responses())
    aaa = next(r for r in q["aggregate"] if r["gene"] == "AAA")
    assert aaa["n_tracks"] == 7  # 3 from row 0 and 4 from row 2
    assert aaa["max_drop_log2fc"] == -0.9  # row 2's value, with nothing in the row saying so
    assert sum(q["_model"]["gene_merge"]["merged"][0]["values_emitted"]) == aaa["n_tracks"]


def test_rows_carrying_no_symbol_are_counted_rather_than_silently_absent(stub_alphagenome, tmp_path):
    """The second response has no gene_name column, so the effects loop iterates an empty list and the
    response contributes nothing at all. The record says so instead of leaving it to be inferred."""
    q = _frozen_quantities(new, responses())
    second = q["_model"]["gene_axis_outputs"][1]
    assert second["gene_name_present"] is False
    assert len(second["rows"]) == 2
    assert q["_model"]["gene_merge"]["rows_without_key"] == 2
    assert q["_model"]["gene_merge"]["rows"] == 6  # 4 from the first response, 2 from the second


def test_the_audit_reports_no_merge_when_every_symbol_is_distinct():
    axes = [
        {
            "output": 0,
            "rows": [
                {"row": 0, "gene_name": "A", "gene_id": "E1", "values_emitted": 1},
                {"row": 1, "gene_name": "B", "gene_id": "E2", "values_emitted": 2},
            ],
        }
    ]
    audit = new.gene_merge_audit(axes)
    assert audit["merged_keys"] == 0 and audit["merged"] == [] and audit["merged_rows"] == 0
    assert audit["rows"] == 2 and audit["keys"] == 2


# --- the record and the effects name tracks the same way ---------------------------------------------


def test_tissue_names_matches_the_expression_the_effects_loop_uses(stub_alphagenome, tmp_path):
    """`tissue_names` is a factored copy of the effects loop's own rule, so a divergence between them
    would make the record describe a different response than the effects came from. Every branch of the
    rule is present in the stub: a missing gtex tissue, the string "nan", an empty string, and a real
    tissue that overrides the biosample name."""
    for adata in responses():
        names = list(adata.var.get("biosample_name", adata.var.index))
        gtex = list(adata.var.get("gtex_tissue", [None] * len(names)))
        inline = [
            str(g) if g and str(g) not in ("nan", "") else str(n) for g, n in zip(gtex, names, strict=False)
        ]
        assert new.tissue_names(adata) == inline
    assert new.tissue_names(responses()[0]) == ["K562", "K562", "HepG2", "Liver", "GM12878"]


def test_the_recorded_cells_are_the_cells_the_consumer_keeps():
    """Held separately from `enhancer_target.CELLS` so the adapter does not import its consumer; equal,
    so the record covers exactly the cells a closure test can read."""
    assert new.CELL_TRACKS == et.CELLS


def test_the_documented_gene_axis_columns_are_the_pinned_clients_own():
    """Read from alphagenome 0.9.0's `_construct_anndata_from_proto`, which sets gene_id unconditionally
    and the rest only when the gene proto carries the field."""
    assert new.GENE_AXIS_COLUMNS_090 == (
        "gene_id",
        "strand",
        "gene_name",
        "gene_type",
        "junction_Start",
        "junction_End",
    )
    assert "gene_id" in new.GENE_AXIS_COLUMNS_090


def test_a_non_finite_value_is_recorded_as_absent_rather_than_as_a_number():
    assert new._jsonable(float("nan")) is None
    assert new._jsonable(float("inf")) is None
    assert new._jsonable(np.float64(1.5)) == 1.5
    assert new._jsonable(np.int64(3)) == 3
    assert new._jsonable(None) is None


def test_the_whole_record_is_json_writable(stub_alphagenome, tmp_path):
    """The record is written into every cached answer, so it has to survive json.dumps with no NaN."""
    q = _frozen_quantities(new, responses())
    text = json.dumps({k: v for k, v in q["_model"].items()}, allow_nan=False)
    assert json.loads(text)["gene_axis_schema"] == new.GENE_AXIS_SCHEMA


# --- why the registered quantity cannot have moved ---------------------------------------------------
#
# The 20-pair K562 identity check (scripts/astroreg_register.py, `identity_check`) annotates pairs through
# `genomeos.attribution.crispri`. Running it would show one sample of agreement; the import closure shows
# why there is nothing to disagree about, for every pair and not just twenty. These tests hold that
# structural claim so that it fails the moment it stops being true.


def _genomeos_closure(start: Path, root: Path = REPO) -> tuple[set[str], list[str]]:
    """Every `genomeos` module reachable from `start` by following imports, and any import the walker
    could not resolve. Relative imports are resolved against the file's own package rather than skipped:
    `genomeos/ir/__init__.py` uses them, so a walker that ignored them would report an empty result that
    meant "the walk stopped here" while reading as "there is nothing here".
    """
    import ast

    seen: set[Path] = set()
    modules: set[str] = set()
    unresolved: list[str] = []
    stack = [start]
    while stack:
        path = stack.pop()
        if path in seen or not path.exists():
            continue
        seen.add(path)
        try:
            pkg = list(path.resolve().relative_to(root.resolve()).parts[:-1])
        except ValueError:
            pkg = None
        for node in ast.walk(ast.parse(path.read_text())):
            names: list[str] = []
            if isinstance(node, ast.ImportFrom):
                parts: list[str] = []
                if node.level:
                    if pkg is None or node.level - 1 > len(pkg):
                        unresolved.append(f"{path}:{node.lineno}")
                        continue
                    parts = pkg[: len(pkg) - (node.level - 1)]
                if node.module:
                    parts = parts + node.module.split(".")
                base = ".".join(parts)
                if not base:
                    unresolved.append(f"{path}:{node.lineno}")
                    continue
                # the package, and each imported name that is itself a module: `from p import m` is the
                # only spelling of a submodule import, so a walker that took only `p` would stop short
                names = [base] + [
                    f"{base}.{a.name}"
                    for a in node.names
                    if (root / base.replace(".", "/") / f"{a.name}.py").exists()
                    or (root / base.replace(".", "/") / a.name / "__init__.py").exists()
                ]
            elif isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            for mod in names:
                if not mod.startswith("genomeos"):
                    continue
                modules.add(mod)
                stem = mod.replace(".", "/")
                stack.append(root / f"{stem}.py")
                stack.append(root / stem / "__init__.py")
    return modules, unresolved


#: The module that computes the frozen deletion feature: `FEATURES`, `annotate`, `_annotate_one`,
#: `deletion_for`. The 20-pair K562 identity check annotates its pairs through this and nothing else.
FROZEN_FEATURE = "genomeos/attribution/crispri.py"


def test_the_frozen_feature_cannot_reach_the_module_this_lane_changed():
    """The load-bearing claim. Nothing in the module that computes the frozen deletion feature imports
    `genomeos.predict`, so no change in this lane's module is reachable from the registered quantity --
    for every pair, not just the twenty the identity check samples. If this fails, the identity of the
    registered quantity stops following from the code's shape and has to be measured again.

    Scoped to this module deliberately. The scripts that run the check (scripts/astroreg_register.py,
    scripts/astroreg2_endpoint.py) do reach this lane's module in their import graph, through nine
    unrelated files -- attribution/compile.py, direction_v2.py, eqtl.py, executor.py,
    not_open_profile.py, vista.py, genome/regulation.py, predict/__init__.py and
    predict/individual_effects.py. An import graph is not a call graph, so that breadth says nothing
    either way; what the feature computes from is this module's own closure, which is where the claim is
    made and where it can fail.
    """
    modules, unresolved = _genomeos_closure(REPO / FROZEN_FEATURE)
    assert not unresolved, f"imports the walker could not resolve: {unresolved}"
    assert modules, "the walker followed nothing, so its silence means nothing"
    assert not [m for m in modules if m.startswith("genomeos.predict")]


def test_the_closure_walker_finds_a_planted_import(tmp_path):
    """The counterfactual for the four tests above. A file that does reach the changed module must be
    reported, or their silence would be the walker's and not the code's."""
    planted = tmp_path / "planted.py"
    planted.write_text("from genomeos.predict.alphagenome_adapter import AlphaGenomeAdapter\n")
    modules, relative = _genomeos_closure(planted)
    assert not relative
    assert "genomeos.predict.alphagenome_adapter" in modules


def test_the_planted_import_is_found_through_a_second_file(tmp_path):
    """And transitively, since the claim above is about a closure and not about one file's import list.
    A planted tree whose entry file reaches the changed module only through a second file must still be
    reported, or a real indirect dependency could sit behind a passing test."""
    (tmp_path / "genomeos" / "predict").mkdir(parents=True)
    (tmp_path / "genomeos" / "__init__.py").write_text("")
    (tmp_path / "genomeos" / "predict" / "__init__.py").write_text("")
    (tmp_path / "genomeos" / "predict" / "alphagenome_adapter.py").write_text("")
    (tmp_path / "genomeos" / "leaf.py").write_text("from genomeos.predict import alphagenome_adapter\n")
    entry = tmp_path / "entry.py"
    entry.write_text("import genomeos.leaf\n")
    modules, unresolved = _genomeos_closure(entry, root=tmp_path)
    assert not unresolved
    assert "genomeos.leaf" in modules
    assert "genomeos.predict" in modules  # reached only through genomeos/leaf.py
    assert [m for m in modules if m.startswith("genomeos.predict")]


def test_a_relative_import_is_followed_and_not_skipped(tmp_path):
    """The walker's own blind spot, planted: a module reached by `from . import x` must be followed, or
    the four tests above could pass because the walk stopped at a relative import."""
    (tmp_path / "genomeos" / "attribution").mkdir(parents=True)
    (tmp_path / "genomeos" / "__init__.py").write_text("")
    (tmp_path / "genomeos" / "attribution" / "__init__.py").write_text("from . import frozen\n")
    (tmp_path / "genomeos" / "attribution" / "frozen.py").write_text(
        "from ..predict import alphagenome_adapter\n"
    )
    (tmp_path / "genomeos" / "predict").mkdir()
    (tmp_path / "genomeos" / "predict" / "__init__.py").write_text("")
    entry = tmp_path / "genomeos" / "attribution" / "__init__.py"
    modules, unresolved = _genomeos_closure(entry, root=tmp_path)
    assert not unresolved
    assert "genomeos.attribution.frozen" in modules  # found through `from . import frozen`
    assert "genomeos.predict" in modules  # found through `from ..predict import ...`


def test_the_frozen_consumer_reads_the_collapsed_value_and_not_the_new_record():
    """Why fixing the pooling was forbidden, in the consumer's own source: `crispri.deletion_values`
    takes `by_cell[cell]`, the one value per cell name that survived the collapse. It does not read
    `by_cell_summary` and it does not read this lane's record, so neither can move the feature; it also
    means the feature's input is still the collapsed value, which this lane does not change."""
    src = (REPO / "genomeos/attribution/crispri.py").read_text()
    assert '(g.get("by_cell") or {}).get(cell)' in src
    assert "by_cell_summary" not in src
    assert "gene_axis" not in src and "gene_merge" not in src


# --- the number a result can state ------------------------------------------------------------------


def _answer(audit: dict | None) -> dict:
    return {"id": "E1", "model": {} if audit is None else {"gene_merge": audit}}


def test_the_merge_report_counts_the_pooled_rows_across_answers(stub_alphagenome, tmp_path):
    """The number AstroREG-2's result needs: how many gene rows were pooled across the answers it used,
    traceable back to the symbols they were pooled under."""
    q = _frozen_quantities(new, responses())
    answers = [{"id": "E1", "model": q["_model"]}, {"id": "E2", "model": q["_model"]}]
    report = new.merge_report(answers)
    assert report["answers"] == 2
    assert report["answers_with_record"] == 2
    assert report["answers_without_record"] == 0
    assert report["answers_with_a_merge"] == 2
    assert report["merged_symbols"] == 2  # one per answer
    assert report["merged_rows"] == 4  # two rows pooled in each
    assert report["pooled_symbols"] == {"AAA": 2}
    assert report["rows_without_symbol"] == 4  # the second response's two rows, in each answer


def test_an_answer_with_no_record_is_not_counted_as_an_answer_with_no_merge():
    """The distinction the whole report turns on: the 963,406 answers already written carry no axis, so
    they are silent on pooling rather than negative on it. Folding them in as zeros would turn a missing
    record into evidence that nothing merged."""
    report = new.merge_report([_answer(None), _answer(None)])
    assert report["answers"] == 2
    assert report["answers_without_record"] == 2
    assert report["answers_with_record"] == 0
    assert report["answers_with_a_merge"] == 0
    assert report["merged_rows"] == 0
    assert "not evidence" in report["no_record_note"]


def test_an_answer_that_recorded_no_merge_is_distinguished_from_one_with_no_record():
    recorded_none = _answer(new.gene_merge_audit([{"output": 0, "rows": [{"row": 0, "gene_name": "A"}]}]))
    report = new.merge_report([recorded_none, _answer(None)])
    assert report["answers_with_record"] == 1
    assert report["answers_without_record"] == 1
    assert report["answers_with_a_merge"] == 0
    assert report["gene_rows"] == 1  # only the answer that recorded an axis contributes rows


def test_the_merge_report_survives_an_answer_with_no_model_at_all():
    assert new.merge_report([{"id": "E1"}, {}])["answers_without_record"] == 2
