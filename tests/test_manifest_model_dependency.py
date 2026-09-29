# SPDX-License-Identifier: AGPL-3.0-or-later
"""R9 follow-up: what the disk says about the AlphaGenome model behind the all-element sweep, carried by
the headline manifests, and pinned by the scorer from now on. Stubs only: no request, no network."""

import json
from pathlib import Path

import pytest

from genomeos import manifest as mf
from genomeos.predict import alphagenome_adapter as ag
from genomeos.predict.enhancer_target import score_element

RESULTS = Path("data/results")
CARRIERS = (
    "node_containment_audit",
    "node_containment_measured",
    "constrained_unknown_targets",
    "crispri_published",
    "unknown_coverage",
)


# --- negatives first -------------------------------------------------------------------------------


def test_an_unknown_model_has_no_dependency_block():
    with pytest.raises(KeyError):
        mf.model_dependency("some-other-model")


def test_the_sweep_never_claims_a_model_version_it_did_not_request():
    d = mf.model_dependency("alphagenome")
    assert d["model_version"] is None
    assert d["unpinned"].startswith("unpinned: ")
    assert d["track_metadata"]["stored"] is False and d["track_metadata"]["sha256"] is None
    assert d["track_metadata"]["unpinned"].startswith("unpinned: ")
    # the documented default is named as a document, never as the version that answered
    assert "not verified" in d["documented_default"]


def test_validate_names_a_model_dependency_that_neither_pins_nor_says_why():
    base = {k: "n/a: test" for k in ("assembly", "coordinates", "partitions")}
    base |= {
        "sources": [{"accession": "a", "version": "1"}],
        "inputs": [{"path": "p", "sha256": "x", "partition": None}],
    }
    base |= {"code": {"git_sha": "abc", "dirty": False}, "parameters": {}, "exclusions": []}
    assert mf.validate(base) == []  # the block is optional: 955 results never had one
    bad = base | {"model_dependencies": [{"name": "AlphaGenome", "model_version": None}]}
    assert any("model dependency" in p for p in mf.validate(bad))
    assert any("model dependency" in p for p in mf.validate(base | {"model_dependencies": {"name": "x"}}))
    good = base | {"model_dependencies": [mf.model_dependency("alphagenome")]}
    assert mf.validate(good) == []


def test_lock_version_is_none_for_a_package_the_lock_does_not_have(tmp_path):
    lock = tmp_path / "uv.lock"
    lock.write_text(
        '[[package]]\nname = "numpy"\nversion = "2.0.0"\n\n'
        '[[package]]\nname = "alphagenomex"\nversion = "9"\n'
    )
    assert mf.lock_version("alphagenome", lock) is None
    assert mf.lock_version("numpy", lock) == "2.0.0"
    assert mf.lock_version("alphagenome", tmp_path / "missing.lock") is None


def test_a_scorer_without_run_metadata_writes_no_model_record(tmp_path):
    def bare(chrom, pos, ref, alt):
        return [("G", "liver", -0.3)]

    hit = score_element(bare, lambda locus: "ACGTACGTAC", "chr21", "E1", 10, 14, cache=tmp_path)
    assert "model" not in hit
    assert "model" not in json.loads((tmp_path / "chr21" / "E1.json").read_text())


def test_run_metadata_records_an_unrequested_version_as_unrequested():
    class Client:
        _model_version = None

    m = ag.run_metadata(Client(), date="2026-09-28")
    assert m["model_version"] is None and m["model_version_note"].startswith("unrequested")
    assert m["date"] == "2026-09-28"


# --- positives -------------------------------------------------------------------------------------


def test_the_sweep_block_records_what_the_disk_establishes():
    d = mf.model_dependency("alphagenome")
    assert d["client"]["package"] == "alphagenome" and d["client"]["version"] == "0.9.0"
    assert d["api"]["service"] == "google.gdm.gdmscience.alphagenome.v1main.DnaModelService"
    assert d["run_dates"]["first"] == "2026-09-12" and d["run_dates"]["last"] == "2026-09-16"
    obs = d["track_metadata"]["observed"]
    assert obs["tracks_per_element_modal"] == 371 and obs["tissue_names"] == 316
    assert len(obs["tissue_names_sha256"]) == 64
    # the live lock still holds the version the sweep ran with; an upgrade must restate the block
    assert d["client"]["version_in_lock_now"] == mf.lock_version("alphagenome")


def test_the_lock_on_disk_pins_the_client_the_sweep_ran_with():
    assert mf.lock_version("alphagenome") == mf.model_dependency("alphagenome")["client"]["version"]


def test_track_fingerprint_is_order_free_and_counts_modal_tracks():
    answers = [
        {"tracks": 371, "genes": [{"max_drop_tissue": "liver", "max_rise_tissue": "brain"}]},
        {"tracks": 371, "genes": [{"max_drop_tissue": "brain", "max_rise_tissue": None}]},
        {"tracks": 370, "genes": []},
    ]
    a = mf.track_fingerprint(answers)
    b = mf.track_fingerprint(list(reversed(answers)))
    assert a == b
    assert a["elements_scanned"] == 3 and a["tracks_per_element_modal"] == 371 and a["tissue_names"] == 2


def test_the_decorator_adds_the_block_to_whatever_the_writer_returns():
    @mf.depends_on_models("alphagenome")
    def manifest(x):
        return {"parameters": {"x": x}}

    out = manifest(3)
    assert out["parameters"] == {"x": 3} and out["model_dependencies"] == [mf.model_dependency("alphagenome")]
    assert manifest.__name__ == "manifest"


def test_with_model_dependencies_adds_and_never_removes():
    m = {"sources": [1], "parameters": {"x": 1}}
    out = mf.with_model_dependencies(m, "alphagenome")
    assert m == {"sources": [1], "parameters": {"x": 1}}  # the caller's dict is not changed
    assert out["sources"] == [1] and out["parameters"] == {"x": 1}
    assert out["model_dependencies"][0]["name"] == "AlphaGenome"


@pytest.mark.parametrize("name", CARRIERS)
def test_the_headline_results_carry_the_sweep_block(name):
    m = json.loads((RESULTS / f"{name}.json").read_text())[mf.KEY]
    (dep,) = m["model_dependencies"]
    assert dep == mf.model_dependency("alphagenome")


@pytest.mark.parametrize(
    "script",
    [
        "node_containment_audit",
        "constrained_unknown_targets",
        "crispri_published",
        "unknown_coverage",
        "crispri_direction",
    ],
)
def test_the_writers_attach_the_block_through_one_helper(script):
    src = Path(f"scripts/{script}.py").read_text()
    assert '@mf.depends_on_models("alphagenome")' in src
    assert src.index('@mf.depends_on_models("alphagenome")') < src.index("def manifest(")


def test_a_live_scorer_pins_client_version_model_and_tracks_without_a_request():
    pytest.importorskip("alphagenome")
    import numpy as np
    import pandas as pd

    class Scores:
        obs = pd.DataFrame({"gene_name": ["G1", "G2"]})
        var = pd.DataFrame(
            {"biosample_name": ["liver", "K562"], "gtex_tissue": ["Liver", ""]}, index=["t1", "t2"]
        )
        X = np.array([[-0.4, 0.0], [0.01, 0.2]])

    class StubClient:  # dna_client.DnaClient's shape: a requested version name or None, and score_variant
        _model_version = "ALL_FOLDS"
        calls = 0

        def score_variant(self, **kw):
            StubClient.calls += 1
            return [Scores()]

    a = ag.AlphaGenomeAdapter(api_key="stub")
    a._client = StubClient()
    scorer = a._live_scorer(threshold=0.05)
    assert scorer.model["model_version"] == "ALL_FOLDS" and scorer.model["client_version"] == "0.9.0"
    assert "tracks_sha256" not in scorer.model  # nothing has answered yet
    out = scorer("chr21", 100, "ACG", "A")
    assert StubClient.calls == 1 and ("G1", "Liver", -0.4) in out
    assert scorer.model["tracks"] == 2 and len(scorer.model["tracks_sha256"]) == 64
    assert scorer.model["api"]["service"] == "google.gdm.gdmscience.alphagenome.v1main.DnaModelService"


def test_score_element_writes_the_scorers_run_metadata_into_the_cache(tmp_path):
    def scorer(chrom, pos, ref, alt):
        return [("G", "liver", -0.3)]

    scorer.model = {"client_version": "0.9.0", "model_version": None, "date": "2026-09-28", "tracks": 1}
    hit = score_element(scorer, lambda locus: "ACGTACGTAC", "chr21", "E2", 10, 14, cache=tmp_path)
    cached = json.loads((tmp_path / "chr21" / "E2.json").read_text())
    assert cached["model"] == scorer.model and hit["model"] == scorer.model
    scorer.model["tracks"] = 99  # a later answer on the same scorer does not rewrite this one's record
    assert cached["model"]["tracks"] == 1 and hit["model"]["tracks"] == 1
