# SPDX-License-Identifier: AGPL-3.0-or-later
"""The third increment's count: the rules it must not bend, and the audit it must record."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

from genomeos.attribution import cell2
from genomeos.attribution import measured as ms
from genomeos.attribution import targets as tg

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "data/results/response_map_increment3_count.json"


def _load(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def count():
    return _load("_rm3_count", "scripts/response_map_increment3_count.py")


def test_the_reporter_assays_are_the_three_the_measured_layer_holds_beside_crispri(count):
    """The kinds are not a list this script chose: they are the layer's own four, less the one."""
    assert set(ms.ASSAYS) - {"crispri"} == {"lentimpra", "vista", "satmut"}
    assert "crispri" in ms.OUTCOME_KIND


def test_the_eqtl_margin_is_the_one_the_retained_set_was_distilled_under(count):
    """A margin chosen here would make the attachment rule this lane's, not the record's."""
    eqtl_targets = _load("_eqtl_targets_for_test", "scripts/eqtl_targets.py")
    assert count.EQTL_MARGIN == eqtl_targets.MARGIN


def test_hits_in_is_the_same_window_the_retained_set_was_kept_under(count):
    """`eqtl.Intervals.at` keeps a 1-based variant at `s < pos <= e` over a half-open interval.

    `hits_in` must apply exactly that window around the element widened by the margin, or this count
    would admit or drop a variant the distillation did not, under a rule of this lane's own making.
    """
    from genomeos.attribution import eqtl

    start, end = 1000, 1000
    lo, hi = start - count.EQTL_MARGIN, end + count.EQTL_MARGIN
    iv = eqtl.Intervals()
    iv.add("chr1", lo, hi, "E")
    iv.freeze()
    positions = (lo - 1, lo, lo + 1, end, hi, hi + 1)
    rows = [{"pos": p} for p in positions]
    got = {r["pos"] for r in count.hits_in(rows, start, end)}
    kept = {p for p in positions if iv.at("chr1", p)}
    assert got == kept
    assert got == {lo + 1, end, hi}  # lo itself is outside: s < pos, not s <= pos


def test_covered_by_finds_an_overlap_and_refuses_a_gap(count):
    index = [(100, 200), (5000, 5100)]
    assert count.covered_by(index, 150, 160)
    assert count.covered_by(index, 50, 150)
    assert not count.covered_by(index, 300, 400)
    assert not count.covered_by([], 1, 2)


def test_the_cache_audit_records_an_open_it_was_never_told_about(count, tmp_path):
    """The point of the audit: the open is recorded without the caller declaring it."""
    root = tmp_path / "elements"
    root.mkdir()
    f = root / "chr21.json"
    f.write_text("{}")
    audit = count._CacheAudit(root)
    with audit:
        with open(f) as fh:  # noqa: PTH123 - the patched builtin is what is under test
            fh.read()
        with open(tmp_path / "elsewhere.txt", "w") as fh:
            fh.write("x")
    assert audit.block()["opens"] == 1
    assert f.resolve().as_posix() in audit.block()["files"]
    assert open is audit._real  # restored


def test_the_cache_audit_root_is_the_real_cache(count):
    assert count._CacheAudit().root == Path(tg.ELEMENT_CACHE).resolve()


def test_the_no_measured_gene_token_suppresses_only_the_gene_arm(count):
    """Two reporter keys 10 Mb apart must not be one locus, and two 10 kb apart must be."""
    far = [
        cell2.LocusKey("lentimpra", "chr1", 1_000, 1_200, "lentimpra|A"),
        cell2.LocusKey("lentimpra", "chr1", 20_000_000, 20_000_200, "lentimpra|B"),
    ]
    near = [
        cell2.LocusKey("lentimpra", "chr1", 1_000, 1_200, "lentimpra|A"),
        cell2.LocusKey("lentimpra", "chr1", 11_000, 11_200, "lentimpra|B"),
    ]
    assert cell2.count_loci(far) == 2
    assert cell2.count_loci(near) == 1


def test_reporter_outcomes_splits_by_the_cell_the_chain_was_measured_in(count):
    row = {
        "assays": ["crispri", "lentimpra"],
        "measured": {
            "lentimpra": {
                "label_by_cell": {"K562": "silent", "HepG2": "active"},
                "cells_active": ["HepG2"],
                "cells_conflicting": [],
            }
        },
    }
    got = count.reporter_outcomes({"E1"}, {"E1": row}, {"E1": {"K562"}}, "lentimpra")
    assert got["by_element"] == {"at_least_one_cell_active": 1}
    assert got["cell_labels"] == {"silent": 1, "active": 1}
    # the chain was measured in K562, where the reporter called the sequence silent
    assert got["in_the_same_cell_as_the_chains_crispri_observation"] == {"silent": 1}
    assert got["elements_whose_reporter_cells_include_none_of_the_chains_cells"] == 0


def test_reporter_outcomes_counts_an_element_with_no_label_under_its_own_name(count):
    row = {"assays": ["lentimpra"], "measured": {"lentimpra": {}}}
    got = count.reporter_outcomes({"E1"}, {"E1": row}, {}, "lentimpra")
    assert got["by_element"] == {"no_cell_label_in_the_record": 1}
    assert got["elements_whose_reporter_cells_include_none_of_the_chains_cells"] == 1


def test_eqtl_outcomes_never_calls_a_shared_gene_agreement(count):
    hits = [{"gene_id": "ENSG1", "tissue": "Whole_Blood", "slope": -0.3, "pval": 1e-9}]
    got = count.eqtl_outcomes({"E1": hits}, {"E1": {"HBE1"}}, {"ENSG1": {"HBE1"}})
    assert got["elements_where_a_measured_gene_is_also_a_retained_egene"] == 1
    assert got["per_element"][0]["of_those_also_a_retained_egene"] == ["HBE1"]
    assert "not agreement" in got["that_is_not_agreement"]


def test_eqtl_outcomes_cannot_match_an_unresolved_ensembl_id(count):
    hits = [{"gene_id": "ENSG9", "tissue": "Lung", "slope": 0.1, "pval": 1e-8}]
    got = count.eqtl_outcomes({"E1": hits}, {"E1": {"HBE1"}}, {})
    assert got["elements_where_a_measured_gene_is_also_a_retained_egene"] == 0


def test_increment_2s_limit_is_carried_word_for_word(count):
    assert "23.26% of one of the measured layer's four assays" in count.INCREMENT_2_LIMIT
    assert (
        "they read a sequence, or the bases inside it, not the native locus after a perturbation of it"
        in count.INCREMENT_2_LIMIT
    )


@pytest.fixture(scope="module")
def payload():
    if not RESULT.exists():
        pytest.skip("the genome-wide count has not been run here")
    with RESULT.open() as fh:
        return json.load(fh)


@pytest.mark.skipif(not RESULT.exists(), reason="the genome-wide count has not been run here")
class TestTheCommittedCount:
    def test_it_reproduces_the_110_increment_2_committed(self, payload):
        assert payload["reconciliation"]["holds"] is True
        assert payload["reconciliation"]["addable_loci"] == 110

    def test_no_reporter_assay_contributes_a_native_locus_observation(self, payload):
        for k in payload["kinds_of_evidence"]:
            if k["kind"] in ("lentimpra", "vista", "satmut"):
                assert k["an_observation_of_the_native_locus"] is False
                assert k["native_locus_observations_on_the_existing_110"] == 0
                assert k["native_locus_loci_on_the_existing_110"] == 0

    def test_every_kind_states_what_its_observation_is(self, payload):
        for k in payload["kinds_of_evidence"]:
            assert k["what_the_observation_is"]
            assert k["source"]

    def test_the_eqtl_count_names_its_own_frame(self, payload):
        k = next(k for k in payload["kinds_of_evidence"] if k["kind"] == "eqtl")
        assert (
            k["shown_elements_with_a_retained_hit"]
            <= (k["shown_elements_the_retained_set_could_have_seen_at_all"])
        )
        assert "UNASSESSED" in k["the_frame_is_the_limit"]

    def test_the_cross_tabulation_partitions_the_110(self, payload):
        b = payload["kinds_per_locus"]
        assert sum(b["beyond_those_two"].values()) == 110
        assert b["loci_carrying_at_least_one_further_kind"] + b["loci_carrying_nothing_further"] == 110

    def test_no_model_request_and_no_per_element_cache_open(self, payload):
        assert payload["alphagenome_requests"] == 0
        assert payload["per_element_response_cache"]["opens"] == 0

    def test_it_says_what_it_cannot_establish(self, payload):
        assert len(payload["cannot_establish"]) >= 5
