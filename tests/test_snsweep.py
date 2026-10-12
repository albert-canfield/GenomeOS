# SPDX-License-Identifier: AGPL-3.0-or-later
"""`genomeos/attribution/snsweep.py`: the adopted rule is lane-pseudobulk's, unchanged.

THE WHOLE PROTECTION IN THIS LANE IS THAT IT DID NOT MOVE THE RULE IT WAS GIVEN, AND THAT ITS ONE
CHANGE -- the modality -- ONLY EVER NARROWS. lane-pbrun's agreement rests on four cell types of
which one is single-NUCLEUS, and this lane decides whether three more can be found beside it. A
lane that could reach four by relaxing a clause, by counting a single-CELL experiment, or by
letting a near-spelling into the headline would have answered nothing. So these tests pin the
numbers, the imported objects and every clause that could have been relaxed, and they plant the
cases that would catch a relaxation:

  - the minimum is lane-pseudobulk's OBJECT, not a copy of its numbers, so a drift is impossible;
  - the two cross-archive readings are lane-scsweep's OBJECTS, for the same reason;
  - a single-CELL experiment is NOT single-nucleus in any of the three archives' vocabularies,
    which is where the three extra cell types would have had to come from;
  - the ATAC arm of a multiome is not an RNA library, though its sample says nuclei;
  - a label that is not the identical string counts in the wider count and NOT in the headline,
    which is the whole of the SK-N-SH case;
  - this lane's verdict never returns the adopted rule's `QUALIFIES`, because it does not count
    the 1,000 genes;
  - the committed reading is RECOMPUTED by the decision function, so the count and the verdict in
    the record cannot disagree with each other or change silently.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import tracked_paths as tp

from genomeos.attribution import pseudobulk as pb
from genomeos.attribution import scsweep as sc
from genomeos.attribution import snsweep as sn


def _sweep():
    """`scripts/snsweep_sweep.py` loaded as a module. It is a script, as its siblings are."""
    spec = importlib.util.spec_from_file_location(
        "snsweep_sweep", Path(__file__).resolve().parents[1] / "scripts" / "snsweep_sweep.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_minimum_is_the_predecessor_s_own_object_and_not_a_copy():
    """IMPORTED, not copied: a drift between the lanes is impossible by construction."""
    assert sn.MINIMUM is pb.MINIMUM
    assert sn.MINIMUM["cell_types"] == 4
    assert sn.MINIMUM["genes"] == 1000
    assert sn.ADOPTED_RULE_VERSION == pb.SELECTION_RULE["version"] == "pseudobulk-1"


def test_the_two_cross_archive_readings_are_the_predecessor_s_own_objects():
    """A near-spelling gets into the headline only by rewriting lane-scsweep's policy, not this."""
    assert sn.BIOSAMPLE_IDENTITY_POLICY is sc.BIOSAMPLE_IDENTITY_POLICY
    assert sn.TERM_MATCH_POLICY is sc.TERM_MATCH_POLICY
    assert sn.WORDING_LIMIT is pb.WORDING_LIMIT


def test_the_modality_counted_is_the_single_NUCLEUS_one_and_the_others_are_named():
    """The narrowing, stated as objects: the adopted rule's three titles split one from two."""
    assert sn.SINGLE_NUCLEUS_ASSAY == "snRNA-seq"
    assert sn.SINGLE_CELL_ASSAYS_NOT_COUNTED == ("scRNA-seq", "long read scRNA-seq")
    assert set(pb.SINGLE_CELL_ASSAYS) == {sn.SINGLE_NUCLEUS_ASSAY} | set(sn.SINGLE_CELL_ASSAYS_NOT_COUNTED)


def test_PLANTED_a_single_CELL_experiment_is_not_counted_in_any_archive_s_vocabulary():
    """The relaxation that would have reached 4, refused in all three vocabularies.

    ENCODE's GM12878 scRNA-seq, CELLxGENE's Tabula Sapiens testis and the GEO K562 h5ad are the
    three cell types lane-pbrun read beside ovary. Each is single-CELL, and each is refused here.
    """
    assert sn.is_single_nucleus_encode("snRNA-seq") is True
    for title in ("scRNA-seq", "long read scRNA-seq", "polyA plus RNA-seq", ""):
        assert sn.is_single_nucleus_encode(title) is False
    assert sn.is_single_nucleus_cellxgene(["nucleus"]) is True
    assert sn.is_single_nucleus_cellxgene(["cell"]) is False
    assert sn.is_single_nucleus_cellxgene([]) is False
    assert sn.is_single_nucleus_geo("RNA-Seq", "Fresh_K562 | droplet single-cell RNA-seq of cells") is False


def test_PLANTED_the_ATAC_arm_of_a_multiome_is_not_an_RNA_library():
    """GSM8159280 and GSM8159281 are the same nuclei; only one of them is RNA."""
    nuclei = "SK.N.SH | Single nuclei | 10X Genomics Cell Ranger pipeline"
    assert sn.is_single_nucleus_geo("RNA-Seq", nuclei) is True
    assert sn.is_single_nucleus_geo("ATAC-seq", nuclei) is False


def test_PLANTED_the_bare_word_nuclei_is_not_a_modality():
    """`nuclei` matches nuclear extraction and nuclear RNA in bulk assays, so it is not enough."""
    assert sn.is_single_nucleus_geo("RNA-Seq", "nuclei were isolated for nuclear RNA-seq") is False
    assert sn.is_single_nucleus_geo("RNA-Seq", "single nuclei RNA-seq") is True
    assert sn.is_single_nucleus_geo("RNA-Seq", "whole cell and single nucleus libraries") is False


def test_below_four_is_the_predecessor_s_own_complete_answer():
    for n in (0, 1, 2, 3):
        assert sn.verdict(n) == "NO_SOURCE_QUALIFIES"


def test_at_four_the_cell_type_half_is_reached_and_QUALIFIES_is_never_claimed():
    """The adopted `QUALIFIES` needs BOTH numbers, and this lane counts no genes."""
    for n in (4, 5, 13):
        assert sn.verdict(n) == "CELL_TYPE_MINIMUM_REACHED"
        assert sn.verdict(n) != "QUALIFIES"
    assert pb.verdict(13, None) == "NO_SOURCE_QUALIFIES"


def test_the_covered_set_is_the_thirteen_and_is_not_widened():
    terms = sn.covered_terms()
    assert len(terms) == 13
    assert terms["ovary"] == "ovary"
    assert terms["SK_N_SH"] == "SK-N-SH"
    assert set(terms) == {c["cell"] for c in pb.covered_cells()}


def test_the_predecessor_count_is_split_by_modality_from_the_committed_records():
    """Where GM12878, K562 and testis leave the count: they are single-CELL, nothing else."""
    tp.must_be_present(
        "data/results/pseudobulk_inventory.json",
        "data/results/scsweep_inventory.json",
        was="the predecessor lanes' inventories are committed records of this repository",
    )
    single_nucleus = sn.predecessor_single_nucleus_cells()
    single_cell_only = sn.predecessor_single_cell_only_cells()
    assert single_nucleus == ["ovary"]
    assert "GM12878" in single_cell_only
    assert not set(single_nucleus) & set(single_cell_only)


def test_every_named_candidate_carries_its_modality_evidence_and_its_query():
    for candidate in sn.CANDIDATES:
        assert candidate["modality_evidence"].strip()
        assert candidate["found_by"].strip()
        assert candidate["archive"] in ("ENCODE", "CELLxGENE", "GEO")
        assert candidate["cell"] in sn.covered_terms()


def test_every_rejection_names_the_clause_it_fails():
    clauses = {r["clause_failed"] for r in sn.REJECTED}
    assert any("biosample identity" in c for c in clauses)
    assert any("not treated" in c for c in clauses)
    assert any("no single-nucleus experiment" in c for c in clauses)
    for rejection in sn.REJECTED:
        assert rejection["why"].strip()
        assert rejection["accession"].strip()


def test_the_entry_script_declares_relative_input_paths():
    """Rule (d): `input_entry` records what it is given, so a worktree-absolute path is unrebuildable."""
    sweep = _sweep()
    for path in sweep.DECLARED_INPUTS:
        assert not path.startswith("/"), path
        assert path.startswith("data/"), path
        assert (Path(__file__).resolve().parents[1] / path).exists(), path


def test_the_geo_query_terms_do_not_include_the_bare_word_nuclei():
    """The query that would have inflated every count, absent from the committed script."""
    sweep = _sweep()
    joined = " ".join(sweep.GEO_SN_TERMS).lower()
    assert "snrna" in joined
    assert "single-nucleus" in joined
    assert '"nuclei"' not in joined


def test_the_committed_result_says_what_it_counted_and_declares_relative_inputs():
    """The result as committed. `must_be_present` RAISES by name instead of skipping.

    A skip keyed on this file's absence would go inert the moment the file was committed, which is
    the defect `tests/test_inert_skip_guards.py` exists for (`d427f7d`).
    """
    tp.must_be_present(
        "data/results/snsweep_inventory.json", was="the sweep has not been run on this machine"
    )
    record = json.loads(sn.INVENTORY_PATH.read_text())
    assert record["answer"] in ("NO_SOURCE_QUALIFIES", "CELL_TYPE_MINIMUM_REACHED")
    # The verdict is RECOMPUTED, never reprinted: the count and the answer cannot disagree.
    assert record["answer"] == sn.verdict(record["exact_label_count"])
    assert record["exact_label_count"] == len(record["exact_label_cells"])
    assert record["same_biosample_count"] == len(record["same_biosample_cells"])
    assert set(record["exact_label_cells"]) <= set(record["same_biosample_cells"])
    assert record["minimum_cell_types"] == 4 == sn.MINIMUM["cell_types"]
    assert record["rule_adopted_verbatim"]["commit"] == "2b9cd0a"
    assert record["rule_adopted_verbatim"]["minimum"] == sn.MINIMUM
    assert record["alphagenome_requests"] == 0
    assert record["covered_cell_count"] == 13
    # Rule (d): every declared input path is RELATIVE to the repository.
    for entry in record["result_manifest"]["inputs"]:
        assert not entry["path"].startswith("/"), entry["path"]
        assert entry["path"].startswith("data/"), entry["path"]
    # No file was fetched: every candidate file says so of itself.
    for candidate in record["candidates"]:
        for f in candidate["files"]:
            assert f["fetched"] is False


def test_the_committed_count_is_BELOW_the_minimum_and_the_lane_closes():
    """THE READING ITSELF, pinned: 1 of 13 on the headline count, 2 on the wider one.

    If a later run finds a fourth single-nucleus cell type this test FAILS LOUDLY, which is what
    it is for: the answer Albert is given rests on this number, and a number that can move without
    a test moving is a number no one checked.
    """
    tp.must_be_present(
        "data/results/snsweep_inventory.json", was="the sweep has not been run on this machine"
    )
    record = json.loads(sn.INVENTORY_PATH.read_text())
    assert record["exact_label_count"] == 1
    assert record["exact_label_cells"] == ["ovary"]
    assert record["same_biosample_count"] == 2
    assert record["same_biosample_cells"] == ["SK_N_SH", "ovary"]
    assert record["exact_label_count"] < sn.MINIMUM["cell_types"]
    assert record["same_biosample_count"] < sn.MINIMUM["cell_types"]
    assert record["answer"] == "NO_SOURCE_QUALIFIES"
    assert record["answer_means"] == sn.CLOSES_WITHOUT_A_LIFT
    assert "SINGLE-CELL PSEUDOBULK ONLY" in sn.CLOSES_WITHOUT_A_LIFT
    assert "UNINTERPRETABLE" in sn.CLOSES_WITHOUT_A_LIFT
    assert "CLOSES WITHOUT A LIFT" in sn.CLOSES_WITHOUT_A_LIFT


def test_the_one_GEO_candidate_counts_in_the_WIDER_count_and_not_in_the_headline():
    """SK-N-SH is the whole case for the imported TERM_MATCH_POLICY, so it is pinned here.

    `SK.N.SH` is the same biosample as `SK-N-SH` and is not the identical string. Admitting it to
    the headline would be the one relaxation that could have moved this lane's number, and it
    would still not have reached 4.
    """
    tp.must_be_present(
        "data/results/snsweep_inventory.json", was="the sweep has not been run on this machine"
    )
    record = json.loads(sn.INVENTORY_PATH.read_text())
    geo = [c for c in record["candidates"] if c["archive"] == "GEO"]
    assert len(geo) == 1
    candidate = geo[0]
    assert candidate["cell"] == "SK_N_SH"
    assert candidate["label_in_archive"] == "SK.N.SH"
    assert candidate["label_is_identical_string"] is False
    assert candidate["admitted_on_every_clause_but_the_label"] is True
    assert "SK_N_SH" not in record["exact_label_cells"]
    assert "SK_N_SH" in record["same_biosample_cells"]
    # Its modality came from the sample's own words, and its matrix was never downloaded.
    modality = next(c for c in candidate["clauses"] if c["clause"] == "SINGLE-NUCLEUS modality")
    assert modality["satisfied"] is True
    assert all(f["fetched"] is False for f in candidate["files"])


def test_ENCODE_has_no_single_nucleus_experiment_for_twelve_of_the_thirteen():
    """The negative that decides the lane, and it is a NEGATIVE and not an unsearched gap."""
    tp.must_be_present(
        "data/results/snsweep_inventory.json", was="the sweep has not been run on this machine"
    )
    record = json.loads(sn.INVENTORY_PATH.read_text())
    per_assay = record["encode_covered_terms_per_assay"]
    assert per_assay[sn.SINGLE_NUCLEUS_ASSAY] == ["ovary"]
    # The facet was ANSWERED for every assay title: a zero here is a count and not a failure.
    for facet in record["encode_assay_facets"]:
        assert facet["http_status"] in (200, 404)
        if facet["http_status"] == 404:
            assert facet["total"] == 0
    with_sn = [c for c in record["per_cell"] if c["encode_snRNA_seq_experiments"]]
    assert [c["cell"] for c in with_sn] == ["ovary"]


def test_a_rate_limited_search_is_NOT_recorded_as_a_count_of_zero():
    """lane-scsweep's first run carried four HTTP 429s into its record as null counts."""
    tp.must_be_present(
        "data/results/snsweep_inventory.json", was="the sweep has not been run on this machine"
    )
    record = json.loads(sn.INVENTORY_PATH.read_text())
    assert len(record["geo_search_counts"]) == 26  # thirteen cell types, samples and series
    for search in record["geo_search_counts"]:
        if search["record_count"] is None:
            assert search["http_status"] != 200
            assert "NOT a count of zero" in search["null_count_means"]
        else:
            assert search["http_status"] == 200
