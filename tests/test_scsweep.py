# SPDX-License-Identifier: AGPL-3.0-or-later
"""`genomeos/attribution/scsweep.py`: the adopted rule is lane-pseudobulk's, unchanged.

THE WHOLE PROTECTION IN THIS LANE IS THAT IT DID NOT MOVE THE RULE IT WAS GIVEN. Its predecessor
answered `NO_SOURCE_QUALIFIES` at 2 of 13 against a minimum of 4, and this lane searched the three
archives that predecessor named as unsearched. A lane that can change the count by relaxing a
clause has answered nothing, so these tests pin the clauses that could have been relaxed and the
two numbers, and they plant the cases that would catch a relaxation:

  - the minimum is lane-pseudobulk's OBJECT, not a copy of its numbers, so a drift is impossible;
  - a Perturb-seq or CRISPRi non-targeting control is named MODIFIED, which is where the lift was
    expected to come from and where under this rule it cannot;
  - a covered cell type that appears only as a `cell_type` annotation inside a tissue dataset is
    NOT admitted, which is the cell-composition decision the adopted rule does not license;
  - a label that is not the identical string counts in the wider count and NOT in the headline;
  - this lane's verdict never returns the adopted rule's `QUALIFIES`, because it does not count
    the 1,000 genes.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import tracked_paths as tp

from genomeos.attribution import pseudobulk as pb
from genomeos.attribution import scsweep as sc


def _sweep():
    """`scripts/scsweep_sweep.py` loaded as a module. It is a script, as its siblings are."""
    spec = importlib.util.spec_from_file_location(
        "scsweep_sweep", Path(__file__).resolve().parents[1] / "scripts" / "scsweep_sweep.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_minimum_is_the_predecessor_s_own_object_and_not_a_copy():
    """IMPORTED, not copied: a drift between the two lanes is impossible by construction."""
    assert sc.MINIMUM is pb.MINIMUM
    assert sc.MINIMUM["cell_types"] == 4
    assert sc.MINIMUM["genes"] == 1000
    assert sc.ADOPTED_RULE_VERSION == pb.SELECTION_RULE["version"] == "pseudobulk-1"


def test_below_four_is_the_predecessor_s_own_complete_answer():
    for n in (0, 1, 2, 3):
        assert sc.verdict(n) == "NO_SOURCE_QUALIFIES"


def test_at_four_the_cell_type_half_is_reached_and_QUALIFIES_is_never_claimed():
    """The adopted `QUALIFIES` needs BOTH numbers, and this lane counts no genes."""
    for n in (4, 5, 13):
        assert sc.verdict(n) == "CELL_TYPE_MINIMUM_REACHED"
        assert sc.verdict(n) != "QUALIFIES"
    # The predecessor's own function, on a gene count this lane cannot supply.
    assert pb.verdict(4, None) == "NO_SOURCE_QUALIFIES"


def test_the_covered_set_is_thirteen_and_is_not_widened():
    terms = sc.covered_terms()
    assert len(terms) == 13
    assert terms == {c["cell"]: c[pb.MATCH_FIELD] for c in pb.covered_cells()}
    for candidate in sc.CANDIDATES:
        assert candidate["cell"] in terms, candidate["cell"]


def test_every_candidate_names_the_query_that_found_it_and_its_archive_label():
    for candidate in sc.CANDIDATES:
        assert candidate["found_by"].strip()
        assert candidate["label_in_archive"].strip()
        assert candidate["archive"] in sc.ARCHIVES


def test_PLANTED_a_perturb_seq_non_targeting_control_is_named_MODIFIED():
    """Where the lift was expected to come from, and where this rule does not let it.

    A non-targeting control cell carries a guide vector and a Cas9 effector. The adopted clause
    excludes a genetically modified biosample in the same words it excludes CRISPRi RNA-seq, so an
    edit that admitted these would have moved the rule, and it fails here by name.
    """
    sweep = _sweep()
    planted = (
        "cell line: K562; genotype: Zim3-dCas9 expressing; "
        "treatment: non-targeting control sgRNA, lentiviral transduction"
    ).lower()
    assert [w for w in sweep.MODIFIED_WORDS if w in planted], (
        "a non-targeting control must be read as a genetically modified biosample"
    )
    for required in ("cas9", "sgrna", "lentivir"):
        assert required in sweep.MODIFIED_WORDS
    # The negative control of the same instrument: an unmodified culture is not named modified.
    clean = "cell line: k562, l1210; cell type: lymphoblast; tissue: human and mouse blood"
    assert [w for w in sweep.MODIFIED_WORDS if w in clean] == []


def test_a_pfa_fixed_preparation_is_named_TREATED():
    """The admitted K562 file is the FRESH one, and the instrument must see why."""
    sweep = _sweep()
    assert [w for w in sweep.TREATED_WORDS if w in "pfa fixed cells"], (
        "fixation is a treatment of the cells and must be read as one"
    )
    assert "genotype" not in " ".join(sweep.MODIFIED_WORDS)


def test_the_rejections_name_the_two_the_supervisor_expected():
    """Perturb-seq controls and HepG2, each with the clause it fails or the queries that missed it."""
    blob = json.dumps(sc.REJECTED)
    assert "Perturb-seq" in blob
    assert "HepG2" in blob
    perturb = [r for r in sc.REJECTED if "Perturb-seq" in r["archive"]]
    assert perturb and perturb[0]["clause_failed"] == "not genetically modified"
    hepg2 = [r for r in sc.REJECTED if r["cell"] == "HepG2"]
    assert hepg2 and "these queries" in hepg2[0]["why"]


def test_the_cell_type_inside_a_tissue_is_NOT_admitted():
    """The composition decision the adopted rule does not license, kept shut.

    CELLxGENE carries `astrocyte`, `hepatocyte`, `keratinocyte`, `cardiac muscle cell` and
    `CD14-positive monocyte` as CL annotations inside tissue datasets. Admitting one would widen
    the rule, so the record must name that exclusion and no candidate may rest on it.
    """
    rejected = [r for r in sc.REJECTED if r["clause_failed"] == "biosample identity"]
    assert rejected, "the exclusion must be recorded, not merely observed"
    assert "cell-composition decision" in rejected[0]["why"]
    for cell in ("astrocyte", "hepatocyte", "keratinocyte", "cardiac_muscle_cell"):
        assert cell not in [c["cell"] for c in sc.CANDIDATES]


def test_a_label_that_is_not_the_identical_string_is_not_in_the_headline_count():
    """Both counts are reported and the headline is the stricter one."""
    assert "exact_label_count" in sc.TERM_MATCH_POLICY
    assert "same_biosample_count" in sc.TERM_MATCH_POLICY
    assert "STRICTER" in sc.TERM_MATCH_POLICY
    terms = sc.covered_terms()
    differing = [c for c in sc.CANDIDATES if c["label_in_archive"] != terms[c["cell"]]]
    assert differing, "the fork is real in this lane's own candidates, not hypothetical"
    for candidate in differing:
        assert "same_biosample_count" in candidate["notes"]


def test_this_lane_adds_no_biosample_and_reads_no_value():
    assert "adds NO biosample" in sc.__doc__
    assert "READS NO RNA VALUE" in sc.__doc__
    assert "stop list of 2026-09-27" in sc.__doc__


def test_the_declared_inputs_are_relative_to_the_repository():
    """Rule (d): `input_entry` records what it is given verbatim (`d72f959`)."""
    sweep = _sweep()
    for path in sweep.DECLARED_INPUTS:
        assert not path.startswith("/"), path
        assert path.startswith("data/"), path
        assert (Path(__file__).resolve().parents[1] / path).exists(), path


def test_the_committed_result_says_what_it_counted_and_declares_relative_inputs():
    """The result as committed. `must_be_present` RAISES by name instead of skipping.

    A skip keyed on this file's absence would go inert the moment the file was committed, which is
    the defect `tests/test_inert_skip_guards.py` exists for (`d427f7d`): eight guards went inert
    that way, not because anyone edited them but because the lanes recorded their measurements.
    """
    tp.must_be_present(
        "data/results/scsweep_inventory.json", was="the sweep has not been run on this machine"
    )
    record = json.loads(sc.INVENTORY_PATH.read_text())
    assert record["answer"] in ("NO_SOURCE_QUALIFIES", "CELL_TYPE_MINIMUM_REACHED")
    assert record["answer"] == sc.verdict(record["exact_label_count"])
    assert record["exact_label_count"] == len(record["exact_label_cells"])
    assert record["same_biosample_count"] == len(record["same_biosample_cells"])
    assert set(record["exact_label_cells"]) <= set(record["same_biosample_cells"])
    assert set(record["already_qualifying_from_encode"]) <= set(record["exact_label_cells"])
    assert record["minimum_cell_types"] == 4 == sc.MINIMUM["cell_types"]
    assert record["rule_adopted_verbatim"]["commit"] == "2b9cd0a"
    assert record["alphagenome_requests"] == 0
    # Rule (d): every declared input path is RELATIVE to the repository, never worktree-absolute.
    for entry in record["result_manifest"]["inputs"]:
        assert not entry["path"].startswith("/"), entry["path"]
        assert entry["path"].startswith("data/"), entry["path"]
    # No file was fetched: every candidate file says so of itself.
    for candidate in record["candidates"]:
        for f in candidate["files"]:
            assert f["fetched"] is False


def test_a_rate_limited_search_is_NOT_recorded_as_a_count_of_zero():
    """The first run's defect, pinned: four of thirteen NCBI searches answered 429.

    A null `sample_count` carried a rate limit into the record where a reader could take it for a
    negative about that cell type. Every search in the committed result is either answered or says
    in its own words that it was not.
    """
    tp.must_be_present(
        "data/results/scsweep_inventory.json", was="the sweep has not been run on this machine"
    )
    record = json.loads(sc.INVENTORY_PATH.read_text())
    for search in record["geo_search_counts"]:
        if search["sample_count"] is None:
            assert search["http_status"] != 200
            assert "NOT a count of zero" in search["null_count_means"]
        else:
            assert search["http_status"] == 200
