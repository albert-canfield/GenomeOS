# SPDX-License-Identifier: AGPL-3.0-or-later
"""The selection rule of lane-oodinventory, pinned so a later search cannot bend it.

The lane's whole value is that the rule was fixed before any public source was looked at. These
tests hold the two parts that a search would be tempted to loosen: what counts as absent (a
child of a covered parent is NOT absent, and nor is a parent of a covered child), and the
refusal to call a candidate absent whose ontology relation was never checked.
"""

from __future__ import annotations

import json
from pathlib import Path

from genomeos.attribution import oodinventory as ood

COVERED = {"UBERON:0002048", "CL:0000084", "EFO:0002067"}


def _candidate(**kw) -> ood.Candidate:
    base = {
        "curie": "CL:0002063",
        "label": "type II pneumocyte",
        "source": "test",
        "accessibility": "ENCSR000AAA",
        "rna": "ENCSR000BBB",
        "licence": "open",
        "relation_checked": True,
    }
    base.update(kw)
    return ood.Candidate(**base)


def test_covered_set_is_the_rna_rows_of_the_metadata_we_hold():
    rows = ood.rna_rows()
    assert rows, "the metadata copy has rna_seq rows"
    assert all(r["output"] == "rna_seq" for r in rows)
    curies = ood.covered_curies()
    # Every RNA row carries a CURIE, which is why the rule can match on it at all.
    assert all(r["ontology_curie"].strip() for r in rows)
    assert curies == {r["ontology_curie"].strip() for r in rows}


def test_no_rna_row_carries_an_encode_accession():
    # lane-trainexp's finding, re-measured here: accession-level matching is impossible.
    assert not [r for r in ood.rna_rows() if r["name"].startswith("ENC")]


def test_a_curie_in_the_metadata_is_not_absent():
    row = ood.verdict(_candidate(curie="CL:0000084"), COVERED)
    assert row["curie_in_metadata"] is True
    assert row["absent_status"] == "curie_present"
    assert row["absent_by_rule"] is False
    assert row["qualifies"] is False


def test_a_child_of_a_covered_parent_is_not_absent():
    row = ood.verdict(_candidate(related_covered=("UBERON:0002048",), notes="lung contains it"), COVERED)
    assert row["curie_in_metadata"] is False
    assert row["related_covered"] == ["UBERON:0002048"]
    assert row["absent_status"] == "related_to_covered"
    assert row["absent_by_rule"] is False
    assert row["qualifies"] is False


def test_a_parent_of_a_covered_child_is_not_absent_either():
    row = ood.verdict(_candidate(curie="CL:0000542", related_covered=("CL:0000084",)), COVERED)
    assert row["absent_status"] == "related_to_covered"
    assert row["absent_by_rule"] is False


def test_an_unchecked_relation_is_never_counted_absent():
    row = ood.verdict(_candidate(relation_checked=False), COVERED)
    assert row["absent_status"] == "relation_unchecked"
    assert row["absent_by_rule"] is False
    assert row["qualifies"] is False


def test_absent_plus_both_modalities_plus_open_licence_qualifies():
    row = ood.verdict(_candidate(), COVERED)
    assert row["absent_status"] == "absent"
    assert row["absent_by_rule"] is True
    assert row["qualifies"] is True


def test_absent_but_missing_a_modality_or_closed_licence_does_not_qualify():
    assert ood.verdict(_candidate(accessibility=None), COVERED)["qualifies"] is False
    assert ood.verdict(_candidate(rna=None), COVERED)["qualifies"] is False
    assert ood.verdict(_candidate(licence=None), COVERED)["qualifies"] is False
    assert ood.verdict(_candidate(disqualified="dbGaP, apply"), COVERED)["qualifies"] is False


def test_empty_or_all_disqualified_reads_no_source_qualifies():
    assert ood.inventory([], COVERED)["answer"] == "NO_SOURCE_QUALIFIES"
    out = ood.inventory([_candidate(curie="CL:0000084")], COVERED)
    assert out["answer"] == "NO_SOURCE_QUALIFIES"
    assert out["qualifying_count"] == 0
    assert out["covered_curie_count"] == 3


def test_one_qualifying_candidate_flips_the_answer():
    out = ood.inventory([_candidate()], COVERED)
    assert out["answer"] == "QUALIFIES"
    assert out["qualifying_count"] == 1


def test_the_rule_states_what_absent_from_the_metadata_does_not_mean():
    limit = ood.SELECTION_RULE["wording_limit"]
    assert "not 'absent from the training set'" in limit
    assert "what we hold, not what AlphaGenome saw" in limit
    assert "adds no biosample" in ood.SELECTION_RULE["no_biosample"]
    assert "2026-09-27 stop list" in ood.SELECTION_RULE["no_biosample"]


def test_committed_rule_file_matches_the_module_and_names_its_set():
    path = Path("data/results/oodinventory_rule.json")
    if not path.exists():  # pragma: no cover - written by scripts/oodinventory_curies.py
        return
    payload = json.loads(path.read_text())
    assert payload["selection_rule"] == ood.SELECTION_RULE
    assert payload["covered_curie_count"] == len(payload["covered_curies"])
    assert payload["covered_curie_digest"] == ood.curie_digest(payload["covered_curies"])
    assert payload["rna_rows_with_accession"] == 0
