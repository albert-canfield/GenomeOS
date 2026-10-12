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

import pytest

from genomeos.attribution import oodinventory as ood

#: `needs_local_data`, in the shape tests/test_trainexp.py and tests/test_astrorun.py already use for
#: this same file. `ood.rna_rows()` and `ood.covered_curies()` read AlphaGenome's recorded track
#: metadata, which lives in `data/cache` -- git-ignored machine-local data by the data boundary, so a
#: fresh clone and CI cannot have it and the tests below FAILED there with FileNotFoundError. They now
#: skip BY NAME, which is the difference between "not run here" and "passed". No assertion is weakened:
#: the covered set is still read from the file and never from a constant.
needs_track_metadata = pytest.mark.needs_local_data(
    str(ood.METADATA_PATH), how="scripts/entex_feasibility.py writes it (AG_METADATA)"
)

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


@needs_track_metadata
def test_covered_set_is_the_rna_rows_of_the_metadata_we_hold():
    rows = ood.rna_rows()
    assert rows, "the metadata copy has rna_seq rows"
    assert all(r["output"] == "rna_seq" for r in rows)
    curies = ood.covered_curies()
    # Every RNA row carries a CURIE, which is why the rule can match on it at all.
    assert all(r["ontology_curie"].strip() for r in rows)
    assert curies == {r["ontology_curie"].strip() for r in rows}


@needs_track_metadata
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


# --- the search, whose rule is the one above ---------------------------------------------------
# These hold the parts of scripts/oodinventory_search.py that decide what the list contains: the
# assay sets (a perturbation RNA-seq or a CAGE experiment must never count), the hand-declared
# relations OLS cannot serve, and the sibling disclosure.


def test_the_assay_sets_do_not_overlap_and_exclude_the_perturbations():
    from scripts import oodinventory_search as search

    assert not set(search.ACCESSIBILITY_ASSAYS) & set(search.RNA_ASSAYS)
    for assay in search.ACCESSIBILITY_ASSAYS + search.RNA_ASSAYS:
        assert assay not in search.EXCLUDED_ASSAYS
    for assay in ("shRNA RNA-seq", "CRISPR RNA-seq", "siRNA RNA-seq", "CRISPRi RNA-seq"):
        assert "genetically modified" in search.EXCLUDED_ASSAYS[assay]
    for assay in ("CAGE", "RAMPAGE", "PRO-cap", "small RNA-seq", "microRNA-seq"):
        assert assay in search.EXCLUDED_ASSAYS


@needs_track_metadata
def test_every_hand_declared_relation_names_a_covered_curie():
    from scripts import oodinventory_search as search

    covered = ood.covered_curies()
    for curie, (parents, evidence) in search.DECLARED_BY_HAND.items():
        assert parents, curie
        assert all(p in covered for p in parents), (curie, parents)
        assert len(evidence) > 40, curie
        # A hand declaration is only for a term OLS cannot answer: an NTR term, or an id
        # migration where the covered set carries the same name under the old id.
        assert curie.startswith("NTR:") or curie == "CL:4006001"


@needs_track_metadata
def test_a_hand_declared_term_is_not_absent_and_an_undeclared_ntr_term_is_unchecked():
    from scripts import oodinventory_search as search

    covered = ood.covered_curies()
    row = ood.verdict(
        ood.Candidate(
            curie="CL:4006001",
            label="fibroblast of skin of scalp",
            source="ENCODE",
            related_covered=search.DECLARED_BY_HAND["CL:4006001"][0],
            relation_checked=True,
        ),
        covered,
    )
    assert row["related_covered"] == ["NTR:0000524"]
    assert row["absent_by_rule"] is False
    unknown = ood.verdict(ood.Candidate(curie="NTR:9999999", label="invented", source="ENCODE"), covered)
    assert unknown["absent_status"] == "relation_unchecked"
    assert unknown["absent_by_rule"] is False


@needs_track_metadata
def test_the_sibling_flag_fires_on_a_donor_line_and_a_cell_line_and_not_on_a_tissue():
    from scripts import oodinventory_search as search

    covered = ood.covered_curies()
    assert search._sibling_flag("EFO:0010835", "GM19023", "cell line", covered)
    assert search._sibling_flag("EFO:0001247", "T47D", "cell line", covered)
    assert search._sibling_flag("UBERON:0002394", "bile duct", "tissue", covered) is None
    # A covered term is not a candidate at all, so it carries no flag.
    assert search._sibling_flag("EFO:0002784", "GM12878", "cell line", covered) is None


def test_a_named_donor_line_is_a_leaf_and_a_tissue_is_not():
    from scripts import oodinventory_search as search

    assert search._donor_line("GM19023", "cell line") is True
    assert search._donor_line("HG02571", "cell line") is True
    assert search._donor_line("T47D", "cell line") is False
    assert search._donor_line("bile duct", "tissue") is False


def test_the_ols_url_knows_efo_from_the_obo_ontologies_and_refuses_ntr():
    from scripts import oodinventory_search as search

    assert "www.ebi.ac.uk%252Fefo%252FEFO_0001247" in search.ols_url("EFO:0001247", "hierarchicalAncestors")
    assert "/ontologies/cl/" in search.ols_url("CL:0000451", "hierarchicalAncestors")
    assert "/ontologies/uberon/" in search.ols_url("UBERON:0002394", "hierarchicalAncestors")
    assert "/ontologies/clo/" in search.ols_url("CLO:0024877", "hierarchicalAncestors")
    assert search.ols_url("NTR:0000497", "hierarchicalAncestors") is None


def test_the_other_sources_each_name_the_clause_they_fail():
    from scripts import oodinventory_search as search

    named = {s["source"] for s in search.OTHER_SOURCES}
    assert {"GTEx", "Tabula Sapiens", "GEO / SRA"} <= named
    assert any("Human Cell Atlas" in s for s in named)
    for entry in search.OTHER_SOURCES:
        assert entry["clause_failed"] and len(entry["detail"]) > 60


@needs_track_metadata
def test_the_two_sibling_kinds_are_told_apart_and_a_tissue_has_neither():
    from scripts import oodinventory_search as search

    covered = ood.covered_curies()
    assert (
        search._sibling_kind("EFO:0010835", "GM19023", "cell line", covered)
        == "donor_variant_of_a_covered_lymphoblastoid_line"
    )
    assert (
        search._sibling_kind("EFO:0001247", "T47D", "cell line", covered)
        == "immortalised_or_cancer_line_with_covered_neighbours"
    )
    assert search._sibling_kind("UBERON:0002394", "bile duct", "tissue", covered) == ""
    # Each named line's flag names the covered neighbour it rests on, not a generic sentence.
    for curie, text in search.LINE_NEIGHBOURS.items():
        assert search._sibling_flag(curie, "x", "cell line", covered) == text
        assert ":" in text and "covered set holds" in text


def test_a_single_nucleus_only_row_gets_the_pseudobulk_caveat_and_a_bulk_row_does_not():
    from scripts import oodinventory_search as search

    rows = [
        {
            "curie": "UBERON:0002394",
            "label": "bile duct",
            "absent_by_rule": True,
            "rna": "ENCSR060SEA (snRNA-seq), ENCSR769KWM (snRNA-seq)",
        },
        {
            "curie": "EFO:0001247",
            "label": "T47D",
            "absent_by_rule": True,
            "rna": "ENCSR000BYD (polyA plus RNA-seq)",
        },
        {
            "curie": "EFO:0002067",
            "label": "K562",
            "absent_by_rule": False,
            "rna": "ENCSR000AAA (snRNA-seq)",
        },
    ]
    got = search._shortlist_caveats(rows)
    assert [c["curie"] for c in got] == ["UBERON:0002394"]
    assert "pseudobulk" in got[0]["caveat"]
