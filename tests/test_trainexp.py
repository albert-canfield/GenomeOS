# SPDX-License-Identifier: AGPL-3.0-or-later
"""What `genomeos.attribution.trainexp` must not be allowed to say.

This lane's one risk is a wording risk: an "unknown" read as "not a training track", or a
biosample+assay match reported as if it identified an experiment. The tests below are about that
first, and about the mechanics second.

The tests that need AlphaGenome's recorded track metadata carry `needs_local_data`: the file is
git-ignored machine-local data, so in a bare worktree or in CI they SKIP BY NAME rather than fail.
No assertion is weakened by the marker; it moves where a test runs, never what it claims.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genomeos.attribution import context_evidence as ce
from genomeos.attribution import freshrna as fr
from genomeos.attribution import trainexp as tx

needs_track_metadata = pytest.mark.needs_local_data(
    str(tx.METADATA), how="scripts/entex_feasibility.py writes it (AG_METADATA)"
)

RESULT = Path("data/results/trainexp_independence.json")


# ---- (a) what the wording may never do ---------------------------------------------------------


def test_unknown_is_never_described_as_independence():
    """`CLASS_UNKNOWN` must say in its own text that it is not evidence of independence."""
    text = tx.CLASSES[tx.CLASS_UNKNOWN].lower()
    assert "not a training track" in text
    assert "absence of a basis" in text


def test_probable_calls_its_own_method_suggestive_and_matched_calls_its_own_strong():
    assert "suggestive only" in tx.CLASSES[tx.CLASS_PROBABLE].lower()
    assert "suggestive only" in tx.METHOD_BIOSAMPLE_ASSAY.lower()
    assert "strong" in tx.CLASSES[tx.CLASS_MATCHED].lower()
    assert "strong" in tx.METHOD_ACCESSION.lower()


def test_probable_does_not_claim_the_experiment_is_the_one_behind_the_track():
    """A biosample+assay match may say the model saw SOME such RNA, never WHICH experiment."""
    text = tx.CLASSES[tx.CLASS_PROBABLE]
    assert "not that this experiment is the one behind that track" in text.lower()


def test_every_quantity_is_labelled_post_hoc_and_decides_nothing():
    for text in (tx.POST_HOC, tx.ELEVEN_CELL_NOTE):
        low = text.lower()
        assert "post hoc" in low
        assert "nothing" in low
    assert "decides nothing" in tx.POST_HOC.lower()


def test_the_training_set_record_states_it_is_not_documented_at_track_level():
    low = tx.TRAINING_SET_RECORD.lower()
    assert "not documented at track level" in low
    assert "does not say which experiment" in low


def test_the_eleven_cell_note_does_not_claim_to_test_training_independence():
    assert "not a test of independence from alphagenome's training data" in tx.ELEVEN_CELL_NOTE.lower()


# ---- (b) the sentence carried beside every M2i "better" -----------------------------------------


def _entry(klass: str) -> dict:
    return {"class": klass}


def test_sentence_says_unknown_when_nothing_is_settled_by_accession():
    entries = [_entry(tx.CLASS_PROBABLE)] * 13
    s = tx.independence_sentence(entries)
    assert s.startswith("fresh to this project; its independence from AlphaGenome's training data is")
    assert "unknown for 13 of 13" in s
    assert "same biosample and assay exists for 13 of 13" in s
    assert "settles nothing either way" in s


def test_sentence_counts_only_accession_matches_as_settled():
    entries = [_entry(tx.CLASS_MATCHED)] * 4 + [_entry(tx.CLASS_PROBABLE)] * 9
    assert "unknown for 9 of 13" in tx.independence_sentence(entries)


def test_sentence_never_says_independent():
    for entries in (
        [_entry(tx.CLASS_UNKNOWN)] * 13,
        [_entry(tx.CLASS_PROBABLE)] * 13,
        [_entry(tx.CLASS_MATCHED)] * 13,
    ):
        assert "is independent" not in tx.independence_sentence(entries).lower()


# ---- (c) the classification mechanics ----------------------------------------------------------


def test_accession_match_wins_and_names_its_method():
    out = tx.classify_one(
        "K562",
        "ENCSR000AEP",
        "K562",
        "EFO:0002067",
        "polyA plus RNA-seq",
        {"ENCSR000AEP": ["rna_seq"]},
        {},
        {},
    )
    assert out["class"] == tx.CLASS_MATCHED
    assert out["method"] == tx.METHOD_ACCESSION


def test_biosample_assay_match_is_probable_not_matched():
    out = tx.classify_one(
        "K562",
        "ENCSR000AEP",
        "K562",
        "EFO:0002067",
        "polyA plus RNA-seq",
        {},
        {("EFO:0002067", "polyA plus RNA-seq"): ["EFO:0002067 polyA plus RNA-seq"]},
        {},
    )
    assert out["class"] == tx.CLASS_PROBABLE
    assert out["method"] == tx.METHOD_BIOSAMPLE_ASSAY


def test_same_biosample_other_assay_is_unknown_and_says_which_assays_exist():
    out = tx.classify_one(
        "astrocyte",
        "ENCSR920UAO",
        "astrocyte",
        "CL:0000127",
        "polyA plus RNA-seq",
        {},
        {("CL:0000127", "total RNA-seq"): ["CL:0000127 total RNA-seq"]},
        {"CL:0000127": ["total RNA-seq"]},
    )
    assert out["class"] == tx.CLASS_UNKNOWN
    assert "total RNA-seq" in out["evidence"]


def test_a_biosample_absent_altogether_is_unknown_and_says_so():
    out = tx.classify_one("x", "ENCSR000XXX", "x", "CL:9999999", "total RNA-seq", {}, {}, {})
    assert out["class"] == tx.CLASS_UNKNOWN
    assert "none at all" in out["evidence"]


def test_tally_counts_experiments_and_files_per_class():
    entries = [
        {"class": tx.CLASS_PROBABLE, "files": 2},
        {"class": tx.CLASS_PROBABLE, "files": 3},
        {"class": tx.CLASS_UNKNOWN, "files": 1},
    ]
    t = tx.tally(entries)
    assert t[tx.CLASS_PROBABLE] == {"experiments": 2, "files": 5}
    assert t[tx.CLASS_UNKNOWN] == {"experiments": 1, "files": 1}
    assert t[tx.CLASS_MATCHED] == {"experiments": 0, "files": 0}


def test_cell_terms_covers_every_freshrna_cell_and_matches_the_authority():
    terms = tx.cell_terms()
    assert set(fr.CELLS) <= set(terms)
    for cell, (name, curie) in terms.items():
        assert ce.ident(name) == cell
        assert ce.READER_TERMS[name] == curie


def test_absent_metadata_reports_rather_than_fetching(tmp_path):
    with pytest.raises(tx.MetadataAbsentError) as e:
        tx.read_tracks(tmp_path / "nowhere.csv")
    assert "no independence is claimed" in str(e.value)


# ---- (d) the eleven-cell subset ----------------------------------------------------------------


def test_subset_drops_exactly_the_two_cells_m2i_disclosed():
    keep = tx.subset_indices(fr.CELLS, tx.SHARED_BIOSAMPLE_CELLS)
    dropped = tx.excluded_cells(fr.CELLS, tx.SHARED_BIOSAMPLE_CELLS)
    assert len(keep) == 11
    assert dropped == ["hepatocyte", "cardiac_muscle_cell"] or set(dropped) == {
        "hepatocyte",
        "cardiac_muscle_cell",
    }


def test_the_two_dropped_cells_are_the_ones_the_sources_file_names():
    sources = json.loads(fr.SOURCES.read_text())
    assert sorted(sources["non_overlap"]["cells_sharing_a_read_biosample"]) == sorted(
        tx.SHARED_BIOSAMPLE_CELLS
    )


def test_nothing_in_this_lane_may_fit_anything():
    with pytest.raises(fr.NoFitError):
        fr.fit_scale()


# ---- (e) the committed result says what this module says ---------------------------------------


@needs_track_metadata
@pytest.mark.skipif(not RESULT.exists(), reason="data/results/trainexp_independence.json not written yet")
def test_the_result_reports_zero_accession_matches_for_rna_seq():
    d = json.loads(RESULT.read_text())
    m = d["matching"]
    assert m["rna_seq_rows_carrying_an_accession"] == 0
    assert m["accession_matching_available_for_rna_seq"] is False
    assert m["fresh_thirteen"]["tally"][tx.CLASS_MATCHED]["experiments"] == 0
    assert m["already_read_twenty_seven"]["tally"][tx.CLASS_MATCHED]["experiments"] == 0


@needs_track_metadata
@pytest.mark.skipif(not RESULT.exists(), reason="data/results/trainexp_independence.json not written yet")
def test_the_result_covers_both_whole_sets():
    d = json.loads(RESULT.read_text())
    fresh = d["matching"]["fresh_thirteen"]["tally"]
    read = d["matching"]["already_read_twenty_seven"]["tally"]
    assert sum(v["experiments"] for v in fresh.values()) == 13
    assert sum(v["files"] for v in fresh.values()) == 22
    assert sum(v["experiments"] for v in read.values()) == 13
    assert sum(v["files"] for v in read.values()) == 27


@needs_track_metadata
@pytest.mark.skipif(not RESULT.exists(), reason="data/results/trainexp_independence.json not written yet")
def test_the_eleven_cell_subset_reproduces_m2i_on_thirteen_and_fits_nothing():
    d = json.loads(RESULT.read_text())["eleven_cell_subset"]
    assert d["thirteen_reproduces_the_committed_result"] is True
    assert d["scale"] == fr.SCALE
    assert len(d["cells_kept"]) == 11
    assert sorted(d["cells_dropped"]) == sorted(tx.excluded_cells(fr.CELLS, tx.SHARED_BIOSAMPLE_CELLS))
    assert "POST HOC" in d["reading"]
