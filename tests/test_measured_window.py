# SPDX-License-Identifier: AGPL-3.0-or-later
"""measured.py moved onto the window, against its own registered falsifier.

Registered in `genomeos/attribution/onetarget2.py` (24adf33, aed8ae9). measured.py's falsifier:
"an element whose measured regulated gene is at the bar in the window but is not the head."

That question is asked of the SWEEP's window and is not new evidence about the element. It cannot
turn a disagreement into an agreement of the compiled claim -- the compiled claim names one gene --
and the only thing it can establish is whether the claim was narrower than the model it came from.

Every test here is synthetic. The real-data arm is `scripts/onetarget2_run.py`.
"""

from __future__ import annotations

import pytest

from genomeos.attribution import measured as me

#: The control below is a reading of the REAL rows: a planted layer would make it assert against a
#: fiction. data/knowledge is git-ignored machine-local data, so in a fresh checkout it skips BY NAME.
needs_the_response_cache = pytest.mark.needs_local_data(
    "data/knowledge/alphagenome/elements",
    how="scripts/enhancer_targets_all.py writes it (worker_scorer, threshold=0.0); "
    "data/knowledge is git-ignored machine-local data, so a fresh checkout cannot have it",
)


def _crispri(regulated=(), not_regulated=(), well_powered=None):
    c = {"genes_regulated": list(regulated), "genes_not_regulated": list(not_regulated)}
    if well_powered is not None:
        c["genes_no_effect_well_powered"] = list(well_powered)
    return {"crispri": c}


# ---------------------------------------------------------------------------
# the control: no committed field may move
# ---------------------------------------------------------------------------


@needs_the_response_cache
def test_with_no_reader_a_row_does_not_gain_the_key_at_all_on_the_real_chromosome():
    """`responses=None` is the control every consumer keeps: not one field added, not one changed.

    Asserted on the real chr21 rows, because the claim is about this module's committed output and a
    synthetic layer would only prove it about a fiction. The genome-wide form of the same check is
    `scripts/onetarget2_run.py`'s own control.
    """
    from pathlib import Path

    from genomeos.attribution.targets import ElementResponses
    from genomeos.genome import Annotation

    ch = "chr21"
    coding = {
        g.symbol
        for g in Annotation.from_gff3(
            Path("data/reference") / f"gencode_v50_{ch}.gff3.gz", {ch}
        ).genes.values()
        if g.type == "protein_coding"
    }
    plain = me.rows(ch)
    assert plain, "chr21 has measured rows on this machine"
    assert all("window" not in r for r in plain)
    windowed = me.rows(ch, responses=ElementResponses(), coding=coding)
    assert len(windowed) == len(plain)
    for a, b in zip(plain, windowed, strict=True):
        assert {k: v for k, v in b.items() if k != "window"} == a, a["id"]
        assert "window" in b


def test_a_not_cached_window_is_named_rather_than_empty():
    m = _crispri(regulated=["TOP"])
    assert me.window_agreement(m, None)["not_cached"] is True


def test_the_compiled_claim_is_untouched_by_the_window_question():
    """agreement() reads ONE gene and window_agreement() reads the window; neither changes the other."""
    m = _crispri(regulated=["SECOND"], not_regulated=["TOP"])
    a = me.agreement("TOP", m)
    assert a["crispri"] == me.DISAGREES
    w = me.window_agreement(m, [("TOP", -0.5), ("SECOND", -0.3)])
    assert w["regulated_gene_in_window"] is True
    assert me.agreement("TOP", m) == a  # asking the window changed nothing about the claim


# ---------------------------------------------------------------------------
# the falsifier
# ---------------------------------------------------------------------------


def test_the_registered_falsifier_FIRES_when_the_measured_gene_is_at_the_bar_but_not_the_head():
    m = _crispri(regulated=["SECOND"], not_regulated=["TOP"])
    w = me.window_agreement(m, [("TOP", -0.50), ("SECOND", -0.30)])
    assert me.agreement("TOP", m)["crispri"] == me.DISAGREES
    assert w["regulated_gene_in_window"] is True
    assert w["regulated_genes_in_window"] == ["SECOND"]
    assert w["best_rank_of_a_regulated_gene"] == 2  # rank 1 is the head; 2 is what the table dropped


def test_rank_one_means_the_compact_table_already_had_it_and_nothing_was_dropped():
    m = _crispri(regulated=["TOP"])
    w = me.window_agreement(m, [("TOP", -0.50), ("SECOND", -0.30)])
    assert w["best_rank_of_a_regulated_gene"] == 1


def test_a_window_with_no_regulated_gene_says_False_and_not_None():
    m = _crispri(regulated=["ELSEWHERE"])
    w = me.window_agreement(m, [("TOP", -0.5)])
    assert w["regulated_gene_in_window"] is False and w["regulated_genes_in_window"] == []
    assert w["best_rank_of_a_regulated_gene"] is None


def test_an_element_the_cache_does_not_hold_is_None_and_NEVER_False():
    """A missing window is not an empty one. False here would be a claim nobody measured."""
    w = me.window_agreement(_crispri(regulated=["X"]), None)
    assert w["regulated_gene_in_window"] is None
    assert w["not_cached"] is True and w["genes_at_bar"] is None


def test_a_row_with_no_crispri_measurement_is_None_rather_than_False():
    """lentiMPRA and VISTA ask whether the sequence acts, not which gene; the question does not apply."""
    w = me.window_agreement({"lentimpra": {"cells_active": ["K562"]}}, [("TOP", -0.5)])
    assert w["regulated_gene_in_window"] is None
    assert w["genes_at_bar"] == 1 and w["not_cached"] is False


def test_several_regulated_genes_in_one_window_are_all_named_and_the_best_rank_is_the_nearest_head():
    m = _crispri(regulated=["THIRD", "SECOND"])
    w = me.window_agreement(m, [("TOP", -0.9), ("SECOND", -0.5), ("THIRD", -0.2)])
    assert w["regulated_genes_in_window"] == ["SECOND", "THIRD"]  # window order, strongest first
    assert w["best_rank_of_a_regulated_gene"] == 2
