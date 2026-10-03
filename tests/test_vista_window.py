# SPDX-License-Identifier: AGPL-3.0-or-later
"""vista.py asks whether any gene at the bar agrees with the element's observed tissue groups.

Wave 2 of section 5 item 10, under the registration at `genomeos/attribution/onetarget2.py`
(24adf33), and the first module of the wave to need an AMENDMENT to it: `onetarget2.AMENDMENTS`
records that the shared reader gains the TRACK beside the gene at the bar
(`attribution/targets.tracks_at_bar`), because vista's whole claim is about the predicted track
GROUP and `WindowReading` drops the track.

The amendment's own risk is that it becomes a SECOND SIZE RULE. The test below asserts that
`tracks_at_bar`'s (gene, signed) projection is EQUAL to `genes_at_bar`'s output -- on fixtures and
on a real chromosome's cache -- so the two cannot drift apart unnoticed, and
`window_tissue_agreement` raises if they ever disagree on a row it is reading.

vista.py's registered falsifier: "a VISTA row whose `tissue_agrees` is False on the head and True
for some other gene at the bar. One such row falsifies the stored agreement rate as an agreement
rate of the ELEMENT; zero of them leaves the stored rate standing unchanged."
"""

from __future__ import annotations

import pytest

from genomeos.attribution import onetarget2 as ot
from genomeos.attribution import vista
from genomeos.attribution.targets import genes_at_bar, tracks_at_bar

#: Two real AlphaGenome track names, one in each of two VISTA groups, so the fixture varies on the
#: axis `group_of_track` turns on rather than reusing one tissue everywhere.
NEURAL = "brain"
LIMB = "muscle of trunk"


def _g(gene: str, drop: float, rise: float, drop_tissue=NEURAL, rise_tissue=LIMB) -> dict:
    return {
        "gene": gene,
        "max_drop_log2fc": drop,
        "max_rise_log2fc": rise,
        "max_drop_tissue": drop_tissue,
        "max_rise_tissue": rise_tissue,
    }


class _Reader:
    def __init__(self, records):
        self.records = records

    def element(self, chrom, element_id):
        return self.records.get(element_id)


#: R1's head reaches the bar by its FALL on a neural track and a second gene reaches it by its RISE
#: on a limb track; R2's window names only its head; R3 is not cached; R4 was skipped as over-length.
#: R6 holds a gene whose FALL AND RISE ARE EQUAL, which is the axis the size rule's TIE ORDER turns
#: on. Without it the `fall >= up` tie clause is unreached: mutating it to `fall > up` left the
#: equality test green, which is how it was found.
RECORDS = {
    "R1": {"genes": [_g("APP", -0.8, 0.05), _g("SOD1", -0.1, 0.4)]},
    "R2": {"genes": [_g("APP", -0.5, 0.01)]},
    "R3": None,
    "R6": {"genes": [_g("TIED", -0.4, 0.4)]},
}
READER = _Reader(RECORDS)


def _row(rid, head, head_group, agrees, groups, skipped=None, tissues=None):
    #: The four fields `summarise` reads on EVERY row and the committed vista_chr*.json files do
    #: not carry: enhancer_like, any_ccre, domain and constrained_fraction. They are here so the
    #: summarise control can be asserted on fixtures, which is the only place it can be -- the
    #: stored rows were written from richer in-memory rows and summarise cannot be re-run on them.
    r = {
        "id": rid,
        "status": "positive",
        "predicted": {"gene": head, "tissue": NEURAL, "strength": "strong"} if head else None,
        "predicted_group": head_group,
        "tissue_agrees": agrees,
        "skipped": skipped,
        "enhancer_like": 1,
        "any_ccre": 1,
        "domain": "chr21:D1",
        "constrained_fraction": 0.3,
        "verdict_coding": "agrees with nearest TSS in domain",
    }
    if groups is not None:
        r["groups"] = groups
    if tissues is not None:
        r["tissues"] = tissues
    return r


# --------------------------------------------------------------------------------------------
# the amendment is not a second size rule
# --------------------------------------------------------------------------------------------


def test_tracks_at_bar_projects_exactly_onto_genes_at_bar_on_fixtures():
    for rec in RECORDS.values():
        bar = tracks_at_bar(rec)
        plain = genes_at_bar(rec)
        if rec is None:
            assert bar is None and plain is None
            continue
        assert [(g, v) for g, v, _ in bar] == plain


def test_a_tie_between_the_fall_and_the_rise_goes_to_the_fall_in_both_readers():
    """`predict_target`'s rule is "ties go to activation", so a tie takes the FALL and the fall's
    track. Both readers must take the same side or they are two rules."""
    bar = tracks_at_bar(RECORDS["R6"])
    assert bar == [("TIED", -0.4, NEURAL)]
    assert genes_at_bar(RECORDS["R6"]) == [("TIED", -0.4)]


def test_tracks_at_bar_keeps_the_track_of_whichever_side_won():
    bar = tracks_at_bar(RECORDS["R1"])
    assert bar[0] == ("APP", -0.8, NEURAL)  # the fall won, so the fall's track
    assert bar[1] == ("SOD1", 0.4, LIMB)  # the rise won, so the rise's track


def test_tracks_at_bar_projects_exactly_onto_genes_at_bar_on_real_chromosome_data(needs_chr21=None):
    """The same equality on the real cache, where the fixtures cannot reach: the measurement, not
    the assumption. Skipped by name when the git-ignored response cache is not on this machine."""
    from pathlib import Path

    from genomeos.attribution.targets import ELEMENT_CACHE, ElementResponses

    if not (Path(ELEMENT_CACHE) / "chr21.json.gz").exists():
        pytest.skip("the sweep's response cache is git-ignored and is not in this checkout")
    responses = ElementResponses()
    responses._load("chr21")
    checked = 0
    for eid in list(responses._archive)[:2000]:
        rec = responses.element("chr21", eid)
        assert [(g, v) for g, v, _ in tracks_at_bar(rec)] == genes_at_bar(rec), eid
        checked += 1
    assert checked > 0


# --------------------------------------------------------------------------------------------
# the control and the populations kept apart
# --------------------------------------------------------------------------------------------


def test_without_a_reader_there_is_no_arm_and_no_new_key():
    assert vista.window_tissue_agreement([], "chr21") == {}
    rows = [
        _row("R1", "APP", "neural", False, ["limb and mesenchyme"]),
    ]
    both = [dict(r) for r in rows] + [dict(r, status="negative") for r in rows]
    off = vista.summarise([dict(r) for r in both])
    on = vista.summarise([dict(r) for r in both], "chr21", READER)
    assert "window" not in off
    for k in off:
        assert off[k] == on[k], k
    assert set(on) - set(off) == {"window"}


def test_a_reader_without_the_chromosome_refuses():
    """PLANT: the rows do not carry the chromosome and the cache is keyed by it, so a window could
    be read from the wrong one. The match is on 'from the wrong one', only this refusal's phrase."""
    rows = [_row("R1", "APP", "neural", False, ["limb and mesenchyme"])]
    with pytest.raises(ValueError, match="from the wrong one"):
        vista.summarise(rows, None, READER)


def test_a_skipped_row_and_an_unjudged_row_are_counted_apart():
    rows = [
        _row("R4", None, None, None, ["neural"], skipped="longer than 100000 bp"),
        _row("R2", "APP", "neural", None, []),  # no observed groups at all
        _row("R2", "APP", "neural", True, ["neural"]),
    ]
    w = vista.window_tissue_agreement(rows, "chr21", READER)["window"]
    assert w["rows_skipped_as_over_length"] == 1
    assert w["rows_with_no_observed_groups"] == 1
    assert w["rows_judged_under_the_window"] == 1
    assert w["rows_not_cached"] == 0
    assert "never as an element whose window named nothing" in w["skipped_is_not_unpredicted"]


def test_a_row_the_cache_does_not_hold_is_counted_by_name():
    rows = [_row("R3", "APP", "neural", False, ["neural"])]
    w = vista.window_tissue_agreement(rows, "chr21", READER)["window"]
    assert w["rows_not_cached"] == 1 and w["rows_judged_under_the_window"] == 0
    assert w["rows_agreeing_on_a_gene_the_head_dropped"] == 0


def test_the_observed_groups_of_a_committed_row_are_rederived_from_its_tissues():
    """A committed vista_chr*.json row carries `tissues` and not `groups`. Deriving them with the
    module's own rule is what makes the falsifier askable from the stored files at all."""
    tissues = sorted(vista.GROUP_OF)[:1]
    expected = sorted(vista.tissue_groups(tissues))
    rows = [_row("R2", "APP", expected[0], True, None, tissues=tissues)]
    w = vista.window_tissue_agreement(rows, "chr21", READER)["window"]
    assert w["rows_with_no_observed_groups"] == 0
    assert w["rows_judged_under_the_window"] == 1


# --------------------------------------------------------------------------------------------
# the falsifier, and the refusals
# --------------------------------------------------------------------------------------------


def test_a_row_that_agrees_on_a_gene_the_head_dropped_fires():
    """R1's head APP sits on a neural track; SOD1 at the bar sits on a limb track, which is what
    the element actually expresses. The head said no and the element says yes."""
    rows = [_row("R1", "APP", "neural", False, ["limb and mesenchyme"])]
    w = vista.window_tissue_agreement(rows, "chr21", READER)["window"]
    assert w["rows_agreeing_on_a_gene_the_head_dropped"] == 1
    f = w["flipped"][0]
    assert f["head"] == "APP" and f["gene_at_the_bar"] == "SOD1"
    assert f["its_group"] == "limb and mesenchyme" and f["its_log2_fold_change"] == 0.4
    assert "FIRES" in w["verdict_against_its_own_falsifier"]


def test_a_row_whose_head_already_agrees_is_not_looked_at():
    rows = [_row("R1", "APP", "neural", True, ["neural"])]
    w = vista.window_tissue_agreement(rows, "chr21", READER)["window"]
    assert w["rows_agreeing_on_a_gene_the_head_dropped"] == 0
    assert "does NOT fire" in w["verdict_against_its_own_falsifier"]


def test_a_second_gene_whose_track_has_no_vista_group_contributes_nothing():
    reader = _Reader({"R5": {"genes": [_g("APP", -0.8, 0.05), _g("SOD1", -0.1, 0.4, rise_tissue="K562")]}})
    assert vista.group_of_track("K562") is None
    rows = [_row("R5", "APP", "neural", False, ["limb and mesenchyme"])]
    w = vista.window_tissue_agreement(rows, "chr21", reader)["window"]
    assert w["rows_agreeing_on_a_gene_the_head_dropped"] == 0


def test_the_two_readers_disagreeing_refuses_rather_than_reporting(monkeypatch):
    """PLANT: the amendment has become a second size rule. The match is on 'a second size rule',
    a phrase only this refusal uses."""
    monkeypatch.setattr(vista, "window_tissue_agreement", vista.window_tissue_agreement)  # keep the real one
    import genomeos.attribution.targets as tg

    monkeypatch.setattr(tg, "tracks_at_bar", lambda record, min_effect=0.1: [("WRONG", -9.9, NEURAL)])
    rows = [_row("R1", "APP", "neural", False, ["limb and mesenchyme"])]
    with pytest.raises(ValueError, match="a second size rule"):
        vista.window_tissue_agreement(rows, "chr21", READER)


def test_an_any_gene_head_disagreement_stops_the_wave():
    """PLANT: the stored head is not the window's strongest gene."""
    rows = [_row("R1", "SOD1", "limb and mesenchyme", False, ["neural"])]
    with pytest.raises(ot.InvariantRefutedError, match="THE WAVE STOPS"):
        vista.window_tissue_agreement(rows, "chr21", READER)
