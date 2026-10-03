# SPDX-License-Identifier: AGPL-3.0-or-later
"""closure.py moved onto the window, against its own registered falsifier.

Registered in `genomeos/attribution/onetarget2.py` (24adf33, aed8ae9). closure.py's falsifier:
"a gene whose closure gains or loses an element when the elements are grouped by every gene at the
bar instead of by the head alone."

The control comes first and is the thing that must hold: `attributed_elements` is NOT changed, and
`window_elements(chrom, None)` reproduces it key for key and value for value, so no committed
closure figure can move. Only then is the window grouping asked anything.
"""

from __future__ import annotations

import gzip
import json

import pytest

from genomeos.attribution import closure as cl
from genomeos.attribution import onetarget2 as ot
from genomeos.attribution.targets import ElementResponses


def _g(name, drop, rise=None, cells=None):
    g = {"gene": name, "max_drop_log2fc": drop, "by_cell": dict(cells or {})}
    if rise is not None:
        g["max_rise_log2fc"] = rise
    return g


def _element(eid, start, end, head, log2, action="activates", by_cell=None, any_gene=None):
    """One compact row. `any_gene` is the `predicted` head when it is NOT the coding head.

    A fixture that set `predicted` and `predicted_coding` to the same gene while the window held a
    larger non-coding fall would describe a table the sweep cannot write -- `predicted` is the
    maximum over EVERY gene in the window and `predicted_coding` over the coding ones. The first
    draft of this file did exactly that and `check_head_invariant` refused it, which is the refusal
    working on its author before it worked on anything else.
    """
    pc = {"gene": head, "log2_fold_change": log2, "action": action, "tissue": "K562"}
    p = dict(pc) if any_gene is None else {**pc, "gene": any_gene[0], "log2_fold_change": any_gene[1]}
    return {
        "id": eid,
        "start": start,
        "end": end,
        "predicted": p,
        "predicted_coding": pc,
        "predicted_coding_by_cell": dict(by_cell or {}),
    }


@pytest.fixture
def world(tmp_path, monkeypatch):
    """Three elements: two cached with a second coding gene at the bar, one not cached at all."""
    els = [
        _element("E1", 100, 200, "TOP", -0.55, by_cell={"K562": -0.40}, any_gene=("NONCODING", -0.90)),
        _element("E2", 300, 400, "SOLO", -0.30),
        _element("E3", 500, 600, "ONLY", 0.42, action="represses"),
    ]
    monkeypatch.setattr(cl, "attributed_elements", cl.attributed_elements)
    monkeypatch.setattr("genomeos.attribution.targets.attributed", lambda c, *a, **k: els)
    archive = {
        "E1": {
            "id": "E1",
            "genes": [
                _g("TOP", -0.55, cells={"K562": -0.40}),
                _g("SECOND", -0.33, cells={"K562": -0.21}),
                _g("NONCODING", -0.90),
                _g("BELOW", -0.02),
            ],
        },
        "E2": {"id": "E2", "genes": [_g("SOLO", -0.30), _g("TINY", -0.01)]},
    }
    with gzip.open(tmp_path / "chrT.json.gz", "wt") as fh:
        json.dump(archive, fh)
    return els, ElementResponses(tmp_path)


CODING = {"TOP", "SECOND", "SOLO", "ONLY", "BELOW", "TINY"}


# ---------------------------------------------------------------------------
# the control: no committed closure figure may move
# ---------------------------------------------------------------------------


def test_the_window_grouping_with_no_reader_IS_the_committed_grouping(world):
    """`responses=None` opens nothing and reproduces attributed_elements exactly."""
    _els, _r = world
    old = cl.attributed_elements("chrT")
    new, census = cl.window_elements("chrT", None, CODING)
    assert new == old
    assert census["elements"] == 3 and census["not_cached"] == 3
    assert census["memberships"] == 3


def test_the_two_readings_compute_the_head_effect_to_the_bit(world):
    """effect_of(signed) and sign * abs(log2_fold_change) are the same arithmetic, not nearly."""
    _els, r = world
    old = cl.attributed_elements("chrT")
    new, _ = cl.window_elements("chrT", r, CODING)
    for gene, rows in old.items():
        by_id = {x["id"]: x for x in new.get(gene, [])}
        for row in rows:
            assert row["id"] in by_id, (gene, row["id"])
            assert by_id[row["id"]]["effect"] == pytest.approx(row["effect"], abs=1e-12)
    assert cl.effect_of(-0.55) == 0.55  # a fall on deletion: the element activates
    assert cl.effect_of(0.42) == -0.42  # a rise on deletion: it represses


def test_an_element_the_cache_does_not_hold_keeps_its_compact_reading_and_is_counted(world):
    """Never a zero and never dropped: E3 is in no archive and still reaches ONLY's closure."""
    _els, r = world
    new, census = cl.window_elements("chrT", r, CODING)
    assert census["not_cached"] == 1
    assert [x["id"] for x in new["ONLY"]] == ["E3"]
    assert new["ONLY"][0]["effect"] == pytest.approx(-0.42)


def test_without_a_coding_set_the_coding_window_is_not_formed_and_says_so(world):
    """Forming it against an uncoded window would compare two populations."""
    _els, r = world
    new, census = cl.window_elements("chrT", r, None)
    assert census["no_coding_window"] == 2 and census["not_cached"] == 1
    assert new == cl.attributed_elements("chrT")


# ---------------------------------------------------------------------------
# the falsifier
# ---------------------------------------------------------------------------


def test_the_registered_falsifier_FIRES_a_gene_gains_an_element(world):
    """SECOND is at the bar at E1 and is not its head, so E1 joins a closure it was absent from."""
    _els, r = world
    old = cl.attributed_elements("chrT")
    new, census = cl.window_elements("chrT", r, CODING)
    assert "SECOND" not in old
    assert [x["id"] for x in new["SECOND"]] == ["E1"]
    assert census["memberships"] > census["elements"]
    gained = {g: sorted({x["id"] for x in new[g]} - {x["id"] for x in old.get(g, [])}) for g in new}
    assert {g: v for g, v in gained.items() if v} == {"SECOND": ["E1"]}


def test_no_gene_LOSES_an_element_because_the_head_is_always_at_its_own_bar(world):
    """The other half of the falsifier, and the half that must not fire."""
    _els, r = world
    old = cl.attributed_elements("chrT")
    new, _ = cl.window_elements("chrT", r, CODING)
    for gene, rows in old.items():
        assert {x["id"] for x in rows} <= {x["id"] for x in new.get(gene, [])}, gene


def test_a_noncoding_gene_at_the_bar_does_NOT_enter_a_coding_closure(world):
    """NONCODING has the largest fall at E1 and is outside the coding set, so it is excluded."""
    _els, r = world
    new, _ = cl.window_elements("chrT", r, CODING)
    assert "NONCODING" not in new


def test_a_gene_below_the_bar_does_not_enter_either(world):
    _els, r = world
    new, _ = cl.window_elements("chrT", r, CODING)
    assert "BELOW" not in new and "TINY" not in new


def test_a_non_head_genes_per_cell_effects_come_from_its_own_tracks_and_absence_is_not_zero(world):
    _els, r = world
    new, _ = cl.window_elements("chrT", r, CODING)
    assert new["SECOND"][0]["by_cell"] == {"K562": pytest.approx(0.21)}
    assert new["SOLO"][0]["by_cell"] == {}  # the sweep scored no cell track for it: absent, not 0.0


# ---------------------------------------------------------------------------
# the one invariant the wave stands on
# ---------------------------------------------------------------------------


def test_this_move_does_not_refute_the_invariant(world):
    _els, r = world
    _new, census = cl.window_elements("chrT", r, CODING)
    assert census["head_disagrees"] == 0
    ot.check_head_invariant(census, "closure.window_elements")


def test_PLANTED_a_head_the_window_disagrees_with_STOPS_THE_WAVE(tmp_path, monkeypatch):
    """A tree where the table's head is not the window's: the wave must stop, not report.

    Mutation: delete the raise in onetarget2.check_head_invariant and this test fails.
    """
    els = [_element("E1", 100, 200, "CLAIMED", -0.10)]
    monkeypatch.setattr("genomeos.attribution.targets.attributed", lambda c, *a, **k: els)
    with gzip.open(tmp_path / "chrT.json.gz", "wt") as fh:
        json.dump({"E1": {"genes": [_g("OTHER", -0.80), _g("CLAIMED", -0.10)]}}, fh)
    _new, census = cl.window_elements("chrT", ElementResponses(tmp_path), {"CLAIMED", "OTHER"})
    assert census["head_disagrees"] == 1
    with pytest.raises(ot.InvariantRefutedError):
        ot.check_head_invariant(census, "closure.window_elements")


def test_a_coding_head_disagreement_is_counted_by_name_and_does_not_stop_the_wave(tmp_path, monkeypatch):
    """The registered limit: 10 of 4,794 elements had no predicted_coding where the cache named one."""
    e = _element("E1", 100, 200, "HEAD", -0.50, any_gene=("CODINGBIGGER", -0.70))
    monkeypatch.setattr("genomeos.attribution.targets.attributed", lambda c, *a, **k: [e])
    with gzip.open(tmp_path / "chrT.json.gz", "wt") as fh:
        json.dump({"E1": {"genes": [_g("HEAD", -0.50), _g("CODINGBIGGER", -0.70)]}}, fh)
    r = ElementResponses(tmp_path)
    _new, census = cl.window_elements("chrT", r, {"HEAD", "CODINGBIGGER"})
    # the any-gene head agrees, so the wave does not stop; the CODING head does not, and is counted
    assert census["head_disagrees"] == 0
    assert census["coding_head_disagrees"] == 1
    ot.check_head_invariant(census, "closure.window_elements")
