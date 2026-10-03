# SPDX-License-Identifier: AGPL-3.0-or-later
"""rearrangements.py asks whether the recipient gene is anywhere in the window, not only at a head.

Wave 2 of section 5 item 10, under the registration at `genomeos/attribution/onetarget2.py`
(24adf33). Its registered falsifier: "a rearrangement whose `names_the_recipient` goes False ->
True. This is a yes/no about a NAMED gene rather than about any gene, so THE_INVARIANT does not
protect it, and one flip establishes that the benchmark's negative was a property of the projection.
`genes_named`'s truncation to 8 must be widened or declared before the count is read, because a
longer list silently truncated would hide the gain."

The truncation is DECLARED and not widened: widening it would change a committed benchmark field,
while `genes_named_total` beside it says what the eight are a sample of. That is asserted below on a
fixture with MORE than eight named genes, which is the axis the demand turns on.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from genomeos.attribution import onetarget2 as ot
from genomeos.benchmark import rearrangements as rr


def _g(gene: str, drop: float, rise: float) -> dict:
    return {
        "gene": gene,
        "max_drop_log2fc": drop,
        "max_rise_log2fc": rise,
        "max_drop_tissue": "K562",
        "max_rise_tissue": "occipital lobe",
    }


class _Reader:
    def __init__(self, records):
        self.records = records

    def element(self, chrom, element_id):
        return self.records.get(element_id)

    def at_bar(self, chrom, element_id, min_effect=0.1, genes=None):
        from genomeos.attribution.targets import genes_at_bar

        got = genes_at_bar(self.element(chrom, element_id), min_effect)
        if got is None or genes is None:
            return got
        return [(g, v) for g, v in got if g in genes]


CH = SimpleNamespace(chrom="chr21")
CODING = {f"G{i}" for i in range(20)} | {"PAX3", "EPHA4"}


def _patch(monkeypatch, rows, records):
    monkeypatch.setattr(rr, "named_by_a_derived_layer", rr.named_by_a_derived_layer)
    import genomeos.benchmark.loci as loci

    monkeypatch.setattr(loci, "_deletion_rows", lambda chrom, lo, hi, rd: rows)


def _el(eid, head):
    return {"id": eid, "predicted_coding": {"gene": head}, "predicted": {"gene": head}}


def test_without_a_reader_the_answer_is_the_committed_one(monkeypatch, tmp_path):
    rows = [_el("E1", "EPHA4")]
    _patch(monkeypatch, rows, {"E1": {"genes": [_g("EPHA4", -0.8, 0.05), _g("PAX3", -0.1, 0.4)]}})
    out = rr.named_by_a_derived_layer(CH, tmp_path, 1000, "PAX3")
    assert out == {
        "elements_scored": 1,
        "genes_named": ["EPHA4"],
        "names_the_recipient": False,
        "pending": None,
    }


def test_with_a_reader_the_committed_keys_are_untouched_and_four_are_added(monkeypatch, tmp_path):
    rows = [_el("E1", "EPHA4")]
    reader = _Reader({"E1": {"genes": [_g("EPHA4", -0.8, 0.05), _g("PAX3", -0.1, 0.4)]}})
    _patch(monkeypatch, rows, None)
    off = rr.named_by_a_derived_layer(CH, tmp_path, 1000, "PAX3")
    on = rr.named_by_a_derived_layer(CH, tmp_path, 1000, "PAX3", reader, CODING)
    for k in off:
        assert off[k] == on[k], k
    assert set(on) - set(off) == {
        "genes_named_total",
        "genes_named_truncated_to",
        "window",
        "names_the_recipient_in_the_window",
        "the_claim_flipped",
    }


def test_the_claim_flips_false_to_true_and_says_so(monkeypatch, tmp_path):
    rows = [_el("E1", "EPHA4")]
    reader = _Reader({"E1": {"genes": [_g("EPHA4", -0.8, 0.05), _g("PAX3", -0.1, 0.4)]}})
    _patch(monkeypatch, rows, None)
    out = rr.named_by_a_derived_layer(CH, tmp_path, 1000, "PAX3", reader, CODING)
    assert out["names_the_recipient"] is False
    assert out["names_the_recipient_in_the_window"] is True
    assert out["the_claim_flipped"] is True
    assert out["window"]["genes_at_the_bar_and_not_a_head"] == ["PAX3"]
    assert "a genuine result and not a defect" in out["window"]["limit"]


def test_a_recipient_already_at_a_head_does_not_flip(monkeypatch, tmp_path):
    rows = [_el("E1", "PAX3")]
    reader = _Reader({"E1": {"genes": [_g("PAX3", -0.8, 0.05)]}})
    _patch(monkeypatch, rows, None)
    out = rr.named_by_a_derived_layer(CH, tmp_path, 1000, "PAX3", reader, CODING)
    assert out["names_the_recipient"] is True and out["the_claim_flipped"] is False


def test_the_truncation_is_declared_rather_than_widened(monkeypatch, tmp_path):
    """TWELVE named genes, so the committed `genes_named` IS truncated and the total says so."""
    rows = [_el(f"E{i}", f"G{i}") for i in range(12)]
    reader = _Reader({f"E{i}": {"genes": [_g(f"G{i}", -0.8, 0.05)]} for i in range(12)})
    _patch(monkeypatch, rows, None)
    out = rr.named_by_a_derived_layer(CH, tmp_path, 1000, "PAX3", reader, CODING)
    assert len(out["genes_named"]) == 8  # the committed field is NOT widened
    assert out["genes_named_total"] == 12
    assert out["genes_named_truncated_to"] == rr.GENES_NAMED_TRUNCATION == 8
    assert out["window"]["genes_at_the_bar_total"] == 12
    assert out["window"]["genes_at_the_bar_and_not_a_head_total"] == 0


def test_an_element_the_cache_does_not_hold_is_counted_by_name(monkeypatch, tmp_path):
    rows = [_el("E1", "EPHA4"), _el("E2", "EPHA4")]
    reader = _Reader({"E1": {"genes": [_g("EPHA4", -0.8, 0.05)]}, "E2": None})
    _patch(monkeypatch, rows, None)
    out = rr.named_by_a_derived_layer(CH, tmp_path, 1000, "PAX3", reader, CODING)
    assert out["window"]["elements_not_cached"] == 1


def test_no_scored_element_keeps_its_pending_reason(monkeypatch, tmp_path):
    _patch(monkeypatch, [], None)
    out = rr.named_by_a_derived_layer(CH, tmp_path, 1000, "PAX3", _Reader({}), CODING)
    assert out["pending"] == "no deletion has been scored within 100 kb of the recipient"
    assert out["names_the_recipient_in_the_window"] is False


def test_an_any_gene_head_disagreement_stops_the_wave(monkeypatch, tmp_path):
    """PLANT: the element's head is not its window's strongest gene."""
    rows = [_el("E1", "PAX3")]
    reader = _Reader({"E1": {"genes": [_g("EPHA4", -0.9, 0.05), _g("PAX3", -0.1, 0.02)]}})
    _patch(monkeypatch, rows, None)
    with pytest.raises(ot.InvariantRefutedError, match="THE WAVE STOPS"):
        rr.named_by_a_derived_layer(CH, tmp_path, 1000, "PAX3", reader, CODING)
