# SPDX-License-Identifier: AGPL-3.0-or-later
"""variation.py reads the coding window beside the one coding gene it kept.

Wave 2 of section 5 item 10, under the registration at `genomeos/attribution/onetarget2.py`
(24adf33). variation.py's falsifier is registered as EXPECTED NOT TO FIRE: "any row's
log2_fold_change changing. The head's magnitude is the window maximum by construction, so this
number must NOT move; if it does, predict_target's size rule and genes_at_bar's are not the same
rule and every magnitude in the census is in doubt." A move is therefore a refutation that stops
the wave, never a finding, and `window_elements_of` raises rather than reporting it.

THE REGISTERED LIMIT IS LOAD-BEARING HERE and is asserted, not only written down: this module reads
`predicted_coding`, and the measured invariant is an ANY-GENE one. A coding-head disagreement is
expected at the 10-in-4,794 rate the 2026-09-27 result found, is counted by name, and must NOT
refuse -- there is a test below that plants one and asserts the function returns.
"""

from __future__ import annotations

import json

import pytest

from genomeos.attribution import onetarget2 as ot
from genomeos.attribution import variation


def _g(gene: str, drop: float, rise: float) -> dict:
    return {
        "gene": gene,
        "max_drop_log2fc": drop,
        "max_rise_log2fc": rise,
        "max_drop_tissue": "K562",
        "max_rise_tissue": "occipital lobe",
    }


class _Reader:
    """A response cache of exactly these records, with the real `genes_at_bar` behind it."""

    def __init__(self, records: dict[str, dict | None]) -> None:
        self.records = records

    def element(self, chrom, element_id):
        return self.records.get(element_id)

    def at_bar(self, chrom, element_id, min_effect=0.1, genes=None):
        from genomeos.attribution.targets import genes_at_bar

        got = genes_at_bar(self.element(chrom, element_id), min_effect)
        if got is None or genes is None:
            return got
        return [(g, v) for g, v in got if g in genes]


#: Varied on the axes the guards turn on: E1's coding window names a second coding gene; E2's names
#: only its head; E3 is not cached; E4 has NO coding head although its window names a coding gene,
#: which is the 2026-09-27 exception and must be counted rather than refused.
CODING = {"APP", "SOD1", "CYYR1", "RUNX1"}


def _results(tmp_path, elements):
    d = tmp_path
    (d / "constrained_targets_chr21.json").write_text(json.dumps({"elements": elements}))
    (d / "enhancer_targets_chr21.json").write_text(json.dumps({"elements": []}))
    return d


def _el(eid, start, gene, l2, coding_gene=True):
    e = {"id": eid, "start": start, "end": start + 200, "constrained_fraction": 0.3}
    if gene is not None:
        e["predicted"] = {"gene": gene, "log2_fold_change": l2}
        if coding_gene:
            e["predicted_coding"] = {"gene": gene, "log2_fold_change": l2}
    return e


def test_without_a_reader_the_rows_and_the_arm_are_the_committed_call(tmp_path):
    d = _results(tmp_path, [_el("E1", 1000, "APP", -0.8)])
    rows, arm = variation.window_elements_of("chr21", d)
    assert arm == {}
    assert rows == variation.elements_of("chr21", d)
    assert all("window" not in r for r in rows)


def test_with_a_reader_every_committed_field_is_identical_and_window_is_the_only_new_key(tmp_path):
    d = _results(tmp_path, [_el("E1", 1000, "APP", -0.8), _el("E2", 5000, "SOD1", -0.5)])
    reader = _Reader(
        {
            "E1": {"genes": [_g("APP", -0.8, 0.05), _g("CYYR1", -0.1, 0.4)]},
            "E2": {"genes": [_g("SOD1", -0.5, 0.01)]},
        }
    )
    off = variation.elements_of("chr21", d)
    rows, arm = variation.window_elements_of("chr21", d, reader, CODING)
    for a, b in zip(off, rows, strict=True):
        assert set(b) - set(a) == {"window"}
        assert all(a[k] == b[k] for k in a)
    assert arm["magnitudes_compared"] == 2
    assert arm["rows_gaining_a_coding_gene"] == 1 and arm["coding_genes_gained"] == 1
    assert rows[0]["window"]["extra"] == ["CYYR1"]
    assert rows[1]["window"]["extra"] == []
    assert "does NOT fire, as registered" in arm["verdict_against_its_own_falsifier"]


def test_an_element_the_cache_does_not_hold_is_counted_and_not_compared(tmp_path):
    d = _results(tmp_path, [_el("E3", 9000, "APP", -0.8)])
    rows, arm = variation.window_elements_of("chr21", d, _Reader({"E3": None}), CODING)
    assert arm["rows_not_cached"] == 1 and arm["magnitudes_compared"] == 0
    assert rows[0]["window"] == {"not_cached": True, "at_bar": [], "extra": []}


def test_a_magnitude_that_moved_refuses_and_names_the_row(tmp_path):
    """PLANT: the stored magnitude is not the window's value for the same gene.

    The sweep cannot write this, so it is planted. The match is on 'THE MAGNITUDE MOVED', a phrase
    only this guard's message uses.
    """
    d = _results(tmp_path, [_el("E1", 1000, "APP", -0.8)])
    reader = _Reader({"E1": {"genes": [_g("APP", -0.6, 0.05)]}})
    with pytest.raises(variation.MagnitudeMovedError, match="THE MAGNITUDE MOVED"):
        variation.window_elements_of("chr21", d, reader, CODING)


def test_a_stored_gene_absent_from_its_own_window_refuses_too(tmp_path):
    """PLANT: the head gene is not at the bar at all, so there is no value to compare."""
    d = _results(tmp_path, [_el("E1", 1000, "APP", -0.8)])
    reader = _Reader({"E1": {"genes": [_g("CYYR1", -0.4, 0.02)]}})
    with pytest.raises(variation.MagnitudeMovedError, match="THE MAGNITUDE MOVED"):
        variation.window_elements_of("chr21", d, reader, CODING)


def test_an_any_gene_head_disagreement_stops_the_wave(tmp_path):
    """PLANT: the any-gene head is not the window's strongest gene. One refusal, for one invariant.

    The row's coding head and magnitude still agree with the coding window, so this test turns on
    the ANY-GENE axis alone -- which is the axis `check_head_invariant` watches.
    """
    e = _el("E1", 1000, "APP", -0.8)
    e["predicted"] = {"gene": "APP", "log2_fold_change": -0.8}
    d = _results(tmp_path, [e])
    #: CYYR1 is the strongest gene in the window, so the any-gene head should be CYYR1 and is APP.
    reader = _Reader({"E1": {"genes": [_g("APP", -0.8, 0.05), _g("CYYR1", -0.9, 0.02)]}})
    with pytest.raises(ot.InvariantRefutedError, match="THE WAVE STOPS"):
        variation.window_elements_of("chr21", d, reader, CODING)


def test_a_coding_head_disagreement_is_counted_by_name_and_does_not_refuse(tmp_path):
    """THE REGISTERED LIMIT, asserted. An element with no `predicted_coding` where the window names a
    coding gene is the 2026-09-27 exception: 10 of 4,794. It is counted and the function RETURNS."""
    e = _el("E4", 1000, "APP", -0.8, coding_gene=False)  # predicted, but no predicted_coding
    d = _results(tmp_path, [e])
    reader = _Reader({"E4": {"genes": [_g("APP", -0.8, 0.05)]}})
    rows, arm = variation.window_elements_of("chr21", d, reader, CODING)
    assert arm["coding_head_disagreements"] == 1
    assert arm["any_gene_head_disagreements"] == 0
    assert arm["magnitudes_compared"] == 0  # there is no stored coding target to compare
    assert rows[0]["target"] is None
    assert "coding-head disagreement is expected" in arm["limit"]


def test_the_limit_and_the_falsifier_are_carried_in_the_arm(tmp_path):
    d = _results(tmp_path, [_el("E2", 5000, "SOD1", -0.5)])
    reader = _Reader({"E2": {"genes": [_g("SOD1", -0.5, 0.01)]}})
    _, arm = variation.window_elements_of("chr21", d, reader, CODING)
    assert arm["falsifier"] == variation.WINDOW_FALSIFIER
    assert arm["limit"] == variation.CODING_HEAD_IS_NOT_INVARIANT
