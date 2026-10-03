# SPDX-License-Identifier: AGPL-3.0-or-later
"""decompile.py shows the elements whose window names a gene, not only those whose head does.

Wave 2 of section 5 item 10, under the registration at `genomeos/attribution/onetarget2.py`
(24adf33). decompile.py's registered falsifier: "a gene whose decompiled element set is empty today
and non-empty under the window. That is the strongest form of this defect -- a gene the decompiler
says has no scored element when the sweep scored one at the bar for it -- and one instance
establishes it; none means the selection never excluded anything."

THE REGISTERED LIMIT IS LOAD-BEARING HERE: `_elements` selects on `predicted_coding`, so this module
depends on the CODING head specifically, and the measured invariant is an ANY-GENE one. That is why
`window_elements` REFUSES to form an any-gene window beside a coding head: comparing the two would
be comparing different populations, which is the defect the whole wave exists to stop.
"""

from __future__ import annotations

import json

import pytest

from genomeos import decompile as dc


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


CODING = {"APP", "CYYR1", "SOD1", "RUNX1"}


def _results(tmp_path):
    """Two sampled runs. E1's coding head is APP and its window also names CYYR1; E2's head is SOD1
    and its window names only SOD1; E3 is not in the response cache at all."""
    constrained = [
        {
            "id": "E1",
            "start": 1000,
            "end": 1200,
            "predicted_coding": {"gene": "APP", "log2_fold_change": -0.8},
        },
        {
            "id": "E3",
            "start": 9000,
            "end": 9200,
            "predicted_coding": {"gene": "APP", "log2_fold_change": -0.3},
        },
    ]
    uniform = [
        {
            "id": "E2",
            "start": 5000,
            "end": 5200,
            "predicted_coding": {"gene": "SOD1", "log2_fold_change": -0.5},
        },
    ]
    (tmp_path / "constrained_targets_chr21.json").write_text(json.dumps({"elements": constrained}))
    (tmp_path / "enhancer_targets_chr21.json").write_text(json.dumps({"elements": uniform}))
    return tmp_path


READER = _Reader(
    {
        "E1": {"genes": [_g("APP", -0.8, 0.05), _g("CYYR1", -0.1, 0.4)]},
        "E2": {"genes": [_g("SOD1", -0.5, 0.01)]},
        "E3": None,
    }
)


def test_without_a_reader_it_is_the_committed_function(tmp_path):
    d = _results(tmp_path)
    for symbol in ("APP", "SOD1", "CYYR1"):
        assert dc.window_elements(symbol, "chr21", d) == dc._elements(symbol, "chr21", d)


def test_a_reader_without_the_coding_set_refuses_rather_than_widening(tmp_path):
    """PLANT: asking for a window without the coding symbols. The match is on 'two different
    populations', a phrase only this refusal uses."""
    d = _results(tmp_path)
    with pytest.raises(ValueError, match="two different populations"):
        dc.window_elements("APP", "chr21", d, READER)


def test_a_gene_the_projection_hid_entirely_gains_its_element(tmp_path):
    """The registered falsifier's strongest form: CYYR1 has NO decompiled element today."""
    d = _results(tmp_path)
    assert dc._elements("CYYR1", "chr21", d) == []
    got = dc.window_elements("CYYR1", "chr21", d, READER, CODING)
    assert [e["id"] for e in got] == ["E1"]
    e = got[0]
    assert e["source"] == dc.SOURCE_WINDOW
    assert e["head_gene"] == "APP" and e["rank_in_the_coding_window"] == 2
    assert e["at_bar_log2_fold_change"] == 0.4


def test_an_appended_row_carries_no_head_derived_field(tmp_path):
    """`action`, `log2_fold_change`, `tissue`, `confidence` and `certainty` belong to the head's
    chosen gene and say nothing about a second mover, so they are ABSENT rather than filled."""
    d = _results(tmp_path)
    e = dc.window_elements("CYYR1", "chr21", d, READER, CODING)[0]
    for k in ("action", "log2_fold_change", "tissue", "confidence", "certainty"):
        assert k not in e, k


def test_a_gene_whose_head_elements_already_hold_it_gains_nothing(tmp_path):
    d = _results(tmp_path)
    head = dc._elements("APP", "chr21", d)
    assert [e["id"] for e in head] == ["E1", "E3"]
    got = dc.window_elements("APP", "chr21", d, READER, CODING)
    assert got == head  # APP is its own head at both, and E3 is not cached


def test_a_gene_named_by_nobody_is_still_empty(tmp_path):
    d = _results(tmp_path)
    assert dc.window_elements("RUNX1", "chr21", d, READER, CODING) == []


def test_the_registered_limit_is_stated_in_the_module():
    assert "coding head" in dc.CODING_HEAD_IS_NOT_INVARIANT
    assert "4,794" in dc.CODING_HEAD_IS_NOT_INVARIANT
