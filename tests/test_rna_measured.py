# SPDX-License-Identifier: AGPL-3.0-or-later
"""Measured RNA as a filter on parser candidates (no network: the bigWig pass is stubbed)."""

from genomeos.coords import Strand
from genomeos.genome.rna_measured import MeasuredRna
from genomeos.genome.segments import Prediction


def test_filter_keeps_candidates_with_signal_on_their_exons():
    m = MeasuredRna("K562", "t", tracks={"tracks": {"+": {"accession": "a"}, "-": {"accession": "b"}}})

    def fake_prepare(exons_by_strand, progress=None):
        for strand, exons in exons_by_strand.items():
            for a, b in exons:
                m.fractions[(strand, a, b)] = 1.0 if strand == "+" and a < 300 else 0.0

    m.prepare = fake_prepare
    made = Prediction("t", Strand.PLUS, 100, 250, [(100, 150), (200, 250)], 1.0)
    silent = Prediction("t", Strand.PLUS, 400, 500, [(400, 450), (460, 500)], 1.0)
    minus = Prediction("t", Strand.MINUS, 100, 250, [(100, 150)], 1.0)
    kept = m.filter([made, silent, minus], 0.5)
    assert [(p.strand.value, p.start) for p in kept] == [("+", 100)]
    assert m.summary()["exons_measured"] == 5
