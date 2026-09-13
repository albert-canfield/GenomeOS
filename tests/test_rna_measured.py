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


def test_orient_swaps_tracks_labelled_by_the_read(monkeypatch):
    from genomeos.genome import rna_measured

    tracks = {"tracks": {"+": {"accession": "P", "href": "plus"}, "-": {"accession": "M", "href": "minus"}}}
    m = MeasuredRna("IMR-90", "t", tracks=tracks)
    exons = {"+": [(100, 200), (300, 400)], "-": [(500, 600)]}

    def covered_by_strand(href, chrom, ex, signal, progress=None, with_bytes=False):
        # the file named "plus" covers the minus-strand exons (500..600), "minus" covers the plus ones
        out = [1.0 if ((a >= 500) == (href == "plus")) else 0.0 for a, b in ex]
        return (out, 0) if with_bytes else out

    monkeypatch.setattr(rna_measured, "_covered", covered_by_strand)
    m.prepare(exons)
    assert m.orientation == {"IMR-90": "swapped"}
    assert m.summary()["tracks"]["IMR-90"] == {"+": "M", "-": "P"}
    assert m.covered_fraction("+", 100, 200) == 1.0 and m.covered_fraction("-", 500, 600) == 1.0


def test_orient_keeps_tracks_labelled_by_the_transcript(monkeypatch):
    from genomeos.genome import rna_measured

    def covered(href, chrom, ex, signal, progress=None, with_bytes=False):
        out = [1.0 if ((a < 500) == (href == "plus")) else 0.0 for a, b in ex]
        return (out, 0) if with_bytes else out

    monkeypatch.setattr(rna_measured, "_covered", covered)
    tracks = {"tracks": {"+": {"accession": "P", "href": "plus"}, "-": {"accession": "M", "href": "minus"}}}
    m = MeasuredRna("K562", "t", tracks=tracks)
    m.prepare({"+": [(100, 200)], "-": [(500, 600)]})
    assert m.orientation == {"K562": "as labelled"} and m.summary()["tracks"]["K562"] == {"+": "P", "-": "M"}
