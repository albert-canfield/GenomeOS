# SPDX-License-Identifier: AGPL-3.0-or-later
"""Predicted RNA runs as a filter on parser candidates (no network)."""

from genomeos.coords import Strand
from genomeos.genome.segments import Prediction
from genomeos.predict.rna_tracks import RnaCoverage, filter_predictions


def test_covered_fraction_and_filter(tmp_path):
    cov = RnaCoverage("t", 10_000, cache=tmp_path)
    cov.runs = {"+": [(100, 200), (300, 400)], "-": [(1000, 2000)]}
    cov._starts = {s: [a for a, _ in r] for s, r in cov.runs.items()}
    assert cov.covered_fraction("+", 100, 200) == 1.0
    assert cov.covered_fraction("+", 150, 350) == 0.5  # 50 of 200 bases in each run
    assert cov.covered_fraction("+", 500, 600) == 0.0 and cov.covered_fraction("-", 1500, 1600) == 1.0
    made = Prediction("t", Strand.PLUS, 100, 400, [(100, 200), (300, 400)], 10.0)
    silent = Prediction("t", Strand.PLUS, 500, 700, [(500, 600), (650, 700)], 10.0)
    minus = Prediction("t", Strand.MINUS, 1000, 2000, [(1000, 1100), (1900, 2000)], 10.0)
    kept = filter_predictions([made, silent, minus], cov)
    assert [p.start for p in kept] == [100, 1000]
