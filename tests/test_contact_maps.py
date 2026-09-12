# SPDX-License-Identifier: AGPL-3.0-or-later
"""Insulation minima from a contact map and their comparison with inferred boundaries (no network)."""

import numpy as np

from genomeos.predict.contact_maps import compare, insulation, local_minima


def test_insulation_dips_at_a_boundary_and_compare_counts():
    n = 200
    m = np.zeros((n, n, 2), dtype=np.float32)
    # two domains, 0..100 and 100..200, with contacts inside and none across
    m[:100, :100, :] = 1.0
    m[100:, 100:, :] = 1.0
    prof = insulation(m, w=10)
    assert prof[50] == 1.0 and prof[100] == 0.0  # across the boundary nothing touches
    mins = local_minima(prof, min_depth=0.3, spacing=5)
    assert mins == [100]
    cmp = compare([205_824, 50_000], [210_824, 900_000], tolerance=20_000)
    assert cmp["inferred_on_a_predicted_boundary"] == 1 and cmp["fraction_inferred_supported"] == 0.5
    assert cmp["predicted_on_an_inferred_boundary"] == 1 and cmp["fraction_predicted_matched"] == 0.5
    assert compare([205_824], [210_824])["median_distance_inferred_to_predicted"] == 5_000
