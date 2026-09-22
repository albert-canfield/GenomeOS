"""The statistics behind the depth banding of the reader's per-biosample family.

The claim these support is that `enhancers_active` is substantially assay depth and that
what survives depth replicates across disjoint halves of the same genome, so the helpers
have to be right about rank correlation, least-squares residuals and the permutation
test before the banding means anything.
"""

import importlib.util
import os

_SPEC = importlib.util.spec_from_file_location(
    "reader_depth_family",
    os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts", "reader_depth_family.py"
    ),
)
rdf = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(rdf)


def test_spearman_is_rank_based_and_signed():
    x = [1, 2, 3, 4, 5]
    assert rdf.spearman(x, [10, 20, 30, 40, 50]) == 1.0
    assert rdf.spearman(x, [50, 40, 30, 20, 10]) == -1.0
    # monotone but wildly non-linear: Pearson would not see 1.0, Spearman must
    assert rdf.spearman(x, [1, 2, 4, 8, 10_000]) == 1.0
    # ties share the average rank rather than an arbitrary order
    assert rdf.spearman([1, 1, 2, 2], [5, 5, 9, 9]) == 1.0


def test_fit_residuals_recovers_the_line_and_leaves_nothing():
    x = [1.0, 2.0, 3.0, 4.0]
    a, b, resid = rdf.fit_residuals(x, [3.0, 5.0, 7.0, 9.0])
    assert round(a, 6) == 1.0 and round(b, 6) == 2.0
    assert all(abs(r) < 1e-9 for r in resid)
    # an offset point shows up in the residual and the residuals still sum to zero
    _, _, resid = rdf.fit_residuals(x, [3.0, 5.0, 7.0, 12.0])
    assert abs(sum(resid)) < 1e-9 and max(resid) == resid[-1]


def test_permutation_p_separates_a_real_association_from_a_shuffled_one():
    x = [float(i) for i in range(13)]
    y = [float(i) for i in range(13)]
    rho = rdf.spearman(x, y)
    # 13 labels: only the identity permutation reaches rho 1.0, so p is the floor
    assert rdf.permutation_p(x, y, rho, 2000, seed=1) < 0.01
    # a covariate that carries no order cannot beat its own observed rho often either way,
    # but a low bar (rho 0) must be met about half the time
    assert 0.3 < rdf.permutation_p(x, y, 0.0, 2000, seed=1) < 0.7


def test_nodes_open_style_median_split_cannot_distinguish_two_samples():
    """The defect this sweep found, stated as a property: calling a node open when it is at
    or above the median density *of that sample* returns half the nodes whatever the depth."""
    for scale in (1, 10, 100):
        densities = [i * scale for i in range(1, 1001)]
        median = sorted(densities)[len(densities) // 2]
        assert sum(1 for d in densities if d >= median) == 500
