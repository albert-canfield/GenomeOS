# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""N1b's frozen code, on made-up inputs only, plus two checks against the downloaded file's obs fields.

Every test below runs before any byte of `.X` is read. The three that need the git-ignored store carry
`needs_local_data` and SKIP BY NAME where it is absent, because CI runs on a fresh checkout; none of
them is a bare assert on a path.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from genomeos.attribution import n1b_calibration as n1b

ROOT = Path(__file__).resolve().parent.parent
NORM_IDENTITIES = "data/cache/n1/K562_gwps_normalized_bulk_01.identities.json"
RAW_FILE = "data/cache/n1/K562_gwps_raw_bulk_01.h5ad"

#: INJECTED, not read live. A live RSS reading is valid only in a process that has done nothing
#: substantial first, and a guard that reads a process-wide figure measures the test process rather
#: than the work -- which is why one such guard is being removed from this project. The real run
#: records its own live reading beside this caveat.
RSS_CEILING_BYTES = 512 * 1024 * 1024


# --- the registered constants, pinned ----------------------------------------------------------------


def test_the_registered_bands_are_the_ones_the_study_was_opened_with() -> None:
    assert (n1b.VAR_LO, n1b.VAR_HI) == (0.80, 1.25)
    assert (n1b.TAIL_LO, n1b.TAIL_HI) == (0.04, 0.06)
    assert n1b.TAIL_Z == 1.96
    assert n1b.STRATA == 5
    assert n1b.NOMINAL_VAR == 1.0
    assert math.isclose(n1b.NOMINAL_TAIL, 0.04999579029644087, rel_tol=1e-12)


def test_an_inconclusive_reading_says_the_data_cannot_tell_and_never_no_information() -> None:
    assert n1b.INCONCLUSIVE_WORDING == "the data cannot tell"
    src = (ROOT / "genomeos/attribution/n1b_calibration.py").read_text()
    assert "no information" not in src


def test_the_module_reads_nothing_of_n1_but_the_file_identity() -> None:
    src = (ROOT / "genomeos/attribution/n1b_calibration.py").read_text()
    assert "from genomeos.attribution.n1_perturb_response import SOURCE" in src
    for forbidden in ("run_v2", "run_v3", "calibration_gate", "pooled_t", "candidate_ledger", "RESPONDS_T"):
        assert forbidden not in src, f"{forbidden} would couple N1b to N1's closed study"


def test_the_module_takes_no_live_memory_reading() -> None:
    src = (ROOT / "genomeos/attribution/n1b_calibration.py").read_text()
    for forbidden in ("getrusage", "psutil", "ru_maxrss", "/proc/self/status"):
        assert forbidden not in src
    assert RSS_CEILING_BYTES == 512 * 1024 * 1024


# --- the usable rows ---------------------------------------------------------------------------------


def test_a_row_is_dropped_only_for_the_cell_count_the_scale_factor_needs() -> None:
    got = n1b.usable_rows([10.0, float("nan"), 0.0, 3.0, float("inf")])
    assert got["rows_used"] == [0, 3]
    assert [d["row"] for d in got["rows_dropped"]] == [1, 2, 4]
    assert "not finite" in got["rows_dropped"][0]["reason"]
    assert "below" in got["rows_dropped"][1]["reason"]
    assert got["rows_given"] == 5


# --- the strata --------------------------------------------------------------------------------------


def test_quintiles_are_balanced_and_stratum_zero_is_the_least_expressed() -> None:
    expr = [7.0, 0.1, 3.0, 0.2, 5.0, 0.3, 9.0, 0.4, 1.0, 2.0, 11.0]
    st = n1b.expression_quintiles(expr)
    sizes = [st.count(s) for s in range(5)]
    assert max(sizes) - min(sizes) <= 1
    assert st[expr.index(min(expr))] == 0
    assert st[expr.index(max(expr))] == 4


def test_no_gene_has_a_stratum_when_there_are_no_genes() -> None:
    assert n1b.expression_quintiles([]) == []


# --- T itself ----------------------------------------------------------------------------------------


def test_t_is_x_times_the_square_root_of_the_row_cell_count() -> None:
    t = n1b.t_from_rows([[1.0, -2.0], [0.5, 0.25]], [100.0, 4.0])
    assert np.allclose(t, [[10.0, -20.0], [1.0, 0.5]])


def test_t_refuses_a_cell_count_that_does_not_match_the_rows() -> None:
    with pytest.raises(ValueError, match="cell counts"):
        n1b.t_from_rows([[1.0], [2.0]], [3.0])


def test_a_gene_with_a_non_finite_value_or_expression_cannot_be_calibrated() -> None:
    t = np.array([[1.0, 2.0, 3.0, 4.0], [1.0, float("nan"), 3.0, 4.0]])
    assert n1b.uncalibratable_genes(t, [1.0, 1.0, float("nan"), 1.0]) == [1, 2]


# --- the two figures and the bands -------------------------------------------------------------------


def _calibrated(rows: int = 60, genes: int = 400, seed: int = 7) -> np.ndarray:
    """T drawn exactly standard normal, by way of X and n: X = z / sqrt(n), so T = X * sqrt(n) = z."""
    rng = np.random.default_rng(seed)
    n = rng.integers(11, 562, size=rows).astype(float)
    z = rng.standard_normal((rows, genes))
    x = z / np.sqrt(n)[:, None]
    return n1b.t_from_rows(x, n)


def test_a_standard_normal_t_passes_both_bands_overall_and_in_every_quintile() -> None:
    t = _calibrated()
    fig = n1b.gene_figures(t)
    strata = n1b.expression_quintiles(list(range(t.shape[1])))
    scopes = {"overall": list(range(t.shape[1]))}
    recs = []
    for s in range(5):
        cols = [j for j, st in enumerate(strata) if st == s]
        recs.append({"stratum": s, **n1b.summarise(fig, cols, t)})
    overall = n1b.summarise(fig, scopes["overall"], t)
    verdict = n1b.decide(overall, recs)
    assert verdict["passed"], verdict["failed_scopes"]
    assert n1b.VAR_LO <= overall["var_primary"] <= n1b.VAR_HI
    assert n1b.TAIL_LO <= overall["tail_primary"] <= n1b.TAIL_HI
    assert "calibrated" in verdict["reading"]


def test_an_inflated_t_fails_the_variance_band_and_the_reading_calls_it_unusable() -> None:
    t = _calibrated() * 1.6
    fig = n1b.gene_figures(t)
    overall = n1b.summarise(fig, range(t.shape[1]), t)
    assert overall["var_primary"] > n1b.VAR_HI
    verdict = n1b.decide(overall, [])
    assert not verdict["passed"]
    assert "UNUSABLE" in verdict["reading"]
    assert "Albert's decision" in verdict["reading"]


def test_a_deflated_t_fails_the_variance_band() -> None:
    t = _calibrated() * 0.7
    overall = n1b.summarise(n1b.gene_figures(t), range(t.shape[1]), t)
    assert overall["var_primary"] < n1b.VAR_LO
    assert not n1b.decide(overall, [])["passed"]


def test_the_two_bands_are_independent_a_unit_variance_t_can_still_miss_the_tail() -> None:
    """T = +/-1 alternating: the variance is m/(m-1) -- the registered denominator is ddof=1, so a
    balanced two-point T gives 40/39 here and 514/513 = 1.00195 on the real rows, which is inside the
    band and is a bias the band carries, not a figure that moves it -- and the tail share is exactly 0.
    A rule that read the variance alone would call this calibrated, and the normal approximation is
    plainly false here."""
    m = 40
    t = np.where((np.add.outer(np.arange(m), np.arange(200)) % 2) == 0, 1.0, -1.0)
    overall = n1b.summarise(n1b.gene_figures(t), range(t.shape[1]), t)
    assert math.isclose(overall["var_primary"], m / (m - 1), rel_tol=1e-12)
    assert n1b.VAR_LO <= overall["var_primary"] <= n1b.VAR_HI
    assert overall["tail_primary"] == 0.0
    verdict = n1b.decide(overall, [])
    assert not verdict["passed"]
    assert verdict["failed_scopes"] == [{"scope": "overall", "legs": ["tail"]}]


def test_the_variance_can_miss_while_the_tail_holds_so_neither_leg_is_redundant() -> None:
    """A gene whose T is 0 on 95% of rows and +/-5.0966 on 5%: variance 1.299 (outside 1.25) while the
    tail share is exactly 0.05 (inside the band)."""
    rows, genes = 200, 100
    t = np.zeros((rows, genes))
    spike = math.sqrt(1.2987 / 0.05)
    for j in range(genes):
        hit = (np.arange(rows) + j) % 20 == 0  # exactly 10 of 200 rows, i.e. 5%
        t[hit, j] = spike * np.where(np.arange(hit.sum()) % 2 == 0, 1.0, -1.0)
    overall = n1b.summarise(n1b.gene_figures(t), range(genes), t)
    assert overall["var_primary"] > n1b.VAR_HI
    assert n1b.TAIL_LO <= overall["tail_primary"] <= n1b.TAIL_HI
    assert n1b.decide(overall, [])["failed_scopes"] == [{"scope": "overall", "legs": ["variance"]}]


def test_one_failing_quintile_fails_the_whole_rule_so_an_overall_pass_cannot_hide_it() -> None:
    good = {"var_primary": 1.0, "tail_primary": 0.05}
    strata = [{"stratum": s, **good} for s in range(5)]
    strata[0] = {"stratum": 0, "var_primary": 0.31, "tail_primary": 0.004}
    verdict = n1b.decide(good, strata)
    assert not verdict["passed"]
    assert verdict["failed_scopes"] == [{"scope": "quintile 0", "legs": ["variance", "tail"]}]
    assert verdict["legs"][0]["scope"] == "overall"
    assert verdict["legs"][0]["variance"]["in_band"]


def test_a_figure_that_cannot_be_formed_fails_its_band_and_says_why() -> None:
    v = n1b.band_verdict(None, n1b.VAR_LO, n1b.VAR_HI)
    assert not v["in_band"]
    assert v["reason"] == "the figure could not be formed"


# --- every ratio beside its absolute level and its control -------------------------------------------


def test_a_ratio_is_never_emitted_without_its_absolute_level_and_its_control() -> None:
    r = n1b.ratio_report(0.06, n1b.NOMINAL_TAIL, "share of |T| > 1.96")
    assert set(r) == {"name", "absolute", "control", "control_kind", "ratio_to_control"}
    assert r["absolute"] == 0.06
    assert r["control"] == n1b.NOMINAL_TAIL
    assert "not an empirical control group" in r["control_kind"]
    assert math.isclose(r["ratio_to_control"], 0.06 / n1b.NOMINAL_TAIL)


def test_every_summarised_scope_carries_both_ratios_with_their_levels() -> None:
    t = _calibrated(rows=20, genes=50)
    rec = n1b.summarise(n1b.gene_figures(t), range(50), t)
    assert len(rec["ratios"]) == 2
    for r in rec["ratios"]:
        assert r["absolute"] is not None and r["control"] is not None
    assert rec["var_pooled_over_values"] is not None
    assert rec["second_moment_mean_over_genes"] is not None
    assert rec["var_mean_over_genes"] is not None


def test_an_empty_scope_reports_no_figure_rather_than_a_zero() -> None:
    t = _calibrated(rows=10, genes=20)
    rec = n1b.summarise(n1b.gene_figures(t), [], t)
    assert rec["var_primary"] is None and rec["tail_primary"] is None
    assert rec["reason"] == "no gene in this scope"


# --- the intervals -----------------------------------------------------------------------------------


def test_a_constant_bootstrap_is_degenerate_decides_nothing_and_routes_to_a_bound() -> None:
    iv = n1b.bootstrap_interval([0.0] * 500, 0.0, degenerate_fallback=n1b.wilson_interval(0, 514))
    assert iv["identical_resample_share"] == 1.0
    assert iv["degenerate"] is True
    assert iv["decides"] == "nothing: every resample returned the identical value"
    assert iv["fallback"]["kind"] == "wilson"
    assert iv["fallback"]["low"] == 0.0
    assert 0.0 < iv["fallback"]["high"] < 0.02


def test_the_wilson_bound_counts_clusters_and_never_the_row_by_gene_values() -> None:
    by_cluster = n1b.wilson_interval(0, 514)
    by_value = n1b.wilson_interval(0, 514 * 8248)
    assert by_cluster["clusters"] == 514
    assert by_cluster["high"] > by_value["high"] * 100
    assert "not on the (row, gene) share" in by_cluster["note"]


def test_a_wilson_bound_needs_a_cluster() -> None:
    assert n1b.wilson_interval(0, 0)["low"] is None


def test_every_interval_reports_its_identical_resample_share() -> None:
    iv = n1b.bootstrap_interval([0.1, 0.2, 0.3, 0.2], 0.2)
    assert iv["identical_resample_share"] == 0.5
    assert iv["low"] <= 0.2 <= iv["high"]
    assert n1b.bootstrap_interval([], 1.0)["identical_resample_share"] is None


def test_the_cluster_bootstrap_resamples_rows_and_brackets_the_point_estimate() -> None:
    t = _calibrated(rows=40, genes=120, seed=11)
    fig = n1b.gene_figures(t)
    scopes = {"overall": list(range(120))}
    points = {"overall": n1b.summarise(fig, scopes["overall"], t)}
    iv = n1b.cluster_bootstrap(t, scopes, points)["overall"]
    assert iv["cluster"] == "non-targeting row" and iv["clusters"] == 40
    for leg in ("variance", "tail"):
        assert iv[leg]["low"] <= points["overall"][f"{leg[:4] if leg == 'tail' else 'var'}_primary"]
        assert iv[leg]["identical_resample_share"] is not None
    assert any("UNCONFIRMED" in a for a in iv["assumptions"])
    assert len(iv["assumptions"]) == 4


def test_a_zero_exceedance_scope_routes_its_tail_interval_to_the_wilson_bound() -> None:
    t = np.where((np.add.outer(np.arange(30), np.arange(60)) % 2) == 0, 1.0, -1.0)
    fig = n1b.gene_figures(t)
    scopes = {"overall": list(range(60))}
    points = {"overall": n1b.summarise(fig, scopes["overall"], t)}
    assert points["overall"]["rows_with_any_exceedance"] == 0
    iv = n1b.cluster_bootstrap(t, scopes, points)["overall"]["tail"]
    assert iv["identical_resample_share"] == 1.0
    assert iv["fallback"]["successes"] == 0 and iv["fallback"]["clusters"] == 30


def test_the_variance_has_no_exact_fallback_and_the_module_says_the_data_cannot_tell() -> None:
    t = np.ones((5, 4)) * np.array([1.0, -1.0, 1.0, -1.0])[None, :]
    fig = n1b.gene_figures(t)
    scopes = {"overall": [0, 1, 2, 3]}
    points = {"overall": n1b.summarise(fig, scopes["overall"], t)}
    iv = n1b.cluster_bootstrap(t, scopes, points)["overall"]["variance"]
    assert iv["identical_resample_share"] == 1.0
    assert n1b.INCONCLUSIVE_WORDING in iv["fallback"]["note"]


def test_a_variance_about_the_genes_own_mean_needs_two_rows() -> None:
    with pytest.raises(ValueError, match="at least two rows"):
        n1b.gene_figures(np.zeros((1, 5)))


# --- the layout probe --------------------------------------------------------------------------------


class _FakeDatasetId:
    def __init__(self, offset):
        self._offset = offset

    def get_offset(self):
        return self._offset


class _FakeDataset:
    def __init__(self, offset, shape, dtype):
        self.id = _FakeDatasetId(offset)
        self.shape = shape
        self.dtype = np.dtype(dtype)


def test_the_layout_probe_gives_the_row_stride_of_a_contiguous_float32_matrix() -> None:
    got = n1b.x_layout({"X": _FakeDataset(2048, (11258, 8248), "float32")})
    assert got == {
        "dataset": "X",
        "offset": 2048,
        "rows": 11258,
        "cols": 8248,
        "dtype": "float32",
        "row_stride_bytes": 32992,
        "layout": "contiguous",
    }


def test_the_layout_probe_refuses_a_chunked_matrix_because_a_row_has_no_single_range() -> None:
    with pytest.raises(RuntimeError, match="not contiguous"):
        n1b.x_layout({"X": _FakeDataset(None, (10, 10), "float32")})


def test_the_layout_probe_refuses_a_width_that_would_make_the_stride_wrong() -> None:
    with pytest.raises(RuntimeError, match="not float32"):
        n1b.x_layout({"X": _FakeDataset(2048, (10, 10), "float64")})


def test_the_whole_file_md5_is_declared_unverifiable_rather_than_claimed() -> None:
    assert "NOT verified" in n1b.MD5_NOT_VERIFIABLE
    assert "374,587,922" in n1b.MD5_NOT_VERIFIABLE
    assert n1b.source_file("normalized")["md5"] == "a3dfaa94ea8724217f5ecb1e14a5f0c8"
    assert n1b.source_file("normalized")["bytes"] == 374_587_922


# --- the range log -----------------------------------------------------------------------------------


def test_the_range_reader_logs_a_url_a_range_and_a_sha256_for_every_read() -> None:
    rf = n1b.RangeFile("https://example.invalid/f", 100)
    payload = b"abcdefghij"

    def fake_get(start, end, purpose):
        rf.fetched += end - start + 1
        rf.log.append(
            {
                "url": rf.url,
                "range": f"bytes={start}-{end}",
                "bytes": end - start + 1,
                "sha256": __import__("hashlib").sha256(payload).hexdigest(),
                "purpose": purpose,
            }
        )
        return payload

    rf._get = fake_get  # noqa: SLF001
    got = rf.fetch_rows(2048, 10, [0, 3], "X")
    assert sorted(got) == [0, 3]
    assert [e["range"] for e in rf.log] == ["bytes=2048-2057", "bytes=2078-2087"]
    for entry in rf.log:
        assert entry["url"] and entry["sha256"] and entry["bytes"] == 10 and "row" in entry["purpose"]
    assert rf.fetched == 20


# --- the downloaded store: the row counts, measured ---------------------------------------------------


@pytest.mark.needs_local_data(NORM_IDENTITIES, how="N1 fetched the normalized file's index by range")
def test_the_non_targeting_row_count_is_measured_from_the_files_own_index() -> None:
    idx = json.loads((ROOT / NORM_IDENTITIES).read_text())["obs_index"]
    assert len(idx) == 11258
    nt = [i for i, s in enumerate(idx) if "non-targeting" in s]
    assert len(nt) == 585
    assert all(s.split("_")[1:] == ["non-targeting"] * 3 for s in (idx[i] for i in nt))


@pytest.mark.needs_local_data(RAW_FILE, how="N1 downloaded the raw pseudobulk file into data/cache/n1")
def test_the_core_set_and_the_unusable_rows_are_the_same_71_rows() -> None:
    h5py = pytest.importorskip("h5py")
    idx = json.loads((ROOT / NORM_IDENTITIES).read_text())["obs_index"]
    nt = np.array([i for i, s in enumerate(idx) if "non-targeting" in s])
    with h5py.File(ROOT / RAW_FILE, "r") as h:
        core = h["obs/core_control"][()]
        n = h["obs/num_cells_filtered"][()]
        ce = h["obs/control_expr"][()]
    assert int(core.sum()) == 514
    assert int(core[nt].sum()) == 514, "core_control is true nowhere outside the non-targeting rows"
    got = n1b.usable_rows([float(v) for v in n[nt]])
    assert len(got["rows_used"]) == 514
    assert len(got["rows_dropped"]) == 71
    assert set(nt[got["rows_used"]].tolist()) == set(nt[core[nt]].tolist())
    assert not np.isfinite(ce[nt]).any(), "obs/control_expr is NaN on every non-targeting row"


# --- what a pass licenses, and what cannot be established ---------------------------------------------


def test_a_verdict_cannot_be_emitted_without_the_pass_licence_attached() -> None:
    """The caveat has to TRAVEL with the figures. A reader who sees only `decide`'s output must still
    see that a pass is about non-targeting subsets of the normalisation reference."""
    for passing in (True, False):
        overall = {"var_primary": 1.0 if passing else 2.0, "tail_primary": 0.05}
        v = n1b.decide(overall, [])
        assert v["passed"] is passing
        assert v["pass_licence"] == n1b.PASS_LICENCE
        assert v["cannot_establish"] == list(n1b.CANNOT_ESTABLISH)


def test_the_pass_licence_names_the_assumption_and_the_step_that_must_restate_it() -> None:
    lic = n1b.PASS_LICENCE
    assert "non-targeting subsets of the normalisation reference" in lic
    assert "ASSUMPTION" in lic
    assert "step 2" in lic and "beside every perturbed-row figure" in lic


def test_a_passing_reading_states_the_licence_in_its_own_sentence() -> None:
    v = n1b.decide({"var_primary": 1.0, "tail_primary": 0.05}, [])
    assert v["passed"]
    assert n1b.PASS_LICENCE in v["reading"]


def test_what_cannot_be_established_includes_the_perturbed_rows_and_the_unconfirmed_divisor() -> None:
    joined = " ".join(n1b.CANNOT_ESTABLISH)
    assert "UNCONFIRMED" in joined
    assert "PER-FACTOR COUNTS ARE NOT READ IN N1b, AT ALL" in joined
    assert "fitness phenotype and not metadata" in joined
    assert "No factor row is read on a pass or a fail" in joined
    assert n1b.PASS_LICENCE in n1b.CANNOT_ESTABLISH


def test_the_licence_moves_no_band() -> None:
    """The addition bounds what a pass LICENSES, not what counts as one. The bands are the same
    objects they were when the study was frozen, and a figure outside one still fails."""
    assert (n1b.VAR_LO, n1b.VAR_HI, n1b.TAIL_LO, n1b.TAIL_HI, n1b.STRATA) == (0.80, 1.25, 0.04, 0.06, 5)
    assert not n1b.decide({"var_primary": 1.26, "tail_primary": 0.05}, [])["passed"]
    assert not n1b.decide({"var_primary": 1.0, "tail_primary": 0.0601}, [])["passed"]


# --- the cell-count distribution: quantiles only, no identity -------------------------------------------


def test_the_cell_count_distribution_persists_no_row_and_no_index() -> None:
    from scripts import n1b_fetch

    counts = [float(v) for v in range(10, 110)]
    got = n1b_fetch.cell_count_distribution(counts, list(range(100)))
    assert set(got) == {
        "rows",
        "with_a_finite_count",
        "quantiles",
        "mean",
        "cells_in_total",
        "what_is_not_here",
    }
    assert got["rows"] == 100 and got["with_a_finite_count"] == 100
    assert sorted(got["quantiles"]) == sorted(str(q) for q in n1b_fetch.CELL_COUNT_QUANTILES)
    assert got["quantiles"]["0.0"] == 10.0 and got["quantiles"]["1.0"] == 109.0
    assert "no per-row count and no row index" in got["what_is_not_here"]


def test_the_distribution_drops_non_finite_counts_and_says_how_many_it_kept() -> None:
    from scripts import n1b_fetch

    got = n1b_fetch.cell_count_distribution([1.0, float("nan"), 3.0], [0, 1, 2])
    assert got["rows"] == 3 and got["with_a_finite_count"] == 2
    assert got["quantiles"]["0.5"] == 2.0


def test_a_scope_with_no_finite_count_reports_no_quantiles_rather_than_a_zero() -> None:
    from scripts import n1b_fetch

    got = n1b_fetch.cell_count_distribution([float("nan")], [0])
    assert got["quantiles"] is None and got["with_a_finite_count"] == 0


def test_the_disclosure_wording_is_the_registered_one() -> None:
    from scripts import n1b_fetch

    assert n1b_fetch.DISCLOSURE == "perturbed-row cell-count distribution read, no identity, no expression"


@pytest.mark.needs_local_data(
    "data/cache/n1b/probe.json", how="scripts/n1b_fetch.py probe fetches the file's index by range"
)
def test_the_probe_stores_a_perturbed_distribution_and_no_perturbed_row() -> None:
    probe = json.loads((ROOT / "data/cache/n1b/probe.json").read_text())
    cells = probe["cell_counts"]
    assert cells["disclosure"] == "perturbed-row cell-count distribution read, no identity, no expression"
    assert cells["perturbed"]["rows"] == probe["rows_total"] - probe["non_targeting_count"]
    assert "per_row" not in json.dumps(cells)
    # the only per-row obs the probe keeps is for the NON-TARGETING rows, and it keeps 585 of each
    for column, values in probe["obs"].items():
        assert len(values) == probe["non_targeting_count"], column
