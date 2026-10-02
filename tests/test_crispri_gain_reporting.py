"""A gain is refused where the feature was never available, and an interval says how it was made.

The committed `crispri_published.json` of 2026-09-27 reports a deletion gain with an interval for three
strata whose `deletion_available` is false (HCT116 +0.0004, Jurkat +0.0038, WTC11 -0.0049). Those cells are
outside `MODEL_CELLS`, so no deletion value was ever read for them: the difference between the two models
there is an artefact of the model form. These tests fail on that behaviour.
"""

from genomeos.attribution import crispri


def test_an_unavailable_stratum_reports_no_number():
    out = crispri.gain_unavailable()
    assert out["gain"] is None
    assert out["ci95"] is None
    assert "not available" in out["unavailable"]


def test_the_gain_is_not_even_computed_where_the_feature_is_missing():
    """Not merely blanked afterwards: the computation must not run, so no number can leak through."""
    calls = []

    def gain() -> dict[str, object]:
        calls.append(1)
        return {"gain": -0.0049, "ci95": [-0.0174, 0.0], "resamples": 200}

    refused = crispri.gain_where_available(False, gain)
    assert calls == []
    assert refused["gain"] is None
    assert "unavailable" in refused

    allowed = crispri.gain_where_available(True, gain)
    assert calls == [1]
    assert allowed["gain"] == -0.0049


def test_an_available_stratum_keeps_its_number():
    assert crispri.gain_where_available(True, lambda: {"gain": 0.1361})["gain"] == 0.1361


def test_the_reason_names_the_artefact_rather_than_calling_it_zero():
    """ "No gain" and "a gain of about zero" are different claims; only the first is true here."""
    reason = crispri.UNAVAILABLE_GAIN
    assert "artefact" in reason
    assert "not a measurement" in reason


def test_an_interval_says_how_many_clusters_and_draws_it_rests_on():
    made = crispri._interval_provenance(clusters=24, requested=2000, kept=1850)
    assert made["clusters"] == 24
    assert made["draws_requested"] == 2000
    assert made["draws_dropped"] == 150
    assert made["resamples"] == 1850


def test_an_interval_below_the_minimum_says_so_rather_than_staying_silent():
    """200 draws put each tail on a handful of values, which the digits of a bound do not show."""
    assert crispri._interval_provenance(24, 200, 200)["met_minimum"] is False
    assert crispri._interval_provenance(24, 2000, 1500)["met_minimum"] is True


def test_the_default_resample_count_clears_the_stated_minimum():
    assert crispri.BOOTSTRAPS >= crispri.MIN_RESAMPLES
    assert crispri.MIN_RESAMPLES >= 1000


def test_the_estimators_report_their_provenance():
    """Both estimators, so neither can be quoted without saying what the interval rests on."""
    pairs = [
        crispri.Pair(
            chrom=f"chr{1 + i % 4}",
            start=i * 1000,
            end=i * 1000 + 300,
            gene=f"G{i}",
            cell="K562",
            dataset="test",
            distance=10_000 + i * 100,
            dhs=1.0,
            h3k27ac=1.0,
            regulated=(i % 3 == 0),
        )
        for i in range(24)
    ]
    a = [0.9 if p.regulated else 0.1 for p in pairs]
    b = [0.5 for _ in pairs]
    for out in (
        crispri.gain_interval(a, b, pairs, n=50),
        crispri.weighted_gain(a, b, pairs, weighted=False, n=50),
    ):
        for key in ("clusters", "draws_requested", "draws_dropped", "met_minimum", "resamples"):
            assert key in out, key
        assert out["draws_requested"] == 50
        assert out["met_minimum"] is False


# --- an interval needs enough clusters to mean 95% --------------------------------------------------------


def _pairs_on(chroms: int, per: int = 8) -> list:
    out = []
    for c in range(chroms):
        for i in range(per):
            out.append(
                crispri.Pair(
                    chrom=f"chr{c + 1}",
                    start=i * 1000,
                    end=i * 1000 + 300,
                    gene=f"G{c}_{i}",
                    cell="K562",
                    dataset="test",
                    distance=10_000 + i * 100,
                    dhs=1.0,
                    h3k27ac=1.0,
                    regulated=(i % 4 == 0),
                )
            )
    return out


def test_an_interval_from_too_few_clusters_is_withheld_not_reported():
    """GM12878 reported [-0.0953, 0.1916] on 7 chromosomes. A ci95 field is read as a 95% interval
    whatever is written beside it, so below the minimum there is no ci95 at all."""
    pairs = _pairs_on(7)
    a = [0.9 if p.regulated else 0.1 for p in pairs]
    b = [0.5 for _ in pairs]
    for out in (
        crispri.gain_interval(a, b, pairs, n=200),
        crispri.weighted_gain(a, b, pairs, weighted=False, n=200),
    ):
        assert out["clusters"] == 7
        assert out["enough_clusters"] is False
        assert out["ci95"] is None
        assert "interval unreliable: 7 clusters" in out["interval_unreliable"]


def test_the_point_estimate_survives_a_withheld_interval():
    """Only the interval depends on the resampling, so the gain is still reported."""
    pairs = _pairs_on(5)
    a = [0.9 if p.regulated else 0.1 for p in pairs]
    b = [0.5 for _ in pairs]
    out = crispri.gain_interval(a, b, pairs, n=200)
    assert out["gain"] is not None
    assert out["ci95"] is None


def test_enough_clusters_keeps_its_interval():
    pairs = _pairs_on(crispri.MIN_CLUSTERS_FOR_AN_INTERVAL + 2)
    a = [0.9 if p.regulated else 0.1 for p in pairs]
    b = [0.5 for _ in pairs]
    out = crispri.gain_interval(a, b, pairs, n=200)
    assert out["enough_clusters"] is True
    assert out["ci95"] is not None
    assert "interval_unreliable" not in out


def test_the_minimum_is_ten_clusters():
    assert crispri.MIN_CLUSTERS_FOR_AN_INTERVAL == 10
