# SPDX-License-Identifier: AGPL-3.0-or-later
"""The re-banding: a band is allowed only when the population's own rate can carry it."""

from __future__ import annotations

import gzip
import json

import pytest

from genomeos.attribution import crispri
from genomeos.attribution import target_calibration as tc


def test_bands_touched_covers_an_interval_inside_one_band() -> None:
    assert tc.bands_touched(0.92, 0.99) == ["0.9-1"]
    assert tc.bands_touched(0.26, 0.49) == ["0.25-0.5"]


def test_bands_touched_reports_every_band_a_wide_interval_reaches() -> None:
    assert tc.bands_touched(0.43, 0.72) == ["0.25-0.5", "0.5-0.75"]
    assert tc.bands_touched(0.0, 1.0) == list(tc.BAND_LABELS)


def test_a_tight_interval_on_enough_pairs_keeps_a_band() -> None:
    v = tc.observed_band(46, 46)
    assert v["band"] == "0.9-1"
    assert v["ci95"][0] > 0.9
    assert v["bands_the_interval_touches"] == ["0.9-1"]


def test_an_interval_spanning_two_bands_is_not_calibrated() -> None:
    """The registered outcome for the drop = 0 stratum: 25 of 43 is 0.58 [0.43, 0.72]."""
    v = tc.observed_band(25, 43)
    assert v["band"] == tc.NOT_CALIBRATED
    assert "spans 2 published bands" in v["why"]
    assert v["observed_rate"] == 0.5814


def test_the_count_is_checked_before_the_interval() -> None:
    """46 of 46 has an interval inside 0.9-1, and a population too thin still buys no band."""
    assert tc.observed_band(46, 46)["band"] == "0.9-1"
    thin = tc.observed_band(46, 46, min_pooled=100)
    assert thin["bands_the_interval_touches"] == ["0.9-1"]
    assert thin["band"] == tc.NOT_CALIBRATED
    assert thin["why"] == "fewer than 100 pooled pairs"


def test_a_perfect_run_of_29_still_cannot_name_a_band() -> None:
    """The count rule is not the only guard: at 29 of 29 the interval reaches 0.75-0.9 as well."""
    v = tc.observed_band(29, 29, min_pooled=1)
    assert v["bands_the_interval_touches"] == ["0.75-0.9", "0.9-1"]
    assert v["band"] == tc.NOT_CALIBRATED


def test_an_empty_population_says_so_rather_than_banding_zero() -> None:
    v = tc.observed_band(0, 0)
    assert v["band"] == tc.NOT_CALIBRATED
    assert v["pairs"] == 0
    assert "observed_rate" not in v


def test_move_of_names_the_direction_and_treats_a_lost_band_apart() -> None:
    assert tc.move_of("0.5-0.75", "0.9-1") == "up"
    assert tc.move_of("0.9-1", "0.25-0.5") == "down"
    assert tc.move_of("0.5-0.75", "0.5-0.75") == "unchanged"
    assert tc.move_of("0.9-1", tc.NOT_CALIBRATED) == "lost its band"


def _pair(gene: str, regulated: bool, drop: float) -> crispri.Pair:
    p = crispri.Pair(
        chrom="chr1",
        start=100,
        end=200,
        gene=gene,
        cell=tc.CELL,
        dataset="D",
        distance=1000.0,
        dhs=1.0,
        h3k27ac=1.0,
        regulated=regulated,
    )
    p.covered = True
    p.features = {"deletion_drop": drop, "top_target": 1.0, "scored": 1.0}
    return p


def test_population_bands_splits_and_bands_each_side_from_its_own_pairs() -> None:
    rows = [_pair(f"G{i}", True, 0.6) for i in range(40)]
    rows += [_pair(f"H{i}", i < 20, 0.0) for i in range(40)]
    t = tc.population_bands(rows, tc.drop_stratum)
    assert t["0.5 < drop <= inf"]["band"] == "0.9-1"
    assert t["0.5 < drop <= inf"]["pairs"] == 40
    # half of forty is 0.50 [0.35, 0.65]: two bands, so no number
    assert t["= 0"]["band"] == tc.NOT_CALIBRATED
    assert t["= 0"]["observed_rate"] == 0.5


def test_the_registration_names_the_rule_it_is_judged_by() -> None:
    reg = tc.PREREGISTERED_BANDS
    assert "MIN_POOLED_FOR_A_BAND = 30" in reg
    assert "not calibrated here" in reg
    assert "THE FALSIFIER." in reg
    assert "612,323" in reg and "593,765" in reg


@pytest.fixture
def world(tmp_path):
    """Four elements on chr1, a registry class file and a GENCODE file that agree with them."""
    els = tmp_path / "elements"
    els.mkdir()
    rows = [
        {
            "id": f"E{i}",
            "start": 1000 + i * 10_000,
            "end": 1200 + i * 10_000,
            "predicted": {"gene": "AAA", "tissue": "K562"},
            "predicted_coding": {"gene": "AAA", "tissue": "placenta"},
            "predicted_by_cell": {tc.CELL: -0.6 if i else -0.01},
            "predicted_coding_by_cell": {tc.CELL: -0.6 if i else -0.01},
        }
        for i in range(4)
    ]
    (els / "chr1.json").write_text(json.dumps(rows))
    ref = tmp_path / "reference"
    ref.mkdir()
    with gzip.open(ref / "gencode_v50_chr1.gff3.gz", "wt") as fh:
        fh.write("##gff-version 3\n")
        fh.write("chr1\tHAVANA\tgene\t500000\t500900\t.\t+\t.\tID=g;gene_name=AAA\n")
    return {"elements": els, "reference": ref, "results": tmp_path}


def test_reband_chromosome_keeps_the_published_band_beside_the_new_one(world) -> None:
    """Every target is counted once under each of the two bands, and the moves add up."""
    els, ref, tmp_path = world["elements"], world["reference"], world["results"]
    weights = {"target": [0.0, 1.0, 0.0, 0.0], "sweep": [0.0] * 10}
    stratum_band = {
        "0.5 < drop <= inf": {"band": "0.9-1"},
        "0 < drop <= 0.1": {"band": tc.NOT_CALIBRATED},
    }
    out = tc.reband_chromosome("chr1", weights, stratum_band, elements=els, reference=ref, results=tmp_path)
    any_gene = out["any gene"]
    assert any_gene["targets"] == 4
    assert sum(any_gene["published_bands"].values()) == 4
    assert sum(any_gene["rebanded"].values()) == 4
    assert sum(any_gene["moves"].values()) == 4
    assert any_gene["rebanded"][tc.NOT_CALIBRATED] == 1  # the one element with a drop of 0.01
    assert any_gene["moves"]["lost its band"] == 1
    # the alternative axis is counted per arm, on that arm's own target
    assert any_gene["named_on_track"] == {tc.CELL: 4}
    assert out["coding gene"]["named_on_track"] == {"another track": 4}


def test_reband_chromosome_refuses_a_chromosome_it_has_no_elements_for(tmp_path) -> None:
    out = tc.reband_chromosome("chr9", {"target": [], "sweep": []}, {}, elements=tmp_path)
    assert "refused" in out


def test_a_stratum_the_benchmark_never_reaches_is_banded_rather_than_dropped(world) -> None:
    """A stratum absent from `stratum_band` must not silently vanish: it is `not calibrated here`."""
    els, ref, tmp_path = world["elements"], world["reference"], world["results"]
    weights = {"target": [0.0, 1.0, 0.0, 0.0], "sweep": [0.0] * 10}
    out = tc.reband_chromosome("chr1", weights, {}, elements=els, reference=ref, results=tmp_path)
    any_gene = out["any gene"]
    assert any_gene["targets"] == 4
    assert any_gene["rebanded"] == {tc.NOT_CALIBRATED: 4}
    assert any_gene["moves"] == {"lost its band": 4}
