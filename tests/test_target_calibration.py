# SPDX-License-Identifier: AGPL-3.0-or-later
"""The CRISPRi calibration: intervals, reliability, per-bin coverage and the transfer to the sweep.

Everything here is synthetic and local: a handful of pairs, a deletion table of three elements, a
registry class file and a GENCODE file written into tmp_path. No network, no cached data.
"""

import gzip
import io
import json
import math

import pytest

from genomeos.attribution import crispri
from genomeos.attribution import target_calibration as tc

HEADER = (
    "chrom\tchromStart\tchromEnd\tname\tEffectSize\tmeasuredGeneSymbol\tValidConnection\tCellType\t"
    "Regulated\tDataset\tdistanceToTSS\tDHS.RPM\tH3K27ac.RPM\tstartTSS\n"
)


def row(chrom, start, end, gene, regulated, distance, cell="K562", dhs=1.0, k27=1.0):
    return (
        f"{chrom}\t{start}\t{end}\t{gene}|x\t-0.1\t{gene}\tTRUE\t{cell}\t"
        f"{'TRUE' if regulated else 'FALSE'}\tD\t{distance}\t{dhs}\t{k27}\t{start + int(distance)}\n"
    )


def element(eid, start, end, gene, k562):
    pred = {"gene": gene, "log2_fold_change": k562, "tissue": "K562"}
    return {
        "id": eid,
        "start": start,
        "end": end,
        "inferred": {"gene": gene, "distance": 1000, "basis": "nearest TSS in domain"},
        "predicted": pred,
        "predicted_by_cell": {"K562": k562},
        "predicted_coding": pred,
        "predicted_coding_by_cell": {"K562": k562},
    }


@pytest.fixture
def world(tmp_path):
    """A deletion table, a registry class file and a GENCODE file that agree with each other."""
    elements = tmp_path / "elements"
    elements.mkdir()
    (elements / "chr1.json").write_text(
        json.dumps(
            [
                element("EH1", 1_000, 1_500, "AAA", -0.4),
                element("EH2", 10_000, 10_500, "BBB", -0.05),
                element("EH3", 20_000, 20_500, "CCC", 0.3),  # a predicted rise: not a drop
                element("EH4", 30_000, 30_500, "NOSUCHGENE", -0.4),  # GENCODE does not know it
            ]
        )
    )
    results = tmp_path / "results"
    results.mkdir()
    with gzip.open(results / "ccres_chr1.bed.gz", "wt") as fh:
        fh.write("# synthetic\n")
        fh.write("chr1\t1000\t1500\tEH1\tdELS\t0\n")
        fh.write("chr1\t10000\t10500\tEH2\tpELS\t1\n")
        fh.write("chr1\t20000\t20500\tEH3\tPLS\t0\n")
        fh.write("chr1\t30000\t30500\tEH4\tCTCF-only\t0\n")
    reference = tmp_path / "reference"
    reference.mkdir()
    with gzip.open(reference / "gencode_v50_chr1.gff3.gz", "wt") as fh:
        fh.write("##gff-version 3\n")
        for name, start in (("AAA", 31_000), ("BBB", 41_000), ("CCC", 51_000)):
            fh.write(f"chr1\tHAVANA\tgene\t{start}\t{start + 900}\t.\t+\t.\tID=g;gene_name={name}\n")
    return {"elements": elements, "results": results, "reference": reference}


def test_wilson_interval_brackets_the_rate():
    lo, hi = tc.wilson(5, 10)
    assert lo < 0.5 < hi
    assert tc.wilson(0, 0) == (0.0, 1.0)
    wide, narrow = tc.wilson(1, 10), tc.wilson(100, 1000)
    assert wide[1] - wide[0] > narrow[1] - narrow[0]  # more pairs, a tighter interval
    assert tc.wilson(0, 100)[0] == 0.0


def test_chi2_survival_matches_known_points():
    assert tc.chi2_sf(3.841, 1) == pytest.approx(0.05, abs=1e-3)
    assert tc.chi2_sf(2.0, 2) == pytest.approx(math.exp(-1.0), abs=1e-9)
    assert tc.chi2_sf(15.507, 8) == pytest.approx(0.05, abs=1e-3)
    assert tc.chi2_sf(0.0, 8) == 1.0


def test_equal_count_bins_split_evenly_and_in_order():
    p = [0.9, 0.1, 0.5, 0.2, 0.8, 0.4]
    groups = tc.equal_count_bins(p, 3)
    assert [len(g) for g in groups] == [2, 2, 2]
    assert [round(min(p[i] for i in g), 1) for g in groups] == [0.1, 0.4, 0.8]
    assert tc.equal_count_bins([], 3) == []


def test_reliability_is_honest_about_a_wrong_probability():
    # 200 pairs at a true rate of 0.5, a model that says 0.5, and one that says 0.05
    labels = [i % 2 == 0 for i in range(200)]
    good = tc.reliability([0.5] * 200, labels, bins=2)
    bad = tc.reliability([0.05] * 200, labels, bins=2)
    assert good["bins_consistent"] == 2 and good["ece"] < 0.02
    assert bad["bins_consistent"] == 0 and bad["ece"] == pytest.approx(0.45, abs=0.01)
    assert bad["hosmer_lemeshow"]["p"] < 0.01 < good["hosmer_lemeshow"]["p"]


def test_class_features_use_dels_as_the_reference_level():
    assert set(tc.class_features("dELS").values()) == {0.0}
    assert tc.class_features("PLS")["class_PLS"] == 1.0
    assert tc.class_features("nonsense")["class_unknown"] == 0.0  # only a literal unknown is the level
    assert tc.class_features("unknown")["class_unknown"] == 1.0


def test_matched_element_prefers_the_element_that_predicts_the_gene():
    els = [element("EHa", 0, 100, "ZZZ", -0.9), element("EHb", 0, 100, "AAA", -0.1)]
    assert tc.matched_element(els, "AAA")["id"] == "EHb"
    assert tc.matched_element(els, "QQQ")["id"] == "EHa"  # no match: the largest magnitude
    assert tc.matched_element([], "AAA") is None


def test_tss_distance_floors_and_reports_the_missing_value():
    assert tc.tss_distance(1_000, 1_200) == float(crispri.MIN_DISTANCE)
    assert tc.tss_distance(1_000, 31_000) == 30_000.0
    assert tc.tss_distance(1_000, None) is None


def pairs_for(world, text):
    pairs = crispri.parse(io.StringIO(text))
    table = crispri.DeletionTable(world["elements"])
    crispri.annotate(pairs, table)
    counts = tc.add_features(pairs, table, world["reference"], world["results"])
    return pairs, table, counts


def test_add_features_counts_every_missing_value_per_arm(world):
    text = (
        HEADER
        + row("chr1", 1_000, 1_500, "AAA", True, 30_000)
        + row("chr1", 10_000, 10_500, "BBB", False, 31_000)
        + row("chr1", 1_000, 1_500, "NOSUCHGENE", True, 5_000)  # no GENCODE entry: counted, not zero
        + row("chr1", 90_000, 90_500, "AAA", False, 1_000)  # no deleted element under it
    )
    pairs, _, counts = pairs_for(world, text)
    by_gene = {p.gene: p for p in pairs if p.start != 90_000}
    assert by_gene["AAA"].features["class_unknown"] == 0.0  # dELS, the reference level
    assert by_gene["AAA"].features["log_tss_distance"] == pytest.approx(math.log(30_999 - 1_250))
    assert by_gene["BBB"].features["class_pELS"] == 1.0
    assert by_gene["NOSUCHGENE"].features["scored"] == 0.0
    assert counts["no_gencode_tss_regulated"] == 1
    assert counts["no_deleted_element_not regulated"] == 1
    assert counts["scored_regulated"] == 1 and counts["pairs_regulated"] == 2
    assert len(tc.scored(pairs)) == 2
    # the GENCODE distance is compared with the benchmark's own column, never silently replaced
    agreement = tc.distance_agreement([p for p in pairs if p.covered])
    assert agreement["pairs"] == 2 and agreement["median_bp"] >= 0


def test_stratum_coverage_and_bin_coverage_show_the_arm_that_is_missing(world):
    text = (
        HEADER
        + row("chr1", 1_000, 1_500, "AAA", True, 30_000)
        # the same stratum (top target, the same drop) in the other arm, and no GENCODE TSS for it
        + row("chr1", 30_000, 30_500, "NOSUCHGENE", False, 30_000)
    )
    pairs, _, _ = pairs_for(world, text)
    cov = tc.stratum_coverage(pairs)
    strata = [tc.stratum(p) for p in pairs]
    assert strata[0] == strata[1] == "top=1|drop=0.2-0.5"
    # the unscored pair is in the other arm of a stratum: the bin must be able to say so
    b = tc.bin_coverage([strata[0]], cov)
    assert b["coverage_regulated"][2] == 1.0
    assert b["coverage_not_regulated"] == [0, 1, 0.0]


def test_by_drop_band_reports_coverage_before_the_rate(world):
    text = (
        HEADER
        + row("chr1", 1_000, 1_500, "AAA", True, 30_000)
        + row("chr1", 10_000, 10_500, "BBB", False, 31_000)
        + row("chr1", 20_000, 20_500, "CCC", False, 31_000)
    )
    pairs, _, _ = pairs_for(world, text)
    bands = tc.by_drop_band([p for p in pairs if p.covered and p.features["top_target"]])
    keys = list(bands[0])
    assert keys.index("coverage_regulated") < keys.index("rate")
    zero, big = bands[0], bands[3]
    assert zero["drop"] == "= 0" and zero["pairs"] == 1  # the predicted rise lands at a drop of zero
    assert big["pairs"] == 1 and big["rate"] == 1.0 and big["ci95"][0] < 1.0
    assert sum(b["pairs"] for b in bands) == 3


def test_element_row_and_sweep_band_every_element(world):
    els = json.loads((world["elements"] / "chr1.json").read_text())
    starts = {"AAA": 31_000, "BBB": 41_000}
    classes = {"EH1": ("dELS", False), "EH2": ("pELS", True)}
    assert len(els) == 4
    got = tc.element_row(els[0], "predicted", starts, classes)
    assert got["cls"] == "dELS" and got["drop"] == pytest.approx(0.4)
    assert got["features"]["top_target"] == 1.0 and got["features"]["node_target"] == 1.0
    assert tc.element_row(els[2], "predicted", starts, classes)["missing"] == (
        "no_gencode_tss_for_the_predicted_gene"
    )
    assert tc.element_row({"start": 0, "end": 1}, "predicted", starts, classes) is None

    weights = {"sweep": [0.0] * (len(tc.SWEEP_FEATURES) + 1), "target": [0.0, 1.0, 0.0, 0.0]}
    out = tc.sweep_chromosome("chr1", weights, world["elements"], world["reference"], world["results"])
    assert out["elements"] == 4
    coding = out["coding gene"]
    assert coding["targets"] == 3  # EH4's gene is not in GENCODE
    assert coding["missing"]["no_gencode_tss_for_the_predicted_gene"] == 1
    assert sum(coding["drop_bands"].values()) == 3
    assert list(coding["drop_bands"]) == [b for b in tc.DROP_LABELS if b in coding["drop_bands"]]
    assert sum(coding["predicted_target_bands"].values()) == 3
    assert sum(sum(v.values()) for v in coding["by_registry_class"].values()) == 3
    assert tc.sweep_chromosome("chr9", weights, world["elements"], world["reference"], world["results"]) == {
        "refused": "no deleted elements cached for this chromosome"
    }


def test_drop_band_labels_line_up_with_the_measured_table():
    assert tc.drop_band_label(0.0) == "= 0"
    assert tc.drop_band_label(0.05) == "0 < drop <= 0.1"
    assert tc.drop_band_label(3.0) == "0.5 < drop <= inf"
    assert tc.DROP_LABELS[0] == "= 0" and len(tc.DROP_LABELS) == len(tc.DROP_BANDS)


def test_band_of_covers_the_unit_interval():
    assert tc.band_of(0.0) == "0-0.02"
    assert tc.band_of(0.03) == "0.02-0.05"
    assert tc.band_of(0.6) == "0.5-0.75"
    assert tc.band_of(0.99) == "0.9-1"
    assert tc.band_of(1.0) == "0.9-1"


def test_score_end_to_end_is_deterministic_and_judges_the_preregistration(world):
    """A synthetic screen where the deletion is the truth: the fit must find it and judge itself."""
    lines = [HEADER]
    for i in range(60):
        chrom = "chr1"
        # every third pair is a real target of EH1 (a strong drop), the rest are not
        if i % 3 == 0:
            lines.append(row(chrom, 1_000, 1_500, "AAA", True, 30_000))
        else:
            lines.append(row(chrom, 10_000, 10_500, "QQQ", False, 31_000))
    text = "".join(lines)
    training = crispri.parse(io.StringIO(text))
    heldout = crispri.parse(io.StringIO(text))
    table = crispri.DeletionTable(world["elements"])
    result = tc.score(training, heldout, table, world["reference"], world["results"])
    assert result["preregistered"] == tc.PREREGISTERED_CALIBRATION
    assert result["verdict"] in ("passed", "failed")
    assert set(result["checks"]) >= {"reliable", "ece_better", "brier_better"}
    held = result["heldout_k562"]["calibrated, sweep features"]
    assert held["pairs"] == 20  # only the pairs whose gene GENCODE knows are scored
    assert held["reliability"]["rows"][0]["coverage_regulated"][1] >= 1
    assert result["populations"]["predicted target"]["training_pairs"] == 20
    assert "K562 pooled (both, after the held-out test)" in result["measured_by_drop_band"]
    assert tc.SCOPE in result["scope"] and result["falsifies_the_transfer"]
    again = tc.score(
        crispri.parse(io.StringIO(text)),
        crispri.parse(io.StringIO(text)),
        crispri.DeletionTable(world["elements"]),
        world["reference"],
        world["results"],
    )
    assert again["heldout_k562"] == result["heldout_k562"]  # no randomness anywhere


def test_no_network_is_touched(world, monkeypatch):
    def refuse(*a, **k):
        raise AssertionError("the calibration must not reach the network")

    monkeypatch.setattr("urllib.request.urlopen", refuse)
    text = HEADER + row("chr1", 1_000, 1_500, "AAA", True, 30_000)
    pairs, _, _ = pairs_for(world, text)
    assert tc.predict([0.0] * (len(tc.SWEEP_FEATURES) + 1), tc.scored(pairs), tc.SWEEP_FEATURES) == [0.5]


def test_prevalence_shift_describes_a_level_error_without_fixing_the_verdict():
    # a curve right in shape, wrong in level: the true rate is 0.2, the model says about 0.1
    labels = [i % 5 == 0 for i in range(200)]
    p = [0.05 if i % 2 else 0.15 for i in range(200)]
    before = tc.calibration(p, labels, bins=2)
    shifted = tc.prevalence_shift(p, labels, bins=2)
    assert shifted["observed_rate"] == pytest.approx(0.2)
    assert shifted["log_odds_shift"] > 0
    assert shifted["ece"] < before["reliability"]["ece"]
    assert shifted["brier"] < before["brier"]
    assert "not a test" in shifted["note"]
    assert tc.prevalence_shift([], [])["pairs"] == 0
