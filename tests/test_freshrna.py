# SPDX-License-Identifier: AGPL-3.0-or-later
"""M2i (lane-freshrna): the registered form, the fixed scale and the refusal to fit.

The load-bearing test is `test_fit_scale_refuses`: this lane's whole claim rests on nothing being
fitted, so the code must raise when anything asks it for a parameter.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import pytest

from genomeos.attribution import combine2 as c2
from genomeos.attribution import context_evidence as ce
from genomeos.attribution import freshrna as fr
from genomeos.attribution import magcal as mc
from genomeos.attribution import multrna as mr

RESULT = Path("data/results/freshrna_result.json")
SOURCES = Path("data/results/freshrna_sources.json")
REGISTRATION = Path("data/results/freshrna_registration.json")


# ---- nothing is fitted -----------------------------------------------------------------------------
def test_fit_scale_refuses() -> None:
    with pytest.raises(fr.NoFitError):
        fr.fit_scale([{"P": [1.0], "D": [1.0]}])
    with pytest.raises(fr.NoFitError):
        fr.fit_scale()


def test_scale_is_the_constant_of_ac5580b() -> None:
    assert fr.SCALE == 0.2145911381525326
    assert "NOTHING IS FITTED IN THIS LANE" in fr.NOTHING_IS_FITTED
    assert "ac5580b" in fr.SCALE_ORIGIN


def test_the_run_script_never_calls_fit_scale() -> None:
    for p in ("scripts/freshrna_run.py", "scripts/freshrna_register.py"):
        text = Path(p).read_text()
        assert "fit_scale(" not in text, f"{p} calls fit_scale"


# ---- what is carried verbatim ----------------------------------------------------------------------
def test_no_reader_biosample_is_added() -> None:
    assert fr.READER_TERMS is ce.READER_TERMS
    assert len(fr.CELLS) == 13
    assert list(fr.CARRIED["cells"]) == list(c2.CELLS)


def test_endpoint_and_baseline_are_carried_word_for_word() -> None:
    assert fr.MEASURED == c2.MEASURED
    assert c2.PREDICTIONS["no_rule_baseline"] == fr.NO_RULE_BASELINE
    assert fr.PSEUDOCOUNT == c2.PSEUDOCOUNT
    assert fr.MIN_GENES == c2.MIN_GENES
    assert mr.PREDICTIONS["form"] == fr.PREDICTIONS["form"]
    assert fr.FORM_INTERVAL == mc.FORM_INTERVAL
    assert fr.PERMUTATIONS == mc.PERMUTATIONS
    assert c2.BOOTSTRAP["seed"] == fr.BOOTSTRAP["seed"]
    assert c2.BOOTSTRAP["resamples"] == fr.BOOTSTRAP["resamples"]


def test_the_registration_is_complete_and_bound() -> None:
    fields = fr.registered_fields()
    assert set(fields) == set(fr.BOUND)
    for name in fr.BOUND.values():
        assert getattr(fr, name) is not None
    assert "replicates on a fresh measurement of the same cells" in fr.SCOPE.lower()
    assert "2% of the variance" in fr.SIZE.lower()
    for outcome in ("BETTER", "WORSE", "INCONCLUSIVE", "NOT_TESTABLE", "HARNESS_DEFECT"):
        assert outcome in fr.OUTCOMES


# ---- the metric, at the fixed scale ----------------------------------------------------------------
def test_delta_mse_is_combine2s_metric_with_the_scaled_prediction() -> None:
    rows = [{"P": [1.0, -2.0, 0.5], "D": [0.5, -1.0, 0.25]}]
    a = fr.SCALE
    want = sum(2 * a * p * d - (a * p) ** 2 for p, d in zip(rows[0]["P"], rows[0]["D"], strict=True)) / 3
    assert fr.delta_mse(rows) == pytest.approx(want)
    scaled = [{"P": [a * p for p in rows[0]["P"]], "D": rows[0]["D"]}]
    assert fr.delta_mse(rows) == pytest.approx(c2.dmse(scaled))


def test_delta_mse_at_scale_one_is_the_uncalibrated_form() -> None:
    rows = [{"P": [1.0, -2.0], "D": [1.0, -1.0]}]
    assert fr.delta_mse(rows, scale=1.0) == pytest.approx(c2.dmse(rows))


def test_a_perfect_prediction_at_the_right_scale_beats_the_baseline() -> None:
    rows = [{"P": [1.0 / fr.SCALE, -1.0 / fr.SCALE], "D": [1.0, -1.0]}]
    assert fr.delta_mse(rows) == pytest.approx(1.0)


def test_gene_terms_sum_to_the_delta_mse() -> None:
    rows = [{"P": [1.0, -2.0], "D": [0.3, 0.7]}, {"P": [0.5, 0.5], "D": [-0.2, 0.9]}]
    s, n = fr.gene_terms(rows)
    assert sum(s) / sum(n) == pytest.approx(fr.delta_mse(rows))


def test_decide_is_the_three_way_table() -> None:
    assert fr.decide(0.01, 0.05) == "BETTER"
    assert fr.decide(-0.05, -0.01) == "WORSE"
    assert fr.decide(-0.01, 0.05) == "INCONCLUSIVE"
    assert fr.decide(0.0, 0.05) == "INCONCLUSIVE"
    assert "CALIBRATION_REMOVES" not in fr.OUTCOMES


def test_the_split_is_m2hs() -> None:
    for symbol in ("TP53", "OCA2", "A3GALT2", "HERC2"):
        want = int(hashlib.sha256(symbol.encode("utf-8")).hexdigest()[:8], 16) % 2 == 1
        assert fr.is_held_out(symbol) is want
        assert fr.is_held_out(symbol) is mc.is_held_out(symbol)


def test_locus_clusters_group_by_the_megabase() -> None:
    rows = [
        {"chrom": "chr1", "start": 10, "P": [1.0], "D": [1.0]},
        {"chrom": "chr1", "start": 999_999, "P": [1.0], "D": [1.0]},
        {"chrom": "chr1", "start": 1_000_001, "P": [1.0], "D": [1.0]},
        {"chrom": "chr2", "start": 10, "P": [1.0], "D": [1.0]},
    ]
    groups = fr.locus_clusters(rows)
    assert len(groups) == 3
    assert groups[("chr1", 0)] == [0, 1]


def test_the_clustered_bootstrap_brackets_the_point_estimate() -> None:
    rows = [{"chrom": "chr1", "start": i * 100_000, "P": [1.0, -1.0], "D": [0.3, -0.2]} for i in range(60)]
    lo, hi, n = fr.bootstrap_delta_clustered(rows, 200, 1)
    assert n == 6
    assert lo - 1e-9 <= fr.delta_mse(rows) <= hi + 1e-9


def test_the_control_permutes_within_the_gene_and_keeps_its_values() -> None:
    rows = [{"P": [1.0, 2.0, 3.0] + [0.0] * 10, "D": [1.0, 2.0, 3.0] + [0.0] * 10} for _ in range(5)]
    vals = fr.permuted_control(rows, 20, 7)
    assert len(vals) == 20
    assert all(math.isfinite(v) for v in vals)
    # the permutation cannot change sum(P^2), so no permutation can beat the aligned prediction
    assert max(vals) <= fr.delta_mse(rows) + 1e-12


def test_full_reading_carries_the_scope_the_size_and_the_non_overlap() -> None:
    text = fr.full_reading("x: better")
    assert fr.SCOPE in text
    assert fr.SIZE in text
    assert fr.NON_OVERLAP in text


# ---- the committed record --------------------------------------------------------------------------
@pytest.mark.skipif(not SOURCES.exists(), reason="the sources file is not committed yet")
def test_the_fresh_sources_are_disjoint_from_what_the_project_has_read() -> None:
    d = json.loads(SOURCES.read_text())
    ov = d["non_overlap"]
    assert ov["intersection"]["experiments"] == []
    assert ov["intersection"]["files"] == []
    assert ov["intersection"]["libraries"] == []
    assert ov["intersection_with_the_whole_read_set"] == []
    assert set(d["biosamples"]) == set(fr.CELLS)
    assert d["values_read"].startswith("none")
    assert all(a in fr.RNA_ASSAYS for a in d["assays"].values())


@pytest.mark.skipif(not REGISTRATION.exists(), reason="the registration is not committed yet")
def test_the_registration_pins_the_scale_and_the_population() -> None:
    d = json.loads(REGISTRATION.read_text())
    assert d["scale"] == fr.SCALE
    assert d["harness_check"]["passed"] is True
    assert d["harness_check"]["scale_matches"] is True
    assert d["harness_check"]["delta_mse_reproduced"] is True
    assert d["primary_population"]["genes"] == d["harness_check"]["m2h_held_out_genes_committed"]


@pytest.mark.skipif(not RESULT.exists(), reason="the run has not happened yet")
def test_the_result_carries_the_registered_wording() -> None:
    d = json.loads(RESULT.read_text())
    assert d["outcome"] in fr.OUTCOMES
    assert d["scale"] == fr.SCALE
    assert fr.SCOPE in d["reading"]
    assert "NOTHING IS FITTED IN THIS LANE" in d["nothing_is_fitted"]
    assert d["alphagenome_requests"] == 0
    if d["outcome"] in ("BETTER", "WORSE", "INCONCLUSIVE"):
        # SCOPE quotes "the model predicts expression" in order to FORBID it, and `reading` carries
        # SCOPE, so the claim check belongs on the verdict sentence alone
        d = {**d, "reading": d["verdict_sentence"]}
        assert fr.decide(*d["primary"]["ci"]) == d["outcome"]
        assert "the model predicts expression" not in d["reading"]
