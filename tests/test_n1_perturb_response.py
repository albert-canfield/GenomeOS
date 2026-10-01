# SPDX-License-Identifier: AGPL-3.0-or-later
"""N1's frozen rules on synthetic inputs only (genomeos/attribution/n1_perturb_response.py).

No test reads the Perturb-seq files or any measurement: every row, count and score below is made up here,
so the rules are checked before any outcome exists.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from genomeos.attribution import n1_perturb_response as n1

# --- a synthetic world -------------------------------------------------------------------------------------

N_CONTROL = 585
N_GENES = 200
N_CELLS = 100.0  # every synthetic row has 100 cells, so X = T / 10


def gid(j: int) -> str:
    return f"ENSG{j:011d}"


GENES = [gid(j) for j in range(N_GENES)]


def control_t(row: int, gene: int, hits: int) -> float:
    """A deterministic control T: `hits` of the 585 rows reach 3.5 for each gene, the rest stay small."""
    return 3.5 if (row + gene) % N_CONTROL < hits else 0.1 * (((row * 7 + gene * 3) % 11) - 5) / 5


class FakeReader:
    """Serves synthetic rows and records what was read, in order."""

    def __init__(self, world: dict):
        self.w = world
        self.calls: list[str] = []

    def identities(self):
        self.calls.append("identities")
        ident = {"obs_index": self.w["obs_index"], "var_gene_id": GENES}
        raw = self.w.get("raw_identity", ident)
        return {"normalized": ident, "raw": raw}

    def control_rows(self, rows, cols):
        self.calls.append("control_rows")
        assert all(self.w["obs_index"][i].split("_")[1] == n1.NON_TARGETING for i in rows)
        return {
            "x": [[self.w["x"][i][c] for c in cols] for i in rows],
            "raw_x": [[self.w["raw"][i][c] for c in cols] for i in rows],
            "n": [self.w["n"][i] for i in rows],
        }

    def knockdown_fields(self, rows):
        self.calls.append("knockdown_fields")
        return {i: dict(self.w["kd"][i]) for i in rows}

    def response_rows(self, rows, cols):
        self.calls.append("response_rows")
        return {i: [self.w["x"][i][c] for c in cols] for i in rows}


def make_world(
    n_factors: int = 40,
    low_stratum_hits: int = 2,
    clusters: int = 12,
    fold: float = 0.2,
    attribution_better: bool = True,
) -> tuple[dict, dict]:
    """585 control rows and `n_factors` factor rows over 200 genes, with their frozen plan."""
    obs, x, raw, n, kd = [], {}, {}, {}, {}
    for r in range(N_CONTROL):
        i = len(obs)
        obs.append(f"{20000 + r}_non-targeting_non-targeting_non-targeting")
        x[i] = [
            control_t(r, j, low_stratum_hits if j < 20 else 2) / math.sqrt(N_CELLS) for j in range(N_GENES)
        ]
        raw[i] = [0.01 * (j + 1) for j in range(N_GENES)]  # gene j's control expression rises with j
        n[i] = N_CELLS
    cands = []
    for k in range(n_factors):
        i = len(obs)
        fgene = gid(1000 + k)
        obs.append(f"{k}_F{k}_P1P2_{fgene}")
        responders = {GENES[(k * 5 + q) % N_GENES] for q in range(20)}
        x[i] = [(5.0 if g in responders else 0.0) / math.sqrt(N_CELLS) for g in GENES]
        n[i] = N_CELLS
        kd[i] = {"num_cells_filtered": N_CELLS, "control_expr": 1.0, "fold_expr": fold, "pct_expr": fold - 1}
        hit = sorted(responders)[:15]
        miss = [g for g in GENES if g not in responders][:15]
        good, poor = {g: 0.9 for g in hit + miss[:5]}, {g: 0.9 for g in hit[:5] + miss}
        cands.append(
            {
                "factor": f"F{k}",
                "cluster": f"C{k % clusters}",
                "gene_id": fgene,
                "rows": [{"index": i, "label": obs[i], "tss": "P1P2"}],
                "excluded_genes": [fgene],
                "scores": {
                    "attribution": good if attribution_better else poor,
                    "proximity": poor if attribution_better else good,
                },
            }
        )
    world = {"obs_index": obs, "x": x, "raw": raw, "n": n, "kd": kd}
    plan = {
        "identity_digest": n1.identity_digest(obs, GENES),
        "universe": GENES,
        "control_rows": list(range(N_CONTROL)),
        "candidates": cands,
    }
    return world, plan


# --- the calibration gate ---------------------------------------------------------------------------------


def gate_inputs(low_hits: int) -> tuple[list, list, list]:
    xs = [[control_t(r, j, low_hits if j < 20 else 2) / 10 for j in range(N_GENES)] for r in range(N_CONTROL)]
    raw = [[0.01 * (j + 1) for j in range(N_GENES)] for _ in range(N_CONTROL)]
    return xs, [N_CELLS] * N_CONTROL, raw


def test_gate_passes_when_control_tails_are_near_nominal():
    g = n1.calibration_gate(*gate_inputs(low_hits=2))  # 2 of 585 per gene: 1.27x nominal everywhere
    assert g["passed"]
    assert g["overall"]["ratio"] == pytest.approx((2 / 585) / n1.NOMINAL_RATE)
    assert len(g["strata"]) == 10 and all(s["passed"] for s in g["strata"])
    assert g["strata"][0]["expression_max"] < g["strata"][9]["expression_min"]


def test_gate_fails_when_only_the_least_expressed_decile_is_inflated():
    g = n1.calibration_gate(*gate_inputs(low_hits=4))  # stratum 0: 2.53x; overall stays under 1.5x
    assert g["overall"]["passed"]
    assert not g["strata"][0]["passed"] and all(s["passed"] for s in g["strata"][1:])
    assert not g["passed"]


def test_gate_fails_a_stratum_too_thin_to_test():
    xs, n, raw = gate_inputs(low_hits=2)
    g = n1.calibration_gate(
        [r[:100] for r in xs], n, [r[:100] for r in raw]
    )  # 10 genes a decile: 15.8 expected
    assert g["strata"][0]["expected"] < n1.GATE_MIN_EXPECTED
    assert not g["passed"]


def test_gate_sets_aside_uncalibratable_genes_and_bad_rows():
    xs, n, raw = gate_inputs(low_hits=2)
    xs[3][7] = float("nan")
    n[5] = float("nan")
    g = n1.calibration_gate(xs, n, raw)
    assert g["genes_uncalibrated"] == [7] and g["rows_excluded_bad_cell_count"] == 1


def test_a_failed_gate_ends_the_run_before_any_factor_row_is_read():
    world, plan = make_world(low_stratum_hits=4)
    reader = FakeReader(world)
    out = n1.run(plan, reader)
    assert out["status"] == "gate_failed"
    assert reader.calls == ["identities", "control_rows"]
    assert [r["step"] for r in out["reads"]] == ["identities", "calibration gate"]


# --- predictions: missing is not zero, and the shared universe ---------------------------------------------


def synthetic_predictions() -> n1.Predictions:
    results = {
        "chrA": {
            "genes": {
                "G1": {"requires": [{"factor": "F", "score": 0.95}]},
                "G2": {"requires": [{"factor": "OTHER", "score": 0.9}]},  # scanned; F absent: a zero
                "G5": {"requires": []},  # scanned, nothing enriched: a zero for every factor
            },
            "elements": [
                {"id": "E1", "target": "G1", "requires": [{"factor": "F", "score": 0.88}]},
                {"id": "E2", "target": "G1", "requires": [{"factor": "F", "score": 0.93}]},
                {"id": "E3", "target": "G2", "requires": []},  # scanned element with no factor: a zero
                {"id": "E4", "target": None, "requires": [{"factor": "F", "score": 0.99}]},
                {"id": "E5", "target": "G3", "requires": [{"factor": "F", "score": 0.97}]},
                {"id": "E6", "target": "G5", "requires": []},
            ],
        }
    }
    index = {"chrA": {"G1": ["ID1"], "G2": ["ID2"], "G3": ["ID3"], "G5": ["ID5"]}}
    return n1.predictions(results, index)


def test_missing_is_not_zero():
    pred = synthetic_predictions()
    assert pred.promoter["F"] == {"ID1": 0.95}
    assert pred.element["F"] == {"ID1": 0.93, "ID3": 0.97}  # the gene's best element hit
    assert "ID2" in pred.promoter_scanned and "ID2" in pred.element_attributed  # a zero in both arms
    assert "ID3" not in pred.promoter_scanned  # G3's promoter was never scanned: missing
    assert pred.dropped["element_without_target"] == 1
    u = n1.shared_universe({"ID1", "ID2", "ID3", "ID4", "ID5"}, pred)
    assert u["genes"] == ["ID1", "ID2", "ID5"]  # ID3 is missing, not a zero; ID4 has neither
    sc = n1.factor_scores("F", pred, u["genes"])
    assert sc == {"attribution": {"ID1": 0.93}, "proximity": {"ID1": 0.95}}
    cand = {"factor": "F", "cluster": "C", "excluded_genes": [], "scores": sc}
    res = n1.factor_result(cand, u["genes"], [4.0, 0.0, 0.0])
    assert res["genes"] == 3  # ID2 and ID5 enter as zeros; ID3 does not enter at all


def test_shared_universe_counts_every_exclusion():
    pred = synthetic_predictions()
    u = n1.shared_universe({"ID1", "ID2", "ID3", "ID4", "ID5", "ID9"}, pred)
    assert u["coverage"]["universe"] == 3 and u["coverage"]["measured"] == 6
    ex = u["exclusions"]
    assert ex["measured_without_scanned_promoter"] == 3  # ID3, ID4, ID9
    assert ex["measured_with_element_without_promoter"] == 1  # ID3
    assert ex["measured_with_promoter_without_sampled_element"] == 0


# --- exclusions, aliases, the ledger ---------------------------------------------------------------------


def test_the_perturbed_gene_and_its_crispri_neighbours_are_excluded():
    tss = {
        "F": ("chr1", 1_000_000),
        "A": ("chr1", 1_010_000),
        "B": ("chr1", 1_010_001),
        "C": ("chr2", 1_000_000),
    }
    assert n1.excluded_genes("F", tss, ["A", "B", "C"]) == ["A", "F"]
    scores = {"attribution": {"F": 0.9, "B": 0.9}, "proximity": {"C": 0.9}}
    cand = {"factor": "F", "cluster": "F", "excluded_genes": ["F", "A"], "scores": scores}
    res = n1.factor_result(cand, ["F", "A", "B", "C"], [30.0, 30.0, 4.0, 0.0])
    assert res["genes"] == 2 and res["responders"] == 1  # the knocked-down gene's own drop is not counted
    assert res["auroc_attribution"] == 1.0 and res["auroc_proximity"] == 0.0


def test_alias_rule_and_heterodimers():
    assert n1.monomer_name("Max::Myc") is None and n1.monomer_name("Spi1") == "SPI1"
    sym = {"SPI1": {"E1"}, "DUP": {"E2", "E3"}, "LOST": {"E9"}}
    obs = {"OLDNAME": {"E5"}, "SPI1": {"E7"}}
    assert n1.resolve_factor("SPI1", sym, {"E1"}, obs) == {
        "status": "resolved",
        "gene_id": "E1",
        "via": "gencode_v50_name",
    }
    assert n1.resolve_factor("DUP", sym, {"E2"}, obs)["status"] == "ambiguous_gencode_name"
    assert n1.resolve_factor("LOST", sym, set(), obs)["status"] == "not_perturbed"
    assert n1.resolve_factor("OLDNAME", sym, {"E5"}, obs)["via"] == "perturbation_index_symbol"
    assert n1.resolve_factor("NONE", sym, set(), obs)["status"] == "unresolved_name"
    # SPI1 resolves in GENCODE to E1, so the index's own SPI1 row (E7) is never borrowed
    assert n1.resolve_factor("SPI1", sym, {"E7"}, obs)["status"] == "not_perturbed"


def test_candidate_ledger_steps_and_pending_measurement():
    pred = n1.Predictions()
    universe = [gid(j) for j in range(30)]
    pred.promoter_scanned |= set(universe)
    pred.element_attributed |= set(universe)
    for j in range(12):
        pred.promoter["GOOD"][gid(j)] = 0.9
        pred.element["GOOD"][gid(j + 12)] = 0.9
        pred.promoter["THIN"][gid(j)] = 0.9
    pred.element["THIN"][gid(0)] = 0.9
    rows = [
        {**n1.parse_row(f"{i}_{s}_P1_{e}"), "index": i}
        for i, (s, e) in enumerate([("GOOD", "ENSGGOOD"), ("THIN", "ENSGTHIN"), ("GOOD", "ENSGGOOD")])
    ]
    rows.append({**n1.parse_row("9_non-targeting_non-targeting_non-targeting"), "index": 9})
    sym = {"GOOD": {"ENSGGOOD"}, "THIN": {"ENSGTHIN"}, "NOTHERE": {"ENSGNOTHERE"}}
    ledger = n1.candidate_ledger(
        ["GOOD", "THIN", "NOTHERE", "A::B", "nosuch"], sym, rows, pred, universe, {}, {"GOOD": "Fam"}
    )
    by = {line["jaspar_name"]: line for line in ledger}
    assert by["A::B"]["status"] == "excluded" and "heterodimer" in by["A::B"]["steps"]["monomer"]
    assert by["NOSUCH"]["steps"]["resolved"] == "excluded: unresolved_name"
    assert by["NOTHERE"]["steps"]["perturbed"].startswith("excluded")
    assert by["THIN"]["status"] == "excluded" and by["THIN"]["positives"] == {
        "attribution": 1,
        "proximity": 12,
    }
    good = by["GOOD"]
    assert good["status"] == "candidate" and good["cluster"] == "Fam" and len(good["rows"]) == 2
    assert good["steps"]["knockdown_eligible"] == n1.PENDING and good["steps"]["auroc_defined"] == n1.PENDING


# --- knockdown eligibility, field semantics, the multi-row rule --------------------------------------------


def test_knockdown_rules_in_order():
    base = {"num_cells_filtered": 100.0, "control_expr": 0.5, "fold_expr": 0.2, "pct_expr": -0.8}
    assert n1.knockdown_status(base) == {"eligible": True, "reason": None}
    assert n1.knockdown_status({**base, "pct_expr": -80.0})["eligible"]  # the same reading, in percent
    assert n1.knockdown_status({**base, "num_cells_filtered": 24.0})["reason"] == "too_few_cells"
    low = {**base, "control_expr": 0.09}  # 100 cells x 0.09 = 9 expected target UMIs
    assert n1.knockdown_status(low)["reason"] == "knockdown_unassessable"
    assert n1.knockdown_status({**base, "control_expr": float("nan")})["reason"] == "knockdown_unassessable"
    assert (
        n1.knockdown_status({**base, "fold_expr": 0.41, "pct_expr": -0.59})["reason"]
        == "insufficient_knockdown"
    )
    assert n1.knockdown_status({**base, "pct_expr": 0.2})["reason"] == "semantics_mismatch"


def test_a_semantics_mismatch_stops_the_run_before_responses():
    world, plan = make_world()
    first = plan["candidates"][0]["rows"][0]["index"]
    world["kd"][first]["pct_expr"] = world["kd"][first]["fold_expr"]  # pct_expr is not fold - 1
    reader = FakeReader(world)
    out = n1.run(plan, reader)
    assert out["status"] == "semantics_check_failed" and out["rows_mismatched"] == [first]
    assert "response_rows" not in reader.calls


def test_multi_row_rule_pools_eligible_rows_by_cells():
    t = n1.pooled_t([([0.3, -0.1], 100.0), ([0.1, 0.1], 300.0)])
    assert t[0] == pytest.approx((0.3 * 100 + 0.1 * 300) / 400 * 20)
    assert t[1] == pytest.approx((-0.1 * 100 + 0.1 * 300) / 400 * 20)
    assert n1.pooled_t([([0.3], 100.0)]) == [pytest.approx(3.0)]


def test_multi_row_factor_uses_every_eligible_row_and_no_ineligible_one():
    world, plan = make_world()
    cand = plan["candidates"][0]
    i0 = cand["rows"][0]["index"]
    # a second row for the same factor whose knockdown fails, and which responds far more strongly
    i1 = len(world["obs_index"])
    world["obs_index"].append(f"999_F0_P2_{cand['gene_id']}")
    world["x"][i1] = [9.0 for _ in GENES]
    world["kd"][i1] = {"num_cells_filtered": 300.0, "control_expr": 1.0, "fold_expr": 0.9, "pct_expr": -0.1}
    cand["rows"].append({"index": i1, "label": world["obs_index"][i1], "tss": "P2"})
    plan["identity_digest"] = n1.identity_digest(world["obs_index"], GENES)
    out = n1.run(plan, FakeReader(world))
    f0 = next(f for f in out["factors"] if f["factor"] == "F0")
    assert f0["rows_pooled"] == 1 and f0["responders"] == 20  # the stronger, ineligible row is not used
    row_status = next(e for e in out["eligibility"] if e["factor"] == "F0")["rows"]
    assert [r["reason"] for r in row_status] == [None, "insufficient_knockdown"]
    # now both rows pass: they are pooled by cells, whatever their responses
    world["kd"][i1].update(fold_expr=0.1, pct_expr=-0.9)
    world["x"][i1] = list(world["x"][i0])
    out = n1.run(plan, FakeReader(world))
    f0 = next(f for f in out["factors"] if f["factor"] == "F0")
    assert f0["rows_pooled"] == 2 and f0["responders"] == 20


# --- AUROC, equal weighting, the bootstrap, the criteria, the floor ----------------------------------------


def test_undefined_auroc_is_reported_never_dropped():
    assert n1.auroc([0.1, 0.2], [False, False]) == (None, "no_responders")
    assert n1.auroc([0.1, 0.2], [True, True]) == (None, "all_responders")
    assert n1.auroc([0.0, 0.0, 0.0], [True, False, False]) == (None, "no_score_variation")
    assert n1.auroc([0.0, 0.9, 0.9, 0.0], [False, True, False, False]) == (pytest.approx(5 / 6), None)
    rows = [result(f"D{k}", 0.6, 0.5) for k in range(30)]
    rows += [
        {
            **result("N1", None, None),
            "undefined_attribution": "no_responders",
            "undefined_proximity": "no_responders",
        },
        {**result("N2", None, 0.5), "undefined_attribution": "no_score_variation"},
        {**result("N3", 0.5, None), "undefined_proximity": "no_score_variation"},
    ]
    out = n1.paired_analysis(rows)
    cov = out["coverage"]
    assert cov["factors_analysed"] == 33 and cov["defined"] == 30 and cov["undefined"] == 3
    assert cov["undefined_reasons"] == {
        "no_responders": 1,
        "no_score_variation_attribution": 1,
        "no_score_variation_proximity": 1,
    }


def result(name: str, att: float | None, prox: float | None, cluster: str | None = None) -> dict:
    return {
        "factor": name,
        "cluster": cluster or name,
        "auroc_attribution": att,
        "auroc_proximity": prox,
        "undefined_attribution": None,
        "undefined_proximity": None,
        "difference": att - prox if att is not None and prox is not None else None,
        "genes": 100,
    }


def test_equal_weighting_across_factors():
    rows = [result(f"S{k}", 0.55, 0.5) for k in range(29)]
    rows.append({**result("BIG", 0.9, 0.5), "genes": 10_000, "responders": 5_000})
    out = n1.paired_analysis(rows)
    assert out["estimate"] == pytest.approx((29 * 0.05 + 0.4) / 30)  # one factor, one vote


def test_cluster_bootstrap_resamples_families_not_factors():
    # 20 factors of one family gain 1.0; 11 singleton families gain 0. Resampling factors would keep the
    # lower bound far above 0; resampling families leaves the big family out of about a third of draws.
    rows = [result(f"M{k}", 1.0, 0.0, cluster="BIG") for k in range(20)]
    rows += [result(f"S{k}", 0.5, 0.5) for k in range(11)]
    out = n1.paired_analysis(rows)
    assert out["estimate"] == pytest.approx(20 / 31) and out["clusters"] == 12
    assert out["interval_95"][0] == 0.0
    again = n1.paired_analysis(rows)
    assert again["interval_95"] == out["interval_95"]  # the seed is fixed


def test_too_few_clusters_leave_the_lower_bound_not_assessable():
    rows = [result(f"M{k}", 0.6, 0.5, cluster=f"C{k % 9}") for k in range(30)]
    out = n1.paired_analysis(rows)
    assert out["clusters"] == 9 and out["interval_95"] is None
    assert out["criteria"]["lower_bound"] == {
        "rule": "95% cluster-bootstrap lower bound > 0",
        "value": None,
        "met": None,
        "assessable": False,
    }
    assert out["criteria"]["point_gain"]["met"]


def test_the_two_criteria_are_kept_separate():
    # a large mean carried by one family: the gain is met, the lower bound is not
    rows = [result(f"M{k}", 1.0, 0.0, cluster="BIG") for k in range(20)]
    rows += [result(f"S{k}", 0.5, 0.5) for k in range(11)]
    c = n1.paired_analysis(rows)["criteria"]
    assert c["point_gain"]["met"] and not c["lower_bound"]["met"]
    # a small, consistent gain: the lower bound is above 0, the gain is below +0.02
    rows = [result(f"S{k}", 0.51, 0.5) for k in range(30)]
    c = n1.paired_analysis(rows)["criteria"]
    assert c["lower_bound"]["met"] and not c["point_gain"]["met"]
    assert c["point_gain"]["value"] == pytest.approx(0.01)


def test_insufficient_coverage_stops_without_an_estimate():
    rows = [result(f"S{k}", 0.6, 0.5) for k in range(29)]
    rows.append(
        {
            **result("U", None, None),
            "undefined_attribution": "no_responders",
            "undefined_proximity": "no_responders",
        }
    )
    out = n1.paired_analysis(rows)
    assert out["status"] == "insufficient_coverage" and out["estimate"] is None and out["criteria"] is None


def test_too_few_eligible_factors_stop_before_any_response_row_is_read():
    world, plan = make_world()
    for cand in plan["candidates"][:11]:  # 29 of 40 stay eligible
        world["kd"][cand["rows"][0]["index"]].update(fold_expr=0.8, pct_expr=-0.2)
    reader = FakeReader(world)
    out = n1.run(plan, reader)
    assert out["status"] == "insufficient_coverage"
    assert reader.calls == ["identities", "control_rows", "knockdown_fields"]


def test_identity_mismatch_stops_the_run():
    world, plan = make_world()
    plan["identity_digest"] = "0" * 64
    reader = FakeReader(world)
    assert n1.run(plan, reader)["status"] == "identity_mismatch" and reader.calls == ["identities"]
    world, plan = make_world()
    world["raw_identity"] = {"obs_index": world["obs_index"], "var_gene_id": list(reversed(GENES))}
    assert n1.run(plan, FakeReader(world))["status"] == "identity_mismatch"


def test_a_full_synthetic_run_in_the_registered_order():
    world, plan = make_world()
    reader = FakeReader(world)
    out = n1.run(plan, reader)
    assert reader.calls == ["identities", "control_rows", "knockdown_fields", "response_rows"]
    assert out["status"] == "analysed" and out["gate"]["passed"]
    a = out["analysis"]
    assert a["coverage"]["defined"] == 40 and a["clusters"] == 12
    assert a["estimate"] > n1.MIN_GAIN and a["criteria"]["lower_bound"]["met"]
    world, plan = make_world(attribution_better=False)
    a = n1.run(plan, FakeReader(world))["analysis"]
    assert (
        a["estimate"] < 0
        and not a["criteria"]["point_gain"]["met"]
        and not a["criteria"]["lower_bound"]["met"]
    )


def test_no_rule_is_a_parameter():
    import inspect

    for fn in (
        n1.calibration_gate,
        n1.knockdown_status,
        n1.factor_result,
        n1.paired_analysis,
        n1.cluster_bootstrap,
        n1.run,
    ):
        params = set(inspect.signature(fn).parameters)
        assert not params & {"threshold", "alpha", "seed", "floor", "b", "max_fold", "min_cells"}, fn.__name__


# --- the file reader, on a synthetic pair of files --------------------------------------------------------


def write_pair(tmp_path, world: dict) -> dict:
    """A synthetic pair of files in the downloaded files' layout, with a statistic the reader must skip."""
    h5py = pytest.importorskip("h5py")
    import numpy as np

    rows = len(world["obs_index"])
    paths = {}
    for kind in ("normalized", "raw"):
        p = tmp_path / f"{kind}.h5ad"
        with h5py.File(p, "w") as h:
            h.create_dataset("obs/gene_transcript", data=np.array(world["obs_index"], dtype="S"))
            h.create_dataset("var/gene_id", data=np.array(GENES, dtype="S"))
            src = world["x"] if kind == "normalized" else world["raw"]
            h.create_dataset(
                "X", data=np.array([src[i] if i in src else [0.0] * N_GENES for i in range(rows)], "f4")
            )
            for f in n1.H5adPseudobulk.KNOCKDOWN:
                col = [world["kd"].get(i, {}).get(f, float("nan")) for i in range(rows)]
                if f == "num_cells_filtered":
                    col = [world["n"][i] for i in range(rows)]
                h.create_dataset(f"obs/{f}", data=np.array(col, "f8"))
            h.create_dataset("obs/energy_test_p_value", data=np.zeros(rows))  # present, never read
        paths[kind] = p
    return paths


def test_h5ad_reader_checks_md5_and_reads_only_named_datasets(tmp_path):
    world, plan = make_world()
    paths = write_pair(tmp_path, world)
    expected = {k: {"md5": n1.file_digests(p)[0]} for k, p in paths.items()}
    reader = n1.H5adPseudobulk(paths["normalized"], paths["raw"], expected=expected)
    assert set(reader.sha256) == {"normalized", "raw"}
    out = n1.run(plan, reader)
    assert out["status"] == "analysed" and out["analysis"]["coverage"]["defined"] == 40
    assert not any("energy" in d for d in reader.datasets_read)
    assert reader.datasets_read[:4] == [
        "normalized:obs/gene_transcript",
        "normalized:var/gene_id",
        "raw:obs/gene_transcript",
        "raw:var/gene_id",
    ]
    with pytest.raises(ValueError, match="md5"):
        n1.H5adPseudobulk(paths["normalized"], paths["raw"], expected={**expected, "raw": {"md5": "0" * 32}})


@pytest.mark.skip(
    reason="amendment 1: scripts/n1_run.py now refuses and points to scripts/n1_run_v2.py; "
    "test_the_original_run_script_refuses_and_points_to_v2 checks that"
)
def test_the_run_script_needs_an_authorisation_and_runs_once(tmp_path):
    import importlib.util

    from genomeos.results import save_result

    spec = importlib.util.spec_from_file_location(
        "n1_run", Path(__file__).resolve().parents[1] / "scripts/n1_run.py"
    )
    run_script = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(run_script)
    world, plan = make_world()
    paths = write_pair(tmp_path, world)
    expected = {k: {"md5": n1.file_digests(p)[0]} for k, p in paths.items()}
    results = tmp_path / "results"
    results.mkdir()
    with pytest.warns(UserWarning):
        save_result("n1_registration", {"plan": plan}, results)  # outside the registry: warns, writes
    args = (paths["normalized"], paths["raw"])
    with pytest.raises(SystemExit, match="authorisation"):
        run_script.execute(*args, "  ", results, expected)
    out = run_script.execute(*args, "synthetic test", results, expected)
    assert out["status"] == "analysed" and set(out["files_sha256"]) == {"normalized", "raw"}
    assert (results / "n1_result.json").exists()
    with pytest.raises(SystemExit, match="applied once"):
        run_script.execute(*args, "synthetic test", results, expected)


def load_script(name: str):
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        name, Path(__file__).resolve().parents[1] / f"scripts/{name}.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_original_run_script_refuses_and_points_to_v2(tmp_path):
    # amendment 1 supersedes the run of f9a9130: the original runner refuses before reading anything
    run_script = load_script("n1_run")
    with pytest.raises(SystemExit, match="n1_run_v2.py"):
        run_script.execute(tmp_path / "absent.h5ad", tmp_path / "absent_raw.h5ad", "synthetic test", tmp_path)


# --- amendment 1 -------------------------------------------------------------------------------------------


def brute_auroc(scores, labels):
    pos = [s for s, x in zip(scores, labels, strict=True) if x]
    neg = [s for s, x in zip(scores, labels, strict=True) if not x]
    wins = sum(1.0 if a > b else 0.5 if a == b else 0.0 for a in pos for b in neg)
    return wins / (len(pos) * len(neg))


def test_auroc_v2_counts_ties_half_and_a_constant_score_gives_half():
    assert n1.auroc_v2([0.0, 0.0, 0.0], [True, False, False]) == (0.5, None)
    assert n1.auroc_v2([0.9, 0.9], [True, False]) == (0.5, None)
    assert n1.auroc_v2([0.1, 0.2], [False, False]) == (None, "no_responders")
    assert n1.auroc_v2([0.1, 0.2], [True, True]) == (None, "all_responders")
    for k in range(40):  # deterministic cases with many ties, against pair counting
        scores = [((j * 7 + k) % 4) * 0.3 for j in range(15)]
        labels = [((j * 5 + k) % 3) == 0 for j in range(15)]
        assert n1.auroc_v2(scores, labels)[0] == pytest.approx(brute_auroc(scores, labels))
        old = n1.auroc(scores, labels)[0]
        if old is not None:  # where the original was defined, the two agree
            assert old == pytest.approx(n1.auroc_v2(scores, labels)[0])


def test_stratified_auroc_compares_within_expression_strata():
    # the score tracks expression and so do the responders: unstratified it looks informative, within
    # strata it is constant, so the stratified AUROC is 0.5
    strata = [0] * 10 + [1] * 10
    scores = [0.0] * 10 + [1.0] * 10
    labels = [j == 0 for j in range(10)] + [j < 5 for j in range(10)]
    assert n1.auroc_v2(scores, labels)[0] == pytest.approx(62 / 84)
    assert n1.stratified_auroc(scores, labels, strata) == (0.5, None, 2)
    # strata weighted by their responder x non-responder pairs: 9 pairs at 1.0 and 25 at 0.0
    scores = [1.0] + [0.0] * 9 + [0.0] * 5 + [1.0] * 5
    value, _, used = n1.stratified_auroc(scores, labels, strata)
    assert used == 2 and value == pytest.approx(9 / 34)
    # a stratum with one class contributes nothing; none with both is undefined, never dropped silently
    assert n1.stratified_auroc([0.1, 0.2, 0.3, 0.4], [True, True, False, False], [0, 0, 1, 1]) == (
        None,
        "no_stratum_with_both_classes",
        0,
    )


def test_run_v2_stops_before_the_knockdown_rule():
    world, plan = make_world()
    reader = FakeReader(world)
    out = n1.run_v2(plan, reader)
    assert out["status"] == "eligibility_unsupported" and out["gate"]["passed"]
    assert reader.calls == ["identities", "control_rows"]  # no candidate row is read
    world, plan = make_world(low_stratum_hits=4)
    out = n1.run_v2(plan, FakeReader(world))
    assert out["status"] == "gate_failed"


def interleaved_world():
    """make_world with control expression permuted, so each decile holds genes from every index range."""
    world, plan = make_world()
    for i in plan["control_rows"]:
        world["raw"][i] = [0.01 * (((j * 37) % N_GENES) + 1) for j in range(N_GENES)]
    return world, plan


def test_the_analysis_after_the_stop_once_a_documented_amendment_supports_it(monkeypatch):
    monkeypatch.setattr(n1, "ELIGIBILITY_SUPPORTED", True)
    world, plan = interleaved_world()
    plan["candidates"][0]["scores"]["proximity"] = {}  # a constant predictor: 0.5, not undefined
    reader = FakeReader(world)
    out = n1.run_v2(plan, reader)
    assert reader.calls == ["identities", "control_rows", "knockdown_fields", "response_rows"]
    assert out["status"] == "analysed"
    f0 = next(f for f in out["factors"] if f["factor"] == "F0")
    assert f0["auroc_proximity"] == 0.5 and f0["difference"] is not None and f0["strata_used"] > 1
    a = out["analysis"]
    assert a["coverage"]["defined"] == 40 and a["coverage"]["undefined_reasons"] == {}
    assert a["criteria"]["point_gain"]["met"] and a["criteria"]["lower_bound"]["met"]
    assert a["secondary_unstratified"]["factors"] == 40


def test_freeze_problems_refuse_any_difference(tmp_path):
    files = {}
    for name in ("original.json", "module.py", "runner.py"):
        files[name] = tmp_path / name
        files[name].write_text(name)
    amendment = {
        "amends": {"sha256": n1.sha256_file(files["original.json"])},
        "frozen_code": {
            "module": {"sha256": n1.sha256_file(files["module.py"])},
            "runner": {"sha256": n1.sha256_file(files["runner.py"])},
        },
        "constants": json.loads(json.dumps(n1.CONSTANTS_V2)),
    }
    args = (files["original.json"], files["module.py"], files["runner.py"])
    assert n1.freeze_problems(amendment, *args) == []
    files["module.py"].write_text("edited")
    assert [p.split(":")[0] for p in n1.freeze_problems(amendment, *args)] == ["analysis module"]
    files["module.py"].write_text("module.py")
    amendment["constants"] = {**amendment["constants"], "responds_t": 2.5}
    assert n1.freeze_problems(amendment, *args) == [
        "constants: CONSTANTS_V2 differ from the registered constants"
    ]


def test_the_v2_run_script_checks_the_freeze_first_then_runs_once(tmp_path):
    from genomeos.results import save_result

    run_v2 = load_script("n1_run_v2")
    world, plan = make_world()
    paths = write_pair(tmp_path, world)
    expected = {k: {"md5": n1.file_digests(p)[0]} for k, p in paths.items()}
    args = (paths["normalized"], paths["raw"])

    def results_with(constants: dict) -> Path:
        d = tmp_path / f"results{len(list(tmp_path.glob('results*')))}"
        d.mkdir()
        with pytest.warns(UserWarning):
            save_result("n1_registration", {"plan": plan}, d)
            amendment = {
                "amends": {"sha256": n1.sha256_file(d / "n1_registration.json")},
                "frozen_code": {
                    "module": {"sha256": n1.sha256_file(run_v2.MODULE)},
                    "runner": {"sha256": n1.sha256_file(run_v2.RUNNER)},
                },
                "constants": constants,
            }
            save_result(n1.AMENDMENT, amendment, d)
        return d

    bad = results_with({**json.loads(json.dumps(n1.CONSTANTS_V2)), "kd_max_fold": 0.5})
    wrong_md5 = {k: {"md5": "0" * 32} for k in expected}
    with pytest.raises(SystemExit, match="differ from the freeze"):  # refused before any file is opened
        run_v2.execute(*args, "synthetic test", bad, wrong_md5)
    good = results_with(json.loads(json.dumps(n1.CONSTANTS_V2)))
    with pytest.raises(SystemExit, match="authorisation"):
        run_v2.execute(*args, " ", good, expected)
    out = run_v2.execute(*args, "synthetic test", good, expected)
    assert out["status"] == "eligibility_unsupported" and set(out["files_sha256"]) == {"normalized", "raw"}
    assert not any("obs/fold_expr" in d for d in out["datasets_read"])
    with pytest.raises(SystemExit, match="applied once"):
        run_v2.execute(*args, "synthetic test", good, expected)


def test_reader_v2_decodes_only_the_selected_columns(tmp_path, monkeypatch):
    h5py = pytest.importorskip("h5py")
    world, plan = make_world()
    paths = write_pair(tmp_path, world)
    expected = {k: {"md5": n1.file_digests(p)[0]} for k, p in paths.items()}
    reader = n1.H5adPseudobulkV2(paths["normalized"], paths["raw"], expected=expected)
    seen = []
    original = h5py.Dataset.__getitem__

    def spy(self, args, *a, **k):
        if self.name == "/X":
            seen.append(args)
        return original(self, args, *a, **k)

    monkeypatch.setattr(h5py.Dataset, "__getitem__", spy)
    cols = [150, 3, 77, 12]  # unsorted on purpose: values come back in the order asked
    rows = plan["control_rows"][:5]
    got = reader.control_rows(rows, cols)
    assert seen and all(
        isinstance(a, tuple) and isinstance(a[0], int) and len(a[1]) == len(cols) for a in seen
    )
    assert got["x"][0] == pytest.approx([world["x"][rows[0]][c] for c in cols], rel=1e-6, abs=1e-7)
    assert got["raw_x"][4] == pytest.approx([world["raw"][rows[4]][c] for c in cols], rel=1e-6)
    assert "normalized:X[5 rows x 4 columns]" in reader.datasets_read
