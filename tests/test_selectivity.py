"""Marker logic: gates, healthy burden, the score and the pre-registered benchmark.

Synthetic atlas and synthetic cell lines throughout, so the scoring rules are
pinned before any real run. No test touches the network.
"""

from __future__ import annotations

from genomeos.therapeutics import atlas, selectivity
from genomeos.therapeutics.providers import hpa_column


def test_essential_columns_are_the_eight_weighed_organs_in_atlas_spelling():
    columns = selectivity.essential_columns()
    assert set(columns) == {
        "heart muscle",
        "cerebral cortex",
        "cerebellum",
        "hypothalamus",
        "bone marrow",
        "lung",
        "liver",
        "kidney",
    }
    assert columns["heart muscle"] == "Tissue RNA - heart muscle [nTPM]"
    assert columns["bone marrow"] == hpa_column("t_RNA_bone_marrow")


def _gates():
    # four lines of one lineage, four of another; TUM is on every tumour line,
    # SHARED too, and HEALTHYMARK only on the second lineage's lines.
    expression = {
        "ACH-1": {"TUM": 6.0, "SHARED": 6.0, "HEALTHYMARK": 0.1},
        "ACH-2": {"TUM": 5.0, "SHARED": 5.0, "HEALTHYMARK": 0.0},
        "ACH-3": {"TUM": 4.0, "SHARED": 0.5, "HEALTHYMARK": 0.0},
        "ACH-4": {"TUM": 0.2, "SHARED": 6.0, "HEALTHYMARK": 0.0},
        "ACH-5": {"TUM": 6.0, "SHARED": 6.0, "HEALTHYMARK": 6.0},
        "ACH-6": {"TUM": 6.0, "SHARED": 6.0, "HEALTHYMARK": 6.0},
        "ACH-7": {"TUM": 6.0, "SHARED": 0.0, "HEALTHYMARK": 6.0},
        "ACH-8": {"TUM": 0.0, "SHARED": 0.0, "HEALTHYMARK": 6.0},
    }
    lineage_of = {m: ("breast" if i < 4 else "lung") for i, m in enumerate(sorted(expression))}
    return selectivity.Gates(expression, lineage_of, min_lines=4)


def test_gates_are_bitmaps_and_coverage_is_per_lineage():
    gates = _gates()
    assert len(gates.lines) == 8
    assert set(gates.lineages) == {"breast", "lung"}
    assert gates.coverage(gates.mask(("TUM",)), "breast") == 0.75
    assert gates.coverage(gates.mask(("TUM", "SHARED")), "breast") == 0.5
    # the NOT marker removes exactly the lines that carry it
    assert gates.coverage(gates.mask(("TUM",), ("HEALTHYMARK",)), "lung") == 0.0
    assert gates.coverage(gates.mask(("TUM",), ("HEALTHYMARK",)), "breast") == 0.75
    assert gates.mask(("ABSENT",)) == 0


def test_a_lineage_with_too_few_lines_is_not_scored():
    expression = {f"ACH-{i}": {"TUM": 6.0} for i in range(3)}
    gates = selectivity.Gates(expression, dict.fromkeys(expression, "tiny"), min_lines=15)
    assert gates.lineages == {}


def test_healthy_burden_is_bounded_by_the_lower_marker_and_vetoed_by_a_not_marker():
    levels = {
        "A": {"liver": 40.0, "lung": 5.0},
        "B": {"liver": 2.0, "lung": 4.0},
        "C": {"liver": 60.0, "lung": 0.0},
    }
    assert selectivity.burden(("A",), (), levels) == (40.0, "liver")
    # an AND gate can only hit a cell carrying both, so the liver falls from 40 towards B's
    # 2 nTPM - but no further than a quarter of 40, because bulk tissue RNA cannot say the two
    # markers are on different cells
    assert selectivity.burden(("A", "B"), (), levels) == (10.0, "liver")
    assert selectivity.burden(("A", "B"), (), levels, co_expression_floor=0.0) == (4.0, "lung")
    # C vetoes the liver, where it is carried above the protection threshold
    assert selectivity.burden(("A",), ("C",), levels) == (5.0, "lung")
    # a gene with no atlas values contributes no tissue rather than a zero
    assert selectivity.burden(("A", "UNKNOWN"), (), levels) == (0.0, None)


def test_the_score_trades_coverage_against_essential_tissue():
    assert selectivity.score(1.0, 0.0) == 1.0
    assert selectivity.score(1.0, 10.0) == 0.5
    assert selectivity.score(0.5, 0.0) > selectivity.score(1.0, 20.0)


def test_the_gate_space_is_singles_and_pairs_and_ordered_and_not_pairs():
    space = selectivity.combinations(["A", "B", "C"])
    kinds = [k for _r, _e, k in space]
    assert kinds.count("A") == 3
    assert kinds.count("A and B") == 3
    assert kinds.count("A and not C") == 6
    assert (("A",), ("B",), "A and not C") in space
    assert (("B",), ("A",), "A and not C") in space


def test_a_gate_that_removes_a_healthy_burden_outranks_the_bare_marker():
    gates = _gates()
    levels = {
        "TUM": {"liver": 30.0, "lung": 1.0},
        "SHARED": {"liver": 0.5, "lung": 0.5},
        "HEALTHYMARK": {"liver": 80.0, "lung": 0.0},
    }
    rows = selectivity.rank_combinations(gates, levels, selectivity.combinations(sorted(levels)))
    bare = next(r for r in rows if r["required"] == ["TUM"] and not r["excluded"])
    gated = next(r for r in rows if r["required"] == ["TUM"] and r["excluded"] == ["HEALTHYMARK"])
    paired = next(r for r in rows if r["required"] == ["SHARED", "TUM"])
    # the bare marker carries the liver; both gates get rid of it, one by vetoing the tissue
    # and one by requiring a second marker the liver does not carry
    assert bare["healthy_burden_ntpm"] == 30.0 and bare["worst_essential_tissue"] == "liver"
    assert gated["healthy_burden_ntpm"] == 1.0 and gated["coverage"] == bare["coverage"]
    assert paired["healthy_burden_ntpm"] == 7.5  # a quarter of TUM's 30, not SHARED's 0.5
    assert gated["score"] > bare["score"] and paired["score"] > bare["score"]
    # the veto is credited in full because a measured presence is evidence; the AND gate is not
    assert gated["score"] > paired["score"]


def test_the_benchmark_claim_is_declared_with_its_thresholds_and_its_drugs():
    assert "0.05" in selectivity.PREREGISTERED_CLAIM
    assert "0.1" in selectivity.PREREGISTERED_CLAIM and "95th percentile" in selectivity.PREREGISTERED_CLAIM
    assert selectivity.MIN_BENCHMARK_DELTA == 0.10
    for gene in ("ERBB2", "CD19", "MSLN", "CEACAM5", "DLL3", "EGFR"):
        assert gene in selectivity.BENCHMARK
    assert "trastuzumab" in selectivity.BENCHMARK["ERBB2"]


def test_the_pool_carries_the_benchmark_and_an_equally_arbitrary_placebo_set():
    genes = [f"G{i:03d}" for i in range(60)] + ["ERBB2", "CD19"]
    expression = {f"ACH-{i}": dict.fromkeys(genes, 6.0) for i in range(20)}
    gates = selectivity.Gates(expression, dict.fromkeys(expression, "breast"), min_lines=4)
    levels = {g: {"liver": 1.0} for g in genes}
    singles = selectivity.rank_combinations(gates, levels, selectivity.single_gates(sorted(genes)))
    pool, placebo = selectivity.build_pool(singles, gates, levels, pool_size=10, placebo_size=5)
    assert len(placebo) == 5
    assert set(placebo) <= set(pool)
    assert not set(placebo) & set(selectivity.BENCHMARK)
    # the benchmark genes are carried in whether or not they scored their way into the top ten
    assert {"ERBB2", "CD19"} <= set(pool)
    assert len(pool) == len({*pool})


def test_the_placebo_null_says_what_a_meaningless_gene_set_scores():
    # scores depend only on the lineage, so no gene set can have an advantage
    rows = [
        {
            "required": [f"G{i}", f"G{(i + 1) % 20}"],
            "excluded": [],
            "logic": "A and B",
            "lineage": "breast",
            "score": round(0.9 - (i % 10) / 100, 3),
        }
        for i in range(60)
    ]
    pool = sorted({g for r in rows for g in r["required"]})
    null = selectivity.placebo_null(rows, pool, size=3, draws=40, sample=20)
    assert null["gene_sets_drawn"] == 40
    assert null["flat"], null
    assert abs(null["mean_delta"]) < 0.2
    # and a planted set is above the null's 95th percentile
    for r in rows:
        if "G1" in r["required"]:
            r["score"] = 1.0
    rows.sort(key=lambda r: -r["score"])
    planted = selectivity.benchmark_test(
        rows, benchmark={"G1": "planted"}, sample=20, pool=pool, null_draws=40
    )
    assert planted["cliffs_delta"] > planted["placebo_null"]["p95_delta"]
    assert planted["delta_above_null_p95"]


def test_the_benchmark_test_detects_both_a_real_and_an_absent_advantage():
    good = [
        {"required": ["ERBB2"], "excluded": [], "logic": "A", "lineage": "breast", "score": 0.9},
        {"required": ["CD19"], "excluded": [], "logic": "A", "lineage": "lymphoid", "score": 0.8},
    ]
    filler = [
        {"required": [f"G{i}"], "excluded": [], "logic": "A", "lineage": "breast", "score": 0.5 - i / 100}
        for i in range(40)
    ]
    supported = selectivity.benchmark_test(good + filler, sample=40)
    assert supported["supported"]
    assert supported["median_percentile_benchmark"] > supported["median_percentile_random"]
    assert supported["best_per_benchmark_gene"]["ERBB2"]["rank"] == 1
    assert supported["best_per_benchmark_gene"]["ERBB2"]["percentile"] > 90
    # a benchmark gene with no gate at all is reported as absent from the pool, not as a zero
    assert "MSLN" in supported["benchmark_genes_not_in_the_pool"]
    assert "MSLN" not in supported["best_per_benchmark_gene"]
    buried = filler[:20] + good + filler[20:]
    for i, row in enumerate(buried):
        row["score"] = 1.0 - i / 100
    not_supported = selectivity.benchmark_test(buried, sample=40)
    assert not not_supported["supported"]


def test_ranked_gates_are_predictions_capped_below_measurement():
    ev = selectivity.evidence_for(
        {
            "required": ["ERBB2"],
            "excluded": ["HEALTHYMARK"],
            "lineage": "breast",
            "coverage": 0.8,
            "healthy_burden_ntpm": 2.0,
            "worst_essential_tissue": "liver",
            "score": 0.9,
        }
    )
    assert ev.source_type == "prediction" and ev.level == "computational"
    assert ev.confidence <= 0.5
    assert "ERBB2 and not HEALTHYMARK" in ev.claim


def test_run_end_to_end_against_a_synthetic_atlas_and_a_synthetic_screen(monkeypatch):
    columns = selectivity.essential_columns()
    table = {
        "genes": {
            "TUM": {"Gene": "TUM", columns["liver"]: 0.2, columns["lung"]: 0.1},
            "SHARED": {"Gene": "SHARED", columns["liver"]: 0.3, columns["lung"]: 0.2},
            "HEALTHYMARK": {"Gene": "HEALTHYMARK", columns["liver"]: 90.0, columns["lung"]: 0.0},
            "ERBB2": {"Gene": "ERBB2", columns["liver"]: 1.0, columns["lung"]: 1.0},
        },
        "meta": {},
    }
    monkeypatch.setattr(atlas, "load", lambda path=None: table)
    expression = {
        f"ACH-{i}": {"TUM": 6.0, "SHARED": 6.0 if i % 2 else 0.0, "HEALTHYMARK": 0.0, "ERBB2": 6.0}
        for i in range(16)
    }
    monkeypatch.setattr(
        selectivity.depmap,
        "matrix",
        lambda *a, **k: {
            "values": expression,
            "genes": ["TUM", "SHARED", "HEALTHYMARK", "ERBB2"],
            "missing": [],
        },
    )
    monkeypatch.setattr(
        selectivity.depmap,
        "models",
        lambda *a, **k: {m: {"OncotreeLineage": "breast"} for m in expression},
    )
    result = selectivity.run(net=False, pool_size=4)
    assert result["data"]["cell_lines"] == 16
    assert result["combinations_scored"] > 0
    assert result["top_combinations"][0]["coverage"] == 1.0
    assert result["benchmark"]["combinations_with_a_benchmark_target"] > 0
    assert "no binder" in result["not_a_design"].lower()
    assert any("cancer cell lines" in limit for limit in result["limits"])
