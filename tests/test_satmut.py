"""Measured grammar from saturation mutagenesis: base calls, AUC, sweeps, loci and pooling."""

from genomeos.knowledge.satmut import auc, base_table, by_element, loci, pool, site_coverage, sweep


def _row(element, pos, alt, coef, p, barcodes=50, chrom="chr7", ref="A"):
    return {
        "element": element,
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "barcodes": barcodes,
        "coef": coef,
        "p": p,
    }


def test_base_table_calls_functional_bases_from_substitutions_only():
    rows = [
        _row("E", 100, "C", -0.8, 1e-9),
        _row("E", 100, "G", -0.2, 1e-3),
        _row("E", 101, "T", 0.05, 0.4),
        _row("E", 102, "-", -1.5, 1e-12),  # a deletion: left out
        _row("E", 103, "G", -2.0, 1e-9, barcodes=3),  # too few barcodes to count
    ]
    t = base_table(rows)
    assert set(t) == {100, 101}  # 102 has no substitution, 103 none with enough barcodes
    assert t[100]["functional"] and t[100]["strong"] and t[100]["effect"] == -0.8 and t[100]["measured"] == 2
    assert not t[101]["functional"] and not t[101]["strong"]
    weak = base_table([_row("E", 5, "C", -0.1, 1e-9)])
    assert weak[5]["functional"] and not weak[5]["strong"]  # significant but under the 0.25 log2 effect


def test_auc_and_loci_grouping():
    assert auc([3, 4], [1, 2]) == 1.0 and auc([1, 2], [3, 4]) == 0.0 and auc([1, 2], [1, 2]) == 0.5
    assert auc([], [1]) is None
    rows = [_row("A", p, "C", 0, 1) for p in (10, 20)]  # spans 10..20
    rows += [_row("A2", p, "C", 0, 1) for p in (12, 13, 14, 15, 16)]  # 12..16, inside A's span
    rows += [_row("B", 900, "C", 0, 1, chrom="chr1")]
    exps = by_element(rows)
    groups = loci(exps)
    assert groups == {"A2": ["A2", "A"], "B": ["B"]}  # overlapping spans are one locus; most measured first


def test_sweep_and_site_coverage():
    bases = {p: {"functional": p in (11, 12), "strong": p == 11, "measured": 3} for p in range(1, 21)}
    hits = [("ETS", 10, 13, 0, 0.95), ("SOX", 0, 3, 1, 0.82)]  # 0-based starts on a sequence at offset 0
    cov = {0.80: site_coverage(hits, 0, 0.80), 0.90: site_coverage(hits, 0, 0.90)}
    assert set(cov[0.90]) == {11, 12, 13} and set(cov[0.80]) == {1, 2, 3, 11, 12, 13}
    rows = {r["threshold"]: r for r in sweep(bases, cov)}
    assert rows[0.90]["recall"] == 1.0 and rows[0.90]["precision"] == round(2 / 3, 4)
    assert rows[0.90]["enrichment"] == round((2 / 3) / (2 / 20), 3) and rows[0.90]["table"] == [2, 1, 0, 17]
    assert rows[0.80]["bases_in_sites"] == 0.3
    strong = {r["threshold"]: r for r in sweep(bases, cov, "strong")}
    assert strong[0.90]["table"] == [1, 2, 0, 17]
    from genomeos.knowledge.satmut import motif_score

    ms = motif_score(hits, 0)
    assert ms[11] == 0.95 and ms[1] == 0.82 and 5 not in ms


def test_pool_counts_each_family_site_once_and_tests_families():
    def result(sites, table):
        block = {
            "bases": 10,
            "phylop_auc": 0.7,
            "motif_auc": 0.55,
            "sweep": [{"table": table}] * 4,
            "conserved_counts": {"in_site": [4, 10], "outside": [2, 40]},
        }
        return {"bases_measured": 100, "functional": block, "strong": block, "sites": sites}

    hox = [{"family": "HOX", "start": s, "functional": True, "strong": s < 5} for s in range(10)]
    dup = [{"family": "HOX", "start": 0, "functional": True, "strong": True}]  # a second matrix, same start
    other = [{"family": "SOX", "start": 100 + s, "functional": False, "strong": False} for s in range(20)]
    out = pool([result(hox + dup + other, [5, 5, 5, 85]), result([], [1, 9, 9, 81])])
    assert out["loci"] == 2 and out["phylop_auc_median"] == 0.7 and out["motif_auc_median"] == 0.55
    assert out["sweep"][0]["recall"] == round(6 / 20, 4) and out["definition"] == "functional"
    assert out["conserved"]["in_site"] == {"functional": 8, "bases": 20, "share": 0.4}
    hox_row = next(r for r in out["families"] if r["family"] == "HOX")
    assert hox_row["sites"] == 10 and hox_row["functional_sites"] == 10 and hox_row["q"] <= 0.05
    assert [r["family"] for r in out["families_q05"]] == ["HOX"]
    strong = pool([result(hox + dup + other, [5, 5, 5, 85])], "strong")
    assert next(r for r in strong["families"] if r["family"] == "HOX")["functional_sites"] == 5
