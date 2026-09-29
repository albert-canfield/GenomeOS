# SPDX-License-Identifier: AGPL-3.0-or-later
"""Audit C checkpoint 1 (lane-entex, 2026-09-29): the EN-TEx feasibility census.
Synthetic rows only: no archive, no network, no allelic table."""

from __future__ import annotations

import ast
import gzip
import hashlib
import importlib.util
import io
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("entex", ROOT / "scripts" / "entex_feasibility.py")
ex = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ex)


def line(hap1=7, hap2=1000, ratio="0.007", p="1e-9", flag="1", assay="ATAC-seq", donor="ENC-001"):
    fields = [
        "chr1",
        "100",
        "400",
        "EH38D1",
        str(hap1),
        str(hap2),
        "ENCSR000AAA",
        donor,
        "thyroid_gland",
        assay,
    ]
    return "\t".join([*fields, ratio, p, flag]) + "\n"


# ---------------------------------------------------------------- the outcome discipline


def test_the_discipline_names_every_column_once():
    named = (
        ex.COLUMNS_READ
        + ex.COLUMNS_SUMMED_NEVER_KEPT
        + ex.COLUMNS_REDUCED_UNSIGNED
        + ex.COLUMNS_NEVER_INDEXED
    )
    assert sorted(named) == sorted(ex.c5.ENTEX_HEADER)
    assert len(named) == len(set(named))


def test_depth_only_keeps_the_sum_and_nothing_signed():
    r = ex.depth_only(line())
    assert r == ("chr1", 100, 400, "EH38D1", "ENCSR000AAA", "ENC-001", "thyroid_gland", "ATAC-seq", 1007, 1)
    assert 7 not in r and 1000 not in r


def test_depth_only_cannot_tell_the_haplotypes_apart():
    assert ex.depth_only(line(7, 1000)) == ex.depth_only(line(1000, 7))
    assert ex.depth_only(line(ratio="0.9", p="0.5")) == ex.depth_only(line(ratio="0.1", p="1e-30"))


def test_depth_only_refuses_an_unknown_flag_without_showing_it():
    with pytest.raises(ValueError, match="value not shown") as e:
        ex.depth_only(line(flag="hap1"))
    assert "hap1" not in str(e.value)
    with pytest.raises(ValueError, match="value not shown"):
        ex.depth_only(line(hap1="x"))
    with pytest.raises(ValueError, match="fields"):
        ex.depth_only("chr1\t1\t2\n")


def test_no_code_indexes_the_ratio_or_p_value_fields():
    tree = ast.parse((ROOT / "scripts" / "entex_feasibility.py").read_text())
    hits = [
        n.lineno
        for n in ast.walk(tree)
        if isinstance(n, ast.Subscript) and isinstance(n.slice, ast.Constant) and n.slice.value in (10, 11)
    ]
    assert hits == []


class FakeSource:
    def __init__(self, text: str) -> None:
        self.url = "http://example.invalid/t.tsv"
        self.bytes = len(text.encode())
        self.sha = hashlib.sha256(text.encode()).hexdigest()

    def record(self, **extra):
        return {"url": self.url, "sha256": self.sha, "bytes": self.bytes, **extra}


def build(tmp_path, monkeypatch, swap: bool) -> bytes:
    rows = {
        "cCREs_default_AS.tsv": [line(7, 1000), line(3, 30, assay="HM-ChIP-seq_CTCF")],
        "genes_default_AS.tsv": [line(11, 500, assay="RNA-seq", flag="0")],
    }

    def fake(name):
        text = "".join(rows[name])
        if swap:  # the same table with the two haplotype columns exchanged
            split = (r.split("\t") for r in text.splitlines(True))
            text = "".join("\t".join(f[:4] + [f[5], f[4]] + f[6:]) for f in split)
        return FakeSource(text), io.StringIO(text)

    d = tmp_path / ("swap" if swap else "plain")
    monkeypatch.setattr(ex, "CACHE", d)
    monkeypatch.setattr(ex, "DEPTH_CACHE", d / "depth_only.tsv.gz")
    monkeypatch.setattr(ex, "DEPTH_SIDECAR", d / "depth_only.json")
    monkeypatch.setattr(ex, "LEDGER", d / "ledger.jsonl")
    monkeypatch.setattr(ex.c5, "_entex_rows", fake)
    side = ex.build_depth_cache(ex.Cost("test"))
    assert side["tables"]["ccre"]["rows"] == 2 and side["tables"]["ccre"]["rows_kept"] == 1
    assert side["experiments"]["ccre"]["HM-ChIP-seq_CTCF"]["ENC-001"]["thyroid_gland"] == ["ENCSR000AAA"]
    return (d / "depth_only.tsv.gz").read_bytes()


def test_the_cache_holds_totals_only_and_is_blind_to_the_haplotype_labels(tmp_path, monkeypatch):
    plain = build(tmp_path, monkeypatch, swap=False)
    assert plain == build(tmp_path, monkeypatch, swap=True)
    rows = [r.split("\t") for r in gzip.decompress(plain).decode().splitlines()]
    assert rows[0] == ex.CACHE_COLUMNS
    assert [r[-2:] for r in rows[1:]] == [["1007", "1"], ["511", "0"]]
    assert not any(v in r for r in rows for v in ("7", "1000", "11", "500", "0.007", "1e-9"))


# ---------------------------------------------------------------- phase, intervals, loci


BLOCKS = [
    (0, 1_000, "NoInfo", "NoInfo"),
    (900, 5_000, "Paternal", "Maternal"),
    (6_000, 9_000, "Maternal", "Paternal"),
]


def test_block_of_names_one_block_or_says_why_not():
    assert ex.block_of(BLOCKS, 100, 200) == (0, "one")
    assert ex.block_of(BLOCKS, 950, 960) == (None, "ambiguous")
    assert ex.block_of(BLOCKS, 4_900, 6_100) == (None, "none")
    assert ex.block_of([], 1, 2) == (None, "none")


def test_orientation_needs_one_block_or_two_known_parents():
    assert ex.orientation(BLOCKS, 1, 1) == "same_block"
    assert ex.orientation(BLOCKS, 1, 2) == "parental_origin"
    assert ex.orientation(BLOCKS, 0, 1) == "not_orientable"
    assert ex.orientation(BLOCKS, None, 1) == "not_orientable"


def test_intervals_find_an_early_long_interval():
    iv = ex.Intervals({"chr1": [(0, 10_000), (100, 200), (300, 400)]})
    assert iv.hits("chr1", 5_000, 5_001)
    assert iv.hits("chr1", 150, 160)
    assert not iv.hits("chr1", 10_000, 10_010)
    assert not iv.hits("chr2", 0, 10)
    assert ex.Intervals({"chr1": [(100, 200)]}).any_of("chr1", [(0, 50), (199, 300)])


def test_loci_are_frozen_on_the_parent_universe_and_a_filter_cannot_split_one():
    tss = {"chr1": {0, 900_000, 1_800_000, 3_000_000}, "chr2": {5}}
    clusters = ex.freeze_clusters(tss)
    assert len(set(clusters.values())) == 3  # 0-1.8 Mb chain, 3 Mb, chr2
    genes = {"chr1": {f"G{t}": {"tss": t} for t in tss["chr1"]}, "chr2": {"G5": {"tss": 5}}}

    def rows(*names):
        return [{"chrom": "chr1", "gene": n, "element": 0, "donor": "ENC-001", "tissue": "t"} for n in names]

    whole = ex.summarize(rows("G0", "G900000", "G1800000"), genes, clusters)
    split = ex.summarize(rows("G0", "G1800000"), genes, clusters)  # recomputing here would give two loci
    assert whole["parent_clusters_touched"] == split["parent_clusters_touched"] == 1
    assert split["parent_clusters_touched_by_donor"] == {"ENC-001": 1}
    assert split["target_gene_loci"] == 2


def test_read_sv_keeps_copy_number_calls_with_a_non_reference_genotype(tmp_path):
    vcf = tmp_path / "sv.vcf"
    vcf.write_text(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts\n"
        "chr1\t1001\ta\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=2000\tGT\t1|0\n"
        "chr1\t5001\tb\tN\t<INS>\t.\tPASS\tSVTYPE=INS;END=5001\tGT\t1|1\n"
        "chr1\t8001\tc\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;END=9000\tGT\t0|0\n"
    )
    iv, types = ex.read_sv(vcf)
    assert iv.n == 1 and iv.hits("chr1", 1500, 1501) and not iv.hits("chr1", 5000, 5002)
    assert types == {"DEL": 1, "INS": 1, "DUP": 1}


# ---------------------------------------------------------------- matching and the rule


def cands(**over):
    base = {
        "gene": np.array(["A", "B", "C", "T"], dtype=object),
        "tss": np.array([12_000, 30_000, 8_000, 10_000]),
        "depth": np.array([100.0, 100.0, 100.0, 100.0]),
        "expr": np.array([3.0, 3.0, 3.0, 3.0]),
        "block": np.array([1, 1, 1, 1]),
        "excluded": np.array([False, False, False, False]),
    }
    return {**base, **over}


ROW = {"gene": "T", "mid": 0, "distance": 10_000, "rna": 100.0, "expr": 3.0}


def test_alternatives_match_on_distance_expression_depth_and_block():
    assert ex.count_alternatives(ROW, cands(), 47, 1) == 2  # A and C; B is 3x as far; T is the target
    assert ex.count_alternatives(ROW, cands(expr=np.array([5.0, 3.0, 3.0, 3.0])), 47, 1) == 1
    assert ex.count_alternatives(ROW, cands(depth=np.array([300.0, 100.0, 20.0, 100.0])), 47, 1) == 0
    assert ex.count_alternatives(ROW, cands(block=np.array([2, 1, 1, 1])), 47, 1) == 1
    assert ex.count_alternatives(ROW, cands(excluded=np.array([True, False, False, False])), 47, 1) == 1
    assert ex.count_alternatives(ROW, cands(), 47, None) == 0
    assert ex.count_alternatives({**ROW, "expr": None}, cands(), 47, 1) == 0
    assert ex.count_alternatives(ROW, None, 47, 1) == 0


def test_balanced_marginals_are_per_side_only():
    rows = [{"atac_flag": 1, "rna_flag": 0}, {"atac_flag": 0, "rna_flag": 0}, {"atac_flag": 0, "rna_flag": 1}]
    assert ex.balanced_marginals(rows) == {"rows": 3, "element_side_balanced": 2, "gene_side_balanced": 2}


def test_the_go_rule_was_fixed_from_the_read_formula():
    assert (ex.GO_FLOOR, ex.GO_LOCI, ex.GO_DONORS, ex.GO_DONOR_LOCI) == (47, 194, 3, 30)
    assert ex.GO_FLOOR in ex.FLOORS


def test_go_decision_needs_loci_donors_and_no_exposure():
    stage = {
        "parent_clusters_touched": 250,
        "parent_clusters_touched_by_donor": {"A": 40, "B": 31, "C": 30, "D": 2},
    }
    assert ex.go_decision(stage, False)["go"]
    assert not ex.go_decision(stage, True)["go"]
    assert not ex.go_decision({**stage, "parent_clusters_touched": 193}, False)["go"]
    few = {**stage, "parent_clusters_touched_by_donor": {"A": 400, "B": 29, "C": 30}}
    assert not ex.go_decision(few, False)["go"]


GTEX_V8_COLUMNS = {
    "Adrenal Gland", "Pancreas", "Breast - Mammary Tissue", "Artery - Coronary", "Esophagus - Muscularis",
    "Esophagus - Mucosa", "Muscle - Skeletal", "Esophagus - Gastroesophageal Junction",
    "Heart - Left Ventricle", "Skin - Sun Exposed (Lower leg)", "Adipose - Visceral (Omentum)", "Ovary",
    "Small Intestine - Terminal Ileum",
    "Prostate", "Heart - Atrial Appendage", "Liver", "Colon - Sigmoid", "Spleen", "Stomach",
    "Adipose - Subcutaneous", "Skin - Not Sun Exposed (Suprapubic)", "Testis", "Artery - Aorta", "Thyroid",
    "Artery - Tibial", "Nerve - Tibial", "Colon - Transverse", "Lung", "Uterus", "Vagina",
}  # fmt: skip


def test_every_tissue_maps_to_a_gtex_column():
    assert set(ex.TISSUE_TO_GTEX.values()) <= GTEX_V8_COLUMNS
    assert ex.norm_tissue("Peyer's patch") == "peyer's_patch"
    assert ex.norm_gtex("Skin - Sun Exposed (Lower leg)") == ex.norm_gtex("Skin_Sun_Exposed_Lower_leg")
