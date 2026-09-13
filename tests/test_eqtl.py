# SPDX-License-Identifier: AGPL-3.0-or-later
"""GTEx eQTL archive read member by member, pairs kept inside elements, targets judged against eGenes."""

import gzip
import io
import tarfile

from genomeos.attribution import eqtl

PAIRS = (
    "variant_id\tgene_id\ttss_distance\tma_samples\tma_count\tmaf\tpval_nominal\tslope\tslope_se\n"
    "chr21_1000_A_G_b38\tENSG00000001.5\t-500\t10\t12\t0.1\t1e-8\t0.5\t0.05\n"
    "chr21_1500_C_T_b38\tENSG00000002.1\t900\t10\t12\t0.1\t1e-6\t-0.3\t0.05\n"
    "chr21_9000_C_T_b38\tENSG00000003.1\t900\t10\t12\t0.1\t1e-6\t-0.3\t0.05\n"
)


def _tar(tmp_path):
    p = tmp_path / "eqtl.tar"
    with tarfile.open(p, "w") as tar:
        for tissue in ("Liver", "Brain_Cerebellum"):
            data = gzip.compress(PAIRS.encode())
            info = tarfile.TarInfo(f"GTEx_Analysis_v8_eQTL/{tissue}.v8.signif_variant_gene_pairs.txt.gz")
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
            small = gzip.compress(b"gene_id\n")
            info = tarfile.TarInfo(f"GTEx_Analysis_v8_eQTL/{tissue}.v8.egenes.txt.gz")
            info.size = len(small)
            tar.addfile(info, io.BytesIO(small))
    return str(p)


def test_tar_members_and_streaming(tmp_path):
    src = _tar(tmp_path)
    members = eqtl.tar_members(src)
    assert [m[0].rsplit("/", 1)[-1] for m in members] == [
        "Liver.v8.signif_variant_gene_pairs.txt.gz",
        "Liver.v8.egenes.txt.gz",
        "Brain_Cerebellum.v8.signif_variant_gene_pairs.txt.gz",
        "Brain_Cerebellum.v8.egenes.txt.gz",
    ]
    name, off, size = members[0]
    rows = list(eqtl.stream_pairs(src, off, size))
    assert rows[0] == ("chr21", 1000, "A", "G", "ENSG00000001", "0.5", "1e-8")
    assert eqtl.tissue_of(name) == "Liver"


def test_intervals_find_covering_elements():
    iv = eqtl.Intervals()
    iv.add("chr21", 900, 1200, "E1")  # covers 1-based 901..1200
    iv.add("chr21", 1000, 1600, "E2")
    iv.add("chr22", 0, 10, "X")
    iv.freeze()
    assert sorted(iv.at("chr21", 1000)) == ["E1"]
    assert sorted(iv.at("chr21", 1100)) == ["E1", "E2"]
    assert iv.at("chr21", 1500) == ["E2"] and iv.at("chr21", 9000) == [] and iv.at("chr1", 5) == []


def test_distil_keeps_only_hits_and_is_resumable(tmp_path):
    src = _tar(tmp_path)
    iv = eqtl.Intervals()
    iv.add("chr21", 900, 1200, "E1")
    iv.add("chr21", 1000, 1600, "E2")
    iv.freeze()
    know = tmp_path / "gtex"
    s = eqtl.distil(iv, know, src, tissues=["Liver"])
    assert s["tissues"] == 1 and s["pairs_scanned_this_run"] == 3 and s["hits_this_run"] == 2
    hits = eqtl.load_hits(know)
    assert set(hits) == {"E1", "E2"}
    assert [h["gene_id"] for h in hits["E1"]] == ["ENSG00000001"]
    assert [h["gene_id"] for h in hits["E2"]] == ["ENSG00000002"]
    s2 = eqtl.distil(iv, know, src)  # Liver is skipped, Brain streamed
    assert s2["tissues"] == 2 and s2["pairs_scanned_this_run"] == 3
    assert {h["tissue"] for h in eqtl.load_hits(know)["E2"]} == {"Liver", "Brain_Cerebellum"}
    assert eqtl.load_hits(know, "chr22") == {}


def test_tissue_matching():
    assert eqtl.tissue_matches("Brain_Cerebellum", "Brain_Cerebellum") is True
    assert eqtl.tissue_matches("cerebellum", "Brain_Cerebellum") is None  # not a GTEx-style track name
    assert eqtl.tissue_matches("Brain_Cortex", "Brain_Cerebellum") is True  # same organ
    assert eqtl.tissue_matches("Whole_Blood", "Liver") is False
    assert eqtl.tissue_matches("K562", "Liver") is None and eqtl.tissue_matches(None, "Liver") is None


def test_score_and_summarise():
    symbols = {"ENSG00000001": "A", "ENSG00000002": "B"}
    hits = [
        {"tissue": "Liver", "gene_id": "ENSG00000001"},
        {"tissue": "Brain_Cortex", "gene_id": "ENSG00000001"},
        {"tissue": "Liver", "gene_id": "ENSG00000002"},
    ]
    e1 = {
        "id": "E1",
        "predicted": {"gene": "A", "tissue": "Liver"},
        "predicted_coding": {"gene": "A"},
        "inferred": {"gene": "C"},
    }
    e2 = {"id": "E2", "predicted": {"gene": "Z", "tissue": "K562"}, "inferred": {"gene": "B"}}
    e3 = {"id": "E3", "predicted": None, "inferred": {"gene": "B"}}
    r1 = eqtl.score_element(e1, hits, symbols)
    assert r1["egenes"] == ["A", "B"] and r1["predicted_in_egenes"] and r1["inferred_in_egenes"] is False
    assert r1["predicted_tissue_is_an_eqtl_tissue"] is True
    r2 = eqtl.score_element(e2, hits[2:], symbols)
    assert r2["predicted_in_egenes"] is False and r2["inferred_in_egenes"] is True
    assert r2["predicted_tissue_is_an_eqtl_tissue"] is None
    r3 = eqtl.score_element(e3, [], symbols)
    assert r3["n_egenes"] == 0 and r3["predicted_in_egenes"] is None
    s = eqtl.summarise([r1, r2, r3])
    assert s["elements"] == 3 and s["with_eqtl"] == 2 and s["fraction_with_eqtl"] == 0.667
    assert s["predicted_target_is_an_egene"] == 0.5 and s["inferred_target_is_an_egene"] == 0.5
    assert s["when_they_disagree"] == {
        "elements": 2,
        "predicted_right": 1,
        "inferred_right": 1,
        "both_right": 0,
        "neither_right": 0,
    }
    assert s["both_judged"] == 2 and s["either_target_is_an_egene"] == 1.0
    assert s["coding_judged"] == 1 and s["when_coding_disagrees"] == {
        "elements": 1,
        "predicted_right": 1,
        "inferred_right": 0,
    }
    assert s["predicted_tissue_is_an_eqtl_tissue"] == 1.0 and s["tissue_judged"] == 1
