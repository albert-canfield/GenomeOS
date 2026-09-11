# SPDX-License-Identifier: AGPL-3.0-or-later
"""Enhancer to gene, predicted: the logic runs on an injected scorer, no network."""

from genomeos.coords import Locus
from genomeos.predict.enhancer_target import aggregate, compare, predict_target, score_element, summarise


def fake_scorer(chrom, pos, ref, alt):
    # three genes, two tracks: GENE_A drops hard in liver, GENE_B rises a little, GENE_C is flat
    return [
        ("GENE_A", "liver", -0.45),
        ("GENE_A", "brain", -0.02),
        ("GENE_B", "liver", 0.12),
        ("GENE_B", "brain", 0.01),
        ("GENE_C", "liver", 0.0),
        ("GENE_C", "brain", -0.01),
    ]


def test_aggregate_and_predict_target():
    rows = aggregate(fake_scorer("chr21", 1, "AC", "A"))
    assert rows[0]["gene"] == "GENE_A" and rows[0]["max_drop_tissue"] == "liver"
    p = predict_target(rows)
    assert p["gene"] == "GENE_A" and p["action"] == "activates" and p["strength"] == "strong"
    assert p["confidence"] == 0.45 and p["tissue"] == "liver"
    # below the threshold nothing is named
    assert predict_target(rows, min_effect=0.5) is None
    # a silencer: only the rise survives a threshold above the drop... use a rows set where rise wins
    rows2 = aggregate([("GENE_B", "liver", 0.2), ("GENE_B", "brain", -0.05)])
    assert predict_target(rows2)["action"] == "represses"


def test_compare_verdicts():
    inferred = [{"gene": "GENE_A", "distance": 10_000}, {"gene": "GENE_B", "distance": 50_000}]
    assert compare(None, inferred) == "no predicted effect"
    assert compare({"gene": "GENE_A"}, inferred) == "agrees with nearest TSS in domain"
    assert compare({"gene": "GENE_B"}, inferred) == "another gene in the same domain"
    assert compare({"gene": "LNC1"}, inferred, domain_genes=["LNC1"]) == "another gene in the same domain"
    assert compare({"gene": "FAR"}, inferred) == "gene outside the domain"


def test_score_element_caches_and_summarises(tmp_path):
    calls = []

    def scorer(chrom, pos, ref, alt):
        calls.append((chrom, pos, ref, alt))
        return fake_scorer(chrom, pos, ref, alt)

    def fetch(locus: Locus) -> str:
        return "G" + "ACGT" * 5  # anchor base then a 20 bp element

    inferred = [{"gene": "GENE_A", "distance": 12_345}]
    r = score_element(scorer, fetch, "chr21", "EH1", 100, 120, inferred, cache=tmp_path)
    assert calls == [("chr21", 100, "GACGTACGTACGTACGTACGT", "G")]
    assert r["predicted"]["gene"] == "GENE_A" and r["verdict"] == "agrees with nearest TSS in domain"
    assert r["length"] == 20 and r["genes_in_window"] == 3
    # second call answers from the cache
    r2 = score_element(scorer, fetch, "chr21", "EH1", 100, 120, inferred, cache=tmp_path)
    assert len(calls) == 1 and r2["predicted"] == r["predicted"]
    s = summarise([r, {"predicted": None, "verdict": "no predicted effect"}])
    assert s["elements_scored"] == 2 and s["with_predicted_target"] == 1
    assert s["fraction_agreeing_with_nearest"] == 1.0 and s["top_tissues"] == {"liver": 1}
    assert s["median_distance_when_agreeing"] == 12_345


def test_cached_prediction_reads_only_the_cache(tmp_path):
    from genomeos.predict.enhancer_target import cached_prediction

    assert cached_prediction("chr21", "EH-missing", cache=tmp_path) is None
    score_element(fake_scorer, lambda loc: "G" + "A" * 10, "chr21", "EH2", 50, 60, cache=tmp_path)
    p = cached_prediction("chr21", "EH2", cache=tmp_path)
    assert p["gene"] == "GENE_A" and p["action"] == "activates"
