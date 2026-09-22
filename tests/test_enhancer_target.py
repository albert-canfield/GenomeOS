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


def test_aggregate_keeps_the_cell_lines_own_tracks():
    from genomeos.predict.enhancer_target import aggregate, by_cell_of, has_cells

    rows = aggregate([("G", "K562", -0.5), ("G", "liver", -0.9), ("G", "IMR-90", 0.1), ("H", "HepG2", 0.2)])
    g = next(r for r in rows if r["gene"] == "G")
    assert g["by_cell"] == {"K562": -0.5, "IMR-90": 0.1} and g["max_drop_tissue"] == "liver"
    assert by_cell_of(rows, {"gene": "G"}) == {"K562": -0.5, "IMR-90": 0.1}
    assert by_cell_of(rows, {"gene": "Z"}) is None and by_cell_of(rows, None) is None
    assert has_cells({"genes": rows}) and not has_cells({"genes": [{"gene": "old"}]}) and not has_cells(None)


def test_pack_folds_a_chromosome_and_load_cached_reads_the_archive(tmp_path):
    import json

    from genomeos.predict.enhancer_target import archive_path, cache_path, load_cached, pack

    for eid, gene in (("E1", "APP"), ("E2", "SOD1")):
        p = cache_path("chr21", eid, tmp_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps({"id": eid, "genes": [{"gene": gene, "by_cell": {}}]}))
    r = pack("chr21", tmp_path, remove=True)
    assert r["elements"] == 2 and r["packed"] == 2 and r["removed"] == 2
    assert archive_path("chr21", tmp_path).exists() and not cache_path("chr21", "E1", tmp_path).exists()
    assert load_cached("chr21", "E1", tmp_path)["genes"][0]["gene"] == "APP"
    assert load_cached("chr21", "E9", tmp_path) is None
    # a per-element file written after the pack still wins, and a second pack absorbs it
    cache_path("chr21", "E3", tmp_path).parent.mkdir(parents=True, exist_ok=True)
    cache_path("chr21", "E3", tmp_path).write_text(json.dumps({"id": "E3", "genes": []}))
    assert load_cached("chr21", "E3", tmp_path)["id"] == "E3"
    assert pack("chr21", tmp_path, remove=True)["elements"] == 3
