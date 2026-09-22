"""CTCF motif orientation at inferred node boundaries (genome/domains.py)."""

from __future__ import annotations

from genomeos.genome.domains import (
    boundaries_from_ccres,
    boundary_clusters,
    domain_edges,
    infer_domains,
    orient_boundaries,
)
from genomeos.genome.regulatory import CCRE


def _ctcf(pos: int, name: str) -> CCRE:
    return CCRE("chrT", pos - 100, pos + 100, name, "CTCF-only", True)


def test_clusters_and_edges_reproduce_infer_domains():
    ccres = [_ctcf(p, f"E{i}") for i, p in enumerate([100_000, 103_000, 200_000, 220_000, 400_000, 600_000])]
    clusters = boundary_clusters(ccres)
    assert [p for p, _ in clusters] == boundaries_from_ccres(ccres)
    assert clusters[0][1] == ["E0", "E1"]  # merged within 5 kb
    edges = domain_edges(800_000, [p for p, _ in clusters])
    doms = infer_domains("chrT", 800_000, ccres)
    assert edges == [doms[0].start, *[d.end for d in doms]]
    assert 220_000 not in edges  # closer than the minimum domain to the previous edge


def test_convergent_divergent_and_unmarked_edges():
    # edges at 200k, 400k, 600k, 800k on a 1 Mb chromosome
    ccres = [_ctcf(p, f"B{p // 1000}") for p in (200_000, 400_000, 600_000, 800_000)]
    strands = {"B200": {"+"}, "B400": {"-"}, "B600": {"-"}, "B800": {"+"}}
    rows = {r["position"]: r for r in orient_boundaries(1_000_000, ccres, strands)}
    # domain 200k-400k: + then - -> convergent; domain 600k-800k: - then + -> divergent
    assert rows[200_000]["class"] == "convergent" and rows[400_000]["class"] == "convergent"
    assert rows[600_000]["class"] == "divergent" and rows[800_000]["class"] == "divergent"
    strands["B600"] = set()
    rows = {r["position"]: r for r in orient_boundaries(1_000_000, ccres, strands)}
    assert rows[600_000]["class"] == "none" and rows[800_000]["class"] == "motif"
    assert rows[200_000]["sites"] == "+" and rows[600_000]["sites"] == "none"
    both = {**strands, "B400": {"+", "-"}}
    rows = {r["position"]: r for r in orient_boundaries(1_000_000, ccres, both)}
    assert rows[400_000]["sites"] == "both" and rows[400_000]["class"] == "convergent"


def test_oriented_boundaries_place_reverse_to_forward_flips():
    from genomeos.genome.domains import ORIENTED_EVIDENCE, oriented_boundaries

    sites = [(100_000, "-"), (104_000, "+"), (300_000, "+"), (320_000, "-"), (500_000, "-"), (560_000, "+")]
    ccres = [CCRE("chrT", p - 100, p + 100, f"S{p // 1000}", "dELS", True) for p, _ in sites]
    ccres.append(CCRE("chrT", 700_000, 700_200, "N700", "dELS", False))  # no CTCF support: never a site
    strands = {f"S{p // 1000}": {s} for p, s in sites}
    strands["N700"] = {"+"}
    assert oriented_boundaries(ccres, strands) == [102_000, 530_000]
    strands["S104"] = {"+", "-"}  # both strands: cannot say which way it points, skipped
    assert oriented_boundaries(ccres, strands) == [
        200_000,
        530_000,
    ]  # 100 kb reverse now meets 300 kb forward
    single = {f"S{p // 1000}": {s} for p, s in sites}
    doms = infer_domains("chrT", 1_000_000, ccres, orientation=single)
    assert [d.start for d in doms] == [0, 102_000, 530_000]  # the original strands, both flips
    assert all(d.to_dict()["evidence"] == ORIENTED_EVIDENCE for d in doms)
    assert infer_domains("chrT", 1_000_000, ccres)[0].to_dict()["evidence"] != ORIENTED_EVIDENCE


def test_motif_strands_scan_and_cache(tmp_path):
    from genomeos.genome import motifs as mo
    from genomeos.genome.domains import ctcf_motif_strands

    if not mo.JASPAR_PATH.exists():
        import pytest

        pytest.skip("JASPAR profiles not cached locally")
    forward = "AAAAAGCCACCAGGGGGCGCAAAAA"
    reverse = "TTTTTGCGCCCCCTGGTGGCTTTTT"
    seqs = {0: forward, 1000: reverse, 2000: "A" * 25}
    ccres = [CCRE("chrT", s, s + 25, f"E{s}", "CTCF-only", True) for s in seqs]
    fetch = lambda s, e: seqs[s]  # noqa: E731
    got = ctcf_motif_strands("chrT", ccres, fetch, cache=tmp_path)
    assert got == {"E0": {"+"}, "E1000": {"-"}, "E2000": set()}
    again = ctcf_motif_strands("chrT", ccres, lambda s, e: "", cache=tmp_path)  # read from the cache
    assert again == got


def test_site_scores_serve_every_variant_of_the_call(tmp_path):
    """The stricter site call (registered 2026-09-21): one scan, four ways of reading it."""
    from genomeos.genome import motifs as mo
    from genomeos.genome.domains import STRICT_RELATIVE, ctcf_motif_sites, strand_variant

    if not mo.JASPAR_PATH.exists():
        import pytest

        pytest.skip("JASPAR profiles not cached locally")
    forward = "AAAAAGCCACCAGGGGGCGCAAAAA"
    reverse = "TTTTTGCGCCCCCTGGTGGCTTTTT"
    # a strong forward core with a degenerate reverse-strand hit beside it: the case the loose call
    # reads as both-strand and drops, and the stricter call orients by its best hit
    mixed = "GCCACCAGGGGGCGCTGCGCCCCCTGGTAGC"
    seqs = {0: forward, 1000: reverse, 2000: "A" * 25, 3000: mixed}
    ccres = [CCRE("chrT", s, s + len(q), f"E{s}", "CTCF-only", True) for s, q in seqs.items()]
    sites = ctcf_motif_sites("chrT", ccres, lambda s, e: seqs[s], cache=tmp_path)
    assert ctcf_motif_sites("chrT", ccres, lambda s, e: "", cache=tmp_path) == sites  # cached

    loose = strand_variant(sites)
    assert loose["E0"] == {"+"} and loose["E1000"] == {"-"} and loose["E2000"] == set()
    best = strand_variant(sites, best_hit=True)
    # an element the loose call cannot orient is orientable by its best hit, and never the other way
    assert all(len(v) <= 1 for v in best.values())
    assert {c for c, v in loose.items() if v} == {c for c, v in best.items() if v}
    assert len(loose["E3000"]) == 2 and best["E3000"] == {"+"}

    strict = strand_variant(sites, relative=STRICT_RELATIVE, best_hit=True)
    assert all(strict[c] <= best[c] for c in sites)  # strictness only ever drops a site
    for cid, ((fr, _fs), (rr, _rs)) in sites.items():
        assert bool(strict[cid]) == (max(fr, rr) >= STRICT_RELATIVE)
