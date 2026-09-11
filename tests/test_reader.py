from types import SimpleNamespace

from genomeos.coords import Locus, Strand
from genomeos.genome.domains import Domain
from genomeos.genome.reader import PeakIndex, compare, load_peaks, read_chromosome, slug
from genomeos.genome.regulatory import CCRE


class _Ann:
    def __init__(self, genes):
        self.genes = {g.id: g for g in genes}


def _gene(gid, sym, start, end, strand=Strand.PLUS, typ="protein_coding"):
    return SimpleNamespace(id=gid, symbol=sym, type=typ, locus=Locus("c", start, end, strand))


def test_peak_index_and_reader(tmp_path, monkeypatch):
    from genomeos.genome import reader

    monkeypatch.setattr(reader, "RESULTS", tmp_path)
    import gzip

    with gzip.open(tmp_path / f"dnase_{slug('K562')}_c.bed.gz", "wt") as fh:
        fh.write("# test\n900\t1200\t5.0\n50000\t50400\t2.0\n119900\t120300\t9.0\n")
    peaks = load_peaks("K562", "c")
    idx = PeakIndex(peaks)
    assert idx.covered_bp(1000, 1100) == 100 and idx.overlapping(2000, 3000) == []
    ann = _Ann(
        [
            _gene("g1", "A", 1000, 5000),
            _gene("g2", "B", 30000, 40000, Strand.MINUS),
            _gene("g3", "C", 119000, 125000),
            _gene("g4", "NC", 1000, 2000, typ="lncRNA"),
        ]
    )
    domains = [
        Domain("c:D1", "c", 0, 60000, coding_genes=2),
        Domain("c:D2", "c", 60000, 100000),
        Domain("c:D3", "c", 100000, 200000, coding_genes=1),
    ]
    ccres = [CCRE("c", 50100, 50300, "e1", "dELS", False), CCRE("c", 70000, 70200, "e2", "pELS", False)]
    r = read_chromosome("K562", "c", ann, domains, ccres)
    assert r["genes_read"] == 2 and r["genes_silent"] == 1 and r["silent_genes"] == ["B"]
    assert r["top_read"][0]["gene"] == "C" and r["top_read"][0]["signal"] == 9.0
    assert r["enhancers_active"] == 1 and r["enhancers"] == 2
    assert r["nodes_silent"] == 1 and r["silent_node_ids"] == ["c:D2"]
    other = dict(
        r,
        cell_type="HepG2",
        top_read=[{"gene": "B", "signal": 1.0}],
        _read_all=["B"],
        silent_genes=["A", "C"],
    )
    c = compare(r, other)
    assert c["read_in_a_only"] == ["A", "C"] and c["read_in_b_only"] == ["B"] and c["read_in_both"] == 0
