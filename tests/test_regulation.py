from types import SimpleNamespace

from genomeos.coords import Locus, Strand
from genomeos.genome.domains import Domain
from genomeos.genome.regulation import assign_targets, regulation_of
from genomeos.genome.regulatory import CCRE
from genomeos.lang.parser import parse


def _gene(gid, symbol, start, end, strand=Strand.PLUS, typ="protein_coding"):
    return SimpleNamespace(
        id=gid, symbol=symbol, type=typ, locus=Locus("c", start, end, strand), transcripts={}
    )


class _Ann:
    def __init__(self, genes):
        self.genes = {g.id: g for g in genes}
        self.by_symbol = {g.symbol: g.id for g in genes}

    def gene(self, s):
        return self.genes[self.by_symbol.get(s, s)]


def _setup():
    genes = [
        _gene("g1", "A", 10_000, 30_000),  # TSS 10,000, domain 1
        _gene("g2", "B", 60_000, 90_000, Strand.MINUS),  # TSS 89,999, domain 1
        _gene("g3", "C", 150_000, 170_000),  # TSS 150,000, domain 2
        _gene("g4", "NC", 12_000, 13_000, typ="lncRNA"),
    ]
    ann = _Ann(genes)
    ccres = [
        CCRE("c", 9_600, 9_900, "E_prom", "PLS", False),  # promoter of A
        CCRE("c", 40_000, 40_300, "E_enh1", "dELS", False),  # domain 1, nearest A (30 kb) vs B (50 kb)
        CCRE("c", 20_000, 20_300, "E_enh2", "pELS", False),  # inside A
        CCRE("c", 99_800, 100_200, "E_ctcf", "CTCF-only", True),  # boundary
        CCRE("c", 140_000, 140_300, "E_enh3", "dELS", False),  # domain 2: reaches C only
        CCRE("c", 149_000, 149_300, "E_open", "DNase-H3K4me3", False),  # open chromatin at C's TSS
    ]
    domains = [
        Domain("c:D1", "c", 0, 100_000),
        Domain("c:D2", "c", 100_000, 200_000, genes=["C"], coding_genes=1),
    ]
    return ann, ccres, domains


def test_targets_respect_domain_boundaries():
    ann, ccres, domains = _setup()
    els = {e.id: e for e in assign_targets("c", ccres, ann, domains)}
    assert els["E_prom"].cls == "promoter" and els["E_prom"].targets[0]["gene"] == "A"
    assert els["E_prom"].targets[0]["basis"] == "promoter of"
    enh1 = els["E_enh1"].targets
    assert [t["gene"] for t in enh1] == ["A", "B"] and enh1[0]["basis"] == "nearest TSS in domain"
    assert enh1[1]["basis"] == "same domain" and enh1[1]["confidence"] < enh1[0]["confidence"]
    assert "C" not in {t["gene"] for t in enh1}  # other side of the CTCF boundary
    assert [t["gene"] for t in els["E_enh3"].targets] == ["C"]
    assert els["E_ctcf"].cls == "insulator" and els["E_ctcf"].targets == []
    assert els["E_open"].targets[0]["gene"] == "C" and els["E_open"].cls == "open_chromatin"
    assert els["E_enh1"].domain == "c:D1" and els["E_enh3"].domain == "c:D2"


def test_regulation_of_a_gene():
    ann, ccres, _ = _setup()
    r = regulation_of("A", "c", ccres, ann, 200_000)
    assert [p["id"] for p in r["promoters"]] == ["E_prom"]
    assert r["enhancers_in_domain"] == 2 and r["enhancers_inside_gene"] == 1
    assert r["enhancers"][0]["id"] == "E_enh2" and r["enhancers"][0]["intragenic"]
    assert r["domain"] and r["domain"]["coding_genes"] == 2
    assert "B" in r["competing_genes"] and "C" not in r["competing_genes"]
    assert len(r["insulators_bounding"]) == 1


def test_biolang_element_block():
    m = parse(
        "module t\nelement APP_enh1 { class: enhancer; locus: chr21:25900000-25900400(+); "
        'targets: APP, CYYR1; domain: chr21:D89; evidence: inferred "same CTCF domain"; confidence: 0.4 }\n'
    )
    e = m.entities["APP_enh1"]
    assert e.kind == "regulatory_element" and e.cls == "enhancer" and e.domain == "chr21:D89"
    assert [t["gene"] for t in e.targets] == ["APP", "CYYR1"] and e.locus.start == 25900000
