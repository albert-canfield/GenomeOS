from pathlib import Path

import pytest

from genomeos.ir import Action
from genomeos.knowledge import LigandReceptorTable

INTER = Path("data/knowledge/cellphonedb_interaction_input.csv")
GENES = Path("data/knowledge/cellphonedb_gene_input.csv")


@pytest.mark.skipif(not (INTER.exists() and GENES.exists()), reason="download CellPhoneDB files first")
def test_cellphonedb_pairs_and_bioir():
    t = LigandReceptorTable.from_cellphonedb(INTER, GENES)
    assert len(t) > 1500
    pairs = {(a, b) for a, b, _, _ in t.pairs} | {(b, a) for a, b, _, _ in t.pairs}
    for a, b in (("CXCL12", "CXCR4"), ("EGF", "EGFR"), ("DLL4", "NOTCH1"), ("TGFB1", "TGFBR2")):
        assert (a, b) in pairs, (a, b)
    m = t.to_module()
    assert all(r.action is Action.BIND for r in m.rules)
    assert len(m.rules) > 1500
    assert m.confidence_report()["rule"] == 0.8
    assert any(p == "CXCR4" for p, _ in t.partners("CXCL12"))
