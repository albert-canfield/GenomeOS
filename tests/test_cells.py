from pathlib import Path

import pytest

from genomeos.knowledge import CellTypes

CL = Path("data/knowledge/cl-basic.obo")


@pytest.mark.skipif(not CL.exists(), reason="download cl-basic.obo first")
def test_cell_ontology_loads_and_compiles():
    cells = CellTypes.from_obo(CL)
    assert len(cells) > 2500
    hits = cells.search("cardiac muscle cell")
    assert hits and hits[0].id.startswith("CL:")
    cm = next(t for t in hits if t.name == "cardiac muscle cell")
    lineage = cells.lineage(cm.id)
    assert any("muscle cell" in x for x in lineage) and any("cell" in x for x in lineage[-1:])
    m = cells.to_module(limit=100)
    assert len(m.entities) == 100
    assert all(e.kind == "cell_type" for e in m.entities.values())
