"""Task 1.5: GO / Reactome import and data-driven library membership."""

import pytest

from genomeos.knowledge import Ontology
from genomeos.lib import LIBRARIES, KnowledgeBase

needs = pytest.mark.skipif(not KnowledgeBase.available(), reason="download knowledge files first")


def test_obo_parser_on_small_ontology(tmp_path):
    obo = tmp_path / "t.obo"
    obo.write_text(
        "format-version: 1.2\n\n[Term]\nid: A:1\nname: root\n\n"
        "[Term]\nid: A:2\nname: child\nis_a: A:1 ! root\n\n"
        "[Term]\nid: A:3\nname: part\nrelationship: part_of A:2\n\n[Typedef]\nid: part_of\n"
    )
    o = Ontology.from_obo(obo)
    assert len(o) == 3
    assert o.ancestors("A:3") == {"A:3", "A:2", "A:1"}
    assert o.descendants("A:1") == {"A:1", "A:2", "A:3"}
    assert o.search("chi")[0].id == "A:2"


@needs
def test_every_catalogue_go_term_exists_and_is_current():
    kb = KnowledgeBase()
    for lib in LIBRARIES.values():
        for t in lib.go_terms:
            assert t in kb.go, (lib.id, t)
            assert not kb.go.terms[t].obsolete, (lib.id, t, kb.go.name(t))


@needs
def test_catalogue_genes_are_members_by_data():
    kb = KnowledgeBase()
    results = kb.verify_all()
    assert len(results) >= 40
    missing = [(v.library, v.missing) for v in results if v.missing]
    total = sum(len(v.checked) for v in results)
    miss = sum(len(v.missing) for v in results)
    # Agreement is ~95%: the misses are developmental master regulators and stem-cell markers,
    # which human GO annotation records only as "regulation of transcription" (see docs/PROGRESS.md).
    assert miss / total < 0.06, missing
    assert not any(v.unknown_terms for v in results)


@needs
def test_membership_counts_are_plausible():
    kb = KnowledgeBase()
    assert 150 < len(kb.members(LIBRARIES["core.replication"])) < 600
    assert len(kb.members(LIBRARIES["blueprint.signalling_toolkit"])) > 1000
    assert "TP53" in kb.libraries_of("TP53") or "core.dna_repair" in kb.libraries_of("TP53")
    assert "core.dna_repair" in kb.libraries_of("BRCA1")
    assert "timer.circadian" in kb.libraries_of("PER2")
