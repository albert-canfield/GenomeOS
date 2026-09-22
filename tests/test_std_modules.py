"""The standard library modules added for human development: Carnegie stages, signalling, turnover; and
the protein import from the packaged proteome."""

from pathlib import Path

import pytest

from genomeos.knowledge.ligand_receptor import DEVELOPMENTAL, LigandReceptorTable, to_std_signalling
from genomeos.lang import BioLangError, parse
from genomeos.lang.parser import IMPORT_RESOLVERS
from genomeos.organism.human import to_std_turnover
from genomeos.results import load_result


def test_human_stages_cover_the_embryonic_period_in_order():
    m = parse("module t\nimport bio.std.human_stages\n")
    names = [s.name for s in m.stages]
    assert names[:3] == ["CS1", "CS2", "CS3"] and names[-2:] == ["CS23", "Fetal"] and len(names) == 24
    assert m.stage_at(16 * 1440) == "CS7" and m.stage_at(100 * 1440) == "Fetal"
    assert all(s.evidence.source for s in m.stages) and m.timer("human_segmentation_clock").minutes() == 300


def test_signalling_module_is_generated_from_cellphonedb():
    m = parse("module t\nimport bio.std.signalling\n")
    sigs = m.signals()
    assert len(sigs) > 100 and all(s.mode == "contact" and s.ligand and s.receptor and s.sets for s in sigs)
    ids = {s.id for s in sigs}
    assert "DLL4_NOTCH1" in ids and "WNT3A_FZD1" in ids and "EPO_EPOR" in ids and "NODAL_ACVR2B" in ids
    notch = next(s for s in sigs if s.id == "DLL4_NOTCH1")
    assert (
        notch.sets == "Notch_signal"
        and "CellPhoneDB" in notch.evidence.source
        and "Bray" in notch.evidence.source
    )
    files = ("data/knowledge/cellphonedb_interaction_input.csv", "data/knowledge/cellphonedb_gene_input.csv")
    if all(Path(f).exists() for f in files):
        table = LigandReceptorTable.from_cellphonedb(*files)
        assert to_std_signalling(table) == Path("genomeos/std/signalling.bio").read_text()
    assert all(v[1] for v in DEVELOPMENTAL.values())  # every pathway carries its review


def test_turnover_module_is_generated_from_sender_and_milo():
    m = parse("module t\nimport bio.std.human_turnover\n")
    timers = {t.name: t for t in m.timers}
    assert (
        timers["Erythrocytes_lifespan"].minutes() == 119 * 1440
        and timers["Neutrophils_lifespan"].duration == 6.6
    )
    assert all(t.evidence.kind.value == "experimental" for t in m.timers)
    r = load_result("human_cell_turnover")
    if r is not None:
        assert to_std_turnover(r) == Path("genomeos/std/human_turnover.bio").read_text()


def test_protein_import_resolves_from_the_packaged_proteome():
    from genomeos.lib.biolang import register

    register()
    assert "protein" in IMPORT_RESOLVERS
    m = parse(
        "module t\nimport protein:TP53\nimport protein:GATA1\nimport protein:TP53\n"
    )  # repeated import merges once
    p = m.entities["TP53"]
    assert p.kind == "protein" and p.accession == "P04637" and p.domains and p.pathways
    with pytest.raises(ValueError, match="not in the packaged proteome"):
        parse("module t\nimport protein:NOSUCHGENE\n")
    with pytest.raises(BioLangError, match="no resolver registered"):
        parse("module t\nimport nothing:here\n")
