"""The Ma 2021 transcription-factor atlas: distil, call presence per factor, generate the reader, check it
against textbook tissue relations; and the `express` action that gives cells their measured factors."""

import io
import json
import zipfile
from pathlib import Path

import pytest

from genomeos.lang import parse
from genomeos.organism.reference import ReferenceLineage, parse_wormweb
from genomeos.organism.reference import distil as distil_lineage
from genomeos.organism.tf_atlas import (
    CELLS_FILE,
    cells_expressing,
    distil,
    load_cells,
    save_cells,
    summary,
    textbook_check,
    to_bio_reader,
)
from genomeos.runtime.body import Body


def _archive() -> bytes:
    """Two factors, one with two strains, over three cells and three frames."""
    rows = {
        "ELT-2_S1.csv": [("Ea", 1, 0.0), ("Ea", 2, 50.0), ("Ea", 3, 120.0), ("Ep", 3, 30.0), ("ABa", 3, 5.0)],
        "ELT-2_S2.csv": [("Ea", 1, 0.0), ("Ep", 2, 90.0), ("ABa", 2, 0.0)],
        "HLH-1_S9.csv": [("ABa", 1, 2.0), ("Ea", 2, 0.3), ("Ep", 1, 0.0)],
    }
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        for name, rs in rows.items():
            text = "TF_strain,Time,Cell-name,Raw-expression,Adjustment-expression\n" + "".join(
                f"{name[:-4]},{t},{c},{v},{v}\n" for c, t, v in rs
            )
            z.writestr(name, text)
    return buf.getvalue()


def test_distil_calls_presence_per_factor_and_merges_strains():
    table = distil(_archive())
    assert set(table) == {"ELT-2", "HLH-1"} and table["ELT-2"]["strains"] == ["ELT-2_S1", "ELT-2_S2"]
    assert table["ELT-2"]["max"] == 120.0 and table["ELT-2"]["threshold"] == 24.0
    assert set(table["ELT-2"]["cells"]) == {"Ea", "Ep"}  # ABa at 5 is under a fifth of the maximum
    assert table["ELT-2"]["cells"]["Ea"] == [120.0, 2] and table["ELT-2"]["cells"]["Ep"] == [90.0, 2]
    assert set(table["HLH-1"]["cells"]) == {"ABa"} and table["HLH-1"]["cells"]["ABa"] == [2.0, 1]
    assert cells_expressing(table) == {"ABa": ["HLH-1"], "Ea": ["ELT-2"], "Ep": ["ELT-2"]}
    bio = to_bio_reader(table)
    assert "decision read_Ea { action: express; when: cell = Ea; sets: ELT-2;" in bio
    assert 'evidence: experimental "Ma et al. 2021' in bio and bio.count("action: express") == 3


def test_compact_table_round_trips(tmp_path):
    table = distil(_archive())
    path = save_cells(table, tmp_path / "cells.json")
    data = json.loads(path.read_text())
    assert data["factors"] == ["ELT-2", "HLH-1"] and data["cells"]["Ea"] == [[0, 2]]
    assert load_cells(path) == {"ABa": ["HLH-1"], "Ea": ["ELT-2"], "Ep": ["ELT-2"]}
    s = summary(table, [])
    assert s["factors"] == 2 and s["strains"] == 3 and s["cells_with_a_factor"] == 3


def test_textbook_check_scores_precision_and_recall_against_the_lineage():
    toy = {
        "id": "p0",
        "name": "P0",
        "did": "P0",
        "data": {"levelDistance": 0, "totalDistance": 0, "deathDistance": 0, "type": ""},
        "children": [
            {
                "id": "e",
                "name": "E",
                "did": "E",
                "data": {"levelDistance": 10, "totalDistance": 10, "deathDistance": 0, "type": ""},
                "children": [
                    {
                        "id": "ea",
                        "name": "Ea",
                        "did": "Ea",
                        "data": {
                            "levelDistance": 5,
                            "totalDistance": 15,
                            "deathDistance": 0,
                            "type": "intestine",
                        },
                    },
                    {
                        "id": "ep",
                        "name": "Ep",
                        "did": "Ep",
                        "data": {
                            "levelDistance": 5,
                            "totalDistance": 15,
                            "deathDistance": 0,
                            "type": "intestine",
                        },
                    },
                ],
            },
            {
                "id": "aba",
                "name": "ABa",
                "did": "ABa",
                "data": {"levelDistance": 10, "totalDistance": 10, "deathDistance": 0, "type": "muscle"},
            },
        ],
    }
    ref = ReferenceLineage.from_dict(distil_lineage(parse_wormweb("var json = " + json.dumps(toy) + ";")))
    checks = textbook_check(distil(_archive()), ref, None)
    by = {c["factor"]: c for c in checks}
    assert by["ELT-2"]["precision"] == 1.0 and by["ELT-2"]["recall"] == 1.0
    assert by["HLH-1"]["precision"] == 1.0 and by["HLH-1"]["recall"] == 1.0
    assert by["PHA-4"]["status"] == "factor not in atlas"


def test_express_gives_cells_their_measured_factors_and_drops_inherited_ones():
    src = """
module toy.express
import bio.std.development
organism E { root: R; cell_type: Zygote; factors: M }
cell_type Gut { parent: Blastomere }
decision first { action: divide; when: cell = R; daughters: A, B; after: 1 min }
decision read_A { action: express; when: cell = A; sets: ELT-2, END-1 }
decision read_B { action: express; when: cell = B; sets: HLH-1 }
decision second { action: divide; when: cell = A; daughters: Aa, Ap; after: 1 min }
decision read_Aa { action: express; when: cell = Aa; sets: ELT-2 }
decision gut { action: differentiate; when: ELT-2 = present, END-1 = absent; to: Gut }
"""
    b = Body(parse(src)).run(until=10)
    assert sorted(b.cells["A"].factors) == ["ELT-2", "END-1", "M"] and b.cells["B"].factors.get("HLH-1")
    assert (
        sorted(b.cells["Aa"].factors) == ["ELT-2", "M"] and b.cells["Aa"].cell_type == "Gut"
    )  # END-1 dropped
    assert sorted(b.cells["Ap"].factors) == ["ELT-2", "END-1", "M"]  # not covered by the reader: inherits
    k = Body(parse(src), knockouts={"ELT-2"}).run(until=10)
    assert "ELT-2" not in k.cells["Aa"].factors and k.cells["Aa"].cell_type == "Zygote"
    bad = "module t\nimport bio.std.development\norganism O { root: R; cell_type: Zygote }\n"
    with pytest.raises(Exception, match="express needs"):
        parse(bad + "decision x { action: express; when: cell = R }")


@pytest.mark.skipif(
    not Path(CELLS_FILE).exists(), reason="run `genomeos data distil --only celegans_tf_atlas`"
)
def test_distilled_atlas_agrees_with_the_textbook():
    from genomeos.results import load_result

    r = load_result("celegans_tf_atlas")
    assert r["factors"] == 266 and r["strains"] == 291 and r["cells_with_a_factor"] > 1100
    by = {c["factor"]: c for c in r["textbook"]}
    assert by["ELT-2"]["precision"] > 0.9 and by["ELT-2"]["recall"] == 1.0
    assert (
        by["ELT-7"]["precision"] > 0.95
        and by["HLH-1"]["precision"] > 0.9
        and by["UNC-120"]["precision"] > 0.9
    )
    assert (
        0.4 < by["HLH-1"]["recall"] < 0.7
    )  # HLH-1 marks body-wall muscle, not the pharyngeal and other muscles
    assert by["PHA-4"]["precision"] > 0.8
    cells = load_cells()
    # END-1 comes on at the 4E stage, ELT-2 at 8E, as the genetics says; E itself carries neither, nor PHA-4
    assert all("END-1" in cells[c] for c in ("Eal", "Ear", "Epl", "Epr"))
    assert all("ELT-2" in cells[c] for c in ("Ealaa", "Eprpp")) and "PHA-4" not in cells["E"]
    assert Path("data/organisms/celegans/reader.bio").exists()
