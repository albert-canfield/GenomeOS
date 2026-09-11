# SPDX-License-Identifier: AGPL-3.0-or-later
"""Node comparison across species by symbol (offline)."""

import json

from genomeos.genome.mouse import compare_nodes, human_nodes_by_symbol


def test_compare_nodes_by_symbol(tmp_path):
    (tmp_path / "domains_chr21.json").write_text(
        json.dumps(
            {
                "domains": [
                    {"id": "chr21:D1", "genes": ["APP", "APP-DT"]},
                    {"id": "chr21:D2", "genes": ["SOD1", "SFRS15"]},
                ]
            }
        )
    )
    index, n = human_nodes_by_symbol(tmp_path)
    assert n == 2 and index["APP"] == "chr21:D1" and index["SOD1"] == "chr21:D2"
    index["FAR"] = "chr5:D9"
    mouse = [
        {"id": "chr16:D1", "start": 0, "end": 10, "coding_symbols": ["App", "Sod1"]},  # adjacent human nodes
        {"id": "chr16:D2", "start": 10, "end": 20, "coding_symbols": ["Sod1", "Sfrs15"]},  # one node
        {"id": "chr16:D3", "start": 20, "end": 30, "coding_symbols": ["App", "Gm1234"]},  # one match only
        {"id": "chr16:D4", "start": 30, "end": 40, "coding_symbols": []},
        {"id": "chr16:D5", "start": 40, "end": 50, "coding_symbols": ["App", "Far"]},  # another chromosome
    ]
    c = compare_nodes(mouse, index)
    assert (c["conserved"], c["split_adjacent"], c["split_scattered"], c["unmapped"]) == (1, 1, 1, 2)
    assert c["tested"] == 3 and c["fraction_conserved"] == 0.333 and c["fraction_same_neighbourhood"] == 0.667
    verdicts = [r["verdict"] for r in c["rows"]]
    assert verdicts == ["split_adjacent", "conserved", "unmapped", "unmapped", "split_scattered"]
    assert c["orthology"].startswith("inferred")
    # a curated orthology maps a renamed mouse gene onto its human orthologue
    c2 = compare_nodes(mouse, index, {"Gm1234": ["APP-DT"], "Far": ["SOD1"]})
    assert c2["orthology"].startswith("curated")
    assert [r["verdict"] for r in c2["rows"]] == [
        "split_adjacent",
        "conserved",
        "conserved",
        "unmapped",
        "split_adjacent",
    ]
