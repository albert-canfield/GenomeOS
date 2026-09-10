"""Phase 5: the early C. elegans lineage from a zygote bootstrap and cited cycle times."""

from genomeos.organism import run_lineage


def test_early_lineage_matches_sulston():
    lin = run_lineage(until_min=150)
    assert {c.name for c in lin.alive_at(0.5)} == {"AB", "P1"}  # 2-cell stage
    four = {c.name for c in lin.alive_at(23)}
    assert four == {"ABa", "ABp", "EMS", "P2"}, four  # 4-cell stage after ~22 min
    founders = {"AB", "MS", "E", "C", "D", "P4"}
    assert founders <= set(lin.cells)
    assert lin.cells["E"].parent == "EMS" and lin.cells["MS"].parent == "EMS"
    assert lin.cells["P4"].parent == "P3" and lin.cells["D"].parent == "P3"
    # 26–28 cells around 100 min (Sulston: 28-cell stage ≈ 100 min at 20 °C)
    assert 22 <= lin.count_at(100) <= 34, lin.count_at(100)
    # AB divides fastest, E slowest: at 120 min AB descendants outnumber E descendants many-fold
    assert lin.lineage_count_at("AB", 120) >= 8 * lin.lineage_count_at("E", 120)
    assert lin.lineage_count_at("E", 120) in (2, 4)
    assert lin.module is not None and lin.module.parameters["ab_cycle_min"].evidence.source.startswith(
        "Sulston"
    )
    tree = lin.tree(depth=3)
    assert tree[0].startswith("P0") and any(line.strip().startswith("EMS") for line in tree)
