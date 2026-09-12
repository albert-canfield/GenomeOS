"""Phylogenetic profiling between libraries from the presence matrix (area J step 5)."""

from genomeos.knowledge.profiling import build, gain_curve, pearson, profile, unpack


def test_unpack_profile_pearson_and_gain_curve():
    assert unpack("7", 3) == [1, 1, 1] and unpack("1", 3) == [0, 0, 1] and unpack("", 3) == [0, 0, 0]
    assert profile(["7", "1"], 3) == [0.5, 0.5, 1.0] and profile([], 2) == [0.0, 0.0]
    assert pearson([1, 2, 3], [2, 4, 6]) == 1.0 and pearson([1, 2, 3], [3, 2, 1]) == -1.0
    assert pearson([1, 1, 1], [1, 2, 3]) is None and pearson([1, 2], [1, 2]) is None
    gc = gain_curve(["Eukaryota", "Metazoa", "Metazoa", "Mammalia"])
    assert (
        gc["Eukaryota"] == 0.25 and gc["Metazoa"] == 0.75 and gc["Vertebrata"] == 0.75 and gc["Homo"] == 1.0
    )


def test_build_finds_libraries_with_the_same_history():
    # species deep to shallow: yeast, fly, fish, mouse; masks are 4 bits
    presence = {
        "species": ["yeast", "fly", "fish", "mouse"],
        "strata": {
            "yeast": "Opisthokonta",
            "fly": "Bilateria",
            "fish": "Euteleostomi",
            "mouse": "Euarchontoglires",
        },
        "genes": {
            **{f"old{i}": "f" for i in range(8)},  # present everywhere
            **{f"vert{i}": "3" for i in range(8)},  # fish and mouse only
            **{f"mam{i}": "1" for i in range(8)},  # mouse only
        },
    }
    origin = {
        "genes": {
            **{f"old{i}": {"origin": "Opisthokonta"} for i in range(8)},
            **{f"vert{i}": {"origin": "Euteleostomi"} for i in range(8)},
            **{f"mam{i}": {"origin": "Euarchontoglires"} for i in range(8)},
        }
    }
    members = {
        "core.a": [f"old{i}" for i in range(8)],
        "core.b": [f"old{i}" for i in range(4)] + [f"old{i}" for i in range(4, 8)],
        "systems.v": [f"vert{i}" for i in range(8)],
        "systems.m": [f"mam{i}" for i in range(8)],
        "tiny": ["old0"],
    }
    s = build(presence, origin, members, min_members=8)
    assert "tiny" not in s["libraries"] and s["species"] == 4
    assert s["libraries"]["core.a"]["profile"] == [1.0, 1.0, 1.0, 1.0]
    assert s["libraries"]["systems.v"]["gain_curve"]["Euteleostomi"] == 1.0
    assert s["libraries"]["systems.v"]["gain_curve"]["Bilateria"] == 0.0
    top = s["top_pairs"][0]
    assert set(top["libraries"]) == {"core.a", "core.b"} and top["residual_correlation"] == 1.0
    bottom = s["bottom_pairs"][-1]
    assert bottom["residual_correlation"] < 0
