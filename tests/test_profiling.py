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


def test_species_masks_and_mask_profile_match_the_plain_profile():
    from genomeos.knowledge.profiling import mask_profile, set_mask, species_masks

    presence = {"species": ["a", "b", "c", "d"], "genes": {"G1": "f", "G2": "3", "G3": "1"}, "strata": {}}
    species, masks, index = species_masks(presence)
    assert species == ["a", "b", "c", "d"] and set(index) == {"G1", "G2", "G3"}
    s_ = set_mask(["G1", "G2", "G3"], index)
    assert [round(x, 4) for x in mask_profile(s_, masks)] == profile(["f", "3", "1"], 4)  # profile() rounds
    assert mask_profile(s_, masks, keep=[2, 3]) == [2 / 3, 1.0]


def test_pair_null_separates_shared_history_from_shared_age():
    """Two libraries whose genes were lost together in one clade stand out against age-matched random
    sets of the same stratum; two libraries of ordinary genes of the same age do not."""
    import random

    from genomeos.knowledge.profiling import pair_null

    rng = random.Random(1)
    n_species = 40
    species = [f"s{i}" for i in range(n_species)]
    strata = {sp: ("Life" if i < 10 else "Metazoa") for i, sp in enumerate(species)}

    def mask(present):
        bits = "".join("1" if i in present else "0" for i in range(n_species))
        return format(int(bits, 2), "010x")

    genes, origin = {}, {}
    # 400 ordinary animal genes, present in a random 80% of animal species
    for k in range(400):
        present = {i for i in range(10, n_species) if rng.random() < 0.8}
        genes[f"N{k}"] = mask(present)
        origin[f"N{k}"] = {"origin": "Metazoa"}
    # 30 genes lost together in species 25..39 (one clade), split between two libraries
    for k in range(30):
        present = {i for i in range(10, 25) if rng.random() < 0.95}
        genes[f"L{k}"] = mask(present)
        origin[f"L{k}"] = {"origin": "Metazoa"}
    presence = {"species": species, "strata": strata, "genes": genes}
    members = {
        "lost.a": [f"L{k}" for k in range(15)],
        "lost.b": [f"L{k}" for k in range(15, 30)],
        "plain.a": [f"N{k}" for k in range(15)],
        "plain.b": [f"N{k}" for k in range(15, 30)],
        "copy.a": [f"N{k}" for k in range(40, 55)],
        "copy.b": [f"N{k}" for k in range(40, 55)] + ["N60"],
    }
    s = pair_null(presence, {"genes": origin}, members, min_members=8, draws=150, seed=3)
    by = {tuple(r["libraries"]): r for r in s["pairs"]}
    lost = by[("lost.a", "lost.b")]
    assert lost["z"] > 3 and lost["q"] <= 0.05 and lost["r_exclusive_eukaryotes"] is not None
    plain = by[("plain.a", "plain.b")]
    assert abs(plain["z"]) < 3 and plain["q"] > 0.05
    copy = by[("copy.a", "copy.b")]
    assert "note" in copy and copy["shared_members"] == 15  # the same genes under two names: not tested
    assert ("lost.a", "lost.b") in {tuple(r["libraries"]) for r in s["significant_q05"]}
