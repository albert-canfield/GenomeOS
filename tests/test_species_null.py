# SPDX-License-Identifier: AGPL-3.0-or-later
"""The species-tree null for area J: the tree, the two nulls, and the outcome rules."""

from __future__ import annotations

import importlib.util
import random
import re
from pathlib import Path

import numpy as np
import pytest

from genomeos.knowledge import profiling as pf
from genomeos.knowledge import species_null as sn

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("profiling_genes", ROOT / "scripts" / "profiling_genes.py")
assert spec and spec.loader
pg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pg)


# ------------------------------------------------------------------------------------------- tree


def test_newick_names_lengths_and_internal_labels() -> None:
    root = sn.parse_newick("((A:1,B:2.5)AB:0.5,'C/x=1':3)top;")
    assert root.name == "top"
    ab, c = root.children
    assert (ab.name, ab.length) == ("AB", 0.5)
    assert [(t.name, t.length) for t in ab.children] == [("A", 1.0), ("B", 2.5)]
    assert c.name == "C/x=1" and c.length == 3.0
    assert [t.name for t in sn.leaves(root)] == ["A", "B", "C/x=1"]


@pytest.mark.parametrize("bad", ["((A,B);", "(A,B));", "A,B;"])
def test_newick_refuses_unbalanced_strings(bad: str) -> None:
    with pytest.raises(ValueError):
        sn.parse_newick(bad)


PAN = (
    "((Bact_one,Bact_two)Bacteria,(Saccharomyces_cerevisiae_S288C,"
    "(Drosophila_melanogaster,(Ciona_savignyi,(Homo_sapiens,Danio_rerio)Euteleostomi)Chordata)Bilateria)"
    "Opisthokonta)root;"
)
VERT = (
    "((((((Homo_sapiens,Mus_musculus)Mammalia,Gallus_gallus_reference_breed)Amniota,Danio_rerio)"
    "Euteleostomi,Petromyzon_marinus)Vertebrata,(Ciona_savignyi,Ciona_intestinalis)Ciona)Chordata,"
    "Drosophila_melanogaster)Bilateria;"
)


def test_the_vertebrate_clade_comes_from_the_vertebrate_tree() -> None:
    merged = sn.merge_trees(sn.parse_newick(PAN), sn.parse_newick(VERT))
    tips = {sn.norm(t.name): t for t in sn.leaves(merged)}
    assert {"mus_musculus", "gallus_gallus_reference_breed", "petromyzon_marinus"} <= set(tips)
    vert = sn.mrca([tips["homo_sapiens"], tips["petromyzon_marinus"]])
    assert vert.name == "Vertebrata"
    assert {sn.norm(t.name) for t in sn.leaves(vert)} == {
        "homo_sapiens",
        "mus_musculus",
        "gallus_gallus_reference_breed",
        "danio_rerio",
        "petromyzon_marinus",
    }
    # the outgroup the pan tree lacks is placed where the vertebrate tree places it
    ciona = sn.mrca([tips["ciona_savignyi"], tips["ciona_intestinalis"]])
    assert {sn.norm(t.name) for t in sn.leaves(ciona)} == {"ciona_savignyi", "ciona_intestinalis"}
    assert sn.mrca([tips["homo_sapiens"], tips["bact_one"]]) is merged
    assert all(v.length is None for v in sn.preorder(merged))


def test_a_merge_that_would_need_a_guess_raises() -> None:
    arthropods = "(Strigamia_maritima,Tetranychus_urticae)Arthropoda)Bilateria"
    vert = VERT.replace("Drosophila_melanogaster)Bilateria", arthropods)
    with pytest.raises(ValueError, match="no relative"):
        sn.merge_trees(sn.parse_newick(PAN), sn.parse_newick(vert))


def test_species_match_one_tip_each_and_refuse_to_choose() -> None:
    tips = [
        "Gallus_gallus_reference_breed",
        "Neogale_vison",
        "Aedes_aegypti_strain_LVP_AGWG",
        "Pyrococcus_horikoshii_OT3",
        "Canis_lupus_dingo",
        "Canis_lupus_familiaris",
    ]
    got = sn.match_species(
        [
            "gallus_gallus",
            "neovison_vison",
            "aedes_aegypti_lvpagwg",
            "pyrococcus_horikoshii_ot3_gca_000011105",
        ],
        tips,
    )
    assert got == {
        "gallus_gallus": "gallus_gallus_reference_breed",
        "neovison_vison": "neogale_vison",
        "aedes_aegypti_lvpagwg": "aedes_aegypti_strain_lvp_agwg",
        "pyrococcus_horikoshii_ot3_gca_000011105": "pyrococcus_horikoshii_ot3",
    }
    with pytest.raises(ValueError, match="exactly one"):
        sn.match_species(["canis_lupus"], tips)  # two candidates is an error, never a choice
    with pytest.raises(ValueError, match="two species on one tip"):
        sn.match_species(["gallus_gallus", "gallus_gallus_reference_breed"], tips)


def test_pruning_collapses_the_nodes_left_with_one_child() -> None:
    root = sn.parse_newick("(((A,B)ab,C)abc,(D,E)de)r;")
    out = sn.prune(root, {"a", "c", "d"})
    assert sn.to_newick(out) == "((a,c),d);"


def small_tree():
    root = sn.parse_newick("((a:1,b:2)ab:3,c:4)r;")
    return root, sn.arrays(root, ["c", "a", "b"])


def test_covariance_is_the_shared_path_from_the_root() -> None:
    root, t = small_tree()
    cov = sn.covariance(t, sn.newick_lengths(root))
    # columns c, a, b: c alone at 4; a and b share the 3 above ab
    assert cov.tolist() == [[4.0, 0.0, 0.0], [0.0, 4.0, 3.0], [0.0, 3.0, 5.0]]


def test_grafen_lengths_are_ultrametric_with_the_root_at_one() -> None:
    root = sn.parse_newick("(((a,b),c,d),((e,f),g));")
    t = sn.arrays(root, list("abcdefg"))
    cov = sn.covariance(t, sn.grafen_lengths(t))
    assert np.allclose(np.diag(cov), 1.0)
    assert cov[0, 1] > cov[0, 2] > cov[0, 4] == 0.0
    assert np.all(np.linalg.eigvalsh(cov) > 0)


def test_the_brownian_bank_has_the_tree_s_covariance() -> None:
    root = sn.parse_newick("(((a,b),c,d),((e,f),g));")
    t = sn.arrays(root, list("abcdefg"))
    cov = sn.covariance(t, sn.grafen_lengths(t))
    bm = sn.Brownian(cov, sims=40_000, seed=1)
    # rows are centred, so compare with the centred covariance
    h = np.eye(7) - 1 / 7
    assert np.allclose(np.cov(bm.w1, rowvar=False), h @ cov @ h, atol=0.03)
    assert abs(float(np.corrcoef(bm.w1[:, 0], bm.w2[:, 0])[0, 1])) < 0.03  # the banks are independent


def test_the_one_pass_pearson_is_the_row_by_row_pearson() -> None:
    root = sn.parse_newick("(((a,b),c,d),((e,f),g));")
    t = sn.arrays(root, list("abcdefg"))
    bm = sn.Brownian(sn.covariance(t, sn.grafen_lengths(t)), sims=50, seed=2)
    rng = np.random.default_rng(3)
    fa, fb = rng.normal(size=7), rng.normal(size=7)
    fast = bm.null_r(fa, fb, 0.7, 1.3)
    slow = [sn.pearson(fa + 0.7 * bm.w1[i], fb + 1.3 * bm.w2[i]) for i in range(50)]
    assert np.allclose(fast, slow)


def clade_tree(n: int, stem: float, tip: float) -> tuple[sn.TreeArrays, np.ndarray]:
    a = ",".join(f"a{i}:{tip}" for i in range(n))
    b = ",".join(f"b{i}:{tip}" for i in range(n))
    root = sn.parse_newick(f"(({a}):{stem},({b}):{stem});")
    t = sn.arrays(root, [f"a{i}" for i in range(n)] + [f"b{i}" for i in range(n)])
    return t, sn.covariance(t, sn.newick_lengths(root))


def test_shared_ancestry_alone_makes_two_clade_patterns_correlate() -> None:
    """The reason for the null: one deep split looks like forty agreeing species on a star."""
    n = 20
    _, cov_deep = clade_tree(n, stem=0.99, tip=0.01)
    rng = np.random.default_rng(4)
    split = np.r_[np.ones(n), -np.ones(n)]
    da = split + rng.normal(scale=0.3, size=2 * n)
    db = split + rng.normal(scale=0.3, size=2 * n)
    r_obs = sn.pearson(da, db)
    assert r_obs > 0.8
    zero = np.zeros(2 * n)
    deep = sn.Brownian(cov_deep, sims=20_000, seed=5)
    star = sn.Brownian(np.eye(2 * n), sims=20_000, seed=5)
    p_deep = float((deep.null_r(zero, zero, 1.0, 1.0) >= r_obs).mean())
    p_star = float((star.null_r(zero, zero, 1.0, 1.0) >= r_obs).mean())
    assert p_star < 0.001  # forty independent species: a striking agreement
    assert p_deep > 0.2  # two clades: one coincidence of sign, close to a coin toss


def test_bh_is_the_procedure_section_14_used() -> None:
    rng = random.Random(6)
    ps = [rng.random() ** 3 for _ in range(200)]
    q = sn.bh(ps)
    order = sorted(range(len(ps)), key=lambda i: ps[i])
    prev = 1.0
    want = [0.0] * len(ps)
    for k in range(len(ps) - 1, -1, -1):
        prev = min(prev, ps[order[k]] * len(ps) / (k + 1))
        want[order[k]] = prev
    assert q == want


def test_presence_bits_read_the_hex_as_profiling_does() -> None:
    rng = random.Random(7)
    n = 13
    width = (n + 3) // 4
    genes = {f"g{i}": format(rng.getrandbits(n), f"0{width}x") for i in range(30)}
    _, names, bits = sn.presence_bits({"species": [f"s{i}" for i in range(n)], "genes": genes})
    for i, g in enumerate(names):
        assert bits[i].astype(int).tolist() == pf.unpack(genes[g], n)


def synthetic_pairs() -> sn.PairInputs:
    rng = np.random.default_rng(8)
    n_sp, n_genes = 12, 300
    bits = rng.random((n_genes, n_sp)) < np.linspace(0.2, 0.9, n_sp)
    genes = [f"g{i}" for i in range(n_genes)]
    hexes = {
        g: format(int("".join("1" if b else "0" for b in row), 2), "03x")
        for g, row in zip(genes, bits, strict=True)
    }
    presence = {"species": [f"s{i}" for i in range(n_sp)], "genes": hexes}
    origin = {"genes": {g: {"origin": ("Life", "Metazoa", "Mammalia")[i % 3]} for i, g in enumerate(genes)}}
    return sn.pair_inputs(presence, origin)


def test_the_age_expectation_is_the_mean_of_section_14_s_matched_draws() -> None:
    d = synthetic_pairs()
    rows = np.array([0, 1, 2, 3, 6, 9, 12, 15], dtype=np.int64)  # mixed strata
    rng = np.random.default_rng(9)
    draws = np.mean([sn.profile(d, sn.matched_rows(d, rows, rng)) for _ in range(4000)], axis=0)
    assert np.allclose(draws, sn.age_expectation(d, rows), atol=0.01)
    one = sn.matched_rows(d, rows, rng)
    assert sorted(d.origin_of_row[one].tolist()) == sorted(d.origin_of_row[rows].tolist())


def test_testable_pairs_are_those_with_enough_members_on_each_side() -> None:
    sets = {"x": set(range(20)), "y": set(range(15, 40)), "z": set(range(3, 18))}
    got = {(a, b): (len(oa), len(ob)) for a, b, oa, ob in sn.testable_pairs(sets)}
    # x - z keeps 0..2 and 18,19: five members, under the floor of eight
    assert got == {("x", "y"): (15, 20), ("y", "z"): (22, 12)}


# ---------------------------------------------------------------------------------------- Dollo


def dollo_toy():
    # columns: a, b, c, d (species) and human last
    root = sn.parse_newick("(((homo_sapiens,a)ha,b)hab,(c,d)cd)r;")
    t = sn.arrays(root, ["a", "b", "c", "d", "homo_sapiens"])
    bits = np.array(
        [
            [1, 0, 0, 0],  # born at ha
            [0, 1, 0, 0],  # born at hab, lost on a
            [0, 0, 1, 0],  # born at the root, lost on a, b and d
            [1, 1, 1, 1],  # born at the root, never lost
            [0, 0, 0, 0],  # human only
            [1, 0, 1, 1],  # born at the root, lost on b
        ],
        dtype=bool,
    )
    return t, bits, sn.Dollo(t, bits, human_column=4)


def test_birth_is_the_common_ancestor_of_human_and_the_carriers() -> None:
    t, _, dl = dollo_toy()
    names = [t.names[int(dl.path[j])] for j in dl.birth]
    assert names == ["ha", "hab", "r", "r", "homo_sapiens", "r"]


def test_losses_are_counted_where_the_gene_leaves_a_clade() -> None:
    t, _, dl = dollo_toy()
    idx = {n: i for i, n in enumerate(t.names)}
    # human carries every gene, so ha and hab hold every gene born at or above them: gene 1 is lost
    # on a; gene 2 on a, b and d; gene 5 on b
    assert dl.losses[idx["a"]] == 2
    assert dl.losses[idx["b"]] == 2
    assert dl.losses[idx["d"]] == 1
    assert dl.losses[idx["hab"]] == 0
    assert dl.losses[idx["cd"]] == 0
    assert all(dl.p[dl.on_path] == 0)  # nothing is ever lost on the way to human


def test_the_fitted_rate_gives_each_gene_its_own_number_of_species() -> None:
    _, bits, dl = dollo_toy()
    for row in range(len(bits)):
        assert dl.expected_carriers(row) == pytest.approx(bits[row].sum(), abs=1e-6)


def test_a_simulated_gene_stays_in_its_birth_clade_and_keeps_its_count_on_average() -> None:
    _, _, dl = dollo_toy()
    sims = dl.simulate(4000, seed=10)
    assert not sims[0][:, 1:].any()  # born at ha: never in b, c or d
    assert not sims[1][:, 2:].any()  # born at hab: never in c or d
    assert sims[3].all()  # never lost
    assert not sims[4].any()  # human only
    assert sims[5].sum(axis=1).mean() == pytest.approx(3, abs=0.08)


def test_mean_pair_jaccard_is_the_mean_of_section_19_s_pairs() -> None:
    rng = np.random.default_rng(11)
    x = rng.random((25, 30)) < 0.4
    x[3] = False  # two empty rows: their union is empty and the pair scores 0
    x[4] = False
    masks = [int("".join("1" if b else "0" for b in row), 2) for row in x]
    want = np.mean([pg.jaccard(masks[i], masks[j]) for i in range(25) for j in range(i + 1, 25)])
    assert sn.mean_pair_jaccard(x) == pytest.approx(want, abs=1e-6)


def test_birth_matched_sets_keep_each_birth_node_s_count() -> None:
    _, _, dl = dollo_toy()
    rng = np.random.default_rng(12)
    rows = np.array([2, 5], dtype=np.int64)  # both born at the root
    for _ in range(20):
        got = sn.birth_matched(dl, rows, rng)
        assert sorted(dl.birth[got].tolist()) == sorted(dl.birth[rows].tolist())
        assert len(set(got.tolist())) == len(got)


def test_classes_reduce_to_section_19_s_on_twenty_draws() -> None:
    assert [sn.classify(k, 20) for k in (0, 1, 9, 10, 20)] == [
        "clear",
        "marginal",
        "marginal",
        "at_or_below",
        "at_or_below",
    ]
    assert [sn.classify(k, 200) for k in (9, 10, 99, 100)] == ["clear", "marginal", "marginal", "at_or_below"]
    assert sn.chance_ceiling(42) == 5  # binomial(42, 0.05): P(<= 4) < 0.95 <= P(<= 5)


# -------------------------------------------------------------------------------------- outcomes


def placeholders(text: str) -> set[str]:
    return set(re.findall(r"{(\w+)}", text))


def stored_14(n_pairs: int = 4, called: int = 2) -> dict:
    pairs = [{"libraries": [f"l{i}", f"m{i}"], "r_exclusive": 0.1 * i} for i in range(n_pairs)]
    return {"pairs": pairs, "significant_q05": pairs[:called]}


def null_a_rows(qs: list[float], cal_share: float = 0.04, halves_share: float = 0.8) -> dict:
    return {
        "pairs": [
            {"libraries": [f"l{i}", f"m{i}"], "r_exclusive": 0.1 * i, "p_tree": q, "q_tree": q}
            for i, q in enumerate(qs)
        ],
        "calibration": {
            "pairs": 100,
            "at_or_below_005": round(cal_share * 100),
            "share_at_or_below_005": cal_share,
        },
        "halves": {
            "libraries": 10,
            "at_or_below_005": round(halves_share * 10),
            "share_at_or_below_005": halves_share,
        },
    }


def test_null_a_s_rule_is_applied_in_its_registered_order() -> None:
    s = stored_14()
    assert sn.select_a(null_a_rows([0.5, 0.5, 0.5, 0.5]), None, s)[0] == "none_survive"
    assert sn.select_a(null_a_rows([0.5] * 4, halves_share=0.4), None, s)[0] == "none_survive_low_power"
    key, text, f = sn.select_a(null_a_rows([0.01, 0.5, 0.02, 0.5]), null_a_rows([0.5] * 4), s)
    assert key == "some_survive"
    assert (f["n_surv"], f["k_14"], f["k_new"], f["lost_14"]) == (2, 1, 1, 1)
    assert "Under equal branch lengths" in text and "{" not in text
    # an uncalibrated null decides before any survivor is read
    assert sn.select_a(null_a_rows([0.01] * 4, cal_share=0.2), None, s)[0] == "null_not_calibrated"
    bad = null_a_rows([0.5] * 4)
    bad["pairs"][0]["r_exclusive"] = 0.3
    assert sn.select_a(bad, None, s)[0] == "not_reproduced"


def b_rows(classes: list[str]) -> list[dict]:
    return [
        {"library": f"l{i}", "members_with_a_profile": 10, "mean_pair_jaccard": 0.5, "class": c}
        for i, c in enumerate(classes)
    ]


def stored_19(n: int, clear: int) -> dict:
    return {
        "libraries": [
            {
                "library": f"l{i}",
                "members_with_a_profile": 10,
                "mean_pair_jaccard": 0.5,
                "matched_draws_at_or_above_observed": 0 if i < clear else 5,
            }
            for i in range(n)
        ]
    }


def test_null_b_s_rule_is_applied_in_its_registered_order() -> None:
    s = stored_19(42, 21)
    assert sn.select_b(b_rows(["marginal"] * 42), s)[0] == "none_clear"
    assert sn.select_b(b_rows(["clear"] * 5 + ["marginal"] * 37), s)[0] == "at_chance"
    key, text, f = sn.select_b(b_rows(["clear"] * 6 + ["at_or_below"] * 36), s)
    assert key == "some_clear" and "{" not in text
    assert (f["n_clear"], f["k_19"], f["n_below"]) == (6, 6, 36)
    off = b_rows(["clear"] * 42)
    off[0]["mean_pair_jaccard"] = 0.51
    assert sn.select_b(off, s)[0] == "not_reproduced"


def test_every_wording_states_only_figures_the_rule_fills() -> None:
    a = {
        "n_tested", "n_tested_14", "n_mismatch", "n_14", "n_surv", "k_14", "k_new", "lost_14",
        "cal_hits", "n_cal", "cal_share", "halves_hits", "n_halves",
    }  # fmt: skip
    for text in sn.A_OUTCOMES.values():
        assert placeholders(text) <= a
    assert placeholders(sn.A_SENSITIVITY) <= a | {"n_surv_eq", "k_14_eq", "cal_hits_eq"}
    b = {"n_lib", "n_mismatch", "n_sets", "n_clear", "n_marg", "n_below", "n_19", "k_19", "k_new"}
    for text in sn.B_OUTCOMES.values():
        assert placeholders(text) <= b | {"expected_chance", "ceiling"}


def test_no_wording_calls_a_pair_a_dependency_without_saying_it_is_not() -> None:
    allowed = re.compile(
        r"no library pair is called a dependency|calls nothing a dependency|nothing is called a "
        r"dependency|the condition the area set before a pair could be called a dependency",
        re.IGNORECASE,
    )
    for text in (*sn.A_OUTCOMES.values(), *sn.B_OUTCOMES.values()):
        if "dependency" in text:
            assert allowed.search(text), text


# ------------------------------------------------------------------------- on the project's tree

PAN_FILE = "data/knowledge/species_tree/pan_default_protein-trees_default.nh"
VERT_FILE = "data/knowledge/species_tree/vert_default_protein-trees_default.nh"
TREES_HOW = "run scripts/species_null_register.py --fetch (two Ensembl species trees, about 13 KB)"


@pytest.mark.needs_local_data(PAN_FILE, VERT_FILE, how=TREES_HOW)
def test_the_project_s_species_sit_where_the_taxonomy_puts_them() -> None:
    from genomeos.results import load_result

    species = load_result("origin_presence_genome_wide")["species"]
    root, _ = sn.species_tree(species, with_human=True)
    tips = {t.name: t for t in sn.leaves(root)}
    assert len(tips) == len(species) + 1 == 356
    expect = {
        "drosophila_melanogaster": "Bilateria",
        "saccharomyces_cerevisiae": "Opisthokonta",
        "arabidopsis_thaliana": "Eukaryota",
        "danio_rerio": "Euteleostomi",
        "gallus_gallus": "Amniota",
        "pan_troglodytes": "Homininae",
    }
    for other, clade in expect.items():
        assert sn.mrca([tips["homo_sapiens"], tips[other]]).name == clade
    assert sn.sha256(Path(PAN_FILE)) == sn.TREE_SHA256[sn.PAN_TREE.name]
    assert sn.sha256(Path(VERT_FILE)) == sn.TREE_SHA256[sn.VERT_TREE.name]


# ------------------------------------------------------------------- the two nulls end to end


TREE12 = "(((s0,s1,s2),(s3,(s4,s5))),((s6,s7),(s8,(s9,(s10,s11)))));"


def test_null_a_runs_every_testable_pair_with_its_calibration_and_halves() -> None:
    d = synthetic_pairs()
    t = sn.arrays(sn.parse_newick(TREE12), d.species)
    bm = sn.Brownian(sn.covariance(t, sn.grafen_lengths(t)), sims=500, seed=13)
    genes = d.genes
    sets = {"x": set(genes[0:40]), "y": set(genes[30:80]), "z": set(genes[100:140])}
    out = sn.null_a(d, bm, sets)
    assert out["tested"] == 3
    for row in out["pairs"]:
        assert 1 / 501 <= row["p_tree"] <= 1 and row["q_tree"] >= row["p_tree"]
        assert row["r_exclusive"] == pytest.approx(
            pf.pearson(
                list(
                    sn.profile(d, sn.rows_of(d, sets[row["libraries"][0]] - sets[row["libraries"][1]]))
                    - d.genome
                ),
                list(
                    sn.profile(d, sn.rows_of(d, sets[row["libraries"][1]] - sets[row["libraries"][0]]))
                    - d.genome
                ),
            ),
            abs=1e-4,
        )
    assert out["calibration"]["pairs"] == 3
    assert out["halves"]["libraries"] == 3


def test_null_b_reports_every_library_against_its_random_sets() -> None:
    rng = np.random.default_rng(14)
    species = [f"s{i}" for i in range(12)]
    newick = TREE12.replace("(s10,s11)", "(s10,(s11,homo_sapiens))")
    t = sn.arrays(sn.parse_newick(newick), [*species, "homo_sapiens"])
    bits = rng.random((200, 12)) < 0.5
    bits[:30, 6:] = True  # thirty genes that agree on half the tree
    bits[:30, :6] = False
    dl = sn.Dollo(t, bits, human_column=12)
    sims = dl.simulate(10, seed=15)
    rows = sn.null_b(
        dl, bits, sims, {"together": np.arange(30), "random": np.arange(100, 140)}, sets=20, seed=16
    )
    by = {r["library"]: r for r in rows}
    assert set(by) == {"together", "random"}
    assert all(r["class"] in ("clear", "marginal", "at_or_below") for r in rows)
    assert by["together"]["mean_pair_jaccard"] == pytest.approx(1.0)
