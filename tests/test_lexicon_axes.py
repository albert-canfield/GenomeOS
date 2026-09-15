# SPDX-License-Identifier: AGPL-3.0-or-later
"""The lexicon's base-resolution axes: the phyloP byte scale, the trinucleotide classes, the indexed
motif scan against the reference scanner, the decoys, and the tests' arithmetic."""

import random

import numpy as np
import pytest

from genomeos.attribution import lexicon_axes as A  # noqa: N812 - the module is read as a namespace


def test_phylop_bytes_put_the_threshold_on_a_bin_edge():
    p = np.array([2.27, 2.2699, 4.0, -3.0, np.nan, 50.0, -50.0], dtype=np.float32)
    q = A.quantise_phylop(p)
    assert q[0] >= A.PHYLOP_ANCHOR and q[1] < A.PHYLOP_ANCHOR
    assert q[4] == A.PHYLOP_MISSING and q[5] == 127 and q[6] == -127
    back = A.dequantise_phylop(q)
    assert abs(back[2] - 4.0) <= 0.1 and np.isnan(back[4])
    bins = A.phylop_bin(q)
    assert bins[4] == 0 and bins[0] > bins[1] and bins[2] == 6


def test_trinucleotide_classes_collapse_strands():
    fwd = A.trinucleotide_class(A.encode("ACGTN"))
    rc = A.trinucleotide_class(A.encode("NACGT"))  # reverse complement of ACGTN
    # ACG (centre C) and its reverse complement CGT (centre G) are one class
    assert fwd[1] == rc[3]
    assert fwd[0] == 255 and fwd[3] == 255  # a neighbour off the end or an N
    assert set(A.trinucleotide_class(A.encode("ACGTACGTTTGCA"))[1:-1]) <= set(range(32))


def test_indexed_scan_matches_the_reference_scanner():
    from genomeos.genome.motifs import all_hits, load_motifs

    rng = random.Random(3)
    seq = (
        "".join(rng.choice("ACGT") for _ in range(1500))
        + "N" * 20
        + "".join(rng.choice("ACGT") for _ in range(400))
    )
    motifs = load_motifs(relative=0.95)[:80]
    codes = A.encode(seq)
    order, offsets = A.kmer_order(codes)
    mine = set()
    for m in motifs:
        for strand_motif, strand in ((m, 0), (A.reverse_complement_motif(m), 1)):
            for s in A.scan_profile(strand_motif, codes, order, offsets):
                mine.add((m.name.upper(), int(s), strand))
    ref = {(f, st, sd) for f, st, _, sd, _ in all_hits(motifs, seq)}
    assert mine == ref


def test_decoys_keep_composition_and_dedupe_merges_identical_site_sets():
    from genomeos.genome.motifs import load_motifs

    m = load_motifs(relative=0.95)[5]
    d = A.decoy_motifs([m])[0]
    assert sorted(map(tuple, zip(*m.counts, strict=True))) == sorted(map(tuple, zip(*d.counts, strict=True)))
    kept, aliases = A.dedupe_sites({"B": [(1, 5), (9, 12)], "A": [(1, 5), (9, 12)], "C": [(2, 6)]})
    assert set(kept) == {"A", "C"} and aliases["A"] == ["B"]


def fake_arrays(length=10_000, variable_every=100):
    q = np.full(length, 0, dtype=np.int8)
    q[2000:3000] = 40  # a constrained stretch
    informative = np.ones(length, dtype=bool)
    variable = np.zeros(length, dtype=bool)
    variable[::variable_every] = True
    ctx = np.zeros(length, dtype=np.uint8)
    return {
        "length": length,
        "ctx": ctx,
        "q": q,
        "measured": q != A.PHYLOP_MISSING,
        "constrained": q >= A.PHYLOP_ANCHOR,
        "cell_m": np.zeros(length, dtype=np.uint16),
        "rate_m": np.array([0.1, 0.0], dtype=np.float32),
        "cell_ms": np.zeros(length, dtype=np.uint8),
        "rate_ms": np.array([0.1, 0.0], dtype=np.float32),
        "gnocchi": None,
        "human": True,
        "informative": informative,
        "variable": variable,
        "cell0": np.zeros(length, dtype=np.uint16),
        "rate0": np.array([1 / variable_every, 0.0], dtype=np.float32),
        "cell1": np.zeros(length, dtype=np.uint16),
        "rate1": np.array([1 / variable_every, 0.0], dtype=np.float32),
        "cell_t": np.zeros(length, dtype=np.uint16),
        "rate_t": np.array([1 / variable_every, 0.0], dtype=np.float32),
    }


def test_occurrence_sums_and_rows():
    arr = fake_arrays()
    starts = np.array([2000, 2500, 5000, -5, 9_990])
    ends = np.array([2100, 2600, 5100, 10, 10_100])
    s = A.occurrence_sums(starts, ends, arr)
    assert list(s["bases"]) == [100, 100, 100, 10, 10]
    assert list(s["m_obs"][:3]) == [100, 100, 0]
    assert s["m_exp"][0] == pytest.approx(10.0)
    assert s["h_obs"][0] == 1 and s["h_e0"][0] == pytest.approx(1.0)
    positions = [(2000 + 40 * i, 2020 + 40 * i) for i in range(25)]
    rows = A.unit_table("u", "test", positions, arr)
    star = rows[0]
    assert star["constrained_fraction"] == 1.0 and star["z_mammals"] == pytest.approx(5.0)
    assert star["constrained_flank_fraction"] > 0.5  # the flanks sit in the same constrained stretch


def test_sandwich_z_floors_at_poisson():
    obs = np.array([1.0, 0.0, 1.0, 0.0])
    exp = np.array([0.5, 0.5, 0.5, 0.5])
    assert A.sandwich_z(obs, exp) == pytest.approx(0.0)
    assert A.sandwich_z(np.array([5.0]), np.array([1.0])) == pytest.approx(4 / 4.0)
    assert A.sandwich_z(obs, np.zeros(4)) is None


def test_classify_applies_fixed_thresholds_to_a_control_family():
    def row(pm, ph, z1=-1.0):
        return {
            "context": "intron",
            "p_mammals": pm,
            "p_human_depleted": ph,
            "p_human_given_phylop": 0.5,
            "z_human_given_phylop": z1,
            "p_mammals_flank": pm,
            "p_human_flank_depleted": ph,
            "p_human_flank_two_sided": 0.5,
            "z_human_flank": -1.0,
        }

    real = A.classify([row(1e-9, 1e-9), row(0.5, 0.5), row(1e-9, 0.9)])
    assert real["counts"]["syntax"] == 1 and real["counts"]["held_mammals"] == 2
    control = A.classify([row(0.01, 0.01), row(0.4, 0.4)], fixed=real["tests"])
    assert control["counts"]["syntax"] == 0


def test_sign_test():
    assert A._sign_p(5, 10) == pytest.approx(1.0)
    assert A._sign_p(0, 10) < 0.01
    assert A._sign_p(0, 0) is None


def test_motif_against_its_decoy_and_the_conditional_share():
    def row(kind, factor, site_oe, flank_oe, cons, flank_cons, held_m=True, held_p=False):
        return {
            "kind": kind,
            "unit": f"{kind}:{factor}",
            "context": "intron",
            "human_oe": site_oe,
            "human_oe_flank": flank_oe,
            "human_expected_trinucleotide": 50.0,
            "constrained_fraction": cons,
            "constrained_flank_fraction": flank_cons,
            "held_mammals_vs_flanks": held_m,
            "held_people_vs_flanks": held_p,
        }

    rows = []
    for i in range(12):
        rows.append(row("jaspar95", f"F{i}", 0.8, 1.0, 0.04, 0.02, held_p=i < 6))
        rows.append(row("decoy95", f"F{i}", 1.0, 1.0, 0.02, 0.02, held_p=i < 1))
    out = A.motif_against_decoy(rows)["intron"]
    assert out["people_motif_less_variable"] == 12 and out["people_sign_p"] < 0.001
    assert out["mammal_median_difference"] == pytest.approx(1.0)
    cond = A.conditional_people(rows)
    assert cond["jaspar95"]["share"] == 0.5 and cond["decoy95"]["share"] == pytest.approx(1 / 12, abs=1e-3)
