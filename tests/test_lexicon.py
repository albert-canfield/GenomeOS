# SPDX-License-Identifier: AGPL-3.0-or-later
"""The lexicon: its statistics, its per-base maps and strata, the vocabulary found from sequence, and
the two questions' nulls on synthetic chromosomes where the answer is known."""

import math
import random

import pytest

from genomeos.attribution import lexicon as L  # noqa: N812 - the module is read as a namespace


def poisson_sf(observed, expected):
    return 1.0 - sum(math.exp(-expected) * expected**k / math.factorial(k) for k in range(observed))


@pytest.mark.parametrize("observed,expected", [(1, 0.5), (5, 0.47), (10, 10.0), (12, 2.35), (3, 20.0)])
def test_poisson_tail_matches_the_direct_sum(observed, expected):
    assert L.poisson_tail(observed, expected) == pytest.approx(
        poisson_sf(observed, expected), rel=1e-6, abs=1e-12
    )


def test_poisson_tail_does_not_underflow_on_large_expectations():
    # exp(-1000) underflows to zero; a naive lower-tail sum then calls every count unremarkable
    assert L.poisson_tail(1200, 1000.0) < 1e-8
    assert L.poisson_tail(900, 1000.0) > 0.99
    assert L.poisson_tail(0, 5.0) == 1.0
    assert L.poisson_tail(3, 0.0) == 0.0


def test_benjamini_hochberg_threshold_and_count():
    thr, passing = L.benjamini_hochberg([0.001, 0.008, 0.039, 0.041, 0.9], fdr=0.05)
    assert passing == 2 and thr == 0.008
    assert L.benjamini_hochberg([]) == (0.0, 0)


def test_cosine():
    assert L.cosine({"a": 1, "b": 0}, {"a": 2}) == pytest.approx(1.0)
    assert L.cosine({"a": 1}, {"b": 1}) == 0.0


def fake_maps(monkeypatch, length, paints=(), repeats=(), conserved=(), measured=(), human=()):
    monkeypatch.setattr(L.Maps, "_fill", lambda self, results_dir: None)
    m = L.Maps("chrT", length)
    for a, b, ctx in paints:
        m._paint(a, b, ctx)
    for arr, spans in (
        (m.repeat, repeats),
        (m.conserved, conserved),
        (m.human_measured, measured),
        (m.human, human),
    ):
        for a, b in spans:
            m._flag(arr, a, b)
    return m


def test_windows_strata_and_context_matched_conservation(monkeypatch):
    rng = random.Random(1)
    length = 40_000
    # four windows: two GC-poor, two GC-rich; the second half is repeat-covered
    seq = "".join(
        "".join(
            rng.choice("AT" if w < 2 else "GC") if rng.random() < 0.8 else rng.choice("ACGT")
            for _ in range(10_000)
        )
        for w in range(4)
    )
    m = fake_maps(
        monkeypatch,
        length,
        paints=[(0, 1_000, "cds"), (20_000, 21_000, "cds")],
        repeats=[(20_000, 40_000)],
        conserved=[(0, 1_000), (5_000, 5_100)],
    )
    w = L.Windows(seq, m)
    assert len(w.rows) == 4
    assert w.rows[0]["context_bp"]["cds"] == 1_000
    assert w.rows[0]["context_conserved"] == {"cds": 1_000, "other": 100}
    st = w.rows[0]["stratum"]
    assert w.stratum_at(15) == st and w.stratum_at(10**9) is None
    assert w.conserved_share_in_context[(st, "cds")] > w.conserved_share_in_context[(st, "other")]


def test_context_enrichment_is_flat_when_the_unit_ignores_context(monkeypatch):
    m = fake_maps(monkeypatch, 20_000, paints=[(0, 5_000, "intron")])
    w = L.Windows("ACGT" * 5_000, m)
    rng = random.Random(3)
    occ = []
    for _ in range(400):
        p = rng.randrange(20_000)
        occ.append((p, m.context_at(p), w.stratum_at(p)))
    out = L.context_enrichment(occ, w)
    assert out["intron"]["expected"] == pytest.approx(100, rel=0.01)
    assert out["intron"]["p"] > 0.001
    # a unit that sits only in introns is enriched there
    only = [(p, "intron", w.stratum_at(p)) for p in range(0, 5_000, 50)]
    assert L.context_enrichment(only, w)["intron"]["ratio"] > 3


def test_markov_expectation_reproduces_a_markov_sequence():
    # a sequence with no word structure beyond its trinucleotides: 6-mers fall near expectation
    rng = random.Random(5)
    seq = "".join(rng.choice("ACGT") for _ in range(200_000))
    counts = {}
    for i in range(len(seq) - 5):
        counts[seq[i : i + 6]] = counts.get(seq[i : i + 6], 0) + 1
    from collections import Counter

    di, tri = Counter(), Counter()
    for word, n in counts.items():
        for x in range(5):
            di[word[x : x + 2]] += n
        for x in range(4):
            tri[word[x : x + 3]] += n
    total = sum(counts.values())
    exp = L._markov2_expectation("ACGTAC", di, tri, total)
    assert exp == pytest.approx(total / 4096, rel=0.1)


def test_seeds_find_a_planted_repeat_and_extension_measures_it(monkeypatch):
    monkeypatch.setattr(L, "SEED_OCCURRENCE_CAP", 30)
    rng = random.Random(7)
    motif = "".join(rng.choice("ACGT") for _ in range(40))
    parts = []
    for _ in range(80):
        parts.append("".join(rng.choice("ACGT") for _ in range(300)))
        parts.append(motif)
    seq = "".join(parts)
    found = L.seeds(seq, k=16, minimum=50, top=100)
    assert found, "the planted 40-base segment gives 25 seeds of 16 occurring 80 times"
    word, pos = next(iter(found.items()))
    assert -pos[-1] == 80 and len(pos) - 1 == 30  # true count kept, occurrences sampled to the cap
    assert word in motif
    ext = L.extend(seq, pos[:-1], k=16)
    assert ext["consensus_length"] == 40
    families = L.seed_families({w: p[:-1] for w, p in found.items()}, k=16)
    assert len(families) == 1 and families[0]["span"] == 40


def test_fossil_similarity_sees_proximity_not_the_partition():
    # composition drifts smoothly along the chromosome and ignores the node boundaries
    rng = random.Random(11)
    blocks = []
    for i in range(80):
        x = i / 80
        comp = {"L1:A": 1000 * (1 - x) + rng.random() * 50, "L2:B": 1000 * x + rng.random() * 50}
        blocks.append(
            {"start": i * 10_000, "end": i * 10_000 + 2_000, "composition": comp, "repeat_bp": 1000 + i}
        )
    nodes = [{"id": f"n{k}", "start": k * 50_000, "end": (k + 1) * 50_000} for k in range(16)]
    out = L.fossil_node_similarity(blocks, nodes)
    assert out["within_node_mean_cosine"] > out["between_node_mean_cosine"]
    assert out["p_label_permutation"] < 0.05  # neighbours are alike
    assert out["p_boundary_shift"] > 0.05  # but moving the boundaries keeps it: no partition effect


def test_syntax_classes_and_the_own_element_null():
    def unit(name, kind, cons, exp, human, cp, **kw):
        return {
            "unit": name,
            "kind": kind,
            "occurrences": 50,
            "bp": 400,
            "contexts": {"enhancer_like": 50},
            "conserved_fraction": cons,
            "conserved_expected_fraction": exp,
            "human_constrained_fraction": human,
            "conserved_p": cp,
            "human_p": 0.5,
            "block_tiers": {"regulatory": 5},
            "block_cases": {"syntax": 5},
            **kw,
        }

    index = {
        "units": [
            unit(
                "jaspar:HELD",
                "jaspar",
                0.5,
                0.05,
                0.5,
                1e-9,
                host_conserved_p=1e-9,
                host_conserved_expected_fraction=0.1,
            ),
            unit(
                "jaspar:BORROWED",
                "jaspar",
                0.5,
                0.05,
                0.5,
                1e-9,
                host_conserved_p=0.6,
                host_conserved_expected_fraction=0.5,
            ),
            unit("repeat:X:Y", "repeat", 0.02, 0.03, 0.05, 0.9),
            unit("seed:ACGT", "seed", 0.02, 0.03, None, 0.9),
        ]
    }
    out = L.syntax_classes(index)
    assert out["by_case"] == {"syntax": 2, "tolerant": 1, "unmeasured": 1}
    assert out["syntax_survival"]["jaspar_lost_to_own_element_null"] == 1
    assert [u["unit"] for u in out["syntax_candidates"]] == ["jaspar:HELD"]


def test_locus_sampling_is_uniform_and_shared_across_seeds():
    a = list(range(0, 1_000_000, 500))  # 2,000 occurrences of one seed
    b = [p + 20 for p in a]  # its neighbour in the same segment
    sa, sb = L.sample_loci(a, cap=200), L.sample_loci(b, cap=200)
    assert len(sa) == 200 and sa == sorted(sa)
    shared = len({p + 20 for p in sa} & set(sb))
    assert shared > 100  # the same loci are kept for both, most of the time
    assert min(sa) < 250_000 and max(sa) > 750_000  # spread along the chromosome, not the first ones


def test_tandem_words_and_self_defined_contexts():
    assert L.is_tandem("CACACACACACACACA") and L.is_tandem("GATGGATGGATGGATG")
    assert not L.is_tandem("GACAGAAGCATTCTCA")
    assert "promoter_like" in L.SELF_DEFINED["ccre:PLS"]
