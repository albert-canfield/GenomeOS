"""The compression probe: every code length checked against a brute-force count on small sequences."""

import math
import random

import numpy as np
import pytest

from genomeos.attribution import compress as cz


def _seq(n, seed=1, alphabet="ACGT"):
    rng = random.Random(seed)
    return "".join(rng.choice(alphabet) for _ in range(n))


def _codes(text):
    return cz.Chromosome.from_text("t", text).s


def _ctx(s, i, k):
    return tuple(int(s[i - j]) if i - j >= 0 else 0 for j in range(1, k + 1))


def test_contexts_match_a_loop():
    s = _codes(_seq(50))
    for k in (0, 1, 3, 5):
        x = cz.contexts(s, k)
        for i in range(len(s)):
            want = sum(v << (2 * j) for j, v in enumerate(_ctx(s, i, k)))
            assert int(x[i]) == want
    rc = cz.rc_contexts(s, 3)
    for j in range(len(s) - 3):
        want = sum((3 - int(s[j + m])) << (2 * (m - 1)) for m in range(1, 4))
        assert int(rc[j]) == want


@pytest.mark.parametrize("revcomp", [False, True])
def test_adaptive_counts_equal_sequential_counting(revcomp, monkeypatch):
    s = _codes(_seq(400, seed=3, alphabet="AACGTTT"))
    k = 3
    for parts in (1 << 25, 64):  # one part, and the hash space split into many
        monkeypatch.setattr(cz, "COUNT_GROUP", parts)
        n_ca, n_c = cz.adaptive_counts(s, k, revcomp=revcomp)
        seen: dict = {}
        n = len(s)
        for i in range(n):
            c = _ctx(s, i, k)
            assert n_ca[i] == seen.get((c, int(s[i])), 0)
            assert n_c[i] == sum(seen.get((c, a), 0) for a in range(4))
            seen[(c, int(s[i]))] = seen.get((c, int(s[i])), 0) + 1
            j = i - k  # the reverse-strand event at j is known once base j + k is
            if revcomp and j >= 0 and j + k < n - 1:
                rc = tuple(3 - int(s[j + m]) for m in range(1, k + 1))
                key = (rc, 3 - int(s[j]))
                seen[key] = seen.get(key, 0) + 1


def test_static_counts_equal_counting_both_strands_of_the_other_sequence():
    test = _codes(_seq(120, seed=5))
    train = _codes(_seq(300, seed=6, alphabet="ACGTTG"))
    k = 2
    state_t = np.array([i % 2 for i in range(len(test))], dtype=np.uint8)
    state_r = np.array([i % 3 == 0 for i in range(len(train))], dtype=np.uint8)
    n_ca, n_c = cz.static_counts(
        test, train, k, test_state=state_t, train_state=state_r, train_rc_state=state_r, state_bits=1
    )
    table: dict = {}
    for i in range(len(train)):
        key = (_ctx(train, i, k), int(state_r[i]), int(train[i]))
        table[key] = table.get(key, 0) + 1
    for j in range(len(train) - k):
        rc = tuple(3 - int(train[j + m]) for m in range(1, k + 1))
        key = (rc, int(state_r[j]), 3 - int(train[j]))
        table[key] = table.get(key, 0) + 1
    for i in range(len(test)):
        c, st = _ctx(test, i, k), int(state_t[i])
        assert n_ca[i] == table.get((c, st, int(test[i])), 0)
        assert n_c[i] == sum(table.get((c, st, a), 0) for a in range(4))


def test_kt_closed_form_equals_the_sequential_code():
    s = _codes(_seq(300, seed=9))
    k = 2
    n_ca, n_c = cz.adaptive_counts(s, k, revcomp=False)
    sequential = -np.log2((n_ca + 0.5) / (n_c + 2.0)).sum()
    key = cz.contexts(s, k).astype(np.int64) * 4 + s
    counts = np.bincount(key, minlength=4 ** (k + 1)).reshape(-1, 4)
    assert cz.kt_bits(counts, 0.5) == pytest.approx(sequential, rel=1e-9)


def test_lgamma_fast_is_exact_enough():
    x = np.array([0.5, 3.0, 255.5, 256.0, 1e3, 1e7])
    want = [math.lgamma(v) for v in x]
    assert np.allclose(cz._lgamma_fast(x), want, rtol=1e-12, atol=1e-9)


def test_mixture_is_a_distribution_and_causal():
    rng = np.random.default_rng(0)
    n, m = 200, 3
    probs = rng.dirichlet(np.ones(4), size=(m, n))  # every model's full distribution per base
    s = rng.integers(0, 4, n)
    base = [cz.quantise(-np.log2(probs[j, np.arange(n), s])) for j in range(m)]
    out = cz.mix(base, window=16, beta=1.0, chunk=64)
    i = 150
    total = 0.0
    for a in range(4):  # same history, each possible base at i
        models = []
        for j in range(m):
            x = base[j].copy()
            x[i] = cz.quantise(-np.log2(probs[j, i, a]))
            models.append(x)
        alt = cz.mix(models, window=16, beta=1.0, chunk=64)
        assert np.array_equal(alt[:i], out[:i])  # the past does not see the base
        total += 2.0 ** -float(alt[i])
    assert total == pytest.approx(1.0, abs=2e-3)  # quantisation of the model codes only


def test_mixture_follows_the_better_model():
    n = 5000
    good = np.full(n, cz.quantise(0.1), dtype=np.uint16)
    bad = np.full(n, cz.quantise(4.0), dtype=np.uint16)
    bits = cz.mix([bad, good], window=32, beta=1.0)
    assert bits[100:].mean() < 0.15


def test_annotation_costs():
    assert cz.delta_bits([1]) == 1
    assert cz.delta_bits([2, 3]) == 8  # 4 bits each
    seq = [0, 0, 1, 0]
    counts = [0, 0]
    bits = 0.0
    for x in seq:
        bits -= math.log2((counts[x] + 0.5) / (sum(counts) + 1.0))
        counts[x] += 1
    assert cz.categorical_bits(seq, 2) == pytest.approx(bits)
    assert cz.interval_bits([10, 30], [20, 35]) > 0
    assert cz.names_bits(["AluY"]) == 40


def test_copy_model_predicts_a_duplicated_segment():
    rng = random.Random(4)
    unit = _seq(3000, seed=11)
    mutated = "".join(c if rng.random() > 0.02 else rng.choice("ACGT") for c in unit)
    t, p = _codes(mutated), _codes(unit)
    off = cz.copy_offsets(t, p)
    usable = off >= 0
    x, bits, hits = cz.copy_codes(t, p, off, usable)
    assert len(x) > 2900
    assert hits / len(x) > 0.95
    assert bits.mean() < 0.5
    # a reversed-complement partner is read on its own strand by the caller
    rc = np.where(p < 4, 3 - p, 4).astype(np.uint8)[::-1]
    back = np.where(rc < 4, 3 - rc, 4).astype(np.uint8)[::-1]
    assert np.array_equal(back, p)


def test_gc_state_reads_only_the_past():
    s = _codes(_seq(3000, seed=2))
    st = cz.gc_state(s, window=100)
    s2 = s.copy()
    s2[2000:] = 1  # change the future
    st2 = cz.gc_state(s2, window=100)
    assert np.array_equal(st[:2001], st2[:2001])
    assert st.min() >= 1 and st.max() <= 18
    islands = _codes("CG" * 300)
    assert cz.gc_state(islands, window=100)[-1] == 5 * 3 + 2 + 1  # GC-rich, CpG island-like


def test_coding_states_phase_and_strand():
    c = cz.Chromosome.from_text("t", "N" * 5 + _seq(40))
    fwd, rc = cz.coding_states(c, [(False, [(5, 11, 0)]), (True, [(20, 26, 1)])])
    assert fwd[:6].tolist() == [1, 2, 3, 1, 2, 3]
    # minus strand, phase 1: counted from the segment's high end, one base skipped first
    assert fwd[15:21].tolist()[::-1] == [4 + (o - 1) % 3 for o in range(6)]
    assert rc[:3].tolist() == [4, 5, 6]


def test_pack_and_general_compressors():
    s = _codes(_seq(4000, seed=8))
    assert len(cz.pack2(s)) == 1000
    out = cz.general_compressors(s, ascii_xz=False)
    assert 1.5 < out["bzip2_packed"]["bits_per_base"] < 3.0


def test_models_on_a_repetitive_sequence_beat_two_bits():
    unit = _seq(500, seed=21)
    text = _seq(2000, seed=22) + unit + _seq(1000, seed=23) + unit
    s = _codes(text)
    adaptive = cz.adaptive_model(s, 12).astype(float) / cz.SCALE
    tail = adaptive[3500 + 12 :].mean()
    assert tail < 1.0  # the second copy is predicted
    train = _codes(_seq(5000, seed=24))
    naive = cz.static_model(s, train, 2).astype(float) / cz.SCALE
    assert naive.mean() == pytest.approx(2.0, abs=0.1)  # random sequence: nothing to learn


def test_plan_stacks_names_every_question():
    plan = cz.plan_stacks(["repeats", "coding"])
    assert ("naive", "adaptive", "coding") in plan["over_generic"]
    assert plan["over_repeat_aware"] == [("naive", "adaptive", "repeats", "coding")]
    assert plan["leave_one_out"][0] == ("naive", "adaptive", "coding")


def test_direct_and_sorted_counts_agree(monkeypatch):
    test = _codes(_seq(500, seed=31))
    train = _codes(_seq(900, seed=32, alphabet="AACGT"))
    st_t = (np.arange(len(test)) % 3).astype(np.uint8)
    st_r = (np.arange(len(train)) % 3).astype(np.uint8)
    kw = {"test_state": st_t, "train_state": st_r, "train_rc_state": st_r, "state_bits": 2}
    direct = cz.static_counts(test, train, 4, **kw)
    monkeypatch.setattr(cz, "DIRECT_BITS", 0)
    sorted_ = cz.static_counts(test, train, 4, **kw)
    assert np.array_equal(direct[0], sorted_[0]) and np.array_equal(direct[1], sorted_[1])


def test_an_inactive_model_takes_no_weight():
    rng = np.random.default_rng(1)
    a = cz.quantise(rng.uniform(0.5, 3.0, 3000))
    b = cz.quantise(rng.uniform(0.5, 3.0, 3000))
    off = np.full(3000, cz.INACTIVE, dtype=np.uint16)
    assert np.allclose(cz.mix([a, b, off], 16, 0.5, chunk=700), cz.mix([a, b], 16, 0.5, chunk=700))
    part = b.copy()
    part[:1500] = cz.INACTIVE
    bits = cz.mix([a, part], 16, 0.5)
    assert np.allclose(bits[:1500], cz.mix([a], 16, 0.5)[:1500])


def test_integer_code_is_proper_and_cheaper_than_delta_on_similar_values():
    values = [300, 310, 290, 305] * 50
    assert cz.integer_bits(values) < cz.delta_bits(values)
    assert cz.integer_bits([1]) == pytest.approx(-math.log2(0.5 / 32))


def test_training_set_paints_each_part():
    a = cz.Chromosome.from_text("a", "ACGTNNACGT")
    b = cz.Chromosome.from_text("b", "GGGG")
    train = cz.Chromosome.training_set([a, b])
    assert train.n == 12 and [p.name for p in train.parts] == ["a", "b"]

    def state(c):
        st = np.full(c.n, 1 if c.name == "a" else 2, dtype=np.uint8)
        return st, st + 1

    t, r, rc = cz._states(a, train, state)
    assert r.tolist() == [1] * 8 + [2] * 4 and rc.tolist() == [2] * 8 + [3] * 4


def test_priming_counts_the_held_sequence_before_the_first_base():
    held = _codes(_seq(600, seed=41))
    s = _codes(_seq(300, seed=42) + _seq(600, seed=41)[:300])
    k = 12
    plain = cz.adaptive_model(s, k).astype(float) / cz.SCALE
    primed = cz.adaptive_model(s, k, prime=held).astype(float) / cz.SCALE
    assert primed[312:].mean() < 0.5 < plain[312:].mean()
    plan = cz.plan_stacks(["duplications"], control=True)
    assert plan["control"][1] == ("naive", "adaptive_primed", "duplications")


def _loop_contexts(s, k):
    return [sum(int(s[i - j]) << (2 * (j - 1)) for j in range(1, k + 1) if i - j >= 0) for i in range(len(s))]


def test_doubling_contexts_match_a_loop_at_high_orders():
    s = _codes(_seq(120, seed=51))
    for k in (7, 13, 24):
        assert cz.contexts(s, k).tolist() == _loop_contexts(s, k)
        rc = cz.rc_contexts(s, k)
        for j in range(len(s) - k):
            want = sum((3 - int(s[j + m])) << (2 * (m - 1)) for m in range(1, k + 1))
            assert int(rc[j]) == want


def test_a_segment_reads_the_same_counts_as_the_whole(monkeypatch, tmp_path):
    s = _codes(_seq(3000, seed=52, alphabet="AACGTT"))
    train = _codes(_seq(2000, seed=53))
    k, tb = 5, cz._time_bits(len(s))
    monkeypatch.setattr(cz, "DIRECT_BITS", 0)  # the sorted path
    whole = cz._CountSink(len(s))
    cz.count_spans(
        cz.Span(s, 0, len(s), "query", k), cz._static_events(train, k, None, None), 0, tb, False, whole
    )
    part = cz._CountSink(1000)
    cz.count_spans(
        cz.Span(s, 1500, 2500, "query", k), cz._static_events(train, k, None, None), 0, tb, False, part
    )
    assert np.array_equal(part.n_ca, whole.n_ca[1500:2500]) and np.array_equal(part.n_c, whole.n_c[1500:2500])
    # adaptive, with the count rows waiting in temporary files and many hash groups
    monkeypatch.setattr(cz, "SPILL_BYTES", 1)
    monkeypatch.setattr(cz, "COUNT_GROUP", 500)
    a_ca, a_c = cz.adaptive_counts(s, 12)
    out = np.zeros(len(s), dtype=np.uint16)
    cz.adaptive_model(s, 12, out=out, spill=tmp_path)
    assert np.array_equal(out, cz._codes(a_ca, a_c, cz.alpha_for(12)))
    assert not any(tmp_path.iterdir())  # temporaries gone


def test_view_models_equal_the_whole_chromosome_slice():
    whole = cz.Chromosome.from_text("t", "NN" + _seq(4000, seed=54) + "NNNN" + _seq(2000, seed=55))
    train = cz.Chromosome.from_text("r", _seq(3000, seed=56))
    view = whole.view(2500, 4500)
    assert view.owns(int(whole.gpos[2500])) and not view.owns(int(whole.gpos[4500]))
    for k in (3, 10):
        a = cz._state_models(whole, train, "m", (None, None, None), 0, (k,))[f"m{k}"]
        b = cz._state_models(view, train, "m", (None, None, None), 0, (k,))[f"m{k}"]
        assert np.array_equal(a[2500:4500], b)
    assert np.array_equal(cz._gc_with_lead(view), cz.gc_state(whole.s)[2500:4500])
    idx, ok = view.stream_index([int(whole.gpos[2600]), 0])
    assert ok.tolist() == [True, False] and idx[0] == 100


def test_one_pass_mixes_every_stack_as_separate_mixtures_would():
    rng = np.random.default_rng(3)
    n = 3000
    models = [cz.quantise(rng.uniform(0.2, 3.5, n)) for _ in range(5)]
    models[3][:1200] = cz.INACTIVE
    groups = {"a": models[:2], "b": models[2:4], "c": models[4:]}
    stacks = {"a": ("a",), "ab": ("a", "b"), "abc": ("a", "b", "c"), "bc": ("b", "c")}
    got = {k: np.zeros(n) for k in stacks}
    for lo, hi, bits in cz.mix_stacks(groups, stacks, 16, 0.5, chunk=700):
        for k in stacks:
            got[k][lo:hi] = bits[k]
    for k, key in stacks.items():
        want = cz.mix([m for g in key for m in groups[g]], 16, 0.5, chunk=700)
        assert np.allclose(got[k], want, atol=1e-4)


def test_disk_guard_refuses_to_cross_the_floor(tmp_path):
    guard = cz.DiskGuard(tmp_path / "not" / "yet", floor_gb=0.0)
    free = guard.check()
    assert free > 0 and guard.min_free_gb == free
    with pytest.raises(cz.DiskFloorError):
        cz.DiskGuard(tmp_path, floor_gb=free + 1e6).check()
    with pytest.raises(cz.DiskFloorError):
        guard.check(need_bytes=(free + 1) * 1e9)


def _fixture(tmp_path, monkeypatch):
    import gzip
    import json

    ref, res, cache = tmp_path / "ref", tmp_path / "res", tmp_path / "cache"
    for d in (ref, res, cache):
        d.mkdir()
    alu = _seq(300, seed=60)
    parts = []
    for i in range(20):
        parts.append(_seq(1500, seed=61 + i))
        parts.append(alu)
    text = "N" * 500 + "".join(parts)
    (ref / "chrA.fa").write_text(">chrA\n" + text + "\n")
    (ref / "chrB.fa").write_text(">chrB\n" + "".join(_seq(1500, seed=90 + i) + alu for i in range(12)) + "\n")
    for name, length, step in (("chrA", len(text), 1800), ("chrB", 12 * 1800, 1800)):
        off = 500 if name == "chrA" else 0
        rows = [
            f"{off + i * step + 1500}\t{off + i * step + 1800}\tSINE\tAlu\tAluY\t0.1"
            for i in range(length // step)
        ]
        with gzip.open(res / f"rmsk_{name}.bed.gz", "wt") as fh:
            fh.write("\n".join(rows) + "\n")
        with gzip.open(res / f"ccres_{name}.bed.gz", "wt") as fh:
            fh.write(f"{name}\t{off + 100}\t{off + 400}\tE1\tdELS\t0\n")
        with gzip.open(cache / f"cpgIslandExt_{name}.bed.gz", "wt") as fh:
            fh.write(f"{name}\t{off + 3000}\t{off + 3400}\n")
        blocks = [
            {"start": off, "end": off + 9000, "guess": {"tier": "neutral"}, "class": "unique_intergenic"},
            {
                "start": off + 9000,
                "end": off + 18000,
                "guess": {"tier": "constrained_unknown"},
                "class": "unique_intergenic",
            },
            {
                "start": off + 18000,
                "end": off + 30000,
                "guess": {"tier": "fossil"},
                "class": "interspersed_repeat_SINE",
            },
        ]
        (res / f"budget_{name}.json").write_text(json.dumps({"blocks": blocks}))
        unknown = [{"start": b["start"], "end": b["end"], "class": b["class"]} for b in blocks]
        (res / f"unknown_{name}.json").write_text(json.dumps({"blocks": unknown}))
    monkeypatch.setattr(cz, "CACHE", cache)
    return ref, res


def test_the_pass_does_not_depend_on_the_segment_size(tmp_path, monkeypatch):
    ref, res = _fixture(tmp_path, monkeypatch)
    kw = {"results_dir": res, "reference": ref, "spill": tmp_path / "spill", "floor_gb": 0.0}
    whole = cz.run_pass("chrA", "chrB", segment=10**9, **kw)
    split = cz.run_pass("chrA", "chrB", segment=14_000, **kw)
    assert split["cost"]["segments"] == 3 and whole["cost"]["segments"] == 1
    a = {r["stack"]: r["bits_per_base"] for r in whole["stacks"]}
    b = {r["stack"]: r["bits_per_base"] for r in split["stacks"]}
    assert a.keys() == b.keys() == cz.all_stacks(list(cz.LAYER_ORDER)).keys()
    for k in a:
        assert b[k] == pytest.approx(a[k], abs=0.01)  # only the mixing windows restart at a boundary
    assert a["naive+adaptive"] < a["naive"]  # the repeated segment is found
    full = whole["key_stacks"]["full"]
    tiers = {r["region"]: r for r in whole["regions"]["rows"]}
    chrom = tiers["chromosome"]
    assert chrom["interspersed_bits_per_base"][full] < 1.0 < chrom["unique_bits_per_base"][full]
    assert not (tmp_path / "spill").exists() or not any((tmp_path / "spill").iterdir())
    for x in (whole, split):
        (res / f"compress_pass_{x['chrom']}.json").write_text(__import__("json").dumps(x))
    genome = cz.rollup(["chrA", "chrZ"], res)
    assert genome["chromosomes"] == ["chrA"] and genome["missing"] == ["chrZ"]
    assert genome["regions"]["unique_bootstrap"][full]["blocks"]["neutral"] == 1
